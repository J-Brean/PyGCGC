"""Durable MS-aware chromatographic components.

Candidate rectangles are temporary calculation windows.  The objects in this
module are the scientific result: overlapping peak contributions, each with its
own spectrum, amount, position and confidence.  They are JSON serialisable so a
reviewed fit can survive a project reload without inventing hard pixel borders.
"""

from __future__ import annotations

import csv
import os
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import dataclass, field, asdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import ndimage as ndi
from skimage.measure import find_contours
from threadpoolctl import threadpool_limits

from .candidates import (CandidateRegion, CandidateSet, CandidateSettings, candidate_cube,
                         candidate_regions, find_candidates)
from .modelselect import select_components, separation_in_widths
from .overlap import (BlobCube, OverlapFit, _shape, assess_components, component_profiles,
                      fit_overlap, fit_shapes_fast, initialise_overlap, refit_spectra,
                      shape_params)
from .peaks import PeakSettings, make_frame


SCHEMA_VERSION = 1
_WORKER_THREAD_LIMITER = None


def _limit_component_worker_threads() -> None:
    """Prevent every fit process from starting a full BLAS thread pool."""
    global _WORKER_THREAD_LIMITER
    _WORKER_THREAD_LIMITER = threadpool_limits(limits=1)


@dataclass(frozen=True)
class ComponentDetectionSettings:
    candidates: CandidateSettings = field(default_factory=CandidateSettings)
    link_t2_s: float = 0.25
    link_modulations: int = 1
    pad_t2_s: float = 0.4
    pad_modulations: int = 4
    model: str = "gaussian"
    max_components: int = 4
    screen_component_counts: bool = True
    n_starts: int = 2
    background_order: int = 1
    screening_background_order: int = 0
    low_confidence_ratio: float = 0.30
    min_bic_gap: float = 10.0
    min_explained_fraction: float = 0.80
    max_restart_position_spread: float = 1.0
    screening_ions: int = 8
    max_ions: int = 30
    fast_max_sigma_row: float = 3.0
    fast_max_sigma_col: float = 1.2
    fast_max_nfev: int = 24
    # Model-order selection. "pursuit" grows each window's model from one peak and
    # accepts extra peaks only when spectra and held-out data justify them (see
    # modelselect.py); "seeds" is the older behaviour in which the number of
    # candidate apexes in a window caps the number of fitted peaks.
    model_selection: str = "pursuit"
    pursuit_max_sigma_row: float = 7.0
    pursuit_max_sigma_col: float = 2.5
    residual_min_score: float = 4.0
    max_seed_proposals: int = 2
    max_fitted_proposals: int = 2
    distinct_cosine: float = 0.90
    resolved_sigma: float = 2.5
    min_relative_volume: float = 0.02
    cv_folds: int = 3
    cv_min_gain: float = 0.02
    merge_windows: bool = True
    merge_link_sigma: float = 3.0
    duplicate_sigma: float = 1.0
    edge_margin_sigma: float = 1.0
    truncated_cosine: float = 0.70
    # What to do with a peak a calculation window cut through: "recentre" refits it in a
    # window built around it, "drop" discards it when another window covers the area,
    # "none" keeps it. Compared on the benchmark before choosing a default.
    edge_policy: str = "drop"
    merge_max_rows: int = 160
    merge_max_cols: int = 14
    workers: int = max(1, min(4, os.cpu_count() or 1))
    refine_selected: bool = False
    fast_proposals: bool = False
    automatically_separate: bool = True

    def to_dict(self) -> dict:
        result = asdict(self)
        result["candidates"] = self.candidates.to_dict()
        return result

    @classmethod
    def from_dict(cls, data: dict) -> "ComponentDetectionSettings":
        data = dict(data or {})
        data["candidates"] = CandidateSettings.from_dict(data.get("candidates", {}))
        return cls(**{key: value for key, value in data.items()
                      if key in cls.__dataclass_fields__})


@dataclass
class Component:
    id: str
    group_id: str
    t1: float
    t2: float
    volume: float
    mz: np.ndarray
    spectrum: np.ndarray
    model: str
    row0: float
    col0: float
    sigma_row: float
    sigma_col: float
    tau: float | None = None
    tau_row: float | None = None
    rho: float | None = None
    row0_se: float = np.nan
    col0_se: float = np.nan
    relative_to_largest: float = 1.0
    confidence: str = "provisional"
    reason: str = ""
    profile: np.ndarray = field(default_factory=lambda: np.array([]))

    def to_dict(self) -> dict:
        return {
            "id": self.id, "group_id": self.group_id, "t1": self.t1, "t2": self.t2,
            "volume": self.volume, "mz": self.mz.tolist(), "spectrum": self.spectrum.tolist(),
            "model": self.model, "row0": self.row0, "col0": self.col0,
            "sigma_row": self.sigma_row, "sigma_col": self.sigma_col,
            "tau": self.tau, "tau_row": self.tau_row, "rho": self.rho,
            "row0_se": self.row0_se, "col0_se": self.col0_se,
            "relative_to_largest": self.relative_to_largest,
            "confidence": self.confidence, "reason": self.reason,
            "profile": self.profile.tolist() if self.profile.size else [],
        }

    @classmethod
    def from_dict(cls, data: dict) -> "Component":
        values = dict(data)
        for key in ("mz", "spectrum", "profile"):
            values[key] = np.asarray(values.get(key, []), dtype=float)
        return cls(**{key: value for key, value in values.items()
                      if key in cls.__dataclass_fields__})


@dataclass
class ComponentGroup:
    id: str
    region_index: int
    components: list[Component]
    model: str
    row0: int
    col0: int
    n_rows: int
    n_cols: int
    pitch: float
    period: float
    t1_start: float
    bic: float
    weighted_residual: float
    explained_fraction: float
    converged: bool
    unresolved: bool
    warning: str = ""
    template_source: str | None = None
    n_starts: int = 1
    bic_gap: float = np.inf
    restart_position_spread: float = np.nan

    def to_dict(self) -> dict:
        return {
            "id": self.id, "region_index": self.region_index,
            "components": [component.to_dict() for component in self.components],
            "model": self.model, "row0": self.row0, "col0": self.col0,
            "n_rows": self.n_rows, "n_cols": self.n_cols, "pitch": self.pitch,
            "period": self.period, "t1_start": self.t1_start, "bic": self.bic,
            "weighted_residual": self.weighted_residual,
            "explained_fraction": self.explained_fraction,
            "converged": self.converged, "unresolved": self.unresolved,
            "warning": self.warning, "template_source": self.template_source,
            "n_starts": self.n_starts,
            "bic_gap": self.bic_gap,
            "restart_position_spread": self.restart_position_spread,
        }

    @classmethod
    def from_dict(cls, data: dict) -> "ComponentGroup":
        values = dict(data)
        values["components"] = [Component.from_dict(item)
                                for item in values.get("components", [])]
        return cls(**{key: value for key, value in values.items()
                      if key in cls.__dataclass_fields__})

    def contribution_map(self, component_index: int) -> np.ndarray:
        component = self.components[component_index]
        if component.profile.size:
            profile = component.profile
        else:
            profile = component_profiles_from_durable(component, (self.n_rows, self.n_cols))
        return profile * float(component.spectrum.sum())

    def contour(self, component_index: int, fraction: float = 0.10
                ) -> tuple[np.ndarray, np.ndarray]:
        """Global-time outline of one component at `fraction` of its own peak height.

        The outline is taken from the fitted shape itself, not traced on the pixel
        grid: one pixel is a whole modulation in t1, so a traced outline had only a
        few points across and looked angular. A Gaussian's outline is an exact
        ellipse. Other analytic shapes are traced on a grid 20 times finer in t1
        and 4 times finer in t2, and a measured template is interpolated onto that
        grid first. The outline is not cut off at the calculation window's edge.
        Returned as NaN-separated (start, end, NaN) triples, which the drawing code
        moves edge by edge across the t2 wrap.
        """
        component = self.components[component_index]
        if float(np.sum(component.spectrum)) <= 0 or not 0 < fraction < 1:
            return np.array([]), np.array([])
        curves = []                                   # (rows, cols) in window pixel units
        if self.model == "gaussian" and not component.profile.size:
            radius = np.sqrt(-2.0 * np.log(fraction))
            angle = np.linspace(0.0, 2.0 * np.pi, 97)
            curves.append((component.row0 + radius * component.sigma_row * np.sin(angle),
                           component.col0 + radius * component.sigma_col * np.cos(angle)))
        else:
            step_r, step_c = 0.25, 0.05
            if component.profile.size:
                profile = ndi.zoom(component.profile, (1 / step_r, 1 / step_c), order=3)
                offset = (0.0, 0.0)
                scale = ((profile.shape[0] - 1) / max(component.profile.shape[0] - 1, 1),
                         (profile.shape[1] - 1) / max(component.profile.shape[1] - 1, 1))
            else:
                rows = np.arange(-4.0, self.n_rows + 4.0, step_r)
                cols = np.arange(-2.0, self.n_cols + 2.0, step_c)
                rr, cc = np.meshgrid(rows, cols, indexing="ij")
                profile = _shape(self.model, rr.ravel(), cc.ravel(),
                                 shape_params(self.model, component)).reshape(rr.shape)
                offset, scale = (rows[0], cols[0]), (1 / step_r, 1 / step_c)
            maximum = float(profile.max()) if profile.size else 0.0
            if maximum <= 0:
                return np.array([]), np.array([])
            for trace in find_contours(profile, fraction * maximum):
                curves.append((offset[0] + trace[:, 0] / scale[0],
                               offset[1] + trace[:, 1] / scale[1]))
        xs, ys = [], []
        for rows, cols in curves:
            if len(rows) < 2:
                continue
            x = self.t1_start + (self.col0 + cols + 0.5) * self.period   # mid-modulation, as BlobCube.t1_of
            y = (self.row0 + rows) * self.pitch
            gap = np.full(len(x) - 1, np.nan)
            xs.append(np.column_stack((x[:-1], x[1:], gap)).ravel())
            ys.append(np.column_stack((y[:-1], y[1:], gap)).ravel())
        if not xs:
            return np.array([]), np.array([])
        return np.concatenate(xs), np.concatenate(ys)


@dataclass
class ComponentSet:
    groups: list[ComponentGroup] = field(default_factory=list)
    schema_version: int = SCHEMA_VERSION

    @property
    def components(self) -> list[Component]:
        return [component for group in self.groups for component in group.components]

    def group_for_region(self, region_index: int) -> ComponentGroup | None:
        return next((group for group in self.groups if group.region_index == region_index), None)

    def replace_group(self, group: ComponentGroup) -> None:
        self.groups = [old for old in self.groups if old.region_index != group.region_index]
        self.groups.append(group)
        self.groups.sort(key=lambda item: item.region_index)

    def to_dict(self) -> dict:
        return {"schema_version": self.schema_version,
                "groups": [group.to_dict() for group in self.groups]}

    @classmethod
    def from_dict(cls, data: dict) -> "ComponentSet":
        data = data or {}
        return cls(groups=[ComponentGroup.from_dict(item) for item in data.get("groups", [])],
                   schema_version=int(data.get("schema_version", SCHEMA_VERSION)))


def component_rows(component_set: ComponentSet, run_name: str = "") -> list[dict]:
    """Return one auditable, flat export row for every fitted component."""
    rows = []
    for group in component_set.groups:
        for component in group.components:
            order = np.argsort(component.spectrum)[::-1]
            top_ions = [f"{component.mz[index]:g}"
                        for index in order[:5] if component.spectrum[index] > 0]
            rows.append({
                "run": run_name,
                "component_id": component.id,
                "group_id": group.id,
                "group_region": group.region_index,
                "t1_s": component.t1,
                "t1_min": component.t1 / 60.0,
                "t2_s": component.t2,
                "estimated_intensity": component.volume,
                "relative_to_group_largest": component.relative_to_largest,
                "diagnostic_ions": ";".join(top_ions),
                "confidence": component.confidence,
                "component_reason": component.reason,
                "group_unresolved": group.unresolved,
                "group_warning": group.warning,
                "model": group.model,
                "explained_fraction": group.explained_fraction,
                "bic": group.bic,
                "bic_gap": group.bic_gap,
                "fit_converged": group.converged,
                "restart_position_spread_pixels": group.restart_position_spread,
            })
    return rows


def write_components_csv(component_set: ComponentSet, path,
                         run_name: str = "") -> Path:
    """Write peaks without hiding unresolved-fit and confidence warnings."""
    path = Path(path)
    rows = component_rows(component_set, run_name)
    fields = list(rows[0]) if rows else [
        "run", "component_id", "group_id", "group_region", "t1_s", "t1_min", "t2_s",
        "estimated_intensity", "relative_to_group_largest", "diagnostic_ions", "confidence",
        "component_reason", "group_unresolved", "group_warning", "model",
        "explained_fraction", "bic", "bic_gap", "fit_converged",
        "restart_position_spread_pixels"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path


def component_profiles_from_durable(component: Component,
                                    shape: tuple[int, int]) -> np.ndarray:
    """Reconstruct one normalised analytic profile from stored parameters."""
    from .overlap import OverlapComponent

    fitted = OverlapComponent(
        row0=component.row0, col0=component.col0,
        sigma_row=component.sigma_row, sigma_col=component.sigma_col,
        tau=component.tau, tau_row=component.tau_row, rho=component.rho,
        t1=component.t1, t2=component.t2, volume=component.volume,
        spectrum=component.spectrum, row0_se=component.row0_se,
        col0_se=component.col0_se)
    return component_profiles(component.model, [fitted], shape)[:, :, 0]


def group_from_fit(run_name: str, region_index: int, fit: OverlapFit,
                   cube: BlobCube, low_confidence_ratio: float = 0.30,
                   bic_gap: float = np.inf, min_bic_gap: float = 0.0,
                   min_explained_fraction: float = -np.inf,
                   max_restart_position_spread: float = np.inf
                   ) -> ComponentGroup:
    """Convert a local numerical fit into a durable component group."""
    assessments = assess_components(fit, low_confidence_ratio)
    group_id = f"{run_name}:region:{int(region_index)}"
    components = []
    for index, (fitted, assessment) in enumerate(zip(fit.components, assessments)):
        profile = (fit.profiles[:, :, index].copy()
                   if fit.model == "template" and fit.profiles.size else np.array([]))
        components.append(Component(
            id=f"{group_id}:component:{index}", group_id=group_id,
            t1=float(fitted.t1), t2=float(fitted.t2), volume=float(fitted.volume),
            mz=np.asarray(cube.ions, dtype=float).copy(),
            spectrum=np.asarray(fitted.spectrum, dtype=float).copy(), model=fit.model,
            row0=float(fitted.row0), col0=float(fitted.col0),
            sigma_row=float(fitted.sigma_row), sigma_col=float(fitted.sigma_col),
            tau=None if fitted.tau is None else float(fitted.tau),
            tau_row=None if fitted.tau_row is None else float(fitted.tau_row),
            rho=None if fitted.rho is None else float(fitted.rho),
            row0_se=float(fitted.row0_se), col0_se=float(fitted.col0_se),
            relative_to_largest=assessment.relative_to_largest,
            confidence=assessment.confidence, reason=assessment.reason, profile=profile))
    observed = cube.data[cube.mask]
    # Match the fitter's count-weighted objective.  An unweighted variance
    # fraction can become negative on real cleaned data because a few abundant
    # ions dominate it even though the fit deliberately downweights their shot
    # noise.  Zero signal is the auditable null model here; signed background
    # remains part of the fitted prediction.
    denominator = float(np.sum(observed ** 2 / np.clip(observed, 1.0, None)))
    explained = (1.0 - float(fit.weighted_rss) / denominator
                 if denominator > 0 else np.nan)
    weighted_residual = float(np.sqrt(fit.weighted_rss / max(observed.size, 1)))
    low = [component for component in components if component.confidence == "low"]
    ambiguous_order = np.isfinite(bic_gap) and bic_gap < min_bic_gap
    poor_explanation = np.isfinite(explained) and explained < min_explained_fraction
    unstable = (np.isfinite(fit.restart_position_spread)
                and fit.restart_position_spread > max_restart_position_spread)
    unresolved = bool(low or not fit.converged or ambiguous_order or poor_explanation or unstable)
    reasons = []
    if not fit.converged:
        reasons.append("fit did not fully converge")
    if low:
        reasons.append(f"{len(low)} weak component{'s' if len(low) != 1 else ''}")
    if ambiguous_order:
        reasons.append(f"component count ambiguous (BIC gap {bic_gap:.1f})")
    if poor_explanation:
        reasons.append(f"fit explains only {100*explained:.0f}% of local variation")
    if unstable:
        reasons.append(
            f"repeated fits move a component by {fit.restart_position_spread:.1f} pixels")
    return ComponentGroup(
        id=group_id, region_index=int(region_index), components=components,
        model=fit.model, row0=int(cube.row0), col0=int(cube.col0),
        n_rows=int(cube.data.shape[0]), n_cols=int(cube.data.shape[1]),
        pitch=float(cube.pitch), period=float(cube.period), t1_start=float(cube.t1_start),
        bic=float(fit.bic), weighted_residual=weighted_residual,
        explained_fraction=float(explained), converged=bool(fit.converged),
        unresolved=unresolved, warning="; ".join(reasons),
        template_source=getattr(fit, "template_source", None), n_starts=int(fit.n_starts),
        bic_gap=float(bic_gap), restart_position_spread=float(fit.restart_position_spread))


def fit_candidate_region(run, fold, candidates: CandidateSet, region: CandidateRegion,
                         region_index: int, run_name: str,
                         settings: ComponentDetectionSettings = ComponentDetectionSettings()
                         ) -> ComponentGroup:
    """Fit and assess all plausible overlapping components in one region."""
    ordered = sorted(region.candidate_indices,
                     key=lambda index: float(candidates.score[index]), reverse=True)
    hints = [SimpleNamespace(row=float(candidates.row[index] - region.row0),
                             col=float(candidates.col[index] - region.col0))
             for index in ordered]
    cube = candidate_cube(run, fold, candidates, region, max_ions=settings.max_ions)
    screening_cube = _subset_cube(cube, settings.screening_ions)
    return _fit_candidate_cubes(
        screening_cube, cube, hints, region_index, run_name, settings)


def _subset_cube(cube: BlobCube, max_ions: int) -> BlobCube:
    """Use the leading diagnostic ions without folding a region twice."""
    count = min(max(int(max_ions), 1), cube.data.shape[2])
    return BlobCube(cube.data[:, :, :count], cube.mask, cube.ions[:count],
                    cube.row0, cube.col0, cube.pitch, cube.period, cube.t1_start)


def _fit_selected(screening_cube: BlobCube, cube: BlobCube, hints: list,
                  region_index: int, run_name: str,
                  settings: ComponentDetectionSettings) -> ComponentGroup:
    """Fit one window with the number of peaks chosen from spectra and held-out data.

    ``hints`` are candidate apexes. They are proposals to the model search, not
    its model order: the search may keep fewer of them or add peaks that no
    candidate suggested.
    """
    selection = select_components(screening_cube, cube, hints, settings)
    group = group_from_fit(
        run_name, region_index, selection.fit, cube, settings.low_confidence_ratio,
        bic_gap=selection.bic_gap, min_bic_gap=settings.min_bic_gap,
        min_explained_fraction=settings.min_explained_fraction,
        max_restart_position_spread=settings.max_restart_position_spread)
    notes = list(selection.notes)
    if any(c.sigma_row >= settings.pursuit_max_sigma_row - 0.01 or
           c.sigma_col >= settings.pursuit_max_sigma_col - 0.01
           for c in selection.fit.components):
        group.unresolved = True
        notes.append("a fitted peak reached the width limit; its amount is unreliable")
    if notes:
        group.unresolved = group.unresolved or any(
            n.startswith(("an extra peak", "two components")) for n in notes)
        group.warning = "; ".join(filter(None, [group.warning] + notes))
    return group


def _fit_candidate_cubes(screening_cube: BlobCube, cube: BlobCube, hints: list,
                         region_index: int, run_name: str,
                         settings: ComponentDetectionSettings) -> ComponentGroup:
    """CPU-heavy, process-safe part of one region fit."""
    if (settings.model_selection == "pursuit" and settings.model == "gaussian"
            and not settings.refine_selected):
        return _fit_selected(screening_cube, cube, hints, region_index, run_name, settings)
    max_k = max(1, min(int(settings.max_components), len(hints)))
    if not settings.refine_selected:
        if max_k == 1:
            # There is no overlap to separate.  The candidate apex and its
            # local half-height widths already define the only peak; solve its
            # spectrum and background directly instead of launching a
            # nonlinear optimiser merely to make a sub-pixel adjustment.
            fit = refit_spectra(
                initialise_overlap(cube, 1, settings.model, hints), cube,
                include_background=True, background_order=settings.background_order)
            return group_from_fit(
                run_name, region_index, fit, cube, settings.low_confidence_ratio,
                bic_gap=np.inf, min_bic_gap=settings.min_bic_gap,
                min_explained_fraction=settings.min_explained_fraction,
                max_restart_position_spread=settings.max_restart_position_spread)
        component_counts = (range(1, max_k + 1)
                            if settings.screen_component_counts else (max_k,))
        screened_fits = []
        for k in component_counts:
            proposal = fit_shapes_fast(
                screening_cube, k, model=settings.model, apex_hint=hints,
                background_order=settings.screening_background_order,
                max_nfev=settings.fast_max_nfev,
                max_sigma_row=settings.fast_max_sigma_row,
                max_sigma_col=settings.fast_max_sigma_col)
            screened_fits.append(refit_spectra(
                proposal, screening_cube, include_background=True,
                background_order=settings.screening_background_order))
        screened_fits.sort(key=lambda result: result.bic)
        screened = screened_fits[0]
        bic_gap = (float(screened_fits[1].bic - screened.bic)
                   if len(screened_fits) > 1 else np.inf)
        fit = refit_spectra(screened, cube, include_background=True,
                            background_order=settings.background_order)
        group = group_from_fit(
            run_name, region_index, fit, cube, settings.low_confidence_ratio,
            bic_gap=bic_gap, min_bic_gap=settings.min_bic_gap,
            min_explained_fraction=settings.min_explained_fraction,
            max_restart_position_spread=settings.max_restart_position_spread)
        if settings.model == "gaussian" and any(
                component.sigma_row >= settings.fast_max_sigma_row - 0.01 or
                component.sigma_col >= settings.fast_max_sigma_col - 0.01
                for component in fit.components):
            group.unresolved = True
            note = ("a fitted peak reached the fast width limit; inspect or refit this "
                    "group before trusting its amount")
            group.warning = f"{group.warning}; {note}" if group.warning else note
        note = ("fast diagnostic-ion shape fit; open this group to run repeated-start "
                "refinement")
        group.warning = f"{group.warning}; {note}" if group.warning else note
        return group
    component_counts = (range(1, max_k + 1)
                        if settings.screen_component_counts else (max_k,))
    fits = [fit_overlap(screening_cube, k, model=settings.model, apex_hint=hints,
                        include_background=True,
                        background_order=settings.screening_background_order,
                        n_starts=1, seed=region_index + k)
            for k in component_counts]
    fits.sort(key=lambda result: result.bic)
    screened = fits[0]
    bic_gap = float(fits[1].bic - screened.bic) if len(fits) > 1 else np.inf
    # Model-order screening needs only the strongest diagnostic ions.  Once k
    # is chosen, refit that model on the wider ion set so stored spectra and
    # amounts are not restricted to the cheap screening subset.
    refined_hints = [SimpleNamespace(row=component.row0, col=component.col0)
                     for component in screened.components]
    fit = fit_overlap(cube, screened.k, model=settings.model, apex_hint=refined_hints,
                      include_background=True, background_order=settings.background_order,
                      n_starts=(1 if screened.k == 1 else settings.n_starts),
                      seed=region_index + 1000)
    group = group_from_fit(run_name, region_index, fit, cube,
                          settings.low_confidence_ratio, bic_gap=bic_gap,
                          min_bic_gap=settings.min_bic_gap,
                          min_explained_fraction=settings.min_explained_fraction,
                          max_restart_position_spread=settings.max_restart_position_spread)
    return group


def propose_candidate_group(run, fold, candidates: CandidateSet, region: CandidateRegion,
                            region_index: int, run_name: str,
                            settings: ComponentDetectionSettings) -> ComponentGroup:
    """Create cheap overlapping proposals; detailed separation remains on demand."""
    cube = candidate_cube(run, fold, candidates, region, max_ions=settings.max_ions)
    background = np.percentile(cube.data, 20.0, axis=(0, 1))
    sigma_row = max(settings.candidates.smooth_t2_s / max(cube.pitch, 1e-9), 1.0)
    sigma_col = max(settings.candidates.smooth_modulations, 0.5)
    group_id = f"{run_name}:region:{int(region_index)}"
    components = []
    for order, candidate_index in enumerate(region.candidate_indices):
        row = float(candidates.row[candidate_index] - region.row0)
        col = float(candidates.col[candidate_index] - region.col0)
        row_index = int(np.clip(round(row), 0, cube.data.shape[0] - 1))
        col_index = int(np.clip(round(col), 0, cube.data.shape[1] - 1))
        spectrum = np.maximum(cube.data[row_index, col_index] - background, 0.0)
        provisional = Component(
            id=f"{group_id}:component:{order}", group_id=group_id,
            t1=float(candidates.t1[candidate_index]),
            t2=float(candidates.t2[candidate_index]), volume=0.0,
            mz=cube.ions.copy(), spectrum=spectrum, model="gaussian",
            row0=row, col0=col, sigma_row=sigma_row, sigma_col=sigma_col,
            confidence="provisional",
            reason="candidate estimate; open the group for joint spectral separation")
        profile = component_profiles_from_durable(provisional, cube.mask.shape)
        provisional.volume = float(profile.sum() * spectrum.sum())
        components.append(provisional)
    largest = max((component.volume for component in components), default=0.0)
    for component in components:
        component.relative_to_largest = component.volume / largest if largest > 0 else 0.0
    return ComponentGroup(
        id=group_id, region_index=region_index, components=components,
        model="gaussian", row0=cube.row0, col0=cube.col0,
        n_rows=cube.data.shape[0], n_cols=cube.data.shape[1],
        pitch=cube.pitch, period=cube.period, t1_start=cube.t1_start,
        bic=np.nan, weighted_residual=np.nan, explained_fraction=np.nan,
        converged=False, unresolved=True,
        warning=("rapid candidate proposals with approximate intensities; open this group "
                 "to jointly separate overlapping spectra"),
        n_starts=0, bic_gap=np.nan, restart_position_spread=np.nan)


def _prepare_window(run, fold, candidates, region, index, run_name, settings):
    """Build the cubes and seed hints for one window (the parent process does this)."""
    ordered = sorted(region.candidate_indices,
                     key=lambda item: float(candidates.score[item]), reverse=True)
    hints = [SimpleNamespace(row=float(candidates.row[item] - region.row0),
                             col=float(candidates.col[item] - region.col0))
             for item in ordered]
    cube = candidate_cube(run, fold, candidates, region, max_ions=settings.max_ions)
    screening_cube = _subset_cube(cube, settings.screening_ions)
    return screening_cube, cube, hints, index, run_name, settings


def _fit_windows(run, fold, candidates, windows: list, run_name: str,
                 settings: ComponentDetectionSettings, cancel=None, progress=None
                 ) -> list[ComponentGroup]:
    """Fit ``windows`` (a list of (region_index, region)) and return their groups.

    At most a few dense cubes exist at once: a long acquisition can hold
    thousands of windows, and building every cube before fitting starts would
    create a large and unnecessary memory spike.
    """
    workers = min(max(1, int(settings.workers)), max(len(windows), 1))
    groups = []
    if workers == 1 or len(windows) < 2:
        for done, (index, region) in enumerate(windows, start=1):
            if cancel is not None:
                cancel()
            groups.append(_fit_candidate_cubes(
                *_prepare_window(run, fold, candidates, region, index, run_name, settings)))
            if progress is not None:
                progress(done, len(windows))
        return groups
    with ProcessPoolExecutor(
            max_workers=workers, initializer=_limit_component_worker_threads) as executor:
        pending, next_index, done = set(), 0, 0
        while next_index < len(windows) or pending:
            if cancel is not None:
                cancel()
            while next_index < len(windows) and len(pending) < 2 * workers:
                index, region = windows[next_index]
                pending.add(executor.submit(_fit_candidate_cubes, *_prepare_window(
                    run, fold, candidates, region, index, run_name, settings)))
                next_index += 1
            finished, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in finished:
                groups.append(future.result())
                done += 1
                if progress is not None:
                    progress(done, len(windows))
    return groups


def _reindex(group: ComponentGroup, index: int, run_name: str) -> None:
    group.region_index = int(index)
    group.id = f"{run_name}:region:{int(index)}"
    for number, component in enumerate(group.components):
        component.group_id = group.id
        component.id = f"{group.id}:component:{number}"


def _window_margin(group: ComponentGroup, component: Component, image_shape: tuple) -> float:
    """Distance of a peak centre from the nearest calculation-window edge, in widths.

    Edges that are the edge of the chromatogram itself do not count: only an
    artificial cut through the data makes a peak there unreliable.
    """
    top, left = group.row0, group.col0
    bottom, right = top + group.n_rows - 1, left + group.n_cols - 1
    distances = []
    if top > 0:
        distances.append(component.row0 / max(component.sigma_row, 1e-6))
    if bottom < image_shape[0] - 1:
        distances.append((group.n_rows - 1 - component.row0) / max(component.sigma_row, 1e-6))
    if left > 0:
        distances.append(component.col0 / max(component.sigma_col, 1e-6))
    if right < image_shape[1] - 1:
        distances.append((group.n_cols - 1 - component.col0) / max(component.sigma_col, 1e-6))
    return float(min(distances)) if distances else float("inf")


def _aligned_cosine(a: Component, b: Component) -> float:
    """Spectral cosine of two components that may have been fitted on different ions.

    Every window chooses its own ions, so spectra from different windows are only
    comparable after placing them on a common m/z axis (absent ions count as zero).
    """
    axis = np.union1d(a.mz, b.mz)
    left, right = np.zeros(axis.size), np.zeros(axis.size)
    left[np.searchsorted(axis, a.mz)] = np.clip(a.spectrum, 0.0, None)
    right[np.searchsorted(axis, b.mz)] = np.clip(b.spectrum, 0.0, None)
    return float(np.dot(left, right) / max(np.linalg.norm(left) * np.linalg.norm(right), 1e-12))

def _drop_duplicates(entries: list, image_shape: tuple,
                     settings: ComponentDetectionSettings) -> list:
    """Keep each peak once when several windows fitted it.

    A peak wider than its calculation window, or one whose ion apexes seeded
    several windows, is fitted more than once. The copy whose centre is furthest
    from its window's edge is the best conditioned, so it is kept. "Same peak"
    means the same place (within ``duplicate_sigma`` widths) and the same spectrum.
    """
    items = []
    for number, (_region, group) in enumerate(entries):
        for component in group.components:
            items.append((_window_margin(group, component, image_shape), component.volume,
                          number, component))
    items.sort(key=lambda item: (-item[0], -item[1]))
    kept, dropped = [], set()
    for margin, _volume, number, component in items:
        group = entries[number][1]
        here = SimpleNamespace(row0=group.row0 + component.row0, col0=group.col0 + component.col0,
                               sigma_row=component.sigma_row, sigma_col=component.sigma_col)
        duplicate = False
        for other_margin, other_number, other, other_component in kept:
            if other_number == number:
                continue
            # A truncated fit underestimates its own width, so measure the
            # distance in the wider of the two.
            separation = float(np.hypot(
                (here.row0 - other.row0) / max(here.sigma_row, other.sigma_row, 1e-6),
                (here.col0 - other.col0) / max(here.sigma_col, other.sigma_col, 1e-6)))
            # The spectrum of a fit cut by its window edge is contaminated by its
            # neighbours, so it is held to a looser standard than a well-centred one.
            needed = (settings.truncated_cosine if min(margin, other_margin) < settings.edge_margin_sigma
                      else settings.distinct_cosine)
            if (separation < settings.duplicate_sigma
                    and _aligned_cosine(component, other_component) >= needed):
                duplicate = True
                break
            if (settings.edge_policy == "drop" and margin < settings.edge_margin_sigma <= other_margin
                    and separation < settings.merge_link_sigma):
                duplicate = True
                break
        if duplicate:
            dropped.add(id(component))
        else:
            kept.append((margin, number, here, component))
    out = []
    for region, group in entries:
        group.components = [c for c in group.components if id(c) not in dropped]
        if group.components:
            largest = max(c.volume for c in group.components)
            for component in group.components:
                component.relative_to_largest = component.volume / largest if largest > 0 else 0.0
            out.append((region, group))
    return out


def _recentre_edge_components(run, fold, candidates, entries: list, image_shape: tuple,
                              run_name: str, settings: ComponentDetectionSettings,
                              cancel=None) -> list:
    """Refit peaks that a calculation window cut through, in a window centred on them.

    A peak whose centre lies within ``edge_margin_sigma`` widths of an artificial
    window edge was fitted from a truncated view of itself. Its position, width and
    spectrum are unreliable whether it is a fragment of a larger neighbour or a
    peak in its own right, so it is removed from that window and the same model
    search is run in a window built around it. A fragment becomes a duplicate of
    its parent and is discarded by :func:`_drop_duplicates`; a real peak now gets
    a proper fit.
    """
    jobs, centres = [], []
    for _region, group in entries:
        for component in group.components:
            if _window_margin(group, component, image_shape) >= settings.edge_margin_sigma:
                continue
            row, col = group.row0 + component.row0, group.col0 + component.col0
            if any(abs(row - r) < 0.5 * h and abs(col - c) < 0.5 * w for r, c, h, w in centres):
                continue        # an already-planned recentred window covers it
            half_rows = int(max(np.ceil(settings.pad_t2_s / group.pitch),
                                np.ceil(4.0 * component.sigma_row)))
            half_cols = int(max(settings.pad_modulations, np.ceil(3.0 * component.sigma_col)))
            row0, row1 = max(int(row) - half_rows, 0), min(int(row) + half_rows + 1, image_shape[0])
            col0, col1 = max(int(col) - half_cols, 0), min(int(col) + half_cols + 1, image_shape[1])
            inside = np.flatnonzero((candidates.row >= row0) & (candidates.row < row1)
                                    & (candidates.col >= col0) & (candidates.col < col1))
            jobs.append(CandidateRegion(tuple(int(i) for i in inside), row0, row1, col0, col1))
            centres.append((row, col, row1 - row0, col1 - col0))
    if not jobs:
        return entries
    refitted = _fit_windows(run, fold, candidates, list(enumerate(jobs)), run_name,
                            settings, cancel=cancel)
    refitted.sort(key=lambda group: group.region_index)
    survivors = []
    for region, group in entries:
        group.components = [c for c in group.components
                            if _window_margin(group, c, image_shape) >= settings.edge_margin_sigma]
        if group.components:
            survivors.append((region, group))
    return survivors + list(zip(jobs, refitted))

def _merge_interacting_windows(run, fold, candidates, entries: list, run_name: str,
                               settings: ComponentDetectionSettings, cancel=None,
                               max_rounds: int = 2) -> list:
    """Refit together pairs of windows whose distinct peaks reach into each other.

    Two peaks straddling the boundary between windows are otherwise fitted
    without knowledge of each other. Windows are joined pairwise, strongest
    interaction first, and only while the union stays small, so a dense
    chromatogram cannot chain into one enormous problem. The enlarged window is
    solved by the ordinary model search.

    ``entries`` is a list of (region, group); the merged list is returned.
    """
    for _ in range(max_rounds):
        strength = []
        for i in range(len(entries)):
            for j in range(i + 1, len(entries)):
                (region_i, group_i), (region_j, group_j) = entries[i], entries[j]
                if (region_i.col1 < region_j.col0 - 1 or region_j.col1 < region_i.col0 - 1
                        or region_i.row1 < region_j.row0 - 1 or region_j.row1 < region_i.row0 - 1):
                    continue
                closest = min(
                    separation_in_widths(
                        SimpleNamespace(row0=group_i.row0 + a.row0, col0=group_i.col0 + a.col0,
                                        sigma_row=a.sigma_row, sigma_col=a.sigma_col),
                        SimpleNamespace(row0=group_j.row0 + b.row0, col0=group_j.col0 + b.col0,
                                        sigma_row=b.sigma_row, sigma_col=b.sigma_col))
                    for a in group_i.components for b in group_j.components)
                if closest < settings.merge_link_sigma:
                    strength.append((closest, i, j))
        strength.sort()
        used, jobs = set(), []
        for _closest, i, j in strength:
            if i in used or j in used:
                continue
            regions = [entries[i][0], entries[j][0]]
            union = CandidateRegion(
                tuple(sorted({n for r in regions for n in r.candidate_indices})),
                min(r.row0 for r in regions), max(r.row1 for r in regions),
                min(r.col0 for r in regions), max(r.col1 for r in regions))
            if (union.row1 - union.row0 > settings.merge_max_rows
                    or union.col1 - union.col0 > settings.merge_max_cols):
                continue
            used.update((i, j))
            jobs.append(union)
        if not jobs:
            break
        merged = _fit_windows(run, fold, candidates, list(enumerate(jobs)), run_name,
                              settings, cancel=cancel)
        merged.sort(key=lambda group: group.region_index)
        entries = ([entry for n, entry in enumerate(entries) if n not in used]
                   + list(zip(jobs, merged)))
        entries.sort(key=lambda item: (item[0].col0, item[0].row0))
    return entries

def detect_components(run, fold, run_name: str,
                      settings: ComponentDetectionSettings = ComponentDetectionSettings(),
                      floor: float | None = None, progress=None, cancel=None) -> tuple:
    """Find candidates and return automatically fitted overlapping components.

    Candidate rectangles are calculation windows: they bound each local fit but
    neither cap the number of fitted peaks nor appear in the result. The
    returned ``regions`` are the final windows after any merging, in the order
    that ``ComponentGroup.region_index`` refers to.
    """
    candidates = find_candidates(run, fold, settings.candidates, floor=floor)
    frame = make_frame(fold, PeakSettings(wrap_below_s=settings.candidates.wrap_below_s))
    regions = candidate_regions(
        candidates, (frame.rows, len(frame.mid_t1)),
        link_t2_s=settings.link_t2_s, link_modulations=settings.link_modulations,
        pad_t2_s=settings.pad_t2_s, pad_modulations=settings.pad_modulations)
    component_set = ComponentSet()
    if not settings.automatically_separate:
        return candidates, regions, component_set
    if settings.fast_proposals:
        for done, region in enumerate(regions, start=1):
            component_set.replace_group(propose_candidate_group(
                run, fold, candidates, region, done - 1, run_name, settings))
            if cancel is not None:
                cancel()
            if progress is not None:
                progress(done, len(regions))
        return candidates, regions, component_set
    groups = _fit_windows(run, fold, candidates, list(enumerate(regions)), run_name,
                          settings, cancel=cancel, progress=progress)
    groups.sort(key=lambda group: group.region_index)
    entries = list(zip(regions, groups))
    if (settings.merge_windows and settings.model_selection == "pursuit"
            and settings.model == "gaussian" and not settings.refine_selected):
        shape = (frame.rows, len(frame.mid_t1))
        entries = _drop_duplicates(entries, shape, settings)
        entries = _merge_interacting_windows(run, fold, candidates, entries, run_name,
                                             settings, cancel=cancel)
        entries = _drop_duplicates(entries, shape, settings)
        if settings.edge_policy == "recentre":
            entries = _recentre_edge_components(run, fold, candidates, entries, shape, run_name,
                                                settings, cancel=cancel)
        entries = _drop_duplicates(entries, shape, settings)
        entries.sort(key=lambda item: (item[0].col0, item[0].row0))
        for index, (_region, group) in enumerate(entries):
            _reindex(group, index, run_name)
    for _region, group in entries:
        component_set.replace_group(group)
    return candidates, [region for region, _group in entries], component_set