"""Durable MS-aware chromatographic components.

Candidate rectangles are temporary calculation windows.  The objects in this
module are the scientific result: overlapping peak contributions, each with its
own spectrum, amount, position and confidence.  They are JSON serialisable so a
reviewed fit can survive a project reload without inventing hard pixel borders.
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from types import SimpleNamespace

import numpy as np
from skimage.measure import find_contours

from .candidates import (CandidateRegion, CandidateSet, CandidateSettings, candidate_cube,
                         candidate_regions, find_candidates)
from .overlap import BlobCube, OverlapFit, assess_components, component_profiles, fit_overlap
from .peaks import PeakSettings, make_frame


SCHEMA_VERSION = 1


@dataclass(frozen=True)
class ComponentDetectionSettings:
    candidates: CandidateSettings = field(default_factory=CandidateSettings)
    link_t2_s: float = 0.25
    link_modulations: int = 1
    pad_t2_s: float = 0.4
    pad_modulations: int = 4
    model: str = "gaussian"
    max_components: int = 4
    n_starts: int = 3
    background_order: int = 1
    low_confidence_ratio: float = 0.30
    min_bic_gap: float = 10.0
    min_explained_fraction: float = 0.80
    screening_ions: int = 20
    max_ions: int = 60
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
        """Return NaN-separated global-time contours for one component."""
        contribution = self.contribution_map(component_index)
        maximum = float(contribution.max()) if contribution.size else 0.0
        if maximum <= 0:
            return np.array([]), np.array([])
        xs, ys = [], []
        for contour in find_contours(contribution, fraction * maximum):
            x = self.t1_start + (self.col0 + contour[:, 1]) * self.period
            y = (self.row0 + contour[:, 0]) * self.pitch
            for index in range(max(len(x) - 1, 0)):
                xs.extend((x[index], x[index + 1], np.nan))
                ys.extend((y[index], y[index + 1], np.nan))
        return np.asarray(xs), np.asarray(ys)


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
                   min_explained_fraction: float = -np.inf
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
    unresolved = bool(low or not fit.converged or ambiguous_order or poor_explanation)
    reasons = []
    if not fit.converged:
        reasons.append("fit did not fully converge")
    if low:
        reasons.append(f"{len(low)} weak component{'s' if len(low) != 1 else ''}")
    if ambiguous_order:
        reasons.append(f"component count ambiguous (BIC gap {bic_gap:.1f})")
    if poor_explanation:
        reasons.append(f"fit explains only {100*explained:.0f}% of local variation")
    return ComponentGroup(
        id=group_id, region_index=int(region_index), components=components,
        model=fit.model, row0=int(cube.row0), col0=int(cube.col0),
        n_rows=int(cube.data.shape[0]), n_cols=int(cube.data.shape[1]),
        pitch=float(cube.pitch), period=float(cube.period), t1_start=float(cube.t1_start),
        bic=float(fit.bic), weighted_residual=weighted_residual,
        explained_fraction=float(explained), converged=bool(fit.converged),
        unresolved=unresolved, warning="; ".join(reasons),
        template_source=getattr(fit, "template_source", None), n_starts=int(fit.n_starts),
        bic_gap=float(bic_gap))


def fit_candidate_region(run, fold, candidates: CandidateSet, region: CandidateRegion,
                         region_index: int, run_name: str,
                         settings: ComponentDetectionSettings = ComponentDetectionSettings()
                         ) -> ComponentGroup:
    """Fit and assess all plausible overlapping components in one region."""
    screening_cube = candidate_cube(
        run, fold, candidates, region,
        max_ions=min(settings.screening_ions, settings.max_ions))
    ordered = sorted(region.candidate_indices,
                     key=lambda index: float(candidates.score[index]), reverse=True)
    hints = [SimpleNamespace(row=float(candidates.row[index] - region.row0),
                             col=float(candidates.col[index] - region.col0))
             for index in ordered]
    max_k = max(1, min(int(settings.max_components), len(hints)))
    fits = [fit_overlap(screening_cube, k, model=settings.model, apex_hint=hints,
                        include_background=True, background_order=settings.background_order,
                        n_starts=1, seed=region_index + k)
            for k in range(1, max_k + 1)]
    fits.sort(key=lambda result: result.bic)
    screened = fits[0]
    bic_gap = float(fits[1].bic - screened.bic) if len(fits) > 1 else np.inf
    # Model-order screening needs only the strongest diagnostic ions.  Once k
    # is chosen, refit that model on the wider ion set so stored spectra and
    # amounts are not restricted to the cheap screening subset.
    cube = candidate_cube(run, fold, candidates, region, max_ions=settings.max_ions)
    refined_hints = [SimpleNamespace(row=component.row0, col=component.col0)
                     for component in screened.components]
    fit = fit_overlap(cube, screened.k, model=settings.model, apex_hint=refined_hints,
                      include_background=True, background_order=settings.background_order,
                      n_starts=settings.n_starts, seed=region_index + 1000)
    return group_from_fit(run_name, region_index, fit, cube,
                          settings.low_confidence_ratio, bic_gap=bic_gap,
                          min_bic_gap=settings.min_bic_gap,
                          min_explained_fraction=settings.min_explained_fraction)


def detect_components(run, fold, run_name: str,
                      settings: ComponentDetectionSettings = ComponentDetectionSettings(),
                      floor: float | None = None, progress=None) -> tuple:
    """Find candidates and return automatically fitted overlapping components."""
    candidates = find_candidates(run, fold, settings.candidates, floor=floor)
    frame = make_frame(fold, PeakSettings(wrap_below_s=settings.candidates.wrap_below_s))
    regions = candidate_regions(
        candidates, (frame.rows, len(frame.mid_t1)),
        link_t2_s=settings.link_t2_s, link_modulations=settings.link_modulations,
        pad_t2_s=settings.pad_t2_s, pad_modulations=settings.pad_modulations)
    component_set = ComponentSet()
    if settings.automatically_separate:
        for index, region in enumerate(regions):
            component_set.replace_group(
                fit_candidate_region(run, fold, candidates, region, index, run_name, settings))
            if progress is not None:
                progress(index + 1, len(regions))
    return candidates, regions, component_set
