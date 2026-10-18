"""Automatic choice of how many overlapping peaks a local window contains.

A calculation window is only a bounded piece of the chromatogram. How many
fitted peaks live in it is decided here, from the data, and never by how many
apex seeds the candidate finder happened to put there. Seeds are demoted to
*proposals*, on the same footing as proposals drawn from the spectral residual.

The search starts from one peak and grows one peak at a time. A proposed peak
is accepted only when every one of these holds:

1. The larger model beats the current one on BIC while the current one is
   allowed to be broad. A single wide peak is therefore the null hypothesis.
2. The peak carries a mass spectrum different from every peak it overlaps
   (cosine similarity below ``distinct_cosine``), or it is spatially resolved from
   them (centres more than ``resolved_sigma`` widths apart). Two overlapping
   components with indistinguishable spectra are not sufficient evidence for
   two compounds: one non-Gaussian profile can also explain that signal.
3. It carries at least ``min_relative_volume`` of the window's signal.
4. Refitting on held-out spatial blocks improves reconstruction in most folds.

A backward pass then removes any peak made redundant by the others. Whatever
was refused or merged for want of spectral evidence is reported as a note, so
an unresolved shoulder is visible to the user rather than silently absorbed.

Model comparisons run on a handful of diagnostic ions (the screening cube). The
spectra used to judge distinctness and reported to the user come from the wider
ion set.

With adaptive shapes enabled, each count is compared as independent Gaussians
and correlated Gaussians. Correlation must earn the additional BIC penalty;
neighbouring compounds are never forced to have identical widths.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from types import SimpleNamespace

import numpy as np
from scipy import ndimage as ndi

from .overlap import (BlobCube, OverlapComponent, OverlapFit, fit_shapes_fast, predict_overlap,
                      refit_spectra)


@dataclass
class Selection:
    fit: OverlapFit                 # final fit on the wide ion set
    screening_fit: OverlapFit       # the same shapes on the screening ions
    bic_gap: float                  # BIC distance to the nearest alternative model order
    notes: list = field(default_factory=list)       # what the user should be told
    rejections: list = field(default_factory=list)  # every refused proposal, for diagnosis
    n_fits: int = 0                 # nonlinear fits run, a direct cost measure
    n_proposals: int = 0


def spectral_cosine(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.clip(a, 0.0, None), np.clip(b, 0.0, None)
    return float(np.dot(a, b) / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-12))


def separation_in_widths(a, b) -> float:
    """Distance in mean peak widths, accounting for fitted axis correlation."""
    sigma_row = max(0.5 * (a.sigma_row + b.sigma_row), 1e-6)
    sigma_col = max(0.5 * (a.sigma_col + b.sigma_col), 1e-6)
    dr, dc = (a.row0-b.row0)/sigma_row, (a.col0-b.col0)/sigma_col
    rho = np.clip(.5*((getattr(a, 'rho', None) or 0.0)
                         +(getattr(b, 'rho', None) or 0.0)), -.95, .95)
    return float(np.sqrt(max((dr*dr+dc*dc-2*rho*dr*dc)/(1-rho*rho), 0)))


def has_valley(a, b, depth: float = 0.95) -> bool:
    """True when the summed profile of two fitted peaks dips between their centres.

    We accept two same-spectrum peaks as separate components only when their
    combined profile shows two maxima. A tail gives a second component without
    a dip, however far away it sits. Heights follow from volume and widths.
    """
    def amplitude(p):
        rho = getattr(p, 'rho', None) or 0.0
        return p.volume / max(2.0 * np.pi * p.sigma_row * p.sigma_col
                              * np.sqrt(1-rho*rho), 1e-12)

    t = np.linspace(0.0, 1.0, 41)
    rows = a.row0 + t * (b.row0 - a.row0)
    cols = a.col0 + t * (b.col0 - a.col0)
    total = np.zeros_like(rows)
    for p in (a, b):
        dr, dc = (rows-p.row0)/p.sigma_row, (cols-p.col0)/p.sigma_col
        rho = getattr(p, 'rho', None) or 0.0
        total += amplitude(p)*np.exp(-.5*(dr*dr+dc*dc-2*rho*dr*dc)/(1-rho*rho))
    return bool(total[1:-1].min() < depth * min(total[0], total[-1]))


def _hint(row: float, col: float):
    return SimpleNamespace(row=float(row), col=float(col))


def _centres(fit: OverlapFit) -> list:
    return [_hint(c.row0, c.col0) for c in fit.components]


class _Search:
    """One window's model search; counts nonlinear fits as its cost measure."""

    def __init__(self, screen: BlobCube, full: BlobCube, seeds: list, settings):
        self.screen, self.full, self.seeds, self.s = screen, full, seeds, settings
        self.caps = (float(settings.pursuit_max_sigma_row), float(settings.pursuit_max_sigma_col))
        self.n_fits = 0
        self.rejections: list = []
        self.n_proposals = 0

    # ------------------------------------------------------------------ fitting
    def fit(self, cube: BlobCube, hints: list, slack=(1.0, 0.6)) -> OverlapFit:
        fits = []
        models = ("gaussian", "tilt") if self.s.adaptive_shapes else ("gaussian",)
        for model in models:
            self.n_fits += 1
            shaped = fit_shapes_fast(
                cube, len(hints), model=model, apex_hint=hints,
                background_order=self.s.background_order,
                max_nfev=self.s.fast_max_nfev, position_rows=slack[0], position_cols=slack[1],
                max_sigma_row=self.caps[0], max_sigma_col=self.caps[1])
            fits.append(refit_spectra(shaped, cube, include_background=True,
                                      background_order=self.s.background_order))
        # Retain the simpler shape unless correlation earns its extra parameters.
        if len(fits) == 2 and fits[1].bic + self.s.min_bic_gap < fits[0].bic:
            return fits[1]
        return fits[0]

    def with_shapes(self, fit: OverlapFit, keep: list, extra=None) -> OverlapFit:
        """The same fitted shapes with some peaks dropped and/or one more added.

        Spectra are not solved here (see :meth:`linear`). Screening a proposal or a
        removal this way costs one small linear solve instead of a nonlinear fit.
        """
        components = [fit.components[i] for i in keep]
        if extra is not None:
            sigma_row = float(np.mean([c.sigma_row for c in fit.components]))
            sigma_col = float(np.mean([c.sigma_col for c in fit.components]))
            components.append(OverlapComponent(
                row0=float(extra[0]), col0=float(extra[1]), sigma_row=sigma_row,
                sigma_col=sigma_col, tau=None, t1=0.0, t2=0.0, volume=0.0,
                spectrum=np.zeros(self.screen.data.shape[2]), row0_se=np.nan, col0_se=np.nan,
                rho=0.0 if fit.model == "tilt" else None))
        k = len(components)
        return replace(fit, k=k, components=components,
                       n_shape_parameters=k * (5 if fit.model == "tilt" else 4))

    def linear(self, shaped: OverlapFit) -> OverlapFit:
        return refit_spectra(shaped, self.screen, include_background=True,
                             background_order=self.s.background_order)
    def on_full(self, fit: OverlapFit) -> OverlapFit:
        return refit_spectra(fit, self.full, include_background=True,
                             background_order=self.s.background_order)

    def summit(self) -> object:
        """Brightest smoothed point away from the rim: where a lone peak starts."""
        image = ndi.gaussian_filter(self.screen.data.sum(axis=2), (1.5, 0.5), mode="nearest")
        rows, cols = image.shape
        margin_r, margin_c = min(2, rows // 4), min(1, cols // 4)
        inner = np.zeros_like(image, dtype=bool)
        inner[margin_r:rows - margin_r, margin_c:cols - margin_c] = True
        row, col = np.unravel_index(int(np.argmax(np.where(inner, image, -np.inf))), image.shape)
        return _hint(row, col)

    # ---------------------------------------------------------------- proposals
    def residual_proposals(self, fit: OverlapFit, n: int = 2) -> list:
        """Peaks suggested by the positive residual, from its leading spectral factor.

        The residual is noise-weighted and spatially smoothed. Its first
        singular pair gives one spatial profile and one spectrum that explain
        the most leftover signal jointly over the ions, which is more sensitive
        to a weak distinct compound than the largest single pixel.
        """
        data = self.screen.data
        predicted = predict_overlap(fit, self.screen)
        z = np.maximum(data - predicted, 0.0) / np.sqrt(np.clip(data, 1.0, None))
        z = ndi.gaussian_filter(z, (1.0, 0.5, 0.0), mode="nearest")
        rows, cols, n_ions = z.shape
        if n_ions < 2:
            return []
        u, s, _vt = np.linalg.svd(z.reshape(-1, n_ions), full_matrices=False)
        profile = np.abs(u[:, 0]).reshape(rows, cols) * s[0] / np.sqrt(n_ions)
        profile[:2] = profile[-2:] = 0.0
        profile[:, :1] = profile[:, -1:] = 0.0
        if profile.max() < self.s.residual_min_score:
            return []
        proposals, work = [], profile.copy()
        sigma_row = np.mean([c.sigma_row for c in fit.components])
        sigma_col = np.mean([c.sigma_col for c in fit.components])
        for _ in range(n):
            if work.max() < self.s.residual_min_score:
                break
            row, col = np.unravel_index(int(np.argmax(work)), work.shape)
            proposals.append(_hint(row, col))
            rr, cc = np.indices(work.shape, dtype=float)
            work[np.hypot((rr - row) / sigma_row, (cc - col) / sigma_col) < 1.5] = 0.0
        return proposals

    def proposals(self, fit: OverlapFit) -> list:
        """Unused seeds first (they already have ion-level support), then residual peaks."""
        chosen = []
        for seed in self.seeds:
            if all(np.hypot((seed.row - c.row0) / max(c.sigma_row, 1e-6),
                            (seed.col - c.col0) / max(c.sigma_col, 1e-6)) > 1.0
                   for c in fit.components):
                chosen.append(seed)
        chosen = chosen[:self.s.max_seed_proposals]
        return chosen + self.residual_proposals(fit)

    # ---------------------------------------------------------------- decisions
    def heldout_gain(self, small: list, big: list) -> tuple[float, int]:
        """Held-out weighted-RMSE gain of the larger model, per spatial fold."""
        cube, folds = self.screen, self.s.cv_folds
        rr, cc = np.indices(cube.mask.shape)
        groups = ((rr // 5) + cc) % folds
        gains = []
        for held_out in range(folds):
            test = cube.mask & (groups == held_out)
            train = cube.mask & ~test
            if test.sum() < 4 or train.sum() < 30:
                continue
            training = BlobCube(cube.data, train, cube.ions, cube.row0, cube.col0,
                                cube.pitch, cube.period, cube.t1_start)
            observed = cube.data[test]
            weight = 1.0 / np.sqrt(np.clip(observed, 1.0, None))
            errors = []
            for hints in (small, big):
                prediction = predict_overlap(self.fit(training, hints), cube)[test]
                errors.append(float(np.sqrt(np.mean(((prediction - observed) * weight) ** 2))))
            gains.append(1.0 - errors[1] / max(errors[0], 1e-12))
        if not gains:
            return 0.0, 0
        return float(np.median(gains)), int(np.sum(np.asarray(gains) > 0))

    def add_is_justified(self, current: OverlapFit, bigger: OverlapFit) -> tuple[bool, str]:
        s = self.s
        # `current` already may be as broad as the caps allow, so a lone wide
        # peak is the null model that the extra component has to beat.
        if bigger.bic + s.min_bic_gap >= current.bic:
            return False, f"no BIC support (gain {current.bic - bigger.bic:.0f})"
        volumes = np.asarray([c.volume for c in bigger.components])
        if volumes.min() < s.min_relative_volume * max(volumes.sum(), 1e-12):
            return False, f"component too small ({100 * volumes.min() / volumes.sum():.1f}%)"
        full = self.on_full(bigger)
        for i in range(bigger.k):
            for j in range(i + 1, bigger.k):
                if not self.pair_is_distinct(full.components[i], full.components[j]):
                    return False, (f"spectrum cosine {spectral_cosine(full.components[i].spectrum, full.components[j].spectrum):.2f} "
                                   f"at {separation_in_widths(bigger.components[i], bigger.components[j]):.1f} widths")
        gain, positive = self.heldout_gain(_centres(current), _centres(bigger))
        if positive * 3 < 2 * s.cv_folds or gain < s.cv_min_gain:
            return False, f"held-out gain {100 * gain:.1f}%"
        return True, ""

    def pair_is_distinct(self, a, b) -> bool:
        """Different spectra, or resolved: two maxima with a dip between them."""
        return (spectral_cosine(a.spectrum, b.spectrum) < self.s.distinct_cosine
                or (separation_in_widths(a, b) >= self.s.resolved_sigma and has_valley(a, b)))


def select_components(screen: BlobCube, full: BlobCube, seeds: list, settings) -> Selection:
    """Choose the number, positions and spectra of the peaks in one window."""
    search = _Search(screen, full, seeds, settings)
    max_k = max(1, int(settings.max_components))
    # Start from the window's brightest point, not from a seed, so a broad peak
    # scattered with ion apexes begins as one component.
    fit = search.fit(screen, [search.summit()], slack=(3.0, 1.0))
    rejected_bic: list = []
    while fit.k < max_k:
        proposals = search.proposals(fit)
        if not proposals:
            break
        search.n_proposals += len(proposals)
        everyone = list(range(fit.k))
        ranked = sorted(proposals, key=lambda p: search.linear(
            search.with_shapes(fit, everyone, (p.row, p.col))).bic)
        accepted = None
        for proposal in ranked[:settings.max_fitted_proposals]:
            bigger = search.fit(screen, _centres(fit) + [proposal], slack=(2.0, 0.8))
            ok, reason = search.add_is_justified(fit, bigger)
            if ok:
                accepted = bigger
                break
            search.rejections.append(reason)
            rejected_bic.append(bigger.bic)
        if accepted is None:
            break
        fit = accepted
    fit, merged = _prune(search, fit)
    gap = _order_gap(search, fit, rejected_bic)
    notes = _summarise_notes(search.rejections, merged)
    return Selection(fit=search.on_full(fit), screening_fit=fit, bic_gap=gap, notes=notes,
                     rejections=search.rejections, n_fits=search.n_fits,
                     n_proposals=search.n_proposals)


def _prune(search: _Search, fit: OverlapFit) -> tuple[OverlapFit, list]:
    """Drop peaks that duplicate another peak's spectrum without being resolved from it."""
    merged, changed = [], False
    while fit.k > 1:
        full = search.on_full(fit)
        redundant = None
        for i in range(fit.k):
            for j in range(i + 1, fit.k):
                if redundant is None and not search.pair_is_distinct(
                        full.components[i], full.components[j]):
                    weaker = i if fit.components[i].volume < fit.components[j].volume else j
                    redundant = (i, j, weaker)
        if redundant is not None:
            i, j, weaker = redundant
            cos = spectral_cosine(full.components[i].spectrum, full.components[j].spectrum)
            merged.append(f"two components with near-identical spectra (cosine {cos:.2f}) "
                          f"{separation_in_widths(fit.components[i], fit.components[j]):.1f} widths "
                          "apart were merged")
            fit = search.linear(search.with_shapes(
                fit, [n for n in range(fit.k) if n != weaker]))
            changed = True
            continue
        # A peak the others make unnecessary earns no BIC when removed.
        smaller = None
        for i in range(fit.k):
            trial = search.linear(search.with_shapes(fit, [n for n in range(fit.k) if n != i]))
            if trial.bic <= fit.bic + search.s.min_bic_gap:
                smaller = trial
                break
        if smaller is None:
            break
        merged.append("component removed: no BIC support once the others were fitted")
        fit, changed = smaller, True
    if changed:
        # Survivors were shaped jointly with the peaks just removed.
        fit = search.fit(search.screen, _centres(fit), slack=(3.0, 1.0))
    return fit, merged


def _order_gap(search: _Search, fit: OverlapFit, rejected_bic: list) -> float:
    """BIC distance between the chosen model order and its nearest alternative."""
    if fit.k > 1:
        weakest = int(np.argmin([c.volume for c in fit.components]))
        smaller = search.linear(search.with_shapes(
            fit, [n for n in range(fit.k) if n != weakest]))
        return float(smaller.bic - fit.bic)
    return float(min(rejected_bic) - fit.bic) if rejected_bic else float("inf")

def _summarise_notes(rejections: list, merged: list) -> list:
    """Keep the messages a user can act on: spectral merges and near-misses."""
    notes = list(dict.fromkeys(merged))
    spectral = [r for r in rejections if r.startswith("spectrum cosine")]
    if spectral:
        notes.append("an extra peak was proposed but not accepted (" + spectral[0]
                     + "); if two compounds are present they cannot be separated by their spectra")
    return notes
