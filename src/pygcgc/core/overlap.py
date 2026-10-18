"""✂️ Overlap toolkit: resolving compounds sharing one TIC blob.

Two peaks that share a footprint (`BlobSet.co_elution` flags the
candidates) get a single summed spectrum and a single volume from 🔍,
which is wrong for both of them. This module offers five ways to pull
them apart, cheapest first, each buying more separation power at the
cost of more assumptions:

1. **Unique-ion apex map** (`unique_ion_apexes`). No fitting at all: find
   where each of the blob's own significant ions actually peaks inside
   the footprint, and cluster those positions. One cluster means the
   ions agree (probably one compound); more than one is the same
   heuristic `BlobSet.co_elution` already flags, but here every ion's
   own apex is reported, which is what you need to draw a sensible cut.
2. **Geometric splits** (`split_by_line`, `split_by_watershed`,
   `split_by_valley`). Cut the footprint in two (or more) without
   modelling anything: a straight line through two points, a fresh
   watershed seeded wherever you say the separate apexes are, or the
   lowest-intensity path between two points (the "valley" between two
   peaks). `split_volumes` then just sums each region.
3. **Joint peak-shape fit** (`fit_overlap`, `fit_template_overlap`, `choose_k`). A proper
   quantitative deconvolution: k two-dimensional peak shapes (Gaussian,
   or an EMG along t1 for a tailing peak) each multiplied by their own
   mass spectrum, fitted jointly to the blob's own ion cube. The shape
   parameters (position, width, tailing) are the only nonlinear
   unknowns -- at every trial set of shapes, each component's spectrum
   is solved directly by (non-negative) least squares, a standard
   separable/variable-projection trick that keeps the nonlinear search
   small. The same approximate Poisson weights are used in the inner and
   outer fits. A measured-template alternative shifts copies of a trusted real
   peak profile instead of imposing an analytic curve. An optional signed local
   background is fitted per ion while component spectra remain non-negative,
   and repeated starts reduce dependence on one initial guess. `choose_k`
   fits k = 1..max_k and keeps whichever minimises BIC, counting both shape
   and spectral coefficients. Parameter uncertainty comes from the Jacobian at the
   solution (`OverlapComponent.row0_se`/`col0_se`), the standard
   asymptotic covariance for a least-squares fit.
4. **MCR-ALS** (`mcr_als`, `mcr_restarts`). No peak shape assumed at all: alternating
   non-negative least squares recovers a free-form elution profile and
   a free-form spectrum for each of k components directly from the
   data. `unimodal=True` additionally nudges each component's profile
   toward a single peak (`_enforce_unimodal`) -- a cheap capping-toward-
   the-peak approximation to a proper isotonic-regression projection,
   good enough to stop a component inventing a second bump, not an
   exact least-squares-optimal unimodal fit. Repeated starts report residual,
   explained variance, convergence and matched-spectrum stability.
5. **PARAFAC across runs** (`parafac`, `core_consistency`,
   `split_half_stability`, `build_run_tensor`). The other four methods
   all work on one run's blob in isolation. PARAFAC instead resolves a
   compound using the same feature matched across *several* runs at
   once (`build_run_tensor` builds the (run, elution position, ion)
   array), which is more identifiable in principle (Kruskal's
   trilinearity result: a genuine three-way structure has a
   near-unique decomposition without needing a peak shape or a
   non-negativity assumption to pin it down) but only if the same
   compounds are genuinely present with the same spectrum across those
   runs. `core_consistency` (Bro and Kiers, 2003) and
   `split_half_stability` are the two standard sanity checks on whether
   a chosen number of components is actually trustworthy, not just a
   number that happened to fit. PARAFAC2 (which would additionally
   tolerate each run having its own elution-time shift) is **not
   built**: real multi-run data to tell whether the shift is worth the
   extra machinery hasn't been available, the same reasoning that has
   kept 🫧's ROC threshold and ⛓️'s qualifier-ratio gate provisional.

Fitting has been compared in 32 windows from four calibration acquisitions.
Known-area injections on measured backgrounds test numerical recovery;
identified co-eluting standards are still needed for chemical validation.
See docs/validation/peak-fitting-2026-10-18.md and ROADMAP.md.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import least_squares, linear_sum_assignment, lsq_linear, nnls
from scipy.special import erf, log_ndtr
from skimage.graph import route_through_array
from skimage.segmentation import watershed

from .fold import Fold
from .peaks import BlobSet, Frame
from .run import Run


@dataclass(frozen=True)
class OverlapSettings:
    ion_share_min: float = 0.03      # an ion needs at least this share of the blob's spectrum to be used
    min_ions: int = 3                # fall back to the top this-many ions if too few clear the share
    max_components: int = 4          # BIC search ceiling for peak-model fitting
    model: str = "gaussian"          # or "emg"

    def to_dict(self) -> dict:
        return dict(ion_share_min=self.ion_share_min, min_ions=self.min_ions,
                    max_components=self.max_components, model=self.model)

    @classmethod
    def from_dict(cls, d: dict) -> "OverlapSettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


# ------------------------------------------------------------------ the cube
@dataclass
class BlobCube:
    """A small (rows, cols, ions) window around one blob, its own
    footprint mask, and enough geometry to read pixels back as (t1, t2).
    Everything else in this module is built from one of these, so a run
    and a fold are only ever read once per blob."""
    data: np.ndarray             # (rows, cols, n_ions)
    mask: np.ndarray             # (rows, cols) bool -- this blob's own footprint
    ions: np.ndarray             # m/z of each ion in `data`'s last axis
    row0: int                    # this window's offset into the extended frame
    col0: int
    pitch: float
    period: float
    t1_start: float

    def t1_of(self, col) -> np.ndarray:
        """t1 of a (fractional) column: a column sits at the middle of its modulation,
        the convention of `Frame.mid_t1`, candidates, legacy blobs and quantification."""
        return self.t1_start + (self.col0 + np.asarray(col, dtype=float) + 0.5) * self.period

    def t2_of(self, row) -> np.ndarray:
        return (self.row0 + np.asarray(row, dtype=float)) * self.pitch


def blob_cube(run: Run, fold: Fold, blobs: BlobSet, i: int,
             settings: OverlapSettings = OverlapSettings()) -> BlobCube:
    """The ion cube around blob `i`, for whichever of its own ions
    actually carry its signal: everything at least `ion_share_min` of
    the blob's summed spectrum, or the top `min_ions` if fewer than
    that clear the bar (a faint or very co-eluted blob can have no ion
    individually over the share threshold)."""
    rows, cols = blobs.footprint(i)
    if rows.size == 0:
        raise ValueError("blob has an empty footprint")
    r0, r1, c0, c1 = int(rows.min()), int(rows.max()), int(cols.min()), int(cols.max())
    spectrum = blobs.spectra[i]
    total = float(spectrum.sum())
    share = spectrum / total if total > 0 else spectrum
    idx = np.flatnonzero(share >= settings.ion_share_min)
    if idx.size < max(settings.min_ions, 2):
        order = np.argsort(spectrum)[::-1]
        idx = order[spectrum[order] > 0][:max(settings.min_ions, 2)]
    ions = blobs.mz[idx]
    frame = Frame(fold=fold, wrap_rows=blobs.wrap_rows)
    cube = frame.extend(fold.cube(run, mz=ions, dtype=np.float64))
    data = cube[r0:r1 + 1, c0:c1 + 1, :]
    mask = blobs.labels[r0:r1 + 1, c0:c1 + 1] == i + 1
    return BlobCube(data=data, mask=mask, ions=ions, row0=r0, col0=c0,
                    pitch=blobs.pitch, period=blobs.period, t1_start=blobs.t1_start)


# ------------------------------------------------------------- 1. unique ions
@dataclass
class IonApex:
    ion: float
    row: float
    col: float
    t1: float
    t2: float
    cluster: int


@dataclass
class UniqueIonMap:
    apexes: list
    n_components: int


def unique_ion_apexes(cube: BlobCube, gap: float = 1.2) -> UniqueIonMap:
    """Where each of the cube's own ions actually peaks inside the
    footprint, clustered (single-linkage, `gap` in units of one
    modulation along t1 and two rows along t2, matching
    `PeakSettings.coelution_gap`). More than one cluster is the first,
    cheapest sign of more than one compound under the footprint."""
    from scipy.cluster.hierarchy import fclusterdata
    pts, hits = [], []
    for k in range(cube.data.shape[2]):
        vals = np.where(cube.mask, cube.data[:, :, k], -np.inf)
        if not np.any(np.isfinite(vals)) or float(np.max(vals)) <= 0:
            continue
        r, c = np.unravel_index(int(np.argmax(vals)), vals.shape)
        pts.append((float(c), float(r) / 2.0))
        hits.append((float(cube.ions[k]), r, c))
    if len(pts) < 2:
        clusters = np.ones(len(pts), dtype=int)
    else:
        clusters = fclusterdata(np.asarray(pts), t=gap, criterion="distance")
    apexes = [IonApex(ion=ion, row=float(r), col=float(c), t1=float(cube.t1_of(c)), t2=float(cube.t2_of(r)),
                      cluster=int(cl))
             for (ion, r, c), cl in zip(hits, clusters)]
    n = int(max((a.cluster for a in apexes), default=1))
    return UniqueIonMap(apexes=apexes, n_components=n)


# ------------------------------------------------------------- 2. geometric splits
def split_by_line(cube: BlobCube, p1: tuple, p2: tuple) -> np.ndarray:
    """Split the footprint by a straight line through two (row, col)
    points (fractional is fine). Returns a label image: 0 outside the
    footprint, 1 or 2 either side of the line."""
    (r1, c1), (r2, c2) = p1, p2
    rows, cols = np.mgrid[0:cube.mask.shape[0], 0:cube.mask.shape[1]]
    side = (c2 - c1) * (rows - r1) - (r2 - r1) * (cols - c1)
    labels = np.where(side >= 0, 1, 2)
    labels[~cube.mask] = 0
    return labels.astype(np.int32)


def split_by_watershed(cube: BlobCube, seeds: list) -> np.ndarray:
    """Split the footprint with a fresh watershed grown from the given
    seed pixels [(row, col), ...], on the cube's own summed image --
    the same downhill-growth idea `peaks.find_blobs` uses for a whole
    run, confined to one footprint and seeded by hand (e.g. from
    `unique_ion_apexes`) instead of automatically."""
    img = cube.data.sum(axis=2)
    markers = np.zeros(cube.mask.shape, dtype=np.int32)
    for k, (r, c) in enumerate(seeds, start=1):
        markers[int(round(r)), int(round(c))] = k
    return watershed(-img, markers, mask=cube.mask).astype(np.int32)


def split_by_valley(cube: BlobCube, p1: tuple, p2: tuple) -> tuple:
    """The lowest-intensity route between two points inside the
    footprint (e.g. two ion apexes) -- the "valley" a watershed would
    grow into on its own, made explicit so a person can see or override
    it. Returns the path's (rows, cols)."""
    img = cube.data.sum(axis=2).astype(np.float64)
    cost = np.where(cube.mask, img, np.inf)
    indices, _cost = route_through_array(cost, tuple(int(round(v)) for v in p1),
                                         tuple(int(round(v)) for v in p2), fully_connected=True)
    rows = np.array([r for r, c in indices])
    cols = np.array([c for r, c in indices])
    return rows, cols


def split_volumes(cube: BlobCube, labels: np.ndarray) -> dict:
    """Summed spectrum, volume and centroid position for each labelled
    region (label 0 is ignored). Volumes are on the cube's own ions
    only (see `blob_cube`), not the run's full spectrum -- fine as long
    as `ion_share_min` was low enough to capture most of the blob's
    real signal; check `cube.ions` if a number looks low."""
    out = {}
    for lab in sorted(int(x) for x in np.unique(labels) if x > 0):
        rows, cols = np.nonzero(labels == lab)
        spectrum = cube.data[rows, cols, :].sum(axis=0)
        weight = cube.data[rows, cols, :].sum(axis=1)
        wsum = float(weight.sum())
        r_c = float(np.average(rows, weights=weight)) if wsum > 0 else float(rows.mean())
        c_c = float(np.average(cols, weights=weight)) if wsum > 0 else float(cols.mean())
        out[lab] = dict(volume=float(spectrum.sum()), spectrum=spectrum,
                        t1=float(cube.t1_of(c_c)), t2=float(cube.t2_of(r_c)))
    return out


# ------------------------------------------------------------- 3. peak-shape fit
MODELS = ("gaussian", "emg", "emg_t2", "emg2", "tilt", "emg_t2_tilt", "cell", "emg_t2_cell")
# gaussian: centre and width in each direction. emg: an exponential tail along
# t1 (the columns). emg_t2: a tail along t2 (the rows), the usual tailing
# direction for a second-dimension column. emg2: a tail in both. tilt: a
# Gaussian whose t2 apex drifts with t1 (a correlation `rho` between the axes),
# which real peaks do because each modulation is a separate injection.
# emg_t2_tilt: a t2 tail on a tilted peak. cell and emg_t2_cell: as gaussian and
# emg_t2 but with the t1 profile integrated over each column (see
# `_cell_gaussian_1d`), for peaks as narrow as a modulation or narrower. Extra
# parameters follow the four common ones, in the order listed here.
_EXTRAS = {"gaussian": (), "emg": ("tau",), "emg_t2": ("tau",), "emg2": ("tau", "tau"),
           "tilt": ("rho",), "emg_t2_tilt": ("tau", "rho"), "cell": (), "emg_t2_cell": ("tau",)}
_N_SHAPE = {m: 4 + len(e) for m, e in _EXTRAS.items()}


def _n_shape_params(model: str) -> int:
    if model not in _N_SHAPE:
        raise ValueError(f"unknown peak model {model!r}; choose from {MODELS}")
    return _N_SHAPE[model]


def _gaussian_1d(x: np.ndarray, x0: float, s: float) -> np.ndarray:
    return np.exp(-0.5 * ((x - x0) / max(s, 1e-6)) ** 2)


def _cell_gaussian_1d(x: np.ndarray, x0: float, s: float) -> np.ndarray:
    """A Gaussian integrated over each unit-wide cell centred on `x`, scaled to
    a maximum of about one. Each column of the folded image is one whole
    modulation, i.e. the first-dimension peak integrated over that interval, so
    a peak narrower than a column (common: alkane bands are 5 to 10 s wide
    against a 10 s modulation) is not a point-sampled Gaussian. For a width of
    several columns this is the same as `_gaussian_1d`."""
    s = max(s, 1e-6)
    z = np.asarray(x, dtype=float) - x0
    cell = 0.5 * (erf((z + 0.5) / (np.sqrt(2.0) * s)) - erf((z - 0.5) / (np.sqrt(2.0) * s)))
    return cell / max(float(erf(0.5 / (np.sqrt(2.0) * s))), 1e-12)


def _emg_1d(x: np.ndarray, x0: float, sigma: float, tau: float) -> np.ndarray:
    """Exponentially modified Gaussian (a Gaussian elution profile with
    an exponential tail, `tau` its time constant), in the numerically
    stable log-CDF form. Scaling depends only on the parameters, so rendering a
    subset of pixels gives exactly the same values as rendering the full image.
    The fitted linear coefficient absorbs the scale."""
    sigma, tau = max(sigma, 1e-6), max(tau, 1e-6)
    dx = np.asarray(x, dtype=float) - x0
    log_profile = (-np.log(tau) + sigma ** 2 / (2.0 * tau ** 2) - dx / tau +
                   log_ndtr(dx / sigma - sigma / tau))
    log_profile += np.log(sigma * np.sqrt(2.0 * np.pi))
    return np.exp(np.clip(log_profile, -745.0, 0.0))


def _shape(model: str, rr: np.ndarray, cc: np.ndarray, params: np.ndarray) -> np.ndarray:
    if model == "emg":
        r0, c0, sr, sc, tau = params
        return _gaussian_1d(rr, r0, sr) * _emg_1d(cc, c0, sc, tau)
    if model == "emg_t2":
        r0, c0, sr, sc, tau_r = params
        return _emg_1d(rr, r0, sr, tau_r) * _gaussian_1d(cc, c0, sc)
    if model == "emg2":
        r0, c0, sr, sc, tau_c, tau_r = params
        return _emg_1d(rr, r0, sr, tau_r) * _emg_1d(cc, c0, sc, tau_c)
    if model == "cell":
        r0, c0, sr, sc = params[:4]
        return _gaussian_1d(rr, r0, sr) * _cell_gaussian_1d(cc, c0, sc)
    if model == "emg_t2_cell":
        r0, c0, sr, sc, tau_r = params
        return _emg_1d(rr, r0, sr, tau_r) * _cell_gaussian_1d(cc, c0, sc)
    if model == "tilt":
        r0, c0, sr, sc, rho = params
        zr, zc = (rr - r0) / max(sr, 1e-6), (cc - c0) / max(sc, 1e-6)
        return np.exp(-(zr ** 2 + zc ** 2 - 2.0 * rho * zr * zc) / (2.0 * (1.0 - rho ** 2)))
    if model == "emg_t2_tilt":
        r0, c0, sr, sc, tau_r, rho = params
        # The t2 axis is sheared along t1 first (the apex drifts by rho * sr / sc rows per column).
        shifted = rr - rho * (sr / max(sc, 1e-6)) * (cc - c0)
        return _emg_1d(shifted, r0, sr, tau_r) * _gaussian_1d(cc, c0, sc)
    r0, c0, sr, sc = params[:4]
    return np.exp(-0.5 * ((rr - r0) / max(sr, 1e-6)) ** 2 - 0.5 * ((cc - c0) / max(sc, 1e-6)) ** 2)


def _shapes_from_theta(theta: np.ndarray, k: int, model: str, rr: np.ndarray, cc: np.ndarray) -> np.ndarray:
    p = _n_shape_params(model)
    mats = np.empty((rr.size, k))
    for j in range(k):
        mats[:, j] = _shape(model, rr, cc, theta[j * p:(j + 1) * p])
    return mats


def _solve_linear(X: np.ndarray, D: np.ndarray, non_negative: bool = True,
                  weights: np.ndarray | None = None,
                  n_nonnegative: int | None = None) -> np.ndarray:
    """Coefficients solving the same weighted objective as the outer fit.

    ``n_nonnegative`` constrains only the first columns.  This lets compound
    spectra remain physical while a fitted background follows signed residuals
    left by preprocessing.
    """
    k = X.shape[1]
    coef = np.empty((k, D.shape[1]))
    for ion in range(D.shape[1]):
        Xi, yi = X, D[:, ion]
        if weights is not None:
            wi = weights[:, ion]
            Xi, yi = X * wi[:, None], yi * wi
        if non_negative and n_nonnegative is not None and n_nonnegative < k:
            lower = np.r_[np.zeros(n_nonnegative),
                          np.full(k - n_nonnegative, -np.inf)]
            coef[:, ion] = lsq_linear(
                Xi, yi, bounds=(lower, np.full(k, np.inf)), method="bvls").x
        elif non_negative:
            coef[:, ion], _ = nnls(Xi, yi)
        else:
            coef[:, ion], *_ = np.linalg.lstsq(Xi, yi, rcond=None)
    return coef


def _initial_guess(cube: BlobCube, k: int, model: str, apex_hint=None) -> np.ndarray:
    if apex_hint and len(apex_hint) >= k:
        centres = [(a.row, a.col) for a in apex_hint[:k]]
    else:
        img = np.where(cube.mask, cube.data.sum(axis=2), -np.inf)
        order = np.argsort(img.ravel())[::-1]
        centres = []
        for flat in order:
            r, c = np.unravel_index(int(flat), img.shape)
            if all(abs(r - rr) + abs(c - cc) > 1 for rr, cc in centres):
                centres.append((float(r), float(c)))
            if len(centres) == k:
                break
        while len(centres) < k:
            centres.append(centres[-1] if centres else (cube.mask.shape[0] / 2.0, cube.mask.shape[1] / 2.0))
    # Start each width at what the data show: the full width at half maximum of
    # the summed-ion profile through the centre, along each axis. A fixed start
    # of one pixel is far too narrow for real 20 ms rows (peaks are several rows
    # wide) and the optimiser then never leaves it.
    img = np.where(cube.mask, np.clip(cube.data.sum(axis=2), 0, None), 0.0)

    def sigma_from(profile: np.ndarray, i: int) -> float:
        i = int(np.clip(round(i), 0, profile.size - 1))
        half = profile[i] / 2.0
        if half <= 0:
            return 1.0
        lo = hi = i
        while lo > 0 and profile[lo - 1] >= half:
            lo -= 1
        while hi < profile.size - 1 and profile[hi + 1] >= half:
            hi += 1
        return max((hi - lo + 1) / 2.355, 0.7)

    guesses = []
    for r0, c0 in centres:
        ri = int(np.clip(round(r0), 0, img.shape[0] - 1))
        ci = int(np.clip(round(c0), 0, img.shape[1] - 1))
        sr, sc = sigma_from(img[:, ci], ri), sigma_from(img[ri, :], ci)
        start = {"emg": [sc], "emg_t2": [sr], "emg2": [sc, sr], "tilt": [0.0], "emg_t2_tilt": [sr, 0.0],
                 "emg_t2_cell": [sr]}
        guesses += [r0, c0, sr, sc] + start.get(model, [])
    return np.array(guesses, dtype=float)


def _bounds(k: int, model: str, shape: tuple) -> tuple:
    p = _n_shape_params(model)
    lo, hi = [], []
    for _ in range(k):
        lo += [-2.0, -2.0, 0.05, 0.05]
        hi += [shape[0] + 1.0, shape[1] + 1.0, float(max(shape)), float(max(shape))]
        for kind in _EXTRAS[model]:
            lo += [-0.95 if kind == "rho" else 0.01]
            hi += [0.95 if kind == "rho" else float(max(shape) * 3)]
    return np.array(lo[:k * p]), np.array(hi[:k * p])


@dataclass
class OverlapComponent:
    row0: float
    col0: float
    sigma_row: float
    sigma_col: float
    tau: float | None
    t1: float
    t2: float
    volume: float
    spectrum: np.ndarray
    row0_se: float
    col0_se: float
    tau_row: float | None = None
    rho: float | None = None


@dataclass
class OverlapFit:
    k: int
    model: str
    components: list
    bic: float
    rss: float
    converged: bool
    weighted_rss: float = np.nan
    n_parameters: int = 0
    n_shape_parameters: int = 0
    background: np.ndarray = field(default_factory=lambda: np.array([]))
    background_surface: np.ndarray = field(default_factory=lambda: np.array([]))
    n_starts: int = 1
    restart_position_spread: float = np.nan
    profiles: np.ndarray = field(default_factory=lambda: np.array([]))


@dataclass(frozen=True)
class ComponentAssessment:
    component: int
    relative_to_largest: float
    confidence: str
    reason: str


def assess_components(fit: OverlapFit, low_confidence_ratio: float = 0.30
                      ) -> list[ComponentAssessment]:
    """Give fitted amounts a deliberately simple, auditable confidence label.

    Recovery on the supplied acquisition currently supports components at least
    30% as large as their neighbour under suitable backgrounds.  Smaller ones
    are retained, as requested, but labelled low confidence.  No component is
    called high confidence until identified standards validate the boundary.
    """
    if low_confidence_ratio <= 0 or low_confidence_ratio > 1:
        raise ValueError("low_confidence_ratio must be in (0, 1]")
    largest = max((component.volume for component in fit.components), default=0.0)
    out = []
    for j, component in enumerate(fit.components):
        relative = float(component.volume / largest) if largest > 0 else 0.0
        reasons = []
        if not fit.converged:
            reasons.append("fit did not fully converge")
        if relative < low_confidence_ratio:
            reasons.append(f"only {100*relative:.0f}% of the largest component")
        confidence = "low" if reasons else "provisional"
        reason = "; ".join(reasons) if reasons else "above the current recovery limit"
        out.append(ComponentAssessment(j, relative, confidence, reason))
    return out


def component_profiles(model: str, components: list[OverlapComponent],
                       shape: tuple[int, int], profiles: np.ndarray | None = None
                       ) -> np.ndarray:
    """Render normalised spatial profiles for fitted components."""
    if model == "template" and profiles is not None and np.asarray(profiles).size:
        return np.asarray(profiles, dtype=float)
    maps = np.zeros(shape + (len(components),), dtype=float)
    rr, cc = np.indices(shape, dtype=float)
    for j, component in enumerate(components):
        maps[:, :, j] = _shape(model, rr.ravel(), cc.ravel(),
                               shape_params(model, component)).reshape(shape)
    return maps


def shape_params(model: str, component) -> np.ndarray:
    """The parameter vector `_shape` takes, read from a fitted component's fields."""
    params = [component.row0, component.col0, component.sigma_row, component.sigma_col]
    if model == "emg":
        params.append(component.tau)
    elif model in ("emg_t2", "emg_t2_cell"):
        params.append(component.tau_row)
    elif model == "emg2":
        params.extend((component.tau, component.tau_row))
    elif model == "tilt":
        params.append(component.rho)
    elif model == "emg_t2_tilt":
        params.extend((component.tau_row, component.rho))
    return np.asarray(params, dtype=float)


def initialise_overlap(cube: BlobCube, k: int, model: str = "gaussian",
                       apex_hint=None) -> OverlapFit:
    """Build peak shapes from candidate positions and local half-height widths."""
    theta = _initial_guess(cube, k, model, apex_hint)
    p = _n_shape_params(model)
    components = []
    for index in range(k):
        params = theta[index * p:(index + 1) * p]
        row0, col0, sigma_row, sigma_col = params[:4]
        tau = float(params[4]) if model in ("emg", "emg2") else None
        tau_row = (float(params[4]) if model in ("emg_t2", "emg_t2_tilt", "emg_t2_cell")
                   else float(params[5]) if model == "emg2" else None)
        rho = (float(params[4]) if model == "tilt"
               else float(params[5]) if model == "emg_t2_tilt" else None)
        components.append(OverlapComponent(
            row0=float(row0), col0=float(col0), sigma_row=float(sigma_row),
            sigma_col=float(sigma_col), tau=tau, tau_row=tau_row, rho=rho,
            t1=float(cube.t1_of(col0)), t2=float(cube.t2_of(row0)), volume=0.0,
            spectrum=np.zeros(cube.data.shape[2]), row0_se=np.nan, col0_se=np.nan))
    return OverlapFit(k=k, model=model, components=components, bic=np.nan, rss=np.nan,
                      converged=True, weighted_rss=np.nan,
                      n_parameters=k * p, n_starts=0,
                      restart_position_spread=np.nan)


def fit_shapes_fast(cube: BlobCube, k: int, model: str = "gaussian",
                    apex_hint=None, background_order: int = 0,
                    max_nfev: int = 12, position_rows: float = 1.0,
                    position_cols: float = 0.6, max_sigma_row: float = 3.0,
                    max_sigma_col: float = 1.2,
                    shared_widths: bool = False) -> OverlapFit:
    """Fit local shapes with the same constrained objective as the final spectra.

    Spectra stay non-negative, while the local background can be signed after
    cleaning. Widths are independent unless explicitly shared. Projected
    derivatives include the change in the re-solved spectra; holding spectra
    fixed gave an inconsistent Jacobian and frequent false non-convergence.
    """
    rows, cols = np.nonzero(cube.mask)
    rr, cc = rows.astype(float), cols.astype(float)
    data = cube.data[rows, cols].astype(float)
    weights = 1.0 / np.sqrt(np.clip(data, 1.0, None))
    if background_order == 1:
        rn = rr / max(cube.mask.shape[0] - 1, 1)
        cn = cc / max(cube.mask.shape[1] - 1, 1)
        background = np.column_stack(((1-rn)*(1-cn), (1-rn)*cn,
                                       rn*(1-cn), rn*cn))
    elif background_order == 0:
        background = np.ones((rows.size, 1))
    else:
        raise ValueError("background_order must be 0 or 1")
    initial = _initial_guess(cube, k, model, apex_hint)
    lo, hi = _bounds(k, model, cube.mask.shape)
    p = _n_shape_params(model)
    initial_matrix = initial.reshape(k, p)
    lo_matrix, hi_matrix = lo.reshape(k, p), hi.reshape(k, p)
    lo_matrix[:, 0] = np.maximum(lo_matrix[:, 0], initial_matrix[:, 0] - position_rows)
    hi_matrix[:, 0] = np.minimum(hi_matrix[:, 0], initial_matrix[:, 0] + position_rows)
    lo_matrix[:, 1] = np.maximum(lo_matrix[:, 1], initial_matrix[:, 1] - position_cols)
    hi_matrix[:, 1] = np.minimum(hi_matrix[:, 1], initial_matrix[:, 1] + position_cols)
    hi_matrix[:, 2] = np.minimum(hi_matrix[:, 2], max_sigma_row)
    hi_matrix[:, 3] = np.minimum(hi_matrix[:, 3], max_sigma_col)
    initial = np.clip(initial, lo, hi)

    # Optional width sharing removes two parameters per additional Gaussian.
    # Independent widths are the default: enforcing equality biased recovery
    # when a narrow peak overlapped a broader one in the calibration windows.
    share = bool(shared_widths and model == "gaussian" and k > 1)
    if share:
        initial_matrix = initial.reshape(k, p)
        reduced_initial = np.r_[initial_matrix[:, :2].ravel(),
                                np.median(initial_matrix[:, 2], axis=0),
                                np.median(initial_matrix[:, 3], axis=0)]
        reduced_lo = np.r_[lo_matrix[:, :2].ravel(),
                           np.max(lo_matrix[:, 2]), np.max(lo_matrix[:, 3])]
        reduced_hi = np.r_[hi_matrix[:, :2].ravel(),
                           np.min(hi_matrix[:, 2]), np.min(hi_matrix[:, 3])]
        reduced_initial = np.clip(reduced_initial, reduced_lo, reduced_hi)

        def expand(theta):
            matrix = np.empty((k, p), dtype=float)
            matrix[:, :2] = theta[:2*k].reshape(k, 2)
            matrix[:, 2] = theta[-2]
            matrix[:, 3] = theta[-1]
            return matrix.ravel()
    else:
        reduced_initial, reduced_lo, reduced_hi = initial, lo, hi

        def expand(theta):
            return theta

    cached = {}

    def solve(theta):
        if 'theta' not in cached or not np.array_equal(cached['theta'], theta):
            peaks = _shapes_from_theta(expand(theta), k, model, rr, cc)
            design = np.column_stack((peaks, background))
            coefficients = _solve_linear(design, data, True, weights, n_nonnegative=k)
            cached.update(theta=theta.copy(), peaks=peaks, design=design,
                          coefficients=coefficients, residual=(design@coefficients-data)*weights)
        return cached

    def residual(theta):
        return solve(theta)['residual'].ravel()

    def jacobian(theta):
        solved = solve(theta)
        derivatives = _gaussian_shape_derivatives(
            expand(theta).reshape(k, p), solved['peaks'], rr, cc)
        return _projected_shape_jacobian(
            solved['design'], solved['coefficients'], weights,
            solved['residual'], derivatives, shared_widths=share)

    result = least_squares(residual, reduced_initial,
                           bounds=(reduced_lo, reduced_hi), max_nfev=max_nfev,
                           ftol=1e-5, xtol=1e-5, gtol=1e-5,
                           jac=jacobian if model in ('gaussian', 'tilt') else '2-point',
                           x_scale="jac")
    result_theta = expand(result.x)
    fitted = initialise_overlap(cube, k, model, apex_hint)
    for index, component in enumerate(fitted.components):
        params = result_theta[index * p:(index + 1) * p]
        component.row0, component.col0 = map(float, params[:2])
        component.sigma_row, component.sigma_col = map(float, params[2:4])
        component.t1 = float(cube.t1_of(component.col0))
        component.t2 = float(cube.t2_of(component.row0))
        if model in ("emg", "emg2"):
            component.tau = float(params[4])
        if model in ("emg_t2", "emg_t2_tilt", "emg_t2_cell"):
            component.tau_row = float(params[4])
        elif model == "emg2":
            component.tau_row = float(params[5])
        if model == "tilt":
            component.rho = float(params[4])
        elif model == "emg_t2_tilt":
            component.rho = float(params[5])
    fitted.converged = bool(result.success)
    fitted.n_shape_parameters = int(2 * k + 2 if share else k * p)
    fitted.n_starts = 1
    fitted.restart_position_spread = np.nan  # one start cannot establish stability
    return fitted


def _gaussian_shape_derivatives(parameters, peaks, rr, cc):
    """Profile derivatives for independent or correlated Gaussians."""
    derivatives = np.zeros((len(rr), *parameters.shape))
    for j, values in enumerate(parameters):
        r0, c0, sr, sc = values[:4]
        rho = values[4] if len(values) == 5 else 0.0
        zr, zc = (rr-r0)/sr, (cc-c0)/sc
        den = 1-rho*rho
        derivatives[:, j, 0] = peaks[:, j]*(zr-rho*zc)/(sr*den)
        derivatives[:, j, 1] = peaks[:, j]*(zc-rho*zr)/(sc*den)
        derivatives[:, j, 2] = derivatives[:, j, 0]*zr
        derivatives[:, j, 3] = derivatives[:, j, 1]*zc
        if len(values) == 5:
            derivatives[:, j, 4] = peaks[:, j]*(zr*zc*(1+rho*rho)-rho*(zr*zr+zc*zc))/den**2
    return derivatives


def _projected_shape_jacobian(design, coefficients, weights, residual, derivatives,
                               shared_widths=False):
    """Differentiate the constrained separable solve on its current active set.

    Both terms from differentiating the linear optimum are retained. The
    projection removes changes which the spectra/background can absorb; the
    second term accounts for a nonzero residual. Inactive spectral coefficients
    stay zero locally. This agrees with finite differences away from active-set
    transitions and avoids re-solving every ion for every shape parameter.
    """
    n_pixels, k, p = derivatives.shape
    n_parameters = 2*k+2 if shared_widths else k*p
    jac = np.zeros((n_pixels, coefficients.shape[1], n_parameters))
    for ion in range(coefficients.shape[1]):
        active_peaks = np.flatnonzero(coefficients[:k, ion] > 0)
        active = np.r_[active_peaks, np.arange(k, design.shape[1])]
        matrix = design[:, active]*weights[:, ion, None]
        inverse = np.linalg.pinv(matrix)
        direct = np.zeros((n_pixels, n_parameters))
        gradient = np.zeros((len(active), n_parameters))
        for index, j in enumerate(active_peaks):
            weighted = derivatives[:, j]*weights[:, ion, None]
            columns = [2*j, 2*j+1, 2*k, 2*k+1] if shared_widths else list(range(j*p, (j+1)*p))
            direct[:, columns] += weighted*coefficients[j, ion]
            gradient[index, columns] = weighted.T@residual[:, ion]
        jac[:, ion] = direct-matrix@(inverse@direct)-inverse.T@gradient
    return jac.reshape(-1, n_parameters)


def component_maps(fit: OverlapFit, cube: BlobCube) -> np.ndarray:
    """Return each fitted component's TIC contribution without hard ownership.

    The result has shape ``(rows, columns, components)``.  Contributions may be
    positive at the same pixel, which is the visual and quantitative replacement
    for assigning each pixel to one blob border.
    """
    profiles = component_profiles(fit.model, fit.components, cube.mask.shape, fit.profiles)
    maps = np.zeros_like(profiles)
    for j, component in enumerate(fit.components):
        maps[:, :, j] = profiles[:, :, j] * component.spectrum.sum()
    maps[~cube.mask] = 0.0
    return maps


def refit_spectra(fit: OverlapFit, cube: BlobCube,
                  include_background: bool = True, background_order: int = 1) -> OverlapFit:
    """Keep fitted shapes fixed and solve spectra/background on another ion cube.

    This is the inexpensive second half of automatic detection: diagnostic ions
    locate and shape the peaks, then a wider ion set supplies their final spectra
    without repeating the nonlinear optimisation.
    """
    if background_order not in (0, 1):
        raise ValueError("background_order must be 0 (constant) or 1 (bilinear)")
    rows, cols = np.nonzero(cube.mask)
    data = cube.data[rows, cols].astype(float)
    weights = 1.0 / np.sqrt(np.clip(data, 1.0, None))
    profiles = component_profiles(fit.model, fit.components, cube.mask.shape, fit.profiles)
    peak_design = profiles[rows, cols]
    rr, cc = rows.astype(float), cols.astype(float)
    if include_background and background_order == 1:
        rn = rr / max(cube.mask.shape[0] - 1, 1)
        cn = cc / max(cube.mask.shape[1] - 1, 1)
        background_basis = np.column_stack(((1-rn)*(1-cn), (1-rn)*cn,
                                             rn*(1-cn), rn*cn))
    elif include_background:
        background_basis = np.ones((rows.size, 1))
    else:
        background_basis = np.zeros((rows.size, 0))
    design = (np.column_stack((peak_design, background_basis))
              if include_background else peak_design)
    coefficients = _solve_linear(
        design, data, True, weights,
        n_nonnegative=fit.k if include_background else None)
    prediction = design @ coefficients
    rss = float(np.sum((prediction - data) ** 2))
    weighted_rss = float(np.sum(((prediction - data) * weights) ** 2))
    components = [replace(
        component, spectrum=coefficients[index].copy(),
        volume=float(profiles[:, :, index].sum() * coefficients[index].sum()))
        for index, component in enumerate(fit.components)]
    n_background = background_basis.shape[1]
    n_shape_parameters = (fit.n_shape_parameters or
                          (2 * fit.k if fit.model == 'template'
                           else fit.k * _n_shape_params(fit.model)))
    n_parameters = int(n_shape_parameters
                       + (fit.k + n_background) * data.shape[1])
    n_obs = data.size
    bic = (n_obs * np.log(max(weighted_rss / max(n_obs, 1), 1.0))
           + n_parameters * np.log(max(n_obs, 1)))
    background_coefficients = coefficients[fit.k:]
    background_surface = np.zeros(cube.data.shape, dtype=float)
    if n_background:
        background_surface[rows, cols] = background_basis @ background_coefficients
    return OverlapFit(
        k=fit.k, model=fit.model, components=components, bic=float(bic), rss=rss,
        converged=fit.converged, weighted_rss=weighted_rss,
        n_parameters=n_parameters,
        n_shape_parameters=int(n_shape_parameters),
        background=(background_coefficients[0].copy() if n_background == 1
                    else background_coefficients.copy()),
        background_surface=background_surface if n_background else np.array([]),
        n_starts=fit.n_starts, restart_position_spread=fit.restart_position_spread,
        profiles=profiles if fit.model == "template" else np.array([]))


def predict_overlap(fit: OverlapFit, cube: BlobCube) -> np.ndarray:
    """Reconstruct every pixel, including positions held out during fitting."""
    rr, cc = np.indices(cube.mask.shape, dtype=float)
    prediction = np.zeros(cube.data.shape, dtype=float)
    for j, component in enumerate(fit.components):
        if fit.model == "template" and fit.profiles.size:
            prediction += fit.profiles[:, :, j, None] * component.spectrum[None, None, :]
            continue
        params = [component.row0, component.col0,
                  component.sigma_row, component.sigma_col]
        if fit.model == "emg":
            params.append(component.tau)
        elif fit.model in ("emg_t2", "emg_t2_cell"):
            params.append(component.tau_row)
        elif fit.model == "emg2":
            params.extend((component.tau, component.tau_row))
        elif fit.model == "tilt":
            params.append(component.rho)
        elif fit.model == "emg_t2_tilt":
            params.extend((component.tau_row, component.rho))
        shape = _shape(fit.model, rr.ravel(), cc.ravel(), np.asarray(params)).reshape(
            cube.mask.shape)
        prediction += shape[:, :, None] * component.spectrum[None, None, :]
    if fit.background.size:
        if fit.background.ndim == 1:
            prediction += fit.background[None, None, :]
        elif fit.background.shape[0] == 4:
            rn = rr / max(cube.mask.shape[0] - 1, 1)
            cn = cc / max(cube.mask.shape[1] - 1, 1)
            basis = np.stack(((1-rn)*(1-cn), (1-rn)*cn, rn*(1-cn), rn*cn), axis=2)
            prediction += np.einsum("rqb,bi->rqi", basis, fit.background)
        else:
            raise ValueError("unknown fitted background basis")
    return prediction


def _render_template(template: np.ndarray, shape: tuple[int, int],
                     row: float, col: float) -> np.ndarray:
    """Move a measured peak template to a sub-pixel centre on ``shape``."""
    template = np.maximum(np.asarray(template, dtype=float), 0.0)
    if template.ndim != 2 or template.sum() <= 0:
        raise ValueError("template must be a non-empty two-dimensional positive array")
    template = template / template.sum()
    centre = ndi.center_of_mass(template)
    rr, cc = np.indices(shape, dtype=float)
    moved = ndi.map_coordinates(
        template, [rr - row + centre[0], cc - col + centre[1]],
        order=1, mode="constant", cval=0.0, prefilter=False)
    total = moved.sum()
    return moved / total if total > 0 else moved


def fit_template_overlap(cube: BlobCube, template: np.ndarray, k: int,
                         apex_hint=None, include_background: bool = False,
                         background_order: int = 0, n_starts: int = 1,
                         seed: int = 0) -> OverlapFit:
    """Fit shifted copies of a measured chromatographic profile.

    Only each copy's row and column are nonlinear.  Its non-negative mass
    spectrum and an optional constant or bilinear background are solved with
    the same count-weighted variable-projection objective as ``fit_overlap``.
    The template is normalised to unit area, so each component's reported
    volume is the sum of its fitted spectrum.
    """
    if k < 1:
        raise ValueError("k must be at least one")
    if background_order not in (0, 1):
        raise ValueError("background_order must be 0 (constant) or 1 (bilinear)")
    rows, cols = np.nonzero(cube.mask)
    D = cube.data[rows, cols, :].astype(float)
    weights = 1.0 / np.sqrt(np.clip(D, 1.0, None))
    initial = _initial_guess(cube, k, "gaussian", apex_hint).reshape(k, 4)[:, :2].ravel()
    lo = np.tile([-2.0, -2.0], k)
    hi = np.tile([cube.mask.shape[0] + 1.0, cube.mask.shape[1] + 1.0], k)
    initial = np.clip(initial, lo, hi)

    rr = rows.astype(float) / max(cube.mask.shape[0] - 1, 1)
    cc = cols.astype(float) / max(cube.mask.shape[1] - 1, 1)
    if include_background and background_order == 1:
        background_basis = np.column_stack(((1-rr)*(1-cc), (1-rr)*cc,
                                             rr*(1-cc), rr*cc))
    elif include_background:
        background_basis = np.ones((rows.size, 1))
    else:
        background_basis = np.zeros((rows.size, 0))

    def rendered(theta):
        return np.stack([_render_template(template, cube.mask.shape,
                                          theta[2*j], theta[2*j+1])
                         for j in range(k)], axis=2)

    def design(theta):
        peak_profiles = rendered(theta)[rows, cols, :]
        return (np.column_stack((peak_profiles, background_basis))
                if include_background else peak_profiles)

    def residuals(theta):
        X = design(theta)
        coef = _solve_linear(X, D, True, weights,
                             n_nonnegative=k if include_background else None)
        prediction = X @ coef
        if not np.all(np.isfinite(prediction)):
            return np.full(D.size, 1e20)
        return ((prediction - D) * weights).ravel()

    rng = np.random.default_rng(seed)
    starts = [initial]
    for _ in range(1, max(int(n_starts), 1)):
        starts.append(np.clip(initial + rng.normal(0.0, 0.75, initial.size), lo, hi))
    results = [least_squares(residuals, start, bounds=(lo, hi)) for start in starts]
    result = min(results, key=lambda item: float(np.sum(item.fun ** 2)))
    profiles = rendered(result.x)
    X = design(result.x)
    coef = _solve_linear(X, D, True, weights,
                         n_nonnegative=k if include_background else None)
    prediction = X @ coef
    rss = float(np.sum((prediction - D) ** 2))
    weighted_rss = float(np.sum(((prediction - D) * weights) ** 2))

    dof = max(result.fun.size - result.x.size, 1)
    variance = float(np.sum(result.fun ** 2) / dof)
    try:
        covariance = variance * np.linalg.inv(result.jac.T @ result.jac)
        errors = np.sqrt(np.clip(np.diag(covariance), 0, None))
    except np.linalg.LinAlgError:
        errors = np.full(result.x.size, np.nan)
    components = []
    for j in range(k):
        row, col = result.x[2*j:2*j+2]
        components.append(OverlapComponent(
            row0=float(row), col0=float(col), sigma_row=np.nan, sigma_col=np.nan,
            tau=None, t1=float(cube.t1_of(col)), t2=float(cube.t2_of(row)),
            volume=float(coef[j].sum()), spectrum=coef[j].copy(),
            row0_se=float(errors[2*j]), col0_se=float(errors[2*j+1])))
    n_background = background_basis.shape[1]
    n_parameters = int(result.x.size + (k + n_background) * D.shape[1])
    n_obs = D.size
    bic = (n_obs * np.log(max(weighted_rss / n_obs, 1.0)) +
           n_parameters * np.log(n_obs))
    background_coefficients = coef[k:]
    background_surface = np.zeros(cube.data.shape, dtype=float)
    if n_background:
        background_surface[rows, cols] = background_basis @ background_coefficients
    return OverlapFit(
        k=k, model="template", components=components, bic=float(bic), rss=rss,
        converged=bool(result.success), weighted_rss=weighted_rss,
        n_parameters=n_parameters,
        background=(background_coefficients[0].copy() if n_background == 1
                    else background_coefficients.copy()),
        background_surface=background_surface if n_background else np.array([]),
        n_starts=len(starts), profiles=profiles)


@dataclass
class OverlapCV:
    k: int
    model: str
    fold_rmse: np.ndarray
    mean_rmse: float
    converged: np.ndarray


def cross_validate_overlap(cube: BlobCube, k: int, folds: int = 4,
                           block_rows: int = 3, **fit_kw) -> OverlapCV:
    """Fit on spatial blocks and score noise-weighted error on held-out pixels."""
    if folds < 2:
        raise ValueError("cross-validation needs at least two folds")
    rows, cols = np.indices(cube.mask.shape)
    groups = ((rows // max(int(block_rows), 1)) + cols) % folds
    scores, converged = [], []
    for held_out in range(folds):
        test = cube.mask & (groups == held_out)
        train = cube.mask & ~test
        if not test.any() or train.sum() < max(k * _n_shape_params(fit_kw.get("model", "gaussian")), 1):
            continue
        training_cube = BlobCube(cube.data, train, cube.ions, cube.row0, cube.col0,
                                 cube.pitch, cube.period, cube.t1_start)
        fit = fit_overlap(training_cube, k, **fit_kw)
        prediction = predict_overlap(fit, cube)
        observed = cube.data[test]
        weights = 1.0 / np.sqrt(np.clip(observed, 1.0, None))
        scores.append(float(np.sqrt(np.mean(((prediction[test] - observed) * weights) ** 2))))
        converged.append(fit.converged)
    if not scores:
        raise ValueError("cube is too small for the requested spatial folds")
    values = np.asarray(scores)
    return OverlapCV(k=k, model=fit_kw.get("model", "gaussian"), fold_rmse=values,
                     mean_rmse=float(values.mean()), converged=np.asarray(converged, dtype=bool))


def fit_overlap(cube: BlobCube, k: int, model: str = "gaussian", non_negative: bool = True,
                apex_hint=None, include_background: bool = False, n_starts: int = 1,
                seed: int = 0, background_order: int = 0,
                tolerance: float | None = None, max_nfev: int | None = None) -> OverlapFit:
    """A joint fit of k two-dimensional peak shapes times mass spectra
    to the cube (see the module docstring: variable projection, only
    the shape parameters are the nonlinear unknowns). Weights
    approximate Poisson counting noise (`1/sqrt(count)`, floored at 1
    count so an empty pixel doesn't get an infinite weight)."""
    mask = cube.mask
    rows, cols = np.nonzero(mask)
    rr, cc = rows.astype(float), cols.astype(float)
    D = cube.data[rows, cols, :].astype(float)
    weights = 1.0 / np.sqrt(np.clip(D, 1.0, None))

    theta0 = _initial_guess(cube, k, model, apex_hint)
    lo, hi = _bounds(k, model, mask.shape)
    theta0 = np.clip(theta0, lo, hi)

    if background_order not in (0, 1):
        raise ValueError("background_order must be 0 (constant) or 1 (bilinear)")
    if include_background and background_order == 1:
        rn = rr / max(cube.mask.shape[0] - 1, 1)
        cn = cc / max(cube.mask.shape[1] - 1, 1)
        background_basis = np.column_stack(((1-rn)*(1-cn), (1-rn)*cn,
                                             rn*(1-cn), rn*cn))
    elif include_background:
        background_basis = np.ones((rr.size, 1))
    else:
        background_basis = np.zeros((rr.size, 0))

    def design(theta):
        X = _shapes_from_theta(theta, k, model, rr, cc)
        return np.column_stack([X, background_basis]) if include_background else X

    def residuals(theta):
        X = design(theta)
        coef = _solve_linear(X, D, non_negative, weights,
                             n_nonnegative=k if include_background and non_negative else None)
        with np.errstate(over="ignore", invalid="ignore"):
            prediction = X @ coef
        if not np.all(np.isfinite(prediction)):
            return np.full(D.size, 1e20)
        return ((prediction - D) * weights).ravel()

    rng = np.random.default_rng(seed)
    starts = [theta0]
    p_shape = _n_shape_params(model)
    for _ in range(1, max(int(n_starts), 1)):
        trial = theta0.copy().reshape(k, p_shape)
        trial[:, :2] += rng.normal(0.0, 1.0, size=(k, 2))
        trial[:, 2:4] *= np.exp(rng.normal(0.0, 0.25, size=(k, 2)))
        for e, kind in enumerate(_EXTRAS[model]):
            if kind == "tau":
                trial[:, 4 + e] *= np.exp(rng.normal(0.0, 0.25, size=k))
            else:
                trial[:, 4 + e] += rng.normal(0.0, 0.15, size=k)
        starts.append(np.clip(trial.ravel(), lo, hi))
    solver_options = {}
    if tolerance is not None:
        solver_options.update(ftol=float(tolerance), xtol=float(tolerance),
                              gtol=float(tolerance))
    if max_nfev is not None:
        solver_options["max_nfev"] = int(max_nfev)
    results = [least_squares(residuals, start, bounds=(lo, hi), **solver_options)
               for start in starts]
    result = min(results, key=lambda r: float(np.sum(r.fun ** 2)))
    theta = result.x
    best_positions = theta.reshape(k, p_shape)[:, :2]
    aligned_positions = []
    for candidate in results:
        positions = candidate.x.reshape(k, p_shape)[:, :2]
        distance = np.linalg.norm(best_positions[:, None, :] - positions[None, :, :], axis=2)
        best_index, candidate_index = linear_sum_assignment(distance)
        aligned = np.empty_like(positions)
        aligned[best_index] = positions[candidate_index]
        aligned_positions.append(aligned)
    position_stack = np.stack(aligned_positions)
    restart_position_spread = float(np.max(np.sqrt(np.var(position_stack, axis=0).sum(axis=1))))
    X = design(theta)
    coef = _solve_linear(X, D, non_negative, weights,
                         n_nonnegative=k if include_background and non_negative else None)
    prediction = X @ coef
    if not np.all(np.isfinite(prediction)):
        raise RuntimeError("overlap fit produced a non-finite prediction")
    rss = float(np.sum((prediction - D) ** 2))
    weighted_rss = float(np.sum(((prediction - D) * weights) ** 2))

    J = result.jac
    dof = max(result.fun.size - theta.size, 1)
    s_sq = float(np.sum(result.fun ** 2) / dof)
    try:
        cov = s_sq * np.linalg.inv(J.T @ J)
        se = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        se = np.full(theta.size, np.nan)

    p = _n_shape_params(model)
    peak_shapes = X[:, :k]
    components = []
    for j in range(k):
        params, errs = theta[j * p:(j + 1) * p], se[j * p:(j + 1) * p]
        r0, c0, sr, sc = params[:4]
        tau = float(params[4]) if model in ("emg", "emg2") else None          # along t1
        tau_row = (float(params[4]) if model in ("emg_t2", "emg_t2_tilt", "emg_t2_cell") else
                   float(params[5]) if model == "emg2" else None)              # along t2
        rho = (float(params[4]) if model == "tilt" else
               float(params[5]) if model == "emg_t2_tilt" else None)
        components.append(OverlapComponent(
            row0=float(r0), col0=float(c0), sigma_row=float(sr), sigma_col=float(sc), tau=tau,
            tau_row=tau_row, rho=rho,
            t1=float(cube.t1_of(c0)), t2=float(cube.t2_of(r0)),
            volume=float(peak_shapes[:, j].sum() * coef[j].sum()), spectrum=coef[j].copy(),
            row0_se=float(errs[0]), col0_se=float(errs[1])))
    n_obs = D.size
    # Spectra are estimated parameters even though variable projection solves
    # them analytically.  Omitting them made extra components artificially cheap.
    n_background = background_basis.shape[1]
    n_linear = k + n_background
    n_parameters = int(theta.size + n_linear * D.shape[1])
    # The weights assert roughly unit residual variance.  Once a model fits
    # below that noise scale, rewarding ever smaller numerical residuals lets
    # redundant components fit interpolation error in otherwise noiseless test
    # data.  Cap the likelihood improvement at the stated noise floor.
    bic = (n_obs * np.log(max(weighted_rss / n_obs, 1.0)) +
           n_parameters * np.log(n_obs))
    background_coefficients = coef[k:]
    background_surface = np.zeros(cube.data.shape, dtype=float)
    if n_background:
        background_surface[rows, cols] = background_basis @ background_coefficients
    return OverlapFit(k=k, model=model, components=components, bic=float(bic), rss=rss,
                      converged=bool(result.success), weighted_rss=weighted_rss,
                      n_parameters=n_parameters,
                      background=(background_coefficients[0].copy() if n_background == 1
                                  else background_coefficients.copy()),
                      background_surface=background_surface if n_background else np.array([]),
                      n_starts=len(starts), restart_position_spread=restart_position_spread)


def choose_k(cube: BlobCube, max_k: int = 4, model: str = "gaussian", apex_hint=None,
             include_background: bool = False, n_starts: int = 1, seed: int = 0,
             background_order: int = 0) -> OverlapFit:
    """Fit k = 1..max_k (never more than the footprint has pixels to
    support) and keep whichever minimises BIC."""
    n_pixels = int(np.count_nonzero(cube.mask))
    max_k = max(1, min(max_k, n_pixels // max(_n_shape_params(model), 1)))
    fits = [fit_overlap(cube, k, model, apex_hint=apex_hint,
                        include_background=include_background, n_starts=n_starts,
                        seed=seed + k, background_order=background_order)
            for k in range(1, max_k + 1)]
    return min(fits, key=lambda f: f.bic)


# ------------------------------------------------------------- 4. MCR-ALS
def _unimodal_1d(v: np.ndarray) -> np.ndarray:
    """Nearest non-negative sequence sharing v's own peak location that
    never rises again after it starts falling, found by capping each
    value at its neighbour closer to the peak. A cheap, well-defined
    approximation to isotonic regression, not the exact least-squares
    projection -- see the module docstring."""
    v = np.clip(np.asarray(v, dtype=float), 0, None)
    if v.size <= 2:
        return v
    peak = int(np.argmax(v))
    out = v.copy()
    for i in range(peak - 1, -1, -1):
        out[i] = min(out[i], out[i + 1])
    for i in range(peak + 1, out.size):
        out[i] = min(out[i], out[i - 1])
    return out


def _enforce_unimodal(component_map: np.ndarray, iterations: int = 2) -> np.ndarray:
    m = component_map.copy()
    for _ in range(iterations):
        for r in range(m.shape[0]):
            m[r, :] = _unimodal_1d(m[r, :])
        for c in range(m.shape[1]):
            m[:, c] = _unimodal_1d(m[:, c])
    return m


@dataclass
class MCRResult:
    concentrations: np.ndarray   # (rows, cols, k), zero outside the footprint
    spectra: np.ndarray           # (k, n_ion)
    ions: np.ndarray
    n_iter: int
    converged: bool
    resid: float


def mcr_als(cube: BlobCube, k: int, max_iter: int = 200, tol: float = 1e-6,
           unimodal: bool = True, seed: int = 0) -> MCRResult:
    """Multivariate curve resolution by alternating non-negative least
    squares: recover k free-form elution profiles and k free-form
    spectra directly from the cube, no peak shape assumed (see the
    module docstring)."""
    mask = cube.mask
    shape = mask.shape
    idx = np.flatnonzero(mask.ravel())
    D = cube.data.reshape(-1, cube.data.shape[-1])[idx].astype(float)
    n_pix, n_ion = D.shape
    if k < 1 or k > n_pix:
        raise ValueError("need at least as many footprint pixels as components")

    rng = np.random.default_rng(seed)
    picks = rng.choice(n_pix, size=k, replace=False)
    S = np.clip(D[picks].astype(float).copy(), 1e-9, None)

    prev_ssq, converged, n_iter, C = np.inf, False, 0, np.zeros((n_pix, k))
    for it in range(max_iter):
        n_iter = it + 1
        A = S.T
        for p in range(n_pix):
            C[p], _ = nnls(A, D[p])
        if unimodal:
            full = np.zeros(shape + (k,))
            full.reshape(-1, k)[idx] = C
            for j in range(k):
                full[:, :, j] = _enforce_unimodal(full[:, :, j])
            C = full.reshape(-1, k)[idx]
        S = np.zeros((k, n_ion))
        for ion in range(n_ion):
            S[:, ion], _ = nnls(C, D[:, ion])
        resid = D - C @ S
        ssq = float(np.sum(resid ** 2))
        if abs(prev_ssq - ssq) < tol * max(prev_ssq, 1.0):
            converged = True
            prev_ssq = ssq
            break
        prev_ssq = ssq

    full_c = np.zeros(shape + (k,))
    full_c.reshape(-1, k)[idx] = C
    return MCRResult(concentrations=full_c, spectra=S, ions=cube.ions, n_iter=n_iter,
                     converged=converged, resid=float(np.sqrt(prev_ssq / max(D.size, 1))))


@dataclass
class MCRStability:
    """Repeated-start evidence; low residual alone does not make a solution unique."""
    best: MCRResult
    residuals: np.ndarray
    explained_variance: np.ndarray
    mean_spectral_similarity: float
    minimum_spectral_similarity: float
    converged: np.ndarray


def mcr_restarts(cube: BlobCube, k: int, n_starts: int = 5, seed: int = 0,
                 **mcr_kw) -> MCRStability:
    """Run MCR-ALS repeatedly and match components by spectral cosine.

    Similarities compare every restart with the first solution after optimal
    one-to-one component matching. Values near one support a reproducible
    factorisation; they do not establish identity or component count.
    """
    if n_starts < 2:
        raise ValueError("MCR stability needs at least two starts")
    runs = [mcr_als(cube, k, seed=seed + i, **mcr_kw) for i in range(n_starts)]
    reference = runs[0].spectra / np.maximum(
        np.linalg.norm(runs[0].spectra, axis=1)[:, None], 1e-12)
    similarities = []
    for result in runs[1:]:
        spectra = result.spectra / np.maximum(
            np.linalg.norm(result.spectra, axis=1)[:, None], 1e-12)
        cosine = reference @ spectra.T
        rows, cols = linear_sum_assignment(-cosine)
        similarities.append(float(np.mean(cosine[rows, cols])))
    residuals = np.asarray([result.resid for result in runs])
    data = cube.data[cube.mask]
    data_rms = float(np.sqrt(np.mean(data ** 2))) if data.size else 0.0
    explained = 1.0 - (residuals / max(data_rms, 1e-12)) ** 2
    return MCRStability(
        best=runs[int(np.argmin(residuals))], residuals=residuals,
        explained_variance=explained,
        mean_spectral_similarity=float(np.mean(similarities)),
        minimum_spectral_similarity=float(np.min(similarities)),
        converged=np.asarray([result.converged for result in runs], dtype=bool))


# ------------------------------------------------------------- 5. PARAFAC
@dataclass
class PARAFACResult:
    factors: tuple      # (A, B, C): (runs, r), (elution positions, r), (ions, r)
    r: int
    n_iter: int
    converged: bool
    resid: float


def _khatri_rao(U: np.ndarray, V: np.ndarray) -> np.ndarray:
    """Column-wise Khatri-Rao product: column j is kron(U[:, j], V[:, j])."""
    return np.stack([np.kron(U[:, j], V[:, j]) for j in range(U.shape[1])], axis=1)


def _mode_update(Xn: np.ndarray, KR: np.ndarray, non_negative: bool) -> np.ndarray:
    n, r = Xn.shape[0], KR.shape[1]
    M = np.empty((n, r))
    for i in range(n):
        if non_negative:
            M[i], _ = nnls(KR, Xn[i])
        else:
            M[i], *_ = np.linalg.lstsq(KR, Xn[i], rcond=None)
    return M


def parafac(tensor: np.ndarray, r: int, max_iter: int = 200, tol: float = 1e-8,
           non_negative: bool = True, seed: int = 0) -> PARAFACResult:
    """A basic trilinear (CP) decomposition by alternating least
    squares: `tensor[i, j, k] ~= sum_r A[i, r] * B[j, r] * C[k, r]`.
    Non-negativity is enforced by solving each mode's update with
    `nnls` directly (row by row), not a plain least squares followed by
    clipping."""
    I, J, K = tensor.shape
    rng = np.random.default_rng(seed)
    A, B, C = rng.random((I, r)), rng.random((J, r)), rng.random((K, r))
    prev, converged, n_iter = np.inf, False, 0
    for it in range(max_iter):
        n_iter = it + 1
        A = _mode_update(tensor.reshape(I, J * K), _khatri_rao(B, C), non_negative)
        B = _mode_update(np.transpose(tensor, (1, 0, 2)).reshape(J, I * K), _khatri_rao(A, C), non_negative)
        C = _mode_update(np.transpose(tensor, (2, 0, 1)).reshape(K, I * J), _khatri_rao(A, B), non_negative)
        recon = np.einsum('ir,jr,kr->ijk', A, B, C)
        ssq = float(np.sum((tensor - recon) ** 2))
        if abs(prev - ssq) < tol * max(prev, 1.0):
            converged = True
            prev = ssq
            break
        prev = ssq
    return PARAFACResult(factors=(A, B, C), r=r, n_iter=n_iter, converged=converged,
                         resid=float(np.sqrt(prev / max(tensor.size, 1))))


def core_consistency(tensor: np.ndarray, factors: tuple) -> float:
    """Bro and Kiers' (2003) core consistency diagnostic: how closely
    the best Tucker core the factors imply resembles the ideal
    superdiagonal one a genuine trilinear model would have. Near 100
    says the chosen number of components is trustworthy; low or
    negative says it's over-factored."""
    A, B, C = factors
    r = A.shape[1]
    Ainv, Binv, Cinv = (np.linalg.pinv(M) for M in (A, B, C))
    G = np.einsum('ijk,ri->rjk', tensor, Ainv)
    G = np.einsum('rjk,sj->rsk', G, Binv)
    G = np.einsum('rsk,tk->rst', G, Cinv)
    ideal = np.zeros((r, r, r))
    for i in range(r):
        ideal[i, i, i] = 1.0
    return float(100.0 * (1.0 - np.sum((G - ideal) ** 2) / np.sum(ideal ** 2)))


def split_half_stability(tensor: np.ndarray, r: int, axis: int = 0, seed: int = 0, **parafac_kw) -> float:
    """Run PARAFAC on two random halves of the tensor split along
    `axis` (the runs mode, by default) and report the mean correlation
    between the two halves' recovered profiles on the two modes that
    survive the split. Components are matched between halves by best
    correlation on the first surviving mode (a greedy match, not a full
    assignment) before averaging."""
    n = tensor.shape[axis]
    if n < 4:
        raise ValueError("need at least 4 slices along the split axis for a meaningful split-half test")
    rng = np.random.default_rng(seed)
    idx = rng.permutation(n)
    half = n // 2
    fit_a = parafac(np.take(tensor, idx[:half], axis=axis), r, seed=seed, **parafac_kw)
    fit_b = parafac(np.take(tensor, idx[half:], axis=axis), r, seed=seed, **parafac_kw)
    surviving = [m for m in range(3) if m != axis]

    Fa0, Fb0 = fit_a.factors[surviving[0]], fit_b.factors[surviving[0]]
    corr = np.abs(np.corrcoef(Fa0.T, Fb0.T)[:r, r:])
    order = np.argmax(corr, axis=1)

    scores = []
    for m in surviving:
        Fa, Fb = fit_a.factors[m], fit_b.factors[m]
        for j, matched in enumerate(order):
            c = np.corrcoef(Fa[:, j], Fb[:, matched])[0, 1]
            scores.append(abs(c) if np.isfinite(c) else 0.0)
    return float(np.mean(scores)) if scores else 0.0


def build_run_tensor(prepped: dict, ref_t1: float, ref_t2: float, mz: list,
                     t1_window_s: float = 12.0, t2_window_s: float = 0.3,
                     n_t1: int = 9, n_t2: int = 9) -> tuple:
    """A small (run, elution position, ion) tensor around one matched
    feature (e.g. a `batch.Feature`'s own `ref_t1`/`ref_t2`), for
    `parafac`. Every run in `prepped` (a name -> object with `.clean`
    and `.fold`, such as `batch.PrepResult`) contributes an
    `n_t1` x `n_t2` grid of each ion's intensity, read directly from
    that run's own folded ion image by bilinear interpolation
    (`scipy.ndimage.map_coordinates`) at the same absolute (t1, t2)
    window -- this assumes retention times are already comparable
    across runs (true of a batch that's already phase-corrected to one
    method), not a full geometric alignment the way `align.py` would
    do for two runs that genuinely drift apart.

    Returns (tensor, run_names, ions) with `tensor.shape ==
    (len(run_names), n_t1 * n_t2, len(ions))`.
    """
    from scipy.ndimage import map_coordinates

    ions = np.asarray(mz, dtype=float)
    t1_grid = np.linspace(ref_t1 - t1_window_s, ref_t1 + t1_window_s, n_t1)
    t2_grid = np.linspace(max(ref_t2 - t2_window_s, 0.0), ref_t2 + t2_window_s, n_t2)
    tt1, tt2 = np.meshgrid(t1_grid, t2_grid, indexing="ij")

    run_names = list(prepped)
    tensor = np.zeros((len(run_names), n_t1 * n_t2, len(ions)))
    for i, name in enumerate(run_names):
        p = prepped[name]
        fold_p = p.fold
        col = (tt1 - fold_p.t1[0]) / fold_p.spec.period
        row = tt2 / fold_p.row_pitch
        for k, ion in enumerate(ions):
            img = fold_p.ion_image(p.clean, ion)
            vals = map_coordinates(img, [row.ravel(), col.ravel()], order=1, mode="constant", cval=0.0)
            tensor[i, :, k] = np.clip(vals, 0, None)
    return tensor, run_names, ions
