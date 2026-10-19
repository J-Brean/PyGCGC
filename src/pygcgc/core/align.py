"""Aligning one run onto another.

The idea in plain terms
------------------------
Run to run, retention drifts a little: a few seconds in t1, a few
hundredths of a second in t2. To compare two chromatograms pixel for
pixel (blink, difference, overlay) or tile several in a grid, one run's
image is resampled onto another's (the reference's) time grid.

You give this a handful of anchors: the same peak's true position in the
reference run (`ref_t1`, `ref_t2`) and in the run being aligned (`run_t1`,
`run_t2`). A transform is fitted through them and used two ways: forward,
run-to-reference, to judge the fit; and inverse, reference-to-run, to
resample the run's image onto the reference grid (for every reference
pixel, look up where that point came from in the run, and read it there
by bilinear interpolation).

How many anchors buys what
---------------------------
- **One anchor**: a plain shift in t1 and in t2 (the median shift, if
  more than one anchor is given at this tier).
- **Three or more**: a straight-line stretch and shift in t1
  (`ref_t1 = a1 + b1 * run_t1`), plus a t2 offset allowed to vary
  linearly with t1 (`ref_t2 = run_t2 + a2 + b2 * run_t1`) — useful when
  the second dimension drifts differently early and late in the run.
  Both fits are Huber-weighted, so one mis-clicked anchor doesn't wreck
  the rest.
- **Six or more**: a smooth map instead of a straight line — a monotone
  cubic (PCHIP) through the t1 anchors, and a smoothing spline for the t2
  offset against t1. This follows a retention curve that isn't linear
  without chasing every last wiggle.

A transform that comes out with its t1 anchors in the wrong order (the
run's peaks arrive in a different sequence than the reference's) cannot
be a genuine time correction; rather than raise, `monotonic` is set to
`False` and `warning` explains why, so the caller can decide what to do
with a transform that isn't to be trusted.

`fit_alignment(points, kind="auto")` picks the simplest tier that already
explains the anchors well, only reaching for a more flexible one if the
simpler one leaves the anchors clearly off; ties go to the simpler tier.
`cross_validate` refits with each anchor held out in turn and predicts
it, an honest error in seconds rather than the always-optimistic in-
sample residual — used by `AlignmentTransform.summary()`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy import interpolate as ip
from scipy import ndimage as ndi

from .fold import Fold
from .phase import Anchor, locate_anchor
from .run import Run

FIT_MIN = {"shift": 1, "affine": 3, "smooth": 6}
AUTO_MIN = {"shift": 1, "affine": 4, "smooth": 7}
TIER_ORDER = ("shift", "affine", "smooth")
AUTO_T1_TOL = 1.0     # seconds: an in-sample fit this good (or better) is "good enough" for auto
AUTO_T2_TOL = 0.02
DEFAULT_EXTEND_FACTOR = 1.15


class AlignmentError(ValueError):
    pass


@dataclass(frozen=True)
class AnchorPoint:
    """The same feature's true position in the reference run and in the
    run being aligned. `ion` is the m/z it was found or clicked on (0 if
    unknown); `name` labels it (a compound, or a plain "P1" for a
    manually clicked point)."""
    name: str
    ref_t1: float
    ref_t2: float
    run_t1: float
    run_t2: float
    ion: int = 0


def refine_anchor_position(image: np.ndarray, fold: Fold, t1: float, t2: float,
                           columns: int = 1, t2_radius: float = .15) -> tuple[float, float] | None:
    """Find the brightest nearby pixel, staying inside a small search box.

    Search at most ``columns`` modulations either side and ``t2_radius``
    seconds vertically. At a modulation boundary, use the neighbouring
    column's rows, just as an extended chromatogram does. Return coordinates
    in the clicked peak's unwrapped frame, so a boundary crossing stays a
    small adjustment. Ties favour the nearest pixel. Empty, flat or entirely
    non-positive neighbourhoods have no supported maximum and return None.
    This is an intensity maximum, not a centroid or chemical identification.
    """
    if columns < 0 or int(columns) != columns or not np.isfinite(t2_radius) or t2_radius <= 0:
        raise ValueError("Search limits must be non-negative whole columns and positive seconds")
    if image.shape != fold.shape or not np.isfinite([t1, t2]).all():
        raise ValueError("An anchor needs finite coordinates and an unextended folded image")
    centre = int(np.floor((t1-fold.t1[0])/fold.spec.period))
    cols = np.arange(max(0, centre-columns), min(fold.n_cols, centre+columns+1), dtype=int)
    rows = np.arange(int(np.ceil((t2-t2_radius)/fold.row_pitch-1e-9)),
                     int(np.floor((t2+t2_radius)/fold.row_pitch+1e-9))+1)
    rr, cc = np.meshgrid(rows, cols, indexing="ij")
    actual_col = cc + rr // fold.n_rows
    valid = (actual_col >= 0) & (actual_col < fold.n_cols)
    rr, cc, actual_col = rr[valid], cc[valid], actual_col[valid]
    values = image[rr % fold.n_rows, actual_col]
    good = np.isfinite(values)
    if not good.any() or values[good].max() <= 0 or np.ptp(values[good]) == 0:
        return None
    best = np.flatnonzero(good & (values == values[good].max()))
    distance = ((fold.t1[cc[best]]+fold.spec.period/2-t1)/fold.spec.period)**2
    distance += ((rr[best]*fold.row_pitch-t2)/t2_radius)**2
    chosen = best[np.argmin(distance)]
    return float(fold.t1[cc[chosen]]+fold.spec.period/2), float(rr[chosen]*fold.row_pitch)


@dataclass
class AlignmentTransform:
    """A fitted map between a run's true coordinates and a reference
    run's. `forward(run_t1, run_t2)` gives the reference coordinates;
    `inverse(ref_t1, ref_t2)` gives the run's. Both take and return
    scalars or numpy arrays."""
    kind: str
    anchors: list[AnchorPoint]
    forward: Callable = field(repr=False, compare=False)
    inverse: Callable = field(repr=False, compare=False)
    monotonic: bool = True
    warning: str = ""
    dt1: float | None = None
    dt2: float | None = None
    a1: float | None = None
    b1: float | None = None
    a2: float | None = None
    b2: float | None = None

    def residuals(self) -> list[dict]:
        """How well the fitted transform predicts each anchor it was
        built from (in-sample; see `cross_validate` for an honest,
        out-of-sample error)."""
        out = []
        for p in self.anchors:
            t1p, t2p = self.forward(p.run_t1, p.run_t2)
            out.append(dict(name=p.name, ref_t1=p.ref_t1, ref_t2=p.ref_t2,
                            run_t1=p.run_t1, run_t2=p.run_t2,
                            resid_t1=float(t1p) - p.ref_t1, resid_t2=float(t2p) - p.ref_t2))
        return out

    def summary(self) -> str:
        n = len(self.anchors)
        s = f"{self.kind} from {n} anchor{'s' if n != 1 else ''}"
        if not self.monotonic:
            s += f"; {self.warning}"
        r1, r2 = cross_validate(self.anchors, self.kind)
        if np.isfinite(r1):
            s += f"; leave-one-out error {r1 * 1e3:.0f} ms in t1, {r2 * 1e3:.0f} ms in t2"
        else:
            s += "; too few anchors to cross-validate"
        return s


# ------------------------------------------------------------------ fitting
def _huber_wls(X: np.ndarray, y: np.ndarray, iters: int = 8) -> np.ndarray:
    """Iteratively reweighted least squares with a Huber weight, so one
    anchor far from the rest doesn't dominate the fit."""
    w = np.ones(y.size)
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    for _ in range(iters):
        sw = np.sqrt(w)
        beta = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        resid = y - X @ beta
        mad = 1.4826 * np.median(np.abs(resid - np.median(resid)))
        if mad <= 1e-12:
            break
        c = 1.345 * mad
        w = np.minimum(1.0, c / np.maximum(np.abs(resid), 1e-12))
    return beta


def _fit_shift(points: list[AnchorPoint]):
    t1s = np.array([p.run_t1 for p in points])
    t2s = np.array([p.run_t2 for p in points])
    t1r = np.array([p.ref_t1 for p in points])
    t2r = np.array([p.ref_t2 for p in points])
    dt1 = float(np.median(t1r - t1s))
    dt2 = float(np.median(t2r - t2s))

    def fwd(t1, t2):
        return np.asarray(t1, dtype=float) + dt1, np.asarray(t2, dtype=float) + dt2

    def inv(t1, t2):
        return np.asarray(t1, dtype=float) - dt1, np.asarray(t2, dtype=float) - dt2

    return fwd, inv, {"dt1": dt1, "dt2": dt2, "monotonic": True, "warning": ""}


def _fit_affine(points: list[AnchorPoint]):
    t1s = np.array([p.run_t1 for p in points])
    t2s = np.array([p.run_t2 for p in points])
    t1r = np.array([p.ref_t1 for p in points])
    t2r = np.array([p.ref_t2 for p in points])
    X1 = np.column_stack([np.ones_like(t1s), t1s])
    a1, b1 = (float(v) for v in _huber_wls(X1, t1r))
    c_d = _huber_wls(X1, t2r - t2s)
    a2, b2 = float(c_d[0]), float(c_d[1])
    monotonic = b1 > 0
    warning = "" if monotonic else ("this run's t1 anchors are not increasing in the same order as "
                                    "the reference's (the fitted stretch is negative)")
    b1_safe = b1 if abs(b1) > 1e-12 else (1e-12 if b1 >= 0 else -1e-12)

    def fwd(t1, t2):
        t1 = np.asarray(t1, dtype=float)
        t2 = np.asarray(t2, dtype=float)
        return a1 + b1 * t1, t2 + a2 + b2 * t1

    def inv(t1, t2):
        t1 = np.asarray(t1, dtype=float)
        t2 = np.asarray(t2, dtype=float)
        t1s_ = (t1 - a1) / b1_safe
        return t1s_, t2 - (a2 + b2 * t1s_)

    return fwd, inv, {"a1": a1, "b1": b1, "a2": a2, "b2": b2, "monotonic": monotonic, "warning": warning}


def _fit_smooth(points: list[AnchorPoint]):
    order = np.argsort([p.run_t1 for p in points])
    t1s = np.array([points[i].run_t1 for i in order])
    t2s = np.array([points[i].run_t2 for i in order])
    t1r = np.array([points[i].ref_t1 for i in order])
    t2r = np.array([points[i].ref_t2 for i in order])
    monotonic = bool(np.all(np.diff(t1s) > 0) and np.all(np.diff(t1r) > 0))
    warning = "" if monotonic else ("the t1 anchors are not increasing together in both runs; "
                                    "the smooth curve is not increasing and cannot be trusted")
    if not monotonic:
        # A safe (if approximate) fallback: sort and drop ties so something can still be built.
        keep = np.concatenate([[True], np.diff(t1s) > 0])
        t1s, t2s, t1r, t2r = t1s[keep], t2s[keep], t1r[keep], t2r[keep]
        keep2 = np.concatenate([[True], np.diff(t1r) > 0])
        t1s, t2s, t1r, t2r = t1s[keep2], t2s[keep2], t1r[keep2], t2r[keep2]
    if len(t1s) < 2:
        fwd, inv, extra = _fit_shift(points)
        extra = {**extra, "monotonic": False, "warning": warning or extra["warning"]}
        return fwd, inv, extra
    t1_curve = ip.PchipInterpolator(t1s, t1r, extrapolate=True)
    try:
        off_curve = ip.make_smoothing_spline(t1s, t2r - t2s)
    except Exception:  # noqa: BLE001  (older scipy without make_smoothing_spline, or too few points)
        k = min(3, len(t1s) - 1)
        off_curve = ip.UnivariateSpline(t1s, t2r - t2s, k=k, s=len(t1s))
    lo, hi = t1s[0], t1s[-1]
    grid = np.linspace(lo, hi, max(400, 40 * len(t1s)))
    grid_ref = t1_curve(grid)
    if np.any(np.diff(grid_ref) <= 0):
        monotonic = False
        warning = warning or "the smooth t1 curve is not increasing across its range"
        order2 = np.argsort(grid_ref, kind="stable")
        grid, grid_ref = grid[order2], grid_ref[order2]
        keep3 = np.concatenate([[True], np.diff(grid_ref) > 0])
        grid, grid_ref = grid[keep3], grid_ref[keep3]
    inv_curve = ip.PchipInterpolator(grid_ref, grid, extrapolate=True)

    def fwd(t1, t2):
        t1 = np.asarray(t1, dtype=float)
        t2 = np.asarray(t2, dtype=float)
        return t1_curve(t1), t2 + off_curve(t1)

    def inv(t1, t2):
        t1 = np.asarray(t1, dtype=float)
        t2 = np.asarray(t2, dtype=float)
        t1s_ = inv_curve(t1)
        return t1s_, t2 - off_curve(t1s_)

    return fwd, inv, {"monotonic": monotonic, "warning": warning}


_FITTERS: dict[str, Callable] = {"shift": _fit_shift, "affine": _fit_affine, "smooth": _fit_smooth}


def _build(kind: str, points: list[AnchorPoint]) -> AlignmentTransform:
    fwd, inv, extra = _FITTERS[kind](points)
    return AlignmentTransform(kind=kind, anchors=list(points), forward=fwd, inverse=inv, **extra)


def fit_alignment(points: list[AnchorPoint], kind: str = "auto") -> AlignmentTransform:
    if kind not in ("auto", *FIT_MIN):
        raise AlignmentError(f"unknown alignment kind '{kind}'")
    n = len(points)
    if kind != "auto":
        if n < FIT_MIN[kind]:
            raise AlignmentError(f"{kind} needs at least {FIT_MIN[kind]} anchors, got {n}")
        return _build(kind, points)
    if n < FIT_MIN["shift"]:
        raise AlignmentError("need at least one anchor")
    # Judge each eligible tier by its LEAVE-ONE-OUT error, not how well it fits
    # the anchors it was built from: a smooth (PCHIP) curve passes exactly
    # through every anchor by construction, so an in-sample check would
    # always call it "good enough" even when it wildly overshoots between
    # sparse, noisy anchors. Cross-validation catches that overfitting;
    # in-sample fit quality would not.
    candidates = []
    for tier in TIER_ORDER:
        if n < AUTO_MIN[tier]:
            continue
        tf = _build(tier, points)
        r1, r2 = cross_validate(points, tier)
        candidates.append((tier, tf, r1, r2))
        if r1 <= AUTO_T1_TOL and r2 <= AUTO_T2_TOL:
            return tf

    def relative_error(c):
        _tier, _tf, r1, r2 = c
        r1 = r1 / AUTO_T1_TOL if np.isfinite(r1) else np.inf
        r2 = r2 / AUTO_T2_TOL if np.isfinite(r2) else np.inf
        return r1 + r2

    return min(candidates, key=relative_error)[1]


def cross_validate(points: list[AnchorPoint], kind: str) -> tuple[float, float]:
    """Leave-one-out RMS error (seconds) of `kind` at predicting an
    anchor it wasn't fitted on. NaN if fewer than `kind`'s own minimum
    anchors are left over in any fold (so no fold could be computed)."""
    if kind not in FIT_MIN:
        raise AlignmentError(f"unknown alignment kind '{kind}'")
    need = FIT_MIN[kind]
    e1, e2 = [], []
    for i in range(len(points)):
        rest = points[:i] + points[i + 1:]
        if len(rest) < need:
            continue
        try:
            tf = _build(kind, rest)
        except AlignmentError:
            continue
        t1p, t2p = tf.forward(points[i].run_t1, points[i].run_t2)
        e1.append(float(t1p) - points[i].ref_t1)
        e2.append(float(t2p) - points[i].ref_t2)
    if not e1:
        return float("nan"), float("nan")
    e1, e2 = np.array(e1), np.array(e2)
    return float(np.sqrt(np.mean(e1 ** 2))), float(np.sqrt(np.mean(e2 ** 2)))


# --------------------------------------------------------------- propagation
def propagate_anchors(anchors: list[Anchor], run: Run, fold: Fold,
                      min_height: float = 50.0) -> dict[str, AnchorPoint | None]:
    """Look for each reference anchor's compound in `run`, on its own ion
    within its own window (the same search `locate_anchor` uses for
    phase), and report where it actually is here. A name maps to `None`
    if it wasn't found above `min_height`. The reported `run_t1` is the
    midpoint of whichever modulation the peak's apex falls in — the
    natural resolution for "where in the run this is", one modulation
    wide either way."""
    out: dict[str, AnchorPoint | None] = {}
    for a in anchors:
        t_obs, _height = locate_anchor(run, a, min_height=min_height)
        if t_obs is None:
            out[a.name] = None
            continue
        col, t2 = fold.to_fold(t_obs)
        run_t1 = float(fold.t1[int(col)] + 0.5 * fold.spec.period)
        out[a.name] = AnchorPoint(a.name, a.t1, a.t2, run_t1, float(t2), int(a.mz))
    return out


# ------------------------------------------------------------------ resampling
def blank_mapping(sample_transform: AlignmentTransform | None, blank_transform: AlignmentTransform | None):
    """A function taking a sample-run coordinate (t1, t2) to the blank
    run's, by way of the shared reference: sample -> reference with the
    sample's own transform, reference -> blank with the blank's inverse. A
    `None` transform means that run is the reference (or has none), so its
    coordinates are used as they are. `None` for both is the identity."""
    def to_blank(t1, t2):
        if sample_transform is not None:
            t1, t2 = sample_transform.forward(t1, t2)
        if blank_transform is not None:
            t1, t2 = blank_transform.inverse(t1, t2)
        return t1, t2
    return to_blank


def resample_between(image: np.ndarray, src_fold: Fold, dst_fold: Fold, to_src: Callable) -> np.ndarray:
    """`image` (folded by `src_fold`) read onto `dst_fold`'s pixel grid:
    every destination pixel is looked up at `to_src(t1, t2)` (its coordinates
    in the source run) by bilinear interpolation. t2 wraps round the source's
    period; pixels outside the source's t1 range come back as 0. The
    destination and source may differ in phase and start time."""
    t1_grid, t2_grid = np.meshgrid(dst_fold.t1, np.arange(dst_fold.n_rows) * dst_fold.row_pitch)
    t1_src, t2_src = to_src(t1_grid.ravel(), t2_grid.ravel())
    col_src = (np.asarray(t1_src) - src_fold.t1[0]) / src_fold.spec.period
    row_src = np.asarray(t2_src) / src_fold.row_pitch
    k = np.floor(row_src / src_fold.n_rows)
    row_src = row_src - k * src_fold.n_rows
    col_src = col_src + k
    out = ndi.map_coordinates(np.nan_to_num(image, nan=0.0), [row_src, col_src], order=1,
                              mode="constant", cval=0.0)
    bad = (col_src < 0) | (col_src > src_fold.n_cols - 1)
    return np.where(bad, 0.0, out).reshape(dst_fold.n_rows, dst_fold.n_cols)


def resample_to_reference(image: np.ndarray, run_fold: Fold, ref_fold: Fold,
                          transform: AlignmentTransform, extended: bool = False,
                          extend_factor: float = DEFAULT_EXTEND_FACTOR) -> np.ndarray:
    """`image` (as folded by `run_fold`) resampled onto `ref_fold`'s grid:
    every reference pixel is read from wherever `transform` says it came
    from in the run, by bilinear interpolation. With `extended`, the
    reference t2 axis reaches past one period (as `Fold.extend` does for
    display) so a wrapped peak shows up whole; pixels with no counterpart
    in the run (outside its own t1 range) come back as 0."""
    n_cols = ref_fold.n_cols
    n_rows = ref_fold.n_rows
    if extended:
        n_rows = n_rows + int(round((extend_factor - 1.0) * ref_fold.n_rows))
    t1_ref = ref_fold.t1
    t2_ref = np.arange(n_rows) * ref_fold.row_pitch
    t1_grid, t2_grid = np.meshgrid(t1_ref, t2_ref)
    t1_run, t2_run = transform.inverse(t1_grid.ravel(), t2_grid.ravel())
    col_src = (t1_run - run_fold.t1[0]) / run_fold.spec.period
    row_src = t2_run / run_fold.row_pitch
    k = np.floor(row_src / run_fold.n_rows)
    row_src = row_src - k * run_fold.n_rows
    col_src = col_src + k
    out = ndi.map_coordinates(np.nan_to_num(image, nan=0.0), [row_src, col_src], order=1,
                              mode="constant", cval=0.0)
    bad = (col_src < 0) | (col_src > run_fold.n_cols - 1)
    return np.where(bad, 0.0, out).reshape(n_rows, n_cols)
