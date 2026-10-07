"""The n-alkane ladder, homologous series, and group-type stencils.

The idea in plain terms
------------------------
A GC×GC hydrocarbon chromatogram is usually a big diagonal "smudge": every
homologous series (n-alkanes, branched alkanes, cycloalkanes, ...) climbs
across the plot together, one series per rough band of t2, because adding
a CH2 group moves a compound almost the same distance in t1 whichever
series it belongs to. The n-alkanes are the one series you can name with
confidence (they're what the standard is run for), so everything else in
this module is described *relative to them*: a carbon-number-equivalent
position along t1 (Van den Dool interpolation between known alkane
apexes), and an offset in t2 measured from the alkane line itself, not
from zero. That coordinate system is what makes a stencil drawn on one
run (typically the alkane standard, where the ladder is clean) transfer
onto another run's own smudge: the stencil is recorded in that relative
space, and "lifting" it onto a different run just re-reads its vertices
against *that* run's own alkane ladder.

Finding the ladder (`find_alkane_series`, `propagate_alkane_series`)
---------------------------------------------------------------------
Given two confirmed points (carbon number, t1 and t2 of where you clicked),
the ladder grows outward one carbon at a time. Each next carbon is predicted
from the straight line through its two nearest confirmed neighbours (in t1
and in t2; see `_predict` for why not a quadratic) and searched for on the alkane quantifier
ion (m/z 57 by default). The search is two-dimensional: it looks only at scans
whose t2 lies within `t2_window_s` of the predicted t2 and takes the tallest
pulse whose qualifier ratios (m/z 71/57 and 85/57 by default) match the
expected n-alkane pattern. This matters on real hydrocarbon samples, where a
brighter m/z 57 pulse from a branched or cyclic compound at another t2 is
almost always nearby: a time-only search (as phase estimation uses for
internal standards) locks onto that instead. The expected ratios can be edited,
or measured from a peak you click (see `pixelspec.spectrum_at`).
`propagate_alkane_series` does the same search on another run, centred on
the reference run's own found positions rather than a prediction, which
is the same simplification `align.propagate_anchors` makes for internal
standards: retention drifts a little run to run, not by a whole modulation.

Stencils (`Stencil`, `stencil_from_polygon`, `Stencil.to_absolute`)
--------------------------------------------------------------------
A stencil is a polygon in (carbon, t2 offset from the alkane line) space.
Draw one on a run with a known `AlkaneSeries` (`stencil_from_polygon`
records it in that space); `to_absolute` reads it back out against *any*
run's own `AlkaneSeries`, including one that's never seen the polygon
drawn on it. `stencil_totals` then sums the TIC volume of every blob
whose apex falls inside the resulting real-coordinate polygon.

Homologous series (`find_homologous_series`) and markers (`find_markers`)
---------------------------------------------------------------------------
Other series are found the cheap way: candidate blobs carrying a strong
share of a diagnostic ion pair are placed in carbon-equivalent space (the
alkane ladder's own t1-to-carbon map, so the shrinking spacing through a
temperature programme is handled) and must sit one carbon apart at a
consistent t2 offset from the alkane line. Every candidate is tried in turn
as the reference phase and the phase explaining the most carbons wins (a
small, exhaustive stand-in for RANSAC). Markers (pristane, phytane) are
looked up at a measured carbon-equivalent offset from their neighbouring
n-alkane (see `MARKER_OFFSET_CARBONS`).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.signal import find_peaks

from .fold import Fold
from .peaks import BlobSet
from .pixelspec import ratios_to, spectrum_at
from .run import Run

MARKER_NEAR_ALKANE: dict[str, int] = {"pristane": 17, "phytane": 18}
# Carbon-equivalent offset of each marker from its neighbouring n-alkane, measured on the
# supplied real acquisition from the t1 centroid of the m/z 183 excess (a branched-alkane
# ion) against the alkane's own m/z 57 centroid: pristane +0.006 C, phytane +0.08 C. Both
# are within one 10 s modulation (about 0.045 C) of the alkane, i.e. they co-elute with
# n-C17 and n-C18 in the first dimension; literature RIs on a 5% phenyl phase (1708 and
# 1812) give +0.08 and +0.12. The earlier guess of -0.3 C put the search a whole
# modulation column or more before the peak.
MARKER_OFFSET_CARBONS: dict[str, float] = {"pristane": 0.01, "phytane": 0.08}


@dataclass(frozen=True)
class HydrocarbonSettings:
    quant_ion: int = 57
    qual_ions: tuple[int, ...] = (71, 85)
    # Expected n-alkane ratios to the quantifier, in the order of `qual_ions`.
    # Measured on the earlier acquisition (docs/validation/hydrocarbons-2026-09-20.md). Re-measured on
    # the September 2026 runs the n-alkane medians are lower (0.60, 0.38), but a gate centred there
    # admits branched isomers with lower ratios still and gave worse ladders (14 eV molecular-ion
    # check: 25 carbons confirmed and 3 contradicted, against 36 and 1), so these values stay
    # (benchmarks/real_data_audit/alkane_ratios.py).
    qual_ratios: tuple[float, ...] = (0.66, 0.44)
    qual_tol: float = 0.3             # allowed relative deviation of each ratio
    enforce_qualifiers: bool = True   # reject a candidate whose ratios are outside the tolerance
    carbon_range: tuple[int, int] = (8, 40)
    t1_window_s: float = 45.0         # per-carbon search window in t1 (widens with the gap to the nearest hit)
    t2_window_s: float = 0.25         # ...and in t2, around the predicted alkane line
    seed_window_s: float = 6.0        # search window when refining a clicked seed
    seed_t2_window_s: float = 0.5
    propagate_window_s: float = 15.0  # search window when propagating to another run
    max_gap: int = 2                  # carbons that may be missed in a row before that direction stops
    min_height: float = 50.0
    refine_passes: int = 2
    diagnostic_frac_min: float = 0.1  # a homologous-series candidate needs this share on its major ion
    pair_min_ratio: float = 0.02      # a candidate needs its minor diagnostic ion at least this share of the major
    series_tol_carbon: float = 0.15   # how far from a whole-carbon slot (in carbon-equivalents, 0.15 = 15 RI) still counts
    series_t2_tol_s: float = 0.5      # allowed spread of a series' t2 offset from the alkane line
    min_series_members: int = 4
    ladder_shape: str = "piecewise"   # path drawn through the alkanes: see `LADDER_SHAPES`
    ladder_degree: int = 3            # polynomial degree when `ladder_shape` is "polynomial"
    assign_carbon_tol: float = 0.15   # a blob this close (carbon-equivalents) to an alkane can be that n-alkane
    assign_t2_tol_s: float = 0.3      # ...and this close in t2 to the alkane line

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "HydrocarbonSettings":
        d = {k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__}
        for key in ("qual_ions", "qual_ratios", "carbon_range"):
            if key in d:
                d[key] = tuple(d[key])
        return cls(**d)


# ------------------------------------------------------------------ ladder
@dataclass
class AlkaneHit:
    carbon: int
    t1: float                       # s, mid-modulation (fold convention)
    t2: float                       # s
    height: float
    qual_ratios: tuple[float, ...]
    seeded: bool = False            # a user-given point, not found by search
    qual_ok: bool = True            # ratios within tolerance of the expected n-alkane pattern


def _interp_extrap(x: float, xp: np.ndarray, fp: np.ndarray) -> float:
    """Linear interpolation between known points, linear extrapolation
    (from the nearest segment's slope) beyond them."""
    if x < xp[0]:
        slope = (fp[1] - fp[0]) / (xp[1] - xp[0]) if xp.size > 1 else 0.0
        return float(fp[0] + slope * (x - xp[0]))
    if x > xp[-1]:
        slope = (fp[-1] - fp[-2]) / (xp[-1] - xp[-2]) if xp.size > 1 else 0.0
        return float(fp[-1] + slope * (x - xp[-1]))
    return float(np.interp(x, xp, fp))


LADDER_SHAPES = {"piecewise": "Straight between points", "smooth": "Smooth curve through points",
                 "polynomial": "Polynomial fit"}


def _curve_value(x: float, xp: np.ndarray, fp: np.ndarray, shape: str, degree: int) -> float:
    """Value at `x` of the curve through (`xp`, `fp`), `xp` increasing.

    piecewise: straight lines between the points (the original behaviour).
    smooth: a monotone cubic (PCHIP) through every point, so a ladder that bends or
        wobbles is followed without overshooting between points.
    polynomial: a least-squares polynomial of `degree`, which averages over a wonky
        point instead of following it.
    Beyond the first or last point every shape continues along the end slope, because
    a polynomial or spline extrapolated freely runs away within a few carbons."""
    if xp.size < 3 or shape == "piecewise" or np.any(np.diff(xp) <= 0):    # repeated positions: no curve to fit
        return _interp_extrap(x, xp, fp)
    if shape == "polynomial":
        centre, scale = float(xp.mean()), float(xp.std()) or 1.0      # scaled so a high degree stays well conditioned
        poly = np.polynomial.Polynomial.fit((xp - centre) / scale, fp, min(degree, xp.size - 1))
        inside = lambda v: float(poly((v - centre) / scale))          # noqa: E731
        slope = lambda v: float(poly.deriv()((v - centre) / scale) / scale)    # noqa: E731
    else:
        spline = PchipInterpolator(xp, fp)
        inside = lambda v: float(spline(v))                           # noqa: E731
        slope = lambda v: float(spline.derivative()(v))               # noqa: E731
    if x < xp[0]:
        return inside(xp[0]) + slope(xp[0]) * (x - xp[0])
    if x > xp[-1]:
        return inside(xp[-1]) + slope(xp[-1]) * (x - xp[-1])
    return inside(x)


@dataclass
class AlkaneSeries:
    """The n-alkane ladder found in one run: confirmed carbon numbers and
    where each one actually is. `shape` and `degree` say what path is drawn
    through them (see `_curve_value`); they change how positions between and
    beyond the hits are read, never the hits themselves."""
    hits: dict[int, AlkaneHit] = field(default_factory=dict)
    shape: str = "piecewise"
    degree: int = 3
    _inverse: tuple | None = field(default=None, init=False, repr=False, compare=False)

    def __len__(self) -> int:
        return len(self.hits)

    def carbons(self) -> np.ndarray:
        return np.array(sorted(self.hits), dtype=np.int64)

    def _arrays(self):
        cs = self.carbons()
        if cs.size < 2:
            raise ValueError("need at least two confirmed alkanes to interpolate")
        t1s = np.array([self.hits[int(c)].t1 for c in cs])
        t2s = np.array([self.hits[int(c)].t2 for c in cs])
        return cs.astype(float), t1s, t2s

    def t1_of(self, carbon: float) -> float:
        cs, t1s, _ = self._arrays()
        return _curve_value(float(carbon), cs, t1s, self.shape, self.degree)

    def t2_of(self, carbon: float) -> float:
        """The alkane line's own t2 at this carbon number (interpolated)."""
        cs, _, t2s = self._arrays()
        return _curve_value(float(carbon), cs, t2s, self.shape, self.degree)

    def carbon_of(self, t1: float) -> float:
        cs, t1s, _ = self._arrays()
        if self.shape == "piecewise":
            order = np.argsort(t1s)
            return _interp_extrap(float(t1), t1s[order], cs[order])
        # A curve fitted to t1 against carbon is inverted numerically: fitting carbon against t1
        # separately would not be the same curve and drifts away from the hits.
        key = (self.shape, self.degree, tuple(cs), tuple(t1s))
        if self._inverse is None or self._inverse[0] != key:      # built once per set of hits, not per blob
            grid = np.linspace(cs[0] - 3.0, cs[-1] + 3.0, 600)
            times = np.maximum.accumulate([_curve_value(c, cs, t1s, self.shape, self.degree) for c in grid])
            keep = np.concatenate(([True], np.diff(times) > 0))   # a wobbling curve must still be invertible
            self._inverse = (key, times[keep], grid[keep])
        return _interp_extrap(float(t1), self._inverse[1], self._inverse[2])

    def to_relative(self, t1: float, t2: float) -> tuple[float, float]:
        c = self.carbon_of(t1)
        return c, float(t2) - self.t2_of(c)

    def to_absolute(self, carbon: float, dt2: float) -> tuple[float, float]:
        return self.t1_of(carbon), self.t2_of(carbon) + float(dt2)

    def step_s(self) -> float:
        """Average per-carbon t1 spacing across the confirmed range."""
        cs, t1s, _ = self._arrays()
        return float((t1s[-1] - t1s[0]) / (cs[-1] - cs[0]))


def qualifier_ratios_at(run: Run, fold: Fold, t1: float, t2: float,
                        s: HydrocarbonSettings = HydrocarbonSettings()) -> tuple[float, ...]:
    """The alkane qualifier ions' intensities relative to the quantifier at a
    point (baseline-subtracted, summed over a few rows). This is what
    'use the peak I clicked as the reference' measures."""
    spec = spectrum_at(run, fold, t1, t2)
    return ratios_to(spec, run.mz, s.quant_ion, s.qual_ions)


def ratios_match(ratios, s: HydrocarbonSettings) -> bool:
    """True if every ratio is within `qual_tol` (relative) of its expected value.
    No expected ratios means nothing to check."""
    if not s.qual_ratios:
        return True
    exp = np.asarray(s.qual_ratios, dtype=float)
    got = np.asarray(ratios, dtype=float)
    if got.size != exp.size:
        return True
    return bool(np.all(np.abs(got - exp) <= s.qual_tol * exp))


def _apex_time(x: np.ndarray, times: np.ndarray, i: int) -> float:
    """Sub-scan apex of the pulse whose top scan is `i`: a parabola through
    the log of the three highest points (a Gaussian is a parabola in log space)."""
    if i <= 0 or i >= x.size - 1:
        return float(times[i])
    a, b, c = x[i - 1], x[i], x[i + 1]
    if a > 0 and c > 0:
        a, b, c = np.log(a), np.log(b), np.log(c)
    den = a - 2 * b + c
    delta = float(np.clip(0.5 * (a - c) / den, -0.5, 0.5)) if den != 0 else 0.0
    step = times[i + 1] - times[i] if delta > 0 else times[i] - times[i - 1]
    return float(times[i] + delta * step)


def _search_carbon(run: Run, fold: Fold, carbon: int, t1_guess: float, t2_guess: float,
                   t1_window: float, t2_window: float, s: HydrocarbonSettings,
                   enforce: bool | None = None, avoid=(), avoid_tol: float = 0.0) -> AlkaneHit | None:
    """Tallest quantifier-ion pulse within `t1_window` seconds and `t2_window`
    of (t1_guess, t2_guess) (t1 is the middle of the modulation, as everywhere
    else) whose qualifier ratios match. With `enforce` False the tallest pulse
    is returned even if they don't (flagged `qual_ok=False`). A pulse within
    `avoid_tol` seconds of any t1 in `avoid` (an alkane already assigned to
    another carbon) is skipped."""
    enforce = s.enforce_qualifiers if enforce is None else enforce
    period = fold.spec.period
    col = int(np.floor((t1_guess - fold.spec.t0 - fold.spec.phase) / period))
    t_expect = float(fold.to_time(col, t2_guess))
    lo, hi = np.searchsorted(run.times, [t_expect - t1_window, t_expect + t1_window])
    if hi - lo < 5:
        return None
    x = run.xic(s.quant_ion)[lo:hi]
    _, t2s = fold.to_fold(run.times[lo:hi])
    y = np.where(np.abs(t2s - t2_guess) <= t2_window, x, 0.0)
    peaks, _ = find_peaks(y, height=s.min_height)
    if peaks.size == 0:
        return None
    first = None
    for i in peaks[np.argsort(y[peaks])[::-1]][:5]:
        t_obs = _apex_time(x, run.times[lo:hi], int(i))
        c, t2 = fold.to_fold(t_obs)
        t1_mid = float(fold.t1[int(c)] + 0.5 * period)
        if any(abs(t1_mid - a) < avoid_tol for a in avoid):
            continue
        ratios = qualifier_ratios_at(run, fold, t1_mid, float(t2), s)
        ok = ratios_match(ratios, s)
        hit = AlkaneHit(carbon=carbon, t1=t1_mid, t2=float(t2), height=float(x[i]),
                        qual_ratios=ratios, qual_ok=ok)
        if ok:
            return hit
        first = first or hit
    return None if enforce else first


def _predict(series: AlkaneSeries, carbon: int) -> tuple[float, float, float]:
    """Expected (t1, t2) of `carbon` and the local per-carbon t1 spacing, from
    the two nearest confirmed alkanes: the local straight line (secant)
    through them. Alkane spacing shrinks steadily through a programme, so a
    line through more distant points is biased; a quadratic would follow the
    curvature but t1 is only known to a modulation (10 s here) and a quadratic
    amplifies that noise when extrapolating. The search window is wide enough
    to absorb both."""
    cs, t1s, t2s = series._arrays()
    order = np.argsort(np.abs(cs - carbon))[:2]
    x = cs[order]
    t1_fit, t2_fit = np.polyfit(x, t1s[order], 1), np.polyfit(x, t2s[order], 1)
    return float(np.polyval(t1_fit, carbon)), float(np.polyval(t2_fit, carbon)), abs(float(t1_fit[0]))


def find_alkane_series(run: Run, fold: Fold, seeds: list[tuple[int, float, float]],
                       settings: HydrocarbonSettings = HydrocarbonSettings()) -> AlkaneSeries:
    """Find the n-alkane ladder from two (or more) confirmed points.

    `seeds`: (carbon number, t1 (mid-modulation, s), t2 (s)) for each
    compound you've already identified (typically two clicks). Each seed
    is snapped to the true local apex (its qualifier ratios are reported
    but not enforced: you said it is an alkane), then the ladder grows
    outward from the confirmed carbons one step at a time, each prediction
    coming from the nearest confirmed neighbours. A carbon that can't be found
    widens the next prediction's window; after `max_gap` misses in a row
    that direction stops. `refine_passes` rounds retry any carbon still
    missing once its neighbours are known.
    """
    s = settings
    if len(seeds) < 2:
        raise ValueError("need at least two seed points to start a series")
    series = AlkaneSeries()
    for carbon, t1_guess, t2_guess in seeds:
        hit = _search_carbon(run, fold, int(carbon), t1_guess, t2_guess, s.seed_window_s,
                             s.seed_t2_window_s, s, enforce=False)
        if hit is None:
            raise ValueError(f"couldn't confirm the seed at C{carbon} near t1={t1_guess:.1f} s")
        hit.seeded = True
        series.hits[int(carbon)] = hit
    # A ladder grown from seeds in the wrong order runs backwards through the run.
    cs = series.carbons()
    t1s = [series.hits[int(c)].t1 for c in cs]
    if len(cs) < len(seeds) or np.any(np.diff(t1s) <= 0):
        pairs = ", ".join(f"C{c} at {t / 60:.2f} min" for c, t in zip(cs, t1s))
        raise ValueError(f"the seeds must be different peaks that elute in carbon order ({pairs})")

    lo, hi = s.carbon_range
    for _ in range(max(s.refine_passes, 1)):
        tried: set[int] = set()
        found_any = False
        while True:
            todo = [(min(abs(c - k) for k in series.hits), c) for c in range(lo, hi + 1)
                    if c not in series.hits and c not in tried]
            todo = [(d, c) for d, c in todo if d <= s.max_gap + 1]
            if not todo:
                break
            d, c = min(todo)
            tried.add(c)
            t1_pred, t2_pred, step = _predict(series, c)
            taken = [h.t1 for h in series.hits.values()]
            # Never reach as far as the neighbouring carbon, however tight the spacing.
            window = min(s.t1_window_s, 0.45 * step) * d
            hit = _search_carbon(run, fold, c, t1_pred, t2_pred, window, s.t2_window_s * d, s,
                                 avoid=taken, avoid_tol=0.4 * step)
            if hit is not None:
                series.hits[c] = hit
                found_any = True
        if not found_any:
            break
    return series


def propagate_alkane_series(reference: AlkaneSeries, run: Run, fold: Fold,
                            settings: HydrocarbonSettings = HydrocarbonSettings()) -> AlkaneSeries:
    """The same carbons `reference` confirmed, searched for in `run` near
    the reference's own positions (retention drifts a little between
    runs, not by a whole modulation -- the same assumption
    `align.propagate_anchors` makes for internal standards)."""
    series = AlkaneSeries()
    for carbon, ref_hit in reference.hits.items():
        hit = _search_carbon(run, fold, carbon, ref_hit.t1, ref_hit.t2, settings.propagate_window_s,
                             settings.t2_window_s, settings)
        if hit is not None:
            series.hits[carbon] = hit
    return series


# ------------------------------------------------------------- homologues
@dataclass
class HomologSeries:
    name: str
    members: list[dict] = field(default_factory=list)   # blob, carbon, t1, t2, height


def find_homologous_series(blobs: BlobSet, ion_pair: tuple[int, int], series: AlkaneSeries,
                           name: str = "series",
                           settings: HydrocarbonSettings = HydrocarbonSettings()) -> HomologSeries | None:
    """A homologous series among blobs carrying both ions of `ion_pair` (the
    major ion at a share of the blob's volume of at least
    `diagnostic_frac_min`, the minor at least `pair_min_ratio` of the major).

    Members are spaced one carbon apart in *carbon-equivalent* space (the
    alkane ladder's own t1-to-carbon map, so the shrinking spacing through a
    temperature programme is handled), share a t2 offset from the alkane line
    within `series_t2_tol_s`, and there is at most one member per carbon (the
    tallest): a vertical streak of blobs at one t1, or several isomers at one
    carbon, cannot pass for a series. Every candidate is tried in turn as the
    reference phase and the phase explaining the most carbons wins. `None` if
    fewer than `min_series_members` line up."""
    s = settings
    major, minor = ion_pair
    mz0 = int(blobs.mz[0])
    n = len(blobs)
    if n == 0:
        return None
    ia, ib = major - mz0, minor - mz0
    if not (0 <= ia < blobs.spectra.shape[1] and 0 <= ib < blobs.spectra.shape[1]):
        return None
    majors = np.clip(blobs.spectra[:, ia], 0, None)
    minors = np.clip(blobs.spectra[:, ib], 0, None)
    share = majors / np.clip(blobs.volume, 1e-12, None)
    candidates = np.flatnonzero((majors > 0) & (share >= s.diagnostic_frac_min)
                                & (minors >= s.pair_min_ratio * majors))
    if candidates.size < s.min_series_members:
        return None
    carbon = np.array([series.carbon_of(float(t)) for t in blobs.t1[candidates]])
    dt2 = blobs.t2[candidates] - np.array([series.t2_of(c) for c in carbon])
    height = blobs.height[candidates]
    best_score, best_members = (0, 0.0), None
    for seed in range(candidates.size):
        phase = carbon[seed] % 1.0
        offset = (carbon - phase + 0.5) % 1.0 - 0.5          # signed distance to the nearest whole-carbon slot
        near = np.abs(offset) <= s.series_tol_carbon
        if near.sum() < s.min_series_members:
            continue
        near &= np.abs(dt2 - np.median(dt2[near])) <= s.series_t2_tol_s
        tallest: dict[int, int] = {}
        for i in np.flatnonzero(near):
            slot = int(np.round(carbon[i] - phase))
            if slot not in tallest or height[i] > height[tallest[slot]]:
                tallest[slot] = int(i)
        members = list(tallest.values())
        score = (len(members), float(height[members].sum()))
        if score > best_score:
            best_score, best_members = score, members
    if best_members is None or len(best_members) < s.min_series_members:
        return None
    members = [dict(blob=int(candidates[i]), carbon=float(carbon[i]), t1=float(blobs.t1[candidates[i]]),
                    t2=float(blobs.t2[candidates[i]]), height=float(height[i])) for i in best_members]
    members.sort(key=lambda m: m["t1"])
    return HomologSeries(name=name, members=members)


def find_markers(blobs: BlobSet, series: AlkaneSeries) -> dict[str, int | None]:
    """The blob nearest to where pristane and phytane are expected, a
    measured offset (`MARKER_OFFSET_CARBONS`) from n-C17 and n-C18. At a 10 s
    modulation they share a column with their alkane, so the nearest blob is
    often the alkane's own: confirm with the m/z 183 and 113 excess."""
    out: dict[str, int | None] = {}
    for name, near in MARKER_NEAR_ALKANE.items():
        try:
            t1, t2 = series.to_absolute(near + MARKER_OFFSET_CARBONS[name], 0.0)
        except ValueError:
            out[name] = None
            continue
        out[name] = blobs.nearest(t1, t2)
    return out


# --------------------------------------------------------------- stencils
@dataclass
class Stencil:
    name: str
    vertices: list[tuple[float, float]]   # (carbon, t2 offset from the alkane line)
    response_factor: float = 1.0

    def to_absolute(self, series: AlkaneSeries) -> list[tuple[float, float]]:
        return [series.to_absolute(c, d) for c, d in self.vertices]

    def to_dict(self) -> dict:
        return dict(name=self.name, vertices=[list(v) for v in self.vertices],
                    response_factor=self.response_factor)

    @classmethod
    def from_dict(cls, d: dict) -> "Stencil":
        return cls(name=d["name"], vertices=[tuple(v) for v in d["vertices"]],
                   response_factor=float(d.get("response_factor", 1.0)))


def stencil_from_polygon(vertices_t1t2: list[tuple[float, float]], series: AlkaneSeries,
                         name: str, response_factor: float = 1.0) -> Stencil:
    """Record a polygon drawn in real (t1, t2) on one run as a `Stencil`
    in that run's own alkane-relative coordinates -- the space any other
    run's own `AlkaneSeries` can translate it back out of."""
    return Stencil(name=name, vertices=[series.to_relative(t1, t2) for t1, t2 in vertices_t1t2],
                   response_factor=response_factor)


def _point_in_polygon(px: np.ndarray, py: np.ndarray, poly_x: np.ndarray, poly_y: np.ndarray) -> np.ndarray:
    """Ray-casting point-in-polygon test, vectorised over (px, py)."""
    px, py = np.atleast_1d(px), np.atleast_1d(py)
    inside = np.zeros(px.shape, dtype=bool)
    n = len(poly_x)
    j = n - 1
    for i in range(n):
        xi, yi, xj, yj = poly_x[i], poly_y[i], poly_x[j], poly_y[j]
        straddles = (yi > py) != (yj > py)
        with np.errstate(divide="ignore", invalid="ignore"):
            x_cross = xj + (py - yj) * (xi - xj) / (yi - yj)
        inside ^= straddles & (px < x_cross)
        j = i
    return inside


def stencil_totals(polygon_t1t2: list[tuple[float, float]], blobs: BlobSet) -> tuple[float, int, list[int]]:
    """Summed TIC volume, blob count, and blob indices for every blob
    whose apex falls inside a (real-coordinate) stencil polygon."""
    if len(blobs) == 0 or len(polygon_t1t2) < 3:
        return 0.0, 0, []
    xs = np.array([p[0] for p in polygon_t1t2])
    ys = np.array([p[1] for p in polygon_t1t2])
    inside = _point_in_polygon(blobs.t1, blobs.t2, xs, ys)
    idx = np.flatnonzero(inside)
    return float(blobs.volume[idx].sum()), int(idx.size), idx.tolist()


@dataclass
class StencilSet:
    """A named frame of stencils, built against one run's `AlkaneSeries`
    and lifted onto any other run's own."""
    name: str
    stencils: list[Stencil] = field(default_factory=list)

    def lift(self, series: AlkaneSeries) -> dict[str, list[tuple[float, float]]]:
        return {s.name: s.to_absolute(series) for s in self.stencils}

    def totals(self, series: AlkaneSeries, blobs: BlobSet) -> dict[str, dict]:
        lifted = self.lift(series)
        stencil_by_name = {s.name: s for s in self.stencils}
        out = {}
        for name, poly in lifted.items():
            volume, n, idx = stencil_totals(poly, blobs)
            rf = stencil_by_name[name].response_factor
            out[name] = dict(volume=volume, n_blobs=n, blob_indices=idx,
                             amount=volume * rf, polygon=poly)
        return out

    def to_dict(self) -> dict:
        return dict(name=self.name, stencils=[s.to_dict() for s in self.stencils])

    @classmethod
    def from_dict(cls, d: dict) -> "StencilSet":
        return cls(name=d.get("name", ""), stencils=[Stencil.from_dict(s) for s in d.get("stencils", [])])
