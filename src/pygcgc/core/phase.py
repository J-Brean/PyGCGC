"""Finding the modulation phase.

Two methods, which answer slightly different questions.

The void method finds a good place to cut. It folds the run with no phase
correction, adds up all the columns to get the typical second-dimension
profile, finds the widest stretch of that profile with almost no signal,
and puts the fold boundary just inside its far edge, so signal starts near
the bottom of the image. Nothing gets cut in half, which is what you want
for display and peak finding. It cannot say where the modulator's cycle
truly starts, because any cut inside the quiet stretch looks equally good.
A solvent front is left out of the profile, and columns are weighted so
that a few very bright ones cannot decide the answer. A `contrast` value
says how much quieter the cut is than the typical signal; it is low when
wrap-around leaves no quiet stretch.

The anchor method finds the true phase. Given a few compounds whose real
second-dimension times are known (internal standards are ideal), it finds
each one in the raw stream, measures when its brightest pulse arrives,
and works out which phase would put that pulse at the known t2. The
answers from the different anchors are combined with a circular mean
(see circstats), after dropping any anchor that disagrees wildly with the
rest. Note that a run-wide shift in second-dimension retention looks
exactly like a phase change, so on shifted samples this method absorbs
that shift too. For aligning to expected positions that is what you want.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import gaussian_filter1d, uniform_filter1d

from . import circstats as cs
from .fold import Fold, FoldSpec
from .run import Run
from .solvent import find_solvent_front


class PhaseEstimationError(RuntimeError):
    pass


@dataclass
class Anchor:
    """A compound with a known second-dimension time.

    mz        : ion(s) used to find it
    t1        : expected first-dimension retention (absolute seconds)
    t2        : expected second-dimension time (seconds, may exceed the period)
    t1_window : half-width of the search window around t1 (seconds)
    """
    name: str
    mz: int | list[int]
    t1: float
    t2: float
    t1_window: float = 12.0


@dataclass
class AnchorHit:
    name: str
    apex_time: float | None
    height: float
    phase: float | None
    used: bool = False


@dataclass
class PhaseEstimate:
    phase: float
    method: str
    spread: float = float("nan")      # seconds; circular SD (anchors) or NaN
    quiet_width: float = float("nan")  # seconds; width of the void (void method)
    hits: list[AnchorHit] = field(default_factory=list)
    profile: np.ndarray | None = None  # smoothed t2 profile (void method)
    contrast: float = float("nan")     # 1 - signal at the cut / typical signal (void method)
    confidence: str = ""               # high, medium or low (void method)
    solvent: tuple[float, float] | None = None   # excluded solvent front, seconds

    def summary(self) -> str:
        s = f"Phase {self.phase:.3f} s ({self.method})"
        if self.method == "anchors":
            used = sum(h.used for h in self.hits)
            s += f", {used}/{len(self.hits)} anchors, spread {self.spread * 1e3:.1f} ms"
        elif self.method == "void":
            s += f", quiet band {self.quiet_width:.2f} s wide, {self.confidence} confidence"
            if self.solvent:
                s += f", solvent front {self.solvent[0]:.0f}-{self.solvent[1]:.0f} s excluded"
        elif self.method == "fallback":
            s += ", borrowed from other runs in the sequence"
        return s


# --------------------------------------------------------------------- void
def _background_corrected(img: np.ndarray, pct: float = 20.0) -> np.ndarray:
    # Bleed and baseline are roughly constant down each column; peaks aren't.
    base = np.percentile(img, pct, axis=0, keepdims=True)
    return np.clip(img - base, 0.0, None)


def _longest_circular_run(mask: np.ndarray) -> tuple[int, int]:
    """Start index and length of the longest run of True, allowing wrap-round."""
    n = mask.size
    if mask.all():
        return 0, n
    if not mask.any():
        return 0, 0
    # Rotate so the array starts on a False, then runs can't straddle the end.
    shift = int(np.flatnonzero(~mask)[0])
    m = np.roll(mask, -shift)
    best_start, best_len, start = 0, 0, None
    for i, v in enumerate(np.append(m, False)):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if i - start > best_len:
                best_start, best_len = start, i - start
            start = None
    return (best_start + shift) % n, best_len


def phase_from_void(run: Run, period: float, t0: float | None = None,
                    n_rows: int | None = None, smooth_s: float = 0.05,
                    quiet_fraction: float = 0.08, window_s: float = 0.3,
                    edge_fraction: float = 0.05, skip_solvent: bool = True) -> PhaseEstimate:
    """Put the cut where the least signal is split across it.

    The cost of a cut is the signal within `window_s` of it, about one peak
    width. The cut goes in the longest quiet stretch, `edge_fraction` of a
    period inside its far edge. Where wrap-around leaves no quiet stretch, it
    still goes at the least damaging place, and a low `contrast` says the
    choice is weak. Columns are weighted so that
    no few very bright ones (a solvent front, say) decide the answer, and a
    detected solvent front is left out altogether."""
    fold0 = Fold(run.times, FoldSpec(period=period, phase=0.0, t0=t0, n_rows=n_rows))
    img = _background_corrected(fold0.tic_image(run))
    solvent = find_solvent_front(run.times, run.tic) if skip_solvent else None
    use = fold0.complete.copy()
    if solvent:
        use &= ~((fold0.t1 + period > solvent[0]) & (fold0.t1 < solvent[1]))
    img = img[:, use]
    col = img.sum(axis=0)
    if not np.any(col > 0):
        raise PhaseEstimationError("no signal to find the phase from")
    # Columns brighter than typical count equally; dimmer ones count less.
    img = img / np.maximum(col, np.median(col[col > 0]))
    profile = img.sum(axis=1)
    sigma = max(1.0, smooth_s / fold0.row_pitch)
    prof = gaussian_filter1d(profile, sigma, mode="wrap")
    size = max(1, int(round(window_s / fold0.row_pitch)))
    cost = uniform_filter1d(prof, size, mode="wrap")
    # Signal scale: the brightest tenth of the profile. The median would be
    # wrong here, since more than half the period is often empty.
    lo, hi = cost.min(), float(np.sort(cost)[-max(1, cost.size // 10):].mean())
    quiet = cost <= lo + quiet_fraction * (hi - lo)
    start, length = _longest_circular_run(quiet)
    # Cut just inside the far edge of the quiet stretch, so signal starts near
    # the bottom of the image; the middle if the stretch is narrow.
    inset = min(length / 2.0, (edge_fraction * period + window_s / 2.0) / fold0.row_pitch)
    cut_row = (start + length - 1 - inset) % fold0.n_rows
    contrast = float(1.0 - cost[int(round(cut_row)) % fold0.n_rows] / hi) if hi > lo else 0.0
    confidence = "high" if contrast >= 0.8 else "medium" if contrast >= 0.5 else "low"
    return PhaseEstimate(phase=float(cut_row * fold0.row_pitch) % period, method="void",
                         quiet_width=(length - 1) * fold0.row_pitch + window_s, profile=prof,
                         contrast=contrast, confidence=confidence, solvent=solvent)


# ------------------------------------------------------------------ anchors
def locate_anchor(run: Run, anchor: Anchor, min_height: float = 50.0) -> tuple[float | None, float]:
    """Apex time and height of the brightest pulse of an anchor's ion(s)
    inside its search window. The apex is refined between scans by fitting
    a parabola to the log of the three highest points (a Gaussian peak is a
    parabola in log space)."""
    lo = np.searchsorted(run.times, anchor.t1 - anchor.t1_window)
    hi = np.searchsorted(run.times, anchor.t1 + anchor.t1_window)
    if hi - lo < 3:
        return None, 0.0
    x = run.xic(anchor.mz)[lo:hi]
    i = int(np.argmax(x))
    height = float(x[i])
    if height < min_height or i == 0 or i == x.size - 1:
        return None, height
    a, b, c = x[i - 1], x[i], x[i + 1]
    if a > 0 and c > 0:
        la, lb, lc = np.log(a), np.log(b), np.log(c)
        den = la - 2 * lb + lc
    else:
        la, lb, lc = a, b, c
        den = a - 2 * b + c
    delta = 0.5 * (la - lc) / den if den != 0 else 0.0
    delta = float(np.clip(delta, -0.5, 0.5))
    t = run.times[lo + i]
    step = run.times[lo + i + 1] - t if delta > 0 else t - run.times[lo + i - 1]
    return float(t + delta * step), height


def phase_from_anchors(run: Run, anchors: list[Anchor], period: float,
                       t0: float | None = None, max_deviation: float = 0.15,
                       min_height: float = 50.0) -> PhaseEstimate:
    t0 = float(run.times[0]) if t0 is None else float(t0)
    hits: list[AnchorHit] = []
    for a in anchors:
        t_obs, h = locate_anchor(run, a, min_height=min_height)
        phi = None if t_obs is None else float(((t_obs - t0) - a.t2) % period)
        hits.append(AnchorHit(a.name, t_obs, h, phi))

    phis = np.array([h.phase for h in hits if h.phase is not None])
    if phis.size == 0:
        raise PhaseEstimationError("none of the anchors were found in this run")

    med = cs.circ_median(phis, period)
    for h in hits:
        if h.phase is not None:
            h.used = abs(float(cs.circ_diff(h.phase, med, period))) <= max_deviation
    kept = np.array([h.phase for h in hits if h.used])
    return PhaseEstimate(phase=cs.circ_mean(kept, period), method="anchors",
                         spread=cs.circ_std(kept, period) if kept.size > 1 else float("nan"),
                         hits=hits)


def auto_phase(run: Run, period: float, anchors: list[Anchor] | None = None,
               min_anchors: int = 2, fallback: float | None = None, **kw) -> PhaseEstimate:
    """Anchors when enough are found (a table with a single anchor needs just
    that one); otherwise a known fallback phase (for example from other runs
    in the same sequence); otherwise the void method."""
    if anchors:
        try:
            est = phase_from_anchors(run, anchors, period, t0=kw.get("t0"))
            if sum(h.used for h in est.hits) >= min(min_anchors, len(anchors)):
                return est
        except PhaseEstimationError:
            pass
    if fallback is not None:
        return PhaseEstimate(phase=float(fallback) % period, method="fallback")
    return phase_from_void(run, period, **kw)


def phase_for_sequence(runs: list[Run], period: float, anchors: list[Anchor] | None = None,
                       min_anchors: int = 2) -> list[PhaseEstimate]:
    """Phase for every run in a sequence.

    The phase is a property of the instrument timing, so it barely changes
    between runs. Runs where the anchors work get their own estimate; runs
    where they don't (tube blanks, say) borrow the circular median of the
    others. The void method is used only if no run has usable anchors.
    """
    first = [auto_phase(r, period, anchors, min_anchors, fallback=float("nan")) for r in runs]
    good = [e.phase for e in first if e.method == "anchors"]
    if not good:
        voids = []
        for r in runs:
            try:
                voids.append(phase_from_void(r, period))
            except PhaseEstimationError:
                voids.append(None)     # e.g. a run with no signal
        solid = [v.phase for v in voids if v is not None and v.confidence != "low"]
        if not solid:
            return [v if v is not None else PhaseEstimate(phase=0.0, method="fallback") for v in voids]
        shared = cs.circ_median(np.array(solid), period)
        return [v if v is not None and v.confidence != "low"
                else PhaseEstimate(phase=shared, method="fallback") for v in voids]
    shared = cs.circ_median(np.array(good), period)
    return [e if e.method == "anchors" else PhaseEstimate(phase=shared, method="fallback")
            for e in first]
