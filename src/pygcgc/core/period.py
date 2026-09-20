"""Checking the modulation period.

If the period used for folding is slightly wrong, every modulation is cut
a little early or late, and the error builds up: peaks drift steadily up or
down the image as the run goes on. With a 4 s period, an error of 4 ms is
invisible at the start but has moved peaks by nearly 2 s after 450
modulations.

Tracking individual peaks cannot give the period: a genuine drift of
second-dimension retention through the run looks exactly like a period
error. The stream as a whole can. Every eluting compound, and any baseline
disturbance from the modulator, repeats once per modulation, so the total
signal is a pulse train at the modulation frequency. That is the period
folding needs, whatever the retention does.

`estimate_period`: finds that periodicity with no prior knowledge. A
spectrum of the compressed, high-passed TIC gives a coarse answer, using a
comb that rewards power at f, 2f, 3f and penalises power halfway between,
so half or double the true period is not mistaken for it. Then the run is
folded at trial periods and the period that gives the sharpest average
profile is kept, first coarsely and then finely. Repeating the fine search
on each fifth of the run separately gives a standard error. Discrete
compounds are visible for only a few modulations each, so they pin the
period down to a few milliseconds at best; structure locked to the
modulator across the whole run (column bleed, an unresolved hump) does much
better. The result therefore reports both whether a period is established
(`confidence`) and how far peaks could drift over the run (`drift`). A
solvent front is masked first. A
genuine drift of one compound's t2 is still invisible to this method; only
anchors can separate that from a period error.

`period_from_anchors`: each anchor gives its own phase estimate. If the
period is off by d, the estimate from an anchor in modulation k is off by
k × d, so the estimates fall on a sloping line against modulation number.
Fitting that line gives the period correction (slope) and the phase
(intercept), with a standard error on the correction so you can see
whether it is real. Because the anchors' true t2 values are known, a
genuine t2 drift through the run doesn't fool it. The fit first searches
for the slope that makes the corrected estimates agree best on the circle,
which copes with large errors that wrap round, then assigns each anchor
its exact modulation number and fits an ordinary straight line.

`period_from_points` is the same fit for points you click on the image
along a series you expect to be flat (for example the n-alkanes), where
the true t2 values are unknown but assumed equal.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import uniform_filter1d

from . import circstats as cs
from .phase import Anchor, PhaseEstimationError, locate_anchor
from .run import Run
from .solvent import find_solvent_front


PRECISE_DRIFT = 0.25   # seconds; PeriodEstimate.drift above this means the period should be checked against anchors


@dataclass
class PeriodEstimate:
    period: float
    method: str
    phase: float | None = None
    residual_sd: float = float("nan")   # seconds
    correction_se: float = float("nan")  # standard error of the period correction (s)
    nominal: float = float("nan")
    n_points: int = 0
    confidence: str = ""                 # high, medium or low (method "auto" only)
    z: float = float("nan")              # spectral comb score against the rest of the search band
    drift: float = float("nan")          # seconds: one standard error of the period times the modulations in the run
    raw_period: float = float("nan")     # before snapping to a round value
    acf_period: float = float("nan")     # independent autocorrelation estimate
    solvent: tuple[float, float] | None = None   # masked solvent front, seconds
    prior_agrees: bool | None = None     # None when no prior was given

    @property
    def correction(self) -> float:
        return self.period - self.nominal

    @property
    def precise(self) -> bool:
        """False if peaks at the end of the run could be displaced by more than
        a quarter of a second by the remaining period uncertainty."""
        return bool(np.isfinite(self.drift) and self.drift <= PRECISE_DRIFT)

    @property
    def significant(self) -> bool:
        """True if the correction is more than three standard errors from zero."""
        return bool(np.isfinite(self.correction_se) and
                    abs(self.correction) > 3 * self.correction_se)

    def summary(self) -> str:
        if self.method == "auto":
            s = f"Period {self.period:.4f} s (estimated from the run, {self.confidence} confidence"
            if np.isfinite(self.correction_se):
                s += f", ± {self.correction_se * 1e3:.2f} ms"
            s += f", {self.n_points} modulations)"
            if not self.precise:
                s += (f"; peaks could drift by about {self.drift:.1f} s over the run, "
                      "so check the period against anchors")
            if self.solvent:
                s += f"; solvent front {self.solvent[0]:.0f}-{self.solvent[1]:.0f} s excluded"
            if self.prior_agrees is False:
                s += f"; DISAGREES with the expected {self.nominal:.4f} s"
            return s
        s = f"Period {self.period:.5f} s ({self.method})"
        if self.phase is not None:
            s += f", phase {self.phase:.3f} s"
        if self.n_points:
            s += f", {self.n_points} points, residual SD {self.residual_sd * 1e3:.1f} ms"
        if np.isfinite(self.correction_se):
            s += (f"; correction {self.correction * 1e3:+.3f} ± {self.correction_se * 1e3:.3f} ms"
                  f" ({'significant' if self.significant else 'not significant'})")
        return s


def _fit_drift(abs_times: np.ndarray, t2: np.ndarray, period: float, t0: float,
               rel_range: float, weights: np.ndarray | None):
    """Joint fit of phase and period from points with known t2.

    Model: t_obs = t0 + phase + k * true_period + t2, where k is the
    modulation number. Returns (true period, phase, residual SD).
    """
    rel = abs_times - t0 - t2
    w = np.ones_like(rel) if weights is None else weights
    k = np.floor(rel / period)
    phi = rel - k * period                      # per-point phase at the nominal period

    # Search period errors d for the best agreement on the circle.
    d_max = rel_range * period
    grid = np.linspace(-d_max, d_max, 4001)
    R = np.array([cs.resultant_length(phi - k * d, period, w) for d in grid])
    d0 = grid[int(np.argmax(R))]
    centre = cs.circ_mean(phi - k * d0, period, w)

    # With a good first guess, assign each point its modulation number exactly,
    # after which the problem is an ordinary straight-line fit.
    k = np.round((rel - centre) / (period + d0))
    y = rel - k * period                        # = phase + k * d
    A = np.vstack([np.ones_like(k), k]).T
    sw = np.sqrt(w)
    (phase, d), *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
    resid = y - (phase + d * k)
    true_period = float(period + d)
    n = rel.size
    if n > 2:
        sd = float(np.sqrt(np.sum(w * resid ** 2) / np.sum(w) * n / (n - 2)))
        kbar = np.average(k, weights=w)
        sxx = float(np.sum(w * (k - kbar) ** 2) / np.mean(w))
        se = sd / np.sqrt(sxx) if sxx > 0 else float("nan")
    else:
        sd = se = float("nan")
    return true_period, float(phase % true_period), sd, float(se)


def period_from_anchors(run: Run, anchors: list[Anchor], nominal: float,
                        t0: float | None = None, rel_range: float = 0.005,
                        min_height: float = 50.0) -> PeriodEstimate:
    t0 = float(run.times[0]) if t0 is None else float(t0)
    times, t2s, heights = [], [], []
    for a in anchors:
        t_obs, h = locate_anchor(run, a, min_height=min_height)
        if t_obs is not None:
            times.append(t_obs), t2s.append(a.t2), heights.append(h)
    if len(times) < 3:
        raise PhaseEstimationError(f"need at least 3 anchors to fit the period, found {len(times)}")
    times, t2s = np.array(times), np.array(t2s)
    period, phase, sd, se = _fit_drift(times, t2s, nominal, t0, rel_range, None)
    return PeriodEstimate(period=period, method="anchors", phase=phase, residual_sd=sd,
                          correction_se=se, nominal=nominal, n_points=len(times))


def period_from_points(abs_times, nominal: float, t0: float,
                       rel_range: float = 0.005) -> PeriodEstimate:
    """Period from clicked points along a series assumed flat in true t2.
    Only the period is meaningful here; the phase is not reported."""
    abs_times = np.asarray(abs_times, dtype=np.float64)
    if abs_times.size < 2:
        raise ValueError("need at least two points")
    period, _, sd, se = _fit_drift(abs_times, np.zeros_like(abs_times), nominal, t0, rel_range, None)
    return PeriodEstimate(period=period, method="points", residual_sd=sd, correction_se=se,
                          nominal=nominal, n_points=abs_times.size)


# ------------------------------------------------------------------- blind
def estimate_period_for_runs(runs: list[Run], max_runs: int = 3, **kw) -> tuple[PeriodEstimate | None, list[PeriodEstimate]]:
    """Blind estimate from a few evenly spaced runs of a sequence (all of them
    would be slow, and they share one modulator). Returns the most precise
    estimate that is not low confidence (None if there is none), and every
    estimate that could be made."""
    pick = np.unique(np.linspace(0, len(runs) - 1, min(max_runs, len(runs))).astype(int))
    found = []
    for i in pick:
        try:
            found.append(estimate_period(runs[i], **kw))
        except PhaseEstimationError:
            pass                                    # no signal, or too short
    usable = [f for f in found if f.confidence != "low"]
    best = min(usable, key=lambda f: f.drift if np.isfinite(f.drift) else np.inf) if usable else None
    return best, found


def _fold_stat(tn: np.ndarray, y: np.ndarray, period: float, nbins: int) -> float:
    """Between-bin sum of squares of y folded at `period`. Large when the
    average profile is sharp, which happens only at the right period."""
    ph = tn / period
    ph -= np.floor(ph)
    idx = (ph * nbins).astype(np.intp)
    n = np.bincount(idx, minlength=nbins)
    s = np.bincount(idx, weights=y, minlength=nbins)
    return float(np.sum(s ** 2 / np.maximum(n, 1)) - s.sum() ** 2 / y.size)


def _fold_scan(tn, y, centre, half_range, step, nbins) -> tuple[float, bool, float]:
    """Best period on a grid around `centre`, refined by a parabola.
    Also says whether the best point sat on the edge of the grid, and how
    far the best fold statistic stands above the rest of the grid (robust
    z-score: large for one sharp peak, small for a flat noisy landscape)."""
    grid = centre + np.arange(-half_range, half_range + 0.5 * step, step)
    q = np.array([_fold_stat(tn, y, p, nbins) for p in grid])
    i = int(np.argmax(q))
    edge = i == 0 or i == grid.size - 1
    best = float(grid[i])
    if not edge:
        a, b, c = q[i - 1], q[i], q[i + 1]
        den = a - 2 * b + c
        if den < 0:
            best += 0.5 * (a - c) / den * step
    mad = float(np.median(np.abs(q - np.median(q))))
    sharp = float((q[i] - np.median(q)) / (1.4826 * mad)) if mad > 0 else float("inf")
    return best, edge, sharp


def _parabolic(values: np.ndarray, i: int) -> float:
    """Fractional index of the peak of `values` near integer index i."""
    if i <= 0 or i >= values.size - 1:
        return float(i)
    a, b, c = values[i - 1], values[i], values[i + 1]
    den = a - 2 * b + c
    return float(i + (0.5 * (a - c) / den if den != 0 else 0.0))


def estimate_period(run: Run, prior: float | None = None, prior_tol: float = 0.03,
                    search: tuple[float, float] = (1.0, 20.0), n_blocks: int = 5,
                    block_half: float = 0.002, snap: float = 0.01) -> PeriodEstimate:
    """Modulation period from the run itself. With `prior`, the search is
    limited to prior ± prior_tol and the result says whether it agrees.

    The result is snapped to a multiple of `snap` seconds (methods use round
    periods) only if that is within two standard errors of the estimate."""
    pmin, pmax = search
    pick_min, pick_max = pmin, pmax
    if prior is not None:
        pick_min, pick_max = max(pmin, prior * (1 - prior_tol)), min(pmax, prior * (1 + prior_tol))

    times, tic = run.times, run.tic
    active = np.flatnonzero(tic > 0)
    if active.size < 100:
        raise PhaseEstimationError("no signal in this run")
    dt = float(np.median(np.diff(times)))
    solvent = find_solvent_front(times, tic)

    # Trim empty stretches at either end, compress the dynamic range so a few
    # huge peaks do not dominate, and remove slow trends (window: three of the
    # longest periods searched) so only structure near the modulation remains.
    sl = slice(active[0], active[-1] + 1)
    t, x = times[sl], tic[sl]
    y = np.log1p(np.clip(x, 0.0, None))
    y -= uniform_filter1d(y, max(3, int(round(3 * search[1] / dt))), mode="nearest")
    keep = np.ones(t.size, dtype=bool)
    if solvent:
        keep = ~((t >= solvent[0]) & (t <= solvent[1]))
    y[~keep] = 0.0
    y[keep] -= y[keep].mean()
    tn = t - t[0]
    span = float(tn[-1])
    if span < 12 * pmin:
        raise PhaseEstimationError(f"run too short to find the period ({span:.0f} s)")

    # Coarse: spectrum on a regular grid, scored with the alternating comb.
    n = t.size
    yg = np.interp(t[0] + dt * np.arange(n), t, y) * np.hanning(n)
    nfft = 1 << int(np.ceil(np.log2(4 * n)))
    spec = np.fft.rfft(yg, nfft)
    power = np.abs(spec) ** 2
    freq = np.fft.rfftfreq(nfft, dt)
    n_harm = int(max(2, min(6, np.floor(0.9 * (0.5 / dt) * pmin))))
    in_band = (freq >= 1 / pmax) & (freq <= 1 / pmin)
    cand = freq[in_band]
    score = np.zeros(cand.size)
    for h in range(1, n_harm + 1):
        score += np.interp(h * cand, freq, power) - np.interp((h - 0.5) * cand, freq, power)
    # The z-score always uses the whole search band; a prior only limits
    # where the peak may be picked.
    allowed = (cand >= 1 / pick_max) & (cand <= 1 / pick_min)
    i_best = int(np.flatnonzero(allowed)[np.argmax(score[allowed])])
    f0 = cand[0] + _parabolic(score, i_best) * (cand[1] - cand[0])
    p0 = 1.0 / f0
    if span < 10 * p0:
        raise PhaseEstimationError(f"run too short to find the period ({span / p0:.0f} modulations)")
    mad = float(np.median(np.abs(score - np.median(score))))
    z = float((score[i_best] - np.median(score)) / (1.4826 * mad)) if mad > 0 else float("inf")

    # Independent check: the strongest autocorrelation peak in the picked range.
    acf = np.fft.irfft(power)[: n // 2]
    lags = np.arange(acf.size) * dt
    lo, hi = np.searchsorted(lags, pick_min), np.searchsorted(lags, pick_max)
    j = lo + int(np.argmax(acf[lo:hi])) if hi > lo else lo
    acf_period = float(lags[j] + (_parabolic(acf, j) - j) * dt)
    near = (lags[lo:hi] > 0.98 * p0) & (lags[lo:hi] < 1.02 * p0)
    acf_agrees = abs(acf_period - p0) < 0.02 * p0 or (
        near.any() and acf[lo:hi][near].max() >= 0.5 * acf[j])

    # Fine: fold at trial periods, first in wide bins over the coarse
    # uncertainty, then in narrower ones. A period error e smears a feature
    # by e * span / period seconds, so each step moves the fold by a quarter
    # of a bin.
    n_rows = max(1, int(round(p0 / dt)))
    bins = sorted({32, min(128, max(32, n_rows // 2)), max(64, n_rows // 2)})
    dp_fft = p0 ** 2 / span
    p, half = p0, 1.5 * dp_fft
    for nb in bins:
        step = 0.25 * (p0 / nb) * p0 / span
        p, _, _ = _fold_scan(tn[keep], y[keep], p, half, step, nb)
        half = 4 * step
    period_raw = p

    # Standard error: estimate the period again from each block of the run on
    # its own. Blocks that share most of their data (leave-one-out) would agree
    # even where the period is poorly constrained; separate blocks do not.
    # Discrete compounds are only visible for a few modulations and carry
    # little information about the period; what constrains it is structure
    # that stays locked to the modulator across the whole run.
    nb_b = bins[1] if len(bins) > 1 else bins[0]
    block_span = span / n_blocks
    step_b = 0.25 * (period_raw / nb_b) * period_raw / block_span
    edges = np.searchsorted(tn, np.linspace(0, span, n_blocks + 1))
    edges[-1] = t.size
    block_est = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = keep[a:b]
        if m.sum() < 100:
            continue
        pj, _, _ = _fold_scan(tn[a:b][m], y[a:b][m], period_raw, block_half * period_raw,
                              step_b, nb_b)
        block_est.append(pj)
    block_est = np.array(block_est)
    se = float(block_est.std(ddof=1) / np.sqrt(block_est.size)) if block_est.size > 2 else float("nan")

    # Folding the whole run at once lets different compounds overlap by chance
    # at some trial periods, which biases the answer by a few milliseconds when
    # nothing stays locked to the modulator all run. Blocks are free of that
    # bias, so if the whole-run answer disagrees with them, trust them.
    if np.isfinite(se) and abs(period_raw - np.median(block_est)) > 2 * se:
        period_raw = float(np.median(block_est))

    period = period_raw
    snapped = round(period_raw / snap) * snap
    if np.isfinite(se) and abs(snapped - period_raw) <= max(2 * se, 5e-4):
        period = snapped

    agrees_prior = None if prior is None else bool(abs(period - prior) <= max(3 * se, 0.001))
    n_mod = int(span / period)
    # Confidence says whether a modulation period is established at all: a clear
    # spectral line, an independent autocorrelation peak at the same place, and
    # enough modulations. How finely it is pinned down is a separate matter,
    # given by the drift: how far peaks at the end of the run could sit from
    # where they should be (standard error times number of modulations).
    drift = se * n_mod
    # A prior must not steer the search onto a weak alias (a wrong prior can
    # catch the second and fourth harmonics of the real line): the picked peak
    # has to be comparable with the strongest one anywhere in the band.
    dominant = score[i_best] >= 0.5 * score.max()
    if z >= 8 and acf_agrees and dominant and n_mod >= 30:
        confidence = "high"
    elif z >= 5 and acf_agrees and dominant and n_mod >= 15:
        confidence = "medium"
    else:
        confidence = "low"
    return PeriodEstimate(period=float(period), method="auto", correction_se=se,
                          nominal=float("nan") if prior is None else float(prior),
                          n_points=n_mod, confidence=confidence, z=z, drift=drift,
                          raw_period=float(period_raw),
                          acf_period=acf_period, solvent=solvent, prior_agrees=agrees_prior)
