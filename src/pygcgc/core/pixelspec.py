"""The mass spectrum under one clicked point of the 2D chromatogram.

In plain terms: to use a peak you clicked as a reference (a siloxane
spectrum, an alkane's qualifier ratios), we need its mass spectrum with the
chemical background taken off. The spectrum of the clicked pixel is summed
over a few neighbouring rows (a peak is spread over several rows of its
modulation), and the spectrum of that same modulation's own baseline, the
median over every row of the column, is subtracted. The baseline is the
right thing to remove here because column bleed, air and carrier-gas ions sit
at a steady level along the whole modulation while the clicked peak occupies
only a few rows of it.
"""

from __future__ import annotations

import numpy as np

from .fold import Fold
from .run import Run


def spectrum_at(run: Run, fold: Fold, t1: float, t2: float, half_rows: int = 2,
                subtract_baseline: bool = True) -> np.ndarray:
    """Dense spectrum (length `run.n_mz`) at a point. `t1` is the absolute
    time of the modulation's middle (the convention `Fold.t1 + period / 2`
    used everywhere else), `t2` the second-dimension time in seconds.

    The rows within `half_rows` of the point are summed; with
    `subtract_baseline`, the same number of rows' worth of the modulation's
    median spectrum is removed and negatives are clipped to zero."""
    col = int(np.floor((t1 - fold.spec.t0 - fold.spec.phase) / fold.spec.period))
    start = fold.spec.t0 + fold.spec.phase + col * fold.spec.period
    lo = np.searchsorted(run.times, start)
    hi = np.searchsorted(run.times, start + fold.spec.period)
    if hi - lo < 3:
        return np.zeros(run.n_mz)
    block = run.data[lo:hi].toarray().astype(np.float64)
    tt = run.times[lo:hi] - start
    near = np.abs(tt - t2) <= (half_rows + 0.5) * fold.row_pitch
    if not near.any():
        return np.zeros(run.n_mz)
    spec = block[near].sum(axis=0)
    if subtract_baseline:
        spec = spec - np.median(block, axis=0) * near.sum()
    return np.clip(spec, 0.0, None)


def ratios_to(spectrum: np.ndarray, mz: np.ndarray, quant: int, ions) -> tuple[float, ...]:
    """Each of `ions`' intensity relative to `quant` (0 where `quant` is 0)."""
    base = float(spectrum[int(quant) - int(mz[0])])
    out = []
    for m in ions:
        j = int(m) - int(mz[0])
        v = float(spectrum[j]) if 0 <= j < spectrum.size else 0.0
        out.append(v / base if base > 0 else 0.0)
    return tuple(out)
