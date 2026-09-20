"""Scan timing check.

Folding assumes scans arrive at a steady rate. This module checks that.
It places every scan on the nearest slot of an evenly spaced grid and
fits a straight line of scan time against slot number: the slope is the
real scan interval, and the leftovers (residuals) show jitter. A skipped
slot is a gap where a scan was dropped; `fill_gaps` inserts empty scans there so the raw stream
and the fold stay in step.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp

from .run import Run

JITTER_WARN_FRACTION = 0.1  # warn if RMS jitter exceeds 10% of the interval


@dataclass
class TimingReport:
    interval: float                 # seconds per scan, from the fit
    rate: float                     # scans per second
    rms_jitter: float               # seconds
    max_residual: float             # seconds, largest absolute residual
    gaps: list[tuple[int, int]] = field(default_factory=list)  # (scan index before gap, scans missing)
    n_scans: int = 0
    collisions: int = 0   # scans that land on an already-used slot (very irregular timing)

    @property
    def missing_scans(self) -> int:
        return sum(n for _, n in self.gaps)

    @property
    def ok(self) -> bool:
        return (not self.gaps and not self.collisions
                and self.rms_jitter <= JITTER_WARN_FRACTION * self.interval)

    def summary(self) -> str:
        lines = [
            f"Scans: {self.n_scans:,}",
            f"Rate: {self.rate:.4f} Hz (interval {self.interval * 1e3:.4f} ms)",
            f"Jitter: {self.rms_jitter * 1e3:.3f} ms RMS, worst {self.max_residual * 1e3:.3f} ms",
        ]
        if self.gaps:
            lines.append(f"Gaps: {len(self.gaps)} ({self.missing_scans} scans missing)")
        else:
            lines.append("Gaps: none")
        if self.collisions:
            lines.append(f"Irregular: {self.collisions} scans too close to their neighbours "
                         "to place on a regular grid")
        lines.append("Timing OK" if self.ok else "Timing needs attention")
        return "\n".join(lines)


def check_timing(times: np.ndarray) -> TimingReport:
    t = np.asarray(times, dtype=np.float64)
    n = t.size
    if n < 3:
        raise ValueError("need at least three scans to check timing")

    # Put every scan on its nearest slot of an ideal, evenly spaced grid.
    # Dropped scans show up as skipped slots. The interval has to be accurate
    # enough that its error, added up over the run, stays well under half a
    # scan. So we start with a rough interval (median span of 10 scans, which
    # jitter and the odd gap barely affect), fit a line to the first few
    # hundred scans, and keep extending the stretch fourfold, refitting each
    # time, until it covers the whole run.
    lag = min(10, n - 1)
    slope = float(np.median((t[lag:] - t[:-lag]) / lag))
    intercept = float(t[0])
    end = min(n, 500)
    while True:
        slot = np.rint((t[:end] - intercept) / slope)
        slope, intercept = np.polyfit(slot, t[:end], 1)
        if end == n:
            break
        end = min(n, end * 4)
    slot = np.rint((t - intercept) / slope).astype(np.int64)
    slot -= slot[0]
    jumps = np.diff(slot)
    collisions = int(np.sum(jumps < 1))
    slope, intercept = np.polyfit(slot, t, 1)
    resid = t - (intercept + slope * slot)
    gaps = [(int(i), int(jumps[i] - 1)) for i in np.flatnonzero(jumps > 1)]

    return TimingReport(interval=float(slope), rate=float(1.0 / slope),
                        rms_jitter=float(np.sqrt(np.mean(resid ** 2))),
                        max_residual=float(np.max(np.abs(resid))),
                        gaps=gaps, n_scans=n, collisions=collisions)


def fill_gaps(run: Run, report: TimingReport | None = None) -> Run:
    """Return a run with empty scans inserted wherever scans were dropped.
    Existing scans keep their recorded times; inserted ones are spaced
    evenly across the gap."""
    report = report or check_timing(run.times)
    if report.collisions:
        raise ValueError("timing too irregular to fill gaps safely; inspect the file")
    if not report.gaps:
        return run
    n_old = run.n_scans
    extra = np.zeros(n_old, dtype=np.int64)
    for i, missing in report.gaps:
        extra[i + 1:] += missing
    new_index = np.arange(n_old) + extra
    n_new = n_old + report.missing_scans

    times = np.empty(n_new)
    times[new_index] = run.times
    filled = np.zeros(n_new, dtype=bool)
    filled[new_index] = True
    for i, missing in report.gaps:
        a, b = run.times[i], run.times[i + 1]
        times[new_index[i] + 1:new_index[i] + 1 + missing] = np.linspace(a, b, missing + 2)[1:-1]

    # Re-point each old row at its new row; inserted rows stay empty.
    coo = run.data.tocoo()
    data = sp.csr_matrix((coo.data, (new_index[coo.row], coo.col)), shape=(n_new, run.n_mz))
    out = Run(times=times, mz=run.mz, data=data, name=run.name, source=run.source,
              meta=dict(run.meta))
    out.meta["inserted_scans"] = int(report.missing_scans)
    return out
