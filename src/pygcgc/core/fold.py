"""Folding: turning the raw scan stream into a 2D chromatogram.

The detector records one long stream. Every modulation period the
modulator releases a new slice of the first-column effluent onto the
second column, so the stream is really a sequence of short second-dimension
chromatograms laid end to end. Folding cuts the stream at each modulation
and stacks the pieces side by side: column k of the image is modulation k
(first-dimension time), and the vertical position within a column is the
time since that modulation started (second-dimension time).

Two details matter.

Phase. The modulator's cycle rarely starts exactly when the acquisition
does. `phase` is the offset, in seconds, between the start of acquisition
(`t0`) and the first modulation boundary. Get it wrong and peaks sit too
high or too low, and some wrap over the top of the image.

Fractional folding. If the period is not a whole number of scans (4 s at
49.98 Hz, say), cutting fixed blocks of scans makes peaks creep diagonally
through the run. So each scan is placed by its time: its position within
its modulation, in rows, usually falls between two rows, and its intensity
is shared between those two rows in proportion to closeness. Totals are
conserved exactly.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
import scipy.sparse as sp

from .run import Run

_SNAP = 1e-6  # row positions within a millionth of a row are treated as exact


@dataclass(frozen=True)
class FoldSpec:
    """How to fold a run.

    period : modulation period in seconds.
    phase  : seconds from t0 to the first modulation boundary.
    t0     : reference time; defaults to the first scan time.
    n_rows : rows per modulation; defaults to period / scan interval.
    """
    period: float
    phase: float = 0.0
    t0: float | None = None
    n_rows: int | None = None

    def resolved(self, times: np.ndarray) -> "FoldSpec":
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or not times.size or not np.isfinite(times).all():
            raise ValueError("folding needs finite scan times in a non-empty one-dimensional array")
        if not np.isfinite(self.period) or self.period <= 0:
            raise ValueError("modulation period must be finite and positive")
        if not np.isfinite(self.phase):
            raise ValueError("modulation phase must be finite")
        t0 = float(times[0]) if self.t0 is None else float(self.t0)
        if not np.isfinite(t0):
            raise ValueError("fold time origin must be finite")
        n_rows = self.n_rows
        if n_rows is None:
            if times.size < 2 or np.any(np.diff(times) <= 0):
                raise ValueError("estimating fold rows needs at least two increasing scan times")
            n_rows = int(round(self.period / float(np.median(np.diff(times)))))
        if not np.isfinite(n_rows) or n_rows < 1 or int(n_rows) != n_rows:
            raise ValueError("fold rows must be a positive whole number; check period and scan spacing")
        return replace(self, t0=t0, n_rows=int(n_rows))


class Fold:
    """Where every scan lands in the folded image, for one run and one FoldSpec.

    Building a Fold is cheap (it only works on scan times), so changing the
    phase or period just means building a new one.
    """

    def __init__(self, times: np.ndarray, spec: FoldSpec):
        times = np.asarray(times, dtype=np.float64)
        spec = spec.resolved(times)
        self.spec = spec
        self.times = times
        P, n_rows = spec.period, spec.n_rows
        self.row_pitch = P / n_rows

        rel = times - spec.t0 - spec.phase
        col = np.floor(rel / P).astype(np.int64)
        pos = (rel / P - col) * n_rows
        pos = np.round(pos / _SNAP) * _SNAP
        row = np.floor(pos).astype(np.int64)
        w = pos - row
        w[w < _SNAP] = 0.0
        # Carry rows that rounded up to n_rows into the next modulation.
        col = col + row // n_rows
        row = row % n_rows
        row_up = row + 1
        col_up = col + row_up // n_rows
        row_up = row_up % n_rows

        valid = col >= 0
        has_up = valid & (w > 0)
        if not np.any(valid):
            raise ValueError("no scans fall after the first modulation boundary")
        self.n_cols = int(max(col[valid].max(), col_up[has_up].max() if has_up.any() else 0)) + 1

        self._lo = np.where(valid, col * n_rows + row, -1)
        self._up = np.where(has_up, col_up * n_rows + row_up, -1)
        self._w = w
        self._valid = valid
        self._has_up = has_up

    # ------------------------------------------------------------------ axes
    @property
    def n_rows(self) -> int:
        return self.spec.n_rows

    @property
    def shape(self) -> tuple[int, int]:
        return (self.n_rows, self.n_cols)

    @property
    def t1(self) -> np.ndarray:
        """Start time (s) of each modulation, i.e. the first-dimension axis."""
        s = self.spec
        return s.t0 + s.phase + np.arange(self.n_cols) * s.period

    @property
    def t2(self) -> np.ndarray:
        """Second-dimension time (s) of each row."""
        return np.arange(self.n_rows) * self.row_pitch

    @property
    def complete(self) -> np.ndarray:
        """True for modulations fully covered by the acquisition."""
        dt = float(np.median(np.diff(self.times)))
        start = self.t1
        return (start >= self.times[0] - 0.5 * dt) & (start + self.spec.period <= self.times[-1] + dt)

    # ------------------------------------------------------------ coordinates
    def to_fold(self, t_abs):
        """Absolute time(s) -> (column index, t2 in seconds)."""
        s = self.spec
        rel = np.asarray(t_abs, dtype=np.float64) - s.t0 - s.phase
        col = np.floor(rel / s.period).astype(np.int64)
        return col, rel - col * s.period

    def to_time(self, col, t2):
        """(column index, t2 in seconds) -> absolute time. t2 may exceed the
        period, which is how wrapped peaks are described."""
        s = self.spec
        return s.t0 + s.phase + np.asarray(col) * s.period + np.asarray(t2)

    # ---------------------------------------------------------------- images
    def image(self, signal: np.ndarray) -> np.ndarray:
        """Fold one value per scan into an (n_rows, n_cols) image."""
        signal = np.asarray(signal, dtype=np.float64)
        if signal.shape != self.times.shape:
            raise ValueError("signal must have one value per scan")
        size = self.n_rows * self.n_cols
        v, u = self._valid, self._has_up
        img = np.bincount(self._lo[v], weights=signal[v] * (1.0 - self._w[v]), minlength=size)
        img += np.bincount(self._up[u], weights=signal[u] * self._w[u], minlength=size)
        return img.reshape(self.n_cols, self.n_rows).T

    def tic_image(self, run: Run) -> np.ndarray:
        return self.image(run.tic)

    def ion_image(self, run: Run, mz) -> np.ndarray:
        """Image of one m/z, or the sum of several."""
        return self.image(run.xic(mz))

    def cube(self, run: Run, mz=None, dtype=np.float32) -> np.ndarray:
        """(n_rows, n_cols, n_ions) array for the chosen ions (all if None).
        All ions for a whole run is roughly 170 MB, so pass a subset when you can."""
        if run.n_scans != self.times.size:
            raise ValueError("run does not match this fold")
        cols = np.arange(run.n_mz) if mz is None else run.mz_index(mz)
        sub = run.data[:, cols] if mz is not None else run.data
        coo = sp.coo_matrix(sub)
        n_ions = cols.size
        size = self.n_rows * self.n_cols * n_ions
        out = np.zeros(size, dtype=dtype)
        scan, ion, val = coo.row, coo.col, coo.data.astype(np.float64)
        v = self._valid[scan]
        np.add.at(out, self._lo[scan[v]] * n_ions + ion[v], (val[v] * (1.0 - self._w[scan[v]])).astype(dtype))
        u = self._has_up[scan]
        np.add.at(out, self._up[scan[u]] * n_ions + ion[u], (val[u] * self._w[scan[u]]).astype(dtype))
        return np.ascontiguousarray(out.reshape(self.n_cols, self.n_rows, n_ions).transpose(1, 0, 2))

    # -------------------------------------------------------------- wrapping
    def extend(self, img: np.ndarray, factor: float = 1.3) -> tuple[np.ndarray, np.ndarray]:
        """Extend the t2 axis beyond one period so wrapped peaks look whole.

        A compound released in modulation k with t2 longer than the period
        elutes during modulation k+1, near its bottom. The extra band on top
        of column k is therefore the bottom of column k+1.
        Returns the extended image and its t2 axis.
        """
        extra = int(round((factor - 1.0) * self.n_rows))
        if extra <= 0:
            return img, self.t2
        extra = min(extra, self.n_rows)
        top = np.zeros((extra,) + img.shape[1:], dtype=img.dtype)
        top[:, :-1] = img[:extra, 1:]
        ext = np.concatenate([img, top], axis=0)
        t2 = np.arange(ext.shape[0]) * self.row_pitch
        return ext, t2


def fold(run: Run, period: float, phase: float = 0.0, **kw) -> Fold:
    """Shortcut: Fold for a run."""
    return Fold(run.times, FoldSpec(period=period, phase=phase, **kw))
