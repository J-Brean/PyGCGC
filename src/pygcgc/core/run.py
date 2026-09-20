"""The Run: one acquisition held in memory.

In plain terms, a run is a long list of mass spectra, one per scan, each
stamped with the time it was recorded. Most (scan, m/z) cells are empty
because the instrument only stores signals above a threshold, so the data
are kept as a sparse matrix with scans as rows and unit m/z as columns.
A 31-minute run at 50 Hz has about a million stored values, which is
roughly 12 MB this way, against about 170 MB if every cell were stored.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path
from typing import Iterable

import numpy as np
import scipy.sparse as sp


@dataclass(eq=False)   # runs compare by identity; comparing arrays is meaningless here
class Run:
    """One GC×GC-MS acquisition.

    Attributes
    ----------
    times:
        Scan acquisition times in seconds, shape (n_scans,).
    mz:
        Unit m/z axis, shape (n_mz,), increasing by one.
    data:
        Sparse intensity matrix, shape (n_scans, n_mz), CSR, float32.
    name:
        Short name, usually the file stem.
    source:
        Path of the file the run was read from, if any.
    meta:
        Anything else worth keeping from the file header.
    """

    times: np.ndarray
    mz: np.ndarray
    data: sp.csr_matrix
    name: str = ""
    source: Path | None = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.times = np.asarray(self.times, dtype=np.float64)
        self.mz = np.asarray(self.mz, dtype=np.int64)
        if not (sp.issparse(self.data) and self.data.format == "csr"):
            self.data = sp.csr_matrix(self.data)
        if self.data.dtype != np.float32:
            self.data = self.data.astype(np.float32)
        n_scans, n_mz = self.data.shape
        if n_scans != self.times.size:
            raise ValueError(f"{n_scans} rows in data but {self.times.size} scan times")
        if n_mz != self.mz.size:
            raise ValueError(f"{n_mz} columns in data but {self.mz.size} m/z values")
        if self.mz.size > 1 and not np.all(np.diff(self.mz) == 1):
            raise ValueError("m/z axis must be contiguous unit masses")
        if self.times.size > 1 and np.any(np.diff(self.times) <= 0):
            raise ValueError("scan times must be strictly increasing")

    # ------------------------------------------------------------------ shape
    @property
    def n_scans(self) -> int:
        return self.times.size

    @property
    def n_mz(self) -> int:
        return self.mz.size

    @property
    def nnz(self) -> int:
        return int(self.data.nnz)

    @property
    def duration(self) -> float:
        return float(self.times[-1] - self.times[0]) if self.n_scans else 0.0

    # ---------------------------------------------------------------- signals
    @cached_property
    def tic(self) -> np.ndarray:
        """Total ion signal per scan (float64)."""
        return np.asarray(self.data.sum(axis=1), dtype=np.float64).ravel()

    @cached_property
    def _csc(self) -> sp.csc_matrix:
        # Column slicing is much faster on CSC; built once, on first use.
        return self.data.tocsc()

    def mz_index(self, mz: int | Iterable[int]) -> np.ndarray:
        """Column indices for one or more unit m/z values."""
        m = np.atleast_1d(np.asarray(mz, dtype=np.int64))
        idx = m - self.mz[0]
        bad = (idx < 0) | (idx >= self.n_mz)
        if np.any(bad):
            raise KeyError(f"m/z {m[bad].tolist()} outside {self.mz[0]}-{self.mz[-1]}")
        return idx

    def xic(self, mz: int | Iterable[int]) -> np.ndarray:
        """Extracted ion signal per scan, summed over the given m/z values."""
        cols = self.mz_index(mz)
        sub = self._csc[:, cols]
        return np.asarray(sub.sum(axis=1), dtype=np.float64).ravel()

    def spectrum(self, scans: int | slice | np.ndarray) -> np.ndarray:
        """Dense spectrum (length n_mz), summed over the requested scans."""
        rows = self.data[scans]
        if rows.ndim == 1 or rows.shape[0] == 1:
            return np.asarray(rows.toarray(), dtype=np.float64).ravel()
        return np.asarray(rows.sum(axis=0), dtype=np.float64).ravel()

    def scan_at(self, t: float) -> int:
        """Index of the scan nearest to time t (seconds)."""
        i = int(np.searchsorted(self.times, t))
        if i <= 0:
            return 0
        if i >= self.n_scans:
            return self.n_scans - 1
        return i if (self.times[i] - t) < (t - self.times[i - 1]) else i - 1

    def with_times(self, times: np.ndarray) -> "Run":
        """Copy of this run with new scan times (data shared)."""
        return Run(times=np.asarray(times, float), mz=self.mz, data=self.data,
                   name=self.name, source=self.source, meta=dict(self.meta))

    def __repr__(self) -> str:
        return (f"Run(name={self.name!r}, scans={self.n_scans}, "
                f"t={self.times[0]:.2f}-{self.times[-1]:.2f} s, "
                f"m/z={self.mz[0]}-{self.mz[-1]}, nnz={self.nnz:,})")
