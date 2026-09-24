"""Known-answer peaks for benchmarking the automatic component finder.

Semi-synthetic means: real background, real spectra, invented peak positions.
Peaks are written into a run at scan level (the way the instrument would sample
them), so the complete automatic pipeline sees them: candidate detection, local
windows, model selection and spectrum fitting. Nothing is told where a peak is.

Truth is given in fold-image coordinates (``row`` in second-dimension rows,
``col`` in modulations), which is what fitted components report. ``area`` is the
pixel-summed signal over all ions before any window or ion selection.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy import sparse
from scipy.optimize import linear_sum_assignment

from .fold import Fold
from .overlap import _emg_1d, _gaussian_1d
from .run import Run


@dataclass(frozen=True)
class TruePeak:
    row: float
    col: float
    area: float
    spectrum: np.ndarray          # non-negative, one value per run.mz channel
    sigma_row: float = 3.0
    sigma_col: float = 0.8
    tau_row: float = 0.0          # exponential tail along t2 (rows); 0 = Gaussian
    tau_col: float = 0.0          # exponential tail along t1 (modulations)

    def profile(self, row, col) -> np.ndarray:
        gr = (_emg_1d(row, self.row, self.sigma_row, self.tau_row) if self.tau_row > 0
              else _gaussian_1d(row, self.row, self.sigma_row))
        gc = (_emg_1d(col, self.col, self.sigma_col, self.tau_col) if self.tau_col > 0
              else _gaussian_1d(col, self.col, self.sigma_col))
        return gr * gc

    def apex(self) -> tuple[float, float]:
        """Position of the profile maximum; differs from the nominal centre when tailing."""
        rows = np.arange(self.row - 6 * self.sigma_row - 8 * self.tau_row,
                         self.row + 6 * self.sigma_row + 8 * self.tau_row, 0.1)
        cols = np.arange(self.col - 4 * self.sigma_col - 6 * self.tau_col,
                         self.col + 4 * self.sigma_col + 6 * self.tau_col, 0.05)
        r = rows[np.argmax(self.profile(rows, np.full_like(rows, self.col)))]
        c = cols[np.argmax(self.profile(np.full_like(cols, self.row), cols))]
        return float(r), float(c)


@dataclass
class TruthCase:
    """One planted scenario: what was injected and how it should be scored."""
    name: str
    kind: str                                  # single, single_broad, single_tail, pair, pair_tail
    peaks: list
    tags: dict = field(default_factory=dict)   # ratio, separation, cosine, ... for grouping
    tic_maxima: int = 0                        # local maxima of the true summed profile


def inject_peaks(run: Run, fold: Fold, peaks: list[TruePeak], seed: int = 0,
                 poisson: bool = True) -> Run:
    """Return ``run`` with the peaks added at scan level.

    A scan's row is its fractional position inside its modulation and its column
    is its modulation index, exactly as :class:`Fold` places it. Counting noise
    with variance equal to the added signal is applied to the added part only.
    """
    rng = np.random.default_rng(seed)
    n_rows = fold.n_rows
    lo, weight = fold._lo, fold._w
    valid = lo >= 0
    scan_col = np.where(valid, lo // n_rows, -1).astype(float)
    scan_row = np.where(valid, lo % n_rows + weight, -1.0)
    blocks = []
    for peak in peaks:
        shape = np.where(valid, peak.profile(scan_row, scan_col), 0.0)
        keep = np.flatnonzero(shape > 1e-4)
        if not keep.size:
            continue
        spectrum = np.clip(np.asarray(peak.spectrum, dtype=float), 0.0, None)
        # area is the pixel sum: sum(shape over scans) x sum(spectrum) x scale
        scale = peak.area / (shape[keep].sum() * spectrum.sum())
        values = scale * shape[keep][:, None] * spectrum[None, :]
        if poisson:
            values = values + rng.normal(0.0, np.sqrt(np.clip(values, 0.0, None)))
        blocks.append(sparse.csr_matrix((values.ravel(), (np.repeat(keep, spectrum.size),
                                                  np.tile(np.arange(spectrum.size), keep.size))),
                                        shape=run.data.shape))
    data = run.data.copy()
    for block in blocks:
        data = data + block
    return replace(run, data=sparse.csr_matrix(data))


def real_spectrum_library(clean: Run, fold: Fold, candidates, min_apex_tic: float = 8e4,
                          min_ions: int = 6, isolation_rows: float = 8.0,
                          isolation_cols: float = 1.0) -> np.ndarray:
    """Spectra read from the apexes of strong real peaks of the run (rows: spectra, sum 1).

    The local median is removed, so background chemistry does not become part of
    a "compound". A crowded chromatogram has almost no candidate free of
    neighbours, so only apexes with no other candidate within a few pixels are used;
    the spectra are real but not guaranteed pure, which is adequate as planted truth.
    """
    cube = fold.cube(clean, dtype=np.float32)
    tic = cube.sum(axis=2)
    spectra = []
    for row, col in zip(candidates.row, candidates.col):
        if row < 30 or row >= fold.n_rows - 30 or col < 5 or col >= fold.n_cols - 5:
            continue
        near = (np.abs(candidates.row - row) <= isolation_rows) & (
            np.abs(candidates.col - col) <= isolation_cols)
        if near.sum() > 1 or tic[row, col] < min_apex_tic:
            continue
        apex = cube[row - 1:row + 2, col].mean(axis=0)
        ring = cube[row - 30:row + 31, col - 4:col + 5].reshape(-1, cube.shape[2])
        spectrum = np.clip(apex - np.median(ring, axis=0), 0.0, None)
        if spectrum.sum() <= 0 or np.count_nonzero(spectrum >= 0.03 * spectrum.sum()) < min_ions:
            continue
        spectra.append(spectrum / spectrum.sum())
    return np.asarray(spectra)


def cosine(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.dot(a, b) / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-12))


def spectral_pairs(library: np.ndarray, lo: float, hi: float, rng, n: int) -> list:
    """Random library pairs whose cosine similarity lies in [lo, hi)."""
    found, tries = [], 0
    while len(found) < n and tries < 20000:
        i, j = rng.choice(len(library), size=2, replace=False)
        if lo <= cosine(library[i], library[j]) < hi:
            found.append((int(i), int(j)))
        tries += 1
    return found


def free_sites(candidates, fold: Fold, rng, n: int, clearance_rows: float = 25.0,
               clearance_cols: float = 3.0, spacing_rows: float = 30.0,
               spacing_cols: float = 6.0, row_range=(45, 330), edge_cols: int = 9) -> list:
    """Random (row, col) positions away from every real candidate and from each other."""
    sites = []
    for _ in range(20000):
        if len(sites) == n:
            break
        row = float(rng.uniform(*row_range))
        col = float(rng.uniform(edge_cols, fold.n_cols - edge_cols))
        if np.any((np.abs(candidates.row - row) <= clearance_rows) &
                  (np.abs(candidates.col - col) <= clearance_cols)):
            continue
        if any(abs(row - r) <= spacing_rows and abs(col - c) <= spacing_cols for r, c in sites):
            continue
        sites.append((row, col))
    return sites


def count_tic_maxima(peaks: list[TruePeak]) -> int:
    """How many separate maxima the noise-free summed profile has (1 means hidden overlap)."""
    r0 = min(p.row for p in peaks) - 25
    c0 = min(p.col for p in peaks) - 5
    rows, cols = np.mgrid[r0:r0 + 70, c0:c0 + 14]
    tic = sum(p.area * p.profile(rows.ravel().astype(float), cols.ravel().astype(float))
              .reshape(rows.shape) / (2 * np.pi * p.sigma_row * p.sigma_col) for p in peaks)
    from scipy import ndimage as ndi
    peaks_found = (tic == ndi.maximum_filter(tic, size=(3, 3), mode="constant")) & (tic > 0.05 * tic.max())
    return int(ndi.label(peaks_found, structure=np.ones((3, 3)))[1])


# ----------------------------------------------------------------------- scoring
def _components_global(component_set) -> list:
    """(row, col, component) with fold-image coordinates instead of window-local ones."""
    out = []
    for group in component_set.groups:
        for component in group.components:
            out.append((group.row0 + component.row0, group.col0 + component.col0, component))
    return out


def score_case(case: TruthCase, component_set, library_ions: np.ndarray,
               ignore: list | None = None, radius_rows: float = 15.0,
               radius_cols: float = 3.0, gate_rows: float = 8.0, gate_cols: float = 2.0,
               pitch: float = 0.02, period: float = 10.0, faint_ratio: float = 0.05) -> dict:
    """Match fitted components near the truth to the planted peaks and measure errors.

    ``ignore`` lists (row, col) positions of components that already existed in
    the un-injected run; they are not counted as extras. An unmatched component
    below ``faint_ratio`` of its group's largest is a faint extra, not a material one.
    """
    truths = case.peaks
    apexes = [p.apex() for p in truths]
    centre_row = float(np.mean([a[0] for a in apexes]))
    centre_col = float(np.mean([a[1] for a in apexes]))
    nearby = []
    for row, col, component in _components_global(component_set):
        if abs(row - centre_row) > radius_rows or abs(col - centre_col) > radius_cols:
            continue
        if ignore and any(abs(row - r) < 1.5 and abs(col - c) < 0.75 for r, c in ignore):
            continue
        nearby.append((row, col, component))
    n_true, n_found = len(truths), len(nearby)
    cost = np.full((n_true, max(n_found, 1)), 1e6)
    for i, (tr, tc) in enumerate(apexes):
        for j, (row, col, _c) in enumerate(nearby):
            if abs(row - tr) <= gate_rows and abs(col - tc) <= gate_cols:
                cost[i, j] = np.hypot((row - tr) / 3.0, (col - tc) / 1.0)
    matches = []
    if n_found:
        for i, j in zip(*linear_sum_assignment(cost)):
            if cost[i, j] < 1e5:
                matches.append((int(i), int(j)))
    peak_rows = []
    for i, j in matches:
        row, col, component = nearby[j]
        truth = truths[i]
        ions = np.asarray(component.mz, dtype=int) - int(library_ions[0])
        truth_spectrum = np.clip(truth.spectrum[ions], 0.0, None)
        truth_area = truth.area * truth_spectrum.sum() / max(np.clip(truth.spectrum, 0, None).sum(), 1e-12)
        peak_rows.append({
            "truth_index": i, "ratio_to_largest": truth.area / max(p.area for p in truths),
            "row_error_s": (row - apexes[i][0]) * pitch,
            "col_error_s": (col - apexes[i][1]) * period,
            "cosine": cosine(truth_spectrum, np.clip(component.spectrum, 0.0, None)),
            "area_error": (component.volume - truth_area) / truth_area,
            "confidence": component.confidence})
    matched_truth = {i for i, _ in matches}
    matched_found = {j for _, j in matches}
    faint = sum(1 for j, (_r, _c, comp) in enumerate(nearby)
                if j not in matched_found and comp.relative_to_largest < faint_ratio)
    return {"name": case.name, "kind": case.kind, "tags": case.tags,
            "tic_maxima": case.tic_maxima, "n_true": n_true, "n_found": n_found,
            "n_matched": len(matches), "count_correct": n_found == n_true and len(matches) == n_true,
            "missed": [i for i in range(n_true) if i not in matched_truth],
            "extra": n_found - len(matches), "extra_faint": faint,
            "extra_material": n_found - len(matches) - faint, "peaks": peak_rows}
