"""Background removal, despiking and the noise map.

The idea in plain terms
-----------------------
Background changes slowly along the first dimension (column bleed rises
with oven temperature), and sometimes carries a faint pattern that repeats
every modulation (modulator pulses, second-column bleed). Peaks are sharp
and sparse. So, ion by ion:

1. In every modulation, take a low value down the t2 axis (the 20th
   percentile by default). Peaks occupy only a few rows, so this lands on
   the background. That gives a rough baseline along t1.
2. Smooth it with arPLS, an asymmetric smoother that hugs the bottom of a
   signal and ignores anything poking up (Baek et al., 2015).
3. Optionally remove the repeating t2 pattern: for each t2 row, take a low
   rolling percentile across many neighbouring modulations. A stripe that
   runs through the whole run survives; a peak spanning a few modulations
   does not.
4. Subtract the background from every stored data point. Small negative
   values are kept, so noise stays unbiased; only the display clips them.

Before that, single-scan spikes are removed: a point far above all four of
its neighbouring scans cannot be a real peak, because real peaks span
several scans.

Optionally, an already-cleaned blank is then subtracted, scaled per ion,
and never taking a point below zero. The scale is worked out separately
for each of the ion's own artefacts (its significant pixels, grouped into
connected regions), so two artefacts sharing an ion but varying
independently between tubes no longer share one compromise scale factor;
whatever significant signal isn't part of any of those regions still
falls back to one scale for the ion as a whole, exactly as before.

Finally the noise map records, for every ion and every stretch of t1, the
spread of what's left away from peaks. Every signal-to-noise ratio later
on comes from it. The instrument doesn't store values below its threshold,
so where most pixels are empty the noise can't be measured; there it is
set to a third of the storage threshold.

Only ions that actually carry background are processed in steps 1 to 3
(an ion whose pixels are mostly empty has a background of zero), which
keeps a whole run to well under a second.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace

import numpy as np
import scipy.sparse as sp
from scipy.linalg import solveh_banded
from scipy.ndimage import label as _connected_components
from scipy.ndimage import percentile_filter
from scipy.stats import norm

from .align import resample_between
from .fold import Fold
from .run import Run


REFERENCE_COLUMNS = 465   # a 31-minute run at a 4 s period

# ----------------------------------------------------------------- settings
@dataclass(frozen=True)
class BackgroundSettings:
    percentile: float = 20.0          # low value taken down each modulation
    stiffness: float = 1e4            # arPLS lambda: larger = stiffer baseline (1e4 bends over roughly 10 modulations)
    remove_pattern: bool = True       # remove repeating t2 structure
    pattern_window: int = 51          # modulations in the rolling window (odd)
    pattern_percentile: float = 25.0
    despike: bool = True
    spike_threshold: float = 8.0      # robust SDs above the neighbouring scans
    noise_block: int = 20             # modulations per noise-map block
    subtract_blank: bool = False      # the blank itself is chosen by the caller
    scale_blank: bool = True          # scale the blank per ion to each run

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "BackgroundSettings":
        names = cls.__dataclass_fields__
        return cls(**{k: v for k, v in (d or {}).items() if k in names})


# -------------------------------------------------------------------- arPLS
def _second_difference_bands(n: int) -> np.ndarray:
    """Upper bands of DᵀD for the second-difference matrix D (shape 3 × n)."""
    D = sp.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n))
    DtD = (D.T @ D).todia()
    ab = np.zeros((3, n))
    for k in (0, 1, 2):
        diag = DtD.diagonal(k)
        ab[2 - k, k:] = diag
    return ab


def arpls(y: np.ndarray, stiffness: float = 1e3, max_iter: int = 50, tol: float = 1e-3) -> np.ndarray:
    """Asymmetrically reweighted penalised least squares baseline.

    Points above the current baseline by more than the noise get little
    weight, so the curve follows the lower envelope of the signal.
    """
    y = np.asarray(y, dtype=np.float64)
    n = y.size
    if n < 3:
        return y.copy()
    H = _second_difference_bands(n) * stiffness
    w = np.ones(n)
    z = y
    for _ in range(max_iter):
        ab = H.copy()
        ab[2] += w
        z = solveh_banded(ab, w * y, lower=False, check_finite=False)
        d = y - z
        neg = d[d < 0]
        if neg.size < 2:
            break
        m, s = neg.mean(), neg.std()
        if s <= 0:
            break
        arg = np.clip(2.0 * (d - (2.0 * s - m)) / s, -50, 50)
        w_new = 1.0 / (1.0 + np.exp(arg))
        if np.linalg.norm(w - w_new) / max(np.linalg.norm(w), 1e-12) < tol:
            w = w_new
            break
        w = w_new
    return z


# --------------------------------------------------------------- helpers
def storage_threshold(run: Run) -> float:
    """The instrument's apparent storage threshold (low end of stored values)."""
    v = run.data.data
    v = v[v > 0]
    if v.size == 0:
        return 1.0
    return float(np.percentile(v, 0.1))


def _scan_positions(fold: Fold) -> tuple[np.ndarray, np.ndarray]:
    """Column and nearest row for every scan (scans before the first
    boundary are treated as belonging to column 0)."""
    col, t2 = fold.to_fold(fold.times)
    row = np.rint(t2 / fold.row_pitch).astype(np.int64)
    col = col + row // fold.n_rows
    row = row % fold.n_rows
    col = np.clip(col, 0, fold.n_cols - 1)
    return col, row


def _occupancy(run: Run, col_of_scan: np.ndarray, n_cols: int) -> np.ndarray:
    """Stored points per (ion, modulation)."""
    coo = run.data.tocoo()
    flat = coo.col.astype(np.int64) * n_cols + col_of_scan[coo.row]
    return np.bincount(flat, minlength=run.n_mz * n_cols).reshape(run.n_mz, n_cols)


def _fill_incomplete(b: np.ndarray, complete: np.ndarray) -> np.ndarray:
    if complete.all() or not complete.any():
        return b
    idx = np.arange(b.size)
    b = b.copy()
    b[~complete] = np.interp(idx[~complete], idx[complete], b[complete])
    return b


# ------------------------------------------------------------------ spikes
def despike(run: Run, threshold: float = 8.0, floor: float | None = None) -> tuple[Run, int]:
    """Remove single-scan spikes.

    A stored point is a spike when it stands more than `threshold` noise
    units above the largest of its four neighbouring scans (two either
    side, same ion) and is also more than four times that neighbour.
    Noise is modelled as counting noise plus a floor.
    """
    floor = storage_threshold(run) / 3.0 if floor is None else floor
    csc = run.data.tocsc(copy=True)
    csc.sort_indices()
    removed = 0
    for j in range(run.n_mz):
        a, b = csc.indptr[j], csc.indptr[j + 1]
        if b - a == 0:
            continue
        rows = csc.indices[a:b]
        vals = csc.data[a:b].astype(np.float64)
        neigh = np.zeros((4, rows.size))
        for k, off in enumerate((-2, -1, 1, 2)):
            pos = np.clip(np.searchsorted(rows, rows + off), 0, rows.size - 1)
            hit = rows[pos] == rows + off
            neigh[k] = np.where(hit, vals[pos], 0.0)
        nb = neigh.max(axis=0)
        sigma = np.sqrt(np.clip(nb, 0, None) + floor ** 2)
        spike = ((vals - nb) > threshold * sigma) & (vals > 4.0 * np.maximum(nb, floor))
        if spike.any():
            # Replace with what the neighbours say is there (usually nothing).
            csc.data[a:b][spike] = np.median(neigh[:, spike], axis=0)
            removed += int(spike.sum())
    csc.eliminate_zeros()
    out = Run(times=run.times, mz=run.mz, data=csc.tocsr(), name=run.name,
              source=run.source, meta=dict(run.meta))
    return out, removed


# ------------------------------------------------------------------ noise
@dataclass
class NoiseMap:
    """Noise (one standard deviation) per ion and per stretch of t1."""
    mz: np.ndarray
    edges: np.ndarray              # absolute times of block boundaries, length n_blocks + 1
    sigma: np.ndarray              # (n_mz, n_blocks)
    floor: float

    def block_of(self, t) -> np.ndarray:
        return np.clip(np.searchsorted(self.edges, t, side="right") - 1, 0, self.sigma.shape[1] - 1)

    def at(self, mz: int, t) -> np.ndarray:
        return self.sigma[int(mz) - int(self.mz[0]), self.block_of(t)]

    def for_ions(self, ions, t) -> np.ndarray:
        """Noise of a summed signal at times t (independent ions add in quadrature)."""
        blocks = self.block_of(t)
        if ions is None:
            s = self.sigma
        else:
            s = self.sigma[np.asarray(list(ions), dtype=np.int64) - int(self.mz[0])]
        return np.sqrt(np.sum(s[:, blocks] ** 2, axis=0))


def _robust_sigma(vals: np.ndarray, floor: float, iterations: int = 4) -> float:
    """Spread of the background pixels, ignoring peaks (iterative clipping)."""
    v = vals[np.isfinite(vals)]
    if v.size < 10 or np.mean(v != 0) < 0.5:
        return floor
    med = np.median(v)
    sig = 1.4826 * np.median(np.abs(v - med))
    for _ in range(iterations):
        if sig <= 0:
            break
        keep = v < med + 3.0 * sig
        med = np.median(v[keep])
        sig = 1.4826 * np.median(np.abs(v[keep] - med))
    return float(max(sig, floor))


def estimate_noise(run: Run, fold: Fold, block: int = 20, floor: float | None = None,
                   ions: np.ndarray | None = None) -> NoiseMap:
    floor = storage_threshold(run) / 3.0 if floor is None else floor
    n_blocks = max(1, int(np.ceil(fold.n_cols / block)))
    edges = fold.t1[0] + np.arange(n_blocks + 1) * block * fold.spec.period
    sigma = np.full((run.n_mz, n_blocks), floor)
    col_of_scan, _ = _scan_positions(fold)
    occ = _occupancy(run, col_of_scan, fold.n_cols)
    busy = np.flatnonzero((occ / fold.n_rows).max(axis=1) >= 0.5) if ions is None else ions
    complete = fold.complete
    for j in busy:
        img = fold.image(run.xic(int(run.mz[j])))
        for b in range(n_blocks):
            c0, c1 = b * block, min((b + 1) * block, fold.n_cols)
            cols = np.arange(c0, c1)[complete[c0:c1]]
            if cols.size:
                sigma[j, b] = _robust_sigma(img[:, cols].ravel(), floor)
    return NoiseMap(mz=run.mz.copy(), edges=edges, sigma=sigma, floor=floor)


# ------------------------------------------------------------------ result
@dataclass
class BackgroundResult:
    run: Run                                   # the corrected run
    noise: NoiseMap
    settings: BackgroundSettings
    removed_by_ion: np.ndarray                 # counts removed per ion (background + blank)
    total_before: float
    spikes_removed: int = 0
    ions_processed: list[int] = field(default_factory=list)
    baselines: dict[int, np.ndarray] = field(default_factory=dict)   # ion -> baseline per modulation
    blank_name: str | None = None
    blank_scale: dict[int, float] = field(default_factory=dict)

    @property
    def removed_fraction(self) -> float:
        return float(self.removed_by_ion.sum() / self.total_before) if self.total_before else 0.0

    def top_ions(self, n: int = 5) -> list[tuple[int, float]]:
        order = np.argsort(self.removed_by_ion)[::-1][:n]
        tot = max(self.removed_by_ion.sum(), 1e-12)
        return [(int(self.run.mz[i]), float(self.removed_by_ion[i] / tot))
                for i in order if self.removed_by_ion[i] > 0]

    def summary(self) -> str:
        parts = [f"Removed {self.removed_fraction:.1%} of the signal"]
        top = self.top_ions(4)
        if top:
            parts.append("mostly m/z " + ", ".join(str(m) for m, _ in top))
        if self.settings.despike:
            parts.append(f"{self.spikes_removed} spike{'s' if self.spikes_removed != 1 else ''}")
        if self.blank_name:
            ks = list(self.blank_scale.values())
            if ks:
                parts.append(f"blank {self.blank_name} scaled {min(ks):.2f} to {max(ks):.2f}")
            else:
                parts.append(f"blank {self.blank_name} (nothing to subtract)")
        return "; ".join(parts)

    def residual_pixels(self, fold: Fold, mz: int) -> tuple[np.ndarray, float]:
        """Background pixels left after correction (peaks excluded) and the
        noise-map SD they should match. Used by the diagnostics plots."""
        img = fold.image(self.run.xic(mz))[:, fold.complete]
        cols_t = fold.t1[fold.complete]
        sig = self.noise.at(mz, cols_t)
        keep = np.abs(img) < 3.0 * sig[None, :]
        return img[keep], float(np.median(sig))


# ------------------------------------------------------------------ main
def _background_surface(img: np.ndarray, complete: np.ndarray, s: BackgroundSettings):
    """Background level per modulation, plus the repeating t2 pattern (or None).

    A low percentile sits below the middle of the noise, so on its own it
    would leave the corrected background slightly positive. The noise spread
    is measured from the gap between the chosen percentile and one a quarter
    as high (both well below any peaks), and the percentile is then raised by
    the matching amount for normally distributed noise. The same correction
    is applied to the repeating pattern.
    """
    p = float(s.percentile)
    q_lo, q_p = np.percentile(img, [p / 4.0, p], axis=0)
    gap = norm.ppf(p / 100.0) - norm.ppf(p / 400.0)
    sigma = np.clip((q_p - q_lo) / gap, 0, None)
    centre = q_p - norm.ppf(p / 100.0) * sigma
    centre = _fill_incomplete(centre, complete)
    # The same stiffness should mean the same smoothness in minutes whatever
    # the run length, and the penalty grows with the fourth power of length.
    base = arpls(centre, s.stiffness * (centre.size / REFERENCE_COLUMNS) ** 4)
    pattern = None
    good = np.flatnonzero(complete)
    if s.remove_pattern and good.size >= 5:
        # Work on complete modulations only: a partly recorded one has empty
        # rows that would look like a deep dip.
        resid = img[:, good] - base[None, good]
        win = min(int(s.pattern_window), good.size)
        if win % 2 == 0:
            win -= 1
        win = max(3, win)
        pp = float(s.pattern_percentile)
        low = percentile_filter(resid, pp, size=(1, win), mode="nearest")
        sig = np.broadcast_to(sigma[good], resid.shape).copy()
        sig_s = percentile_filter(sig, 50, size=(1, win), mode="nearest")
        pat_good = low - norm.ppf(pp / 100.0) * sig_s
        # Only keep a pattern that stands clear of the noise.
        pat_good[np.abs(pat_good) < 2.0 * sig_s / np.sqrt(win)] = 0.0
        # Incomplete modulations borrow the pattern of the nearest complete one.
        nearest = good[np.clip(np.searchsorted(good, np.arange(img.shape[1])), 0, good.size - 1)]
        lookup = np.searchsorted(good, nearest)
        pattern = pat_good[:, lookup]
    return base, pattern


def remove_background(run: Run, fold: Fold, settings: BackgroundSettings = BackgroundSettings(),
                      blank: "BackgroundResult | None" = None, blank_fold: Fold | None = None,
                      blank_name: str | None = None, progress=None, to_blank=None) -> BackgroundResult:
    """Despike, remove background, optionally subtract a cleaned blank, and
    build the noise map. `progress(done, total, message)` is called as it goes.
    `to_blank(t1, t2)` maps this run's coordinates onto the blank's (see
    `align.blank_mapping`); left out, the two are assumed to line up pixel for pixel."""
    s = settings
    total_before = float(run.data.sum())
    floor = storage_threshold(run) / 3.0
    step = (lambda d, t, m: progress(d, t, m)) if progress else (lambda *a: None)

    spikes = 0
    work = run
    if s.despike:
        step(0, 4, "Removing spikes")
        work, spikes = despike(run, s.spike_threshold, floor)

    step(1, 4, "Estimating background")
    col_of_scan, row_of_scan = _scan_positions(fold)
    occ = _occupancy(work, col_of_scan, fold.n_cols)
    # An ion has background where most of a modulation's rows hold data.
    need = (100.0 - s.percentile) / 100.0
    candidates = np.flatnonzero((occ / fold.n_rows >= need).any(axis=1))

    csc = work.data.tocsc(copy=True)
    removed = np.zeros(run.n_mz)
    removed_spikes = np.asarray(run.data.sum(axis=0)).ravel() - np.asarray(work.data.sum(axis=0)).ravel()
    baselines: dict[int, np.ndarray] = {}
    complete = fold.complete
    for k, j in enumerate(candidates):
        a, b_ = csc.indptr[j], csc.indptr[j + 1]
        rows = csc.indices[a:b_]
        xic = np.zeros(run.n_scans)
        xic[rows] = csc.data[a:b_]
        img = fold.image(xic)
        base, pattern = _background_surface(img, complete, s)
        bg = base[col_of_scan[rows]]
        if pattern is not None:
            bg = bg + pattern[row_of_scan[rows], col_of_scan[rows]]
        csc.data[a:b_] -= bg.astype(csc.data.dtype)
        removed[j] += float(bg.sum())
        baselines[int(run.mz[j])] = base

    blank_scale: dict[int, float] = {}
    if s.subtract_blank and blank is not None and blank_fold is not None:
        step(2, 4, "Subtracting blank")
        # The blank is read onto this run's own pixel grid. Without an
        # alignment mapping that is a whole-modulation shift of t1 (the old
        # behaviour); with one (`to_blank`, from the fitted alignment of both
        # runs to the reference) it follows the retention drift between them.
        if to_blank is None:
            period = fold.spec.period
            shift = int(np.rint((fold.t1[0] - blank_fold.t1[0]) / period))
            offset = blank_fold.t1[0] + shift * period - fold.t1[0]

            def to_blank(t1, t2, offset=offset):
                return np.asarray(t1) + offset, np.asarray(t2)
        bcsc = blank.run.data.tocsc()
        all_cols = col_of_scan
        for j in range(run.n_mz):
            m = int(run.mz[j])
            bj = m - int(blank.run.mz[0])
            if bj < 0 or bj >= blank.run.n_mz or bcsc.indptr[bj + 1] == bcsc.indptr[bj]:
                continue
            bimg_own = blank_fold.image(blank.run.xic(m))
            bsig = blank.noise.at(m, blank_fold.t1)
            bimg = resample_between(bimg_own, blank_fold, fold, to_blank)
            significant = resample_between((bimg_own > 10.0 * bsig[None, :]).astype(float),
                                           blank_fold, fold, to_blank) > 0.5
            if significant.sum() < 3:
                continue
            a, b_ = csc.indptr[j], csc.indptr[j + 1]
            rows = csc.indices[a:b_]
            cols = all_cols[rows]
            k_scan = np.ones(rows.size)
            k_global = 1.0
            if s.scale_blank:
                xic = np.zeros(run.n_scans)
                xic[rows] = csc.data[a:b_]
                simg = fold.image(xic)
                # Compare sample and blank on the blank's significant pixels.
                sv, bv = simg[significant], bimg[significant]
                good = bv > 0
                k_global = float(np.clip(np.median(sv[good] / bv[good]), 0.0, 10.0)) if good.sum() >= 10 else 1.0
                # This ion's artefacts, taken separately: each connected
                # region of significant blank pixels gets its own ratio, so
                # one artefact scaling badly doesn't drag another's with it.
                groups, n_groups = _connected_components(significant)
                if n_groups > 1:
                    k_by_group = np.full(n_groups + 1, k_global)
                    for g in range(1, n_groups + 1):
                        gmask = groups == g
                        sv_g, bv_g = simg[gmask], bimg[gmask]
                        good_g = bv_g > 0
                        if good_g.sum() >= 10:
                            k_by_group[g] = float(np.clip(np.median(sv_g[good_g] / bv_g[good_g]), 0.0, 10.0))
                    k_scan[:] = k_by_group[groups[row_of_scan[rows], cols]]
                else:
                    k_scan[:] = k_global
            sub = k_scan * bimg[row_of_scan[rows], cols]
            # Never take a point below zero: blank peaks rarely line up exactly,
            # and over-subtraction would dig holes that look like negative peaks.
            sub = np.clip(np.minimum(sub, csc.data[a:b_]), 0, None)
            csc.data[a:b_] -= sub.astype(csc.data.dtype)
            removed[j] += float(sub.sum())
            blank_scale[m] = k_global

    step(3, 4, "Building noise map")
    corrected = Run(times=run.times, mz=run.mz, data=csc.tocsr(), name=run.name,
                    source=run.source, meta=dict(run.meta, background=s.to_dict()))
    noise = estimate_noise(corrected, fold, s.noise_block, floor, ions=candidates)
    step(4, 4, "Done")
    return BackgroundResult(run=corrected, noise=noise, settings=s,
                            removed_by_ion=removed + np.clip(removed_spikes, 0, None),
                            total_before=total_before, spikes_removed=spikes,
                            ions_processed=[int(run.mz[j]) for j in candidates],
                            baselines=baselines, blank_name=blank_name if blank_scale or blank else None,
                            blank_scale=blank_scale)
