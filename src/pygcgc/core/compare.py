"""⚖️ Comparing two runs: images, pixel spectra and peaks.

Run A is the one everything is drawn on; run B is read onto A's grid.
Three questions, three tools:

- Images. Where does the signal differ? B's ion image is resampled onto
  A's pixels, then shown in turn (blink), side by side, subtracted or as a
  two-colour overlay.
- Pixel spectra. Where does the chemistry differ, whatever the amounts?
  At every pixel the full mass spectrum of A is scored against B's at the
  same place. A pixel where both runs hold the same compound scores near 1
  even if one run holds ten times more of it; a pixel where a different
  compound elutes scores low.
- Peaks. Which detected peaks are in both runs, which in only one, and by
  how much do the shared ones differ? Peaks are paired one to one: first
  by identity where both runs carry the same confident name (a 🏷️
  internal standard or a tier 2a/2b library identification), then by
  position and spectrum together.

Coordinates. When both runs have a 📍 transform (or one is the reference
and the other has a transform), B is read through them. Otherwise B is
read at the same (t1, t2) as A, which is fine for runs acquired close
together whose retention drifts by less than a peak width, and is what the
position windows of the peak matching allow for.

Spectral similarity is the weighted cosine that 🏷️ library search uses
(intensity to the power 0.6 times m/z cubed, after Stein and Scott, 1994),
squared, so it runs from 0 (nothing shared) to 1 (identical) and means the
same as the forward term of a library match score.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy import ndimage as ndi

from .align import blank_mapping
from .fold import Fold
from .run import Run

WEIGHT_A = 0.6          # intensity exponent, as IdentifySettings.weight_a
WEIGHT_B = 3.0          # m/z exponent, as IdentifySettings.weight_b
MAX_DT1 = 12.0          # s: 1.5 modulations at 8 s, and 📍's own anchor search window
MAX_DT2 = 0.2           # s: about one second-dimension peak width
MIN_SIMILARITY = 0.7    # the score 🏷️ asks of a confident library match


# ------------------------------------------------------------------ coordinates
def pair_mapping(transform_a, transform_b) -> Callable:
    """A function taking A's (t1, t2) to B's, through the shared 📍
    reference. `None` for a transform means that run is the reference or is
    compared without alignment; `None` for both is the identity. The same
    composition 🧹 uses to read a blank through the alignment."""
    return blank_mapping(transform_a, transform_b)


def source_pixels(src_fold: Fold, dst_fold: Fold, to_src: Callable) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Where every pixel of `dst_fold`'s grid falls in `src_fold`'s image, as
    fractional (row, col) arrays of the destination's shape, plus a mask of
    destination pixels outside the source run. Wrapping in t2 is followed
    into the next modulation, as `align.resample_between` does."""
    t1_grid, t2_grid = np.meshgrid(dst_fold.t1, np.arange(dst_fold.n_rows) * dst_fold.row_pitch)
    t1_src, t2_src = to_src(t1_grid.ravel(), t2_grid.ravel())
    col = (np.asarray(t1_src, dtype=float) - src_fold.t1[0]) / src_fold.spec.period
    row = np.asarray(t2_src, dtype=float) / src_fold.row_pitch
    k = np.floor(row / src_fold.n_rows)
    row = row - k * src_fold.n_rows
    col = col + k
    outside = (col < 0) | (col > src_fold.n_cols - 1)
    shape = (dst_fold.n_rows, dst_fold.n_cols)
    return row.reshape(shape), col.reshape(shape), outside.reshape(shape)


def resample(image: np.ndarray, pixels: tuple[np.ndarray, np.ndarray, np.ndarray]) -> np.ndarray:
    """`image` read at the fractional pixels from `source_pixels`, by bilinear
    interpolation; zero outside the source run."""
    row, col, outside = pixels
    out = ndi.map_coordinates(np.nan_to_num(image, nan=0.0), [row.ravel(), col.ravel()], order=1,
                              mode="constant", cval=0.0).reshape(row.shape)
    return np.where(outside, 0.0, out)


# ------------------------------------------------------------------ spectra
def similarity(mz: np.ndarray, u: np.ndarray, v: np.ndarray,
               weight_a: float = WEIGHT_A, weight_b: float = WEIGHT_B) -> float:
    """Weighted cosine, squared, of two spectra on the same m/z axis."""
    w = np.power(np.clip(np.asarray(mz, dtype=float), 1, None), weight_b)
    wu = np.power(np.clip(u, 0, None), weight_a) * w
    wv = np.power(np.clip(v, 0, None), weight_a) * w
    denom = float(np.sqrt(float((wu ** 2).sum()) * float((wv ** 2).sum())))
    return 0.0 if denom <= 0 else float((wu * wv).sum() / denom) ** 2


def common_axis(mz_a: np.ndarray, mz_b: np.ndarray) -> tuple[np.ndarray, slice, slice]:
    """The m/z values both runs recorded, and the slices of each run's axis that give them."""
    lo, hi = max(int(mz_a[0]), int(mz_b[0])), min(int(mz_a[-1]), int(mz_b[-1]))
    if hi < lo:
        return np.array([], dtype=np.int64), slice(0, 0), slice(0, 0)
    return (np.arange(lo, hi + 1), slice(lo - int(mz_a[0]), hi - int(mz_a[0]) + 1),
            slice(lo - int(mz_b[0]), hi - int(mz_b[0]) + 1))


def pixel_spectrum(run: Run, fold: Fold, t1: float, t2: float, half_width: int = 1) -> np.ndarray:
    """The spectrum at one point of the 2D image: the scan nearest to it and
    `half_width` scans either side, summed. `t1` may fall between
    modulation starts (a mapped coordinate); the nearest modulation is used."""
    col = int(round((t1 - fold.t1[0]) / fold.spec.period))
    scan = run.scan_at(float(fold.t1[0] + col * fold.spec.period + t2))
    return run.spectrum(slice(max(scan - half_width, 0), min(scan + half_width + 1, run.n_scans)))


@dataclass
class SimilarityMap:
    similarity: np.ndarray      # on A's grid; NaN where neither run has any signal
    tic_a: np.ndarray           # A's summed signal over the shared m/z range
    tic_b: np.ndarray           # B's, resampled onto A's grid
    mz: np.ndarray              # the shared m/z range that was compared


def similarity_map(run_a: Run, fold_a: Fold, run_b: Run, fold_b: Fold, to_b: Callable,
                   weight_a: float = WEIGHT_A, weight_b: float = WEIGHT_B,
                   progress=None, cancel=None) -> SimilarityMap:
    """Spectral similarity of A and B at every pixel of A's grid. Built one
    ion at a time (fold, resample, accumulate the three sums of the cosine),
    so memory stays at a few images however many ions there are."""
    mz, cols_a, cols_b = common_axis(run_a.mz, run_b.mz)
    pixels = source_pixels(fold_b, fold_a, to_b)
    dot = np.zeros(fold_a.shape)
    norm_a = np.zeros(fold_a.shape)
    norm_b = np.zeros(fold_a.shape)
    tic_a = np.zeros(fold_a.shape)
    tic_b = np.zeros(fold_a.shape)
    csc_a = run_a.data.tocsc()[:, cols_a]
    csc_b = run_b.data.tocsc()[:, cols_b]
    for k, m in enumerate(mz):
        if cancel is not None:
            cancel()
        if progress is not None and k % 20 == 0:
            progress(k, mz.size, f"m/z {m}")
        xa = np.asarray(csc_a[:, k].toarray(), dtype=np.float64).ravel()
        xb = np.asarray(csc_b[:, k].toarray(), dtype=np.float64).ravel()
        if not xa.any() and not xb.any():
            continue
        img_a = fold_a.image(xa)
        img_b = resample(fold_b.image(xb), pixels)
        tic_a += img_a
        tic_b += img_b
        w = float(m) ** weight_b
        wa = np.power(np.clip(img_a, 0, None), weight_a) * w
        wb = np.power(np.clip(img_b, 0, None), weight_a) * w
        dot += wa * wb
        norm_a += wa ** 2
        norm_b += wb ** 2
    denom = np.sqrt(norm_a * norm_b)
    with np.errstate(invalid="ignore", divide="ignore"):
        sim = np.where(denom > 0, (dot / denom) ** 2, np.nan)
    sim[(norm_a > 0) != (norm_b > 0)] = 0.0       # signal in one run only: nothing shared
    return SimilarityMap(similarity=sim, tic_a=tic_a, tic_b=tic_b, mz=mz)


# ------------------------------------------------------------------ peaks
def peak_spectra(peaks) -> np.ndarray:
    """Each peak's own spectrum: the apex (single-scan) one where it was
    computed, the footprint sum otherwise, as 🏷️ identification does."""
    apex = getattr(peaks, "apex_spectra", np.zeros((0, 0)))
    return apex if apex.shape[0] == len(peaks) and apex.size else peaks.spectra


def attach_labels(t1: np.ndarray, t2: np.ndarray, labels: list[tuple[float, float, str]],
                  max_dt1: float, max_dt2: float) -> list[str]:
    """Hang each (t1, t2, text) label on the nearest peak within the window
    (labels were found on this run's own peaks, so they sit on one; the
    window only absorbs a peak list rebuilt since). Returns one string per
    peak, several labels on one peak joined by "; "."""
    out: list[list[str]] = [[] for _ in range(len(t1))]
    for lt1, lt2, text in labels:
        if len(t1) == 0 or not np.isfinite(lt1) or not np.isfinite(lt2):
            continue
        d1, d2 = np.abs(t1 - lt1), np.abs(t2 - lt2)
        ok = (d1 <= max_dt1) & (d2 <= max_dt2)
        if ok.any():
            i = int(np.argmin(np.where(ok, np.hypot(d1 / max_dt1, d2 / max_dt2), np.inf)))
            if text not in out[i]:
                out[i].append(text)
    return ["; ".join(x) for x in out]


@dataclass
class PeakPair:
    """One row of the peak comparison: a peak in both runs, or in one only
    (`a` or `b` is None). Positions are in each run's own coordinates;
    `dt1` and `dt2` are B's position, mapped onto A, minus A's."""
    a: int | None
    b: int | None
    t1_a: float = np.nan
    t2_a: float = np.nan
    t1_b: float = np.nan
    t2_b: float = np.nan
    dt1: float = np.nan
    dt2: float = np.nan
    similarity: float = np.nan
    volume_a: float = 0.0
    volume_b: float = 0.0
    name_a: str = ""
    name_b: str = ""
    class_a: str = ""
    class_b: str = ""
    how: str = ""               # "name", "position and spectrum", or "" when unpaired

    @property
    def status(self) -> str:
        return "both" if self.a is not None and self.b is not None else ("only A" if self.a is not None else "only B")

    @property
    def log2_ratio(self) -> float:
        """log2 of B's volume over A's; ±inf for a peak in one run only."""
        if self.volume_a <= 0 and self.volume_b <= 0:
            return np.nan
        if self.volume_a <= 0:
            return np.inf
        if self.volume_b <= 0:
            return -np.inf
        return float(np.log2(self.volume_b / self.volume_a))

    @property
    def name(self) -> str:
        return self.name_a or self.name_b


def match_peaks(a, b, to_a: Callable, names_a: list[str] | None = None, names_b: list[str] | None = None,
                max_dt1: float = MAX_DT1, max_dt2: float = MAX_DT2, min_similarity: float = MIN_SIMILARITY,
                weight_a: float = WEIGHT_A, weight_b: float = WEIGHT_B) -> list[PeakPair]:
    """Pair the peaks of runs A and B one to one.

    `a` and `b` are peak lists shaped like `BlobSet` (t1, t2, volume, mz,
    spectra). `to_a` maps B's (t1, t2) onto A's. `names_a` and `names_b` are
    confident identities per peak ("" for none); several names on one peak
    are "; "-separated.

    First pass: a name carried by one peak in each run pairs them, whatever
    their distance, because an identity already made the match. Second pass:
    every remaining A peak and B peak within the windows (after mapping) is
    scored on spectral similarity, and pairs are accepted best first, each
    peak used once, down to `min_similarity`. Unpaired peaks become "only A"
    or "only B" rows. Rows come back in A's t1 order (B's for B-only rows)."""
    n_a, n_b = len(a), len(b)
    names_a = names_a or [""] * n_a
    names_b = names_b or [""] * n_b
    mz, cols_a, cols_b = common_axis(np.asarray(a.mz), np.asarray(b.mz))
    spec_a = np.asarray(peak_spectra(a), dtype=float)[:, cols_a] if n_a else np.zeros((0, mz.size))
    spec_b = np.asarray(peak_spectra(b), dtype=float)[:, cols_b] if n_b else np.zeros((0, mz.size))
    t1_b_on_a, t2_b_on_a = (np.asarray(v, dtype=float) for v in to_a(np.asarray(b.t1, dtype=float),
                                                                       np.asarray(b.t2, dtype=float)))
    paired: dict[int, tuple[int, str]] = {}
    used_b: set[int] = set()

    def split(s: str) -> set[str]:
        return {x.strip() for x in s.split(";") if x.strip()}

    # Pass 1: identities.
    where_b: dict[str, list[int]] = {}
    for j, s in enumerate(names_b):
        for nm in split(s):
            where_b.setdefault(nm, []).append(j)
    for i, s in enumerate(names_a):
        for nm in sorted(split(s)):
            js = [j for j in where_b.get(nm, []) if j not in used_b]
            if i in paired or not js:
                continue
            j = min(js, key=lambda j: np.hypot((t1_b_on_a[j] - a.t1[i]) / max_dt1, (t2_b_on_a[j] - a.t2[i]) / max_dt2))
            paired[i] = (j, "name")
            used_b.add(j)

    # Pass 2: position and spectrum.
    candidates = []
    for i in range(n_a):
        if i in paired:
            continue
        near = np.flatnonzero((np.abs(t1_b_on_a - a.t1[i]) <= max_dt1) & (np.abs(t2_b_on_a - a.t2[i]) <= max_dt2))
        for j in near:
            if int(j) in used_b:
                continue
            s = similarity(mz, spec_a[i], spec_b[j], weight_a, weight_b)
            if s >= min_similarity:
                d = np.hypot((t1_b_on_a[j] - a.t1[i]) / max_dt1, (t2_b_on_a[j] - a.t2[i]) / max_dt2)
                candidates.append((-s, d, i, int(j)))
    for _neg, _d, i, j in sorted(candidates):
        if i not in paired and j not in used_b:
            paired[i] = (j, "position and spectrum")
            used_b.add(j)

    def row(i: int | None, j: int | None, how: str = "") -> PeakPair:
        p = PeakPair(a=i, b=j, how=how)
        if i is not None:
            p.t1_a, p.t2_a, p.volume_a, p.name_a = float(a.t1[i]), float(a.t2[i]), float(a.volume[i]), names_a[i]
        if j is not None:
            p.t1_b, p.t2_b, p.volume_b, p.name_b = float(b.t1[j]), float(b.t2[j]), float(b.volume[j]), names_b[j]
        if i is not None and j is not None:
            p.dt1, p.dt2 = float(t1_b_on_a[j] - a.t1[i]), float(t2_b_on_a[j] - a.t2[i])
            p.similarity = similarity(mz, spec_a[i], spec_b[j], weight_a, weight_b)
        return p

    rows = [row(i, paired[i][0], paired[i][1]) if i in paired else row(i, None) for i in range(n_a)]
    rows += [row(None, j) for j in range(n_b) if j not in used_b]
    return sorted(rows, key=lambda p: p.t1_a if p.a is not None else float(t1_b_on_a[p.b]))


def summarise_pairs(pairs: list[PeakPair]) -> dict:
    """Counts and the median volume ratio of the shared peaks, for a one-line summary."""
    both = [p for p in pairs if p.status == "both"]
    ratios = np.array([p.log2_ratio for p in both if np.isfinite(p.log2_ratio)])
    return {"both": len(both), "only A": sum(p.status == "only A" for p in pairs),
            "only B": sum(p.status == "only B" for p in pairs),
            "named": sum(bool(p.name) for p in pairs),
            "median_log2_ratio": float(np.median(ratios)) if ratios.size else np.nan}


PAIR_FIELDS = ["status", "matched_by", "name_a", "name_b", "class_a", "class_b", "t1_a_min", "t2_a_s",
               "t1_b_min", "t2_b_s", "dt1_s", "dt2_s", "spectral_similarity", "volume_a", "volume_b",
               "log2_b_over_a"]


def write_pairs_csv(pairs: list[PeakPair], path, run_a: str = "", run_b: str = "") -> Path:
    """One row per peak pair or unpaired peak, with the two run names in a header comment."""
    path = Path(path)

    def f(x, fmt):
        return "" if x is None or not np.isfinite(x) else format(x, fmt)

    with open(path, "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# A = {run_a}; B = {run_b}\n")
        w = csv.writer(fh)
        w.writerow(PAIR_FIELDS)
        for p in pairs:
            r = p.log2_ratio
            w.writerow([p.status, p.how, p.name_a, p.name_b, p.class_a, p.class_b,
                        f(p.t1_a / 60.0, ".4f"), f(p.t2_a, ".3f"), f(p.t1_b / 60.0, ".4f"), f(p.t2_b, ".3f"),
                        f(p.dt1, ".2f"), f(p.dt2, ".3f"), f(p.similarity, ".3f"),
                        f(p.volume_a, ".6g"), f(p.volume_b, ".6g"),
                        "inf" if r == np.inf else "-inf" if r == -np.inf else f(r, ".3f")])
    return path


def class_shares(summary_a: list, summary_b: list) -> list[tuple]:
    """Join two `ClassResult.summary` lists into (class, n A, n B, % A, % B,
    B minus A in percentage points), largest combined share first. Shares,
    not volumes, so runs of different total amount compare directly."""
    a = {name: (n, pct) for name, n, _vol, pct in summary_a}
    b = {name: (n, pct) for name, n, _vol, pct in summary_b}
    rows = []
    for name in set(a) | set(b):
        na, pa = a.get(name, (0, 0.0))
        nb, pb = b.get(name, (0, 0.0))
        rows.append((name, na, nb, pa, pb, pb - pa))
    return sorted(rows, key=lambda r: (-(r[3] + r[4]), r[0]))
