"""Finding and measuring peaks ("blobs").

The folded image and wrap-around
--------------------------------
A compound whose second-dimension time is longer than the modulation
period arrives near the bottom of the *next* modulation. We work on an
image extended upwards by `wrap_below_s` (the bottom rows of the next
modulation repeated on top), and only look at rows from `wrap_below_s` up
to one period above it. Every recorded sample appears exactly once in
that band. Peaks below `wrap_below_s` in a modulation are therefore
treated as wrapped from the modulation before, and are reported with
their true t2 (longer than the period). Set it inside the quiet band just
above zero, where nothing is retained; 0.5 s suits the dummy data.

TIC blobs (`find_blobs`)
------------------------
1. Smooth the total-ion image a little (an elongated Gaussian, narrow
   along t2 and about one modulation along t1), and turn it into
   signal-to-noise by dividing by its noise measured in stretches of t1
   (robust spread of the pixels, ignoring peaks).
2. Do the same, ion by ion, for every ion with enough stored points to be
   worth the trouble (`min_ion_points`), and sum the per-ion S/N images.
   A compound spread thinly but consistently across many ions can stand
   out here even when it's unremarkable in the raw TIC, which is
   dominated by whichever handful of ions happen to be brightest. The
   detection surface used for seeding and the watershed is the pixel-wise
   maximum of this sum and the plain TIC S/N, so nothing the old,
   TIC-only method already found is lost.
3. Local maxima above `snr_min` (on that combined surface) are peak
   seeds.
4. A watershed from the seeds, inside the area above `floor_snr`, gives
   each peak its footprint: it grows downhill from each apex until it
   meets a neighbour or drops into the noise. Each footprint is then
   trimmed to `blob_rel_floor` of its own TIC apex (1% by default), so a
   big peak doesn't spread across low-level noise, and flooded again so
   every pixel stays connected to its apex.
5. Every stored data point is assigned to the footprint it falls in, so
   each blob gets a summed spectrum and a volume (the sum of everything
   inside it), plus the spectrum of the scan(s) nearest its apex pixel.
   Blobs spanning fewer than `min_modulations` are dropped.

On a long run, step 4 (seeding and the watershed) is split into
overlapping stretches of t1 (`chunk_modulations`, `chunk_overlap`) and run
across a thread pool, one peak processed exactly once by whichever chunk
owns its seed; correctness relies on the overlap being wider than any
peak actually seen (comfortably true of the dummy data's default). Short
runs (the common case so far) fall into a single chunk and behave exactly
as the unchunked code did.

Co-elution (`BlobSet.co_elution`). Every ion holding a decent share of a
blob's spectrum (`coelution_frac`) gets its own local apex pixel found
inside the footprint. If those apexes fall into more than one cluster
(single-linkage, `coelution_gap`, in units of one modulation along t1 and
two rows along t2) the blob is flagged as a probable co-elution: its
ions don't agree on where the compound actually peaks. This is a cheap
heuristic, not a separation; the overlap toolkit (milestone 8) resolves
what it finds.

Target peaks (`measure_ion_peak`)
---------------------------------
A target is measured on its own quantifier ion, which is far more
sensitive than the TIC and less disturbed by neighbours. Inside the
target's window, local maxima of the smoothed quantifier-ion image are
candidates. Each gets a footprint by the same watershed (seeded with every
maximum nearby, so neighbours keep their own signal), and the quantifier
and qualifier ions are summed over it. The best candidate is chosen in
`quant.py`.

Positions
---------
t1 is the intensity-weighted mean of the modulation mid-times across the
footprint, which estimates the first-dimension apex even though the peak
was sampled only every few seconds. t2 is the apex refined by a parabola
through the log of the three highest points.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field, replace

import numpy as np
from scipy import ndimage as ndi
from scipy.cluster.hierarchy import fclusterdata
from skimage.segmentation import watershed

from .background import NoiseMap, _robust_sigma, _scan_positions, storage_threshold
from .fold import Fold
from .run import Run


@dataclass(frozen=True)
class PeakSettings:
    snr_min: float = 10.0            # TIC blob seeds
    target_snr_min: float = 5.0      # quantifier-ion candidates
    floor_snr: float = 2.0           # footprints stop here
    rel_floor: float = 0.02          # ...or at this fraction of the apex, whichever is higher (targets)
    blob_rel_floor: float = 0.01     # the same for TIC blobs
    rim_baseline: bool = True        # subtract the level around a target peak (e.g. a UCM hump)
    grow_rows: int = 3               # widen target footprints by this many rows...
    grow_modulations: int = 1        # ...and modulations, within the peak's own basin
    smooth_t2_s: float = 0.03        # Gaussian SD along t2
    smooth_modulations: float = 0.6  # Gaussian SD along t1
    wrap_below_s: float = 0.5        # apexes below this t2 belong to the previous modulation
    min_modulations: int = 2
    min_pixels: int = 6
    noise_block: int = 20            # modulations per noise block (TIC image)
    ion_snr: bool = True             # also detect on the sum of per-ion S/N images, not the TIC alone
    min_ion_points: int = 15         # an ion needs at least this many stored points to help with detection
    coelution_frac: float = 0.05     # an ion needs at least this share of a blob's volume to test for co-elution
    coelution_gap: float = 1.2       # ion apexes further apart than this (see module docstring) are separate clusters
    chunk_modulations: int = 300     # split the watershed pass into pieces this wide, for speed on long runs
    chunk_overlap: int = 15          # ...with this much overlap either side, discarded at the seams
    chunk_workers: int = 4

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "PeakSettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


# ------------------------------------------------------------------ frame
@dataclass
class Frame:
    """The extended image geometry shared by everything in this module."""
    fold: Fold
    wrap_rows: int

    @property
    def n_rows(self) -> int:
        return self.fold.n_rows

    @property
    def rows(self) -> int:           # rows in the extended image
        return self.fold.n_rows + self.wrap_rows

    @property
    def pitch(self) -> float:
        return self.fold.row_pitch

    @property
    def period(self) -> float:
        return self.fold.spec.period

    @property
    def mid_t1(self) -> np.ndarray:
        return self.fold.t1 + 0.5 * self.period

    def extend(self, img: np.ndarray) -> np.ndarray:
        if self.wrap_rows == 0:
            return img
        return self.fold.extend(img, 1.0 + (self.wrap_rows + 0.25) / self.n_rows)[0][: self.rows]

    def canonical(self) -> np.ndarray:
        """Rows that hold each recorded sample exactly once."""
        m = np.zeros((self.rows, self.fold.n_cols), dtype=bool)
        m[self.wrap_rows:self.n_rows + self.wrap_rows] = True
        return m

    def scan_pixels(self) -> tuple[np.ndarray, np.ndarray]:
        """Extended-frame pixel (row, col) of every scan, in the canonical band."""
        col, row = _scan_positions(self.fold)
        wrapped = row < self.wrap_rows
        r = np.where(wrapped, row + self.n_rows, row)
        c = np.where(wrapped, col - 1, col)
        return r, c

    def time_to_pixel(self, t1_mid: float, t2: float) -> tuple[int, float]:
        """(column, extended row position) for a true t1 and t2."""
        col = int(np.floor((t1_mid - self.fold.t1[0]) / self.period))
        rowpos = t2 / self.pitch
        # A t2 beyond the extended band belongs to the next modulation's rows.
        while rowpos >= self.rows:
            rowpos -= self.n_rows
            col += 1
        while rowpos < self.wrap_rows:
            rowpos += self.n_rows
            col -= 1
        return col, rowpos


def make_frame(fold: Fold, settings: PeakSettings) -> Frame:
    wrap_rows = int(np.ceil(max(settings.wrap_below_s, 0.0) / fold.row_pitch))
    return Frame(fold=fold, wrap_rows=min(wrap_rows, fold.n_rows // 2))


def _smooth(img: np.ndarray, frame: Frame, s: PeakSettings) -> np.ndarray:
    sig_r = max(s.smooth_t2_s / frame.pitch, 0.0)
    sig_c = max(s.smooth_modulations, 0.0)
    if sig_r == 0 and sig_c == 0:
        return img
    return ndi.gaussian_filter(img, (sig_r, sig_c), mode="nearest")


def _parabolic(y0, y1, y2) -> float:
    if y0 > 0 and y1 > 0 and y2 > 0:
        y0, y1, y2 = np.log(y0), np.log(y1), np.log(y2)
    den = y0 - 2 * y1 + y2
    return float(np.clip(0.5 * (y0 - y2) / den, -0.5, 0.5)) if den < 0 else 0.0


def _column_noise(img: np.ndarray, frame: Frame, block: int, floor: float) -> np.ndarray:
    """Noise per column of an (extended) image, from its quiet pixels."""
    n_cols = img.shape[1]
    out = np.empty(n_cols)
    complete = frame.fold.complete
    for b0 in range(0, n_cols, block):
        cols = np.arange(b0, min(b0 + block, n_cols))
        use = cols[complete[cols]] if complete[cols].any() else cols
        out[cols] = _robust_sigma(img[:frame.n_rows, use].ravel(), floor)
    return out


# ------------------------------------------------------------------ blobs
@dataclass
class BlobSet:
    """Every TIC blob in one run."""
    mz: np.ndarray
    t1: np.ndarray          # s, first-dimension apex (intensity-weighted)
    t2: np.ndarray          # s, true second-dimension apex (may exceed the period)
    col: np.ndarray         # modulation of the apex
    row: np.ndarray         # extended-frame row of the apex
    height: np.ndarray      # smoothed TIC at the apex
    snr: np.ndarray
    volume: np.ndarray      # summed signal (all ions) over the footprint
    n_pixels: np.ndarray
    n_modulations: np.ndarray
    spectra: np.ndarray     # (n_blobs, n_mz) summed over the footprint
    labels: np.ndarray      # extended-frame label image; blob i has label i + 1
    wrap_rows: int
    pitch: float
    period: float
    t1_start: float         # start time of modulation 0
    settings: PeakSettings = field(default_factory=PeakSettings)
    width_t1: np.ndarray = field(default_factory=lambda: np.array([]))     # s, FWHM along t1
    width_t2: np.ndarray = field(default_factory=lambda: np.array([]))     # s, FWHM along t2
    apex_spectra: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))  # (n_blobs, n_mz), nearest-scan only
    co_elution: np.ndarray = field(default_factory=lambda: np.array([], dtype=bool))
    co_elution_clusters: np.ndarray = field(default_factory=lambda: np.array([], dtype=np.int32))

    def __len__(self) -> int:
        return int(self.t1.size)

    def ion_volume(self, i: int, ions) -> float:
        idx = np.asarray(list(ions), dtype=np.int64) - int(self.mz[0])
        return float(self.spectra[i, idx].sum())

    def footprint(self, i: int) -> tuple[np.ndarray, np.ndarray]:
        """Extended-frame (rows, cols) of blob i."""
        return np.nonzero(self.labels == i + 1)

    def outline_times(self, i: int) -> tuple[np.ndarray, np.ndarray]:
        """Outline of blob i as (t1 in s, t2 in s) points of pixel edges, for drawing."""
        rows, cols = self.footprint(i)
        if rows.size == 0:
            return np.array([]), np.array([])
        r0, r1, c0, c1 = rows.min(), rows.max(), cols.min(), cols.max()
        sub = (self.labels[r0:r1 + 1, c0:c1 + 1] == i + 1)
        pad = np.pad(sub, 1)
        xs, ys = [], []
        # Horizontal and vertical pixel edges where inside meets outside.
        dr = pad[1:, 1:-1] != pad[:-1, 1:-1]
        for rr, cc in zip(*np.nonzero(dr)):
            y = (r0 + rr - 0.5) * self.pitch
            xs += [c0 + cc, c0 + cc + 1, np.nan]
            ys += [y, y, np.nan]
        dc = pad[1:-1, 1:] != pad[1:-1, :-1]
        for rr, cc in zip(*np.nonzero(dc)):
            x = c0 + cc
            ys += [(r0 + rr - 0.5) * self.pitch, (r0 + rr + 0.5) * self.pitch, np.nan]
            xs += [x, x, np.nan]
        t1 = self.t1_start + np.asarray(xs, dtype=float) * self.period
        return t1, np.asarray(ys, dtype=float)

    def at_pixel(self, row: int, col: int) -> int | None:
        if 0 <= row < self.labels.shape[0] and 0 <= col < self.labels.shape[1]:
            lab = int(self.labels[row, col])
            return lab - 1 if lab > 0 else None
        return None

    def nearest(self, t1: float, t2: float, t1_scale: float = 60.0, t2_scale: float = 0.5) -> int | None:
        if len(self) == 0:
            return None
        d = np.hypot((self.t1 - t1) / t1_scale, (self.t2 - t2) / t2_scale)
        return int(np.argmin(d))


def _reindexed(arr: np.ndarray, idx: np.ndarray, n: int) -> np.ndarray:
    """`arr[idx]`, unless `arr` is one of the optional per-blob fields that
    an older or hand-built BlobSet may not carry (empty rather than length
    `n`), in which case it's left alone."""
    return arr[idx] if arr.shape[0] == n else arr


def delete_blobs(blobs: BlobSet, indices) -> BlobSet:
    """A copy of `blobs` with the given blobs removed and their footprint
    pixels cleared, everything else renumbered to stay contiguous."""
    idx = {int(i) for i in indices}
    keep = np.array([i not in idx for i in range(len(blobs))], dtype=bool)
    return _subset_blobs(blobs, keep)


def _subset_blobs(blobs: BlobSet, keep: np.ndarray) -> BlobSet:
    n = len(blobs)
    kept = np.flatnonzero(keep)
    remap = np.zeros(n + 1, dtype=np.int32)
    remap[kept + 1] = np.arange(1, kept.size + 1)
    return replace(blobs, t1=blobs.t1[kept], t2=blobs.t2[kept], col=blobs.col[kept], row=blobs.row[kept],
                    height=blobs.height[kept], snr=blobs.snr[kept], volume=blobs.volume[kept],
                    n_pixels=blobs.n_pixels[kept], n_modulations=blobs.n_modulations[kept],
                    spectra=blobs.spectra[kept], labels=remap[blobs.labels],
                    width_t1=_reindexed(blobs.width_t1, kept, n), width_t2=_reindexed(blobs.width_t2, kept, n),
                    apex_spectra=_reindexed(blobs.apex_spectra, kept, n),
                    co_elution=_reindexed(blobs.co_elution, kept, n),
                    co_elution_clusters=_reindexed(blobs.co_elution_clusters, kept, n))


def merge_blobs(blobs: BlobSet, indices) -> BlobSet:
    """Combine two or more blobs into one.

    Volume, pixel count, modulations spanned and the summed spectrum are
    exact sums: footprints never overlap, so nothing is double-counted.
    The combined apex position and width are exact too, from each blob's
    own position and width weighted by its volume (the parallel-axis
    theorem: combining two known distributions doesn't need the original
    pixels back). The apex, height and S/N are taken from whichever
    contributing blob was tallest. Co-elution isn't re-tested; run 🔍
    again to refresh that flag on a merged peak.
    """
    idx = np.array(sorted({int(i) for i in indices}))
    if idx.size < 2:
        raise ValueError("need at least two blobs to merge")
    n = len(blobs)
    keep = np.ones(n, dtype=bool)
    keep[idx] = False
    rest = _subset_blobs(blobs, keep)

    vol = np.clip(blobs.volume[idx], 1e-12, None)
    t1 = float((vol * blobs.t1[idx]).sum() / vol.sum())
    t2 = float((vol * blobs.t2[idx]).sum() / vol.sum())

    def combined_width(pos: np.ndarray, width: np.ndarray, centre: float) -> float:
        if width.shape[0] != n:
            return 0.0
        var = (width[idx] / 2.3548) ** 2
        total = (vol * (var + (pos[idx] - centre) ** 2)).sum() / vol.sum()
        return float(2.3548 * np.sqrt(np.clip(total, 0, None)))

    width_t1 = combined_width(blobs.t1, blobs.width_t1, t1)
    width_t2 = combined_width(blobs.t2, blobs.width_t2, t2)
    apex = idx[int(np.argmax(blobs.height[idx]))]

    union = np.isin(blobs.labels, (idx + 1).tolist())
    labels = rest.labels.copy()
    labels[union] = len(rest) + 1

    def append(arr, extra):
        return np.concatenate([arr, [extra]]) if arr.shape[0] == len(rest) else arr

    def append_row(arr, extra):
        return np.vstack([arr, extra[None, :]]) if arr.shape[0] == len(rest) else arr

    return replace(rest, t1=np.append(rest.t1, t1), t2=np.append(rest.t2, t2),
                    col=np.append(rest.col, blobs.col[apex]), row=np.append(rest.row, blobs.row[apex]),
                    height=np.append(rest.height, blobs.height[apex]), snr=np.append(rest.snr, blobs.snr[apex]),
                    volume=np.append(rest.volume, float(blobs.volume[idx].sum())),
                    n_pixels=np.append(rest.n_pixels, int(union.sum())),
                    n_modulations=np.append(rest.n_modulations, int(np.unique(np.nonzero(union)[1]).size)),
                    spectra=np.vstack([rest.spectra, blobs.spectra[idx].sum(axis=0, keepdims=True)]),
                    labels=labels, width_t1=append(rest.width_t1, width_t1),
                    width_t2=append(rest.width_t2, width_t2),
                    apex_spectra=append_row(rest.apex_spectra, blobs.apex_spectra[apex])
                    if blobs.apex_spectra.shape[0] == n else rest.apex_spectra,
                    co_elution=append(rest.co_elution, False), co_elution_clusters=append(rest.co_elution_clusters, 1))


def add_blob_at(blobs: BlobSet, run: Run, fold: Fold, t1: float, t2: float) -> BlobSet:
    """A copy of `blobs` with one more, grown from the clicked pixel itself
    (refined to the tallest TIC pixel nearby) rather than a local maximum
    found automatically, for a real peak that automatic detection missed
    (too faint, or shadowed by a bigger neighbour). Raises `ValueError` if
    there's nothing there, or the click landed on an existing blob."""
    s = blobs.settings
    frame = make_frame(fold, s)
    raw = frame.extend(fold.tic_image(run))
    sm = _smooth(raw, frame, s)
    col, rpos = frame.time_to_pixel(t1, t2)
    row = int(round(np.clip(rpos, 0, sm.shape[0] - 1)))
    col = int(np.clip(col, 0, sm.shape[1] - 1))
    r0, r1 = max(row - 3, 0), min(row + 4, sm.shape[0])
    c0, c1 = max(col - 1, 0), min(col + 2, sm.shape[1])
    rr, cc = np.unravel_index(int(np.argmax(sm[r0:r1, c0:c1])), (r1 - r0, c1 - c0))
    row, col = r0 + rr, c0 + cc
    if blobs.labels.shape == sm.shape and blobs.labels[row, col] != 0:
        raise ValueError("that pixel already belongs to a blob")
    peak = float(sm[row, col])
    if peak <= 0:
        raise ValueError("nothing there to measure")
    canon = frame.canonical()
    markers = np.zeros(sm.shape, dtype=np.int32)
    markers[row, col] = 1
    free = blobs.labels == 0 if blobs.labels.shape == sm.shape else np.ones(sm.shape, dtype=bool)
    mask = canon & free & (sm >= s.blob_rel_floor * peak)
    fp = watershed(-sm, markers, mask=mask) == 1
    if not fp.any():
        raise ValueError("nothing there to measure")

    n = len(blobs)
    new_labels = blobs.labels.copy() if blobs.labels.shape == sm.shape else np.zeros(sm.shape, dtype=np.int32)
    new_labels[fp] = n + 1
    spectra, _ = _spectra_by_label(run, new_labels, frame, n + 1)
    apex_labels = np.zeros_like(new_labels)
    apex_labels[row, col] = n + 1
    apex_spectra, _ = _spectra_by_label(run, apex_labels, frame, n + 1)

    mid = frame.mid_t1
    colvol = np.bincount(np.nonzero(fp)[1], weights=np.clip(raw, 0, None)[fp], minlength=sm.shape[1])
    colsum = max(float(colvol.sum()), 1e-12)
    new_t1 = float((colvol * mid).sum() / colsum)
    new_width_t1 = 2.3548 * float(np.sqrt(np.clip((colvol * (mid - new_t1) ** 2).sum() / colsum, 0, None)))
    d = _parabolic(sm[row - 1, col], sm[row, col], sm[row + 1, col]) if 0 < row < sm.shape[0] - 1 else 0.0
    new_t2 = (row + d) * frame.pitch
    t2_pos = np.arange(sm.shape[0]) * frame.pitch
    rowvol = np.bincount(np.nonzero(fp)[0], weights=np.clip(raw, 0, None)[fp], minlength=sm.shape[0])
    rowsum = max(float(rowvol.sum()), 1e-12)
    new_width_t2 = 2.3548 * float(np.sqrt(np.clip((rowvol * (t2_pos - new_t2) ** 2).sum() / rowsum, 0, None)))

    def append(arr, extra):
        return np.append(arr, extra) if arr.shape[0] == n else arr

    return replace(blobs, t1=np.append(blobs.t1, new_t1), t2=np.append(blobs.t2, new_t2),
                    col=np.append(blobs.col, col), row=np.append(blobs.row, row),
                    height=np.append(blobs.height, peak), snr=np.append(blobs.snr, peak),
                    volume=np.append(blobs.volume, float(spectra[-1].sum())),
                    n_pixels=np.append(blobs.n_pixels, int(fp.sum())),
                    n_modulations=np.append(blobs.n_modulations, int((colvol > 0).sum())),
                    spectra=spectra.astype(np.float32), labels=new_labels,
                    width_t1=append(blobs.width_t1, new_width_t1), width_t2=append(blobs.width_t2, new_width_t2),
                    apex_spectra=(np.vstack([blobs.apex_spectra, apex_spectra[-1:].astype(np.float32)])
                                 if blobs.apex_spectra.shape[0] == n else blobs.apex_spectra),
                    co_elution=append(blobs.co_elution, False), co_elution_clusters=append(blobs.co_elution_clusters, 1))


def _spectra_by_label(run: Run, labels: np.ndarray, frame: Frame, n: int):
    r, c = frame.scan_pixels()
    ok = (c >= 0) & (c < labels.shape[1])
    scan_label = np.zeros(run.n_scans, dtype=np.int64)
    scan_label[ok] = labels[r[ok], c[ok]]
    coo = run.data.tocoo()
    lab = scan_label[coo.row]
    keep = lab > 0
    flat = (lab[keep] - 1) * run.n_mz + coo.col[keep]
    spectra = np.bincount(flat, weights=coo.data[keep].astype(np.float64),
                          minlength=n * run.n_mz).reshape(n, run.n_mz)
    return spectra, scan_label


def _robust_sigma_batch(vals: np.ndarray, floor: float, iterations: int = 4) -> np.ndarray:
    """`_robust_sigma`, for every column of a 2D array at once (rows are the
    pixels of one t1 block, columns are independent ions)."""
    v = vals.astype(np.float64)
    n = v.shape[0]
    frac_nonzero = np.mean(v != 0, axis=0)
    med = np.median(v, axis=0)
    sig = 1.4826 * np.median(np.abs(v - med[None, :]), axis=0)
    for _ in range(iterations):
        bad = sig <= 0
        keep = v < (med[None, :] + 3.0 * sig[None, :])
        any_keep = keep.any(axis=0)
        # A column with nothing kept would be all-NaN below; numpy warns
        # about that even though such a column's result is discarded via
        # `bad`/`any_keep` regardless, so fall back to the raw column there
        # instead (this runs in worker threads, where a racy warnings
        # filter wouldn't reliably suppress it anyway).
        vv = np.where(keep | ~any_keep[None, :], v, np.nan)
        new_med = np.nanmedian(vv, axis=0)
        new_sig = 1.4826 * np.nanmedian(np.abs(vv - new_med[None, :]), axis=0)
        skip = bad | ~any_keep
        med = np.where(skip, med, new_med)
        sig = np.where(skip, sig, new_sig)
    sig = np.where(np.isnan(sig) | (n < 10) | (frac_nonzero < 0.5), floor, sig)
    return np.maximum(sig, floor)


def _column_noise_batch(cube: np.ndarray, frame: Frame, block: int, floor: float) -> np.ndarray:
    """`_column_noise`, for every ion in a cube at once. Returns (n_cols, n_ions)."""
    n_cols, n_ions = cube.shape[1], cube.shape[2]
    out = np.empty((n_cols, n_ions))
    complete = frame.fold.complete
    for b0 in range(0, n_cols, block):
        cols = np.arange(b0, min(b0 + block, n_cols))
        use = cols[complete[cols]] if complete[cols].any() else cols
        vals = cube[:frame.n_rows][:, use, :].reshape(-1, n_ions)
        out[cols] = _robust_sigma_batch(vals, floor)
    return out


def _active_ion_cube(run: Run, fold: Fold, frame: Frame, s: PeakSettings):
    """Extended-frame cube of every ion with enough stored points to be
    worth including in the detection surface (see the module docstring)."""
    counts = np.diff(run.data.tocsc().indptr)
    idx = np.flatnonzero(counts >= max(int(s.min_ion_points), 1))
    if idx.size == 0:
        return np.array([], dtype=np.int64), None
    ions = run.mz[idx]
    cube = fold.cube(run, mz=ions, dtype=np.float32)
    return ions, frame.extend(cube)


def _ion_snr_cube(cube: np.ndarray | None, frame: Frame, s: PeakSettings, floor: float):
    """Smoothed per-ion cube and the sum of its per-ion S/N images, or
    (None, None) if there are no ions to work with."""
    if cube is None or cube.shape[-1] == 0:
        return None, None
    sig_r = max(s.smooth_t2_s / frame.pitch, 0.0)
    sig_c = max(s.smooth_modulations, 0.0)
    sm = ndi.gaussian_filter(cube, (sig_r, sig_c, 0.0), mode="nearest") if (sig_r or sig_c) else cube
    noise = _column_noise_batch(sm, frame, s.noise_block, floor)
    total = (sm / noise[None, :, :]).sum(axis=2)
    return sm, total


def _chunk_bounds(n_cols: int, size: int, overlap: int):
    """(padded_lo, padded_hi, core_lo, core_hi) for every piece of a long
    image: the core columns partition the run exactly once; the padding
    either side gives a chunk's watershed enough context that no real peak
    is cut off at the seam. One chunk (no padding) if the run is short."""
    if n_cols <= size:
        yield 0, n_cols, 0, n_cols
        return
    start = 0
    while start < n_cols:
        core_hi = min(start + size, n_cols)
        yield max(start - overlap, 0), min(core_hi + overlap, n_cols), start, core_hi
        start = core_hi


def _segment_chunk(det: np.ndarray, sm: np.ndarray, canon: np.ndarray, s: PeakSettings):
    """Seeds and a label image for one piece of the detection surface.
    Exactly what a single, unchunked run of `find_blobs` used to do."""
    peaks = (det == ndi.maximum_filter(det, size=(5, 3), mode="nearest")) & (det >= s.snr_min) & canon
    seeds = np.argwhere(peaks)
    markers = np.zeros(det.shape, dtype=np.int32)
    markers[seeds[:, 0], seeds[:, 1]] = np.arange(1, len(seeds) + 1)
    mask = canon & (det >= s.floor_snr)
    labels = watershed(-det, markers, mask=mask) if len(seeds) else np.zeros(det.shape, np.int32)
    if len(seeds) and s.blob_rel_floor > 0:
        # Trim each footprint to a fraction of its own TIC apex, so big peaks
        # don't spread across low-level noise, then flood again so every kept
        # pixel stays connected to its seed.
        heights = np.concatenate([[np.inf], sm[seeds[:, 0], seeds[:, 1]]])
        mask = mask & (sm >= s.blob_rel_floor * heights[labels])
        labels = watershed(-det, markers, mask=mask)
        labels = _trim_labels_to_own_height(labels, seeds, sm, s.blob_rel_floor)
    return seeds, labels


def _trim_labels_to_own_height(labels: np.ndarray, seeds: np.ndarray,
                               intensity: np.ndarray, relative_floor: float) -> np.ndarray:
    """Enforce each final basin's own relative floor and connectivity.

    A second watershed can assign a pixel to a different marker from the one
    whose threshold admitted it to the mask.  Check the final ownership, then
    retain only the connected part containing that marker.  This is deliberately
    the last operation: another global flood could recreate the same ownership
    mismatch.
    """
    if relative_floor <= 0 or len(seeds) == 0:
        return labels
    out = np.zeros_like(labels)
    structure = ndi.generate_binary_structure(2, 1)
    for j, (row, col) in enumerate(seeds, start=1):
        allowed = (labels == j) & (intensity >= relative_floor * intensity[row, col])
        components, _ = ndi.label(allowed, structure=structure)
        component = int(components[row, col])
        if component:
            out[components == component] = j
    return out


def _coelution_flags(labels: np.ndarray, seeds: np.ndarray, spectra: np.ndarray, run: Run,
                     ion_mz: np.ndarray, sm_cube: np.ndarray | None, s: PeakSettings):
    """Per-blob co-elution flag and cluster count. See the module docstring."""
    n = len(seeds)
    flags = np.zeros(n, dtype=bool)
    clusters = np.ones(n, dtype=np.int32)
    if sm_cube is None or n == 0 or ion_mz.size < 2:
        return flags, clusters
    mz0 = int(run.mz[0])
    for i in range(n):
        rows, cols = np.nonzero(labels == i + 1)
        if rows.size == 0:
            continue
        r0, r1, c0, c1 = rows.min(), rows.max(), cols.min(), cols.max()
        subfp = (labels[r0:r1 + 1, c0:c1 + 1] == i + 1)
        total = float(spectra[i].sum())
        if total <= 0:
            continue
        pts = []
        for k, ion in enumerate(ion_mz):
            idx = int(ion) - mz0
            if not (0 <= idx < spectra.shape[1]) or spectra[i, idx] / total < s.coelution_frac:
                continue
            sub = sm_cube[r0:r1 + 1, c0:c1 + 1, k]
            vals = np.where(subfp, sub, -np.inf)
            if not np.any(np.isfinite(vals)) or float(np.max(vals)) <= 0:
                continue
            rr, cc = np.unravel_index(int(np.argmax(vals)), vals.shape)
            # Scaled so one modulation along t1 and two rows along t2 both
            # count as one unit of distance.
            pts.append((float(cc + c0), float(rr + r0) / 2.0))
        if len(pts) >= 2:
            cl = fclusterdata(np.asarray(pts), t=s.coelution_gap, criterion="distance")
            clusters[i] = int(cl.max())
            flags[i] = clusters[i] > 1
    return flags, clusters


def find_blobs(run: Run, fold: Fold, settings: PeakSettings = PeakSettings(),
               floor: float | None = None) -> BlobSet:
    """TIC blobs of a (preferably cleaned) run. See the module docstring."""
    s = settings
    frame = make_frame(fold, s)
    floor = storage_threshold(run) / 3.0 if floor is None else floor
    raw = frame.extend(fold.tic_image(run))
    sm = _smooth(raw, frame, s)
    tic_noise = _column_noise(sm, frame, s.noise_block, floor)
    det = sm / tic_noise[None, :]
    canon = frame.canonical()

    ion_mz, sm_cube = np.array([], dtype=np.int64), None
    if s.ion_snr:
        ion_mz, cube = _active_ion_cube(run, fold, frame, s)
        sm_cube, ion_sum = _ion_snr_cube(cube, frame, s, floor)
        if ion_sum is not None:
            det = np.maximum(det, ion_sum)

    bounds = list(_chunk_bounds(det.shape[1], max(int(s.chunk_modulations), 1), max(int(s.chunk_overlap), 0)))
    if len(bounds) == 1:
        pieces = [_segment_chunk(det, sm, canon, s)]
    else:
        with ThreadPoolExecutor(max_workers=max(1, int(s.chunk_workers))) as ex:
            pieces = list(ex.map(
                lambda b: _segment_chunk(det[:, b[0]:b[1]], sm[:, b[0]:b[1]], canon[:, b[0]:b[1]], s), bounds))

    # Stitch the chunks together: each peak is kept exactly once, by
    # whichever chunk owns the column its seed falls in.
    labels = np.zeros(det.shape, dtype=np.int32)
    seeds_list = []
    next_id = 0
    for (lo, _hi, core_lo, core_hi), (loc_seeds, loc_labels) in zip(bounds, pieces):
        for j in range(len(loc_seeds)):
            r, c = int(loc_seeds[j, 0]), int(loc_seeds[j, 1])
            gc = c + lo
            if not (core_lo <= gc < core_hi):
                continue
            fr, fc = np.nonzero(loc_labels == j + 1)
            next_id += 1
            labels[fr, fc + lo] = next_id
            seeds_list.append((r, gc))
    seeds = np.array(seeds_list, dtype=np.int64).reshape(-1, 2) if seeds_list else np.zeros((0, 2), dtype=np.int64)

    n = len(seeds)
    n_rows_img, n_cols_img = labels.shape
    npix = np.bincount(labels.ravel(), minlength=n + 1)[1:]
    cols_idx = np.broadcast_to(np.arange(n_cols_img), labels.shape)
    colvol = np.bincount((labels * n_cols_img + cols_idx).ravel(), weights=np.clip(raw, 0, None).ravel(),
                         minlength=(n + 1) * n_cols_img).reshape(n + 1, n_cols_img)[1:]
    rows_idx = np.broadcast_to(np.arange(n_rows_img)[:, None], labels.shape)
    rowvol = np.bincount((labels * n_rows_img + rows_idx).ravel(), weights=np.clip(raw, 0, None).ravel(),
                         minlength=(n + 1) * n_rows_img).reshape(n + 1, n_rows_img)[1:]
    nmod = (colvol > 0).sum(axis=1)
    keep = (nmod >= s.min_modulations) & (npix >= s.min_pixels)
    # Relabel so kept blobs are 1..k
    new_id = np.zeros(n + 1, dtype=np.int32)
    new_id[1:][keep] = np.arange(1, keep.sum() + 1)
    labels = new_id[labels]
    seeds, colvol, rowvol, npix, nmod = seeds[keep], colvol[keep], rowvol[keep], npix[keep], nmod[keep]
    k = len(seeds)

    spectra, _ = _spectra_by_label(run, labels, frame, k)
    apex_labels = np.zeros_like(labels)
    apex_labels[seeds[:, 0], seeds[:, 1]] = np.arange(1, k + 1)
    apex_spectra, _ = _spectra_by_label(run, apex_labels, frame, k)

    mid = frame.mid_t1
    colsum = np.maximum(colvol.sum(axis=1), 1e-12)
    t1 = (colvol * mid[None, :]).sum(axis=1) / colsum
    width_t1 = 2.3548 * np.sqrt(np.clip((colvol * (mid[None, :] - t1[:, None]) ** 2).sum(axis=1) / colsum, 0, None))

    t2 = np.empty(k)
    for j, (r, c) in enumerate(seeds):
        col = sm[:, c]
        d = _parabolic(col[r - 1], col[r], col[r + 1]) if 0 < r < len(col) - 1 else 0.0
        t2[j] = (r + d) * frame.pitch
    t2_pos = np.arange(n_rows_img) * frame.pitch
    rowsum = np.maximum(rowvol.sum(axis=1), 1e-12)
    width_t2 = 2.3548 * np.sqrt(np.clip((rowvol * (t2_pos[None, :] - t2[:, None]) ** 2).sum(axis=1) / rowsum, 0, None))

    co_elution, co_clusters = _coelution_flags(labels, seeds, spectra, run, ion_mz, sm_cube, s)

    return BlobSet(mz=run.mz.copy(), t1=t1, t2=t2, col=seeds[:, 1].copy(), row=seeds[:, 0].copy(),
                   height=sm[seeds[:, 0], seeds[:, 1]], snr=det[seeds[:, 0], seeds[:, 1]],
                   volume=spectra.sum(axis=1), n_pixels=npix, n_modulations=nmod,
                   spectra=spectra.astype(np.float32), labels=labels, wrap_rows=frame.wrap_rows,
                   pitch=frame.pitch, period=frame.period, t1_start=float(fold.t1[0]), settings=s,
                   width_t1=width_t1, width_t2=width_t2, apex_spectra=apex_spectra.astype(np.float32),
                   co_elution=co_elution, co_elution_clusters=co_clusters)


# ------------------------------------------------------------- ion peaks
@dataclass
class IonPeak:
    """One candidate peak on a target's quantifier ion."""
    t1: float
    t2: float
    col: int
    row: int
    snr: float
    quant_volume: float
    qual_volumes: tuple[float, ...]
    max_raw: float                 # largest raw quantifier value inside the footprint
    n_pixels: int
    n_modulations: int
    pixels: tuple[np.ndarray, np.ndarray]   # extended-frame (rows, cols)
    baseline: float = 0.0          # level subtracted per pixel (from the rim)
    touches_edge: bool = False     # footprint reaches the wrap boundary, so part may be cut off
    touches_larger: bool = False   # shares a border with a bigger peak, so the split is uncertain

    @property
    def qual_ratios(self) -> tuple[float, ...]:
        q = self.quant_volume
        return tuple(v / q if q > 0 else np.nan for v in self.qual_volumes)


def measure_ion_peak(run: Run, raw_run: Run, fold: Fold, noise: NoiseMap | None, quant_ion: int,
                     qual_ions, t1: float, t2: float, t1_window: float, t2_window: float,
                     settings: PeakSettings = PeakSettings(),
                     frame: Frame | None = None, cache: dict | None = None) -> list[IonPeak]:
    """Candidate peaks of one ion inside a target's window. See the module docstring."""
    s = settings
    frame = frame or make_frame(fold, s)
    cache = {} if cache is None else cache

    def ext_image(r: Run, ion: int, key: str):
        k = (key, ion)
        if k not in cache:
            cache[k] = frame.extend(fold.image(r.xic(ion)))
        return cache[k]

    img = ext_image(run, quant_ion, "clean")
    c_exp, rpos = frame.time_to_pixel(t1, t2)
    half_c = int(np.ceil(t1_window / frame.period)) + 1
    half_r = int(np.ceil(t2_window / frame.pitch))
    pad_c, pad_r = half_c + 4, half_r + int(np.ceil(0.8 / frame.pitch))
    r_lo, r_hi = max(int(rpos) - pad_r, 0), min(int(rpos) + pad_r + 1, frame.rows)
    c_lo, c_hi = max(c_exp - pad_c, 0), min(c_exp + pad_c + 1, img.shape[1])
    if r_hi - r_lo < 3 or c_hi - c_lo < 1:
        return []
    sub = img[r_lo:r_hi, c_lo:c_hi]
    sm = _smooth(sub, frame, s)
    if "floor" not in cache:
        cache["floor"] = storage_threshold(raw_run) / 3.0
    floor = cache["floor"]
    if noise is not None:
        sig = noise.at(quant_ion, frame.mid_t1[c_lo:c_hi])
    else:
        sig = np.full(c_hi - c_lo, floor)
    sig = np.maximum(sig, 1e-9)
    snr = sm / sig[None, :]
    canon = frame.canonical()[r_lo:r_hi, c_lo:c_hi]
    maxima = (sm == ndi.maximum_filter(sm, size=(5, 3), mode="nearest")) & (snr >= s.target_snr_min) & canon
    seeds = np.argwhere(maxima)
    if len(seeds) == 0:
        return []
    markers = np.zeros(sm.shape, dtype=np.int32)
    markers[seeds[:, 0], seeds[:, 1]] = np.arange(1, len(seeds) + 1)
    labels = watershed(-sm, markers, mask=canon & (snr >= s.floor_snr))
    basins = watershed(-sm, markers, mask=canon)

    # Candidates are seeds inside the window proper.
    t1_mid = frame.mid_t1
    quals = [ext_image(run, q, "clean")[r_lo:r_hi, c_lo:c_hi] for q in qual_ions]
    out = []
    for j, (r, c) in enumerate(seeds):
        cc, rr = c + c_lo, r + r_lo
        if abs(cc - c_exp) > half_c or abs(rr - rpos) > half_r + 1:
            continue
        fp = (labels == j + 1) & (sm >= s.rel_floor * sm[r, c])
        fp_lab = ndi.label(fp)[0]
        fp = fp_lab == fp_lab[r, c] if fp_lab[r, c] else fp
        # The noise-level cut trims the tails of small peaks. Cleaned noise
        # averages to zero, so widening the footprint a little (within the
        # peak's own basin) recovers the tails without adding bias.
        if s.grow_rows or s.grow_modulations:
            st = np.ones((2 * s.grow_rows + 1, 2 * s.grow_modulations + 1), bool)
            fp = ndi.binary_dilation(fp, structure=st) & (basins == j + 1)
        base_q = _rim_level(sub, fp, canon) if s.rim_baseline else 0.0
        colsum = np.clip(sub - base_q, 0, None)
        colsum[~fp] = 0
        cv = colsum.sum(axis=0)
        if cv.sum() <= 0:
            continue
        pt1 = float((cv * t1_mid[c_lo:c_hi]).sum() / cv.sum())
        col = sm[:, c]
        d = _parabolic(col[r - 1], col[r], col[r + 1]) if 0 < r < len(col) - 1 else 0.0
        pt2 = (rr + d) * frame.pitch
        qv = tuple(float((q[fp] - (_rim_level(q, fp, canon) if s.rim_baseline else 0.0)).sum())
                   for q in quals)
        rows, cols = np.nonzero(fp)
        max_raw = float(_max_raw(raw_run, frame, quant_ion, rows + r_lo, cols + c_lo, cache)) \
            if raw_run is not None else float("nan")
        out.append(IonPeak(t1=pt1, t2=pt2, col=int(cc), row=int(rr), snr=float(snr[r, c]),
                           quant_volume=float((sub[fp] - base_q).sum()), qual_volumes=qv, max_raw=max_raw,
                           baseline=float(base_q), touches_edge=bool(_touches_edge(fp, canon)),
                           touches_larger=bool(_touches_larger(fp, labels, seeds, sm, j)),
                           n_pixels=int(fp.sum()), n_modulations=int((cv > 0).sum()),
                           pixels=(rows + r_lo, cols + c_lo)))
    return out


def measure_ion_polygon(run: Run, raw_run: Run, fold: Fold, noise: NoiseMap | None, quant_ion: int,
                        qual_ions, polygon, settings: PeakSettings = PeakSettings(),
                        frame: Frame | None = None, cache: dict | None = None) -> IonPeak | None:
    """Integrate one ion inside a hand-drawn polygon instead of the automatic
    footprint. `polygon` is a list of (t1, t2) vertices in seconds, t1 in the
    2D view's convention (the start of a modulation is `fold.t1[c]`, so the
    modulation's middle is half a period later). The volume is the sum over the
    pixels inside, less the same rim baseline the automatic footprints use.
    `None` if the polygon holds no pixels or no signal."""
    from skimage.draw import polygon as _fill
    s = settings
    frame = frame or make_frame(fold, s)
    cache = {} if cache is None else cache

    def ext_image(r: Run, ion: int, key: str):
        k = (key, ion)
        if k not in cache:
            cache[k] = frame.extend(fold.image(r.xic(ion)))
        return cache[k]

    rr, cc = [], []
    for t1, t2 in polygon:
        c = (t1 - frame.fold.t1[0]) / frame.period - 0.5
        r = t2 / frame.pitch
        if r < frame.wrap_rows:            # the wrapped tail sits below the next modulation's rows
            r += frame.n_rows
            c -= 1.0
        rr.append(r)
        cc.append(c)
    img = ext_image(run, quant_ion, "clean")
    canon = frame.canonical()
    fp = np.zeros(img.shape, dtype=bool)
    ri, ci = _fill(np.asarray(rr), np.asarray(cc), shape=img.shape)
    fp[ri, ci] = True
    fp &= canon
    if not fp.any():
        return None
    base = _rim_level(img, fp, canon) if s.rim_baseline else 0.0
    vals = np.clip(img - base, 0, None) * fp
    if vals.sum() <= 0:
        return None
    rows, cols = np.nonzero(fp)
    cv = vals.sum(axis=0)
    t1 = float((cv * frame.mid_t1[: cv.size]).sum() / cv.sum())
    rv = vals.sum(axis=1)
    t2 = float((rv * np.arange(rv.size)).sum() / rv.sum() * frame.pitch)
    apex = np.unravel_index(int(np.argmax(np.where(fp, img, -np.inf))), img.shape)
    quals = tuple(float((ext_image(run, q, "clean")[fp]
                         - (_rim_level(ext_image(run, q, "clean"), fp, canon) if s.rim_baseline else 0.0)).sum())
                  for q in qual_ions)
    if "floor" not in cache:
        cache["floor"] = storage_threshold(raw_run) / 3.0
    sig = float(noise.at(quant_ion, frame.mid_t1[[apex[1]]])[0]) if noise is not None else cache["floor"]
    max_raw = float(_max_raw(raw_run, frame, quant_ion, rows, cols, cache)) if raw_run is not None else float("nan")
    return IonPeak(t1=t1, t2=t2, col=int(apex[1]), row=int(apex[0]), snr=float(img[apex] / max(sig, 1e-9)),
                   quant_volume=float((img[fp] - base).sum()), qual_volumes=quals, max_raw=max_raw,
                   n_pixels=int(fp.sum()), n_modulations=int((cv > 0).sum()), pixels=(rows, cols),
                   baseline=float(base), touches_edge=bool(_touches_edge(fp, canon)), touches_larger=False)


def _rim_level(img: np.ndarray, fp: np.ndarray, canon: np.ndarray) -> float:
    """Level around a footprint: a low percentile of the ring of pixels just
    outside it. A low percentile ignores neighbouring peaks on the ring;
    clipping at zero avoids inflating peaks on cleaned data, where the ring
    is only noise."""
    ring = ndi.binary_dilation(fp, structure=np.ones((5, 3), bool)) & ~fp & canon
    if ring.sum() < 4:
        return 0.0
    return max(float(np.percentile(img[ring], 25)), 0.0)


def _touches_larger(fp, labels, seeds, sm, j) -> bool:
    ring = ndi.binary_dilation(fp, structure=np.ones((3, 3), bool)) & ~fp
    neighbours = np.unique(labels[ring])
    own = sm[seeds[j, 0], seeds[j, 1]]
    for lab in neighbours:
        if lab > 0 and lab != j + 1:
            r, c = seeds[lab - 1]
            if sm[r, c] > own:
                return True
    return False


def _touches_edge(fp: np.ndarray, canon: np.ndarray) -> bool:
    outside = ndi.binary_dilation(fp, structure=np.ones((3, 1), bool)) & ~canon
    return bool(outside.any())


def _max_raw(raw_run: Run, frame: Frame, ion: int, rows, cols, cache) -> float:
    """Largest raw value of one ion among the scans inside a footprint."""
    key = ("scanpix",)
    if key not in cache:
        cache[key] = frame.scan_pixels()
    r, c = cache[key]
    kx = ("rawxic", ion)
    if kx not in cache:
        cache[kx] = raw_run.xic(ion)
    x = cache[kx]
    cmin, cmax = cols.min(), cols.max()
    sel = np.flatnonzero((c >= cmin) & (c <= cmax))
    if sel.size == 0:
        return 0.0
    fp = set(zip(rows.tolist(), cols.tolist()))
    inside = [i for i in sel if (r[i], c[i]) in fp]
    return float(x[inside].max()) if inside else 0.0


def detect_ceiling(run: Run, min_count: int = 3) -> float | None:
    """The detector's saturation value, if the data show one: several stored
    points sharing exactly the largest value."""
    v = run.data.data
    if v.size == 0:
        return None
    top = float(v.max())
    return top if int(np.sum(v >= top * (1 - 1e-7))) >= min_count else None
