"""Regions of a chromatogram found from the pixels themselves, with no peak detection.

The idea in plain terms
------------------------
`cluster.py` groups detected peaks. This groups the chromatogram: every
scan is one pixel of the folded image and carries a mass spectrum, so the
image is a stack of spectra. Pixels whose spectra look alike are put in the
same region, and the result is a label image the same shape as the
chromatogram (`cluster.plot_cluster_chromatogram` draws it). Because it
never asks where a peak starts or ends, it works on the parts blob detection
handles badly: an unresolved hump, a saturated solvent front, overlapping peaks.

Method (`fit_regions`)
-----------------------
1. Scans brighter than a percentile of the run's total signal are kept (the
   rest are baseline). Their spectra are cut to the usable m/z range,
   weighted (intensity^a * m/z^b) and scaled to unit length, exactly as for
   peaks.
2. The spectra are reduced by truncated SVD to `n_components` directions
   (scikit-learn), scaled to unit length again, and grouped by k-means into
   `n_regions` (Lloyd, k-means++ starts, fixed seed).
3. Each kept scan's region is painted on its pixel (the most common region
   where a pixel holds several scans), then a majority filter of radius
   `smooth` pixels removes speckle without growing the regions into the
   unlabelled baseline.
4. Regions are numbered by the signal they hold, largest first.

The fitted model labels any other run (`SegmentResult.label_image`), so a
region means the same thing in every run, and `region_shares` gives each
run's share of its signal in each region, a small table to compare runs by.
Fit on several runs at once by giving `fit_regions` all of them.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import scipy.sparse as sp
from scipy import ndimage as ndi

from .background import _scan_positions
from .fold import Fold
from .run import Run


@dataclass(frozen=True)
class SegmentSettings:
    n_regions: int = 12
    mz_min: int = 35                   # below this, air, water, CO2 and solvent dominate the spectrum
    exclude_ions: tuple = ()           # e.g. SILOXANE_IONS
    weight_a: float = 0.5              # as ClusterSettings
    weight_b: float = 1.0
    keep_percentile: float = 70.0      # scans below this percentile of the run's total signal are baseline
    max_scans: int = 150_000           # fit on at most this many scans (the brightest), across all runs given
    n_components: int = 20             # SVD directions kept before k-means
    smooth: int = 1                    # radius in pixels of the majority filter; 0 leaves the raw labelling
    seed: int = 0

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "SegmentSettings":
        d = {k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__}
        if "exclude_ions" in d:
            d["exclude_ions"] = tuple(d["exclude_ions"])
        return cls(**d)


@dataclass
class RegionSummary:
    label: int
    n_pixels: int
    percent: float            # of the signal in labelled pixels
    t1_median: float          # s
    t2_median: float          # s
    top_ions: list            # m/z with the largest mean share of the ion current


@dataclass
class SegmentResult:
    settings: SegmentSettings
    mz: np.ndarray            # m/z axis of the run(s) it was fitted on
    mz_used: np.ndarray       # the ions the distance is built from
    svd: object               # the fitted TruncatedSVD
    centroids: np.ndarray     # (k, n_components) unit-length region centres
    mean_spectra: np.ndarray  # (k, len(mz)) mean share of ion current, for reading a region

    def __len__(self) -> int:
        return int(self.centroids.shape[0])

    def top_ions(self, k: int, n: int = 6) -> list:
        return [int(m) for m in self.mz[np.argsort(self.mean_spectra[k])[::-1][:n]]]

    def label_image(self, run: Run, fold: Fold) -> np.ndarray:
        """Region of every pixel of `run`, (rows, modulations); -1 for baseline pixels."""
        M, idx, _ = _scan_spectra(run, self.settings, self.mz_used, self.settings.keep_percentile)
        labels = self._nearest(M)
        return _paint(labels, idx, fold, len(self), self.settings.smooth)

    def _nearest(self, M) -> np.ndarray:
        Z = _reduce(self.svd, M)
        return np.argmax(Z @ self.centroids.T, axis=1)

    def summary(self, labels: np.ndarray, run: Run, fold: Fold) -> list:
        """One `RegionSummary` per region, from a label image of `run`."""
        signal = fold.tic_image(run)
        total = float(signal[labels >= 0].sum()) or 1.0
        rows, cols = np.indices(labels.shape)
        out = []
        for k in range(len(self)):
            sel = labels == k
            out.append(RegionSummary(
                label=k, n_pixels=int(sel.sum()), percent=100.0 * float(signal[sel].sum()) / total,
                t1_median=float(np.median(fold.t1[cols[sel]] + 0.5 * fold.spec.period)) if sel.any() else float("nan"),
                t2_median=float(np.median(rows[sel] * fold.row_pitch)) if sel.any() else float("nan"),
                top_ions=self.top_ions(k)))
        return out


# --------------------------------------------------------------- internals
def _scan_spectra(run: Run, s: SegmentSettings, mz_use: np.ndarray | None, percentile: float, limit: int | None = None):
    """Weighted, unit-length spectra of the brightest scans (sparse), and which scans they are."""
    from sklearn.preprocessing import normalize
    if mz_use is None:
        keep = (run.mz >= s.mz_min) & ~np.isin(run.mz, s.exclude_ions)
        mz_use = run.mz[keep]
    cols = run.mz_index(mz_use)
    tic = run.tic
    positive = tic[tic > 0]
    if positive.size == 0:
        return sp.csr_matrix((0, mz_use.size)), np.zeros(0, dtype=np.int64), mz_use
    idx = np.flatnonzero(tic >= np.percentile(positive, percentile))
    if limit is not None and idx.size > limit:
        idx = idx[np.argpartition(tic[idx], -limit)[-limit:]]
    idx = np.sort(idx)
    M = run.data[idx][:, cols].astype(np.float64)
    M.data = np.clip(M.data, 0, None) ** s.weight_a
    M = normalize(M @ sp.diags(mz_use.astype(np.float64) ** s.weight_b))
    nonempty = np.asarray(M.getnnz(axis=1)).ravel() > 0
    return M[nonempty], idx[nonempty], mz_use


def _reduce(svd, M) -> np.ndarray:
    Z = svd.transform(M)
    return Z / np.maximum(np.linalg.norm(Z, axis=1, keepdims=True), 1e-12)


def _paint(labels: np.ndarray, scans: np.ndarray, fold: Fold, k: int, smooth: int) -> np.ndarray:
    """Put each kept scan's region on its pixel (the commonest where several share one), then smooth."""
    col, row = _scan_positions(fold)
    pixel = col[scans] * fold.n_rows + row[scans]
    n_pix = fold.n_rows * fold.n_cols
    counts = np.bincount(labels * n_pix + pixel, minlength=k * n_pix).reshape(k, n_pix)
    seen = counts.sum(axis=0) > 0
    grid = np.where(seen, counts.argmax(axis=0), -1).reshape(fold.n_cols, fold.n_rows).T
    if smooth > 0:
        size = 2 * smooth + 1
        votes = np.stack([ndi.uniform_filter((grid == j).astype(np.float32), size=size, mode="constant")
                          for j in range(k)])
        # Never grow into baseline: a pixel that had no scan stays unlabelled.
        grid = np.where(grid >= 0, votes.argmax(axis=0), -1)
    return grid


# ------------------------------------------------------------------- fitting
def fit_regions(pairs, settings: SegmentSettings = SegmentSettings()) -> SegmentResult:
    """Fit the regions on one or more runs. `pairs` is a list of `(run, fold)`; the runs must share an m/z axis
    (the range they all cover is used). The model does not depend on the folds, which are only needed to
    paint a run's regions onto its own pixel grid (`SegmentResult.label_image`)."""
    from sklearn.cluster import KMeans
    from sklearn.decomposition import TruncatedSVD
    pairs = list(pairs)
    if not pairs:
        raise ValueError("no runs to fit")
    lo, hi = max(int(r.mz[0]) for r, _ in pairs), min(int(r.mz[-1]) for r, _ in pairs)
    if hi < lo:
        raise ValueError("the runs share no m/z range")
    axis = np.arange(lo, hi + 1)
    mz_use = axis[(axis >= settings.mz_min) & ~np.isin(axis, settings.exclude_ions)]
    per_run = max(settings.max_scans // len(pairs), 1)
    chosen = [(r, *_scan_spectra(r, settings, mz_use, settings.keep_percentile, per_run)[:2]) for r, _ in pairs]
    M = sp.vstack([m for _, m, _ in chosen if m.shape[0]]).tocsr()
    if M.shape[0] < max(settings.n_regions, 2):
        raise ValueError("too few bright scans to find regions")
    n_comp = int(min(settings.n_components, M.shape[1] - 1, M.shape[0] - 1))
    svd = TruncatedSVD(n_components=n_comp, random_state=settings.seed).fit(M)
    Z = _reduce(svd, M)
    km = KMeans(n_clusters=settings.n_regions, n_init=3, random_state=settings.seed).fit(Z)
    labels = km.labels_
    centres = km.cluster_centers_ / np.maximum(np.linalg.norm(km.cluster_centers_, axis=1, keepdims=True), 1e-12)

    # Number by the signal held, largest first.
    signal = np.concatenate([r.tic[idx] for r, m, idx in chosen if m.shape[0]])
    held = np.bincount(labels, weights=signal, minlength=settings.n_regions)
    order = np.argsort(held)[::-1]
    rank = np.empty(order.size, dtype=int)
    rank[order] = np.arange(order.size)
    labels = rank[labels]
    centres = centres[order]

    # Mean share of ion current per region, over the full m/z axis, for reading a region by its ions.
    from sklearn.preprocessing import normalize
    raw = sp.vstack([r.data[idx][:, r.mz_index(mz_use)].astype(np.float64) for r, m, idx in chosen if m.shape[0]]).tocsr()
    share = normalize(raw, norm="l1")
    indicator = sp.csr_matrix((np.ones(labels.size), (labels, np.arange(labels.size))), shape=(settings.n_regions, labels.size))
    counts = np.maximum(np.asarray(indicator.sum(axis=1)).ravel(), 1)
    mean_used = (indicator @ share).toarray() / counts[:, None]
    mean_spectra = np.zeros((settings.n_regions, axis.size))
    mean_spectra[:, np.searchsorted(axis, mz_use)] = mean_used
    return SegmentResult(settings=settings, mz=axis, mz_used=mz_use, svd=svd, centroids=centres, mean_spectra=mean_spectra)


def segment_pixels(run: Run, fold: Fold, settings: SegmentSettings = SegmentSettings()):
    """Fit regions on one run and label it: `(SegmentResult, label image)`."""
    result = fit_regions([(run, fold)], settings)
    return result, result.label_image(run, fold)


def region_shares(labels: np.ndarray, run: Run, fold: Fold, n_regions: int) -> np.ndarray:
    """Percent of `run`'s signal, over its labelled pixels, that falls in each region (length `n_regions`)."""
    signal = fold.tic_image(run)
    inside = labels >= 0
    total = float(signal[inside].sum()) or 1.0
    return 100.0 * np.bincount(labels[inside], weights=signal[inside], minlength=n_regions) / total
