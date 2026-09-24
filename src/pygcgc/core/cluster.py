"""Unsupervised grouping of peaks, within one chromatogram and across many.

The idea in plain terms
------------------------
`chemclass.py` says which kind of compound a peak is, from rules you write
down in advance. This asks a different question with no rules: which peaks
look alike? Peaks are grouped by the shape of their mass spectra alone, so
an alkane series, an alkylbenzene band or a family nobody has written a
rule for each fall out as a cluster. Plotted at their retention times, the
clusters divide the chromatogram into chemically coherent regions.

By default retention time is deliberately not part of the distance. Two
runs never line up perfectly, and a clustering that used position would turn
drift into false differences. Spectra do not drift, so clusters found this way
mean the same thing in every run and no alignment is needed to compare runs.
`ClusterSettings(method="density")` is the exception: it adds retention and
clusters by HDBSCAN (scikit-learn), to find bands and series that spectra
alone cannot tell apart, at the price of leaving more peaks unassigned and
needing runs aligned before they are pooled. `segment.py` goes the other way
and divides the chromatogram into regions with no peaks at all.

Method (`cluster_peaks`)
-------------------------
1. Each peak's spectrum is cut to the usable m/z range, weighted
   (intensity^a * m/z^b) and scaled to unit length. Peaks with too few ions
   to describe a shape are left unassigned rather than forced into a group.
2. The largest `max_fit` peaks are clustered by hierarchical clustering
   (average linkage on cosine distance) and the tree is cut at
   `distance_threshold`, or at `n_clusters`. Big peaks have the cleanest
   spectra, and the cost of the tree grows with the square of the number of
   peaks, so the rest are not put in it.
3. Every other peak joins the cluster with the nearest centroid, if that
   is within `assign_max_distance`, and is otherwise left unassigned (-1).
   The same rule (`ClusterResult.assign`) places a new run into an existing
   segmentation.
4. Clusters are numbered by total volume, largest first.

Comparing runs
---------------
Pool the peaks of several runs into one `PeakTable`, cluster once, and each
run is described by how much volume it has in each cluster
(`cluster_features`). That is a small, comparable table, and it is returned
as the same `Feature` list the batch pipeline produces, so `pca`, `volcano`
and `plot_pca` work on it unchanged. `permanova` tests whether groups of
runs differ overall.

Peaks come from `PeakTable.from_blobs` (the legacy detector). Any other
detector only needs its own constructor.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, replace
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform

from .batch import Feature
from .peaks import BlobSet

LINKAGES = ("average", "complete", "weighted", "single")
METHODS = ("tree", "density")

# Okabe and Ito (2008) without black, which disappears on a dark background. A chromatogram can have a
# hundred clusters, more than colours can tell apart, so only the first few get a colour and the rest share one.
CLUSTER_COLOURS = ((230, 159, 0), (86, 180, 233), (0, 158, 115), (240, 228, 66),
                   (0, 114, 178), (213, 94, 0), (204, 121, 167))
OTHER_COLOUR = (190, 190, 190)


def cluster_colour(k: int) -> tuple:
    return CLUSTER_COLOURS[k] if 0 <= k < len(CLUSTER_COLOURS) else OTHER_COLOUR


@dataclass(frozen=True)
class ClusterSettings:
    mz_min: int = 35                  # below this, air, water, CO2 and solvent dominate the spectrum
    exclude_ions: tuple = ()          # e.g. SILOXANE_IONS, so column bleed does not group by itself
    weight_a: float = 0.5             # intensity exponent; gentler than the library-search 0.6 so weak ions still count
    weight_b: float = 1.0             # m/z exponent; library search uses 3, which would cluster on the high-mass end alone
    min_ions: int = 5                 # a spectrum with fewer non-zero ions is left unassigned (as in ClassSettings)
    linkage: str = "average"
    n_clusters: int | None = None     # if set, cut the tree here instead of at distance_threshold
    distance_threshold: float = 0.3   # cosine distance; 0.3 is a mean similarity of 0.7 within a merged group
    min_cluster_size: int = 3         # smaller groups are left unassigned
    max_fit: int = 5000               # peaks in the tree; 5000 needs about 200 MB
    assign_max_distance: float = 0.4  # a peak further than this from every centroid is left unassigned
    soft_weight: float = 0.0          # weight of the 14 eV spectrum against the 70 eV one (0 = ignore it; 0.5 = half as much)
    method: str = "tree"              # "tree": spectra only, hierarchical. "density": spectra and retention, HDBSCAN
    retention_weight: float = 1.0     # density only: how much retention counts against spectrum (0 = not at all)
    t1_scale: float = 120.0           # density only: seconds of t1 that count as one unit of spectral distance
    t2_scale: float = 1.0             # density only: seconds of t2 that count as one unit of spectral distance
    density_min_size: int = 15        # density only: smallest group of peaks called a cluster
    density_min_samples: int = 3      # density only: neighbours a peak needs to count as dense; lower gives more, smaller clusters
    density_selection: str = "leaf"   # density only: "leaf" (finest stable clusters) or "eom" (the largest; on a real
                                      # run it chains almost every peak into one cluster because density is continuous)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "ClusterSettings":
        d = {k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__}
        if "exclude_ions" in d:
            d["exclude_ions"] = tuple(d["exclude_ions"])
        return cls(**d)


# ------------------------------------------------------------------- peaks
@dataclass
class PeakTable:
    """What clustering needs to know about a set of peaks, whichever run or
    detector they came from."""
    run: np.ndarray          # (n,) name of the run each peak is from
    blob: np.ndarray         # (n,) its index within that run's detector output
    t1: np.ndarray           # s
    t2: np.ndarray           # s
    volume: np.ndarray
    co_elution: np.ndarray   # (n,) bool: spectrum is probably a mixture
    spectra: np.ndarray      # (n, len(mz))
    mz: np.ndarray           # contiguous unit masses
    soft_spectra: np.ndarray | None = None   # (n, len(soft_mz)): the same peaks at 14 eV, if a partner run was loaded
    soft_mz: np.ndarray | None = None

    def __len__(self) -> int:
        return int(self.t1.size)

    @classmethod
    def from_blobs(cls, blobs: BlobSet, run: str = "", prefer_apex: bool = True, min_ions: int = 5,
                   soft: tuple | None = None) -> "PeakTable":
        """Blobs from either detector. The apex spectrum is used where it has
        enough ions, because a footprint sum picks up a co-eluting neighbour
        (the same choice `classify_blobs` makes). `soft` is `(spectra, mz)`
        from `soft_spectra_for_blobs`: the same blobs at a lower ionisation energy."""
        n = len(blobs)
        if n == 0:
            return cls(np.array([], dtype=object), np.array([], dtype=int), np.array([]), np.array([]), np.array([]),
                       np.array([], dtype=bool), np.zeros((0, blobs.mz.size), dtype=np.float32), blobs.mz)
        spectra = blobs.spectra
        if prefer_apex and blobs.apex_spectra.shape[0] == n:
            thin = (blobs.apex_spectra > 0).sum(axis=1) < min_ions
            spectra = np.where(thin[:, None], blobs.spectra, blobs.apex_spectra)
        co = blobs.co_elution if blobs.co_elution.size == n else np.zeros(n, dtype=bool)
        soft_spectra, soft_mz = (None, None) if soft is None else (np.asarray(soft[0], dtype=np.float32), np.asarray(soft[1]))
        if soft_spectra is not None and soft_spectra.shape[0] != n:
            raise ValueError(f"{soft_spectra.shape[0]} soft spectra for {n} blobs")
        return cls(run=np.full(n, run, dtype=object), blob=np.arange(n), t1=np.asarray(blobs.t1, float),
                   t2=np.asarray(blobs.t2, float), volume=np.asarray(blobs.volume, float),
                   co_elution=np.asarray(co, dtype=bool), spectra=np.asarray(spectra, dtype=np.float32),
                   mz=np.asarray(blobs.mz), soft_spectra=soft_spectra, soft_mz=soft_mz)

    @classmethod
    def concat(cls, tables: list) -> "PeakTable":
        """Pool several tables. Runs can differ in m/z range; the pool keeps
        the range they all cover. The 14 eV spectra are kept only if every
        run has them (a run without would otherwise look like a run of
        featureless soft spectra)."""
        tables = [t for t in tables if len(t)]
        if not tables:
            raise ValueError("no peaks to pool")

        def trimmed(get_spectra, get_mz):
            lo = max(int(get_mz(t)[0]) for t in tables)
            hi = min(int(get_mz(t)[-1]) for t in tables)
            if hi < lo:
                return None, None
            return (np.vstack([get_spectra(t)[:, lo - int(get_mz(t)[0]):hi - int(get_mz(t)[0]) + 1] for t in tables]),
                    np.arange(lo, hi + 1))

        spectra, mz = trimmed(lambda t: t.spectra, lambda t: t.mz)
        if spectra is None:
            raise ValueError("the runs share no m/z range")
        soft_spectra = soft_mz = None
        if all(t.soft_spectra is not None for t in tables):
            soft_spectra, soft_mz = trimmed(lambda t: t.soft_spectra, lambda t: t.soft_mz)
        cat = lambda name: np.concatenate([getattr(t, name) for t in tables])  # noqa: E731
        return cls(run=cat("run"), blob=cat("blob"), t1=cat("t1"), t2=cat("t2"), volume=cat("volume"),
                   co_elution=cat("co_elution"), spectra=spectra, mz=mz, soft_spectra=soft_spectra, soft_mz=soft_mz)

    def subset(self, index) -> "PeakTable":
        return PeakTable(self.run[index], self.blob[index], self.t1[index], self.t2[index], self.volume[index],
                         self.co_elution[index], self.spectra[index], self.mz,
                         None if self.soft_spectra is None else self.soft_spectra[index], self.soft_mz)


# ------------------------------------------------------------------ result
@dataclass
class ClusterSummary:
    label: int
    n_peaks: int
    n_runs: int
    volume: float
    percent: float            # of the volume of every peak in the table, unassigned included
    t1_median: float          # s
    t2_median: float          # s
    top_ions: list            # m/z with the largest mean share of the ion current
    silhouette: float         # mean over this cluster's peaks in the tree; nan if it had none


@dataclass
class ClusterResult:
    labels: np.ndarray        # (n,) cluster of each peak, -1 for unassigned
    settings: ClusterSettings
    mz: np.ndarray            # m/z axis of the table that was clustered
    mz_used: np.ndarray       # the ions the distance is built from
    centroids: np.ndarray     # (k, len(mz_used)) unit-length weighted spectra
    mean_spectra: np.ndarray  # (k, len(mz)) mean share of ion current, for reading a cluster
    fit_index: np.ndarray     # rows of the table that went into the tree
    tree: np.ndarray          # scipy linkage matrix of those rows
    silhouette: np.ndarray    # (len(fit_index),) nan where unassigned
    soft_used: np.ndarray | None = None   # the 14 eV ions in the distance, if it used them (centroids then hold both blocks)

    def __len__(self) -> int:
        return int(self.centroids.shape[0])

    def top_ions(self, k: int, n: int = 6) -> list:
        return [int(m) for m in self.mz[np.argsort(self.mean_spectra[k])[::-1][:n]]]

    def assign(self, table: PeakTable) -> np.ndarray:
        """Cluster of every peak in `table`, e.g. a run that was not in the fit."""
        X, ok, _, _ = _features(table, self.settings, self.mz_used, self.soft_used)
        return _nearest_centroid(X, ok, self.centroids, self.settings.assign_max_distance)

    def summary(self, table: PeakTable) -> list:
        """One `ClusterSummary` per cluster, largest first. `table` is the one that was clustered."""
        total = float(table.volume.sum()) or 1.0
        fit_labels = self.labels[self.fit_index]
        rows = []
        for k in range(len(self)):
            sel = self.labels == k
            s = self.silhouette[fit_labels == k]
            s = s[~np.isnan(s)]              # peaks placed by centroid rather than by the tree have none
            rows.append(ClusterSummary(
                label=k, n_peaks=int(sel.sum()), n_runs=len(set(table.run[sel])), volume=float(table.volume[sel].sum()),
                percent=100.0 * float(table.volume[sel].sum()) / total, t1_median=float(np.median(table.t1[sel])),
                t2_median=float(np.median(table.t2[sel])), top_ions=self.top_ions(k),
                silhouette=float(s.mean()) if s.size else float("nan")))
        return rows


# --------------------------------------------------------------- internals
def _block(spectra: np.ndarray, mz: np.ndarray, mz_use: np.ndarray, s: ClusterSettings):
    """Weighted spectra on the ions in `mz_use`, their lengths, and the unweighted intensities."""
    cols = mz_use - int(mz[0])
    if cols.min() < 0 or cols.max() >= mz.size:
        raise ValueError(f"the table covers m/z {mz[0]}-{mz[-1]}, the clusters were built on {mz_use[0]}-{mz_use[-1]}")
    raw = np.clip(spectra[:, cols], 0, None)
    X = raw ** s.weight_a * mz_use.astype(np.float32) ** s.weight_b
    return X, np.linalg.norm(X, axis=1), raw


def _features(table: PeakTable, s: ClusterSettings, mz_use: np.ndarray | None = None,
              soft_use: np.ndarray | None = None):
    """Weighted, unit-length spectra of every peak, which rows had enough ions
    to be used, and the ions used (given, or chosen from the settings). With
    `soft_weight` above zero and 14 eV spectra in the table, each spectrum
    is the 70 eV block followed by the 14 eV block scaled by the square root
    of that weight, so the cosine between two peaks is the weighted mean of
    their two cosines. `mz_use` and `soft_use` from an earlier result make a
    new table comparable with it."""
    keep = lambda axis: axis[(axis >= s.mz_min) & ~np.isin(axis, s.exclude_ions)]   # noqa: E731
    if mz_use is None:
        mz_use = keep(table.mz)
    X, norm, raw = _block(table.spectra, table.mz, mz_use, s)
    ok = ((raw > 0).sum(axis=1) >= s.min_ions) & (norm > 0)
    X[ok] /= norm[ok, None]
    X[~ok] = 0.0
    if soft_use is not None and table.soft_spectra is None:
        raise ValueError("these clusters were built with 14 eV spectra; the peaks have none")
    if soft_use is None and s.soft_weight > 0 and table.soft_spectra is not None:
        soft_use = keep(table.soft_mz)
    if soft_use is not None and soft_use.size:
        Xs, ns, _ = _block(table.soft_spectra, table.soft_mz, soft_use, s)
        has = ns > 0                       # a peak with no 14 eV signal just contributes nothing from that block
        Xs[has] /= ns[has, None]
        Xs[~has] = 0.0
        X = np.hstack([X, np.sqrt(s.soft_weight) * Xs])
        X /= np.maximum(np.linalg.norm(X, axis=1, keepdims=True), 1e-12)
    else:
        soft_use = None
    return X.astype(np.float32), ok, mz_use, soft_use


def _nearest_centroid(X: np.ndarray, ok: np.ndarray, centroids: np.ndarray, max_distance: float) -> np.ndarray:
    labels = np.full(X.shape[0], -1, dtype=int)
    if centroids.shape[0] == 0 or not ok.any():
        return labels
    dist = 1.0 - X[ok] @ centroids.T
    best = dist.argmin(axis=1)
    labels[ok] = np.where(dist[np.arange(best.size), best] <= max_distance, best, -1)
    return labels


def _silhouette(dist: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Silhouette of every clustered row of a square distance matrix
    (Rousseeuw, 1987); nan for unassigned rows, and everywhere if there is
    only one cluster, where it is undefined."""
    out = np.full(labels.size, np.nan)
    inside = labels >= 0
    k = int(labels.max()) + 1 if inside.any() else 0
    if k < 2:
        return out
    onehot = np.zeros((int(inside.sum()), k))
    onehot[np.arange(onehot.shape[0]), labels[inside]] = 1.0
    size = onehot.sum(axis=0)
    total = dist[np.ix_(inside, inside)] @ onehot                 # summed distance from each peak to each cluster
    own = np.take_along_axis(total, labels[inside][:, None], 1)[:, 0] / np.maximum(size[labels[inside]] - 1, 1)
    others = total / size
    others[np.arange(onehot.shape[0]), labels[inside]] = np.inf
    nearest = others.min(axis=1)
    s = (nearest - own) / np.maximum(np.maximum(own, nearest), 1e-12)
    s[size[labels[inside]] == 1] = 0.0                            # the convention for a cluster of one
    out[inside] = s
    return out


# ------------------------------------------------------------------ cluster
def cluster_peaks(table: PeakTable, settings: ClusterSettings = ClusterSettings()) -> ClusterResult:
    if settings.linkage not in LINKAGES:
        raise ValueError(f"linkage must be one of {LINKAGES}, not {settings.linkage!r}")
    if settings.method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, not {settings.method!r}")
    X, ok, mz_used, soft_used = _features(table, settings)
    if settings.method == "density":
        return _cluster_density(table, settings, X, ok, mz_used, soft_used)

    usable = np.flatnonzero(ok)
    fit_index = np.sort(usable[np.argsort(table.volume[usable])[::-1][:settings.max_fit]])
    labels = np.full(len(table), -1, dtype=int)
    centroids = np.zeros((0, X.shape[1]), dtype=np.float32)
    tree = np.zeros((0, 4))
    silhouette = np.full(fit_index.size, np.nan)

    if fit_index.size >= 2:
        Xf = X[fit_index]
        dist = np.clip(1.0 - Xf @ Xf.T, 0.0, None).astype(np.float64)   # unit-length, non-negative spectra: 0 to 1
        np.fill_diagonal(dist, 0.0)
        tree = linkage(squareform(dist, checks=False), method=settings.linkage)
        cut = (fcluster(tree, settings.n_clusters, criterion="maxclust") if settings.n_clusters
               else fcluster(tree, settings.distance_threshold, criterion="distance"))
        big = [c for c in np.unique(cut) if (cut == c).sum() >= settings.min_cluster_size]
        fit_labels = np.full(fit_index.size, -1, dtype=int)
        for new, c in enumerate(big):
            fit_labels[cut == c] = new
        silhouette = _silhouette(dist, fit_labels)

        if big:
            centroids = np.array([Xf[fit_labels == k].mean(axis=0) for k in range(len(big))])
            centroids /= np.linalg.norm(centroids, axis=1, keepdims=True)
            labels = _nearest_centroid(X, ok, centroids, settings.assign_max_distance)
            # The tree's decision stands for peaks it put in a real cluster. Peaks it left in tiny groups
            # (typically small, noisy ones) get the same nearest-centroid chance as peaks outside the tree.
            in_cluster = fit_labels >= 0
            labels[fit_index[in_cluster]] = fit_labels[in_cluster]

            # Number by volume, largest first.
            vol = np.array([table.volume[labels == k].sum() for k in range(len(big))])
            rank = np.empty(len(big), dtype=int)
            rank[np.argsort(vol)[::-1]] = np.arange(len(big))
            labels = np.where(labels >= 0, rank[np.clip(labels, 0, None)], -1)
            centroids = centroids[np.argsort(vol)[::-1]]

    return ClusterResult(labels=labels, settings=settings, mz=table.mz, mz_used=mz_used, centroids=centroids,
                         mean_spectra=_mean_spectra(table, labels, centroids.shape[0]), fit_index=fit_index,
                         tree=tree, silhouette=silhouette, soft_used=soft_used)


def _cluster_density(table: PeakTable, s: ClusterSettings, X: np.ndarray, ok: np.ndarray, mz_used, soft_used) -> ClusterResult:
    """Spectra and retention together, by HDBSCAN (Campello et al., 2013) as in scikit-learn. The unit-length
    spectra are reduced to their ten main directions (the Euclidean distance between unit vectors is
    sqrt(2 (1 - cosine)), so the ordering is the cosine distance's) and retention time, scaled by `t1_scale`,
    `t2_scale` and `retention_weight`, is appended. HDBSCAN needs no cluster count, copes with clusters of
    different density, and leaves what belongs to none as noise (-1). Retention is compared as recorded, so
    pooling runs needs them aligned first. The centroids are spectral only: `ClusterResult.assign` places a new
    peak by its spectrum alone. There is no tree and no silhouette."""
    try:
        from sklearn.cluster import HDBSCAN
    except ImportError as exc:                       # pragma: no cover - depends on the environment
        raise ImportError("retention-aware clustering needs scikit-learn: pip install 'pygcgc[cluster]'") from exc
    usable = np.flatnonzero(ok)
    labels = np.full(len(table), -1, dtype=int)
    centroids = np.zeros((0, X.shape[1]), dtype=np.float32)
    if usable.size >= max(s.density_min_size, 2) + 1:
        centred = X[usable] - X[usable].mean(axis=0)
        k = min(10, centred.shape[1], usable.size - 1)
        _, _, vt = np.linalg.svd(centred, full_matrices=False)
        embedding = centred @ vt[:k].T
        retention = np.column_stack([table.t1[usable] / s.t1_scale, table.t2[usable] / s.t2_scale]) * s.retention_weight
        options = dict(min_cluster_size=max(s.density_min_size, 2), min_samples=s.density_min_samples,
                       cluster_selection_method=s.density_selection)
        try:
            model = HDBSCAN(copy=True, **options)          # copy: its default changes in scikit-learn 1.10
        except TypeError:                                  # older scikit-learn has no such argument
            model = HDBSCAN(**options)
        raw = model.fit(np.hstack([embedding, retention])).labels_
        found = [c for c in np.unique(raw) if c >= 0]
        if found:
            vol = np.array([table.volume[usable][raw == c].sum() for c in found])
            order = np.argsort(vol)[::-1]                                    # number by volume, largest first
            for new, c in enumerate(np.array(found)[order]):
                labels[usable[raw == c]] = new
            centroids = np.array([X[labels == j].mean(axis=0) for j in range(len(found))])
            centroids /= np.maximum(np.linalg.norm(centroids, axis=1, keepdims=True), 1e-12)
    return ClusterResult(labels=labels, settings=s, mz=table.mz, mz_used=mz_used, centroids=centroids.astype(np.float32),
                         mean_spectra=_mean_spectra(table, labels, centroids.shape[0]), fit_index=usable,
                         tree=np.zeros((0, 4)), silhouette=np.full(usable.size, np.nan), soft_used=soft_used)


def _mean_spectra(table: PeakTable, labels: np.ndarray, k: int) -> np.ndarray:
    """Share of ion current, so a cluster is read as "these ions, in these proportions"."""
    share = table.spectra / np.maximum(table.spectra.clip(0).sum(axis=1, keepdims=True), 1e-12)
    return np.array([share[labels == j].mean(axis=0) for j in range(k)]).reshape(k, table.mz.size)


# ------------------------------------------------------------- across runs
def apply_pqn(features: list, run_names: list) -> list:
    """Probabilistic quotient normalisation (Dieterle et al., 2006) of a run-by-feature table, returned as new
    `Feature`s. Percent-of-total scaling assumes the total is comparable between runs; when one big cluster
    changes, it is not, and every other cluster appears to change the other way. PQN instead takes the median
    profile over runs as a reference, finds each run's median quotient against it (the dilution-like factor
    most features agree on, so a few changed features do not move it) and divides the run by that.
    Apply it to percent-scaled values."""
    matrix = np.array([[f.volumes.get(r, 0.0) for f in features] for r in run_names], dtype=float)
    reference = np.median(matrix, axis=0)
    use = reference > 0
    quotient = np.median(matrix[:, use] / reference[use], axis=1) if use.any() else np.ones(len(run_names))
    quotient = np.where(quotient > 0, quotient, 1.0)
    return [Feature(id=f.id, ref_t1=f.ref_t1, ref_t2=f.ref_t2,
                    volumes={r: float(matrix[i, j] / quotient[i]) for i, r in enumerate(run_names)})
            for j, f in enumerate(features)]


def cluster_features(table: PeakTable, result: ClusterResult, run_names: list | None = None,
                     relative: bool = False, pqn: bool = False) -> list:
    """Each cluster as a `Feature` whose volumes are the run totals in that
    cluster, in the form `pca`, `volcano` and `feature_table_matrix` take.
    A run with nothing in a cluster gets 0.0. `run_names` lists the runs to
    include (default: those in the table, in order of appearance); name a
    run that found no peaks at all and it appears as zeros. `relative`
    gives percent of the run's total peak volume instead, which removes
    differences in overall loading; normalise to an internal standard first
    if you want absolute amounts. `pqn` then applies `apply_pqn` (use it
    with `relative`)."""
    if run_names is None:
        run_names = list(dict.fromkeys(table.run))
    totals = {r: float(table.volume[table.run == r].sum()) for r in run_names}
    features = []
    for k in range(len(result)):
        in_k = result.labels == k
        vols = {}
        for r in run_names:
            v = float(table.volume[in_k & (table.run == r)].sum())
            vols[r] = 100.0 * v / totals[r] if relative and totals[r] > 0 else v
        features.append(Feature(id=k, ref_t1=float(np.median(table.t1[in_k])), ref_t2=float(np.median(table.t2[in_k])),
                                volumes=vols))
    return apply_pqn(features, run_names) if pqn and features else features


def permanova(X: np.ndarray, groups, metric: str = "braycurtis", permutations: int = 999, seed: int = 0) -> tuple:
    """Permutational MANOVA (Anderson, 2001): do the groups of rows of `X`
    (runs, by cluster volume) sit further apart than chance allows? Returns
    the pseudo-F and a permutation p-value. Bray-Curtis suits abundances
    such as `cluster_features(..., relative=True)`; with metric="euclidean"
    on one column the pseudo-F is the ordinary ANOVA F."""
    groups = np.asarray(groups)
    _, g = np.unique(groups, return_inverse=True)
    n, a = g.size, int(g.max()) + 1
    if n <= a or a < 2:
        raise ValueError("need at least two groups and more runs than groups")
    d2 = squareform(pdist(np.asarray(X, dtype=float), metric=metric)) ** 2
    ss_total = d2.sum() / (2 * n)

    def pseudo_f(assignment):
        onehot = np.zeros((n, a))
        onehot[np.arange(n), assignment] = 1.0
        ss_within = (np.diag(onehot.T @ d2 @ onehot) / (2 * onehot.sum(axis=0))).sum()
        if ss_within <= 1e-12 * max(ss_total, 1e-300):       # identical runs within every group: perfectly separated
            return np.inf
        return ((ss_total - ss_within) / (a - 1)) / (ss_within / (n - a))

    observed = pseudo_f(g)
    rng = np.random.default_rng(seed)
    beaten = sum(pseudo_f(rng.permutation(g)) >= observed for _ in range(permutations))
    return float(observed), (1 + beaten) / (1 + permutations)


def write_clusters_csv(table: PeakTable, result: ClusterResult, path, names: dict | None = None) -> Path:
    """One row per peak: where it is, how big, and which cluster. `names` maps a cluster to a name you gave it;
    if any are given there is a `cluster_name` column."""
    path = Path(path)
    names = names or {}
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["run", "blob", "t1_min", "t2_s", "volume", "cluster", "cluster_top_ions"] + (["cluster_name"] if names else []))
        for i in range(len(table)):
            k = int(result.labels[i])
            w.writerow([table.run[i], int(table.blob[i]), f"{table.t1[i] / 60:.3f}", f"{table.t2[i]:.3f}",
                        f"{table.volume[i]:.6g}", k if k >= 0 else "", " ".join(map(str, result.top_ions(k, 5))) if k >= 0 else ""]
                       + ([names.get(k, "") if k >= 0 else ""] if names else []))
    return path


# ------------------------------------------------------- reading the result
def clusters_of_run(table: PeakTable, result: ClusterResult, run: str, n_blobs: int) -> np.ndarray:
    """Cluster of each of `run`'s `n_blobs` detector blobs, -1 where a blob has none."""
    out = np.full(n_blobs, -1, dtype=int)
    mine = table.run == run
    out[table.blob[mine].astype(int)] = result.labels[mine]
    return out


def cluster_label_image(blobs: BlobSet, cluster_of_blob) -> np.ndarray:
    """The cluster of every pixel of the plain folded image, as (rows,
    modulations): the cluster of the blob whose footprint covers the pixel,
    or -1 outside every clustered blob. The detector works on an image
    extended above the top edge (the bottom of the next modulation, so
    wrapped peaks are whole); this puts that band back where it came from,
    so the result lines up with `Fold.tic_image`."""
    lab = np.concatenate([[-1], np.asarray(cluster_of_blob, dtype=int)])   # label 0 is "no blob"
    ext, w = blobs.labels, blobs.wrap_rows
    n_rows = ext.shape[0] - w
    plain = np.zeros((n_rows, ext.shape[1]), dtype=ext.dtype)
    plain[w:] = ext[w:n_rows]
    if w:
        plain[:w, 1:] = ext[n_rows:, :-1]     # extension row n_rows + j of column c is plain row j of column c + 1
    return lab[plain]


def class_composition(labels, classes) -> dict:
    """For each cluster, the classes (from `classify_blobs`) of its peaks as
    [(class, fraction of the cluster's peaks)], largest first. `classes` is one
    name per peak; None or "" counts as "unclassified". Compares the two
    approaches: a cluster that is 95% one class is that class found without rules."""
    labels = np.asarray(labels, dtype=int)
    names = np.array([c or "unclassified" for c in classes], dtype=object)
    out = {}
    for k in np.unique(labels[labels >= 0]):
        found, counts = np.unique(names[labels == k], return_counts=True)
        out[int(k)] = sorted(((str(n), float(c) / counts.sum()) for n, c in zip(found, counts)), key=lambda x: -x[1])
    return out


def merge_clusters(table: PeakTable, result: ClusterResult, merge) -> ClusterResult:
    """A copy of `result` with the clusters in `merge` joined into one and
    everything renumbered by volume. The distances are not recomputed: the
    joined cluster's centroid is the size-weighted mean of its parts, and its
    peaks lose their silhouette (which described the old split)."""
    merge = sorted({int(k) for k in merge})
    if len(merge) < 2 or merge[0] < 0 or merge[-1] >= len(result):
        raise ValueError("name at least two existing clusters to merge")
    target = np.arange(len(result))
    target[merge] = merge[0]
    keep = sorted(set(target.tolist()))
    new_of = np.array([keep.index(t) for t in target])
    labels = np.where(result.labels >= 0, new_of[np.clip(result.labels, 0, None)], -1)
    counts = np.array([(result.labels == j).sum() for j in range(len(result))])
    centroids = []
    for new in range(len(keep)):
        parts = np.flatnonzero(new_of == new)
        c = (result.centroids[parts] * counts[parts, None]).sum(axis=0)
        centroids.append(c / max(float(np.linalg.norm(c)), 1e-12))
    centroids = np.array(centroids, dtype=np.float32)
    order = np.argsort([-table.volume[labels == j].sum() for j in range(len(keep))])
    rank = np.empty(len(keep), dtype=int)
    rank[order] = np.arange(len(keep))
    labels = np.where(labels >= 0, rank[np.clip(labels, 0, None)], -1)
    silhouette = np.where(np.isin(result.labels[result.fit_index], merge), np.nan, result.silhouette)
    return replace(result, labels=labels, centroids=centroids[order], silhouette=silhouette,
                   mean_spectra=_mean_spectra(table, labels, len(keep)))


# ---------------------------------------------------- how much to trust it
def adjusted_rand(a, b) -> float:
    """Adjusted Rand index (Hubert and Arabie, 1985) between two labellings,
    over the peaks both assign: 1 is identical partitions, about 0 is chance."""
    a, b = np.asarray(a, dtype=int), np.asarray(b, dtype=int)
    both = (a >= 0) & (b >= 0)
    a, b = np.unique(a[both], return_inverse=True)[1], np.unique(b[both], return_inverse=True)[1]
    n = a.size
    if n < 2:
        return float("nan")
    table = np.zeros((a.max() + 1, b.max() + 1))
    np.add.at(table, (a, b), 1)
    pairs = lambda x: x * (x - 1) / 2   # noqa: E731
    within, rows, cols = pairs(table).sum(), pairs(table.sum(axis=1)).sum(), pairs(table.sum(axis=0)).sum()
    expected = rows * cols / pairs(n)
    denom = 0.5 * (rows + cols) - expected
    return 1.0 if denom == 0 else float((within - expected) / denom)


def threshold_scan(table: PeakTable, settings: ClusterSettings, thresholds) -> list:
    """Cluster at each merge distance in `thresholds` to see how the answer
    changes: clusters, share of the volume left unassigned, mean silhouette."""
    rows = []
    for t in thresholds:
        r = cluster_peaks(table, replace(settings, distance_threshold=float(t), n_clusters=None))
        s = r.silhouette[~np.isnan(r.silhouette)]
        rows.append(dict(threshold=float(t), n_clusters=len(r),
                         unassigned_percent=100.0 * float(table.volume[r.labels < 0].sum()) / (float(table.volume.sum()) or 1.0),
                         silhouette=float(s.mean()) if s.size else float("nan")))
    return rows


def stability(table: PeakTable, settings: ClusterSettings, n_boot: int = 20, fraction: float = 0.8,
              seed: int = 0) -> dict:
    """Would a different sample of the same peaks give the same clusters?
    Clusters a random `fraction` of the usable peaks `n_boot` times, places
    every peak by the nearest centroid, and compares each labelling with the
    full one by the adjusted Rand index. Returns the mean, SD and every value.
    Near 1 means the grouping does not depend on which peaks were seen."""
    reference = cluster_peaks(table, settings)
    if len(reference) < 2:
        raise ValueError("fewer than two clusters: nothing to compare")
    ref_labels = reference.assign(table)
    _, ok, _, _ = _features(table, settings, reference.mz_used, reference.soft_used)
    usable = np.flatnonzero(ok)
    rng = np.random.default_rng(seed)
    scores = []
    for _ in range(n_boot):
        pick = np.sort(rng.choice(usable, size=max(int(fraction * usable.size), 2), replace=False))
        scores.append(adjusted_rand(ref_labels, cluster_peaks(table.subset(pick), settings).assign(table)))
    return dict(mean=float(np.nanmean(scores)), sd=float(np.nanstd(scores)), values=scores, n_boot=n_boot)
