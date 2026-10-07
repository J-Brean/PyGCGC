"""MS-aware chromatographic peak candidates without quantitative borders.

The legacy blob detector reduces every ion to one summed surface before it
looks for peaks.  On a busy real chromatogram, small positive residuals from
hundreds of ions accumulate and make nearly the whole image look significant.
This module keeps the ion evidence separate: it finds a local chromatographic
maximum on each ion, then groups ions whose maxima agree in position.

A candidate is an apex hypothesis, not an integrated compound.  Several
components may later share one analysis region and one pixel may contribute to
several fitted components.  This distinction is important for genuine
co-elution and is why this module deliberately does not return hard borders.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
import os

import numpy as np
from scipy import ndimage as ndi
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .background import storage_threshold
from .fold import Fold
from .peaks import Frame, make_frame, PeakSettings
from .run import Run


def _gaussian_filter_ions(cube: np.ndarray, sigma: tuple[float, float, float]) -> np.ndarray:
    """Filter independent ion planes concurrently without changing results."""
    n_ions = cube.shape[2]
    workers = min(4, max(1, os.cpu_count() or 1), max(n_ions, 1))
    if workers == 1 or n_ions < 8:
        return ndi.gaussian_filter(cube, sigma, mode="nearest")
    chunks = [chunk for chunk in np.array_split(np.arange(n_ions), workers) if chunk.size]

    def filtered(indices):
        return indices, ndi.gaussian_filter(
            cube[:, :, indices], sigma, mode="nearest")

    out = np.empty_like(cube)
    with ThreadPoolExecutor(max_workers=workers) as executor:
        for indices, values in executor.map(filtered, chunks):
            out[:, :, indices] = values
    return out


@dataclass(frozen=True)
class CandidateSettings:
    seed_snr: float = 6.0
    single_ion_snr: float = 15.0
    min_support_ions: int = 2
    tic_prominence_snr: float = 2.0
    smooth_t2_s: float = 0.02
    smooth_modulations: float = 0.35
    broad_t2_s: float = 0.18
    broad_modulations: float = 2.0
    apex_rows: int = 2
    apex_modulations: int = 1
    min_ion_points: int = 15
    noise_block: int = 20
    wrap_below_s: float = 0.5
    poisson_noise: bool = True

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "CandidateSettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


@dataclass
class CandidateSet:
    """Local apex hypotheses and the ions independently supporting each one."""
    t1: np.ndarray
    t2: np.ndarray
    row: np.ndarray
    col: np.ndarray
    score: np.ndarray
    strongest_snr: np.ndarray
    support_count: np.ndarray
    support_ions: list[np.ndarray] = field(default_factory=list)
    support_snr: list[np.ndarray] = field(default_factory=list)
    pitch: float = 0.0
    period: float = 0.0
    t1_start: float = 0.0
    wrap_rows: int = 0

    def __len__(self) -> int:
        return int(self.t1.size)


@dataclass(frozen=True)
class CandidateRegion:
    """A padded local rectangle which may contain several components."""
    candidate_indices: tuple[int, ...]
    row0: int
    row1: int                 # exclusive
    col0: int
    col1: int                 # exclusive


def _block_location_scale(cube: np.ndarray, frame: Frame, block: int,
                          floor: float) -> tuple[np.ndarray, np.ndarray]:
    """Robust local centre and spread for every column and ion.

    Only the upper tail is clipped during refinement, since chromatographic
    peaks are positive.  The raw-file storage floor remains the lower bound on
    scale when cleaned data contain many tiny residuals or exact zeros.
    """
    n_cols, n_ions = cube.shape[1:]
    centre = np.empty((n_cols, n_ions), dtype=np.float32)
    scale = np.empty((n_cols, n_ions), dtype=np.float32)
    complete = frame.fold.complete
    for b0 in range(0, n_cols, max(int(block), 1)):
        cols = np.arange(b0, min(b0 + max(int(block), 1), n_cols))
        use = cols[complete[cols]] if complete[cols].any() else cols
        values = cube[:frame.n_rows, use].reshape(-1, n_ions).astype(np.float64)
        loc = np.median(values, axis=0)
        spread = 1.4826 * np.median(np.abs(values - loc), axis=0)
        for _ in range(3):
            keep = values <= loc + 3.0 * np.maximum(spread, floor)
            clipped = np.where(keep, values, np.nan)
            loc = np.nanmedian(clipped, axis=0)
            spread = 1.4826 * np.nanmedian(np.abs(clipped - loc), axis=0)
        centre[cols] = np.nan_to_num(loc, nan=0.0).astype(np.float32)
        scale[cols] = np.maximum(np.nan_to_num(spread, nan=floor), floor).astype(np.float32)
    return centre, scale


def find_candidates(run: Run, fold: Fold,
                    settings: CandidateSettings = CandidateSettings(),
                    floor: float | None = None) -> CandidateSet:
    """Find local per-ion peak evidence and group agreeing apexes.

    When ``run`` is processed data, pass the raw acquisition's
    ``storage_threshold(raw) / 3`` as ``floor``.  Falling back to the supplied
    run is convenient for raw and synthetic data but cannot recover provenance
    lost by a processing layer.
    """
    s = settings
    peak_geometry = PeakSettings(wrap_below_s=s.wrap_below_s)
    frame = make_frame(fold, peak_geometry)
    floor = storage_threshold(run) / 3.0 if floor is None else float(floor)

    counts = np.diff(run.data.tocsc().indptr)
    idx = np.flatnonzero(counts >= max(int(s.min_ion_points), 1))
    if idx.size == 0:
        return _empty(frame)
    ions = run.mz[idx]
    cube = frame.extend(fold.cube(run, mz=ions, dtype=np.float32))

    fine = (max(s.smooth_t2_s / frame.pitch, 0.0),
            max(s.smooth_modulations, 0.0), 0.0)
    smoothed = _gaussian_filter_ions(cube, fine) if fine[:2] != (0.0, 0.0) else cube
    centre, scale = _block_location_scale(smoothed, frame, s.noise_block, floor)

    broad = (max(s.broad_t2_s / frame.pitch, fine[0] + 1e-6),
             max(s.broad_modulations, fine[1] + 1e-6), 0.0)
    local_background = _gaussian_filter_ions(smoothed, broad)
    denominator = scale[None, :, :]
    if s.poisson_noise:
        # A fixed storage floor badly overstates significance on an abundant
        # background ion.  For count-like detector signal, shot noise grows
        # approximately with sqrt(intensity).  The empirical spread still wins
        # wherever it is larger.  This remains an assumption to validate against
        # replicate/blank data, so it is explicit in the settings.
        denominator = np.maximum(denominator, np.sqrt(np.maximum(local_background, 0.0)))
    evidence = np.maximum(smoothed - local_background, 0.0) / denominator

    # Ion agreement must also make a locally peak-shaped contribution to the
    # total signal.  This rejects broad spectral background whose many ions can
    # share small extrema without forming a chromatographic peak in the TIC.
    tic = smoothed.sum(axis=2, dtype=np.float64).astype(np.float32)
    tic_centre, tic_scale = _block_location_scale(tic[:, :, None], frame,
                                                   s.noise_block, floor)
    tic_background = ndi.gaussian_filter(tic, broad[:2], mode="nearest")
    tic_denominator = tic_scale[None, :, 0]
    if s.poisson_noise:
        tic_denominator = np.maximum(tic_denominator,
                                     np.sqrt(np.maximum(tic_background, 0.0)))
    tic_evidence = np.maximum(tic - tic_background, 0.0) / tic_denominator
    tic_support = ndi.maximum_filter(
        tic_evidence,
        size=(2 * max(int(s.apex_rows), 0) + 1,
              2 * max(int(s.apex_modulations), 0) + 1), mode="constant")

    canonical = frame.canonical()
    local_max = evidence == ndi.maximum_filter(evidence, size=(5, 3, 1), mode="nearest")
    ion_seeds = local_max & (evidence >= s.seed_snr) & canonical[:, :, None]
    if not ion_seeds.any():
        return _empty(frame)

    neighbourhood = (2 * max(int(s.apex_rows), 0) + 1,
                     2 * max(int(s.apex_modulations), 0) + 1, 1)
    nearby = ndi.maximum_filter(ion_seeds.astype(np.uint8), size=neighbourhood,
                                mode="constant") > 0
    support = nearby.sum(axis=2).astype(np.int32)
    seeded_strength = np.where(ion_seeds, evidence, 0.0)
    strength = ndi.maximum_filter(seeded_strength, size=neighbourhood,
                                  mode="constant").max(axis=2)
    # A grouped candidate must remain on an observed ion apex.  The expanded
    # neighbourhood is only for counting agreeing ions; allowing every nearby
    # pixel to compete would shift a flat score plateau by up to the grouping
    # radius.
    has_seed = ion_seeds.any(axis=2)
    eligible = (has_seed & canonical & (tic_support >= s.tic_prominence_snr) &
                ((support >= s.min_support_ions) | (strength >= s.single_ion_snr)))
    score_surface = np.where(eligible, strength + 2.0 * np.maximum(support - 1, 0), 0.0)
    maxima = eligible & (score_surface == ndi.maximum_filter(
        score_surface, size=(2 * s.apex_rows + 1, 2 * s.apex_modulations + 1), mode="constant"))

    # Collapse equal-score plateaux deterministically to their strongest ion
    # seed (then top-left as a stable final tie-break).
    plateau, n_plateaux = ndi.label(maxima, structure=ndi.generate_binary_structure(2, 2))
    points = []
    for label in range(1, n_plateaux + 1):
        rows, cols = np.nonzero(plateau == label)
        if rows.size == 0:
            continue
        values = score_surface[rows, cols]
        best = int(np.argmax(values))
        points.append((int(rows[best]), int(cols[best])))
    if not points:
        return _empty(frame)

    rows = np.asarray([p[0] for p in points], dtype=np.int64)
    cols = np.asarray([p[1] for p in points], dtype=np.int64)
    order = np.lexsort((rows, cols))
    rows, cols = rows[order], cols[order]
    support_ions, support_snr = [], []
    counts_out, strongest, scores = [], [], []
    ar, ac = max(int(s.apex_rows), 0), max(int(s.apex_modulations), 0)
    for row, col in zip(rows, cols):
        r0, r1 = max(row-ar, 0), min(row+ar+1, evidence.shape[0])
        c0, c1 = max(col-ac, 0), min(col+ac+1, evidence.shape[1])
        sub_seeds = ion_seeds[r0:r1, c0:c1]
        supported = np.flatnonzero(sub_seeds.any(axis=(0, 1)))
        values = np.array([evidence[r0:r1, c0:c1, k][sub_seeds[:, :, k]].max()
                           for k in supported], dtype=float)
        support_ions.append(ions[supported].copy())
        support_snr.append(values)
        counts_out.append(len(supported))
        strongest.append(float(values.max()))
        scores.append(float(values.max() + 2.0 * max(len(supported)-1, 0)))

    return CandidateSet(
        t1=frame.mid_t1[cols].astype(float), t2=(rows * frame.pitch).astype(float),
        row=rows, col=cols, score=np.asarray(scores), strongest_snr=np.asarray(strongest),
        support_count=np.asarray(counts_out, dtype=np.int32), support_ions=support_ions,
        support_snr=support_snr, pitch=frame.pitch, period=frame.period,
        t1_start=float(fold.t1[0]), wrap_rows=frame.wrap_rows)


def _empty(frame: Frame) -> CandidateSet:
    return CandidateSet(*(np.array([], dtype=float) for _ in range(2)),
                        *(np.array([], dtype=np.int64) for _ in range(2)),
                        *(np.array([], dtype=float) for _ in range(2)),
                        np.array([], dtype=np.int32), [], [], frame.pitch,
                        frame.period, float(frame.fold.t1[0]), frame.wrap_rows)


def candidate_regions(candidates: CandidateSet, image_shape: tuple[int, int],
                      link_t2_s: float = 0.25, link_modulations: int = 1,
                      pad_t2_s: float = 0.15, pad_modulations: int = 1) -> list[CandidateRegion]:
    """Join nearby apex hypotheses and give each group a padded fit rectangle."""
    n = len(candidates)
    if n == 0:
        return []
    link_rows = max(int(np.ceil(link_t2_s / candidates.pitch)), 0)
    # Complete linkage prevents a long chain of individually-near candidates
    # from becoming one enormous region.  Such chains previously produced, for
    # example, a 365-row/64-candidate window and were both slow and scientifically
    # incompatible with a local-background model.  Chebyshev distance preserves
    # the old rectangular t1/t2 limits while requiring the whole group to fit.
    if n == 1:
        labels = np.ones(1, dtype=int)
    else:
        scaled = np.column_stack((
            candidates.col / max(int(link_modulations), 1),
            candidates.row / max(int(link_rows), 1)))
        # A complete-linkage cluster at t=1 has every pair within distance 1, so
        # it lies inside one connected component of the "within 1" graph. Cluster
        # each component alone: same partition, but linkage's O(n^2) distance
        # matrix is only ever built for one small component, not all candidates.
        pairs = cKDTree(scaled).query_pairs(1.0, p=np.inf, output_type="ndarray")
        graph = coo_matrix((np.ones(len(pairs), dtype=np.int8), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
        n_comp, comp = connected_components(graph, directed=False)
        labels = np.zeros(n, dtype=int)
        next_label = 1
        for c in range(n_comp):
            members = np.flatnonzero(comp == c)
            if len(members) == 1:
                labels[members] = next_label
                next_label += 1
                continue
            tree = linkage(scaled[members], method="complete", metric="chebyshev")
            sub = fcluster(tree, t=1.0, criterion="distance")
            labels[members] = sub + next_label - 1
            next_label += int(sub.max())
    groups = {}
    for i, label in enumerate(labels):
        groups.setdefault(int(label), []).append(i)
    pad_rows = max(int(np.ceil(pad_t2_s / candidates.pitch)), 0)
    n_rows, n_cols = image_shape[:2]
    out = []
    for indices in groups.values():
        rows, cols = candidates.row[indices], candidates.col[indices]
        out.append(CandidateRegion(
            tuple(indices), max(int(rows.min()) - pad_rows, 0),
            min(int(rows.max()) + pad_rows + 1, n_rows),
            max(int(cols.min()) - pad_modulations, 0),
            min(int(cols.max()) + pad_modulations + 1, n_cols)))
    return sorted(out, key=lambda r: (r.col0, r.row0))


def candidate_cube(run: Run, fold: Fold, candidates: CandidateSet,
                   region: CandidateRegion, max_ions: int = 60, min_ions: int = 8):
    """Build an overlap-compatible cube without inheriting a legacy border.

    Ions which independently supported a candidate are retained first.  Any
    remaining capacity is filled by ions with the largest local variation, so
    weak diagnostic structure is not excluded merely because it forms a small
    fraction of the region's total spectrum.
    """
    from .overlap import BlobCube

    supported_mz = {}
    for i in region.candidate_indices:
        for mz, snr in zip(candidates.support_ions[i], candidates.support_snr[i]):
            supported_mz[int(mz)] = max(float(snr), supported_mz.get(int(mz), 0.0))
    ranked_support = sorted(supported_mz, key=supported_mz.get, reverse=True)
    supported_idx = [int(m - run.mz[0]) for m in ranked_support
                     if run.mz[0] <= m <= run.mz[-1]]
    # Rank supplemental ions from the sparse scans in these modulation columns.
    # The old implementation folded every ion over the entire acquisition for
    # every region, then discarded all but 60.  That repeated dense allocation
    # dominated full-run time and memory.  Summing the local sparse scan block
    # selects useful abundant ions without constructing the all-ion cube.
    t1_lo = candidates.t1_start + region.col0 * candidates.period
    t1_hi = candidates.t1_start + region.col1 * candidates.period
    scan0, scan1 = np.searchsorted(run.times, (t1_lo, t1_hi))
    if scan1 > scan0:
        local_total = np.asarray(run.data[scan0:scan1].sum(axis=0)).ravel()
    else:
        local_total = np.zeros(run.n_mz, dtype=float)
    order = np.argsort(local_total)[::-1]
    chosen = list(dict.fromkeys(supported_idx + order.tolist()))
    positive = [i for i in chosen if local_total[i] > 0 or i in supported_idx]
    limit = max(int(max_ions), int(min_ions))
    chosen = np.asarray(positive[:limit], dtype=np.int64)
    if chosen.size < min(min_ions, int(np.count_nonzero(local_total > 0))):
        chosen = order[local_total[order] > 0][:min_ions]
    data = _fold_region(run, fold, chosen, region, candidates.wrap_rows).astype(
        np.float64, copy=False)
    mask = np.ones(data.shape[:2], dtype=bool)
    return BlobCube(data=data, mask=mask, ions=run.mz[chosen],
                    row0=region.row0, col0=region.col0, pitch=candidates.pitch,
                    period=candidates.period, t1_start=candidates.t1_start)


def _fold_region(run: Run, fold: Fold, chosen: np.ndarray, region: CandidateRegion,
                 wrap_rows: int) -> np.ndarray:
    """Fold selected ions directly into one extended local region."""
    n_rows = fold.n_rows
    out_rows = region.row1 - region.row0
    out_cols = region.col1 - region.col0
    out = np.zeros((out_rows, out_cols, len(chosen)), dtype=np.float32)
    if not len(chosen):
        return out

    def scan_hits(flat_pixels: np.ndarray) -> np.ndarray:
        valid = flat_pixels >= 0
        rows = np.where(valid, flat_pixels % n_rows, -1)
        cols = np.where(valid, flat_pixels // n_rows, -1)
        regular = (valid & (cols >= region.col0) & (cols < region.col1) &
                   (rows >= region.row0) & (rows < min(region.row1, n_rows)))
        extended = (valid & (rows < wrap_rows) &
                    (cols - 1 >= region.col0) & (cols - 1 < region.col1) &
                    (n_rows + rows >= region.row0) &
                    (n_rows + rows < region.row1))
        return regular | extended

    # Only scans in these modulation columns can land in the local rectangle.
    # Testing every scan in the acquisition for every region made run-wide
    # fitting scale roughly as regions x acquisition length.  Keep one extra
    # modulation on either side for interpolation and the extended wrap rows.
    period = fold.spec.period
    time_lo = fold.t1[0] + (region.col0 - 1) * period
    time_hi = fold.t1[0] + (region.col1 + 1) * period
    scan0, scan1 = np.searchsorted(run.times, (time_lo, time_hi))
    nearby_scans = np.arange(max(int(scan0), 0), min(int(scan1), run.n_scans))
    scan_mask = (scan_hits(fold._lo[nearby_scans]) |
                 scan_hits(fold._up[nearby_scans]))
    scan_indices = nearby_scans[scan_mask]
    if not scan_indices.size:
        return out
    sub = run.data[scan_indices][:, chosen].tocoo()
    scans = scan_indices[sub.row]
    values = sub.data.astype(np.float64)

    def add(flat_pixels: np.ndarray, weights: np.ndarray) -> None:
        flat = flat_pixels[scans]
        valid = flat >= 0
        rows = flat % n_rows
        cols = flat // n_rows
        regular = (valid & (cols >= region.col0) & (cols < region.col1) &
                   (rows >= region.row0) & (rows < min(region.row1, n_rows)))
        if regular.any():
            np.add.at(out, (rows[regular] - region.row0,
                            cols[regular] - region.col0, sub.col[regular]),
                      (values[regular] * weights[scans[regular]]).astype(np.float32))
        extended_rows = n_rows + rows
        extended_cols = cols - 1
        extended = (valid & (rows < wrap_rows) &
                    (extended_cols >= region.col0) & (extended_cols < region.col1) &
                    (extended_rows >= region.row0) & (extended_rows < region.row1))
        if extended.any():
            np.add.at(out, (extended_rows[extended] - region.row0,
                            extended_cols[extended] - region.col0, sub.col[extended]),
                      (values[extended] * weights[scans[extended]]).astype(np.float32))

    add(fold._lo, 1.0 - fold._w)
    add(fold._up, fold._w)
    return out
