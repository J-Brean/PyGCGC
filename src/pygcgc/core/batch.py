"""Batch processing across a whole sequence: a cached pipeline, a feature
table matched across runs, QC control charts, and simple cross-run
statistics.

The pipeline (`run_batch`)
---------------------------
The expensive, embarrassingly parallel part of processing a sequence is
per run: fold, 🧹 clean, optionally 🫧 and 🔍. Those steps run in a
process pool (real multiprocessing -- the GIL means a thread pool, used
elsewhere in this codebase for GUI responsiveness, doesn't actually
parallelise CPU-bound work like this), one worker per run, each step's
result cached to disk keyed by a hash of its inputs and settings. A step
whose inputs and settings haven't changed since the last run loads
straight from the cache instead of recomputing, which is also the whole
resumption story: an interrupted batch, re-run with the same cache
directory, just picks up wherever it left off, because everything already
cached loads instantly and everything not yet cached still needs doing.
One run's failure doesn't stop the others; it's reported and skipped.

Calibration and quantification (`quantify_prepared`, in `quant.py`) is
run once, in the main process, on the results of that parallel prep --
it needs every run together (to learn positions from the calibration
runs and fit calibration curves), and is fast enough on top of an
already-built set of images that there's nothing to gain by
parallelising it too.

Cross-run outputs
------------------
`build_feature_table` matches every other run's blobs onto one reference
run's, by the Hungarian algorithm on the same scaled 2D distance used
elsewhere (`BlobSet.nearest`'s convention); a run with no match within
`max_distance` for a given feature is "gap-filled" instead -- a small box
integrated directly at the feature's expected position in that run's own
cleaned image, on the assumption that something is very probably there
even though automatic detection missed it.

`westgard_flags` applies the standard multi-rule QC scheme (1-2s, 1-3s,
2-2s, R-4s, 4-1s, 10x) to one series of values (an internal standard's
volume across a sequence, say).

`pca` is a plain SVD-based PCA (no new dependency) on the feature table,
log-transformed and autoscaled first, the standard preprocessing for this
kind of data. `volcano` is a two-group Welch's t-test per feature (unequal
variance, since there's no reason to assume two sample groups share one),
with Benjamini-Hochberg FDR correction across every feature tested.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
import pickle
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import ttest_ind

from .background import BackgroundResult, BackgroundSettings, NoiseMap, remove_background, storage_threshold
from .fold import Fold, FoldSpec
from .io import read_many, discover_runs
from .peaks import BlobSet, PeakSettings, find_blobs
from .period import estimate_period_for_runs
from .phase import phase_for_sequence
from .quant import QuantReport, QuantSettings, SequenceRun, quantify_prepared, read_sequence_extras
from .run import Run
from .siloxane import SiloxaneResult, SiloxaneSettings, remove_siloxane_blobs, remove_siloxanes
from .targets import load_anchors, load_levels, load_targets


@dataclass(frozen=True)
class BatchSettings:
    background: BackgroundSettings = field(default_factory=BackgroundSettings)
    run_siloxane: bool = False
    siloxane: SiloxaneSettings = field(default_factory=SiloxaneSettings)
    run_blobs: bool = True
    peaks: PeakSettings = field(default_factory=PeakSettings)
    run_quant: bool = True
    quant: QuantSettings = field(default_factory=QuantSettings)

    def to_dict(self) -> dict:
        return dict(background=self.background.to_dict(), run_siloxane=self.run_siloxane,
                   siloxane=self.siloxane.to_dict(), run_blobs=self.run_blobs,
                   peaks=self.peaks.to_dict(), run_quant=self.run_quant, quant=self.quant.to_dict())

    @classmethod
    def from_dict(cls, d: dict) -> "BatchSettings":
        d = d or {}
        return cls(background=BackgroundSettings.from_dict(d.get("background", {})),
                  run_siloxane=bool(d.get("run_siloxane", False)),
                  siloxane=SiloxaneSettings.from_dict(d.get("siloxane", {})),
                  run_blobs=bool(d.get("run_blobs", True)),
                  peaks=PeakSettings.from_dict(d.get("peaks", {})),
                  run_quant=bool(d.get("run_quant", True)),
                  quant=QuantSettings.from_dict(d.get("quant", {})))


@dataclass(frozen=True)
class BatchConfig:
    """The saved pipeline: settings plus which runs go where for the
    cross-run outputs. `period=None` reads `method.json`."""
    period: float | None = None
    settings: BatchSettings = field(default_factory=BatchSettings)
    reference: str = ""                  # run for the feature table; "" picks one automatically
    group_a: tuple[str, ...] = ()        # for the volcano plot
    group_b: tuple[str, ...] = ()

    def to_dict(self) -> dict:
        return dict(period=self.period, settings=self.settings.to_dict(), reference=self.reference,
                   group_a=list(self.group_a), group_b=list(self.group_b))

    @classmethod
    def from_dict(cls, d: dict) -> "BatchConfig":
        d = d or {}
        return cls(period=d.get("period"), settings=BatchSettings.from_dict(d.get("settings", {})),
                  reference=d.get("reference", ""), group_a=tuple(d.get("group_a", [])),
                  group_b=tuple(d.get("group_b", [])))

    def save(self, path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=2))

    @classmethod
    def load(cls, path) -> "BatchConfig":
        return cls.from_dict(json.loads(Path(path).read_text()))


# ------------------------------------------------------------------ cache
def _hash_key(*parts) -> str:
    h = hashlib.sha256()
    for p in parts:
        h.update(repr(p).encode())
    return h.hexdigest()[:20]


class StepCache:
    """One pickled file per cache key, under `directory`."""

    def __init__(self, directory: Path):
        self.dir = Path(directory)
        self.dir.mkdir(parents=True, exist_ok=True)

    def get(self, key: str):
        p = self.dir / f"{key}.pkl"
        if not p.exists():
            return None
        try:
            with open(p, "rb") as fh:
                return pickle.load(fh)
        except (pickle.PickleError, EOFError, OSError):
            return None

    def set(self, key: str, value) -> None:
        p = self.dir / f"{key}.pkl"
        tmp = p.with_suffix(".tmp")
        with open(tmp, "wb") as fh:
            pickle.dump(value, fh)
        tmp.replace(p)


# ------------------------------------------------------------------- prep
@dataclass
class PrepResult:
    """One run's output from the parallel, cached part of the pipeline."""
    name: str
    run: Run
    clean: Run
    fold: Fold
    noise: NoiseMap
    background: BackgroundResult
    siloxane: SiloxaneResult | None = None
    blobs: BlobSet | None = None


def _run_cache_base(run: Run, period: float, phase: float) -> str:
    stat = run.source.stat() if run.source and Path(run.source).exists() else None
    fingerprint = (stat.st_mtime_ns, stat.st_size) if stat else (run.nnz, float(run.data.sum()))
    return _hash_key(str(run.source), fingerprint, period, phase,
                     {k:run.meta.get(k) for k in ('source_sha256', 'parser_version',
                      'energy_ev', 'processing_state', 'time_origin_policy', 'time_origin_seconds')})


def _prep_one(run: Run, period: float, phase: float, settings: BatchSettings,
             cache_dir: str | None) -> PrepResult:
    """The per-run pipeline step: what a process-pool worker (or, with
    `max_workers=1`, the main process directly) actually runs."""
    fold = Fold(run.times, FoldSpec(period=period, phase=phase))
    cache = StepCache(cache_dir) if cache_dir else None
    base = _run_cache_base(run, period, phase)

    bg_key = _hash_key(base, "background", settings.background.to_dict())
    bg = cache.get(bg_key) if cache else None
    if bg is None:
        bg = remove_background(run, fold, settings.background)
        if cache:
            cache.set(bg_key, bg)
    clean, chain_key = bg.run, bg_key

    sil = None
    if settings.run_siloxane:
        sil_key = _hash_key(chain_key, "siloxane", settings.siloxane.to_dict(),
                            settings.peaks.to_dict() if settings.siloxane.method == "blobs" else None)
        sil = cache.get(sil_key) if cache else None
        if sil is None:
            if settings.siloxane.method == "blobs":
                # The blobs to block out have to be found before the siloxane is gone.
                pre_blobs = find_blobs(clean, fold, settings.peaks, floor=storage_threshold(run) / 3.0)
                sil = remove_siloxane_blobs(clean, fold, pre_blobs, settings.siloxane)
            else:
                sil = remove_siloxanes(clean, fold, settings.siloxane)
            if cache:
                cache.set(sil_key, sil)
        clean, chain_key = sil.run, sil_key

    blobs = None
    if settings.run_blobs:
        raw_floor = storage_threshold(run) / 3.0
        blob_key = _hash_key(chain_key, "blobs-v2-raw-floor", settings.peaks.to_dict(), raw_floor)
        blobs = cache.get(blob_key) if cache else None
        if blobs is None:
            blobs = find_blobs(clean, fold, settings.peaks, floor=raw_floor)
            if cache:
                cache.set(blob_key, blobs)

    return PrepResult(name=run.name, run=run, clean=clean, fold=fold, noise=bg.noise,
                      background=bg, siloxane=sil, blobs=blobs)


# --------------------------------------------------------------- pipeline
@dataclass
class BatchResult:
    prepped: dict                 # run name -> PrepResult
    errors: dict                  # run name -> error message
    report: QuantReport | None
    feature_table: list | None
    reference: str = ""
    period: float | None = None   # the period used, including one estimated from the runs
    period_estimated: bool = False


def run_batch(folder, config: BatchConfig = BatchConfig(), cache_dir=None, max_workers: int | None = None,
             progress=None) -> BatchResult:
    """Run the whole pipeline on a folder in the dummy-pack layout
    (`sequence.csv`, `target_table.csv`, calibration levels, `method.json`
    -- see README.md). One run's failure is reported, not fatal.
    `max_workers=1` skips the process pool entirely (for tests, and any
    environment where spawning processes is unwelcome); otherwise
    `cores - 1` by default.
    """
    folder = Path(folder)
    period = config.period
    if period is None:
        method = folder / "method.json"
        if method.exists():
            period = float(json.loads(method.read_text())["modulation_period_s"])

    targets_path = folder / "target_table.csv"
    seq_path = folder / "sequence.csv"
    levels_path = next((folder / n for n in ("calibration_levels.csv", "calibration_levels_ng.csv")
                       if (folder / n).exists()), None)
    targets = load_targets(targets_path) if targets_path.exists() else []
    levels = load_levels(levels_path) if levels_path else {}
    anchors = load_anchors(targets_path) if targets_path.exists() else None
    extras = read_sequence_extras(seq_path) if seq_path.exists() else {}

    if seq_path.exists():
        rows = list(csv.DictReader(open(seq_path, newline="", encoding="utf-8-sig")))
        files = []
        for row in rows:
            hits = list(folder.rglob(row["file"]))
            if len(hits) != 1:
                raise ValueError(f"expected one match for {row['file']}; found {len(hits)}")
            files.append(hits[0])
    else:
        files = discover_runs(folder)
        rows = [{"file": p.name, "type": "sample", "level": ""} for p in files]

    step = progress or (lambda *a: None)
    step(0, len(files) + 1, "Reading")
    raw_runs = read_many(files)
    estimated = period is None
    if estimated:
        best, found = estimate_period_for_runs(raw_runs)
        if best is None:
            guess = f" (best guess {found[0].period:.4f} s, low confidence)" if found else ""
            raise ValueError(f"no period given, no method.json, and none could be estimated{guess}")
        period = best.period
        step(0, len(files) + 1, best.summary())
    phases = phase_for_sequence(raw_runs, period, anchors)

    settings = config.settings
    cache_str = str(cache_dir) if cache_dir else None
    prepped: dict = {}
    errors: dict = {}
    n = len(raw_runs)
    max_workers = (1 if any(r.meta.get('format') == 'Markes LSC' for r in raw_runs)
                   else max(1, (os.cpu_count() or 2) - 1)) if max_workers is None else max_workers

    if max_workers <= 1:
        for k, run in enumerate(raw_runs):
            step(k + 1, n + 1, f"Preparing {run.name}")
            try:
                prepped[run.name] = _prep_one(run, period, phases[k].phase, settings, cache_str)
            except Exception as exc:  # noqa: BLE001
                errors[run.name] = str(exc)
    else:
        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            futures = {ex.submit(_prep_one, raw_runs[k], period, phases[k].phase, settings, cache_str): k
                      for k in range(n)}
            done = 0
            for fut in as_completed(futures):
                k = futures[fut]
                done += 1
                step(done, n + 1, f"Prepared {raw_runs[k].name}")
                try:
                    prepped[raw_runs[k].name] = fut.result()
                except Exception as exc:  # noqa: BLE001
                    errors[raw_runs[k].name] = str(exc)

    report = None
    if settings.run_quant and targets and levels:
        prepared = []
        info = []
        for row in rows:
            name = Path(row["file"]).stem
            if name not in prepped:
                continue
            p = prepped[name]
            prepared.append((p.run, p.clean, p.fold, p.noise))
            ex_ = extras.get(name, {})
            info.append(SequenceRun(name, row.get("type", ""), (row.get("level") or "").strip(),
                                    ex_.get("cal_set", ""), ex_.get("volume")))
        if prepared:
            report = quantify_prepared(prepared, info, targets, levels, settings.quant)

    reference = config.reference or next((n for n in prepped if prepped[n].blobs is not None), "")
    feature_table = None
    if settings.run_blobs and reference and reference in prepped:
        feature_table = build_feature_table(prepped, reference)

    step(n + 1, n + 1, "Done")
    return BatchResult(prepped=prepped, errors=errors, report=report, feature_table=feature_table,
                       reference=reference, period=period, period_estimated=estimated)


# ------------------------------------------------------------ feature table
@dataclass
class Feature:
    id: int
    ref_t1: float
    ref_t2: float
    volumes: dict = field(default_factory=dict)     # run name -> volume
    filled: dict = field(default_factory=dict)       # run name -> gap-filled, not matched


def _box_sum(run: Run, fold: Fold, t1: float, t2: float, t1_half_s: float = 6.0,
            t2_half_s: float = 0.15) -> float:
    """A small window's TIC sum around an expected position, for filling
    in a feature that a run's own blob detection missed."""
    img = fold.tic_image(run)
    col = int(round((t1 - fold.t1[0]) / fold.spec.period - 0.5))
    row = int(round(t2 / fold.row_pitch)) % img.shape[0]
    c_half = max(int(round(t1_half_s / fold.spec.period)), 1)
    r_half = max(int(round(t2_half_s / fold.row_pitch)), 1)
    r0, r1 = max(row - r_half, 0), min(row + r_half + 1, img.shape[0])
    c0, c1 = max(col - c_half, 0), min(col + c_half + 1, img.shape[1])
    if r1 <= r0 or c1 <= c0:
        return 0.0
    return float(np.clip(img[r0:r1, c0:c1], 0, None).sum())


def build_feature_table(prepped: dict, reference: str, t1_scale: float = 60.0, t2_scale: float = 0.5,
                        max_distance: float = 1.5) -> list:
    """Every blob in `reference`'s run becomes a feature; every other
    run's own blobs are matched onto them by the Hungarian algorithm (a
    scaled 2D distance, `BlobSet.nearest`'s convention), and a run with no
    match within `max_distance` is gap-filled from its own cleaned image
    instead of left blank."""
    ref = prepped[reference].blobs
    if ref is None or len(ref) == 0:
        return []
    features = [Feature(id=i, ref_t1=float(ref.t1[i]), ref_t2=float(ref.t2[i])) for i in range(len(ref))]
    for name, r in prepped.items():
        if name == reference:
            for i, feat in enumerate(features):
                feat.volumes[name] = float(ref.volume[i])
                feat.filled[name] = False
            continue
        blobs = r.blobs
        matched: dict = {}
        if blobs is not None and len(blobs) > 0:
            cost = np.hypot((ref.t1[:, None] - blobs.t1[None, :]) / t1_scale,
                            (ref.t2[:, None] - blobs.t2[None, :]) / t2_scale)
            row_ind, col_ind = linear_sum_assignment(cost)
            for ri, ci in zip(row_ind, col_ind):
                if cost[ri, ci] <= max_distance:
                    matched[int(ri)] = int(ci)
        for i, feat in enumerate(features):
            if i in matched:
                feat.volumes[name] = float(blobs.volume[matched[i]])
                feat.filled[name] = False
            else:
                feat.volumes[name] = _box_sum(r.clean, r.fold, feat.ref_t1, feat.ref_t2)
                feat.filled[name] = True
    return features


def feature_table_matrix(feature_table: list, run_names: list) -> np.ndarray:
    """(n_runs, n_features) volume matrix, in `run_names`' order."""
    return np.array([[f.volumes.get(name, 0.0) for f in feature_table] for name in run_names], dtype=float)


# ------------------------------------------------------------------ QC
@dataclass
class QCFlag:
    index: int
    rule: str
    message: str


def westgard_flags(values, mean: float | None = None, sd: float | None = None) -> list:
    """The standard multi-rule QC scheme (Westgard et al., 1981):
    1-2s/1-3s (one point beyond 2 or 3 SD), 2-2s (two consecutive points
    beyond 2 SD on the same side), R-4s (a 4 SD range between two
    consecutive points), 4-1s (four consecutive points beyond 1 SD on the
    same side), 10x (ten consecutive points on the same side of the
    mean). `mean`/`sd` default to the series' own; pass a target mean/SD
    to judge against an external reference instead."""
    values = np.asarray(values, dtype=float)
    mean = float(np.mean(values)) if mean is None else mean
    sd = float(np.std(values, ddof=1)) if sd is None else sd
    flags: list = []
    if sd <= 0 or values.size == 0:
        return flags
    z = (values - mean) / sd
    for i in range(z.size):
        if abs(z[i]) > 3:
            flags.append(QCFlag(i, "1-3s", f"point {i} is {z[i]:+.1f} SD from the mean"))
        elif abs(z[i]) > 2:
            flags.append(QCFlag(i, "1-2s", f"point {i} is {z[i]:+.1f} SD from the mean (warning)"))
    for i in range(1, z.size):
        if (z[i] > 2 and z[i - 1] > 2) or (z[i] < -2 and z[i - 1] < -2):
            flags.append(QCFlag(i, "2-2s", f"points {i - 1} and {i} are both beyond 2 SD, same side"))
        if abs(z[i] - z[i - 1]) > 4:
            flags.append(QCFlag(i, "R-4s", f"points {i - 1} and {i} span more than 4 SD"))
    for i in range(3, z.size):
        w = z[i - 3:i + 1]
        if np.all(w > 1) or np.all(w < -1):
            flags.append(QCFlag(i, "4-1s", f"points {i - 3}-{i} are all beyond 1 SD, same side"))
    for i in range(9, z.size):
        w = z[i - 9:i + 1]
        if np.all(w > 0) or np.all(w < 0):
            flags.append(QCFlag(i, "10x", f"points {i - 9}-{i} are all on the same side of the mean"))
    return flags


# ------------------------------------------------------------ statistics
@dataclass
class PCAResult:
    scores: np.ndarray               # (n_runs, n_components)
    loadings: np.ndarray             # (n_features, n_components)
    explained_variance: np.ndarray   # (n_components,), fraction of total
    run_names: list
    feature_ids: list


def pca(feature_table: list, run_names: list, n_components: int = 2, log: bool = True) -> PCAResult:
    """Plain SVD-based PCA (`X = U S Vt`; scores = `U S`, loadings = `Vt.T`)
    on the feature table, log-transformed and autoscaled (mean-centred,
    unit variance per feature) first -- the usual preprocessing for this
    kind of data, so that a handful of huge peaks don't dominate purely
    because of their scale."""
    X = feature_table_matrix(feature_table, run_names)
    if log:
        X = np.log1p(np.clip(X, 0, None))
    mean = X.mean(axis=0)
    sd = X.std(axis=0)
    Xc = (X - mean) / np.where(sd > 0, sd, 1.0)
    u, s, vt = np.linalg.svd(Xc, full_matrices=False)
    k = min(n_components, s.size)
    scores = u[:, :k] * s[:k]
    loadings = vt[:k].T
    total = float((s ** 2).sum())
    explained = (s[:k] ** 2) / total if total > 0 else np.zeros(k)
    return PCAResult(scores=scores, loadings=loadings, explained_variance=explained,
                     run_names=list(run_names), feature_ids=[f.id for f in feature_table])


def benjamini_hochberg(p_values) -> np.ndarray:
    """FDR-adjusted p-values (Benjamini & Hochberg, 1995)."""
    p_values = np.asarray(p_values, dtype=float)
    n = p_values.size
    if n == 0:
        return p_values
    order = np.argsort(p_values)
    ranks = np.arange(1, n + 1)
    adjusted = p_values[order] * n / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out = np.empty(n)
    out[order] = adjusted
    return out


@dataclass
class VolcanoPoint:
    feature_id: int
    ref_t1: float
    ref_t2: float
    log2_fold_change: float
    p_value: float
    q_value: float


def volcano(feature_table: list, group_a: list, group_b: list, min_group: int = 2) -> list:
    """Welch's t-test (unequal variance) per feature between two groups of
    runs, log2 fold change of the (offset) group means, and
    Benjamini-Hochberg FDR correction across every feature actually
    tested. A feature needs at least `min_group` non-zero values in each
    group to be tested at all."""
    raw = []
    for f in feature_table:
        a = np.array([f.volumes.get(n, 0.0) for n in group_a])
        b = np.array([f.volumes.get(n, 0.0) for n in group_b])
        if int((a > 0).sum()) < min_group or int((b > 0).sum()) < min_group:
            continue
        _t, p = ttest_ind(a, b, equal_var=False)
        fc = float(np.log2((a.mean() + 1.0) / (b.mean() + 1.0)))
        raw.append((f, fc, float(p)))
    if not raw:
        return []
    q = benjamini_hochberg(np.array([r[2] for r in raw]))
    return [VolcanoPoint(feature_id=f.id, ref_t1=f.ref_t1, ref_t2=f.ref_t2, log2_fold_change=fc,
                         p_value=p, q_value=float(qi)) for (f, fc, p), qi in zip(raw, q)]
