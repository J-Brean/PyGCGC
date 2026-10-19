"""🧩 Peak-clustering state for the app.

Holds the settings, and the latest clustering: one `PeakTable` of the peaks
that went in (one run, or every 70 eV run pooled), the `ClusterResult`, and
the blobs each peak came from, so the 2D view can outline them. Results are
recomputed on demand and not saved, like 🔍's blobs and 🧬's classes; the
settings are saved. The science is all in `core/cluster.py`.

Peaks come from the legacy detector's blobs. A run with no 🔍 blobs (an
MS-aware user) gets the legacy detector run here, privately, on the same
cleaned run, as 🧬 does.
"""

from __future__ import annotations

from dataclasses import replace

from PySide6 import QtCore

import numpy as np

from ..core.background import remove_background, storage_threshold
from ..core.batch import Feature, feature_table_matrix, volcano
from ..core.chemclass import ClassSettings, soft_spectra_for_blobs
from ..core.cluster import (ClusterSettings, ClusterSummary, PeakTable, apply_pqn, class_composition, cluster_features,
                            cluster_label_image, cluster_peaks, clusters_of_run, merge_clusters, permanova, stability,
                            threshold_scan)
from ..core.fold import Fold
from ..core.peaks import find_blobs
from ..core.peaktable import respectrum_components
from ..core.segment import SegmentSettings, fit_regions, region_shares
from ..core.siloxane import SILOXANE_IONS
from .blob_state import MS_AWARE

HARD_ENERGIES = (None, 70.0)
SCOPES = ("run", "all")
UNITS = ("peaks", "regions")


class ClusterState(QtCore.QObject):
    changed = QtCore.Signal()          # a result, the settings or its staleness changed
    runningChanged = QtCore.Signal(bool)

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.settings = ClusterSettings()
        self.skip_siloxane = False            # leave the siloxane ions out of the distance
        self.scope = "run"                    # "run": the selected run; "all": every 70 eV run, pooled
        self.running = False
        self._input_revision = 0
        self.unit = "peaks"                   # "peaks": cluster detected peaks; "regions": segment the pixels themselves
        self.n_regions = 12
        self.normalise = "percent"            # run table: "percent" of the run's own total, or "pqn" (see apply_pqn)
        self.names = {"peaks": {}, "regions": {}}   # cluster or region -> the name the user gave it (not saved: numbers change)
        self._stale = {"peaks": False, "regions": False}
        self.regions = None                   # SegmentResult (fitted on every run in `region_runs`)
        self.region_labels: dict = {}         # run name -> label image
        self.region_shares: dict = {}         # run name -> percent of its signal in each region
        self.region_summaries: dict = {}      # run name -> [RegionSummary]
        self.region_runs: list = []
        self.region_notes = ""
        self.table: PeakTable | None = None
        self.result = None                    # ClusterResult for `table`
        self.summary: list = []               # ClusterSummary rows, largest first
        self.blobs_by_run: dict = {}          # run name -> the BlobSet its peaks came from
        self.run_names: list = []             # runs in the result, in order
        self.notes = ""
        self.use_soft = False                 # add the 14 eV spectrum of a tandem partner run to the distance
        self.groups: dict = {}                # run name -> the user's group label, for the run comparison
        self.check_text = ""                  # result of the last stability or threshold check

    # The staleness of whichever unit is showing; setting it sets that unit's.
    @property
    def stale(self) -> bool:
        return self._stale[self.unit]

    @stale.setter
    def stale(self, value: bool) -> None:
        self._stale[self.unit] = bool(value)

    # ------------------------------------------------------------ queries
    def has_result(self) -> bool:
        return (self.regions is not None) if self.unit == "regions" else (self.result is not None)

    def n_groups(self) -> int:
        """Clusters or regions in the result being shown."""
        if self.unit == "regions":
            return 0 if self.regions is None else len(self.regions)
        return 0 if self.result is None else len(self.result)

    def top_ions(self, k: int, n: int = 3) -> list:
        source = self.regions if self.unit == "regions" else self.result
        return source.top_ions(k, n)

    def active_runs(self) -> list:
        return self.region_runs if self.unit == "regions" else self.run_names

    def label_image(self, index: int):
        """Cluster or region of every pixel of run `index` on the plain folded grid (-1 outside), or None."""
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        if self.unit == "regions":
            return self.region_labels.get(p.entries[index].name)
        blobs, per_blob = self.blobs_for(index), self.blob_clusters(index)
        return None if blobs is None or per_blob is None else cluster_label_image(blobs, per_blob)

    def current_summary(self, index: int) -> list:
        """`ClusterSummary` rows for the table: the peak clusters, or run `index`'s regions in the same shape
        (pixels stand in for peaks)."""
        if self.unit == "peaks":
            return self.summary
        p = self.project
        rows = self.region_summaries.get(p.entries[index].name, []) if 0 <= index < len(p.entries) else []
        return [ClusterSummary(label=r.label, n_peaks=r.n_pixels, n_runs=len(self.region_runs), volume=float("nan"),
                               percent=r.percent, t1_median=r.t1_median, t2_median=r.t2_median, top_ions=r.top_ions,
                               silhouette=float("nan")) for r in rows]

    def energy_of(self, index: int):
        return (self.project.entries[index].import_meta or {}).get("energy_ev")

    def eligible(self) -> list[int]:
        """Loaded runs whose spectra are 70 eV (or of unknown energy): clusters are read from the fragmentation pattern."""
        p = self.project
        return [i for i, e in enumerate(p.entries) if e.run is not None and self.energy_of(i) in HARD_ENERGIES]

    def why_not(self, index: int) -> str:
        """Empty if the current settings can be run, otherwise the reason."""
        p = self.project
        if self.scope == "all":
            return "" if self.eligible() else "Open at least one 70 eV run."
        if not 0 <= index < len(p.entries) or p.entries[index].run is None:
            return "Select a loaded run."
        energy = self.energy_of(index)
        if energy not in HARD_ENERGIES:
            return f"This is the {energy:g} eV channel. Select the 70 eV run."
        return ""

    def effective_settings(self) -> ClusterSettings:
        return replace(self.settings, exclude_ions=tuple(SILOXANE_IONS) if self.skip_siloxane else (),
                       soft_weight=0.5 if self.use_soft else 0.0)

    def rows_for(self, index: int):
        """(blob index, cluster) arrays for run `index`'s peaks in the current result, or None."""
        p = self.project
        if self.result is None or self.table is None or not 0 <= index < len(p.entries):
            return None
        mine = self.table.run == p.entries[index].name
        if not mine.any():
            return None
        return self.table.blob[mine].astype(int), self.result.labels[mine]

    def blobs_for(self, index: int):
        p = self.project
        return self.blobs_by_run.get(p.entries[index].name) if 0 <= index < len(p.entries) else None

    def blob_clusters(self, index: int):
        """Cluster of every blob of run `index` (-1 for none), or None if it is not in the current result."""
        blobs = self.blobs_for(index)
        if blobs is None or self.result is None or self.table is None:
            return None
        return clusters_of_run(self.table, self.result, self.project.entries[index].name, len(blobs))

    def fingerprint(self):
        """(features, percent matrix) or None: each run's share in each cluster (of its peak volume) or region (of
        its signal). Rows of the matrix are runs, in `active_runs()` order; columns are clusters or regions."""
        if self.unit == "regions":
            if self.regions is None or not self.region_runs:
                return None
            first = self.region_summaries[self.region_runs[0]]
            feats = [Feature(id=k, ref_t1=first[k].t1_median, ref_t2=first[k].t2_median,
                             volumes={n: float(self.region_shares[n][k]) for n in self.region_runs})
                     for k in range(len(self.regions))]
            if self.normalise == "pqn":
                feats = apply_pqn(feats, self.region_runs)
            return feats, feature_table_matrix(feats, self.region_runs)
        if self.result is None or self.table is None or len(self.result) == 0:
            return None
        feats = cluster_features(self.table, self.result, self.run_names, relative=True, pqn=self.normalise == "pqn")
        return feats, feature_table_matrix(feats, self.run_names)

    def name_map(self) -> dict:
        """cluster or region -> name, for those the user has named."""
        return {k: v for k, v in self.names[self.unit].items() if v}

    def set_name(self, label: int, name: str) -> None:
        self.names[self.unit][int(label)] = name.strip()

    def forget_runs(self, entries) -> None:
        names = {e.name for e in entries}
        if names & set(self.run_names):
            self.table = self.result = None
            self.summary, self.blobs_by_run, self.run_names, self.notes = [], {}, [], ""
            self._stale["peaks"] = False
            self.changed.emit()
        if names & set(self.region_runs):
            self.regions = None
            self.region_labels, self.region_shares, self.region_summaries = {}, {}, {}
            self.region_runs, self.region_notes = [], ""
            self._stale["regions"] = False
            self.changed.emit()

    # ----------------------------------------------------------- settings
    def set_settings(self, settings: ClusterSettings, skip_siloxane: bool, scope: str, use_soft: bool = False) -> None:
        if (settings, skip_siloxane, scope, use_soft) != (self.settings, self.skip_siloxane, self.scope, self.use_soft):
            self.settings, self.skip_siloxane, self.scope, self.use_soft = settings, bool(skip_siloxane), scope, bool(use_soft)
            self.mark_stale()
            self.changed.emit()

    def set_unit(self, unit: str, n_regions: int, normalise: str = "percent") -> None:
        """Which result is shown (peaks or pixel regions), how many regions to look for, and how the run table
        is normalised."""
        unit = unit if unit in UNITS else "peaks"
        normalise = normalise if normalise in ("percent", "pqn") else "percent"
        if (unit, n_regions, normalise) != (self.unit, self.n_regions, self.normalise):
            changed_regions = n_regions != self.n_regions
            self.unit, self.n_regions, self.normalise = unit, int(n_regions), normalise
            if changed_regions and self.regions is not None:
                self._stale["regions"] = True
            self.changed.emit()

    def soft_partner(self, index: int):
        """A loaded 14 eV run from the same acquisition as run `index` (found the way 🧬 finds it), or None."""
        return self.project.classes.soft_partner(index)

    def mark_stale(self, *_) -> None:
        self._input_revision += 1
        have = {"peaks": self.result is not None, "regions": self.regions is not None}
        newly = [u for u in UNITS if have[u] and not self._stale[u]]
        for u in newly:
            self._stale[u] = True
        if newly:
            self.changed.emit()

    # ---------------------------------------------------------- clustering
    def _cancelled(self) -> None:
        self.running = False
        self.runningChanged.emit(False)
        self.changed.emit()

    def _job_inputs(self):
        return (self.project.analysis_inputs(), self._input_revision, self.settings, self.skip_siloxane,
                self.scope, self.use_soft, self.unit, self.n_regions,
                self.project.blobs.settings, self.project.blobs.component_settings,
                self.project.blobs.method)

    def cluster(self, index: int, on_done=None) -> bool:
        p = self.project
        reason = self.why_not(index)
        if reason:
            p.message.emit(f"🧩 {reason}")
            return False
        if self.running:
            return False
        if self.unit == "regions":
            return self._segment(index, on_done)
        targets = self.eligible() if self.scope == "all" else [index]
        # Everything the worker needs is read here, on the main thread.
        ms_aware = p.blobs.method == MS_AWARE
        def soft_input(i):
            j = self.soft_partner(i) if self.use_soft else None
            return None if j is None else (p.base_run(j), p.entries[j].background is not None)

        inputs = [(p.entries[i], p.blobs.blobs_for(i), p.base_run(i), p.fold(i), p.entries[i].run,
                   p.blobs.components_for(i) if ms_aware else None, soft_input(i)) for i in targets]
        blob_settings, settings = p.blobs.settings, self.effective_settings()
        n = len(inputs)
        captured = self._job_inputs()
        self.running = True
        self.runningChanged.emit(True)

        def work(token, progress):
            tables, blobs_by_run, found_here, with_soft = [], {}, 0, 0
            for k, (entry, blobs, base_run, fold, raw, component_set, soft_src) in enumerate(inputs):
                progress(k, n + 1, f"Peaks of {entry.name}")
                if blobs is None and component_set is not None and component_set.groups:
                    # fitted overlapping peaks, spectra re-solved on every ion: no detector run
                    blobs, _model, _footprint = respectrum_components(base_run, fold, component_set)
                elif blobs is None:
                    blobs = find_blobs(base_run, fold, blob_settings, floor=storage_threshold(raw) / 3.0)
                    found_here += 1
                token.check()
                soft = None
                if soft_src is not None and len(blobs):
                    src, cleaned = soft_src
                    if not cleaned:
                        progress(k, n + 1, f"Removing background from the 14 eV partner of {entry.name}")
                        src = remove_background(src, Fold(src.times, fold.spec)).run
                    try:
                        soft = soft_spectra_for_blobs(src, fold, blobs, ClassSettings())
                        with_soft += 1
                    except (ValueError, AttributeError):
                        soft = None             # not the same acquisition: use the 70 eV spectrum alone for this run
                blobs_by_run[entry.name] = blobs
                tables.append(PeakTable.from_blobs(blobs, entry.name, min_ions=settings.min_ions, soft=soft))
            progress(n, n + 1, "Clustering")
            table = PeakTable.concat(tables)      # raises ValueError if no run has any peaks
            result = cluster_peaks(table, settings)
            return table, result, blobs_by_run, [e.name for e, *_ in inputs], found_here, with_soft

        def done(out):
            if captured != self._job_inputs():
                self._cancelled()
                p.message.emit("Clustering discarded: inputs changed; cluster again")
                return
            table, result, blobs_by_run, names, found_here, with_soft = out
            self.table, self.result, self.blobs_by_run, self.run_names = table, result, blobs_by_run, names
            self.summary = result.summary(table)
            unassigned = 100.0 * float(table.volume[result.labels < 0].sum()) / (float(table.volume.sum()) or 1.0)
            parts = [f"{len(result)} clusters from {len(table)} peaks in {len(names)} run{'s' * (len(names) != 1)}",
                     f"{unassigned:.1f}% of the volume unassigned"]
            if found_here:
                parts.append(f"legacy detector run here on {found_here} run{'s' * (found_here != 1)} with no 🔍 blobs")
            if self.use_soft:
                parts.append(f"14 eV spectra used ({with_soft} of {len(names)} runs had a partner)" if result.soft_used is not None
                             else "14 eV spectra not used: not every run has a partner run loaded")
            self.notes = "; ".join(parts)
            self.check_text = ""
            self.names["peaks"] = {}                 # a new clustering renumbers everything
            self._stale["peaks"] = False
            self.running = False
            self.runningChanged.emit(False)
            self.changed.emit()
            if on_done:
                on_done(result)

        def failed(msg):
            self.running = False
            self.runningChanged.emit(False)
            p.message.emit(f"🧩 {msg}")
            self.changed.emit()

        p.jobs.submit("🧩 Clustering peaks", work, done, failed, on_cancel=self._cancelled)
        return True

    def _segment(self, index: int, on_done=None) -> bool:
        """Pixel regions: fit on the selected run, or on every 70 eV run at once, and label each."""
        p = self.project
        captured = self._job_inputs()
        targets = self.eligible() if self.scope == "all" else [index]
        pairs = [(p.entries[i].name, p.base_run(i), p.fold(i)) for i in targets]
        settings = SegmentSettings(n_regions=self.n_regions,
                                   exclude_ions=tuple(SILOXANE_IONS) if self.skip_siloxane else ())
        n = len(pairs)
        self.running = True
        self.runningChanged.emit(True)

        def work(token, progress):
            progress(0, n + 1, "Fitting regions")
            model = fit_regions([(run, fold) for _, run, fold in pairs], settings)
            token.check()
            labels, shares, summaries = {}, {}, {}
            for k, (name, run, fold) in enumerate(pairs):
                progress(k + 1, n + 1, f"Labelling {name}")
                labels[name] = model.label_image(run, fold)
                shares[name] = region_shares(labels[name], run, fold, len(model))
                summaries[name] = model.summary(labels[name], run, fold)
                token.check()
            return model, labels, shares, summaries, [name for name, _, _ in pairs]

        def done(out):
            if captured != self._job_inputs():
                self._cancelled()
                p.message.emit("Pixel regions discarded: inputs changed; fit regions again")
                return
            model, labels, shares, summaries, names = out
            self.regions, self.region_labels, self.region_shares, self.region_summaries = model, labels, shares, summaries
            self.region_runs = names
            self.names["regions"] = {}
            self.region_notes = (f"{len(model)} regions fitted on {len(names)} run{'s' * (len(names) != 1)}; "
                                 "regions are found from the pixels' spectra, with no peak detection")
            self._stale["regions"] = False
            self.running = False
            self.runningChanged.emit(False)
            self.changed.emit()
            if on_done:
                on_done(model)

        def failed(msg):
            self.running = False
            self.runningChanged.emit(False)
            p.message.emit(f"🧩 {msg}")
            self.changed.emit()

        p.jobs.submit("🧩 Finding pixel regions", work, done, failed, on_cancel=self._cancelled)
        return True

    # ------------------------------------------------------- reading a result
    def merge(self, labels) -> None:
        """Join the listed clusters into one (for a family the tree split). Raises ValueError for fewer than two."""
        if self.unit != "peaks" or self.result is None:
            raise ValueError("only peak clusters can be merged")
        self.result = merge_clusters(self.table, self.result, labels)
        self.summary = self.result.summary(self.table)
        self.check_text = ""
        self.names["peaks"] = {}                     # merging renumbers by volume
        self.changed.emit()

    def composition_for(self, index: int) -> dict:
        """cluster -> (main compound class, its share of the cluster's peaks in run `index`), from 🧬, or {} if
        that run has not been classified on the same blobs."""
        p = self.project
        if self.unit != "peaks":
            return {}
        per_blob, blobs = self.blob_clusters(index), self.blobs_for(index)
        classes = p.classes.result_for(index)
        used = p.classes.blobs_classified(index)
        if per_blob is None or classes is None or used is None or len(used) != len(blobs) or not np.array_equal(used.t1, blobs.t1):
            return {}
        comp = class_composition(per_blob, [c.cls for c in classes.calls])
        return {k: rows[0] for k, rows in comp.items()}

    def set_group(self, run: str, label: str) -> None:
        label = label.strip()
        if label:
            self.groups[run] = label
        else:
            self.groups.pop(run, None)

    def compare_groups(self) -> dict:
        """Do the user's groups of runs differ in how their peak volume is spread over the clusters? PERMANOVA
        (Bray-Curtis on the percent shares), plus, for exactly two groups, a Welch t-test per cluster with
        Benjamini-Hochberg correction. Raises ValueError with the reason if it cannot be done."""
        fp = self.fingerprint()
        names = self.active_runs()
        if fp is None or len(names) < 4:
            raise ValueError("Pool at least four runs first (scope: every 70 eV run).")
        feats, matrix = fp
        tagged = [(n, self.groups[n]) for n in names if n in self.groups]
        labels = sorted({g for _, g in tagged})
        counts = {g: sum(1 for _, x in tagged if x == g) for g in labels}
        if len(labels) < 2 or min(counts.values()) < 2:
            raise ValueError("Give at least two groups, each of at least two runs, in the Group column.")
        rows = [names.index(n) for n, _ in tagged]
        f, pval = permanova(matrix[rows], [g for _, g in tagged], seed=0)
        points = []
        if len(labels) == 2:
            a, b = ([n for n, g in tagged if g == lab] for lab in labels)
            points = volcano(feats, a, b)
        return dict(f=f, p=pval, groups=labels, counts=counts, points=points)

    # ------------------------------------------------------ how much to trust
    def run_check(self, kind: str, on_done=None) -> bool:
        """`stability` (bootstrap adjusted Rand index) or `scan` (clusters at a range of thresholds), as a job."""
        if self.unit != "peaks" or self.result is None or self.table is None or self.running:
            return False
        table, settings = self.table, self.effective_settings()

        def work(token, progress):
            if kind == "stability":
                s = stability(table, settings)
                return (f"Stability: adjusted Rand index {s['mean']:.2f} ± {s['sd']:.2f} over {s['n_boot']} resamples "
                        "of 80% of the peaks (1 means the same clusters every time).")
            rows = threshold_scan(table, settings, (0.1, 0.2, 0.3, 0.4, 0.5))
            return "Merge distance → clusters, % volume unassigned, mean silhouette:\n" + "\n".join(
                f"  {r['threshold']:.1f} → {r['n_clusters']}, {r['unassigned_percent']:.1f}%, "
                f"{'n/a' if r['silhouette'] != r['silhouette'] else format(r['silhouette'], '.2f')}" for r in rows)

        def done(text):
            self.check_text = text
            self.changed.emit()
            if on_done:
                on_done(text)

        self.project.jobs.submit("🧩 Checking clusters", work, done,
                                 lambda msg: self.project.message.emit(f"🧩 {msg}"))
        return True

    # -------------------------------------------------------------- files
    def to_dict(self) -> dict:
        # Results are not stored: cluster again after reopening.
        return dict(settings=self.settings.to_dict(), skip_siloxane=self.skip_siloxane, scope=self.scope,
                    use_soft=self.use_soft, groups=dict(self.groups), unit=self.unit, n_regions=self.n_regions,
                    normalise=self.normalise)

    def from_dict(self, d: dict) -> None:
        d = d or {}
        self.settings = ClusterSettings.from_dict(d.get("settings", {}))
        self.skip_siloxane = bool(d.get("skip_siloxane", False))
        self.scope = d.get("scope", "run") if d.get("scope") in SCOPES else "run"
        self.use_soft = bool(d.get("use_soft", False))
        self.groups = {str(k): str(v) for k, v in (d.get("groups") or {}).items()}
        self.unit = d.get("unit", "peaks") if d.get("unit") in UNITS else "peaks"
        self.n_regions = int(d.get("n_regions", 12))
        self.normalise = d.get("normalise", "percent") if d.get("normalise") in ("percent", "pqn") else "percent"
        self.names = {"peaks": {}, "regions": {}}
        self.table = self.result = None
        self.summary, self.blobs_by_run, self.run_names, self.notes = [], {}, [], ""
        self.regions, self.region_labels, self.region_shares, self.region_summaries = None, {}, {}, {}
        self.region_runs, self.region_notes = [], ""
        self._stale = {"peaks": False, "regions": False}
        self.check_text = ""
        self.changed.emit()
