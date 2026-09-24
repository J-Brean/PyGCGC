"""The project: which runs are open, how they are folded, and where the
supporting files (sequence list, target table, method) live.

Views never fold anything themselves; they ask the project for a Fold,
so every view agrees on the period and phase of every run.
"""

from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass, replace
from datetime import datetime
from pathlib import Path

import numpy as np
from PySide6 import QtCore


from ..core import (BackgroundResult, BackgroundSettings, Fold, FoldSpec, SiloxaneResult,
                    SiloxaneSettings, check_timing, estimate_period_for_runs, load_anchors,
                    phase_for_sequence, read_run, remove_background, remove_siloxane_blobs,
                    remove_siloxanes)
from ..core.align import blank_mapping
from ..core.phase import Anchor, PhaseEstimate
from ..core.run import Run
from ..core.io import discover_runs
from .jobs import JobManager, Cancelled
from .align_state import AlignState
from .blob_state import BlobState
from .class_state import ClassState
from .cluster_state import ClusterState
from .hydrocarbon_state import HydrocarbonState
from .istd_state import IstdState
from .quantstate import QuantState

PROJECT_VERSION = 1
SIDECAR_NAME = "project.pygcgc.json"     # default project file dropped alongside a folder of runs
AUTOSAVE_DELAY_MS = 1000


@dataclass(eq=False)
class RunEntry:
    path: Path
    name: str
    kind: str = ""
    level: str = ""
    run: Run | None = None
    phase: float = 0.0
    phase_method: str = "unset"      # unset, anchors, void, fallback, manual
    phase_detail: str = ""
    timing: str = ""
    error: str = ""
    import_meta: dict | None = None
    background: BackgroundResult | None = None   # applied 🧹 layer
    siloxane: SiloxaneResult | None = None        # applied 🫧 layer, on top of 🧹 if both are on


class Project(QtCore.QObject):
    runsChanged = QtCore.Signal()                 # entries added, removed or metadata changed
    loadingChanged = QtCore.Signal(bool)
    foldChanged = QtCore.Signal(int)              # run index, or -1 for every run
    anchorsChanged = QtCore.Signal()
    layersChanged = QtCore.Signal()               # a processing layer was added, removed or updated
    message = QtCore.Signal(str)                  # short status-bar messages
    autosaved = QtCore.Signal()                    # processing was written to disk without being asked

    def __init__(self, jobs: JobManager, parent=None):
        super().__init__(parent)
        self.jobs = jobs
        self.entries: list[RunEntry] = []
        self.period: float = 4.0
        self.period_confirmed = False      # set by hand or by the method file; never by an estimate
        self.period_method = "default"     # default, file, manual, auto
        self.period_detail = ""
        self.anchors: list[Anchor] | None = None
        self.targets_path: Path | None = None
        self.sequence_path: Path | None = None
        self.path: Path | None = None
        self.processed_at: str | None = None       # when this project was last saved, for "processed on ..."
        self._data_folder: Path | None = None       # folder opened from, for the default autosave location
        self._loading_project = False                # suppresses autosave while a project is being restored
        self.view_state: dict = {}
        self._folds: dict[tuple, Fold] = {}
        self.loading = False
        # 🧹 background layer
        self.bg_settings = BackgroundSettings()
        self.bg_blank: str | None = None
        self.show_processed = True
        self.preview: tuple[RunEntry, BackgroundResult] | None = None
        self.foldChanged.connect(self._fold_changed)
        # 🫧 siloxane layer
        self.siloxane_settings = SiloxaneSettings()
        # 🧪 quantification
        self.quant = QuantState(self, self)
        self.foldChanged.connect(self.quant.mark_stale)
        self.layersChanged.connect(self._layers_for_quant)
        # 📍 alignment
        self.align = AlignState(self, self)
        self.align.changed.connect(self._alignment_changed_under_blank)
        # 🔍 blob detection
        self.blobs = BlobState(self, self)
        self.foldChanged.connect(self.blobs.mark_stale)
        self.layersChanged.connect(self._layers_for_blobs)
        # ⛓️ hydrocarbon and isomer mapper
        self.hydrocarbons = HydrocarbonState(self, self)
        # 🧬 compound classes (after ⛓️: the ladder gives it a polarity axis)
        self.classes = ClassState(self, self)
        self.foldChanged.connect(self.classes.mark_stale)
        self.blobs.changed.connect(self.classes.mark_stale)
        self.hydrocarbons.changed.connect(self.classes.mark_stale)
        # 🧩 peak clustering (from the same legacy blobs; no ladder needed)
        self.clusters = ClusterState(self, self)
        self.foldChanged.connect(self.clusters.mark_stale)
        self.blobs.changed.connect(self.clusters.mark_stale)
        # 🏷️ internal standards (click-labelled, matched by position and spectrum)
        self.istd = IstdState(self, self)
        # autosave: written silently once a step actually finishes, never mid-preview
        self._autosave_timer = QtCore.QTimer(self)
        self._autosave_timer.setSingleShot(True)
        self._autosave_timer.setInterval(AUTOSAVE_DELAY_MS)
        self._autosave_timer.timeout.connect(self._autosave)
        self.align.changed.connect(self.request_autosave)
        self.istd.changed.connect(self.request_autosave)
        for state in (self.blobs, self.quant, self.classes, self.clusters, self.hydrocarbons):
            state.runningChanged.connect(self._run_finished)

    # ------------------------------------------------------------ queries
    def __len__(self) -> int:
        return len(self.entries)

    def loaded(self) -> list[int]:
        return [i for i, e in enumerate(self.entries) if e.run is not None]

    def fold(self, index: int) -> Fold | None:
        e = self.entries[index]
        if e.run is None:
            return None
        key = (id(e.run), self.period, e.phase)
        f = self._folds.get(key)
        if f is None:
            f = Fold(e.run.times, FoldSpec(period=self.period, phase=e.phase))
            if len(self._folds) > 64:
                self._folds.clear()
            self._folds[key] = f
        return f

    def _layers_for_quant(self) -> None:
        # Quantification uses the applied 🧹 layers; if they changed, results are out of date.
        snap = {e: e.background for e in self.entries}
        if getattr(self, "_last_bg", None) is not None and snap != self._last_bg:
            self.quant.mark_stale()
        self._last_bg = snap

    def _layers_for_blobs(self) -> None:
        # Blobs are found on `base_run`, so only a change beneath 🫧 makes them out of date.
        snap = (self.show_processed, id(self.preview),
                tuple((id(e.run), id(e.background)) for e in self.entries))
        if getattr(self, "_last_blob_layers", None) is not None and snap != self._last_blob_layers:
            self.blobs.mark_stale()
        self._last_blob_layers = snap

    def index_of(self, name: str) -> int | None:
        for i, e in enumerate(self.entries):
            if e.name == name:
                return i
        return None

    # ------------------------------------------------------------ layers
    def any_background(self) -> bool:
        return any(e.background is not None for e in self.entries)

    def any_siloxane(self) -> bool:
        return any(e.siloxane is not None for e in self.entries)

    def result_for(self, index: int) -> BackgroundResult | SiloxaneResult | None:
        """The processed result a view should use for this run, if any:
        🫧 (built on whatever was displayed when it ran) if applied,
        otherwise 🧹, otherwise raw."""
        e = self.entries[index]
        if self.preview is not None and self.preview[0] is e:
            return self.preview[1]
        if self.show_processed and e.siloxane is not None:
            return e.siloxane
        if self.show_processed and e.background is not None:
            return e.background
        return None

    def display_run(self, index: int) -> Run | None:
        e = self.entries[index]
        res = self.result_for(index)
        return res.run if res is not None else e.run

    def base_run(self, index: int) -> Run | None:
        """What is displayed, minus any 🫧 layer: 🧹 if shown, otherwise raw.
        🔍 detects blobs on this, so a blob blocked out by 🫧 is still found
        (and can be scored, kept or un-removed), and 🫧 is applied to this, so
        applying it again replaces the earlier layer instead of stacking on it."""
        e = self.entries[index]
        if self.preview is not None and self.preview[0] is e:
            return self.preview[1].run
        if self.show_processed and e.background is not None:
            return e.background.run
        return e.run

    def set_show_processed(self, on: bool) -> None:
        if on != self.show_processed:
            self.show_processed = bool(on)
            self.layersChanged.emit()

    def set_preview(self, entry: RunEntry | None, result: BackgroundResult | None) -> None:
        self.preview = None if entry is None or result is None else (entry, result)
        self.layersChanged.emit()

    def background_snapshot(self) -> dict:
        return {e: e.background for e in self.entries}

    def set_backgrounds(self, results: dict, settings: BackgroundSettings | None = None,
                        blank: str | None = "__keep__") -> None:
        for e, res in results.items():
            if e in self.entries:
                e.background = res
                e.siloxane = None     # was fit against the old 🧹 layer; re-apply 🫧 if still wanted
        if settings is not None:
            self.bg_settings = settings
        if blank != "__keep__":
            self.bg_blank = blank
        self.layersChanged.emit()
        self.runsChanged.emit()
        self.request_autosave()

    def _alignment_changed_under_blank(self) -> None:
        """Blank subtraction reads the blank through the runs' alignment, so
        an applied blank subtraction is out of date once an alignment changes.
        Nothing is recomputed automatically (a background job per run per
        anchor click would be far too heavy); the message says what to do."""
        if self.bg_blank and self.bg_settings.subtract_blank and self.any_background():
            self.message.emit("Alignment changed: apply 🧹 again to redo the blank subtraction with it")

    def compute_background(self, indices: list[int], settings: BackgroundSettings,
                           blank: str | None, on_done=None, preview: bool = False) -> None:
        """Run 🧹 on the given runs in background jobs. If a blank is to be
        subtracted, it is cleaned first (without subtracting itself) and the
        other runs start once it is ready. `on_done(results)` receives a dict
        of entry -> result when every run has finished."""
        entries = [self.entries[i] for i in indices if self.entries[i].run is not None]
        if not entries:
            return
        blank_i = self.index_of(blank) if (settings.subtract_blank and blank) else None
        blank_entry = self.entries[blank_i] if blank_i is not None else None
        if blank_entry is not None and blank_entry.run is None:
            blank_entry = None
        results: dict = {}
        remaining = set(entries)
        latest = "bg-preview" if preview else None

        def start(entry, blank_res, blank_fold):
            idx = self.entries.index(entry)
            run, fold = entry.run, self.fold(idx)
            s = settings
            use_blank = blank_res is not None and entry is not blank_entry
            if not use_blank:
                s = replace(settings, subtract_blank=False)
            # Blank subtraction follows the alignment where there is one: each
            # run's transform to the reference (none for the reference itself).
            t_sample = self.align.transforms.get(entry.name)
            t_blank = self.align.transforms.get(blank) if blank else None
            to_blank = blank_mapping(t_sample, t_blank) if (t_sample or t_blank) else None

            def work(token, progress):
                return remove_background(run, fold, s,
                                         blank=blank_res if use_blank else None,
                                         blank_fold=blank_fold if use_blank else None,
                                         blank_name=blank if use_blank else None,
                                         progress=progress, to_blank=to_blank if use_blank else None)

            def done(res):
                results[entry] = res
                remaining.discard(entry)
                if not remaining and on_done:
                    on_done(results)

            name = f"🧹 {'Preview' if preview else 'Background'}: {entry.name}"
            self.jobs.submit(name, work, done, lambda t: self.message.emit(f"{name} failed: {t}"),
                             latest_key=latest, quiet=preview)

        if blank_entry is None:
            for e in entries:
                start(e, None, None)
            return

        bi = self.entries.index(blank_entry)
        b_run, b_fold = blank_entry.run, self.fold(bi)
        b_settings = replace(settings, subtract_blank=False)

        def blank_work(token, progress):
            return remove_background(b_run, b_fold, b_settings, progress=progress)

        def blank_done(bres):
            if blank_entry in remaining:
                results[blank_entry] = bres
                remaining.discard(blank_entry)
            for e in entries:
                if e is not blank_entry:
                    start(e, bres, b_fold)
            if not remaining and on_done:
                on_done(results)

        self.jobs.submit(f"🧹 Cleaning blank {blank_entry.name}", blank_work, blank_done,
                         lambda t: self.message.emit(f"Blank cleaning failed: {t}"),
                         latest_key=("bg-preview-blank" if preview else None), quiet=preview)

    # ------------------------------------------------------------ 🫧 layer
    def siloxane_snapshot(self) -> dict:
        return {e: e.siloxane for e in self.entries}

    def set_siloxanes(self, results: dict, settings: SiloxaneSettings | None = None) -> None:
        for e, res in results.items():
            if e in self.entries:
                e.siloxane = res
        if settings is not None:
            self.siloxane_settings = settings
        self.layersChanged.emit()
        self.runsChanged.emit()
        self.request_autosave()

    def compute_siloxane(self, indices: list[int], settings: SiloxaneSettings, on_done=None) -> None:
        """Run 🫧 on the given runs (their currently displayed data, 🧹
        cleaned if that's been applied) in background jobs. `on_done(results)`
        receives a dict of entry -> result when every run has finished. Not
        wired into 🧪's own measurement, which still reads `entry.background`
        directly -- 🫧 affects the 2D/1D view and 🔍 blob detection, both of
        which go through `display_run`/`result_for`."""
        entries = [self.entries[i] for i in indices if self.entries[i].run is not None]
        by_blobs = settings.method == "blobs"
        if by_blobs:
            # The blob method blocks out the 🔍 blobs, so a run needs current ones.
            missing = [e for e in entries if self.blobs.by_run.get(e) is None or e in self.blobs.stale]
            if missing:
                self.message.emit("🫧 Find blobs first: " + ", ".join(e.name for e in missing))
            entries = [e for e in entries if e not in missing]
        if not entries:
            if on_done:
                on_done({})
            return
        results: dict = {}
        remaining = set(entries)
        for e in entries:
            idx = self.entries.index(e)
            run, fold = self.base_run(idx), self.fold(idx)
            blobs = self.blobs.by_run.get(e) if by_blobs else None

            def work(token, progress, run=run, fold=fold, blobs=blobs):
                if blobs is not None:
                    return remove_siloxane_blobs(run, fold, blobs, settings)
                return remove_siloxanes(run, fold, settings)

            def done(res, e=e):
                results[e] = res
                remaining.discard(e)
                if not remaining and on_done:
                    on_done(results)

            name = f"🫧 Siloxanes: {e.name}"
            self.jobs.submit(name, work, done, lambda t, name=name: self.message.emit(f"{name} failed: {t}"))

    def _fold_changed(self, index: int) -> None:
        """A phase or period change moves the fold, so applied background
        layers are recomputed with their own settings. A 🫧 layer, built on
        top of one of them, isn't recomputed the same way (it would mean
        chaining another background job after this one); it's simply
        cleared, so a stale fit is never shown silently -- re-apply 🫧
        after moving the fold."""
        if index == -1:
            affected = [i for i, e in enumerate(self.entries) if e.background is not None]
            siloxane_affected = [i for i, e in enumerate(self.entries) if e.siloxane is not None]
        else:
            affected = [index] if self.entries[index].background is not None else []
            blank_i = self.index_of(self.bg_blank) if self.bg_blank else None
            if blank_i == index:
                affected = [i for i, e in enumerate(self.entries) if e.background is not None]
            siloxane_affected = [index] if self.entries[index].siloxane is not None else []
        for i in siloxane_affected:
            self.entries[i].siloxane = None
        if siloxane_affected:
            self.layersChanged.emit()
        if not affected:
            return
        settings = self.entries[affected[0]].background.settings
        blank = self.bg_blank

        def done(results):
            for e, res in results.items():
                if e.background is not None:     # skip runs whose layer was removed meanwhile
                    e.background = res
            self.layersChanged.emit()
            self.request_autosave()

        self.compute_background(affected, settings, blank, on_done=done)

    # ------------------------------------------------------------ editing
    def set_phase(self, index: int, phase: float, method: str = "manual", detail: str = "") -> None:
        e = self.entries[index]
        e.phase = float(phase) % self.period
        e.phase_method = method
        e.phase_detail = detail
        self.foldChanged.emit(index)
        self.runsChanged.emit()

    def set_phases(self, phases: dict[int, tuple[float, str, str]]) -> None:
        for i, (p, method, detail) in phases.items():
            e = self.entries[i]
            e.phase, e.phase_method, e.phase_detail = float(p) % self.period, method, detail
        self.foldChanged.emit(-1)
        self.runsChanged.emit()

    def set_period(self, period: float, method: str = "manual", detail: str = "") -> None:
        if period <= 0:
            raise ValueError("period must be positive")
        self.period_confirmed = method != "auto"
        self.period_method, self.period_detail = method, detail
        self.period = float(period)
        for e in self.entries:
            e.phase %= self.period
        self.foldChanged.emit(-1)
        self.runsChanged.emit()

    def phase_snapshot(self) -> dict[int, tuple[float, str, str]]:
        return {i: (e.phase, e.phase_method, e.phase_detail) for i, e in enumerate(self.entries)}

    def fold_snapshot(self) -> dict:
        return {"period": self.period, "confirmed": self.period_confirmed, "method": self.period_method,
                "detail": self.period_detail, "phases": self.phase_snapshot()}

    def restore_fold(self, snap: dict) -> None:
        self.period, self.period_confirmed = snap["period"], snap["confirmed"]
        self.period_method, self.period_detail = snap["method"], snap["detail"]
        self.set_phases(snap["phases"])

    # ------------------------------------------------------------ loading
    def add_files(self, paths, on_loaded=None) -> None:
        self._loading_project = False      # a fresh user action; any earlier stuck load guard is over
        paths = [Path(p) for p in paths]
        known = {e.path.resolve() for e in self.entries}
        new = [RunEntry(path=p, name=p.stem) for p in paths if p.resolve() not in known]
        if not new:
            return
        self._apply_sequence(new)
        self.entries.extend(new)
        self.runsChanged.emit()
        self._read(new, on_loaded)

    def _read(self, entries: list[RunEntry], on_loaded=None, estimate: bool = True) -> None:
        self.loading = True
        self.loadingChanged.emit(True)

        def work(token, progress):
            out = []
            for k, e in enumerate(entries):
                progress(k, len(entries), f"Reading {e.path.name}")
                try:
                    run = read_run(e.path, check_cancel=token.check,
                                   progress=lambda done, total, message: progress(
                                       k * 1000 + int(1000 * done / max(total, 1)),
                                       len(entries) * 1000, message))
                    if e.import_meta and e.import_meta.get('source_sha256') != run.meta.get('source_sha256'):
                        raise ValueError('Source file changed since this project was saved; reopen it as a new run')
                    out.append((e, run, check_timing(run.times), ""))
                except Cancelled:
                    raise
                except Exception as exc:  # noqa: BLE001
                    out.append((e, None, None, str(exc)))
            progress(len(entries), len(entries), "Read")
            return out

        def done(results):
            for e, run, rep, err in results:
                e.run, e.error = run, err
                if run is not None and run.meta.get('format') == 'Markes LSC':
                    e.import_meta = {key: run.meta[key] for key in (
                        'source_sha256', 'parser_version', 'acquisition_id', 'energy_ev',
                        'processing_state', 'time_origin_seconds', 'time_origin_policy')}
                    e.phase_detail = 'Native LSC: verify modulation period and vendor retention-time origin.'
                e.timing = "" if rep is None else ("OK" if rep.ok else "Check timing")
            self.loading = False
            self.loadingChanged.emit(False)
            self.runsChanged.emit()
            self.foldChanged.emit(-1)
            bad = [e.name for e, _, _, err in results if err]
            if bad:
                self.message.emit(f"Could not read: {', '.join(bad)}")
            if estimate:
                new = [i for i, x in enumerate(self.entries) if x in entries and x.run is not None]
                self.auto_fold(apply_to=new)
            if on_loaded:
                on_loaded()

        def failed(text):
            self.loading = False
            self.loadingChanged.emit(False)
            self._loading_project = False      # a project load that dies here would otherwise wedge autosave off
            self.message.emit(f"Loading failed: {text}")

        n = len(entries)
        def cancelled():
            for entry in entries:
                if entry.run is None:
                    entry.error = 'Import cancelled; remove and reopen to retry'
            self.loading = False
            self.loadingChanged.emit(False)
            self._loading_project = False      # a project load that dies here would otherwise wedge autosave off
            self.runsChanged.emit()

        self.jobs.submit(f"📂 Reading {n} run{'s' if n > 1 else ''}", work, done, failed,
                         on_cancel=cancelled)

    def open_folder(self, folder, on_loaded=None) -> None:
        """Discover acquisition alternatives off the GUI thread, then load."""
        self._loading_project = False      # a fresh user action; any earlier stuck load guard is over
        self._data_folder = Path(folder)
        self.loading = True
        self.loadingChanged.emit(True)

        def work(token, progress):
            progress(0, 1, 'Discovering native channels and run files')
            return discover_runs(folder, check_cancel=token.check)

        def finish():
            self.loading = False
            self.loadingChanged.emit(False)

        def done(files):
            finish()
            self._open_folder_selected(folder, files, on_loaded)

        def failed(message):
            finish()
            self.message.emit(message)

        self.jobs.submit('Discovering runs', work, done, failed, on_cancel=finish)

    def _open_folder_selected(self, folder, files, on_loaded=None) -> None:
        """Open a folder of runs. Picks up sequence.csv, target_table.csv and
        method.json next to the data if they exist (the dummy pack layout)."""
        folder = Path(folder)
        if (folder / 'sequence.csv').exists():
            with (folder / 'sequence.csv').open(newline='', encoding='utf-8-sig') as fh:
                requested = [row['file'] for row in csv.DictReader(fh)]
            files = []
            for name in requested:
                matches = list(folder.rglob(name))
                if len(matches) != 1:
                    self.message.emit(f'Expected one file matching {name}; found {len(matches)}')
                    return
                files.append(matches[0])
        if not files:
            self.message.emit(f"No supported run files in {folder}")
            return
        method = folder / "method.json"
        if method.exists():
            try:
                self.period = float(json.loads(method.read_text())["modulation_period_s"])
                self.period_confirmed = True
                self.period_method, self.period_detail = "file", f"From {method.name}"
            except (KeyError, ValueError):
                pass
        if (folder / "sequence.csv").exists():
            self.load_sequence(folder / "sequence.csv")
        if (folder / "target_table.csv").exists():
            self.load_targets(folder / "target_table.csv", estimate=False)
        for name in ("calibration_levels.csv", "calibration_levels_ng.csv"):
            if (folder / name).exists():
                self.quant.set_levels(folder / name)
                break
        self.add_files(files, on_loaded)

    def load_sequence(self, path) -> None:
        self.sequence_path = Path(path)
        self.quant.set_sequence(path)
        self._apply_sequence(self.entries)
        self.runsChanged.emit()

    def _apply_sequence(self, entries) -> None:
        if not self.sequence_path or not self.sequence_path.exists():
            return
        info = {}
        with open(self.sequence_path, newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                info[Path(row.get("file", "")).stem] = row
        for e in entries:
            row = info.get(e.name)
            if row:
                e.kind = row.get("type", "")
                e.level = row.get("level", "")

    def load_targets(self, path, estimate: bool = True) -> None:
        self.targets_path = Path(path)
        self.anchors = load_anchors(path)
        self.quant.set_targets(path)
        self.anchorsChanged.emit()
        self.message.emit(f"{len(self.anchors)} anchors from {self.targets_path.name}")
        if estimate and self.loaded():
            self.estimate_phases()

    def estimate_phases(self, on_done=None, apply_to: list[int] | None = None) -> None:
        """Estimate phases using every loaded run (so runs without anchors can
        borrow), but only change the runs in `apply_to` (default: all)."""
        idx = self.loaded()
        if not idx:
            return
        ents = [self.entries[i] for i in idx]
        chosen = ents if apply_to is None else [self.entries[i] for i in apply_to]
        runs = [e.run for e in ents]
        period, anchors = self.period, self.anchors

        def work(token, progress):
            progress(0, 1, "Estimating phase")
            return phase_for_sequence(runs, period, anchors)

        def done(estimates: list[PhaseEstimate]):
            before = self.phase_snapshot()
            # Runs may have been removed meanwhile, so match by identity.
            self.set_phases({self.entries.index(e): (est.phase, est.method, est.summary())
                             for e, est in zip(ents, estimates)
                             if e in chosen and e in self.entries})
            if on_done:
                on_done(before, self.phase_snapshot())

        self.jobs.submit("🔄 Estimating phase", work, done)

    def auto_fold(self, apply_to: list[int] | None = None, on_done=None, force_period: bool = False) -> None:
        """Estimate the modulation period, then the phase of every run.

        A period set by hand or by the method file is kept (the run is only
        checked against it), unless `force_period`. An estimated period is
        applied at once and flagged, but only if a period is really
        established; otherwise the current one stays. Phases are then
        estimated at the period in use, as in `estimate_phases`."""
        idx = self.loaded()
        if not idx:
            return
        ents = [self.entries[i] for i in idx]
        chosen = ents if apply_to is None else [self.entries[i] for i in apply_to]
        runs = [e.run for e in ents]
        period, anchors = self.period, self.anchors
        keep = self.period_confirmed and not force_period

        def work(token, progress):
            progress(0, 2, "Estimating the modulation period")
            best, found = estimate_period_for_runs(runs)   # blind, even when checking a confirmed period
            use = period if keep or best is None else best.period
            progress(1, 2, "Estimating phase")
            return best, found, use, phase_for_sequence(runs, use, anchors)

        def done(result):
            best, found, use, estimates = result
            if self.period != period or any(e not in self.entries for e in ents):
                self.message.emit("Runs or period changed while estimating; estimate discarded")
                return
            before = self.fold_snapshot()
            apply = chosen
            if keep:
                se = best.correction_se if best is not None and np.isfinite(best.correction_se) else 0.0
                if best is not None and abs(best.period - period) > max(3 * se, 0.001):
                    self.message.emit(f"The run suggests a period of {best.period:.4f} s, not the "
                                      f"{period:.4f} s set for it. Use 'Estimate period' to apply it.")
            elif best is not None:
                if use != period:
                    apply = ents      # phases found at the old period no longer mean anything
                self.period = use
                self.period_confirmed = False
                self.period_method, self.period_detail = "auto", best.summary()
                self.message.emit(f"Modulation period estimated: {best.summary()}")
            else:
                guess = f" (best guess {found[0].period:.4f} s, low confidence)" if found else ""
                self.period_detail = f"Could not be estimated from the run{guess}; the default is unchanged"
                self.message.emit(f"Modulation period {self.period_detail}")
            self.set_phases({self.entries.index(e): (est.phase, est.method, est.summary())
                             for e, est in zip(ents, estimates) if e in apply})
            self.runsChanged.emit()
            if on_done:
                on_done(before, self.fold_snapshot())

        self.jobs.submit("🔄 Estimating period and phase", work, done)

    def remove(self, indices) -> None:
        removed = [self.entries[i] for i in set(indices)]
        self.quant.forget_runs(removed)
        self.blobs.forget_runs(removed)
        self.hydrocarbons.forget_runs(removed)
        self.classes.forget_runs(removed)
        self.clusters.forget_runs(removed)
        for i in sorted(set(indices), reverse=True):
            if self.preview is not None and self.preview[0] is self.entries[i]:
                self.preview = None
            if self.bg_blank == self.entries[i].name:
                self.bg_blank = None
            del self.entries[i]
        self._folds.clear()
        self.runsChanged.emit()
        self.foldChanged.emit(-1)

    # ------------------------------------------------------------ autosave
    def _run_finished(self, running: bool) -> None:
        if not running:
            self.request_autosave()

    def request_autosave(self) -> None:
        """Write the project soon, silently. A job finishing and the layer it
        produces arrive together, so requests are coalesced into one write."""
        if not self._loading_project:
            self._autosave_timer.start()

    def flush_autosave(self) -> None:
        if self._autosave_timer.isActive():
            self._autosave_timer.stop()
            self._autosave()

    def _default_autosave_path(self) -> Path | None:
        folder = self._data_folder or (self.entries[0].path.parent if self.entries else None)
        return folder / SIDECAR_NAME if folder is not None else None

    def _autosave(self) -> None:
        # Skipped while a project is being restored: that would overwrite the
        # file being read and reset `processed_at` to now.
        if self._loading_project or not self.entries:
            return
        path = self.path or self._default_autosave_path()
        if path is None:
            return
        try:
            self.save(path, autosave=True)
        except OSError as exc:
            self.message.emit(f"Autosave to {path} failed: {exc}")

    # ------------------------------------------------------------ files
    def to_dict(self, base: Path | None = None) -> dict:
        def rel(p):
            if p is None:
                return None
            p = Path(p)
            if base is not None:
                try:
                    return str(p.resolve().relative_to(base.resolve()))
                except ValueError:
                    pass
            return str(p)

        return {
            "pygcgc_project": PROJECT_VERSION,
            "processed_at": self.processed_at,
            "period": self.period,
            "period_confirmed": self.period_confirmed,
            "period_method": self.period_method,
            "period_detail": self.period_detail,
            "sequence": rel(self.sequence_path),
            "targets": rel(self.targets_path),
            "runs": [{"path": rel(e.path), "phase": e.phase, "phase_method": e.phase_method,
                      "import_meta": e.import_meta,
                      "background": e.background is not None, "siloxane": e.siloxane is not None}
                     for e in self.entries],
            "background": {"settings": self.bg_settings.to_dict(), "blank": self.bg_blank,
                           "show_processed": self.show_processed},
            "siloxane": {"settings": self.siloxane_settings.to_dict()},
            "quant": self.quant.to_dict(rel),
            "align": self.align.to_dict(),
            "blobs": self.blobs.to_dict(),
            "hydrocarbons": self.hydrocarbons.to_dict(),
            "classes": self.classes.to_dict(),
            "clusters": self.clusters.to_dict(),
            "istd": self.istd.to_dict(),
            "view": self.view_state,
        }

    def save(self, path, autosave: bool = False) -> None:
        path = Path(path)
        self._loading_project = False      # an explicit save is proof loading (if any was stuck) has ended
        previous = self.processed_at
        self.processed_at = datetime.now().isoformat(timespec="seconds")
        # Written beside the target, then swapped in, so a crash mid-write
        # never leaves the only copy of the processing half-written.
        tmp = path.with_name(path.name + ".tmp")
        try:
            tmp.write_text(json.dumps(self.to_dict(path.parent), indent=2))
            os.replace(tmp, path)
        except OSError:
            self.processed_at = previous
            raise
        self.path = path
        if autosave:
            self.autosaved.emit()
        else:
            self.message.emit(f"Saved {path.name}")

    def load(self, path, on_loaded=None) -> None:
        path = Path(path)
        data = json.loads(path.read_text())
        if data.get("pygcgc_project") != PROJECT_VERSION:
            raise ValueError("not a PyGCGC project file, or a newer version")
        base = path.parent
        self._loading_project = True      # restoring settings triggers the same calls processing does
        self._autosave_timer.stop()        # a write still pending belongs to the project being replaced

        def finished():
            self._loading_project = False
            if on_loaded:
                on_loaded()

        def resolve(p):
            if p is None:
                return None
            p = Path(p)
            return p if p.is_absolute() else base / p

        self.entries.clear()
        self._folds.clear()
        self.path = path
        self.processed_at = data.get("processed_at")
        self.period = float(data["period"])
        self.period_confirmed = bool(data.get('period_confirmed', False))
        self.period_method = data.get("period_method", "manual" if self.period_confirmed else "default")
        self.period_detail = data.get("period_detail", "")
        self.view_state = data.get("view", {})
        self.sequence_path = resolve(data.get("sequence"))
        self.anchors, self.targets_path = None, None
        if data.get("targets") and resolve(data["targets"]).exists():
            self.targets_path = resolve(data["targets"])
            self.anchors = load_anchors(self.targets_path)
            self.quant.set_targets(self.targets_path)
            self.anchorsChanged.emit()
        if self.sequence_path is not None and self.sequence_path.exists():
            self.quant.set_sequence(self.sequence_path)
        self.quant.from_dict(data.get("quant", {}), resolve)
        self.align.from_dict(data.get("align", {}))
        self.blobs.from_dict(data.get("blobs", {}))
        self.hydrocarbons.from_dict(data.get("hydrocarbons", {}))
        self.classes.from_dict(data.get("classes", {}))
        self.clusters.from_dict(data.get("clusters", {}))
        self.istd.from_dict(data.get("istd", {}))
        bg =data.get("background", {})
        self.bg_settings = BackgroundSettings.from_dict(bg.get("settings", {}))
        self.bg_blank = bg.get("blank")
        self.show_processed = bool(bg.get("show_processed", True))
        self.siloxane_settings = SiloxaneSettings.from_dict(data.get("siloxane", {}).get("settings", {}))
        self.preview = None
        entries, wants_bg, wants_sil = [], [], []
        for r in data["runs"]:
            p = resolve(r["path"])
            e = RunEntry(path=p, name=p.stem, phase=float(r.get("phase", 0.0)),
                         phase_method=r.get("phase_method", "unset"), import_meta=r.get('import_meta'))
            entries.append(e)
            wants_bg.append(bool(r.get("background")))
            wants_sil.append(bool(r.get("siloxane")))
        self._apply_sequence(entries)
        self.entries = entries
        self.runsChanged.emit()
        self.layersChanged.emit()

        def loaded():
            # Saved phases are kept; background layers are rebuilt from their settings.
            todo = [i for i, w in enumerate(wants_bg) if w and self.entries[i].run is not None]
            if todo:
                def done(results):
                    self.set_backgrounds(results)
                    resiloxane()
                self.compute_background(todo, self.bg_settings, self.bg_blank, on_done=done)
            else:
                resiloxane()

        def resiloxane():
            # 🫧 isn't stored either; rebuilt on top of whatever 🧹 just restored.
            todo = [i for i, w in enumerate(wants_sil) if w and self.entries[i].run is not None]
            if todo:
                def done(results):
                    self.set_siloxanes(results)
                    requant()

                def go():
                    self.compute_siloxane(todo, self.siloxane_settings, on_done=done)

                # The blob method needs the 🔍 blobs, found on the layer beneath 🫧.
                if self.siloxane_settings.method == "blobs":
                    if not self.blobs.run_all(on_done=go, method="legacy"):
                        requant()
                else:
                    go()
            else:
                requant()

        def reblob():
            # Blobs aren't stored either; rebuilt from the saved settings and edit log
            # (unless 🫧 already needed them and found them just now).
            if self.blobs.method == "ms-aware":
                fresh = all(self.blobs.candidates_by_run.get(e) is not None
                            and e not in self.blobs.candidate_stale
                            for e in self.entries if e.run is not None)
                wanted = self.blobs.candidates_were_run
            else:
                fresh = all(self.blobs.by_run.get(e) is not None and e not in self.blobs.stale
                            for e in self.entries if e.run is not None)
                wanted = self.blobs.was_run
            if wanted and self.loaded() and not fresh:
                self.blobs.run_all(on_done=finished)
            else:
                finished()

        def requant():
            # Results aren't stored; they are rebuilt from the saved settings and decisions.
            if self.quant.was_run and self.quant.ready()[0]:
                self.quant.run(on_done=reblob)
            else:
                reblob()

        self._read(entries, loaded, estimate=False)
