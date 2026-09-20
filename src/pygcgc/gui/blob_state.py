"""🔍 Peak-detection state for the app.

Runs either MS-aware candidate/region detection or legacy blob detection on
every loaded run. Both read the 🧹 layer where available, raw otherwise, and
never the 🫧 layer. The selected method and both settings objects are saved;
derived results are rebuilt when a project opens.

Delete, merge and add are kept as a small edit log per run name, replayed
in order on top of a fresh `find_blobs` result whenever blobs are
recomputed (a fold or 🧹 change) -- the same idea as 🧪's review decisions
sitting on top of a fresh quantification. This is a best-effort
mechanism: it's positional (which blob is "index 3" depends on
`find_blobs` producing blobs in the same order), which holds up well
across a fold or background change but not across a change to the 🔍
settings themselves, so changing settings clears the log rather than risk
silently replaying it against the wrong blobs.
"""

from __future__ import annotations

from dataclasses import replace

from PySide6 import QtCore

from ..core.background import storage_threshold
from ..core.candidates import CandidateSettings
from ..core.components import (ComponentDetectionSettings, ComponentGroup, ComponentSet,
                               detect_components)
from ..core.peaks import (BlobSet, PeakSettings, add_blob_at, delete_blobs, find_blobs,
                          merge_blobs)


LEGACY = "legacy"
MS_AWARE = "ms-aware"


class BlobState(QtCore.QObject):
    changed = QtCore.Signal()          # results, staleness or settings changed
    runningChanged = QtCore.Signal(bool)

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.method = LEGACY
        self.settings = PeakSettings()
        self.component_settings = ComponentDetectionSettings()
        self.candidate_settings = self.component_settings.candidates
        self.by_run: dict = {}          # RunEntry -> BlobSet
        self.candidates_by_run: dict = {}
        self.regions_by_run: dict = {}
        self.component_sets_by_name: dict[str, ComponentSet] = {}
        self.accepted_groups_by_name: dict[str, dict[int, ComponentGroup]] = {}
        self.edits: dict[str, list[tuple]] = {}
        self.stale: set = set()         # RunEntry
        self.candidate_stale: set = set()
        self.component_stale_names: set[str] = set()
        self.running = False
        self.was_run = False
        self.candidates_were_run = False

    def forget_runs(self, entries) -> None:
        for e in entries:
            self.by_run.pop(e, None)
            self.candidates_by_run.pop(e, None)
            self.regions_by_run.pop(e, None)
            self.component_sets_by_name.pop(e.name, None)
            self.accepted_groups_by_name.pop(e.name, None)
            self.stale.discard(e)
            self.candidate_stale.discard(e)
            self.component_stale_names.discard(e.name)

    # ------------------------------------------------------------ settings
    def set_settings(self, settings: PeakSettings) -> None:
        self.settings = settings
        self.edits = {}                 # a new detection surface invalidates positional edits
        if self.by_run:
            self.stale = set(self.by_run)
        self.changed.emit()

    def set_method(self, method: str) -> None:
        if method not in (LEGACY, MS_AWARE):
            raise ValueError(method)
        if method != self.method:
            self.method = method
            self.changed.emit()

    def set_candidate_settings(self, settings: CandidateSettings) -> None:
        self.candidate_settings = settings
        self.component_settings = replace(self.component_settings, candidates=settings)
        if self.candidates_by_run:
            self.candidate_stale = set(self.candidates_by_run)
        self.component_stale_names.update(self.component_sets_by_name)
        self.changed.emit()

    def set_component_settings(self, settings: ComponentDetectionSettings) -> None:
        self.component_settings = settings
        self.candidate_settings = settings.candidates
        if self.candidates_by_run:
            self.candidate_stale = set(self.candidates_by_run)
        self.component_stale_names.update(self.component_sets_by_name)
        self.changed.emit()

    def mark_stale(self, *_):
        if (self.by_run or self.candidates_by_run or self.component_sets_by_name) and not self.running:
            self.stale = set(self.by_run)
            self.candidate_stale = set(self.candidates_by_run)
            self.component_stale_names.update(self.component_sets_by_name)
            self.changed.emit()

    def blobs_for(self, index: int) -> BlobSet | None:
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        return self.by_run.get(p.entries[index])

    def candidates_for(self, index: int):
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        return self.candidates_by_run.get(p.entries[index])

    def regions_for(self, index: int):
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        return self.regions_by_run.get(p.entries[index])

    def components_for(self, index: int) -> ComponentSet | None:
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        return self.component_sets_by_name.get(p.entries[index].name)

    def accept_component_group(self, index: int, group: ComponentGroup) -> None:
        p = self.project
        if not 0 <= index < len(p.entries):
            return
        name = p.entries[index].name
        component_set = self.component_sets_by_name.setdefault(name, ComponentSet())
        component_set.replace_group(group)
        self.accepted_groups_by_name.setdefault(name, {})[group.region_index] = group
        self.component_stale_names.discard(name)
        self.changed.emit()

    def selected_stale(self, entry) -> bool:
        return entry in (self.candidate_stale if self.method == MS_AWARE else self.stale)

    # ------------------------------------------------------------ running
    def run_all(self, on_done=None, method: str | None = None) -> bool:
        entries = [self.project.entries[i] for i in self.project.loaded()]
        if not entries or self.running:
            return False
        method = self.method if method is None else method
        if method not in (LEGACY, MS_AWARE):
            raise ValueError(method)
        self.running = True
        self.runningChanged.emit(True)
        remaining = set(entries)
        settings = self.settings
        component_settings = self.component_settings

        for e in entries:
            idx = self.project.entries.index(e)
            run, fold = self.project.base_run(idx), self.project.fold(idx)
            raw_floor = storage_threshold(e.run) / 3.0
            edits = list(self.edits.get(e.name, []))

            def work(token, progress, run=run, fold=fold, edits=edits,
                     raw_floor=raw_floor, run_name=e.name):
                if method == MS_AWARE:
                    def report(done, total):
                        progress(done, total, "Separating overlapping components")
                    return detect_components(run, fold, run_name, component_settings,
                                             floor=raw_floor, progress=report)
                blobs = find_blobs(run, fold, settings, floor=raw_floor)
                return _replay(blobs, run, fold, edits)

            def done(result, e=e):
                if method == MS_AWARE:
                    candidates, regions, automatic = result
                    self.candidates_by_run[e], self.regions_by_run[e] = candidates, regions
                    if component_settings.automatically_separate:
                        for accepted in self.accepted_groups_by_name.get(e.name, {}).values():
                            automatic.replace_group(accepted)
                        self.component_sets_by_name[e.name] = automatic
                        self.component_stale_names.discard(e.name)
                    self.candidate_stale.discard(e)
                else:
                    self.by_run[e] = result
                    self.stale.discard(e)
                remaining.discard(e)
                if not remaining:
                    self.running = False
                    if method == MS_AWARE:
                        self.candidates_were_run = True
                    else:
                        self.was_run = True
                    self.runningChanged.emit(False)
                    self.changed.emit()
                    if on_done:
                        on_done()

            def failed(text, e=e):
                self.project.message.emit(f"🔍 {e.name} failed: {text}")
                remaining.discard(e)
                if not remaining:
                    self.running = False
                    self.runningChanged.emit(False)
                    self.changed.emit()
                    if on_done:
                        on_done()

            label = "MS-aware peaks" if method == MS_AWARE else "legacy blobs"
            self.project.jobs.submit(f"🔍 Finding {label}: {e.name}", work, done, failed)
        return True

    # ------------------------------------------------------------ editing
    def delete(self, index: int, indices) -> None:
        e, blobs = self._entry_and_blobs(index)
        if blobs is None:
            return
        self.by_run[e] = delete_blobs(blobs, indices)
        self.edits.setdefault(e.name, []).append(("delete", tuple(sorted(int(i) for i in indices))))
        self.changed.emit()

    def merge(self, index: int, indices) -> None:
        e, blobs = self._entry_and_blobs(index)
        if blobs is None:
            return
        self.by_run[e] = merge_blobs(blobs, indices)
        self.edits.setdefault(e.name, []).append(("merge", tuple(sorted(int(i) for i in indices))))
        self.changed.emit()

    def add(self, index: int, t1: float, t2: float) -> None:
        e, blobs = self._entry_and_blobs(index)
        if blobs is None:
            return
        run, fold = self.project.base_run(index), self.project.fold(index)
        self.by_run[e] = add_blob_at(blobs, run, fold, t1, t2)
        self.edits.setdefault(e.name, []).append(("add", (float(t1), float(t2))))
        self.changed.emit()

    def _entry_and_blobs(self, index: int):
        p = self.project
        if not 0 <= index < len(p.entries):
            return None, None
        e = p.entries[index]
        return e, self.by_run.get(e)

    # ------------------------------------------------------------------ undo
    def snapshot_run(self, index: int):
        """This run's `BlobSet` and edit log, for `EditBlobs` to restore."""
        e, blobs = self._entry_and_blobs(index)
        if e is None:
            return None
        return (e, blobs, list(self.edits.get(e.name, [])))

    def restore_run(self, snapshot) -> None:
        if snapshot is None:
            return
        e, blobs, edits = snapshot
        if blobs is None:
            self.by_run.pop(e, None)
        else:
            self.by_run[e] = blobs
        if edits:
            self.edits[e.name] = edits
        else:
            self.edits.pop(e.name, None)
        self.changed.emit()

    def snapshot_components(self, index: int):
        if not 0 <= index < len(self.project.entries):
            return None
        name = self.project.entries[index].name
        component_set = self.component_sets_by_name.get(name)
        accepted = self.accepted_groups_by_name.get(name, {})
        return (name,
                None if component_set is None else ComponentSet.from_dict(component_set.to_dict()),
                {region: ComponentGroup.from_dict(group.to_dict())
                 for region, group in accepted.items()})

    def restore_components(self, snapshot) -> None:
        if snapshot is None:
            return
        name, component_set, accepted = snapshot
        if component_set is None:
            self.component_sets_by_name.pop(name, None)
        else:
            self.component_sets_by_name[name] = ComponentSet.from_dict(component_set.to_dict())
        if accepted:
            self.accepted_groups_by_name[name] = {
                region: ComponentGroup.from_dict(group.to_dict())
                for region, group in accepted.items()}
        else:
            self.accepted_groups_by_name.pop(name, None)
        self.changed.emit()

    # -------------------------------------------------------------- files
    def to_dict(self) -> dict:
        return dict(method=self.method, settings=self.settings.to_dict(),
                    candidate_settings=self.candidate_settings.to_dict(),
                    component_settings=self.component_settings.to_dict(),
                    component_sets={name: result.to_dict()
                                    for name, result in self.component_sets_by_name.items()},
                    accepted_groups={name: {str(region): group.to_dict()
                                            for region, group in groups.items()}
                                     for name, groups in self.accepted_groups_by_name.items()},
                    edits={r: [list(op) for op in v] for r, v in self.edits.items()},
                    was_run=self.was_run, candidates_were_run=self.candidates_were_run)

    def from_dict(self, d: dict) -> None:
        d = d or {}
        self.method = d.get("method", LEGACY)
        if self.method not in (LEGACY, MS_AWARE):
            self.method = LEGACY
        self.settings = PeakSettings.from_dict(d.get("settings", {}))
        if d.get("component_settings"):
            self.component_settings = ComponentDetectionSettings.from_dict(
                d.get("component_settings", {}))
        else:
            candidate_settings = CandidateSettings.from_dict(d.get("candidate_settings", {}))
            self.component_settings = ComponentDetectionSettings(candidates=candidate_settings)
        self.candidate_settings = self.component_settings.candidates
        self.edits = {r: [tuple(op) for op in v] for r, v in d.get("edits", {}).items()}
        self.by_run, self.stale = {}, set()
        self.candidates_by_run, self.regions_by_run, self.candidate_stale = {}, {}, set()
        self.component_sets_by_name = {
            name: ComponentSet.from_dict(value)
            for name, value in d.get("component_sets", {}).items()}
        self.accepted_groups_by_name = {
            name: {int(region): ComponentGroup.from_dict(group)
                   for region, group in groups.items()}
            for name, groups in d.get("accepted_groups", {}).items()}
        self.component_stale_names = set()
        self.was_run = bool(d.get("was_run"))
        self.candidates_were_run = bool(d.get("candidates_were_run"))
        self.changed.emit()


def _replay(blobs: BlobSet, run, fold, edits: list[tuple]) -> BlobSet:
    for op, arg in edits:
        try:
            if op == "delete":
                blobs = delete_blobs(blobs, arg)
            elif op == "merge":
                blobs = merge_blobs(blobs, arg)
            elif op == "add":
                blobs = add_blob_at(blobs, run, fold, *arg)
        except (ValueError, IndexError):
            continue   # no longer applies to this run of find_blobs; skip it rather than fail the lot
    return blobs
