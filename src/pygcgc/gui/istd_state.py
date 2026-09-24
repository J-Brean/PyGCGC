"""🏷️ Internal-standard state for the app.

Holds the IS library, the peak labelled as each IS on a reference run, and
the per-run match results found by searching for them elsewhere. The
matching itself is all in `core/istd.py`; this just wires it to the
project (for runs, folds and peaks) and to background jobs, the same shape
as 📍's `AlignState`.

Peaks: like 🧬's `ClassState`, this wants full spectra. If the legacy
detector's blobs are there (`blobs_for`), they're used directly; otherwise,
for the MS-aware default, a run's components are re-solved on every ion via
`core.peaktable.respectrum_components`. Nothing is cached across calls --
each label or propagate re-resolves the one run it needs, inside its own
background job, the same cost model `ClassState.classify` already accepts.
"""

from __future__ import annotations

import numpy as np
from PySide6 import QtCore

from ..core.istd import (IstdLibraryEntry, IstdMatch, IstdMatchSettings, IstdReference,
                         load_library, match_run, peak_spectrum, save_library)
from ..core.peaktable import respectrum_components
from .blob_state import MS_AWARE


class IstdState(QtCore.QObject):
    changed = QtCore.Signal()   # library, a reference, or match results changed

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.library: list[IstdLibraryEntry] = load_library()
        self.reference_run: str | None = None
        self.references: dict[str, IstdReference] = {}
        self.matches: dict[str, dict[str, IstdMatch]] = {}
        self.settings = IstdMatchSettings()
        project.runsChanged.connect(self._prune)

    def _prune(self) -> None:
        names = {e.name for e in self.project.entries}
        if self.reference_run not in names:
            self.reference_run = None
        for name in [n for n in self.matches if n not in names]:
            self.matches.pop(name, None)

    # ------------------------------------------------------------- library
    def add_custom(self, name: str, cas: str = "", product_code: str = "") -> None:
        self.library = self.library + [IstdLibraryEntry(name=name, cas=cas, product_code=product_code, custom=True)]
        save_library(self.library)
        self.changed.emit()

    def remove_from_library(self, name: str) -> None:
        self.library = [e for e in self.library if e.name != name]
        save_library(self.library)
        self.changed.emit()

    # ------------------------------------------------------------- peaks
    def _peaks_for(self, index: int):
        """This run's peaks with full spectra, however they were found --
        `None` if peak detection hasn't been run on it at all. Only ever
        called from inside a background job's `work()`."""
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        blobs = p.blobs.blobs_for(index)
        if blobs is not None:
            return blobs
        if p.blobs.method == MS_AWARE:
            component_set = p.blobs.components_for(index)
            if component_set is not None and component_set.groups:
                base_run, fold = p.base_run(index), p.fold(index)
                peaks, _model, _footprint = respectrum_components(base_run, fold, component_set)
                return peaks
        return None

    @staticmethod
    def _resolve_index(peaks, peak_ref) -> int | None:
        """`peak_ref` is a stable component id (MS-aware) or a positional
        blob index (legacy, the same positional trade-off `BlobState`
        itself already accepts for its edit log)."""
        if peaks is None or len(peaks) == 0:
            return None
        if isinstance(peak_ref, str):
            ids = getattr(peaks, "component_id", None)
            if not ids or peak_ref not in ids:
                return None
            return ids.index(peak_ref)
        i = int(peak_ref)
        return i if 0 <= i < len(peaks) else None

    # ------------------------------------------------------------- labelling
    def label(self, name: str, index: int, peak_ref, on_done=None) -> None:
        p = self.project
        if not 0 <= index < len(p.entries):
            return
        entry = p.entries[index]
        cas = next((e.cas for e in self.library if e.name == name), "")

        def work(token, progress):
            return self._peaks_for(index)

        def done(peaks) -> None:
            i = self._resolve_index(peaks, peak_ref)
            if i is None:
                p.message.emit(f"🏷️ Couldn't read that peak on {entry.name} -- run peak detection first")
                return
            self.references[name] = IstdReference(
                name=name, cas=cas, run=entry.name, t1=float(peaks.t1[i]), t2=float(peaks.t2[i]),
                mz0=int(peaks.mz[0]), spectrum=peak_spectrum(peaks, i))
            if self.reference_run is None:
                self.reference_run = entry.name
            self.matches.pop(entry.name, None)   # a fresh label invalidates any match already found on its own run
            self.changed.emit()
            if on_done:
                on_done()

        p.jobs.submit(f"🏷️ Labelling {name} on {entry.name}", work, done)

    def remove_reference(self, name: str) -> None:
        self.references.pop(name, None)
        for per_run in self.matches.values():
            per_run.pop(name, None)
        self.changed.emit()

    # ------------------------------------------------------------- batch
    def propagate_to(self, indices: list[int], on_done=None) -> None:
        if not self.references:
            self.project.message.emit("🏷️ Label at least one internal standard first")
            return
        p = self.project
        entries = [p.entries[i] for i in indices
                  if p.entries[i].run is not None and p.entries[i].name != self.reference_run]
        if not entries:
            return
        remaining = set(entries)
        not_found: list[str] = []
        references = list(self.references.values())
        settings = self.settings

        def start(entry) -> None:
            idx = p.entries.index(entry)

            def work(token, progress):
                peaks = self._peaks_for(idx)
                if peaks is None:
                    return None
                return match_run(references, entry.name, peaks, settings)

            def done(results) -> None:
                if results is None:
                    not_found.append(f"{entry.name} (no peaks detected)")
                else:
                    self.matches[entry.name] = {m.name: m for m in results}
                    missing = [m.name for m in results if m.tier == "not found"]
                    if missing:
                        not_found.append(f"{entry.name} ({', '.join(missing)})")
                remaining.discard(entry)
                if not remaining:
                    self.changed.emit()
                    if not_found:
                        self.project.message.emit(
                            "🏷️ Couldn't find every internal standard in: " + "; ".join(not_found))
                    if on_done:
                        on_done()

            p.jobs.submit(f"🏷️ Finding internal standards: {entry.name}", work, done)

        for e in entries:
            start(e)

    # ------------------------------------------------------------- queries
    def all_matches(self) -> list[IstdMatch]:
        return [m for per_run in self.matches.values() for m in per_run.values()]

    def matches_for(self, index: int) -> dict[str, IstdMatch]:
        p = self.project
        return self.matches.get(p.entries[index].name, {}) if 0 <= index < len(p.entries) else {}

    # ------------------------------------------------------------- files
    def to_dict(self) -> dict:
        return dict(
            library=[e.to_dict() for e in self.library],
            reference_run=self.reference_run,
            references={name: dict(name=r.name, cas=r.cas, run=r.run, t1=r.t1, t2=r.t2, mz0=r.mz0,
                                   spectrum=r.spectrum.tolist())
                       for name, r in self.references.items()},
            matches={run: {name: dict(run=m.run, name=m.name, peak_id=m.peak_id, t1=m.t1, t2=m.t2,
                                      dt1=m.dt1, dt2=m.dt2, position_score=m.position_score,
                                      spectral_score=m.spectral_score, forward=m.forward, reverse=m.reverse,
                                      combined_score=m.combined_score, tier=m.tier)
                           for name, m in per_run.items()}
                    for run, per_run in self.matches.items()},
            settings=self.settings.to_dict())

    def from_dict(self, d: dict) -> None:
        d = d or {}
        if d.get("library"):
            self.library = [IstdLibraryEntry.from_dict(e) for e in d["library"]]
        self.reference_run = d.get("reference_run")
        self.references = {
            name: IstdReference(name=r["name"], cas=r.get("cas", ""), run=r["run"], t1=r["t1"], t2=r["t2"],
                               mz0=r["mz0"], spectrum=np.asarray(r["spectrum"], dtype=float))
            for name, r in d.get("references", {}).items()}
        self.matches = {
            run: {name: IstdMatch(**m) for name, m in per_run.items()}
            for run, per_run in d.get("matches", {}).items()}
        self.settings = IstdMatchSettings.from_dict(d.get("settings", {}))
        self.changed.emit()
