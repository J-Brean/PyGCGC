"""📚 Spectral library matching for the app: the loaded library, each run's matches and the user's decisions.

The science is in `core/libmatch.py` (scoring, retention-index gate, flags) and `core/nist_sdf.py` (reading
the NIST files). This holds what the window needs between uses: the settings, the library (loaded once, it is
large), the results of matching each run, and the decisions the user made about each peak.

Results are not saved: they are recomputed in a minute from the saved settings. Decisions are saved, each with
a *content signature* of the inputs it was made on (the run's fitted components, the library, the settings,
the ladder, the fold). A decision is shown only while the run's current matches have the same signature, so
a decision never silently attaches to a peak that has since moved, been refitted or been matched against
different evidence. A decision never changes the confidence level shown beside it.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

from PySide6 import QtCore

from ..core.compare import pair_mapping
from ..core.identify import load_msp
from ..core.libmatch import MatchSettings, match_queries, queries_from_peaks
from ..core.nist_sdf import compact_from_entries, open_library
from .blob_state import MS_AWARE

ACTIONS = ("accepted", "rejected", "chosen", "unknown")


@dataclass
class RunMatches:
    matches: list
    content_signature: str
    runtime_signature: tuple
    notes: list = field(default_factory=list)
    ladder_note: str = ""
    mz_range: tuple = (0, 0)
    seconds: float = 0.0


class MatchState(QtCore.QObject):
    changed = QtCore.Signal()               # results, decisions or settings changed
    libraryChanged = QtCore.Signal()
    runningChanged = QtCore.Signal(bool)

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.settings = MatchSettings()
        self.source = ""                    # a NIST folder or an MSP file
        self.library = None
        self.library_status = "No library loaded."
        self.loading = False
        self.ladder_run = ""                # "" means each run's own ladder
        self.results: dict[str, RunMatches] = {}
        self.decisions: dict[str, dict[str, dict]] = {}
        self.running = False

    # ------------------------------------------------------------- library
    def load_library(self, source, on_done=None) -> bool:
        """Read a NIST folder (with its retention indices for the chosen column class) or an MSP file in a job."""
        source = Path(source)
        if self.loading:
            return False
        if not source.exists():
            self.library_status = f"Not found: {source}"
            self.libraryChanged.emit()
            return False
        column_class = self.settings.column_class
        self.loading = True
        self.library_status = f"Reading {source.name}…"
        self.libraryChanged.emit()

        def work(token, progress):
            def report(fraction, text=""):
                token.check()
                progress(int(1000 * fraction), 1000, text or "Reading the library")
            if source.is_dir():
                return open_library(source, column_class, progress=report, cancel=token.check)
            entries = load_msp(source)
            if not entries:
                raise ValueError(f"{source.name} parsed to zero entries: check it is a plain-text MSP export")
            return compact_from_entries(entries, label=source.name)

        def done(library):
            self.loading = False
            self.library = library
            self.source = str(source)
            self.results = {}
            n_ri = int((library.ri_exp_n > 0).sum()) if len(library) else 0
            self.library_status = (f"{library.label}: {len(library):,} spectra; literature retention index for "
                                   f"{n_ri:,} ({column_class.lower()})" if source.is_dir() else
                                   f"{library.label}: {len(library):,} spectra; {n_ri:,} with a retention index")
            self.libraryChanged.emit()
            self.changed.emit()
            self.project.request_autosave()
            if on_done:
                on_done()

        def failed(message):
            self.loading = False
            self.library_status = f"Could not read {source.name}: {message}"
            self.libraryChanged.emit()

        def cancelled():
            self.loading = False
            self.library_status = "Reading the library was cancelled."
            self.libraryChanged.emit()

        self.project.jobs.submit("📚 Reading the spectral library", work, done, failed, on_cancel=cancelled)
        return True

    def set_settings(self, settings: MatchSettings) -> None:
        old = self.settings
        self.settings = settings
        if settings.column_class != old.column_class and self.source and Path(self.source).is_dir():
            self.load_library(self.source)          # the retention indices depend on the column class
        self.changed.emit()
        self.project.request_autosave()

    def set_ladder_run(self, name: str) -> None:
        self.ladder_run = name or ""
        self.changed.emit()
        self.project.request_autosave()

    # ------------------------------------------------------------- inputs
    def ladder_for(self, index: int):
        """`(series, borrowed, t1_map, note)` for run `index`: its own ladder, or the chosen run's through the alignment."""
        p = self.project
        entry = p.entries[index]
        own = p.hydrocarbons.series_for(index)
        if own is not None and len(own) >= 2 and self.ladder_run in ("", entry.name):
            carbons = own.carbons()
            return own, False, None, f"own alkane ladder, C{int(carbons.min())} to C{int(carbons.max())}"
        other_index = p.index_of(self.ladder_run) if self.ladder_run else None
        if other_index is not None:
            other = p.hydrocarbons.series_for(other_index)
            if other is not None and len(other) >= 2:
                tf_this, tf_other = p.align.transforms.get(entry.name), p.align.transforms.get(self.ladder_run)
                reference = p.align.reference
                aligned = (tf_this is not None or entry.name == reference) and (tf_other is not None or self.ladder_run == reference)
                mapping = pair_mapping(tf_this, tf_other)
                note = (f"alkane ladder borrowed from {self.ladder_run}, through the alignment" if aligned else
                        f"alkane ladder borrowed from {self.ladder_run} with no alignment between the runs: "
                        "the same first-dimension times are used")
                return other, True, (lambda t1, t2: mapping(t1, t2)[0]), note
        return None, False, None, "no alkane ladder: matching rests on the spectrum alone"

    def ladder_signature(self, index: int) -> tuple:
        series, borrowed, _map, _note = self.ladder_for(index)
        if series is None:
            return ()
        return (borrowed, self.ladder_run, tuple((int(c), round(float(series.hits[int(c)].t1), 2)) for c in series.carbons()))

    def content_signature(self, index: int) -> str:
        """What a decision was made on; stable between sessions (no object identities)."""
        p, b = self.project, self.project.blobs
        entry = p.entries[index]
        library = (self.library.label, len(self.library)) if self.library is not None else None
        parts = (b.component_signature(index) if b.method == MS_AWARE else len(b.blobs_for(index) or ()),
                 library, json.dumps(self.settings.to_dict(), sort_keys=True), self.ladder_signature(index),
                 round(float(p.period), 6), round(float(entry.phase), 6), b.method)
        return hashlib.sha256(repr(parts).encode()).hexdigest()[:16]

    def runtime_signature(self, index: int) -> tuple:
        p = self.project
        return (self.content_signature(index), id(p.base_run(index)), id(self.library))

    def why_not(self, index: int) -> list[str]:
        """What is missing before run `index` can be matched, as a list of sentences."""
        p, out = self.project, []
        if not 0 <= index < len(p.entries) or p.entries[index].run is None:
            return ["Select a loaded run."]
        if self.library is None:
            out.append("Load a library: the NIST folder, or an MSP file." if not self.loading else "The library is still loading.")
        peaks_ready = (p.blobs.components_for(index) is not None if p.blobs.method == MS_AWARE
                       else p.blobs.blobs_for(index) is not None)
        if not peaks_ready:
            out.append("Find peaks first (step 4): the matcher identifies the detected peaks.")
        energy = (p.entries[index].import_meta or {}).get("energy_ev")
        if energy is not None and energy != 70:
            out.append(f"This is the {energy:g} eV channel; library spectra are 70 eV. Select the 70 eV run.")
        return out

    # ------------------------------------------------------------- running
    def current_for(self, index: int) -> RunMatches | None:
        """The run's matches if they were made on the inputs it has now, else None."""
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        found = self.results.get(p.entries[index].name)
        return found if found is not None and found.runtime_signature == self.runtime_signature(index) else None

    def results_for(self, index: int) -> RunMatches | None:
        """The run's last matches, even if its inputs have changed since."""
        p = self.project
        return self.results.get(p.entries[index].name) if 0 <= index < len(p.entries) else None

    def match_runs(self, indices: list[int], on_done=None) -> bool:
        p = self.project
        indices = [i for i in indices if 0 <= i < len(p.entries) and not self.why_not(i)]
        if not indices or self.running:
            return False
        settings, library = self.settings, self.library
        before = {i: self.runtime_signature(i) for i in indices}
        prepared = {i: self.ladder_for(i) for i in indices}
        labels = {i: self._class_labels(i) for i in indices}
        self.running = True
        self.runningChanged.emit(True)

        def work(token, progress):
            import time
            out = {}
            for number, index in enumerate(indices):
                token.check()
                started = time.time()
                name = p.entries[index].name
                progress(number, len(indices), f"Matching {name}")
                peaks = p.blobs.full_spectrum_peaks(index)
                if peaks is None or len(peaks) == 0:
                    out[index] = RunMatches([], "", (), ["no peaks"])
                    continue
                series, borrowed, t1_map, note = prepared[index]
                queries = queries_from_peaks(peaks, series, period=p.period, settings=settings, borrowed=borrowed,
                                             t1_map=t1_map, class_labels=labels[index])
                run = p.base_run(index)
                mz_range = (int(run.mz[0]), int(run.mz[-1]))
                matches = match_queries(queries, library, settings, mz_range=mz_range, cancel=token.check)
                out[index] = RunMatches(matches, "", (), [], note, mz_range, time.time() - started)
            return out

        def finish():
            self.running = False
            self.runningChanged.emit(False)

        def done(out):
            finish()
            if any(self.runtime_signature(i) != before[i] for i in indices):
                p.message.emit("Library matching discarded: the peaks, library or settings changed while it ran. Match again.")
                return
            for index, found in out.items():
                found.content_signature = self.content_signature(index)
                found.runtime_signature = self.runtime_signature(index)
                self.results[p.entries[index].name] = found
            self.changed.emit()
            self.project.request_autosave()
            if on_done:
                on_done()

        def failed(message):
            finish()
            p.message.emit(f"Library matching failed: {message}")

        p.jobs.submit("📚 Matching peaks to the library", work, done, failed, on_cancel=finish)
        return True

    def _class_labels(self, index: int) -> list | None:
        """Each peak's 🧬 class where one was found with confidence, in the order of the peaks being matched."""
        p = self.project
        result, used = p.classes.result_for(index), p.classes.blobs_classified(index)
        peaks_ready = p.blobs.components_for(index) if p.blobs.method == MS_AWARE else p.blobs.blobs_for(index)
        if result is None or used is None or p.classes.is_stale(index) or peaks_ready is None:
            return None
        n = len(getattr(peaks_ready, "components", peaks_ready))
        if len(result.calls) != n:
            return None
        return [c.cls if c.cls and c.tier in ("confident", "probable") else "" for c in result.calls]

    # ------------------------------------------------------------- decisions
    def decide(self, index: int, label: str, action: str, name: str = "", cas: str = "") -> None:
        if action not in ACTIONS:
            raise ValueError(f"unknown decision '{action}'")
        found = self.current_for(index)
        if found is None:
            return
        run = self.project.entries[index].name
        self.decisions.setdefault(run, {})[label] = dict(action=action, name=name, cas=cas,
                                                          signature=found.content_signature)
        self.changed.emit()
        self.project.request_autosave()

    def clear_decision(self, index: int, label: str) -> None:
        run = self.project.entries[index].name
        if self.decisions.get(run, {}).pop(label, None) is not None:
            self.changed.emit()
            self.project.request_autosave()

    def decisions_for(self, index: int) -> dict:
        """Decisions that still stand: made on the same inputs the current matches were."""
        found = self.current_for(index)
        if found is None:
            return {}
        run = self.project.entries[index].name
        return {label: d for label, d in self.decisions.get(run, {}).items() if d.get("signature") == found.content_signature}

    def forget_runs(self, entries) -> None:
        for e in entries:
            self.results.pop(e.name, None)
            self.decisions.pop(e.name, None)

    # ------------------------------------------------------------- files
    def to_dict(self) -> dict:
        return dict(settings=self.settings.to_dict(), source=self.source, ladder_run=self.ladder_run,
                    decisions={run: dict(d) for run, d in self.decisions.items() if d})

    def from_dict(self, d: dict) -> None:
        d = d or {}
        self.settings = MatchSettings.from_dict(d.get("settings", {}))
        self.source = d.get("source", "")
        self.ladder_run = d.get("ladder_run", "")
        self.decisions = {run: dict(v) for run, v in d.get("decisions", {}).items()}
        self.results = {}
        self.library = None
        self.library_status = "No library loaded." if not self.source else f"Library to load: {self.source}"
        self.libraryChanged.emit()
        self.changed.emit()
