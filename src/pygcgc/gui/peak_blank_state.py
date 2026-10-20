"""➖ Subtracting a blank run's peaks from every other run: the settings and the results.

The pairing and arithmetic are in `core/peak_blank.py`. This holds the user's choices (which run is the
blank, the factor, the matching windows), the peak-by-peak results for each run, and the link to
quantification: while "apply to quantified amounts" is on, the same blank run and factor are written into
`QuantSettings`, so 🧪 reports amounts as sample minus blank and marks its results out of date when they change.

Results are not saved (they are recomputed in seconds from the saved settings); the settings are.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from PySide6 import QtCore

from ..core.compare import pair_mapping
from ..core.peak_blank import PeakBlankSettings, subtract_blank_peaks, suggest_factor, summarise
from .blob_state import MS_AWARE


@dataclass
class RunSubtraction:
    rows: list
    summary: dict
    suggestion: dict
    runtime_signature: tuple = ()
    peaks: object = field(default=None, repr=False)     # the sample's peak table, for display


class PeakBlankState(QtCore.QObject):
    changed = QtCore.Signal()
    runningChanged = QtCore.Signal(bool)

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.settings = PeakBlankSettings()
        self.apply_to_quant = True
        self.results: dict[str, RunSubtraction] = {}
        self.running = False

    # ------------------------------------------------------------- settings
    def set_settings(self, settings: PeakBlankSettings, apply_to_quant: bool | None = None) -> None:
        if apply_to_quant is not None:
            self.apply_to_quant = bool(apply_to_quant)
        self.settings = settings
        self._push_to_quant()
        self.changed.emit()
        self.project.request_autosave()

    def _push_to_quant(self) -> None:
        """Keep 🧪's blank settings equal to these while the link is on, and empty when it is off."""
        quant = self.project.quant
        s, q = self.settings, quant.settings
        wanted = (dict(blank_run=s.blank_run, blank_factor=s.factor, blank_max_dt1=s.max_dt1, blank_max_dt2=s.max_dt2,
                       blank_min_similarity=s.min_similarity, blank_significance=s.significance)
                  if self.apply_to_quant else dict(blank_run=""))
        new = replace(q, **wanted)
        if new != q:
            quant.set_settings(new)

    # ------------------------------------------------------------- inputs
    def blank_index(self) -> int | None:
        return self.project.index_of(self.settings.blank_run) if self.settings.blank_run else None

    def why_not(self, index: int) -> list[str]:
        """What is missing before run `index` can have the blank taken from it."""
        p, out = self.project, []
        b = self.blank_index()
        if b is None:
            out.append("Choose the blank run.")
        elif index == b:
            out.append("This is the blank run itself.")
        if not 0 <= index < len(p.entries) or p.entries[index].run is None:
            return out + ["Select a loaded run."]

        def peaks_ready(i):
            return (p.blobs.components_for(i) is not None if p.blobs.method == MS_AWARE
                    else p.blobs.blobs_for(i) is not None)

        if not peaks_ready(index):
            out.append("Find peaks in this run first (step 4).")
        if b is not None and not peaks_ready(b):
            out.append("Find peaks in the blank run first (step 4).")
        if b is not None and index != b:
            e1 = (p.entries[index].import_meta or {}).get("energy_ev")
            e2 = (p.entries[b].import_meta or {}).get("energy_ev")
            if e1 != e2:
                out.append(f"The blank is {e2} eV and this run {e1} eV: use runs of the same ionisation energy.")
        return out

    def runtime_signature(self, index: int) -> tuple:
        p, b = self.project, self.blank_index()
        sample, blank = p.entries[index], p.entries[b] if b is not None else None
        tf = lambda name: (lambda t: None if t is None else (t.kind, tuple(t.anchors)))(p.align.transforms.get(name))
        sig = lambda i: (p.blobs.component_signature(i) if p.blobs.method == MS_AWARE else len(p.blobs.blobs_for(i) or ()),
                         id(p.base_run(i)))
        return (self.settings, sig(index), sig(b) if b is not None else None, tf(sample.name),
                tf(blank.name) if blank else None, round(float(p.period), 6))

    def current_for(self, index: int) -> RunSubtraction | None:
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        found = self.results.get(p.entries[index].name)
        if found is None or self.why_not(index) or found.runtime_signature != self.runtime_signature(index):
            return None
        return found

    def results_for(self, index: int) -> RunSubtraction | None:
        p = self.project
        return self.results.get(p.entries[index].name) if 0 <= index < len(p.entries) else None

    # ------------------------------------------------------------- running
    def subtract_runs(self, indices: list[int] | None = None, on_done=None) -> bool:
        p = self.project
        b = self.blank_index()
        if indices is None:
            indices = [i for i in p.loaded() if i != b]
        indices = [i for i in indices if i != b and not self.why_not(i)]
        if not indices or b is None or self.running:
            return False
        settings = self.settings
        before = {i: self.runtime_signature(i) for i in indices}
        maps = {i: pair_mapping(p.align.transforms.get(settings.blank_run), p.align.transforms.get(p.entries[i].name))
                for i in indices}                                  # the blank's coordinates onto each sample's
        self.running = True
        self.runningChanged.emit(True)

        def work(token, progress):
            blank = p.blobs.full_spectrum_peaks(b)
            out = {}
            for number, index in enumerate(indices):
                token.check()
                progress(number, len(indices), f"Subtracting the blank from {p.entries[index].name}")
                peaks = p.blobs.full_spectrum_peaks(index)
                rows = subtract_blank_peaks(peaks, blank, settings, maps[index])
                out[index] = RunSubtraction(rows, summarise(rows, settings), suggest_factor(rows), peaks=peaks)
            return out

        def finish():
            self.running = False
            self.runningChanged.emit(False)

        def done(out):
            finish()
            if any(self.runtime_signature(i) != before[i] for i in indices):
                p.message.emit("Blank subtraction discarded: the peaks, alignment or settings changed while it ran. Run it again.")
                return
            for index, found in out.items():
                found.runtime_signature = self.runtime_signature(index)
                self.results[p.entries[index].name] = found
            self.changed.emit()
            if on_done:
                on_done()

        def failed(message):
            finish()
            p.message.emit(f"Blank subtraction failed: {message}")

        p.jobs.submit("➖ Subtracting the blank's peaks", work, done, failed, on_cancel=finish)
        return True

    def forget_runs(self, entries) -> None:
        names = {e.name for e in entries}
        for name in names:
            self.results.pop(name, None)
        if self.settings.blank_run in names:
            self.settings = replace(self.settings, blank_run="")
            self._push_to_quant()

    # ------------------------------------------------------------- files
    def to_dict(self) -> dict:
        return dict(settings=self.settings.to_dict(), apply_to_quant=self.apply_to_quant)

    def from_dict(self, d: dict) -> None:
        d = d or {}
        self.settings = PeakBlankSettings.from_dict(d.get("settings", {}))
        self.apply_to_quant = bool(d.get("apply_to_quant", True))
        self.results = {}
        self.changed.emit()
