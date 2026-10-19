"""🧬 Compound-class state for the app.

Keeps one `ClassResult` per run, the rule set in use (the shipped
`presets/chemical_classes.toml` or a file of the user's own), which classes
the user has switched off, and whether to use the 14 eV channel of a tandem
acquisition when one is loaded alongside. Results are recomputed on demand
and not saved, like 🔍's blobs and ⛓️'s ladder; the settings are saved.
The science is all in `core/chemclass.py`.

Blobs: the classifier wants full spectra, so it uses the legacy detector's
blobs. If the user has run only the MS-aware detector (the default), the
legacy detector is run here, privately, on the same cleaned run.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from PySide6 import QtCore

from ..core.background import remove_background, storage_threshold
from ..core.chemclass import (ClassSettings, RuleSet, classify_blobs, load_rules, soft_spectra_for_blobs)
from ..core.fold import Fold
from ..core.peaks import find_blobs
from ..core.peaktable import respectrum_components
from .blob_state import MS_AWARE

HARD_ENERGIES = (None, 70.0)


class ClassState(QtCore.QObject):
    changed = QtCore.Signal()          # results, rules or settings changed
    runningChanged = QtCore.Signal(bool)

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.settings = ClassSettings()
        self.rules_path: Path | None = None          # None: the shipped rules
        self.rules: RuleSet = load_rules()
        self.enabled: dict[str, bool] = {}           # the user's overrides of each class's `enabled`
        self.use_soft = True
        self.results: dict = {}                      # RunEntry -> ClassResult
        self.blobs_used: dict = {}                   # RunEntry -> the BlobSet that was classified
        self.notes: dict = {}                        # RunEntry -> what the last run used, for the panel
        self.own_blobs: dict = {}                    # RunEntry -> BlobSet found here (MS-aware users)
        self.stale: set = set()
        self.running = False
        self._input_revision = 0

    # ------------------------------------------------------------ queries
    def forget_runs(self, entries) -> None:
        for e in entries:
            for d in (self.results, self.blobs_used, self.notes, self.own_blobs):
                d.pop(e, None)
            self.stale.discard(e)

    def result_for(self, index: int):
        p = self.project
        return self.results.get(p.entries[index]) if 0 <= index < len(p.entries) else None

    def blobs_classified(self, index: int):
        p = self.project
        return self.blobs_used.get(p.entries[index]) if 0 <= index < len(p.entries) else None

    def is_stale(self, index: int) -> bool:
        p = self.project
        return 0 <= index < len(p.entries) and p.entries[index] in self.stale

    def rules_in_use(self) -> RuleSet:
        """The loaded rules with the user's on/off choices applied."""
        classes = [replace(c, enabled=self.enabled.get(c.name, c.enabled)) for c in self.rules.classes]
        return replace(self.rules, classes=classes)

    def energy_of(self, index: int):
        meta = self.project.entries[index].import_meta or {}
        return meta.get("energy_ev")

    def soft_partner(self, index: int) -> int | None:
        """A loaded run from the same acquisition at a lower ionisation energy (the 14 eV channel of a tandem file)."""
        p = self.project
        meta = p.entries[index].import_meta or {}
        acq, energy = meta.get("acquisition_id"), meta.get("energy_ev")
        if acq is None or energy is None:
            return None
        for j, e in enumerate(p.entries):
            m = e.import_meta or {}
            if (j != index and e.run is not None and m.get("acquisition_id") == acq
                    and m.get("energy_ev") is not None and m["energy_ev"] < energy
                    and m.get("processing_state", "instrument") == "instrument"):
                return j
        return None

    def why_not(self, index: int) -> str:
        """Empty if run `index` can be classified, otherwise the reason."""
        p = self.project
        if not 0 <= index < len(p.entries) or p.entries[index].run is None:
            return "Select a loaded run."
        energy = self.energy_of(index)
        if energy not in HARD_ENERGIES:
            return (f"This is the {energy:g} eV channel. Classes are read from the 70 eV spectrum; select the "
                    "70 eV run (the 14 eV run is used automatically as its partner).")
        return ""

    def _set_running(self, on: bool) -> None:
        self.running = on
        self.runningChanged.emit(on)

    # -------------------------------------------------------------- rules
    def set_rules_path(self, path) -> None:
        """Load a rule file. A mistake raises ValueError with the class and criterion it is in."""
        rules = load_rules(path)
        self.rules, self.rules_path, self.enabled = rules, Path(path), {}
        self._rules_changed()

    def reset_rules(self) -> None:
        self.rules, self.rules_path, self.enabled = load_rules(), None, {}
        self._rules_changed()

    def reload_rules(self) -> None:
        self.rules = load_rules(self.rules_path)
        self._rules_changed()

    def set_enabled(self, name: str, on: bool) -> None:
        self.enabled[name] = bool(on)
        self._rules_changed()

    def set_use_soft(self, on: bool) -> None:
        self.use_soft = bool(on)
        self._rules_changed()

    def set_settings(self, settings: ClassSettings) -> None:
        self.settings = settings
        self._rules_changed()

    def _rules_changed(self) -> None:
        self._input_revision += 1
        self.stale = set(self.results)
        self.changed.emit()

    def mark_stale(self, *_) -> None:
        self._input_revision += 1
        if self.results:
            self.stale = set(self.results)
            self.changed.emit()

    # ---------------------------------------------------------- classifying
    def classify(self, index: int, on_done=None) -> bool:
        p = self.project
        reason = self.why_not(index)
        if reason:
            p.message.emit(f"🧬 {reason}")
            return False
        if self.running:
            return False
        e = p.entries[index]
        fold = p.fold(index)
        base_run = p.base_run(index)
        blobs = p.blobs.blobs_for(index)
        component_set = p.blobs.components_for(index) if p.blobs.method == MS_AWARE else None
        series = p.hydrocarbons.series_for(index)
        rules, settings = self.rules_in_use(), self.settings
        partner = self.soft_partner(index) if self.use_soft else None
        soft_entry = p.entries[partner] if partner is not None else None
        soft_run = p.base_run(partner) if partner is not None else None
        soft_cleaned = soft_entry is not None and soft_entry.background is not None
        blob_settings = p.blobs.settings
        raw = e.run
        inputs = p.analysis_inputs()
        configuration = (self._input_revision, self.settings, self.rules_in_use(), dict(self.enabled), self.use_soft,
                         p.blobs.settings, p.blobs.component_settings, p.blobs.method)
        self._set_running(True)

        def work(token, progress):
            found_here = blobs is None and component_set is None
            b = blobs
            from_components = b is None and component_set is not None and bool(component_set.groups)
            if from_components:
                # The fitted overlapping peaks, with spectra re-solved on every ion. No detector runs here.
                progress(0, 3, "Fitting full spectra of the peaks")
                b, _model, _footprint = respectrum_components(base_run, fold, component_set)
            elif b is None:
                progress(0, 3, "Finding blobs")
                b = find_blobs(base_run, fold, blob_settings, floor=storage_threshold(raw) / 3.0)
            token.check()
            soft = None
            if soft_run is not None and not from_components:
                src = soft_run
                if not soft_cleaned:
                    progress(1, 3, "Removing background from the 14 eV run")
                    src = remove_background(soft_run, Fold(soft_run.times, fold.spec)).run
                token.check()
                progress(2, 3, "Reading the 14 eV spectra")
                soft = soft_spectra_for_blobs(src, fold, b, settings)
            progress(3, 3, "Classifying")
            result = classify_blobs(b, rules, series, settings, soft)
            return b, result, (found_here, from_components)

        def done(out):
            if (inputs != p.analysis_inputs() or configuration != (
                    self._input_revision, self.settings, self.rules_in_use(), self.enabled, self.use_soft,
                    p.blobs.settings, p.blobs.component_settings, p.blobs.method)):
                self._set_running(False)
                p.message.emit("Classification discarded: inputs changed; classify again")
                return
            b, result, (found_here, from_components) = out
            self.results[e], self.blobs_used[e] = result, b
            if found_here:
                self.own_blobs[e] = b
            source = (" (fitted overlapping peaks, spectra on all ions)" if from_components else
                      " (found here: legacy detector)" if found_here else " (from 🔍)")
            parts = [f"{len(b)} peaks" + source]
            parts.append("ladder used" if series is not None else "no ladder: t2 evidence skipped")
            if soft_entry is not None:
                parts.append(f"14 eV from {soft_entry.name}" + ("" if soft_cleaned else " (background removed here)"))
            elif self.use_soft:
                parts.append("no 14 eV run loaded")
            if from_components and soft_entry is not None:
                parts.append("14 eV spectra are not used with fitted peaks yet (they need pixel footprints)")
            self.notes[e] = "; ".join(parts)
            self.stale.discard(e)
            self._set_running(False)
            self.changed.emit()
            if on_done:
                on_done(result)

        def failed(msg):
            self._set_running(False)
            p.message.emit(f"🧬 {e.name}: {msg}")

        p.jobs.submit(f"🧬 Classifying blobs: {e.name}", work, done, failed,
                      on_cancel=lambda: self._set_running(False))
        return True

    # -------------------------------------------------------------- files
    def to_dict(self) -> dict:
        # Results are not stored: classify again after reopening.
        return dict(settings=self.settings.to_dict(), rules_path=str(self.rules_path) if self.rules_path else None,
                    enabled=dict(self.enabled), use_soft=self.use_soft)

    def from_dict(self, d: dict) -> None:
        d = d or {}
        self.settings = ClassSettings.from_dict(d.get("settings", {}))
        self.use_soft = bool(d.get("use_soft", True))
        self.rules_path, self.rules = None, load_rules()
        path = d.get("rules_path")
        if path and Path(path).exists():
            try:
                self.rules, self.rules_path = load_rules(path), Path(path)
            except ValueError as exc:
                self.project.message.emit(f"🧬 Saved rule file {Path(path).name} is invalid ({exc}); using the shipped rules")
        self.enabled = {k: bool(v) for k, v in (d.get("enabled") or {}).items() if k in {c.name for c in self.rules.classes}}
        self.results, self.blobs_used, self.notes, self.own_blobs, self.stale = {}, {}, {}, {}, set()
        self.changed.emit()
