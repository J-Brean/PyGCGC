"""🧪 Quantification state for the app.

Holds the full target table, calibration levels, settings and review
decisions, runs the per-run measurements as background jobs, and keeps the
resulting QuantReport. The science is all in `core/quant.py`.

Flow of 🧪 Quantify:
1. Measure the calibration runs with the table positions (in parallel).
2. Learn expected positions from them (if enabled).
3. Measure every run with those positions (in parallel).
4. Fit the curves and quantify, which is instant.

Review decisions (reject, accept, exclude a point, change a curve's
settings) only need step 4 again. Picking a different peak by hand
re-measures that one run.
"""

from __future__ import annotations

from pathlib import Path
from dataclasses import asdict, replace

from PySide6 import QtCore

from ..core import (QuantDecisions, QuantReport, QuantSettings, SequenceRun, finish_sequence,
                    learn_positions, learn_qual_ratios, load_levels, load_targets, quantify_run,
                    with_positions)
from ..core.targets import Target


class QuantState(QtCore.QObject):
    changed = QtCore.Signal()            # report, decisions, targets or staleness changed
    runningChanged = QtCore.Signal(bool)

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.targets = None
        self.targets_error = ""
        self._inline_targets = False
        self.levels: dict | None = None
        self.levels_path: Path | None = None
        self.amount_unit = "calibration units"
        self.calibration_setup: dict = {}
        self.settings = QuantSettings()
        self.decisions = QuantDecisions()
        self.report: QuantReport | None = None
        self.positions: dict = {}
        self.ratios: dict = {}           # qualifier ratios learned where the table leaves them blank
        self.runquants: dict = {}        # RunEntry -> RunQuant
        self.stale = False
        self.running = False
        self.was_run = False
        self._input_revision = 0

    # ------------------------------------------------------------ inputs
    def set_targets(self, path) -> None:
        self._inline_targets = False
        self.calibration_setup = {}
        try:
            self.targets = load_targets(path)
            self.targets_error = ""
        except (ValueError, KeyError) as exc:
            self.targets, self.targets_error = None, str(exc)
            self.project.message.emit(f"Target table can't be used for 🧪: {exc}")
        self._invalidate()

    def set_levels(self, path) -> None:
        levels = load_levels(path)
        self.calibration_setup = {}
        self.amount_unit = "calibration units"
        self.levels_path = Path(path)
        self.levels = levels
        self._input_revision += 1
        self.project.message.emit(f"Calibration levels from {self.levels_path.name}: "
                                  f"{', '.join(self.levels)}")
        if self.running:
            self.stale = True
            self.changed.emit()
        elif self.runquants:
            self.refit()
        else:
            self.changed.emit()

    def set_settings(self, settings: QuantSettings) -> None:
        self.settings = settings
        self._invalidate()

    def ready(self) -> tuple[bool, str]:
        if not self.project.loaded():
            return False, "Open some runs first"
        if hasattr(self.project, 'istd') and self.project.istd.outdated_reviews():
            return False, "Reviewed standard choices are out of date; search Standards and review the new matches"
        runs = [self.project.entries[i].run for i in self.project.loaded()]
        energies = {r.meta.get('energy_ev') for r in runs if r.meta.get('energy_ev') is not None}
        if len(energies) > 1:
            return False, 'Quantification cannot mix ionisation energies; select one energy channel'
        if any(r.meta.get('format') == 'Markes LSC' for r in runs) and not self.project.period_confirmed:
            return False, 'Confirm the modulation period in the Fold panel before native-data quantification'
        if self.targets is None:
            return False, self.targets_error or "Load a target table first (File menu)"
        if not any(t.is_istd for t in self.targets):
            return False, "The target table has no internal standards (role = istd)"
        return True, ""

    # ------------------------------------------------------------ queries
    def info(self) -> list[SequenceRun]:
        out = []
        for e in self.project.entries:
            out.append(SequenceRun(e.name, e.kind, e.level, e.cal_set, e.volume, self.runquants.get(e)))
        return out

    def run_quant(self, index: int):
        if not 0 <= index < len(self.project.entries):
            return None
        return self.runquants.get(self.project.entries[index])

    def used_targets(self):
        if self.targets is None:
            return None
        return with_positions(self.targets, self.positions, self.ratios)

    def calibration_entries(self) -> list:
        info = self.info()
        return [e for e, i in zip(self.project.entries, info) if i.is_calibration and e.run is not None]

    # ------------------------------------------------------------ running
    def _invalidate(self) -> None:
        self._input_revision += 1
        if self.report is not None or self.runquants or self.running:
            self.stale = True
        self.changed.emit()

    def mark_stale(self, *_):
        if self.runquants or self.running:
            self._input_revision += 1
            self.stale = True
            self.changed.emit()

    def _set_running(self, on: bool) -> None:
        self.running = on
        self.runningChanged.emit(on)

    def _measure(self, entries, targets, label, on_done) -> None:
        """Measure some runs in parallel; call on_done() when all have finished."""
        p = self.project
        remaining = set(entries)
        unsuccessful = set()
        if not entries:
            on_done()
            return
        settings = self.settings

        def settled(entry, success):
            if not success:
                unsuccessful.add(entry)
            remaining.discard(entry)
            if remaining:
                return
            if unsuccessful:
                self.stale = True
                self._set_running(False)
                self.changed.emit()
            else:
                on_done()
        for e in entries:
            idx = p.entries.index(e)
            bg = e.background
            clean = bg.run if bg is not None else e.run
            noise = bg.noise if bg is not None else None
            raw, fold = e.run, p.fold(idx)
            assigned = {c: v for (r, c), v in self.decisions.assigned.items() if r == e.name}
            reviewed = p.istd.reviewed_matches(e.name)
            assigned.update({name: (match.t1, match.t2) for name, match in reviewed.items() if match.found})
            missing_standards = [name for name, match in reviewed.items() if not match.found]
            polys = {c: v for (r, c), v in self.decisions.footprints.items() if r == e.name}

            def work(token, progress, raw=raw, clean=clean, fold=fold, noise=noise, assigned=assigned,
                     polys=polys, missing_standards=missing_standards):
                result = quantify_run(raw, clean, fold, noise, targets, settings, assigned, progress,
                                      footprints=polys)
                for name in missing_standards:
                    if name in result.hits:
                        result.hits[name] = replace(result.hits[name], found=False,
                            quant_volume=float('nan'), tic_volume=float('nan'), manual=True)
                return result

            def done(rq, e=e):
                self.runquants[e] = rq
                settled(e, True)

            def failed(text, e=e):
                p.message.emit(f"🧪 {e.name} failed: {text}")
                settled(e, False)

            p.jobs.submit(f"🧪 {label}: {e.name}", work, done, failed,
                          on_cancel=lambda e=e: settled(e, False))

    def run(self, on_done=None) -> bool:
        ok, why = self.ready()
        if not ok or self.running:
            if not ok:
                self.project.message.emit(why)
            return False
        self._set_running(True)
        input_revision = self._input_revision
        self.runquants = {}
        self.report = None
        entries = [self.project.entries[i] for i in self.project.loaded()]
        targets = self.targets

        def second_pass():
            used = self.used_targets()
            self._measure(entries, used, "Measuring", finish)

        def finish():
            if input_revision != self._input_revision:
                self.stale = True
                self._set_running(False)
                self.changed.emit()
                self.project.message.emit("Quantification inputs changed during processing. Rerun Quantify.")
                return
            self.finish()
            self._set_running(False)
            self.was_run = True
            if on_done:
                on_done()

        if self.settings.learn_positions:
            cal = self.calibration_entries() or entries

            def learned():
                info = [i for i in self.info() if i.quant is not None]
                self.positions = learn_positions(info, targets, self.settings)
                self.ratios = learn_qual_ratios(info, targets)
                second_pass()

            self.positions, self.ratios = {}, {}
            self._measure(cal, targets, "Learning positions", learned)
        else:
            self.positions, self.ratios = {}, {}
            second_pass()
        return True

    def finish(self) -> None:
        targets = self.used_targets()
        if targets is None:
            return
        self.report = finish_sequence(self.info(), targets, self.levels or {}, self.settings,
                                      self.decisions, self.positions)
        self.stale = False
        self.changed.emit()

    def refit(self) -> None:
        if self.runquants:
            self.finish()

    def requantify(self, entries, on_done=None) -> None:
        targets = self.used_targets()
        if targets is None:
            return
        self._set_running(True)

        def done():
            self.finish()
            self._set_running(False)
            if on_done:
                on_done()

        self._measure(list(entries), targets, "Re-measuring", done)

    def set_decisions(self, decisions: QuantDecisions, requantify_runs=()) -> None:
        self.decisions = decisions
        runs = [e for e in self.project.entries if e.name in set(requantify_runs) and e.run is not None]
        if runs:
            self.requantify(runs)
        else:
            self.refit()
            self.changed.emit()

    def forget_runs(self, entries) -> None:
        for e in entries:
            self.runquants.pop(e, None)
        if self.runquants:
            self.finish()
        else:
            self.report = None
            self.changed.emit()

    # ------------------------------------------------------------ files
    def to_dict(self, rel) -> dict:
        return {"settings": self.settings.to_dict(), "decisions": self.decisions.to_dict(),
                "levels": rel(self.levels_path), "levels_data": self.levels if self.levels_path is None else None,
                "targets_data": [asdict(t) for t in self.targets] if self._inline_targets else None,
                "amount_unit": self.amount_unit, "calibration_setup": self.calibration_setup,
                "was_run": bool(self.runquants)}

    def from_dict(self, d: dict, resolve) -> None:
        d = d or {}
        self.settings = QuantSettings.from_dict(d.get("settings", {}))
        self.decisions = QuantDecisions.from_dict(d.get("decisions", {}))
        self.amount_unit = d.get("amount_unit", "calibration units")
        self.calibration_setup = d.get("calibration_setup", {})
        self._inline_targets = d.get("targets_data") is not None
        if d.get("targets_data") is not None:
            self.targets = [Target(**{**t, "qual_ions": tuple(t.get("qual_ions", ())),
                                     "qual_ratios": tuple(t.get("qual_ratios", ()))})
                            for t in d["targets_data"]]
            self.targets_error = ""
        elif self.project.targets_path is None:
            self.targets, self.targets_error = None, ""
        self.report, self.runquants, self.positions, self.ratios = None, {}, {}, {}
        self.stale = False
        self.levels, self.levels_path = None, None
        lv = resolve(d.get("levels"))
        if lv is not None and Path(lv).exists():
            self.levels_path = Path(lv)
            self.levels = load_levels(lv)
        if d.get("levels_data") is not None:
            self.levels = d["levels_data"]
        self.was_run = bool(d.get("was_run"))
        self.changed.emit()
