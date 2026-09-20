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

from PySide6 import QtCore

from ..core import (QuantDecisions, QuantReport, QuantSettings, SequenceRun, finish_sequence,
                    learn_positions, load_levels, load_targets, quantify_run, read_sequence_extras,
                    with_positions)


class QuantState(QtCore.QObject):
    changed = QtCore.Signal()            # report, decisions, targets or staleness changed
    runningChanged = QtCore.Signal(bool)

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.targets = None
        self.targets_error = ""
        self.levels: dict | None = None
        self.levels_path: Path | None = None
        self.settings = QuantSettings()
        self.decisions = QuantDecisions()
        self.report: QuantReport | None = None
        self.positions: dict = {}
        self.runquants: dict = {}        # RunEntry -> RunQuant
        self.extras: dict = {}           # run name -> {cal_set, volume}
        self.stale = False
        self.running = False
        self.was_run = False

    # ------------------------------------------------------------ inputs
    def set_targets(self, path) -> None:
        try:
            self.targets = load_targets(path)
            self.targets_error = ""
        except (ValueError, KeyError) as exc:
            self.targets, self.targets_error = None, str(exc)
            self.project.message.emit(f"Target table can't be used for 🧪: {exc}")
        self._invalidate()

    def set_levels(self, path) -> None:
        self.levels_path = Path(path)
        self.levels = load_levels(path)
        self.project.message.emit(f"Calibration levels from {self.levels_path.name}: "
                                  f"{', '.join(self.levels)}")
        if self.runquants:
            self.refit()
        else:
            self.changed.emit()

    def set_sequence(self, path) -> None:
        try:
            self.extras = read_sequence_extras(path)
        except (OSError, KeyError, ValueError):
            self.extras = {}

    def set_settings(self, settings: QuantSettings) -> None:
        self.settings = settings
        self._invalidate()

    def ready(self) -> tuple[bool, str]:
        if not self.project.loaded():
            return False, "Open some runs first"
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
            ex = self.extras.get(e.name, {})
            out.append(SequenceRun(e.name, e.kind, e.level, ex.get("cal_set", ""),
                                   ex.get("volume"), self.runquants.get(e)))
        return out

    def run_quant(self, index: int):
        if not 0 <= index < len(self.project.entries):
            return None
        return self.runquants.get(self.project.entries[index])

    def used_targets(self):
        if self.targets is None:
            return None
        return with_positions(self.targets, self.positions) if self.positions else self.targets

    def calibration_entries(self) -> list:
        info = self.info()
        return [e for e, i in zip(self.project.entries, info) if i.is_calibration and e.run is not None]

    # ------------------------------------------------------------ running
    def _invalidate(self) -> None:
        if self.report is not None or self.runquants:
            self.stale = True
        self.changed.emit()

    def mark_stale(self, *_):
        if self.runquants and not self.running:
            self.stale = True
            self.changed.emit()

    def _set_running(self, on: bool) -> None:
        self.running = on
        self.runningChanged.emit(on)

    def _measure(self, entries, targets, label, on_done) -> None:
        """Measure some runs in parallel; call on_done() when all have finished."""
        p = self.project
        remaining = set(entries)
        if not entries:
            on_done()
            return
        settings = self.settings
        for e in entries:
            idx = p.entries.index(e)
            bg = e.background
            clean = bg.run if bg is not None else e.run
            noise = bg.noise if bg is not None else None
            raw, fold = e.run, p.fold(idx)
            assigned = {c: v for (r, c), v in self.decisions.assigned.items() if r == e.name}
            polys = {c: v for (r, c), v in self.decisions.footprints.items() if r == e.name}

            def work(token, progress, raw=raw, clean=clean, fold=fold, noise=noise, assigned=assigned,
                     polys=polys):
                return quantify_run(raw, clean, fold, noise, targets, settings, assigned, progress,
                                    footprints=polys)

            def done(rq, e=e):
                self.runquants[e] = rq
                remaining.discard(e)
                if not remaining:
                    on_done()

            def failed(text, e=e):
                p.message.emit(f"🧪 {e.name} failed: {text}")
                remaining.discard(e)
                if not remaining:
                    on_done()

            p.jobs.submit(f"🧪 {label}: {e.name}", work, done, failed)

    def run(self, on_done=None) -> bool:
        ok, why = self.ready()
        if not ok or self.running:
            if not ok:
                self.project.message.emit(why)
            return False
        self._set_running(True)
        self.runquants = {}
        self.report = None
        entries = [self.project.entries[i] for i in self.project.loaded()]
        targets = self.targets

        def second_pass():
            used = self.used_targets()
            self._measure(entries, used, "Measuring", finish)

        def finish():
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
                second_pass()

            self.positions = {}
            self._measure(cal, targets, "Learning positions", learned)
        else:
            self.positions = {}
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
                "levels": rel(self.levels_path), "was_run": bool(self.runquants)}

    def from_dict(self, d: dict, resolve) -> None:
        d = d or {}
        self.settings = QuantSettings.from_dict(d.get("settings", {}))
        self.decisions = QuantDecisions.from_dict(d.get("decisions", {}))
        self.report, self.runquants, self.positions = None, {}, {}
        self.stale = False
        self.levels, self.levels_path = None, None
        lv = resolve(d.get("levels"))
        if lv is not None and Path(lv).exists():
            self.levels_path = Path(lv)
            self.levels = load_levels(lv)
        self.was_run = bool(d.get("was_run"))
        self.changed.emit()
