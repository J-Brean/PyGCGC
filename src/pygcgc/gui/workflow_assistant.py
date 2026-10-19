"""Readiness and freshness guidance derived from the actual project state."""

from __future__ import annotations

import math

from PySide6 import QtCore

from .blob_state import MS_AWARE
from ..core.calibration_inputs import calibration_coverage_errors


def calibration_problems(project):
    q = project.quant
    if not q.targets or not any(not t.is_istd for t in q.targets):
        return ["Add named targets and their internal standards in calibration setup."]
    if not q.levels:
        return ["Enter known calibration amounts and their unit in calibration setup."]
    sets = {}
    for i in project.loaded():
        e = project.entries[i]
        if 'cal' in e.kind.lower():
            if not e.level or e.level not in q.levels:
                return [f"{e.name}: choose a calibration level with known amounts."]
            sets.setdefault(e.cal_set, set()).add(e.level)
    if not sets:
        return [f"Choose calibration runs and known levels in calibration setup "
                f"(current default minimum: {max(2, q.settings.calibration.min_levels)} usable levels)."]
    for cal_set, levels in sets.items():
        for target in q.targets:
            if target.is_istd:
                continue
            values = [q.levels.get(lv, {}).get(target.compound) for lv in levels]
            valid = [v for v in values if v is not None and math.isfinite(v) and v >= 0]
            if len(valid) != len(values) or len(set(valid)) < 2:
                return [f"{target.compound}, set {cal_set or '(default)'}: supply at least two distinct "
                        "known amounts covering every calibration level."]
    runs = {project.entries[i].name: dict(kind=project.entries[i].kind, level=project.entries[i].level,
                                        cal_set=project.entries[i].cal_set) for i in project.loaded()}
    return calibration_coverage_errors(q.targets, q.levels, runs, q.settings.calibration, q.decisions)


class WorkflowAssistant(QtCore.QObject):
    """Refresh once per event-loop turn and keep stale completion banners honest."""

    def __init__(self, window):
        super().__init__(window)
        self.window = window
        self.timer = QtCore.QTimer(self)
        self.timer.setSingleShot(True)
        self.timer.timeout.connect(self.refresh)
        p = window.project
        for signal in (p.runsChanged, p.foldChanged, p.layersChanged, p.blobs.changed,
                       p.blobs.runningChanged, p.istd.changed, p.quant.changed, p.quant.runningChanged,
                       p.align.changed, window.workflow.progressChanged):
            signal.connect(self.schedule)
        p.align.refinementRunning.connect(self.schedule)
        p.foldChanged.connect(self._fold_changed)
        # Labelling a standard must enable Go before the emitting action returns.
        # Coalesced refreshes still handle the noisier rendering/preparation signals.
        p.istd.changed.connect(self.refresh)
        self.refresh()

    def schedule(self, *_):
        self.timer.start(0)

    def _clear(self, keys):
        workflow = self.window.workflow
        for key in keys:
            if workflow.is_done(key):
                workflow.clear_done(key)
                if workflow.banner.property('for') == key:
                    workflow.banner.hide()

    def _fold_changed(self, *_):
        self._clear(('background', 'blobs', 'istd', 'assign', 'quant', 'export'))

    def refresh(self):
        w, p = self.window, self.window.project
        wf, b, q, state = w.workflow, p.blobs, p.quant, p.istd
        loaded = p.loaded()
        if not loaded:
            for key in ('fold', 'align', 'background', 'blobs', 'istd', 'assign', 'quant', 'export'):
                wf.set_go_enabled(key, False, "Load chromatograms first.")
                wf.set_guidance(key, "Load chromatograms before this step.", 'load', 'Open Load data')
            return
        for key in ('fold', 'align', 'background', 'blobs', 'assign', 'export'):
            wf.set_go_enabled(key, True)
            wf.set_guidance(key, '')
        if len(loaded) < 2 or p.align.refining:
            wf.set_go_enabled('align', False)
            wf.set_guidance('align', 'Wait for anchor refinement to finish.' if p.align.refining else
                             'Alignment needs at least two loaded chromatograms. Use Skip alignment for a single run.')
        elif not p.align.reference:
            wf.set_guidance('align', 'Start picking to use the first run as the reference, or choose one above. '
                             'Each click advances to the next chromatogram.')
        missing = []
        stale = []
        for i in loaded:
            e = p.entries[i]
            result = b.components_for(i) if b.method == MS_AWARE else b.blobs_for(i)
            if result is None:
                missing.append(e.name)
            elif b.selected_stale(e):
                stale.append(e.name)
        peaks_current = not missing and not stale and not b.running
        if missing or stale:
            self._clear(('blobs', 'istd', 'assign', 'export'))
        raw_count = sum(p.entries[i].background is None for i in loaded)
        note = (f"{len(stale)} run(s) need new peak fits after settings or preparation changed. " if stale else '')
        note += (f"{len(missing)} run(s) have no peak results. " if missing else '')
        if raw_count:
            note += f"{raw_count} run(s) have no background layer; fitting uses their raw signal."
        wf.set_guidance('blobs', note, 'background' if raw_count else '',
                        'Open Background removal' if raw_count else '')
        w.blobs_panel.preview_btn.setEnabled(not b.running)
        wf.set_busy('blobs', b.running)
        if not peaks_current:
            wf.set_go_enabled('istd', False)
            wf.set_guidance('istd', 'Refresh Peaks on every loaded run before searching for standards.',
                            'blobs', 'Open Peaks')
            wf.set_go_enabled('assign', False)
            wf.set_guidance('assign', 'Refresh Peaks before reviewing assignments.', 'blobs', 'Open Peaks')
        elif not state.references:
            wf.set_go_enabled('istd', False)
            wf.set_guidance('istd', 'Select an internal standard above, then label its peak on one reference run. '
                            'Search all runs after labelling.')
        else:
            wf.set_go_enabled('istd', not state.running)
            exceptions = state.exceptions()
            stale_matches = sum(bool(state.matches.get(p.entries[i].name)) and state.is_stale(p.entries[i].name)
                                for i in loaded)
            if stale_matches:
                self._clear(('istd', 'export'))
            wf.set_guidance('istd',
                f"{stale_matches} run(s) need a new standard search. " +
                f"{len(exceptions)} missing, uncertain or unsearched standard/run pair(s).",
                'review_istd' if exceptions and not stale_matches and not state.running else '',
                'Review standard exceptions')
        wf.set_busy('istd', state.running)
        w.istd_panel.find_btn.setEnabled(peaks_current and bool(state.references) and not state.running)
        if q.stale:
            self._clear(('quant', 'export'))
        ready, why = q.ready()
        problems = calibration_problems(p) if ready else []
        ok = ready and not problems and not q.running
        outdated_reviews = state.outdated_reviews()
        if outdated_reviews:
            ok = False
        wf.set_go_enabled('quant', ok, why)
        action = ''
        label = ''
        if outdated_reviews:
            note = 'Previously reviewed standard choices are out of date. Search standards again and review the new matches.'
            action, label = 'istd', 'Open Standards'
        elif not ready:
            action = 'fold' if 'modulation period' in why else 'load' if 'energies' in why else 'setup'
            label = {'fold': 'Check modulation period', 'load': 'Open Load data', 'setup': 'Set up calibration and units…'}[action]
            note = why
        elif problems:
            note, action, label = problems[0], 'setup', 'Set up calibration and units…'
        else:
            note = f"Amounts: {q.amount_unit}. Concentrations: {q.amount_unit}/m³ when volume is supplied."
            if q.stale:
                note = 'Quantitative results are out of date; rerun Quantify. ' + note
            exceptions = state.exceptions() if state.references else []
            if exceptions:
                note += f" {len(exceptions)} standard exception(s) remain for review."
                action, label = 'review_istd', 'Review standard exceptions'
        wf.set_guidance('quant', note, action, label)
        wf.set_busy('quant', q.running)
        w.quant_step.unit_label.setText(f"Amount unit: {q.amount_unit}; concentration: {q.amount_unit}/m³.")
        export_notes = []
        if q.report is None:
            export_notes.append('No quantitative results yet; export will contain available analysis and inputs.')
        elif q.stale:
            export_notes.append('Quantitative results are out of date and will be marked stale in the manifest.')
        if stale:
            export_notes.append(f"Peak results are out of date in {len(stale)} run(s).")
        wf.set_guidance('export', ' '.join(export_notes), 'quant' if q.report is not None and q.stale else '',
                        'Open Quantify')
