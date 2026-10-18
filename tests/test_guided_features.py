"""Exercise guided input, exception review, prerequisites and isolated previews as user routes."""

import json
import os
from dataclasses import replace

import numpy as np
import pytest

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
pytest.importorskip('PySide6')
from PySide6 import QtCore, QtWidgets

from pygcgc.core.istd import IstdReference
from pygcgc.core.targets import Target
from pygcgc.gui.blob_state import MS_AWARE
from tests.test_gui import TIMEOUT, window  # noqa: F401


def set_cell(table, row, column, value):
    table.setItem(row, column, QtWidgets.QTableWidgetItem(str(value)))


def guide_targets(w):
    w.quant.targets = [Target('IS', role='istd', quant_ion=50, t1=40, t2=1.1),
                       Target('A', quant_ion=51, t1=90, t2=1.9, istd='IS')]
    w.quant.targets_error = ''


def configured_wizard(qtbot, w):
    guide_targets(w)
    # This fixture has two calibration runs; use the explicit two-point policy
    # for input/persistence checks without changing the production three-level default.
    w.quant.set_settings(replace(w.quant.settings, calibration=replace(w.quant.settings.calibration, min_levels=2)))
    w.show_calibration_setup()
    wizard = w.calibration_setup_wizard
    qtbot.addWidget(wizard)
    assert wizard.currentId() == 0
    assert wizard.button(QtWidgets.QWizard.WizardButton.NextButton).isVisibleTo(wizard)
    button = wizard.button(QtWidgets.QWizard.WizardButton.NextButton)
    assert wizard.rect().contains(button.mapTo(wizard, button.rect().center()))
    button.click()
    assert wizard.currentId() == 1
    wizard.runs_table.cellWidget(1, 1).setCurrentText('calibration')
    set_cell(wizard.runs_table, 1, 2, 'L2')
    button.click()
    assert wizard.currentId() == 2
    assert not wizard._amount_widgets['IS'].text()  # no guessed one-ng spike
    set_cell(wizard.amounts_table, 0, 2, 1000)
    set_cell(wizard.amounts_table, 1, 2, 4000)
    wizard._amount_widgets['IS'].setText('1000')
    return wizard


def test_guided_setup_is_atomic_converts_units_and_survives_project_reload(qtbot, window, tmp_path):
    w = window
    originals = {p: p.read_bytes() for p in (w.project.targets_path, w.project.sequence_path)}
    wizard = configured_wizard(qtbot, w)
    old_levels = w.quant.levels
    wizard.unit_box.setCurrentText('ug')
    assert wizard.draft(True).levels['L1']['A'] == pytest.approx(1)
    assert wizard.draft(True).internal_standard_amounts['IS'] == pytest.approx(1)
    assert w.quant.levels is old_levels  # changing a draft must not mutate active inputs
    wizard.next()
    assert wizard.currentId() == 3
    assert wizard.button(QtWidgets.QWizard.WizardButton.FinishButton).isVisibleTo(wizard)
    wizard.accept()
    assert wizard.isVisible()  # confirmation is required before the application is changed
    wizard.confirm.setChecked(True)
    wizard.accept()
    assert w.quant.amount_unit == 'ug'
    assert w.quant.levels['L2']['A'] == pytest.approx(4)
    assert w.project.entries[1].kind == 'calibration'
    assert all(path.read_bytes() == before for path, before in originals.items())
    out = tmp_path / 'guided.pygcgc.json'
    w.project.save(out)
    from pygcgc.gui.mainwindow import MainWindow
    restored = MainWindow()
    qtbot.addWidget(restored)
    restored.show()
    try:
        loaded = []
        restored.project.load(out, on_loaded=lambda: loaded.append(True))
        qtbot.waitUntil(lambda: bool(loaded) and restored.jobs.active_count() == 0, timeout=TIMEOUT)
        assert restored.quant.amount_unit == 'ug'
        assert [t.compound for t in restored.quant.targets] == ['IS', 'A']
        assert restored.quant.levels == w.quant.levels
        assert restored.quant.calibration_setup['confirmed']
        from pygcgc.gui.export import export_project
        destination = tmp_path / 'export'
        export_project(restored.project, destination)
        manifest = json.loads((destination / 'manifest.json').read_text())
        assert manifest['amount_unit'] == 'ug'
        assert json.loads((destination / 'calibration_setup.json').read_text())['internal_standard_amounts']['IS'] == 1
    finally:
        restored.jobs.cancel_all()
        restored.jobs.wait(5000)
        restored.undo.clear()


def test_cancelling_setup_preserves_inputs_and_run_metadata(qtbot, window):
    wizard = configured_wizard(qtbot, window)
    before = window.project.to_dict()
    wizard.reject()
    assert window.project.to_dict() == before


def test_typing_a_mass_unit_converts_only_the_completed_unit_and_unknown_is_stays_unknown(qtbot, window):
    wizard = configured_wizard(qtbot, window)
    edit = wizard.unit_box.lineEdit()
    edit.setFocus()
    edit.selectAll()
    qtbot.keyClicks(edit, 'ug')
    assert wizard.draft(True).levels['L1']['A'] == pytest.approx(1)
    assert wizard.draft(True).internal_standard_amounts['IS'] == pytest.approx(1)
    wizard._amount_widgets['IS'].clear()
    wizard.next()
    assert wizard.currentId() == 2
    assert 'IS: enter' in wizard.amount_error.text()
    wizard.reject()


def test_prerequisite_actions_and_invalidation_are_visible(qtbot, window):
    w = window
    w.workflow.show_step('istd')
    w.workflow_assistant.refresh()
    assert not w.workflow.go_btn.isEnabled()
    assert 'Peaks' in w.workflow.go_note.text()
    assert w.workflow.prerequisite_btn.isVisibleTo(w.workflow)
    w.workflow.prerequisite_btn.click()
    assert w.workflow.current() == 'blobs'
    w.workflow.show_step('quant')
    w.workflow_assistant.refresh()
    assert not w.workflow.go_btn.isEnabled()
    assert 'calibration' in w.workflow.go_note.text().lower() or 'named targets' in w.workflow.go_note.text().lower()
    w.workflow.prerequisite_btn.click()
    assert w.calibration_setup_wizard.isVisible()
    w.calibration_setup_wizard.reject()
    for key in ('blobs', 'istd', 'assign', 'quant', 'export'):
        w.workflow.mark_done(key)
    w.project.foldChanged.emit(0)
    assert all(not w.workflow.is_done(key) for key in ('blobs', 'istd', 'assign', 'quant', 'export'))


def prepare_standard(qtbot, w, ms_aware=False):
    if ms_aware:
        w.blobs.set_method(MS_AWARE)
        w.blobs.set_component_settings(replace(w.blobs.component_settings, workers=1))
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.jobs.active_count() == 0, timeout=TIMEOUT)
    index = w.project.index_of('RUN_A')
    peaks = w.blobs.components_for(index).components if ms_aware else w.blobs.blobs_for(index)
    if ms_aware:
        component = min(peaks, key=lambda c: abs(c.t1-40)+abs(c.t2-1.1))
        peak_id = component.id
    else:
        peak_id = int(np.argmin(abs(peaks.t1-40)+abs(peaks.t2-1.1)))
    w.istd.label('IS', index, peak_id)
    qtbot.waitUntil(lambda: 'IS' in w.istd.references, timeout=TIMEOUT)
    w.istd.propagate_to(w.project.loaded())
    qtbot.waitUntil(lambda: not w.istd.running and len(w.istd.matches) == 3, timeout=TIMEOUT)
    return index


def test_review_decision_is_per_run_persists_and_invalidates_after_threshold_change(qtbot, window):
    w = window
    prepare_standard(qtbot, w)
    reference = w.istd.references['IS']
    original = reference.spectrum.copy()
    match = w.istd.matches['RUN_A']['IS']
    w.istd.matches['RUN_A']['IS'] = replace(match, n_passing=2, peak_unresolved=True)
    assert ('RUN_A', 'IS') in w.istd.exceptions()
    w.istd.review_match('RUN_A', 'IS')
    assert ('RUN_A', 'IS') not in w.istd.exceptions()
    assert w.istd.matches['RUN_A']['IS'].peak_unresolved  # acknowledgement never removes scientific warnings
    saved = json.loads(json.dumps(w.istd.to_dict()))
    before_signature = w.istd.input_signature(w.project.index_of('RUN_A'))
    before_key = w.istd._signature_cache['RUN_A'][0]
    w.istd.from_dict(saved)
    after_signature = w.istd.input_signature(w.project.index_of('RUN_A'))
    assert before_signature == after_signature, (before_key, w.istd._signature_cache['RUN_A'][0])
    assert w.istd.reviewed('RUN_A', 'IS'), (saved['review_decisions'], w.istd.matches['RUN_A']['IS'])
    np.testing.assert_array_equal(w.istd.references['IS'].spectrum, original)
    w.istd.set_settings(replace(w.istd.settings, min_spectral=.81))
    assert w.istd.is_stale('RUN_A')
    assert not w.istd.reviewed('RUN_A', 'IS')
    with pytest.raises(ValueError, match='out of date'):
        w.istd.review_match('RUN_A', 'IS')


def test_exceptions_queue_loads_candidate_spectra_and_real_fitted_profiles(qtbot, window):
    w = window
    prepare_standard(qtbot, w, ms_aware=True)
    match = w.istd.matches['RUN_A']['IS']
    w.istd.matches['RUN_A']['IS'] = replace(match, n_passing=2)
    w.show_standard_review()
    dialog = w.standard_review
    qtbot.waitUntil(lambda: dialog._data is not None and bool(dialog._candidates), timeout=TIMEOUT)
    assert dialog._pair == ('RUN_A', 'IS')
    assert dialog._candidates[0]['profile'] is not None
    assert dialog.spectra[0].listDataItems()
    assert dialog.profiles[0].listDataItems()
    w.istd.running = True
    w.istd.changed.emit()
    assert not dialog.accept_btn.isEnabled() and not dialog.missing_btn.isEnabled()
    with pytest.raises(ValueError, match='search to finish'):
        w.istd.review_match('RUN_A', 'IS')
    w.istd.running = False
    w.istd.changed.emit()
    before = len(w.istd.exceptions())
    dialog.ack_btn.click()
    assert len(w.istd.exceptions()) == before-1
    dialog.close()
    assert dialog._job is None


def test_preview_uses_current_run_without_mutating_accepted_results(qtbot, window):
    w = window
    w.blobs.set_method(MS_AWARE)
    before = w.blobs.to_dict()
    settings = w.blobs.component_settings
    w.view2d.plot.vb.setRange(xRange=(.3, 1.5), yRange=(0, 4), padding=0)
    w.show_peak_preview()
    dialog = w.peak_preview
    dialog.run_btn.click()
    qtbot.waitUntil(lambda: dialog.apply_btn.isEnabled() and w.jobs.active_count() == 0, timeout=TIMEOUT)
    assert dialog.map._data is not None
    assert np.isfinite(dialog.map._data['model']).all()
    assert w.blobs.to_dict() == before and w.blobs.component_settings == settings
    w.blobs.set_component_settings(replace(settings, fast_max_nfev=80))
    assert not dialog.apply_btn.isEnabled()
    assert 'changed' in dialog.status.text().lower()
    dialog.close()


def test_closing_or_cancelling_preview_discards_late_results(qtbot, window, monkeypatch):
    w = window
    w.show_peak_preview()
    dialog = w.peak_preview
    submitted = []
    monkeypatch.setattr(w.jobs, 'submit', lambda *args, **kw: submitted.append((args, kw)) or 987)
    monkeypatch.setattr(w.jobs, 'cancel', lambda job: None)
    dialog.run_btn.click()
    assert dialog._job == 987
    dialog.cancel_btn.click()
    assert dialog._job is None and not dialog.apply_btn.isEnabled()
    submitted[0][0][2]({})  # late callback is ignored before accessing discarded data
    assert dialog.map._data is None
    dialog.run_btn.click()
    dialog.close()
    submitted[1][0][2]({})
    assert dialog.map._data is None


def test_reviewed_positions_and_confirmed_missing_standards_reach_quantification(qtbot, window, monkeypatch):
    w = window
    prepare_standard(qtbot, w)
    guide_targets(w)
    w.quant.set_settings(replace(w.quant.settings, learn_positions=False))
    chosen = w.istd.matches['RUN_A']['IS']
    w.istd.review_match('RUN_A', 'IS', chosen)
    w.istd.review_match('RUN_B', 'IS', missing=True)
    assert w.istd.reviewed('RUN_B', 'IS')
    import pygcgc.gui.quantstate as module
    original = module.quantify_run
    assigned_inputs = []
    def record(*args, **kwargs):
        assigned_inputs.append(dict(args[6]))
        return original(*args, **kwargs)
    monkeypatch.setattr(module, 'quantify_run', record)
    assert w.quant.run()
    qtbot.waitUntil(lambda: not w.quant.running and w.jobs.active_count() == 0, timeout=TIMEOUT)
    assert any(values.get('IS') == (chosen.t1, chosen.t2) for values in assigned_inputs)
    result = w.quant.run_quant(w.project.index_of('RUN_B'))
    assert not result.hits['IS'].found
    assert np.isnan(result.response('A'))
    assert all(np.isnan(row.amount) for row in w.quant.report.results if row.run == 'RUN_B')


def test_confirmed_missing_with_no_candidate_survives_serialization(qtbot, window):
    w = window
    prepare_standard(qtbot, w)
    w.istd.references['Absent'] = IstdReference('Absent', '', 'RUN_A', 1000., 1., 50, np.array([1.]))
    w.istd.propagate_to(w.project.loaded())
    qtbot.waitUntil(lambda: not w.istd.running, timeout=TIMEOUT)
    assert np.isnan(w.istd.matches['RUN_A']['Absent'].t1)
    w.istd.review_match('RUN_A', 'Absent', missing=True)
    w.istd.from_dict(json.loads(json.dumps(w.istd.to_dict())))
    assert w.istd.reviewed('RUN_A', 'Absent')


def test_quantification_does_not_publish_success_if_inputs_change_during_jobs(window, monkeypatch):
    w = window
    guide_targets(w)
    w.quant.set_settings(replace(w.quant.settings, learn_positions=False))
    measured = []
    monkeypatch.setattr(w.quant, '_measure', lambda entries, targets, label, on_done: measured.append(on_done))
    finished = []
    assert w.quant.run(on_done=lambda: finished.append(True))
    w.quant.mark_stale()
    measured[0]()
    assert not w.quant.running and w.quant.stale
    assert not finished and w.quant.report is None


def test_cancelled_standard_search_resets_running_without_success(qtbot, window, monkeypatch):
    w = window
    prepare_standard(qtbot, w)
    submitted = []
    monkeypatch.setattr(w.jobs, 'submit', lambda *args, **kwargs: submitted.append((args, kwargs)))
    complete = []
    w.istd.propagate_to(w.project.loaded(), on_done=lambda: complete.append(True))
    assert w.istd.running
    for args, kwargs in submitted:
        kwargs['on_cancel']()
    assert not w.istd.running and not complete


def test_guided_units_produce_known_quantities_and_exported_unit_metadata(qtbot, cal_folder, tmp_path):
    from pygcgc.gui.mainwindow import MainWindow
    from pygcgc.gui.export import export_project
    from tests.conftest import SAMPLE_AMOUNTS
    w = MainWindow()
    qtbot.addWidget(w)
    w.show()
    try:
        loaded = []
        w.open_folder(cal_folder, on_loaded=lambda: loaded.append(True))
        qtbot.waitUntil(lambda: bool(loaded) and w.jobs.active_count() == 0, timeout=TIMEOUT)
        w.show_calibration_setup()
        wizard = w.calibration_setup_wizard
        wizard.next()
        wizard.next()
        assert wizard.currentId() == 2
        wizard._amount_widgets['IS'].setText('1')
        wizard.unit_box.setCurrentText('ug')
        wizard.next()
        assert wizard.currentId() == 3
        wizard.confirm.setChecked(True)
        wizard.accept()
        assert w.quant.amount_unit == 'ug'
        w.workflow.go_btn.click()
        qtbot.waitUntil(lambda: not w.quant.running and w.quant.report is not None
                       and w.jobs.active_count() == 0, timeout=TIMEOUT)
        for row in w.quant.report.results:
            if row.run in SAMPLE_AMOUNTS and row.compound in SAMPLE_AMOUNTS[row.run]:
                assert row.amount == pytest.approx(SAMPLE_AMOUNTS[row.run][row.compound]/1000, rel=.15)
        destination = tmp_path / 'units_export'
        export_project(w.project, destination)
        assert json.loads((destination/'manifest.json').read_text())['amount_unit'] == 'ug'
        assert 'ug' in w.results.model.headerData(w.results.model.keys.index('amount'), QtCore.Qt.Orientation.Horizontal)
    finally:
        w.jobs.cancel_all()
        w.jobs.wait(5000)
        w.undo.clear()


def test_theme_switch_is_safe_with_open_and_closed_review_graphs(qtbot, window):
    from pathlib import Path
    from PySide6 import QtGui
    from pygcgc.gui.theme import apply_theme
    for name in ('arial.ttf', 'arialbd.ttf'):
        path = Path('C:/Windows/Fonts')/name
        if path.exists():
            QtGui.QFontDatabase.addApplicationFont(str(path))
    w = window
    prepare_standard(qtbot, w, ms_aware=True)
    w.show_standard_review()
    dialog = w.standard_review
    dialog.only_exceptions.setChecked(False)
    qtbot.waitUntil(lambda: dialog._data is not None and bool(dialog._candidates), timeout=TIMEOUT)
    apply_theme(w, False)
    qtbot.wait(50)
    apply_theme(w, True)
    qtbot.wait(50)
    dialog.close()
    apply_theme(w, False)
    qtbot.wait(50)
    assert not dialog.isVisible() and not dialog._theme_timer.isActive()


def test_applying_a_fitting_preview_explicitly_starts_the_full_route(qtbot, window, monkeypatch):
    w = window
    w.blobs.set_method(MS_AWARE)
    requested = []
    monkeypatch.setattr(w, 'run_blobs', lambda: requested.append(True))
    w.show_peak_preview()
    dialog = w.peak_preview
    dialog.settings = replace(dialog.settings, fast_max_nfev=80)
    dialog.run_btn.click()
    qtbot.waitUntil(lambda: dialog.apply_btn.isEnabled() and w.jobs.active_count() == 0, timeout=TIMEOUT)
    assert not requested and w.blobs.component_settings.fast_max_nfev == 60
    dialog.apply_btn.click()
    assert requested == [True]
    assert w.blobs.component_settings.fast_max_nfev == 80
    assert not dialog.isVisible()


def test_standard_review_remains_valid_after_full_project_reload(qtbot, window, tmp_path):
    from pygcgc.gui.mainwindow import MainWindow
    w = window
    prepare_standard(qtbot, w)
    w.istd.review_match('RUN_A', 'IS', w.istd.matches['RUN_A']['IS'])
    path = tmp_path/'reviewed.pygcgc.json'
    w.project.save(path)
    restored = MainWindow()
    qtbot.addWidget(restored)
    restored.show()
    try:
        loaded = []
        restored.project.load(path, on_loaded=lambda: loaded.append(True))
        qtbot.waitUntil(lambda: bool(loaded) and restored.jobs.active_count() == 0, timeout=TIMEOUT)
        assert restored.istd.reviewed('RUN_A', 'IS')
        assert restored.istd.matches['RUN_A']['IS'].peak_id == w.istd.matches['RUN_A']['IS'].peak_id
    finally:
        restored.jobs.cancel_all()
        restored.jobs.wait(5000)
        restored.undo.clear()


@pytest.mark.parametrize('ms_aware', [True, False])
def test_one_selected_peak_can_be_added_to_the_calibration_draft(qtbot, window, ms_aware):
    w = window
    index = prepare_standard(qtbot, w, ms_aware=ms_aware)
    w.set_run(index)
    if ms_aware:
        peak = w.blobs.components_for(index).components[0]
        w.blobs_panel.select_component(peak.id)
        expected = (peak.t1, peak.t2)
    else:
        w.blobs_panel.select_row(0)
        peaks = w.blobs.blobs_for(index)
        expected = (peaks.t1[0], peaks.t2[0])
    w.show_calibration_setup()
    wizard = w.calibration_setup_wizard
    assert wizard.selected_peak_btn.isEnabled()
    before = wizard.targets_table.rowCount()
    wizard.selected_peak_btn.click()
    assert wizard.targets_table.rowCount() == before+1
    selected = wizard.targets()[-1]
    assert (selected.t1, selected.t2) == pytest.approx(expected)
    assert selected.compound == '' and selected.quant_ion > 0
    wizard.reject()
