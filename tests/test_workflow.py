"""User routes, default detection and recovery from interrupted processing."""

import csv
import json
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from pygcgc.gui.blob_state import MS_AWARE
from tests.conftest import GUI_PEAKS
from tests.test_gui import TIMEOUT, qwin, window  # noqa: F401


def test_optional_alignment_can_be_skipped_and_revisited(window):
    w = window
    assert list(w.step_actions) == ["load", "align", "background", "blobs", "istd", "assign", "quant", "export"]
    w.step_actions["align"].trigger()
    w.workflow.skip_btn.click()
    assert w.workflow.current() == "background"
    assert w.workflow.is_skipped("align")
    assert not w.workflow.is_done("align")
    w.workflow.mark_done("export")
    assert w.workflow.next_key() is None  # finishing must never send a user back to optional preparation
    w.step_actions["align"].trigger()
    assert w.align_panel.isVisibleTo(w.workflow)


def test_peak_menu_does_not_leave_duplicate_buttons_on_the_panel(window):
    w = window
    w.workflow.show_step("blobs")
    for method in (MS_AWARE, "legacy"):
        w.blobs_panel.set_method(method)
        for name in ("add_btn", "merge_btn", "delete_btn", "export_btn", "account_btn"):
            assert getattr(w.blobs_panel, name).isHidden()
    assert len(w.blobs_panel.tools_menu.actions()) == 8


@pytest.mark.parametrize("ending", ["failed_first", "failed_last", "cancelled"])
def test_incomplete_detection_never_claims_success_and_can_be_retried(window, monkeypatch, ending):
    w = window
    submitted = []
    monkeypatch.setattr(w.jobs, "submit", lambda *args, **kwargs: submitted.append((args, kwargs)))
    complete = []
    assert w.blobs.run_all(on_done=lambda: complete.append(True))
    assert w.blobs.running
    if ending == "failed_first":
        submitted[0][0][3]("bad input")
        for args, kwargs in submitted[1:]:
            args[2](None)
    elif ending == "failed_last":
        for args, kwargs in submitted[:-1]:
            args[2](None)
        submitted[-1][0][3]("bad input")
    else:
        for args, kwargs in submitted:
            kwargs["on_cancel"]()
    assert not w.blobs.running
    assert not complete
    assert not w.blobs.was_run
    assert w.blobs.run_all()  # the step is usable again after failure or cancellation
    for args, kwargs in submitted[len(submitted) // 2:]:
        kwargs["on_cancel"]()


def test_default_detector_reaches_standards_assignment_and_export(qtbot, window, tmp_path):
    w = window
    w.blobs.set_method(MS_AWARE)
    w.step_actions["background"].trigger()
    w.workflow.go_btn.click()
    qtbot.waitUntil(lambda: w.workflow.is_done("background") and w.jobs.active_count() == 0, timeout=TIMEOUT)
    w.workflow.banner_next.click()
    assert w.workflow.current() == "blobs"
    w.workflow.go_btn.click()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.candidates_were_run, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    components = w.blobs.components_for(i).components
    assert w.blobs.component_settings.model == "tilt"
    assert components and all(c.model == "tilt" for c in components)
    for t1, t2, ion, _ in GUI_PEAKS:
        nearby = [c for c in components if abs(c.t1 - t1) < 4 and abs(c.t2 - t2) < .15]
        assert nearby, f"Missing planted peak on m/z {ion}"
    component = min(components, key=lambda c: abs(c.t1 - 40) + abs(c.t2 - 1.1))
    w.workflow.banner_next.click()
    assert w.workflow.current() == "istd"
    name = w.istd.library[0].name
    w.istd.label(name, i, component.id)
    qtbot.waitUntil(lambda: name in w.istd.references, timeout=TIMEOUT)
    w.workflow.go_btn.click()
    qtbot.waitUntil(lambda: "RUN_B" in w.istd.matches, timeout=TIMEOUT)
    assert w.istd.matches["RUN_B"][name].found
    w.step_actions["assign"].trigger()
    w.workflow.go_btn.click()
    assert w.workflow.current() == "quant"
    assert w.workflow.is_done("assign")
    qtbot.waitUntil(lambda: w.jobs.active_count() == 0, timeout=TIMEOUT)
    original_path, original_sequence = w.project.path, w.project.sequence_path
    w.step_actions["export"].trigger()
    w.export_everything(tmp_path / "export")
    qtbot.waitUntil(lambda: w.workflow.is_done("export"), timeout=TIMEOUT)
    out = tmp_path / "export"
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["quantification_available"] is False
    assert len(manifest["runs"]) == 3
    assert len(list(csv.DictReader((out / "peaks_001.csv").open()))) == len(components)
    assert (out / "internal_standard_matches.csv").exists()
    assert w.project.path == original_path and w.project.sequence_path == original_sequence


def test_export_preserves_inputs_and_declares_stale_results(window, tmp_path):
    from pygcgc.gui.export import export_project
    w = window
    w.quant.levels = {"L1": {"std0": .25}, "L2": {"std0": 1.0}}
    w.quant.stale = True
    original = w.project.to_dict()
    export_project(w.project, tmp_path)
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["quantification_stale"]
    snapshot = json.loads((tmp_path / "project.pygcgc.json").read_text())
    assert snapshot["targets"] == "target_table.csv"
    assert snapshot["quant"]["levels"] == "calibration_levels.csv"
    assert w.project.to_dict() == original


def test_opening_identifier_does_not_claim_peaks_are_identified(window):
    w = window
    w.workflow.show_step("identify")
    w.workflow.go_btn.click()
    assert not w.workflow.is_done("identify")
    assert not w.workflow.is_done("assign")


def test_uncertain_internal_standard_fit_uses_the_review_colour():
    from pygcgc.core.istd import FOUND, IstdMatch
    from pygcgc.gui.panels.istd import _colour_key
    match = IstdMatch(run="sample", name="IS", status=FOUND, n_passing=1,
                      spectral_score=.999, peak_confidence="low", peak_unresolved=True)
    assert _colour_key(match) == "ambiguous"


def test_loading_assignments_preserves_the_prepared_coordinates(window, monkeypatch):
    from PySide6 import QtWidgets
    w = window
    phases = w.project.phase_snapshot()
    called = []
    monkeypatch.setattr(w, "fix_phase", lambda: called.append(True))
    monkeypatch.setattr(QtWidgets.QFileDialog, "getOpenFileName",
                        lambda *args, **kwargs: (str(w.project.targets_path), "CSV (*.csv)"))
    w.assign_step.buttons["targets"].click()
    assert not called
    assert w.project.phase_snapshot() == phases


@pytest.mark.parametrize("ending", ["failed", "cancelled"])
def test_quantification_failure_does_not_finish_a_partial_sequence(window, monkeypatch, ending):
    w = window
    submitted = []
    monkeypatch.setattr(w.jobs, "submit", lambda *args, **kwargs: submitted.append((args, kwargs)))
    complete = []
    assert w.quant.run(on_done=lambda: complete.append(True))
    for args, kwargs in submitted:
        if ending == "failed":
            args[3]("could not measure")
        else:
            kwargs["on_cancel"]()
    assert not w.quant.running
    assert not complete
    assert w.quant.report is None
    assert w.quant.stale


def test_export_quantified_amounts_round_trips_inputs(qwin, tmp_path):
    from pygcgc.gui.export import export_project
    from pygcgc.core.targets import load_levels, load_targets
    from pygcgc.core.istd import FOUND, IstdMatch
    w = qwin
    w.istd.matches = {"S1": {"IS": IstdMatch(run="S1", name="IS", status=FOUND)}}
    paths = export_project(w.project, tmp_path)
    assert len(paths) == len(set(paths))
    manifest = json.loads((tmp_path / "manifest.json").read_text())
    assert manifest["quantification_available"] and not manifest["quantification_stale"]
    assert load_levels(tmp_path / "calibration_levels.csv") == w.quant.levels
    assert load_targets(tmp_path / "target_table.csv") == w.quant.targets
    comparison = w.quant.report.write(tmp_path / "comparison")
    for path in comparison:
        assert path.name in manifest["files"]
        assert path.read_bytes() == (tmp_path / path.name).read_bytes()
    assert (tmp_path / "internal_standard_matches.csv").exists()


def test_complete_workflow_recovers_known_sample_amounts(qtbot, cal_folder, monkeypatch, tmp_path):
    from dataclasses import replace
    from PySide6 import QtWidgets
    from pygcgc.gui.mainwindow import MainWindow
    from tests.conftest import SAMPLE_AMOUNTS

    w = MainWindow()
    qtbot.addWidget(w)
    w.show()
    assert w.blobs.method == MS_AWARE
    w.blobs.set_component_settings(replace(w.blobs.component_settings, workers=1))
    monkeypatch.setattr(QtWidgets.QFileDialog, "getExistingDirectory",
                        lambda *args, **kwargs: str(cal_folder))
    try:
        w.workflow.go_btn.click()  # Load
        qtbot.waitUntil(lambda: len(w.project.loaded()) == 7 and w.jobs.active_count() == 0, timeout=TIMEOUT)
        w.workflow.banner_next.click()
        assert w.workflow.current() == "align"
        w.workflow.skip_btn.click()
        assert w.workflow.current() == "background"
        w.workflow.go_btn.click()
        qtbot.waitUntil(lambda: w.workflow.is_done("background") and w.jobs.active_count() == 0, timeout=TIMEOUT)
        w.workflow.banner_next.click()
        w.workflow.go_btn.click()  # Peaks
        qtbot.waitUntil(lambda: w.blobs.candidates_were_run and not w.blobs.running, timeout=TIMEOUT)
        w.workflow.banner_next.click()
        assert w.workflow.current() == "istd"
        i = w.project.index_of("CAL_L3")
        component = min(w.blobs.components_for(i).components,
                        key=lambda c: abs(c.t1 - 60) + abs(c.t2 - 1.2))
        w.istd.label("IS", i, component.id)
        qtbot.waitUntil(lambda: "IS" in w.istd.references, timeout=TIMEOUT)
        w.workflow.go_btn.click()
        qtbot.waitUntil(lambda: len(w.istd.matches) == 7 and w.jobs.active_count() == 0, timeout=TIMEOUT)
        assert all(row["IS"].found for row in w.istd.matches.values())
        w.step_actions["assign"].trigger()
        w.workflow.go_btn.click()  # The folder supplied the named target table
        assert w.workflow.current() == "quant"
        w.workflow.go_btn.click()
        qtbot.waitUntil(lambda: w.workflow.is_done("quant") and w.jobs.active_count() == 0, timeout=TIMEOUT)
        for row in w.quant.report.results:
            if row.run in SAMPLE_AMOUNTS and row.compound in SAMPLE_AMOUNTS[row.run]:
                assert row.amount == pytest.approx(SAMPLE_AMOUNTS[row.run][row.compound], rel=.15)
        w.workflow.banner_next.click()
        assert w.workflow.current() == "export"
        destination = tmp_path / "complete_export"
        monkeypatch.setattr(QtWidgets.QFileDialog, "getExistingDirectory",
                            lambda *args, **kwargs: str(destination))
        w.workflow.go_btn.click()
        qtbot.waitUntil(lambda: w.workflow.is_done("export"), timeout=TIMEOUT)
        assert (destination / "targets.csv").exists()
        assert (destination / "internal_standard_matches.csv").exists()
        assert (destination / "internal_standards.csv").exists()
    finally:
        w.jobs.cancel_all()
        w.jobs.wait(5000)
        w.undo.clear()
