"""The ➖ Subtract a run's peaks step: settings, results, link to quantification and persistence."""

import os
from dataclasses import replace

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from pygcgc.core.peak_blank import PeakBlankSettings  # noqa: E402
from tests.test_gui import TIMEOUT, window  # noqa: E402,F401


def peaks_found(qtbot, w):
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    w.set_run(w.project.index_of("RUN_A"))


def subtracted(qtbot, w, blank="RUN_B", factor=1.0):
    peaks_found(qtbot, w)
    state = w.project.blank
    state.set_settings(PeakBlankSettings(blank_run=blank, factor=factor))
    w.workflow.show_step("blank")
    w.workflow.go_btn.click()
    i = w.project.index_of("RUN_A")
    qtbot.waitUntil(lambda: state.current_for(i) is not None, timeout=TIMEOUT)
    return state, i


def test_the_step_sits_under_assign_and_lists_the_runs(qtbot, window):
    w = window
    info = w.workflow.step_info("blank")
    assert info["parent"] == "assign" and "peaks" in info["title"].lower()
    w.blank_panel.refresh_runs(w.project.entries)
    assert w.blank_panel.blank_box.count() == 1 + len(w.project.loaded())
    assert "blank" in w.blank_panel.blank_box.itemText(0).lower()


def test_nothing_runs_until_a_blank_is_chosen_and_peaks_exist(qtbot, window):
    w = window
    state = w.project.blank
    assert state.why_not(0) == ["Choose the blank run.", "Find peaks in this run first (step 4)."] or "Choose the blank run." in state.why_not(0)
    w.run_blank_subtraction()
    assert "Choose the blank run" in w.statusBar().currentMessage()
    state.set_settings(PeakBlankSettings(blank_run="RUN_B"))
    assert any("Find peaks" in m for m in state.why_not(w.project.index_of("RUN_A")))
    assert state.subtract_runs() is False


def test_the_blank_run_is_subtracted_from_the_others_peak_by_peak(qtbot, window):
    w = window
    state, i = subtracted(qtbot, w)
    found = state.results_for(i)
    peaks = w.blobs.full_spectrum_peaks(i)
    assert len(found.rows) == len(peaks) > 0
    assert found.summary["matched"] > 0                                  # the same compounds are in both runs
    for row in found.rows:
        assert 0.0 <= row.corrected <= row.volume + 1e-9 and row.removed == pytest.approx(row.volume - row.corrected)
    assert state.results_for(w.project.index_of("RUN_B")) is None        # the blank itself is not subtracted from
    assert w.workflow.is_done("blank")


def test_the_table_shows_each_peak_against_the_blank_in_positive_terms(qtbot, window):
    w = window
    state, i = subtracted(qtbot, w)
    panel = w.blank_panel
    panel.set_shown_run(i)
    panel.refresh()
    assert panel.table.rowCount() == len(state.results_for(i).rows)
    words = {panel.table.item(r, 9).text() for r in range(panel.table.rowCount())}
    assert all(("above blank" in x) or ("near blank" in x) for x in words)
    volumes = [float(panel.table.item(r, 3).text()) for r in range(panel.table.rowCount())]
    assert volumes == sorted(volumes, reverse=True)
    assert "peaks:" in panel.status.text() and "of the run's signal removed" in panel.status.text()
    assert "Shared peaks" in panel.suggestion.text() and panel.use_suggestion.isEnabled()


def test_using_the_suggestion_sets_the_factor_but_only_when_asked(qtbot, window):
    w = window
    state, i = subtracted(qtbot, w, factor=1.0)
    panel = w.blank_panel
    panel.set_shown_run(i)
    panel.refresh()
    assert panel.factor.value() == 1.0 and state.settings.factor == 1.0
    median = state.results_for(i).suggestion["median"]
    panel.use_suggestion.click()
    assert panel.factor.value() == pytest.approx(round(median, 3))
    assert state.settings.factor == pytest.approx(round(median, 3))        # the panel pushes the new factor to the state


def test_changing_the_factor_makes_the_results_out_of_date_not_silently_wrong(qtbot, window):
    w = window
    state, i = subtracted(qtbot, w)
    state.set_settings(replace(state.settings, factor=0.5))
    assert state.current_for(i) is None and state.results_for(i) is not None
    w.blank_panel.set_shown_run(i)
    w.blank_panel.refresh()
    assert "out of date" in w.blank_panel.status.text()
    state.subtract_runs([i])
    qtbot.waitUntil(lambda: state.current_for(i) is not None, timeout=TIMEOUT)


def test_the_blank_and_factor_are_shared_with_quantification_while_linked(qtbot, window):
    w = window
    state = w.project.blank
    state.set_settings(PeakBlankSettings(blank_run="RUN_B", factor=0.7, max_dt1=9.0))
    q = w.project.quant.settings
    assert q.blank_run == "RUN_B" and q.blank_factor == 0.7 and q.blank_max_dt1 == 9.0
    state.set_settings(state.settings, apply_to_quant=False)
    assert w.project.quant.settings.blank_run == ""
    state.set_settings(state.settings, apply_to_quant=True)
    assert w.project.quant.settings.blank_run == "RUN_B"


def test_the_panel_follows_settings_changed_elsewhere_such_as_a_loaded_project(qtbot, window):
    w = window
    w.blank_panel.refresh_runs(w.project.entries)
    state = w.project.blank
    state.from_dict({"settings": {"blank_run": "RUN_B", "factor": 0.35, "max_dt1": 120.0}, "apply_to_quant": False})
    panel = w.blank_panel
    assert panel.blank_box.currentData() == "RUN_B" and panel.factor.value() == pytest.approx(0.35)
    assert panel.dt1.value() == 2.0 and not panel.apply_quant.isChecked()
    assert state.settings.factor == 0.35                    # syncing the widgets did not write anything back


def test_the_position_limits_are_set_in_minutes_and_seconds_and_turn_away_distant_peaks(qtbot, window):
    w = window
    state, i = subtracted(qtbot, w)
    panel = w.blank_panel
    panel.set_shown_run(i)
    panel.refresh()
    assert panel.dt1.isVisibleTo(panel) and panel.dt2.isVisibleTo(panel)        # always on show, not tucked away
    assert state.settings.max_dt1 == panel.dt1.value() * 60.0
    loose = state.results_for(i).summary["matched"]
    assert loose > 0
    panel.dt1.setValue(0.01)                                                    # 0.6 s in t1
    panel.dt2.setValue(0.01)                                                    # 0.01 s in t2
    assert state.settings.max_dt1 == pytest.approx(0.6) and state.settings.max_dt2 == pytest.approx(0.01)
    assert state.current_for(i) is None                                         # the old pairing is out of date
    state.subtract_runs([i])
    qtbot.waitUntil(lambda: state.current_for(i) is not None, timeout=TIMEOUT)
    tight = state.results_for(i)
    assert tight.summary["matched"] < loose
    assert tight.summary["too_far"] > 0 and "beyond the position limits" in w.blank_panel.status.text()
    assert any("outside the limits, not used" in panel.table.item(r, 9).text() for r in range(panel.table.rowCount()))
    for row in tight.rows:                                                      # a peak turned away loses nothing
        if row.blank_index is None:
            assert row.corrected == max(row.volume, 0.0) and row.removed == 0.0


def test_settings_are_saved_and_the_results_are_recomputed(qtbot, window):
    w = window
    state = w.project.blank
    state.set_settings(PeakBlankSettings(blank_run="RUN_B", factor=0.4, significance=5.0), apply_to_quant=False)
    saved = w.project.to_dict()["peak_blank"]
    assert saved["settings"]["blank_run"] == "RUN_B" and saved["apply_to_quant"] is False
    state.from_dict(saved)
    assert state.settings.factor == 0.4 and state.results == {} and state.apply_to_quant is False


def test_the_main_export_writes_the_subtraction_and_its_settings(qtbot, window, tmp_path):
    import csv
    import json
    from pygcgc.gui.export import export_project
    w = window
    state, i = subtracted(qtbot, w)
    out = tmp_path / "export"
    export_project(w.project, out)
    rows = list(csv.DictReader((out / f"peak_blank_{i + 1:03d}.csv").open(encoding="utf-8")))
    assert len(rows) == len(state.results_for(i).rows) and rows[0]["blank_run"] == "RUN_B"
    assert {r["against_blank"] for r in rows} <= {"above blank", "near blank"}
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["runs"][i]["peak_blank"].startswith("peak_blank_") and manifest["runs"][i]["peak_blank_stale"] is False
    provenance = json.loads((out / "provenance.json").read_text())
    assert provenance["settings"]["peak_blank"]["blank_run"] == "RUN_B"
    assert provenance["settings"]["peak_blank"]["applied_to_quantification"] is True


def test_different_phases_without_an_alignment_are_called_out_because_peaks_would_not_pair(qtbot, window):
    w = window
    state = w.project.blank
    a, b = w.project.index_of("RUN_A"), w.project.index_of("RUN_B")
    w.project.set_phase(a, 0.2)
    w.project.set_phase(b, 0.2)
    state.set_settings(PeakBlankSettings(blank_run="RUN_B"))
    assert state.phase_note(a) == ""
    w.project.set_phase(b, 0.9)
    note = state.phase_note(a)
    assert "differ by 0.70 s" in note and "Use for all runs" in note
    w.blank_panel.set_shown_run(a)
    w.blank_panel.refresh()
    assert "differ by 0.70 s" in w.blank_panel.status.text()
    assert state.phase_note(b) == ""                                           # the blank itself needs no warning
    # an alignment of both runs absorbs the difference
    w.set_reference("RUN_A")
    w.align.add_anchor("RUN_B", "P1", 100.0, 1.0, 100.0, 0.3, 50)
    assert state.phase_note(a) == ""


def test_removing_the_blank_run_clears_the_choice_everywhere(qtbot, window):
    w = window
    state = w.project.blank
    state.set_settings(PeakBlankSettings(blank_run="RUN_C_BLANK"))
    w.remove_runs([w.project.index_of("RUN_C_BLANK")])
    assert state.settings.blank_run == "" and w.project.quant.settings.blank_run == ""


def test_results_are_discarded_if_the_alignment_changes_while_subtracting(qtbot, window, monkeypatch):
    from pygcgc.gui.jobs import CancelToken
    w = window
    peaks_found(qtbot, w)
    state = w.project.blank
    state.set_settings(PeakBlankSettings(blank_run="RUN_B"))
    submitted = []
    monkeypatch.setattr(w.jobs, "submit", lambda *a, **k: submitted.append((a, k)))
    assert state.subtract_runs() and state.running
    args, _kwargs = submitted[-1]
    result = args[1](CancelToken(), lambda *_: None)
    w.project.period = w.project.period + 0.5
    args[2](result)
    assert not state.running and state.results_for(w.project.index_of("RUN_A")) is None and np.isfinite(0)
