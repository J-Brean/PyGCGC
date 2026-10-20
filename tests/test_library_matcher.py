"""The library matcher window and its state."""

import csv
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from pygcgc.core.libmatch import FLAGS, MatchSettings  # noqa: E402
from tests.nist_fixture import write_mini_nist  # noqa: E402
from tests.test_gui import TIMEOUT, _write_test_msp, window  # noqa: E402,F401


def matched_window(qtbot, w, tmp_path):
    """The synthetic project with peaks found, an MSP library loaded and RUN_A matched."""
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    w.set_run(w.project.index_of("RUN_A"))
    w.show_library_matcher()
    state = w.project.matcher
    state.load_library(_write_test_msp(tmp_path / "lib.msp"))
    qtbot.waitUntil(lambda: state.library is not None, timeout=TIMEOUT)
    w.library_matcher.match_btn.click()
    qtbot.waitUntil(lambda: state.current_for(w.index) is not None, timeout=TIMEOUT)
    return w.library_matcher, state


def test_missing_library_and_peaks_are_named_and_matching_is_disabled(qtbot, window):
    w = window
    w.set_run(w.project.index_of("RUN_A"))
    w.show_library_matcher()
    win = w.library_matcher
    text = win.requirements.text()
    assert "Load a library" in text and "Find peaks first" in text
    assert not win.match_btn.isEnabled()
    assert "No library loaded" in win.library_label.text()


def test_no_alkane_ladder_is_stated_and_does_not_block_matching(qtbot, window, tmp_path):
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    assert "no alkane ladder" in win.requirements.text().lower()
    assert win.match_btn.isEnabled()
    match = state.current_for(w.index).matches[0]
    assert "RI_MISSING" in match.flags and match.level != "confident"


def test_the_known_peak_is_named_and_the_tables_are_largest_first(qtbot, window, tmp_path):
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    found = state.current_for(w.index)
    target = min(found.matches, key=lambda m: abs(m.query.t1 - 40.0))
    assert target.best.name == "Compound50"
    shares = [float(win.table.item(r, 3).text()) for r in range(win.table.rowCount())]
    assert shares == sorted(shares, reverse=True) and win.table.rowCount() == len(found.matches)
    assert not hasattr(win, "summary")                        # one table, largest first: no separate "largest peaks"
    volumes = sorted((m.query.volume for m in found.matches), reverse=True)
    top_t1 = next(m.query.t1 for m in found.matches if m.query.volume == volumes[0])
    assert win.table.item(0, 1).text() == f"{top_t1 / 60:.3f}"


def test_selecting_a_peak_shows_the_spectra_candidates_and_every_uncertainty(qtbot, window, tmp_path):
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    win.table.selectRow(0)
    match = win._shown
    assert match is not None and win.candidates.rowCount() == len(match.candidates) >= 1
    assert "Spectrum: match factor" in win.detail_text.text()
    for flag in match.flags:
        assert flag in win.detail_text.text() and FLAGS[flag].split(".")[0][:20] in win.detail_text.text()
    assert len(win.up_bars.opts["x"]) > 0 and len(win.down_bars.opts["x"]) > 0
    assert win.up_bars.opts["height"].max() == pytest.approx(100.0)
    assert win.down_bars.opts["height"].min() == pytest.approx(-100.0)


def test_the_table_says_where_each_peak_stands_against_a_subtracted_blank(qtbot, window, tmp_path):
    from pygcgc.core.peak_blank import PeakBlankSettings
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    assert all(win.table.item(r, 14).text() == "" for r in range(win.table.rowCount()))      # no subtraction yet
    w.project.blank.set_settings(PeakBlankSettings(blank_run="RUN_B"))
    w.project.blank.subtract_runs([w.index])
    qtbot.waitUntil(lambda: w.project.blank.current_for(w.index) is not None, timeout=TIMEOUT)
    win._fill_results()
    words = [win.table.item(r, 14).text() for r in range(win.table.rowCount())]
    assert words and all(x.startswith(("● above blank", "○ near blank")) for x in words)


def test_filters_narrow_the_table_without_losing_the_selection_logic(qtbot, window, tmp_path):
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    total = win.table.rowCount()
    win.name_filter.setText("compound50")
    assert 0 < win.table.rowCount() <= total
    assert all("compound50" in win.table.item(r, 4).text().lower() for r in range(win.table.rowCount()))
    win.name_filter.setText("")
    win.level_filter.setCurrentIndex(win.level_filter.findData("good"))
    assert win.table.rowCount() == sum(m.level in ("confident", "probable") for m in state.current_for(w.index).matches)
    win.level_filter.setCurrentIndex(0)
    win.min_share.setValue(99.0)
    assert win.table.rowCount() <= 1
    win.min_share.setValue(0.0)
    assert win.table.rowCount() == total


def test_a_decision_is_recorded_beside_the_level_and_never_changes_it(qtbot, window, tmp_path):
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    win.table.selectRow(0)
    label, level = win._shown.query.label, win._shown.level
    win.accept_btn.click()
    d = state.decisions_for(w.index)[label]
    assert d["action"] == "accepted" and d["name"] == win._shown.best.name
    assert win._shown.level == level
    assert "accepted" in win.table.item(0, 13).text()
    win.reject_btn.click()
    assert state.decisions_for(w.index)[label]["action"] == "rejected"
    win.clear_btn.click()
    assert label not in state.decisions_for(w.index)


def test_decisions_survive_saving_only_while_the_inputs_are_the_same(qtbot, window, tmp_path):
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    win.table.selectRow(0)
    label = win._shown.query.label
    win.accept_btn.click()
    saved = state.to_dict()
    assert saved["decisions"]["RUN_A"][label]["action"] == "accepted"
    state.from_dict(saved)
    assert state.results == {} and state.library is None            # results are recomputed, not saved
    state.load_library(tmp_path / "lib.msp")
    qtbot.waitUntil(lambda: state.library is not None, timeout=TIMEOUT)
    state.match_runs([w.index])
    qtbot.waitUntil(lambda: state.current_for(w.index) is not None, timeout=TIMEOUT)
    assert label in state.decisions_for(w.index)                    # same inputs: the decision stands
    from dataclasses import replace
    state.set_settings(replace(state.settings, ri_gate_z=2.5))
    state.match_runs([w.index])
    qtbot.waitUntil(lambda: state.current_for(w.index) is not None, timeout=TIMEOUT)
    assert label not in state.decisions_for(w.index)               # different settings: it does not silently carry over


def test_changing_the_settings_makes_the_results_out_of_date_not_wrong(qtbot, window, tmp_path):
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    from dataclasses import replace
    state.set_settings(replace(state.settings, mf_confident=900.0))
    assert state.current_for(w.index) is None and state.results_for(w.index) is not None
    assert "changed since this run was matched" in win.requirements.text()
    assert "out of date" in win.table_note.text()


def test_results_are_discarded_if_the_peaks_change_while_matching(qtbot, window, tmp_path, monkeypatch):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    w.set_run(w.project.index_of("RUN_A"))
    state = w.project.matcher
    state.load_library(_write_test_msp(tmp_path / "lib.msp"))
    qtbot.waitUntil(lambda: state.library is not None, timeout=TIMEOUT)
    submitted = []
    monkeypatch.setattr(w.jobs, "submit", lambda *a, **k: submitted.append((a, k)))
    assert state.match_runs([w.index]) and state.running
    args, kwargs = submitted[-1]
    result = args[1](__import__("pygcgc.gui.jobs", fromlist=["CancelToken"]).CancelToken(), lambda *_: None)
    w.project.period = w.project.period + 0.5                      # the fold moved under the job
    args[2](result)
    assert not state.running and state.results_for(w.index) is None


def test_export_writes_every_doubt_and_the_decision(qtbot, window, tmp_path):
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    win.table.selectRow(0)
    win.accept_btn.click()
    out = tmp_path / "ids.csv"
    win.export_to(out)
    rows = list(csv.DictReader(out.open(encoding="utf-8")))
    assert len(rows) == len(state.current_for(w.index).matches)
    assert float(rows[0]["percent_of_signal"]) >= float(rows[-1]["percent_of_signal"])
    assert "RI_MISSING" in rows[0]["flags"] and FLAGS["RI_MISSING"] in rows[0]["flag_text"]
    assert rows[0]["decision"] == "accepted" and rows[0]["column_class"] == "Semi-standard non-polar"


def test_main_export_includes_identifications_and_the_library_settings(qtbot, window, tmp_path):
    import json
    from pygcgc.gui.export import export_project
    w = window
    matched_window(qtbot, w, tmp_path)
    out = tmp_path / "export"
    export_project(w.project, out)
    index = w.index + 1
    assert (out / f"identifications_{index:03d}.csv").exists()
    manifest = json.loads((out / "manifest.json").read_text())
    row = manifest["runs"][w.index]
    assert row["identifications"] == f"identifications_{index:03d}.csv" and row["identifications_stale"] is False
    provenance = json.loads((out / "provenance.json").read_text())
    assert provenance["settings"]["library_match"]["column_class"] == "Semi-standard non-polar"


def test_promoting_a_match_appends_a_target_where_the_peak_was_found(qtbot, window, tmp_path):
    from pygcgc.core.targets import load_targets
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    match = min(state.current_for(w.index).matches, key=lambda m: abs(m.query.t1 - 40.0))
    path = tmp_path / "targets.csv"
    win.promote_to_file(match, match.best, path)
    targets = load_targets(path)
    assert targets[0].compound == "Compound50" and targets[0].t1 == pytest.approx(match.query.t1)


def test_the_nist_folder_loads_with_literature_indices_for_the_chosen_column(qtbot, window, tmp_path):
    w = window
    write_mini_nist(tmp_path / "nist")
    w.set_run(w.project.index_of("RUN_A"))
    w.show_library_matcher()
    state = w.project.matcher
    state.load_library(tmp_path / "nist")
    qtbot.waitUntil(lambda: state.library is not None, timeout=TIMEOUT)
    assert "NIST" in state.library_status and "literature retention index" in state.library_status
    assert len(state.library) == 7 and (state.library.ri_exp_n > 0).sum() >= 5
    # a polar column has far fewer literature values in the fixture
    from dataclasses import replace
    state.set_settings(replace(state.settings, column_class="Standard polar"))
    qtbot.waitUntil(lambda: not state.loading and "standard polar" in state.library_status, timeout=TIMEOUT)
    assert (state.library.ri_exp_n > 0).sum() < 5


def test_a_borrowed_ladder_is_used_when_a_run_has_none(qtbot, window):
    from pygcgc.core.hydrocarbons import AlkaneHit, AlkaneSeries
    w = window
    series = AlkaneSeries()
    for c, t in ((10, 100.0), (11, 200.0), (12, 300.0)):
        series.hits[c] = AlkaneHit(carbon=c, t1=t, t2=2.0, height=1.0, qual_ratios=())
    w.project.hydrocarbons.series["RUN_B"] = series
    state = w.project.matcher
    i_a = w.project.index_of("RUN_A")
    assert state.ladder_for(i_a)[0] is None
    state.set_ladder_run("RUN_B")
    got, borrowed, t1_map, note = state.ladder_for(i_a)
    assert got is series and borrowed and "RUN_B" in note and "no alignment" in note
    assert t1_map(123.0, 2.0) == pytest.approx(123.0)
    assert state.ladder_for(w.project.index_of("RUN_B"))[1] is False        # its own ladder is its own
    assert MatchSettings().column_class == "Semi-standard non-polar"
    assert np.isfinite(t1_map(1.0, 1.0))
