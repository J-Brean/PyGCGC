"""🏷️ Internal standards (GUI): click to label, the three thresholds, the
search of every chromatogram and the results table, driven through
MainWindow the way a user would."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
pytest.importorskip("pyqtgraph")

from tests.test_gui import TIMEOUT, _peak_pixel, window  # noqa: E402,F401


def _label_first_library_row(w, run_index, ion):
    w.set_run(run_index)
    w.istd_panel.library_table.selectRow(0)
    name = w.istd_panel.selected_name()
    assert name
    col, t2 = _peak_pixel(w, run_index, ion)
    w.start_label_istd()
    assert w._istd_label_pending == name
    w._point_clicked(col, t2)
    assert w._istd_label_pending is None
    return name


def test_label_click_captures_the_reference_peak(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    i_a = w.project.index_of("RUN_A")
    name = _label_first_library_row(w, i_a, 50)
    qtbot.waitUntil(lambda: name in w.istd.references, timeout=TIMEOUT)
    ref = w.istd.references[name]
    assert ref.run == "RUN_A"
    assert ref.t1 == pytest.approx(40.0, abs=2.0)   # GUI_PEAKS[0]: t1=40s, ion 50


def test_propagate_finds_the_shifted_peak_on_the_other_run(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    i_a = w.project.index_of("RUN_A")
    name = _label_first_library_row(w, i_a, 50)
    qtbot.waitUntil(lambda: name in w.istd.references, timeout=TIMEOUT)

    w.propagate_istd()
    qtbot.waitUntil(lambda: "RUN_B" in w.istd.matches and name in w.istd.matches["RUN_B"], timeout=TIMEOUT)
    match = w.istd.matches["RUN_B"][name]
    # RUN_B's peaks sit 1.5 s later in t1 than RUN_A's (synthetic_folder fixture) with the
    # same ion and a scaled amplitude: inside the default thresholds on all three.
    assert match.found
    assert match.dt1 == pytest.approx(1.5, abs=1.0)

    # Every chromatogram is searched, the labelled one included, and the table shows it:
    # the standard down the side, the runs across the top, a tick where it passed.
    table = w.istd_results
    assert set(w.istd.matches) == {e.name for e in w.project.entries if e.run is not None}
    assert table.standards() == [name]
    col = table.runs().index("RUN_B")
    assert table.table.item(0, col).text().startswith("✓")
    assert "Spectral match" in table.table.item(0, col).toolTip()


def test_thresholds_set_in_the_panel_decide_what_is_found(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    name = _label_first_library_row(w, w.project.index_of("RUN_A"), 50)
    qtbot.waitUntil(lambda: name in w.istd.references, timeout=TIMEOUT)
    # RUN_B's standard is 1.5 s later in t1; with only one standard there is no drift
    # to follow, so a 0.5 s first-dimension limit must reject it and say why.
    w.istd_panel.max_dt1.setValue(0.5)
    assert w.istd.settings.max_dt1 == pytest.approx(0.5)
    w.go("istd")
    qtbot.waitUntil(lambda: "RUN_B" in w.istd.matches and "RUN_A" in w.istd.matches, timeout=TIMEOUT)
    assert w.istd.matches["RUN_A"][name].found                 # the labelled peak itself, at Δ 0
    miss = w.istd.matches["RUN_B"][name]
    assert not miss.found and "no peak within 0.5 s in t1" in miss.reason
    col = w.istd_results.runs().index("RUN_B")
    cell = w.istd_results.table.item(0, col)
    assert cell.text().startswith("✗") and "no peak within" in cell.toolTip()
    # Clicking the cell opens that run.
    w.istd_results.cellChosen.emit("RUN_B", name)
    assert w.project.entries[w.index].name == "RUN_B"


def test_cancel_pick_clears_pending_label(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    i_a = w.project.index_of("RUN_A")
    w.set_run(i_a)
    w.istd_panel.library_table.selectRow(0)
    w.start_label_istd()
    assert w._istd_label_pending is not None
    w.cancel_pick()
    assert w._istd_label_pending is None


def test_labelling_without_peak_detection_refuses(qtbot, window):
    w = window
    i_a = w.project.index_of("RUN_A")
    w.set_run(i_a)
    w.istd_panel.library_table.selectRow(0)
    w.start_label_istd()
    assert w._istd_label_pending is None   # refused: no peaks detected yet on this run


def test_ms_aware_labelling_round_trips_through_respectrum(qtbot, window):
    """The default MS-aware detector has no per-component spectrum until
    `respectrum_components` resolves it -- labelling and matching must both
    go through that path, not just the legacy BlobSet one."""
    from pygcgc.gui.blob_state import MS_AWARE

    w = window
    w.blobs_panel.method_box.setCurrentIndex(w.blobs_panel.method_box.findData(MS_AWARE))
    assert w.blobs.method == MS_AWARE
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.candidates_were_run, timeout=TIMEOUT)
    i_a = w.project.index_of("RUN_A")
    name = _label_first_library_row(w, i_a, 50)
    qtbot.waitUntil(lambda: name in w.istd.references, timeout=TIMEOUT)

    w.propagate_istd()
    qtbot.waitUntil(lambda: "RUN_B" in w.istd.matches and name in w.istd.matches["RUN_B"], timeout=TIMEOUT)
    assert w.istd.matches["RUN_B"][name].found


def test_go_on_the_internal_standard_step_searches_every_run(qtbot, window):
    # Go on step 11 used to raise KeyError: the step had no entry in MainWindow.go.
    w = window
    w.workflow.show_step("istd")
    w.go("istd")                                   # nothing labelled yet: says so, does not crash
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    name = _label_first_library_row(w, w.project.index_of("RUN_A"), 50)
    qtbot.waitUntil(lambda: name in w.istd.references, timeout=TIMEOUT)
    w.go("istd")
    qtbot.waitUntil(lambda: "RUN_B" in w.istd.matches, timeout=TIMEOUT)
    assert "istd" in [a for a in w.step_actions]


def test_labelled_standards_become_the_target_table(qtbot, window, tmp_path):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    name = _label_first_library_row(w, w.project.index_of("RUN_A"), 50)
    qtbot.waitUntil(lambda: name in w.istd.references, timeout=TIMEOUT)
    before = {t.compound for t in w.quant.targets}          # the synthetic folder's own table
    w.istd_to_targets(path=tmp_path / "targets.csv")
    rows = {t.compound: t for t in w.quant.targets}
    assert rows[name].role == "istd" and rows[name].quant_ion == 50
    assert before - {name} <= set(rows)                      # the other rows are kept
    assert name in [a.name for a in w.project.anchors]