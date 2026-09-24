"""🏷️ Internal standards (GUI): click-to-label, cross-run propagation, and
the confidence report, driven through MainWindow the way a user would --
select a library row, click a peak, then find it on every other run."""

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
    # RUN_B's peaks sit 1.5s later in t1 than RUN_A's (synthetic_folder fixture) with the
    # same ion and a scaled amplitude -- close enough in position and spectrum for a
    # confident match, not "not found".
    assert match.tier in ("high", "medium")
    assert match.dt1 == pytest.approx(1.5, abs=1.0)

    # The confidence report table reflects it.
    w.istd_panel.refresh()
    rows = w.istd_panel.report_table.rowCount()
    assert rows >= 1
    runs_in_table = {w.istd_panel.report_table.item(r, 0).text() for r in range(rows)}
    assert "RUN_B" in runs_in_table


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
    assert w.istd.matches["RUN_B"][name].tier in ("high", "medium")


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