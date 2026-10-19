"""The anchor route must collect corresponding clicks across the whole sequence."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from pygcgc.gui.jobs import CancelToken
from tests.test_gui import TIMEOUT, _peak_pixel, _wait_render, window  # noqa: F401


def click_peak(qtbot, w, ion=50):
    name = w.project.entries[w.index].name
    _wait_render(qtbot, w, name)
    qtbot.waitUntil(lambda: w.view2d.result.expr == w.renderer.expr and not w.view2d.result.aligned,
                    timeout=TIMEOUT)
    col, t2 = _peak_pixel(w, w.index, ion)
    w._point_clicked(col, t2)


def test_mouse_clicks_and_visible_controls_complete_a_round(qtbot, window):
    from PySide6 import QtCore
    w = window
    w.workflow.show_step("align")
    qtbot.mouseClick(w.workflow.go_btn, QtCore.Qt.MouseButton.LeftButton)
    for name in ("RUN_A", "RUN_B"):
        _wait_render(qtbot, w, name)
        col, t2 = _peak_pixel(w, w.index, 50)
        f = w.project.fold(w.index)
        point = QtCore.QPointF((f.t1[col]+f.spec.period/2)/60, t2)
        scene_point = w.view2d.plot.vb.mapViewToScene(point)
        position = w.view2d.canvas.mapFromScene(scene_point)
        qtbot.mouseClick(w.view2d.canvas.viewport(), QtCore.Qt.MouseButton.LeftButton, pos=position)
        assert w.align.points[name][0].name == "P1"
    assert w.project.entries[w.index].name == "RUN_C_BLANK"
    qtbot.mouseClick(w.align_panel.skip_btn, QtCore.Qt.MouseButton.LeftButton)
    assert w.project.entries[w.index].name == "RUN_A" and w._anchor_pending["name"] == "P2"
    qtbot.mouseClick(w.workflow.go_btn, QtCore.Qt.MouseButton.LeftButton)
    assert w._anchor_stage is None and w.align_panel.refine_btn.isEnabled()
    running = []
    w.align.refinementRunning.connect(running.append)
    qtbot.mouseClick(w.align_panel.refine_btn, QtCore.Qt.MouseButton.LeftButton)
    qtbot.waitUntil(lambda: not w.align.refining, timeout=TIMEOUT)
    assert running == [True, False]
    assert w.align_panel.table.isColumnHidden(2) and not w.align_panel.table.isColumnHidden(4)
    qtbot.mouseClick(w.align_panel.more_btn, QtCore.Qt.MouseButton.LeftButton)
    assert not w.align_panel.table.isColumnHidden(2)


def test_two_anchor_rounds_start_on_first_run_and_keep_a_shared_table(qtbot, window):
    w = window
    w.set_ion_expression("50")
    assert w.align.reference is None
    w.go("align")
    assert w.align.reference == "RUN_A" and w._anchor_stage == "ref"
    assert w.align_panel.add_btn.isVisibleTo(w.workflow)
    assert "Finish picking" in w.workflow.go_btn.text()
    for anchor in ("P1", "P2"):
        assert w._anchor_pending["name"] == anchor
        for run in ("RUN_A", "RUN_B", "RUN_C_BLANK"):
            assert w.project.entries[w.index].name == run
            click_peak(qtbot, w)
            assert any(p.name == anchor for p in w.align.points[run])
        assert w.project.entries[w.index].name == "RUN_A"
        assert w._anchor_stage == "ref"
    assert w.align_panel.table.rowCount() == 6
    assert w._anchor_pending["name"] == "P3"
    w.go("align")  # the same primary action finishes the guided loop
    assert w._anchor_stage is None and not w.align_panel.skip_btn.isVisible()
    saved = w.align.to_dict()
    w.align.from_dict(saved)
    assert w.align.to_dict() == saved
    _wait_render(qtbot, w, "RUN_A")
    w._update_anchor_overlay()
    assert len(w.view2d.anchor_markers.points()) == 2
    w.workflow.show_step("background")
    assert len(w.view2d.anchor_markers.points()) == 0


def test_custom_reference_visits_other_runs_in_sequence_and_can_skip(qtbot, window):
    w = window
    w.set_reference("RUN_B")
    w.start_add_anchor()
    assert w.project.entries[w.index].name == "RUN_B"
    click_peak(qtbot, w)
    assert w.project.entries[w.index].name == "RUN_A"
    w.align_panel.skip_btn.click()
    assert w.project.entries[w.index].name == "RUN_C_BLANK"
    w.align_panel.skip_btn.click()
    assert w.project.entries[w.index].name == "RUN_B"
    assert w._anchor_pending["name"] == "P2"
    assert "RUN_A" not in w.align.points
    w.cancel_pick()
    assert w.align.points["RUN_B"][0].name == "P1"


def test_click_on_previous_image_is_ignored_and_manual_navigation_finishes(qtbot, window):
    w = window
    w.start_add_anchor()
    click_peak(qtbot, w)
    assert w.view2d.result.name == "RUN_A" and w.project.entries[w.index].name == "RUN_B"
    before = w.align.snapshot()
    w._point_clicked(0, .1)
    assert w.align.snapshot() == before  # RUN_B has not drawn yet
    w.set_run(w.project.index_of("RUN_C_BLANK"))
    assert w._anchor_stage is None
    assert w.align.snapshot() == before


def test_changing_reference_or_leaving_alignment_finishes_collection(qtbot, window):
    w = window
    w.start_add_anchor()
    w.set_reference("RUN_B")
    assert w._anchor_stage is None
    assert w.align_panel.reference.currentData() == "RUN_B"
    w.undo.undo()
    assert w.align_panel.reference.currentData() == "RUN_A"
    w.undo.redo()
    assert w.align_panel.reference.currentData() == "RUN_B"
    w.start_add_anchor()
    w.workflow.show_step("background")
    assert w._anchor_stage is None


def test_undoing_reference_click_finishes_instead_of_recording_an_orphan_pair(qtbot, window):
    w = window
    w.start_add_anchor()
    click_peak(qtbot, w)
    assert w._anchor_stage == "run"
    w.undo.undo()
    assert w._anchor_stage is None and not w.align.points
    w.undo.redo()
    assert w.align.points["RUN_A"][0].name == "P1"


def test_matching_click_uses_the_reference_ion_and_input_changes_finish(qtbot, window):
    w = window
    w.set_ion_expression("50")
    w.start_add_anchor()
    click_peak(qtbot, w)
    w.set_ion_expression("51")
    _wait_render(qtbot, w, "RUN_B")
    qtbot.waitUntil(lambda: w.view2d.result.expr == w.renderer.expr, timeout=TIMEOUT)
    w._point_clicked(10, 1.5)
    assert "RUN_B" not in w.align.points
    w.project.foldChanged.emit(-1)
    assert w._anchor_stage is None


def test_table_removes_the_selected_runs_point_not_the_displayed_run(qtbot, window):
    w = window
    w.start_add_anchor()
    click_peak(qtbot, w)
    click_peak(qtbot, w)
    w.cancel_pick()  # currently the blank, but select RUN_B's row
    table = w.align_panel.table
    row = next(r for r in range(table.rowCount()) if table.item(r, 1).text() == "RUN_B")
    table.selectRow(row)
    w.align_panel._remove_selected()
    assert "RUN_B" not in w.align.points and "RUN_A" in w.align.points
    w.undo.undo()
    assert "RUN_B" in w.align.points
    w.remove_anchor("P1", "RUN_A")
    assert not w.align.points  # removing a reference removes its corresponding pairs
    w.undo.undo()
    assert "RUN_B" in w.align.points


def seed_offset_clicks(w, *, reference_row=True):
    w.set_reference("RUN_A")
    positions = {}
    for name in ("RUN_A", "RUN_B"):
        index = w.project.index_of(name)
        col, t2 = _peak_pixel(w, index, 50)
        f = w.project.fold(index)
        positions[name] = (float(f.t1[col] + f.spec.period/2), t2)
    ref = positions["RUN_A"]
    if reference_row:
        w.align.add_anchor("RUN_A", "P1", ref[0], ref[1]+.1, ref[0], ref[1]+.1, 50)
    target = positions["RUN_B"]
    w.align.add_anchor("RUN_B", "P1", ref[0], ref[1]+.1, target[0], target[1]-.1, 50)
    return positions


def test_visible_add_delete_and_restart_controls_are_undoable(qtbot, window):
    from PySide6 import QtCore
    w = window
    w.workflow.show_step("align")
    panel = w.align_panel
    assert panel.add_btn.isVisibleTo(w.workflow) and panel.delete_btn.isVisibleTo(w.workflow)
    assert panel.restart_btn.isVisibleTo(w.workflow) and not panel.delete_btn.isEnabled()
    qtbot.mouseClick(panel.add_btn, QtCore.Qt.MouseButton.LeftButton)
    click_peak(qtbot, w)
    click_peak(qtbot, w)
    w.cancel_pick()
    before = w.align.snapshot()
    panel.table.selectAll()
    count = w.undo.count()
    qtbot.mouseClick(panel.delete_btn, QtCore.Qt.MouseButton.LeftButton)
    assert not w.align.points and w.undo.count() == count+1
    w.undo.undo()
    assert w.align.snapshot() == before
    w.align.set_kind("RUN_B", "shift")
    before = w.align.snapshot()
    w.workflow.mark_done("align", "Old completion")
    count = w.undo.index()
    qtbot.mouseClick(panel.restart_btn, QtCore.Qt.MouseButton.LeftButton)
    assert not w.align.points and not w.align.kind and not w.align.transforms
    assert w.align.reference == "RUN_A" and w._anchor_pending["name"] == "P1"
    assert w._anchor_stage == "ref" and w.index == w.project.index_of("RUN_A")
    assert not w.workflow.is_done("align") and not w.workflow.banner.isVisible()
    assert w.undo.index() == count+1
    w.undo.undo()
    assert w.align.snapshot() == before and w._anchor_stage is None
    w.undo.redo()
    assert not w.align.points and not w.align.kind


@pytest.mark.parametrize("reference_row", [True, False])
def test_refinement_updates_reference_and_pairs_in_one_undo_without_editing_data(qtbot, window, reference_row):
    w = window
    positions = seed_offset_clicks(w, reference_row=reference_row)
    before = w.align.snapshot()
    raw = [e.run.data.copy() for e in w.project.entries]
    history = w.undo.count()
    w.refine_anchors()
    assert w.align.refining and not w.align_panel.refine_btn.isEnabled()
    qtbot.waitUntil(lambda: not w.align.refining, timeout=TIMEOUT)
    for run, points in w.align.points.items():
        for p in points:
            assert (p.ref_t1, p.ref_t2) == pytest.approx(positions["RUN_A"])
            assert (p.run_t1, p.run_t2) == pytest.approx(positions[run])
    after = w.align.snapshot()
    assert after != before and w.undo.count() == history+1
    for entry, original in zip(w.project.entries, raw):
        assert (entry.run.data != original).nnz == 0
    w.undo.undo()
    assert w.align.snapshot() == before
    w.undo.redo()
    assert w.align.snapshot() == after


@pytest.mark.parametrize("change", ["anchors", "fold", "layers", "cancel", "failure"])
def test_refinement_discards_stale_cancelled_or_failed_jobs(window, monkeypatch, change):
    w = window
    seed_offset_clicks(w)
    submitted = []
    monkeypatch.setattr(w.jobs, "submit", lambda *args, **kwargs: submitted.append((args, kwargs)))
    w.refine_anchors()
    args, kwargs = submitted[-1]
    assert w.align.refining
    if change == "cancel":
        kwargs["on_cancel"]()
    elif change == "failure":
        args[3]("unreadable signal")
    else:
        result = args[1](CancelToken(), lambda *_: None)
        if change == "anchors":
            w.align.remove_anchor("RUN_B", "P1")
        elif change == "fold":
            w.project.foldChanged.emit(-1)
        else:
            w.project.layersChanged.emit()
        before = w.align.snapshot()
        args[2](result)
        assert w.align.snapshot() == before
    assert not w.align.refining
