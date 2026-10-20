"""Drawing the alkane ladder by hand, and the hydrocarbon tools working with either peak finder."""

import os
from types import SimpleNamespace

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from PySide6 import QtCore, QtWidgets  # noqa: E402

from pygcgc import HydrocarbonSettings, fold  # noqa: E402
from pygcgc.core.components import Component, ComponentSet  # noqa: E402
from pygcgc.core.hydrocarbons import (MIN_STEP_S, AlkaneHit, AlkaneSeries, control_points_of,  # noqa: E402
                                      ladder_from_control_points, sample_ladder, snap_ladder)
from pygcgc.core.peaktable import component_positions  # noqa: E402
from tests.test_gui import TIMEOUT, window  # noqa: E402,F401
from tests.test_hydrocarbons import CARBONS, PERIOD, STEP, _build_run, _seeded_series  # noqa: E402


# ---------------------------------------------------------------- the geometry
def test_a_ladder_is_the_smooth_curve_through_its_control_points_with_a_hit_on_every_carbon():
    points = [(10, 100.0, 1.0), (15, 300.0, 1.4), (20, 450.0, 1.2), (30, 700.0, 1.8)]
    series = ladder_from_control_points(points)
    assert series.carbons().tolist() == list(range(10, 31))
    for carbon, t1, t2 in points:
        assert series.hits[carbon].t1 == pytest.approx(t1) and series.hits[carbon].t2 == pytest.approx(t2)
    t1s = [series.hits[c].t1 for c in range(10, 31)]
    assert all(b > a for a, b in zip(t1s, t1s[1:]))                    # runs forwards in time
    assert all(h.seeded and h.height == 0 and h.qual_ratios == () for h in series.hits.values())
    t1, t2 = sample_ladder(points, [12.5, 25.0])
    spline = AlkaneSeries(hits={c: AlkaneHit(c, a, b, 0.0, ()) for c, a, b in points}, shape="smooth")
    assert t1[0] == pytest.approx(spline.t1_of(12.5)) and t2[1] == pytest.approx(spline.t2_of(25.0))


def test_an_s_shaped_ladder_follows_its_bends_and_a_straight_one_does_not():
    s_bend = [(10, 100.0, 1.0), (14, 260.0, 2.0), (18, 420.0, 1.0), (22, 580.0, 2.0)]
    t1, t2 = sample_ladder(s_bend, [12.0, 16.0, 20.0])
    assert t2[0] > 1.0 and t2[1] < 2.0 and t2[2] > 1.0                  # it wanders between the control values
    assert t2[1] < t2[0] or t2[1] < 2.0
    line = ladder_from_control_points([(10, 100.0, 1.0), (30, 900.0, 2.0)])
    assert line.hits[20].t2 == pytest.approx(1.5) and line.hits[20].t1 == pytest.approx(500.0)


def test_bad_control_points_are_refused_with_a_reason():
    with pytest.raises(ValueError, match="at least two"):
        ladder_from_control_points([(10, 1.0, 1.0)])
    with pytest.raises(ValueError, match="different whole carbon"):
        ladder_from_control_points([(10, 1.0, 1.0), (10, 50.0, 1.0)])
    with pytest.raises(ValueError, match="different whole carbon"):
        ladder_from_control_points([(10, 1.0, 1.0), (10.5, 50.0, 1.0)])
    with pytest.raises(ValueError, match="in carbon order"):
        ladder_from_control_points([(10, 100.0, 1.0), (12, 100.0 + MIN_STEP_S / 2, 1.0)])
    with pytest.raises(ValueError, match="in carbon order"):
        ladder_from_control_points([(10, 100.0, 1.0), (12, 50.0, 1.0)])


def test_control_points_of_a_series_keep_the_ends_and_spread_evenly():
    series = ladder_from_control_points([(10, 100.0, 1.0), (30, 900.0, 2.0)])
    assert [c for c, _t1, _t2 in control_points_of(series)] == list(range(10, 31))
    few = control_points_of(series, 5)
    assert len(few) == 5 and few[0][0] == 10 and few[-1][0] == 30
    assert control_points_of(series, 50) == control_points_of(series)


def test_snapping_moves_a_rough_ladder_onto_the_real_alkanes_and_leaves_the_rest():
    run = _build_run()
    f = fold(run, PERIOD, 0.7)
    found = _seeded_series(run, f)                              # where the alkanes really are, as the program finds them
    truth = {c: (found.hits[c].t1, found.hits[c].t2) for c in CARBONS}
    drawn = ladder_from_control_points([(10, truth[10][0] + 8.0, truth[10][1] + 0.08),
                                        (30, truth[30][0] - 7.0, truth[30][1] - 0.06),
                                        (34, truth[30][0] + 4 * STEP, truth[30][1])])   # C31 to C34 are not in the run
    snapped, moved = snap_ladder(run, f, drawn, HydrocarbonSettings(carbon_range=(8, 40)))
    assert set(moved) >= set(range(11, 30)) and 31 not in moved and 34 not in moved
    for c in moved:
        assert abs(snapped.hits[c].t1 - truth[c][0]) < PERIOD
        assert snapped.hits[c].height > 0 and not snapped.hits[c].seeded
    assert snapped.hits[33].t1 == pytest.approx(drawn.hits[33].t1)      # no peak there: kept as drawn
    assert snapped.hits[33].seeded
    with pytest.raises(ValueError):
        snap_ladder(run, f, AlkaneSeries(), HydrocarbonSettings())


def test_fitted_components_give_positions_and_volumes_without_any_spectra():
    def component(i, t1, t2, volume):
        return Component(id=f"c{i}", group_id="g", t1=t1, t2=t2, volume=volume, mz=np.array([57.0]),
                         spectrum=np.array([1.0]), model="gaussian", row0=0, col0=0, sigma_row=1, sigma_col=1)
    parts = [component(0, 100.0, 1.0, 5.0), component(1, 300.0, 2.0, 7.0)]
    peaks = component_positions(ComponentSet(groups=[SimpleNamespace(components=parts)]))
    assert len(peaks) == 2 and peaks.volume.tolist() == [5.0, 7.0] and peaks.component_id == ["c0", "c1"]
    assert peaks.nearest(290.0, 2.1) == 1 and peaks.nearest(110.0, 0.9) == 0
    assert component_positions(ComponentSet(groups=[])).nearest(0, 0) is None


# ---------------------------------------------------------------- dragging
def drag(overlay, which, t1_s, t2):
    """Put a dot where the mouse would have taken it and let go."""
    item = overlay._mover if which < 0 else overlay._handles[which]
    overlay._busy = True
    item.setPos((t1_s / 60.0, t2))
    overlay._busy = False
    overlay._dragged(which)
    overlay._finished()


def open_editor(w, index=None):
    w.show_ladder_editor(w.project.index_of("RUN_A") if index is None else index)
    return w.ladder_editor


def test_the_editor_opens_with_a_rough_ladder_when_there_is_none_and_saves_nothing_until_a_dot_moves(qtbot, window):
    w = window
    w.act_aligned.setEnabled(True)
    w.act_aligned.setChecked(True)
    ed = open_editor(w)
    i = w.project.index_of("RUN_A")
    assert ed.active and ed.panel.isVisible() and not w.act_aligned.isChecked()
    assert len(ed.overlay.points) == 2 and w.hydrocarbons.series_for(i) is None
    assert "no ladder yet" in ed.panel.status.text() and ed.panel.title.text().endswith("RUN_A")
    assert ed.overlay.curve.isVisible() and all(h.isVisible() for h in ed.overlay._handles)
    n = w.undo.count()
    drag(ed.overlay, 0, ed.overlay.points[0][1] + 20.0, ed.overlay.points[0][2] + 0.1)
    series = w.hydrocarbons.series_for(i)
    assert series is not None and series.carbons().tolist() == list(range(10, 31))
    assert w.undo.count() == n + 1 and "saved" in ed.panel.status.text()


def test_an_end_dot_stretches_the_ladder_a_middle_dot_bends_it_and_ctrl_moves_only_the_end(qtbot, window, monkeypatch):
    w = window
    ed = open_editor(w)
    ov = ed.overlay
    ov.load([(10, 100.0, 1.0), (15, 300.0, 1.2), (20, 500.0, 1.4), (25, 700.0, 1.6), (30, 900.0, 1.8)])
    before = [list(p) for p in ov.points]
    drag(ov, 4, 1100.0, 2.8)                                  # pull the last end out by 200 s and 1 s
    assert ov.points[0][1:] == before[0][1:]                  # the far end stays
    assert ov.points[4][1] == pytest.approx(1100.0) and ov.points[4][2] == pytest.approx(2.8)
    assert ov.points[2][1] == pytest.approx(before[2][1] + 0.5 * 200.0)        # the middle follows by its share
    assert ov.points[2][2] == pytest.approx(before[2][2] + 0.5 * 1.0)
    assert [p[0] for p in ov.points] == [10, 15, 20, 25, 30]
    ov.load(before)
    drag(ov, 2, 520.0, 2.4)                                   # a middle dot bends: nothing else moves
    assert ov.points[2][1:] == [520.0, 2.4] and ov.points[1][1:] == before[1][1:] and ov.points[3][1:] == before[3][1:]
    ov.load(before)
    monkeypatch.setattr(QtWidgets.QApplication, "keyboardModifiers", staticmethod(lambda: QtCore.Qt.KeyboardModifier.ControlModifier))
    drag(ov, 4, 1100.0, 2.8)                                  # Ctrl: only that end moves
    assert ov.points[4][1:] == [1100.0, 2.8] and all(ov.points[k][1:] == before[k][1:] for k in range(4))


def test_dragging_cannot_make_the_ladder_run_backwards(qtbot, window):
    w = window
    ov = open_editor(w).overlay
    ov.load([(10, 100.0, 1.0), (15, 300.0, 1.2), (20, 500.0, 1.4)])
    drag(ov, 1, 900.0, 1.2)                                   # a middle dot dragged past its neighbour
    assert ov.points[1][1] == pytest.approx(500.0 - MIN_STEP_S)
    drag(ov, 0, 5000.0, 1.0)                                  # an end dragged right over the others
    assert all(b[1] - a[1] >= MIN_STEP_S - 1e-9 for a, b in zip(ov.points, ov.points[1:]))
    ov.series()                                               # and the series can still be built


def test_the_diamond_moves_the_whole_ladder(qtbot, window):
    w = window
    ov = open_editor(w).overlay
    ov.load([(10, 100.0, 1.0), (20, 500.0, 1.4), (30, 900.0, 1.8)])
    before = [list(p) for p in ov.points]
    mover = ov._mover.pos()
    drag(ov, -1, mover.x() * 60.0 + 30.0, mover.y() + 0.2)
    for p, b in zip(ov.points, before):
        assert p[1] == pytest.approx(b[1] + 30.0) and p[2] == pytest.approx(b[2] + 0.2)


def test_double_click_adds_a_dot_and_the_menu_actions_remove_and_renumber_with_checks(qtbot, window):
    w = window
    ed = open_editor(w)
    ov = ed.overlay
    ov.load([(10, 30.0, 1.0), (20, 120.0, 1.4), (30, 210.0, 1.8)])        # inside the synthetic run's 4 minutes
    messages = []
    ov.inconsistent.connect(messages.append)
    t1, t2 = sample_ladder(ov.points, [15])
    scene = ov.plot.vb.mapViewToScene(QtCore.QPointF(float(t1[0]) / 60.0, float(t2[0])))

    class Click:
        accepted = False

        def double(self):
            return True

        def button(self):
            return QtCore.Qt.MouseButton.LeftButton

        def scenePos(self):
            return scene + QtCore.QPointF(3, 2)

        def accept(self):
            self.accepted = True
    click = Click()
    ov._scene_clicked(click)
    assert click.accepted and [p[0] for p in ov.points] == [10, 15, 20, 30]
    ov._scene_clicked(Click())                                # the same place again: C15 already has its dot
    assert len(ov.points) == 4
    assert not ov.renumber(1, 21) and "out of carbon order" in messages[-1]
    assert ov.renumber(1, 16) and [p[0] for p in ov.points] == [10, 16, 20, 30]
    assert ov.remove(1) and len(ov.points) == 3
    ov.remove(0)
    assert not ov.remove(0) and "at least two" in messages[-1]
    assert not ov.set_ends(20, 15) and "lower than the last" in messages[-1]
    assert ov.set_ends(8, 34) and (ov.points[0][0], ov.points[-1][0]) == (8, 34)


def test_the_end_carbons_in_the_panel_relabel_the_ends(qtbot, window):
    w = window
    ed = open_editor(w)
    drag(ed.overlay, 0, ed.overlay.points[0][1] + 5.0, ed.overlay.points[0][2])
    ed.panel.first.setValue(12)
    ed.panel.last.setValue(26)
    ed._ends_edited()
    assert [p[0] for p in ed.overlay.points] == [12, 26]
    assert w.hydrocarbons.series_for(w.project.index_of("RUN_A")).carbons().tolist() == list(range(12, 27))
    ed.panel.first.setValue(40)                               # not below the last: refused and put back
    ed._ends_edited()
    assert ed.panel.first.value() == 12


def test_a_dot_on_every_carbon_and_back_keeps_the_line_where_it_is(qtbot, window):
    w = window
    ed = open_editor(w)
    ov = ed.overlay
    ov.load([(10, 100.0, 1.0), (15, 300.0, 1.5), (20, 500.0, 1.2), (30, 900.0, 1.8)])
    t1, t2 = sample_ladder(ov.points, np.arange(10, 31))
    ov.set_fine(True)
    assert len(ov.points) == 21
    ov.set_fine(False)
    assert 2 <= len(ov.points) <= 5 and ov.points[0][0] == 10 and ov.points[-1][0] == 30
    again = sample_ladder(ov.points, np.arange(10, 31))
    assert np.allclose(again[0][[0, -1]], t1[[0, -1]]) and np.allclose(again[1][[0, -1]], t2[[0, -1]])


def test_every_edit_is_one_undo_step_and_undo_and_redo_restore_the_ladder(qtbot, window):
    w = window
    ed = open_editor(w)
    i = w.project.index_of("RUN_A")
    w.undo.clear()
    drag(ed.overlay, 0, ed.overlay.points[0][1] + 20.0, ed.overlay.points[0][2])
    first = w.hydrocarbons.series_for(i).copy()
    drag(ed.overlay, 1, ed.overlay.points[1][1] - 30.0, ed.overlay.points[1][2] + 0.3)
    second = w.hydrocarbons.series_for(i).copy()
    assert w.undo.count() == 2 and second.hits[30].t1 != first.hits[30].t1
    w.undo.undo()
    assert w.hydrocarbons.series_for(i).hits[30].t1 == pytest.approx(first.hits[30].t1)
    assert ed.overlay.points[-1][1] == pytest.approx(first.hits[30].t1)       # the overlay follows (a dot per carbon)
    w.undo.undo()
    assert w.hydrocarbons.series_for(i) is None                               # back to no ladder
    w.undo.redo()
    w.undo.redo()
    assert w.hydrocarbons.series_for(i).hits[30].t1 == pytest.approx(second.hits[30].t1)


def test_the_drawn_ladder_is_what_the_matcher_and_the_other_steps_read(qtbot, window, tmp_path):
    from tests.test_library_matcher import matched_window
    w = window
    win, state = matched_window(qtbot, w, tmp_path)
    i = w.index
    assert state.current_for(i) is not None and "no alkane ladder" in state.ladder_for(i)[3]
    ed = open_editor(w, i)
    drag(ed.overlay, 0, ed.overlay.points[0][1] + 10.0, ed.overlay.points[0][2])
    series, borrowed, _map, note = state.ladder_for(i)
    assert series is w.hydrocarbons.series_for(i) and not borrowed and "own alkane ladder, C10 to C30" in note
    assert state.current_for(i) is None and state.results_for(i) is not None          # matches are out of date, not lost
    assert "changed since this run was matched" in win.requirements.text()
    assert w.hydrocarbons.can_assign(i) == ""                                         # assignment needs only the ladder and peaks


def test_stored_ladders_survive_a_save_with_their_drawn_flags(qtbot, window):
    w = window
    ed = open_editor(w)
    drag(ed.overlay, 0, ed.overlay.points[0][1] + 10.0, ed.overlay.points[0][2])
    saved = w.hydrocarbons.to_dict()
    w.hydrocarbons.from_dict(saved)
    series = w.hydrocarbons.series_for(w.project.index_of("RUN_A"))
    assert series.carbons().tolist() == list(range(10, 31)) and all(h.seeded for h in series.hits.values())


def test_switching_runs_edits_the_new_runs_ladder(qtbot, window):
    w = window
    ed = open_editor(w)
    drag(ed.overlay, 0, ed.overlay.points[0][1] + 10.0, ed.overlay.points[0][2])
    w.set_run(w.project.index_of("RUN_B"))
    assert ed.index == w.index and ed.panel.title.text().endswith("RUN_B") and "no ladder yet" in ed.panel.status.text()
    ed.close()
    assert not ed.active and not ed.overlay.curve.isVisible() and not ed.panel.isVisible()


def test_the_hydrocarbon_panel_button_opens_the_editor_for_the_open_run(qtbot, window):
    w = window
    w.set_run(w.project.index_of("RUN_A"))
    w.hydro_panel.draw_ladder_btn.click()
    assert w.ladder_editor is not None and w.ladder_editor.active and w.ladder_editor.index == w.index


# ---------------------------------------------------------------- snapping in the app
def test_snap_to_peaks_uses_a_job_and_is_one_undo_step(qtbot, window, monkeypatch):
    # The synthetic run has no m/z 57, so the search itself is checked above; here the wiring is.
    from pygcgc.gui import hydrocarbon_state
    w = window
    ed = open_editor(w)
    i = w.project.index_of("RUN_A")
    drag(ed.overlay, 0, ed.overlay.points[0][1] + 5.0, ed.overlay.points[0][2])
    drawn = w.hydrocarbons.series_for(i).copy()

    def fake_snap(run, fold, series, settings, cancel=None):
        out = series.copy()
        for c in (12, 13):
            out.hits[c] = AlkaneHit(c, series.hits[c].t1 + 2.0, series.hits[c].t2, 9.0, (0.66, 0.44), seeded=False)
        return out, [12, 13]
    monkeypatch.setattr(hydrocarbon_state, "snap_ladder", fake_snap)
    n = w.undo.count()
    ed.panel.snap_btn.click()
    assert not ed.panel.snap_btn.isEnabled()
    qtbot.waitUntil(lambda: not w.hydrocarbons.running and w.undo.count() > n, timeout=TIMEOUT)
    assert w.undo.count() == n + 1 and "Snap 2 of 21" in w.undo.text(w.undo.index() - 1)
    assert ed.panel.fine.isChecked() and "2 of 21 carbons were moved onto a peak" in ed.panel.status.text()
    assert ed.panel.snap_btn.isEnabled()
    assert w.hydrocarbons.series_for(i).hits[12].t1 == pytest.approx(drawn.hits[12].t1 + 2.0)
    assert len(ed.overlay.found.data["x"]) == 2                        # the measured apexes are marked
    w.undo.undo()
    assert w.hydrocarbons.series_for(i).hits[12].t1 == pytest.approx(drawn.hits[12].t1)


def test_a_snap_that_cannot_run_says_so_and_leaves_the_ladder_alone(qtbot, window):
    w = window
    ed = open_editor(w)
    i = w.project.index_of("RUN_A")
    drag(ed.overlay, 0, ed.overlay.points[0][1] + 5.0, ed.overlay.points[0][2])
    drawn = w.hydrocarbons.series_for(i).copy()
    n = w.undo.count()
    ed.panel.snap_btn.click()                                           # the synthetic run has no m/z 57
    qtbot.waitUntil(lambda: ed.panel.snap_btn.isEnabled(), timeout=TIMEOUT)
    assert "Could not snap the ladder" in ed.panel.status.text() and w.undo.count() == n
    assert w.hydrocarbons.series_for(i).hits[20].t1 == pytest.approx(drawn.hits[20].t1)


# ---------------------------------------------------------------- MS-aware peaks
def test_assignment_markers_and_stencil_totals_work_with_ms_aware_peaks(qtbot, window):
    from pygcgc.gui.blob_state import MS_AWARE
    w = window
    w.blobs.set_method(MS_AWARE)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.components_for(i) is not None, timeout=TIMEOUT)
    components = w.blobs.components_for(i).components
    assert components and w.blobs.blobs_for(i) is None                        # no legacy blobs exist
    hydro = w.hydrocarbons
    positions = hydro.positions_for(i)
    assert len(positions) == len(components)
    # a ladder drawn through the synthetic peaks' region, then everything that reads it
    ed = open_editor(w, i)
    drag(ed.overlay, 0, ed.overlay.points[0][1] + 5.0, ed.overlay.points[0][2])
    assert hydro.can_assign(i) == ""
    hydro.add_stencil(i, [(0.0, 0.0), (3000.0, 0.0), (3000.0, 10.0), (0.0, 10.0)], "everything")
    totals = hydro.totals(i)
    assert totals["everything"]["n_blobs"] == len(components)
    assert totals["everything"]["volume"] == pytest.approx(sum(c.volume for c in components))
    assert set(hydro.markers(i)) == {"pristane", "phytane"}
    assert hydro.calls_for(i) is None and hydro.spectra_for(i) is None        # spectra are not built until asked
    done = []
    w.hydro_panel.assignRequested.emit()
    qtbot.waitUntil(lambda: hydro.calls_for(i) is not None, timeout=TIMEOUT)
    assert len(hydro.calls_for(i)) == len(components)
    assert w.hydro_panel.assign_table.rowCount() == len(components)
    done.append(1)
    # moving the ladder redoes the assignments at once, now that the spectra are ready
    before = [c.carbon_eq for c in hydro.calls_for(i)]
    drag(ed.overlay, 0, ed.overlay.points[0][1] + 60.0, ed.overlay.points[0][2])
    after = hydro.calls_for(i)
    assert after is not None and len(after) == len(components) and [c.carbon_eq for c in after] != before
    w.goto_hydrocarbon_peak(0)                                                          # and a row can be visited
    assert done
