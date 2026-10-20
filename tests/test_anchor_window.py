"""The anchoring window: phase check, picking, ANCHOR, and the before-and-after check."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from PySide6 import QtCore  # noqa: E402

from pygcgc.core.align import (AUTO_MIN, FIT_MIN, AnchorPoint, anchor_offsets, anchor_requirements,  # noqa: E402
                               image_similarity)
from tests.test_gui import TIMEOUT, _peak_pixel, _wait_render, confirm_phases, window  # noqa: E402,F401


# ------------------------------------------------------------------ pure functions
def test_requirements_follow_the_library_thresholds_not_a_copy_of_them():
    assert anchor_requirements(0)["auto"] is None
    assert anchor_requirements(1)["auto"] == "shift"
    assert anchor_requirements(AUTO_MIN["affine"] - 1)["auto"] == "shift"
    assert anchor_requirements(AUTO_MIN["affine"])["auto"] == "affine"
    assert anchor_requirements(AUTO_MIN["smooth"])["auto"] == "smooth"
    fits = {f["kind"]: f for f in anchor_requirements(3)["fits"]}
    assert fits["affine"]["least"] == FIT_MIN["affine"] and fits["affine"]["auto"] == AUTO_MIN["affine"]


def test_offsets_are_the_median_distance_between_reference_and_run():
    points = [AnchorPoint("P1", 100.0, 1.0, 97.0, 1.10), AnchorPoint("P2", 200.0, 1.5, 198.0, 1.30),
              AnchorPoint("P3", 300.0, 2.0, 296.0, 2.20)]
    dt1, dt2 = anchor_offsets(points)
    assert dt1 == pytest.approx(3.0) and dt2 == pytest.approx(0.2)
    assert np.isnan(anchor_offsets([])[0])


def test_image_similarity_rewards_matching_peaks_and_ignores_empty_baseline():
    rng = np.random.default_rng(0)
    ref = rng.random((60, 80)) * 0.01
    for r, c in ((10, 10), (30, 40), (50, 70), (20, 60)):
        ref[r - 1:r + 2, c - 1:c + 2] += 100.0
    same = ref + rng.random(ref.shape) * 0.01
    shifted = np.roll(ref, (5, 4), axis=(0, 1))
    assert image_similarity(ref, same) > 0.95
    assert image_similarity(ref, shifted) < image_similarity(ref, same) - 0.5
    assert np.isnan(image_similarity(np.zeros((20, 20)), np.zeros((20, 20))))
    with pytest.raises(ValueError):
        image_similarity(ref, ref[:10])


# ------------------------------------------------------------------ state
def test_clicks_are_only_recorded_while_holding_and_anchor_fits_them(window):
    w = window
    w.set_reference("RUN_A")
    messages = []
    w.project.message.connect(messages.append)
    w.align.begin_hold()
    w.align.add_anchor("RUN_A", "P1", 100.0, 1.0, 100.0, 1.0, 50)
    w.align.add_anchor("RUN_B", "P1", 100.0, 1.0, 98.0, 1.1, 50)
    assert not w.align.transforms and w.align.pending_runs() == ["RUN_B"]
    assert not any("Alignment changed" in m for m in messages)
    w.align.apply_hold()
    assert "RUN_B" in w.align.transforms and not w.align.pending_runs() and not w.align.holding


def test_other_edits_end_the_hold_and_apply(window):
    w = window
    w.set_reference("RUN_A")
    w.align.begin_hold()
    w.align.add_anchor("RUN_B", "P1", 100.0, 1.0, 98.0, 1.1, 50)
    w.align.restore(w.align.snapshot())                  # what undo, refinement and deleting rows all do
    assert not w.align.holding and "RUN_B" in w.align.transforms
    w.align.begin_hold()
    w.set_reference("RUN_B")
    assert not w.align.holding and not w.align.points


def test_a_phase_change_after_anchoring_marks_the_run_and_clears_its_tick(window):
    w = window
    w.set_reference("RUN_A")
    confirm_phases(w)
    w.align.add_anchor("RUN_B", "P1", 100.0, 1.0, 98.0, 1.1, 50)
    assert w.align.stale_phase_runs() == []
    index = w.project.index_of("RUN_B")
    w.project.set_phase(index, w.project.entries[index].phase + 0.5)
    assert w.align.stale_phase_runs() == ["RUN_B"]
    assert "RUN_B" not in w.align.phase_checked and "RUN_A" in w.align.phase_checked
    w.project.set_phases({index: (w.project.entries[index].phase, "manual", "")})
    assert w.align.phase_checked == set()                  # every run's tick goes when all phases are set at once


def test_phase_record_and_ticks_survive_saving_and_old_projects_load_unchanged(window):
    w = window
    w.set_reference("RUN_A")
    confirm_phases(w)
    w.align.add_anchor("RUN_B", "P1", 100.0, 1.0, 98.0, 1.1, 50)
    saved = w.align.to_dict()
    assert saved["phase_at_anchor"].keys() == {"RUN_B"} and set(saved["phase_checked"]) == {e.name for e in w.project.entries}
    w.align.from_dict(saved)
    assert w.align.phase_at_anchor == saved["phase_at_anchor"] and w.align.phase_checked == set(saved["phase_checked"])
    old = {k: v for k, v in saved.items() if k not in ("phase_at_anchor", "phase_checked")}
    w.align.from_dict(old)
    assert w.align.stale_phase_runs() == [] and w.align.phase_checked == set()
    assert "phase_at_anchor" not in w.align.to_dict() or w.align.to_dict()["phase_at_anchor"] == {}


# ------------------------------------------------------------------ the window
def test_picking_cannot_start_until_every_phase_is_ticked(qtbot, window):
    w = window
    w.set_ion_expression("50")
    w.set_reference("RUN_A")
    w.start_add_anchor()
    win = w.anchor_window
    assert win is not None and win.isVisible() and w._anchor_stage is None
    assert win.pages.currentIndex() == 0 and not win.continue_btn.isEnabled()
    assert "phase" in win.phase_heading.text().lower()
    assert win.phase_table.rowCount() == 3
    win._tick_all()
    assert win.continue_btn.isEnabled()
    qtbot.mouseClick(win.continue_btn, QtCore.Qt.MouseButton.LeftButton)
    assert w._anchor_stage == "ref" and win.pages.currentIndex() == 1 and w.align.holding


def test_unticking_one_run_blocks_continuing(window):
    w = window
    w.set_reference("RUN_A")
    w.start_add_anchor()
    win = w.anchor_window
    win._tick_all()
    win.align.set_phase_checked("RUN_B", False)
    assert not win.continue_btn.isEnabled() and "Tick every run" in win.phase_warning.text()


def test_header_names_the_run_and_its_role_and_the_requirement_box_counts(qtbot, window):
    w = window
    w.set_ion_expression("50")
    confirm_phases(w)
    w.start_add_anchor()
    win = w.anchor_window
    assert win.run_heading.text() == "Run 1 of 3: RUN_A" and win.role_badge.text() == "REFERENCE"
    assert "reference" in win.instruction.text().lower() and "P1" in win.instruction.text()
    assert win.skip_btn.isHidden() and not win.finish_btn.isHidden()
    for k, name in enumerate(("RUN_A", "RUN_B")):
        index = w.project.index_of(name)
        col, t2 = _peak_pixel(w, index, 50)
        w.anchor_pane_click(name, col, t2)
    assert win.run_heading.text() == "Run 3 of 3: RUN_C_BLANK" and win.role_badge.text() == "SAMPLE"
    assert not win.skip_btn.isHidden()
    assert [c.toolTip() for c in win.chips[:3]] == ["has this anchor", "has this anchor", "being picked now"]
    w.skip_anchor_run()                                  # the blank has no such peak
    assert w._anchor_pending["name"] == "P2"
    # RUN_C_BLANK was skipped, so no anchor is complete on every sample yet
    assert "complete on every sample: 0" in win.requirements.text()
    assert "a shift: 1 or more" in win.requirements.text()


def test_a_click_on_the_windows_own_image_is_converted_to_the_runs_pixel(qtbot, window):
    w = window
    w.set_ion_expression("50")
    confirm_phases(w)
    w.start_add_anchor()
    win = w.anchor_window
    qtbot.waitUntil(lambda: win._shown_run == "RUN_A", timeout=TIMEOUT)
    index = w.project.index_of("RUN_A")
    col, t2 = _peak_pixel(w, index, 50)
    fold = w.project.fold(index)

    class Click:
        def button(self):
            return QtCore.Qt.MouseButton.LeftButton

        def scenePos(self):
            point = QtCore.QPointF((fold.t1[col] + fold.spec.period / 2) / 60.0, t2)
            return win.pane.plot.vb.mapViewToScene(point)

    win._pane_clicked(Click())
    assert w.align.points["RUN_A"][0].ion == 50
    assert w.align.points["RUN_A"][0].ref_t2 == pytest.approx(t2, abs=fold.row_pitch)
    assert w._anchor_stage == "run"
    # a click while the window still shows the previous run is refused, not recorded on the wrong run
    win._shown_run = "RUN_A"
    before = w.align.snapshot()
    win._pane_clicked(Click())
    assert w.align.snapshot() == before


def test_the_reference_position_is_drawn_on_sample_runs(qtbot, window):
    w = window
    w.set_ion_expression("50")
    confirm_phases(w)
    w.start_add_anchor()
    win = w.anchor_window
    col, t2 = _peak_pixel(w, w.project.index_of("RUN_A"), 50)
    w.anchor_pane_click("RUN_A", col, t2)
    qtbot.waitUntil(lambda: win._shown_run == "RUN_B", timeout=TIMEOUT)
    assert len(win.ghost.data) == 1 and win.ghost.data["y"][0] == pytest.approx(t2)
    xs, ys = win.box.getData()
    assert len(xs) == 5 and ys.max() - ys.min() == pytest.approx(0.3)


def test_anchor_aligns_and_the_result_page_reports_before_and_after(qtbot, window):
    w = window
    w.set_ion_expression("50")
    confirm_phases(w)
    w.start_add_anchor()
    win = w.anchor_window
    for name in ("RUN_A", "RUN_B"):
        col, t2 = _peak_pixel(w, w.project.index_of(name), 50)
        w.anchor_pane_click(name, col, t2)
    w.skip_anchor_run()
    assert "RUN_B" not in w.align.transforms and "not been applied" in win.pending_note.text()
    assert win.anchor_btn.isEnabled()
    qtbot.mouseClick(win.anchor_btn, QtCore.Qt.MouseButton.LeftButton)
    assert w._anchor_stage is None and "RUN_B" in w.align.transforms and win.pages.currentIndex() == 2
    qtbot.waitUntil(lambda: "RUN_B" in win._report, timeout=TIMEOUT)
    assert win.result_table.rowCount() == 1 and win.result_table.item(0, 0).text() == "RUN_B"
    assert win.result_table.item(0, 1).text() == "shift" and win.result_table.item(0, 2).text() == "1"
    before, after = win._report["RUN_B"]["before"], win._report["RUN_B"]["after"]
    assert np.isfinite(before) and np.isfinite(after)
    # One anchor can only shift; the status must say honestly whether that helped.
    expected = "Worse after alignment" if after + 0.02 < before else "OK"
    assert win.result_table.item(0, 7).text().startswith(expected)
    # both comparison views are on the reference's grid, with and without the transform
    assert win.before_view.index_a == w.project.index_of("RUN_A") and win.after_view.index_b == w.project.index_of("RUN_B")
    assert win.before_view.image_source(win.before_view.index_b)[0].shape == win.after_view.image_source(win.after_view.index_b)[0].shape
    win._back_to_anchors()
    assert win.pages.currentIndex() == 1


def test_anchor_without_a_matching_click_says_so_and_changes_nothing(window):
    w = window
    w.set_reference("RUN_A")
    w.align.add_anchor("RUN_A", "P1", 100.0, 1.0, 100.0, 1.0, 50)
    before = w.align.snapshot()
    w.apply_anchors()
    assert w.align.snapshot() == before and "RUN_B" not in w.align.transforms


def test_closing_the_window_ends_picking_but_keeps_the_clicks(qtbot, window):
    w = window
    w.set_ion_expression("50")
    confirm_phases(w)
    w.start_add_anchor()
    col, t2 = _peak_pixel(w, w.project.index_of("RUN_A"), 50)
    w.anchor_pane_click("RUN_A", col, t2)
    w.anchor_window.close()
    assert w._anchor_stage is None and w.align.points["RUN_A"][0].name == "P1"


def test_leaving_a_changed_phase_unticked_shows_which_run_was_anchored_differently(qtbot, window):
    w = window
    w.set_ion_expression("50")
    confirm_phases(w)
    w.start_add_anchor()
    for name in ("RUN_A", "RUN_B"):
        col, t2 = _peak_pixel(w, w.project.index_of(name), 50)
        w.anchor_pane_click(name, col, t2)
    w.skip_anchor_run()
    w.apply_anchors()
    index = w.project.index_of("RUN_B")
    w.project.set_phase(index, w.project.entries[index].phase + 0.7)
    win = w.anchor_window
    win.pages.setCurrentIndex(0)
    win.refresh()
    labels = [win.phase_table.item(r, 1).text() for r in range(win.phase_table.rowCount())]
    assert any("RUN_B" in t and "different phase" in t for t in labels)
    assert "Anchored with a different phase" in win.phase_warning.text()


def test_the_align_panel_reminds_about_phase_and_opens_the_window(qtbot, window):
    w = window
    w.workflow.show_step("align")
    assert "phase" in w.align_panel.phase_note.text().lower()
    qtbot.mouseClick(w.align_panel.window_btn, QtCore.Qt.MouseButton.LeftButton)
    assert w.anchor_window is not None and w.anchor_window.isVisible()
