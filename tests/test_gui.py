"""GUI tests (pytest-qt). They run headless with QT_QPA_PLATFORM=offscreen."""

import os
import time

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
pytest.importorskip("pyqtgraph")

from pygcgc.core import circstats as cs  # noqa: E402
from pygcgc.gui.jobs import JobManager  # noqa: E402
from tests.conftest import GUI_PEAKS, GUI_PHASE  # noqa: E402

TIMEOUT = 20000


# ------------------------------------------------------------------- jobs
def test_job_runs_and_reports(qtbot):
    jm = JobManager()
    out = []
    jm.submit("add", lambda tok, prog: 2 + 2, on_done=out.append)
    qtbot.waitUntil(lambda: out == [4], timeout=TIMEOUT)


def test_job_failure_reported(qtbot):
    jm = JobManager()
    errors = []

    def boom(tok, prog):
        raise RuntimeError("nope")

    jm.submit("boom", boom, on_error=errors.append)
    qtbot.waitUntil(lambda: bool(errors), timeout=TIMEOUT)
    assert "nope" in errors[0]


def test_job_cancel(qtbot):
    jm = JobManager()
    results = []

    def slow(tok, prog):
        for i in range(200):
            prog(i, 200, "working")
            time.sleep(0.005)
        return "finished"

    jid = jm.submit("slow", slow, on_done=results.append)
    qtbot.wait(50)
    jm.cancel(jid)
    qtbot.waitUntil(lambda: jm.jobs.get(jid) is None, timeout=TIMEOUT)
    assert results == []


def test_latest_request_wins(qtbot):
    jm = JobManager(max_threads=4)
    got = []

    def make(v, delay):
        def fn(tok, prog):
            time.sleep(delay)
            return v
        return fn

    jm.submit("old", make("old", 0.3), on_done=got.append, latest_key="k", quiet=True)
    jm.submit("new", make("new", 0.0), on_done=got.append, latest_key="k", quiet=True)
    qtbot.waitUntil(lambda: jm.pool.activeThreadCount() == 0, timeout=TIMEOUT)
    qtbot.wait(100)
    assert got == ["new"]


# ------------------------------------------------------------------ window
@pytest.fixture
def window(qtbot, synthetic_folder):
    from pygcgc.gui.mainwindow import MainWindow
    from pygcgc.gui.blob_state import LEGACY
    w = MainWindow()
    w.blobs.set_method(LEGACY)  # older GUI flows below intentionally exercise legacy blobs
    w.jobs_panel.clear_after_ms = 10
    qtbot.addWidget(w)
    w.resize(1400, 900)
    w.show()
    loaded = []
    w.open_folder(synthetic_folder, on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and w.jobs.active_count() == 0 and w.view2d.result is not None,
                    timeout=TIMEOUT)
    qtbot.waitUntil(lambda: all(e.phase_method != "unset" for e in w.project.entries), timeout=TIMEOUT)
    yield w
    w.jobs.cancel_all()
    w.jobs.wait(5000)
    w.undo.clear()   # no "unsaved changes" prompt on close


def _wait_render(qtbot, w, name):
    qtbot.waitUntil(lambda: w.view2d.result is not None and w.view2d.result.name == name,
                    timeout=TIMEOUT)


def test_folder_opens_with_sequence_and_phases(window):
    p = window.project
    assert [e.name for e in p.entries] == ["RUN_A", "RUN_B", "RUN_C_BLANK"]
    assert [e.kind for e in p.entries] == ["calibration", "sample", "tube blank"]
    assert p.period == 4.0
    a, b, blank = p.entries
    assert a.phase_method == "anchors" and b.phase_method == "anchors"
    assert abs(cs.circ_diff(a.phase, GUI_PHASE, 4.0)) < 0.01
    assert blank.phase_method == "fallback"
    r = window.view2d.result
    assert r.image.shape[0] == 200
    assert window.runs_panel.table.rowCount() == 3


def test_flick_keeps_view_and_levels(qtbot, window):
    window.view2d.plot.setRange(xRange=(1.0, 2.5), yRange=(0.5, 3.0), padding=0)
    qtbot.wait(50)
    view_before = window.view2d.plot.vb.viewRange()
    levels_before = tuple(window.view2d.image.getLevels())
    window.act_next.trigger()
    _wait_render(qtbot, window, "RUN_B")
    assert window.index == 1
    np.testing.assert_allclose(window.view2d.plot.vb.viewRange(), view_before)
    assert tuple(window.view2d.image.getLevels()) == pytest.approx(levels_before)
    window.act_prev.trigger()
    _wait_render(qtbot, window, "RUN_A")


def test_ion_expression(qtbot, window):
    assert window.set_ion_expression("51")
    qtbot.waitUntil(lambda: window.view2d.result.expr.text() == "51", timeout=TIMEOUT)
    img = window.view2d.result.image
    r, c = np.unravel_index(np.argmax(img), img.shape)
    assert abs(r * window.view2d.result.pitch - 1.90) < 0.03
    assert not window.set_ion_expression("banana")
    assert not window.set_ion_expression("999")          # outside the mass range
    assert window.renderer.expr.text() == "51"
    assert window.set_ion_expression("50/51")
    assert window.view2d.scale == "linear"               # ratios aren't shown on a log scale


def test_spectrum_click_and_ion_pick(qtbot, window):
    f = window.project.fold(0)
    t1, t2, ion, _ = GUI_PEAKS[1]
    col, _ = f.to_fold(t1)
    window.show_spectrum_at(float(f.to_time(col, t2)))
    spec = window.spectrum
    assert int(spec.mz[np.argmax(spec.intensity)]) == ion
    window.ion_clicked(52, False)
    assert window.renderer.expr.text() == "52"
    window.ion_clicked(53, True)
    assert window.renderer.expr.text() == "52+53"


def test_region_spectrum(qtbot, window):
    t1, t2, ion, _ = GUI_PEAKS[2]
    window.show_region_spectrum(t1 - 10, t1 + 10, t2 - 0.3, t2 + 0.3)
    spec = window.spectrum
    assert int(spec.mz[np.argmax(spec.intensity)]) == ion
    window.act_region.setChecked(True)
    assert window.view2d.region() is not None


def test_phase_edit_undo_and_merge(qtbot, window):
    e = window.project.entries[0]
    original = e.phase
    window.fold_panel.phase.setValue(1.20)
    window.fold_panel.phase.setValue(1.30)       # merges with the previous edit
    assert e.phase == pytest.approx(1.30)
    assert e.phase_method == "manual"
    assert window.undo.count() == 1
    window.undo.undo()
    assert e.phase == pytest.approx(original)
    assert e.phase_method == "anchors"
    qtbot.waitUntil(lambda: abs(window.view2d.result.phase - original) < 1e-9, timeout=TIMEOUT)


def test_copy_then_fix_phase_undo(qtbot, window):
    window.fold_panel.phase.setValue(2.0)
    window.copy_phase_to_all()
    assert all(x.phase == pytest.approx(2.0) for x in window.project.entries)
    n = window.undo.count()
    window.fix_phase()
    qtbot.waitUntil(lambda: window.undo.count() == n + 1, timeout=TIMEOUT)
    assert window.project.entries[0].phase_method == "anchors"
    window.undo.undo()
    assert all(x.phase == pytest.approx(2.0) for x in window.project.entries)
    window.undo.redo()
    assert window.project.entries[0].phase_method == "anchors"


def test_period_edit_undo(qtbot, window):
    window.fold_panel.period.setValue(4.01)
    assert window.project.period == pytest.approx(4.01)
    window.undo.undo()
    assert window.project.period == pytest.approx(4.0)
    assert window.project.entries[0].phase_method == "anchors"


def test_band_sync(qtbot, window):
    window.view1d.set_band(60.0, 120.0, emit=True)
    a, b = window.view2d.band.getRegion()
    assert (a, b) == pytest.approx((1.0, 2.0))
    window.view2d.band.setRegion((2.0, 3.0))
    window.view2d._band_moved()
    assert window.view1d.band_range() == pytest.approx((120.0, 180.0))


def test_extended_view(qtbot, window):
    window.act_extend.setChecked(True)
    qtbot.waitUntil(lambda: window.view2d.result.extended, timeout=TIMEOUT)
    r = window.view2d.result
    assert r.image.shape[0] == round(200 * 1.15)
    assert window.view2d.period_line.isVisible()


def test_scale_cycle(window):
    assert window.view2d.scale == "log"
    window.act_scale.trigger()
    assert window.view2d.scale == "sqrt"
    assert window.scale_box.currentData() == "sqrt"


def test_project_roundtrip(qtbot, window, tmp_path):
    from pygcgc.gui.mainwindow import MainWindow
    window.fold_panel.phase.setValue(1.11)
    window.set_ion_expression("51")
    window.act_extend.setChecked(True)
    path = tmp_path / "proj.pygcgc.json"
    window.project.path = path
    window.save_project()
    assert window.undo.isClean()

    w2 = MainWindow()
    qtbot.addWidget(w2)
    loaded = []
    w2.open_project(path, on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and w2.view2d.result is not None, timeout=TIMEOUT)
    e = w2.project.entries[0]
    assert e.phase == pytest.approx(1.11)
    assert e.phase_method == "manual"            # saved phases are not re-estimated
    assert w2.project.entries[1].phase_method == "anchors"
    assert w2.renderer.expr.text() == "51"
    assert w2.renderer.extended
    w2.jobs.wait(5000)


def test_remove_runs(qtbot, window):
    window.fold_panel.phase.setValue(1.0)
    window.remove_runs([2])
    assert len(window.project) == 2
    assert window.undo.count() == 0


# ------------------------------------------------------------ 🧹 background
def _apply_bg(qtbot, w):
    n = w.undo.count()
    w.act_background.trigger()
    qtbot.waitUntil(lambda: w.undo.count() == n + 1, timeout=TIMEOUT)
    qtbot.waitUntil(lambda: w.view2d.result is not None and w.view2d.result.processed, timeout=TIMEOUT)


def test_background_apply_toggle_undo(qtbot, window):
    w = window
    assert not w.act_processed.isEnabled()
    _apply_bg(qtbot, w)
    assert all(e.background is not None for e in w.project.entries)
    assert w.act_processed.isEnabled() and w.act_processed.isChecked()
    assert w.layers_panel.background.isChecked()
    assert 54 in w.project.entries[0].background.ions_processed
    raw_total = w.project.entries[0].run.xic(54).sum()
    cleaned_total = w.project.display_run(0).xic(54).sum()
    assert raw_total > 4 * 2e5
    assert cleaned_total == pytest.approx(2e5, rel=0.15)     # only the real m/z 54 peak is left
    # Space shows raw again.
    w.act_processed.trigger()
    qtbot.waitUntil(lambda: not w.view2d.result.processed, timeout=TIMEOUT)
    assert w.project.display_run(0) is w.project.entries[0].run
    w.act_processed.trigger()
    qtbot.waitUntil(lambda: w.view2d.result.processed, timeout=TIMEOUT)
    # Undo removes the layer everywhere; redo brings it back.
    w.undo.undo()
    assert not w.project.any_background()
    assert not w.act_processed.isEnabled()
    w.undo.redo()
    assert w.project.any_background()


def test_background_spectrum_uses_cleaned_data(qtbot, window):
    w = window
    _apply_bg(qtbot, w)
    f = w.project.fold(0)
    col, _ = f.to_fold(150.0)
    w.show_spectrum_at(float(f.to_time(col, 0.3)))     # a quiet spot: only bleed on m/z 54
    spec = w.spectrum
    assert 54 not in spec.mz[spec.intensity > 30]


def test_background_preview_and_split(qtbot, window):
    w = window
    w.bg_panel.preview.setChecked(True)
    qtbot.waitUntil(lambda: w.project.preview is not None, timeout=TIMEOUT)
    qtbot.waitUntil(lambda: w.view2d.result.preview, timeout=TIMEOUT)
    assert not w.project.any_background()               # preview isn't applied
    assert "Preview" in w.bg_panel.summary.text()
    w.bg_panel.split.setChecked(True)
    qtbot.waitUntil(lambda: w.view2d.result.image_before is not None, timeout=TIMEOUT)
    assert w.view2d.divider.isVisible()
    # Changing a setting recomputes the preview.
    old = w.project.preview[1]
    w.bg_panel.percentile.setValue(30)
    qtbot.waitUntil(lambda: w.project.preview is not None and w.project.preview[1] is not old,
                    timeout=TIMEOUT)
    assert w.project.preview[1].settings.percentile == 30
    w.bg_panel.preview.setChecked(False)
    assert w.project.preview is None


def test_background_diagnostics_and_snr(qtbot, window):
    w = window
    _apply_bg(qtbot, w)
    w.bg_panel.diag_ion.setValue(54)
    text = w.bg_panel.diag_text.text()
    assert "Residual mean" in text and "Centred" in text
    assert w.act_snr.isEnabled()
    w.set_ion_expression("54")
    w.act_snr.setChecked(True)
    qtbot.waitUntil(lambda: w.view2d.result.snr, timeout=TIMEOUT)
    img = w.view2d.result.image[:, 5:-5]
    spread = 1.4826 * np.nanmedian(np.abs(img - np.nanmedian(img)))
    assert spread == pytest.approx(1.0, rel=0.35)      # pure noise has an S/N spread of about 1


def test_background_follows_phase_change(qtbot, window):
    w = window
    _apply_bg(qtbot, w)
    old = w.project.entries[0].background
    w.fold_panel.phase.setValue(1.5)
    qtbot.waitUntil(lambda: w.project.entries[0].background is not old, timeout=TIMEOUT)
    assert w.project.entries[0].background is not None


def test_background_with_blank(qtbot, window):
    w = window
    i = w.bg_panel.blank.findData("RUN_C_BLANK")
    assert i > 0
    assert w.bg_panel.blank.itemText(1).startswith("RUN_C_BLANK")     # blanks listed first
    w.bg_panel.blank.setCurrentIndex(i)
    _apply_bg(qtbot, w)
    a = w.project.entries[0].background
    blank = w.project.entries[2].background
    assert a.blank_name == "RUN_C_BLANK"
    assert blank.blank_name is None                     # the blank isn't subtracted from itself
    assert w.project.bg_blank == "RUN_C_BLANK"


def test_background_saved_and_rebuilt(qtbot, window, tmp_path):
    from pygcgc.gui.mainwindow import MainWindow
    w = window
    w.bg_panel.percentile.setValue(25)
    _apply_bg(qtbot, w)
    path = tmp_path / "bg.pygcgc.json"
    w.project.path = path
    w.save_project()
    w2 = MainWindow()
    qtbot.addWidget(w2)
    done = []
    w2.open_project(path, on_loaded=lambda: done.append(True))
    qtbot.waitUntil(lambda: bool(done), timeout=TIMEOUT)
    assert all(e.background is not None for e in w2.project.entries)
    assert w2.project.bg_settings.percentile == 25
    assert w2.bg_panel.settings().percentile == 25
    qtbot.waitUntil(lambda: w2.view2d.result is not None and w2.view2d.result.processed, timeout=TIMEOUT)
    w2.jobs.wait(5000)


def test_background_remove_all(qtbot, window):
    w = window
    _apply_bg(qtbot, w)
    w.remove_background_all()
    assert not w.project.any_background()
    w.undo.undo()
    assert w.project.any_background()


# ------------------------------------------------------------ 🫧 siloxanes
def test_siloxane_apply_undo_and_layers_panel(qtbot, window):
    # The synthetic_folder runs only cover m/z 50-54, well below any
    # siloxane ion, so this exercises the plumbing (job, undo, layers
    # panel, display_run) rather than the fit itself -- that's covered on
    # synthetic data built for it in tests/test_siloxane.py.
    w = window
    n = w.undo.count()
    w.act_siloxane.trigger()
    qtbot.waitUntil(lambda: w.undo.count() == n + 1, timeout=TIMEOUT)
    assert all(e.siloxane is not None for e in w.project.entries)
    assert all(e.siloxane.ions_used == () for e in w.project.entries)   # nothing in range
    assert w.project.any_siloxane()
    assert w.layers_panel.siloxane.isChecked()
    assert w.act_processed.isEnabled() and w.act_processed.isChecked()
    # The (unchanged) siloxane run is what's now displayed.
    assert w.project.display_run(0) is w.project.entries[0].siloxane.run
    w.undo.undo()
    assert not w.project.any_siloxane()
    w.undo.redo()
    assert w.project.any_siloxane()


def test_siloxane_on_top_of_background_and_remove_all(qtbot, window):
    w = window
    _apply_bg(qtbot, w)
    w.act_siloxane.trigger()
    qtbot.waitUntil(lambda: w.project.any_siloxane(), timeout=TIMEOUT)
    # Built from the cleaned run, not raw.
    e = w.project.entries[0]
    assert e.siloxane.total_before == pytest.approx(e.background.run.data.sum())
    w.remove_siloxane_all()
    assert not w.project.any_siloxane()
    assert w.project.any_background()          # 🧹 itself is untouched
    w.undo.undo()
    assert w.project.any_siloxane()


def test_siloxane_flags_shown_in_blobs_panel(qtbot, window):
    w = window
    w.act_siloxane.trigger()
    qtbot.waitUntil(lambda: w.project.any_siloxane(), timeout=TIMEOUT)
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    col = w.blobs_panel.table.columnCount() - 1
    assert w.blobs_panel.table.horizontalHeaderItem(col).text() == "🫧 siloxane?"
    # None of GUI_PEAKS' ions are siloxane ions, so every score is 0% and nothing is flagged.
    values = [w.blobs_panel.table.item(r, col).text() for r in range(w.blobs_panel.table.rowCount())]
    assert values and all(v == "0%" for v in values)


def test_blob_siloxane_removal_blocks_blobs_keeps_them_findable_and_undoes(qtbot, window):
    w = window
    w.sil_panel.threshold.setValue(0)             # certainty >= 0%: every blob counts as siloxane
    w.sil_panel.block_all.setChecked(True)        # the synthetic runs carry no siloxane ions to block
    assert w.project.siloxane_settings.score_threshold == 0.0
    assert w.project.siloxane_settings.method == "blobs"
    n_before = w.undo.count()
    w.act_siloxane.trigger()                      # finds the blobs first, then blocks them out
    qtbot.waitUntil(lambda: w.project.any_siloxane() and w.undo.count() == n_before + 1, timeout=TIMEOUT)
    e = w.project.entries[0]
    assert e.siloxane.n_blobs == len(w.blobs.by_run[e]) > 0
    assert e.siloxane.run.data.nnz < e.run.data.nnz
    assert w.project.display_run(0) is e.siloxane.run
    # Blobs were found beneath the layer, so they are still there and not out of date.
    assert e not in w.blobs.stale and len(w.blobs.by_run[e]) == e.siloxane.n_blobs
    # Applying again replaces the layer; it does not stack on it.
    first = e.siloxane.run
    w.act_siloxane.trigger()
    qtbot.waitUntil(lambda: w.undo.count() == n_before + 2, timeout=TIMEOUT)
    assert w.project.entries[0].siloxane.run.data.nnz == first.data.nnz
    w.undo.undo()
    w.undo.undo()
    assert not w.project.any_siloxane()
    assert w.project.display_run(0) is e.run


# ------------------------------------------------------------ 🧪 quantification
@pytest.fixture
def qwin(qtbot, cal_folder):
    from pygcgc.gui.mainwindow import MainWindow
    from pygcgc.gui.blob_state import LEGACY
    w = MainWindow()
    w.blobs.set_method(LEGACY)
    w.jobs_panel.clear_after_ms = 10
    qtbot.addWidget(w)
    w.resize(1400, 900)
    w.show()
    loaded = []
    w.open_folder(cal_folder, on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and w.jobs.active_count() == 0, timeout=TIMEOUT)
    qtbot.waitUntil(lambda: all(e.phase_method != "unset" for e in w.project.entries), timeout=TIMEOUT)
    w.run_quant()
    qtbot.waitUntil(lambda: w.quant.report is not None and not w.quant.running, timeout=TIMEOUT)
    yield w
    w.jobs.cancel_all()
    w.jobs.wait(5000)
    w.undo.clear()


def _row(w, run, compound):
    return next(r for r in w.quant.report.results if r.run == run and r.compound == compound)


def test_quant_runs_and_fills_panels(qtbot, qwin):
    w = qwin
    r = w.quant.report
    fits = [f for d in r.fits.values() for k, f in d.items() if k == "quant"]
    assert len(fits) == 3 and all(f.ok for f in fits)
    s1a = _row(w, "S1", "A")
    assert s1a.amount == pytest.approx(3.0, rel=0.05)
    assert s1a.concentration == pytest.approx(1.5, rel=0.05)
    assert w.results.model.rowCount() == len(r.results)
    assert w.calibration_panel.compounds.count() == 3
    w.tabs.setCurrentWidget(w.review)
    w.review.show_compound("A")
    assert len(w.review.tiles) == len(r.seq)
    # Every run reports its peaks against the nearest internal standard.
    s1 = [x for x in r.signals if x.run == "S1"]
    assert s1 and all(x.istd == "IS" for x in s1 if x.istd)
    assert {x.compound for x in r.istd} == {"IS"}


def test_quant_overlay_and_goto(qtbot, qwin):
    w = qwin
    i = w.project.index_of("S1")
    w.set_run(i)
    qtbot.waitUntil(lambda: w.view2d.result is not None and w.view2d.result.name == "S1", timeout=TIMEOUT)
    labels = {t.toPlainText() for t in w.view2d.overlay_labels}
    assert {"IS (IS)", "A", "B", "C"} <= labels
    w.goto_peak("S2", "B")
    qtbot.waitUntil(lambda: w.view2d.result is not None and w.view2d.result.name == "S2"
                    and w.view2d.result.expr.text() == "54", timeout=TIMEOUT)
    qtbot.wait(150)
    spec = w.spectrum
    assert int(spec.mz[np.argmax(spec.intensity)]) == 54
    (x0, x1), _ = w.view2d.plot.vb.viewRange()
    assert x0 < (170.0 + 1.5) / 60 < x1


def test_quant_decisions_undo(qtbot, qwin):
    w = qwin
    w.decide([("S1", "A")], "rejected")
    assert _row(w, "S1", "A").status == "rejected"
    w.undo.undo()
    assert _row(w, "S1", "A").status != "rejected"
    w.undo.redo()
    assert _row(w, "S1", "A").status == "rejected"
    # Review keys
    w.tabs.setCurrentWidget(w.review)
    w.review.show_compound("B")
    w.review.set_current(0)
    first = w.review.tiles[0]
    from PySide6 import QtCore
    qtbot.keyClick(w.review, QtCore.Qt.Key.Key_Y)
    assert _row(w, first.run, "B").status == "accepted"


def test_quant_exclude_and_override(qtbot, qwin):
    w = qwin
    fit = w.quant.report.fit("A")
    n = fit.n
    w.exclude_point("A", "CAL_L5", "pipetting error")
    assert w.quant.report.fit("A").n == n - 1
    w.override_curve("A", {"weights": "1"})
    assert w.quant.report.fit("A").weights == "1"
    w.undo.undo()
    w.undo.undo()
    assert w.quant.report.fit("A").n == n
    assert "A" not in w.quant.decisions.overrides


def test_quant_pick_by_hand(qtbot, qwin):
    w = qwin
    w.start_pick("S2", "C")
    qtbot.waitUntil(lambda: w.view2d.result is not None and w.view2d.result.name == "S2", timeout=TIMEOUT)
    hit = w.quant.run_quant(w.index).hits["C"]
    fold = w.project.fold(w.index)
    col = int(np.floor((hit.t1 - fold.t1[0]) / fold.spec.period))
    w._point_clicked(col, float(hit.t2))
    assert ("S2", "C") in w.quant.decisions.assigned
    qtbot.waitUntil(lambda: not w.quant.running, timeout=TIMEOUT)
    assert w.quant.run_quant(w.index).hits["C"].manual
    assert _row(w, "S2", "C").amount == pytest.approx(1.5, rel=0.1)
    w.undo.undo()
    qtbot.waitUntil(lambda: not w.quant.running, timeout=TIMEOUT)
    assert ("S2", "C") not in w.quant.decisions.assigned


def test_quant_reintegrate_by_drawn_footprint_and_undo(qtbot, qwin):
    w = qwin
    w.start_reintegrate("S2", "C")
    qtbot.waitUntil(lambda: w._reintegrate == ("S2", "C") and w.view2d.footprint_roi.isVisible(), timeout=TIMEOUT)
    assert w.act_accept_footprint.isEnabled()
    auto = w.quant.run_quant(w.index).hits["C"]
    poly = w.view2d.footprint_vertices()
    assert len(poly) == 4
    # Shrink the polygon to its middle third in t1: less of the peak is inside, so the volume falls.
    t1a, t1b = poly[0][0], poly[1][0]
    third = (t1b - t1a) / 3.0
    shrunk = [(t1a + third, poly[0][1]), (t1b - third, poly[1][1]), (t1b - third, poly[2][1]), (t1a + third, poly[3][1])]
    w.view2d.start_footprint_edit(shrunk)
    w.accept_footprint()
    assert w._reintegrate is None and not w.act_accept_footprint.isEnabled()
    assert ("S2", "C") in w.quant.decisions.footprints
    qtbot.waitUntil(lambda: not w.quant.running, timeout=TIMEOUT)
    hand = w.quant.run_quant(w.index).hits["C"]
    assert hand.integrated_by_hand and hand.quant_volume < auto.quant_volume
    w.undo.undo()
    qtbot.waitUntil(lambda: not w.quant.running, timeout=TIMEOUT)
    assert ("S2", "C") not in w.quant.decisions.footprints
    assert not w.quant.run_quant(w.index).hits["C"].integrated_by_hand


def test_reintegrate_can_be_cancelled_and_reset(qtbot, qwin):
    w = qwin
    w.start_reintegrate("S2", "C")
    qtbot.waitUntil(lambda: w._reintegrate is not None, timeout=TIMEOUT)
    w.cancel_pick()
    assert w._reintegrate is None and not w.view2d.footprint_roi.isVisible()
    assert ("S2", "C") not in w.quant.decisions.footprints


def test_quant_settings_stale_and_export(qtbot, qwin, tmp_path):
    from dataclasses import replace
    w = qwin
    w.export_results(str(tmp_path / "out"))
    assert sorted(p.name for p in (tmp_path / "out").iterdir()) == [
        "all_signals.csv", "calibrations.csv", "internal_standards.csv", "targets.csv"]
    w.quant.set_settings(replace(w.quant.settings, qual_tolerance=10.0))
    assert w.quant.stale
    assert w.results.banner.isVisible() or not w.results.isVisible()
    w.results.view.setCurrentIndex(w.results.view.findData("istd"))
    assert not w.results.chart.isHidden()


def test_quant_saved_and_rerun(qtbot, qwin, tmp_path):
    from pygcgc.gui.mainwindow import MainWindow
    w = qwin
    w.decide([("S1", "B")], "rejected")
    w.exclude_point("C", "CAL_L1", "low")
    path = tmp_path / "q.pygcgc.json"
    w.project.path = path
    w.save_project()
    w2 = MainWindow()
    qtbot.addWidget(w2)
    done = []
    w2.open_project(path, on_loaded=lambda: done.append(True))
    qtbot.waitUntil(lambda: bool(done) and w2.quant.report is not None, timeout=TIMEOUT)
    assert ("S1", "B") in w2.quant.decisions.rejected
    assert w2.quant.decisions.excluded["C"]["CAL_L1"] == "low"
    assert next(r for r in w2.quant.report.results if r.run == "S1" and r.compound == "B").status == "rejected"
    w2.jobs.wait(5000)


# ------------------------------------------------------------------ 📍 align
def _peak_pixel(w, index, ion):
    f = w.project.fold(index)
    run = w.project.display_run(index)
    img = f.ion_image(run, ion)
    r, c = np.unravel_index(np.argmax(img), img.shape)
    return int(c), float(r * f.row_pitch)


def _setup_alignment(w):
    """Reference RUN_A, one good shift anchor on RUN_B (its peaks sit 1.5 s
    later in t1 than RUN_A's, per the synthetic_folder fixture)."""
    w.set_reference("RUN_A")
    i_a, i_b = w.project.index_of("RUN_A"), w.project.index_of("RUN_B")
    col_ref, t2_ref = _peak_pixel(w, i_a, 50)
    col_run, t2_run = _peak_pixel(w, i_b, 50)
    f_ref, f_run = w.project.fold(i_a), w.project.fold(i_b)
    ref_t1 = float(f_ref.t1[col_ref] + 0.5 * f_ref.spec.period)
    run_t1 = float(f_run.t1[col_run] + 0.5 * f_run.spec.period)
    w.align.add_anchor("RUN_B", "P1", ref_t1, t2_ref, run_t1, t2_run, 50)


def test_align_reference_and_anchor_click_flow(qtbot, window):
    w = window
    w.set_ion_expression("50")
    w.set_reference("RUN_A")
    assert w.align.reference == "RUN_A"
    assert w.align_panel.reference.currentData() == "RUN_A"   # combo reflects it, not just the state
    i_a, i_b = w.project.index_of("RUN_A"), w.project.index_of("RUN_B")
    w.set_run(i_b)
    w.start_add_anchor()
    assert w._anchor_stage == "ref"
    assert w.index == i_a
    _wait_render(qtbot, w, "RUN_A")
    col_ref, t2_ref = _peak_pixel(w, i_a, 50)
    qtbot.waitUntil(lambda: w.view2d.result.expr == w.renderer.expr, timeout=TIMEOUT)
    w._point_clicked(col_ref, t2_ref)
    assert w._anchor_stage == "run"
    assert w.index == i_b
    _wait_render(qtbot, w, "RUN_B")
    col_run, t2_run = _peak_pixel(w, i_b, 50)
    w._point_clicked(col_run, t2_run)
    assert w._anchor_stage == "run" and w.project.entries[w.index].name == "RUN_C_BLANK"
    w.align_panel.skip_btn.click()
    assert w._anchor_stage == "ref" and w.index == i_a
    assert w._anchor_pending["name"] == "P2"
    assert "Choose a new peak" in w.align_panel.pick_status.text()
    assert w.align_panel.table.rowCount() == 2
    assert w.align.points["RUN_A"][0].name == "P1"
    assert [p.name for p in w.align.points["RUN_B"]] == ["P1"]
    assert w.align.points["RUN_B"][0].ion == 50
    tf = w.align.transforms["RUN_B"]
    assert tf.kind == "shift"
    # The fit should exactly reproduce the two points that were actually
    # clicked (each click snaps to its modulation's mid-time, so with only
    # one anchor this need not match the compounds' true 1.5 s offset).
    p = w.align.points["RUN_B"][0]
    t1f, t2f = tf.forward(p.run_t1, p.run_t2)
    assert (t1f, t2f) == pytest.approx((p.ref_t1, p.ref_t2), abs=1e-6)
    w.cancel_pick()


def test_align_add_anchor_undo(qtbot, window):
    w = window
    w.set_ion_expression("50")
    _setup_alignment(w)
    w.set_run(w.project.index_of("RUN_B"))
    n = w.undo.count()
    assert n == 1        # only set_reference went through the mainwindow (add_anchor was called directly)
    w.remove_anchor("P1")
    assert "RUN_B" not in w.align.points
    assert w.undo.count() == n + 1
    w.undo.undo()
    assert "RUN_B" in w.align.transforms   # the direct add is still there; only the remove is undone
    w.undo.redo()
    assert "RUN_B" not in w.align.points


def test_align_kind_override_and_error(qtbot, window):
    w = window
    w.set_ion_expression("50")
    _setup_alignment(w)
    w.set_run(w.project.index_of("RUN_B"))
    w.set_align_kind("affine")               # only one anchor: affine can't be fit
    assert "RUN_B" in w.align.errors
    assert "RUN_B" not in w.align.transforms
    w.undo.undo()
    assert w.align.kind.get("RUN_B", "auto") == "auto"
    assert "RUN_B" in w.align.transforms


def test_align_remove_and_clear_with_undo(qtbot, window):
    w = window
    w.set_ion_expression("50")
    _setup_alignment(w)
    w.set_run(w.project.index_of("RUN_B"))
    w.clear_run_anchors()
    assert "RUN_B" not in w.align.points
    assert "RUN_B" not in w.align.transforms
    w.undo.undo()
    assert "RUN_B" in w.align.transforms


def test_align_propagate_finds_and_flags(qtbot, window):
    w = window
    w.set_ion_expression("50")
    _setup_alignment(w)
    messages = []
    w.project.message.connect(messages.append)
    w.propagate_anchors()
    qtbot.waitUntil(lambda: w.jobs.active_count() == 0 and "RUN_B" in w.align.transforms, timeout=TIMEOUT)
    assert "RUN_C_BLANK" not in w.align.points
    assert any("RUN_C_BLANK" in m for m in messages)   # its peak is nowhere near P1's window


def test_aligned_view_shows_resampled_image(qtbot, window):
    w = window
    w.set_ion_expression("50")
    _setup_alignment(w)
    i_a, i_b = w.project.index_of("RUN_A"), w.project.index_of("RUN_B")
    w.set_run(i_b)
    qtbot.waitUntil(lambda: w.act_aligned.isEnabled(), timeout=TIMEOUT)
    w.act_aligned.setChecked(True)
    qtbot.waitUntil(lambda: w.view2d.result is not None and w.view2d.result.aligned, timeout=TIMEOUT)
    ref_img = w.project.fold(i_a).ion_image(w.project.display_run(i_a), 50)
    _, c_ref = np.unravel_index(np.argmax(ref_img), ref_img.shape)
    _, c_aligned = np.unravel_index(np.argmax(w.view2d.result.image), w.view2d.result.image.shape)
    assert abs(int(c_ref) - int(c_aligned)) <= 1
    w.act_aligned.setChecked(False)
    qtbot.waitUntil(lambda: not w.view2d.result.aligned, timeout=TIMEOUT)


def _open_compare(w, qtbot, a="RUN_A", b="RUN_B"):
    w.act_compare.trigger()
    cw = w.compare_window
    cw.set_pair(a, b)
    qtbot.wait(80)
    return cw


def test_compare_view_modes(qtbot, window):
    w = window
    w.set_ion_expression("50")
    _setup_alignment(w)
    cw = _open_compare(w, qtbot)
    view = cw.images
    assert cw.isVisible() and "Aligned" in cw.align_label.text()
    assert view.current_mode() == "blink"
    assert view.info.text() == ""
    view.mode.setCurrentIndex(view.mode.findData("overlay"))
    assert view.single.image.image is not None
    view.mode.setCurrentIndex(view.mode.findData("difference"))
    assert "robust brightness match" in view.info.text()
    view.mode.setCurrentIndex(view.mode.findData("grid"))
    assert view.grid_canvas.isVisible()
    assert not view.blink_timer.isActive()


def test_compare_works_without_alignment(qtbot, window):
    # The Compare tab used to show nothing until a 📍 reference was set and the run aligned.
    w = window
    w.set_ion_expression("50")
    assert w.align.reference is None
    cw = _open_compare(w, qtbot)
    assert "Not aligned" in cw.align_label.text()
    assert cw.images.info.text() == ""
    cw.show_images("overlay")
    rgb = cw.images.single.image.image
    assert rgb is not None and rgb.ndim == 3 and rgb[..., 1].max() > 0.5   # B (green) is there


def test_compare_uses_alignment_when_both_runs_have_it(qtbot, window):
    w = window
    w.set_ion_expression("50")
    _setup_alignment(w)
    cw = _open_compare(w, qtbot)
    i_a, i_b = w.project.index_of("RUN_A"), w.project.index_of("RUN_B")
    img_a, img_b = cw._image(i_a)[0], cw._image(i_b)[0]
    _, c_a = np.unravel_index(np.argmax(img_a), img_a.shape)
    _, c_b = np.unravel_index(np.argmax(img_b), img_b.shape)
    assert abs(int(c_a) - int(c_b)) <= 1


def test_blink_keeps_one_colour_scale(qtbot, window):
    # Each blink frame used to rescale itself, hiding a difference in amount (RUN_B holds 0.7 of RUN_A).
    w = window
    w.set_ion_expression("50")
    cw = _open_compare(w, qtbot)
    cw.show_images("blink")
    view = cw.images
    first = tuple(view.single.image.levels)
    view._blink_tick()
    assert tuple(view.single.image.levels) == pytest.approx(first)


def test_compare_overlay_after_difference_is_not_washed_out(qtbot, window):
    # Difference sets scalar (min, max) levels on the shared ImageItem; switching
    # to overlay afterwards must not leave those applied to the new RGB image
    # (which flattened it to uniform grey before this was fixed).
    w = window
    w.set_ion_expression("50")
    _setup_alignment(w)
    cw = _open_compare(w, qtbot)
    view = cw.images
    view.mode.setCurrentIndex(view.mode.findData("difference"))
    view.mode.setCurrentIndex(view.mode.findData("overlay"))
    rgb = view.single.image.image
    assert rgb is not None and rgb.ndim == 3
    assert rgb[..., :2].std() > 0.01   # varies; not flattened to one grey value


def test_blink_shortcut_opens_compare_in_blink_mode(qtbot, window):
    w = window
    w.set_ion_expression("50")
    from PySide6 import QtWidgets
    assert "⚖️ Compare" in [a.text() for a in w.findChild(QtWidgets.QToolBar, "toolbar_main").actions()]
    w.act_blink.trigger()
    assert w.compare_window.isVisible()
    assert w.compare_window.images.current_mode() == "blink"
    assert w.tabs.count() == 4                     # the old Compare tab is gone from the main window


def test_compare_pixel_spectra_and_similarity_map(qtbot, window):
    w = window
    cw = _open_compare(w, qtbot)
    cw.tabs.setCurrentIndex(1)
    i_a = w.project.index_of("RUN_A")
    f = w.project.fold(i_a)
    t1, t2, ion, _ = GUI_PEAKS[1]
    col = int((t1 - f.t1[0]) // f.spec.period)
    assert cw.show_pixel(float(f.t1[col]), t2) == pytest.approx(1.0, abs=0.05)   # same compound in both
    cw.compute_similarity()
    qtbot.waitUntil(lambda: cw.sim is not None, timeout=TIMEOUT)
    row = int(round(t2 / f.row_pitch))
    assert cw.sim.similarity[row, col] > 0.9
    assert cw.sim_image.image is not None and cw.sim_bar.isVisible()


def test_compare_class_shares(qtbot, window):
    w = window
    cw = _open_compare(w, qtbot)
    cw.tabs.setCurrentIndex(3)
    assert "Classify both runs" in cw.class_status.text()
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    for name in ("RUN_A", "RUN_B"):
        i = w.project.index_of(name)
        w.set_run(i)
        w.start_classify()
        qtbot.waitUntil(lambda: w.classes.result_for(i) is not None and not w.classes.running, timeout=TIMEOUT)
    cw.tabs.setCurrentIndex(0)
    cw.tabs.setCurrentIndex(3)
    t = cw.class_table
    assert t.rowCount() >= 1
    assert sum(float(t.item(r, 3).text()) for r in range(t.rowCount())) == pytest.approx(100.0, abs=0.5)


def test_compare_pairs_peaks_by_identity_and_position(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    b = w.blobs.blobs_for(w.project.index_of("RUN_B"))
    k = int(np.argmax(b.volume))
    w.project.identifications["RUN_B"] = [(float(b.t1[k]), float(b.t2[k]), "compound X")]
    a = w.blobs.blobs_for(w.project.index_of("RUN_A"))
    ka = int(np.argmin(np.hypot((a.t1 - b.t1[k]) / 60.0, (a.t2 - b.t2[k]) / 0.5)))
    w.project.identifications["RUN_A"] = [(float(a.t1[ka]), float(a.t2[ka]), "compound X")]
    cw = _open_compare(w, qtbot)
    cw.tabs.setCurrentIndex(2)
    cw.match()
    qtbot.waitUntil(lambda: bool(cw._current_pairs()), timeout=TIMEOUT)
    pairs = cw._current_pairs()
    named = [p for p in pairs if p.name == "compound X"]
    assert len(named) == 1 and named[0].how == "name"
    # The fixture's real peaks (the bleed on m/z 54 also gives many small noise blobs)
    big = [p for p in pairs if p.status == "both" and p.volume_a > 1e5]
    assert len(big) >= len(GUI_PEAKS) - 1
    for p in big:
        assert abs(p.dt1) < 4.0 and abs(p.dt2) < 0.05 and p.similarity > 0.99
        assert np.log2(0.6) < p.log2_ratio < 0     # RUN_B was made at 0.7 × RUN_A, plus the same noise
    assert cw.peak_table.rowCount() == len(pairs)
    cw.filter_box.setCurrentIndex(cw.filter_box.findData("named"))
    assert cw.peak_table.rowCount() == 1


def test_align_saved_and_reloaded(qtbot, window, tmp_path):
    from pygcgc.gui.mainwindow import MainWindow
    w = window
    w.set_ion_expression("50")
    _setup_alignment(w)
    path = tmp_path / "align.pygcgc.json"
    w.project.path = path
    w.save_project()

    w2 = MainWindow()
    qtbot.addWidget(w2)
    done = []
    w2.open_project(path, on_loaded=lambda: done.append(True))
    qtbot.waitUntil(lambda: bool(done), timeout=TIMEOUT)
    assert w2.align.reference == "RUN_A"
    assert "RUN_B" in w2.align.transforms
    assert w2.align.points["RUN_B"][0].ion == 50
    w2.jobs.wait(5000)


def test_reference_change_clears_and_is_undoable(qtbot, window):
    w = window
    w.set_ion_expression("50")
    _setup_alignment(w)
    assert w.align.points
    w.set_reference("RUN_C_BLANK")
    assert w.align.reference == "RUN_C_BLANK"
    assert w.align.points == {}
    w.undo.undo()
    assert w.align.reference == "RUN_A"
    assert "RUN_B" in w.align.transforms


# ------------------------------------------------------------------ 🔍 blobs
def test_find_blobs_runs_and_fills_panel(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    for e in w.project.entries:
        assert e in w.blobs.by_run
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    blobs = w.blobs.blobs_for(i)
    assert len(blobs) >= 3                       # GUI_PEAKS has 5 well-separated compounds
    assert w.blobs_panel.table.rowCount() == len(blobs)


def test_ms_aware_detection_is_a_first_class_peak_method(qtbot, window, tmp_path):
    from pygcgc.gui.blob_state import MS_AWARE

    w = window
    w.blobs_panel.method_box.setCurrentIndex(w.blobs_panel.method_box.findData(MS_AWARE))
    assert w.blobs.method == MS_AWARE
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.candidates_were_run,
                    timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    candidates = w.blobs.candidates_for(i)
    regions = w.blobs.regions_for(i)
    assert candidates is not None and len(candidates) > 0
    assert regions
    assert "overlapping blobs" in w.blobs_panel.status.text()

    w.act_blob_overlay.setChecked(True)
    automatic = w.blobs.components_for(i)
    assert automatic is not None and automatic.groups
    assert w.view2d.component_items
    w.blobs_panel.select_row(0)
    w.resolve_region(0)
    assert w.overlap_dialog._source == "candidates"
    assert w.overlap_dialog.cube is not None

    dialog = w.overlap_dialog
    dialog.auto_k.setChecked(False)
    dialog.k_spin.setValue(1)
    dialog._run_fit()
    qtbot.waitUntil(lambda: dialog.accept_fit_btn.isEnabled(), timeout=TIMEOUT)
    dialog.accept_fit_btn.click()
    accepted = w.blobs.components_for(i)
    assert accepted is not None and accepted.group_for_region(0) is not None
    assert len(accepted.group_for_region(0).components) == 1
    assert w.view2d.component_items
    saved = w.blobs.to_dict()
    assert saved["component_sets"]["RUN_A"]["groups"][0]["components"]
    from pygcgc.gui.blob_state import BlobState
    restored = BlobState(w.project)
    restored.from_dict(saved)
    assert restored.component_sets_by_name["RUN_A"].group_for_region(0) is not None
    assert 0 in restored.accepted_groups_by_name["RUN_A"]
    assert 0 in w.blobs.accepted_groups_by_name["RUN_A"]
    w.undo.undo()
    assert 0 not in w.blobs.accepted_groups_by_name.get("RUN_A", {})
    w.undo.redo()
    assert 0 in w.blobs.accepted_groups_by_name["RUN_A"]
    export_path = tmp_path / "components.csv"
    w.export_components(export_path)
    assert export_path.exists()
    assert "group_unresolved" in export_path.read_text(encoding="utf-8")


def test_clicking_an_ms_aware_component_shows_its_own_full_ion_spectrum(qtbot, window):
    from pygcgc.gui.blob_state import MS_AWARE

    w = window
    w.blobs_panel.method_box.setCurrentIndex(w.blobs_panel.method_box.findData(MS_AWARE))
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.components_for(i) is not None
                    and bool(w.blobs.components_for(i).components), timeout=TIMEOUT)
    w.act_blob_overlay.setChecked(True)
    component_set = w.blobs.components_for(i)
    component = component_set.components[0]
    # click on the blob itself, as the 2D view reports a click: its position plus pixel
    w.view2d.last_click = (component.t1 / 60.0, component.t2)
    row, col = w.view2d.pixel_at(*w.view2d.last_click)
    w._point_clicked(int(col), row * w.project.fold(i).row_pitch)
    assert w.blobs_panel.selected_component_ids() == {component.id}
    qtbot.waitUntil(lambda: w.spectrum.mz.size > 0, timeout=TIMEOUT)
    assert "fitted peak" in w.spectrum.title
    assert f"{component.t1 / 60:.3f}" in w.spectrum.title
    # re-solved on every run ion, not truncated to the component's own ~30 diagnostic ions
    assert w.spectrum.mz.size >= component.mz.size
    run = w.project.base_run(i)
    assert set(w.spectrum.mz.astype(int).tolist()) <= set(run.mz.astype(int).tolist())


def test_new_project_defaults_to_automatic_overlapping_blobs(window):
    from pygcgc.gui.blob_state import BlobState, MS_AWARE
    from pygcgc.gui.panels.blobs import BlobsPanel

    fresh = BlobState(window.project)
    assert fresh.method == MS_AWARE
    assert fresh.component_settings.automatically_separate
    fresh.from_dict({})
    assert fresh.method == MS_AWARE
    assert BlobsPanel().method == MS_AWARE


def test_peak_detection_method_is_saved(window):
    from pygcgc.gui.blob_state import MS_AWARE

    window.blobs.set_method(MS_AWARE)
    saved = window.blobs.to_dict()
    window.blobs.set_method("legacy")
    window.blobs.from_dict(saved)
    assert window.blobs.method == MS_AWARE


def test_blob_overlay_and_click_selects_row(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.act_blob_overlay.setChecked(True)
    blobs = w.blobs.blobs_for(i)
    x, _y = w.view2d._blob_outline
    assert len(x) > 0                            # something is actually drawn
    col, t2 = _peak_pixel(w, i, 50)
    w._point_clicked(col, t2)
    row = int(np.argmin(np.hypot(
        (blobs.t1 - (w.project.fold(i).t1[col] + 0.5 * w.project.fold(i).spec.period)) / 60.0,
        (blobs.t2 - t2) / 0.5)))
    assert w.blobs_panel.selected_indices() == [row]
    # the clicked blob's own spectrum is shown, not just whatever pixel was under the cursor
    expected = blobs.apex_spectra[row] if (
        blobs.apex_spectra.shape[0] == len(blobs)
        and (blobs.apex_spectra[row] > 0).sum() >= 5) else blobs.spectra[row]
    np.testing.assert_array_equal(w.spectrum.mz, blobs.mz[expected > 0])
    np.testing.assert_array_equal(w.spectrum.intensity, expected[expected > 0])
    assert f"{blobs.t1[row] / 60:.3f}" in w.spectrum.title


def test_blob_delete_is_undoable(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    n = len(w.blobs.blobs_for(i))
    w.delete_blobs([0])
    assert len(w.blobs.blobs_for(i)) == n - 1
    w.undo.undo()
    assert len(w.blobs.blobs_for(i)) == n
    w.undo.redo()
    assert len(w.blobs.blobs_for(i)) == n - 1


def test_blob_merge_is_exact_and_undoable(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    before = w.blobs.blobs_for(i)
    n = len(before)
    total_volume = float(before.volume[[0, 1]].sum())
    w.merge_blobs([0, 1])
    after = w.blobs.blobs_for(i)
    assert len(after) == n - 1
    assert after.volume[-1] == pytest.approx(total_volume)
    w.undo.undo()
    assert len(w.blobs.blobs_for(i)) == n


def test_ctrl_click_on_footprints_builds_a_selection_to_merge(qtbot, window):
    from PySide6 import QtCore
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.view2d.blob_overlay_on = True
    blobs = w.blobs.blobs_for(i)
    n = len(blobs)
    # Click a footprint pixel well away from each apex, as on a tail fragment.
    clicks = []
    for b in (0, 1):
        rows, cols = blobs.footprint(b)
        far = int(np.argmax(np.hypot(rows - blobs.row[b], cols - blobs.col[b])))
        clicks.append((int(cols[far]), float(rows[far] * blobs.pitch)))
    w.view2d.last_click_modifiers = QtCore.Qt.KeyboardModifier.NoModifier
    w._select_blob_at(*clicks[0])
    w.view2d.last_click_modifiers = QtCore.Qt.KeyboardModifier.ControlModifier
    w._select_blob_at(*clicks[1])
    assert w.blobs_panel.selected_indices() == [0, 1]
    w.merge_selected_blobs()
    assert len(w.blobs.blobs_for(i)) == n - 1


def test_blob_add_recovers_a_deleted_peak(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    blobs = w.blobs.blobs_for(i)
    n = len(blobs)
    col, t2 = _peak_pixel(w, i, 50)              # a real, well-separated peak
    fold = w.project.fold(i)
    t1_mid = float(fold.t1[col] + 0.5 * fold.spec.period)
    w.delete_blobs([blobs.nearest(t1_mid, t2)])
    assert len(w.blobs.blobs_for(i)) == n - 1
    w.start_add_blob()
    assert w._blob_add_pending
    w._point_clicked(col, t2)
    assert not w._blob_add_pending
    assert len(w.blobs.blobs_for(i)) == n
    w.undo.undo()
    assert len(w.blobs.blobs_for(i)) == n - 1


def test_blobs_saved_and_reloaded(qtbot, window, tmp_path):
    from pygcgc.gui.mainwindow import MainWindow
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    n = len(w.blobs.blobs_for(i))
    w.delete_blobs([0])
    path = tmp_path / "blobs.pygcgc.json"
    w.project.path = path
    w.save_project()

    w2 = MainWindow()
    qtbot.addWidget(w2)
    done = []
    w2.open_project(path, on_loaded=lambda: done.append(True))
    qtbot.waitUntil(lambda: bool(done) and w2.blobs.was_run and not w2.blobs.running, timeout=TIMEOUT)
    i2 = w2.project.index_of("RUN_A")
    assert len(w2.blobs.blobs_for(i2)) == n - 1     # the hand deletion survived the reload
    w2.jobs.wait(5000)


# ------------------------------------------------------------ ⛓️ hydrocarbons
@pytest.fixture
def hydro_window(qtbot, alkane_folder):
    from pygcgc.gui.mainwindow import MainWindow
    w = MainWindow()
    # Hydrocarbon mapping still consumes legacy BlobSet features, so these tests
    # choose that method explicitly now that fresh projects start MS-aware.
    w.blobs.set_method("legacy")
    w.jobs_panel.clear_after_ms = 10
    qtbot.addWidget(w)
    w.resize(1400, 900)
    w.show()
    loaded = []
    w.open_folder(alkane_folder, on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and w.jobs.active_count() == 0 and w.view2d.result is not None,
                    timeout=TIMEOUT)
    # No target table (no anchors), so phase isn't estimated automatically
    # on load; the quiet-band (void) method needs an explicit 🔄 Fix phase.
    w.fix_phase()
    qtbot.waitUntil(lambda: all(e.phase_method != "unset" for e in w.project.entries), timeout=TIMEOUT)
    yield w
    w.jobs.cancel_all()
    w.jobs.wait(5000)
    w.undo.clear()


def _click_carbon(w, monkeypatch, t1_abs, t2_abs, carbon):
    from PySide6 import QtWidgets
    # The fixture's t2 values are in the modulator's own frame, so add its true
    # phase to get the peak's absolute time; the fold frame depends on the void cut.
    col, t2 = w.project.fold(w.index).to_fold(t1_abs + GUI_PHASE + t2_abs)
    monkeypatch.setattr(QtWidgets.QInputDialog, "getInt",
                        staticmethod(lambda *a, **k: (carbon, True)))
    w._point_clicked(int(col), float(t2))


def test_find_alkane_series_two_click_flow(qtbot, hydro_window, monkeypatch):
    from tests.conftest import ALKANE_CARBONS, ALKANE_STEP, ALKANE_T1_START
    w = hydro_window
    i_std = w.project.index_of("STD")
    w.set_run(i_std)
    w.start_find_alkanes()
    assert w._alkane_seed_pending == []
    i12, i20 = ALKANE_CARBONS.index(12), ALKANE_CARBONS.index(20)
    _click_carbon(w, monkeypatch, ALKANE_T1_START + i12 * ALKANE_STEP, 1.0 + 0.01 * i12, 12)
    assert w._alkane_seed_pending is not None and len(w._alkane_seed_pending) == 1
    _click_carbon(w, monkeypatch, ALKANE_T1_START + i20 * ALKANE_STEP, 1.0 + 0.01 * i20, 20)
    qtbot.waitUntil(lambda: not w.hydrocarbons.running, timeout=TIMEOUT)
    assert w._alkane_seed_pending is None
    series = w.hydrocarbons.series_for(i_std)
    assert series is not None and len(series) == len(ALKANE_CARBONS)
    assert w.hydro_panel.series_table.rowCount() == len(ALKANE_CARBONS)


def _find_std_series(w, qtbot, monkeypatch):
    from tests.conftest import ALKANE_CARBONS, ALKANE_STEP, ALKANE_T1_START
    i_std = w.project.index_of("STD")
    w.set_run(i_std)
    w.start_find_alkanes()
    for c in (12, 20):
        i = ALKANE_CARBONS.index(c)
        _click_carbon(w, monkeypatch, ALKANE_T1_START + i * ALKANE_STEP, 1.0 + 0.01 * i, c)
    qtbot.waitUntil(lambda: not w.hydrocarbons.running, timeout=TIMEOUT)
    return i_std


def _click_on_hit(w, hit):
    """A 2D-view click exactly on a found alkane's apex pixel."""
    fold = w.project.fold(w.index)
    col = int(np.floor((hit.t1 - fold.spec.t0 - fold.spec.phase) / fold.spec.period))
    w._point_clicked(col, float(hit.t2))


def test_siloxane_panel_presets_and_edits_reach_the_project(qtbot, hydro_window):
    from pygcgc.core.siloxane import BLEED_REFERENCE, CYCLIC_REFERENCE
    w = hydro_window
    assert w.project.siloxane_settings.reference_dict() == BLEED_REFERENCE        # measured default
    w.sil_panel.preset.setCurrentText("Cyclic siloxane peak (measured)")
    w.sil_panel._preset_chosen(0)
    assert w.project.siloxane_settings.reference_dict() == CYCLIC_REFERENCE
    # Edit one intensity by hand: the preset flips to Custom and the project follows.
    w.sil_panel.table.item(0, 1).setText("77")
    assert w.sil_panel.preset.currentText() == "Custom"
    assert w.project.siloxane_settings.reference_dict()[min(CYCLIC_REFERENCE)] == 77.0
    w.sil_panel.threshold.setValue(55)             # the slider is in per cent
    assert w.project.siloxane_settings.score_threshold == pytest.approx(0.55)
    w.sil_panel.method.setCurrentIndex(1)
    w.sil_panel.method.activated.emit(1)
    assert w.project.siloxane_settings.method == "pixels"


def test_siloxane_reference_can_be_picked_from_the_2d_view(qtbot, hydro_window, monkeypatch):
    w = hydro_window
    i_std = _find_std_series(w, qtbot, monkeypatch)
    hit = w.hydrocarbons.series_for(i_std).hits[15]
    w.sil_panel.pick_btn.setChecked(True)
    assert w._siloxane_pick_pending
    _click_on_hit(w, hit)
    assert not w._siloxane_pick_pending and not w.sil_panel.pick_btn.isChecked()
    ref = w.project.siloxane_settings.reference_dict()
    assert set(ref) == {57, 71, 85}                          # the ions the synthetic alkane carries
    assert max(ref.values()) == 100.0 and w.sil_panel.preset.currentText() == "Custom"
    assert w.sil_panel.table.rowCount() == 3


def test_esc_cancels_a_siloxane_pick(qtbot, hydro_window):
    w = hydro_window
    w.set_run(0)
    w.sil_panel.pick_btn.setChecked(True)
    assert w._siloxane_pick_pending
    w.cancel_pick()
    assert not w._siloxane_pick_pending and not w.sil_panel.pick_btn.isChecked()


def test_alkane_qualifier_ratios_are_editable_and_pickable(qtbot, hydro_window, monkeypatch):
    w = hydro_window
    assert w.hydrocarbons.settings.qual_ratios == pytest.approx((0.66, 0.44), abs=0.02)   # measured default
    assert w.hydro_panel.qualifiers()[0] == pytest.approx((0.66, 0.44), abs=0.02)
    w.hydro_panel._qual_spins[0].setValue(0.5)
    assert w.hydrocarbons.settings.qual_ratios[0] == pytest.approx(0.5)
    w.hydro_panel.qual_enforce.setChecked(False)
    assert not w.hydrocarbons.settings.enforce_qualifiers
    # Pick from a clicked alkane: the synthetic peaks carry the measured ratios, so they come back.
    w.hydro_panel.qual_enforce.setChecked(True)
    w.hydro_panel._qual_spins[0].setValue(0.66)
    i_std = _find_std_series(w, qtbot, monkeypatch)
    hit = w.hydrocarbons.series_for(i_std).hits[15]
    w.hydro_panel._qual_spins[0].setValue(0.1)
    w.hydro_panel.qual_pick.setChecked(True)
    assert w._qual_pick_pending
    _click_on_hit(w, hit)
    assert not w._qual_pick_pending
    assert w.hydrocarbons.settings.qual_ratios == pytest.approx((0.66, 0.44), abs=0.04)
    assert w.hydro_panel.qualifiers()[0] == pytest.approx((0.66, 0.44), abs=0.04)
    # ...and the reset button restores the shipped defaults.
    w.hydro_panel._qual_spins[0].setValue(0.2)
    w.hydro_panel.qual_reset.click()
    assert w.hydrocarbons.settings.qual_ratios == pytest.approx((0.66, 0.44), abs=0.02)


def test_series_table_shows_the_qualifier_ratios_per_hit(qtbot, hydro_window, monkeypatch):
    w = hydro_window
    _find_std_series(w, qtbot, monkeypatch)
    assert w.hydro_panel.series_table.columnCount() == 3 + len(w.hydrocarbons.settings.qual_ions)
    assert w.hydro_panel.series_table.horizontalHeaderItem(3).text() == "71/57"
    assert float(w.hydro_panel.series_table.item(0, 3).text()) == pytest.approx(0.66, abs=0.05)


def test_propagate_stencil_lift_and_undo(qtbot, hydro_window, monkeypatch):
    w = hydro_window
    i_std = _find_std_series(w, qtbot, monkeypatch)
    i_sample = w.project.index_of("SAMPLE")

    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)

    w.propagate_alkanes()
    qtbot.waitUntil(lambda: not w.hydrocarbons.running, timeout=TIMEOUT)
    std_series = w.hydrocarbons.series_for(i_std)
    sample_series = w.hydrocarbons.series_for(i_sample)
    assert sample_series is not None and len(sample_series) == len(std_series)

    box = [std_series.to_absolute(15, 0.15), std_series.to_absolute(25, 0.15),
          std_series.to_absolute(25, 0.45), std_series.to_absolute(15, 0.45)]
    w.hydro_panel.draw_btn.setChecked(True)
    w.view2d.stencil_roi.setPoints([[t1 / 60.0, t2] for t1, t2 in box])
    w.view2d.stencil_roi.setPos([0, 0])
    n = w.undo.count()
    w.hydro_panel.saveStencilRequested.emit("branched C15-25", 2.0)
    assert w.undo.count() == n + 1
    assert [s.name for s in w.hydrocarbons.stencils.stencils] == ["branched C15-25"]

    totals_std = w.hydrocarbons.totals(i_std)["branched C15-25"]
    assert totals_std["n_blobs"] == 10

    w.set_run(i_sample)
    totals_sample = w.hydrocarbons.totals(i_sample)["branched C15-25"]
    assert totals_sample["n_blobs"] == totals_std["n_blobs"]         # the lift found the same compounds
    assert totals_sample["volume"] == pytest.approx(totals_std["volume"], rel=0.1)
    assert totals_sample["amount"] == pytest.approx(totals_sample["volume"] * 2.0)

    w.undo.undo()
    assert w.hydrocarbons.stencils.stencils == []
    w.undo.redo()
    assert [s.name for s in w.hydrocarbons.stencils.stencils] == ["branched C15-25"]


def test_stencil_shape_can_be_edited_and_replaced(qtbot, hydro_window, monkeypatch):
    w = hydro_window
    i_std = _find_std_series(w, qtbot, monkeypatch)
    series = w.hydrocarbons.series_for(i_std)
    box = [series.to_absolute(15, 0.15), series.to_absolute(25, 0.15),
           series.to_absolute(25, 0.45), series.to_absolute(15, 0.45)]
    w.hydro_panel.draw_btn.setChecked(True)
    w.view2d.stencil_roi.setPoints([[t1 / 60.0, t2] for t1, t2 in box])
    w.view2d.stencil_roi.setPos([0, 0])
    w.hydro_panel.saveStencilRequested.emit("band", 2.0)
    before = list(w.hydrocarbons.stencils.stencils[0].vertices)
    # Edit: the polygon comes back into the ROI at the same place, and a wider one replaces it.
    w.hydro_panel.editStencilRequested.emit("band")
    assert w.view2d.stencil_vertices() == pytest.approx(box, abs=1e-6)
    assert w.hydro_panel.rf_spin.value() == 2.0
    wider = [series.to_absolute(15, 0.10), series.to_absolute(25, 0.10),
             series.to_absolute(25, 0.60), series.to_absolute(15, 0.60)]
    w.view2d.stencil_roi.setPoints([[t1 / 60.0, t2] for t1, t2 in wider])
    w.view2d.stencil_roi.setPos([0, 0])
    w.hydro_panel.saveStencilRequested.emit("band", 2.0)
    assert [s.name for s in w.hydrocarbons.stencils.stencils] == ["band"]      # replaced, not duplicated
    assert w.hydrocarbons.stencils.stencils[0].vertices != before


def test_homologous_series_button_finds_the_branched_band(qtbot, hydro_window, monkeypatch):
    w = hydro_window
    i_std = _find_std_series(w, qtbot, monkeypatch)
    w.hydro_panel.find_homolog_btn.click()
    assert "blobs" in w.hydro_panel.homolog_status.text()          # needs 🔍 first
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    w.hydro_panel.pair_combo.setCurrentText("43/57 branched")
    w.hydro_panel.find_homolog_btn.click()
    assert w.hydro_panel.homolog_table.rowCount() >= 4
    assert w.hydro_panel.homolog_status.text().startswith(f"{w.hydro_panel.homolog_table.rowCount()} members")
    w.hydro_panel.pair_combo.setCurrentText("bad")
    w.hydro_panel.find_homolog_btn.click()
    assert "major/minor" in w.hydro_panel.homolog_status.text()


def test_settings_dialog_edits_any_settings_dataclass(qtbot):
    from pygcgc.core import HydrocarbonSettings, IdentifySettings, PeakSettings
    from pygcgc.core.batch import BatchSettings
    from pygcgc.gui.dialogs import SettingsDialog
    for cls in (PeakSettings, HydrocarbonSettings, IdentifySettings):
        base = cls()
        dlg = SettingsDialog("t", base)
        qtbot.addWidget(dlg)
        assert dlg.settings() == base                                # untouched form round-trips exactly
    dlg = SettingsDialog("peaks", PeakSettings())
    qtbot.addWidget(dlg)
    getters = dict(dlg._getters)
    assert ("snr_min",) in getters and ("ion_snr",) in getters
    dlg.findChildren(__import__("PySide6").QtWidgets.QDoubleSpinBox)[0].setValue(33.0)
    assert dlg.settings().snr_min == 33.0
    nested = SettingsDialog("batch", BatchSettings())                # nested dataclasses get their own paths
    qtbot.addWidget(nested)
    assert ("peaks", "snr_min") in dict(nested._getters) and ("quant", "peaks", "snr_min") in dict(nested._getters)
    assert nested.settings() == BatchSettings()
    hs = SettingsDialog("h", HydrocarbonSettings(qual_ratios=(0.5, 0.3)))
    qtbot.addWidget(hs)
    assert hs.settings().qual_ratios == (0.5, 0.3)                    # tuples survive as text


def test_blob_and_hydro_settings_dialogs_reach_state(qtbot, hydro_window, monkeypatch):
    from pygcgc.gui import dialogs
    w = hydro_window

    class Fake:
        def __init__(self, *a, **k):
            self.base = a[1]

        def exec(self):
            return True

        def settings(self):
            from dataclasses import replace
            return replace(self.base, **({"snr_min": 7.0} if hasattr(self.base, "snr_min") else {"max_gap": 5}))

    monkeypatch.setattr("pygcgc.gui.mainwindow.SettingsDialog", Fake)
    w.edit_blob_settings()
    assert w.blobs.settings.snr_min == 7.0
    w.edit_hydro_settings()
    assert w.hydrocarbons.settings.max_gap == 5


def test_export_figure_writes_png_and_caption(qtbot, hydro_window, tmp_path):
    w = hydro_window
    w.set_run(0)
    qtbot.waitUntil(lambda: w.view2d.result is not None, timeout=TIMEOUT)
    out = tmp_path / "fig_test.png"
    w.save_view_figure(out, 89.0)
    assert out.exists() and out.stat().st_size > 2000
    caption = out.with_suffix(".txt").read_text(encoding="utf-8")
    assert "modulation period" in caption and "GC" in caption


def test_toolbar_follows_the_working_order(qtbot, hydro_window):
    w = hydro_window
    # Step buttons carry a ✓ (done) or ▸ (suggested next) prefix.
    tb = [a.text().lstrip("✓▸ ") for a in w.findChildren(__import__("PySide6").QtWidgets.QToolBar)[0].actions()
          if a.text()]
    order = ["1 Load", "2 Align", "3 🧹 Background", "4 🔍 Peaks", "5 Standards",
             "6 Assign", "7 🧪 Quantify", "8 Export"]
    pos = [tb.index(x) for x in order if x in tb]
    assert len(pos) == len(order) and pos == sorted(pos)
    # Siloxane and hydrocarbons open from the Peaks panel, not the toolbar.
    assert not any("Siloxane" in x or "Hydrocarbon" in x for x in tb)


def test_siloxane_and_hydrocarbon_tools_open_from_peaks(qtbot, window):
    w = window
    wf = w.workflow
    w.step_actions["blobs"].trigger()
    w.blobs_panel.siloxane_btn.click()
    assert wf.current() == "siloxane" and w.step_actions["blobs"].isChecked()
    assert wf.title.text() == "Siloxane filterer" and wf.back_btn.isVisibleTo(wf)
    wf.back_btn.click()
    assert wf.current() == "blobs" and not wf.back_btn.isVisibleTo(wf)
    w.blobs_panel.hydrocarbon_btn.click()
    assert wf.current() == "hydro" and wf.title.text() == "Hydrocarbon identifier"
    # A finished sub-menu offers the next numbered step, never itself or its sibling.
    wf.mark_done("hydro", "done")
    assert wf.next_key() == "quant"


def test_left_tabs_and_reopening_closed_panels(qtbot, window):
    w = window
    w.show()
    tabs = [w.left_tabs.tabText(i) for i in range(w.left_tabs.count())]
    # The bottom bar became tabs here, so results, standards, history and log live in the workspace.
    assert tabs[0] == "Runs" and tabs[2] == "Layers"
    assert tabs[3:] == ["Results", "Standards", "History", "Log"]
    # Opening a step brings the Tool tab forward and names it after the step.
    w.left_tabs.setCurrentIndex(0)
    w.step_actions["background"].trigger()
    assert w.left_tabs.currentIndex() == w.tool_tab
    assert w.left_tabs.tabText(w.tool_tab) == "Background"      # icons replace the emoji on tabs
    # Closed panels can be brought back from the toolbar's own toggle actions.
    for dock in (w.workspace_dock, w.spectrum_dock):
        dock.close()
        assert not dock.isVisible()
        dock.toggleViewAction().trigger()
        assert dock.isVisible()
    w.workspace_dock.close()
    w.reset_layout()
    assert w.workspace_dock.isVisible()


def test_every_step_says_what_and_why(qtbot, window):
    wf = window.workflow
    for key in wf.keys():
        info = wf.step_info(key)
        assert info["blurb"] and info["why"], key
        if not info["parent"]:          # sub-menus have no toolbar button
            assert info["why"] in window.step_actions[key].toolTip()
    assert all(s["why"] in wf.guide_html() for s in wf._steps)


def test_workflow_step_opens_settings_and_go_runs_it(qtbot, window):
    w = window
    wf = w.workflow
    assert wf.current() == "load"
    # A toolbar step button only opens the settings; nothing runs.
    n = w.undo.count()
    w.step_actions["background"].trigger()
    assert wf.current() == "background"
    assert w.undo.count() == n
    assert w.bg_panel.isVisibleTo(wf)
    assert not wf.is_done("background")
    # Go runs it, ticks the step and offers the next one.
    wf.go_btn.click()
    qtbot.waitUntil(lambda: w.undo.count() == n + 1, timeout=TIMEOUT)
    assert wf.is_done("background")
    assert wf.banner.isVisibleTo(wf)
    assert "Peaks" in wf.banner_next.text()
    assert w.step_actions["background"].text().startswith("✓")
    wf.banner_next.click()
    assert wf.current() == "blobs"
    assert w.step_actions["blobs"].isChecked()
    # Removing the layer un-ticks the step.
    w.remove_background_all()
    assert not wf.is_done("background")


def test_hydrocarbons_saved_and_reloaded(qtbot, hydro_window, monkeypatch, tmp_path):
    from pygcgc.gui.mainwindow import MainWindow
    w = hydro_window
    i_std = _find_std_series(w, qtbot, monkeypatch)
    std_series = w.hydrocarbons.series_for(i_std)
    box = [std_series.to_absolute(15, 0.15), std_series.to_absolute(25, 0.15),
          std_series.to_absolute(25, 0.45), std_series.to_absolute(15, 0.45)]
    w.hydro_panel.draw_btn.setChecked(True)
    w.view2d.stencil_roi.setPoints([[t1 / 60.0, t2] for t1, t2 in box])
    w.view2d.stencil_roi.setPos([0, 0])
    w.hydro_panel.saveStencilRequested.emit("branched C15-25", 3.0)

    path = tmp_path / "hydro.pygcgc.json"
    w.project.path = path
    w.save_project()

    w2 = MainWindow()
    qtbot.addWidget(w2)
    done = []
    w2.open_project(path, on_loaded=lambda: done.append(True))
    qtbot.waitUntil(lambda: bool(done), timeout=TIMEOUT)
    assert [s.name for s in w2.hydrocarbons.stencils.stencils] == ["branched C15-25"]
    assert w2.hydrocarbons.stencils.stencils[0].response_factor == 3.0
    saved = w2.hydrocarbons.series_for(w2.project.index_of("STD"))          # the ladder is kept with the project
    assert saved is not None and sorted(saved.hits) == sorted(std_series.hits)
    assert saved.hits[15].t1 == pytest.approx(std_series.hits[15].t1)
    w2.jobs.wait(5000)


# ------------------------------------------------------------------ 🗂️ batch
def test_batch_dialog_runs_and_fills_tabs(qtbot, window):
    w = window
    w.show_batch_dialog()
    dlg = w.batch_dialog
    assert dlg.isVisible()
    folder = w.project.entries[0].path.parent.parent     # the synthetic_folder root
    dlg.folder_edit.setText(str(folder))
    dlg.workers.setValue(1)
    dlg._run_clicked()
    qtbot.waitUntil(lambda: dlg.result is not None, timeout=TIMEOUT)
    assert dlg.result.errors == {}
    assert set(dlg.result.prepped) == {"RUN_A", "RUN_B", "RUN_C_BLANK"}
    assert dlg.feature_table.rowCount() > 0
    assert dlg.result.report is None       # no calibration_levels.csv in this fixture: quant is skipped
    # PCA (3 runs) ran automatically after the batch finished and drew something.
    assert len(dlg.pca_plot.listDataItems()) > 0
    dlg.group_a_list.setCurrentRow(0)
    dlg.group_a_list.item(0).setSelected(True)
    dlg.group_a_list.item(1).setSelected(True)
    dlg.group_b_list.item(2).setSelected(True)
    dlg._draw_volcano()      # too few in group B (only 1): should not raise, just do nothing


def test_batch_navigate_switches_run_when_open_else_says_so(qtbot, window):
    w = window
    w.show_batch_dialog()
    dlg = w.batch_dialog
    other = next(n for n in ("RUN_A", "RUN_B", "RUN_C_BLANK") if n != w.project.entries[w.index].name)

    dlg.navigateRequested.emit(other, 100.0, 1.0)
    assert w.project.entries[w.index].name == other

    dlg.navigateRequested.emit("NOT_A_LOADED_RUN", 0.0, 0.0)
    assert "NOT_A_LOADED_RUN" in w.statusBar().currentMessage()


def test_batch_config_persists_settings_choices():
    from pygcgc import BatchConfig, BatchSettings
    cfg = BatchConfig(settings=BatchSettings(run_siloxane=True, run_blobs=False))
    assert cfg.settings.run_siloxane and not cfg.settings.run_blobs
    assert BatchConfig.from_dict(cfg.to_dict()) == cfg


# ------------------------------------------------------------------ 📤 report export
def test_batch_dialog_exports_xlsx_html_pdf(qtbot, window, cal_folder, tmp_path):
    w = window
    w.show_batch_dialog()
    dlg = w.batch_dialog
    dlg.folder_edit.setText(str(cal_folder))
    dlg.workers.setValue(1)
    assert not dlg.export_btn.isEnabled()
    dlg._run_clicked()
    qtbot.waitUntil(lambda: dlg.result is not None, timeout=TIMEOUT)
    assert dlg.result.report is not None            # cal_folder has calibration_levels.csv
    assert dlg.export_btn.isEnabled()

    dlg._do_export(tmp_path / "report.xlsx", "xlsx")
    assert (tmp_path / "report.xlsx").exists()

    dlg.group_a_list.item(0).setSelected(True)
    dlg.group_a_list.item(1).setSelected(True)
    dlg.group_b_list.item(2).setSelected(True)
    dlg.group_b_list.item(3).setSelected(True)
    dlg._do_export(tmp_path / "report.html", "html")
    html = (tmp_path / "report.html").read_text(encoding="utf-8")
    assert "Volcano" in html and "PCA" in html

    dlg._do_export(tmp_path / "report.pdf", "pdf")
    assert (tmp_path / "report.pdf").read_bytes().startswith(b"%PDF")
    assert "Wrote" in dlg.status.text()


def test_batch_dialog_export_disabled_until_a_run_completes(window):
    w = window
    w.show_batch_dialog()
    dlg = w.batch_dialog
    assert not dlg.export_btn.isEnabled()


# ------------------------------------------------------------------ ✂️ overlap
def test_overlap_dialog_lists_blobs_and_fills_the_unique_ion_map(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.show_overlap_dialog()
    dlg = w.overlap_dialog
    assert dlg.isVisible()
    blobs = w.blobs.blobs_for(i)
    assert dlg.blob_combo.count() == len(blobs)
    assert dlg.cube is not None
    assert dlg.apex_table.rowCount() > 0

    # GUI_PEAKS[0] (t1=40, ion 50): pick it out by t1 rather than assuming
    # index 0. Its own apex must be ion 50 -- m/z 54's deliberate rising
    # baseline artefact (see conftest.py) can bleed a spurious second
    # cluster into any blob's bounding box, so that count isn't asserted.
    target = int(np.argmin(np.abs(blobs.t1 - 40.0)))
    combo_row = next(r for r in range(dlg.blob_combo.count()) if dlg.blob_combo.itemData(r) == target)
    dlg.blob_combo.setCurrentIndex(combo_row)
    assert any(a.ion == 50.0 for a in dlg._apex_result.apexes)


def test_overlap_dialog_fit_and_mcr_run_as_background_jobs(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.show_overlap_dialog()
    dlg = w.overlap_dialog

    dlg.auto_k.setChecked(False)
    dlg.k_spin.setValue(1)
    dlg._run_fit()
    qtbot.waitUntil(lambda: dlg.fit_table.rowCount() == 1, timeout=TIMEOUT)
    assert "BIC" in dlg.fit_status.text()
    assert dlg.fit_table.columnCount() == 7
    assert dlg.fit_table.item(0, 6).text() in ("provisional", "low")
    assert dlg.component_view._items

    dlg.mcr_k.setValue(1)
    dlg._run_mcr()
    qtbot.waitUntil(lambda: dlg.mcr_table.rowCount() == 1, timeout=TIMEOUT)
    assert "residual" in dlg.mcr_status.text()


def test_overlap_dialog_switches_blob_and_refreshes_on_reshow(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.show_overlap_dialog()
    dlg = w.overlap_dialog
    first_cube = dlg.cube
    if dlg.blob_combo.count() > 1:
        dlg.blob_combo.setCurrentIndex(1)
        assert dlg.cube is not first_cube

    j = w.project.index_of("RUN_B")
    w.set_run(j)
    w.show_overlap_dialog()
    assert dlg._name == "RUN_B"


def test_overlap_dialog_can_use_border_free_candidate_regions(qtbot, window):
    w = window
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.show_overlap_dialog()
    dlg = w.overlap_dialog
    dlg._find_candidate_regions()
    qtbot.waitUntil(lambda: dlg.find_regions_btn.isEnabled(), timeout=TIMEOUT)
    assert dlg._source == "candidates"
    assert dlg.blob_combo.count() > 0
    assert dlg.cube is not None
    assert dlg.cube.mask.all()


def _write_test_msp(path):
    from pygcgc.core.identify import LibraryEntry, write_msp
    entries = [LibraryEntry(name="Compound50", mz=np.array([50]), intensity=np.array([999.0])),
              LibraryEntry(name="Compound51", mz=np.array([51]), intensity=np.array([999.0]))]
    return write_msp(entries, path)


def test_identify_dialog_loads_library_and_matches_a_blob(qtbot, window, tmp_path):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.show_identify_dialog()
    dlg = w.identify_dialog
    assert dlg.isVisible()
    assert not dlg.identify_btn.isEnabled()          # no library loaded yet

    lib_path = _write_test_msp(tmp_path / "lib.msp")
    dlg.load_library_file(lib_path)
    assert dlg.library is not None
    assert "2 compound" in dlg.library_label.text()
    assert dlg.identify_btn.isEnabled()

    dlg._run_identify()
    qtbot.waitUntil(lambda: len(dlg.results) > 0, timeout=TIMEOUT)
    assert dlg.table.rowCount() == len(dlg.results)

    blobs = w.blobs.blobs_for(i)
    target_blob = int(np.argmin(np.abs(blobs.t1 - 40.0)))     # GUI_PEAKS[0]: ion 50
    ident = next(r for r in dlg.results if r.blob_index == target_blob)
    assert ident.best is not None and ident.best.name == "Compound50"


def test_library_browser_filters_edits_deletes_and_saves(qtbot, window, tmp_path):
    from pygcgc.core.identify import load_msp
    from pygcgc.gui.identify_dialog import LibraryBrowser
    w = window
    w.show_identify_dialog()
    dlg = w.identify_dialog
    dlg.load_library_file(_write_test_msp(tmp_path / "lib.msp"))
    n_all = len(dlg.library.entries)
    browser = LibraryBrowser(dlg.library)
    qtbot.addWidget(browser)
    assert browser.table.rowCount() == 2                         # decoys hidden by default
    browser.show_decoys.setChecked(True)
    assert browser.table.rowCount() == n_all
    browser.show_decoys.setChecked(False)
    browser.filter.setText("51")
    assert browser.table.rowCount() == 1
    browser.filter.setText("")
    # Rename and set an RI in place; the rebuilt library reflects both.
    browser.table.item(0, 0).setText("Renamed50")
    browser.table.item(0, 1).setText("1234")
    assert browser.library.entries[browser._rows[0]].name == "Renamed50"
    assert browser.library.entries[browser._rows[0]].ri == 1234.0
    browser.table.selectRow(1)
    browser._delete()
    assert browser.table.rowCount() == 1
    out = tmp_path / "edited.msp"
    browser.save_to(out)
    saved = load_msp(out)
    assert [e.name for e in saved] == ["Renamed50"]
    assert saved[0].ri == 1234.0


def test_identify_settings_reach_the_search(qtbot, window, tmp_path):
    from dataclasses import replace
    w = window
    w.show_identify_dialog()
    dlg = w.identify_dialog
    assert dlg.settings.min_score == 0.2
    dlg.settings = replace(dlg.settings, min_score=0.99)
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    w.set_run(w.project.index_of("RUN_A"))
    w.show_identify_dialog()
    dlg.load_library_file(_write_test_msp(tmp_path / "lib.msp"))
    dlg._run_identify()
    qtbot.waitUntil(lambda: len(dlg.results) > 0, timeout=TIMEOUT)
    assert dlg.settings.min_score == 0.99


def test_identify_dialog_promote_and_inhouse_write_files(qtbot, window, tmp_path):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.show_identify_dialog()
    dlg = w.identify_dialog
    dlg.load_library_file(_write_test_msp(tmp_path / "lib.msp"))
    dlg._run_identify()
    qtbot.waitUntil(lambda: len(dlg.results) > 0, timeout=TIMEOUT)

    blobs = w.blobs.blobs_for(i)
    target_blob = int(np.argmin(np.abs(blobs.t1 - 40.0)))
    ident = next(r for r in dlg.results if r.blob_index == target_blob)
    assert ident.best is not None

    target_csv = tmp_path / "target_table.csv"
    dlg.promote_to_file(ident, target_csv)
    from pygcgc.core.targets import load_targets
    loaded = load_targets(target_csv)
    assert loaded[0].compound == "Compound50" and loaded[0].quant_ion == 50
    before = target_csv.read_bytes()
    dlg.promote_to_file(ident, target_csv)
    assert "Could not write target table" in dlg.status.text()
    assert target_csv.read_bytes() == before

    inhouse_path = tmp_path / "inhouse.msp"
    dlg.add_to_inhouse_file(ident, "MyStandard", inhouse_path)
    from pygcgc.core.identify import load_msp
    entries = load_msp(inhouse_path)
    assert entries and entries[0].name == "MyStandard"
    assert 50 in entries[0].mz.tolist()


def test_identify_dialog_switches_run_on_reshow(qtbot, window):
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    w.show_identify_dialog()
    dlg = w.identify_dialog
    first_name = dlg._name

    j = w.project.index_of("RUN_B") if dlg._name != "RUN_B" else w.project.index_of("RUN_C_BLANK")
    w.set_run(j)
    w.show_identify_dialog()
    assert dlg._name != first_name


def test_surface3d_view_handles_empty_data_without_raising(qtbot):
    from pygcgc.gui.views.surface3d import Surface3DView
    view = Surface3DView()
    qtbot.addWidget(view)
    view.set_data(np.zeros((0, 0)))
    view.set_data(np.random.default_rng(0).random((6, 7)))
    assert view.available


def test_runs_table_sets_type_level_set_and_volume_and_saves_them(qtbot, window, tmp_path):
    from pygcgc.gui.mainwindow import MainWindow
    w = window
    i_b = w.project.index_of("RUN_B")
    # Typed into the table: RUN_B becomes a second calibration level.
    table = w.runs_panel.table
    table.item(i_b, 1).setText("calibration")
    table.item(i_b, 2).setText("L2")
    assert w.project.entries[i_b].kind == "calibration" and w.project.entries[i_b].level == "L2"
    table.item(i_b, 4).setText("not a number")                    # refused, the old value comes back
    assert w.project.entries[i_b].volume is None and table.item(i_b, 4).text() == ""
    # Several runs at once, as the right-click menu does.
    w.runs_panel.infoEdited.emit([0, i_b], "cal_set", "mon")
    w.runs_panel.infoEdited.emit([i_b], "volume", "1.5")
    info = {r.name: r for r in w.quant.info()}
    assert info["RUN_B"].is_calibration and info["RUN_B"].cal_set == "mon" and info["RUN_B"].volume == 1.5
    # Written out as a sequence.csv that loads back the same way...
    w.save_sequence_dialog(tmp_path / "sequence.csv")
    text = (tmp_path / "sequence.csv").read_text()
    assert "RUN_B.cdf,calibration,L2,mon,1.5" in text
    # ...and kept in the project file, winning over any sequence.csv on reopening.
    path = tmp_path / "proj.pygcgc.json"
    w.project.save(path)
    w2 = MainWindow()
    qtbot.addWidget(w2)
    done = []
    w2.open_project(path, on_loaded=lambda: done.append(True))
    qtbot.waitUntil(lambda: bool(done), timeout=TIMEOUT)
    e = w2.project.entries[w2.project.index_of("RUN_B")]
    assert (e.kind, e.level, e.cal_set, e.volume) == ("calibration", "L2", "mon", 1.5)
    w2.jobs.wait(5000)

def test_reopening_a_project_reuses_cached_background_and_blobs(qtbot, window, tmp_path, monkeypatch):
    from pygcgc.gui import blob_state, project as project_module
    from pygcgc.gui.mainwindow import MainWindow
    w = window
    w.apply_background()
    qtbot.waitUntil(lambda: all(e.background is not None for e in w.project.entries), timeout=TIMEOUT)
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    n_blobs = {e.name: len(w.blobs.by_run[e]) for e in w.project.entries}
    assert all(getattr(e.background, "cache_key", None) for e in w.project.entries)
    path = tmp_path / "proj.pygcgc.json"
    w.project.save(path)
    w.jobs.wait(5000)

    # Nothing may be recomputed on reopening: both steps now fail if called.
    def boom(*a, **k):
        raise AssertionError("recomputed instead of read from the cache")
    monkeypatch.setattr(project_module, "remove_background", boom)
    monkeypatch.setattr(blob_state, "find_blobs", boom)
    w2 = MainWindow()
    w2.blobs.set_method(blob_state.LEGACY)
    qtbot.addWidget(w2)
    done = []
    w2.open_project(path, on_loaded=lambda: done.append(True))
    qtbot.waitUntil(lambda: bool(done) and all(e.background is not None for e in w2.project.entries), timeout=TIMEOUT)
    w2.run_blobs()
    qtbot.waitUntil(lambda: not w2.blobs.running and w2.blobs.was_run, timeout=TIMEOUT)
    assert {e.name: len(w2.blobs.by_run[e]) for e in w2.project.entries} == n_blobs
    w2.jobs.wait(5000)
