"""The explanatory animation must stay local, consistent and stoppable."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from pygcgc.gui.fitting_help import FittingHelpDialog


def test_play_pause_seek_replay_and_close_stop_the_timer(qtbot):
    dialog = FittingHelpDialog()
    qtbot.addWidget(dialog)
    assert not dialog.timer.isActive()
    dialog.show()
    assert dialog.timer.isActive()
    qtbot.waitUntil(lambda: dialog.position_ms > 0)
    dialog.play_button.click()
    paused = dialog.position_ms
    qtbot.wait(90)
    assert dialog.position_ms == paused and not dialog.timer.isActive()
    dialog.slider.setValue(9000)
    assert dialog.heading.text().startswith("3.")
    dialog.replay_button.click()
    assert dialog.position_ms == 0 and dialog.timer.isActive()
    dialog.set_position(dialog.DURATION_MS)
    assert not dialog.timer.isActive() and dialog.play_button.text() == "Play"
    dialog.play_button.click()
    assert dialog.position_ms == 0 and dialog.timer.isActive()
    dialog.close()
    assert not dialog.timer.isActive()
    dialog.show()
    assert dialog.position_ms == 0 and dialog.timer.isActive()
    dialog.hide()
    assert not dialog.timer.isActive()


def test_1d_and_2d_share_the_same_contributions_and_final_sum(qtbot):
    dialog = FittingHelpDialog()
    qtbot.addWidget(dialog)
    assert not dialog.sum_curve.isVisible()
    dialog.set_position(dialog.DURATION_MS)
    assert all(curve.isVisible() for curve in dialog.component_curves)
    x, total = dialog.sum_curve.getData()
    np.testing.assert_allclose(x, dialog._detector_time)
    # Folding the detector stream must recover the displayed map, scan for scan.
    np.testing.assert_allclose(total.reshape(len(dialog._x), len(dialog._y)).T,
                               dialog.image.image)
    contributions = sum(curve.getData()[1] for curve in dialog.component_curves)
    np.testing.assert_allclose(contributions + .025, total)
    # Each 2D contour is an equal-height line on the same fitted Gaussian.
    for p, contours in zip(dialog._truth, dialog.contours):
        cx, cy, sx, sy, rho, _height = p
        for curve, level in zip(contours, (.10, .50)):
            x, y = curve.getData()
            dx, dy = (x-cx)/sx, (y-cy)/sy
            np.testing.assert_allclose(np.exp(-.5*(dx*dx+dy*dy-2*rho*dx*dy)/(1-rho*rho)), level)


def test_detector_trace_has_250_ms_pulses_at_modulation_intervals(qtbot):
    """A broad first-dimension envelope must never stand in for the raw trace."""
    from scipy.signal import find_peaks, peak_widths
    dialog = FittingHelpDialog()
    qtbot.addWidget(dialog)
    dialog.set_position(dialog.DURATION_MS)
    for curve in dialog.component_curves:
        time, signal = curve.getData()
        peaks, _ = find_peaks(signal, height=.03)
        assert len(peaks) >= 5
        widths = peak_widths(signal, peaks, rel_height=.5)[0] * dialog._dt
        np.testing.assert_allclose(widths, .250, atol=.001)
        # Tilt moves the pulse slightly within each 10 s cycle.
        np.testing.assert_allclose(np.diff(time[peaks]), 10., atol=.06)
        assert np.all(np.diff(time) > 0)
    assert "detector time" in dialog.trace.getAxis("bottom").labelText.lower()


@pytest.mark.parametrize("position", [0, 4500, 9000, 16000])
def test_highlight_and_enlargement_follow_the_same_modulation(qtbot, position):
    dialog = FittingHelpDialog()
    qtbot.addWidget(dialog)
    dialog.set_position(position)
    selected = dialog._selected_cycle
    start = dialog._starts[selected]
    for band in dialog.cycle_bands:
        np.testing.assert_allclose(band.getRegion(), [start, start+10])
    np.testing.assert_allclose(dialog.zoom_observed.getData()[1],
                               dialog.image.image[:, selected])
    full_sum = dialog.sum_curve.getData()[1].reshape(len(dialog._x), len(dialog._y))
    np.testing.assert_allclose(dialog.zoom_sum.getData()[1], full_sum[selected])
    for zoom, full in zip(dialog.zoom_components, dialog.component_curves):
        np.testing.assert_allclose(zoom.getData()[1],
                                   full.getData()[1].reshape(full_sum.shape)[selected])


def test_peak_help_opens_without_finding_peaks_and_reuses_one_window(qtbot):
    from pygcgc.gui.panels.blobs import BlobsPanel
    panel = BlobsPanel()
    qtbot.addWidget(panel)
    requested = []
    panel.findRequested.connect(lambda: requested.append(True))
    panel.show()
    assert panel.fitting_help_btn.isVisible()
    panel.fitting_help_btn.click()
    help_window = panel._fitting_help
    assert help_window.isVisible() and help_window.timer.isActive()
    panel.fitting_help_btn.click()
    assert panel._fitting_help is help_window and not requested
    help_window.close()
    assert not help_window.timer.isActive()
    panel.set_method("legacy")
    assert panel.fitting_help_btn.isHidden()
    panel.set_method("ms-aware")
    assert panel.fitting_help_btn.isVisible()


def test_manual_refit_starts_with_tilt_and_offers_the_same_help(qtbot):
    from pygcgc.gui.overlap_dialog import OverlapDialog
    dialog = OverlapDialog(None)
    qtbot.addWidget(dialog)
    assert dialog.model_combo.currentText() == "tilt"
    dialog.fitting_help_btn.click()
    assert dialog._fitting_help.isVisible()
    dialog._fitting_help.close()


@pytest.mark.parametrize("dark", [False, True])
def test_help_uses_the_current_theme_without_replacing_it(qtbot, qapp, dark):
    from pygcgc.gui.theme import DARK, LIGHT, apply_theme
    from PySide6 import QtWidgets
    host = QtWidgets.QWidget()
    host.view2d = host.view1d = host.spectrum = None
    qtbot.addWidget(host)
    apply_theme(host, dark)
    stylesheet = qapp.styleSheet()
    dialog = FittingHelpDialog(host)
    qtbot.addWidget(dialog)
    dialog.show()
    plots = (dialog.trace, dialog.zoom, dialog.map)
    assert all(plot.backgroundBrush().color().name() == (DARK if dark else LIGHT)["base"].lower()
               for plot in plots)
    assert qapp.styleSheet() == stylesheet
    apply_theme(host, not dark)
    qtbot.waitUntil(lambda: all(plot.backgroundBrush().color().name()
                               == (LIGHT if dark else DARK)["base"].lower() for plot in plots))
    dialog.close()
