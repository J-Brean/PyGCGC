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
    np.testing.assert_allclose(x, dialog._x)
    np.testing.assert_allclose(total, dialog.image.image.sum(axis=0))
    contributions = sum(curve.getData()[1] for curve in dialog.component_curves)
    np.testing.assert_allclose(contributions + .025*len(dialog._y), total)
    # Each 2D contour is an equal-height line on the same fitted Gaussian.
    for p, contours in zip(dialog._truth, dialog.contours):
        cx, cy, sx, sy, rho, _height = p
        for curve, level in zip(contours, (.10, .50)):
            x, y = curve.getData()
            dx, dy = (x-cx)/sx, (y-cy)/sy
            np.testing.assert_allclose(np.exp(-.5*(dx*dx+dy*dy-2*rho*dx*dy)/(1-rho*rho)), level)


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
    assert dialog.trace.backgroundBrush().color().name() == (DARK if dark else LIGHT)["base"].lower()
    assert qapp.styleSheet() == stylesheet
    apply_theme(host, not dark)
    qtbot.waitUntil(lambda: dialog.trace.backgroundBrush().color().name()
                   == (LIGHT if dark else DARK)["base"].lower())
    dialog.close()
