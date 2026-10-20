"""Every top-level window has a background: only the main window paints its own."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")


@pytest.mark.parametrize("theme", ["LIGHT", "DARK"])
def test_main_windows_get_an_opaque_background(theme):
    # A transparent QMainWindow rule left the 📈 Calibrate window black on screen, with
    # dark text on it; MainWindow paints its gradient in paintEvent and needs no rule.
    from pygcgc.gui import theme as t
    colours = getattr(t, theme)
    qss = t._QSS.format(font=t.FONT, **colours, **t.DENSITY["standard"])
    rule = next(line for line in qss.splitlines() if line.startswith("QMainWindow {"))
    assert "transparent" not in rule
    assert colours["window"] in rule


def test_opening_another_window_does_not_replace_the_application_theme(qtbot, qapp):
    from PySide6 import QtWidgets
    from pygcgc.gui.theme import apply_theme

    def window():
        widget = QtWidgets.QWidget()
        widget.view2d = widget.view1d = widget.spectrum = None
        qtbot.addWidget(widget)
        return widget

    first, second = window(), window()
    apply_theme(first, False)
    style, font, palette, stylesheet = qapp.style(), qapp.font(), qapp.palette(), qapp.styleSheet()
    apply_theme(second, False)
    assert qapp.style() is style
    assert qapp.font() == font and qapp.palette() == palette
    assert qapp.styleSheet() == stylesheet
    apply_theme(second, True)
    assert qapp.styleSheet() != stylesheet
    apply_theme(first, False)
    assert qapp.styleSheet() == stylesheet


def test_colour_theme_switch_does_not_replace_the_shared_fusion_style(qtbot, qapp, monkeypatch):
    from PySide6 import QtWidgets
    from pygcgc.gui.theme import apply_theme
    widget = QtWidgets.QWidget()
    widget.view2d = widget.view1d = widget.spectrum = None
    qtbot.addWidget(widget)
    apply_theme(widget, False)
    replaced = []
    original = qapp.setStyle
    monkeypatch.setattr(qapp, 'setStyle', lambda style: (replaced.append(style), original(style))[1])
    apply_theme(widget, True)
    apply_theme(widget, False)
    assert replaced == []
