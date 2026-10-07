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
    qss = t._QSS.format(font=t.FONT, **colours)
    rule = next(line for line in qss.splitlines() if line.startswith("QMainWindow {"))
    assert "transparent" not in rule
    assert colours["window"] in rule
