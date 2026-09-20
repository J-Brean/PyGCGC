"""Light and dark themes (Qt Fusion style with a matching palette).
Plot backgrounds follow the theme; the chromatogram colour map does not.

The light theme is the default and borrows the Financial Times palette: a salmon
paper background, near-black text and claret highlights."""

from __future__ import annotations

import pyqtgraph as pg
from PySide6 import QtGui, QtWidgets

DARK = {
    "window": "#16121F", "base": "#1D1829", "alt": "#231D31", "text": "#ECE8F3",
    "mid": "#8E86A0", "button": "#241E33", "highlight": "#7A3D8C", "plot": "#0E0B15",
}
LIGHT = {
    "window": "#FFF1E5", "base": "#FFF9F5", "alt": "#F7E7D8", "text": "#33302E",
    "mid": "#807973", "button": "#F2DFCE", "highlight": "#990F3D", "plot": "#FFF9F5",
}


def apply_theme(window, dark: bool) -> None:
    c = DARK if dark else LIGHT
    app = QtWidgets.QApplication.instance()
    app.setStyle("Fusion")
    pal = QtGui.QPalette()
    R = QtGui.QPalette.ColorRole
    col = QtGui.QColor
    pal.setColor(R.Window, col(c["window"]))
    pal.setColor(R.WindowText, col(c["text"]))
    pal.setColor(R.Base, col(c["base"]))
    pal.setColor(R.AlternateBase, col(c["alt"]))
    pal.setColor(R.Text, col(c["text"]))
    pal.setColor(R.Button, col(c["button"]))
    pal.setColor(R.ButtonText, col(c["text"]))
    pal.setColor(R.ToolTipBase, col(c["base"]))
    pal.setColor(R.ToolTipText, col(c["text"]))
    pal.setColor(R.PlaceholderText, col(c["mid"]))
    pal.setColor(R.Mid, col(c["mid"]))
    pal.setColor(R.Highlight, col(c["highlight"]))
    pal.setColor(R.HighlightedText, col("#FFFFFF"))
    disabled = QtGui.QPalette.ColorGroup.Disabled
    for role in (R.WindowText, R.Text, R.ButtonText):
        pal.setColor(disabled, role, col(c["mid"]))
    app.setPalette(pal)
    bg = c["plot"]
    for view in (window.view2d, window.view1d, window.spectrum,
                 getattr(window, "calibration_panel", None), getattr(window, "batch_dialog", None),
                 getattr(window, "overlap_dialog", None)):
        if view is not None:
            view.set_background(bg)
    if getattr(window, "results", None) is not None:
        window.results.set_chart_background(bg)
    # Axis lines and labels must contrast with the plot background: pyqtgraph's
    # defaults are light grey, which vanishes on paper.
    fg = col(c["text"])
    for gv in window.findChildren(pg.GraphicsView):
        for item in (gv.scene().items() if gv.scene() is not None else []):
            if isinstance(item, pg.AxisItem):
                item.setPen(fg)
                item.setTextPen(fg)
