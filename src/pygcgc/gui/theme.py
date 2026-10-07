"""Light and dark themes (Qt Fusion style, a matching palette and a flat stylesheet).
Plot backgrounds follow the theme; the chromatogram colour map does not.

Swiss-style: Helvetica/Arial, flat surfaces with no frames, generous spacing, one
green accent. The dark theme is near-black with a Spotify-like feel; the light theme
is white paper with the same layout. `FadeTabWidget` and `fade_in` give the gentle
transition when a tab or workflow step changes."""

from __future__ import annotations

import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

# Depth comes from layers: a recessed window, raised "card" panels (base) with a hairline
# border, and inset fields (alt). `tint` is the accent at low opacity, for active states.
DARK = {
    "window": "#0A0A0A", "base": "#161616", "alt": "#1F1F1F", "text": "#FFFFFF",
    "mid": "#A1A1A1", "button": "#242424", "hover": "#2E2E2E", "highlight": "#1DB954",
    "plot": "#161616", "line": "#303030", "onhl": "#000000", "gohover": "#1ED760",
    "tint": "rgba(29, 185, 84, 38)", "card": "rgba(22, 22, 22, 135)", "inset": "rgba(22, 22, 22, 70)", "stripe": "rgba(255, 255, 255, 10)", "plot_alpha": 135, "glow": "#10301F",
}
LIGHT = {
    "window": "#F4F5F7", "base": "#FFFFFF", "alt": "#F7F7F8", "text": "#111318",
    "mid": "#6B7080", "button": "#FFFFFF", "hover": "#EEF0F3", "highlight": "#4F46E5",
    "plot": "#FFFFFF", "line": "#E2E4E9", "onhl": "#FFFFFF", "gohover": "#4338CA",
    "tint": "rgba(79, 70, 229, 26)", "card": "rgba(255, 255, 255, 120)", "inset": "rgba(255, 255, 255, 70)", "stripe": "rgba(17, 19, 24, 9)", "plot_alpha": 120, "glow": "#E6E5FB",
}

FONT = '"Helvetica Neue", Helvetica, Arial, sans-serif'

_QSS = """
* {{ font-family: {font}; }}
QWidget {{ color: {text}; }}
QWidget:disabled {{ color: {mid}; }}
QMainWindow {{ background: {window}; }}  /* the main window paints its own gradient instead: see MainWindow.paintEvent */
QDialog {{ background: {window}; }}
QDockWidget {{ background: transparent; }}
QMessageBox {{ background: {window}; }}
QLabel, QCheckBox, QRadioButton {{ background: transparent; }}
QAbstractScrollArea {{ background: transparent; border: none; }}
QScrollArea, QScrollArea > QWidget, QScrollArea > QWidget > QWidget {{ background: transparent; border: none; }}
QCheckBox, QRadioButton {{ spacing: 8px; }}
QToolTip {{ background: {button}; color: {text}; border: none; padding: 4px 6px; }}

QWidget#card {{ background: {card}; border: 1px solid {line}; border-radius: 10px; }}
QTabWidget::pane {{ border: none; border-top: 1px solid {line}; background: transparent; }}
QTabWidget#main_tabs::pane {{ border: 1px solid {line}; border-radius: 10px; background: {card}; }}
QTabBar {{ background: transparent; }}
QMainWindow > QTabBar {{ margin-left: 12px; }}
QTabBar::tab:first {{ padding-left: 0; margin-left: 0; }}
QTabBar::tab {{ background: transparent; color: {mid}; padding: 8px 14px 7px 14px; font-weight: 500;
    margin-right: 2px; border: none; border-bottom: 2px solid transparent; }}
QTabBar::tab:hover {{ color: {text}; }}
QTabBar::tab:selected {{ color: {text}; border-bottom: 2px solid {highlight}; }}

QDockWidget {{ color: {mid}; }}
QDockWidget::title {{ background: transparent; padding: 8px 12px; text-align: left; font-weight: 600; }}

QPushButton, QToolButton {{ background: {button}; color: {text}; border: 1px solid {line};
    border-radius: 8px; padding: 6px 16px; min-height: 18px; }}
QPushButton:hover, QToolButton:hover {{ background: {hover}; }}
QPushButton:pressed, QToolButton:pressed {{ background: {line}; }}
QPushButton:disabled, QToolButton:disabled {{ color: {mid}; background: {alt}; }}
QPushButton:default, QPushButton:checked {{ background: {highlight}; color: {onhl}; border: 1px solid {highlight}; }}
QPushButton#go {{ background: {highlight}; color: {onhl}; border: none; border-radius: 16px; min-height: 20px; font-weight: bold; }}
QPushButton#go:hover {{ background: {gohover}; }}
QPushButton#go:disabled {{ background: {button}; color: {mid}; }}
QToolBar {{ background: transparent; border: none; spacing: 2px; padding: 4px 4px; }}
QFrame#banner {{ background: {tint}; border: 1px solid {tint}; border-radius: 10px; }}
QToolBar QLabel {{ color: {mid}; padding: 0 8px; }}
QToolBar::separator {{ background: {line}; width: 1px; margin: 4px 6px; }}
QToolButton {{ border: 1px solid transparent; border-radius: 8px; padding: 6px 10px; margin: 0 1px; background: transparent; min-height: 18px; }}
QToolButton:hover {{ background: {hover}; }}
QToolButton:checked {{ background: {tint}; color: {highlight}; border: 1px solid {tint}; }}

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QPlainTextEdit, QTextEdit {{
    background: {alt}; color: {text}; border: 1px solid {line}; border-radius: 6px;
    padding: 5px 8px; min-height: 18px; selection-background-color: {highlight}; selection-color: {onhl}; }}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{ border: 1px solid {highlight}; }}
QComboBox::drop-down {{ border: none; width: 18px; }}
QComboBox QAbstractItemView {{ background: {base}; selection-background-color: {highlight};
    selection-color: {onhl}; border: 1px solid {line}; outline: none; }}

QTableView, QTreeView, QListView, QTableWidget, QListWidget, QTreeWidget {{
    background: {inset}; alternate-background-color: {stripe}; border: 1px solid {line}; border-radius: 6px;
    gridline-color: transparent;
    selection-background-color: {highlight}; selection-color: {onhl}; outline: none; }}
QHeaderView::section {{ background: {window}; color: {mid}; border: none;
    border-bottom: 1px solid {line}; padding: 5px 6px; }}

QGroupBox {{ border: none; border-top: 1px solid {line}; margin-top: 14px; padding-top: 8px; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 0px; color: {text}; font-weight: 600; }}

QMenuBar {{ background: transparent; color: {text}; padding: 2px 4px; }}
QMenuBar::item {{ background: transparent; color: {text}; padding: 4px 8px; }}
QMenuBar::item:selected {{ background: {button}; color: {text}; border-radius: 4px; }}
QMenu {{ background: {base}; color: {text}; border: 1px solid {line}; padding: 4px; }}
QMenu::item {{ color: {text}; padding: 5px 22px; border-radius: 4px; }}
QMenu::item:disabled {{ color: {mid}; }}
QMenu::section {{ color: {mid}; padding: 6px 10px 2px 10px; }}
QMenu::item:selected {{ background: {highlight}; color: {onhl}; }}
QMenu::separator {{ height: 1px; background: {line}; margin: 4px 6px; }}

QStatusBar {{ background: transparent; color: {mid}; }}
QStatusBar::item {{ border: none; }}
QStatusBar QLabel {{ padding: 2px 12px; }}

QScrollBar:vertical {{ background: transparent; width: 10px; margin: 0; }}
QScrollBar:horizontal {{ background: transparent; height: 10px; margin: 0; }}
QScrollBar::handle {{ background: {hover}; border-radius: 5px; min-height: 24px; min-width: 24px; }}
QScrollBar::handle:hover {{ background: {mid}; }}
QScrollBar::add-line, QScrollBar::sub-line {{ width: 0; height: 0; }}
QScrollBar::add-page, QScrollBar::sub-page {{ background: transparent; }}

QProgressBar {{ background: {button}; border: none; border-radius: 3px; height: 6px;
    text-align: center; color: transparent; }}
QProgressBar::chunk {{ background: {highlight}; border-radius: 3px; }}
QSlider::groove:horizontal {{ background: {hover}; height: 4px; border-radius: 2px; }}
QSlider::sub-page:horizontal {{ background: {highlight}; border-radius: 2px; }}
QSlider::handle:horizontal {{ background: {text}; width: 12px; margin: -4px 0; border-radius: 6px; }}
QSplitter::handle {{ background: transparent; }}
QMainWindow::separator {{ background: transparent; width: 8px; height: 8px; }}
"""


def fade_in(widget: QtWidgets.QWidget, ms: int = 180) -> None:
    """Fade a widget from transparent to opaque. Skipped while it is not on screen,
    and the effect is removed afterwards so plots render normally again."""
    if widget is None or not widget.isVisible():
        return
    old = getattr(widget, "_fade_anim", None)
    if old is not None:
        old.stop()
    effect = QtWidgets.QGraphicsOpacityEffect(widget)
    widget.setGraphicsEffect(effect)
    anim = QtCore.QPropertyAnimation(effect, b"opacity", widget)
    anim.setDuration(ms)
    anim.setStartValue(0.0)
    anim.setEndValue(1.0)
    anim.setEasingCurve(QtCore.QEasingCurve.Type.OutCubic)
    anim.finished.connect(lambda: widget.setGraphicsEffect(None))
    widget._fade_anim = anim
    anim.start()


class FadeTabWidget(QtWidgets.QTabWidget):
    """QTabWidget whose newly shown page fades in."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.currentChanged.connect(lambda _i: fade_in(self.currentWidget()))


def apply_theme(window, dark: bool) -> None:
    c = DARK if dark else LIGHT
    app = QtWidgets.QApplication.instance()
    app.setStyle("Fusion")
    app.setFont(QtGui.QFont("Helvetica Neue" if "Helvetica Neue" in QtGui.QFontDatabase.families()
                            else "Arial", 9))
    app.setStyleSheet(_QSS.format(font=FONT, **c))
    # Qt's default layout spacing is cramped; give every panel room to breathe.
    for lay in window.findChildren(QtWidgets.QLayout):
        if isinstance(lay, (QtWidgets.QBoxLayout, QtWidgets.QGridLayout, QtWidgets.QFormLayout)):
            lay.setSpacing(max(lay.spacing(), 10))
    # Panels grey their hints with "palette(mid)", which Qt resolves to a washed-out tint
    # under a stylesheet. Swap in the theme's own grey, keeping the original for the next switch.
    for w in window.findChildren(QtWidgets.QWidget):
        original = w.property("_ss_original")
        if original is None:
            original = w.styleSheet()
            if "palette(mid)" not in original:
                continue
            w.setProperty("_ss_original", original)
        w.setStyleSheet(original.replace("palette(mid)", c["mid"] + "; font-size: 8pt"))
    if hasattr(window, "refresh_icons"):
        window.refresh_icons(c)
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
    pal.setColor(R.HighlightedText, col(c["onhl"]))
    disabled = QtGui.QPalette.ColorGroup.Disabled
    for role in (R.WindowText, R.Text, R.ButtonText):
        pal.setColor(disabled, role, col(c["mid"]))
    app.setPalette(pal)
    bg = QtGui.QColor(c["plot"])
    bg.setAlpha(c["plot_alpha"])
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
