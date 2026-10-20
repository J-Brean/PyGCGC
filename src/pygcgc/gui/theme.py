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

# How much room text and controls take. Compact is for laptop screens (1366 x 768, or 1080p
# scaled to 150 %); the default follows the screen unless the user chose one (View > Text size).
DENSITY = {
    "compact": dict(base_pt=8.0, hint_pt=7.0, spacing=4, margin=6, title_scale=1.3, blurb_pt=8.5,
                    tb_pad="2px 2px", btn_pad="3px 10px", btn_min="14px", tool_pad="3px 6px", tab_pad="5px 10px 4px 10px",
                    field_pad="2px 6px", field_min="14px", menu_pad="3px 16px", header_pad="3px 5px"),
    "standard": dict(base_pt=9.0, hint_pt=8.0, spacing=8, margin=10, title_scale=1.5, blurb_pt=9.5,
                     tb_pad="3px 4px", btn_pad="5px 14px", btn_min="18px", tool_pad="5px 8px", tab_pad="7px 12px 6px 12px",
                     field_pad="4px 8px", field_min="16px", menu_pad="4px 20px", header_pad="4px 6px"),
    "large": dict(base_pt=10.5, hint_pt=9.0, spacing=10, margin=12, title_scale=1.9, blurb_pt=10.5,
                  tb_pad="4px 4px", btn_pad="6px 16px", btn_min="18px", tool_pad="6px 10px", tab_pad="8px 14px 7px 14px",
                  field_pad="5px 8px", field_min="18px", menu_pad="5px 22px", header_pad="5px 6px"),
}

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
QTabBar::tab {{ background: transparent; color: {mid}; padding: {tab_pad}; font-weight: 500;
    margin-right: 2px; border: none; border-bottom: 2px solid transparent; }}
QTabBar::tab:hover {{ color: {text}; }}
QTabBar::tab:selected {{ color: {text}; border-bottom: 2px solid {highlight}; }}

QDockWidget {{ color: {mid}; }}
QDockWidget::title {{ background: transparent; padding: {tool_pad}; text-align: left; font-weight: 600; }}

QPushButton, QToolButton {{ background: {button}; color: {text}; border: 1px solid {line};
    border-radius: 8px; padding: {btn_pad}; min-height: {btn_min}; }}
QPushButton:hover, QToolButton:hover {{ background: {hover}; }}
QPushButton:pressed, QToolButton:pressed {{ background: {line}; }}
QPushButton:disabled, QToolButton:disabled {{ color: {mid}; background: {alt}; }}
QPushButton:default, QPushButton:checked {{ background: {highlight}; color: {onhl}; border: 1px solid {highlight}; }}
QPushButton#go {{ background: {highlight}; color: {onhl}; border: none; border-radius: 16px; min-height: 20px; font-weight: bold; }}
QPushButton#go:hover {{ background: {gohover}; }}
QPushButton#go:disabled {{ background: {button}; color: {mid}; }}
QToolBar {{ background: transparent; border: none; spacing: 2px; padding: {tb_pad}; }}
QFrame#banner {{ background: {tint}; border: 1px solid {tint}; border-radius: 10px; }}
QToolBar QLabel {{ color: {mid}; padding: 0 8px; }}
QToolBar::separator {{ background: {line}; width: 1px; margin: 4px 6px; }}
QToolButton {{ border: 1px solid transparent; border-radius: 8px; padding: {tool_pad}; margin: 0 1px; background: transparent; min-height: {btn_min}; }}
QToolButton:hover {{ background: {hover}; }}
QToolButton:checked {{ background: {tint}; color: {highlight}; border: 1px solid {tint}; }}

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QPlainTextEdit, QTextEdit {{
    background: {alt}; color: {text}; border: 1px solid {line}; border-radius: 6px;
    padding: {field_pad}; min-height: {field_min}; selection-background-color: {highlight}; selection-color: {onhl}; }}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {{ border: 1px solid {highlight}; }}
QComboBox::drop-down {{ border: none; width: 18px; }}
QComboBox QAbstractItemView {{ background: {base}; selection-background-color: {highlight};
    selection-color: {onhl}; border: 1px solid {line}; outline: none; }}

QTableView, QTreeView, QListView, QTableWidget, QListWidget, QTreeWidget {{
    background: {inset}; alternate-background-color: {stripe}; border: 1px solid {line}; border-radius: 6px;
    gridline-color: transparent;
    selection-background-color: {highlight}; selection-color: {onhl}; outline: none; }}
QTableView::indicator, QListView::indicator {{ width: 14px; height: 14px; border: 1px solid {mid}; border-radius: 3px;
    background: {base}; }}
QTableView::indicator:checked, QListView::indicator:checked {{ background: {highlight}; border: 1px solid {highlight}; }}
QHeaderView::section {{ background: {window}; color: {mid}; border: none;
    border-bottom: 1px solid {line}; padding: {header_pad}; }}

QGroupBox {{ border: none; border-top: 1px solid {line}; margin-top: 14px; padding-top: 8px; }}
QGroupBox::title {{ subcontrol-origin: margin; left: 0px; color: {text}; font-weight: 600; }}

QMenuBar {{ background: transparent; color: {text}; padding: 2px 4px; }}
QMenuBar::item {{ background: transparent; color: {text}; padding: 4px 8px; }}
QMenuBar::item:selected {{ background: {button}; color: {text}; border-radius: 4px; }}
QMenu {{ background: {base}; color: {text}; border: 1px solid {line}; padding: 4px; menu-scrollable: 1; }}
QMenu::item {{ color: {text}; padding: {menu_pad}; border-radius: 4px; }}
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


def current_density() -> str:
    """The text size in use: the one set on this application, else the saved or screen-based default."""
    app = QtWidgets.QApplication.instance()
    name = app.property("_pygcgc_density") if app is not None else None
    if name in DENSITY:
        return name
    from . import prefs
    return prefs.density()


def set_density(window, name: str, save: bool = True) -> None:
    """Change the text size of the whole application at once, and remember the choice."""
    if name not in DENSITY:
        raise ValueError(f"unknown text size '{name}'")
    if save:
        from . import prefs
        prefs.save_density(name)
    apply_theme(window, getattr(window, "dark", False), density=name)


def fit_to_screen(widget: QtWidgets.QWidget, width: int, height: int, fraction: float = 0.92,
                  centre: bool = False) -> None:
    """Resize a window to the size it wants, but never larger than `fraction` of the screen it is
    on. Dialogs with a parent are positioned by Qt; only top-level windows are centred on request."""
    screen = widget.screen() or QtGui.QGuiApplication.primaryScreen()
    if screen is None:
        widget.resize(width, height)
        return
    avail = screen.availableGeometry()
    w, h = min(width, int(avail.width() * fraction)), min(height, int(avail.height() * fraction))
    widget.resize(w, h)
    if centre:
        widget.move(avail.x() + (avail.width() - w) // 2, avail.y() + (avail.height() - h) // 2)


def apply_theme(window, dark: bool, density: str | None = None) -> None:
    c = DARK if dark else LIGHT
    app = QtWidgets.QApplication.instance()
    if density is not None:
        if density not in DENSITY:
            raise ValueError(f"unknown text size '{density}'")
        app.setProperty("_pygcgc_density", density)
    D = DENSITY[current_density()]
    stylesheet = _QSS.format(font=FONT, **c, **D)
    changed = app.styleSheet() != stylesheet
    # New windows inherit the application theme. Resetting it for each dialog
    # repolishes every existing widget, making repeated peak inspection slow.
    # Fusion is shared by both colour themes. Replacing the style on every
    # light/dark switch resets Qt's default font and needlessly traverses every
    # hidden graphics window; choose it once for the application instead.
    if not app.property('_pygcgc_fusion_selected'):
        app.setStyle("Fusion")
        app.setProperty('_pygcgc_fusion_selected', True)
    font = QtGui.QFont("Helvetica Neue" if "Helvetica Neue" in QtGui.QFontDatabase.families()
                        else "Arial")
    font.setPointSizeF(D["base_pt"])
    # A stylesheet can rewrite the family as a fallback list even though the
    # effective font is unchanged. Resetting that equivalent font recursively
    # propagates native FontChange events through hidden pyqtgraph windows.
    # Compare the resolved face/size/style instead of its request representation.
    current_font, requested_font = QtGui.QFontInfo(app.font()), QtGui.QFontInfo(font)
    if (abs(app.font().pointSizeF() - font.pointSizeF()) > 1e-6
            or current_font.family() != requested_font.family()
            or current_font.pointSizeF() != requested_font.pointSizeF()
            or current_font.weight() != requested_font.weight()
            or current_font.italic() != requested_font.italic()):
        app.setFont(font)
    if changed:
        app.setStyleSheet(stylesheet)
    # Qt's default layout spacing is cramped; give every panel room to breathe.
    # The spacing a layout asked for is remembered, so a switch to a tighter size can go back down.
    for lay in window.findChildren(QtWidgets.QLayout):
        if isinstance(lay, (QtWidgets.QBoxLayout, QtWidgets.QGridLayout, QtWidgets.QFormLayout)):
            original = lay.property("_orig_spacing")
            if original is None:
                original = lay.spacing()
                lay.setProperty("_orig_spacing", original)
            lay.setSpacing(max(int(original), D["spacing"]))
    # Panels grey their hints with "palette(mid)", which Qt resolves to a washed-out tint
    # under a stylesheet. Swap in the theme's own grey, keeping the original for the next switch.
    for w in window.findChildren(QtWidgets.QWidget):
        original = w.property("_ss_original")
        if original is None:
            original = w.styleSheet()
            if "palette(mid)" not in original:
                continue
            w.setProperty("_ss_original", original)
        w.setStyleSheet(original.replace("palette(mid)", c["mid"] + f"; font-size: {D['hint_pt']:g}pt"))
    if hasattr(window, "apply_text_scale"):
        window.apply_text_scale(D)
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
    if app.palette() != pal:
        app.setPalette(pal)
    bg = QtGui.QColor(c["plot"])
    bg.setAlpha(c["plot_alpha"])
    for view in (window.view2d, window.view1d, window.spectrum,
                 getattr(window, "calibration_panel", None), getattr(window, "batch_dialog", None),
                 getattr(window, "align_panel", None),
                 getattr(window, "overlap_dialog", None), getattr(window, "compare_window", None),
                 getattr(window, "anchor_window", None), getattr(window, "library_matcher", None)):
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
