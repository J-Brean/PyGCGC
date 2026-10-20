"""The window must work on a laptop screen: smaller text on request, nothing larger than the screen,
menus that fit, and drop-downs that do not live in a toolbar that can overflow."""

from __future__ import annotations

import pytest
from PySide6 import QtGui, QtWidgets

from pygcgc.gui import prefs
from pygcgc.gui.theme import DENSITY, current_density, fit_to_screen, set_density


@pytest.fixture
def main_window(qtbot):
    from pygcgc.gui.mainwindow import MainWindow
    w = MainWindow()
    qtbot.addWidget(w)
    yield w
    set_density(w, "standard", save=False)         # the application keeps the choice between tests


def all_menus(window):
    return window.menuBar().findChildren(QtWidgets.QMenu)


def test_compact_window_can_shrink_below_a_laptop_screen(main_window):
    set_density(main_window, "compact", save=False)
    hint = main_window.minimumSizeHint()
    assert hint.width() <= 1000 and hint.height() <= 680


def test_no_menu_is_taller_than_a_small_screen(main_window):
    set_density(main_window, "compact", save=False)
    tallest = max(m.sizeHint().height() for m in all_menus(main_window))
    assert tallest < 600


def test_process_menu_is_one_submenu_per_numbered_step_in_order(main_window):
    process = next(m for m in main_window.menuBar().findChildren(QtWidgets.QMenu) if m.title() == "&Process")
    titles = [a.text() for a in process.actions() if a.menu() is not None]
    assert [t.split()[0] for t in titles] == ["1", "2", "3", "4", "5", "6", "7", "8"]


def test_tools_and_panels_are_in_the_menu_bar_not_in_a_toolbar(main_window):
    assert main_window.tools_menu.menuAction() in main_window.menuBar().actions()
    assert main_window.act_compare in main_window.tools_menu.actions()
    for bar in main_window.findChildren(QtWidgets.QToolBar):
        popups = [b for b in bar.findChildren(QtWidgets.QToolButton)
                  if b.menu() is not None and b.text() in ("Tools", "Panels ▾")]
        assert not popups


def test_text_size_changes_the_application_font_and_the_workflow_title(main_window):
    set_density(main_window, "large", save=False)
    large = QtWidgets.QApplication.font().pointSizeF(), main_window.workflow.title.font().pointSizeF()
    set_density(main_window, "compact", save=False)
    compact = QtWidgets.QApplication.font().pointSizeF(), main_window.workflow.title.font().pointSizeF()
    assert compact[0] == DENSITY["compact"]["base_pt"] and large[0] == DENSITY["large"]["base_pt"]
    assert compact[1] < large[1]
    assert current_density() == "compact"


def test_layout_spacing_goes_back_up_after_a_compact_spell(main_window):
    layouts = [lay for lay in main_window.findChildren(QtWidgets.QBoxLayout)]
    set_density(main_window, "compact", save=False)
    compact = [lay.spacing() for lay in layouts]
    set_density(main_window, "large", save=False)
    large = [lay.spacing() for lay in layouts]
    assert sum(compact) < sum(large) and all(a <= b for a, b in zip(compact, large))


def test_changing_the_text_size_does_not_reselect_the_style_or_reset_an_unchanged_font(main_window):
    app = QtWidgets.QApplication.instance()
    set_density(main_window, "compact", save=False)
    style = app.style()
    set_density(main_window, "compact", save=False)
    assert app.style() is style and app.property("_pygcgc_fusion_selected")


def test_the_view_menu_lists_the_three_sizes_and_checks_the_current_one(main_window):
    set_density(main_window, "compact", save=False)
    main_window.density_actions["compact"].setChecked(True)
    assert set(main_window.density_actions) == {"compact", "standard", "large"}
    main_window.density_actions["large"].trigger()
    assert current_density() == "large" and main_window.density_actions["large"].isChecked()
    assert prefs.saved_density() == "large"
    prefs.save_density("standard")


def test_density_defaults_to_the_screen_and_the_saved_choice_wins(monkeypatch):
    class Screen:
        def __init__(self, h):
            self.h = h

        def availableGeometry(self):
            from PySide6 import QtCore
            return QtCore.QRect(0, 0, 1366, self.h)
    assert prefs.default_density(Screen(728)) == "compact"
    assert prefs.default_density(Screen(1040)) == "standard"
    with pytest.raises(ValueError):
        prefs.save_density("enormous")


def test_windows_are_never_larger_than_the_screen_they_open_on(qtbot):
    dialog = QtWidgets.QDialog()
    qtbot.addWidget(dialog)
    screen = QtGui.QGuiApplication.primaryScreen().availableGeometry()
    fit_to_screen(dialog, 5000, 4000)
    assert dialog.width() <= screen.width() and dialog.height() <= screen.height()
    fit_to_screen(dialog, 300, 200)
    assert (dialog.width(), dialog.height()) == (300, 200)


def test_main_window_geometry_is_remembered_and_clamped(qtbot):
    from pygcgc.gui.mainwindow import MainWindow
    w = MainWindow()
    qtbot.addWidget(w)
    w.resize(640, 480)
    w.close()
    assert prefs.saved_geometry("main") is not None
    again = MainWindow()
    qtbot.addWidget(again)
    screen = QtGui.QGuiApplication.primaryScreen().availableGeometry()
    assert again.width() <= screen.width() and again.height() <= screen.height()
