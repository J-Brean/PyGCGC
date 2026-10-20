"""Peak names on the 2D image: hover, optional labels for every peak, and renaming by double-click."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from PySide6 import QtCore, QtWidgets  # noqa: E402

from pygcgc.gui.peak_naming import PeakNameDialog, label_for, rename_peak  # noqa: E402
from tests.test_gui import TIMEOUT, _wait_render, window  # noqa: E402,F401
from tests.test_library_matcher import matched_window  # noqa: E402


def ready(qtbot, w, tmp_path):
    win, state = matched_window(qtbot, w, tmp_path)
    _wait_render(qtbot, w, "RUN_A")
    qtbot.waitUntil(lambda: w.view2d.result is not None and w.view2d._label_geometry is not None, timeout=TIMEOUT)
    return win, state


class Event:
    def __init__(self, view, x, y, double=False):
        self._view, self._x, self._y, self._double = view, x, y, double

    def button(self):
        return QtCore.Qt.MouseButton.LeftButton

    def double(self):
        return self._double

    def scenePos(self):
        return self._view.plot.vb.mapViewToScene(QtCore.QPointF(self._x, self._y))

    def modifiers(self):
        return QtCore.Qt.KeyboardModifier.NoModifier


def peak_position(view, k=0):
    x, y, _rx, _ry = view._label_geometry
    return float(x[k]), float(y[k])


def test_hover_labels_are_on_by_default_and_every_peak_label_is_off(qtbot, window, tmp_path):
    w = window
    ready(qtbot, w, tmp_path)
    assert w.act_peak_labels.isChecked() and not w.act_peak_labels_all.isChecked()
    assert w.view2d.peak_label_hover and not w.view2d.peak_label_all
    n = len(w.project.matcher.results_for(w.index).matches)
    assert len(w.view2d._peak_labels) == n > 0 and w.view2d._label_items == []       # nothing is drawn until you hover
    assert not w.view2d.hover_label.isVisible()


def test_the_label_appears_over_a_peak_and_goes_away_off_it(qtbot, window, tmp_path):
    w = window
    ready(qtbot, w, tmp_path)
    view = w.view2d
    x, y = peak_position(view, 0)
    view._mouse_moved(Event(view, x, y).scenePos())
    key = view.hovered_peak_key()
    assert key is not None and view._peak_labels[0]["key"] == key
    assert view.hover_label.toPlainText().startswith(view._peak_labels[0]["text"].split("\n")[0])
    assert "double-click to rename" in view.hover_label.toPlainText()
    view._mouse_moved(Event(view, x + 5.0, y).scenePos())                           # 5 minutes along, still on the image
    assert view.hovered_peak_key() is None and not view.hover_label.isVisible()
    view._mouse_moved(Event(view, x, y).scenePos())
    assert view.hovered_peak_key() is not None
    view._mouse_moved(QtCore.QPointF(-500.0, -500.0))                                # off the plot altogether
    assert view.hovered_peak_key() is None


def test_hover_labelling_can_be_switched_off_and_every_peak_labelled_instead(qtbot, window, tmp_path):
    w = window
    ready(qtbot, w, tmp_path)
    view = w.view2d
    x, y = peak_position(view, 0)
    w.act_peak_labels.setChecked(False)
    view._mouse_moved(Event(view, x, y).scenePos())
    assert view.hovered_peak_key() is None                                           # off means off
    w.act_peak_labels_all.setChecked(True)
    assert len(view._label_items) == len(view._peak_labels) > 0
    w.act_peak_labels_all.setChecked(False)
    assert view._label_items == []
    w.act_peak_labels.setChecked(False)
    assert view._peak_labels == []


def test_a_double_click_on_a_peak_asks_to_rename_it_and_elsewhere_does_nothing(qtbot, window, tmp_path, monkeypatch):
    w = window
    ready(qtbot, w, tmp_path)
    monkeypatch.setattr("pygcgc.gui.peak_naming.rename_peak", lambda *a: False)      # the real dialog would block
    view = w.view2d
    heard = []
    view.peakLabelActivated.connect(heard.append)
    x, y = peak_position(view, 0)
    view._mouse_clicked(Event(view, x + 6.0, y, double=True))
    assert heard == []
    view._mouse_clicked(Event(view, x, y, double=True))
    assert heard == [view._peak_labels[0]["key"]]


def test_the_main_window_opens_the_naming_dialog_unless_something_is_being_picked(qtbot, window, tmp_path, monkeypatch):
    w = window
    ready(qtbot, w, tmp_path)
    opened = []
    monkeypatch.setattr("pygcgc.gui.peak_naming.rename_peak", lambda parent, project, index, key: opened.append((index, key)) or True)
    key = w.view2d._peak_labels[0]["key"]
    w._rename_peak_label(key)
    assert opened == [(w.index, key)]
    w._blob_add_pending = True
    w._rename_peak_label(key)
    assert len(opened) == 1                                                          # a click while picking is a pick, not a rename
    w._blob_add_pending = False


def dialog_returning(monkeypatch, typed, accept=True):
    def run(self):
        self.edit.setText(typed)
        return 1 if accept else 0
    monkeypatch.setattr(PeakNameDialog, "exec", run)


def target_label(w, state):
    found = state.current_for(w.index)
    return min(found.matches, key=lambda m: abs(m.query.t1 - 40.0)).query.label


def test_a_typed_name_is_checked_and_saved_with_the_library_spelling(qtbot, window, tmp_path, monkeypatch):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    label = target_label(w, state)
    dialog_returning(monkeypatch, "  COMPOUND 51 ")
    assert rename_peak(w, w.project, w.index, label)
    d = state.decisions_for(w.index)[label]
    assert d["action"] == "named" and d["name"] == "Compound51" and d["typed"] == "COMPOUND 51"
    match = next(m for m in state.current_for(w.index).matches if m.query.label == label)
    text, _rgb = label_for(match, d)
    assert text.startswith("Compound51") and "named by you" in text
    qtbot.waitUntil(lambda: next(i for i in w.view2d._peak_labels if i["key"] == label)["text"].startswith("Compound51"), timeout=5000)
    assert match.level == next(m for m in state.current_for(w.index).matches if m.query.label == label).level


def test_a_name_the_library_does_not_have_is_kept_as_typed_and_a_clash_is_noticed(qtbot, window, tmp_path, monkeypatch):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    found = state.current_for(w.index)
    first, second = found.matches[0].query.label, found.matches[1].query.label
    dialog_returning(monkeypatch, "my   own  compound")
    rename_peak(w, w.project, w.index, first)
    assert state.decisions_for(w.index)[first]["name"] == "my own compound" and state.decisions_for(w.index)[first]["cas"] == ""
    seen = {}
    original = PeakNameDialog._changed

    def spy(self, *a):
        original(self, *a)
        seen["notes"] = list(self.check.notes)

    monkeypatch.setattr(PeakNameDialog, "_changed", spy)
    dialog_returning(monkeypatch, "MY OWN COMPOUND")
    rename_peak(w, w.project, w.index, second)
    assert any("already called" in n for n in seen["notes"])


def test_cancelling_changes_nothing_and_removing_my_name_goes_back_to_the_match(qtbot, window, tmp_path, monkeypatch):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    label = target_label(w, state)
    dialog_returning(monkeypatch, "Something else", accept=False)
    assert rename_peak(w, w.project, w.index, label) is False and label not in state.decisions_for(w.index)
    dialog_returning(monkeypatch, "Compound51")
    rename_peak(w, w.project, w.index, label)
    assert state.decisions_for(w.index)[label]["action"] == "named"
    monkeypatch.setattr(PeakNameDialog, "exec", lambda self: (setattr(self, "cleared", True), 1)[1])
    assert rename_peak(w, w.project, w.index, label) and label not in state.decisions_for(w.index)


def test_renaming_before_any_matches_explains_instead_of_failing(qtbot, window, monkeypatch):
    w = window
    said = []
    monkeypatch.setattr(QtWidgets.QMessageBox, "information", lambda *a, **k: said.append(a[2]))
    assert rename_peak(w, w.project, 0, "x") is False
    assert said and "Match this run's peaks" in said[0]


def test_a_named_peak_shows_in_the_matcher_table_and_the_export(qtbot, window, tmp_path, monkeypatch):
    import csv
    w = window
    win, state = ready(qtbot, w, tmp_path)
    label = target_label(w, state)
    dialog_returning(monkeypatch, "compound-51")
    rename_peak(w, w.project, w.index, label)
    win._fill_results()
    names = [win.table.item(r, 4).text() for r in range(win.table.rowCount())]
    assert "Compound51" in names
    out = tmp_path / "ids.csv"
    win.export_to(out)
    rows = {r["peak"]: r for r in csv.DictReader(out.open(encoding="utf-8"))}
    assert rows[label]["decision"] == "named" and rows[label]["decision_name"] == "Compound51"


def test_a_named_peak_stays_named_only_while_the_matches_are_the_same(qtbot, window, tmp_path, monkeypatch):
    from dataclasses import replace
    w = window
    win, state = ready(qtbot, w, tmp_path)
    label = target_label(w, state)
    dialog_returning(monkeypatch, "Compound51")
    rename_peak(w, w.project, w.index, label)
    assert label in state.decisions_for(w.index)
    state.set_settings(replace(state.settings, ri_gate_z=2.5))                       # different evidence
    assert label not in state.decisions_for(w.index)


def test_the_dialog_explains_what_it_checked_and_blocks_an_empty_name(qtbot, window, tmp_path):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    found = state.current_for(w.index)
    match = found.matches[0]
    dialog = PeakNameDialog(state.library, match, "RUN_A", 0.1, [], None)
    qtbot.addWidget(dialog)
    dialog.edit.setText("compound 50")
    assert "Library spelling used" in dialog.report.text() and "Compound50" in dialog.report.text()
    dialog.edit.setText("Compound50")
    assert "In the library" in dialog.report.text()
    dialog.edit.setText("Compound5")
    assert "Not in the library" in dialog.report.text() and "Did you mean" in dialog.report.text()
    dialog.edit.setText("   ")
    assert not dialog.buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Save).isEnabled()
