"""Autosave and the "processed on" record (pytest-qt, headless)."""

import json
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
pytest.importorskip("pyqtgraph")

from PySide6 import QtWidgets  # noqa: E402

from pygcgc.core import BackgroundSettings  # noqa: E402
from pygcgc.gui.project import AUTOSAVE_DELAY_MS, SIDECAR_NAME  # noqa: E402

TIMEOUT = 20000


def _new_window(qtbot):
    from pygcgc.gui.mainwindow import MainWindow
    w = MainWindow()
    w.jobs_panel.clear_after_ms = 10
    qtbot.addWidget(w)
    return w


def _idle(qtbot, w, loaded):
    qtbot.waitUntil(lambda: bool(loaded) and w.jobs.active_count() == 0, timeout=TIMEOUT)


def _choose(monkeypatch, label):
    """Answer the "saved processing found" box without showing it."""
    asked = []

    def fake_exec(box):
        asked.append(box.text())
        next(b for b in box.buttons() if b.text() == label).click()
        return 0

    monkeypatch.setattr(QtWidgets.QMessageBox, "exec", fake_exec)
    return asked


@pytest.fixture
def processed_folder(qtbot, synthetic_folder):
    """The synthetic folder after one window has cleaned run 0 and autosaved."""
    w = _new_window(qtbot)
    loaded = []
    w.open_folder(synthetic_folder, on_loaded=lambda: loaded.append(True))
    _idle(qtbot, w, loaded)
    w.project.compute_background([0], BackgroundSettings(), None, on_done=w.project.set_backgrounds)
    sidecar = synthetic_folder / SIDECAR_NAME
    qtbot.waitUntil(sidecar.exists, timeout=TIMEOUT)
    w.jobs.wait(5000)
    w.project.flush_autosave()
    w.undo.clear()
    w.close()
    return synthetic_folder


def test_processing_step_writes_sidecar(processed_folder):
    sidecar = processed_folder / SIDECAR_NAME
    data = json.loads(sidecar.read_text())
    assert data["processed_at"]
    assert data["runs"][0]["background"] is True
    assert not sidecar.with_name(sidecar.name + ".tmp").exists()


def test_resume_restores_processing_and_date(qtbot, monkeypatch, processed_folder):
    stamp = json.loads((processed_folder / SIDECAR_NAME).read_text())["processed_at"]
    asked = _choose(monkeypatch, "Resume processing")
    w = _new_window(qtbot)
    loaded = []
    w.open_folder(processed_folder, on_loaded=lambda: loaded.append(True))
    _idle(qtbot, w, loaded)
    assert asked and asked[0].startswith("This chromatogram was processed on ")
    assert w.project.path == processed_folder / SIDECAR_NAME
    assert w.project.entries[0].background is not None
    assert "processed on" in w.log.toPlainText()
    # Rebuilding the layers on load must not count as processing them again.
    qtbot.wait(AUTOSAVE_DELAY_MS + 300)
    assert w.project.processed_at == stamp
    assert json.loads((processed_folder / SIDECAR_NAME).read_text())["processed_at"] == stamp
    w.undo.clear()


def test_start_fresh_ignores_saved_processing(qtbot, monkeypatch, processed_folder):
    _choose(monkeypatch, "Start fresh")
    w = _new_window(qtbot)
    loaded = []
    w.open_folder(processed_folder, on_loaded=lambda: loaded.append(True))
    _idle(qtbot, w, loaded)
    assert w.project.path is None
    assert all(e.background is None for e in w.project.entries)
    w.undo.clear()


def test_cancel_opens_nothing(qtbot, monkeypatch, processed_folder):
    _choose(monkeypatch, "Cancel")
    w = _new_window(qtbot)
    w.open_folder(processed_folder)
    assert len(w.project) == 0
