"""Period and phase estimated when runs are opened (pytest-qt, headless)."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
pytest.importorskip("pyqtgraph")

from pygcgc.core import circstats as cs  # noqa: E402
from tests.conftest import write_run_andi  # noqa: E402
from tests.test_period import _run  # noqa: E402

TIMEOUT = 30000
TRUE_PERIOD, TRUE_PHASE = 6.0, 1.7


def _folder(tmp_path, method_period=None):
    folder = tmp_path / "pack"
    (folder / "data").mkdir(parents=True)
    write_run_andi(_run(TRUE_PERIOD, phase=TRUE_PHASE, duration=600.0), folder / "data" / "RUN.cdf")
    if method_period is not None:
        (folder / "method.json").write_text(f'{{"modulation_period_s": {method_period}}}')
    return folder


def _open(qtbot, folder):
    from pygcgc.gui.mainwindow import MainWindow
    w = MainWindow()
    w.jobs_panel.clear_after_ms = 10
    qtbot.addWidget(w)
    w.resize(1400, 900)
    w.show()
    messages = []
    w.project.message.connect(messages.append)
    loaded = []
    w.open_folder(folder, on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and w.jobs.active_count() == 0, timeout=TIMEOUT)
    qtbot.waitUntil(lambda: all(e.phase_method != "unset" for e in w.project.entries), timeout=TIMEOUT)
    return w, messages


def _close(w):
    w.jobs.cancel_all()
    w.jobs.wait(5000)
    w.undo.clear()


def test_period_and_phase_estimated_on_open_without_a_method_file(qtbot, tmp_path):
    w, messages = _open(qtbot, _folder(tmp_path))
    try:
        p = w.project
        assert p.period == pytest.approx(TRUE_PERIOD, abs=0.01)
        assert p.period_method == "auto" and not p.period_confirmed
        assert "estimated from the run" in p.period_detail
        assert any("Modulation period estimated" in m for m in messages)
        e = p.entries[0]
        assert e.phase_method == "void"
        # The cut lies in the quiet stretch, so it is not near a peak's t2 (1.1 to 3.4 s, shifted by the true phase).
        t2_unshifted = [(t2 + TRUE_PHASE) % TRUE_PERIOD for t2 in (0.6, 3.4)]
        assert min(abs(cs.circ_diff(e.phase, t, TRUE_PERIOD)) for t in t2_unshifted) > 0.1
        assert "estimated from the run" in w.fold_panel.period_status.text()
    finally:
        _close(w)


def test_confirmed_period_is_kept_and_a_disagreement_is_reported(qtbot, tmp_path):
    w, messages = _open(qtbot, _folder(tmp_path, method_period=4.0))
    try:
        p = w.project
        assert p.period == 4.0 and p.period_method == "file" and p.period_confirmed
        qtbot.waitUntil(lambda: any("The run suggests a period" in m for m in messages), timeout=TIMEOUT)
    finally:
        _close(w)


def test_estimate_button_replaces_a_confirmed_period_and_can_be_undone(qtbot, tmp_path):
    w, _ = _open(qtbot, _folder(tmp_path, method_period=4.0))
    try:
        p = w.project
        before = p.phase_snapshot()
        w.estimate_period()
        qtbot.waitUntil(lambda: p.period_method == "auto", timeout=TIMEOUT)
        assert p.period == pytest.approx(TRUE_PERIOD, abs=0.01)
        w.undo.undo()
        assert p.period == 4.0 and p.period_method == "file" and p.period_confirmed
        assert p.phase_snapshot() == before
        w.undo.redo()
        assert p.period == pytest.approx(TRUE_PERIOD, abs=0.01) and p.period_method == "auto"
    finally:
        _close(w)


def test_period_estimate_survives_save_and_load(qtbot, tmp_path):
    w, _ = _open(qtbot, _folder(tmp_path))
    try:
        path = tmp_path / "proj.json"
        w.project.save(path)
        import json
        data = json.loads(path.read_text())
        assert data["period_method"] == "auto" and data["period_confirmed"] is False
    finally:
        _close(w)
