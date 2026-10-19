"""Deterministic checks for worker results already queued when Cancel is clicked."""

import os
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
from pygcgc.gui.jobs import JobManager


@pytest.mark.parametrize("fails", [False, True])
def test_cancel_after_worker_finished_still_suppresses_queued_result(qtbot, fails):
    manager = JobManager()
    results, errors, cancellations = [], [], []

    def work(token, progress):
        if fails:
            raise ValueError("finished before cancellation")
        return 42

    job_id = manager.submit("race", work, results.append, errors.append,
                            on_cancel=lambda: cancellations.append(True))
    assert manager.wait(5000)  # Signals are queued; the GUI has not received them.
    manager.cancel(job_id)
    qtbot.waitUntil(lambda: manager.active_count() == 0)
    assert not results and not errors
    assert cancellations == [True]


def test_superseded_queued_result_notifies_cancellation(qtbot):
    manager = JobManager()
    results, cancellations = [], []
    manager.submit("old", lambda *_: 1, results.append, latest_key="test",
                   on_cancel=lambda: cancellations.append("old"))
    assert manager.wait(5000)
    manager.submit("new", lambda *_: 2, results.append, latest_key="test")
    qtbot.waitUntil(lambda: manager.active_count() == 0)
    assert results == [2] and cancellations == ["old"]


def test_late_progress_cannot_reactivate_a_settled_job(qtbot):
    manager = JobManager()
    job_id = manager.submit("complete", lambda *_: 42)
    qtbot.waitUntil(lambda: manager.active_count() == 0)
    manager._on_progress(job_id, 1, 2, "late")
    assert manager.active_count() == 0


@pytest.mark.parametrize("ending", ["failed", "cancelled"])
def test_job_settles_even_if_an_error_or_cancellation_callback_raises(qtbot, ending):
    manager = JobManager()
    def bad_callback(*args):
        raise RuntimeError("callback failed")
    job_id = manager.submit("callback", lambda *_: 42, on_error=bad_callback, on_cancel=bad_callback)
    assert manager.wait(5000)
    with pytest.raises(RuntimeError, match="callback failed"):
        if ending == "failed":
            manager._on_failed(job_id, "worker failed")
        else:
            manager._on_cancelled(job_id)
    assert manager.active_count() == 0
    qtbot.wait(10)  # The already queued success signal must now be ignored.
