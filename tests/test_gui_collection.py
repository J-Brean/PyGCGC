"""Qt cycles must be destroyed on the application thread, not a worker."""

import gc
import os
import weakref

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
from PySide6 import QtCore

from pygcgc.gui.garbage_collection import ensure_gui_collection
from pygcgc.gui.jobs import JobManager


def test_worker_allocations_leave_widget_cycles_for_the_gui_thread(qtbot, qapp):
    jobs = JobManager()
    collector = ensure_gui_collection()
    assert collector is not None and not gc.isenabled()
    destroyed_on = []

    class Object(QtCore.QObject):
        def __del__(self):
            destroyed_on.append(QtCore.QThread.currentThread())

    obj = Object()
    obj.cycle = obj
    ref = weakref.ref(obj)
    del obj
    finished = []

    def work(token, progress):
        # Enough allocations to trigger automatic collection if it were enabled.
        values = [[i] for i in range(2000)]
        return len(values)

    jobs.submit("Allocate", work, finished.append)
    qtbot.waitUntil(lambda: bool(finished), timeout=10000)
    collector.collect()
    assert ref() is None
    assert destroyed_on == [qapp.thread()]


def test_collection_controller_is_shared_by_job_managers(qtbot, qapp):
    first = ensure_gui_collection()
    JobManager()
    JobManager()
    assert ensure_gui_collection() is first
