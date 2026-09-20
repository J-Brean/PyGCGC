"""Background jobs.

The rule for the whole GUI: the main thread only draws and handles input.
Anything that might take noticeable time is a Job. A job runs a plain
Python function in a thread pool, reports progress through Qt signals
(which Qt delivers safely back on the main thread), and can be cancelled.

Jobs submitted with the same `latest_key` replace each other: when a new
one arrives, older ones still waiting are cancelled, and the result of an
older one that has already started is thrown away. This is what keeps
live previews smooth while you drag a slider.
"""

from __future__ import annotations

import itertools
import math
import os
import time
import traceback
from dataclasses import dataclass, field
from typing import Any, Callable

from PySide6 import QtCore, QtGui, QtWidgets


class Cancelled(Exception):
    """Raised inside a job function when its token has been cancelled."""


class CancelToken:
    def __init__(self) -> None:
        self._cancelled = False

    def cancel(self) -> None:
        self._cancelled = True

    @property
    def cancelled(self) -> bool:
        return self._cancelled

    def check(self) -> None:
        if self._cancelled:
            raise Cancelled()


class _Signals(QtCore.QObject):
    progress = QtCore.Signal(int, int, int, str)      # job id, done, total, message
    finished = QtCore.Signal(int, object)             # job id, result
    failed = QtCore.Signal(int, str)                  # job id, error text
    cancelled = QtCore.Signal(int)


class _Runnable(QtCore.QRunnable):
    def __init__(self, job: "JobInfo", signals: _Signals):
        super().__init__()
        self.job = job
        self.signals = signals
        self.setAutoDelete(True)

    def run(self) -> None:
        job = self.job
        if job.token.cancelled:
            self.signals.cancelled.emit(job.id)
            return
        job.started = time.monotonic()

        def progress(done: int, total: int, message: str = "") -> None:
            job.token.check()
            self.signals.progress.emit(job.id, int(done), int(total), message)

        try:
            result = job.fn(job.token, progress)
        except Cancelled:
            self.signals.cancelled.emit(job.id)
        except Exception as exc:  # noqa: BLE001  (report everything to the user)
            detail = "".join(traceback.format_exception_only(type(exc), exc)).strip()
            job.traceback = traceback.format_exc()
            self.signals.failed.emit(job.id, detail)
        else:
            if job.token.cancelled:
                self.signals.cancelled.emit(job.id)
            else:
                self.signals.finished.emit(job.id, result)


@dataclass
class JobInfo:
    id: int
    name: str
    fn: Callable[[CancelToken, Callable], Any]
    on_done: Callable[[Any], None] | None = None
    on_error: Callable[[str], None] | None = None
    latest_key: str | None = None
    quiet: bool = False                    # quiet jobs don't appear in the Jobs panel
    on_cancel: Callable[[], None] | None = None
    token: CancelToken = field(default_factory=CancelToken)
    state: str = "queued"                  # queued, running, done, failed, cancelled
    done: int = 0
    total: int = 0
    message: str = ""
    started: float = 0.0
    traceback: str = ""


class JobManager(QtCore.QObject):
    jobAdded = QtCore.Signal(object)
    jobChanged = QtCore.Signal(object)
    jobRemoved = QtCore.Signal(int)
    busyChanged = QtCore.Signal(int)       # number of active visible jobs

    def __init__(self, parent=None, max_threads: int | None = None):
        super().__init__(parent)
        self.pool = QtCore.QThreadPool(self)
        n = max_threads or max(1, (os.cpu_count() or 2) - 1)
        self.pool.setMaxThreadCount(n)
        self._ids = itertools.count(1)
        self.jobs: dict[int, JobInfo] = {}
        self._latest: dict[str, int] = {}
        self._signals = _Signals()
        self._signals.progress.connect(self._on_progress)
        self._signals.finished.connect(self._on_finished)
        self._signals.failed.connect(self._on_failed)
        self._signals.cancelled.connect(self._on_cancelled)

    # ----------------------------------------------------------------- API
    def submit(self, name: str, fn: Callable[[CancelToken, Callable], Any],
               on_done: Callable[[Any], None] | None = None,
               on_error: Callable[[str], None] | None = None,
               latest_key: str | None = None, quiet: bool = False,
               on_cancel: Callable[[], None] | None = None) -> int:
        job = JobInfo(next(self._ids), name, fn, on_done, on_error, latest_key, quiet)
        job.on_cancel = on_cancel
        if latest_key is not None:
            old = self._latest.get(latest_key)
            if old is not None and old in self.jobs:
                self.jobs[old].token.cancel()
            self._latest[latest_key] = job.id
        self.jobs[job.id] = job
        if not quiet:
            self.jobAdded.emit(job)
            self._emit_busy()
        self.pool.start(_Runnable(job, self._signals))
        return job.id

    def cancel(self, job_id: int) -> None:
        job = self.jobs.get(job_id)
        if job:
            job.token.cancel()

    def cancel_all(self) -> None:
        for job in self.jobs.values():
            job.token.cancel()

    def active_count(self) -> int:
        return sum(1 for j in self.jobs.values() if not j.quiet and j.state in ("queued", "running"))

    def wait(self, msecs: int = 30000) -> bool:
        """Block until the pool is idle (used by tests and on exit)."""
        return self.pool.waitForDone(msecs)

    # ------------------------------------------------------------ handlers
    def _emit_busy(self) -> None:
        self.busyChanged.emit(self.active_count())

    def _on_progress(self, job_id, done, total, message) -> None:
        job = self.jobs.get(job_id)
        if job is None:
            return
        job.state, job.done, job.total, job.message = "running", done, total, message
        if not job.quiet:
            self.jobChanged.emit(job)

    def _is_stale(self, job: JobInfo) -> bool:
        return job.latest_key is not None and self._latest.get(job.latest_key) != job.id

    def _finish(self, job: JobInfo, state: str) -> None:
        job.state = state
        if not job.quiet:
            self.jobChanged.emit(job)
            self._emit_busy()
        if job.latest_key is not None and self._latest.get(job.latest_key) == job.id:
            del self._latest[job.latest_key]
        # Keep finished visible jobs around briefly so the panel can show them.
        if job.quiet or state == "cancelled":
            self.jobs.pop(job.id, None)
            if not job.quiet:
                self.jobRemoved.emit(job.id)

    def _on_finished(self, job_id, result) -> None:
        job = self.jobs.get(job_id)
        if job is None:
            return
        if self._is_stale(job):
            self._finish(job, "cancelled")
            return
        # Callback before state change so callers can chain jobs.
        try:
            if job.on_done:
                job.on_done(result)
        finally:
            self._finish(job, "done")

    def _on_failed(self, job_id, text) -> None:
        job = self.jobs.get(job_id)
        if job is None:
            return
        job.message = text
        stale = self._is_stale(job)
        if job.on_error and not stale:
            job.on_error(text)
        self._finish(job, "cancelled" if stale else "failed")

    def _on_cancelled(self, job_id) -> None:
        job = self.jobs.get(job_id)
        if job is not None:
            if job.on_cancel:
                job.on_cancel()
            self._finish(job, "cancelled")

    def forget(self, job_id: int) -> None:
        if self.jobs.pop(job_id, None) is not None:
            self.jobRemoved.emit(job_id)


class _JobRow(QtWidgets.QWidget):
    def __init__(self, job: JobInfo, manager: JobManager, parent=None):
        super().__init__(parent)
        self.job_id = job.id
        self.manager = manager
        lay = QtWidgets.QGridLayout(self)
        lay.setContentsMargins(4, 2, 4, 2)
        self.title = QtWidgets.QLabel(job.name)
        self.detail = QtWidgets.QLabel("")
        self.detail.setStyleSheet("color: palette(mid);")
        self.bar = GCAnimation()
        self.bar.start()
        self.button = QtWidgets.QPushButton("Cancel")
        self.button.setFlat(True)
        self.button.clicked.connect(self._clicked)
        lay.addWidget(self.title, 0, 0)
        lay.addWidget(self.button, 0, 1, 2, 1)
        lay.addWidget(self.bar, 1, 0)
        lay.addWidget(self.detail, 2, 0, 1, 2)
        self.update_from(job)

    def _clicked(self) -> None:
        job = self.manager.jobs.get(self.job_id)
        if job is None or job.state in ("done", "failed", "cancelled"):
            self.manager.forget(self.job_id)
        else:
            self.manager.cancel(self.job_id)

    def update_from(self, job: JobInfo) -> None:
        if job.total:
            self.bar.set_progress(job.done / job.total)
        elif job.state in ("queued", "running"):
            self.bar.set_progress(None)   # busy indicator
        state_text = {"queued": "Waiting", "running": job.message or "Running",
                      "done": "Done", "failed": f"Failed: {job.message}",
                      "cancelled": "Cancelled"}[job.state]
        self.detail.setText(state_text)
        if job.state in ("done", "failed"):
            self.bar.set_progress(1.0 if job.state == "done" else 0.0)
            self.bar.stop()
            self.button.setText("Clear")


class GCAnimation(QtWidgets.QWidget):
    """Pixel-art progress bar: molecules race along a coiled GC column in a warming oven.

    Drawn on an 8-row pixel grid, as wide as the widget, and scaled up with hard
    edges. The bottom rail of the oven lights up left to right with progress;
    with no progress value (None) it is lit throughout. The timer only runs
    between start() and stop(), so finished rows cost nothing.
    """

    GRID_H = 8
    HEAT = [(60, 110, 200), (70, 170, 190), (230, 200, 70), (240, 140, 50), (220, 70, 50)]
    # (speed in pixels per frame, phase offset in pixels, colour): light species outrun heavy ones
    MOLECULES = [(0.90, 0, (255, 255, 255)), (0.62, 20, (120, 220, 255)),
                 (0.42, 45, (255, 120, 200)), (0.27, 8, (160, 255, 120))]
    COIL_PERIOD = 12

    def __init__(self, parent=None, scale: int = 2):
        super().__init__(parent)
        self.scale = scale
        self.frame = 0
        self.progress: float | None = None
        self.setFixedHeight(self.GRID_H * scale)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(70)
        self.timer.timeout.connect(self._tick)

    def start(self) -> None:
        if not self.timer.isActive():
            self.timer.start()

    def stop(self) -> None:
        self.timer.stop()
        self.update()

    def set_progress(self, fraction: float | None) -> None:
        self.progress = None if fraction is None else min(1.0, max(0.0, fraction))
        self.update()

    def sizeHint(self) -> QtCore.QSize:
        return QtCore.QSize(200, self.GRID_H * self.scale)

    def _tick(self) -> None:
        self.frame += 1
        self.update()

    def _coil_y(self, x: float) -> int:
        return int(round(3.5 + 2.0 * math.sin(2 * math.pi * x / self.COIL_PERIOD)))

    def _heat(self, t: float) -> QtGui.QColor:
        # t in [0, 1) walks up the ramp and back down, so the oven warms and cools in a loop
        t = 1 - abs(2 * (t % 1) - 1)
        pos = t * (len(self.HEAT) - 1)
        i = min(int(pos), len(self.HEAT) - 2)
        f = pos - i
        a, b = self.HEAT[i], self.HEAT[i + 1]
        return QtGui.QColor(*(int(a[k] + (b[k] - a[k]) * f) for k in range(3)))

    def paintEvent(self, _event) -> None:
        p = QtGui.QPainter(self)
        p.scale(self.scale, self.scale)

        def px(x: int, y: int, colour) -> None:
            p.fillRect(int(x), int(y), 1, 1, colour)

        grid_w = max(16, self.width() // self.scale)
        heat = self._heat(self.frame / 160)
        dim = QtGui.QColor(heat)
        dim.setAlpha(90)
        off = QtGui.QColor(128, 128, 128, 60)
        lit = grid_w if self.progress is None else int(round(self.progress * grid_w))
        # Oven walls: top and bottom rails with a flickering edge; the bottom one is the progress readout
        for x in range(grid_w):
            px(x, 0, heat if (x + self.frame // 3) % 4 else dim)
            if x < lit:
                px(x, self.GRID_H - 1, heat if (x - self.frame // 3) % 4 else dim)
            else:
                px(x, self.GRID_H - 1, off)
        # Column coil
        coil = QtGui.QColor(150, 150, 160, 200)
        for x in range(1, grid_w - 1):
            px(x, self._coil_y(x), coil)
        # Molecules
        span = grid_w - 2
        for speed, phase, rgb in self.MOLECULES:
            x = (phase + speed * self.frame) % span + 1
            y = self._coil_y(x)
            colour = QtGui.QColor(*rgb)
            px(x, y, colour)
            px(x + 1, y, colour)
            trail = QtGui.QColor(*rgb, 90)
            px(x - 1, self._coil_y(x - 1), trail)
        p.end()


class JobsPanel(QtWidgets.QWidget):
    """List of visible jobs with progress bars and cancel buttons.
    Finished jobs clear themselves after a few seconds; failures stay until cleared."""

    def __init__(self, manager: JobManager, parent=None, clear_after_ms: int = 4000):
        super().__init__(parent)
        self.manager = manager
        self.clear_after_ms = clear_after_ms
        self.rows: dict[int, _JobRow] = {}
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        top = QtWidgets.QHBoxLayout()
        self.summary = QtWidgets.QLabel("No jobs running")
        cancel_all = QtWidgets.QPushButton("Cancel all")
        cancel_all.clicked.connect(manager.cancel_all)
        top.addWidget(self.summary)
        top.addStretch(1)
        top.addWidget(cancel_all)
        outer.addLayout(top)
        self.list = QtWidgets.QVBoxLayout()
        self.list.addStretch(1)
        holder = QtWidgets.QWidget()
        holder.setLayout(self.list)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(holder)
        outer.addWidget(scroll)
        manager.jobAdded.connect(self._added)
        manager.jobChanged.connect(self._changed)
        manager.jobRemoved.connect(self._removed)
        manager.busyChanged.connect(self._busy)

    def _added(self, job: JobInfo) -> None:
        row = _JobRow(job, self.manager)
        self.rows[job.id] = row
        self.list.insertWidget(self.list.count() - 1, row)

    def _changed(self, job: JobInfo) -> None:
        row = self.rows.get(job.id)
        if row:
            row.update_from(job)
        if job.state == "done":
            # Tie the timer to this panel so it dies with it.
            QtCore.QTimer.singleShot(self.clear_after_ms, self, lambda i=job.id: self.manager.forget(i))

    def _removed(self, job_id: int) -> None:
        row = self.rows.pop(job_id, None)
        if row:
            row.setParent(None)
            row.deleteLater()

    def _busy(self, n: int) -> None:
        self.summary.setText("No jobs running" if n == 0 else f"{n} job{'s' if n > 1 else ''} running")