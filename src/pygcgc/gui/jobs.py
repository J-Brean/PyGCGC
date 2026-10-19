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

from .icons import plain
from .garbage_collection import ensure_gui_collection


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
        ensure_gui_collection()
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
        if job is None or job.state not in ("queued", "running") or job.token.cancelled or self._is_stale(job):
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
        if job is None or job.state not in ("queued", "running"):
            return
        if job.token.cancelled or self._is_stale(job):
            self._on_cancelled(job_id)
            return
        # Callback before state change so callers can chain jobs.
        try:
            if job.on_done:
                job.on_done(result)
        finally:
            self._finish(job, "done")

    def _on_failed(self, job_id, text) -> None:
        job = self.jobs.get(job_id)
        if job is None or job.state not in ("queued", "running"):
            return
        if job.token.cancelled or self._is_stale(job):
            self._on_cancelled(job_id)
            return
        job.message = text
        try:
            if job.on_error:
                job.on_error(text)
        finally:
            self._finish(job, "failed")

    def _on_cancelled(self, job_id) -> None:
        job = self.jobs.get(job_id)
        if job is not None and job.state in ("queued", "running"):
            try:
                if job.on_cancel:
                    job.on_cancel()
            finally:
                self._finish(job, "cancelled")

    def forget(self, job_id: int) -> None:
        if self.jobs.pop(job_id, None) is not None:
            self.jobRemoved.emit(job_id)


def _duration(seconds: float) -> str:
    seconds = max(0.0, seconds)
    if seconds < 1:
        return "under 1 s"
    if seconds < 90:
        return f"{seconds:.0f} s"
    return f"{int(seconds // 60)} min {int(seconds % 60):02d} s"


class GCAnimation(QtWidgets.QWidget):
    """A gas chromatography column in miniature: sample molecules race along a coiled column.

    Light molecules outrun heavy ones, as they do in the oven. A thin track underneath fills
    left to right with progress, or sweeps back and forth when the length of the job is not
    known (progress None). Drawn with the palette's text and highlight colours, so it follows
    the theme. The timer only runs between start() and stop(), so finished rows cost nothing.
    """

    HEIGHT = 38
    # (share of the column covered per second, starting position): light species first
    MOLECULES = [(0.42, 0.00), (0.30, 0.37), (0.21, 0.66), (0.14, 0.12)]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.frame = 0
        self.progress: float | None = None
        self.setFixedHeight(self.HEIGHT)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(33)
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
        return QtCore.QSize(200, self.HEIGHT)

    def _tick(self) -> None:
        self.frame += 1
        self.update()

    def paintEvent(self, _event) -> None:
        pal = self.palette()
        accent = pal.color(QtGui.QPalette.ColorRole.Highlight)
        ink = pal.color(QtGui.QPalette.ColorRole.Text)
        grey = pal.color(QtGui.QPalette.ColorRole.Mid)
        p = QtGui.QPainter(self)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        w = self.width()
        # Coil: a helix seen from the side, so each turn is a loop. x runs left to right.
        loops = max(4, w // 30)
        x0, x1, cy, amp_y, amp_x = 4.0, w - 4.0, 15.0, 9.0, 7.0
        turns = 2 * math.pi * loops

        def point(u: float) -> QtCore.QPointF:
            theta = u * turns
            k = (x1 - x0 - 2 * amp_x) / turns
            return QtCore.QPointF(x0 + amp_x + k * theta - amp_x * math.sin(theta),
                                  cy + amp_y * math.cos(theta))

        path = QtGui.QPainterPath(point(0))
        steps = loops * 14
        for i in range(1, steps + 1):
            path.lineTo(point(i / steps))
        grey_line = QtGui.QColor(grey)
        grey_line.setAlpha(150)
        p.setPen(QtGui.QPen(grey_line, 1.3, QtCore.Qt.PenStyle.SolidLine, QtCore.Qt.PenCapStyle.RoundCap))
        p.drawPath(path)
        # Molecules, each with a short fading trail
        colours = [accent, ink, accent.lighter(135), grey]
        t = self.frame * self.timer.interval() / 1000
        p.setPen(QtCore.Qt.PenStyle.NoPen)
        for (speed, start), colour in zip(self.MOLECULES, colours):
            u = (start + speed * t) % 1.0
            for j in range(4, 0, -1):
                trail = QtGui.QColor(colour)
                trail.setAlpha(28 * (5 - j))
                p.setBrush(trail)
                p.drawEllipse(point((u - j * 0.007) % 1.0), 1.5, 1.5)
            p.setBrush(colour)
            p.drawEllipse(point(u), 2.7, 2.7)
        # Progress track
        track = QtCore.QRectF(2, self.HEIGHT - 5, w - 4, 2.5)
        faint = QtGui.QColor(grey)
        faint.setAlpha(55)
        p.setBrush(faint)
        p.drawRoundedRect(track, 1.25, 1.25)
        p.setBrush(accent)
        if self.progress is None:
            swing = abs(((t * 0.5) % 1.0) * 2 - 1)
            seg = track.width() * 0.25
            p.drawRoundedRect(QtCore.QRectF(track.left() + swing * (track.width() - seg), track.top(),
                                            seg, track.height()), 1.25, 1.25)
        else:
            p.drawRoundedRect(QtCore.QRectF(track.left(), track.top(), track.width() * self.progress,
                                            track.height()), 1.25, 1.25)
        p.end()


class _JobRow(QtWidgets.QWidget):
    """One running task: its name and "done/total", the animation, and what is left to do."""

    def __init__(self, job: JobInfo, manager: JobManager, parent=None):
        super().__init__(parent)
        self.job_id = job.id
        self.manager = manager
        self._samples: list[tuple[float, int]] = []     # (time, items done) at each progress report
        lay = QtWidgets.QGridLayout(self)
        lay.setContentsMargins(0, 2, 0, 2)
        lay.setVerticalSpacing(4)
        self.title = QtWidgets.QLabel(plain(job.name))
        font = self.title.font()
        font.setWeight(QtGui.QFont.Weight.DemiBold)
        self.title.setFont(font)
        self.count = QtWidgets.QLabel("")
        self.detail = QtWidgets.QLabel("")
        self.detail.setStyleSheet("color: palette(mid); font-size: 8pt;")
        self.detail.setWordWrap(True)
        self.bar = GCAnimation()
        self.bar.start()
        self.button = QtWidgets.QPushButton("Cancel")
        self.button.setFlat(True)
        self.button.clicked.connect(self._clicked)
        lay.addWidget(self.title, 0, 0)
        lay.addWidget(self.count, 0, 1)
        lay.addWidget(self.button, 0, 2)
        lay.addWidget(self.bar, 1, 0, 1, 3)
        lay.addWidget(self.detail, 2, 0, 1, 3)
        lay.setColumnStretch(0, 1)
        # ETA text changes with the clock, not only when the job reports progress
        self.clock = QtCore.QTimer(self)
        self.clock.setInterval(500)
        self.clock.timeout.connect(self._refresh)
        self.clock.start()
        self.update_from(job)

    def _clicked(self) -> None:
        job = self.manager.jobs.get(self.job_id)
        if job is None or job.state in ("done", "failed", "cancelled"):
            self.manager.forget(self.job_id)
        else:
            self.manager.cancel(self.job_id)

    def _refresh(self) -> None:
        job = self.manager.jobs.get(self.job_id)
        if job is not None:
            self.update_from(job)

    def _remaining(self, job: JobInfo, now: float) -> str:
        """Time left, from the speed over the last 20 s of progress reports.

        Not elapsed / done: many jobs spend a long first stage before their first report, and
        later items cost more or less than earlier ones, so the whole-job average is far off.
        """
        if not self._samples or self._samples[-1][1] != job.done:
            self._samples.append((now, job.done))
        while len(self._samples) > 2 and now - self._samples[1][0] > 20:
            self._samples.pop(0)
        (t0, d0), (t1, d1) = self._samples[0], self._samples[-1]
        if d1 <= d0 or t1 - t0 < 1.0:
            return "estimating time"
        left = (job.total - job.done) * (t1 - t0) / (d1 - d0) - (now - t1)
        left = max(left, 0.0)
        if left > 30:
            left = round(left / 5) * 5
        return f"about {_duration(left)} remaining"

    def update_from(self, job: JobInfo) -> None:
        active = job.state in ("queued", "running")
        if job.total:
            self.bar.set_progress(job.done / job.total)
            self.count.setText(f"{job.done}/{job.total}")
        elif active:
            self.bar.set_progress(None)   # length unknown: sweeping indicator
        if job.state == "queued":
            text = "Waiting"
        elif job.state == "running":
            text = job.message or "Running"
            now = time.monotonic()
            if job.total and job.done >= job.total:
                text += "  ·  finishing"
            elif job.total:
                text += f"  ·  {self._remaining(job, now)}"
            elif job.started:
                text += f"  ·  {_duration(now - job.started)} so far"
        else:
            text = {"done": "Done", "failed": f"Failed: {job.message}", "cancelled": "Cancelled"}[job.state]
        self.detail.setText(text)
        if job.state in ("done", "failed"):
            self.bar.set_progress(1.0 if job.state == "done" else 0.0)
            self.count.setText("")
            self.bar.stop()
            self.clock.stop()
            self.button.setText("Clear")


class JobsPanel(QtWidgets.QWidget):
    """The running tasks, shown at the foot of the workspace panel and hidden when there are none.
    Finished jobs clear themselves after a few seconds; failures stay until cleared."""

    def __init__(self, manager: JobManager, parent=None, clear_after_ms: int = 4000):
        super().__init__(parent)
        self.manager = manager
        self.clear_after_ms = clear_after_ms
        self.rows: dict[int, _JobRow] = {}
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 8, 0, 0)
        self.cancel_all = QtWidgets.QPushButton("Cancel all")
        self.cancel_all.clicked.connect(manager.cancel_all)
        self.cancel_all.setVisible(False)
        outer.addWidget(self.cancel_all, 0, QtCore.Qt.AlignmentFlag.AlignRight)
        self.list = QtWidgets.QVBoxLayout()
        self.list.setSpacing(10)
        holder = QtWidgets.QWidget()
        holder.setLayout(self.list)
        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.scroll.setMaximumHeight(190)
        self.scroll.setWidget(holder)
        outer.addWidget(self.scroll)
        self.setVisible(False)
        manager.jobAdded.connect(self._added)
        manager.jobChanged.connect(self._changed)
        manager.jobRemoved.connect(self._removed)

    def _sync(self) -> None:
        self.setVisible(bool(self.rows))
        self.cancel_all.setVisible(len(self.rows) > 1)
        # Grow with the rows up to the cap, rather than leaving a tall empty strip
        self.scroll.setFixedHeight(min(190, sum(r.sizeHint().height() + 10 for r in self.rows.values())))

    def _added(self, job: JobInfo) -> None:
        row = _JobRow(job, self.manager)
        self.rows[job.id] = row
        self.list.addWidget(row)
        self._sync()

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
        self._sync()
