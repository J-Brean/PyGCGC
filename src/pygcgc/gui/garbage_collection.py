"""Keep cyclic garbage collection, and Qt object destruction, on the GUI thread.

Python can start automatic collection in any allocating thread. With PySide,
collecting widget cycles from a worker can deadlock against deletion on the
GUI thread. A timer collects those cycles on the application's thread instead.
Reference-counted objects are still released normally.
"""

from __future__ import annotations

import gc

from PySide6 import QtCore


class GuiGarbageCollector(QtCore.QObject):
    def __init__(self, app):
        super().__init__(app)
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(5000)
        self.timer.timeout.connect(self.collect)
        self._was_enabled = gc.isenabled()
        gc.disable()
        self.timer.start()
        app.aboutToQuit.connect(self.restore)

    @QtCore.Slot()
    def collect(self):
        gc.collect()

    @QtCore.Slot()
    def restore(self):
        self.timer.stop()
        if self._was_enabled:
            gc.enable()


def ensure_gui_collection():
    """Install once per Qt application, preserving deliberately disabled GC."""
    app = QtCore.QCoreApplication.instance()
    if app is None or QtCore.QThread.currentThread() != app.thread():
        return None
    collector = getattr(app, "_pygcgc_garbage_collector", None)
    if collector is None and gc.isenabled():
        collector = GuiGarbageCollector(app)
        app._pygcgc_garbage_collector = collector
    return collector
