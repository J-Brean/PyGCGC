"""Start the PyGCGC desktop app.

    pygcgc gui [FOLDER_OR_PROJECT]
    pygcgc-gui [FOLDER_OR_PROJECT]
"""

from __future__ import annotations

import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    try:
        from PySide6 import QtWidgets
    except ImportError:
        print("The GUI needs PySide6 and pyqtgraph:  uv pip install -e \".[gui]\"")
        return 1
    from PySide6 import QtCore, QtGui

    from .mainwindow import MainWindow
    from .theme import LIGHT
    from ..core.io import READERS

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    app.setApplicationName("PyGCGC")
    app.setOrganizationName("PyGCGC")

    # Splash: name in large type, author below, in the light theme's colours.
    art = QtGui.QPixmap(560, 280)
    art.fill(QtGui.QColor(LIGHT["window"]))
    p = QtGui.QPainter(art)
    p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
    p.setPen(QtGui.QColor(LIGHT["highlight"]))
    p.drawRect(0, 0, art.width() - 1, art.height() - 1)
    p.setPen(QtGui.QColor(LIGHT["text"]))
    p.setFont(QtGui.QFont("Georgia", 54, QtGui.QFont.Weight.Bold))
    p.drawText(QtCore.QRect(0, 60, 560, 100), QtCore.Qt.AlignmentFlag.AlignCenter, "PyGCGC")
    p.setPen(QtGui.QColor(LIGHT["mid"]))
    p.setFont(QtGui.QFont("Georgia", 16))
    p.drawText(QtCore.QRect(0, 165, 560, 40), QtCore.Qt.AlignmentFlag.AlignCenter, "James Brean")
    p.end()
    splash = QtWidgets.QSplashScreen(art)
    splash.show()
    app.processEvents()

    win = MainWindow()
    win.show()
    # Keep the splash up for a moment even when the window builds quickly.
    QtCore.QTimer.singleShot(1200, lambda: splash.finish(win))
    if argv:
        target = Path(argv[0])
        if target.is_dir():
            win.open_folder(target)
        elif target.suffix.lower() == ".json":
            win.open_project(target)
        elif target.suffix.lower() in READERS:
            win.project.add_files([target])
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
