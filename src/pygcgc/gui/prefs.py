"""Per-user window preferences (text size, window geometry).

Kept in a small INI file in the user's application-data folder, never in the
registry and never in a project file: a project travels between computers,
a screen size does not. `PYGCGC_SETTINGS` points at another file, which the
tests use so they never touch the real preferences.
"""

from __future__ import annotations

import os

from PySide6 import QtCore

DENSITIES = ("compact", "standard", "large")


def settings() -> QtCore.QSettings:
    path = os.environ.get("PYGCGC_SETTINGS")
    if path:
        return QtCore.QSettings(path, QtCore.QSettings.Format.IniFormat)
    return QtCore.QSettings(QtCore.QSettings.Format.IniFormat, QtCore.QSettings.Scope.UserScope, "PyGCGC", "PyGCGC")


def saved_density() -> str | None:
    """The text size the user chose, or None if they never did."""
    try:
        value = str(settings().value("ui/density", ""))
    except Exception:                                   # noqa: BLE001  (an unreadable file just means no preference)
        return None
    return value if value in DENSITIES else None


def save_density(name: str) -> None:
    if name not in DENSITIES:
        raise ValueError(f"unknown text size '{name}'")
    s = settings()
    s.setValue("ui/density", name)
    s.sync()


def default_density(screen=None) -> str:
    """Compact on a screen under 900 logical pixels high (a 1366 x 768 laptop, or a
    1080p screen scaled to 150 %), otherwise standard."""
    from PySide6 import QtGui
    screen = screen or QtGui.QGuiApplication.primaryScreen()
    if screen is None:
        return "standard"
    return "compact" if screen.availableGeometry().height() < 900 else "standard"


def density() -> str:
    return saved_density() or default_density()


def saved_geometry(key: str):
    try:
        value = settings().value(f"geometry/{key}")
    except Exception:                                   # noqa: BLE001
        return None
    return value if isinstance(value, QtCore.QByteArray) and not value.isEmpty() else None


def save_geometry(key: str, data: QtCore.QByteArray) -> None:
    s = settings()
    s.setValue(f"geometry/{key}", data)
    s.sync()
