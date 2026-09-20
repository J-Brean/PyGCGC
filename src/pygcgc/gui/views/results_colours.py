"""Status colours shared by the results table and review gallery."""

from PySide6 import QtGui

STATUS_HEX = {
    "ok": "#2E9E6B", "accepted": "#2E9E6B", "check": "#E0A020",
    "rejected": "#C2410C", "not found": "#8E86A0",
}
STATUS_QCOLOURS = {k: QtGui.QColor(v) for k, v in STATUS_HEX.items()}
