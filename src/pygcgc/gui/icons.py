"""One-colour line icons, drawn in code so they follow the theme.

The workflow vocabulary uses emoji (🧹, 🔍, ...) in the docs, tooltips and menus. On screen,
colour emoji clash with a flat interface, so toolbar buttons and tabs show these line icons
in the text colour instead (`plain` strips the emoji from a label). Each icon is a list of
shapes on a 24 x 24 grid:

    ("M", x, y)                  start a new line
    ("L", x, y)                  straight segment
    ("C", x1, y1, x2, y2, x, y)  cubic curve
    ("Z",)                       close the shape
    ("O", cx, cy, rx, ry)        ellipse outline
    ("R", x, y, w, h, radius)    rounded rectangle outline
    ("F", cx, cy, r)             filled dot
"""

from __future__ import annotations

import re

from PySide6 import QtCore, QtGui

# Emoji, pictographs and variation selectors, but not the ✓ used for "done" marks.
_EMOJI = re.compile("[\U0001F000-\U0001FAFF☀-✒✔-➿⭐⭕⌀-⏿️‍]")


def plain(text: str) -> str:
    """The label without emoji, for toolbar buttons, tabs and banners."""
    return re.sub(r"\s{2,}", " ", _EMOJI.sub("", text)).strip()


ICONS: dict[str, list[tuple]] = {
    "open": [("M", 3, 6), ("L", 9, 6), ("L", 11, 9), ("L", 21, 9), ("L", 21, 19), ("L", 3, 19), ("Z",)],
    "fold": [("M", 3, 8), ("L", 9, 16), ("L", 15, 8), ("L", 21, 16)],
    "align": [("O", 12, 12, 5, 5), ("M", 12, 3), ("L", 12, 7), ("M", 12, 17), ("L", 12, 21),
              ("M", 3, 12), ("L", 7, 12), ("M", 17, 12), ("L", 21, 12)],
    "background": [("M", 3, 19), ("L", 21, 19), ("M", 3, 15), ("L", 8, 15), ("L", 12, 5), ("L", 16, 15), ("L", 21, 15)],
    "peaks": [("M", 3, 19), ("C", 8, 19, 8, 5, 12, 5), ("C", 16, 5, 16, 19, 21, 19)],
    "siloxane": [("O", 9, 14, 5, 5), ("O", 17, 7, 3, 3), ("O", 18, 17, 2, 2)],
    "hydro": [("M", 2, 15), ("L", 7, 9), ("L", 12, 15), ("L", 17, 9), ("L", 22, 15)],
    "classify": [("R", 4, 4, 7, 7, 2), ("R", 13, 4, 7, 7, 2), ("R", 4, 13, 7, 7, 2), ("R", 13, 13, 7, 7, 2)],
    "identify": [("O", 10, 10, 6, 6), ("M", 15, 15), ("L", 21, 21)],
    "quant": [("M", 5, 20), ("L", 5, 13), ("M", 10, 20), ("L", 10, 5), ("M", 15, 20), ("L", 15, 10),
              ("M", 20, 20), ("L", 20, 7)],
    "cluster": [("O", 7, 8, 3, 3), ("O", 17, 8, 3, 3), ("O", 12, 17, 3, 3)],
    "blank": [("M", 3, 19), ("C", 8, 19, 8, 5, 12, 5), ("C", 16, 5, 16, 19, 21, 19), ("M", 3, 12), ("L", 21, 12)],
    "istd": [("O", 12, 12, 8, 8), ("F", 12, 12, 2)],
    "overlap": [("O", 9, 12, 6, 6), ("O", 15, 12, 6, 6)],
    "batch": [("R", 4, 3, 16, 4, 1.5), ("R", 4, 10, 16, 4, 1.5), ("R", 4, 17, 16, 4, 1.5)],
    "export": [("M", 12, 4), ("L", 12, 15), ("M", 7, 9), ("L", 12, 4), ("L", 17, 9),
               ("M", 5, 16), ("L", 5, 20), ("L", 19, 20), ("L", 19, 16)],
    "outlines": [("O", 12, 12, 8, 6), ("F", 12, 12, 1.5)],
    "detected": [("O", 12, 12, 7, 7), ("M", 12, 2), ("L", 12, 5), ("M", 12, 19), ("L", 12, 22),
                 ("M", 2, 12), ("L", 5, 12), ("M", 19, 12), ("L", 22, 12)],
    "stencil": [("M", 12, 3), ("L", 20, 8), ("L", 20, 16), ("L", 12, 21), ("L", 4, 16), ("L", 4, 8), ("Z",)],
    "clean": [("M", 12, 3), ("C", 12, 9, 15, 12, 21, 12), ("C", 15, 12, 12, 15, 12, 21),
              ("C", 12, 15, 9, 12, 3, 12), ("C", 9, 12, 12, 9, 12, 3)],
    "lock": [("R", 6, 11, 12, 9, 2), ("M", 8, 11), ("L", 8, 8), ("C", 8, 3, 16, 3, 16, 8), ("L", 16, 11)],
    "auto": [("O", 12, 12, 8, 8), ("M", 12, 4), ("L", 12, 20)],
    "extend": [("M", 12, 4), ("L", 12, 20), ("M", 7, 9), ("L", 12, 4), ("L", 17, 9),
               ("M", 7, 15), ("L", 12, 20), ("L", 17, 15)],
    "region": [("R", 5, 5, 14, 14, 2)],
    "reset": [("M", 4, 9), ("L", 4, 4), ("L", 9, 4), ("M", 15, 4), ("L", 20, 4), ("L", 20, 9),
              ("M", 20, 15), ("L", 20, 20), ("L", 15, 20), ("M", 9, 20), ("L", 4, 20), ("L", 4, 15)],
    "prev": [("M", 14, 5), ("L", 7, 12), ("L", 14, 19)],
    "next": [("M", 10, 5), ("L", 17, 12), ("L", 10, 19)],
    "review": [("M", 9, 3), ("L", 15, 3), ("M", 10, 3), ("L", 10, 10), ("L", 5, 19), ("L", 19, 19),
               ("L", 14, 10), ("L", 14, 3)],
}


def make_icon(name: str, normal: str, active: str, disabled: str, size: int = 18) -> QtGui.QIcon:
    """A QIcon with Normal (off), Normal (on, for checked buttons) and Disabled pixmaps."""
    shapes = ICONS[name]
    icon = QtGui.QIcon()
    for colour, mode, state in ((normal, QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off),
                                (active, QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.On),
                                (disabled, QtGui.QIcon.Mode.Disabled, QtGui.QIcon.State.Off)):
        scale = 2          # draw at 2x, so the icon stays sharp on high-density screens
        pix = QtGui.QPixmap(size * scale, size * scale)
        pix.setDevicePixelRatio(scale)
        pix.fill(QtCore.Qt.GlobalColor.transparent)
        p = QtGui.QPainter(pix)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing)
        p.scale(size / 24, size / 24)
        pen = QtGui.QPen(QtGui.QColor(colour), 1.7)
        pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(QtCore.Qt.PenJoinStyle.RoundJoin)
        p.setPen(pen)
        p.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        path = QtGui.QPainterPath()
        for s in shapes:
            op = s[0]
            if op == "M":
                path.moveTo(s[1], s[2])
            elif op == "L":
                path.lineTo(s[1], s[2])
            elif op == "C":
                path.cubicTo(*s[1:7])
            elif op == "Z":
                path.closeSubpath()
            elif op == "O":
                path.addEllipse(QtCore.QPointF(s[1], s[2]), s[3], s[4])
            elif op == "R":
                path.addRoundedRect(QtCore.QRectF(s[1], s[2], s[3], s[4]), s[5], s[5])
            elif op == "F":
                p.setBrush(QtGui.QColor(colour))
                p.drawEllipse(QtCore.QPointF(s[1], s[2]), s[3], s[3])
                p.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        p.drawPath(path)
        p.end()
        icon.addPixmap(pix, mode, state)
    return icon
