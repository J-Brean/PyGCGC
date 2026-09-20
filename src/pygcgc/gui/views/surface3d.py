"""A 3D surface view (pyqtgraph OpenGL), for looking at a small array --
a co-eluting blob's own footprint, in the ✂️ Overlap dialog -- from an
angle a flat 2D image can't show, especially two peaks close enough in
height and position to be hard to tell apart from directly above.

Needs PyOpenGL (the `gui` extra), which `pyqtgraph.opengl` imports at its
own import time, not lazily -- so importing it here is wrapped in a
plain `try`/`except ImportError`, and the widget falls back to an
explanatory label rather than failing to open the dialog it lives in.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets


class Surface3DView(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        self._gl = None
        self._item = None
        self.view = None
        try:
            import pyqtgraph.opengl as gl
            self._gl = gl
            self.view = gl.GLViewWidget()
            lay.addWidget(self.view)
        except ImportError:
            msg = QtWidgets.QLabel('3D view needs PyOpenGL:  uv pip install -e ".[gui]"')
            msg.setWordWrap(True)
            msg.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            lay.addWidget(msg)

    @property
    def available(self) -> bool:
        return self.view is not None

    def set_background(self, colour) -> None:
        if self.view is not None:
            self.view.setBackgroundColor(colour)

    def set_data(self, z: np.ndarray) -> None:
        """Show a 2D array as a height surface, coloured by height with
        the app's own magma colour map. Height is rescaled so the tallest
        point always stands a sensible amount above the (row, col) extent,
        regardless of the data's own units."""
        if self.view is None:
            return
        z = np.asarray(z, dtype=float)
        if self._item is not None:
            self.view.removeItem(self._item)
            self._item = None
        if z.size == 0:
            return
        zmax = float(z.max())
        colours = _height_colours(z, zmax)
        scaled = z * (max(z.shape) * 0.6 / zmax) if zmax > 0 else z
        self._item = self._gl.GLSurfacePlotItem(z=scaled, colors=colours, shader="shaded", smooth=False)
        self.view.addItem(self._item)
        self.view.opts["distance"] = max(z.shape) * 2.2
        self.view.update()


def _height_colours(z: np.ndarray, zmax: float) -> np.ndarray:
    cm = pg.colormap.get("magma")
    t = np.clip(z / zmax, 0.0, 1.0) if zmax > 0 else np.zeros_like(z)
    return cm.map(t.ravel(), mode="float").reshape(z.shape + (4,))
