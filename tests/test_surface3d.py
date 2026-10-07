"""The ✂️ dialog's 3D surface: its colours must be indexable by the mesh's faces."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
gl = pytest.importorskip("pyqtgraph.opengl")


def test_surface_colours_index_by_face():
    # pyqtgraph 0.14 raised IndexError on drawing when colours came as (rows, cols, 4).
    from pygcgc.gui.views.surface3d import surface_item
    z = np.outer(np.hanning(9), np.hanning(6)) * 1e5
    item = surface_item(gl, z)
    colours = item._meshdata.vertexColors(indexed="faces")
    assert colours.shape == (item._meshdata.faceCount(), 3, 4)
