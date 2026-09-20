"""The 1D views.

Top: the reconstructed 1D chromatogram (each modulation summed down t2),
which looks like an ordinary GC trace. Drag the shaded band to pick a t1
range; the same band appears on the 2D view, and the raw stream below
zooms to it.

Bottom: the raw detector stream, one point per scan. When zoomed in to a
few minutes, the modulation boundaries are drawn so you can see how each
peak is sliced (and spot wrap-around).
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from ..render import Rendered

BOUNDARY_PEN = pg.mkPen((126, 232, 250, 90), width=1)
BAND_BRUSH = pg.mkBrush(126, 232, 250, 40)
MAX_BOUNDARY_SPAN_MIN = 3.0


class OneDView(QtWidgets.QWidget):
    bandChanged = QtCore.Signal(float, float)     # t1 range (s)
    timeClicked = QtCore.Signal(float)            # absolute time (s)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.result: Rendered | None = None
        self.canvas = pg.GraphicsLayoutWidget()
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.canvas)

        self.recon = self.canvas.addPlot(row=0, col=0)
        self.recon.setLabel("bottom", "t₁ (min)")
        self.recon.setLabel("left", "reconstructed 1D")
        self.recon.setMenuEnabled(False)
        self.recon.hideButtons()
        self.recon_curve = self.recon.plot(pen=pg.mkPen((209, 149, 224), width=1.3))
        self.band = pg.LinearRegionItem(brush=BAND_BRUSH, pen=pg.mkPen((126, 232, 250, 120)))
        self.band.setVisible(False)
        self.recon.addItem(self.band)
        self.band.sigRegionChangeFinished.connect(self._band_moved)
        self.recon.scene().sigMouseClicked.connect(self._clicked)

        self.raw = self.canvas.addPlot(row=1, col=0)
        self.raw.setLabel("bottom", "acquisition time (min)")
        self.raw.setLabel("left", "raw stream")
        self.raw.setMenuEnabled(False)
        self.raw.hideButtons()
        self.raw.setDownsampling(auto=True, mode="peak")
        self.raw.setClipToView(True)
        self.raw.setAutoVisible(y=True)
        self.raw_curve = self.raw.plot(pen=pg.mkPen((126, 232, 250), width=1))
        self.boundaries = pg.PlotDataItem(pen=BOUNDARY_PEN, connect="pairs")
        self.raw.addItem(self.boundaries, ignoreBounds=True)
        self.raw.sigXRangeChanged.connect(self._update_boundaries)
        self._first = True

    def set_background(self, colour) -> None:
        self.canvas.setBackground(colour)

    def show_result(self, result: Rendered) -> None:
        self.result = result
        x = result.t1 / 60.0 + result.period / 120.0
        self.recon_curve.setData(x, result.profile)
        self.raw_curve.setData(result.times / 60.0, result.signal)
        label = result.expr.label() if not result.expr.is_ratio else \
            f"m/z {'+'.join(map(str, result.expr.numerator or ()))} (numerator)"
        self.raw.setTitle(f"{result.name}   {label}", size="9pt")
        if self._first:
            self._first = False
            self.recon.enableAutoRange()
            self.raw.enableAutoRange()
            if not self.band.isVisible():
                mid = x[np.nanargmax(result.profile)] if np.isfinite(result.profile).any() else x[len(x) // 2]
                self.set_band(mid * 60 - 60, mid * 60 + 60, emit=True)
        self._update_boundaries()

    def set_band(self, t1a: float, t1b: float, emit: bool = False) -> None:
        self.band.blockSignals(True)
        self.band.setRegion((t1a / 60.0, t1b / 60.0))
        self.band.blockSignals(False)
        self.band.setVisible(True)
        self.raw.setXRange(t1a / 60.0, t1b / 60.0, padding=0.02)
        self.raw.enableAutoRange(axis="y")
        if emit:
            self.bandChanged.emit(t1a, t1b)

    def band_range(self):
        if not self.band.isVisible():
            return None
        a, b = self.band.getRegion()
        return a * 60.0, b * 60.0

    def _band_moved(self) -> None:
        a, b = self.band.getRegion()
        self.raw.setXRange(a, b, padding=0.02)
        self.raw.enableAutoRange(axis="y")
        self.bandChanged.emit(a * 60.0, b * 60.0)

    def _clicked(self, ev) -> None:
        if self.result is None or not ev.double():
            return
        for plot in (self.recon, self.raw):
            if plot.sceneBoundingRect().contains(ev.scenePos()):
                x = plot.vb.mapSceneToView(ev.scenePos()).x()
                self.timeClicked.emit(float(x * 60.0))

    def _update_boundaries(self) -> None:
        r = self.result
        if r is None:
            return
        (x0, x1), (y0, y1) = self.raw.vb.viewRange()
        if x1 - x0 > MAX_BOUNDARY_SPAN_MIN:
            self.boundaries.setData([], [])
            return
        starts = r.t1 / 60.0
        sel = starts[(starts >= x0) & (starts <= x1)]
        xs = np.repeat(sel, 2)
        lo, hi = float(np.nanmin(r.signal)), float(np.nanmax(r.signal))
        ys = np.tile([lo, hi], sel.size)
        self.boundaries.setData(xs, ys)
