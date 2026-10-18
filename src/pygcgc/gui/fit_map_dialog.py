"""📊 Fit and residual map: what the fitted blobs explain, and what is left over.

Four linked views of one run, on the same axes as the main chromatogram:

- Observed: the cleaned total-ion image the blobs were fitted to.
- Fitted: the sum of every fitted blob (each blob's shape times its full spectrum), that is,
  the integrated area the blobs assign to each pixel.
- Residual: observed minus fitted. Red is signal no blob explains; blue is where the
  blobs put more signal than was measured. The colour scale is a signed logarithm
  scaled to the residual's own noise level, so both faint leftovers and large misfits
  stay visible: 0 is noise-sized, 1 is ten times the noise, 2 is a hundred times.

- Identified fraction: fitted divided by observed in each pixel, 0 where the blobs explain
  nothing and 1 where they explain all of it. Pixels within 3 noise levels of zero are blank.

Observed and fitted share one logarithmic colour scale, so they can be compared by eye.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from ..core.peaktable import format_accounting


def _diverging_map() -> pg.ColorMap:
    try:
        return pg.colormap.get("CET-D1")                     # blue - white - red
    except Exception:                                        # a build without the CET maps
        return pg.ColorMap([0.0, 0.5, 1.0], [(40, 80, 200), (235, 235, 235), (200, 40, 40)])


def signed_log(residual: np.ndarray, noise: float) -> np.ndarray:
    """sign(r) * log10(1 + |r| / noise): about linear near zero, logarithmic far from it."""
    return np.sign(residual) * np.log10(1.0 + np.abs(residual) / max(noise, 1e-12))


class FitMapWindow(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent, QtCore.Qt.WindowType.Window)
        self.setWindowTitle("📊 Fit and residual map")
        self.resize(1200, 900)
        self._data = None

        self.summary = QtWidgets.QLabel("")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.outlines_box = QtWidgets.QCheckBox("Blob outlines")
        self.outlines_box.setChecked(True)
        self.outlines_box.toggled.connect(self._draw_outlines)
        self.cursor = QtWidgets.QLabel("")
        key = QtWidgets.QLabel("Residual: red = signal no blob explains, blue = blobs exceed the data. "
                               "Scale: 0 is no residual, 1 is about 10× the noise, 2 about 100×.")
        key.setWordWrap(True)

        self.canvas = pg.GraphicsLayoutWidget()
        self.plots, self.images = [], []
        for row, title in enumerate(("Observed (cleaned)", "Fitted blobs (integrated area per pixel)",
                                     "Residual (observed − fitted)",
                                     "Identified fraction (fitted ÷ observed)")):
            plot = self.canvas.addPlot(row=row, col=0)
            plot.setTitle(title, size="9pt")
            plot.setLabel("left", "t₂ (s)")
            plot.setMenuEnabled(False)
            plot.hideButtons()
            if row == 3:
                plot.setLabel("bottom", "t₁ (min)")
            if row:
                plot.setXLink(self.plots[0])
                plot.setYLink(self.plots[0])
            image = pg.ImageItem(axisOrder="row-major")
            plot.addItem(image)
            self.plots.append(plot)
            self.images.append(image)
            plot.scene().sigMouseMoved.connect(lambda pos, p=plot: self._hover(p, pos))
        self.bars = [pg.ColorBarItem(values=(0, 1), colorMap=pg.colormap.get("turbo"), interactive=False,
                                     width=12, label="log₁₀ TIC"),
                     pg.ColorBarItem(values=(0, 1), colorMap=pg.colormap.get("turbo"), interactive=False,
                                     width=12, label="log₁₀ TIC"),
                     pg.ColorBarItem(values=(-1, 1), colorMap=_diverging_map(), interactive=False,
                                     width=12, label="signed log₁₀ (residual / noise)"),
                     pg.ColorBarItem(values=(0, 1), colorMap=pg.colormap.get("viridis"), interactive=False,
                                     width=12, label="fitted / observed")]
        for bar, image, plot in zip(self.bars, self.images, self.plots):
            bar.setImageItem(image, insert_in=plot)
        self._outline_items: list = []

        top = QtWidgets.QHBoxLayout()
        top.addWidget(self.outlines_box)
        top.addStretch(1)
        top.addWidget(self.cursor)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(self.summary)
        lay.addWidget(key)
        lay.addLayout(top)
        lay.addWidget(self.canvas, 1)

    # ------------------------------------------------------------------ data
    def set_data(self, run_name: str, tic: np.ndarray, model: np.ndarray, t1_start: float,
                 period: float, pitch: float, report: dict, outlines: list) -> None:
        """tic and model: (rows, modulations) on one modulation's rows; outlines: (t1 s, t2 s) triples."""
        self.setWindowTitle(f"📊 Fit and residual map: {run_name}")
        residual = tic - model
        centre = float(np.median(residual))
        noise = 1.4826 * float(np.median(np.abs(residual - centre)))
        noise = noise if noise > 0 else float(np.std(residual)) or 1.0
        # Share of each pixel's signal that the blobs account for. Pixels within 3 noise
        # levels of zero are left blank: their ratio is noise divided by noise.
        measured = tic > 3.0 * noise
        fraction = np.full(tic.shape, np.nan)
        fraction[measured] = np.clip(model[measured], 0.0, None) / tic[measured]
        self._data = dict(tic=tic, model=model, residual=residual, fraction=fraction, t1_start=t1_start,
                          period=period, pitch=pitch, outlines=outlines, noise=noise)
        positive = tic[tic > 0]
        floor = max(float(np.percentile(positive, 1.0)), 1e-12) if positive.size else 1.0
        observed_log = np.log10(np.clip(tic, floor, None))
        fitted_log = np.log10(np.clip(model, floor, None))
        lo, hi = np.percentile(observed_log, [90.0, 99.9])       # cleaned data sit mostly at noise
        hi = hi if hi > lo else lo + 1.0
        shown = signed_log(residual, noise)
        span = max(float(np.percentile(np.abs(shown), 99.5)), 0.5)
        rect = QtCore.QRectF(t1_start / 60.0, -0.5 * pitch, tic.shape[1] * period / 60.0, tic.shape[0] * pitch)
        for image, values, bar, levels in ((self.images[0], observed_log, self.bars[0], (lo, hi)),
                                           (self.images[1], fitted_log, self.bars[1], (lo, hi)),
                                           (self.images[2], shown, self.bars[2], (-span, span)),
                                           (self.images[3], fraction, self.bars[3], (0.0, 1.0))):
            image.setImage(values, autoLevels=False)
            image.setRect(rect)
            bar.setLevels(levels)
        self.summary.setText(format_accounting(report) + f" Residual noise level: {noise:,.3g} per pixel.")
        self._draw_outlines()
        self.plots[0].autoRange()

    def _draw_outlines(self) -> None:
        for plot, item in self._outline_items:
            plot.removeItem(item)
        self._outline_items = []
        if self._data is None or not self.outlines_box.isChecked() or not self._data["outlines"]:
            return
        d = self._data
        x = np.concatenate([o[0] for o in d["outlines"]])
        y = np.concatenate([o[1] for o in d["outlines"]])
        n = x.size // 3 * 3
        x, y = x[:n].reshape(-1, 3).copy(), y[:n].reshape(-1, 3).copy()
        # An edge above one modulation was recorded at the bottom of the next one.
        top = d["tic"].shape[0] * d["pitch"] - 0.5 * d["pitch"]
        wrapped = y[:, :2].mean(axis=1) >= top
        y[wrapped, :2] -= d["period"]
        x[wrapped, :2] += d["period"]
        x, y = x / 60.0, y
        x[:, 2] = y[:, 2] = np.nan
        for plot in self.plots:
            item = pg.PlotDataItem(x.ravel(), y.ravel(), connect="finite",
                                   pen=pg.mkPen((255, 255, 255, 150), width=1))   # 1 px: wider pens are slow
            item.setZValue(10)
            plot.addItem(item, ignoreBounds=True)
            self._outline_items.append((plot, item))

    def _hover(self, plot, pos) -> None:
        if self._data is None or not plot.sceneBoundingRect().contains(pos):
            return
        d = self._data
        p = plot.vb.mapSceneToView(pos)
        col = int(np.floor((p.x() * 60.0 - d["t1_start"]) / d["period"]))
        row = int(np.rint(p.y() / d["pitch"]))
        if not (0 <= row < d["tic"].shape[0] and 0 <= col < d["tic"].shape[1]):
            self.cursor.setText("")
            return
        share = d["fraction"][row, col]
        identified = f"   identified {share:.0%}" if np.isfinite(share) else ""
        self.cursor.setText(
            f"t₁ {p.x():.3f} min   t₂ {row * d['pitch']:.2f} s   observed {d['tic'][row, col]:,.4g}   "
            f"fitted {d['model'][row, col]:,.4g}   residual {d['residual'][row, col]:,.4g} "
            f"({d['residual'][row, col] / d['noise']:+.1f}× noise)" + identified)
