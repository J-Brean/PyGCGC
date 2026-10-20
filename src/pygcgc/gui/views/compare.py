"""⚖️ Image comparison: blink, grid, difference and overlay.

All four work on images already read onto run A's grid (through 📍 when
both runs are aligned, at the same retention times otherwise; see
`core/compare.py`), so the same pixel means the same (t1, t2) everywhere.
Nothing here refits or resamples; it only displays what it is given.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

MODES = ("blink", "grid", "difference", "overlay")
MODE_NAMES = {"blink": "Blink", "grid": "Grid", "difference": "Difference", "overlay": "Overlay"}
BLINK_MS = 650
GRID_MAX_TILES = 12


def _levels(img: np.ndarray, pct: float = 99.5) -> tuple[float, float]:
    v = img[np.isfinite(img) & (img > 0)]
    if v.size == 0:
        return 0.0, 1.0
    lo = np.percentile(v, 5.0)
    hi = np.percentile(v, pct)
    return (float(lo), float(hi)) if hi > lo else (0.0, float(v.max()) or 1.0)


def _norm(img: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.clip((np.sqrt(np.clip(img, 0, None)) - np.sqrt(max(lo, 0))) /
                   max(np.sqrt(hi) - np.sqrt(max(lo, 0)), 1e-12), 0, 1)


def robust_scale(ref: np.ndarray, run: np.ndarray) -> float:
    """A robust ref/run intensity scale from their brighter shared pixels,
    so a run that is simply more concentrated doesn't look like a wall of
    difference. The median ratio over the top decile of pixels, in place
    of a full Theil-Sen fit."""
    mask = (ref > np.percentile(ref, 90)) & (run > np.percentile(run, 90)) & (run > 0)
    if mask.sum() < 5:
        mask = run > 0
    if not mask.any():
        return 1.0
    ratio = ref[mask] / run[mask]
    ratio = ratio[np.isfinite(ratio)]
    return float(np.median(ratio)) if ratio.size else 1.0


class _ImagePane(QtWidgets.QWidget):
    """One pyqtgraph image with a title, used by blink/difference/overlay."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.canvas = pg.GraphicsLayoutWidget()
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.canvas)
        self.plot = self.canvas.addPlot()
        self.plot.setLabel("bottom", "t₁ (min)")
        self.plot.setLabel("left", "t₂ (s)")
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.image = pg.ImageItem()
        self.plot.addItem(self.image)
        self.hist = pg.HistogramLUTItem(image=self.image, orientation="vertical")
        self.hist.setMaximumWidth(90)
        self.canvas.addItem(self.hist)

    def set_background(self, colour) -> None:
        self.canvas.setBackground(colour)

    def show_mono(self, img: np.ndarray, t1, pitch: float, period: float, cmap: str,
                  levels: tuple[float, float] | None = None, diverging: bool = False) -> None:
        self.hist.setVisible(True)
        self.hist.gradient.setColorMap(pg.colormap.get(cmap))
        disp = img if diverging else np.log10(np.clip(img, max(np.nanmax(img), 1e-9) * 1e-4, None))
        self.image.setImage(disp, autoLevels=False)
        self._rect(t1, pitch, img.shape)
        if levels is None:
            levels = _levels(disp if not diverging else np.abs(disp), 99.5) if not diverging else \
                (-np.nanpercentile(np.abs(disp), 99), np.nanpercentile(np.abs(disp), 99))
        self.image.setLevels(levels)
        self.hist.setLevels(*levels)

    def show_rgb(self, rgb: np.ndarray, t1, pitch: float) -> None:
        self.hist.setVisible(False)
        # rgb is already float in [0, 1]; an explicit (0, 1) level range stops
        # whatever scalar levels a previous mono/diverging view left behind
        # (still stored on this same ImageItem) from being applied here too,
        # which otherwise washes an RGB image out to flat grey.
        self.image.setImage(rgb, autoLevels=False, levels=(0.0, 1.0))
        self._rect(t1, pitch, rgb.shape[:2])

    def _rect(self, t1, pitch, shape) -> None:
        x0 = t1[0] / 60.0
        w = shape[1] * (t1[1] - t1[0] if len(t1) > 1 else 4.0) / 60.0
        h = shape[0] * pitch
        self.image.setRect(QtCore.QRectF(x0, -0.5 * pitch, w, h))


class CompareView(QtWidgets.QWidget):
    """`image_source(index) -> (image, t1, pitch, period, label) | None` is
    how this view gets a run's image on A's grid without knowing anything
    about folds or alignment itself."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.image_source = None       # set by the compare window
        self.index_a = -1
        self.index_b = -1
        self.grid_indices: list[int] = []
        self._blink_show_ref = True

        self.mode = QtWidgets.QComboBox()
        for m in MODES:
            self.mode.addItem(MODE_NAMES[m], m)
        self.mode.currentIndexChanged.connect(self._mode_changed)
        self.info = QtWidgets.QLabel("Choose two runs to compare")
        self.info.setWordWrap(True)

        self.single = _ImagePane()
        self.blink_timer = QtCore.QTimer(self)
        self.blink_timer.setInterval(BLINK_MS)
        self.blink_timer.timeout.connect(self._blink_tick)
        self.blink_label = pg.TextItem("", anchor=(0, 1), color=(255, 255, 255))
        self.blink_label.setZValue(20)
        self.single.plot.addItem(self.blink_label, ignoreBounds=True)

        self.grid_canvas = pg.GraphicsLayoutWidget()
        self.grid_canvas.setVisible(False)

        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel("Mode"))
        top.addWidget(self.mode)
        top.addStretch(1)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addLayout(top)
        lay.addWidget(self.info)
        lay.addWidget(self.single, 1)
        lay.addWidget(self.grid_canvas, 1)

    def set_background(self, colour) -> None:
        self.single.set_background(colour)
        self.grid_canvas.setBackground(colour)

    # ------------------------------------------------------------- driving
    def set_runs(self, index_a: int, index_b: int, grid_indices: list[int] | None = None) -> None:
        """A (the grid everything is drawn on), B, and the runs the grid mode tiles."""
        self.index_a = index_a
        self.index_b = index_b
        self.grid_indices = list(grid_indices or [])
        self.refresh()

    def _mode_changed(self) -> None:
        self.blink_timer.stop()
        self.refresh()

    def current_mode(self) -> str:
        return self.mode.currentData()

    def refresh(self) -> None:
        mode = self.current_mode()
        self.single.setVisible(mode != "grid")
        self.grid_canvas.setVisible(mode == "grid")
        if self.image_source is None or self.index_a < 0:
            self.info.setText("Choose two runs to compare")
            self.blink_timer.stop()
            return
        if mode == "grid":
            self.blink_timer.stop()
            self._refresh_grid()
            return
        ref = self.image_source(self.index_a)
        cur = self.image_source(self.index_b)
        if ref is None or cur is None:
            missing = "run B" if ref is not None else "run A"
            self.info.setText(f"No image for {missing} (is it loaded?)")
            self.blink_timer.stop()
            return
        self.info.setText("")
        if mode == "blink":
            self._show_blink(ref, cur)
        elif mode == "difference":
            self._show_difference(ref, cur)
        else:
            self._show_overlay(ref, cur)

    # ------------------------------------------------------------- blink
    def _show_blink(self, ref, cur) -> None:
        self._blink_ref, self._blink_cur = ref, cur
        self._blink_show_ref = True
        self._blink_tick()
        if not self.blink_timer.isActive():
            self.blink_timer.start()

    def _blink_tick(self) -> None:
        img, t1, pitch, period, label = self._blink_ref if self._blink_show_ref else self._blink_cur
        # One colour scale for both frames, from A, or a change of amount would
        # be hidden by each frame rescaling itself.
        ref = self._blink_ref[0]
        levels = _levels(np.log10(np.clip(ref, max(np.nanmax(ref), 1e-9) * 1e-4, None)))
        self.single.show_mono(img, t1, pitch, period, "magma", levels=levels)
        self.blink_label.setText(label)
        self.blink_label.setPos(t1[0] / 60.0, img.shape[0] * pitch)
        self._blink_show_ref = not self._blink_show_ref

    # -------------------------------------------------------- difference
    def _show_difference(self, ref, cur) -> None:
        self.blink_timer.stop()
        ref_img, t1, pitch, period, ref_label = ref
        cur_img, *_, cur_label = cur
        k = robust_scale(ref_img, cur_img)
        diff = ref_img - k * cur_img
        self.single.show_mono(diff, t1, pitch, period, "CET-D1", diverging=True)
        self.info.setText(f"{ref_label} minus {cur_label} × {k:.3g} (a robust brightness match, "
                          "not a concentration ratio)")

    # ------------------------------------------------------------ overlay
    def _show_overlay(self, ref, cur) -> None:
        self.blink_timer.stop()
        ref_img, t1, pitch, period, ref_label = ref
        cur_img, *_, cur_label = cur
        lo_r, hi_r = _levels(ref_img)
        lo_c, hi_c = _levels(cur_img)
        r = _norm(ref_img, lo_r, hi_r)
        g = _norm(cur_img, lo_c, hi_c)
        rgb = np.zeros(r.shape + (3,), dtype=np.float32)
        rgb[..., 0] = r          # reference: magenta (red + blue)
        rgb[..., 2] = r
        rgb[..., 1] = g          # this run: green
        self.single.show_rgb(np.clip(rgb, 0, 1), t1, pitch)
        self.info.setText(f"{ref_label} in magenta, {cur_label} in green; shared peaks look white")

    # ----------------------------------------------------------------- grid
    def _refresh_grid(self) -> None:
        self.grid_canvas.clear()
        indices = list(dict.fromkeys([self.index_a, *self.grid_indices]))[:GRID_MAX_TILES]
        tiles = [(i, src) for i in indices if (src := self.image_source(i)) is not None]
        self.info.setText(f"{len(tiles)} run{'s' if len(tiles) != 1 else ''} on A's grid, same colour scale"
                          + (f" (the first {GRID_MAX_TILES})" if len(self.grid_indices) >= GRID_MAX_TILES else ""))
        peak = max((float(np.nanmax(src[0])) for _i, src in tiles), default=1.0)
        lo, hi = 0.0, np.sqrt(peak) or 1.0
        cols = int(np.ceil(np.sqrt(len(tiles))))
        for k, (i, src) in enumerate(tiles):
            r, c = divmod(k, cols)
            plot = self.grid_canvas.addPlot(row=r, col=c)
            plot.hideAxis("left")
            plot.hideAxis("bottom")
            plot.setMenuEnabled(False)
            plot.hideButtons()
            plot.setTitle(src[4] + (" (A)" if i == self.index_a else ""), size="8pt")
            img, t1, pitch, period, _ = src
            item = pg.ImageItem(np.sqrt(np.clip(img, 0, None)))
            item.setColorMap(pg.colormap.get("magma"))
            item.setLevels((lo, hi))
            item.setRect(QtCore.QRectF(t1[0] / 60.0, 0, img.shape[1] * period / 60.0, img.shape[0] * pitch))
            plot.addItem(item)
            if k == 0:
                first = plot
            else:
                plot.setXLink(first)
                plot.setYLink(first)
