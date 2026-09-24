"""The 2D chromatogram view.

Layout (one pyqtgraph canvas):

    +-------------------------------+---------+
    | 2D image (t1 across, t2 up)   | t2      |
    |            [colour key, top right]      |
    +-------------------------------+---------+
    | reconstructed 1D (t1)         |
    +-------------------------------+

The slice follows the cursor column and the 1D trace shares the t1 axis.
Clicking a pixel asks for its spectrum; with the region tool on, a
draggable box gives the summed spectrum of everything inside it.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from ..render import Rendered

SCALES = ("log", "sqrt", "linear")
SCALE_NAMES = {"log": "Log", "sqrt": "Square root", "linear": "Linear"}
DEFAULT_CMAP = "turbo"

CROSS_PEN = pg.mkPen((255, 255, 255, 110), width=1)
PERIOD_PEN = pg.mkPen((255, 255, 255, 170), width=1, style=QtCore.Qt.PenStyle.DashLine)
BAND_BRUSH = pg.mkBrush(126, 232, 250, 40)


def display_transform(img: np.ndarray, scale: str) -> np.ndarray:
    img = np.asarray(img, dtype=np.float64)
    if scale == "log":
        pos = img[np.isfinite(img) & (img > 0)]
        floor = max(np.percentile(pos, 1.0), 1e-12) if pos.size else 1.0
        out = np.log10(np.clip(img, floor, None))
        out[~np.isfinite(img)] = np.nan
        return out
    if scale == "sqrt":
        return np.sqrt(np.clip(img, 0, None))
    return img


def _merge_collinear(x: np.ndarray, y: np.ndarray, tol: float = 1e-7) -> tuple[np.ndarray, np.ndarray]:
    """Join touching horizontal or vertical (start, end, NaN) segments into single segments.

    A legacy blob outline is one segment per pixel edge, so a straight 30-row side is
    30 segments; the full run was ~250,000 of them, rebuilt into a painter path on
    every zoom step. Merging only ever joins segments that touch or overlap on the
    same line, so the drawing is unchanged. Other segments pass through untouched.
    Run this on display coordinates (after the wrap shift), never before it.
    """
    n = len(x) // 3 * 3
    x, y = np.asarray(x, dtype=float)[:n], np.asarray(y, dtype=float)[:n]
    if n == 0:
        return x, y
    x0, x1, y0, y1 = x[0::3], x[1::3], y[0::3], y[1::3]
    vertical = (x0 == x1) & (y0 != y1)
    horizontal = (y0 == y1) & (x0 != x1)
    other = ~(vertical | horizontal)

    def runs(line, a, b):
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        order = np.lexsort((lo, line))
        line, lo, hi = line[order], lo[order], hi[order]
        new = np.ones(line.size, dtype=bool)
        new[1:] = (line[1:] != line[:-1]) | (lo[1:] > hi[:-1] + tol)
        start = np.flatnonzero(new)
        return line[start], np.minimum.reduceat(lo, start), np.maximum.reduceat(hi, start)

    gap = np.nan
    xs, ys = [np.column_stack((x0[other], x1[other], np.full(other.sum(), gap)))], \
             [np.column_stack((y0[other], y1[other], np.full(other.sum(), gap)))]
    if vertical.any():
        at, lo, hi = runs(x0[vertical], y0[vertical], y1[vertical])
        xs.append(np.column_stack((at, at, np.full(at.size, gap))))
        ys.append(np.column_stack((lo, hi, np.full(at.size, gap))))
    if horizontal.any():
        at, lo, hi = runs(y0[horizontal], x0[horizontal], x1[horizontal])
        xs.append(np.column_stack((lo, hi, np.full(at.size, gap))))
        ys.append(np.column_stack((at, at, np.full(at.size, gap))))
    return np.concatenate(xs).ravel(), np.concatenate(ys).ravel()


def auto_levels(disp: np.ndarray, scale: str, processed: bool = False) -> tuple[float, float]:
    vals = disp[np.isfinite(disp)]
    if vals.size == 0:
        return 0.0, 1.0
    # Cleaned data are mostly noise around zero; start the colour scale above it.
    lo_pct = (90.0 if processed else 20.0) if scale == "log" else 1.0
    lo, hi = np.percentile(vals, [lo_pct, 99.9])
    if hi <= lo:
        hi = lo + 1.0
    return float(lo), float(hi)


CMAP_CHOICES = ("turbo", "viridis", "inferno", "magma", "plasma", "cividis", "CET-L1")


class ColourKey(pg.GraphicsObject):
    """A small colour key pinned over the corner of the image.

    The grey track is the range of the data. The coloured part is the range
    currently mapped onto the colour map: the left end is the lower limit and
    the right end the upper limit, both printed underneath. Drag either end to
    move that limit, drag the coloured part to slide both, double-click to type
    exact limits, right-click to change the colour map.
    """

    sigLevelsChanged = QtCore.Signal(float, float)
    sigEditRequested = QtCore.Signal()
    sigColormapRequested = QtCore.Signal(str)

    W, H, PAD, BAR_Y, BAR_H = 150.0, 46.0, 8.0, 19.0, 10.0

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAcceptedMouseButtons(QtCore.Qt.MouseButton.LeftButton | QtCore.Qt.MouseButton.RightButton)
        self.lo, self.hi = 0.0, 1.0
        self.span = (0.0, 1.0)
        self.title = ""
        self._mode = None
        self._last_x = 0.0
        self.log = False    # levels are log10 values, but are labelled as the values themselves
        self.set_colormap(DEFAULT_CMAP)
        self.setToolTip("Drag an end to move that colour limit, drag the coloured part to slide both.\n"
                        "Double-click to type exact limits. Right-click to change the colour map.")

    def boundingRect(self) -> QtCore.QRectF:
        return QtCore.QRectF(0, 0, self.W, self.H)

    def set_colormap(self, name: str) -> None:
        lut = pg.colormap.get(name).getLookupTable(0.0, 1.0, 256, alpha=False)
        self._lut = np.asarray(lut, dtype=np.uint8)
        self.update()

    def set_state(self, lo: float, hi: float, span: tuple[float, float] | None = None) -> None:
        """Levels, and optionally the data range. The track always includes the levels."""
        self.lo, self.hi = float(lo), float(hi)
        if span is not None:
            self.span = (float(span[0]), float(span[1]))
        self.span = (min(self.span[0], self.lo), max(self.span[1], self.hi))
        if self.span[1] <= self.span[0]:
            self.span = (self.span[0], self.span[0] + 1.0)
        self.update()

    def _x(self, v: float) -> float:
        a, b = self.span
        return self.PAD + (v - a) / (b - a) * (self.W - 2 * self.PAD)

    def paint(self, p, *args) -> None:
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, False)
        p.setPen(pg.mkPen(200, 200, 210, 160))
        p.setBrush(pg.mkBrush(18, 18, 26, 215))
        p.drawRect(QtCore.QRectF(0.5, 0.5, self.W - 1, self.H - 1))
        font = p.font()
        font.setPointSizeF(7.5)
        p.setFont(font)
        p.setPen(pg.mkPen(220, 220, 230))
        p.drawText(QtCore.QRectF(self.PAD, 2, self.W - 2 * self.PAD, 14),
                   QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter, self.title)
        track = QtCore.QRectF(self._x(self.span[0]), self.BAR_Y, self._x(self.span[1]) - self._x(self.span[0]), self.BAR_H)
        p.fillRect(track, QtGui.QColor(70, 70, 80))
        x0, x1 = self._x(self.lo), self._x(self.hi)
        if x1 > x0:
            # Slices of the colour map: simple and avoids handing Qt a raw numpy buffer
            n = 64
            edges = np.linspace(x0, x1, n + 1)
            for i in range(n):
                r, g, b = self._lut[int((i + 0.5) / n * 255)]
                p.fillRect(QtCore.QRectF(edges[i], self.BAR_Y, edges[i + 1] - edges[i] + 0.6, self.BAR_H),
                           QtGui.QColor(int(r), int(g), int(b)))
        p.setPen(pg.mkPen(255, 255, 255, width=2))
        for x in (x0, x1):
            p.drawLine(QtCore.QPointF(x, self.BAR_Y - 3), QtCore.QPointF(x, self.BAR_Y + self.BAR_H + 3))
        p.setPen(pg.mkPen(220, 220, 230))
        row = QtCore.QRectF(self.PAD, self.BAR_Y + self.BAR_H + 4, self.W - 2 * self.PAD, 13)
        if self.log:
            # Decade ticks along the track, so the axis reads as logarithmic
            p.setPen(pg.mkPen(200, 200, 210, 200))
            for d in range(int(np.ceil(self.span[0])), int(np.floor(self.span[1])) + 1):
                p.drawLine(QtCore.QPointF(self._x(d), self.BAR_Y + self.BAR_H),
                           QtCore.QPointF(self._x(d), self.BAR_Y + self.BAR_H + 3))
            p.setPen(pg.mkPen(220, 220, 230))
        lo_txt, hi_txt = (f"{10.0 ** self.lo:.3g}", f"{10.0 ** self.hi:.3g}") if self.log else (f"{self.lo:.3g}", f"{self.hi:.3g}")
        p.drawText(row, QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignVCenter, lo_txt)
        p.drawText(row, QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignVCenter, hi_txt)

    def mouseDragEvent(self, ev) -> None:
        if ev.button() != QtCore.Qt.MouseButton.LeftButton:
            ev.ignore()
            return
        ev.accept()
        if ev.isStart():
            x = ev.buttonDownPos().x()
            d_lo, d_hi = abs(x - self._x(self.lo)), abs(x - self._x(self.hi))
            if min(d_lo, d_hi) <= 7:
                self._mode = "lo" if d_lo <= d_hi else "hi"
            elif self._x(self.lo) < x < self._x(self.hi):
                self._mode = "both"
            else:
                self._mode = None
        if self._mode is None:
            return
        a, b = self.span
        # lastPos on the first event is stale in pyqtgraph, so track the previous x ourselves
        last = ev.buttonDownPos().x() if ev.isStart() else self._last_x
        self._last_x = ev.pos().x()
        dv = (ev.pos().x() - last) / (self.W - 2 * self.PAD) * (b - a)
        gap = 1e-3 * (b - a)
        lo, hi = self.lo, self.hi
        if self._mode in ("lo", "both"):
            lo += dv
        if self._mode in ("hi", "both"):
            hi += dv
        if hi - lo < gap:
            return
        self.set_state(lo, hi)
        self.sigLevelsChanged.emit(self.lo, self.hi)

    def hoverEvent(self, ev) -> None:
        # Claim clicks and drags from the image's ViewBox underneath
        if ev.isExit():
            return
        ev.acceptDrags(QtCore.Qt.MouseButton.LeftButton)
        ev.acceptClicks(QtCore.Qt.MouseButton.LeftButton)
        ev.acceptClicks(QtCore.Qt.MouseButton.RightButton)

    def mouseClickEvent(self, ev) -> None:
        ev.accept()
        if ev.button() == QtCore.Qt.MouseButton.LeftButton and ev.double():
            self.sigEditRequested.emit()
        elif ev.button() == QtCore.Qt.MouseButton.RightButton:
            menu = QtWidgets.QMenu()
            for name in CMAP_CHOICES:
                menu.addAction(name)
            act = menu.exec(QtGui.QCursor.pos())
            if act is not None:
                self.sigColormapRequested.emit(act.text())


class ChromatogramView(QtWidgets.QWidget):
    hovered = QtCore.Signal(float, float, float)      # t1 (s), t2 (s), value
    pixelClicked = QtCore.Signal(float)               # absolute time (s)
    pointClicked = QtCore.Signal(int, float)          # modulation index, t2 (s, may exceed the period)
    regionChanged = QtCore.Signal(float, float, float, float)   # t1a, t1b (s), t2a, t2b (s)
    bandChanged = QtCore.Signal(float, float)         # t1 range (s)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.result: Rendered | None = None
        self.scale = "log"
        self.levels_locked = False
        self._levels_key = None
        self._first_image = True
        self._follow_cursor = True
        self._slice_col = None

        pg.setConfigOptions(imageAxisOrder="row-major", antialias=True)
        self.canvas = pg.GraphicsLayoutWidget()
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.canvas)

        # Main image
        self.plot = self.canvas.addPlot(row=0, col=0)
        self.plot.setLabel("bottom", "t₁ (min)")
        self.plot.setLabel("left", "t₂ (s)")
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.image = pg.ImageItem()
        self.image.setColorMap(pg.colormap.get(DEFAULT_CMAP))
        self.plot.addItem(self.image)
        # 🧩 cluster colours: RGBA laid over the image, in the same place
        self.cluster_layer = pg.ImageItem()
        self.cluster_layer.setZValue(2)
        self.cluster_layer.setVisible(False)
        self.plot.addItem(self.cluster_layer)
        self._cluster_rgba = None
        self.vline = pg.InfiniteLine(angle=90, pen=CROSS_PEN)
        self.hline = pg.InfiniteLine(angle=0, pen=CROSS_PEN)
        self.period_line = pg.InfiniteLine(angle=0, pen=PERIOD_PEN)
        self.period_line.setVisible(False)
        for item in (self.vline, self.hline, self.period_line):
            self.plot.addItem(item, ignoreBounds=True)
        self.band = pg.LinearRegionItem(orientation="vertical", brush=BAND_BRUSH, movable=True,
                                        pen=pg.mkPen((126, 232, 250, 120)))
        self.band.setZValue(5)
        self.band.setVisible(False)
        self.plot.addItem(self.band, ignoreBounds=True)
        self.band.sigRegionChangeFinished.connect(self._band_moved)
        # Before/after divider for the 🧹 split view
        self.divider = pg.InfiniteLine(angle=90, movable=True,
                                       pen=pg.mkPen((255, 255, 255, 220), width=2))
        self.divider.setZValue(20)
        self.divider.setVisible(False)
        self.divider.sigPositionChanged.connect(self._compose)
        self.plot.addItem(self.divider, ignoreBounds=True)
        self.before_label = pg.TextItem("before", anchor=(1.1, 0), color=(255, 255, 255))
        self.after_label = pg.TextItem("after", anchor=(-0.1, 0), color=(255, 255, 255))
        for t in (self.before_label, self.after_label):
            t.setZValue(21)
            t.setVisible(False)
            self.plot.addItem(t, ignoreBounds=True)
        self._disp_after = None
        self._disp_before = None
        self.roi = pg.RectROI([0, 0], [1, 1], pen=pg.mkPen((126, 232, 250), width=1.5))
        self.roi.setVisible(False)
        self.plot.addItem(self.roi, ignoreBounds=True)
        self.roi.sigRegionChangeFinished.connect(self._roi_moved)

        # ⛓️ stencil drawing (an editable polygon) and the lifted-stencil overlay
        self.stencil_roi = pg.PolyLineROI([[0, 0], [1, 0], [1, 1], [0, 1]], closed=True,
                                          pen=pg.mkPen((94, 234, 212), width=1.5))
        self.stencil_roi.setVisible(False)
        self.plot.addItem(self.stencil_roi, ignoreBounds=True)
        self.stencil_overlay_on = False
        self._stencil_polys: dict = {}
        self._stencil_items: list = []

        # 🧪 hand-drawn integration footprint: an editable polygon, orange to tell it from a stencil
        self.footprint_roi = pg.PolyLineROI([[0, 0], [1, 0], [1, 1], [0, 1]], closed=True,
                                            pen=pg.mkPen((251, 146, 60), width=2))
        self.footprint_roi.setVisible(False)
        self.plot.addItem(self.footprint_roi, ignoreBounds=True)

        # t2 slice at the cursor column
        self.slice = self.canvas.addPlot(row=0, col=1)
        self.slice.setMaximumWidth(130)
        self.slice.setYLink(self.plot)
        self.slice.hideAxis("left")
        self.slice.setLabel("bottom", "slice")
        self.slice.getAxis("bottom").setStyle(showValues=False)
        self.slice.setMenuEnabled(False)
        self.slice.hideButtons()
        self.slice_before = self.slice.plot(pen=pg.mkPen((140, 140, 150), width=1,
                                                         style=QtCore.Qt.PenStyle.DotLine))
        self.slice_before.setVisible(False)
        self.slice_curve = self.slice.plot(pen=pg.mkPen((126, 232, 250), width=1.2))
        self.slice_hline = pg.InfiniteLine(angle=0, pen=CROSS_PEN)
        self.slice.addItem(self.slice_hline, ignoreBounds=True)

        # Colour bar: a small boxed key floating over the top-right corner of the image.
        # Drag its lines or the band between them to move the levels; double-click for exact limits.
        vb = self.plot.vb
        self.key = ColourKey()
        self.key.setParentItem(vb)
        self.key.setZValue(30)
        self.key.sigLevelsChanged.connect(self.image_levels_from_key)
        self.key.sigEditRequested.connect(self.edit_levels)
        self.key.sigColormapRequested.connect(self.set_colormap)
        vb.sigResized.connect(self._place_colourbar)

        # Reconstructed 1D profile
        self.profile = self.canvas.addPlot(row=1, col=0)
        self.profile.setMaximumHeight(130)
        self.profile.setXLink(self.plot)
        self.profile.setLabel("left", "Σ t₂")
        self.profile.setMenuEnabled(False)
        self.profile.hideButtons()
        self.profile_curve = self.profile.plot(pen=pg.mkPen((209, 149, 224), width=1.2))
        self.profile_vline = pg.InfiniteLine(angle=90, pen=CROSS_PEN)
        self.profile.addItem(self.profile_vline, ignoreBounds=True)
        self.canvas.ci.layout.setRowStretchFactor(0, 5)
        self.canvas.ci.layout.setRowStretchFactor(1, 1)
        self.canvas.ci.layout.setColumnStretchFactor(0, 8)
        self.canvas.ci.layout.setColumnStretchFactor(1, 1)

        self.plot.scene().sigMouseMoved.connect(self._mouse_moved)
        self.plot.scene().sigMouseClicked.connect(self._mouse_clicked)

        # 🧪 overlay: peak outlines and labelled markers for targets and internal standards
        self.overlay_on = True
        self._overlay = ([], (np.array([]), np.array([])))
        self.overlay_outline = pg.PlotDataItem(connect="finite", pen=pg.mkPen((126, 232, 250, 150), width=1))
        self.overlay_markers = pg.ScatterPlotItem(size=11, pxMode=True)
        self.overlay_labels: list[pg.TextItem] = []
        for item in (self.overlay_outline, self.overlay_markers):
            item.setZValue(15)
            self.plot.addItem(item, ignoreBounds=True)

        # 🏷️ overlay: internal standards labelled or matched on this run, coloured by confidence tier
        self._istd_overlay: list[dict] = []
        self.istd_markers = pg.ScatterPlotItem(size=12, pxMode=True)
        self.istd_labels: list[pg.TextItem] = []
        self.istd_markers.setZValue(17)
        self.plot.addItem(self.istd_markers, ignoreBounds=True)

        # 🔍 overlay: every blob found on this run, with the selected ones picked out
        self.blob_overlay_on = False
        self._blob_outline = (np.array([]), np.array([]))
        self._blob_highlight = (np.array([]), np.array([]))
        self._blob_siloxane = (np.array([]), np.array([]))
        self._candidate_regions = (np.array([]), np.array([]))
        self._candidate_highlight = (np.array([]), np.array([]))
        self._candidate_centres = (np.array([]), np.array([]))
        self._component_overlays = []
        self.component_items = []
        self._outline_cache: dict = {}
        self.blob_outline = pg.PlotDataItem(connect="finite", pen=pg.mkPen((250, 204, 21, 130), width=1))
        self.blob_highlight = pg.PlotDataItem(connect="finite", pen=pg.mkPen((250, 204, 21), width=2.2))
        # 🫧 blobs at or above the threshold: the ones removal would block out
        # Width 1: any wider pen uses Qt's slow wide-line stroker, and hundreds of bleed blobs can be flagged.
        self.blob_siloxane = pg.PlotDataItem(connect="finite", pen=pg.mkPen((255, 64, 160), width=1, style=QtCore.Qt.PenStyle.DashLine))
        self.candidate_regions = pg.PlotDataItem(
            connect="finite", pen=pg.mkPen((34, 211, 238, 150), width=1,
                                            style=QtCore.Qt.PenStyle.DashLine))
        self.candidate_highlight = pg.PlotDataItem(
            connect="finite", pen=pg.mkPen((34, 211, 238), width=2.2))
        self.candidate_centres = pg.ScatterPlotItem(
            size=7, symbol="x", pen=pg.mkPen((255, 255, 255), width=1.5), brush=None)
        for item in (self.blob_outline, self.blob_siloxane, self.blob_highlight,
                     self.candidate_regions, self.candidate_highlight, self.candidate_centres):
            item.setZValue(14)
            self.plot.addItem(item, ignoreBounds=True)

        self.empty = pg.TextItem("", anchor=(0.5, 0.5), color=(128, 121, 115))
        self.plot.addItem(self.empty, ignoreBounds=True)
        self.set_empty_message("Open a folder of LSC or CDF runs to start  (Ctrl+O)")

    # ----------------------------------------------------------- appearance
    def set_background(self, colour) -> None:
        self.canvas.setBackground(colour)

    def set_empty_message(self, text: str) -> None:
        self.empty.setText(text)
        self.empty.setVisible(bool(text))
        if text:
            self.empty.setPos(0.5, 0.5)
            self.plot.setRange(xRange=(0, 1), yRange=(0, 1), padding=0)

    def set_image_visible(self, on: bool) -> None:
        """Hide the base chromatogram image (raw or 🧹 cleaned) while leaving overlays on show."""
        self.image.setVisible(on)

    def set_image_opacity(self, alpha: float) -> None:
        self.image.setOpacity(alpha)

    # ------------------------------------------------------------- display
    def show_result(self, result: Rendered) -> None:
        self.set_empty_message("")
        levels_key = (result.expr.text(), result.expr.is_ratio, self.scale, result.snr)
        keep_levels = self.levels_locked and levels_key == self._levels_key and self.result is not None
        old_levels = self.image.getLevels() if keep_levels else None
        self.result = result
        disp = display_transform(result.image, self.scale)
        self._disp_after = disp
        self._disp_before = (display_transform(result.image_before, self.scale)
                             if result.image_before is not None else None)
        split = self._disp_before is not None
        if split and not self.divider.isVisible():
            xr = self.plot.vb.viewRange()[0]
            self.divider.blockSignals(True)
            self.divider.setPos(0.5 * (xr[0] + xr[1]))
            self.divider.blockSignals(False)
        for item in (self.divider, self.before_label, self.after_label):
            item.setVisible(split)
        self.image.setImage(self._composite(), autoLevels=False)
        x0 = result.t1[0] / 60.0
        w = result.image.shape[1] * result.period / 60.0
        h = result.image.shape[0] * result.pitch
        self.image.setRect(QtCore.QRectF(x0, -0.5 * result.pitch, w, h))
        self._place_cluster_layer()
        if old_levels is not None:
            self.image.setLevels(old_levels)
        else:
            self.image.setLevels(auto_levels(disp, self.scale, result.processed))
        finite = disp[np.isfinite(disp)]
        span = tuple(np.percentile(finite, [0.05, 99.99])) if finite.size else (0.0, 1.0)
        self.key.set_state(*self.image.getLevels(), span=span)
        self._levels_key = levels_key
        quantity = "S/N" if result.snr else ("ratio" if result.expr.is_ratio else "intensity")
        self.key.log = self.scale == "log"
        self.key.title = {"log": f"{quantity} (log)", "sqrt": f"√ {quantity}", "linear": quantity}[self.scale]
        self.key.update()
        self.period_line.setVisible(result.extended)
        self.period_line.setPos(result.period)
        state = ("🧹 preview" if result.preview else "🧹 cleaned") if result.processed else "raw"
        extra = "   S/N" if result.snr else ""
        extra += "   📍 aligned" if result.aligned else ""
        self.plot.setTitle(f"{result.name}   {result.expr.label()}{extra}   {state}   "
                           f"phase {result.phase:.3f} s", size="9pt")
        self._place_split_labels()
        self._draw_overlay()
        self._draw_blob_overlay()
        self._draw_istd_overlay()
        prof = result.profile
        self.profile_curve.setData(result.t1 / 60.0 + result.period / 120.0, prof)
        if self._first_image:
            self._first_image = False
            self.reset_view()
        self._update_slice(self._slice_col if self._slice_col is not None else int(np.nanargmax(prof)) if np.isfinite(prof).any() else 0)

    def set_cluster_layer(self, rgba) -> None:
        """RGBA uint8 (rows, modulations, 4), or None to remove it. It is laid over the image
        only while it fits: not on an aligned or extended image, whose pixels are elsewhere."""
        self._cluster_rgba = rgba
        self._place_cluster_layer()

    def _place_cluster_layer(self) -> None:
        r, rgba = self.result, self._cluster_rgba
        fits = (rgba is not None and r is not None and rgba.shape[:2] == r.image.shape
                and not r.aligned and not r.extended)
        self.cluster_layer.setVisible(bool(fits))
        if fits:
            self.cluster_layer.setImage(rgba, autoLevels=False, levels=(0, 255))
            self.cluster_layer.setRect(QtCore.QRectF(r.t1[0] / 60.0, -0.5 * r.pitch,
                                                     r.image.shape[1] * r.period / 60.0, r.image.shape[0] * r.pitch))

    def _composite(self) -> np.ndarray:
        after, before = self._disp_after, self._disp_before
        if before is None or self.result is None:
            return after
        r = self.result
        c = int(np.clip(np.floor((self.divider.value() * 60.0 - r.t1[0]) / r.period), 0, after.shape[1]))
        out = after.copy()
        out[:, :c] = before[:, :c]
        return out

    def _compose(self) -> None:
        if self._disp_before is None:
            return
        levels = self.image.getLevels()
        self.image.setImage(self._composite(), autoLevels=False)
        self.image.setLevels(levels)
        self._place_split_labels()

    def _place_split_labels(self) -> None:
        if not self.divider.isVisible():
            return
        y = self.plot.vb.viewRange()[1][1]
        x = self.divider.value()
        self.before_label.setPos(x, y)
        self.after_label.setPos(x, y)

    def split_position(self) -> float:
        return float(self.divider.value())

    # ------------------------------------------------------------- overlay
    def set_overlay(self, markers: list[dict], outline: tuple) -> None:
        """markers: dicts with t1, t2 (s), label, colour (RGB) and symbol.
        outline: (t1 array, t2 array) in seconds, NaN-separated."""
        self._overlay = (markers, outline)
        self._draw_overlay()

    def set_overlay_visible(self, on: bool) -> None:
        self.overlay_on = bool(on)
        self._draw_overlay()

    def _to_display(self, t1, t2):
        """Where a true (t1, t2) appears on the current image. Anything above
        the displayed height is drawn where it was recorded: one modulation
        later, one period lower."""
        r = self.result
        t1 = np.asarray(t1, dtype=float).copy()
        t2 = np.asarray(t2, dtype=float).copy()
        top = r.image.shape[0] * r.pitch - 0.5 * r.pitch
        for _ in range(3):
            over = np.isfinite(t2) & (t2 >= top)
            if not over.any():
                break
            t2[over] -= r.period
            t1[over] += r.period
        return t1 / 60.0, t2

    def _wrap_outline(self, ot1, ot2) -> tuple[np.ndarray, np.ndarray]:
        """(t1, t2) NaN-separated (start, end, NaN) edge triples, moved as
        whole edges onto the displayed (possibly extended) image, so an
        edge is never split across the wrap."""
        if not len(ot1):
            return np.array([]), np.array([])
        n = len(ot1) // 3 * 3
        seg1 = np.asarray(ot1[:n], float).reshape(-1, 3)
        seg2 = np.asarray(ot2[:n], float).reshape(-1, 3)
        mid1 = seg1[:, :2].mean(axis=1)
        mid2 = seg2[:, :2].mean(axis=1)
        x_mid, y_mid = self._to_display(mid1, mid2)
        dx = x_mid - mid1 / 60.0
        dy = y_mid - mid2
        x = seg1 / 60.0 + dx[:, None]
        y = seg2 + dy[:, None]
        x[:, 2] = np.nan
        y[:, 2] = np.nan
        return x.ravel(), y.ravel()

    def _display_outline(self, outline) -> tuple[np.ndarray, np.ndarray]:
        """`_wrap_outline` then `_merge_collinear`, cached per outline array and display geometry.

        Selecting a blob redraws every overlay, but the big outlines (all blobs, all
        component contours) are the same array objects each time, so they are only
        converted once.
        """
        ot1, ot2 = outline
        if not len(ot1):
            return np.array([]), np.array([])
        r = self.result
        key = (id(ot1), len(ot1), r.image.shape[0], r.pitch, r.period)
        hit = self._outline_cache.get(key)
        if hit is not None and hit[0] is ot1:
            return hit[1]
        drawn = _merge_collinear(*self._wrap_outline(ot1, ot2))
        if len(self._outline_cache) > 20000:
            self._outline_cache.clear()
        self._outline_cache[key] = (ot1, drawn)     # the reference keeps id(ot1) from being reused
        return drawn

    def _draw_overlay(self) -> None:
        for t in self.overlay_labels:
            self.plot.removeItem(t)
        self.overlay_labels = []
        markers, (ot1, ot2) = self._overlay
        if self.result is None or not self.overlay_on:
            self.overlay_outline.setData([], [])
            self.overlay_markers.setData([])
            return
        x, y = self._wrap_outline(ot1, ot2)
        self.overlay_outline.setData(x, y)
        spots = []
        for m in markers:
            x, y = self._to_display([m["t1"]], [m["t2"]])
            colour = m.get("colour", (255, 255, 255))
            spots.append({"pos": (float(x[0]), float(y[0])), "symbol": m.get("symbol", "o"),
                          "pen": pg.mkPen(colour, width=1.6), "brush": None, "data": m.get("key")})
            if m.get("label") and len(self.overlay_labels) < 150:
                t = pg.TextItem(m["label"], anchor=(0, 1), color=colour)
                t.setZValue(16)
                t.setPos(float(x[0]), float(y[0]))
                self.plot.addItem(t, ignoreBounds=True)
                self.overlay_labels.append(t)
        self.overlay_markers.setData(spots)

    # --------------------------------------------------------- 🏷️ internal-standard overlay
    TIER_COLOURS = {"high": (60, 200, 110), "medium": (230, 180, 40), "low": (230, 110, 70),
                    "not found": (150, 150, 150)}

    def set_istd_overlay(self, markers: list[dict]) -> None:
        """markers: dicts with t1, t2 (s), label and tier ("high"/"medium"/"low"/"not found")."""
        self._istd_overlay = markers
        self._draw_istd_overlay()

    def _draw_istd_overlay(self) -> None:
        for t in self.istd_labels:
            self.plot.removeItem(t)
        self.istd_labels = []
        if self.result is None or not self._istd_overlay:
            self.istd_markers.setData([])
            return
        spots = []
        for m in self._istd_overlay:
            x, y = self._to_display([m["t1"]], [m["t2"]])
            colour = self.TIER_COLOURS.get(m.get("tier"), (255, 255, 255))
            spots.append({"pos": (float(x[0]), float(y[0])), "symbol": "t", "size": 12,
                          "pen": pg.mkPen(colour, width=1.8), "brush": None})
            if m.get("label") and len(self.istd_labels) < 50:
                t = pg.TextItem(m["label"], anchor=(0.5, 1.4), color=colour)
                t.setZValue(18)
                t.setPos(float(x[0]), float(y[0]))
                self.plot.addItem(t, ignoreBounds=True)
                self.istd_labels.append(t)
        self.istd_markers.setData(spots)

    # --------------------------------------------------------- 🔍 blob overlay
    def set_blob_overlay(self, outline: tuple, highlight: tuple, siloxane: tuple | None = None) -> None:
        """outline: every blob's edges. highlight: the selected ones', drawn
        bolder on top. siloxane: the ones at or above the 🫧 threshold, dashed
        magenta. All (t1 array, t2 array), NaN-separated."""
        self._blob_outline, self._blob_highlight = outline, highlight
        self._blob_siloxane = siloxane if siloxane is not None else (np.array([]), np.array([]))
        self._candidate_regions = self._candidate_highlight = (np.array([]), np.array([]))
        self._candidate_centres = (np.array([]), np.array([]))
        self._component_overlays = []
        self._draw_blob_overlay()

    def set_candidate_overlay(self, regions: tuple, centres: tuple,
                              highlight: tuple | None = None) -> None:
        """Show MS-aware centres and dashed analysis rectangles.

        The rectangles deliberately do not imply pixel ownership: fitted
        component maps may overlap inside them.
        """
        self._candidate_regions = regions
        self._candidate_centres = centres
        self._candidate_highlight = (highlight if highlight is not None else
                                     (np.array([]), np.array([])))
        self._blob_outline = self._blob_highlight = (np.array([]), np.array([]))
        self._blob_siloxane = (np.array([]), np.array([]))
        self._draw_blob_overlay(components=False)     # contours are untouched by this call

    def set_component_overlay(self, overlays: list[dict]) -> None:
        """Set overlapping component contours.

        Each item contains ``outline`` as global-time arrays, ``colour`` as an
        RGB tuple, and optional ``width``/``style`` values.
        """
        self._component_overlays = list(overlays)
        self._draw_components()

    def set_blob_overlay_visible(self, on: bool) -> None:
        self.blob_overlay_on = bool(on)
        self._draw_blob_overlay()

    def _draw_blob_overlay(self, components: bool = True) -> None:
        if components:
            self._draw_components()
        if self.result is None or not self.blob_overlay_on:
            self.blob_outline.setData([], [])
            self.blob_siloxane.setData([], [])
            self.blob_highlight.setData([], [])
            self.candidate_regions.setData([], [])
            self.candidate_highlight.setData([], [])
            self.candidate_centres.setData([])
            return
        x, y = self._display_outline(self._blob_outline)
        self.blob_outline.setData(x, y)
        x, y = self._display_outline(self._blob_siloxane)
        self.blob_siloxane.setData(x, y)
        x, y = self._display_outline(self._blob_highlight)
        self.blob_highlight.setData(x, y)
        x, y = self._wrap_outline(*self._candidate_regions)
        self.candidate_regions.setData(x, y)
        x, y = self._wrap_outline(*self._candidate_highlight)
        self.candidate_highlight.setData(x, y)
        t1, t2 = self._candidate_centres
        if len(t1):
            x, y = self._to_display(t1, t2)
            self.candidate_centres.setData(x=x, y=y)
        else:
            self.candidate_centres.setData([])

    def _draw_components(self) -> None:
        """Draw the fitted-component contours as a handful of plot items.

        Adding and removing one pyqtgraph item per contour cost about five seconds
        for the ~1,600 contours of a full run. Contours that share a pen are joined
        (NaN-separated) into one item, so the item count is the number of distinct
        pens, and nothing is redrawn unless the overlays themselves change.
        """
        for item in self.component_items:
            self.plot.removeItem(item)
        self.component_items = []
        if self.result is None or not self.blob_overlay_on or not self._component_overlays:
            return
        pens: dict = {}
        for overlay in self._component_overlays:
            x, y = self._display_outline(overlay["outline"])
            if not x.size:
                continue
            key = (tuple(overlay.get("colour", (0, 220, 255))),
                   overlay.get("style", QtCore.Qt.PenStyle.SolidLine),
                   float(overlay.get("width", 1.0)))
            pens.setdefault(key, []).append((x, y))
        for (colour, style, width), outlines in pens.items():
            x, y = np.concatenate([o[0] for o in outlines]), np.concatenate([o[1] for o in outlines])
            item = pg.PlotDataItem(x, y, connect="finite",
                                   pen=pg.mkPen(colour, width=width, style=style))
            item.setZValue(18)
            self.plot.addItem(item, ignoreBounds=True)
            self.component_items.append(item)

    def centre_on(self, t1: float, t2: float, half_t1_s: float = 90.0, half_t2_s: float = 0.8) -> None:
        if self.result is None:
            return
        x, y = self._to_display([t1], [t2])
        self.plot.setRange(xRange=(x[0] - half_t1_s / 60.0, x[0] + half_t1_s / 60.0),
                           yRange=(max(y[0] - half_t2_s, 0), y[0] + half_t2_s), padding=0)

    def reset_view(self) -> None:
        r = self.result
        if r is None:
            return
        x0 = r.t1[0] / 60.0
        self.plot.setRange(xRange=(x0, x0 + r.image.shape[1] * r.period / 60.0),
                           yRange=(0, r.image.shape[0] * r.pitch), padding=0)
        self.profile.enableAutoRange(axis="y")

    def set_scale(self, scale: str) -> None:
        if scale not in SCALES:
            raise ValueError(scale)
        self.scale = scale
        if self.result is not None:
            self.show_result(self.result)

    def cycle_scale(self) -> str:
        self.set_scale(SCALES[(SCALES.index(self.scale) + 1) % len(SCALES)])
        return self.scale

    def auto_levels(self) -> None:
        if self.result is None:
            return
        disp = self._disp_after if self._disp_after is not None else self.image.image
        self.image.setLevels(auto_levels(disp, self.scale, self.result.processed))
        self.key.set_state(*self.image.getLevels())

    def pixel_value(self, row: int, col: int) -> float:
        """Value behind a pixel, respecting the before/after split."""
        r = self.result
        if r.image_before is not None:
            c = np.floor((self.divider.value() * 60.0 - r.t1[0]) / r.period)
            if col < c:
                return float(r.image_before[row, col])
        return float(r.image[row, col])

    def image_levels_from_key(self, lo: float, hi: float) -> None:
        self.image.setLevels((lo, hi))

    def _place_colourbar(self) -> None:
        # Pinned to the top-right corner of the image area
        self.key.setPos(self.plot.vb.width() - self.key.W - 6, 6)

    def _colourbar_hit(self, scene_pos) -> bool:
        return self.key.sceneBoundingRect().contains(scene_pos)

    def set_colormap(self, name: str) -> None:
        self.image.setColorMap(pg.colormap.get(name))
        self.key.set_colormap(name)

    def set_levels(self, lo: float, hi: float) -> None:
        """Set the colour limits, in displayed units (log10, sqrt or linear)."""
        self.image.setLevels((lo, hi))
        self.key.set_state(lo, hi)

    def edit_levels(self) -> None:
        """Ask for exact lower and upper colour limits (double-click on the colour bar)."""
        if self.result is None:
            return
        lo, hi = self.image.getLevels()
        is_log = self.scale == "log"
        # The log scale is an axis transform: the user types intensities (1e3), not exponents (3)
        shown = (10.0 ** lo, 10.0 ** hi) if is_log else (lo, hi)
        dlg = QtWidgets.QDialog(self)
        dlg.setWindowTitle("Colour limits")
        form = QtWidgets.QFormLayout(dlg)
        form.addRow(QtWidgets.QLabel(f"{self.key.title}" + (" (log scale)" if is_log else "")))
        edits = []
        for name, val in zip(("Lower", "Upper"), shown):
            ed = QtWidgets.QLineEdit(f"{val:.6g}")
            form.addRow(name, ed)
            edits.append(ed)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                             | QtWidgets.QDialogButtonBox.StandardButton.Cancel
                                             | QtWidgets.QDialogButtonBox.StandardButton.RestoreDefaults)
        buttons.button(QtWidgets.QDialogButtonBox.StandardButton.RestoreDefaults).setText("Reset to default")
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        reset = []
        def _reset():
            reset.append(True)
            dlg.accept()
        buttons.button(QtWidgets.QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(_reset)
        form.addRow(buttons)
        if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted:
            return
        if reset:
            self.auto_levels()
            return
        try:
            new_lo, new_hi = (float(e.text()) for e in edits)
        except ValueError:
            QtWidgets.QMessageBox.warning(self, "Colour limits", "Enter numbers, for example 1500 or 1e3.")
            return
        if is_log:
            if new_lo <= 0 or new_hi <= 0:
                QtWidgets.QMessageBox.warning(self, "Colour limits", "Limits must be positive on a log scale.")
                return
            new_lo, new_hi = np.log10(new_lo), np.log10(new_hi)
        if new_hi <= new_lo:
            QtWidgets.QMessageBox.warning(self, "Colour limits", "The upper limit must be above the lower limit.")
            return
        self.set_levels(new_lo, new_hi)

    # ------------------------------------------------------------ geometry
    def pixel_at(self, x_min: float, y_s: float):
        r = self.result
        if r is None:
            return None
        col = int(np.floor((x_min * 60.0 - r.t1[0]) / r.period))
        row = int(np.rint(y_s / r.pitch))
        if 0 <= col < r.image.shape[1] and 0 <= row < r.image.shape[0]:
            return row, col
        return None

    def absolute_time(self, x_min: float, y_s: float) -> float | None:
        px = self.pixel_at(x_min, y_s)
        if px is None:
            return None
        row, col = px
        r = self.result
        return float(r.t1[col] + row * r.pitch)

    # ------------------------------------------------------------- mouse
    def _mouse_moved(self, pos) -> None:
        if self.result is None or not self.plot.sceneBoundingRect().contains(pos):
            return
        p = self.plot.vb.mapSceneToView(pos)
        self.vline.setPos(p.x())
        self.hline.setPos(p.y())
        self.slice_hline.setPos(p.y())
        self.profile_vline.setPos(p.x())
        px = self.pixel_at(p.x(), p.y())
        if px is None:
            return
        row, col = px
        if self._follow_cursor:
            self._update_slice(col)
        r = self.result
        self.hovered.emit(float(r.t1[col]), float(row * r.pitch), self.pixel_value(row, col))

    def _mouse_clicked(self, ev) -> None:
        if self.result is None or ev.button() != QtCore.Qt.MouseButton.LeftButton:
            return
        if not self.plot.sceneBoundingRect().contains(ev.scenePos()) or self._colourbar_hit(ev.scenePos()):
            return
        p = self.plot.vb.mapSceneToView(ev.scenePos())
        t = self.absolute_time(p.x(), p.y())
        if t is None:
            return
        px = self.pixel_at(p.x(), p.y())
        self._slice_col = px[1]
        self._update_slice(px[1])
        # pixelClicked (the raw scan spectrum) fires first so that, when the click also
        # lands on a blob, pointClicked's own spectrum can overwrite it; a click that
        # misses every blob leaves the raw spectrum showing.
        self.pixelClicked.emit(t)
        self.pointClicked.emit(int(px[1]), float(px[0] * self.result.pitch))

    def _update_slice(self, col: int) -> None:
        r = self.result
        if r is None or not (0 <= col < r.image.shape[1]):
            return
        self.slice_curve.setData(np.nan_to_num(r.image[:, col]), r.t2)
        if r.image_before is not None:
            self.slice_before.setData(np.nan_to_num(r.image_before[:, col]), r.t2)
            self.slice_before.setVisible(True)
        else:
            self.slice_before.setVisible(False)

    # ------------------------------------------------------------- tools
    def set_region_tool(self, on: bool) -> None:
        r = self.result
        if on and r is not None:
            xr, yr = self.plot.vb.viewRange()
            w, h = (xr[1] - xr[0]) * 0.08, (yr[1] - yr[0]) * 0.15
            self.roi.setPos((xr[0] + (xr[1] - xr[0]) / 2 - w / 2, yr[0] + (yr[1] - yr[0]) / 2 - h / 2))
            self.roi.setSize((w, h))
        self.roi.setVisible(on and r is not None)
        if on and r is not None:
            self._roi_moved()

    def region(self):
        if not self.roi.isVisible():
            return None
        x, y = self.roi.pos()
        w, h = self.roi.size()
        return x * 60.0, (x + w) * 60.0, y, y + h

    def _roi_moved(self) -> None:
        reg = self.region()
        if reg is not None:
            self.regionChanged.emit(*reg)

    # --------------------------------------------------------- ⛓️ stencils
    def start_stencil_draw(self) -> None:
        """Show the stencil ROI as a small default rectangle near the
        centre of the current view, ready to drag into shape."""
        xr, yr = self.plot.vb.viewRange()
        cx, cy = (xr[0] + xr[1]) / 2, (yr[0] + yr[1]) / 2
        w, h = (xr[1] - xr[0]) * 0.06, (yr[1] - yr[0]) * 0.15
        pts = [[cx - w / 2, cy - h / 2], [cx + w / 2, cy - h / 2],
              [cx + w / 2, cy + h / 2], [cx - w / 2, cy + h / 2]]
        self.stencil_roi.setPoints(pts)
        self.stencil_roi.setPos([0, 0])
        self.stencil_roi.setVisible(True)

    def start_footprint_edit(self, vertices: list[tuple[float, float]]) -> None:
        """Show the integration polygon, (t1 in s, t2 in s), ready to drag."""
        self.footprint_roi.setPoints([[t1 / 60.0, t2] for t1, t2 in vertices])
        self.footprint_roi.setPos([0, 0])
        self.footprint_roi.setVisible(True)

    def footprint_vertices(self) -> list[tuple[float, float]]:
        out = []
        for h in self.footprint_roi.getHandles():
            p = self.footprint_roi.mapToView(h.pos())
            out.append((p.x() * 60.0, p.y()))
        return out

    def cancel_footprint_edit(self) -> None:
        self.footprint_roi.setVisible(False)

    def load_stencil_polygon(self, vertices: list[tuple[float, float]]) -> None:
        """Put an existing polygon, (t1 in s, t2 in s), into the editable
        ROI, to drag its corners into a new shape."""
        self.stencil_roi.setPoints([[t1 / 60.0, t2] for t1, t2 in vertices])
        self.stencil_roi.setPos([0, 0])
        self.stencil_roi.setVisible(True)

    def cancel_stencil_draw(self) -> None:
        self.stencil_roi.setVisible(False)

    def stencil_vertices(self) -> list[tuple[float, float]]:
        """The stencil ROI's current polygon, as (t1 in s, t2 in s)."""
        out = []
        for h in self.stencil_roi.getHandles():
            p = self.stencil_roi.mapToView(h.pos())
            out.append((p.x() * 60.0, p.y()))
        return out

    def set_stencil_overlay(self, polygons: dict) -> None:
        """polygons: name -> [(t1 in s, t2 in s), ...], already lifted
        onto the current run's own alkane series."""
        self._stencil_polys = polygons
        self._draw_stencil_overlay()

    def set_stencil_overlay_visible(self, on: bool) -> None:
        self.stencil_overlay_on = bool(on)
        self._draw_stencil_overlay()

    def _draw_stencil_overlay(self) -> None:
        for item in self._stencil_items:
            self.plot.removeItem(item)
        self._stencil_items = []
        if self.result is None or not self.stencil_overlay_on:
            return
        for name, poly in self._stencil_polys.items():
            if len(poly) < 3:
                continue
            t1s = [p[0] for p in poly] + [poly[0][0]]
            t2s = [p[1] for p in poly] + [poly[0][1]]
            x, y = self._to_display(t1s, t2s)
            outline = pg.PlotDataItem(x, y, pen=pg.mkPen((94, 234, 212), width=1.5))
            outline.setZValue(13)
            self.plot.addItem(outline, ignoreBounds=True)
            self._stencil_items.append(outline)
            label = pg.TextItem(name, anchor=(0, 1), color=(94, 234, 212))
            label.setZValue(13)
            label.setPos(float(x[0]), float(y[0]))
            self.plot.addItem(label, ignoreBounds=True)
            self._stencil_items.append(label)

    def set_band(self, t1a: float | None, t1b: float | None) -> None:
        if t1a is None:
            self.band.setVisible(False)
            return
        self.band.blockSignals(True)
        self.band.setRegion((t1a / 60.0, t1b / 60.0))
        self.band.blockSignals(False)
        self.band.setVisible(True)

    def _band_moved(self) -> None:
        a, b = self.band.getRegion()
        self.bandChanged.emit(a * 60.0, b * 60.0)

    def grab_png(self, path) -> None:
        self.canvas.grab().save(str(path))
