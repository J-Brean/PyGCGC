"""📚 The plot explorer beside the peak table: small plots that help decide whether a match is right.

The first plot, and the default, is a miniature of the run's 2D chromatogram with the selected peak boxed, so
you can see at once whether it is a big isolated peak or a small one on a busy shoulder. The others answer the
questions that come next: does the retention index of this peak fit its neighbours, are the large peaks the
well-matched ones, how much of the run's signal is actually identified, and how do this peak's candidates
compare on spectrum and on retention index.

Every dot is one fitted peak. Colour *and* shape give the confidence (so the plots read in greyscale), and the
size follows the share of the run's signal. Clicking a dot selects that peak in the table.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

LEVEL_ORDER = ("confident", "probable", "tentative", "unknown")
# symbol, colour, word, plot marker (a different shape per level, so the plots survive greyscale)
LEVEL_STYLE = {"confident": ("●", "#009E73", "Confident", "o"), "probable": ("◐", "#0072B2", "Probable", "s"),
               "tentative": ("○", "#E69F00", "Tentative", "t"), "unknown": ("✕", "#777777", "Unknown", "x")}
SELECTED = "#D55E00"
PLOTS = (("map", "Where the peak sits in the run (2D map)"),
         ("zoom", "Close-up of the peak's neighbourhood (2D map)"),
         ("ri", "Retention index against t₂"),
         ("mf", "Match factor against signal"),
         ("pareto", "How much of the signal is identified"),
         ("candidates", "This peak's candidates: spectrum against retention index"))
ZOOM_T1_MIN = 1.0           # half-width of the close-up
ZOOM_T2_S = 0.8


@dataclass
class PeakDot:
    """The few numbers a plot needs from one matched peak."""
    label: str
    t1: float                   # minutes
    t2: float                   # seconds
    share: float                # fraction of the run's signal; NaN if unknown
    level: str
    mf: float = float("nan")
    ri: float = float("nan")
    name: str = ""
    shown: bool = True          # False for peaks the table filters hide: drawn faintly


@dataclass
class RunMap:
    """A TIC image for the 2D map: log-scaled, rows along t₂, columns along t₁."""
    image: np.ndarray
    x0: float                   # minutes at the left edge
    width: float                # minutes
    pitch: float                # seconds per row
    levels: tuple
    t2_max: float
    period: float               # seconds
    key: tuple = field(default_factory=tuple)


def tic_map(image: np.ndarray, t1_start: float, period: float, pitch: float, key: tuple = ()) -> RunMap:
    """Log-scale a folded TIC image the way the main image does and pick grey levels from it."""
    from .views.image2d import display_transform
    shown = display_transform(image, "log")
    finite = shown[np.isfinite(shown)]
    levels = tuple(np.percentile(finite, [5.0, 99.8])) if finite.size else (0.0, 1.0)
    return RunMap(shown, t1_start / 60.0, image.shape[1] * period / 60.0, pitch, levels, image.shape[0] * pitch,
                  period, key)


def _size(share: float, biggest: float) -> float:
    if not np.isfinite(share) or biggest <= 0:
        return 6.0
    return 5.0 + 15.0 * float(np.sqrt(max(share, 0.0) / biggest))


class MatchExplorer(QtWidgets.QWidget):
    peakClicked = QtCore.Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.peaks: list[PeakDot] = []
        self.run_map: RunMap | None = None
        self.map_source = None
        self.match = None                     # the selected PeakMatch, for the candidates plot
        self.selected = ""
        self.mf_lines: tuple = ()
        self._shown_map_key = None
        self._by_label: dict = {}

        self.kind = QtWidgets.QComboBox()
        for key, text in PLOTS:
            self.kind.addItem(text, key)
        self.kind.currentIndexChanged.connect(self._kind_changed)
        self.note = QtWidgets.QLabel("")
        self.note.setWordWrap(True)

        # The two map views share one plot; everything else is drawn on a second, cleared each time.
        self.map_plot = pg.PlotWidget()
        self.plot = pg.PlotWidget()
        for p in (self.map_plot, self.plot):
            p.setMenuEnabled(False)
            p.hideButtons()
            p.setMinimumHeight(190)
        self.map_plot.setLabel("bottom", "t₁ (min)")
        self.map_plot.setLabel("left", "t₂ (s)")
        self.image = pg.ImageItem()
        self.image.setZValue(-10)
        ink = pg.ColorMap([0.0, 1.0], [(250, 250, 250), (20, 20, 20)])
        self.image.setColorMap(ink)
        self.map_plot.addItem(self.image)
        self.map_legend = self.map_plot.addLegend(offset=(-4, 4), labelTextSize="8pt")
        self.dots: dict[str, pg.ScatterPlotItem] = {}
        for level in LEVEL_ORDER:
            _s, colour, word, marker = LEVEL_STYLE[level]
            item = pg.ScatterPlotItem(symbol=marker, brush=pg.mkBrush(colour), pen=pg.mkPen("k", width=0.6), name=word,
                                      tip=self._tip, hoverable=True)
            item.sigClicked.connect(self._dot_clicked)
            self.map_plot.addItem(item)
            self.dots[level] = item
        self.box = QtWidgets.QGraphicsRectItem()
        self.box.setPen(pg.mkPen(SELECTED, width=1.5))
        self.box.setBrush(pg.mkBrush(213, 94, 0, 50))
        self.box.setZValue(10)
        self.ring = pg.ScatterPlotItem(symbol="o", size=18, brush=None, pen=pg.mkPen(SELECTED, width=2.2))
        self.ring.setZValue(11)
        self.vline = pg.InfiniteLine(angle=90, pen=pg.mkPen(SELECTED, width=0.8, style=QtCore.Qt.PenStyle.DashLine))
        self.hline = pg.InfiniteLine(angle=0, pen=pg.mkPen(SELECTED, width=0.8, style=QtCore.Qt.PenStyle.DashLine))
        for item in (self.box, self.ring, self.vline, self.hline):
            self.map_plot.addItem(item, ignoreBounds=True)
        self.plot.addLegend(offset=(-4, 4), labelTextSize="8pt")

        self.stack = QtWidgets.QStackedWidget()
        self.stack.addWidget(self.map_plot)
        self.stack.addWidget(self.plot)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel("<b>Plot</b>"))
        top.addWidget(self.kind, 1)
        lay.addLayout(top)
        lay.addWidget(self.stack, 1)
        lay.addWidget(self.note)
        self._tip_names: dict = {}

    # ================================================================ appearance
    def set_background(self, colour) -> None:
        self.map_plot.setBackground(colour)
        self.plot.setBackground(colour)

    def current_kind(self) -> str:
        return self.kind.currentData()

    # ================================================================ data in
    def set_data(self, peaks: list[PeakDot], map_source=None, match=None, mf_lines: tuple = ()) -> None:
        """`map_source` is a function returning the run's `RunMap`, called only when a map is on show (the image
        takes a moment to build)."""
        self.peaks = list(peaks)
        self._by_label = {p.label: p for p in self.peaks}
        self.map_source = map_source
        self.match = match
        self.mf_lines = tuple(mf_lines)
        self._draw()

    def set_selected(self, label: str, match=None) -> None:
        self.selected, self.match = label or "", match
        if self.current_kind() in ("map", "zoom"):
            self._mark_selected()
        else:
            self._draw()

    # ================================================================ drawing
    def _kind_changed(self) -> None:
        self._draw()

    def _draw(self) -> None:
        kind = self.current_kind()
        if kind in ("map", "zoom"):
            self.stack.setCurrentWidget(self.map_plot)
            self._draw_map(zoom=kind == "zoom")
        else:
            self.stack.setCurrentWidget(self.plot)
            {"ri": self._draw_ri, "mf": self._draw_mf, "pareto": self._draw_pareto,
             "candidates": self._draw_candidates}[kind]()

    def _tip(self, x, y, data) -> str:
        p = self._by_label.get(data)
        if p is None:
            return ""
        share = f", {100 * p.share:.2f} % of the run" if np.isfinite(p.share) else ""
        return f"{p.name or '(no candidate)'}\n{p.t1:.3f} min, {p.t2:.2f} s{share}\n{LEVEL_STYLE[p.level][2]}"

    def _dot_clicked(self, _item, points, *_a) -> None:
        if len(points):
            self.peakClicked.emit(str(points[0].data()))

    def _fresh_plot(self, x_label: str, y_label: str):
        """Clear the second plot (its legend entries go with the items) and put its axes in order."""
        self.plot.clear()
        self.plot.setLogMode(x=False, y=False)
        self.plot.setLabel("bottom", x_label)
        self.plot.setLabel("left", y_label)
        self.plot.enableAutoRange()

    def _message(self, text: str) -> None:
        self.note.setText(text)

    # ---- the 2D map
    def _draw_map(self, zoom: bool) -> None:
        rm = self.run_map = self.map_source() if self.map_source is not None and self.peaks else None
        biggest = max((p.share for p in self.peaks if np.isfinite(p.share)), default=0.0)
        if rm is None:
            self.image.setVisible(False)
        else:
            if rm.key != self._shown_map_key:
                self.image.setVisible(True)
                self.image.setImage(rm.image, autoLevels=False, levels=rm.levels)
                self.image.setRect(QtCore.QRectF(rm.x0, -0.5 * rm.pitch, rm.width, rm.image.shape[0] * rm.pitch))
                self._shown_map_key = rm.key
                self.map_plot.setRange(xRange=(rm.x0, rm.x0 + rm.width), yRange=(0, rm.t2_max), padding=0.01)
            self.image.setVisible(True)
        for level, item in self.dots.items():
            mine = [p for p in self.peaks if p.level == level]
            if not mine:
                item.setData([], [])
                continue
            item.setData(x=[p.t1 for p in mine], y=[p.t2 for p in mine], data=[p.label for p in mine],
                         size=[_size(p.share, biggest) if p.shown else 4.0 for p in mine],
                         brush=[pg.mkBrush(LEVEL_STYLE[level][1]) if p.shown else pg.mkBrush(150, 150, 150, 70) for p in mine],
                         pen=[pg.mkPen("k", width=0.6) if p.shown else pg.mkPen(None) for p in mine])
        self._mark_selected(zoom=zoom)
        shown = sum(p.shown for p in self.peaks)
        self._message(("Each dot is a fitted peak: colour and shape give the confidence, size the share of the run's signal. "
                       f"{shown} of {len(self.peaks)} peaks pass the table's filters (the others are small grey dots). "
                       "Click a dot to select it. The orange box marks the selected peak, a few modulations wide.")
                      if self.peaks else "No matches yet. Press Match peaks.")

    def _mark_selected(self, zoom: bool | None = None) -> None:
        zoom = self.current_kind() == "zoom" if zoom is None else zoom
        p = self._by_label.get(self.selected)
        for item in (self.box, self.ring, self.vline, self.hline):
            item.setVisible(p is not None)
        if p is None:
            return
        period = self.run_map.period / 60.0 if self.run_map is not None else 0.07
        half_x, half_y = 2 * period, 0.35
        self.box.setRect(p.t1 - half_x, p.t2 - half_y, 2 * half_x, 2 * half_y)
        self.ring.setData(x=[p.t1], y=[p.t2])
        self.vline.setPos(p.t1)
        self.hline.setPos(p.t2)
        if zoom:
            self.map_plot.setRange(xRange=(p.t1 - ZOOM_T1_MIN, p.t1 + ZOOM_T1_MIN),
                                   yRange=(p.t2 - ZOOM_T2_S, p.t2 + ZOOM_T2_S), padding=0)

    # ---- plots on the second widget
    def _scatter_by_level(self, xs_of, ys_of) -> None:
        biggest = max((p.share for p in self.peaks if np.isfinite(p.share)), default=0.0)
        for level in LEVEL_ORDER:
            mine = [p for p in self.peaks if p.level == level and np.isfinite(xs_of(p)) and np.isfinite(ys_of(p))]
            _s, colour, word, marker = LEVEL_STYLE[level]
            item = pg.ScatterPlotItem(x=[xs_of(p) for p in mine], y=[ys_of(p) for p in mine], data=[p.label for p in mine],
                                      symbol=marker, size=[_size(p.share, biggest) if p.shown else 4.0 for p in mine],
                                      brush=[pg.mkBrush(colour) if p.shown else pg.mkBrush(150, 150, 150, 70) for p in mine],
                                      pen=pg.mkPen("k", width=0.6), name=word, tip=self._tip, hoverable=True)
            item.sigClicked.connect(self._dot_clicked)
            self.plot.addItem(item)

    def _ring_selected(self, x, y) -> None:
        if np.isfinite(x) and np.isfinite(y):
            ring = pg.ScatterPlotItem(x=[x], y=[y], symbol="o", size=20, brush=None, pen=pg.mkPen(SELECTED, width=2.2))
            ring.setZValue(20)
            self.plot.addItem(ring)

    def _draw_ri(self) -> None:
        self._fresh_plot("Retention index", "t₂ (s)")
        have = [p for p in self.peaks if np.isfinite(p.ri)]
        if not have:
            self._message("No peak has a retention index: there is no alkane ladder for this run, so matching rests on the spectrum alone.")
            return
        self._scatter_by_level(lambda p: p.ri, lambda p: p.t2)
        sel = self._by_label.get(self.selected)
        if sel is not None:
            self._ring_selected(sel.ri, sel.t2)
        self._message("Each dot is a fitted peak at its retention index (alkane scale, 100 per carbon) and second-dimension time. "
                      "Compounds of one chemical family line up along a band; a peak far off its neighbours' band is "
                      "more likely a wrong match. The ringed dot is the selected peak.")

    def _draw_mf(self) -> None:
        self._fresh_plot("Share of the run's signal (%)", "Match factor")
        self.plot.setLogMode(x=True, y=False)
        self._scatter_by_level(lambda p: np.log10(max(100 * p.share, 1e-4)) if np.isfinite(p.share) else np.nan,
                               lambda p: p.mf)
        for value, text in self.mf_lines:
            line = pg.InfiniteLine(pos=value, angle=0, pen=pg.mkPen((90, 90, 90), width=0.8, style=QtCore.Qt.PenStyle.DashLine),
                                   label=text, labelOpts=dict(position=0.04, color=(90, 90, 90), movable=False))
            self.plot.addItem(line)
        sel = self._by_label.get(self.selected)
        if sel is not None and np.isfinite(sel.share):
            self._ring_selected(np.log10(max(100 * sel.share, 1e-4)), sel.mf)
        self._message("Large peaks at the left of the 'good' lines (dashed: the confident, probable and tentative thresholds) are "
                      "the ones most worth checking by hand. Peaks with no candidate are not shown.")

    def _draw_pareto(self) -> None:
        self._fresh_plot("Peaks, largest first", "Cumulative share of the run's signal (%)")
        ranked = sorted((p for p in self.peaks if np.isfinite(p.share)), key=lambda p: p.share, reverse=True)
        if not ranked:
            self._message("No signal totals to show.")
            return
        cumulative = 100 * np.cumsum([p.share for p in ranked])
        x = np.arange(1, len(ranked) + 1, dtype=float)
        self.plot.addItem(pg.PlotDataItem(x, cumulative, pen=pg.mkPen((120, 120, 120), width=1.2)))
        for level in LEVEL_ORDER:
            idx = [i for i, p in enumerate(ranked) if p.level == level]
            _s, colour, word, marker = LEVEL_STYLE[level]
            item = pg.ScatterPlotItem(x=x[idx], y=cumulative[idx], data=[ranked[i].label for i in idx], symbol=marker, size=7,
                                      brush=pg.mkBrush(colour), pen=pg.mkPen("k", width=0.6), name=word, tip=self._tip,
                                      hoverable=True)
            item.sigClicked.connect(self._dot_clicked)
            self.plot.addItem(item)
        sel = self._by_label.get(self.selected)
        if sel in ranked:
            k = ranked.index(sel)
            self._ring_selected(x[k], cumulative[k])
        total = sum(p.share for p in ranked)
        good = sum(p.share for p in ranked if p.level in ("confident", "probable")) / total
        unknown = sum(p.share for p in ranked if p.level == "unknown") / total
        self._message(f"Confident or probable matches carry {100 * good:.0f} % of the run's signal; unknown peaks carry "
                      f"{100 * unknown:.0f} %. A curve that rises steeply through grey crosses means the biggest peaks are not named.")

    def _draw_candidates(self) -> None:
        self._fresh_plot("Library RI minus peak RI", "Match factor")
        match = self.match
        if match is None or not match.candidates:
            self._message("Select a peak that has candidates.")
            return
        q = match.query
        has_ri = np.isfinite(q.ri)
        if not has_ri:
            self.plot.setLabel("bottom", "Candidate rank")
        for in_gate, name, brush, pen in ((True, "Inside the RI gate", pg.mkBrush("#0072B2"), pg.mkPen("k", width=0.6)),
                                          (False, "Outside the RI gate", None, pg.mkPen(SELECTED, width=1.6))):
            xs, ys, labels = [], [], []
            for rank, c in enumerate(match.candidates, start=1):
                if bool(c.in_gate) != in_gate:
                    continue
                if has_ri:
                    if not np.isfinite(c.ri_lib):
                        continue
                    xs.append(c.ri_lib - q.ri)
                else:
                    xs.append(float(rank))
                ys.append(c.mf)
                labels.append(rank)
            item = pg.ScatterPlotItem(x=xs, y=ys, symbol="o", size=11, brush=brush, pen=pen, name=name)
            self.plot.addItem(item)
            for xv, yv, rank in zip(xs, ys, labels):
                text = pg.TextItem(str(rank), anchor=(0.5, 1.4), color=(60, 60, 60))
                text.setPos(xv, yv)
                self.plot.addItem(text)
        if has_ri:
            self.plot.addItem(pg.InfiniteLine(pos=0, angle=90, pen=pg.mkPen((90, 90, 90), width=0.8, style=QtCore.Qt.PenStyle.DashLine)))
        missing = sum(not np.isfinite(c.ri_lib) for c in match.candidates) if has_ri else 0
        self._message("Numbers are the candidates' ranks. The best candidates are high (good spectrum) and near the dashed line "
                      "(library RI close to the peak's)." + (f" {missing} candidate(s) with no library RI are not drawn." if missing else "")
                      if has_ri else "This peak has no retention index, so the candidates are ranked by the spectrum alone.")
