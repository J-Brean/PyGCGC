"""🧪 Chemical space: the named peaks of a run plotted by what the compounds are, not when they eluted.

Each peak that has a name *and a formula* is placed by its formula: H/C against O/C (a van Krevelen diagram),
H/C, O/C and double-bond equivalents against carbon number, retention index against carbon number, the Kendrick
mass defect, the 2D chromatogram coloured by chemistry, and the share of the run's signal by heteroatom family
and by carbon number. Dot size is the peak's share of the run's signal.

Peaks left out are counted in the note under the plot, with the reason: no name, a name that is no library's
(typed by the user, so no formula), or a formula the plots cannot read. Nothing is guessed. A name you typed
that matches a library compound keeps that compound's formula.

Clicking a dot selects that peak in the matcher. The arithmetic is `core/chemspace.py`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from ..core.chemspace import Space, space_of
from .peak_naming import identity_of
from .theme import fit_to_screen

OKABE = ["#009E73", "#0072B2", "#E69F00", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442"]
MARKERS = ["o", "s", "t", "d", "p", "h", "star"]
LEVEL_COLOURS = {"confident": "#009E73", "probable": "#0072B2", "tentative": "#E69F00", "unknown": "#777777",
                 "named": "#CC79A7"}
LEVEL_MARKERS = {"confident": "o", "probable": "s", "tentative": "t", "unknown": "x", "named": "d"}
LEVEL_WORDS = {"confident": "Confident", "probable": "Probable", "tentative": "Tentative", "unknown": "Unknown",
               "named": "Named by you"}
SHOW = (("all", "Every named peak with a formula"), ("good", "Confident and probable only"), ("named", "Named by you only"))
COLOURS = (("level", "Confidence"), ("family", "Heteroatom family"), ("oc", "O/C"), ("hc", "H/C"), ("c", "Carbon number"),
           ("dbe", "Double-bond equivalents"))
PLOTS = (("vk", "Van Krevelen: H/C against O/C"),
         ("hc_c", "H/C against carbon number"),
         ("oc_c", "O/C against carbon number"),
         ("dbe_c", "Double-bond equivalents against carbon number"),
         ("ri_c", "Retention index against carbon number"),
         ("kmd", "Kendrick mass defect (CH₂)"),
         ("chrom", "2D chromatogram coloured by chemistry"),
         ("family", "Signal by heteroatom family"),
         ("carbon", "Signal by carbon number"))
AXES = {"vk": ("oc", "hc"), "hc_c": ("c", "hc"), "oc_c": ("c", "oc"), "dbe_c": ("c", "dbe"), "ri_c": ("c", "ri100"),
        "kmd": ("kmass", "kmd"), "chrom": ("t1", "t2")}
LABELS = {"oc": "O/C (atomic)", "hc": "H/C (atomic)", "c": "Carbon number", "dbe": "Double-bond equivalents",
          "ri100": "Retention index / 100", "kmass": "Kendrick mass (CH₂ = 14)", "kmd": "Kendrick mass defect",
          "t1": "t₁ (min)", "t2": "t₂ (s)"}


@dataclass
class SpacePeak:
    label: str
    t1: float                   # minutes
    t2: float                   # seconds
    share: float
    level: str                  # confident / probable / tentative / unknown, or "named" for a name the user typed
    name: str
    formula: str
    space: Space
    ri: float

    def value(self, key: str) -> float:
        s = self.space
        return {"oc": s.oc, "hc": s.hc, "c": float(s.c), "dbe": s.dbe, "ri100": self.ri / 100.0, "kmass": s.kendrick_mass,
                "kmd": s.kmd, "t1": self.t1, "t2": self.t2}[key]


def space_rows(found, decisions: dict) -> tuple[list[SpacePeak], dict]:
    """The peaks that can be placed in chemical space, and a count of those that cannot, by reason."""
    skipped = {"no name": 0, "no formula (a name no library holds)": 0, "formula not readable": 0}
    rows: list[SpacePeak] = []
    if found is None:
        return rows, skipped
    total = sum(max(m.query.volume, 0.0) for m in found.matches)
    for m in found.matches:
        d = decisions.get(m.query.label)
        name, formula, _cas = identity_of(m, d)
        if not name:
            skipped["no name"] += 1
            continue
        if not formula:
            skipped["no formula (a name no library holds)"] += 1
            continue
        sp = space_of(formula)
        if sp is None:
            skipped["formula not readable"] += 1
            continue
        level = "named" if d and d.get("action") == "named" else m.level
        rows.append(SpacePeak(m.query.label, m.query.t1 / 60.0, m.query.t2, m.query.volume / total if total > 0 else float("nan"),
                              level, name, formula, sp, float(m.query.ri)))
    return rows, skipped


def _size(share: float, biggest: float) -> float:
    return 6.0 if not np.isfinite(share) or biggest <= 0 else 6.0 + 16.0 * float(np.sqrt(max(share, 0.0) / biggest))


class ChemSpaceWindow(QtWidgets.QDialog):
    peakClicked = QtCore.Signal(str)

    def __init__(self, matcher_window, parent=None):
        super().__init__(parent or matcher_window)
        self.setWindowTitle("🧪 Chemical space")
        self.setWindowFlag(QtCore.Qt.WindowType.WindowMaximizeButtonHint, True)
        fit_to_screen(self, 900, 720)
        self.matcher = matcher_window
        self.rows: list[SpacePeak] = []
        self.skipped: dict = {}
        self.selected = ""
        self._points: dict = {}               # label -> (x, y) on the plot now showing
        self._ring = None                     # the ring round the selected peak

        self.kind = QtWidgets.QComboBox()
        for key, text in PLOTS:
            self.kind.addItem(text, key)
        self.colour = QtWidgets.QComboBox()
        for key, text in COLOURS:
            self.colour.addItem(text, key)
        self.show_box = QtWidgets.QComboBox()
        for key, text in SHOW:
            self.show_box.addItem(text, key)
        for box in (self.kind, self.colour, self.show_box):
            box.currentIndexChanged.connect(self._draw)
        self.plot = pg.PlotWidget()
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.plot.setMinimumHeight(360)
        self.legend = self.plot.addLegend(offset=(-4, 4), labelTextSize="8pt")
        self.bar = QtWidgets.QLabel("")
        self.bar.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
        self.note = QtWidgets.QLabel("")
        self.note.setWordWrap(True)
        self.about = QtWidgets.QLabel("")
        self.about.setWordWrap(True)
        self.about.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.run_label = QtWidgets.QLabel("")

        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel("<b>Plot</b>"))
        top.addWidget(self.kind, 1)
        top.addWidget(QtWidgets.QLabel("Colour by"))
        top.addWidget(self.colour)
        top.addWidget(QtWidgets.QLabel("Show"))
        top.addWidget(self.show_box)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(self.run_label)
        lay.addLayout(top)
        lay.addWidget(self.plot, 1)
        lay.addWidget(self.bar)
        lay.addWidget(self.about)
        lay.addWidget(self.note)

    # ================================================================ data in
    def set_background(self, colour) -> None:
        self.plot.setBackground(colour)

    def refresh(self) -> None:
        """Rebuild from the matcher's current run and decisions."""
        m = self.matcher
        found = m._current_matches()
        index = m.run_box.currentData()
        self.run_label.setText(f"Run: <b>{m.project.entries[int(index)].name}</b>" if index is not None else "")
        self.rows, self.skipped = space_rows(found, m._decisions())
        self.selected = m._shown.query.label if m._shown is not None else self.selected
        self._draw()

    def set_selected(self, label: str) -> None:
        self.selected = label or ""
        self._draw_ring()

    # ================================================================ drawing
    def _visible(self) -> list[SpacePeak]:
        mode = self.show_box.currentData()
        if mode == "good":
            return [r for r in self.rows if r.level in ("confident", "probable")]
        if mode == "named":
            return [r for r in self.rows if r.level == "named"]
        return list(self.rows)

    def _colouring(self, rows: list[SpacePeak]):
        """`(groups, legend text)`: groups of `(name, colours, symbol, rows)`; continuous colour is one group per dot."""
        key = self.colour.currentData()
        if key == "level":
            groups = [(LEVEL_WORDS[lv], LEVEL_COLOURS[lv], LEVEL_MARKERS[lv], [r for r in rows if r.level == lv])
                      for lv in ("confident", "probable", "tentative", "unknown", "named")]
            return [g for g in groups if g[3]], ""
        if key == "family":
            totals: dict = {}
            for r in rows:
                totals[r.space.family] = totals.get(r.space.family, 0.0) + (r.share if np.isfinite(r.share) else 0.0)
            order = sorted(totals, key=lambda f: -totals[f])
            groups = []
            for i, fam in enumerate(order[:len(OKABE)]):
                groups.append((fam, OKABE[i], MARKERS[i], [r for r in rows if r.space.family == fam]))
            rest = [r for r in rows if r.space.family not in order[:len(OKABE)]]
            if rest:
                groups.append(("other", "#777777", "x", rest))
            return groups, ""
        values = np.array([r.value(key) for r in rows], dtype=float)
        lo, hi = (float(values.min()), float(values.max())) if values.size else (0.0, 1.0)
        span = (hi - lo) or 1.0
        cmap = pg.colormap.get("viridis")
        colours = cmap.map((values - lo) / span, mode="byte") if values.size else []
        brushes = [pg.mkBrush(*(int(c) for c in colours[i])) for i in range(len(rows))]
        return [("", brushes, "o", rows)], f"Colour: {dict(COLOURS)[key]}, {lo:.2f} (purple) to {hi:.2f} (yellow)"

    def _draw(self, *_a) -> None:
        kind = self.kind.currentData()
        self.plot.clear()
        self._ring = None
        self.plot.setLogMode(x=False, y=False)
        self.plot.getPlotItem().getAxis("bottom").setTicks(None)
        self.plot.getPlotItem().getAxis("left").setTicks(None)
        self.plot.enableAutoRange()
        self._points = {}
        rows = self._visible()
        skipped = "; ".join(f"{n} with {why}" for why, n in self.skipped.items() if n)
        self.note.setText(f"{len(rows)} of {len(self.rows) + sum(self.skipped.values())} peaks plotted"
                          + (f". Left out: {skipped}." if skipped else ".")
                          + (" Press Match peaks in the matcher first." if not self.rows and not skipped else ""))
        self.bar.setText("")
        self.about.setText(self._explain(kind))
        if not rows:
            return
        if kind == "family":
            self._draw_family(rows)
        elif kind == "carbon":
            self._draw_carbon(rows)
        else:
            self._draw_scatter(kind, rows)
        self._draw_ring()

    def _draw_scatter(self, kind: str, rows: list[SpacePeak]) -> None:
        xk, yk = AXES[kind]
        if kind == "ri_c":
            rows = [r for r in rows if np.isfinite(r.ri)]
            if not rows:
                self.note.setText(self.note.text() + " None of them has a retention index (no alkane ladder).")
                return
        self.plot.setLabel("bottom", LABELS[xk])
        self.plot.setLabel("left", LABELS[yk])
        biggest = max((r.share for r in rows if np.isfinite(r.share)), default=0.0)
        groups, colour_note = self._colouring(rows)
        for name, colour, symbol, members in groups:
            brush = colour if isinstance(colour, list) else pg.mkBrush(colour)
            item = pg.ScatterPlotItem(x=[r.value(xk) for r in members], y=[r.value(yk) for r in members],
                                      data=[r.label for r in members], size=[_size(r.share, biggest) for r in members],
                                      symbol=symbol, brush=brush, pen=pg.mkPen("k", width=0.6), name=name or None,
                                      tip=self._tip, hoverable=True)
            item.sigClicked.connect(self._clicked)
            self.plot.addItem(item)
        for r in rows:
            self._points[r.label] = (r.value(xk), r.value(yk))
        self._guides(kind, rows)
        self.bar.setText(colour_note)

    def _guides(self, kind: str, rows: list[SpacePeak]) -> None:
        grey = pg.mkPen((110, 110, 110), width=0.8, style=QtCore.Qt.PenStyle.DashLine)
        opts = dict(color=(90, 90, 90), movable=False)

        def line(pos, angle, text, at=0.97):
            self.plot.addItem(pg.InfiniteLine(pos=pos, angle=angle, pen=grey, label=text,
                                              labelOpts=dict(position=at, **opts)), ignoreBounds=True)
        if kind in ("vk", "hc_c"):
            line(2.0, 0, "H/C = 2: saturated hydrocarbons", 0.7 if kind == "vk" else 0.3)
            line(1.0, 0, "H/C = 1: benzene-like", 0.7 if kind == "vk" else 0.3)
        if kind == "ri_c":
            xs = [r.value("c") for r in rows]
            lo, hi = min(xs) - 1, max(xs) + 1
            diag = pg.PlotDataItem([lo, hi], [lo, hi], pen=grey)
            self.plot.addItem(diag, ignoreBounds=True)
            text = pg.TextItem("n-alkanes", color=(90, 90, 90), anchor=(1, 0))
            text.setPos(hi, hi)
            self.plot.addItem(text, ignoreBounds=True)

    def _family_totals(self, rows):
        totals: dict = {}
        for r in rows:
            totals[r.space.family] = totals.get(r.space.family, 0.0) + (r.share if np.isfinite(r.share) else 0.0)
        return sorted(totals.items(), key=lambda kv: -kv[1])

    def _draw_family(self, rows: list[SpacePeak]) -> None:
        self.plot.setLabel("bottom", "Heteroatom family (elements besides carbon and hydrogen)")
        self.plot.setLabel("left", "Share of the run's signal (%)")
        totals = self._family_totals(rows)
        self.plot.addItem(pg.BarGraphItem(x=np.arange(len(totals)), height=[100 * v for _f, v in totals], width=0.7,
                                          brush=pg.mkBrush("#0072B2")))
        self.plot.getPlotItem().getAxis("bottom").setTicks([[(i, f) for i, (f, _v) in enumerate(totals)]])

    def _draw_carbon(self, rows: list[SpacePeak]) -> None:
        self.plot.setLabel("bottom", "Carbon number")
        self.plot.setLabel("left", "Share of the run's signal (%)")
        carbons = sorted({r.space.c for r in rows})

        def signal(c, hydrocarbon):
            return 100 * sum(r.share for r in rows if r.space.c == c and (r.space.family == "CH") == hydrocarbon
                             and np.isfinite(r.share))
        hydro = np.array([signal(c, True) for c in carbons])
        other = np.array([signal(c, False) for c in carbons])
        self.plot.addItem(pg.BarGraphItem(x=carbons, height=hydro, width=0.8, brush=pg.mkBrush("#0072B2"),
                                          name="Hydrocarbons (C and H only)"))
        self.plot.addItem(pg.BarGraphItem(x=carbons, y0=hydro, height=other, width=0.8, brush=pg.mkBrush("#E69F00"),
                                          name="With other elements"))

    def _draw_ring(self) -> None:
        if self._ring is not None:
            self.plot.removeItem(self._ring)
            self._ring = None
        at = self._points.get(self.selected)
        if at is not None:
            self._ring = pg.ScatterPlotItem(x=[at[0]], y=[at[1]], symbol="o", size=22, brush=None,
                                            pen=pg.mkPen("#D55E00", width=2.4))
            self._ring.setZValue(20)
            self.plot.addItem(self._ring)

    def _tip(self, x, y, data) -> str:
        r = next((r for r in self.rows if r.label == data), None)
        if r is None:
            return ""
        s = r.space
        return (f"{r.name}\n{r.formula}   H/C {s.hc:.2f}   O/C {s.oc:.2f}   DBE {s.dbe:g}\n"
                f"{r.t1:.3f} min, {r.t2:.2f} s   {LEVEL_WORDS[r.level]}")

    def _clicked(self, _item, points, *_a) -> None:
        if len(points):
            label = str(points[0].data())
            self.selected = label
            self._draw_ring()
            self.peakClicked.emit(label)

    @staticmethod
    def _explain(kind: str) -> str:
        return {
            "vk": "Each dot is a named peak at its H/C and O/C. Alkanes sit at the top left (H/C = 2, no oxygen), aromatics lower, "
                  "and oxidised compounds (acids, alcohols, ketones) to the right.",
            "hc_c": "Saturation against size. Alkanes lie on H/C = 2; aromatics and PAHs fall as the rings grow.",
            "oc_c": "Oxidation against size. Hydrocarbons lie on the axis.",
            "dbe_c": "Rings plus double bonds against size. A benzene ring is 4, naphthalene 7, pyrene 12.",
            "ri_c": "Retention index/100 against carbon number. Compounds lie on the diagonal if they elute like the n-alkane of "
                    "their size; polar and aromatic compounds sit above it.",
            "kmd": "Homologous series (compounds differing by CH₂) share a Kendrick mass defect, so each series is a horizontal row.",
            "chrom": "The chromatogram with each named peak placed at its t₁ and t₂, coloured by chemistry: chemical families "
                     "occupy bands in the 2D plane.",
            "family": "Share of the run's signal by the elements in the compound's formula.",
            "carbon": "Share of the run's signal by carbon number: blue is hydrocarbons, orange compounds with other elements.",
        }[kind]
