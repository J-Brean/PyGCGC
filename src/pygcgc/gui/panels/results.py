"""🧪 Results tables.

Four views of the same report: target results, all signals (every peak,
normalised to its nearest internal standard), internal standards across
the sequence, and the calibration summary. Double-click a row to jump to
that peak; with target rows selected, Y accepts and N rejects them.
"""

from __future__ import annotations

import math
from dataclasses import asdict

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from ...core.quant import calibration_table

from ..views.results_colours import STATUS_HEX as STATUS_COLOURS

VIEWS = {
    "targets": ("Targets", ["run", "compound", "status", "flags", "amount", "amount_se",
                            "concentration", "response", "t1", "t2", "cal_set", "note"]),
    "signals": ("All signals", ["run", "blob", "t1", "t2", "snr", "volume", "istd", "ratio",
                                "target", "surrogate", "amount", "concentration", "note"]),
    "istd": ("Internal standards", ["run", "compound", "found", "volume", "rel_to_median", "flag",
                                    "t1", "t2"]),
    "calibrations": ("Calibrations", ["cal_set", "compound", "response", "passes", "model", "weights",
                                      "r2", "worst_rel_error", "lod", "loq", "n", "mandel_p",
                                      "intercept_p", "excluded"]),
}

HEADERS = {
    "t1": "t₁ (min)", "t2": "t₂ (s)", "amount_se": "± SE", "concentration": "conc.",
    "rel_to_median": "vs median (%)", "worst_rel_error": "worst error (%)", "snr": "S/N",
    "cal_set": "cal set", "mandel_p": "Mandel p", "intercept_p": "intercept p", "r2": "R²",
}


def _fmt(key, v):
    if v is None:
        return ""
    if isinstance(v, float):
        if math.isnan(v):
            return ""
        if key == "t1":
            return f"{v / 60:.3f}"
        if key == "t2":
            return f"{v:.3f}"
        if key in ("rel_to_median", "worst_rel_error"):
            return f"{v:+.1f}" if key == "rel_to_median" else f"{v:.1f}"
        if key in ("r2",):
            return f"{v:.4f}"
        if key in ("mandel_p", "intercept_p"):
            return f"{v:.3g}"
        return f"{v:.4g}"
    if isinstance(v, bool):
        return "yes" if v else ""
    if isinstance(v, (list, tuple)):
        return "; ".join(map(str, v))
    return str(v)


class RowsModel(QtCore.QAbstractTableModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.rows: list[dict] = []
        self.keys: list[str] = []
        self.amount_unit = "calibration units"

    def set_rows(self, rows: list[dict], keys: list[str]) -> None:
        self.beginResetModel()
        self.rows, self.keys = rows, keys
        self.endResetModel()

    def rowCount(self, parent=QtCore.QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return 0 if parent.isValid() else len(self.keys)

    def headerData(self, section, orientation, role=QtCore.Qt.ItemDataRole.DisplayRole):
        if role == QtCore.Qt.ItemDataRole.DisplayRole and orientation == QtCore.Qt.Orientation.Horizontal:
            k = self.keys[section]
            if k in ('amount', 'amount_se', 'lod', 'loq'):
                return f"{HEADERS.get(k, k)} ({self.amount_unit})"
            if k == 'concentration':
                return f"Concentration ({self.amount_unit}/m³)"
            return HEADERS.get(k, k.replace("_", " "))
        return None

    def data(self, index, role=QtCore.Qt.ItemDataRole.DisplayRole):
        row = self.rows[index.row()]
        key = self.keys[index.column()]
        v = row.get(key)
        if role == QtCore.Qt.ItemDataRole.DisplayRole:
            return _fmt(key, v)
        if role == QtCore.Qt.ItemDataRole.UserRole:        # sort value
            if isinstance(v, float) and math.isnan(v):
                return -math.inf
            return v if isinstance(v, (int, float)) else str(v)
        if role == QtCore.Qt.ItemDataRole.ForegroundRole:
            if key == "status" and v in STATUS_COLOURS:
                return QtGui.QBrush(QtGui.QColor(STATUS_COLOURS[v]))
            if key in ("flag",) and v:
                return QtGui.QBrush(QtGui.QColor(STATUS_COLOURS["check"]))
            if key == "passes":
                return QtGui.QBrush(QtGui.QColor(STATUS_COLOURS["ok" if v else "rejected"]))
        if role == QtCore.Qt.ItemDataRole.ToolTipRole and key in ("flags", "note", "excluded"):
            return _fmt(key, v)
        return None


class SortProxy(QtCore.QSortFilterProxyModel):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setSortRole(QtCore.Qt.ItemDataRole.UserRole)
        self.run = ""
        self.compound = ""
        self.only_check = False

    def lessThan(self, left, right):
        a, b = left.data(self.sortRole()), right.data(self.sortRole())
        try:
            return a < b
        except TypeError:
            return str(a) < str(b)

    def filterAcceptsRow(self, source_row, source_parent):
        row = self.sourceModel().rows[source_row]
        if self.run and row.get("run") not in (None, self.run):
            return False
        if self.compound and self.compound not in (row.get("compound"), row.get("target"),
                                                   row.get("surrogate")):
            return False
        if self.only_check:
            st = row.get("status")
            if st is not None:
                return st in ("check", "not found")
            if "flag" in row:
                return bool(row["flag"])
            if "passes" in row:
                return not row["passes"]
        return True


class ResultsPanel(QtWidgets.QWidget):
    gotoRequested = QtCore.Signal(dict)               # a row
    decisionRequested = QtCore.Signal(list, str)       # [(run, compound)], accepted/rejected/clear
    exportRequested = QtCore.Signal()
    settingsRequested = QtCore.Signal()
    runRequested = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.report = None
        self.view = QtWidgets.QComboBox()
        for key, (label, _) in VIEWS.items():
            self.view.addItem(label, key)
        self.run_box = QtWidgets.QComboBox()
        self.compound_box = QtWidgets.QComboBox()
        self.check_box = QtWidgets.QCheckBox("Only to check")
        self.count = QtWidgets.QLabel("")
        self.banner = QtWidgets.QPushButton("Results are out of date: 🧪 Quantify again")
        self.banner.setStyleSheet("color: #E0A020; text-align: left;")
        self.banner.setFlat(True)
        self.banner.setVisible(False)
        self.banner.clicked.connect(self.runRequested)
        accept = QtWidgets.QPushButton("✅ Accept (Y)")
        reject = QtWidgets.QPushButton("Reject (N)")
        clear = QtWidgets.QPushButton("Clear")
        accept.clicked.connect(lambda: self._decide("accepted"))
        reject.clicked.connect(lambda: self._decide("rejected"))
        clear.clicked.connect(lambda: self._decide("clear"))
        self.decision_buttons = (accept, reject, clear)
        settings = QtWidgets.QPushButton("⚙️")
        settings.setToolTip("Quantification settings")
        settings.clicked.connect(self.settingsRequested)
        export = QtWidgets.QPushButton("Export CSV…")
        export.clicked.connect(self.exportRequested)

        self.model = RowsModel(self)
        self.proxy = SortProxy(self)
        self.proxy.setSourceModel(self.model)
        self.table = QtWidgets.QTableView()
        self.table.setModel(self.proxy)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.doubleClicked.connect(self._double)
        for key, status in (("Y", "accepted"), ("N", "rejected")):
            sc = QtGui.QShortcut(QtGui.QKeySequence(key), self.table)
            sc.setContext(QtCore.Qt.ShortcutContext.WidgetShortcut)
            sc.activated.connect(lambda s=status: self._decide(s))

        top = QtWidgets.QHBoxLayout()
        for w in (QtWidgets.QLabel("Show"), self.view, QtWidgets.QLabel("Run"), self.run_box,
                  QtWidgets.QLabel("Compound"), self.compound_box, self.check_box):
            top.addWidget(w)
        top.addStretch(1)
        top.addWidget(self.count)
        for b in self.decision_buttons:
            top.addWidget(b)
        top.addWidget(settings)
        top.addWidget(export)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addLayout(top)
        lay.addWidget(self.banner)
        # Internal standards across the sequence (shown with that view)
        self.chart = pg.PlotWidget()
        self.chart.setMenuEnabled(False)
        self.chart.hideButtons()
        self.chart.setMinimumHeight(130)
        self.chart.setLabel("left", "IS volume vs median (%)")
        self.chart.setLabel("bottom", "run in sequence order")
        self.chart.setVisible(False)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        split.addWidget(self.chart)
        split.addWidget(self.table)
        lay.addWidget(split)

        self.view.currentIndexChanged.connect(self.refresh)
        self.run_box.currentIndexChanged.connect(self._filters)
        self.compound_box.currentIndexChanged.connect(self._filters)
        self.check_box.toggled.connect(self._filters)

    # ---------------------------------------------------------------- data
    def set_report(self, report, stale: bool = False, running: bool = False) -> None:
        self.report = report
        if running:
            self.banner.setText("🧪 Quantifying…")
            self.banner.setVisible(True)
        else:
            self.banner.setText("Results are out of date: 🧪 Quantify again")
            self.banner.setVisible(stale)
        self._fill_boxes()
        self.refresh()

    def set_amount_unit(self, unit):
        self.model.amount_unit = unit
        if self.model.columnCount():
            self.model.headerDataChanged.emit(QtCore.Qt.Orientation.Horizontal, 0, self.model.columnCount()-1)

    def _fill_boxes(self) -> None:
        r = self.report
        runs = [] if r is None else [x.name for x in r.seq if x.quant is not None]
        comps = [] if r is None else [t.compound for t in r.targets]
        for box, items in ((self.run_box, runs), (self.compound_box, comps)):
            cur = box.currentData()
            box.blockSignals(True)
            box.clear()
            box.addItem("All", "")
            for it in items:
                box.addItem(it, it)
            i = box.findData(cur)
            box.setCurrentIndex(max(i, 0))
            box.blockSignals(False)

    def current_view(self) -> str:
        return self.view.currentData()

    def rows_for(self, key: str) -> list[dict]:
        r = self.report
        if r is None:
            return []
        if key == "targets":
            return [asdict(x) for x in r.results]
        if key == "signals":
            return [asdict(x) for x in r.signals]
        if key == "istd":
            return [asdict(x) for x in r.istd]
        return calibration_table(r.fits)

    def refresh(self) -> None:
        key = self.current_view()
        rows = self.rows_for(key)
        self.model.set_rows(rows, VIEWS[key][1])
        for b in self.decision_buttons:
            b.setEnabled(key == "targets")
        self.chart.setVisible(key == "istd")
        if key == "istd":
            self._draw_istd_chart()
        self._filters()
        self.table.resizeColumnsToContents()

    def set_chart_background(self, colour) -> None:
        self.chart.setBackground(colour)

    def _draw_istd_chart(self) -> None:
        """One line per internal standard: its volume relative to its own
        median across the sequence, with the flagging limits dashed."""
        self.chart.clear()
        r = self.report
        if r is None or not r.istd:
            return
        runs = [x.name for x in r.seq if x.quant is not None]
        order = {n: i for i, n in enumerate(runs)}
        names = sorted({x.compound for x in r.istd})
        cmap = pg.colormap.get("CET-C6")
        for k, name in enumerate(names):
            rows = [x for x in r.istd if x.compound == name and x.run in order]
            xs = np.array([order[x.run] for x in rows], dtype=float)
            ys = np.array([x.rel_to_median if x.found else np.nan for x in rows], dtype=float)
            colour = cmap.map((k + 0.5) / len(names), mode="qcolor")
            self.chart.plot(xs, ys, pen=pg.mkPen(colour, width=1.4), symbol="o", symbolSize=5,
                            symbolBrush=colour, symbolPen=None, connect="finite")
            finite = np.isfinite(ys)
            if finite.any():
                label = pg.TextItem(name, color=colour, anchor=(0, 0.5))
                label.setPos(xs[finite][-1] + 0.1, ys[finite][-1])
                self.chart.addItem(label)
        tol = r.settings.is_tolerance
        for y in (-tol, tol):
            self.chart.addItem(pg.InfiniteLine(angle=0, pos=y, pen=pg.mkPen((142, 134, 160),
                                                                         style=QtCore.Qt.PenStyle.DashLine)))
        self.chart.addItem(pg.InfiniteLine(angle=0, pos=0, pen=pg.mkPen((142, 134, 160))))
        ticks = [(i, n) for i, n in enumerate(runs)]
        self.chart.getAxis("bottom").setTicks([ticks])

    def _filters(self) -> None:
        # Qt 6.10 replaced invalidateFilter() with begin/endFilterChange().
        new_api = hasattr(self.proxy, "beginFilterChange")
        if new_api:
            self.proxy.beginFilterChange()
        self.proxy.run = self.run_box.currentData() or ""
        self.proxy.compound = self.compound_box.currentData() or ""
        self.proxy.only_check = self.check_box.isChecked()
        if new_api:
            self.proxy.endFilterChange()
        else:
            self.proxy.invalidateFilter()
        self.count.setText(f"{self.proxy.rowCount()} of {self.model.rowCount()} rows")

    def select_target(self, run: str, compound: str) -> None:
        if self.current_view() != "targets":
            self.view.setCurrentIndex(self.view.findData("targets"))
        for i in range(self.proxy.rowCount()):
            row = self.model.rows[self.proxy.mapToSource(self.proxy.index(i, 0)).row()]
            if row.get("run") == run and row.get("compound") == compound:
                self.table.selectRow(i)
                self.table.scrollTo(self.proxy.index(i, 0))
                return

    def selected_rows(self) -> list[dict]:
        idx = {self.proxy.mapToSource(i).row() for i in self.table.selectionModel().selectedRows()}
        return [self.model.rows[i] for i in sorted(idx)]

    def _double(self, index) -> None:
        row = self.model.rows[self.proxy.mapToSource(index).row()]
        self.gotoRequested.emit(row)

    def _decide(self, status: str) -> None:
        if self.current_view() != "targets":
            return
        pairs = [(r["run"], r["compound"]) for r in self.selected_rows()]
        if pairs:
            self.decisionRequested.emit(pairs, status)
