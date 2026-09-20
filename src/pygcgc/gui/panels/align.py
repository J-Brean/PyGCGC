"""📍 Align panel.

The reference run, a table of anchors on whichever run is currently
selected (with residuals once there's a transform), a kind override, and
batch anchor-finding. Picking an anchor itself is a two-click flow driven
by the main window (click the peak on the reference, then the same peak
on the target run); this panel starts and reports on it.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

KIND_LABELS = [("Auto", "auto"), ("Shift", "shift"), ("Affine", "affine"), ("Smooth", "smooth")]


class AlignPanel(QtWidgets.QWidget):
    referenceChanged = QtCore.Signal(str)              # run name, or "" for none
    addAnchorRequested = QtCore.Signal()
    removeAnchorRequested = QtCore.Signal(str)          # anchor name, on the current run
    clearRunRequested = QtCore.Signal()
    kindChanged = QtCore.Signal(str)                    # "auto" / "shift" / "affine" / "smooth"
    propagateRequested = QtCore.Signal()
    showAlignedToggled = QtCore.Signal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.align = None
        self.current_run: str | None = None
        self._updating = False

        self.reference = QtWidgets.QComboBox()
        self.reference.setToolTip("Every other run is aligned to this one")
        self.reference.currentIndexChanged.connect(self._reference_changed)

        self.add_btn = QtWidgets.QPushButton("➕ Add anchor  (A)")
        self.add_btn.setToolTip("Click a peak on the reference, then the same peak here")
        self.add_btn.clicked.connect(self.addAnchorRequested)
        self.propagate_btn = QtWidgets.QPushButton("📍 Auto-find on all runs")
        self.propagate_btn.setToolTip("Search every other run for the anchors you've named so far")
        self.propagate_btn.clicked.connect(self.propagateRequested)
        self.clear_btn = QtWidgets.QPushButton("Clear anchors on this run")
        self.clear_btn.clicked.connect(self.clearRunRequested)

        self.kind = QtWidgets.QComboBox()
        for label, v in KIND_LABELS:
            self.kind.addItem(label, v)
        self.kind.currentIndexChanged.connect(
            lambda _i: (not self._updating) and self.kindChanged.emit(self.kind.currentData()))

        self.table = QtWidgets.QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["Name", "ref t1", "ref t2", "run t1", "run t2", "resid"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        remove_row = QtWidgets.QPushButton("Remove selected anchor")
        remove_row.clicked.connect(self._remove_selected)

        self.summary = QtWidgets.QLabel("No reference set yet")
        self.summary.setWordWrap(True)
        self.show_aligned = QtWidgets.QCheckBox("Show aligned in the 2D/1D view")
        self.show_aligned.toggled.connect(self.showAlignedToggled)

        self.resid_plot = pg.PlotWidget()
        self.resid_plot.setMenuEnabled(False)
        self.resid_plot.hideButtons()
        self.resid_plot.setLabel("bottom", "run t1 (min)")
        self.resid_plot.setLabel("left", "residual (s)")
        self.resid_plot.setMinimumHeight(120)
        self.resid_t1 = self.resid_plot.plot([], [], pen=None, symbol="o", symbolSize=7,
                                             symbolBrush=(209, 149, 224), name="t1")
        self.resid_t2 = self.resid_plot.plot([], [], pen=None, symbol="s", symbolSize=6,
                                             symbolBrush=(126, 232, 250), name="t2 (×10)")
        self.resid_plot.addLine(y=0, pen=pg.mkPen((142, 134, 160)))
        self.resid_plot.addLegend()

        form = QtWidgets.QFormLayout()
        form.addRow("Reference", self.reference)
        form.addRow("Transform", self.kind)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addLayout(form)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.add_btn)
        row.addWidget(self.propagate_btn)
        lay.addLayout(row)
        lay.addWidget(self.table, 1)
        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(remove_row)
        row2.addWidget(self.clear_btn)
        lay.addLayout(row2)
        lay.addWidget(self.summary)
        lay.addWidget(self.resid_plot)
        lay.addWidget(self.show_aligned)

    def _reference_changed(self, _i) -> None:
        if not self._updating:
            self.referenceChanged.emit(self.reference.currentData() or "")

    def set_background(self, colour) -> None:
        self.resid_plot.setBackground(colour)

    # ------------------------------------------------------------- filling
    def set_align(self, align) -> None:
        self.align = align

    def refresh_runs(self, entries) -> None:
        cur = self.reference.currentData()
        self._updating = True
        self.reference.clear()
        self.reference.addItem("None", "")
        for e in entries:
            if e.run is not None:
                self.reference.addItem(e.name, e.name)
        i = self.reference.findData(cur if cur else (self.align.reference if self.align else ""))
        self.reference.setCurrentIndex(max(i, 0))
        self._updating = False

    def set_current_run(self, name: str | None) -> None:
        self.current_run = name
        self.refresh()

    def refresh(self) -> None:
        a = self.align
        if a is None:
            return
        is_ref = self.current_run is not None and self.current_run == a.reference
        self.add_btn.setEnabled(bool(a.reference) and self.current_run is not None and not is_ref)
        self.clear_btn.setEnabled(bool(self.current_run and a.points.get(self.current_run)))
        self._updating = True
        self.kind.setCurrentIndex(max(self.kind.findData(a.kind.get(self.current_run, "auto")), 0))
        self._updating = False
        self.show_aligned.setEnabled(bool(a.reference))

        pts = a.points.get(self.current_run, []) if self.current_run else []
        tf = a.transforms.get(self.current_run) if self.current_run else None
        resid = {r["name"]: r for r in tf.residuals()} if tf is not None else {}
        self.table.setRowCount(len(pts))
        for i, p in enumerate(pts):
            r = resid.get(p.name)
            resid_text = f"{r['resid_t1']:+.2f}s / {r['resid_t2'] * 1e3:+.0f}ms" if r else ""
            values = (p.name, f"{p.ref_t1 / 60:.3f}", f"{p.ref_t2:.3f}",
                     f"{p.run_t1 / 60:.3f}", f"{p.run_t2:.3f}", resid_text)
            for c, v in enumerate(values):
                item = self.table.item(i, c)
                if item is None:
                    item = QtWidgets.QTableWidgetItem()
                    self.table.setItem(i, c, item)
                item.setText(v)

        if is_ref:
            self.summary.setText(f"{self.current_run} is the reference: shown as itself, never aligned.")
        elif not a.reference:
            self.summary.setText("Pick a reference run above, then ➕ Add anchor.")
        elif self.current_run in a.errors:
            self.summary.setText(f"{self.current_run}: {a.errors[self.current_run]}")
        elif tf is not None:
            self.summary.setText(f"{self.current_run}: {tf.summary()}")
        elif self.current_run:
            self.summary.setText(f"{self.current_run}: no anchors yet.")
        else:
            self.summary.setText("")

        if tf is not None and pts:
            rows = tf.residuals()
            x = np.array([r["run_t1"] for r in rows]) / 60.0
            self.resid_t1.setData(x, [r["resid_t1"] for r in rows])
            self.resid_t2.setData(x, [r["resid_t2"] * 10 for r in rows])
        else:
            self.resid_t1.setData([], [])
            self.resid_t2.setData([], [])

    def _remove_selected(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if rows and self.current_run:
            name = self.table.item(rows[0], 0).text()
            self.removeAnchorRequested.emit(name)
