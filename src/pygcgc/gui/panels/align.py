"""📍 Align panel.

The reference run, a shared table of clicks, guided picking across runs,
bounded local-maximum refinement and transform diagnostics.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

KIND_LABELS = [("Auto", "auto"), ("Shift", "shift"), ("Affine", "affine"), ("Smooth", "smooth")]


class AlignPanel(QtWidgets.QWidget):
    referenceChanged = QtCore.Signal(str)              # run name, or "" for none
    addAnchorRequested = QtCore.Signal()
    removeAnchorRequested = QtCore.Signal(str, str)     # anchor name, run
    skipAnchorRequested = QtCore.Signal()
    refineRequested = QtCore.Signal()
    clearRunRequested = QtCore.Signal()
    kindChanged = QtCore.Signal(str)                    # "auto" / "shift" / "affine" / "smooth"
    propagateRequested = QtCore.Signal()
    showAlignedToggled = QtCore.Signal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.align = None
        self.current_run: str | None = None
        self._updating = False
        self._picking = False

        self.reference = QtWidgets.QComboBox()
        self.reference.setToolTip("Every other run is aligned to this one")
        self.reference.currentIndexChanged.connect(self._reference_changed)

        self.add_btn = QtWidgets.QPushButton("Pick anchors across runs (A)")
        self.add_btn.setToolTip("Click on the reference, then every loaded run; repeat for the next anchor")
        self.add_btn.clicked.connect(self.addAnchorRequested)
        self.skip_btn = QtWidgets.QPushButton("Skip peak in this run")
        self.skip_btn.clicked.connect(self.skipAnchorRequested)
        self.skip_btn.hide()
        self.pick_status = QtWidgets.QLabel("Start picking to visit each chromatogram in turn. Esc finishes.")
        self.pick_status.setWordWrap(True)
        self.refine_btn = QtWidgets.QPushButton("Refine anchors to local maxima")
        self.refine_btn.setToolTip("Adjust every saved click, including the reference, using its saved ion or TIC. Undo restores all positions.")
        self.refine_btn.clicked.connect(self.refineRequested)
        self.refine_columns = QtWidgets.QSpinBox()
        self.refine_columns.setRange(0, 3)
        self.refine_columns.setValue(1)
        self.refine_columns.setSuffix(" modulation(s)")
        self.refine_t2 = QtWidgets.QDoubleSpinBox()
        self.refine_t2.setRange(.01, 1.)
        self.refine_t2.setSingleStep(.025)
        self.refine_t2.setDecimals(3)
        self.refine_t2.setValue(.15)
        self.refine_t2.setSuffix(" s")
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

        self.table = QtWidgets.QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(["Anchor", "Run", "ref t1 (min)", "ref t2 (s)",
                                             "t1 (min)", "t2 (s)", "residual"])
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setMinimumHeight(160)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.Interactive)
        self.table.setColumnWidth(1, 100)
        for column in (2, 3, 6):
            self.table.setColumnHidden(column, True)
        remove_row = QtWidgets.QPushButton("Remove selected anchor")
        remove_row.setToolTip("Remove the selected table row. A reference row removes its matching points on every run.")
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
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addLayout(form)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.add_btn)
        lay.addLayout(row)
        lay.addWidget(self.pick_status)
        lay.addWidget(self.skip_btn)
        lay.addWidget(self.table, 1)
        lay.addWidget(self.summary)
        lay.addWidget(self.refine_btn)
        lay.addWidget(self.show_aligned)
        self.more_btn = QtWidgets.QToolButton()
        self.more_btn.setText("More alignment options")
        self.more_btn.setCheckable(True)
        self.more_btn.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.more_btn.setArrowType(QtCore.Qt.ArrowType.RightArrow)
        self.details = QtWidgets.QWidget()
        details = QtWidgets.QVBoxLayout(self.details)
        details.setContentsMargins(0, 0, 0, 0)
        limits = QtWidgets.QFormLayout()
        limits.addRow("Search t1 ±", self.refine_columns)
        limits.addRow("Search t2 ±", self.refine_t2)
        limits.addRow("Transform", self.kind)
        details.addLayout(limits)
        details.addWidget(self.propagate_btn)
        details.addWidget(remove_row)
        details.addWidget(self.clear_btn)
        details.addWidget(self.resid_plot)
        self.details.hide()
        self.more_btn.toggled.connect(self._show_details)
        lay.addWidget(self.more_btn)
        lay.addWidget(self.details)

    def _show_details(self, visible: bool) -> None:
        self.details.setVisible(visible)
        self.more_btn.setArrowType(QtCore.Qt.ArrowType.DownArrow if visible else QtCore.Qt.ArrowType.RightArrow)
        for column in (2, 3, 6):
            self.table.setColumnHidden(column, not visible)

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
        desired = (self.align.reference or "") if self.align is not None else cur
        i = self.reference.findData(desired)
        self.reference.setCurrentIndex(max(i, 0))
        self._updating = False

    def set_current_run(self, name: str | None) -> None:
        self.current_run = name
        self.refresh()

    def set_picking(self, message: str = "", can_skip: bool = False) -> None:
        self._picking = bool(message)
        self.pick_status.setText(message or "Start picking to visit each chromatogram in turn. Esc finishes.")
        self.add_btn.setText("Finish picking (Esc)" if message else "Pick anchors across runs (A)")
        self.skip_btn.setVisible(bool(message) and can_skip)
        self.refresh()

    def refresh(self) -> None:
        a = self.align
        if a is None:
            return
        is_ref = self.current_run is not None and self.current_run == a.reference
        self.add_btn.setEnabled(self.current_run is not None and not a.refining)
        self.reference.setEnabled(not self._picking and not a.refining)
        self.refine_btn.setEnabled(bool(a.reference and any(a.points.values())) and not a.refining and not self._picking)
        self.refine_btn.setText("Refining anchors..." if a.refining else "Refine anchors to local maxima")
        self.propagate_btn.setEnabled(bool(a.reference) and not self._picking and not a.refining)
        self.clear_btn.setEnabled(bool(self.current_run and a.points.get(self.current_run)) and not is_ref)
        self._updating = True
        self.kind.setCurrentIndex(max(self.kind.findData(a.kind.get(self.current_run, "auto")), 0))
        self._updating = False
        self.show_aligned.setEnabled(bool(a.reference) and not self._picking)

        pts = a.points.get(self.current_run, []) if self.current_run else []
        tf = a.transforms.get(self.current_run) if self.current_run else None
        all_points = [(run, p) for run, points in a.points.items() for p in points]
        self.table.setRowCount(len(all_points))
        for i, (run, p) in enumerate(all_points):
            transform = a.transforms.get(run)
            residuals = {r["name"]: r for r in transform.residuals()} if transform else {}
            r = residuals.get(p.name)
            resid_text = f"{r['resid_t1']:+.2f}s / {r['resid_t2'] * 1e3:+.0f}ms" if r else ""
            values = (p.name, run, f"{p.ref_t1 / 60:.3f}", f"{p.ref_t2:.3f}",
                     f"{p.run_t1 / 60:.3f}", f"{p.run_t2:.3f}", resid_text)
            for c, v in enumerate(values):
                item = self.table.item(i, c)
                if item is None:
                    item = QtWidgets.QTableWidgetItem()
                    self.table.setItem(i, c, item)
                item.setText(v)
                item.setToolTip(f"{run}: {p.name}\nReference: {p.ref_t1/60:.5f} min, {p.ref_t2:.5f} s\n"
                                f"This run: {p.run_t1/60:.5f} min, {p.run_t2:.5f} s")
                item.setData(QtCore.Qt.ItemDataRole.UserRole, run)
                font = item.font()
                font.setBold(run == self.current_run)
                item.setFont(font)

        if is_ref:
            self.summary.setText(f"{self.current_run} is the reference: shown as itself, never aligned.")
        elif not a.reference:
            self.summary.setText("Start picking to use the first loaded run as the reference, or choose one above.")
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
            run = self.table.item(rows[0], 0).data(QtCore.Qt.ItemDataRole.UserRole)
            self.removeAnchorRequested.emit(name, run)
