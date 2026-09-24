"""🏷️ Internal-standards panel.

A library of internal standards (the 7 shipped defaults, plus whatever a
user has added), a button to label the row selected in that table on
whichever run is open in the 2D view, a button to search every other
loaded run for the internal standards labelled so far, and the resulting
confidence report. Labelling itself is a single click driven by the main
window (like 🔍's "add a blob"); this panel starts and reports on it.
"""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

TIER_COLOURS = {
    "high": QtGui.QColor(60, 160, 90),
    "medium": QtGui.QColor(210, 160, 40),
    "low": QtGui.QColor(200, 90, 60),
    "not found": QtGui.QColor(140, 140, 140),
}


def _offset(value: float, expected: float, fmt: str) -> str:
    d = value - expected
    return format(d, fmt) if d == d else ""


class IstdPanel(QtWidgets.QWidget):
    labelRequested = QtCore.Signal(str)          # IS name, on whichever run is open
    findAllRequested = QtCore.Signal()
    removeReferenceRequested = QtCore.Signal(str)
    addCustomRequested = QtCore.Signal(str, str, str)   # name, cas, product_code
    removeFromLibraryRequested = QtCore.Signal(str)
    exportRequested = QtCore.Signal()
    rowActivated = QtCore.Signal(str, str)        # run, IS name -- re-label that one run by hand

    def __init__(self, parent=None):
        super().__init__(parent)
        self.istd = None
        self.current_run: str | None = None

        self.library_table = QtWidgets.QTableWidget(0, 3)
        self.library_table.setHorizontalHeaderLabels(["Name", "CAS", "Product code"])
        self.library_table.verticalHeader().setVisible(False)
        self.library_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.library_table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.library_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.library_table.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.library_table.itemSelectionChanged.connect(self._selection_changed)

        self.name_field = QtWidgets.QLineEdit()
        self.name_field.setPlaceholderText("Name")
        self.cas_field = QtWidgets.QLineEdit()
        self.cas_field.setPlaceholderText("CAS")
        self.code_field = QtWidgets.QLineEdit()
        self.code_field.setPlaceholderText("Product code (optional)")
        self.add_btn = QtWidgets.QPushButton("➕ Add to library")
        self.add_btn.clicked.connect(self._add_custom)
        self.remove_library_btn = QtWidgets.QPushButton("Remove from library")
        self.remove_library_btn.clicked.connect(self._remove_from_library)

        self.label_btn = QtWidgets.QPushButton("🏷️ Label peak on this run")
        self.label_btn.setToolTip("Click a peak for the selected internal standard on the run open in the 2D view")
        self.label_btn.clicked.connect(self._label_clicked)
        self.find_btn = QtWidgets.QPushButton("🔎 Find on all other chromatograms")
        self.find_btn.clicked.connect(self.findAllRequested)
        self.export_btn = QtWidgets.QPushButton("Export report (CSV)…")
        self.export_btn.clicked.connect(self.exportRequested)

        self.summary = QtWidgets.QLabel("Label an internal standard, then find it on every other run.")
        self.summary.setWordWrap(True)

        self.report_table = QtWidgets.QTableWidget(0, 11)
        self.report_table.setHorizontalHeaderLabels(
            ["Run", "IS", "t1 (min)", "t2 (s)", "Δt1 (s)", "Δt2 (s)", "off drift t1 (s)", "off drift t2 (s)",
             "spectral", "combined", "tier"])
        self.report_table.setToolTip(
            "Δt: offset from the reference peak. Off drift: offset from where the other internal "
            "standards on the same run predict this one should be; this is what the position score uses.")
        self.report_table.verticalHeader().setVisible(False)
        self.report_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.report_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.report_table.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.report_table.cellDoubleClicked.connect(self._row_activated)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(QtWidgets.QLabel("Internal standard library"))
        lay.addWidget(self.library_table, 1)
        form = QtWidgets.QHBoxLayout()
        form.addWidget(self.name_field)
        form.addWidget(self.cas_field)
        form.addWidget(self.code_field)
        lay.addLayout(form)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.add_btn)
        row.addWidget(self.remove_library_btn)
        lay.addLayout(row)
        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(self.label_btn)
        row2.addWidget(self.find_btn)
        lay.addLayout(row2)
        lay.addWidget(self.summary)
        lay.addWidget(QtWidgets.QLabel("Confidence report"))
        lay.addWidget(self.report_table, 2)
        lay.addWidget(self.export_btn)

    # ------------------------------------------------------------- wiring
    def set_istd(self, istd) -> None:
        self.istd = istd

    def set_current_run(self, name: str | None) -> None:
        self.current_run = name
        self.refresh()

    def _selection_changed(self) -> None:
        self.label_btn.setEnabled(bool(self.library_table.selectedItems()) and self.current_run is not None)

    def selected_name(self) -> str | None:
        rows = sorted({i.row() for i in self.library_table.selectedIndexes()})
        if not rows:
            return None
        item = self.library_table.item(rows[0], 0)
        return item.text() if item else None

    def _label_clicked(self) -> None:
        name = self.selected_name()
        if name:
            self.labelRequested.emit(name)

    def _add_custom(self) -> None:
        name = self.name_field.text().strip()
        if not name:
            return
        self.addCustomRequested.emit(name, self.cas_field.text().strip(), self.code_field.text().strip())
        self.name_field.clear()
        self.cas_field.clear()
        self.code_field.clear()

    def _remove_from_library(self) -> None:
        name = self.selected_name()
        if name:
            self.removeFromLibraryRequested.emit(name)

    def _row_activated(self, row: int, _col: int) -> None:
        run = self.report_table.item(row, 0).text()
        name = self.report_table.item(row, 1).text()
        self.rowActivated.emit(run, name)

    # ------------------------------------------------------------- filling
    def refresh(self) -> None:
        if self.istd is None:
            return
        selected = self.selected_name()
        self.library_table.setRowCount(len(self.istd.library))
        for i, entry in enumerate(self.istd.library):
            for c, v in enumerate((entry.name, entry.cas, entry.product_code)):
                item = self.library_table.item(i, c)
                if item is None:
                    item = QtWidgets.QTableWidgetItem()
                    self.library_table.setItem(i, c, item)
                item.setText(v)
            referenced = entry.name in self.istd.references
            for c in range(3):
                self.library_table.item(i, c).setForeground(
                    QtGui.QBrush(TIER_COLOURS["high"]) if referenced else QtGui.QBrush())
        if selected:
            hit = self.library_table.findItems(selected, QtCore.Qt.MatchFlag.MatchExactly)
            if hit:
                self.library_table.selectRow(hit[0].row())
        self.label_btn.setEnabled(bool(self.selected_name()) and self.current_run is not None)
        self.find_btn.setEnabled(bool(self.istd.references))

        matches = self.istd.all_matches()
        self.report_table.setRowCount(len(matches))
        for i, m in enumerate(sorted(matches, key=lambda m: (m.name, m.run))):
            values = (m.run, m.name, f"{m.t1 / 60:.3f}" if m.t1 == m.t1 else "",
                     f"{m.t2:.3f}" if m.t2 == m.t2 else "",
                     f"{m.dt1:+.2f}" if m.dt1 == m.dt1 else "", f"{m.dt2:+.3f}" if m.dt2 == m.dt2 else "",
                     _offset(m.t1, m.expected_t1, "+.2f"), _offset(m.t2, m.expected_t2, "+.3f"),
                     f"{m.spectral_score:.2f}", f"{m.combined_score:.2f}", m.tier)
            for c, v in enumerate(values):
                item = self.report_table.item(i, c)
                if item is None:
                    item = QtWidgets.QTableWidgetItem()
                    self.report_table.setItem(i, c, item)
                item.setText(v)
                item.setForeground(QtGui.QBrush(TIER_COLOURS.get(m.tier, QtGui.QColor("black"))))

        if not self.istd.references:
            self.summary.setText("Select an internal standard above, then 🏷️ label its peak on this run.")
        elif not matches:
            self.summary.setText(
                f"{len(self.istd.references)} internal standard(s) labelled on {self.istd.reference_run}. "
                "🔎 Find on all other chromatograms to match them elsewhere.")
        else:
            not_found = sum(1 for m in matches if m.tier == "not found")
            low = sum(1 for m in matches if m.tier == "low")
            self.summary.setText(
                f"{len(matches)} match(es) across {len(self.istd.matches)} run(s): "
                f"{not_found} not found, {low} low confidence. Double-click a row to re-label it by hand.")
