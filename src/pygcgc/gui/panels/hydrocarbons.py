"""⛓️ Hydrocarbons panel.

Two independent workflows share this panel: finding the n-alkane ladder
(a two-click flow, like 📍 add anchor), and building and lifting group-
type stencils, drawn once as an editable polygon on the 2D view and
carried onto any other run through its own alkane series.
"""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from ...core.hydrocarbons import HydrocarbonSettings, ratios_match

SERIES_COLUMNS = ["Cn", "t1 (min)", "t2 (s)"]
STENCIL_COLUMNS = ["Stencil", "response factor", "volume", "n blobs", "amount"]
DEFAULT_QUALIFIERS = HydrocarbonSettings()
PAIR_PRESETS = ["57/71 n-alkane check", "91/105 alkylbenzenes", "82/83 alkylcyclohexanes",
                "142/156 alkylnaphthalenes", "58/59"]


class HydrocarbonsPanel(QtWidgets.QWidget):
    settingsRequested = QtCore.Signal()
    editStencilRequested = QtCore.Signal(str)
    findHomologRequested = QtCore.Signal(int, int, str)     # major, minor, name
    homologActivated = QtCore.Signal(int)                    # blob index
    qualifiersChanged = QtCore.Signal()
    pickQualifiersToggled = QtCore.Signal(bool)
    findRequested = QtCore.Signal()
    propagateRequested = QtCore.Signal()
    drawRequested = QtCore.Signal()
    saveStencilRequested = QtCore.Signal(str, float)          # name, response factor
    cancelDrawRequested = QtCore.Signal()
    removeStencilsRequested = QtCore.Signal(list)              # names
    responseFactorEdited = QtCore.Signal(str, float)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_run: str | None = None
        self._updating = False

        self.find_btn = QtWidgets.QPushButton("➕ Find alkane series")
        self.find_btn.setToolTip("Click two n-alkanes in the 2D view and type their carbon numbers")
        self.find_btn.clicked.connect(self.findRequested)
        self.propagate_btn = QtWidgets.QPushButton("⛓️ Propagate to other runs")
        self.propagate_btn.setToolTip("Search every other loaded run for this run's confirmed alkanes")
        self.propagate_btn.clicked.connect(self.propagateRequested)
        self._qual_ions: tuple[int, ...] = DEFAULT_QUALIFIERS.qual_ions
        self._quant_ion = DEFAULT_QUALIFIERS.quant_ion
        self._settings = DEFAULT_QUALIFIERS
        self._qual_spins: list[QtWidgets.QDoubleSpinBox] = []
        self.qual_box = QtWidgets.QGroupBox("Alkane qualifier ratios")
        self.qual_form = QtWidgets.QFormLayout()
        self.qual_tol = QtWidgets.QDoubleSpinBox()
        self.qual_tol.setRange(0.05, 1.0)
        self.qual_tol.setSingleStep(0.05)
        self.qual_tol.setDecimals(2)
        self.qual_tol.setToolTip("Allowed relative deviation of each ratio from its expected value")
        self.qual_tol.valueChanged.connect(self._qual_edited)
        self.qual_enforce = QtWidgets.QCheckBox("Reject candidates outside the tolerance")
        self.qual_enforce.setToolTip("Off: the search still prefers a matching pulse but falls back to the "
                                     "tallest one, flagged as a mismatch")
        self.qual_enforce.toggled.connect(self._qual_edited)
        self.qual_pick = QtWidgets.QPushButton("🖱 Pick from 2D view")
        self.qual_pick.setCheckable(True)
        self.qual_pick.setToolTip("Click an n-alkane in the 2D view; its measured ratios become the expected ones")
        self.qual_pick.toggled.connect(self.pickQualifiersToggled)
        self.qual_reset = QtWidgets.QPushButton("Reset to measured default")
        self.qual_reset.clicked.connect(self._qual_reset)
        qlay = QtWidgets.QVBoxLayout(self.qual_box)
        qlay.addLayout(self.qual_form)
        qrow = QtWidgets.QFormLayout()
        qrow.addRow("Tolerance (relative)", self.qual_tol)
        qlay.addLayout(qrow)
        qlay.addWidget(self.qual_enforce)
        qbtn = QtWidgets.QHBoxLayout()
        qbtn.addWidget(self.qual_pick)
        qbtn.addWidget(self.qual_reset)
        qlay.addLayout(qbtn)

        self.series_table = QtWidgets.QTableWidget(0, len(SERIES_COLUMNS))
        self.series_table.setHorizontalHeaderLabels(SERIES_COLUMNS)
        self.series_table.verticalHeader().setVisible(False)
        self.series_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.series_table.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.series_table.setMaximumHeight(140)
        self.series_status = QtWidgets.QLabel("No alkane series on this run yet.")
        self.series_status.setWordWrap(True)

        self.draw_btn = QtWidgets.QPushButton("✏️ Draw stencil")
        self.draw_btn.setCheckable(True)
        self.draw_btn.setToolTip("Drag a polygon onto the 2D view, over the group you want to lump")
        self.draw_btn.toggled.connect(self._draw_toggled)
        self.name_edit = QtWidgets.QLineEdit()
        self.name_edit.setPlaceholderText("Stencil name, e.g. C15-C25 branched")
        self.rf_spin = QtWidgets.QDoubleSpinBox()
        self.rf_spin.setRange(0.001, 1000.0)
        self.rf_spin.setValue(1.0)
        self.rf_spin.setDecimals(3)
        self.rf_spin.setToolTip("Surrogate response factor: amount = volume × this")
        self.save_btn = QtWidgets.QPushButton("💾 Save")
        self.save_btn.setEnabled(False)
        self.save_btn.clicked.connect(self._save_clicked)

        self.stencil_table = QtWidgets.QTableWidget(0, len(STENCIL_COLUMNS))
        self.stencil_table.setHorizontalHeaderLabels(STENCIL_COLUMNS)
        self.stencil_table.verticalHeader().setVisible(False)
        self.stencil_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.stencil_table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.stencil_table.setEditTriggers(
            QtWidgets.QAbstractItemView.EditTrigger.DoubleClicked
            | QtWidgets.QAbstractItemView.EditTrigger.SelectedClicked)
        self.stencil_table.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.stencil_table.itemChanged.connect(self._item_changed)
        self.delete_btn = QtWidgets.QPushButton("Delete selected")
        self.delete_btn.clicked.connect(self._delete_selected)
        self.edit_btn = QtWidgets.QPushButton("✎ Edit shape")
        self.edit_btn.setToolTip("Load the selected stencil into the editable polygon, drag its corners, "
                                 "then 💾 Save (same name) to replace it")
        self.edit_btn.clicked.connect(self._edit_selected)

        self.markers_label = QtWidgets.QLabel("")
        self.markers_label.setWordWrap(True)

        self.settings_btn = QtWidgets.QPushButton("⚙️ All settings…")
        self.settings_btn.setToolTip("Every hydrocarbon setting: search windows, thresholds, series tolerances")
        self.settings_btn.clicked.connect(self.settingsRequested)

        self.pair_combo = QtWidgets.QComboBox()
        self.pair_combo.setEditable(True)
        for text in PAIR_PRESETS:
            self.pair_combo.addItem(text)
        self.pair_combo.setToolTip("Major/minor diagnostic ions, e.g. 91/105. Needs 🔍 blobs and an alkane series "
                                   "on this run. Unvalidated on real data except 57/71.")
        self.find_homolog_btn = QtWidgets.QPushButton("🔗 Find series")
        self.find_homolog_btn.clicked.connect(self._find_homolog_clicked)
        self.homolog_table = QtWidgets.QTableWidget(0, 4)
        self.homolog_table.setHorizontalHeaderLabels(["Carbon eq.", "t1 (min)", "t2 (s)", "Blob"])
        self.homolog_table.verticalHeader().setVisible(False)
        self.homolog_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.homolog_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.homolog_table.setMaximumHeight(150)
        self.homolog_table.cellDoubleClicked.connect(self._homolog_activated)
        self.homolog_status = QtWidgets.QLabel("")
        self.homolog_status.setWordWrap(True)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.find_btn)
        row.addWidget(self.propagate_btn)
        lay.addLayout(row)
        lay.addWidget(self.qual_box)
        lay.addWidget(self.series_table)
        lay.addWidget(self.series_status)
        lay.addWidget(_hline())
        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(self.draw_btn)
        row2.addWidget(QtWidgets.QLabel("Name"))
        row2.addWidget(self.name_edit, 1)
        row2.addWidget(QtWidgets.QLabel("RF"))
        row2.addWidget(self.rf_spin)
        row2.addWidget(self.save_btn)
        lay.addLayout(row2)
        lay.addWidget(self.stencil_table, 1)
        brow = QtWidgets.QHBoxLayout()
        brow.addWidget(self.edit_btn)
        brow.addWidget(self.delete_btn)
        lay.addLayout(brow)
        lay.addWidget(self.markers_label)
        lay.addWidget(_hline())
        hrow = QtWidgets.QHBoxLayout()
        hrow.addWidget(QtWidgets.QLabel("Homologous series, ions"))
        hrow.addWidget(self.pair_combo, 1)
        hrow.addWidget(self.find_homolog_btn)
        lay.addLayout(hrow)
        lay.addWidget(self.homolog_table)
        lay.addWidget(self.homolog_status)
        lay.addWidget(self.settings_btn)
        self._refresh_buttons()

    # ------------------------------------------------------------- filling
    # ------------------------------------------------------ homologous series
    def _find_homolog_clicked(self) -> None:
        text = self.pair_combo.currentText().split()[0]
        try:
            major, minor = (int(x) for x in text.split("/"))
        except ValueError:
            self.homolog_status.setText("Type the ions as major/minor, e.g. 91/105.")
            return
        self.findHomologRequested.emit(major, minor, self.pair_combo.currentText().strip())

    def set_homologs(self, series, message: str = "") -> None:
        """Show a `HomologSeries` (or clear the table for `None`)."""
        members = [] if series is None else series.members
        self.homolog_table.setRowCount(len(members))
        for r, m in enumerate(members):
            for c, v in enumerate((f"{m['carbon']:.2f}", f"{m['t1'] / 60:.3f}", f"{m['t2']:.2f}", str(m["blob"]))):
                self.homolog_table.setItem(r, c, QtWidgets.QTableWidgetItem(v))
        self.homolog_status.setText(message)

    def _homolog_activated(self, row: int, _col: int) -> None:
        item = self.homolog_table.item(row, 3)
        if item is not None:
            self.homologActivated.emit(int(item.text()))

    # ---------------------------------------------------------- qualifiers
    def set_qualifiers(self, s: HydrocarbonSettings) -> None:
        """Show `s`'s expected ratios, tolerance and gate without emitting a change."""
        self._settings = s
        self._updating = True
        if tuple(s.qual_ions) != self._qual_ions or len(self._qual_spins) != len(s.qual_ions):
            while self.qual_form.rowCount():
                self.qual_form.removeRow(0)
            self._qual_spins = []
            for ion in s.qual_ions:
                sp = QtWidgets.QDoubleSpinBox()
                sp.setRange(0.0, 10.0)
                sp.setSingleStep(0.01)
                sp.setDecimals(3)
                sp.valueChanged.connect(self._qual_edited)
                self._qual_spins.append(sp)
                self.qual_form.addRow(f"m/z {ion} / {s.quant_ion}", sp)
        self._qual_ions, self._quant_ion = tuple(s.qual_ions), s.quant_ion
        for sp, v in zip(self._qual_spins, s.qual_ratios or (0.0,) * len(self._qual_spins)):
            sp.setValue(float(v))
        self.qual_tol.setValue(s.qual_tol)
        self.qual_enforce.setChecked(s.enforce_qualifiers)
        self._updating = False

    def qualifiers(self) -> tuple[tuple[float, ...], float, bool]:
        return (tuple(sp.value() for sp in self._qual_spins), self.qual_tol.value(),
                self.qual_enforce.isChecked())

    def set_picking_qualifiers(self, on: bool) -> None:
        self.qual_pick.blockSignals(True)
        self.qual_pick.setChecked(on)
        self.qual_pick.blockSignals(False)

    def _qual_edited(self, *_a) -> None:
        if not self._updating:
            self.qualifiersChanged.emit()

    def _qual_reset(self) -> None:
        d = DEFAULT_QUALIFIERS
        self.set_qualifiers(HydrocarbonSettings(qual_ratios=d.qual_ratios, qual_tol=d.qual_tol,
                                                enforce_qualifiers=d.enforce_qualifiers))
        self.qualifiersChanged.emit()

    # ------------------------------------------------------------- filling
    def set_current_run(self, name: str | None, series, running: bool) -> None:
        self.current_run = name
        self._updating = True
        n = 0 if series is None else len(series)
        headers = SERIES_COLUMNS + [f"{i}/{self._quant_ion}" for i in self._qual_ions]
        self.series_table.setColumnCount(len(headers))
        self.series_table.setHorizontalHeaderLabels(headers)
        self.series_table.setRowCount(n)
        bad = QtGui.QBrush(QtGui.QColor(226, 160, 60, 90))
        if series is not None:
            for i, c in enumerate(sorted(series.hits)):
                h = series.hits[c]
                values = [f"C{c}", f"{h.t1 / 60:.3f}", f"{h.t2:.3f}"] + [f"{r:.2f}" for r in h.qual_ratios]
                # A seed is trusted even if its ratios disagree; other hits are flagged when they do.
                off = not ratios_match(h.qual_ratios, self._settings)
                for col, v in enumerate(values):
                    item = self.series_table.item(i, col)
                    if item is None:
                        item = QtWidgets.QTableWidgetItem()
                        self.series_table.setItem(i, col, item)
                    item.setText(v)
                    item.setBackground(bad if (off and col >= len(SERIES_COLUMNS)) else QtGui.QBrush())
                    if col >= len(SERIES_COLUMNS) and off:
                        item.setToolTip("Outside the expected qualifier ratios")
        self._updating = False
        if running:
            self.series_status.setText(f"⛓️ Working on {name}…" if name else "⛓️ Working…")
        elif series is None:
            self.series_status.setText("No alkane series on this run yet." if name else "")
        else:
            self.series_status.setText(f"{len(series)} confirmed, C{min(series.hits)} to C{max(series.hits)}.")
        self.propagate_btn.setEnabled(series is not None and len(series) > 0)
        self._refresh_buttons()

    def set_stencils(self, stencils, totals: dict, name_field_enabled: bool) -> None:
        self._updating = True
        self.stencil_table.setRowCount(len(stencils))
        for i, s in enumerate(stencils):
            t = totals.get(s.name, {})
            values = (s.name, f"{s.response_factor:g}",
                     f"{t.get('volume', 0):,.0f}" if t else "", str(t.get("n_blobs", "")) if t else "",
                     f"{t.get('amount', 0):,.3g}" if t else "")
            for c, v in enumerate(values):
                item = self.stencil_table.item(i, c)
                if item is None:
                    item = QtWidgets.QTableWidgetItem()
                    self.stencil_table.setItem(i, c, item)
                item.setText(v)
                editable = c == 1
                flags = item.flags()
                item.setFlags((flags | QtCore.Qt.ItemFlag.ItemIsEditable) if editable
                              else (flags & ~QtCore.Qt.ItemFlag.ItemIsEditable))
        self._updating = False

    def set_markers(self, markers: dict) -> None:
        if not markers:
            self.markers_label.setText("")
            return
        parts = [f"{name}: blob #{i}" if i is not None else f"{name}: not found"
                for name, i in markers.items()]
        self.markers_label.setText("Markers -- " + ", ".join(parts))

    def set_drawing(self, on: bool) -> None:
        self.draw_btn.blockSignals(True)
        self.draw_btn.setChecked(on)
        self.draw_btn.blockSignals(False)
        self.save_btn.setEnabled(on)

    def selected_stencil_names(self) -> list[str]:
        rows = sorted({i.row() for i in self.stencil_table.selectedIndexes()})
        return [self.stencil_table.item(r, 0).text() for r in rows]

    def _refresh_buttons(self) -> None:
        self.find_btn.setEnabled(self.current_run is not None)

    def _draw_toggled(self, on: bool) -> None:
        self.save_btn.setEnabled(on)
        (self.drawRequested if on else self.cancelDrawRequested).emit()

    def _save_clicked(self) -> None:
        name = self.name_edit.text().strip()
        if not name:
            QtWidgets.QMessageBox.warning(self, "Name needed", "Give the stencil a name first.")
            return
        self.saveStencilRequested.emit(name, self.rf_spin.value())
        self.name_edit.clear()

    def _edit_selected(self) -> None:
        names = self.selected_stencil_names()
        if len(names) != 1:
            return
        self.name_edit.setText(names[0])
        self.editStencilRequested.emit(names[0])

    def _delete_selected(self) -> None:
        names = self.selected_stencil_names()
        if names:
            self.removeStencilsRequested.emit(names)

    def _item_changed(self, item: QtWidgets.QTableWidgetItem) -> None:
        if self._updating or item.column() != 1:
            return
        name = self.stencil_table.item(item.row(), 0).text()
        try:
            value = float(item.text())
        except ValueError:
            return
        self.responseFactorEdited.emit(name, value)


def _hline() -> QtWidgets.QFrame:
    line = QtWidgets.QFrame()
    line.setFrameShape(QtWidgets.QFrame.Shape.HLine)
    line.setFrameShadow(QtWidgets.QFrame.Shadow.Sunken)
    return line
