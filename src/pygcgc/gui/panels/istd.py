"""🏷️ Internal-standards panel and results table.

`IstdPanel` (the workflow step on the left): the internal standards, with
where each was labelled and on how many runs it was found; the three
thresholds a peak must pass (first-dimension and second-dimension distance
from the expected position, and spectral match); and the buttons to label a
standard by clicking its peak and to search every chromatogram. Go runs the
search.

`IstdResultsTable` (a dock along the bottom): one row per standard, one
column per chromatogram. Each cell says whether the standard passed all
three thresholds there (✓ or ✗, and its colour) and shows the quantity
chosen above the table; its tooltip gives every number and, for a miss,
which threshold failed. Clicking a cell opens that run on the peak;
double-clicking re-labels the standard there by hand.
"""

from __future__ import annotations

from dataclasses import replace

from PySide6 import QtCore, QtGui, QtWidgets

from ...core.istd import FOUND, NO_PEAKS, IstdMatchSettings

STATUS_COLOURS = {
    FOUND: QtGui.QColor(46, 160, 90),
    "ambiguous": QtGui.QColor(215, 160, 30),     # found, but more than one peak passed
    "not found": QtGui.QColor(215, 85, 60),
    NO_PEAKS: QtGui.QColor(140, 140, 140),
}


def _fmt(value: float, fmt: str) -> str:
    return format(value, fmt) if value == value else ""


# What a results cell can show: label -> text of a match.
SHOW = {
    "Spectral match": lambda m: f"{m.spectral_score:.2f}" if m.t1 == m.t1 else "",
    "t1 (min)": lambda m: _fmt(m.t1 / 60.0, ".2f"),
    "t2 (s)": lambda m: _fmt(m.t2, ".2f"),
    "Δt1 from expected (s)": lambda m: _fmt(m.off_t1, "+.1f"),
    "Δt2 from expected (s)": lambda m: _fmt(m.off_t2, "+.2f"),
    "Δt1 from reference (s)": lambda m: _fmt(m.dt1, "+.1f"),
    "Δt2 from reference (s)": lambda m: _fmt(m.dt2, "+.2f"),
}


def _colour_key(m) -> str:
    if m.found:
        return "ambiguous" if m.n_passing > 1 else FOUND
    return m.status


def _tooltip(m, settings: IstdMatchSettings) -> str:
    head = {FOUND: "Found", NO_PEAKS: "No peaks"}.get(m.status, "Not found")
    lines = [f"{m.name} on {m.run}: {head}"]
    if m.t1 == m.t1:
        lines += [f"t1 {m.t1 / 60:.3f} min, t2 {m.t2:.3f} s",
                  f"Δ from expected: t1 {m.off_t1:+.1f} s (limit ±{settings.max_dt1:g}), "
                  f"t2 {m.off_t2:+.2f} s (limit ±{settings.max_dt2:g})",
                  f"Δ from reference: t1 {m.dt1:+.1f} s, t2 {m.dt2:+.2f} s",
                  f"Spectral match {m.spectral_score:.2f} (minimum {settings.min_spectral:.2f}; "
                  f"forward {m.forward:.2f}, reverse {m.reverse:.2f})"]
    if m.expected_t1 == m.expected_t1:
        lines.append(f"Expected at t1 {m.expected_t1 / 60:.3f} min, t2 {m.expected_t2:.3f} s"
                     + (f", from the drift of {m.n_drift_anchors} other standard(s)" if m.n_drift_anchors
                        else ", the reference position"))
    if m.reason:
        lines.append(m.reason)
    return "\n".join(lines)


class IstdPanel(QtWidgets.QWidget):
    labelRequested = QtCore.Signal(str)          # IS name, on whichever run is open
    findAllRequested = QtCore.Signal()
    removeReferenceRequested = QtCore.Signal(str)
    addCustomRequested = QtCore.Signal(str, str, str)   # name, cas, product_code
    removeFromLibraryRequested = QtCore.Signal(str)
    exportRequested = QtCore.Signal()
    targetsRequested = QtCore.Signal()            # write the labelled standards as target-table rows
    settingsChanged = QtCore.Signal(object)       # IstdMatchSettings

    def __init__(self, parent=None):
        super().__init__(parent)
        self.istd = None
        self.current_run: str | None = None
        self._updating = False

        self.library_table = QtWidgets.QTableWidget(0, 6)
        self.library_table.setHorizontalHeaderLabels(
            ["Name", "CAS", "Labelled on", "t1 (min)", "t2 (s)", "Found in"])
        self.library_table.verticalHeader().setVisible(False)
        self.library_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.library_table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.library_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.library_table.horizontalHeader().setSectionResizeMode(
            QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.library_table.setToolTip("Every internal standard in the library. Select one, then 🏷️ label its "
                                      "peak; its reference position and the runs it was found in appear here.")
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
        self.unlabel_btn = QtWidgets.QPushButton("Clear label")
        self.unlabel_btn.setToolTip("Forget the labelled peak of the selected standard, and its results")
        self.unlabel_btn.clicked.connect(self._unlabel)

        self.label_btn = QtWidgets.QPushButton("🏷️ Label peak on this run")
        self.label_btn.setToolTip("Click the peak of the selected internal standard on the run open in the 2D view. "
                                  "Its position and spectrum become the reference every run is searched against.")
        self.label_btn.clicked.connect(self._label_clicked)

        # The three thresholds, one per dimension.
        self.max_dt1 = QtWidgets.QDoubleSpinBox()
        self.max_dt1.setRange(0.5, 900.0)
        self.max_dt1.setSingleStep(1.0)
        self.max_dt1.setDecimals(1)
        self.max_dt1.setSuffix(" s")
        self.max_dt1.setToolTip("First dimension: how far (either side) a peak may sit from the expected t1")
        self.max_dt2 = QtWidgets.QDoubleSpinBox()
        self.max_dt2.setRange(0.01, 20.0)
        self.max_dt2.setSingleStep(0.05)
        self.max_dt2.setDecimals(2)
        self.max_dt2.setSuffix(" s")
        self.max_dt2.setToolTip("Second dimension: how far (either side) a peak may sit from the expected t2")
        self.min_spectral = QtWidgets.QDoubleSpinBox()
        self.min_spectral.setRange(0.0, 1.0)
        self.min_spectral.setSingleStep(0.05)
        self.min_spectral.setDecimals(2)
        self.min_spectral.setToolTip("Mass spectrum: the lowest match to the labelled peak's spectrum that counts "
                                     "(weighted dot product, 1 = identical)")
        self.follow_drift = QtWidgets.QCheckBox("Expect each standard where the others moved to")
        self.follow_drift.setToolTip("On each run, move the expected position by the median offset of the other "
                                     "standards found there with a confident spectrum, so a run-wide retention "
                                     "shift does not break the position thresholds. Off: the labelled position.")
        for w in (self.max_dt1, self.max_dt2, self.min_spectral):
            w.valueChanged.connect(self._settings_edited)
        self.follow_drift.toggled.connect(self._settings_edited)
        thresholds = QtWidgets.QGroupBox("Thresholds: a peak must pass all three")
        form = QtWidgets.QFormLayout(thresholds)
        form.addRow("Max Δt1", self.max_dt1)
        form.addRow("Max Δt2", self.max_dt2)
        form.addRow("Min spectral match", self.min_spectral)
        form.addRow(self.follow_drift)

        self.find_btn = QtWidgets.QPushButton("🔎 Find in every chromatogram")
        self.find_btn.clicked.connect(self.findAllRequested)
        self.export_btn = QtWidgets.QPushButton("Export results (CSV)…")
        self.export_btn.clicked.connect(self.exportRequested)
        self.targets_btn = QtWidgets.QPushButton("📋 Use as target table…")
        self.targets_btn.setToolTip("Write the labelled standards as internal-standard rows of a target table and "
                                    "load it, so they anchor the phase and normalise 🧪 quantification")
        self.targets_btn.clicked.connect(self.targetsRequested)

        self.summary = QtWidgets.QLabel("")
        self.summary.setWordWrap(True)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(self.library_table, 1)
        add_row = QtWidgets.QHBoxLayout()
        add_row.addWidget(self.name_field)
        add_row.addWidget(self.cas_field)
        add_row.addWidget(self.code_field)
        lay.addLayout(add_row)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.add_btn)
        row.addWidget(self.remove_library_btn)
        lay.addLayout(row)
        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(self.label_btn)
        row2.addWidget(self.unlabel_btn)
        lay.addLayout(row2)
        lay.addWidget(thresholds)
        lay.addWidget(self.find_btn)
        lay.addWidget(self.summary)
        row3 = QtWidgets.QHBoxLayout()
        row3.addWidget(self.export_btn)
        row3.addWidget(self.targets_btn)
        lay.addLayout(row3)

    # ------------------------------------------------------------- wiring
    def set_istd(self, istd) -> None:
        self.istd = istd
        self.set_settings(istd.settings)

    def set_settings(self, s: IstdMatchSettings) -> None:
        self._updating = True
        self.max_dt1.setValue(s.max_dt1)
        self.max_dt2.setValue(s.max_dt2)
        self.min_spectral.setValue(s.min_spectral)
        self.follow_drift.setChecked(s.follow_drift)
        self._updating = False

    def settings(self) -> IstdMatchSettings:
        base = self.istd.settings if self.istd is not None else IstdMatchSettings()
        return replace(base, max_dt1=float(self.max_dt1.value()), max_dt2=float(self.max_dt2.value()),
                       min_spectral=float(self.min_spectral.value()), follow_drift=self.follow_drift.isChecked())

    def _settings_edited(self, *_) -> None:
        if not self._updating:
            self.settingsChanged.emit(self.settings())

    def set_current_run(self, name: str | None) -> None:
        self.current_run = name
        self.refresh()

    def _selection_changed(self) -> None:
        name = self.selected_name()
        self.label_btn.setEnabled(bool(name) and self.current_run is not None)
        self.unlabel_btn.setEnabled(bool(name) and self.istd is not None and name in self.istd.references)

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

    def _unlabel(self) -> None:
        name = self.selected_name()
        if name:
            self.removeReferenceRequested.emit(name)

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

    # ------------------------------------------------------------- filling
    def refresh(self) -> None:
        if self.istd is None:
            return
        istd = self.istd
        selected = self.selected_name()
        n_runs = len(istd.matches)
        self.library_table.setRowCount(len(istd.library))
        for i, entry in enumerate(istd.library):
            ref = istd.references.get(entry.name)
            found = sum(1 for per_run in istd.matches.values()
                        if entry.name in per_run and per_run[entry.name].found)
            values = (entry.name, entry.cas, ref.run if ref else "",
                      f"{ref.t1 / 60:.2f}" if ref else "", f"{ref.t2:.2f}" if ref else "",
                      f"{found}/{n_runs}" if ref and n_runs else "")
            for c, v in enumerate(values):
                item = self.library_table.item(i, c)
                if item is None:
                    item = QtWidgets.QTableWidgetItem()
                    self.library_table.setItem(i, c, item)
                item.setText(v)
                item.setForeground(QtGui.QBrush(STATUS_COLOURS[FOUND]) if ref else QtGui.QBrush())
        if selected:
            hit = [r for r in range(self.library_table.rowCount())
                   if self.library_table.item(r, 0).text() == selected]
            if hit:
                self.library_table.selectRow(hit[0])
        self._selection_changed()
        self.find_btn.setEnabled(bool(istd.references))
        self.targets_btn.setEnabled(bool(istd.references))
        self.export_btn.setEnabled(bool(istd.matches))

        matches = istd.all_matches()
        if not istd.references:
            self.summary.setText("Select an internal standard, then 🏷️ label its peak on the run open in the "
                                 "2D view. Then set the thresholds and press Go to search every chromatogram.")
        elif not matches:
            self.summary.setText(f"{len(istd.references)} standard(s) labelled. Press Go to search every "
                                 "chromatogram with these thresholds.")
        else:
            found = sum(m.found for m in matches)
            self.summary.setText(f"{found} of {len(matches)} standard-run pairs pass all three thresholds, over "
                                 f"{n_runs} chromatogram(s). The full table is in the 🏷️ Internal standards dock.")


class IstdResultsTable(QtWidgets.QWidget):
    """Standards down, chromatograms across."""
    cellChosen = QtCore.Signal(str, str)          # run, standard: open the run on it
    relabelRequested = QtCore.Signal(str, str)    # run, standard: label it there by hand

    def __init__(self, parent=None):
        super().__init__(parent)
        self.istd = None
        self.show_box = QtWidgets.QComboBox()
        self.show_box.addItems(list(SHOW))
        self.show_box.currentIndexChanged.connect(lambda _i: self.refresh())
        self.legend = QtWidgets.QLabel("✓ passed all three thresholds (amber: more than one peak did)   "
                                       "✗ did not: hover for why   · no peaks detected.   "
                                       "Click: open that run on it.  Double-click: label it there by hand.")
        self.legend.setWordWrap(True)
        self.table = QtWidgets.QTableWidget(0, 0)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.cellClicked.connect(lambda r, c: self._emit(self.cellChosen, r, c))
        self.table.cellDoubleClicked.connect(lambda r, c: self._emit(self.relabelRequested, r, c))
        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel("Show"))
        top.addWidget(self.show_box)
        top.addWidget(self.legend, 1)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addLayout(top)
        lay.addWidget(self.table)

    def set_istd(self, istd) -> None:
        self.istd = istd

    def standards(self) -> list[str]:
        return [self.table.verticalHeaderItem(r).text() for r in range(self.table.rowCount())]

    def runs(self) -> list[str]:
        return [self.table.horizontalHeaderItem(c).text() for c in range(self.table.columnCount())]

    def _emit(self, signal, row: int, col: int) -> None:
        signal.emit(self.runs()[col], self.standards()[row])

    def refresh(self) -> None:
        if self.istd is None:
            return
        istd = self.istd
        names = [e.name for e in istd.library if e.name in istd.references]
        names += [n for n in istd.references if n not in names]
        runs = [e.name for e in istd.project.entries if e.name in istd.matches]
        self.table.setRowCount(len(names))
        self.table.setColumnCount(len(runs))
        self.table.setVerticalHeaderLabels(names)
        self.table.setHorizontalHeaderLabels(runs)
        show = SHOW[self.show_box.currentText()]
        for r, name in enumerate(names):
            for c, run in enumerate(runs):
                m = istd.matches.get(run, {}).get(name)
                item = self.table.item(r, c)
                if item is None:
                    item = QtWidgets.QTableWidgetItem()
                    item.setTextAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
                    self.table.setItem(r, c, item)
                if m is None:
                    item.setText("")
                    item.setToolTip("Not searched yet: press Go")
                    item.setForeground(QtGui.QBrush())
                    continue
                mark = "✓" if m.found else ("·" if m.status == NO_PEAKS else "✗")
                value = show(m)
                item.setText(f"{mark} {value}" if value else mark)
                item.setToolTip(_tooltip(m, istd.settings))
                item.setForeground(QtGui.QBrush(STATUS_COLOURS[_colour_key(m)]))
