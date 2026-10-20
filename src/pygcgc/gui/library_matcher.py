"""📚 Library matcher: name the detected peaks from the spectrum and the retention index.

Read top to bottom:

1. **Set-up strip.** Which library (the NIST folder, or an MSP file), which column class (5 % phenyl is
   "semi-standard non-polar"), which run or runs, and where the alkane ladder comes from. Anything missing is
   one red sentence saying what to do. Column dimensions and the oven programme are *not* needed: a
   retention index is relative to the alkanes run on the same method.
2. **Largest peaks.** The ten peaks with the most fitted signal, each with its share of the run, its name and
   how sure the match is, so the peaks that matter are never buried in a long table.
3. **All peaks**, largest first, with filters. Confidence is a symbol and a word (never colour alone).
4. **The selected peak.** The measured spectrum against the library's, head to tail; the top candidates; the
   retention index of the peak against each candidate's; and every uncertainty as a sentence.

The matching is `core/libmatch.py`. A decision made here (accept, reject, choose another candidate, mark
unknown) is saved with the project but never raises the confidence shown beside it: a human decision is
recorded, not turned into evidence.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from ..core.libmatch import FLAGS, describe_candidate_rank, largest_peaks, write_identifications_csv
from ..core.nist_sdf import COLUMN_CLASSES
from ..core.identify import promote_to_target
from ..core.targets import write_targets_csv
from .match_state import NAMING_ACTIONS
from .theme import fit_to_screen

LEVELS = {"confident": ("●", "#009E73", "Confident"), "probable": ("◐", "#0072B2", "Probable"),
          "tentative": ("○", "#E69F00", "Tentative"), "unknown": ("✕", "#777777", "Unknown")}
COLUMNS = ["#", "t₁ (min)", "t₂ (s)", "% signal", "Name", "Formula", "Confidence", "MF", "RMF", "RI peak",
           "RI library", "z", "Flags", "Decision", "Against the blank"]
CANDIDATE_COLUMNS = ["Name", "Formula", "MF", "RMF", "RI library", "Source", "z", "RI gate"]
COLUMN_LABELS = {"Standard non-polar": "Standard non-polar (100 % dimethyl)",
                 "Semi-standard non-polar": "Semi-standard non-polar (5 % phenyl)",
                 "Standard polar": "Standard polar (wax)"}
UP, DOWN, WARN = (0, 114, 178), (230, 159, 0), "#D55E00"


class _Number(QtWidgets.QTableWidgetItem):
    """A cell shown formatted but sorted by value."""

    def __init__(self, value: float, text: str):
        super().__init__(text)
        self.value = float(value) if value is not None and np.isfinite(value) else -np.inf

    def __lt__(self, other):
        return self.value < other.value if isinstance(other, _Number) else super().__lt__(other)


def level_text(level: str) -> str:
    symbol, _colour, word = LEVELS[level]
    return f"{symbol} {word}"


class LibraryMatcherWindow(QtWidgets.QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.setWindowTitle("📚 Library matcher")
        self.setWindowFlag(QtCore.Qt.WindowType.WindowMaximizeButtonHint, True)
        fit_to_screen(self, 1250, 860)
        self.mw = main_window
        self.project = main_window.project
        self.state = main_window.project.matcher
        self._rows: list = []                 # PeakMatch per results-table row
        self._shown = None                    # the PeakMatch in the detail pane
        self._updating = False

        # ---- set-up strip
        self.library_label = QtWidgets.QLabel("")
        self.library_label.setWordWrap(True)
        nist = QtWidgets.QPushButton("NIST folder…")
        nist.setToolTip("The folder holding mainlib.sdf (and replib.sdf, nist_ri.sdf), or its SDF subfolder")
        nist.clicked.connect(self._choose_nist)
        msp = QtWidgets.QPushButton("MSP file…")
        msp.clicked.connect(self._choose_msp)
        self.column_box = QtWidgets.QComboBox()
        for key in COLUMN_CLASSES:
            self.column_box.addItem(COLUMN_LABELS[key], key)
        self.column_box.setToolTip("Literature retention indices differ by hundreds of units between column classes. "
                                   "Nothing else about the column, and not the oven programme, is needed.")
        self.column_box.currentIndexChanged.connect(self._column_changed)
        self.run_box = QtWidgets.QComboBox()
        self.run_box.currentIndexChanged.connect(self._refresh_all)
        self.all_runs = QtWidgets.QCheckBox("All 70 eV runs")
        self.ladder_box = QtWidgets.QComboBox()
        self.ladder_box.setToolTip("Each run's own alkane ladder, or borrow another run's (through the alignment if "
                                   "both runs are aligned)")
        self.ladder_box.currentIndexChanged.connect(self._ladder_changed)
        self.match_btn = QtWidgets.QPushButton("📚 Match peaks")
        self.match_btn.setObjectName("go")
        self.match_btn.clicked.connect(self._match)
        settings_btn = QtWidgets.QPushButton("⚙️ Thresholds…")
        settings_btn.clicked.connect(self._edit_settings)
        self.requirements = QtWidgets.QLabel("")
        self.requirements.setWordWrap(True)
        self.requirements.setTextFormat(QtCore.Qt.TextFormat.RichText)

        row1 = QtWidgets.QHBoxLayout()
        row1.addWidget(QtWidgets.QLabel("Library"))
        row1.addWidget(self.library_label, 1)
        row1.addWidget(nist)
        row1.addWidget(msp)
        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(QtWidgets.QLabel("Column"))
        row2.addWidget(self.column_box)
        row2.addSpacing(10)
        row2.addWidget(QtWidgets.QLabel("Run"))
        row2.addWidget(self.run_box, 1)
        row2.addWidget(self.all_runs)
        row2.addSpacing(10)
        row2.addWidget(QtWidgets.QLabel("Alkane ladder"))
        row2.addWidget(self.ladder_box)
        row3 = QtWidgets.QHBoxLayout()
        row3.addWidget(self.match_btn)
        row3.addWidget(settings_btn)
        row3.addStretch(1)

        # ---- largest peaks
        self.summary = QtWidgets.QTableWidget(0, 5)
        self.summary.setHorizontalHeaderLabels(["Rank", "Share of the run", "t₁ (min)", "Best match", "Confidence"])
        self.summary.verticalHeader().setVisible(False)
        self.summary.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.summary.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.summary.horizontalHeader().setSectionResizeMode(3, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.summary.setMaximumHeight(230)
        self.summary.itemSelectionChanged.connect(self._summary_selected)

        # ---- all peaks
        self.level_filter = QtWidgets.QComboBox()
        for key, text in (("all", "All peaks"), ("good", "Confident and probable"), ("attention", "Needs attention"),
                          ("unknown", "Unknown only")):
            self.level_filter.addItem(text, key)
        self.flag_filter = QtWidgets.QCheckBox("Flagged only")
        self.name_filter = QtWidgets.QLineEdit()
        self.name_filter.setPlaceholderText("Filter by name or formula")
        self.min_share = QtWidgets.QDoubleSpinBox()
        self.min_share.setRange(0, 100)
        self.min_share.setDecimals(2)
        self.min_share.setSuffix(" % or more")
        for widget in (self.level_filter, self.flag_filter, self.name_filter, self.min_share):
            signal = (widget.currentIndexChanged if isinstance(widget, QtWidgets.QComboBox) else
                      widget.toggled if isinstance(widget, QtWidgets.QCheckBox) else
                      widget.textChanged if isinstance(widget, QtWidgets.QLineEdit) else widget.valueChanged)
            signal.connect(self._fill_results)
        filters = QtWidgets.QHBoxLayout()
        filters.addWidget(self.level_filter)
        filters.addWidget(self.flag_filter)
        filters.addWidget(self.name_filter, 1)
        filters.addWidget(QtWidgets.QLabel("Signal"))
        filters.addWidget(self.min_share)
        self.table = QtWidgets.QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(4, QtWidgets.QHeaderView.ResizeMode.Interactive)
        self.table.setColumnWidth(4, 220)
        self.table.itemSelectionChanged.connect(self._row_selected)
        self.table_note = QtWidgets.QLabel("")
        self.table_note.setWordWrap(True)

        # ---- selected peak
        self.detail_title = QtWidgets.QLabel("Select a peak")
        self.detail_title.setWordWrap(True)
        self.detail_title.setStyleSheet("font-size: 12pt; font-weight: bold;")
        self.detail_text = QtWidgets.QLabel("")
        self.detail_text.setWordWrap(True)
        self.detail_text.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.detail_text.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)
        self.spectrum_plot = pg.PlotWidget()
        self.spectrum_plot.setMenuEnabled(False)
        self.spectrum_plot.hideButtons()
        self.spectrum_plot.setLabel("bottom", "m/z")
        self.spectrum_plot.setLabel("left", "Relative intensity (%): peak above, library below")
        self.spectrum_plot.setMinimumHeight(190)
        self.spectrum_plot.addLine(y=0, pen=pg.mkPen((128, 128, 128)))
        self.up_bars = pg.BarGraphItem(x=[], height=[], width=0.8, brush=UP)
        self.down_bars = pg.BarGraphItem(x=[], height=[], width=0.8, brush=DOWN)
        self.spectrum_plot.addItem(self.up_bars)
        self.spectrum_plot.addItem(self.down_bars)
        self.ion_labels: list = []
        self.candidates = QtWidgets.QTableWidget(0, len(CANDIDATE_COLUMNS))
        self.candidates.setHorizontalHeaderLabels(CANDIDATE_COLUMNS)
        self.candidates.verticalHeader().setVisible(False)
        self.candidates.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.candidates.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.candidates.setMaximumHeight(170)
        self.candidates.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.candidates.itemSelectionChanged.connect(self._candidate_selected)
        self.ri_plot = pg.PlotWidget()
        self.ri_plot.setMenuEnabled(False)
        self.ri_plot.hideButtons()
        self.ri_plot.setLabel("bottom", "Retention index")
        self.ri_plot.hideAxis("left")
        self.ri_plot.setMaximumHeight(120)
        self.measured_band = pg.LinearRegionItem(values=(0, 0), movable=False, brush=pg.mkBrush(0, 114, 178, 60))
        self.ri_plot.addItem(self.measured_band)
        self.ri_errors = pg.ErrorBarItem(x=[], y=[], left=[], right=[], beam=0.25, pen=pg.mkPen(DOWN, width=2))
        self.ri_points = pg.ScatterPlotItem(size=9, brush=pg.mkBrush(DOWN), pen=pg.mkPen("k"))
        self.ri_plot.addItem(self.ri_errors)
        self.ri_plot.addItem(self.ri_points)

        self.accept_btn = QtWidgets.QPushButton("✔ Accept best")
        self.choose_btn = QtWidgets.QPushButton("Use selected candidate")
        self.reject_btn = QtWidgets.QPushButton("✘ Reject")
        self.rename_btn = QtWidgets.QPushButton("Suggest another name…")
        self.rename_btn.setToolTip("Type a different name for this peak; capitals, commas, dashes and spaces are checked against the library")
        self.unknown_btn = QtWidgets.QPushButton("Mark unknown")
        self.clear_btn = QtWidgets.QPushButton("Clear decision")
        self.promote_btn = QtWidgets.QPushButton("Promote to target…")
        self.export_btn = QtWidgets.QPushButton("Export identifications…")
        self.accept_btn.clicked.connect(lambda: self._decide("accepted"))
        self.choose_btn.clicked.connect(lambda: self._decide("chosen"))
        self.reject_btn.clicked.connect(lambda: self._decide("rejected"))
        self.unknown_btn.clicked.connect(lambda: self._decide("unknown"))
        self.rename_btn.clicked.connect(self._rename)
        self.clear_btn.clicked.connect(self._clear_decision)
        self.promote_btn.clicked.connect(self._promote)
        self.export_btn.clicked.connect(self._export)
        buttons = QtWidgets.QHBoxLayout()
        for b in (self.accept_btn, self.choose_btn, self.rename_btn, self.reject_btn, self.unknown_btn, self.clear_btn):
            buttons.addWidget(b)
        buttons.addStretch(1)
        buttons.addWidget(self.promote_btn)
        buttons.addWidget(self.export_btn)

        left = QtWidgets.QWidget()
        left_lay = QtWidgets.QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 0, 0)
        left_lay.addWidget(QtWidgets.QLabel("<b>Largest peaks</b>"))
        left_lay.addWidget(self.summary)
        left_lay.addWidget(QtWidgets.QLabel("<b>All peaks</b>, largest first"))
        left_lay.addLayout(filters)
        left_lay.addWidget(self.table, 1)
        left_lay.addWidget(self.table_note)
        right = QtWidgets.QWidget()
        right_lay = QtWidgets.QVBoxLayout(right)
        right_lay.setContentsMargins(0, 0, 0, 0)
        right_lay.addWidget(self.detail_title)
        right_lay.addWidget(self.detail_text)
        right_lay.addWidget(self.spectrum_plot, 1)
        right_lay.addWidget(self.candidates)
        right_lay.addWidget(self.ri_plot)
        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        splitter.addWidget(left)
        splitter.addWidget(right)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)

        lay = QtWidgets.QVBoxLayout(self)
        lay.addLayout(row1)
        lay.addLayout(row2)
        lay.addLayout(row3)
        lay.addWidget(self.requirements)
        lay.addWidget(splitter, 1)
        lay.addLayout(buttons)

        self.state.changed.connect(self._refresh_all)
        self.state.libraryChanged.connect(self._refresh_all)
        self.state.runningChanged.connect(self._refresh_all)
        self.project.runsChanged.connect(self._fill_runs)
        self.project.blobs.changed.connect(self._refresh_all)
        self.project.hydrocarbons.changed.connect(self._fill_ladders)
        self._fill_runs()

    # ================================================================ opening
    def open_window(self, index: int | None = None) -> None:
        at = self.run_box.findData(index) if index is not None else -1
        if at >= 0:
            self.run_box.setCurrentIndex(at)
        if self.state.source and self.state.library is None and not self.state.loading:
            self.state.load_library(self.state.source)       # quick when the library is cached
        self._fill_ladders()
        self._refresh_all()
        self.show()
        self.raise_()
        self.activateWindow()

    def set_background(self, colour) -> None:
        self.spectrum_plot.setBackground(colour)
        self.ri_plot.setBackground(colour)

    # ================================================================ set-up
    def _fill_runs(self) -> None:
        current = self.run_box.currentData()
        self._updating = True
        self.run_box.clear()
        for i in self.project.loaded():
            self.run_box.addItem(self.project.entries[i].name, i)
        at = self.run_box.findData(current if current is not None else self.mw.index)
        self.run_box.setCurrentIndex(max(at, 0))
        self._updating = False
        self._fill_ladders()
        self._refresh_all()

    def _fill_ladders(self) -> None:
        self._updating = True
        self.ladder_box.clear()
        self.ladder_box.addItem("Each run's own", "")
        for i in self.project.loaded():
            series = self.project.hydrocarbons.series_for(i)
            if series is not None and len(series) >= 2:
                self.ladder_box.addItem(f"Borrow from {self.project.entries[i].name}", self.project.entries[i].name)
        at = self.ladder_box.findData(self.state.ladder_run)
        self.ladder_box.setCurrentIndex(max(at, 0))
        self._updating = False
        key = self.column_box.findData(self.state.settings.column_class)
        self.column_box.blockSignals(True)
        self.column_box.setCurrentIndex(max(key, 0))
        self.column_box.blockSignals(False)

    def _choose_nist(self) -> None:
        start = self.state.source if self.state.source else ""
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "NIST library folder (holding mainlib.sdf)", start)
        if folder:
            self.state.load_library(folder)

    def _choose_msp(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load an MSP library", "", "MSP library (*.msp)")
        if path:
            self.state.load_library(path)

    def _column_changed(self) -> None:
        if self._updating:
            return
        from dataclasses import replace
        self.state.set_settings(replace(self.state.settings, column_class=self.column_box.currentData()))

    def _ladder_changed(self) -> None:
        if not self._updating:
            self.state.set_ladder_run(self.ladder_box.currentData() or "")

    def _edit_settings(self) -> None:
        from .dialogs import SettingsDialog
        tips = {"ri_gate_z": "A library compound whose RI is further than this many standard deviations from the peak is demoted",
                "sigma_systematic": "RI units allowed for a GC x GC first dimension against one-dimensional literature columns",
                "mf_confident": "Provisional. Dot-product match factors resemble NIST MS Search's but are not identical",
                "exclude_ions": "Ions (column bleed, for example) left out of both spectra before comparing"}
        dialog = SettingsDialog("📚 Matching thresholds (provisional)", self.state.settings, tips=tips, parent=self)
        if dialog.exec():
            new = dialog.settings()
            if new.column_class not in COLUMN_CLASSES:
                QtWidgets.QMessageBox.warning(self, "Column class", "Column class must be one of: " + "; ".join(COLUMN_CLASSES))
                return
            self.state.set_settings(new)
            self._fill_ladders()

    def _indices(self) -> list[int]:
        if self.all_runs.isChecked():
            return [i for i in self.project.loaded()
                    if ((self.project.entries[i].import_meta or {}).get("energy_ev") in (None, 70))]
        index = self.run_box.currentData()
        return [] if index is None else [int(index)]

    def _match(self) -> None:
        self.state.match_runs(self._indices())

    # ================================================================ refresh
    def _refresh_all(self, *_a) -> None:
        if self._updating:
            return
        state = self.state
        index = self.run_box.currentData()
        self.library_label.setText(state.library_status)
        self.library_label.setStyleSheet(f"color: {WARN};" if state.library is None else "")
        lines: list[tuple[str, bool]] = []                 # (sentence, is it a problem?)
        missing: list = []
        if index is None:
            lines.append(("Open and select a run.", True))
        else:
            missing = state.why_not(int(index))
            lines += [(sentence, True) for sentence in missing]
            series, _borrowed, _map, note = state.ladder_for(int(index))
            lines.append((("Ladder: " + note) if series is not None else note.capitalize() + ".", series is None))
            if state.results_for(int(index)) is not None and state.current_for(int(index)) is None:
                lines.append(("The peaks, library or settings have changed since this run was matched: match again.", True))
        self.requirements.setText("<br>".join(f"<b style='color:{WARN}'>{text}</b>" if bad else text for text, bad in lines))
        ready = index is not None and not (missing if not self.all_runs.isChecked() else [])
        self.match_btn.setEnabled(bool(ready) and not state.running and not state.loading)
        self.match_btn.setText("Matching…" if state.running else "📚 Match peaks")
        self._fill_summary()
        self._fill_results()

    def _current_matches(self):
        index = self.run_box.currentData()
        return None if index is None else self.state.results_for(int(index))

    def _decisions(self) -> dict:
        index = self.run_box.currentData()
        return {} if index is None else self.state.decisions_for(int(index))

    def _fill_summary(self) -> None:
        found = self._current_matches()
        top = largest_peaks(found.matches, 10) if found is not None else []
        decisions = self._decisions()
        self.summary.setRowCount(len(top))
        biggest = max((s for _m, s in top if np.isfinite(s)), default=1.0) or 1.0
        for r, (match, share) in enumerate(top):
            bar = "█" * max(1, int(round(14 * share / biggest))) if np.isfinite(share) else ""
            best = match.best
            name = self._display_name(match, decisions)
            cells = (str(r + 1), f"{bar}  {share:.1%}" if np.isfinite(share) else "", f"{match.query.t1 / 60:.3f}",
                     name if best else "(no candidate)", level_text(match.level))
            for c, text in enumerate(cells):
                item = QtWidgets.QTableWidgetItem(text)
                item.setData(QtCore.Qt.ItemDataRole.UserRole, match.query.label)
                if c == 4:
                    item.setForeground(QtGui.QBrush(QtGui.QColor(LEVELS[match.level][1])))
                self.summary.setItem(r, c, item)
        self.summary.resizeColumnToContents(1)

    @staticmethod
    def _display_name(match, decisions) -> str:
        d = decisions.get(match.query.label)
        best = match.best
        if d and d["action"] in NAMING_ACTIONS and d.get("name"):
            return d["name"]
        return best.name if best else ""

    def _passes(self, match, share, decisions) -> bool:
        mode = self.level_filter.currentData()
        if mode == "good" and match.level not in ("confident", "probable"):
            return False
        if mode == "unknown" and match.level != "unknown":
            return False
        if mode == "attention" and not (match.flags or match.level in ("tentative", "unknown")):
            return False
        if self.flag_filter.isChecked() and not match.flags:
            return False
        text = self.name_filter.text().strip().lower()
        best = match.best
        if text and not (best and (text in best.name.lower() or text in best.formula.lower())):
            return False
        return not (np.isfinite(share) and 100 * share < self.min_share.value())

    def _fill_results(self, *_a) -> None:
        found = self._current_matches()
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        self._rows = []
        if found is None:
            self.table_note.setText("No matches yet. Press Match peaks.")
            self._show_detail(None)
            return
        decisions = self._decisions()
        blank = self.state.project.blank.current_for(int(self.run_box.currentData()))
        position = {id(m): k for k, m in enumerate(found.matches)}
        total = sum(max(m.query.volume, 0.0) for m in found.matches)
        ranked = sorted(found.matches, key=lambda m: m.query.volume, reverse=True)
        keep = []
        for rank, m in enumerate(ranked, start=1):
            share = m.query.volume / total if total > 0 else float("nan")
            if self._passes(m, share, decisions):
                keep.append((rank, m, share))
        self.table.setRowCount(len(keep))
        for r, (rank, m, share) in enumerate(keep):
            best = m.best
            d = decisions.get(m.query.label)
            values = [(str(rank), rank), (f"{m.query.t1 / 60:.3f}", m.query.t1), (f"{m.query.t2:.3f}", m.query.t2),
                      (f"{100 * share:.2f}" if np.isfinite(share) else "", share),
                      (self._display_name(m, decisions), None), (best.formula if best else "", None),
                      (level_text(m.level), ("confident", "probable", "tentative", "unknown").index(m.level)),
                      (f"{best.mf:.0f}" if best else "", best.mf if best else None),
                      (f"{best.rmf:.0f}" if best else "", best.rmf if best else None),
                      (f"{m.query.ri:.0f} ± {m.query.ri_sigma:.0f}" if np.isfinite(m.query.ri) else "none", m.query.ri),
                      (f"{best.ri_lib:.0f}" if best and np.isfinite(best.ri_lib) else "", best.ri_lib if best else None),
                      (f"{best.z:+.1f}" if best and np.isfinite(best.z) else "", best.z if best else None),
                      (f"{len(m.flags)}" if m.flags else "", len(m.flags)),
                      (f"{d['action']}: {d['name']}" if d and d.get("name") else (d["action"] if d else ""), None),
                      (self._blank_text(blank, position.get(id(m))), None)]
            for c, (text, key) in enumerate(values):
                item = _Number(key, text) if key is not None else QtWidgets.QTableWidgetItem(text)
                if c == 6:
                    item.setForeground(QtGui.QBrush(QtGui.QColor(LEVELS[m.level][1])))
                if c == 12 and m.flags:
                    item.setForeground(QtGui.QBrush(QtGui.QColor(WARN)))
                    item.setToolTip("\n".join(FLAGS.get(f, f) for f in m.flags))
                item.setData(QtCore.Qt.ItemDataRole.UserRole, m.query.label)
                self.table.setItem(r, c, item)
            self._rows.append(m)
        self.table.setSortingEnabled(True)
        self.table.sortItems(3, QtCore.Qt.SortOrder.DescendingOrder)
        counts = {lv: sum(m.level == lv for m in found.matches) for lv in LEVELS}
        age = "" if self.state.current_for(int(self.run_box.currentData())) is not None else "  (out of date)"
        self.table_note.setText(
            f"{len(keep)} of {len(found.matches)} peaks shown. " + ", ".join(f"{counts[lv]} {lv}" for lv in LEVELS)
            + f". {found.ladder_note}.{age}")
        if self.table.rowCount():
            if self._shown is None or not any(self._shown.query.label == m.query.label for m in self._rows):
                self.table.selectRow(0)
            else:
                self._select_label(self._shown.query.label)
        else:
            self._show_detail(None)

    def _blank_text(self, blank, k) -> str:
        """Where the peak stands against the subtracted blank run, in positive terms; empty if no subtraction was made."""
        if blank is None or k is None or k >= len(blank.rows):
            return ""
        from ..core.peak_blank import blank_status
        row = blank.rows[k]
        text = ("● " if blank_status(row, self.project.blank.settings) == "above blank" else "○ ") + blank_status(row, self.project.blank.settings)
        return text + (f", {100 * row.fraction_removed:.0f}% removed" if row.removed else "")

    def _match_for_label(self, label):
        found = self._current_matches()
        return next((m for m in found.matches if m.query.label == label), None) if found else None

    def _select_label(self, label) -> None:
        for r in range(self.table.rowCount()):
            item = self.table.item(r, 0)
            if item is not None and item.data(QtCore.Qt.ItemDataRole.UserRole) == label:
                self.table.selectRow(r)
                return

    def _summary_selected(self) -> None:
        rows = sorted({i.row() for i in self.summary.selectedIndexes()})
        if rows:
            item = self.summary.item(rows[0], 0)
            if item is not None:
                self._select_label(item.data(QtCore.Qt.ItemDataRole.UserRole))

    def _row_selected(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if not rows:
            return
        item = self.table.item(rows[0], 0)
        self._show_detail(self._match_for_label(item.data(QtCore.Qt.ItemDataRole.UserRole)))

    # ================================================================ detail
    def _show_detail(self, match, candidate_row: int = 0) -> None:
        self._shown = match
        self.candidates.blockSignals(True)
        self.candidates.setRowCount(0)
        self.candidates.blockSignals(False)
        for item in self.ion_labels:
            self.spectrum_plot.removeItem(item)
        self.ion_labels = []
        self.up_bars.setOpts(x=[], height=[])
        self.down_bars.setOpts(x=[], height=[])
        self.ri_errors.setData(x=np.array([]), y=np.array([]), left=np.array([]), right=np.array([]))
        self.ri_points.setData([], [])
        for b in (self.accept_btn, self.choose_btn, self.rename_btn, self.reject_btn, self.unknown_btn, self.clear_btn,
                  self.promote_btn):
            b.setEnabled(match is not None and match.best is not None)
        self.export_btn.setEnabled(self._current_matches() is not None)
        if match is None or match.best is None:
            self.detail_title.setText("Select a peak" if match is None else "No candidate for this peak")
            self.detail_text.setText("" if match is None else self._flag_html(match))
            return
        best = match.best
        q = match.query
        self.detail_title.setText(f"{best.name}   {best.formula}   CAS {best.cas or 'unknown'}")
        self.candidates.blockSignals(True)
        self.candidates.setRowCount(len(match.candidates))
        for r, c in enumerate(match.candidates):
            values = (c.name, c.formula, f"{c.mf:.0f}", f"{c.rmf:.0f}",
                      f"{c.ri_lib:.0f} ± {c.ri_sigma_lib:.0f}" if np.isfinite(c.ri_lib) else "none",
                      c.ri_kind + (f" (n={c.n_ri})" if c.n_ri else ""), f"{c.z:+.1f}" if np.isfinite(c.z) else "",
                      "inside" if c.in_gate else "OUTSIDE")
            for k, text in enumerate(values):
                item = QtWidgets.QTableWidgetItem(text)
                if not c.in_gate:
                    item.setForeground(QtGui.QBrush(QtGui.QColor(WARN)))
                item.setToolTip(describe_candidate_rank(match, r))
                self.candidates.setItem(r, k, item)
        self.candidates.selectRow(min(candidate_row, len(match.candidates) - 1))
        self.candidates.blockSignals(False)
        self._draw_candidate(match, min(candidate_row, len(match.candidates) - 1))
        self.detail_text.setText(
            f"<b>{level_text(match.level)}</b> (Schymanski level {match.tier or '3 or lower'}). "
            f"{self._evidence_html(match)}<br>{self._flag_html(match)}"
            + (f"<br>Decision: <b>{self._decision_text(match)}</b>" if self._decision_text(match) else ""))

    def _evidence_html(self, match) -> str:
        best, q = match.best, match.query
        parts = [f"Spectrum: match factor {best.mf:.0f}, reverse {best.rmf:.0f}."]
        if np.isfinite(q.ri):
            lib = (f"library {best.ri_lib:.0f} ± {best.ri_sigma_lib:.0f} ({best.ri_kind}" + (f", {best.n_ri} values" if best.n_ri else "")
                   + ")") if np.isfinite(best.ri_lib) else "no library RI"
            z = f", z = {best.z:+.1f}" if np.isfinite(best.z) else ""
            parts.append(f"Retention index: peak {q.ri:.0f} ± {q.ri_sigma:.0f}, {lib}{z}.")
        else:
            parts.append("Retention index: not available.")
        if np.isfinite(match.q_value):
            parts.append(f"Estimated false discovery rate {match.q_value:.1%}.")
        return " ".join(parts)

    @staticmethod
    def _flag_html(match) -> str:
        if not match.flags:
            return "<span>No uncertainties flagged.</span>"
        items = "".join(f"<li><b style='color:{WARN}'>{f}</b>: {FLAGS.get(f, f)}</li>" for f in match.flags)
        return f"<b>Uncertainties:</b><ul style='margin:2px 0'>{items}</ul>"

    def _decision_text(self, match) -> str:
        d = self._decisions().get(match.query.label)
        return "" if not d else (f"{d['action']}: {d['name']}" if d.get("name") else d["action"])

    def _candidate_selected(self) -> None:
        rows = sorted({i.row() for i in self.candidates.selectedIndexes()})
        if rows and self._shown is not None and rows[0] < len(self._shown.candidates):
            self._draw_candidate(self._shown, rows[0])

    def _draw_candidate(self, match, k: int) -> None:
        q, c = match.query, match.candidates[k]
        entry = self.state.library.entry(c.index) if self.state.library is not None else None
        for item in self.ion_labels:
            self.spectrum_plot.removeItem(item)
        self.ion_labels = []
        keep = np.asarray(q.intensity) > 0
        mz_q, i_q = np.asarray(q.mz)[keep], np.asarray(q.intensity)[keep]
        i_q = 100.0 * i_q / i_q.max() if i_q.size else i_q
        self.up_bars.setOpts(x=mz_q, height=i_q, width=0.8)
        if entry is not None and entry.mz.size:
            i_l = 100.0 * entry.intensity / entry.intensity.max()
            self.down_bars.setOpts(x=entry.mz, height=-i_l, width=0.8)
            for mz, h, sign in [(m, h, 1) for m, h in sorted(zip(mz_q, i_q), key=lambda t: -t[1])[:8]] + \
                               [(m, h, -1) for m, h in sorted(zip(entry.mz, i_l), key=lambda t: -t[1])[:8]]:
                label = pg.TextItem(str(int(mz)), anchor=(0.5, 1.0 if sign > 0 else 0.0), color=(90, 90, 90))
                label.setPos(float(mz), sign * float(h))
                self.spectrum_plot.addItem(label)
                self.ion_labels.append(label)
        self.spectrum_plot.setYRange(-115, 115)
        # Retention index of the peak against each candidate's, in standard deviations of the combined uncertainty.
        self.measured_band.setVisible(bool(np.isfinite(q.ri)))
        xs, ys, sig = [], [], []
        if np.isfinite(q.ri):
            self.measured_band.setRegion((q.ri - q.ri_sigma, q.ri + q.ri_sigma))
        for r, cand in enumerate(match.candidates):
            if np.isfinite(cand.ri_lib):
                xs.append(cand.ri_lib)
                ys.append(len(match.candidates) - r)
                sig.append(cand.ri_sigma_lib)
        self.ri_errors.setData(x=np.array(xs), y=np.array(ys, dtype=float), left=np.array(sig), right=np.array(sig))
        self.ri_points.setData(xs, ys)
        if xs or np.isfinite(q.ri):
            allx = xs + ([q.ri] if np.isfinite(q.ri) else [])
            spread = max(max(sig, default=0), q.ri_sigma if np.isfinite(q.ri_sigma) else 0, 30)
            self.ri_plot.setXRange(min(allx) - 2 * spread, max(allx) + 2 * spread)
            self.ri_plot.setYRange(0, len(match.candidates) + 1)

    # ================================================================ decisions and files
    def _selected_candidate(self):
        rows = sorted({i.row() for i in self.candidates.selectedIndexes()})
        if self._shown is None or not rows or rows[0] >= len(self._shown.candidates):
            return None
        return self._shown.candidates[rows[0]]

    def _decide(self, action: str) -> None:
        match, index = self._shown, self.run_box.currentData()
        if match is None or index is None:
            return
        cand = self._selected_candidate() if action == "chosen" else match.best
        if action in ("accepted", "chosen") and cand is None:
            return
        name, cas = (cand.name, cand.cas) if action in ("accepted", "chosen") else ("", "")
        self.state.decide(int(index), match.query.label, action, name, cas)

    def _rename(self) -> None:
        index = self.run_box.currentData()
        if self._shown is not None and index is not None:
            from .peak_naming import rename_peak
            rename_peak(self, self.project, int(index), self._shown.query.label)

    def _clear_decision(self) -> None:
        if self._shown is not None and self.run_box.currentData() is not None:
            self.state.clear_decision(int(self.run_box.currentData()), self._shown.query.label)

    def _promote(self) -> None:
        match = self._shown
        d = self._decisions().get(match.query.label) if match else None
        cand = self._selected_candidate() or (match.best if match else None)
        if match is None or cand is None or self.state.library is None:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Append to a target table", str(self.project.targets_path) if self.project.targets_path else "target_table.csv",
            "CSV (*.csv)")
        if path:
            self.promote_to_file(match, cand, path, d)

    def promote_to_file(self, match, cand, path, decision=None) -> None:
        """The write, separate from the file dialog. The target sits where the peak was found in this run."""
        from pathlib import Path
        entry = self.state.library.entry(cand.index)
        target = promote_to_target(entry, SimpleNamespace(ri=None if not np.isfinite(match.query.ri) else float(match.query.ri)),
                                   match.query.t1, match.query.t2)
        append = Path(path).exists()
        try:
            write_targets_csv([target], path, append=append)
        except (OSError, ValueError) as exc:
            self.table_note.setText(f"Could not write the target table: {exc}")
            return
        self.table_note.setText(f"{'Appended' if append else 'Wrote'} {target.compound} to {path}. Reload the target table to use it.")

    def _export(self) -> None:
        index, found = self.run_box.currentData(), self._current_matches()
        if index is None or found is None:
            return
        from pathlib import Path
        name = self.project.entries[int(index)].name
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export identifications", f"{Path(name).stem}_identifications.csv",
                                                        "CSV (*.csv)")
        if path:
            self.export_to(path)

    def export_to(self, path) -> None:
        index, found = int(self.run_box.currentData()), self._current_matches()
        write_identifications_csv(path, self.project.entries[index].name, found.matches, self._decisions(),
                                  self.state.library.label if self.state.library is not None else "", self.state.settings)
        self.table_note.setText(f"Wrote {len(found.matches)} peaks to {path}")
