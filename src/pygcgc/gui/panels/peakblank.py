"""➖ Subtract a run's peaks: the settings side of the step, and the peak-by-peak result.

Choose the blank run and the factor, press Go, and every other run has the blank's partner of each peak taken
away from its volume. The table shows, for the run you pick, each peak with the blank's volume, what was
removed and what is left, and says in positive terms whether the peak stands **above blank** or is **near
blank**. While the link to quantification is on, the same blank and factor make 🧪 report amounts as sample
minus blank.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from ...core.peak_blank import ABOVE, NEAR, REDUCED, REMOVED, blank_status

COLUMNS = ["#", "t₁ (min)", "t₂ (s)", "Volume", "Blank volume", "Removed", "Left", "Removed (%)", "Spectrum match",
           "Against the blank"]
GOOD, WEAK = "#009E73", "#E69F00"


class _Number(QtWidgets.QTableWidgetItem):
    def __init__(self, value: float, text: str):
        super().__init__(text)
        self.value = float(value) if value is not None and np.isfinite(value) else -np.inf

    def __lt__(self, other):
        return self.value < other.value if isinstance(other, _Number) else super().__lt__(other)


class PeakBlankPanel(QtWidgets.QWidget):
    settingsChanged = QtCore.Signal(object, bool)          # PeakBlankSettings, apply to quantification

    def __init__(self, parent=None):
        super().__init__(parent)
        self.state = None
        self._updating = False
        self.blank_box = QtWidgets.QComboBox()
        self.blank_box.setToolTip("The run whose peaks are taken away from every other run")
        self.factor = QtWidgets.QDoubleSpinBox()
        self.factor.setRange(0.0, 50.0)
        self.factor.setDecimals(3)
        self.factor.setSingleStep(0.05)
        self.factor.setValue(1.0)
        self.factor.setToolTip("How much of the blank each sample carries: 1 if it saw exactly the blank's "
                               "contamination, less if it was sampled for a shorter time or diluted")
        self.suggestion = QtWidgets.QLabel("")
        self.suggestion.setWordWrap(True)
        self.use_suggestion = QtWidgets.QPushButton("Use the median")
        self.use_suggestion.clicked.connect(self._use_suggestion)
        self.use_suggestion.setEnabled(False)
        self.apply_quant = QtWidgets.QCheckBox("Also take the blank's amounts off the quantified amounts")
        self.apply_quant.setChecked(True)
        self.apply_quant.setToolTip("Quantification then reports sample minus factor times blank, target by target, and the "
                                    "unknowns' semi-quantitative amounts from their corrected volumes")
        self.more = QtWidgets.QToolButton()
        self.more.setText("Matching windows")
        self.more.setCheckable(True)
        self.more.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.more.setArrowType(QtCore.Qt.ArrowType.RightArrow)
        self.details = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(self.details)
        form.setContentsMargins(0, 0, 0, 0)
        self.dt1 = self._spin(0.5, 120.0, 0.5, 1, " s")
        self.dt2 = self._spin(0.01, 2.0, 0.01, 2, " s")
        self.similarity = self._spin(0.0, 1.0, 0.05, 2, "")
        self.significance = self._spin(1.0, 100.0, 0.5, 1, " x")
        form.addRow("A blank peak may sit within t₁", self.dt1)
        form.addRow("and t₂", self.dt2)
        form.addRow("Spectra must match at least", self.similarity)
        form.addRow("Above blank from", self.significance)
        self.details.hide()
        self.more.toggled.connect(self._toggle_details)

        self.show_box = QtWidgets.QComboBox()
        self.show_box.currentIndexChanged.connect(self._fill_table)
        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)
        self.table = QtWidgets.QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.setMinimumHeight(180)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.factor)
        row.addWidget(self.use_suggestion)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(QtWidgets.QLabel("Blank run"))
        lay.addWidget(self.blank_box)
        lay.addWidget(QtWidgets.QLabel("Factor"))
        lay.addLayout(row)
        lay.addWidget(self.suggestion)
        lay.addWidget(self.apply_quant)
        lay.addWidget(self.more)
        lay.addWidget(self.details)
        lay.addWidget(QtWidgets.QLabel("Show run"))
        lay.addWidget(self.show_box)
        lay.addWidget(self.status)
        lay.addWidget(self.table, 1)
        for widget in (self.blank_box,):
            widget.currentIndexChanged.connect(self._edited)
        for widget in (self.factor, self.dt1, self.dt2, self.similarity, self.significance):
            widget.valueChanged.connect(self._edited)
        self.apply_quant.toggled.connect(self._edited)

    @staticmethod
    def _spin(lo, hi, step, decimals, suffix):
        w = QtWidgets.QDoubleSpinBox()
        w.setRange(lo, hi)
        w.setSingleStep(step)
        w.setDecimals(decimals)
        w.setSuffix(suffix)
        return w

    def _toggle_details(self, visible: bool) -> None:
        self.details.setVisible(visible)
        self.more.setArrowType(QtCore.Qt.ArrowType.DownArrow if visible else QtCore.Qt.ArrowType.RightArrow)

    # ------------------------------------------------------------- state
    def set_state(self, state) -> None:
        self.state = state

    def refresh_runs(self, entries) -> None:
        """Fill the run lists from the loaded runs and show the state's settings."""
        if self.state is None:
            return
        s = self.state.settings
        self._updating = True
        current_show = self.show_box.currentData()
        self.blank_box.clear()
        self.show_box.clear()
        self.blank_box.addItem("(choose the blank run)", "")
        for i, e in enumerate(entries):
            if e.run is not None:
                self.blank_box.addItem(e.name, e.name)
                self.show_box.addItem(e.name, i)
        self.blank_box.setCurrentIndex(max(self.blank_box.findData(s.blank_run), 0))
        self.show_box.setCurrentIndex(max(self.show_box.findData(current_show), 0))
        self.factor.setValue(s.factor)
        self.dt1.setValue(s.max_dt1)
        self.dt2.setValue(s.max_dt2)
        self.similarity.setValue(s.min_similarity)
        self.significance.setValue(s.significance)
        self.apply_quant.setChecked(self.state.apply_to_quant)
        self._updating = False
        self._fill_table()

    def set_shown_run(self, index: int) -> None:
        at = self.show_box.findData(index)
        if at >= 0 and at != self.show_box.currentIndex():
            self.show_box.setCurrentIndex(at)

    def current_settings(self):
        return replace(self.state.settings, blank_run=self.blank_box.currentData() or "", factor=self.factor.value(),
                       max_dt1=self.dt1.value(), max_dt2=self.dt2.value(), min_similarity=self.similarity.value(),
                       significance=self.significance.value())

    def _edited(self, *_a) -> None:
        if self._updating or self.state is None:
            return
        self.settingsChanged.emit(self.current_settings(), self.apply_quant.isChecked())

    def _use_suggestion(self) -> None:
        found = self._shown_result()
        if found is not None and np.isfinite(found.suggestion["median"]):
            self.factor.setValue(round(found.suggestion["median"], 3))

    def _shown_result(self):
        index = self.show_box.currentData()
        return None if index is None or self.state is None else self.state.results_for(int(index))

    # ------------------------------------------------------------- the table
    def refresh(self) -> None:
        self._fill_table()

    def _fill_table(self, *_a) -> None:
        if self._updating or self.state is None:
            return
        index = self.show_box.currentData()
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        self.use_suggestion.setEnabled(False)
        if index is None:
            self.status.setText("Open some runs first.")
            self.suggestion.setText("")
            return
        index = int(index)
        problems = self.state.why_not(index)
        found = self.state.results_for(index)
        warning = self.state.phase_note(index)
        if found is None:
            self.status.setText(" ".join(problems + ([warning] if warning else [])) if (problems or warning)
                                else "Press Subtract to take the blank off this run.")
            self.suggestion.setText("")
            return
        current = self.state.current_for(index) is not None
        s = found.summary
        self.status.setText(
            f"{s['peaks']} peaks: {s['above']} above blank, {s['near']} near blank; {s['matched']} have a partner in the blank, "
            f"{s['removed']} removed entirely, {s['reduced']} reduced; {s['signal_removed']:.1%} of the run's signal removed."
            + ("" if current else " <out of date: the peaks, alignment or settings have changed since>")
            + (" " + warning if warning else ""))
        hint = found.suggestion
        if hint["n"]:
            self.suggestion.setText(
                f"Shared peaks: {hint['n']}. Sample-to-blank volume ratio: median {hint['median']:.2f}, middle half "
                f"{hint['q1']:.2f} to {hint['q3']:.2f}. Only the peaks that are pure contamination show the true factor; "
                "real compounds also in the blank push it up, so the lower end is the safer guide.")
            self.use_suggestion.setEnabled(True)
        else:
            self.suggestion.setText("No peak was shared with the blank, so there is nothing to suggest a factor from.")
        rows = found.rows
        self.table.setRowCount(len(rows))
        for r, row in enumerate(rows):
            status = blank_status(row, self.state.settings)
            words = ("● above blank" if status == ABOVE else "○ near blank")
            if row.blank_index is None:
                words += " (not in the blank)"
            values = [(str(row.index + 1), row.index), (f"{row.t1 / 60:.3f}", row.t1), (f"{row.t2:.3f}", row.t2),
                      (f"{row.volume:.4g}", row.volume), (f"{row.blank_volume:.4g}" if row.blank_index is not None else "", row.blank_volume),
                      (f"{row.removed:.4g}" if row.removed else "", row.removed), (f"{row.corrected:.4g}", row.corrected),
                      (f"{100 * row.fraction_removed:.0f}" if row.removed else "", row.fraction_removed),
                      (f"{row.similarity:.2f}" if np.isfinite(row.similarity) else "", row.similarity), (words, None)]
            for c, (text, key) in enumerate(values):
                item = _Number(key, text) if key is not None else QtWidgets.QTableWidgetItem(text)
                if c == 9:
                    item.setForeground(QtGui.QBrush(QtGui.QColor(GOOD if status == ABOVE else WEAK)))
                if row.match in (REMOVED, REDUCED) and c in (5, 6):
                    item.setForeground(QtGui.QBrush(QtGui.QColor(WEAK if row.match == REMOVED else GOOD)))
                self.table.setItem(r, c, item)
        self.table.setSortingEnabled(True)
        self.table.sortItems(3, QtCore.Qt.SortOrder.DescendingOrder)

    def set_background(self, colour) -> None:      # the table needs no plot background
        pass
