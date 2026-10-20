"""📈 Calibration builder: a separate window for a set of standard runs.

Four tabs, worked left to right, each saying what to do on it:

1. Standards: the runs selected in the main window (or a folder, or files),
   with each run's level amount (read from its name where it says, as in
   "L6_5ng"). Loading cleans them (🧹), finds their blobs (🔍) and lines
   their phases up, with no button to press.
2. Compounds: 🔎 suggested from the standards alone, loaded from a table, or
   clicked on the chromatogram. Tick IS for internal standards. Each
   compound is matched to a blob in every run and divided by its internal
   standard as soon as there is one.
3. Curves: one per compound. Click a point to select it; leave it out or put
   it back, retype its amount, or "Show in chromatogram" and click the right
   blob to rematch it ("No peak here" removes it). Everything updates at once.
4. Apply to samples: the main window's other runs, cleaned and searched the
   same way, each compound found by position and spectrum and read off its
   curve, with flags; exported as CSV.

The window keeps its own project (runs, fold, layers, blobs), apart from the
main window's, and never autosaves over the main project's file. A session
(runs, fold, settings, compounds, hand matches, exclusions and curve choices)
is saved explicitly as JSON; the results are exported as CSV tables. The
science is in `core/calbuilder.py`.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets
from .theme import fit_to_screen

from ..core import BackgroundSettings, CalibrationFit, PeakSettings
from ..core.io import discover_runs
from ..core.phase import Anchor
from ..core.background import remove_background, storage_threshold
from ..core.batch import _hash_key
from ..core.calbuilder import (CalBuilderSettings, CalCompound, CalSession, SuggestSettings,
                               amount_from_name, blob_heights, compound_templates, curve_rows,
                               load_compounds, load_session, match_blob, normalised_blob_rows, point_rows,
                               quantify_samples, sample_rows, save_session, suggest_compounds, t2_offsets,
                               usable_blobs, write_compounds)
from ..core.peaks import find_blobs
from .project import RunEntry
from ..core.siloxane import SiloxaneSettings, flag_siloxane_blobs
from ..core.quant import write_dicts
from .blob_state import LEGACY, MS_AWARE
from ..core.calbuilder import ComponentBlobs, peak_coverage, peaks_from_components

# Share of a run's 300 largest legacy blobs that MS-aware components must sit on before
# components taken from the main window are trusted (fresh detection on the real standards: 0.51;
# the bad shared-cache entry of 7 October 2026: 0.21).
MIN_COVERAGE = 0.35
from ..core.components import detect_components
from ..core.peaktable import respectrum_components
from .dialogs import SettingsDialog
from .jobs import JobManager
from .project import Project
from .render import Renderer
from .theme import DARK, LIGHT
from .views.image2d import ChromatogramView
from .views.spectrum import SpectrumView

OK = (46, 158, 107)
BAD = (224, 160, 32)
GREY = (142, 134, 160)
ISTD_COLOUR = (60, 200, 110)
TARGET_COLOUR = (126, 232, 250)
SELECTED_COLOUR = (251, 146, 60)
MISSING_COLOUR = (150, 150, 150)
CONTOUR_COLOUR = (255, 60, 255)      # magenta: absent from the turbo colour scale
RUN_FILES = "Runs (*.cdf *.CDF *.lsc *.LSC);;All files (*)"
FIXED_COLUMNS = ["Compound", "t₁ (s)", "t₂ (s)", "IS", "Quant ion", "Divide by"]
CLICK_MODES = [("pick", "put this peak's t₁, t₂ into the selected compound (a new row if none is selected)"),
               ("assign", "match the selected compound to this peak in this run only"),
               ("new", "add a new compound here"),
               ("new_is", "add an internal standard here"),
               ("none", "nothing")]


def _fmt(v, digits=4) -> str:
    return "" if v is None or not np.isfinite(v) else f"{v:.{digits}g}"


def all_outlines(blobs) -> tuple[np.ndarray, np.ndarray]:
    """Every blob's outline at once, in `BlobSet.outline_times`' (start, end, NaN) format.

    The same pixel edges as calling `outline_times` blob by blob, found in one pass
    over the label image (an edge shared by two blobs is drawn once). Blob by blob
    took about 10 s for the 5,000 blobs of a real run; this takes a fraction of one.
    MS-aware components have no footprints: their fitted contours are drawn instead."""
    if isinstance(blobs, ComponentBlobs):
        parts = [blobs.outline_times(i) for i in range(len(blobs))]
        parts = [p for p in parts if len(p[0])]
        if not parts:
            return np.array([]), np.array([])
        return np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts])
    pad = np.pad(blobs.labels, 1)
    rr, cc = np.nonzero(pad[1:, 1:-1] != pad[:-1, 1:-1])          # edges along t1, between rows
    y = (rr - 0.5) * blobs.pitch
    h_x = np.column_stack([cc, cc + 1, np.full(cc.size, np.nan)])
    h_y = np.column_stack([y, y, np.full(y.size, np.nan)])
    rr, cc = np.nonzero(pad[1:-1, 1:] != pad[1:-1, :-1])          # edges along t2, between columns
    v_x = np.column_stack([cc, cc, np.full(cc.size, np.nan)]).astype(float)
    v_y = np.column_stack([(rr - 0.5) * blobs.pitch, (rr + 0.5) * blobs.pitch, np.full(rr.size, np.nan)])
    x = np.concatenate([h_x.ravel(), v_x.ravel()])
    return blobs.t1_start + x * blobs.period, np.concatenate([h_y.ravel(), v_y.ravel()])


class _CalProject(Project):
    """The standards' own project. Never autosaved: its runs may share a folder
    with the main project, whose `project.pygcgc.json` must not be overwritten."""

    def request_autosave(self) -> None:
        pass


# ======================================================================== curves
class CurvePanel(QtWidgets.QWidget):
    """One compound's curve, its residuals, and a table of its points."""
    toggleRequested = QtCore.Signal(str, str)              # compound, run: leave out or put back
    amountEdited = QtCore.Signal(str, str, object)         # compound, run, float or None
    showRequested = QtCore.Signal(str, str)                # compound, run
    matchRequested = QtCore.Signal(str, str, str)          # compound, run, "none" or "auto"
    overrideRequested = QtCore.Signal(str, dict)           # compound, overrides ({} clears)
    amountRatioToggled = QtCore.Signal(bool)

    COLUMNS = ["Use", "Run", "Amount", "Signal", "Max point", "IS", "IS signal", "Ratio", "Back-calc.",
               "Error (%)", "Match"]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.session: CalSession | None = None
        self.points: dict = {}
        self.fits: dict = {}
        self.normalised = False
        self.selected: str | None = None     # run
        self._updating = False

        self.compounds = QtWidgets.QListWidget()
        self.compounds.setMaximumWidth(240)
        self.compounds.currentRowChanged.connect(lambda _r: self._show())
        self.weights = QtWidgets.QComboBox()
        for label, v in (("Auto", "auto"), ("None", "1"), ("1/x", "1/x"), ("1/x²", "1/x2")):
            self.weights.addItem(label, v)
        self.model = QtWidgets.QComboBox()
        for label, v in (("Auto (Mandel)", "auto"), ("Linear", "linear"), ("Quadratic", "quadratic")):
            self.model.addItem(label, v)
        self.zero = QtWidgets.QCheckBox("Through zero")
        self.log_axes = QtWidgets.QCheckBox("Log axes")
        self.ratio_x = QtWidgets.QCheckBox("x = amount / IS amount")
        self.ratio_x.setToolTip("Use when the internal standard was not added at the same amount to every run")
        for w in (self.weights, self.model):
            w.currentIndexChanged.connect(lambda _i: self._override())
        self.zero.toggled.connect(lambda _b: self._override())
        self.log_axes.toggled.connect(self._log)
        self.ratio_x.toggled.connect(lambda on: None if self._updating else self.amountRatioToggled.emit(on))

        self.plot = pg.PlotWidget()
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.plot.setLabel("left", "response / IS")
        self.curve = self.plot.plot(pen=pg.mkPen((153, 15, 61), width=1.8))
        self.used = pg.ScatterPlotItem(size=10, hoverable=True)
        self.left_out = pg.ScatterPlotItem(size=10, pen=pg.mkPen(GREY, width=2), brush=None, symbol="x")
        self.ring = pg.ScatterPlotItem(size=18, pen=pg.mkPen(SELECTED_COLOUR, width=2), brush=None)
        for item in (self.used, self.left_out, self.ring):
            self.plot.addItem(item)
        self.used.sigClicked.connect(self._point_clicked)
        self.left_out.sigClicked.connect(self._point_clicked)
        self.resid = pg.PlotWidget()
        self.resid.setMenuEnabled(False)
        self.resid.hideButtons()
        self.resid.setMaximumHeight(130)
        self.resid.setLabel("left", "error (%)")
        self.resid_points = pg.ScatterPlotItem(size=8)
        self.resid.addItem(self.resid_points)
        self.limits = [pg.InfiniteLine(angle=0, pen=pg.mkPen(GREY, style=QtCore.Qt.PenStyle.DashLine))
                       for _ in range(2)]
        for ln in self.limits:
            self.resid.addItem(ln)
        self.resid.addItem(pg.InfiniteLine(angle=0, pos=0, pen=pg.mkPen(GREY)))
        self.resid.setXLink(self.plot)
        self.stats = QtWidgets.QLabel("")
        self.stats.setWordWrap(True)
        self.stats.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)

        self.table = QtWidgets.QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.itemChanged.connect(self._item_changed)
        self.table.itemSelectionChanged.connect(self._row_selected)

        self.toggle_btn = QtWidgets.QPushButton("Leave out")
        self.toggle_btn.setToolTip("Leave the selected point out of the curve, or put it back (Delete)")
        self.show_btn = QtWidgets.QPushButton("Show in chromatogram")
        self.show_btn.setToolTip("Open this run on the Chromatograms tab, ready to click the right blob")
        self.none_btn = QtWidgets.QPushButton("No peak here")
        self.none_btn.setToolTip("This run has no peak for the compound: it gives no point")
        self.auto_btn = QtWidgets.QPushButton("Automatic match")
        self.auto_btn.setToolTip("Forget the hand match in this run and match by position again")
        self.toggle_btn.clicked.connect(self._toggle)
        self.show_btn.clicked.connect(lambda: self._emit_selected(self.showRequested))
        self.none_btn.clicked.connect(lambda: self._emit_match("none"))
        self.auto_btn.clicked.connect(lambda: self._emit_match("auto"))
        QtGui.QShortcut(QtGui.QKeySequence(QtCore.Qt.Key.Key_Delete), self, self._toggle,
                        context=QtCore.Qt.ShortcutContext.WidgetWithChildrenShortcut)

        controls = QtWidgets.QHBoxLayout()
        for w in (QtWidgets.QLabel("Weights"), self.weights, QtWidgets.QLabel("Model"), self.model,
                  self.zero, self.ratio_x, self.log_axes):
            controls.addWidget(w)
        controls.addStretch(1)
        buttons = QtWidgets.QHBoxLayout()
        for b in (self.toggle_btn, self.show_btn, self.none_btn, self.auto_btn):
            buttons.addWidget(b)
        buttons.addStretch(1)
        plots = QtWidgets.QWidget()
        pl = QtWidgets.QVBoxLayout(plots)
        pl.setContentsMargins(0, 0, 0, 0)
        pl.addWidget(self.plot, 3)
        pl.addWidget(self.resid, 1)
        lower = QtWidgets.QWidget()
        ll = QtWidgets.QVBoxLayout(lower)
        ll.setContentsMargins(0, 0, 0, 0)
        ll.addLayout(buttons)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.table, 3)
        self.stats.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop | QtCore.Qt.AlignmentFlag.AlignLeft)
        self.stats.setMinimumWidth(320)
        row.addWidget(self.stats, 2)
        ll.addLayout(row, 1)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        split.addWidget(plots)
        split.addWidget(lower)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 2)
        right = QtWidgets.QVBoxLayout()
        right.addLayout(controls)
        right.addWidget(split, 1)
        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addWidget(self.compounds)
        lay.addLayout(right, 1)
        self._update_buttons()

    def set_background(self, colour) -> None:
        self.plot.setBackground(colour)
        self.resid.setBackground(colour)

    # ------------------------------------------------------------------ data
    def set_data(self, session: CalSession, points: dict, fits: dict, normalised: bool) -> None:
        self.session, self.points, self.fits, self.normalised = session, points, fits, normalised
        cur = self.current_compound()
        self._updating = True
        self.ratio_x.setChecked(session.settings.amount_ratio)
        self.compounds.blockSignals(True)
        self.compounds.clear()
        for c in session.analytes():
            fit = fits.get(c.name)
            mark = "✅ " if isinstance(fit, CalibrationFit) and fit.ok else ("⚠️ " if isinstance(fit, CalibrationFit) else "— ")
            item = QtWidgets.QListWidgetItem(mark + c.name)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, c.name)
            item.setToolTip(fit.summary() if isinstance(fit, CalibrationFit) else str(fit or ""))
            self.compounds.addItem(item)
        names = [c.name for c in session.analytes()]
        if names:
            self.compounds.setCurrentRow(names.index(cur) if cur in names else 0)
        self.compounds.blockSignals(False)
        self._updating = False
        self._show()

    def current_compound(self) -> str | None:
        item = self.compounds.currentItem()
        return item.data(QtCore.Qt.ItemDataRole.UserRole) if item else None

    def select(self, compound: str, run: str | None = None) -> None:
        for i in range(self.compounds.count()):
            if self.compounds.item(i).data(QtCore.Qt.ItemDataRole.UserRole) == compound:
                self.selected = run
                self.compounds.setCurrentRow(i)
                self._show()
                return

    def _show(self) -> None:
        if self._updating:
            return
        c = self.current_compound()
        s = self.session
        for item in (self.curve, self.used, self.left_out, self.ring, self.resid_points):
            item.setData([], [])
        self._fill_table(c)
        self._sync_overrides(c)
        unit = f" ({s.unit})" if (s and s.unit and not s.settings.amount_ratio) else ""
        self.plot.setLabel("bottom", ("amount / IS amount" if s and s.settings.amount_ratio else "amount") + unit)
        self.resid.setLabel("bottom", "")
        if s is None or c is None:
            self.stats.setText("Add compounds (not internal standards) to see curves")
            self._update_buttons()
            return
        if not self.normalised:
            self.stats.setText("Press ➗ Normalise to IS to divide each blob by its internal standard")
            self._update_buttons()
            return
        fit = self.fits.get(c)
        excluded = s.excluded.get(c, {})
        passes = {}
        if isinstance(fit, CalibrationFit):
            passes = dict(zip(fit.labels, fit.passes()))
        ux, uy, ulab, ubr, usym, ox, oy, olab = [], [], [], [], [], [], [], []
        for run in s.runs:
            p = self.points.get((c, run))
            if p is None or not (np.isfinite(p.x) and np.isfinite(p.ratio)):
                continue
            if run in excluded:
                ox.append(p.x)
                oy.append(p.ratio)
                olab.append(run)
            else:
                ux.append(p.x)
                uy.append(p.ratio)
                ulab.append(run)
                ubr.append(pg.mkBrush(*(OK if passes.get(run, True) else BAD)) if passes else pg.mkBrush(*GREY))
                usym.append("t" if p.saturated else "o")      # a triangle: possibly saturated
        self.used.setData(ux, uy, brush=ubr, symbol=usym, data=ulab)
        self.left_out.setData(ox, oy, data=olab)
        sel = self.points.get((c, self.selected)) if self.selected else None
        if sel is not None and np.isfinite(sel.x) and np.isfinite(sel.ratio):
            self.ring.setData([sel.x], [sel.ratio])
        std = s.istd_for(c)
        if isinstance(fit, CalibrationFit):
            top = max(ux + ox) if (ux or ox) else 1.0
            xs = np.linspace(0, top * 1.05, 200)
            self.curve.setData(xs, fit.response(xs))
            self.resid_points.setData(fit.x, fit.rel_error,
                                      brush=[pg.mkBrush(*(OK if p else BAD)) for p in fit.passes()])
            self.limits[0].setPos(fit.settings.accept_re)
            self.limits[1].setPos(-fit.settings.accept_re)
            infl = fit.influential()
            verdict = "passes" if fit.ok else "<b>fails</b> the back-calculation check"
            ex = "".join(f"<br>Left out {k}: {v}" for k, v in fit.excluded.items())
            self.stats.setText(
                f"<b>{c}</b> ÷ {std.name if std else '(no internal standard)'}: {verdict}.<br>{fit.summary()}<br>"
                f"Mandel p = {fit.mandel_p:.3g}; intercept p = {fit.intercept_p:.3g}"
                + (f"<br>Worth a second look: {', '.join(infl)}" if infl else "") + ex)
        else:
            why = fit if fit else "no points"
            if std is None:
                why = "there is no internal standard to divide by (tick IS for one in the compound table)"
            self.stats.setText(f"<b>{c}</b>: no curve yet, {why}.")
        self._update_buttons()

    def _fill_table(self, c: str | None) -> None:
        s = self.session
        self._updating = True
        self.table.setRowCount(0)
        if s is None or c is None:
            self._updating = False
            return
        excluded = s.excluded.get(c, {})
        fit = self.fits.get(c)
        editable = QtCore.Qt.ItemFlag.ItemIsSelectable | QtCore.Qt.ItemFlag.ItemIsEnabled
        self.table.setRowCount(len(s.runs))
        for r, run in enumerate(s.runs):
            p = self.points.get((c, run))
            back = err = float("nan")
            if p is not None and isinstance(fit, CalibrationFit) and np.isfinite(p.ratio):
                back = float(fit.amount(np.array([p.ratio]))[0])
                err = 100.0 * (back - p.x) / p.x if (np.isfinite(p.x) and p.x > 0) else float("nan")
            use = QtWidgets.QTableWidgetItem()
            use.setFlags(editable | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            use.setCheckState(QtCore.Qt.CheckState.Unchecked if run in excluded else QtCore.Qt.CheckState.Checked)
            if run in excluded:
                use.setToolTip(excluded[run])
            match = "none" if p is None or p.blob is None else ("by hand" if p.pinned else "auto")
            if p is not None and p.pinned and p.blob is None:
                match = "no peak (by hand)"
            cells = [use, QtWidgets.QTableWidgetItem(run),
                     QtWidgets.QTableWidgetItem(_fmt(p.amount if p else None, 6))]
            if not self.normalised:
                ratio = istd = istd_sig = back_s = err_s = ""
            else:
                ratio, istd, istd_sig = _fmt(p.ratio), p.istd, _fmt(p.istd_signal)
                back_s, err_s = _fmt(back), (f"{err:+.1f}" if np.isfinite(err) else "")
            cells += [QtWidgets.QTableWidgetItem(x) for x in
                      (_fmt(p.signal), _fmt(p.height, 3), istd, istd_sig, ratio, back_s, err_s, match)]
            if p.saturated:
                cells[4].setForeground(QtGui.QBrush(QtGui.QColor(*BAD)))
                cells[4].setToolTip("At or above the saturation setting (Matching…): the detector may be "
                                    "near its limit, so this point may sit below the curve")
            for k, item in enumerate(cells):
                if k == 2:
                    item.setFlags(editable | QtCore.Qt.ItemFlag.ItemIsEditable)
                    item.setToolTip("Type the amount in this run (blank for none)")
                elif k > 0:
                    item.setFlags(editable)
                item.setData(QtCore.Qt.ItemDataRole.UserRole, run)
                self.table.setItem(r, k, item)
            if run == self.selected:
                self.table.selectRow(r)
        self._updating = False

    def _sync_overrides(self, compound) -> None:
        ov = self.session.overrides.get(compound, {}) if (self.session and compound) else {}
        self._updating = True
        self.weights.setCurrentIndex(max(self.weights.findData(ov.get("weights", "auto")), 0))
        self.model.setCurrentIndex(max(self.model.findData(ov.get("model", "auto")), 0))
        self.zero.setChecked(bool(ov.get("force_zero", False)))
        self._updating = False

    def _override(self) -> None:
        c = self.current_compound()
        if self._updating or c is None:
            return
        ov = {}
        if self.weights.currentData() != "auto":
            ov["weights"] = self.weights.currentData()
        if self.model.currentData() != "auto":
            ov["model"] = self.model.currentData()
        if self.zero.isChecked():
            ov["force_zero"] = True
        self.overrideRequested.emit(c, ov)

    def _log(self, on: bool) -> None:
        self.plot.setLogMode(x=on, y=on)
        self.resid.setLogMode(x=on, y=False)

    # ---------------------------------------------------------------- points
    def _point_clicked(self, _item, points, *_):
        if points is None or len(points) == 0:
            return
        self.selected = points[0].data()
        self._show()

    def _row_selected(self) -> None:
        if self._updating:
            return
        rows = self.table.selectionModel().selectedRows()
        if not rows:
            return
        item = self.table.item(rows[0].row(), 1)
        run = item.data(QtCore.Qt.ItemDataRole.UserRole) if item else None
        if run and run != self.selected:
            self.selected = run
            c = self.current_compound()
            p = self.points.get((c, run))
            if p is not None and np.isfinite(p.x) and np.isfinite(p.ratio):
                self.ring.setData([p.x], [p.ratio])
            else:
                self.ring.setData([], [])
            self._update_buttons()

    def _item_changed(self, item) -> None:
        if self._updating:
            return
        c = self.current_compound()
        run = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if c is None or not run:
            return
        # Signals go out after this handler returns: acting on them rebuilds the table.
        if item.column() == 0:
            self.selected = run
            QtCore.QTimer.singleShot(0, self, lambda: self.toggleRequested.emit(c, run))
        elif item.column() == 2:
            text = item.text().strip()
            try:
                value = float(text) if text else None
            except ValueError:
                QtCore.QTimer.singleShot(0, self, self._show)    # put the old value back
                return
            self.selected = run
            QtCore.QTimer.singleShot(0, self, lambda: self.amountEdited.emit(c, run, value))

    def _update_buttons(self) -> None:
        c, run = self.current_compound(), self.selected
        has = bool(c and run and self.session and run in self.session.runs)
        for b in (self.toggle_btn, self.show_btn, self.none_btn, self.auto_btn):
            b.setEnabled(has)
        left_out = has and run in self.session.excluded.get(c, {})
        self.toggle_btn.setText("Put back" if left_out else "Leave out")
        self.auto_btn.setEnabled(has and (c, run) in self.session.pins)

    def _toggle(self) -> None:
        c, run = self.current_compound(), self.selected
        if c and run:
            self.toggleRequested.emit(c, run)

    def _emit_selected(self, signal) -> None:
        c, run = self.current_compound(), self.selected
        if c and run:
            signal.emit(c, run)

    def _emit_match(self, how: str) -> None:
        c, run = self.current_compound(), self.selected
        if c and run:
            self.matchRequested.emit(c, run, how)


# ================================================================ pick runs dialog
class _PickRunsDialog(QtWidgets.QDialog):
    """Tick which of the main window's runs to use as standards."""

    def __init__(self, entries, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Runs from the main window")
        self.list = QtWidgets.QListWidget()
        for e in entries:
            item = QtWidgets.QListWidgetItem(f"{e.name}   {e.kind} {e.level}".rstrip())
            item.setData(QtCore.Qt.ItemDataRole.UserRole, str(e.path))
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.CheckState.Checked if "cal" in e.kind.lower()
                               else QtCore.Qt.CheckState.Unchecked)
            self.list.addItem(item)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                             | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(QtWidgets.QLabel("Runs typed as calibration are ticked. They are read again from disk."))
        lay.addWidget(self.list)
        lay.addWidget(buttons)

    def paths(self) -> list[str]:
        return [self.list.item(i).data(QtCore.Qt.ItemDataRole.UserRole) for i in range(self.list.count())
                if self.list.item(i).checkState() == QtCore.Qt.CheckState.Checked]


# ======================================================================== window
class SuggestDialog(QtWidgets.QDialog):
    """🔎 Find compounds from the standards alone: type each run's level amount, Find,
    tick what to keep. Calibrants grow with the level; constant peaks are what an
    internal standard looks like (so is a steady contaminant: check before keeping)."""
    COLUMNS = ["Use", "Kind", "t₁ (min)", "t₂ (s)", "Base peak", "Slope", "R²", "CV (%)", "Runs", "Worst cosine"]

    def __init__(self, runs: list[str], amounts: dict, blobs_by_run: dict,
                 settings: SuggestSettings = SuggestSettings(), parent=None, usable: dict | None = None):
        super().__init__(parent)
        self.setWindowTitle("🔎 Suggest compounds from the standards")
        fit_to_screen(self, 900, 640)
        self.runs, self.blobs_by_run, self.settings, self.usable = runs, blobs_by_run, settings, usable
        self.found = []
        self.levels = QtWidgets.QTableWidget(len(runs), 2)
        self.levels.setHorizontalHeaderLabels(["Run", "Level amount"])
        self.levels.verticalHeader().setVisible(False)
        self.levels.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        for r, run in enumerate(runs):
            name = QtWidgets.QTableWidgetItem(run)
            name.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
            self.levels.setItem(r, 0, name)
            v = amounts.get(run)
            self.levels.setItem(r, 1, QtWidgets.QTableWidgetItem("" if v is None else f"{v:g}"))
        self.levels.setMaximumHeight(40 + 30 * min(len(runs), 8))
        find = QtWidgets.QPushButton("🔎 Find")
        find.clicked.connect(self.find)
        options = QtWidgets.QPushButton("Settings…")
        options.clicked.connect(self._edit_settings)
        self.table = QtWidgets.QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSortingEnabled(False)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.summary = QtWidgets.QLabel("Type the amount of each level (read from the run names where they say), "
                                        "then Find. Siloxane-like blobs are left out.")
        self.summary.setWordWrap(True)
        buttons = QtWidgets.QDialogButtonBox()
        self.add_btn = buttons.addButton("Add ticked compounds", QtWidgets.QDialogButtonBox.ButtonRole.AcceptRole)
        buttons.addButton(QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        self.add_btn.setEnabled(False)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(find)
        row.addWidget(options)
        row.addStretch(1)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(self.levels)
        lay.addLayout(row)
        lay.addWidget(self.summary)
        lay.addWidget(self.table, 1)
        lay.addWidget(buttons)

    def amounts(self) -> dict:
        out = {}
        for r, run in enumerate(self.runs):
            item = self.levels.item(r, 1)
            try:
                v = float(item.text()) if item and item.text().strip() else None
            except ValueError:
                v = None
            if v is not None and v > 0:
                out[run] = v
        return out

    def _edit_settings(self) -> None:
        dlg = SettingsDialog("Suggestion settings", self.settings, parent=self)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.settings = dlg.settings()

    def find(self) -> None:
        amounts = self.amounts()
        if len(amounts) < 3:
            self.summary.setText("Give an amount for at least three runs")
            return
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            self.found = suggest_compounds(self.blobs_by_run, amounts, self.settings, self.usable)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self.table.setRowCount(len(self.found))
        for k, g in enumerate(self.found):
            use = QtWidgets.QTableWidgetItem()
            use.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            use.setCheckState(QtCore.Qt.CheckState.Checked if (g.kind == "calibrant" or g.name)
                              else QtCore.Qt.CheckState.Unchecked)
            kind = "calibrant" if g.kind == "calibrant" else (f"IS: {g.name}" if g.name else "constant (IS?)")
            cells = [use, kind, f"{g.t1 / 60:.2f}",
                     f"{g.t2:.3f}", str(g.base_peak), f"{g.slope:.2f}", f"{g.r2:.3f}", f"{100 * g.cv:.0f}",
                     f"{len(g.blobs)}/{len(amounts)}", f"{g.min_cosine:.2f}"]
            for c, v in enumerate(cells):
                item = v if isinstance(v, QtWidgets.QTableWidgetItem) else QtWidgets.QTableWidgetItem(v)
                if c:
                    item.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled | QtCore.Qt.ItemFlag.ItemIsSelectable)
                self.table.setItem(k, c, item)
        n_cal = sum(g.kind == "calibrant" for g in self.found)
        n_is = sum(bool(g.name) for g in self.found)
        self.summary.setText(
            f"{n_cal} calibrants and {n_is} recognised internal standards (ticked), and "
            f"{len(self.found) - n_cal - n_is} other constant peaks. Tick any other internal standards: "
            "constant peaks are added as IS. Calibrant names are placeholders; rename them in the "
            "compound table.")
        self.add_btn.setEnabled(bool(self.found))

    def chosen(self) -> list:
        return [g for k, g in enumerate(self.found)
                if self.table.item(k, 0).checkState() == QtCore.Qt.CheckState.Checked]


def _banner(text: str) -> QtWidgets.QLabel:
    """The "what to do here" line at the top of each tab."""
    label = QtWidgets.QLabel(text)
    label.setWordWrap(True)
    label.setTextFormat(QtCore.Qt.TextFormat.RichText)
    label.setStyleSheet("QLabel { padding: 6px; border-radius: 4px; background: rgba(153, 15, 61, 28); }")
    return label


class ApplyPanel(QtWidgets.QWidget):
    """4 Apply to samples: tick runs, quantify, read and export the results."""
    quantifyRequested = QtCore.Signal(list)          # run names
    exportRequested = QtCore.Signal()
    showRequested = QtCore.Signal(str, float, float)  # run, t1, t2 (in the main window)

    MODES = [("Amount", "amount"), ("Concentration (per m³)", "concentration"), ("Ratio to IS", "ratio")]

    def __init__(self, parent=None):
        super().__init__(parent)
        self.results = []
        self.volumes: dict = {}
        self.unit = ""
        self.banner = _banner("")
        self.samples = QtWidgets.QListWidget()
        self.samples.setMaximumWidth(320)
        self.samples.setToolTip("The main window's runs, other than the standards. Tick the ones to quantify.")
        all_btn = QtWidgets.QPushButton("Tick all")
        none_btn = QtWidgets.QPushButton("Tick none")
        all_btn.clicked.connect(lambda: self._tick_all(True))
        none_btn.clicked.connect(lambda: self._tick_all(False))
        self.go = QtWidgets.QPushButton("▶ Quantify ticked runs")
        self.go.setStyleSheet("QPushButton { font-weight: bold; padding: 6px; }")
        self.go.clicked.connect(lambda: self.quantifyRequested.emit(self.ticked()))
        self.mode = QtWidgets.QComboBox()
        for label, key in self.MODES:
            self.mode.addItem(label, key)
        self.mode.currentIndexChanged.connect(lambda _i: self._fill())
        self.table = QtWidgets.QTableWidget(0, 0)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.cellDoubleClicked.connect(self._show)
        self.table.setToolTip("Amber: flagged (hover for why). n.d.: not found. Double-click to see it in the "
                              "main window.")
        self.export_btn = QtWidgets.QPushButton("📤 Export results…")
        self.export_btn.clicked.connect(self.exportRequested.emit)
        self.export_btn.setEnabled(False)
        self.summary = QtWidgets.QLabel("")
        self.summary.setWordWrap(True)

        left = QtWidgets.QVBoxLayout()
        left.addWidget(QtWidgets.QLabel("Sample runs (from the main window)"))
        left.addWidget(self.samples, 1)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(all_btn)
        row.addWidget(none_btn)
        left.addLayout(row)
        left.addWidget(self.go)
        right = QtWidgets.QVBoxLayout()
        bar = QtWidgets.QHBoxLayout()
        bar.addWidget(QtWidgets.QLabel("Show"))
        bar.addWidget(self.mode)
        bar.addStretch(1)
        bar.addWidget(self.export_btn)
        right.addLayout(bar)
        right.addWidget(self.table, 1)
        right.addWidget(self.summary)
        body = QtWidgets.QHBoxLayout()
        body.addLayout(left)
        body.addLayout(right, 1)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addWidget(self.banner)
        lay.addLayout(body, 1)

    def set_samples(self, entries, standards: set) -> None:
        """The main window's runs; ticked by default unless they are standards or calibration runs."""
        before = {self.samples.item(i).data(QtCore.Qt.ItemDataRole.UserRole):
                  self.samples.item(i).checkState() for i in range(self.samples.count())}
        self.samples.clear()
        for e in entries:
            if str(e.path) in standards:
                continue
            item = QtWidgets.QListWidgetItem(f"{e.name}   {e.kind}".rstrip())
            item.setData(QtCore.Qt.ItemDataRole.UserRole, e.name)
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            default = (QtCore.Qt.CheckState.Unchecked if ("cal" in e.kind.lower() or e.run is None)
                       else QtCore.Qt.CheckState.Checked)
            item.setCheckState(before.get(e.name, default))
            if e.run is None:
                item.setToolTip("Not loaded in the main window")
            self.samples.addItem(item)

    def ticked(self) -> list[str]:
        return [self.samples.item(i).data(QtCore.Qt.ItemDataRole.UserRole) for i in range(self.samples.count())
                if self.samples.item(i).checkState() == QtCore.Qt.CheckState.Checked]

    def _tick_all(self, on: bool) -> None:
        for i in range(self.samples.count()):
            self.samples.item(i).setCheckState(QtCore.Qt.CheckState.Checked if on else QtCore.Qt.CheckState.Unchecked)

    def set_results(self, results, volumes: dict, unit: str) -> None:
        self.results, self.volumes, self.unit = results, volumes, unit
        self.export_btn.setEnabled(bool(results))
        self._fill()

    def _fill(self) -> None:
        runs = list(dict.fromkeys(r.run for r in self.results))
        compounds = list(dict.fromkeys(r.compound for r in self.results))
        by = {(r.run, r.compound): r for r in self.results}
        mode = self.mode.currentData()
        self.table.clear()
        self.table.setRowCount(len(runs))
        self.table.setColumnCount(len(compounds))
        self.table.setHorizontalHeaderLabels(compounds)
        self.table.setVerticalHeaderLabels(runs)
        flagged = missing = 0
        for i, run in enumerate(runs):
            vol = self.volumes.get(run)
            for j, c in enumerate(compounds):
                r = by.get((run, c))
                if r is None:
                    continue
                if not r.found:
                    text = "n.d."
                    missing += 1
                elif mode == "ratio":
                    text = _fmt(r.ratio)
                elif mode == "concentration":
                    text = _fmt(r.amount / vol) if (vol and np.isfinite(r.amount)) else "no volume" if r.found else ""
                else:
                    text = _fmt(r.amount)
                item = QtWidgets.QTableWidgetItem(text)
                item.setData(QtCore.Qt.ItemDataRole.UserRole, (run, r.t1, r.t2))
                tip = [f"{c} in {run}", f"amount {_fmt(r.amount)} ± {_fmt(r.se, 2)} {self.unit}".rstrip(),
                       f"ratio to {r.istd or 'IS'} {_fmt(r.ratio)}", f"spectral match {_fmt(r.cosine, 2)}"]
                if r.flags:
                    tip.append("⚠️ " + "; ".join(r.flags))
                    if r.found:
                        item.setBackground(QtGui.QBrush(QtGui.QColor(*BAD, 90)))
                        flagged += 1
                if not r.found:
                    item.setForeground(QtGui.QBrush(QtGui.QColor(*MISSING_COLOUR)))
                item.setToolTip("\n".join(tip))
                self.table.setItem(i, j, item)
        self.table.resizeColumnsToContents()
        unit = {"amount": self.unit or "standards' unit", "concentration": f"{self.unit or 'amount'} per m³",
                "ratio": "response / IS"}[mode]
        self.summary.setText(
            f"{len(runs)} runs × {len(compounds)} compounds, in {unit}. {flagged} results flagged (amber, hover "
            f"for why), {missing} not found. Double-click a cell to see the peak in the main window."
            if runs else "")

    def _show(self, row: int, col: int) -> None:
        item = self.table.item(row, col)
        if item is not None:
            run, t1, t2 = item.data(QtCore.Qt.ItemDataRole.UserRole)
            if np.isfinite(t1):
                self.showRequested.emit(run, t1, t2)


class CalibrationWindow(QtWidgets.QMainWindow):
    """📈 Build calibration curves from a separate set of standard runs."""
    changed = QtCore.Signal()

    def __init__(self, main_project: Project | None = None, dark: bool = False, parent=None):
        super().__init__(parent)
        self.setWindowTitle("📈 Calibration builder")
        self.setWindowFlag(QtCore.Qt.WindowType.Window, True)
        fit_to_screen(self, 1400, 900)
        self.main_project = main_project
        self.jobs = JobManager(self)
        self.project = _CalProject(self.jobs, self)
        if main_project is not None:
            # The main window's period is only a starting point: standards may have been run with
            # another method (the August 2026 standards at 10 s, the exhaust runs at 8 s), so the
            # period is estimated from the standards themselves unless a method.json sits beside them.
            self.project.period = main_project.period
            self.project.period_confirmed = False
            self.project.bg_settings = main_project.bg_settings
            self.project.blobs.settings = main_project.blobs.settings
            self.project.blobs.component_settings = main_project.blobs.component_settings
            self.project.blobs.candidate_settings = main_project.blobs.candidate_settings
        # Peaks the same way as the main window's 🔍 step: MS-aware components unless it is set to
        # legacy blobs (Settings menu here). Same settings and fold, so the cache is shared.
        self.peak_method = main_project.blobs.method if main_project is not None else MS_AWARE
        self.project.blobs.method = self.peak_method
        self.ms_peaks: dict = {}          # RunEntry -> ComponentBlobs, when peaks are MS-aware
        self.processing = False           # the standards are being cleaned and searched
        self.renderer = Renderer(self.project, self.jobs, parent=self)
        self.session = CalSession()
        self.suggest_settings = SuggestSettings()
        self.heights: dict = {}           # run name -> largest raw point of each blob
        self.normalised = False
        self.resolved: dict = {}
        self.points: dict = {}
        self.fits: dict = {}
        self._undo: list[dict] = []
        self._redo: list[dict] = []
        self._updating = False
        self._outline_cache = None
        self.session_path: Path | None = None
        # Loading runs cleans them and finds their blobs straight away (tests switch this off to
        # set phases by hand first).
        self.auto_process = True
        self.main_window = parent if hasattr(parent, "runs_panel") else None
        self.sample_blobs: dict = {}      # main-window run name -> BlobSet, found with this window's settings
        self.sample_offsets: dict = {}    # ... -> t2 offset from the standards
        self.sample_results: list = []

        # ---------------------------------------------------------- widgets
        self.view = ChromatogramView()
        self.view.set_empty_message("Choose the standard runs first (tab 1)")
        self.view.set_blob_overlay_visible(True)
        self.renderer.ready.connect(self.view.show_result)
        self.renderer.ready.connect(self._view_levels)
        self.renderer.ready.connect(lambda _r: self._draw_overlay())
        self.scale = QtWidgets.QComboBox()
        for label, key in (("Log", "log"), ("Square root", "sqrt"), ("Linear", "linear")):
            self.scale.addItem(label, key)
        self.scale.setToolTip("Colour scale of the chromatogram")
        self.scale.currentIndexChanged.connect(lambda _i: self._set_scale())
        self.view.pointClicked.connect(self._view_clicked)
        self.view.hovered.connect(lambda t1, t2, v: self.statusBar().showMessage(
            f"t₁ {t1 / 60:.2f} min   t₂ {t2:.3f} s   {v:.4g}", 3000))
        # The spectrum of the peak clicked, or of the selected compound's peak in the run shown.
        self.spectrum = SpectrumView()
        self.spectrum.header.setText("Click a peak, or select a compound, to see its spectrum")
        self.spectrum.ionClicked.connect(self._spectrum_ion)

        self.runs = QtWidgets.QListWidget()
        self.runs.setMaximumWidth(260)
        self.runs.currentRowChanged.connect(self._run_selected)
        self.ion = QtWidgets.QLineEdit()
        self.ion.setPlaceholderText("TIC")
        self.ion.setToolTip("Ion or ion expression to show (blank for the TIC), e.g. 57 or 57+71")
        self.ion.setMaximumWidth(120)
        self.ion.editingFinished.connect(self._set_ion)
        self.click_mode = QtWidgets.QComboBox()
        for key, label in CLICK_MODES:
            self.click_mode.addItem(label, key)
        self.click_mode.setToolTip("What a click on the chromatogram does. By default the clicked peak's t₁ and t₂ "
                                   "go into the selected compound, then the next compound without a position is "
                                   "selected; with none selected, each click adds a row.")
        self.all_blobs = QtWidgets.QCheckBox("All blob outlines")
        self.all_blobs.setChecked(True)
        self.all_blobs.setToolTip("Outline every blob, not only the matched ones (bold)")
        self.all_blobs.toggled.connect(lambda _on: self._draw_overlay())
        self.none_here = QtWidgets.QPushButton("No peak in this run")
        self.none_here.setToolTip("The selected compound has no peak in the run shown")
        self.auto_here = QtWidgets.QPushButton("Automatic match")
        self.auto_here.setToolTip("Forget the hand match of the selected compound in the run shown")
        self.none_here.clicked.connect(lambda: self._match_current("none"))
        self.auto_here.clicked.connect(lambda: self._match_current("auto"))

        self.table = QtWidgets.QTableWidget(0, len(FIXED_COLUMNS))
        self.table.setHorizontalHeaderLabels(FIXED_COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.itemChanged.connect(self._table_edited)
        self.table.itemSelectionChanged.connect(self._compound_selected)
        self.add_row = QtWidgets.QPushButton("➕ Compound")
        self.add_row.setToolTip("Add an empty row (or click a peak with 'add a compound here')")
        self.del_row = QtWidgets.QPushButton("🗑 Remove")
        self.fill_btn = QtWidgets.QPushButton("Same amount in every run…")
        self.add_row.clicked.connect(self._add_empty_row)
        self.del_row.clicked.connect(self._remove_selected)
        self.fill_btn.clicked.connect(self._fill_amount)
        self.unit = QtWidgets.QLineEdit()
        self.unit.setPlaceholderText("unit, e.g. ng")
        self.unit.setMaximumWidth(110)
        self.unit.editingFinished.connect(self._set_unit)

        self.curves = CurvePanel()
        self.curves.toggleRequested.connect(self._toggle_point)
        self.curves.amountEdited.connect(self._set_amount)
        self.curves.showRequested.connect(self.show_in_chromatogram)
        self.curves.matchRequested.connect(self._set_match)
        self.curves.overrideRequested.connect(self._set_override)
        self.curves.amountRatioToggled.connect(self._set_amount_ratio)

        self.guide = QtWidgets.QLabel()
        self.guide.setWordWrap(True)
        self.guide.setContentsMargins(6, 2, 6, 2)

        # 1 Standards
        self.standards_banner = _banner("")
        self.run_table = QtWidgets.QTableWidget(0, 4)
        self.run_table.setHorizontalHeaderLabels(["Standard run", "Level amount", "Status", "Phase (s)"])
        self.run_table.verticalHeader().setVisible(False)
        self.run_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.run_table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.run_table.itemChanged.connect(self._level_edited)
        self.from_main_btn = QtWidgets.QPushButton("Use the runs selected in the main window")
        self.from_main_btn.setStyleSheet("QPushButton { font-weight: bold; padding: 6px; }")
        self.from_main_btn.clicked.connect(self.runs_from_main)
        self.open_folder_btn = QtWidgets.QPushButton("📂 Open a folder…")
        self.open_folder_btn.clicked.connect(self.open_folder_dialog)
        self.open_runs_btn = QtWidgets.QPushButton("Open files…")
        self.open_runs_btn.clicked.connect(self.open_runs_dialog)
        remove_btn = QtWidgets.QPushButton("Remove selected")
        remove_btn.clicked.connect(self._remove_selected_runs)
        std_bar = QtWidgets.QHBoxLayout()
        for w in (self.from_main_btn, self.open_folder_btn, self.open_runs_btn, remove_btn):
            std_bar.addWidget(w)
        std_bar.addStretch(1)
        std_bar.addWidget(QtWidgets.QLabel("Amounts in"))
        std_bar.addWidget(self.unit)
        standards = QtWidgets.QWidget()
        sl = QtWidgets.QVBoxLayout(standards)
        sl.setContentsMargins(4, 4, 4, 4)
        sl.addWidget(self.standards_banner)
        sl.addLayout(std_bar)
        sl.addWidget(self.run_table, 1)
        self._std_layout = sl

        # 2 Compounds: how to add them, then the chromatogram and the table
        self.compounds_banner = _banner("")
        self.suggest_btn = QtWidgets.QPushButton("🔎 Suggest from the standards…")
        self.suggest_btn.setStyleSheet("QPushButton { font-weight: bold; padding: 6px; }")
        self.suggest_btn.setToolTip("No table needed: compounds that grow with the level, and the internal "
                                    "standards, found in the runs")
        self.suggest_btn.clicked.connect(self.suggest_dialog)
        load_btn = QtWidgets.QPushButton("📄 Load a compound table…")
        load_btn.setToolTip("A CSV of compound, t1_s, t2_s, role (istd), quant_ion, istd and amounts; "
                            "a target_table.csv also works")
        load_btn.clicked.connect(self.load_compounds_dialog)
        save_btn = QtWidgets.QPushButton("💾 Save the table…")
        save_btn.clicked.connect(self.save_compounds_dialog)
        add_bar = QtWidgets.QHBoxLayout()
        for w in (self.suggest_btn, load_btn, save_btn):
            add_bar.addWidget(w)
        add_bar.addStretch(1)

        # ---------------------------------------------------------- layout
        view_bar = QtWidgets.QHBoxLayout()
        view_bar.addWidget(QtWidgets.QLabel("Ion"))
        view_bar.addWidget(self.ion)
        view_bar.addWidget(QtWidgets.QLabel("Scale"))
        view_bar.addWidget(self.scale)
        view_bar.addWidget(self.all_blobs)
        view_bar.addSpacing(12)
        view_bar.addWidget(QtWidgets.QLabel("A click on the chromatogram:"))
        view_bar.addWidget(self.click_mode)
        view_bar.addWidget(self.none_here)
        view_bar.addWidget(self.auto_here)
        view_bar.addStretch(1)
        table_bar = QtWidgets.QHBoxLayout()
        for w in (self.add_row, self.del_row, self.fill_btn):
            table_bar.addWidget(w)
        table_bar.addStretch(1)
        upper = QtWidgets.QWidget()
        ul = QtWidgets.QVBoxLayout(upper)
        ul.setContentsMargins(0, 0, 0, 0)
        ul.addLayout(view_bar)
        beside = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        beside.addWidget(self.view)
        beside.addWidget(self.spectrum)
        beside.setStretchFactor(0, 3)
        beside.setStretchFactor(1, 1)
        beside.setSizes([900, 340])
        ul.addWidget(beside, 1)
        lower = QtWidgets.QWidget()
        ll = QtWidgets.QVBoxLayout(lower)
        ll.setContentsMargins(0, 0, 0, 0)
        ll.addLayout(table_bar)
        ll.addWidget(self.table, 1)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        split.addWidget(upper)
        split.addWidget(lower)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 2)
        chrom = QtWidgets.QWidget()
        outer = QtWidgets.QVBoxLayout(chrom)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.addWidget(self.compounds_banner)
        outer.addLayout(add_bar)
        cl = QtWidgets.QHBoxLayout()
        left = QtWidgets.QVBoxLayout()
        left.addWidget(QtWidgets.QLabel("Run shown"))
        left.addWidget(self.runs, 1)
        cl.addLayout(left)
        cl.addWidget(split, 1)
        outer.addLayout(cl, 1)
        self.curves_banner = _banner("")
        curves_page = QtWidgets.QWidget()
        cpl = QtWidgets.QVBoxLayout(curves_page)
        cpl.setContentsMargins(4, 4, 4, 0)
        cpl.addWidget(self.curves_banner)
        cpl.addWidget(self.curves, 1)
        self.apply = ApplyPanel()
        self.apply.quantifyRequested.connect(self.quantify_samples)
        self.apply.exportRequested.connect(self.export_results_dialog)
        self.apply.showRequested.connect(self._show_in_main)
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(standards, "1  Standards")
        self.tabs.addTab(chrom, "2  Compounds")
        self.tabs.addTab(curves_page, "3  Curves")
        self.tabs.addTab(self.apply, "4  Apply to samples")
        self.tabs.currentChanged.connect(self._tab_changed)
        central = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(central)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.guide)
        lay.addWidget(self.tabs, 1)
        self.setCentralWidget(central)
        self._build_toolbar()
        self.busy = QtWidgets.QLabel("")
        self.statusBar().addPermanentWidget(self.busy)

        self.project.runsChanged.connect(self._runs_changed)
        self.project.message.connect(lambda m: self.statusBar().showMessage(m, 8000))
        self.project.blobs.changed.connect(self.recompute)
        self.project.blobs.runningChanged.connect(lambda _b: self._update_guide())
        self.project.foldChanged.connect(self._fold_changed)
        self.jobs.jobChanged.connect(self._job_changed)
        self.jobs.busyChanged.connect(lambda n: self.busy.setText("" if n == 0 else f"⏳ {n} job(s)"))
        QtGui.QShortcut(QtGui.QKeySequence.StandardKey.Undo, self, self.undo)
        QtGui.QShortcut(QtGui.QKeySequence("Ctrl+Shift+Z"), self, self.redo)
        self.set_dark(dark)
        self.recompute()

    def _build_toolbar(self) -> None:
        """A short toolbar: the session, the settings, undo. The work itself is on the tabs."""
        tb = self.addToolBar("Calibration")
        tb.setMovable(False)
        tb.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextOnly)

        def act(menu, text, slot, tip=""):
            a = QtGui.QAction(text, self)
            a.setToolTip(tip)
            a.setStatusTip(tip)
            a.triggered.connect(lambda _checked=False: slot())
            menu.addAction(a)
            return a

        def drop(text, tip):
            button = QtWidgets.QToolButton()
            button.setText(text)
            button.setToolTip(tip)
            button.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
            menu = QtWidgets.QMenu(button)
            menu.setToolTipsVisible(True)
            button.setMenu(menu)
            tb.addWidget(button)
            return menu

        session = drop("Session ▾", "Save or open a session, export the tables")
        self.act_save = act(session, "Save session…", self.save_session_dialog, "Everything needed to carry on later")
        self.act_load = act(session, "Open session…", self.open_session_dialog, "Carry on with a saved session")
        session.addSeparator()
        self.act_export = act(session, "📤 Export curves and points…", self.export_dialog,
                              "The curves, points, compounds and every normalised blob as CSV tables")
        self.act_export_results = act(session, "📤 Export sample results…", self.export_results_dialog,
                                      "What tab 4 shows, as CSV")
        self.act_save_compounds = act(session, "Save the compound table…", self.save_compounds_dialog)
        self.act_compounds = act(session, "Load a compound table…", self.load_compounds_dialog)

        settings = drop("⚙️ Settings ▾", "Background, blobs, matching and curve rules")
        self.act_bg_settings = act(settings, "🧹 Background removal…", self.edit_background_settings)
        self.act_peak_settings = act(settings, "🔍 Blob detection…", self.edit_peak_settings)
        self.act_match_settings = act(settings, "Matching and saturation…", self.edit_match_settings,
                                      "Search window for matching, the nearest-IS distance scale, saturation flag")
        self.act_curve_settings = act(settings, "Curve rules…", self.edit_curve_settings,
                                      "Default weights, model and acceptance limits for every curve")
        settings.addSeparator()
        self.act_line_up = act(settings, "Line up phases", lambda: None,
                               "After finding blobs, move any run whose second dimension sits apart from the "
                               "others (judged from blobs of the same spectrum) and find its blobs again")
        self.act_line_up.setCheckable(True)
        self.act_line_up.setChecked(True)
        settings.addSection("Peaks")
        self.act_ms_aware = act(settings, "MS-aware peaks (as the main window)", lambda: self.set_peak_method(MS_AWARE),
                                "Overlapping components fitted from intensity and spectra: step 4's default")
        self.act_legacy = act(settings, "Legacy blobs", lambda: self.set_peak_method(LEGACY),
                              "Hard-border TIC blobs: faster, and a click adds one where detection missed a peak")
        for a in (self.act_ms_aware, self.act_legacy):
            a.setCheckable(True)
        self.act_ms_aware.setChecked(self.peak_method == MS_AWARE)
        self.act_legacy.setChecked(self.peak_method == LEGACY)
        settings.addSeparator()
        self.act_run = act(settings, "▶ Clean and find peaks again", self.run_pipeline,
                           "Needed only after changing a setting; loading runs does it")
        self.act_normalise = act(settings, "➗ Normalise to IS", self.normalise,
                                 "Done as soon as there are blobs and an internal standard")
        tb.addSeparator()
        undo = QtGui.QAction("↶ Undo", self)
        undo.setToolTip("Undo (Ctrl+Z)")
        undo.triggered.connect(lambda _c=False: self.undo())
        redo = QtGui.QAction("↷ Redo", self)
        redo.setToolTip("Redo (Ctrl+Shift+Z)")
        redo.triggered.connect(lambda _c=False: self.redo())
        tb.addAction(undo)
        tb.addAction(redo)

        # The fold lives with the runs, out of the way.
        self.period = QtWidgets.QDoubleSpinBox()
        self.period.setRange(0.5, 30.0)
        self.period.setDecimals(4)
        self.period.setSingleStep(0.01)
        self.period.setSuffix(" s")
        self.period.setValue(self.project.period)
        self.period.setToolTip("Modulation period of the standard runs (estimated from them when loaded)")
        self.period.editingFinished.connect(self._set_period)
        phase_btn = QtWidgets.QPushButton("🔄 Estimate phases again")
        phase_btn.setToolTip("From the internal standards that have a quantifier ion, otherwise from the quiet band")
        phase_btn.clicked.connect(lambda: self.project.estimate_phases())
        process_btn = QtWidgets.QPushButton("▶ Process the standards")
        process_btn.setToolTip("Clean, find peaks and line up phases; loading runs does this by itself, so it "
                               "is needed only after changing a setting")
        process_btn.clicked.connect(lambda: self.run_pipeline())
        fold_row = QtWidgets.QHBoxLayout()
        fold_row.addWidget(process_btn)
        fold_row.addSpacing(20)
        fold_row.addWidget(QtWidgets.QLabel("Modulation period"))
        fold_row.addWidget(self.period)
        fold_row.addWidget(phase_btn)
        fold_row.addStretch(1)
        self._std_layout.addLayout(fold_row)

    def set_dark(self, dark: bool) -> None:
        c = DARK if dark else LIGHT
        self.view.set_background(c["plot"])
        self.spectrum.set_background(c["plot"])
        self.curves.set_background(c["plot"])
        fg = QtGui.QColor(c["text"])
        for gv in self.findChildren(pg.GraphicsView):
            for item in (gv.scene().items() if gv.scene() is not None else []):
                if isinstance(item, pg.AxisItem):
                    item.setPen(fg)
                    item.setTextPen(fg)

    def say(self, text: str) -> None:
        self.statusBar().showMessage(text, 8000)

    def _job_changed(self, job) -> None:
        if job.state == "running":
            self.busy.setText(f"⏳ {job.name} {job.message}")
        elif job.state == "failed":
            self.say(f"{job.name} failed: {job.message}")
        if job.state in ("done", "failed") and self.processing:
            self._fill_run_table()           # each run's status as it finishes
            self._update_guide()

    # ============================================================ runs
    def open_runs_dialog(self) -> None:
        paths, _ = QtWidgets.QFileDialog.getOpenFileNames(self, "Standard runs", "", RUN_FILES)
        if paths:
            self.add_runs(paths)

    def open_folder_dialog(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Folder of standard runs", self._folder())
        if folder:
            self.open_folder(folder)

    def open_folder(self, folder, on_loaded=None) -> None:
        """Every run in a folder, the untreated 70 eV channel of each native acquisition
        (as the main window's folder opening chooses)."""
        QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
        try:
            paths = discover_runs(folder)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        if not paths:
            self.say(f"No supported runs in {folder}")
            return
        self.add_runs(paths, on_loaded)

    def _fold_changed(self, _index: int) -> None:
        if abs(self.period.value() - self.project.period) > 1e-9:
            self.period.blockSignals(True)
            self.period.setValue(self.project.period)
            self.period.blockSignals(False)
        self._update_guide()

    def main_selection(self) -> list:
        """The main window's runs that are selected in its Runs table (two or more), else []."""
        try:
            rows = self.main_window.runs_panel.selected_rows() if self.main_window is not None else []
        except (AttributeError, RuntimeError):
            rows = []
        entries = self.main_project.entries if self.main_project is not None else []
        picked = [entries[i] for i in rows if 0 <= i < len(entries)]
        return picked if len(picked) >= 2 else []

    def runs_from_main(self) -> None:
        """The runs selected in the main window, or (with fewer than two selected) a list to tick,
        calibration runs ticked. Runs already read there are not read again."""
        entries = self.main_project.entries if self.main_project is not None else []
        if not entries:
            self.say("No runs are open in the main window: open your folder there first, or use 📂 here")
            return
        picked = self.main_selection()
        if not picked:
            dlg = _PickRunsDialog(entries, self)
            if dlg.exec() != QtWidgets.QDialog.DialogCode.Accepted or not dlg.paths():
                return
            chosen = set(dlg.paths())
            picked = [e for e in entries if str(e.path) in chosen]
        self.add_main_entries(picked)

    def add_main_entries(self, picked: list, on_loaded=None) -> None:
        """Standards taken from the main window: the same Run objects (never changed by anyone),
        at the main window's period and phases, then processed like any other."""
        known = {e.path.resolve() for e in self.project.entries}
        main = self.main_project
        loaded = [e for e in picked if e.run is not None and e.path.resolve() not in known]
        unread = [e.path for e in picked if e.run is None and e.path.resolve() not in known]
        if loaded:
            if not self.project.entries:
                self.project.period, self.project.period_confirmed = main.period, main.period_confirmed
                self.project.period_method, self.project.period_detail = main.period_method, main.period_detail
                self.period.setValue(self.project.period)
            same_period = abs(self.project.period - main.period) < 1e-6
            new = [RunEntry(path=e.path, name=e.name, kind=e.kind, level=e.level, run=e.run,
                            phase=e.phase if same_period else 0.0,
                            phase_method=e.phase_method if same_period else "unset",
                            phase_detail=e.phase_detail, timing=e.timing, import_meta=e.import_meta)
                   for e in loaded]
            self.project.entries.extend(new)
            self.project.runsChanged.emit()
            self.project.foldChanged.emit(-1)
            if same_period:
                self._loaded(on_loaded)
            else:
                idx = [self.project.entries.index(e) for e in new]
                self.project.auto_fold(apply_to=idx, on_done=lambda *_: self._loaded(on_loaded))
        if unread:
            self.add_runs(unread, on_loaded if not loaded else None)
        self.say(f"{len(loaded) + len(unread)} standard runs from the main window")

    def add_runs(self, paths, on_loaded=None) -> None:
        """Read runs from disk, estimate the period (unless a method.json says it) and each
        phase, then clean them and find their blobs."""
        paths = [Path(p) for p in paths]
        if not self.project.entries:
            for folder in (paths[0].parent, paths[0].parent.parent):
                method = folder / "method.json"
                if method.exists():
                    try:
                        self.project.period = float(json.loads(method.read_text())["modulation_period_s"])
                        self.project.period_confirmed = True
                        self.project.period_method, self.project.period_detail = "file", f"From {method.name}"
                        break
                    except (KeyError, ValueError, OSError):
                        pass
            self.period.setValue(self.project.period)
        known = {e.path.resolve() for e in self.project.entries}
        new = [RunEntry(path=p, name=p.stem) for p in paths if p.resolve() not in known]
        if not new:
            return
        self.project.entries.extend(new)
        self.project.runsChanged.emit()

        def read():
            idx = [i for i, e in enumerate(self.project.entries) if e in new and e.run is not None]
            if not idx:
                return
            self.project.auto_fold(apply_to=idx, on_done=lambda *_: self._loaded(on_loaded))

        self.project._read(new, read, estimate=False)

    def _loaded(self, on_loaded=None) -> None:
        """Runs are read and folded: process them unless told not to."""
        self._runs_changed()
        if on_loaded:
            on_loaded()
        if self.auto_process and not self.project.blobs.running:
            self.say("Cleaning the standards and finding their blobs (a minute or two a run)…")
            self.run_pipeline(on_ready=self._processed)

    def _processed(self) -> None:
        """Standards processed: on to the compounds, which say what to do next."""
        if self.tabs.currentIndex() == 0:
            self.tabs.setCurrentIndex(1)

    def _runs_changed(self) -> None:
        names = [e.name for e in self.project.entries]
        if names != self.session.runs:
            self.session.runs = names
            # A new run's level amount, from its name where it says ("L6_5ng") or the main window's level.
            for e in self.project.entries:
                if e.name not in self.session.levels:
                    v = amount_from_name(e.name)
                    if v is None:
                        v = amount_from_name(e.level or "")
                    if v is not None:
                        self.session.levels[e.name] = v
            self._outline_cache = None
            self.recompute()
        self._fill_run_table()
        cur = self.runs.currentRow()
        self.runs.blockSignals(True)
        self.runs.clear()
        for e in self.project.entries:
            b = self.peaks.get(e)
            state = ("⚠️ " + e.error) if e.error else ("reading…" if e.run is None else
                                                          f"{len(b)} blobs" if b is not None else
                                                          "🧹 cleaned" if e.background is not None else "raw")
            item = QtWidgets.QListWidgetItem(f"{e.name}\n   {state}")
            item.setToolTip(f"{e.path}\nphase {e.phase:.3f} s ({e.phase_method})")
            self.runs.addItem(item)
        self.runs.blockSignals(False)
        if self.project.entries:
            row = cur if 0 <= cur < len(self.project.entries) else 0
            self.runs.setCurrentRow(row)
            if self.renderer.index != row:
                self.renderer.set_index(row)
        self._fill_table()
        self._update_guide()

    def _fill_run_table(self) -> None:
        p = self.project
        self._updating = True
        self.run_table.setRowCount(len(p.entries))
        running = self.processing
        for r, e in enumerate(p.entries):
            b = self.peaks.get(e)
            if e.error:
                status = "⚠️ " + e.error
            elif e.run is None:
                status = "reading…"
            elif b is not None:
                status = f"✓ {len(b)} {'MS-aware peaks' if self.peak_method == MS_AWARE else 'blobs'}"
            elif running:
                status = "processing…"
            else:
                status = "waiting"
            level = self.session.levels.get(e.name)
            cells = [e.name, "" if level is None else f"{level:g}", status,
                     f"{e.phase:.3f} ({e.phase_method})"]
            for k, text in enumerate(cells):
                item = QtWidgets.QTableWidgetItem(text)
                item.setData(QtCore.Qt.ItemDataRole.UserRole, e.name)
                flags = QtCore.Qt.ItemFlag.ItemIsSelectable | QtCore.Qt.ItemFlag.ItemIsEnabled
                if k == 1:
                    flags |= QtCore.Qt.ItemFlag.ItemIsEditable
                    item.setToolTip("The amount of every compound in this run, unless a compound has its "
                                    "own (compound table). Read from the run name where it says.")
                    if level is None:
                        item.setBackground(QtGui.QBrush(QtGui.QColor(*BAD, 70)))
                item.setFlags(flags)
                if k == 0:
                    item.setToolTip(str(e.path))
                self.run_table.setItem(r, k, item)
        self._updating = False

    def _level_edited(self, item) -> None:
        if self._updating or item.column() != 1:
            return
        run, text = item.data(QtCore.Qt.ItemDataRole.UserRole), item.text().strip()
        try:
            value = float(text) if text else None
        except ValueError:
            self.say(f"Could not read '{text}' as an amount")
            QtCore.QTimer.singleShot(0, self, self._fill_run_table)
            return

        def change(s):
            if value is None:
                s.levels.pop(run, None)
            else:
                s.levels[run] = value
        QtCore.QTimer.singleShot(0, self, lambda: self._edit(change))

    def set_level(self, run: str, value: float | None) -> None:
        def change(s):
            if value is None:
                s.levels.pop(run, None)
            else:
                s.levels[run] = float(value)
        self._edit(change)

    def _remove_selected_runs(self) -> None:
        rows = sorted({i.row() for i in self.run_table.selectedIndexes()})
        if rows:
            self.project.remove(rows)

    def _run_selected(self, row: int) -> None:
        if 0 <= row < len(self.project.entries):
            self.renderer.set_index(row)
            self._draw_overlay()
            self._show_compound_spectrum()

    def show_peak_spectrum(self, blobs, k: int, title: str) -> None:
        """Show peak `k` of `blobs` (full-ion spectrum) in the spectrum panel."""
        self.spectrum.set_spectrum(blobs.mz, np.asarray(blobs.spectra[k], dtype=float),
                                   f"{title}: t₁ {blobs.t1[k] / 60:.2f} min, t₂ {blobs.t2[k]:.2f} s")

    def _show_compound_spectrum(self) -> None:
        name, run = self.selected_compound(), self.current_run()
        if name is None or run is None:
            return
        i = self.current_index()
        blobs = self.peaks.get(self.project.entries[i])
        p = self.points.get((name, run))
        if blobs is None or p is None or p.blob is None or not 0 <= p.blob < len(blobs):
            self.spectrum.clear()
            self.spectrum.header.setText(f"{name}: no peak matched in {run}")
            return
        self.show_peak_spectrum(blobs, p.blob, f"{name} in {run}")

    def _spectrum_ion(self, mz: int, add: bool) -> None:
        """An ion clicked in the spectrum is imaged; Shift adds it to the ions shown."""
        text = self.ion.text().strip()
        self.ion.setText(f"{text}+{mz}" if (add and text) else str(mz))
        self._set_ion()

    def current_index(self) -> int:
        return self.renderer.index

    def current_run(self) -> str | None:
        i = self.current_index()
        return self.project.entries[i].name if 0 <= i < len(self.project.entries) else None

    def _remove_run(self) -> None:
        i = self.runs.currentRow()
        if 0 <= i < len(self.project.entries):
            self.project.remove([i])

    def _set_period(self) -> None:
        v = float(self.period.value())
        if abs(v - self.project.period) > 1e-9:
            self.project.set_period(v)
            self.project.estimate_phases()

    def _set_scale(self) -> None:
        self.view.set_scale(self.scale.currentData())
        if self.view.result is not None:
            self._view_levels(self.view.result)

    def _view_levels(self, result) -> None:
        """Colour limits that show the peaks against the cleaned background, not just the tallest.

        The shared view starts a cleaned image's log scale at its 90th percentile (most cleaned
        pixels are noise near zero), which on these sparse calibration runs leaves only the largest
        peaks in colour and looks like a linear scale. Here the limits run from the median of the
        pixels with any signal to the 99.9th percentile."""
        if not result.processed or self.view.scale != "log":
            return
        img = np.asarray(result.image, dtype=float)
        pos = img[np.isfinite(img) & (img > 0)]
        if pos.size < 10:
            return
        lo, hi = np.log10(np.percentile(pos, [50.0, 99.9]))
        if hi > lo:
            self.view.set_levels(float(lo), float(hi))

    def _set_ion(self) -> None:
        try:
            self.renderer.set_expression(self.ion.text().strip())
        except Exception as exc:  # noqa: BLE001  (an unreadable expression is the user's to fix)
            self.say(f"Ion: {exc}")

    # ======================================================== processing
    def edit_background_settings(self) -> None:
        dlg = SettingsDialog("🧹 Background settings", self.project.bg_settings, parent=self)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.project.bg_settings = dlg.settings()
            self.say("🧹 settings changed: press ▶ Clean and find blobs to use them")

    def edit_peak_settings(self) -> None:
        dlg = SettingsDialog("🔍 Blob settings", self.project.blobs.settings, parent=self)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.project.blobs.set_settings(dlg.settings())
            self.say("🔍 settings changed: press ▶ Clean and find blobs to use them")

    def edit_match_settings(self) -> None:
        dlg = SettingsDialog("Matching", self.session.settings, tips={
            "t1_window_s": "A compound matches the blob holding its position, else the nearest apex this close in t₁…",
            "t2_window_s": "…and this close in t₂",
            "t1_scale": "For the nearest internal standard: this many seconds of t₁…",
            "t2_scale": "…count the same as this many seconds of t₂",
            "amount_ratio": "x = amount / IS amount",
            "saturation": "Flag a point (triangle, amber Max point) whose largest raw value reaches this. "
                          "0 turns it off. The August 2026 standards levelled off near 3e5 counts."},
            parent=self)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            new = dlg.settings()
            self._edit(lambda s: setattr(s, "settings", new))

    def edit_curve_settings(self) -> None:
        dlg = SettingsDialog("Curve rules", self.session.calibration, parent=self)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            new = dlg.settings()
            self._edit(lambda s: setattr(s, "calibration", new))

    def run_pipeline(self, on_ready=None, line_up: bool | None = None) -> None:
        """🧹 then 🔍 on every run (both come from the cache when nothing changed), then,
        if "Line up phases" is on, the phases are lined up (`line_up_phases`) and any
        run whose phase moved is cleaned and searched again."""
        idx = self.project.loaded()
        if not idx:
            self.say("Open the standard runs first")
            return
        if self.project.blobs.running:
            self.say("Peak detection is already running")
            return
        line_up = self.act_line_up.isChecked() if line_up is None else line_up
        settings = self.project.bg_settings
        method = self.peak_method
        self.processing = True
        self._stage = "🧹 cleaning"
        self._update_guide()
        self._fill_run_table()

        def finished():
            self.processing = False
            self._blobs_done()
            if on_ready:
                on_ready()

        def components():
            # MS-aware: each run's fitted components (the main window's own, when it has them for
            # this run at this fold and cleaning), re-solved on every ion as the main window's
            # classify and identify steps read them, become the builder's peak list.
            p = self.project
            self.ms_peaks = {}
            todo = [e for e in p.entries if e.run is not None]
            remaining = set(todo)
            self._stage = "🔍 finding MS-aware peaks (a few minutes a run, less if the main window already did)"
            self._update_guide()
            for e in todo:
                i = p.entries.index(e)

                def done(result, e=e):
                    self.ms_peaks[e] = result
                    remaining.discard(e)
                    self._fill_run_table()
                    self._update_guide()
                    if not remaining:
                        finished()

                def failed(text, e=e):
                    self.say(f"{e.name}: MS-aware peaks failed: {text}")
                    remaining.discard(e)
                    if not remaining:
                        finished()
                self.jobs.submit(f"📈 MS-aware peaks: {e.name}",
                                 self._component_work(p, i, self._main_components(e)), done, failed)
            if not todo:
                finished()

        def peaks():
            if method == MS_AWARE:
                components()
            elif line_up:
                finished()                  # the legacy blobs used to line up are the peaks
            else:
                self._stage = "🔍 finding blobs"
                self._update_guide()
                self.project.blobs.run_all(on_done=finished, method=LEGACY)

        def lined_up():
            # Phases are lined up on legacy blobs, which take a minute a run and were checked on
            # real standards; MS-aware components then need finding only once, at the final fold.
            if self.line_up_phases():
                self.run_pipeline(on_ready, line_up=False)
                return
            peaks()

        def cleaned(results):
            self.project.set_backgrounds(results, settings, None)
            if line_up:
                self._stage = "🔍 finding blobs to line the phases up"
                self._update_guide()
                self.project.blobs.run_all(on_done=lined_up, method=LEGACY)
            else:
                peaks()

        self.say("🧹 Cleaning the standard runs…")
        self.project.compute_background(idx, settings, None, on_done=cleaned)

    def _main_components(self, entry):
        """The MS-aware peaks the main window already found on this run, as (component set, the
        run they were fitted on), if they are current and at this fold; else None. These are the
        peaks the person has seen and perhaps separated by hand, so they are used as they are,
        on the main window's own data (cleaned there or not)."""
        main = self.main_project
        if main is None or self.peak_method != MS_AWARE:
            return None
        i = main.index_of(entry.name)
        if i is None:
            return None
        m = main.entries[i]
        cs = main.blobs.component_sets_by_name.get(m.name)
        if cs is None or m.run is None or m.name in main.blobs.component_stale_names:
            return None
        # A sample is the main window's own run, at its own fold: always current. A standard is a
        # copy at the builder's fold, which must be the main window's. (Comparing a sample with the
        # builder's period sent every sample to a fresh search of several minutes.)
        if m is not entry and (Path(m.path).resolve() != Path(entry.path).resolve()
                               or abs(m.phase - entry.phase) > 1e-9
                               or abs(main.period - self.project.period) > 1e-9):
            return None
        return cs, main.base_run(i)

    def _component_work(self, project, i: int, preset=None, raw=None, fold=None, clean=None, bg_key=None):
        """A job making one run's MS-aware peak list: the given component set, or components found
        (cached under the builder's own key), re-solved on every ion. Standards pass their builder project index; samples pass the main project's
        index with their own raw run, fold and cleaning."""
        comp = self.project.blobs.component_settings
        wrap = comp.candidates.wrap_below_s
        e = project.entries[i]
        raw = e.run if raw is None else raw
        fold = project.fold(i) if fold is None else fold
        base = None if clean is not None else project.base_run(i)
        bg_key = project.base_key(i) if (bg_key is None and clean is None) else bg_key
        floor = storage_threshold(raw) / 3.0
        name = e.name
        # Legacy blobs of this run at this fold, when the line-up step found them: a check on
        # components taken from elsewhere.
        legacy = project.blobs.by_run.get(e) if project is self.project else None

        def work(token, progress):
            if preset is not None:
                cs, run = preset                 # the main window's own peaks, on its own data
                if legacy is not None:
                    peaks, _model, _footprint = respectrum_components(run, fold, cs, wrap_below_s=wrap)
                    taken = peaks_from_components(peaks, fold, wrap)
                    if peak_coverage(taken, legacy) >= MIN_COVERAGE:
                        return taken
                    # They do not sit on this run's peaks (found at another fold or by other code):
                    # find them again here.
                    cs, run = None, base if clean is None else clean()
            else:
                cs, run = None, clean() if clean is not None else base
            if cs is None:
                # The builder's own key: on 7 October 2026 the entry shared with the main window's 🔍
                # step for one real standard held components that matched its peaks far worse than a
                # fresh detection (64 of its 300 largest blobs against 153), origin unknown.
                key = _hash_key(bg_key, MS_AWARE, comp.to_dict(), floor, "calbuilder-1") if bg_key else None
                cs = project.cached(key, lambda: detect_components(
                    run, fold, name, comp, floor=floor, cancel=token.check,
                    progress=lambda done, total: progress(done, total, "Separating overlapping components")))[2]
            token.check()
            peaks, _model, _footprint = respectrum_components(run, fold, cs, wrap_below_s=wrap)
            return peaks_from_components(peaks, fold, wrap)
        return work

    @property
    def peaks(self) -> dict:
        """RunEntry -> the peak list the builder works on: MS-aware components or legacy blobs."""
        return self.ms_peaks if self.peak_method == MS_AWARE else self.project.blobs.by_run

    def set_peak_method(self, method: str) -> None:
        """MS-aware components or legacy blobs; the standards are searched again."""
        if method == self.peak_method:
            return
        self.peak_method = method
        self.project.blobs.method = method
        self.sample_blobs, self.sample_offsets, self.sample_results = {}, {}, []
        for a in (self.act_ms_aware, self.act_legacy):
            a.setChecked((a is self.act_ms_aware) == (method == MS_AWARE))
        self._outline_cache = None
        self.recompute()
        if self.project.loaded():
            self.run_pipeline()

    def line_up_phases(self, min_shift: float = 0.1, min_pairs: int = 10) -> bool:
        """Give every standard the same second-dimension position for the same compound.

        Each run's phase is estimated on its own and can come out tenths of a second
        apart (0.73 s for one of the four August 2026 standards), so blobs no longer
        line up between runs. `core.calbuilder.t2_offsets` measures each run's offset
        from blobs of the same spectrum; the runs are moved by their offset less the
        median offset, so the majority stays where it is. Only moves larger than
        `min_shift` are made. Returns whether any phase moved (blobs must then be found
        again)."""
        p = self.project
        # Always on legacy blobs (MS-aware components broke this on the real standards: every run
        # was moved, by up to 1.2 s), so a run's phase is settled before components are fitted.
        legacy = {e.name: p.blobs.by_run.get(e) for e in p.entries}
        names = [n for n, b in legacy.items() if b is not None]
        if len(names) < 2:
            return False
        sil = SiloxaneSettings()
        usable = {n: usable_blobs(legacy[n], p.period) & ~flag_siloxane_blobs(legacy[n], sil) for n in names}
        offsets = t2_offsets({n: legacy[n] for n in names}, names[0], usable, min_pairs=min_pairs)
        good = {n: o for n, (o, _pairs) in offsets.items() if np.isfinite(o)}
        if len(good) < 2:
            return False
        centre = float(np.median(list(good.values())))
        moves = {n: o - centre for n, o in good.items() if abs(o - centre) > min_shift}
        if not moves:
            return False
        if len(moves) * 2 > len(good):
            # The point is to bring a stray run or two back to the rest; when most runs disagree,
            # there is no "rest", and moving them would only make things worse.
            self.say("Phases not lined up: the runs disagree about where their peaks sit (" +
                     ", ".join(f"{n} {o:+.2f} s" for n, o in good.items()) + "). Check the phases in tab 1.")
            return False
        phases = {}
        for n, move in moves.items():
            i = p.index_of(n)
            e = p.entries[i]
            # A run whose t2 must rise by `move` starts its modulations `move` earlier.
            phases[i] = ((e.phase - move) % p.period, "lined up",
                         f"{e.phase_detail}; moved {-move:+.3f} s to line up with the other standards "
                         f"({offsets[n][1]} blob pairs)")
        for e in p.entries:
            e.background = None          # recomputed by the next pass, not by the fold change as well
        p.blobs.by_run = {}
        self.ms_peaks = {}
        p.set_phases(phases)
        self.processing = True               # set_phases refreshes the guide; this pass is not over
        self.say("Phases lined up: " + ", ".join(f"{n} {-m:+.2f} s" for n, m in moves.items())
                 + "; cleaning and finding blobs again")
        return True

    def usable_masks(self) -> dict:
        """Run name -> blobs worth matching across runs: no streaks, solvent front or siloxane."""
        p = self.project
        out = {}
        sil = SiloxaneSettings()
        for e in p.entries:
            b = self.peaks.get(e)
            if b is None or e.run is None:
                continue
            # No solvent-front window: on short runs a large early peak passes for a solvent front,
            # and on the real standards the streak rules alone left the same blobs.
            out[e.name] = usable_blobs(b, p.period) & ~flag_siloxane_blobs(b, sil)
        return out

    def _blobs_done(self) -> None:
        self._outline_cache = None
        self.heights = {}
        n = sum(len(b) for b in self.peaks.values())
        self.say(f"🔍 {n} blobs in {len(self.peaks)} runs; compounds matched")
        self._runs_changed()
        self.recompute()
        # Largest raw point in every blob, for the saturation check (a second or two a run).
        p = self.project
        for e in p.entries:
            b = self.peaks.get(e)
            if b is None or e.run is None:
                continue
            run, fold, name = e.run, p.fold(p.entries.index(e)), e.name

            def done(h, name=name, b=b):
                if self.blobs_by_run().get(name) is b:
                    self.heights[name] = h
                    self.recompute(refresh_table=False)
            self.jobs.submit(f"Peak heights: {name}", lambda _t, _p, run=run, fold=fold, b=b:
                             blob_heights(run, fold, b), done, quiet=True)

    def suggest_dialog(self) -> None:
        if not self.peaks:
            self.say("The standards are still being processed (tab 1)")
            return
        runs = self.session.runs
        # Level amounts: tab 1's, else what the compounds or run names say.
        known = {r: next((c.amounts[r] for c in self.session.analytes() if r in c.amounts), None) for r in runs}
        amounts = {r: self.session.levels.get(r, known[r] if known[r] is not None else amount_from_name(r))
                   for r in runs}
        dlg = SuggestDialog(runs, amounts, self.blobs_by_run(), self.suggest_settings, self, self.usable_masks())
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.suggest_settings = dlg.settings
            self.add_suggestions(dlg.chosen(), dlg.amounts())

    def add_suggestions(self, chosen: list, amounts: dict) -> None:
        """Add suggested compounds, each matched by hand to the blobs the suggestion
        found (by spectrum, which follows drift between days better than position)."""
        blobs = self.blobs_by_run()

        def change(s):
            # The amounts typed in the dialog are the runs' level amounts (tab 1), which every
            # suggested calibrant takes; an internal standard gets 1 in every run.
            s.levels.update({r: float(v) for r, v in amounts.items()})
            for g in chosen:
                is_istd = g.kind == "constant"
                c = s.add_compound(CalCompound(
                    name=("IS " if is_istd and not g.name else "") + g.label, t1=g.t1, t2=g.t2, is_istd=is_istd,
                    amounts={r: 1.0 for r in s.runs} if is_istd else {}))
                for run, j in g.blobs.items():
                    s.pin_blob(c.name, run, blobs[run], j)
        self._edit(change)
        n_is = sum(g.kind == "constant" for g in chosen)
        self.say(f"Added {len(chosen) - n_is} calibrants and {n_is} internal standards")

    def normalise(self) -> None:
        if not self.peaks:
            self.say("▶ Clean and find blobs first")
            return
        if not self.session.istds():
            self.say("Mark at least one compound as an internal standard (the IS column) first")
            return
        self.normalised = True
        self.recompute()
        self.tabs.setCurrentIndex(2)

    def blobs_by_run(self) -> dict:
        return {e.name: self.peaks.get(e) for e in self.project.entries}

    # ======================================================== state and undo
    def _edit(self, change, refresh_table: bool = True) -> None:
        """Apply a change to the session, keeping the old state for undo."""
        self._undo.append(self.session.to_dict())
        del self._undo[:-100]
        self._redo.clear()
        change(self.session)
        self.recompute(refresh_table)

    def undo(self) -> None:
        if not self._undo:
            self.say("Nothing to undo")
            return
        self._redo.append(self.session.to_dict())
        self._restore(self._undo.pop())

    def redo(self) -> None:
        if not self._redo:
            self.say("Nothing to redo")
            return
        self._undo.append(self.session.to_dict())
        self._restore(self._redo.pop())

    def _restore(self, d: dict) -> None:
        runs = self.session.runs
        self.session = CalSession.from_dict(d)
        self.session.runs = runs
        self.recompute()

    def recompute(self, refresh_table: bool = True) -> None:
        # Internal standards with a quantifier ion anchor the phase (🔄 Phase), as in the main window.
        anchors = [Anchor(c.name, c.quant_ion, c.t1, c.t2) for c in self.session.istds() if c.quant_ion]
        self.project.anchors = anchors or None
        blobs = self.blobs_by_run()
        # Divided by the internal standards as soon as there are blobs and an internal standard.
        if self.peaks and self.session.istds():
            self.normalised = True
        self.resolved = self.session.resolve(blobs)
        self.points = self.session.points(blobs, self.resolved, self.heights)
        self.fits = self.session.fits(self.points) if self.normalised else {}
        if refresh_table:
            self._fill_table()
        self.curves.set_data(self.session, self.points, self.fits, self.normalised)
        self._draw_overlay()
        if self.sample_blobs:
            self._requantify()
        self._update_guide()
        self.changed.emit()

    def _update_guide(self) -> None:
        """Tab titles tick off what is done; each tab's banner says what to do on it next."""
        p, s = self.project, self.session
        n_runs = len(p.loaded())
        processing = self.processing
        processed = bool(self.peaks) and not processing
        n_done = sum(1 for e in p.entries if self.peaks.get(e) is not None)
        kind = "MS-aware peaks" if self.peak_method == MS_AWARE else "legacy blobs"
        levels = sum(r in s.levels for r in s.runs)
        n_cmp, n_is = len(s.analytes()), len(s.istds())
        fits = [f for f in self.fits.values() if isinstance(f, CalibrationFit)]
        n_ok = sum(f.ok for f in fits)
        done = [processed and levels >= 3, n_cmp > 0 and n_is > 0, bool(fits), bool(self.sample_results)]
        for k, title in enumerate(["Standards", "Compounds", "Curves", "Apply to samples"]):
            self.tabs.setTabText(k, f"{'✓' if done[k] else k + 1}  {title}")

        # 1 Standards
        main_sel = len(self.main_selection())
        if not p.entries:
            hint = ("<b>Choose your calibration runs.</b> Select them in the main window's Runs table and press "
                    "<i>Use the runs selected in the main window</i>"
                    + (f" ({main_sel} selected now)" if main_sel else "")
                    + ", or open a folder or files here. They are cleaned and searched for blobs straight away.")
        elif processing:
            hint = (f"<b>⏳ Processing {n_runs} standards: {getattr(self, '_stage', 'working')}</b> "
                    f"({n_done} of {n_runs} done; the status column follows each run). "
                    + ("Fill in the missing level amounts (amber) meanwhile." if levels < len(s.runs) else
                       "Tab 2 opens when it finishes."))
        elif not processed:
            hint = ("<b>The standards are not processed.</b> Press <i>▶ Process the standards</i> "
                    "(below the table).")
        elif levels < len(s.runs):
            hint = ("<b>Give each run its level amount</b> (amber cells): the amount of every compound in that run, "
                    "unless a compound has its own in the compound table. Then go to <b>2 Compounds</b>.")
        else:
            hint = f"✓ Standards ready ({kind}). Go to <b>2 Compounds</b>."
        self.standards_banner.setText(hint)
        self.from_main_btn.setText(f"Use the {main_sel} runs selected in the main window" if main_sel
                                   else "Use runs from the main window…")
        self.from_main_btn.setEnabled(self.main_project is not None and bool(self.main_project.entries))

        # 2 Compounds
        if processing:
            hint = (f"⏳ The standards are being processed ({getattr(self, '_stage', 'working')}, {n_done} of "
                    f"{n_runs} runs done; MS-aware peaks take a few minutes a run). This tab fills in when it "
                    "finishes.")
        elif not processed:
            hint = "The standards are not processed yet: see tab 1."
        elif not s.compounds:
            hint = ("<b>Say which compounds are in the standards</b>: <i>🔎 Suggest from the standards</i> "
                    "finds them for you; or load a table; or click a peak with "
                    "<i>A click on the chromatogram: add a compound here</i>.")
        elif not n_is:
            hint = ("<b>Mark at least one internal standard</b>: tick its IS box in the table (or add one by clicking "
                    "with <i>add an internal standard here</i>).")
        else:
            missing = sum(1 for (c, r), pt in self.points.items() if pt.blob is None and s.compound(c)
                          and not s.compound(c).is_istd)
            hint = (f"✓ {n_cmp} compounds and {n_is} internal standards. Amber names are not found in the run shown; "
                    "select one and click its peak to fix it. Go to <b>3 Curves</b>."
                    + (f" ({missing} compound-run pairs not found)" if missing else ""))
        self.compounds_banner.setText(hint)

        # 3 Curves
        if not fits:
            hint = "No curves yet: they appear once tab 2 has compounds, an internal standard and level amounts."
        else:
            hint = (f"<b>{n_ok} of {len(fits)} curves pass</b> the back-calculation check (✅ in the list). For a ⚠️ one, "
                    "click its worst point: <i>Leave out</i> (Delete) drops it, <i>Show in chromatogram</i> lets you "
                    "click the right blob. Then go to <b>4 Apply to samples</b>.")
        self.curves_banner.setText(hint)

        # 4 Apply to samples
        if self.main_project is None or not self.main_project.entries:
            hint = "Open your samples in the main window (📂 there), then come back to this tab."
        elif not fits:
            hint = "Build the curves first (tabs 1 to 3)."
        else:
            warn = ""
            if abs(self.main_project.period - p.period) > 0.01:
                warn = (f" ⚠️ The main window's period ({self.main_project.period:.3f} s) differs from the standards' "
                        f"({p.period:.3f} s): the samples must be run with the same method.")
            if self.sample_results:
                hint = ("✓ Results below follow every change to the curves. Hover an amber cell for why it is "
                        "flagged; <i>📤 Export results</i> writes them as CSV. Tick more runs and ▶ Quantify to add "
                        "them." + warn)
            else:
                hint = ("<b>Tick the sample runs and press ▶ Quantify.</b> Each is cleaned and searched like the "
                        "standards; every compound is found by position and spectrum and read off its curve." + warn)
        self.apply.banner.setText(hint)

        state = "processing the standards…" if processing else (
            f"{n_runs} standards · {n_cmp} compounds · {n_is} IS · {n_ok}/{len(fits)} curves pass"
            + (f" · {len({r.run for r in self.sample_results})} samples quantified" if self.sample_results else ""))
        self.guide.setText(state)

    def event(self, e) -> bool:
        # Coming back from the main window: its run selection may have changed.
        if e.type() == QtCore.QEvent.Type.WindowActivate:
            self._update_guide()
        return super().event(e)

    def _tab_changed(self, index: int) -> None:
        if index == 3:
            self._refresh_samples()
        self._update_guide()

    # ======================================================== 4 apply to samples
    def _refresh_samples(self) -> None:
        main = self.main_project
        standards = {str(e.path) for e in self.project.entries}
        self.apply.set_samples(main.entries if main is not None else [], standards)

    def quantify_samples(self, names: list[str], on_done=None) -> None:
        """Clean the chosen main-window runs and find their blobs with this window's settings (from
        the shared cache where the main window already did the same), then read every compound off
        its curve. Results update by themselves when a curve changes afterwards."""
        main = self.main_project
        if main is None or not names:
            self.say("Tick some sample runs first")
            return
        if not self.fits:
            self.say("Build the curves first")
            return
        bg_settings, peak_settings = self.project.bg_settings, self.project.blobs.settings
        todo = []
        for name in names:
            i = main.index_of(name)
            if i is None or main.entries[i].run is None:
                continue
            if name in self.sample_blobs:
                continue
            todo.append(i)
        remaining = set(todo)

        def finish():
            self._requantify()
            self._update_guide()
            self.tabs.setCurrentWidget(self.apply)
            if on_done:
                on_done()

        if not todo:
            finish()
            return
        ref = next((n for n, b in self.blobs_by_run().items() if b is not None), None)
        ref_blobs = self.blobs_by_run().get(ref)
        ref_usable = self.usable_masks().get(ref)
        # In the main window's job list (bottom left), where the person is looking.
        jobs = getattr(main, "jobs", None) or self.jobs
        # A fresh MS-aware search already uses every core, so those runs go one at a time; runs
        # whose peaks the main window found are only re-solved and go at once.
        searches: list = []

        def next_search():
            if searches:
                jobs.submit(*searches.pop(0))

        for i in todo:
            e = main.entries[i]
            raw, fold = e.run, main.fold(i)
            raw_key = main.raw_key(i)
            bg_key = _hash_key(raw_key, "background", bg_settings.to_dict()) if raw_key else None
            floor = storage_threshold(raw) / 3.0

            ms_work, preset = None, None
            if self.peak_method == MS_AWARE:
                preset = self._main_components(e)

                def clean_fn(raw=raw, fold=fold, bg_key=bg_key):
                    return main.cached(bg_key, lambda: remove_background(raw, fold, bg_settings)).run
                ms_work = self._component_work(main, i, preset, raw=raw, fold=fold,
                                               clean=clean_fn, bg_key=bg_key)
            search = ms_work is not None and preset is None

            def work(token, progress, raw=raw, fold=fold, bg_key=bg_key, floor=floor, name=e.name, ms_work=ms_work):
                if ms_work is not None:
                    progress(0, 3, "MS-aware peaks")
                    blobs = ms_work(token, progress)
                else:
                    progress(0, 3, "Cleaning")
                    clean = main.cached(bg_key, lambda: remove_background(raw, fold, bg_settings))
                    token.check()
                    progress(1, 3, "Finding blobs")
                    blob_key = _hash_key(bg_key, LEGACY, peak_settings.to_dict(), floor) if bg_key else None
                    blobs = main.cached(blob_key, lambda: find_blobs(clean.run, fold, peak_settings, floor=floor))
                progress(2, 3, "Lining up with the standards")
                usable = {name: usable_blobs(blobs, fold.spec.period)}
                if ref_usable is not None:
                    usable[ref] = ref_usable
                off = t2_offsets({ref: ref_blobs, name: blobs}, ref, usable)[name][0] if ref_blobs is not None \
                    else float("nan")
                return blobs, off

            def done(result, name=e.name, i=i, search=search):
                self.sample_blobs[name], self.sample_offsets[name] = result
                remaining.discard(i)
                if search:
                    next_search()
                if not remaining:
                    finish()

            def failed(text, name=e.name, i=i, search=search):
                self.say(f"{name} failed: {text}")
                remaining.discard(i)
                if search:
                    next_search()
                if not remaining:
                    finish()

            if search:
                searches.append((f"📈 Quantify {e.name} (finding MS-aware peaks)", work, done, failed))
            else:
                jobs.submit(f"📈 Quantify {e.name}", work, done, failed)
        n_search = len(searches)
        next_search()
        self.say(f"Quantifying {len(todo)} sample runs; progress is in the main window's job list"
                 + (f". {n_search} have no MS-aware peaks from the main window's 🔍 step and are searched one at a "
                    "time (minutes each): run 🔍 there first to skip this" if n_search else ""))

    def _requantify(self) -> None:
        """Read the processed samples off the current curves (fast; run after every change)."""
        main = self.main_project
        if not self.sample_blobs or main is None or not self.fits:
            return
        ticked = set(self.apply.ticked()) or set(self.sample_blobs)
        samples = {n: b for n, b in self.sample_blobs.items() if n in ticked}
        templates = compound_templates(self.session, self.blobs_by_run(), self.resolved)
        self.sample_results = quantify_samples(self.session, self.fits, templates, samples, self.sample_offsets)
        volumes = {e.name: e.volume for e in main.entries if e.volume}
        self.apply.set_results(self.sample_results, volumes, self.session.unit)

    def export_results(self, path) -> Path:
        main = self.main_project
        volumes = {e.name: e.volume for e in main.entries if e.volume} if main is not None else {}
        return write_dicts(sample_rows(self.sample_results, self.session.unit, volumes), path)

    def export_results_dialog(self) -> None:
        if not self.sample_results:
            self.say("Quantify some samples first (tab 4)")
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export sample results",
                                                        str(Path(self._folder() or ".") / "calibrated_results.csv"),
                                                        "CSV (*.csv)")
        if path:
            self.say(f"Wrote {self.export_results(path).name}")

    def _show_in_main(self, run: str, t1: float, t2: float) -> None:
        mw = self.main_window
        go = getattr(mw, "_batch_navigate", None)
        if go is None:
            self.say("Open the main window to see it")
            return
        go(run, t1, t2)
        mw.raise_()
        mw.activateWindow()

    # ======================================================== compound table
    def _fill_table(self) -> None:
        s = self.session
        cur = self.selected_compound()
        self._updating = True
        self.table.setColumnCount(len(FIXED_COLUMNS) + len(s.runs))
        self.table.setHorizontalHeaderLabels(FIXED_COLUMNS + [f"{r}\namount" for r in s.runs])
        self.table.setRowCount(len(s.compounds))
        flags = QtCore.Qt.ItemFlag.ItemIsSelectable | QtCore.Qt.ItemFlag.ItemIsEnabled
        run = self.current_run()
        std_names = [c.name for c in s.istds()]
        for r, c in enumerate(s.compounds):
            p = self.points.get((c.name, run)) if run else None
            name = QtWidgets.QTableWidgetItem(c.name)
            state = "not matched in this run" if (p is None or p.blob is None) else (
                "matched by hand in this run" if p.pinned else "matched in this run")
            name.setToolTip(state)
            if p is not None and p.blob is None and self.peaks:
                name.setForeground(QtGui.QBrush(QtGui.QColor(*BAD)))
            items = [name, QtWidgets.QTableWidgetItem(f"{c.t1:.2f}"), QtWidgets.QTableWidgetItem(f"{c.t2:.3f}")]
            is_item = QtWidgets.QTableWidgetItem()
            is_item.setFlags(flags | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            is_item.setCheckState(QtCore.Qt.CheckState.Checked if c.is_istd else QtCore.Qt.CheckState.Unchecked)
            items.append(is_item)
            ion = QtWidgets.QTableWidgetItem("" if c.quant_ion is None else str(c.quant_ion))
            ion.setToolTip("Blank: the whole blob (every ion)")
            items.append(ion)
            for k, item in enumerate(items):
                if k != 3:
                    item.setFlags(flags | QtCore.Qt.ItemFlag.ItemIsEditable)
                item.setData(QtCore.Qt.ItemDataRole.UserRole, c.name)
                self.table.setItem(r, k, item)
            # Divide by: the nearest internal standard unless one is chosen. The combo in
            # each row is reused, so a refresh never leaves replaced widgets awaiting deletion.
            holder = QtWidgets.QTableWidgetItem()
            holder.setFlags(flags)
            holder.setData(QtCore.Qt.ItemDataRole.UserRole, c.name)
            self.table.setItem(r, 5, holder)
            combo = self.table.cellWidget(r, 5)
            if not isinstance(combo, QtWidgets.QComboBox):
                combo = QtWidgets.QComboBox()
                combo.currentIndexChanged.connect(
                    lambda _i, cb=combo: self._set_istd(cb.property("compound"), cb.currentData()))
                self.table.setCellWidget(r, 5, combo)
            combo.blockSignals(True)
            combo.clear()
            combo.setProperty("compound", c.name)
            combo.setEnabled(not c.is_istd)
            if c.is_istd:
                combo.addItem("(itself)", "")
            else:
                nearest = s.istd_for(c.name) if not c.istd else None
                combo.addItem(f"nearest ({nearest.name})" if nearest else "nearest", "")
                for n in std_names:
                    combo.addItem(n, n)
                combo.setCurrentIndex(max(combo.findData(c.istd), 0))
            combo.blockSignals(False)
            for j, run_name in enumerate(s.runs):
                v = c.amounts.get(run_name)
                level = s.levels.get(run_name) if not c.is_istd else None
                if v is not None:
                    item = QtWidgets.QTableWidgetItem(f"{v:g}")
                    item.setToolTip("This compound's own amount in this run (blank it to use the run's level)")
                elif level is not None:
                    # The run's level amount, shown grey: typing over it gives the compound its own.
                    item = QtWidgets.QTableWidgetItem(f"{level:g}")
                    item.setForeground(QtGui.QBrush(QtGui.QColor(*MISSING_COLOUR)))
                    item.setToolTip("The run's level amount (tab 1); type to give this compound its own")
                else:
                    item = QtWidgets.QTableWidgetItem("")
                item.setFlags(flags | QtCore.Qt.ItemFlag.ItemIsEditable)
                item.setData(QtCore.Qt.ItemDataRole.UserRole, c.name)
                self.table.setItem(r, len(FIXED_COLUMNS) + j, item)
            if c.name == cur:
                self.table.selectRow(r)
        self._updating = False

    def selected_compound(self) -> str | None:
        rows = self.table.selectionModel().selectedRows() if self.table.selectionModel() else []
        if not rows:
            return None
        item = self.table.item(rows[0].row(), 0)
        return item.data(QtCore.Qt.ItemDataRole.UserRole) if item else None

    def select_compound(self, name: str) -> None:
        for r in range(self.table.rowCount()):
            item = self.table.item(r, 0)
            if item is not None and item.data(QtCore.Qt.ItemDataRole.UserRole) == name:
                self.table.selectRow(r)
                return

    def _compound_selected(self) -> None:
        if not self._updating:
            self._draw_overlay()
            self._show_compound_spectrum()

    def _table_edited(self, item) -> None:
        if self._updating:
            return
        name = item.data(QtCore.Qt.ItemDataRole.UserRole)
        col, text = item.column(), item.text().strip()
        checked = item.checkState() == QtCore.Qt.CheckState.Checked
        # Applied once this handler has returned, because applying it rebuilds the table.
        QtCore.QTimer.singleShot(0, self, lambda: self.edit_cell(name, col, text, checked))

    def edit_cell(self, name: str, col: int, text: str, checked: bool = False) -> None:
        """Apply an edit of one cell of the compound table."""
        if self.session.compound(name) is None:
            return
        try:
            if col == 0:
                new = [name]
                self._edit(lambda s: new.__setitem__(0, s.rename(name, text)))
                self.select_compound(new[0])
            elif col in (1, 2):
                value = float(text)
                self._edit(lambda s: setattr(s.compound(name), "t1" if col == 1 else "t2", value))
            elif col == 3:
                self._edit(lambda s: setattr(s.compound(name), "is_istd", checked))
            elif col == 4:
                ion = int(round(float(text))) if text else None
                self._edit(lambda s: setattr(s.compound(name), "quant_ion", ion))
            elif col >= len(FIXED_COLUMNS):
                run = self.session.runs[col - len(FIXED_COLUMNS)]
                self._set_amount(name, run, float(text) if text else None)
        except (ValueError, IndexError):
            self.say(f"Could not read '{text}'")
            self._fill_table()

    def _set_istd(self, name: str, std: str) -> None:
        if not self._updating and name:
            QtCore.QTimer.singleShot(0, self, lambda: self._edit(
                lambda s: setattr(s.compound(name), "istd", std or "")))

    def _set_amount(self, name: str, run: str, value) -> None:
        def change(s):
            c = s.compound(name)
            if value is None:
                c.amounts.pop(run, None)
            else:
                c.amounts[run] = float(value)
        self._edit(change)

    def _fill_amount(self) -> None:
        name = self.selected_compound()
        if name is None:
            self.say("Select a compound in the table first")
            return
        value, ok = QtWidgets.QInputDialog.getDouble(self, "Same amount in every run",
                                                     f"Amount of {name} in every run", 1.0, 0.0, 1e12, 6)
        if ok:
            self._edit(lambda s: s.compound(name).amounts.update({r: float(value) for r in s.runs}))

    def _set_unit(self) -> None:
        text = self.unit.text().strip()
        if text != self.session.unit:
            self._edit(lambda s: setattr(s, "unit", text), refresh_table=False)

    def _add_empty_row(self) -> None:
        def change(s):
            s.add_compound(CalCompound(name=f"compound {len(s.compounds) + 1}", t1=0.0, t2=0.0))
        self._edit(change)
        self.select_compound(self.session.compounds[-1].name)

    def _remove_selected(self) -> None:
        name = self.selected_compound()
        if name is not None:
            self._edit(lambda s: s.remove_compound(name))

    # ======================================================== chromatogram clicks
    def _view_clicked(self, col: int, t2: float) -> None:
        i = self.current_index()
        mode = self.click_mode.currentData()
        if i < 0 or self.project.entries[i].run is None:
            return
        fold = self.project.fold(i)
        t1 = float(fold.t1[col] + 0.5 * fold.spec.period)
        run = self.project.entries[i].name
        blobs = self.peaks.get(self.project.entries[i])
        tight = CalBuilderSettings(t1_window_s=0.75 * fold.spec.period, t2_window_s=3 * fold.row_pitch)
        hit = match_blob(blobs, t1, t2, tight)
        snap = (float(blobs.t1[hit]), float(blobs.t2[hit])) if hit is not None else (t1, float(t2))
        if hit is not None:
            self.show_peak_spectrum(blobs, hit, f"Peak in {run}")
        if mode == "none":
            return
        if mode in ("new", "new_is"):
            is_istd = mode == "new_is"
            default = self.session.unique_name("IS" if is_istd else "compound")
            name, ok = QtWidgets.QInputDialog.getText(
                self, "New internal standard" if is_istd else "New compound", "Name", text=default)
            if not ok:
                return
            self._edit(lambda s: s.add_compound(CalCompound(name=name.strip() or default, t1=snap[0], t2=snap[1],
                                                            is_istd=is_istd)))
            self.select_compound(self.session.compounds[-1].name)
            return
        name = self.selected_compound()
        if mode == "pick" and name is None:
            # Nothing selected: a new row, named after the peak until it is renamed.
            label = "compound"
            if hit is not None:
                label = f"m/z {int(blobs.mz[int(np.argmax(blobs.spectra[hit]))])} @ {snap[0] / 60:.2f} min"
            self._edit(lambda s: s.add_compound(CalCompound(name=label, t1=snap[0], t2=snap[1])))
            # Left unselected, so the next click adds another row rather than moving this one.
            self.table.clearSelection()
            self.say(f"Added {self.session.compounds[-1].name}: t₁ {snap[0]:.1f} s, t₂ {snap[1]:.3f} s. "
                     "Rename it in the table; tick IS for an internal standard.")
            return
        if name is None:
            self.say("Select a compound in the table first, then click its peak")
            return
        if mode == "pick":
            def move(s):
                c = s.compound(name)
                c.t1, c.t2 = snap
                s.unpin(name, run)
            self._edit(move)
            # Step on to the next compound still without a position, so a loaded name list can be
            # placed by clicking its peaks in turn.
            names = self.session.names()
            later = names[names.index(name) + 1:]
            todo = [n for n in later if self.session.compound(n).t1 <= 0]
            if todo:
                self.select_compound(todo[0])
            else:
                self.table.clearSelection()
            self.say(f"{name}: t₁ {snap[0]:.1f} s, t₂ {snap[1]:.3f} s (from the peak clicked)"
                     + (f". Next: click the peak for {todo[0]}" if todo else ""))
            return
        # assign: a click away from every blob adds a blob there first
        if blobs is None:
            self.say("▶ Clean and find blobs first")
            return
        if hit is None and self.peak_method == MS_AWARE:
            # No component there: take the nearest one within the usual matching window.
            hit = match_blob(blobs, t1, float(t2), self.session.settings)
            if hit is None:
                self.say("No MS-aware peak near there. Switch to legacy blobs (Settings) to add one by clicking.")
                return
        if hit is None:
            try:
                self.project.blobs.add(i, t1, float(t2))
            except ValueError as exc:
                self.say(f"No blob added: {exc}")
                return
            self._outline_cache = None
            blobs = self.peaks.get(self.project.entries[i])
            hit = len(blobs) - 1
        self._edit(lambda s: s.pin_blob(name, run, blobs, hit))
        self.say(f"{name} in {run}: matched by hand")

    def _match_current(self, how: str) -> None:
        name, run = self.selected_compound(), self.current_run()
        if name is None or run is None:
            self.say("Select a compound in the table first")
            return
        self._set_match(name, run, how)

    def _set_match(self, name: str, run: str, how: str) -> None:
        if how == "none":
            self._edit(lambda s: s.pin(name, run, None))
        else:
            self._edit(lambda s: s.unpin(name, run))

    def _toggle_point(self, name: str, run: str) -> None:
        if run in self.session.excluded.get(name, {}):
            self._edit(lambda s: s.include(name, run), refresh_table=False)
        else:
            self._edit(lambda s: s.exclude(name, run, "left out by hand"), refresh_table=False)

    def _set_override(self, name: str, ov: dict) -> None:
        def change(s):
            if ov:
                s.overrides[name] = dict(ov)
            else:
                s.overrides.pop(name, None)
        self._edit(change, refresh_table=False)

    def _set_amount_ratio(self, on: bool) -> None:
        self._edit(lambda s: setattr(s, "settings", replace(s.settings, amount_ratio=bool(on))), refresh_table=False)

    def show_in_chromatogram(self, name: str, run: str) -> None:
        i = self.project.index_of(run)
        if i is None:
            return
        self.tabs.setCurrentIndex(1)
        self.runs.setCurrentRow(i)
        self.select_compound(name)
        self.click_mode.setCurrentIndex(self.click_mode.findData("assign"))
        p = self.points.get((name, run))
        c = self.session.compound(name)
        t1, t2 = (p.t1, p.t2) if (p is not None and p.blob is not None) else (c.t1, c.t2)
        QtCore.QTimer.singleShot(200, self, lambda: self.view.centre_on(t1, t2))
        self.say(f"Click the blob that is {name} in {run}, or 'No peak in this run'")

    # ======================================================== overlay
    def _draw_overlay(self) -> None:
        i = self.current_index()
        empty = (np.array([]), np.array([]))
        if not 0 <= i < len(self.project.entries):
            return
        e = self.project.entries[i]
        blobs = self.peaks.get(e)
        run = e.name
        selected = self.selected_compound()
        markers, matched = [], set()
        # With many compounds (a suggested set runs to dozens) only the selected one and
        # the internal standards are labelled; the rest keep their markers.
        crowded = len(self.session.compounds) > 25
        for c in self.session.compounds:
            p = self.points.get((c.name, run))
            colour = ISTD_COLOUR if c.is_istd else TARGET_COLOUR
            if c.name == selected:
                colour = SELECTED_COLOUR
            show = not crowded or c.is_istd or c.name == selected
            if p is not None and p.blob is not None:
                matched.add(p.blob)
                markers.append(dict(t1=p.t1, t2=p.t2, colour=colour, symbol="t" if c.is_istd else "o",
                                    label=(c.name + (" ✋" if p.pinned else "")) if show else ""))
            else:
                label = c.name + " (no match)" if blobs is not None else c.name
                markers.append(dict(t1=c.t1, t2=c.t2, colour=MISSING_COLOUR if c.name != selected else colour,
                                    symbol="x", label=label if show else ""))
        self.view.set_overlay(markers, empty)
        if blobs is None or len(blobs) == 0:
            self.view.set_blob_overlay(empty, empty)
            return
        key = (id(blobs), id(blobs.labels), len(blobs))
        if self._outline_cache is None or self._outline_cache[0] != key:
            self._outline_cache = (key, all_outlines(blobs), {})
        _key, outline, per_blob = self._outline_cache
        sel = sorted(k for k in matched if 0 <= k < len(blobs))
        for k in sel:
            if k not in per_blob:
                per_blob[k] = blobs.outline_times(k)
        highlight = ((np.concatenate([per_blob[k][0] for k in sel]), np.concatenate([per_blob[k][1] for k in sel]))
                     if sel else empty)
        if isinstance(blobs, ComponentBlobs):
            # Fitted contours drawn in magenta: the view's faint yellow outline vanished against
            # the yellow and green of the colour scale, so components looked missing or misplaced.
            self.view.set_blob_overlay(empty, highlight)
            if self.all_blobs.isChecked():
                self.view.set_component_overlay([dict(outline=outline, colour=CONTOUR_COLOUR, width=1.3)])
        else:
            self.view.set_blob_overlay(outline if self.all_blobs.isChecked() else empty, highlight)

    # ======================================================== files
    def load_compounds_dialog(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Compound table", self._folder(), "CSV (*.csv)")
        if path:
            self.load_compounds(path)

    def load_compounds(self, path) -> None:
        try:
            comps = load_compounds(path, self.session.runs)
        except (OSError, ValueError, KeyError) as exc:
            QtWidgets.QMessageBox.warning(self, "Compound table", f"Couldn't read {Path(path).name}: {exc}")
            return

        def change(s):
            for c in comps:
                old = s.compound(c.name)
                if old is not None:
                    s.remove_compound(c.name)
                s.compounds.append(c)
        self._edit(change)
        n_is = sum(c.is_istd for c in comps)
        self.say(f"{len(comps)} compounds ({n_is} internal standards) from {Path(path).name}")

    def save_compounds_dialog(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save compound table",
                                                        str(Path(self._folder()) / "calibration_compounds.csv"),
                                                        "CSV (*.csv)")
        if path:
            write_compounds(path, self.session.compounds, self.session.runs)
            self.say(f"Wrote {Path(path).name}")

    def _folder(self) -> str:
        if self.session_path is not None:
            return str(self.session_path.parent)
        return str(self.project.entries[0].path.parent) if self.project.entries else ""

    def export(self, folder) -> list[Path]:
        folder = Path(folder)
        blobs = self.blobs_by_run()
        fits = self.session.fits(self.points)
        paths = [write_compounds(folder / "calibration_compounds.csv", self.session.compounds, self.session.runs),
                 write_dicts(curve_rows(self.session, fits), folder / "calibration_curves.csv"),
                 write_dicts(point_rows(self.session, self.points, fits), folder / "calibration_points.csv")]
        rows = normalised_blob_rows(self.session, blobs, self.resolved)
        if rows:
            paths.append(write_dicts(rows, folder / "calibration_normalised_blobs.csv"))
        return paths

    def export_dialog(self) -> None:
        if not self.normalised:
            self.say("➗ Normalise to IS first")
            return
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Export calibration tables to", self._folder())
        if folder:
            paths = self.export(folder)
            self.say(f"Wrote {', '.join(p.name for p in paths)}")

    def session_extra(self) -> dict:
        p = self.project
        return dict(runs=[dict(path=str(e.path), phase=e.phase, phase_method=e.phase_method) for e in p.entries],
                    period=p.period, period_confirmed=p.period_confirmed,
                    background=p.bg_settings.to_dict(), peaks=p.blobs.settings.to_dict(),
                    blob_edits={r: [list(op) for op in v] for r, v in p.blobs.edits.items()},
                    normalised=self.normalised, peak_method=self.peak_method,
                    components=p.blobs.component_settings.to_dict())

    def save_session(self, path) -> Path:
        self.session_path = save_session(path, self.session, self.session_extra())
        self.say(f"Saved {self.session_path.name}")
        return self.session_path

    def save_session_dialog(self) -> None:
        default = self.session_path or Path(self._folder() or ".") / "calibration_session.json"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save calibration session", str(default), "JSON (*.json)")
        if path:
            self.save_session(path)

    def open_session_dialog(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open calibration session", self._folder(),
                                                        "JSON (*.json)")
        if path:
            self.open_session(path)

    def open_session(self, path, on_ready=None) -> None:
        """Read the runs again, restore their fold, then clean and find blobs (both come from the cache
        when nothing changed) and replay the session on top."""
        try:
            session, d = load_session(path)
        except (OSError, ValueError, KeyError) as exc:
            QtWidgets.QMessageBox.warning(self, "Calibration session", f"Couldn't open {Path(path).name}: {exc}")
            return
        self.session_path = Path(path)
        p = self.project
        if p.entries:
            p.remove(list(range(len(p.entries))))
        p.period = float(d.get("period", p.period))
        p.period_confirmed = bool(d.get("period_confirmed", True))
        self.period.setValue(p.period)
        p.bg_settings = BackgroundSettings.from_dict(d.get("background", {}))
        p.blobs.settings = PeakSettings.from_dict(d.get("peaks", {}))
        p.blobs.edits = {r: [tuple(op) for op in v] for r, v in d.get("blob_edits", {}).items()}
        if d.get("components"):
            from ..core.components import ComponentDetectionSettings
            p.blobs.component_settings = ComponentDetectionSettings.from_dict(d["components"])
            p.blobs.candidate_settings = p.blobs.component_settings.candidates
        method = d.get("peak_method", LEGACY)          # sessions saved before MS-aware peaks used blobs
        self.peak_method, p.blobs.method = method, method
        self.act_ms_aware.setChecked(method == MS_AWARE)
        self.act_legacy.setChecked(method == LEGACY)
        self.session = session
        self.normalised = bool(d.get("normalised", False))
        self._undo.clear()
        self._redo.clear()
        self.unit.setText(session.unit)
        runs = d.get("runs", [])
        phases = {Path(r["path"]).resolve(): (float(r["phase"]), r.get("phase_method", "manual")) for r in runs}
        missing = [r["path"] for r in runs if not Path(r["path"]).exists()]
        if missing:
            self.say(f"Missing: {', '.join(Path(m).name for m in missing)}")
        paths = [Path(r["path"]) for r in runs if Path(r["path"]).exists()]
        if not paths:
            self.recompute()
            return
        from .project import RunEntry
        new = [RunEntry(path=q, name=q.stem) for q in paths]
        p.entries.extend(new)
        p.runsChanged.emit()

        def loaded():
            for e in p.entries:
                ph = phases.get(e.path.resolve())
                if ph is not None:
                    e.phase, e.phase_method = ph[0] % p.period, ph[1]
            p.foldChanged.emit(-1)
            p.runsChanged.emit()
            self.run_pipeline(on_ready)

        p._read(new, loaded, estimate=False)

    def closeEvent(self, event) -> None:
        self.jobs.cancel_all()
        super().closeEvent(event)
