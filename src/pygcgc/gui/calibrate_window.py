"""📈 Calibration builder: a separate window for a set of standard runs.

Work from left to right along its toolbar:

1. 📂 Open the standard runs (or take them from the main window). They are
   folded at the main window's period, or the one in a `method.json` beside
   them, and each run's phase is estimated.
2. Say which compounds are in them: load a compound table, or click peaks
   in the chromatogram ("Click: add a compound" or "add an internal
   standard"), then type each compound's amount per run in the table.
3. ▶ Clean and find blobs: 🧹 background removal and 🔍 TIC blobs on every
   run. Each compound is matched to a blob in each run automatically.
4. ➗ Normalise to IS: every compound's blob is divided by its internal
   standard's (the named one, otherwise the nearest) in the same run.
5. Curves: one per compound. Click a point to select it; leave it out or put
   it back, retype its amount, or "Show in chromatogram" and click the right
   blob to rematch it ("No peak here" removes it). Everything updates at once.

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

from ..core import BackgroundSettings, CalibrationFit, PeakSettings
from ..core.io import discover_runs
from ..core.phase import Anchor
from ..core.calbuilder import (CalBuilderSettings, CalCompound, CalSession, SuggestSettings,
                               amount_from_name, blob_heights, curve_rows, load_compounds, load_session,
                               match_blob, normalised_blob_rows, point_rows, save_session,
                               suggest_compounds, t2_offsets, usable_blobs, write_compounds)
from ..core.solvent import find_solvent_front
from ..core.siloxane import SiloxaneSettings, flag_siloxane_blobs
from ..core.quant import write_dicts
from .blob_state import LEGACY
from .dialogs import SettingsDialog
from .jobs import JobManager
from .project import Project
from .render import Renderer
from .theme import DARK, LIGHT
from .views.image2d import ChromatogramView

OK = (46, 158, 107)
BAD = (224, 160, 32)
GREY = (142, 134, 160)
ISTD_COLOUR = (60, 200, 110)
TARGET_COLOUR = (126, 232, 250)
SELECTED_COLOUR = (251, 146, 60)
MISSING_COLOUR = (150, 150, 150)
RUN_FILES = "Runs (*.cdf *.CDF *.lsc *.LSC);;All files (*)"
FIXED_COLUMNS = ["Compound", "t₁ (s)", "t₂ (s)", "IS", "Quant ion", "Divide by"]
CLICK_MODES = [("assign", "match the selected compound to this blob (this run)"),
               ("move", "move the selected compound here (every run)"),
               ("new", "add a compound here"),
               ("new_is", "add an internal standard here"),
               ("none", "nothing")]


def _fmt(v, digits=4) -> str:
    return "" if v is None or not np.isfinite(v) else f"{v:.{digits}g}"


def all_outlines(blobs) -> tuple[np.ndarray, np.ndarray]:
    """Every blob's outline at once, in `BlobSet.outline_times`' (start, end, NaN) format.

    The same pixel edges as calling `outline_times` blob by blob, found in one pass
    over the label image (an edge shared by two blobs is drawn once). Blob by blob
    took about 10 s for the 5,000 blobs of a real run; this takes a fraction of one."""
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
        self.resize(900, 640)
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


class CalibrationWindow(QtWidgets.QMainWindow):
    """📈 Build calibration curves from a separate set of standard runs."""
    changed = QtCore.Signal()

    def __init__(self, main_project: Project | None = None, dark: bool = False, parent=None):
        super().__init__(parent)
        self.setWindowTitle("📈 Calibration builder")
        self.setWindowFlag(QtCore.Qt.WindowType.Window, True)
        self.resize(1400, 900)
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
        self.project.blobs.method = LEGACY
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

        # ---------------------------------------------------------- widgets
        self.view = ChromatogramView()
        self.view.set_empty_message("Open the standard runs to start (📂 in this window's toolbar)")
        self.view.set_blob_overlay_visible(True)
        self.renderer.ready.connect(self.view.show_result)
        self.renderer.ready.connect(lambda _r: self._draw_overlay())
        self.view.pointClicked.connect(self._view_clicked)
        self.view.hovered.connect(lambda t1, t2, v: self.statusBar().showMessage(
            f"t₁ {t1 / 60:.2f} min   t₂ {t2:.3f} s   {v:.4g}", 3000))

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
        self.click_mode.setToolTip("What a click on the chromatogram does. Matching needs a compound selected "
                                   "in the table; a click away from every blob adds a blob there.")
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

        # ---------------------------------------------------------- layout
        view_bar = QtWidgets.QHBoxLayout()
        view_bar.addWidget(QtWidgets.QLabel("Ion"))
        view_bar.addWidget(self.ion)
        view_bar.addWidget(self.all_blobs)
        view_bar.addSpacing(12)
        view_bar.addWidget(QtWidgets.QLabel("Click:"))
        view_bar.addWidget(self.click_mode)
        view_bar.addWidget(self.none_here)
        view_bar.addWidget(self.auto_here)
        view_bar.addStretch(1)
        table_bar = QtWidgets.QHBoxLayout()
        for w in (self.add_row, self.del_row, self.fill_btn, QtWidgets.QLabel("Amounts in"), self.unit):
            table_bar.addWidget(w)
        table_bar.addStretch(1)
        upper = QtWidgets.QWidget()
        ul = QtWidgets.QVBoxLayout(upper)
        ul.setContentsMargins(0, 0, 0, 0)
        ul.addLayout(view_bar)
        ul.addWidget(self.view, 1)
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
        cl = QtWidgets.QHBoxLayout(chrom)
        cl.setContentsMargins(4, 4, 4, 4)
        left = QtWidgets.QVBoxLayout()
        left.addWidget(QtWidgets.QLabel("Standard runs"))
        left.addWidget(self.runs, 1)
        self.remove_run_btn = QtWidgets.QPushButton("Remove run")
        self.remove_run_btn.clicked.connect(self._remove_run)
        left.addWidget(self.remove_run_btn)
        cl.addLayout(left)
        cl.addWidget(split, 1)
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(chrom, "Chromatograms and compounds")
        self.tabs.addTab(self.curves, "Curves")
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
        tb = self.addToolBar("Calibration")
        tb.setMovable(False)
        tb.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextOnly)

        def act(text, slot, tip):
            a = QtGui.QAction(text, self)
            a.setToolTip(tip)
            a.triggered.connect(lambda _checked=False: slot())
            tb.addAction(a)
            return a

        self.act_open_folder = act("📂 Open folder…", self.open_folder_dialog,
                                   "Every run in a folder (the untreated 70 eV channel of native files)")
        self.act_open = act("Open runs…", self.open_runs_dialog, "Choose the standard runs one by one")
        self.act_from_main = act("From main window…", self.runs_from_main,
                                 "Use runs already open in the main window (calibration runs are ticked)")
        tb.addSeparator()
        tb.addWidget(QtWidgets.QLabel(" Period "))
        self.period = QtWidgets.QDoubleSpinBox()
        self.period.setRange(0.5, 30.0)
        self.period.setDecimals(4)
        self.period.setSingleStep(0.01)
        self.period.setSuffix(" s")
        self.period.setValue(self.project.period)
        self.period.setToolTip("Modulation period of the standard runs")
        self.period.editingFinished.connect(self._set_period)
        tb.addWidget(self.period)
        self.act_fold = act("🔄 Phase", lambda: self.project.estimate_phases(),
                            "Estimate each run's phase again: from the internal standards that have a "
                            "quantifier ion, otherwise from the quiet band")
        tb.addSeparator()
        self.act_compounds = act("Load compounds…", self.load_compounds_dialog,
                                 "A CSV of compound, t1_s, t2_s, role (istd), quant_ion, istd, amount "
                                 "and one amount column per run; a target_table.csv also works")
        self.act_save_compounds = act("Save compounds…", self.save_compounds_dialog,
                                      "Write the compound table, with the amounts, as a CSV")
        tb.addSeparator()
        self.act_bg_settings = act("🧹 Settings…", self.edit_background_settings, "Background removal settings")
        self.act_peak_settings = act("🔍 Settings…", self.edit_peak_settings, "Blob detection settings")
        self.act_line_up = act("Line up phases", lambda: None,
                               "After finding blobs, move any run whose second dimension sits apart from the "
                               "others (judged from blobs of the same spectrum) and find its blobs again")
        self.act_line_up.setCheckable(True)
        self.act_line_up.setChecked(True)
        self.act_run = act("▶ Clean and find blobs", self.run_pipeline,
                           "🧹 background removal then 🔍 TIC blobs on every run")
        self.act_suggest = act("🔎 Suggest compounds…", self.suggest_dialog,
                               "Find the calibrants (grow with the level) and constant peaks (internal "
                               "standards?) from the standards alone, with no compound table")
        self.act_normalise = act("➗ Normalise to IS", self.normalise,
                                 "Divide every compound's blob by its internal standard's in the same run")
        tb.addSeparator()
        self.act_match_settings = act("Matching…", self.edit_match_settings,
                                      "Search window for automatic matching and the nearest-IS distance scale")
        self.act_curve_settings = act("Curve rules…", self.edit_curve_settings,
                                      "Default weights, model and acceptance limits for every curve")
        tb.addSeparator()
        self.act_export = act("📤 Export…", self.export_dialog,
                              "Write the curves, points, compounds and every normalised blob as CSV tables")
        self.act_save = act("Save session…", self.save_session_dialog, "Everything needed to carry on later")
        self.act_load = act("Open session…", self.open_session_dialog, "Carry on with a saved session")

    def set_dark(self, dark: bool) -> None:
        c = DARK if dark else LIGHT
        self.view.set_background(c["plot"])
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

    def runs_from_main(self) -> None:
        entries = self.main_project.entries if self.main_project is not None else []
        if not entries:
            self.say("No runs are open in the main window")
            return
        dlg = _PickRunsDialog(entries, self)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted and dlg.paths():
            self.add_runs(dlg.paths())

    def add_runs(self, paths, on_loaded=None) -> None:
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
        self.project.add_files(paths, on_loaded)

    def _runs_changed(self) -> None:
        names = [e.name for e in self.project.entries]
        if names != self.session.runs:
            self.session.runs = names
            self._outline_cache = None
            self.recompute()
        cur = self.runs.currentRow()
        self.runs.blockSignals(True)
        self.runs.clear()
        for e in self.project.entries:
            b = self.project.blobs.by_run.get(e)
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

    def _run_selected(self, row: int) -> None:
        if 0 <= row < len(self.project.entries):
            self.renderer.set_index(row)
            self._draw_overlay()

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
            self.say("Blob detection is already running")
            return
        line_up = self.act_line_up.isChecked() if line_up is None else line_up
        settings = self.project.bg_settings

        def cleaned(results):
            self.project.set_backgrounds(results, settings, None)

            def found():
                if line_up and self.line_up_phases():
                    self.run_pipeline(on_ready, line_up=False)
                    return
                self._blobs_done()
                if on_ready:
                    on_ready()
            self.project.blobs.run_all(on_done=found, method=LEGACY)

        self.say("🧹 Cleaning the standard runs…")
        self.project.compute_background(idx, settings, None, on_done=cleaned)

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
        names = [e.name for e in p.entries if p.blobs.by_run.get(e) is not None]
        if len(names) < 2:
            return False
        blobs = self.blobs_by_run()
        offsets = t2_offsets({n: blobs[n] for n in names}, names[0], self.usable_masks(), min_pairs=min_pairs)
        good = {n: o for n, (o, _pairs) in offsets.items() if np.isfinite(o)}
        if len(good) < 2:
            return False
        centre = float(np.median(list(good.values())))
        moves = {n: o - centre for n, o in good.items() if abs(o - centre) > min_shift}
        if not moves:
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
        p.set_phases(phases)
        self.say("Phases lined up: " + ", ".join(f"{n} {-m:+.2f} s" for n, m in moves.items())
                 + "; cleaning and finding blobs again")
        return True

    def usable_masks(self) -> dict:
        """Run name -> blobs worth matching across runs: no streaks, solvent front or siloxane."""
        p = self.project
        out = {}
        sil = SiloxaneSettings()
        for e in p.entries:
            b = p.blobs.by_run.get(e)
            if b is None or e.run is None:
                continue
            solvent = find_solvent_front(e.run.times, np.asarray(e.run.tic, dtype=float))
            out[e.name] = usable_blobs(b, p.period, solvent) & ~flag_siloxane_blobs(b, sil)
        return out

    def _blobs_done(self) -> None:
        self._outline_cache = None
        self.heights = {}
        n = sum(len(b) for b in self.project.blobs.by_run.values())
        self.say(f"🔍 {n} blobs in {len(self.project.blobs.by_run)} runs; compounds matched")
        self._runs_changed()
        self.recompute()
        # Largest raw point in every blob, for the saturation check (a second or two a run).
        p = self.project
        for e in p.entries:
            b = p.blobs.by_run.get(e)
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
        if not self.project.blobs.by_run:
            self.say("▶ Clean and find blobs first")
            return
        runs = self.session.runs
        # Level amounts: what the compounds already say, else what the run names say.
        known = {r: next((c.amounts[r] for c in self.session.analytes() if r in c.amounts), None) for r in runs}
        amounts = {r: known[r] if known[r] is not None else amount_from_name(r) for r in runs}
        dlg = SuggestDialog(runs, amounts, self.blobs_by_run(), self.suggest_settings, self, self.usable_masks())
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.suggest_settings = dlg.settings
            self.add_suggestions(dlg.chosen(), dlg.amounts())

    def add_suggestions(self, chosen: list, amounts: dict) -> None:
        """Add suggested compounds, each matched by hand to the blobs the suggestion
        found (by spectrum, which follows drift between days better than position)."""
        blobs = self.blobs_by_run()

        def change(s):
            for g in chosen:
                is_istd = g.kind == "constant"
                c = s.add_compound(CalCompound(
                    name=("IS " if is_istd and not g.name else "") + g.label, t1=g.t1, t2=g.t2, is_istd=is_istd,
                    amounts={r: 1.0 for r in s.runs} if is_istd else dict(amounts)))
                for run, j in g.blobs.items():
                    s.pin_blob(c.name, run, blobs[run], j)
        self._edit(change)
        n_is = sum(g.kind == "constant" for g in chosen)
        self.say(f"Added {len(chosen) - n_is} calibrants and {n_is} internal standards")

    def normalise(self) -> None:
        if not self.project.blobs.by_run:
            self.say("▶ Clean and find blobs first")
            return
        if not self.session.istds():
            self.say("Mark at least one compound as an internal standard (the IS column) first")
            return
        self.normalised = True
        self.recompute()
        self.tabs.setCurrentWidget(self.curves)

    def blobs_by_run(self) -> dict:
        return {e.name: self.project.blobs.by_run.get(e) for e in self.project.entries}

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
        self.resolved = self.session.resolve(blobs)
        self.points = self.session.points(blobs, self.resolved, self.heights)
        self.fits = self.session.fits(self.points) if self.normalised else {}
        if refresh_table:
            self._fill_table()
        self.curves.set_data(self.session, self.points, self.fits, self.normalised)
        self._draw_overlay()
        self._update_guide()
        self.changed.emit()

    def _update_guide(self) -> None:
        loaded = bool(self.project.loaded())
        compounds = bool(self.session.compounds)
        istd = bool(self.session.istds())
        blobs = bool(self.project.blobs.by_run)
        stale = bool(self.project.blobs.stale)
        amounts = any(c.amounts for c in self.session.analytes())

        def step(done, text):
            return f"<span style='color:#2E9E6B'>✓ {text}</span>" if done else f"<b>{text}</b>"

        parts = [step(loaded, "1 📂 open the standard runs"),
                 step(compounds and istd, "2 add compounds and mark internal standards (table or clicks)"),
                 step(amounts, "3 type the amounts"),
                 step(blobs and not stale, "4 ▶ clean and find blobs" + (" (out of date)" if stale else "")),
                 step(self.normalised, "5 ➗ normalise to IS"),
                 "6 check each curve"]
        if self.project.blobs.running:
            parts[3] = "4 finding blobs…"
        self.guide.setText("  ·  ".join(parts))

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
            if p is not None and p.blob is None and self.project.blobs.by_run:
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
                item = QtWidgets.QTableWidgetItem("" if v is None else f"{v:g}")
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
        if i < 0 or mode == "none" or self.project.entries[i].run is None:
            return
        fold = self.project.fold(i)
        t1 = float(fold.t1[col] + 0.5 * fold.spec.period)
        run = self.project.entries[i].name
        blobs = self.project.blobs.by_run.get(self.project.entries[i])
        tight = CalBuilderSettings(t1_window_s=0.75 * fold.spec.period, t2_window_s=3 * fold.row_pitch)
        hit = match_blob(blobs, t1, t2, tight)
        snap = (float(blobs.t1[hit]), float(blobs.t2[hit])) if hit is not None else (t1, float(t2))
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
        if name is None:
            self.say("Select a compound in the table first, then click its peak")
            return
        if mode == "move":
            def move(s):
                c = s.compound(name)
                c.t1, c.t2 = snap
                s.unpin(name, run)
            self._edit(move)
            self.say(f"{name} moved to t₁ {snap[0]:.1f} s, t₂ {snap[1]:.3f} s")
            return
        # assign: a click away from every blob adds a blob there first
        if blobs is None:
            self.say("▶ Clean and find blobs first")
            return
        if hit is None:
            try:
                self.project.blobs.add(i, t1, float(t2))
            except ValueError as exc:
                self.say(f"No blob added: {exc}")
                return
            self._outline_cache = None
            blobs = self.project.blobs.by_run.get(self.project.entries[i])
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
        self.tabs.setCurrentIndex(0)
        self.runs.setCurrentRow(i)
        self.select_compound(name)
        self.click_mode.setCurrentIndex(0)
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
        blobs = self.project.blobs.by_run.get(e)
        run = e.name
        selected = self.selected_compound()
        markers, matched = [], set()
        for c in self.session.compounds:
            p = self.points.get((c.name, run))
            colour = ISTD_COLOUR if c.is_istd else TARGET_COLOUR
            if c.name == selected:
                colour = SELECTED_COLOUR
            if p is not None and p.blob is not None:
                matched.add(p.blob)
                markers.append(dict(t1=p.t1, t2=p.t2, colour=colour, symbol="t" if c.is_istd else "o",
                                    label=c.name + (" ✋" if p.pinned else "")))
            else:
                markers.append(dict(t1=c.t1, t2=c.t2, colour=MISSING_COLOUR if c.name != selected else colour,
                                    symbol="x", label=c.name + " (no match)" if blobs is not None else c.name))
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
                    normalised=self.normalised)

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
