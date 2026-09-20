"""📈 Calibration panel.

Pick a compound to see its curve (response ratio against amount), the
back-calculation errors per level, and the fit statistics. Click a point
to select it; then exclude it (a reason is required and recorded),
include it again, or open its peak. The weighting, model and zero-intercept
choices can be overridden per compound. Every change is undoable.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from ...core import CalibrationFit

OK = (46, 158, 107)
BAD = (224, 160, 32)
GREY = (142, 134, 160)


class CalibrationPanel(QtWidgets.QWidget):
    overrideRequested = QtCore.Signal(str, dict)          # compound, overrides ({} clears)
    excludeRequested = QtCore.Signal(str, str, str)       # compound, run, reason ("" includes again)
    gotoRequested = QtCore.Signal(str, str)               # run, compound

    def __init__(self, parent=None):
        super().__init__(parent)
        self.report = None
        self.decisions = None
        self.selected: str | None = None       # run name of the selected point
        self._updating = False

        self.compounds = QtWidgets.QListWidget()
        self.compounds.setMaximumWidth(260)
        self.compounds.currentRowChanged.connect(lambda _r: self._show())
        self.cal_set = QtWidgets.QComboBox()
        self.cal_set.setToolTip("The curve for one calibration set. A sample without its own set "
                                "is quantified from a blend of the sets either side of it in the "
                                "sequence (see its 'cal set' column in Results), not shown here.")
        self.kind = QtWidgets.QComboBox()
        self.kind.addItem("Quantifier ion / IS", "quant")
        self.kind.addItem("Whole peak (TIC) / IS", "tic")
        self.kind.setToolTip("The TIC curve is used for semi-quantitative amounts of nearby peaks")
        self.weights = QtWidgets.QComboBox()
        for label, v in (("Auto", "auto"), ("None", "1"), ("1/x", "1/x"), ("1/x²", "1/x2")):
            self.weights.addItem(label, v)
        self.model = QtWidgets.QComboBox()
        for label, v in (("Auto (Mandel)", "auto"), ("Linear", "linear"), ("Quadratic", "quadratic")):
            self.model.addItem(label, v)
        self.zero = QtWidgets.QCheckBox("Through zero")
        for w in (self.cal_set, self.kind):
            w.currentIndexChanged.connect(lambda _i: self._show())
        for w in (self.weights, self.model):
            w.currentIndexChanged.connect(lambda _i: self._override())
        self.zero.toggled.connect(lambda _b: self._override())

        self.plot = pg.PlotWidget()
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.plot.setLabel("bottom", "amount")
        self.plot.setLabel("left", "response ratio")
        self.curve = self.plot.plot(pen=pg.mkPen((209, 149, 224), width=1.5))
        self.points = pg.ScatterPlotItem(size=9, hoverable=True)
        self.points.sigClicked.connect(self._point_clicked)
        self.plot.addItem(self.points)
        self.excl = pg.ScatterPlotItem(size=9, pen=pg.mkPen(GREY), brush=None, symbol="x")
        self.plot.addItem(self.excl)
        self.samples = pg.ScatterPlotItem(size=7, pen=None, brush=pg.mkBrush(126, 232, 250, 160),
                                          symbol="d")
        self.plot.addItem(self.samples)
        self.resid = pg.PlotWidget()
        self.resid.setMenuEnabled(False)
        self.resid.hideButtons()
        self.resid.setMaximumHeight(140)
        self.resid.setLabel("left", "error (%)")
        self.resid.setLabel("bottom", "amount")
        self.resid_points = pg.ScatterPlotItem(size=8)
        self.resid.addItem(self.resid_points)
        self.limits = [pg.InfiniteLine(angle=0, pen=pg.mkPen(GREY, style=QtCore.Qt.PenStyle.DashLine))
                       for _ in range(2)]
        for ln in self.limits:
            self.resid.addItem(ln)
        self.resid.addItem(pg.InfiniteLine(angle=0, pos=0, pen=pg.mkPen(GREY)))
        self.resid.setXLink(self.plot)
        self.stats = QtWidgets.QLabel("Run 🧪 Quantify to see calibrations")
        self.stats.setWordWrap(True)
        self.stats.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.point_label = QtWidgets.QLabel("Click a point to select it")
        self.exclude_btn = QtWidgets.QPushButton("Exclude…")
        self.include_btn = QtWidgets.QPushButton("Include")
        self.open_btn = QtWidgets.QPushButton("Open peak")
        self.exclude_btn.clicked.connect(self._exclude)
        self.include_btn.clicked.connect(self._include)
        self.open_btn.clicked.connect(self._open)
        self.log_axes = QtWidgets.QCheckBox("Log axes")
        self.log_axes.toggled.connect(self._log)

        controls = QtWidgets.QHBoxLayout()
        for w in (QtWidgets.QLabel("Set"), self.cal_set, self.kind, QtWidgets.QLabel("Weights"),
                  self.weights, QtWidgets.QLabel("Model"), self.model, self.zero, self.log_axes):
            controls.addWidget(w)
        controls.addStretch(1)
        prow = QtWidgets.QHBoxLayout()
        prow.addWidget(self.point_label, 1)
        for b in (self.exclude_btn, self.include_btn, self.open_btn):
            prow.addWidget(b)
        right = QtWidgets.QVBoxLayout()
        right.addLayout(controls)
        right.addWidget(self.plot, 3)
        right.addWidget(self.resid, 1)
        right.addLayout(prow)
        right.addWidget(self.stats)
        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addWidget(self.compounds)
        lay.addLayout(right, 1)
        self._update_point_buttons()

    def set_background(self, colour) -> None:
        self.plot.setBackground(colour)
        self.resid.setBackground(colour)

    # ---------------------------------------------------------------- data
    def set_report(self, report, decisions) -> None:
        self.report, self.decisions = report, decisions
        cur = self.current_compound()
        sets = sorted({cs for cs, _ in report.fits}) if report else []
        self._updating = True
        cs_cur = self.cal_set.currentData()
        self.cal_set.clear()
        for cs in sets:
            self.cal_set.addItem(cs or "(one set)", cs)
        self.cal_set.setCurrentIndex(max(self.cal_set.findData(cs_cur), 0))
        self.compounds.clear()
        if report:
            for t in report.targets:
                if t.is_istd:
                    continue
                fit = report.fit(t.compound, sets[0] if sets else "")
                ok = isinstance(fit, CalibrationFit) and fit.ok
                item = QtWidgets.QListWidgetItem(("✅ " if ok else "⚠️ ") + t.compound)
                item.setData(QtCore.Qt.ItemDataRole.UserRole, t.compound)
                self.compounds.addItem(item)
        self._updating = False
        names = [self.compounds.item(i).data(QtCore.Qt.ItemDataRole.UserRole)
                 for i in range(self.compounds.count())]
        if names:
            self.compounds.setCurrentRow(names.index(cur) if cur in names else 0)
        else:
            self._show()

    def current_compound(self) -> str | None:
        item = self.compounds.currentItem()
        return item.data(QtCore.Qt.ItemDataRole.UserRole) if item else None

    def select_compound(self, compound: str) -> None:
        for i in range(self.compounds.count()):
            if self.compounds.item(i).data(QtCore.Qt.ItemDataRole.UserRole) == compound:
                self.compounds.setCurrentRow(i)
                return

    def current_fit(self):
        c = self.current_compound()
        if self.report is None or c is None:
            return None
        return self.report.fit(c, self.cal_set.currentData() or "", self.kind.currentData())

    def _show(self) -> None:
        if self._updating:
            return
        c = self.current_compound()
        fit = self.current_fit()
        self._sync_overrides(c)
        for item in (self.curve, self.points, self.excl, self.samples, self.resid_points):
            item.setData([], [])
        if fit is None:
            self.stats.setText("Run 🧪 Quantify to see calibrations" if self.report is None else "")
            return
        if isinstance(fit, str):
            self.stats.setText(f"No curve: {fit}")
            return
        passes = fit.passes()
        brushes = [pg.mkBrush(*(OK if p else BAD)) for p in passes]
        self.points.setData(fit.x, fit.y, brush=brushes, data=fit.labels,
                            pen=[pg.mkPen("w", width=2) if lab == self.selected else None
                                 for lab in fit.labels])
        self.excl.setData(fit.x_excluded, fit.y_excluded)
        xs = np.linspace(0, max(fit.x.max(), fit.x_excluded.max() if fit.x_excluded.size else 0) * 1.05, 200)
        self.curve.setData(xs, fit.response(xs))
        self.resid_points.setData(fit.x, fit.rel_error, brush=brushes)
        lim = fit.settings.accept_re
        self.limits[0].setPos(lim)
        self.limits[1].setPos(-lim)
        # Samples, placed on the curve at their measured response.
        kind = self.kind.currentData()
        sx, sy = [], []
        for r in self.report.seq:
            if r.quant is None or r.is_calibration:
                continue
            y = r.quant.response(c, kind)
            if np.isfinite(y):
                x = float(fit.amount(np.array([y]))[0])
                if np.isfinite(x):
                    sx.append(x)
                    sy.append(y)
        self.samples.setData(sx, sy)
        ex = "".join(f"<br>Excluded {k}: {v}" for k, v in fit.excluded.items())
        infl = fit.influential()
        warn = f"<br>Worth a second look: {', '.join(infl)}" if infl else ""
        verdict = "passes" if fit.ok else "<b>fails</b> the back-calculation check"
        self.stats.setText(
            f"<b>{c}</b>: {verdict}.<br>{fit.summary()}<br>"
            f"Mandel p = {fit.mandel_p:.3g} (quadratic only if &lt; {fit.settings.alpha}); "
            f"intercept p = {fit.intercept_p:.3g}; weight scores "
            + ", ".join(f"{k}: {v:.1f}" for k, v in fit.weight_scores.items()) + ex + warn)

    def _sync_overrides(self, compound) -> None:
        ov = (self.decisions.overrides.get(compound, {}) if (self.decisions and compound) else {})
        self._updating = True
        self.weights.setCurrentIndex(max(self.weights.findData(ov.get("weights", "auto")), 0))
        self.model.setCurrentIndex(max(self.model.findData(ov.get("model", "auto")), 0))
        self.zero.setChecked(bool(ov.get("force_zero", False)))
        self._updating = False

    def _override(self) -> None:
        if self._updating:
            return
        c = self.current_compound()
        if c is None:
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

    # --------------------------------------------------------------- points
    def _point_clicked(self, _item, points, *_):
        if points is None or len(points) == 0:
            return
        self.selected = points[0].data()
        self._update_point_buttons()
        self._show()

    def _update_point_buttons(self) -> None:
        c = self.current_compound()
        excluded = (self.decisions.excluded.get(c, {}) if (self.decisions and c) else {})
        has = self.selected is not None
        self.point_label.setText(f"Selected: {self.selected}" if has else "Click a point to select it")
        self.exclude_btn.setEnabled(has and self.selected not in excluded)
        self.include_btn.setEnabled(bool(excluded))
        self.open_btn.setEnabled(has)

    def _exclude(self) -> None:
        c = self.current_compound()
        if not (c and self.selected):
            return
        reason, ok = QtWidgets.QInputDialog.getText(
            self, "Exclude point", f"Why leave {self.selected} out of the {c} curve?")
        if ok and reason.strip():
            self.excludeRequested.emit(c, self.selected, reason.strip())

    def _include(self) -> None:
        c = self.current_compound()
        excluded = self.decisions.excluded.get(c, {}) if self.decisions else {}
        run = self.selected if self.selected in excluded else next(iter(excluded), None)
        if c and run:
            self.excludeRequested.emit(c, run, "")

    def _open(self) -> None:
        c = self.current_compound()
        if c and self.selected:
            self.gotoRequested.emit(self.selected, c)
