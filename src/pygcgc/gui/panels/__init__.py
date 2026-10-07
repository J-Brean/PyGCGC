"""Dockable side panels: the run list, fold settings and layers."""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from ...core import period_from_anchors
from ...core.phase import PhaseEstimationError
from ..project import Project

METHOD_TEXT = {"anchors": "anchors", "void": "quiet band", "bleed": "bleed onset", "fallback": "borrowed",
               "manual": "by hand", "unset": "not set"}


class RunsPanel(QtWidgets.QWidget):
    """Table of runs. Selecting a row shows that run.

    Type, level, calibration set and volume can be typed straight into the
    table (double-click, or F2), or set for several selected runs at once
    from the right-click menu, which can also save them as a sequence.csv.
    A type containing "cal" plus a level makes a calibration run; a type
    containing "blank" is offered as a blank."""
    currentChanged = QtCore.Signal(int)
    removeRequested = QtCore.Signal(list)
    infoEdited = QtCore.Signal(list, str, object)      # run indices, field, value
    saveSequenceRequested = QtCore.Signal()

    COLUMNS = ("Run", "Type", "Level", "Cal set", "Volume (m³)", "Phase (s)", "Phase from", "Timing")
    EDITABLE = {1: "kind", 2: "level", 3: "cal_set", 4: "volume"}
    TYPES = ("calibration", "sample", "blank", "tube blank", "QC")

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self.table = QtWidgets.QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.DoubleClicked
                                   | QtWidgets.QAbstractItemView.EditTrigger.EditKeyPressed)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.table.customContextMenuRequested.connect(self._menu)
        self.table.currentCellChanged.connect(lambda r, *_: r >= 0 and self.currentChanged.emit(r))
        self.table.itemChanged.connect(self._item_changed)
        self.hint = QtWidgets.QLabel("PgUp / PgDn to flick between runs. Double-click Type, Level, Cal set or "
                                     "Volume to edit; right-click to set several runs at once.")
        self.hint.setWordWrap(True)
        self.hint.setStyleSheet("color: palette(mid);")
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addWidget(self.table)
        lay.addWidget(self.hint)
        project.runsChanged.connect(self.refresh)

    def refresh(self) -> None:
        current = self.table.currentRow()
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.project.entries))
        editable = QtCore.Qt.ItemFlag.ItemIsEditable
        for i, e in enumerate(self.project.entries):
            loaded = e.run is not None
            phase = f"{e.phase:.3f}" if loaded else ""
            status = e.error and "Unreadable" or (e.timing if loaded else "Loading…")
            volume = "" if e.volume is None else f"{e.volume:g}"
            values = (e.name, e.kind, e.level, e.cal_set, volume, phase,
                      METHOD_TEXT.get(e.phase_method, e.phase_method), status)
            if e.import_meta:
                energy = e.import_meta.get('energy_ev')
                processing = 'vendor DBC' if e.import_meta.get('processing_state') == 'vendor_dbc' else 'instrument'
                label = f"{e.name} [{energy:g} eV; {processing}]" if energy is not None else f"{e.name} [{processing}]"
                values = (label, *values[1:])
                if loaded and not self.project.period_confirmed:
                    values = (*values[:7], 'Check modulation period')
            for c, v in enumerate(values):
                item = self.table.item(i, c)
                if item is None:
                    item = QtWidgets.QTableWidgetItem()
                    self.table.setItem(i, c, item)
                item.setText(str(v))
                item.setFlags(item.flags() | editable if c in self.EDITABLE else item.flags() & ~editable)
                tip = e.phase_detail if c in (5, 6) else (e.error or str(e.path) if c == 0 else "")
                item.setToolTip(tip)
                item.setForeground(QtGui.QBrush(QtGui.QColor("#FB923C")) if (e.error or status == "Check timing")
                                   and c == 7 else QtGui.QBrush())
        if 0 <= current < self.table.rowCount():
            self.table.setCurrentCell(current, 0)
        self.table.blockSignals(False)

    def select(self, index: int) -> None:
        if 0 <= index < self.table.rowCount() and index != self.table.currentRow():
            self.table.setCurrentCell(index, 0)

    def selected_rows(self) -> list[int]:
        return sorted({i.row() for i in self.table.selectedIndexes()})

    def _item_changed(self, item) -> None:
        field = self.EDITABLE.get(item.column())
        if field is None:
            return
        text = item.text().strip()
        if field == "volume" and text:
            try:
                float(text)
            except ValueError:
                self.refresh()            # not a number: put the old value back
                return
        self.infoEdited.emit([item.row()], field, text)

    def _menu(self, pos) -> None:
        rows = self.selected_rows()
        if not rows:
            return
        n = f"{len(rows)} run{'s' if len(rows) > 1 else ''}"
        menu = QtWidgets.QMenu(self)
        types = menu.addMenu(f"Set type of {n}")
        for kind in self.TYPES:
            types.addAction(kind, lambda k=kind: self.infoEdited.emit(rows, "kind", k))
        other = [("level", "Set level…", "Calibration level (a row of the levels file):"),
                 ("cal_set", "Set calibration set…", "Calibration set (leave empty for none):"),
                 ("volume", "Set sampled volume…", "Sampled air volume in m³ (leave empty for none):")]
        for field, label, prompt in other:
            menu.addAction(label, lambda f=field, p=prompt, l=label: self._ask(rows, f, l, p))
        menu.addSeparator()
        menu.addAction("Save types and levels as sequence.csv…", self.saveSequenceRequested)
        menu.addSeparator()
        remove = menu.addAction(f"Remove {n} from project")
        chosen = menu.exec(self.table.viewport().mapToGlobal(pos))
        if chosen is remove:
            self.removeRequested.emit(rows)

    def _ask(self, rows, field, title, prompt) -> None:
        text, ok = QtWidgets.QInputDialog.getText(self, title.rstrip("…"), prompt)
        if not ok:
            return
        text = text.strip()
        if field == "volume" and text:
            try:
                float(text)
            except ValueError:
                return
        self.infoEdited.emit(rows, field, text)


class FoldPanel(QtWidgets.QWidget):
    """Modulation period and phase for the current run."""
    phaseEdited = QtCore.Signal(float)
    periodEdited = QtCore.Signal(float)
    fixPhaseRequested = QtCore.Signal()
    copyPhaseRequested = QtCore.Signal()
    periodApplyRequested = QtCore.Signal(float)
    estimatePeriodRequested = QtCore.Signal()

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self.index = -1
        self._suggested_period: float | None = None

        self.period = QtWidgets.QDoubleSpinBox()
        self.period.setDecimals(4)
        self.period.setRange(0.5, 30.0)
        self.period.setSingleStep(0.001)
        self.period.setSuffix(" s")
        self.period.setKeyboardTracking(False)
        self.period.valueChanged.connect(self.periodEdited)
        self.period_status = QtWidgets.QLabel("")
        self.period_status.setWordWrap(True)
        self.estimate_btn = QtWidgets.QPushButton("🔄 Estimate period")
        self.estimate_btn.setToolTip("Estimate the modulation period and every phase from the runs themselves, "
                                     "replacing the current period even if it was set by hand")
        self.estimate_btn.clicked.connect(self.estimatePeriodRequested)

        self.phase = QtWidgets.QDoubleSpinBox()
        self.phase.setDecimals(3)
        self.phase.setSingleStep(0.01)
        self.phase.setSuffix(" s")
        self.phase.setKeyboardTracking(False)
        self.phase.setWrapping(True)
        self.phase.valueChanged.connect(self._phase_spin)

        self.slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.slider.valueChanged.connect(self._phase_slider)

        self.method = QtWidgets.QLabel("")
        self.method.setWordWrap(True)
        fix = self.fix_btn = QtWidgets.QPushButton("🔄 Fix phase")
        fix.setToolTip("Estimate the phase of every run (internal standards if a target table is loaded)")
        fix.clicked.connect(self.fixPhaseRequested)
        copy = QtWidgets.QPushButton("Use for all runs")
        copy.setToolTip("Give every run this run's phase")
        copy.clicked.connect(self.copyPhaseRequested)
        self.check = QtWidgets.QPushButton("Check period")
        self.check.setToolTip("Fit the period from this run's internal standards")
        self.check.clicked.connect(self._check_period)
        self.apply_period = QtWidgets.QPushButton("Apply")
        self.apply_period.setVisible(False)
        self.apply_period.clicked.connect(self._apply_period)
        self.period_result = QtWidgets.QLabel("")
        self.period_result.setWordWrap(True)
        self.anchors = QtWidgets.QLabel("No target table loaded, so phases come from the quiet band")
        self.anchors.setWordWrap(True)
        self.anchors.setStyleSheet("color: palette(mid);")

        form = QtWidgets.QFormLayout()
        form.addRow("Modulation period", self.period)
        form.addRow("", self.period_status)
        form.addRow("Phase (this run)", self.phase)
        form.addRow("", self.slider)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addLayout(form)
        lay.addWidget(self.method)
        lay.addWidget(self.estimate_btn)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(fix)
        row.addWidget(copy)
        lay.addLayout(row)
        row2 = QtWidgets.QHBoxLayout()
        row2.addWidget(self.check)
        row2.addWidget(self.apply_period)
        row2.addStretch(1)
        lay.addLayout(row2)
        lay.addWidget(self.period_result)
        lay.addWidget(self.anchors)
        lay.addStretch(1)

        project.runsChanged.connect(self.refresh)
        project.foldChanged.connect(lambda _i: self.refresh())
        project.anchorsChanged.connect(self._anchors_changed)
        self.refresh()

    def set_index(self, index: int) -> None:
        self.index = index
        self.period_result.setText("")
        self.apply_period.setVisible(False)
        self.refresh()

    def refresh(self) -> None:
        p = self.project
        widgets = (self.period, self.phase, self.slider)
        for w in widgets:
            w.blockSignals(True)
        self.period.setValue(p.period)
        self.phase.setRange(0.0, p.period)
        self.slider.setRange(0, int(round(p.period * 100)))
        entry = p.entries[self.index] if 0 <= self.index < len(p.entries) else None
        enabled = entry is not None and entry.run is not None
        self.phase.setEnabled(enabled)
        self.slider.setEnabled(enabled)
        self.check.setEnabled(enabled and bool(p.anchors))
        self.estimate_btn.setEnabled(bool(p.loaded()))
        status = {"default": "Default: not estimated yet", "manual": "Set by hand"}.get(p.period_method, "")
        if p.period_method in ("file", "auto") or (p.period_method == "default" and p.period_detail):
            status = p.period_detail
        self.period_status.setText(status)
        # Flag an estimate whose precision is limited, or one that could not be made.
        flagged = "check the period" in status or "Could not" in status
        self.period_status.setStyleSheet("color: #b36b00;" if flagged else "color: palette(mid);")
        if entry is not None:
            self.phase.setValue(entry.phase)
            self.slider.setValue(int(round(entry.phase * 100)))
            text = {"anchors": "From internal standards", "void": "From the quiet band (no anchors found)",
                    "bleed": "From the start of the column bleed (no anchors found)",
                    "fallback": "Borrowed from the other runs (no anchors found)",
                    "manual": "Set by hand", "unset": "Not estimated yet"}.get(entry.phase_method, "")
            self.method.setText(text + (f"\n{entry.phase_detail}" if entry.phase_detail and
                                        entry.phase_method in ("anchors", "void", "bleed") else ""))
        else:
            self.method.setText("")
        for w in widgets:
            w.blockSignals(False)

    def _anchors_changed(self) -> None:
        a = self.project.anchors
        if a:
            names = ", ".join(x.name for x in a[:4]) + (" …" if len(a) > 4 else "")
            self.anchors.setText(f"{len(a)} anchors: {names}")
        self.refresh()

    def _phase_spin(self, v: float) -> None:
        self.slider.blockSignals(True)
        self.slider.setValue(int(round(v * 100)))
        self.slider.blockSignals(False)
        self.phaseEdited.emit(v)

    def _phase_slider(self, v: int) -> None:
        self.phase.blockSignals(True)
        self.phase.setValue(v / 100.0)
        self.phase.blockSignals(False)
        self.phaseEdited.emit(v / 100.0)

    def _check_period(self) -> None:
        p = self.project
        e = p.entries[self.index]
        try:
            est = period_from_anchors(e.run, p.anchors, p.period)
        except PhaseEstimationError as exc:
            self.period_result.setText(str(exc).capitalize())
            self.apply_period.setVisible(False)
            return
        corr = est.correction * 1e3
        se = est.correction_se * 1e3
        if est.significant:
            self.period_result.setText(f"Best fit {est.period:.5f} s ({corr:+.3f} ± {se:.3f} ms). "
                                       "That differs from the current period.")
            self._suggested_period = est.period
            self.apply_period.setVisible(True)
        else:
            self.period_result.setText(f"Period confirmed ({corr:+.3f} ± {se:.3f} ms, not significant).")
            self.apply_period.setVisible(False)

    def _apply_period(self) -> None:
        if self._suggested_period is not None:
            self.periodApplyRequested.emit(self._suggested_period)
            self.apply_period.setVisible(False)


class LayersPanel(QtWidgets.QWidget):
    """Processing layers, plus every other on/off overlay the 2D view can draw,
    gathered in one place. 🧹 background removal is switched on and off here (or
    with Space). The blob, stencil and cluster checkboxes mirror the toggles in
    the top bar and the Cluster panel: either one moves the other.
    Later layers are listed so it's clear where they will appear."""

    FUTURE = []

    def __init__(self, project: Project, processed_action: QtGui.QAction,
                aligned_action: QtGui.QAction, parent=None):
        super().__init__(parent)
        self.project = project
        self.action = processed_action
        self.aligned_action = aligned_action
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        self.raw = QtWidgets.QCheckBox("Raw data")
        self.raw.setChecked(True)
        self.raw.setToolTip("The base chromatogram image (raw or 🧹 cleaned, per the toggle below). "
                            "Hiding it leaves any overlays on show.")
        lay.addWidget(self.raw)
        alpha_row = QtWidgets.QHBoxLayout()
        alpha_row.setContentsMargins(22, 0, 0, 0)
        alpha_row.addWidget(QtWidgets.QLabel("Alpha"))
        self.raw_alpha = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.raw_alpha.setRange(0, 100)
        self.raw_alpha.setValue(100)
        self.raw_alpha.setToolTip("Opacity of the chromatogram image")
        alpha_row.addWidget(self.raw_alpha)
        lay.addLayout(alpha_row)
        self.background = QtWidgets.QCheckBox("🧹 Background removed   (Space)")
        self.background.setToolTip("Show the cleaned data. Apply 🧹 first.")
        self.background.toggled.connect(processed_action.setChecked)
        processed_action.toggled.connect(self._sync)
        processed_action.changed.connect(self._sync)
        lay.addWidget(self.background)
        self.count = QtWidgets.QLabel("")
        self.count.setStyleSheet("color: palette(mid); margin-left: 22px;")
        lay.addWidget(self.count)
        self.aligned = QtWidgets.QCheckBox("📍 Aligned   (see the Align panel)")
        self.aligned.toggled.connect(aligned_action.setChecked)
        aligned_action.toggled.connect(self._sync)
        aligned_action.changed.connect(self._sync)
        lay.addWidget(self.aligned)
        self.siloxane = QtWidgets.QCheckBox("🫧 Siloxanes removed")
        self.siloxane.setEnabled(False)          # shown under the same raw/cleaned toggle as 🧹
        self.siloxane.setToolTip("On top of 🧹, where it's been applied (toolbar 🫧 button)")
        lay.addWidget(self.siloxane)
        self.siloxane_count = QtWidgets.QLabel("")
        self.siloxane_count.setStyleSheet("color: palette(mid); margin-left: 22px;")
        lay.addWidget(self.siloxane_count)
        self.blob_overlay = QtWidgets.QCheckBox("🔍 Detected peaks")
        self.blob_overlay.setToolTip("Outline every detected peak (also the top-bar 🔍 toggle)")
        lay.addWidget(self.blob_overlay)
        self.stencil_overlay = QtWidgets.QCheckBox("⛓️ Stencil outlines")
        self.stencil_overlay.setToolTip("Every hydrocarbon stencil, lifted onto this run's own alkane series "
                                        "(also the top-bar ⛓️ toggle)")
        lay.addWidget(self.stencil_overlay)
        self.cluster_fill = QtWidgets.QCheckBox("🧩 Cluster fill")
        self.cluster_fill.setToolTip("Fill peaks with their cluster colour (also in the Cluster panel)")
        lay.addWidget(self.cluster_fill)
        self.cluster_outline = QtWidgets.QCheckBox("🧩 Cluster outline")
        self.cluster_outline.setToolTip("Outline peaks by cluster colour (also in the Cluster panel)")
        lay.addWidget(self.cluster_outline)
        for name, when in self.FUTURE:
            box = QtWidgets.QCheckBox(name)
            box.setEnabled(False)
            box.setToolTip(f"Arrives in {when}")
            lay.addWidget(box)
        lay.addStretch(1)
        project.layersChanged.connect(self._sync)
        project.runsChanged.connect(self._sync)
        self._sync()

    def _sync(self, *_):
        self.background.blockSignals(True)
        self.background.setEnabled(self.action.isEnabled())
        self.background.setChecked(self.action.isChecked())
        self.background.blockSignals(False)
        n = sum(e.background is not None for e in self.project.entries)
        total = len(self.project.entries)
        self.count.setText(f"on {n} of {total} runs" if n else "")
        self.aligned.blockSignals(True)
        self.aligned.setEnabled(self.aligned_action.isEnabled())
        self.aligned.setChecked(self.aligned_action.isChecked())
        self.aligned.blockSignals(False)
        n_sil = sum(e.siloxane is not None for e in self.project.entries)
        self.siloxane.setChecked(n_sil > 0)
        self.siloxane_count.setText(f"on {n_sil} of {total} runs" if n_sil else "")
