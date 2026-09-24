"""🫧 Siloxane panel.

The reference spectrum that flagging and removal fit against: a preset (the
measured column bleed is the default), edited by hand in the table, or taken
from any peak you click in the 2D view. The bar chart shows the reference as
it will be used. Blobs are also scored against a second reference (the
measured cyclic siloxane peak by default) and keep the larger score, because
bleed and discrete Dn peaks fragment differently.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from ...core.siloxane import PRESETS, SiloxaneSettings

CUSTOM = "Custom"


class SiloxanePanel(QtWidgets.QWidget):
    settingsChanged = QtCore.Signal()
    pickToggled = QtCore.Signal(bool)
    applyRequested = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._updating = False

        self.preset = QtWidgets.QComboBox()
        self.preset.addItems(list(PRESETS) + [CUSTOM])
        self.preset.setToolTip("Column bleed is the default, measured on the supplied real run. "
                               "Use the cyclic preset for a discrete Dn contaminant peak.")
        self.preset.activated.connect(self._preset_chosen)
        self.second = QtWidgets.QComboBox()
        self.second.addItem("None", None)
        for name in PRESETS:
            self.second.addItem(name, name)
        self.second.setToolTip("Blobs are scored against this reference too and keep the larger score. Column "
                               "bleed and cyclic siloxane peaks look different, so scoring against both catches "
                               "both kinds. Not used by the per-pixel method.")
        self.second.activated.connect(self._widgets_edited)
        self._second_custom: dict[int, float] = {}

        self.pick_btn = QtWidgets.QPushButton("🖱 Pick from 2D view")
        self.pick_btn.setCheckable(True)
        self.pick_btn.setToolTip("Click a siloxane peak (or the bleed) in the 2D view; its spectrum "
                                 "becomes the reference. Esc cancels.")
        self.pick_btn.toggled.connect(self.pickToggled)
        self.baseline_box = QtWidgets.QCheckBox("Subtract local baseline")
        self.baseline_box.setChecked(True)
        self.baseline_box.setToolTip("Take the modulation's median spectrum off the clicked point. "
                                     "On for a discrete peak; off when you click plain column bleed, "
                                     "which the baseline would otherwise remove.")

        self.table = QtWidgets.QTableWidget(0, 2)
        self.table.setHorizontalHeaderLabels(["m/z", "Relative intensity"])
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.table.itemChanged.connect(self._table_edited)
        add = QtWidgets.QPushButton("＋ Ion")
        add.clicked.connect(self._add_row)
        remove = QtWidgets.QPushButton("－ Remove selected")
        remove.clicked.connect(self._remove_rows)

        self.plot = pg.PlotWidget()
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.plot.setMouseEnabled(x=False, y=False)
        self.plot.setLabel("bottom", "m/z")
        self.plot.setLabel("left", "relative (%)")
        self.plot.setFixedHeight(130)
        self.bars = pg.BarGraphItem(x=[], height=[], width=1.2, brush=(209, 149, 224))
        self.plot.addItem(self.bars)

        self.method = QtWidgets.QComboBox()
        self.method.addItem("Block out siloxane blobs", "blobs")
        self.method.addItem("Fit and subtract per pixel", "pixels")
        self.method.setToolTip("Blobs: every 🔍 blob at or above the certainty below is zeroed. "
                               "Pixels: the older fit of the reference pattern at every pixel.")
        self.method.activated.connect(self._widgets_edited)

        # Certainty that a blob is siloxane, 0 to 100 %. Blobs at or above it are outlined
        # in dashed magenta in the 2D view (🔍 Blob outlines) and marked in the 🔍 table.
        self.threshold = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.threshold.setRange(0, 100)
        self.threshold.setPageStep(5)
        self.threshold.setToolTip("A blob's certainty is (share of its volume on these ions) × "
                                  "(cosine similarity to the reference). Blobs at or above the slider "
                                  "are outlined in dashed magenta and are what 🫧 blocks out. Provisional scale.")
        self.threshold_label = QtWidgets.QLabel("15%")
        self.threshold_label.setMinimumWidth(36)
        self.threshold.valueChanged.connect(lambda v: self.threshold_label.setText(f"{v}%"))
        self.blob_count = QtWidgets.QLabel("")
        self.blob_count.setStyleSheet("color: palette(mid);")
        self.block_all = QtWidgets.QCheckBox("Block every ion in the blob")
        self.block_all.setToolTip("On: the whole signal under a siloxane blob is removed. "
                                  "Off (default): only the ions of the two references are, so anything "
                                  "co-eluting inside the blob survives.")
        self.block_all.toggled.connect(self._widgets_edited)
        self.min_ions = QtWidgets.QSpinBox()
        self.min_ions.setRange(1, 30)
        self.min_ions.setToolTip("A blob or pixel needs at least this many of the ions to show signal")
        self.min_cosine = QtWidgets.QDoubleSpinBox()
        self.min_cosine.setRange(0.0, 1.0)
        self.min_cosine.setSingleStep(0.05)
        self.min_cosine.setDecimals(2)
        self.min_cosine.setToolTip("Removal only touches a pixel whose own ion ratios resemble the "
                                   "reference at least this closely")
        for w in (self.threshold, self.min_cosine):
            w.valueChanged.connect(self._widgets_edited)
        self.min_ions.valueChanged.connect(self._widgets_edited)

        self.apply_btn = QtWidgets.QPushButton("🫧 Apply to all runs")
        self.apply_btn.clicked.connect(self.applyRequested)
        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)
        self.status.setStyleSheet("color: palette(mid);")

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Reference"))
        row.addWidget(self.preset, 1)
        lay.addLayout(row)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Also score against"))
        row.addWidget(self.second, 1)
        lay.addLayout(row)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.pick_btn)
        row.addWidget(self.baseline_box)
        lay.addLayout(row)
        lay.addWidget(self.plot)
        lay.addWidget(self.table, 1)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(add)
        row.addWidget(remove)
        lay.addLayout(row)
        form = QtWidgets.QFormLayout()
        form.addRow("Method", self.method)
        srow = QtWidgets.QHBoxLayout()
        srow.addWidget(self.threshold, 1)
        srow.addWidget(self.threshold_label)
        form.addRow("Siloxane certainty ≥", srow)
        form.addRow("", self.blob_count)
        form.addRow("", self.block_all)
        form.addRow("Min. ions present", self.min_ions)
        form.addRow("Min. fit cosine (pixels)", self.min_cosine)
        lay.addLayout(form)
        lay.addWidget(self.apply_btn)
        lay.addWidget(self.status)
        self.set_settings(SiloxaneSettings())

    # ------------------------------------------------------------- reading
    def settings(self) -> SiloxaneSettings:
        ref = self._read_table()
        base = SiloxaneSettings.from_reference(ref) if ref else SiloxaneSettings()
        name = self.second.currentData()
        second = PRESETS[name] if name in PRESETS else (self._second_custom if name == CUSTOM else None)
        return SiloxaneSettings(ions=base.ions, reference=base.reference,
                                score_threshold=self.threshold.value() / 100.0,
                                min_ions_present=self.min_ions.value(),
                                min_fit_cosine=self.min_cosine.value(),
                                method=self.method.currentData(),
                                block_all_ions=self.block_all.isChecked()).with_second(second)

    def subtract_baseline(self) -> bool:
        return self.baseline_box.isChecked()

    def _read_table(self) -> dict[int, float]:
        out: dict[int, float] = {}
        for r in range(self.table.rowCount()):
            try:
                mz = int(round(float(self.table.item(r, 0).text())))
                v = float(self.table.item(r, 1).text())
            except (AttributeError, ValueError):
                continue
            if mz > 0 and v > 0:
                out[mz] = v
        return out

    # ------------------------------------------------------------- writing
    def set_settings(self, s: SiloxaneSettings, label: str | None = None) -> None:
        """Show `s` without emitting a change."""
        self._updating = True
        self._fill_table(s.reference_dict())
        self.threshold.setValue(int(round(100 * s.score_threshold)))
        self.min_ions.setValue(s.min_ions_present)
        self.min_cosine.setValue(s.min_fit_cosine)
        self.method.setCurrentIndex(max(0, self.method.findData(s.method)))
        self.block_all.setChecked(s.block_all_ions)
        second = s.second_dict()
        name = next((n for n, ref in PRESETS.items() if ref == second), CUSTOM if second else None)
        if name == CUSTOM:
            self._second_custom = second
            if self.second.findData(CUSTOM) < 0:
                self.second.addItem(CUSTOM, CUSTOM)
        self.second.setCurrentIndex(max(0, self.second.findData(name)))
        name = label or next((n for n, ref in PRESETS.items() if ref == s.reference_dict()), CUSTOM)
        self.preset.setCurrentText(name if name in PRESETS or name == CUSTOM else CUSTOM)
        self._draw(s.reference_dict())
        self._updating = False

    def set_reference(self, ref: dict[int, float], note: str = "") -> None:
        """A picked or edited reference: show it, mark the preset Custom and emit a change."""
        self._updating = True
        self._fill_table(ref)
        self.preset.setCurrentText(CUSTOM)
        self._draw(ref)
        self._updating = False
        self.status.setText(note)
        self.settingsChanged.emit()

    def set_picking(self, on: bool) -> None:
        self.pick_btn.blockSignals(True)
        self.pick_btn.setChecked(on)
        self.pick_btn.blockSignals(False)

    def set_blob_count(self, n_flagged: int | None, n_total: int) -> None:
        """How many of the current run's 🔍 blobs the slider currently marks as siloxane."""
        self.blob_count.setText("" if n_flagged is None else f"{n_flagged} of {n_total} blobs on this run")

    def set_stale_note(self, n_layers: int) -> None:
        self.status.setText(f"{n_layers} run{'s' if n_layers != 1 else ''} carry a 🫧 layer built with "
                            "earlier settings. Apply again to refresh." if n_layers else "")

    def _fill_table(self, ref: dict[int, float]) -> None:
        self.table.blockSignals(True)
        self.table.setRowCount(len(ref))
        for r, mz in enumerate(sorted(ref)):
            self.table.setItem(r, 0, QtWidgets.QTableWidgetItem(str(mz)))
            self.table.setItem(r, 1, QtWidgets.QTableWidgetItem(f"{ref[mz]:g}"))
        self.table.blockSignals(False)

    def _draw(self, ref: dict[int, float]) -> None:
        mz = np.array(sorted(ref), dtype=float)
        h = np.array([ref[int(m)] for m in mz])
        top = h.max() if h.size else 1.0
        self.bars.setOpts(x=mz, height=100.0 * h / top if h.size else [], width=1.2)

    # -------------------------------------------------------------- events
    def _preset_chosen(self, _i: int) -> None:
        name = self.preset.currentText()
        if name in PRESETS:
            self._updating = True
            self._fill_table(dict(PRESETS[name]))
            self._draw(dict(PRESETS[name]))
            self._updating = False
            self.status.setText("")
            self.settingsChanged.emit()

    def _table_edited(self, _item) -> None:
        if self._updating:
            return
        self.preset.setCurrentText(CUSTOM)
        self._draw(self._read_table())
        self.settingsChanged.emit()

    def _widgets_edited(self, *_a) -> None:
        if not self._updating:
            self.settingsChanged.emit()

    def _add_row(self) -> None:
        self._updating = True
        r = self.table.rowCount()
        self.table.insertRow(r)
        self.table.setItem(r, 0, QtWidgets.QTableWidgetItem(""))
        self.table.setItem(r, 1, QtWidgets.QTableWidgetItem(""))
        self._updating = False
        self.table.editItem(self.table.item(r, 0))

    def _remove_rows(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True)
        if not rows:
            return
        self.table.blockSignals(True)
        for r in rows:
            self.table.removeRow(r)
        self.table.blockSignals(False)
        self._table_edited(None)
