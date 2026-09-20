"""🧹 Background panel: settings, preview, apply, and diagnostics.

Changing any setting with Preview on recomputes the current run in the
background and shows it, so you can judge the effect before applying it
to every run. The diagnostics show what's left of the background for one
ion: after a good correction it should look like noise centred on zero,
with a spread matching the noise map (histogram), and lie along the
straight line in the normal quantile plot.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets
from scipy.stats import norm

from ...core import BackgroundSettings

STIFFNESS_STEPS = [1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]


class BackgroundPanel(QtWidgets.QWidget):
    settingsChanged = QtCore.Signal()
    applyAllRequested = QtCore.Signal()
    applyThisRequested = QtCore.Signal()
    removeAllRequested = QtCore.Signal()
    previewToggled = QtCore.Signal(bool)
    splitToggled = QtCore.Signal(bool)
    diagnosticIonChanged = QtCore.Signal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._building = True
        d = BackgroundSettings()

        self.percentile = self._dspin(1, 50, 1, 0, " %", d.percentile,
                                      "Low value taken down each modulation. Raise it if the "
                                      "background is left too high; lower it in very crowded regions.")
        self.stiffness = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.stiffness.setRange(0, len(STIFFNESS_STEPS) - 1)
        self.stiffness.setValue(STIFFNESS_STEPS.index(d.stiffness))
        self.stiffness.setToolTip("How stiff the baseline is along t₁. Stiffer follows only slow drift.")
        self.stiffness_label = QtWidgets.QLabel()
        self.pattern = QtWidgets.QCheckBox("Remove repeating t₂ pattern")
        self.pattern.setChecked(d.remove_pattern)
        self.pattern.setToolTip("Stripes that repeat every modulation (modulator pulses, "
                                "second-column bleed). Peaks spanning a few modulations are left alone.")
        self.window = self._spin(5, 401, 2, d.pattern_window, " modulations",
                                 "Width of the rolling window. It should be much wider than any peak.")
        self.despike = QtWidgets.QCheckBox("Remove single-scan spikes")
        self.despike.setChecked(d.despike)
        self.spike_threshold = self._dspin(3, 50, 0.5, 1, " SD", d.spike_threshold,
                                           "How far above its neighbouring scans a point must stand.")
        self.noise_block = self._spin(5, 200, 5, d.noise_block, " modulations",
                                      "Length of each stretch of t₁ in the noise map.")
        self.blank = QtWidgets.QComboBox()
        self.blank.setToolTip("A blank run to subtract after cleaning (tube or system blank).")
        self.scale_blank = QtWidgets.QCheckBox("Scale blank to each run (per ion)")
        self.scale_blank.setChecked(d.scale_blank)
        self.scale_blank.setToolTip("Artefacts vary between tubes. Per-ion scaling can over- or "
                                    "under-correct when two artefacts share an ion; per-peak "
                                    "scaling arrives with blob detection.")

        self.preview = QtWidgets.QCheckBox("Preview on this run")
        self.preview.setToolTip("Recompute the current run whenever a setting changes")
        self.preview.toggled.connect(self.previewToggled)
        self.split = QtWidgets.QCheckBox("Before / after split")
        self.split.setToolTip("Drag the white line on the image to compare")
        self.split.toggled.connect(self.splitToggled)
        apply_all = self.apply_all_btn = QtWidgets.QPushButton("🧹 Apply to all runs")
        apply_all.clicked.connect(self.applyAllRequested)
        apply_this = QtWidgets.QPushButton("Apply to this run")
        apply_this.clicked.connect(self.applyThisRequested)
        remove_all = QtWidgets.QPushButton("Remove from all runs")
        remove_all.clicked.connect(self.removeAllRequested)
        reset = QtWidgets.QPushButton("Reset settings")
        reset.clicked.connect(self.reset)
        self.summary = QtWidgets.QLabel("Not applied yet")
        self.summary.setWordWrap(True)

        form = QtWidgets.QFormLayout()
        form.addRow("Low percentile", self.percentile)
        srow = QtWidgets.QHBoxLayout()
        srow.addWidget(self.stiffness, 1)
        srow.addWidget(self.stiffness_label)
        form.addRow("Baseline stiffness", srow)
        form.addRow(self.pattern)
        form.addRow("Pattern window", self.window)
        form.addRow(self.despike)
        form.addRow("Spike threshold", self.spike_threshold)
        form.addRow("Noise map blocks", self.noise_block)
        form.addRow("Subtract blank", self.blank)
        blank_note = QtWidgets.QLabel('The blank is read through the alignment (step 2) when both runs have '
                                      'anchors; otherwise original coordinates are used. Align first.')
        blank_note.setWordWrap(True)
        form.addRow(blank_note)
        form.addRow(self.scale_blank)

        # Diagnostics
        self.diag_ion = QtWidgets.QSpinBox()
        self.diag_ion.setRange(1, 2000)
        self.diag_ion.setValue(207)
        self.diag_ion.setToolTip("Ion to check. The worst-affected ions are listed below.")
        self.diag_ion.valueChanged.connect(self.diagnosticIonChanged)
        self.diag_text = QtWidgets.QLabel("")
        self.diag_text.setWordWrap(True)
        self.hist_plot = pg.PlotWidget()
        self.hist_plot.setMenuEnabled(False)
        self.hist_plot.hideButtons()
        self.hist_plot.setLabel("bottom", "residual")
        self.hist_plot.setMinimumHeight(140)
        self.hist_bars = pg.PlotDataItem(stepMode="center", fillLevel=0,
                                         brush=pg.mkBrush(209, 149, 224, 150), pen=None)
        self.hist_fit = pg.PlotDataItem(pen=pg.mkPen((126, 232, 250), width=1.5))
        self.hist_plot.addItem(self.hist_bars)
        self.hist_plot.addItem(self.hist_fit)
        self.qq_plot = pg.PlotWidget()
        self.qq_plot.setMenuEnabled(False)
        self.qq_plot.hideButtons()
        self.qq_plot.setLabel("bottom", "normal quantile")
        self.qq_plot.setLabel("left", "residual / SD")
        self.qq_plot.setMinimumHeight(140)
        self.qq_points = pg.ScatterPlotItem(size=3, pen=None, brush=pg.mkBrush(209, 149, 224, 160))
        self.qq_line = pg.PlotDataItem([-3, 3], [-3, 3], pen=pg.mkPen((126, 232, 250), width=1.2))
        self.qq_plot.addItem(self.qq_points)
        self.qq_plot.addItem(self.qq_line)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addLayout(form)
        toggles = QtWidgets.QHBoxLayout()
        toggles.addWidget(self.preview)
        toggles.addWidget(self.split)
        lay.addLayout(toggles)
        grid = QtWidgets.QGridLayout()
        grid.addWidget(apply_all, 0, 0)
        grid.addWidget(apply_this, 0, 1)
        grid.addWidget(remove_all, 1, 0)
        grid.addWidget(reset, 1, 1)
        lay.addLayout(grid)
        lay.addWidget(self.summary)
        box = QtWidgets.QGroupBox("Diagnostics")
        dl = QtWidgets.QVBoxLayout(box)
        irow = QtWidgets.QHBoxLayout()
        irow.addWidget(QtWidgets.QLabel("Ion m/z"))
        irow.addWidget(self.diag_ion)
        irow.addStretch(1)
        dl.addLayout(irow)
        plots = QtWidgets.QHBoxLayout()
        plots.addWidget(self.hist_plot)
        plots.addWidget(self.qq_plot)
        dl.addLayout(plots)
        dl.addWidget(self.diag_text)
        lay.addWidget(box, 1)

        for w in (self.percentile, self.window, self.spike_threshold, self.noise_block):
            w.valueChanged.connect(self._changed)
        self.stiffness.valueChanged.connect(self._changed)
        for w in (self.pattern, self.despike, self.scale_blank):
            w.toggled.connect(self._changed)
        self.blank.currentIndexChanged.connect(self._changed)
        self._building = False
        self._update_enabled()

    # -------------------------------------------------------------- widgets
    @staticmethod
    def _spin(lo, hi, step, value, suffix, tip):
        w = QtWidgets.QSpinBox()
        w.setRange(lo, hi)
        w.setSingleStep(step)
        w.setValue(value)
        w.setSuffix(suffix)
        w.setToolTip(tip)
        w.setKeyboardTracking(False)
        return w

    @staticmethod
    def _dspin(lo, hi, step, decimals, suffix, value, tip):
        w = QtWidgets.QDoubleSpinBox()
        w.setRange(lo, hi)
        w.setSingleStep(step)
        w.setDecimals(decimals)
        w.setValue(value)
        w.setSuffix(suffix)
        w.setToolTip(tip)
        w.setKeyboardTracking(False)
        return w

    def _update_enabled(self) -> None:
        self.window.setEnabled(self.pattern.isChecked())
        self.spike_threshold.setEnabled(self.despike.isChecked())
        self.scale_blank.setEnabled(self.blank.currentIndex() > 0)
        self.stiffness_label.setText(f"{STIFFNESS_STEPS[self.stiffness.value()]:.0e}")

    def _changed(self, *_):
        self._update_enabled()
        if not self._building:
            self.settingsChanged.emit()

    # -------------------------------------------------------------- values
    def settings(self) -> BackgroundSettings:
        win = self.window.value()
        return BackgroundSettings(
            percentile=float(self.percentile.value()),
            stiffness=float(STIFFNESS_STEPS[self.stiffness.value()]),
            remove_pattern=self.pattern.isChecked(),
            pattern_window=win if win % 2 else win + 1,
            despike=self.despike.isChecked(),
            spike_threshold=float(self.spike_threshold.value()),
            noise_block=int(self.noise_block.value()),
            subtract_blank=self.blank.currentIndex() > 0,
            scale_blank=self.scale_blank.isChecked(),
        )

    def blank_name(self) -> str | None:
        return self.blank.currentData() if self.blank.currentIndex() > 0 else None

    def set_values(self, s: BackgroundSettings, blank: str | None) -> None:
        self._building = True
        self.percentile.setValue(s.percentile)
        nearest = int(np.argmin([abs(np.log10(v) - np.log10(s.stiffness)) for v in STIFFNESS_STEPS]))
        self.stiffness.setValue(nearest)
        self.pattern.setChecked(s.remove_pattern)
        self.window.setValue(s.pattern_window)
        self.despike.setChecked(s.despike)
        self.spike_threshold.setValue(s.spike_threshold)
        self.noise_block.setValue(s.noise_block)
        self.scale_blank.setChecked(s.scale_blank)
        i = self.blank.findData(blank) if (blank and s.subtract_blank) else 0
        self.blank.setCurrentIndex(max(i, 0))
        self._building = False
        self._update_enabled()

    def reset(self) -> None:
        self.set_values(BackgroundSettings(), None)
        self.settingsChanged.emit()

    def set_run_names(self, entries) -> None:
        current = self.blank_name()
        self._building = True
        self.blank.clear()
        self.blank.addItem("None", None)
        blanks = [e for e in entries if "blank" in (e.kind or "").lower()]
        others = [e for e in entries if e not in blanks]
        for e in blanks + others:
            label = e.name + (f"  ({e.kind})" if e.kind else "")
            self.blank.addItem(label, e.name)
        i = self.blank.findData(current) if current else 0
        self.blank.setCurrentIndex(max(i, 0))
        self._building = False
        self._update_enabled()

    # ---------------------------------------------------------- reporting
    def show_summary(self, text: str) -> None:
        self.summary.setText(text)

    def show_diagnostics(self, residual: np.ndarray | None, sigma: float | None,
                         top_ions: list[tuple[int, float]] | None = None) -> None:
        worst = ""
        if top_ions:
            worst = "Most signal removed at m/z " + ", ".join(
                f"{m} ({f:.0%})" for m, f in top_ions[:5]) + "."
        if residual is None or residual.size < 50 or not sigma:
            self.hist_bars.setData([], [])
            self.hist_fit.setData([], [])
            self.qq_points.setData([], [])
            msg = ("Apply or preview 🧹 to see diagnostics." if residual is None else
                   "This ion is mostly below the storage threshold, so there's no background to check.")
            self.diag_text.setText(msg + (f"\n{worst}" if worst else ""))
            return
        v = residual[np.isfinite(residual)]
        lo, hi = np.percentile(v, [0.5, 99.5])
        edges = np.linspace(lo, hi, 41)
        counts, _ = np.histogram(v, bins=edges)
        self.hist_bars.setData(edges, counts)
        x = np.linspace(lo, hi, 200)
        width = edges[1] - edges[0]
        self.hist_fit.setData(x, v.size * width * norm.pdf(x, 0.0, sigma))
        sample = np.sort(np.random.default_rng(0).choice(v, size=min(v.size, 1500), replace=False))
        q = norm.ppf((np.arange(sample.size) + 0.5) / sample.size)
        self.qq_points.setData(q, sample / sigma)
        mean, sd = float(v.mean()), float(v.std())
        verdict = ("Centred and matches the noise map." if abs(mean) < 0.5 * sigma and 0.7 < sd / sigma < 1.4
                   else "Background may be over- or under-corrected; try another percentile.")
        self.diag_text.setText(f"Residual mean {mean:+.2f}, SD {sd:.2f}; noise map SD {sigma:.2f}. "
                               f"{verdict}" + (f"\n{worst}" if worst else ""))
