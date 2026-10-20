"""Cancellable peak-fitting previews isolated from accepted project results."""

from __future__ import annotations

from dataclasses import replace

from PySide6 import QtCore, QtWidgets
from .theme import fit_to_screen

from ..core.background import storage_threshold
from ..core.components import detect_components
from ..core.fold import fold
from ..core.io import crop_run
from ..core.peaktable import model_on_frame, respectrum_components, signal_accounting
from .dialogs import SettingsDialog
from .fit_map_dialog import FitMapWindow


def preview_fit(run, raw, original_fold, settings, name, bounds, token, progress):
    """Use the production detector on a run or time interval, preserving fold origin.

    Every second-dimension row is retained, so a rectangular selection does not
    cut through a peak across modulations. The display can zoom to its t2 bounds.
    """
    token.check()
    if bounds is not None:
        start, end, _low, _high = bounds
        if end <= start:
            raise ValueError("Preview interval must have positive width")
        run = crop_run(run, start-original_fold.spec.period, end+original_fold.spec.period)
        if run.n_scans < 2:
            raise ValueError("Choose a preview interval containing at least two scans.")
        f = fold(run, original_fold.spec.period, original_fold.spec.phase, t0=original_fold.spec.t0)
    else:
        f = original_fold
    phase = ["Finding candidate peaks"]
    def stage(message):
        token.check()
        phase[0] = message
        progress(0, 0, message)
    candidates, regions, components = detect_components(
        run, f, name, replace(settings, workers=1), floor=storage_threshold(raw)/3,
        progress=lambda done, total: progress(done, total, phase[0]), cancel=token.check, stage=stage)
    stage("Solving full spectra and building residual maps")
    _peaks, model, footprint = respectrum_components(
        run, f, components, wrap_below_s=settings.candidates.wrap_below_s,
        progress=lambda done, total: progress(done, total, phase[0]), cancel=token.check)
    token.check()
    tic = f.tic_image(run)
    accounting = signal_accounting(tic, model, footprint, f.row_pitch, f.spec.period, float(f.t1[0]))
    model, _mask = model_on_frame(model, footprint, tic.shape[0])
    outlines = [group.contour(j, .10) for group in components.groups for j in range(len(group.components))]
    return dict(tic=tic, model=model, report=accounting, outlines=outlines,
                start=float(f.t1[0]), period=f.spec.period, pitch=f.row_pitch,
                count=len(components.components), unresolved=sum(g.unresolved for g in components.groups),
                bounds=bounds)


class PeakPreviewDialog(QtWidgets.QDialog):
    fitAllRequested = QtCore.Signal()

    def __init__(self, project, index, bounds, parent=None):
        super().__init__(parent)
        self.project, self.entry, self.bounds = project, project.entries[index], bounds
        self.settings = project.blobs.component_settings
        self._job = None
        self._generation = 0
        self._basis = None
        self.setWindowTitle(f"Preview peak fitting: {self.entry.name}")
        fit_to_screen(self, 1150, 850)
        self.scope = QtWidgets.QComboBox()
        self.scope.addItem("Visible first-time interval (all second-time rows)", "region")
        self.scope.addItem("Current run", "run")
        self.scope.setToolTip("Retains full second-dimension profiles to avoid clipping overlapping peaks.")
        self.scope.currentIndexChanged.connect(self._invalidate)
        self.settings_btn = QtWidgets.QPushButton("Preview settings…")
        self.settings_btn.clicked.connect(self.edit_settings)
        self.run_btn = QtWidgets.QPushButton("Run preview")
        self.run_btn.clicked.connect(self.run_preview)
        self.cancel_btn = QtWidgets.QPushButton("Cancel preview")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self.cancel)
        self.apply_btn = QtWidgets.QPushButton("Use settings and fit all runs")
        self.apply_btn.setEnabled(False)
        self.apply_btn.clicked.connect(self.apply_settings)
        self.progress = QtWidgets.QProgressBar()
        self.status = QtWidgets.QLabel("Preview is separate from accepted peaks and quantitative results. "
                                     "An interval uses the existing background layer and may differ at its edges.")
        self.status.setWordWrap(True)
        self.map = FitMapWindow(self)
        self.map.setWindowFlags(QtCore.Qt.WindowType.Widget)
        for plot in self.map.plots:
            plot.getAxis('bottom').enableAutoSIPrefix(False)
        row = QtWidgets.QHBoxLayout()
        for widget in (self.scope, self.settings_btn, self.run_btn, self.cancel_btn, self.apply_btn):
            row.addWidget(widget)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(row)
        layout.addWidget(self.status)
        layout.addWidget(self.progress)
        layout.addWidget(self.map, 1)
        for signal in (project.foldChanged, project.layersChanged, project.runsChanged, project.blobs.changed):
            signal.connect(self._inputs_changed)
        project.jobs.jobChanged.connect(self._job_progress)

    def _job_progress(self, job):
        if job.id != self._job:
            return
        self.status.setText(job.message)
        self.progress.setRange(0, job.total if job.total else 0)
        if job.total:
            self.progress.setValue(job.done)

    def _signature(self):
        if self.entry not in self.project.entries:
            return None
        i = self.project.entries.index(self.entry)
        return (id(self.project.base_run(i)), self.project.period, self.entry.phase,
                self.project.blobs.component_settings.to_dict(), self.scope.currentData())

    def _invalidate(self, *_):
        self._generation += 1
        if self._job is not None:
            self.project.jobs.cancel(self._job)
        self._job = None
        self.apply_btn.setEnabled(False)
        self.run_btn.setEnabled(self.entry in self.project.entries and self.entry.run is not None)
        self.cancel_btn.setEnabled(False)
        self.progress.setRange(0, 1)
        self.progress.setValue(0)
        self.status.setText("Preview settings, scope or inputs changed. Run a fresh preview before applying settings.")

    def _inputs_changed(self, *_):
        if self._basis is not None and self._basis != self._signature():
            self._invalidate()
            self.status.setText("Inputs changed. Run a fresh preview before applying settings.")

    def edit_settings(self):
        dialog = SettingsDialog("Preview peak-fitting settings", self.settings, parent=self)
        if dialog.exec():
            self.settings = dialog.settings()
            self._invalidate()
            self.status.setText("Preview settings changed; accepted project settings are unchanged.")

    def run_preview(self):
        self._invalidate()
        if self.entry not in self.project.entries or self.entry.run is None:
            self.status.setText("This run is no longer loaded.")
            return
        index = self.project.entries.index(self.entry)
        run, raw, f = self.project.base_run(index), self.entry.run, self.project.fold(index)
        source = ('background preview' if self.project.preview is not None and self.project.preview[0] is self.entry
                  else 'existing background layer' if self.project.show_processed and self.entry.background is not None
                  else 'raw signal; no background layer in use')
        settings = self.settings
        bounds = self.bounds if self.scope.currentData() == "region" else None
        self._basis = self._signature()
        generation = self._generation
        self.run_btn.setEnabled(False)
        self.apply_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.progress.setRange(0, 0)
        self.status.setText("Finding candidate peaks… Progress is also available in Jobs.")

        def work(token, progress):
            def report(done, total, message):
                progress(done, total, message)
            return preview_fit(run, raw, f, settings, self.entry.name, bounds, token, report)

        def done(result):
            if generation != self._generation or self._basis != self._signature() or not self.isVisible():
                return
            self._job = None
            self.run_btn.setEnabled(True)
            self.cancel_btn.setEnabled(False)
            self.apply_btn.setEnabled(True)
            self.progress.setRange(0, 1)
            self.progress.setValue(1)
            self.status.setText(f"Preview from {source}: {result['count']} components; {result['unresolved']} unresolved groups. "
                                "Inspect the fitted sum and residual before applying settings to all runs.")
            self.map.set_data(self.entry.name + " (preview)", result['tic'], result['model'],
                              result['start'], result['period'], result['pitch'], result['report'], result['outlines'])
            self.map.plots[0].setTitle(f"Observed ({source})")
            self.map.plots[3].setTitle("Fitted fraction (component contribution ÷ observed)")
            if bounds:
                x0, x1, y0, y1 = bounds
                self.map.plots[0].setRange(xRange=(x0/60, x1/60), yRange=(y0, y1), padding=0)

        def ended(message):
            if generation == self._generation:
                self._invalidate()
                self.status.setText(message)
        self._job = self.project.jobs.submit(f"Peak preview: {self.entry.name}", work, done,
                    lambda message: ended(f"Preview failed: {message}"),
                    on_cancel=lambda: ended("Preview cancelled. Accepted results are unchanged."))

    def cancel(self):
        self._invalidate()
        self.status.setText("Preview cancelled. Accepted results are unchanged.")

    def apply_settings(self):
        if not self.apply_btn.isEnabled() or self._basis != self._signature():
            self._inputs_changed()
            return
        self.project.blobs.set_component_settings(self.settings)
        self.fitAllRequested.emit()
        self.close()

    def closeEvent(self, event):
        self._invalidate()
        super().closeEvent(event)

    def showEvent(self, event):
        super().showEvent(event)
        self.map.canvas.setBackground(QtWidgets.QApplication.palette().color(
            self.backgroundRole()))
