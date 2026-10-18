"""Sequence-wide exceptions queue with candidate spectra and fitted-profile comparisons."""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from ..core.components import component_profiles_from_durable
from ..core.istd import peak_spectrum, review_candidates
from ..core.peaks import PeakSettings, make_frame
from .blob_state import MS_AWARE


class StandardReviewDialog(QtWidgets.QDialog):
    gotoRequested = QtCore.Signal(str, float, float)
    searchRequested = QtCore.Signal()

    def __init__(self, state, parent=None):
        super().__init__(parent)
        self.state = state
        self.setWindowTitle("Review standard exceptions")
        self.resize(1150, 850)
        self._pair = None
        self._generation = 0
        self._job = None
        self._candidates = []
        self._data = None
        self._theme_timer = QtCore.QTimer(self)
        self._theme_timer.setSingleShot(True)
        self._theme_timer.timeout.connect(self._draw_candidates)
        self.only_exceptions = QtWidgets.QCheckBox("Only missing, uncertain or stale matches")
        self.only_exceptions.setChecked(True)
        self.only_exceptions.toggled.connect(self.refresh)
        self.counter = QtWidgets.QLabel()
        self.next_btn = QtWidgets.QPushButton("Next exception")
        self.next_btn.clicked.connect(self.next_exception)
        search = QtWidgets.QPushButton("Search standards again")
        search.clicked.connect(self.searchRequested)
        top = QtWidgets.QHBoxLayout()
        for widget in (self.only_exceptions, self.counter, self.next_btn, search):
            top.addWidget(widget)
        self.queue = QtWidgets.QTableWidget(0, 3)
        self.queue.setHorizontalHeaderLabels(["Run", "Standard", "Why it needs review"])
        self.queue.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.queue.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.queue.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.queue.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.queue.horizontalHeader().setStretchLastSection(True)
        self.queue.itemSelectionChanged.connect(self._pair_selected)
        self.queue.setMaximumHeight(180)
        self.queue.setMinimumHeight(100)
        self.candidates = QtWidgets.QTableWidget(0, 6)
        self.candidates.setHorizontalHeaderLabels(["Candidate", "t1 (min)", "t2 (s)", "Spectrum score",
                                                   "Automatic thresholds", "Fit quality"])
        self.candidates.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.candidates.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.candidates.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.candidates.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.candidates.horizontalHeader().setStretchLastSection(True)
        self.candidates.itemSelectionChanged.connect(self._candidate_selected)
        self.candidates.setMaximumHeight(160)
        self.candidates.setMinimumHeight(110)
        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        self.spectra, self.profiles, self.notes = [], [], []
        comparisons = QtWidgets.QHBoxLayout()
        for title in ("Selected candidate", "Alternative candidate"):
            column = QtWidgets.QVBoxLayout()
            spectrum = pg.PlotWidget(title=title + ": spectrum")
            spectrum.setLabel("bottom", "m/z")
            spectrum.setLabel("left", "Relative intensity")
            profile = pg.PlotWidget(title="Local group: observed, fitted contributions and residual")
            profile.setLabel("bottom", "First retention time (min)")
            profile.setLabel("left", "Summed signal")
            spectrum.getAxis('bottom').enableAutoSIPrefix(False)
            profile.getAxis('bottom').enableAutoSIPrefix(False)
            note = QtWidgets.QLabel()
            note.setWordWrap(True)
            column.addWidget(spectrum, 1)
            column.addWidget(profile, 1)
            column.addWidget(note)
            comparisons.addLayout(column, 1)
            self.spectra.append(spectrum)
            self.profiles.append(profile)
            self.notes.append(note)
        self.override = QtWidgets.QCheckBox("Use this candidate despite failed automatic thresholds")
        self.override.toggled.connect(self._buttons)
        self.accept_btn = QtWidgets.QPushButton("Use selected candidate and next")
        self.accept_btn.clicked.connect(self.accept_candidate)
        self.ack_btn = QtWidgets.QPushButton("Keep current match; mark reviewed")
        self.ack_btn.clicked.connect(self.acknowledge)
        self.missing_btn = QtWidgets.QPushButton("Confirm missing and next")
        self.missing_btn.clicked.connect(self.confirm_missing)
        goto = QtWidgets.QPushButton("Show selected peak")
        goto.clicked.connect(self.goto_candidate)
        self.refit_btn = QtWidgets.QPushButton("Inspect/refit this group")
        self.refitRequested = None
        self.refit_btn.clicked.connect(self.inspect_group)
        self.refit_btn.setToolTip("Open the existing manual refit tool for this candidate's group.")
        self.decision_note = QtWidgets.QLabel("Review records your decision per run. The reference spectrum is unchanged. "
            "Unresolved-fit warnings and missing-standard flags remain in results and exports. "
            "Spectrum: grey = reference, colour = candidate. Profile: neutral = observed, orange = "
            "all fitted contributions, colour = candidate, dashed red = residual (including local background).")
        self.decision_note.setWordWrap(True)
        buttons = QtWidgets.QHBoxLayout()
        for widget in (self.accept_btn, self.ack_btn, self.missing_btn, goto, self.refit_btn):
            buttons.addWidget(widget)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addLayout(top)
        layout.addWidget(self.queue)
        layout.addWidget(self.status)
        layout.addWidget(self.candidates)
        layout.addLayout(comparisons, 1)
        layout.addWidget(self.override)
        layout.addLayout(buttons)
        layout.addWidget(self.decision_note)
        state.changed.connect(self.refresh)
        state.project.quant.runningChanged.connect(self._buttons)
        for signal in (state.project.blobs.changed, state.project.runsChanged,
                       state.project.layersChanged, state.project.foldChanged):
            signal.connect(self.refresh)
        self.refresh()

    def refresh(self, *_):
        previous = self._pair
        previous_row = max(self.queue.currentRow(), 0)
        exceptions = self.state.exceptions()
        pairs = exceptions if self.only_exceptions.isChecked() else [
            (self.state.project.entries[i].name, name) for i in self.state.project.loaded()
            for name in self.state.references]
        self.queue.blockSignals(True)
        self.queue.setRowCount(len(pairs))
        for row, (run, name) in enumerate(pairs):
            match = self.state.matches.get(run, {}).get(name)
            why = ("Out of date: rerun standard search" if self.state.is_stale(run)
                   else "Not searched" if match is None
                   else "Reviewed (fit warnings retained)" if self.state.reviewed(run, name)
                   else match.reason or ("Uncertain fit" if match.needs_review else "Found"))
            for column, value in enumerate((run, name, why)):
                self.queue.setItem(row, column, QtWidgets.QTableWidgetItem(value))
        self._pairs = pairs
        self.queue.blockSignals(False)
        self.counter.setText(f"{len(exceptions)} decision(s) remaining")
        self.next_btn.setEnabled(bool(exceptions))
        if pairs:
            row = pairs.index(previous) if previous in pairs else min(previous_row, len(pairs)-1)
            self.queue.selectRow(row)
            self._pair_selected()
        else:
            self._clear()
            self.status.setText("Label at least one internal standard to create the review queue."
                if not self.state.references else
                "Load chromatograms before reviewing standards." if not self.state.project.loaded() else
                "No exceptions remain. Reviewed uncertain fits still retain their quality warnings.")

    def _clear(self):
        self._generation += 1
        if self._job is not None:
            self.state.project.jobs.cancel(self._job)
        self._job = None
        self._pair, self._data, self._candidates = None, None, []
        self.candidates.setRowCount(0)
        for spectrum, profile, note in zip(self.spectra, self.profiles, self.notes):
            spectrum.clear()
            profile.clear()
            note.clear()
        self._buttons()

    def _pair_selected(self):
        if not self.isVisible():
            return
        row = self.queue.currentRow()
        if not 0 <= row < len(self._pairs):
            return
        pair = self._pairs[row]
        p = self.state.project
        # A refresh can invalidate a displayed comparison even when its row stays selected.
        if pair == self._pair and self._data is not None:
            i = p.index_of(pair[0])
            if i is not None and self._data['signature'] == self.state.input_signature(i) and not self.state.is_stale(pair[0]):
                self._buttons()
                return
        self._clear()
        self._pair = pair
        run, name = pair
        self.override.setChecked(False)
        if self.state.is_stale(run):
            self.status.setText("These matches are out of date or not searched. Search standards again "
                                "after refreshing Peaks; decisions are disabled until results are current.")
            self._buttons()
            return
        i = p.index_of(run)
        signature = self.state.input_signature(i)
        generation = self._generation
        reference = self.state.references[name]
        match = self.state.matches.get(run, {}).get(name)
        expected = ((match.expected_t1, match.expected_t2) if match and np.isfinite(match.expected_t1)
                    else (reference.t1, reference.t2))
        self.status.setText(f"Loading candidate spectra and profiles for {name} on {run}…")
        base, f = p.base_run(i), p.fold(i)
        components = p.blobs.components_for(i) if p.blobs.method == MS_AWARE else None
        settings = self.state.settings

        def work(token, progress):
            peaks = self.state._peaks_for(i)
            token.check()
            candidates = review_candidates(reference, run, peaks, settings, expected)
            details = []
            frame = make_frame(f, PeakSettings(wrap_below_s=p.blobs.component_settings.candidates.wrap_below_s))
            observed = frame.extend(f.tic_image(base)) if components is not None else None
            for candidate in candidates:
                token.check()
                index = self.state._resolve_index(peaks, candidate.peak_id)
                detail = dict(match=candidate, mz=peaks.mz, spectrum=peak_spectrum(peaks, index), profile=None,
                              group_index=None)
                if components is not None:
                    group = next((g for g in components.groups if any(c.id == candidate.peak_id for c in g.components)), None)
                    if group is not None:
                        profiles = [component_profiles_from_durable(c, (group.n_rows, group.n_cols))*c.spectrum.sum()
                                    for c in group.components]
                        selected = next(j for j, c in enumerate(group.components) if c.id == candidate.peak_id)
                        local = observed[group.row0:group.row0+group.n_rows, group.col0:group.col0+group.n_cols]
                        x = (group.t1_start + (group.col0+np.arange(group.n_cols)+.5)*group.period)/60
                        detail['profile'] = (x, local.sum(axis=0), sum(profiles).sum(axis=0), profiles[selected].sum(axis=0))
                        detail['warning'] = group.warning
                        detail['group_index'] = group.region_index
                details.append(detail)
            return dict(reference=reference, candidates=details, signature=signature)

        def done(result):
            if generation != self._generation or self._pair != pair or self.state.is_stale(run):
                return
            self._job = None
            self._data = result
            self._candidates = result['candidates']
            self.candidates.blockSignals(True)
            self.candidates.setRowCount(len(self._candidates))
            for row, detail in enumerate(self._candidates):
                c = detail['match']
                values = [str(c.peak_id), f"{c.t1/60:.3f}", f"{c.t2:.3f}", f"{c.spectral_score:.3f}",
                          "pass" if c.found else c.reason,
                          "unresolved" if c.peak_unresolved else c.peak_confidence or "legacy"]
                for col, value in enumerate(values):
                    self.candidates.setItem(row, col, QtWidgets.QTableWidgetItem(value))
            self.candidates.blockSignals(False)
            self.status.setText(f"{name} on {run}: compare alternatives, then keep a match or confirm it is missing.")
            if self._candidates:
                self.candidates.selectRow(0)
            self._draw_candidates()

        def failed(message):
            if generation == self._generation:
                self.status.setText(f"Could not load comparison: {message}")
        self._job = p.jobs.submit(f"Review standard candidates: {name} on {run}", work, done, failed,
                                 on_cancel=lambda: failed("cancelled"))

    def _buttons(self, *_):
        current = bool(self._pair and not self.state.running and not self.state.project.quant.running
                       and not self.state.is_stale(self._pair[0]))
        row = self.candidates.currentRow()
        selected = self._candidates[row] if 0 <= row < len(self._candidates) else None
        self.accept_btn.setEnabled(current and selected is not None
                                   and (selected['match'].found or self.override.isChecked()))
        self.ack_btn.setEnabled(current and self.state.matches.get(self._pair[0], {}).get(self._pair[1]) is not None
                                if self._pair else False)
        self.missing_btn.setEnabled(current)
        self.refit_btn.setEnabled(current and selected is not None and selected.get('group_index') is not None)

    def _draw_candidates(self):
        if not self.isVisible():
            return
        row = self.candidates.currentRow()
        alternative = next((j for j in range(len(self._candidates)) if j != row), -1)
        palette = QtWidgets.QApplication.palette()
        foreground = palette.color(QtGui.QPalette.ColorRole.WindowText)
        for j, (spectrum, profile, note, index) in enumerate(zip(self.spectra, self.profiles, self.notes, (row, alternative))):
            spectrum.clear()
            profile.clear()
            spectrum.setBackground(palette.color(QtGui.QPalette.ColorRole.Base))
            profile.setBackground(palette.color(QtGui.QPalette.ColorRole.Base))
            spectrum.setTitle(('Selected' if j == 0 else 'Alternative') + ' candidate: spectrum', color=foreground.name())
            profile.setTitle('Local group: observed, fitted contributions and residual', color=foreground.name())
            for plot in (spectrum, profile):
                for axis in ('left', 'bottom'):
                    plot.getAxis(axis).setTextPen(foreground)
                    plot.getAxis(axis).setPen(foreground)
            if not 0 <= index < len(self._candidates):
                note.setText("No alternative candidate")
                continue
            d = self._candidates[index]
            c, ref = d['match'], self._data['reference']
            colour = '#169ab4' if j == 0 else '#d34b99'
            spectrum.plot(np.repeat(ref.mz0+np.arange(len(ref.spectrum)), 2),
                          np.column_stack((np.zeros(len(ref.spectrum)), ref.spectrum/max(ref.spectrum.max(), 1e-12))).ravel(),
                          pen=pg.mkPen('#888888', width=3), connect='pairs')
            spectrum.plot(np.repeat(d['mz'], 2),
                          np.column_stack((np.zeros(len(d['spectrum'])), d['spectrum']/max(d['spectrum'].max(), 1e-12))).ravel(),
                          pen=pg.mkPen(colour, width=1.5), connect='pairs')
            if d['profile'] is not None:
                x, observed, fitted, contribution = d['profile']
                profile.plot(x, observed, pen=pg.mkPen(foreground, width=1.5))
                profile.plot(x, fitted, pen=pg.mkPen('#d39426', width=1.5))
                profile.plot(x, contribution, pen=pg.mkPen(colour, width=2))
                profile.plot(x, observed-fitted, pen=pg.mkPen('#d85a5a', style=QtCore.Qt.PenStyle.DashLine))
            note.setText(f"{c.peak_id}: score {c.spectral_score:.3f}; t1 {c.t1/60:.3f} min, t2 {c.t2:.3f} s. "
                         + (d.get('warning') or c.reason or "No automatic fit warning"))
        self._buttons()

    def _candidate_selected(self):
        self.override.setChecked(False)
        self._draw_candidates()

    def _decide(self, **kwargs):
        if self._pair is None:
            return
        try:
            self.state.review_match(*self._pair, **kwargs)
        except ValueError as exc:
            self.status.setText(str(exc))

    def accept_candidate(self):
        row = self.candidates.currentRow()
        if self.accept_btn.isEnabled() and 0 <= row < len(self._candidates):
            self._decide(match=self._candidates[row]['match'])

    def acknowledge(self):
        self._decide()

    def confirm_missing(self):
        self._decide(missing=True)

    def next_exception(self):
        pending = set(self.state.exceptions())
        for offset in range(1, len(self._pairs)+1):
            row = (self.queue.currentRow()+offset) % len(self._pairs)
            if self._pairs[row] in pending:
                self.queue.selectRow(row)
                return

    def goto_candidate(self):
        row = self.candidates.currentRow()
        if self._pair and 0 <= row < len(self._candidates):
            c = self._candidates[row]['match']
            self.gotoRequested.emit(self._pair[0], c.t1, c.t2)

    def inspect_group(self):
        row = self.candidates.currentRow()
        if self.refitRequested and self._pair and 0 <= row < len(self._candidates):
            self.refitRequested(self._pair[0], self._candidates[row]['group_index'])

    def closeEvent(self, event):
        self._theme_timer.stop()
        self._clear()
        super().closeEvent(event)

    def showEvent(self, event):
        super().showEvent(event)
        self.refresh()

    def changeEvent(self, event):
        super().changeEvent(event)
        if (event.type() in (QtCore.QEvent.Type.PaletteChange, QtCore.QEvent.Type.ApplicationPaletteChange)
                and hasattr(self, '_theme_timer') and self.isVisible()):
            # Clearing/rebuilding graphics during a global Qt repolish can leave
            # its native font/palette traversal holding deleted scene objects.
            self._theme_timer.start(0)
