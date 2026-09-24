"""Inspect and refit one legacy blob or one MS-aware component group.

The joint fit displays overlapping contribution maps rather than exclusive
pixel ownership. Fits made from an MS-aware analysis region can be accepted
into project state; acceptance is persistent and undoable.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets
from skimage.measure import find_contours

from ..core.background import storage_threshold
from ..core.candidates import candidate_cube, candidate_regions, find_candidates
from ..core.components import group_from_fit
from ..core.overlap import (MODELS, OverlapSettings, assess_components, blob_cube,
                            choose_k, component_maps, fit_overlap, fit_template_overlap, mcr_als,
                            split_by_watershed, split_volumes, unique_ion_apexes)
from .views.surface3d import Surface3DView


class ComponentMapView(QtWidgets.QWidget):
    """Observed TIC with overlapping fitted-component contributions."""

    COLOURS = ((0, 220, 255), (255, 105, 180), (255, 205, 60),
               (120, 235, 120), (190, 130, 255), (255, 140, 70))

    def __init__(self, parent=None):
        super().__init__(parent)
        self.plot = pg.PlotWidget()
        self.plot.setLabel("bottom", "t₁ (min)")
        self.plot.setLabel("left", "t₂ (s)")
        self.plot.setAspectLocked(False)
        self.image = pg.ImageItem(axisOrder="row-major")
        self.image.setColorMap(pg.ColorMap([0.0, 1.0], [(0, 0, 0), (220, 220, 220)]))
        self.plot.addItem(self.image)
        self.note = QtWidgets.QLabel(
            "Colours show fitted contributions, not ownership. Solid contours are 10% "
            "of that component's maximum; dotted contours are 50%. Overlap is expected.")
        self.note.setWordWrap(True)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(self.plot, 1)
        lay.addWidget(self.note)
        self._items = []

    def clear_fit(self) -> None:
        for item in self._items:
            self.plot.removeItem(item)
        self._items.clear()
        self.image.clear()

    def set_fit(self, cube, fit) -> None:
        self.clear_fit()
        if cube is None or fit is None:
            return
        observed = np.log1p(np.maximum(cube.data.sum(axis=2), 0.0))
        self.image.setImage(observed, autoLevels=True)
        rect = QtCore.QRectF(cube.t1_of(-0.5) / 60.0, cube.t2_of(-0.5),
                             cube.data.shape[1] * cube.period / 60.0,
                             cube.data.shape[0] * cube.pitch)
        self.image.setRect(rect)
        maps = component_maps(fit, cube)
        assessments = assess_components(fit)
        for j, (component, assessment) in enumerate(zip(fit.components, assessments)):
            contribution = maps[:, :, j]
            maximum = float(contribution.max())
            if maximum <= 0:
                continue
            colour = self.COLOURS[j % len(self.COLOURS)]
            rgba = np.zeros(contribution.shape + (4,), dtype=np.uint8)
            rgba[..., :3] = colour
            rgba[..., 3] = np.asarray(
                115 * np.sqrt(np.clip(contribution / maximum, 0, 1)), dtype=np.uint8)
            layer = pg.ImageItem(rgba, axisOrder="row-major")
            layer.setRect(rect)
            layer.setZValue(5 + j)
            self.plot.addItem(layer)
            self._items.append(layer)
            for fraction, style in ((0.10, QtCore.Qt.PenStyle.SolidLine),
                                    (0.50, QtCore.Qt.PenStyle.DotLine)):
                for contour in find_contours(contribution, fraction * maximum):
                    x = cube.t1_of(contour[:, 1]) / 60.0
                    y = cube.t2_of(contour[:, 0])
                    curve = pg.PlotDataItem(x, y, pen=pg.mkPen(colour, width=2, style=style))
                    curve.setZValue(20 + j)
                    self.plot.addItem(curve)
                    self._items.append(curve)
            label = pg.TextItem(
                f"C{j+1}  {component.volume:,.0f}\n{assessment.confidence}",
                anchor=(0.5, 1.1), color=colour)
            label.setPos(component.t1 / 60.0, component.t2)
            label.setZValue(30 + j)
            self.plot.addItem(label)
            self._items.append(label)
        self.plot.autoRange()


class OverlapDialog(QtWidgets.QDialog):
    fitAccepted = QtCore.Signal(str, object)

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.setWindowTitle("✂️ Overlap")
        self.resize(880, 640)
        self.project = project
        self._name = None
        self._run = None
        self._fold = None
        self._blobs = None
        self._series = None
        self._candidates = None
        self._regions = []
        self._source = "blobs"
        self.cube = None
        self._last_fit = None

        self.blob_combo = QtWidgets.QComboBox()
        self.blob_combo.currentIndexChanged.connect(self._blob_chosen)

        top = QtWidgets.QHBoxLayout()
        self.source_label = QtWidgets.QLabel("Legacy blob:")
        top.addWidget(self.source_label)
        top.addWidget(self.blob_combo, 1)
        self.find_regions_btn = QtWidgets.QPushButton("Find MS-aware candidate regions")
        self.find_regions_btn.setToolTip(
            "Find peak centres ion by ion and build wider analysis rectangles without hard ownership borders")
        self.find_regions_btn.clicked.connect(self._find_candidate_regions)
        top.addWidget(self.find_regions_btn)

        self.tabs = QtWidgets.QTabWidget()
        self.surface = Surface3DView()
        self.tabs.addTab(self.surface, "3D surface")
        self.component_view = ComponentMapView()
        self.tabs.addTab(self.component_view, "Component contributions")
        self.tabs.addTab(self._build_apex_page(), "Unique ions / split")
        self.tabs.addTab(self._build_fit_page(), "Peak-shape fit")
        self.tabs.addTab(self._build_mcr_page(), "MCR-ALS")

        lay = QtWidgets.QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(self.tabs, 1)

    def set_background(self, colour) -> None:
        self.surface.set_background(colour)

    # ------------------------------------------------------------------ pages
    def _build_apex_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(page)
        self.apex_table = QtWidgets.QTableWidget(0, 4)
        self.apex_table.setHorizontalHeaderLabels(["ion", "t1 (min)", "t2 (s)", "cluster"])
        self.apex_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        lay.addWidget(self.apex_table)
        self.apex_summary = QtWidgets.QLabel("")
        self.apex_summary.setWordWrap(True)
        lay.addWidget(self.apex_summary)
        split_btn = QtWidgets.QPushButton("Split at each cluster's own apex (watershed)")
        split_btn.clicked.connect(self._run_split)
        lay.addWidget(split_btn)
        self.split_table = QtWidgets.QTableWidget(0, 3)
        self.split_table.setHorizontalHeaderLabels(["region", "t1 (min) / t2 (s)", "volume"])
        self.split_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        lay.addWidget(self.split_table)
        return page

    def _build_fit_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(page)
        controls = QtWidgets.QHBoxLayout()
        self.model_combo = QtWidgets.QComboBox()
        self.model_combo.addItems(list(MODELS) + ["measured alkane template"])
        self.model_combo.setToolTip("gaussian: axis-aligned. emg: tail along t1. emg_t2: tail along t2. "
                                    "emg2: both. tilt: t2 apex drifts with t1. cell: t1 profile integrated "
                                    "over each modulation (peaks narrower than a column). See "
                                    "docs/validation/hydrocarbons-2026-09-20.md: on real alkane peaks the gain "
                                    "over gaussian was at most a few per cent.")
        self.auto_k = QtWidgets.QCheckBox("Choose k automatically (BIC)")
        self.auto_k.setChecked(True)
        self.k_spin = QtWidgets.QSpinBox()
        self.k_spin.setRange(1, 6)
        self.k_spin.setValue(2)
        self.k_spin.setEnabled(False)
        self.auto_k.toggled.connect(lambda on: self.k_spin.setEnabled(not on))
        fit_btn = QtWidgets.QPushButton("Fit peak shapes")
        fit_btn.clicked.connect(self._run_fit)
        self.accept_fit_btn = QtWidgets.QPushButton("Accept components")
        self.accept_fit_btn.setEnabled(False)
        self.accept_fit_btn.setToolTip(
            "Save these overlapping component shapes, spectra and estimated amounts in the project")
        self.accept_fit_btn.clicked.connect(self._accept_fit)
        for w in (QtWidgets.QLabel("Model:"), self.model_combo, self.auto_k,
                 QtWidgets.QLabel("k:"), self.k_spin, fit_btn, self.accept_fit_btn):
            controls.addWidget(w)
        controls.addStretch(1)
        lay.addLayout(controls)
        self.fit_table = QtWidgets.QTableWidget(0, 7)
        self.fit_table.setHorizontalHeaderLabels(["component", "t1 (min)", "t2 (s)",
                                                  "position SE (px, row/col)", "estimated amount",
                                                  "relative size", "confidence"])
        self.fit_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        lay.addWidget(self.fit_table)
        self.fit_status = QtWidgets.QLabel("")
        self.fit_status.setWordWrap(True)
        lay.addWidget(self.fit_status)
        return page

    def _build_mcr_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(page)
        controls = QtWidgets.QHBoxLayout()
        self.mcr_k = QtWidgets.QSpinBox()
        self.mcr_k.setRange(1, 6)
        self.mcr_k.setValue(2)
        self.mcr_unimodal = QtWidgets.QCheckBox("Unimodal")
        self.mcr_unimodal.setChecked(True)
        mcr_btn = QtWidgets.QPushButton("Run MCR-ALS")
        mcr_btn.clicked.connect(self._run_mcr)
        for w in (QtWidgets.QLabel("k:"), self.mcr_k, self.mcr_unimodal, mcr_btn):
            controls.addWidget(w)
        controls.addStretch(1)
        lay.addLayout(controls)
        self.mcr_table = QtWidgets.QTableWidget(0, 3)
        self.mcr_table.setHorizontalHeaderLabels(["component", "dominant ion", "volume"])
        self.mcr_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        lay.addWidget(self.mcr_table)
        self.mcr_status = QtWidgets.QLabel("")
        self.mcr_status.setWordWrap(True)
        lay.addWidget(self.mcr_status)
        return page

    # ------------------------------------------------------------------ context
    def set_run(self, name, run, fold, blobs, series=None, *, candidates=None,
                regions=None, candidate_region=None) -> None:
        """Called by the main window with whichever run is open. Lists
        co-eluting blobs (🔍's own flag) first, since they're the ones
        this dialog exists for, but every blob is available."""
        self._name, self._run, self._fold, self._blobs, self._series = name, run, fold, blobs, series
        self._candidates, self._regions, self._source = candidates, regions or [], "blobs"
        self.source_label.setText("Legacy blob:")
        self.find_regions_btn.setEnabled(run is not None and fold is not None)
        self.find_regions_btn.setText("Find MS-aware candidate regions")
        self.setWindowTitle(f"✂️ Overlap -- {name}" if name else "✂️ Overlap")
        self.blob_combo.blockSignals(True)
        self.blob_combo.clear()
        if candidates is not None and regions is not None:
            self._source = "candidates"
            self.source_label.setText("MS-aware region:")
            self._populate_candidate_regions()
            if candidate_region is not None and 0 <= int(candidate_region) < self.blob_combo.count():
                self.blob_combo.setCurrentIndex(int(candidate_region))
        elif blobs is not None and len(blobs):
            co = blobs.co_elution if blobs.co_elution.shape[0] == len(blobs) else None
            order = sorted(range(len(blobs)), key=lambda i: 0 if (co is not None and co[i]) else 1)
            for i in order:
                flag = "  (co-elution?)" if (co is not None and co[i]) else ""
                self.blob_combo.addItem(f"#{i}  t1={blobs.t1[i] / 60:.3f} min  t2={blobs.t2[i]:.3f} s{flag}", i)
        self.blob_combo.blockSignals(False)
        self._blob_chosen()

    def _populate_candidate_regions(self) -> None:
        self.blob_combo.clear()
        for region_index, region in enumerate(self._regions):
            idx = np.asarray(region.candidate_indices, dtype=int)
            strongest = int(idx[np.argmax(self._candidates.score[idx])])
            self.blob_combo.addItem(
                f"#{region_index}  t1={self._candidates.t1[strongest]/60:.3f} min  "
                f"t2={self._candidates.t2[strongest]:.3f} s  "
                f"({len(idx)} candidate{'s' if len(idx) != 1 else ''})", region_index)

    def _find_candidate_regions(self) -> None:
        if self._run is None or self._fold is None:
            return
        run, folded, name = self._run, self._fold, self._name
        index = self.project.index_of(name) if name is not None else None
        raw = self.project.entries[index].run if index is not None else run
        floor = storage_threshold(raw) / 3.0
        self.find_regions_btn.setEnabled(False)
        self.find_regions_btn.setText("Finding candidate regions…")

        def work(token, progress):
            candidates = find_candidates(run, folded, floor=floor)
            regions = candidate_regions(
                candidates, (folded.n_rows + candidates.wrap_rows, folded.n_cols),
                pad_t2_s=0.4, pad_modulations=4)
            return candidates, regions

        def done(result):
            if name != self._name:
                return
            self._candidates, self._regions = result
            self._source = "candidates"
            self.source_label.setText("MS-aware region:")
            self.blob_combo.blockSignals(True)
            self._populate_candidate_regions()
            self.blob_combo.blockSignals(False)
            self.find_regions_btn.setEnabled(True)
            self.find_regions_btn.setText(f"Refresh candidate regions ({len(self._regions)})")
            self._blob_chosen()

        def failed(text):
            self.find_regions_btn.setEnabled(True)
            self.find_regions_btn.setText("Find MS-aware candidate regions")
            self.fit_status.setText(f"Candidate detection failed: {text}")

        self.project.jobs.submit("✂️ Finding MS-aware candidate regions", work, done, failed)

    def _blob_chosen(self) -> None:
        self.cube = None
        self._last_fit = None
        self.accept_fit_btn.setEnabled(False)
        self.component_view.clear_fit()
        for table in (self.apex_table, self.split_table, self.fit_table, self.mcr_table):
            table.setRowCount(0)
        self.apex_summary.setText("")
        self.fit_status.setText("")
        self.mcr_status.setText("")
        i = self.blob_combo.currentData()
        if i is None or self._run is None:
            self.surface.set_data(np.zeros((0, 0)))
            return
        try:
            if self._source == "candidates":
                self.cube = candidate_cube(self._run, self._fold, self._candidates,
                                           self._regions[int(i)], max_ions=60)
            elif self._blobs is not None:
                self.cube = blob_cube(self._run, self._fold, self._blobs, int(i),
                                      OverlapSettings())
        except ValueError:
            self.cube = None
            return
        self.surface.set_data(self.cube.data.sum(axis=2))
        self._fill_apex_map()

    # ------------------------------------------------------------------ 1/2
    def _fill_apex_map(self) -> None:
        result = unique_ion_apexes(self.cube)
        self.apex_table.setRowCount(len(result.apexes))
        for r, a in enumerate(result.apexes):
            values = [f"{a.ion:g}", f"{a.t1 / 60:.3f}", f"{a.t2:.3f}", str(a.cluster)]
            for c, v in enumerate(values):
                self.apex_table.setItem(r, c, QtWidgets.QTableWidgetItem(v))
        if result.n_components <= 1:
            self.apex_summary.setText("Every ion peaks in the same place -- probably one compound.")
        else:
            self.apex_summary.setText(f"{result.n_components} apex clusters: these ions don't agree on "
                                      "where the compound peaks -- probably more than one.")
        self._apex_result = result

    def _run_split(self) -> None:
        self.split_table.setRowCount(0)
        if self.cube is None or getattr(self, "_apex_result", None) is None:
            return
        result = self._apex_result
        if result.n_components <= 1:
            self.apex_summary.setText(self.apex_summary.text() + "  Nothing to split.")
            return
        seeds = []
        for cluster in sorted(set(a.cluster for a in result.apexes)):
            members = [a for a in result.apexes if a.cluster == cluster]
            best = max(members, key=lambda a: self.cube.data[int(round(a.row)), int(round(a.col)),
                                                              int(np.flatnonzero(self.cube.ions == a.ion)[0])])
            seeds.append((best.row, best.col))
        labels = split_by_watershed(self.cube, seeds)
        totals = split_volumes(self.cube, labels)
        self.split_table.setRowCount(len(totals))
        for r, (lab, info) in enumerate(sorted(totals.items())):
            self.split_table.setItem(r, 0, QtWidgets.QTableWidgetItem(str(lab)))
            self.split_table.setItem(r, 1, QtWidgets.QTableWidgetItem(
                f"{info['t1'] / 60:.3f} min / {info['t2']:.3f} s"))
            self.split_table.setItem(r, 2, QtWidgets.QTableWidgetItem(f"{info['volume']:,.0f}"))

    # ------------------------------------------------------------------ 3
    def _measured_alkane_template(self):
        """Build a local real-peak profile from a confirmed alkane ladder."""
        if self._series is None or not self._series.hits:
            raise ValueError("Find the alkane series on this run before using a measured template")
        choices = [carbon for carbon in self._series.hits if carbon != 17]
        if not choices:
            raise ValueError("The alkane series has no non-C17 reference peak")
        carbon = 16 if 16 in choices else min(choices, key=lambda value: abs(value - 16))
        hit = self._series.hits[carbon]
        wanted = [57, 71, 85, 14 * carbon + 2]
        ions = [ion for ion in wanted if self._run.mz[0] <= ion <= self._run.mz[-1]]
        images = np.stack([self._fold.image(self._run.xic(ion)) for ion in ions], axis=2)
        column = int(np.argmin(abs((self._fold.t1 + 0.5 * self._fold.spec.period) - hit.t1)))
        row = int(round(hit.t2 / self._fold.row_pitch))
        r0, r1 = max(row - 12, 0), min(row + 13, images.shape[0])
        c0, c1 = max(column - 2, 0), min(column + 3, images.shape[1])
        window = images[r0:r1, c0:c1].astype(float)
        rim = np.zeros(window.shape[:2], dtype=bool)
        rim[[0, -1], :] = True
        rim[:, [0, -1]] = True
        baseline = np.median(window[rim], axis=0)
        profile = np.maximum(window - baseline[None, None, :], 0.0).sum(axis=2)
        if profile.sum() <= 0:
            raise ValueError(f"C{carbon} has no positive profile after its local baseline")
        return profile / profile.sum(), f"n-C{carbon}"

    def _run_fit(self) -> None:
        if self.cube is None:
            return
        cube = self.cube
        model = self.model_combo.currentText()
        auto = self.auto_k.isChecked()
        k = self.k_spin.value()
        self.fit_status.setText("Fitting…")

        def work(token, progress):
            if model == "measured alkane template":
                template, source = self._measured_alkane_template()
                if auto:
                    fits = [fit_template_overlap(
                        cube, template, trial_k, include_background=True,
                        background_order=1, n_starts=5, seed=trial_k)
                            for trial_k in range(1, 5)]
                    fit = min(fits, key=lambda result: result.bic)
                else:
                    fit = fit_template_overlap(
                        cube, template, k, include_background=True,
                        background_order=1, n_starts=5)
                fit.template_source = source
                return fit
            if auto:
                return choose_k(cube, max_k=4, model=model, include_background=True,
                                background_order=1, n_starts=5)
            return fit_overlap(cube, k, model=model, include_background=True,
                               background_order=1, n_starts=5)

        def done(fit):
            self._last_fit = fit
            self.accept_fit_btn.setEnabled(self._source == "candidates")
            assessments = assess_components(fit)
            self.fit_table.setRowCount(len(fit.components))
            for r, (comp, assessment) in enumerate(zip(fit.components, assessments)):
                values = [str(r), f"{comp.t1 / 60:.3f}", f"{comp.t2:.3f}",
                         f"{comp.row0_se:.2f} / {comp.col0_se:.2f}", f"{comp.volume:,.0f}",
                         f"{100*assessment.relative_to_largest:.0f}%", assessment.confidence]
                for c, v in enumerate(values):
                    item = QtWidgets.QTableWidgetItem(v)
                    item.setToolTip(assessment.reason)
                    if assessment.confidence == "low":
                        item.setBackground(QtGui.QColor(255, 160, 40, 75))
                    self.fit_table.setItem(r, c, item)
            self.component_view.set_fit(cube, fit)
            source = (f", shape={fit.template_source}" if
                      hasattr(fit, "template_source") else "")
            low = sum(assessment.confidence == "low" for assessment in assessments)
            group_note = (f" Unresolved group: {low} weak estimate{'s' if low != 1 else ''} retained."
                          if low else "")
            self.fit_status.setText(f"k={fit.k}, BIC={fit.bic:,.1f}{source}"
                                    f"{' (converged)' if fit.converged else ' (did not fully converge)'}. "
                                    "Amounts below 30% of the largest component are retained but marked low confidence."
                                    + group_note)

        def failed(text):
            self.fit_status.setText(f"Failed: {text}")

        self.project.jobs.submit("✂️ Fitting peak shapes", work, done, failed)

    def _accept_fit(self) -> None:
        region_index = self.blob_combo.currentData()
        if (self._source != "candidates" or self._last_fit is None or self.cube is None
                or self._name is None or region_index is None):
            return
        group = group_from_fit(self._name, int(region_index), self._last_fit, self.cube)
        self.fitAccepted.emit(self._name, group)
        self.fit_status.setText(self.fit_status.text() + " Accepted into the project.")
        self.accept_fit_btn.setEnabled(False)

    # ------------------------------------------------------------------ 4
    def _run_mcr(self) -> None:
        if self.cube is None:
            return
        cube = self.cube
        k = self.mcr_k.value()
        unimodal = self.mcr_unimodal.isChecked()
        self.mcr_status.setText("Running…")

        def work(token, progress):
            return mcr_als(cube, k, unimodal=unimodal)

        def done(result):
            self.mcr_table.setRowCount(result.spectra.shape[0])
            for r in range(result.spectra.shape[0]):
                dom = int(np.argmax(result.spectra[r]))
                volume = float(result.spectra[r].sum() * result.concentrations[..., r].sum())
                values = [str(r), f"{result.ions[dom]:g}", f"{volume:,.0f}"]
                for c, v in enumerate(values):
                    self.mcr_table.setItem(r, c, QtWidgets.QTableWidgetItem(v))
            self.mcr_status.setText(f"{result.n_iter} iteration(s)"
                                    f"{' (converged)' if result.converged else ' (did not fully converge)'}, "
                                    f"residual {result.resid:.3g}")

        def failed(text):
            self.mcr_status.setText(f"Failed: {text}")

        self.project.jobs.submit("✂️ MCR-ALS", work, done, failed)
