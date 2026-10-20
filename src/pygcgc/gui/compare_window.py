"""⚖️ Compare window: two runs, by image, by pixel spectrum and by peak.

Opened from the toolbar (⚖️ Compare, key 5). Choose run A and run B at the
top; every tab uses that pair. A is the grid everything is drawn on.

- Images: blink, grid, difference and overlay of the current ion expression.
- Spectra: a map of how alike the two runs' mass spectra are at every
  pixel, and the two spectra at any pixel clicked, head to tail.
- Peaks: the detected peaks of A and B paired one to one, by identity
  first (🏷️ internal standards, confident 🏷️ library hits) and then by
  position and spectrum, with the volume ratio of each shared peak and the
  peaks found in one run only.
- Classes: each 🧬 compound class's share of the two runs.

📍 alignment is used when both runs have it, and is not needed otherwise:
B is then read at A's own retention times. The science is in
`core/compare.py`.

A window of its own, not a step or a docked panel, because comparing
changes nothing: like 🏷️ and 🗂️ there is nothing to undo or save
(docs/architecture.md section 5).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets
from .theme import fit_to_screen

from ..core.compare import (MAX_DT1, MAX_DT2, MIN_SIMILARITY, class_shares, common_axis, match_peaks,
                            attach_labels, pair_mapping, pixel_spectrum, resample, similarity,
                            similarity_map, source_pixels, summarise_pairs, write_pairs_csv)
from ..core.peaktable import respectrum_components
from .blob_state import MS_AWARE
from .views.compare import CompareView

PEAK_COLUMNS = ["In", "Name", "Class A", "Class B", "t₁ A (min)", "t₂ A (s)", "Δt₁ (s)", "Δt₂ (s)",
                "Similarity", "Volume A", "Volume B", "log₂ B/A", "Paired by"]
FILTERS = {"all": "All peaks", "both": "In both", "only A": "Only in A", "only B": "Only in B",
           "named": "Identified only"}
CLASS_COLUMNS = ["Class", "Peaks A", "Peaks B", "% of A", "% of B", "B − A (points)"]
COLOUR_A = (230, 159, 0)        # Okabe-Ito orange
COLOUR_B = (0, 114, 178)        # Okabe-Ito blue


def _num(x: float, fmt: str) -> str:
    return "" if x is None or not np.isfinite(x) else format(x, fmt)


class _NumberItem(QtWidgets.QTableWidgetItem):
    """A table cell shown formatted but sorted by value (blanks last)."""

    def __init__(self, value: float, fmt: str):
        super().__init__(_num(value, fmt))
        self.value = float(value) if value is not None and np.isfinite(value) else np.inf

    def __lt__(self, other):
        return self.value < other.value if isinstance(other, _NumberItem) else super().__lt__(other)


class CompareWindow(QtWidgets.QDialog):
    navigateRequested = QtCore.Signal(str, float, float)     # run name, t1 (s), t2 (s)

    def __init__(self, project, expr_source, parent=None):
        super().__init__(parent)
        self.setWindowTitle("⚖️ Compare runs")
        self.setWindowFlag(QtCore.Qt.WindowType.WindowMaximizeButtonHint, True)
        fit_to_screen(self, 1200, 820)
        self.project = project
        self.align = project.align
        self.expr_source = expr_source       # () -> the main window's IonExpr
        self._peaks_cache: dict = {}         # (run name, id of its detection result) -> peaks with spectra
        self.pairs: list = []
        self._pairs_for: tuple | None = None
        self.sim = None
        self._sim_for: tuple | None = None

        self.a_box = QtWidgets.QComboBox()
        self.b_box = QtWidgets.QComboBox()
        for box in (self.a_box, self.b_box):
            box.setMinimumWidth(220)
            box.currentIndexChanged.connect(self._pair_changed)
        swap = QtWidgets.QPushButton("⇄")
        swap.setToolTip("Swap A and B")
        swap.setFixedWidth(34)
        swap.clicked.connect(self.swap)
        self.align_label = QtWidgets.QLabel("")
        self.align_label.setWordWrap(True)
        top = QtWidgets.QHBoxLayout()
        top.addWidget(QtWidgets.QLabel("A"))
        top.addWidget(self.a_box)
        top.addWidget(swap)
        top.addWidget(QtWidgets.QLabel("B"))
        top.addWidget(self.b_box)
        top.addWidget(self.align_label, 1)

        self.images = CompareView()
        self.images.image_source = self._image
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self.images, "Images")
        self.tabs.addTab(self._build_spectra_tab(), "Spectra")
        self.tabs.addTab(self._build_peaks_tab(), "Peaks")
        self.tabs.addTab(self._build_classes_tab(), "Classes")
        self.tabs.currentChanged.connect(lambda _i: self._refresh_tab())

        lay = QtWidgets.QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(self.tabs, 1)

        project.runsChanged.connect(self._runs_changed)
        project.layersChanged.connect(self._data_changed)
        project.foldChanged.connect(lambda _i: self._data_changed())
        self.align.changed.connect(self._data_changed)
        self._runs_changed()

    # ================================================================ building
    def _build_spectra_tab(self) -> QtWidgets.QWidget:
        self.sim_btn = QtWidgets.QPushButton("Map spectral similarity")
        self.sim_btn.setToolTip("Score A's mass spectrum against B's at every pixel (all shared m/z). "
                                "Takes a few seconds per run pair.")
        self.sim_btn.clicked.connect(self.compute_similarity)
        self.floor = QtWidgets.QDoubleSpinBox()
        self.floor.setRange(0.0, 50.0)
        self.floor.setDecimals(2)
        # 5 %: on the real exhaust runs, pixels below it are mostly noise, whose spectra
        # score low against anything and would read as a chemical difference.
        self.floor.setValue(5.0)
        self.floor.setSuffix(" %")
        self.floor.setToolTip("Hide pixels whose signal, in the brighter of the two runs, is below this share "
                              "of the brightest pixel (99.9th percentile): noise has no meaningful spectrum")
        self.floor.valueChanged.connect(lambda _v: self._show_similarity())
        self.sim_status = QtWidgets.QLabel("Click the image to compare the two spectra at that point.")
        self.sim_status.setWordWrap(True)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.sim_btn)
        row.addWidget(QtWidgets.QLabel("Hide below"))
        row.addWidget(self.floor)
        row.addWidget(self.sim_status, 1)

        self.sim_canvas = pg.GraphicsLayoutWidget()
        self.sim_plot = self.sim_canvas.addPlot()
        self.sim_plot.setLabel("bottom", "t₁ (min)")
        self.sim_plot.setLabel("left", "t₂ (s)")
        self.sim_plot.setMenuEnabled(False)
        self.sim_plot.hideButtons()
        self.sim_image = pg.ImageItem()
        self.sim_plot.addItem(self.sim_image)
        self.sim_bar = pg.ColorBarItem(values=(0.0, 1.0), colorMap=pg.colormap.get("viridis"),
                                       label="similarity", interactive=False, width=12)
        self.sim_bar.setImageItem(self.sim_image, insert_in=self.sim_plot)
        self.sim_bar.setVisible(False)
        self.sim_marker = pg.ScatterPlotItem(size=12, symbol="+", pen=pg.mkPen((255, 255, 255), width=2))
        self.sim_plot.addItem(self.sim_marker)
        self.sim_canvas.scene().sigMouseClicked.connect(self._sim_clicked)

        self.mirror = pg.PlotWidget()
        self.mirror.setMenuEnabled(False)
        self.mirror.hideButtons()
        self.mirror.setLabel("bottom", "m/z")
        self.mirror.setLabel("left", "relative intensity (%)")
        self.mirror.addLine(y=0, pen=pg.mkPen((128, 128, 128)))
        self.sticks_a = pg.PlotDataItem([], [], connect="pairs", pen=pg.mkPen(COLOUR_A, width=2))
        self.sticks_b = pg.PlotDataItem([], [], connect="pairs", pen=pg.mkPen(COLOUR_B, width=2))
        self.mirror.addItem(self.sticks_a)
        self.mirror.addItem(self.sticks_b)
        self.mirror_label = QtWidgets.QLabel("")
        self.mirror_label.setWordWrap(True)
        self.mirror_label.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)

        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        split.addWidget(self.sim_canvas)
        lower = QtWidgets.QWidget()
        ll = QtWidgets.QVBoxLayout(lower)
        ll.setContentsMargins(0, 0, 0, 0)
        ll.addWidget(self.mirror_label)
        ll.addWidget(self.mirror, 1)
        split.addWidget(lower)
        split.setSizes([520, 260])
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addLayout(row)
        lay.addWidget(split, 1)
        return w

    def _build_peaks_tab(self) -> QtWidgets.QWidget:
        self.max_dt1 = QtWidgets.QDoubleSpinBox()
        self.max_dt1.setRange(0.5, 300.0)
        self.max_dt1.setValue(MAX_DT1)
        self.max_dt1.setSuffix(" s")
        self.max_dt1.setToolTip("Largest first-dimension distance (after 📍, if used) at which two peaks can pair")
        self.max_dt2 = QtWidgets.QDoubleSpinBox()
        self.max_dt2.setRange(0.01, 5.0)
        self.max_dt2.setDecimals(2)
        self.max_dt2.setSingleStep(0.05)
        self.max_dt2.setValue(MAX_DT2)
        self.max_dt2.setSuffix(" s")
        self.max_dt2.setToolTip("Largest second-dimension distance at which two peaks can pair")
        self.min_sim = QtWidgets.QDoubleSpinBox()
        self.min_sim.setRange(0.0, 1.0)
        self.min_sim.setDecimals(2)
        self.min_sim.setSingleStep(0.05)
        self.min_sim.setValue(MIN_SIMILARITY)
        self.min_sim.setToolTip("Least spectral similarity (weighted cosine squared, as 🏷️ uses) for two "
                                "unnamed peaks to pair. Peaks with the same identity pair regardless.")
        match_btn = QtWidgets.QPushButton("Pair peaks")
        match_btn.setToolTip("Pair A's detected peaks with B's (needs 🔍 peak detection on both)")
        match_btn.clicked.connect(self.match)
        self.match_btn = match_btn
        self.filter_box = QtWidgets.QComboBox()
        for key, text in FILTERS.items():
            self.filter_box.addItem(text, key)
        self.filter_box.currentIndexChanged.connect(lambda _i: self._fill_peaks())
        export_btn = QtWidgets.QPushButton("Export CSV…")
        export_btn.clicked.connect(self._export_pairs)
        controls = QtWidgets.QHBoxLayout()
        for label, w in (("Δt₁ ≤", self.max_dt1), ("Δt₂ ≤", self.max_dt2), ("similarity ≥", self.min_sim)):
            controls.addWidget(QtWidgets.QLabel(label))
            controls.addWidget(w)
        controls.addWidget(match_btn)
        controls.addStretch(1)
        controls.addWidget(QtWidgets.QLabel("Show"))
        controls.addWidget(self.filter_box)
        controls.addWidget(export_btn)
        self.peak_status = QtWidgets.QLabel("Press Pair peaks.")
        self.peak_status.setWordWrap(True)
        self.peak_status.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)

        self.peak_table = QtWidgets.QTableWidget(0, len(PEAK_COLUMNS))
        self.peak_table.setHorizontalHeaderLabels(PEAK_COLUMNS)
        self.peak_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.peak_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.peak_table.setSortingEnabled(True)
        self.peak_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.peak_table.setToolTip("Double-click a row to show that peak in the main window")
        self.peak_table.cellDoubleClicked.connect(self._peak_double_clicked)

        self.scatter = pg.PlotWidget()
        self.scatter.setMenuEnabled(False)
        self.scatter.hideButtons()
        self.scatter.setLabel("bottom", "log₁₀ volume in A")
        self.scatter.setLabel("left", "log₁₀ volume in B")
        self.scatter_points = pg.ScatterPlotItem(size=6, pen=pg.mkPen((40, 40, 40), width=0.6),
                                                 brush=pg.mkBrush(*COLOUR_B, 160))
        self.scatter_named = pg.ScatterPlotItem(size=8, symbol="s", pen=pg.mkPen((40, 40, 40), width=0.6),
                                                brush=pg.mkBrush(*COLOUR_A, 220))
        self.one_to_one = pg.PlotDataItem([], [], pen=pg.mkPen((128, 128, 128), style=QtCore.Qt.PenStyle.DashLine))
        for item in (self.one_to_one, self.scatter_points, self.scatter_named):
            self.scatter.addItem(item)
        self.scatter.setMinimumWidth(280)

        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        split.addWidget(self.peak_table)
        split.addWidget(self.scatter)
        split.setSizes([820, 340])
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addLayout(controls)
        lay.addWidget(self.peak_status)
        lay.addWidget(split, 1)
        return w

    def _build_classes_tab(self) -> QtWidgets.QWidget:
        self.class_status = QtWidgets.QLabel("")
        self.class_status.setWordWrap(True)
        self.class_table = QtWidgets.QTableWidget(0, len(CLASS_COLUMNS))
        self.class_table.setHorizontalHeaderLabels(CLASS_COLUMNS)
        self.class_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.class_table.setSortingEnabled(True)
        self.class_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        w = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addWidget(self.class_status)
        lay.addWidget(self.class_table, 1)
        return w

    def set_background(self, colour) -> None:
        self.images.set_background(colour)
        self.sim_canvas.setBackground(colour)
        self.mirror.setBackground(colour)
        self.scatter.setBackground(colour)

    # ================================================================ runs
    def _runs_changed(self) -> None:
        p = self.project
        names = [p.entries[i].name for i in p.loaded()]
        keep_a, keep_b = self.a_box.currentData(), self.b_box.currentData()
        for box in (self.a_box, self.b_box):
            box.blockSignals(True)
            box.clear()
            for n in names:
                box.addItem(n, n)
            box.blockSignals(False)
        a = keep_a if keep_a in names else (self.align.reference if self.align.reference in names else
                                            (names[0] if names else None))
        b = keep_b if keep_b in names else next((n for n in names if n != a), a)
        self.set_pair(a, b)

    def set_pair(self, a: str | None, b: str | None) -> None:
        """Choose runs A and B by name (the main window and tests call this)."""
        for box, name in ((self.a_box, a), (self.b_box, b)):
            i = box.findData(name) if name is not None else -1
            box.blockSignals(True)
            box.setCurrentIndex(i)
            box.blockSignals(False)
        self._pair_changed()

    def swap(self) -> None:
        self.set_pair(self.b_box.currentData(), self.a_box.currentData())

    def index_a(self) -> int:
        i = self.project.index_of(self.a_box.currentData()) if self.a_box.currentData() else None
        return -1 if i is None else i

    def index_b(self) -> int:
        i = self.project.index_of(self.b_box.currentData()) if self.b_box.currentData() else None
        return -1 if i is None else i

    def _pair_changed(self, *_a) -> None:
        self._update_align_label()
        if self.isVisible():
            self._refresh_tab()

    def showEvent(self, event) -> None:
        super().showEvent(event)
        self._refresh_tab()

    def _data_changed(self, *_a) -> None:
        """Processing, folds or alignment changed: what was computed for this pair may be out of date."""
        self._sim_for = None
        self._pairs_for = None
        self._peaks_cache.clear()
        self._update_align_label()
        if self.isVisible():
            self._refresh_tab()

    def expression_changed(self) -> None:
        if self.isVisible() and self.tabs.currentWidget() is self.images:
            self.images.refresh()

    def show_images(self, mode: str | None = None) -> None:
        self.tabs.setCurrentWidget(self.images)
        if mode is not None:
            i = self.images.mode.findData(mode)
            if i >= 0:
                self.images.mode.setCurrentIndex(i)

    # ================================================================ alignment
    def _transform(self, name: str):
        return None if name == self.align.reference else self.align.transforms.get(name)

    def aligned(self, name_a: str | None = None, name_b: str | None = None) -> bool:
        """Whether both runs sit in 📍's frame: the reference, or a run with a transform."""
        name_a = name_a or self.a_box.currentData()
        name_b = name_b or self.b_box.currentData()
        runs = self.align.aligned_runs()
        return bool(self.align.reference) and name_a in runs and name_b in runs

    def mapping(self, name_from: str, name_to: str):
        """(t1, t2) in `name_from`'s coordinates -> `name_to`'s, through 📍 if both are aligned."""
        if name_from != name_to and self.aligned(name_from, name_to):
            return pair_mapping(self._transform(name_from), self._transform(name_to))
        return pair_mapping(None, None)

    def _update_align_label(self) -> None:
        a, b = self.a_box.currentData(), self.b_box.currentData()
        if not a or not b:
            text = "Open at least one run."
        elif a == b:
            text = "A and B are the same run."
        elif self.aligned(a, b):
            text = f"📍 Aligned: B is read through the alignment to reference {self.align.reference}."
        else:
            text = ("Not aligned: B is read at A's own retention times. Fine for small drift; "
                    "📍 Align both runs if peaks look doubled in the overlay.")
        self.align_label.setText(text)

    # ================================================================ images
    def _image(self, index: int):
        """Run `index`'s image of the main window's ion expression, on A's grid."""
        p, ia = self.project, self.index_a()
        if not (0 <= index < len(p.entries)) or ia < 0 or p.entries[index].run is None:
            return None
        fold_a, expr = p.fold(ia), self.expr_source()
        img = expr.image(p.display_run(index), p.fold(index))
        if index != ia:
            name, name_a = p.entries[index].name, p.entries[ia].name
            img = resample(img, source_pixels(p.fold(index), fold_a, self.mapping(name_a, name)))
        return img, fold_a.t1, fold_a.row_pitch, fold_a.spec.period, p.entries[index].name

    # ================================================================ tabs
    def _refresh_tab(self) -> None:
        w = self.tabs.currentIndex()
        if w == 0:
            self.images.set_runs(self.index_a(), self.index_b(), self.project.loaded())
        else:
            self.images.blink_timer.stop()
        if w == 1:
            self._show_similarity()
        elif w == 2:
            self._fill_peaks()
        elif w == 3:
            self._fill_classes()

    def hideEvent(self, event) -> None:
        self.images.blink_timer.stop()
        super().hideEvent(event)

    # ================================================================ spectra
    def _pair_key(self) -> tuple:
        return (self.a_box.currentData(), self.b_box.currentData(), self.aligned())

    def compute_similarity(self) -> None:
        p, ia, ib = self.project, self.index_a(), self.index_b()
        if ia < 0 or ib < 0:
            self.sim_status.setText("Choose two loaded runs.")
            return
        run_a, run_b, fold_a, fold_b = p.display_run(ia), p.display_run(ib), p.fold(ia), p.fold(ib)
        to_b, key = self.mapping(p.entries[ia].name, p.entries[ib].name), self._pair_key()
        self.sim_btn.setEnabled(False)
        self.sim_status.setText("Mapping spectral similarity…")

        def work(token, progress):
            return similarity_map(run_a, fold_a, run_b, fold_b, to_b, progress=progress, cancel=token.check)

        def done(result):
            self.sim_btn.setEnabled(True)
            self.sim, self._sim_for = result, key
            self._show_similarity()

        def failed(text):
            self.sim_btn.setEnabled(True)
            self.sim_status.setText(f"Failed: {text}")

        p.jobs.submit(f"⚖️ Spectral similarity: {key[0]} vs {key[1]}", work, done, failed)

    def _show_similarity(self) -> None:
        p, ia = self.project, self.index_a()
        if ia < 0 or p.entries[ia].run is None:
            self.sim_image.clear()
            return
        fold = p.fold(ia)
        if self.sim is not None and self._sim_for == self._pair_key():
            s = self.sim
            bright = np.maximum(s.tic_a, s.tic_b)
            top = float(np.percentile(bright, 99.9)) if bright.size else 0.0
            shown = np.where(bright >= self.floor.value() / 100.0 * top, s.similarity, np.nan)
            self.sim_image.setImage(shown, autoLevels=False)
            self.sim_bar.setVisible(True)
            self.sim_bar.setColorMap(pg.colormap.get("viridis"))
            self.sim_bar.setLevels((0.0, 1.0))
            finite = shown[np.isfinite(shown)]
            weights = bright[np.isfinite(shown)]
            mean = float(np.average(finite, weights=weights)) if finite.size and weights.sum() > 0 else np.nan
            self.sim_status.setText(
                f"Spectral similarity of A and B at every pixel ({s.mz[0]}–{s.mz[-1]} m/z; 1 = same spectrum). "
                f"Signal-weighted mean {_num(mean, '.2f')}. Click a point to compare its spectra.")
        else:
            tic = fold.tic_image(p.display_run(ia))
            self.sim_image.setImage(np.sqrt(np.clip(tic, 0, None)), autoLevels=True)
            self.sim_image.setColorMap(pg.colormap.get("CET-L1"))
            self.sim_bar.setVisible(False)
            self.sim_status.setText("Run A's signal. Click a point to compare the two spectra there, "
                                    "or map the similarity everywhere.")
        self.sim_image.setRect(QtCore.QRectF(fold.t1[0] / 60.0, -0.5 * fold.row_pitch,
                                             fold.n_cols * fold.spec.period / 60.0, fold.n_rows * fold.row_pitch))

    def _sim_clicked(self, ev) -> None:
        vb = self.sim_plot.vb
        if ev.button() != QtCore.Qt.MouseButton.LeftButton or not vb.sceneBoundingRect().contains(ev.scenePos()):
            return
        pt = vb.mapSceneToView(ev.scenePos())
        self.show_pixel(pt.x() * 60.0, pt.y())

    def show_pixel(self, t1: float, t2: float) -> float | None:
        """Head-to-tail spectra of A and B at A's (t1 in s, t2 in s); returns their similarity."""
        p, ia, ib = self.project, self.index_a(), self.index_b()
        if ia < 0 or ib < 0:
            return None
        run_a, run_b = p.display_run(ia), p.display_run(ib)
        name_a, name_b = p.entries[ia].name, p.entries[ib].name
        t1_b, t2_b = self.mapping(name_a, name_b)(np.array([t1]), np.array([t2]))
        spec_a = pixel_spectrum(run_a, p.fold(ia), t1, t2)
        spec_b = pixel_spectrum(run_b, p.fold(ib), float(t1_b[0]), float(t2_b[0]))
        mz, sa, sb = common_axis(run_a.mz, run_b.mz)
        spec_a, spec_b = spec_a[sa], spec_b[sb]
        s = similarity(mz, spec_a, spec_b)
        for item, spec, sign in ((self.sticks_a, spec_a, 1.0), (self.sticks_b, spec_b, -1.0)):
            top = float(spec.max())
            y = sign * 100.0 * spec / top if top > 0 else np.zeros_like(spec)
            keep = np.abs(y) >= 0.5
            item.setData(np.repeat(mz[keep].astype(float), 2), np.column_stack([np.zeros(keep.sum()), y[keep]]).ravel())
        self.sim_marker.setData([t1 / 60.0], [t2])
        self.mirror_label.setText(
            f"t₁ {t1 / 60:.3f} min, t₂ {t2:.2f} s. A ({name_a}, up, orange; total {spec_a.sum():.3g}) against "
            f"B ({name_b}, down, blue; total {spec_b.sum():.3g}): similarity {s:.2f}. Each spectrum scaled to its "
            "own largest ion; three scans summed.")
        return s

    # ================================================================ peaks
    def _peaks_for(self, index: int):
        """(cache key, a callable giving this run's peaks with full spectra), or None without 🔍 results.
        The callable does the slow part and runs inside the job."""
        p = self.project
        blobs = p.blobs.blobs_for(index)
        name = p.entries[index].name
        if blobs is not None:
            return (name, id(blobs)), lambda: blobs
        if p.blobs.method == MS_AWARE:
            component_set = p.blobs.components_for(index)
            if component_set is not None and component_set.groups:
                base_run, fold = p.base_run(index), p.fold(index)
                return (name, id(component_set)), lambda: respectrum_components(base_run, fold, component_set)[0]
        return None

    def _labels(self, index: int) -> tuple[list, list]:
        """Confident identities and 🧬 class calls on one run, as (t1, t2, text) lists."""
        p = self.project
        name = p.entries[index].name
        names = [(m.t1, m.t2, m.name) for m in p.istd.matches_for(index).values() if m.found]
        names += [(t1, t2, text) for t1, t2, text in p.identifications.get(name, [])]
        classes = []
        result, blobs = p.classes.result_for(index), p.classes.blobs_classified(index)
        if result is not None and blobs is not None and len(blobs) == len(result):
            classes = [(float(blobs.t1[i]), float(blobs.t2[i]), c.cls) for i, c in enumerate(result.calls) if c.cls]
        return names, classes

    def match(self) -> None:
        p, ia, ib = self.project, self.index_a(), self.index_b()
        if ia < 0 or ib < 0 or ia == ib:
            self.peak_status.setText("Choose two different loaded runs.")
            return
        src_a, src_b = self._peaks_for(ia), self._peaks_for(ib)
        missing = [p.entries[i].name for i, s in ((ia, src_a), (ib, src_b)) if s is None]
        if missing:
            self.peak_status.setText(f"Run 🔍 peak detection on {' and '.join(missing)} first.")
            return
        (key_a, get_a), (key_b, get_b) = src_a, src_b
        cache = self._peaks_cache
        names_a, classes_a = self._labels(ia)
        names_b, classes_b = self._labels(ib)
        tol_a = (p.fold(ia).spec.period, 0.1)
        tol_b = (p.fold(ib).spec.period, 0.1)
        to_a = self.mapping(p.entries[ib].name, p.entries[ia].name)
        settings = (self.max_dt1.value(), self.max_dt2.value(), self.min_sim.value())
        key = self._pair_key() + settings
        self.match_btn.setEnabled(False)
        self.peak_status.setText("Pairing peaks…")

        def work(token, progress):
            progress(0, 3, "Peaks of A")
            peaks_a = cache[key_a] if key_a in cache else get_a()
            token.check()
            progress(1, 3, "Peaks of B")
            peaks_b = cache[key_b] if key_b in cache else get_b()
            token.check()
            progress(2, 3, "Pairing")
            pairs = match_peaks(peaks_a, peaks_b, to_a,
                                attach_labels(peaks_a.t1, peaks_a.t2, names_a, *tol_a),
                                attach_labels(peaks_b.t1, peaks_b.t2, names_b, *tol_b),
                                max_dt1=settings[0], max_dt2=settings[1], min_similarity=settings[2])
            cls_a = attach_labels(peaks_a.t1, peaks_a.t2, classes_a, *tol_a)
            cls_b = attach_labels(peaks_b.t1, peaks_b.t2, classes_b, *tol_b)
            for pair in pairs:
                pair.class_a = cls_a[pair.a] if pair.a is not None else ""
                pair.class_b = cls_b[pair.b] if pair.b is not None else ""
            return peaks_a, peaks_b, pairs

        def done(result):
            peaks_a, peaks_b, pairs = result
            cache[key_a], cache[key_b] = peaks_a, peaks_b
            self.pairs, self._pairs_for = pairs, key
            self.match_btn.setEnabled(True)
            self._fill_peaks()

        def failed(text):
            self.match_btn.setEnabled(True)
            self.peak_status.setText(f"Failed: {text}")

        p.jobs.submit(f"⚖️ Pairing peaks: {key[0]} vs {key[1]}", work, done, failed)

    def _current_pairs(self) -> list:
        settings = (self.max_dt1.value(), self.max_dt2.value(), self.min_sim.value())
        return self.pairs if self._pairs_for == self._pair_key() + settings else []

    def _fill_peaks(self) -> None:
        pairs = self._current_pairs()
        want = self.filter_box.currentData()
        shown = [(k, q) for k, q in enumerate(pairs)
                 if want == "all" or (want == "named" and q.name) or q.status == want]
        t = self.peak_table
        t.setSortingEnabled(False)
        t.setRowCount(len(shown))
        for r, (k, q) in enumerate(shown):
            in_a = q.a is not None
            name = (q.name_a if not q.name_b or q.name_b == q.name_a else
                    q.name_b if not q.name_a else f"{q.name_a} / {q.name_b}")
            cells = [(q.status, None), (name, None), (q.class_a, None), (q.class_b, None),
                     ((q.t1_a if in_a else q.t1_b) / 60.0, ".3f"), (q.t2_a if in_a else q.t2_b, ".2f"),
                     (q.dt1, ".1f"), (q.dt2, ".3f"), (q.similarity, ".2f"),
                     (q.volume_a if in_a else np.nan, ".3g"), (q.volume_b if q.b is not None else np.nan, ".3g"),
                     (q.log2_ratio, "+.2f"), (q.how, None)]
            for c, (v, fmt) in enumerate(cells):
                item = QtWidgets.QTableWidgetItem(v) if fmt is None else _NumberItem(v, fmt)
                if c == 0:
                    item.setData(QtCore.Qt.ItemDataRole.UserRole, k)
                t.setItem(r, c, item)
        t.setSortingEnabled(True)

        both = [q for q in pairs if q.status == "both" and q.volume_a > 0 and q.volume_b > 0]
        named = [q for q in both if q.name]
        plain = [q for q in both if not q.name]
        self.scatter_points.setData([np.log10(q.volume_a) for q in plain], [np.log10(q.volume_b) for q in plain])
        self.scatter_named.setData([np.log10(q.volume_a) for q in named], [np.log10(q.volume_b) for q in named])
        if both:
            v = np.log10([x for q in both for x in (q.volume_a, q.volume_b)])
            self.one_to_one.setData([v.min(), v.max()], [v.min(), v.max()])
        else:
            self.one_to_one.setData([], [])
        if not pairs:
            self.peak_status.setText("Press Pair peaks (needs 🔍 peak detection on both runs).")
            return
        s = summarise_pairs(pairs)
        med = s["median_log2_ratio"]
        self.peak_status.setText(
            f"{s['both']} peaks in both, {s['only A']} only in A, {s['only B']} only in B; {s['named']} identified "
            f"(🏷️ internal standards and tier 2a/2b library hits). Median log₂ B/A of the shared peaks "
            f"{_num(med, '+.2f')} (B holds {_num(2 ** med if np.isfinite(med) else np.nan, '.2f')} times A). "
            "Squares on the plot are identified peaks; the dashed line is 1:1.")

    def _peak_double_clicked(self, row: int, _col: int) -> None:
        item = self.peak_table.item(row, 0)
        pairs = self._current_pairs()
        if item is None or not pairs:
            return
        q = pairs[int(item.data(QtCore.Qt.ItemDataRole.UserRole))]
        if q.a is not None:
            self.navigateRequested.emit(self.a_box.currentData(), q.t1_a, q.t2_a)
        else:
            self.navigateRequested.emit(self.b_box.currentData(), q.t1_b, q.t2_b)

    def _export_pairs(self) -> None:
        pairs = self._current_pairs()
        if not pairs:
            self.peak_status.setText("Pair the peaks first.")
            return
        a, b = self.a_box.currentData(), self.b_box.currentData()
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export peak comparison",
                                                        f"{Path(a).stem}_vs_{Path(b).stem}_peaks.csv", "CSV (*.csv)")
        if path:
            write_pairs_csv(pairs, path, a, b)
            self.peak_status.setText(f"Wrote {Path(path).name}")

    # ================================================================ classes
    def _fill_classes(self) -> None:
        p, ia, ib = self.project, self.index_a(), self.index_b()
        summaries = []
        for i in (ia, ib):
            result = p.classes.result_for(i) if i >= 0 else None
            blobs = p.classes.blobs_classified(i) if i >= 0 else None
            summaries.append(result.summary(blobs.volume) if result is not None and blobs is not None else None)
        t = self.class_table
        t.setSortingEnabled(False)
        if any(s is None for s in summaries):
            t.setRowCount(0)
            self.class_status.setText("Classify both runs (🧬 step 7) to compare their compound classes.")
            return
        rows = class_shares(*summaries)
        t.setRowCount(len(rows))
        for r, (name, na, nb, pa, pb, diff) in enumerate(rows):
            t.setItem(r, 0, QtWidgets.QTableWidgetItem(name))
            for c, (v, fmt) in enumerate(((na, "d"), (nb, "d"), (pa, ".1f"), (pb, ".1f"), (diff, "+.1f")), start=1):
                t.setItem(r, c, _NumberItem(v, fmt))
        t.setSortingEnabled(True)
        self.class_status.setText("Each class's share of the run's total peak volume, so runs of different total "
                                  "amount compare directly.")
