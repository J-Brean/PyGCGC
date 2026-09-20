"""🗂️ Batch dialog.

A one-shot, whole-folder analysis, deliberately not part of the per-run
project state the other tools live in (there's nothing to undo, and its
result isn't tied to which run happens to be open): pick a folder,
choose which steps to run, and go. Everything runs through the existing
`JobManager`, so the window stays responsive; the pipeline itself
(`core/batch.run_batch`) does the real parallel work in a process pool.
"""

from __future__ import annotations

import os
from dataclasses import replace
from pathlib import Path

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from .. import __version__
from ..core.batch import BatchConfig, BatchSettings, pca, run_batch, volcano, westgard_flags
from ..core.report import build_provenance, write_excel, write_html_report, write_pdf_report


class BatchDialog(QtWidgets.QDialog):
    navigateRequested = QtCore.Signal(str, float, float)   # run name, t1 (s), t2 (s)

    def __init__(self, project, default_folder=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("🗂️ Batch")
        self.resize(900, 620)
        self.project = project
        self.result = None
        self._active_job_id: int | None = None
        self._last_config: BatchConfig | None = None
        project.jobs.jobChanged.connect(self._job_changed)

        self.folder_edit = QtWidgets.QLineEdit(str(default_folder or ""))
        browse = QtWidgets.QPushButton("Browse…")
        browse.clicked.connect(self._browse)
        self.run_siloxane = QtWidgets.QCheckBox("🫧 Siloxane removal")
        self.run_blobs = QtWidgets.QCheckBox("🔍 Blob detection")
        self.run_blobs.setChecked(True)
        self.run_quant = QtWidgets.QCheckBox("🧪 Quantify")
        self.run_quant.setChecked(True)
        self.workers = QtWidgets.QSpinBox()
        self.workers.setRange(1, max(1, os.cpu_count() or 1))
        self.workers.setValue(max(1, (os.cpu_count() or 2) - 1))
        self.batch_settings = BatchSettings()
        self.reference_edit = QtWidgets.QLineEdit()
        self.reference_edit.setPlaceholderText("automatic (first run with blobs)")
        self.reference_edit.setToolTip("Run whose blobs define the feature table; the file stem, as in the run list")
        settings_btn = QtWidgets.QPushButton("⚙️ Settings…")
        settings_btn.setToolTip("Every setting of every step: 🧹 background, 🫧, 🔍 peaks and 🧪 quantification")
        settings_btn.clicked.connect(self._edit_settings)
        self.run_btn = QtWidgets.QPushButton("▶ Run batch")
        self.run_btn.clicked.connect(self._run_clicked)
        self.export_combo = QtWidgets.QComboBox()
        self.export_combo.addItems(["Excel (.xlsx)", "HTML report (.html)", "PDF report (.pdf)"])
        self.export_btn = QtWidgets.QPushButton("📤 Export…")
        self.export_btn.setEnabled(False)
        self.export_btn.setToolTip("Write the tables and charts currently shown, plus provenance "
                                   "(version, period, settings, any manual edits given), to a file.")
        self.export_btn.clicked.connect(self._export_clicked)
        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)

        top = QtWidgets.QFormLayout()
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.folder_edit, 1)
        row.addWidget(browse)
        top.addRow("Folder", row)
        steps = QtWidgets.QHBoxLayout()
        for w in (self.run_siloxane, self.run_blobs, self.run_quant):
            steps.addWidget(w)
        top.addRow("Steps (🧹 background always runs)", steps)
        top.addRow("Worker processes", self.workers)
        top.addRow("Feature-table reference", self.reference_edit)
        top.addRow("", settings_btn)

        self.tabs = QtWidgets.QTabWidget()
        self.feature_table = QtWidgets.QTableWidget(0, 0)
        self.feature_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.feature_table.cellDoubleClicked.connect(self._feature_row_activated)
        self.tabs.addTab(self.feature_table, "Feature table")

        qc_page = QtWidgets.QWidget()
        qc_lay = QtWidgets.QVBoxLayout(qc_page)
        self.qc_combo = QtWidgets.QComboBox()
        self.qc_combo.currentIndexChanged.connect(self._draw_qc)
        qc_lay.addWidget(self.qc_combo)
        self.qc_plot = pg.PlotWidget()
        self.qc_plot.setMenuEnabled(False)
        self.qc_plot.setLabel("left", "volume")
        self.qc_plot.setLabel("bottom", "run order")
        qc_lay.addWidget(self.qc_plot)
        self.qc_flags = QtWidgets.QLabel("")
        self.qc_flags.setWordWrap(True)
        qc_lay.addWidget(self.qc_flags)
        self.tabs.addTab(qc_page, "QC chart")

        self.pca_plot = pg.PlotWidget()
        self.pca_plot.setMenuEnabled(False)
        self.pca_plot.setLabel("bottom", "PC1")
        self.pca_plot.setLabel("left", "PC2")
        self.tabs.addTab(self.pca_plot, "PCA")

        volcano_page = QtWidgets.QWidget()
        v_lay = QtWidgets.QVBoxLayout(volcano_page)
        groups = QtWidgets.QHBoxLayout()
        self.group_a_list = QtWidgets.QListWidget()
        self.group_a_list.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.MultiSelection)
        self.group_b_list = QtWidgets.QListWidget()
        self.group_b_list.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.MultiSelection)
        groups.addWidget(_labelled("Group A", self.group_a_list))
        groups.addWidget(_labelled("Group B", self.group_b_list))
        v_lay.addLayout(groups)
        self.volcano_btn = QtWidgets.QPushButton("Compute volcano plot")
        self.volcano_btn.clicked.connect(self._draw_volcano)
        v_lay.addWidget(self.volcano_btn)
        self.volcano_plot = pg.PlotWidget()
        self.volcano_plot.setMenuEnabled(False)
        self.volcano_plot.setLabel("bottom", "log2 fold change (A / B)")
        self.volcano_plot.setLabel("left", "-log10(q)")
        v_lay.addWidget(self.volcano_plot)
        self.tabs.addTab(volcano_page, "Volcano")

        run_row = QtWidgets.QHBoxLayout()
        run_row.addWidget(self.run_btn, 1)
        run_row.addWidget(self.export_combo)
        run_row.addWidget(self.export_btn)

        lay = QtWidgets.QVBoxLayout(self)
        lay.addLayout(top)
        lay.addLayout(run_row)
        lay.addWidget(self.status)
        lay.addWidget(self.tabs, 1)

    def set_background(self, colour) -> None:
        for plot in (self.qc_plot, self.pca_plot, self.volcano_plot):
            plot.setBackground(colour)

    # ------------------------------------------------------------------ run
    def _edit_settings(self) -> None:
        from .dialogs import SettingsDialog
        dlg = SettingsDialog("🗂️ Batch settings", self.batch_settings, parent=self)
        if dlg.exec():
            self.batch_settings = dlg.settings()
            self.run_siloxane.setChecked(self.batch_settings.run_siloxane)
            self.run_blobs.setChecked(self.batch_settings.run_blobs)
            self.run_quant.setChecked(self.batch_settings.run_quant)

    def _browse(self) -> None:
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Folder to batch-process")
        if folder:
            self.folder_edit.setText(folder)

    def _run_clicked(self) -> None:
        folder = Path(self.folder_edit.text())
        if not folder.is_dir():
            self.status.setText("Not a folder.")
            return
        settings = replace(self.batch_settings, run_siloxane=self.run_siloxane.isChecked(),
                           run_blobs=self.run_blobs.isChecked(), run_quant=self.run_quant.isChecked())
        config = BatchConfig(settings=settings, reference=self.reference_edit.text().strip())
        self._last_config = config
        cache_dir = folder / ".pygcgc_cache"
        self.run_btn.setEnabled(False)
        self.export_btn.setEnabled(False)
        self.status.setText("Starting…")

        def work(token, progress):
            return run_batch(folder, config, cache_dir=cache_dir, max_workers=self.workers.value(),
                             progress=progress)

        def done(result):
            self._active_job_id = None
            self.run_btn.setEnabled(True)
            self.result = result
            n_err = len(result.errors)
            n_ok = len(result.prepped)
            msg = f"{n_ok} run{'s' if n_ok != 1 else ''} done"
            if n_err:
                msg += f", {n_err} failed: {', '.join(result.errors)}"
            self.status.setText(msg)
            self._fill_feature_table(result)
            self._fill_qc(result)
            self._fill_pca(result)
            self._fill_groups(result)
            self.export_btn.setEnabled(True)

        def failed(text):
            self._active_job_id = None
            self.run_btn.setEnabled(True)
            self.status.setText(f"Failed: {text}")

        self._active_job_id = self.project.jobs.submit("🗂️ Batch", work, done, failed)

    def _job_changed(self, job) -> None:
        if job.id == self._active_job_id:
            self.status.setText(f"{job.message} ({job.done}/{job.total})")

    # -------------------------------------------------------------- export
    def _export_clicked(self) -> None:
        kind = ["xlsx", "html", "pdf"][self.export_combo.currentIndex()]
        filters = {"xlsx": "Excel (*.xlsx)", "html": "HTML (*.html)", "pdf": "PDF (*.pdf)"}
        default = str(Path(self.folder_edit.text()) / f"report.{kind}")
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Export report", default, filters[kind])
        if path:
            self._do_export(Path(path), kind)

    def _do_export(self, path: Path, kind: str) -> None:
        """Write whatever the tabs are currently showing (the QC compound
        selected in the combo, the volcano groups selected, if any -- so
        the file matches what's on screen) plus provenance. The batch
        dialog has no interactive edit log of its own to add (see
        docs/architecture.md section 5: 🗂️ is deliberately not part of the per-run
        Project state that accumulates one), so `edits` is always empty
        here; only period, method and settings are recorded."""
        if self.result is None:
            return
        folder = Path(self.folder_edit.text())
        method_path = folder / "method.json"
        method = {}
        if method_path.exists():
            import json
            method = json.loads(method_path.read_text())
        period = (self._last_config.period if self._last_config else None) or method.get("modulation_period_s")
        settings = self._last_config.settings.to_dict() if self._last_config else {}
        prov = build_provenance(__version__, period=period, method=method, settings=settings,
                                runs=[p.run for p in self.result.prepped.values()])

        run_names = list(self.result.prepped)
        qc_values, qc_label = None, ""
        if self.qc_combo.count():
            compound = self.qc_combo.currentData()
            rows = [r for r in self.result.report.istd if r.compound == compound] \
                if self.result.report is not None else []
            if rows:
                qc_values = np.array([r.volume for r in rows])
                qc_label = self.qc_combo.currentText()
        pca_result = None
        if self.result.feature_table and len(run_names) >= 2:
            pca_result = pca(self.result.feature_table, run_names, n_components=2)
        volcano_points = None
        group_a = [i.text() for i in self.group_a_list.selectedItems()]
        group_b = [i.text() for i in self.group_b_list.selectedItems()]
        if self.result.feature_table and len(group_a) >= 2 and len(group_b) >= 2:
            volcano_points = volcano(self.result.feature_table, group_a, group_b)

        try:
            if kind == "xlsx":
                write_excel(path, report=self.result.report, feature_table=self.result.feature_table,
                           run_names=run_names, provenance=prov)
            elif kind == "html":
                write_html_report(path, prov, report=self.result.report,
                                 feature_table=self.result.feature_table, run_names=run_names,
                                 qc_values=qc_values, qc_label=qc_label, pca_result=pca_result,
                                 volcano_points=volcano_points)
            else:
                write_pdf_report(path, prov, report=self.result.report,
                                feature_table=self.result.feature_table, qc_values=qc_values,
                                qc_label=qc_label, pca_result=pca_result, volcano_points=volcano_points)
        except ImportError as exc:
            self.status.setText(f"Export needs an extra that isn't installed: {exc}")
            return
        self.status.setText(f"Wrote {path}")

    # ------------------------------------------------------------ feature table
    def _fill_feature_table(self, result) -> None:
        ft = result.feature_table or []
        run_names = list(result.prepped)
        cols = ["#", "t1 (min)", "t2 (s)"] + run_names
        self.feature_table.setColumnCount(len(cols))
        self.feature_table.setHorizontalHeaderLabels(cols)
        self.feature_table.setRowCount(len(ft))
        for i, f in enumerate(ft):
            values = [str(f.id), f"{f.ref_t1 / 60:.3f}", f"{f.ref_t2:.3f}"]
            for name in run_names:
                v = f.volumes.get(name)
                values.append("" if v is None else (f"{v:,.0f}*" if f.filled.get(name) else f"{v:,.0f}"))
            for c, v in enumerate(values):
                self.feature_table.setItem(i, c, QtWidgets.QTableWidgetItem(v))
        self.feature_table.resizeColumnsToContents()

    def _feature_row_activated(self, row: int, _col: int) -> None:
        if not self.result or not self.result.feature_table:
            return
        f = self.result.feature_table[row]
        # Jump to whichever run's own column was double-clicked, if it's one
        # of the run columns; otherwise the reference run.
        col_names = [self.feature_table.horizontalHeaderItem(c).text()
                    for c in range(self.feature_table.columnCount())]
        run_name = col_names[_col] if _col >= 3 else self.result.reference
        if run_name not in self.result.prepped:
            run_name = self.result.reference
        self.navigateRequested.emit(run_name, f.ref_t1, f.ref_t2)

    # ------------------------------------------------------------------- QC
    def _fill_qc(self, result) -> None:
        self.qc_combo.blockSignals(True)
        self.qc_combo.clear()
        if result.report is not None:
            for row in result.report.istd:
                label = f"{row.compound} (IS)"
                if self.qc_combo.findText(label) < 0:
                    self.qc_combo.addItem(label, row.compound)
        self.qc_combo.blockSignals(False)
        self._draw_qc()

    def _draw_qc(self) -> None:
        self.qc_plot.clear()
        self.qc_flags.setText("")
        if self.result is None or self.result.report is None or self.qc_combo.count() == 0:
            return
        compound = self.qc_combo.currentData()
        rows = [r for r in self.result.report.istd if r.compound == compound]
        if not rows:
            return
        values = np.array([r.volume for r in rows], dtype=float)
        x = np.arange(len(values))
        self.qc_plot.plot(x, values, pen=None, symbol="o", symbolSize=7, symbolBrush=(126, 232, 250))
        mean, sd = float(values.mean()), float(values.std(ddof=1)) if len(values) > 1 else 0.0
        for k, colour in ((0, (150, 150, 150)), (2, (224, 160, 32)), (3, (194, 65, 12))):
            if k == 0:
                self.qc_plot.addLine(y=mean, pen=pg.mkPen(colour, style=QtCore.Qt.PenStyle.DashLine))
            else:
                self.qc_plot.addLine(y=mean + k * sd, pen=pg.mkPen(colour, style=QtCore.Qt.PenStyle.DotLine))
                self.qc_plot.addLine(y=mean - k * sd, pen=pg.mkPen(colour, style=QtCore.Qt.PenStyle.DotLine))
        flags = westgard_flags(values, mean=mean, sd=sd)
        if flags:
            names = [r.run for r in rows]
            text = "; ".join(f"{names[f.index]}: {f.rule}" for f in flags[:8])
            self.qc_flags.setText(f"{len(flags)} flag(s) -- {text}")
        else:
            self.qc_flags.setText("No Westgard rule violations.")

    # ------------------------------------------------------------------ PCA
    def _fill_pca(self, result) -> None:
        self.pca_plot.clear()
        ft = result.feature_table
        run_names = list(result.prepped)
        if not ft or len(run_names) < 2:
            return
        res = pca(ft, run_names, n_components=2)
        scatter = pg.ScatterPlotItem(pos=res.scores[:, :2], size=11, pen=pg.mkPen((126, 232, 250)),
                                     brush=pg.mkBrush(126, 232, 250, 150))
        self.pca_plot.addItem(scatter)
        for (x, y), name in zip(res.scores[:, :2], res.run_names):
            t = pg.TextItem(name, anchor=(0, 1), color=(200, 200, 210))
            t.setPos(float(x), float(y))
            self.pca_plot.addItem(t)
        ev = res.explained_variance
        self.pca_plot.setLabel("bottom", f"PC1 ({ev[0]:.0%})" if ev.size > 0 else "PC1")
        self.pca_plot.setLabel("left", f"PC2 ({ev[1]:.0%})" if ev.size > 1 else "PC2")
        scatter.sigClicked.connect(lambda _item, points, _ev: self._pca_point_clicked(points, res))

    def _pca_point_clicked(self, points, res) -> None:
        if not points:
            return
        i = int(points[0].index())
        name = res.run_names[i]
        p = self.result.prepped.get(name)
        if p is not None:
            self.navigateRequested.emit(name, float(p.run.times[0]), 0.0)

    # -------------------------------------------------------------- volcano
    def _fill_groups(self, result) -> None:
        for widget in (self.group_a_list, self.group_b_list):
            widget.clear()
            widget.addItems(list(result.prepped))

    def _draw_volcano(self) -> None:
        self.volcano_plot.clear()
        if self.result is None or not self.result.feature_table:
            return
        group_a = [i.text() for i in self.group_a_list.selectedItems()]
        group_b = [i.text() for i in self.group_b_list.selectedItems()]
        if len(group_a) < 2 or len(group_b) < 2:
            return
        points = volcano(self.result.feature_table, group_a, group_b)
        if not points:
            return
        x = np.array([p.log2_fold_change for p in points])
        y = -np.log10(np.clip(np.array([p.q_value for p in points]), 1e-300, None))
        sig = np.array([p.q_value < 0.05 for p in points])
        colours = np.where(sig, "#e08420", "#7e8ea6")
        scatter = pg.ScatterPlotItem(x=x, y=y, size=9,
                                     brush=[pg.mkBrush(c) for c in colours])
        self.volcano_plot.addItem(scatter)
        self.volcano_plot.addLine(y=-np.log10(0.05), pen=pg.mkPen((150, 150, 150),
                                                                   style=QtCore.Qt.PenStyle.DashLine))
        scatter.sigClicked.connect(lambda _item, pts, _ev: self._volcano_point_clicked(pts, points))

    def _volcano_point_clicked(self, pts, points) -> None:
        if not pts or self.result is None:
            return
        i = int(pts[0].index())
        vp = points[i]
        self.navigateRequested.emit(self.result.reference, vp.ref_t1, vp.ref_t2)


def _labelled(text: str, widget: QtWidgets.QWidget) -> QtWidgets.QWidget:
    box = QtWidgets.QWidget()
    lay = QtWidgets.QVBoxLayout(box)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.addWidget(QtWidgets.QLabel(text))
    lay.addWidget(widget)
    return box
