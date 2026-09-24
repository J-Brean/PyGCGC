"""🧩 Peak-clustering panel.

Settings on top (which runs, how fine the grouping is), the result
underneath: one row per cluster, and, when several runs were pooled, how
much of each run's peak volume falls in each cluster. The panel holds no
state of its own; `ClusterState` does.

A chromatogram can have a hundred clusters, more than colours can tell
apart, so only the largest few are coloured (Okabe-Ito, colour-blind safe)
and the rest share a grey.
"""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from ...core.cluster import CLUSTER_COLOURS, OTHER_COLOUR, ClusterSettings, cluster_colour  # noqa: F401  (re-exported)
from .classes import _Num

RUN_TABLE_CLUSTERS = 12       # cluster columns in the run table; the CSV has all of them


class ClustersPanel(QtWidgets.QWidget):
    settingsEdited = QtCore.Signal()
    clusterRequested = QtCore.Signal()
    overlayToggled = QtCore.Signal(bool)
    clusterSelected = QtCore.Signal(int)          # -1 when nothing is selected
    exportPeaksRequested = QtCore.Signal()
    exportRunsRequested = QtCore.Signal()
    figureRequested = QtCore.Signal()
    mergeRequested = QtCore.Signal()
    checkRequested = QtCore.Signal(str)           # "stability" or "scan"
    groupEdited = QtCore.Signal(str, str)         # run name, group label
    testGroupsRequested = QtCore.Signal()
    pcaFigureRequested = QtCore.Signal()
    nameEdited = QtCore.Signal(int, str)          # cluster or region, the name typed for it

    COLUMNS = ("Cluster", "Peaks", "Runs", "Volume %", "t1 (min)", "t2 (s)", "Top ions", "Silhouette", "Made of", "Name")
    NAME_COLUMN = 9

    def __init__(self, parent=None):
        super().__init__(parent)
        self._filling = False
        d = ClusterSettings()

        self.status = QtWidgets.QLabel("Select a run, then press Go.")
        self.status.setWordWrap(True)
        self.notes = QtWidgets.QLabel("")
        self.notes.setWordWrap(True)
        self.notes.setStyleSheet("color: palette(mid);")

        self.scope = QtWidgets.QComboBox()
        self.scope.addItem("This run", "run")
        self.scope.addItem("Every 70 eV run, pooled", "all")
        self.scope.setToolTip("Pooling clusters the peaks of all runs together, so a cluster means the same thing "
                              "in every run and the runs can be compared cluster by cluster.")
        self.unit_box = QtWidgets.QComboBox()
        self.unit_box.addItem("Detected peaks", "peaks")
        self.unit_box.addItem("Pixel regions (no peak detection)", "regions")
        self.unit_box.setToolTip("Peaks: group the blobs from peak detection by spectrum. Pixel regions: group every bright "
                                 "scan's spectrum and paint the result on the chromatogram, so unresolved humps, saturated "
                                 "peaks and overlaps are divided up too. Regions are found by k-means, so you choose how many.")
        self.n_regions = QtWidgets.QSpinBox()
        self.n_regions.setRange(2, 60)
        self.n_regions.setValue(12)
        self.n_regions.setToolTip("How many regions to divide the chromatogram into (pixel regions only)")
        self.method_box = QtWidgets.QComboBox()
        self.method_box.addItem("Spectra only (tree)", "tree")
        self.method_box.addItem("Spectra and retention (density)", "density")
        self.method_box.setToolTip("Tree: peaks with similar spectra group together wherever they elute, so clusters mean "
                                   "the same in every run and need no alignment. Density (HDBSCAN): peaks must also be "
                                   "close in retention, which finds bands and series but leaves more peaks unassigned; "
                                   "pooled runs should be aligned first.")
        self.retention = QtWidgets.QDoubleSpinBox()
        self.retention.setRange(0.0, 5.0)
        self.retention.setSingleStep(0.25)
        self.retention.setDecimals(2)
        self.retention.setValue(d.retention_weight)
        self.retention.setToolTip("Density only: how much retention counts against spectrum. 0 ignores it; 120 s of t1 or "
                                  "1 s of t2 counts as one unit of spectral distance at weight 1.")
        self.threshold = QtWidgets.QDoubleSpinBox()
        self.threshold.setRange(0.05, 0.80)
        self.threshold.setSingleStep(0.05)
        self.threshold.setDecimals(2)
        self.threshold.setValue(d.distance_threshold)
        self.threshold.setKeyboardTracking(False)
        self.threshold.setToolTip("Spectral distance (1 minus the cosine similarity) at which groups stop being merged. "
                                  "Lower gives more, tighter clusters (isomers and carbon numbers apart); higher gives "
                                  "fewer, broader ones (compound classes).")
        self.fixed_box = QtWidgets.QCheckBox("Cut at a set number of clusters")
        self.fixed_box.setToolTip("Ignore the distance and cut the tree into this many groups. Groups smaller than the "
                                  "minimum size are dropped, so you may get fewer.")
        self.fixed_n = QtWidgets.QSpinBox()
        self.fixed_n.setRange(2, 500)
        self.fixed_n.setValue(20)
        self.fixed_n.setEnabled(False)
        self.min_size = QtWidgets.QSpinBox()
        self.min_size.setRange(1, 100)
        self.min_size.setValue(d.min_cluster_size)
        self.min_size.setToolTip("A group with fewer peaks than this is not called a cluster.")
        self.siloxane_box = QtWidgets.QCheckBox("Leave out siloxane ions")
        self.siloxane_box.setToolTip("Drops the polydimethylsiloxane ions (73, 147, 207, 281 …) from the comparison so "
                                     "column bleed and septum peaks do not form a cluster of their own.")
        self.soft_box = QtWidgets.QCheckBox("Add the 14 eV spectrum when a partner run is loaded")
        self.soft_box.setToolTip("Peaks are compared on their 70 eV fragmentation pattern and, at half the weight, their "
                                 "14 eV spectrum, where the molecular ion survives. It can separate peaks whose 70 eV "
                                 "spectra look alike. Needs the 14 eV file of the same acquisition loaded alongside "
                                 "(File → Open files…), as for 🧬. Runs without a partner fall back to 70 eV alone, "
                                 "unless pooling, where every run needs one.")
        self.soft_box.toggled.connect(lambda _on: self.settingsEdited.emit())
        for w in (self.scope, self.unit_box, self.method_box):
            w.currentIndexChanged.connect(lambda _i: self._edited())
        for w in (self.threshold, self.fixed_n, self.min_size, self.n_regions, self.retention):
            w.valueChanged.connect(lambda _v: self._edited())
        self.fixed_box.toggled.connect(self.fixed_n.setEnabled)
        self.fixed_box.toggled.connect(lambda _on: self.settingsEdited.emit())
        self.siloxane_box.toggled.connect(lambda _on: self.settingsEdited.emit())

        self.fill_box = QtWidgets.QCheckBox("Fill peaks with their cluster colour in the 2D view")
        self.fill_box.setChecked(True)
        self.fill_box.setToolTip("The chromatogram coloured by cluster (also in the Layers tab). Select rows in the "
                                 "table (Ctrl or Shift for several) to colour just those; with none selected the "
                                 "seven largest are coloured and the rest are light grey.")
        self.fill_box.toggled.connect(self.overlayToggled)
        self.overlay_box = QtWidgets.QCheckBox("Outline peaks by cluster in the 2D view")
        self.overlay_box.setChecked(False)
        self.overlay_box.setToolTip("Coloured outlines instead of (or as well as) the fill (also in the Layers tab).")
        self.overlay_box.toggled.connect(self.overlayToggled)

        self.table = QtWidgets.QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.DoubleClicked
                                   | QtWidgets.QAbstractItemView.EditTrigger.EditKeyPressed)   # only the Name column is editable
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSortingEnabled(True)
        self.table.currentCellChanged.connect(self._row_changed)
        self.table.itemChanged.connect(self._table_item_changed)

        self.merge_btn = QtWidgets.QPushButton("Merge selected")
        self.merge_btn.setToolTip("Join the selected clusters (Ctrl or Shift for several) into one, for a family the tree split")
        self.merge_btn.clicked.connect(self.mergeRequested)
        self.stability_btn = QtWidgets.QPushButton("Check stability")
        self.stability_btn.setToolTip("Recluster random 80% samples of the peaks and see how well the groups reproduce "
                                 "(adjusted Rand index; 1 is identical every time)")
        self.stability_btn.clicked.connect(lambda: self.checkRequested.emit("stability"))
        self.scan_btn = QtWidgets.QPushButton("Scan distances")
        self.scan_btn.setToolTip("How many clusters, how much unassigned volume and what mean silhouette at merge distances 0.1 to 0.5")
        self.scan_btn.clicked.connect(lambda: self.checkRequested.emit("scan"))
        self.check_label = QtWidgets.QLabel("")
        self.check_label.setWordWrap(True)
        self.check_label.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)

        self.norm_box = QtWidgets.QComboBox()
        self.norm_box.addItem("Percent of each run's own total", "percent")
        self.norm_box.addItem("Percent, then PQN", "pqn")
        self.norm_box.setToolTip("How the run table is scaled before PCA, volcano and PERMANOVA. Percent of the run's own total "
                                 "assumes totals are comparable between runs. PQN (probabilistic quotient normalisation) "
                                 "also divides each run by its median quotient against the median run, so one large cluster "
                                 "changing does not make every other cluster appear to change the other way.")
        self.norm_box.currentIndexChanged.connect(lambda _i: self._edited())
        self.run_label = QtWidgets.QLabel("Share of each run's peak volume, %")
        self.run_table = QtWidgets.QTableWidget(0, 0)
        self.run_table.verticalHeader().setVisible(False)
        self.run_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.DoubleClicked
                                       | QtWidgets.QAbstractItemView.EditTrigger.EditKeyPressed)
        self.run_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.run_table.setMaximumHeight(170)
        self.run_table.itemChanged.connect(self._run_item_changed)
        self.run_label.setVisible(False)
        self.run_table.setVisible(False)
        test_btn = QtWidgets.QPushButton("Test groups")
        test_btn.setToolTip("Type a group name for each run in the Group column (double-click), then test whether the "
                            "groups differ overall (PERMANOVA) and, for two groups, which clusters differ")
        test_btn.clicked.connect(self.testGroupsRequested)
        pca_btn = QtWidgets.QPushButton("Run PCA figure…")
        pca_btn.setToolTip("PCA of the runs on their cluster shares, as a figure with a caption")
        pca_btn.clicked.connect(self.pcaFigureRequested)
        self.group_label = QtWidgets.QLabel("")
        self.group_label.setWordWrap(True)
        self.group_label.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.group_box = QtWidgets.QWidget()
        gl = QtWidgets.QVBoxLayout(self.group_box)
        gl.setContentsMargins(0, 0, 0, 0)
        grow = QtWidgets.QHBoxLayout()
        grow.addWidget(test_btn)
        grow.addWidget(pca_btn)
        gl.addLayout(grow)
        gl.addWidget(self.group_label)
        self.group_box.setVisible(False)

        figure_btn = QtWidgets.QPushButton("Save coloured chromatogram…")
        figure_btn.setToolTip("A publication figure of this run coloured by cluster (89 or 183 mm, 300 dpi) with a "
                              "caption file. Colours the selected clusters, or the seven largest.")
        figure_btn.clicked.connect(self.figureRequested)
        peaks_btn = QtWidgets.QPushButton("Export peaks CSV…")
        peaks_btn.clicked.connect(self.exportPeaksRequested)
        runs_btn = QtWidgets.QPushButton("Export run table CSV…")
        runs_btn.setToolTip("Every cluster's share of every run's peak volume (%): the input to PCA, volcano plots "
                            "and PERMANOVA")
        runs_btn.clicked.connect(self.exportRunsRequested)

        form = QtWidgets.QFormLayout()
        form.addRow("Cluster", self.scope)
        form.addRow("Group", self.unit_box)
        form.addRow("Regions", self.n_regions)
        form.addRow("Method", self.method_box)
        form.addRow("Retention weight", self.retention)
        form.addRow("Merge up to distance", self.threshold)
        form.addRow(self.fixed_box, self.fixed_n)
        form.addRow("Smallest cluster (peaks)", self.min_size)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(self.status)
        lay.addWidget(self.notes)
        lay.addLayout(form)
        lay.addWidget(self.siloxane_box)
        lay.addWidget(self.soft_box)
        lay.addWidget(self.fill_box)
        lay.addWidget(self.overlay_box)
        lay.addWidget(self.table, 1)
        tools = QtWidgets.QHBoxLayout()
        for b in (self.merge_btn, self.stability_btn, self.scan_btn):
            tools.addWidget(b)
        lay.addLayout(tools)
        lay.addWidget(self.check_label)
        lay.addWidget(self.run_label)
        lay.addWidget(self.norm_box)
        lay.addWidget(self.run_table)
        lay.addWidget(self.group_box)
        lay.addWidget(figure_btn)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(peaks_btn)
        row.addWidget(runs_btn)
        lay.addLayout(row)
        self._update_enabled()

    # -------------------------------------------------------------- setting
    def values(self) -> tuple:
        """(settings, leave out siloxane ions, scope, add the 14 eV spectrum) as shown. The ion list and the
        14 eV weight are added by the state."""
        base = ClusterSettings(distance_threshold=self.threshold.value(), min_cluster_size=self.min_size.value(),
                               n_clusters=self.fixed_n.value() if self.fixed_box.isChecked() else None,
                               method=self.method_box.currentData(), retention_weight=self.retention.value())
        return base, self.siloxane_box.isChecked(), self.scope.currentData(), self.soft_box.isChecked()

    def unit_values(self) -> tuple:
        """(peaks or regions, number of regions, how the run table is normalised)"""
        return self.unit_box.currentData(), self.n_regions.value(), self.norm_box.currentData()

    def _edited(self) -> None:
        self._update_enabled()
        self.settingsEdited.emit()

    def _update_enabled(self) -> None:
        """Grey out what does not apply: pixel regions ignore the peak settings, the tree ignores retention."""
        regions = self.unit_box.currentData() == "regions"
        density = self.method_box.currentData() == "density"
        self.n_regions.setEnabled(regions)
        for w in (self.method_box, self.fixed_box, self.min_size, self.soft_box, self.merge_btn, self.stability_btn,
                  self.scan_btn):
            w.setEnabled(not regions)
        self.retention.setEnabled(not regions and density)
        self.threshold.setEnabled(not regions and not density)
        self.fixed_n.setEnabled(not regions and not density and self.fixed_box.isChecked())
        self.table.horizontalHeaderItem(1).setText("Pixels" if regions else "Peaks")

    def set_values(self, settings: ClusterSettings, skip_siloxane: bool, scope: str, use_soft: bool = False,
                   unit: str = "peaks", n_regions: int = 12, normalise: str = "percent") -> None:
        """Show saved settings without announcing them as an edit."""
        widgets = (self.threshold, self.fixed_n, self.min_size, self.fixed_box, self.siloxane_box, self.scope,
                   self.soft_box, self.unit_box, self.n_regions, self.method_box, self.retention, self.norm_box)
        for w in widgets:
            w.blockSignals(True)
        self.norm_box.setCurrentIndex(max(self.norm_box.findData(normalise), 0))
        self.soft_box.setChecked(use_soft)
        self.unit_box.setCurrentIndex(max(self.unit_box.findData(unit), 0))
        self.n_regions.setValue(n_regions)
        self.method_box.setCurrentIndex(max(self.method_box.findData(settings.method), 0))
        self.retention.setValue(settings.retention_weight)
        self.threshold.setValue(settings.distance_threshold)
        self.min_size.setValue(settings.min_cluster_size)
        self.fixed_box.setChecked(settings.n_clusters is not None)
        self.fixed_n.setValue(settings.n_clusters or self.fixed_n.value())
        self.fixed_n.setEnabled(settings.n_clusters is not None)
        self.siloxane_box.setChecked(skip_siloxane)
        self.scope.setCurrentIndex(max(self.scope.findData(scope), 0))
        for w in widgets:
            w.blockSignals(False)
        self._update_enabled()

    def set_current(self, status: str, note: str, summary: list, stale: bool, runs=None, composition=None,
                    check: str = "", groups=None, names=None) -> None:
        """`summary` is the `ClusterSummary` list (or empty); `runs` is (run names, (features, percent matrix)) or None;
        `composition` maps a cluster to (main compound class, its share); `groups` maps a run to its group label;
        `names` maps a cluster to the name the user gave it."""
        self.status.setText(status + ("  Out of date: settings, blobs or runs changed. Press Go." if stale and summary else ""))
        self.notes.setText(note)
        self.check_label.setText(check)
        self._fill_table(summary, composition or {}, names or {})
        self._fill_runs(runs, groups or {})

    def set_group_text(self, text: str) -> None:
        self.group_label.setText(text)

    # ---------------------------------------------------------------- fills
    def _fill_table(self, summary: list, composition: dict, names: dict) -> None:
        self._filling = True
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(summary))
        for r, s in enumerate(summary):
            sil = "" if s.silhouette != s.silhouette else f"{s.silhouette:.2f}"       # nan: placed by centroid only
            label = _Num(f"■ {s.label}", s.label)
            label.setForeground(QtGui.QBrush(QtGui.QColor(*cluster_colour(s.label))))
            made = composition.get(s.label)
            cells = (label, _Num(str(s.n_peaks), s.n_peaks), _Num(str(s.n_runs), s.n_runs),
                     _Num(f"{s.percent:.1f}", s.percent), _Num(f"{s.t1_median / 60.0:.1f}", s.t1_median),
                     _Num(f"{s.t2_median:.2f}", s.t2_median), QtWidgets.QTableWidgetItem(" ".join(map(str, s.top_ions))),
                     _Num(sil, float("-inf") if not sil else s.silhouette),
                     QtWidgets.QTableWidgetItem(f"{made[0]} ({100 * made[1]:.0f}%)" if made else ""),
                     QtWidgets.QTableWidgetItem(names.get(s.label, "")))
            cells[6].setToolTip("m/z with the largest mean share of the ion current in this cluster, largest first")
            cells[8].setToolTip("The compound class most of this cluster's peaks were given by 🧬, in the current run. "
                                "Blank until that run is classified.")
            cells[self.NAME_COLUMN].setToolTip("Double-click to name this cluster. The name is used in the figure key and "
                                               "the peaks CSV. Names are lost when you cluster again.")
            for c, item in enumerate(cells):
                if c != self.NAME_COLUMN:
                    item.setFlags(item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, c, item)
        self.table.setSortingEnabled(True)
        self._filling = False

    def _fill_runs(self, runs, groups: dict) -> None:
        show = runs is not None and len(runs[0]) > 1
        self.run_label.setVisible(show)
        self.norm_box.setVisible(show)
        self.run_table.setVisible(show)
        self.group_box.setVisible(show)
        if not show:
            return
        names, (feats, matrix) = runs[0], runs[1]
        k = min(len(feats), RUN_TABLE_CLUSTERS)
        self._filling = True
        self.run_table.setColumnCount(k + 3)
        self.run_table.setHorizontalHeaderLabels(["Run", "Group"] + [f"c{f.id}" for f in feats[:k]] + ["unassigned"])
        self.run_table.setRowCount(len(names))
        read_only = ~QtCore.Qt.ItemFlag.ItemIsEditable
        for r, name in enumerate(names):
            cells = [QtWidgets.QTableWidgetItem(name), QtWidgets.QTableWidgetItem(groups.get(name, ""))]
            cells += [QtWidgets.QTableWidgetItem(f"{matrix[r, c]:.1f}") for c in range(k)]
            cells.append(QtWidgets.QTableWidgetItem(f"{max(100.0 - matrix[r].sum(), 0.0):.1f}"))
            for c, item in enumerate(cells):
                if c != 1:
                    item.setFlags(item.flags() & read_only)
                self.run_table.setItem(r, c, item)
        self._filling = False
        if len(feats) > k:
            self.run_label.setText(f"Share of each run's peak volume, %  (largest {k} of {len(feats)} clusters; the CSV has all)")
        else:
            self.run_label.setText("Share of each run's peak volume, %")

    def select_cluster(self, label: int) -> None:
        for r in range(self.table.rowCount()):
            item = self.table.item(r, 0)
            if item is not None and int(item.data(QtCore.Qt.ItemDataRole.UserRole)) == label:
                self.table.setCurrentCell(r, 0)
                return

    def selected_clusters(self) -> list:
        """Clusters of the selected rows, in ascending order."""
        labels = set()
        for r in {i.row() for i in self.table.selectedIndexes()}:
            item = self.table.item(r, 0)
            if item is not None:
                labels.add(int(item.data(QtCore.Qt.ItemDataRole.UserRole)))
        return sorted(labels)

    def set_colours(self, colours: dict) -> None:
        """Paint each row's swatch with the colour it has in the 2D view (grey if it is not coloured there)."""
        for r in range(self.table.rowCount()):
            item = self.table.item(r, 0)
            if item is not None:
                k = int(item.data(QtCore.Qt.ItemDataRole.UserRole))
                item.setForeground(QtGui.QBrush(QtGui.QColor(*colours.get(k, OTHER_COLOUR))))

    def selected_cluster(self) -> int:
        item = self.table.item(self.table.currentRow(), 0) if self.table.currentRow() >= 0 else None
        return -1 if item is None else int(item.data(QtCore.Qt.ItemDataRole.UserRole))

    def _row_changed(self, row: int, *_a) -> None:
        if not self._filling:
            self.clusterSelected.emit(self.selected_cluster() if row >= 0 else -1)

    def _table_item_changed(self, item) -> None:
        if not self._filling and item.column() == self.NAME_COLUMN:
            label = self.table.item(item.row(), 0).data(QtCore.Qt.ItemDataRole.UserRole)
            self.nameEdited.emit(int(label), item.text())

    def _run_item_changed(self, item) -> None:
        if not self._filling and item.column() == 1:
            self.groupEdited.emit(self.run_table.item(item.row(), 0).text(), item.text())
