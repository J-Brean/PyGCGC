"""Peak-detection panel for both legacy borders and MS-aware regions."""

from __future__ import annotations

import numpy as np
from PySide6 import QtCore, QtWidgets

from ..blob_state import LEGACY, MS_AWARE

LEGACY_COLUMNS = ["#", "t1 (min)", "t2 (s)", "volume", "S/N", "width t1 (s)",
                  "width t2 (s)", "n mod", "pixels", "co-elution?", "🫧 siloxane?"]
REGION_COLUMNS = ["#", "t1 (min)", "t2 (s)", "candidate centres", "accepted components",
                  "supporting ions", "strongest evidence", "t1 span (s)", "t2 span (s)"]
COMPONENT_COLUMNS = ["#", "t1 (min)", "t2 (s)", "estimated intensity", "relative size",
                     "diagnostic ions", "confidence", "group"]


class BlobsPanel(QtWidgets.QWidget):
    findRequested = QtCore.Signal()
    methodChanged = QtCore.Signal(str)
    resolveRequested = QtCore.Signal(int)
    addRequested = QtCore.Signal()
    deleteRequested = QtCore.Signal(list)
    mergeRequested = QtCore.Signal(list)
    rowActivated = QtCore.Signal(int)
    selectionChanged = QtCore.Signal(list)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.blobs = None
        self.candidates = None
        self.regions = None
        self.component_set = None
        self.method = LEGACY
        self.current_run: str | None = None
        self._updating = False
        self._row_regions: list[int] = []
        self._row_components: list = []

        self.method_box = QtWidgets.QComboBox()
        self.method_box.addItem("MS-aware separation", MS_AWARE)
        self.method_box.addItem("Legacy blob borders", LEGACY)
        self.method_box.setToolTip(
            "MS-aware separation keeps ion evidence separate and permits overlapping compounds. "
            "Legacy detection assigns every pixel to one hard-edged blob.")
        self.method_box.currentIndexChanged.connect(self._method_changed)
        method_row = QtWidgets.QHBoxLayout()
        method_row.addWidget(QtWidgets.QLabel("Detection method:"))
        method_row.addWidget(self.method_box, 1)

        self.explanation = QtWidgets.QLabel()
        self.explanation.setWordWrap(True)
        self.find_btn = QtWidgets.QPushButton("Find peaks")
        self.find_btn.clicked.connect(self.findRequested)
        self.resolve_btn = QtWidgets.QPushButton("Resolve selected region…")
        self.resolve_btn.clicked.connect(self._resolve_selected)
        self.add_btn = QtWidgets.QPushButton("➕ Add")
        self.add_btn.clicked.connect(self.addRequested)
        self.merge_btn = QtWidgets.QPushButton("Merge  (G)")
        self.merge_btn.clicked.connect(self._merge_selected)
        self.delete_btn = QtWidgets.QPushButton("Delete")
        self.delete_btn.clicked.connect(self._delete_selected)

        self.table = QtWidgets.QTableWidget(0, len(LEGACY_COLUMNS))
        self.table.setHorizontalHeaderLabels(LEGACY_COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.itemSelectionChanged.connect(self._selection_changed)
        self.table.cellDoubleClicked.connect(lambda r, _c: self.rowActivated.emit(r))
        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)

        buttons = QtWidgets.QHBoxLayout()
        for widget in (self.find_btn, self.resolve_btn, self.add_btn, self.merge_btn, self.delete_btn):
            buttons.addWidget(widget)
        buttons.addStretch(1)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addLayout(method_row)
        layout.addWidget(self.explanation)
        layout.addLayout(buttons)
        layout.addWidget(self.table, 1)
        layout.addWidget(self.status)
        self.set_method(LEGACY)

    def set_method(self, method: str) -> None:
        self.method = method
        self.method_box.blockSignals(True)
        index = self.method_box.findData(method)
        if index >= 0:
            self.method_box.setCurrentIndex(index)
        self.method_box.blockSignals(False)
        ms_aware = method == MS_AWARE
        self.explanation.setText(
            "Finds peak centres from ions that independently agree, then groups nearby centres into "
            "shared analysis regions. The rectangles are work areas, not claims that one compound owns "
            "every enclosed pixel." if ms_aware else
            "Finds peaks on one combined intensity surface and assigns each pixel to one hard-edged blob. "
            "Kept for compatibility with existing siloxane, ladder and library-search tools.")
        columns = REGION_COLUMNS if ms_aware else LEGACY_COLUMNS
        self.table.setColumnCount(len(columns))
        self.table.setHorizontalHeaderLabels(columns)
        self._refresh_buttons()

    def set_current_run(self, name: str | None, blobs, stale: bool, running: bool,
                        siloxane_scores=None, siloxane_threshold: float = 0.0,
                        *, method: str = LEGACY, candidates=None, regions=None,
                        component_set=None) -> None:
        self.current_run, self.blobs = name, blobs
        self.candidates, self.regions = candidates, regions
        self.component_set = component_set
        self._updating = True
        self.set_method(method)
        if method == MS_AWARE:
            self._fill_regions(stale, running)
        else:
            self._fill_blobs(stale, running, siloxane_scores, siloxane_threshold)
        self._updating = False
        self._refresh_buttons()

    def _fill_regions(self, stale: bool, running: bool) -> None:
        candidates, regions = self.candidates, self.regions
        if self.component_set is not None and self.component_set.components:
            self._fill_components(stale, running)
            return
        self.table.setColumnCount(len(REGION_COLUMNS))
        self.table.setHorizontalHeaderLabels(REGION_COLUMNS)
        self._row_regions = list(range(0 if regions is None else len(regions)))
        self._row_components = [None] * len(self._row_regions)
        count = 0 if regions is None else len(regions)
        self.table.setRowCount(count)
        for row, region in enumerate(regions or []):
            indices = np.asarray(region.candidate_indices, dtype=int)
            strongest = int(indices[np.argmax(candidates.score[indices])])
            support = {int(mz) for i in indices for mz in candidates.support_ions[i]}
            accepted = (self.component_set.group_for_region(row)
                        if self.component_set is not None else None)
            values = (str(row), f"{candidates.t1[strongest] / 60:.3f}",
                      f"{candidates.t2[strongest]:.3f}", str(len(indices)),
                      str(len(accepted.components)) if accepted is not None else "",
                      str(len(support)), f"{candidates.score[strongest]:.1f}",
                      f"{(region.col1-region.col0)*candidates.period:.1f}",
                      f"{(region.row1-region.row0)*candidates.pitch:.3f}")
            self._set_row(row, values)
        if running:
            self.status.setText(f"Finding MS-aware peak regions on {self.current_run}…")
        elif regions is None:
            self.status.setText("No MS-aware detection has been run on this data yet." if self.current_run else "")
        elif stale:
            self.status.setText(f"{count} regions, out of date. Find peaks to refresh.")
        else:
            accepted = (len(self.component_set.components)
                        if self.component_set is not None else 0)
            self.status.setText(
                f"{count} analysis regions from {len(candidates)} candidate centres; "
                f"{accepted} accepted overlapping component{'s' if accepted != 1 else ''}. ")

    def _fill_components(self, stale: bool, running: bool) -> None:
        self.table.setColumnCount(len(COMPONENT_COLUMNS))
        self.table.setHorizontalHeaderLabels(COMPONENT_COLUMNS)
        rows = [(group, component) for group in self.component_set.groups
                for component in group.components]
        self._row_regions = [group.region_index for group, _component in rows]
        self._row_components = [component for _group, component in rows]
        self.table.setRowCount(len(rows))
        for row, (group, component) in enumerate(rows):
            order = np.argsort(component.spectrum)[::-1]
            ions = [f"{component.mz[index]:g}" for index in order[:3]
                    if component.spectrum[index] > 0]
            confidence = component.confidence
            if group.unresolved:
                confidence += " / unresolved group"
            values = (str(row), f"{component.t1 / 60:.3f}", f"{component.t2:.3f}",
                      f"{component.volume:,.0f}",
                      f"{100 * component.relative_to_largest:.0f}%", ", ".join(ions),
                      confidence, str(group.region_index))
            self._set_row(row, values)
        if running:
            self.status.setText(f"Separating overlapping components on {self.current_run}…")
        else:
            unresolved = sum(group.unresolved for group in self.component_set.groups)
            stale_note = " Results are out of date." if stale else ""
            self.status.setText(
                f"{len(rows)} overlapping components in {len(self.component_set.groups)} groups; "
                f"{unresolved} unresolved group{'s' if unresolved != 1 else ''}.{stale_note}")

    def _fill_blobs(self, stale, running, siloxane_scores, siloxane_threshold) -> None:
        blobs = self.blobs
        self._row_regions = []
        self._row_components = []
        count = 0 if blobs is None else len(blobs)
        self.table.setRowCount(count)
        for i in range(count):
            has_width = blobs.width_t1.shape[0] == count
            has_co = blobs.co_elution.shape[0] == count
            has_sil = siloxane_scores is not None and len(siloxane_scores) == count
            values = (str(i), f"{blobs.t1[i] / 60:.3f}", f"{blobs.t2[i]:.3f}",
                      f"{blobs.volume[i]:,.0f}", f"{blobs.snr[i]:.1f}",
                      f"{blobs.width_t1[i]:.2f}" if has_width else "",
                      f"{blobs.width_t2[i]:.3f}" if has_width else "",
                      str(int(blobs.n_modulations[i])), str(int(blobs.n_pixels[i])),
                      ("yes" if blobs.co_elution[i] else "") if has_co else "",
                      (f"{'yes ' if siloxane_scores[i] >= siloxane_threshold else ''}"
                       f"{100 * siloxane_scores[i]:.0f}%") if has_sil else "")
            self._set_row(i, values)
        if running:
            self.status.setText(f"Finding legacy blobs on {self.current_run}…")
        elif blobs is None:
            self.status.setText("No legacy blobs found on this run yet." if self.current_run else "")
        elif stale:
            self.status.setText(f"{count} blobs, out of date. Find peaks to refresh.")
        else:
            self.status.setText(f"{count} legacy blobs.")

    def _set_row(self, row, values) -> None:
        for column, value in enumerate(values):
            item = self.table.item(row, column)
            if item is None:
                item = QtWidgets.QTableWidgetItem()
                self.table.setItem(row, column, item)
            item.setText(value)

    def select_row(self, index: int) -> None:
        self._updating = True
        self.table.clearSelection()
        if 0 <= index < self.table.rowCount():
            self.table.selectRow(index)
            self.table.scrollToItem(self.table.item(index, 0))
        self._updating = False
        self._refresh_buttons()

    def selected_indices(self) -> list[int]:
        return sorted({index.row() for index in self.table.selectedIndexes()})

    def region_for_row(self, row: int) -> int | None:
        return self._row_regions[row] if 0 <= row < len(self._row_regions) else None

    def component_for_row(self, row: int):
        return self._row_components[row] if 0 <= row < len(self._row_components) else None

    def selected_region_indices(self) -> list[int]:
        return sorted({self._row_regions[row] for row in self.selected_indices()
                       if 0 <= row < len(self._row_regions)})

    def select_region(self, region_index: int) -> None:
        row = next((row for row, value in enumerate(self._row_regions)
                    if value == region_index), -1)
        self.select_row(row)

    def select_component(self, component_id: str) -> None:
        row = next((row for row, component in enumerate(self._row_components)
                    if component is not None and component.id == component_id), -1)
        self.select_row(row)

    def _method_changed(self) -> None:
        if not self._updating:
            method = self.method_box.currentData()
            self.set_method(method)
            self.methodChanged.emit(method)

    def _selection_changed(self) -> None:
        self._refresh_buttons()
        if not self._updating:
            self.selectionChanged.emit(self.selected_indices())

    def _refresh_buttons(self) -> None:
        selected = len(self.selected_indices())
        legacy = self.method == LEGACY
        self.resolve_btn.setEnabled(not legacy and selected == 1)
        self.resolve_btn.setVisible(not legacy)
        self.add_btn.setEnabled(legacy and self.current_run is not None)
        self.add_btn.setVisible(legacy)
        self.merge_btn.setEnabled(legacy and selected >= 2)
        self.merge_btn.setVisible(legacy)
        self.delete_btn.setEnabled(legacy and selected >= 1)
        self.delete_btn.setVisible(legacy)

    def _resolve_selected(self) -> None:
        rows = self.selected_indices()
        if len(rows) == 1:
            region = self.region_for_row(rows[0])
            if region is not None:
                self.resolveRequested.emit(region)

    def _delete_selected(self) -> None:
        rows = self.selected_indices()
        if rows:
            self.deleteRequested.emit(rows)

    def _merge_selected(self) -> None:
        rows = self.selected_indices()
        if len(rows) >= 2:
            self.mergeRequested.emit(rows)
