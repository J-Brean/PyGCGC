"""Peak-detection panel for both legacy borders and MS-aware regions."""

from __future__ import annotations

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from ..blob_state import LEGACY, MS_AWARE
from ..fitting_help import fitting_help_button

LEGACY_COLUMNS = ["#", "t1 (min)", "t2 (s)", "volume", "S/N", "width t1 (s)",
                  "width t2 (s)", "n mod", "pixels", "co-elution?", "🫧 siloxane?"]
REGION_COLUMNS = ["#", "t1 (min)", "t2 (s)", "candidate centres", "accepted components",
                  "supporting ions", "strongest evidence", "t1 span (s)", "t2 span (s)"]
COMPONENT_COLUMNS = ["#", "t1 (min)", "t2 (s)", "estimated intensity", "relative size",
                     "diagnostic ions", "confidence", "group"]


class BlobsPanel(QtWidgets.QWidget):
    findRequested = QtCore.Signal()
    previewRequested = QtCore.Signal()
    methodChanged = QtCore.Signal(str)
    resolveRequested = QtCore.Signal(int)
    addRequested = QtCore.Signal()
    deleteRequested = QtCore.Signal(list)
    mergeRequested = QtCore.Signal(list)
    rowActivated = QtCore.Signal(int)
    selectionChanged = QtCore.Signal(list)
    exportRequested = QtCore.Signal()
    accountingRequested = QtCore.Signal()
    siloxaneRequested = QtCore.Signal()
    hydrocarbonRequested = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.blobs = None
        self.candidates = None
        self.regions = None
        self.component_set = None
        self.method = MS_AWARE
        self.current_run: str | None = None
        self._updating = False
        self._row_regions: list[int] = []
        self._row_components: list = []
        self._row_groups: list = []

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
        self.preview_btn = QtWidgets.QPushButton("Preview fitting…")
        self.preview_btn.setToolTip("Try fitting the selected region, visible interval or current run before processing all runs")
        self.preview_btn.clicked.connect(self.previewRequested)
        self.fitting_help_btn = fitting_help_button(self)
        method_row.addWidget(self.fitting_help_btn)
        self.resolve_btn = QtWidgets.QPushButton("Inspect/refit selected blob…")
        self.resolve_btn.clicked.connect(self._resolve_selected)
        self.add_btn = QtWidgets.QPushButton("➕ Add")
        self.add_btn.clicked.connect(self.addRequested)
        self.merge_btn = QtWidgets.QPushButton("Merge  (G)")
        self.merge_btn.clicked.connect(self._merge_selected)
        self.delete_btn = QtWidgets.QPushButton("Delete")
        self.delete_btn.clicked.connect(self._delete_selected)
        self.export_btn = QtWidgets.QPushButton("Export peaks CSV…")
        self.export_btn.clicked.connect(self.exportRequested)
        self.account_btn = QtWidgets.QPushButton("📊 Fit and residual map…")
        self.account_btn.setToolTip("Opens a window with the observed chromatogram, the sum of the fitted blobs and "
                                    "what is left over, and reports how much of the signal the blobs explain.")
        self.account_btn.clicked.connect(self.accountingRequested)
        self.accounting = QtWidgets.QLabel("")
        self.accounting.setWordWrap(True)
        self.accounting.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        # Tools that work on the peaks found here; each opens its own sub-menu in the Tool tab.
        self.siloxane_btn = QtWidgets.QPushButton("🫧 Siloxane filterer…")
        self.siloxane_btn.setToolTip("Score every peak against a siloxane reference spectrum and filter out the ones "
                                     "above a certainty you set")
        self.siloxane_btn.clicked.connect(self.siloxaneRequested)
        self.hydrocarbon_btn = QtWidgets.QPushButton("⛓️ Hydrocarbon identifier…")
        self.hydrocarbon_btn.setToolTip("Anchor the n-alkane ladder and assign every peak to the hydrocarbon it is")
        self.hydrocarbon_btn.clicked.connect(self.hydrocarbonRequested)
        more = QtWidgets.QToolButton()
        more.setText("Peak tools…")
        more.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QtWidgets.QMenu(more)
        self.tools_menu = menu
        menu_buttons = (self.resolve_btn, self.add_btn, self.merge_btn, self.delete_btn, self.export_btn,
                        self.account_btn, self.siloxane_btn, self.hydrocarbon_btn)
        actions = []
        for button in menu_buttons:
            button.setParent(self)
            button.hide()
            action = menu.addAction(button.text(), button.click)
            action.setToolTip(button.toolTip())
            actions.append((action, button))
        menu.aboutToShow.connect(lambda: [a.setEnabled(b.isEnabled()) for a, b in actions])
        more.setMenu(menu)

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
        for widget in (self.find_btn, self.preview_btn, more):
            buttons.addWidget(widget)
        buttons.addStretch(1)
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addLayout(method_row)
        layout.addWidget(self.explanation)
        layout.addLayout(buttons)
        layout.addWidget(self.table, 1)
        layout.addWidget(self.accounting)
        layout.addWidget(self.status)
        self.set_method(MS_AWARE)

    def set_method(self, method: str) -> None:
        self.method = method
        self.method_box.blockSignals(True)
        index = self.method_box.findData(method)
        if index >= 0:
            self.method_box.setCurrentIndex(index)
        self.method_box.blockSignals(False)
        ms_aware = method == MS_AWARE
        self.fitting_help_btn.setVisible(ms_aware)
        self.preview_btn.setVisible(ms_aware)
        self.explanation.setText(
            "Fits overlapping tilted Gaussian peaks using both position and mass spectrum. Each row and "
            "colour is one proposed peak; unresolved estimates remain visible with a warning." if ms_aware else
            "Finds peaks on one combined intensity surface and assigns each pixel to one hard-edged blob. "
            "Kept for compatibility with existing siloxane, ladder and library-search tools.")
        columns = COMPONENT_COLUMNS if ms_aware else LEGACY_COLUMNS
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
        self.table.setColumnCount(len(COMPONENT_COLUMNS))
        self.table.setHorizontalHeaderLabels(COMPONENT_COLUMNS)
        self._row_regions = []
        self._row_components = []
        self._row_groups = []
        self.table.setRowCount(0)
        if running:
            self.status.setText(f"Finding overlapping peaks on {self.current_run}…")
        elif regions is None:
            self.status.setText("Press Find peaks to detect overlapping blobs automatically."
                                if self.current_run else "")
        elif stale:
            self.status.setText("Peak results are out of date. Press Find peaks to refresh.")
        else:
            self.status.setText("No overlapping blobs were fitted. Enable automatic separation "
                                "in detection settings and press Find peaks.")

    def _fill_components(self, stale: bool, running: bool) -> None:
        self.table.setColumnCount(len(COMPONENT_COLUMNS))
        self.table.setHorizontalHeaderLabels(COMPONENT_COLUMNS)
        rows = [(group, component) for group in self.component_set.groups
                for component in group.components]
        self._row_regions = [group.region_index for group, _component in rows]
        self._row_components = [component for _group, component in rows]
        self._row_groups = [group for group, _component in rows]
        self.table.setRowCount(len(rows))
        for row, (group, component) in enumerate(rows):
            order = np.argsort(component.spectrum)[::-1]
            ions = [f"{component.mz[index]:g}" for index in order[:3]
                    if component.spectrum[index] > 0]
            confidence = ("candidate proposal" if group.warning.startswith("rapid candidate proposals")
                          else component.confidence)
            if group.unresolved:
                confidence += " / unresolved group"
            values = (str(row), f"{component.t1 / 60:.3f}", f"{component.t2:.3f}",
                      f"{component.volume:,.0f}",
                      f"{100 * component.relative_to_largest:.0f}%", ", ".join(ions),
                      confidence, str(group.region_index))
            self._set_row(row, values)
            confidence_item = self.table.item(row, 6)
            if component.confidence == "low":
                confidence_item.setForeground(QtGui.QColor("#ff9f43"))
            detail = group.warning or "No automatic quality warning"
            confidence_item.setToolTip(detail)
            explained_text = (f"{100 * group.explained_fraction:.1f}%"
                              if np.isfinite(group.explained_fraction) else "not measured")
            bic_text = (f"{group.bic_gap:.1f}"
                        if np.isfinite(group.bic_gap) else "not applicable")
            spread_text = (f"{group.restart_position_spread:.2f} pixels"
                           if np.isfinite(group.restart_position_spread) else "not checked")
            self.table.item(row, 7).setToolTip(
                f"Explained fraction: {explained_text}\n"
                f"BIC gap: {bic_text}\n"
                f"Restart position spread: {spread_text}\n"
                f"{detail}")
        if running:
            self.status.setText(f"Separating overlapping components on {self.current_run}…")
        else:
            unresolved = sum(group.unresolved for group in self.component_set.groups)
            stale_note = " Results are out of date." if stale else ""
            proposals = sum(group.warning.startswith("rapid candidate proposals")
                            for group in self.component_set.groups)
            if proposals:
                self.status.setText(
                    f"{len(rows)} rapid peak proposals in {len(self.component_set.groups)} groups. "
                    "Open a group to jointly separate its overlapping spectra."
                    f"{stale_note}")
            else:
                self.status.setText(
                    f"{len(rows)} overlapping blobs found automatically; "
                    f"{unresolved} unresolved group{'s' if unresolved != 1 else ''}.{stale_note}")

    def _fill_blobs(self, stale, running, siloxane_scores, siloxane_threshold) -> None:
        blobs = self.blobs
        self._row_regions = []
        self._row_components = []
        self._row_groups = []
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

    def toggle_row(self, index: int) -> None:
        """Add a row to the selection, or drop it if already selected (Ctrl+click on the image)."""
        if not 0 <= index < self.table.rowCount():
            return
        self._updating = True
        model = self.table.selectionModel()
        model.select(self.table.model().index(index, 0),
                     QtCore.QItemSelectionModel.SelectionFlag.Toggle
                     | QtCore.QItemSelectionModel.SelectionFlag.Rows)
        self.table.scrollToItem(self.table.item(index, 0))
        self._updating = False
        self._refresh_buttons()

    def show_accounting(self, text: str) -> None:
        self.accounting.setText(text)

    def selected_indices(self) -> list[int]:
        return sorted({index.row() for index in self.table.selectedIndexes()})

    def region_for_row(self, row: int) -> int | None:
        return self._row_regions[row] if 0 <= row < len(self._row_regions) else None

    def component_for_row(self, row: int):
        return self._row_components[row] if 0 <= row < len(self._row_components) else None

    def group_for_row(self, row: int):
        return self._row_groups[row] if 0 <= row < len(self._row_groups) else None

    def selected_component_ids(self) -> set[str]:
        """Ids of the selected fitted blobs: each one is its own thing, even inside one fitted group."""
        return {self._row_components[row].id for row in self.selected_indices()
                if 0 <= row < len(self._row_components) and self._row_components[row] is not None}

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
        self.resolve_btn.hide()
        self.add_btn.setEnabled(legacy and self.current_run is not None)
        self.add_btn.hide()
        self.merge_btn.setEnabled(legacy and selected >= 2)
        self.merge_btn.hide()
        self.delete_btn.setEnabled(legacy and selected >= 1)
        self.delete_btn.hide()
        has_components = (self.component_set is not None
                          and bool(self.component_set.components))
        self.export_btn.hide()
        self.export_btn.setEnabled(not legacy and has_components)
        self.account_btn.hide()
        self.account_btn.setEnabled(not legacy and has_components)
        self.accounting.setVisible(not legacy)

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
