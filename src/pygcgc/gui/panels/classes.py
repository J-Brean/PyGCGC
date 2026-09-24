"""🧬 Compound-class panel.

Settings on top (which rules, which classes are on, whether to use the 14 eV
channel), the result underneath: a volume summary by class, a table of every
blob, and, for the selected blob, the plain-text reason it got its class.
The panel holds no state of its own; `ClassState` does.
"""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

GROUP_COLOURS = {
    "hydrocarbon": (0, 114, 178), "aromatic": (230, 159, 0), "oxygen": (0, 158, 115),
    "nitrogen": (213, 94, 0), "sulfur": (204, 121, 167), "contaminant": (153, 153, 153),
}
DEFAULT_COLOUR = (200, 200, 200)
TIER_ORDER = {"confident": 3, "probable": 2, "possible": 1, "unclassified": 0}


def group_colour(group: str | None) -> tuple:
    return GROUP_COLOURS.get(group or "", DEFAULT_COLOUR)


class _Num(QtWidgets.QTableWidgetItem):
    """A cell that displays text but sorts on a number."""

    def __init__(self, text: str, key: float):
        super().__init__(text)
        self.setData(QtCore.Qt.ItemDataRole.UserRole, float(key))

    def __lt__(self, other) -> bool:
        return self.data(QtCore.Qt.ItemDataRole.UserRole) < other.data(QtCore.Qt.ItemDataRole.UserRole)


class ClassesPanel(QtWidgets.QWidget):
    classifyRequested = QtCore.Signal()
    softToggled = QtCore.Signal(bool)
    loadRulesRequested = QtCore.Signal()
    copyRulesRequested = QtCore.Signal()
    reloadRulesRequested = QtCore.Signal()
    shippedRulesRequested = QtCore.Signal()
    classToggled = QtCore.Signal(str, bool)
    overlayToggled = QtCore.Signal(bool)
    blobSelected = QtCore.Signal(int)
    blobActivated = QtCore.Signal(int)
    exportRequested = QtCore.Signal()

    COLUMNS = ("#", "t1 (min)", "t2 (s)", "Class", "Tier", "Score", "Runner-up", "Held back by")

    def __init__(self, parent=None):
        super().__init__(parent)
        self._filling = False
        self._result = None

        self.status = QtWidgets.QLabel("Select a run, then press Go.")
        self.status.setWordWrap(True)
        self.notes = QtWidgets.QLabel("")
        self.notes.setWordWrap(True)
        self.notes.setStyleSheet("color: palette(mid);")

        self.soft_box = QtWidgets.QCheckBox("Use the 14 eV channel when a partner run is loaded")
        self.soft_box.setChecked(True)
        self.soft_box.setToolTip("A tandem acquisition stores its 70 eV and 14 eV scans as two files. At 14 eV the "
                                 "molecular ion survives, which gives the nitrogen rule and the homologous series "
                                 "(M mod 14) as extra evidence. Open the 14 eV file with File → Open files… "
                                 "alongside the 70 eV one.")
        self.soft_box.toggled.connect(self.softToggled)
        self.soft_status = QtWidgets.QLabel("")
        self.soft_status.setWordWrap(True)

        self.rules_label = QtWidgets.QLabel("Rules: shipped")
        self.rules_label.setWordWrap(True)
        load = QtWidgets.QPushButton("Load rules…")
        load.setToolTip("Use your own rule file (TOML, same format as the shipped one)")
        load.clicked.connect(self.loadRulesRequested)
        copy = QtWidgets.QPushButton("Save an editable copy…")
        copy.setToolTip("Write the shipped rules where you can edit them, and use that copy")
        copy.clicked.connect(self.copyRulesRequested)
        reload_ = QtWidgets.QPushButton("Reload")
        reload_.setToolTip("Read the rule file again after editing it")
        reload_.clicked.connect(self.reloadRulesRequested)
        shipped = QtWidgets.QPushButton("Shipped")
        shipped.setToolTip("Go back to the shipped rules")
        shipped.clicked.connect(self.shippedRulesRequested)

        self.class_list = QtWidgets.QListWidget()
        self.class_list.setMaximumHeight(150)
        self.class_list.setToolTip("Untick a class to leave it out. Hover for what it looks for.")
        self.class_list.itemChanged.connect(self._class_toggled)

        self.overlay_box = QtWidgets.QCheckBox("Colour blob outlines by group in the 2D view")
        self.overlay_box.setChecked(True)
        self.overlay_box.toggled.connect(self.overlayToggled)
        self.legend = QtWidgets.QLabel("")
        self.legend.setTextFormat(QtCore.Qt.TextFormat.RichText)

        self.summary = QtWidgets.QTableWidget(0, 3)
        self.summary.setHorizontalHeaderLabels(["Class", "Blobs", "Volume %"])
        self.summary.verticalHeader().setVisible(False)
        self.summary.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.summary.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.summary.setMaximumHeight(170)

        self.hide_un = QtWidgets.QCheckBox("Hide unclassified")
        self.hide_un.toggled.connect(lambda _on: self._fill_table())
        export = QtWidgets.QPushButton("Export CSV…")
        export.clicked.connect(self.exportRequested)
        self.table = QtWidgets.QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSortingEnabled(True)
        self.table.currentCellChanged.connect(self._row_changed)
        self.table.cellDoubleClicked.connect(self._row_activated)

        self.explain = QtWidgets.QPlainTextEdit()
        self.explain.setReadOnly(True)
        self.explain.setMaximumHeight(150)
        self.explain.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont))
        self.explain.setPlaceholderText("Select a blob to see every criterion, its value and its score.")

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(self.status)
        lay.addWidget(self.notes)
        lay.addWidget(self.soft_box)
        lay.addWidget(self.soft_status)
        lay.addWidget(self.rules_label)
        row = QtWidgets.QHBoxLayout()
        for b in (load, copy, reload_, shipped):
            row.addWidget(b)
        lay.addLayout(row)
        lay.addWidget(QtWidgets.QLabel("Classes"))
        lay.addWidget(self.class_list)
        lay.addWidget(self.overlay_box)
        lay.addWidget(self.legend)
        lay.addWidget(self.summary)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Blobs"))
        row.addStretch(1)
        row.addWidget(self.hide_un)
        row.addWidget(export)
        lay.addLayout(row)
        lay.addWidget(self.table, 1)
        lay.addWidget(QtWidgets.QLabel("Why this class"))
        lay.addWidget(self.explain)

    # -------------------------------------------------------------- setting
    def set_rules(self, rules, enabled: dict, label: str) -> None:
        """Show the class list (tick = on) for `rules` with the user's `enabled` overrides."""
        self._filling = True
        self.rules_label.setText(f"Rules: {label}")
        self.class_list.clear()
        for c in rules.classes:
            item = QtWidgets.QListWidgetItem(f"{c.name}   [{c.group}]")
            item.setFlags(item.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(QtCore.Qt.CheckState.Checked if enabled.get(c.name, c.enabled)
                               else QtCore.Qt.CheckState.Unchecked)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, c.name)
            item.setForeground(QtGui.QBrush(QtGui.QColor(*group_colour(c.group))))
            item.setToolTip(c.note or c.name)
            self.class_list.addItem(item)
        self._filling = False
        groups = sorted({c.group for c in rules.classes})
        self.legend.setText("  ".join(
            f'<span style="color: rgb{group_colour(g)}">■</span> {g}' for g in groups))

    def set_soft(self, checked: bool, text: str) -> None:
        self.soft_box.blockSignals(True)
        self.soft_box.setChecked(checked)
        self.soft_box.blockSignals(False)
        self.soft_status.setText(text)

    def set_current(self, status: str, note: str, result, blobs, volume, stale: bool) -> None:
        """The current run's result (or None). `volume` is one number per blob, for the summary."""
        self._result, self._blobs, self._volume = result, blobs, volume
        self.status.setText(status + ("  Out of date: settings, blobs or the ladder changed. Press Go." if stale and result else ""))
        self.notes.setText(note)
        self._fill_summary()
        self._fill_table()

    # ---------------------------------------------------------------- fills
    def _fill_summary(self) -> None:
        res = self._result
        rows = res.summary(self._volume) if res is not None else []
        self.summary.setRowCount(len(rows))
        groups = {c.name: c.group for c in res.classes} if res is not None else {}
        for r, (name, n, _v, pct) in enumerate(rows):
            for c, text in enumerate((name, str(n), f"{pct:.1f}")):
                item = QtWidgets.QTableWidgetItem(text)
                if c == 0:
                    item.setForeground(QtGui.QBrush(QtGui.QColor(*group_colour(groups.get(name)))))
                self.summary.setItem(r, c, item)

    def _fill_table(self) -> None:
        res, blobs = self._result, getattr(self, "_blobs", None)
        self._filling = True
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        self.explain.clear()
        if res is None or blobs is None:
            self._filling = False
            self.table.setSortingEnabled(True)
            return
        show = [i for i, c in enumerate(res.calls) if not (self.hide_un.isChecked() and c.cls is None)]
        self.table.setRowCount(len(show))
        for r, i in enumerate(show):
            c = res.calls[i]
            cells = (_Num(str(i), i), _Num(f"{blobs.t1[i] / 60.0:.2f}", blobs.t1[i]), _Num(f"{blobs.t2[i]:.2f}", blobs.t2[i]),
                     QtWidgets.QTableWidgetItem(c.cls or "unclassified"),
                     _Num(c.tier + (" *" if c.co_elution else ""), TIER_ORDER.get(c.tier, 0)),
                     _Num(f"{c.score:.2f}", c.score),
                     QtWidgets.QTableWidgetItem(f"{c.runner_up} ({c.runner_up_score:.2f})" if c.runner_up else ""),
                     QtWidgets.QTableWidgetItem(c.limiting or ""))
            cells[3].setForeground(QtGui.QBrush(QtGui.QColor(*group_colour(c.group))))
            if c.co_elution:
                cells[4].setToolTip("* co-eluting: the spectrum is a mixture, so the tier is lowered")
            for col, item in enumerate(cells):
                self.table.setItem(r, col, item)
        self.table.setSortingEnabled(True)
        self._filling = False

    def select_blob(self, index: int) -> None:
        """Select blob `index` in the table (and show why it got its class)."""
        for r in range(self.table.rowCount()):
            item = self.table.item(r, 0)
            if item is not None and int(item.data(QtCore.Qt.ItemDataRole.UserRole)) == index:
                self.table.setCurrentCell(r, 0)
                return

    # --------------------------------------------------------------- events
    def _blob_at(self, row: int) -> int | None:
        item = self.table.item(row, 0)
        return None if item is None else int(item.data(QtCore.Qt.ItemDataRole.UserRole))

    def _row_changed(self, row: int, *_a) -> None:
        if self._filling or row < 0 or self._result is None:
            return
        i = self._blob_at(row)
        if i is not None:
            self.explain.setPlainText(self._result.explain(i))
            self.blobSelected.emit(i)

    def _row_activated(self, row: int, _col: int) -> None:
        i = self._blob_at(row)
        if i is not None:
            self.blobActivated.emit(i)

    def _class_toggled(self, item) -> None:
        if not self._filling:
            self.classToggled.emit(item.data(QtCore.Qt.ItemDataRole.UserRole),
                                   item.checkState() == QtCore.Qt.CheckState.Checked)
