"""📍 The anchoring window: check the phase, click matching peaks run by run, press Anchor, check the result.

Three pages, in order.

1. **Phase check.** An anchor records where a peak sits inside its modulation. Change a run's phase
   afterwards and every peak moves in the second dimension, so its anchors are wrong. The window therefore
   asks the user to look at each run's phase first and tick it; a tick is cleared whenever that run's phase
   or the period changes (`AlignState.phase_checked`).
2. **Pick anchors.** A large header names the run being clicked ("Run 3 of 8"), whether it is the reference
   or a sample, a strip shows which runs have this anchor, and a box says how many anchors the alignment
   needs and how many are complete. The image is the run's own, so a click here is a click on the run named
   at the top. The picking sequence itself (reference first, then every other run, skip, Esc) is the
   main window's `start_add_anchor` machinery, unchanged; this page is its face.
3. **Check the result.** Pressing **Anchor** fits every run's transform (`AlignState.apply_hold`). Until
   then the clicks are only recorded, so the images and any blank subtraction do not move under the user's
   mouse. The result page shows, for each run, the reference and the run overlaid before and after
   alignment, with numbers: how far the anchors were apart before, the leave-one-out error after, and how
   much more alike the images became.

The window changes nothing but through `AlignState` and the main window's existing, undoable alignment
edits.
"""

from __future__ import annotations

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from ..core.align import anchor_offsets, anchor_requirements, cross_validate, image_similarity
from ..core.compare import resample, source_pixels
from .theme import fit_to_screen
from .views.compare import CompareView, _ImagePane

PAGE_PHASE, PAGE_PICK, PAGE_RESULT = 0, 1, 2
CHIP_COLOURS = {"done": "#009E73", "current": "#0072B2", "skipped": "#E69F00", "pending": "#999999"}
SEARCH_COLUMNS, SEARCH_T2 = 1, 0.150          # the refinement defaults, drawn as the box a click can snap within


def _heading(text: str, size_pt: float = 14.0) -> QtWidgets.QLabel:
    label = QtWidgets.QLabel(text)
    label.setWordWrap(True)
    label.setStyleSheet(f"font-size: {size_pt:g}pt; font-weight: bold;")
    return label


class AnchorWindow(QtWidgets.QDialog):
    def __init__(self, main_window, parent=None):
        super().__init__(parent)
        self.setWindowTitle("📍 Anchor the runs")
        self.setWindowFlag(QtCore.Qt.WindowType.WindowMaximizeButtonHint, True)
        fit_to_screen(self, 1100, 760)
        self.mw = main_window
        self.project = main_window.project
        self.align = main_window.align
        self._shown_run: str | None = None         # the run whose image the pick pane holds
        self._shown_key: tuple | None = None
        self._pane_geometry: tuple | None = None   # (t1 of column 0, period, pitch) of that image
        self._report: dict = {}
        self._thumbs: dict = {}
        self.pages = QtWidgets.QStackedWidget()
        self.pages.addWidget(self._build_phase_page())
        self.pages.addWidget(self._build_pick_page())
        self.pages.addWidget(self._build_result_page())
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(self.pages)
        self.align.changed.connect(self.refresh)
        self.align.phaseStatusChanged.connect(self.refresh)
        self.project.runsChanged.connect(self.refresh)

    # ================================================================ pages
    def _build_phase_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        self.phase_heading = _heading("Check the phase of every run before anchoring.")
        note = QtWidgets.QLabel(
            "An anchor records where a peak sits inside its modulation. If a run's phase is changed afterwards, "
            "every peak in it moves in the second dimension and its anchors are wrong. Look at each run below: "
            "a correct phase keeps every peak whole, a wrong one cuts peaks across the top and bottom edges. "
            "Tick each run you are happy with. A tick is removed if that run's phase or the period changes.")
        note.setWordWrap(True)
        self.phase_table = QtWidgets.QTableWidget(0, 5)
        self.phase_table.setHorizontalHeaderLabels(["Checked", "Run", "Phase (s)", "Found by", "Folded TIC"])
        self.phase_table.verticalHeader().setVisible(False)
        self.phase_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.phase_table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        self.phase_table.horizontalHeader().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.phase_table.setIconSize(QtCore.QSize(200, 56))
        self.phase_table.verticalHeader().setDefaultSectionSize(62)
        self.phase_table.itemChanged.connect(self._phase_item_changed)
        self.phase_warning = QtWidgets.QLabel("")
        self.phase_warning.setWordWrap(True)
        self.phase_warning.setStyleSheet("color: #D55E00; font-weight: bold;")
        tick_all = QtWidgets.QPushButton("Tick all runs")
        tick_all.clicked.connect(self._tick_all)
        open_fold = QtWidgets.QPushButton("Open Fold and phase…")
        open_fold.clicked.connect(lambda: self.mw.workflow.show_step("fold"))
        self.continue_btn = QtWidgets.QPushButton("Continue to anchoring")
        self.continue_btn.setObjectName("go")
        self.continue_btn.clicked.connect(self._continue)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(tick_all)
        row.addWidget(open_fold)
        row.addStretch(1)
        row.addWidget(self.continue_btn)
        lay = QtWidgets.QVBoxLayout(page)
        lay.addWidget(self.phase_heading)
        lay.addWidget(note)
        lay.addWidget(self.phase_table, 1)
        lay.addWidget(self.phase_warning)
        lay.addLayout(row)
        return page

    def _build_pick_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        self.run_heading = _heading("", 16.0)
        self.role_badge = QtWidgets.QLabel("")
        self.role_badge.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.role_badge.setMinimumWidth(110)
        self.instruction = QtWidgets.QLabel("")
        self.instruction.setWordWrap(True)
        self.instruction.setStyleSheet("font-size: 11pt;")
        head = QtWidgets.QHBoxLayout()
        head.addWidget(self.run_heading, 1)
        head.addWidget(self.role_badge)
        self.chip_row = QtWidgets.QHBoxLayout()
        self.chip_row.setSpacing(4)
        self.chips: list[QtWidgets.QLabel] = []
        self.requirements = QtWidgets.QLabel("")
        self.requirements.setWordWrap(True)
        self.requirements.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        self.requirements.setMargin(6)
        self.coverage = QtWidgets.QLabel("")
        self.coverage.setWordWrap(True)

        self.pane = _ImagePane()
        self.pane.plot.scene().sigMouseClicked.connect(self._pane_clicked)
        self.ghost = self._marker_item((213, 94, 0), "x", 16)
        self.box = self.pane.plot.plot([], [], pen=self._pen((213, 94, 0), 1.2, dash=True))
        self.placed = self._marker_item((0, 158, 115), "o", 12)
        for item in (self.ghost, self.placed):
            self.pane.plot.addItem(item)
        self.pane_note = QtWidgets.QLabel("")
        self.pane_note.setWordWrap(True)

        self.ion_edit = QtWidgets.QLineEdit()
        self.ion_edit.setPlaceholderText("TIC or one ion, e.g. 178")
        self.ion_edit.setMaximumWidth(150)
        self.ion_edit.returnPressed.connect(self._set_ion)
        ion_btn = QtWidgets.QPushButton("Set signal")
        ion_btn.clicked.connect(self._set_ion)
        self.start_btn = QtWidgets.QPushButton("Start picking (A)")
        self.start_btn.clicked.connect(self.mw.start_add_anchor)
        self.skip_btn = QtWidgets.QPushButton("Skip this peak in this run")
        self.skip_btn.clicked.connect(self.mw.skip_anchor_run)
        self.finish_btn = QtWidgets.QPushButton("Finish picking (Esc)")
        self.finish_btn.clicked.connect(self.mw.cancel_pick)
        back = QtWidgets.QPushButton("‹ Phase check")
        back.clicked.connect(lambda: self.pages.setCurrentIndex(PAGE_PHASE))
        self.anchor_btn = QtWidgets.QPushButton("ANCHOR")
        self.anchor_btn.setObjectName("go")
        self.anchor_btn.setMinimumHeight(40)
        self.anchor_btn.setStyleSheet("font-size: 13pt;")
        self.anchor_btn.clicked.connect(self._anchor)
        self.pending_note = QtWidgets.QLabel("")
        self.pending_note.setWordWrap(True)

        controls = QtWidgets.QHBoxLayout()
        controls.addWidget(QtWidgets.QLabel("Signal"))
        controls.addWidget(self.ion_edit)
        controls.addWidget(ion_btn)
        controls.addSpacing(12)
        controls.addWidget(self.start_btn)
        controls.addWidget(self.skip_btn)
        controls.addWidget(self.finish_btn)
        controls.addStretch(1)
        lay = QtWidgets.QVBoxLayout(page)
        lay.addLayout(head)
        lay.addWidget(self.instruction)
        lay.addLayout(self.chip_row)
        lay.addWidget(self.requirements)
        lay.addWidget(self.coverage)
        lay.addWidget(self.pane, 1)
        lay.addWidget(self.pane_note)
        lay.addLayout(controls)
        bottom = QtWidgets.QHBoxLayout()
        bottom.addWidget(back)
        bottom.addWidget(self.pending_note, 1)
        bottom.addWidget(self.anchor_btn)
        lay.addLayout(bottom)
        return page

    def _build_result_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        self.result_heading = _heading("Check what the alignment did.")
        self.result_note = QtWidgets.QLabel(
            "Left: the reference in magenta and the run in green as they were, shared peaks look white only where "
            "they already agreed. Right: after alignment. Peaks that line up turn white.")
        self.result_note.setWordWrap(True)
        self.result_table = QtWidgets.QTableWidget(0, 8)
        self.result_table.setHorizontalHeaderLabels(
            ["Run", "Transform", "Anchors", "Apart before (t₁ s / t₂ ms)", "Leave-one-out error (t₁ s / t₂ ms)",
             "Image match before", "Image match after", "Status"])
        self.result_table.verticalHeader().setVisible(False)
        self.result_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.result_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.result_table.setMaximumHeight(190)
        self.result_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
        self.result_table.itemSelectionChanged.connect(self._result_run_changed)
        self.before_view = CompareView()
        self.after_view = CompareView()
        self.before_view.image_source = self._before_image
        self.after_view.image_source = self._after_image
        for view in (self.before_view, self.after_view):
            view.mode.setCurrentIndex(view.mode.findData("overlay"))
        before_box = QtWidgets.QGroupBox("Before")
        QtWidgets.QVBoxLayout(before_box).addWidget(self.before_view)
        after_box = QtWidgets.QGroupBox("After")
        QtWidgets.QVBoxLayout(after_box).addWidget(self.after_view)
        views = QtWidgets.QHBoxLayout()
        views.addWidget(before_box)
        views.addWidget(after_box)
        back = QtWidgets.QPushButton("‹ Back to anchors")
        back.clicked.connect(self._back_to_anchors)
        done = QtWidgets.QPushButton("Done")
        done.setObjectName("go")
        done.clicked.connect(self.close)
        undo = QtWidgets.QPushButton("Undo last alignment edit")
        undo.clicked.connect(lambda: self.mw.undo.undo())
        bottom = QtWidgets.QHBoxLayout()
        bottom.addWidget(back)
        bottom.addWidget(undo)
        bottom.addStretch(1)
        bottom.addWidget(done)
        lay = QtWidgets.QVBoxLayout(page)
        lay.addWidget(self.result_heading)
        lay.addWidget(self.result_note)
        lay.addWidget(self.result_table)
        lay.addLayout(views, 1)
        lay.addLayout(bottom)
        return page

    @staticmethod
    def _pen(rgb, width, dash=False):
        import pyqtgraph as pg
        return pg.mkPen(rgb, width=width, style=QtCore.Qt.PenStyle.DashLine if dash else QtCore.Qt.PenStyle.SolidLine)

    @staticmethod
    def _marker_item(rgb, symbol, size):
        import pyqtgraph as pg
        return pg.ScatterPlotItem(pen=pg.mkPen(rgb, width=2), brush=None, symbol=symbol, size=size)

    def set_background(self, colour) -> None:
        self.pane.set_background(colour)
        self.before_view.set_background(colour)
        self.after_view.set_background(colour)

    # ================================================================ opening
    def open_window(self) -> None:
        """Show the window on the page that fits: phase check first, then picking."""
        if self.pages.currentIndex() == PAGE_RESULT and not self.mw._anchor_stage:
            pass
        elif self._all_checked():
            self.pages.setCurrentIndex(PAGE_PICK)
        else:
            self.pages.setCurrentIndex(PAGE_PHASE)
        self.refresh()
        self.show()
        self.raise_()
        self.activateWindow()

    def go_to_pick_page(self) -> None:
        self.pages.setCurrentIndex(PAGE_PICK)
        self.refresh()
        self.show()
        self.raise_()

    def runs_to_align(self) -> list[str]:
        return [self.project.entries[i].name for i in self.project.loaded()]

    def _all_checked(self) -> bool:
        names = self.runs_to_align()
        return bool(names) and all(n in self.align.phase_checked for n in names)

    # ================================================================ phase page
    def _fill_phase_table(self) -> None:
        names = self.runs_to_align()
        stale = set(self.align.stale_phase_runs())
        self.phase_table.blockSignals(True)
        self.phase_table.setRowCount(len(names))
        for row, name in enumerate(names):
            index = self.project.index_of(name)
            entry = self.project.entries[index]
            check = QtWidgets.QTableWidgetItem()
            check.setFlags(QtCore.Qt.ItemFlag.ItemIsUserCheckable | QtCore.Qt.ItemFlag.ItemIsEnabled)
            check.setCheckState(QtCore.Qt.CheckState.Checked if name in self.align.phase_checked
                                else QtCore.Qt.CheckState.Unchecked)
            check.setData(QtCore.Qt.ItemDataRole.UserRole, name)
            self.phase_table.setItem(row, 0, check)
            label = name + ("   (anchored with a different phase)" if name in stale else "")
            item = QtWidgets.QTableWidgetItem(label)
            if name in stale:
                item.setForeground(QtGui.QBrush(QtGui.QColor("#D55E00")))
            self.phase_table.setItem(row, 1, item)
            self.phase_table.setItem(row, 2, QtWidgets.QTableWidgetItem(f"{entry.phase:.3f}"))
            self.phase_table.setItem(row, 3, QtWidgets.QTableWidgetItem(str(entry.phase_method)))
            thumb = QtWidgets.QTableWidgetItem()
            icon = self._thumbnail(index)
            if icon is not None:
                thumb.setIcon(icon)
            self.phase_table.setItem(row, 4, thumb)
        self.phase_table.blockSignals(False)
        ready = self._all_checked()
        self.continue_btn.setEnabled(ready)
        self.phase_warning.setText(
            ("Anchored with a different phase: " + ", ".join(sorted(stale)) + ". Pick their anchors again or press "
             "Anchor to refit.") if stale else ("" if ready else "Tick every run to continue."))

    def _thumbnail(self, index: int):
        """A small picture of the folded TIC, so a phase that cuts peaks in two can be seen."""
        entry = self.project.entries[index]
        key = (id(entry.run), self.project.period, entry.phase)
        if key in self._thumbs:
            return self._thumbs[key]
        icon = self._draw_thumbnail(index)
        self._thumbs[key] = icon
        return icon

    def _draw_thumbnail(self, index: int):
        try:
            fold = self.project.fold(index)
            run = self.project.display_run(index)
            image = np.sqrt(np.clip(np.nan_to_num(fold.image(run.tic)), 0, None))
        except Exception:                           # noqa: BLE001  (a missing picture must not stop the page)
            return None
        if image.size == 0 or image.max() <= 0:
            return None
        step_r, step_c = max(image.shape[0] // 56, 1), max(image.shape[1] // 200, 1)
        small = image[::step_r, ::step_c]
        scaled = (255 * np.clip(small / np.percentile(small[small > 0], 99.5) if (small > 0).any() else small, 0, 1)).astype(np.uint8)
        scaled = np.ascontiguousarray(scaled[::-1])
        qimage = QtGui.QImage(scaled.data, scaled.shape[1], scaled.shape[0], scaled.shape[1],
                              QtGui.QImage.Format.Format_Grayscale8).copy()
        return QtGui.QIcon(QtGui.QPixmap.fromImage(qimage))

    def _phase_item_changed(self, item) -> None:
        if item.column() != 0:
            return
        name = item.data(QtCore.Qt.ItemDataRole.UserRole)
        self.align.set_phase_checked(name, item.checkState() == QtCore.Qt.CheckState.Checked)

    def _tick_all(self) -> None:
        for name in self.runs_to_align():
            self.align.phase_checked.add(name)
        self.align.phaseStatusChanged.emit()

    def _continue(self) -> None:
        if not self._all_checked():
            return
        self.pages.setCurrentIndex(PAGE_PICK)
        self.mw.start_add_anchor()

    # ================================================================ pick page
    def _samples(self) -> list[str]:
        return [n for n in self.runs_to_align() if n != self.align.reference]

    def _complete_anchors(self) -> list[str]:
        samples = self._samples()
        return [name for name in self.align.anchor_names()
                if samples and all(any(p.name == name for p in self.align.points.get(run, [])) for run in samples)]

    def _refresh_pick_page(self) -> None:
        mw, picking = self.mw, self.mw._anchor_stage is not None
        pending = mw._anchor_pending if picking else {}
        order = pending.get("order", [self.align.reference] + self._samples() if self.align.reference else [])
        cursor = pending.get("cursor", 0)
        name = pending.get("name") or self.align.next_anchor_name()
        run = order[cursor] if picking and order else None
        reference = self.align.reference

        if picking and run:
            is_ref = run == reference
            self.run_heading.setText(f"Run {cursor + 1} of {len(order)}: {run}")
            self.role_badge.setText("REFERENCE" if is_ref else "SAMPLE")
            self.role_badge.setStyleSheet(
                f"background: {'#0072B2' if is_ref else '#009E73'}; color: white; font-weight: bold; padding: 6px;"
                " border-radius: 4px;")
            self.instruction.setText(
                (f"Anchor {name}. Click a clear peak on the reference. Every other run will be matched to it."
                 if is_ref else
                 f"Anchor {name}. Click the same peak as on the reference (the orange cross, inside the dashed box "
                 "where a click snaps to the brightest nearby pixel), or skip if this run does not have it."))
        else:
            self.run_heading.setText("Not picking")
            self.role_badge.setText("")
            self.role_badge.setStyleSheet("")
            self.instruction.setText("Press Start picking to click matching peaks run by run, or press ANCHOR "
                                     "to align from the clicks already recorded.")
        self.start_btn.setVisible(not picking)
        self.skip_btn.setVisible(picking and run is not None and run != reference)
        self.finish_btn.setVisible(picking)

        # Chips: one per run, for the anchor being picked.
        while len(self.chips) < len(order):
            chip = QtWidgets.QLabel("")
            chip.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
            self.chip_row.addWidget(chip)
            self.chips.append(chip)
        for k, chip in enumerate(self.chips):
            chip.setVisible(k < len(order))
            if k >= len(order):
                continue
            state = ("current" if picking and k == cursor else
                     "done" if any(p.name == name for p in self.align.points.get(order[k], [])) else
                     "skipped" if picking and k < cursor else "pending")
            chip.setText(f"{k + 1}  {order[k]}")
            chip.setToolTip({"done": "has this anchor", "current": "being picked now", "skipped": "skipped",
                             "pending": "still to come"}[state])
            chip.setStyleSheet(f"background: {CHIP_COLOURS[state]}; color: white; padding: 3px 6px;"
                               " border-radius: 4px; font-size: 8pt;")

        # How many anchors are needed.
        complete = self._complete_anchors()
        need = anchor_requirements(len(complete))
        kinds = "; ".join(f"{f['label']}: {f['least']} or more" for f in need["fits"])
        auto = {"shift": "a shift", "affine": "a stretch and shift", "smooth": "a smooth curve", None: "nothing yet"}[need["auto"]]
        self.requirements.setText(
            f"<b>Anchors complete on every sample: {len(complete)}.</b> Needed: {kinds}. "
            f"With {len(complete)}, automatic choice may use <b>{auto}</b>. "
            f"Recommended: {need['recommended'][0]} to {need['recommended'][1]}, spread from early to late in the run.")
        self.coverage.setText(self._coverage_text())

        if run is not None and picking:
            self._show_run_image(run)
        else:
            self._update_markers(None)
        n_pending = len(self.align.pending_runs())
        self.anchor_btn.setEnabled(bool(complete) or bool(self.align.points) and any(
            r != reference for r in self.align.points))
        self.pending_note.setText(
            f"Clicks recorded for {n_pending} run(s) have not been applied. Press ANCHOR to align them." if n_pending
            else ("All recorded clicks are applied." if self.align.transforms else ""))

    def _coverage_text(self) -> str:
        reference = self.align.reference
        points = self.align.points.get(reference, []) if reference else []
        index = self.project.index_of(reference) if reference else None
        if not points or index is None:
            return ""
        fold = self.project.fold(index)
        lo, hi = float(fold.t1[0]), float(fold.t1[-1])
        times = [p.ref_t1 for p in points]
        span = (max(times) - min(times)) / max(hi - lo, 1e-9)
        text = f"Anchors cover {min(times) / 60:.1f} to {max(times) / 60:.1f} min of {lo / 60:.1f} to {hi / 60:.1f} min ({span:.0%})."
        if len(points) >= 2 and span < 1 / 3:
            text += " <b>They sit in one third of the run or less: add anchors early and late.</b>"
        return text

    # ---- image of the run being clicked
    def _show_run_image(self, run: str) -> None:
        expr = self.mw.renderer.expr
        key = (run, expr.text(), id(self.project.display_run(self.project.index_of(run))))
        if key == self._shown_key:
            self._update_markers(run)
            return
        index = self.project.index_of(run)
        if index is None or self.project.entries[index].run is None:
            return
        self._shown_key = None
        self._shown_run = None
        self.pane_note.setText(f"Drawing {run}…")
        data, fold = self.project.display_run(index), self.project.fold(index)

        def work(token, progress):
            return expr.image(data, fold)

        def done(image):
            if self.mw._anchor_stage is None and self.pages.currentIndex() != PAGE_PICK:
                return
            self.pane.show_mono(image, fold.t1, fold.row_pitch, fold.spec.period, "magma")
            self._shown_run, self._shown_key = run, key
            self._pane_geometry = (float(fold.t1[0]), float(fold.spec.period), float(fold.row_pitch), fold.n_cols)
            self.pane_note.setText(f"{run}: signal {expr.label()}. Click a peak to place the anchor.")
            self._update_markers(run)

        def failed(message):
            self.pane_note.setText(f"Could not draw {run}: {message}")

        self.project.jobs.submit(f"📍 Drawing {run}", work, done, failed, latest_key="anchor-window-image", quiet=True)

    def _update_markers(self, run: str | None) -> None:
        """Orange cross and dashed box where the reference peak is; green rings for this run's anchors."""
        ghost, box, placed = ([], []), ([], []), ([], [])
        reference = self.align.reference
        if run is not None and reference:
            name = self.mw._anchor_pending.get("name") if self.mw._anchor_stage else None
            ref_point = next((p for p in self.align.points.get(reference, []) if p.name == name), None)
            if ref_point is not None and run != reference and self._pane_geometry is not None:
                period = self._pane_geometry[1]
                x0, x1 = (ref_point.ref_t1 - SEARCH_COLUMNS * period) / 60.0, (ref_point.ref_t1 + SEARCH_COLUMNS * period) / 60.0
                y0, y1 = ref_point.ref_t2 - SEARCH_T2, ref_point.ref_t2 + SEARCH_T2
                ghost = ([ref_point.ref_t1 / 60.0], [ref_point.ref_t2])
                box = ([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0])
            mine = self.align.points.get(run, [])
            placed = ([p.run_t1 / 60.0 for p in mine], [p.run_t2 for p in mine])
        self.ghost.setData(x=ghost[0], y=ghost[1])
        self.box.setData(*box)
        self.placed.setData(x=placed[0], y=placed[1])

    def _pane_clicked(self, event) -> None:
        if event.button() != QtCore.Qt.MouseButton.LeftButton or self._shown_run is None or self._pane_geometry is None:
            return
        if not self.pane.plot.sceneBoundingRect().contains(event.scenePos()):
            return
        point = self.pane.plot.vb.mapSceneToView(event.scenePos())
        t1_first, period, pitch, n_cols = self._pane_geometry
        col = int(np.floor((point.x() * 60.0 - t1_first) / period))
        if not (0 <= col < n_cols):
            return
        self.mw.anchor_pane_click(self._shown_run, col, float(np.rint(point.y() / pitch) * pitch))

    def _set_ion(self) -> None:
        text = self.ion_edit.text().strip()
        if text:
            self.mw.set_ion_expression(text)
            self._shown_key = None
            self.refresh()

    # ================================================================ anchoring
    def _anchor(self) -> None:
        self.mw.apply_anchors()
        self.show_results()

    def _back_to_anchors(self) -> None:
        self.pages.setCurrentIndex(PAGE_PICK)
        self.refresh()

    def show_results(self) -> None:
        self.pages.setCurrentIndex(PAGE_RESULT)
        self.refresh()
        self._compute_report()

    def _result_runs(self) -> list[str]:
        return [n for n in self._samples() if n in self.align.points]

    def _fill_result_table(self) -> None:
        runs = self._result_runs()
        self.result_table.blockSignals(True)
        self.result_table.setRowCount(len(runs))
        for row, run in enumerate(runs):
            pts = self.align.points.get(run, [])
            tf = self.align.transforms.get(run)
            dt1, dt2 = anchor_offsets(pts)
            if tf is not None:
                e1, e2 = cross_validate(pts, tf.kind)
                loo = f"{e1:.2f} / {e2 * 1e3:.0f}" if np.isfinite(e1) else "too few to check"
            else:
                loo = ""
            rep = self._report.get(run, {})
            before, after = rep.get("before"), rep.get("after")
            status = self._status(run, tf, before, after)
            values = (run, tf.kind if tf else "none", str(len(pts)), f"{dt1:.2f} / {dt2 * 1e3:.0f}", loo,
                      "" if before is None or not np.isfinite(before) else f"{before:.2f}",
                      "" if after is None or not np.isfinite(after) else f"{after:.2f}", status)
            for c, v in enumerate(values):
                item = QtWidgets.QTableWidgetItem(v)
                if c == 7:
                    ok = status == "OK"
                    item.setForeground(QtGui.QBrush(QtGui.QColor("#009E73" if ok else "#D55E00")))
                self.result_table.setItem(row, c, item)
        self.result_table.blockSignals(False)
        if runs and not self.result_table.selectedItems():
            self.result_table.selectRow(0)

    def _status(self, run, tf, before, after) -> str:
        if run in self.align.errors:
            return self.align.errors[run]
        if tf is None:
            return "Not aligned yet"
        if not tf.monotonic:
            return "Not monotonic: " + tf.warning
        if before is not None and after is not None and np.isfinite(before) and np.isfinite(after) and after + 0.02 < before:
            return "Worse after alignment: add anchors across the run"
        if run in self.align.stale_phase_runs():
            return "Phase changed since anchoring"
        if len(self.align.points.get(run, [])) < 2:
            return "OK (one anchor: a shift only)"
        return "OK"

    def _compute_report(self) -> None:
        """Image similarity before and after, for each run, in a job."""
        reference = self.align.reference
        runs = self._result_runs()
        ref_index = self.project.index_of(reference) if reference else None
        if ref_index is None or not runs:
            self._fill_result_table()
            return
        expr = self.mw.renderer.expr
        project, align = self.project, self.align

        def work(token, progress):
            ref_fold = project.fold(ref_index)
            ref_img = expr.image(project.display_run(ref_index), ref_fold)
            out = {}
            for k, run in enumerate(runs):
                token.check()
                progress(k, len(runs), f"Comparing {run}")
                index = project.index_of(run)
                fold = project.fold(index)
                raw = expr.image(project.display_run(index), fold)
                identity = resample(raw, source_pixels(fold, ref_fold, lambda t1, t2: (t1, t2)))
                aligned = align.image_for(index, expr)
                out[run] = dict(before=image_similarity(ref_img, identity),
                                after=image_similarity(ref_img, aligned) if aligned is not None else float("nan"))
            return out

        def done(report):
            self._report = report
            self._fill_result_table()

        self._fill_result_table()
        self.project.jobs.submit("📍 Comparing before and after", work, done, lambda m: self.mw.say(f"Comparison failed: {m}"),
                                 latest_key="anchor-window-report", quiet=True)
        self._result_run_changed()

    def _result_run_changed(self) -> None:
        rows = sorted({i.row() for i in self.result_table.selectedIndexes()})
        runs = self._result_runs()
        reference = self.align.reference
        if not rows or rows[0] >= len(runs) or not reference:
            return
        ref_index, index = self.project.index_of(reference), self.project.index_of(runs[rows[0]])
        for view in (self.before_view, self.after_view):
            view.set_runs(ref_index, index, [])

    def _before_image(self, index: int):
        """Run `index` on the reference's grid with no alignment applied."""
        ref_index = self.project.index_of(self.align.reference)
        p, expr = self.project, self.mw.renderer.expr
        if ref_index is None or not (0 <= index < len(p.entries)) or p.entries[index].run is None:
            return None
        ref_fold = p.fold(ref_index)
        img = expr.image(p.display_run(index), p.fold(index))
        if index != ref_index:
            img = resample(img, source_pixels(p.fold(index), ref_fold, lambda t1, t2: (t1, t2)))
        return img, ref_fold.t1, ref_fold.row_pitch, ref_fold.spec.period, p.entries[index].name

    def _after_image(self, index: int):
        ref_index = self.project.index_of(self.align.reference)
        if ref_index is None or not (0 <= index < len(self.project.entries)):
            return None
        img = self.align.image_for(index, self.mw.renderer.expr)
        if img is None:
            return None
        ref_fold = self.project.fold(ref_index)
        return img, ref_fold.t1, ref_fold.row_pitch, ref_fold.spec.period, self.project.entries[index].name

    # ================================================================ refresh
    def refresh(self) -> None:
        if not self.project.entries:
            return
        if self.pages.currentIndex() == PAGE_PHASE:
            self._fill_phase_table()
        self._refresh_pick_page()
        if self.pages.currentIndex() == PAGE_RESULT:
            self._fill_result_table()

    def hideEvent(self, event) -> None:
        for view in (self.before_view, self.after_view):
            view.blink_timer.stop()
        super().hideEvent(event)

    def closeEvent(self, event) -> None:
        # Closing the window ends picking, as pressing Esc does; the clicks are kept.
        if self.mw._anchor_stage is not None:
            self.mw.cancel_pick()
        super().closeEvent(event)

