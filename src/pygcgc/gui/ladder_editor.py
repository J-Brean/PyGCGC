"""📏 Draw the alkane ladder on the 2D chromatogram and pull it onto the alkanes.

The ladder is the line of n-alkanes (C10, C11, C12, ...) that every retention index is measured against. Here it is
drawn over the 2D image as a curve with a dot at each of a few *control carbons*, and you move it by hand:

- the dots at the two ends stretch, shorten or turn the whole ladder (the dots between follow in proportion);
  hold Ctrl to move only that end;
- a dot in the middle bends the ladder there, which is how an S-shaped ladder is followed;
- the diamond moves everything at once;
- double-click the line, near a tick, to add a dot; right-click a dot to remove it or change its carbon number;
- "A dot on every carbon" gives each carbon its own dot, for fine work;
- "Snap to peaks" then moves each carbon onto the nearest real m/z 57 apex (the search the program uses when it
  finds the ladder itself), so a rough drawing is enough.

Each time you let go, the ladder is stored for the run (one undo step), and everything that reads it follows: the
retention indices of the library matcher (whose results are then marked out of date), the hydrocarbon
assignments, the stencils and the markers. This works with either peak finder, because it needs only the run.
The geometry is in `core/hydrocarbons.py` (`ladder_from_control_points`, `snap_ladder`); this file is the drawing.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

from ..core.hydrocarbons import MIN_STEP_S, control_points_of, ladder_from_control_points, sample_ladder
from .commands import EditLadder
from .theme import fit_to_screen

FEW_HANDLES = 5                      # control carbons when the ladder is drawn "simply"
HANDLE_RGB, HOVER_RGB = (255, 120, 0), (255, 255, 255)
MOVER_RGB = (0, 200, 255)
DEFAULT_RANGE = (10, 30)             # carbons offered when there is no ladder yet
HINT = ("<b>Drag the dots onto the alkanes.</b> The dots at the ends stretch or turn the whole ladder (hold Ctrl to move "
        "just that end); the dots between bend it; the blue diamond moves it all. Double-click the line near a tick to "
        "add a dot; right-click a dot to remove it or change its carbon number.")


class _Handle(pg.TargetItem):
    """A draggable dot that also reports a right-click."""
    contextRequested = QtCore.Signal(object, object)           # this handle, the screen position

    def mouseClickEvent(self, ev):
        if ev.button() == QtCore.Qt.MouseButton.RightButton:
            ev.accept()
            self.contextRequested.emit(self, ev.screenPos())
        else:
            ev.ignore()


class LadderOverlay(QtCore.QObject):
    """The ladder as items on the 2D view. It holds control points `[carbon, t1 (s), t2 (s)]` and nothing else."""
    edited = QtCore.Signal(str)                                 # a change is finished; the text names it for undo
    inconsistent = QtCore.Signal(str)                           # a change was refused, and why

    def __init__(self, view, parent=None):
        super().__init__(parent)
        self.view, self.plot = view, view.plot
        self.points: list[list[float]] = []
        self.active = False
        self._handles: list[_Handle] = []
        self._mover: pg.TargetItem | None = None
        self._kind = ""                                         # what the drag in progress is, for the undo list
        self._start = None                                    # the control points when a drag began
        self._busy = False                                      # True while the dots are being placed by code
        self._labels: list[pg.TextItem] = []
        self._mover_lift = 0.0
        self.halo = pg.PlotDataItem(pen=pg.mkPen(0, 0, 0, 170, width=5))
        self.curve = pg.PlotDataItem(pen=pg.mkPen(255, 255, 255, 240, width=2.2))
        self.ticks = pg.ScatterPlotItem(size=7, symbol="+", pen=pg.mkPen(255, 255, 255, 255, width=1.5), brush=None)
        self.found = pg.ScatterPlotItem(size=9, symbol="x", pen=pg.mkPen(0, 0, 0, 220, width=1.6), brush=None)
        for item in (self.halo, self.curve, self.ticks, self.found):
            item.setZValue(30)
            item.setVisible(False)
            self.plot.addItem(item, ignoreBounds=True)
        self.plot.scene().sigMouseClicked.connect(self._scene_clicked)

    # ------------------------------------------------------------ what is shown
    def load(self, points, found=()) -> None:
        """Show control points `(carbon, t1, t2)`; `found` are `(t1, t2)` of measured apexes, drawn as small crosses."""
        self.points = [[float(c), float(t1), float(t2)] for c, t1, t2 in sorted(points)]
        self.found.setData(x=[t1 / 60.0 for t1, _ in found], y=[t2 for _, t2 in found])
        self._rebuild_handles()
        self._redraw()

    def show(self, on: bool) -> None:
        self.active = bool(on)
        for item in (self.halo, self.curve, self.ticks, self.found, *self._handles, self._mover, *self._labels):
            if item is not None:
                item.setVisible(self.active)

    def series(self):
        """The ladder as an `AlkaneSeries` (raises `ValueError` if the control points are not usable)."""
        return ladder_from_control_points(self.points, "smooth")

    # ------------------------------------------------------------ building
    def _rebuild_handles(self) -> None:
        for h in self._handles:
            self.plot.removeItem(h)
        if self._mover is not None:
            self.plot.removeItem(self._mover)
        self._handles = []
        yr = self.plot.vb.viewRange()[1]
        self._mover_lift = 0.06 * (yr[1] - yr[0])
        for k, (_carbon, t1, t2) in enumerate(self.points):
            h = _Handle(pos=(t1 / 60.0, t2), size=15, symbol="o", movable=True,
                        pen=pg.mkPen(*HANDLE_RGB, width=2), hoverPen=pg.mkPen(*HOVER_RGB, width=2.5),
                        brush=pg.mkBrush(*HANDLE_RGB, 90), hoverBrush=pg.mkBrush(*HANDLE_RGB, 200))
            h.setZValue(40)
            h.sigPositionChanged.connect(lambda _item, k=k: self._dragged(k))
            h.sigPositionChangeFinished.connect(lambda _item: self._finished())
            h.contextRequested.connect(self._context)
            self.plot.addItem(h, ignoreBounds=True)
            h.setVisible(self.active)
            self._handles.append(h)
        self._mover = pg.TargetItem(pos=(0, 0), size=19, symbol="d", movable=True, label="move all",
                                    pen=pg.mkPen(*MOVER_RGB, width=2), hoverPen=pg.mkPen(*HOVER_RGB, width=2.5),
                                    brush=pg.mkBrush(*MOVER_RGB, 90), hoverBrush=pg.mkBrush(*MOVER_RGB, 200))
        self._mover.setZValue(40)
        self._mover.sigPositionChanged.connect(lambda _item: self._dragged(-1))
        self._mover.sigPositionChangeFinished.connect(lambda _item: self._finished())
        self.plot.addItem(self._mover, ignoreBounds=True)
        self._mover.setVisible(self.active)
        self._place_handles()

    def _place_handles(self) -> None:
        """Put every dot where the control points say, without those moves counting as drags."""
        self._busy = True
        try:
            for h, (_c, t1, t2) in zip(self._handles, self.points):
                h.setPos((t1 / 60.0, t2))
            if self._mover is not None and len(self.points) >= 2:
                mid = 0.5 * (self.points[0][0] + self.points[-1][0])
                t1, t2 = sample_ladder(self.points, [mid])
                self._mover.setPos((float(t1[0]) / 60.0, float(t2[0]) + self._mover_lift))
        finally:
            self._busy = False

    def _redraw(self) -> None:
        if len(self.points) < 2:
            for item in (self.halo, self.curve, self.ticks):
                item.setVisible(False)
            return
        try:
            first, last = self.points[0][0], self.points[-1][0]
            c = np.linspace(first, last, 200)
            t1, t2 = sample_ladder(self.points, c)
            whole = np.arange(int(first), int(last) + 1)
            tx, ty = sample_ladder(self.points, whole)
        except ValueError:
            return                                           # mid-drag and momentarily out of order: keep what is drawn
        self.halo.setData(t1 / 60.0, t2)
        self.curve.setData(t1 / 60.0, t2)
        self.ticks.setData(x=tx / 60.0, y=ty)
        for item in (self.halo, self.curve, self.ticks):
            item.setVisible(self.active)
        # Every carbon is named (every fifth if there are many); the carbons that have a dot are named in orange.
        controls = {int(p[0]) for p in self.points}
        wanted = [(int(cc), x, y, int(cc) in controls) for cc, x, y in zip(whole, tx, ty)
                  if int(cc) in controls or len(whole) <= 30 or int(cc) % 5 == 0]
        while len(self._labels) < len(wanted):
            label = pg.TextItem("", anchor=(0.5, 1.5), color=(255, 255, 255))
            label.setZValue(31)
            self.plot.addItem(label, ignoreBounds=True)
            self._labels.append(label)
        for i, label in enumerate(self._labels):
            if i < len(wanted):
                cc, x, y, control = wanted[i]
                label.fill = pg.mkBrush(*HANDLE_RGB, 215) if control else pg.mkBrush(0, 0, 0, 130)
                label.setText(f"C{cc}")
                label.setPos(x / 60.0, y)
                label.setVisible(self.active)
            else:
                label.setVisible(False)

    # ------------------------------------------------------------ dragging
    def _dragged(self, which: int) -> None:
        """A dot is being dragged: `which` is its control point, or -1 for the diamond."""
        if self._busy or len(self.points) < 2:
            return
        if self._start is None:
            self._start = [list(p) for p in self.points]
        start, n = self._start, len(self.points)
        item = self._mover if which < 0 else self._handles[which]
        x, y = item.pos().x() * 60.0, item.pos().y()
        ends_follow = (which in (0, n - 1) and n > 2
                       and not (QtWidgets.QApplication.keyboardModifiers() & QtCore.Qt.KeyboardModifier.ControlModifier))
        self._kind = "Move the ladder" if which < 0 else (
            "Stretch the ladder" if ends_follow else f"Bend the ladder at C{int(self.points[which][0])}")
        if which < 0:
            home = sample_ladder(start, [0.5 * (start[0][0] + start[-1][0])])
            d1, d2 = x - float(home[0][0]), y - self._mover_lift - float(home[1][0])
            for p, s in zip(self.points, start):
                p[1], p[2] = s[1] + d1, s[2] + d2
        elif ends_follow:
            # An end: the rest of the ladder follows in proportion to its distance from the far end.
            far = n - 1 if which == 0 else 0
            span = start[which][0] - start[far][0]
            f = [(s[0] - start[far][0]) / span for s in start]
            d1, d2 = x - start[which][1], y - start[which][2]
            lo, hi = -np.inf, np.inf                           # the ladder must keep running forwards in t1
            for j in range(n - 1):
                gap, df = start[j + 1][1] - start[j][1] - MIN_STEP_S, f[j + 1] - f[j]
                if df > 0:
                    lo = max(lo, -gap / df)
                elif df < 0:
                    hi = min(hi, -gap / df)
            d1 = float(np.clip(d1, lo, hi))
            for p, s, fj in zip(self.points, start, f):
                p[1], p[2] = s[1] + fj * d1, s[2] + fj * d2
        else:
            lo = start[which - 1][1] + MIN_STEP_S if which > 0 else -np.inf
            hi = start[which + 1][1] - MIN_STEP_S if which < n - 1 else np.inf
            self.points[which][1] = float(np.clip(x, lo, hi))
            self.points[which][2] = y
        self._place_handles()
        self._redraw()

    def _finished(self) -> None:
        if self._start is None:
            return
        start, self._start = self._start, None
        if all(abs(a[1] - b[1]) < 1e-9 and abs(a[2] - b[2]) < 1e-9 for a, b in zip(start, self.points)):
            return
        self.edited.emit(self._kind)

    # ------------------------------------------------------------ adding, removing, renumbering
    def _scene_clicked(self, ev) -> None:
        """Double-click near a tick of the line: add a dot at that carbon."""
        if not self.active or not ev.double() or ev.button() != QtCore.Qt.MouseButton.LeftButton or len(self.points) < 2:
            return
        if not self.plot.sceneBoundingRect().contains(ev.scenePos()):
            return
        first, last = int(self.points[0][0]), int(self.points[-1][0])
        whole = np.arange(first, last + 1)
        try:
            tx, ty = sample_ladder(self.points, whole)
        except ValueError:
            return
        vb = self.plot.vb
        best, best_d = None, 14.0                               # within 14 pixels of a tick
        for c, x, y in zip(whole, tx, ty):
            p = vb.mapViewToScene(QtCore.QPointF(x / 60.0, y))
            d = np.hypot(p.x() - ev.scenePos().x(), p.y() - ev.scenePos().y())
            if d < best_d and int(c) not in {int(q[0]) for q in self.points}:
                best, best_d = (int(c), float(x), float(y)), d
        if best is not None:
            self.points = sorted(self.points + [[float(best[0]), best[1], best[2]]])
            self._rebuild_handles()
            self._redraw()
            ev.accept()
            self.edited.emit(f"Add a ladder dot at C{best[0]}")

    def _context(self, handle, screen_pos) -> None:
        k = self._handles.index(handle)
        menu = QtWidgets.QMenu()
        remove = menu.addAction(f"Remove the dot at C{int(self.points[k][0])}")
        remove.setEnabled(len(self.points) > 2)
        number = menu.addAction("Change its carbon number…")
        chosen = menu.exec(screen_pos.toPoint() if hasattr(screen_pos, "toPoint") else screen_pos)
        if chosen is remove:
            self.remove(k)
        elif chosen is number:
            current = int(self.points[k][0])
            value, ok = QtWidgets.QInputDialog.getInt(None, "Carbon number", "Which n-alkane is this dot on?", current, 1, 100)
            if ok:
                self.renumber(k, value)

    def remove(self, k: int) -> bool:
        if len(self.points) <= 2:
            self.inconsistent.emit("A ladder needs at least two dots.")
            return False
        carbon = int(self.points[k][0])
        del self.points[k]
        self._rebuild_handles()
        self._redraw()
        self.edited.emit(f"Remove the ladder dot at C{carbon}")
        return True

    def renumber(self, k: int, carbon: int) -> bool:
        """Give dot `k` another carbon number, which must keep the dots in carbon order."""
        carbon = int(carbon)
        if carbon == int(self.points[k][0]):
            return True
        below = self.points[k - 1][0] if k > 0 else -np.inf
        above = self.points[k + 1][0] if k < len(self.points) - 1 else np.inf
        if not below < carbon < above:
            self.inconsistent.emit(f"C{carbon} would put the dots out of carbon order (it must lie between "
                                   f"C{int(below) if np.isfinite(below) else '…'} and C{int(above) if np.isfinite(above) else '…'}).")
            return False
        old = int(self.points[k][0])
        self.points[k][0] = float(carbon)
        self._rebuild_handles()
        self._redraw()
        self.edited.emit(f"Renumber the ladder dot C{old} as C{carbon}")
        return True

    def set_ends(self, first: int, last: int) -> bool:
        """Relabel the end dots, which sets the carbons the ladder spans."""
        if first >= last:
            self.inconsistent.emit("The first carbon must be lower than the last.")
            return False
        n = len(self.points)
        if n > 2 and not (first < self.points[1][0] and last > self.points[-2][0]):
            self.inconsistent.emit("The end carbons must lie outside the dots between them.")
            return False
        if int(first) == int(self.points[0][0]) and int(last) == int(self.points[-1][0]):
            return True
        self.points[0][0], self.points[-1][0] = float(first), float(last)
        self._rebuild_handles()
        self._redraw()
        self.edited.emit(f"Set the ladder to C{first} to C{last}")
        return True

    def set_fine(self, on: bool) -> None:
        """A dot on every carbon, or only a few. Both keep the line where it is."""
        if len(self.points) < 2:
            return
        c = np.arange(int(self.points[0][0]), int(self.points[-1][0]) + 1)
        t1, t2 = sample_ladder(self.points, c)
        every = [(int(cc), float(a), float(b)) for cc, a, b in zip(c, t1, t2)]
        if on:
            chosen = every
        else:
            keep = np.unique(np.round(np.linspace(0, len(every) - 1, min(FEW_HANDLES, len(every)))).astype(int))
            chosen = [every[i] for i in keep]
        self.points = [list(p) for p in chosen]
        self._rebuild_handles()
        self._redraw()


class LadderPanel(QtWidgets.QDialog):
    """The small window that floats beside the 2D view while the ladder is being drawn."""

    def __init__(self, parent):
        super().__init__(parent)
        self.setWindowTitle("📏 Alkane ladder")
        self.setWindowFlag(QtCore.Qt.WindowType.Tool, True)
        fit_to_screen(self, 380, 400)
        self.title = QtWidgets.QLabel("")
        self.title.setWordWrap(True)
        self.title.setStyleSheet("font-weight: bold;")
        self.hint = QtWidgets.QLabel(HINT)
        self.hint.setWordWrap(True)
        self.hint.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.first = QtWidgets.QSpinBox()
        self.last = QtWidgets.QSpinBox()
        for box in (self.first, self.last):
            box.setRange(1, 100)
            box.setPrefix("C")
        self.first.setToolTip("The n-alkane under the first dot")
        self.last.setToolTip("The n-alkane under the last dot")
        self.fine = QtWidgets.QCheckBox("A dot on every carbon")
        self.fine.setToolTip("Off: a few dots, which move the ladder as a smooth curve. On: each carbon can be placed alone.")
        self.snap_btn = QtWidgets.QPushButton("🎯 Snap to peaks")
        self.snap_btn.setToolTip("Move each carbon onto the nearest real alkane peak (m/z 57 with the alkane ion ratios). "
                                 "Carbons with no such peak nearby stay where you put them.")
        self.back_btn = QtWidgets.QPushButton("↩ Back to the saved ladder")
        self.done_btn = QtWidgets.QPushButton("Done")
        self.done_btn.setObjectName("go")
        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)
        self.status.setTextFormat(QtCore.Qt.TextFormat.RichText)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("From"))
        row.addWidget(self.first)
        row.addWidget(QtWidgets.QLabel("to"))
        row.addWidget(self.last)
        row.addStretch(1)
        buttons = QtWidgets.QHBoxLayout()
        buttons.addWidget(self.snap_btn)
        buttons.addWidget(self.back_btn)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(self.title)
        lay.addWidget(self.hint)
        lay.addLayout(row)
        lay.addWidget(self.fine)
        lay.addLayout(buttons)
        lay.addWidget(self.status, 1)
        lay.addWidget(self.done_btn)


class LadderEditor(QtCore.QObject):
    """Ties the overlay, the panel, the stored ladder and the undo stack together for one run at a time."""

    def __init__(self, main_window):
        super().__init__(main_window)
        self.mw = main_window
        self.project = main_window.project
        self.hydro = main_window.project.hydrocarbons
        self.overlay = LadderOverlay(main_window.view2d, self)
        self.panel = LadderPanel(main_window)
        self.index = -1
        self.saved = None                                       # the stored ladder when last looked at
        self._updating = False
        self._internal = False                                  # True while this editor is the one changing the store
        self.overlay.edited.connect(self._edited)
        self.overlay.inconsistent.connect(lambda text: self.say(f"<span style='color:#D55E00'>{text}</span>"))
        self.panel.done_btn.clicked.connect(self.close)
        self.panel.snap_btn.clicked.connect(self.snap)
        self.panel.back_btn.clicked.connect(self.reload)
        self.panel.fine.toggled.connect(self._fine_toggled)
        self.panel.first.editingFinished.connect(self._ends_edited)
        self.panel.last.editingFinished.connect(self._ends_edited)
        self.panel.finished.connect(lambda _r: self.close())
        self.hydro.changed.connect(self._store_changed)

    @property
    def active(self) -> bool:
        return self.overlay.active

    def say(self, text: str) -> None:
        self.panel.status.setText(text)

    # ------------------------------------------------------------ opening and closing
    def open(self, index: int) -> None:
        """Draw run `index`'s ladder over the 2D view (a first, rough one if it has none)."""
        self.index = int(index)
        self.mw.act_aligned.setChecked(False)                  # the ladder is in this run's own coordinates
        self.mw.tabs.setCurrentIndex(0)
        self.overlay.show(True)
        self.reload()
        self.panel.show()
        self.panel.raise_()

    def close(self) -> None:
        self.overlay.show(False)
        self.panel.hide()
        self.index = -1

    def run_changed(self) -> None:
        """The main window moved to another run: edit that one's ladder instead."""
        if self.active and self.mw.index >= 0 and self.mw.index != self.index:
            self.open(self.mw.index)

    def reload(self) -> None:
        """Read the stored ladder into the overlay (or a rough starting ladder if the run has none)."""
        self._updating = True
        series = self.hydro.series_for(self.index)
        name = self.project.entries[self.index].name
        self.panel.title.setText(f"Alkane ladder of {name}")
        if series is not None and len(series) >= 2:
            self.saved = series.copy()
            # Every carbon gets a dot, so what is drawn is exactly what is stored; "A dot on every carbon" off simplifies it.
            self.overlay.load(control_points_of(series), [(h.t1, h.t2) for h in series.hits.values() if h.height > 0])
            carbons = series.carbons()
            self.say(f"{len(series)} carbons, C{int(carbons.min())} to C{int(carbons.max())}. Crosses mark apexes "
                     "found in the data; the white line is the ladder.")
        else:
            self.saved = None
            view = self.mw.view2d.plot.vb.viewRange()
            (x0, x1), (y0, y1) = view
            c0, c1 = DEFAULT_RANGE
            self.overlay.load([(c0, (x0 + 0.25 * (x1 - x0)) * 60.0, y0 + 0.25 * (y1 - y0)),
                               (c1, (x0 + 0.75 * (x1 - x0)) * 60.0, y0 + 0.55 * (y1 - y0))])
            self.say("This run has no ladder yet. Drag the dots onto your alkanes; it is saved as soon as you move one. "
                     "Then press <b>Snap to peaks</b>.")
        points = self.overlay.points
        self.panel.first.setValue(int(points[0][0]))
        self.panel.last.setValue(int(points[-1][0]))
        self.panel.fine.blockSignals(True)
        self.panel.fine.setChecked(len(points) > FEW_HANDLES)
        self.panel.fine.blockSignals(False)
        self.panel.snap_btn.setEnabled(not self.hydro.running)
        self._updating = False

    # ------------------------------------------------------------ editing
    def _edited(self, text: str) -> None:
        """A drag or an add/remove finished: store the ladder, as one undo step."""
        try:
            series = self.overlay.series()
        except ValueError as exc:
            self.say(f"<span style='color:#D55E00'>{exc}</span>")
            return
        self._store(series, text)

    def _store(self, series, text: str) -> None:
        before = self.saved
        self._internal = True
        try:
            self.hydro.set_series(self.index, series)
        finally:
            self._internal = False
        after = self.hydro.snapshot_series(self.index)
        self.saved = after.copy()
        self.mw.undo.push(EditLadder(self.hydro, self.index, before, after, f"📏 {text}"))
        carbons = after.carbons()
        self._updating = True
        self.panel.first.setValue(int(carbons.min()))
        self.panel.last.setValue(int(carbons.max()))
        self._updating = False
        self.say(f"{text}: saved (C{int(carbons.min())} to C{int(carbons.max())}, {len(after)} carbons). "
                 "Library matches made on the earlier ladder are now out of date: match again.")

    def _store_changed(self) -> None:
        """The stored ladder changed from outside (the alkane search, propagation, a different shape): follow it."""
        if self.active and not self._internal and 0 <= self.index < len(self.project.entries):
            now = self.hydro.series_for(self.index)
            if now is not self.saved and (now is None or self.saved is None or
                                          {c: (h.t1, h.t2) for c, h in now.hits.items()} !=
                                          {c: (h.t1, h.t2) for c, h in self.saved.hits.items()}):
                self.reload()

    def _fine_toggled(self, on: bool) -> None:
        if not self._updating:
            self.overlay.set_fine(on)

    def _ends_edited(self) -> None:
        if self._updating:
            return
        if not self.overlay.set_ends(self.panel.first.value(), self.panel.last.value()):
            self._updating = True
            self.panel.first.setValue(int(self.overlay.points[0][0]))
            self.panel.last.setValue(int(self.overlay.points[-1][0]))
            self._updating = False

    def snap(self) -> None:
        if self.saved is None:
            self.say("Move a dot first, so the ladder is saved; then snap it to the peaks.")
            return
        self.panel.snap_btn.setEnabled(False)
        self.say("Looking for the alkane peaks…")

        def done(series, snapped):
            self.panel.snap_btn.setEnabled(True)
            self._store(series, f"Snap {len(snapped)} of {len(series)} carbons to peaks")
            self.panel.fine.blockSignals(True)
            self.panel.fine.setChecked(True)
            self.panel.fine.blockSignals(False)
            self.reload()
            self.say(f"{len(snapped)} of {len(series)} carbons were moved onto a peak (crosses); the other "
                     f"{len(series) - len(snapped)} stay where you drew them. Check them against the image.")

        def failed(message):
            self.panel.snap_btn.setEnabled(True)
            self.say(f"<span style='color:#D55E00'>Could not snap the ladder: {message}. The ladder is as you drew it.</span>")

        if not self.hydro.snap(self.index, on_done=done, on_fail=failed):
            self.panel.snap_btn.setEnabled(True)
