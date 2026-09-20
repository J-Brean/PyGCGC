"""🧪 Review gallery.

One tile per peak: the quantifier-ion image around it, its outline, the
amount and any flags, with a coloured border (green fine, amber to check,
red rejected, grey not found). Show one compound across every run, or
every compound in one run.

Keys: arrows move, Y accepts, N rejects, U clears the decision,
Enter (or double-click) opens the peak in the 2D view, P picks a
different peak for it by clicking in the 2D view.
"""

from __future__ import annotations

import math

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from .results_colours import STATUS_QCOLOURS

TILE_W, TILE_H = 220, 170
IMG_H = 110


class _Tile:
    def __init__(self, run, compound, row, hit, image, outline, apex):
        self.run, self.compound, self.row, self.hit = run, compound, row, hit
        self.image, self.outline, self.apex = image, outline, apex


class ReviewView(QtWidgets.QWidget):
    gotoRequested = QtCore.Signal(str, str)               # run, compound
    decisionRequested = QtCore.Signal(list, str)          # [(run, compound)], status
    pickRequested = QtCore.Signal(str, str)               # run, compound
    reintegrateRequested = QtCore.Signal(str, str)        # run, compound (key I)

    def __init__(self, image_source, parent=None):
        """image_source(run_name, ion) -> (extended image, Frame) or None."""
        super().__init__(parent)
        self.image_source = image_source
        self.report = None
        self.tiles: list[_Tile] = []
        self.current = 0
        self._lut = pg.colormap.get("magma").getLookupTable(nPts=256, alpha=False)

        self.mode = QtWidgets.QComboBox()
        self.mode.addItem("One compound, every run", "compound")
        self.mode.addItem("One run, every compound", "run")
        self.which = QtWidgets.QComboBox()
        self.only_check = QtWidgets.QCheckBox("Only to check")
        self.info = QtWidgets.QLabel("")
        help_ = QtWidgets.QLabel("Y accept, N reject, U clear, Enter open, P pick another peak")
        help_.setStyleSheet("color: palette(mid);")
        self.mode.currentIndexChanged.connect(self._fill_which)
        self.which.currentIndexChanged.connect(self.rebuild)
        self.only_check.toggled.connect(self.rebuild)

        self.canvas = _Canvas(self)
        self.scroll = QtWidgets.QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.canvas)
        top = QtWidgets.QHBoxLayout()
        for w in (self.mode, self.which, self.only_check):
            top.addWidget(w)
        top.addStretch(1)
        top.addWidget(help_)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addLayout(top)
        lay.addWidget(self.scroll, 1)
        lay.addWidget(self.info)
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)

    # ------------------------------------------------------------------ data
    def set_report(self, report) -> None:
        self.report = report
        self._fill_which()

    def _fill_which(self) -> None:
        cur = self.which.currentData()
        self.which.blockSignals(True)
        self.which.clear()
        if self.report is not None:
            if self.mode.currentData() == "compound":
                for t in self.report.targets:
                    self.which.addItem(t.compound + ("  (IS)" if t.is_istd else ""), t.compound)
            else:
                for r in self.report.seq:
                    if r.quant is not None:
                        self.which.addItem(r.name, r.name)
        i = self.which.findData(cur)
        self.which.setCurrentIndex(max(i, 0))
        self.which.blockSignals(False)
        self.rebuild()

    def show_compound(self, compound: str) -> None:
        self.mode.setCurrentIndex(0)
        i = self.which.findData(compound)
        if i >= 0:
            self.which.setCurrentIndex(i)

    def rebuild(self) -> None:
        keep = self.tiles[self.current].run if self.tiles and 0 <= self.current < len(self.tiles) else None
        keep_c = self.tiles[self.current].compound if keep else None
        self.tiles = []
        r = self.report
        if r is not None and self.which.count():
            key = self.which.currentData()
            rows = {(x.run, x.compound): x for x in r.results}
            targets = {t.compound: t for t in r.targets}
            if self.mode.currentData() == "compound":
                pairs = [(s.name, key) for s in r.seq if s.quant is not None]
            else:
                pairs = [(key, t.compound) for t in r.targets]
            for run, comp in pairs:
                row = rows.get((run, comp))
                if row is None:
                    continue
                if self.only_check.isChecked() and row.status not in ("check", "not found"):
                    continue
                q = next(s.quant for s in r.seq if s.name == run)
                self.tiles.append(self._tile(run, comp, row, q.hits.get(comp), targets[comp]))
        self.current = 0
        for i, t in enumerate(self.tiles):
            if (t.run, t.compound) == (keep, keep_c):
                self.current = i
        self.canvas.relayout()
        self._update_info()

    def _tile(self, run, comp, row, hit, target) -> _Tile:
        src = self.image_source(run, target.quant_ion)
        if src is None:
            return _Tile(run, comp, row, hit, None, None, None)
        img, frame = src
        if hit is not None and hit.found:
            r0, c0 = hit.apex
        else:
            c0, rp = frame.time_to_pixel(*(hit.expected if hit is not None and hit.expected
                                           else (target.t1, target.t2)))
            r0 = int(rp)
        half_r = max(int(round(0.6 / frame.pitch)), 10)
        half_c = 6
        rs, cs = slice(max(r0 - half_r, 0), r0 + half_r + 1), slice(max(c0 - half_c, 0), c0 + half_c + 1)
        crop = img[rs, cs]
        if crop.size == 0:
            return _Tile(run, comp, row, hit, None, None, None)
        v = np.sqrt(np.clip(crop, 0, None))
        top = v.max() if v.max() > 0 else 1.0
        idx = np.clip((v / top) * 255, 0, 255).astype(np.uint8)
        rgb = self._lut[idx][::-1]      # flip so t2 increases upwards
        rgb = np.ascontiguousarray(rgb)
        h, w, _ = rgb.shape
        qimg = QtGui.QImage(rgb.data, w, h, 3 * w, QtGui.QImage.Format.Format_RGB888).copy()
        outline = None
        if hit is not None and hit.found and len(hit.pixels) == 2:
            mask = np.zeros(crop.shape, bool)
            rr, cc = hit.pixels
            sel = (rr >= rs.start) & (rr < rs.start + h) & (cc >= cs.start) & (cc < cs.start + w)
            mask[rr[sel] - rs.start, cc[sel] - cs.start] = True
            outline = mask[::-1]
        apex = (h - 1 - (r0 - rs.start), c0 - cs.start)
        return _Tile(run, comp, row, hit, qimg, outline, apex)

    # ---------------------------------------------------------------- keys
    def keyPressEvent(self, ev) -> None:
        k = ev.key()
        n = len(self.tiles)
        if not n:
            return super().keyPressEvent(ev)
        cols = self.canvas.columns()
        moves = {QtCore.Qt.Key.Key_Right: 1, QtCore.Qt.Key.Key_Left: -1,
                 QtCore.Qt.Key.Key_Down: cols, QtCore.Qt.Key.Key_Up: -cols}
        t = self.tiles[self.current]
        if k in moves:
            self.set_current(self.current + moves[k])
        elif k == QtCore.Qt.Key.Key_Y:
            self.decisionRequested.emit([(t.run, t.compound)], "accepted")
            self.set_current(self.current + 1)
        elif k == QtCore.Qt.Key.Key_N:
            self.decisionRequested.emit([(t.run, t.compound)], "rejected")
            self.set_current(self.current + 1)
        elif k == QtCore.Qt.Key.Key_U:
            self.decisionRequested.emit([(t.run, t.compound)], "clear")
        elif k in (QtCore.Qt.Key.Key_Return, QtCore.Qt.Key.Key_Enter):
            self.gotoRequested.emit(t.run, t.compound)
        elif k == QtCore.Qt.Key.Key_P:
            self.pickRequested.emit(t.run, t.compound)
        elif k == QtCore.Qt.Key.Key_I:
            self.reintegrateRequested.emit(t.run, t.compound)
        else:
            super().keyPressEvent(ev)

    def set_current(self, i: int) -> None:
        if not self.tiles:
            return
        self.current = int(np.clip(i, 0, len(self.tiles) - 1))
        rect = self.canvas.tile_rect(self.current)
        self.scroll.ensureVisible(rect.center().x(), rect.center().y(), TILE_W // 2, TILE_H // 2)
        self.canvas.update()
        self._update_info()

    def _update_info(self) -> None:
        if not self.tiles:
            self.info.setText("Nothing to review" if self.report else "Run 🧪 Quantify first")
            return
        t = self.tiles[self.current]
        n_check = sum(x.row.status == "check" for x in self.tiles)
        self.info.setText(f"{self.current + 1} of {len(self.tiles)}   ({n_check} to check)   "
                          f"{t.run}, {t.compound}: {t.row.status}"
                          + (f", {', '.join(t.row.flags)}" if t.row.flags else ""))


class _Canvas(QtWidgets.QWidget):
    def __init__(self, view: ReviewView):
        super().__init__()
        self.view = view
        self.setMouseTracking(False)

    def columns(self) -> int:
        return max(1, self.view.scroll.viewport().width() // TILE_W)

    def relayout(self) -> None:
        n = len(self.view.tiles)
        rows = math.ceil(n / self.columns()) if n else 0
        self.setMinimumHeight(rows * TILE_H)
        self.update()

    def resizeEvent(self, ev) -> None:
        self.relayout()
        super().resizeEvent(ev)

    def tile_rect(self, i: int) -> QtCore.QRect:
        cols = self.columns()
        return QtCore.QRect((i % cols) * TILE_W, (i // cols) * TILE_H, TILE_W, TILE_H)

    def tile_at(self, pos) -> int | None:
        cols = self.columns()
        c, r = pos.x() // TILE_W, pos.y() // TILE_H
        i = r * cols + c
        return i if (c < cols and 0 <= i < len(self.view.tiles)) else None

    def mousePressEvent(self, ev) -> None:
        i = self.tile_at(ev.position().toPoint())
        if i is not None:
            self.view.set_current(i)
            self.view.setFocus()

    def mouseDoubleClickEvent(self, ev) -> None:
        i = self.tile_at(ev.position().toPoint())
        if i is not None:
            t = self.view.tiles[i]
            self.view.gotoRequested.emit(t.run, t.compound)

    def paintEvent(self, ev) -> None:
        p = QtGui.QPainter(self)
        pal = self.palette()
        text_col = pal.color(QtGui.QPalette.ColorRole.WindowText)
        mid = pal.color(QtGui.QPalette.ColorRole.Mid)
        for i, t in enumerate(self.view.tiles):
            rect = self.tile_rect(i)
            if not ev.rect().intersects(rect):
                continue
            inner = rect.adjusted(6, 6, -6, -6)
            colour = STATUS_QCOLOURS.get(t.row.status, mid)
            pen = QtGui.QPen(colour, 4 if i == self.view.current else 2)
            p.setPen(pen)
            p.drawRoundedRect(inner, 6, 6)
            img_rect = QtCore.QRect(inner.left() + 6, inner.top() + 6, inner.width() - 12, IMG_H)
            if t.image is not None:
                p.drawImage(img_rect, t.image)
                sx = img_rect.width() / t.image.width()
                sy = img_rect.height() / t.image.height()
                if t.outline is not None:
                    p.setPen(QtGui.QPen(QtGui.QColor(126, 232, 250), 1))
                    m = t.outline
                    ys, xs = np.nonzero(m & ~np.pad(m, 1)[2:, 1:-1])     # bottom edges
                    for y, x in zip(ys, xs):
                        p.drawLine(QtCore.QPointF(img_rect.left() + x * sx, img_rect.top() + (y + 1) * sy),
                                   QtCore.QPointF(img_rect.left() + (x + 1) * sx, img_rect.top() + (y + 1) * sy))
                    ys, xs = np.nonzero(m & ~np.pad(m, 1)[:-2, 1:-1])    # top edges
                    for y, x in zip(ys, xs):
                        p.drawLine(QtCore.QPointF(img_rect.left() + x * sx, img_rect.top() + y * sy),
                                   QtCore.QPointF(img_rect.left() + (x + 1) * sx, img_rect.top() + y * sy))
                    ys, xs = np.nonzero(m & ~np.pad(m, 1)[1:-1, :-2])    # left edges
                    for y, x in zip(ys, xs):
                        p.drawLine(QtCore.QPointF(img_rect.left() + x * sx, img_rect.top() + y * sy),
                                   QtCore.QPointF(img_rect.left() + x * sx, img_rect.top() + (y + 1) * sy))
                    ys, xs = np.nonzero(m & ~np.pad(m, 1)[1:-1, 2:])     # right edges
                    for y, x in zip(ys, xs):
                        p.drawLine(QtCore.QPointF(img_rect.left() + (x + 1) * sx, img_rect.top() + y * sy),
                                   QtCore.QPointF(img_rect.left() + (x + 1) * sx, img_rect.top() + (y + 1) * sy))
                if t.apex is not None and (t.hit is None or not t.hit.found):
                    ay, ax = t.apex
                    p.setPen(QtGui.QPen(QtGui.QColor(255, 255, 255, 160), 1, QtCore.Qt.PenStyle.DashLine))
                    cx, cy = img_rect.left() + (ax + 0.5) * sx, img_rect.top() + (ay + 0.5) * sy
                    p.drawEllipse(QtCore.QPointF(cx, cy), 8, 8)
            else:
                p.setPen(mid)
                p.drawText(img_rect, QtCore.Qt.AlignmentFlag.AlignCenter, "no image")
            p.setPen(text_col)
            title = t.run if self.view.mode.currentData() == "compound" else t.compound
            f = p.font()
            f.setBold(True)
            p.setFont(f)
            p.drawText(QtCore.QRect(inner.left() + 6, img_rect.bottom() + 2, inner.width() - 12, 16),
                       QtCore.Qt.AlignmentFlag.AlignLeft, _elide(title, 30))
            f.setBold(False)
            p.setFont(f)
            amt = t.row.amount
            line = f"{amt:.3g} ± {t.row.amount_se:.2g}" if np.isfinite(amt) else (
                f"response {t.row.response:.3g}" if np.isfinite(t.row.response) else t.row.status)
            if t.hit is not None and t.hit.manual:
                line += "  (picked by hand)"
            p.drawText(QtCore.QRect(inner.left() + 6, img_rect.bottom() + 18, inner.width() - 12, 16),
                       QtCore.Qt.AlignmentFlag.AlignLeft, line)
            p.setPen(colour)
            flags = ", ".join(t.row.flags) if t.row.flags else t.row.status
            p.drawText(QtCore.QRect(inner.left() + 6, img_rect.bottom() + 33, inner.width() - 12, 16),
                       QtCore.Qt.AlignmentFlag.AlignLeft, _elide(flags, 34))
        p.end()


def _elide(s: str, n: int) -> str:
    return s if len(s) <= n else s[: n - 1] + "…"
