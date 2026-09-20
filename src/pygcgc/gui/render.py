"""Builds what the chromatogram views show, off the main thread.

The views say what they want (which run, which ion expression, extended
or not); the renderer checks its cache, otherwise folds the data in a
background job, and announces the result. Only the newest request counts,
so flicking quickly through runs never shows a stale image. After each
image it quietly prepares the neighbouring runs so the next flick is
instant.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
from PySide6 import QtCore

from ..core import IonExpr, IonExpressionError, parse_ion_expr
from .constants import EXTEND_FACTOR
from .jobs import JobManager
from .project import Project


@dataclass
class Rendered:
    index: int
    name: str
    expr: IonExpr
    extended: bool
    image: np.ndarray        # (rows, modulations)
    t1: np.ndarray           # modulation start times (s)
    t2: np.ndarray           # row times (s)
    pitch: float
    period: float
    phase: float
    n_rows: int              # rows in one period (image may be taller if extended)
    times: np.ndarray        # scan times (s)
    signal: np.ndarray       # per-scan signal behind the image (numerator for ratios)
    processed: bool = False  # showing a 🧹 layer (or its preview)
    preview: bool = False
    snr: bool = False        # image is signal-to-noise
    image_before: np.ndarray | None = None   # raw image for the before/after split
    aligned: bool = False     # image is on the reference run's grid, not this run's own

    @property
    def profile(self) -> np.ndarray:
        """Reconstructed 1D trace: each modulation summed down t2 (one period only)."""
        img = self.image[: self.n_rows]
        return np.nansum(img, axis=0)


class Renderer(QtCore.QObject):
    ready = QtCore.Signal(object)       # Rendered
    failed = QtCore.Signal(str)

    def __init__(self, project: Project, jobs: JobManager, cache_size: int = 48, parent=None):
        super().__init__(parent)
        self.project = project
        self.jobs = jobs
        self.cache: OrderedDict[tuple, Rendered] = OrderedDict()
        self.cache_size = cache_size
        self.index: int = -1
        self.expr: IonExpr = IonExpr(None)
        self.extended: bool = False
        self.split: bool = False
        self.snr: bool = False
        self.aligned: bool = False
        self.align = None
        self.preload = True
        self._pending: set[tuple] = set()
        self._timer = QtCore.QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(15)
        self._timer.timeout.connect(self._render)
        project.foldChanged.connect(self._fold_changed)
        project.runsChanged.connect(self._runs_changed)
        project.layersChanged.connect(self.request)

    def set_align(self, align) -> None:
        self.align = align
        align.changed.connect(self.request)

    def set_aligned(self, on: bool) -> None:
        self.aligned = bool(on)
        self.request()

    def aligned_available(self, index: int | None = None) -> bool:
        i = self.index if index is None else index
        a = self.align
        if a is None or not a.reference or not (0 <= i < len(self.project.entries)):
            return False
        e = self.project.entries[i]
        return e.run is not None and (e.name == a.reference or e.name in a.transforms)

    # ------------------------------------------------------------- requests
    def set_index(self, index: int) -> None:
        self.index = index
        self.request()

    def set_expression(self, text: str) -> IonExpr:
        expr = parse_ion_expr(text)
        e = self._entry()
        if e is not None and e.run is not None:
            try:
                expr.validate(e.run)
            except KeyError as exc:
                raise IonExpressionError(str(exc).strip("'\"")) from None
        self.expr = expr
        self.request()
        return expr

    def set_extended(self, on: bool) -> None:
        self.extended = bool(on)
        self.request()

    def _use_aligned(self, index: int) -> bool:
        return self.aligned and self.aligned_available(index)

    def set_split(self, on: bool) -> None:
        self.split = bool(on)
        self.request()

    def set_snr(self, on: bool) -> None:
        self.snr = bool(on)
        self.request()

    def snr_available(self, index: int | None = None) -> bool:
        i = self.index if index is None else index
        return self._entry(i) is not None and self.project.result_for(i) is not None

    def request(self) -> None:
        self._timer.start()

    # ------------------------------------------------------------- internals
    def _entry(self, index: int | None = None):
        i = self.index if index is None else index
        if 0 <= i < len(self.project.entries):
            return self.project.entries[i]
        return None

    def _key(self, index: int):
        e = self._entry(index)
        if e is None or e.run is None:
            return None
        res = self.project.result_for(index)
        use_aligned = self._use_aligned(index)
        split = self.split and res is not None and not use_aligned
        snr = self.snr and e.background is not None and not self.expr.is_ratio and not use_aligned
        ref_id = id(self.align.transforms.get(e.name)) if (use_aligned and self.align) else None
        return (id(e.run), id(res) if res is not None else None, self.project.period,
                round(e.phase, 9), self.expr.text(), self.expr.is_ratio, self.extended, split, snr,
                use_aligned, self.align.reference if use_aligned else None, ref_id)

    def _fold_changed(self, index: int) -> None:
        if index == -1 or index == self.index:
            self.request()

    def _runs_changed(self) -> None:
        if self.index >= len(self.project.entries):
            self.index = len(self.project.entries) - 1
        self.request()

    def _remember(self, key, result: Rendered) -> None:
        self.cache[key] = result
        self.cache.move_to_end(key)
        while len(self.cache) > self.cache_size:
            self.cache.popitem(last=False)

    def _compute(self, index: int):
        """Return a function that builds the image for run `index` in a worker."""
        e = self._entry(index)
        raw = e.run
        res = self.project.result_for(index)
        run = res.run if res is not None else raw
        preview = self.project.preview is not None and self.project.preview[0] is e
        fold = self.project.fold(index)
        expr, extended = self.expr, self.extended
        use_aligned = self._use_aligned(index)
        split = self.split and res is not None and not use_aligned
        # S/N always comes from 🧹's own noise map, whether or not 🫧 also
        # ran on top of it (🫧 doesn't compute one of its own).
        noise = e.background.noise if e.background is not None else None
        snr = self.snr and noise is not None and not expr.is_ratio and not use_aligned
        period, phase, name = self.project.period, e.phase, e.name

        if use_aligned:
            align = self.align
            ref_i = self.project.index_of(align.reference)
            ref_fold = self.project.fold(ref_i)

            def work(token, progress):
                img = align.image_for(index, expr, extended=extended)
                token.check()
                t2 = ref_fold.t2
                if extended:
                    _, t2 = ref_fold.extend(np.zeros((ref_fold.n_rows, 1)), EXTEND_FACTOR)
                signal = expr.signal(run, expr.numerator)
                return Rendered(index=index, name=name, expr=expr, extended=extended,
                                image=img, t1=ref_fold.t1, t2=t2, pitch=ref_fold.row_pitch,
                                period=ref_fold.spec.period, phase=phase, n_rows=ref_fold.n_rows,
                                times=raw.times, signal=signal, processed=res is not None,
                                preview=preview, aligned=True)
            return work

        def build(r):
            img = expr.image(r, fold)
            if snr:
                sigma = noise.for_ions(expr.numerator, fold.t1)
                img = img / sigma[None, :]
            if extended:
                img, _ = fold.extend(img, EXTEND_FACTOR)
            return img

        def work(token, progress):
            img = build(run)
            token.check()
            before = build(raw) if split else None
            t2 = fold.extend(np.zeros((fold.n_rows, 1)), EXTEND_FACTOR)[1] if extended else fold.t2
            signal = expr.signal(run, expr.numerator)
            return Rendered(index=index, name=name, expr=expr, extended=extended,
                            image=img, t1=fold.t1, t2=t2, pitch=fold.row_pitch,
                            period=period, phase=phase, n_rows=fold.n_rows,
                            times=run.times, signal=signal, processed=res is not None,
                            preview=preview, snr=snr, image_before=before)
        return work

    def _render(self) -> None:
        key = self._key(self.index)
        if key is None:
            return
        hit = self.cache.get(key)
        if hit is not None:
            self.cache.move_to_end(key)
            hit.index = self.index
            self.ready.emit(hit)
            self._preload()
            return
        index = self.index

        def done(result: Rendered):
            self._remember(key, result)
            if index == self.index and self._key(index) == key:
                self.ready.emit(result)
                self._preload()

        self.jobs.submit("Drawing", self._compute(index), done, self.failed.emit,
                         latest_key="render", quiet=True)

    def _preload(self) -> None:
        if not self.preload:
            return
        for j in (self.index + 1, self.index - 1):
            key = self._key(j)
            if key is None or key in self.cache or key in self._pending:
                continue
            self._pending.add(key)

            def done(result, key=key):
                self._pending.discard(key)
                self._remember(key, result)

            def err(_text, key=key):
                self._pending.discard(key)

            self.jobs.submit("Preloading", self._compute(j), done, err, quiet=True)
