"""📍 Alignment state for the app.

Holds which run is the reference, the anchor points collected for every
other run, and the transform fitted from them. The maths is all in
`core/align.py`; this just wires it to the project (for runs and folds)
and to background jobs (for batch anchor-finding).
"""

from __future__ import annotations

from PySide6 import QtCore

from ..core.align import AlignmentError, AnchorPoint, fit_alignment, propagate_anchors, resample_to_reference
from ..core.phase import Anchor
from .constants import EXTEND_FACTOR

DEFAULT_WINDOW = 12.0


class AlignState(QtCore.QObject):
    changed = QtCore.Signal()   # reference, an anchor, a transform, or a kind override changed

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.reference: str | None = None
        self.points: dict[str, list[AnchorPoint]] = {}
        self.kind: dict[str, str] = {}
        self.transforms: dict = {}
        self.errors: dict[str, str] = {}
        project.runsChanged.connect(self._prune)

    def _prune(self) -> None:
        names = {e.name for e in self.project.entries}
        if self.reference not in names:
            self.reference = None
        for d in (self.points, self.kind, self.transforms, self.errors):
            for k in [k for k in d if k not in names]:
                d.pop(k)

    # ------------------------------------------------------------- editing
    def snapshot(self) -> dict:
        return dict(reference=self.reference, kind=dict(self.kind),
                    points={k: list(v) for k, v in self.points.items()})

    def restore(self, snap: dict) -> None:
        self.reference = snap.get("reference")
        self.kind = dict(snap.get("kind", {}))
        self.points = {k: list(v) for k, v in snap.get("points", {}).items()}
        self._recompute_all()
        self.changed.emit()

    def set_reference(self, name: str | None) -> None:
        self.reference = name
        self.points, self.transforms, self.errors = {}, {}, {}
        self.changed.emit()

    def set_kind(self, run: str, kind: str) -> None:
        self.kind[run] = kind
        self._recompute(run)
        self.changed.emit()

    def add_anchor(self, run: str, name: str, ref_t1: float, ref_t2: float,
                   run_t1: float, run_t2: float, ion: int = 0) -> None:
        pts = self.points.setdefault(run, [])
        pts[:] = [p for p in pts if p.name != name]
        pts.append(AnchorPoint(name, float(ref_t1), float(ref_t2), float(run_t1), float(run_t2), int(ion)))
        self._recompute(run)
        self.changed.emit()

    def remove_anchor(self, run: str, name: str) -> None:
        pts = self.points.get(run, [])
        pts[:] = [p for p in pts if p.name != name]
        if not pts:
            self.points.pop(run, None)
            self.transforms.pop(run, None)
            self.errors.pop(run, None)
        else:
            self._recompute(run)
        self.changed.emit()

    def _recompute(self, run: str) -> None:
        pts = self.points.get(run)
        if not pts:
            self.transforms.pop(run, None)
            self.errors.pop(run, None)
            return
        try:
            self.transforms[run] = fit_alignment(pts, self.kind.get(run, "auto"))
            self.errors.pop(run, None)
        except AlignmentError as exc:
            self.transforms.pop(run, None)
            self.errors[run] = str(exc)

    def _recompute_all(self) -> None:
        self.transforms, self.errors = {}, {}
        for run in list(self.points):
            self._recompute(run)

    # ------------------------------------------------------------- queries
    def anchor_names(self) -> list[str]:
        names: list[str] = []
        for pts in self.points.values():
            for p in pts:
                if p.name not in names:
                    names.append(p.name)
        return names

    def next_anchor_name(self) -> str:
        existing = set(self.anchor_names())
        i = 1
        while f"P{i}" in existing:
            i += 1
        return f"P{i}"

    def reference_anchors(self) -> list[Anchor]:
        """One `Anchor` per distinct name that has a known ion, built from
        wherever it was first clicked (the reference position is the same
        wherever a name appears, by construction). Used for batch search."""
        seen: dict[str, Anchor] = {}
        for pts in self.points.values():
            for p in pts:
                if p.name not in seen and p.ion:
                    seen[p.name] = Anchor(p.name, p.ion, p.ref_t1, p.ref_t2, t1_window=DEFAULT_WINDOW)
        return list(seen.values())

    def aligned_runs(self) -> set[str]:
        """Names usable as-is: the reference itself, plus every run with a
        transform (whatever its own errors dict says about *other* runs)."""
        out = set(self.transforms)
        if self.reference:
            out.add(self.reference)
        return out

    def image_for(self, index: int, expr, extended: bool = False):
        """This run's ion image, resampled onto the reference's grid. None
        if there's no reference, or this run (not the reference) has no
        usable transform yet."""
        p = self.project
        if self.reference is None or not (0 <= index < len(p.entries)):
            return None
        entry = p.entries[index]
        ref_i = p.index_of(self.reference)
        if ref_i is None or entry.run is None or p.entries[ref_i].run is None:
            return None
        ref_fold = p.fold(ref_i)
        if entry.name == self.reference:
            img = expr.image(p.display_run(index), ref_fold)
            if extended:
                img, _ = ref_fold.extend(img, EXTEND_FACTOR)
            return img
        tf = self.transforms.get(entry.name)
        if tf is None:
            return None
        run_fold = p.fold(index)
        img = expr.image(p.display_run(index), run_fold)
        return resample_to_reference(img, run_fold, ref_fold, tf, extended=extended)

    # ------------------------------------------------------------- batch
    def propagate_to(self, indices: list[int], on_done=None) -> None:
        anchors = self.reference_anchors()
        if not anchors:
            self.project.message.emit("No named anchors with a known ion to search for yet")
            return
        entries = [self.project.entries[i] for i in indices
                  if self.project.entries[i].run is not None and self.project.entries[i].name != self.reference]
        if not entries:
            return
        remaining = set(entries)
        failed: list[str] = []

        def start(entry) -> None:
            idx = self.project.entries.index(entry)
            run, fold = entry.run, self.project.fold(idx)

            def work(token, progress):
                return propagate_anchors(anchors, run, fold)

            def done(found) -> None:
                missing = [n for n, v in found.items() if v is None]
                for name, pt in found.items():
                    if pt is not None:
                        pts = self.points.setdefault(entry.name, [])
                        pts[:] = [p for p in pts if p.name != name] + [pt]
                if missing:
                    failed.append(f"{entry.name} ({', '.join(missing)})")
                self._recompute(entry.name)
                remaining.discard(entry)
                if not remaining:
                    self.changed.emit()
                    if failed:
                        self.project.message.emit("Couldn't find every anchor in: " + "; ".join(failed))
                    if on_done:
                        on_done()

            self.project.jobs.submit(f"📍 Finding anchors: {entry.name}", work, done)

        for e in entries:
            start(e)

    # ------------------------------------------------------------- files
    def to_dict(self) -> dict:
        return dict(reference=self.reference, kind=dict(self.kind),
                    points={run: [dict(name=p.name, ref_t1=p.ref_t1, ref_t2=p.ref_t2,
                                       run_t1=p.run_t1, run_t2=p.run_t2, ion=p.ion) for p in pts]
                           for run, pts in self.points.items()})

    def from_dict(self, d: dict) -> None:
        d = d or {}
        self.reference = d.get("reference")
        self.kind = dict(d.get("kind", {}))
        self.points = {run: [AnchorPoint(**p) for p in pts] for run, pts in d.get("points", {}).items()}
        self._recompute_all()
        self.changed.emit()
