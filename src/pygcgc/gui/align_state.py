"""📍 Alignment state for the app.

Holds which run is the reference, the anchor points collected for every
other run, and the transform fitted from them. The maths is all in
`core/align.py`; this just wires it to the project (for runs and folds)
and to background jobs (for batch anchor-finding).
"""

from __future__ import annotations

from dataclasses import replace

from PySide6 import QtCore

from ..core.align import (AlignmentError, AnchorPoint, fit_alignment, propagate_anchors,
                          refine_anchor_position, resample_to_reference)
from ..core.phase import Anchor
from .constants import EXTEND_FACTOR

DEFAULT_WINDOW = 12.0


class AlignState(QtCore.QObject):
    changed = QtCore.Signal()   # reference, an anchor, a transform, or a kind override changed
    refinementRunning = QtCore.Signal(bool)
    phaseStatusChanged = QtCore.Signal()    # a run's phase was ticked as checked, or changed since it was anchored

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.reference: str | None = None
        self.points: dict[str, list[AnchorPoint]] = {}
        self.kind: dict[str, str] = {}
        self.transforms: dict = {}
        self.errors: dict[str, str] = {}
        self.refining = False
        # While anchors are being picked, clicks are recorded but the transforms are left alone until the
        # user presses Anchor: a transform refitted after every click would shift the images mid-pick and
        # invalidate any blank-subtracted layer each time. `applied_points` is what the transforms were fitted from.
        self.holding = False
        self.applied_points: dict[str, list[AnchorPoint]] = {}
        self.phase_at_anchor: dict[str, float] = {}      # each run's phase when its transform was fitted
        self.phase_checked: set[str] = set()             # runs whose phase the user has confirmed
        self._revision = 0
        for signal in (self.changed, project.runsChanged, project.foldChanged, project.layersChanged):
            signal.connect(self._invalidate_refinement)
        project.runsChanged.connect(self._prune)
        project.foldChanged.connect(self._phase_may_have_changed)

    def _invalidate_refinement(self, *_args):
        self._revision += 1

    def _prune(self) -> None:
        names = {e.name for e in self.project.entries}
        if self.reference not in names:
            self.reference = None
        for d in (self.points, self.kind, self.transforms, self.errors, self.applied_points, self.phase_at_anchor):
            for k in [k for k in d if k not in names]:
                d.pop(k)
        self.phase_checked &= names

    # ------------------------------------------------------------- phase
    def _entry(self, name: str):
        return next((e for e in self.project.entries if e.name == name), None)

    def _phase_may_have_changed(self, index: int = -1) -> None:
        """A run's phase or period changed: its confirmation no longer stands."""
        entries = self.project.entries
        names = [e.name for e in entries] if not (0 <= index < len(entries)) else [entries[index].name]
        before = set(self.phase_checked)
        self.phase_checked -= set(names)
        if self.phase_checked != before or self.stale_phase_runs():
            self.phaseStatusChanged.emit()

    def set_phase_checked(self, name: str, checked: bool) -> None:
        (self.phase_checked.add if checked else self.phase_checked.discard)(name)
        self.phaseStatusChanged.emit()

    def stale_phase_runs(self) -> list[str]:
        """Runs whose phase differs from what it was when their transform was fitted. Their anchors
        record positions inside a modulation, so a changed phase makes them wrong."""
        out = []
        for name, phase in self.phase_at_anchor.items():
            entry = self._entry(name)
            if entry is None or name not in self.transforms:
                continue
            period = float(self.project.period) if getattr(self.project, "period", None) else 0.0
            delta = abs(float(entry.phase) - phase)
            if period:
                delta = min(delta % period, period - delta % period)
            if delta > 1e-3:
                out.append(name)
        return out

    # ------------------------------------------------------------- holding
    def begin_hold(self) -> None:
        """Start recording clicks without refitting."""
        self.holding = True

    def pending_runs(self) -> list[str]:
        """Runs whose clicks differ from the ones their transform was fitted from."""
        return [run for run in set(self.points) | set(self.applied_points)
                if run != self.reference and self.points.get(run, []) != self.applied_points.get(run, [])]

    def apply_hold(self) -> None:
        """The Anchor button: fit every run's transform from its clicks."""
        self.holding = False
        self._recompute_all()
        self.changed.emit()

    # ------------------------------------------------------------- editing
    def snapshot(self) -> dict:
        return dict(reference=self.reference, kind=dict(self.kind),
                    points={k: list(v) for k, v in self.points.items()})

    def restore(self, snap: dict) -> None:
        self.holding = False
        self.reference = snap.get("reference")
        self.kind = dict(snap.get("kind", {}))
        self.points = {k: list(v) for k, v in snap.get("points", {}).items()}
        self._recompute_all()
        self.changed.emit()

    def set_reference(self, name: str | None) -> None:
        self.holding = False
        self.reference = name
        self.points, self.transforms, self.errors = {}, {}, {}
        self.applied_points, self.phase_at_anchor = {}, {}
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
        if self.holding:
            return
        pts = self.points.get(run)
        if not pts:
            self.transforms.pop(run, None)
            self.errors.pop(run, None)
            self.applied_points.pop(run, None)
            self.phase_at_anchor.pop(run, None)
            return
        self.applied_points[run] = list(pts)
        entry = self._entry(run)
        if entry is not None:
            self.phase_at_anchor[run] = float(entry.phase)
        try:
            self.transforms[run] = fit_alignment(pts, self.kind.get(run, "auto"))
            self.errors.pop(run, None)
        except AlignmentError as exc:
            self.transforms.pop(run, None)
            self.errors[run] = str(exc)

    def _recompute_all(self) -> None:
        self.transforms, self.errors, self.applied_points, self.phase_at_anchor = {}, {}, {}, {}
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
    def refine_all(self, columns: int = 1, t2_radius: float = .15, on_done=None) -> None:
        """Refine saved clicks in a job, then apply all runs as one undoable edit.

        Uses each anchor's saved ion, or TIC for ion 0, on the current processing
        layer. Reference positions are searched once and shared by every pair.
        Input changes while the job runs invalidate its result.
        """
        if self.refining or not self.reference or not any(self.points.values()):
            return
        before, revision = self.snapshot(), self._revision
        sources = {e.name: (self.project.display_run(i), self.project.fold(i))
                   for i, e in enumerate(self.project.entries) if e.run is not None}
        reference = self.reference
        if reference not in sources:
            return
        self.refining = True
        self.refinementRunning.emit(True)

        def work(token, progress):
            originals = {}
            # Prefer explicit reference clicks; old projects only stored pairs.
            for run in [reference] + [r for r in before["points"] if r != reference]:
                for p in before["points"].get(run, []):
                    originals.setdefault(p.name, (p.ref_t1, p.ref_t2, p.ion))
            missing, moved = [], 0

            def search(run, specs):
                nonlocal moved
                found, images = {}, {}
                if run not in sources:
                    missing.extend(f"{run}/{name}" for name in specs)
                    return found
                data, folded = sources[run]
                for name, (t1, t2, ion) in specs.items():
                    token.check()
                    try:
                        if ion not in images:
                            images[ion] = folded.image(data.xic(ion) if ion else data.tic)
                        point = refine_anchor_position(images[ion], folded, t1, t2, columns, t2_radius)
                    except KeyError:
                        point = None
                    if point is None:
                        missing.append(f"{run}/{name}")
                    else:
                        found[name] = point
                        moved += int(abs(point[0]-t1) > 1e-9 or abs(point[1]-t2) > 1e-9)
                return found

            ref_points = search(reference, originals)
            after = dict(before, points={})
            runs = list(before["points"])
            for number, run in enumerate(runs):
                progress(number, len(runs), f"Refining anchors: {run}")
                points = before["points"][run]
                found = ref_points if run == reference else search(
                    run, {p.name: (p.run_t1, p.run_t2, p.ion) for p in points})
                after["points"][run] = [replace(
                    p, ref_t1=ref_points.get(p.name, (p.ref_t1, p.ref_t2))[0],
                    ref_t2=ref_points.get(p.name, (p.ref_t1, p.ref_t2))[1],
                    run_t1=found.get(p.name, (p.run_t1, p.run_t2))[0],
                    run_t2=found.get(p.name, (p.run_t1, p.run_t2))[1]) for p in points]
            return after, moved, missing

        def finish():
            self.refining = False
            self.refinementRunning.emit(False)

        def done(result):
            finish()
            if revision != self._revision:
                self.project.message.emit("Anchor refinement discarded because the inputs changed. Run it again.")
                return
            after, moved, missing = result
            self.restore(after)
            if on_done:
                on_done(before, after)
            text = f"Refined anchors: {moved} positions adjusted"
            if missing:
                text += f"; {len(missing)} without a supported maximum, kept unchanged: " + ", ".join(missing)
            self.project.message.emit(text)

        def failed(error):
            finish()
            self.project.message.emit(f"Anchor refinement failed: {error}")

        self.project.jobs.submit("Refine anchor positions", work, done, failed, on_cancel=finish)

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
        before = self.snapshot()
        inputs = self.project.analysis_inputs()
        found_by_run = {}
        unsuccessful = set()

        def settled(entry, success):
            if not success:
                unsuccessful.add(entry)
            remaining.discard(entry)
            if remaining:
                return
            if unsuccessful:
                return
            if inputs != self.project.analysis_inputs() or self.snapshot() != before:
                self.project.message.emit("Anchor search discarded: inputs changed; search again")
                return
            points = {run: list(pts) for run, pts in before["points"].items()}
            for run, found in found_by_run.items():
                for name, point in found.items():
                    if point is not None:
                        points[run] = [p for p in points.get(run, []) if p.name != name] + [point]
            self.restore(dict(before, points=points))
            if failed:
                self.project.message.emit("Couldn't find every anchor in: " + "; ".join(failed))
            if on_done:
                on_done()

        def start(entry) -> None:
            idx = self.project.entries.index(entry)
            run, fold = entry.run, self.project.fold(idx)

            def work(token, progress):
                return propagate_anchors(anchors, run, fold)

            def done(found) -> None:
                missing = [n for n, v in found.items() if v is None]
                found_by_run[entry.name] = found
                if missing:
                    failed.append(f"{entry.name} ({', '.join(missing)})")
                settled(entry, True)

            def error(message):
                self.project.message.emit(f"Anchor search failed for {entry.name}: {message}")
                settled(entry, False)

            self.project.jobs.submit(f"📍 Finding anchors: {entry.name}", work, done, error,
                                     on_cancel=lambda: settled(entry, False))

        for e in entries:
            start(e)

    # ------------------------------------------------------------- files
    def to_dict(self) -> dict:
        out = dict(reference=self.reference, kind=dict(self.kind),
                   points={run: [dict(name=p.name, ref_t1=p.ref_t1, ref_t2=p.ref_t2,
                                      run_t1=p.run_t1, run_t2=p.run_t2, ion=p.ion) for p in pts]
                           for run, pts in self.points.items()})
        # Only written when there is something to record, so projects without it are unchanged.
        if self.phase_at_anchor:
            out["phase_at_anchor"] = dict(self.phase_at_anchor)
        if self.phase_checked:
            out["phase_checked"] = sorted(self.phase_checked)
        return out

    def from_dict(self, d: dict) -> None:
        d = d or {}
        self.holding = False
        self.reference = d.get("reference")
        self.kind = dict(d.get("kind", {}))
        self.points = {run: [AnchorPoint(**p) for p in pts] for run, pts in d.get("points", {}).items()}
        self._recompute_all()
        # Projects saved before this record how each run was anchored have none: nothing is called stale.
        self.phase_at_anchor = {k: float(v) for k, v in d.get("phase_at_anchor", {}).items() if k in self.transforms}
        self.phase_checked = set(d.get("phase_checked", []))
        self.changed.emit()
        self.phaseStatusChanged.emit()
