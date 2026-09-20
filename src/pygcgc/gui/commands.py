"""Undoable edits.

Every change to how data are processed goes through a command on the undo
stack, so Ctrl+Z always works. View settings (zoom, colour levels, which
ion is shown) are not processing and are not on the stack.
"""

from __future__ import annotations

from PySide6.QtGui import QUndoCommand

from .project import Project

_PHASE_ID = 1001
_PERIOD_ID = 1002


class SetPhase(QUndoCommand):
    """Phase change for one run. Consecutive edits to the same run (slider
    drags, spin-box clicks) merge into one undo step."""

    def __init__(self, project: Project, index: int, new: float):
        e = project.entries[index]
        super().__init__(f"Phase of {e.name}: {new:.3f} s")
        self.project, self.index = project, index
        self.old = (e.phase, e.phase_method, e.phase_detail)
        self.new = (float(new), "manual", "set by hand")

    def id(self) -> int:
        return _PHASE_ID

    def mergeWith(self, other) -> bool:
        if not isinstance(other, SetPhase) or other.index != self.index:
            return False
        self.new = other.new
        self.setText(other.text())
        return True

    def redo(self) -> None:
        self.project.set_phase(self.index, *self.new)

    def undo(self) -> None:
        self.project.set_phase(self.index, *self.old)


class SetPhases(QUndoCommand):
    """Phase change for many runs at once (automatic estimate, copy to all)."""

    def __init__(self, project: Project, before: dict, after: dict, text: str):
        super().__init__(text)
        self.project = project
        self.before, self.after = before, after
        self._first = True

    def redo(self) -> None:
        # The change has usually been applied already when the command is
        # pushed (the estimate arrives from a background job), so skip once.
        if self._first:
            self._first = False
            return
        self.project.set_phases(self.after)

    def undo(self) -> None:
        self.project.set_phases(self.before)


class SetPeriod(QUndoCommand):
    def __init__(self, project: Project, new: float):
        super().__init__(f"Modulation period: {new:.4f} s")
        self.project = project
        self.old = project.fold_snapshot()
        self.new = float(new)

    def id(self) -> int:
        return _PERIOD_ID

    def mergeWith(self, other) -> bool:
        if not isinstance(other, SetPeriod):
            return False
        self.new = other.new
        self.setText(other.text())
        return True

    def redo(self) -> None:
        self.project.set_period(self.new)

    def undo(self) -> None:
        self.project.restore_fold(self.old)


class SetFold(QUndoCommand):
    """Period and every phase, estimated in the background and already
    applied when this is pushed."""

    def __init__(self, project: Project, before: dict, after: dict, text: str):
        super().__init__(text)
        self.project = project
        self.before, self.after = before, after
        self._first = True

    def redo(self) -> None:
        if self._first:
            self._first = False
            return
        self.project.restore_fold(self.after)

    def undo(self) -> None:
        self.project.restore_fold(self.before)


class ApplyBackground(QUndoCommand):
    """🧹 layers added, changed or removed on some runs. Results are computed
    in the background and already applied when this is pushed."""

    def __init__(self, project: Project, before: dict, after: dict,
                 before_settings, after_settings, before_blank, after_blank, text: str):
        super().__init__(text)
        self.project = project
        self.before, self.after = before, after
        self.bs, self.as_ = before_settings, after_settings
        self.bb, self.ab = before_blank, after_blank
        self._first = True

    def redo(self) -> None:
        if self._first:
            self._first = False
            return
        self.project.set_backgrounds(self.after, self.as_, self.ab)

    def undo(self) -> None:
        self.project.set_backgrounds(self.before, self.bs, self.bb)


class ApplySiloxane(QUndoCommand):
    """🫧 layers added or removed on some runs. Results are computed in the
    background and already applied when this is pushed."""

    def __init__(self, project: Project, before: dict, after: dict, before_settings, after_settings, text: str):
        super().__init__(text)
        self.project = project
        self.before, self.after = before, after
        self.bs, self.as_ = before_settings, after_settings
        self._first = True

    def redo(self) -> None:
        if self._first:
            self._first = False
            return
        self.project.set_siloxanes(self.after, self.as_)

    def undo(self) -> None:
        self.project.set_siloxanes(self.before, self.bs)


class EditDecisions(QUndoCommand):
    """A review decision: accept, reject, exclude a calibration point,
    override a curve's settings, or pick a peak by hand."""

    def __init__(self, quant, before, after, text: str, requantify=()):
        super().__init__(text)
        self.quant = quant
        self.before, self.after = before, after
        self.requantify = tuple(requantify)

    def redo(self) -> None:
        self.quant.set_decisions(self.after.copy(), self.requantify)

    def undo(self) -> None:
        self.quant.set_decisions(self.before.copy(), self.requantify)


class EditAlignment(QUndoCommand):
    """An anchor added or removed, or a run's transform kind changed.
    Applied already when pushed (the point was recorded, or the kind
    picked, before this command exists); redo is a no-op the first time."""

    def __init__(self, align, before: dict, after: dict, text: str):
        super().__init__(text)
        self.align, self.before, self.after = align, before, after
        self._first = True

    def redo(self) -> None:
        if self._first:
            self._first = False
            return
        self.align.restore(self.after)

    def undo(self) -> None:
        self.align.restore(self.before)


class EditBlobs(QUndoCommand):
    """A blob deleted, merged or added by hand on one run. Applied already
    when pushed; redo is a no-op the first time."""

    def __init__(self, blobs, before, after, text: str):
        super().__init__(text)
        self.blobs, self.before, self.after = blobs, before, after
        self._first = True

    def redo(self) -> None:
        if self._first:
            self._first = False
            return
        self.blobs.restore_run(self.after)

    def undo(self) -> None:
        self.blobs.restore_run(self.before)


class EditComponents(QUndoCommand):
    """An accepted MS-aware component-group fit, applied before pushing."""

    def __init__(self, peaks, before, after, text: str):
        super().__init__(text)
        self.peaks, self.before, self.after = peaks, before, after
        self._first = True

    def redo(self) -> None:
        if self._first:
            self._first = False
            return
        self.peaks.restore_components(self.after)

    def undo(self) -> None:
        self.peaks.restore_components(self.before)


class EditHydrocarbons(QUndoCommand):
    """A stencil added, removed or reweighted. Applied already when
    pushed; redo is a no-op the first time."""

    def __init__(self, hydro, before, after, text: str):
        super().__init__(text)
        self.hydro, self.before, self.after = hydro, before, after
        self._first = True

    def redo(self) -> None:
        if self._first:
            self._first = False
            return
        self.hydro.restore_stencils(self.after)

    def undo(self) -> None:
        self.hydro.restore_stencils(self.before)


class SetReference(QUndoCommand):
    """Choosing a reference run clears every anchor (they were relative to
    the old one), so this is its own command rather than an EditAlignment."""

    def __init__(self, align, before: dict, after_name: str | None, text: str):
        super().__init__(text)
        self.align, self.before, self.after_name = align, before, after_name
        self._first = True

    def redo(self) -> None:
        if self._first:
            self._first = False
            return
        self.align.set_reference(self.after_name)

    def undo(self) -> None:
        self.align.restore(self.before)
