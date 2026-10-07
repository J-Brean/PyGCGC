"""⛓️ Hydrocarbon and isomer mapper state for the app.

Keeps one `AlkaneSeries` per run name (found by two clicks, or propagated
from another run's; saved with the project), and a single `StencilSet` --
the "frame" -- drawn once
against whichever run's series was current, and lifted onto any other
run's own series on demand. Finding/propagating the ladder isn't
undoable (like 🔍's `find_blobs`, it's a detection result, not an edit);
adding, removing or reweighting a stencil is (like 🔍's delete/merge/add).
The science is all in `core/hydrocarbons.py`.
"""

from __future__ import annotations

from PySide6 import QtCore

from ..core.hcassign import HydrocarbonCall, assign_hydrocarbons
from ..core.hydrocarbons import (AlkaneHit, AlkaneSeries, HydrocarbonSettings, StencilSet,
                                 find_alkane_series, find_markers, propagate_alkane_series,
                                 stencil_from_polygon)


class HydrocarbonState(QtCore.QObject):
    changed = QtCore.Signal()          # series, stencils or settings changed
    runningChanged = QtCore.Signal(bool)

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.settings = HydrocarbonSettings()
        self.series: dict[str, AlkaneSeries] = {}     # run name -> ladder
        self.stencils = StencilSet(name="frame")
        self.calls: dict[str, list[HydrocarbonCall]] = {}   # run name -> hydrocarbon of every blob (not saved: cheap to redo)
        self.running = False

    def forget_runs(self, entries) -> None:
        for e in entries:
            self.series.pop(e.name, None)
            self.calls.pop(e.name, None)

    def set_settings(self, settings: HydrocarbonSettings) -> None:
        """Store new settings and redraw every ladder with the chosen shape."""
        self.settings = settings
        self._apply_shape()
        self.changed.emit()

    def _apply_shape(self) -> None:
        for s in self.series.values():
            s.shape, s.degree = self.settings.ladder_shape, self.settings.ladder_degree

    # ------------------------------------------------------------ assigning
    def calls_for(self, index: int) -> list[HydrocarbonCall] | None:
        """The run's assignments, or None if there are none or the blobs have changed since."""
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        calls = self.calls.get(p.entries[index].name)
        blobs = p.blobs.blobs_for(index)
        return calls if calls is not None and blobs is not None and len(calls) == len(blobs) else None

    def assign(self, indices) -> tuple[int, list[str]]:
        """Assign the blobs of each run that has both a ladder (two alkanes or more) and blobs.
        Returns the number of runs done and the names of the runs that lacked something."""
        p = self.project
        done, skipped = 0, []
        for i in indices:
            name = p.entries[i].name
            series, blobs = self.series.get(name), p.blobs.blobs_for(i)
            if series is None or len(series) < 2 or blobs is None or len(blobs) == 0:
                skipped.append(name)
                continue
            self.calls[name] = assign_hydrocarbons(blobs, series, self.settings)
            done += 1
        if done:
            self.changed.emit()
        return done, skipped

    def series_for(self, index: int):
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        return self.series.get(p.entries[index].name)

    def _set_running(self, on: bool) -> None:
        self.running = on
        self.runningChanged.emit(on)

    # ------------------------------------------------------------ finding
    def find_series(self, index: int, seeds, on_done=None) -> None:
        p = self.project
        e = p.entries[index]
        run, fold = p.display_run(index), p.fold(index)
        settings = self.settings
        self._set_running(True)

        def work(token, progress):
            return find_alkane_series(run, fold, seeds, settings)

        def done(series):
            series.shape, series.degree = settings.ladder_shape, settings.ladder_degree
            self.series[e.name] = series
            self._set_running(False)
            p.request_autosave()
            self.changed.emit()
            if on_done:
                on_done(series)

        def failed(msg):
            self._set_running(False)
            p.message.emit(f"⛓️ {e.name}: {msg}")

        p.jobs.submit(f"⛓️ Finding alkanes: {e.name}", work, done, failed)

    def propagate_series(self, from_index: int, to_indices, on_done=None) -> None:
        p = self.project
        ref_series = self.series_for(from_index)
        if ref_series is None:
            p.message.emit("Find the alkane series on this run first")
            return
        entries = [p.entries[i] for i in to_indices
                  if p.entries[i].run is not None and i != from_index]
        if not entries:
            if on_done:
                on_done()
            return
        remaining = set(entries)
        settings = self.settings
        missing = []
        self._set_running(True)
        for e in entries:
            idx = p.entries.index(e)
            run, fold = p.display_run(idx), p.fold(idx)

            def work(token, progress, run=run, fold=fold):
                return propagate_alkane_series(ref_series, run, fold, settings)

            def done(series, e=e):
                if len(series):
                    series.shape, series.degree = settings.ladder_shape, settings.ladder_degree
                    self.series[e.name] = series
                else:
                    missing.append(e.name)
                remaining.discard(e)
                if not remaining:
                    self._set_running(False)
                    self.changed.emit()
                    p.request_autosave()
                    if missing:
                        p.message.emit("No alkanes found in: " + ", ".join(missing))
                    if on_done:
                        on_done()

            p.jobs.submit(f"⛓️ Propagating alkanes: {e.name}", work, done)

    # ------------------------------------------------------------ stencils
    def snapshot_stencils(self) -> StencilSet:
        return StencilSet(name=self.stencils.name, stencils=list(self.stencils.stencils))

    def restore_stencils(self, snap: StencilSet) -> None:
        self.stencils = StencilSet(name=snap.name, stencils=list(snap.stencils))
        self.changed.emit()

    def add_stencil(self, index: int, vertices_t1t2, name: str, response_factor: float = 1.0) -> bool:
        series = self.series_for(index)
        if series is None:
            self.project.message.emit("Find the alkane series on this run first")
            return False
        stencil = stencil_from_polygon(vertices_t1t2, series, name, response_factor)
        self.stencils.stencils = [s for s in self.stencils.stencils if s.name != name] + [stencil]
        self.changed.emit()
        return True

    def remove_stencils(self, names) -> None:
        names = set(names)
        self.stencils.stencils = [s for s in self.stencils.stencils if s.name not in names]
        self.changed.emit()

    def set_response_factor(self, name: str, value: float) -> None:
        for s in self.stencils.stencils:
            if s.name == name:
                s.response_factor = float(value)
        self.changed.emit()

    def totals(self, index: int) -> dict:
        series = self.series_for(index)
        blobs = self.project.blobs.blobs_for(index)
        if series is None or blobs is None or not self.stencils.stencils:
            return {}
        return self.stencils.totals(series, blobs)

    def markers(self, index: int) -> dict:
        series = self.series_for(index)
        blobs = self.project.blobs.blobs_for(index)
        if series is None or blobs is None:
            return {}
        return find_markers(blobs, series)

    # -------------------------------------------------------------- files
    def to_dict(self) -> dict:
        # Ladders are small and cost clicks to make, so they are kept with the project.
        series = {name: [dict(carbon=h.carbon, t1=h.t1, t2=h.t2, height=h.height, qual_ratios=list(h.qual_ratios),
                              seeded=h.seeded, qual_ok=h.qual_ok) for h in s.hits.values()]
                  for name, s in self.series.items()}
        return dict(settings=self.settings.to_dict(), stencils=self.stencils.to_dict(), series=series)

    def from_dict(self, d: dict) -> None:
        d = d or {}
        self.settings = HydrocarbonSettings.from_dict(d.get("settings", {}))
        self.stencils = StencilSet.from_dict(d.get("stencils") or {"name": "frame", "stencils": []})
        self.series = {name: AlkaneSeries(hits={int(h["carbon"]): AlkaneHit(
                           carbon=int(h["carbon"]), t1=float(h["t1"]), t2=float(h["t2"]), height=float(h["height"]),
                           qual_ratios=tuple(h.get("qual_ratios", ())), seeded=bool(h.get("seeded")),
                           qual_ok=bool(h.get("qual_ok", True))) for h in hits})
                       for name, hits in (d.get("series") or {}).items()}
        self.calls = {}
        self._apply_shape()
        self.changed.emit()
