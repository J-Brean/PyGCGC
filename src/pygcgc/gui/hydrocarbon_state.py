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
                                 find_alkane_series, find_markers, propagate_alkane_series, snap_ladder,
                                 stencil_from_polygon)
from ..core.peaktable import component_positions
from .blob_state import MS_AWARE


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

    # ------------------------------------------------------------ the peaks
    # Both peak finders work here. Legacy blobs carry everything. MS-aware components give positions and
    # volumes at once; their all-ion spectra (needed to say what kind of hydrocarbon a peak is) take a job,
    # which `prepare` runs and `spectra_for` then reads from the cache.
    def positions_for(self, index: int):
        """The run's peaks with positions and volumes: its legacy blobs, or its fitted components. None if neither."""
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        if p.blobs.method == MS_AWARE:
            components = p.blobs.components_for(index)
            return component_positions(components) if components is not None else None
        return p.blobs.blobs_for(index)

    def spectra_for(self, index: int):
        """The run's peaks with a mass spectrum on every ion, if they are ready (never starts any work)."""
        return self.project.blobs.cached_full_spectrum_peaks(index)

    def peak_word(self) -> str:
        return "MS-aware peaks" if self.project.blobs.method == MS_AWARE else "legacy blobs"

    # ------------------------------------------------------------ assigning
    def calls_for(self, index: int) -> list[HydrocarbonCall] | None:
        """The run's assignments, or None if there are none or the peaks have changed since."""
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        calls = self.calls.get(p.entries[index].name)
        peaks = self.positions_for(index)
        return calls if calls is not None and peaks is not None and len(calls) == len(peaks) else None

    def can_assign(self, index: int) -> str:
        """Empty if run `index` can be assigned now, else what is missing, in a sentence."""
        series = self.series_for(index)
        if series is None or len(series) < 2:
            return "no alkane ladder (two alkanes or more)"
        peaks = self.positions_for(index)
        if peaks is None or len(peaks) == 0:
            return f"no peaks (find the {self.peak_word()} first)"
        return ""

    def assign(self, indices) -> tuple[int, list[str]]:
        """Assign the peaks of each run that has a ladder (two alkanes or more), peaks and, for MS-aware peaks, their
        all-ion spectra already prepared. Returns the number of runs done and the names of the runs that lacked
        something; `assign_async` prepares the spectra first."""
        p = self.project
        done, skipped = 0, []
        for i in indices:
            name = p.entries[i].name
            peaks = self.spectra_for(i)
            if self.can_assign(i) or peaks is None or len(peaks) == 0:
                skipped.append(name)
                continue
            self.calls[name] = assign_hydrocarbons(peaks, self.series[name], self.settings)
            done += 1
        if done:
            self.changed.emit()
        return done, skipped

    def assign_async(self, indices, on_done=None) -> int:
        """`assign` for runs whose spectra may not be ready: prepare each in a job, then assign. `on_done(done, skipped)`
        runs at the end. Returns the number of runs being worked on."""
        p = self.project
        wanted = [i for i in indices if not self.can_assign(i)]
        skipped = [p.entries[i].name for i in indices if self.can_assign(i)]
        pending = set(wanted)
        failed: list[str] = []

        def settle(i, ok):
            pending.discard(i)
            if not ok:
                failed.append(p.entries[i].name)
            if pending:
                return
            done, short = self.assign([j for j in wanted if p.entries[j].name not in failed])
            if on_done:
                on_done(done, skipped + failed + short)

        if not wanted and on_done:
            on_done(0, skipped)
        for i in wanted:
            if not p.blobs.prepare_full_spectrum_peaks(i, on_done=lambda i=i: settle(i, True),
                                                       on_fail=lambda _m, i=i: settle(i, False)):
                settle(i, False)
        return len(wanted)

    def refresh_calls(self, index: int) -> None:
        """Re-assign a run whose assignments exist and whose spectra are ready (the ladder moved under them)."""
        name = self.project.entries[index].name
        peaks = self.spectra_for(index)
        if name in self.calls and peaks is not None and len(peaks) and not self.can_assign(index):
            self.calls[name] = assign_hydrocarbons(peaks, self.series[name], self.settings)
        else:
            self.calls.pop(name, None)

    def series_for(self, index: int):
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        return self.series.get(p.entries[index].name)

    # ------------------------------------------------------------ a ladder drawn or snapped by hand
    def snapshot_series(self, index: int):
        """A run's ladder as it is now (a copy, or None), for undo."""
        series = self.series_for(index)
        return None if series is None else series.copy()

    def set_series(self, index: int, series: AlkaneSeries | None) -> None:
        """Replace (or, with None, remove) a run's ladder. Everything that reads the ladder follows: the assignments
        are redone if their spectra are ready, the stencils are lifted onto the new line and the library matcher's
        results are marked out of date, because the retention indices moved."""
        p = self.project
        name = p.entries[index].name
        if series is None:
            self.series.pop(name, None)
            self.calls.pop(name, None)
        else:
            series.shape, series.degree = self.settings.ladder_shape, self.settings.ladder_degree
            self.series[name] = series
            self.refresh_calls(index)
        p.request_autosave()
        self.changed.emit()

    def snap(self, index: int, on_done=None, on_fail=None) -> bool:
        """Move a drawn ladder onto the real alkane apexes, in a job. `on_done(series, snapped_carbons)` receives
        the result, which is not stored: the caller decides (so it can be undone as one step). `on_fail(message)`
        hears of a search that could not run, for example a run without the alkane ion."""
        if self.running:
            return False
        p = self.project
        series = self.series_for(index)
        if series is None or len(series) < 2:
            p.message.emit("Draw a ladder of at least two carbons first")
            return False
        run, fold = p.display_run(index), p.fold(index)
        settings, name = self.settings, p.entries[index].name
        inputs = p.analysis_inputs()
        self._set_running(True)

        def work(token, progress):
            progress(0, 0, "Snapping the ladder to the alkane apexes")
            return snap_ladder(run, fold, series, settings, cancel=token.check)

        def done(result):
            self._set_running(False)
            if inputs != p.analysis_inputs() or settings != self.settings or self.series_for(index) is not series:
                p.message.emit("Snapping discarded: the run, settings or ladder changed while it ran")
                return
            if on_done:
                on_done(*result)

        def failed(message):
            self._set_running(False)
            p.message.emit(f"⛓️ {name}: {message}")
            if on_fail:
                on_fail(message)

        p.jobs.submit(f"⛓️ Snapping the ladder: {name}", work, done, failed, on_cancel=lambda: self._set_running(False))
        return True

    def _set_running(self, on: bool) -> None:
        self.running = on
        self.runningChanged.emit(on)

    # ------------------------------------------------------------ finding
    def find_series(self, index: int, seeds, on_done=None) -> None:
        if self.running:
            return
        p = self.project
        e = p.entries[index]
        run, fold = p.display_run(index), p.fold(index)
        settings = self.settings
        inputs = p.analysis_inputs()
        self._set_running(True)

        def work(token, progress):
            return find_alkane_series(run, fold, seeds, settings)

        def done(series):
            if inputs != p.analysis_inputs() or settings != self.settings:
                self._set_running(False)
                p.message.emit("Alkane search discarded: inputs changed; search again")
                return
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

        p.jobs.submit(f"⛓️ Finding alkanes: {e.name}", work, done, failed,
                      on_cancel=lambda: self._set_running(False))

    def propagate_series(self, from_index: int, to_indices, on_done=None) -> None:
        if self.running:
            return
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
        inputs = p.analysis_inputs()
        missing = []
        unsuccessful = set()
        self._set_running(True)

        def settled(entry, success):
            if not success:
                unsuccessful.add(entry)
            remaining.discard(entry)
            if not remaining:
                self._set_running(False)
                self.changed.emit()
                p.request_autosave()
                if missing:
                    p.message.emit("No alkanes found in: " + ", ".join(missing))
                if not unsuccessful and on_done:
                    on_done()
        for e in entries:
            idx = p.entries.index(e)
            run, fold = p.display_run(idx), p.fold(idx)

            def work(token, progress, run=run, fold=fold):
                return propagate_alkane_series(ref_series, run, fold, settings)

            def done(series, e=e):
                if (inputs != p.analysis_inputs() or settings != self.settings
                        or self.series_for(from_index) is not ref_series):
                    p.message.emit(f"Alkane propagation for {e.name} discarded: inputs changed")
                    settled(e, False)
                    return
                if len(series):
                    series.shape, series.degree = settings.ladder_shape, settings.ladder_degree
                    self.series[e.name] = series
                else:
                    missing.append(e.name)
                settled(e, True)

            def failed(message, e=e):
                p.message.emit(f"Alkane propagation for {e.name} failed: {message}")
                settled(e, False)

            p.jobs.submit(f"⛓️ Propagating alkanes: {e.name}", work, done, failed,
                          on_cancel=lambda e=e: settled(e, False))

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
        peaks = self.positions_for(index)
        if series is None or len(series) < 2 or peaks is None or not self.stencils.stencils:
            return {}
        return self.stencils.totals(series, peaks)

    def markers(self, index: int) -> dict:
        series = self.series_for(index)
        peaks = self.positions_for(index)
        if series is None or len(series) < 2 or peaks is None:
            return {}
        return find_markers(peaks, series)

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
