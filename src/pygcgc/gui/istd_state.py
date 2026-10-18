"""🏷️ Internal-standard state for the app.

Holds the IS library, the peak labelled as each IS on a reference run, and
the per-run match results found by searching for them elsewhere. The
matching itself is all in `core/istd.py`; this just wires it to the
project (for runs, folds and peaks) and to background jobs, the same shape
as 📍's `AlignState`.

Peaks: like 🧬's `ClassState`, this wants full spectra. If the legacy
detector's blobs are there (`blobs_for`), they're used directly; otherwise,
for the MS-aware default, a run's components are re-solved on every ion via
`core.peaktable.respectrum_components`. Nothing is cached across calls --
each label or propagate re-resolves the one run it needs, inside its own
background job, the same cost model `ClassState.classify` already accepts.
"""

from __future__ import annotations

from dataclasses import asdict, replace
import hashlib
import json

import numpy as np
from PySide6 import QtCore

from ..core.istd import (NO_PEAKS, IstdLibraryEntry, IstdMatch, IstdMatchSettings, IstdReference,
                         load_library, match_run, peak_spectrum, save_library)
from ..core.peaktable import respectrum_components
from .blob_state import MS_AWARE


class IstdState(QtCore.QObject):
    changed = QtCore.Signal()   # library, a reference, or match results changed

    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.project = project
        self.library: list[IstdLibraryEntry] = load_library()
        self.reference_run: str | None = None
        self.references: dict[str, IstdReference] = {}
        self.matches: dict[str, dict[str, IstdMatch]] = {}
        self.settings = IstdMatchSettings()
        self.running = False
        self.match_inputs: dict[str, str] = {}
        self.review_decisions: dict[str, dict] = {}
        self._signature_cache = {}
        project.runsChanged.connect(self._prune)
        self.changed.connect(self._invalidate_reviewed_quant)
        for signal in (project.blobs.changed, project.layersChanged, project.foldChanged):
            signal.connect(self._invalidate_reviewed_quant)

    def _prune(self) -> None:
        names = {e.name for e in self.project.entries}
        if self.reference_run not in names:
            self.reference_run = None
        for name in [n for n in self.matches if n not in names]:
            self.matches.pop(name, None)

    def input_signature(self, index: int) -> str:
        """Stable result/input identity; uncacheable layers conservatively need re-review on reload."""
        p, b = self.project, self.project.blobs
        e = p.entries[index]
        result = b.components_for(index) if b.method == MS_AWARE else b.blobs_for(index)
        key = (p.base_key(index) or f"session:{id(p.base_run(index))}", p.period, e.phase,
               b.method, repr(b.component_settings if b.method == MS_AWARE else b.settings),
               tuple(id(g) for g in result.groups) if b.method == MS_AWARE and result is not None else id(result),
               repr(self.settings), tuple((r.name, r.t1, r.t2, hashlib.sha256(np.asarray(r.spectrum, dtype=np.float64).tobytes()).hexdigest())
                                         for r in self.references.values()))
        cached = self._signature_cache.get(e.name)
        if cached is not None and cached[0] == key:
            return cached[1]
        if result is None:
            payload = None
        elif b.method == MS_AWARE:
            payload = result.to_dict()
        else:
            payload = (result.t1.tolist(), result.t2.tolist(), result.volume.tolist(),
                       hashlib.sha256(result.spectra.tobytes()).hexdigest())
        # The object identities above only cache this computation; persisted signatures
        # contain the reproducible source/settings/result values instead.
        stable = (key[:5], key[6:], payload)
        signature = hashlib.sha256(json.dumps(stable, sort_keys=True, default=str).encode()).hexdigest()
        self._signature_cache[e.name] = (key, signature)
        return signature

    def is_stale(self, run: str) -> bool:
        index = self.project.index_of(run)
        if index is None or self.project.entries[index].run is None:
            return True
        return (self.project.blobs.running or self.project.blobs.selected_stale(self.project.entries[index])
                or self.match_inputs.get(run) != self.input_signature(index))

    def _review_key(self, run, name):
        return json.dumps([run, name])

    def reviewed(self, run, name) -> bool:
        decision = self.review_decisions.get(self._review_key(run, name))
        match = self.matches.get(run, {}).get(name)
        return bool(decision and match and not self.is_stale(run)
                    and decision['signature'] == self.match_inputs.get(run)
                    and json.dumps(decision['match'], sort_keys=True) == json.dumps(asdict(match), sort_keys=True))

    def reviewed_matches(self, run):
        return {name: match for name, match in self.matches.get(run, {}).items() if self.reviewed(run, name)}

    def exceptions(self):
        """Missing, uncertain, unsearched or stale standard/run pairs, in sequence order."""
        result = []
        for i in self.project.loaded():
            run = self.project.entries[i].name
            for name in self.references:
                match = self.matches.get(run, {}).get(name)
                if (match is None or self.is_stale(run) or not match.found or match.needs_review
                        or self._review_key(run, name) in self.review_decisions) and not self.reviewed(run, name):
                    result.append((run, name))
        return result

    def outdated_reviews(self):
        return [(self.project.entries[i].name, name) for i in self.project.loaded() for name in self.references
                if self._review_key(self.project.entries[i].name, name) in self.review_decisions
                and not self.reviewed(self.project.entries[i].name, name)]

    def _invalidate_reviewed_quant(self, *_):
        q = self.project.quant
        if self.review_decisions and (q.runquants or q.running) and not q.stale and self.outdated_reviews():
            q.mark_stale()

    def review_match(self, run, name, match=None, *, missing=False):
        """Save a per-run choice/acknowledgement without changing the global reference."""
        if self.project.quant.running:
            raise ValueError("Wait for quantification to finish before changing standard decisions.")
        if self.running:
            raise ValueError("Wait for the standard search to finish before changing review decisions.")
        if self.is_stale(run):
            raise ValueError("Standard matches are out of date; search this sequence again first.")
        current = self.matches.get(run, {}).get(name)
        if current is None:
            raise ValueError("Search for internal standards before recording a review.")
        if missing:
            match = replace(current, status="not found", peak_id=None,
                            reason="Reviewed: internal standard is absent; quantitative results must retain missing-IS flags")
        elif match is not None:
            match = replace(match, status="found", reason=("Manual per-run selection. " + match.reason).strip())
        else:
            match = current
        self.matches.setdefault(run, {})[name] = match
        self.review_decisions[self._review_key(run, name)] = dict(
            signature=self.match_inputs[run], match=asdict(match), missing=missing,
            manual=match is not current)
        self.project.quant.mark_stale()
        self.changed.emit()

    # ------------------------------------------------------------- library
    def add_custom(self, name: str, cas: str = "", product_code: str = "") -> None:
        self.library = self.library + [IstdLibraryEntry(name=name, cas=cas, product_code=product_code, custom=True)]
        save_library(self.library)
        self.changed.emit()

    def remove_from_library(self, name: str) -> None:
        self.library = [e for e in self.library if e.name != name]
        save_library(self.library)
        self.changed.emit()

    # ------------------------------------------------------------- peaks
    def _peaks_for(self, index: int):
        """This run's peaks with full spectra, however they were found --
        `None` if peak detection hasn't been run on it at all. Only ever
        called from inside a background job's `work()`."""
        p = self.project
        if not 0 <= index < len(p.entries):
            return None
        # The selected method decides, as it does for the click that picks the peak: legacy
        # blobs also exist on an MS-aware project once 🫧 or ⛓️ has found them, and a
        # component id cannot be looked up among them.
        if p.blobs.method == MS_AWARE:
            component_set = p.blobs.components_for(index)
            if component_set is None:
                return None
            base_run, fold = p.base_run(index), p.fold(index)
            peaks, _model, _footprint = respectrum_components(base_run, fold, component_set)
            return peaks
        return p.blobs.blobs_for(index)

    @staticmethod
    def _resolve_index(peaks, peak_ref) -> int | None:
        """`peak_ref` is a stable component id (MS-aware) or a positional
        blob index (legacy, the same positional trade-off `BlobState`
        itself already accepts for its edit log)."""
        if peaks is None or len(peaks) == 0:
            return None
        if isinstance(peak_ref, str):
            ids = getattr(peaks, "component_id", None)
            if not ids or peak_ref not in ids:
                return None
            return ids.index(peak_ref)
        i = int(peak_ref)
        return i if 0 <= i < len(peaks) else None

    # ------------------------------------------------------------- labelling
    def label(self, name: str, index: int, peak_ref, on_done=None) -> None:
        p = self.project
        if not 0 <= index < len(p.entries):
            return
        entry = p.entries[index]
        cas = next((e.cas for e in self.library if e.name == name), "")

        def work(token, progress):
            return self._peaks_for(index)

        def done(peaks) -> None:
            i = self._resolve_index(peaks, peak_ref)
            if i is None:
                p.message.emit(f"🏷️ Couldn't read that peak on {entry.name} -- run peak detection first")
                return
            self.references[name] = IstdReference(
                name=name, cas=cas, run=entry.name, t1=float(peaks.t1[i]), t2=float(peaks.t2[i]),
                mz0=int(peaks.mz[0]), spectrum=peak_spectrum(peaks, i))
            if self.reference_run is None:
                self.reference_run = entry.name
            self.matches.pop(entry.name, None)   # a fresh label invalidates any match already found on its own run
            self.changed.emit()
            if on_done:
                on_done()

        p.jobs.submit(f"🏷️ Labelling {name} on {entry.name}", work, done)

    def remove_reference(self, name: str) -> None:
        self.references.pop(name, None)
        for per_run in self.matches.values():
            per_run.pop(name, None)
        self.changed.emit()

    # ------------------------------------------------------------- batch
    def set_settings(self, settings: IstdMatchSettings) -> None:
        """New thresholds. Results already found stay until the search is run again."""
        self.settings = settings
        self.changed.emit()

    def propagate_to(self, indices: list[int], on_done=None) -> None:
        """Search every run in `indices` (the labelling run included) for every
        labelled standard, with the current thresholds."""
        if self.running:
            return
        if not self.references:
            self.project.message.emit("🏷️ Label at least one internal standard first")
            return
        p = self.project
        entries = [p.entries[i] for i in indices if p.entries[i].run is not None]
        if not entries:
            return
        remaining = set(entries)
        missing: list[str] = []
        references = list(self.references.values())
        settings = self.settings
        self.running = True
        unsuccessful = []
        self.changed.emit()

        def settled(entry, success=True):
            if not success:
                unsuccessful.append(entry.name)
            remaining.discard(entry)
            if not remaining:
                self.running = False
                self.changed.emit()
                if missing:
                    self.project.message.emit("Internal standards need review in: " + "; ".join(missing))
                if not unsuccessful and on_done:
                    on_done()

        def start(entry) -> None:
            idx = p.entries.index(entry)
            signature = self.input_signature(idx)

            def work(token, progress):
                peaks = self._peaks_for(idx)
                if peaks is None:
                    return [IstdMatch(run=entry.name, name=r.name, status=NO_PEAKS,
                                      reason="peak detection has not been run on this chromatogram")
                            for r in references]
                return match_run(references, entry.name, peaks, settings)

            def done(results) -> None:
                if entry not in p.entries or signature != self.input_signature(p.entries.index(entry)):
                    settled(entry, False)
                    return
                for j, match in enumerate(results):
                    decision = self.review_decisions.get(self._review_key(entry.name, match.name))
                    if decision and decision['signature'] == signature and decision.get('manual'):
                        results[j] = IstdMatch.from_dict(decision['match'])
                self.matches[entry.name] = {m.name: m for m in results}
                self.match_inputs[entry.name] = signature
                lost = [m.name for m in results if not m.found]
                if lost:
                    missing.append(f"{entry.name} ({', '.join(lost)})")
                settled(entry)

            p.jobs.submit(f"🏷️ Finding internal standards: {entry.name}", work, done,
                          lambda message, e=entry: (p.message.emit(message), settled(e, False)),
                          on_cancel=lambda e=entry: settled(e, False))

        for e in entries:
            start(e)

    # ------------------------------------------------------------- queries
    def all_matches(self) -> list[IstdMatch]:
        return [m for per_run in self.matches.values() for m in per_run.values()]

    def matches_for(self, index: int) -> dict[str, IstdMatch]:
        p = self.project
        return self.matches.get(p.entries[index].name, {}) if 0 <= index < len(p.entries) else {}

    # ------------------------------------------------------------- files
    def to_dict(self) -> dict:
        return dict(
            library=[e.to_dict() for e in self.library],
            reference_run=self.reference_run,
            references={name: dict(name=r.name, cas=r.cas, run=r.run, t1=r.t1, t2=r.t2, mz0=r.mz0,
                                   spectrum=r.spectrum.tolist())
                       for name, r in self.references.items()},
            matches={run: {name: asdict(m) for name, m in per_run.items()}
                     for run, per_run in self.matches.items()},
            settings=self.settings.to_dict(), match_inputs=self.match_inputs,
            review_decisions=self.review_decisions)

    def from_dict(self, d: dict) -> None:
        d = d or {}
        if d.get("library"):
            self.library = [IstdLibraryEntry.from_dict(e) for e in d["library"]]
        self.reference_run = d.get("reference_run")
        self.references = {
            name: IstdReference(name=r["name"], cas=r.get("cas", ""), run=r["run"], t1=r["t1"], t2=r["t2"],
                               mz0=r["mz0"], spectrum=np.asarray(r["spectrum"], dtype=float))
            for name, r in d.get("references", {}).items()}
        # Results saved before the three-threshold search (no "status") are dropped; Go re-finds them.
        self.matches = {
            run: {name: IstdMatch.from_dict(m) for name, m in per_run.items() if "status" in m}
            for run, per_run in d.get("matches", {}).items()}
        self.settings = IstdMatchSettings.from_dict(d.get("settings", {}))
        self.match_inputs = d.get("match_inputs", {})
        self.review_decisions = d.get("review_decisions", {})
        self._signature_cache = {}
        self.running = False
        self.changed.emit()
