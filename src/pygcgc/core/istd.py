"""🏷️ Internal standards: a user-labelled library, and cross-run matching.

The idea in plain terms
------------------------
Today an internal standard (IS) is only ever a hand-typed target-table row
(`core/targets.py`, `role="istd"`), located on each run by a windowed apex
search on one declared ion (`core/phase.py:locate_anchor`) -- no spectral
comparison, no confidence score. This module is the alternative: a peak
*clicked* on one reference chromatogram (its own position and spectrum,
captured on the spot) is searched for on every other chromatogram using
both position and mass-spectral similarity, with a confidence score and
tier reported per run so a weak match can be spotted and corrected by hand.

1. **Library** (`IstdLibraryEntry`, `load_library`/`save_library`). A small
   list of name/CAS/product-code rows -- the product code is a lab
   catalogue number, not a retention time, so it carries no position
   information. Seeded from the shipped `presets/istd_library.json` into a
   user-writable file the first time it's needed, so custom entries a user
   adds persist across projects and sessions.
2. **Reference** (`IstdReference`). One labelled peak: which IS, which run
   it was clicked on, its t1/t2, and its own spectrum -- captured exactly
   the way `identify.py` captures an accepted standard's spectrum
   (`identify.build_inhouse_entry`'s docstring: `BlobSet.apex_spectra[i]`,
   `mz0 = int(blobs.mz[0])`).
3. **Matching** (`find_istd_on_run`). Candidates are every peak on the
   target run within a generous position window of the reference (position
   alone is not trusted -- that's the whole point of adding a spectral
   term: it has to survive the retention drift a position-only search
   can't). Each candidate is scored two ways and combined: a position
   score (a scaled-distance falloff, the same idea as `BlobSet.nearest`)
   and a spectral score (`identify.search`'s Stein and Scott (1994)
   dot-product scoring, candidate as the query, the reference as a
   one-entry library -- exactly the role `identify.py` normally has query
   and library play, so its asymmetric reverse score means what it always
   means: ignore channels the *reference* doesn't have, i.e. don't
   penalise a candidate for a co-eluting neighbour's extra ions).
4. **Confidence** (`IstdMatchSettings` tiers). A simple threshold on the
   combined score -- "high"/"medium"/"low", or "not found" when nothing
   at all falls inside the position window.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .align import AlignmentError, AnchorPoint, fit_alignment
from .identify import IdentifySettings, build_inhouse_entry, build_library, search

PRESET_LIBRARY_PATH = Path(__file__).resolve().parents[1] / "presets" / "istd_library.json"

# The candidate/reference search itself should never drop or re-threshold a
# hit -- that's `find_istd_on_run`'s own job, via the combined score and its
# tiers -- so every prefilter `identify.search` would otherwise apply is
# switched off here.
_SEARCH_SETTINGS = IdentifySettings(min_score=0.0, min_shared_top_ions=0)


# ------------------------------------------------------------------ library
@dataclass
class IstdLibraryEntry:
    """One internal standard a user can label peaks against. `product_code`
    is a lab catalogue number (not a retention time) kept only for display."""
    name: str
    cas: str = ""
    product_code: str = ""
    custom: bool = False

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "IstdLibraryEntry":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


def default_user_library_path() -> Path:
    """Where a user's IS library lives, seeded from the shipped defaults on
    first use. No user-config-directory convention exists elsewhere in this
    codebase (no appdirs/platformdirs/QStandardPaths) -- this is a new one,
    a plain dotfolder under the user's home directory."""
    return Path.home() / ".pygcgc" / "istd_library.json"


def load_library(path: Path | None = None) -> list[IstdLibraryEntry]:
    path = Path(path) if path is not None else default_user_library_path()
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(PRESET_LIBRARY_PATH, path)
    with open(path, "r", encoding="utf-8") as f:
        rows = json.load(f)
    return [IstdLibraryEntry.from_dict(r) for r in rows]


def save_library(entries: list[IstdLibraryEntry], path: Path | None = None) -> None:
    path = Path(path) if path is not None else default_user_library_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump([e.to_dict() for e in entries], f, indent=2)


def add_custom_entry(entries: list[IstdLibraryEntry], name: str, cas: str = "",
                     product_code: str = "", path: Path | None = None) -> list[IstdLibraryEntry]:
    """`entries` plus one new custom entry, saved immediately so it survives
    across projects and sessions."""
    updated = list(entries) + [IstdLibraryEntry(name=name, cas=cas, product_code=product_code, custom=True)]
    save_library(updated, path)
    return updated


# ------------------------------------------------------------------ reference & matches
@dataclass
class IstdReference:
    """One IS as labelled on the reference run: its own position and
    spectrum, captured from the clicked peak."""
    name: str
    cas: str
    run: str
    t1: float
    t2: float
    mz0: int
    spectrum: np.ndarray


@dataclass
class IstdMatchSettings:
    t1_window: float = 90.0     # s either side of the reference, for the candidate search box
    t2_window: float = 1.5      # s either side
    t1_scale: float = 60.0      # scaled-distance denominator, same idea as BlobSet.nearest's default
    t2_scale: float = 0.5
    position_weight: float = 0.4
    spectral_weight: float = 0.6
    tier_high: float = 0.75
    tier_medium: float = 0.45
    # First pass (`match_run`): a wider box and a spectrum-only criterion for the
    # confident matches whose drift predicts where every other IS should be.
    anchor_t1_window: float = 180.0
    anchor_t2_window: float = 2.5
    anchor_min_spectral: float = 0.75

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "IstdMatchSettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


@dataclass
class IstdMatch:
    run: str
    name: str
    peak_id: object                     # component id (MS-aware) or blob index (legacy), or None if not found
    t1: float
    t2: float
    dt1: float
    dt2: float
    position_score: float
    spectral_score: float
    forward: float
    reverse: float
    combined_score: float
    tier: str                           # "high", "medium", "low", "not found"
    expected_t1: float = float("nan")   # where the other ISs' drift predicts it; the reference position if none
    expected_t2: float = float("nan")
    n_drift_anchors: int = 0            # how many other ISs that prediction came from


def position_score(dt1: float, dt2: float, settings: IstdMatchSettings = IstdMatchSettings()) -> float:
    return float(np.exp(-np.hypot(dt1 / settings.t1_scale, dt2 / settings.t2_scale)))


def _tier(score: float, settings: IstdMatchSettings) -> str:
    if score >= settings.tier_high:
        return "high"
    if score >= settings.tier_medium:
        return "medium"
    return "low"


def peak_spectrum(peaks, i: int) -> np.ndarray:
    """A peak's own spectrum: the apex (single-scan) spectrum where one was
    computed, the footprint-summed spectrum otherwise -- the same fallback
    `identify.py` uses (`identify.py:424`, `use_apex`). Public: also used to
    capture a reference peak's spectrum from a GUI click (`gui/istd_state.py`)."""
    use_apex = peaks.apex_spectra.shape[0] == len(peaks) and peaks.apex_spectra.size
    return peaks.apex_spectra[i] if use_apex else peaks.spectra[i]


def _peak_id(peaks, i: int):
    component_id = getattr(peaks, "component_id", None)
    return component_id[i] if component_id else int(i)


def _spectral_scores(reference: IstdReference, peaks, indices) -> list[tuple[float, float, float]]:
    """(score, forward, reverse) of each candidate against the reference:
    candidate as the query, reference as a one-entry library."""
    ref_library = build_library([build_inhouse_entry(reference.name, reference.spectrum, reference.mz0,
                                                     cas=reference.cas)])
    mz0 = int(peaks.mz[0])
    out = []
    for i in indices:
        cand = build_inhouse_entry(f"candidate_{i}", peak_spectrum(peaks, int(i)), mz0)
        hits = search(cand.mz, cand.intensity, ref_library, query_ri=None, settings=_SEARCH_SETTINGS)
        out.append((hits[0].score, hits[0].forward, hits[0].reverse) if hits else (0.0, 0.0, 0.0))
    return out


def find_istd_on_run(reference: IstdReference, run_name: str, peaks,
                     settings: IstdMatchSettings = IstdMatchSettings(),
                     expected: tuple[float, float] | None = None, n_drift_anchors: int = 0) -> IstdMatch:
    """The best-matching peak for `reference` on `peaks` (a `BlobSet` or a
    `ComponentPeaks` -- anything with the same `t1`/`t2`/`mz`/`spectra`
    shape), by position and spectral similarity combined. The search box
    and the position score are centred on `expected` (where this IS is
    predicted to sit on this run), or on the reference position if none is
    given. `dt1`/`dt2` are always reported against the reference. `tier="not
    found"` if nothing at all falls inside the box."""
    e1, e2 = expected if expected is not None else (reference.t1, reference.t2)
    not_found = IstdMatch(run=run_name, name=reference.name, peak_id=None,
                          t1=float("nan"), t2=float("nan"), dt1=float("nan"), dt2=float("nan"),
                          position_score=0.0, spectral_score=0.0, forward=0.0, reverse=0.0,
                          combined_score=0.0, tier="not found", expected_t1=float(e1), expected_t2=float(e2),
                          n_drift_anchors=n_drift_anchors)
    if len(peaks) == 0:
        return not_found
    t1 = np.asarray(peaks.t1, dtype=float)
    t2 = np.asarray(peaks.t2, dtype=float)
    r1, r2 = t1 - e1, t2 - e2
    in_window = np.flatnonzero((np.abs(r1) <= settings.t1_window) & (np.abs(r2) <= settings.t2_window))
    if in_window.size == 0:
        return not_found

    best = None
    for i, (spectral, forward, reverse) in zip(in_window, _spectral_scores(reference, peaks, in_window)):
        pos = position_score(float(r1[i]), float(r2[i]), settings)
        combined = settings.position_weight * pos + settings.spectral_weight * spectral
        if best is None or combined > best.combined_score:
            best = IstdMatch(run=run_name, name=reference.name, peak_id=_peak_id(peaks, int(i)),
                             t1=float(t1[i]), t2=float(t2[i]),
                             dt1=float(t1[i] - reference.t1), dt2=float(t2[i] - reference.t2),
                             position_score=float(pos), spectral_score=float(spectral),
                             forward=float(forward), reverse=float(reverse),
                             combined_score=float(combined), tier=_tier(combined, settings),
                             expected_t1=float(e1), expected_t2=float(e2), n_drift_anchors=n_drift_anchors)
    return best


def match_run(references: list[IstdReference], run_name: str, peaks,
              settings: IstdMatchSettings = IstdMatchSettings()) -> list[IstdMatch]:
    """Every reference on one run, allowing for the run's own retention drift.

    Internal standards drift together, so an absolute position penalty
    would mark every correct match on a shifted run as doubtful. First,
    each IS's best spectral match inside a wide box is taken as an anchor
    if its spectral score alone clears `anchor_min_spectral`. Then each IS
    is scored against where the *other* anchors predict it should be: their
    median t1/t2 offset (`core.align.fit_alignment(..., "shift")`, fitted
    leave-one-out, so an IS never vouches for its own position). With no
    other anchor the reference position is used."""
    anchors = []
    if len(peaks):
        t1 = np.asarray(peaks.t1, dtype=float)
        t2 = np.asarray(peaks.t2, dtype=float)
        for ref in references:
            box = np.flatnonzero((np.abs(t1 - ref.t1) <= settings.anchor_t1_window)
                                 & (np.abs(t2 - ref.t2) <= settings.anchor_t2_window))
            if box.size == 0:
                continue
            scores = [s for s, _f, _r in _spectral_scores(ref, peaks, box)]
            k = int(np.argmax(scores))
            if scores[k] >= settings.anchor_min_spectral:
                anchors.append(AnchorPoint(ref.name, ref.t1, ref.t2, float(t1[box[k]]), float(t2[box[k]])))

    out = []
    for ref in references:
        others = [a for a in anchors if a.name != ref.name]
        expected = None
        if others:
            try:
                # A median shift, not "auto": with a handful of anchors, one of them
                # possibly mislabelled, an affine or smooth fit extrapolates wildly.
                tf = fit_alignment(others, "shift")
                p1, p2 = tf.inverse(ref.t1, ref.t2)
                expected = (float(p1), float(p2))
            except AlignmentError:
                expected = None
        out.append(find_istd_on_run(ref, run_name, peaks, settings, expected=expected,
                                    n_drift_anchors=len(others) if expected is not None else 0))
    return out


# ------------------------------------------------------------------ report
def build_report(matches: list[IstdMatch]) -> list[dict]:
    return [dict(run=m.run, name=m.name, t1=m.t1, t2=m.t2, dt1=m.dt1, dt2=m.dt2,
                position_score=m.position_score, spectral_score=m.spectral_score,
                forward=m.forward, reverse=m.reverse, combined_score=m.combined_score, tier=m.tier,
                expected_t1=m.expected_t1, expected_t2=m.expected_t2, n_drift_anchors=m.n_drift_anchors)
            for m in matches]


def write_istd_report_csv(matches: list[IstdMatch], path) -> Path:
    import csv
    path = Path(path)
    rows = build_report(matches)
    fieldnames = ["run", "name", "t1", "t2", "dt1", "dt2", "position_score", "spectral_score",
                 "forward", "reverse", "combined_score", "tier", "expected_t1", "expected_t2", "n_drift_anchors"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path
