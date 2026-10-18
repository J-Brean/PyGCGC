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
3. **Matching** (`find_istd_on_run`, `match_run`). Three thresholds, one
   per dimension, all of which a peak must pass to count as the standard:
   its first-dimension time within `max_dt1` seconds of where the standard
   is expected, its second-dimension time within `max_dt2`, and its mass
   spectrum matching the reference at `min_spectral` or better. Of the peaks
   that pass, the best spectral match is taken. The spectral match is
   `identify.search`'s Stein and Scott (1994) weighted dot-product score,
   candidate as the query and the reference as a one-entry library, so its
   reverse score ignores ions the reference does not have (a co-eluting
   neighbour's extra ions are not held against the candidate).
4. **Drift** (`match_run`). Internal standards drift together, so the
   expected position of each standard on a run is, by default, its
   reference position moved by the median offset of the *other* standards
   found there with a confident spectrum. A whole run shifted by a minute
   therefore still passes tight position thresholds, and a standard never
   vouches for its own position.
5. **Target table** (`references_to_targets`). The labelled standards become
   `role="istd"` rows, so the same peaks anchor the phase and normalise 🧪.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from .align import AlignmentError, AnchorPoint, fit_alignment
from .identify import IdentifySettings, build_inhouse_entry, build_library, search
from .targets import Target

PRESET_LIBRARY_PATH = Path(__file__).resolve().parents[1] / "presets" / "istd_library.json"

# The candidate/reference search itself should never drop a hit -- thresholding
# is `find_istd_on_run`'s own job -- so every prefilter `identify.search` would
# otherwise apply is switched off here.
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


FOUND = "found"
NOT_FOUND = "not found"
NO_PEAKS = "no peaks"


@dataclass
class IstdMatchSettings:
    """The three thresholds a peak must pass to be taken as the standard.

    Defaults, from the four real runs of September 2026 (labelled on one run,
    drift followed): correct matches sat within 1.5 s in t1 and 0.1 s in t2
    of the expected position, with spectral matches of 0.97 to 0.99."""
    max_dt1: float = 30.0               # s either side of the expected first-dimension time
    max_dt2: float = 0.5                # s either side of the expected second-dimension time
    min_spectral: float = 0.80          # weighted dot-product match to the reference spectrum, 0 to 1
    follow_drift: bool = True           # expect each standard where the other standards' offset puts it
    # Only for estimating that offset: a wide box, and a spectrum confident on its own.
    drift_t1_window: float = 180.0
    drift_t2_window: float = 2.5
    drift_min_spectral: float = 0.90

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "IstdMatchSettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


@dataclass
class IstdMatch:
    """One standard on one run. When it is not found, t1/t2 and the scores
    describe the closest miss (the best spectrum inside the position box), if
    there was one, and `reason` says which threshold it failed."""
    run: str
    name: str
    status: str = NOT_FOUND             # FOUND, NOT_FOUND or NO_PEAKS
    peak_id: object = None              # component id (MS-aware) or blob index (legacy)
    t1: float = float("nan")
    t2: float = float("nan")
    dt1: float = float("nan")           # from the reference peak
    dt2: float = float("nan")
    expected_t1: float = float("nan")   # the reference position, moved by the drift if followed
    expected_t2: float = float("nan")
    spectral_score: float = 0.0
    forward: float = 0.0
    reverse: float = 0.0
    n_passing: int = 0                  # peaks passing all three thresholds; more than one is worth a look
    n_drift_anchors: int = 0            # other standards the expected position came from
    reason: str = ""
    peak_confidence: str = ""
    peak_unresolved: bool = False

    @property
    def found(self) -> bool:
        return self.status == FOUND

    @property
    def needs_review(self) -> bool:
        """A passing spectral match can still belong to an uncertain peak fit."""
        return self.found and (self.n_passing > 1 or self.peak_unresolved or self.peak_confidence == "low")

    @property
    def off_t1(self) -> float:
        """Distance from the expected first-dimension time: what `max_dt1` is tested against."""
        return self.t1 - self.expected_t1

    @property
    def off_t2(self) -> float:
        return self.t2 - self.expected_t2

    @classmethod
    def from_dict(cls, d: dict) -> "IstdMatch":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


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
    """`reference` on one run's `peaks` (a `BlobSet` or a `ComponentPeaks`:
    anything with `t1`, `t2`, `mz` and `spectra`), by the three thresholds.

    The position box is centred on `expected`, or on the reference position if
    none is given. Of the peaks inside it whose spectrum passes `min_spectral`,
    the best spectral match is taken. If none passes, the result is NOT_FOUND
    and carries the closest miss so the table can say why."""
    e1, e2 = expected if expected is not None else (reference.t1, reference.t2)
    base = dict(run=run_name, name=reference.name, expected_t1=float(e1), expected_t2=float(e2),
                n_drift_anchors=n_drift_anchors)
    if len(peaks) == 0:
        return IstdMatch(status=NO_PEAKS, reason="no peaks detected on this run", **base)
    t1 = np.asarray(peaks.t1, dtype=float)
    t2 = np.asarray(peaks.t2, dtype=float)
    box = np.flatnonzero((np.abs(t1 - e1) <= settings.max_dt1) & (np.abs(t2 - e2) <= settings.max_dt2))
    if box.size == 0:
        return IstdMatch(reason=f"no peak within {settings.max_dt1:g} s in t1 and {settings.max_dt2:g} s in t2",
                         **base)
    scores = _spectral_scores(reference, peaks, box)
    best = int(np.argmax([s for s, _f, _r in scores]))
    i = int(box[best])
    spectral, forward, reverse = scores[best]
    n_passing = sum(s >= settings.min_spectral for s, _f, _r in scores)
    if n_passing:
        status = FOUND
        reason = "" if n_passing == 1 else f"{n_passing} peaks pass; the best spectral match was taken"
    else:
        status = NOT_FOUND
        reason = (f"{box.size} peak(s) in position, best spectral match {spectral:.2f} "
                  f"below {settings.min_spectral:.2f}")
    confidence = getattr(peaks, "confidence", [])
    unresolved = getattr(peaks, "unresolved", [])
    peak_confidence = str(confidence[i]) if len(confidence) > i else ""
    peak_unresolved = bool(unresolved[i]) if len(unresolved) > i else False
    if status == FOUND and (peak_unresolved or peak_confidence == "low"):
        warning = "Peak fit is unresolved or has low confidence; review this match"
        reason = f"{reason}. {warning}" if reason else warning
    return IstdMatch(status=status, peak_id=_peak_id(peaks, i), t1=float(t1[i]), t2=float(t2[i]),
                     dt1=float(t1[i] - reference.t1), dt2=float(t2[i] - reference.t2),
                     spectral_score=float(spectral), forward=float(forward), reverse=float(reverse),
                     n_passing=int(n_passing), reason=reason, peak_confidence=peak_confidence,
                     peak_unresolved=peak_unresolved, **base)


def match_run(references: list[IstdReference], run_name: str, peaks,
              settings: IstdMatchSettings = IstdMatchSettings()) -> list[IstdMatch]:
    """Every reference on one run.

    With `follow_drift`, each standard is first looked for in a wide box, and
    kept as a drift anchor if its best spectral match alone reaches
    `drift_min_spectral`. Each standard is then expected at its reference
    position moved by the median t1/t2 offset of the *other* anchors
    (`core.align.fit_alignment(..., "shift")`; a median shift, not an affine
    fit, because with a handful of anchors, one possibly wrong, anything more
    flexible extrapolates wildly). With no other anchor, or `follow_drift`
    off, it is expected at its reference position."""
    anchors = []
    if settings.follow_drift and len(peaks):
        t1 = np.asarray(peaks.t1, dtype=float)
        t2 = np.asarray(peaks.t2, dtype=float)
        for ref in references:
            box = np.flatnonzero((np.abs(t1 - ref.t1) <= settings.drift_t1_window)
                                 & (np.abs(t2 - ref.t2) <= settings.drift_t2_window))
            if box.size == 0:
                continue
            scores = [s for s, _f, _r in _spectral_scores(ref, peaks, box)]
            k = int(np.argmax(scores))
            if scores[k] >= settings.drift_min_spectral:
                anchors.append(AnchorPoint(ref.name, ref.t1, ref.t2, float(t1[box[k]]), float(t2[box[k]])))

    out = []
    for ref in references:
        others = [a for a in anchors if a.name != ref.name]
        expected = None
        if others:
            try:
                p1, p2 = fit_alignment(others, "shift").inverse(ref.t1, ref.t2)
                expected = (float(p1), float(p2))
            except AlignmentError:
                expected = None
        out.append(find_istd_on_run(ref, run_name, peaks, settings, expected=expected,
                                    n_drift_anchors=len(others) if expected is not None else 0))
    return out


# ------------------------------------------------------------------ to the target table
def references_to_targets(references: list[IstdReference], n_quals: int = 2, min_mz: int = 50) -> list[Target]:
    """Target-table rows (role `istd`) for labelled internal standards, so the
    same peaks anchor the phase and normalise quantification.

    Ions are ranked by intensity times m/z, so a strong high-mass ion (a
    molecular ion) beats a slightly stronger low-mass fragment that many
    compounds share; ions below `min_mz` are left out. The best is the
    quantifier, the next `n_quals` the qualifiers, with their measured ratios.
    Positions are the labelled peak's own. Check the ions before relying on
    them: an ion shared with a co-eluting compound makes a poor quantifier."""
    out = []
    for ref in references:
        spec = np.clip(np.asarray(ref.spectrum, dtype=float), 0.0, None)
        mz = ref.mz0 + np.arange(spec.size)
        order = [j for j in np.argsort(spec * mz)[::-1] if mz[j] >= min_mz and spec[j] > 0][:n_quals + 1]
        if not order:
            continue
        q, quals = order[0], order[1:]
        out.append(Target(compound=ref.name, role="istd", quant_ion=int(mz[q]),
                          qual_ions=tuple(int(mz[j]) for j in quals),
                          qual_ratios=tuple(round(float(spec[j] / spec[q]), 3) for j in quals),
                          t1=round(float(ref.t1), 1), t2=round(float(ref.t2), 3)))
    return out


# ------------------------------------------------------------------ report
REPORT_FIELDS = ["run", "name", "status", "t1_min", "t2_s", "dt1_s", "dt2_s", "off_expected_t1_s",
                 "off_expected_t2_s", "spectral_match", "forward", "reverse", "n_passing", "n_drift_anchors",
                 "max_dt1_s", "max_dt2_s", "min_spectral", "reason",
                 "peak_confidence", "peak_unresolved", "needs_review"]


def build_report(matches: list[IstdMatch], settings: IstdMatchSettings = IstdMatchSettings()) -> list[dict]:
    """One row per standard per run, with the thresholds it was judged against."""
    return [dict(run=m.run, name=m.name, status=m.status, t1_min=m.t1 / 60.0, t2_s=m.t2, dt1_s=m.dt1, dt2_s=m.dt2,
                 off_expected_t1_s=m.off_t1, off_expected_t2_s=m.off_t2, spectral_match=m.spectral_score,
                 forward=m.forward, reverse=m.reverse, n_passing=m.n_passing, n_drift_anchors=m.n_drift_anchors,
                 max_dt1_s=settings.max_dt1, max_dt2_s=settings.max_dt2, min_spectral=settings.min_spectral,
                 reason=m.reason, peak_confidence=m.peak_confidence,
                 peak_unresolved=m.peak_unresolved, needs_review=m.needs_review)
            for m in matches]


def write_istd_report_csv(matches: list[IstdMatch], path, settings: IstdMatchSettings = IstdMatchSettings()) -> Path:
    import csv
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=REPORT_FIELDS)
        writer.writeheader()
        writer.writerows(build_report(matches, settings))
    return path
