"""Targeted quantification across a sequence.

The idea in plain terms
-----------------------
1. In every run, find every target and internal standard on its own
   quantifier ion (`quantify_run`). Several candidate peaks may sit in a
   target's window; each is scored on how close it is to the expected
   position, how well its qualifier-ion ratios match the standard's, and
   how clearly it stands out of the noise. The best one wins. If the
   runner-up scores almost as well, the result is flagged as ambiguous.
   Every TIC blob in the run is found as well.
2. Each target is divided by its internal standard (the one named in the
   target table, otherwise the nearest one), giving a response ratio.
3. The calibration runs give, for each compound, pairs of known amount and
   response ratio, and a curve is fitted to them (`calibration.py`).
   Calibration runs can be grouped into sets (a `cal_set` column in the
   sequence), for example one set at the start of a week and another at
   the end. A sample with its own `cal_set` uses that set exactly. Any
   other sample sits between two calibration blocks in run order, and
   (`interpolate_cal_sets`, on by default) its amount is a weighted blend
   of what each bracketing set's curve would give, the weight set by how
   far along the sample sits between them (`bracket_cal_sets`). This
   assumes drift is roughly linear in run order between calibrations,
   which is the usual reason to bracket a sequence in the first place. A
   sample outside every calibration set (before the first, or after the
   last) is flagged as extrapolated and uses the nearest set alone. With
   only one set, or with interpolation turned off, a sample uses its own
   set if given, otherwise the most recent calibration set before it.
4. Samples are read back through the curves to give amounts, and amounts
   divided by the sampled air volume give concentrations.
5. Every other peak ("all signals") is reported as its volume divided by
   the volume of the nearest internal standard in the same run. If a
   calibrated target lies nearby, the peak is also given a
   semi-quantitative amount using that target's TIC calibration. This
   assumes the unknown responds like its neighbour, which is only roughly
   true within a compound class; such amounts are labelled as
   "as <target>".

"Nearest" uses a scaled distance in which `t1_scale` seconds along the
first dimension count the same as `t2_scale` seconds along the second.

Two response kinds are calibrated for every target. The quantifier-ion
ratio is used for the target itself. The TIC ratio (whole-peak volume over
whole-IS volume) is used when the curve stands in for nearby unknowns,
because unknowns have no quantifier ion of their own.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from .background import NoiseMap
from .calibration import CalibrationError, CalibrationFit, CalibrationSettings, fit_calibration
from .background import storage_threshold
from .fold import Fold
from .peaks import (BlobSet, PeakSettings, _max_raw, detect_ceiling, find_blobs, make_frame,
                    measure_ion_peak, measure_ion_polygon)
from .run import Run
from .targets import Target


@dataclass(frozen=True)
class QuantSettings:
    peaks: PeakSettings = field(default_factory=PeakSettings)
    calibration: CalibrationSettings = field(default_factory=CalibrationSettings)
    qual_tolerance: float = 20.0      # % relative deviation allowed on each qualifier ratio
    t1_tolerance: float = 4.0         # s; how far from the expected t1 a peak may sit (1 SD in scoring)
    t2_tolerance: float = 0.04        # s; the same along t2
    use_istd_offsets: bool = True     # move expected positions with the paired internal standard
    learn_positions: bool = True      # take expected positions from the calibration runs
    ambiguity: float = 0.8            # runner-up/best score above this is flagged
    t1_scale: float = 60.0            # s along t1 that count as one distance unit
    t2_scale: float = 0.5             # s along t2 that count as one distance unit
    surrogate_max_distance: float = 3.0
    saturation: float | None = None   # detector ceiling; None = find it in the data
    is_tolerance: float = 50.0        # % from the median before an IS area is flagged
    interpolate_cal_sets: bool = True  # blend bracketing calibration sets by sequence position
    # OPTIONAL, off by default. When a target's quantifier ion is saturated, quantify it from the
    # first unsaturated qualifier ion instead: volume(qualifier) / expected ratio, the ratio
    # coming from the target table (put an M+1 or other isotope ion in the qualifier list, with
    # its ratio, to use it). The hit is flagged as quantified on that ion.
    saturated_qualifier: bool = False

    def to_dict(self) -> dict:
        d = asdict(self)
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "QuantSettings":
        d = dict(d or {})
        peaks = PeakSettings.from_dict(d.pop("peaks", {}))
        cal = CalibrationSettings.from_dict(d.pop("calibration", {}))
        return cls(peaks=peaks, calibration=cal,
                   **{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# ------------------------------------------------------------------ decisions
@dataclass
class QuantDecisions:
    """Everything a person decided during review. Saved with the project."""
    rejected: set = field(default_factory=set)        # {(run, compound)}
    accepted: set = field(default_factory=set)        # {(run, compound)} checked and fine
    excluded: dict = field(default_factory=dict)      # {compound: {run: reason}} left out of the curve
    overrides: dict = field(default_factory=dict)     # {compound: {weights/model/force_zero}}
    assigned: dict = field(default_factory=dict)      # {(run, compound): (t1, t2)} peak chosen by hand
    footprints: dict = field(default_factory=dict)    # {(run, compound): [(t1, t2), ...]} hand-drawn integration polygon

    def status_of(self, run: str, compound: str) -> str | None:
        if (run, compound) in self.rejected:
            return "rejected"
        if (run, compound) in self.accepted:
            return "accepted"
        return None

    def to_dict(self) -> dict:
        return {
            "rejected": sorted(list(x) for x in self.rejected),
            "accepted": sorted(list(x) for x in self.accepted),
            "excluded": self.excluded,
            "overrides": self.overrides,
            "assigned": [[r, c, t1, t2] for (r, c), (t1, t2) in self.assigned.items()],
            "footprints": [[r, c, [list(v) for v in poly]] for (r, c), poly in self.footprints.items()],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "QuantDecisions":
        d = d or {}
        return cls(rejected={tuple(x) for x in d.get("rejected", [])},
                   accepted={tuple(x) for x in d.get("accepted", [])},
                   excluded={k: dict(v) for k, v in d.get("excluded", {}).items()},
                   overrides={k: dict(v) for k, v in d.get("overrides", {}).items()},
                   assigned={(r, c): (float(t1), float(t2)) for r, c, t1, t2 in d.get("assigned", [])},
                   footprints={(r, c): [(float(a), float(b)) for a, b in poly]
                               for r, c, poly in d.get("footprints", [])})

    def copy(self) -> "QuantDecisions":
        return QuantDecisions.from_dict(self.to_dict())


# ------------------------------------------------------------------ per run
@dataclass
class TargetHit:
    compound: str
    role: str
    found: bool
    t1: float = float("nan")
    t2: float = float("nan")
    snr: float = float("nan")
    quant_volume: float = float("nan")
    tic_volume: float = float("nan")
    qual_ratios: tuple = ()
    qual_expected: tuple = ()
    qual_ok: bool | None = None
    saturated: bool = False
    ambiguous: bool = False
    shared_with: tuple = ()
    touches_edge: bool = False
    touches_larger: bool = False
    off_position: bool = False
    expected: tuple = ()           # (t1, t2) used for the search, after corrections
    score: float = 0.0
    n_candidates: int = 0
    manual: bool = False
    integrated_by_hand: bool = False   # a hand-drawn polygon replaced the automatic footprint
    quant_ion_used: int | None = None  # set when a saturated quantifier was replaced by a qualifier (optional)
    blob: int | None = None
    pixels: tuple = ()             # extended-frame (rows, cols) of the footprint
    apex: tuple = ()               # extended-frame (row, col)

    def flags(self) -> list[str]:
        if not self.found:
            return ["not found"]
        f = []
        if self.qual_ok is False:
            f.append("qualifier ratio")
        if self.saturated:
            f.append("saturated" if self.quant_ion_used is None
                     else f"saturated (quantified on m/z {self.quant_ion_used})")
        if self.integrated_by_hand:
            f.append("integrated by hand")
        if self.ambiguous:
            f.append("ambiguous")
        if self.shared_with:
            f.append("shares peak with " + ", ".join(self.shared_with))
        if self.touches_edge:
            f.append("cut by wrap boundary")
        if self.touches_larger:
            f.append("on the flank of a larger peak")
        if self.off_position:
            f.append("off position")
        return f


@dataclass
class RunQuant:
    run: str
    hits: dict                      # compound -> TargetHit
    blobs: BlobSet
    istd_for: dict                  # compound -> internal standard name
    ceiling: float | None

    def istd_hit(self, compound: str) -> TargetHit | None:
        name = self.istd_for.get(compound)
        return self.hits.get(name) if name else None

    def response(self, compound: str, kind: str = "quant") -> float:
        hit = self.hits.get(compound)
        ref = self.istd_hit(compound)
        if hit is None or ref is None or not hit.found or not ref.found:
            return float("nan")
        attr = "quant_volume" if kind == "quant" else "tic_volume"
        num, den = getattr(hit, attr), getattr(ref, attr)
        return float(num / den) if den > 0 else float("nan")

    def nearest_istd(self, t1: float, t2: float, t1_scale: float, t2_scale: float) -> TargetHit | None:
        found = [h for h in self.hits.values() if h.role == "istd" and h.found]
        if not found:
            return None
        d = [np.hypot((h.t1 - t1) / t1_scale, (h.t2 - t2) / t2_scale) for h in found]
        return found[int(np.argmin(d))]


def _istd_pairs(targets: list[Target], s: QuantSettings) -> dict:
    stds = [t for t in targets if t.is_istd]
    out = {}
    for t in targets:
        if t.is_istd:
            out[t.compound] = t.compound
        elif t.istd:
            out[t.compound] = t.istd
        elif stds:
            d = [np.hypot((i.t1 - t.t1) / s.t1_scale, (i.t2 - t.t2) / s.t2_scale) for i in stds]
            out[t.compound] = stds[int(np.argmin(d))].compound
    return out


def quantify_run(raw: Run, clean: Run, fold: Fold, noise: NoiseMap | None, targets: list[Target],
                 settings: QuantSettings = QuantSettings(), assigned: dict | None = None,
                 progress=None, footprints: dict | None = None) -> RunQuant:
    """Find every target, internal standard and TIC blob in one run.
    `assigned` maps compound -> (t1, t2) for peaks chosen by hand; `footprints`
    maps compound -> a hand-drawn polygon [(t1, t2), ...] to integrate instead
    of the automatic footprint."""
    s = settings
    assigned = assigned or {}
    footprints = footprints or {}
    step = progress or (lambda *a: None)
    step(0, 3, "Finding peaks")
    # The storage threshold (and so the noise floor) is a property of the raw file.
    blobs = find_blobs(clean, fold, s.peaks, floor=storage_threshold(raw) / 3.0)
    frame = make_frame(fold, s.peaks)
    cache: dict = {}
    tic_ext = frame.extend(fold.tic_image(clean))
    ceiling = s.saturation if s.saturation is not None else detect_ceiling(raw)
    step(1, 3, "Measuring targets")
    hits = {}
    pairs = _istd_pairs(targets, s)
    offsets: dict = {}
    # Internal standards first, so targets can follow their retention shifts.
    ordered = [t for t in targets if t.is_istd] + [t for t in targets if not t.is_istd]
    for t in ordered:
        manual = t.compound in assigned
        e1, e2 = t.t1, t.t2
        if not t.is_istd and s.use_istd_offsets and pairs.get(t.compound) in offsets:
            d1, d2 = offsets[pairs[t.compound]]
            e1, e2 = e1 + d1, e2 + d2
        c1, c2 = assigned.get(t.compound, (e1, e2))
        w1 = t.t1_window if not manual else max(fold.spec.period, 1.0)
        w2 = t.t2_window if not manual else max(3 * fold.row_pitch, 0.06)
        tol1 = s.t1_tolerance if not t.is_istd else max(s.t1_tolerance, t.t1_window / 2)
        tol2 = s.t2_tolerance if not t.is_istd else max(s.t2_tolerance, t.t2_window / 2)
        poly = footprints.get(t.compound)
        if poly:
            manual = True
            p = measure_ion_polygon(clean, raw, fold, noise, t.quant_ion, t.qual_ions, poly, s.peaks,
                                    frame=frame, cache=cache)
            cands = [p] if p is not None else []
            if p is not None:
                c1, c2 = p.t1, p.t2          # the drawn peak is where it is: no position penalty
        else:
            cands = measure_ion_peak(clean, raw, fold, noise, t.quant_ion, t.qual_ions, c1, c2,
                                     w1, w2, s.peaks, frame=frame, cache=cache)
        if not cands:
            hits[t.compound] = TargetHit(t.compound, t.role, False, manual=manual,
                                         qual_expected=t.qual_ratios, expected=(c1, c2))
            continue
        scores = []
        for p in cands:
            pos = np.exp(-0.5 * (((p.t1 - c1) / tol1) ** 2 + ((p.t2 - c2) / tol2) ** 2))
            qual = 1.0
            for obs, exp in zip(p.qual_ratios, t.qual_ratios):
                if exp > 0 and np.isfinite(obs):
                    dev = 100.0 * abs(obs / exp - 1.0)
                    qual *= np.exp(-0.5 * (dev / s.qual_tolerance) ** 2)
            size = 1.0 - np.exp(-p.snr / 10.0)
            # A floor keeps scores comparable when every candidate is far off.
            scores.append(max(pos, 1e-12) * max(qual, 1e-6) * size)
        order = np.argsort(scores)[::-1]
        best = cands[order[0]]
        ambiguous = (not manual and len(order) > 1
                     and scores[order[1]] >= s.ambiguity * scores[order[0]])
        qual_ok = None
        if t.qual_ratios:
            devs = [100.0 * abs(o / e - 1.0) for o, e in zip(best.qual_ratios, t.qual_ratios)
                    if e > 0 and np.isfinite(o)]
            qual_ok = bool(devs) and max(devs) <= s.qual_tolerance
        off = (not manual) and (abs(best.t1 - c1) > 3 * tol1 or abs(best.t2 - c2) > 3 * tol2)
        rows, cols = best.pixels
        saturated = bool(ceiling is not None and best.max_raw >= ceiling * 0.999)
        quant_volume, used_ion = best.quant_volume, None
        if saturated and s.saturated_qualifier:
            # Optional fallback: the first qualifier that is itself unsaturated and has an expected ratio.
            for ion, qvol, expected in zip(t.qual_ions, best.qual_volumes, t.qual_ratios):
                if expected > 0 and qvol > 0 and _max_raw(raw, frame, ion, rows, cols, cache) < ceiling * 0.999:
                    quant_volume, used_ion = float(qvol / expected), int(ion)
                    break
        hits[t.compound] = TargetHit(
            compound=t.compound, role=t.role, found=True, t1=best.t1, t2=best.t2, snr=best.snr,
            quant_volume=quant_volume, tic_volume=float(tic_ext[rows, cols].sum()),
            qual_ratios=tuple(best.qual_ratios), qual_expected=t.qual_ratios, qual_ok=qual_ok,
            saturated=saturated, quant_ion_used=used_ion, integrated_by_hand=bool(poly),
            ambiguous=bool(ambiguous), touches_edge=best.touches_edge,
            touches_larger=best.touches_larger, off_position=bool(off), expected=(c1, c2),
            score=float(scores[order[0]]), n_candidates=len(cands), manual=manual,
            blob=blobs.at_pixel(best.row, best.col), pixels=(rows, cols), apex=(best.row, best.col))
        if t.is_istd and not off:
            offsets[t.compound] = (best.t1 - t.t1, best.t2 - t.t2)

    # Targets on the same ion that picked the same peak cannot both be right.
    step(2, 3, "Checking overlaps")
    by_peak: dict = {}
    for t in targets:
        h = hits[t.compound]
        if h.found:
            by_peak.setdefault((t.quant_ion, h.apex), []).append(t.compound)
    for names in by_peak.values():
        if len(names) > 1:
            for n in names:
                h = hits[n]
                h.shared_with = tuple(x for x in names if x != n)
    step(3, 3, "Done")
    hits = {t.compound: hits[t.compound] for t in targets}     # table order
    return RunQuant(run=raw.name, hits=hits, blobs=blobs, istd_for=pairs, ceiling=ceiling)


def learn_positions(seq: list["SequenceRun"], targets: list[Target],
                    settings: QuantSettings = QuantSettings()) -> dict:
    """Expected positions learned from the data.

    Internal standards: the median position over every run where they were
    found cleanly. Targets: the median over calibration runs, after removing
    each run's internal-standard shift, so the result is in the same frame as
    the table. Only clean hits count (no flags). Returns {compound: (t1, t2)}.
    """
    out = {}
    pairs = _istd_pairs(targets, settings)
    for t in targets:
        pts = []
        for r in seq:
            q = r.quant
            if q is None:
                continue
            h = q.hits.get(t.compound)
            if h is None or not h.found or h.flags():
                continue
            if t.is_istd:
                pts.append((h.t1, h.t2))
            elif r.is_calibration:
                ref = q.hits.get(pairs.get(t.compound, ""))
                d1 = d2 = 0.0
                if settings.use_istd_offsets and ref is not None and ref.found:
                    it = next(x for x in targets if x.compound == ref.compound)
                    d1, d2 = ref.t1 - it.t1, ref.t2 - it.t2
                pts.append((h.t1 - d1, h.t2 - d2))
        if pts:
            arr = np.array(pts)
            out[t.compound] = (float(np.median(arr[:, 0])), float(np.median(arr[:, 1])))
    return out


def with_positions(targets: list[Target], positions: dict) -> list[Target]:
    from dataclasses import replace
    return [replace(t, t1=positions[t.compound][0], t2=positions[t.compound][1])
            if t.compound in positions else t for t in targets]


# ------------------------------------------------------------- sequence
@dataclass
class SequenceRun:
    name: str
    kind: str = ""
    level: str = ""
    cal_set: str = ""
    volume: float | None = None         # sampled air volume (m³)
    quant: RunQuant | None = None

    @property
    def is_calibration(self) -> bool:
        return "cal" in (self.kind or "").lower() and bool(self.level)


def read_sequence_extras(path) -> dict:
    """cal_set and volume_m3 columns of a sequence CSV, keyed by run name."""
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            name = Path(row.get("file", "")).stem
            vol = row.get("volume_m3")
            out[name] = {"cal_set": (row.get("cal_set") or "").strip(),
                         "volume": float(vol) if vol not in (None, "") else None}
    return out


def cal_set_for(seq: list[SequenceRun], index: int) -> str:
    """The single calibration set for a run: its own, otherwise the most
    recent calibration set before it, otherwise the first one there is.
    Used when interpolation is off, or as the label for an exact match."""
    run = seq[index]
    if run.cal_set:
        return run.cal_set
    for j in range(index, -1, -1):
        if seq[j].is_calibration and seq[j].cal_set:
            return seq[j].cal_set
    sets = [r.cal_set for r in seq if r.is_calibration]
    return sets[0] if sets else ""


def _set_positions(seq: list[SequenceRun]) -> dict[str, float]:
    """Each calibration set's position in the sequence: the mean index of
    its own calibration runs (levels usually run together, so this is close
    to the middle of that block)."""
    idx: dict[str, list[int]] = {}
    for i, r in enumerate(seq):
        if r.is_calibration and r.cal_set:
            idx.setdefault(r.cal_set, []).append(i)
    return {cs: float(np.mean(v)) for cs, v in idx.items()}


@dataclass
class CalBracket:
    """Which calibration set (or two) a run should be read against."""
    before: str             # set to use alone, or the earlier of a pair
    after: str               # equal to `before` when only one set applies
    frac: float = 0.0        # 0 uses `before` alone, 1 uses `after` alone
    extrapolated: bool = False  # the run sits outside every calibration set

    @property
    def blended(self) -> bool:
        return self.before != self.after and 0.0 < self.frac < 1.0

    def label(self) -> str:
        if not self.before:
            return ""
        if self.before == self.after:
            return self.before + (" (extrapolated)" if self.extrapolated else "")
        return f"{self.before}→{self.after} ({self.frac * 100:.0f}%)"


def bracket_cal_sets(seq: list[SequenceRun], index: int, interpolate: bool = True) -> CalBracket:
    """Which calibration set(s) apply to the run at `index`.

    A run with its own `cal_set` uses that set exactly. Otherwise, with two
    or more calibration sets and interpolation on, the two nearest sets by
    sequence position (one at or before, one at or after) are returned with
    a blend fraction; a run before the first set or after the last uses
    that single nearest set, flagged as extrapolated. With interpolation
    off, or only one set, this matches `cal_set_for`.
    """
    run = seq[index]
    positions = _set_positions(seq)
    if run.cal_set and run.cal_set in positions:
        return CalBracket(run.cal_set, run.cal_set, 0.0)
    if not positions:
        return CalBracket("", "")
    if not interpolate or len(positions) == 1:
        cs = cal_set_for(seq, index)
        return CalBracket(cs, cs, 0.0)
    items = sorted(positions.items(), key=lambda kv: kv[1])
    before = [it for it in items if it[1] <= index]
    after = [it for it in items if it[1] >= index]
    if before and after and before[-1][0] != after[0][0]:
        (b_set, b_pos), (a_set, a_pos) = before[-1], after[0]
        frac = float(np.clip((index - b_pos) / (a_pos - b_pos), 0.0, 1.0))
        return CalBracket(b_set, a_set, frac)
    if before:
        return CalBracket(before[-1][0], before[-1][0], 0.0, extrapolated=not after)
    return CalBracket(after[0][0], after[0][0], 0.0, extrapolated=True)


def _predict(fit, resp: float) -> tuple[float, float]:
    if not isinstance(fit, CalibrationFit) or not np.isfinite(resp):
        return float("nan"), float("nan")
    a, e = fit.amount_with_se(np.array([resp]))
    return float(a[0]), float(e[0])


def blended_amount(fits: dict, compound: str, kind: str, bracket: CalBracket, resp: float):
    """Amount and its standard error at one calibration bracket.

    A single set is read directly. A pair is read from both curves and
    blended by `frac` (linear interpolation of the two amounts, and their
    errors combined in quadrature with the same weights); this is the
    bracketing-standard approach used to correct for drift between two
    calibration runs. Extra flags note when the blend leans entirely on
    one side because the other set has no fit for this compound, and the
    LOD/LOQ/top-of-curve checks use whichever set carries more weight.
    Returns (amount, se, flags, note).
    """
    fit_b = fits.get((bracket.before, compound), {}).get(kind)
    flags: list[str] = []
    note = ""
    if bracket.before == bracket.after:
        amount, se = _predict(fit_b, resp)
        if isinstance(fit_b, str):
            note = fit_b
        elif isinstance(fit_b, CalibrationFit) and np.isfinite(amount):
            flags += _range_flags(fit_b, amount)
        if bracket.extrapolated and bracket.before:
            flags.append("extrapolated calibration")
        return amount, se, flags, note
    fit_a = fits.get((bracket.after, compound), {}).get(kind)
    ab, eb = _predict(fit_b, resp)
    aa, ea = _predict(fit_a, resp)
    if not (np.isfinite(ab) or np.isfinite(aa)):
        note = fit_b if isinstance(fit_b, str) else (fit_a if isinstance(fit_a, str) else "")
        return float("nan"), float("nan"), flags, note
    if not np.isfinite(ab):
        flags.append(f"only the {bracket.after} calibration has this compound")
        return aa, ea, flags + _range_flags(fit_a, aa), note
    if not np.isfinite(aa):
        flags.append(f"only the {bracket.before} calibration has this compound")
        return ab, eb, flags + _range_flags(fit_b, ab), note
    f = bracket.frac
    amount = (1.0 - f) * ab + f * aa
    se = float(np.hypot((1.0 - f) * eb, f * ea))
    ref = fit_b if f <= 0.5 else fit_a
    return amount, se, _range_flags(ref, amount), note


def _range_flags(fit: CalibrationFit, amount: float) -> list[str]:
    flags = []
    if amount < fit.loq:
        flags.append("below LOQ" if amount >= fit.lod else "below LOD")
    if amount > fit.x.max() * 1.05:
        flags.append("above top level")
    return flags


def calibration_points(seq: list[SequenceRun], compound: str, levels: dict, cal_set: str,
                       kind: str, decisions: QuantDecisions):
    x, y, labels = [], [], []
    for r in seq:
        if not r.is_calibration or r.quant is None or (r.cal_set or "") != cal_set:
            continue
        amount = levels.get(r.level, {}).get(compound)
        if amount is None:
            continue
        resp = r.quant.response(compound, kind)
        if decisions.status_of(r.name, compound) == "rejected":
            resp = float("nan")
        x.append(float(amount))
        y.append(resp)
        labels.append(r.name)
    return np.array(x), np.array(y), labels


def build_calibrations(seq: list[SequenceRun], targets: list[Target], levels: dict,
                       settings: QuantSettings, decisions: QuantDecisions) -> dict:
    """{(cal_set, compound): {"quant": fit or error text, "tic": fit or error text}}"""
    sets = sorted({r.cal_set or "" for r in seq if r.is_calibration})
    out = {}
    for cs in sets:
        for t in targets:
            if t.is_istd:
                continue
            base = settings.calibration
            ov = decisions.overrides.get(t.compound, {})
            cal = CalibrationSettings.from_dict({**base.to_dict(), **ov})
            entry = {}
            for kind in ("quant", "tic"):
                x, y, labels = calibration_points(seq, t.compound, levels, cs, kind, decisions)
                try:
                    entry[kind] = fit_calibration(t.compound, x, y, labels,
                                                  decisions.excluded.get(t.compound, {}), cal)
                except (CalibrationError, np.linalg.LinAlgError) as exc:
                    entry[kind] = str(exc)
            out[(cs, t.compound)] = entry
    return out


@dataclass
class SampleResult:
    run: str
    compound: str
    role: str
    found: bool
    status: str               # ok, check, rejected, accepted, not found
    flags: list
    t1: float
    t2: float
    response: float
    amount: float
    amount_se: float
    concentration: float
    cal_set: str
    note: str = ""


def quantify_sequence(seq: list[SequenceRun], targets: list[Target], fits: dict,
                      decisions: QuantDecisions,
                      settings: QuantSettings = QuantSettings()) -> list[SampleResult]:
    rows = []
    for i, r in enumerate(seq):
        if r.quant is None:
            continue
        bracket = bracket_cal_sets(seq, i, settings.interpolate_cal_sets)
        for t in targets:
            h = r.quant.hits.get(t.compound)
            if h is None:
                continue
            flags = h.flags()
            decision = decisions.status_of(r.name, t.compound)
            if decision == "rejected":
                status = "rejected"
            elif not h.found:
                status = "not found"
            elif decision == "accepted":
                status = "accepted"
            else:
                status = "check" if flags else "ok"
            resp = r.quant.response(t.compound, "quant")
            amount = se = float("nan")
            note = ""
            if not t.is_istd and status not in ("rejected", "not found") and np.isfinite(resp):
                amount, se, extra, note = blended_amount(fits, t.compound, "quant", bracket, resp)
                flags = flags + extra
            if status == "ok" and len(flags) > len(h.flags()):
                status = "check"
            conc = amount / r.volume if (r.volume and np.isfinite(amount)) else float("nan")
            rows.append(SampleResult(r.name, t.compound, t.role, h.found, status, flags,
                                     h.t1, h.t2, resp, amount, se, conc, bracket.label(), note))
    return rows


@dataclass
class SignalRow:
    run: str
    blob: int
    t1: float
    t2: float
    snr: float
    volume: float
    istd: str
    istd_volume: float
    ratio: float                  # volume / nearest IS volume
    target: str                   # target whose peak this is, if any
    surrogate: str                # calibrated target used for the semi-quantitative amount
    distance: float
    amount: float                 # "as <surrogate>"
    concentration: float
    note: str = ""


def all_signals(seq: list[SequenceRun], targets: list[Target], fits: dict,
                settings: QuantSettings) -> list[SignalRow]:
    s = settings
    rows = []
    calibrated = [t for t in targets if not t.is_istd]
    for i, r in enumerate(seq):
        q = r.quant
        if q is None:
            continue
        bracket = bracket_cal_sets(seq, i, s.interpolate_cal_sets)
        owner = {h.blob: h.compound for h in q.hits.values() if h.found and h.blob is not None}
        for b in range(len(q.blobs)):
            t1, t2, vol = float(q.blobs.t1[b]), float(q.blobs.t2[b]), float(q.blobs.volume[b])
            ref = q.nearest_istd(t1, t2, s.t1_scale, s.t2_scale)
            ratio = vol / ref.tic_volume if (ref is not None and ref.tic_volume > 0) else float("nan")
            sur, dist, amount = "", float("nan"), float("nan")
            best_d = np.inf
            for t in calibrated:
                h = q.hits.get(t.compound)
                tt1, tt2 = (h.t1, h.t2) if (h is not None and h.found) else (t.t1, t.t2)
                d = float(np.hypot((tt1 - t1) / s.t1_scale, (tt2 - t2) / s.t2_scale))
                fit = fits.get((bracket.before, t.compound), {}).get("tic")
                if fit is None or not isinstance(fit, CalibrationFit):
                    fit = fits.get((bracket.after, t.compound), {}).get("tic")
                if d < best_d and isinstance(fit, CalibrationFit):
                    best_d, sur = d, t.compound
            if sur and best_d <= s.surrogate_max_distance and ref is not None:
                # The surrogate curve is in units of "target TIC / its IS TIC".
                sref = q.istd_hit(sur)
                if sref is not None and sref.found and sref.tic_volume > 0:
                    amount, _se, _flags, _note = blended_amount(fits, sur, "tic", bracket,
                                                                vol / sref.tic_volume)
                dist = best_d
            else:
                sur = ""
            note = ""
            if np.isfinite(amount) and amount < 0:
                amount, note = float("nan"), "below the surrogate curve"
            conc = amount / r.volume if (r.volume and np.isfinite(amount)) else float("nan")
            rows.append(SignalRow(r.name, b, t1, t2, float(q.blobs.snr[b]), vol,
                                  ref.compound if ref else "", ref.tic_volume if ref else float("nan"),
                                  ratio, owner.get(b, ""), sur, dist, amount, conc, note))
    return rows


@dataclass
class IstdRow:
    run: str
    compound: str
    found: bool
    t1: float
    t2: float
    volume: float
    rel_to_median: float           # %
    flag: bool


def istd_table(seq: list[SequenceRun], targets: list[Target], settings: QuantSettings) -> list[IstdRow]:
    rows = []
    for t in (x for x in targets if x.is_istd):
        vols = [r.quant.hits[t.compound].quant_volume for r in seq
                if r.quant is not None and r.quant.hits.get(t.compound) is not None
                and r.quant.hits[t.compound].found]
        med = float(np.median(vols)) if vols else float("nan")
        for r in seq:
            if r.quant is None or t.compound not in r.quant.hits:
                continue
            h = r.quant.hits[t.compound]
            rel = 100.0 * (h.quant_volume / med - 1.0) if (h.found and med > 0) else float("nan")
            flag = (not h.found and "blank" not in (r.kind or "").lower()) or \
                (np.isfinite(rel) and abs(rel) > settings.is_tolerance)
            rows.append(IstdRow(r.name, t.compound, h.found, h.t1, h.t2, h.quant_volume, rel, bool(flag)))
    return rows


def write_csv(rows, path) -> Path:
    path = Path(path)
    if not rows:
        path.write_text("")
        return path
    fields = list(asdict(rows[0]).keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            d = asdict(r)
            for k, v in d.items():
                if isinstance(v, (list, tuple)):
                    d[k] = "; ".join(map(str, v))
            w.writerow(d)
    return path


def calibration_table(fits: dict) -> list[dict]:
    out = []
    for (cs, compound), entry in sorted(fits.items()):
        for kind, fit in entry.items():
            if isinstance(fit, CalibrationFit):
                out.append(dict(cal_set=cs, compound=compound, response=kind, model=fit.model,
                                weights=fit.weights, zero=fit.zero, intercept=fit.intercept,
                                slope=fit.slope, curvature=fit.curvature, r2=fit.r2,
                                mandel_p=fit.mandel_p, intercept_p=fit.intercept_p,
                                worst_rel_error=float(np.max(np.abs(fit.rel_error))),
                                passes=fit.ok, lod=fit.lod, loq=fit.loq, n=fit.n,
                                excluded="; ".join(f"{k}: {v}" for k, v in fit.excluded.items())))
            else:
                out.append(dict(cal_set=cs, compound=compound, response=kind, model="",
                                weights="", zero="", intercept="", slope="", curvature="", r2="",
                                mandel_p="", intercept_p="", worst_rel_error="", passes=False,
                                lod="", loq="", n=0, excluded=str(fit)))
    return out


def write_dicts(rows: list[dict], path) -> Path:
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        if rows:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
    return path


# ------------------------------------------------------------ whole sequence
@dataclass
class QuantReport:
    """Everything quantification produced for a sequence."""
    seq: list
    targets: list                 # with learned positions, if any
    positions: dict               # learned {compound: (t1, t2)}
    fits: dict
    results: list
    signals: list
    istd: list
    settings: QuantSettings

    def fit(self, compound: str, cal_set: str = "", kind: str = "quant"):
        return self.fits.get((cal_set, compound), {}).get(kind)

    def write(self, outdir) -> list[Path]:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        return [
            write_csv(self.results, outdir / "targets.csv"),
            write_dicts(calibration_table(self.fits), outdir / "calibrations.csv"),
            write_csv(self.signals, outdir / "all_signals.csv"),
            write_csv(self.istd, outdir / "internal_standards.csv"),
        ]


def finish_sequence(seq: list, targets: list, levels: dict, settings: QuantSettings,
                    decisions: QuantDecisions, positions: dict | None = None) -> QuantReport:
    """Calibrate and quantify once every run has its RunQuant."""
    fits = build_calibrations(seq, targets, levels, settings, decisions)
    return QuantReport(seq=seq, targets=targets, positions=positions or {}, fits=fits,
                       results=quantify_sequence(seq, targets, fits, decisions, settings),
                       signals=all_signals(seq, targets, fits, settings),
                       istd=istd_table(seq, targets, settings), settings=settings)


def quantify_prepared(prepared: list, info: list, targets: list, levels: dict,
                      settings: QuantSettings = QuantSettings(),
                      decisions: QuantDecisions | None = None, progress=None) -> QuantReport:
    """Quantify runs that are already folded and cleaned.

    prepared: [(raw Run, cleaned Run, Fold, NoiseMap or None)]
    info:     [SequenceRun without quant] in the same order
    With `learn_positions`, the calibration runs are measured first,
    expected positions are learned from them, and then every run is
    measured with the learned positions.
    """
    energies = {raw.meta.get('energy_ev') for raw, *_ in prepared if raw.meta.get('energy_ev') is not None}
    if len(energies) > 1:
        raise ValueError('Quantification cannot mix ionisation energies; select one energy channel for the sequence')
    decisions = decisions or QuantDecisions()
    step = progress or (lambda *a: None)

    def measure(tlist, which):
        for k, i in enumerate(which):
            raw, clean, fold, noise = prepared[i]
            step(k, len(which), f"Measuring {info[i].name}")
            assigned = {c: v for (r, c), v in decisions.assigned.items() if r == info[i].name}
            polys = {c: v for (r, c), v in decisions.footprints.items() if r == info[i].name}
            info[i].quant = quantify_run(raw, clean, fold, noise, tlist, settings, assigned, footprints=polys)

    positions = {}
    if settings.learn_positions:
        cal = [i for i, r in enumerate(info) if r.is_calibration] or list(range(len(info)))
        measure(targets, cal)
        positions = learn_positions([info[i] for i in cal], targets, settings)
        targets = with_positions(targets, positions)
    measure(targets, list(range(len(info))))
    return finish_sequence(info, targets, levels, settings, decisions, positions)
