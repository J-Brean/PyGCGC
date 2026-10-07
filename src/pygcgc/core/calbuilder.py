"""📈 Calibration builder: curves from a set of standard runs, made by hand.

The idea in plain terms
-----------------------
A calibration set is a handful of runs of standard mixtures. For each
compound we know roughly where it elutes (t1, t2) and how much of it went
into each run. The compounds come from a table or from clicks on the
chromatogram, and some of them are internal standards.

1. Every run is cleaned (🧹) and its blobs found (🔍, TIC blobs).
2. Each compound is matched to one blob per run: the blob whose footprint
   contains the compound's position, otherwise the nearest blob apex inside
   a search ellipse (`t1_window_s` by `t2_window_s`). A match can be changed
   by hand on any run ("pinned" to another blob), or removed, in which case
   that run gives no point for the compound.
3. A compound's signal is its blob volume (every ion summed over the
   footprint), or the volume on its own quantifier ion when one is given.
4. Dividing by the signal of an internal standard in the same run gives the
   IS-normalised response, which cancels run-to-run changes in injection
   volume and detector sensitivity. The internal standard is the one named
   for the compound, otherwise the nearest one by a scaled distance in which
   `t1_scale` seconds of t1 count the same as `t2_scale` seconds of t2 (the
   same rule and defaults as 🧪).
5. Each compound's curve is response ratio against amount (or against the
   amount ratio to its internal standard, when the internal standard was not
   added at the same amount to every run), fitted by
   `calibration.fit_calibration`: weights chosen by back-calculation error,
   straight or quadratic by Mandel's test. Points can be left out of a curve
   and put back again; the reason is recorded.

Nothing here changes a run or a blob set; a session is a small record
(compounds, pins, exclusions, curve choices) replayed on top of whatever
detection produced, so it survives recomputing the blobs.
"""

from __future__ import annotations

import csv
import json
import re
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np

from .calibration import CalibrationError, CalibrationFit, CalibrationSettings, fit_calibration
from .peaks import BlobSet, PeakSettings, make_frame


@dataclass
class ComponentBlobs(BlobSet):
    """MS-aware components (`peaktable.ComponentPeaks`) as the builder reads peaks.

    Components overlap, so there is no footprint image: `labels` is all zero, a
    compound is matched to the nearest component apex, and outlines are the fitted
    10% contours. Everything else (t1, t2, volume, full-ion spectra) reads as a blob's."""
    contours: list = field(default_factory=list)

    def outline_times(self, i: int):
        if i < len(self.contours) and self.contours[i] is not None:
            return self.contours[i]
        return np.array([]), np.array([])


def peaks_from_components(peaks, fold, wrap_below_s: float = 0.5) -> ComponentBlobs:
    """`ComponentPeaks` from `peaktable.respectrum_components` wrapped as `ComponentBlobs`."""
    frame = make_frame(fold, PeakSettings(wrap_below_s=wrap_below_s))
    t1 = np.asarray(peaks.t1, dtype=float)
    t2 = np.asarray(peaks.t2, dtype=float)
    n = t1.size
    col = np.floor((t1 - frame.fold.t1[0]) / frame.period).astype(int)
    row = np.round(t2 / frame.pitch).astype(int)
    volume = np.asarray(peaks.volume, dtype=float)
    return ComponentBlobs(
        mz=np.asarray(peaks.mz), t1=t1, t2=t2, col=col, row=row, height=volume.copy(), snr=np.zeros(n),
        volume=volume, n_pixels=np.zeros(n, dtype=int), n_modulations=np.zeros(n, dtype=int),
        spectra=np.asarray(peaks.spectra, dtype=np.float32),
        labels=np.zeros((frame.rows, fold.n_cols), dtype=np.int32), wrap_rows=frame.wrap_rows,
        pitch=frame.pitch, period=frame.period, t1_start=float(fold.t1[0]),
        co_elution=np.asarray(peaks.co_elution, dtype=bool), contours=list(peaks.contours))

ISTD_WORDS = {"istd", "is", "internal standard", "internal_standard", "surrogate"}
TRUE_WORDS = {"1", "true", "yes", "y", "x"}


@dataclass(frozen=True)
class CalBuilderSettings:
    t1_window_s: float = 20.0      # search half-width along t1 for an automatic match
    t2_window_s: float = 0.25      # ...and along t2
    t1_scale: float = 60.0         # s of t1 that count as one unit of "nearest IS" distance
    t2_scale: float = 0.5          # s of t2 that count as one unit
    amount_ratio: bool = False     # x = amount / IS amount, for an IS not spiked equally into every run
    saturation: float = 0.0        # a point whose largest raw value reaches this is flagged (0: never)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "CalBuilderSettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


@dataclass
class CalCompound:
    name: str
    t1: float                       # s, first-dimension apex (the blob table's convention)
    t2: float                       # s, true second-dimension time
    is_istd: bool = False
    quant_ion: int | None = None    # None: the whole blob (every ion)
    istd: str = ""                  # internal standard to divide by; blank for the nearest
    amounts: dict[str, float] = field(default_factory=dict)   # run name -> amount

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "CalCompound":
        d = dict(d)
        d["amounts"] = {str(k): float(v) for k, v in (d.get("amounts") or {}).items()}
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class CalPoint:
    """One compound in one run."""
    compound: str
    run: str
    amount: float                   # NaN when none is given for this run
    blob: int | None                # index into that run's BlobSet
    pinned: bool                    # matched by hand
    t1: float                       # the matched blob's apex (NaN if none)
    t2: float
    signal: float
    istd: str                       # internal standard divided by ("" for none)
    istd_signal: float
    ratio: float                    # signal / istd_signal
    x: float                        # amount, or amount / IS amount
    height: float = float("nan")    # largest raw single-ion value in the blob, when known
    saturated: bool = False         # height at or above the saturation setting


# ------------------------------------------------------------------ matching
def blob_pixel(blobs: BlobSet, t1: float, t2: float) -> tuple[int, int]:
    """(row, column) of the extended-frame pixel at a true (t1, t2), as `Frame.time_to_pixel`."""
    n_rows = blobs.labels.shape[0] - blobs.wrap_rows
    col = int(np.floor((t1 - blobs.t1_start) / blobs.period))
    rowpos = t2 / blobs.pitch
    for _ in range(4):
        if rowpos >= blobs.labels.shape[0] - 0.5:
            rowpos -= n_rows
            col += 1
        elif rowpos < blobs.wrap_rows - 0.5:
            rowpos += n_rows
            col -= 1
        else:
            break
    return int(round(rowpos)), col


def match_blob(blobs: BlobSet | None, t1: float, t2: float,
               settings: CalBuilderSettings = CalBuilderSettings()) -> int | None:
    """The blob at (t1, t2): the one whose footprint holds that point, otherwise
    the nearest apex inside the search ellipse, otherwise None."""
    if blobs is None or len(blobs) == 0 or not (np.isfinite(t1) and np.isfinite(t2)):
        return None
    if blobs.labels.size:
        hit = blobs.at_pixel(*blob_pixel(blobs, t1, t2))
        if hit is not None:
            return hit
    d = _wrapped_distance(blobs.t1, blobs.t2, t1, t2, settings.t1_window_s, settings.t2_window_s, blobs.period)
    i = int(np.argmin(d))
    return i if d[i] <= 1.0 else None


def blob_signal(blobs: BlobSet, i: int, quant_ion: int | None = None) -> float:
    """Blob volume over every ion, or over one quantifier ion."""
    if quant_ion is None:
        return float(blobs.volume[i])
    k = int(quant_ion) - int(blobs.mz[0])
    if not 0 <= k < blobs.spectra.shape[1]:
        return float("nan")
    return float(blobs.spectra[i, k])


# ------------------------------------------------------------------ session
@dataclass
class CalSession:
    """Everything a person decided: compounds, amounts, hand matches, exclusions, curve choices."""
    runs: list[str] = field(default_factory=list)
    compounds: list[CalCompound] = field(default_factory=list)
    pins: dict = field(default_factory=dict)          # (compound, run) -> (t1, t2) or None (no peak)
    excluded: dict = field(default_factory=dict)      # compound -> {run: reason}
    overrides: dict = field(default_factory=dict)     # compound -> {weights, model, force_zero}
    settings: CalBuilderSettings = field(default_factory=CalBuilderSettings)
    calibration: CalibrationSettings = field(default_factory=CalibrationSettings)
    unit: str = ""                                    # of the amounts, for labels and exports
    levels: dict = field(default_factory=dict)        # run -> level amount, for every compound without its own

    def amount(self, compound: "CalCompound", run: str) -> float:
        """A compound's amount in a run: its own, else (not for an internal standard) the run's level."""
        if run in compound.amounts:
            return float(compound.amounts[run])
        if not compound.is_istd and run in self.levels:
            return float(self.levels[run])
        return float("nan")

    # ----------------------------------------------------------- compounds
    def compound(self, name: str) -> CalCompound | None:
        return next((c for c in self.compounds if c.name == name), None)

    def names(self) -> list[str]:
        return [c.name for c in self.compounds]

    def istds(self) -> list[CalCompound]:
        return [c for c in self.compounds if c.is_istd]

    def analytes(self) -> list[CalCompound]:
        return [c for c in self.compounds if not c.is_istd]

    def unique_name(self, base: str) -> str:
        names, name, k = set(self.names()), base, 2
        while name in names:
            name, k = f"{base} {k}", k + 1
        return name

    def add_compound(self, compound: CalCompound) -> CalCompound:
        compound.name = self.unique_name(compound.name.strip() or "compound")
        self.compounds.append(compound)
        return compound

    def remove_compound(self, name: str) -> None:
        self.compounds = [c for c in self.compounds if c.name != name]
        self.pins = {k: v for k, v in self.pins.items() if k[0] != name}
        self.excluded.pop(name, None)
        self.overrides.pop(name, None)
        for c in self.compounds:
            if c.istd == name:
                c.istd = ""

    def rename(self, old: str, new: str) -> str:
        c = self.compound(old)
        new = new.strip()
        if c is None or not new or new == old:
            return old
        new = self.unique_name(new)
        c.name = new
        self.pins = {((new if k[0] == old else k[0]), k[1]): v for k, v in self.pins.items()}
        for d in (self.excluded, self.overrides):
            if old in d:
                d[new] = d.pop(old)
        for other in self.compounds:
            if other.istd == old:
                other.istd = new
        return new

    def istd_for(self, name: str) -> CalCompound | None:
        """The internal standard a compound is divided by: its own choice, else the nearest."""
        c = self.compound(name)
        if c is None or c.is_istd:
            return None
        stds = self.istds()
        if c.istd:
            named = next((s for s in stds if s.name == c.istd), None)
            if named is not None:
                return named
        if not stds:
            return None
        s = self.settings
        d = [np.hypot((i.t1 - c.t1) / s.t1_scale, (i.t2 - c.t2) / s.t2_scale) for i in stds]
        return stds[int(np.argmin(d))]

    # -------------------------------------------------------------- matches
    def pin(self, compound: str, run: str, t1: float | None, t2: float | None = None) -> None:
        """Match a compound to the blob at (t1, t2) in one run; t1=None means no peak there."""
        self.pins[(compound, run)] = None if t1 is None else (float(t1), float(t2))

    def pin_blob(self, compound: str, run: str, blobs: BlobSet, i: int) -> None:
        """Match a compound to blob `i` of one run. The position kept is the centre of
        the blob's apex pixel, which always lies in its own footprint (its weighted
        t1 may not, for a skewed peak), so it finds the same blob again after the
        blobs are recomputed."""
        if isinstance(blobs, ComponentBlobs):
            # No footprints: the exact apex, matched back as the nearest component.
            self.pin(compound, run, float(blobs.t1[i]), float(blobs.t2[i]))
            return
        t1 = blobs.t1_start + (float(blobs.col[i]) + 0.5) * blobs.period
        self.pin(compound, run, t1, float(blobs.row[i]) * blobs.pitch)

    def unpin(self, compound: str, run: str) -> None:
        self.pins.pop((compound, run), None)

    def resolve(self, blobs_by_run: dict) -> dict:
        """(compound, run) -> blob index or None, for every compound and run."""
        out = {}
        for c in self.compounds:
            for run in self.runs:
                blobs = blobs_by_run.get(run)
                key = (c.name, run)
                if key in self.pins:
                    p = self.pins[key]
                    out[key] = None if p is None else match_blob(blobs, p[0], p[1], self.settings)
                else:
                    out[key] = match_blob(blobs, c.t1, c.t2, self.settings)
        return out

    def points(self, blobs_by_run: dict, resolved: dict | None = None, heights: dict | None = None) -> dict:
        """(compound, run) -> CalPoint for every compound and run. `heights` maps run
        name to `blob_heights` of its blobs, for the saturation flag."""
        resolved = self.resolve(blobs_by_run) if resolved is None else resolved
        raw = {}
        for c in self.compounds:
            for run in self.runs:
                i = resolved.get((c.name, run))
                blobs = blobs_by_run.get(run)
                if i is None or blobs is None:
                    raw[(c.name, run)] = (None, np.nan, np.nan, np.nan)
                else:
                    raw[(c.name, run)] = (i, float(blobs.t1[i]), float(blobs.t2[i]),
                                          blob_signal(blobs, i, c.quant_ion))
        out = {}
        nan = float("nan")
        for c in self.compounds:
            std = c if c.is_istd else self.istd_for(c.name)
            for run in self.runs:
                i, t1, t2, sig = raw[(c.name, run)]
                amount = self.amount(c, run)
                std_sig = raw[(std.name, run)][3] if std is not None else nan
                ratio = sig / std_sig if (np.isfinite(std_sig) and std_sig > 0) else nan
                x = amount
                if self.settings.amount_ratio and std is not None:
                    std_amount = float(std.amounts.get(run, nan))
                    x = amount / std_amount if std_amount > 0 else nan
                h = heights.get(run) if heights else None
                height = float(h[i]) if (h is not None and i is not None and i < len(h)) else nan
                out[(c.name, run)] = CalPoint(
                    compound=c.name, run=run, amount=amount, blob=i, pinned=(c.name, run) in self.pins,
                    t1=t1, t2=t2, signal=sig, istd=std.name if std is not None else "",
                    istd_signal=std_sig, ratio=ratio, x=x, height=height,
                    saturated=bool(self.settings.saturation > 0 and height >= self.settings.saturation))
        return out

    # ---------------------------------------------------------------- curves
    def exclude(self, compound: str, run: str, reason: str = "left out by hand") -> None:
        self.excluded.setdefault(compound, {})[run] = reason or "left out by hand"

    def include(self, compound: str, run: str) -> None:
        ex = self.excluded.get(compound, {})
        ex.pop(run, None)
        if not ex:
            self.excluded.pop(compound, None)

    def curve_settings(self, compound: str) -> CalibrationSettings:
        ov = self.overrides.get(compound, {})
        return replace(self.calibration, **{k: v for k, v in ov.items()
                                            if k in CalibrationSettings.__dataclass_fields__})

    def fit(self, compound: str, points: dict) -> CalibrationFit | str:
        """The curve for one compound, or the reason there is none."""
        pts = [points[(compound, r)] for r in self.runs if (compound, r) in points]
        if not pts:
            return "no points"
        x = np.array([p.x for p in pts])
        y = np.array([p.ratio for p in pts])
        labels = [p.run for p in pts]
        try:
            return fit_calibration(compound, x, y, labels=labels,
                                   excluded=self.excluded.get(compound, {}),
                                   settings=self.curve_settings(compound))
        except (CalibrationError, np.linalg.LinAlgError, ValueError) as exc:
            return str(exc)

    def fits(self, points: dict) -> dict:
        return {c.name: self.fit(c.name, points) for c in self.analytes()}

    # ----------------------------------------------------------------- files
    def to_dict(self) -> dict:
        return dict(runs=list(self.runs), compounds=[c.to_dict() for c in self.compounds],
                    pins=[[k[0], k[1], None if v is None else list(v)] for k, v in self.pins.items()],
                    excluded={k: dict(v) for k, v in self.excluded.items()},
                    overrides={k: dict(v) for k, v in self.overrides.items()},
                    settings=self.settings.to_dict(), calibration=self.calibration.to_dict(),
                    unit=self.unit, levels={k: float(v) for k, v in self.levels.items()})

    @classmethod
    def from_dict(cls, d: dict) -> "CalSession":
        d = d or {}
        return cls(runs=list(d.get("runs", [])),
                   compounds=[CalCompound.from_dict(c) for c in d.get("compounds", [])],
                   pins={(c, r): (None if v is None else (float(v[0]), float(v[1])))
                         for c, r, v in d.get("pins", [])},
                   excluded={k: dict(v) for k, v in d.get("excluded", {}).items()},
                   overrides={k: dict(v) for k, v in d.get("overrides", {}).items()},
                   settings=CalBuilderSettings.from_dict(d.get("settings", {})),
                   calibration=CalibrationSettings.from_dict(d.get("calibration", {})),
                   unit=str(d.get("unit", "")),
                   levels={k: float(v) for k, v in d.get("levels", {}).items()})

    def copy(self) -> "CalSession":
        return CalSession.from_dict(json.loads(json.dumps(self.to_dict())))


# --------------------------------------------------------------- the table
def _key(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", s.lower())


def _run_columns(header: list[str], runs: list[str]) -> dict[str, str]:
    """Which amount column belongs to which run: the column named after the run
    (ignoring case and punctuation), otherwise a column named after one part of
    the run name (e.g. "L3" for CAL_L3_mon) that no other run shares."""
    out = {}
    cols = {_key(h): h for h in header}
    for run in runs:
        if _key(run) in cols:
            out[run] = cols[_key(run)]
    for run in runs:
        if run in out:
            continue
        tokens = {_key(t) for t in re.split(r"[_\-\s.]+", run) if t}
        hits = [h for k, h in cols.items() if k and k in tokens]
        hits = [h for h in hits if sum(_key(h) in {_key(t) for t in re.split(r"[_\-\s.]+", r)}
                                       for r in runs) == 1]
        if len(hits) == 1:
            out[run] = hits[0]
    return out


def _float(v) -> float | None:
    try:
        v = str(v).strip()
        return float(v) if v else None
    except ValueError:
        return None


def load_compounds(path, runs: list[str]) -> list[CalCompound]:
    """Read a compound table (CSV).

    Columns (only `compound`, t1 and t2 are needed): `compound`; `t1_s` (or
    `t1_expected_s`, or `t1_min` in minutes); `t2_s` (or `t2_expected_s`);
    `role` (`istd` marks an internal standard; a `target_table.csv` reads
    as it is); `quant_ion` (blank for the whole blob); `istd` (the internal
    standard to divide by, blank for the nearest); `amount` (the same amount
    in every run); then one column per run holding that run's amount, named
    after the run (or after the part of its name that differs, such as L3).
    """
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        header = [h for h in (reader.fieldnames or []) if h]
        rows = list(reader)
    low = {h.strip().lower(): h for h in header}

    def col(*names):
        return next((low[n] for n in names if n in low), None)

    c_name = col("compound", "name")
    c_t1 = col("t1_s", "t1_expected_s", "t1")
    c_t1m = col("t1_min", "t1_expected_min")
    c_t2 = col("t2_s", "t2_expected_s", "t2")
    if c_name is None or (c_t1 is None and c_t1m is None) or c_t2 is None:
        raise ValueError("the table needs a compound column, a t1 column (t1_s or t1_min) and a t2 column (t2_s)")
    c_role, c_is = col("role", "type"), col("is_istd", "internal_standard")
    c_ion, c_std, c_amount = col("quant_ion", "ion"), col("istd"), col("amount")
    known = {c for c in (c_name, c_t1, c_t1m, c_t2, c_role, c_is, c_ion, c_std, c_amount) if c}
    run_cols = _run_columns([h for h in header if h not in known], runs)
    out = []
    for row in rows:
        name = (row.get(c_name) or "").strip()
        if not name:
            continue
        t1 = _float(row.get(c_t1)) if c_t1 else None
        if t1 is None and c_t1m:
            m = _float(row.get(c_t1m))
            t1 = None if m is None else 60.0 * m
        t2 = _float(row.get(c_t2))
        if t1 is None or t2 is None:
            raise ValueError(f"{name}: t1 and t2 are needed")
        is_istd = ((row.get(c_role) or "").strip().lower() in ISTD_WORDS if c_role else False) or \
                  ((row.get(c_is) or "").strip().lower() in TRUE_WORDS if c_is else False)
        ion = _float(row.get(c_ion)) if c_ion else None
        amounts = {}
        default = _float(row.get(c_amount)) if c_amount else None
        for run in runs:
            v = _float(row.get(run_cols[run])) if run in run_cols else None
            v = default if v is None else v
            if v is not None:
                amounts[run] = v
        out.append(CalCompound(name=name, t1=t1, t2=t2, is_istd=is_istd,
                               quant_ion=None if ion is None else int(round(ion)),
                               istd=(row.get(c_std) or "").strip() if c_std else "", amounts=amounts))
    return out


def write_compounds(path, compounds: list[CalCompound], runs: list[str]) -> Path:
    """Write the compound table in the format `load_compounds` reads."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["compound", "t1_s", "t2_s", "role", "quant_ion", "istd"] + list(runs))
        for c in compounds:
            w.writerow([c.name, f"{c.t1:.3f}", f"{c.t2:.4f}", "istd" if c.is_istd else "target",
                        "" if c.quant_ion is None else c.quant_ion, c.istd]
                       + ["" if r not in c.amounts else f"{c.amounts[r]:g}" for r in runs])
    return path


# ---------------------------------------------------------------- exports
def point_rows(session: CalSession, points: dict, fits: dict) -> list[dict]:
    """Every compound in every run, with its normalised response and back-calculation."""
    rows = []
    for c in session.compounds:
        fit = fits.get(c.name)
        for run in session.runs:
            p = points.get((c.name, run))
            if p is None:
                continue
            back = err = float("nan")
            if isinstance(fit, CalibrationFit) and np.isfinite(p.ratio):
                back = float(fit.amount(np.array([p.ratio]))[0])
                if np.isfinite(p.x) and p.x > 0:
                    err = 100.0 * (back - p.x) / p.x
            rows.append(dict(compound=c.name, role="istd" if c.is_istd else "target", run=run,
                             amount=p.amount, x=p.x, matched=p.blob is not None, by_hand=p.pinned,
                             t1_s=p.t1, t2_s=p.t2, signal=p.signal, istd=p.istd,
                             istd_signal=p.istd_signal, ratio=p.ratio, max_point=p.height,
                             saturated=p.saturated,
                             in_curve=(not c.is_istd and run not in session.excluded.get(c.name, {})
                                       and np.isfinite(p.x) and np.isfinite(p.ratio)),
                             excluded_because=session.excluded.get(c.name, {}).get(run, ""),
                             back_calculated=back, error_pct=err))
    return rows


def curve_rows(session: CalSession, fits: dict) -> list[dict]:
    rows = []
    for name, fit in fits.items():
        std = session.istd_for(name)
        base = dict(compound=name, istd=std.name if std is not None else "",
                    x="amount / IS amount" if session.settings.amount_ratio else "amount",
                    unit=session.unit)
        if isinstance(fit, CalibrationFit):
            rows.append(dict(base, model=fit.model, weights=fit.weights, through_zero=fit.zero,
                             intercept=fit.intercept, slope=fit.slope, curvature=fit.curvature,
                             r2=fit.r2, mandel_p=fit.mandel_p, intercept_p=fit.intercept_p,
                             worst_error_pct=float(np.max(np.abs(fit.rel_error))) if fit.n else float("nan"),
                             passes=fit.ok, lod=fit.lod, loq=fit.loq, n=fit.n,
                             excluded="; ".join(f"{k}: {v}" for k, v in fit.excluded.items()),
                             note=""))
        else:
            rows.append(dict(base, model="", weights="", through_zero="", intercept="", slope="",
                             curvature="", r2="", mandel_p="", intercept_p="", worst_error_pct="",
                             passes=False, lod="", loq="", n=0, excluded="", note=str(fit)))
    return rows


def normalised_blob_rows(session: CalSession, blobs_by_run: dict, resolved: dict) -> list[dict]:
    """Every blob in every run divided by its nearest internal standard in that run.

    Both are whole-blob volumes (every ion), since an unknown blob has no
    quantifier ion of its own. The nearest internal standard is judged from
    the internal standards' own positions in that run (their matched blobs),
    with the same scaled distance used for compounds."""
    s = session.settings
    names = {(run, i): c for (c, run), i in resolved.items() if i is not None}
    rows = []
    for run in session.runs:
        blobs = blobs_by_run.get(run)
        if blobs is None or len(blobs) == 0:
            continue
        stds = [(std, resolved.get((std.name, run))) for std in session.istds()]
        stds = [(std, i) for std, i in stds if i is not None]
        for i in range(len(blobs)):
            name, ratio, sig = names.get((run, i), ""), float("nan"), float("nan")
            if stds:
                d = [np.hypot((blobs.t1[j] - blobs.t1[i]) / s.t1_scale, (blobs.t2[j] - blobs.t2[i]) / s.t2_scale)
                     for _, j in stds]
                std, j = stds[int(np.argmin(d))]
                sig = float(blobs.volume[j])
                ratio = float(blobs.volume[i]) / sig if sig > 0 else float("nan")
            rows.append(dict(run=run, blob=i, t1_s=float(blobs.t1[i]), t2_s=float(blobs.t2[i]),
                             volume=float(blobs.volume[i]), snr=float(blobs.snr[i]), compound=name,
                             istd=std.name if stds else "", istd_signal=sig, ratio=ratio))
    return rows


# ------------------------------------------------------------- suggestions
@dataclass(frozen=True)
class SuggestSettings:
    """How compounds are found from the standards alone (`suggest_compounds`).

    The calibrant slope range is wide on purpose. On the August 2026 standards
    (0.25 to 25 ng, a 100-fold range) the standard compounds grew only 5 to 15-fold:
    the tubes carry a background of the same compounds at the low end and the
    detector saturates at the top, giving log-log slopes of 0.4 to 0.9."""
    max_blobs: int = 1500          # the largest blobs of the highest level considered
    t1_window_s: float = 30.0      # the same compound in another run: this close in t1...
    t2_window_s: float = 0.3       # ...and t2...
    min_cosine: float = 0.85       # ...with spectra at least this alike
    slope_min: float = 0.3         # a calibrant: rises with every level, log volume against log amount
    slope_max: float = 1.5         # ...at a slope in this range...
    min_r2: float = 0.85           # ...along a fairly straight line
    constant_slope: float = 0.15   # an internal-standard candidate: |slope| below this...
    constant_cv: float = 0.35      # ...and volumes within this coefficient of variation
    min_volume_fraction: float = 1e-3   # blobs smaller than this share of the largest are noise, not compounds

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "SuggestSettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


@dataclass
class Suggestion:
    kind: str                       # "calibrant" or "constant"
    t1: float                       # s, in the reference (highest-level) run
    t2: float
    base_peak: int                  # most intense ion of the reference blob
    blobs: dict                     # run -> blob index
    volumes: dict                   # run -> blob volume
    slope: float                    # of log volume against log amount
    r2: float
    cv: float                       # coefficient of variation of the volumes
    min_cosine: float               # worst spectral match to the reference blob
    name: str = ""                  # a library internal standard its spectrum matches, if any

    @property
    def label(self) -> str:
        return self.name or f"m/z {self.base_peak} @ {self.t1 / 60:.2f} min, {self.t2:.2f} s"


# Spectral signatures of the shipped internal standards (presets/istd_library.json):
# quantifier ion, qualifier ions and their ratios to it, measured on the five standards
# found in the September 2026 exhaust runs (benchmarks/real_data_audit/internal_standards.py).
ISTD_SIGNATURES = {
    "acenaphthene-d10": (164, (162, 160), (0.963, 0.492)),
    "phenanthrene-d10": (188, (184, 160), (0.186, 0.13)),
    "n-tetracosane-d50": (66, (82, 98), (0.723, 0.512)),
    "chrysene-d12": (240, (236, 120), (0.277, 0.15)),
    "perylene-d12": (264, (260, 132), (0.199, 0.187)),
}


def recognise_istd(spectrum: np.ndarray, mz: np.ndarray, tolerance: float = 0.5) -> str:
    """The shipped internal standard whose signature this spectrum carries, or "".

    The quantifier ion must be among the spectrum's five largest, and each qualifier
    ratio within `tolerance` (relative) of the measured one: loose enough for a blob
    that shares its footprint with something else, tight enough that an undeuterated
    compound (whose ions are 8 to 12 mass units lower) never qualifies."""
    spectrum = np.asarray(spectrum, dtype=float)
    top = set(np.asarray(mz)[np.argsort(-spectrum)[:5]].tolist())
    lo = int(mz[0])
    for name, (quant, quals, ratios) in ISTD_SIGNATURES.items():
        if quant not in top or not lo <= quant <= int(mz[-1]):
            continue
        q = spectrum[quant - lo]
        ok = all(lo <= ion <= int(mz[-1]) and abs(spectrum[ion - lo] / q - r) <= tolerance * r
                 for ion, r in zip(quals, ratios))
        if q > 0 and ok:
            return name
    return ""


def usable_blobs(blobs: BlobSet, period: float, solvent: tuple[float, float] | None = None,
                 max_width_fraction: float = 0.06) -> np.ndarray:
    """Blobs worth matching across runs: not part of a streak along t2 and not inside
    the solvent front, `(start, end)` in seconds.

    A streak is a solvent tail or an overloaded compound smeared up the second
    dimension. It shows as one blob wider than `max_width_fraction` of the period at
    half height, or as three or more blobs with the same base peak within one
    modulation of each other spread over more than a third of the period (on the
    August 2026 standards, toluene-like m/z 91 streaks 9 to 14 min after injection)."""
    n = len(blobs)
    keep = np.ones(n, dtype=bool)
    if blobs.width_t2.size == n:
        keep &= blobs.width_t2 < max_width_fraction * period
    if solvent is not None:
        keep &= ~((blobs.t1 >= solvent[0]) & (blobs.t1 <= solvent[1] + 2 * period))
    if n:
        unit = _unit(blobs.spectra)
        col = np.round((blobs.t1 - blobs.t1_start) / period).astype(int)
        for i in range(n):
            near = np.flatnonzero(np.abs(col - col[i]) <= 1)
            same = near[unit[near] @ unit[i] > 0.98]
            # The same spectrum in three or more blobs spread over more than a second up the
            # second dimension is a smear. A big peak whose tail breaks into pieces keeps its
            # main blob (half or more of the group's signal); a smear has no such blob.
            if (same.size >= 3 and np.ptp(blobs.t2[same]) > 1.0
                    and blobs.volume[i] < 0.5 * blobs.volume[same].sum()):
                keep[i] = False
    return keep


def blob_heights(run, fold, blobs: BlobSet) -> np.ndarray:
    """The largest single-ion value recorded inside each blob, from the raw run.

    A detector near its limit stops growing with the amount: on the August 2026
    standards the largest point of the calibrants rose from about 2.1e5 counts at
    5 ng to only 2.5e5 at 25 ng (none above 3e5), and their curves bent over at the
    top. The height lets such points be seen and left out. MS-aware components have no
    footprint to look inside, so they get NaN."""
    if isinstance(blobs, ComponentBlobs):
        return np.full(len(blobs), np.nan)
    frame = make_frame(fold, blobs.settings)
    r, c = frame.scan_pixels()
    ok = (c >= 0) & (c < blobs.labels.shape[1]) & (r >= 0) & (r < blobs.labels.shape[0])
    lab = np.zeros(run.n_scans, dtype=np.int64)
    lab[ok] = blobs.labels[r[ok], c[ok]]
    top = np.asarray(run.data.max(axis=1).todense()).ravel()
    out = np.zeros(len(blobs) + 1)
    np.maximum.at(out, lab, top)
    return out[1:]


def t2_offsets(blobs_by_run: dict, ref: str, usable: dict | None = None, t1_window_s: float = 30.0,
               min_cosine: float = 0.95, top: int = 400, min_pairs: int = 10,
               min_volume_fraction: float = 1e-3) -> dict:
    """How far each run's second dimension sits from the reference run's: {run: (offset s, pairs)}.

    The same compound in two runs of one method should have the same t2. The phase of
    each run is estimated on its own, and on the August 2026 standards it came out 0.73 s
    apart for one of four runs (the bleed method), so blobs no longer lined up. Here the
    largest blobs of each run are paired with the reference blob of the most similar
    spectrum (cosine at least `min_cosine`) within `t1_window_s`, whatever its t2; the
    offset is the most common t2 difference (reference minus run, the median of the
    pairs within 0.1 s of the mode of a 0.05 s histogram). Adding it to the run's t2,
    by moving the run's phase back by as much, lines the run up with the reference.
    Fewer than `min_pairs` pairs gives no offset (NaN)."""
    rb = blobs_by_run[ref]
    ru = _unit(rb.spectra)
    period = rb.period
    out = {}
    for run, b in blobs_by_run.items():
        if b is None or run == ref:
            out[run] = (0.0, len(rb)) if run == ref else (float("nan"), 0)
            continue
        ok = usable.get(run) if usable else None
        ok_ref = usable.get(ref) if usable else None
        u = _unit(b.spectra)
        # Only real peaks: single-ion noise blobs pair up at random t2 and can outvote them.
        floor = min_volume_fraction * float(b.volume.max()) if len(b) else 0.0
        ref_floor = min_volume_fraction * float(rb.volume.max()) if len(rb) else 0.0
        order = [i for i in np.argsort(-b.volume) if (ok is None or ok[i]) and b.volume[i] >= floor][:top]
        dts = []
        for i in order:
            cand = np.flatnonzero((np.abs(rb.t1 - b.t1[i]) < t1_window_s) & (rb.volume >= ref_floor) &
                                  (True if ok_ref is None else ok_ref))
            if cand.size == 0:
                continue
            c = ru[cand] @ u[i]
            j = int(np.argmax(c))
            if c[j] >= min_cosine:
                dts.append((rb.t2[cand[j]] - b.t2[i] + period / 2) % period - period / 2)
        if len(dts) < min_pairs:
            out[run] = (float("nan"), len(dts))
            continue
        dts = np.asarray(dts)
        hist, edges = np.histogram(dts, bins=np.arange(-period / 2, period / 2 + 0.05, 0.05))
        mode = edges[int(np.argmax(hist))] + 0.025
        near = dts[np.abs(dts - mode) < 0.1]
        out[run] = (float(np.median(near)), int(near.size))
    return out


def amount_from_name(name: str) -> float | None:
    """The amount written into a run name, as in "L6_5ng_on_tube" (5) or "L2_5uL_025ng"
    (0.25: a leading zero marks a decimal, as written on the vials). None if absent."""
    m = re.search(r"(?<![\d.])(\d+(?:\.\d+)?)\s*(?:ng|pg|ug|µg)", name, re.IGNORECASE)
    if m is None:
        return None
    text = m[1]
    if len(text) > 1 and text.startswith("0") and "." not in text:
        text = "0." + text[1:]
    return float(text)


def _wrapped_distance(t1, t2, t1_ref: float, t2_ref: float, w1: float, w2: float, period: float) -> np.ndarray:
    """Scaled distance from (t1_ref, t2_ref), the smaller of reading each blob as it is or
    across the wrap boundary (one modulation earlier with t2 one period longer, or the
    reverse). A peak with t2 near the boundary is recorded on either side of it from run
    to run, by a few hundredths of a second of phase."""
    dt1 = np.asarray(t1) - t1_ref
    dt2 = np.asarray(t2) - t2_ref
    return np.minimum.reduce([np.hypot(dt1 / w1, dt2 / w2),
                              np.hypot((dt1 + period) / w1, (dt2 - period) / w2),
                              np.hypot((dt1 - period) / w1, (dt2 + period) / w2)])


def _unit(spectra: np.ndarray) -> np.ndarray:
    s = np.asarray(spectra, dtype=float)
    n = np.linalg.norm(s, axis=-1, keepdims=True)
    return s / np.where(n > 0, n, 1.0)


def suggest_compounds(blobs_by_run: dict, amounts: dict,
                      settings: SuggestSettings = SuggestSettings(),
                      usable: dict | None = None) -> list[Suggestion]:
    """Compounds that behave like calibrants, or like internal standards, across a set of standards.

    The idea in plain terms: in a calibration series every standard compound grows
    with the level and an internal standard stays the same, while most other peaks
    (bleed, contaminants, noise) do neither in a tidy way. Each of the largest blobs
    of the highest level is looked for in every other run (nearest apex inside the
    search window whose spectrum matches it best, if alike enough; each blob is used
    once, the largest reference blobs claiming first). Its volumes are then compared
    with the amounts on a log-log line. A calibrant is found in at least three runs
    and all but one, grows with every level, and lies on a fairly straight log-log
    line of positive slope. A constant peak is found in every run with a slope near
    zero and little scatter: that is what an internal standard looks like (so is a
    steady contaminant, so constant peaks are only named and ticked when their
    spectrum carries a shipped internal standard's signature, `recognise_istd`).

    Runs must share one phase (see `t2_offsets`). `amounts` maps run name to its
    level amount; `usable` maps run name to a boolean mask of blobs worth matching
    (`usable_blobs`, and not siloxane-like, say).
    """
    s = settings
    runs = [r for r in amounts if blobs_by_run.get(r) is not None and np.isfinite(amounts[r]) and amounts[r] > 0]
    if len(runs) < 3:
        return []
    ref = max(runs, key=lambda r: amounts[r])
    rb = blobs_by_run[ref]
    usable = {r: (np.asarray(usable[r], dtype=bool) if usable and r in usable else
                  np.ones(len(blobs_by_run[r]), dtype=bool)) for r in runs}
    units = {r: _unit(blobs_by_run[r].spectra) for r in runs}
    used = {r: np.zeros(len(blobs_by_run[r]), dtype=bool) for r in runs}
    floor = s.min_volume_fraction * float(rb.volume[usable[ref]].max()) if usable[ref].any() else 0.0
    order = [i for i in np.argsort(-rb.volume) if usable[ref][i] and rb.volume[i] >= floor][:s.max_blobs]
    x_all = {r: np.log10(amounts[r]) for r in runs}
    out, named = [], set()
    for i in order:
        found, cos = {ref: int(i)}, [1.0]
        for r in runs:
            if r == ref:
                continue
            b = blobs_by_run[r]
            d = _wrapped_distance(b.t1, b.t2, rb.t1[i], rb.t2[i], s.t1_window_s, s.t2_window_s, rb.period)
            cand = np.flatnonzero((d <= 1.0) & usable[r] & ~used[r])
            if cand.size == 0:
                continue
            c = units[r][cand] @ units[ref][i]
            k = int(np.argmax(c))
            if c[k] >= s.min_cosine:
                found[r] = int(cand[k])
                cos.append(float(c[k]))
        vols = {r: float(blobs_by_run[r].volume[j]) for r, j in found.items()}
        x = np.array([x_all[r] for r in found])
        y = np.log10(np.array([max(vols[r], 1e-12) for r in found]))
        if len(found) < 3 or np.ptp(x) == 0:
            continue
        slope, intercept = np.polyfit(x, y, 1)
        resid = y - (intercept + slope * x)
        ss = float(np.sum((y - y.mean()) ** 2))
        r2 = 1.0 - float(np.sum(resid ** 2)) / ss if ss > 0 else 0.0
        v = np.array(list(vols.values()))
        cv = float(v.std(ddof=1) / v.mean()) if v.mean() > 0 else float("inf")
        by_amount = [vols[r] for r in sorted(found, key=lambda r: amounts[r])]
        rising = bool(np.all(np.diff(by_amount) > 0))
        name = recognise_istd(rb.spectra[i], rb.mz)
        if name in named:
            name = ""              # a smaller piece of a standard already found (largest blobs come first)
        kind = None
        if name and len(found) == len(runs):
            # A known internal standard is kept whatever its trend: on the August 2026
            # standards they lost about half their signal in the 25 ng run.
            kind = "constant"
        elif (len(found) >= max(3, len(runs) - 1) and rising and s.slope_min <= slope <= s.slope_max
                and r2 >= s.min_r2):
            kind = "calibrant"
        elif len(found) == len(runs) and abs(slope) < s.constant_slope and cv <= s.constant_cv:
            kind = "constant"
        if kind is None:
            continue
        for r, j in found.items():
            used[r][j] = True
        if kind == "constant" and name:
            named.add(name)
        out.append(Suggestion(kind=kind, t1=float(rb.t1[i]), t2=float(rb.t2[i]),
                              base_peak=int(rb.mz[int(np.argmax(rb.spectra[i]))]), blobs=found,
                              volumes=vols, slope=float(slope), r2=float(r2), cv=cv, min_cosine=min(cos),
                              name=name if kind == "constant" else ""))
    out.sort(key=lambda g: g.t1)
    return out


# ------------------------------------------------------- applying the curves
@dataclass
class SampleResult:
    run: str
    compound: str
    found: bool
    t1: float
    t2: float
    signal: float
    istd: str
    istd_signal: float
    ratio: float
    amount: float                   # in the standards' unit; NaN when it cannot be read off the curve
    se: float
    cosine: float                   # spectral match to the compound in the standards
    flags: list = field(default_factory=list)


@dataclass
class CompoundTemplate:
    """What a compound looks like in the standards: where, and with what spectrum."""
    t1: float
    t2: float
    spectrum: np.ndarray            # unit vector over the standards' m/z axis
    mz0: int


def compound_templates(session: CalSession, blobs_by_run: dict, resolved: dict | None = None) -> dict:
    """Each compound's mean position and mean (unit) spectrum over its matched blobs in the standards."""
    resolved = session.resolve(blobs_by_run) if resolved is None else resolved
    out = {}
    for c in session.compounds:
        pos, spec, mz0 = [], [], None
        for run in session.runs:
            i, b = resolved.get((c.name, run)), blobs_by_run.get(run)
            if i is None or b is None:
                continue
            pos.append((b.t1[i], b.t2[i]))
            spec.append(_unit(b.spectra[i]))
            mz0 = int(b.mz[0])
        if pos:
            p = np.array(pos)
            out[c.name] = CompoundTemplate(t1=float(p[:, 0].mean()), t2=float(p[:, 1].mean()),
                                           spectrum=_unit(np.mean(spec, axis=0)), mz0=mz0)
    return out


def _on_axis(spectrum: np.ndarray, mz0: int, mz: np.ndarray) -> np.ndarray:
    """A template spectrum moved onto another run's m/z axis (zero where it has no ion)."""
    out = np.zeros(mz.size)
    src = np.arange(spectrum.size) + mz0
    keep = (src >= mz[0]) & (src <= mz[-1])
    out[src[keep] - int(mz[0])] = spectrum[keep]
    return out


def find_in_sample(template: CompoundTemplate, blobs: BlobSet, t2_offset: float = 0.0,
                   settings: CalBuilderSettings = CalBuilderSettings(), min_cosine: float = 0.8):
    """(blob index, cosine) of a compound in a sample run, or (None, best cosine).

    A sample holds far more than a standard, so position alone is not enough: among
    the blobs inside the search window around the compound's position in the
    standards (moved by the run's `t2_offset` from the standards, see `t2_offsets`),
    the one whose spectrum is most like the compound's is taken, if at least
    `min_cosine` alike."""
    if blobs is None or len(blobs) == 0:
        return None, float("nan")
    d = _wrapped_distance(blobs.t1, blobs.t2, template.t1, template.t2 - t2_offset,
                          settings.t1_window_s, settings.t2_window_s, blobs.period)
    cand = np.flatnonzero(d <= 1.0)
    if cand.size == 0:
        return None, float("nan")
    ref = _on_axis(template.spectrum, template.mz0, blobs.mz)
    cos = _unit(blobs.spectra[cand]) @ ref
    k = int(np.argmax(cos))
    return (int(cand[k]) if cos[k] >= min_cosine else None), float(cos[k])


def quantify_samples(session: CalSession, fits: dict, templates: dict, samples: dict,
                     t2_offsets_by_run: dict | None = None, min_cosine: float = 0.8) -> list[SampleResult]:
    """Read every compound in every sample run off its curve.

    Each compound and its internal standard are found by `find_in_sample`; the ratio
    of their volumes (or quantifier-ion volumes) is turned into an amount by inverse
    prediction on the compound's curve, with its standard error. With "x = amount /
    IS amount", the result is multiplied by the internal standard's median amount in
    the standards, which assumes the samples were spiked like them. Flags say when a
    result cannot be trusted as it stands: not found, internal standard not found,
    no curve, a curve that fails its own check, or an amount outside the range of the
    standards that made the curve (extrapolated)."""
    out = []
    for run, blobs in samples.items():
        offset = (t2_offsets_by_run or {}).get(run, 0.0)
        offset = 0.0 if not np.isfinite(offset) else offset
        hits = {}
        for c in session.compounds:
            tpl = templates.get(c.name)
            hits[c.name] = find_in_sample(tpl, blobs, offset, session.settings, min_cosine) if tpl else (None, np.nan)
        for c in session.analytes():
            i, cos = hits[c.name]
            std = session.istd_for(c.name)
            j = hits[std.name][0] if std is not None else None
            flags = []
            sig = blob_signal(blobs, i, c.quant_ion) if i is not None else float("nan")
            std_sig = blob_signal(blobs, j, std.quant_ion) if j is not None else float("nan")
            if i is None:
                flags.append("not found")
            if std is None:
                flags.append("no internal standard")
            elif j is None:
                flags.append(f"{std.name} not found")
            ratio = sig / std_sig if (np.isfinite(std_sig) and std_sig > 0 and np.isfinite(sig)) else float("nan")
            fit = fits.get(c.name)
            amount = se = float("nan")
            if not isinstance(fit, CalibrationFit):
                flags.append("no curve")
            else:
                if not fit.ok:
                    flags.append("curve fails its check")
                if np.isfinite(ratio):
                    a, e = fit.amount_with_se(np.array([ratio]))
                    amount, se = float(a[0]), float(e[0])
                    if session.settings.amount_ratio and std is not None:
                        std_amounts = [v for v in std.amounts.values() if np.isfinite(v)]
                        k = float(np.median(std_amounts)) if std_amounts else float("nan")
                        amount, se = amount * k, se * k
                        lo, hi = fit.x.min() * k, fit.x.max() * k
                    else:
                        lo, hi = fit.x.min(), fit.x.max()
                    if not np.isfinite(amount):
                        flags.append("beyond the curve")
                    elif amount < lo:
                        flags.append("below the lowest standard")
                    elif amount > hi:
                        flags.append("above the highest standard")
            out.append(SampleResult(
                run=run, compound=c.name, found=i is not None,
                t1=float(blobs.t1[i]) if i is not None else float("nan"),
                t2=float(blobs.t2[i]) if i is not None else float("nan"),
                signal=sig, istd=std.name if std is not None else "", istd_signal=std_sig, ratio=ratio,
                amount=amount, se=se, cosine=cos, flags=flags))
    return out


def sample_rows(results: list[SampleResult], unit: str = "", volumes: dict | None = None) -> list[dict]:
    """One row per sample and compound, with a concentration where the run has a volume (m³)."""
    rows = []
    for r in results:
        vol = (volumes or {}).get(r.run)
        conc = r.amount / vol if (vol and np.isfinite(r.amount)) else float("nan")
        rows.append(dict(run=r.run, compound=r.compound, amount=r.amount, se=r.se, unit=unit,
                         volume_m3=vol if vol else "", concentration=conc,
                         concentration_unit=f"{unit} m-3" if unit else "per m3",
                         ratio=r.ratio, istd=r.istd, signal=r.signal, istd_signal=r.istd_signal,
                         t1_s=r.t1, t2_s=r.t2, spectral_match=r.cosine, flags="; ".join(r.flags)))
    return rows


def save_session(path, session: CalSession, extra: dict | None = None) -> Path:
    path = Path(path)
    path.write_text(json.dumps(dict(kind="pygcgc-calibration", version=1, session=session.to_dict(),
                                    **(extra or {})), indent=1), encoding="utf-8")
    return path


def load_session(path) -> tuple[CalSession, dict]:
    d = json.loads(Path(path).read_text(encoding="utf-8"))
    if d.get("kind") != "pygcgc-calibration":
        raise ValueError(f"{Path(path).name} is not a calibration session")
    return CalSession.from_dict(d.get("session", {})), d
