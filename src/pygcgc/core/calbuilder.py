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
from .peaks import BlobSet

ISTD_WORDS = {"istd", "is", "internal standard", "internal_standard", "surrogate"}
TRUE_WORDS = {"1", "true", "yes", "y", "x"}


@dataclass(frozen=True)
class CalBuilderSettings:
    t1_window_s: float = 20.0      # search half-width along t1 for an automatic match
    t2_window_s: float = 0.25      # ...and along t2
    t1_scale: float = 60.0         # s of t1 that count as one unit of "nearest IS" distance
    t2_scale: float = 0.5          # s of t2 that count as one unit
    amount_ratio: bool = False     # x = amount / IS amount, for an IS not spiked equally into every run

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
    d = np.hypot((blobs.t1 - t1) / settings.t1_window_s, (blobs.t2 - t2) / settings.t2_window_s)
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

    def points(self, blobs_by_run: dict, resolved: dict | None = None) -> dict:
        """(compound, run) -> CalPoint for every compound and run."""
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
                amount = float(c.amounts.get(run, nan))
                std_sig = raw[(std.name, run)][3] if std is not None else nan
                ratio = sig / std_sig if (np.isfinite(std_sig) and std_sig > 0) else nan
                x = amount
                if self.settings.amount_ratio and std is not None:
                    std_amount = float(std.amounts.get(run, nan))
                    x = amount / std_amount if std_amount > 0 else nan
                out[(c.name, run)] = CalPoint(
                    compound=c.name, run=run, amount=amount, blob=i, pinned=(c.name, run) in self.pins,
                    t1=t1, t2=t2, signal=sig, istd=std.name if std is not None else "",
                    istd_signal=std_sig, ratio=ratio, x=x)
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
                    unit=self.unit)

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
                   unit=str(d.get("unit", "")))

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
                             istd_signal=p.istd_signal, ratio=p.ratio,
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
