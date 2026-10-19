"""Target tables and calibration levels.

Target table columns: compound, role (target or istd), cls, quant_ion,
qual_ions and qual_ratios (semicolon-separated), t1_expected_s,
t1_window_s, t2_expected_s, t2_window_s, istd, stock_factor, ri. Only
compound, quant_ion, t1_expected_s and t2_expected_s are required.

Rows with role 'istd' are the internal standards. They make the default
phase anchors and the reference peaks every signal is normalised to.
"""

from __future__ import annotations

import csv
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .phase import Anchor


def load_anchors(path: str | Path, roles: tuple[str, ...] = ("istd",),
                 default_window: float = 12.0) -> list[Anchor]:
    anchors = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            if row.get("role", "").strip().lower() not in roles:
                continue
            window = row.get("t1_window_s") or default_window
            anchors.append(Anchor(
                name=row["compound"].strip(),
                mz=int(float(row["quant_ion"])),
                t1=float(row["t1_expected_s"]),
                t2=float(row["t2_expected_s"]),
                t1_window=float(window),
            ))
    if not anchors:
        raise ValueError(f"no rows with role in {roles} found in {path}")
    return anchors


# ------------------------------------------------------------- full targets
@dataclass
class Target:
    """One row of the target table.

    t1 is the absolute first-dimension retention (s); t2 is the true
    second-dimension time (s, may exceed the modulation period).
    `istd` names the internal standard this target is normalised to; if it
    is blank, the nearest internal standard is used.
    """
    compound: str
    role: str = "target"                 # target or istd
    cls: str = ""
    quant_ion: int = 0
    qual_ions: tuple[int, ...] = ()
    qual_ratios: tuple[float, ...] = ()  # qualifier / quantifier, from a clean standard
    t1: float = 0.0
    t1_window: float = 12.0              # half-width (s)
    t2: float = 0.0
    t2_window: float = 0.25              # half-width (s)
    istd: str = ""
    stock_factor: float = 1.0
    ri: float | None = None

    @property
    def is_istd(self) -> bool:
        return self.role == "istd"


def _split(text, cast):
    text = (text or "").strip()
    return tuple(cast(x) for x in text.replace(",", ";").split(";") if x.strip()) if text else ()


def _mass(value) -> int:
    mass = float(value)
    if not math.isfinite(mass) or mass <= 0 or not mass.is_integer():
        raise ValueError(f"Ion masses must be positive whole numbers, got {value!r}")
    return int(mass)


def _reader(fh, required):
    reader = csv.DictReader(fh)
    fields = [name.strip() for name in (reader.fieldnames or [])]
    if len(fields) != len(set(fields)):
        raise ValueError("CSV contains duplicate column names")
    missing = set(required) - set(fields)
    if missing:
        raise ValueError("CSV is missing required columns: " + ", ".join(sorted(missing)))
    reader.fieldnames = fields
    return reader


def load_targets(path: str | Path) -> list[Target]:
    """Read a target table CSV (see README for the columns)."""
    out = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in _reader(fh, ("compound", "quant_ion", "t1_expected_s", "t2_expected_s")):
            name = (row.get("compound") or "").strip()
            if not name:
                continue
            quals = _split(row.get("qual_ions"), _mass)
            ratios = _split(row.get("qual_ratios"), float)
            if ratios and len(ratios) != len(quals):
                raise ValueError(f"{name}: {len(quals)} qualifier ions but {len(ratios)} ratios")
            ri = row.get("ri")
            out.append(Target(
                compound=name,
                role=(row.get("role") or "target").strip().lower(),
                cls=(row.get("cls") or "").strip(),
                quant_ion=_mass(row["quant_ion"]),
                qual_ions=quals,
                qual_ratios=ratios,
                t1=float(row["t1_expected_s"]),
                t1_window=float(row.get("t1_window_s") or 12.0),
                t2=float(row["t2_expected_s"]),
                t2_window=float(row.get("t2_window_s") or 0.25),
                istd=(row.get("istd") or "").strip(),
                stock_factor=float(row.get("stock_factor") or 1.0),
                ri=float(ri) if ri not in (None, "") else None,
            ))
    names = [t.compound for t in out]
    dupes = {n for n in names if names.count(n) > 1}
    if dupes:
        raise ValueError(f"duplicate compounds in target table: {', '.join(sorted(dupes))}")
    for t in out:
        if t.role not in ("target", "istd"):
            raise ValueError(f"{t.compound}: role must be target or istd")
        if not all(math.isfinite(v) for v in (t.t1, t.t2, t.t1_window, t.t2_window, t.stock_factor, *t.qual_ratios)):
            raise ValueError(f"{t.compound}: retention positions, windows, factors and ratios must be finite")
        if t.t1_window <= 0 or t.t2_window <= 0 or t.stock_factor <= 0 or any(v < 0 for v in t.qual_ratios):
            raise ValueError(f"{t.compound}: windows and stock factors must be positive; ratios cannot be negative")
        if t.ri is not None and (not math.isfinite(t.ri) or t.ri <= 0):
            raise ValueError(f"{t.compound}: RI must be finite and positive")
        if t.istd and not any(ref.compound == t.istd and ref.is_istd for ref in out):
            raise ValueError(f"{t.compound}: internal standard '{t.istd}' is not in the table")
    if not out:
        raise ValueError("Target table contains no named targets")
    return out


_TARGET_COLUMNS = ["compound", "role", "cls", "quant_ion", "qual_ions", "qual_ratios", "t1_expected_s",
                  "t1_window_s", "t2_expected_s", "t2_window_s", "istd", "stock_factor", "ri"]


def write_targets_csv(targets: list[Target], path: str | Path, append: bool = False) -> Path:
    """Write (or append to) a target table CSV in the format `load_targets`
    reads -- how a 🏷️ search hit promoted to a target actually reaches the
    rest of the app: written here, then the target table is reloaded."""
    path = Path(path)
    rows, columns = [], list(_TARGET_COLUMNS)
    if append and path.exists():
        existing = load_targets(path)
        names = [t.compound for t in existing] + [t.compound for t in targets]
        if len(names) != len(set(names)):
            raise ValueError("Cannot append duplicate compound names to a target table")
        with path.open(newline="", encoding="utf-8-sig") as fh:
            reader = _reader(fh, _TARGET_COLUMNS[:1])
            columns = list(reader.fieldnames)
            columns.extend(c for c in _TARGET_COLUMNS if c not in columns)
            rows = list(reader)
    for t in targets:
        values = [t.compound, t.role, t.cls, t.quant_ion, ";".join(map(str, t.qual_ions)),
                  ";".join(f"{r:g}" for r in t.qual_ratios), t.t1, t.t1_window, t.t2, t.t2_window,
                  t.istd, t.stock_factor, "" if t.ri is None else t.ri]
        rows.append(dict(zip(_TARGET_COLUMNS, values)))
    # A minimal or reordered input table needs its own header respected. Replace
    # only after the complete new table has been written successfully.
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", newline="", encoding="utf-8",
                                         dir=path.parent, suffix=".csv", delete=False) as fh:
            temporary = Path(fh.name)
            writer = csv.DictWriter(fh, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return path


def load_levels(path: str | Path) -> dict[str, dict[str, float]]:
    """Amount of every compound at every calibration level.

    The CSV has a `level` column and one column per compound (amounts in ng,
    or whatever unit you use consistently). Returns {level: {compound: amount}}.
    """
    out: dict[str, dict[str, float]] = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in _reader(fh, ("level",)):
            level = (row.pop("level", "") or "").strip()
            if level:
                if level in out:
                    raise ValueError(f"Duplicate calibration level: {level}")
                amounts = {k.strip(): float(v) for k, v in row.items() if k and v not in ("", None) and v.strip()}
                if any(not math.isfinite(v) or v < 0 for v in amounts.values()):
                    raise ValueError(f"{level}: calibration amounts must be finite and non-negative")
                out[level] = amounts
    return out
