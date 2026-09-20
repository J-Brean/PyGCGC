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


def load_targets(path: str | Path) -> list[Target]:
    """Read a target table CSV (see README for the columns)."""
    out = []
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            name = (row.get("compound") or "").strip()
            if not name:
                continue
            quals = _split(row.get("qual_ions"), lambda x: int(float(x)))
            ratios = _split(row.get("qual_ratios"), float)
            if ratios and len(ratios) != len(quals):
                raise ValueError(f"{name}: {len(quals)} qualifier ions but {len(ratios)} ratios")
            ri = row.get("ri")
            out.append(Target(
                compound=name,
                role=(row.get("role") or "target").strip().lower(),
                cls=(row.get("cls") or "").strip(),
                quant_ion=int(float(row["quant_ion"])),
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
        if t.istd and t.istd not in names:
            raise ValueError(f"{t.compound}: internal standard '{t.istd}' is not in the table")
    return out


_TARGET_COLUMNS = ["compound", "role", "cls", "quant_ion", "qual_ions", "qual_ratios", "t1_expected_s",
                  "t1_window_s", "t2_expected_s", "t2_window_s", "istd", "stock_factor", "ri"]


def write_targets_csv(targets: list[Target], path: str | Path, append: bool = False) -> Path:
    """Write (or append to) a target table CSV in the format `load_targets`
    reads -- how a 🏷️ search hit promoted to a target actually reaches the
    rest of the app: written here, then the target table is reloaded."""
    path = Path(path)
    write_header = not (append and path.exists())
    mode = "a" if append and path.exists() else "w"
    with open(path, mode, newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        if write_header:
            w.writerow(_TARGET_COLUMNS)
        for t in targets:
            w.writerow([t.compound, t.role, t.cls, t.quant_ion, ";".join(map(str, t.qual_ions)),
                       ";".join(f"{r:g}" for r in t.qual_ratios), t.t1, t.t1_window, t.t2, t.t2_window,
                       t.istd, t.stock_factor, "" if t.ri is None else t.ri])
    return path


def load_levels(path: str | Path) -> dict[str, dict[str, float]]:
    """Amount of every compound at every calibration level.

    The CSV has a `level` column and one column per compound (amounts in ng,
    or whatever unit you use consistently). Returns {level: {compound: amount}}.
    """
    out: dict[str, dict[str, float]] = {}
    with open(path, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            level = (row.pop("level", "") or "").strip()
            if level:
                out[level] = {k.strip(): float(v) for k, v in row.items() if k and v not in ("", None)}
    return out
