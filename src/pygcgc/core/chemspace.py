"""Chemical-space coordinates of a compound from its formula.

The idea in plain terms
------------------------
A molecular formula says how oxidised (O/C), how saturated (H/C) and how big (carbon number) a compound is.
Plotting every named peak in those coordinates sorts a chromatogram by chemistry instead of by retention:
alkanes sit at H/C = 2 and O/C = 0, aromatic hydrocarbons and PAHs low in H/C, acids and alcohols up the O/C
axis. This module only turns a formula into numbers; the plots are in `gui/chemspace_window.py`.

Formulas are the simple kind a library gives (`C10H22`, `C6H18O3Si3`, `C14D10`). Deuterium counts as hydrogen.
Anything with brackets, charges or unknown symbols returns `None` rather than a guess.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

# Monoisotopic masses (u) of the elements a GC-EI library contains.
MASS = {"H": 1.00782503207, "C": 12.0, "N": 14.0030740048, "O": 15.99491461956, "F": 18.99840322, "Si": 27.9769265325,
        "P": 30.97376163, "S": 31.972071, "Cl": 34.96885268, "Br": 78.9183371, "I": 126.904473, "B": 11.0093054,
        "Se": 79.9165218, "Ge": 73.9211778, "Sn": 119.9021947, "As": 74.9215965}
HALOGENS = ("F", "Cl", "Br", "I")
_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")
CH2 = 14.01565                                      # the Kendrick scale: CH2 = 14.00000


def parse_formula(formula: str) -> dict | None:
    """`{"C": 10, "H": 22}` for `C10H22`; `None` if the text is not a plain formula."""
    text = (formula or "").strip()
    if not text or not re.fullmatch(r"(?:[A-Z][a-z]?\d*)+", text):
        return None
    counts: dict = {}
    for symbol, n in _TOKEN.findall(text):
        symbol = "H" if symbol == "D" else symbol
        if symbol not in MASS:
            return None
        counts[symbol] = counts.get(symbol, 0) + (int(n) if n else 1)
    return counts or None


@dataclass(frozen=True)
class Space:
    """The chemical-space coordinates of one formula."""
    formula: str
    c: int
    h: int
    o: int
    n: int
    s: int
    si: int
    hal: int
    hc: float
    oc: float
    dbe: float
    mass: float
    kmd: float                  # Kendrick mass defect, CH2 scale
    kendrick_mass: float
    family: str                 # "CH", "CHO", "CHN", "CHOS", ...


def space_of(formula: str) -> Space | None:
    """Coordinates for a formula, or `None` if it cannot be read or has no carbon."""
    counts = parse_formula(formula)
    if counts is None or counts.get("C", 0) <= 0:
        return None
    c, h = counts.get("C", 0), counts.get("H", 0)
    o, n, s, si = counts.get("O", 0), counts.get("N", 0), counts.get("S", 0), counts.get("Si", 0)
    p = counts.get("P", 0)
    hal = sum(counts.get(x, 0) for x in HALOGENS)
    mass = sum(MASS[k] * v for k, v in counts.items())
    dbe = (c + si) - (h + hal) / 2.0 + (n + p) / 2.0 + 1.0
    km = mass * 14.0 / CH2
    parts = ["C", "H"] + [sym for sym, count in (("O", o), ("N", n), ("S", s), ("Si", si), ("P", p)) if count]
    if hal:
        parts.append("X")
    other = [k for k in counts if k not in ("C", "H", "O", "N", "S", "Si", "P") and k not in HALOGENS]
    family = "".join(parts) + ("+" if other else "")
    return Space(formula=formula.strip(), c=c, h=h, o=o, n=n, s=s, si=si, hal=hal, hc=h / c, oc=o / c, dbe=dbe,
                 mass=mass, kmd=round(km) - km, kendrick_mass=km, family=family)
