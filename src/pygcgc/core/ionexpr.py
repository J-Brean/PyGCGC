"""Ion expressions: what the image shows.

Typed in the ion box (M key) or built by clicking ions in the spectrum.

    TIC          total ion signal
    178          one ion
    178+202+228  sum of ions
    200:260      every ion from 200 to 260 inclusive
    57/71        ratio image (useful for telling compound classes apart)
    57+71/85     ratio of a sum to a single ion (the slash splits the two sides)

Ratio pixels where the denominator is below `ratio_floor` are set to NaN,
so noise doesn't produce huge ratios.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np

from .fold import Fold
from .run import Run

_TOKEN = re.compile(r"^\s*(\d+)\s*(?::\s*(\d+))?\s*$")


class IonExpressionError(ValueError):
    pass


@dataclass(frozen=True)
class IonExpr:
    numerator: tuple[int, ...] | None          # None means TIC
    denominator: tuple[int, ...] | None = None

    @property
    def is_tic(self) -> bool:
        return self.numerator is None and self.denominator is None

    @property
    def is_ratio(self) -> bool:
        return self.denominator is not None

    def text(self) -> str:
        def side(ions):
            if ions is None:
                return "TIC"
            # Write runs of three or more consecutive ions as a:b.
            parts, i = [], 0
            while i < len(ions):
                j = i
                while j + 1 < len(ions) and ions[j + 1] == ions[j] + 1:
                    j += 1
                if j - i >= 2:
                    parts.append(f"{ions[i]}:{ions[j]}")
                else:
                    parts.extend(str(x) for x in ions[i:j + 1])
                i = j + 1
            return "+".join(parts)
        s = side(self.numerator)
        return f"{s}/{side(self.denominator)}" if self.is_ratio else s

    def label(self) -> str:
        return "TIC" if self.is_tic else f"m/z {self.text()}"

    def ions(self) -> tuple[int, ...]:
        return tuple(self.numerator or ()) + tuple(self.denominator or ())

    def signal(self, run: Run, side: tuple[int, ...] | None) -> np.ndarray:
        return run.tic if side is None else run.xic(list(side))

    def image(self, run: Run, fold: Fold, ratio_floor: float = 20.0) -> np.ndarray:
        num = fold.image(self.signal(run, self.numerator))
        if not self.is_ratio:
            return num
        den = fold.image(self.signal(run, self.denominator))
        out = np.full_like(num, np.nan)
        ok = den > ratio_floor
        out[ok] = num[ok] / den[ok]
        return out

    def validate(self, run: Run) -> None:
        for m in self.ions():
            run.mz_index(m)   # raises KeyError with a helpful message


def _parse_side(text: str) -> tuple[int, ...] | None:
    text = text.strip()
    if text.upper() == "TIC":
        return None
    ions: list[int] = []
    for part in text.split("+"):
        m = _TOKEN.match(part)
        if not m:
            raise IonExpressionError(f"can't read '{part.strip()}': use numbers, '+', 'a:b' or TIC")
        a = int(m.group(1))
        b = int(m.group(2)) if m.group(2) else a
        if b < a:
            a, b = b, a
        ions.extend(range(a, b + 1))
    # Keep order of first appearance, drop repeats.
    seen: dict[int, None] = {}
    for i in ions:
        seen.setdefault(i, None)
    return tuple(seen)


def parse_ion_expr(text: str) -> IonExpr:
    if text is None or not text.strip():
        return IonExpr(None)
    parts = text.split("/")
    if len(parts) > 2:
        raise IonExpressionError("only one '/' is allowed")
    num = _parse_side(parts[0])
    den = _parse_side(parts[1]) if len(parts) == 2 else None
    if len(parts) == 2 and num is None and den is None:
        raise IonExpressionError("TIC/TIC is always 1")
    return IonExpr(num, den) if len(parts) == 2 else IonExpr(num)
