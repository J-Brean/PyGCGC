"""Assign each blob to the hydrocarbon it is.

The idea in plain terms
------------------------
Two things say what a hydrocarbon blob is: where it sits against the n-alkane
ladder, and what its mass spectrum looks like.

Position. The ladder (`hydrocarbons.AlkaneSeries`) turns a blob's t1 into a
carbon-equivalent number (retention index / 100) and its t2 into an offset from
the alkane line. The ladder can bend and wobble; the shape chosen on the series
(straight between points, smooth curve, or polynomial) decides how that
conversion follows it, so a curved ladder gives the right carbon number
everywhere along it.

Spectrum. Each hydrocarbon family has ions that dominate its spectrum. The
share of a blob's total signal carried by a family's ions (`TYPE_RULES`) says
which family it is closest to.

A blob is an n-alkane, labelled "n-C17", only if all of these hold: it sits on a
whole carbon number of the ladder (within `assign_carbon_tol`), within
`assign_t2_tol_s` of the alkane line in t2, its spectrum is alkane-like, and its
57/71/85 ion ratios match the alkane pattern. A blob that is alkane-like but fails
the position or ratio test is a branched alkane at that carbon-equivalent number
("branched C17-eq"). Other families keep their carbon-equivalent number too.
Blobs that match no family stay "unassigned" rather than being forced into one.

The ion sets and minimum shares in `TYPE_RULES` are starting values, not
validated on real data. They are plain numbers in this file so they can be
edited.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .hydrocarbons import AlkaneSeries, HydrocarbonSettings, ratios_match
from .peaks import BlobSet
from .pixelspec import ratios_to

# family: (ions that characterise it, smallest share of the blob's signal those ions must carry)
TYPE_RULES: dict[str, tuple[tuple[int, ...], float]] = {
    "alkane": ((57, 71, 85, 99), 0.30),
    "cycloalkane": ((55, 69, 83, 97, 82, 96), 0.30),
    "alkylbenzene": ((91, 92, 105, 119, 133), 0.25),
    "alkylnaphthalene": ((128, 141, 142, 155, 156, 170), 0.25),
}
UNASSIGNED = "unassigned"


@dataclass
class HydrocarbonCall:
    blob: int
    kind: str                # "n-alkane", "branched alkane", "cycloalkane", "alkylbenzene", ... or "unassigned"
    carbon: int | None       # whole carbon number for an n-alkane, else the nearest to the carbon-equivalent
    carbon_eq: float         # position on the ladder, in carbons
    dt2: float               # t2 above the alkane line (s)
    share: float             # signal share of the family's ions (0 when unassigned)

    @property
    def label(self) -> str:
        if self.kind == UNASSIGNED:
            return UNASSIGNED
        if self.kind == "n-alkane":
            return f"n-C{self.carbon}"
        return f"{self.kind} C{self.carbon}-eq"


def assign_hydrocarbons(blobs: BlobSet, series: AlkaneSeries,
                        settings: HydrocarbonSettings = HydrocarbonSettings()) -> list[HydrocarbonCall]:
    """One `HydrocarbonCall` per blob, in blob order. Needs a ladder of at least two alkanes."""
    mz = np.asarray(blobs.mz)
    spectra = np.asarray(blobs.spectra, dtype=float)
    lo, hi = settings.carbon_range
    # Without the alkane quantifier ion in range the ratio test cannot be made, so it cannot pass.
    can_check = mz.size > 0 and mz[0] <= settings.quant_ion <= mz[-1]
    calls = []
    for i in range(len(blobs)):
        carbon_eq = series.carbon_of(float(blobs.t1[i]))
        dt2 = float(blobs.t2[i]) - series.t2_of(carbon_eq)
        spec = spectra[i]
        total = float(spec.sum())
        shares = {name: (float(spec[np.isin(mz, ions)].sum()) / total if total > 0 else 0.0)
                  for name, (ions, _) in TYPE_RULES.items()}
        # The family whose share clears its own minimum by the widest margin.
        fits = {n: shares[n] / TYPE_RULES[n][1] for n in shares if shares[n] >= TYPE_RULES[n][1]}
        nearest = int(round(carbon_eq))
        if not fits:
            calls.append(HydrocarbonCall(i, UNASSIGNED, None, carbon_eq, dt2, 0.0))
            continue
        family = max(fits, key=fits.get)
        if family == "alkane":
            on_ladder = (lo <= nearest <= hi and abs(carbon_eq - nearest) <= settings.assign_carbon_tol
                         and abs(dt2) <= settings.assign_t2_tol_s
                         and can_check
                         and ratios_match(ratios_to(spec, mz, settings.quant_ion, settings.qual_ions), settings))
            kind = "n-alkane" if on_ladder else "branched alkane"
        else:
            kind = family
        calls.append(HydrocarbonCall(i, kind, nearest, carbon_eq, dt2, shares[family]))
    return calls
