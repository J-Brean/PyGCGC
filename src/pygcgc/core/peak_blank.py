"""Taking one run's peaks away from another's: a blank subtracted peak by peak, not pixel by pixel.

The idea in plain terms
------------------------
Subtracting a blank pixel by pixel subtracts noise from noise: where the blank's baseline is messy, the
difference is messier than either run. A fitted peak, though, is a number summed over many pixels, and
the same contaminant peak (a siloxane, a plasticiser, a column-bleed hump's neighbour) sits in the blank
and in the sample at nearly the same place with nearly the same mass spectrum. So:

1. For every peak in the sample, look for its partner in the blank: within a small window in both retention
   times (after mapping the blank onto the sample through the alignment, if there is one) and with a similar
   full mass spectrum. Partners are paired one to one, best spectral match first (`core.compare.match_peaks`).
2. Take the blank peak's volume times a factor the user chooses away from the sample peak's volume, never
   going below zero. The factor is the amount of blank the sample carries: 1 if the sample saw exactly the
   blank's contamination, less if the sample was sampled for less time or diluted.
3. A sample peak with no partner is untouched: it is not in the blank.

Every sample peak also says how it stands against the blank, in positive terms: **above blank** when its
volume is at least `significance` times the blank's contribution (or the blank has none), otherwise **near
blank**, i.e. mostly what the blank would explain.

`suggest_factor` offers a number to start from (the median sample-to-blank volume ratio of the shared peaks)
but never applies it: whether the shared peaks are pure contamination is the user's judgement.

The same pairing is used to correct quantified amounts (`core.quant`): the amount of a target in the blank,
times the factor, is taken from its amount in the sample.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from .compare import MAX_DT1, MAX_DT2, MIN_SIMILARITY, match_peaks

ABOVE, NEAR = "above blank", "near blank"
UNTOUCHED, REDUCED, REMOVED = "untouched", "reduced", "removed"


@dataclass(frozen=True)
class PeakBlankSettings:
    blank_run: str = ""
    factor: float = 1.0
    max_dt1: float = MAX_DT1              # s: how far a blank peak may sit from the sample's in t1
    max_dt2: float = MAX_DT2              # s: and in t2
    min_similarity: float = MIN_SIMILARITY    # weighted cosine squared of the two mass spectra, as ⚖️ Compare uses
    significance: float = 3.0             # a peak this many times the blank's contribution is "above blank"

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "PeakBlankSettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


@dataclass
class PeakBlankRow:
    index: int                            # the sample peak's position in its table
    t1: float
    t2: float
    volume: float
    blank_index: int | None = None
    blank_volume: float = 0.0
    similarity: float = float("nan")
    dt1: float = float("nan")
    dt2: float = float("nan")
    removed: float = 0.0
    corrected: float = 0.0
    ratio: float = float("inf")           # sample volume over the blank's contribution (factor x blank volume)

    @property
    def match(self) -> str:
        if self.blank_index is None or self.removed <= 0:
            return UNTOUCHED
        return REMOVED if self.corrected <= 0 else REDUCED

    @property
    def fraction_removed(self) -> float:
        return self.removed / self.volume if self.volume > 0 else 0.0


def subtract_blank_peaks(sample, blank, settings: PeakBlankSettings, to_sample=None) -> list:
    """One `PeakBlankRow` per sample peak, in the sample's order.

    `sample` and `blank` are peak tables shaped like `BlobSet` (t1, t2, volume, mz, spectra).
    `to_sample` takes the blank's `(t1, t2)` onto the sample's; leave it out when the runs are not aligned
    (the windows then absorb ordinary drift)."""
    n = len(sample)
    rows = [PeakBlankRow(i, float(sample.t1[i]), float(sample.t2[i]), float(sample.volume[i]),
                         corrected=max(float(sample.volume[i]), 0.0)) for i in range(n)]
    if n == 0 or blank is None or len(blank) == 0:
        return rows
    identity = to_sample or (lambda t1, t2: (np.asarray(t1, dtype=float), np.asarray(t2, dtype=float)))
    pairs = match_peaks(sample, blank, identity, max_dt1=settings.max_dt1, max_dt2=settings.max_dt2,
                        min_similarity=settings.min_similarity)
    for pair in pairs:
        if pair.a is None or pair.b is None:
            continue
        row = rows[pair.a]
        vb = max(float(pair.volume_b), 0.0)
        row.blank_index, row.blank_volume = int(pair.b), vb
        row.similarity, row.dt1, row.dt2 = float(pair.similarity), float(pair.dt1), float(pair.dt2)
        contribution = settings.factor * vb
        row.removed = float(min(contribution, max(row.volume, 0.0)))
        row.corrected = float(max(row.volume, 0.0) - row.removed)
        row.ratio = row.volume / contribution if contribution > 0 else float("inf")
    return rows


def blank_status(row: PeakBlankRow, settings: PeakBlankSettings) -> str:
    """"above blank" or "near blank", in the positive sense described in the module docstring."""
    return ABOVE if row.ratio >= settings.significance else NEAR


def corrected_spectra(sample, rows: list) -> np.ndarray:
    """The sample peaks' spectra scaled to their corrected volumes (a peak that was reduced keeps its shape)."""
    spectra = np.asarray(sample.spectra, dtype=float).copy()
    for row in rows:
        if row.volume > 0:
            spectra[row.index] *= row.corrected / row.volume
        else:
            spectra[row.index] *= 0.0
    return spectra


def suggest_factor(rows: list) -> dict:
    """A starting factor: the median sample-to-blank volume ratio of the peaks that have a partner.

    Returns `{"median", "q1", "q3", "n"}` (NaN and 0 when nothing is shared). Offered, never applied: if the
    shared peaks include real sample compounds the ratios run high, and the lower quartile is the better guide."""
    ratios = np.array([r.volume / r.blank_volume for r in rows
                       if r.blank_index is not None and r.blank_volume > 0 and r.volume > 0])
    if ratios.size == 0:
        return dict(median=float("nan"), q1=float("nan"), q3=float("nan"), n=0)
    q1, med, q3 = np.percentile(ratios, [25, 50, 75])
    return dict(median=float(med), q1=float(q1), q3=float(q3), n=int(ratios.size))


def summarise(rows: list, settings: PeakBlankSettings) -> dict:
    """Counts and the share of the run's signal removed."""
    total = sum(max(r.volume, 0.0) for r in rows)
    removed = sum(r.removed for r in rows)
    return dict(peaks=len(rows), matched=sum(r.blank_index is not None for r in rows),
                removed=sum(r.match == REMOVED for r in rows), reduced=sum(r.match == REDUCED for r in rows),
                untouched=sum(r.match == UNTOUCHED for r in rows),
                above=sum(blank_status(r, settings) == ABOVE for r in rows),
                near=sum(blank_status(r, settings) == NEAR for r in rows),
                signal_removed=removed / total if total > 0 else 0.0)

FIELDS = ["run", "peak", "t1_min", "t2_s", "volume", "blank_volume", "removed", "corrected", "percent_removed",
          "spectrum_match", "match", "against_blank", "blank_run", "factor"]


def write_rows_csv(path, run: str, rows: list, settings: PeakBlankSettings):
    """One line per sample peak: its volume, the blank's partner, what was removed and what is left."""
    import csv
    from pathlib import Path
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for r in rows:
            writer.writerow(dict(
                run=run, peak=r.index + 1, t1_min=f"{r.t1 / 60:.3f}", t2_s=f"{r.t2:.3f}", volume=f"{r.volume:.6g}",
                blank_volume=f"{r.blank_volume:.6g}" if r.blank_index is not None else "", removed=f"{r.removed:.6g}",
                corrected=f"{r.corrected:.6g}", percent_removed=f"{100 * r.fraction_removed:.1f}",
                spectrum_match=f"{r.similarity:.3f}" if np.isfinite(r.similarity) else "", match=r.match,
                against_blank=blank_status(r, settings), blank_run=settings.blank_run, factor=f"{settings.factor:g}"))
    return path
