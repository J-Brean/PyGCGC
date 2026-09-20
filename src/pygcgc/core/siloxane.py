"""Siloxane removal.

Ageing or overheated GC columns (and septa) bleed dimethylsiloxane
fragments, and cyclic siloxane contaminants (Dn: octamethylcyclotetra-,
decamethylcyclopenta-, dodecamethylcyclohexasiloxane and so on) turn up as
their own peaks. Both fragment under EI onto a shared, recognisable set of
ions, dominated by m/z 73 ((CH3)3Si+).

Flagging (`siloxane_score`, `flag_siloxane_blobs`)
---------------------------------------------------
Every TIC blob's summed spectrum (`core/peaks.BlobSet.spectra`) is scored:
the fraction of its volume sitting on the siloxane ions, times its cosine
similarity to the reference pattern. The default reference is the column
bleed measured in the supplied real acquisition (`BLEED_REFERENCE`); it can
be replaced with a preset, edited, or taken from any peak clicked in the 2D
view (`reference_from_spectrum`).
A blob needs at least `min_ions_present` of those ions to show any signal
at all before it's scored; otherwise the score is 0, not a spurious high
cosine from one or two ions. A blob at or above `score_threshold` is
flagged as probable siloxane. That threshold is provisional (a reasonable
guess, not yet chosen from a ROC curve on labelled data, because no real
or dummy data has been available to build one from); see HANDOVER.md.

Blob removal (`remove_siloxane_blobs`), the default method
-----------------------------------------------------------
Every blob is scored as above, and every blob at or above the threshold
(the "% certainty" slider in the 🫧 panel) is blocked out: each scan that
falls inside its footprint is set to zero, on every ion (or only on the
siloxane ions, if `block_all_ions` is off). Nothing is fitted. A blob is
either siloxane or it is not, and the threshold decides which. Because the
blobs are the ones the 🔍 detector already found, deleting a blob in the
🔍 panel keeps it out of the removal, and the outlines drawn in the 2D view
show exactly what the current threshold would block. The removal is a layer,
so it is undoable and nothing in the raw data changes.

Pixel removal (`remove_siloxanes`), the older method
-----------------------------------------------------
Independent of flagging, every pixel of every siloxane ion's image is fit
against the same reference pattern by non-negative least squares. With a
single reference (one generic siloxane profile, not one per homologue),
that reduces to a closed-form clipped projection:

    c = max(0, (observed . reference) / (reference . reference))

for every pixel, where `observed` is that pixel's value on each siloxane
ion and `c` is how much of the reference pattern best explains it. Only
`min(c * reference, observed)` is then subtracted from each ion, and
nothing is subtracted from a pixel with fewer than `min_ions_present`
siloxane ions present (there's nothing to fit a ratio to) or whose own
ion ratios don't actually resemble the reference (cosine similarity below
`min_fit_cosine`) -- a handful of ions rising together isn't enough on its
own; they need to rise together in roughly the right proportions. A real
analyte that happens to share one siloxane ion (behenic acid on m/z 73,
unit mass matching (CH3)3Si+) doesn't bring the *other* siloxane ions up
at all, let alone in the right ratio, so nothing is removed there; column
bleed and Dn peaks, which do show that ratio across many ions at once,
are subtracted almost entirely. Ions are never zeroed outright, only the
fitted siloxane component.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
from scipy import sparse

from .background import _scan_positions
from .fold import Fold
from .peaks import BlobSet, Frame
from .run import Run

# Default reference: the column bleed measured in the supplied real acquisition
# (`real_data/`, untreated 70 eV, m/z 38 to 300). It is the spectrum that rises
# with the oven temperature over the last 30 minutes of the run: each ion's
# 20 s median signal was regressed on m/z 207's, and the slope (207 = 100) is
# the reference intensity; only ions with R^2 >= 0.9 and slope >= 2 are kept,
# which drops the constant chemical background (m/z 39, 40, 44 and so on). It is a
# low-bleed 5%-phenyl-type pattern: m/z 207 dominates and m/z 73 is only about
# 7%, unlike the textbook cyclic Dn pattern the earlier default assumed.
# The instrument's mass range ends at 300, so ions above that are absent.
BLEED_REFERENCE: dict[int, float] = {
    73: 6.7, 75: 2.5, 96: 9.6, 97: 2.8, 119: 2.0, 133: 7.3, 147: 2.1, 177: 3.2, 191: 8.1,
    192: 2.6, 193: 6.0, 207: 100.0, 208: 17.2, 209: 10.2, 267: 2.3, 281: 18.1,
    282: 4.3, 283: 2.8,
}

# A discrete cyclic siloxane peak measured in the same run (29.3 min, t2 0.5 s):
# baseline-subtracted spectrum, m/z 73 = 100. Use it (or click a peak of your own)
# where a Dn contaminant, not column bleed, is what needs removing.
CYCLIC_REFERENCE: dict[int, float] = {
    73: 100.0, 147: 50.4, 207: 6.6, 221: 7.0, 267: 1.8, 281: 56.2, 282: 15.9, 283: 9.8,
}

# The pre-2026-09-20 generic pattern, kept for comparison. It was reasoned from the
# usual shape of a Dn spectrum, not measured, and does not resemble this
# instrument's column bleed (it puts m/z 73 first).
GENERIC_DN_REFERENCE: dict[int, float] = {
    73: 100.0, 147: 45.0, 207: 20.0, 221: 8.0, 281: 25.0, 295: 8.0,
    327: 6.0, 341: 5.0, 355: 18.0, 401: 4.0, 415: 4.0, 429: 4.0,
}

PRESETS: dict[str, dict[int, float]] = {
    "Column bleed (measured)": BLEED_REFERENCE,
    "Cyclic siloxane peak (measured)": CYCLIC_REFERENCE,
    "Generic Dn (old default)": GENERIC_DN_REFERENCE,
}

SILOXANE_REFERENCE: dict[int, float] = BLEED_REFERENCE
SILOXANE_IONS: tuple[int, ...] = tuple(sorted(SILOXANE_REFERENCE))


def reference_from_spectrum(spectrum: np.ndarray, mz: np.ndarray, min_rel: float = 0.03,
                            max_ions: int = 20, min_mz: int = 45) -> dict[int, float]:
    """A siloxane reference (base ion = 100) from a mass spectrum, such as
    the one under a peak clicked in the 2D view. Keeps at most `max_ions`
    ions of at least `min_rel` of the base above `min_mz`; low masses are
    left out because water, nitrogen, oxygen and argon sit there in every
    spectrum. Edit the result afterwards: which ions belong is your call."""
    spec = np.clip(np.asarray(spectrum, dtype=np.float64), 0, None)
    spec = np.where(np.asarray(mz) >= min_mz, spec, 0.0)
    if spec.max() <= 0:
        return {}
    keep = np.flatnonzero(spec >= min_rel * spec.max())
    keep = keep[np.argsort(spec[keep])[::-1][:max_ions]]
    top = float(spec.max())
    return {int(mz[j]): round(float(100.0 * spec[j] / top), 1) for j in sorted(keep)}


@dataclass(frozen=True)
class SiloxaneSettings:
    ions: tuple[int, ...] = SILOXANE_IONS
    reference: tuple[float, ...] = tuple(SILOXANE_REFERENCE[m] for m in SILOXANE_IONS)
    score_threshold: float = 0.3     # provisional; see the module docstring
    min_ions_present: int = 3        # a blob or pixel needs at least this many siloxane ions with signal
    min_fit_cosine: float = 0.8      # ...and its own ion ratios need to resemble the reference this closely (pixel method only)
    method: str = "blobs"            # "blobs": block out every blob scoring >= score_threshold; "pixels": per-pixel fit
    block_all_ions: bool = True      # blobs method: zero every ion in the footprint, not only the siloxane ions

    @classmethod
    def from_reference(cls, reference: dict[int, float], **kw) -> "SiloxaneSettings":
        """Settings whose ions and reference intensities come from a dict
        {m/z: intensity}, keeping the other settings at their defaults
        unless overridden in `kw`."""
        ions = tuple(sorted(int(m) for m in reference))
        return cls(ions=ions, reference=tuple(float(reference[m]) for m in ions), **kw)

    def reference_dict(self) -> dict[int, float]:
        return {int(m): float(r) for m, r in zip(self.ions, self.reference)}

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "SiloxaneSettings":
        d = {k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__}
        if "ions" in d:
            d["ions"] = tuple(d["ions"])
        if "reference" in d:
            d["reference"] = tuple(d["reference"])
        return cls(**d)


# ------------------------------------------------------------------ flag
def siloxane_score(spectrum: np.ndarray, mz0: int, s: SiloxaneSettings = SiloxaneSettings()) -> float:
    """How much one blob's summed spectrum looks like siloxane: the
    fraction of its volume on the siloxane ions, times its cosine
    similarity to the reference pattern. 0 if fewer than `min_ions_present`
    of them carry any signal."""
    idx = np.asarray(s.ions, dtype=np.int64) - mz0
    in_range = (idx >= 0) & (idx < spectrum.size)
    vals = np.zeros(len(s.ions))
    vals[in_range] = np.clip(spectrum[idx[in_range]], 0, None)
    total = float(np.clip(spectrum, 0, None).sum())
    if total <= 0 or int((vals > 0).sum()) < s.min_ions_present:
        return 0.0
    frac = float(vals.sum() / total)
    ref = np.asarray(s.reference, dtype=np.float64)
    denom = float(np.linalg.norm(vals) * np.linalg.norm(ref))
    cosine = float(np.dot(vals, ref) / denom) if denom > 0 else 0.0
    return frac * cosine


def siloxane_scores(blobs: BlobSet, s: SiloxaneSettings = SiloxaneSettings()) -> np.ndarray:
    """`siloxane_score` for every blob in a `BlobSet`."""
    if len(blobs) == 0:
        return np.array([])
    mz0 = int(blobs.mz[0])
    return np.array([siloxane_score(blobs.spectra[i], mz0, s) for i in range(len(blobs))])


def flag_siloxane_blobs(blobs: BlobSet, s: SiloxaneSettings = SiloxaneSettings()) -> np.ndarray:
    """Boolean flag per blob: probable siloxane, at or above `score_threshold`."""
    return siloxane_scores(blobs, s) >= s.score_threshold


# --------------------------------------------------------------- removal
@dataclass
class SiloxaneResult:
    run: Run
    removed_by_ion: dict = field(default_factory=dict)   # ion -> counts removed
    total_before: float = 0.0
    settings: SiloxaneSettings = field(default_factory=SiloxaneSettings)
    ions_used: tuple = ()                                 # siloxane ions actually within this run's m/z range
    n_blobs: int = 0                                      # blobs method: how many blobs were blocked out

    @property
    def removed_fraction(self) -> float:
        return float(sum(self.removed_by_ion.values()) / self.total_before) if self.total_before else 0.0

    def summary(self) -> str:
        if not self.ions_used:
            return "No siloxane ions in range; nothing to remove"
        removed = sum(self.removed_by_ion.values())
        if self.settings.method == "blobs":
            if self.n_blobs == 0:
                return "No blob reached the siloxane threshold; nothing removed"
            return (f"Blocked out {self.n_blobs} siloxane blob{'s' if self.n_blobs != 1 else ''} "
                    f"({self.removed_fraction:.1%} of the total signal)")
        if removed <= 0:
            return f"Nothing removed ({len(self.ions_used)} siloxane ions checked)"
        top = sorted(self.removed_by_ion.items(), key=lambda kv: -kv[1])[:4]
        ions_txt = ", ".join(str(m) for m, v in top if v > 0)
        return f"Removed {self.removed_fraction:.1%} of the signal on the siloxane ions; mostly m/z {ions_txt}"


def remove_siloxanes(run: Run, fold: Fold, settings: SiloxaneSettings = SiloxaneSettings()) -> SiloxaneResult:
    """Fit and subtract the siloxane component pixel by pixel. See the
    module docstring. A run with fewer than `min_ions_present` siloxane
    ions in its m/z range is returned unchanged."""
    s = settings
    total_before = float(run.data.sum())
    pairs = [(m, r) for m, r in zip(s.ions, s.reference) if run.mz[0] <= m <= run.mz[-1]]
    if len(pairs) < s.min_ions_present:
        return SiloxaneResult(run=run, total_before=total_before, settings=s, ions_used=())
    present_ions = np.array([m for m, _ in pairs], dtype=np.int64)
    ref = np.array([r for _, r in pairs], dtype=np.float64)

    cube = fold.cube(run, mz=present_ions, dtype=np.float64)   # (rows, cols, K)
    n_present = (cube > 0).sum(axis=2)
    num = (cube * ref[None, None, :]).sum(axis=2)
    den = float((ref ** 2).sum())
    cube_norm = np.sqrt((cube ** 2).sum(axis=2))
    cosine = np.divide(num, cube_norm * np.linalg.norm(ref), out=np.zeros_like(num),
                       where=cube_norm > 0)
    fit_here = (n_present >= s.min_ions_present) & (cosine >= s.min_fit_cosine)
    c = np.where(fit_here, np.clip(num / den, 0, None), 0.0)
    sub = np.clip(np.minimum(c[:, :, None] * ref[None, None, :], cube), 0, None)

    col_of_scan, row_of_scan = _scan_positions(fold)
    csc = run.data.tocsc(copy=True)
    idx_map = run.mz_index(present_ions)
    removed_by_ion: dict = {}
    for k, j in enumerate(idx_map):
        a, b = csc.indptr[j], csc.indptr[j + 1]
        if a == b:
            continue
        rows = csc.indices[a:b]
        take = sub[row_of_scan[rows], col_of_scan[rows], k]
        take = np.clip(np.minimum(take, csc.data[a:b]), 0, None)
        csc.data[a:b] -= take.astype(csc.data.dtype)
        removed_by_ion[int(present_ions[k])] = float(take.sum())
    csc.eliminate_zeros()
    corrected = Run(times=run.times, mz=run.mz, data=csc.tocsr(), name=run.name,
                    source=run.source, meta=dict(run.meta, siloxane=s.to_dict()))
    return SiloxaneResult(run=corrected, removed_by_ion=removed_by_ion, total_before=total_before,
                          settings=s, ions_used=tuple(int(m) for m in present_ions))


def remove_siloxane_blobs(run: Run, fold: Fold, blobs: BlobSet,
                          settings: SiloxaneSettings = SiloxaneSettings()) -> SiloxaneResult:
    """Block out every blob scoring at or above `settings.score_threshold`.
    `blobs` must have been found on this run (or one with the same scans)
    with this fold, on data that still contains the siloxane; see the
    module docstring."""
    s = settings
    total_before = float(run.data.sum())
    ions_used = tuple(int(m) for m in s.ions if run.mz[0] <= m <= run.mz[-1])
    chosen = np.flatnonzero(siloxane_scores(blobs, s) >= s.score_threshold) if len(blobs) else np.array([], int)
    if chosen.size == 0:
        return SiloxaneResult(run=run, total_before=total_before, settings=s, ions_used=ions_used)

    # Each scan sits in exactly one pixel of the blob label image (the same mapping
    # that gave each blob its spectrum), so a blob's scans are the ones to zero.
    r, c = Frame(fold=fold, wrap_rows=blobs.wrap_rows).scan_pixels()
    ok = (c >= 0) & (c < blobs.labels.shape[1])
    scan_label = np.zeros(run.n_scans, dtype=np.int64)
    scan_label[ok] = blobs.labels[r[ok], c[ok]]
    blocked = np.isin(scan_label, chosen + 1)

    coo = run.data.tocoo()
    hit = blocked[coo.row]
    if not s.block_all_ions:
        hit &= np.isin(run.mz[coo.col], np.asarray(s.ions))
    removed = np.bincount(coo.col[hit], weights=coo.data[hit].astype(np.float64), minlength=run.n_mz)
    keep = ~hit
    data = sparse.csr_matrix((coo.data[keep], (coo.row[keep], coo.col[keep])), shape=run.data.shape)
    corrected = Run(times=run.times, mz=run.mz, data=data, name=run.name,
                    source=run.source, meta=dict(run.meta, siloxane=s.to_dict()))
    removed_by_ion = {int(run.mz[j]): float(removed[j]) for j in np.flatnonzero(removed)}
    return SiloxaneResult(run=corrected, removed_by_ion=removed_by_ion, total_before=total_before,
                          settings=s, ions_used=ions_used, n_blobs=int(chosen.size))
