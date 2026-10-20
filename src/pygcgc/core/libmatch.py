"""Matching peaks to a spectral library using the spectrum and the retention index.

The idea in plain terms
------------------------
A peak is described by two things: its mass spectrum and where it came out
of the first column. A library compound is described the same way (a
spectrum and a retention index, RI). A good assignment needs both to agree.

- **Retention index.** The RI of a peak is 100 times the carbon number an
  n-alkane would have at the same first-dimension time (from the run's
  alkane ladder). Different columns give different RIs for the same
  compound, so the library RI is taken for one *column class* (here, by
  default, semi-standard non-polar, i.e. 5 % phenyl). Nothing else about the
  column, and not the oven programme, is needed: an RI is relative to the
  alkanes run on the same method, which cancels most of that. How far the
  peak is from the library value is expressed in standard deviations (`z`),
  combining three sources of uncertainty: where the peak's first-dimension
  time is known to (a fraction of a modulation), how well the library value
  is known (the spread of published values), and a fixed allowance for the
  difference between a GC x GC first dimension and the one-dimensional
  columns of the literature.
- **Spectrum.** The match factor (MF) is 999 times the squared cosine of
  the two intensity- and mass-weighted spectra (Stein and Scott, 1994). The
  reverse match factor (RMF) ignores ions the library entry does not have.
  These are dot-product match factors like NIST MS Search's, not identical to
  them, so NIST's thresholds are not quoted here as equivalent.
- **Combining them.** The RI is a gate and the spectrum ranks. Candidates
  whose RI is within `ri_gate_z` standard deviations (or that have no RI at
  all, which is flagged) are ranked by MF. A strong spectral match outside
  the gate is never silently dropped: it is shown as a demoted candidate and
  the peak is flagged `RI_CONFLICT`.

Every doubt is a named flag (`FLAGS`), and each peak gets one confidence
level: confident, probable, tentative or unknown. The thresholds are
provisional (see `MatchSettings`). False discovery rates come from
target-decoy competition: each library spectrum has a decoy with the same
ions and RI but shuffled intensities, and a peak whose best match is
no better than its best decoy is not believed.

Nothing here reaches Schymanski level 1: that needs a co-injected standard.
"""

from __future__ import annotations

import csv
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from .nist_sdf import (DEFAULT_COLUMN_CLASS, ESTIMATE_BIAS, ESTIMATE_MAX_ERROR, SINGLE_SIGMA, CompactLibrary)

FLAGS = {
    "RI_MISSING": "No alkane ladder for this peak, so the match rests on the spectrum alone.",
    "RI_ESTIMATED": "The library RI is NIST's estimate from structure, not a measurement.",
    "RI_LIBRARY_MISSING": "The library entry has no retention index, so it could not be checked.",
    "RI_EXTRAPOLATED": "The peak lies outside the carbon range of the alkane ladder.",
    "RI_BORROWED": "The alkane ladder was taken from another run.",
    "RI_TENSION": "The RI disagrees with the library by more than 2 standard deviations.",
    "RI_CONFLICT": "A good spectral match exists only outside the RI gate.",
    "ISOMER_AMBIGUOUS": "Another isomer (same formula) fits nearly as well.",
    "CLOSE_RUNNER_UP": "A different compound fits nearly as well.",
    "FEW_IONS": "Fewer than 5 ions above 1 % of the base peak.",
    "COELUTION": "The peak shares its fit with another, or is unresolved.",
    "HIGH_FDR": "The estimated false discovery rate is above the limit.",
    "CLASS_DISAGREES": "The compound class from the spectrum conflicts with the library formula.",
}
LEVELS = ("confident", "probable", "tentative", "unknown")
SCHYMANSKI = {"confident": "2a", "probable": "2b", "tentative": "3", "unknown": ""}


@dataclass(frozen=True)
class MatchSettings:
    """All thresholds are provisional and unvalidated; see docs/validation."""
    column_class: str = DEFAULT_COLUMN_CLASS
    weight_a: float = 0.6               # Stein and Scott (1994) intensity exponent
    weight_b: float = 3.0               # ... and m/z exponent
    ri_gate_z: float = 3.0
    ri_caution_z: float = 2.0
    sigma_systematic: float = 15.0      # RI units; GC x GC first dimension against 1D literature columns
    sigma_t1_fraction: float = 0.25     # of the modulation period, for fitted components
    mf_confident: float = 850.0
    mf_probable: float = 750.0
    mf_tentative: float = 650.0
    runner_up_mf: float = 30.0
    q_max: float = 0.05
    fdr_min_queries: int = 20           # fewer peaks than this cannot give a meaningful q-value
    few_ions: int = 5
    few_ion_fraction: float = 0.01
    use_estimated_ri: bool = True
    use_decoys: bool = True
    max_candidates: int = 5
    exclude_ions: tuple = ()
    decoy_seed: int = 0

    def to_dict(self) -> dict:
        d = asdict(self)
        d["exclude_ions"] = list(self.exclude_ions)
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "MatchSettings":
        d = {k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__}
        if "exclude_ions" in d:
            d["exclude_ions"] = tuple(d["exclude_ions"])
        return cls(**d)


# ------------------------------------------------------------------ the query
@dataclass
class Query:
    """One peak to identify: a full-ion spectrum and, if there is a ladder, an RI."""
    mz: np.ndarray
    intensity: np.ndarray
    ri: float = float("nan")
    ri_sigma: float = float("nan")
    flags: tuple = ()                   # RI_EXTRAPOLATED, RI_BORROWED, COELUTION, ... known before searching
    t1: float = float("nan")
    t2: float = float("nan")
    volume: float = 0.0
    label: str = ""                     # component id, for display and decisions
    class_label: str = ""               # the 🧬 compound class from the spectrum, if one was found with confidence


# The elements a compound of each 🧬 class must (or, for hydrocarbons, may only) contain. A library formula that breaks
# this contradicts the class, which usually means the library hit is a look-alike.
HYDROCARBON_CLASSES = {"alkane", "alkene", "cycloalkane", "alkylbenzene", "indane / indene / tetralin", "naphthalene", "PAH"}
CLASS_NEEDS = {"N-heteroaromatic": "N", "aliphatic amine": "N", "alcohol / ether": "O", "ketone / aldehyde": "O",
               "carboxylic acid": "O", "ester": "O", "aromatic oxygenate": "O", "S-heteroaromatic": "S",
               "siloxane": "Si", "phthalate": "O"}


def formula_elements(formula: str) -> set:
    """The element symbols in a formula such as `C10H22O` or `C6H18O3Si3`."""
    return set(re.findall(r"[A-Z][a-z]?", formula or ""))


def class_conflict(class_label: str, formula: str) -> bool:
    """True when the library formula cannot belong to the 🧬 class (a hydrocarbon class with a heteroatom, or
    a heteroatom class without that heteroatom). An empty class or formula never conflicts."""
    if not class_label or not formula:
        return False
    elements = formula_elements(formula)
    if class_label in HYDROCARBON_CLASSES:
        return bool(elements - {"C", "H"})
    need = CLASS_NEEDS.get(class_label)
    return bool(need) and need not in elements


def peak_ri(series, t1: float, sigma_t1: float) -> tuple[float, float, bool]:
    """`(RI, uncertainty, extrapolated)` for a first-dimension time, from an alkane ladder.

    The uncertainty is the time uncertainty times the local slope of RI
    against time. `extrapolated` is True outside the ladder's carbon range."""
    carbon = float(series.carbon_of(float(t1)))
    h = max(sigma_t1, 1e-3)
    slope = 100.0 * abs(float(series.carbon_of(float(t1) + h)) - float(series.carbon_of(float(t1) - h))) / (2 * h)
    carbons = series.carbons()
    extrapolated = carbon < float(min(carbons)) - 1e-9 or carbon > float(max(carbons)) + 1e-9
    return 100.0 * carbon, slope * sigma_t1, bool(extrapolated)


def queries_from_peaks(peaks, series=None, period: float = 8.0, settings: MatchSettings = MatchSettings(),
                       borrowed: bool = False, labels: list | None = None, t1_map=None,
                       class_labels: list | None = None) -> list:
    """`Query` objects for a `ComponentPeaks`-shaped object (t1, t2, volume, mz, spectra, co_elution).

    `series` is an `AlkaneSeries` or None. `borrowed` marks a ladder taken from another run; `t1_map`, if
    given, takes a peak's first-dimension time in this run to the time it would have in the ladder's run
    (through the alignment). `class_labels` gives each peak's 🧬 class, or "" for none."""
    out = []
    sigma_t1 = settings.sigma_t1_fraction * period
    has_ladder = series is not None and len(series) >= 2
    unresolved = getattr(peaks, "unresolved", None)
    ids = list(getattr(peaks, "component_id", []))
    for i in range(len(peaks)):
        flags = []
        ri = sigma = float("nan")
        if has_ladder:
            try:
                ladder_t1 = float(peaks.t1[i]) if t1_map is None else float(t1_map(float(peaks.t1[i]), float(peaks.t2[i])))
                ri, sigma, extrapolated = peak_ri(series, ladder_t1, sigma_t1)
                if extrapolated:
                    flags.append("RI_EXTRAPOLATED")
                if borrowed:
                    flags.append("RI_BORROWED")
            except ValueError:
                ri = sigma = float("nan")
        if not np.isfinite(ri):
            flags.append("RI_MISSING")
        if bool(peaks.co_elution[i]) or (unresolved is not None and len(unresolved) > i and bool(unresolved[i])):
            flags.append("COELUTION")
        out.append(Query(mz=np.asarray(peaks.mz, dtype=np.int64), intensity=np.asarray(peaks.spectra[i], dtype=np.float64),
                         ri=ri, ri_sigma=sigma, flags=tuple(flags), t1=float(peaks.t1[i]), t2=float(peaks.t2[i]),
                         volume=float(peaks.volume[i]),
                         label=labels[i] if labels is not None else (str(ids[i]) if i < len(ids) else str(i)),
                         class_label=class_labels[i] if class_labels is not None and i < len(class_labels) else ""))
    return out


# ------------------------------------------------------------------ results
@dataclass
class Candidate:
    index: int
    name: str
    formula: str
    cas: str
    mf: float                           # 0 to 999
    rmf: float
    ri_lib: float = float("nan")
    ri_sigma_lib: float = float("nan")
    ri_kind: str = "none"               # "experimental", "estimated", "library" (MSP file) or "none"
    n_ri: int = 0
    z: float = float("nan")
    in_gate: bool = True
    reason: str = ""


@dataclass
class PeakMatch:
    query: Query
    candidates: list = field(default_factory=list)      # best first: in-gate by MF, then demoted
    level: str = "unknown"
    flags: list = field(default_factory=list)
    q_value: float = float("nan")
    n_ions: int = 0
    decoy_mf: float = float("nan")

    @property
    def best(self):
        return self.candidates[0] if self.candidates else None

    @property
    def tier(self) -> str:
        return SCHYMANSKI[self.level]


# ------------------------------------------------------------------ library arrays
def library_ri(lib: CompactLibrary, settings: MatchSettings):
    """`(ri, sigma, kind)` per entry: the experimental index if there is one, else NIST's estimate.

    `kind` is 0 for none, 1 for experimental (or an MSP file's own value), 2 for estimated. NIST's estimate
    describes a standard non-polar column; it is corrected to the chosen column class using the bias and spread
    measured against the experimental values (`nist_sdf.ESTIMATE_BIAS`), and is not used for a polar column or
    where NIST's own stated error is 120 or more."""
    key = ("ri", settings.column_class, settings.use_estimated_ri)
    hit = lib._cache.get(key)
    if hit is not None:
        return hit
    ri = lib.ri_exp.astype(np.float64).copy()
    sigma = lib.ri_exp_sigma.astype(np.float64).copy()
    sigma[np.isnan(sigma) & np.isfinite(ri)] = SINGLE_SIGMA    # an MSP file's index has no stated spread
    kind = np.where(np.isfinite(ri), 1, 0).astype(np.int8)
    bias = ESTIMATE_BIAS.get(settings.column_class)
    if settings.use_estimated_ri and bias is not None:
        offset, factor = bias
        err = lib.ri_est_err.astype(np.float64)
        est = np.isnan(ri) & np.isfinite(lib.ri_est) & np.isfinite(err) & (err > 0) & (err < ESTIMATE_MAX_ERROR)
        ri[est] = lib.ri_est[est] - offset
        sigma[est] = err[est] * factor
        kind[est] = 2
    out = (ri, sigma, kind)
    lib._cache[key] = out
    return out


def _weight(data: np.ndarray, columns: np.ndarray, a: float, b: float, exclude: tuple) -> np.ndarray:
    """Stein and Scott weighting: intensity**a times m/z**b; excluded ions get zero."""
    w = np.power(data.astype(np.float64), a) * np.power(np.clip(columns, 1, None).astype(np.float64), b)
    if exclude:
        w[np.isin(columns, np.asarray(exclude))] = 0.0
    return w


def _weighted_library(lib: CompactLibrary, lo: int, hi: int, a: float, b: float, exclude: tuple, decoy: int | None):
    """The library's weighted spectra for m/z `lo`..`hi`, row-major, with squared row norms. Cached.

    Row-major, multiplied by a block of queries at once, is the fastest layout measured on the full library
    (3.6 ms per query for 32 at a time, against 47 ms one at a time). Only the latest range is kept per kind
    (genuine and decoy), because each copy is ~200 MB."""
    key = ("w", lo, hi, a, b, tuple(exclude), decoy)
    hit = lib._cache.get(key)
    if hit is not None:
        return hit
    for stale in [k for k in lib._cache if k[0] == "w" and k[-1] == decoy]:
        del lib._cache[stale]
    m = lib.matrix if decoy is None else _decoy_matrix(lib, decoy)
    hi_c = min(hi, m.shape[1] - 1)
    if lo <= hi_c:
        csr = m[:, lo:hi_c + 1].tocsr().astype(np.float32)
    else:
        csr = sp.csr_matrix((m.shape[0], 0), dtype=np.float32)
    w = _weight(csr.data, lo + csr.indices, a, b, exclude)
    csr.data = w.astype(np.float32)
    row_of = np.repeat(np.arange(csr.shape[0]), np.diff(csr.indptr))
    norm2 = np.bincount(row_of, weights=w * w, minlength=csr.shape[0]).astype(np.float64)
    out = (csr, norm2)
    lib._cache[key] = out
    return out


def _decoy_matrix(lib: CompactLibrary, seed: int) -> sp.csr_matrix:
    """Same ions in each entry, intensities shuffled among them (built once, not kept)."""
    m = lib.matrix
    rng = np.random.default_rng(seed)
    row_of = np.repeat(np.arange(m.shape[0]), np.diff(m.indptr))
    order = np.lexsort((rng.random(m.nnz), row_of))
    return sp.csr_matrix((m.data[order], m.indices.copy(), m.indptr.copy()), shape=m.shape, dtype=np.float32)


def _weighted_rows(lib: CompactLibrary, ids, lo: int, hi: int, a: float, b: float, exclude: tuple):
    """A few library rows, weighted the same way, as a row-major matrix with their squared norms."""
    rows = lib.matrix[np.asarray(ids, dtype=np.int64)]
    hi_c = min(hi, rows.shape[1] - 1)
    sub = rows[:, lo:hi_c + 1].tocsr().astype(np.float32) if lo <= hi_c else sp.csr_matrix((len(ids), 0), dtype=np.float32)
    w = _weight(sub.data, lo + sub.indices, a, b, exclude)
    sub.data = w.astype(np.float32)
    row_of = np.repeat(np.arange(sub.shape[0]), np.diff(sub.indptr))
    return sub, np.bincount(row_of, weights=w * w, minlength=sub.shape[0]).astype(np.float64)


# ------------------------------------------------------------------ scoring
def _forward_from_dots(dots: np.ndarray, norm2: np.ndarray, q2: float) -> np.ndarray:
    """Forward score (squared weighted cosine) against every library entry, from the dot products."""
    d = dots.astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where((norm2 > 0) & (q2 > 0), d * d / (norm2 * q2), 0.0)
    return np.clip(f, 0.0, 1.0)


def _reverse(sub: sp.csr_matrix, norm2: np.ndarray, qw: np.ndarray) -> np.ndarray:
    """Like the forward score, but the query is cut down to the ions the library entry has."""
    pattern = sub.copy()
    pattern.data = np.ones_like(pattern.data)
    q2_rev = np.asarray(pattern @ (qw.astype(np.float64) ** 2).astype(np.float32), dtype=np.float64).ravel()
    dots = np.asarray(sub @ qw.astype(np.float32), dtype=np.float64).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where((norm2 > 0) & (q2_rev > 0), dots ** 2 / (norm2 * q2_rev), 0.0)
    return np.clip(r, 0.0, 1.0)


def q_at_scores(win_scores, win_is_decoy, scores) -> np.ndarray:
    """False-discovery q-values, by target-decoy competition, at arbitrary scores.

    `win_*` describe each query's winner (its best target or its best decoy,
    whichever scored higher). The raw rate at a threshold is decoys over
    targets at or above it; the q-value is the lowest raw rate at any
    threshold at or below the score, so q never improves going the wrong way."""
    win_scores = np.asarray(win_scores, dtype=np.float64)
    win_is_decoy = np.asarray(win_is_decoy, dtype=bool)
    scores = np.atleast_1d(np.asarray(scores, dtype=np.float64))
    if win_scores.size == 0:
        return np.full(scores.shape, np.nan)
    order = np.argsort(win_scores)                      # ascending
    s_asc, d_asc = win_scores[order], win_is_decoy[order]
    n = s_asc.size
    # decoys and targets at or above each ascending position
    decoys_above = np.cumsum(d_asc[::-1])[::-1]
    targets_above = np.cumsum((~d_asc)[::-1])[::-1]
    raw = decoys_above / np.maximum(targets_above, 1)
    q_asc = np.minimum.accumulate(raw)                  # lowest raw rate at or below each score
    pos = np.searchsorted(s_asc, scores, side="right") - 1
    out = np.where(pos >= 0, q_asc[np.clip(pos, 0, n - 1)], 1.0)
    return np.clip(out, 0.0, 1.0)


def _collapse(lib: CompactLibrary, idx: np.ndarray, f: np.ndarray, keep: int) -> list:
    """Best-scoring entry per compound, best first, at most `keep` of them."""
    if f.size > 600:                    # replicates of one compound are few; no need to sort a quarter of a million
        top = np.argpartition(-f, 600)[:600]
        order = top[np.argsort(-f[top])]
    else:
        order = np.argsort(-f)
    seen, out = set(), []
    for k in order:
        i = int(idx[k])
        key = lib.compound_key(i)
        if key in seen:
            continue
        seen.add(key)
        out.append((i, float(f[k])))
        if len(out) >= keep:
            break
    return out


BLOCK = 32          # queries multiplied against the library together (see _weighted_library)


def match_queries(queries: list, lib: CompactLibrary, settings: MatchSettings = MatchSettings(),
                  mz_range: tuple | None = None, progress=None, cancel=None) -> list:
    """Match every query peak to the library; one `PeakMatch` each, in order.

    `mz_range` is the run's acquired `(lowest, highest)` unit mass; library
    ions outside it are dropped before comparing, so a library spectrum is
    not penalised for ions the instrument never recorded. It defaults to
    the queries' own m/z range."""
    n_q = len(queries)
    if n_q == 0 or len(lib) == 0:
        return [PeakMatch(query=q) for q in queries]
    if mz_range is None:
        mz_range = (int(min(q.mz.min() for q in queries if q.mz.size)), int(max(q.mz.max() for q in queries if q.mz.size)))
    lo, hi = int(mz_range[0]), int(mz_range[1])
    a, b = settings.weight_a, settings.weight_b
    exclude = tuple(settings.exclude_ions)
    sub, norm2 = _weighted_library(lib, lo, hi, a, b, exclude, None)
    dsub = dnorm2 = None
    if settings.use_decoys:
        dsub, dnorm2 = _weighted_library(lib, lo, hi, a, b, exclude, settings.decoy_seed)
    ri_lib, sig_lib, kind = library_ri(lib, settings)
    n_cols = sub.shape[1]
    has_ri = np.isfinite(ri_lib)
    best_target, best_decoy = np.zeros(n_q), np.zeros(n_q)

    def make(i, fwd, rev, gated, z):
        if kind[i] == 1 and lib.ri_link[i] == "library":
            ri_kind = "library"                       # an MSP file's own value
        else:
            ri_kind = ("none", "experimental", "estimated")[int(kind[i])]
        return Candidate(index=i, name=lib.names[i], formula=lib.formula[i], cas=lib.cas[i], mf=999.0 * fwd,
                         rmf=999.0 * rev, ri_lib=float(ri_lib[i]) if has_ri[i] else float("nan"),
                         ri_sigma_lib=float(sig_lib[i]) if has_ri[i] else float("nan"),
                         ri_kind=ri_kind, n_ri=int(lib.ri_exp_n[i]) if kind[i] == 1 else 0, z=z, in_gate=gated)

    # Prepare every query: its ions inside the acquired range, and its weighted vector.
    results, prepared = [], []
    for q in queries:
        mz = np.asarray(q.mz, dtype=np.int64)
        inten = np.asarray(q.intensity, dtype=np.float64)
        ok = (mz >= lo) & (mz <= hi) & (inten > 0)
        if exclude:
            ok &= ~np.isin(mz, np.asarray(exclude))
        mz, inten = mz[ok], inten[ok]
        match = PeakMatch(query=q, flags=list(q.flags))
        if inten.size:
            match.n_ions = int((inten > settings.few_ion_fraction * inten.max()).sum())
        results.append(match)
        qw_full = np.zeros(hi - lo + 1)
        if inten.size:
            qw_full[mz - lo] = _weight(inten, mz, a, b, ())
        q2 = float((qw_full ** 2).sum())       # the whole query counts, even ions beyond the library's heaviest
        prepared.append((qw_full, q2))

    for start in range(0, n_q, BLOCK):
        if cancel is not None:
            cancel()
        if progress is not None:
            progress(start / n_q, f"Matching peak {start + 1} of {n_q}")
        block = range(start, min(start + BLOCK, n_q))
        Q = np.zeros((n_cols, len(block)), dtype=np.float32)
        for k, qi in enumerate(block):
            Q[:, k] = prepared[qi][0][:n_cols]
        dots = np.asarray(sub @ Q)
        dots_d = np.asarray(dsub @ Q) if dsub is not None else None

        for k, qi in enumerate(block):
            q, match = queries[qi], results[qi]
            qw_full, q2 = prepared[qi]
            if q2 <= 0:
                continue
            qw = qw_full[:n_cols]
            f_all = _forward_from_dots(dots[:, k], norm2, q2)
            use_ri = bool(np.isfinite(q.ri))
            if use_ri:
                sigma_tot = np.sqrt(np.where(has_ri, sig_lib, 0.0) ** 2 + q.ri_sigma ** 2 + settings.sigma_systematic ** 2)
                with np.errstate(invalid="ignore"):
                    z_all = np.where(has_ri, (q.ri - ri_lib) / sigma_tot, np.nan)
                gate = ~has_ri | (np.abs(z_all) <= settings.ri_gate_z)
                idx = np.flatnonzero(gate)
            else:
                z_all = None
                idx = np.arange(len(lib))

            def z_of(i):
                return float(z_all[i]) if z_all is not None and has_ri[i] else float("nan")

            if idx.size:
                keep = _collapse(lib, idx, f_all[idx], settings.max_candidates)
                ids = [i for i, _ in keep]
                rows, nrm = _weighted_rows(lib, ids, lo, hi, a, b, exclude)
                rmf = _reverse(rows, nrm, qw)
                match.candidates = [make(i, fwd, rev, True, z_of(i)) for (i, fwd), rev in zip(keep, rmf)]
                best_target[qi] = keep[0][1] if keep else 0.0
                if dots_d is not None:
                    fd = _forward_from_dots(dots_d[idx, k], dnorm2[idx], q2)
                    best_decoy[qi] = float(fd.max())
                    match.decoy_mf = 999.0 * best_decoy[qi]

            # A strong match outside the RI gate is shown, demoted, when nothing inside it is convincing.
            if use_ri and 999.0 * best_target[qi] < settings.mf_probable:
                out_idx = np.flatnonzero(~gate)
                if out_idx.size and 999.0 * f_all[out_idx].max() >= settings.mf_probable:
                    outside = _collapse(lib, out_idx, f_all[out_idx], 3)
                    ids = [i for i, _ in outside]
                    rows, nrm = _weighted_rows(lib, ids, lo, hi, a, b, exclude)
                    for (i, fwd), rev in zip(outside, _reverse(rows, nrm, qw)):
                        c = make(i, fwd, rev, False, z_of(i))
                        c.reason = "outside the RI gate"
                        match.candidates.append(c)
                    match.flags.append("RI_CONFLICT")

    # False discovery rate by target-decoy competition across all the peaks.
    q_values = np.full(n_q, np.nan)
    if settings.use_decoys and n_q >= settings.fdr_min_queries:
        used = (best_target > 0) | (best_decoy > 0)
        q_values = q_at_scores(np.maximum(best_target, best_decoy)[used], (best_decoy > best_target)[used], best_target)
    for qi, match in enumerate(results):
        match.q_value = float(q_values[qi]) if match.candidates else float("nan")
        _grade(match, lib, settings)
    if progress is not None:
        progress(1.0, "Done")
    return results


def _grade(match: PeakMatch, lib: CompactLibrary, settings: MatchSettings) -> None:
    """Set the confidence level and the remaining flags of one peak."""
    flags = match.flags
    best = match.best
    order = list(FLAGS)

    def rank(flag):
        return order.index(flag) if flag in order else len(order)

    if best is None:
        match.level = "unknown"
        match.flags = sorted(set(flags), key=rank)
        return
    in_gate = [c for c in match.candidates if c.in_gate]
    if match.n_ions < settings.few_ions:
        flags.append("FEW_IONS")
    if best.ri_kind == "none" and np.isfinite(match.query.ri):
        flags.append("RI_LIBRARY_MISSING")
    if best.ri_kind == "estimated":
        flags.append("RI_ESTIMATED")
    if best.in_gate and np.isfinite(best.z) and abs(best.z) > settings.ri_caution_z:
        flags.append("RI_TENSION")
    runner = [c for c in in_gate[1:] if lib.compound_key(c.index) != lib.compound_key(best.index)
              and best.mf - c.mf <= settings.runner_up_mf]
    if runner:
        same_formula = any(c.formula and c.formula == best.formula for c in runner)
        flags.append("ISOMER_AMBIGUOUS" if same_formula else "CLOSE_RUNNER_UP")
    if class_conflict(match.query.class_label, best.formula):
        flags.append("CLASS_DISAGREES")
    high_fdr = np.isfinite(match.q_value) and match.q_value > settings.q_max
    if high_fdr:
        flags.append("HIGH_FDR")

    ri_evidence = np.isfinite(match.query.ri) and best.ri_kind != "none" and np.isfinite(best.z)
    mf = best.mf
    if not best.in_gate:
        level = "unknown"
    elif not ri_evidence:
        level = "probable" if mf >= settings.mf_confident else ("tentative" if mf >= settings.mf_tentative else "unknown")
    elif mf >= settings.mf_confident and best.ri_kind in ("experimental", "library") and abs(best.z) <= settings.ri_caution_z:
        level = "confident"
    elif mf >= settings.mf_probable and abs(best.z) <= settings.ri_gate_z:
        level = "probable"
    elif mf >= settings.mf_tentative:
        level = "tentative"
    else:
        level = "unknown"
    if level == "confident" and (runner or high_fdr or "FEW_IONS" in flags):
        level = "probable"
    match.level = level
    match.flags = sorted(set(flags), key=rank)


def describe_candidate_rank(match: PeakMatch, k: int) -> str:
    """One line saying why the k-th candidate ranks where it does."""
    c = match.candidates[k]
    parts = [f"MF {c.mf:.0f}", f"RMF {c.rmf:.0f}"]
    if np.isfinite(c.z):
        parts.append(f"RI {c.ri_lib:.0f} ({c.ri_kind}, z = {c.z:+.1f})")
    else:
        parts.append("no RI check")
    if not c.in_gate:
        parts.append("outside the RI gate")
    return "; ".join(parts)

# ------------------------------------------------------------------ summaries and export
def largest_peaks(matches: list, n: int = 10) -> list:
    """The `n` peaks with the most fitted signal, largest first, as `(match, share)` where `share` is the
    fraction of all the peaks' volume."""
    total = sum(max(m.query.volume, 0.0) for m in matches)
    ranked = sorted(matches, key=lambda m: m.query.volume, reverse=True)[:n]
    return [(m, (m.query.volume / total) if total > 0 else float("nan")) for m in ranked]


IDENTIFICATION_FIELDS = ["run", "rank", "peak", "t1_min", "t2_s", "volume", "percent_of_signal", "name", "formula", "cas",
                         "confidence", "schymanski_level", "mf", "rmf", "ri_measured", "ri_measured_sigma",
                         "ri_library", "ri_library_sigma", "ri_kind", "z", "q_value", "flags", "flag_text",
                         "runner_up", "decision", "decision_name", "decision_formula", "decision_cas", "library", "column_class"]


def write_identifications_csv(path, run: str, matches: list, decisions: dict | None = None, library: str = "",
                              settings: MatchSettings = MatchSettings()) -> Path:
    """One row per peak, largest first: the best candidate, every doubt as a code and as a sentence, and the
    user's decision. `decisions` maps a peak's label to `{"action": ..., "name": ...}`."""
    decisions = decisions or {}
    path = Path(path)
    total = sum(max(m.query.volume, 0.0) for m in matches)
    ranked = sorted(matches, key=lambda m: m.query.volume, reverse=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=IDENTIFICATION_FIELDS)
        writer.writeheader()
        for rank, m in enumerate(ranked, start=1):
            best = m.best
            runner = m.candidates[1] if len(m.candidates) > 1 and m.candidates[1].in_gate else None
            decision = decisions.get(m.query.label, {})

            def num(x, fmt):
                return "" if x is None or not np.isfinite(x) else format(x, fmt)
            writer.writerow(dict(
                run=run, rank=rank, peak=m.query.label, t1_min=num(m.query.t1 / 60.0, ".3f"), t2_s=num(m.query.t2, ".3f"),
                volume=num(m.query.volume, ".6g"), percent_of_signal=num(100.0 * m.query.volume / total if total > 0 else float("nan"), ".3f"),
                name=best.name if best else "", formula=best.formula if best else "", cas=best.cas if best else "",
                confidence=m.level, schymanski_level=m.tier,
                mf=num(best.mf if best else None, ".0f"), rmf=num(best.rmf if best else None, ".0f"),
                ri_measured=num(m.query.ri, ".1f"), ri_measured_sigma=num(m.query.ri_sigma, ".1f"),
                ri_library=num(best.ri_lib if best else None, ".1f"), ri_library_sigma=num(best.ri_sigma_lib if best else None, ".1f"),
                ri_kind=best.ri_kind if best else "", z=num(best.z if best else None, ".2f"), q_value=num(m.q_value, ".3f"),
                flags=";".join(m.flags), flag_text=" ".join(FLAGS.get(f, f) for f in m.flags),
                runner_up=runner.name if runner else "", decision=decision.get("action", ""),
                decision_name=decision.get("name", ""), decision_formula=decision.get("formula") or "",
                decision_cas=decision.get("cas") or "", library=library, column_class=settings.column_class))
    return path
