"""🏷️ MS identification: MSP libraries, search, and confidence.

The idea in plain terms
------------------------
1. **Retention index.** Already built in milestone 10: `AlkaneSeries.carbon_of(t1)`
   is the Van den Dool carbon-number-equivalent at a t1; multiplying by
   100 gives the conventional Kovats-style retention index. Nothing new
   needed here -- `retention_index` below is a one-line wrapper, kept for
   a readable name at the call site.
2. **Library** (`load_msp`/`write_msp`/`build_library`). A library is a
   list of `LibraryEntry` (name, unit-mass spectrum, retention index,
   optional metadata), imported from any MSP file -- NIST via `lib2nist`,
   MoNA, MassBank, or one built here from accepted IDs on real standards
   (`build_inhouse_entry`). `build_library` packs every entry's spectrum
   into one sparse (entries, unit mass) matrix and precomputes each
   entry's own top ions, so searching doesn't re-walk every spectrum from
   scratch each time.
3. **Search** (`search`). Two cheap prefilters first -- a retention-index
   window and at least one shared "top ion" -- keep a real library from
   needing a full spectral comparison against every entry for every
   query. What's left is scored three ways: a forward and a reverse
   weighted dot product (Stein and Scott, 1994 intensity/mass weighting;
   reverse ignores mass channels the library entry itself has nothing
   on, which matters when the query spectrum carries a co-eluting
   neighbour's peaks too) and a peak-overlap ratio, combined and then
   reduced by a retention-index penalty that grows with how far the
   query's own RI sits from the library entry's.
4. **Confidence** (`assign_qvalues`, `classify`). `assign_qvalues` is
   target-decoy competition, the standard way to turn a raw score into
   an estimated false discovery rate without needing to know the true
   answer: every query's single best hit competes against the *whole*
   library, real entries and `generate_decoys`' shuffled ones alike, and
   the decoy-to-target ratio above each score gives that score's q-value.
   `classify` gives each best hit a simplified confidence tier, loosely
   after Schymanski and colleagues (2014) -- only tiers 2a ("probable
   structure", library match with retention-index confirmation), 2b (a
   library match alone, missing or ambiguous RI evidence) and 3
   ("tentative candidate", a weaker match) are reachable from a
   unit-mass EI spectrum plus RI: level 1 needs a co-injected reference
   standard (see `promote_to_target`) and level 4 needs an exact-mass
   molecular formula, neither available from this kind of data.
5. **Promote to target** (`promote_to_target`). Turn a confident hit into
   a `Target` row `core/quant.py` can measure and calibrate against, at
   the position it was actually found, with a quantifier and qualifier
   ions read straight from the library entry's own top peaks.

None of this has been checked against a real library yet -- there's no
NIST licence in hand, and no dummy pack in this environment either
(`tests/test_identify.py` has one test marked `@pytest.mark.dummy`,
against the dummy pack's own `synthetic_library.msp`, ready for when
either turns up). Everything else here is tested against small,
hand-built synthetic libraries with a known right answer. See HANDOVER.md.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from .hydrocarbons import AlkaneSeries
from .peaks import BlobSet
from .targets import Target


def retention_index(series: AlkaneSeries, t1: float) -> float:
    """The conventional (Kovats-style) retention index at `t1`, from an
    already-found n-alkane ladder: 100 x the Van den Dool carbon-number
    equivalent."""
    return 100.0 * series.carbon_of(float(t1))


# ------------------------------------------------------------------ settings
@dataclass(frozen=True)
class IdentifySettings:
    ri_window: float = 50.0           # RI units either side of the query, for prefiltering and the penalty
    ri_penalty_full: float = 100.0    # RI difference at which the penalty term saturates at ri_penalty_weight
    top_ions: int = 5                 # ions used for the cheap "shared top ion" prefilter
    min_shared_top_ions: int = 1
    weight_a: float = 0.6             # Stein and Scott (1994) intensity exponent
    weight_b: float = 3.0             # ... and m/z exponent
    forward_weight: float = 0.5
    reverse_weight: float = 0.3
    ratio_weight: float = 0.2
    ri_penalty_weight: float = 0.3
    min_score: float = 0.2            # hits below this aren't reported at all
    confident_score: float = 0.7      # tier 2a/2b needs at least this
    ambiguity: float = 0.02           # a runner-up within this fraction of the best score makes it ambiguous

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> "IdentifySettings":
        return cls(**{k: v for k, v in (d or {}).items() if k in cls.__dataclass_fields__})


# ------------------------------------------------------------------ library
@dataclass
class LibraryEntry:
    name: str
    mz: np.ndarray                    # unit masses actually stored (any order)
    intensity: np.ndarray             # matching intensities, same length as mz
    ri: float | None = None
    formula: str = ""
    cas: str = ""
    is_decoy: bool = False
    meta: dict = field(default_factory=dict)

    def top_ions(self, n: int) -> np.ndarray:
        if self.mz.size == 0:
            return np.array([], dtype=np.int64)
        order = np.argsort(self.intensity)[::-1]
        return self.mz[order[:n]].astype(np.int64)


@dataclass
class Library:
    """Every entry's spectrum in one sparse (entries, unit mass) matrix,
    plus each entry's own precomputed top ions -- built once by
    `build_library` so `search` never re-walks a spectrum from scratch."""
    entries: list
    mz_min: int
    mz_max: int
    matrix: sp.csr_matrix
    top_ions: list

    def __len__(self) -> int:
        return len(self.entries)


def build_library(entries: list, top_ions: int = 5) -> Library:
    entries = list(entries)
    with_peaks = [e for e in entries if e.mz.size]
    if not with_peaks:
        return Library(entries=entries, mz_min=0, mz_max=0, matrix=sp.csr_matrix((len(entries), 0)),
                       top_ions=[np.array([], dtype=np.int64) for _ in entries])
    mz_min = int(min(int(e.mz.min()) for e in with_peaks))
    mz_max = int(max(int(e.mz.max()) for e in with_peaks))
    n_mz = mz_max - mz_min + 1
    rows, cols, vals, tops = [], [], [], []
    for i, e in enumerate(entries):
        if e.mz.size:
            idx = e.mz.astype(np.int64) - mz_min
            rows.extend([i] * idx.size)
            cols.extend(idx.tolist())
            vals.extend(np.asarray(e.intensity, dtype=np.float64).tolist())
        tops.append(e.top_ions(top_ions))
    matrix = sp.csr_matrix((vals, (rows, cols)), shape=(len(entries), n_mz))
    return Library(entries=entries, mz_min=mz_min, mz_max=mz_max, matrix=matrix, top_ions=tops)


# ------------------------------------------------------------------ MSP import/export
_RI_KEYS = {"ri", "retentionindex", "retention_index", "rikovats", "kovatsri", "rt2ri"}
_SKIP_KEYS = {"name", "numpeaks", "num peaks", "formula", "cas", "cas#"}


def load_msp(path, is_decoy: bool = False) -> list:
    """Parse an MSP library file (the common NIST/MoNA/MassBank text
    interchange format: `Key: value` metadata lines, a `Num Peaks:` line,
    then that many `mz intensity` peak lines, blocks separated by a blank
    line or the next `Name:`). Unrecognised metadata keys are kept in
    `LibraryEntry.meta` rather than dropped. Peak lines may pack more
    than one `mz intensity` pair, semicolon-separated, as some exports do.

    Real exports vary more than this handles (quoted per-peak comments,
    unusual encodings); nothing here has been checked against a real
    NIST or MoNA export yet -- see HANDOVER.md.
    """
    entries: list = []
    name = ri = formula = cas = None
    meta: dict = {}
    peaks: list = []

    def flush():
        if name is None:
            return
        if peaks:
            mz = np.array([p[0] for p in peaks], dtype=np.int64)
            inten = np.array([p[1] for p in peaks], dtype=np.float64)
            order = np.argsort(mz)
            mz, inten = mz[order], inten[order]
        else:
            mz, inten = np.array([], dtype=np.int64), np.array([], dtype=np.float64)
        entries.append(LibraryEntry(name=name, mz=mz, intensity=inten, ri=ri, formula=formula or "",
                                    cas=cas or "", is_decoy=is_decoy, meta=dict(meta)))

    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n\r")
            if not line.strip():
                continue
            if ":" in line:
                key, _, value = line.partition(":")
                key_norm = key.strip().lower().replace(" ", "")
                value = value.strip()
                if key_norm == "name":
                    flush()
                    name, ri, formula, cas, meta, peaks = value, None, None, None, {}, []
                elif key_norm in _RI_KEYS:
                    try:
                        ri = float(value)
                    except ValueError:
                        pass
                elif key_norm == "formula":
                    formula = value
                elif key_norm in ("cas", "cas#"):
                    cas = value
                elif key_norm == "numpeaks":
                    pass
                else:
                    meta[key.strip()] = value
            else:
                for chunk in line.split(";"):
                    parts = chunk.split()
                    if len(parts) >= 2:
                        try:
                            peaks.append((float(parts[0]), float(parts[1])))
                        except ValueError:
                            continue
    flush()
    return entries


def write_msp(entries: list, path) -> Path:
    """Write library entries back out as MSP -- the way to persist an
    in-house library built from accepted IDs (`build_inhouse_entry`)."""
    path = Path(path)
    lines = []
    for e in entries:
        lines.append(f"Name: {e.name}")
        if e.ri is not None:
            lines.append(f"RI: {e.ri:g}")
        if e.formula:
            lines.append(f"Formula: {e.formula}")
        if e.cas:
            lines.append(f"CAS#: {e.cas}")
        for k, v in e.meta.items():
            lines.append(f"{k}: {v}")
        lines.append(f"Num Peaks: {e.mz.size}")
        for mz, inten in zip(e.mz, e.intensity):
            lines.append(f"{int(round(mz))} {inten:g}")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def build_inhouse_entry(name: str, spectrum: np.ndarray, mz0: int, ri: float | None = None,
                        formula: str = "", cas: str = "") -> LibraryEntry:
    """A `LibraryEntry` from one accepted standard's own apex spectrum
    (`BlobSet.apex_spectra[i]`, `mz0 = int(blobs.mz[0])`) -- the "in-house
    library built from accepted IDs" the spec asks for."""
    spectrum = np.asarray(spectrum, dtype=np.float64)
    nz = np.flatnonzero(spectrum > 0)
    return LibraryEntry(name=name, mz=(nz + mz0).astype(np.int64), intensity=spectrum[nz],
                        ri=ri, formula=formula, cas=cas, is_decoy=False)


def generate_decoys(entries: list, seed: int = 0) -> list:
    """One shuffled decoy per real entry: the same set of mass channels,
    with intensities randomly reassigned among them -- a plausible set of
    ions, an implausible pattern, the standard target-decoy strategy for
    estimating a library search's own false discovery rate without
    needing to already know the right answer. Kept at the same RI as the
    real entry, so it competes fairly within the same RI window."""
    rng = np.random.default_rng(seed)
    decoys = []
    for e in entries:
        if e.mz.size == 0:
            continue
        shuffled = rng.permutation(e.intensity)
        decoys.append(LibraryEntry(name=f"DECOY_{e.name}", mz=e.mz.copy(), intensity=shuffled,
                                   ri=e.ri, formula="", cas="", is_decoy=True))
    return decoys


# ------------------------------------------------------------------ scoring
def _weighted(mz_axis: np.ndarray, intensity: np.ndarray, a: float, b: float) -> np.ndarray:
    return np.power(np.clip(intensity, 0, None), a) * np.power(np.clip(mz_axis, 1, None), b)


def _cosine_sq(u: np.ndarray, v: np.ndarray) -> float:
    denom = float(np.sqrt(float((u ** 2).sum()) * float((v ** 2).sum())))
    if denom <= 0:
        return 0.0
    return float((u * v).sum() / denom) ** 2


def _score_one(query_mz: np.ndarray, query_intensity: np.ndarray, entry: LibraryEntry,
              settings: IdentifySettings) -> tuple:
    if entry.mz.size == 0:
        return 0.0, 0.0, 0.0
    lo = int(min(query_mz.min(), entry.mz.min()))
    hi = int(max(query_mz.max(), entry.mz.max()))
    axis = np.arange(lo, hi + 1)
    q = np.zeros(axis.size)
    q[query_mz.astype(np.int64) - lo] = query_intensity
    lib = np.zeros(axis.size)
    lib[entry.mz.astype(np.int64) - lo] = entry.intensity
    wq = _weighted(axis, q, settings.weight_a, settings.weight_b)
    wl = _weighted(axis, lib, settings.weight_a, settings.weight_b)
    forward = _cosine_sq(wq, wl)
    mask = lib > 0
    reverse = _cosine_sq(wq[mask], wl[mask]) if mask.any() else 0.0
    ratio = float((q[mask] > 0).sum() / mask.sum()) if mask.any() else 0.0
    return forward, reverse, ratio


@dataclass
class SearchHit:
    name: str
    score: float
    forward: float
    reverse: float
    ratio: float
    ri: float | None
    ri_diff: float | None
    is_decoy: bool
    entry_index: int
    q_value: float | None = None
    tier: str | None = None


def search(query_mz, query_intensity, library: Library, query_ri: float | None = None,
          settings: IdentifySettings = IdentifySettings()) -> list:
    """Every library entry that clears the prefilters, scored and sorted,
    best first. `query_ri` is optional -- without it (no alkane ladder
    found yet on this run), the RI prefilter and penalty are both skipped
    and everything rests on the spectrum alone."""
    query_mz = np.asarray(query_mz, dtype=np.int64)
    query_intensity = np.asarray(query_intensity, dtype=np.float64)
    if len(library) == 0 or query_mz.size == 0:
        return []
    top_order = np.argsort(query_intensity)[::-1]
    my_top = set(int(m) for m in query_mz[top_order[:settings.top_ions]])

    hits = []
    for i, entry in enumerate(library.entries):
        if query_ri is not None and entry.ri is not None and abs(entry.ri - query_ri) > settings.ri_window:
            continue
        shared = len(my_top & set(int(m) for m in library.top_ions[i]))
        if shared < settings.min_shared_top_ions:
            continue
        forward, reverse, ratio = _score_one(query_mz, query_intensity, entry, settings)
        ri_diff = abs(entry.ri - query_ri) if (entry.ri is not None and query_ri is not None) else None
        penalty = (min(1.0, ri_diff / settings.ri_penalty_full) * settings.ri_penalty_weight
                  if ri_diff is not None else 0.0)
        score = (settings.forward_weight * forward + settings.reverse_weight * reverse
                + settings.ratio_weight * ratio - penalty)
        if score < settings.min_score:
            continue
        hits.append(SearchHit(name=entry.name, score=float(score), forward=float(forward), reverse=float(reverse),
                              ratio=float(ratio), ri=entry.ri, ri_diff=ri_diff, is_decoy=entry.is_decoy,
                              entry_index=i))
    hits.sort(key=lambda h: h.score, reverse=True)
    return hits


# ------------------------------------------------------------------ confidence
def assign_qvalues(best_hits: list) -> None:
    """Target-decoy competition: rank every query's own best hit (`None`
    entries are skipped) by score, and at each rank set q = (decoys at or
    above this score) / (targets at or above), then take a running
    minimum from the best score down so q never *improves* going the
    wrong way. Mutates each hit's `q_value` in place."""
    real = [h for h in best_hits if h is not None]
    if not real:
        return
    order = sorted(range(len(real)), key=lambda i: real[i].score, reverse=True)
    n_target = n_decoy = 0
    raw = []
    for i in order:
        if real[i].is_decoy:
            n_decoy += 1
        else:
            n_target += 1
        raw.append(n_decoy / max(n_target, 1))
    running_min = 1.0
    qvals = [0.0] * len(order)
    for pos in range(len(order) - 1, -1, -1):
        running_min = min(running_min, raw[pos])
        qvals[pos] = running_min
    for pos, i in enumerate(order):
        real[i].q_value = qvals[pos]


def classify(hits: list, settings: IdentifySettings = IdentifySettings()) -> str | None:
    """A simplified confidence tier for the best of one query's hits
    (`hits[0]`), loosely after Schymanski and colleagues (2014) -- see
    the module docstring for which levels are and aren't reachable here."""
    if not hits:
        return None
    best = hits[0]
    ambiguous = len(hits) > 1 and hits[1].score >= best.score * (1.0 - settings.ambiguity)
    has_ri = best.ri_diff is not None
    if best.score >= settings.confident_score and has_ri and not ambiguous:
        return "2a"
    if best.score >= settings.confident_score:
        return "2b"
    return "3"


@dataclass
class Identification:
    blob_index: int
    hits: list

    @property
    def best(self):
        return self.hits[0] if self.hits else None


def identify_blobs(blobs: BlobSet, library: Library, series: AlkaneSeries | None = None,
                   settings: IdentifySettings = IdentifySettings()) -> list:
    """Search every blob's own apex spectrum (the cleanest one available
    -- footprint sums pick up a co-eluting neighbour) against a library,
    with the retention index from `series` if one's been found on this
    run. Assigns q-values and a confidence tier across the whole batch of
    best hits, the same run's blobs competing together the way
    target-decoy search is meant to."""
    mz0 = int(blobs.mz[0]) if blobs.mz.size else 0
    use_apex = blobs.apex_spectra.shape[0] == len(blobs) and blobs.apex_spectra.size
    out = []
    for i in range(len(blobs)):
        spectrum = blobs.apex_spectra[i] if use_apex else blobs.spectra[i]
        nz = np.flatnonzero(spectrum > 0)
        if nz.size == 0:
            out.append(Identification(blob_index=i, hits=[]))
            continue
        query_ri = None
        if series is not None and len(series) >= 2:
            try:
                query_ri = retention_index(series, float(blobs.t1[i]))
            except ValueError:
                query_ri = None
        hits = search(nz + mz0, spectrum[nz], library, query_ri, settings)
        out.append(Identification(blob_index=i, hits=hits))
    assign_qvalues([ident.best for ident in out])
    for ident in out:
        if ident.best is not None:
            ident.best.tier = classify(ident.hits, settings)
    return out


def identify_components(component_set, library: Library, series: AlkaneSeries | None = None,
                        settings: IdentifySettings = IdentifySettings()) -> list:
    """Search separated component spectra instead of mixed blob spectra."""
    out = []
    for index, component in enumerate(component_set.components):
        spectrum = np.asarray(component.spectrum, dtype=float)
        mz = np.asarray(component.mz, dtype=np.int64)
        keep = spectrum > 0
        if not np.any(keep):
            out.append(Identification(blob_index=index, hits=[]))
            continue
        query_ri = None
        if series is not None and len(series) >= 2:
            try:
                query_ri = retention_index(series, float(component.t1))
            except ValueError:
                query_ri = None
        out.append(Identification(
            blob_index=index,
            hits=search(mz[keep], spectrum[keep], library, query_ri, settings)))
    assign_qvalues([ident.best for ident in out])
    for ident in out:
        if ident.best is not None:
            ident.best.tier = classify(ident.hits, settings)
    return out


# ------------------------------------------------------------------ promote
def promote_to_target(entry: LibraryEntry, hit: SearchHit, t1: float, t2: float, n_qualifiers: int = 2,
                      role: str = "target", cls: str = "") -> Target:
    """Turn a confident hit into a `Target` row, at the position it was
    actually found on this run (not wherever the library entry's own
    reference run had it), with quantifier/qualifier ions read from the
    library entry's own strongest peaks."""
    order = np.argsort(entry.intensity)[::-1] if entry.mz.size else np.array([], dtype=np.int64)
    if order.size == 0:
        raise ValueError(f"{entry.name}: library entry has no peaks to promote")
    quant_ion = int(entry.mz[order[0]])
    quant_height = float(entry.intensity[order[0]])
    qual_idx = order[1:1 + n_qualifiers]
    qual_ions = tuple(int(m) for m in entry.mz[qual_idx])
    qual_ratios = tuple(float(entry.intensity[j] / quant_height) if quant_height > 0 else 0.0 for j in qual_idx)
    return Target(compound=entry.name, role=role, cls=cls, quant_ion=quant_ion, qual_ions=qual_ions,
                 qual_ratios=qual_ratios, t1=float(t1), t2=float(t2), ri=hit.ri if hit is not None else entry.ri)
