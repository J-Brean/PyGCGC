"""Reading the NIST mass spectral library from its SDF text files.

The idea in plain terms
------------------------
NIST ships its 2011 library as plain-text SDF files. Each record is a
chemical-structure block (which we skip), then labelled text fields such as
the name, formula, CAS number and the mass spectrum as `m/z intensity`
pairs, then a line of four dollar signs. Three files matter here:

- `mainlib.sdf` and `replib.sdf`: electron-ionisation spectra (the main and
  the replicate library). Many records also carry an *estimated* Kovats
  retention index and its estimation error.
- `nist_ri.sdf`: about 350,000 *experimental* retention indices from the
  literature, one record per measurement, each saying which kind of column
  (standard non-polar, semi-standard non-polar, standard polar), which index
  scale and which temperature programme it came from.

Reading 700 MB of text takes a minute or two, so the first read is written
to a cache of compact arrays (`save_compact` / `load_compact`) that reloads
in seconds and is rebuilt only when a source file changes.

A library is held as one sparse (entries, unit mass) matrix, not one Python
object per spectrum: a quarter of a million small objects would be slow and
large. `CompactLibrary.entry(i)` builds a `LibraryEntry` on demand, so the
existing display, promotion and export code keeps working.

The NIST files are licensed to the user. Nothing here copies, ships or
caches them anywhere but the user's own local application-data folder.
"""

from __future__ import annotations

import hashlib
import os
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from .identify import LibraryEntry

PARSER_VERSION = 3
_FIELD = re.compile(r"^>\s*<([^>]+)>")

EI_FIELDS = ("NAME", "FORMULA", "MW", "CASNO", "NISTNO", "MASS SPECTRAL PEAKS",
             "ESTIMATED KOVATS RI", "RI ESTIMATION ERROR")
RI_FIELDS = ("NAME", "FORMULA", "CASNO", "NISTNO", "KOVATS INDEX", "COLUMN CLASS", "DATA TYPE",
             "PROGRAM TYPE")


class Cancelled(Exception):
    """Raised by a caller's `cancel` callback to stop a long read."""


# ------------------------------------------------------------------ parsing
def iter_sdf_records(path, fields=None, progress=None, cancel=None):
    """Yield `(title, {FIELD: text})` for each record of an SDF file.

    `fields` limits which labelled fields are kept (names are upper-cased);
    leaving the rest out keeps memory small, because the synonym lists are
    long. Multi-line values are joined with newlines. The file is read as
    Latin-1 so a stray byte can never stop the read, and universal newlines
    mean CRLF files parse the same as LF ones. `progress(fraction)` is
    called every few thousand records; `cancel()` may raise to stop.
    """
    wanted = None if fields is None else {f.upper() for f in fields}
    size = max(os.path.getsize(path), 1)
    position = 0
    title: str | None = None
    record: dict[str, str] = {}
    current: str | None = None
    buffer: list[str] = []
    count = 0

    def close_field():
        nonlocal current, buffer
        if current is not None and (wanted is None or current in wanted):
            record[current] = "\n".join(buffer)
        current, buffer = None, []

    with open(path, encoding="latin-1", newline=None) as handle:
        for line in handle:
            position += len(line)
            line = line.rstrip("\n")
            if title is None:
                title = line.strip()
                continue
            if line.startswith("$$$$"):
                close_field()
                yield title, record
                title, record = None, {}
                count += 1
                if count % 4000 == 0:
                    if cancel is not None:
                        cancel()
                    if progress is not None:
                        progress(min(position / size, 1.0))
                continue
            if line.startswith(">"):
                match = _FIELD.match(line)
                if match:
                    close_field()
                    current = match.group(1).strip().upper()
                    continue
            if current is not None:
                if line.strip() == "":
                    close_field()
                else:
                    buffer.append(line)
    if title is not None and record:             # a last record with no terminator
        close_field()
        yield title, record
    if progress is not None:
        progress(1.0)


def parse_peaks(text: str) -> tuple[np.ndarray, np.ndarray]:
    """`(mz, intensity)` from a MASS SPECTRAL PEAKS field.

    Reads all the numbers and pairs them, so one pair per line and several
    pairs per line (space or semicolon separated) both work. Returns unit
    masses, rounded."""
    tokens = text.replace(";", " ").replace(",", " ").split()
    if len(tokens) < 2:
        return np.zeros(0, dtype=np.int32), np.zeros(0, dtype=np.float32)
    values = np.array(tokens[: len(tokens) // 2 * 2], dtype=np.float64)
    mz = np.rint(values[0::2]).astype(np.int32)
    intensity = values[1::2].astype(np.float32)
    keep = (mz > 0) & (intensity > 0)
    return mz[keep], intensity[keep]


def clean_cas(text: str) -> str:
    """Digits only, so `7782-39-0` and NIST's `7782390` are the same key."""
    return re.sub(r"\D", "", text or "")


def clean_name(text: str) -> str:
    """Lower case with punctuation and spaces collapsed, for name matching."""
    return re.sub(r"[^a-z0-9]+", " ", (text or "").lower()).strip()


def _float(text: str | None) -> float:
    try:
        return float((text or "").split()[0])
    except (ValueError, IndexError):
        return float("nan")


# ------------------------------------------------------------------ the library
@dataclass(eq=False)
class CompactLibrary:
    """Every entry's spectrum in one sparse matrix, with columns of metadata.

    `matrix[i, m]` is entry `i`'s intensity at unit mass `m` (NIST scale,
    0 to 9999). `ri_est` and `ri_est_err` are NIST's own estimate and its
    error where it gave one (NaN otherwise). `ri_exp`, `ri_exp_sigma` and
    `ri_exp_n` are filled by `attach_experimental_ri`: the median of the
    literature indices for the chosen column class, a robust spread, and how
    many measurements gave it (0 when there is none).
    """
    names: list
    formula: list
    cas: list                   # digits only
    nistno: np.ndarray          # int64, -1 where absent
    mw: np.ndarray
    source: np.ndarray          # int8: 0 main library, 1 replicate library, 2 other
    ri_est: np.ndarray
    ri_est_err: np.ndarray
    matrix: sp.csr_matrix
    ri_exp: np.ndarray = field(default=None)
    ri_exp_sigma: np.ndarray = field(default=None)
    ri_exp_n: np.ndarray = field(default=None)
    ri_link: list = field(default=None)     # how the experimental RI was linked: "", "cas", "nistno" or "name"
    label: str = ""
    key: str = ""
    _cache: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        n = len(self.names)
        if self.ri_exp is None:
            self.ri_exp = np.full(n, np.nan, dtype=np.float32)
            self.ri_exp_sigma = np.full(n, np.nan, dtype=np.float32)
            self.ri_exp_n = np.zeros(n, dtype=np.int32)
            self.ri_link = [""] * n

    def __len__(self) -> int:
        return len(self.names)

    @property
    def mz_max(self) -> int:
        return self.matrix.shape[1] - 1

    def entry(self, i: int) -> LibraryEntry:
        row = self.matrix.getrow(int(i))
        ri = self.ri_exp[i] if np.isfinite(self.ri_exp[i]) else self.ri_est[i]
        return LibraryEntry(name=self.names[i], mz=row.indices.astype(np.int64),
                            intensity=row.data.astype(np.float64),
                            ri=None if not np.isfinite(ri) else float(ri),
                            formula=self.formula[i], cas=self.cas[i], is_decoy=False,
                            meta={"NISTNO": str(self.nistno[i])} if self.nistno[i] >= 0 else {})

    def compound_key(self, i: int) -> str:
        """What makes two entries the same compound: CAS number if there is one, otherwise the name."""
        return self.cas[i] if self.cas[i] else "n:" + clean_name(self.names[i])


def compact_from_entries(entries: list, label: str = "") -> CompactLibrary:
    """Pack a list of `LibraryEntry` (an MSP library) into a `CompactLibrary`.

    An entry's `ri` is treated as a measured index of unknown precision."""
    entries = [e for e in entries if not e.is_decoy]
    n = len(entries)
    mz_max = max((int(e.mz.max()) for e in entries if e.mz.size), default=0)
    indptr = np.zeros(n + 1, dtype=np.int64)
    indices, data = [], []
    for i, e in enumerate(entries):
        order = np.argsort(e.mz)
        indices.append(np.asarray(e.mz, dtype=np.int32)[order])
        data.append(np.asarray(e.intensity, dtype=np.float32)[order])
        indptr[i + 1] = indptr[i] + e.mz.size
    matrix = sp.csr_matrix(
        (np.concatenate(data) if data else np.zeros(0, np.float32),
         np.concatenate(indices) if indices else np.zeros(0, np.int32), indptr),
        shape=(n, mz_max + 1), dtype=np.float32)
    lib = CompactLibrary(
        names=[e.name for e in entries], formula=[e.formula for e in entries],
        cas=[clean_cas(e.cas) for e in entries], nistno=np.full(n, -1, dtype=np.int64),
        mw=np.full(n, np.nan, dtype=np.float32), source=np.full(n, 2, dtype=np.int8),
        ri_est=np.full(n, np.nan, dtype=np.float32), ri_est_err=np.full(n, np.nan, dtype=np.float32),
        matrix=matrix, label=label)
    for i, e in enumerate(entries):
        if e.ri is not None:
            lib.ri_exp[i], lib.ri_exp_sigma[i], lib.ri_exp_n[i], lib.ri_link[i] = e.ri, np.nan, 1, "library"
    return lib


def read_ei_library(paths, progress=None, cancel=None) -> CompactLibrary:
    """Parse `mainlib.sdf` (and `replib.sdf`, if given) into a `CompactLibrary`.

    `progress(fraction, text)` is called as the files are read."""
    paths = [Path(p) for p in paths]
    names, formula, cas, nistno, mw, source, ri_est, ri_err = [], [], [], [], [], [], [], []
    ptr, idx_parts, dat_parts = [0], [], []
    for number, path in enumerate(paths):
        kind = 1 if "rep" in path.stem.lower() else 0

        def step(fraction, number=number, path=path):
            if progress is not None:
                progress((number + fraction) / len(paths), f"Reading {path.name}")

        for title, rec in iter_sdf_records(path, EI_FIELDS, step, cancel):
            mz, intensity = parse_peaks(rec.get("MASS SPECTRAL PEAKS", ""))
            if mz.size == 0:
                continue
            names.append(rec.get("NAME", title).strip() or title)
            formula.append(rec.get("FORMULA", "").strip())
            cas.append(clean_cas(rec.get("CASNO", "")))
            try:
                nistno.append(int(rec.get("NISTNO", "").split()[0]))
            except (ValueError, IndexError):
                nistno.append(-1)
            mw.append(_float(rec.get("MW")))
            source.append(kind)
            ri_est.append(_float(rec.get("ESTIMATED KOVATS RI")))
            ri_err.append(_float(rec.get("RI ESTIMATION ERROR")))
            idx_parts.append(mz)
            dat_parts.append(intensity)
            ptr.append(ptr[-1] + mz.size)
    if not names:
        raise ValueError("no spectra found in " + ", ".join(p.name for p in paths))
    mz_max = max(int(p.max()) for p in idx_parts)
    matrix = sp.csr_matrix((np.concatenate(dat_parts), np.concatenate(idx_parts), np.array(ptr, dtype=np.int64)),
                           shape=(len(names), mz_max + 1), dtype=np.float32)
    matrix.sum_duplicates()
    return CompactLibrary(
        names=names, formula=formula, cas=cas, nistno=np.array(nistno, dtype=np.int64),
        mw=np.array(mw, dtype=np.float32), source=np.array(source, dtype=np.int8),
        ri_est=np.array(ri_est, dtype=np.float32), ri_est_err=np.array(ri_err, dtype=np.float32),
        matrix=matrix, label=", ".join(p.name for p in paths))


# ------------------------------------------------------------------ experimental indices
@dataclass(eq=False)
class RIRecords:
    """The literature retention indices, one array element per measurement."""
    value: np.ndarray
    cas: list
    nistno: np.ndarray
    name: list
    column_class: list
    data_type: list
    program_type: list

    def __len__(self) -> int:
        return int(self.value.size)


def read_ri_records(path, progress=None, cancel=None) -> RIRecords:
    """Parse `nist_ri.sdf`. Records without a numeric index are dropped."""
    value, cas, nistno, name, column_class, data_type, program_type = [], [], [], [], [], [], []
    for title, rec in iter_sdf_records(path, RI_FIELDS, (lambda f: progress(f, "Reading " + Path(path).name))
                                       if progress else None, cancel):
        v = _float(rec.get("KOVATS INDEX"))
        if not np.isfinite(v):
            continue
        value.append(v)
        cas.append(clean_cas(rec.get("CASNO", "")))
        try:
            nistno.append(int(rec.get("NISTNO", "").split()[0]))
        except (ValueError, IndexError):
            nistno.append(-1)
        name.append(clean_name(rec.get("NAME", title)))
        column_class.append(rec.get("COLUMN CLASS", "").strip())
        data_type.append(rec.get("DATA TYPE", "").strip())
        program_type.append(rec.get("PROGRAM TYPE", "").strip())
    return RIRecords(np.array(value, dtype=np.float64), cas, np.array(nistno, dtype=np.int64), name,
                     column_class, data_type, program_type)


COLUMN_CLASSES = ("Standard non-polar", "Semi-standard non-polar", "Standard polar")
DEFAULT_COLUMN_CLASS = "Semi-standard non-polar"       # 5% phenyl; confirmed by James, 20 October 2026
DEFAULT_DATA_TYPES = ("Linear RI", "Normal alkane RI")
DEFAULT_PROGRAM_TYPES = ("Ramp", "Complex")
MIN_SIGMA = 5.0             # RI units: no experimental index is known better than this
SINGLE_SIGMA = 15.0         # RI units: one or two measurements give no usable spread
DEFINITION_SIGMA = 1.0      # RI units: an n-alkane's index is fixed by the definition of the scale
ALL_COLUMNS = "All column types"
DEFINITION_TYPE = "Normal alkane RI value specified by scale definition"

# NIST's estimated RI (mainlib `ESTIMATED KOVATS RI`) describes a standard non-polar column. Checked against
# the experimental values in nist_ri.sdf on 20 October 2026 (45,219 and 18,850 compounds with both): on standard
# non-polar the estimate minus the measurement has median 0 and a robust spread equal to NIST's stated error; on
# semi-standard non-polar (5 % phenyl) it reads 21 units low with a spread 1.25 x NIST's stated error; on a polar
# column it is useless (median -454). Estimates with a stated error of 120 or more were off by -319 (median): not used.
# column class -> (estimate minus measurement, spread as a multiple of NIST's stated error)
ESTIMATE_BIAS = {"Standard non-polar": (0.0, 1.0), "Semi-standard non-polar": (-21.0, 1.25)}
ESTIMATE_MAX_ERROR = 120.0


def attach_experimental_ri(lib: CompactLibrary, records: RIRecords,
                           column_class: str = DEFAULT_COLUMN_CLASS,
                           data_types=DEFAULT_DATA_TYPES, program_types=DEFAULT_PROGRAM_TYPES) -> dict:
    """Give each library entry the median literature RI for `column_class`.

    Linking tries the CAS number, then NIST's own number, then the
    normalised name. The spread is `1.4826 x MAD` of the measurements (at
    least `MIN_SIGMA`), or `SINGLE_SIGMA` when there are fewer than three.
    The n-alkanes are special: NIST records their index as fixed by the
    definition of the scale (100 x carbon number), for every column, so
    those entries take that value with a spread of `DEFINITION_SIGMA`.
    Returns how many entries were linked by each route, for the report."""
    keep = np.array([c == column_class and d in data_types and p in program_types
                     for c, d, p in zip(records.column_class, records.data_type, records.program_type)], dtype=bool)
    definition = np.array([c == ALL_COLUMNS and d.startswith(DEFINITION_TYPE)
                           for c, d in zip(records.column_class, records.data_type)], dtype=bool)
    # per route: key -> measurements, and key -> values fixed by definition
    tables = {"cas": ({}, {}), "nistno": ({}, {}), "name": ({}, {})}
    for k in np.flatnonzero(keep | definition):
        v = float(records.value[k])
        slot = 1 if definition[k] else 0
        if records.cas[k]:
            tables["cas"][slot].setdefault(records.cas[k], []).append(v)
        if records.nistno[k] >= 0:
            tables["nistno"][slot].setdefault(int(records.nistno[k]), []).append(v)
        tables["name"][slot].setdefault(records.name[k], []).append(v)
    n = len(lib)
    lib.ri_exp = np.full(n, np.nan, dtype=np.float32)
    lib.ri_exp_sigma = np.full(n, np.nan, dtype=np.float32)
    lib.ri_exp_n = np.zeros(n, dtype=np.int32)
    lib.ri_link = [""] * n
    counts = {"cas": 0, "nistno": 0, "name": 0, "none": 0}
    for i in range(n):
        keys = (("cas", lib.cas[i] or None), ("nistno", int(lib.nistno[i]) if lib.nistno[i] >= 0 else None),
                ("name", clean_name(lib.names[i])))
        found = None
        for how, key in keys:
            if key is None:
                continue
            measured, fixed = tables[how]
            if key in fixed or key in measured:
                found = (how, key)
                break
        if found is None:
            counts["none"] += 1
            continue
        how, key = found
        measured, fixed = tables[how]
        if key in fixed:
            arr = np.asarray(fixed[key], dtype=np.float64)
            median, sigma = float(np.median(arr)), DEFINITION_SIGMA
        else:
            arr = np.asarray(measured[key], dtype=np.float64)
            median = float(np.median(arr))
            sigma = (max(1.4826 * float(np.median(np.abs(arr - median))), MIN_SIGMA)
                     if arr.size >= 3 else SINGLE_SIGMA)
        lib.ri_exp[i], lib.ri_exp_sigma[i], lib.ri_exp_n[i], lib.ri_link[i] = median, sigma, arr.size, how
        counts[how] += 1
    lib._cache.clear()
    return counts

# ------------------------------------------------------------------ cache
def library_cache_dir() -> Path:
    """Local application data, never the (possibly synced) data folder.
    `PYGCGC_LIBRARY_CACHE` overrides it, as the tests do."""
    if os.environ.get("PYGCGC_LIBRARY_CACHE"):
        return Path(os.environ["PYGCGC_LIBRARY_CACHE"])
    base = os.environ.get("LOCALAPPDATA") or str(Path.home() / ".cache")
    return Path(base) / "PyGCGC" / "libraries"


def source_key(paths) -> str:
    """A key that changes when a source file's path, size or modification time does."""
    h = hashlib.sha256(str(PARSER_VERSION).encode())
    for p in sorted(Path(p) for p in paths):
        st = p.stat()
        h.update(f"{p.resolve()}|{st.st_size}|{int(st.st_mtime)}".encode())
    return h.hexdigest()[:16]


def _pack(strings: list) -> np.ndarray:
    return np.frombuffer("\n".join(s.replace("\n", " ") for s in strings).encode("utf-8"), dtype=np.uint8)


def _unpack(array: np.ndarray, n: int) -> list:
    if n == 0:
        return []
    parts = array.tobytes().decode("utf-8").split("\n")
    return parts if len(parts) == n else (parts + [""] * n)[:n]


def save_compact(lib: CompactLibrary, directory) -> None:
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    np.savez(directory / "meta.npz", names=_pack(lib.names), formula=_pack(lib.formula), cas=_pack(lib.cas),
             nistno=lib.nistno, mw=lib.mw, source=lib.source, ri_est=lib.ri_est, ri_est_err=lib.ri_est_err,
             n=np.array(len(lib)), label=_pack([lib.label]), shape=np.array(lib.matrix.shape))
    m = lib.matrix
    np.savez(directory / "spectra.npz", data=m.data, indices=m.indices, indptr=m.indptr, shape=np.array(m.shape))


def load_compact(directory) -> CompactLibrary:
    directory = Path(directory)
    with np.load(directory / "meta.npz") as meta, np.load(directory / "spectra.npz") as spec:
        n = int(meta["n"])
        matrix = sp.csr_matrix((spec["data"], spec["indices"], spec["indptr"]),
                               shape=tuple(int(x) for x in spec["shape"]), dtype=np.float32)
        return CompactLibrary(
            names=_unpack(meta["names"], n), formula=_unpack(meta["formula"], n), cas=_unpack(meta["cas"], n),
            nistno=meta["nistno"], mw=meta["mw"], source=meta["source"], ri_est=meta["ri_est"],
            ri_est_err=meta["ri_est_err"], matrix=matrix, label=_unpack(meta["label"], 1)[0])


def open_ei_library(paths, progress=None, cancel=None, cache_root=None) -> CompactLibrary:
    """The EI library for these SDF files: from the local cache if the files
    are unchanged, otherwise parsed (a minute or two) and then cached."""
    paths = [Path(p) for p in paths]
    key = source_key(paths)
    folder = Path(cache_root or library_cache_dir()) / ("ei_" + key)
    if (folder / "meta.npz").exists() and (folder / "spectra.npz").exists():
        try:
            lib = load_compact(folder)
            lib.key = key
            return lib
        except Exception:                        # noqa: BLE001  (a damaged cache is just rebuilt)
            pass
    lib = read_ei_library(paths, progress, cancel)
    lib.key = key
    try:
        save_compact(lib, folder)
    except OSError:
        pass                                     # a full disk must not lose the library just read
    return lib


def find_nist_files(folder) -> dict:
    """Locate `mainlib`, `replib` and `nist_ri` SDF files in `folder` (or its `SDF` subfolder)."""
    folder = Path(folder)
    found = {}
    for base in (folder, folder / "SDF"):
        if not base.is_dir():
            continue
        for p in base.iterdir():
            stem = p.stem.lower()
            if p.suffix.lower() == ".sdf" and stem in ("mainlib", "replib", "nist_ri"):
                found[stem] = p
    return found
