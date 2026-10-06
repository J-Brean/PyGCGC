"""File import. Readers are looked up by file extension, so new formats
(mzML, vendor CSV) can be added later without touching anything else."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Iterable

from ..run import Run
from .andi import crop_run, read_andi, write_andi
from .markes_lsc import LSCFile, read_lsc

READERS: dict[str, Callable[[Path], Run]] = {
    ".cdf": read_andi,
    ".nc": read_andi,
    ".lsc": read_lsc,
}


def register_reader(extension: str, reader: Callable[[Path], Run]) -> None:
    READERS[extension.lower()] = reader


def read_run(path: str | Path, *, check_cancel=None, progress=None) -> Run:
    path = Path(path)
    reader = READERS.get(path.suffix.lower())
    if reader is None:
        known = ", ".join(sorted(READERS))
        raise ValueError(f"No reader for '{path.suffix}' files (known: {known}).")
    if check_cancel:
        check_cancel()
    if reader is read_lsc:
        return reader(path, check_cancel=check_cancel, progress=progress)
    return reader(path)


def discover_runs(folder, *, include_alternatives=False, check_cancel=None):
    """Discover files, preferring untreated 70 eV within unambiguous acquisitions.

    Explicit file opening bypasses this preference. Unknown/ambiguous identities
    remain visible; corrupt files remain visible so loading can report the error.
    """
    files = sorted(p for p in Path(folder).rglob('*') if p.is_file() and p.suffix.lower() in READERS)
    if include_alternatives:
        return files
    groups, selected = {}, []
    for path in files:
        if check_cancel:
            check_cancel()
        if path.suffix.lower() != '.lsc':
            selected.append(path)
            continue
        try:
            with LSCFile(path, check_cancel, validate_index=False) as source:
                identity = source.identity()
        except (ValueError, OSError) as exc:
            # A TOF-DS re-export ("Open TOF Data File", SD01 spectra) is a copy in a format
            # we don't read, saved beside the original; it is left out, not reported as broken.
            if 'no SD04 spectra' not in str(exc):
                selected.append(path)
            continue
        groups.setdefault(identity['acquisition_id'], []).append((path, identity))
    for group in groups.values():
        preferred = [p for p,m in group if m['energy_ev'] == 70 and m['processing_state'] == 'instrument']
        kinds = [(m['energy_ev'],m['processing_state']) for p,m in group]
        unambiguous = len(set(kinds)) == len(kinds) and all(energy is not None for energy,state in kinds)
        selected.extend(preferred if len(preferred) == 1 and unambiguous else [p for p,m in group])
    return sorted(selected)


def file_filter():
    return 'Supported runs (' + ' '.join('*'+ext for ext in sorted(READERS)) + ')'


def read_many(paths: Iterable[str | Path], workers: int = 4,
              progress: Callable[[int, int, Path], None] | None = None) -> list[Run]:
    """Read several files at once. Import is mostly disk and decompression
    work, so a small thread pool speeds it up without extra processes.
    Results come back in the order given."""
    paths = [Path(p) for p in paths]
    if any(p.suffix.lower() == '.lsc' for p in paths):
        workers = 1  # A native run can need hundreds of MB before processing caches.
    out: list[Run | None] = [None] * len(paths)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(read_run, p): i for i, p in enumerate(paths)}
        for done, fut in enumerate(as_completed(futures), start=1):
            i = futures[fut]
            out[i] = fut.result()
            if progress:
                progress(done, len(paths), paths[i])
    return out  # type: ignore[return-value]


__all__ = ["read_run", "read_many", "read_andi", "read_lsc", "write_andi", "crop_run", "register_reader", "READERS", "discover_runs", "file_filter"]
