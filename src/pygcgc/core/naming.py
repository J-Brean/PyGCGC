"""Checking a compound name someone types against the library's spelling.

The idea in plain terms
------------------------
The same compound is written many ways: `Toluene`, `toluene`, `TOLUENE`; `1,2-dimethylbenzene`,
`1,2 dimethylbenzene`, `1 2 dimethylbenzene`; NIST's own `Benzene, 1,2-dimethyl-`. When a person
renames a peak, the name should not become a new, unrelated string just because of capitals,
commas, dashes or spaces. So a typed name is compared with the library's names by a *key*: lower
case, with everything that is not a letter or a digit removed. Names with the same key are the same
compound written differently, and the library's spelling is offered (and used by default). Names
whose comma-separated parts are the same in another order (`Benzene, 1,2-dimethyl-` and
`1,2-dimethyl-, benzene`) also match. Digits stay in the key, so `1,2-` and `1,3-` are different
compounds. A name the library does not know is kept as typed, with whitespace tidied, and called
out as new; close library names are suggested in case it was a slip.

Nothing here knows about Qt.
"""

from __future__ import annotations

import difflib
import re
from dataclasses import dataclass, field


def format_cas(digits: str) -> str:
    """`108883` as `108-88-3` (the library keeps CAS numbers as digits only)."""
    d = re.sub(r"\D", "", digits or "")
    return f"{d[:-3]}-{d[-3:-1]}-{d[-1]}" if len(d) >= 5 else d


def tidy(name: str) -> str:
    """Whitespace collapsed, leading and trailing spaces and stray commas removed. Case is kept."""
    return re.sub(r"\s+", " ", (name or "").strip()).strip(" ,;")


def name_key(name: str) -> str:
    """Lower case letters and digits only: capitals, spaces, commas, dashes and brackets do not matter."""
    return re.sub(r"[^a-z0-9]", "", (name or "").lower())


def part_key(name: str) -> str:
    """Like `name_key`, but the comma-separated parts may come in any order."""
    parts = [name_key(p) for p in re.split(r"\s*,\s+|\s*;\s*", name or "") if name_key(p)]
    return "|".join(sorted(parts))


@dataclass
class NameCheck:
    status: str                         # "empty", "exact", "respelled" or "new"
    name: str = ""                      # what would be saved
    cas: str = ""
    formula: str = ""
    matches: list = field(default_factory=list)      # library spellings with the same key: (index, name, cas, formula)
    suggestions: list = field(default_factory=list)  # close library names, when the name is new
    notes: list = field(default_factory=list)
    problems: list = field(default_factory=list)     # things that stop it being saved


def _index(library) -> tuple:
    """`(key -> entry indices, part key -> entry indices)` for a library, built once and kept on it."""
    cache = getattr(library, "_cache", None)
    hit = cache.get("name_index") if cache is not None else None
    if hit is not None:
        return hit
    by_key: dict = {}
    by_part: dict = {}
    for i, name in enumerate(library.names):
        by_key.setdefault(name_key(name), []).append(i)
        by_part.setdefault(part_key(name), []).append(i)
    if cache is not None:
        cache["name_index"] = (by_key, by_part)
    return by_key, by_part


def check_name(typed: str, library=None, candidates=(), taken=()) -> NameCheck:
    """Compare a typed name with the library and with the names already used.

    `candidates` are names the peak's own matches suggested (offered as spelling hints when the name is new);
    `taken` is `(peak label, name)` for the other peaks already named in this run."""
    name = tidy(typed)
    if not name_key(name):
        return NameCheck("empty", problems=["Type a name."])
    check = NameCheck("new", name=name)
    key = name_key(name)
    if library is None or len(library) == 0:
        check.notes.append("No library is loaded, so the spelling could not be checked; it is saved as typed.")
    else:
        by_key, by_part = _index(library)
        hits = by_key.get(key) or by_part.get(part_key(name)) or []
        if hits:
            spellings = []
            for i in hits:
                entry = (i, library.names[i], library.cas[i], library.formula[i])
                if entry[1] not in [s[1] for s in spellings]:
                    spellings.append(entry)
            check.matches = spellings
            best = spellings[0]
            check.name, check.cas, check.formula = best[1], best[2], best[3]
            if best[1] == name:
                check.status = "exact"
                check.notes.append(f"In the library as typed (CAS {format_cas(best[2]) or 'unknown'}).")
            else:
                check.status = "respelled"
                check.notes.append(f"The same compound as the library's \"{best[1]}\" (they differ only in capitals, "
                                   "punctuation, spacing or the order of comma-separated parts): the library spelling is used.")
            if len(spellings) > 1:
                others = "; ".join(f"{s[1]} (CAS {format_cas(s[2]) or 'unknown'})" for s in spellings[1:4])
                check.notes.append(f"The library also has: {others}. Different CAS numbers mean different compounds: check which you mean.")
        else:
            check.notes.append("Not in the library: it will be saved as typed (whitespace tidied). Check the spelling if it should be there.")
            pool = {name_key(c): c for c in candidates if name_key(c)}
            close = [pool[k] for k in difflib.get_close_matches(key, list(pool), n=3, cutoff=0.75)]
            starts = []
            if len(key) >= 4:
                starts = [library.names[v[0]] for k, v in by_key.items() if k.startswith(key) and k != key]
                starts = sorted(set(starts), key=lambda s: (len(s), s))[:5]
            check.suggestions = list(dict.fromkeys(close + starts))
    for label, other in taken:
        if name_key(other) == name_key(check.name) or name_key(other) == key:
            check.notes.append(f"Another peak in this run is already called \"{other}\" (peak {label}).")
            break
    return check
