"""Checking typed compound names against the library's spelling."""

import numpy as np
import pytest

from pygcgc.core.identify import LibraryEntry
from pygcgc.core.naming import check_name, format_cas, name_key, part_key, tidy
from pygcgc.core.nist_sdf import compact_from_entries


def library():
    names = [("Toluene", "108-88-3", "C7H8"), ("Benzene, 1,2-dimethyl-", "95-47-6", "C8H10"),
             ("Benzene, 1,3-dimethyl-", "108-38-3", "C8H10"), ("Decane", "124-18-5", "C10H22"),
             ("Naphthalene, 2-methyl-", "91-57-6", "C11H10"), ("Naphthalene, 1-methyl-", "90-12-0", "C11H10")]
    lib = compact_from_entries([LibraryEntry(n, np.array([91]), np.array([999.0]), cas=c, formula=f) for n, c, f in names])
    return lib


def test_capitals_commas_dashes_and_spaces_do_not_make_a_new_name():
    assert name_key("Benzene, 1,2-dimethyl-") == name_key("BENZENE 1 2 dimethyl") == "benzene12dimethyl"
    assert name_key("1,2-Dimethylbenzene") == name_key("1 2 dimethylbenzene") == name_key("1,2dimethylbenzene")
    assert tidy("  Toluene ,  ") == "Toluene" and tidy("a   b\t c") == "a b c"
    assert name_key("1,2-dimethyl") != name_key("1,3-dimethyl")                 # locants are part of the name


def test_cas_numbers_are_shown_with_hyphens():
    assert format_cas("108883") == "108-88-3" and format_cas("7782390") == "7782-39-0" and format_cas("") == ""


def test_comma_separated_parts_can_come_in_any_order():
    assert part_key("Benzene, 1,2-dimethyl-") == part_key("1,2-dimethyl-, benzene")
    assert part_key("Benzene, 1,2-dimethyl-") != part_key("Benzene, 1,3-dimethyl-")


@pytest.mark.parametrize("typed", ["toluene", "TOLUENE", " Toluene ", "tolu-ene", "Tolu ene"])
def test_the_library_spelling_replaces_a_respelling(typed):
    lib = library()
    c = check_name(typed, lib)
    assert c.name == "Toluene" and c.cas == "108883" and c.formula == "C7H8"
    assert c.status in ("respelled", "exact")
    assert (c.status == "exact") == (typed.strip() == "Toluene")      # padding is tidied away, so it is still exact


def test_exact_names_are_reported_as_such_and_other_orders_are_recognised():
    lib = library()
    assert check_name("Toluene", lib).status == "exact"
    c = check_name("1,2-dimethyl-, benzene", lib)
    assert c.status == "respelled" and c.name == "Benzene, 1,2-dimethyl-" and c.cas == "95476"
    c = check_name("benzene 1 3 dimethyl", lib)
    assert c.name == "Benzene, 1,3-dimethyl-"                                     # not the 1,2 compound


def test_a_new_name_is_kept_as_typed_tidied_and_called_out_with_close_suggestions():
    lib = library()
    c = check_name("  Naphthalene,   2-ethyl-", lib, candidates=["Naphthalene, 2-methyl-", "Decane"])
    assert c.status == "new" and c.name == "Naphthalene, 2-ethyl-" and c.cas == ""
    assert any("Not in the library" in n for n in c.notes)
    c = check_name("Decan", lib)
    assert c.status == "new" and "Decane" in c.suggestions
    typo = check_name("Tolune", lib, candidates=["Toluene"])
    assert typo.status == "new" and "Toluene" in typo.suggestions


def test_an_empty_name_is_refused_and_a_missing_library_is_said_not_assumed():
    assert check_name("   ", library()).status == "empty" and check_name(" , ", None).problems
    c = check_name("Whatever", None)
    assert c.status == "new" and c.name == "Whatever" and "No library is loaded" in c.notes[0]


def test_a_name_already_used_by_another_peak_is_pointed_out():
    lib = library()
    c = check_name("toluene", lib, taken=[("peak 4", "Decane"), ("peak 7", "TOLUENE")])
    assert any("peak 7" in n and "already called" in n for n in c.notes)
    assert not any("already called" in n for n in check_name("toluene", lib, taken=[("peak 4", "Decane")]).notes)


def test_two_library_entries_with_one_key_but_different_cas_are_flagged():
    lib = compact_from_entries([LibraryEntry("Xylene", np.array([91]), np.array([1.0]), cas="1330-20-7"),
                                LibraryEntry("XYLENE", np.array([91]), np.array([1.0]), cas="95-47-6")])
    c = check_name("xylene", lib)
    assert len(c.matches) == 2 and any("Different CAS numbers" in n for n in c.notes)


def test_the_index_is_built_once_per_library():
    lib = library()
    check_name("toluene", lib)
    first = lib._cache["name_index"]
    check_name("decane", lib)
    assert lib._cache["name_index"] is first
