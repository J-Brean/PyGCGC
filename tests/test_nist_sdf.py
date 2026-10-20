"""Reading NIST-style SDF files: parsing, the compact library, the RI link and the cache."""

import os
from pathlib import Path

import numpy as np
import pytest

from tests.nist_fixture import MAIN, REP, RI, write_mini_nist
from pygcgc.core.nist_sdf import (DEFAULT_COLUMN_CLASS, CompactLibrary, Cancelled, attach_experimental_ri,
                                  clean_cas, clean_name, compact_from_entries, find_nist_files, iter_sdf_records,
                                  load_compact, open_ei_library, parse_peaks, read_ei_library, read_ri_records,
                                  save_compact, source_key)
from pygcgc.core.identify import LibraryEntry


@pytest.fixture
def nist(tmp_path):
    return write_mini_nist(tmp_path)


def test_records_survive_crlf_odd_headers_and_missing_fields(nist):
    records = list(iter_sdf_records(nist["mainlib"]))
    assert len(records) == len(MAIN)
    names = [rec["NAME"] for _t, rec in records]
    assert names[0] == "Decane" and names[-1] == "Compound without CAS"
    assert records[1][1]["CASNO"] == "1120214"            # "> <CASNO>" with one space
    assert "CASNO" not in records[-1][1]                  # empty fields are omitted
    assert records[0][1]["SYNONYMS"] == "a synonym\nanother synonym"
    assert "\r" not in records[0][1]["MASS SPECTRAL PEAKS"]


def test_field_filter_keeps_only_the_wanted_fields(nist):
    _title, rec = next(iter_sdf_records(nist["mainlib"], fields=("name", "formula")))
    assert set(rec) == {"NAME", "FORMULA"}


def test_peaks_read_from_one_pair_per_line_and_several_pairs_per_line():
    mz, inten = parse_peaks("43 999\n57 850\n71 500")
    assert mz.tolist() == [43, 57, 71] and inten.tolist() == [999, 850, 500]
    mz2, inten2 = parse_peaks("43 999 57 850\n71 500 85 300")
    assert mz2.tolist() == [43, 57, 71, 85] and inten2.tolist() == [999, 850, 500, 300]
    assert parse_peaks("")[0].size == 0
    assert parse_peaks("5 0\n6 10")[0].tolist() == [6]    # zero intensities dropped


def test_key_helpers_ignore_punctuation():
    assert clean_cas("7782-39-0") == "7782390"
    assert clean_name("Benzene, 1-chloro-4-nitro-") == clean_name("BENZENE 1 chloro 4 nitro")


def test_ei_library_has_one_row_per_spectrum_and_keeps_replicates(nist):
    lib = read_ei_library([nist["mainlib"], nist["replib"]])
    assert len(lib) == len(MAIN) + len(REP)
    assert lib.source.tolist() == [0] * len(MAIN) + [1]
    assert lib.cas[-2] == "" and lib.cas[2] == "108883"
    assert np.isnan(lib.ri_est[-2]) and lib.ri_est[0] == 1000 and lib.ri_est_err[0] == 20
    entry = lib.entry(0)
    assert entry.name == "Decane" and entry.mz.tolist() == [43, 57, 71, 85, 142]
    assert entry.intensity.max() == 999 and entry.ri == 1000          # falls back to the estimate
    # replicates of the same compound share a key
    assert lib.compound_key(2) == lib.compound_key(len(MAIN))


def test_cancelling_stops_a_read(nist):
    def stop():
        raise Cancelled
    big = nist["mainlib"].parent / "big.sdf"
    big.write_bytes(nist["mainlib"].read_bytes() * 1000)
    with pytest.raises(Cancelled):
        read_ei_library([big], cancel=stop)


def test_experimental_ri_uses_the_chosen_column_class_and_programme(nist):
    lib = read_ei_library([nist["mainlib"]])
    records = read_ri_records(nist["nist_ri"])
    assert len(records) == len(RI)
    counts = attach_experimental_ri(lib, records, DEFAULT_COLUMN_CLASS)
    decane = 0
    # three usable semi-standard ramp/complex records: 1000, 1001, 999 (the isothermal 1040 and polar 1000 are out)
    assert lib.ri_exp_n[decane] == 3 and lib.ri_exp[decane] == pytest.approx(1000.0)
    assert lib.ri_exp_sigma[decane] >= 5.0
    toluene = 2
    assert lib.ri_exp_n[toluene] == 2 and lib.ri_exp[toluene] == pytest.approx(768.0)
    assert lib.ri_exp_sigma[toluene] == 15.0                          # fewer than three: no usable spread
    assert lib.ri_link[toluene] == "cas"
    assert counts["cas"] >= 4
    # the entry with no CAS is linked by its normalised name, whatever the punctuation
    nocas = 5
    assert lib.ri_link[nocas] == "name" and lib.ri_exp[nocas] == pytest.approx(880.0)
    polar = attach_experimental_ri(lib, records, "Standard polar")
    assert lib.ri_exp_n[decane] == 1 and np.isnan(lib.ri_exp[toluene])
    assert polar["none"] > 0


def test_n_alkane_indices_come_from_the_scale_definition_for_any_column(nist):
    lib = read_ei_library([nist["mainlib"]])
    records = read_ri_records(nist["nist_ri"])
    undecane = 1
    for column in ("Semi-standard non-polar", "Standard polar"):
        attach_experimental_ri(lib, records, column)
        assert lib.ri_exp[undecane] == pytest.approx(1100.0) and lib.ri_exp_sigma[undecane] == pytest.approx(1.0)


def test_a_nist_number_links_when_the_cas_number_does_not(nist):
    lib = read_ei_library([nist["mainlib"]])
    records = read_ri_records(nist["nist_ri"])
    records.cas = ["" for _ in records.cas]                             # no CAS in the RI file
    counts = attach_experimental_ri(lib, records)
    assert counts["nistno"] >= 3 and lib.ri_link[0] == "nistno"


def test_cache_round_trip_and_invalidation(nist, tmp_path):
    cache = tmp_path / "cache"
    lib = open_ei_library([nist["mainlib"]], cache_root=cache)
    assert (cache / ("ei_" + lib.key) / "spectra.npz").exists()
    again = open_ei_library([nist["mainlib"]], cache_root=cache)
    assert again.names == lib.names and (again.matrix != lib.matrix).nnz == 0
    assert again.cas == lib.cas and np.array_equal(again.nistno, lib.nistno)
    before = source_key([nist["mainlib"]])
    nist["mainlib"].write_bytes(nist["mainlib"].read_bytes() + b"\r\n")
    os.utime(nist["mainlib"], (1, 1))
    assert source_key([nist["mainlib"]]) != before


def test_a_damaged_cache_is_rebuilt(nist, tmp_path):
    cache = tmp_path / "cache"
    lib = open_ei_library([nist["mainlib"]], cache_root=cache)
    (cache / ("ei_" + lib.key) / "spectra.npz").write_bytes(b"not a zip")
    assert len(open_ei_library([nist["mainlib"]], cache_root=cache)) == len(MAIN)


def test_msp_entries_pack_into_a_compact_library():
    entries = [LibraryEntry("A", np.array([40, 50]), np.array([10.0, 999.0]), ri=900.0, cas="12-34-5"),
               LibraryEntry("B", np.array([41]), np.array([999.0]), ri=None),
               LibraryEntry("DECOY_A", np.array([40]), np.array([1.0]), is_decoy=True)]
    lib = compact_from_entries(entries)
    assert len(lib) == 2 and lib.cas[0] == "12345" and lib.ri_link[0] == "library"
    assert lib.ri_exp[0] == 900 and np.isnan(lib.ri_exp[1])


def test_folder_discovery_finds_the_three_files(nist, tmp_path):
    found = find_nist_files(tmp_path)
    assert set(found) == {"mainlib", "replib", "nist_ri"}
    assert find_nist_files(tmp_path / "SDF").keys() == found.keys()


@pytest.mark.real_data
def test_real_library_counts():
    folder = os.environ.get("PYGCGC_NIST_SDF")
    if not folder or not Path(folder).exists():
        pytest.skip("set PYGCGC_NIST_SDF to the folder holding mainlib.sdf to run this")
    found = find_nist_files(folder)
    lib = read_ei_library([found["mainlib"]])
    assert len(lib) == 212961 and isinstance(lib, CompactLibrary)
