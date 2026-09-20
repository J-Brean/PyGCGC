"""🏷️ MS identification: MSP libraries, search, confidence, promote-to-target."""

from types import SimpleNamespace

import numpy as np
import pytest

from pygcgc.core.hydrocarbons import AlkaneHit, AlkaneSeries
from pygcgc.core.identify import (IdentifySettings, LibraryEntry, assign_qvalues,
                                  build_inhouse_entry, build_library, classify, generate_decoys,
                                  identify_blobs, identify_components, load_msp, promote_to_target, retention_index,
                                  search, write_msp)
from pygcgc.core.peaks import BlobSet
from pygcgc.core.targets import load_targets, write_targets_csv


def _entry(name, mz, intensity, ri=None):
    return LibraryEntry(name=name, mz=np.array(mz, dtype=np.int64),
                        intensity=np.array(intensity, dtype=np.float64), ri=ri)


# ------------------------------------------------------------------ MSP
def test_msp_round_trips_through_write_and_load(tmp_path):
    entries = [_entry("Phenanthrene", [178, 176, 152, 89], [999, 120, 30, 15], ri=1771.0),
              _entry("No RI compound", [100, 101], [999, 50])]
    path = write_msp(entries, tmp_path / "lib.msp")
    loaded = load_msp(path)
    assert [e.name for e in loaded] == ["Phenanthrene", "No RI compound"]
    assert loaded[0].ri == pytest.approx(1771.0)
    assert loaded[1].ri is None
    np.testing.assert_array_equal(np.sort(loaded[0].mz), sorted([178, 176, 152, 89]))


def test_load_msp_handles_field_name_variants_and_multi_pair_lines(tmp_path):
    text = (
        "Name: Anthracene\n"
        "Synon: some synonym\n"
        "RetentionIndex: 1774\n"
        "Formula: C14H10\n"
        "CAS#: 120-12-7\n"
        "Num Peaks: 3\n"
        "178 999; 176 80\n"
        "152 40\n"
        "\n"
        "Name: NoPeaksListed\n"
        "RI: 1200\n"
        "Num Peaks: 0\n"
    )
    path = tmp_path / "variants.msp"
    path.write_text(text, encoding="utf-8")
    entries = load_msp(path)
    assert len(entries) == 2
    a = entries[0]
    assert a.name == "Anthracene" and a.ri == pytest.approx(1774.0) and a.formula == "C14H10"
    assert a.cas == "120-12-7"
    assert a.meta.get("Synon") == "some synonym"
    assert sorted(a.mz.tolist()) == [152, 176, 178]
    assert entries[1].mz.size == 0


def test_load_msp_marks_decoys_when_asked(tmp_path):
    path = write_msp([_entry("X", [50], [999])], tmp_path / "x.msp")
    entries = load_msp(path, is_decoy=True)
    assert entries[0].is_decoy


# ------------------------------------------------------------------ library / scoring
def test_build_library_packs_spectra_and_top_ions():
    entries = [_entry("A", [50, 55, 60], [100, 999, 50]), _entry("B", [55, 70], [999, 200])]
    lib = build_library(entries, top_ions=2)
    assert len(lib) == 2
    assert lib.mz_min == 50 and lib.mz_max == 70
    assert lib.matrix.shape == (2, 21)
    assert set(lib.top_ions[0].tolist()) == {55, 50}      # 999 then 100, not 60 (50 > 50... tie broken by order)
    assert set(lib.top_ions[1].tolist()) == {55, 70}


def test_search_finds_its_own_spectrum_as_the_top_hit():
    entries = [_entry("Target", [50, 91, 105], [200, 999, 400], ri=1500.0),
              _entry("Decoy-ish", [51, 92, 106], [150, 900, 300], ri=1520.0)]
    lib = build_library(entries)
    hits = search([50, 91, 105], [200, 999, 400], lib, query_ri=1500.0)
    assert hits and hits[0].name == "Target"
    assert hits[0].forward == pytest.approx(1.0, abs=1e-6)
    assert hits[0].reverse == pytest.approx(1.0, abs=1e-6)


def test_search_ranks_within_ri_window_only():
    entries = [_entry("InWindow", [50, 60], [999, 500], ri=1500.0),
              _entry("FarAway", [50, 60], [999, 500], ri=2000.0)]
    lib = build_library(entries)
    hits = search([50, 60], [999, 500], lib, query_ri=1490.0,
                 settings=IdentifySettings(ri_window=50.0))
    assert [h.name for h in hits] == ["InWindow"]


def test_ri_penalty_separates_near_identical_isomers_by_retention_index():
    # Same spectrum (indistinguishable by mass spec alone), well-separated RI --
    # exactly the phenanthrene/anthracene, BbF/BkF acceptance case the spec asks for.
    entries = [_entry("Phenanthrene", [178, 176, 152], [999, 120, 30], ri=1771.0),
              _entry("Anthracene", [178, 176, 152], [999, 120, 30], ri=1774.0)]
    lib = build_library(entries)
    settings = IdentifySettings(ri_window=50.0, ri_penalty_full=10.0, ri_penalty_weight=0.9)
    hits = search([178, 176, 152], [999, 120, 30], lib, query_ri=1774.0, settings=settings)
    assert hits[0].name == "Anthracene"
    hits2 = search([178, 176, 152], [999, 120, 30], lib, query_ri=1771.0, settings=settings)
    assert hits2[0].name == "Phenanthrene"


def test_search_without_a_retention_index_still_works():
    entries = [_entry("A", [50, 60], [999, 500], ri=1500.0)]
    lib = build_library(entries)
    hits = search([50, 60], [999, 500], lib, query_ri=None)
    assert hits and hits[0].name == "A" and hits[0].ri_diff is None


def test_search_returns_nothing_for_an_empty_library_or_query():
    lib = build_library([])
    assert search([50], [999], lib) == []
    assert search([], [], build_library([_entry("A", [50], [999])])) == []


# ------------------------------------------------------------------ decoys / q-values
def test_generate_decoys_keeps_the_same_ions_and_ri_but_shuffles_intensity():
    entries = [_entry("A", [50, 60, 70], [999, 500, 100], ri=1500.0)]
    decoys = generate_decoys(entries, seed=0)
    assert len(decoys) == 1
    d = decoys[0]
    assert d.is_decoy and d.ri == 1500.0
    assert sorted(d.mz.tolist()) == [50, 60, 70]
    assert sorted(d.intensity.tolist()) == sorted([999.0, 500.0, 100.0])


def test_assign_qvalues_gives_genuine_hits_low_q_and_is_monotone():
    rng = np.random.default_rng(0)
    entries = [_entry(f"C{i}", [50 + i, 60 + i, 70 + i], [999, 500, 200], ri=1000.0 + i * 100)
              for i in range(8)]
    decoys = generate_decoys(entries, seed=1)
    lib = build_library(entries + decoys)

    best_hits = []
    for i, e in enumerate(entries):
        hits = search(e.mz, e.intensity, lib, query_ri=e.ri)
        best_hits.append(hits[0] if hits else None)
    assign_qvalues(best_hits)
    assert all(h is not None and h.name == entries[i].name for i, h in enumerate(best_hits))
    assert all(h.q_value is not None and h.q_value <= 0.05 for h in best_hits)

    # monotonicity: sorted by score descending, q-value never decreases going down the list
    ordered = sorted(best_hits, key=lambda h: h.score, reverse=True)
    qs = [h.q_value for h in ordered]
    assert all(qs[i] <= qs[i + 1] + 1e-12 for i in range(len(qs) - 1))


def test_assign_qvalues_handles_no_hits_gracefully():
    assign_qvalues([None, None])   # must not raise


# ------------------------------------------------------------------ confidence tiers
def test_classify_tiers():
    settings = IdentifySettings(confident_score=0.7, ambiguity=0.02)
    confident = [_hit("A", 0.9, ri_diff=1.0), _hit("B", 0.5, ri_diff=2.0)]
    assert classify(confident, settings) == "2a"

    no_ri = [_hit("A", 0.9, ri_diff=None), _hit("B", 0.4, ri_diff=None)]
    assert classify(no_ri, settings) == "2b"

    ambiguous = [_hit("A", 0.90, ri_diff=1.0), _hit("B", 0.895, ri_diff=1.0)]
    assert classify(ambiguous, settings) == "2b"

    weak = [_hit("A", 0.3, ri_diff=1.0)]
    assert classify(weak, settings) == "3"

    assert classify([], settings) is None


def _hit(name, score, ri_diff):
    from pygcgc.core.identify import SearchHit
    return SearchHit(name=name, score=score, forward=score, reverse=score, ratio=1.0,
                     ri=1500.0 if ri_diff is not None else None, ri_diff=ri_diff, is_decoy=False,
                     entry_index=0)


# ------------------------------------------------------------------ end-to-end identify_blobs
def _make_blob_set(compounds, t1_for_carbon):
    """A hand-built BlobSet (see overlap's own test precedent) carrying
    only what `identify_blobs` reads: mz, apex_spectra, t1."""
    mz_min = min(int(mz.min()) for _n, mz, _i, _c in compounds)
    mz_max = max(int(mz.max()) for _n, mz, _i, _c in compounds)
    mz_axis = np.arange(mz_min, mz_max + 1)
    n = len(compounds)
    apex = np.zeros((n, mz_axis.size))
    t1 = np.zeros(n)
    for i, (_name, mz, intensity, carbon) in enumerate(compounds):
        apex[i, mz - mz_min] = intensity
        t1[i] = t1_for_carbon(carbon)
    zeros = np.zeros(n)
    return BlobSet(mz=mz_axis, t1=t1, t2=zeros + 1.0, col=zeros.astype(int), row=zeros.astype(int),
                  height=zeros + 1.0, snr=zeros + 10.0, volume=apex.sum(axis=1), n_pixels=zeros + 10,
                  n_modulations=zeros + 2, spectra=apex.copy(), labels=np.zeros((1, 1), dtype=np.int32),
                  wrap_rows=0, pitch=0.02, period=4.0, t1_start=0.0, apex_spectra=apex)


def test_identify_blobs_recovers_at_least_90_percent_of_synthetic_compounds():
    rng = np.random.default_rng(0)
    carbons = list(range(12, 32, 2))                  # 10 compounds, RI 1200..3100
    entries = []
    for k, carbon in enumerate(carbons):
        mz = rng.choice(np.arange(45, 250), size=6, replace=False)
        intensity = rng.uniform(50, 999, size=6)
        intensity[0] = 999.0                          # a clear base peak
        entries.append(_entry(f"Compound{k}", mz, intensity, ri=float(carbon * 100)))
    library = build_library(entries + generate_decoys(entries, seed=1))

    series = AlkaneSeries(hits={10: AlkaneHit(carbon=10, t1=100.0, t2=1.0, height=1e5, qual_ratios=()),
                                30: AlkaneHit(carbon=30, t1=300.0, t2=1.0, height=1e5, qual_ratios=())})
    t1_for_carbon = lambda c: 100.0 + (c - 10) * 10.0    # matches the series' own step (10 s/carbon)

    compounds = []
    for k, (carbon, e) in enumerate(zip(carbons, entries)):
        intensity = e.intensity.copy()
        if k % 3 == 0:                                # perturb a third of them: one extra noise peak
            intensity = intensity * (1.0 + rng.normal(0, 0.05, intensity.size))
        compounds.append((e.name, e.mz, np.clip(intensity, 0, None), carbon))
    blobs = _make_blob_set(compounds, t1_for_carbon)

    results = identify_blobs(blobs, library, series)
    correct = sum(1 for i, r in enumerate(results) if r.best is not None and r.best.name == entries[i].name)
    assert correct / len(entries) >= 0.9
    assert all(r.best is None or r.best.tier in ("2a", "2b", "3") for r in results)


def test_identify_blobs_works_without_an_alkane_series():
    entries = [_entry("Solo", [50, 60], [999, 500])]
    library = build_library(entries)
    blobs = _make_blob_set([("Solo", np.array([50, 60]), np.array([999.0, 500.0]), 15)],
                           lambda c: 150.0)
    results = identify_blobs(blobs, library, series=None)
    assert results[0].best is not None and results[0].best.name == "Solo"
    assert results[0].best.ri_diff is None


def test_identify_components_uses_each_separated_spectrum():
    library = build_library([
        _entry("A", [50, 51], [999, 100]),
        _entry("B", [50, 51], [80, 999])])
    component_set = SimpleNamespace(components=[
        SimpleNamespace(mz=np.array([50, 51]), spectrum=np.array([999.0, 100.0]),
                        t1=100.0, t2=1.0),
        SimpleNamespace(mz=np.array([50, 51]), spectrum=np.array([80.0, 999.0]),
                        t1=100.0, t2=1.1)])
    results = identify_components(component_set, library)
    assert [result.best.name for result in results] == ["A", "B"]


# ------------------------------------------------------------------ promote / in-house
def test_promote_to_target_reads_quant_and_qualifier_ions_from_the_entry():
    entry = _entry("Fluoranthene", [202, 101, 200, 88], [999, 200, 150, 50], ri=1610.0)
    hit = search(entry.mz, entry.intensity, build_library([entry]), query_ri=1610.0)[0]
    target = promote_to_target(entry, hit, t1=250.0, t2=1.8, n_qualifiers=2)
    assert target.compound == "Fluoranthene"
    assert target.quant_ion == 202
    assert set(target.qual_ions) == {101, 200}
    assert target.t1 == 250.0 and target.t2 == 1.8
    assert target.ri == pytest.approx(1610.0)


def test_promote_to_target_rejects_an_empty_entry():
    empty = LibraryEntry(name="Empty", mz=np.array([], dtype=np.int64), intensity=np.array([]))
    with pytest.raises(ValueError):
        promote_to_target(empty, None, t1=0.0, t2=0.0)


def test_write_targets_csv_round_trips_through_load_targets(tmp_path):
    entry = _entry("Fluoranthene", [202, 101, 200, 88], [999, 200, 150, 50], ri=1610.0)
    hit = search(entry.mz, entry.intensity, build_library([entry]), query_ri=1610.0)[0]
    target = promote_to_target(entry, hit, t1=250.0, t2=1.8)
    path = write_targets_csv([target], tmp_path / "target_table.csv")
    loaded = load_targets(path)
    assert len(loaded) == 1
    assert loaded[0].compound == "Fluoranthene"
    assert loaded[0].quant_ion == 202
    assert loaded[0].ri == pytest.approx(1610.0)


def test_write_targets_csv_can_append_to_an_existing_table(tmp_path):
    first = _entry("A", [50], [999], ri=1000.0)
    second = _entry("B", [60], [999], ri=1100.0)
    path = tmp_path / "target_table.csv"
    write_targets_csv([promote_to_target(first, None, t1=100.0, t2=1.0)], path)
    write_targets_csv([promote_to_target(second, None, t1=110.0, t2=1.0)], path, append=True)
    loaded = load_targets(path)
    assert [t.compound for t in loaded] == ["A", "B"]


def test_build_inhouse_entry_extracts_nonzero_channels():
    spectrum = np.zeros(20)
    spectrum[5] = 999.0
    spectrum[12] = 300.0
    entry = build_inhouse_entry("Standard1", spectrum, mz0=40, ri=1234.0)
    assert sorted(entry.mz.tolist()) == [45, 52]
    assert entry.ri == 1234.0
    assert not entry.is_decoy


# ------------------------------------------------------------------ retention index
def test_retention_index_matches_100_times_carbon_of():
    series = AlkaneSeries(hits={10: AlkaneHit(carbon=10, t1=100.0, t2=1.0, height=1.0, qual_ratios=()),
                                20: AlkaneHit(carbon=20, t1=200.0, t2=1.0, height=1.0, qual_ratios=())})
    assert retention_index(series, 150.0) == pytest.approx(100.0 * series.carbon_of(150.0))
    assert retention_index(series, 150.0) == pytest.approx(1500.0)


# ------------------------------------------------------------------ dummy data
@pytest.mark.dummy
def test_dummy_library_search_finds_known_targets(dummy_dir, dummy_runs, dummy_truth):
    """The milestone's own acceptance bar (>=90% correct top hits) against
    the dummy pack's own `synthetic_library.msp`, on whichever of the
    already-confirmed target names (see test_siloxane.py's own dummy
    test) also turn up in the library -- library compound-name spelling
    hasn't been checked against the truth files' own naming, so this
    skips rather than guesses if too few line up. RI isn't used here
    (finding a real alkane series blind, without the dummy pack to check
    naming/positions against, risked a fragile test); a full RI-aware
    check is worth adding once this file is actually to hand."""
    from pygcgc.core.background import remove_background
    from pygcgc.core.fold import fold
    from pygcgc.core.peaks import find_blobs
    from tests.conftest import PERIOD, TRUE_PHASE

    library_path = dummy_dir / "synthetic_library.msp"
    if not library_path.exists():
        pytest.skip("synthetic_library.msp not in this dummy pack")
    entries = load_msp(library_path)
    assert entries, "synthetic_library.msp parsed to zero entries; check the MSP format assumptions"
    by_name = {e.name.lower(): e for e in entries}

    run = dummy_runs["CAL_L4"]
    f = fold(run, PERIOD, TRUE_PHASE)
    clean = remove_background(run, f).run
    blobs = find_blobs(clean, f)

    candidates = ("naphthalene", "fluorene", "pyrene", "1-nitropyrene", "behenic acid")
    known = [c for c in candidates if ("CAL_L4", c) in dummy_truth and c.lower() in by_name]
    if len(known) < 3:
        pytest.skip("too few confirmed target names match the library's own naming to test this "
                    "meaningfully; update the naming once the dummy pack is to hand")

    library = build_library(entries + generate_decoys(entries, seed=0))
    correct = 0
    for name in known:
        t = dummy_truth[("CAL_L4", name)]
        col, t2 = f.to_fold(float(t["apex_scan_time_s"]))
        t1 = float(f.t1[col] + 0.5 * f.spec.period)
        i = blobs.nearest(t1, t2)
        spectrum = blobs.apex_spectra[i]
        mz0 = int(blobs.mz[0])
        nz = np.flatnonzero(spectrum > 0)
        hits = search(nz + mz0, spectrum[nz], library, query_ri=None)
        if hits and hits[0].name.lower() == name:
            correct += 1
    assert correct / len(known) >= 0.9
