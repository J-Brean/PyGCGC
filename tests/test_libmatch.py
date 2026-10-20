"""Library matching by spectrum and retention index."""

import numpy as np
import pytest

from tests.nist_fixture import write_mini_nist
from pygcgc.core.hydrocarbons import AlkaneHit, AlkaneSeries
from pygcgc.core.identify import LibraryEntry
from pygcgc.core.libmatch import (FLAGS, MatchSettings, Query, match_queries, peak_ri, q_at_scores,
                                  queries_from_peaks)
from pygcgc.core.nist_sdf import (attach_experimental_ri, compact_from_entries, read_ei_library, read_ri_records)


@pytest.fixture
def lib(tmp_path):
    files = write_mini_nist(tmp_path)
    library = read_ei_library([files["mainlib"], files["replib"]])
    attach_experimental_ri(library, read_ri_records(files["nist_ri"]))
    return library


def query_from(lib, i, ri=float("nan"), sigma=10.0, flags=()):
    row = lib.matrix.getrow(i)
    return Query(mz=row.indices.astype(np.int64), intensity=row.data.astype(np.float64), ri=ri, ri_sigma=sigma,
                 flags=tuple(flags), t1=100.0, t2=2.0, volume=1.0)


def ladder():
    """C10 at 100 s, C11 at 200 s, C12 at 300 s: one carbon per 100 s."""
    s = AlkaneSeries()
    for c, t in ((10, 100.0), (11, 200.0), (12, 300.0)):
        s.hits[c] = AlkaneHit(carbon=c, t1=t, t2=2.0, height=1.0, qual_ratios=())
    return s


# ------------------------------------------------------------------ retention index of a peak
def test_peak_ri_follows_the_ladder_and_carries_the_time_uncertainty():
    ri, sigma, extrapolated = peak_ri(ladder(), 150.0, sigma_t1=2.0)
    assert ri == pytest.approx(1050.0) and sigma == pytest.approx(2.0) and not extrapolated    # 100 RI per 100 s
    assert peak_ri(ladder(), 350.0, 2.0)[2] and peak_ri(ladder(), 50.0, 2.0)[2]


def test_queries_from_peaks_flag_missing_ladders_borrowed_ladders_and_coelution():
    class Peaks:
        t1 = np.array([150.0, 400.0]); t2 = np.array([2.0, 2.0]); volume = np.array([5.0, 1.0])
        mz = np.arange(40, 60); spectra = np.ones((2, 20)); co_elution = np.array([False, True])
        component_id = ["a", "b"]; unresolved = np.array([False, False])

        def __len__(self):
            return 2
    plain = queries_from_peaks(Peaks(), ladder(), period=8.0)
    assert plain[0].ri == pytest.approx(1050.0) and plain[0].flags == ()
    assert plain[1].flags == ("RI_EXTRAPOLATED", "COELUTION")
    assert plain[0].ri_sigma == pytest.approx(0.25 * 8.0)             # a quarter of a modulation, 1 RI per second
    assert queries_from_peaks(Peaks(), ladder(), 8.0, borrowed=True)[0].flags == ("RI_BORROWED",)
    none = queries_from_peaks(Peaks(), None)
    assert np.isnan(none[0].ri) and "RI_MISSING" in none[0].flags


# ------------------------------------------------------------------ q-values
def test_q_values_come_from_decoys_and_never_improve_going_the_wrong_way():
    scores = [0.95, 0.9, 0.85, 0.6, 0.5, 0.4]
    decoy = [False, False, False, True, False, True]
    q = q_at_scores(scores, decoy, scores)
    assert q[0] == q[1] == q[2] == 0.0
    # lowest raw rate at or below each score: 2/4 at 0.4, 1/4 at 0.5, then 1/3 at 0.6 (kept at the lower 1/4)
    assert q[3] == pytest.approx(1 / 4)
    assert q[4] == pytest.approx(1 / 4) and q[5] == pytest.approx(2 / 4)
    assert np.all(np.diff(q) >= -1e-12)                              # worse scores never have a lower q
    assert np.isnan(q_at_scores([], [], [0.5])[0])


# ------------------------------------------------------------------ search
def test_ri_separates_isomers_the_spectrum_cannot(lib):
    """Propylbenzene (RI 950) and its look-alike (RI 1050) have almost identical spectra."""
    settings = MatchSettings(use_decoys=False)
    propyl = 3
    q = query_from(lib, propyl, ri=952.0)
    match = match_queries([q], lib, settings)[0]
    assert match.best.name == "Propylbenzene" and match.best.in_gate
    assert match.best.mf > 990
    assert not any(c.name.startswith("Isopropylbenzene") and c.in_gate for c in match.candidates)
    # the other way round: at RI 1050 the look-alike wins
    other = match_queries([query_from(lib, propyl, ri=1049.0)], lib, settings)[0]
    assert other.best.name.startswith("Isopropylbenzene")
    # with no RI at all the two are ranked by spectrum and flagged as hard to tell apart
    blind = match_queries([query_from(lib, propyl, flags=("RI_MISSING",))], lib, settings)[0]
    assert "ISOMER_AMBIGUOUS" in blind.flags and "RI_MISSING" in blind.flags and blind.level != "confident"


def test_a_clean_experimental_match_is_confident_and_sentences_exist_for_every_flag(lib):
    match = match_queries([query_from(lib, 0, ri=1001.0)], lib, MatchSettings(use_decoys=False))[0]
    assert match.best.name == "Decane" and match.best.ri_kind == "experimental" and match.best.n_ri == 3
    assert match.level == "confident" and match.tier == "2a" and abs(match.best.z) < 1
    assert match.flags == []
    assert all(isinstance(text, str) and text for text in FLAGS.values())


def test_replicate_spectra_collapse_to_one_candidate_per_compound(lib):
    match = match_queries([query_from(lib, 2, ri=768.0)], lib, MatchSettings(use_decoys=False))[0]
    assert [c.name for c in match.candidates].count("Toluene") == 1


def test_an_estimated_ri_is_flagged_and_cannot_make_a_peak_confident(lib):
    # Remove the experimental values so only NIST's estimate (Decane 1000, error 20) is left.
    lib.ri_exp[:] = np.nan; lib.ri_exp_sigma[:] = np.nan; lib.ri_exp_n[:] = 0; lib._cache.clear()
    match = match_queries([query_from(lib, 0, ri=1002.0)], lib, MatchSettings(use_decoys=False))[0]
    assert match.best.ri_kind == "estimated" and "RI_ESTIMATED" in match.flags
    assert match.level == "probable"
    polar = match_queries([query_from(lib, 0, ri=1002.0)], lib,
                          MatchSettings(use_decoys=False, column_class="Standard polar"))[0]
    assert polar.best.ri_kind == "none"                             # the estimate describes non-polar columns only


def test_the_estimate_is_corrected_for_the_column_and_dropped_when_nist_doubts_it(lib):
    """NIST's estimate describes a standard non-polar column; a 5 % phenyl column reads 21 units higher."""
    from pygcgc.core.libmatch import library_ri
    lib.ri_exp[:] = np.nan; lib.ri_exp_sigma[:] = np.nan; lib.ri_exp_n[:] = 0; lib._cache.clear()
    ri, sigma, kind = library_ri(lib, MatchSettings())
    assert kind[0] == 2 and ri[0] == pytest.approx(1021.0) and sigma[0] == pytest.approx(20 * 1.25)
    ri_std, sigma_std, _ = library_ri(lib, MatchSettings(column_class="Standard non-polar"))
    assert ri_std[0] == pytest.approx(1000.0) and sigma_std[0] == pytest.approx(20.0)
    lib.ri_est_err[0] = 150.0; lib._cache.clear()                    # NIST itself says this one is poor
    assert library_ri(lib, MatchSettings())[2][0] == 0
    assert library_ri(lib, MatchSettings(use_estimated_ri=False))[2][1] == 0


def test_a_good_spectrum_outside_the_ri_gate_is_shown_demoted_and_flagged(lib):
    q = query_from(lib, 0, ri=1400.0)                               # Decane's spectrum, 400 RI units too late
    match = match_queries([q], lib, MatchSettings(use_decoys=False))[0]
    assert "RI_CONFLICT" in match.flags and match.level == "unknown"
    assert any((not c.in_gate) and c.name == "Decane" for c in match.candidates)
    assert not match.best.in_gate and match.best.mf > 990      # nothing inside the gate: the demoted match is still shown
    assert "RI_TENSION" not in match.flags                     # RI_CONFLICT already says it


def test_a_library_entry_with_no_ri_is_still_a_candidate_but_flagged(lib):
    nocas = 5
    lib.ri_exp[nocas] = np.nan; lib.ri_exp_n[nocas] = 0; lib.ri_est[nocas] = np.nan; lib._cache.clear()
    match = match_queries([query_from(lib, nocas, ri=880.0)], lib, MatchSettings(use_decoys=False))[0]
    assert match.best.name == "Compound without CAS" and "RI_LIBRARY_MISSING" in match.flags
    assert match.level == "probable"


def test_peaks_with_few_ions_or_shared_fits_are_flagged(lib):
    q = Query(mz=np.array([91, 92]), intensity=np.array([999.0, 600.0]), ri=768.0, ri_sigma=5.0, flags=("COELUTION",))
    match = match_queries([q], lib, MatchSettings(use_decoys=False))[0]
    assert {"FEW_IONS", "COELUTION"} <= set(match.flags)
    assert match.level != "confident"


def test_ions_outside_the_acquired_range_do_not_count_against_a_library_spectrum(lib):
    """Decane has its molecular ion at m/z 142. A run that stops at 100 never saw it."""
    row = lib.matrix.getrow(0)
    keep = row.indices < 100
    q = Query(mz=row.indices[keep].astype(np.int64), intensity=row.data[keep].astype(np.float64), ri=1000.0,
              ri_sigma=5.0)
    settings = MatchSettings(use_decoys=False)
    short = match_queries([q], lib, settings, mz_range=(35, 100))[0].best.mf
    wide = match_queries([q], lib, settings, mz_range=(35, 200))[0].best.mf
    assert short > wide and short > 990


def test_excluded_ions_are_ignored_in_both_spectra(lib):
    row = lib.matrix.getrow(0)
    data = row.data.astype(np.float64).copy()
    data[row.indices == 57] = 1.0                                   # column bleed or similar at 57
    q = Query(mz=row.indices.astype(np.int64), intensity=data, ri=1000.0, ri_sigma=5.0)
    plain = match_queries([q], lib, MatchSettings(use_decoys=False))[0].best.mf
    cleaned = match_queries([q], lib, MatchSettings(use_decoys=False, exclude_ions=(57,)))[0].best.mf
    assert cleaned > plain


def test_vectorised_search_equals_a_slow_reference(tmp_path):
    rng = np.random.default_rng(3)
    entries = []
    for i in range(300):
        mz = np.unique(rng.integers(30, 250, size=int(rng.integers(6, 40))))
        entries.append(LibraryEntry(f"compound {i}", mz, rng.integers(1, 1000, mz.size).astype(float),
                                    ri=float(rng.integers(600, 2200))))
    lib = compact_from_entries(entries)
    lib.ri_exp_sigma[:] = 10.0
    mz = np.unique(rng.integers(30, 250, size=30))
    inten = rng.integers(1, 1000, mz.size).astype(float)
    q = Query(mz=mz, intensity=inten, ri=float("nan"))
    got = match_queries([q], lib, MatchSettings(use_decoys=False, max_candidates=5), mz_range=(30, 249))[0]

    def reference(e):
        axis = np.arange(30, 250)
        u = np.zeros(axis.size); v = np.zeros(axis.size)
        u[mz - 30] = inten; v[e.mz - 30] = e.intensity
        wu = u ** 0.6 * axis ** 3.0; wv = v ** 0.6 * axis ** 3.0
        return 999.0 * (wu @ wv) ** 2 / ((wu @ wu) * (wv @ wv))

    scores = np.array([reference(e) for e in entries])
    top = np.argsort(-scores)[:5]
    assert [c.name for c in got.candidates] == [entries[i].name for i in top]
    assert np.allclose([c.mf for c in got.candidates], scores[top], rtol=1e-3)


def test_reverse_match_ignores_library_gaps(lib):
    """A query with extra ions from a neighbour lowers MF but not RMF."""
    row = lib.matrix.getrow(2)                                      # Toluene
    mz = np.concatenate([row.indices, [120, 150]]).astype(np.int64)
    inten = np.concatenate([row.data, [900.0, 900.0]]).astype(np.float64)
    m = match_queries([Query(mz=mz, intensity=inten, ri=768.0, ri_sigma=5.0)], lib, MatchSettings(use_decoys=False))[0]
    assert m.best.name == "Toluene" and m.best.rmf > m.best.mf + 100


def test_decoys_give_q_values_for_many_peaks_and_none_for_few(lib):
    queries = [query_from(lib, i % 6, ri=r) for i, r in zip(range(30), [1000, 1100, 768, 950, 1050, 880] * 5)]
    many = match_queries(queries, lib, MatchSettings())
    assert all(np.isfinite(m.q_value) for m in many if m.candidates)
    few = match_queries(queries[:5], lib, MatchSettings())
    assert all(np.isnan(m.q_value) for m in few)
    assert "HIGH_FDR" not in {f for m in few for f in m.flags}


def test_progress_and_cancel_are_called(lib):
    seen = []
    match_queries([query_from(lib, 0, ri=1000.0)] * 3, lib, MatchSettings(use_decoys=False),
                  progress=lambda f, t: seen.append(f))
    assert seen and seen[-1] == 1.0

    def stop():
        raise RuntimeError("stop")
    with pytest.raises(RuntimeError):
        match_queries([query_from(lib, 0)], lib, cancel=stop)


def test_settings_round_trip():
    s = MatchSettings(column_class="Standard polar", exclude_ions=(73, 207), mf_confident=900.0)
    assert MatchSettings.from_dict(s.to_dict()) == s
