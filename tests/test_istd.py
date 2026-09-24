"""🏷️ Internal-standard library and cross-run matching."""

import numpy as np
import pytest

from pygcgc.core.istd import (IstdLibraryEntry, IstdMatchSettings, IstdReference, add_custom_entry,
                              build_report, find_istd_on_run, load_library, position_score,
                              save_library, write_istd_report_csv)
from pygcgc.core.peaks import BlobSet


def _blob_set(peaks):
    """A hand-built BlobSet carrying only what `find_istd_on_run` reads:
    mz, t1, t2, spectra/apex_spectra (see test_identify.py's own precedent)."""
    if not peaks:
        empty = np.zeros(0)
        return BlobSet(mz=np.arange(0, 1), t1=empty, t2=empty, col=empty.astype(int), row=empty.astype(int),
                      height=empty, snr=empty, volume=empty, n_pixels=empty, n_modulations=empty,
                      spectra=np.zeros((0, 1)), labels=np.zeros((1, 1), dtype=np.int32),
                      wrap_rows=0, pitch=0.02, period=4.0, t1_start=0.0, apex_spectra=np.zeros((0, 1)))
    mz_min = min(int(mz.min()) for _t1, _t2, mz, _i in peaks)
    mz_max = max(int(mz.max()) for _t1, _t2, mz, _i in peaks)
    mz_axis = np.arange(mz_min, mz_max + 1)
    n = len(peaks)
    apex = np.zeros((n, mz_axis.size))
    t1 = np.zeros(n)
    t2 = np.zeros(n)
    for i, (pt1, pt2, mz, intensity) in enumerate(peaks):
        apex[i, mz - mz_min] = intensity
        t1[i], t2[i] = pt1, pt2
    zeros = np.zeros(n)
    return BlobSet(mz=mz_axis, t1=t1, t2=t2, col=zeros.astype(int), row=zeros.astype(int),
                  height=zeros + 1.0, snr=zeros + 10.0, volume=apex.sum(axis=1), n_pixels=zeros + 10,
                  n_modulations=zeros + 2, spectra=apex.copy(), labels=np.zeros((1, 1), dtype=np.int32),
                  wrap_rows=0, pitch=0.02, period=4.0, t1_start=0.0, apex_spectra=apex)


# ------------------------------------------------------------------ library
def test_load_library_seeds_the_user_path_from_the_shipped_defaults(tmp_path):
    path = tmp_path / "istd_library.json"
    entries = load_library(path)
    assert path.exists()
    names = [e.name for e in entries]
    assert "Acenaphthene-d10" in names
    assert len(entries) == 7
    assert all(not e.custom for e in entries)
    entry = next(e for e in entries if e.name == "Acenaphthene-d10")
    assert entry.cas == "15067-26-2"
    assert entry.product_code == "1524.12"


def test_load_library_does_not_reseed_an_existing_user_library(tmp_path):
    path = tmp_path / "istd_library.json"
    entries = load_library(path)
    save_library(entries[:2], path)
    reloaded = load_library(path)
    assert len(reloaded) == 2


def test_add_custom_entry_round_trips(tmp_path):
    path = tmp_path / "istd_library.json"
    entries = load_library(path)
    updated = add_custom_entry(entries, "My-own-IS", cas="000-00-0", product_code="", path=path)
    assert any(e.name == "My-own-IS" and e.custom for e in updated)
    reloaded = load_library(path)
    assert any(e.name == "My-own-IS" and e.custom for e in reloaded)


def test_library_entry_to_dict_round_trip():
    entry = IstdLibraryEntry(name="Chrysene-d12", cas="1719-03-5", product_code="1024.18")
    assert IstdLibraryEntry.from_dict(entry.to_dict()) == entry


# ------------------------------------------------------------------ matching
def _reference(t1, t2, mz, intensity, name="Acenaphthene-d10", cas="15067-26-2", run="ref"):
    mz = np.asarray(mz)
    mz0 = int(mz.min())
    spectrum = np.zeros(int(mz.max()) - mz0 + 1)
    spectrum[mz - mz0] = intensity
    return IstdReference(name=name, cas=cas, run=run, t1=t1, t2=t2, mz0=mz0, spectrum=spectrum)


def test_find_istd_prefers_the_true_match_over_a_spectrally_similar_but_distant_decoy():
    ref_mz = [95, 96, 97, 98]
    ref_intensity = [100, 999, 300, 50]
    reference = _reference(500.0, 2.0, ref_mz, ref_intensity)
    peaks = _blob_set([
        (500.5, 2.05, np.array(ref_mz), np.array(ref_intensity)),   # true match: near position, same spectrum
        (503.0, 8.0, np.array(ref_mz), np.array(ref_intensity)),     # decoy: identical spectrum, far outside window
        (500.4, 1.95, np.array([60, 61, 62]), np.array([999, 200, 50])),  # decoy: close position, unrelated spectrum
    ])
    match = find_istd_on_run(reference, "sample", peaks)
    assert match.tier in ("high", "medium")
    np.testing.assert_allclose(match.t1, 500.5)
    np.testing.assert_allclose(match.t2, 2.05)


def test_find_istd_rejects_a_nearby_but_spectrally_dissimilar_decoy():
    ref_mz = [95, 96, 97, 98]
    ref_intensity = [100, 999, 300, 50]
    reference = _reference(500.0, 2.0, ref_mz, ref_intensity)
    peaks = _blob_set([
        (500.2, 2.02, np.array([60, 61, 62]), np.array([999, 200, 50])),
    ])
    match = find_istd_on_run(reference, "sample", peaks)
    assert match.spectral_score < 0.3
    assert match.tier == "low"


def test_find_istd_reports_not_found_when_nothing_is_in_the_position_window():
    reference = _reference(500.0, 2.0, [95, 96], [100, 999])
    settings = IstdMatchSettings(t1_window=5.0, t2_window=0.1)
    peaks = _blob_set([(600.0, 2.0, np.array([95, 96]), np.array([100, 999]))])
    match = find_istd_on_run(reference, "sample", peaks, settings)
    assert match.tier == "not found"
    assert match.peak_id is None


def test_find_istd_reports_not_found_on_an_empty_run():
    reference = _reference(500.0, 2.0, [95, 96], [100, 999])
    peaks = _blob_set([])
    match = find_istd_on_run(reference, "sample", peaks)
    assert match.tier == "not found"


# ------------------------------------------------------------------ scoring / tiers
def test_position_score_decays_with_distance():
    settings = IstdMatchSettings()
    assert position_score(0.0, 0.0, settings) == pytest.approx(1.0)
    assert position_score(60.0, 0.0, settings) < position_score(10.0, 0.0, settings)


@pytest.mark.parametrize("dt1,expected", [(0.0, "high"), (10.0, "high"), (30.0, "medium"), (70.0, "low")])
def test_tier_thresholds_follow_the_combined_score(dt1, expected):
    """With spectral_weight=0, the combined score is just the position score,
    so its exponential decay against `dt1` gives predictable tier boundaries
    (tier_high=0.75 at dt1~17.3s, tier_medium=0.45 at dt1~48.0s, for the
    default t1_scale=60)."""
    reference = _reference(0.0, 0.0, [95, 96], [100, 999])
    peaks = _blob_set([(dt1, 0.0, np.array([95, 96]), np.array([100, 999]))])
    settings = IstdMatchSettings(position_weight=1.0, spectral_weight=0.0)
    match = find_istd_on_run(reference, "sample", peaks, settings)
    assert match.tier == expected


# ------------------------------------------------------------------ drift-aware matching
_SPECTRA = {"A": ([95, 96, 97], [100, 999, 300]), "B": ([60, 61, 62], [999, 200, 50]),
            "C": ([120, 121, 122], [300, 999, 100]), "D": ([80, 82, 84], [999, 500, 200])}
_REF_POS = {"A": (500.0, 2.0), "B": (900.0, 1.5), "C": (1300.0, 3.0), "D": (1700.0, 2.5)}


def _references():
    return [_reference(*_REF_POS[k], *_SPECTRA[k], name=k) for k in _REF_POS]


def _shifted_run(dt1, dt2):
    return _blob_set([(t1 + dt1, t2 + dt2, np.array(_SPECTRA[k][0]), np.array(_SPECTRA[k][1]))
                      for k, (t1, t2) in _REF_POS.items()])


def test_match_run_scores_position_against_the_other_standards_drift():
    """A whole run shifted by 45 s and 0.7 s: every IS still matches, and with
    the shift predicted from the other standards each is high confidence,
    where scoring against the raw reference position would not be."""
    from pygcgc.core.istd import match_run

    peaks = _shifted_run(45.0, 0.7)
    matches = {m.name: m for m in match_run(_references(), "shifted", peaks)}
    for name, m in matches.items():
        assert m.tier == "high", (name, m)
        assert m.dt1 == pytest.approx(45.0)
        assert m.n_drift_anchors == 3
        assert m.expected_t1 == pytest.approx(_REF_POS[name][0] + 45.0, abs=1e-6)
    raw = find_istd_on_run(_references()[0], "shifted", peaks)
    assert raw.combined_score < matches["A"].combined_score


def test_match_run_does_not_let_a_standard_vouch_for_its_own_position():
    """With only one confident standard there are no *other* anchors, so it
    falls back to its reference position rather than predicting itself."""
    from pygcgc.core.istd import match_run

    peaks = _shifted_run(45.0, 0.7)
    (m,) = match_run(_references()[:1], "shifted", peaks)
    assert m.n_drift_anchors == 0
    assert m.expected_t1 == pytest.approx(_REF_POS["A"][0])


# ------------------------------------------------------------------ report
def test_build_report_and_csv_round_trip(tmp_path):
    reference = _reference(500.0, 2.0, [95, 96], [100, 999])
    peaks = _blob_set([(500.1, 2.0, np.array([95, 96]), np.array([100, 999]))])
    match = find_istd_on_run(reference, "sample", peaks)
    rows = build_report([match])
    assert rows[0]["run"] == "sample"
    assert rows[0]["name"] == "Acenaphthene-d10"
    path = write_istd_report_csv([match], tmp_path / "report.csv")
    text = path.read_text(encoding="utf-8")
    assert "sample" in text and "Acenaphthene-d10" in text
