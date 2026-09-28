"""🏷️ Internal-standard library and cross-run matching."""

import numpy as np
import pytest

from pygcgc.core.istd import (FOUND, NO_PEAKS, NOT_FOUND, IstdLibraryEntry, IstdMatch, IstdMatchSettings,
                              IstdReference, add_custom_entry, build_report, find_istd_on_run, load_library,
                              match_run, references_to_targets, save_library, write_istd_report_csv)
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
    mz_min = min(int(np.min(mz)) for _t1, _t2, mz, _i in peaks)
    mz_max = max(int(np.max(mz)) for _t1, _t2, mz, _i in peaks)
    mz_axis = np.arange(mz_min, mz_max + 1)
    n = len(peaks)
    apex = np.zeros((n, mz_axis.size))
    t1 = np.zeros(n)
    t2 = np.zeros(n)
    for i, (pt1, pt2, mz, intensity) in enumerate(peaks):
        apex[i, np.asarray(mz) - mz_min] = intensity
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


# ------------------------------------------------------------------ the three thresholds
REF_MZ, REF_I = [95, 96, 97, 98], [100, 999, 300, 50]
OTHER_MZ, OTHER_I = [60, 61, 62], [999, 200, 50]


def _reference(t1, t2, mz, intensity, name="Acenaphthene-d10", cas="15067-26-2", run="ref"):
    mz = np.asarray(mz)
    mz0 = int(mz.min())
    spectrum = np.zeros(int(mz.max()) - mz0 + 1)
    spectrum[mz - mz0] = intensity
    return IstdReference(name=name, cas=cas, run=run, t1=t1, t2=t2, mz0=mz0, spectrum=spectrum)


def test_a_peak_passing_all_three_thresholds_is_found():
    reference = _reference(500.0, 2.0, REF_MZ, REF_I)
    peaks = _blob_set([(505.0, 2.1, REF_MZ, REF_I)])
    m = find_istd_on_run(reference, "sample", peaks)
    assert m.status == FOUND and m.found
    assert (m.t1, m.t2) == pytest.approx((505.0, 2.1))
    assert (m.dt1, m.dt2) == pytest.approx((5.0, 0.1))
    assert m.spectral_score > 0.99
    assert m.n_passing == 1 and m.reason == ""


@pytest.mark.parametrize("t1,t2", [(531.0, 2.0), (500.0, 2.6)], ids=["t1", "t2"])
def test_a_peak_just_outside_either_position_threshold_is_not_found(t1, t2):
    """First and second dimension are each tested on their own: 31 s is outside
    max_dt1 = 30 s even at the right t2, and 0.6 s outside max_dt2 = 0.5 s even
    at the right t1, however good the spectrum."""
    reference = _reference(500.0, 2.0, REF_MZ, REF_I)
    m = find_istd_on_run(reference, "sample", _blob_set([(t1, t2, REF_MZ, REF_I)]),
                         IstdMatchSettings(max_dt1=30.0, max_dt2=0.5, follow_drift=False))
    assert m.status == NOT_FOUND
    assert "no peak within 30 s in t1 and 0.5 s in t2" in m.reason


def test_a_peak_in_position_with_the_wrong_spectrum_is_not_found_and_says_why():
    reference = _reference(500.0, 2.0, REF_MZ, REF_I)
    m = find_istd_on_run(reference, "sample", _blob_set([(500.2, 2.02, OTHER_MZ, OTHER_I)]))
    assert m.status == NOT_FOUND
    assert m.spectral_score < 0.3
    assert (m.t1, m.t2) == pytest.approx((500.2, 2.02))     # the closest miss is reported
    assert "below 0.80" in m.reason


def test_min_spectral_decides_a_borderline_spectrum():
    """A spectrum that matches partly passes a lenient threshold and fails a strict one."""
    reference = _reference(500.0, 2.0, REF_MZ, REF_I)
    partial = _blob_set([(500.0, 2.0, REF_MZ, [999, 300, 100, 600])])    # same ions, other ratios
    score = find_istd_on_run(reference, "s", partial, IstdMatchSettings(min_spectral=0.0)).spectral_score
    assert 0.2 < score < 0.98
    assert find_istd_on_run(reference, "s", partial, IstdMatchSettings(min_spectral=score - 0.01)).found
    assert not find_istd_on_run(reference, "s", partial, IstdMatchSettings(min_spectral=score + 0.01)).found


def test_of_the_peaks_that_pass_the_best_spectrum_is_taken_and_counted():
    reference = _reference(500.0, 2.0, REF_MZ, REF_I)
    peaks = _blob_set([(500.5, 2.0, REF_MZ, [100, 999, 300, 150]),     # passes, but a slightly worse match
                       (510.0, 2.2, REF_MZ, REF_I),                      # passes, identical spectrum
                       (500.1, 2.0, OTHER_MZ, OTHER_I)])                 # nearest, but fails the spectrum
    m = find_istd_on_run(reference, "sample", peaks, IstdMatchSettings(min_spectral=0.8))
    assert m.found and m.t1 == pytest.approx(510.0)
    assert m.n_passing == 2 and "2 peaks pass" in m.reason


def test_an_empty_run_reports_no_peaks():
    reference = _reference(500.0, 2.0, REF_MZ, REF_I)
    m = find_istd_on_run(reference, "sample", _blob_set([]))
    assert m.status == NO_PEAKS and not m.found


# ------------------------------------------------------------------ drift
_SPECTRA = {"A": ([95, 96, 97], [100, 999, 300]), "B": ([60, 61, 62], [999, 200, 50]),
            "C": ([120, 121, 122], [300, 999, 100]), "D": ([80, 82, 84], [999, 500, 200])}
_REF_POS = {"A": (500.0, 2.0), "B": (900.0, 1.5), "C": (1300.0, 3.0), "D": (1700.0, 2.5)}


def _references():
    return [_reference(*_REF_POS[k], *_SPECTRA[k], name=k) for k in _REF_POS]


def _shifted_run(dt1, dt2):
    return _blob_set([(t1 + dt1, t2 + dt2, _SPECTRA[k][0], _SPECTRA[k][1]) for k, (t1, t2) in _REF_POS.items()])


def test_following_drift_finds_every_standard_on_a_shifted_run():
    """A whole run 45 s and 0.7 s later than the labelled one is outside tight
    position thresholds, but the other standards moved by the same amount, so
    with drift followed every standard is expected in the right place."""
    peaks = _shifted_run(45.0, 0.7)
    tight = IstdMatchSettings(max_dt1=5.0, max_dt2=0.2)
    matches = {m.name: m for m in match_run(_references(), "shifted", peaks, tight)}
    for name, m in matches.items():
        assert m.found, (name, m)
        assert m.dt1 == pytest.approx(45.0) and m.off_t1 == pytest.approx(0.0, abs=1e-6)
        assert m.n_drift_anchors == 3
    fixed = match_run(_references(), "shifted", peaks, IstdMatchSettings(max_dt1=5.0, max_dt2=0.2,
                                                                          follow_drift=False))
    assert not any(m.found for m in fixed)


def test_a_standard_does_not_vouch_for_its_own_position():
    """With only one standard there are no *other* anchors, so it is expected at
    its labelled position, not wherever it happens to be found."""
    peaks = _shifted_run(45.0, 0.7)
    (m,) = match_run(_references()[:1], "shifted", peaks)
    assert m.n_drift_anchors == 0
    assert m.expected_t1 == pytest.approx(_REF_POS["A"][0])
    assert m.status == NOT_FOUND                           # 0.7 s in t2 is outside max_dt2 = 0.5 s


# ------------------------------------------------------------------ report and saving
def test_report_rows_carry_the_thresholds_and_the_reason(tmp_path):
    reference = _reference(500.0, 2.0, REF_MZ, REF_I)
    found = find_istd_on_run(reference, "run1", _blob_set([(500.1, 2.0, REF_MZ, REF_I)]))
    missed = find_istd_on_run(reference, "run2", _blob_set([(500.1, 2.0, OTHER_MZ, OTHER_I)]))
    settings = IstdMatchSettings()
    rows = build_report([found, missed], settings)
    assert [r["status"] for r in rows] == [FOUND, NOT_FOUND]
    assert rows[0]["max_dt1_s"] == settings.max_dt1 and rows[0]["min_spectral"] == settings.min_spectral
    assert "below" in rows[1]["reason"]
    text = write_istd_report_csv([found, missed], tmp_path / "report.csv", settings).read_text(encoding="utf-8")
    assert "run1" in text and "run2" in text and "spectral_match" in text


def test_a_match_round_trips_through_a_dict_and_ignores_old_keys():
    from dataclasses import asdict
    reference = _reference(500.0, 2.0, REF_MZ, REF_I)
    m = find_istd_on_run(reference, "run1", _blob_set([(500.1, 2.0, REF_MZ, REF_I)]))
    assert IstdMatch.from_dict(asdict(m)) == m
    assert IstdMatch.from_dict({**asdict(m), "tier": "high", "combined_score": 0.9}) == m


# ------------------------------------------------------------------ to the target table
def test_labelled_standards_become_istd_rows_of_a_target_table(tmp_path):
    from pygcgc.core.targets import load_targets, write_targets_csv
    # phenanthrene-d10's measured spectrum on the real runs: 188 base, then 189, 184, 187 (and 80 below 100)
    ref = _reference(2708.2, 4.60, [80, 94, 184, 187, 188, 189], [30, 20, 18, 15, 100, 16],
                     name="phenanthrene-d10")
    (t,) = references_to_targets([ref])
    assert t.role == "istd" and t.quant_ion == 188
    assert t.qual_ions == (184, 189) or t.qual_ions == (189, 184)
    assert t.qual_ratios == pytest.approx(tuple({184: 0.18, 189: 0.16}[m] for m in t.qual_ions))
    assert (t.t1, t.t2) == (2708.2, 4.6)
    # The rows round-trip through the target-table file and become phase anchors.
    path = write_targets_csv([t], tmp_path / "target_table.csv")
    assert load_targets(path)[0].quant_ion == 188
