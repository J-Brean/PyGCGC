"""⛓️ The n-alkane ladder, homologous series, stencils and markers."""

import numpy as np
import pytest

from pygcgc import (HydrocarbonSettings, PeakSettings, find_alkane_series, find_blobs,
                    find_homologous_series, find_markers, fold, propagate_alkane_series,
                    qualifier_ratios_at, ratios_match, remove_background, stencil_from_polygon,
                    stencil_totals)
from tests.conftest import make_run

PERIOD = 4.0
RATE = 50.0
CARBONS = list(range(10, 31))
T1_START, STEP = 200.0, 40.0


def _fold_pos(f, t1_abs, t2_abs):
    col, t2 = f.to_fold(t1_abs + t2_abs)
    return float(f.t1[col] + 0.5 * f.spec.period), float(t2)


def _build_run(shift_t1=0.0, shift_t2=0.0, branched=True, seed=0):
    peaks = []
    for i in range(len(CARBONS)):
        t1 = T1_START + i * STEP + shift_t1
        t2 = 1.0 + 0.01 * i + shift_t2
        peaks += [(t1, t2, 57, 2e5), (t1, t2, 71, 1.32e5), (t1, t2, 85, 8.8e4)]   # the measured n-alkane ratios
    if branched:
        for i in range(len(CARBONS)):
            t1 = T1_START + i * STEP + 5.0 + shift_t1
            t2 = 1.3 + 0.01 * i + shift_t2
            peaks += [(t1, t2, 43, 1.5e5), (t1, t2, 57, 3e4)]
    return make_run(peaks, phase=0.7, rate=RATE, mz=(40, 90), duration=1100.0, seed=seed)


def _seeded_series(run, f, carbons=(12, 20), **settings_kw):
    seeds = []
    for c in carbons:
        i = CARBONS.index(c)
        seeds.append((c, *_fold_pos(f, T1_START + i * STEP, 1.0 + 0.01 * i)))
    settings = HydrocarbonSettings(carbon_range=(8, 32), **settings_kw)
    return find_alkane_series(run, f, seeds, settings)


def test_find_alkane_series_from_two_seeds():
    run = _build_run()
    f = fold(run, PERIOD, 0.7)
    series = _seeded_series(run, f)
    assert len(series) == len(CARBONS)          # every carbon in range confirmed
    assert series.step_s() == pytest.approx(STEP, abs=0.1)
    # The true t2 drift (0.01 s per carbon) is recovered.
    t2_15, t2_25 = series.t2_of(15), series.t2_of(25)
    assert (t2_25 - t2_15) == pytest.approx(0.10, abs=0.01)


def test_carbon_and_t1_are_inverses():
    run = _build_run()
    f = fold(run, PERIOD, 0.7)
    series = _seeded_series(run, f)
    for c in (11.5, 18.0, 27.3):
        assert series.carbon_of(series.t1_of(c)) == pytest.approx(c, abs=1e-6)
    # Extrapolation beyond the confirmed range still follows the same line.
    assert series.t1_of(35) == pytest.approx(series.t1_of(30) + 5 * STEP, abs=0.5)


def test_needs_at_least_two_seeds():
    run = _build_run()
    f = fold(run, PERIOD, 0.7)
    with pytest.raises(ValueError):
        find_alkane_series(run, f, [(12, 300.0, 1.0)])


def test_propagate_alkane_series_follows_a_shifted_run():
    ref_run = _build_run()
    ref_f = fold(ref_run, PERIOD, 0.7)
    ref_series = _seeded_series(ref_run, ref_f)

    sample_run = _build_run(shift_t1=6.0, shift_t2=0.03, seed=5)
    sample_f = fold(sample_run, PERIOD, 0.7)
    sample_series = propagate_alkane_series(ref_series, sample_run, sample_f)
    assert len(sample_series) == len(ref_series)
    # Every carbon moved by roughly the true shift (one modulation of
    # quantisation noise either way, per the propagate_anchors convention).
    for c in ref_series.hits:
        assert sample_series.hits[c].t1 - ref_series.hits[c].t1 == pytest.approx(6.0, abs=PERIOD)


def test_stencil_lifts_correctly_onto_a_shifted_run():
    ref_run = _build_run()
    ref_f = fold(ref_run, PERIOD, 0.7)
    ref_series = _seeded_series(ref_run, ref_f)
    clean_ref = remove_background(ref_run, ref_f).run
    blobs_ref = find_blobs(clean_ref, ref_f, PeakSettings(snr_min=8), floor=1.0)

    # A box around the branched-alkane band, C15-C25, in the reference's
    # own alkane-relative coordinates.
    box = [ref_series.to_absolute(15, 0.15), ref_series.to_absolute(25, 0.15),
          ref_series.to_absolute(25, 0.45), ref_series.to_absolute(15, 0.45)]
    stencil = stencil_from_polygon(box, ref_series, "branched C15-25", response_factor=2.0)
    vol_ref, n_ref, _ = stencil_totals(stencil.to_absolute(ref_series), blobs_ref)
    assert n_ref == 10          # C15 through C24 inclusive (C25 sits exactly on the box edge)

    sample_run = _build_run(shift_t1=6.0, shift_t2=0.03, seed=7)
    sample_f = fold(sample_run, PERIOD, 0.7)
    sample_series = propagate_alkane_series(ref_series, sample_run, sample_f)
    clean_sample = remove_background(sample_run, sample_f).run
    blobs_sample = find_blobs(clean_sample, sample_f, PeakSettings(snr_min=8), floor=1.0)

    lifted = stencil.to_absolute(sample_series)
    vol_s, n_s, _ = stencil_totals(lifted, blobs_sample)
    assert n_s == n_ref                                          # the same compounds, despite the shift
    assert vol_s == pytest.approx(vol_ref, rel=0.05)
    assert stencil.response_factor == 2.0


def test_stencil_lift_finds_nothing_when_polygon_sits_outside_the_data():
    ref_run = _build_run()
    ref_f = fold(ref_run, PERIOD, 0.7)
    ref_series = _seeded_series(ref_run, ref_f)
    clean_ref = remove_background(ref_run, ref_f).run
    blobs_ref = find_blobs(clean_ref, ref_f, PeakSettings(snr_min=8), floor=1.0)
    box = [ref_series.to_absolute(15, 5.0), ref_series.to_absolute(25, 5.0),
          ref_series.to_absolute(25, 5.3), ref_series.to_absolute(15, 5.3)]
    vol, n, idx = stencil_totals(box, blobs_ref)
    assert n == 0 and vol == 0.0 and idx == []


def test_find_homologous_series_recovers_the_branched_alkanes():
    run = _build_run()
    f = fold(run, PERIOD, 0.7)
    series = _seeded_series(run, f)
    clean = remove_background(run, f).run
    blobs = find_blobs(clean, f, PeakSettings(snr_min=8), floor=1.0)
    hs = find_homologous_series(blobs, (43, 57), series, name="branched")
    assert hs is not None
    assert len(hs.members) == len(CARBONS)
    carbons_found = sorted(m["carbon"] for m in hs.members)
    assert carbons_found[0] == pytest.approx(10, abs=0.5)
    assert carbons_found[-1] == pytest.approx(30, abs=0.5)


def test_find_homologous_series_none_when_ion_pair_absent():
    run = _build_run(branched=False)
    f = fold(run, PERIOD, 0.7)
    series = _seeded_series(run, f)
    clean = remove_background(run, f).run
    blobs = find_blobs(clean, f, PeakSettings(snr_min=8), floor=1.0)
    assert find_homologous_series(blobs, (43, 57), series, name="branched") is None


def test_find_markers_returns_a_blob_near_each_expected_position():
    run = _build_run()
    f = fold(run, PERIOD, 0.7)
    series = _seeded_series(run, f)
    clean = remove_background(run, f).run
    blobs = find_blobs(clean, f, PeakSettings(snr_min=8), floor=1.0)
    markers = find_markers(blobs, series)
    assert set(markers) == {"pristane", "phytane"}
    assert all(0 <= i < len(blobs) for i in markers.values())


def test_settings_round_trip():
    s = HydrocarbonSettings(quant_ion=71, min_series_members=6)
    assert HydrocarbonSettings.from_dict(s.to_dict()) == s
    assert HydrocarbonSettings.from_dict({"min_height": 99, "unknown": 1}).min_height == 99
    t = HydrocarbonSettings(qual_ratios=(0.5, 0.3), enforce_qualifiers=False)
    assert HydrocarbonSettings.from_dict(t.to_dict()) == t


def test_default_qualifier_ratios_are_the_measured_ones():
    s = HydrocarbonSettings()
    assert s.qual_ions == (71, 85)
    assert s.qual_ratios == pytest.approx((0.66, 0.44), abs=0.03)


def test_ratios_match_uses_the_relative_tolerance():
    s = HydrocarbonSettings(qual_ratios=(0.6, 0.4), qual_tol=0.25)
    assert ratios_match((0.6, 0.4), s)
    assert ratios_match((0.7, 0.35), s)
    assert not ratios_match((0.9, 0.4), s)          # 50% high
    assert ratios_match((0.0, 9.0), HydrocarbonSettings(qual_ratios=()))   # nothing to check


def test_qualifier_ratios_at_a_clicked_peak_recover_the_planted_ratios():
    run = _build_run(branched=False)
    f = fold(run, PERIOD, 0.7)
    i = CARBONS.index(15)
    t1, t2 = _fold_pos(f, T1_START + i * STEP, 1.0 + 0.01 * i)
    series = _seeded_series(run, f)
    hit = series.hits[15]
    ratios = qualifier_ratios_at(run, f, hit.t1, hit.t2)
    assert ratios == pytest.approx((0.66, 0.44), abs=0.04)


def _run_with_impostors(offset_s=8.0, missing=()):
    """Every alkane has a taller neighbour on the same t2 whose spectrum is
    not an n-alkane's (57 and 43, no 71 or 85). Carbons in `missing` have
    only the impostor."""
    peaks = []
    for i, c in enumerate(CARBONS):
        t1, t2 = T1_START + i * STEP, 1.0 + 0.01 * i
        if c not in missing:
            peaks += [(t1, t2, 57, 2e5), (t1, t2, 71, 1.32e5), (t1, t2, 85, 8.8e4)]
        peaks += [(t1 + offset_s, t2, 57, 6e5), (t1 + offset_s, t2, 43, 6e5)]
    return make_run(peaks, phase=0.7, rate=RATE, mz=(40, 90), duration=1100.0, seed=3)


def test_the_qualifier_gate_ignores_a_taller_impostor_at_the_same_t2():
    run = _run_with_impostors()
    f = fold(run, PERIOD, 0.7)
    gated = _seeded_series(run, f)
    for c in (14, 18, 24):
        real_t1 = T1_START + CARBONS.index(c) * STEP
        assert gated.hits[c].t1 == pytest.approx(real_t1, abs=PERIOD)
        assert gated.hits[c].qual_ok


def test_a_carbon_with_only_an_impostor_is_rejected_unless_the_gate_is_switched_off():
    run = _run_with_impostors(missing=(18,))
    f = fold(run, PERIOD, 0.7)
    assert 18 not in _seeded_series(run, f).hits
    ungated = _seeded_series(run, f, enforce_qualifiers=False)
    assert 18 in ungated.hits and not ungated.hits[18].qual_ok        # reported, but flagged


def test_the_search_is_two_dimensional_a_bright_pulse_at_another_t2_is_ignored():
    peaks = []
    for i in range(len(CARBONS)):
        t1, t2 = T1_START + i * STEP, 1.0 + 0.01 * i
        peaks += [(t1, t2, 57, 2e5), (t1, t2, 71, 1.32e5), (t1, t2, 85, 8.8e4)]
        peaks += [(t1 + 3.0, t2 + 1.5, 57, 8e5), (t1 + 3.0, t2 + 1.5, 71, 5.3e5), (t1 + 3.0, t2 + 1.5, 85, 3.5e5)]
    run = make_run(peaks, phase=0.7, rate=RATE, mz=(40, 90), duration=1100.0, seed=4)
    f = fold(run, PERIOD, 0.7)
    series = _seeded_series(run, f)
    for c in (13, 19, 27):
        # the alkane line (about 0.9 s after the fold phase), not the 4x brighter pulse 1.5 s above it
        assert series.hits[c].t2 == pytest.approx(0.85 + 0.01 * CARBONS.index(c), abs=0.15)


class _FakeBlobs:
    """The few fields `find_homologous_series` reads."""
    def __init__(self, rows, mz0=40, n_mz=100):
        self.mz = np.arange(mz0, mz0 + n_mz)
        self.t1 = np.array([r["t1"] for r in rows])
        self.t2 = np.array([r["t2"] for r in rows])
        self.height = np.array([r.get("h", 1000.0) for r in rows])
        self.volume = np.full(len(rows), 1000.0)
        self.spectra = np.zeros((len(rows), n_mz))
        for k, r in enumerate(rows):
            self.spectra[k, 91 - mz0] = 400.0
            self.spectra[k, 105 - mz0] = 100.0

    def __len__(self):
        return len(self.t1)


def _ladder():
    from pygcgc import AlkaneHit, AlkaneSeries
    return AlkaneSeries(hits={c: AlkaneHit(carbon=c, t1=100.0 + (c - 10) * (60.0 - 2.0 * (c - 10)),
                                           t2=1.0 + 0.02 * (c - 10), height=1e4, qual_ratios=(0.66, 0.44))
                              for c in range(10, 20)})


def test_a_streak_of_blobs_at_one_t1_is_not_a_series_but_an_evenly_spaced_run_is():
    ladder = _ladder()
    # 12 blobs at a single t1 with t2 all over the place (a solvent-front streak) ...
    rows = [dict(t1=ladder.t1_of(12.3) + 0.2 * k, t2=1.0 + 0.6 * k) for k in range(12)]
    assert find_homologous_series(_FakeBlobs(rows), (91, 105), ladder) is None
    # ... and a real series: one member per carbon at 10.4, 11.4, ... on a steady t2 offset. The spacing in
    # seconds shrinks with carbon, which only a carbon-equivalent model (not a fixed step in seconds) survives.
    rows = [dict(t1=ladder.t1_of(c + 0.4), t2=ladder.t2_of(c + 0.4) + 0.4, h=1000.0 + c) for c in range(10, 18)]
    hs = find_homologous_series(_FakeBlobs(rows + rows[:2]), (91, 105), ladder)      # duplicates collapse to one per carbon
    assert hs is not None and len(hs.members) == 8
    assert [round(m["carbon"] % 1, 1) for m in hs.members] == [0.4] * 8


def test_a_series_needs_both_diagnostic_ions():
    ladder = _ladder()
    rows = [dict(t1=ladder.t1_of(c + 0.4), t2=ladder.t2_of(c + 0.4) + 0.4) for c in range(10, 18)]
    blobs = _FakeBlobs(rows)
    blobs.spectra[:, 105 - 40] = 0.0        # the minor ion is absent everywhere
    assert find_homologous_series(blobs, (91, 105), ladder) is None


def test_marker_offsets_are_measured_and_small():
    from pygcgc.core.hydrocarbons import MARKER_OFFSET_CARBONS
    assert set(MARKER_OFFSET_CARBONS) == {"pristane", "phytane"}
    assert all(abs(v) < 0.2 for v in MARKER_OFFSET_CARBONS.values())      # within a column of the alkane, not -0.3
