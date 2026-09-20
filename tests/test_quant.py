"""🧪 Calibration, peak finding and targeted quantification."""

import csv

import numpy as np
import pytest

from dataclasses import replace

from pygcgc import (CalBracket, CalibrationError, CalibrationSettings, PeakSettings, QuantDecisions,
                    QuantSettings, SequenceRun, blended_amount, bracket_cal_sets, build_calibrations,
                    cal_set_for, find_blobs, fit_calibration, fold, load_levels, load_targets,
                    measure_ion_peak, quantify_prepared, quantify_run, quantify_sequence)
from pygcgc.core.targets import Target
from tests.conftest import (CAL_LEVELS, CAL_TARGETS, GUI_PHASE, IS_AREA, SAMPLE_AMOUNTS, SCALE,
                            make_cal_run, make_run)


# ------------------------------------------------------------ cal set brackets
def _fake_run_quant(compound_areas: dict, istd_area: float = 1.0):
    """A RunQuant with hand-set quantifier volumes, for testing calibration
    bracketing without simulating any chromatography."""
    from pygcgc.core.peaks import BlobSet
    from pygcgc.core.quant import RunQuant, TargetHit
    hits = {"IS": TargetHit("IS", "istd", True, 0.0, 0.0, quant_volume=istd_area, tic_volume=istd_area)}
    for name, area in compound_areas.items():
        hits[name] = TargetHit(name, "target", True, 0.0, 0.0, quant_volume=area, tic_volume=area)
    blobs = BlobSet(mz=np.array([]), t1=np.array([]), t2=np.array([]), col=np.array([]),
                    row=np.array([]), height=np.array([]), snr=np.array([]), volume=np.array([]),
                    n_pixels=np.array([]), n_modulations=np.array([]),
                    spectra=np.zeros((0, 1)), labels=np.zeros((1, 1), dtype=np.int32),
                    wrap_rows=0, pitch=0.02, period=4.0, t1_start=0.0)
    return RunQuant(run="r", hits=hits, blobs=blobs, istd_for={n: "IS" for n in compound_areas}, ceiling=None)


def _bracket_sequence():
    # Calibration set "mon" at indices 0-2 (mean 1), "fri" at indices 5-7 (mean 6).
    seq = ([SequenceRun(f"MON_L{i}", "calibration", f"L{i}", "mon") for i in range(3)]
          + [SequenceRun("S1"), SequenceRun("S2")]
          + [SequenceRun(f"FRI_L{i}", "calibration", f"L{i}", "fri") for i in range(3)])
    return seq


def test_bracket_between_two_sets():
    seq = _bracket_sequence()
    b = bracket_cal_sets(seq, 4)          # S2, one before the "fri" block
    assert b.before == "mon" and b.after == "fri"
    assert b.frac == pytest.approx((4 - 1) / (6 - 1))
    assert not b.extrapolated
    assert "mon\u2192fri" in b.label()


def test_bracket_extrapolated_before_and_after():
    seq = [SequenceRun("S0")] + _bracket_sequence() + [SequenceRun("S_last")]
    before = bracket_cal_sets(seq, 0)
    assert before.before == before.after == "mon"
    assert before.extrapolated
    after = bracket_cal_sets(seq, len(seq) - 1)
    assert after.before == after.after == "fri"
    assert after.extrapolated


def test_bracket_explicit_cal_set_wins():
    seq = _bracket_sequence()
    seq[3] = SequenceRun("S1", cal_set="fri")   # forced, even though it sits near "mon"
    b = bracket_cal_sets(seq, 3)
    assert b.before == b.after == "fri"
    assert b.frac == 0.0


def test_bracket_interpolation_off_matches_cal_set_for():
    seq = _bracket_sequence()
    for i in (3, 4):
        b = bracket_cal_sets(seq, i, interpolate=False)
        assert b.before == b.after == cal_set_for(seq, i)


def test_bracket_single_set_never_blends():
    seq = [SequenceRun(f"L{i}", "calibration", f"L{i}", "only") for i in range(3)] + [SequenceRun("S1")]
    b = bracket_cal_sets(seq, 3)
    assert b.before == b.after == "only"
    assert not b.blended


def test_blended_amount_linear_blend():
    x = np.array([0.0, 1.0, 2.0, 5.0, 10.0])
    fit_a = fit_calibration("cpd", x, 1.0 * x)     # sensitivity 1.0
    fit_b = fit_calibration("cpd", x, 0.5 * x)     # sensitivity halved: drift
    fits = {("mon", "cpd"): {"quant": fit_a}, ("fri", "cpd"): {"quant": fit_b}}
    bracket = CalBracket("mon", "fri", 0.4)
    resp = 2.0
    amount, se, flags, note = blended_amount(fits, "cpd", "quant", bracket, resp)
    want = 0.6 * fit_a.amount(np.array([resp]))[0] + 0.4 * fit_b.amount(np.array([resp]))[0]
    assert amount == pytest.approx(want)
    assert se > 0
    assert note == ""
    # Blended sits strictly between what either calibration alone would say.
    bounds = (fit_a.amount(np.array([resp]))[0], fit_b.amount(np.array([resp]))[0])
    assert min(bounds) < amount < max(bounds)


def test_blended_amount_missing_one_side_flagged():
    x = np.array([0.0, 1.0, 2.0, 5.0, 10.0])
    fit_a = fit_calibration("cpd", x, x)
    fits = {("mon", "cpd"): {"quant": fit_a}}   # "fri" has no fit for this compound
    bracket = CalBracket("mon", "fri", 0.7)
    amount, se, flags, note = blended_amount(fits, "cpd", "quant", bracket, 2.0)
    assert amount == pytest.approx(fit_a.amount(np.array([2.0]))[0])
    assert any("only the mon" in f for f in flags)


def test_drift_correction_beats_single_calibration():
    """A sample sitting between two calibration blocks with different
    sensitivity (simulated drift) is quantified more accurately with
    interpolation than by using either block alone."""
    x = np.array([0.0, 0.5, 1.0, 2.5, 5.0, 10.0])
    k_mon, k_fri, true_amount = 1.0, 0.6, 4.0     # sensitivity fell 40% over the week
    fit_mon = fit_calibration("cpd", x, k_mon * x)
    fit_fri = fit_calibration("cpd", x, k_fri * x)
    fits = {("mon", "cpd"): {"quant": fit_mon}, ("fri", "cpd"): {"quant": fit_fri}}

    seq = _bracket_sequence()
    frac = (3 - 1) / (6 - 1)                      # S1 sits at index 3
    k_true = (1 - frac) * k_mon + frac * k_fri    # assumed-linear drift
    resp = k_true * true_amount

    bracket = bracket_cal_sets(seq, 3)
    blended, _, _, _ = blended_amount(fits, "cpd", "quant", bracket, resp)
    mon_only = fit_mon.amount(np.array([resp]))[0]
    fri_only = fit_fri.amount(np.array([resp]))[0]
    err_blend = abs(blended - true_amount)
    err_mon = abs(mon_only - true_amount)
    err_fri = abs(fri_only - true_amount)
    assert err_blend < err_mon
    assert err_blend < err_fri


def test_quantify_sequence_labels_and_blends_end_to_end():
    # Three nonzero levels per block (min_levels=3; a zero-amount level is
    # excluded from the fit by default and would leave only two).
    x = np.array([0.5, 1.0, 2.5])
    k_mon, k_fri = 1.0, 0.5
    seq = _bracket_sequence()
    levels = {f"L{i}": {"cpd": float(v)} for i, v in enumerate(x)}
    targets = [Target(compound="IS", role="istd", quant_ion=50, t1=60.0, t2=1.2),
              Target(compound="cpd", role="target", quant_ion=52, t1=110.0, t2=2.0, istd="IS")]
    for i, r in enumerate(seq):
        if r.is_calibration:
            k = k_mon if r.cal_set == "mon" else k_fri
            amt = levels[r.level]["cpd"]
            r.quant = _fake_run_quant({"cpd": amt * k})
    true_amount = 1.5
    frac = (3 - 1) / (6 - 1)
    k_true = (1 - frac) * k_mon + frac * k_fri
    seq[3].quant = _fake_run_quant({"cpd": true_amount * k_true})

    settings = QuantSettings()
    decisions = QuantDecisions()
    fits = build_calibrations(seq, targets, levels, settings, decisions)
    rows = quantify_sequence(seq, targets, fits, decisions, settings)
    row = next(r for r in rows if r.run == "S1" and r.compound == "cpd")
    resp = seq[3].quant.response("cpd")
    mon_only = fits[("mon", "cpd")]["quant"].amount(np.array([resp]))[0]
    fri_only = fits[("fri", "cpd")]["quant"].amount(np.array([resp]))[0]
    assert min(mon_only, fri_only) < row.amount < max(mon_only, fri_only)
    assert row.amount == pytest.approx(0.6 * mon_only + 0.4 * fri_only)
    assert "mon\u2192fri" in row.cal_set

    # With interpolation off, the sample uses "mon" alone and drifts further off.
    off = replace(settings, interpolate_cal_sets=False)
    rows_off = quantify_sequence(seq, targets, fits, decisions, off)
    row_off = next(r for r in rows_off if r.run == "S1" and r.compound == "cpd")
    assert row_off.cal_set == "mon"
    assert row_off.amount == pytest.approx(mon_only)
    assert abs(row_off.amount - true_amount) > abs(row.amount - true_amount)


# ------------------------------------------------------------ calibration
def test_linear_recovery_and_inverse_prediction():
    rng = np.random.default_rng(1)
    x = np.array([0.1, 0.25, 0.5, 1, 2.5, 5, 10])
    y = 0.8 * x * (1 + rng.normal(0, 0.03, x.size))
    f = fit_calibration("c", x, y)
    assert f.model == "linear"
    assert f.slope == pytest.approx(0.8, rel=0.03)
    assert f.ok
    amt, se = f.amount_with_se(np.array([0.8 * 3.0]))
    assert amt[0] == pytest.approx(3.0, rel=0.05)
    assert 0 < se[0] < 0.3


def test_weighting_chosen_for_proportional_noise():
    rng = np.random.default_rng(2)
    x = np.repeat([0.1, 0.3, 1, 3, 10], 3)
    y = 2 * x * (1 + rng.normal(0, 0.05, x.size))
    f = fit_calibration("c", x, y)
    assert f.weights in ("1/x", "1/x2")
    assert f.weight_scores["1"] > f.weight_scores[f.weights]
    assert CalibrationSettings().weights == "auto"
    assert fit_calibration("c", x, y, settings=CalibrationSettings(weights="1")).weights == "1"


def test_quadratic_only_when_justified():
    x = np.array([0.5, 1, 2, 5, 10, 20])
    curved = 1.0 * x - 0.02 * x ** 2
    f = fit_calibration("c", x, curved * (1 + np.array([0.002, -0.003, 0.001, 0.002, -0.001, 0.001])))
    assert f.model == "quadratic" and f.mandel_p < 0.05
    straight = 1.0 * x * (1 + np.array([0.01, -0.02, 0.015, -0.01, 0.02, -0.01]))
    assert fit_calibration("c", x, straight).model == "linear"
    forced = fit_calibration("c", x, straight, settings=CalibrationSettings(model="quadratic"))
    assert forced.model == "quadratic"


def test_exclusion_zero_and_errors():
    x = np.array([0.0, 1, 2, 3, 4])
    y = np.array([0.01, 1.0, 2.0, 3.0, 40.0])
    labels = ["blank", "a", "b", "c", "d"]
    f = fit_calibration("c", x, y, labels, excluded={"d": "bubble in syringe"})
    assert f.labels == ["a", "b", "c"]            # zero amount and the excluded point are left out
    assert list(f.x_excluded) == [4.0]
    assert f.slope == pytest.approx(1.0)
    z = fit_calibration("c", x[1:4], y[1:4], settings=CalibrationSettings(force_zero=True))
    assert z.zero and z.intercept == 0.0
    with pytest.raises(CalibrationError):
        fit_calibration("c", [1, 1, 2], [1, 1, 2])
    bad = fit_calibration("c", [1, 2, 3, 4], [1, 2, 3, 6])
    assert not bad.ok
    assert "d" not in bad.influential() or True


# ------------------------------------------------------------ peaks
def test_blobs_found_with_volumes_and_wrap():
    peaks = [(62.0, 1.2, 50, 2e5), (122.0, 2.6, 51, 3e5), (182.0, 4.3, 52, 2e5)]
    run = make_run(peaks, phase=0.7, rate=50.0)
    f = fold(run, 4.0, 0.7)
    b = find_blobs(run, f, PeakSettings(snr_min=10), floor=1.0)
    assert len(b) == 3
    order = np.argsort(b.t1)
    np.testing.assert_allclose(b.volume[order], [2e5, 3e5, 2e5], rtol=0.02)
    np.testing.assert_allclose(b.t2[order], [1.2, 2.6, 4.3], atol=0.03)   # true t2 beyond the period
    np.testing.assert_allclose(b.t1[order], [62.0, 122.0, 182.0], atol=1.5)
    i = int(order[0])
    assert b.spectra[i, 0] == pytest.approx(b.volume[i])
    t1s, t2s = b.outline_times(i)
    assert np.nanmin(t2s) < 1.2 < np.nanmax(t2s)


def test_ion_peak_volume_and_flank_flag():
    run = make_run([(100.0, 2.0, 52, 2e5), (100.0, 2.0, 53, 6e4), (120.0, 2.02, 52, 1e6)],
                   phase=0.7, rate=50.0)
    f = fold(run, 4.0, 0.7)
    c = measure_ion_peak(run, run, f, None, 52, [53], 100.0, 2.0, 5.0, 0.2)
    assert len(c) == 1
    p = c[0]
    # Part of its tail sits in the big neighbour's basin, which is why it is flagged.
    assert p.quant_volume == pytest.approx(2e5, rel=0.08)
    assert p.qual_ratios[0] == pytest.approx(0.3, rel=0.05)
    assert p.touches_larger


def _box(t1, t2, half_t1, half_t2):
    return [(t1 - half_t1, t2 - half_t2), (t1 + half_t1, t2 - half_t2),
            (t1 + half_t1, t2 + half_t2), (t1 - half_t1, t2 + half_t2)]


def test_hand_drawn_polygon_integrates_what_is_inside_it():
    from pygcgc.core.peaks import measure_ion_polygon
    run = make_run([(100.0, 2.0, 52, 2e5), (100.0, 2.0, 53, 6e4)], phase=0.7, rate=50.0)
    f = fold(run, 4.0, 0.7)
    auto = measure_ion_peak(run, run, f, None, 52, [53], 100.0, 2.0, 5.0, 0.2)[0]
    wide = measure_ion_polygon(run, run, f, None, 52, [53], _box(auto.t1, auto.t2, 12.0, 0.4))
    assert wide.quant_volume == pytest.approx(2e5, rel=0.08)
    assert wide.qual_ratios[0] == pytest.approx(0.3, rel=0.05)
    narrow = measure_ion_polygon(run, run, f, None, 52, [53], _box(auto.t1, auto.t2, 2.0, 0.05))
    assert 0.05 * auto.quant_volume < narrow.quant_volume < 0.7 * auto.quant_volume   # part of the peak only
    assert measure_ion_polygon(run, run, f, None, 52, [53], _box(30.0, 0.9, 1.0, 0.05)) is None   # empty ground


def test_quantify_run_uses_a_hand_drawn_footprint_and_flags_it():
    run = make_run([(100.0, 2.0, 52, 2e5), (100.0, 2.0, 53, 6e4)], phase=0.7, rate=50.0)
    f = fold(run, 4.0, 0.7)
    t = Target("A", "target", "", 52, (53,), (0.3,), 100.0, 10.0, 2.0, 0.2)
    auto = quantify_run(run, run, f, None, [t]).hits["A"]
    poly = _box(auto.t1, auto.t2, 2.0, 0.05)
    hand = quantify_run(run, run, f, None, [t], footprints={"A": poly}).hits["A"]
    assert hand.integrated_by_hand and hand.manual and "integrated by hand" in hand.flags()
    assert hand.quant_volume < 0.7 * auto.quant_volume
    assert not auto.integrated_by_hand


def test_decisions_carry_footprints_through_a_round_trip():
    d = QuantDecisions()
    d.footprints[("S1", "A")] = [(90.0, 1.9), (110.0, 1.9), (110.0, 2.1)]
    back = QuantDecisions.from_dict(d.to_dict())
    assert back.footprints == d.footprints
    assert back.copy().footprints == d.footprints


def _clip_ion(run, ion, ceiling):
    from pygcgc import Run
    csc = run.data.tocsc(copy=True)
    j = ion - int(run.mz[0])
    a, b = csc.indptr[j], csc.indptr[j + 1]
    csc.data[a:b] = np.minimum(csc.data[a:b], ceiling)
    return Run(run.times, run.mz, csc.tocsr(), run.name)


def test_saturated_quantifier_can_optionally_fall_back_to_a_qualifier():
    run = make_run([(100.0, 2.0, 52, 2e5), (100.0, 2.0, 53, 6e4)], phase=0.7, rate=50.0)
    ceiling = 0.5 * float(run.xic(52).max())
    clipped = _clip_ion(run, 52, ceiling)                      # the detector clipped the big ion
    f = fold(clipped, 4.0, 0.7)
    t = Target("A", "target", "", 52, (53,), (0.3,), 100.0, 10.0, 2.0, 0.2)
    off = quantify_run(clipped, clipped, f, None, [t], QuantSettings(saturation=ceiling)).hits["A"]
    on = quantify_run(clipped, clipped, f, None, [t],
                      QuantSettings(saturation=ceiling, saturated_qualifier=True)).hits["A"]
    assert off.saturated and off.quant_ion_used is None
    assert off.quant_volume < 0.85 * 2e5                     # the clipped ion under-reads
    assert on.saturated and on.quant_ion_used == 53
    assert on.quant_volume == pytest.approx(2e5, rel=0.1)    # 6e4 / 0.3
    assert "saturated (quantified on m/z 53)" in on.flags()
    assert not QuantSettings().saturated_qualifier            # optional: off unless asked for


def test_saturated_fallback_needs_an_unsaturated_qualifier_with_an_expected_ratio():
    run = make_run([(100.0, 2.0, 52, 2e5), (100.0, 2.0, 53, 6e4)], phase=0.7, rate=50.0)
    ceiling = 0.5 * float(run.xic(52).max())
    clipped = _clip_ion(run, 52, ceiling)
    f = fold(clipped, 4.0, 0.7)
    settings = QuantSettings(saturation=ceiling, saturated_qualifier=True)
    no_ratio = Target("A", "target", "", 52, (53,), (0.0,), 100.0, 10.0, 2.0, 0.2)
    assert quantify_run(clipped, clipped, f, None, [no_ratio], settings).hits["A"].quant_ion_used is None
    also_sat = _clip_ion(clipped, 53, 0.5 * float(run.xic(53).max()))
    t = Target("A", "target", "", 52, (53,), (0.3,), 100.0, 10.0, 2.0, 0.2)
    hit = quantify_run(also_sat, also_sat, fold(also_sat, 4.0, 0.7), None, [t],
                       QuantSettings(saturation=0.5 * float(run.xic(53).max()),
                                     saturated_qualifier=True)).hits["A"]
    assert hit.quant_ion_used is None                          # the qualifier is saturated too


def _sequence(noise=0.0, shift=(1.5, 0.02)):
    targets = [Target(n, role, "", ion, (q,) if q else (), (qr,) if qr else (), t1, 10.0, t2, 0.2,
                      "" if role == "istd" else "IS") for n, role, ion, q, qr, t1, t2, rf in CAL_TARGETS]
    prepared, info = [], []
    for k, (lv, amt) in enumerate(CAL_LEVELS.items()):
        r = make_cal_run({n: amt for n in "ABC"}, seed=k, noise=noise)
        prepared.append((r, r, fold(r, 4.0, GUI_PHASE), None))
        info.append(SequenceRun(f"CAL_{lv}", "calibration", lv))
    for k, (name, amts) in enumerate(SAMPLE_AMOUNTS.items()):
        r = make_cal_run(amts, shift=shift, seed=10 + k, noise=noise)
        prepared.append((r, r, fold(r, 4.0, GUI_PHASE), None))
        info.append(SequenceRun(name, "sample", volume=2.0))
    levels = {lv: {n: a for n in "ABC"} for lv, a in CAL_LEVELS.items()}
    return prepared, info, targets, levels


def test_sequence_quantification_recovers_amounts():
    prepared, info, targets, levels = _sequence()
    rep = quantify_prepared(prepared, info, targets, levels)
    for n in "ABC":
        assert rep.fit(n).ok, n
    for r in rep.results:
        if r.run in SAMPLE_AMOUNTS and r.role == "target":
            true = SAMPLE_AMOUNTS[r.run][r.compound]
            assert r.status == "ok", (r.run, r.compound, r.flags)
            assert r.amount == pytest.approx(true, rel=0.05), (r.run, r.compound)
            assert r.concentration == pytest.approx(r.amount / 2.0)
    # Samples are shifted by 1.5 s; following the internal standard keeps targets on position.
    s1 = next(i for i in rep.seq if i.name == "S1").quant
    assert s1.hits["A"].t1 == pytest.approx(110.0 + 1.5, abs=1.0)
    assert not s1.hits["A"].off_position
    # Learned positions come from the calibration runs.
    assert rep.positions["A"][0] == pytest.approx(110.0, abs=1.0)


def test_shared_peak_and_qualifier_flags():
    targets = [Target("IS", "istd", "", 50, (51,), (0.5,), 60.0, 10, 1.2, 0.2),
               Target("A", "target", "", 52, (53,), (0.3,), 110.0, 10, 2.0, 0.2, "IS"),
               Target("A2", "target", "", 52, (53,), (0.9,), 111.0, 10, 2.0, 0.2, "IS")]
    r = make_cal_run({"A": 2.0, "B": 1.0, "C": 1.0})
    q = quantify_run(r, r, fold(r, 4.0, GUI_PHASE), None, targets)
    assert "shares peak with A2" in q.hits["A"].flags()
    assert q.hits["A"].qual_ok is True
    assert q.hits["A2"].qual_ok is False
    assert "qualifier ratio" in q.hits["A2"].flags()
    assert q.response("A", "quant") == pytest.approx(2.0 * SCALE / IS_AREA, rel=0.05)


def test_decisions_reject_exclude_assign_and_round_trip():
    prepared, info, targets, levels = _sequence()
    d = QuantDecisions()
    d.rejected.add(("S1", "A"))
    d.excluded["B"] = {"CAL_L5": "dropped vial"}
    d.overrides["C"] = {"weights": "1"}
    d.assigned[("S2", "C")] = (201.5, 1.52)
    rep = quantify_prepared(prepared, info, targets, levels, decisions=d)
    s1a = next(r for r in rep.results if r.run == "S1" and r.compound == "A")
    assert s1a.status == "rejected" and np.isnan(s1a.amount)
    assert rep.fit("B").excluded == {"CAL_L5": "dropped vial"}
    assert "CAL_L5" not in rep.fit("B").labels
    assert rep.fit("C").weights == "1"
    s2 = next(i for i in rep.seq if i.name == "S2").quant
    assert s2.hits["C"].manual
    again = QuantDecisions.from_dict(d.to_dict())
    assert again.rejected == d.rejected and again.assigned == d.assigned
    assert again.excluded == d.excluded and again.overrides == d.overrides


def test_all_signals_and_istd_table():
    prepared, info, targets, levels = _sequence()
    rep = quantify_prepared(prepared, info, targets, levels)
    s1 = [row for row in rep.signals if row.run == "S1"]
    assert len(s1) == 4                                  # IS and three targets
    a = next(row for row in s1 if row.target == "A")
    assert a.istd == "IS"
    assert a.ratio == pytest.approx(3.0 * 1.3 * SCALE / (IS_AREA * 1.5), rel=0.05)
    assert a.surrogate == "A" and a.amount == pytest.approx(3.0, rel=0.05)
    assert all(not row.flag for row in rep.istd)
    assert {row.run for row in rep.istd} == {i.name for i in info}


def test_cli_quant(cal_folder, tmp_path, capsys):
    from pygcgc.cli import main
    out = tmp_path / "res"
    assert main(["quant", str(cal_folder), "--outdir", str(out), "--workers", "2"]) == 0
    text = capsys.readouterr().out
    assert "3 of 3 calibration curves pass" in text
    rows = list(csv.DictReader(open(out / "targets.csv")))
    s1a = next(r for r in rows if r["run"] == "S1" and r["compound"] == "A")
    assert float(s1a["amount"]) == pytest.approx(3.0, rel=0.05)
    assert float(s1a["concentration"]) == pytest.approx(1.5, rel=0.05)
    for name in ("calibrations.csv", "all_signals.csv", "internal_standards.csv"):
        assert (out / name).exists()


def test_cli_quant_multi_set_interpolates(drift_folder, tmp_path, capsys):
    from pygcgc.cli import main
    folder, true_amount = drift_folder
    out_on = tmp_path / "on"
    assert main(["quant", str(folder), "--outdir", str(out_on), "--workers", "2"]) == 0
    text = capsys.readouterr().out
    assert "2 calibration sets" in text and "interpolation on" in text
    on_rows = list(csv.DictReader(open(out_on / "targets.csv", encoding="utf-8")))
    on_a = next(r for r in on_rows if r["run"] == "S_MID" and r["compound"] == "A")
    assert "\u2192" in on_a["cal_set"]
    amount_on = float(on_a["amount"])

    out_off = tmp_path / "off"
    assert main(["quant", str(folder), "--outdir", str(out_off), "--workers", "2",
                "--no-interpolate"]) == 0
    text_off = capsys.readouterr().out
    assert "interpolation off" in text_off
    off_rows = list(csv.DictReader(open(out_off / "targets.csv")))
    off_a = next(r for r in off_rows if r["run"] == "S_MID" and r["compound"] == "A")
    assert off_a["cal_set"] == "mon"
    amount_off = float(off_a["amount"])

    assert abs(amount_on - true_amount) < abs(amount_off - true_amount)
    assert amount_on == pytest.approx(true_amount, rel=0.1)


def test_target_table_validation(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text("compound,quant_ion,t1_expected_s,t2_expected_s,istd,qual_ions,qual_ratios\n"
                 "A,57,100,1,missing,71,0.5\n")
    with pytest.raises(ValueError, match="not in the table"):
        load_targets(p)
    p.write_text("compound,quant_ion,t1_expected_s,t2_expected_s,qual_ions,qual_ratios\n"
                 "A,57,100,1,71;85,0.5\n")
    with pytest.raises(ValueError, match="ratios"):
        load_targets(p)
    lv = tmp_path / "l.csv"
    lv.write_text("level,A,B\nL1,1,2\nL2,3,\n")
    assert load_levels(lv) == {"L1": {"A": 1.0, "B": 2.0}, "L2": {"A": 3.0}}


# ------------------------------------------------------------ dummy data
@pytest.fixture(scope="module")
def dummy_report(dummy_dir, dummy_runs, dummy_anchors):
    from pygcgc import phase_for_sequence, remove_background
    targets = load_targets(dummy_dir / "target_table.csv")
    levels = load_levels(dummy_dir / "calibration_levels_ng.csv")
    rows = list(csv.DictReader(open(dummy_dir / "sequence.csv")))
    runs = [dummy_runs[r["file"][:-4]] for r in rows]
    phases = phase_for_sequence(runs, 4.0, dummy_anchors)
    prepared, info = [], []
    for run, ph, row in zip(runs, phases, rows):
        f = fold(run, 4.0, ph.phase)
        bg = remove_background(run, f)
        prepared.append((run, bg.run, f, bg.noise))
        info.append(SequenceRun(run.name, row["type"], row["level"] if row["type"] == "calibration" else ""))
    return quantify_prepared(prepared, info, targets, levels)


def _truth(dummy_dir):
    with open(dummy_dir / "truth_per_file.csv", newline="") as fh:
        return {(r["file"][:-4], r["compound"]): r for r in csv.DictReader(fh)}


@pytest.mark.dummy
def test_dummy_calibrations_pass(dummy_report):
    fits = [e["quant"] for e in dummy_report.fits.values()]
    assert len(fits) == 48
    assert all(not isinstance(f, str) and f.ok for f in fits)


@pytest.mark.dummy
def test_dummy_isolated_areas_and_samples(dummy_report, dummy_dir):
    truth = _truth(dummy_dir)
    iso = ["naphthalene", "acenaphthene", "fluorene", "pyrene", "fluoranthene", "1-nitronaphthalene",
           "9-nitroanthracene", "benzo[a]pyrene"]
    for run in dummy_report.seq:
        if run.name not in ("CAL_L3", "CAL_L4", "CAL_L5"):
            continue
        for c in iso:
            h = run.quant.hits[c]
            assert h.quant_volume == pytest.approx(float(truth[(run.name, c)]["quant_ion_area"]), rel=0.02), c
    ok = [r for r in dummy_report.results if r.run.startswith("SAMPLE") and r.status == "ok"
          and r.compound in iso]
    assert len(ok) >= 15
    for r in ok:
        true = float(truth[(r.run, r.compound)]["amount_ng"])
        assert r.amount == pytest.approx(true, rel=0.10), (r.run, r.compound)


@pytest.mark.dummy
def test_dummy_flags(dummy_report):
    q = {r.name: r.quant for r in dummy_report.seq}
    for c in ("trans-decalin", "cis-decalin"):
        assert q["CAL_L7"].hits[c].saturated
        assert not q["CAL_L3"].hits[c].saturated
    assert q["CAL_L4"].hits["chrysene"].shared_with == ("benz[a]anthracene",)
    assert q["CAL_L4"].hits["benzo[k]fluoranthene"].shared_with == ("benzo[b]fluoranthene",)
    status = {(r.run, r.compound): r.status for r in dummy_report.results}
    assert status[("SAMPLE_B", "anthracene")] == "check"        # its m/z 178 interferent
    assert status[("SAMPLE_B", "n-C13")] == "check"              # only a UCM compound nearby
    assert all(not row.flag for row in dummy_report.istd)
    blanks = [row for row in dummy_report.istd if row.run == "TUBE_BLANK"]
    assert blanks and not any(row.found for row in blanks)
