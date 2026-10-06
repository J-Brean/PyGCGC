"""📈 Calibration builder: matching, IS normalisation, curves, tables and sessions."""

from __future__ import annotations

import numpy as np
import pytest

from pygcgc import fold, find_blobs
from pygcgc.core.calbuilder import (CalBuilderSettings, CalCompound, CalSession, blob_pixel,
                                    curve_rows, load_compounds, load_session, match_blob,
                                    normalised_blob_rows, point_rows, save_session, write_compounds)
from pygcgc.core.calibration import CalibrationFit

from tests.conftest import CAL_LEVELS, GUI_PHASE, IS_AREA, SCALE, make_cal_run

RUNS = [f"CAL_{lv}" for lv in CAL_LEVELS]


@pytest.fixture(scope="module")
def blobs_by_run():
    out = {}
    for k, (lv, amt) in enumerate(CAL_LEVELS.items()):
        run = make_cal_run({n: amt for n in ("A", "B", "C")}, seed=k)
        out[f"CAL_{lv}"] = find_blobs(run, fold(run, 4.0, GUI_PHASE))
    return out


def _session(**kw) -> CalSession:
    amounts = {f"CAL_{lv}": amt for lv, amt in CAL_LEVELS.items()}
    return CalSession(runs=list(RUNS), compounds=[
        CalCompound("IS", 60.0, 1.20, is_istd=True, amounts={r: 1.0 for r in RUNS}),
        CalCompound("A", 110.0, 2.00, amounts=dict(amounts)),
        CalCompound("B", 170.0, 2.60, amounts=dict(amounts)),
        CalCompound("C", 200.0, 1.50, amounts=dict(amounts)),
    ], **kw)


def test_every_compound_is_matched_to_its_own_blob(blobs_by_run):
    s = _session()
    resolved = s.resolve(blobs_by_run)
    for (name, run), i in resolved.items():
        assert i is not None, (name, run)
        c = s.compound(name)
        b = blobs_by_run[run]
        assert abs(b.t1[i] - c.t1) < 4.0 and abs(b.t2[i] - c.t2) < 0.05


def test_a_point_inside_a_footprint_matches_that_blob(blobs_by_run):
    b = blobs_by_run["CAL_L3"]
    i = match_blob(b, 110.0, 2.0)
    row, col = blob_pixel(b, float(b.t1[i]), float(b.t2[i]))
    assert b.at_pixel(row, col) == i
    # Nothing near an empty stretch of the run, even with a generous window.
    assert match_blob(b, 140.0, 3.5, CalBuilderSettings(t1_window_s=5, t2_window_s=0.1)) is None


def test_normalised_curve_recovers_the_known_response(blobs_by_run):
    s = _session()
    points = s.points(blobs_by_run)
    fits = s.fits(points)
    assert set(fits) == {"A", "B", "C"}
    # A: quantifier 52 plus a qualifier at 0.3 on 53, response 1 per unit; IS: 50 plus 0.5 on 51.
    fit = fits["A"]
    assert isinstance(fit, CalibrationFit) and fit.ok
    expected = 1.0 * SCALE * 1.3 / (IS_AREA * 1.5)
    assert fit.slope == pytest.approx(expected, rel=0.03)
    for run in RUNS:
        assert points[("IS", run)].ratio == pytest.approx(1.0)
        assert points[("A", run)].istd == "IS"


def test_quantifier_ion_signal_is_used_when_given(blobs_by_run):
    s = _session()
    s.compound("A").quant_ion = 52
    fit = s.fits(s.points(blobs_by_run))["A"]
    assert fit.slope == pytest.approx(1.0 * SCALE / (IS_AREA * 1.5), rel=0.03)


def test_pins_move_or_remove_a_point_and_unpin_restores_it(blobs_by_run):
    s = _session()
    auto = s.resolve(blobs_by_run)
    s.pin("A", "CAL_L3", 200.0, 1.50)                 # hand-match A onto C's blob in one run
    moved = s.resolve(blobs_by_run)
    assert moved[("A", "CAL_L3")] == auto[("C", "CAL_L3")]
    assert moved[("A", "CAL_L2")] == auto[("A", "CAL_L2")]
    s.pin("A", "CAL_L4", None)                        # no peak in this run
    points = s.points(blobs_by_run)
    assert points[("A", "CAL_L4")].blob is None and np.isnan(points[("A", "CAL_L4")].ratio)
    assert points[("A", "CAL_L3")].pinned
    assert s.fit("A", points).n == 4
    s.unpin("A", "CAL_L3")
    s.unpin("A", "CAL_L4")
    assert s.resolve(blobs_by_run) == auto


def test_pinning_any_blob_finds_that_blob_again(blobs_by_run):
    s = _session()
    b = blobs_by_run["CAL_L2"]
    for i in range(len(b)):
        s.pin_blob("A", "CAL_L2", b, i)
        assert s.resolve(blobs_by_run)[("A", "CAL_L2")] == i


def test_excluding_and_including_a_point(blobs_by_run):
    s = _session()
    points = s.points(blobs_by_run)
    s.exclude("B", "CAL_L5", "")
    fit = s.fit("B", points)
    assert fit.n == 4 and "CAL_L5" in fit.excluded and fit.x_excluded.size == 1
    s.include("B", "CAL_L5")
    assert s.fit("B", points).n == 5 and "B" not in s.excluded
    s.overrides["B"] = {"model": "quadratic", "force_zero": True}
    fit = s.fit("B", points)
    assert fit.model == "quadratic" and fit.zero


def test_too_few_points_gives_a_reason_not_a_curve(blobs_by_run):
    s = _session()
    for run in RUNS[:3]:
        s.pin("C", run, None)
    assert isinstance(s.fit("C", s.points(blobs_by_run)), str)


def test_nearest_or_named_internal_standard():
    s = CalSession(runs=["r"], compounds=[
        CalCompound("IS1", 60.0, 1.0, is_istd=True), CalCompound("IS2", 600.0, 2.0, is_istd=True),
        CalCompound("X", 100.0, 1.1), CalCompound("Y", 550.0, 1.9, istd="IS1")])
    assert s.istd_for("X").name == "IS1"
    assert s.istd_for("Y").name == "IS1"
    s.compound("Y").istd = ""
    assert s.istd_for("Y").name == "IS2"
    assert s.istd_for("IS1") is None


def test_amount_ratio_divides_by_the_internal_standard_amount(blobs_by_run):
    s = _session(settings=CalBuilderSettings(amount_ratio=True))
    for r in RUNS:
        s.compound("IS").amounts[r] = 2.0
    points = s.points(blobs_by_run)
    assert points[("A", "CAL_L5")].x == pytest.approx(CAL_LEVELS["L5"] / 2.0)


def test_rename_and_remove_keep_the_session_consistent():
    s = _session()
    s.pin("A", "CAL_L1", 110.0, 2.0)
    s.exclude("A", "CAL_L2")
    s.compound("B").istd = "IS"
    s.rename("IS", "IS-d10")
    assert s.compound("B").istd == "IS-d10"
    s.rename("A", "anthracene")
    assert ("anthracene", "CAL_L1") in s.pins and "anthracene" in s.excluded
    assert s.rename("B", "C") == "C 2"            # names stay unique
    s.remove_compound("anthracene")
    assert not s.pins and not s.excluded


def test_compound_table_reads_run_columns_level_tokens_and_defaults(tmp_path):
    runs = ["CAL_L1_mon", "CAL_L2_mon", "CAL_L10_mon"]
    p = tmp_path / "compounds.csv"
    p.write_text("compound,t1_min,t2_s,role,quant_ion,istd,amount,L1,L2,CAL_L10_mon\n"
                 "phen-d10,45.1,4.6,istd,188,,2.5,,,\n"
                 "phenanthrene,45.2,4.55,target,178,phen-d10,,0.1,0.2,1.0\n", encoding="utf-8")
    comps = load_compounds(p, runs)
    std, phen = comps
    assert std.is_istd and std.quant_ion == 188 and std.t1 == pytest.approx(2706.0)
    assert std.amounts == {r: 2.5 for r in runs}
    assert phen.istd == "phen-d10" and not phen.is_istd
    assert phen.amounts == {"CAL_L1_mon": 0.1, "CAL_L2_mon": 0.2, "CAL_L10_mon": 1.0}
    out = write_compounds(tmp_path / "out.csv", comps, runs)
    again = load_compounds(out, runs)
    assert [c.to_dict() for c in again] == [c.to_dict() for c in comps]


def test_a_target_table_reads_as_a_compound_table(tmp_path):
    p = tmp_path / "target_table.csv"
    p.write_text("compound,role,quant_ion,qual_ions,t1_expected_s,t2_expected_s,istd\n"
                 "IS,istd,50,51,60,1.2,\nA,target,52,53,110,2.0,IS\n", encoding="utf-8")
    comps = load_compounds(p, ["r1"])
    assert [c.name for c in comps] == ["IS", "A"] and comps[0].is_istd and comps[1].quant_ion == 52
    with pytest.raises(ValueError):
        (tmp_path / "bad.csv").write_text("compound,t2_s\nX,1\n", encoding="utf-8")
        load_compounds(tmp_path / "bad.csv", ["r1"])


def test_amounts_read_from_run_names():
    from pygcgc.core.calbuilder import amount_from_name
    assert amount_from_name("2026-08-25_Cal_L2_5uL_025ng_on_tube_1_70eV") == 0.25
    assert amount_from_name("2026-08-24_L4_1ng_on_tube_1_70eV") == 1.0
    assert amount_from_name("2026-08-24_L8_25ng_on_tube_1_70eV") == 25.0
    assert amount_from_name("CAL_L1_2.5 ng") == 2.5
    assert amount_from_name("SAMPLE_A") is None


def test_suggestions_find_calibrants_and_the_internal_standard(blobs_by_run):
    from pygcgc.core.calbuilder import SuggestSettings, suggest_compounds
    amounts = {f"CAL_{lv}": v for lv, v in CAL_LEVELS.items()}
    found = suggest_compounds(blobs_by_run, amounts, SuggestSettings())
    kinds = {(round(g.t1 / 10) * 10, g.kind) for g in found}
    # IS at 60 s (constant), A 110, B 170, C 200 (one unit per unit amount)
    assert (60, "constant") in kinds
    for t1 in (110, 170, 200):
        assert (t1, "calibrant") in kinds, kinds
    cal = next(g for g in found if g.kind == "calibrant")
    assert cal.slope == pytest.approx(1.0, abs=0.05) and len(cal.blobs) == len(RUNS)
    # Too few levels to say anything.
    assert suggest_compounds(blobs_by_run, dict(list(amounts.items())[:2])) == []


def test_t2_offsets_find_a_run_folded_at_the_wrong_phase():
    from pygcgc.core.calbuilder import t2_offsets
    blobs = {}
    for k, (lv, amt) in enumerate(CAL_LEVELS.items()):
        run = make_cal_run({n: amt for n in ("A", "B", "C")}, seed=k)
        # CAL_L4 folded 0.3 s late: everything in it reads 0.3 s low in t2
        phase = GUI_PHASE + (0.3 if lv == "L4" else 0.0)
        blobs[f"CAL_{lv}"] = find_blobs(run, fold(run, 4.0, phase))
    off = t2_offsets(blobs, "CAL_L1", min_cosine=0.9, min_pairs=3)
    assert off["CAL_L4"][0] == pytest.approx(0.3, abs=0.03)
    for r in ("CAL_L2", "CAL_L3", "CAL_L5"):
        assert abs(off[r][0]) < 0.03
    assert np.isnan(t2_offsets(blobs, "CAL_L1")["CAL_L4"][0])      # too few pairs to say


def test_streaks_and_the_solvent_front_are_not_matched(blobs_by_run):
    from dataclasses import replace
    from pygcgc.core.calbuilder import usable_blobs
    b = blobs_by_run["CAL_L3"]
    assert usable_blobs(b, 4.0).all()
    near_a = np.abs(b.t1 - 110) < 10
    assert not usable_blobs(b, 4.0, solvent=(100.0, 105.0))[near_a].any()
    # three blobs of one base peak stacked up a modulation: a smear
    k = int(np.flatnonzero(near_a)[np.argmax(b.volume[near_a])])
    idx = np.array([k, k, k])
    smear = replace(b, t1=np.r_[b.t1, b.t1[idx]], t2=np.r_[b.t2, [0.6, 2.0, 3.4]],
                    spectra=np.vstack([b.spectra, b.spectra[idx]]), width_t2=np.r_[b.width_t2, [0.1] * 3],
                    volume=np.r_[b.volume, b.volume[idx]])
    keep = usable_blobs(smear, 4.0)
    assert not keep[-3:].any() and keep[np.flatnonzero(~near_a)].all()


def test_shipped_internal_standards_are_recognised_by_spectrum():
    from pygcgc.core.calbuilder import recognise_istd
    mz = np.arange(38, 301)
    s = np.zeros(mz.size)
    s[[188 - 38, 184 - 38, 160 - 38, 80 - 38]] = [100, 19, 13, 8]       # phenanthrene-d10
    assert recognise_istd(s, mz) == "phenanthrene-d10"
    s2 = np.zeros(mz.size)
    s2[[178 - 38, 176 - 38, 152 - 38]] = [100, 20, 8]                    # phenanthrene, undeuterated
    assert recognise_istd(s2, mz) == ""


def test_blob_heights_are_the_largest_raw_point():
    from pygcgc.core.calbuilder import blob_heights
    run = make_cal_run({"A": 5.0, "B": 5.0, "C": 5.0}, seed=3)
    f = fold(run, 4.0, GUI_PHASE)
    b = find_blobs(run, f)
    h = blob_heights(run, f, b)
    assert h.shape == (len(b),) and np.all(h > 0)
    assert h.max() == pytest.approx(run.data.max())


def test_saturated_points_are_flagged(blobs_by_run):
    from dataclasses import replace
    s = _session()
    heights = {r: np.full(len(b), 1e5) for r, b in blobs_by_run.items()}
    i = s.resolve(blobs_by_run)[("A", "CAL_L5")]
    heights["CAL_L5"][i] = 3e5
    s.settings = replace(s.settings, saturation=2.5e5)
    p = s.points(blobs_by_run, heights=heights)
    assert p[("A", "CAL_L5")].saturated and p[("A", "CAL_L5")].height == 3e5
    assert not p[("A", "CAL_L4")].saturated


def test_session_round_trip_and_exports(tmp_path, blobs_by_run):
    s = _session()
    s.pin("A", "CAL_L2", 110.0, 2.0)
    s.pin("B", "CAL_L1", None)
    s.exclude("C", "CAL_L5", "saturated")
    s.overrides["C"] = {"weights": "1/x"}
    path = save_session(tmp_path / "cal.json", s, extra={"period": 4.0})
    back, raw = load_session(path)
    assert back.to_dict() == s.to_dict() and raw["period"] == 4.0
    points = back.points(blobs_by_run)
    fits = back.fits(points)
    rows = point_rows(back, points, fits)
    assert len(rows) == len(RUNS) * 4
    row = next(r for r in rows if r["compound"] == "C" and r["run"] == "CAL_L5")
    assert not row["in_curve"] and row["excluded_because"] == "saturated"
    good = next(r for r in rows if r["compound"] == "A" and r["run"] == "CAL_L3")
    assert good["in_curve"] and abs(good["error_pct"]) < 5
    curves = curve_rows(back, fits)
    assert {r["compound"] for r in curves} == {"A", "B", "C"}
    assert next(r for r in curves if r["compound"] == "C")["weights"] == "1/x"
    blob_rows = normalised_blob_rows(back, blobs_by_run, back.resolve(blobs_by_run))
    is_rows = [r for r in blob_rows if r["compound"] == "IS"]
    assert len(is_rows) == len(RUNS) and all(r["ratio"] == pytest.approx(1.0) for r in is_rows)
