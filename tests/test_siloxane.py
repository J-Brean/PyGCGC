"""🫧 Siloxane flagging and removal."""

import numpy as np
import pytest

from pygcgc import (SiloxaneSettings, fold, find_blobs, flag_siloxane_blobs, remove_background,
                    remove_siloxane_blobs, remove_siloxanes, siloxane_scores)
from pygcgc.core.siloxane import (BLEED_REFERENCE, CYCLIC_REFERENCE, PRESETS, SILOXANE_IONS,
                                  SILOXANE_REFERENCE, reference_from_spectrum)
from tests.conftest import make_run


def _box(run, f, t1, t2, ion, half_rows=0.4, cols=4):
    img = f.ion_image(run, ion)
    c = int((t1 - f.spec.phase) // f.spec.period)
    r = int(round(t2 / f.row_pitch))
    h = int(half_rows / f.row_pitch)
    return img[r - h:r + h + 1, c - cols:c + cols + 1].sum()


def _siloxane_and_analyte_run():
    """One blob shaped exactly like the reference siloxane pattern (every
    reference ion in the right ratio) and one isolated analyte peak on m/z 73
    alone (behenic acid, in spirit), well separated in t1."""
    scale = 2000.0
    siloxane = [(100.0, 1.5, m, scale * SILOXANE_REFERENCE[m]) for m in SILOXANE_IONS]
    analyte = (200.0, 2.0, 73, 2e5)
    run = make_run(siloxane + [analyte], phase=0.7, rate=50.0, mz=(73, 429), duration=260.0)
    f = fold(run, 4.0, 0.7)
    return run, f


def test_siloxane_peak_scores_high_and_analyte_scores_zero():
    run, f = _siloxane_and_analyte_run()
    clean = remove_background(run, f).run
    blobs = find_blobs(clean, f, floor=1.0)
    assert len(blobs) == 2
    order = np.argsort(blobs.t1)          # siloxane peak (t1=100) first, analyte (t1=200) second
    scores = siloxane_scores(blobs)[order]
    assert scores[0] > 0.9                # matches the reference pattern closely
    assert scores[1] == 0.0               # only one ion present: not scored at all
    flags = flag_siloxane_blobs(blobs)[order]
    assert flags[0] and not flags[1]


def test_remove_siloxanes_strips_the_pattern_but_keeps_the_isolated_analyte():
    run, f = _siloxane_and_analyte_run()
    clean = remove_background(run, f).run
    res = remove_siloxanes(clean, f)
    sil_before = _box(clean, f, 100.0, 1.5, 73)
    sil_after = _box(res.run, f, 100.0, 1.5, 73)
    beh_before = _box(clean, f, 200.0, 2.0, 73)
    beh_after = _box(res.run, f, 200.0, 2.0, 73)
    assert sil_after < 0.1 * sil_before                  # the siloxane peak is mostly gone
    assert beh_after == pytest.approx(beh_before, rel=0.05)   # behenic acid within 5% (the milestone's own bar)
    assert res.run.data.data.min() >= 0                  # never taken below zero


def test_blob_removal_blocks_the_siloxane_blob_and_leaves_the_analyte_alone():
    run, f = _siloxane_and_analyte_run()
    clean = remove_background(run, f).run
    blobs = find_blobs(clean, f, floor=1.0)
    res = remove_siloxane_blobs(clean, f, blobs)
    assert res.n_blobs == 1
    # Blocked out on every ion; the footprint stops at 1% of the apex, so a faint edge tail remains.
    assert _box(res.run, f, 100.0, 1.5, 207) < 0.01 * _box(clean, f, 100.0, 1.5, 207)
    assert _box(res.run, f, 100.0, 1.5, 73) < 0.01 * _box(clean, f, 100.0, 1.5, 73)
    assert _box(res.run, f, 200.0, 2.0, 73) == pytest.approx(_box(clean, f, 200.0, 2.0, 73))   # analyte untouched
    assert res.run.data.nnz < clean.data.nnz
    assert clean.data.sum() > 0 and res.removed_fraction > 0     # the input run itself is unchanged
    assert "1 siloxane blob" in res.summary()


def test_blob_removal_follows_the_threshold_and_can_keep_the_other_ions():
    run, f = _siloxane_and_analyte_run()
    clean = remove_background(run, f).run
    blobs = find_blobs(clean, f, floor=1.0)
    strict = remove_siloxane_blobs(clean, f, blobs, SiloxaneSettings(score_threshold=1.01))
    assert strict.n_blobs == 0 and strict.run is clean
    # Only the siloxane ions: a non-siloxane ion inside the footprint survives.
    run2 = make_run([(100.0, 1.5, m, 2000.0 * SILOXANE_REFERENCE[m]) for m in SILOXANE_IONS]
                    + [(100.0, 1.5, 300, 5e4)], phase=0.7, rate=50.0, mz=(73, 429), duration=260.0)
    f2 = fold(run2, 4.0, 0.7)
    clean2 = remove_background(run2, f2).run
    blobs2 = find_blobs(clean2, f2, floor=1.0)
    res = remove_siloxane_blobs(clean2, f2, blobs2, SiloxaneSettings(block_all_ions=False))
    assert _box(res.run, f2, 100.0, 1.5, 207) < 0.01 * _box(clean2, f2, 100.0, 1.5, 207)
    assert _box(res.run, f2, 100.0, 1.5, 300) == pytest.approx(_box(clean2, f2, 100.0, 1.5, 300))


def test_isolated_single_ion_pixel_is_never_touched():
    # Fewer than min_ions_present siloxane ions present at all: the fitted
    # component is forced to zero regardless of how the one ion looks.
    run = make_run([(100.0, 1.5, 73, 1e5)], phase=0.7, rate=50.0, mz=(73, 300))
    f = fold(run, 4.0, 0.7)
    res = remove_siloxanes(run, f)
    assert sum(res.removed_by_ion.values()) == 0.0


def test_a_real_coincidence_in_the_wrong_ratio_is_not_removed():
    # Three siloxane-listed ions rising together, but nowhere near the
    # reference proportions (equal amounts, not 45:20:25) -- a coincidence,
    # not siloxane. The per-pixel cosine gate should leave it alone.
    run = make_run([(100.0, 1.5, 147, 5e4), (100.0, 1.5, 207, 5e4), (100.0, 1.5, 281, 5e4)],
                   phase=0.7, rate=50.0, mz=(73, 300))
    f = fold(run, 4.0, 0.7)
    clean = remove_background(run, f).run
    res = remove_siloxanes(clean, f)
    for ion in (147, 207, 281):
        before = _box(clean, f, 100.0, 1.5, ion)
        after = _box(res.run, f, 100.0, 1.5, ion)
        assert after == pytest.approx(before, rel=0.02)


def test_no_siloxane_ions_in_range_is_a_no_op():
    run = make_run([(100.0, 1.5, 73, 1e5)], phase=0.7, rate=50.0, mz=(70, 80))
    f = fold(run, 4.0, 0.7)
    res = remove_siloxanes(run, f)
    assert res.ions_used == ()
    assert res.run is run
    assert res.summary()


def test_settings_round_trip():
    s = SiloxaneSettings(score_threshold=0.5, min_ions_present=4)
    assert SiloxaneSettings.from_dict(s.to_dict()) == s
    assert SiloxaneSettings.from_dict({"score_threshold": 0.6, "unknown": 1}).score_threshold == 0.6


def test_default_reference_is_the_measured_column_bleed():
    s = SiloxaneSettings()
    assert s.reference_dict() == BLEED_REFERENCE
    assert max(BLEED_REFERENCE, key=BLEED_REFERENCE.get) == 207      # 207 is the base ion, not 73
    assert BLEED_REFERENCE[73] < 10
    assert max(CYCLIC_REFERENCE, key=CYCLIC_REFERENCE.get) == 73
    assert set(PRESETS) >= {"Column bleed (measured)", "Cyclic siloxane peak (measured)"}


def test_settings_from_reference_and_round_trip():
    s = SiloxaneSettings.from_reference({281: 50.0, 73: 100.0, 147: 20.0}, min_ions_present=2)
    assert s.ions == (73, 147, 281) and s.reference == (100.0, 20.0, 50.0)
    assert s.min_ions_present == 2
    assert SiloxaneSettings.from_dict(s.to_dict()) == s


def test_reference_from_spectrum_keeps_strong_ions_and_drops_low_masses():
    mz = np.arange(38, 301)
    spec = np.zeros(mz.size)
    spec[73 - 38], spec[207 - 38], spec[281 - 38] = 100.0, 40.0, 1.0     # 281 is below min_rel
    spec[44 - 38] = 500.0                                                 # CO2: below min_mz, must be ignored
    ref = reference_from_spectrum(spec, mz, min_rel=0.02)
    assert ref == {73: 100.0, 207: 40.0}
    assert reference_from_spectrum(np.zeros(mz.size), mz) == {}


def test_a_picked_spectrum_removes_a_peak_the_default_would_leave():
    """Column bleed is the default reference; a discrete Dn peak (73 with 147
    and 281) doesn't resemble it, so it is left alone until the reference is
    replaced by that peak's own spectrum."""
    run = make_run([(100.0, 1.5, m, 2000.0 * v) for m, v in CYCLIC_REFERENCE.items()],
                   phase=0.7, rate=50.0, mz=(73, 429), duration=260.0)
    f = fold(run, 4.0, 0.7)
    clean = remove_background(run, f).run
    before = _box(clean, f, 100.0, 1.5, 73)
    kept = _box(remove_siloxanes(clean, f).run, f, 100.0, 1.5, 73)
    gone = _box(remove_siloxanes(clean, f, SiloxaneSettings.from_reference(CYCLIC_REFERENCE)).run, f, 100.0, 1.5, 73)
    assert kept == pytest.approx(before, rel=0.02)
    assert gone < 0.1 * before


# ------------------------------------------------------------------ dummy data
@pytest.mark.dummy
def test_dummy_siloxane_score_separates_real_targets_from_planted_siloxanes(dummy_runs, dummy_truth):
    """The provisional `score_threshold` (see the module docstring) hasn't
    been chosen from a real ROC curve yet -- this is that curve's first
    real check. If it starts failing once the dummy (or real) data are
    available, that's the signal to actually run the ROC analysis the
    milestone asked for and update the default, not to loosen the test."""
    from tests.conftest import PERIOD, TRUE_PHASE
    run = dummy_runs["CAL_L4"]
    f = fold(run, PERIOD, TRUE_PHASE)
    clean = remove_background(run, f).run
    blobs = find_blobs(clean, f)
    scores = siloxane_scores(blobs)

    def score_near(apex_scan_time_s):
        col, t2 = f.to_fold(float(apex_scan_time_s))
        t1 = float(f.t1[col] + 0.5 * f.spec.period)
        i = blobs.nearest(t1, t2)
        return float(scores[i]) if i is not None else 0.0

    siloxanes = [c for c in ("siloxane D4", "siloxane D5", "siloxane D6")
                if ("CAL_L4", c) in dummy_truth]
    targets = [c for c in ("naphthalene", "fluorene", "pyrene", "1-nitropyrene", "behenic acid")
              if ("CAL_L4", c) in dummy_truth]
    assert siloxanes and targets, "dummy truth layout changed; update the compound names above"
    sil_scores = [score_near(dummy_truth[("CAL_L4", c)]["apex_scan_time_s"]) for c in siloxanes]
    target_scores = [score_near(dummy_truth[("CAL_L4", c)]["apex_scan_time_s"]) for c in targets]
    assert min(sil_scores) > max(target_scores)
