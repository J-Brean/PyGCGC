"""🧹 background removal, despiking, blank subtraction and the noise map."""

import csv

import numpy as np
import pytest
import scipy.sparse as sp

from pygcgc import (BackgroundSettings, Run, arpls, despike, estimate_noise, fold,
                    remove_background)
from tests.conftest import PERIOD, TRUE_PHASE, make_run

PEAKS = [(62.0, 1.20, 50, 2e5), (122.0, 2.60, 51, 3e5), (182.0, 2.0, 53, 2e5)]


def test_blank_subtraction_matches_mass_not_column_position():
    from pygcgc.core.io.andi import crop_run
    sample = make_run([(22.,1.,51,2e5),(22.,1.,53,2e5)],duration=48,mz=(50,55))
    blank = make_run([(22.,1.,51,1e5)],duration=48,mz=(50,55))
    narrow = crop_run(blank,0,48,mz_min=51,mz_max=51)
    settings = BackgroundSettings(despike=False,remove_pattern=False,subtract_blank=True,scale_blank=False)
    sf,bf,nf = fold(sample,4),fold(blank,4),fold(narrow,4)
    full = remove_background(blank,bf)
    small = remove_background(narrow,nf)
    a = remove_background(sample,sf,settings,blank=full,blank_fold=bf)
    b = remove_background(sample,sf,settings,blank=small,blank_fold=nf)
    np.testing.assert_allclose(a.run.data.toarray(),b.run.data.toarray())
    assert b.run.xic(51).sum() < sample.xic(51).sum()


def with_background(run, ion=53, level=(40.0, 120.0), noise=4.0, stripe=0.0,
                    stripe_t2=3.0, spikes=0, seed=0, phase=0.0, period=4.0):
    """Add a rising baseline (plus optional stripe and spikes) to one ion."""
    rng = np.random.default_rng(seed)
    dense = run.data.toarray().astype(np.float64)
    j = ion - run.mz[0]
    frac = (run.times - run.times[0]) / (run.times[-1] - run.times[0])
    base = level[0] + (level[1] - level[0]) * frac ** 2
    t2 = (run.times - run.times[0] - phase) % period
    base = base + stripe * np.exp(-0.5 * ((t2 - stripe_t2) / 0.05) ** 2)
    dense[:, j] += base + rng.normal(0, noise, run.n_scans)
    spike_at = []
    if spikes:
        spike_at = rng.choice(np.arange(10, run.n_scans - 10), spikes, replace=False)
        dense[spike_at, 50 - run.mz[0]] += rng.uniform(500, 2000, spikes)
    out = Run(times=run.times, mz=run.mz, data=sp.csr_matrix(dense), name="bg")
    return out, base, np.sort(spike_at)


def box_sum(run, f, t1, t2, ion, half_rows=0.4, cols=4):
    img = f.ion_image(run, ion)
    c = int((t1 - f.spec.phase) // f.spec.period)
    r = int(round(t2 / f.row_pitch))
    h = int(half_rows / f.row_pitch)
    return img[r - h:r + h + 1, c - cols:c + cols + 1].sum()


# ------------------------------------------------------------------ arPLS
def test_arpls_follows_baseline_under_peaks():
    rng = np.random.default_rng(0)
    x = np.linspace(0, 1, 465)
    base = 50 + 100 * x ** 2
    y = base + rng.normal(0, 2, x.size)
    for c in (80, 200, 330):
        y[c - 3:c + 4] += 600 * np.exp(-0.5 * ((np.arange(7) - 3) / 1.5) ** 2)
    z = arpls(y, 1e5)
    assert np.max(np.abs(z - base)) < 3.0


def test_arpls_short_inputs():
    np.testing.assert_allclose(arpls(np.array([1.0, 2.0])), [1.0, 2.0])
    assert arpls(np.array([1.0, 2.0, 3.0])).shape == (3,)


# ------------------------------------------------------------ synthetic runs
def test_background_centred_and_peak_kept():
    run, base, _ = with_background(make_run(PEAKS, rate=50.0))
    f = fold(run, 4.0)
    res = remove_background(run, f, BackgroundSettings(despike=False))
    assert 53 in res.ions_processed
    img = f.ion_image(res.run, 53)[:, f.complete]
    sig = res.noise.at(53, f.t1[f.complete])
    quiet = np.abs(img) < 3 * sig[None, :]
    assert abs(img[quiet].mean()) < 0.5
    assert res.noise.at(53, 120.0) == pytest.approx(4.0, rel=0.2)
    clean = make_run(PEAKS, rate=50.0)
    expected = box_sum(clean, f, 182.0, 2.0, 53)
    assert box_sum(res.run, f, 182.0, 2.0, 53) == pytest.approx(expected, rel=0.02)
    # Ions without background are untouched.
    np.testing.assert_allclose(res.run.xic(51), run.xic(51))


def test_repeating_pattern_removed_only_when_asked():
    run, _, _ = with_background(make_run(PEAKS, rate=50.0), stripe=60.0, stripe_t2=3.0)
    f = fold(run, 4.0)
    row = int(round(3.0 / f.row_pitch))
    on = remove_background(run, f, BackgroundSettings(despike=False))
    off = remove_background(run, f, BackgroundSettings(despike=False, remove_pattern=False))
    stripe_on = f.ion_image(on.run, 53)[row, f.complete].mean()
    stripe_off = f.ion_image(off.run, 53)[row, f.complete].mean()
    assert stripe_off > 20
    assert abs(stripe_on) < 3
    # The real peak on this ion survives pattern removal.
    clean = make_run(PEAKS, rate=50.0)
    assert box_sum(on.run, f, 182.0, 2.0, 53) == pytest.approx(box_sum(clean, f, 182.0, 2.0, 53), rel=0.03)


def test_despike_removes_spikes_not_peaks():
    run, _, spikes = with_background(make_run(PEAKS, rate=50.0), spikes=25, seed=3)
    cleaned, n = despike(run, 8.0, floor=3.0)
    assert n == 25
    np.testing.assert_allclose(cleaned.xic(50)[spikes], make_run(PEAKS, rate=50.0).xic(50)[spikes],
                               atol=1.0)
    f = fold(run, 4.0)
    assert box_sum(cleaned, f, 62.0, 1.2, 50) == pytest.approx(
        box_sum(make_run(PEAKS, rate=50.0), f, 62.0, 1.2, 50), rel=1e-4)


def test_blank_subtraction_removes_artefact_and_keeps_analyte():
    artefact = (100.0, 1.5, 52, 1e5)
    sample = make_run(PEAKS + [artefact], rate=50.0)
    blank = make_run([(100.0, 1.5, 52, 5e4)], rate=50.0)       # half as much in the blank
    f = fold(sample, 4.0)
    bres = remove_background(blank, f)
    s = BackgroundSettings(subtract_blank=True)
    res = remove_background(sample, f, s, blank=bres, blank_fold=f, blank_name="blank")
    assert res.blank_scale[52] == pytest.approx(2.0, rel=0.05)
    assert box_sum(res.run, f, 100.0, 1.5, 52) < 0.05 * box_sum(sample, f, 100.0, 1.5, 52)
    assert box_sum(res.run, f, 62.0, 1.2, 50) == pytest.approx(box_sum(sample, f, 62.0, 1.2, 50))
    assert res.run.data.data.min() >= 0          # never below zero from the blank step
    unscaled = remove_background(sample, f, BackgroundSettings(subtract_blank=True, scale_blank=False),
                                 blank=bres, blank_fold=f, blank_name="blank")
    assert unscaled.blank_scale[52] == 1.0


def test_blank_subtraction_after_alignment_follows_retention_drift():
    """The blank's artefact sits 12 s later in t1 and 0.2 s higher in t2 than
    the sample's. Read pixel for pixel it is missed; with the mapping built
    from the two runs' alignment to the reference it is removed."""
    from pygcgc.core.align import AnchorPoint, blank_mapping, fit_alignment
    sample = make_run(PEAKS + [(100.0, 1.5, 52, 1e5)], rate=50.0, duration=260.0)
    blank = make_run([(112.0, 1.7, 52, 5e4)], rate=50.0, duration=260.0)
    f = fold(sample, 4.0)
    bres = remove_background(blank, f)
    s = BackgroundSettings(subtract_blank=True)
    naive = remove_background(sample, f, s, blank=bres, blank_fold=f, blank_name="blank")
    assert box_sum(naive.run, f, 100.0, 1.5, 52) > 0.8 * box_sum(sample, f, 100.0, 1.5, 52)
    # The sample is the reference; the blank's one anchor says (112, 1.7) there is (100, 1.5) here.
    t_blank = fit_alignment([AnchorPoint("art", 100.0, 1.5, 112.0, 1.7, 52)], "auto")
    aligned = remove_background(sample, f, s, blank=bres, blank_fold=f, blank_name="blank",
                                to_blank=blank_mapping(None, t_blank))
    assert box_sum(aligned.run, f, 100.0, 1.5, 52) < 0.1 * box_sum(sample, f, 100.0, 1.5, 52)
    assert box_sum(aligned.run, f, 62.0, 1.2, 50) == pytest.approx(box_sum(sample, f, 62.0, 1.2, 50))


def test_blank_subtraction_scales_two_artefacts_on_one_ion_independently():
    # Two artefacts on the same ion, well apart in t1, that vary
    # independently between tubes: the sample has twice the blank's amount
    # of the first, and half the blank's amount of the second. A single,
    # shared scale factor (the old behaviour) would only ever get the
    # average right, leaving one of them under- or over-corrected.
    sample = make_run(PEAKS + [(100.0, 1.5, 52, 1e5), (200.0, 1.5, 52, 3e4)], rate=50.0, duration=260.0)
    blank = make_run([(100.0, 1.5, 52, 5e4), (200.0, 1.5, 52, 6e4)], rate=50.0, duration=260.0)
    f = fold(sample, 4.0)
    bres = remove_background(blank, f)
    res = remove_background(sample, f, BackgroundSettings(subtract_blank=True), blank=bres, blank_fold=f,
                            blank_name="blank")
    assert box_sum(res.run, f, 100.0, 1.5, 52) < 0.05 * box_sum(sample, f, 100.0, 1.5, 52)
    assert box_sum(res.run, f, 200.0, 1.5, 52) < 0.05 * box_sum(sample, f, 200.0, 1.5, 52)
    assert box_sum(res.run, f, 62.0, 1.2, 50) == pytest.approx(box_sum(sample, f, 62.0, 1.2, 50))


def test_settings_round_trip():
    s = BackgroundSettings(percentile=15, stiffness=1e6, remove_pattern=False)
    assert BackgroundSettings.from_dict(s.to_dict()) == s
    assert BackgroundSettings.from_dict({"percentile": 30, "unknown": 1}).percentile == 30


def test_noise_floor_for_sparse_ions():
    run = make_run(PEAKS, rate=50.0)
    nm = estimate_noise(run, fold(run, 4.0), floor=3.0)
    assert np.all(nm.sigma == 3.0)
    assert nm.for_ions([50, 51], np.array([10.0]))[0] == pytest.approx(3.0 * np.sqrt(2))


# ------------------------------------------------------------------ dummy data
def _truth(dummy_dir):
    with open(dummy_dir / "truth_per_file.csv", newline="") as fh:
        return {(r["file"].replace(".cdf", ""), r["compound"]): r for r in csv.DictReader(fh)}


def _box(run, f, t, mz, before=0.4, after=0.6, cb=5, ca=6):
    col, t2 = f.to_fold(float(t["apex_scan_time_s"]))
    img, _ = f.extend(f.ion_image(run, mz), 1.99)
    r0 = int(round(t2 / f.row_pitch))
    return img[max(r0 - int(before / f.row_pitch), 0):r0 + int(after / f.row_pitch),
               col - cb:col + ca].sum()


@pytest.fixture(scope="module")
def cleaned_dummy(dummy_runs):
    out = {}
    for name in ("CAL_L1", "CAL_L4", "SAMPLE_A", "TUBE_BLANK"):
        run = dummy_runs[name]
        f = fold(run, PERIOD, TRUE_PHASE)
        out[name] = (run, f, remove_background(run, f))
    return out


@pytest.mark.dummy
def test_dummy_bleed_ions_centred(cleaned_dummy):
    run, f, res = cleaned_dummy["CAL_L4"]
    for mz in (40, 44, 73, 207, 281):
        raw = f.ion_image(run, mz)[:, f.complete]
        img = f.ion_image(res.run, mz)[:, f.complete]
        sig = res.noise.at(mz, f.t1[f.complete])
        quiet = np.abs(img) < 3 * sig[None, :]
        assert np.median(raw) > 15, mz
        assert abs(img[quiet].mean()) < 0.25 * np.median(sig), mz
        assert img[quiet].std() == pytest.approx(np.median(sig), rel=0.3), mz


@pytest.mark.dummy
def test_dummy_areas(cleaned_dummy, dummy_dir):
    truth = _truth(dummy_dir)
    for name in ("CAL_L4", "SAMPLE_A"):
        run, f, res = cleaned_dummy[name]
        for compound in ("naphthalene", "fluorene", "pyrene", "1-nitropyrene"):
            t = truth[(name, compound)]
            mz = int(t["quant_ion"])
            assert _box(res.run, f, t, mz) == pytest.approx(_box(run, f, t, mz), rel=0.005), compound
        d4 = truth[(name, "siloxane D4")]
        assert _box(res.run, f, d4, 281) == pytest.approx(float(d4["quant_ion_area"]), rel=0.01)
    # Behenic acid is quantified on m/z 73, a bleed ion; its tail needs a wide box.
    for name in ("CAL_L1", "CAL_L4", "SAMPLE_A"):
        run, f, res = cleaned_dummy[name]
        t = truth[(name, "behenic acid")]
        area = _box(res.run, f, t, 73, before=0.4, after=2.2, cb=4, ca=12)
        raw = _box(run, f, t, 73, before=0.4, after=2.2, cb=4, ca=12)
        assert area == pytest.approx(float(t["quant_ion_area"]), rel=0.05), name
        assert raw > 1.05 * float(t["quant_ion_area"]), name


@pytest.mark.dummy
def test_dummy_spikes(cleaned_dummy):
    for name in ("CAL_L4", "SAMPLE_A"):
        _, _, res = cleaned_dummy[name]
        # The generator adds about 86 spikes; a few land on peaks and are left.
        assert 75 <= res.spikes_removed <= 90, name


@pytest.mark.dummy
def test_dummy_blank_subtraction(cleaned_dummy, dummy_dir, dummy_runs):
    truth = _truth(dummy_dir)
    brun, bf, bres = cleaned_dummy["TUBE_BLANK"]
    run, f, plain = cleaned_dummy["SAMPLE_A"]
    res = remove_background(run, f, BackgroundSettings(subtract_blank=True), blank=bres,
                            blank_fold=bf, blank_name="TUBE_BLANK")
    arts = ("benzaldehyde", "phenol", "acetophenone", "benzoic acid", "siloxane D5", "siloxane D6")
    before = sum(_box(plain.run, f, truth[("SAMPLE_A", a)], int(truth[("SAMPLE_A", a)]["quant_ion"]))
                 for a in arts)
    after = sum(_box(res.run, f, truth[("SAMPLE_A", a)], int(truth[("SAMPLE_A", a)]["quant_ion"]))
                for a in arts)
    assert after < 0.35 * before
    for compound in ("naphthalene", "fluorene", "pyrene"):
        t = truth[("SAMPLE_A", compound)]
        mz = int(t["quant_ion"])
        assert _box(res.run, f, t, mz) == pytest.approx(_box(plain.run, f, t, mz), rel=1e-6)


@pytest.mark.dummy
def test_dummy_speed(dummy_runs):
    import time
    run = dummy_runs["SAMPLE_B"]
    f = fold(run, PERIOD, TRUE_PHASE)
    t = time.perf_counter()
    remove_background(run, f)
    assert time.perf_counter() - t < 5.0
