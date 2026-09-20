import numpy as np
import pytest

from pygcgc import Fold, FoldSpec, fold
from tests.conftest import PERIOD, TRUE_PHASE, make_run

PEAKS = [(62.0, 1.20, 50, 1e5), (122.0, 2.60, 51, 2e5), (182.0, 4.30, 52, 1.5e5)]


def _valid_total(run, f):
    rel = run.times - f.spec.t0 - f.spec.phase
    return run.tic[rel >= 0].sum()


@pytest.mark.parametrize("phase", [0.0, 0.7, 2.35])
def test_intensity_conserved(phase):
    run = make_run(PEAKS, phase=phase)
    f = fold(run, 4.0, phase)
    assert f.tic_image(run).sum() == pytest.approx(_valid_total(run, f), rel=1e-9)


def test_intensity_conserved_with_fractional_period():
    run = make_run(PEAKS, rate=20.37)
    f = fold(run, 4.0, 0.0)
    assert (4.0 / f.row_pitch) == f.n_rows
    assert f.tic_image(run).sum() == pytest.approx(_valid_total(run, f), rel=1e-9)


def test_peak_lands_where_expected():
    phase = 0.7
    run = make_run(PEAKS, phase=phase)
    f = fold(run, 4.0, phase)
    img = f.ion_image(run, 50)
    r, c = np.unravel_index(np.argmax(img), img.shape)
    assert r == round(1.20 / f.row_pitch)
    assert c == int((62.0 - phase) // 4.0)


def test_wrapped_peak_and_extension():
    phase = 0.7
    run = make_run(PEAKS, phase=phase)
    f = fold(run, 4.0, phase)
    img = f.ion_image(run, 52)
    k = int((182.0 - phase) // 4.0)
    r, c = np.unravel_index(np.argmax(img), img.shape)
    # A 4.3 s compound released in modulation k arrives 0.3 s into modulation k+1.
    assert (r, c) == (round(0.30 / f.row_pitch), k + 1)
    ext, t2 = f.extend(img, 1.3)
    # In the extended image the same peak also appears whole on top of column k.
    r_ext = round(4.30 / f.row_pitch)
    assert int(np.argmax(ext[:, k])) == r_ext
    assert ext[r_ext, k] == img[round(0.30 / f.row_pitch), k + 1]
    assert t2[r_ext] == pytest.approx(4.30)


def test_no_diagonal_creep_with_fractional_period():
    # 4 s at 20.37 Hz is not a whole number of scans; peaks must stay level.
    run = make_run([(120.0, 2.0, 50, 1e6)], rate=20.37, sigma1=10.0, duration=240.0)
    f = fold(run, 4.0, 0.0)
    img = f.ion_image(run, 50)
    cols = np.flatnonzero(img.sum(axis=0) > 0.05 * img.sum(axis=0).max())
    centres = [(img[:, c] * f.t2).sum() / img[:, c].sum() for c in cols]
    assert len(cols) >= 5
    assert np.ptp(centres) < 0.3 * f.row_pitch
    assert np.mean(centres) == pytest.approx(2.0, abs=0.5 * f.row_pitch)


def test_coordinate_round_trip():
    run = make_run(PEAKS)
    f = Fold(run.times, FoldSpec(period=4.0, phase=0.55, t0=300.0 - 300.0))
    t = np.array([10.0, 57.3, 199.99])
    col, t2 = f.to_fold(t)
    np.testing.assert_allclose(f.to_time(col, t2), t)
    assert np.all((t2 >= 0) & (t2 < 4.0))


def test_cube_matches_ion_images():
    run = make_run(PEAKS, phase=0.3)
    f = fold(run, 4.0, 0.3)
    cube = f.cube(run, [50, 52])
    np.testing.assert_allclose(cube[:, :, 0], f.ion_image(run, 50), rtol=1e-5, atol=1e-3)
    np.testing.assert_allclose(cube[:, :, 1], f.ion_image(run, 52), rtol=1e-5, atol=1e-3)
    full = f.cube(run)
    np.testing.assert_allclose(full.sum(axis=2), f.tic_image(run), rtol=1e-5, atol=1e-3)


def test_complete_flags_partial_columns():
    run = make_run(PEAKS, duration=41.0)
    f = fold(run, 4.0, 0.5)
    assert f.complete[0] and not f.complete[-1]


# ------------------------------------------------------------------ dummy data
ISOLATED = ["naphthalene", "acenaphthene-d10", "fluorene", "pyrene", "perylene-d12",
            "1-nitronaphthalene", "naphthalene-d8"]


@pytest.mark.dummy
def test_dummy_apex_positions(dummy_runs, dummy_truth):
    run = dummy_runs["CAL_L4"]
    f = fold(run, PERIOD, TRUE_PHASE)
    for name in ISOLATED:
        t = dummy_truth[("CAL_L4", name)]
        col, t2 = f.to_fold(float(t["apex_scan_time_s"]))
        img = f.ion_image(run, int(t["quant_ion"]))
        r0 = int(round(t2 / f.row_pitch))
        win = img[max(r0 - 10, 0):r0 + 11, col - 2:col + 3]
        dr, dc = np.unravel_index(np.argmax(win), win.shape)
        obs_r, obs_c = max(r0 - 10, 0) + dr, col - 2 + dc
        assert abs(obs_c - col) <= 1, name
        assert abs(obs_r - r0) <= 1, name
        # With the true phase, t2 of the apex equals the compound's true t2.
        assert t2 == pytest.approx(float(t["t2_true_s"]) % PERIOD, abs=0.021), name


@pytest.mark.dummy
def test_dummy_naive_fold_matches_truth(dummy_runs, dummy_truth):
    run = dummy_runs["SAMPLE_B"]
    f = fold(run, PERIOD, 0.0)
    for (file, name), t in dummy_truth.items():
        if file != "SAMPLE_B":
            continue
        _, t2 = f.to_fold(float(t["apex_scan_time_s"]))
        assert t2 == pytest.approx(float(t["folded_t2_s"]), abs=1e-3), name


@pytest.mark.dummy
def test_dummy_wraparound_nitropyrene(dummy_runs, dummy_truth):
    run = dummy_runs["CAL_L4"]
    t = dummy_truth[("CAL_L4", "1-nitropyrene")]
    assert t["wraps"] == "True"
    f = fold(run, PERIOD, TRUE_PHASE)
    col, t2 = f.to_fold(float(t["apex_scan_time_s"]))
    assert t2 < 1.0                       # appears near the bottom of the next modulation
    ext, _ = f.extend(f.ion_image(run, int(t["quant_ion"])), 1.3)
    # Extended, the whole peak sits above the period line, in the modulation that released it.
    r = int(round((t2 + PERIOD) / f.row_pitch))
    assert abs(int(np.argmax(ext[:, col - 1])) - r) <= 1


@pytest.mark.dummy
def test_dummy_conservation(dummy_runs):
    run = dummy_runs["SAMPLE_A"]
    f = fold(run, PERIOD, TRUE_PHASE)
    assert f.tic_image(run).sum() == pytest.approx(_valid_total(run, f), rel=1e-9)
