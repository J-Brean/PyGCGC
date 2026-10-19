"""📍 Alignment: transform fitting, cross-validation, propagation, resampling."""

import numpy as np
import pytest

from pygcgc import fold
from pygcgc.core.align import (AlignmentError, AnchorPoint, cross_validate, fit_alignment,
                               propagate_anchors, refine_anchor_position, resample_to_reference)
from pygcgc.core.phase import Anchor
from tests.conftest import make_run


def test_refinement_finds_only_the_bounded_local_maximum():
    f = fold(make_run([], rate=100), 4.)
    image = np.zeros(f.shape)
    image[151, 24] = 40
    image[180, 25] = 1000  # brighter, outside the vertical limit
    image[150, 28] = 2000  # outside the modulation limit
    assert refine_anchor_position(image, f, 102., 1.45) == pytest.approx((98., 1.51))
    assert refine_anchor_position(image, f, 102., 1.45, columns=0) is None


@pytest.mark.parametrize("height", [0., -3., np.nan, np.inf, 2.])
def test_refinement_leaves_unsupported_neighbourhoods_unchanged(height):
    f = fold(make_run([]), 4.)
    assert refine_anchor_position(np.full(f.shape, height), f, 102., 1.5) is None


def test_refinement_ties_favour_nearest_pixel_and_do_not_drift_on_repeat():
    f = fold(make_run([], rate=100), 4.)
    image = np.zeros(f.shape)
    image[150, 24:27] = 40
    point = refine_anchor_position(image, f, 102., 1.49)
    assert point == pytest.approx((102., 1.50))
    assert refine_anchor_position(image, f, *point) == point


@pytest.mark.parametrize("clicked_t2,row,column,expected_t2", [
    (3.98, 2, 26, 4.02), (.02, 398, 24, -.02), (4.05, 8, 26, 4.08)])
def test_refinement_tracks_wrapped_pulses_without_a_large_coordinate_jump(clicked_t2, row, column, expected_t2):
    f = fold(make_run([], rate=100), 4.)
    image = np.zeros(f.shape)
    image[row, column] = 10
    assert refine_anchor_position(image, f, 102., clicked_t2, columns=0) == pytest.approx((102., expected_t2))


def test_refinement_at_acquisition_edges_does_not_wrap_columns():
    f = fold(make_run([], rate=100), 4.)
    image = np.zeros(f.shape)
    image[-2, -1] = 1000
    image[3, 0] = 10
    assert refine_anchor_position(image, f, 2., .01, columns=0) == pytest.approx((2., .03))


@pytest.mark.parametrize("columns,radius", [(-1, .1), (1.5, .1), (1, 0.), (1, np.inf)])
def test_refinement_rejects_invalid_search_limits(columns, radius):
    f = fold(make_run([]), 4.)
    with pytest.raises(ValueError):
        refine_anchor_position(np.zeros(f.shape), f, 102., 1.5, columns, radius)


def _affine_points(n, a1=5.0, b1=1.02, a2=0.01, b2=0.0002, seed=0, t2_noise=0.0):
    rng = np.random.default_rng(seed)
    ref_t1 = np.linspace(300, 2000, n)
    ref_t2 = rng.uniform(0.5, 3.5, n)
    run_t1 = (ref_t1 - a1) / b1
    run_t2 = ref_t2 - (a2 + b2 * run_t1) + rng.normal(0, t2_noise, n)
    return [AnchorPoint(f"a{i}", ref_t1[i], ref_t2[i], run_t1[i], run_t2[i]) for i in range(n)]


# ------------------------------------------------------------------- fitting
def test_shift_recovers_a_pure_shift():
    pts = [AnchorPoint("p", 110.0, 1.55, 100.0, 1.50)]
    tf = fit_alignment(pts, kind="shift")
    assert tf.dt1 == pytest.approx(10.0)
    assert tf.dt2 == pytest.approx(0.05)
    t1, t2 = tf.forward(100.0, 1.50)
    assert (t1, t2) == pytest.approx((110.0, 1.55))
    t1i, t2i = tf.inverse(110.0, 1.55)
    assert (t1i, t2i) == pytest.approx((100.0, 1.50))


def test_shift_uses_median_against_one_outlier():
    pts = [AnchorPoint(f"p{i}", 100.0 + i, 1.0, 90.0 + i, 1.0) for i in range(5)]
    pts.append(AnchorPoint("bad", 500.0, 1.0, 90.0, 1.0))   # wildly wrong dt1
    tf = fit_alignment(pts, kind="shift")
    assert tf.dt1 == pytest.approx(10.0, abs=0.6)


def test_affine_recovers_stretch_and_shift():
    pts = _affine_points(5)
    tf = fit_alignment(pts, kind="affine")
    assert tf.a1 == pytest.approx(5.0, abs=1e-6)
    assert tf.b1 == pytest.approx(1.02, abs=1e-6)
    assert tf.a2 == pytest.approx(0.01, abs=1e-6)
    assert tf.b2 == pytest.approx(0.0002, abs=1e-8)
    assert tf.monotonic
    for p in pts:
        t1, t2 = tf.forward(p.run_t1, p.run_t2)
        assert (t1, t2) == pytest.approx((p.ref_t1, p.ref_t2), abs=1e-6)


def test_affine_robust_to_one_bad_anchor():
    pts = _affine_points(6)
    bad = AnchorPoint("bad", pts[0].ref_t1, pts[0].ref_t2, pts[0].run_t1 + 200.0, pts[0].run_t2)
    good = pts[1:]
    tf = fit_alignment(good + [bad], kind="affine")
    huber_err = max(abs(tf.forward(p.run_t1, p.run_t2)[0] - p.ref_t1) for p in good)

    # A plain (non-robust) least-squares fit through the same six points,
    # for comparison: one gross outlier among six leverage points still
    # pulls a bounded-influence fit a little, but nothing like this.
    x = np.array([p.run_t1 for p in good] + [bad.run_t1])
    y = np.array([p.ref_t1 for p in good] + [bad.ref_t1])
    A = np.vstack([np.ones_like(x), x]).T
    a0, b0 = np.linalg.lstsq(A, y, rcond=None)[0]
    ols_err = max(abs((a0 + b0 * p.run_t1) - p.ref_t1) for p in good)

    assert huber_err < 0.2 * ols_err


def test_affine_flags_non_monotonic():
    # A run whose t1 order is the reverse of the reference's: no honest
    # straight line through this is increasing.
    pts = [AnchorPoint("a", 300.0, 1.0, 90.0, 1.0), AnchorPoint("b", 200.0, 1.0, 190.0, 1.0),
          AnchorPoint("c", 100.0, 1.0, 290.0, 1.0)]
    tf = fit_alignment(pts, kind="affine")
    assert not tf.monotonic
    assert tf.b1 < 0
    assert "not increasing" in tf.warning


def test_smooth_matches_nonlinear_drift():
    rng = np.random.default_rng(1)
    run_t1 = np.linspace(300, 2000, 8)
    true_t1 = run_t1 + 5.0 + 0.01 * (run_t1 - 1000) + 2e-6 * (run_t1 - 1000) ** 2
    ref_t2_base = rng.uniform(0.5, 3.5, 8)
    true_t2 = ref_t2_base
    run_t2 = true_t2 - (0.02 * np.sin((run_t1 - 300) / 1700 * np.pi))
    pts = [AnchorPoint(f"a{i}", true_t1[i], true_t2[i], run_t1[i], run_t2[i]) for i in range(8)]
    tf = fit_alignment(pts, kind="smooth")
    assert tf.monotonic
    for p in pts:
        t1, t2 = tf.forward(p.run_t1, p.run_t2)
        assert t1 == pytest.approx(p.ref_t1, abs=0.05)
        assert t2 == pytest.approx(p.ref_t2, abs=0.01)
    # Round trip through the inverse should land back close to the start.
    t1i, t2i = tf.inverse(*tf.forward(run_t1[3], run_t2[3]))
    assert (t1i, t2i) == pytest.approx((run_t1[3], run_t2[3]), abs=1e-6)


def test_too_few_anchors_raises():
    with pytest.raises(AlignmentError):
        fit_alignment([], kind="shift")
    with pytest.raises(AlignmentError):
        fit_alignment(_affine_points(2), kind="affine")
    with pytest.raises(AlignmentError):
        fit_alignment(_affine_points(5), kind="smooth")
    with pytest.raises(AlignmentError):
        fit_alignment(_affine_points(3), kind="not-a-kind")


# --------------------------------------------------------------- auto/cross-val
def test_auto_picks_shift_when_that_is_enough():
    pts = [AnchorPoint(f"p{i}", 100.0 + i * 50, 1.0 + 0.01 * i, 90.0 + i * 50, 1.0 + 0.01 * i)
          for i in range(6)]   # a pure shift, but enough points to try affine/smooth too
    tf = fit_alignment(pts, kind="auto")
    assert tf.kind == "shift"


def test_auto_picks_affine_for_affine_data():
    tf = fit_alignment(_affine_points(9), kind="auto")
    assert tf.kind == "affine"


def test_auto_falls_back_to_shift_below_auto_minimum():
    # 3 anchors: affine can be *fit* (FIT_MIN=3) but not cross-validated
    # (AUTO_MIN=4), so auto should not offer it.
    tf = fit_alignment(_affine_points(3), kind="auto")
    assert tf.kind == "shift"


def test_cross_validate_improves_with_more_anchors():
    few = _affine_points(4, t2_noise=0.02, seed=2)
    many = _affine_points(10, t2_noise=0.02, seed=2)
    r1_few, _ = cross_validate(few, "affine")
    r1_many, _ = cross_validate(many, "affine")
    assert np.isfinite(r1_few) and np.isfinite(r1_many)
    assert r1_many <= r1_few * 1.5    # more anchors should not make it clearly worse


def test_cross_validate_nan_when_not_enough_spare_anchors():
    r1, r2 = cross_validate(_affine_points(3), "affine")
    assert np.isnan(r1) and np.isnan(r2)


# ------------------------------------------------------------------ propagate
def test_propagate_anchors_finds_shifted_peaks():
    ref_anchors = [Anchor("IS1", 50, 100.0, 1.2, t1_window=15.0),
                  Anchor("IS2", 51, 300.0, 2.4, t1_window=15.0)]
    shift = 6.0
    run = make_run([(100.0 + shift, 1.2, 50, 2e5), (300.0 + shift, 2.4, 51, 3e5)], rate=50.0, duration=350.0)
    f = fold(run, 4.0, 0.0)
    found = propagate_anchors(ref_anchors, run, f)
    assert set(found) == {"IS1", "IS2"}
    for name, a in zip(("IS1", "IS2"), ref_anchors):
        pt = found[name]
        assert pt is not None
        assert pt.run_t1 == pytest.approx(a.t1 + shift, abs=0.5)
        assert pt.run_t2 == pytest.approx(a.t2, abs=0.03)
        assert pt.ref_t1 == a.t1 and pt.ref_t2 == a.t2


def test_propagate_anchors_flags_missing():
    ref_anchors = [Anchor("Absent", 50, 900.0, 1.0, t1_window=10.0)]
    run = make_run([(100.0, 1.2, 51, 2e5)], rate=50.0)
    f = fold(run, 4.0, 0.0)
    found = propagate_anchors(ref_anchors, run, f)
    assert found["Absent"] is None


# ------------------------------------------------------------------ resample
def test_resample_shift_recovers_peak_and_conserves_volume():
    run = make_run([(100.0, 1.5, 50, 2e5)], phase=0.5, rate=50.0)
    f_run = fold(run, 4.0, 0.5)
    ref = make_run([(105.0, 1.52, 50, 2e5)], phase=0.5, rate=50.0)
    f_ref = fold(ref, 4.0, 0.5)
    tf = fit_alignment([AnchorPoint("p", 105.0, 1.52, 100.0, 1.5)], kind="shift")
    img = f_run.ion_image(run, 50)
    aligned = resample_to_reference(img, f_run, f_ref, tf)
    ref_img = f_ref.ion_image(ref, 50)
    assert aligned.shape == ref_img.shape
    r_a, c_a = np.unravel_index(np.argmax(aligned), aligned.shape)
    r_r, c_r = np.unravel_index(np.argmax(ref_img), ref_img.shape)
    assert abs(r_a - r_r) <= 1 and abs(c_a - c_r) <= 1
    assert aligned.sum() == pytest.approx(ref_img.sum(), rel=0.02)


def test_resample_identity_gives_back_the_same_image():
    run = make_run([(100.0, 1.5, 50, 2e5), (200.0, 2.5, 51, 1e5)], phase=0.3, rate=50.0)
    f = fold(run, 4.0, 0.3)
    tf = fit_alignment([AnchorPoint("p", 0.0, 0.0, 0.0, 0.0)], kind="shift")   # dt1=dt2=0
    img = f.ion_image(run, 50)
    aligned = resample_to_reference(img, f, f, tf)
    np.testing.assert_allclose(aligned, img, atol=img.max() * 1e-6)


def test_resample_extended_is_taller():
    run = make_run([(100.0, 1.5, 50, 2e5)], phase=0.3, rate=50.0)
    f = fold(run, 4.0, 0.3)
    tf = fit_alignment([AnchorPoint("p", 0.0, 0.0, 0.0, 0.0)], kind="shift")
    img = f.ion_image(run, 50)
    plain = resample_to_reference(img, f, f, tf, extended=False)
    extended = resample_to_reference(img, f, f, tf, extended=True)
    assert extended.shape[0] > plain.shape[0]
    assert extended.shape[1] == plain.shape[1]
