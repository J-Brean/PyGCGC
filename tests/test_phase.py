import numpy as np
import pytest

from pygcgc import (Anchor, PhaseEstimationError, auto_phase, fold, locate_anchor,
                    phase_for_sequence, phase_from_anchors, phase_from_void)
from pygcgc.core import circstats as cs
from tests.conftest import PERIOD, TRUE_PHASE, make_run

PEAKS = [(40.0, 1.10, 50, 1e5), (90.0, 1.90, 51, 2e5), (150.0, 2.70, 52, 1.5e5),
         (200.0, 1.40, 53, 1e5), (70.0, 2.20, 54, 1e5)]


def _anchors():
    return [Anchor(f"c{i}", ion, t1, t2) for i, (t1, t2, ion, _) in enumerate(PEAKS)]


@pytest.mark.parametrize("phase", [0.0, 0.73, 3.95])
def test_anchor_phase_synthetic(phase):
    run = make_run(PEAKS, phase=phase, rate=50.0)
    est = phase_from_anchors(run, _anchors(), 4.0)
    assert abs(cs.circ_diff(est.phase, phase, 4.0)) < 0.005
    assert all(h.used for h in est.hits)


def test_anchor_outlier_rejected():
    run = make_run(PEAKS, phase=0.73, rate=50.0)
    anchors = _anchors()
    anchors[0].t2 += 0.8                  # a wrong expected t2
    est = phase_from_anchors(run, anchors, 4.0)
    assert not est.hits[0].used
    assert abs(cs.circ_diff(est.phase, 0.73, 4.0)) < 0.005


def test_missing_anchor_reported():
    run = make_run(PEAKS, phase=0.73)
    assert locate_anchor(run, Anchor("absent", 54, 400.0, 1.0))[0] is None
    with pytest.raises(PhaseEstimationError):
        phase_from_anchors(run, [Anchor("absent", 54, 400.0, 1.0)], 4.0)


def _boundary_fraction(run, phase, half_width=0.1):
    f = fold(run, 4.0, phase)
    prof = f.tic_image(run)
    prof = np.clip(prof - np.percentile(prof, 20, axis=0, keepdims=True), 0, None).sum(axis=1)
    n = int(round(half_width / f.row_pitch))
    near = np.r_[prof[:n], prof[-n:]].sum()
    return near / prof.sum()


def test_void_cut_avoids_peaks():
    run = make_run(PEAKS, phase=0.73, rate=50.0)
    est = phase_from_void(run, 4.0)
    assert est.method == "void"
    assert _boundary_fraction(run, est.phase) < 1e-3
    # A deliberately bad cut straight through a peak row is much worse.
    assert _boundary_fraction(run, (0.73 + 1.9) % 4.0) > 0.05


def _first_signal_row(run, phase, period=4.0):
    f = fold(run, period, phase)
    prof = f.tic_image(run)
    prof = np.clip(prof - np.percentile(prof, 20, axis=0, keepdims=True), 0, None).sum(axis=1)
    return f.t2[np.flatnonzero(prof > 0.01 * prof.max())[0]]


def test_void_cut_sits_just_before_the_signal():
    run = make_run(PEAKS, phase=0.73, rate=50.0)
    est = phase_from_void(run, 4.0)
    assert est.confidence == "high" and est.contrast > 0.8
    # Signal starts soon after the cut, not half a period away.
    assert 0.05 < _first_signal_row(run, est.phase) < 0.9


def test_void_confidence_is_low_when_wraparound_leaves_no_gap():
    # Compounds spread over the whole period, so there is nowhere quiet to cut.
    t2s = np.linspace(0.1, 3.9, 14)
    peaks = [(40.0 + 12 * i, t2, 50 + i % 5, 2e5) for i, t2 in enumerate(t2s)]
    est = phase_from_void(make_run(peaks, phase=0.73, rate=50.0, duration=240.0, sigma1=30.0), 4.0)
    assert est.confidence != "high"


def test_void_ignores_a_solvent_front():
    from tests.test_period import _add_solvent_front
    run = _add_solvent_front(make_run(PEAKS, phase=0.73, rate=50.0), centre=20.0, width=8.0)
    est = phase_from_void(run, 4.0)
    assert est.solvent is not None and est.solvent[0] < 20.0 < est.solvent[1]
    assert _boundary_fraction(make_run(PEAKS, phase=0.73, rate=50.0), est.phase) < 1e-3


def test_void_refuses_a_run_with_no_signal():
    import scipy.sparse as sp
    from pygcgc import Run
    empty = Run(times=np.arange(2000) / 50.0, mz=np.arange(50, 52), data=sp.csr_matrix((2000, 2)))
    with pytest.raises(PhaseEstimationError):
        phase_from_void(empty, 4.0)


def test_sequence_without_anchors_lends_a_confident_phase_to_a_weak_run():
    good = make_run(PEAKS, phase=0.73, rate=50.0)
    blank = make_run([(120.0, 0.9, 50, 3e3)], phase=0.73, rate=50.0, noise=30.0, seed=4)
    ests = phase_for_sequence([good, blank], 4.0, None)
    assert ests[0].method == "void"
    assert ests[1].method == "fallback"
    assert ests[1].phase == pytest.approx(ests[0].phase)


def test_auto_prefers_anchors_then_fallback_then_void():
    run = make_run(PEAKS, phase=0.73, rate=50.0)
    assert auto_phase(run, 4.0, _anchors()).method == "anchors"
    absent = [Anchor("x", 54, 400.0, 1.0)]
    assert auto_phase(run, 4.0, absent, fallback=0.5).method == "fallback"
    assert auto_phase(run, 4.0, absent).method == "void"


# ------------------------------------------------------------------ dummy data
@pytest.mark.dummy
def test_dummy_calibration_runs(dummy_runs, dummy_anchors):
    for name, run in dummy_runs.items():
        if name.startswith("CAL"):
            est = phase_from_anchors(run, dummy_anchors, PERIOD)
            assert abs(cs.circ_diff(est.phase, TRUE_PHASE, PERIOD)) <= 0.02, name
            assert sum(h.used for h in est.hits) == len(dummy_anchors), name


@pytest.mark.dummy
def test_dummy_samples(dummy_runs, dummy_anchors):
    # Samples carry run-wide t2 shifts of up to 0.04 s, which the phase absorbs.
    for name in ("SAMPLE_A", "SAMPLE_B", "SAMPLE_C"):
        est = phase_from_anchors(dummy_runs[name], dummy_anchors, PERIOD)
        assert abs(cs.circ_diff(est.phase, TRUE_PHASE, PERIOD)) <= 0.05, name


@pytest.mark.dummy
def test_dummy_sequence_fallback(dummy_runs, dummy_anchors):
    names = sorted(dummy_runs)
    runs = [dummy_runs[n] for n in names]
    ests = dict(zip(names, phase_for_sequence(runs, PERIOD, dummy_anchors)))
    assert ests["TUBE_BLANK"].method == "fallback"
    assert abs(cs.circ_diff(ests["TUBE_BLANK"].phase, TRUE_PHASE, PERIOD)) <= 0.02
    assert all(e.method == "anchors" for n, e in ests.items() if n != "TUBE_BLANK")


@pytest.mark.dummy
@pytest.mark.parametrize("name", ["CAL_L4", "CAL_L7", "SAMPLE_B"])
def test_dummy_void_cut_is_quiet(dummy_runs, name):
    run = dummy_runs[name]
    est = phase_from_void(run, PERIOD)
    assert _boundary_fraction(run, est.phase, 0.05) < 0.005
