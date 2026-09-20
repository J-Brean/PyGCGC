import numpy as np
import pytest
import scipy.sparse as sp

from pygcgc import (Anchor, PhaseEstimationError, Run, estimate_period, period_from_anchors,
                    period_from_points)
from pygcgc.core import circstats as cs
from tests.conftest import PERIOD, TRUE_PHASE, make_run

PEAKS = [(30.0, 1.10, 50, 1e5), (90.0, 1.90, 51, 2e5), (170.0, 2.70, 52, 1.5e5),
         (230.0, 1.40, 53, 1e5), (60.0, 2.20, 54, 1e5), (205.0, 3.1, 50, 1e5)]


def _anchors():
    return [Anchor(f"c{i}", ion, t1, t2) for i, (t1, t2, ion, _) in enumerate(PEAKS)]


@pytest.mark.parametrize("true_period", [4.0, 4.004, 3.99])
def test_synthetic_period_and_phase(true_period):
    run = make_run(PEAKS, period=true_period, phase=0.6, rate=50.0, duration=260.0)
    est = period_from_anchors(run, _anchors(), 4.0)
    assert est.period == pytest.approx(true_period, abs=2e-4)
    assert abs(cs.circ_diff(est.phase, 0.6, true_period)) < 0.01
    assert est.significant == (abs(true_period - 4.0) > 1e-3)


def test_points_along_flat_series():
    t0, phase, true_p, c = 0.0, 0.4, 4.003, 1.25
    k = np.array([10, 60, 150, 300, 420])
    clicks = t0 + phase + k * true_p + c
    est = period_from_points(clicks, 4.0, t0)
    assert est.period == pytest.approx(true_p, abs=1e-6)


def test_too_few_anchors():
    run = make_run(PEAKS[:2], rate=50.0)
    with pytest.raises(PhaseEstimationError):
        period_from_anchors(run, _anchors()[:2], 4.0)


# ------------------------------------------------------------ blind estimate
def _many_peaks(duration, n=14, seed=5):
    rng = np.random.default_rng(seed)
    t1 = np.linspace(0.08, 0.92, n) * duration
    return [(float(t), float(rng.uniform(0.6, 3.4)), 50 + i % 5, float(rng.uniform(5e4, 3e5)))
            for i, t in enumerate(t1)]


def _add_solvent_front(run, centre=40.0, width=12.0, height=2e6):
    dense = run.data.toarray()
    dense[:, 0] += height * np.exp(-0.5 * ((run.times - centre) / width) ** 2) * (1 + 0.3 * np.sin(2 * np.pi * run.times / 4.0))
    return Run(times=run.times, mz=run.mz, data=sp.csr_matrix(dense), name=run.name)


def _add_bleed(run, period, phase, t2=1.5, height=3e3, width=0.15):
    """Column bleed: one pulse per modulation, all run long, locked to the modulator."""
    dense = run.data.toarray()
    ph = (run.times - run.times[0] - phase) % period
    dense[:, -1] += height * (1 + run.times / run.times[-1]) * np.exp(-0.5 * ((ph - t2) / width) ** 2)
    return Run(times=run.times, mz=run.mz, data=sp.csr_matrix(dense), name=run.name)


def _run(period, phase=0.6, duration=500.0, bleed=True, sigma1=25.0, noise=20.0, peaks=None):
    run = make_run(_many_peaks(duration) if peaks is None else peaks, period=period, phase=phase,
                   rate=50.0, duration=duration, noise=noise, seed=3, sigma1=sigma1)
    return _add_bleed(run, period, phase) if bleed else run


@pytest.mark.parametrize("true_period,phase,duration", [(4.0, 0.6, 500.0), (4.0037, 2.9, 500.0),
                                                        (10.0, 1.8, 1500.0), (10.0017, 0.4, 1500.0),
                                                        (2.5, 1.1, 300.0)])
def test_blind_period_with_locked_bleed(true_period, phase, duration):
    est = estimate_period(_run(true_period, phase, duration))
    assert est.method == "auto"
    assert est.raw_period == pytest.approx(true_period, abs=1e-3)
    assert est.confidence == "high" and est.precise


@pytest.mark.parametrize("true_period,duration", [(4.0, 500.0), (4.0037, 500.0), (10.0017, 1500.0)])
def test_compounds_alone_give_the_period_to_a_few_ms(true_period, duration):
    # Discrete compounds carry little long-range information, and the result says so.
    est = estimate_period(_run(true_period, duration=duration, bleed=False))
    assert est.confidence in ("high", "medium")
    assert abs(est.raw_period - true_period) < max(4 * est.correction_se, 2e-3)
    assert est.raw_period == pytest.approx(true_period, abs=0.012)


def test_blind_period_snaps_only_when_consistent():
    assert estimate_period(_run(4.0)).period == pytest.approx(4.0, abs=1e-9)
    est = estimate_period(_run(4.0037))
    assert est.period == pytest.approx(4.0037, abs=1e-3) and est.period != 4.0


def test_blind_period_with_solvent_front():
    run = _add_solvent_front(_run(4.0037, duration=600.0))
    est = estimate_period(run)
    assert est.solvent is not None and est.solvent[0] < 40 < est.solvent[1]
    assert est.period == pytest.approx(4.0037, abs=1e-3)
    assert est.confidence == "high"


def test_blind_period_with_leading_empty_stretch():
    run = _run(4.0037, peaks=_many_peaks(500.0)[3:])
    dense = run.data.toarray()
    dense[run.times < 100.0] = 0.0     # the detector was not recording yet
    run = Run(times=run.times, mz=run.mz, data=sp.csr_matrix(dense))
    assert estimate_period(run).period == pytest.approx(4.0037, abs=1e-3)


def test_prior_narrows_the_search_and_flags_disagreement():
    run = _run(4.0037)
    disagree = estimate_period(run, prior=4.0)
    assert disagree.prior_agrees is False and disagree.period == pytest.approx(4.0037, abs=1e-3)
    assert estimate_period(run, prior=4.0037).prior_agrees
    # A wrong prior cannot pull the answer outside its window, and the result is not trusted.
    far = estimate_period(run, prior=6.0, prior_tol=0.03)
    assert 5.8 <= far.period <= 6.2 and far.confidence == "low"


def test_blind_period_noise_only_is_low_confidence():
    rng = np.random.default_rng(0)
    times = np.arange(20000) / 50.0
    run = Run(times=times, mz=np.arange(50, 52), data=sp.csr_matrix(rng.random((20000, 2)) * 100))
    assert estimate_period(run).confidence == "low"


def test_blind_period_refuses_bad_input():
    with pytest.raises(PhaseEstimationError):
        estimate_period(make_run([(20.0, 1.0, 50, 1e5)], rate=50.0, duration=20.0))   # too short
    empty = Run(times=np.arange(5000) / 50.0, mz=np.arange(50, 52), data=sp.csr_matrix((5000, 2)))
    with pytest.raises(PhaseEstimationError):
        estimate_period(empty)


def test_estimate_period_for_runs_picks_the_most_precise_usable_estimate():
    from pygcgc import estimate_period_for_runs
    good = _run(4.0037)
    empty = Run(times=good.times, mz=good.mz, data=sp.csr_matrix(good.data.shape))
    best, found = estimate_period_for_runs([empty, good, empty])
    assert best is not None and best.period == pytest.approx(4.0037, abs=1e-3)
    assert len(found) == 1                       # the empty runs could not be estimated
    assert estimate_period_for_runs([empty])[0] is None


def test_cli_estimates_the_period_when_none_is_given(tmp_path, capsys):
    from pygcgc.cli import main
    from tests.conftest import write_run_andi
    write_run_andi(_run(6.0, phase=1.7, duration=600.0), tmp_path / "RUN.cdf")
    npz = tmp_path / "out.npz"
    assert main(["fold", str(tmp_path / "RUN.cdf"), "--npz", str(npz)]) == 0
    assert "estimated from the run" in capsys.readouterr().out
    assert float(np.load(npz)["period"]) == pytest.approx(6.0, abs=0.01)
    assert main(["phase", str(tmp_path / "RUN.cdf")]) == 0


def test_cli_stops_when_no_period_can_be_found(tmp_path, capsys):
    from pygcgc.cli import main
    from tests.conftest import write_run_andi
    rng = np.random.default_rng(0)
    noise = Run(times=np.arange(20000) / 50.0, mz=np.arange(50, 52),
                data=sp.csr_matrix(rng.random((20000, 2)) * 100))
    write_run_andi(noise, tmp_path / "NOISE.cdf")
    assert main(["fold", str(tmp_path / "NOISE.cdf")]) == 2
    assert "give --period" in capsys.readouterr().out


@pytest.mark.dummy
def test_dummy_nominal_period_confirmed(dummy_runs, dummy_anchors):
    for name, run in dummy_runs.items():
        if name == "TUBE_BLANK":
            continue
        est = period_from_anchors(run, dummy_anchors, PERIOD)
        assert est.period == pytest.approx(PERIOD, abs=1e-4), name
        assert abs(cs.circ_diff(est.phase, TRUE_PHASE, PERIOD)) <= 0.05, name


@pytest.mark.dummy
@pytest.mark.parametrize("eps", [0.001, -0.0005])
def test_dummy_detects_clock_mismatch(dummy_runs, dummy_anchors, eps):
    # Pretend the scan clock runs slightly fast or slow: the apparent period changes.
    run = dummy_runs["CAL_L4"]
    scaled = run.with_times(run.times[0] + (run.times - run.times[0]) * (1 + eps))
    est = period_from_anchors(scaled, dummy_anchors, PERIOD)
    assert est.period == pytest.approx(PERIOD * (1 + eps), abs=1e-4)
    assert est.significant
