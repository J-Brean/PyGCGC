"""The baseline-subtracted spectrum under a clicked point."""

import numpy as np
import pytest

from pygcgc import fold, ratios_to, spectrum_at
from tests.conftest import make_run

PERIOD = 4.0


def _run_with_bleed_and_a_peak():
    # A peak (ions 100, 150, 200 in 5:2:1) on top of a steady bleed on m/z 207 and 44
    # that fills the whole modulation.
    peaks = [(100.0, 1.5, 100, 5e5), (100.0, 1.5, 150, 2e5), (100.0, 1.5, 200, 1e5)]
    run = make_run(peaks, phase=0.7, rate=50.0, mz=(40, 260), duration=260.0)
    times = run.times
    data = run.data.tolil()
    bleed_cols = [207 - run.mz[0], 44 - run.mz[0]]
    for j in bleed_cols:
        data[:, j] = 300.0
    from pygcgc import Run
    return Run(times, run.mz, data.tocsr(), "bleed")


def test_spectrum_at_a_peak_recovers_the_ion_ratios_and_drops_the_steady_background():
    run = _run_with_bleed_and_a_peak()
    f = fold(run, PERIOD, 0.7)
    apex_time = run.times[int(np.argmax(run.xic(100)))]
    col, apex_t2 = f.to_fold(apex_time)
    t1 = float(f.t1[int(col)] + 0.5 * PERIOD)
    apex_t2 = float(apex_t2)
    with_baseline = spectrum_at(run, f, t1, apex_t2, subtract_baseline=False)
    clean = spectrum_at(run, f, t1, apex_t2)
    assert with_baseline[207 - 40] > 0 and with_baseline[44 - 40] > 0
    assert clean[207 - 40] == pytest.approx(0.0, abs=1.0) and clean[44 - 40] == pytest.approx(0.0, abs=1.0)
    assert clean.min() >= 0
    r = ratios_to(clean, run.mz, 100, (150, 200))
    assert r == pytest.approx((0.4, 0.2), abs=0.03)


def test_spectrum_at_a_point_with_no_scans_is_zero():
    run = _run_with_bleed_and_a_peak()
    f = fold(run, PERIOD, 0.7)
    spec = spectrum_at(run, f, t1=1e6, t2=1.0)
    assert spec.shape == (run.n_mz,) and not spec.any()


def test_ratios_to_a_missing_quantifier_are_zero():
    mz = np.arange(40, 60)
    assert ratios_to(np.zeros(20), mz, 50, (51, 52)) == (0.0, 0.0)
    spec = np.zeros(20)
    spec[10] = 4.0
    spec[11] = 2.0
    assert ratios_to(spec, mz, 50, (51, 99)) == (0.5, 0.0)          # 99 is outside the range
