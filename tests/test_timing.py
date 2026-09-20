import numpy as np
import pytest

from pygcgc import check_timing, fill_gaps
from tests.conftest import make_run


def test_regular_timing():
    t = 300 + np.arange(5000) / 50.0
    rep = check_timing(t)
    assert rep.rate == pytest.approx(50.0, rel=1e-9)
    assert rep.gaps == []
    assert rep.ok


def test_gaps_found_and_filled(synthetic_peaks):
    run = make_run(synthetic_peaks)
    keep = np.ones(run.n_scans, bool)
    keep[1000:1003] = False          # drop three scans
    keep[3000] = False               # and one more
    gappy = type(run)(times=run.times[keep], mz=run.mz, data=run.data[keep], name="gappy")
    rep = check_timing(gappy.times)
    assert rep.missing_scans == 4
    assert len(rep.gaps) == 2
    assert rep.rate == pytest.approx(20.0, rel=1e-9)   # gaps must not bias the rate
    filled = fill_gaps(gappy, rep)
    assert filled.n_scans == run.n_scans
    np.testing.assert_allclose(filled.times, run.times, atol=1e-9)
    assert filled.data[1000:1003].nnz == 0
    np.testing.assert_allclose(filled.tic[keep], gappy.tic)


@pytest.mark.parametrize("sd, ok", [(0.001, True), (0.0022, False)])
def test_jitter_reported(sd, ok):
    rng = np.random.default_rng(1)
    t = np.arange(10000) / 50.0 + rng.normal(0, sd, 10000)
    rep = check_timing(t)
    assert rep.rms_jitter == pytest.approx(sd, rel=0.1)
    assert rep.gaps == []          # jitter must not be mistaken for dropped scans
    assert rep.ok is ok            # the warning level is 10% of the interval


@pytest.mark.dummy
def test_dummy_timing(dummy_runs):
    for run in dummy_runs.values():
        rep = check_timing(run.times)
        assert rep.ok
        assert rep.rate == pytest.approx(50.0, rel=1e-6)
