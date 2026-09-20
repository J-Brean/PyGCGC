import numpy as np
import pytest

from pygcgc.core import circstats as cs


def test_mean_across_the_wrap():
    # 3.9 s and 0.1 s are neighbours on a 4 s circle; their mean is 0.0 s, not 2.0 s.
    m = cs.circ_mean([3.9, 0.1], 4.0)
    assert min(m, 4.0 - m) == pytest.approx(0.0, abs=1e-9)


def test_mean_plain_case():
    assert cs.circ_mean([1.0, 1.2, 1.4], 4.0) == pytest.approx(1.2, abs=1e-6)


def test_std_small_for_tight_cluster_and_large_for_spread():
    tight = cs.circ_std([3.98, 3.99, 0.0, 0.01, 0.02], 4.0)
    spread = cs.circ_std(np.linspace(0, 4, 50, endpoint=False), 4.0)
    assert tight < 0.02
    assert spread > 1.0


def test_diff_is_shortest_signed():
    assert cs.circ_diff(0.1, 3.9, 4.0) == pytest.approx(0.2)
    assert cs.circ_diff(3.9, 0.1, 4.0) == pytest.approx(-0.2)


def test_median_ignores_outlier():
    assert 0.50 <= cs.circ_median([0.50, 0.52, 0.51, 2.5], 4.0) <= 0.52
    m = cs.circ_median([3.95, 0.05, 0.0, 0.02, 2.0], 4.0)
    assert abs(cs.circ_diff(m, 0.0, 4.0)) <= 0.05
