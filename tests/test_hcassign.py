"""Hydrocarbon assignment and the shape of the alkane ladder."""

import numpy as np
import pytest

from pygcgc.core.hcassign import assign_hydrocarbons
from pygcgc.core.hydrocarbons import AlkaneHit, AlkaneSeries, HydrocarbonSettings

MZ = np.arange(40, 201)


def _curved_ladder(carbons, shape="piecewise", degree=3):
    # Spacing between alkanes shrinks steadily with carbon number, so t1 is a curve in carbon, not a line.
    hits = {c: AlkaneHit(carbon=c, t1=100.0 + (c - 10) * (80.0 - 4.0 * (c - 10)), t2=1.0 + 0.03 * (c - 10),
                         height=1e4, qual_ratios=(0.66, 0.44)) for c in carbons}
    return AlkaneSeries(hits=hits, shape=shape, degree=degree)


class _Blobs:
    def __init__(self, rows):
        self.mz = MZ
        self.t1 = np.array([r["t1"] for r in rows])
        self.t2 = np.array([r["t2"] for r in rows])
        self.spectra = np.zeros((len(rows), MZ.size))
        for k, r in enumerate(rows):
            for ion, v in r["ions"].items():
                self.spectra[k, ion - MZ[0]] = v

    def __len__(self):
        return len(self.t1)


ALKANE = {57: 100, 71: 66, 85: 44, 43: 60, 41: 20}
BENZENE = {91: 100, 105: 40, 119: 20, 92: 15, 77: 10}
NOISE = {i: 10 for i in range(41, 200, 4)}


def test_a_smooth_ladder_follows_a_curve_that_straight_segments_cut_across():
    every_other = range(10, 21, 2)
    truth = _curved_ladder(range(10, 21))
    for shape in ("smooth", "polynomial"):
        series = _curved_ladder(every_other, shape)
        straight = _curved_ladder(every_other)
        err_curve = abs(series.t1_of(15.0) - truth.t1_of(15.0))
        err_straight = abs(straight.t1_of(15.0) - truth.t1_of(15.0))
        assert err_curve < err_straight, shape
        assert series.carbon_of(series.t1_of(13.5)) == pytest.approx(13.5, abs=0.05)


def test_every_shape_passes_through_the_hits_and_extrapolates_sensibly():
    for shape in ("piecewise", "smooth"):
        series = _curved_ladder(range(10, 20), shape)
        for c in (10, 14, 19):
            assert series.t1_of(c) == pytest.approx(series.hits[c].t1)
        assert series.t1_of(8) < series.t1_of(10) and series.t1_of(22) > series.t1_of(19)
    poly = _curved_ladder(range(10, 20), "polynomial", 3)
    assert poly.t1_of(25) > poly.t1_of(19) and np.isfinite(poly.t2_of(5))


def test_a_polynomial_averages_over_one_wonky_point():
    series = _curved_ladder(range(10, 20), "polynomial", 2)
    clean = series.t1_of(15)
    series.hits[15].t1 += 12.0                     # one alkane found a modulation early or late
    wonky_straight = _curved_ladder(range(10, 20))
    wonky_straight.hits[15].t1 += 12.0
    assert abs(series.t1_of(15) - clean) < 12.0 * 0.5
    assert wonky_straight.t1_of(15) - clean == pytest.approx(12.0)


def test_blobs_are_assigned_by_position_and_spectrum():
    series = _curved_ladder(range(10, 20))
    on = series.t1_of(14.0), series.t2_of(14.0)
    rows = [dict(t1=on[0], t2=on[1], ions=ALKANE),                                           # on the ladder, alkane spectrum
            dict(t1=series.t1_of(14.5), t2=series.t2_of(14.5), ions=ALKANE),                  # between two carbons
            dict(t1=on[0], t2=on[1] + 0.8, ions=ALKANE),                                      # far above the alkane line
            dict(t1=series.t1_of(13.0), t2=series.t2_of(13.0) + 0.5, ions=BENZENE),
            dict(t1=series.t1_of(12.0), t2=series.t2_of(12.0), ions=NOISE),
            dict(t1=on[0], t2=on[1], ions={57: 100, 71: 5, 85: 3, 43: 30})]                  # alkane-like, wrong ratios
    calls = assign_hydrocarbons(_Blobs(rows), series)
    assert [c.label for c in calls] == ["n-C14", "branched alkane C14-eq", "branched alkane C14-eq",
                                        "alkylbenzene C13-eq", "unassigned", "branched alkane C14-eq"]
    assert calls[2].dt2 == pytest.approx(0.8, abs=1e-6)


def test_the_ladder_shape_changes_which_blobs_sit_on_a_carbon():
    every_other = range(10, 21, 2)
    truth = _curved_ladder(range(10, 21))
    rows = [dict(t1=truth.t1_of(15.0), t2=truth.t2_of(15.0), ions=ALKANE)]
    settings = HydrocarbonSettings(assign_carbon_tol=0.05)
    straight = assign_hydrocarbons(_Blobs(rows), _curved_ladder(every_other), settings)[0]
    curved = assign_hydrocarbons(_Blobs(rows), _curved_ladder(every_other, "smooth"), settings)[0]
    assert curved.label == "n-C15"
    assert straight.label != "n-C15"           # straight segments put the same blob 0.1 carbon away
