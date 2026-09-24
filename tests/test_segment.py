"""Pixel-level regions: spectra of the scans themselves, no peak detection."""

import numpy as np
import pytest

pytest.importorskip("sklearn")

from pygcgc import fold                                                       # noqa: E402
from pygcgc.core.segment import (SegmentSettings, fit_regions, region_shares,  # noqa: E402
                                 segment_pixels)
from tests.conftest import make_run                                            # noqa: E402

ALKANE_IONS = (57, 71, 85, 99, 113)
AROMATIC_IONS = (65, 77, 91, 105, 119)
WEIGHTS = (1.0, 0.7, 0.4, 0.2, 0.1)
PERIOD, PHASE = 4.0, 0.7


def _run(area_b=2e5, seed=1, shift=0.0):
    """Six alkane-like and six aromatic-like peaks, well apart in t1."""
    peaks, where = [], {"A": [], "B": []}
    for i in range(6):
        for fam, ions in (("A", ALKANE_IONS), ("B", AROMATIC_IONS)):
            t1 = 61.3 + 130.0 * i + (55.0 if fam == "B" else 0.0) + shift
            t2 = 1.207 + 0.25 * i + (0.4 if fam == "B" else 0.0)
            area = area_b if fam == "B" else 2e5
            peaks += [(t1, t2, ion, area * w) for ion, w in zip(ions, WEIGHTS)]
            where[fam].append((t1, t2))
    run = make_run(peaks, phase=PHASE, rate=50.0, mz=(50, 130), duration=900.0, seed=seed)
    return run, fold(run, PERIOD, PHASE), where


def _region_at(labels, f, t1, t2):
    """The commonest region in a 3 x 3 patch around a peak's apex pixel."""
    start = PHASE + np.floor((t1 - PHASE) / PERIOD) * PERIOD
    col = int(np.floor((start + t2 - PHASE) / PERIOD))
    row = int(round(t2 / f.row_pitch))
    patch = labels[row - 1:row + 2, col - 1:col + 2].ravel()
    patch = patch[patch >= 0]
    return int(np.bincount(patch).argmax()) if patch.size else -1


def test_two_spectral_families_fall_into_two_regions():
    run, f, where = _run()
    result, labels = segment_pixels(run, f, SegmentSettings(n_regions=2, keep_percentile=50))
    assert labels.shape == f.shape and len(result) == 2
    a = {_region_at(labels, f, *p) for p in where["A"]}
    b = {_region_at(labels, f, *p) for p in where["B"]}
    assert len(a) == 1 and len(b) == 1 and a != b and -1 not in a | b
    assert result.top_ions(next(iter(a)), 1) == [57] and result.top_ions(next(iter(b)), 1) == [65]


def test_regions_are_numbered_by_signal_and_baseline_stays_unlabelled():
    run, f, where = _run(area_b=1e6)                                    # the aromatic family holds most of the signal
    result, labels = segment_pixels(run, f, SegmentSettings(n_regions=2, keep_percentile=50))
    assert result.top_ions(0, 1) == [65]
    assert (labels == -1).mean() > 0.5                                  # most of the plane is quiet
    shares = region_shares(labels, run, f, 2)
    assert shares.sum() == pytest.approx(100.0) and shares[0] > shares[1]


def test_a_fitted_model_labels_another_run_the_same_way():
    run1, f1, where = _run(seed=1)
    run2, f2, where2 = _run(seed=2, shift=6.0)                           # drifted in t1, as a real second run would be
    result = fit_regions([(run1, f1)], SegmentSettings(n_regions=2, keep_percentile=50))
    l1, l2 = result.label_image(run1, f1), result.label_image(run2, f2)
    for fam in ("A", "B"):
        assert ({_region_at(l1, f1, *p) for p in where[fam]} == {_region_at(l2, f2, *p) for p in where2[fam]}
                and len({_region_at(l2, f2, *p) for p in where2[fam]}) == 1)


def test_fitting_on_several_runs_at_once_and_the_summary():
    run1, f1, where = _run(seed=1)
    run2, f2, _ = _run(seed=2, shift=6.0)
    result = fit_regions([(run1, f1), (run2, f2)], SegmentSettings(n_regions=2, keep_percentile=50))
    labels = result.label_image(run1, f1)
    rows = result.summary(labels, run1, f1)
    assert [r.label for r in rows] == [0, 1] and sum(r.percent for r in rows) == pytest.approx(100.0)
    assert all(r.n_pixels > 0 and r.t1_median > 0 for r in rows)


def test_smoothing_removes_speckle_without_growing_into_baseline():
    run, f, _ = _run()
    raw = segment_pixels(run, f, SegmentSettings(n_regions=2, keep_percentile=50, smooth=0))[1]
    smooth = segment_pixels(run, f, SegmentSettings(n_regions=2, keep_percentile=50, smooth=2))[1]
    assert ((smooth >= 0) == (raw >= 0)).all()                          # the labelled area is unchanged


def test_settings_round_trip_and_bad_input():
    s = SegmentSettings(n_regions=7, exclude_ions=(73, 207))
    assert SegmentSettings.from_dict(__import__("json").loads(__import__("json").dumps(s.to_dict()))) == s
    with pytest.raises(ValueError):
        fit_regions([])
    run, f, _ = _run()
    with pytest.raises(ValueError, match="too few"):
        fit_regions([(run, f)], SegmentSettings(n_regions=5000, keep_percentile=99.99))
