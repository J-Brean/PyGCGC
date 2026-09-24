"""MS-aware candidate perception before component fitting or integration."""

import numpy as np

from pygcgc import (CandidateSet, CandidateSettings, candidate_cube,
                    candidate_regions, find_candidates, fold)
from pygcgc.core.peaks import Frame
from tests.conftest import make_run


def test_groups_several_ions_supporting_one_peak():
    peaks = [(100.0, 2.0, ion, 2e4) for ion in (60, 61, 62, 63)]
    run = make_run(peaks, phase=.7, rate=50., mz=(58, 66), noise=2., seed=4)
    candidates = find_candidates(run, fold(run, 4., .7), floor=1.)
    near = np.flatnonzero((abs(candidates.t1 - 100.) < 3.) &
                          (abs(candidates.t2 - 2.) < .1))
    assert len(near) == 1
    assert candidates.support_count[near[0]] >= 4
    assert set((60, 61, 62, 63)).issubset(candidates.support_ions[near[0]])


def test_distinct_ion_apexes_make_distinct_candidates_without_hard_borders():
    run = make_run([(100., 1.8, 60, 2e5), (100., 2.15, 61, 2e5)],
                   phase=.7, rate=50., mz=(58, 63), noise=1., seed=2)
    settings = CandidateSettings(single_ion_snr=10.)
    candidates = find_candidates(run, fold(run, 4., .7), settings, floor=1.)
    assert np.any(abs(candidates.t2 - 1.8) < .08)
    assert np.any(abs(candidates.t2 - 2.15) < .08)


def test_unstructured_noise_does_not_accumulate_across_ions():
    run = make_run([], phase=.7, rate=50., mz=(50, 100), noise=2., seed=7)
    candidates = find_candidates(run, fold(run, 4., .7), floor=1.)
    assert len(candidates) <= 2


def test_candidate_regions_join_neighbours_and_keep_distant_peaks_separate():
    candidates = CandidateSet(
        t1=np.array([10., 14., 50.]), t2=np.array([1., 1.1, 2.]),
        row=np.array([10, 12, 20]), col=np.array([2, 3, 10]),
        score=np.ones(3), strongest_snr=np.ones(3), support_count=np.ones(3, int),
        support_ions=[np.array([50]), np.array([51]), np.array([52])],
        support_snr=[np.ones(1)] * 3, pitch=.05, period=4., t1_start=0.)
    regions = candidate_regions(candidates, (100, 20), link_t2_s=.15,
                                link_modulations=1, pad_t2_s=.1)
    assert [r.candidate_indices for r in regions] == [(0, 1), (2,)]
    assert regions[0].row0 == 8 and regions[0].row1 == 15


def test_candidate_cube_uses_local_support_without_a_legacy_footprint():
    run = make_run([(100., 2., 60, 2e5), (100., 2.1, 61, 1e5)],
                   phase=.7, rate=50., mz=(58, 70))
    folded = fold(run, 4., .7)
    candidates = find_candidates(run, folded, CandidateSettings(single_ion_snr=8.), floor=1.)
    regions = candidate_regions(candidates, (folded.n_rows + 20, folded.n_cols))
    region = min(regions, key=lambda r: abs(candidates.t1[r.candidate_indices[0]] - 100.))
    cube = candidate_cube(run, folded, candidates, region)
    assert cube.mask.all()
    assert 60 in cube.ions and 61 in cube.ions
    assert cube.data.shape[:2] == (region.row1-region.row0, region.col1-region.col0)
    expected = Frame(folded, candidates.wrap_rows).extend(
        folded.cube(run, mz=cube.ions, dtype=np.float32))[
            region.row0:region.row1, region.col0:region.col1]
    assert np.allclose(cube.data, expected)
