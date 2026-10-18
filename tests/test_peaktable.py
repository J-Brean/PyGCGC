"""Fitted components as a downstream peak list, and the signal accounting."""

from dataclasses import replace

import numpy as np
import pytest
from scipy import sparse

from pygcgc import ComponentDetectionSettings, Run, detect_components, fold
from pygcgc.core.peaktable import (ComponentPeaks, format_accounting, respectrum_components,
                                   signal_accounting)
from pygcgc.core.synthetic import TruePeak, cosine, inject_peaks

MZ = np.arange(40, 100)                      # 60 ions: more than a window's 30, so truncation is real


def _spectrum(seed, n_ions=60):
    rng = np.random.default_rng(seed)
    s = np.zeros(n_ions)
    s[rng.choice(n_ions, 14, replace=False)] = rng.uniform(300, 3000, 14)
    return s


@pytest.fixture(scope="module")
def scene():
    period, pitch, n_scans = 10.0, 0.02, 12000                     # 20 modulations of 500 rows
    times = np.arange(n_scans) * pitch
    empty = Run(times, MZ, sparse.csr_matrix((n_scans, MZ.size)), "synthetic")
    folded = fold(empty, period, 0.0)
    spectra = [_spectrum(1), _spectrum(2)]
    truth = [TruePeak(200.0, 6.0, 4e6, spectra[0], 3.0, 0.8),
             TruePeak(340.0, 13.0, 2.5e6, spectra[1], 3.0, 0.8)]
    run = inject_peaks(empty, folded, truth, seed=3)
    settings = replace(ComponentDetectionSettings(), workers=1)
    _cands, _regions, found = detect_components(run, folded, "S", settings, floor=5.0)
    return run, folded, found, truth, spectra


def test_components_are_found_and_respectrummed_on_every_ion(scene):
    run, folded, found, truth, spectra = scene
    peaks, model, footprint = respectrum_components(run, folded, found)
    assert isinstance(peaks, ComponentPeaks) and len(peaks) == len(found.components) >= 2
    assert peaks.spectra.shape == (len(peaks), MZ.size)
    for planted, spectrum in zip(truth, spectra):
        i = int(np.argmin(np.hypot((peaks.t2 - planted.row * 0.02) / 0.1,
                                   (peaks.t1 - (planted.col * 10.0 + 5.0)) / 10.0)))
        assert cosine(peaks.spectra[i], spectrum) > 0.97          # all 60 ions, not just the window's 30
        assert peaks.volume[i] == pytest.approx(planted.area, rel=0.2)


def test_adapter_offers_what_classify_and_cluster_read(scene):
    run, folded, found, *_ = scene
    peaks, *_ = respectrum_components(run, folded, found)
    n = len(peaks)
    assert peaks.mz.size == peaks.spectra.shape[1]
    assert peaks.co_elution.shape == (n,) and peaks.t1.shape == (n,) and peaks.volume.shape == (n,)
    assert peaks.apex_spectra.shape[0] != n                       # so they fall back to `spectra`
    x, y = peaks.outline_times(0)
    assert len(x) == len(y) and len(x) % 3 == 0 and np.isnan(x[2])   # (start, end, NaN) triples
    from pygcgc.core.chemclass import classify_blobs
    from pygcgc.core.cluster import PeakTable
    assert len(classify_blobs(peaks).calls) == n
    assert len(PeakTable.from_blobs(peaks, "S")) == n


def test_signal_accounting_explains_the_planted_signal_and_reports_the_rest(scene):
    run, folded, found, *_ = scene
    _peaks, model, footprint = respectrum_components(run, folded, found)
    report = signal_accounting(folded.tic_image(run), model, footprint, folded.row_pitch,
                               folded.spec.period, float(folded.t1[0]))
    assert report["percent_explained"] > 85.0
    assert report["percent_unexplained_outside_all_peaks"] < 10.0
    assert "explain" in format_accounting(report)


def test_signal_accounting_flags_a_peak_the_model_does_not_contain(scene):
    run, folded, found, *_ = scene
    _peaks, model, footprint = respectrum_components(run, folded, found)
    tic = folded.tic_image(run).copy()
    tic[100:130, 15:17] += 4e5                                    # a peak nothing fitted
    report = signal_accounting(tic, model, footprint, folded.row_pitch, folded.spec.period,
                               float(folded.t1[0]))
    top = report["largest_unexplained_patches"][0]
    assert not top["inside_a_fitted_footprint"]
    assert top["t2_s"] == pytest.approx(115 * 0.02, abs=0.4)
    assert report["percent_unexplained_outside_all_peaks"] > 5.0


def test_downstream_spectra_use_the_same_constrained_solve_as_the_fitter(scene):
    from pygcgc.core.candidates import CandidateRegion, _fold_region
    from pygcgc.core.components import component_profiles_from_durable
    from pygcgc.core.overlap import BlobCube, initialise_overlap, refit_spectra
    from pygcgc.core.peaks import PeakSettings, make_frame

    run, folded, found, *_ = scene
    peaks, *_ = respectrum_components(run, folded, found)
    offset = 0
    frame = make_frame(folded, PeakSettings())
    for group in found.groups:
        region = CandidateRegion((), group.row0, group.row0 + group.n_rows,
                                  group.col0, group.col0 + group.n_cols)
        data = _fold_region(run, folded, np.arange(run.n_mz), region, frame.wrap_rows)
        cube = BlobCube(data, np.ones(data.shape[:2], bool), run.mz,
                         group.row0, group.col0, group.pitch, group.period, group.t1_start)
        # Reconstruct shapes through the durable representation, independently
        # of the downstream adapter's coefficients.
        shaped = initialise_overlap(cube, len(group.components))
        shaped.model = 'template'
        shaped.profiles = np.stack([component_profiles_from_durable(c, cube.mask.shape)
                                    for c in group.components], axis=2)
        fit = refit_spectra(shaped, cube)
        for component in fit.components:
            np.testing.assert_allclose(peaks.spectra[offset], component.spectrum, rtol=2e-6, atol=1e-5)
            assert peaks.volume[offset] == pytest.approx(component.volume, rel=2e-6)
            offset += 1
