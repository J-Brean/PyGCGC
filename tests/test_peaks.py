"""🔍 Blob detection: the sum-of-per-ion-S/N surface, widths, apex spectra,
co-elution clustering, chunked/threaded processing and manual editing."""

import numpy as np
import pytest

from pygcgc import PeakSettings, fold, find_blobs
from pygcgc.core.peaks import (_trim_labels_to_own_height, add_blob_at,
                               delete_blobs, merge_blobs)
from tests.conftest import make_run


def test_widths_and_apex_spectrum():
    sigma1, sigma2 = 6.0, 0.12
    run = make_run([(100.0, 2.0, 60, 3e5)], phase=0.7, rate=50.0, mz=(58, 63),
                   sigma1=sigma1, sigma2=sigma2, duration=200.0)
    f = fold(run, 4.0, 0.7)
    b = find_blobs(run, f, PeakSettings(snr_min=8), floor=1.0)
    assert len(b) == 1
    # FWHM of a Gaussian is 2.3548 sigma; the footprint's own column and row
    # profiles (unsmoothed) should recover the true peak shape closely.
    assert b.width_t1[0] == pytest.approx(2.3548 * sigma1, rel=0.1)
    assert b.width_t2[0] == pytest.approx(2.3548 * sigma2, rel=0.1)
    # The apex spectrum is the scan(s) nearest the single apex pixel, not
    # the whole footprint, so it holds far less than the total volume...
    assert 0 < b.apex_spectra.sum() < 0.1 * b.volume[0]
    # ...but is still on the right ion (m/z 60 is index 2 of 58-63).
    assert np.nonzero(b.apex_spectra[0])[0].tolist() == [60 - run.mz[0]]


def test_isolated_peak_is_not_flagged_as_coelution():
    run = make_run([(100.0, 2.0, 60, 2e5), (100.0, 2.05, 61, 2e5)], phase=0.7, rate=50.0, mz=(58, 63))
    f = fold(run, 4.0, 0.7)
    b = find_blobs(run, f, PeakSettings(snr_min=8, coelution_frac=0.05), floor=1.0)
    assert len(b) == 1              # close enough in t2 to share one footprint
    assert not b.co_elution[0]
    assert b.co_elution_clusters[0] == 1


def test_two_ions_apart_in_one_footprint_are_flagged_as_coelution():
    # Same t1, far enough apart in t2 that the two ions peak in visibly
    # different places, but still close enough to watershed into one blob.
    run = make_run([(100.0, 2.0, 60, 2e5), (100.0, 2.15, 61, 2e5)], phase=0.7, rate=50.0, mz=(58, 63))
    f = fold(run, 4.0, 0.7)
    b = find_blobs(run, f, PeakSettings(snr_min=8, coelution_frac=0.05), floor=1.0)
    assert len(b) == 1
    assert b.co_elution[0]
    assert b.co_elution_clusters[0] == 2


def test_a_faint_minor_ion_does_not_trigger_the_flag():
    # The second ion holds only a trace of the blob's volume, below
    # coelution_frac, so it shouldn't be tested at all.
    run = make_run([(100.0, 2.0, 60, 2e5), (100.0, 2.15, 61, 2e3)], phase=0.7, rate=50.0, mz=(58, 63))
    f = fold(run, 4.0, 0.7)
    b = find_blobs(run, f, PeakSettings(snr_min=8, coelution_frac=0.05), floor=1.0)
    assert len(b) == 1
    assert not b.co_elution[0]


def test_chunked_processing_matches_unchunked_on_a_long_run():
    period = 4.0
    peaks = [(t1, 1.5 + 0.3 * (i % 3), 60 + (i % 3), 2e5)
            for i, t1 in enumerate(np.arange(50.0, 3900.0, 250.0))]
    run = make_run(peaks, phase=0.7, rate=50.0, mz=(58, 63), sigma1=3.0, sigma2=0.08, duration=4000.0)
    f = fold(run, period, 0.7)
    whole = find_blobs(run, f, PeakSettings(snr_min=8, chunk_modulations=100000), floor=1.0)
    chunked = find_blobs(run, f, PeakSettings(snr_min=8, chunk_modulations=40, chunk_overlap=15,
                                              chunk_workers=4), floor=1.0)
    assert len(whole) == len(peaks) == len(chunked)
    ow, oc = np.argsort(whole.t1), np.argsort(chunked.t1)
    np.testing.assert_allclose(whole.t1[ow], chunked.t1[oc], atol=0.5)
    np.testing.assert_allclose(whole.volume[ow], chunked.volume[oc], rtol=0.02)


def test_a_peak_too_narrow_for_one_ion_alone_is_still_found_by_the_ion_sum():
    # Several ions each too faint on their own to clear snr_min against the
    # TIC (which is dominated by an unrelated bright neighbour), but which
    # agree closely enough in position that their summed S/N clears it.
    rng = np.random.default_rng(0)
    bright = (40.0, 1.0, 70, 5e5)          # a big, unrelated peak elsewhere
    faint_ions = [71, 72, 73, 74, 75]
    faint = [(120.0, 2.0, ion, 1.2e4) for ion in faint_ions]
    run = make_run([bright] + faint, phase=0.7, rate=50.0, mz=(60, 80), noise=2.0, seed=3)
    f = fold(run, 4.0, 0.7)
    tic_only = find_blobs(run, f, PeakSettings(snr_min=10, ion_snr=False), floor=1.0)
    with_ion_sum = find_blobs(run, f, PeakSettings(snr_min=10, ion_snr=True), floor=1.0)
    assert len(with_ion_sum) >= len(tic_only)
    assert any(abs(t1 - 120.0) < 2.0 for t1 in with_ion_sum.t1)


def test_final_relative_floor_uses_the_final_owner_and_stays_connected():
    # Pixel (1, 2) was admitted by the tall first marker's threshold, but a
    # subsequent watershed assigned it to the shorter second marker.  It is
    # below that final owner's own 50% floor and must not survive.  Pixel
    # (0, 3) then becomes disconnected from marker 2 and must also be removed.
    labels = np.array([[1, 1, 0, 2],
                       [1, 1, 2, 2],
                       [0, 0, 2, 2]], dtype=np.int32)
    seeds = np.array([[1, 0], [2, 2]])
    intensity = np.array([[10., 8., 0., 8.],
                          [20., 12., 4., 4.],
                          [0., 0., 10., 8.]])
    trimmed = _trim_labels_to_own_height(labels, seeds, intensity, .5)
    assert trimmed[1, 2] == 0
    assert trimmed[0, 3] == 0
    assert trimmed[2, 2] == 2 and trimmed[2, 3] == 2
    for j, (row, col) in enumerate(seeds, start=1):
        assert np.all(intensity[trimmed == j] >= .5 * intensity[row, col])


# ------------------------------------------------------------------ editing
def test_delete_blob_removes_it_and_nothing_else():
    peaks = [(60.0, 1.2, 50, 2e5), (120.0, 2.6, 51, 3e5), (180.0, 4.3, 52, 1.5e5)]
    run = make_run(peaks, phase=0.7, rate=50.0)
    f = fold(run, 4.0, 0.7)
    b = find_blobs(run, f, PeakSettings(snr_min=10), floor=1.0)
    order = np.argsort(b.t1)
    middle = int(order[1])
    survivors_t1 = sorted(b.t1[order[[0, 2]]])
    b2 = delete_blobs(b, [middle])
    assert len(b2) == 2
    np.testing.assert_allclose(sorted(b2.t1), survivors_t1)
    # The deleted footprint's pixels are gone, not left dangling under some id.
    assert int((b2.labels > 0).sum()) == int(b.n_pixels[order[[0, 2]]].sum())
    assert b2.labels.max() == 2


def test_merge_two_footprints_recovers_the_split_peak():
    # Two watershed footprints close enough together that they came from
    # one real, slightly bimodal-looking feature.
    peaks = [(100.0, 1.5, 60, 1.5e5), (114.0, 1.5, 60, 1.0e5)]
    run = make_run(peaks, phase=0.7, rate=50.0, mz=(58, 63), sigma1=2.0)
    f = fold(run, 4.0, 0.7)
    b = find_blobs(run, f, PeakSettings(snr_min=8), floor=1.0)
    assert len(b) == 2          # watershed genuinely split it in two
    total_volume = float(b.volume.sum())
    total_pixels = int(b.n_pixels.sum())
    merged = merge_blobs(b, [0, 1])
    assert len(merged) == 1
    assert merged.volume[0] == pytest.approx(total_volume)         # an exact sum: footprints don't overlap
    assert merged.n_pixels[0] == total_pixels
    # A volume-weighted position between the two original apexes.
    assert min(b.t1) < merged.t1[0] < max(b.t1)
    np.testing.assert_array_equal(merged.spectra[0], b.spectra.sum(axis=0))


def test_add_blob_recovers_a_peak_too_faint_for_automatic_detection():
    bright = (40.0, 1.0, 70, 5e5)
    faint = (120.0, 2.0, 71, 3e3)
    run = make_run([bright, faint], phase=0.7, rate=50.0, mz=(65, 75))
    f = fold(run, 4.0, 0.7)
    b = find_blobs(run, f, PeakSettings(snr_min=10), floor=50.0)
    assert len(b) == 1                                              # the faint one was missed
    b2 = add_blob_at(b, run, f, 120.0, 2.0)
    assert len(b2) == 2
    assert b2.volume[-1] == pytest.approx(3e3, rel=0.1)
    assert b2.t1[-1] == pytest.approx(120.0, abs=1.0)
    with pytest.raises(ValueError):
        add_blob_at(b2, run, f, 40.0, 1.0)                          # already part of the bright blob
