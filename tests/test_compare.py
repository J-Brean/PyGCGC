"""⚖️ core/compare.py: resampling between runs, spectral similarity, peak pairing."""

from dataclasses import dataclass

import numpy as np
import pytest

from pygcgc.core.align import AnchorPoint, fit_alignment
from pygcgc.core.compare import (PeakPair, attach_labels, class_shares, match_peaks, pair_mapping, pixel_spectrum,
                                 resample, similarity, similarity_map, source_pixels, summarise_pairs,
                                 write_pairs_csv)
from pygcgc.core.fold import Fold, FoldSpec
from tests.conftest import make_run

PEAKS = [(40.0, 1.1, 50, 2e5), (90.0, 1.9, 51, 3e5), (150.0, 2.7, 52, 2e5)]


@dataclass
class Peaks:
    t1: np.ndarray
    t2: np.ndarray
    volume: np.ndarray
    mz: np.ndarray
    spectra: np.ndarray

    def __len__(self):
        return self.t1.size


def peaks(rows, mz=(50, 54)):
    """rows of (t1, t2, volume, ion): one single-ion spectrum per peak."""
    axis = np.arange(mz[0], mz[1] + 1)
    spectra = np.zeros((len(rows), axis.size))
    for k, (_t1, _t2, _v, ion) in enumerate(rows):
        spectra[k, ion - mz[0]] = 1.0
    return Peaks(np.array([r[0] for r in rows], float), np.array([r[1] for r in rows], float),
                 np.array([r[2] for r in rows], float), axis, spectra)


def fold_of(run):
    return Fold(run.times, FoldSpec(period=4.0))


def test_identity_mapping_resamples_a_run_onto_itself():
    run = make_run(PEAKS)
    f = fold_of(run)
    img = f.tic_image(run)
    out = resample(img, source_pixels(f, f, pair_mapping(None, None)))
    assert np.allclose(out, img)


def test_mapping_follows_both_transforms_through_the_reference():
    # The same peak is at 100 s in A, 103 s in the reference and 105 s in B.
    tf_a = fit_alignment([AnchorPoint("P", 103.0, 1.0, 100.0, 1.0)], "shift")
    tf_b = fit_alignment([AnchorPoint("P", 103.0, 1.0, 105.0, 1.0)], "shift")
    t1, t2 = pair_mapping(tf_a, tf_b)(np.array([100.0]), np.array([1.0]))
    assert t1[0] == pytest.approx(105.0)
    assert t2[0] == pytest.approx(1.0)


def test_similarity_bounds():
    mz = np.arange(50, 55)
    u = np.array([1.0, 0, 5, 0, 2])
    assert similarity(mz, u, 3 * u) == pytest.approx(1.0)     # amount does not matter
    assert similarity(mz, u, np.array([0, 1.0, 0, 1, 0])) == 0.0
    assert similarity(mz, u, np.zeros(5)) == 0.0


def test_similarity_map_is_one_on_shared_peaks_and_zero_where_the_chemistry_differs():
    a = make_run(PEAKS, seed=1)
    # B: the first two peaks at a tenth of the amount, the third swapped for another ion
    b = make_run([(40.0, 1.1, 50, 2e4), (90.0, 1.9, 51, 3e4), (150.0, 2.7, 54, 2e5)], seed=2)
    fa, fb = fold_of(a), fold_of(b)
    s = similarity_map(a, fa, b, fb, pair_mapping(None, None))
    for t1, t2, ion, _ in PEAKS:
        r, c = int(round(t2 / fa.row_pitch)), int((t1 - fa.t1[0]) // 4.0)
        expected = 0.0 if ion == 52 else 1.0
        assert s.similarity[r, c] == pytest.approx(expected, abs=1e-6)
    assert np.isnan(s.similarity[0, 0])                      # no signal in either run


def test_pixel_spectrum_reads_the_peak_ion():
    run = make_run(PEAKS)
    f = fold_of(run)
    col = int((90.0 - f.t1[0]) // 4.0)
    spec = pixel_spectrum(run, f, float(f.t1[col]), 1.9)
    assert run.mz[np.argmax(spec)] == 51


def test_match_pairs_shifted_peaks_and_reports_the_rest():
    a = peaks([(40.0, 1.1, 100.0, 50), (90.0, 1.9, 300.0, 51), (150.0, 2.7, 200.0, 52)])
    b = peaks([(41.5, 1.1, 200.0, 50), (91.5, 1.9, 150.0, 51), (151.5, 2.7, 200.0, 53),   # wrong ion
               (200.0, 1.0, 50.0, 54)])
    pairs = match_peaks(a, b, pair_mapping(None, None))
    both = [p for p in pairs if p.status == "both"]
    assert [(p.a, p.b) for p in both] == [(0, 0), (1, 1)]
    assert both[0].log2_ratio == pytest.approx(1.0)
    assert both[1].log2_ratio == pytest.approx(-1.0)
    assert both[0].dt1 == pytest.approx(1.5)
    assert both[0].how == "position and spectrum"
    s = summarise_pairs(pairs)
    assert (s["both"], s["only A"], s["only B"]) == (2, 1, 2)
    only_b = [p for p in pairs if p.status == "only B"]
    assert all(p.log2_ratio == np.inf for p in only_b)


def test_match_is_one_to_one_best_spectrum_first():
    # Two B peaks inside A's window: the one with A's spectrum wins even though it is further away.
    a = peaks([(90.0, 1.9, 100.0, 51)])
    b = peaks([(90.5, 1.9, 100.0, 52), (95.0, 1.9, 100.0, 51)])
    pairs = match_peaks(a, b, pair_mapping(None, None))
    pair = next(p for p in pairs if p.status == "both")
    assert pair.b == 1


def test_identities_pair_peaks_outside_the_windows():
    a = peaks([(40.0, 1.1, 100.0, 50), (90.0, 1.9, 100.0, 51)])
    b = peaks([(70.0, 1.6, 100.0, 50), (90.0, 1.9, 100.0, 51)])   # 30 s away: no position match
    pairs = match_peaks(a, b, pair_mapping(None, None), ["toluene-d8", ""], ["toluene-d8", ""])
    named = next(p for p in pairs if p.name == "toluene-d8")
    assert (named.a, named.b, named.how) == (0, 0, "name")
    assert named.similarity == pytest.approx(1.0)
    assert summarise_pairs(pairs)["both"] == 2


def test_attach_labels_hangs_each_label_on_the_nearest_peak():
    t1, t2 = np.array([40.0, 90.0, 94.0]), np.array([1.1, 1.9, 1.9])
    out = attach_labels(t1, t2, [(90.5, 1.9, "x"), (93.9, 1.92, "y"), (300.0, 1.0, "far"), (40.0, 1.1, "z"),
                                 (40.1, 1.1, "z2")], 4.0, 0.1)
    assert out == ["z; z2", "x", "y"]


def test_class_shares_joins_both_runs():
    rows = class_shares([("alkane", 3, 60.0, 60.0), ("PAH", 1, 40.0, 40.0)],
                        [("alkane", 2, 10.0, 30.0), ("siloxane", 4, 40.0, 70.0)])
    by = {r[0]: r for r in rows}
    assert by["alkane"] == ("alkane", 3, 2, 60.0, 30.0, -30.0)
    assert by["PAH"][2] == 0 and by["siloxane"][1] == 0
    assert [r[0] for r in rows] == ["alkane", "siloxane", "PAH"]     # largest combined share first


def test_pairs_csv(tmp_path):
    pairs = [PeakPair(a=0, b=1, t1_a=60.0, t2_a=1.0, t1_b=61.0, t2_b=1.0, dt1=1.0, dt2=0.0, similarity=0.9,
                      volume_a=10.0, volume_b=20.0, name_a="x", how="name"),
             PeakPair(a=None, b=2, t1_b=120.0, t2_b=2.0, volume_b=5.0)]
    path = write_pairs_csv(pairs, tmp_path / "p.csv", "RUN_A", "RUN_B")
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "# A = RUN_A; B = RUN_B"
    assert lines[2].startswith("both,name,x,")
    assert lines[2].endswith(",1.000")
    assert lines[3].startswith("only B,") and lines[3].endswith(",inf")
