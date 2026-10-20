"""Subtracting a blank run peak by peak, and from quantified amounts."""

from dataclasses import dataclass, field, replace

import numpy as np
import pytest

from pygcgc import QuantSettings, SequenceRun, fold, quantify_prepared
from pygcgc.core.peak_blank import (ABOVE, NEAR, REDUCED, REMOVED, UNTOUCHED, PeakBlankSettings, blank_status,
                                    corrected_spectra, subtract_blank_peaks, suggest_factor, summarise)
from pygcgc.core.quant import SampleResult, blank_correct_amounts
from pygcgc.core.targets import Target
from tests.conftest import CAL_LEVELS, CAL_TARGETS, GUI_PHASE, SAMPLE_AMOUNTS, make_cal_run
from tests.test_quant import _sequence


@dataclass
class Table:
    """A peak table shaped like a BlobSet, from (t1, t2, volume, {mz: intensity}) rows."""
    peaks: list
    mz: np.ndarray = field(default_factory=lambda: np.arange(40, 220))

    def __post_init__(self):
        self.t1 = np.array([p[0] for p in self.peaks], dtype=float)
        self.t2 = np.array([p[1] for p in self.peaks], dtype=float)
        self.volume = np.array([p[2] for p in self.peaks], dtype=float)
        self.spectra = np.zeros((len(self.peaks), self.mz.size))
        for i, p in enumerate(self.peaks):
            for ion, height in p[3].items():
                self.spectra[i, ion - 40] = height
        self.apex_spectra = np.zeros((0, 0))

    def __len__(self):
        return len(self.peaks)


SILOXANE, PLASTICISER, ALKANE = {73: 100, 147: 60, 207: 40}, {91: 100, 149: 80, 167: 30}, {57: 100, 71: 70, 85: 50}


# ------------------------------------------------------------------ peak by peak
def test_a_shared_contaminant_loses_the_blanks_volume_times_the_factor():
    sample = Table([(100.0, 1.0, 100.0, SILOXANE), (200.0, 1.5, 500.0, ALKANE)])
    blank = Table([(101.0, 1.02, 40.0, SILOXANE)])
    full = subtract_blank_peaks(sample, blank, PeakBlankSettings(factor=1.0))
    half = subtract_blank_peaks(sample, blank, PeakBlankSettings(factor=0.5))
    assert full[0].corrected == pytest.approx(60.0) and full[0].removed == pytest.approx(40.0)
    assert half[0].corrected == pytest.approx(80.0)
    assert full[0].match == REDUCED and full[0].blank_index == 0
    assert full[1].match == UNTOUCHED and full[1].corrected == 500.0 and full[1].blank_index is None


def test_a_peak_is_never_taken_below_zero():
    sample = Table([(100.0, 1.0, 30.0, SILOXANE)])
    blank = Table([(100.0, 1.0, 40.0, SILOXANE)])
    row = subtract_blank_peaks(sample, blank, PeakBlankSettings(factor=3.0))[0]
    assert row.corrected == 0.0 and row.removed == 30.0 and row.match == REMOVED and row.fraction_removed == 1.0


def test_factor_zero_changes_nothing():
    sample = Table([(100.0, 1.0, 100.0, SILOXANE)])
    row = subtract_blank_peaks(sample, Table([(100.0, 1.0, 40.0, SILOXANE)]), PeakBlankSettings(factor=0.0))[0]
    assert row.corrected == 100.0 and row.match == UNTOUCHED


def test_a_neighbour_with_a_different_spectrum_is_not_the_blanks_peak():
    sample = Table([(100.0, 1.0, 100.0, ALKANE)])
    blank = Table([(100.5, 1.01, 40.0, SILOXANE)])                    # in the window, but a different compound
    assert subtract_blank_peaks(sample, blank, PeakBlankSettings())[0].blank_index is None


def test_a_peak_outside_the_windows_is_not_matched():
    sample = Table([(100.0, 1.0, 100.0, SILOXANE)])
    far_t1, far_t2 = Table([(140.0, 1.0, 40.0, SILOXANE)]), Table([(100.0, 1.6, 40.0, SILOXANE)])
    assert subtract_blank_peaks(sample, far_t1, PeakBlankSettings())[0].blank_index is None
    assert subtract_blank_peaks(sample, far_t2, PeakBlankSettings())[0].blank_index is None
    assert subtract_blank_peaks(sample, far_t1, PeakBlankSettings(max_dt1=60.0))[0].blank_index == 0


def test_one_blank_peak_is_used_once_by_the_better_matching_sample_peak():
    close_but_different = {73: 100, 147: 20, 207: 10, 91: 40}
    sample = Table([(100.0, 1.0, 100.0, close_but_different), (101.0, 1.0, 100.0, SILOXANE)])
    blank = Table([(100.5, 1.0, 40.0, SILOXANE)])
    rows = subtract_blank_peaks(sample, blank, PeakBlankSettings(min_similarity=0.5))
    assert rows[1].blank_index == 0 and rows[0].blank_index is None


def test_the_alignment_mapping_is_used_to_read_the_blank_at_the_samples_times():
    sample = Table([(100.0, 1.0, 100.0, SILOXANE)])
    blank = Table([(160.0, 1.0, 40.0, SILOXANE)])                     # the blank ran 60 s later
    assert subtract_blank_peaks(sample, blank, PeakBlankSettings())[0].blank_index is None
    rows = subtract_blank_peaks(sample, blank, PeakBlankSettings(), lambda t1, t2: (t1 - 60.0, t2))
    assert rows[0].blank_index == 0 and rows[0].dt1 == pytest.approx(0.0)


def test_peaks_are_labelled_above_blank_or_near_blank_by_the_ratio_to_the_blanks_contribution():
    sample = Table([(100.0, 1.0, 400.0, SILOXANE), (200.0, 1.5, 100.0, PLASTICISER), (300.0, 2.0, 80.0, ALKANE)])
    blank = Table([(100.0, 1.0, 40.0, SILOXANE), (200.0, 1.5, 60.0, PLASTICISER)])
    s = PeakBlankSettings(factor=1.0, significance=3.0)
    rows = subtract_blank_peaks(sample, blank, s)
    assert [blank_status(r, s) for r in rows] == [ABOVE, NEAR, ABOVE]       # 10x, 1.7x, and not in the blank at all
    summary = summarise(rows, s)
    assert summary["above"] == 2 and summary["near"] == 1 and summary["matched"] == 2 and summary["untouched"] == 1
    assert summary["signal_removed"] == pytest.approx(100.0 / 580.0)


def test_the_suggested_factor_is_the_median_ratio_and_is_never_applied():
    sample = Table([(100.0, 1.0, 100.0, SILOXANE), (200.0, 1.5, 300.0, PLASTICISER), (300.0, 2.0, 80.0, ALKANE)])
    blank = Table([(100.0, 1.0, 50.0, SILOXANE), (200.0, 1.5, 100.0, PLASTICISER)])
    rows = subtract_blank_peaks(sample, blank, PeakBlankSettings())
    hint = suggest_factor(rows)
    assert hint["n"] == 2 and hint["median"] == pytest.approx(2.5) and hint["q1"] < hint["median"] < hint["q3"]
    assert rows[0].corrected == pytest.approx(50.0)                             # the default factor of 1 was used
    assert suggest_factor(subtract_blank_peaks(sample, Table([]), PeakBlankSettings()))["n"] == 0


def test_corrected_spectra_keep_their_shape_and_scale_to_the_corrected_volume():
    sample = Table([(100.0, 1.0, 100.0, SILOXANE)])
    rows = subtract_blank_peaks(sample, Table([(100.0, 1.0, 40.0, SILOXANE)]), PeakBlankSettings())
    scaled = corrected_spectra(sample, rows)
    assert scaled[0, 73 - 40] == pytest.approx(60.0) and scaled[0, 147 - 40] == pytest.approx(36.0)


def test_an_empty_blank_leaves_every_peak_untouched():
    sample = Table([(100.0, 1.0, 100.0, SILOXANE)])
    assert subtract_blank_peaks(sample, None, PeakBlankSettings())[0].corrected == 100.0
    assert subtract_blank_peaks(Table([]), Table([(1.0, 1.0, 1.0, SILOXANE)]), PeakBlankSettings()) == []


def test_settings_round_trip():
    s = PeakBlankSettings(blank_run="B1", factor=0.4, max_dt1=9.0)
    assert PeakBlankSettings.from_dict(s.to_dict()) == s


# ------------------------------------------------------------------ amounts
def result(run, compound, amount, se=0.1, role="target", status="ok"):
    return SampleResult(run, compound, role, True, status, [], 100.0, 1.0, 1.0, amount, se, amount / 2.0, "")


def test_amounts_have_the_blanks_amount_times_the_factor_taken_away():
    seq = [SequenceRun("CAL_1", "calibration", "L1"), SequenceRun("BLANK", "blank"), SequenceRun("S1", "sample", volume=2.0)]
    rows = [result("BLANK", "A", 1.0, 0.05), result("S1", "A", 3.0, 0.2), result("S1", "B", 5.0), result("S1", "IS", 1.0, role="istd"),
            result("CAL_1", "A", 2.0), result("BLANK", "B", float("nan"), status="not found")]
    out = {(r.run, r.compound): r for r in blank_correct_amounts(rows, seq, QuantSettings(blank_run="BLANK", blank_factor=0.5))}
    a = out[("S1", "A")]
    assert a.amount == pytest.approx(2.5) and a.amount_raw == 3.0 and a.blank_amount == 1.0
    assert a.amount_se == pytest.approx(np.hypot(0.2, 0.5 * 0.05)) and a.concentration == pytest.approx(1.25)
    assert a.blank_status == ABOVE and "blank 1 x 0.5 subtracted" in a.note
    assert out[("S1", "B")].amount == 5.0 and out[("S1", "B")].blank_amount == 0.0       # the blank did not show it
    assert out[("S1", "IS")].amount == 1.0 and out[("CAL_1", "A")].amount == 2.0          # standards left alone
    assert out[("BLANK", "A")].amount == 1.0 and out[("BLANK", "A")].amount_raw != out[("BLANK", "A")].amount_raw


def test_an_amount_at_or_below_the_blank_is_zero_and_flagged_not_negative():
    seq = [SequenceRun("BLANK", "blank"), SequenceRun("S1", "sample", volume=2.0)]
    rows = [result("BLANK", "A", 2.0), result("S1", "A", 1.5)]
    got = next(r for r in blank_correct_amounts(rows, seq, QuantSettings(blank_run="BLANK")) if r.run == "S1")
    assert got.amount == 0.0 and "at or below blank" in got.flags and got.status == "check"
    assert got.blank_status == NEAR and got.concentration == 0.0


def test_without_a_blank_run_nothing_changes():
    seq = [SequenceRun("S1", "sample")]
    rows = [result("S1", "A", 3.0)]
    assert blank_correct_amounts(rows, seq, QuantSettings()) is rows
    assert blank_correct_amounts(rows, seq, QuantSettings(blank_run="MISSING")) is rows


def test_blank_settings_are_part_of_the_quantification_settings_and_round_trip():
    s = QuantSettings(blank_run="B", blank_factor=0.7, blank_significance=5.0)
    again = QuantSettings.from_dict(s.to_dict())
    assert again.blank_run == "B" and again.blank_factor == 0.7 and again.blank_settings().significance == 5.0


# ------------------------------------------------------------------ end to end
def _sequence_with_blank(blank_amounts):
    prepared, info, targets, levels = _sequence()
    run = make_cal_run(blank_amounts, shift=(1.5, 0.02), seed=99, noise=0.0)
    prepared.append((run, run, fold(run, 4.0, GUI_PHASE), None))
    info.append(SequenceRun("BLANK", "blank", volume=2.0))
    return prepared, info, targets, levels


def test_quantified_amounts_are_sample_minus_blank_end_to_end():
    prepared, info, targets, levels = _sequence_with_blank({"A": 1.0, "B": 0.5, "C": 0.0})
    plain = quantify_prepared(prepared, info, targets, levels)
    corrected = quantify_prepared(prepared, info, targets, levels, settings=QuantSettings(blank_run="BLANK", blank_factor=1.0))
    half = quantify_prepared(prepared, info, targets, levels, settings=QuantSettings(blank_run="BLANK", blank_factor=0.5))
    blank_a = next(r for r in plain.results if r.run == "BLANK" and r.compound == "A")
    assert blank_a.amount == pytest.approx(1.0, rel=0.1)
    for sample in SAMPLE_AMOUNTS:
        for compound, blank_amount in (("A", 1.0), ("B", 0.5), ("C", 0.0)):
            raw = next(r for r in plain.results if r.run == sample and r.compound == compound)
            full = next(r for r in corrected.results if r.run == sample and r.compound == compound)
            part = next(r for r in half.results if r.run == sample and r.compound == compound)
            assert full.amount == pytest.approx(max(SAMPLE_AMOUNTS[sample][compound] - blank_amount, 0.0), rel=0.1, abs=0.1)
            assert part.amount == pytest.approx(max(SAMPLE_AMOUNTS[sample][compound] - 0.5 * blank_amount, 0.0), rel=0.1, abs=0.1)
            assert full.amount_raw == pytest.approx(raw.amount) and full.blank_status in (ABOVE, NEAR)
    assert next(r for r in corrected.results if r.run == "BLANK" and r.compound == "A").amount == blank_a.amount


def test_all_signals_use_the_peak_level_correction_before_reading_the_curve():
    prepared, info, targets, levels = _sequence_with_blank({"A": 0.5, "B": 0.25, "C": 0.0})
    plain = quantify_prepared(prepared, info, targets, levels)
    corrected = quantify_prepared(prepared, info, targets, levels, settings=QuantSettings(blank_run="BLANK"))
    raw = next(r for r in plain.signals if r.run == "S1" and r.target == "A")
    fixed = next(r for r in corrected.signals if r.run == "S1" and r.target == "A")
    assert fixed.volume == pytest.approx(raw.volume)                              # the raw volume is still reported
    assert fixed.blank_volume > 0 and fixed.corrected_volume == pytest.approx(raw.volume - fixed.blank_volume)
    assert fixed.amount == pytest.approx(SAMPLE_AMOUNTS["S1"]["A"] - 0.5, rel=0.1)
    assert fixed.blank_status == ABOVE and raw.blank_status == "" and np.isnan(raw.corrected_volume)
    # a spiked internal standard is in the blank-corrected sample on purpose: its peak is not reduced
    spiked = next(r for r in corrected.signals if r.run == "S1" and r.target == "IS")
    assert spiked.blank_status == "" and np.isnan(spiked.corrected_volume)
    # the blank's own and the calibration runs' rows are not corrected
    assert all(np.isnan(r.corrected_volume) for r in corrected.signals if r.run == "BLANK" or r.run.startswith("CAL_"))
    assert replace(QuantSettings(), blank_run="BLANK").blank_run == "BLANK" and CAL_LEVELS and CAL_TARGETS and Target
