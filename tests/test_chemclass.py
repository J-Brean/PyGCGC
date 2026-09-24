"""Compound-class rules (core/chemclass.py and presets/chemical_classes.toml).

The reference spectra are typed from the familiar EI patterns of one
representative compound per class (approximate, not library-exact). They test
that each shipped class can fire on what it is meant to describe, not that
its thresholds are right for any real sample; see real_data/CLASS_VALIDATION.md
for that.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from pygcgc.core.chemclass import (ClassSettings, classify_blobs, classify_spectra, compute_features, load_rules,
                                    soft_spectra_for_blobs)
from pygcgc.core.fold import fold
from pygcgc.core.peaks import PeakSettings, make_frame
from pygcgc.core.run import Run
from scipy import sparse
from pygcgc.core.hydrocarbons import AlkaneHit, AlkaneSeries
from pygcgc.core.siloxane import CYCLIC_REFERENCE

MZ = np.arange(38, 301)

REFERENCE = {
    "alkane": {41: 45, 43: 90, 55: 35, 56: 15, 57: 100, 69: 18, 70: 14, 71: 65, 83: 9, 84: 6, 85: 40, 97: 4, 99: 10, 113: 4, 226: 5},
    "alkene": {41: 85, 43: 85, 55: 100, 56: 45, 57: 65, 69: 85, 70: 45, 71: 30, 83: 70, 84: 35, 85: 15, 97: 55, 98: 25, 111: 25, 224: 4},
    "cycloalkane": {41: 55, 55: 75, 56: 30, 57: 20, 67: 20, 69: 30, 81: 15, 82: 40, 83: 100, 84: 8, 97: 5, 140: 12},
    "alkylbenzene": {39: 10, 51: 8, 65: 10, 77: 5, 78: 7, 91: 100, 92: 50, 105: 8, 134: 25, 135: 2},
    "indane / indene / tetralin": {39: 8, 51: 8, 63: 6, 65: 6, 77: 6, 91: 30, 103: 6, 115: 45, 116: 20, 117: 100, 118: 60, 119: 5},
    "naphthalene": {51: 6, 63: 5, 64: 8, 75: 5, 101: 2, 102: 8, 126: 3, 127: 10, 128: 100, 129: 11},
    "PAH": {76: 7, 88: 8, 89: 12, 150: 3, 151: 6, 152: 9, 176: 20, 177: 6, 178: 100, 179: 15},
    "N-heteroaromatic": {50: 5, 51: 6, 62: 5, 63: 8, 89: 35, 90: 33, 116: 9, 117: 100, 118: 9},          # indole
    "aliphatic amine": {44: 100, 58: 40, 72: 25, 86: 15, 43: 20, 57: 18, 41: 15, 55: 10, 199: 3},
    "alcohol / ether": {43: 40, 45: 100, 55: 20, 57: 12, 69: 12, 83: 10, 97: 8},                          # 2-octanol
    "ketone / aldehyde": {41: 20, 43: 100, 55: 12, 57: 22, 58: 65, 59: 25, 71: 22, 85: 8, 96: 8, 226: 3},
    "carboxylic acid": {41: 55, 43: 80, 55: 50, 57: 55, 60: 90, 69: 30, 73: 100, 83: 20, 97: 18, 129: 35, 256: 20},
    "ester": {41: 45, 43: 60, 55: 40, 57: 20, 59: 8, 74: 100, 75: 15, 87: 60, 97: 12, 143: 20, 227: 8, 239: 3, 270: 10},
    "aromatic oxygenate": {39: 20, 40: 10, 55: 6, 65: 25, 66: 40, 94: 100, 95: 7},                        # phenol
    "siloxane": dict(CYCLIC_REFERENCE),
    "phthalate": {57: 20, 71: 12, 104: 10, 113: 8, 149: 100, 150: 10, 167: 30, 279: 8},
}


def matrix(spectra: list) -> np.ndarray:
    out = np.zeros((len(spectra), MZ.size))
    for i, spec in enumerate(spectra):
        for m, v in spec.items():
            out[i, m - 38] = v
    return out


@pytest.fixture(scope="module")
def rules():
    return load_rules()


# ------------------------------------------------------------------ the shipped presets
def test_every_enabled_class_has_a_reference_and_fires(rules):
    enabled = [c.name for c in rules.classes if c.enabled]
    assert set(enabled) == set(REFERENCE), "add a reference spectrum for every enabled class"
    res = classify_spectra(matrix(list(REFERENCE.values())), MZ, rules)
    got = {name: call.cls for name, call in zip(REFERENCE, res.calls)}
    assert got == {name: name for name in REFERENCE}


def test_sulfur_class_ships_disabled(rules):
    sulfur = next(c for c in rules.classes if c.name == "S-heteroaromatic")
    assert not sulfur.enabled
    # ... and a disabled class never scores, even on a spectrum that satisfies it.
    res = classify_spectra(matrix([{134: 100, 135: 9, 136: 5, 89: 15, 63: 8, 108: 6}]), MZ, rules)
    assert res.calls[0].cls != "S-heteroaromatic"
    assert "S-heteroaromatic" not in [c.name for c in res.classes]


def test_indole_and_indane_are_separated_by_the_nitrogen_rule(rules):
    res = classify_spectra(matrix([REFERENCE["N-heteroaromatic"], REFERENCE["indane / indene / tetralin"]]), MZ, rules)
    assert [c.cls for c in res.calls] == ["N-heteroaromatic", "indane / indene / tetralin"]


def test_alkene_with_weak_molecular_ion_is_not_a_cycloalkane(rules):
    res = classify_spectra(matrix([REFERENCE["alkene"]]), MZ, rules)
    assert res.calls[0].cls == "alkene"


def test_unresolvable_class_pairs_are_flagged_not_forced(rules):
    """A ketone-like spectrum with a large m/z 44 (an amine's signature) must not come out as a confident ketone."""
    spec = {43: 30, 44: 60, 58: 100, 72: 40, 86: 20, 57: 15, 41: 15}
    call = classify_spectra(matrix([spec]), MZ, rules).calls[0]
    assert call.cls != "ketone / aldehyde"


def test_too_little_signal_is_unclassified(rules):
    call = classify_spectra(matrix([{57: 5, 43: 3, 71: 2}]), MZ, rules).calls[0]
    assert call.cls is None and call.tier == "unclassified"


# ------------------------------------------------------------------ features
def test_shares_and_missing_ions(rules, tmp_path):
    spec = matrix([{43: 40, 57: 40, 100: 5, 120: 5, 150: 5, 200: 5}])
    assert compute_features(spec, MZ, rules)["alkyl"][0] == pytest.approx(0.8)
    # m/z 30 is below this instrument's range: a set that names it is scored on the ions that exist.
    custom = load_rules(_write(tmp_path, '[ions]\nx = [30, 43]\n[[class]]\nname = "c"\nneed = ["x: 0.1 -> 0.2"]\n'))
    assert compute_features(spec, MZ, custom)["x"][0] == pytest.approx(0.4)


def test_molecular_ion_steps_down_from_the_isotope_and_ignores_background(rules):
    feats = compute_features(matrix([{128: 100, 129: 11, 127: 10, 102: 8, 51: 6}, {44: 100, 57: 12, 71: 9, 85: 6, 43: 9}]), MZ, rules)
    assert feats["M_mass"][0] == 128 and feats["M_odd"][0] == 0 and feats["M_rel"][0] == pytest.approx(1.0)
    assert feats["M_mass"][1] < 60 or np.isnan(feats["M_mass"][1]) or feats["M_mass"][1] != 44


def test_dt2_is_unknown_without_a_ladder_and_outside_it(rules):
    hits = {c: AlkaneHit(c, 60.0 * (10 + 5 * (c - 12)), 0.5 + 0.05 * (c - 12), 0.0, ()) for c in range(12, 20)}
    series = AlkaneSeries(hits)
    spec = matrix([REFERENCE["alkane"]] * 3)
    t1 = np.array([60.0 * 20, 60.0 * 20, 60.0 * 200])            # the third is far past C19
    t2 = np.array([series.t2_of(series.carbon_of(60.0 * 20)) + 1.5] * 2 + [3.0])
    none = compute_features(spec, MZ, rules, t1, t2, None)
    assert np.isnan(none["dt2"]).all()
    with_ladder = compute_features(spec, MZ, rules, t1, t2, series)
    assert with_ladder["dt2"][0] == pytest.approx(1.5, abs=1e-6)
    assert np.isnan(with_ladder["dt2"][2])


def test_unknown_evidence_is_skipped_not_failed(rules):
    """The alkane rule prefers dt2. With no ladder that preference is n/a and the score must not fall."""
    res = classify_spectra(matrix([REFERENCE["alkane"]]), MZ, rules)
    assert res.calls[0].score > 0.9
    assert "n/a" in res.explain(0)


# ------------------------------------------------------------------ result behaviour
def test_co_elution_lowers_the_tier_one_notch(rules):
    spec = matrix([REFERENCE["alkylbenzene"]] * 2)
    plain = classify_spectra(spec, MZ, rules, co_elution=[False, True])
    assert plain.calls[0].tier == "confident" and plain.calls[1].tier == "probable"


def test_explain_lists_every_criterion_of_the_winning_class(rules):
    res = classify_spectra(matrix([REFERENCE["PAH"]]), MZ, rules)
    text = res.explain(0)
    pah = next(c for c in rules.classes if c.name == "PAH")
    assert text.startswith("blob 0: PAH")
    assert all(c.text in text for c in pah.need + pah.prefer)


def test_summary_is_volume_weighted(rules):
    res = classify_spectra(matrix([REFERENCE["alkane"], REFERENCE["alkylbenzene"], {57: 1}]), MZ, rules)
    rows = {name: (n, pct) for name, n, _v, pct in res.summary(np.array([1.0, 3.0, 4.0]))}
    assert rows["alkane"][1] == pytest.approx(12.5) and rows["alkylbenzene"][1] == pytest.approx(37.5)
    assert rows["unclassified"][1] == pytest.approx(50.0)


def test_classify_blobs_falls_back_to_the_footprint_when_the_apex_is_empty(rules):
    footprint = matrix([REFERENCE["alkane"], REFERENCE["alkane"]])
    apex = np.vstack([footprint[0], np.zeros(MZ.size)])
    blobs = SimpleNamespace(mz=MZ, spectra=footprint, apex_spectra=apex, t1=np.array([100.0, 200.0]),
                            t2=np.array([1.0, 1.0]), co_elution=np.array([False, False]))
    res = classify_blobs(blobs, rules)
    assert [c.cls for c in res.calls] == ["alkane", "alkane"]


# ------------------------------------------------------------------ the rule file
def _write(tmp_path, text, name="rules.toml"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


GOOD = """
[ions]
a = [43, 57]
[[class]]
name = "x"
need = ["a: 0.2 -> 0.5"]
"""


def test_a_minimal_custom_rule_file_works(tmp_path):
    rs = load_rules(_write(tmp_path, GOOD))
    res = classify_spectra(matrix([{43: 50, 57: 50, 91: 1, 100: 1, 120: 1}]), MZ, rs)
    assert res.calls[0].cls == "x"


def test_a_utf8_bom_from_a_windows_editor_is_accepted(tmp_path):
    p = tmp_path / "bom.toml"
    p.write_bytes(b"\xef\xbb\xbf" + GOOD.encode())
    assert load_rules(p).classes[0].name == "x"


@pytest.mark.parametrize("bad, message", [
    (GOOD.replace("a: 0.2", "aa: 0.2"), "unknown feature 'aa'"),
    (GOOD.replace("a: 0.2 -> 0.5", "a: 0.2 => 0.5"), "can't read"),
    (GOOD.replace("0.2 -> 0.5", "0.5 -> 0.5"), "same start and end"),
    (GOOD.replace('need = ["a: 0.2 -> 0.5"]', "need = []"), "no 'need'"),
    (GOOD + '\n[[class]]\nname = "x"\nneed = ["a: 0.2 -> 0.5"]\n', "appears twice"),
    (GOOD + '\n[ratios]\nr = ["a", "nope"]\n', "needs two names"),
])
def test_rule_file_mistakes_are_reported_with_context(tmp_path, bad, message):
    with pytest.raises(ValueError, match=message):
        load_rules(_write(tmp_path, bad))


def test_the_unknown_feature_error_suggests_the_right_name(tmp_path):
    with pytest.raises(ValueError, match="Did you mean.*alkyl"):
        load_rules(_write(tmp_path, '[ions]\nalkyl = [43]\n[[class]]\nname = "x"\nneed = ["alky: 0.2 -> 0.5"]\n'))


def test_settings_round_trip():
    s = ClassSettings(confident=0.8)
    assert ClassSettings.from_dict(s.to_dict()) == s


# ------------------------------------------------------------------ the 14 eV channel
def test_soft_features_are_unknown_without_a_soft_spectrum(rules):
    feats = compute_features(matrix([REFERENCE["alkylbenzene"]]), MZ, rules)
    assert np.isnan(feats["soft_M_mass"]).all() and np.isnan(feats["soft_M_mod14"]).all()


def test_soft_molecular_ion_is_found_at_a_lower_threshold(rules):
    """An n-alkane M+ is only ~4% of the base peak at 14 eV: too weak for the 70 eV picker, found by the soft one."""
    soft = {57: 100, 71: 90, 85: 70, 99: 30, 113: 18, 127: 13, 226: 4, 225: 0.3}
    hard = matrix([REFERENCE["alkane"]])
    feats = compute_features(hard, MZ, rules, soft_spectra=matrix([soft]), soft_mz=MZ)
    assert feats["soft_M_mass"][0] == 226 and feats["soft_M_mod14"][0] == 2
    assert feats["M_mass"][0] != 226                          # the 70 eV feature is untouched by the soft spectrum


def test_soft_series_and_parity_change_the_score_only_through_prefer(rules):
    tol = matrix([{39: 10, 51: 8, 65: 10, 77: 5, 78: 7, 91: 100, 92: 50, 105: 8, 134: 25, 135: 2}])
    right = matrix([{91: 100, 105: 20, 119: 30, 134: 60, 135: 6}])          # M 134: 134 mod 14 = 8, an alkylbenzene
    wrong = matrix([{91: 100, 105: 20, 119: 30, 136: 60, 137: 6}])          # M 136: 136 mod 14 = 10
    a = classify_spectra(tol, MZ, rules, soft_spectra=right, soft_mz=MZ).calls[0]
    b = classify_spectra(tol, MZ, rules, soft_spectra=wrong, soft_mz=MZ).calls[0]
    none = classify_spectra(tol, MZ, rules).calls[0]
    assert a.cls == b.cls == none.cls == "alkylbenzene"
    assert a.score > b.score and a.score >= none.score - 1e-9


def test_equality_criteria_parse_and_score(tmp_path):
    p = _write(tmp_path, '[ions]\nx = [43]\n[[class]]\nname = "c"\nneed = ["x: 0.1 -> 0.2"]\nprefer = ["soft_M_mod14 == 8"]\n')
    assert load_rules(p).classes[0].prefer[0].op == "=="


def test_soft_spectra_for_blobs_use_the_same_pixels():
    times = 0.1 * np.arange(400)
    run = Run(times, np.arange(38, 301), sparse.csr_matrix(np.zeros((400, 263))), "hard")
    dense = np.zeros((400, 263))
    dense[:, 91 - 38] = 100.0
    soft_run = Run(times + 0.01, np.arange(38, 301), sparse.csr_matrix(dense), "soft")
    f = fold(run, 2.0, 0.0)
    frame = make_frame(f, PeakSettings())
    labels = np.zeros((frame.rows, f.n_cols), dtype=np.int32)
    labels[10:13, 3:5] = 1
    blobs = SimpleNamespace(settings=PeakSettings(), labels=labels, row=np.array([11]), col=np.array([3]))
    spectra, mz = soft_spectra_for_blobs(soft_run, f, blobs)
    r, c = frame.scan_pixels()
    n_scans = int(((c >= 0) & (c < labels.shape[1]) & (labels[r.clip(0, labels.shape[0] - 1), c.clip(0, labels.shape[1] - 1)] == 1)).sum())
    assert n_scans > 0 and spectra.shape == (1, 263) and spectra[0, 91 - 38] == pytest.approx(100.0 * n_scans)
    bad = SimpleNamespace(settings=PeakSettings(), labels=np.zeros((frame.rows + 3, f.n_cols), dtype=np.int32),
                          row=np.array([0]), col=np.array([0]))
    with pytest.raises(ValueError, match="same acquisition"):
        soft_spectra_for_blobs(soft_run, f, bad)