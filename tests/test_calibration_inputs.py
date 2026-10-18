"""Known amounts, ambiguous filenames and unit declarations cannot silently change truth."""

from dataclasses import replace

import pytest

from pygcgc.core.calibration_inputs import CalibrationInputs, filename_amount
from pygcgc.core.calibration_inputs import calibration_coverage_errors
from pygcgc.core.calibration import CalibrationSettings
from pygcgc.core.quant import QuantDecisions
from pygcgc.core.targets import Target


def valid_setup():
    return CalibrationInputs(
        targets=[Target('IS', role='istd', quant_ion=50, t1=60, t2=1),
                 Target('A', quant_ion=51, t1=100, t2=2, istd='IS')],
        levels={'L1': {'A': 1, 'IS': 1}, 'L2': {'A': 5, 'IS': 1}},
        runs={'cal1': {'kind': 'calibration', 'level': 'L1'},
              'cal2': {'kind': 'calibration', 'level': 'L2'}, 'sample': {'kind': 'sample'}},
        internal_standard_amounts={'IS': 1}, confirmed=True)


@pytest.mark.parametrize('filename,unit,amount', [
    ('Cal_L4_1ng_on_tube', 'ng', 1), ('Cal_0.25ng_on_tube', 'pg', 250),
    ('Cal_25ng_on_tube', 'ug', .025), ('Cal_5ug_on_tube', 'ng', 5000)])
def test_filename_suggestions_convert_explicit_mass_units(filename, unit, amount):
    result, note = filename_amount(filename, unit)
    assert result == pytest.approx(amount)
    assert 'confirm' in note


@pytest.mark.parametrize('filename', ['Cal_L2_025ng_on_tube', '1ng_5ng_mix', 'unknown', 'mix_1e-3ng', 'mix_-1ng'])
def test_ambiguous_or_missing_filename_amount_stays_unknown(filename):
    assert filename_amount(filename)[0] is None


def test_setup_requires_confirmed_actual_is_amount_and_real_calibration_range():
    setup = valid_setup()
    assert not setup.errors()
    assert any('Confirm' in e for e in replace(setup, confirmed=False).errors())
    assert any('constant across' in e for e in replace(setup, internal_standard_amounts={}).errors())
    assert any('distinct' in e for e in replace(setup, levels={'L1': {'A': 1}, 'L2': {'A': 1}}).errors())
    assert any('internal standard from' in e for e in replace(setup,
        targets=[setup.targets[0], replace(setup.targets[1], istd='unknown')]).errors())


@pytest.mark.parametrize('value', [float('nan'), float('inf'), -1])
def test_invalid_amounts_do_not_pass_setup(value):
    setup = valid_setup()
    setup.levels['L1']['A'] = value
    assert setup.errors()


def test_duplicate_names_and_invalid_roles_are_rejected():
    setup = valid_setup()
    setup.targets.append(replace(setup.targets[1], compound='a', role='unknown'))
    assert any('unique' in e for e in setup.errors())
    assert any('valid role' in e for e in setup.errors())


def test_coverage_obeys_actual_curve_minimum_zero_policy_and_exclusions():
    setup = valid_setup()
    policy = CalibrationSettings()
    errors = calibration_coverage_errors(setup.targets, setup.levels, setup.runs, policy)
    assert errors and 'require 3' in errors[0]
    setup.runs['cal0'] = {'kind': 'calibration', 'level': 'L0'}
    setup.levels['L0'] = {'A': 0}
    assert calibration_coverage_errors(setup.targets, setup.levels, setup.runs, policy)
    assert not calibration_coverage_errors(setup.targets, setup.levels, setup.runs, replace(policy, include_zero=True))
    decisions = QuantDecisions(excluded={'A': {'cal2': 'wrong peak'}})
    assert calibration_coverage_errors(setup.targets, setup.levels, setup.runs,
                                       replace(policy, include_zero=True), decisions)
    assert not calibration_coverage_errors(setup.targets, setup.levels, setup.runs,
                                          replace(policy, min_levels=2))
