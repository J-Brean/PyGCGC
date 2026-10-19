"""Changing the declared amount unit must not change the calibration physics."""

import numpy as np
import pytest
from pygcgc.core.calibration import CalibrationSettings, fit_calibration, _invert


@pytest.mark.parametrize("scale", [1., 1e3, 1e6, 1e9])
@pytest.mark.parametrize("model", ["linear", "quadratic"])
@pytest.mark.parametrize("weights", ["1", "1/x", "1/x2"])
@pytest.mark.parametrize("zero", [False, True])
def test_calibration_recovers_the_same_amounts_in_different_units(scale, model, weights, zero):
    x = np.array([.25, 1, 5, 10, 25.])
    y = (0 if zero else 2) + 3*x + (.2*x*x if model == "quadratic" else 0)
    fit = fit_calibration("known", x*scale, y,
                          settings=CalibrationSettings(model=model, weights=weights, force_zero=zero))
    np.testing.assert_allclose(fit.back_calc/scale, x, rtol=1e-9, atol=1e-10)
    assert fit.intercept == pytest.approx(0 if zero else 2, abs=1e-9)
    assert fit.slope*scale == pytest.approx(3, abs=1e-9)
    if model == "quadratic":
        assert fit.curvature*scale**2 == pytest.approx(.2, abs=1e-9)


def test_quadratic_inverse_keeps_precision_close_to_the_intercept():
    result = _invert(0., 1e8, 1., np.array([1., 10.]))
    np.testing.assert_allclose(1e8*result + result*result, [1., 10.], rtol=1e-14)


def test_small_curvature_coefficients_are_not_discarded_when_units_are_large():
    amount = 2e9
    response = 3e-9*amount + 2e-19*amount*amount
    assert float(_invert(0., 3e-9, 2e-19, response)) == pytest.approx(amount, rel=1e-12)


@pytest.mark.parametrize("scale", [1e3, 1e6, 1e9])
def test_inverse_prediction_uncertainty_scales_with_the_amount_unit(scale):
    x = np.array([.25, 1, 5, 10, 25.])
    y = 2+3*x+.2*x*x + np.array([.03, -.02, .02, -.05, .04])
    settings = CalibrationSettings(model="quadratic", weights="1/x")
    ref = fit_calibration("known", x, y, settings=settings)
    actual = fit_calibration("known", x*scale, y, settings=settings)
    amounts, se = actual.amount_with_se(np.array([5., 25.]))
    ref_amounts, ref_se = ref.amount_with_se(np.array([5., 25.]))
    np.testing.assert_allclose(amounts/scale, ref_amounts, rtol=1e-9)
    np.testing.assert_allclose(se/scale, ref_se, rtol=1e-8)
