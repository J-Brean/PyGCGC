"""Known-answer checks for failures exposed by the calibration chromatograms."""

from types import SimpleNamespace

import numpy as np
import pytest

from pygcgc.core.components import ComponentDetectionSettings
from pygcgc.core.modelselect import select_components
from pygcgc.core.overlap import (BlobCube, _emg_1d, _shape, fit_shapes_fast,
                                  predict_overlap, refit_spectra)


def test_old_project_uses_corrected_fitting_defaults_without_losing_custom_budgets():
    upgraded = ComponentDetectionSettings.from_dict({'fast_max_nfev': 24})
    assert upgraded.adaptive_shapes and upgraded.fast_max_nfev == 60
    assert ComponentDetectionSettings.from_dict({'fast_max_nfev': 120}).fast_max_nfev == 120
    explicit = ComponentDetectionSettings.from_dict({'adaptive_shapes': False, 'fast_max_nfev': 24})
    assert not explicit.adaptive_shapes and explicit.fast_max_nfev == 24


def test_emg_rendering_is_independent_of_which_pixels_are_requested():
    x = np.arange(40.0)
    whole = _emg_1d(x, 12.3, 2.1, 5.0)
    subset = np.array([0, 9, 12, 20, 30, 39])  # excludes the mode
    np.testing.assert_allclose(_emg_1d(x[subset], 12.3, 2.1, 5.0), whole[subset],
                               rtol=1e-13)


def test_saved_measured_profile_is_preserved_for_downstream_spectrum_extraction():
    from pygcgc.core.components import Component, component_profiles_from_durable
    profile = np.array([[0., .02, 0.], [.1, .4, .1], [.03, .3, .02], [0., .03, 0.]])
    component = Component(id='c', group_id='g', t1=10., t2=.2, volume=100.,
                          mz=np.array([50., 51.]), spectrum=np.array([80., 20.]),
                          model='template', row0=1.5, col0=1., sigma_row=np.nan,
                          sigma_col=np.nan, profile=profile)
    restored = Component.from_dict(component.to_dict())
    rendered = component_profiles_from_durable(restored, profile.shape)
    np.testing.assert_allclose(rendered, profile)
    assert rendered.sum()*restored.spectrum.sum() == pytest.approx(component.volume)


def test_independent_widths_recover_unequal_overlapping_peak_amounts():
    rr, cc = np.indices((49, 11), dtype=float)
    spectra = np.array([[9000., 2000., 200., 60.], [80., 400., 3500., 1200.]])
    params = [[22., 4.2, 4.8, 1.1], [27., 4.7, 1.4, .55]]
    profiles = [_shape('gaussian', rr, cc, np.array(p)) for p in params]
    data = 25 + sum(p[:, :, None] * s for p, s in zip(profiles, spectra))
    cube = BlobCube(data, np.ones(rr.shape, bool), np.arange(50, 54), 0, 0, .02, 10., 0.)
    hints = [SimpleNamespace(row=p[0], col=p[1]) for p in params]
    shaped = fit_shapes_fast(cube, 2, apex_hint=hints, max_nfev=60,
                             max_sigma_row=7, max_sigma_col=2.5,
                             background_order=1)
    fit = refit_spectra(shaped, cube)
    truth = [p.sum() * s.sum() for p, s in zip(profiles, spectra)]
    assert fit.converged
    assert [c.volume for c in fit.components] == pytest.approx(truth, rel=.03)
    assert [c.sigma_row for c in fit.components] == pytest.approx([4.8, 1.4], rel=.03)


def test_tilted_peak_is_fitted_as_one_component_with_its_full_amount():
    rr, cc = np.indices((55, 11), dtype=float)
    params = np.array([27., 5., 4.2, .9, -.78])
    profile = _shape('tilt', rr, cc, params)
    spectrum = np.array([7000., 2500., 200., 50., 400., 900., 60., 100.])
    data = 20 + profile[:, :, None] * spectrum
    cube = BlobCube(data, np.ones(rr.shape, bool), np.arange(50, 58), 0, 0, .02, 10., 0.)
    seeds = [SimpleNamespace(row=27., col=5.), SimpleNamespace(row=23., col=6.)]
    selected = select_components(cube, cube, seeds, ComponentDetectionSettings())
    assert selected.fit.k == 1
    assert selected.fit.model == 'tilt'
    assert selected.fit.components[0].volume == pytest.approx(profile.sum()*spectrum.sum(), rel=.02)
    assert selected.fit.converged
    assert np.linalg.norm(predict_overlap(selected.fit, cube)-data)/np.linalg.norm(data) < .02


@pytest.mark.parametrize('model', ['gaussian', 'tilt'])
def test_projected_jacobian_matches_finite_differences_with_signed_background(model):
    from scipy.optimize._numdiff import approx_derivative
    from pygcgc.core.overlap import (_gaussian_shape_derivatives,
                                     _projected_shape_jacobian, _solve_linear)
    rr, cc = np.indices((19, 7), dtype=float)
    rr, cc = rr.ravel(), cc.ravel()
    parameters = np.array([[8., 2.8, 2.2, .8], [11., 3.2, 1.4, .6]])
    if model == 'tilt':
        parameters = np.column_stack((parameters, [-.4, .2]))
    spectrum = np.array([[250., 80., 0.], [60., 0., 200.]])
    background = np.column_stack((np.ones(len(rr)), rr/rr.max()))
    profiles = np.column_stack([_shape(model, rr, cc, p) for p in parameters])
    rng = np.random.default_rng(52)
    data = profiles@spectrum + background@np.array([[-2., 4., 5.], [3., -1., 2.]])
    data += rng.normal(0, .5, data.shape)
    weights = 1/np.sqrt(np.maximum(data, 1))

    def residual(theta):
        profiles = np.column_stack([_shape(model, rr, cc, p)
                                    for p in theta.reshape(parameters.shape)])
        design = np.column_stack((profiles, background))
        coef = _solve_linear(design, data, True, weights, n_nonnegative=2)
        return ((design@coef-data)*weights).ravel()

    design = np.column_stack((profiles, background))
    coefficients = _solve_linear(design, data, True, weights, n_nonnegative=2)
    analytic = _projected_shape_jacobian(
        design, coefficients, weights, residual(parameters.ravel()).reshape(data.shape),
        _gaussian_shape_derivatives(parameters, profiles, rr, cc))
    numerical = approx_derivative(residual, parameters.ravel(), method='3-point')
    np.testing.assert_allclose(analytic, numerical, rtol=2e-4, atol=2e-5)
