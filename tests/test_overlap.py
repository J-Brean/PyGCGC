"""✂️ Overlap toolkit: resolving two compounds sharing one footprint."""

from dataclasses import dataclass

import numpy as np
import pytest

from pygcgc import Fold, FoldSpec
from pygcgc.core.overlap import (BlobCube, MCRResult, OverlapFit, OverlapSettings, PARAFACResult,
                                 UniqueIonMap, _unimodal_1d, blob_cube, build_run_tensor,
                                 assess_components, choose_k, component_maps, core_consistency,
                                 cross_validate_overlap, fit_overlap,
                                 fit_template_overlap,
                                 mcr_als, mcr_restarts, parafac, predict_overlap,
                                 split_by_line, split_by_valley, split_by_watershed,
                                 split_half_stability, split_volumes, unique_ion_apexes)
from pygcgc.core.peaks import BlobSet
from tests.conftest import make_run

PERIOD = 4.0
RATE = 50.0
ROW_PITCH = PERIOD / round(PERIOD * RATE)     # 0.02 s

# Two compounds landing in the same modulation column (t1 = 100 and 100.8,
# both // 4 == 25) but 5 rows (0.10 s) apart in t2, on different ions --
# exactly the "one footprint, several ions that don't quite agree" scenario
# the overlap toolkit is for.
PEAK_A = (100.0, 1.50, 50, 1e5)
PEAK_B = (100.8, 1.60, 51, 1e5)


def _hand_blob(run, fold, mask) -> BlobSet:
    """A `BlobSet` with exactly one blob, its footprint given directly
    (bypassing `find_blobs`) so tests control the overlap scenario
    precisely rather than fighting automatic watershed segmentation."""
    cube_all = fold.cube(run, dtype=np.float64)
    spectrum = cube_all[mask].sum(axis=0)
    labels = np.where(mask, 1, 0).astype(np.int32)
    tic = cube_all.sum(axis=2)
    masked = np.where(mask, tic, -np.inf)
    r, c = np.unravel_index(int(np.argmax(masked)), tic.shape)
    rows_idx, cols_idx = np.nonzero(mask)
    return BlobSet(mz=run.mz.copy(), t1=np.array([float(fold.t1[c])]), t2=np.array([r * fold.row_pitch]),
                  col=np.array([c]), row=np.array([r]), height=np.array([float(tic[r, c])]),
                  snr=np.array([10.0]), volume=np.array([float(spectrum.sum())]),
                  n_pixels=np.array([int(mask.sum())]), n_modulations=np.array([int(np.unique(cols_idx).size)]),
                  spectra=spectrum.reshape(1, -1), labels=labels, wrap_rows=0,
                  pitch=fold.row_pitch, period=fold.spec.period, t1_start=float(fold.t1[0]))


@pytest.fixture
def overlap_scene():
    run = make_run([PEAK_A, PEAK_B], period=PERIOD, phase=0.0, rate=RATE, duration=130.0,
                   mz=(48, 55), sigma1=4.0, sigma2=0.06, noise=0.0, seed=0)
    fold = Fold(run.times, FoldSpec(period=PERIOD, phase=0.0))
    mask = np.zeros(fold.shape, dtype=bool)
    mask[60:100, 20:31] = True
    blobs = _hand_blob(run, fold, mask)
    cube = blob_cube(run, fold, blobs, 0, OverlapSettings(ion_share_min=0.01, min_ions=2))
    return run, fold, blobs, cube


# ------------------------------------------------------------------ blob_cube
def test_blob_cube_selects_the_blobs_own_significant_ions(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    assert isinstance(cube, BlobCube)
    assert 50 in cube.ions and 51 in cube.ions
    assert cube.data.shape[:2] == cube.mask.shape
    assert cube.mask.sum() == 40 * 11
    assert cube.data.shape[2] == cube.ions.size


def test_blob_cube_time_conversion_matches_the_fold(overlap_scene):
    _run, fold, _blobs, cube = overlap_scene
    assert cube.t1_of(0) == pytest.approx(fold.t1[cube.col0])
    assert cube.t2_of(0) == pytest.approx(cube.row0 * fold.row_pitch)


# ------------------------------------------------------------- 1. unique ions
def test_unique_ion_apexes_separates_the_two_ions(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    result = unique_ion_apexes(cube)
    assert isinstance(result, UniqueIonMap)
    assert result.n_components == 2
    by_ion = {a.ion: a for a in result.apexes}
    assert by_ion[50.0].cluster != by_ion[51.0].cluster
    # ion 51's own apex should sit a few rows below ion 50's (t2 = 1.60 vs 1.50).
    assert by_ion[51.0].row > by_ion[50.0].row


def test_unique_ion_apexes_reports_one_cluster_for_a_single_compound():
    run = make_run([PEAK_A], period=PERIOD, phase=0.0, rate=RATE, duration=130.0, mz=(48, 55),
                   sigma1=4.0, sigma2=0.06, noise=0.0, seed=1)
    fold = Fold(run.times, FoldSpec(period=PERIOD, phase=0.0))
    mask = np.zeros(fold.shape, dtype=bool)
    mask[60:100, 20:31] = True
    blobs = _hand_blob(run, fold, mask)
    cube = blob_cube(run, fold, blobs, 0, OverlapSettings(ion_share_min=0.01, min_ions=2))
    result = unique_ion_apexes(cube)
    assert result.n_components == 1


# ------------------------------------------------------------- 2. geometric splits
def _apex_pixels(cube):
    apexes = {a.ion: a for a in unique_ion_apexes(cube).apexes}
    return (apexes[50.0].row, apexes[50.0].col), (apexes[51.0].row, apexes[51.0].col)


def test_split_by_line_and_split_volumes_separate_by_ion(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    (r50, c50), (r51, c51) = _apex_pixels(cube)
    r_mid = (r50 + r51) / 2.0
    labels = split_by_line(cube, (r_mid, 0), (r_mid, cube.mask.shape[1] - 1))
    assert set(np.unique(labels)) <= {0, 1, 2}
    totals = split_volumes(cube, labels)
    assert len(totals) == 2
    ion50_idx = int(np.flatnonzero(cube.ions == 50)[0])
    ion51_idx = int(np.flatnonzero(cube.ions == 51)[0])
    # Ion 50 sits at the lower row (t2 = 1.50 < 1.60): whichever region has
    # the lower mean row should be the one its spectrum dominates.
    lower_label = min(totals, key=lambda lab: totals[lab]["t2"])
    upper_label = next(lab for lab in totals if lab != lower_label)
    assert totals[lower_label]["spectrum"][ion50_idx] > totals[lower_label]["spectrum"][ion51_idx]
    assert totals[upper_label]["spectrum"][ion51_idx] > totals[upper_label]["spectrum"][ion50_idx]


def test_split_by_watershed_seeded_at_each_apex_separates_the_ions(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    seeds = list(_apex_pixels(cube))
    labels = split_by_watershed(cube, seeds)
    assert set(np.unique(labels[cube.mask])) == {1, 2}
    totals = split_volumes(cube, labels)
    ion50_idx = int(np.flatnonzero(cube.ions == 50)[0])
    ion51_idx = int(np.flatnonzero(cube.ions == 51)[0])
    label_a, label_b = totals.keys()
    dominant_ion = {label_a: int(np.argmax(totals[label_a]["spectrum"])),
                    label_b: int(np.argmax(totals[label_b]["spectrum"]))}
    assert {dominant_ion[label_a], dominant_ion[label_b]} == {ion50_idx, ion51_idx}


def test_split_by_valley_stays_inside_the_footprint_and_joins_the_endpoints(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    p1, p2 = _apex_pixels(cube)
    rows, cols = split_by_valley(cube, p1, p2)
    assert cube.mask[rows, cols].all()
    assert (rows[0], cols[0]) == (int(round(p1[0])), int(round(p1[1])))
    assert (rows[-1], cols[-1]) == (int(round(p2[0])), int(round(p2[1])))


# ------------------------------------------------------------- 3. peak-shape fit
def test_fit_overlap_recovers_two_components_near_the_true_positions(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    (r50, c50), (r51, c51) = _apex_pixels(cube)
    fit = fit_overlap(cube, k=2, model="gaussian")
    assert fit.weighted_rss >= 0
    assert fit.n_parameters == 2 * 4 + 2 * cube.data.shape[2]
    assert isinstance(fit, OverlapFit)
    assert len(fit.components) == 2
    rows_found = sorted(comp.row0 for comp in fit.components)
    rows_true = sorted([r50, r51])
    assert rows_found[0] == pytest.approx(rows_true[0], abs=2.0)
    assert rows_found[1] == pytest.approx(rows_true[1], abs=2.0)
    for comp in fit.components:
        assert np.isfinite(comp.row0_se) and np.isfinite(comp.col0_se)
        assert comp.volume > 0


def test_every_peak_model_fits_and_reports_its_own_parameters(overlap_scene):
    from pygcgc.core.overlap import MODELS, _N_SHAPE
    _run, _fold, _blobs, cube = overlap_scene
    for model in MODELS:
        fit = fit_overlap(cube, k=1, model=model)
        assert fit.n_parameters == _N_SHAPE[model] + cube.data.shape[2], model
        c = fit.components[0]
        assert np.isfinite(c.volume) and c.volume > 0, model
        assert (c.tau is not None) == (model in ("emg", "emg2"))
        assert (c.tau_row is not None) == (model in ("emg_t2", "emg2", "emg_t2_tilt", "emg_t2_cell"))
        assert (c.rho is not None) == (model in ("tilt", "emg_t2_tilt"))
        assert np.isfinite(predict_overlap(fit, cube)).all(), model
    with pytest.raises(ValueError):
        fit_overlap(cube, k=1, model="banana")


def _synthetic_cube(model, params, shape=(60, 11), ions=(50, 51)):
    from pygcgc.core.overlap import _shape
    rr, cc = np.mgrid[0:shape[0], 0:shape[1]].astype(float)
    g = _shape(model, rr.ravel(), cc.ravel(), np.asarray(params, dtype=float)).reshape(shape)
    data = np.stack([5000.0 * g, 2000.0 * g], axis=2) + 1.0
    return BlobCube(data=data, mask=np.ones(shape, bool), ions=np.array(ions), row0=0, col0=0,
                    pitch=0.02, period=10.0, t1_start=0.0)


def test_a_t2_tailing_peak_is_fitted_better_by_the_t2_tail_model():
    cube = _synthetic_cube("emg_t2", [25.0, 5.0, 2.0, 1.5, 6.0])
    plain = fit_overlap(cube, 1, "gaussian", n_starts=3)
    tailed = fit_overlap(cube, 1, "emg_t2", n_starts=3)
    assert tailed.weighted_rss < 0.2 * plain.weighted_rss
    assert tailed.components[0].tau_row == pytest.approx(6.0, rel=0.25)


def test_a_tilted_peak_is_fitted_better_by_the_tilt_model():
    cube = _synthetic_cube("tilt", [25.0, 5.0, 2.0, 1.5, 0.7])
    plain = fit_overlap(cube, 1, "gaussian", n_starts=3)
    tilted = fit_overlap(cube, 1, "tilt", n_starts=3)
    assert tilted.weighted_rss < 0.2 * plain.weighted_rss
    assert tilted.components[0].rho == pytest.approx(0.7, abs=0.15)


def test_the_cell_model_handles_a_peak_narrower_than_one_column():
    cube = _synthetic_cube("cell", [25.0, 5.0, 2.0, 0.3])
    fit = fit_overlap(cube, 1, "cell", n_starts=3)
    assert fit.components[0].sigma_col == pytest.approx(0.3, abs=0.1)


def test_choose_k_prefers_two_components_for_a_genuine_two_peak_footprint(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    best = choose_k(cube, max_k=3, model="gaussian")
    assert best.k == 2


def test_overlap_fit_can_estimate_a_nonnegative_local_background(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    raised = BlobCube(data=cube.data + 7.0, mask=cube.mask, ions=cube.ions,
                      row0=cube.row0, col0=cube.col0, pitch=cube.pitch,
                      period=cube.period, t1_start=cube.t1_start)
    fit = fit_overlap(raised, k=2, include_background=True, n_starts=2, seed=3)
    assert fit.n_starts == 2
    assert fit.background.shape == (cube.data.shape[2],)
    assert np.all(fit.background >= 0)
    assert np.median(fit.background) == pytest.approx(7.0, abs=1.0)


def test_overlap_fit_can_estimate_a_bilinear_local_background(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    rr, cc = np.indices(cube.mask.shape)
    background = 2.0 + 0.2 * rr + 0.4 * cc
    raised = BlobCube(data=cube.data + background[:, :, None], mask=cube.mask,
                      ions=cube.ions, row0=cube.row0, col0=cube.col0,
                      pitch=cube.pitch, period=cube.period, t1_start=cube.t1_start)
    fit = fit_overlap(raised, k=2, include_background=True, background_order=1,
                      n_starts=2, seed=3)
    assert fit.background.shape == (4, cube.data.shape[2])
    assert fit.background_surface.shape == cube.data.shape
    assert np.median(abs(fit.background_surface.sum(axis=2) -
                         background * cube.data.shape[2])) < 3.0


def test_overlap_background_can_follow_negative_cleaning_residuals(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    rr, cc = np.indices(cube.mask.shape)
    background = -6.0 + 0.03 * rr + 0.08 * cc
    shifted = BlobCube(data=cube.data + background[:, :, None], mask=cube.mask,
                       ions=cube.ions, row0=cube.row0, col0=cube.col0,
                       pitch=cube.pitch, period=cube.period,
                       t1_start=cube.t1_start)
    fit = fit_overlap(shifted, k=2, include_background=True, background_order=1,
                      n_starts=2, seed=8)
    assert fit.background.min() < 0
    reconstructed_background = fit.background_surface.sum(axis=2)
    assert np.median(abs(reconstructed_background -
                         background * cube.data.shape[2])) < 3.0


def test_overlap_prediction_and_spatial_cross_validation_favour_two_peaks(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    fit = fit_overlap(cube, 2)
    prediction = predict_overlap(fit, cube)
    assert prediction.shape == cube.data.shape
    one = cross_validate_overlap(cube, 1, folds=3, n_starts=2)
    two = cross_validate_overlap(cube, 2, folds=3, n_starts=2)
    assert two.mean_rmse < one.mean_rmse
    assert two.fold_rmse.shape == (3,)


def test_known_component_areas_are_recovered_with_a_sloping_background():
    from pygcgc.core.overlap import _shape

    shape = (50, 9)
    rr, cc = np.indices(shape, dtype=float)
    params = ([17.0, 3.0, 2.2, 0.8], [32.0, 5.0, 2.6, 0.9])
    spectra = np.asarray([[6000.0, 250.0, 100.0],
                          [150.0, 4200.0, 300.0]])
    profiles = [_shape("gaussian", rr.ravel(), cc.ravel(), np.asarray(p)).reshape(shape)
                for p in params]
    background = 3.0 + 0.12 * rr + 0.25 * cc
    data = sum(profile[:, :, None] * spectrum[None, None, :]
               for profile, spectrum in zip(profiles, spectra))
    data += background[:, :, None]
    cube = BlobCube(data, np.ones(shape, bool), np.array([50, 51, 52]), 0, 0,
                    0.02, 10.0, 0.0)
    hints = [type("Apex", (), {"row": p[0], "col": p[1]})() for p in params]
    fit = fit_overlap(cube, 2, include_background=True, background_order=1,
                      apex_hint=hints, n_starts=3, seed=4)

    truth = np.asarray([profile.sum() * spectrum.sum()
                        for profile, spectrum in zip(profiles, spectra)])
    recovered = sorted(fit.components, key=lambda component: component.row0)
    measured = np.asarray([component.volume for component in recovered])
    assert measured == pytest.approx(truth, rel=0.02)


def test_measured_template_fit_recovers_shifted_irregular_peaks():
    from pygcgc.core.overlap import _render_template

    tr, tc = np.indices((17, 5), dtype=float)
    template = np.exp(-0.5 * ((tr-7.0)/2.0) ** 2 - 0.5 * ((tc-2.0)/0.7) ** 2)
    template += 0.35 * np.exp(-0.5 * ((tr-11.0)/3.5) ** 2 -
                              0.5 * ((tc-2.5)/0.9) ** 2)
    shape = (45, 9)
    centres = ((16.2, 3.1), (28.4, 5.2))
    profiles = [_render_template(template, shape, *centre) for centre in centres]
    spectra = np.asarray([[7000.0, 200.0, 50.0], [100.0, 4200.0, 500.0]])
    data = sum(profile[:, :, None] * spectrum[None, None, :]
               for profile, spectrum in zip(profiles, spectra))
    cube = BlobCube(data, np.ones(shape, bool), np.array([50, 51, 52]), 0, 0,
                    0.02, 10.0, 0.0)
    hints = [type("Apex", (), {"row": row, "col": col})() for row, col in centres]
    fit = fit_template_overlap(cube, template, 2, apex_hint=hints, n_starts=3, seed=2)

    recovered = sorted(fit.components, key=lambda component: component.row0)
    assert [component.row0 for component in recovered] == pytest.approx(
        [centre[0] for centre in centres], abs=0.15)
    assert [component.col0 for component in recovered] == pytest.approx(
        [centre[1] for centre in centres], abs=0.15)
    assert [component.volume for component in recovered] == pytest.approx(
        spectra.sum(axis=1), rel=0.02)
    assert np.allclose(predict_overlap(fit, cube), data, rtol=0.02, atol=1e-3)


def test_component_maps_overlap_and_weak_amount_is_retained_with_warning(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    fit = fit_overlap(cube, 2)
    maps = component_maps(fit, cube)
    assert maps.shape == cube.mask.shape + (2,)
    assert np.any((maps[:, :, 0] > 0) & (maps[:, :, 1] > 0))
    assert np.allclose(maps.sum(axis=2), predict_overlap(fit, cube).sum(axis=2))

    fit.components[1].volume = 0.1 * fit.components[0].volume
    assessments = assess_components(fit)
    assert assessments[0].confidence == "provisional"
    assert assessments[1].confidence == "low"
    assert assessments[1].relative_to_largest == pytest.approx(0.1)


def test_fixed_shape_spectrum_refit_keeps_positions_and_recovers_signal(overlap_scene):
    from pygcgc.core.overlap import refit_spectra

    _run, _fold, _blobs, cube = overlap_scene
    screened = fit_overlap(cube, 2, include_background=True, background_order=1)
    refitted = refit_spectra(screened, cube, include_background=True, background_order=1)
    assert [(c.row0, c.col0) for c in refitted.components] == pytest.approx(
        [(c.row0, c.col0) for c in screened.components])
    assert all(component.spectrum.shape == cube.ions.shape
               for component in refitted.components)
    assert refitted.weighted_rss <= screened.weighted_rss * 1.001


def test_fast_shared_width_fit_recovers_two_known_overlapping_areas():
    from types import SimpleNamespace
    from pygcgc.core.overlap import fit_shapes_fast, refit_spectra

    rr, cc = np.indices((29, 7), dtype=float)
    centres = [(12.0, 2.2), (17.0, 2.8)]
    spectra = np.asarray([[5000., 1500., 200., 50.],
                          [100., 800., 3700., 400.]])
    profiles = [np.exp(-0.5 * (((rr-r)/3.0)**2 + ((cc-c)/0.8)**2))
                for r, c in centres]
    background = 20.0 + rr * 0.2 + cc * 0.4
    data = sum(profile[:, :, None] * spectrum[None, None, :]
               for profile, spectrum in zip(profiles, spectra))
    data += background[:, :, None]
    cube = BlobCube(data, np.ones(rr.shape, bool), np.array([50, 51, 52, 53]),
                    0, 0, 0.02, 10.0, 0.0)
    hints = [SimpleNamespace(row=r, col=c) for r, c in centres]
    shaped = fit_shapes_fast(cube, 2, apex_hint=hints, max_nfev=30)
    fit = refit_spectra(shaped, cube, background_order=1)
    truth = [profile.sum() * spectrum.sum()
             for profile, spectrum in zip(profiles, spectra)]
    assert shaped.n_shape_parameters == 6
    assert fit.components[0].sigma_row == pytest.approx(fit.components[1].sigma_row)
    assert [component.volume for component in fit.components] == pytest.approx(truth, rel=0.1)


def test_fitted_component_group_round_trips_and_reconstructs_contours(overlap_scene):
    from pygcgc.core.components import (ComponentGroup, ComponentSet, component_rows,
                                        group_from_fit, write_components_csv)

    _run, _fold, _blobs, cube = overlap_scene
    fit = fit_overlap(cube, 2, model="gaussian", include_background=True,
                      background_order=1, n_starts=3)
    group = group_from_fit("RUN", 4, fit, cube)
    assert np.isfinite(group.restart_position_spread)
    restored = ComponentGroup.from_dict(group.to_dict())
    assert restored.id == "RUN:region:4"
    assert len(restored.components) == 2
    assert restored.components[0].volume == pytest.approx(group.components[0].volume)
    assert np.allclose(restored.components[0].spectrum, group.components[0].spectrum)
    for index in range(2):
        x, y = restored.contour(index)
        assert len(x) == len(y) > 0
        assert np.isfinite(x).any() and np.isfinite(y).any()

    component_set = ComponentSet.from_dict(ComponentSet([group]).to_dict())
    assert component_set.group_for_region(4) is not None
    assert len(component_set.components) == 2

    rows = component_rows(component_set, "RUN")
    assert len(rows) == 2
    assert rows[0]["estimated_intensity"] > 0
    assert rows[0]["group_unresolved"] == group.unresolved

    fit.restart_position_spread = 2.5
    unstable = group_from_fit("RUN", 5, fit, cube, max_restart_position_spread=1.0)
    assert unstable.unresolved
    assert "repeated fits move" in unstable.warning


def test_component_csv_includes_quality_warnings(overlap_scene, tmp_path):
    import csv
    from pygcgc.core.components import ComponentSet, group_from_fit, write_components_csv

    _run, _fold, _blobs, cube = overlap_scene
    fit = fit_overlap(cube, 2)
    group = group_from_fit("RUN", 2, fit, cube, min_explained_fraction=2.0)
    path = write_components_csv(ComponentSet([group]), tmp_path / "components.csv", "RUN")
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    assert len(rows) == 2
    assert rows[0]["run"] == "RUN"
    assert rows[0]["group_unresolved"] == "True"
    assert "fit explains only" in rows[0]["group_warning"]


# ------------------------------------------------------------- 4. MCR-ALS
def test_unimodal_1d_caps_toward_the_peak():
    # Peak (argmax) is index 3 (value 7): everything left is capped by its
    # neighbour closer to the peak, everything right likewise.
    out = _unimodal_1d(np.array([0.0, 5.0, 3.0, 7.0, 2.0]))
    assert list(out) == [0.0, 3.0, 3.0, 7.0, 2.0]


def test_mcr_als_recovers_ion_dominant_spectra_for_two_components(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    result = mcr_als(cube, k=2, unimodal=True, seed=0)
    assert isinstance(result, MCRResult)
    assert result.spectra.shape == (2, cube.ions.size)
    dominant = {int(np.argmax(result.spectra[j])) for j in range(2)}
    ion50_idx = int(np.flatnonzero(cube.ions == 50)[0])
    ion51_idx = int(np.flatnonzero(cube.ions == 51)[0])
    assert dominant == {ion50_idx, ion51_idx}
    assert result.concentrations.shape == cube.mask.shape + (2,)


def test_mcr_restarts_reports_stable_two_component_solution(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    stability = mcr_restarts(cube, 2, n_starts=4, max_iter=100)
    assert stability.residuals.shape == (4,)
    assert stability.explained_variance.shape == (4,)
    assert stability.mean_spectral_similarity > .9
    assert stability.minimum_spectral_similarity > .8
    assert stability.best.spectra.shape[0] == 2


def test_mcr_als_rejects_more_components_than_footprint_pixels(overlap_scene):
    _run, _fold, _blobs, cube = overlap_scene
    with pytest.raises(ValueError):
        mcr_als(cube, k=cube.mask.sum() + 1)


# ------------------------------------------------------------- 5. PARAFAC
def _synthetic_trilinear_tensor(rng, I, J, K, r):
    A = rng.random((I, r)) + 0.1
    B = rng.random((J, r)) + 0.1
    C = rng.random((K, r)) + 0.1
    tensor = np.einsum('ir,jr,kr->ijk', A, B, C)
    return tensor, (A, B, C)


def test_parafac_reconstructs_a_clean_trilinear_tensor():
    rng = np.random.default_rng(0)
    tensor, _true = _synthetic_trilinear_tensor(rng, 6, 9, 5, 2)
    fit = parafac(tensor, r=2, non_negative=True, seed=1)
    assert isinstance(fit, PARAFACResult)
    recon = np.einsum('ir,jr,kr->ijk', *fit.factors)
    assert np.allclose(recon, tensor, atol=1e-3, rtol=1e-2)
    assert fit.resid < 1e-3


def test_core_consistency_is_high_at_the_true_rank_and_lower_when_overfactored():
    rng = np.random.default_rng(2)
    tensor, _true = _synthetic_trilinear_tensor(rng, 8, 10, 6, 2)
    fit_true = parafac(tensor, r=2, seed=3)
    cc_true = core_consistency(tensor, fit_true.factors)
    fit_over = parafac(tensor, r=4, seed=3)
    cc_over = core_consistency(tensor, fit_over.factors)
    assert cc_true > 80.0
    assert cc_over < cc_true


def test_split_half_stability_is_high_for_a_clean_trilinear_tensor():
    rng = np.random.default_rng(4)
    tensor, _true = _synthetic_trilinear_tensor(rng, 12, 10, 6, 2)
    score = split_half_stability(tensor, r=2, axis=0, seed=5)
    assert score > 0.8


def test_split_half_stability_needs_at_least_four_slices():
    tensor = np.random.default_rng(0).random((2, 5, 4))
    with pytest.raises(ValueError):
        split_half_stability(tensor, r=1, axis=0)


# ------------------------------------------------------------- build_run_tensor
@dataclass
class _Prepped:
    clean: object
    fold: object


def test_build_run_tensor_reads_each_run_at_the_matched_position():
    prepped = {}
    for i, amp in enumerate((1e5, 2e5, 0.5e5)):
        run = make_run([(100.0, 1.5, 50, amp)], period=PERIOD, phase=0.0, rate=RATE, duration=130.0,
                       mz=(48, 55), sigma1=4.0, sigma2=0.06, noise=0.0, seed=i)
        fold = Fold(run.times, FoldSpec(period=PERIOD, phase=0.0))
        prepped[f"R{i}"] = _Prepped(clean=run, fold=fold)

    tensor, run_names, ions = build_run_tensor(prepped, ref_t1=100.0, ref_t2=1.5, mz=[50],
                                               t1_window_s=8.0, t2_window_s=0.2, n_t1=5, n_t2=5)
    assert run_names == ["R0", "R1", "R2"]
    assert tensor.shape == (3, 25, 1)
    assert np.all(tensor >= 0)
    # The middle run (2x the amplitude) should read consistently higher than the faintest one.
    assert tensor[1].sum() > tensor[2].sum()
