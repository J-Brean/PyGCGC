"""Model-order selection: how many peaks a window holds is decided by the data.

Each scenario is a small noisy cube with a known answer. Seeds (candidate apexes)
are deliberately wrong in number: the selected model order must not follow them.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from pygcgc.core.components import (ComponentDetectionSettings, ComponentGroup,
                                    _drop_duplicates, _fit_selected)
from pygcgc.core.modelselect import has_valley, select_components, spectral_cosine
from pygcgc.core.overlap import BlobCube, _emg_1d

IONS = np.arange(50, 62)
SETTINGS = ComponentDetectionSettings()


def _spectrum(peak_ions, rng=None):
    spectrum = np.full(IONS.size, 30.0)
    for ion, height in peak_ions.items():
        spectrum[ion - IONS[0]] = height
    return spectrum


A = _spectrum({50: 6000, 51: 2500, 52: 400, 53: 100, 54: 60, 60: 30})
B = _spectrum({54: 300, 55: 2800, 56: 4200, 57: 900, 58: 200, 59: 100})
A_LIKE = A * 1.0 + np.linspace(0, 40, IONS.size)     # same compound, different amount of noise floor


def _cube(peaks, shape=(45, 9), seed=0, noise=True):
    """peaks: (row, col, sigma_row, sigma_col, spectrum, tau_row)."""
    rr, cc = np.indices(shape, dtype=float)
    data = np.full(shape + (IONS.size,), 20.0)
    for row, col, sr, sc, spectrum, tau in peaks:
        gr = _emg_1d(rr, row, sr, tau) if tau else np.exp(-0.5 * ((rr - row) / sr) ** 2)
        profile = gr * np.exp(-0.5 * ((cc - col) / sc) ** 2)
        data = data + profile[:, :, None] * spectrum[None, None, :]
    if noise:
        rng = np.random.default_rng(seed)
        data = data + rng.normal(0.0, np.sqrt(np.clip(data, 1.0, None)))
    return BlobCube(data, np.ones(shape, bool), IONS.astype(float), 0, 0, 0.02, 10.0, 0.0)


def _seeds(*positions):
    return [SimpleNamespace(row=float(r), col=float(c)) for r, c in positions]


def _select(cube, seeds):
    screen = BlobCube(cube.data[:, :, :8], cube.mask, cube.ions[:8], 0, 0, 0.02, 10.0, 0.0)
    return select_components(screen, cube, seeds, SETTINGS)


def test_isolated_gaussian_is_one_component_whatever_the_seeds_say():
    cube = _cube([(22.0, 4.0, 3.0, 0.8, A, 0)])
    for seeds in (_seeds((22, 4)),
                  _seeds((22, 4), (16, 3), (29, 5)),                  # spurious apexes on the flanks
                  _seeds((17, 4), (27, 4), (22, 2), (22, 6))):
        assert _select(cube, seeds).fit.k == 1


def test_broad_isolated_peak_with_scattered_seeds_stays_one_component():
    cube = _cube([(22.0, 4.0, 6.0, 1.6, A, 0)], shape=(50, 11))
    selection = _select(cube, _seeds((15, 4), (29, 3), (22, 6), (20, 2)))
    assert selection.fit.k == 1
    assert selection.fit.components[0].sigma_row == pytest.approx(6.0, rel=0.15)


def test_tailing_peak_is_not_split_into_two_same_spectrum_components():
    cube = _cube([(18.0, 4.0, 2.5, 0.8, A, 6.0)], shape=(60, 9), seed=3)
    assert _select(cube, _seeds((20, 4), (30, 4))).fit.k == 1


def test_hidden_shoulder_is_found_without_any_second_seed():
    # One apparent TIC maximum: 3 rows apart at sigma 3, unequal spectra.
    cube = _cube([(20.0, 4.0, 3.0, 0.8, A, 0), (23.0, 4.5, 3.0, 0.8, B * 0.6, 0)])
    selection = _select(cube, _seeds((20, 4)))                       # a single seed
    assert selection.fit.k == 2
    rows = sorted(c.row0 for c in selection.fit.components)
    assert rows == pytest.approx([20.0, 23.0], abs=1.3)


def test_seed_count_never_caps_or_forces_the_model_order():
    two = _cube([(15.0, 3.0, 3.0, 0.8, A, 0), (28.0, 5.0, 3.0, 0.8, B, 0)])
    assert _select(two, _seeds((15, 3))).fit.k == 2                   # one seed, two peaks
    assert _select(two, _seeds((15, 3), (28, 5), (21, 4), (10, 2), (33, 6))).fit.k == 2  # five seeds, two peaks


def test_spectrally_distinct_overlapping_pair_is_separated_with_accurate_spectra():
    cube = _cube([(19.0, 4.0, 3.0, 0.8, A, 0), (25.0, 4.0, 3.0, 0.8, B, 0)], seed=5)
    fit = _select(cube, _seeds((19, 4), (25, 4))).fit
    assert fit.k == 2
    ordered = sorted(fit.components, key=lambda c: c.row0)
    assert spectral_cosine(ordered[0].spectrum, A) > 0.98
    assert spectral_cosine(ordered[1].spectrum, B) > 0.98


def test_same_spectrum_and_unresolved_is_reported_as_one_peak_with_a_note():
    # Two identical-spectrum peaks 4 rows apart at sigma 3 have no dip between them.
    cube = _cube([(20.0, 4.0, 3.0, 0.8, A, 0), (24.0, 4.0, 3.0, 0.8, A_LIKE, 0)], seed=2)
    selection = _select(cube, _seeds((20, 4), (24, 4)))
    assert selection.fit.k == 1


def test_same_spectrum_but_resolved_peaks_stay_two_components():
    # Isomer-like: identical spectra, 16 rows apart (well over 2 sigma) with a clear valley.
    cube = _cube([(14.0, 4.0, 3.0, 0.8, A, 0), (32.0, 4.0, 3.0, 0.8, A_LIKE, 0)], shape=(46, 9), seed=4)
    assert _select(cube, _seeds((14, 4), (32, 4))).fit.k == 2


def test_selection_is_deterministic():
    cube = _cube([(19.0, 4.0, 3.0, 0.8, A, 0), (25.0, 4.0, 3.0, 0.8, B, 0)], seed=9)
    first = _select(cube, _seeds((19, 4)))
    second = _select(cube, _seeds((19, 4)))
    assert [c.row0 for c in first.fit.components] == [c.row0 for c in second.fit.components]
    assert first.bic_gap == second.bic_gap


def test_valley_detection():
    def peak(row, col, volume=1e5, sigma_row=3.0, sigma_col=0.8):
        return SimpleNamespace(row0=row, col0=col, volume=volume, sigma_row=sigma_row,
                               sigma_col=sigma_col)
    assert has_valley(peak(10, 4), peak(24, 4))                        # 4.7 widths, clear dip
    assert not has_valley(peak(10, 4), peak(14, 4))                    # merged shoulder
    assert not has_valley(peak(10, 4), peak(20, 4, volume=1e4))        # tail-like: small neighbour, no dip


def test_fit_selected_flags_a_shape_the_data_want_split_but_spectra_cannot():
    # A strongly tailing peak: the data prefer more Gaussians, the spectra give no
    # reason to call them different compounds. One peak is reported, with a warning.
    cube = _cube([(18.0, 4.0, 2.5, 0.8, A, 8.0)], shape=(70, 9), seed=3)
    screen = BlobCube(cube.data[:, :, :8], cube.mask, cube.ions[:8], 0, 0, 0.02, 10.0, 0.0)
    group = _fit_selected(screen, cube, _seeds((20, 4), (30, 4)), 0, "T", SETTINGS)
    assert len(group.components) == 1
    assert group.unresolved
    assert group.warning


def _fake_group(index, row0, col0, components, n_rows=41, n_cols=9):
    group = ComponentGroup(
        id=f"T:region:{index}", region_index=index, components=components, model="gaussian",
        row0=row0, col0=col0, n_rows=n_rows, n_cols=n_cols, pitch=0.02, period=10.0,
        t1_start=0.0, bic=0.0, weighted_residual=0.0, explained_fraction=1.0,
        converged=True, unresolved=False)
    return group


def _fake_component(row, col, spectrum, volume=1e6, mz=None):
    from pygcgc.core.components import Component
    return Component(id="c", group_id="g", t1=0.0, t2=0.0, volume=volume,
                     mz=IONS.astype(float) if mz is None else np.asarray(mz, dtype=float),
                     spectrum=spectrum, model="gaussian", row0=row, col0=col,
                     sigma_row=3.0, sigma_col=0.8)


def test_duplicate_of_one_peak_from_two_windows_is_kept_once_in_the_better_centred_window():
    centred = _fake_group(0, 100, 10, [_fake_component(20.0, 4.0, A)])
    edge = _fake_group(1, 118, 10, [_fake_component(2.0, 4.0, A * 1.02)])   # same peak, on its window's edge
    entries = _drop_duplicates([(None, centred), (None, edge)], (500, 60), SETTINGS)
    assert len(entries) == 1
    assert entries[0][1] is centred


def test_duplicate_detection_compares_spectra_on_a_common_mz_axis():
    # Each window picks its own ions. The same compound seen through two different
    # ion lists must still be recognised as one peak, and different compounds must not.
    ions_a = np.arange(50, 62)
    ions_b = np.arange(53, 65)                                     # shifted by three ions
    spectrum = np.zeros(80)
    spectrum[[50, 51, 55, 56, 58]] = [900, 400, 3000, 2000, 500]  # keyed by m/z
    other = np.zeros(80)
    other[[52, 57, 59, 60, 63]] = [800, 500, 2500, 1500, 700]
    centred = _fake_group(0, 100, 10, [_fake_component(20.0, 4.0, spectrum[ions_a], mz=ions_a)])
    edge = _fake_group(1, 118, 10, [_fake_component(2.0, 4.0, spectrum[ions_b], mz=ions_b)])
    assert len(_drop_duplicates([(None, centred), (None, edge)], (500, 60), SETTINGS)) == 1
    # a different compound in its own window, well inside it and a few widths away
    stranger = _fake_group(1, 105, 10, [_fake_component(20.0, 4.0, other[ions_b], mz=ions_b)])
    assert len(_drop_duplicates([(None, centred), (None, stranger)], (500, 60), SETTINGS)) == 2


def test_distinct_neighbour_in_another_window_is_not_treated_as_duplicate():
    first = _fake_group(0, 100, 10, [_fake_component(20.0, 4.0, A)])
    second = _fake_group(1, 100, 10, [_fake_component(21.0, 4.5, B)])
    entries = _drop_duplicates([(None, first), (None, second)], (500, 60), SETTINGS)
    assert len(entries) == 2


def test_new_settings_round_trip_and_legacy_keys_are_ignored():
    settings = ComponentDetectionSettings(distinct_cosine=0.93, model_selection="seeds")
    restored = ComponentDetectionSettings.from_dict(settings.to_dict())
    assert restored.distinct_cosine == 0.93 and restored.model_selection == "seeds"
    legacy = settings.to_dict() | {"search_hidden_components": True}
    assert ComponentDetectionSettings.from_dict(legacy).distinct_cosine == 0.93
