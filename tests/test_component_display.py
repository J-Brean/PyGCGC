"""How fitted blobs are drawn and picked: smooth outlines, one blob per click, the residual map."""

import numpy as np
import pytest
from PySide6 import QtCore

from pygcgc.core.components import Component, ComponentGroup, ComponentSet
from pygcgc.core.overlap import _shape
from tests.test_gui import window  # noqa: F401  (the synthetic-folder main window fixture)

PITCH, PERIOD, T1_START = 0.02, 8.0, 100.0


def _component(row0, col0, sigma_row=3.0, sigma_col=0.8, model="gaussian", tau_row=None, name="c0"):
    return Component(id=name, group_id="g", t1=T1_START + (10 + col0 + 0.5) * PERIOD,
                     t2=(200 + row0) * PITCH, volume=1e6, mz=np.array([50.0, 51.0]),
                     spectrum=np.array([900.0, 100.0]), model=model, row0=row0, col0=col0,
                     sigma_row=sigma_row, sigma_col=sigma_col, tau_row=tau_row)


def _group(components, model="gaussian"):
    return ComponentGroup(id="g", region_index=0, components=components, model=model, row0=200, col0=10,
                          n_rows=41, n_cols=9, pitch=PITCH, period=PERIOD, t1_start=T1_START, bic=0.0,
                          weighted_residual=0.0, explained_fraction=1.0, converged=True, unresolved=False)


def _to_pixels(group, x, y):
    """Outline points back into the group's window pixel units."""
    x0, x1, y0, y1 = x[0::3], x[1::3], y[0::3], y[1::3]
    cols = (np.r_[x0, x1] - T1_START) / PERIOD - group.col0 - 0.5
    rows = np.r_[y0, y1] / PITCH - group.row0
    return rows, cols


@pytest.mark.parametrize("fraction", [0.10, 0.50])
def test_gaussian_outline_is_a_smooth_ellipse_at_the_requested_height(fraction):
    group = _group([_component(20.3, 4.4)])
    x, y = group.contour(0, fraction)
    assert x.size // 3 >= 90                                     # many points, not a handful per column
    rows, cols = _to_pixels(group, x, y)
    c = group.components[0]
    height = np.exp(-0.5 * (((rows - c.row0) / c.sigma_row) ** 2 + ((cols - c.col0) / c.sigma_col) ** 2))
    assert height == pytest.approx(fraction, abs=1e-9)


def test_outline_centre_matches_the_reported_position():
    # The outline and the table must agree on where the blob is: both mid-modulation.
    group = _group([_component(20.3, 4.4)])
    x, y = group.contour(0, 0.5)
    c = group.components[0]
    assert np.nanmean(x) == pytest.approx(c.t1, abs=1e-6)
    assert np.nanmean(y) == pytest.approx(c.t2, abs=1e-6)


def test_tailing_outline_follows_the_fitted_shape():
    component = _component(15.0, 4.0, sigma_row=2.5, model="emg_t2", tau_row=6.0)
    group = _group([component], model="emg_t2")
    x, y = group.contour(0, 0.10)
    rows, cols = _to_pixels(group, x, y)
    params = np.array([component.row0, component.col0, component.sigma_row, component.sigma_col,
                       component.tau_row])
    fine_r, fine_c = np.meshgrid(np.arange(-4, 45, 0.25), np.arange(-2, 11, 0.05), indexing="ij")
    peak = _shape("emg_t2", fine_r.ravel(), fine_c.ravel(), params).max()
    height = _shape("emg_t2", rows, cols, params) / peak
    assert height == pytest.approx(0.10, abs=0.01)
    assert rows.max() - component.row0 > component.row0 - rows.min()   # the tail side reaches further


def test_a_blob_with_no_signal_has_no_outline():
    component = _component(20.0, 4.0)
    component.spectrum = np.zeros(2)
    x, y = _group([component]).contour(0, 0.1)
    assert x.size == 0 and y.size == 0


# ---------------------------------------------------------------- GUI
def _two_blob_window(qtbot, window):
    """A run whose only fitted group holds two overlapping blobs."""
    from pygcgc.gui.blob_state import MS_AWARE

    w = window
    w.blobs.set_method(MS_AWARE)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    fold = w.project.fold(i)
    t1_start, period, pitch = float(fold.t1[0]), fold.spec.period, fold.row_pitch
    group = ComponentGroup(id="RUN_A:region:0", region_index=0, components=[], model="gaussian",
                           row0=40, col0=5, n_rows=41, n_cols=9, pitch=pitch, period=period,
                           t1_start=t1_start, bic=0.0, weighted_residual=0.0, explained_fraction=1.0,
                           converged=True, unresolved=False)
    for n, (row0, col0) in enumerate(((18.0, 4.0), (24.0, 4.6))):
        group.components.append(Component(
            id=f"RUN_A:region:0:component:{n}", group_id=group.id,
            t1=t1_start + (group.col0 + col0 + 0.5) * period, t2=(group.row0 + row0) * pitch,
            volume=1e6, mz=np.array([50.0, 51.0]), spectrum=np.array([900.0, 100.0]), model="gaussian",
            row0=row0, col0=col0, sigma_row=3.0, sigma_col=0.8))
    w.blobs.component_sets_by_name["RUN_A"] = ComponentSet(groups=[group])
    w._blobs_changed()
    w.act_blob_overlay.setChecked(True)
    return w, i, group


def _click(w, t1_s, t2_s):
    w.view2d.last_click = (t1_s / 60.0, t2_s)
    col, row = w.view2d.pixel_at(t1_s / 60.0, t2_s)[1], w.view2d.pixel_at(t1_s / 60.0, t2_s)[0]
    w._point_clicked(col, row * w.project.fold(w.index).row_pitch)


def test_clicking_one_blob_selects_and_highlights_only_that_blob(qtbot, window):
    w, _i, group = _two_blob_window(qtbot, window)
    target = group.components[1]
    _click(w, target.t1, target.t2)
    assert w.blobs_panel.selected_component_ids() == {target.id}
    wide = [o for o in w.view2d._component_overlays if o["width"] > 1]
    assert len(wide) == 1                                        # its outline alone, nothing else
    w.act_half_contours.setChecked(True)
    wide = [o for o in w.view2d._component_overlays if o["width"] > 1]
    assert len(wide) == 2                                        # plus the optional half-height ring


def test_each_blob_is_drawn_as_one_outline_unless_half_height_rings_are_asked_for(qtbot, window):
    w, _i, group = _two_blob_window(qtbot, window)
    blobs = [o for o in w.view2d._component_overlays if o["colour"] in
             {(0, 220, 255), (255, 105, 180)}]
    assert len(blobs) == len(group.components)
    assert all(o["style"] == QtCore.Qt.PenStyle.SolidLine for o in blobs)
    w.act_half_contours.setChecked(True)
    blobs = [o for o in w.view2d._component_overlays if o["colour"] in
             {(0, 220, 255), (255, 105, 180)}]
    assert len(blobs) == 2 * len(group.components)


def test_clicking_between_blobs_selects_the_one_whose_outline_contains_the_click(qtbot, window):
    w, _i, group = _two_blob_window(qtbot, window)
    first, second = group.components
    pitch = w.project.fold(w.index).row_pitch
    # 1.2 widths out on the side away from the second blob, which is 3.2 widths away
    _click(w, first.t1, first.t2 - 1.2 * first.sigma_row * pitch)
    assert w.blobs_panel.selected_component_ids() == {first.id}
    _click(w, second.t1, second.t2 + 1.2 * second.sigma_row * pitch)
    assert w.blobs_panel.selected_component_ids() == {second.id}


def test_clicking_empty_chromatogram_selects_no_blob(qtbot, window):
    w, _i, group = _two_blob_window(qtbot, window)
    far = group.components[0]
    _click(w, far.t1 + 40 * w.project.fold(w.index).spec.period, far.t2)
    assert w.blobs_panel.selected_component_ids() == set()


def test_fit_and_residual_map_opens_with_four_matching_images(qtbot, window):
    from pygcgc.gui.blob_state import MS_AWARE

    w = window
    w.blobs.set_method(MS_AWARE)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.components_for(i) is not None
                    and bool(w.blobs.components_for(i).groups), timeout=20000)
    w.run_signal_accounting()
    qtbot.waitUntil(lambda: getattr(w, "fit_map", None) is not None and w.fit_map.isVisible(), timeout=20000)
    shapes = {image.image.shape for image in w.fit_map.images}
    assert len(shapes) == 1
    d = w.fit_map._data
    np.testing.assert_allclose(d["residual"], d["tic"] - d["model"])
    assert len(w.fit_map.images) == 4
    share = d["fraction"][np.isfinite(d["fraction"])]
    assert share.size and share.min() >= 0.0
    measured = d["tic"] > 3.0 * d["noise"]
    np.testing.assert_array_equal(np.isfinite(d["fraction"]), measured)
    assert "explain" in w.fit_map.summary.text()
    w.fit_map.close()
