"""🧩 Peak-clustering step in the GUI (pytest-qt, headless)."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
pytest.importorskip("pyqtgraph")

from tests.conftest import GUI_PHASE, make_run, write_run_andi  # noqa: E402
from tests.test_gui import TIMEOUT  # noqa: E402

ALKANE_IONS = (57, 71, 85, 99, 113)
AROMATIC_IONS = (65, 77, 91, 105, 119)
WEIGHTS = (1.0, 0.7, 0.4, 0.2, 0.1)
N_ALKANE = {"R1": 6, "R2": 6, "R3": 6}
N_AROMATIC = {"R1": 4, "R2": 6, "R3": 8}     # the aromatic share rises from R1 to R3


def _peaks(n_alkane, n_aromatic):
    peaks = []
    for j in range(n_alkane + n_aromatic):
        ions = ALKANE_IONS if j < n_alkane else AROMATIC_IONS
        # Off the sampling grid: an apex exactly between two pixels is a tie and the detector splits it.
        t1, t2 = 61.3 + 55.0 * j, 1.207 + 0.11 * (j % 12)
        peaks += [(t1, t2, ion, 2e5 * w) for ion, w in zip(ions, WEIGHTS)]
    return peaks


@pytest.fixture
def cluster_folder(tmp_path):
    folder = tmp_path / "cluster_pack"
    (folder / "data").mkdir(parents=True)
    for k, name in enumerate(("R1", "R2", "R3")):
        run = make_run(_peaks(N_ALKANE[name], N_AROMATIC[name]), phase=GUI_PHASE, rate=50.0, mz=(50, 130),
                       duration=900.0, seed=k + 1)
        write_run_andi(run, folder / "data" / f"{name}.cdf")
    (folder / "sequence.csv").write_text("file,type\nR1.cdf,sample\nR2.cdf,sample\nR3.cdf,sample\n")
    (folder / "method.json").write_text('{"modulation_period_s": 4.0}')
    return folder


@pytest.fixture
def cwin(qtbot, cluster_folder):
    from pygcgc.gui.mainwindow import MainWindow
    w = MainWindow()
    w.jobs_panel.clear_after_ms = 10
    qtbot.addWidget(w)
    w.resize(1400, 900)
    w.show()
    loaded = []
    w.open_folder(cluster_folder, on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and w.jobs.active_count() == 0 and w.view2d.result is not None,
                    timeout=TIMEOUT)
    qtbot.waitUntil(lambda: all(e.phase_method != "unset" for e in w.project.entries), timeout=TIMEOUT)
    yield w
    w.jobs.cancel_all()
    w.jobs.wait(5000)
    w.undo.clear()


def _go(qtbot, w):
    w.step_actions["cluster"].trigger()
    w.workflow.goRequested.emit("cluster")
    qtbot.waitUntil(lambda: w.clusters.result is not None and not w.clusters.running, timeout=TIMEOUT)


def _go_again(qtbot, w):
    """Cluster again after a result already exists (the first result would satisfy `_go`'s wait at once)."""
    before = w.clusters.result
    w.workflow.goRequested.emit("cluster")
    qtbot.waitUntil(lambda: w.clusters.result is not before and not w.clusters.running, timeout=TIMEOUT)


def _go_regions(qtbot, w):
    w.step_actions["cluster"].trigger()
    w.workflow.goRequested.emit("cluster")
    qtbot.waitUntil(lambda: w.clusters.regions is not None and not w.clusters.running, timeout=TIMEOUT)


def _regions_mode(w, n=2):
    w.cluster_panel.unit_box.setCurrentIndex(1)
    w.cluster_panel.n_regions.setValue(n)


def test_pixel_regions_divide_the_run_without_peak_detection(qtbot, cwin, tmp_path):
    w = cwin
    i = w.project.index_of("R2")
    w.set_run(i)
    _regions_mode(w)
    assert w.clusters.unit == "regions" and w.clusters.n_regions == 2
    assert not w.cluster_panel.merge_btn.isEnabled() and not w.cluster_panel.threshold.isEnabled()   # peak-only controls
    assert w.cluster_panel.table.horizontalHeaderItem(1).text() == "Pixels"
    _go_regions(qtbot, w)
    assert w.clusters.result is None and len(w.clusters.regions) == 2 and w.workflow.is_done("cluster")
    assert w.cluster_panel.table.rowCount() == 2 and "no peak detection" in w.cluster_panel.notes.text()
    assert not w.clusters.blobs_by_run                                        # no blobs were needed
    qtbot.waitUntil(lambda: w.view2d.result is not None and w.view2d.result.name == "R2" and w.view2d.cluster_layer.isVisible(),
                    timeout=TIMEOUT)
    painted = w.view2d._cluster_rgba[..., 3] > 0
    assert painted.any() and not painted.all()
    with pytest.raises(ValueError):
        w.clusters.merge([0, 1])                                              # regions cannot be merged
    assert not w.clusters.run_check("stability")
    path = tmp_path / "fig_R2_regions.png"
    w.save_cluster_figure(path, 89.0)
    assert path.stat().st_size > 5000
    assert "divided into regions of similar mass spectrum" in path.with_suffix(".txt").read_text(encoding="utf-8")


def test_switching_between_peaks_and_regions_keeps_each_result(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R2"))
    _go(qtbot, w)                                                             # peaks
    n_peak_clusters = len(w.clusters.result)
    _regions_mode(w, 3)
    assert w.cluster_panel.table.rowCount() == 0 and "Not clustered yet" in w.cluster_panel.status.text()
    _go_regions(qtbot, w)
    assert w.cluster_panel.table.rowCount() == 3
    w.cluster_panel.unit_box.setCurrentIndex(0)
    assert w.cluster_panel.table.rowCount() == n_peak_clusters                # the peak result is still there
    assert w.cluster_panel.table.horizontalHeaderItem(1).text() == "Peaks"


def test_pooled_regions_give_a_run_table_of_signal_shares(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R1"))
    w.cluster_panel.scope.setCurrentIndex(1)
    _regions_mode(w)
    _go_regions(qtbot, w)
    assert sorted(w.clusters.region_runs) == ["R1", "R2", "R3"]
    assert w.cluster_panel.run_table.isVisible() and w.cluster_panel.run_table.rowCount() == 3
    feats, matrix = w.clusters.fingerprint()
    assert matrix.shape == (3, 2) and np.allclose(matrix.sum(axis=1), 100.0)
    assert "Pool at least four runs" in (_expect_error(w.clusters.compare_groups))


def _expect_error(call):
    try:
        call()
    except ValueError as exc:
        return str(exc)
    return ""


def test_density_method_and_retention_weight_reach_the_clustering(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R2"))
    assert not w.cluster_panel.retention.isEnabled() and w.cluster_panel.threshold.isEnabled()
    w.cluster_panel.method_box.setCurrentIndex(1)
    w.cluster_panel.retention.setValue(2.0)
    assert w.cluster_panel.retention.isEnabled() and not w.cluster_panel.threshold.isEnabled()
    s = w.clusters.effective_settings()
    assert s.method == "density" and s.retention_weight == 2.0
    _go(qtbot, w)
    assert w.clusters.result.settings.method == "density"                    # it ran (12 peaks: too few for a 15-peak cluster)
    saved = w.project.to_dict()["clusters"]
    assert saved["settings"]["method"] == "density" and saved["unit"] == "peaks"


def test_a_cluster_can_be_named_and_the_name_reaches_the_figure_key_and_csv(qtbot, cwin, tmp_path, monkeypatch):
    from PySide6 import QtCore, QtWidgets
    w = cwin
    w.set_run(w.project.index_of("R2"))
    _go(qtbot, w)
    row = next(r for r in range(w.cluster_panel.table.rowCount())
               if w.cluster_panel.table.item(r, 0).data(0x0100) == 0)                  # cluster 0, wherever sorting put it
    w.cluster_panel.table.item(row, w.cluster_panel.NAME_COLUMN).setText("  alkylbenzenes ")
    assert w.clusters.name_map() == {0: "alkylbenzenes"}
    editable = QtCore.Qt.ItemFlag.ItemIsEditable
    assert w.cluster_panel.table.item(row, w.cluster_panel.NAME_COLUMN).flags() & editable
    assert not w.cluster_panel.table.item(row, 1).flags() & editable                    # every other column is read-only
    fig = tmp_path / "f.png"
    w.save_cluster_figure(fig, 89.0)
    assert "c0: alkylbenzenes" in fig.with_suffix(".txt").read_text(encoding="utf-8")
    csv_path = tmp_path / "p.csv"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(csv_path), "")))
    w.export_cluster_peaks()
    assert csv_path.read_text().splitlines()[0].endswith("cluster_name")
    _go_again(qtbot, w)
    assert w.clusters.name_map() == {}                                                  # clustering again renumbers, so names go


def test_pqn_normalisation_reaches_the_run_table_and_is_saved(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R1"))
    w.cluster_panel.scope.setCurrentIndex(1)
    _go(qtbot, w)
    plain = w.clusters.fingerprint()[1].copy()
    w.cluster_panel.norm_box.setCurrentIndex(1)
    assert w.clusters.normalise == "pqn"
    scaled = w.clusters.fingerprint()[1]
    assert scaled.shape == plain.shape and not np.allclose(scaled, plain)
    assert w.project.to_dict()["clusters"]["normalise"] == "pqn"


def test_the_unit_and_region_count_are_saved(qtbot, cwin):
    w = cwin
    _regions_mode(w, 5)
    saved = w.project.to_dict()["clusters"]
    assert saved["unit"] == "regions" and saved["n_regions"] == 5
    w.clusters.from_dict(saved)
    assert w.clusters.unit == "regions" and w.clusters.n_regions == 5 and w.clusters.regions is None


def test_cluster_is_step_ten_after_quantify(qtbot, cwin):
    assert cwin.workflow.keys()[-1] == "cluster"
    assert cwin.step_actions["cluster"].text().lstrip("✓▸ ") == "10 🧩 Cluster"
    assert cwin.step_actions["quant"].text().lstrip("✓▸ ") == "9 🧪 Quantify"     # nothing renumbered


def test_go_clusters_this_run_and_outlines_the_peaks(qtbot, cwin):
    w = cwin
    i = w.project.index_of("R2")
    w.set_run(i)
    _go(qtbot, w)
    assert len(w.clusters.result) == 2 and w.workflow.is_done("cluster")
    assert w.cluster_panel.table.rowCount() == 2
    assert "legacy detector run here" in w.cluster_panel.notes.text()      # no 🔍 blobs: found privately
    assert not w.cluster_panel.run_table.isVisible()                        # one run: no run table
    assert w._cluster_overlays() == []                                      # outlines are off by default
    w.cluster_panel.overlay_box.setChecked(True)
    overlays = w._cluster_overlays()
    assert overlays and all(len(o["outline"][0]) > 0 for o in overlays)
    w.cluster_panel.select_cluster(1)
    assert list(w._cluster_colours()) == [1]                                # selecting colours just that cluster
    w.step_actions["quant"].trigger()
    assert w._cluster_overlays() != []                                      # a layer, not gated to the step
    w.cluster_panel.overlay_box.setChecked(False)
    assert w._cluster_overlays() == []                                      # off once the checkbox is off


def test_the_2d_view_is_filled_with_cluster_colours(qtbot, cwin):
    w = cwin
    i = w.project.index_of("R2")
    w.set_run(i)
    _go(qtbot, w)
    qtbot.waitUntil(lambda: w.view2d.result is not None and w.view2d.result.name == "R2" and w.view2d.cluster_layer.isVisible(),
                    timeout=TIMEOUT)
    rgba = w.view2d._cluster_rgba
    assert rgba.shape == w.view2d.result.image.shape + (4,) and rgba.dtype == np.uint8
    painted = rgba[..., 3] > 0
    assert painted.any() and not painted.all()
    from pygcgc.core.cluster import CLUSTER_COLOURS
    assert {tuple(c) for c in rgba[painted][:, :3]} <= set(CLUSTER_COLOURS)   # both clusters are in the palette
    assert len({tuple(c) for c in rgba[painted][:, :3]}) == 2                  # and they differ
    w.cluster_panel.fill_box.setChecked(False)
    assert not w.view2d.cluster_layer.isVisible()
    w.cluster_panel.fill_box.setChecked(True)
    w.step_actions["quant"].trigger()
    assert w.view2d.cluster_layer.isVisible()                                   # a layer, not gated to the step
    w.cluster_panel.fill_box.setChecked(False)
    assert not w.view2d.cluster_layer.isVisible()                               # off once the checkbox is off


def test_the_coloured_chromatogram_is_saved_with_a_caption(qtbot, cwin, tmp_path):
    w = cwin
    w.set_run(w.project.index_of("R3"))
    _go(qtbot, w)
    path = tmp_path / "fig_R3_clusters.png"
    w.save_cluster_figure(path, 89.0)
    assert path.stat().st_size > 5000
    caption = path.with_suffix(".txt").read_text(encoding="utf-8")
    assert "coloured by the spectral cluster" in caption and "c0:" in caption and "cosine distance of 0.3" in caption


def test_merge_selected_joins_clusters(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R2"))
    _go(qtbot, w)
    assert len(w.clusters.result) == 2
    w.merge_selected_clusters()                                             # nothing selected: refused, nothing changes
    assert len(w.clusters.result) == 2
    w.cluster_panel.table.selectAll()
    w.merge_selected_clusters()
    assert len(w.clusters.result) == 1 and w.cluster_panel.table.rowCount() == 1


def test_stability_and_distance_scan_report_in_the_panel(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R2"))
    _go(qtbot, w)
    w.clusters.run_check("scan")
    qtbot.waitUntil(lambda: bool(w.clusters.check_text), timeout=TIMEOUT)
    assert "Merge distance" in w.cluster_panel.check_label.text() and "0.3 → 2" in w.cluster_panel.check_label.text()
    w.clusters.check_text = ""
    w.clusters.run_check("stability")
    qtbot.waitUntil(lambda: bool(w.clusters.check_text), timeout=TIMEOUT)
    assert "adjusted Rand index 1.00" in w.cluster_panel.check_label.text()


def test_group_comparison_needs_enough_runs_and_reports_permanova(qtbot, cwin, tmp_path, monkeypatch):
    from PySide6 import QtWidgets
    from pygcgc.core.cluster import PeakTable, cluster_peaks
    from tests.test_cluster import _two_groups
    w = cwin
    w.set_run(w.project.index_of("R1"))
    w.cluster_panel.scope.setCurrentIndex(1)
    _go(qtbot, w)
    w.test_cluster_groups()                                                 # three runs, no groups
    assert "at least four runs" in w.cluster_panel.group_label.text()
    assert w.cluster_panel.group_box.isVisible()
    # Swap in a pooled result with two clearly different groups of six runs, to exercise the statistics.
    tables, names, groups = _two_groups(6)
    pooled = PeakTable.concat(tables)
    cs = w.clusters
    cs.table, cs.result, cs.run_names = pooled, cluster_peaks(pooled), names
    cs.summary = cs.result.summary(pooled)
    cs.groups = {}
    w.test_cluster_groups()
    assert "at least two groups" in w.cluster_panel.group_label.text()
    for n, g in zip(names, groups):
        w.clusters.set_group(n, g)
    w.test_cluster_groups()
    text = w.cluster_panel.group_label.text()
    assert "PERMANOVA" in text and "p = 0.00" in text and "low (n = 6) vs high" not in text     # groups sort alphabetically
    assert "high (n = 6) vs low (n = 6)" in text and "Clusters differing at q < 0.05" in text
    assert w.project.to_dict()["clusters"]["groups"] == dict(zip(names, groups))               # and they are saved
    png = tmp_path / "pca.png"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(png), "")))
    w.export_cluster_pca()
    assert png.stat().st_size > 5000 and "Principal component analysis of 12" in png.with_suffix(".txt").read_text(encoding="utf-8")


def test_editing_the_group_column_reaches_the_state(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R1"))
    w.cluster_panel.scope.setCurrentIndex(1)
    _go(qtbot, w)
    item = w.cluster_panel.run_table.item(0, 1)
    item.setText("control")
    assert w.clusters.groups == {w.cluster_panel.run_table.item(0, 0).text(): "control"}
    item.setText("")
    assert w.clusters.groups == {}


def test_the_14_eV_partner_is_added_to_the_distance_when_asked(qtbot, cwin):
    w = cwin
    a, b = w.project.index_of("R1"), w.project.index_of("R2")
    w.project.entries[a].import_meta = {"acquisition_id": "x", "energy_ev": 70.0, "processing_state": "instrument"}
    w.project.entries[b].import_meta = {"acquisition_id": "x", "energy_ev": 14.0, "processing_state": "instrument"}
    w.set_run(a)
    w.cluster_panel.soft_box.setChecked(True)
    assert w.clusters.use_soft and w.clusters.effective_settings().soft_weight == 0.5
    _go(qtbot, w)
    assert w.clusters.result.soft_used is not None and w.clusters.table.soft_spectra is not None
    assert "14 eV spectra used (1 of 1 runs had a partner)" in w.cluster_panel.notes.text()
    w.cluster_panel.soft_box.setChecked(False)
    _go_again(qtbot, w)
    assert w.clusters.result.soft_used is None


def test_made_of_shows_the_class_of_each_cluster_after_classifying(qtbot, cwin):
    w = cwin
    i = w.project.index_of("R2")
    w.set_run(i)
    _go(qtbot, w)
    assert all(w.cluster_panel.table.item(r, 8).text() == "" for r in range(w.cluster_panel.table.rowCount()))
    w.start_classify()
    qtbot.waitUntil(lambda: w.classes.result_for(i) is not None and not w.classes.running, timeout=TIMEOUT)
    made = [w.cluster_panel.table.item(r, 8).text() for r in range(w.cluster_panel.table.rowCount())]
    assert all("(" in m and m.endswith("%)") for m in made)


def test_pooling_runs_gives_a_run_table_that_tracks_the_aromatic_share(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R1"))
    w.cluster_panel.scope.setCurrentIndex(1)                                # every 70 eV run, pooled
    assert w.clusters.scope == "all"
    _go(qtbot, w)
    assert sorted(w.clusters.run_names) == ["R1", "R2", "R3"]
    assert w.cluster_panel.run_table.isVisible() and w.cluster_panel.run_table.rowCount() == 3
    feats, matrix = w.clusters.fingerprint()
    names = w.clusters.run_names
    assert np.allclose(matrix.sum(axis=1), 100.0)                           # every peak was assigned
    aromatic = next(f.id for f in feats if 65 in w.clusters.result.top_ions(f.id, 3))
    share = dict(zip(names, matrix[:, aromatic]))
    assert share["R1"] < share["R2"] < share["R3"]


def test_editing_a_setting_marks_the_result_out_of_date(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R1"))
    _go(qtbot, w)
    assert "Out of date" not in w.cluster_panel.status.text()
    w.cluster_panel.threshold.setValue(0.5)
    assert w.clusters.stale and "Out of date" in w.cluster_panel.status.text()


def test_the_14_eV_run_is_refused(qtbot, cwin):
    w = cwin
    i = w.project.index_of("R1")
    w.project.entries[i].import_meta = {"acquisition_id": "x", "energy_ev": 14.0, "processing_state": "instrument"}
    w.set_run(i)
    assert "70 eV" in w.clusters.why_not(i)
    assert not w.workflow.step_info("cluster")["go_enabled"]
    assert w.clusters.eligible() == [w.project.index_of("R2"), w.project.index_of("R3")]


def test_settings_survive_a_project_round_trip_and_results_do_not(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R1"))
    w.cluster_panel.fixed_box.setChecked(True)
    w.cluster_panel.fixed_n.setValue(7)
    w.cluster_panel.siloxane_box.setChecked(True)
    _go(qtbot, w)
    saved = w.project.to_dict()["clusters"]
    assert saved["settings"]["n_clusters"] == 7 and saved["skip_siloxane"] is True
    w.clusters.from_dict(saved)
    assert w.clusters.result is None and w.clusters.settings.n_clusters == 7 and w.clusters.skip_siloxane


def test_removing_a_clustered_run_forgets_the_result(qtbot, cwin):
    w = cwin
    w.set_run(w.project.index_of("R1"))
    _go(qtbot, w)
    w.remove_runs([w.project.index_of("R1")])
    assert w.clusters.result is None and w.cluster_panel.table.rowCount() == 0


def test_exports_write_the_peak_and_run_tables(qtbot, cwin, tmp_path, monkeypatch):
    from PySide6 import QtWidgets
    w = cwin
    w.set_run(w.project.index_of("R1"))
    w.cluster_panel.scope.setCurrentIndex(1)
    _go(qtbot, w)
    peaks, runs = tmp_path / "p.csv", tmp_path / "r.csv"
    for path, call in ((peaks, w.export_cluster_peaks), (runs, w.export_cluster_runs)):
        monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", staticmethod(lambda *a, p=path, **k: (str(p), "")))
        call()
    assert len(peaks.read_text().splitlines()) == 1 + len(w.clusters.table)
    header = runs.read_text().splitlines()[0].split(",")
    assert header[:3] == ["feature", "t1_min", "t2_s"] and {"R1", "R2", "R3"} <= set(header)
