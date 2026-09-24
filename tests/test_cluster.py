"""Unsupervised peak clustering: within a run, pooled across runs, and the
run-level statistics built on it."""

import numpy as np
import pytest
from scipy.stats import f_oneway

from pygcgc import PeakSettings, feature_table_matrix, fold, find_blobs, pca, volcano
from pygcgc.core.cluster import (ClusterSettings, PeakTable, cluster_features, cluster_peaks, permanova,
                                 write_clusters_csv)
from tests.conftest import make_run

MZ = np.arange(40, 131)
ALKANE_IONS = (57, 71, 85, 99, 113)
AROMATIC_IONS = (65, 77, 91, 105, 119)


def _spectra(kind: np.ndarray, rng, jitter=0.15):
    """Two spectral families with a little run-to-run variation in the ratios."""
    ions = {0: ALKANE_IONS, 1: AROMATIC_IONS}
    spectra = np.zeros((kind.size, MZ.size), dtype=np.float32)
    for i, k in enumerate(kind):
        for ion, w in zip(ions[k], (1.0, 0.7, 0.4, 0.2, 0.1)):
            spectra[i, ion - MZ[0]] = 1e4 * w * (1 + jitter * rng.standard_normal())
    return spectra


def _table(kind, run, seed=0, volume=None):
    rng = np.random.default_rng(seed)
    n = kind.size
    return PeakTable(run=np.full(n, run, dtype=object), blob=np.arange(n), t1=rng.uniform(0, 3000, n),
                     t2=rng.uniform(0.5, 3, n), volume=rng.uniform(1e5, 1e6, n) if volume is None else volume,
                     co_elution=np.zeros(n, dtype=bool), spectra=_spectra(kind, rng), mz=MZ)


def _same_partition(labels, truth):
    """Every truth group maps to one label and no label is shared between groups."""
    mapped = {t: set(labels[truth == t]) for t in np.unique(truth)}
    return all(len(v) == 1 for v in mapped.values()) and len({next(iter(v)) for v in mapped.values()}) == len(mapped)


def test_two_spectral_families_are_recovered_from_real_blobs():
    # Six peaks of each family, well separated in retention so each is its own blob.
    peaks, truth = [], []
    for i in range(6):
        for fam, ions in enumerate((ALKANE_IONS, AROMATIC_IONS)):
            # Off the sampling grid: an apex exactly between two pixels is a tie and the detector splits it.
            t1, t2 = 61.3 + 130.0 * i + 55.0 * fam, 1.207 + 0.25 * i + 0.4 * fam
            peaks += [(t1, t2, ion, 2e5 * w * (1 + 0.1 * i)) for ion, w in zip(ions, (1.0, 0.7, 0.4, 0.2, 0.1))]
            truth.append((t1, fam))
    run = make_run(peaks, phase=0.7, rate=50.0, mz=(50, 130), duration=900.0, seed=1)
    blobs = find_blobs(run, fold(run, 4.0, 0.7), PeakSettings(snr_min=8), floor=1.0)
    assert len(blobs) == 12

    table = PeakTable.from_blobs(blobs, "R1")
    result = cluster_peaks(table)
    kind = np.array([truth[int(np.argmin([abs(t - t1) for t, _ in truth]))][1] for t1 in table.t1])
    assert len(result) == 2 and (result.labels >= 0).all()
    assert _same_partition(result.labels, kind)
    top = {frozenset(result.top_ions(k, 3)) for k in range(2)}
    assert top == {frozenset({57, 71, 85}), frozenset({65, 77, 91})}


def test_tree_on_a_subset_then_assignment_gives_the_same_answer():
    rng = np.random.default_rng(3)
    kind = rng.integers(0, 2, 200)
    table = _table(kind, "R1", seed=3)
    full = cluster_peaks(table)
    small = cluster_peaks(table, ClusterSettings(max_fit=40))
    assert len(small.fit_index) == 40
    assert _same_partition(full.labels, kind) and _same_partition(small.labels, kind)


def test_clusters_are_numbered_by_volume():
    kind = np.array([0] * 30 + [1] * 10)
    table = _table(kind, "R1", volume=np.where(kind == 1, 5e6, 1e5).astype(float))
    result = cluster_peaks(table)
    assert result.top_ions(0, 1) == [65] and result.top_ions(1, 1) == [57]     # each family's strongest ion
    assert [r.label for r in result.summary(table)] == [0, 1]
    assert result.summary(table)[0].percent > result.summary(table)[1].percent


def test_odd_spectra_are_left_unassigned_not_forced_in():
    kind = np.array([0] * 20 + [1] * 20)
    table = _table(kind, "R1")
    table.spectra[0] = 0
    table.spectra[0, 60 - MZ[0]] = 1e4                 # one ion: too few to describe a shape
    table.spectra[1] = 0
    table.spectra[1, [45 - MZ[0], 50 - MZ[0], 55 - MZ[0], 60 - MZ[0], 62 - MZ[0]]] = 1e4   # matches neither family
    result = cluster_peaks(table)
    assert result.labels[0] == -1 and result.labels[1] == -1
    assert (result.labels[2:] >= 0).all()


def test_a_tiny_group_is_dropped_by_min_cluster_size():
    kind = np.array([0] * 30 + [1] * 2)
    result = cluster_peaks(_table(kind, "R1"))
    assert len(result) == 1 and (result.labels[kind == 1] == -1).all()


def test_a_new_run_is_assigned_to_an_existing_segmentation():
    fit = _table(np.random.default_rng(1).integers(0, 2, 100), "R1", seed=1)
    result = cluster_peaks(fit)
    kind = np.random.default_rng(2).integers(0, 2, 60)
    new = _table(kind, "R2", seed=2)
    assert _same_partition(result.assign(new), kind)
    with pytest.raises(ValueError):
        result.assign(PeakTable(**{**new.__dict__, "spectra": new.spectra[:, :30], "mz": MZ[:30]}))


def test_pooling_keeps_the_shared_mz_range():
    a, b = _table(np.zeros(5, int), "A"), _table(np.zeros(5, int), "B")
    b = PeakTable(**{**b.__dict__, "spectra": b.spectra[:, 10:], "mz": MZ[10:]})
    pooled = PeakTable.concat([a, b])
    assert pooled.mz[0] == MZ[10] and pooled.spectra.shape == (10, MZ.size - 10)
    assert list(pooled.run) == ["A"] * 5 + ["B"] * 5


def test_settings_round_trip_through_json_types():
    s = ClusterSettings(exclude_ions=(73, 207), n_clusters=4)
    assert ClusterSettings.from_dict(__import__("json").loads(__import__("json").dumps(s.to_dict()))) == s


# ----------------------------------------------------------- across runs
def _two_groups(n_per_group=5, shift=0.5):
    """Runs of group 'low' hold mostly family 0; runs of group 'high' hold more of family 1."""
    tables, names, groups = [], [], []
    for g, share1 in (("low", 0.5 - shift / 2), ("high", 0.5 + shift / 2)):
        for i in range(n_per_group):
            rng = np.random.default_rng(100 * (g == "high") + i)
            kind = (rng.random(80) < share1).astype(int)
            name = f"{g}{i}"
            tables.append(_table(kind, name, seed=int(rng.integers(1e6))))
            names.append(name)
            groups.append(g)
    return tables, names, groups


def test_pooled_clusters_feed_the_existing_pca_and_volcano():
    tables, names, groups = _two_groups()
    pooled = PeakTable.concat(tables)
    result = cluster_peaks(pooled)
    assert len(result) == 2
    feats = cluster_features(pooled, result, relative=True)
    X = feature_table_matrix(feats, names)
    assert X.shape == (10, 2)
    assert np.allclose(X.sum(axis=1), 100.0)                       # nothing was left unassigned
    assert pca(feats, names).scores.shape == (10, 2)
    pts = volcano(feats, names[:5], names[5:])
    assert len(pts) == 2 and all(p.q_value < 0.05 for p in pts)     # the two families move in opposite directions


def test_a_run_with_no_peaks_still_appears_as_zeros():
    tables, names, _ = _two_groups(2)
    pooled = PeakTable.concat(tables)
    feats = cluster_features(pooled, cluster_peaks(pooled), run_names=names + ["empty"])
    assert all(f.volumes["empty"] == 0.0 for f in feats)


def test_permanova_matches_anova_in_one_dimension():
    rng = np.random.default_rng(0)
    x = np.concatenate([rng.normal(0, 1, 6), rng.normal(1.5, 1, 7), rng.normal(0.5, 1, 5)])
    groups = np.repeat(["a", "b", "c"], [6, 7, 5])
    f, p = permanova(x[:, None], groups, metric="euclidean", permutations=999)
    assert f == pytest.approx(f_oneway(x[groups == "a"], x[groups == "b"], x[groups == "c"]).statistic)
    assert 0.001 <= p <= 1.0


def test_pqn_removes_a_whole_run_dilution_that_percent_scaling_cannot():
    from pygcgc.core.batch import Feature
    from pygcgc.core.cluster import apply_pqn
    # Three runs of the same four-cluster profile; run c holds twice the amount of everything, but a big new
    # cluster (the last) dominates it, so percent-of-total makes the other clusters look lower in run c.
    base = np.array([10.0, 20.0, 30.0, 40.0])
    rows = {"a": base, "b": base * 1.0, "c": np.concatenate([base[:3] * 2.0, [400.0]])}
    total = {r: v.sum() for r, v in rows.items()}
    feats = [Feature(id=j, ref_t1=0.0, ref_t2=0.0, volumes={r: 100 * v[j] / total[r] for r, v in rows.items()})
             for j in range(4)]
    names = list(rows)
    fixed = apply_pqn(feats, names)
    # After PQN the three unchanged clusters agree across runs even though the fourth moved.
    for j in range(3):
        vals = [f for f in fixed if f.id == j][0].volumes
        assert vals["a"] == pytest.approx(vals["b"]) and vals["c"] == pytest.approx(vals["a"], rel=1e-6)
    raw = [f for f in feats if f.id == 0][0].volumes
    assert raw["c"] != pytest.approx(raw["a"], rel=0.05)                # ...which percent scaling alone did not manage


def test_pqn_flag_on_cluster_features_and_names_in_the_peak_csv(tmp_path):
    tables, names, _ = _two_groups(3)
    pooled = PeakTable.concat(tables)
    result = cluster_peaks(pooled)
    plain = cluster_features(pooled, result, relative=True)
    scaled = cluster_features(pooled, result, relative=True, pqn=True)
    assert len(scaled) == len(plain) and set(scaled[0].volumes) == set(names)
    path = write_clusters_csv(pooled, result, tmp_path / "c.csv", names={0: "alkylbenzenes"})
    lines = path.read_text().splitlines()
    assert lines[0].endswith("cluster_name") and any(l.endswith(",alkylbenzenes") for l in lines[1:])


def test_permanova_handles_groups_of_identical_runs_without_dividing_by_zero():
    X = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 1.0]])
    f, p = permanova(X, ["a", "a", "b", "b"], seed=0)
    assert f == np.inf and 0 < p <= 1


def test_permanova_separates_real_group_differences_and_not_shuffled_ones():
    # Six runs a side: five against five has only 252 distinct splits, too few for a small p.
    tables, names, groups = _two_groups(6)
    pooled = PeakTable.concat(tables)
    X = feature_table_matrix(cluster_features(pooled, cluster_peaks(pooled), relative=True), names)
    _, p_real = permanova(X, groups, seed=1)
    assert p_real <= 0.01
    shuffled = np.random.default_rng(7).permutation(groups)
    _, p_null = permanova(X, shuffled, seed=1)
    assert p_null > 0.05
    assert permanova(X, groups, seed=1) == permanova(X, groups, seed=1)   # seeded, so repeatable


# ------------------------------------------------------ reading a clustering
def test_label_image_paints_each_blobs_footprint_and_restores_wrapped_rows():
    from types import SimpleNamespace
    from pygcgc.core.cluster import cluster_label_image
    # 6 plain rows, 4 modulations, 2 wrap rows: extended image has 8 rows.
    ext = np.zeros((8, 4), dtype=int)
    ext[3:5, 1:3] = 1          # blob 1, in the canonical band
    ext[6, 0] = 2              # blob 2, in the wrapped band: plain row 0, column 1
    blobs = SimpleNamespace(labels=ext, wrap_rows=2)
    img = cluster_label_image(blobs, [5, 7])
    assert img.shape == (6, 4)
    assert (img[3:5, 1:3] == 5).all() and img[0, 1] == 7
    assert (img >= 0).sum() == 4 + 1                      # nothing else was painted
    assert cluster_label_image(blobs, [-1, -1]).max() == -1


def test_merge_joins_clusters_and_renumbers_by_volume():
    from pygcgc.core.cluster import merge_clusters
    kind = np.array([0] * 30 + [1] * 20)
    table = _table(kind, "R1")
    result = cluster_peaks(table, ClusterSettings(n_clusters=6, min_cluster_size=1))    # over-split on purpose
    assert len(result) == 6
    merged = merge_clusters(table, result, range(len(result)))
    assert len(merged) == 1 and (merged.labels == 0).sum() == (result.labels >= 0).sum()
    family0 = [k for k in range(len(result)) if result.top_ions(k, 1) == [57]]
    two = merge_clusters(table, result, family0)
    assert len(two.labels[kind == 0][two.labels[kind == 0] >= 0].tolist() and set(two.labels[kind == 0])) >= 1
    assert len(set(two.labels[kind == 0])) == 1 and len(two) == len(result) - len(family0) + 1
    with pytest.raises(ValueError):
        merge_clusters(table, result, [0])


def test_class_composition_lists_the_majority_class_first():
    from pygcgc.core.cluster import class_composition
    comp = class_composition([0, 0, 0, 1, 1, -1], ["alkane", "alkane", None, "PAH", "PAH", "PAH"])
    assert comp[0][0] == ("alkane", pytest.approx(2 / 3)) and comp[0][1][0] == "unclassified"
    assert comp[1] == [("PAH", 1.0)] and -1 not in comp


def test_adjusted_rand_is_one_for_the_same_partition_and_near_zero_for_unrelated():
    from pygcgc.core.cluster import adjusted_rand
    a = np.repeat([0, 1, 2], 100)
    assert adjusted_rand(a, np.repeat([5, 3, 9], 100)) == pytest.approx(1.0)          # names do not matter
    assert abs(adjusted_rand(a, np.random.default_rng(0).permutation(a))) < 0.05


def test_threshold_scan_gives_more_clusters_as_the_threshold_tightens():
    from pygcgc.core.cluster import threshold_scan
    kind = np.random.default_rng(1).integers(0, 2, 120)
    rows = threshold_scan(_table(kind, "R1"), ClusterSettings(), [0.02, 0.3, 0.6])
    assert [r["threshold"] for r in rows] == [0.02, 0.3, 0.6]
    assert rows[1]["n_clusters"] == 2 and rows[1]["unassigned_percent"] == pytest.approx(0.0)
    assert rows[2]["n_clusters"] <= rows[1]["n_clusters"]              # merging further never adds clusters here


def test_stability_is_high_for_clear_groups_and_repeatable():
    from pygcgc.core.cluster import stability
    table = _table(np.random.default_rng(2).integers(0, 2, 150), "R1")
    s = stability(table, ClusterSettings(), n_boot=8, seed=3)
    assert s["n_boot"] == 8 and s["mean"] > 0.95 and len(s["values"]) == 8
    assert stability(table, ClusterSettings(), n_boot=8, seed=3) == s
    with pytest.raises(ValueError):                                     # one family: nothing to compare
        stability(_table(np.zeros(40, dtype=int), "R1"), ClusterSettings(), n_boot=3)


# ------------------------------------------------------- retention-aware (density)
def _bands():
    """Family 0 in two retention bands and family 1 in a third: 40 peaks in each band."""
    kind = np.array([0] * 80 + [1] * 40)
    table = _table(kind, "R1", seed=5)
    rng = np.random.default_rng(6)
    table.t1 = np.concatenate([rng.uniform(300, 500, 40), rng.uniform(2200, 2400, 40), rng.uniform(1200, 1400, 40)])
    table.t2 = rng.uniform(1.0, 2.0, 120)
    return table, np.array([0] * 40 + [1] * 40 + [2] * 40)


def test_density_clustering_uses_retention_to_split_peaks_the_spectra_cannot():
    pytest.importorskip("sklearn")
    table, band = _bands()
    assert len(cluster_peaks(table)) == 2                                        # by spectrum alone: two families
    result = cluster_peaks(table, ClusterSettings(method="density", density_selection="eom"))
    assert len(result) == 3 and (result.labels >= 0).all()
    assert _same_partition(result.labels, band)                                  # ...and the two family-0 bands come apart
    assert result.tree.shape == (0, 4) and np.isnan(result.silhouette).all()


def test_density_clustering_with_no_retention_weight_matches_the_spectral_answer():
    pytest.importorskip("sklearn")
    table, band = _bands()
    result = cluster_peaks(table, ClusterSettings(method="density", density_selection="eom", retention_weight=0.0))
    kind = np.array([0] * 80 + [1] * 40)
    assert len(result) == 2 and _same_partition(result.labels, kind)


def test_density_clusters_are_numbered_by_volume_and_leave_scatter_as_noise():
    pytest.importorskip("sklearn")
    table, band = _bands()
    table.volume = np.where(band == 2, 5e6, 1e5).astype(float)
    scatter = _table(np.zeros(5, dtype=int), "R1", seed=8)                       # family-0 spectra, each far from the rest
    scatter.t1 = np.array([800.0, 1800.0, 2900.0, 3500.0, 4100.0])
    scatter.t2 = np.full(5, 1.5)
    pooled = PeakTable.concat([table, scatter])
    result = cluster_peaks(pooled, ClusterSettings(method="density", density_selection="eom"))
    assert result.top_ions(0, 1) == [65]                                         # the heavy band (family 1) is cluster 0
    assert (result.labels[-3:] == -1).all()                                      # far from every band in retention: noise
    placed = result.assign(table)                                                # by spectrum alone: the two family-0
    assert (placed[80:] == 0).all() and set(placed[:80]) <= {1, 2}               # bands share a spectrum, so either will do


def test_an_unknown_method_is_refused():
    with pytest.raises(ValueError, match="method"):
        cluster_peaks(_table(np.zeros(10, dtype=int), "R1"), ClusterSettings(method="magic"))


# ------------------------------------------------------------------- 14 eV
def _with_soft(table, kind, confusable):
    """Add 14 eV spectra. Where `confusable`, both families share one 70 eV spectrum but keep different 14 eV ones."""
    rng = np.random.default_rng(9)
    spectra = table.spectra.copy()
    soft = np.zeros((len(table), 40), dtype=np.float32)
    soft_mz = np.arange(90, 130)
    for i, k in enumerate(kind):
        soft[i, (100 if k == 0 else 120) - 90] = 1e4 * (1 + 0.1 * rng.standard_normal())       # the molecular ion
        soft[i, (101 if k == 0 else 121) - 90] = 5e3
        soft[i, 95 - 90], soft[i, 110 - 90], soft[i, 125 - 90] = 3e3, 3e3, 3e3
        if confusable:
            spectra[i] = _spectra(np.array([0]), rng)[0]
    return PeakTable(**{**table.__dict__, "spectra": spectra, "soft_spectra": soft, "soft_mz": soft_mz})


def test_the_14_eV_block_separates_peaks_whose_70_eV_spectra_match():
    kind = np.array([0] * 20 + [1] * 20)
    table = _with_soft(_table(kind, "R1"), kind, confusable=True)
    assert len(cluster_peaks(table)) == 1                                              # 70 eV alone cannot tell them apart
    result = cluster_peaks(table, ClusterSettings(soft_weight=1.0))
    assert len(result) == 2 and _same_partition(result.labels, kind)
    assert result.soft_used is not None and result.centroids.shape[1] == result.mz_used.size + result.soft_used.size


def test_soft_spectra_are_needed_to_assign_with_a_soft_clustering():
    kind = np.array([0] * 20 + [1] * 20)
    table = _with_soft(_table(kind, "R1"), kind, confusable=True)
    result = cluster_peaks(table, ClusterSettings(soft_weight=1.0))
    assert _same_partition(result.assign(table), kind)
    with pytest.raises(ValueError, match="14 eV"):
        result.assign(PeakTable(**{**table.__dict__, "soft_spectra": None, "soft_mz": None}))


def test_pooling_keeps_soft_spectra_only_when_every_run_has_them():
    kind = np.zeros(5, dtype=int)
    a = _with_soft(_table(kind, "A"), kind, confusable=False)
    b = _with_soft(_table(kind, "B"), kind, confusable=False)
    assert PeakTable.concat([a, b]).soft_spectra.shape == (10, 40)
    assert PeakTable.concat([a, _table(kind, "C")]).soft_spectra is None


def test_csv_has_one_row_per_peak(tmp_path):
    table = _table(np.array([0] * 10 + [1] * 10), "R1")
    path = write_clusters_csv(table, cluster_peaks(table), tmp_path / "c.csv")
    lines = path.read_text().splitlines()
    assert len(lines) == 21 and lines[0].startswith("run,blob,t1_min")
