"""🗂️ Batch pipeline: caching, the process pool, feature table, QC and
cross-run statistics."""

import numpy as np
import pytest

from pygcgc import (BatchConfig, BatchSettings, Feature, pca, run_batch, volcano, westgard_flags)
import pygcgc.core.batch as batch_mod


def _count_calls(monkeypatch, name):
    calls = {"n": 0}
    original = getattr(batch_mod, name)

    def wrapped(*a, **k):
        calls["n"] += 1
        return original(*a, **k)

    monkeypatch.setattr(batch_mod, name, wrapped)
    return calls


def test_pipeline_runs_and_caches(cal_folder, tmp_path, monkeypatch):
    bg_calls = _count_calls(monkeypatch, "remove_background")
    blob_calls = _count_calls(monkeypatch, "find_blobs")
    cache_dir = tmp_path / "cache"

    result = run_batch(cal_folder, cache_dir=cache_dir, max_workers=1)
    assert result.errors == {}
    assert set(result.prepped) == {"CAL_L1", "CAL_L2", "CAL_L3", "CAL_L4", "CAL_L5", "S1", "S2"}
    assert bg_calls["n"] == 7 and blob_calls["n"] == 7
    assert result.report is not None
    assert len(result.report.results) > 0
    assert result.feature_table is not None and len(result.feature_table) > 0

    bg_calls["n"] = blob_calls["n"] = 0
    result2 = run_batch(cal_folder, cache_dir=cache_dir, max_workers=1)
    assert bg_calls["n"] == 0 and blob_calls["n"] == 0        # everything loaded from cache
    assert set(result2.prepped) == set(result.prepped)


def test_changing_settings_invalidates_only_that_step_onward(cal_folder, tmp_path, monkeypatch):
    bg_calls = _count_calls(monkeypatch, "remove_background")
    blob_calls = _count_calls(monkeypatch, "find_blobs")
    cache_dir = tmp_path / "cache"
    run_batch(cal_folder, cache_dir=cache_dir, max_workers=1)

    bg_calls["n"] = blob_calls["n"] = 0
    from pygcgc import PeakSettings
    changed = BatchConfig(settings=BatchSettings(peaks=PeakSettings(snr_min=12)))
    run_batch(cal_folder, config=changed, cache_dir=cache_dir, max_workers=1)
    assert bg_calls["n"] == 0          # background settings unchanged: still cached
    assert blob_calls["n"] == 7        # blob settings changed: recomputed for every run


def test_batch_blob_detection_uses_the_raw_storage_floor(cal_folder, tmp_path, monkeypatch):
    seen = []
    original = batch_mod.find_blobs

    def capture(clean, fold, settings, **kwargs):
        seen.append(kwargs.get("floor"))
        return original(clean, fold, settings, **kwargs)

    monkeypatch.setattr(batch_mod, "find_blobs", capture)
    result = run_batch(cal_folder, cache_dir=tmp_path / "cache", max_workers=1)
    expected = [batch_mod.storage_threshold(p.run) / 3.0 for p in result.prepped.values()]
    assert sorted(seen) == pytest.approx(sorted(expected))


def test_one_run_failing_does_not_stop_the_batch(cal_folder, tmp_path, monkeypatch):
    def boom(run, *a, **k):
        if run.name == "S1":
            raise RuntimeError("synthetic failure")
        return original(run, *a, **k)

    original = batch_mod.remove_background
    monkeypatch.setattr(batch_mod, "remove_background", boom)
    result = run_batch(cal_folder, cache_dir=tmp_path / "cache", max_workers=1)
    assert "S1" in result.errors
    assert "S2" in result.prepped and "CAL_L1" in result.prepped
    assert len(result.prepped) == 6


def test_batch_config_round_trip(tmp_path):
    cfg = BatchConfig(period=4.0, reference="CAL_L4", group_a=("S1",), group_b=("S2",))
    path = tmp_path / "pipeline.json"
    cfg.save(path)
    loaded = BatchConfig.load(path)
    assert loaded == cfg


# ------------------------------------------------------------------ QC
def test_westgard_1_3s_and_1_2s():
    values = np.array([10.0, 10.1, 9.9, 10.2, 9.8, 10.0, 9.9, 10.1, 10.0, 10.0])
    values_out = np.append(values, 16.0)          # a clear outlier
    flags = westgard_flags(values_out)
    at_outlier = [f for f in flags if f.index == 10]
    assert at_outlier and at_outlier[0].rule in ("1-2s", "1-3s")


def test_westgard_2_2s_rule():
    mean, sd = 10.0, 1.0
    values = np.array([10, 10, 10, 10, 12.5, 12.5, 10, 10])
    flags = westgard_flags(values, mean=mean, sd=sd)
    assert any(f.rule == "2-2s" and f.index == 5 for f in flags)


def test_westgard_10x_rule():
    mean, sd = 10.0, 1.0
    values = np.full(12, 10.3)       # every point on the same (high) side, within 1 SD
    flags = westgard_flags(values, mean=mean, sd=sd)
    assert any(f.rule == "10x" for f in flags)


def test_westgard_quiet_series_mostly_clean():
    rng = np.random.default_rng(0)
    values = rng.normal(10.0, 0.3, 20)
    flags = westgard_flags(values)
    # A little noise can trip 1-2s occasionally; the strict multi-point
    # rules should not fire on genuinely unremarkable data.
    assert not any(f.rule in ("2-2s", "R-4s", "4-1s", "10x") for f in flags)


# ------------------------------------------------------------------ stats
def test_pca_separates_distinct_runs():
    features = [Feature(id=i, ref_t1=0.0, ref_t2=0.0,
                        volumes={"r1": 100 + i, "r2": 100 + i, "r3": 500 + i * 5})
               for i in range(6)]
    res = pca(features, ["r1", "r2", "r3"], n_components=2)
    assert res.scores.shape == (3, 2)
    assert res.explained_variance[0] > 0.99
    assert abs(res.scores[0, 0] - res.scores[1, 0]) < abs(res.scores[0, 0] - res.scores[2, 0])


def test_volcano_finds_the_planted_difference():
    rng = np.random.default_rng(1)
    features = []
    for i in range(20):
        a = {f"a{k}": max(rng.normal(100, 5), 0) for k in range(4)}
        shifted = i < 3
        b = {f"b{k}": max(rng.normal(300 if shifted else 100, 5), 0) for k in range(4)}
        features.append(Feature(id=i, ref_t1=0.0, ref_t2=0.0, volumes={**a, **b}))
    points = volcano(features, [f"a{k}" for k in range(4)], [f"b{k}" for k in range(4)])
    significant = {p.feature_id for p in points if p.q_value < 0.05}
    assert significant == {0, 1, 2}
    for p in points:
        if p.feature_id in significant:
            assert p.log2_fold_change < 0     # group b higher than group a


def test_cli_batch(cal_folder, tmp_path, capsys):
    from pygcgc.cli import main
    out = tmp_path / "batch_out"
    rc = main(["batch", str(cal_folder), "--outdir", str(out), "--workers", "1"])
    text = capsys.readouterr().out
    assert rc == 0
    assert "7 run(s) prepared, 0 failed" in text
    for name in ("targets.csv", "calibrations.csv", "all_signals.csv", "internal_standards.csv",
                "feature_table.csv"):
        assert (out / name).exists()
    assert (cal_folder / ".pygcgc_cache").is_dir()


def test_cli_batch_saves_and_loads_config(cal_folder, tmp_path):
    from pygcgc.cli import main
    config_path = tmp_path / "pipeline.json"
    out = tmp_path / "batch_out"
    assert main(["batch", str(cal_folder), "--outdir", str(out), "--workers", "1",
               "--no-quant", "--save-config", str(config_path)]) == 0
    assert config_path.exists()
    from pygcgc import BatchConfig
    saved = BatchConfig.load(config_path)
    assert saved.settings.run_quant is False

    out2 = tmp_path / "batch_out2"
    assert main(["batch", str(cal_folder), "--outdir", str(out2), "--workers", "1",
               "--config", str(config_path)]) == 0
    assert not (out2 / "targets.csv").exists()          # quant was off in the loaded config


def test_volcano_skips_features_with_too_few_nonzero_values():
    features = [Feature(id=0, ref_t1=0.0, ref_t2=0.0,
                        volumes={"a1": 100.0, "a2": 0.0, "b1": 100.0, "b2": 0.0})]
    assert volcano(features, ["a1", "a2"], ["b1", "b2"], min_group=2) == []


# ------------------------------------------------------------------ 📤 reports
def test_cli_batch_writes_provenance_json_always(cal_folder, tmp_path):
    import json

    from pygcgc.cli import main
    out = tmp_path / "batch_out"
    assert main(["batch", str(cal_folder), "--outdir", str(out), "--workers", "1"]) == 0
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["period"] == 4.0
    assert prov["version"]
    assert prov["edits"] == {}          # a scripted batch run has no interactive review to report


def test_cli_batch_report_all_writes_xlsx_html_pdf(cal_folder, tmp_path):
    from pygcgc.cli import main
    out = tmp_path / "batch_out"
    assert main(["batch", str(cal_folder), "--outdir", str(out), "--workers", "1",
               "--report", "all"]) == 0
    assert (out / "report.xlsx").exists()
    assert (out / "report.html").exists()
    assert (out / "report.pdf").exists()
    assert "PyGCGC" in (out / "report.html").read_text(encoding="utf-8")


def test_cli_batch_report_uses_config_groups_for_volcano(cal_folder, tmp_path):
    from pygcgc.cli import main
    config = BatchConfig(settings=BatchSettings(), group_a=("CAL_L1", "CAL_L2"),
                         group_b=("CAL_L4", "CAL_L5"))
    config_path = tmp_path / "pipeline.json"
    config.save(config_path)
    out = tmp_path / "batch_out"
    assert main(["batch", str(cal_folder), "--outdir", str(out), "--workers", "1",
               "--config", str(config_path), "--report", "html"]) == 0
    assert "Volcano" in (out / "report.html").read_text(encoding="utf-8")


def test_cli_quant_writes_provenance_and_optional_xlsx(cal_folder, tmp_path):
    import json

    from pygcgc.cli import main
    out = tmp_path / "quant_out"
    assert main(["quant", str(cal_folder), "--outdir", str(out), "--workers", "1", "--xlsx"]) == 0
    prov = json.loads((out / "provenance.json").read_text())
    assert prov["period"] == 4.0
    assert (out / "report.xlsx").exists()
