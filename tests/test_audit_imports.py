"""Import failures must not replace working project inputs or lose empty data."""

import os
import numpy as np
import pytest
from netCDF4 import Dataset
from pygcgc.core.targets import Target, load_targets, write_targets_csv
from pygcgc.core.io.andi import read_andi, write_andi
from tests.conftest import make_run
from tests.test_andi import write_andi as write_fixture, SCANS


def test_append_target_respects_minimal_reordered_header_and_extra_columns(tmp_path):
    import csv
    path = tmp_path / "targets.csv"
    path.write_text("note,t2_expected_s,quant_ion,compound,t1_expected_s\nkeep me,1,50,A,20\n")
    added = Target("B", quant_ion=51, t1=30, t2=2, ri=1200)
    write_targets_csv([added], path, append=True)
    actual = load_targets(path)
    assert [t.compound for t in actual] == ["A", "B"]
    assert actual[1] == added
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert rows[0]["note"] == "keep me"


def test_append_duplicate_target_keeps_existing_file_unchanged(tmp_path):
    path = write_targets_csv([Target("A", quant_ion=50)], tmp_path / "targets.csv")
    before = path.read_bytes()
    with pytest.raises(ValueError, match="duplicate"):
        write_targets_csv([Target("A", quant_ion=51)], path, append=True)
    assert path.read_bytes() == before


def test_target_write_failure_keeps_existing_file_unchanged(tmp_path, monkeypatch):
    import pygcgc.core.targets as module
    path = write_targets_csv([Target("A", quant_ion=50)], tmp_path / "targets.csv")
    before = path.read_bytes()
    def fail(*args):
        raise PermissionError("locked file")
    monkeypatch.setattr(module.os, "replace", fail)
    with pytest.raises(PermissionError):
        write_targets_csv([Target("B", quant_ion=51)], path, append=True)
    assert path.read_bytes() == before
    assert list(tmp_path.iterdir()) == [path]


def test_empty_andi_roundtrip_preserves_declared_mass_axis(tmp_path):
    run = make_run([], mz=(50, 54))
    actual = read_andi(write_andi(run, tmp_path/"blank.cdf"))
    np.testing.assert_array_equal(actual.mz, run.mz)
    np.testing.assert_array_equal(actual.times, run.times)
    assert actual.nnz == 0


def test_andi_backward_times_cannot_survive_the_duplicate_scan_filter(tmp_path):
    path = tmp_path/"backward.cdf"
    write_fixture(path, [10., 1., 2.], SCANS)
    with pytest.warns(UserWarning, match="non-increasing"):
        actual = read_andi(path)
    np.testing.assert_array_equal(actual.times, [10.])
    np.testing.assert_allclose(actual.tic, [30.])


def test_andi_respects_scan_offsets_in_a_padded_file(tmp_path):
    path = tmp_path/"padded.cdf"
    write_fixture(path, [0., 1., 2.], SCANS)
    with Dataset(path, "a") as ds:
        ds.variables["point_count"][:] = [1, 1, 1]
        ds.variables["scan_index"][:] = [1, 3, 5]
    actual = read_andi(path)
    np.testing.assert_allclose(actual.tic, [20., 1., 3.])


@pytest.mark.parametrize("field,value", [
    ("quant_ion", "50.5"), ("quant_ion", "0"), ("t1_expected_s", "nan"),
    ("t2_expected_s", "inf"), ("t1_window_s", "-1"), ("t2_window_s", "0"),
    ("stock_factor", "nan"), ("stock_factor", "-2"), ("qual_ratios", "-0.1"),
    ("qual_ions", "51.3"), ("role", "internalstandard"), ("ri", "inf")])
def test_target_import_rejects_invalid_numbers_and_roles(tmp_path, field, value):
    import csv
    row = dict(compound="IS", role="istd", quant_ion="50", t1_expected_s="30", t2_expected_s="1",
               t1_window_s="12", t2_window_s=".25", stock_factor="1", qual_ions="51", qual_ratios=".1", ri="")
    row[field] = value
    path = tmp_path/"bad.csv"
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    with pytest.raises(ValueError):
        load_targets(path)


def test_target_cannot_normalise_to_a_row_that_is_not_an_internal_standard(tmp_path):
    path = write_targets_csv([Target("A", quant_ion=50, t1=20, t2=1),
                              Target("B", quant_ion=51, t1=40, t2=2, istd="A")], tmp_path/"bad.csv")
    with pytest.raises(ValueError, match="internal standard"):
        load_targets(path)


@pytest.mark.parametrize("text", ["level,A\nL1,nan\n", "level,A\nL1,-1\n", "level,A\nL1,1\nL1,2\n",
                                  "level,A,A\nL1,1,2\n"])
def test_levels_cannot_silently_use_invalid_or_overwritten_amounts(tmp_path, text):
    from pygcgc.core.targets import load_levels
    path = tmp_path/"bad.csv"
    path.write_text(text)
    with pytest.raises(ValueError):
        load_levels(path)


@pytest.mark.parametrize("times,masses", [([0., np.nan], [50]), ([0., np.inf], [50]),
                                           ([[0.], [1.]], [50]), ([0., 1.], [50.5])])
def test_run_rejects_axes_that_would_corrupt_folding(times, masses):
    from pygcgc.core.run import Run
    with pytest.raises(ValueError):
        Run(times, masses, np.zeros((2, 1)))


@pytest.mark.parametrize("kwargs", [dict(period=np.nan), dict(period=np.inf), dict(period=0),
                                    dict(period=.001), dict(period=4., phase=np.nan),
                                    dict(period=4., n_rows=0), dict(period=4., n_rows=3.5)])
def test_invalid_fold_settings_fail_before_integer_conversion_or_division(kwargs):
    from pygcgc.core.fold import FoldSpec
    with pytest.raises(ValueError):
        FoldSpec(**kwargs).resolved(np.arange(10)/50)


def test_failed_target_import_preserves_all_previous_inputs(qapp, tmp_path):
    from pygcgc.gui.jobs import JobManager
    from pygcgc.gui.project import Project
    project = Project(JobManager())
    good = write_targets_csv([Target("IS", role="istd", quant_ion=50, t1=30, t2=1)], tmp_path/"good.csv")
    project.load_targets(good, estimate=False)
    before = (project.targets_path, list(project.anchors), list(project.quant.targets))
    bad = tmp_path/"bad.csv"
    bad.write_text("compound,role,quant_ion,t1_expected_s,t2_expected_s,qual_ions,qual_ratios\n"
                   "IS,istd,50,30,1,51;52,0.2\n")
    with pytest.raises(ValueError):
        project.load_targets(bad, estimate=False)
    assert (project.targets_path, project.anchors, project.quant.targets) == before


def test_target_only_table_can_be_imported_before_adding_standards(qapp, tmp_path):
    from pygcgc.gui.jobs import JobManager
    from pygcgc.gui.project import Project
    project = Project(JobManager())
    table = write_targets_csv([Target("Target", quant_ion=50, t1=30, t2=1)], tmp_path/"targets.csv")
    project.load_targets(table, estimate=False)
    assert [t.compound for t in project.quant.targets] == ["Target"]
    assert project.anchors == []
    saved = tmp_path/"project.pygcgc.json"
    project.save(saved)
    restored = Project(JobManager())
    restored.load(saved)
    assert [t.compound for t in restored.quant.targets] == ["Target"]
    assert restored.anchors == []
