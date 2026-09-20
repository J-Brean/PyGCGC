import warnings

import numpy as np
import pytest
from netCDF4 import Dataset

from pygcgc import read_run


def write_andi(path, times, per_scan, time_units="seconds", with_counts=True,
               with_index=True, mass_scale=None):
    """per_scan: list of (masses, intensities) per scan."""
    masses = np.concatenate([np.asarray(m, float) for m, _ in per_scan]) if per_scan else np.array([])
    inten = np.concatenate([np.asarray(i, float) for _, i in per_scan]) if per_scan else np.array([])
    counts = np.array([len(m) for m, _ in per_scan], dtype=np.int32)
    with Dataset(path, "w", format="NETCDF3_CLASSIC") as ds:
        ds.global_mass_min = float(masses.min())
        ds.global_mass_max = float(masses.max())
        ds.createDimension("scan_number", len(times))
        ds.createDimension("point_number", masses.size)
        v = ds.createVariable("scan_acquisition_time", "f8", ("scan_number",))
        v.units = time_units
        v[:] = times
        if with_counts:
            ds.createVariable("point_count", "i4", ("scan_number",))[:] = counts
        if with_index:
            ds.createVariable("scan_index", "i4", ("scan_number",))[:] = np.concatenate([[0], np.cumsum(counts)[:-1]])
        if mass_scale:
            mv = ds.createVariable("mass_values", "i2", ("point_number",))
            mv.scale_factor = mass_scale
            mv[:] = masses
        else:
            ds.createVariable("mass_values", "f4", ("point_number",))[:] = masses
        ds.createVariable("intensity_values", "f4", ("point_number",))[:] = inten


SCANS = [([50.0, 51.0], [10.0, 20.0]), ([50.0], [5.0]), ([52.0, 53.0, 50.0], [1.0, 2.0, 3.0])]


def test_basic_read(tmp_path):
    p = tmp_path / "a.cdf"
    write_andi(p, [1.0, 1.02, 1.04], SCANS)
    r = read_run(p)
    assert r.n_scans == 3
    assert list(r.mz) == [50, 51, 52, 53]
    np.testing.assert_allclose(r.tic, [30, 5, 6])
    np.testing.assert_allclose(r.xic(50), [10, 5, 3])
    assert r.name == "a"


def test_scan_index_only(tmp_path):
    p = tmp_path / "b.cdf"
    write_andi(p, [1.0, 1.02, 1.04], SCANS, with_counts=False)
    np.testing.assert_allclose(read_run(p).tic, [30, 5, 6])


def test_minutes_converted(tmp_path):
    p = tmp_path / "c.cdf"
    write_andi(p, [5.0, 5.001, 5.002], SCANS, time_units="Minutes")
    np.testing.assert_allclose(read_run(p).times, [300.0, 300.06, 300.12])


def test_centroids_rounded_and_merged(tmp_path):
    p = tmp_path / "d.cdf"
    # 50.1 and 49.9 both round to 50 and must be summed.
    write_andi(p, [0.0, 0.1], [([50.1, 49.9, 51.0], [1.0, 2.0, 4.0]), ([51.2], [7.0])])
    r = read_run(p)
    np.testing.assert_allclose(r.xic(50), [3.0, 0.0])
    np.testing.assert_allclose(r.xic(51), [4.0, 7.0])


def test_warns_for_high_resolution_masses(tmp_path):
    p = tmp_path / "e.cdf"
    scans = [([50.45, 128.4, 60.6], [1.0, 1.0, 1.0])] * 3
    write_andi(p, [0.0, 0.1, 0.2], scans)
    with pytest.warns(UserWarning, match="high-resolution"):
        read_run(p)


def test_scaled_integer_masses(tmp_path):
    p = tmp_path / "f.cdf"
    write_andi(p, [1.0, 1.02, 1.04], SCANS, mass_scale=0.05)
    assert list(read_run(p).mz) == [50, 51, 52, 53]


def test_missing_structure_is_an_error(tmp_path):
    p = tmp_path / "g.cdf"
    write_andi(p, [1.0, 1.02, 1.04], SCANS, with_counts=False, with_index=False)
    with pytest.raises(ValueError, match="point_count"):
        read_run(p)


def test_unknown_extension(tmp_path):
    with pytest.raises(ValueError, match="No reader"):
        read_run(tmp_path / "x.unsupported")


@pytest.mark.dummy
def test_dummy_files_consistent(dummy_dir, dummy_runs):
    for name, run in dummy_runs.items():
        with Dataset(dummy_dir / "data" / f"{name}.cdf") as ds:
            n_points = int(ds["point_count"][:].sum())
            tic_file = np.asarray(ds["total_intensity"][:])
        assert run.n_scans == 93000
        assert run.times[0] == pytest.approx(300.0)
        assert run.mz[0] == 35 and run.mz[-1] == 500
        assert run.nnz == n_points
        np.testing.assert_allclose(run.tic, tic_file, rtol=1e-5, atol=1e-2)


def test_write_and_crop_round_trip(tmp_path):
    from pygcgc import crop_run, write_andi
    from tests.conftest import make_run
    run = make_run([(62.0, 1.2, 50, 1e5), (122.0, 2.6, 52, 1e5)], rate=20.0)
    back = read_run(write_andi(run, tmp_path / "w.cdf"))
    np.testing.assert_allclose(back.times, run.times)
    np.testing.assert_allclose(back.data.toarray(), run.data.toarray(), rtol=1e-6, atol=1e-4)
    small = crop_run(run, 100.0, 150.0, 51, 53)
    assert small.times[0] >= 100.0 and small.times[-1] <= 150.0
    assert list(small.mz) == [51, 52, 53]
    np.testing.assert_allclose(small.xic(52), run.xic(52)[(run.times >= 100) & (run.times <= 150)])
    back = read_run(write_andi(small, tmp_path / "s.cdf"))
    assert list(back.mz) == [51, 52, 53]


def test_cli_crop_and_clean(tmp_path, capsys):
    from pygcgc import write_andi
    from pygcgc.cli import main
    from tests.conftest import make_run
    run = make_run([(62.0, 1.2, 50, 1e5)], rate=20.0)
    src = write_andi(run, tmp_path / "r.cdf")
    assert main(["crop", str(src), "--start", "0.5", "--end", "2", "-o", str(tmp_path / "c.cdf")]) == 0
    assert (tmp_path / "c.cdf").exists()
    assert main(["clean", str(src), "--period", "4", "--phase", "0", "--outdir", str(tmp_path / "out")]) == 0
    assert (tmp_path / "out" / "r_clean.cdf").exists()
    assert "Removed" in capsys.readouterr().out
