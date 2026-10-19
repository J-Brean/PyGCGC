"""Reader for ANDI-MS (AIA) netCDF files, the .cdf export from ChromSpace
and most other GC-MS software.

The format stores every scan's centroids end to end in two long arrays
(masses and intensities). `point_count` says how many belong to each scan,
and `scan_index` says where each scan starts. We rebuild the per-scan
structure from those, round masses to unit m/z and pack everything into
a sparse matrix.
"""

from __future__ import annotations

import threading
import warnings
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from netCDF4 import Dataset

from ..run import Run

# Fraction of centroids allowed to sit far (>0.3 Da) from a whole mass before
# we warn that the file probably isn't unit-resolution data.
OFF_NOMINAL_TOLERANCE = 0.02

# The netCDF/HDF5 C libraries are not guaranteed to be thread-safe, so file
# access is serialised. Building the sparse matrix happens outside the lock,
# so reading several files in threads still gains something.
_NETCDF_LOCK = threading.Lock()


def _var(ds, name, required=True):
    if name in ds.variables:
        v = ds.variables[name]
        v.set_auto_mask(False)          # plain arrays; scale factors still applied
        return np.asarray(v[:])
    if required:
        raise ValueError(f"ANDI file is missing the '{name}' variable")
    return None


def read_andi(path: str | Path) -> Run:
    path = Path(path)
    with _NETCDF_LOCK, Dataset(path, "r") as ds:
        times = _var(ds, "scan_acquisition_time").astype(np.float64)
        masses = _var(ds, "mass_values").astype(np.float64)
        inten = _var(ds, "intensity_values").astype(np.float64)
        counts = _var(ds, "point_count", required=False)
        starts = _var(ds, "scan_index", required=False)
        attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
        time_units = getattr(ds.variables["scan_acquisition_time"], "units", "seconds")

    n_scans = times.size
    if times.ndim != 1 or masses.ndim != 1 or inten.shape != masses.shape:
        raise ValueError("ANDI scan times, masses and intensities must be matching one-dimensional arrays")

    # Work out which scan each stored point belongs to.
    if counts is None and starts is None:
        raise ValueError("ANDI file has neither 'point_count' nor 'scan_index'")
    if counts is None:
        ends = np.append(starts[1:], masses.size)
        counts = ends - starts
    if counts.shape != (n_scans,) or not np.isfinite(counts).all() or np.any(counts < 0) or np.any(counts != np.floor(counts)):
        raise ValueError("ANDI point_count must contain one non-negative integer per scan")
    counts = counts.astype(np.int64)
    total = int(counts.sum())
    if starts is not None:
        if (starts.shape != (n_scans,) or not np.isfinite(starts).all()
                or np.any(starts < 0) or np.any(starts != np.floor(starts))):
            raise ValueError("ANDI scan_index must contain one non-negative integer per scan")
        starts = starts.astype(np.int64)
        if np.any(starts + counts > masses.size) or np.any(starts[1:] < (starts+counts)[:-1]):
            raise ValueError("ANDI scan_index and point_count overlap or exceed the stored points")
    elif total != masses.size and not (total == 0 and masses.size == 1):
        raise ValueError(
            f"point counts add up to {total:,} but the file holds {masses.size:,} points")
    scan_of_point = np.repeat(np.arange(n_scans), counts)
    if starts is not None:
        # Respect offsets: some valid files include unused padding between scans.
        offsets = np.arange(total) - np.repeat(np.cumsum(counts)-counts, counts)
        selected = starts[scan_of_point] + offsets
        masses, inten = masses[selected], inten[selected]
    elif total == 0:
        # Classic netCDF needs one padding point for an otherwise empty run.
        masses, inten = masses[:0], inten[:0]

    # Some exports write minutes despite the spec asking for seconds.
    if str(time_units).strip().lower().startswith("min"):
        times = times * 60.0

    # Drop junk points before binning.
    good = np.isfinite(masses) & np.isfinite(inten) & (inten > 0)
    masses, inten, scan_of_point = masses[good], inten[good], scan_of_point[good]

    unit = np.rint(masses).astype(np.int64)
    off = np.abs(masses - unit) > 0.3
    if masses.size and off.mean() > OFF_NOMINAL_TOLERANCE:
        warnings.warn(
            f"{path.name}: {off.mean():.1%} of masses are more than 0.3 Da from a whole "
            "number. This looks like high-resolution data; values were rounded to unit m/z.",
            stacklevel=2)

    if unit.size:
        mz_min, mz_max = int(unit.min()), int(unit.max())
        # Respect the declared scan range if the file gives one and it is wider.
        mz_min = min(mz_min, int(np.floor(attrs.get("global_mass_min", mz_min))))
        mz_max = max(mz_max, int(np.ceil(attrs.get("global_mass_max", mz_max))))
    else:
        mz_min = int(np.floor(attrs.get("global_mass_min", 0)))
        mz_max = int(np.ceil(attrs.get("global_mass_max", mz_min)))
    mz_axis = np.arange(mz_min, mz_max + 1)

    # COO -> CSR sums duplicates, which merges two centroids that round to the
    # same unit mass within one scan.
    data = sp.coo_matrix((inten.astype(np.float32), (scan_of_point, unit - mz_min)),
                         shape=(n_scans, mz_axis.size)).tocsr()
    data.sum_duplicates()

    # Scan times must increase. Drop exact repeats (seen in some exports).
    finite = np.isfinite(times)
    previous = np.concatenate([[-np.inf], np.maximum.accumulate(np.where(finite, times, -np.inf))[:-1]])
    keep = finite & (times > previous) if n_scans else finite
    if not keep.all():
        warnings.warn(f"{path.name}: dropped {np.sum(~keep)} scans with non-increasing times",
                      stacklevel=2)
        times, data = times[keep], data[keep]

    meta = {
        "format": "ANDI-MS",
        "experiment_title": attrs.get("experiment_title", ""),
        "experiment_type": attrs.get("experiment_type", ""),
        "ionisation": attrs.get("test_ionization_mode", ""),
        "detector": attrs.get("test_detector_type", ""),
        "comments": attrs.get("administrative_comments", ""),
        "attributes": attrs,
    }
    return Run(times=times, mz=mz_axis, data=data, name=path.stem, source=path, meta=meta)


def write_andi(run: Run, path: str | Path, comment: str = "") -> Path:
    """Write a run as a minimal ANDI-MS file (enough for PyGCGC and most
    other readers). Used for cropped test fixtures and for exporting
    cleaned runs."""
    path = Path(path)
    coo = run.data.tocoo()
    order = np.lexsort((coo.col, coo.row))
    rows, cols, vals = coo.row[order], coo.col[order], coo.data[order]
    counts = np.bincount(rows, minlength=run.n_scans).astype(np.int32)
    starts = np.concatenate([[0], np.cumsum(counts)[:-1]]).astype(np.int32)
    tic = np.bincount(rows, weights=vals, minlength=run.n_scans)
    with _NETCDF_LOCK, Dataset(path, "w", format="NETCDF3_CLASSIC") as ds:
        ds.dataset_completeness = "C1+C2"
        ds.ms_template_revision = "1.0.1"
        ds.experiment_type = "Centroided Mass Spectrum"
        ds.experiment_title = run.name
        ds.administrative_comments = comment or f"Written by PyGCGC from {run.source or run.name}"
        ds.global_mass_min = float(run.mz[0])
        ds.global_mass_max = float(run.mz[-1])
        ds.createDimension("scan_number", run.n_scans)
        ds.createDimension("point_number", max(int(vals.size), 1))
        v = ds.createVariable("scan_acquisition_time", "f8", ("scan_number",))
        v.units = "seconds"
        v[:] = run.times
        ds.createVariable("total_intensity", "f8", ("scan_number",))[:] = tic
        ds.createVariable("point_count", "i4", ("scan_number",))[:] = counts
        ds.createVariable("scan_index", "i4", ("scan_number",))[:] = starts
        mv = ds.createVariable("mass_values", "f4", ("point_number",))
        iv = ds.createVariable("intensity_values", "f4", ("point_number",))
        if vals.size:
            mv[:] = run.mz[cols].astype(np.float32)
            iv[:] = vals.astype(np.float32)
    return path


def crop_run(run: Run, t_start: float, t_end: float, mz_min: int | None = None,
             mz_max: int | None = None) -> Run:
    """A time (and optionally mass) window of a run, for small test files."""
    keep = (run.times >= t_start) & (run.times <= t_end)
    lo = 0 if mz_min is None else int(mz_min) - int(run.mz[0])
    hi = run.n_mz if mz_max is None else int(mz_max) - int(run.mz[0]) + 1
    lo, hi = max(lo, 0), min(hi, run.n_mz)
    data = run.data[np.flatnonzero(keep)][:, lo:hi]
    return Run(times=run.times[keep], mz=run.mz[lo:hi], data=data,
               name=f"{run.name}_crop", source=run.source, meta=dict(run.meta))
