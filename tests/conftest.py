"""Shared test fixtures.

Two kinds of test data:

* Small synthetic runs built here, where every number is known exactly.
  These always run.
* The dummy data pack (synthetic GC×GC-ToF-MS files with ground truth).
  Point PYGCGC_DUMMY_DATA at the folder containing data/, target_table.csv
  and truth_per_file.csv. Tests that need it are skipped if it is missing.
"""

from __future__ import annotations

import csv
import os
import tempfile
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

from pygcgc import Run

# The GUI's results cache goes to a throwaway folder, never the user's own.
os.environ["PYGCGC_CACHE"] = tempfile.mkdtemp(prefix="pygcgc-test-cache-")

TRUE_PHASE = 0.55      # the dummy pack's hidden phase
PERIOD = 4.0


def make_run(peaks, period=4.0, phase=0.0, rate=20.0, t_start=0.0, duration=240.0,
             mz=(50, 54), sigma1=4.0, sigma2=0.08, noise=0.0, seed=0, name="synthetic"):
    """Build a small modulated run.

    peaks: list of (t1, t2, ion, area). t1 is the first-dimension apex
    (absolute seconds), t2 the true second-dimension time (seconds, may
    exceed the period). Each modulation k starts at t_start + phase + k*period
    and releases the fraction of the first-dimension peak that eluted during
    that modulation.
    """
    rng = np.random.default_rng(seed)
    times = t_start + np.arange(int(round(duration * rate))) / rate
    mz_axis = np.arange(mz[0], mz[1] + 1)
    dense = np.zeros((times.size, mz_axis.size))
    k_all = np.arange(-2, int(duration / period) + 3)
    starts = t_start + phase + k_all * period
    from scipy.special import erf
    for t1, t2, ion, area in peaks:
        cdf = lambda x: 0.5 * (1 + erf((x - t1) / (sigma1 * np.sqrt(2))))
        frac = cdf(starts + period) - cdf(starts)
        for s, f in zip(starts, frac):
            if f < 1e-6:
                continue
            apex = s + t2
            dense[:, ion - mz[0]] += (area * f / (sigma2 * np.sqrt(2 * np.pi)) / rate
                                      * np.exp(-0.5 * ((times - apex) / sigma2) ** 2))
    if noise:
        dense += rng.normal(0, noise, dense.shape).clip(0)
    dense[dense < 1e-9] = 0.0
    return Run(times=times, mz=mz_axis, data=sp.csr_matrix(dense), name=name)


@pytest.fixture
def synthetic_peaks():
    # Three compounds on different ions, one wrapping past the period.
    return [(60.0, 1.20, 50, 1e5), (120.0, 2.60, 51, 2e5), (180.0, 4.30, 52, 1.5e5)]


# ------------------------------------------------------------------ dummy pack
def _dummy_dir() -> Path | None:
    candidates = []
    if os.environ.get("PYGCGC_DUMMY_DATA"):
        candidates.append(Path(os.environ["PYGCGC_DUMMY_DATA"]))
    here = Path(__file__).resolve().parent
    candidates += [here / "dummy_gcxgc", here.parent / "dummy_gcxgc",
                   here.parent.parent / "dummy_gcxgc"]
    for c in candidates:
        if (c / "data").is_dir() and (c / "target_table.csv").exists():
            return c
    return None


@pytest.fixture(scope="session")
def dummy_dir() -> Path:
    d = _dummy_dir()
    if d is None:
        pytest.skip("dummy data pack not found; set PYGCGC_DUMMY_DATA")
    return d


@pytest.fixture(scope="session")
def dummy_runs(dummy_dir):
    from pygcgc import read_many
    files = sorted((dummy_dir / "data").glob("*.cdf"))
    return {r.name: r for r in read_many(files)}


@pytest.fixture(scope="session")
def dummy_anchors(dummy_dir):
    from pygcgc import load_anchors
    return load_anchors(dummy_dir / "target_table.csv")


@pytest.fixture(scope="session")
def dummy_truth(dummy_dir):
    rows = {}
    with open(dummy_dir / "truth_per_file.csv", newline="") as fh:
        for r in csv.DictReader(fh):
            rows[(r["file"].replace(".cdf", ""), r["compound"])] = r
    return rows


# ------------------------------------------------------------ synthetic folder
GUI_PEAKS = [(40.0, 1.10, 50, 2e5), (90.0, 1.90, 51, 3e5), (150.0, 2.70, 52, 2e5),
             (200.0, 1.40, 53, 2e5), (70.0, 2.20, 54, 2e5)]
GUI_PHASE = 0.7


def write_run_andi(run, path):
    """Write a Run as a minimal ANDI file."""
    from netCDF4 import Dataset
    coo = run.data.tocoo()
    order = np.lexsort((coo.col, coo.row))
    rows, cols, vals = coo.row[order], coo.col[order], coo.data[order]
    counts = np.bincount(rows, minlength=run.n_scans).astype(np.int32)
    with Dataset(path, "w", format="NETCDF3_CLASSIC") as ds:
        ds.global_mass_min = float(run.mz[0])
        ds.global_mass_max = float(run.mz[-1])
        ds.administrative_comments = "synthetic"
        ds.createDimension("scan_number", run.n_scans)
        ds.createDimension("point_number", vals.size)
        v = ds.createVariable("scan_acquisition_time", "f8", ("scan_number",))
        v.units = "seconds"
        v[:] = run.times
        ds.createVariable("point_count", "i4", ("scan_number",))[:] = counts
        ds.createVariable("mass_values", "f4", ("point_number",))[:] = run.mz[cols].astype(np.float32)
        ds.createVariable("intensity_values", "f4", ("point_number",))[:] = vals.astype(np.float32)


@pytest.fixture
def synthetic_folder(tmp_path):
    """Folder in the dummy-pack layout: three runs (one without anchors),
    sequence.csv, target_table.csv and method.json."""
    folder = tmp_path / "synthetic_pack"
    (folder / "data").mkdir(parents=True)
    runs = {
        "RUN_A": make_run(GUI_PEAKS, phase=GUI_PHASE, rate=50.0, seed=1),
        "RUN_B": make_run([(t1 + 1.5, t2, ion, a * 0.7) for t1, t2, ion, a in GUI_PEAKS],
                          phase=GUI_PHASE, rate=50.0, seed=2),
        "RUN_C_BLANK": make_run([(120.0, 0.9, 50, 3e3)], phase=GUI_PHASE, rate=50.0, seed=3),
    }
    # A rising "bleed" baseline with noise on m/z 54, so 🧹 has work to do.
    for k, (name, run) in enumerate(runs.items()):
        rng = np.random.default_rng(10 + k)
        dense = run.data.toarray()
        frac = (run.times - run.times[0]) / (run.times[-1] - run.times[0])
        dense[:, 54 - run.mz[0]] += 40 + 60 * frac + rng.normal(0, 3, run.n_scans)
        run = Run(times=run.times, mz=run.mz, data=sp.csr_matrix(dense), name=name)
        write_run_andi(run, folder / "data" / f"{name}.cdf")
    with open(folder / "sequence.csv", "w", newline="") as fh:
        fh.write("file,type,level,level_ng\nRUN_A.cdf,calibration,L1,1\n"
                 "RUN_B.cdf,sample,,\nRUN_C_BLANK.cdf,tube blank,,0\n")
    with open(folder / "target_table.csv", "w", newline="") as fh:
        fh.write("compound,role,quant_ion,t1_expected_s,t2_expected_s,t1_window_s\n")
        for i, (t1, t2, ion, _) in enumerate(GUI_PEAKS):
            fh.write(f"std{i},istd,{ion},{t1},{t2},6\n")
    (folder / "method.json").write_text('{"modulation_period_s": 4.0}')
    return folder


# ------------------------------------------------------------- alkane folder
ALKANE_CARBONS = list(range(10, 26))
ALKANE_T1_START, ALKANE_STEP = 200.0, 40.0


def _alkane_peaks(shift_t1=0.0, shift_t2=0.0, branched=True):
    peaks = []
    for i in range(len(ALKANE_CARBONS)):
        t1 = ALKANE_T1_START + i * ALKANE_STEP + shift_t1
        t2 = 1.0 + 0.01 * i + shift_t2
        peaks += [(t1, t2, 57, 2e5), (t1, t2, 71, 1.32e5), (t1, t2, 85, 8.8e4)]   # the measured n-alkane ratios
    if branched:
        for i in range(len(ALKANE_CARBONS)):
            t1 = ALKANE_T1_START + i * ALKANE_STEP + 5.0 + shift_t1
            t2 = 1.3 + 0.01 * i + shift_t2
            peaks += [(t1, t2, 43, 1.5e5), (t1, t2, 57, 3e4)]
    return peaks


@pytest.fixture
def alkane_folder(tmp_path):
    """Folder with an n-alkane standard ("STD", C10-C25 on m/z 57, plus a
    parallel branched-alkane series on m/z 43) and a shifted "SAMPLE" run,
    for the ⛓️ hydrocarbon mapper. No target table: phase comes from the
    quiet-band (void) method."""
    folder = tmp_path / "alkane_pack"
    (folder / "data").mkdir(parents=True)
    runs = {
        "STD": make_run(_alkane_peaks(), phase=GUI_PHASE, rate=50.0, mz=(40, 90),
                        duration=1100.0, seed=1),
        "SAMPLE": make_run(_alkane_peaks(shift_t1=6.0, shift_t2=0.03), phase=GUI_PHASE, rate=50.0,
                           mz=(40, 90), duration=1100.0, seed=2),
    }
    for name, run in runs.items():
        write_run_andi(run, folder / "data" / f"{name}.cdf")
    with open(folder / "sequence.csv", "w", newline="") as fh:
        fh.write("file,type\nSTD.cdf,calibration\nSAMPLE.cdf,sample\n")
    (folder / "method.json").write_text('{"modulation_period_s": 4.0}')
    return folder


# ------------------------------------------------------- calibration folder
CAL_LEVELS = {"L1": 0.5, "L2": 1.0, "L3": 2.0, "L4": 5.0, "L5": 10.0}
CAL_TARGETS = [
    # name, role, ion, qual ion, qual ratio, t1, t2, response per unit amount
    ("IS", "istd", 50, 51, 0.5, 60.0, 1.20, None),
    ("A", "target", 52, 53, 0.3, 110.0, 2.00, 1.0),
    ("B", "target", 54, None, None, 170.0, 2.60, 0.6),
    ("C", "target", 52, 53, 0.3, 200.0, 1.50, 1.4),
]
SAMPLE_AMOUNTS = {"S1": {"A": 3.0, "B": 0.8, "C": 7.0}, "S2": {"A": 0.7, "B": 6.0, "C": 1.5}}
IS_AREA = 3e5
SCALE = 1e5


def make_cal_run(amounts: dict, shift=(0.0, 0.0), seed=0, noise=0.0):
    peaks = []
    for name, role, ion, qion, qr, t1, t2, rf in CAL_TARGETS:
        area = IS_AREA if role == "istd" else amounts.get(name, 0.0) * rf * SCALE
        if area <= 0:
            continue
        peaks.append((t1 + shift[0], t2 + shift[1], ion, area))
        if qion:
            peaks.append((t1 + shift[0], t2 + shift[1], qion, area * qr))
    return make_run(peaks, phase=GUI_PHASE, rate=50.0, seed=seed, noise=noise)


@pytest.fixture
def cal_folder(tmp_path):
    folder = tmp_path / "cal_pack"
    (folder / "data").mkdir(parents=True)
    rows = []
    for k, (lv, amt) in enumerate(CAL_LEVELS.items()):
        run = make_cal_run({n: amt for n in ("A", "B", "C")}, seed=k)
        write_run_andi(run, folder / "data" / f"CAL_{lv}.cdf")
        rows.append(f"CAL_{lv}.cdf,calibration,{lv},,")
    for k, (name, amts) in enumerate(SAMPLE_AMOUNTS.items()):
        run = make_cal_run(amts, shift=(1.5, 0.02), seed=10 + k)
        write_run_andi(run, folder / "data" / f"{name}.cdf")
        rows.append(f"{name}.cdf,sample,,,2.0")
    (folder / "sequence.csv").write_text("file,type,level,cal_set,volume_m3\n" + "\n".join(rows) + "\n")
    lines = ["compound,role,quant_ion,qual_ions,qual_ratios,t1_expected_s,t1_window_s,t2_expected_s,t2_window_s,istd"]
    for name, role, ion, qion, qr, t1, t2, rf in CAL_TARGETS:
        lines.append(f"{name},{role},{ion},{qion or ''},{qr or ''},{t1},10,{t2},0.2,{'' if role == 'istd' else 'IS'}")
    (folder / "target_table.csv").write_text("\n".join(lines) + "\n")
    lv = ["level,A,B,C"] + [f"{k},{v},{v},{v}" for k, v in CAL_LEVELS.items()]
    (folder / "calibration_levels.csv").write_text("\n".join(lv) + "\n")
    (folder / "method.json").write_text('{"modulation_period_s": 4.0}')
    return folder


@pytest.fixture
def drift_folder(tmp_path):
    """Two calibration blocks ('mon', 'fri') with a target-specific
    sensitivity drift between them (an internal standard alone would not
    show this: it is the target/IS response ratio that has drifted), and
    one sample sitting between them in run order."""
    folder = tmp_path / "drift_pack"
    (folder / "data").mkdir(parents=True)
    rows = []

    def make_scaled(amounts, rf_scale, shift=(0.0, 0.0), seed=0):
        peaks = []
        for name, role, ion, qion, qr, t1, t2, rf in CAL_TARGETS:
            area = IS_AREA if role == "istd" else amounts.get(name, 0.0) * rf * rf_scale * SCALE
            if area <= 0:
                continue
            peaks.append((t1 + shift[0], t2 + shift[1], ion, area))
            if qion:
                peaks.append((t1 + shift[0], t2 + shift[1], qion, area * qr))
        return make_run(peaks, phase=GUI_PHASE, rate=50.0, seed=seed)

    for k, (lv, amt) in enumerate(CAL_LEVELS.items()):
        run = make_scaled({n: amt for n in ("A", "B", "C")}, rf_scale=1.0, seed=k)
        write_run_andi(run, folder / "data" / f"MON_{lv}.cdf")
        rows.append(f"MON_{lv}.cdf,calibration,{lv},mon,")

    true_amount = 4.0
    # Chosen so that blending the two bracketing curves (equal weight, since
    # the sample sits exactly halfway between the two blocks by run order)
    # recovers true_amount exactly: with sensitivities 1.0 and 0.5, that
    # blend is 1.5 * (true response-space amount), so the injected
    # sensitivity is scaled down by 1/1.5 to compensate.
    rf_scale = 1.0 / 1.5
    sample = make_scaled({"A": true_amount, "B": 1.0, "C": 1.0}, rf_scale=rf_scale,
                         shift=(1.5, 0.02), seed=50)
    write_run_andi(sample, folder / "data" / "S_MID.cdf")
    rows.append("S_MID.cdf,sample,,,2.0")

    for k, (lv, amt) in enumerate(CAL_LEVELS.items()):
        run = make_scaled({n: amt for n in ("A", "B", "C")}, rf_scale=0.5, seed=20 + k)
        write_run_andi(run, folder / "data" / f"FRI_{lv}.cdf")
        rows.append(f"FRI_{lv}.cdf,calibration,{lv},fri,")

    (folder / "sequence.csv").write_text("file,type,level,cal_set,volume_m3\n" + "\n".join(rows) + "\n")
    lines = ["compound,role,quant_ion,qual_ions,qual_ratios,t1_expected_s,t1_window_s,t2_expected_s,t2_window_s,istd"]
    for name, role, ion, qion, qr, t1, t2, rf in CAL_TARGETS:
        lines.append(f"{name},{role},{ion},{qion or ''},{qr or ''},{t1},10,{t2},0.2,{'' if role == 'istd' else 'IS'}")
    (folder / "target_table.csv").write_text("\n".join(lines) + "\n")
    lv2 = ["level,A,B,C"] + [f"{k},{v},{v},{v}" for k, v in CAL_LEVELS.items()]
    (folder / "calibration_levels.csv").write_text("\n".join(lv2) + "\n")
    (folder / "method.json").write_text('{"modulation_period_s": 4.0}')
    return folder, true_amount
