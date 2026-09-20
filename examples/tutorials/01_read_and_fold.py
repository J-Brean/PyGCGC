"""Milestone 1 walkthrough: read a sequence, check timing, find the phase,
check the period, fold, and save a figure.

Run from the repository folder:
    python examples/tutorials/01_read_and_fold.py path/to/dummy_gcxgc
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

import pygcgc as g

folder = Path(sys.argv[1] if len(sys.argv) > 1 else "dummy_gcxgc")
period = 4.0  # from the instrument method

# 1. Read every run in the sequence (threaded).
files = sorted((folder / "data").glob("*.cdf"))
runs = g.read_many(files)
print(f"Read {len(runs)} runs")

# 2. Check scan timing for each run.
for run in runs:
    report = g.check_timing(run.times)
    status = "OK" if report.ok else "CHECK"
    print(f"{run.name:<12} {report.rate:.4f} Hz  gaps {report.missing_scans}  {status}")

# 3. Phase for the whole sequence, using the internal standards as anchors.
#    Runs without them (tube blanks) borrow the phase of the others.
anchors = g.load_anchors(folder / "target_table.csv")
phases = g.phase_for_sequence(runs, period, anchors)
for run, est in zip(runs, phases):
    print(f"{run.name:<12} {est.summary()}")

# 4. Check the modulation period on one calibration run.
cal = next(r for r in runs if r.name == "CAL_L4")
print(g.period_from_anchors(cal, anchors, period).summary())

# 5. Fold one sample with and without the phase correction, and plot both.
sample = next(r for r in runs if r.name == "SAMPLE_B")
phase = phases[runs.index(sample)].phase
naive = g.fold(sample, period, 0.0)
fixed = g.fold(sample, period, phase)
img_naive = naive.tic_image(sample)
img_fixed, t2_ext = fixed.extend(fixed.tic_image(sample), 1.15)

fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.0), dpi=300, sharex=True)
for ax, img, f, t2, title in [
    (axes[0], img_naive, naive, naive.t2, "No phase correction"),
    (axes[1], img_fixed, fixed, t2_ext, f"Phase {phase:.3f} s, t2 axis extended"),
]:
    floor = np.percentile(img[img > 0], 5)
    ax.imshow(np.clip(img, floor, None), origin="lower", aspect="auto", cmap="magma",
              norm=LogNorm(vmin=floor, vmax=img.max()),
              extent=[f.t1[0] / 60, (f.t1[-1] + period) / 60, 0, t2[-1] + f.row_pitch])
    ax.set_ylabel("$t_2$ (s)", fontsize=8)
    ax.set_title(title, fontsize=8)
    ax.tick_params(labelsize=7)
axes[1].axhspan(period, t2_ext[-1] + fixed.row_pitch, color="w", alpha=0.12, lw=0)
axes[1].set_xlabel("$t_1$ (min)", fontsize=8)
fig.tight_layout()
fig.savefig("01_sample_B.png")
print("Saved 01_sample_B.png")
