"""Milestone 3 walkthrough: clean one sample, check the result, and plot
raw against cleaned for a bleed ion.

    python examples/tutorials/02_background.py path/to/dummy_gcxgc
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import pygcgc as g

folder = Path(sys.argv[1] if len(sys.argv) > 1 else "dummy_gcxgc")
period = 4.0
ion = 73          # a column-bleed ion that behenic acid is also quantified on

# 1. Read the sequence and find every run's phase from the internal standards.
runs = g.read_many(sorted((folder / "data").glob("*.cdf")))
anchors = g.load_anchors(folder / "target_table.csv")
phases = {r.name: e.phase for r, e in zip(runs, g.phase_for_sequence(runs, period, anchors))}

# 2. Clean the tube blank, then a calibration run with the blank subtracted.
blank = next(r for r in runs if r.name == "TUBE_BLANK")
blank_fold = g.fold(blank, period, phases[blank.name])
blank_clean = g.remove_background(blank, blank_fold)

run = next(r for r in runs if r.name == "CAL_L1")
f = g.fold(run, period, phases[run.name])
settings = g.BackgroundSettings(subtract_blank=True)
res = g.remove_background(run, f, settings, blank=blank_clean, blank_fold=blank_fold,
                          blank_name=blank.name)
print(res.summary())

# 3. Check the leftover background on the bleed ion: it should be centred on
#    zero with a spread matching the noise map.
residual, sigma = res.residual_pixels(f, ion)
print(f"m/z {ion}: residual mean {residual.mean():+.2f}, SD {residual.std():.2f}, noise map {sigma:.2f}")

# 4. Plot raw and cleaned images of the bleed ion side by side.
raw_img = f.ion_image(run, ion)
clean_img = f.ion_image(res.run, ion)
extent = [f.t1[0] / 60, (f.t1[-1] + period) / 60, 0, period]
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.8), dpi=300, sharey=True)
for ax, img, title in [(axes[0], raw_img, "raw"), (axes[1], clean_img, "cleaned")]:
    top = np.percentile(clean_img, 99.9)
    ax.imshow(img, origin="lower", aspect="auto", cmap="magma", vmin=0, vmax=top, extent=extent)
    ax.set_title(f"m/z {ion}, {title}", fontsize=8)
    ax.set_xlabel("$t_1$ (min)", fontsize=8)
    ax.tick_params(labelsize=7)
axes[0].set_ylabel("$t_2$ (s)", fontsize=8)
fig.tight_layout()
fig.savefig("02_background.png")
print("Saved 02_background.png")
