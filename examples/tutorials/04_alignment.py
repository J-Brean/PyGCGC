"""Milestone 5 walkthrough: align two runs on their internal standards and
compare them.

    python examples/tutorials/04_alignment.py path/to/demo_gcxgc
"""

import sys
from pathlib import Path

import pygcgc as g

folder = Path(sys.argv[1] if len(sys.argv) > 1 else "demo_gcxgc")
period = 8.0

anchors = g.load_anchors(folder / "target_table.csv")
ref = g.read_run(folder / "data" / "CAL_L3_20ng.cdf")
run = g.read_run(folder / "data" / "SAMPLE_1_urban_air.cdf")

ref_fold = g.fold(ref, period, g.auto_phase(ref, period, anchors).phase)
run_fold = g.fold(run, period, g.auto_phase(run, period, anchors).phase)

# Find every internal standard in both runs, and pair them up as anchors.
ref_found = g.propagate_anchors(anchors, ref, ref_fold)
run_found = g.propagate_anchors(anchors, run, run_fold)
points = [g.AnchorPoint(name, ref_found[name].run_t1, ref_found[name].run_t2,
                        run_found[name].run_t1, run_found[name].run_t2)
          for name in ref_found if ref_found[name] and run_found.get(name)]
print(f"{len(points)} of {len(anchors)} internal standards found in both runs")

transform = g.fit_alignment(points)   # kind="auto"
print(transform.summary())
for r in transform.residuals():
    print(f"  {r['name']:<20s} resid t1 {r['resid_t1']:+.2f} s   resid t2 {r['resid_t2'] * 1e3:+.0f} ms")

# Resample the urban sample's TIC image onto CAL_L3_20ng's grid, so they can be compared
# pixel for pixel (this is what the ⚖️ Compare window's blink/difference/overlay
# modes are built on).
run_tic = run_fold.tic_image(run)
aligned = g.resample_to_reference(run_tic, run_fold, ref_fold, transform)
ref_tic = ref_fold.tic_image(ref)
print(f"reference image {ref_tic.shape}, aligned image {aligned.shape}")

# A robust brightness-matched difference, the same idea the ⚖️ Compare window's
# Difference mode uses (see gui/views/compare.robust_scale): match on the
# median ratio of the two images' brightest shared pixels, so a run that
# is simply more concentrated doesn't look like a wall of difference.
import numpy as np

mask = (ref_tic > np.percentile(ref_tic, 90)) & (aligned > np.percentile(aligned, 90)) & (aligned > 0)
k = float(np.median(ref_tic[mask] / aligned[mask])) if mask.any() else 1.0
diff = ref_tic - k * aligned
print(f"difference image built (scale {k:.3g}); "
     f"{(abs(diff) > 3 * diff.std()).mean():.1%} of pixels differ by more than 3 SD")
