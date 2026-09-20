"""Milestone 7 walkthrough: flag which blobs look like siloxane, and fit
and subtract the siloxane component from the whole run.

    python examples/milestone7_siloxane.py path/to/dummy_gcxgc
"""

import sys
from pathlib import Path

import pygcgc as g

folder = Path(sys.argv[1] if len(sys.argv) > 1 else "dummy_gcxgc")
period = 4.0

anchors = g.load_anchors(folder / "target_table.csv")
run = g.read_run(folder / "data" / "CAL_L4.cdf")
phase = g.auto_phase(run, period, anchors).phase
fold = g.fold(run, period, phase)

bg = g.remove_background(run, fold)
blobs = g.find_blobs(bg.run, fold)
scores = g.siloxane_scores(blobs)
flagged = g.flag_siloxane_blobs(blobs)
print(f"{flagged.sum()} of {len(blobs)} blobs flagged as probable siloxane "
     f"(threshold {g.SiloxaneSettings().score_threshold})")

order = scores.argsort()[::-1][:8]
for i in order:
    print(f"  t1 {blobs.t1[i] / 60:6.2f} min   t2 {blobs.t2[i]:5.2f} s   "
         f"score {scores[i]:.2f}{'  <-- flagged' if flagged[i] else ''}")

# Fit and subtract the siloxane component everywhere, and check a real
# analyte on a shared ion (behenic acid, m/z 73) survives.
result = g.remove_siloxanes(bg.run, fold)
print(result.summary())


def box_sum(run, f, apex_scan_time_s, ion, before=0.4, after=0.6, cb=5, ca=6):
    col, t2 = f.to_fold(apex_scan_time_s)
    img, _ = f.extend(f.ion_image(run, ion), 1.99)
    r0 = int(round(t2 / f.row_pitch))
    return img[max(r0 - int(before / f.row_pitch), 0):r0 + int(after / f.row_pitch),
               col - cb:col + ca].sum()


behenic_t = None
import csv

with open(folder / "truth_per_file.csv", newline="") as fh:
    for row in csv.DictReader(fh):
        if row["file"] == "CAL_L4.cdf" and row["compound"] == "behenic acid":
            behenic_t = float(row["apex_scan_time_s"])
if behenic_t is not None:
    before = box_sum(bg.run, fold, behenic_t, 73)
    after = box_sum(result.run, fold, behenic_t, 73)
    print(f"behenic acid on m/z 73: {before:,.0f} -> {after:,.0f} "
         f"({100 * after / before:.1f}% remains, milestone bar is >= 95%)")
