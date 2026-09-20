"""Milestone 6 walkthrough: find every blob in a run, look at what makes
each one interesting, and try the hand-editing operations.

    python examples/milestone6_blobs.py path/to/dummy_gcxgc
"""

import sys
from pathlib import Path

import pygcgc as g
from pygcgc.core.peaks import add_blob_at, delete_blobs, merge_blobs

folder = Path(sys.argv[1] if len(sys.argv) > 1 else "dummy_gcxgc")
period = 4.0

anchors = g.load_anchors(folder / "target_table.csv")
run = g.read_run(folder / "data" / "CAL_L4.cdf")
phase = g.auto_phase(run, period, anchors).phase
fold = g.fold(run, period, phase)

clean = g.remove_background(run, fold).run
blobs = g.find_blobs(clean, fold)
print(f"{len(blobs)} blobs found")

# The five biggest, by volume, with what makes each one worth a second look.
order = blobs.volume.argsort()[::-1][:5]
for i in order:
    flag = " CO-ELUTION?" if blobs.co_elution[i] else ""
    print(f"  t1 {blobs.t1[i] / 60:6.2f} min   t2 {blobs.t2[i]:5.2f} s   "
         f"volume {blobs.volume[i]:10,.0f}   S/N {blobs.snr[i]:6.1f}   "
         f"width {blobs.width_t1[i]:4.1f} s x {blobs.width_t2[i] * 1e3:4.0f} ms{flag}")

# Delete is just bookkeeping (zero the footprint, renumber); merge is an
# exact combination of two blobs' volume, spectrum, position and width, with
# no need to re-read the raw data (see the module docstring in core/peaks.py).
i, j = int(order[0]), int(order[1])
before = len(blobs)
trimmed = delete_blobs(blobs, [i])
print(f"deleted the biggest blob: {before} -> {len(trimmed)}")

combined = merge_blobs(blobs, [i, j])
print(f"merged the two biggest into one: volume {combined.volume[-1]:,.0f} "
     f"(exactly {blobs.volume[i] + blobs.volume[j]:,.0f})")

# Add grows a fresh footprint from a clicked pixel, for a real peak that
# automatic detection missed. Here: put one back after deleting it.
recovered = add_blob_at(trimmed, clean, fold, float(blobs.t1[i]), float(blobs.t2[i]))
print(f"added a blob back at the same spot: volume {recovered.volume[-1]:,.0f} "
     f"(originally {blobs.volume[i]:,.0f})")
