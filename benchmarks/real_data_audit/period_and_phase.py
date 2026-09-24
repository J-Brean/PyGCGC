"""Are the modulation period and phase found without being told, on the real runs?

The period is not in the LSC file. `estimate_period` finds it from the total signal's
periodicity alone; the anchor-free phase is a place to cut (the start of the
modulator-locked bleed), not the modulator's proven start. The user set 8.0 s and
2.77 s by hand in the GUI.

Run from the repository root: python benchmarks/real_data_audit/period_and_phase.py
"""

import time

from common import PERIOD, PHASE, RUNS, load

from pygcgc.core.period import estimate_period
from pygcgc.core.phase import phase_unanchored

for key in RUNS:
    run = load(key)
    t = time.perf_counter()
    est = estimate_period(run)
    print(f"{key:5}: {est.summary()} ({time.perf_counter() - t:.0f} s)")
    print(f"{'':7}anchor-free phase at {PERIOD} s: {phase_unanchored(run, PERIOD).summary()}; set by hand: {PHASE} s")
