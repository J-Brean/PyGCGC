"""Milestone 4 walkthrough: clean, calibrate and quantify a whole sequence,
then print a few results and write the four tables.

    python examples/milestone4_quant.py path/to/dummy_gcxgc
"""

import csv
import sys
from pathlib import Path

import pygcgc as g

folder = Path(sys.argv[1] if len(sys.argv) > 1 else "dummy_gcxgc")
period = 4.0

# The sequence file sets the order; the target table gives the internal standards.
rows = list(csv.DictReader(open(folder / "sequence.csv", newline="")))
runs = g.read_many([folder / "data" / r["file"] for r in rows])
targets = g.load_targets(folder / "target_table.csv")
anchors = g.load_anchors(folder / "target_table.csv")
levels = g.load_levels(folder / "calibration_levels_ng.csv")
extras = g.read_sequence_extras(folder / "sequence.csv")
phases = g.phase_for_sequence(runs, period, anchors)

# Fold and clean every run, and describe it for quantification.
prepared, info = [], []
for run, est, row in zip(runs, phases, rows):
    f = g.fold(run, period, est.phase)
    bg = g.remove_background(run, f)
    prepared.append((run, bg.run, f, bg.noise))
    ex = extras.get(run.name, {})
    info.append(g.SequenceRun(run.name, row["type"], row["level"], ex.get("cal_set", ""), ex.get("volume")))
    print(f"cleaned {run.name}")

report = g.quantify_prepared(prepared, info, targets, levels)

fits = [d["quant"] for d in report.fits.values() if isinstance(d.get("quant"), g.CalibrationFit)]
print(f"{sum(f.ok for f in fits)} of {len(fits)} calibration curves pass")
print(report.fit("fluorene").summary())

for r in report.results:
    if r.run == "SAMPLE_A" and r.compound in ("naphthalene", "fluorene", "pyrene", "anthracene"):
        print(f"{r.compound:<12} {r.amount:8.3f} ± {r.amount_se:.3f} ng   {r.status}  {', '.join(r.flags)}")

signals = [s for s in report.signals if s.run == "SAMPLE_A"]
print(f"SAMPLE_A: {len(signals)} peaks, each reported as volume / nearest internal standard")
for s in sorted(signals, key=lambda s: -s.volume)[:5]:
    print(f"  t1 {s.t1 / 60:6.2f} min  t2 {s.t2:4.2f} s  ratio {s.ratio:7.3f} to {s.istd}"
          + (f"  ~{s.amount:.2f} ng as {s.surrogate}" if s.surrogate else ""))

for p in report.write("milestone4_results"):
    print("wrote", p)
