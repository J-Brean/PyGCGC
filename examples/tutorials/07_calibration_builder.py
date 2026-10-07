"""📈 Calibration builder walkthrough: curves from the four standards and a compound
table, normalised to the internal standards, edited, then read off for the samples.
This is what the 📈 Calibrate window does, step by step, without the window.

    python examples/tutorials/07_calibration_builder.py path/to/demo_gcxgc
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np

import pygcgc as g
from pygcgc.core.calbuilder import (CalSession, compound_templates, curve_rows, load_compounds,
                                    quantify_samples, sample_rows, t2_offsets)
from pygcgc.core.quant import write_dicts

folder = Path(sys.argv[1] if len(sys.argv) > 1 else "demo_gcxgc")
period = json.loads((folder / "method.json").read_text())["modulation_period_s"]
rows = list(csv.DictReader(open(folder / "sequence.csv", newline="")))
standards = [Path(r["file"]).stem for r in rows if r["type"] == "calibration"]
samples = [Path(r["file"]).stem for r in rows if r["type"] == "sample"]
volumes = {Path(r["file"]).stem: float(r["volume_m3"]) for r in rows if r["volume_m3"]}

# 1. Every run folded on its internal standards, cleaned, and its blobs found.
runs = g.read_many([folder / "data" / r["file"] for r in rows])
phases = g.phase_for_sequence(runs, period, g.load_anchors(folder / "target_table.csv"))
blobs = {}
for run, est in zip(runs, phases):
    f = g.fold(run, period, est.phase)
    blobs[run.name] = g.find_blobs(g.remove_background(run, f).run, f)
    print(f"{run.name:<26} phase {est.phase:.3f} s   {len(blobs[run.name]):5d} blobs")

# 2. The compounds and their amounts per level come from a table. Each is matched to
#    the blob at its position in every standard and divided by its nearest internal standard.
session = CalSession(runs=standards, compounds=load_compounds(folder / "calibration_compounds.csv", standards),
                     unit="ng")
points = session.points(blobs)
fits = session.fits(points)
ok = [n for n, fit in fits.items() if isinstance(fit, g.CalibrationFit) and fit.ok]
print(f"\n{len(ok)} of {len(fits)} curves pass the back-calculation check")
for name in ("naphthalene", "fluorene", "pyrene", "n-C20"):
    fit = fits[name]
    pts = ", ".join(f"{points[(name, r)].x:g} ng -> {points[(name, r)].ratio:.3f}" for r in standards)
    print(f"  {name:<12} {fit.model} {fit.weights}, R2 {fit.r2:.4f}   ({pts})")

# 3. A point can be left out of one curve, with a reason that is kept. Here the top
#    standard of naphthalene, where the detector starts to bend over.
session.exclude("naphthalene", "CAL_L4_50ng", "detector bending over at 50 ng")
fits = session.fits(points)
print(f"\nnaphthalene without the 50 ng point: {fits['naphthalene'].model}, "
      f"worst back-calculation {np.max(np.abs(fits['naphthalene'].rel_error)):.1f}%")

# 4. The samples: each compound is found by its position and spectrum in the standards
#    (moved by the run's own t2 offset from them), and read off its curve.
offsets = {run: dt for run, (dt, _n) in t2_offsets(blobs, "CAL_L3_20ng").items()}
templates = compound_templates(session, blobs)
results = quantify_samples(session, fits, templates, {s: blobs[s] for s in samples}, offsets)

truth = {(r["run"], r["compound"]): float(r["amount_ng"])
         for r in csv.DictReader(open(folder / "truth" / "truth_amounts.csv", newline=""))}
print("\nrun                        compound             found ng   true ng   flags")
for r in results:
    if r.compound in ("naphthalene", "phenanthrene", "pyrene", "benzo[a]pyrene", "n-C20", "9-fluorenone"):
        true = truth.get((r.run, r.compound), np.nan)
        print(f"{r.run:<26} {r.compound:<18} {r.amount:9.2f} {true:9.2f}   {'; '.join(r.flags)}")

# 5. The same tables the window exports.
write_dicts(curve_rows(session, fits), "07_curves.csv")
write_dicts(sample_rows(results, session.unit, volumes), "07_samples.csv")
print("\nwrote 07_curves.csv and 07_samples.csv")
