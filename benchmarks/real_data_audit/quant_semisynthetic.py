"""Does 🧪 quantification recover known amounts in a real matrix?

No real calibration sequence exists yet, so one is made: two internal standards and
four targets are written into the real runs (cropped to 12-52 min) at scan level,
with counting noise, at five calibration levels, in three samples of known amount
and in a blank. Two targets sit on ions the fuel and exhaust matrix hardly carries
(m/z 136, 212); two sit on busy ions (m/z 57 near the alkanes, m/z 105 among the
alkylbenzenes). Everything planted in a run is multiplied by a random injection
factor (0.8 to 1.2), which the internal standard should cancel, and moved by a
small retention drift. The target table gives rough positions and leaves the
qualifier ratios blank, as a first table typed from a quick look would.

The folder is written in the layout the CLI and GUI expect (data/*.cdf,
sequence.csv, target_table.csv, calibration_levels.csv, method.json) under
real_data/quant_validation/, then `pygcgc quant` runs on it and the amounts are
compared with the truth. Open the folder in the GUI to review it by hand.

Run from the repository root: python benchmarks/real_data_audit/quant_semisynthetic.py
"""

import csv
import json
import shutil

import numpy as np
from common import PERIOD, PHASE, REAL, load

from pygcgc import crop_run, fold, write_andi
from pygcgc.cli import main as cli
from pygcgc.core.synthetic import TruePeak, inject_peaks

rng = np.random.default_rng(20260924)
out = REAL / "quant_validation"
if out.exists():
    shutil.rmtree(out)
(out / "data").mkdir(parents=True)

T_START, T_END = 12 * 60.0, 52 * 60.0
# name, role, quantifier, qualifiers, spectrum {m/z: relative}, t1 (s), t2 (s), response per ng
COMPOUNDS = [
    ("phenanthrene-d10", "istd", 188, (189, 160), {188: 100, 189: 15, 160: 8, 94: 8, 80: 6}, 25.2 * 60, 2.45, None),
    ("chrysene-d12", "istd", 240, (236, 241), {240: 100, 236: 30, 241: 20, 120: 12}, 44.1 * 60, 3.85, None),
    ("T_unique136", "target", 136, (108, 137), {136: 100, 108: 12, 137: 11, 68: 6}, 18.3 * 60, 2.05, 1.0),
    ("T_unique212", "target", 212, (213, 106), {212: 100, 213: 17, 106: 10, 104: 5}, 36.6 * 60, 3.40, 0.7),
    ("T_branched57", "target", 57, (71, 85), {57: 100, 43: 80, 71: 45, 85: 25, 113: 8, 183: 4}, 30.9 * 60, 0.93, 1.3),
    ("T_arom105", "target", 105, (120, 77), {105: 100, 120: 30, 77: 10, 91: 8, 79: 6}, 20.4 * 60, 1.85, 0.9),
]
IS_AREA, SCALE = 1.5e6, 1.5e5          # internal standard area; area per ng per unit response (all ions summed)
LEVELS = {"L1": 1.0, "L2": 2.0, "L3": 5.0, "L4": 10.0, "L5": 20.0}
SAMPLES = {"S1": {"T_unique136": 3.0, "T_unique212": 7.5, "T_branched57": 0.6, "T_arom105": 14.0},
           "S2": {"T_unique136": 12.0, "T_unique212": 1.5, "T_branched57": 4.0, "T_arom105": 2.5},
           "S3": {"T_unique136": 0.4, "T_unique212": 16.0, "T_branched57": 9.0, "T_arom105": 6.0}}
targets = [c[0] for c in COMPOUNDS if c[1] == "target"]
SEQUENCE = ([(f"CAL_{lv}", "calibration", lv, dict.fromkeys(targets, amt), bg)
             for (lv, amt), bg in zip(LEVELS.items(), ["ATWC", "BTWC", "T191", "T60", "ATWC"])]
            + [(s, "sample", "", amts, bg) for (s, amts), bg in zip(SAMPLES.items(), ["BTWC", "T191", "T60"])]
            + [("BLANK", "tube blank", "", {}, "BTWC")])

truth = []
for k, (name, kind, level, amounts, bg) in enumerate(SEQUENCE):
    full = load(bg)
    run = crop_run(full, T_START, T_END)
    # Keep the modulation boundaries of the whole acquisition.
    f = fold(run, PERIOD, (float(full.times[0]) + PHASE - float(run.times[0])) % PERIOD)
    injection = float(rng.uniform(0.8, 1.2))
    peaks = []
    for comp, role, _q, _quals, spec, t1, t2, response in COMPOUNDS:
        area = IS_AREA if role == "istd" else amounts.get(comp, 0.0) * response * SCALE
        if area <= 0:
            continue
        vec = np.zeros(run.n_mz)
        for m, v in spec.items():
            vec[m - run.mz[0]] = v
        t1_run, t2_run = t1 + rng.normal(0, 2.0), t2 + rng.normal(0, 0.01)
        # TruePeak works in fold coordinates; a column's t1 is its middle, as in the target table
        peaks.append(TruePeak(row=t2_run / f.row_pitch, col=(t1_run - f.t1[0]) / PERIOD - 0.5, area=area * injection,
                              spectrum=vec, sigma_row=0.05 / f.row_pitch, sigma_col=0.9, tau_row=0.03 / f.row_pitch))
        if role == "target":
            truth.append((name, comp, amounts.get(comp, 0.0)))
    planted = inject_peaks(run, f, peaks, seed=100 + k)
    planted.name = name
    write_andi(planted, out / "data" / f"{name}.cdf", comment="semi-synthetic quantification check")
    print(f"wrote {name} on {bg} (injection factor {injection:.2f})")

with open(out / "sequence.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["file", "type", "level"])
    w.writerows([f"{name}.cdf", kind, level] for name, kind, level, _a, _b in SEQUENCE)
with open(out / "target_table.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["compound", "role", "quant_ion", "qual_ions", "qual_ratios", "t1_expected_s", "t1_window_s",
                "t2_expected_s", "t2_window_s", "istd"])
    for comp, role, q, quals, _spec, t1, t2, _r in COMPOUNDS:
        w.writerow([comp, role, q, ";".join(map(str, quals)), "", round(t1 + rng.normal(0, 3), 1), 12,
                    round(t2 + rng.normal(0, 0.03), 2), 0.25, ""])
with open(out / "calibration_levels.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["level"] + targets)
    w.writerows([lv] + [amt] * len(targets) for lv, amt in LEVELS.items())
(out / "method.json").write_text(json.dumps({"modulation_period_s": PERIOD}))

cli(["quant", str(out), "--outdir", str(out / "results")])

found = {(r["run"], r["compound"]): r for r in csv.DictReader(open(out / "results" / "targets.csv"))}
print(f"\n{'run':7} {'compound':13} {'truth':>6} {'found':>6} {'error':>7}  status  flags")
errors = {}
for run, comp, amount in truth:
    r = found[(run, comp)]
    a = float(r["amount"]) if r["amount"] not in ("", "nan") else float("nan")
    err = 100 * (a - amount) / amount if amount > 0 else float("nan")
    errors.setdefault(comp, []).append((run, err, r["status"]))
    print(f"{run:7} {comp:13} {amount:6.2f} {a:6.2f} {err:+6.1f}%  {r['status']:6}  {r['flags']}")
print()
for comp, rows in errors.items():
    cal = [abs(e) for run, e, _s in rows if run.startswith("CAL")]
    silent = [f"{run} {e:+.0f}%" for run, e, s in rows if s == "ok" and abs(e) > 20]
    print(f"{comp:13} median calibration error {np.median(cal):4.1f}%; unflagged errors above 20%: {silent or 'none'}")
