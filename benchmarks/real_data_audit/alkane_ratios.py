"""The n-alkane qualifier ratios (71/57, 85/57) on the current instrument tune.

Each run's ladder is grown with the qualifier gate off (so the measurement does not
depend on the value being measured), from the same two clicks as ladder.py. Each
hit's ratios are read at its apex, baseline-subtracted (`qualifier_ratios_at`).
The median over C13 to C28 (where the ratios are steady) comes out at 0.60 and 0.38,
below the shipped 0.66 and 0.44. The second half grows the ladders again with the gate
at each pair and judges them by the 14 eV molecular-ion check of ladder.py: the shipped
values do better, because the confounders are branched alkanes whose ratios are lower
still (for example 0.53 and 0.29 beside n-C18 on T60), and a gate centred on the
n-alkane median lets them in. So the shipped defaults stay.

Run from the repository root: python benchmarks/real_data_audit/alkane_ratios.py
"""

import numpy as np
from common import PERIOD, PHASE, RUNS, load

from pygcgc import HydrocarbonSettings, find_alkane_series, fold, remove_background

SEEDS = [(15, 1953.0, 1.13), (16, 2188.0, 1.12)]
rows, prepared = [], {}
for key in RUNS:
    run = load(key)
    f = fold(run, PERIOD, PHASE)
    clean = remove_background(run, f).run
    prepared[key] = (clean, f, load(key, "14eV"))
    series = find_alkane_series(clean, f, SEEDS, HydrocarbonSettings(enforce_qualifiers=False))
    for c, hit in sorted(series.hits.items()):
        rows.append((key, c, *hit.qual_ratios))

rows = np.array([(c, a, b) for _k, c, a, b in rows])
print(" C   71/57 (median, range)      85/57 (median, range)   runs")
for c in np.unique(rows[:, 0]).astype(int):
    r = rows[rows[:, 0] == c]
    print(f"{c:2d}   {np.median(r[:, 1]):.2f} ({r[:, 1].min():.2f}-{r[:, 1].max():.2f})"
          f"            {np.median(r[:, 2]):.2f} ({r[:, 2].min():.2f}-{r[:, 2].max():.2f})         {len(r)}")
steady = rows[(rows[:, 0] >= 13) & (rows[:, 0] <= 28)]
q71, q85 = np.median(steady[:, 1]), np.median(steady[:, 2])
print(f"\nMedian over C13-C28 ({len(steady)} hits): 71/57 = {q71:.2f}, 85/57 = {q85:.2f}")
within = np.all(np.abs(steady[:, 1:] - [q71, q85]) <= 0.3 * np.array([q71, q85]), axis=1)
print(f"Within the 30% gate of those values: {within.sum()} of {len(steady)} hits (C13-C28); "
      f"{np.sum(np.all(np.abs(rows[:, 1:] - [q71, q85]) <= 0.3 * np.array([q71, q85]), axis=1))} of {len(rows)} overall")
old = np.all(np.abs(steady[:, 1:] - [0.66, 0.44]) <= 0.3 * np.array([0.66, 0.44]), axis=1)
print(f"Within the 30% gate of the old defaults (0.66, 0.44): {old.sum()} of {len(steady)} (C13-C28)")

print("\nLadders with the gate at each pair of ratios, judged by the 14 eV molecular ion (to C21)")
for label, ratios in (("shipped 0.66/0.44", (0.66, 0.44)), (f"measured {q71:.2f}/{q85:.2f}", (round(q71, 2), round(q85, 2)))):
    confirmed = contradicted = 0
    lengths = []
    for key, (clean, f, run14) in prepared.items():
        series = find_alkane_series(clean, f, SEEDS, HydrocarbonSettings(qual_ratios=ratios))
        lengths.append(f"{key} {len(series.hits)}")
        for c, hit in series.hits.items():
            if c > 21:
                continue
            col, _ = f.to_fold(hit.t1 - 0.5 * PERIOD)
            i = run14.scan_at(float(f.to_time(col, hit.t2)))
            spec = run14.spectrum(slice(i - 4, i + 5)) - run14.spectrum(slice(i - 60, i - 40)) * 9 / 20
            m_ions = {n: spec[14 * n + 2 - run14.mz[0]] / spec.max() for n in range(max(c - 3, 7), min(c + 3, 21) + 1)}
            best = max(m_ions, key=m_ions.get)
            if m_ions[best] >= 0.008 and m_ions[best] >= 1.5 * sorted(m_ions.values())[-2]:
                confirmed += best == c
                contradicted += best != c
    print(f"  {label}: carbons {', '.join(lengths)}; confirmed {confirmed}, contradicted {contradicted}")