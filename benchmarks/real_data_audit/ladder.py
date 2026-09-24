"""Does the ⛓️ n-alkane ladder finder work on the real runs?

Each run's ladder is grown from two clicks (C15 and C16, placed a few seconds and a
few hundredths of a second off, as a hand would), with the qualifier gate on (the
default) and off. Every carbon it assigns up to C21 is then checked against the
14 eV channel of the same acquisition, where the n-alkane molecular ion (14n + 2)
is visible after local background subtraction: "confirmed" if 14n + 2 is the
strongest of the alkane molecular ions within three carbons, "contradicted" if another's is, and
"no evidence" if none stands out. Above C21 the molecular ion is beyond the mass
range (300), so only the spacing can be judged. Finally the ATWC ladder is
propagated to the other runs, as the ⛓️ panel's Propagate does.

Run from the repository root: python benchmarks/real_data_audit/ladder.py
"""

import numpy as np
from common import PERIOD, PHASE, RUNS, load

from pygcgc import HydrocarbonSettings, find_alkane_series, fold, propagate_alkane_series, remove_background

SEEDS = [(15, 1953.0, 1.13), (16, 2188.0, 1.12)]     # t1 (s, mid-modulation), t2 (s), near n-C15 and n-C16
ladders = {}
for key in RUNS:
    run70, run14 = load(key, "70eV"), load(key, "14eV")
    f70 = fold(run70, PERIOD, PHASE)
    clean = remove_background(run70, f70).run
    for label, settings in (("gate on", HydrocarbonSettings()), ("gate off", HydrocarbonSettings(enforce_qualifiers=False))):
        series = find_alkane_series(clean, f70, SEEDS, settings)
        verdicts = {"confirmed": [], "contradicted": [], "no evidence": []}
        for c, hit in sorted(series.hits.items()):
            if c > 21:
                continue
            # 14 eV spectrum at the hit, minus the signal 0.4 to 0.6 s earlier in the same modulation
            col, _ = f70.to_fold(hit.t1 - 0.5 * PERIOD)
            i = run14.scan_at(float(f70.to_time(col, hit.t2)))
            spec = run14.spectrum(slice(i - 4, i + 5)) - run14.spectrum(slice(i - 60, i - 40)) * 9 / 20
            # Neighbouring carbons only: lower masses also carry fragment ions of heavier compounds.
            m_ions = {n: spec[14 * n + 2 - run14.mz[0]] / spec.max() for n in range(max(c - 3, 7), min(c + 3, 21) + 1)}
            best = max(m_ions, key=m_ions.get)
            if m_ions[best] < 0.008 or m_ions[best] < 1.5 * sorted(m_ions.values())[-2]:
                verdicts["no evidence"].append(c)
            elif best == c:
                verdicts["confirmed"].append(c)
            else:
                verdicts["contradicted"].append(f"C{c} (14 eV says C{best})")
        cs = sorted(series.hits)
        gaps = [c for c in range(cs[0], cs[-1] + 1) if c not in series.hits]
        spacing = np.diff([series.hits[c].t1 for c in cs]) / np.diff(cs)
        print(f"{key:5} {label:8}: C{cs[0]}-C{cs[-1]}, {len(cs)} carbons, missing {gaps or 'none'}; "
              f"to C21 confirmed {len(verdicts['confirmed'])}, contradicted {verdicts['contradicted'] or 'none'}, "
              f"no 14 eV evidence {verdicts['no evidence'] or 'none'}; spacing {spacing.min():.0f}-{spacing.max():.0f} s")
        if label == "gate on":
            ladders[key] = (series, clean, f70)

print("\nPropagating the ATWC ladder to the other runs")
reference = ladders["ATWC"][0]
for key in RUNS:
    if key == "ATWC":
        continue
    own, clean, f70 = ladders[key]
    prop = propagate_alkane_series(reference, clean, f70)
    agree = [c for c in prop.hits if c in own.hits and abs(prop.hits[c].t1 - own.hits[c].t1) <= PERIOD]
    print(f"{key:5}: {len(prop.hits)} of {len(reference.hits)} carbons found; {len(agree)} agree with the run's own "
          f"ladder to within one modulation; missing C{sorted(set(reference.hits) - set(prop.hits))}")
