"""Are the deuterated internal standards in the real runs, and does 🏷️ find them?

For each standard in the shipped IS library, the brightest peak of a characteristic
ion is located on the cleaned 70 eV run and its spectrum checked. The standard is
then labelled once on ATWC (as a click in the 🏷️ panel would) and matched on the
other runs by position and spectrum (`core.istd.match_run`). Finally 🧪's own peak
finder measures every standard in every run from a target table, and a draft
target table of the standards is written to real_data/target_table_istd.csv.

Di-n-octyl phthalate-d4 and n-dotriacontane-d66 are not found: the phthalate's
ions are absent, and n-C32 elutes after the end of these 76-minute runs.

Run from the repository root: python benchmarks/real_data_audit/internal_standards.py
"""

import numpy as np
from common import PERIOD, PHASE, REAL, RUNS, load

from pygcgc import PeakSettings, QuantSettings, Target, find_blobs, fold, quantify_run, remove_background
from pygcgc.core.background import storage_threshold
from pygcgc.core.istd import IstdReference, match_run, peak_spectrum, references_to_targets
from pygcgc.core.targets import write_targets_csv

# name: (quantifier, qualifiers, first-look t1 (s), t2 (s)). Perylene-d12 wraps: its true t2 is above the period.
STANDARDS = {"acenaphthene-d10": (164, (162, 160), 33.58 * 60, 3.74),
             "phenanthrene-d10": (188, (184, 160), 45.05 * 60, 4.63),
             "n-tetracosane-d50": (66, (82, 98), 60.38 * 60, 1.36),
             "chrysene-d12": (240, (236, 120), 65.58 * 60, 6.35),
             "perylene-d12": (264, (260, 132), 76.12 * 60, 1.16 + PERIOD)}
prepared = {}
for key in RUNS:
    raw = load(key)
    f = fold(raw, PERIOD, PHASE)
    bg = remove_background(raw, f)
    blobs = find_blobs(bg.run, f, PeakSettings(), floor=storage_threshold(raw) / 3.0)
    prepared[key] = (raw, bg, f, blobs)
    print(f"=== {key}")
    for name, (q, quals, _t1, _t2) in STANDARDS.items():
        img = f.ion_image(bg.run, q)
        img[f.t2 < 0.5, :] = 0
        r, c = np.unravel_index(np.argmax(img), img.shape)
        i = bg.run.scan_at(float(f.to_time(c, f.t2[r])))
        spec = bg.run.spectrum(slice(i - 3, i + 4))
        rel = " ".join(f"{m}:{100 * spec[m - bg.run.mz[0]] / spec.max():.0f}" for m in (q, *quals))
        print(f"  {name:18} t1 {(f.t1[c] + PERIOD / 2) / 60:6.2f} min, t2 {f.t2[r]:4.2f} s, "
              f"height {img[r, c]:8.0f} ({img[r, c] / np.median(img[img > 0]):6.0f} x median); relative {rel}")

print("\n🏷️ labelled on ATWC, matched on the other runs")
raw, bg, f, blobs = prepared["ATWC"]
references = []
for name, (q, quals, t1, t2) in STANDARDS.items():
    i = blobs.nearest(t1, t2 % PERIOD)
    references.append(IstdReference(name=name, cas="", run="ATWC", t1=float(blobs.t1[i]), t2=float(blobs.t2[i]),
                                    mz0=int(blobs.mz[0]), spectrum=peak_spectrum(blobs, i)))
for key in ("BTWC", "T191", "T60"):
    for m in match_run(references, key, prepared[key][3]):
        print(f"  {key:5} {m.name:18} {m.tier:9} score {m.combined_score:.2f} (spectral {m.spectral_score:.2f}), "
              f"dt1 {m.dt1:+5.1f} s, dt2 {m.dt2:+.3f} s")
print("  target rows the 🏷️ panel's 'Use as target table' would write:")
for t in references_to_targets(references):
    print(f"    {t.compound:18} quantifier {t.quant_ion}, qualifiers {t.qual_ions} at {t.qual_ratios}")

print("\n🧪 measured from a target table")
targets = [Target(name, "istd", "", q, quals, (), t1, 12.0, t2, 0.25) for name, (q, quals, t1, t2) in STANDARDS.items()]
hits = {}
for key, (raw, bg, f, _blobs) in prepared.items():
    rq = quantify_run(raw, bg.run, f, bg.noise, targets, QuantSettings())
    for name, h in rq.hits.items():
        hits.setdefault(name, []).append(h)
        print(f"  {key:5} {name:18} t1 {h.t1 / 60:6.2f} min, t2 {h.t2:5.2f} s, S/N {h.snr:7.0f}, "
              f"volume {h.quant_volume:9.3g}, flags {h.flags() or 'none'}")
rows = []
for name, (q, quals, _t1, _t2) in STANDARDS.items():
    found = [h for h in hits[name] if h.found]
    rows.append(Target(name, "istd", "", q, quals,
                       tuple(round(float(np.median([h.qual_ratios[k] for h in found])), 3) for k in range(len(quals))),
                       round(float(np.median([h.t1 for h in found])), 1), 12.0,
                       round(float(np.median([h.t2 for h in found])), 3), 0.25))
    volumes = np.array([h.quant_volume for h in found])
    print(f"  {name:18} found in {len(found)}/{len(RUNS)} runs; volume spread {100 * volumes.std() / volumes.mean():.0f}% (CV)")
print(f"\nWrote {write_targets_csv(rows, REAL / 'target_table_istd.csv')}")
