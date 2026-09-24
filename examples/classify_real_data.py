"""Compound-class assignment on the real 70 eV run: reproduces every number in
real_data/CLASS_VALIDATION.md and writes real_data/class_map.png (+ .txt caption).

Uses the same cached decode of the native file as validate_real_data.py (the
first run of that script creates it). Blob detection and the alkane ladder take
about a minute. Run with the project's Python environment:
    python examples/classify_real_data.py
"""
import tempfile
from pathlib import Path

import matplotlib
import numpy as np
from scipy import sparse

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pygcgc.core.background import remove_background, storage_threshold  # noqa: E402
from pygcgc.core.chemclass import classify_blobs, load_rules  # noqa: E402
from pygcgc.core.fold import fold  # noqa: E402
from pygcgc.core.hydrocarbons import HydrocarbonSettings, find_alkane_series  # noqa: E402
from pygcgc.core.peaks import PeakSettings, find_blobs  # noqa: E402
from pygcgc.core.run import Run  # noqa: E402

PERIOD, PHASE, MIN = 10.0, 1.81, 60.0
SOLVENT_FRONT_END_MIN = 11.0        # the second solvent-front streak is at 10.4 min (HYDROCARBON_VALIDATION.md)
root = Path(__file__).resolve().parents[1] / "real_data"
cache = Path(tempfile.gettempdir()) / "pygcgc_real70.npz"
if not cache.exists():
    raise SystemExit("Run examples/validate_real_data.py first: it decodes the native file and caches it.")

z = np.load(cache)
run = Run(z["times"], z["mz"], sparse.csr_matrix((z["values"], z["indices"], z["indptr"]), shape=tuple(z["shape"])), "real70")
f = fold(run, PERIOD, PHASE)
ladder = find_alkane_series(run, f, [(13, 23.45 * MIN, 0.74), (18, 42.95 * MIN, 0.92)], HydrocarbonSettings())
blobs = find_blobs(remove_background(run, f).run, f, PeakSettings(), floor=storage_threshold(run) / 3.0)
res = classify_blobs(blobs, load_rules(), ladder)
print(f"{len(blobs)} blobs, {len(ladder)} alkanes in the ladder (C{min(ladder.hits)} to C{max(ladder.hits)})")

cls = np.array([c.cls or "unclassified" for c in res.calls])
grp = np.array([c.group or "unclassified" for c in res.calls])
vol, t1_min, t2, dt2 = blobs.volume, blobs.t1 / MIN, blobs.t2, res.features["dt2"]
empty_apex = (blobs.apex_spectra > 0).sum(axis=1) < res.settings.min_ions
print(f"{int(empty_apex.sum())} blobs ({100 * vol[empty_apex].sum() / vol.sum():.1f}% of volume) have an empty apex spectrum: footprint used")

# ------------------------------------------------------------------ tables
after = t1_min >= SOLVENT_FRONT_END_MIN
print(f"\nvolume by class, all blobs / after the solvent front (t1 >= {SOLVENT_FRONT_END_MIN:g} min):")
total_after = vol[after].sum()
for name, n, v, pct in res.summary(vol):
    sel = (cls == name) & after
    print(f"  {name:28s} n={n:5d} {pct:5.1f}%   |  n={int(sel.sum()):5d} {100 * vol[sel].sum() / total_after:5.1f}%")
print("tiers:", {t: int(sum(c.tier == t for c in res.calls)) for t in ("confident", "probable", "possible", "unclassified")})

print("\nalkane ladder: blob with the most m/z 57 signal within the search box of each ladder carbon")
i57 = 57 - int(blobs.mz[0])
hits = 0
for c in sorted(ladder.hits):
    h = ladder.hits[c]
    box = np.flatnonzero(np.hypot((blobs.t1 - h.t1) / 6.0, (blobs.t2 - h.t2) / 0.3) < 2.5)
    i = int(box[np.argmax(np.clip(blobs.spectra[box, i57], 0, None))])
    hits += cls[i] == "alkane"
    print(f"  C{c}: blob {i} -> {cls[i]} ({res.calls[i].score:.2f}, {res.calls[i].tier})")
print(f"  {hits} of {len(ladder)} are alkane")

print("\nt2 above the alkane line (s), blobs after the solvent front with the ladder covering them:")
order = ["alkane", "alkene", "cycloalkane", "alkylbenzene", "indane / indene / tetralin", "naphthalene", "PAH",
         "N-heteroaromatic", "aromatic oxygenate", "alcohol / ether", "ketone / aldehyde", "carboxylic acid", "ester",
         "aliphatic amine", "siloxane", "phthalate"]
box_data, box_names = [], []
for name in order:
    sel = (cls == name) & after & ~np.isnan(dt2)
    if sel.sum() >= 5:
        p = np.percentile(dt2[sel], [10, 50, 90])
        print(f"  {name:28s} n={int(sel.sum()):4d}  p10 {p[0]:5.2f}  median {p[1]:5.2f}  p90 {p[2]:5.2f}")
        box_data.append(dt2[sel])
        box_names.append(f"{name} ({int(sel.sum())})")

# ------------------------------------------------------------------ figure
plt.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.6, "xtick.direction": "in",
                     "ytick.direction": "in", "xtick.major.size": 3, "ytick.major.size": 3, "xtick.top": True,
                     "ytick.right": True, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
                     "legend.fontsize": 7, "legend.frameon": False})
style = {"hydrocarbon": ("#0072B2", "o"), "aromatic": ("#E69F00", "s"), "oxygen": ("#009E73", "^"),
         "nitrogen": ("#D55E00", "D"), "contaminant": ("#999999", "x"), "unclassified": ("#CCCCCC", ".")}
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(183 / 25.4, 95 / 25.4), gridspec_kw={"width_ratios": [1.15, 1]})
fig.subplots_adjust(left=0.06, right=0.995, bottom=0.13, top=0.99, wspace=0.62)
size = 3 + 4 * np.clip(np.log10(np.clip(vol, 1, None)) - 4, 0, None)
for g in ("unclassified", "contaminant", "hydrocarbon", "aromatic", "oxygen", "nitrogen"):
    sel = grp == g
    colour, marker = style[g]
    kw = {"linewidths": 0.4} if marker == "x" else {"edgecolors": "k", "linewidths": 0.2}
    ax_a.scatter(t1_min[sel], t2[sel], s=size[sel], c=colour, marker=marker, label=f"{g} ({int(sel.sum())})", **kw)
ax_a.set_xlabel("first-dimension retention time (min)")
ax_a.set_ylabel("second-dimension retention time (s)")
ax_a.set_xlim(5, 72)
ax_a.set_ylim(0, 10)
ax_a.legend(loc="upper right", ncol=2, markerscale=1.2, handletextpad=0.1, columnspacing=0.6, labelspacing=0.2,
            borderaxespad=0.3)
ax_a.text(0.02, 0.03, "(a)", transform=ax_a.transAxes, va="bottom")

bp = ax_b.boxplot(box_data, orientation="horizontal", whis=(10, 90), showfliers=False, widths=0.6,
                  medianprops={"color": "k", "linewidth": 0.8}, boxprops={"linewidth": 0.6}, whiskerprops={"linewidth": 0.6},
                  capprops={"linewidth": 0.6})
ax_b.set_yticks(range(1, len(box_names) + 1))
ax_b.set_yticklabels(box_names)
ax_b.invert_yaxis()
ax_b.set_xlabel("t2 above the n-alkane line (s)")
ax_b.set_xlim(-0.5, 10.5)
ax_b.set_xticks(range(0, 11, 2))
ax_b.axvline(0, color="k", lw=0.4, ls=":")
ax_b.text(0.98, 0.03, "(b)", transform=ax_b.transAxes, ha="right", va="bottom")
out = root / "class_map.png"
fig.savefig(out, dpi=300)
print(f"\nwrote {out}")

caption = f"""Compound-class assignment of {len(blobs)} blobs from the untreated 70 eV run (Repsol3 test3, 1 uL injection, m/z 38 to 300, modulation period {PERIOD:g} s).

(a) Blob apexes in the first- and second-dimension retention plane, coloured and shaped by chemical group. Marker area increases with blob volume (logarithmic, m/z summed over the footprint). (b) Distribution of the second-dimension retention time above the n-alkane line (t2 minus the t2 of the C12 to C28 ladder at the blob's own carbon number) for each class with at least five blobs; the number of blobs is given in brackets. Boxes show the interquartile range, the bar the median and whiskers the 10th and 90th percentiles. The dotted line marks the alkane line. Blobs earlier than {SOLVENT_FRONT_END_MIN:g} min (solvent front) and blobs outside the ladder's carbon range are omitted from (b).

Classes were assigned from each blob's apex mass spectrum (footprint spectrum where the apex spectrum was empty, {int(empty_apex.sum())} blobs) with the rules in chemical_classes.toml. The alkane ladder was found from two seeds (C13 and C18). The t2 offset enters the rules only as supporting evidence, which can raise a class score from 0.5 to 1.0; no class requires it, so panel (b) is a consistency check and not an independent measurement. The alkane line is the ladder itself, so the alkane median lies at zero by construction only if the alkane blobs sit on the ladder. The dense block of blobs at t2 above 4 s and t1 below 20 min in (a) is not separated from the rest; it produces the upper whisker of the alkylbenzene box in (b). The S-heteroaromatic class was disabled. Blobs that met no class at a score of 0.3 are shown as unclassified.
"""
out.with_suffix(".txt").write_text(caption, encoding="utf-8")
print(f"wrote {out.with_suffix('.txt')}")
