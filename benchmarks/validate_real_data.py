"""Reproduce the real-data measurements behind the siloxane, alkane, marker,
diagnostic-pair and phenanthrene/anthracene defaults.

Reads the untreated 70 eV native file in real_data/ (period 10 s, phase 1.81 s, as
supplied), prints every number quoted in docs/validation/hydrocarbons-2026-09-20.md and
changes no production setting. The first run decodes the 446 MB file (a few
minutes); the sparse matrix is cached in the system temp folder afterwards.
Run with the project's Python environment: python benchmarks/validate_real_data.py
"""
import tempfile
from pathlib import Path

import numpy as np
from scipy import sparse

from pygcgc.core.background import remove_background, storage_threshold
from pygcgc.core.fold import fold
from pygcgc.core.hydrocarbons import (HydrocarbonSettings, find_alkane_series, find_homologous_series)
from pygcgc.core.identify import IdentifySettings, LibraryEntry, build_library, classify, retention_index, search
from pygcgc.core.io.markes_lsc import read_lsc
from pygcgc.core.peaks import PeakSettings, find_blobs
from pygcgc.core.pixelspec import spectrum_at
from pygcgc.core.run import Run
from pygcgc.core.siloxane import SiloxaneSettings, remove_siloxanes

PERIOD, PHASE, MIN = 10.0, 1.81, 60.0
root = Path(__file__).resolve().parents[1] / "real_data"
cache = Path(tempfile.gettempdir()) / "pygcgc_real70.npz"

if cache.exists():
    z = np.load(cache)
    run = Run(z["times"], z["mz"], sparse.csr_matrix((z["values"], z["indices"], z["indptr"]),
                                                     shape=tuple(z["shape"])), "real70")
else:
    run = read_lsc(next(root.glob("*_70eV.lsc")))
    d = run.data
    np.savez(cache, times=run.times, mz=run.mz, values=d.data, indices=d.indices, indptr=d.indptr,
             shape=np.array(d.shape))
f = fold(run, PERIOD, PHASE)
t1_mid = f.t1 + 0.5 * PERIOD
print(f"{run.n_scans} scans, m/z {run.mz[0]}-{run.mz[-1]}, {f.n_cols} modulations of {f.n_rows} rows\n")

# ---------------------------------------------------------------- siloxane
print("== column bleed reference (each ion's 20 s median regressed on m/z 207, 40-71 min)")
meds = []
for t in np.arange(40 * MIN, 71 * MIN, 20.0):
    sel = (run.times >= t) & (run.times < t + 20.0)
    meds.append(np.median(run.data[sel].toarray(), axis=0))
meds = np.array(meds)
x = meds[:, 207 - run.mz[0]]
design = np.column_stack([x, np.ones_like(x)])
rows = []
for j, m in enumerate(run.mz):
    coef, *_ = np.linalg.lstsq(design, meds[:, j], rcond=None)
    ss = ((meds[:, j] - meds[:, j].mean()) ** 2).sum()
    r2 = 1 - ((meds[:, j] - design @ coef) ** 2).sum() / ss if ss > 0 else 0.0
    rows.append((int(m), 100 * coef[0], r2))
kept = [(m, s) for m, s, r2 in rows if r2 >= 0.9 and s >= 2.0]
print("  kept (R2 >= 0.9, slope >= 2):", ", ".join(f"{m}:{s:.1f}" for m, s in sorted(kept)))

print("\n== removal on the raw run with the default (bleed) reference")
i207 = f.image(run.xic(207))
res = remove_siloxanes(run, f, SiloxaneSettings())
n207 = f.image(res.run.xic(207))
late = (t1_mid > 55 * MIN) & (t1_mid < 70 * MIN)
print(f"  m/z 207, 55-70 min: {100 * (1 - n207[:, late].sum() / i207[:, late].sum()):.1f}% removed; {res.summary()}")

# ------------------------------------------------------------ alkane ladder
print("\n== n-alkane ladder from two clicked seeds (C13 and C18)")
seeds = [(13, 23.45 * MIN, 0.74), (18, 42.95 * MIN, 0.92)]
ladder = find_alkane_series(run, f, seeds, HydrocarbonSettings())
prev = None
r71, r85 = [], []
for c in sorted(ladder.hits):
    h = ladder.hits[c]
    step = "" if prev is None else f"  step {(h.t1 - prev) / MIN:.2f} min"
    prev = h.t1
    print(f"  C{c}: t1 {h.t1 / MIN:.3f} min, t2 {h.t2:.2f} s, 71/57 {h.qual_ratios[0]:.2f}, 85/57 {h.qual_ratios[1]:.2f}{step}")
    if c >= 13 and c != 17:            # C17 co-elutes with pristane
        r71.append(h.qual_ratios[0])
        r85.append(h.qual_ratios[1])
print(f"  median ratios, C13-C28 without C17: 71/57 {np.median(r71):.3f}, 85/57 {np.median(r85):.3f}")
x57 = run.xic(57)
lo, hi = np.searchsorted(run.times, [42.95 * MIN - 1.5, 42.95 * MIN + 1.5])
i = lo + int(np.argmax(x57[lo:hi]))
print(f"  the old time-only search at the C18 seed (+-1.5 s of the clicked time) takes the tallest m/z 57 pulse, "
      f"which is at t2 {f.to_fold(run.times[i])[1]:.1f} s with height {x57[i]:.0f}, not on the alkane line")

# ------------------------------------------------------ pristane / phytane
print("\n== pristane and phytane: t1 centroid of the m/z 183 excess against the alkane's m/z 57 centroid")
i57, i183 = f.image(run.xic(57)), f.image(run.xic(183))
t2 = np.arange(f.n_rows) * f.row_pitch
for cn in (17, 18):
    h = ladder.hits[cn]
    spacing = (ladder.t1_of(cn + 1) - ladder.t1_of(cn - 1)) / 2
    cols = np.flatnonzero(np.abs(t1_mid - h.t1) <= 40)
    band = (t2 >= h.t2 - 0.25) & (t2 <= h.t2 + 0.25)
    p57, p183 = i57[band][:, cols].sum(axis=0), i183[band][:, cols].sum(axis=0)
    excess = np.clip(p183 - 0.015 * p57, 0, None)
    off = (excess * t1_mid[cols]).sum() / excess.sum() - (p57 * t1_mid[cols]).sum() / p57.sum()
    print(f"  C{cn}: offset {off:+.1f} s = {off / spacing:+.3f} C (one column = {PERIOD / spacing:.3f} C)")

# ------------------------------------------------ homologous diagnostic pairs
print("\n== diagnostic ion pairs on the blobs of the whole cleaned run (default detector)")
clean = remove_background(run, f).run
blobs = find_blobs(clean, f, PeakSettings(), floor=storage_threshold(run) / 3.0)
mz0 = int(blobs.mz[0])
print(f"  {len(blobs)} blobs")
for pair in ((57, 71), (91, 105), (82, 83), (142, 156), (58, 59)):
    share = np.clip(blobs.spectra[:, pair[0] - mz0], 0, None) / np.clip(blobs.volume, 1e-12, None)
    hs = find_homologous_series(blobs, pair, ladder, name=str(pair))
    cands = int((share >= 0.1).sum())
    print(f"  {pair}: {cands} candidate blobs (major share >= 10%) -> "
          + ("no series" if hs is None else f"{len(hs.members)} members, carbon-equivalents "
             + " ".join(f"{m['carbon']:.2f}" for m in hs.members)))

# ------------------------------------------------ phenanthrene / anthracene
print("\n== C14H10 candidates (m/z 178 base) and RI")
peak = spectrum_at(run, f, 42.451 * MIN, 2.98)
ri = retention_index(ladder, 42.451 * MIN)
per_col = 100.0 * PERIOD / (ladder.t1_of(19) - ladder.t1_of(18))
print(f"  clean 178 peak at t1 42.45 min, t2 2.98 s: RI {ri:.1f}; one column = {per_col:.1f} RI units")
big = spectrum_at(run, f, 44.451 * MIN, 3.26)
print(f"  the larger 178 peak at 44.45 min, t2 3.26 s (RI {retention_index(ladder, 44.451 * MIN):.0f}) has m/z 188 at "
      f"{100 * big[188 - run.mz[0]] / big[178 - run.mz[0]]:.0f}% of 178 (a d10 standard co-elutes)")
keep = peak > 0.01 * peak.max()
qmz, qint = run.mz[keep].astype(np.int64), peak[keep]
lib = build_library([LibraryEntry("phenanthrene", qmz, qint.copy(), ri=1777.0),      # RIs assumed from the literature
                     LibraryEntry("anthracene", qmz, qint.copy(), ri=1793.0)])       # for a 5%-phenyl phase
for label, s in (("default", IdentifySettings()),
                 ("ri_penalty_full 30", IdentifySettings(ri_window=40.0, ri_penalty_full=30.0))):
    print(f"  {label}:")
    for q in (1775, 1780, 1785, 1790, 1795):
        hits = search(qmz, qint, lib, query_ri=float(q), settings=s)
        print(f"    query RI {q}: {hits[0].name}, margin {hits[0].score - hits[1].score:.3f}, tier {classify(hits, s)}")
