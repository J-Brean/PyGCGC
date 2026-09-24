"""Does the 🫧 siloxane remover do the right thing on the real runs?

The remover scores the legacy 🔍 blobs, so this does too. Independent labels come
from each blob's own summed spectrum: "siloxane" if at least half its volume sits
on characteristic Si fragments and its base peak is one of them, "clearly not" if
under 10% does, "mixed" otherwise. For both measured references it reports what
the default threshold (0.3) flags, a threshold scan on the clearly labelled blobs,
and how much non-siloxane signal (m/z 57, 91, 128) blocking the flagged footprints
destroys. The per-pixel method is shown for comparison.

Run from the repository root: python benchmarks/real_data_audit/siloxane.py
"""

import time
from dataclasses import replace

import numpy as np
from common import PERIOD, PHASE, load

from pygcgc import PeakSettings, find_blobs, fold, remove_background
from pygcgc.core.background import storage_threshold
from pygcgc.core.siloxane import (CYCLIC_REFERENCE, SiloxaneSettings, remove_siloxane_blobs, remove_siloxanes,
                                  siloxane_scores)

SI_IONS = [73, 96, 133, 147, 191, 193, 207, 208, 209, 221, 267, 281, 282, 283]
SI_BASE = {73, 147, 207, 221, 267, 281}
WATCH = (57, 91, 128, 73, 207, 281)
for key in ("ATWC", "T60"):
    raw = load(key)
    f = fold(raw, PERIOD, PHASE)
    clean = remove_background(raw, f).run
    t = time.perf_counter()
    blobs = find_blobs(clean, f, PeakSettings(), floor=storage_threshold(raw) / 3.0)
    mz0 = int(blobs.mz[0])
    spec = np.clip(blobs.spectra, 0, None)
    si_frac = spec[:, [m - mz0 for m in SI_IONS]].sum(axis=1) / np.maximum(spec.sum(axis=1), 1e-12)
    base = blobs.mz[np.argmax(spec, axis=1)]
    is_si = (si_frac >= 0.5) & np.isin(base, list(SI_BASE))
    is_other = si_frac < 0.10
    print(f"\n=== {key}: {len(blobs)} legacy blobs ({time.perf_counter() - t:.0f} s), covering "
          f"{100 * (blobs.labels > 0).mean():.0f}% of the image; {is_si.sum()} siloxane-like, "
          f"{is_other.sum()} clearly not, {len(blobs) - is_si.sum() - is_other.sum()} mixed")
    before = {m: float(clean.xic(m).clip(0).sum()) for m in WATCH}
    for name, s in (("bleed (default)", SiloxaneSettings()), ("cyclic peak", SiloxaneSettings.from_reference(CYCLIC_REFERENCE))):
        score = siloxane_scores(blobs, s)
        flagged = score >= s.score_threshold
        scan = "; ".join(f"{thr}: {np.sum((score >= thr) & is_si)}/{is_si.sum()} siloxane, "
                         f"{np.sum((score >= thr) & is_other)} not" for thr in (0.05, 0.1, 0.15, 0.2, 0.3, 0.5))
        res = remove_siloxane_blobs(clean, f, blobs, s)
        after = {m: float(res.run.xic(m).clip(0).sum()) for m in WATCH}
        print(f"{name:16}: flags {flagged.sum()} at 0.3 ({np.sum(flagged & is_si)} siloxane-like, "
              f"{np.sum(flagged & ~is_si & ~is_other)} mixed, {np.sum(flagged & is_other)} clearly not)")
        print(f"{'':16}  threshold scan {scan}")
        print(f"{'':16}  blocking removes " + ", ".join(f"m/z {m} {100 * (1 - after[m] / before[m]):.1f}%" for m in WATCH))
    # Recommended in ROADMAP.md: score against both references, keep the larger, block siloxane ions only.
    score = np.maximum(siloxane_scores(blobs, SiloxaneSettings()),
                       siloxane_scores(blobs, SiloxaneSettings.from_reference(CYCLIC_REFERENCE)))
    scan = "; ".join(f"{thr}: {np.sum((score >= thr) & is_si)}/{is_si.sum()} siloxane, "
                     f"{np.sum((score >= thr) & is_other)} not" for thr in (0.1, 0.15, 0.2, 0.3))
    print(f"{'both, larger':16}: threshold scan {scan}")
    # remove_siloxane_blobs scores against one reference, so hand it only the chosen blobs' footprints
    # and a zero threshold.
    picked = replace(blobs, labels=np.where(np.isin(blobs.labels, np.flatnonzero(score >= 0.15) + 1), blobs.labels, 0))
    res = remove_siloxane_blobs(clean, f, picked, SiloxaneSettings(score_threshold=0.0, block_all_ions=False))
    after = {m: float(res.run.xic(m).clip(0).sum()) for m in WATCH}
    print(f"{'':16}  at 0.15, siloxane ions only, removes "
          + ", ".join(f"m/z {m} {100 * (1 - after[m] / before[m]):.1f}%" for m in WATCH))
    res = remove_siloxanes(clean, f, SiloxaneSettings(method="pixels"))
    after = {m: float(res.run.xic(m).clip(0).sum()) for m in WATCH}
    print(f"{'per-pixel fit':16}: removes " + ", ".join(f"m/z {m} {100 * (1 - after[m] / before[m]):.1f}%" for m in WATCH))
