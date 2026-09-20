"""Compare MS-aware candidate perception on the cached real-data segment."""

from dataclasses import replace
from pathlib import Path
import json
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse

from pygcgc import (CandidateSettings, Run, find_candidates, fold,
                    remove_background)
from pygcgc.core.background import storage_threshold


root = Path(__file__).resolve().parents[1]
data_dir = root / "real_data" / "blob_diagnosis"
z = np.load(data_dir / "segment.npz")
meta = json.loads((data_dir / "segment.json").read_text())
run = Run(z["times"], z["mz"], sparse.csr_matrix(
    (z["values"], z["indices"], z["indptr"]), shape=tuple(z["shape"])),
    "candidate benchmark")
folded = fold(run, meta["period"], meta["phase"], t0=meta["fold_t0"])
clean = remove_background(run, folded).run
raw_floor = storage_threshold(run) / 3.0

base = CandidateSettings()
cases = {
    "default": base,
    "strict": replace(base, seed_snr=8.0, single_ion_snr=25.0),
    "three_ion": replace(base, min_support_ions=3, single_ion_snr=25.0),
}
image = folded.tic_image(clean)
t1 = folded.t1 + 0.5 * folded.spec.period
t2 = np.arange(folded.n_rows) * folded.row_pitch
focus_bounds = (14 * 60, 16 * 60, 0.6, 2.0)
background_bounds = (14 * 60, 16 * 60, 5.0, 9.0)

fig, axes = plt.subplots(len(cases), 2, figsize=(8, 8), constrained_layout=True)
metrics = {"floor": raw_floor, "cases": {}}
for row_axes, (name, settings) in zip(axes, cases.items()):
    started = time.perf_counter()
    candidates = find_candidates(clean, folded, settings, floor=raw_floor)
    elapsed = time.perf_counter() - started
    in_focus = ((candidates.t1 >= focus_bounds[0]) & (candidates.t1 <= focus_bounds[1]) &
                (candidates.t2 >= focus_bounds[2]) & (candidates.t2 <= focus_bounds[3]))
    in_background = ((candidates.t1 >= background_bounds[0]) &
                     (candidates.t1 <= background_bounds[1]) &
                     (candidates.t2 >= background_bounds[2]) &
                     (candidates.t2 <= background_bounds[3]))
    metrics["cases"][name] = {
        "seconds": elapsed,
        "total_candidates": len(candidates),
        "focus_candidates": int(in_focus.sum()),
        "background_candidates": int(in_background.sum()),
        "focus_support_counts": candidates.support_count[in_focus].tolist(),
        "focus_strongest_snr": candidates.strongest_snr[in_focus].tolist(),
        "focus_t1_s": candidates.t1[in_focus].tolist(),
        "focus_t2_s": candidates.t2[in_focus].tolist(),
    }
    print(name, metrics["cases"][name], flush=True)
    for ax, bounds, selected in ((row_axes[0], focus_bounds, in_focus),
                                 (row_axes[1], background_bounds, in_background)):
        mesh = ax.pcolormesh(t1 / 60, t2, np.log10(1 + np.maximum(image, 0)),
                             shading="nearest", cmap="magma", rasterized=True)
        ax.scatter(candidates.t1[selected] / 60, candidates.t2[selected],
                   s=14, facecolors="none", edgecolors="cyan", linewidths=.7)
        ax.set(xlim=(bounds[0] / 60, bounds[1] / 60), ylim=bounds[2:],
               xlabel="t1 (min)", ylabel="t2 (s)", title=name)

fig.colorbar(mesh, ax=axes, label="log10(1 + cleaned TIC)", shrink=.5)
fig.savefig(data_dir / "candidate_comparison.png", dpi=200)
plt.close(fig)
(data_dir / "candidate_metrics.json").write_text(json.dumps(metrics, indent=2))
