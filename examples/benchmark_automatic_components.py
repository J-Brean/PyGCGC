"""Run automatic overlapping-component fits in the supplied crowded window."""

from pathlib import Path
import json
import time

import numpy as np
from scipy import sparse

from pygcgc import (ComponentDetectionSettings, Run, candidate_regions,
                    find_candidates, fit_candidate_region, fold, remove_background)
from pygcgc.core.background import storage_threshold
from pygcgc.core.peaks import PeakSettings, make_frame


root = Path(__file__).resolve().parents[1]
data_dir = root / "real_data" / "blob_diagnosis"
z = np.load(data_dir / "segment.npz")
meta = json.loads((data_dir / "segment.json").read_text())
raw = Run(z["times"], z["mz"], sparse.csr_matrix(
    (z["values"], z["indices"], z["indptr"]), shape=tuple(z["shape"])),
    "automatic component benchmark")
folded = fold(raw, meta["period"], meta["phase"], t0=meta["fold_t0"])
clean = remove_background(raw, folded).run
settings = ComponentDetectionSettings()
floor = storage_threshold(raw) / 3.0

started = time.perf_counter()
candidates = find_candidates(clean, folded, settings.candidates, floor=floor)
frame = make_frame(folded, PeakSettings(wrap_below_s=settings.candidates.wrap_below_s))
regions = candidate_regions(
    candidates, (frame.rows, len(frame.mid_t1)), link_t2_s=settings.link_t2_s,
    link_modulations=settings.link_modulations, pad_t2_s=settings.pad_t2_s,
    pad_modulations=settings.pad_modulations)

focus = []
for region_index, region in enumerate(regions):
    indices = np.asarray(region.candidate_indices, dtype=int)
    if np.any((candidates.t1[indices] >= 14 * 60) & (candidates.t1[indices] <= 16 * 60) &
              (candidates.t2[indices] >= 0.6) & (candidates.t2[indices] <= 2.0)):
        focus.append((region_index, region))

groups = []
for done, (region_index, region) in enumerate(focus, start=1):
    group = fit_candidate_region(clean, folded, candidates, region, region_index,
                                 "real-focus", settings)
    groups.append(group)
    print(f"{done}/{len(focus)} region {region_index}: {len(group.components)} component(s), "
          f"unresolved={group.unresolved}, BIC gap={group.bic_gap:.1f}", flush=True)

metrics = {
    "settings": settings.to_dict(),
    "seconds": time.perf_counter() - started,
    "candidate_count_full_segment": len(candidates),
    "focus_region_count": len(groups),
    "focus_component_count": sum(len(group.components) for group in groups),
    "unresolved_group_count": sum(group.unresolved for group in groups),
    "groups": [{
        "region_index": group.region_index,
        "component_count": len(group.components),
        "bic": group.bic,
        "bic_gap": group.bic_gap,
        "weighted_residual": group.weighted_residual,
        "explained_fraction": group.explained_fraction,
        "unresolved": group.unresolved,
        "warning": group.warning,
        "components": [{
            "t1_s": component.t1, "t2_s": component.t2,
            "volume": component.volume,
            "relative_to_largest": component.relative_to_largest,
            "confidence": component.confidence,
            "top_ions": component.mz[np.argsort(component.spectrum)[::-1][:5]].tolist(),
        } for component in group.components],
    } for group in groups],
}
(data_dir / "automatic_component_metrics.json").write_text(json.dumps(metrics, indent=2))
print(json.dumps({key: value for key, value in metrics.items() if key != "groups"}, indent=2))
