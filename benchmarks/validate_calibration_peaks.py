"""Check the default detector around a manually labelled calibration peak.

This is a recovery and cross-level consistency check, not an accuracy estimate
for the whole chromatogram. The reference is a user's saved label, not a vendor
answer key. Raw instrument files and the existing project are read only.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from pygcgc import fold, read_run, remove_background
from pygcgc.core.background import storage_threshold
from pygcgc.core.components import ComponentDetectionSettings, detect_components
from pygcgc.core.io import crop_run
from pygcgc.core.istd import IstdReference, find_istd_on_run
from pygcgc.core.peaks import find_blobs
from pygcgc.core.peaktable import respectrum_components


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=4)
    args = parser.parse_args()
    reference_bytes = (args.folder / "project.pygcgc.json").read_bytes()
    project = json.loads(reference_bytes.decode("utf-8"))
    references = [IstdReference(**{**r, "spectrum": np.asarray(r["spectrum"])})
                  for r in project["istd"]["references"].values()]
    if not references:
        parser.error("The saved project must contain at least one manually labelled standard")
    settings = replace(ComponentDetectionSettings(), workers=1)
    report = {"scope": "120-second crops around saved manual labels; not whole-run accuracy",
              "reference_project_sha256": hashlib.sha256(reference_bytes).hexdigest(),
              "period_s": project["period"], "settings": settings.to_dict(),
              "results": []}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    for entry in project["runs"][:args.limit]:
        path = args.folder / Path(entry["path"]).name
        print(f"Reading {path.name}", flush=True)
        raw = read_run(path)
        full_fold = fold(raw, project["period"], entry["phase"])
        for ref_number, ref in enumerate(references, start=1):
            started = time.perf_counter()
            crop = crop_run(raw, ref.t1 - 60, ref.t1 + 60)
            f = fold(crop, project["period"], entry["phase"], t0=full_fold.spec.t0)
            cleaned = remove_background(crop, f).run
            floor = storage_threshold(raw) / 3
            candidates, regions, components = detect_components(
                cleaned, f, path.stem, settings, floor=floor,
                progress=lambda done, total: print(f"  {done}/{total} fit windows", flush=True))
            component_path = args.out.with_name(f"{path.stem}_reference_{ref_number:02d}_components.json")
            component_path.write_text(json.dumps(components.to_dict()), encoding="utf-8")
            peaks, _, _ = respectrum_components(cleaned, f, components)
            legacy = find_blobs(cleaned, f, floor=floor)
            row = {"run": path.name, "standard": ref.name, "candidate_count": len(candidates),
                   "source_sha256": raw.meta.get("source_sha256"),
                   "components_file": component_path.name,
                   "component_count": len(peaks), "legacy_count": len(legacy),
                   "unresolved_groups": sum(g.unresolved for g in components.groups),
                   "seconds": time.perf_counter() - started}
            near = np.flatnonzero((np.abs(peaks.t1 - ref.t1) <= 30) & (np.abs(peaks.t2 - ref.t2) <= .5))
            row["nearby_components"] = [
                {"id": peaks.component_id[j], "t1": peaks.t1[j], "t2": peaks.t2[j],
                 "volume": peaks.volume[j], "confidence": peaks.confidence[j],
                 "unresolved": bool(peaks.unresolved[j])} for j in near]
            for method, peak_list in (("ms-aware", peaks), ("legacy", legacy)):
                hit = find_istd_on_run(ref, path.stem, peak_list)
                row[method] = {"found": hit.found, "spectral_score": hit.spectral_score,
                               "t1": hit.t1, "t2": hit.t2, "reason": hit.reason,
                               "passing_candidates": hit.n_passing, "needs_review": hit.needs_review,
                               "peak_unresolved": hit.peak_unresolved}
            report["results"].append(row)
            args.out.write_text(json.dumps(report, indent=2), encoding="utf-8")
            print(json.dumps(row), flush=True)
        del raw


if __name__ == "__main__":
    main()
