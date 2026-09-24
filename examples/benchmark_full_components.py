"""Time MS-aware separation on the supplied full 70 eV acquisition."""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import time

from pygcgc import ComponentDetectionSettings, detect_components, fold, remove_background
from pygcgc.core.background import storage_threshold
from pygcgc.core.io.markes_lsc import read_lsc


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int,
                        default=ComponentDetectionSettings().workers)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    source = next((root / "real_data").glob("*_70eV.lsc"))
    timings = {}

    started = time.perf_counter()
    run = read_lsc(source)
    timings["read_seconds"] = time.perf_counter() - started
    print(f"Read {run.n_scans:,} scans in {timings['read_seconds']:.1f} s", flush=True)

    stage = time.perf_counter()
    folded = fold(run, 10.0, 1.81)
    timings["fold_seconds"] = time.perf_counter() - stage

    stage = time.perf_counter()
    clean = remove_background(run, folded).run
    timings["background_seconds"] = time.perf_counter() - stage
    print(f"Removed background in {timings['background_seconds']:.1f} s", flush=True)

    stage = time.perf_counter()
    candidates, regions, components = detect_components(
        clean, folded, source.stem,
        replace(ComponentDetectionSettings(), workers=args.workers),
        floor=storage_threshold(run) / 3.0)
    timings["detection_seconds"] = time.perf_counter() - stage
    timings["total_seconds"] = time.perf_counter() - started

    result = {
        "source": source.name,
        "scans": run.n_scans,
        "mass_channels": run.n_mz,
        "candidates": len(candidates),
        "regions": len(regions),
        "components": len(components.components),
        "unresolved_regions": sum(group.unresolved for group in components.groups),
        "workers": args.workers,
        **timings,
    }
    suffix = "" if args.workers == 1 else f"_workers{args.workers}"
    output = root / "real_data" / "blob_diagnosis" / f"full_component_benchmark{suffix}.json"
    output.write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
