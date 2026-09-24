"""Benchmark automatic overlapping-component fits in the crowded real window."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import replace
import json
from pathlib import Path
from types import SimpleNamespace
import time

import numpy as np
from scipy import sparse

from pygcgc import (ComponentDetectionSettings, Run, candidate_cube, candidate_regions,
                    find_candidates, fold, remove_background)
from pygcgc.core.background import storage_threshold
from pygcgc.core.components import (_fit_candidate_cubes, _limit_component_worker_threads,
                                    _subset_cube, propose_candidate_group)
from pygcgc.core.peaks import PeakSettings, make_frame
from pygcgc.core.overlap import initialise_overlap


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--screening-ions", type=int, default=8)
    parser.add_argument("--max-ions", type=int, default=30)
    parser.add_argument("--max-sigma-row", type=float, default=3.0)
    parser.add_argument("--max-nfev", type=int,
                        default=ComponentDetectionSettings().fast_max_nfev)
    parser.add_argument("--starts", type=int, default=2)
    parser.add_argument("--all-regions", action="store_true")
    parser.add_argument("--proposals", action="store_true")
    parser.add_argument("--screen-counts", action=argparse.BooleanOptionalAction,
                        default=None)
    parser.add_argument("--threads", action="store_true")
    parser.add_argument("--widths", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    data_dir = root / "real_data" / "blob_diagnosis"
    z = np.load(data_dir / "segment.npz")
    meta = json.loads((data_dir / "segment.json").read_text())
    raw = Run(z["times"], z["mz"], sparse.csr_matrix(
        (z["values"], z["indices"], z["indptr"]), shape=tuple(z["shape"])),
        "automatic component benchmark")
    folded = fold(raw, meta["period"], meta["phase"], t0=meta["fold_t0"])
    clean = remove_background(raw, folded).run
    defaults = ComponentDetectionSettings()
    settings = replace(defaults, screening_ions=args.screening_ions,
                       max_ions=args.max_ions, n_starts=args.starts,
                       workers=args.workers,
                       fast_max_sigma_row=args.max_sigma_row,
                       fast_max_nfev=args.max_nfev,
                       screen_component_counts=(defaults.screen_component_counts
                                                if args.screen_counts is None else
                                                args.screen_counts))
    floor = storage_threshold(raw) / 3.0

    started = time.perf_counter()
    candidates = find_candidates(clean, folded, settings.candidates, floor=floor)
    frame = make_frame(folded, PeakSettings(wrap_below_s=settings.candidates.wrap_below_s))
    regions = candidate_regions(
        candidates, (frame.rows, len(frame.mid_t1)), link_t2_s=settings.link_t2_s,
        link_modulations=settings.link_modulations, pad_t2_s=settings.pad_t2_s,
        pad_modulations=settings.pad_modulations)
    candidate_seconds = time.perf_counter() - started

    if args.widths:
        widths = []
        for region in regions:
            if len(region.candidate_indices) != 1:
                continue
            index = region.candidate_indices[0]
            cube = candidate_cube(clean, folded, candidates, region, max_ions=8)
            hint = SimpleNamespace(row=float(candidates.row[index] - region.row0),
                                   col=float(candidates.col[index] - region.col0))
            peak = initialise_overlap(cube, 1, apex_hint=[hint]).components[0]
            widths.append((peak.sigma_row, peak.sigma_col,
                           float(candidates.score[index])))
        widths = np.asarray(widths)
        print(json.dumps({
            "single_candidate_regions": len(widths),
            "sigma_row_quantiles": np.quantile(widths[:, 0], [.1, .25, .5, .75, .9]).tolist(),
            "sigma_col_quantiles": np.quantile(widths[:, 1], [.1, .25, .5, .75, .9]).tolist(),
            "strongest_20": widths[np.argsort(widths[:, 2])[-20:]].tolist(),
        }, indent=2))
        return

    if args.proposals:
        groups = [propose_candidate_group(clean, folded, candidates, region, index,
                                          "real-focus", settings)
                  for index, region in enumerate(regions)]
        region_stats = sorted(({
            "region": index, "candidates": len(region.candidate_indices),
            "rows": region.row1 - region.row0, "cols": region.col1 - region.col0,
            "pixels": (region.row1 - region.row0) * (region.col1 - region.col0),
        } for index, region in enumerate(regions)),
                              key=lambda item: item["pixels"], reverse=True)
        print(json.dumps({
            "seconds": time.perf_counter() - started,
            "candidate_seconds": candidate_seconds,
            "region_count": len(groups),
            "component_count": sum(len(group.components) for group in groups),
            "largest_regions": region_stats[:15],
        }, indent=2))
        return

    focus = []
    for region_index, region in enumerate(regions):
        indices = np.asarray(region.candidate_indices, dtype=int)
        in_focus = np.any((candidates.t1[indices] >= 14 * 60) &
                          (candidates.t1[indices] <= 16 * 60) &
                          (candidates.t2[indices] >= 0.6) &
                          (candidates.t2[indices] <= 2.0))
        if args.all_regions or in_focus:
            ordered = sorted(region.candidate_indices,
                             key=lambda item: float(candidates.score[item]), reverse=True)
            hints = [SimpleNamespace(row=float(candidates.row[item] - region.row0),
                                     col=float(candidates.col[item] - region.col0))
                     for item in ordered]
            full = candidate_cube(clean, folded, candidates, region,
                                  max_ions=settings.max_ions)
            screening = _subset_cube(full, settings.screening_ions)
            focus.append((screening, full, hints, region_index, "real-focus", settings))

    fit_started = time.perf_counter()
    groups = []
    if args.workers <= 1:
        for payload in focus:
            group = _fit_candidate_cubes(*payload)
            groups.append(group)
            _print_group(len(groups), len(focus), group)
    else:
        executor_type = ThreadPoolExecutor if args.threads else ProcessPoolExecutor
        executor_kw = ({} if args.threads else
                       {"initializer": _limit_component_worker_threads})
        with executor_type(max_workers=args.workers, **executor_kw) as executor:
            futures = [executor.submit(_fit_candidate_cubes, *payload) for payload in focus]
            for future in as_completed(futures):
                group = future.result()
                groups.append(group)
                _print_group(len(groups), len(focus), group)
    groups.sort(key=lambda group: group.region_index)
    fit_seconds = time.perf_counter() - fit_started

    metrics = {
        "settings": settings.to_dict(), "seconds": time.perf_counter() - started,
        "candidate_seconds": candidate_seconds, "fit_seconds": fit_seconds,
        "workers": args.workers, "candidate_count_full_segment": len(candidates),
        "focus_region_count": len(groups),
        "focus_component_count": sum(len(group.components) for group in groups),
        "unresolved_group_count": sum(group.unresolved for group in groups),
        "groups": [_group_metrics(group) for group in groups],
    }
    suffix = f"_row{args.max_sigma_row:g}" if args.max_sigma_row != 3.0 else ""
    suffix += f"_nfev{args.max_nfev}" if args.max_nfev != 24 else ""
    (data_dir / f"automatic_component_metrics{suffix}.json").write_text(
        json.dumps(metrics, indent=2))
    print(json.dumps({key: value for key, value in metrics.items() if key != "groups"},
                     indent=2))


def _print_group(done, total, group) -> None:
    print(f"{done}/{total} region {group.region_index}: {len(group.components)} component(s), "
          f"unresolved={group.unresolved}, BIC gap={group.bic_gap:.1f}, "
          f"restart spread={group.restart_position_spread:.2f} px", flush=True)


def _group_metrics(group) -> dict:
    return {
        "region_index": group.region_index, "component_count": len(group.components),
        "bic": group.bic, "bic_gap": group.bic_gap,
        "restart_position_spread_pixels": group.restart_position_spread,
        "weighted_residual": group.weighted_residual,
        "explained_fraction": group.explained_fraction,
        "unresolved": group.unresolved, "warning": group.warning,
        "components": [{
            "t1_s": component.t1, "t2_s": component.t2, "volume": component.volume,
            "relative_to_largest": component.relative_to_largest,
            "confidence": component.confidence,
            "top_ions": component.mz[np.argsort(component.spectrum)[::-1][:5]].tolist(),
        } for component in group.components],
    }


if __name__ == "__main__":
    main()
