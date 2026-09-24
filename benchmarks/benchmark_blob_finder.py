"""Known-answer benchmark of the automatic overlapping-blob finder.

Real cleaned background and real spectra, planted peaks, the whole automatic
pipeline (no positions supplied). Each injection pass plants one case of every
type at non-interacting, candidate-free sites of the cached ten-minute segment.

    python benchmarks/benchmark_blob_finder.py --split tuning --label baseline \
        --set search_hidden_components=False

Splits use different random seeds and disjoint halves of the real spectrum
library (even/odd), so held-out cases share neither positions nor spectra with
tuning cases. Tune on ``tuning``; run ``heldout`` once, on frozen settings.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import time

import numpy as np
from scipy import sparse

from pygcgc import ComponentDetectionSettings, Run, detect_components, fold, remove_background
from pygcgc.core.background import storage_threshold
from pygcgc.core.synthetic import (TruePeak, TruthCase, count_tic_maxima, cosine, free_sites,
                                   inject_peaks, real_spectrum_library, score_case,
                                   spectral_pairs)

SEEDS = {"tuning": list(range(101, 111)), "heldout": list(range(901, 911))}
# (kind, ratio, separation in rows, separation in columns, spectral class)
PAIR_GRID = ([("pair", r, s, c, "distinct") for r in (1.0, 0.3, 0.1)
              for s, c in ((2.0, 0.0), (3.5, 0.5), (6.0, 1.0))] +
             [("pair", r, s, c, "similar") for r in (1.0, 0.3)
              for s, c in ((3.5, 0.5), (6.0, 1.0))])


def load_background(root: Path):
    data_dir = root / "real_data" / "blob_diagnosis"
    z = np.load(data_dir / "segment.npz")
    meta = json.loads((data_dir / "segment.json").read_text())
    raw = Run(z["times"], z["mz"], sparse.csr_matrix(
        (z["values"], z["indices"], z["indptr"]), shape=tuple(z["shape"])), "benchmark")
    folded = fold(raw, meta["period"], meta["phase"], t0=meta["fold_t0"])
    clean = remove_background(raw, folded).run
    return raw, folded, clean, storage_threshold(raw) / 3.0


def apex_to_area(apex_tic: float, sigma_row: float, sigma_col: float) -> float:
    return apex_tic * 2.0 * np.pi * sigma_row * sigma_col


def build_cases(rng, sites, library, half: int) -> list[TruthCase]:
    """One case of every type for one injection pass; ``half`` picks the library half."""
    pool = np.flatnonzero(np.arange(len(library)) % 2 == half)
    lib = library[pool]
    distinct = spectral_pairs(lib, 0.0, 0.5, rng, 200)
    similar = spectral_pairs(lib, 0.6, 0.85, rng, 200)
    specs = ([("single", 0, 0, 0, "-")] * 2 + [("single_broad", 0, 0, 0, "-")] * 2 +
             [("single_tail", 0, 0, 0, "-")] * 2 + PAIR_GRID + [("pair_tail", 0.3, 3.5, 0.5, "distinct")] * 2)
    order = rng.permutation(len(specs))
    cases, counters = [], {"distinct": 0, "similar": 0}
    for site, index in zip(sites, order):
        kind, ratio, sep_rows, sep_cols, spec_class = specs[index]
        row, col = site
        apex = float(np.exp(rng.uniform(np.log(3e5), np.log(1.5e6))))
        sr, sc = 3.0 * rng.uniform(0.85, 1.15), 0.8 * rng.uniform(0.9, 1.1)
        tags = {"apex_tic": apex}
        if kind.startswith("single"):
            tau_row = tau_col = 0.0
            if kind == "single_broad":
                sr, sc = rng.uniform(5.0, 7.0), rng.uniform(1.4, 1.8)
            if kind == "single_tail":
                sr = 2.5
                tau_row, tau_col = rng.choice([(3.0, 0.0), (6.0, 0.0), (0.0, 1.0)])
            spectrum = lib[int(rng.integers(len(lib)))]
            peaks = [TruePeak(row, col, apex_to_area(apex, sr, sc), spectrum, sr, sc,
                              float(tau_row), float(tau_col))]
        else:
            pairs = distinct if spec_class == "distinct" else similar
            i, j = pairs[counters[spec_class] % len(pairs)]
            counters[spec_class] += 1
            first_major = bool(rng.integers(2))
            ratios = (1.0, ratio) if first_major else (ratio, 1.0)
            tail = 4.0 if kind == "pair_tail" else 0.0
            peaks = []
            for n, (spectrum, r) in enumerate(((lib[i], ratios[0]), (lib[j], ratios[1]))):
                s_row = sr * rng.uniform(0.95, 1.05)
                peaks.append(TruePeak(
                    row + n * sep_rows, col + n * sep_cols,
                    apex_to_area(apex * r, s_row, sc), spectrum, s_row, sc,
                    tail if (kind == "pair_tail" and r == 1.0) else 0.0))
            tags.update(ratio=ratio, sep_rows=sep_rows, sep_cols=sep_cols,
                        spectral_class=spec_class,
                        spectral_cosine=cosine(lib[i], lib[j]))
        cases.append(TruthCase(f"{kind}@{row:.0f},{col:.0f}", kind, peaks, tags,
                               count_tic_maxima(peaks)))
    return cases


def summarise(results: list[dict]) -> dict:
    """Aggregate per case kind and, for pairs, per abundance ratio."""
    def stats(rows):
        peaks = [p for r in rows for p in r["peaks"]]
        area = np.abs([p["area_error"] for p in peaks]) if peaks else np.array([np.nan])
        return {
            "cases": len(rows),
            "count_correct": float(np.mean([r["count_correct"] for r in rows])),
            "extra_component_rate": float(np.mean([r["extra"] > 0 for r in rows])),
            "material_extra_rate": float(np.mean([r["extra_material"] > 0 for r in rows])),
            "peak_recall": float(sum(r["n_matched"] for r in rows) / max(sum(r["n_true"] for r in rows), 1)),
            "median_centre_error_s": float(np.median([np.hypot(p["row_error_s"], p["col_error_s"] / 10)
                                                      for p in peaks])) if peaks else np.nan,
            "median_cosine": float(np.median([p["cosine"] for p in peaks])) if peaks else np.nan,
            "area_within_10pc": float(np.mean(area <= 0.10)),
            "area_within_25pc": float(np.mean(area <= 0.25)),
            "median_abs_area_error": float(np.median(area)),
        }
    out = {}
    for kind in ("single", "single_broad", "single_tail", "pair_tail"):
        rows = [r for r in results if r["kind"] == kind]
        if rows:
            out[kind] = stats(rows)
    for ratio in (1.0, 0.3, 0.1):
        rows = [r for r in results if r["kind"] == "pair" and r["tags"]["ratio"] == ratio]
        if rows:
            out[f"pair_ratio_{ratio:g}"] = stats(rows)
            hidden = [r for r in rows if r["tic_maxima"] == 1]
            if hidden:
                out[f"pair_ratio_{ratio:g}_hidden_tic_max"] = stats(hidden)
    for cls in ("distinct", "similar"):
        rows = [r for r in results if r["kind"] == "pair" and r["tags"]["spectral_class"] == cls]
        if rows:
            out[f"pair_{cls}_spectra"] = stats(rows)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split", choices=["tuning", "heldout"], default="tuning")
    parser.add_argument("--label", default="run")
    parser.add_argument("--passes", type=int, default=None, help="injection passes (default: all seeds)")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--set", nargs="*", default=[], metavar="FIELD=VALUE",
                        help="override ComponentDetectionSettings fields")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    raw, folded, clean, floor = load_background(root)
    overrides = {}
    for item in args.set:
        key, value = item.split("=", 1)
        try:
            overrides[key] = json.loads(value.lower() if value in ("True", "False") else value)
        except json.JSONDecodeError:
            overrides[key] = value                    # a bare string such as edge_policy=drop
    settings = replace(ComponentDetectionSettings(workers=args.workers), **overrides)

    # Baseline: what the pipeline finds in the untouched run.
    started = time.perf_counter()
    candidates, _regions, base = detect_components(clean, folded, "base", settings, floor=floor)
    base_seconds = time.perf_counter() - started
    ignore = [(g.row0 + c.row0, g.col0 + c.col0) for g in base.groups for c in g.components]
    library = real_spectrum_library(clean, folded, candidates)
    print(f"background run: {len(candidates)} candidates, {len(ignore)} components, "
          f"{base_seconds:.1f} s; spectrum library {len(library)}", flush=True)

    results, timings = [], []
    seeds = SEEDS[args.split][:args.passes]
    for seed in seeds:
        rng = np.random.default_rng(seed)
        sites = free_sites(candidates, folded, rng, 21)
        cases = build_cases(rng, sites, library, half=0 if args.split == "tuning" else 1)
        planted = [p for case in cases for p in case.peaks]
        injected = inject_peaks(clean, folded, planted, seed=seed)
        started = time.perf_counter()
        _c, regions, found = detect_components(injected, folded, f"seed{seed}", settings, floor=floor)
        elapsed = time.perf_counter() - started
        timings.append({"seed": seed, "seconds": elapsed, "windows": len(regions),
                        "components": len(found.components)})
        rows = [score_case(case, found, clean.mz, ignore=ignore, pitch=folded.row_pitch,
                           period=folded.spec.period) for case in cases]
        # Injecting strong peaks shifts the detector's global noise statistics, so
        # baseline components move by a pixel or two and a few appear or vanish.
        # Positional matching to the baseline would overcount, so report the net
        # change in the number of components outside every case neighbourhood.
        centres = [np.mean([p.apex() for p in case.peaks], axis=0) for case in cases]
        off_target = sum(
            not any(abs(g.row0 + c.row0 - a) <= 15 and abs(g.col0 + c.col0 - b) <= 3
                    for a, b in centres)
            for g in found.groups for c in g.components)
        timings[-1]["off_target_components"] = int(off_target)
        timings[-1]["off_target_net_of_baseline"] = int(off_target - len(ignore))
        results.extend(rows)
        print(f"seed {seed}: {len(cases)} cases, {elapsed:.1f} s, {len(regions)} windows, "
              f"off-target components {off_target} vs baseline {len(ignore)}", flush=True)

    summary = summarise(results)
    header = f"{'group':32s} {'n':>3s} {'count':>6s} {'extra':>6s} {'xmat':>5s} {'recall':>6s} {'cos':>5s} {'<10%':>5s} {'<25%':>5s} {'|dA|':>6s} {'dt s':>6s}"
    print(header)
    for name, s in summary.items():
        print(f"{name:32s} {s['cases']:3d} {s['count_correct']:6.2f} {s['extra_component_rate']:6.2f} {s['material_extra_rate']:5.2f} "
              f"{s['peak_recall']:6.2f} {s['median_cosine']:5.2f} {s['area_within_10pc']:5.2f} "
              f"{s['area_within_25pc']:5.2f} {s['median_abs_area_error']:6.2f} {s['median_centre_error_s']:6.3f}")
    out = {"label": args.label, "split": args.split, "settings": settings.to_dict(),
           "overrides": overrides, "timings": timings,
           "background_run_seconds": base_seconds, "summary": summary, "cases": results}
    path = root / "real_data" / "blob_diagnosis" / f"blob_benchmark_{args.split}_{args.label}.json"
    path.write_text(json.dumps(out, indent=1, default=float))
    print("wrote", path)


if __name__ == "__main__":
    main()
