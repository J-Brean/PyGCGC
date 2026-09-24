"""How faint a peak does the automatic finder detect, and what does more sensitivity cost?

Weak isolated Gaussian peaks with real spectra are planted at known apex intensities
into the cleaned real segment. Detection rate is reported per intensity band for
several candidate-detector settings. As a false-positive proxy, the number of fitted
components in the reviewed background band (t1 14 to 16 min, t2 5 to 9 s, visually
background-dominated but not a certified blank) is reported for the same settings.

    python benchmarks/benchmark_sensitivity.py
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from benchmark_blob_finder import apex_to_area, load_background          # noqa: E402
from pygcgc import ComponentDetectionSettings, detect_components          # noqa: E402
from pygcgc.core.synthetic import (TruePeak, free_sites, inject_peaks,     # noqa: E402
                                   real_spectrum_library)

VARIANTS = {
    "default (seed 6, tic 2)": {},
    "seed 5": {"seed_snr": 5.0},
    "seed 4": {"seed_snr": 4.0},
    "seed 3": {"seed_snr": 3.0},
    "seed 4, tic 1": {"seed_snr": 4.0, "tic_prominence_snr": 1.0},
    "seed 3, tic 1": {"seed_snr": 3.0, "tic_prominence_snr": 1.0},
}
BANDS = ((1.5e4, 4e4), (4e4, 1e5), (1e5, 3e5))
SEEDS = (201, 202, 203)


def in_band(group, component, fold) -> bool:
    row, col = group.row0 + component.row0, group.col0 + component.col0
    t2, t1 = row * fold.row_pitch, fold.t1[0] + col * fold.spec.period
    return 5.0 <= t2 <= 9.0 and 14 * 60 <= t1 <= 16 * 60      # fold.t1 is absolute acquisition time


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    _raw, folded, clean, floor = load_background(root)
    base_settings = replace(ComponentDetectionSettings(), workers=4)
    candidates0, _r, _c = detect_components(clean, folded, "b", base_settings, floor=floor)[:3]
    library = real_spectrum_library(clean, folded, candidates0)
    results = {}
    for name, overrides in VARIANTS.items():
        settings = replace(base_settings, candidates=replace(base_settings.candidates, **overrides))
        started = time.perf_counter()
        _c, _r, background = detect_components(clean, folded, "b", settings, floor=floor)
        band = sum(in_band(g, c, folded) for g in background.groups for c in g.components)
        hits = {band_: [] for band_ in BANDS}
        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            sites = free_sites(candidates0, folded, rng, 21)
            peaks, truth = [], []
            for row, col in sites:
                apex = float(np.exp(rng.uniform(np.log(BANDS[0][0]), np.log(BANDS[-1][1]))))
                sr, sc = 3.0 * rng.uniform(0.85, 1.15), 0.8 * rng.uniform(0.9, 1.1)
                spectrum = library[int(rng.integers(len(library)))]
                peaks.append(TruePeak(row, col, apex_to_area(apex, sr, sc), spectrum, sr, sc))
                truth.append((row, col, apex))
            injected = inject_peaks(clean, folded, peaks, seed=seed)
            _c, _r, found = detect_components(injected, folded, f"s{seed}", settings, floor=floor)
            located = [(g.row0 + c.row0, g.col0 + c.col0) for g in found.groups for c in g.components]
            for row, col, apex in truth:
                hit = any(abs(r - row) <= 8 and abs(k - col) <= 2 for r, k in located)
                for lo, hi in BANDS:
                    if lo <= apex < hi:
                        hits[(lo, hi)].append(hit)
        results[name] = {"background_band_components": int(band),
                         "background_components": len(background.components),
                         "detection_rate": {f"{lo:.0e}-{hi:.0e}": float(np.mean(v)) if v else None
                                            for (lo, hi), v in hits.items()},
                         "n_per_band": {f"{lo:.0e}-{hi:.0e}": len(v) for (lo, hi), v in hits.items()},
                         "seconds": time.perf_counter() - started}
        r = results[name]
        print(f"{name:26s} band comps {band:3d}  total comps {r['background_components']:4d}  "
              + "  ".join(f"{k}: {v:.2f}" for k, v in r["detection_rate"].items())
              + f"  ({r['seconds']:.0f} s)", flush=True)
    path = root / "real_data" / "blob_diagnosis" / "sensitivity_benchmark.json"
    path.write_text(json.dumps(results, indent=1))
    print("wrote", path)


if __name__ == "__main__":
    main()
