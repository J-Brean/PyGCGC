"""Known-area two-component injections into a real background region.

This is the first quantitative gate. It uses real background and real spectral
factors, but ideal Gaussian profiles, so passing it validates bookkeeping and
local fitting under background interference rather than realistic peak shape.
"""

from pathlib import Path
from types import SimpleNamespace
import json

import numpy as np
from scipy import sparse
from scipy.optimize import linear_sum_assignment

from pygcgc import (BlobCube, Run, candidate_cube, candidate_regions,
                    find_candidates, fit_overlap, fold, mcr_restarts,
                    remove_background)
from pygcgc.core.background import storage_threshold
from pygcgc.core.peaks import Frame


root = Path(__file__).resolve().parents[1]
data_dir = root / "real_data" / "blob_diagnosis"
z = np.load(data_dir / "segment.npz")
meta = json.loads((data_dir / "segment.json").read_text())
raw = Run(z["times"], z["mz"], sparse.csr_matrix(
    (z["values"], z["indices"], z["indptr"]), shape=tuple(z["shape"])),
    "injection benchmark")
folded = fold(raw, meta["period"], meta["phase"], t0=meta["fold_t0"])
clean = remove_background(raw, folded).run
floor = storage_threshold(raw) / 3.0
candidates = find_candidates(clean, folded, floor=floor)
regions = candidate_regions(candidates, (folded.n_rows + candidates.wrap_rows,
                                         folded.n_cols))

# Learn two distinct spectral templates from the stable real two-factor solution
# at the reported crowded location. Their identities remain unknown.
target_t1, target_t2 = 15 * 60, 1.2
region = min(regions, key=lambda r: min(
    np.hypot((candidates.t1[list(r.candidate_indices)] - target_t1) / 60.0,
             (candidates.t2[list(r.candidate_indices)] - target_t2) / 0.5)))
source_cube = candidate_cube(clean, folded, candidates, region, max_ions=60)
source = mcr_restarts(source_cube, 2, n_starts=5, max_iter=100).best
ranked = []
for j in range(2):
    ranked.extend(source_cube.ions[np.argsort(source.spectra[j])[::-1][:15]].tolist())
ions = np.asarray(list(dict.fromkeys(ranked)), dtype=int)
source_idx = np.asarray([int(np.flatnonzero(source_cube.ions == ion)[0]) for ion in ions])
spectra = source.spectra[:, source_idx].astype(float)
spectra /= np.maximum(spectra.sum(axis=1)[:, None], 1e-12)

# A patch in the reviewed background-comparison band. It had no candidate centres
# under the new default detector, but retains the real structured background.
frame = Frame(folded, candidates.wrap_rows)
full = frame.extend(folded.cube(clean, mz=ions, dtype=np.float64))
centre_col = int(np.argmin(abs((folded.t1 + 0.5 * folded.spec.period) - target_t1)))
centre_row = int(round(7.0 / folded.row_pitch))
row0, row1 = centre_row - 12, centre_row + 13
col0, col1 = centre_col - 2, centre_col + 3
background = full[row0:row1, col0:col1].copy()
mask = np.ones(background.shape[:2], dtype=bool)
rr, cc = np.indices(mask.shape, dtype=float)


def gaussian(row, col, sigma_row=3.0, sigma_col=0.8):
    shape = np.exp(-0.5 * ((rr-row)/sigma_row) ** 2 -
                   0.5 * ((cc-col)/sigma_col) ** 2)
    return shape / shape.sum()


results = []
base_area = 2.5e7
for separation_rows in (3, 6, 10):
    for minor_ratio in (1.0, 0.3, 0.1):
        centres = [(12.0-separation_rows/2, 1.7),
                   (12.0+separation_rows/2, 2.3)]
        truth = np.asarray([base_area, base_area * minor_ratio])
        injected = background.copy()
        for j, ((row, col), area) in enumerate(zip(centres, truth)):
            injected += gaussian(row, col)[:, :, None] * area * spectra[j][None, None, :]
        cube = BlobCube(injected, mask, ions, row0, col0, folded.row_pitch,
                        folded.spec.period, float(folded.t1[0]))
        hints = [SimpleNamespace(row=row, col=col) for row, col in centres]
        fit = fit_overlap(cube, 2, model="gaussian", include_background=True,
                          background_order=1, n_starts=8, seed=separation_rows,
                          apex_hint=hints)
        recovered_spectra = np.asarray([component.spectrum for component in fit.components])
        truth_unit = spectra / np.maximum(np.linalg.norm(spectra, axis=1)[:, None], 1e-12)
        recovered_unit = recovered_spectra / np.maximum(
            np.linalg.norm(recovered_spectra, axis=1)[:, None], 1e-12)
        cosine = truth_unit @ recovered_unit.T
        truth_idx, recovered_idx = linear_sum_assignment(-cosine)
        recovered = np.asarray([fit.components[j].volume for j in recovered_idx])
        errors = (recovered - truth[truth_idx]) / truth[truth_idx]
        results.append({
            "separation_rows": separation_rows,
            "separation_s": separation_rows * folded.row_pitch,
            "minor_ratio": minor_ratio,
            "truth": truth[truth_idx].tolist(),
            "recovered": recovered.tolist(),
            "relative_error": errors.tolist(),
            "spectral_cosine": cosine[truth_idx, recovered_idx].tolist(),
            "converged": fit.converged,
        })
        print(separation_rows, minor_ratio, np.round(errors, 3), flush=True)

summary = {
    "method": "ideal Gaussian profiles and two real MCR spectral factors injected into real background",
    "background_window": {"t1_columns": [col0, col1], "t2_rows": [row0, row1]},
    "ions": ions.tolist(), "base_area": base_area, "cases": results,
    "within_10_percent": int(sum(np.all(np.abs(case["relative_error"]) <= .1)
                                 for case in results)),
    "total_cases": len(results),
}
(data_dir / "injection_metrics.json").write_text(json.dumps(summary, indent=2))
