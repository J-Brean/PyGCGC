"""Initial component-fit audit near the reported real-data overlap."""

from pathlib import Path
import json
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.optimize import linear_sum_assignment

from pygcgc import (Run, candidate_cube, candidate_regions, choose_k,
                    cross_validate_overlap, find_candidates, fit_overlap, fold, remove_background,
                    mcr_als, unique_ion_apexes)
from pygcgc.core.background import storage_threshold
from pygcgc.core.overlap import _shape


root = Path(__file__).resolve().parents[1]
data_dir = root / "real_data" / "blob_diagnosis"
z = np.load(data_dir / "segment.npz")
meta = json.loads((data_dir / "segment.json").read_text())
raw = Run(z["times"], z["mz"], sparse.csr_matrix(
    (z["values"], z["indices"], z["indptr"]), shape=tuple(z["shape"])),
    "component benchmark")
folded = fold(raw, meta["period"], meta["phase"], t0=meta["fold_t0"])
clean = remove_background(raw, folded).run
candidates = find_candidates(clean, folded, floor=storage_threshold(raw) / 3.0)
regions = candidate_regions(candidates, (folded.n_rows + candidates.wrap_rows, folded.n_cols))

# Pick the region nearest the location supplied by James. Distance scales match
# the project's normal retention-space convention.
target_t1, target_t2 = 15 * 60, 1.2
distance = []
for region in regions:
    idx = np.asarray(region.candidate_indices)
    distance.append(np.min(np.hypot((candidates.t1[idx] - target_t1) / 60.0,
                                    (candidates.t2[idx] - target_t2) / 0.5)))
region = regions[int(np.argmin(distance))]
cube = candidate_cube(clean, folded, candidates, region, max_ions=60)
apex_map = unique_ion_apexes(cube)

fits = []
started = time.perf_counter()
for model in ("gaussian", "emg"):
    for k in range(1, min(4, len(region.candidate_indices) + 2) + 1):
        fits.append(fit_overlap(cube, k, model=model, include_background=True,
                                n_starts=5, seed=100 + k))
elapsed = time.perf_counter() - started
best = min(fits, key=lambda fit: fit.bic)
cross_validation = {}
for k in (1, 2, 3):
    cv = cross_validate_overlap(cube, k, folds=3, block_rows=3,
                                model="gaussian", include_background=True,
                                background_order=1, n_starts=3, seed=200+k)
    cross_validation[str(k)] = {
        "fold_rmse": cv.fold_rmse.tolist(), "mean_rmse": cv.mean_rmse,
        "converged": cv.converged.tolist()}

# Free-form comparison. Repeated starts expose whether MCR has a stable answer;
# a low residual by itself is not enough evidence for a component count.
mcr_results = {}
mcr_saved = {}
data_rms = float(np.sqrt(np.mean(cube.data ** 2)))
for k in range(1, 5):
    runs = [mcr_als(cube, k, max_iter=100, seed=seed) for seed in range(5)]
    reference = runs[0].spectra / np.maximum(np.linalg.norm(runs[0].spectra, axis=1)[:, None], 1e-12)
    similarities = []
    for result in runs[1:]:
        spectra = result.spectra / np.maximum(np.linalg.norm(result.spectra, axis=1)[:, None], 1e-12)
        similarity = reference @ spectra.T
        rows_match, cols_match = linear_sum_assignment(-similarity)
        similarities.append(float(np.mean(similarity[rows_match, cols_match])))
    mcr_results[str(k)] = {
        "residuals": [result.resid for result in runs],
        "explained_variance": [1.0 - (result.resid / max(data_rms, 1e-12)) ** 2
                               for result in runs],
        "converged": [result.converged for result in runs],
        "mean_restart_spectral_similarity": float(np.mean(similarities)),
        "minimum_restart_spectral_similarity": float(np.min(similarities)),
    }
    best_run = min(runs, key=lambda result: result.resid)
    mcr_results[str(k)]["best_component_signal"] = [
        float(best_run.concentrations[:, :, j].sum() * best_run.spectra[j].sum())
        for j in range(k)]
    mcr_results[str(k)]["best_component_top_ions"] = [
        cube.ions[np.argsort(best_run.spectra[j])[::-1][:10]].tolist() for j in range(k)]
    if k in (2, 3):
        mcr_saved[k] = best_run

rr, cc = np.indices(cube.mask.shape, dtype=float)
reconstructed = np.zeros(cube.data.shape, dtype=float)
for component in best.components:
    params = [component.row0, component.col0, component.sigma_row, component.sigma_col]
    if best.model == "emg":
        params.append(component.tau)
    reconstructed += _shape(best.model, rr.ravel(), cc.ravel(), np.asarray(params)).reshape(
        cube.mask.shape)[:, :, None] * component.spectrum[None, None, :]
reconstructed += best.background_surface
observed_tic = cube.data.sum(axis=2)
fitted_tic = reconstructed.sum(axis=2)
residual_tic = observed_tic - fitted_tic

metrics = {
    "requested_t1_s": target_t1,
    "requested_t2_s": target_t2,
    "region": {"candidate_indices": list(region.candidate_indices),
               "row_bounds": [region.row0, region.row1],
               "column_bounds": [region.col0, region.col1],
               "shape": list(cube.data.shape), "ions": cube.ions.tolist()},
    "unique_ion_clusters": apex_map.n_components,
    "fit_seconds": elapsed,
    "fits": [{"k": fit.k, "model": fit.model, "bic": fit.bic,
              "rss": fit.rss, "weighted_rss": fit.weighted_rss,
              "n_parameters": fit.n_parameters, "converged": fit.converged}
             for fit in fits],
    "mcr": mcr_results,
    "spatial_cross_validation": cross_validation,
    "selected": {"k": best.k, "model": best.model, "bic": best.bic,
                 "component_t1_s": [x.t1 for x in best.components],
                 "component_t2_s": [x.t2 for x in best.components],
                 "component_volumes": [x.volume for x in best.components],
                 "residual_fraction": float(np.linalg.norm(residual_tic) /
                                            max(np.linalg.norm(observed_tic), 1e-12))},
}
(data_dir / "component_metrics.json").write_text(json.dumps(metrics, indent=2))
print(json.dumps(metrics["selected"], indent=2), flush=True)

t1 = cube.t1_of(np.arange(cube.data.shape[1])) / 60.0
t2 = cube.t2_of(np.arange(cube.data.shape[0]))
fig, axes = plt.subplots(1, 3, figsize=(8, 3), constrained_layout=True)
panels = ((observed_tic, "Observed"), (fitted_tic, f"Fitted: {best.k} {best.model}"),
          (residual_tic, "Residual"))
limit = np.quantile(np.abs(residual_tic), .99)
for ax, (image, title) in zip(axes, panels):
    if title == "Residual":
        mesh = ax.pcolormesh(t1, t2, image, cmap="coolwarm", vmin=-limit, vmax=limit,
                             shading="nearest")
    else:
        mesh = ax.pcolormesh(t1, t2, np.log10(1 + np.maximum(image, 0)), cmap="magma",
                             shading="nearest")
    fig.colorbar(mesh, ax=ax, shrink=.75)
    ax.set(title=title, xlabel="t1 (min)", ylabel="t2 (s)")
fig.savefig(data_dir / "component_fit.png", dpi=200)
plt.close(fig)

fig, axes = plt.subplots(2, 3, figsize=(8, 5), constrained_layout=True)
for row, k in enumerate((2, 3)):
    result = mcr_saved[k]
    for j in range(3):
        ax = axes[row, j]
        if j >= k:
            ax.axis("off")
            continue
        mesh = ax.pcolormesh(t1, t2, result.concentrations[:, :, j],
                             shading="nearest", cmap="viridis")
        fig.colorbar(mesh, ax=ax, shrink=.7)
        top = cube.ions[np.argsort(result.spectra[j])[::-1][:4]]
        ax.set(title=f"MCR k={k}, C{j+1}; top {','.join(map(str, top))}",
               xlabel="t1 (min)", ylabel="t2 (s)")
fig.savefig(data_dir / "mcr_components.png", dpi=200)
plt.close(fig)
