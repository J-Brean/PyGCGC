"""Recovery test using chromatographic profiles measured in the real run.

An "injection" here is a controlled software experiment: copy a shape measured
in the run and two mass-spectral patterns into a quiet real-data patch at known
amounts, then ask the overlap fitter to recover those amounts.  The ten-minute
crop contains no automatically defensible isolated source, so the selected
shape is taken from the confirmed n-alkane ladder in the full acquisition.  The
alkane identity makes it a much better shape reference, although it is still an
observed sample peak rather than a separately acquired pure standard.
"""

from pathlib import Path
from types import SimpleNamespace
import json
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi, sparse
from scipy.optimize import linear_sum_assignment, nnls

from pygcgc import (BlobCube, Run, candidate_cube, candidate_regions,
                    find_candidates, fit_overlap, fit_template_overlap, fold,
                    mcr_restarts, remove_background)
from pygcgc.core.background import storage_threshold
from pygcgc.core.hydrocarbons import HydrocarbonSettings, find_alkane_series
from pygcgc.core.io.markes_lsc import read_lsc
from pygcgc.core.overlap import MODELS
from pygcgc.core.peaks import Frame


root = Path(__file__).resolve().parents[1]
data_dir = root / "real_data" / "blob_diagnosis"
z = np.load(data_dir / "segment.npz")
meta = json.loads((data_dir / "segment.json").read_text())
raw = Run(z["times"], z["mz"], sparse.csr_matrix(
    (z["values"], z["indices"], z["indptr"]), shape=tuple(z["shape"])),
    "measured-profile benchmark")
folded = fold(raw, meta["period"], meta["phase"], t0=meta["fold_t0"])
clean = remove_background(raw, folded).run
floor = storage_threshold(raw) / 3.0
candidates = find_candidates(clean, folded, floor=floor)
regions = candidate_regions(candidates, (folded.n_rows + candidates.wrap_rows,
                                         folded.n_cols))


def unit_rows(values):
    values = np.asarray(values, dtype=float)
    return values / np.maximum(np.linalg.norm(values, axis=1)[:, None], 1e-12)


# Use n-C16 from the independently established alkane ladder as a trustworthy
# real chromatographic shape.  C17 is deliberately avoided because it co-elutes
# with pristane in this sample.  Loading the existing full-run cache takes a few
# seconds; without it the native file is decoded once and cached in the system
# temporary directory, matching examples/validate_real_data.py.
cache = Path(tempfile.gettempdir()) / "pygcgc_real70.npz"
if cache.exists():
    full_z = np.load(cache)
    full_run = Run(full_z["times"], full_z["mz"], sparse.csr_matrix(
        (full_z["values"], full_z["indices"], full_z["indptr"]),
        shape=tuple(full_z["shape"])), "full real 70 eV")
else:
    full_run = read_lsc(next((root / "real_data").glob("*_70eV.lsc")))
    matrix = full_run.data
    np.savez(cache, times=full_run.times, mz=full_run.mz, values=matrix.data,
             indices=matrix.indices, indptr=matrix.indptr, shape=np.array(matrix.shape))
full_fold = fold(full_run, meta["period"], meta["phase"])
ladder = find_alkane_series(
    full_run, full_fold, [(13, 23.45 * 60, 0.74), (18, 42.95 * 60, 0.92)],
    HydrocarbonSettings())
alkane = ladder.hits[16]
source_ions = np.asarray([57, 71, 85, 226])  # C16 molecular ion is m/z 226.
source_data = np.stack([full_fold.image(full_run.xic(int(ion)))
                        for ion in source_ions], axis=2)
source_col = int(np.argmin(abs((full_fold.t1 + 0.5 * full_fold.spec.period) - alkane.t1)))
source_row = int(round(alkane.t2 / full_fold.row_pitch))
source_window = source_data[source_row-12:source_row+13,
                            source_col-2:source_col+3].astype(float)
rim = np.zeros(source_window.shape[:2], dtype=bool)
rim[[0, -1], :] = True
rim[:, [0, -1]] = True
baseline = np.median(source_window[rim], axis=0)
profile = np.maximum(source_window - baseline[None, None, :], 0.0).sum(axis=2)
profile /= max(profile.sum(), 1e-12)
profile_source = {
    "name": "n-C16", "carbon_number": 16, "t1_s": float(alkane.t1),
    "t2_s": float(alkane.t2), "ions": source_ions.tolist(), "profile": profile,
}

# The two spectra still come from the stable two-factor solution in the
# difficult 15 min / 1.2 s region.  This changes only the chromatographic shape
# relative to benchmark_injection_recovery.py, making the comparison interpretable.
target_t1, target_t2 = 15 * 60, 1.2
target_region = min(regions, key=lambda region: min(
    np.hypot((candidates.t1[list(region.candidate_indices)] - target_t1) / 60.0,
             (candidates.t2[list(region.candidate_indices)] - target_t2) / 0.5)))
spectral_cube = candidate_cube(clean, folded, candidates, target_region, max_ions=60)
spectral_fit = mcr_restarts(spectral_cube, 2, n_starts=5, max_iter=100).best
ranked_ions = []
for j in range(2):
    ranked_ions.extend(spectral_cube.ions[
        np.argsort(spectral_fit.spectra[j])[::-1][:15]].tolist())
all_ions = np.asarray(list(dict.fromkeys(ranked_ions)), dtype=int)
source_idx = np.asarray([int(np.flatnonzero(spectral_cube.ions == ion)[0])
                         for ion in all_ions])
spectra = spectral_fit.spectra[:, source_idx].astype(float)
spectra /= np.maximum(spectra.sum(axis=1)[:, None], 1e-12)

# Reuse a quiet patch from the reviewed high-t2 comparison band.  It contains
# actual structured detector/background signal rather than generated noise.
frame = Frame(folded, candidates.wrap_rows)
full = frame.extend(folded.cube(clean, mz=all_ions, dtype=np.float64))
target_t1 = 15 * 60
centre_col = int(np.argmin(abs((folded.t1 + 0.5 * folded.spec.period) - target_t1)))
centre_row = int(round(7.0 / folded.row_pitch))
row0, row1 = centre_row - 22, centre_row + 23
col0, col1 = centre_col - 4, centre_col + 5
background = full[row0:row1, col0:col1].copy()
# A background estimate available to a real analysis: corresponding pixels in
# the immediately neighbouring five-column windows, averaged.  The actual
# injected patch is not used to construct this estimate.
flank_background = 0.5 * (full[row0:row1, col0-9:col1-9] +
                          full[row0:row1, col0+9:col1+9])
mask = np.ones(background.shape[:2], dtype=bool)


def place(profile, row, col):
    """Centre a measured profile on a larger canvas without changing its area."""
    canvas = np.zeros(mask.shape, dtype=float)
    h, w = profile.shape
    r0 = (canvas.shape[0] - h) // 2
    c0 = (canvas.shape[1] - w) // 2
    canvas[r0:r0+h, c0:c0+w] = profile
    centre = ndi.center_of_mass(canvas)
    moved = ndi.shift(canvas, (row-centre[0], col-centre[1]), order=1,
                      mode="constant", cval=0.0, prefilter=False)
    return moved / max(moved.sum(), 1e-12)


def fitted_errors(fit, truth):
    fitted_spectra = np.asarray([component.spectrum for component in fit.components])
    similarity = unit_rows(spectra) @ unit_rows(fitted_spectra).T
    truth_idx, fitted_idx = linear_sum_assignment(-similarity)
    amounts = np.asarray([fit.components[j].volume for j in fitted_idx])
    return ((amounts - truth[truth_idx]) / truth[truth_idx]).tolist()


results = []
model_screen_case = None
base_area = 2.5e7
for separation_rows in (3, 6, 10):
    for minor_ratio in (1.0, 0.3, 0.1):
        centres = [(22.0-separation_rows/2, 3.7),
                   (22.0+separation_rows/2, 4.3)]
        profiles = [place(profile_source["profile"], row, col) for row, col in centres]
        truth = np.asarray([base_area, base_area * minor_ratio])
        signal_only = np.zeros_like(background)
        for profile, spectrum, area in zip(profiles, spectra, truth):
            signal_only += profile[:, :, None] * area * spectrum[None, None, :]
        observed = background + signal_only
        cube = BlobCube(observed, mask, all_ions, row0, col0, folded.row_pitch,
                        folded.spec.period, float(folded.t1[0]))
        hints = [SimpleNamespace(row=row, col=col) for row, col in centres]
        fit = fit_overlap(cube, 2, model="gaussian", include_background=True,
                          background_order=1, n_starts=8, seed=separation_rows,
                          apex_hint=hints)
        recovered_spectra = np.asarray([component.spectrum for component in fit.components])
        cosine = unit_rows(spectra) @ unit_rows(recovered_spectra).T
        truth_idx, recovered_idx = linear_sum_assignment(-cosine)
        recovered = np.asarray([fit.components[j].volume for j in recovered_idx])
        errors = (recovered - truth[truth_idx]) / truth[truth_idx]
        shifted_template_fit = fit_template_overlap(
            cube, profile_source["profile"], 2, apex_hint=hints,
            include_background=True, background_order=1, n_starts=5,
            seed=separation_rows)
        shifted_template_errors = fitted_errors(shifted_template_fit, truth)

        # Control: use the exact measured shape at the known centres, while
        # still estimating both component spectra and a sloping background.
        # If this succeeds where the Gaussian fit fails, peak-shape mismatch is
        # the cause rather than a bookkeeping or spectral-separation failure.
        rows, cols = np.indices(mask.shape, dtype=float)
        rn = rows.ravel() / max(mask.shape[0] - 1, 1)
        cn = cols.ravel() / max(mask.shape[1] - 1, 1)
        design = np.column_stack((profiles[0].ravel(), profiles[1].ravel(),
                                  (1-rn)*(1-cn), (1-rn)*cn,
                                  rn*(1-cn), rn*cn))
        component_design = design[:, :2]
        flat = observed.reshape(-1, observed.shape[2])
        weights = 1.0 / np.sqrt(np.clip(flat, 1.0, None))
        coefficients = np.column_stack([
            nnls(design * weights[:, ion, None], flat[:, ion] * weights[:, ion])[0]
            for ion in range(flat.shape[1])])
        template_spectra = coefficients[:2]
        template_cosine = unit_rows(spectra) @ unit_rows(template_spectra).T
        template_truth_idx, template_recovered_idx = linear_sum_assignment(-template_cosine)
        template_recovered = coefficients[:2].sum(axis=1)[template_recovered_idx]
        template_errors = ((template_recovered - truth[template_truth_idx]) /
                           truth[template_truth_idx])

        # Estimate the smooth background only where neither proposed component
        # has appreciable template signal.  Holding it fixed afterwards prevents
        # a broad background term from taking intensity away from the peaks.
        peak_exclusion = ((profiles[0] > 0.01 * profiles[0].max()) |
                          (profiles[1] > 0.01 * profiles[1].max())).ravel()
        background_design = design[:, 2:]
        outside = ~peak_exclusion
        masked_background_coefficients = np.column_stack([
            nnls(background_design[outside] * weights[outside, ion, None],
                 flat[outside, ion] * weights[outside, ion])[0]
            for ion in range(flat.shape[1])])
        masked_background = (background_design @ masked_background_coefficients).reshape(
            observed.shape)
        masked_corrected = observed - masked_background
        masked_cube = BlobCube(masked_corrected, mask, all_ions, row0, col0,
                               folded.row_pitch, folded.spec.period,
                               float(folded.t1[0]))
        masked_fit = fit_template_overlap(
            masked_cube, profile_source["profile"], 2, apex_hint=hints,
            n_starts=5, seed=separation_rows)
        masked_errors = fitted_errors(masked_fit, truth)

        flank_corrected = observed - flank_background
        flank_cube = BlobCube(flank_corrected, mask, all_ions, row0, col0,
                              folded.row_pitch, folded.spec.period,
                              float(folded.t1[0]))
        flank_fit = fit_overlap(flank_cube, 2, model="gaussian", n_starts=8,
                                seed=separation_rows, apex_hint=hints)
        flank_spectra = np.asarray([component.spectrum for component in flank_fit.components])
        flank_cosine = unit_rows(spectra) @ unit_rows(flank_spectra).T
        flank_truth_idx, flank_recovered_idx = linear_sum_assignment(-flank_cosine)
        flank_recovered = np.asarray(
            [flank_fit.components[j].volume for j in flank_recovered_idx])
        flank_errors = ((flank_recovered - truth[flank_truth_idx]) /
                        truth[flank_truth_idx])

        flank_flat = flank_corrected.reshape(-1, flank_corrected.shape[2])
        # Noise is set by the observed counts before subtraction, not by signed
        # residuals after the neighbouring background estimate is removed.
        flank_weights = 1.0 / np.sqrt(np.clip(flat, 1.0, None))
        flank_coefficients = np.column_stack([
            nnls(component_design * flank_weights[:, ion, None],
                 flank_flat[:, ion] * flank_weights[:, ion])[0]
            for ion in range(flank_flat.shape[1])])
        flank_template_spectra = flank_coefficients[:2]
        flank_template_cosine = unit_rows(spectra) @ unit_rows(flank_template_spectra).T
        flank_template_truth_idx, flank_template_recovered_idx = linear_sum_assignment(
            -flank_template_cosine)
        flank_template_recovered = flank_coefficients.sum(axis=1)[flank_template_recovered_idx]
        flank_template_errors = ((flank_template_recovered - truth[flank_template_truth_idx]) /
                                 truth[flank_template_truth_idx])

        clean_cube = BlobCube(signal_only, mask, all_ions, row0, col0, folded.row_pitch,
                              folded.spec.period, float(folded.t1[0]))
        clean_fit = fit_overlap(clean_cube, 2, model="gaussian", n_starts=8,
                                seed=separation_rows, apex_hint=hints)
        clean_spectra = np.asarray([component.spectrum for component in clean_fit.components])
        clean_cosine = unit_rows(spectra) @ unit_rows(clean_spectra).T
        clean_truth_idx, clean_recovered_idx = linear_sum_assignment(-clean_cosine)
        clean_recovered = np.asarray(
            [clean_fit.components[j].volume for j in clean_recovered_idx])
        clean_errors = ((clean_recovered - truth[clean_truth_idx]) /
                        truth[clean_truth_idx])
        clean_shifted_template_fit = fit_template_overlap(
            clean_cube, profile_source["profile"], 2, apex_hint=hints,
            n_starts=5, seed=separation_rows)
        clean_shifted_template_errors = fitted_errors(clean_shifted_template_fit, truth)

        clean_flat = signal_only.reshape(-1, signal_only.shape[2])
        clean_weights = 1.0 / np.sqrt(np.clip(clean_flat, 1.0, None))
        clean_coefficients = np.column_stack([
            nnls(component_design * clean_weights[:, ion, None],
                 clean_flat[:, ion] * clean_weights[:, ion])[0]
            for ion in range(clean_flat.shape[1])])
        clean_template_spectra = clean_coefficients[:2]
        clean_template_cosine = unit_rows(spectra) @ unit_rows(clean_template_spectra).T
        clean_template_truth_idx, clean_template_recovered_idx = linear_sum_assignment(
            -clean_template_cosine)
        clean_template_recovered = clean_coefficients.sum(axis=1)[clean_template_recovered_idx]
        clean_template_errors = ((clean_template_recovered - truth[clean_template_truth_idx]) /
                                 truth[clean_template_truth_idx])
        results.append({
            "separation_rows": separation_rows,
            "separation_s": separation_rows * folded.row_pitch,
            "minor_ratio": minor_ratio, "truth": truth[truth_idx].tolist(),
            "recovered": recovered.tolist(), "relative_error": errors.tolist(),
            "spectral_cosine": cosine[truth_idx, recovered_idx].tolist(),
            "shifted_template_relative_error": shifted_template_errors,
            "masked_background_template_relative_error": masked_errors,
            "template_recovered": template_recovered.tolist(),
            "template_relative_error": template_errors.tolist(),
            "template_spectral_cosine": template_cosine[
                template_truth_idx, template_recovered_idx].tolist(),
            "flank_gaussian_relative_error": flank_errors.tolist(),
            "flank_template_relative_error": flank_template_errors.tolist(),
            "clean_gaussian_relative_error": clean_errors.tolist(),
            "clean_shifted_template_relative_error": clean_shifted_template_errors,
            "clean_template_relative_error": clean_template_errors.tolist(),
            "converged": fit.converged,
        })
        if separation_rows == 6 and minor_ratio == 0.3:
            model_screen_case = (signal_only.copy(), observed.copy(), truth.copy(), hints)
        print(separation_rows, minor_ratio, "Gaussian", np.round(errors, 3),
              "shifted template", np.round(shifted_template_errors, 3),
              "masked background", np.round(masked_errors, 3),
              "measured template", np.round(template_errors, 3),
              "flank Gaussian", np.round(flank_errors, 3),
              "flank template", np.round(flank_template_errors, 3),
              "clean", np.round(clean_errors, 3),
              "clean template", np.round(clean_template_errors, 3), flush=True)

# Check whether one of the shape families already present in PyGCGC solves the
# trusted-profile mismatch.  One representative case limits runtime; the winner
# must still pass the full grid before it can become a default.
screen_signal, screen_observed, screen_truth, screen_hints = model_screen_case
model_screen = {}
for model in MODELS:
    clean_model_cube = BlobCube(screen_signal, mask, all_ions, row0, col0,
                                folded.row_pitch, folded.spec.period,
                                float(folded.t1[0]))
    real_model_cube = BlobCube(screen_observed, mask, all_ions, row0, col0,
                               folded.row_pitch, folded.spec.period,
                               float(folded.t1[0]))
    clean_model_fit = fit_overlap(clean_model_cube, 2, model=model, n_starts=3,
                                  seed=17, apex_hint=screen_hints)
    real_model_fit = fit_overlap(real_model_cube, 2, model=model, include_background=True,
                                 background_order=1, n_starts=3, seed=17,
                                 apex_hint=screen_hints)
    model_screen[model] = {
        "clean_relative_error": fitted_errors(clean_model_fit, screen_truth),
        "real_background_relative_error": fitted_errors(real_model_fit, screen_truth),
        "clean_weighted_rss": clean_model_fit.weighted_rss,
        "real_background_weighted_rss": real_model_fit.weighted_rss,
    }

# Repeat one representative recovery case on the quietest candidate-free window
# available near several t2 values.  This tests whether the result is peculiar
# to the original high-t2 background patch.
background_screen = []
for requested_t2 in (1.2, 3.0, 5.0, 7.0, 9.0):
    requested_row = int(round(requested_t2 / folded.row_pitch))
    options = []
    for candidate_col in range(5, folded.n_cols - 5):
        occupied = ((abs(candidates.row - requested_row) <= 22) &
                    (abs(candidates.col - candidate_col) <= 4)).any()
        if occupied:
            continue
        patch = full[requested_row-22:requested_row+23,
                     candidate_col-4:candidate_col+5]
        if patch.shape == screen_signal.shape:
            options.append((float(patch.sum()), candidate_col, patch))
    if not options:
        background_screen.append({"requested_t2_s": requested_t2,
                                  "status": "no candidate-free window"})
        continue
    _, candidate_col, patch = min(options, key=lambda item: item[0])
    trial_cube = BlobCube(patch + screen_signal, mask, all_ions,
                          requested_row-22, candidate_col-4, folded.row_pitch,
                          folded.spec.period, float(folded.t1[0]))
    trial_fit = fit_template_overlap(
        trial_cube, profile_source["profile"], 2, apex_hint=screen_hints,
        include_background=True, background_order=1, n_starts=5, seed=31)
    background_screen.append({
        "requested_t2_s": requested_t2, "status": "fitted",
        "t1_s": float(folded.t1[candidate_col] + 0.5 * folded.spec.period),
        "background_sum": float(patch.sum()),
        "relative_error": fitted_errors(trial_fit, screen_truth),
    })


summary = {
    "method": "measured n-C16 profile, with two real-overlap MCR spectra, copied into real background",
    "source_warning": "n-C16 is confirmed from the alkane ladder but is an observed sample peak, not a separately acquired pure standard",
    "profile_source": {k: v for k, v in profile_source.items()
                       if k not in ("profile", "ions", "spectrum")},
    "source_spectral_cosine": float((unit_rows(spectra) @ unit_rows(spectra).T)[0, 1]),
    "ions": all_ions.tolist(), "base_area": base_area, "cases": results,
    "within_10_percent": int(sum(np.all(np.abs(case["relative_error"]) <= .1)
                                 for case in results)),
    "shifted_template_within_10_percent": int(sum(
        np.all(np.abs(case["shifted_template_relative_error"]) <= .1)
        for case in results)),
    "masked_background_template_within_10_percent": int(sum(
        np.all(np.abs(case["masked_background_template_relative_error"]) <= .1)
        for case in results)),
    "template_within_10_percent": int(sum(
        np.all(np.abs(case["template_relative_error"]) <= .1) for case in results)),
    "flank_gaussian_within_10_percent": int(sum(
        np.all(np.abs(case["flank_gaussian_relative_error"]) <= .1) for case in results)),
    "flank_template_within_10_percent": int(sum(
        np.all(np.abs(case["flank_template_relative_error"]) <= .1) for case in results)),
    "clean_gaussian_within_10_percent": int(sum(
        np.all(np.abs(case["clean_gaussian_relative_error"]) <= .1) for case in results)),
    "clean_shifted_template_within_10_percent": int(sum(
        np.all(np.abs(case["clean_shifted_template_relative_error"]) <= .1)
        for case in results)),
    "clean_template_within_10_percent": int(sum(
        np.all(np.abs(case["clean_template_relative_error"]) <= .1) for case in results)),
    "total_cases": len(results),
    "model_screen_case": {"separation_s": 0.12, "minor_ratio": 0.3},
    "model_screen": model_screen,
    "background_window_screen": background_screen,
}
(data_dir / "measured_profile_metrics.json").write_text(json.dumps(summary, indent=2))

matrices = [np.asarray([[max(abs(x) for x in case[key])
                         for case in results if case["minor_ratio"] == ratio]
                        for ratio in (1.0, 0.3, 0.1)])
            for key in ("relative_error", "shifted_template_relative_error",
                        "masked_background_template_relative_error",
                        "template_relative_error")]
limit = max(1.0, *(float(matrix.max()) for matrix in matrices))
fig, axes = plt.subplots(1, 4, figsize=(14.5, 3.2), constrained_layout=True)
for ax, matrix, title in zip(axes, matrices,
                             ("Gaussian fit", "Movable measured template",
                              "Background fitted outside peaks",
                              "Known-position shape control")):
    image = ax.imshow(matrix, vmin=0, vmax=limit, cmap="magma", aspect="auto")
    for r in range(matrix.shape[0]):
        for c in range(matrix.shape[1]):
            ax.text(c, r, f"{100*matrix[r,c]:.0f}%", ha="center", va="center",
                    color="white" if matrix[r,c] > .35 * limit else "black")
    ax.set(xticks=range(3), xticklabels=("0.06", "0.12", "0.20"),
           yticks=range(3), yticklabels=("100%", "30%", "10%"),
           xlabel="Peak separation in t2 (s)", ylabel="Smaller peak relative size",
           title=title)
fig.colorbar(image, ax=axes, label="largest absolute relative error")
fig.savefig(data_dir / "measured_profile_recovery.png", dpi=200)
plt.close(fig)
