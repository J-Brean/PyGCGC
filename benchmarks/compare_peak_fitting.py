"""Compare eight fitting techniques on fixed calibration windows and held-out pixels.

Raw acquisitions are read only. The manifest fixes window geometry, ions and
initial component counts, so this tests fitting rather than chemical identity
or candidate detection. Three interleaved spatial-block folds test prediction.
The original implementation is loaded from the recorded local Git revision.
Run from the repository root; results and extracted cubes belong in real_data/.
"""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np
from scipy import ndimage as ndi
from scipy.optimize import least_squares
from scipy.optimize import linear_sum_assignment
from threadpoolctl import threadpool_limits

from pygcgc import fold, read_run, remove_background
from pygcgc.core.candidates import CandidateRegion, _fold_region
from pygcgc.core.io import crop_run
from pygcgc.core.components import ComponentDetectionSettings
from pygcgc.core.modelselect import select_components
from pygcgc.core.overlap import (BlobCube, _emg_1d, _solve_linear, fit_overlap,
                                  predict_overlap, _shape)

MANIFEST = Path(__file__).parent / 'data' / 'calibration_fitting_windows.json'
METHODS = ('original', 'gaussian', 'cell', 'emg_t2', 'tilt', 'emg_t2_tilt',
           'modulated', 'modulated_tail')


def load_original(ref: str, module: str):
    """Load the recorded scientific code without changing the working tree."""
    source = subprocess.check_output(['git', 'show', f'{ref}:src/pygcgc/core/{module}.py'])
    name = f'pygcgc.core._benchmark_original_{module}'
    spec = importlib.util.spec_from_loader(name, loader=None)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    exec(compile(source, name, 'exec'), loaded.__dict__)
    return loaded


def modulated_fit(cube, hints, tail=False, max_nfev=60):
    """Experimental shared-spectrum peak with free modulation amplitudes.

    A common second-dimension width and linear apex drift link the slices.
    Each peak is restricted to seven columns around its seed. One amplitude is
    fixed to remove spectral scale ambiguity. This is a benchmark alternative,
    not a production model or a claim that a factor is one chemical compound.
    """
    nr, nc = cube.mask.shape
    rows, cols = np.nonzero(cube.mask)
    rr, cc = np.indices((nr, nc), dtype=float)
    rn, cn = rr.ravel()/max(nr-1, 1), cc.ravel()/max(nc-1, 1)
    background = np.column_stack(((1-rn)*(1-cn), (1-rn)*cn, rn*(1-cn), rn*cn))
    data = cube.data[rows, cols]
    weights = 1/np.sqrt(np.maximum(data, 1))
    # Initialisation uses training pixels only, including when estimating drift.
    image = np.where(cube.mask, np.maximum(cube.data.sum(axis=2), 0), 0)
    initial, lower, upper, pieces = [], [], [], []
    for hint in hints:
        active = np.arange(max(0, round(hint.col)-3), min(nc, round(hint.col)+4))
        local = image[:, active].copy()
        local[np.abs(np.arange(nr)-hint.row) > 12] = 0
        total = max(local.sum(), 1)
        c0 = (local.sum(axis=0)*active).sum()/total
        r0 = (local*np.arange(nr)[:, None]).sum()/total
        variance = (local*(active-c0)**2).sum()/total
        slope = (local*(np.arange(nr)[:, None]-r0)*(active-c0)).sum()/total/max(variance, .01)
        centre = r0+slope*(hint.col-c0)
        width = np.sqrt((local*(np.arange(nr)[:, None]-centre-slope*(active-hint.col))**2).sum()/total)
        envelope = local.sum(axis=0)
        anchor = int(np.argmax(envelope))
        envelope /= max(envelope[anchor], 1)
        start = len(initial)
        initial.extend([centre, max(.5, width), slope])
        lower.extend([max(-1, hint.row-8), .3, -10])
        upper.extend([min(nr, hint.row+8), 12, 10])
        if tail:
            initial.append(max(.1, width*.3)); lower.append(.01); upper.append(15)
        initial.extend(np.delete(envelope, anchor))
        lower.extend([0]*(len(active)-1)); upper.extend([5]*(len(active)-1))
        pieces.append((start, active, anchor, hint.col))

    def profiles(theta):
        out = np.zeros((nr, nc, len(hints)))
        for j, (start, active, anchor, c0) in enumerate(pieces):
            centre, width, slope = theta[start:start+3]
            pos = start+3
            if tail:
                tau = theta[pos]; pos += 1
            envelope = np.insert(theta[pos:pos+len(active)-1], anchor, 1)
            for m, col in enumerate(active):
                mu = centre+slope*(col-c0)
                profile = (_emg_1d(np.arange(nr), mu, width, tau) if tail
                           else np.exp(-.5*((np.arange(nr)-mu)/width)**2))
                out[:, col, j] = profile/max(profile.sum(), 1e-20)*envelope[m]
        return out

    def solve(theta):
        shapes = profiles(theta)
        design = np.column_stack((shapes[rows, cols], background.reshape(nr, nc, -1)[rows, cols]))
        coefficients = _solve_linear(design, data, True, weights, n_nonnegative=len(hints))
        return shapes, coefficients, ((design@coefficients-data)*weights).ravel()

    result = least_squares(lambda t: solve(t)[2],
                           np.clip(initial, np.array(lower)+1e-8, np.array(upper)-1e-8),
                           bounds=(lower, upper), max_nfev=max_nfev,
                           ftol=1e-5, xtol=1e-5, gtol=1e-5)
    shapes, coefficients, _ = solve(result.x)
    prediction = (np.column_stack((shapes.reshape(-1, len(hints)), background))@coefficients)
    volumes = shapes.sum(axis=(0, 1))*coefficients[:len(hints)].sum(axis=1)
    return prediction.reshape(cube.data.shape), bool(result.success), volumes.tolist()


def fit_method(cube, hints, method, original):
    if method == 'original':
        fit = original.refit_spectra(original.fit_shapes_fast(
            cube, len(hints), apex_hint=hints, max_nfev=24,
            max_sigma_row=7, max_sigma_col=2.5), cube)
        return original.predict_overlap(fit, cube), fit.converged, [c.volume for c in fit.components]
    if method.startswith('modulated'):
        return modulated_fit(cube, hints, tail=method.endswith('_tail'))
    fit = fit_overlap(cube, len(hints), model=method, apex_hint=hints,
                       include_background=True, background_order=1,
                       tolerance=1e-5, max_nfev=60)
    return predict_overlap(fit, cube), fit.converged, [c.volume for c in fit.components]


def extract(folder, entry, out):
    geometry = hashlib.sha256(json.dumps(entry, sort_keys=True).encode()).hexdigest()[:10]
    cache = out / (Path(entry['file']).stem + f'_{geometry}_windows.npz')
    if not cache.exists():
        raw = read_run(folder/entry['file'])
        cropped = crop_run(raw, *entry['crop'])
        folded = fold(cropped, entry['period'], entry['phase'], t0=entry['t0'])
        clean = remove_background(cropped, folded).run
        arrays = {}
        for n, window in enumerate(entry['windows']):
            region = CandidateRegion((), window['row0'], window['row0']+window['n_rows'],
                                      window['col0'], window['col0']+window['n_cols'])
            ions = np.searchsorted(clean.mz, window['ions'])
            arrays[f'data{n}'] = _fold_region(clean, folded, ions, region, 0).astype(float)
        arrays['background'] = folded.cube(clean, mz=entry['windows'][0]['ions'],
                                           dtype=np.float64)[325:380, 192:203]
        arrays['source_sha256'] = np.array(entry['sha256'])
        np.savez_compressed(cache, **arrays)
    return np.load(cache)


def recovery_cases():
    a = np.array([10, 5, 2, 1, .3, .1, .5, 1, .2, .4, .7, .8])
    b = a[::-1].copy()
    return {
        'single': [('gaussian', [25, 5, 3, .8], a, 2e6)],
        'tilted': [('tilt', [25, 5, 4, .9, -.75], a, 2e6)],
        'unequal_widths': [('gaussian', [23, 4.5, 4.5, 1.1], a, 2e6),
                           ('gaussian', [28, 5, 1.3, .55], b, 6e5)],
        'hidden_shoulder': [('gaussian', [24, 4.5, 3, .8], a, 2e6),
                            ('gaussian', [27, 5, 3, .8], b, 8e5)],
        'tail': [('emg_t2', [22, 5, 2.5, .8, 6], a, 2e6)],
        'tilted_tail': [('emg_t2_tilt', [22, 5, 2.5, .8, 4, -.6], a, 2e6)],
        'resolved_same_spectrum': [('gaussian', [16, 4.5, 2.8, .8], a, 2e6),
                                   ('gaussian', [34, 5, 2.8, .8], a, 1e6)],
        'tilted_pair': [('tilt', [23, 4.5, 3.5, .9, -.65], a, 2e6),
                        ('tilt', [29, 5, 2, .7, -.4], b, 7e5)],
    }


def recovery(background, ions, run, original_search):
    """Known areas and counts, added to measured background with counting noise.

    Seeds deliberately omit a hidden shoulder and include a false flank seed
    for single peaks. Both algorithms must choose the component count. This
    tests the window fitter and model selection; it bypasses candidate finding.
    """
    rng = np.random.default_rng(19+run)
    rr, cc = np.indices(background.shape[:2], dtype=float)
    settings = ComponentDetectionSettings(workers=1)
    results = []
    for case, peaks in recovery_cases().items():
        signal = np.zeros_like(background)
        truth = []
        for model, parameters, spectrum, volume in peaks:
            profile = _shape(model, rr, cc, np.array(parameters))
            profile /= profile.sum()
            spectrum = spectrum/spectrum.sum()
            signal += profile[:, :, None]*spectrum*volume
            truth.append((spectrum, volume, parameters[:2]))
        data = background+signal+rng.normal(0, np.sqrt(np.maximum(signal, 0)))
        cube = BlobCube(data, np.ones(data.shape[:2], bool), ions, 0, 0, .02, 10, 0)
        screen = replace(cube, data=data[:, :, :8], ions=ions[:8])
        hints = [SimpleNamespace(row=peaks[0][1][0]+1, col=peaks[0][1][1])]
        if len(peaks) == 1:
            hints.append(SimpleNamespace(row=peaks[0][1][0]+8, col=peaks[0][1][1]+1))
        for method, search, options in (
                ('original', original_search, replace(settings, fast_max_nfev=24)),
                ('fixed', select_components, settings)):
            start = time.perf_counter()
            fit = search(screen, cube, hints, options).fit
            cosine = np.array([[np.dot(c.spectrum, s)/max(np.linalg.norm(c.spectrum)*np.linalg.norm(s), 1e-12)
                                for c in fit.components] for s, _, _ in truth])
            # Spatial tie-breaking matters for the resolved identical-spectrum pair.
            distance = np.array([[np.hypot((c.row0-position[0])/10, c.col0-position[1])
                                  for c in fit.components] for _, _, position in truth])
            ii, jj = linear_sum_assignment(1-cosine+.01*distance)
            results.append(dict(run=run, case=case, method=method, true_k=len(peaks),
                                found_k=fit.k, converged=fit.converged, model=fit.model,
                                area_errors=[abs(fit.components[j].volume/truth[i][1]-1)
                                             for i, j in zip(ii, jj)],
                                cosines=[float(cosine[i, j]) for i, j in zip(ii, jj)],
                                seconds=time.perf_counter()-start))
    return results


def measured_profile_recovery(arrays, entry, run, original):
    """Known total areas from copies of an observed, non-analytic peak profile.

    The leading ion of the reference window supplies morphology, not chemical
    identity or a pure-compound spectrum. Its three-column neighbourhood and
    ten-row neighbourhood are copied into measured background with invented
    fixed spectra. Counts are supplied here, isolating shape/area recovery.
    """
    observed = arrays['data0'][:, :, 0]
    apex = np.unravel_index(int(np.argmax(observed)), observed.shape)
    sr, sc = np.indices(observed.shape)
    source = np.maximum(observed-np.median(observed[[0, -1]]), 0)
    source[(abs(sr-apex[0]) > 10) | (abs(sc-apex[1]) > 3)] = 0
    background = arrays['background']
    rr, cc = np.indices(background.shape[:2], dtype=float)
    a = np.array([10, 5, 2, 1, .3, .1, .5, 1, .2, .4, .7, .8])
    b = a[::-1].copy()
    rng = np.random.default_rng(83+run)
    results = []
    for count in (1, 2):
        signal = np.zeros_like(background)
        hints = []
        truth = 0.0
        for r, c, volume, spectrum in [(25, 5, 2e6, a), (29, 5.6, 6e5, b)][:count]:
            profile = ndi.map_coordinates(source, [rr-r+apex[0], cc-c+apex[1]],
                                           order=1, mode='constant', cval=0)
            profile /= profile.sum()
            signal += profile[:, :, None]*(spectrum/spectrum.sum())*volume
            hints.append(SimpleNamespace(row=r, col=c))
            truth += volume
        data = background+signal+rng.normal(0, np.sqrt(np.maximum(signal, 0)))
        cube = BlobCube(data, np.ones(data.shape[:2], bool),
                         np.array(entry['windows'][0]['ions']), 0, 0, .02, 10, 0)
        for method in METHODS:
            _, converged, volumes = fit_method(cube, hints, method, original)
            results.append(dict(run=run, components=count, method=method,
                                converged=bool(converged), total_area_error=abs(sum(volumes)/truth-1)))
    return results


def summarise(out):
    """Write paired comparisons and a standalone scientific figure."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    rows = json.loads((out/'comparison.json').read_text())
    recovered = json.loads((out/'recovery.json').read_text())
    baseline = {(r['run'], r['window']): r for r in rows if r['method'] == 'original'}
    methods = []
    distributions = []
    for method in METHODS:
        selected = [r for r in rows if r['method'] == method]
        ratios = [r['cv_rmse']/baseline[r['run'], r['window']]['cv_rmse'] for r in selected]
        distributions.append(ratios)
        methods.append(dict(method=method, windows=len(selected),
                            median_cv_ratio=float(np.median(ratios)),
                            better_windows=int(np.sum(np.array(ratios) < 1)),
                            converged=sum(r['converged'] for r in selected),
                            median_seconds=float(np.median([r['seconds'] for r in selected]))))
    recovery_summary = []
    errors = []
    for method in ('original', 'fixed'):
        selected = [r for r in recovered if r['method'] == method]
        values = [e for r in selected for e in r['area_errors']]
        errors.append(values)
        recovery_summary.append(dict(method=method, cases=len(selected), matched_peaks=len(values),
                                     correct_count=sum(r['found_k'] == r['true_k'] for r in selected),
                                     converged=sum(r['converged'] for r in selected),
                                     median_area_error=float(np.median(values)),
                                     p90_area_error=float(np.percentile(values, 90))))
    summary = dict(methods=methods, recovery=recovery_summary)
    measured_path = out/'measured_profile_recovery.json'
    if measured_path.exists():
        measured = json.loads(measured_path.read_text())
        summary['measured_profiles'] = [
            dict(method=method,
                 cases=sum(r['method'] == method for r in measured),
                 median_total_area_error=float(np.median([
                     r['total_area_error'] for r in measured if r['method'] == method])))
            for method in METHODS]
    (out/'summary.json').write_text(json.dumps(summary, indent=2))
    labels = ['Original', 'Independent Gaussian', 'Cell-integrated Gaussian', 'Tailing Gaussian',
              'Tilted Gaussian', 'Tailing + tilt', 'Free modulation amplitudes', 'Free amplitudes + tail']
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), layout='constrained')
    axes[0].boxplot(distributions, vert=False, showfliers=False, widths=.55)
    axes[0].set_yticks(np.arange(1, len(labels)+1), labels)
    axes[0].axvline(1, color='#777777', linestyle='--', linewidth=1)
    axes[0].set_xlabel('Held-out error / original error (lower is better)')
    axes[0].set_title(f'{len(baseline)} measured calibration windows')
    axes[0].invert_yaxis()
    for values, label, colour in zip(errors, ['Original', 'Fixed adaptive fitter'], ['#777777', '#0072B2']):
        values = 100*np.sort(values)
        axes[1].plot(values, np.arange(1, len(values)+1)/len(values), label=label,
                      color=colour, linewidth=2)
    axes[1].set_xscale('symlog', linthresh=1)
    axes[1].set_xlabel('Absolute area error (%)')
    axes[1].set_ylabel('Fraction of matched planted peaks')
    axes[1].set_title('Analytic shapes on measured backgrounds')
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(out/'fitting_comparison.png', dpi=180)
    fig.savefig(out/'fitting_comparison.svg')
    plt.close(fig)
    print(json.dumps(summary, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--limit', type=int, default=4)
    parser.add_argument('--summarise-only', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.summarise_only:
        summarise(args.out)
        return
    manifest = json.loads(MANIFEST.read_text())
    original = load_original(manifest['baseline_ref'], 'overlap')
    original_search = load_original(manifest['baseline_ref'], 'modelselect')
    for name in ('fit_shapes_fast', 'refit_spectra', 'predict_overlap'):
        setattr(original_search, name, getattr(original, name))
    results, recovered, measured = [], [], []
    for run, entry in enumerate(manifest['runs'][:args.limit]):
        with (args.folder/entry['file']).open('rb') as source:
            if hashlib.file_digest(source, 'sha256').hexdigest() != entry['sha256']:
                raise ValueError(f"Source differs from the recorded acquisition: {entry['file']}")
        arrays = extract(args.folder, entry, args.out)
        if str(arrays['source_sha256']) != entry['sha256']:
            raise ValueError('Cached windows belong to a different source')
        recovered.extend(recovery(arrays['background'], np.array(entry['windows'][0]['ions']),
                                   run, original_search.select_components))
        (args.out/'recovery.json').write_text(json.dumps(recovered, indent=2))
        measured.extend(measured_profile_recovery(arrays, entry, run, original))
        (args.out/'measured_profile_recovery.json').write_text(json.dumps(measured, indent=2))
        for n, window in enumerate(entry['windows']):
            data = arrays[f'data{n}']
            cube = BlobCube(data, np.ones(data.shape[:2], bool), np.array(window['ions']),
                             window['row0'], window['col0'], .02, entry['period'],
                             entry['t0']+entry['phase'])
            hints = [SimpleNamespace(row=r, col=c) for r, c in window['hints']]
            rr, cc = np.indices(cube.mask.shape)
            groups = (rr//3+cc) % 3
            for method in METHODS:
                started = time.perf_counter()
                prediction, converged, volumes = fit_method(cube, hints, method, original)
                scores = []
                for held_out in range(3):
                    train = replace(cube, mask=groups != held_out)
                    pred, _, _ = fit_method(train, hints, method, original)
                    observed = data[groups == held_out]
                    scores.append(float(np.sqrt(np.mean((pred[groups == held_out]-observed)**2
                                                         / np.maximum(observed, 1)))))
                result = dict(run=run, window=n, method=method, k=len(hints),
                              converged=bool(converged), seconds=time.perf_counter()-started,
                              cv_rmse=float(np.mean(scores)), volumes=volumes,
                              l1=float(np.sum(abs(prediction-data))/np.sum(abs(data))),
                              explained=float(1-np.sum((prediction-data)**2/np.maximum(data, 1))
                                              / np.sum(data**2/np.maximum(data, 1))))
                results.append(result)
                print(json.dumps(result), flush=True)
                (args.out/'comparison.json').write_text(json.dumps(results, indent=2))
    summarise(args.out)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
