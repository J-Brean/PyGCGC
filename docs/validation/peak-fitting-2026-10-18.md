# Peak fitting corrections and comparison, 18 October 2026

The production default is now **tilted Gaussian**, at the user's request after
the comparison below. Overlapping peaks have independent widths and correlations;
zero correlation still represents an axis-aligned peak. Automatic component-count
selection and interacting-window merging remain active. Manual GUI refits also
start with tilt. Ordinary Gaussian and other models remain available explicitly.

The original correction and automatic-recovery measurements below used the
preceding adaptive Gaussian/tilt selector. Those figures describe that measured
version; they are not silently relabelled as results of the new fixed-tilt default.

This is a change to the fit, spectrum extraction and recovered amounts. The
confidence thresholds have not been relaxed.

## Follow-up: tilted Gaussian as the default

After the user selected tilt, the same 32 analytic known-area cases were rerun
on the four cached measured backgrounds, with the same noise seeds and misleading
candidate seeds. The fixed-tilt default chose **32/32 component counts correctly**,
converged in **32/32**, and matched **48/48 injected peaks**. Median absolute
matched-peak area error was **0.41%**, with a 90th percentile of **7.65%**.
This is numerical recovery of analytic injections, not a measured accuracy for
real compound amounts. The preceding adaptive result remains reported separately below.

The full default detector was also rerun in all four 120-second calibration crops,
including candidate detection, component-count pursuit, merging, deduplication,
full-ion spectrum extraction and internal-standard matching:

| Filename level | Components in crop | Matched group components | Match spectrum score | Group converged | Explained fraction | Passing candidates |
|---|---:|---:|---:|---|---:|---:|
| L4, `1ng` | 43 | 2 | 0.9867 | yes | 57.05% | 2 |
| L8, `25ng` | 43 | 2 | 0.9946 | yes | 62.45% | 2 |
| L2, `025ng` | 49 | 1 | 0.9867 | yes | 33.72% | 2 |
| L6, `5ng` | 41 | 1 | 0.9882 | yes | 38.64% | 1 |

All four saved-reference matches were recovered. Their groups still require review:
convergence and high spectral similarity do not remove the unresolved-mixture
limitations. More detected components than the adaptive run do not establish more
correct compounds. Source hashes matched the recorded acquisitions and the saved
private project hash was unchanged. The opt-in four-level regression passed in
248 seconds, with one fitting worker and other tests running concurrently.

Previous automatic defaults upgrade to tilt on settings reload. A serialised
`defaults_version=2` marker preserves deliberately selected alternatives in new
projects without exposing an implementation setting in the GUI. Saved reviewed
fits remain stored until recomputation. Low-level fitting API defaults stay
Gaussian for existing scripts; the automatic detector and manual GUI use tilt.

The Peaks panel and manual fit dialog now offer **▶ How it works**, a local
16-second synthetic explanation with the same profiles in 1D and 2D. It animates
initial guesses, independent widths and tilt, and explains why extra components
need evidence. Pause, replay and scrubbing are available; hiding or closing the
dialog stops its timer. Screenshots were inspected in light and dark themes.
This explanation neither runs an optimiser nor changes sample results.

Follow-up results and logs are under `real_data/peak_fitting/tilt_default_*`,
separate from the original comparison artifacts.

Verification of the follow-up: the complete suite passed **652 tests, with 30
skipped and no failures**, in one run (351 seconds). The optional real-calibration
test passed separately across all four levels. After the final settings-marker
and animation layout/theme refinements, the affected fitting, model-selection,
animation and end-to-end workflow modules passed **45 tests**. The animation tests
cover matching 1D/2D profiles, playback, pause, seek, replay, end-of-clip behaviour,
timer shutdown, window reuse, current themes and theme changes while open.

## What was wrong

1. The fast fitter shared both widths between neighbouring peaks. A narrow
   shoulder overlapping a broad peak could therefore have the wrong area even
   with the right component count and almost perfect spectra.
2. The fast shape fit, final diagnostic-ion spectra and full-ion downstream
   spectra solved different objectives. Two stages clipped negative coefficients
   from unconstrained least squares. Clipping does not give the constrained
   optimum when peak profiles overlap. All production stages now use the same
   count weights and non-negative spectra, with a signed local background.
3. The old fast Jacobian held the solved spectra fixed. The replacement projects
   out changes the spectra/background can absorb and includes the residual term
   from differentiating their optimum, on the current active set. Both independent
   and correlated Gaussian derivatives agree with central finite differences.
   This uses the separable least-squares approach described by
   [Golub and Pereyra](https://epubs.siam.org/doi/abs/10.1137/0710036), with the
   derivative evaluated on the constrained solve's active set.
4. EMG profiles were normalised by the maximum over the requested pixels. A fit
   made on a masked training cube consequently predicted a different scale when
   rendered on the full cube. The scale now depends only on the parameters.
5. Automatic selection only tried axis-aligned Gaussians. The labelled-standard
   region visibly shifts in second-dimension time across successive modulation
   slices. The automatic search now compares axis-aligned and tilted Gaussians
   at each proposed component count. Same-spectrum separation and valley checks
   also account for Gaussian correlation.

The default evaluation budget rises from 24 to 60. Serial fitting, like the
existing process workers, limits BLAS threads to avoid launching a thread pool
for each small matrix solve. Two template failures were also corrected: a measured-template refit has two
shift parameters per component, and downstream spectrum extraction must restore
the saved measured profile rather than render a Gaussian from its placeholder widths.
Template accounting footprints now use the actual profile rather than undefined
Gaussian widths, so their fitted signal is not reported as lying outside every peak.

An earlier check which only raised the original fitter's budget from 24 to 120
made the L4 reference fit converge, but its explained fraction stayed around
49 to 50%. A larger budget alone did not resolve the modelling problem.

## The real-data comparison

Eight techniques were compared in 32 fixed windows from the four native 70 eV
calibration files in `C:\REAL_DATA_FOR_PYGCGC\Cal files for James`. Each run
contributes the saved labelled-standard window and the seven largest remaining
windows from the earlier default-detector results. This is a deliberately
difficult, non-random sample of a 120-second retention-time interval, not a
whole-chromatogram accuracy estimate.

The manifest in `benchmarks/data/calibration_fitting_windows.json` records the
file hashes, fold origin, phase, window geometry, ions and original seeds. The
original scientific code is loaded from Git revision
`115df2b47240f589abf9904ff8ce74d4d89ed158`, without changing the working tree.
Twelve ions and the original component count are held fixed per window to isolate
the fitting comparison. The alternatives have independent widths, the same
per-ion count-weighted objective, a signed bilinear background and a budget of
60 evaluations. The original fast pipeline retains its original 24 evaluations,
shared widths and shape/spectrum weighting mismatch.

Three spatial folds with groups `(row // 3 + column) % 3` hold pixels out before
fitting. Every method is scored on the same held-out count-weighted RMSE.
The experimental modulation models have free amplitudes in up to seven columns,
a shared mass spectrum and a Gaussian or EMG second-dimension profile with
linear apex drift. Their initialisation uses training pixels only.

| Technique | Median held-out error / original | Better windows | Full fits converged |
|---|---:|---:|---:|
| Original fast fitter | 1.000 | 0/32 | 13/32 |
| Independent Gaussian | 0.975 | 26/32 | 32/32 |
| Cell-integrated Gaussian | 0.974 | 25/32 | 32/32 |
| EMG along t2 | 1.001 | 16/32 | 32/32 |
| Tilted Gaussian | 0.946 | 22/32 | 32/32 |
| EMG with tilt | 0.948 | 22/32 | 32/32 |
| Free modulation amplitudes | 0.959 | 23/32 | 31/32 |
| Free amplitudes with EMG | 0.956 | 22/32 | 32/32 |

Tilt has the lowest paired median, a **5.4% reduction** in held-out error.
The EMG-with-tilt result is close at 5.2%; these data do not establish a
universal winner between them. The simpler tilted Gaussian is the practical
choice. The median wall time for a full fit plus its three validation fits was
0.19 seconds for the original and 0.88 seconds for tilt, using the benchmark
solver. These are concurrent-workload timings, not full-run throughput; the
production Gaussian solver additionally uses the projected analytic Jacobian.

The table reports the median of paired error ratios, rather than a ratio of
unpaired medians. Different techniques win individual windows. More flexible
tailing and free-amplitude models do not consistently beat the simpler tilted
Gaussian. The first corrected production choice compared ordinary and tilted
Gaussians. The user subsequently selected tilt as the default; this fits the
correlation parameter for every component, including values near zero.

## Known-answer recovery in measured backgrounds

Eight scenarios were added to a quiet 55-row by 11-modulation patch from each of
the four cleaned runs: an isolated Gaussian, a tilted peak, unequal-width
overlap, a hidden shoulder, tailing, tailing with tilt, resolved identical-spectrum
peaks, and a tilted pair. Known areas and spectra were injected with counting
noise using fixed seeds. The background itself is measured signal.

Both automatic searches start with incomplete or excessive seeds and must choose
the component count. This tests the window fitter and model selection; it
bypasses candidate finding. A one-to-one spectral match with spatial tie-breaking
scores the recovered areas. Missing or extra components are reported separately
through the count result.

| Check | Original | Fixed adaptive fitter |
|---|---:|---:|
| Correct component count | 26/32 | 32/32 |
| Converged | 21/32 | 32/32 |
| Matched planted peaks | 46/48 | 48/48 |
| Median absolute area error, matched peaks | 7.03% | 0.35% |
| 90th-percentile area error, matched peaks | 54.5% | 7.6% |

The original missed two planted peaks in the one-to-one matching; errors above
are for matched peaks. A tilted overlapping pair, for example, changed from
roughly 50% area errors to less than 0.5% in all four measured backgrounds.

These are instrument-signal areas, not nanograms or a calibration of the August
mixture. Filename concentrations and the tentative 1 ng internal-standard amount
were not used as truth. Strong tailing remains a harder case for the Gaussian
family, which is visible in the upper tail of the error distribution.

## Recovery with copied observed peak shapes

To avoid judging area recovery only on analytic shapes, the strongest ion of
each reference window supplies an observed profile. Its local background is
removed, the profile is limited to ten rows and three columns either side of
the apex, then translated and normalised to a known area. One copy or an
overlapping pair with invented fixed spectra is added to measured background
with counting noise. These eight cases supply the correct component count to
each fitter, and score total area rather than individual compound identity.

| Technique | Median absolute total-area error |
|---|---:|
| Original fast fitter | 30.8% |
| Independent Gaussian | 28.6% |
| Cell-integrated Gaussian | 28.4% |
| EMG along t2 | 28.0% |
| Tilted Gaussian | 18.4% |
| EMG with tilt | 19.6% |
| Free modulation amplitudes | 10.2% |
| Free amplitudes with EMG | 9.8% |

This is a material limitation of the Gaussian family: the 0.35% analytic-case
result does not describe arbitrary observed shapes. Free amplitudes with EMG
give the best total-area recovery in this small test. They remain an experimental
benchmark alternative because they did not give the lowest held-out error in
the 32 real windows and have more fitting freedom. Tilt is the practical automatic
choice for reconstruction, not a universal winner for quantitative recovery.

## Default detector on the real crops

The saved reference was found in all four crops. All four selected fit groups
converged; the earlier audit selected a non-converged group in every run.

| Filename level | Passing reference candidates | Spectrum cosine | Local explained fraction | Components in matched group |
|---|---:|---:|---:|---:|
| L2, 025ng | 1 | 0.9867 | 33.3% | 1 |
| L4, 1ng | 2 | 0.9867 | 57.0% | 2 |
| L6, 5ng | 1 | 0.9882 | 38.6% | 1 |
| L8, 25ng | 1 | 0.9946 | 62.5% | 2 |

Filename levels are labels rather than verified mixture concentrations. The
explained fraction is a count-weighted comparison with zero signal, not a
chemical accuracy score; selected group boundaries also change with detection.
The L2 and L6 reference fits still explain only 33% and 39% of their local data.
All four reference matches therefore retain review warnings. This remains a
real modelling limitation, despite convergence and the controlled recovery gains.

The serial 120-second crop checks took 164 to 407 seconds per file on this
machine while other checks ran. More shape models, independent widths and a
larger iteration budget cost computation. The GUI normally uses up to four fit
processes; two real windows were also checked for agreement between serial and
process-pool fitting. The opt-in four-file test permits 30 minutes, without
weakening its recovery or warning assertions.


The manually saved reference spectrum is a user label from this dataset, not an
independently acquired pure spectrum. A changed best match or component count
does not establish chemical accuracy. The real data still contain local shape
and spectral structure which these models cannot explain completely; genuine
unresolved fits retain their warnings. Numerical convergence and controlled area
recovery are substantially stronger evidence here than a high spectral cosine
alone, but neither establishes the identities or absolute amounts of co-eluting
real compounds.

## Verification and reproduction

New regression checks fail on the original implementation for EMG subset
invariance, unequal-width area recovery and fitting one tilted compound as one
component with the right amount. Central finite differences independently check
the projected Jacobian with active non-negative constraints and signed background.
The downstream adapter is checked against the fitted-spectrum contract, and the
default workflow test reaches standard assignment and export.

The broad GUI run exposed repeated application-wide style, font and palette
resets whenever another dialog opened. Theme application now updates those
global values only when they change; per-window plot colours and hint styles
still refresh. A regression checks opening another window with the same theme
and switching between light and dark. This also reduces delays in opening the
peak inspection dialogs.

The focused fitting/workflow run passed 73 cases. After the measured-template
restore correction, the science checks passed 60 cases and component display
checks passed 10 cases. The four-run real-data script completed and its recovery,
warning propagation, source hashes and unchanged project hash were checked.
The complete suite passed in fresh processes per test file: 641 passed and 30
skipped across 46 files. After the final legacy-budget and template-footprint
regressions were added, both affected modules were rerun (13 passed). Combined
coverage is **643 passed, 30 skipped, zero failures**, over 673 unique cases.
The monolithic run was interrupted for slow repeated GUI setup; no test files
were omitted from the isolated run.

```powershell
$env:QT_QPA_PLATFORM = "offscreen"
python -m pytest tests -q -p no:cacheprovider
python benchmarks/compare_peak_fitting.py "C:\REAL_DATA_FOR_PYGCGC\Cal files for James" --out real_data/peak_fitting/methods
python benchmarks/validate_calibration_peaks.py "C:\REAL_DATA_FOR_PYGCGC\Cal files for James" --out real_data/peak_fitting/calibration_peaks.json
```

Measurements, extracted cubes and PNG/SVG figures are ignored under
`real_data/peak_fitting/`; the script, manifest and this report are tracked.
Raw files and the existing project remain read only. Cached detections invalidate
automatically because the changed scientific modules contribute to the existing
code fingerprint. Old projects which saved the former 24-evaluation default
upgrade to 60; other custom budgets, and budgets explicitly saved with adaptive
fitting, remain unchanged. Restart the application and run Peaks again to apply
the corrected fitter to an existing project. Saved reviewed fits are preserved
until recomputation is requested.
