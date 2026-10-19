> **Dated validation record.** Results, timings and test totals describe the version
> and dataset checked here. See the [documentation index](../README.md) for current
> guidance and superseding evidence. Historical measurements are retained.

# Guided workflow implementation and validation, 18 October 2026

The four requested workflow improvements are implemented: calibration setup,
standard-exception review, prerequisite/freshness guidance and isolated peak-fitting
previews. Tilted Gaussian remains the automatic and manual fitting default.

## Calibration setup

Assign and Quantify both open the same four-page wizard. It connects named targets,
quantifier/qualifier ions, internal-standard assignments, run types, calibration
sets, known amounts and units. Exactly one selected MS-aware component or legacy
blob can supply a draft position and suggested quantifier, without inventing a
compound identity. Retention-time precision is retained.

The draft does not alter the active project until Apply. Guided targets and levels
are embedded in project JSON, together with unit and confirmation metadata; original
target/sequence CSVs are not overwritten. Ordinary external CSV loading remains
supported. Common mass units convert the draft numbers; custom units retain them
and require checking before confirmation.

IS amounts are initially unknown unless supplied by existing confirmed inputs.
There is no assumed 1 ng addition. Filename amounts are suggestions; `025ng`,
multiple amounts, unsupported scientific notation and negative tokens are not
silently converted into known amounts. The operator confirms actual amounts and
constant additions before applying. Varying IS additions remain an advanced
amount-ratio calibration-builder workflow.

Setup and readiness use the active calibration policy, including the three-level
default, per-compound overrides, zero-level use and excluded/rejected points.
Scientific requirements are not reduced to make a setup look ready.

![Explicit amounts in the synthetic three-level setup](../images/screenshot_calibration_setup.png)

## Exceptions and quantitative consequences

The Standards panel and grid open an exceptions queue. Unsearched, stale, missing
and uncertain standard/run pairs are listed in sequence order; previous reviewed
choices that become obsolete also return to the queue. Candidate spectra are shown
against the reference, side by side with an alternative. MS-aware candidates also
show observed signal, all local fitted contributions, the candidate contribution
and a residual projection. Residuals here include local background.

Keeping a match records an acknowledgement while retaining fit warnings. Selecting
another candidate records a per-run choice without replacing the global reference;
failed automatic thresholds require an explicit override. Confirming a miss keeps
the standard missing. Reviewed positions guide the independent quantifier-ion
measurement; confirmed missing standards cannot supply its normalisation. Fitted
component intensities do not replace calibrated target areas.

Decisions carry reproducible source, processing, result, reference and matching
signatures. Object identity caches signature computation but is excluded from
persisted signatures. Spectral data are normalised to a consistent numeric type
before signing. Valid decisions survive a complete project reload; changed inputs
require a fresh search and review. Uncacheable processing layers conservatively
require review on reopening. Decision edits and obsolete decisions make dependent
quantitative results stale. Cancelled or failed searches do not report completion.

![Two actual passing candidates in the synthetic screenshot fixture](../images/screenshot_standard_review.png)

## Prerequisites and previews

The main action explains missing inputs and offers a direct prerequisite link.
Standards requires fresh detected peaks and a reference. Alignment requires a
reference and searchable anchors, or an explicit skip. Quantify checks targets,
known amounts, current curve requirements, native-data period confirmation and
energy-channel consistency. Raw peak fitting remains available with a background
warning. Stale completion banners are cleared. Input revisions prevent a running
quantification from reporting success against inputs changed during its jobs.

Peak previews use production detection and full-ion spectrum extraction on a
selected/visible interval or current run, preserving fold origin. Interval previews
retain full second-time profiles and one modulation of padding per side. Aligned
view bounds are mapped back to source-run coordinates. Existing background layers
are reused; crop boundaries can still affect the fit.

Progress distinguishes candidate finding, fitting, merging and full-spectrum
extraction. Observed, fitted, residual and fitted-fraction maps support inspection.
Cancellation, closure and changed inputs discard late results. Accepted components,
project fitting settings and quantitative results are unchanged until the operator
explicitly applies preview settings and starts ordinary all-run fitting.

![Synthetic preview, separate from accepted component state](../images/screenshot_peak_preview.png)

## Evidence

Route tests cover draft cancellation, mass-unit conversion, required confirmation,
project reopening, selected-peak import in both detection modes, per-run review,
missing standards, stale review invalidation, interrupted searches, stale quantitative
jobs, preview isolation and explicit application. Known synthetic sample amounts
are checked after converting the calibration inputs from ng to ug, including export
metadata and Results units. Source CSV bytes remain unchanged.

The complete suite passed **678 tests with 31 skips and zero failures** in one run
(354 seconds). After the final coverage, navigation and interaction refinements,
the affected workflow, fitting, model-selection and theme modules passed **81 tests**.
The final input/GUI checks passed **33 tests**, including importing a selected peak in
both detection modes without rounding its retention time, and typing a completed
mass-unit change without converting intermediate keystrokes. Across the full run,
final affected-module runs and native preview, combined coverage is **687 passed,
30 skipped and zero failures**, over 717 unique cases. The JUnit union is recorded
in `real_data/peak_fitting/guided_test_summary.json`. The workflow implementation is
commit `20a2d18`; the final calibration-curve unit labels in `825ec21` passed the
known-quantity/unit-export regression separately.
The two suite warnings are the existing SciPy precision-loss warnings on nearly
identical synthetic statistical samples.

An opt-in native-data preview passed separately in **55 seconds** on the L4 1 ng
calibration acquisition. It used the original 10 s period and 1.66 s phase, preserving
the full run's fold origin in the 120-second labelled-reference crop. The production
preview returned **43 components and 29 unresolved groups**, with finite observed
and model maps. The saved private project bytes were unchanged. This checks the
real-data preview path, not chemical identity or amount accuracy. The earlier
[fitting comparison](peak-fitting-2026-10-18.md) documents those scientific limits.

Visual QA used a synthetic three-level sequence with two resolved same-spectrum
peaks inside the IS search box, producing a genuine ambiguous standard match.
Wizard pages, side-by-side candidate comparisons and preview maps were inspected
in light and dark themes. It exposed two faults that ordinary route tests missed:
Windows' native wizard style hid the navigation controls offscreen, and reselecting
Fusion on a colour switch led to a native crash during font propagation after a
review window closed. The wizard now uses ModernStyle; Fusion is selected once,
equivalent resolved fonts are retained, and plot redraws after palette changes are
deferred. The full visual route and dedicated theme regression then completed.

## Reproduction

```powershell
$env:QT_QPA_PLATFORM = "offscreen"
python -m pytest tests -q -p no:cacheprovider
python -m pytest tests/test_calibration_inputs.py tests/test_guided_features.py tests/test_workflow.py tests/test_theme_windows.py -q -p no:cacheprovider
$env:PYGCGC_CALIBRATION_DATA = "C:\REAL_DATA_FOR_PYGCGC\Cal files for James"
python -m pytest tests/test_real_calibration.py::test_peak_preview_on_native_calibration_data_preserves_inputs -q -p no:cacheprovider
```

Private-data tests skip without that environment variable. Generated reports,
JUnit records, screenshots and the synthetic visual fixture are ignored under
`real_data/peak_fitting/guided_*`. Only synthetic screenshots are included in the docs.
The real acquisitions and their saved project remain read only.
