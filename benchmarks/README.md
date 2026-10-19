# Benchmarks and real-data checks

For current usage and evidence, start with the [documentation index](../docs/README.md).
The [19 October code audit](../docs/validation/code-audit-2026-10-19.md) records the latest
full-suite result (850 passed, 31 optional-data skips), native checks and wheel verification.

The guided workflow has route-level regressions as well as numerical tests:

```powershell
$env:QT_QPA_PLATFORM = "offscreen"
python -m pytest tests/test_calibration_inputs.py tests/test_guided_features.py tests/test_workflow.py -q -p no:cacheprovider
```

They check known quantities after changing amount units, setup persistence and cancellation,
candidate review, propagation of missing-standard decisions, stale-input rejection, and
preview isolation. These checks exercise the GUI workflow rather than estimating chemical
identification accuracy. See [the workflow validation](../docs/validation/guided-workflow-2026-10-18.md).

Scripts that measure how well the methods work, as opposed to `tests/` (which check that the
code does what it says) and `examples/tutorials/` (which show how to use it). None of them
changes a default. Run them from the repository root with the project's Python environment;
generated outputs belong under the repository's ignored `real_data/`. Private raw data
now lives under `C:\REAL_DATA_FOR_PYGCGC`. The modern calibration scripts take that input
path explicitly; legacy scripts below retain hard-coded original paths and fold settings.

## Calibration peak recovery (October 2026)

`validate_calibration_peaks.py` accepts the data location explicitly and leaves the instrument
files and their saved project untouched. It compares the default MS-aware detector with the
legacy detector in 120-second crops around the saved manually labelled standards. Its JSON
records ambiguity, fit warnings and cross-level recovery. This is not a whole-run accuracy or
amount-recovery benchmark. Findings: [workflow-and-peaks-2026-10-18.md](../docs/validation/workflow-and-peaks-2026-10-18.md).

```powershell
python benchmarks/validate_calibration_peaks.py "C:\REAL_DATA_FOR_PYGCGC\Cal files for James" --out real_data/workflow_audit/calibration_peaks.json
$env:PYGCGC_CALIBRATION_DATA = "C:\REAL_DATA_FOR_PYGCGC\Cal files for James"
python -m pytest tests/test_real_calibration.py -q -p no:cacheprovider
```

The real-data tests are opt-in and skip without this environment variable. They also require
the saved, non-empty manually labelled standard reference expected by the script; the data
folder alone is insufficient. Neither test establishes independently known chemical identity
or spike amounts. Private instrument files are not needed for CI. Older benchmarks below
describe their original datasets and paths; the files have since moved to `C:\REAL_DATA_FOR_PYGCGC` on the development machine.

## Fitting comparison on the August calibration files

`compare_peak_fitting.py` compares the original fitter, independent Gaussians,
cell integration, tailing, tilt, tailing with tilt, and two experimental models
with free modulation amplitudes. It scores three held-out spatial folds in 32
fixed windows, then tests automatic component counts and known areas planted
into four measured backgrounds. The manifest records file hashes, fold origin,
window geometry, ions and the original Git revision; no raw files are modified.

```powershell
python benchmarks/compare_peak_fitting.py "C:\REAL_DATA_FOR_PYGCGC\Cal files for James" --out real_data/peak_fitting/methods
```

It writes JSON measurements and PNG/SVG comparison figures. These are fitting
and recovery tests, not estimates of chemical identification accuracy. See
[the fitting report](../docs/validation/peak-fitting-2026-10-18.md).

Tilted Gaussian is now the production automatic and manual GUI default. The original
32-case analytic recovery figures in the report used the preceding adaptive Gaussian/tilt
selector; a new run uses the current default. Keep outputs in a separate folder when
comparing versions. Low-level fitting functions retain their explicit/default Gaussian
API for existing scripts, so pass `model="tilt"` when calling them directly.

## real_data_audit/: the four tandem acquisitions (September 2026)

Written for the audit of 24 September 2026, on four 70 and 14 eV acquisitions originally
under `real_data/` (8 s modulation, phase 2.77 s). `common.py` still hard-codes that layout
and those fold settings; it does not read `PYGCGC_REAL_DATA`. Its original raw files or
decoded caches must be available in that layout to reproduce the historical results.
Do not point these scripts at the August 10 s calibration acquisitions. Interpretation:
[docs/validation/audit-2026-09-24.md](../docs/validation/audit-2026-09-24.md).
The first script to run decodes the files once into `real_data/.cache/` (about a minute).

| Script | Question | Minutes |
|---|---|---|
| `period_and_phase.py` | Are the period and phase found without being told? | 1 |
| `ladder.py` | Does the ⛓️ n-alkane ladder find the right carbons (checked against 14 eV molecular ions)? | 2 |
| `internal_standards.py` | Are the deuterated standards there, and do 🏷️ and 🧪 find them? Writes a draft target table. | 3 |
| `siloxane.py` | What does 🫧 flag and block, and what else does it remove? | 3 |
| `quant_semisynthetic.py` | Does 🧪 recover known amounts planted in a real matrix? | 6 |

## MS-aware peak finding (the earlier acquisition)

These read `real_data/blob_diagnosis/segment.npz`, a cached ten-minute crop of the earlier
10 s-modulation acquisition (no longer in `real_data/`), so they still reproduce their
published numbers as long as that file is kept. Background:
[docs/history/ms-aware-handover-2026-09-21.md](../docs/history/ms-aware-handover-2026-09-21.md).

| Script | What it measures |
|---|---|
| `benchmark_blob_finder.py` | Planted-peak benchmark of the whole automatic finder: count, position, spectrum, area; tuning and held-out splits |
| `benchmark_sensitivity.py` | Detection rate of weak planted peaks against candidate settings |
| `benchmark_injection_recovery.py` | Known-area two-component injections into real background |
| `benchmark_measured_profiles.py` | The same with a measured n-C16 profile |
| `benchmark_automatic_components.py` | Automatic component fits across the crowded window |
| `benchmark_candidates.py`, `benchmark_components.py` | The first candidate and component audits |

## Written for the earlier acquisition, not yet updated

These take the first `*_70eV.lsc` in `real_data/` and assume a 10 s period and 1.81 s phase,
so they are historical reproduction tools, not general runners for the current private
data folder. Without the original input or crop, they cannot reproduce their reports.

| Script | Report |
|---|---|
| `validate_real_data.py` | [hydrocarbons-2026-09-20.md](../docs/validation/hydrocarbons-2026-09-20.md) (siloxane reference, alkane ratios, markers) |
| `classify_real_data.py` | [classes-2026-09-21.md](../docs/validation/classes-2026-09-21.md) (compound-class rules) |
| `diagnose_real_blobs.py` | [blob-diagnosis-and-plan.md](../docs/history/blob-diagnosis-and-plan.md) (makes `segment.npz`) |
| `benchmark_full_components.py` | full-file timing of MS-aware separation |

## lsc_format/

The investigation that decoded the Markes FT04/SD04 `.lsc` format, and its outputs.
Report: [docs/validation/native-lsc-format.md](../docs/validation/native-lsc-format.md).
