# Local data and generated checks

Private instrument files live outside the repository at `C:\REAL_DATA_FOR_PYGCGC` on the
development machine. `Cal files for James` contains the four calibration levels and a saved
project with a manually labelled internal standard. Read those inputs without overwriting them.

This directory holds generated audit results, fitted components, previews and test logs. These
outputs are ignored by Git; only this README is tracked. Scripts from the September audit may
still assume instrument files or cached crops at the former repository-local paths.

The current automatic fitter uses tilted Gaussians. `peak_fitting/tilt_default_*` holds
follow-up test logs, private-data recovery results and animation previews for that default.
The earlier `peak_fitting/methods/` comparison records the eight-method study and its
preceding adaptive-selector recovery results; see the tracked
[fitting report](../docs/validation/peak-fitting-2026-10-18.md) for their scope.

The current calibration check accepts an explicit input folder:

```powershell
python benchmarks/validate_calibration_peaks.py "C:\REAL_DATA_FOR_PYGCGC\Cal files for James" --out real_data/workflow_audit/calibration_peaks.json
```

To run the optional recovery test, set `PYGCGC_CALIBRATION_DATA` to that input folder and run
`python -m pytest tests/test_real_calibration.py -q`. Private files are not needed for CI.
