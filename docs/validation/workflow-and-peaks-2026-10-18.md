> **Dated validation record.** Results, timings and test totals describe the version
> and dataset checked here. See the [documentation index](../README.md) for current
> guidance and superseding evidence. Historical measurements are retained.

# Workflow and default peak detector, 18 October 2026

The detector findings below describe the initial workflow audit, before fitting
was changed. The subsequent [fitting corrections and eight-method comparison](peak-fitting-2026-10-18.md)
supersede those numerical fitting results; the workflow findings remain current.

## What changed

The main route is now Load, Align (optional), Background, Peaks, Standards, Assign,
Quantify, Export. Folding and phase checks sit under Load. Hydrocarbon assignment,
classification, library search and clustering sit under Assign. Specialist analysis,
peak edits and display toggles remain available through menus. Alignment has an explicit
skip action, and reaching the end no longer sends the user back to unfinished earlier steps.

The export action writes the available quantitative and per-run peak tables, internal-standard
matches, assignments, target and calibration inputs, a project snapshot and provenance.
A manifest records missing quantification and stale results. It leaves the active project's
save paths unchanged and references the original raw files. Quantitative internal-standard
values and internal-standard search matches have different filenames.

Fixed misleading completion after failed peak or quantitative jobs, running flags stuck after
cancellation, identification marked complete merely by opening its window, and phase being
re-estimated when loading target assignments after data preparation.
An empty or absent MS-aware result also no longer falls back to hidden legacy peaks during
internal-standard matching.

A full-suite run exposed a Qt deletion deadlock while applying background removal. Native
stack inspection showed the GUI thread waiting in QObject destruction while a worker tried
to destroy Qt wrappers during cyclic garbage collection. The application now disables
automatic cyclic collection and runs collection from a shared GUI-thread timer. A regression
checks that widget cycles are destroyed on the application thread after worker allocations.

## Real calibration check

Input: the four original untreated 70 eV files under
`C:\REAL_DATA_FOR_PYGCGC\Cal files for James`. The saved project supplies the 10 s period,
1.66 s phase and a manually labelled acenaphthene-d10 reference at t1 = 1971.096 s,
t2 = 2.72254 s. Raw files and their existing project were read only.

For each level, take a 120-second crop around the reference, preserving the full run's fold
origin, remove background with default settings, and run the default MS-aware detector.
One worker replaces four solely for resource control; the detection and fit settings are
unchanged. Re-solve spectra over all ions and match to the saved reference using the usual
position and spectral thresholds. Run the legacy detector on the same cleaned crop.

| Filename level | MS-aware spectrum score | Passing candidates | MS-aware components in crop | Legacy spectrum score |
|---|---:|---:|---:|---:|
| L2, `025ng` | 0.9842 | 2 | 46 | 0.9483 |
| L4, `1ng` | 0.9975 | 2 | 45 | 0.9628 |
| L6, `5ng` | 0.9902 | 2 | 35 | 0.9554 |
| L8, `25ng` | 0.9867 | 3 | 45 | 0.9621 |

Both methods recover a match at every level. **Every best MS-aware match belongs to an
unresolved, low-confidence fit that did not fully converge.** The local groups containing the
best matches explain approximately 34%, 49%, 39% and 58% of local variation at L2, L4, L6 and
L8 respectively. The L2, L4 and L6 groups also warn that another proposed peak has a spectrum
too similar to separate reliably.

This supports recovery of a spectral feature across levels. It does not establish reliable
component separation, correct peak counts or component amounts. A high spectral score alone
must not be read as a reliable fit. The saved reference comes from the same analysis system,
and crop-local background estimation differs from processing a complete run. These are
limitations of this check, not independent validation against vendor output.

Fit confidence and unresolved flags now travel into internal-standard matches and CSV exports.
The Standards grid shows Review for a passing but uncertain or ambiguous match. Thresholds
and detector defaults have not changed.

As a separate sensitivity check, increasing the fitting budget from 24 to 120 function
evaluations on the L4 crop made the standard's group converge. Explained variation changed
only from 48.8% to 49.6%, two candidates still passed the spectral match, and the fit remained
unresolved. More iterations alone did not solve separation on this crop. This experimental
setting was not adopted. Its results are retained locally under `workflow_audit/fit_budget_120/`.

## Compound and amount information

The supplied `AEROSOLS\Standards.xlsx` is an inventory of compounds and CAS numbers, with
classes and formula information on later sheets. Its column named RI contains values matching
catalogue codes in the original list, including 1875.10 for n-decane and 1524.12 for
acenaphthene-d10. It is not usable as chromatographic retention-index ground truth.

The calibration filenames suggest 0.25, 1, 5 and 25 ng on tube; `025ng` remains an interpretation
of the filename. The user tentatively recalled 1 ng of internal standards. No unconfirmed amount
has been written into production calibration settings, and no amount-recovery claim is made.
The inventory does not supply reference spectra or labelled retention coordinates for all compounds.

## Tests and reproduction

The earlier suite passed 613 tests with 29 skips. Several GUI fixtures deliberately select the
legacy detector. Those tests are useful for compatibility but cannot validate the default route.
The new tests cover default MS-aware detection of known planted peaks, standard propagation,
assignment and export navigation, exported quantitative values and calibration inputs, optional
alignment, failure/cancellation recovery and carrying uncertain-fit warnings into results.
An end-to-end test follows the eight-step route using the default detector and recovers known
sample amounts before checking the exported files.

The real-data check is explicit and reproducible:

```powershell
python benchmarks/validate_calibration_peaks.py "C:\REAL_DATA_FOR_PYGCGC\Cal files for James" --out real_data/workflow_audit/calibration_peaks.json
$env:PYGCGC_CALIBRATION_DATA = "C:\REAL_DATA_FOR_PYGCGC\Cal files for James"
python -m pytest tests/test_real_calibration.py -q -p no:cacheprovider
```

The real-data regression passed on all four levels. It skips without the environment variable
so private instrument data are not required in CI. Machine-readable results and fitted components
remain under the ignored `real_data/workflow_audit/` directory. No raw data are committed.

Final verification: the complete suite ran 663 cases, with 632 passes, 30 skips and one
calibration-display assertion that ran before its asynchronous image was available. That test
now waits for the correct rendered run and its contours; the calibration GUI and collection
modules then passed all 19 checks. This verifies all 633 applicable cases across the full run
and the corrected module rerun. The optional real-data check passed separately in 109 s.
Five new failure/confidence regression cases also failed when replayed against the original
implementation in memory, and pass with the fixes. The remaining two suite warnings are
precision-loss warnings from synthetic statistical comparisons.

The next scientific benchmark needs independently labelled peaks across the chromatogram,
including overlapping pairs, missed compounds and impurities. It should score one-to-one
recovery, false detections, duplicate detections, retention errors, spectral contamination and
amount recovery separately. The current synthetic suite and the four-level recovery check
cannot replace that benchmark.
