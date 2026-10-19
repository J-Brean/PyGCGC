# Code audit, 19 October 2026

The user requested a codebase-wide bug check. This pass found and fixed numerical,
file-import and asynchronous-state bugs that were not detected by the previous suite.
The authorised tilted Gaussian default, fitting thresholds, calibration acceptance
rules and minimum number of calibration levels remain unchanged.

## Scope

- Inventoried and parsed all 174 tracked Python files at the start of the audit,
  spanning 51,965 lines in application code, tests, benchmarks and examples. Checked
  for duplicate definitions and ran Pyflakes across all tracked Python files. No concrete
  undefined production names were found; intentional re-exports and unused imports
  were not treated as functional bugs.
- Compiled source, tests, benchmarks and examples, imported all 88 package modules
  with an offscreen QApplication, and checked the CLI help entry point.
- Ran the complete existing test suite before changes: **722 passed, 31 skipped**.
- Manually traced job lifecycle, project loading and removal, data axes and folding,
  ANDI indexing, target/level CSVs, calibration fitting and inversion, alignment,
  background and siloxane jobs, peak detection, classification, clustering,
  hydrocarbon propagation and quantitative measurement.
- Added deterministic regressions for failures and numerical invariants. This was
  whole-tree automated coverage plus targeted manual review, not a line-by-line
  inspection of every module or a new interactive tour of every screen.

## Findings and fixes

| Problem | Consequence | Correction and regression coverage |
|---|---|---|
| Cancellation after a worker finished but before its queued Qt signal was handled | A cancelled calculation could still apply its result or report an error | Check cancellation at delivery, notify cancellation for superseded jobs, ignore terminal-job progress, and settle jobs even if a callback raises |
| Analysis inputs edited while a worker was running | Old peaks, assignments, phase estimates, background layers or measurements could be accepted against new inputs | Capture source identity, fold/layer state and relevant settings; discard obsolete results and settle running flags without reporting success |
| Alignment propagation applied individual results as they arrived | Deleted points or a changed reference could be repopulated by late jobs | Buffer the batch, require unchanged alignment and source inputs, and apply only a complete successful batch |
| Cancelled classification, clustering or alkane jobs | Controls could remain stuck in a running state | Add cancellation settlement, including partial propagation failures |
| Empty alkane series displayed in the GUI | Taking a minimum over an empty hit list raised an exception | Treat an empty series as unavailable in the panel |
| Quantification rebuilt before a removed run left the project | The rebuilt report still included the removed run | Remove project entries before rebuilding quantitative results; invalidate rather than refit an active measurement |
| Failed target import changed phase anchors before discovering invalid qualifiers | Project anchors and quantitative targets could disagree | Validate the complete table before replacing any project input |
| Target table contained no internal-standard rows | Import or project reload failed even though it was a valid target table | Permit target-only imports and reloads; quantification still requires internal standards |
| Fractional ion masses, non-finite numbers, invalid roles or references, duplicate levels | Truncated masses, invalid arithmetic or silent overwriting of input amounts | Validate target and level CSV headers, roles, masses, finite values, windows, factors and internal-standard references |
| Identified compounds appended by positional column order | Minimal or reordered target tables were corrupted; duplicate promotions made a table unreadable | Respect the existing header, preserve extra columns, reject duplicate appends and replace the file only after a complete successful write; show write failures in the dialog |
| Zero-signal ANDI export used one padding point | The application's reader rejected its own empty export or lost the declared mass axis | Recognise the zero-count placeholder and retain the declared mass range |
| ANDI scan offsets or backwards acquisition times | Padded files assigned intensities to the wrong scans; time filtering could retain a backwards step | Honour scan offsets, validate ranges/counts, and compare scan times against the cumulative preceding maximum |
| Invalid Run axes or FoldSpec values | Non-finite/fractional axes or zero rows failed later in unrelated calculations | Reject invalid axes and fold settings with direct ValueErrors |
| Polynomial calibration formed normal equations on unscaled columns | Changing only amount units could substantially change fitted and back-calculated amounts | Solve the weighted, column-scaled design with an SVD pseudoinverse and derive covariance from the same inverse |
| Quadratic inverse used a fixed small-curvature cutoff and subtraction of nearly equal values | Unit-dependent branch selection or loss of precision | Use the exact linear case and a rationalised quadratic root on the rising branch |

### Calibration evidence

For the known curve `response = 2 + 3*x + 0.2*x*x` at amounts
`[0.25, 1, 5, 10, 25]`, the previous quadratic fit recovered amounts accurately
in the original units. Expressing the same amounts with a factor of 1,000 gave
maximum relative back-calculation error **2.13** (213%); a factor of 1,000,000
gave **10.59** (1,059%). The corrected solve kept the maximum error below
`9e-14` across scale factors from 1 to 1,000,000,000 in this example.

Regression tests cover linear and quadratic models, free and zero intercepts,
all three explicit weighting schemes, several unit scales, inverse-root precision
and uncertainty scaling. Model selection and acceptance criteria were not relaxed.
These known-answer checks establish numerical correctness, not chemical accuracy
on an unknown mixture. Re-run Quantify to regenerate numerical results after updating.

## Verification and evidence

The final complete suite passed **850 tests**, with **31 optional-data skips** and
**two existing SciPy precision-loss warnings** in constant-data statistical-report
fixtures. Runtime was 1,396 seconds on this Windows session. Compared with the
722-test baseline, this adds **128 regression cases**. There were no test failures
or Qt callback exceptions. The final run includes every committed source fix.

The skips cover unconfigured dummy data, legacy native audit/acquisition fixtures
and saved manually labelled calibration references. They are not evidence that
those checks passed against the current private data. The separate read-only
native checks below exercise the four available calibration acquisitions.

Compilation, all 88 package imports and CLI help passed. A wheel built successfully;
its extracted package imported all 88 modules using the existing environment's
dependencies, its CLI help ran, and all 91 source package files were present.
This verifies the Python distribution, not a fresh dependency installation or the
separate PyInstaller Windows executable build.

Read-only native smoke checks reloaded all four calibration acquisitions listed in
`benchmarks/data/calibration_fitting_windows.json`, cropped and folded them, and
repeated the bounded m/z 164 anchor-maximum search. Every refined point increased
intensity, from approximately 48,668 to 334,376; 21,560 to 220,409; 101,572 to
316,396; and 86,305 to 290,603. Raw sizes/modification times and hashes of existing
project JSON files were unchanged. These checks exercise reading, Run, Fold and
refinement on real inputs, not compound identification or mass recovery.

All three target CSVs generated for the user in the preceding task still import
through the GUI's Project API, with respectively 3, 4 and 3 targets/anchors. Their
bytes were unchanged. The four-row table retains its explicit unconfirmed-isomer
entry; successful parsing does not establish that identity.

Tests live in `test_audit_jobs.py`, `test_audit_imports.py`,
`test_audit_state_jobs.py` and `test_audit_calibration_units.py`, plus the existing
Identify GUI regression. Cancellation tests deliberately queue a completed worker
result before cancelling it. State tests deliver results after changing inputs,
rather than depending on thread timing. Import tests check previous input/file
preservation on failure, not just that an exception occurs.

Local logs, JUnit XML, source inventory and numerical reproduction data are under
ignored `real_data/code_audit_2026_10_19/`. Private raw data, compound workbooks,
existing projects and generated caches are not included in this report or commit.
The final full run is `release_suite.log` / `release_suite.xml`; `final_suite.*`
records the earlier intermediate run, before the last numerical and file fixes.

## Limits

Passing tests do not establish that all bugs are gone. The suite contains optional
data-dependent skips, and this pass did not add independent chemical truth for
Peak ID or absolute mass recovery. A workbook inventory and tentative IS amount
are insufficient to certify those. Known real-data limits and unvalidated method
paths remain listed in [ROADMAP](../../ROADMAP.md), including legacy-only tools,
blank subtraction validation and spectral-library identification.
