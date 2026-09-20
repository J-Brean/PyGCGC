# Native LSC import: implementation and validation

Native import uses `src/pygcgc/core/io/markes_lsc.py`; forensic evidence remains in
`real_data/NATIVE_LSC_REPORT.md`. Original instrument files have not been changed.

## Implemented

- Strict FT04/DI02/SD04 framing and index validation, count checks, signed-safe
  fractional timestamps, original centroid access and unit-mass CSR construction.
- Duplicate-preserving metadata, source SHA-256, parser revision, channel identity
  and vendor background-processing state.
- Native file opening from Python, CLI and GUI, plus acquisition-aware discovery.
  Folder discovery runs in a GUI worker. Untreated 70 eV is preferred only when
  the alternatives can be grouped unambiguously. Open files and sequence CSVs
  select alternatives explicitly.
- Per-scan import progress, cancellation cleanup, serial native imports and a
  conservative single-process default for native batch processing.
- Project persistence of source selection and fingerprint, including a check for
  changed source contents on reopening. Native identity participates in batch
  cache keys and JSON/Excel/HTML/PDF report provenance.
- Visible energy/processing labels and an unconfirmed-modulation-period indicator.
  Mixed-energy quantification is rejected. The GUI's 70 eV library-identification
  action is disabled for native non-70-eV channels.
- A blank-subtraction mass-axis bug was fixed: blank columns are now selected by
  m/z rather than sample-column position. Alignment transforms are still not
  applied to blanks; the background panel explicitly describes that limitation.

## Evidence

Final regression run, 2026-09-20: **295 passed, 25 skipped, 2 warnings** in
117.44 seconds. Skips require the unavailable dummy data pack. Both warnings are
SciPy precision-loss warnings from existing near-identical-sample statistical
tests. The complete real-file tests and real GUI rendering test ran and passed.

Earlier suite attempts intermittently stalled in the GUI quantification fixture;
the isolated case and the final full suite passed. No cause was established and
no speculative change to the peak-segmentation algorithm was made. If this recurs,
capture thread stacks (`-o faulthandler_timeout=40`) rather than treating native
reader validation as evidence that all GUI threading issues are resolved.

The real-data regression compares all 645,045 scans in the three supplied files
against the earlier independent audit, including file hashes, point counts,
duration, empty scans, total signal, every scan's TIC, and five complete EIC traces.

Untreated 70 eV contains 27,398,105 input points and 27,397,056 unit-mass CSR entries
after 1,049 duplicate bins are summed. Its CSR arrays occupy 220,036,512 bytes.
One observed import took about six seconds; this is an observation on this machine,
not a cross-platform performance guarantee or a peak-memory measurement.

The existing timing checker reports 50.0000 Hz, 0.353 ms RMS jitter, 6.447 ms worst
residual, no gaps and no collisions. A 650–690 s crop exercises signal-conserving
folding, background removal and peak detection while checking that raw arrays and
scan times remain unchanged. Background residuals are deliberately signed: the
existing algorithm keeps small negative values to avoid biasing noise estimates.

Tests also cover malformed/truncated files, timestamp rollovers, mixed channels
inside an unsupported file, empty scans, duplicate unit bins, CDF round trips,
batch selection, project reopening, changed source contents, cancellation and
rendering the complete real 70 eV file in the GUI.

## Not yet independently validated

- Vendor-displayed spectra, absolute intensity normalization and retention-time zero.
- The actual modulation period for this acquisition. The 4 s period used in the
  numerical smoke test is not a determination of the instrument method.
- Analytical recovery, calibration accuracy, identification or separation quality
  on standards. A successful reader regression cannot establish these.
- Other LSC revisions/instruments/acquisition modes or DAT-only decoding.
- Alignment-aware blank subtraction and the scientific defaults listed in ROADMAP.md.

Production times are relative to the whole acquisition-clock second shared by the
energy channels, retaining their offset. Audit CSV times remain relative to each
channel's first scan. This intentional origin difference is explicitly checked.

Use `python -m pytest -q` for the regression suite. Real-file checks look for the
supplied files and audit artifacts in `real_data/`, or under `PYGCGC_REAL_DATA`.
They skip when those local files are unavailable; no instrument binaries were
added to the test fixtures.
