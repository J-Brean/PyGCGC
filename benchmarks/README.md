# Benchmarks and real-data checks

Scripts that measure how well the methods work, as opposed to `tests/` (which check that the
code does what it says) and `examples/tutorials/` (which show how to use it). None of them
changes a default. Run them from the repository root with the project's Python environment;
they read from, and write into, `real_data/` (git-ignored, see `real_data/README.md`).

## real_data_audit/: the current four tandem acquisitions (September 2026)

Written for the audit of 24 September 2026, on the four 70 and 14 eV acquisitions now in
`real_data/` (8 s modulation, phase 2.77 s). Results and interpretation:
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
so on the current files they run but do not reproduce the numbers in their reports.

| Script | Report |
|---|---|
| `validate_real_data.py` | [hydrocarbons-2026-09-20.md](../docs/validation/hydrocarbons-2026-09-20.md) (siloxane reference, alkane ratios, markers) |
| `classify_real_data.py` | [classes-2026-09-21.md](../docs/validation/classes-2026-09-21.md) (compound-class rules) |
| `diagnose_real_blobs.py` | [blob-diagnosis-and-plan.md](../docs/history/blob-diagnosis-and-plan.md) (makes `segment.npz`) |
| `benchmark_full_components.py` | full-file timing of MS-aware separation |

## lsc_format/

The investigation that decoded the Markes FT04/SD04 `.lsc` format, and its outputs.
Report: [docs/validation/native-lsc-format.md](../docs/validation/native-lsc-format.md).
