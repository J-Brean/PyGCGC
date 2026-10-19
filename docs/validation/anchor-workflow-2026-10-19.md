# Anchor workflow and local-maximum refinement, 19 October 2026

The former two-click route paired one target run with the reference and stopped.
The new route stores the reference click, visits all other loaded chromatograms
for the same point, then returns to the reference for the next anchor. A shared
table retains every click, while labelled crosses show their image coordinates.
Missing peaks can be skipped without generating a match. Finishing preserves
completed clicks; the in-progress navigation itself is not saved.

Refinement searches each anchor's saved ion, or TIC, on the current raw or
processed layer. It selects the brightest finite positive pixel within ±1
modulation and ±0.150 s by default. Limits are adjustable. Wrapped neighbourhoods
use the adjacent modulation's rows, and returned coordinates remain close to the
original unwrapped position. Equal maxima favour proximity; flat, non-positive or
non-finite regions keep their original point. This is an intensity maximum, not
an intensity-weighted centroid or independent evidence of compound identity.

Reference coordinates are refined once per name and shared by every matching
pair. Legacy projects without explicit reference rows remain supported. Results
are applied atomically with transform recomputation and one undo command. A job
is discarded after anchor, run, fold or processing-layer changes; cancellation
and failure leave saved positions unchanged.

## Checks

The full suite passed **720 tests, with 31 skips and no failures** in 431 s.
Skips require optional dummy/native data; the two existing SciPy precision-loss
warnings come from statistical tests with nearly identical values. Final
dropdown and interface refinements were followed by **40 GUI/workflow checks**
and then **15 anchor checks**, all passing. The latter include actual mouse
clicks on the plot, Skip, Finish and Refine controls, plus compact-table options.
Logs and JUnit reports are in `real_data/anchor_workflow/full_suite.*`,
`final_gui_checks.xml` and `final_anchor_checks.xml`.

- Core regressions check search bounds, ties, repeated refinement, unsupported
  pixels, wrap continuity, acquisition edges and invalid limits.
- GUI regressions follow two complete anchor rounds, automatic first-reference
  selection, custom references, skipping, finishing, stale-image clicks, manual
  navigation, ion consistency, overlays, project-state round trips, row removal
  and reference changes through undo/redo.
- Refinement regressions cover reference and target movement, legacy pairs,
  one-action undo/redo, unchanged raw data, and late/cancelled/failed jobs.
- Light and dark visual checks exercised picking and refinement on the small
  synthetic GUI fixture. Search limits and diagnostics collapse behind More
  alignment options; a scroll area keeps the table and controls usable in the
  narrow workspace. The main button toggles its picking label.

## Native-data smoke check

Four untreated 70 eV calibration acquisitions from `Cal files for James` were
read, using the fitting benchmark's fixed crop, phase and period. In the
approximately 1971 s, 2.72 s region on m/z 164, a local high point was selected
and the simulated click offset by +80 ms in t2. Refinement found a higher
intensity within the requested bounds on all four acquisitions:

| Calibration filename level | Clicked pixel intensity | Refined intensity |
|---|---:|---:|
| L4, 1 ng | 48,668 | 334,376 |
| L8, 25 ng | 21,560 | 220,409 |
| L2, 025ng token | 101,572 | 316,396 |
| L6, 5 ng | 86,305 | 290,603 |

The amounts above describe filename tokens, not independently confirmed tube
amounts. This checks bounded maximum selection on measured signal; it does not
validate chemical correspondence or alignment accuracy against known real drift.
Raw-file sizes and modification timestamps were unchanged, and the saved JSON
files retained their SHA-256 hashes. The check only reads private inputs.
Generated scripts, screenshots and results are under the ignored
`real_data/anchor_workflow/`, including `native_refinement.json`.
