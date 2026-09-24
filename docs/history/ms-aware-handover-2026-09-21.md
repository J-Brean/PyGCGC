> **Historical document.** Kept for reference; superseded by the README, ROADMAP.md and
> docs/architecture.md. File paths and numbers describe the project when it was written
> (scripts have since moved from examples/ to benchmarks/, and the 10 s acquisition it refers
> to has left real_data/).

# Claude Code handover: MS-aware overlapping blob finder

Date: 21 September 2026

This is an honest handover of the current state of the work. The requested
outcome is an **automatic blob finder that reports overlapping 2D GC×GC peak
shapes**. Each reported peak should have its own fitted mass spectrum and
estimated intensity. It must not require the user to draw or select a region
before it can find peaks.

The project is a dirty worktree. Do **not** reset or discard changes. Besides
the MS-aware work described here, there is a separate, uncommitted chemical
classification/clustering feature (`chemclass.py`, `cluster.py`, their GUI
state/panels, presets, and related tests) that is outside this handover.

## Start here

Read these files before changing the algorithm:

1. `README.md`
2. `docs/history/handover-2026-09-21.md`
3. `docs/history/blob-diagnosis-and-plan.md`
4. `src/pygcgc/core/candidates.py`
5. `src/pygcgc/core/components.py`
6. `src/pygcgc/core/overlap.py`

The real file is in `real_data/`. Its intended fold is a **10.0 s modulation
period and 1.81 s phase**. The user identified the troublesome area as roughly
**t1 = 15 min, t2 = 1.2 s**. A ten-minute cached crop is in
`real_data/blob_diagnosis/segment.npz`; its metadata are in `segment.json`.

Use this Python interpreter on this Windows machine:

```powershell
C:/Users/BreanJ/venvs/main/Scripts/python.exe
```

## What was wrong with the old system

The legacy detector (`core/peaks.py`) makes a summed chromatogram and assigns
every pixel to exactly one watershed-like blob. This is fundamentally incapable
of representing two co-eluting compounds at the same pixel. It also produced
giant merged blobs and false blobs in noisy background.

The replacement work correctly introduced a different data model:

```
local mass-channel evidence -> fitted component 1 + fitted component 2 + ...
                                     each has position, shape, mass spectrum,
                                     estimated intensity and a warning/status
```

The intended visual result is several overlapping coloured contours. It is not
a set of mutually exclusive borders.

## What is implemented and reasonably solid

### Data model and UI

`Component` and `ComponentGroup` in `core/components.py` persist a component's
position, Gaussian/other shape parameters, spectrum, intensity, confidence and
fit diagnostics. A `ComponentSet` is saved in the project.

The GUI can render each component as overlapping 10% solid and 50% dotted
contours, list components, identify fitted spectra, and export them to CSV.
Relevant code:

- `src/pygcgc/core/components.py`
- `src/pygcgc/gui/mainwindow.py` (`_update_blob_overlay`)
- `src/pygcgc/gui/panels/blobs.py`
- `src/pygcgc/gui/overlap_dialog.py`

The manual overlap dialog is a review/refit tool. It must not be treated as the
normal way components are discovered.

### Automatic candidate and fitting pipeline

`detect_components` in `core/components.py` runs automatically across all
loaded data. Its historical pipeline is:

1. `find_candidates` finds local mass-channel maxima and requires agreement
   between ions plus local TIC evidence.
2. `candidate_regions` groups nearby maxima into small padded windows.
3. `candidate_cube` folds only selected local ions/scans into each window.
4. `fit_shapes_fast` fits one or more overlapping Gaussian shapes on eight
   diagnostic ions.
5. `refit_spectra` estimates spectra and amounts on 30 ions.

The windows are an implementation convenience so the solver has bounded data.
They are not legacy hard blobs and should not be presented as the user-facing
answer.

### Performance work

The full supplied 70 eV LSC file has 215,015 scans and 263 mass channels.
Measured with the previous default configuration:

| Stage | Approximate time |
|---|---:|
| Read | 27 s |
| Background correction | 59 s |
| MS-aware fitting, four processes | 68 s |
| Total | 154 s |

The implementation uses a bounded `ProcessPoolExecutor` queue; it must not
build dense cubes for all 1,277 windows at once. Four workers are the default.
Eight workers measured only about two seconds faster in fitting, at more memory
cost. See `benchmarks/benchmark_full_components.py`.

## Critical mismatch with the user's expectation

The code had an architectural contradiction:

- it **does automatically run** the candidate/fitting pipeline over the whole
  chromatogram;
- but it started the GUI in legacy mode and exposed calculation windows as
  “MS-aware candidate regions”.

That made it look like the user had to find/select regions before anything
happened. This is the core reason the user is unhappy.

The latest uncommitted edits attempted to correct that:

- new `BlobState` and `BlobsPanel` instances default to `MS_AWARE`;
- the normal overlay hides candidate centres and padded rectangles, showing only
  fitted component contours;
- the MS-aware table is empty until automatic component fitting has produced
  components, instead of showing a table of internal regions;
- UI wording says “overlapping blobs” and “inspect/refit selected blob”.

Those presentation changes are directionally correct, but should be reviewed
carefully in the GUI. Existing saved projects keep their saved method. The
legacy method remains selectable for compatibility with downstream tools.

## Automatic model selection (21 September 2026, replaces the residual prototype)

The prototype (`_residual_apex`, `_fit_hidden_component`) was benchmarked and
removed: it made false splitting worse (exact count on 1:1 pairs 0.28 to 0.16).

**Design** (`core/modelselect.py`, `core/components.py`). Candidate seeds are
proposals, not the model order. Each calculation window starts from one peak at
its brightest point and grows one peak at a time. Proposals are unused seeds and
peaks suggested by the leading spectral factor of the positive residual; they are
ranked by a cheap fixed-shape linear fit and only the best two get a nonlinear fit.
An extra peak is accepted only if it beats the current model on BIC, carries at
least 2% of the signal, has a spectrum distinct from the peaks it overlaps (cosine
below 0.90) or is resolved from them (two maxima with a valley between), and
improves held-out spatial-block reconstruction. A backward pass removes peaks made
redundant. What was refused for want of spectral evidence is written into the
group warning. Across windows, duplicates (same place, same spectrum compared on a
common m/z axis) are removed, peaks cut by a window edge are dropped when another
window covers them (`edge_policy`), and interacting neighbouring windows are refitted
jointly. `model_selection="seeds"` restores the old behaviour for comparison.

**Benchmark** (`benchmarks/benchmark_blob_finder.py`, `core/synthetic.py`). Peaks
with real spectra are planted at scan level into the cleaned real segment, at
candidate-free sites, and the whole pipeline runs with no positions supplied.
Tuning (seeds 101 to 110, even-indexed spectra) and held-out (seeds 901 to 910,
odd-indexed spectra) are disjoint. Settings were frozen before the held-out run.

Held-out results, new selection against seed-limited baseline (20 to 90 cases per row; the 1:0.1 exact count is worse than the baseline, see below):

| Case | Exact count | Extra components | Area within 10% | Median area error |
|---|---|---|---|---|
| single Gaussian | 0.85 vs 0.25 | 0.15 vs 0.75 | 0.85 vs 0.70 | 0.02 vs 0.03 |
| broad single | 0.90 vs 0.05 | 0.10 vs 0.95 | 0.85 vs 0.00 | 0.01 vs 0.44 |
| tailing single | 0.80 vs 0.10 | 0.20 vs 0.90 | 0.60 vs 0.30 | 0.08 vs 0.28 |
| pair 1:1 | 0.70 vs 0.28 | 0.06 vs 0.52 | 0.71 vs 0.33 | 0.03 vs 0.30 |
| pair 1:0.3 | 0.56 vs 0.40 | 0.10 vs 0.40 | 0.71 vs 0.26 | 0.03 vs 0.57 |
| pair 1:0.1 | 0.17 vs 0.33 | 0.07 vs 0.40 | 0.91 vs 0.31 | 0.03 vs 0.45 |
| tailing pair | 0.20 vs 0.35 | 0.05 vs 0.40 | 0.44 vs 0.15 | 0.12 vs 0.47 |

Spectral cosine of matched components is about 1.00. Peak recall for 1:0.1 is 0.58
(baseline 0.85, which buys it with three to six components per case); for 1:0.3 it is
0.80 (0.88). Both are worse than the baseline and are the main open weakness.

**Measured limits and open work**

- Close (0.7 to 1.2 widths), unequal pairs are usually merged into one component.
  Recall for a 10% minor peak is 0.6. Neither the held-out gate nor the cosine
  threshold is the cause (each disabled or loosened in isolation changed results only
  within noise); estimated spectra of strongly overlapped peaks are biased together.
- Tailing peaks are Gaussian-only in this fitter; tailing pairs are poor.
- A strong peak elsewhere changes the candidate finder's global noise statistics,
  so 60 of 231 background candidates appeared or vanished in one run. Candidate
  normalisation should be made local.
- Runtime on the ten-minute segment is about 19 s against 4.5 s (four workers).
  Full file: peak finding through the GUI took 141 s in the offscreen run; the
  dedicated full-file benchmark has not been rerun.
- Real crowded data: the full file gives 789 components in 436 groups, 417 flagged
  unresolved. No labelled real standard exists, so this is not a validated count.
- Edge policies compared on tuning: `drop` (default) and `none` and `recentre`
  differ little; `recentre` costs about 20% more time for no measured gain.
## Downstream use, speed and sensitivity (21 September 2026, second pass)

- **Speed.** A blob click cost 5.6 s (MS-aware) or 2.9 s (legacy) on the full file, all
  of it pyqtgraph adding and removing about 3,000 plot items, not spectrum sums. Contours
  are now batched by pen, cached per fitted group, and only redrawn when they change:
  0.49 s and 0.3 s.
- **Speed, 24 September.** Panning or zooming with the overlay open still cost 350 to
  430 ms per step. An isolated pyqtgraph benchmark showed the cause is pen width: any
  width above 1 px sends Qt to its general wide-line stroker (1,600 contours: 1.7 s per
  frame at width 1.8, 8 ms at width 1, with or without antialiasing). Bulk overlays
  (component contours, class and cluster outlines, siloxane) are now width 1; only a
  selected group is drawn wide. Legacy outlines (one segment per pixel edge, ~250,000
  segments) are also merged into collinear runs after the wrap shift
  (`image2d._merge_collinear`, cached per outline). On the supplied file, overlay on:
  MS-aware pan 15 ms, zoom 39 ms, hover 13 ms, click 59 ms; legacy pan 10 ms, zoom
  15 ms, click 18 ms. View culling was tried and removed: no measurable gain at width 1.
- **Blob spectrum.** Clicking a blob (in the view or table) shows that blob's spectrum:
  an MS-aware peak's spectrum re-solved on every ion of its window, or a legacy blob's
  apex spectrum (footprint sum if the apex has fewer than 5 ions).
- **Downstream.** Nothing downstream read `ComponentSet`. Classify and Cluster silently
  ran the legacy detector when no legacy blobs existed (about 22 s each). `core/peaktable.py`
  now re-solves every component's spectrum on all ions (window spectra hold only about 30
  ions and can lack the molecular ion) and presents them as a `BlobSet`-shaped
  `ComponentPeaks`; Classify and Cluster use it (about 6 s each, no detector run). Still
  legacy-only: siloxane blocking (needs pixel labels), hydrocarbon mapping, batch, the
  14 eV soft-spectrum path of Classify, and targeted quantification's all-signals table.
- **Signal accounting** (📊 button, `signal_accounting`). On the supplied file the fitted
  peaks explain 37% of the summed signal, 45% lies outside every fitted footprint. The
  faintest pixels hold under 2% of the signal; most unexplained signal sits at 1e4 to 1e6
  per pixel, largely as a broad unresolved hump plus peak-like leftovers (the report
  separates the two).
- **Sensitivity** (`benchmarks/benchmark_sensitivity.py`, weak planted singles). Detection
  rate by apex TIC, default vs `tic_prominence_snr=1`: 2e4 to 4e4: 0.23 vs 0.35;
  4e4 to 1e5: 0.50 vs 0.81; 1e5 to 3e5: 0.95 vs 1.00. `seed_snr` has no effect. The price:
  87 to 152 components in the segment and 0 to 7 in the reviewed background band.
  Default not changed; `tic_prominence_snr` is exposed in the candidate settings.

## Accuracy evidence and current limits

There is no labelled real standard file yet. Do not claim the method is
accurate or send component amounts to production quantification.

The best available known-area test is
`benchmarks/benchmark_injection_recovery.py`. It injects ideal Gaussian profiles
with two real spectral factors into real local background. With the previous
fast shared-width method, both areas were within 10% in **5 of 9** cases. Every
case where the smaller peak was 10% of the larger failed. This is useful as a
regression test, not validation on realistic peak shapes.

On the cached ten-minute segment, the old automatic result was about 185
components in 120 calculation windows. Around 110 windows were flagged
unresolved, largely because the Gaussian reconstruction explained too little
of the local variation. The full run previously had 1,869 fitted components in
1,277 windows and 1,271 unresolved windows.

Attempts already tried:

- wider t2 Gaussian limits: harmed the known-area test (5/9 to 1/9);
- doubling fast optimizer steps (12 to 24): reduced premature
  non-convergence flags but not unresolved groups;
- process/thread parallelism: four processes materially helped; more did not;
- complete linkage in candidate grouping: fixed enormous chained windows;
- local sparse ion selection: removed a major time/memory bottleneck.

## Recommended next work

### 1. Repair the conceptual pipeline before tuning thresholds

Do not make “candidate regions” a user concept. Keep them internal or rename
them `fit_windows` in code if that helps avoid the confusion. The normal UI
flow should be:

```text
Find peaks -> automatic overlapping component contours and component table
             -> select a component only to inspect/refit it
```

Open the GUI on the real file and check this visually. The normal overlay must
show contours, not rectangles or candidate crosses. Confirm that a fresh
project starts in MS-aware mode.

### 2. Replace seed-count-limited model selection

The target algorithm should generate possible components from the data itself,
then let mass spectra distinguish them. Sensible routes to investigate:

- residual-driven peak pursuit: fit one component, inspect the **spectral
  residual tensor**, add the strongest spatially distinct residual factor, and
  stop only when BIC/held-out reconstruction says to stop;
- non-negative matrix/tensor factorisation (MCR-ALS) inside a sliding window,
  followed by spatial peak fitting; this may be better for non-Gaussian shapes;
- local rank tests on the mass-channel matrix to decide whether one spectral
  factor is inadequate before proposing another;
- use the existing `mcr_als` / `mcr_restarts` implementations in
  `core/overlap.py` as a diagnostic baseline.

For every route, compare k=1, k=2, etc. against a single broad peak and use a
penalty/held-out criterion. Do not split every residual fluctuation.

### 3. Make windows independent of initial apex count

The present candidate finder must still provide enough spatial coverage to make
the process fast, but its detected maxima must not dictate the number of output
blobs. Consider tiling/overlapping small windows around meaningful signal, then
deduplicating component centres afterward. Alternatively, grow windows from
TIC/mass-channel signal islands and use model order to decide components.

This is the main technical job. It needs a clean design and benchmarks, rather
than more small threshold changes.

### 4. Build a proper benchmark suite now

The user will have labelled data in a few days. Before then, create reproducible
synthetic and semi-synthetic cases:

- same position / different mass spectrum;
- partial t1/t2 overlap with unequal areas (1:1, 1:0.3, 1:0.1);
- broad single peaks to measure false splitting;
- non-Gaussian/tailing shapes;
- structured real background from the supplied file.

Measure: correct component count, position error, spectral cosine similarity,
area error, false-positive rate, and runtime. Split tuning and held-out cases.

When labelled data arrives, preserve a small non-sensitive crop and its truth
table as an optional test fixture, if James permits it.

### 5. Keep uncertainty visible

Until validation, retain weak/ambiguous components in the GUI with warnings but
do not use their estimated intensities in production quantification. The user
wants an estimate of badly resolved material, so “unresolved” must mean “shown
with an uncertainty warning”, not silently removed.

## Verification state

Before the final UI/default and residual-search edits, the complete test suite
passed:

```text
393 passed, 25 skipped, 2 warnings in 582.40s
```

The skipped tests require the optional dummy-data pack.

After the latest edits, the focused command was started and should be allowed
to finish or rerun:

```powershell
& C:/Users/BreanJ/venvs/main/Scripts/python.exe -m pytest `
  tests/test_overlap.py tests/test_candidates.py `
  tests/test_gui.py::test_ms_aware_detection_is_a_first_class_peak_method `
  tests/test_gui.py::test_new_project_defaults_to_automatic_overlapping_blobs `
  -q -p no:cacheprovider
```

At minimum, rerun the entire suite after deciding what to retain:

```powershell
& C:/Users/BreanJ/venvs/main/Scripts/python.exe -m pytest tests -q -p no:cacheprovider
```

## Useful commands

```powershell
# full real-file timing; writes ignored JSON under real_data/blob_diagnosis
& C:/Users/BreanJ/venvs/main/Scripts/python.exe benchmarks/benchmark_full_components.py --workers 4

# cached ten-minute automatic benchmark
& C:/Users/BreanJ/venvs/main/Scripts/python.exe benchmarks/benchmark_automatic_components.py --all-regions

# known-area injection regression
& C:/Users/BreanJ/venvs/main/Scripts/python.exe benchmarks/benchmark_injection_recovery.py --fast-only

# launch the application
& C:/Users/BreanJ/venvs/main/Scripts/pygcgc-gui.exe
```

## Communication note

James has been clear that he wants **automatically found, overlapping blobs**,
not a region browser with a manual deconvolution tool bolted onto it. Explain
calculation windows only as internal solver plumbing. Be candid about the lack
of labelled validation, but keep moving on the automatic discovery algorithm
while those data are unavailable.
