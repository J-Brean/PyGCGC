# Changelog

Dated record of what was built and what real-data checks found, newest first. Add new
entries at the top; keep the detail that a fresh session needs (decisions, gotchas, what was
not checked). This was section 2 of the old docs/history/handover-2026-09-21.md.

## State of play

**Colleague feedback round, 2026-10-20** (branch `feedback-2026-10-20`; plan in `docs/plans/feedback-2026-10-20.md`).
Four requests, all built.

- *Small screens.* Compact, standard and large text sizes (View → Text size, default from the screen height), windows
  clamped to their screen, a main window that can shrink below 900 px wide (the Review and Calibration pages had been
  forcing 1200 px), the Process menu as one submenu per numbered step (it had two "4" sections and a "6" before the "5"),
  and Tools moved into the menu bar because a drop-down button in a crowded toolbar falls into its overflow arrow.
- *Anchoring window.* Phase check first (ticks cleared by any phase or period change), a header naming the run ("Run 3 of
  8", reference or sample), a live count of anchors against what each transform needs, the run's own image with the
  reference position and snap box drawn on it, an ANCHOR button, and a before/after overlay with offsets, leave-one-out error
  and image match. Clicks are held (`AlignState.holding`) until ANCHOR, so nothing moves under the mouse and blank layers
  are not invalidated by every click. Two existing expectations changed deliberately (a click no longer refits at once;
  picking needs the phases ticked).
- *📚 Library matcher* on the NIST 2011 SDF files: spectrum plus retention index (column class 5 % phenyl), largest peaks
  first, every uncertainty a named flag, decoy q-values, decisions saved with a content signature. Column dimensions and the
  oven programme are not needed. NIST's estimated RI is corrected per column class from the bias measured over 45,000
  compounds. About 20 ms per peak on 244,000 spectra. Evidence and limits in
  [the matcher note](docs/validation/library-matcher-2026-10-20.md); thresholds are provisional.
- *➖ Peak-level blank subtraction* (new Assign sub-step) and **amounts as sample minus factor times blank** in 🧪 (James's
  decision of 20 October); unknowns' volumes corrected peak by peak before the surrogate curve. Warns when the blank and
  the sample have different phases and no alignment, found on the real runs below.
- Test totals: see the end of this entry.

**Documentation reconciliation, 2026-10-19.**
Added a documentation index, historical-record index and complete input-table reference.
README, user guide, demo tour, architecture, benchmark instructions, roadmap and agent notes
now agree on the eight-step route, guided anchor actions, tilted Gaussian default, saved
components/ladders, source-data location, current theme and main-project export/provenance.
The roadmap distinguishes implemented export from remaining rich-report work and the
checked Python wheel from the untested Windows executable. Calibration filename amounts
and tentative IS identities remain explicitly unconfirmed.

Dated reports retain their measurements and test totals, with clear scope/supersession
notes and a dated correction to the original calibration IS interpretation. Repaired
13 broken local references, including links to unshipped local figures. The CSV examples
parse through the native readers, all 13 target columns are documented, Python examples
compile, and local file/heading links were checked. This was documentation-only; the latest
software suite remains the 850-pass, 31-skip code audit below.

**Codebase audit and regression fixes, 2026-10-19.**
The previous full suite passed 722 tests, yet deterministic checks reproduced queued
cancellation races, stale-result publication, stuck analysis controls, inconsistent
failed imports and a zero-signal ANDI round-trip failure. Job delivery now honours
cancellation and supersession. Analysis states guard source/settings snapshots;
alignment propagation buffers a batch before applying it. Failed/cancelled class,
cluster and alkane jobs settle, empty alkane series display safely, and run removal
rebuilds quantitative results after removing the entry without refitting an active job.

Target imports are transactional and accept target-only tables, including project
reload. CSV validation rejects invalid masses, roles, references and numeric values,
and duplicate calibration levels. Library-hit promotion respects arbitrary valid CSV
column order, preserves extra columns, rejects duplicate names and writes atomically.
ANDI loading respects scan offsets and zero-count padding, retains an empty export's
mass axis and rejects invalid counts/ranges. Backwards scan filtering and Run/FoldSpec
validation now prevent invalid axes from reaching downstream calculations.

Weighted polynomial calibration now solves the column-scaled design directly rather
than forming normal equations. A known quadratic formerly had 213% maximum relative
back-calculation error after a 1,000-fold amount-unit conversion; the corrected example
stays below `9e-14` across scale factors through a billion. Stable quadratic inversion
also preserves small curvature and near-intercept precision. Scientific model-selection,
acceptance and minimum-level rules are unchanged. Re-run Quantify to refresh results.
See [the audit report](docs/validation/code-audit-2026-10-19.md) for scope, tests and limits.
The final suite passed **850 tests**, with 31 optional-data skips and the two existing
SciPy precision-loss warnings. Added 128 regression cases. Read-only checks on four
native calibration acquisitions, all three generated target CSVs, 88 module imports
and a built wheel also passed. Raw data and existing private projects were unchanged.

**Visible alignment point controls, 2026-10-19.**
Added visible Add points, Delete selected points and Start alignment over buttons.
Deletion supports multiple table rows in one undoable edit. Restart keeps the chosen
reference, clears all points and transform overrides, cancels pending propagation and
returns to P1. Undo restores the previous alignment and stops the new picking session.
These controls sit outside the advanced options; docs and the alignment screenshot were updated.
All 40 affected anchor, alignment and workflow checks passed, including mouse-driven
multi-row deletion, restart and undo/redo restoration.

**Guided anchor collection and bounded refinement, 2026-10-19.**
Align now records a reference click immediately, visits every other loaded chromatogram
for the same named anchor, then returns to the reference and prompts for the next point.
The shared table retains all runs' clicks; labelled crosses show saved positions.
Skip missing peaks or finish with Esc. The first loaded run is the default reference.
Picking uses original run coordinates, rejects stale renders, and ends safely on manual
navigation or input changes. Saved reference rows use the existing project format.

**Refine anchors to local maxima** adjusts saved reference and target positions on their
recorded ion or TIC and current processing layer. Its bounded search defaults to ±1
modulation and ±0.150 s, handles wrapped peaks and keeps unsupported neighbourhoods
unchanged. Reference updates are shared across pairs; the whole operation is undoable.
Cancellation, failure and changed inputs cannot apply a late result. The primary button
switches to Finish picking; limits, propagation and transform diagnostics are under
More alignment options. Reference selection now follows undo correctly. The panel scrolls
instead of compressing controls, and residual plots follow the current theme.
The full suite passed 720 tests with 31 skips; final GUI/workflow and anchor runs
passed 40 and 15 checks respectively. Bounded refinement increased the selected
pixel intensity on all four native calibration runs in a read-only smoke check.
See [the anchor workflow check](docs/validation/anchor-workflow-2026-10-19.md).

**Correct detector-time fitting animation, 2026-10-19.**
Replaced the misleading broad 1D projection with an unfolded detector trace: 250 ms
component pulses repeated at a 10 s illustrative modulation interval. A moving cycle
highlight links the pulse train to a magnified within-cycle trace and its 2D column.
All three views use the same sampled profiles, including fitted component contributions.
Playback, seeking, themes and timer shutdown remain available. The guide and screenshot
now explain this folding relationship. This changes the explanation only; production
fitting and acquisition settings are unchanged.
The animation, component-display and workflow checks passed **35 tests**. Light/dark
screenshots at the start, fitting stage and final frame were inspected.

**Guided calibration, exceptions review, prerequisites and fitting previews, 2026-10-18.**
All four workflow priorities are implemented. A transactional four-page setup connects
named compounds, quantifier/qualifier ions, internal standards, run types, calibration
sets, levels, amounts and output units. Guided inputs are embedded in project snapshots;
external CSVs and raw data are not overwritten. Common mass-unit changes convert entered
amounts. Unknown IS amounts and ambiguous filename concentrations remain blank, and the
operator confirms actual amounts and constant IS additions before applying.

Standards now has an exceptions queue with side-by-side candidate spectra and local fitted
profiles, missing confirmations and per-run selections. Review decisions persist with
source/settings/result signatures; changed inputs require search and review again.
Reviewed positions feed quantitative remeasurement, and confirmed missing standards cannot
normalise targets. Fit warnings remain visible. Failed/cancelled standard searches settle
their running state without reporting successful completion.

Workflow guidance shows prerequisites, links to the step that resolves them and invalidates
completion banners for stale results. Quantification rejects a run whose inputs changed
during processing. Peak previews run production fitting on an interval or current run,
report detection/fitting/merging stages, show observed/fitted/residual maps, support cancellation,
and discard late or stale results. Applying preview settings is an explicit action that
starts the normal full-sequence fit; preview components never replace accepted results.
Exports include declared units, setup confirmation, standard review decisions and freshness.
Setup honours the current curve minimum (three levels by default), zero-level policy and
exclusions. Unit conversion occurs on a completed edit; selected-peak retention precision
is retained. Visual QA fixed missing native wizard navigation and a Windows theme-switch
crash after review graphs closed. The full suite passed 678 tests with 31 skips; final
affected-module/input runs and a native-data preview bring combined coverage to 687 passes,
30 skips and no failures over 717 unique cases. See
[the guided-workflow validation](docs/validation/guided-workflow-2026-10-18.md).

**Tilted Gaussian default and animated fitting help, 2026-10-18.**
At the user's request, automatic MS-aware separation and manual GUI refits now start with
tilted Gaussians, with independent widths. Automatic component-count pursuit, duplicate
removal and neighbouring-window merging remain active for tilt. The search honours an
explicit model rather than silently fitting Gaussians. Older automatic defaults upgrade on
reload; versioned settings preserve deliberate alternative choices. Accepted manual fits
remain overrides when their calculation regions still match after detection. A contextual
**▶ How it works** button plays a local 16-second
synthetic illustration with matching 1D projections and 2D contours, pause, replay and seek.
The animation explains spectra, tilt and evidence for additional components, distinguishes
fitting from compound identification, and stops its timer when hidden or closed.
All 32 analytic recovery cases on measured backgrounds selected the right counts and
converged; median matched-peak area error was 0.41%. The full suite passed 652 tests
with 30 skipped; the private-data regression passed all four calibration levels.
Final affected-module checks passed 45 tests after the layout and settings refinements.

**Peak fitting corrections and method comparison, 2026-10-18** (same branch).
Automatic fitting now compares independent Gaussian widths with correlated
Gaussian shapes, allowing modulation-to-modulation apex drift when BIC supports
it. Shape fitting and final full-ion spectra use the same constrained, count
weighted solve. The projected Jacobian includes the response of the solved
spectra and background, and agrees with finite differences. EMG rendering now
uses a parameter-dependent scale which stays consistent on held-out pixels.
The evaluation budget is 60, and serial fits also limit BLAS threads. Saved
measured-template shapes now survive downstream spectrum extraction. Eight
techniques are compared on 32 calibration windows, with analytic and copied
observed-profile recovery in four measured backgrounds. Tilt has the lowest
median held-out error; the experimental free-modulation EMG has the lowest
copied-profile total-area error. See [the fitting report](docs/validation/peak-fitting-2026-10-18.md).
Opening another dialog no longer resets the application-wide theme when it is
unchanged, avoiding repeated repolishing of every existing widget.

**Simpler workflow and calibration peak check, 2026-10-18** (branch `workflow-real-data-2026-10-18`).
Eight main steps now follow Load, optional Align, Background, Peaks, Standards, Assign, Quantify,
Export. Modulation checks live under Load, specialist assignment tools under Assign, and the
peak edit and display controls sit in menus. Alignment can be skipped explicitly. Loading a
target table no longer re-estimates phase after data preparation. Export writes the available
peak, standard and quantitative tables, assignments, input tables, project snapshot and provenance
without changing the active project's save paths. Missing and stale results are recorded in a manifest.
Failed or cancelled peak detection and quantification no longer report a successfully finished
sequence or leave their running flags set. Opening library search no longer marks identification done.
An empty MS-aware result no longer silently falls back to previously detected legacy peaks when
searching for internal standards.
The full-suite run also exposed a PySide deadlock during background removal: automatic cyclic
collection in a worker raced GUI-object deletion on the main thread. A shared Qt timer now
collects cycles on the GUI thread; reference-counted objects still release normally.

The default detector recovered a saved acenaphthene-d10 label across four real calibration levels,
but with multiple passing candidates and unresolved, non-converged fits explaining only 34–58%
of local variation. These warnings now reach internal-standard matches, the review grid and CSVs.
No detection defaults changed. New tests exercise the default GUI route, interruption recovery,
exported values and inputs, and opt-in real-data recovery. Evidence and limitations are in
[the validation report](docs/validation/workflow-and-peaks-2026-10-18.md).

**Demonstration data, new screenshots, 2026-10-07** (same branch). `examples/demo_data/generate_demo_data.py`
writes seven synthetic ANDI runs (four standards at 1 to 50 ng, urban air, diesel exhaust, wood smoke;
8 s modulation, 50 Hz, m/z 38 to 400, 48 min) with about 150 named compounds, five deuterated internal
standards, impurities (siloxanes, phthalates, BHT, solvent front, overloaded toluene) and artefacts
(per-run phase, retention drift, bleed, modulator stripe, spikes, storage threshold, a detector that bends
over and saturates), the supporting files, an MSP library and the truth. Spectra and retention indices are
approximations, not library values. The CLI quantifies the samples within a few per cent where peaks are
clean; it flags the overlapping pairs and over-reads alkanes on m/z 57 over the diesel hump, as on real
exhaust. The legacy detector finds about 1500 to 4500 blobs per run, mostly small residuals on minor bleed
ions near the storage threshold, much as on the real runs. `docs/images/make_screenshots.py` drives the app
on that pack and retakes every screenshot; README and user guide now show every step, and the guide has a
📈 Calibrate section. Tutorials 01 to 06 run on the pack; 07 is the calibration builder without the window.
Fixed: the 📈 window drew black on screen (a transparent `QMainWindow` style rule); the ✂️ 3D surface raised
IndexError on drawing with pyqtgraph 0.14 (colours must be one row per vertex); labelling an internal standard failed ("Couldn't read that peak") on an MS-aware project that also held legacy blobs, as it does once 🫧 or ⛓️ has run (`IstdState._peaks_for` now follows the selected method, as the click does).

**Toolbar reordered to nine steps, hydrocarbon assignment, 2026-10-07** (branch `blob-display-2026-09-28`).
Toolbar is now 1 Fold, 2 Align, 3 Background, 4 Peaks, 5 IS, 6 Quantify, 7 Classify, 8 Identify, 9 Cluster.
Siloxane and hydrocarbons are sub-menus of Peaks (buttons under the peak table; `WorkflowPanel.add_step(...,
parent="blobs")`: no number, no toolbar button, a Back link, never offered as the next step). New
`core/hcassign.py`: each legacy blob is an n-alkane (on a ladder carbon, near the alkane line in t2, alkane
spectrum, 57/71/85 ratios pass), a branched alkane, cycloalkane, alkylbenzene, alkylnaphthalene or
unassigned. Its ion sets and minimum shares are starting values, **not validated on real data**; ask
before treating them as defaults. `AlkaneSeries` gains `shape` (piecewise, the old behaviour and still the
default; smooth PCHIP; polynomial) so a curved or wonky ladder is followed, with the inverse (t1 to carbon)
taken numerically from the forward curve. The estimated modulation period is now rounded to the nearest
0.1 s (`Project.auto_fold`); a period set by hand or by the method file is untouched. Not done: editing or
dropping single ladder points, saving assignments with the project, outlining assignments on the 2D view.

**📈 Calibration builder: usability pass and applying the curves, 2026-10-06** (same branch). The
window is now four tabs that say what to do next and tick themselves off: 1 Standards (the runs
selected in the main window, reused without reading them again; loading cleans them, finds blobs and
lines up phases with no button; a level amount per run, read from the name, stands in for every
compound without its own), 2 Compounds (Suggest, load or click; divided by the internal standard as
soon as there is one), 3 Curves, 4 Apply to samples (new: the main window's other runs are cleaned and
searched with the same settings, each compound found by position and spectrum, read off its curve
with its standard error and flags, as amounts, concentrations or ratios; results follow every curve
edit; CSV export). Settings and session moved into two menus. Fixed: with noisy data, hundreds of
single-ion noise blobs were suggested as internal standards and could move a run's phase; blobs
smaller than a thousandth of the largest are now left out of both. The solvent-front window was
dropped from matching: a large early peak passed for one on short runs, and on the real standards the
streak rules alone left the same blobs. On the synthetic pack the whole route recovers the planted
sample amounts within 0.2%.

**📈 Calibration builder on real standards, 2026-10-06** (same branch). Run end to end on four
August 2026 standard levels (0.25 to 25 ng, 10 s method) supplied without a compound list
([report](docs/validation/calibration-standards-2026-10-06.md)). New: 🔎 Suggest compounds (blobs
that grow with every level are calibrants; flat ones are IS candidates, named and ticked when their
spectrum carries a shipped internal standard's signature); level amounts read from run names;
phases lined up from blob pairs after the first blob search (one run's bleed-method phase was
0.73 s off, leaving 6 of 800 blobs matchable, 123 after); solvent front and t2 smears left out of
matching; wrap-aware matching; largest raw point per blob with an optional saturation flag; 📂 Open
folder; the main window's period no longer carried over as confirmed. Fixed: folder opening listed
TOF-DS SD01 re-exports and then failed to load them. Findings: four internal standards found in
every level, about half their signal at 25 ng; calibrant signal rose 5 to 15-fold over the 100-fold
range (tube background, plateau near 3e5 counts); 3 of 59 curves pass on all levels, 24 without
25 ng. 591 tests pass.

**📈 Calibration builder, 2026-09-28** (branch `calibrate-2026-09-28`). A separate window (📈 Calibrate
on the main toolbar) for a set of standard runs with its own project: compounds from a CSV (a
`target_table.csv` also reads) or from clicks on the chromatogram, amounts per run, internal
standards ticked in the table; ▶ runs 🧹 and legacy 🔍 blobs; each compound is matched to the blob
holding its position or the nearest apex within 20 s / 0.25 s; ➗ divides by the named or nearest
internal standard; one curve per compound with points left out, put back, retyped or rematched by
a click (a click off every blob adds one). Hand matches are stored as positions, so they survive
recomputing blobs. Sessions save as JSON; results export as four CSVs. Science in
`core/calbuilder.py`, curves by the existing `fit_calibration`. Internal standards with a
quantifier ion anchor the window's phase. Checked on synthetic standards (24 new tests) and, for
matching only, on two real 70 eV runs with the five labelled internal standards (all found in both;
5,000 to 6,000 blobs per run, outlines drawn in 0.08 s by a vectorised pass instead of 10 s blob by
blob). Not yet run on a real calibration sequence, and the curves do not yet feed 🧪.

**After the audit, 2026-09-24** (agreed with James). 🫧 blobs are scored against both measured
references (bleed and cyclic) and keep the larger score; the threshold is 0.15 and blocking removes
the siloxane ions only (on two runs: every siloxane-like blob caught, no hydrocarbon signal lost;
old projects keep their saved settings). The n-alkane qualifier ratios were re-measured (0.60, 0.38)
and the defaults kept, because a gate centred on them admitted branched isomers. Run type, level,
calibration set and volume are edited in the Runs table, saved with the project and written out as
`sequence.csv`; quantification reads them from the runs (QuantState.extras is gone). Alkane ladders
are saved with the project; background layers and 🔍 detections are cached on local disk, keyed by
file, fold, settings and a fingerprint of the code (`StepCache` gained a size limit; the batch cache
uses the fingerprint too). `read_lsc` decodes with array operations: identical output on all nine
real files, 1.2 to 3.1 s instead of 14 to 18 s. 547 tests pass.

**Audit and reorganisation, 2026-09-24** (branch `audit-2026-09-24`; full report in
[docs/validation/audit-2026-09-24.md](docs/validation/audit-2026-09-24.md)). The previously
uncommitted tree was committed as found, then every feature was run on the four tandem
acquisitions now in `real_data/` (8 s modulation) and on semi-synthetic sequences built from them.
Fixed: targets whose true t2 exceeds the period were always "off position" (perylene-d12);
blank `qual_ratios` are now learned from the calibration runs, as the README promised; alkane
seeds in the wrong carbon order are refused instead of growing a ladder backwards; Go on the
internal-standard step raised `KeyError`; the step was missing from the toolbar; a disabled
"3D, milestone 8" tab was left over. New: the 🏷️ IS step writes its labelled standards into the
target table, so they anchor the phase and normalise 🧪 (before, labelling fed nothing).
Reorganised: benchmark and investigation scripts moved to `benchmarks/` (with new
`real_data_audit/` scripts reproducing every audit number), validation reports to
`docs/validation/`, superseded plans and handovers to `docs/history/`, screenshots to
`docs/images/`; `real_data/` is now entirely git-ignored; HANDOVER.md became CLAUDE.md;
known-limitations.md was merged into a rewritten ROADMAP.md; README rewritten; stale statements
in the user guide and architecture notes corrected. 543 tests pass, 28 skip (dummy pack).

**🧩 Unsupervised clustering, 2026-09-21.** Workflow step 10 (`core/cluster.py`, `core/segment.py`,
`gui/cluster_state.py`, `gui/panels/clusters.py`). Peaks (legacy blobs, or the MS-aware components
through `core/peaktable.py`) are grouped by mass spectrum: cosine distance, average-linkage tree fitted
on the largest 5,000 peaks with the rest placed by nearest centroid, and retention deliberately left out
so runs need no alignment. Options: a 14 eV block, HDBSCAN on spectra plus retention (scikit-learn), and
pixel regions that need no peak detection (truncated SVD and k-means on the brightest scans). Pooled runs
give a cluster or region table (percent of each run, optionally PQN) that feeds `pca`, `volcano` and a
PERMANOVA; there are merge, stability (bootstrap adjusted Rand index) and distance-scan tools, a Name
column, and the cluster-coloured chromatogram as a live 2D fill and an 89 or 183 mm figure with a caption
file. CLI: `pygcgc batch --cluster --regions N`. Real-data check, one 70 eV acquisition (4,068 legacy
blobs): 102 clusters at merge distance 0.3 with 2.2% of the volume unassigned, and an adjusted Rand
index of 0.52 against the rule-based 🧬 classes (clusters of 94% alkane, 100% alkylbenzene, 86% siloxane
and 100% naphthalene appeared without any rules). Pixel regions (12) took 11.6 s on 215,015 scans and
separated the C3 and C4 alkylbenzene bands, the aromatic roof-tile band, the aliphatics and the siloxane
bleed. HDBSCAN with default settings collapsed the run into one cluster of 3,695 peaks, so "leaf"
selection is the default and still leaves about 30% of the volume unassigned. **Not checked:** any real
multi-sample set (the cross-run parts are tested on synthetic runs only), tandem 14 eV data, and the
frozen-app build with scikit-learn. Legacy footprints cover most of the plane on this run, so peak
clusters are only as good as those boundaries. Full suite at the time: 509 passed, 25 skipped, plus three
unrelated `test_gui.py` failures (stencil lift, homolog table, settings dialog) that fail identically
with the cluster code switched off. A session that froze mid-work had left duplicate function
definitions in `cluster.py` (the later one silently wins); they were removed. See
[known limitations](ROADMAP.md) and the [ROADMAP](ROADMAP.md) for what to check first.

**🧬 Compound classes, 2026-09-21.** Built as a readable replacement for GC Image's CLIC, run after blob detection.
`core/chemclass.py` classifies blobs from a TOML rule file (`presets/chemical_classes.toml`, 17 classes, ion-set
shares, fuzzy ramps, need/prefer, a tier and a per-blob reason). Evidence: the apex spectrum, the ⛓️ ladder (t2 above
the alkane line) and, when the 14 eV file of the same acquisition is loaded, soft molecular-ion features
(`soft_M_mass`, `soft_M_mod14`, ...) read on the 70 eV pixel grid. GUI: workflow step 7 (`gui/class_state.py`,
`gui/panels/classes.py`), so Identify is now step 8 and Quantify step 9. Checked on `real_data/`
(`docs/validation/classes-2026-09-21.md`): the alkane ladder recovers, and 14 eV recovers 9 of 10 alkane molecular masses but
changes only about 25 of 4068 calls. Not resolved: alkene against cycloalkane at unit mass; sulfur heteroaromatics
(class disabled, M+2 unreliable); whether the alkylbenzene block at t2 above 4 s is real. MS-aware users get blobs
found privately by the legacy detector. Results are not saved in the project. Three `tests/test_gui.py` tests
(stencil lift, homologous series button, settings dialogs) fail while the MS-aware method is the default in
`blob_state.py`, because they expect legacy blobs; they pass with the legacy default.

**Automatic model selection, downstream use and speed, 2026-09-21:** the number of
fitted peaks per window is now chosen by `core/modelselect.py` (grow from one peak;
accept only on BIC, spectral distinctness or resolution, and held-out gain; merge
same-spectrum unresolved peaks with a warning) and duplicates across windows are removed.
The residual prototype was removed after it made false splitting worse. Held-out planted-peak
benchmark (`benchmarks/benchmark_blob_finder.py`): exact count for 1:1 pairs 0.70 (was 0.28),
false splits of broad single peaks 10% (was 95%), median area error for distinct pairs
0.02 to 0.03 (was 0.3 to 0.6). Open: 1:0.1 pairs, tailing pairs, runtime about 4x on the
ten-minute segment. A fixed bug: spectra from different windows were compared on mismatched
ion lists. Second pass: `core/peaktable.py` gives Classify and Cluster the fitted peaks
(spectra re-solved on all ions), which previously ran the legacy detector silently;
a blob click fell from 5.6 s to 0.49 s (contours batched and cached), and on 24 September
pan/zoom with the overlay open fell from ~400 ms to ~15-40 ms per step (bulk overlay pens
at 1 px, since wider pens use Qt's slow stroker; legacy pixel-edge outlines merged into
runs); clicking a blob now shows its own spectrum; 📊 signal accounting
added (on the supplied file the fitted peaks explain 37% of the signal); the sensitivity
trade-off is in `benchmarks/benchmark_sensitivity.py`. See
[docs/history/ms-aware-handover-2026-09-21.md](docs/history/ms-aware-handover-2026-09-21.md). Amounts are still not
validated against a labelled standard and must not feed quantification.

**Real blob audit, 2026-09-20:** James supplied period 10 s and phase 1.81 s
for the acquisition in `real_data/`, highlighting t1 around 15 min and t2 around
1.2 s. See [diagnosis and phased plan](docs/history/blob-diagnosis-and-plan.md).
`benchmarks/diagnose_real_blobs.py` reproduces six detector comparisons using a
cached ten-minute scan crop and exports a small all-ion focus cube, metrics and
figures under `real_data/blob_diagnosis/`. Default ion-summed detection labels
about 99% of a background-dominated comparison band. GUI/batch derive their
noise floor from processed data, unlike quantification's explicit raw floor.
That caller inconsistency and a final watershed ownership defect are now fixed;
fixing the floor alone did not solve excessive coverage. `core/candidates.py`
adds an MS-aware candidate/region path. On the cached benchmark it
finds 231 candidates, retains 21 in the crowded focus and finds none in the
background comparison band. Legacy remains the initial method for compatibility,
but the MS-aware route is now an equal saved choice in workflow step 4. The first
real component audit finds two ion-apex clusters; corrected BIC and repeatable
MCR support two stable factors, while a third is unstable. A known-area injection
benchmark using those spectra and a real background passes both component areas
within 10% in only 4 of 9 cases: equal peaks pass from 0.06 to 0.20 s separation,
a 30% minor peak passes at 0.20 s, and 10% minor peaks fail throughout. This is
an explicit block on feeding deconvolved components into quantification. The
confirmed full-run n-C16 peak now supplies a measured profile. The new
`fit_template_overlap` moves copies of that profile and passes 9 of 9 signal-only
cases. Allowing signed fitted background while keeping compound spectra
non-negative raises the real-background result to 6 of 9: equal and 30%-minor
peaks pass, but every 10%-minor case fails. Quiet candidate-free windows at t2
5–9 s pass a representative case; t2 3 s fails and no comparable empty window
exists near the crowded t2 1.2 s band. Background suitability and weak-component
rejection remain blockers. The plan
separates noise/candidate repair, spectral deconvolution, identifiability validation
and GUI/quant integration. Chemical identities and signal-dependent recovery
limits remain unvalidated; do not interpret milestone completion as that validation.

**First-class visual integration, 2026-09-20:** workflow step 4 is now Peak
detection and separation, with a saved selector for MS-aware separation or
legacy blob borders. `BlobState` owns candidate results and regions beside its
legacy `BlobSet`; its runner executes the selected method across every loaded
run. The peak table switches schemas, the 2D view draws centre crosses and cyan
dashed analysis rectangles, and Resolve selected region opens that region in
the ✂️ dialog without recomputing it. Candidate rectangles are deliberately not
converted into fake hard-edged `BlobSet` objects. Existing siloxane blocking,
hydrocarbon totals and identification still request legacy blobs explicitly
until they gain component-aware implementations. A Component contributions tab draws
each fitted contribution as a translucent colour with 10%/50% contours; overlaps
are intentional and mean both components contribute there. The fit table retains
every estimated amount, adds relative size and confidence, and marks components
below 30% of the largest as low confidence. If any are weak, the status calls the
set an unresolved group. When the run has an alkane ladder, "measured alkane
template" derives a shape from n-C16 or the nearest usable alkane. This remains
read-only dialog output; no project persistence, identification, export or quant
consumes the components yet.

**Real-data checks of 🫧 and ⛓️, 2026-09-20:** full write-up and numbers in
[docs/validation/hydrocarbons-2026-09-20.md](docs/validation/hydrocarbons-2026-09-20.md), reproduced by
`benchmarks/validate_real_data.py`. What changed in the code:
- 🫧 default reference is now the column bleed measured in the real run (m/z 207 base, 73 only about 7%), not the generic Dn pattern (`BLEED_REFERENCE`). Presets (`SILOXANE_PRESETS`), table editing and click-to-pick live in the new 🫧 Siloxane tab (`gui/panels/siloxane.py`, `MainWindow._finish_siloxane_pick`). The picked spectrum comes from `core/pixelspec.spectrum_at` (baseline-subtracted by the modulation's median spectrum, switchable) and `siloxane.reference_from_spectrum`. Panel edits set `Project.siloxane_settings` immediately (flags follow); applied layers are not refitted, the panel reports how many are stale. `ApplySiloxane`'s undo still stores the same settings before and after; settings are not part of undo.
- ⛓️ alkane search is now two-dimensional (`_search_carbon`: t2 window, qualifier gate, sub-scan apex, duplicate guard; `_predict`: local secant; window capped at 0.45 of the local spacing). The old time-only search failed on the real run. Seeds carry t2 and are trusted (ratios reported, not enforced). `HydrocarbonSettings` gained `qual_ratios` (0.66, 0.44), `qual_tol`, `enforce_qualifiers`, `t2_window_s`, `seed_t2_window_s`, `max_gap`, `pair_min_ratio`, `series_tol_carbon`, `series_t2_tol_s`; `series_tol_s` is gone. The ⛓️ panel edits, resets and picks the ratios; the ladder table shows each hit's ratios. Synthetic test and conftest alkanes now carry the measured ratios.
- `find_homologous_series` works in carbon-equivalent space (one member per carbon, t2-consistent). Its diagnostic pairs are unvalidated: no series on this run apart from the (57, 71) positive control.
- `MARKER_OFFSET_CARBONS` is a dict, {pristane: +0.01, phytane: +0.08}, measured; both are within one 10 s column of their alkane.
- Phenanthrene/anthracene by RI could not be measured directly (no resolved anthracene); the RI penalty scale is documented as too small to decide the pair. 🏷️ thresholds were not changed.
- Gotchas: `real_data` has one acquisition, so `propagate_alkane_series` is untested on real data. The alkane gate ends the ladder at C12 because lighter alkanes have lower 71/57 and 85/57. The shell's Bash tool cannot run external programs here; use PowerShell.

**Blob-based 🫧 removal, 2026-09-20.** James wanted siloxane removed by blocking out whole 🔍 blobs chosen by a "% certainty" slider, undoable, instead of per-pixel fitting. `SiloxaneSettings.method` is `"blobs"` (default) or `"pixels"` (the old `remove_siloxanes`, unchanged). `siloxane.remove_siloxane_blobs(run, fold, blobs, settings)` scores every blob (`siloxane_scores`), takes those at or above `score_threshold` and zeroes the scans inside their footprints (all ions, or only the siloxane ions with `block_all_ions=False`). Points to know:
- **Blobs are now detected on `Project.base_run(i)`** (🧹 layer or raw, never 🫧), and 🫧 is applied to `base_run` too. Without this, applying 🫧 would remove the siloxane blobs from the detection and re-applying would stack layers. `Project._layers_for_blobs` replaces the direct `layersChanged -> blobs.mark_stale` link so a 🫧 change no longer marks blobs stale.
- `MainWindow.apply_siloxane` runs 🔍 first when a run has no fresh blobs, then applies; `Project.compute_siloxane` skips runs without fresh blobs (and calls `on_done({})` if none are left). Project load re-runs 🔍 before rebuilding a blob-method 🫧 layer.
- The slider is `SiloxanePanel.threshold` (0 to 100, stored as `score_threshold` / 100). It redraws the dashed-magenta outlines (`Image2D.set_blob_overlay(..., siloxane)`), the 🔍 table's 🫧 column (now a percentage) and the count in the panel. The undo command is the existing `ApplySiloxane`.
- `core/batch._prep_one` honours `method`: for `"blobs"` it detects blobs on the cleaned run, blocks them out, then runs the normal blob step on the result. The cache key includes the peak settings in that case.
- Footprints stop at 1% of the apex, so a faint edge tail survives (0.4% in the synthetic test). A dilation option is the obvious fix if it matters on real data. The 0.3 threshold is still provisional and has not been tried on a real run: check which blobs it flags on the supplied acquisition before trusting it. Column bleed is a smear, not a blob, so it is still the per-pixel method's job.

**Workflow panel, 2026-09-20.** James found the eight squeezed left-hand tabs hard to read and wanted a clear route through the tool: open a step, set its options on the left, press Go, then be prompted to the next step. Now:
- James's follow-up: the step buttons live only on the toolbar (no duplicate grid in the panel), siloxane comes after blob detection, and alignment comes before background removal so a blank can be subtracted through the alignment. Order: fold, align, background, blobs, siloxane, hydro, identify, quant.
- Further round: the left side is one `Workspace` dock holding `MainWindow.left_tabs` (Runs, Tool, Layers; the Tool tab is renamed by `_step_opened`), replacing the stacked Runs/Layers/workflow docks. `add_step` now takes a `why` string shown as "Why here", repeated in the toolbar tooltips and in Help → Workflow guide (`WorkflowPanel.guide_html`). The Process menu is in workflow order with a section per step. Closed docks come back from the Workspace and Spectrum toggle buttons and the Panels ▾ menu on the main toolbar (`reset_layout` restores `_default_layout`, saved at the end of `__init__`). The default theme is now light, in the FT palette (`MainWindow(dark=False)`). `identify_dialog.py` was edited elsewhere during this work and was not touched; `MainWindow` only calls `set_run(name, run, fold, blobs, series)` on it and now applies the theme after creating it.
- `gui/workflow.py`: `WorkflowPanel` (the open step's settings stack, a Go button, a next-step banner), plus `QuantStepPanel` and `IdentifyStepPanel`. Steps are registered in `MainWindow.__init__` with `add_step(key, label, title, blurb, widget, go_label)`. The toolbar buttons are checkable `MainWindow.step_actions`; `_step_progress` prefixes ✓ (done) or ▸ (suggested next) from `WorkflowPanel.progressChanged`.
- The old `bg_dock`, `sil_dock`, `align_dock`, `blobs_dock`, `hydro_dock` and Fold dock are gone. The `show_*_panel` methods call `workflow.show_step(key)`. Layers is a tab beside Runs.
- Toolbar step buttons (`MainWindow.step_actions`) only open a step. `MainWindow.go(key)` maps Go to the existing methods; the `act_background`, `act_siloxane`, `act_find_blobs` and similar actions still run the tools directly and are what the menus and tests use.
- Each step's main method calls `workflow.mark_done(key, message)` when it finishes (blobs, hydro, quant, align through `on_done`), which ticks the step and shows the banner. `remove_background_all` and `remove_siloxane_all` call `clear_done`. `_quant_step_status` keeps the 🧪 step's checklist and Go enabled state current.
- The panels' own primary buttons (`bg_panel.apply_all_btn`, `sil_panel.apply_btn`, `blobs_panel.find_btn`, `hydro_panel.find_btn`, `align_panel.propagate_btn`, `fold_panel.fix_btn`) are hidden, not removed.
- Checked: the GUI tests for background, siloxane, toolbar and workflow pass; the offscreen screenshot confirmed the layout but has no fonts, so real-font readability has not been seen. The full suite had not finished when this was written. Tick state is not saved with the project. Identify is marked done when its window opens.

**GUI completion round, 2026-09-20** (ROADMAP sections 2, 3, 4, 5; details there and in README "Working order, settings and figures"):
- `gui/dialogs.SettingsDialog` is a generic form for any settings dataclass (`_auto_editor`, `_rebuild`; tuples as text, optional numbers as blank, dicts skipped). Used for peaks, hydrocarbons, identify, batch (nested tabs) and reachable from menus/panels.
- New in `core`: `report.plot_chromatogram` / `chromatogram_caption` (the GUI's 📤 Figure…, `MainWindow.save_view_figure`), `peaks.measure_ion_polygon` and `QuantDecisions.footprints` (hand-drawn integration, applied in `quantify_run(footprints=...)`, redone by the fast path only through `requantify`), `QuantSettings.saturated_qualifier` with `TargetHit.quant_ion_used`, `align.blank_mapping` / `resample_between` and `remove_background(to_blank=...)` (blank subtraction after alignment; `Project.compute_background` builds the mapping from `project.align.transforms`), and eight peak models in `overlap.MODELS`.
- `fit_overlap` used to start every width at one pixel and never left it on real 20 ms rows; initial widths now come from the half-maximum profile. Shapes: no real gain over gaussian (ROADMAP section 5).
- GUI: toolbar reordered; ⛓️ panel has 🔗 Find series, ✎ Edit shape and ⚙️ All settings; 🏷️ dialog has Browse library… (`LibraryBrowser`) and ⚙️; Review key I and Process menu start re-integration (`MainWindow.start_reintegrate`, orange `footprint_roi` on the 2D view; Ctrl+Enter accepts); D toggles 2D/1D.
- Not done or not checked: blank-after-alignment has only synthetic tests (one real acquisition); the batch pipeline and CLI still do not subtract blanks; footprint polygons are not shown in the Review tiles; the M+1 fallback needs the M+1 ratio typed into the target table (it is not computed from formula).

All 12 milestones are done (version 0.12.0 -- the number counts completed milestones; it matches the milestone number now only because there happen to be 12 of each). The 2026-09-20 full run (later that day, after the MS-aware component work) passed 388 and skipped 25 requiring the absent dummy data pack (413 tests in all); real-file checks ran against the supplied acquisition.

| Milestone | Content | State |
|---|---|---|
| 1 | Reader, timing, folding, phase, period | Done |
| 2 | App shell: 2D and 1D views, spectrum, runs, fold panel, jobs, undo, projects | Done |
| 3 | 🧹 Background, despiking, noise map, blank subtraction, S/N view, `clean` and `crop` commands | Done |
| 4 | 🧪 Targets and internal standards in every run, calibration sets and bracketing (drift interpolation), weighted curves, all-signal ratios to the nearest IS, semi-quant amounts, review gallery, `quant` command | Done |
| 5 | 📍 Reference run, anchors (shift/affine/smooth), leave-one-out error, batch search, aligned view, Compare (blink/grid/difference/overlay) | Done |
| 6 | 🔍 Ion-summed S/N detection, widths, co-elution flag, chunked/threaded detection, delete/merge/add editing, per-peak blank scaling | Done |
| 7 | 🫧 Siloxane flagging (score) and per-pixel NNLS-style removal, as its own layer on top of 🧹 | Done |
| 8 | ✂️ Unique-ion apex map, geometric splits, joint peak-shape fit (BIC, Jacobian uncertainty), MCR-ALS, cross-run PARAFAC, a 3D surface view | Done |
| 9 | Retention index (reused from ⛓️) and 🏷️ MS identifier: MSP import/export, weighted dot-product search, target-decoy q-values, confidence tiers, promote to target, in-house library | Done, but **never seen a real library** -- see below |
| 10 | ⛓️ n-Alkane ladder (two clicks, then propagate to other runs), group-type stencils recorded in alkane-relative coordinates and lifted onto any run, homologous series, markers | Done |
| 11 | 🗂️ Cached, process-pool batch pipeline (reusing quant), feature table (Hungarian matching plus gap fill), Westgard QC, PCA, volcano plots | Done |
| 12 | 📤 Provenance, Excel/HTML/PDF reports with flagged rows, publication figures, a PyInstaller spec (untested) | Done |

**Native-data investigation, 2026-09-19:** All three instrument-produced LSC files
in `real_data/` have now been decoded by standalone investigative scripts (all
645,045 spectra / 53,674,922 ion records), with TIC/EIC outputs and structural audits.
See [docs/validation/native-lsc-format.md](docs/validation/native-lsc-format.md) for the current
findings and integration proposal.
LSC stores calibrated centroid masses as integer m/z × 1000 and float64 intensities;
direct import is implemented without DAT decoding or CDF conversion in
`core/io/markes_lsc.py`. Registry, GUI file/folder/launch paths, batch discovery,
project source hashes and report provenance are wired up. Folder discovery prefers
untreated 70 eV for unambiguous acquisitions; explicit files/sequence rows can select
14 eV or DBC. `LSCFile.spectrum()` preserves original centroids; `read_lsc()` rounds
and sums into CSR. Native imports run serially in `read_many`, and batch defaults
to one process for these larger data. Scan decoding supports progress/cancellation.

Times use whole acquisition-clock seconds as a shared origin, preserving channel
offsets; original metadata records the policy. The GUI labels unconfirmed modulation
periods and skips automatic native phase estimation until the period is confirmed.
Vendor-spectrum/time-origin validation remains outstanding. Blank alignment is still
not implemented (now explicitly labelled); a separate mass-axis indexing bug in
blank subtraction was fixed. Mixed-energy quantification is rejected. The GUI
disables its 70 eV library search for other native energies. Scientific defaults
remain provisional. Tests live in `test_markes_lsc.py` and `test_lsc_gui.py`; full-file
regressions find `real_data/` or `PYGCGC_REAL_DATA`, with no binary fixtures added.
