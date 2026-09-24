# ROADMAP: what's left, what's shaky, what's missing

All 12 planned milestones are built (see README.md's Status table and docs/milestone-specs.md for what each one actually delivered). This file pulls every "provisional," "not built," "not yet checked" and "untested" note scattered across docs/known-limitations.md and CHANGELOG.md's per-milestone write-ups into one place, organised by what to *do* about each one, so a future session (or James) can pick a starting point without re-reading both files cover to cover.

**Keep this current.** When a caveat here gets resolved (real data checks out an assumption, a GUI gap gets filled, a threshold gets properly tuned), delete it or move it to "done" rather than leaving it stale — the same discipline HANDOVER.md itself asks for.

The detailed build specification for replacing hard blob ownership with persistent,
overlapping MS-aware components is
[`docs/MS_AWARE_COMPONENT_IMPLEMENTATION_PLAN.md`](docs/MS_AWARE_COMPONENT_IMPLEMENTATION_PLAN.md).
It defines the data model, automatic pipeline, interface, persistence, downstream migration,
validation gates and phased delivery order.

## 0. The one thing that changes almost everything else

**Native FT04/SD04 LSC import is implemented and tested on the three supplied instrument files.**
Full-file counts, timestamps, TIC and five EIC traces match the independent audit.
The GUI can open native files, and project persistence, channel selection and
cancellation have regression tests. Folder discovery prefers untreated 70 eV;
explicit files or sequence rows select alternatives. Vendor-spectrum equivalence,
the retention-time origin, this acquisition's modulation period and analytical
performance on standards remain unvalidated. A real-data crop exercises folding,
background removal and detection numerically; it does not validate scientific
defaults. See `real_data/NATIVE_LSC_REPORT.md` and HANDOVER.md. The assumptions to
check next are:

| Assumption | Lives in |
|---|---|
| Times in seconds, strictly increasing | `core/io/andi.py` |
| Unit-mass centroids; duplicates within a scan summed after rounding | `core/io/andi.py` |
| One modulation period for the whole run | `core/fold.py` |
| Scan rate steady apart from dropped scans (jitter under 10%) | `core/timing.py` |
| Values below the storage threshold are absent; noise floor = threshold / 3 | `core/background.py` |
| Bleed is smooth along t1 over ~10 modulations | `core/background.py` |
| Nothing is retained in the first 0.5 s of t2 (`PeakSettings.wrap_below_s`) | `core/peaks.py` |
| The brightest value occurring ≥3 times is the detector ceiling | `core/peaks.detect_ceiling` |
| A calibration run has "cal" in its type plus a level | `core/quant.SequenceRun` |

Answered by James on 2026-09-19 (HANDOVER.md section 10): the modulation period varies per run (usually 2 to 13 s); "cycloalkanes 1 and 2" are stand-ins for many C10+ cycloalkanes; air volumes range from a few cc/min for 30 min to several days; no NIST library yet; no .exe needed. Native LSC now provides the unfolded spectra directly. Still needed for validated analysis: this acquisition's confirmed modulation period, vendor spectrum/time comparisons, and the real internal-standard list with concentrations.

## 1. Needs a real ROC curve or a real measurement, not a guess

These aren't bugs — they're placeholders the milestone spec always intended to replace once data existed to tune them against. Each has a test already written and marked `@pytest.mark.dummy`, ready to run the moment the dummy pack or real data is available.

- **🔍 Replacement candidate/deconvolution path.** The 2026-09-20 real-data
  audit and first implementation are in `docs/BLOB_DIAGNOSIS_AND_PLAN.md`.
  Per-ion candidates plus TIC prominence remove candidate centres from the
  reviewed background comparison band while retaining 21 in the crowded focus.
  The fast shared-width fit passes 5 of 9 ideal-profile injections into real
  background (the older unrestricted fit passed 4 of 9);
  equal peaks pass down to 0.06 s, but every 10%-minor case fails. Annotate
  known/isolated/weak/held-out windows before making it the default or using
  fitted factors for quantification. It is now a saved first-class method beside
  legacy blob borders in step 4, with its own table, centre/region overlay and
  direct **Resolve selected region…** action. A movable n-C16 measured-template model
  passes 9 of 9 signal-only cases and 6 of 9 in the chosen real-background
  window after signed background support was added; every 10%-minor case still
  fails. A representative 30%-minor case passes in quiet candidate-free windows
  at t2 5–9 s but not at t2 3 s, so background suitability needs a rejection rule. BIC,
  ion apexes and MCR restarts support two stable mathematical factors near
  15 min / 1.2 s; a third is unstable and identities remain unknown.
  The Overlap dialog exposes border-free candidate regions and overlapping
  coloured contribution maps. Durable component groups are fitted automatically,
  persisted, manually replaceable, undoable, listed as peaks, drawn as overlapping
  contours, searched by the identification dialog and exportable to CSV with warnings.
  A measured-alkane-template option appears once that run has an alkane ladder.
  Complete-link grouping now prevents candidate chains from forming huge calculation
  windows. The ten-minute cache has 120 bounded regions and takes about 13 s; the full
  supplied file takes 68 s for detection/separation with four worker processes after
  preparation (154 s including reading and background correction). Most crowded
  focus groups still fail the present explained-signal gate. Next: measured or locally learned
  peak profiles and held-out prediction gates,
  neighbouring-group reconciliation, residual search and run-wide caching. Do not
  send component volumes to production quantification yet.
  Production spectra use 30 ions after eight-ion count/shape screening. Weak components
  remain visible and explicitly low confidence; do not send their volumes to production
  quantification yet. Four worker processes now cut the supplied full-file separation to
  68 s; a bounded job queue avoids retaining every local cube at once. Increasing the
  fast optimizer limit from 12 to 24 steps reduces premature non-convergence flags on
  the ten-minute segment (37 to 10), but still leaves 110 of 120 regions unresolved.
  Widening the Gaussian width limit worsened known-area recovery (5/9 to 1/9), so it
  is not a route to trustworthy weak-component amounts. The next accuracy decision
  needs annotated known peaks or a known co-eluting standard mix.

- **🫧 `SiloxaneSettings.score_threshold` (0.3).** Still provisional. The reference spectrum is now measured (column bleed from the real run, with a cyclic-peak preset, table editing and click-to-pick in the 🫧 panel), but the threshold has not been chosen from a ROC curve: the real run has no labelled siloxane blobs, and with the bleed reference 873 of 4068 blobs score at or above 0.3 (nearly all in the bleed-dominated last 20 minutes). Label some blobs, then choose it. Run `tests/test_siloxane.py::test_dummy_siloxane_score_separates_real_targets_from_planted_siloxanes` when the dummy pack is to hand.
- **⛓️ Alkane qualifier gate** is now enforced, against measured defaults (71/57 = 0.66, 85/57 = 0.44, 30% tolerance), editable in the ⛓️ panel and pickable from a clicked alkane. On the real run the ratios are steady from C13 to C28 (71/57 0.60 to 0.69, 85/57 0.37 to 0.51) but fall for lighter alkanes (C12 0.54/0.32, C11 0.48/0.25, C10 0.36/0.17), so the gate stops the ladder at C12. Lower the expected ratios or switch the gate off to extend it. The ratios also drift up slowly with carbon number, which a single expected value does not capture.
- **⛓️ Alkane search** was time-only and failed on the real run (see `real_data/HYDROCARBON_VALIDATION.md`); it is now two-dimensional (t2 window around the predicted alkane line, local secant prediction). Validated on one acquisition only. Still to check on a second run: `propagate_alkane_series` (only one real acquisition exists, so it is untested on real data), and ladders that start below C12.
- **⛓️ Pristane/phytane offset** measured: +0.006 C (pristane) and +0.08 C (phytane) from n-C17 and n-C18, replacing the guessed -0.3. At a 10 s modulation both share a column with their alkane (one column is about 0.045 C), so the offset is limited by modulation, and the phytane evidence is weak (the m/z 183 excess is small). A shorter modulation period would resolve them.
- **⛓️ Homologous-series diagnostic ion pairs** (91/105, 82/83, 142/156, 58/59) are **not validated**: the real run gave no series for any of them. The finder was also broken on real data (it locked onto solvent-front streaks and assumed a constant spacing in seconds); it now works in carbon-equivalent space, keeps one member per carbon and checks t2 consistency, and as a positive control the (57, 71) pair recovers the n-alkane ladder from blobs alone (14 members, each within 0.1 C of an integer). 142/156 (alkylnaphthalenes) is not a one-carbon series in any case: each methyl adds about 1.3 C. Needs a sample with these series at usable blob quality (blob detection over-labels this run, see `docs/BLOB_DIAGNOSIS_AND_PLAN.md`).
- **🏷️ Every identification threshold** (`IdentifySettings`'s `ri_window`, `ri_penalty_full`, the forward/reverse/ratio/RI-penalty weights, `confident_score`, `ambiguity`) is tuned against small hand-built synthetic libraries only. **No real library has been loaded yet** ("I'll add the NIST library later" — James). The moment one exists: load it, run identification on real blobs, and expect to retune most of these.
- **🏷️ Phenanthrene/anthracene, BbF/BkF/BaP separation by RI** — tested on the real run for phenanthrene/anthracene. One clean 178 peak (RI 1786) is present, plus a larger 178 peak at RI 1845 that co-elutes with a d10 standard (m/z 188); no resolved anthracene, so the pair could not be measured directly. With one modulation column equal to 5 RI units, and assumed library RIs of 1777 and 1793, the real peak sits on the midpoint and the default settings cannot decide (margin 0.006 at RI 1786, tier flips between 1780 and 1790). The default penalty scale (`ri_penalty_full` 100) gives at most a 0.048 score margin for a 16-unit RI difference against an ambiguity threshold of 0.02; `ri_penalty_full` 30 triples it. Any call needs RI good to about 4 units, which needs a co-injected standard to calibrate the RI offset and sub-column t1 centroiding. BbF/BkF/BaP untested.
- **Detector saturation from peak shape.** `detect_ceiling` uses the brightest value seen at least 3 times. James (2026-09-19) wants it estimated from the data and asks whether peak shape (flat-topped apexes) can do it. Try once real saturated peaks exist.
- **Alignment precision.** Run-to-run drift is sub-second while peaks are about 0.25 s wide, so check 📍 residuals against peak width on real data, and whether re-referencing partway through a campaign is needed.
- **🗂️ Comparisons James wants:** similar-sample pairs, samples against blanks, samples against standard mixes. Check the feature-matching radius on these with real data.
- **Internal-standard pairing.** Default is nearest by normalised (t1, t2) distance; consider making the target table's `istd` column optional (auto nearest) once the real internal-standard list exists.
- **✂️ BbF/BkF and benz[a]anthracene/chrysene** (both in the dummy data) are the natural first real test of `fit_overlap`/`choose_k`/`mcr_als` against the milestone's own 10% acceptance target.
- **🧩 Every default, on a real multi-sample set.** Run **Scan distances** and **Check stability** on real samples before trusting the merge distance (0.3), and try pixel regions pooled across runs with PQN on, then **Test groups**, on samples whose groups you already know differ. The cross-run tables, PQN and PERMANOVA have only seen synthetic runs. Also check that an internal standard or a known spiked compound falls in the cluster or region you expect, and that replicate injections cluster together (a replicate PERMANOVA should not separate them).

## 2. Built in `core/`, not reachable from the GUI

Real, tested functions with no button yet. Each is a candidate for "wire this up" once real use shows it's worth the UI, not a sign anything is missing internally.

Done 2026-09-20: `find_homologous_series` (🔗 Find series in the ⛓️ panel: ion-pair presets, members table, double-click to go to a member) and `save_figure` (📤 Figure… on the toolbar and File menu: the 2D view's window, levels and scale as a 89 or 183 mm figure at 300 dpi, PNG or TIFF, with a sibling caption `.txt` to finish).

- **✂️ `split_by_line` / `split_by_valley`** — hand-drawn geometric splits exist and are tested; the Overlap dialog only wires up the one-click watershed-from-apexes split.
- **✂️ `parafac` / `core_consistency` / `split_half_stability` / `build_run_tensor`** — cross-run PARAFAC needs a "pick a feature, pick some runs" flow the single-run Overlap dialog doesn't have. A natural fit for the batch dialog instead, where several runs are already in view.
- **📤 `build_provenance`'s manual-edit-log arguments** (`blob_edits`, `quant_decisions`, `stencils`) — accepted and tested, but every shipped caller (the batch dialog, both CLI commands) passes none of them, since 🗂️ and 📤 are deliberately outside the per-run Project that accumulates that state. A "report the run open right now" action on the main window would be the natural place to wire this in.

## 3. Missing settings dialogs

Done 2026-09-20. `gui/dialogs.SettingsDialog` builds a form from any settings dataclass (a tab per nested one), and is wired to every one: the active peak method's `CandidateSettings` or `PeakSettings` (Process, 🔍 Peak detection settings…), `HydrocarbonSettings` (⛓️ panel, ⚙️ All settings…), `IdentifySettings` (⚙️ in the 🏷️ dialog), 🗂️ batch (⚙️ Settings…: every step's settings, plus a Feature-table reference box that sets `BatchConfig.reference`), `SiloxaneSettings` (🫧 tab) and 🧪 (existing dialog, now with the saturated-qualifier option). Fields the form cannot show (dictionaries) are left as they are. Still fixed constants with no control: the feature-matching radius and gap-fill window in 🗂️.

## 4. Deliberately simplified, not wrong — but worth knowing the limits of

- **🧹 Blank subtraction now follows the alignment** (2026-09-20; unchecked on real data, which has one acquisition). The blank's cleaned per-ion image is read onto each run's own pixel grid through both runs' alignment to the reference (`align.blank_mapping`, `align.resample_between`); with no alignment it is the old whole-modulation shift. The order to use is therefore: align (📍), then apply 🧹 with the blank. Changing an alignment does not redo an applied blank subtraction; the status bar says to apply 🧹 again. A run with no transform is treated as the reference, so give every non-reference run (blank included) anchors. Two artefacts overlapping on one ion still share a single scale factor.
- **📍 Batch anchor-finding** matches each named anchor independently rather than a full assignment-plus-RANSAC pass across every blob; no small-multiples XIC gallery.
- **🔍 Co-elution clustering** (single-linkage on each ion's own local apex) is a heuristic never checked against a real co-eluting pair.
- **✂️ Component-count validation remains incomplete.** The fit now counts
  spectral coefficients in BIC and uses a consistent weighted objective, local
  background and repeated starts. In the first real crowded region, BIC, ion
  apexes and stable MCR restarts support two mathematical components; spatial
  holdout gives only a small further gain at three. Known-area injection into a
  real background passes the 10% area target for equal peaks down to 0.06 s and
  for a 30% minor peak at 0.20 s. The current fast shared-width fit passes 5 of 9
  cases overall, while 10%-minor peaks still fail. Measured isolated profiles, held-out annotations and
  window/noise perturbations are still required before quantification integration.
- **✂️ MCR-ALS unimodality** is a cheap capping-toward-the-peak approximation, not a true isotonic-regression projection.
- **✂️ PARAFAC2 was never built at all** — PARAFAC assumes every run's own elution profile has the same shape; PARAFAC2 would tolerate per-run retention-time drift, but needs materially more machinery, and there's no real multi-run drift yet to say it's worth building.
- **🏷️ Confidence tiers cover only Schymanski levels 2a/2b/3.** Level 1 needs a co-injected reference standard (the closest thing here is `promote_to_target`, which still needs a person to accept the ID); level 4 needs an exact-mass molecular formula this unit-mass detector can't give.
- **🏷️ MSP parsing** has only been checked against a few hand-written files — quoted per-peak comments and unusual encodings in a real NIST/MoNA/MassBank export aren't handled.
- **🗂️ The feature-matching radius and gap-fill window** (`max_distance`, `t1_half_s`/`t2_half_s`) are fixed constants; whether they should scale with real peak width (from 🔍) is unknown.
- **🗂️ Progress reporting** is a single job message, not a per-run progress bar.
- **📤 Reports never carry a real manual-edit log** (see section 2 above — the mechanism exists, nothing populates it yet).
- **📦 The PyInstaller spec (`packaging/pygcgc.spec`) has never been built, on any platform.** There was no environment with PyInstaller available while writing it. Before trusting it: install PyInstaller, build it, run the frozen app against the dummy data pack, and add whatever `hiddenimports` its own build warnings ask for. It now also lists scikit-learn's usual dynamic imports (for 🧩), untested for the same reason.
- **🧩 What is simplified in the clustering.** Retention time is left out of the distance by default so runs need no alignment to compare (the density method adds it, and then pooled runs should be aligned). Peak clusters use average-linkage on the largest 5,000 peaks and place the rest by nearest centroid. Pixel regions are k-means on truncated-SVD spectra, smoothed afterwards, not a spatially aware model, and their number is chosen by hand. The run table is scaled to percent of each run's own total, optionally PQN; there is no internal-standard normalisation, centred log-ratio (negative values break the existing PCA's log transform) or QC drift correction. Cluster results and names are not saved in the project.

## 5. Real gaps against the original design spec (not simplifications — just not done)

Done 2026-09-20 (all of these except the peak-shape verdict below):

- **✂️ Peak-shape options, checked on real data.** Added `emg_t2` (t2 tail), `emg2` (both), `tilt` (t2 apex drifting with t1), `emg_t2_tilt`, `cell` and `emg_t2_cell` (t1 profile integrated over each modulation) beside `gaussian` and `emg`, all in the ✂️ dialog. On eleven isolated real alkane peaks (C16 to C26, five ions, cleaned data) none helps by more than about 2 to 3% in weighted residual per point, and the weighted residual variance per point is 100 to 1500 for every model (1 would be counting noise alone), so the misfit is not the shape family. What did matter: the fitter started every width at one pixel, so it never left the start on 20 ms rows; initial widths now come from the half-maximum profile. Real alkane peaks are often narrower than one 10 s modulation (fitted t1 width 0.1 to 0.3 columns), so a point-sampled Gaussian in t1 is the wrong model there; `cell` fits that case (recovers 0.3 columns on a synthetic peak), but did not lower the real residual either. Likely limits: co-eluting neighbours inside the window and modulation-to-modulation apex jitter. Gaussian stays the default.
- **1D / 2D.** The 1D tab shows the raw chromatogram (reconstructed trace on top, the unfolded detector stream below). Press D to toggle between it and the 2D image (also View menu).
- **🏷️ Library browser and editor** (Browse library… in the 🏷️ dialog): filter, edit name, RI, formula and CAS in place, delete entries, see the spectrum, save the real entries as MSP.
- **⛓️ Stencil shapes are editable**: select a stencil, ✎ Edit shape, drag the corners, 💾 Save with the same name to replace it (undoable).
- **🧪 Hand-drawn integration footprint**: in the Review gallery press I on a peak (or Process → 🧪 Re-integrate peak…); an orange polygon starts as the automatic footprint's bounding box; drag its corners, Ctrl+Enter to integrate inside it (Esc cancels, ↺ resets to automatic). Stored in the decisions (saved with the project, undoable, counted in provenance), flagged "integrated by hand". Rim baseline is still subtracted.
- **🧪 Saturated quantifier fallback, optional and off by default**: with Quantification settings → "Saturated target: quantify on a qualifier", a saturated quantifier is replaced by the first unsaturated qualifier's volume divided by its expected ratio from the target table (list an M+1 ion there with its ratio to use one). Flagged "saturated (quantified on m/z N)". The expected ratio must be right for unsaturated conditions, and a calibration built on a mix of quantifier-ion and qualifier-ion volumes carries that assumption.
- **🧪 No process pool for measurement** — runs are still measured in the thread pool, not the process pool 🗂️ batch later added for its own per-run prep.

## 6. Open questions only James can answer

The full list (15 questions, with James's answers of 2026-09-19 inline) is HANDOVER.md section 10. Real files now import directly; independent vendor comparisons and real standards are the next validation steps. Headline topics:

- Instrument/method basics: internal standards, cycloalkane identities, air volumes, saturation level source.
- Whether small real files can be committed to the repo for regression tests.
- Real-data tuning for 🔍 (S/N thresholds, co-elution), 🫧 (reference spectrum), ⛓️ (diagnostic ion pairs, marker offset, qualifier gate), 🗂️ (matching radius, volcano groups, worker count), ✂️ (which method he'd trust day to day), 🏷️ (library timeline, threshold tuning).
- Whether a standalone installer (PyInstaller) is worth finishing, or a shared `uv` environment is enough.
- 📍 alignment reference-run strategy for a long campaign.

## Suggested order, if starting fresh with real data

1. Work through section 0 (the reader and its assumptions) before anything else — everything downstream depends on it.
2. Get answers to the blocking instrument-basics questions in section 6.
3. Run the dummy-marked tests once the dummy pack is available, even before real data, to catch anything section 0 didn't (siloxane ROC, hydrocarbon qualifier ratios, 🏷️ library search against `synthetic_library.msp`).
4. Load a real library into 🏷️ and see what breaks (section 1's biggest single item).
5. Everything else in sections 1–5 is independent and can be picked up in whatever order real use surfaces as annoying first.
