# HANDOVER: PyGCGC

Notes for whoever picks this up next, most likely Claude Code working with James and his first real data. Read this file and `README.md` first; every core module also opens with a plain-English docstring explaining its method. The detail that used to live here has moved (see the map below). **`ROADMAP.md`** is the short version: every provisional threshold, missing GUI wire-up and open question in this file and README's Known limitations, pulled into one prioritised list with a suggested starting order -- read it first if the question is just "what should I work on."

**This file exists so a brand-new conversation can pick this project up cold.** James's preferred workflow is to clear the conversation between milestones rather than compress it: a compress keeps a lossy summary of a lot of exploratory back-and-forth (scratch scripts, debugging dead ends) that isn't useful going forward, at real token cost on every later turn, whereas a clear plus a fresh read of this file and `README.md` gets equivalent or better context for less — provided this file stays genuinely complete. So: keep it that way. Update it as you go, not just at the end (state of play, decisions, gotchas, open questions), the same way it's been kept up to date through milestone 9.

**Current component architecture, 2026-09-20:** the MS-aware result is now a sequence of
overlapping peak shapes, each with its own mass spectrum and estimated intensity. Candidate
rectangles are local calculation windows and disappear once components exist. Durable
component/group objects, project persistence, undoable manual acceptance, automatic fitting,
overlapping contours, the component table, component-spectrum identification and CSV export
are implemented. The full specification and remaining stages are in
[`docs/MS_AWARE_COMPONENT_IMPLEMENTATION_PLAN.md`](docs/MS_AWARE_COMPONENT_IMPLEMENTATION_PLAN.md).

The automatic real-data audit uses the cached 14–16 minute segment at a 10 s modulation
period and 1.81 s phase. Six crowded groups produce 14 proposed components. All six are
correctly marked unresolved: their Gaussian fits explain 33–78% of the locally weighted
signal and four groups contain weak components. Screening model order on 20 ions and then
refitting the chosen model on 60 ions reduced this benchmark from about 266 s to 152 s
without changing the component counts. This is useful review output, but it is not yet a
validated measurement of 14 compounds. See
`real_data/blob_diagnosis/automatic_component_metrics.json`.

## Where things live

| Question | File |
|---|---|
| What has been built and found, by date | [CHANGELOG.md](CHANGELOG.md) |
| Code layout, rules, decisions, gotchas (old sections 5 to 8) | [docs/architecture.md](docs/architecture.md) |
| Original milestone specs and acceptance tests (old section 9) | [docs/milestone-specs.md](docs/milestone-specs.md) |
| What is left, shaky or missing | [ROADMAP.md](ROADMAP.md) |
| What the app does and how to use it | [docs/user-guide.md](docs/user-guide.md) |
| What has not been checked on real data | [docs/known-limitations.md](docs/known-limitations.md) |

Section numbers quoted elsewhere ("HANDOVER.md section 5") refer to the old single-file
numbering: 2 is the changelog, 5 to 8 are the architecture file, 9 is the milestone specs.
Sections 1, 3, 4 and 10 stay below.

## 1. What this is and who it's for

PyGCGC processes GC×GC-EI-ToF-MS data. It is a Python package (`pygcgc.core`, no GUI dependencies) plus a PySide6/pyqtgraph desktop app (`pygcgc.gui`).

The user is James Brean, an atmospheric chemist at the University of Birmingham.

**Instrument.** An older Markes/SepSolve GC×GC with a 70 eV EI time-of-flight MS (BenchTOF family, unit mass). Samples are thermally desorbed from Tenax tubes. ChromSpace is the vendor software, and data leave it as ANDI/AIA `.cdf`.

**Target compounds** (from James's calibration scheme; the dummy data mirrors it):
- n-alkanes, PAHs, nitro-PAHs and 2-octadecanone, from 100 µg/mL stocks
- oxy-PAHs, from 10 µg/mL stocks
- decalin, two cycloalkanes and behenic acid, from 1000 µg/mL stocks
- deuterated internal standards

The calibration has 8 levels, the first of which is an internal-standard-only blank.

**The internal standard and cycloalkane identities in the dummy data are placeholders.** James's answers (section 10, Q2-3, 2026-09-19): there are many internal standards and which ones apply depends on the sample/method, so a maintained list is needed, and the default IS-pairing rule should be the nearest internal standard by normalised (t1, t2) distance rather than a fixed per-target assignment (see section 7). "Cycloalkanes 1 and 2" are not two specific compounds -- they stand in for a large number of different cycloalkanes from C10 upwards, which doesn't matter for the software's fundamentals (any name/RI works for these rows). The actual internal-standard list (names, concentrations, which targets pair with which) is still needed before building a real `target_table.csv`.

## 3. Setting up

James works on Windows without admin rights, uses **uv with one central environment** shared across projects, writes code in VS Code, and keeps code and data in a OneDrive-synced folder.

```
uv pip install -e ".[gui,test]"      # with the central environment active
pytest                               # set QT_QPA_PLATFORM=offscreen if Qt complains
pygcgc-gui path\to\folder
```

Keep large data folders out of the repository. If they sit in OneDrive, sync may lock files mid-read; a local, unsynced copy is safer for heavy work.

The dummy data pack (`dummy_gcxgc/`) is found automatically if it sits beside or inside the repo, or via the `PYGCGC_DUMMY_DATA` environment variable. It contains:

| Item | Contents |
|---|---|
| `data/*.cdf` | 12 synthetic runs at 50 Hz, m/z 35 to 500, 4 s modulation, 5 to 36 min |
| `sequence.csv`, `target_table.csv`, `method.json` | the supporting files, in the layout the app expects |
| `calibration_levels_ng.csv` | amount of every target at every level |
| `synthetic_library.msp` | synthetic spectra with retention indices, plus 60 decoys |
| `truth_*.csv`, `hidden_truth.json` | ground truth: true phase 0.55 s, saturation ceiling 1.2e5, 86 spikes per run, and more |
| `generate_dummy_data.py` | the generator (fixed seed) |
| `README.md` | lists the deliberate problems |

## 4. First job with real data

Work through this in order and report findings to James before changing defaults.

1. **Summarise the files.** Run `pygcgc info` on a few runs and look at:
   - scan count and rate (expect tens of Hz)
   - start time and units (seconds are expected; the reader converts minutes if the `units` attribute says so)
   - the m/z range
   - the number of stored points
   - the timing report
   - any warnings (a high-resolution warning means masses are not unit)
2. **Check what ChromSpace actually exports.** The reader assumes the *unfolded* detector stream with all masses: one scan per acquisition time, centroids in `mass_values` and `intensity_values`, and `point_count` or `scan_index`. If the export is folded, TIC-only, profile-mode or split across files, the reader in `core/io/andi.py` needs adapting. Open a file with `netCDF4.Dataset` and list its variables and global attributes. Record anything unexpected in this file.
3. **Encoding details.** Check:
   - `mass_values` dtype and any `scale_factor` (auto-applied)
   - whether intensities are integers
   - whether stored values have a clear lower threshold (`core/background.storage_threshold`)
   - whether any stored values are negative or NaN (the reader drops them)
4. **Modulation period.** Get it from James or the method (it isn't in the ANDI file) and put it in `method.json` next to the data.
5. **Make a real target table and levels file.** Open a run in the app, find the internal standards, and note their apex t₁ (absolute seconds) and t₂. Build `target_table.csv` with `role=istd` rows. Then run `pygcgc phase ... --targets target_table.csv` across a sequence. Check that:
   - anchors are found in most runs
   - the phase spread is a few ms
   - the period check is "not significant"

   A significant period correction on real data is interesting: ask James before applying it.

   James plans to run standards of known composition and concentration as calibration sequences. Build `calibration_levels.csv` from his dilution scheme (one row per level, one column per compound, amounts in the unit he wants results in), and add `cal_set` and `volume_m3` columns to `sequence.csv` if there is more than one calibration block or if samples have air volumes.
6. **Visual check against ChromSpace.** Run `pygcgc fold RUN.cdf --period P --targets ... --extend 1.15 --png check.png` and compare with ChromSpace's own 2D plot: orientation, t₂ range, wrap-around and the positions of known compounds.
7. **Background on real data.** In the app, apply 🧹 and look at the diagnostics for bleed ions (73, 207, 281, 355) and air ions (40, 44 if in range). The residuals should be centred with SD close to the noise map. Try the split view on a big tailing peak and on a steep bleed rise. If the defaults (`BackgroundSettings`: percentile 20, stiffness 1e4, pattern window 51) need changing, change them in one place and update the tests and README.
8. **First real quantification.** Run 🧪 Quantify (or `pygcgc quant FOLDER`) and check, in this order:
   - every internal standard is found in every run (Results → Internal standards; the trend chart should be flat within the tolerance)
   - learned positions look sensible (the overlay circles sit on the peaks)
   - the Calibration tab: which curves fail, and why (look at the relative-error plot before touching the weights)
   - the Review tab with Only to check ticked: which flags are common, and whether the defaults in `QuantSettings` suit real peak shapes (tolerances, S/N minima, `wrap_below_s`)
   - the detector ceiling found from the data (`detect_ceiling`) against the instrument's known saturation level
9. **Real regression fixtures.** Use `pygcgc crop RUN.cdf --start A --end B --mz-min X --mz-max Y -o small.cdf` to make a few small files: a few minutes, including an internal standard and some bleed.
   - **Ask James before committing any real data.** If he'd rather not, keep them outside the repo, point a `PYGCGC_REAL_DATA` variable at them, and write tests that skip when it is missing, following the `dummy_dir` fixture in `tests/conftest.py`.
   - Test at least: the reader's invariants, a stable phase estimate, and centred background residuals.
10. **Native `.lsc` files:** James explicitly requested direct instrument-file
    support on 2026-09-19. The standalone investigation now decodes this FT04/SD04
    variant successfully. Follow `real_data/NATIVE_LSC_REPORT.md` for integration,
    remaining validation and acquisition/channel/DBC handling. Native import is
    now registered alongside ANDI; independent vendor validation is still pending.

### Assumptions to verify, with where they live

| Assumption | Code |
|---|---|
| Times in seconds, strictly increasing (repeats are dropped with a warning) | `core/io/andi.py` |
| Unit-mass centroids; duplicates within a scan are summed after rounding | `core/io/andi.py` |
| One modulation period for the whole run | `core/fold.py` |
| Acquisition start (`t0`) is the reference for the phase | `core/fold.py` |
| Scan rate is steady apart from dropped scans; jitter under 10% of the interval is fine | `core/timing.py` |
| Internal standard quantifier ions are distinctive enough that the brightest pulse in a ±12 s window is the standard | `core/phase.locate_anchor` |
| Values below the storage threshold are absent; noise floor = threshold / 3 | `core/background.py` |
| Bleed is smooth along t₁ over about ten modulations | `core/background.py` |
| Nothing is retained in the first 0.5 s of t₂, so peaks there have wrapped (`PeakSettings.wrap_below_s`) | `core/peaks.py` |
| The brightest stored value that occurs at least 3 times is the detector ceiling | `core/peaks.detect_ceiling` |
| A calibration run has "cal" in its type and a level; everything else is a sample or blank | `core/quant.SequenceRun` |
| Unknown peaks respond like their nearest calibrated neighbour (semi-quant only) | `core/quant.all_signals` |

## 10. Questions to ask James

1. What is the modulation period, and does ChromSpace export the raw unfolded stream with all masses?
   **James, 2026-09-19:** The period varies from run to run, usually 2 to 13 s, and is never stored in the datafiles, so it always has to come from the method (`method.json`); do not expect to recover it from the file. Still open: whether the ChromSpace export is the raw unfolded all-mass stream. That needs a real file (section 4).
2. Which internal standards are used, at what concentrations, and which targets pair with which standard?
   **James, 2026-09-19:** Many internal standards, and which apply depends on the sample. The identities are always known beforehand, so keep a list of them. Default pairing: normalise the two retention times and use the internal standard at minimum distance (matches the All signals rule in section 7). Still needed from James: the actual list, with concentrations.
3. Which cycloalkanes are the "cycloalkanes 1 and 2" in the standards?
   **James, 2026-09-19:** They are not two specific compounds. They stand for a large set of cycloalkanes from C10 upwards, which doesn't matter for the fundamentals of the software. Closed.
4. What air volumes are sampled per tube (for concentrations)?
   **James, 2026-09-19:** Varies widely, from a few cc/min for 30 min to several days. `volume_m3` must stay a per-sample column with a wide range. Closed.
5. Is there a NIST (or Wiley) licence for library export yet? (As of milestone 9, not yet — James's own words were "I'll add the NIST library later" — so 🏷️ identification was built and tested against small synthetic libraries only. The moment a real library (`lib2nist`-exported MSP, MoNA, MassBank, or otherwise) exists, load it through 🏷️'s "Load library" and see what actually breaks or needs retuning before trusting any of its default thresholds.)
   **James, 2026-09-19:** Still no library (checked 2026-09-19).
6. May small cropped real files be committed to the repository for tests?
   **James, 2026-09-19:** "Yeah I guess?" A tentative yes. Confirm again before actually committing any file. James also wants tests run against his real data rather than dummy data (Q11, Q14).
7. Should the detector saturation level be taken from the instrument manual, or estimated from the data?
   **James, 2026-09-19:** Estimate it from the data (as `detect_ceiling` does now). He asks whether it can also be got from the peak shape (a flat-topped or clipped apex). Worth trying once real saturated peaks exist; listed in ROADMAP.md.
8. For alignment: which run makes a sensible default reference (a mid-sequence calibration run, or something else), and is a single reference for the whole project enough, or does drift across a long campaign call for re-referencing partway through?
   **James, 2026-09-19:** The reference-run choice is unanswered. The constraint is that drift is sub-second run to run while peaks are only about 0.25 s wide, so alignment precision matters a great deal, and the case for re-referencing across a long campaign should be tested on real data.
9. For 🔍 blob detection on real data: are the default S/N thresholds (`snr_min` 10, `floor_snr` 2) about right, or does a busier real chromatogram need them raised? And does the co-elution flag actually fire sensibly on a known co-eluting pair (e.g. phenanthrene/anthracene), or does it need retuning (`coelution_frac`, `coelution_gap`)?
   **James, 2026-09-19:** Don't know; test on real data.
10. For 🫧 siloxane removal: is there a real reference spectrum for the siloxanes actually seen on this column (a measured Dn standard, or even a representative real bleed spectrum), to replace the generic, reasoned-not-measured `SILOXANE_REFERENCE`? And once real (or the dummy) data are to hand, `tests/test_siloxane.py::test_dummy_siloxane_score_separates_real_targets_from_planted_siloxanes` needs running to actually choose `score_threshold` from a ROC curve, as the milestone always intended, rather than the provisional 0.3.
   **James, 2026-09-19:** No siloxane standard exists. Siloxanes come from many sources and look different each time. There are a few characteristic peaks to look for, but he can't recall the m/z values (the generic list in the milestone 7 spec, 73, 207, 281, 355 and so on, is the starting point). Choose the threshold from real bleed once seen.
11. For ⛓️ the hydrocarbon mapper: what real diagnostic ion pairs does James actually want for the homologous series he sees (the spec named 91/105, 82/83, 142/156, 58/59, but `find_homologous_series` takes any pair, untested against real data)? Is pristane/phytane's assumed position (a fixed carbon-equivalent offset before nC17/nC18) close enough to be useful, or does it need a real measured offset for this column? And is skipping a hard qualifier-ratio gate on the alkane finder (57/71/85 reported but not enforced, since real fragmentation ratios vary and there was nothing but synthetic data to calibrate a threshold against) the right call, or does real data show it needs tightening?
   **James, 2026-09-19:** Test against his real data.
12. For 🗂️ batch: is the fixed feature-matching radius (`max_distance`, in scaled t1/t2 units) and gap-fill window (`t1_half_s`/`t2_half_s` in `_box_sum`) about right for real peak widths, or do they need to scale with something (peak width from 🔍, say)? What does James actually want to compare in a real volcano plot (two sample groups, a treatment vs control, a time series?) — the two-group Welch's-t-test shape was a reasonable default, not something he'd specified. And is cores − 1 the right default worker count on his machine, or does real batch size and file size on a OneDrive-synced folder call for something more conservative?
   **James, 2026-09-19:** Not sure of the matching radius; check on real data. Comparisons he wants to make: pairs of similar samples, samples against blanks, and samples against mixes of standards (so the batch comparison should not assume a single treatment-versus-control design). Worker count unanswered.
13. For 📤 reports and the installer: does James actually want a standalone PyGCGC.exe he can hand to a colleague without a Python environment (worth someone finally installing PyInstaller and building `packaging/pygcgc.spec` to find out what it's missing), or is `uv pip install -e .` on a shared environment good enough for everyone who'd use this? And for reports specifically -- is a report ever generated from something other than a fresh batch run (say, from a project someone has been reviewing by hand for a while), which is the case that would actually need `build_provenance`'s manual-edit-log arguments wired up to something real, rather than always empty as they are now?
   **James, 2026-09-19:** No standalone .exe for now; handing colleagues the Python package is enough. Report-source question unanswered.
14. For ✂️ overlap: BbF/BkF and the benz[a]anthracene/chrysene pair (both in the dummy data) are the obvious first real test once either the dummy pack or real data are available -- does `fit_overlap`/`choose_k` actually separate them within the spec's original 10% target, and does `mcr_als` agree? Is a joint peak-shape fit (parametric, needs a shape assumption) or MCR-ALS (free-form, needs non-negativity and, optionally, unimodality) the one James would actually trust more day to day, or does it depend on the pair? Is PARAFAC across runs (not yet reachable from the GUI -- see section 9) something he'd use often enough to be worth a proper "pick a feature, pick some runs" GUI flow, or is calling it from Python rare enough that it's not worth building yet?
   **James, 2026-09-19:** Use real data, not the dummy data, for this test. Method preference and PARAFAC GUI question unanswered.
15. For 🏷️ identification, once a real library exists: does `search`'s prefilter (RI window plus a shared top ion) actually keep enough true candidates while cutting the rest, or does a real library's spectral diversity need a wider or narrower window? Are the score weights (`forward_weight`/`reverse_weight`/`ratio_weight`/`ri_penalty_weight`) and the tier thresholds (`confident_score`, `ambiguity`) anywhere close to sensible, or do they need a genuine ROC-style tune the way 🫧's `score_threshold` still needs? And does phenanthrene/anthracene (or BbF/BkF/BaP) actually get separated correctly by the RI penalty on a real alkane ladder, the spec's own named acceptance case -- the synthetic test proves the mechanism works, not that the default `ri_window`/`ri_penalty_full` are the right size for a real column's real RI spread.
   **James, 2026-09-19:** No idea; needs a real library and real data.
