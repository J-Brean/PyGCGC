# PyGCGC: notes for coding agents

Read this, then README.md and ROADMAP.md. The user is James Brean (atmospheric chemistry,
University of Birmingham). Code conventions are in [docs/architecture.md](docs/architecture.md).

## Working rules

- **Version control.** Work on a branch, commit in small, described steps, and never `git stash`,
  `reset` or `checkout --` over work you did not write. Other sessions (Claude and Codex) edit
  this tree; after any pause, check `git status` and re-read files before editing them.
- **Environment.** Windows, one shared uv environment at `C:\Users\BreanJ\venvs\main`. Run Python,
  pytest and git from PowerShell; Git Bash here lacks git and most tools.
- **Tests.** `python -m pytest tests -q -p no:cacheprovider` with `QT_QPA_PLATFORM=offscreen`
  (several minutes; run it in the background). Without offscreen, GUI tests can hang on a dialog.
  Optional dummy/native/reference tests skip without configured inputs. Every bug fixed gets
  a test that failed before. Native GUI tests can autosave: use isolated project paths and
  preserve existing private projects when adding real-data coverage.
- **Real data** lives in `C:\REAL_DATA_FOR_PYGCGC`; `Cal files for James` is the main calibration
  resource. Treat its raw files and existing project as read only. Generated audit results go in
  the repository's ignored `real_data/`. Opt-in calibration tests use `PYGCGC_CALIBRATION_DATA`;
  older scripts may still refer to the previous data location. Ask before changing a scientific
  default; report the evidence first.
- **Results cache.** 🧹 layers and 🔍 detections are cached in `%LOCALAPPDATA%\PyGCGC\results`
  (tests use a temporary folder via `PYGCGC_CACHE`), keyed by `core.batch.CODE_FINGERPRINT`, a hash
  of the modules whose results are cached. A new module that feeds a cached result must be added
  to that list, or stale results will be served.
- **Line endings.** The repository stores LF; some editors here write CRLF. Scripted
  find-and-replace on `\n` silently misses CRLF files, so check that an edit applied.
- **Style.** British English, no em dashes, plain-English docstrings for anything methodological,
  scripts that read top to bottom. Figures follow the user's global figure rules.

## Where things are

The [19 October code audit](docs/validation/code-audit-2026-10-19.md) adds deterministic
queued-cancellation, obsolete-result, transactional-import and calibration-unit tests.
Workers must reject changed source/settings before publishing and settle running flags
on cancellation. Calibration must solve the scaled weighted design directly, retaining
unit invariance; do not form normal equations. Target-table appends preserve arbitrary
valid header order and extra columns, reject duplicates, and replace files atomically.

The 20 October feedback round added four things; read the notes in docs/architecture.md section 5 before
touching them. **Text size and screens:** `theme.DENSITY`, `prefs.py`, `fit_to_screen`; tests point
`PYGCGC_SETTINGS` at a temporary INI. **Anchoring** is a window (`gui/anchor_window.py`) over the main
window's unchanged picking machinery; clicks are held (`AlignState.holding`) until ANCHOR, and picking needs
every run's phase ticked. **📚 Library matcher** (`core/nist_sdf.py`, `core/libmatch.py`, `gui/match_state.py`,
`gui/library_matcher.py`): the NIST files are licensed to James, so never commit, copy or ship them or anything
parsed from them (tests use the invented files in `tests/nist_fixture.py`; `PYGCGC_NIST_SDF` and
`PYGCGC_LIBRARY_CACHE` are for the opt-in real checks); the 5 % phenyl column class and the 850/750/650
thresholds are James's choices of 20 October and still provisional. **➖ Peak-level blank subtraction**
(`core/peak_blank.py`, `QuantSettings.blank_*`): James decided quantified amounts must be sample minus factor times
blank. `C:\REAL_DATA_FOR_PYGCGC\background_removal` holds one blank-like run (9003, Dilutor) and one sample (9002);
the autosave writes beside the data, so scripts must switch `Project._autosave` off.

The user authorised **tilted Gaussian as the default** on 18 October 2026, after the
eight-method comparison. `ComponentDetectionSettings` uses `model="tilt"`, independent widths
and automatic component-count pursuit. Both Gaussian models must retain pursuit and window
merging. `defaults_version=2` distinguishes new deliberate model choices from the former
automatic default on project reload. Existing results are not rewritten on load. Accepted
manual fits are reapplied after detection when their saved region geometry still matches;
they need an explicit manual refit to adopt tilt.
Manual refits also start with tilt; low-level fitting APIs keep their Gaussian defaults for
script compatibility. The local, scripted 1D/2D explanation is `gui/fitting_help.py`; it must
not start fitting jobs or leave a timer running when hidden. Its 1D chromatogram must be
an unfolded detector pulse train, never a summed first-dimension projection.

The four workflow priorities are now implemented. `gui/calibration_setup.py` is a draft-only
wizard until Apply; `core/calibration_inputs.py` validates amounts and treats ambiguous filename
tokens as unknown. Constant IS additions are explicit; never fill the tentative private 1 ng
as fact. Guided targets/levels/unit/confirmation live in `QuantState` and project JSON, without
overwriting input CSVs. Legacy external-input loading remains supported.
`gui/standard_review.py` uses per-run `IstdState.review_decisions`, with reproducible input/result
signatures. Acknowledgement retains fit warnings. Reviewed choices guide quantitative positions;
confirmed misses cannot normalise targets. `gui/workflow_assistant.py` derives prerequisites
and freshness; `gui/peak_preview.py` runs cancellable production fitting without storing preview
components. The preview can explicitly apply settings and start the normal all-run job.
Setup and readiness use the actual calibration minimum (three by default), zero-level policy
and curve exclusions. Do not lower scientific requirements to make the wizard pass.
Fusion style is selected once per QApplication; colour switches must not recreate it or reset
an equivalent resolved font. Candidate-plot redraws after palette changes are deferred until
Qt repolishing finishes. Visual QA reproduced a native Windows crash without these guards.

Alignment's primary action now collects one anchor across all loaded runs, returns to the
reference and prompts for the next. Reference clicks are stored as identity AnchorPoints
in `AlignState.points[reference]`; legacy projects with target-only pairs still load.
The table shows all runs. Esc ends picking but keeps clicks; missing peaks can be skipped.
Picking must use unaligned images and reject clicks on a previous run's pending render.
`core.align.refine_anchor_position` searches a bounded intensity maximum, not a centroid;
defaults are ±1 modulation and ±0.150 s in t2. Wrapped rows come from adjacent columns.
`AlignState.refine_all` runs on saved ions/TIC and the current layer in a cancellable job,
updates each reference position once, and rejects results after input edits. Apply as one
undoable edit. Main picking/refinement stay visible; search limits and transform diagnostics
are under More alignment options. The panel scrolls rather than squeezing the anchor table.
Add points, Delete selected points and Start alignment over are visible controls, explicitly
requested by the user. Delete selected rows in one undo command; deleting a reference row
removes its matching pairs. Restart retains the reference but clears points and kind overrides
and resumes P1. Undo must stop that new picker before restored P1 can be overwritten.

| | |
|---|---|
| `src/pygcgc/core/` | all the science, no Qt; one module per method, each opening with an explanation |
| `src/pygcgc/gui/` | PySide6/pyqtgraph app; `mainwindow.py` wires everything, `*_state.py` hold per-step state |
| `src/pygcgc/presets/` | shipped compound-class rules and internal-standard library |
| `tests/` | pytest; synthetic fixtures in `conftest.py` |
| `benchmarks/` | method benchmarks and real-data checks (see its README) |
| `examples/tutorials/` | walkthroughs on the dummy data pack |
| `docs/` | user guide, architecture, validation reports, history (superseded plans and handovers) |

## Current data and documentation

Private instrument data is under `C:\REAL_DATA_FOR_PYGCGC`, including the August 2026
10 s calibration acquisitions in `Cal files for James` and earlier 8 s exhaust/diesel
acquisitions. The repository's ignored `real_data/` holds generated outputs and legacy
crops/caches. Older scripts may hard-code that former raw-data location; check their
paths, fold settings and reference requirements before running them.

The supplied `AEROSOLS/Standards.xlsx` inventory is available. Its RI column contains
catalogue codes. Acenaphthene-d10, phenanthrene-d10 and n-tetracosane-d50 have inventory-
matched calibration starter tables. The m/z 240 deuterated PAH is isomer-unconfirmed;
do not perpetuate the former heuristic chrysene-d12 label as established identity.
Exact per-compound amounts, ambiguous filename amounts and IS additions are unconfirmed.

Use [docs/README.md](docs/README.md) as the documentation map. Keep current instructions
in README, user guide, architecture and ROADMAP. Changelog and dated validation results
record their original revision; add a dated correction or supersession note instead of
silently replacing historical measurements. Preserve archived plans in `docs/history/`.
Update UI labels and step numbers against the code, check relative links, and distinguish
the main Export action from batch reports and a Python wheel from a Windows executable.
