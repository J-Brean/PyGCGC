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
  (about 4 minutes; run it in the background). Without offscreen, GUI tests can hang on a dialog.
  28 tests skip without the dummy data pack. Every bug fixed gets a test that failed before.
- **Real data** lives in `real_data/` and is never committed. Benchmarks read it; tests skip
  without it. Ask before changing a scientific default; report the evidence first.
- **Style.** British English, no em dashes, plain-English docstrings for anything methodological,
  scripts that read top to bottom. Figures follow the user's global figure rules.

## Where things are

| | |
|---|---|
| `src/pygcgc/core/` | all the science, no Qt; one module per method, each opening with an explanation |
| `src/pygcgc/gui/` | PySide6/pyqtgraph app; `mainwindow.py` wires everything, `*_state.py` hold per-step state |
| `src/pygcgc/presets/` | shipped compound-class rules and internal-standard library |
| `tests/` | pytest; synthetic fixtures in `conftest.py` |
| `benchmarks/` | method benchmarks and real-data checks (see its README) |
| `examples/tutorials/` | walkthroughs on the dummy data pack |
| `docs/` | user guide, architecture, validation reports, history (superseded plans and handovers) |

## State on 24 September 2026

- Branch `audit-2026-09-24` holds a full audit: a snapshot commit of the previously uncommitted
  tree, bug fixes, and this reorganisation. Its findings are in
  [docs/validation/audit-2026-09-24.md](docs/validation/audit-2026-09-24.md) and ROADMAP.md.
- `real_data/` now holds four tandem (70 and 14 eV) acquisitions at an 8 s modulation (two
  Repsol exhaust samples, ATWC and BTWC, and two diesel samples, D50T191 and D50T60). The
  earlier 10 s acquisition that most older reports describe has gone; only its cached crop
  (`real_data/blob_diagnosis/segment.npz`) remains.
- Five of the seven library internal standards are in every run and are found by both 🏷️ and
  🧪 (draft table: `benchmarks/real_data_audit/internal_standards.py`).
- No real calibration sequence or spectral library exists yet.
