# PyGCGC

A desktop app and Python package for comprehensive two-dimensional gas chromatography with
electron-ionisation time-of-flight mass spectrometry (GC×GC-EI-ToF-MS). It reads native Markes
FT04/SD04 `.lsc` files (including tandem 70/14 eV acquisitions) and ANDI (AIA) `.cdf` exports.

![The main window: a synthetic diesel-exhaust run after background removal](docs/images/screenshot_main.png)

| | | |
|---|---|---|
| ![Peaks](docs/images/screenshot_peaks.png) MS-aware peaks | ![Hydrocarbons](docs/images/screenshot_hydrocarbons.png) n-alkane ladder | ![Internal standards](docs/images/screenshot_istd.png) internal standards |
| ![Calibration](docs/images/screenshot_calibration.png) calibration curves | ![Classes](docs/images/screenshot_classify.png) compound classes | ![Clusters](docs/images/screenshot_cluster.png) spectral clusters |
| ![Calibrate window](docs/images/screenshot_calibrate_curves.png) 📈 calibration builder | ![Compare](docs/images/screenshot_align_overlay.png) alignment overlay | ![Batch](docs/images/screenshot_batch.png) batch PCA |

Every screenshot is of the synthetic [demonstration data](examples/demo_data/README.md) and is
described in the [user guide](docs/user-guide.md).

## Status

Every step below is built and tested on synthetic data with a known answer (547 tests). The
right-hand column says what has been checked on real runs; the details are in
[docs/validation/audit-2026-09-24.md](docs/validation/audit-2026-09-24.md) and what is left is in
[ROADMAP.md](ROADMAP.md).

| Step | What | On real data |
|---|---|---|
| Read | native `.lsc` (70 eV, 14 eV, vendor DBC) and ANDI | Decodes every scan of the supplied files in 1 to 3 s each; vendor comparison of spectra and times still pending |
| 1 🔄 Fold | period and phase, from the run or from anchors | Period found unaided on four runs (8.0000 s, as set in the method); anchor-free phase within 0.2 s of the hand-set value |
| 2 📍 Align | anchors, shift/affine/smooth transforms, compare modes | Not yet checked on real drift |
| 3 🧹 Background | bleed, stripes, spikes, noise map, blank subtraction | Used on every real run; blank subtraction unchecked |
| 4 🔍 Peaks | MS-aware overlapping components (default) or legacy blobs | Benchmarked with planted peaks in real background; no labelled real standard yet |
| 4a 🫧 Siloxane filterer (in Peaks) | score blobs against bleed and cyclic references; block their siloxane ions, or fit per pixel | Defaults chosen on two runs: every siloxane-like blob caught, none of the hydrocarbon signal removed |
| 4b ⛓️ Hydrocarbon identifier (in Peaks) | n-alkane ladder (straight, smooth or polynomial), peak assignment, stencils, homologous series | Ladder from two clicks reaches C9–C30 on four runs; one or two carbons wrong per run. Assignment and ladder shapes not yet checked on real data |
| 5 🏷️ IS | label internal standards, find them on every run, write them to the target table | Five deuterated standards found with high confidence on all four runs |
| 6 🧪 Quantify | targets, internal standards, weighted calibration, bracketing | Planted standards in real matrix: within a few % on selective ions, 15–60% off on busy ones |
| 7 🧬 Classify | compound classes from readable rules | Checked on one earlier run |
| 8 🏷️ Identify | MSP library search with target-decoy confidence | No real library yet |
| 9 🧩 Cluster | peaks or pixels grouped by spectrum, within and across runs | Checked on one earlier run; not on a sample set |
| 📈 Calibrate | separate window: standard runs, compounds from a table or clicks, IS-normalised blobs, curves edited point by point | Four real levels without a compound list: 59 calibrants and four internal standards suggested; phases lined up; responses bend at 25 ng |
| 🗂️ Batch, 📤 Reports | cached parallel pipeline, feature table, QC, PCA, volcano, Excel/HTML/PDF | Ran end to end on nine real-matrix runs |

## Install

Python 3.11 or newer. From the repository folder, in PowerShell with your environment active:

```powershell
uv pip install -e ".[gui,test]"
```

This installs two commands, `pygcgc-gui` and `pygcgc`. The core alone (`uv pip install -e .`)
needs only numpy, scipy, scikit-image and netCDF4; the extras are `plot` (matplotlib),
`report` (openpyxl), `cluster` (scikit-learn), `gui` and `test`.

## Quick start

```powershell
pygcgc-gui "C:\path\to\folder\of\runs"
```

Work along the numbered toolbar: open a step, set it on the left, press **▶ Go**, then follow the
banner to the next step. A first look needs only 1 🔄 Fold and 3 🧹 Background. Set each run's
type and calibration level in the Runs table (or with a `sequence.csv`). The project saves itself
beside the runs (`project.pygcgc.json`) and offers to resume next time; background layers and
peak detections are cached on local disk, so reopening does not recompute them. The
[user guide](docs/user-guide.md) covers every step.

For native files, confirm the modulation period in the Fold step (it is not stored in the file);
folder opening picks untreated 70 eV, and **Open files…** loads 14 eV or DBC channels.

## Demonstration data

Seven synthetic runs (four calibration standards, three samples) with about 150 compounds,
impurities, artefacts and the true answer for every number:

```powershell
python examples/demo_data/generate_demo_data.py demo_gcxgc
pygcgc-gui demo_gcxgc
```

[examples/demo_data/README.md](examples/demo_data/README.md) lists what is in the pack and gives
a tour of every step, including the 📈 calibration builder. The tutorials in
`examples/tutorials/` run on the same pack.

## Supporting files

All optional for viewing and cleaning; picked up when they sit in the folder you open.
Quantification needs a target table, a sequence and calibration levels.

**`method.json`**: `{"modulation_period_s": 8.0}`.

**`sequence.csv`**: what each run is, in running order (or set it in the Runs table and save it from there).

```
file,type,level,cal_set,volume_m3
CAL_L1_mon.cdf,calibration,L1,monday,
SAMPLE_A.cdf,sample,,,1.2
TUBE_BLANK.cdf,tube blank,,,
```

A type containing "cal" plus a `level` is a calibration run; a type containing "blank" is offered
as a blank. `cal_set` groups calibration runs; a sample between two sets is read against both,
blended by position (turn this off in the quantification settings). `volume_m3` turns amounts into
concentrations.

**`target_table.csv`**: the compounds in your standards; `role` = `istd` marks internal standards,
which also anchor the phase. The 🏷️ IS step can write these rows for you.

```
compound,role,quant_ion,qual_ions,qual_ratios,t1_expected_s,t1_window_s,t2_expected_s,t2_window_s,istd
phenanthrene-d10,istd,188,184;160,,2708,12,4.60,0.25,
phenanthrene,target,178,176;152,,2715,12,4.55,0.25,phenanthrene-d10
```

`t1_expected_s` is the absolute first-dimension time; `t2_expected_s` the true second-dimension
time, which may exceed the period for a compound that wraps. Positions are refined from the
calibration runs, so rough values do. Leave `qual_ratios` blank to learn them from the
calibration runs. Leave `istd` blank to use the nearest internal standard.

**`calibration_levels.csv`**: the amount of every compound at every level, in the unit you want
results in.

```
level,phenanthrene,fluorene
L1,0.1,0.1
L2,0.25,0.25
```

## Command line

```
pygcgc info  RUN [RUN ...]
pygcgc phase RUN [...] [--period 8] [--targets CSV]
pygcgc fold  RUN [--period 8] [--phase 2.77 | --targets CSV] [--mz 57] [--extend 1.15] [--png out.png]
pygcgc clean RUN [...] [--period 8] [--outdir cleaned]
pygcgc crop  RUN --start 10 --end 14 [--mz-min 35 --mz-max 300] [-o small.cdf]
pygcgc quant FOLDER [--outdir results] [--xlsx] [--no-interpolate]
pygcgc batch FOLDER [--siloxane] [--no-blobs] [--no-quant] [--cluster] [--regions N] [--report all]
pygcgc gui   [FOLDER | PROJECT.json]
```

Without `--period` (and no `method.json`) the period is estimated from the runs and printed.
`quant` writes targets, calibrations, all-signals and internal-standard tables plus
`provenance.json`. `batch` runs the same pipeline across a process pool with a disk cache, adds a
feature table, and optionally clusters and reports. Alignment, the alkane ladder, overlap and
identification are interactive and have no command.

## Python

```python
from pathlib import Path
import pygcgc as g

run = g.read_run("real_data/sample_70eV.lsc")
f = g.fold(run, 8.0, 2.77)
clean = g.remove_background(run, f)
blobs = g.find_blobs(clean.run, f)
ladder = g.find_alkane_series(clean.run, f, [(15, 1953.0, 1.13), (16, 2188.0, 1.12)])
print(ladder.carbon_of(2500.0))          # retention as a carbon number
```

Every module in `pygcgc.core` opens with a plain-English account of its method. The tutorials in
`examples/tutorials/` walk through a whole sequence on the demonstration data pack: reading and
folding, background, quantification, alignment, peaks, siloxanes and the calibration builder.

## Tests

```powershell
$env:QT_QPA_PLATFORM = "offscreen"; python -m pytest tests -q
```

About four minutes. 28 tests need the synthetic dummy data pack (`dummy_gcxgc/` beside the
repository, or `PYGCGC_DUMMY_DATA`) and skip without it; one needs a native run in `real_data/`.

## Documentation

| File | What it holds |
|---|---|
| [docs/user-guide.md](docs/user-guide.md) | every step and panel, and how the methods work |
| [examples/demo_data/](examples/demo_data/README.md) | the synthetic demonstration pack, its generator and a tour of the app; `docs/images/make_screenshots.py` retakes every screenshot from it |
| [ROADMAP.md](ROADMAP.md) | known limits, recommended changes and suggested features, in priority order |
| [CHANGELOG.md](CHANGELOG.md) | dated record of what changed and what real-data checks found |
| [docs/architecture.md](docs/architecture.md) | code layout, rules, design decisions, gotchas |
| [docs/validation/](docs/validation/) | real-data reports, newest first: [audit-2026-09-24](docs/validation/audit-2026-09-24.md) |
| [benchmarks/](benchmarks/README.md) | scripts that reproduce every benchmark and real-data number |
| [CLAUDE.md](CLAUDE.md) | notes for coding agents working on the repository |
| [docs/history/](docs/history/) | superseded plans, specs and handovers, kept for reference |
