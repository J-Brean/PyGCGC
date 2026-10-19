# Documentation

Current guidance was reconciled with the code on 19 October 2026. Start with
[installation and quick start](../README.md), then follow **Load → Align (optional)
→ Background → Peaks → Standards → Assign → Quantify → Export**.

## Current reference

| Document | Use it for |
|---|---|
| [User guide](user-guide.md) | workflow, controls, specialist tools, persistence and methods |
| [Input tables](input-tables.md) | CSV columns, units, validation and local starter tables |
| [Demonstration data](../examples/demo_data/README.md) | a reproducible synthetic tour with known amounts |
| [Architecture](architecture.md) | state ownership, jobs, fitting defaults, caches and development conventions |
| [Roadmap](../ROADMAP.md) | outstanding validation and implementation gaps |
| [Benchmarks](../benchmarks/README.md) | reproduction commands and dataset requirements |
| [Agent notes](../CLAUDE.md) | workspace rules and current implementation constraints |
| [Changelog](../CHANGELOG.md) | dated changes and the checks performed at each revision |

## Latest evidence

| Report | What it establishes |
|---|---|
| [Code audit, 19 October](validation/code-audit-2026-10-19.md) | 850 tests passed, 31 optional-data skips; numerical, import and job-state fixes; native and wheel smoke checks |
| [Anchor workflow, 19 October](validation/anchor-workflow-2026-10-19.md) | guided anchor collection, editing, restart, bounded refinement and four-run native smoke checks |
| [Guided workflow, 18 October](validation/guided-workflow-2026-10-18.md) | calibration setup, standard exceptions, prerequisites and fitting previews |
| [Fitting comparison, 18 October](validation/peak-fitting-2026-10-18.md) | eight methods on 32 real windows and known-answer recovery in measured backgrounds; tilted Gaussian default |
| [Initial workflow audit, 18 October](validation/workflow-and-peaks-2026-10-18.md) | workflow/export changes; its initial detector findings are superseded by the fitting comparison |

Read each report's date, settings and validation limits. Earlier test totals describe
earlier revisions. A passing software test, spectral similarity or filename does not
establish chemical identity or an independently known amount. Re-run Quantify after
the 19 October numerical calibration fix to regenerate amounts and uncertainties.

## Earlier validation

These retain the evidence and assumptions of their original datasets. They are not
current instructions or proof that later versions reproduce identical numbers.

| Report | Original scope |
|---|---|
| [Calibration standards, 6 October](validation/calibration-standards-2026-10-06.md) | four August 10 s calibration acquisitions; suggested compounds and provisional amounts/IS identities |
| [Real-data audit, 24 September](validation/audit-2026-09-24.md) | four 8 s tandem exhaust/diesel acquisitions; phase, ladder, siloxanes and planted quantitative recovery |
| [Classes, 21 September](validation/classes-2026-09-21.md) | classification on the earlier 10 s acquisition |
| [Hydrocarbons, 20 September](validation/hydrocarbons-2026-09-20.md) | ladder/marker/reference checks on the earlier 10 s acquisition |
| [Native import](validation/native-import.md) | implementation and scan-level checks for the original native acquisition |
| [Native format investigation](validation/native-lsc-format.md) | forensic format notes and the original integration proposal |

Private raw data now lives under `C:\REAL_DATA_FOR_PYGCGC`; the main calibration
resource is `Cal files for James`. Generated outputs and legacy cached crops remain
under the ignored repository `real_data/`. An old script's hard-coded path or fold
setting must be checked against its original dataset. See the benchmark guide.

## Historical design records

[History](history/README.md) contains superseded plans and handovers. It preserves
past reasoning, limitations and proposed work. Current behaviour is described by
the reference documents above. Local figures not committed with a historical report
are identified as local artefacts rather than displayed as broken images.
