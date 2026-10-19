# Input tables

Use **Quantify → Set up calibration and units…** for a guided setup, or import CSV
files through the advanced controls. CSVs use a header row and UTF-8 encoding; a UTF-8
BOM is accepted. Retention times and search half-widths are in seconds.

## Target table

| Column | Meaning and default |
|---|---|
| `compound` | required, unique compound name; unnamed rows are ignored |
| `quant_ion` | required, positive whole-number quantifier m/z; `188.0` is accepted |
| `t1_expected_s` | required, finite absolute first-dimension retention time |
| `t2_expected_s` | required, finite true second-dimension time, including values beyond the modulation period for wrapped peaks |
| `role` | `target` by default; `istd` marks an internal standard and phase anchor |
| `cls` | optional compound-class label |
| `qual_ions` | optional qualifier m/z values, separated by semicolons |
| `qual_ratios` | optional non-negative qualifier/quantifier ratios, in the same order and number as the qualifier ions; blank ratios can be learned from calibration runs |
| `t1_window_s` | positive search half-width, default 12 s |
| `t2_window_s` | positive search half-width, default 0.25 s |
| `istd` | name of an `istd` row; blank uses the nearest usable internal standard |
| `stock_factor` | positive relative stock multiplier, default 1; not an amount or an IS spike |
| `ri` | optional positive chromatographic retention index; never a catalogue identifier |

Numeric values must be finite. Duplicate headers or compound names, unknown roles,
fractional ion masses, mismatched qualifier ratios and invalid internal-standard
references are rejected. A failed Project target import preserves its working
inputs. Importing a valid table replaces the active target list. A target-only
table can be imported and reopened; quantification still requires internal standards.
GUI target imports preserve the prepared phase. Use **Load → Check modulation period
and phase** to re-estimate explicitly. Opening a data folder with an IS table can
use its anchors during initial folding; the Project API also offers explicit estimation.

```csv
compound,role,quant_ion,qual_ions,qual_ratios,t1_expected_s,t1_window_s,t2_expected_s,t2_window_s,istd
phenanthrene-d10,istd,188,184;160,,2708,12,4.60,0.25,
phenanthrene,target,178,176;152,,2715,12,4.55,0.25,phenanthrene-d10
```

These positions illustrate the format; match them to your acquisition and fold.
Promoting an identified peak into an existing CSV preserves valid column order and
extra columns, rejects duplicate names and replaces the file only after a complete
successful write. Reload that CSV to use the new target.

## Calibration levels and units

```csv
level,phenanthrene
L1,0.1
L2,0.25
L3,1.0
```

One row per unique level, one column per target compound. Amounts are finite,
non-negative amounts on tube in one consistent unit. Blank cells are unknown,
not zero. The default calibration needs three distinct usable levels per target
and calibration set; exclusions, overrides and zero-level policy affect coverage.

Guided setup records the output unit and actual constant IS additions. Switching
between pg, ng, ug and mg converts draft amounts; custom units require checking
the numbers. External level CSVs do not declare a unit and retain the label
“calibration units” until setup declares it. Concentration is amount unit/m³
when sample air volume is present. The table-driven quantitative route assumes
constant IS addition across the sequence; the separate calibration builder has
an amount-ratio route for varying additions. See the [user guide](user-guide.md).

## Sequence and method

```csv
file,type,level,cal_set,volume_m3
CAL_L1.cdf,calibration,L1,set1,
CAL_L2.cdf,calibration,L2,set1,
CAL_L3.cdf,calibration,L3,set1,
SAMPLE_A.cdf,sample,,,1.2
```

A type containing `cal` and a level marks a calibration run; a type containing
`blank` makes it available as a blank. `cal_set` groups curves and supports
bracketing in run order. `volume_m3` is sampled air volume. Saved edits in the
Runs table take priority over `sequence.csv` on project reload.

`method.json` can supply a confirmed period, for example
`{"modulation_period_s": 8.0}`. Use the actual method; 8 s and 10 s target
positions are not interchangeable. Native `.lsc` files do not declare the
modulation period, so confirm it under Load before quantitative analysis.

## Local starter tables supplied on 19 October

The files below are in the development workspace's ignored
`real_data/target_table_exports/`, with their own evidence and README. They are
private-data aids and are not shipped with the open-source package.

| File | Scope |
|---|---|
| `calibration_10s_targets.csv` | three inventory-matched IS rows for `Cal files for James`; start here for the August calibration acquisitions |
| `calibration_10s_targets_with_d12_review.csv` | those rows plus `D12 PAH m/z 240 (isomer unconfirmed)` for review |
| `exhaust_8s_targets.csv` | the three IS names with the earlier September exhaust/diesel positions |

The three names are acenaphthene-d10, phenanthrene-d10 and n-tetracosane-d50.
These tables help locate and review standards; they are neither a full target list
for the supplied inventory nor a calibration-amount file. The old chrysene-d12
assignment cannot distinguish it from benzo[a]anthracene-d12 on this evidence.
Perylene-d12 is omitted from the starter tables. The workbook column called RI
contains catalogue codes. Neither the tentative 1 ng IS recollection nor the
ambiguous `025ng` filename has been treated as a confirmed amount.
