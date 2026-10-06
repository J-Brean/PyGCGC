# 📈 Calibration builder on the August 2026 standards (6 October 2026)

Four standard runs, supplied without a compound list, taken through the 📈 Calibration builder
end to end: open the folder, clean and find blobs, suggest compounds, normalise to the internal
standards, fit the curves.

## Data

| Level | Run | Amount (from the name) |
|---|---|---|
| L2 | `2026-08-25_Cal_L2_5uL_025ng_on_tube_1_70eV.lsc` | 0.25 ng |
| L4 | `2026-08-24_L4_1ng_on_tube_1_70eV.lsc` | 1 ng |
| L6 | `2026-08-25_Cal_L6_5ng_on_tube_1_70eV.lsc` | 5 ng |
| L8 | `2026-08-24_L8_25ng_on_tube_1_70eV.lsc` | 25 ng |

Untreated 70 eV channel of each tandem acquisition. "025ng" was read as 0.25 ng: to be confirmed.
Three runs also have a `<run> 70eV.lsc` copy written by TOF-DS "Open TOF Data File" (SD01
spectra, a format PyGCGC does not read); folder opening listed them and failed to load them, and
now skips them. Period estimated from the runs: 10.0015 s (high confidence), not the 8 s of the
exhaust runs.

## Findings

1. **Phase.** The bleed-method phases came out 1.98, 1.42, 2.10 and 2.12 s (L4, L8, L2, L6).
   Blob pairs of near-identical spectrum put L8 0.73 s off the others; with each run's own phase,
   6 of the 800 largest L8 blobs could be matched in every level, and 123 with L8 moved. The
   builder now lines the phases up itself (L8 moved to 2.078 s, 122 pairs).
2. **Streaks.** The solvent front and an overloaded m/z 91 compound smear up the whole second
   dimension from 6.5 to 14 min and were the first "calibrants" suggested. Blobs in the solvent
   front and pieces of a smear (three or more blobs of the same spectrum spread over more than
   1 s of t2 in one modulation, none holding half the signal) are no longer matched.
3. **Internal standards.** Acenaphthene-d10, phenanthrene-d10, n-tetracosane-d50 and
   chrysene-d12 were found in every level and recognised by spectrum. Perylene-d12 was present
   but too weak and its qualifier ratios too far off to recognise. In the 25 ng run the IS
   signals were about half those of the other levels (acenaphthene-d10 1.6e7 against 2.9 to
   3.7e7; phenanthrene-d10 1.2e7 against 2.1 to 2.5e7).
4. **Calibrants.** 59 suggested: n-alkanes (m/z 57), alkylbenzenes (119, 133), methyl-, dimethyl-
   and trimethylnaphthalenes (142, 156, 170), biphenyl or acenaphthene (154), methylphenanthrenes
   (192) and others. Over the 100-fold range of amount their volumes rose only 5 to 15-fold
   (log-log slopes 0.4 to 0.9). Two causes are visible: a background of the same compounds at
   the low end (an intercept), and a plateau at the top, with the largest raw point of the
   calibrant blobs at a median 2.1e5 counts at 5 ng and 2.5e5 at 25 ng (none above 3.0e5).
5. **Curves.** IS-normalised, automatic weights and model: 3 of 59 pass the back-calculation
   check on all four levels; 24 of 59 with the 25 ng point left out. On 43 calibrants found in
   every level after 14 min, dividing by the internal standard halved the median back-calculation
   error of a linear fit at 1 ng (31 to 24%) and at 25 ng (54 to 28%).

## What is needed

The compound list and amounts of the mix, and the internal-standard amount per tube. With four
levels, Mandel's test has one degree of freedom, so a curve is almost never judged quadratic; a
fifth level between 5 and 25 ng, or a dilution of L8, would show where the response bends.

Reproduce with the builder (📂 Open folder, ▶ Clean and find blobs, 🔎 Suggest compounds,
➗ Normalise to IS); the session and exported tables are not committed (real data).
