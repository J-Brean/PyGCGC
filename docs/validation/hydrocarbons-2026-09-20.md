> Describes the earlier acquisition (10 s modulation, 2025-10-09) that is no longer in
> real_data/. The current runs are covered by [audit-2026-09-24.md](audit-2026-09-24.md).

# Real-data checks: siloxane reference, alkane ladder, markers, diagnostic pairs, phenanthrene/anthracene

20 September 2026. Untreated 70 eV native file, modulation period 10 s, phase 1.81 s, 215,015 scans,
m/z 38 to 300, 72 minutes. Every number below is printed by
`benchmarks/validate_real_data.py`. No data were committed or uploaded, and no threshold outside
`SiloxaneSettings`, `HydrocarbonSettings` and the marker offsets was changed.

## Siloxane reference spectrum

The run has column bleed that rises from about 45 min and dominates m/z 207. The bleed spectrum is
not the textbook cyclic Dn pattern the old default assumed (73 = 100, 147 = 45): m/z 73 is only about 7% of the base ion.

Method: 20 s median spectra from 40 to 71 min, each ion regressed on m/z 207, slope × 100 taken as the
reference intensity. Ions with R² of at least 0.9 and slope of at least 2 were kept, which removes the constant
background (m/z 39, 40, 44).

| m/z | 207 | 281 | 208 | 209 | 96 | 191 | 133 | 73 | 193 | 282 | 177 | 283 | 97 | 192 | 75 | 267 | 147 | 119 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Relative | 100 | 18.1 | 17.2 | 10.2 | 9.6 | 8.1 | 7.3 | 6.7 | 6.0 | 4.3 | 3.2 | 2.8 | 2.8 | 2.6 | 2.5 | 2.3 | 2.1 | 2.0 |

Discrete cyclic peaks also occur (29.3 min and 22.3 min, t2 about 0.5 s): baseline-subtracted, the 29.3 min peak is
73 = 100, 281 = 56, 147 = 50, 282 = 16, 283 = 10, 221 = 7, 207 = 7. It is shipped as a preset.

Removal on the raw run with the bleed reference takes 97.0% of m/z 207 out of 55 to 70 min and leaves m/z 57 unchanged
to five figures. The bleed reference leaves the 29.3 min cyclic peak untouched (0% of its m/z 73 removed); the cyclic
reference removes 85% of it and 0.1% of the siloxane-ion signal elsewhere. A bleed reference flags 873 of 4068 blobs at
the current 0.3 threshold, nearly all in the last 20 minutes; the threshold itself is not validated.

## n-Alkane ladder

The low-t2 ridge in the m/z 57 image is the n-alkane ladder: the molecular ions of C12 to C15, C17 and C18 (m/z 170, 184,
198, 212, 240, 254) peak on it at the ladder positions, and the spacing shrinks smoothly from 4.33 min (C13 to C14) to
2.33 min (C26 to C28). t2 rises from 0.67 s (C12) to 1.38 s (C28).

The old search failed on this run. It took the tallest m/z 57 pulse within a few seconds of the clicked time at any t2;
at the C18 seed that pulse is at t2 5.9 s. The clicked t1 is also the middle of the modulation, up to half a period (5 s)
from the apex time. The search is now two-dimensional, predicts each carbon from its two nearest neighbours, and takes the
tallest pulse whose qualifier ratios match. From seed pairs (C13, C18), (C12, C14) and (C20, C24), and from seeds clicked a
column off, it returns the same ladder, C12 to C28, 17 carbons, every step monotone.

Qualifier ratios along the ladder (baseline-subtracted, summed over five rows at the apex):

| | C12 | C13 | C14 | C15 | C16 | C17 | C18 | C19 to C28 |
|---|---|---|---|---|---|---|---|---|
| 71/57 | 0.54 | 0.63 | 0.60 | 0.64 | 0.64 | 0.81 | 0.63 | 0.66 to 0.69 |
| 85/57 | 0.32 | 0.41 | 0.37 | 0.39 | 0.39 | 0.37 | 0.41 | 0.41 to 0.51 |

C17 co-elutes with pristane. The medians over C13 to C28 without C17 are 0.666 and 0.438, the shipped defaults
(0.66, 0.44, 30% tolerance). Lighter alkanes fall away (C11 0.48/0.25, C10 0.36/0.17), so the gate ends the ladder at C12.

## Pristane and phytane

Comparing the t1 centroid of the m/z 183 excess (183 above 1.5% of 57, a branched-alkane ion) with the alkane's own
m/z 57 centroid: pristane +1.3 s (+0.006 C), phytane +17 s (+0.084 C). One column is 0.047 to 0.049 C, so pristane
co-elutes with n-C17 within a column (its 183 excess falls in one column, 3336 counts) and the phytane offset rests on a
weak signal (an excess of 53 to 397 counts per column across five columns). The old guess of -0.3 C would have searched 65 s before n-C17. New defaults:
+0.01 C and +0.08 C. Literature RIs on a 5% phenyl phase (1708 and 1812) would give +0.08 and +0.12.

## Homologous-series diagnostic ion pairs

The old finder produced 18 to 22 "members" for (91, 105), (57, 71) and (43, 57), all of them blobs at one t1 (the solvent
front at 7.4 and 10.4 min, t2 4 to 10 s). It assumed a constant spacing in seconds (the real ladder spacing halves) and kept several
members per position. It now works in carbon-equivalents from the alkane ladder, keeps one member per carbon and
requires a consistent t2 offset from the alkane line.

| Pair | Candidate blobs (major share of at least 10%) | Series found |
|---|---|---|
| (57, 71) | 86 | 14 members, carbon-equivalents 9.9 to 25.0, all within 0.1 C of an integer: the n-alkane ladder, a positive control |
| (91, 105) | 80 (most in the solvent front) | none |
| (82, 83) | 6 | none |
| (142, 156) | 2 | none |
| (58, 59) | 1 | none |

None of the four spec pairs is validated. The result is a property of this sample and of blob detection on it (4068 blobs,
with the solvent-front streaks and background fragmentation described in `docs/history/blob-diagnosis-and-plan.md`), not proof that the
series are absent. 142/156 would not be a one-carbon series in any case: methylnaphthalene to dimethylnaphthalene adds
about 1.3 C.

## Phenanthrene and anthracene

Two C14H10-type peaks (m/z 178 base, 176 at 17%, 152 at 8%) lie in the RI 1700 to 1950 window. The
peak at t1 42.45 min, t2 2.98 s has RI 1785.7. The larger one at 44.45 min, t2 3.26 s (RI 1845) also carries m/z 188 at 71% of
178, a co-eluting d10 standard, and is 59 RI units away, too far to be anthracene against phenanthrene. No resolved anthracene is
present, so the pair could not be measured directly.

What could be measured is the ladder's RI precision and the behaviour of the RI penalty. One 10 s column is 5.0 RI units at
this position. With library RIs assumed (not measured) at 1777 for phenanthrene and 1793 for anthracene, and the real spectrum used for
both, the default settings give:

| Query RI | 1775 | 1780 | 1785 | 1790 | 1795 |
|---|---|---|---|---|---|
| Default: winner, margin | phenanthrene 0.048 | phenanthrene 0.030 | phenanthrene 0.000 | anthracene 0.030 | anthracene 0.048 |
| `ri_penalty_full` 30: winner, margin | phenanthrene 0.160 | phenanthrene 0.100 | phenanthrene 0.000 | anthracene 0.100 | anthracene 0.160 |

The real peak (1786) sits on the midpoint. A 16-unit RI difference gives at most a 0.048 margin at the defaults against an
ambiguity threshold of 0.02, and a 2a call flips between RI 1780 and 1790, two columns apart. Telling the pair apart by RI
alone needs RI good to about 4 units, so a co-injected standard to calibrate the library offset and sub-column t1 centroiding.
Identification thresholds were left unchanged, as agreed.
