# Library matcher and peak-level blank subtraction, 20 October 2026

Evidence for the 📚 matcher (NIST 2011 SDF files) and the ➖ blank subtraction added on 20 October
2026. Everything here was run on the branch `feedback-2026-10-20` with the NIST files in
`C:\REAL_DATA_FOR_PYGCGC\nist\SDF` (licensed to James; no library content is stored in the
repository or in this note beyond counts and compound names used as examples). All
thresholds the matcher uses are **provisional**: James agreed them on 20 October 2026 as
starting values, and nothing below validates them against a mixture of known composition.

## The files

| File | Records | Used for |
|---|---|---|
| `mainlib.sdf` | 212,961 EI spectra (108,401 with a CAS number; 197,301 with an estimated Kovats RI) | spectra |
| `replib.sdf` | 30,932 replicate spectra | extra spectra of the same compounds, collapsed to the best per compound at search time |
| `nist_ri.sdf` | 346,757 experimental RI records. Column class: 113,379 standard non-polar, 149,090 semi-standard non-polar, 84,243 standard polar. Index scale: 155,748 "Normal alkane RI", 111,778 "Linear RI", 73,379 "Kovats RI", 5,852 "Lee RI". | literature retention indices |

The matcher uses the semi-standard non-polar records (a 5 % phenyl column, confirmed by James) of
"Linear RI" and "Normal alkane RI" type from temperature-programmed ("Ramp" or "Complex")
measurements. Isothermal records and the PAH-based "Lee RI" scale are excluded by default.
n-Alkanes appear in `nist_ri.sdf` once, with the class "All column types" and the type "Normal alkane RI
value specified by scale definition" (for example 1000 for decane); those entries take that value with a spread of 1.

## Reading and linking (timings on this workstation)

| Step | Time |
|---|---|
| First read of `mainlib.sdf` and `replib.sdf` into 243,893 spectra (25.6 million non-zero peaks, m/z up to 1760), written to the local cache | 28 s |
| Reading the cached spectra in a later session | a few seconds |
| Reading `nist_ri.sdf` (346,756 numeric records) | 18 s |
| Linking indices to spectra, semi-standard non-polar | 2 s |
| Building the weighted copy of the library for one m/z range, plus its shuffled-decoy copy (once per session and range) | 10 s |

Linking the semi-standard non-polar records to the spectra (final code, including the n-alkane rule): by CAS
number 18,075 entries, by NIST number a further 25,065, by normalised name 3,064; 197,689 of the 243,893
spectra (81 %) have no measured index for that class. For those the matcher falls back on NIST's
estimate (below) or, if that is also missing, on the spectrum alone and flags `RI_LIBRARY_MISSING`.

## NIST's estimated retention index

`mainlib.sdf` carries an *estimated* Kovats index and an estimation error for 197,301 entries. Compared
with the measured indices of `nist_ri.sdf` for the compounds that have both:

| Measured column class | Compounds with both | Estimate minus measurement, median | Robust spread (1.4826 MAD) | Median stated error |
|---|---|---|---|---|
| Standard non-polar | 18,850 | 0.0 | 62.3 | 72 |
| Semi-standard non-polar | 45,219 | −21.0 | 89.0 | 89 |
| Standard polar | 9,666 | −454 | 262 | 55 |

So the estimate describes a **standard non-polar** column; on a 5 % phenyl column it reads about 21
units low; on a polar column it is unusable. Splitting by NIST's own stated error: with an error under 120 the
spread of (estimate − measurement) is close to the stated error (38 to 106 units for stated
errors of 40 to 120, semi-standard) but estimates with a stated error of **120 or more** were off by −319 (median) with a
spread of 222, so they are not used. The matcher therefore uses the estimate only for the two non-polar
classes, corrected by the measured bias (21 units for 5 % phenyl), with a spread of 1.25 times NIST's stated
error on 5 % phenyl and 1.0 times on standard non-polar, never for a polar column, and
flags it `RI_ESTIMATED`. The inference that the estimate describes the standard non-polar class comes from this comparison;
it is not taken from NIST documentation.

## Search speed

The library is held as a sparse matrix and queries are multiplied against it 32 at a time.

| Method | Time per peak (243,893 spectra, with decoys) |
|---|---|
| One peak at a time, column-wise gather | about 100 ms |
| 32 peaks at a time, row-wise (adopted) | about 25 ms (3.6 ms for the spectral products alone; 2,000 peaks in 50 s) |

A run of 2,000 fitted components therefore takes under a minute once the library is loaded, plus the one-off 10 s above.
The weighted library and its decoys occupy about 200 MB each and the decoy build needs about a gigabyte transiently.

## Self-match check (a best case, not an accuracy estimate)

2,000 library spectra chosen at random, each multiplied by log-normal noise (σ = 0.15) and given a retention index
equal to its own (measured where available, otherwise estimated) plus Gaussian noise (σ = 8 units), were
matched against the whole library with the final code. The right *compound* was top for **99.3 %**. Of the 2,000, 274 were graded
confident, 1,725 probable, 0 tentative and 1 unknown. Most are only "probable" because their library index is an
estimate (which cannot make a peak confident), or an isomer with a nearly identical spectrum and index is within 30
points. This tests the arithmetic and the gates; it says nothing about how well real peaks, with co-elution and
background, will be named.

Named compounds, matched with their own library spectrum and index plus 7 units: decane, naphthalene,
phenanthrene, eicosane and pyrene all returned as top candidate (match factor 999, z between 0.2 and 0.4). Decane and
eicosane are graded confident (their index is exact by definition); the others are probable because of their isomers.
Phenanthrene and anthracene (997), pyrene and fluoranthene (997), and naphthalene and a methylene-indene (981) were flagged
`ISOMER_AMBIGUOUS`, as they should be: their spectra and indices are too close for either to separate them.

## Real runs: a blank-like run and a sample

Files in `C:\REAL_DATA_FOR_PYGCGC\background_removal`: `9003_Test1_RDE_Dilutor` (treated as the blank; James did not say which
of the two was which) and `9002_Test3_RDE(w)_D_bDOC` (the sample). Both native 70 eV runs, 618,000 scans, 103 min, 10 s
modulation, m/z 38 to 400. The blank is far from empty: the sample's TIC is 1.05 to 1.5 times the blank's in
every 10-minute block (for example 1.0 × 10⁹ against 1.5 × 10⁹ at 40 to 50 min). Read only. Because fitting two
103-minute runs takes a long time, the check used a window of 30 to 40 min (60,000 scans each), one fitting process, the default tilted
Gaussian MS-aware fit: 226 components in the blank, 234 in the sample (about 7.5 min each).

**Phases.** The automatic estimates were 0.86 s for the sample (void method) and 1.37 s for the blank (bleed
method), 0.51 s apart on the same method and instrument.

| | Peaks paired with the blank (of 234) | Median spectrum match of pairs | Median offset of pairs | Signal removed, factor 1 |
|---|---|---|---|---|
| Each run at its own estimated phase | 5 | 0.76 | 0.64 s in t₁, 0.047 s in t₂ | 0.1 % |
| Both at 0.86 s | 18 | 0.91 | 2.5 s in t₁, 0.13 s in t₂ | 4.2 % (2.3 % at factor 0.5) |

A phase that differs by 0.5 s moves every peak by 0.5 s in t₂, outside the 0.2 s matching window, so
only 5 pairs were found, by chance. With a common phase 18 pairs were found; four were removed completely and 14 reduced.
223 sample peaks were above blank, 11 near blank. The sample-to-blank volume ratio of the pairs had median 2.2 (middle half 1.6 to 8.6).
The two largest sample peaks (at 36.2 min and 37.5 min, t₂ 5.5 s and 7.2 s) have partners with spectral match 0.95 and lose 11 and 13 %.
The residual 0.13 s median t₂ offset and 2.5 s in t₁ show real retention drift between the runs, which alignment
would absorb. **This is why the blank step now warns when the phases differ and the runs are not aligned**, and
why the anchoring window insists the phases are checked. Only about 8 % of the sample's peaks have a partner in this blank: it
removes a modest share of the signal at this window, and says so.

The largest of the two paired peaks, at 37.5 min, matches anthracene-d₁₀ (probable, match factor 880): very likely a deuterated
internal standard spiked into both runs. Subtracting a blank from a spiked standard would be wrong, so internal standards
are left alone in the quantified amounts and in the all-signals table. The peak-level table itself does not know which peaks are
standards; read those rows accordingly.

**Matcher on the same window** (sample, 234 peaks, no alkane ladder, the full NIST library, 15 s including the one-off setup):
0 confident, 47 probable, 59 tentative, 128 unknown. Nothing can be confident without a retention index, by design; every
peak carries `RI_MISSING`, and most `RI_ESTIMATED` and `COELUTION`. The ten largest peaks mostly return low match factors
(400 to 530) with long ester or acid names, the signature of a mixed spectrum from a co-eluting cluster at t₂ about 5.5 s; the
clear ones are anthracene-d₁₀ (880), dodecyl acrylate (913) and a tetrahydro-phenylnaphthalene (927, isomer ambiguous). A real
ladder and the 🧬 class would be needed to say more. This is a demonstration that the pipeline runs end to end on real data, not
a validation of the names.

## What is not checked

- Two runs only for the blank step, one window, and the sample/blank assignment is assumed. Alignment of these runs was not tried.
- Whole-run fits and the 14 eV channels were not used.

- No mixture of known composition. The August standards' starter names and the supplied inventory are the next test,
  but their amounts and some identities are unconfirmed (ROADMAP section 1).
- The 15-unit allowance for the difference between a GC×GC first dimension and the one-dimensional columns of the literature
  is a starting value; the median measured-minus-library index of the confirmed standards would set it.
- The 850, 750 and 650 match-factor thresholds, the 30-point runner-up margin and the 5 % false-discovery limit.
- Which of 5 % phenyl's literature columns resemble the instrument's first dimension (DB-5, Rxi-5, HP-5 and so on are pooled).
- The estimated-index correction is a single offset; its error probably depends on compound class.
