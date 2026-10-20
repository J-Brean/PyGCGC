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

Linking the semi-standard non-polar records to the spectra, before the n-alkane rule: by CAS
number 17,950 entries, by NIST number a further 25,065, by normalised name 3,064; 197,814 of the 243,893
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
| 32 peaks at a time, row-wise (adopted) | about 20 ms (3.6 ms for the spectral products alone) |

A run of 2,000 fitted components therefore takes about 40 s once the library is loaded, plus the one-off 10 s above.
The weighted library and its decoys occupy about 200 MB each and the decoy build needs about a gigabyte transiently.

## Self-match check (a best case, not an accuracy estimate)

2,000 library spectra chosen at random, each multiplied by log-normal noise (σ = 0.15) and given a retention index
equal to its own (measured where available, otherwise estimated) plus Gaussian noise (σ = 8 units), were
matched against the whole library. The right *compound* was top for **99.3 %**. Of the 2,000, 274 were graded
confident, 1,725 probable, 0 tentative and 1 unknown. Most are only "probable" because their library index is an
estimate (which cannot make a peak confident), or an isomer with a nearly identical spectrum and index is within 30
points. This tests the arithmetic and the gates; it says nothing about how well real peaks, with co-elution and
background, will be named.

Named compounds, matched with their own library spectrum and index plus 7 units: decane, naphthalene,
phenanthrene, eicosane and pyrene all returned as top candidate (match factor 999, z between 0.2 and 0.4).
Phenanthrene and anthracene (997), pyrene and fluoranthene (997), and naphthalene and a methylene-indene (981) were flagged
`ISOMER_AMBIGUOUS`, as they should be: their spectra and indices are too close for either to separate them.

## Real runs: a blank-like run and a sample

*(Filled in from the run of 20 October; see the section below.)*

## What is not checked

- No mixture of known composition. The August standards' starter names and the supplied inventory are the next test,
  but their amounts and some identities are unconfirmed (ROADMAP section 1).
- The 15-unit allowance for the difference between a GC×GC first dimension and the one-dimensional columns of the literature
  is a starting value; the median measured-minus-library index of the confirmed standards would set it.
- The 850, 750 and 650 match-factor thresholds, the 30-point runner-up margin and the 5 % false-discovery limit.
- Which of 5 % phenyl's literature columns resemble the instrument's first dimension (DB-5, Rxi-5, HP-5 and so on are pooled).
- The estimated-index correction is a single offset; its error probably depends on compound class.
