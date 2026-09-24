# Real-data blob diagnosis and development plan

20 September 2026. This is a measured baseline and proposed development sequence.
The new route is now a first-class, saved choice in the main Peak detection step;
legacy blob borders remain alongside it for compatibility.

## Plain-English summary

### What was going wrong

The old detector first combined evidence from hundreds of mass channels into one
image. Small positive leftovers from background removal accumulated, so broad
areas with no obvious chromatographic peak looked extremely strong to the
software. Once a false peak was started, its border was allowed to grow down to
only 1% of its height. Because the supposed "background" was already above that
level, borders spread across large parts of the chromatogram and joined several
real peaks together.

There was also a smaller consistency bug: the GUI and batch path estimated the
instrument's minimum stored signal from already-cleaned data. That made the
minimum about 41 times too small in this file. Fixing it was necessary, but did
not solve the main accumulation problem. Finally, a hard border assumes every
point belongs to exactly one compound. Real overlapping peaks do not obey that
assumption: two compounds can contribute signal at the same time and mass.

### What has been tried

The experimental detector now looks for a peak separately on each mass channel,
then asks whether several channels agree on its position and whether the total
signal also has a local peak there. It records a candidate centre and a small
analysis rectangle instead of immediately claiming every surrounding point for
one compound. On the ten-minute test segment this gives 231 candidates instead
of 612 old blobs. The reviewed background band goes from roughly 99% covered by
blob borders to no candidate centres, while the difficult area still contains
21 candidates.

Inside a difficult rectangle, the code then tries to explain the measured data
as overlapping contributions from several components. Each component has a
chromatographic shape and a mass spectrum. At the supplied location near 15 min
and 1.2 s, several independent checks support two repeatable signal patterns. We
cannot yet say that they are two identified compounds, because this file is not
a known mixture and no reference identities were supplied.

Recovery tests provide a known answer. In these tests, "injection" means adding
a software copy of a peak to a real quiet part of the chromatogram at an amount
chosen in advance. Nothing is injected into the instrument. The program is then
asked to separate the copies and report their amounts. With ideal bell-shaped
peaks, the current fast shared-width fit puts both amounts within 10% in 5 of 9
combinations of spacing and relative size (the older unrestricted fit passed 4).
Equal-sized peaks work well, but a peak only 10% as large as its neighbour does not.

A newer control uses the confirmed n-C16 peak from the full acquisition instead
of a perfect bell curve. The existing bell, tailing, tilted and narrow-cell
models all lose roughly 19% of the larger component on a representative clean
case. A new experimental fitter instead moves copies of the measured n-C16
shape and solves their mass spectra and amounts. It passes all nine cases without
real background. Allowing the fitted background to be signed, as the cleaned
data actually are, raises recovery in the chosen real-background window to 6 of
9: every equal-size and 30%-minor case passes, while every 10%-minor case fails.
This establishes a useful initial operating limit rather than universal success.

### What happens next

1. Extract a shape library from several confirmed n-alkanes and, when available,
   pure standards. Test how the shape changes across retention time.
2. Test the measured-template model against those held-out shapes rather than
   using n-C16 as both the template source and recovery shape.
3. Add an explicit quality rule: weak components and unsuitable background must
   be reported as an unresolved group rather than assigned unreliable amounts.
4. Repeat the known-answer recovery grid over peak spacing, relative size and
   signal strength. This will define when PyGCGC may report two amounts and when
   it must report an unresolved group.
5. Annotate held-out real windows containing clear peaks, weak peaks, background
   and edge cases. The route is already available in the GUI for inspection;
   only after those checks pass should its fitted amounts feed quantitative results.

### What replaces blob borders in the interface

Step 4 is now **Peak detection and separation**, with a saved method selector:
**MS-aware separation** or **Legacy blob borders**. MS-aware detection runs over
every loaded run from the same Go action as legacy detection. Its normal table
has one row per fitted component, with position, estimated intensity, diagnostic
ions and confidence. The 2D view draws overlapping 10% and 50% contours for each
component. Candidate rectangles remain calculation windows and are only shown
when no fit is available or a group is opened for review. **Resolve selected
region…** opens that group's calculation window.

The old yellow borders and add/merge/delete tools appear when Legacy blob borders
is selected. Existing siloxane blocking, hydrocarbon feature totals and library
search still depend on `BlobSet`, so those paths compute or use legacy blobs for
now. They must be migrated to persistent components separately; the GUI does not
mislabel candidate rectangles as quantitative blobs merely to satisfy them.

After fitting, the new Component contributions tab draws each component in a
different translucent colour over the observed signal. Each has a 10% and 50%
contribution contour. Contours may cross and colours may overlap because both
components can contribute to the same measurement. The table reports each
component's estimated amount, size relative to the largest component, and a
confidence label. Following James's decision, a weak component is retained and
reported; below the current 30% recovery boundary it is orange and labelled
low confidence, and the whole set is described as an unresolved group.

When an alkane ladder has been found for the run, the fit model menu also offers
"measured alkane template". It builds a reference profile from n-C16 (or the
nearest usable confirmed alkane) in that run. Component groups now persist in the
project, can be replaced by an accepted manual fit with undo/redo, feed component
spectra to identification, and export to CSV with their warnings. They do not yet
feed siloxane handling, batch matching or production quantification.

## Reproduction and scope

Read HANDOVER.md and README.md, traced the detector, background, overlap and
GUI/batch/quantification callers, and ran six controlled detector variants on
the supplied untreated 70 eV acquisition. The folder contains one acquisition
with untreated 70 eV, 14 eV and vendor DBC channels. This audit uses untreated
70 eV only.

- Period: **10 s**; phase: **1.81 s**, relative to the first scan as in the app.
- Original first scan: 0.2265625 s; first modulation boundary: 2.0365625 s.
- Cached context: 600.2265625–1200.2265625 s, 30,001 scans, 60 modulations,
  500 rows per modulation (20 ms row pitch), all original unit-mass channels.
- Main comparison: **t1 = 14–16 min, t2 = 0.6–2.0 s**.
- Background comparison: t1 = 14–16 min, t2 = 5–9 s. This is a visually
  background-dominated comparison band, **not a known blank or proven zero-signal region**.
- Crop folding uses t0 = original t0 + 600 s, an integer number of periods.
  Resetting phase relative to an arbitrary crop's first scan would move the data.

Run `python examples/diagnose_real_blobs.py` in the project environment.
The first run imports the source; subsequent runs reuse the sparse crop.
All derived real-data files are under `real_data/blob_diagnosis/`:

- `segment.npz`: approximately 8 MB sparse scan data for preprocessing tests.
- `segment.json`: original origin, crop origin, source path/SHA-256 and fold settings.
- `focus_cube.npz`: approximately 1.1 MB, raw and cleaned all-ion arrays for the
  focus window, with t1 midpoint and t2 axes. Use the wider scan crop to expand
  windows when a fitted peak reaches an edge.
- `metrics.json`: machine-readable results and per-blob extents.
- [Border comparisons](../real_data/blob_diagnosis/comparison.png).
- [Individual ion maps](../real_data/blob_diagnosis/ions.png).

The context is deliberately wider than the analysis window because the background
pattern filter spans 51 modulations. This is still a crop experiment: its baseline
fit and noise blocks can differ from a whole-run fit. Full-run preprocessing
followed by extraction, plus wider-context comparisons, is an early validation
task. No real data have been committed or uploaded.

## What the experiments show

Counts below are for the entire ten-minute context. Coverage and median scores
are for the background comparison band. Scores are the detector's current
numbers; they must not be interpreted as calibrated statistical S/N.

| Variant | Blobs | Median background score | Background pixels labelled |
|---|---:|---:|---:|
| Raw, default detector | 550 | 489 | 99.6% |
| Cleaned, default detector | 612 | 660 | 99.0% |
| Cleaned, original raw noise floor | 475 | 408 | 99.6% |
| Cleaned, TIC-only detector | 50 | 1.53 | 0% |
| Cleaned, border cutoff raised from 1% to 10% | 612 | 660 | 98.5% |
| Cleaned, less smoothing | 1,114 | 3,004 | 93.8% |

Each detector-plus-diagnostic case takes roughly 2–4.3 seconds on this machine,
excluding import, background removal, cube export and plotting. Reduced smoothing
uses 0.25 modulations and 0.015 s Gaussian SD, versus 0.6 and 0.03 by default.

In the focus window the cleaned default has 22 intersecting blobs. Their maximum
individual spans are 110 s along t1 and 0.8 s along t2 (these maxima need not
belong to the same blob). Raising the relative cutoff reduces the maximum spans
to 80 s and 0.46 s, but leaves the background-detection problem substantially intact.
TIC-only gives 16 intersecting blobs; that is not evidence it has recovered the
correct number of compounds. Reduced smoothing gives 33 and increases false
structure elsewhere.

The ion maps differ substantially: for example, m/z 71 emphasises the lower
ridge, while 128 and 142 show different structure near the crowded 1.2 s region.
This supports investigating spectral separation. These are mixed ion images,
not identified compounds or a verified component count.

## Diagnosis

**1. The detection score accumulates background across ions.** In
`core/peaks.py::_ion_snr_cube`, each smoothed ion image is divided by its spread,
then all 261 eligible ions are summed. There is no local baseline centring in
this step and no calibration for the number or correlation of ions. Sparse ions
fall back to the storage-derived floor. Residual positive background therefore
accumulates into very large scores. `find_blobs` takes the maximum of that sum
and TIC S/N. Both seed threshold 10 and growth threshold 2 are passed throughout
the background comparison band in all ion-summed variants. Watershed consequently
partitions background fluctuations into apparent peaks.

Simply dividing by sqrt(number of ions) would not fix nonzero means, correlated
fragments or censored measurements. The combined statistic needs a measured null
distribution and evidence of chromatographic peak structure.

**2. The raw noise-floor rule is inconsistently applied.** `find_blobs` defaults
to `storage_threshold(run)/3`. The GUI's `BlobState.run_all` and batch pipeline
pass the displayed/processed run without an explicit floor. Quantification
already passes `storage_threshold(raw)/3`. Here the raw floor is 0.6667 and the
cleaned floor 0.01634, a factor of 40.8. Tiny positive cleaning residuals are not
the instrument's original storage threshold. Fixing this inconsistency is
necessary, but the raw-floor experiment demonstrates that it is insufficient.
The storage floor itself is only a provisional bound, not a complete noise model.

**3. Border growth has no reliable background stop.** The same inflated score
controls both seeding and watershed support. A weak seed has a very low 1%-of-apex
cutoff, so background can remain attached over long distances. Hard watershed
also cannot allocate two compounds' contributions to the same scan/ion.

There is a smaller implementation defect: `_segment_chunk` trims using the first
watershed's per-label thresholds and then runs a global watershed again. Pixels
can change owners, so the final ownership need not obey its own peak's cutoff.
The audit finds 9 such pixels for cleaned defaults and 29 with the 10% cutoff
over the context. This deserves a regression test, but is not the main explanation
for near-total background coverage.

**4. Resolution and evidence are being discarded early.** Default t1 smoothing
has SD 6 s at this period. The fixed maximum filter is five rows by three
modulations, and the detector collapses all spectra to one surface before
seeding. Weak components on strong neighbours can lack a maximum on that surface.
The co-elution flag examines only ions contributing at least 5% of a blob's total
spectrum and just one apex per ion; shared dominant fragments can hide a smaller
component. These are limitations to test, not evidence that every broad blob
contains a known number of peaks.

**5. Existing deconvolution was a starting point, with important gaps.** The
overlap dialog is manual and its results do not replace detection or quantification.
`blob_cube` inherits the existing footprint and normally keeps only ions with
at least 3% of the summed spectrum (with a small fallback set). That can discard
weak diagnostic ions or truncate component tails. In the audited pre-change
code, the shape fit's inner spectral NNLS was unweighted while its outer residual
was weighted, and BIC omitted the estimated spectral coefficients from its
parameter penalty. Solving parameters analytically does not remove their
statistical degrees of freedom. These objective and penalty problems were
corrected during the first implementation pass below. Component
count and Jacobian uncertainty still must not be treated as validated evidence.
MCR uses a simple approximate unimodality operation and needs restart/stability
assessment and an explicit background treatment.

## Phased approach and completion criteria

### Phase 1: Establish the small benchmark and repair preprocessing contracts

The initial crop, ion maps and baseline comparisons above are complete. Next:

- Freeze several windows: the crowded focus, visually isolated peaks, weak peaks,
  background-dominated regions, and a wrap/edge case. Keep some out of tuning.
- Extract the same windows after whole-run preprocessing and compare with crop
  processing. Retain raw data, signed residuals and masks for missing scans.
- Pass raw-derived floor/noise provenance consistently through GUI, batch and
  quantification. Preserve per-ion noise information rather than recalculating
  an instrument floor from cleaned residuals.
- Inspect the dark patches introduced by cleaning and test despiking, baseline
  and repeating-pattern removal independently. Do not assume cleaning preserves
  all weak peaks. Compare vendor DBC only as another processing result, not truth.

**Gate:** reproducible timing and crop/full-run agreement within documented
tolerances, identical detector inputs across callers, and no unexplained loss of
the selected isolated peak traces. Add focused regression tests for those contracts
and final border-cutoff ownership. Keep a single fixed baseline report for comparison.

### Phase 2: Reliable peak candidates and bounded analysis regions

- Estimate local baseline and noise per ion, acknowledging intensity dependence,
  stored-value censoring and correlation introduced by folding/smoothing.
- Detect peaks within each modulation's ion traces using prominence, plausible
  width and support over successive scans. Link evidence across modulations using
  t2 position, profile similarity and spectra. Handle genuinely narrow t1 peaks
  explicitly rather than imposing a universal two-modulation requirement.
- Combine agreeing ion traces into candidates; allow a weak distinctive ion to
  seed a component hidden in TIC. Avoid demanding that every useful ion be unique.
- Separate candidate acceptance from the lower threshold defining an analysis
  region. Use local excess signal to stop growth. Watershed can provide initial
  regions and visual outlines for isolated peaks.
- Choose thresholds from empirical background/blank behaviour. Label uncertain
  candidates and regions that touch an analysis edge.

**Gate:** substantially fewer detections in reviewed background windows while
retaining reviewed weak and isolated peaks. Measure candidate precision/recall
where annotations exist, detections per background area, fragmentation/merging,
and stability to modest threshold/context changes. A lower total blob count alone
does not pass this gate. Pick numerical operating targets after the annotations
are available; a true false-discovery rate cannot be claimed from this one sample.

### Phase 3: Separate overlaps using both GC profiles and MS

Fit small padded rectangles around candidate groups, independent of the old hard
footprints. Start with the region near 15 min / 1.2 s and retain all informative
ions, selected by local evidence rather than global spectral percentage.

Use the model `signal(position, ion) = background + sum(profile_k * spectrum_k)
+ residual`. Component contributions may overlap at every position. Normalise
spectra consistently so integrated profiles have a defined signal-volume meaning.

Compare two implementations on exactly the same windows and weighting:

1. Corrected joint Gaussian/asymmetric shape fits with consistent weighted NNLS,
   constrained widths, multiple starts and explicit background.
2. Noise-weighted MCR-ALS with nonnegative components, soft chromatographic
   smoothness/shape constraints and multiple starts. Preserve signed measurement
   residuals; do not automatically force negative noise to become positive signal.

Account for movement of t2 across modulations if residual maps show it. A fixed
axis-aligned 2D Gaussian may be too restrictive. For difficult regions, compare
fitting original scan-time traces with fitting the interpolated cube.

**Gate:** residuals no longer contain repeatable chromatographic peaks; component
spectra and areas are stable across starts, neighbouring window sizes and reasonable
weighting choices. Test recovery by injecting measured isolated-peak profiles and
spectra into real backgrounds at known separations and abundance ratios. Keep the
untouched real overlap as a separate test, since injection truth is constructed.
An initial area-recovery target of 10% for adequately strong injected components
matches the earlier project ambition; weak components need signal-dependent limits.

### Phase 4: Decide which separations are supportable

Compare component counts using held-out observations, noise-scaled residuals and
stability, with complexity penalties consistent with the actual fitted model.
Account for correlated folded pixels when defining holdouts. Inspect withheld ion
traces and diagnostic-ion ratios as well as TIC reconstruction. Use bootstrap or
noise perturbations that preserve scan/modulation structure.

**Gate:** extra components must explain reproducible signal and retain interpretable
spectra/areas under perturbation. Report unresolved groups when competing solutions
are equally plausible. Exact co-elution with indistinguishable spectra cannot be
uniquely separated from this acquisition alone. Even distinct spectra can be
non-identifiable without sufficient variation, selective ions or reference information.
Do not turn an arbitrary MCR factorisation into confident compound identities.

Known mixtures/standards and independent repeat acquisitions provide stronger
validation when available. The existing 14 eV channel might add information later,
but needs time alignment and energy-specific spectra and response; it is not an
independent replicate or directly interchangeable with 70 eV.

### Phase 5: Integrate validated components into PyGCGC

Introduce component results containing profile, spectrum, integrated signal,
uncertainty/stability and unresolved-group membership. Keep region labels as
navigation/display geometry. Show overlapping contours or contribution maps,
observed-versus-fitted ion traces and residuals. Integrate each fitted component
once, with background and unexplained signal accounted for separately.

Wire the same validated component results into identification, quantification,
batch matching and exports. Add settings/provenance, cancellation and local cube
caching. Run whole-file checks, held-out windows and existing core/GUI regressions;
test invariance to tile seams and wrapping. Preserve the legacy detector as a
comparison mode until the new route passes these checks.

**Gate:** GUI/batch/quant outputs agree, fitted signal accounting is consistent,
weak/isolated peaks survive, and failure/ambiguity is visible. Only then make the
new route the default. PARAFAC2 or shift-aware multi-run methods are a later option
when several independent runs justify their extra assumptions and complexity.

## Methodological references

NIST's [AMDIS overview](https://chemdata.nist.gov/mass-spc/amdis/overview/page5.htm)
supports separating noise analysis, local baseline handling and spectral
deconvolution. It is useful precedent for the proposed architecture, not evidence
that PyGCGC currently implements or matches AMDIS performance.

[Schneide et al., SIML preprint](https://arxiv.org/abs/2412.12114) describes
shift-aware constraints within MCR-ALS for GC×GC-ToF-MS and comparisons against
MCR-ALS/PARAFAC2 variants. It motivates testing t2 movement and fit stability;
no superiority of that method has been established on this acquisition.

Recommended next implementation milestone: annotate held-out weak, isolated,
background and edge windows, then repeat the recovery benchmark with measured
isolated-peak profiles. The present real-data and injection audits establish the
failure mechanism and a first operating envelope; they do not establish chemical
identities or complete sensitivity.

## Implementation progress, 20 September 2026

The first Phase 1–3 pass is implemented, with the legacy detector still the GUI
default:

- GUI and batch blob detection now use the raw acquisition's storage floor, as
  quantification already did. The batch cache key was versioned so old blobs are
  not silently reused.
- Final watershed footprints are checked against their final owner's relative
  cutoff and reduced to the connected part containing that owner. This fixes the
  post-reflood ownership defect.
- `core/candidates.py` provides `find_candidates`, `candidate_regions` and
  `candidate_cube`. It preserves per-ion maxima, requires agreeing ion evidence
  plus local TIC prominence, returns apex hypotheses without quantitative hard
  borders, and builds padded fit regions using locally supported/varying ions.
- On the cached ten-minute segment, default candidates number 231 versus 612
  cleaned legacy blobs. The reviewed comparison band falls from 99.0% labelled
  pixels to **zero candidate centres**, while the crowded focus retains 21.
  This passes the background half of Phase 2's gate; sensitivity cannot pass
  until reviewed/known peaks are annotated.
- The joint shape fit now uses identical weights for its inner spectral NNLS and
  outer nonlinear objective, can fit a signed local background while keeping
  component spectra non-negative, supports
  repeated starts, reports the weighted residual, and counts fitted spectra in
  BIC. Integrated component signal now includes the fitted shape integral.
- `mcr_restarts` reports restart residuals, explained variance, convergence and
  optimally matched spectral stability.

At the candidate region nearest 15 min / 1.2 s, individual-ion apexes form two
clusters. Corrected BIC selects two Gaussians over one or three, and EMG does not
improve the result. The fitted centres are (897.08 s, 1.105 s) and (927.44 s,
1.180 s). The remaining TIC residual norm is 20.9% of the observed TIC, so the
fit is useful as a two-component hypothesis but is not a complete description.
Free-form MCR gives:

| Components | Explained variance | Minimum restart spectral similarity | Interpretation |
|---:|---:|---:|---|
| 1 | 96.46% | ~1.000 | inadequate reconstruction |
| 2 | 99.10% | ~1.000 | stable, materially better; all starts converged |
| 3 | 99.32% best | 0.796 | unstable; four of five starts did not converge |
| 4 | 99.35% best | 0.797 | unstable over-factorisation |

The stable two-factor result has different leading-ion orderings, with one led
by m/z 119/134 and the other by 115/91. These are mathematical factors, not
identified compounds. Three-fold spatial holdout RMSE decreases from 286.2 at
k=1 to 259.3 at k=2 and 250.3 at k=3. The folds are highly imbalanced and k=3
improves only 3.5% over k=2, while its MCR solution is unstable. The combined
evidence therefore supports **two stable factors**, with no defensible third.

The first quantitative recovery gate injects ideal Gaussian profiles carrying
the two real MCR spectra into an actual background-only patch. Equal-area pairs
recover within 10% at 0.06, 0.12 and 0.20 s separation. A 30%-abundance minor
component passes at 0.20 s, is borderline at 0.12 s (major error -10.8%), and
fails at 0.06 s. A 10%-abundance minor component fails at every separation; even
at 0.20 s its area is 31.5% high. That unrestricted fit met the 10% target in
4 of 9 cases; the current shared-width fast fit meets it in 5 of 9. Spectral
cosines remain high in the successful cases, confirming
that the failure is area attribution rather than simple non-convergence.

A second diagnostic uses the confirmed n-C16 profile from the full run. All
existing analytic shapes retain about 19% bias on a representative clean case.
`fit_template_overlap` instead shifts copies of the measured profile while
solving non-negative component spectra. It passes all nine signal-only cases.
The linear solver now permits signed background coefficients while keeping
component spectra non-negative, matching the signed residuals produced by
cleaning. With that correction, the movable template passes 6 of 9 cases in the
chosen real-background window: all equal-size and 30%-minor cases pass and all
10%-minor cases fail. A representative 30%-minor case also passes in the quietest
candidate-free windows at t2 5, 7 and 9 s, but fails at t2 3 s; no candidate-free
window of the tested size exists near t2 1.2 s. Background suitability must
therefore be part of the acceptance decision.

The run-wide route now uses complete-link candidate grouping, shared local Gaussian widths,
eight-ion component-count/shape screening, and one 30-ion production spectrum solve. The
cached ten-minute segment contains 231 candidate centres in 120 bounded regions and takes
about 13 s after startup. On the complete supplied 70 eV acquisition (215,015 scans), file
reading took about 27 s, threaded background correction 60 s, and four-process
detection/separation 68 s, for 154 s total. One-process separation took 118 to 131 s in
separate runs. The same end-to-end benchmark took 344 s before the performance pass.
Most crowded 14–16 min groups remain unresolved because the Gaussian reconstruction
explains too little local variation; faster execution has not been treated as evidence of
chemical correctness.

On the cached ten-minute segment, raising the fast optimizer's step limit from 12 to 24
reduced non-convergence flags from 37 to 10 with less than a second of extra fitting time,
but left 110 of 120 regions unresolved. Widening the second-dimension width limit from
three to six scan rows barely changed the median explained fraction (0.610 to 0.619)
and reduced known-area copy recovery from 5/9 to 1/9. The width cap remains three rows.

The rapid overlapping-Gaussian fit currently limits the second-dimension width to three
scan rows by default. That helps prevent a weak factor from expanding into background in
the known-area copy test, but isolated real peaks have a wider range of observed widths
(median about three rows; the upper quartile is substantially wider). Single-candidate
groups use their own observed width. Any overlapping fit that reaches the rapid width
limit is marked unresolved, and the limit is exposed in MS-aware settings. A library of
measured isolated shapes is needed to replace this conservative global bound.

See `component_metrics.json`, `component_fit.png`, `mcr_components.png`,
`injection_metrics.json`, `measured_profile_metrics.json` and
`measured_profile_recovery.png` under `real_data/blob_diagnosis/`. The automatic audit is
`automatic_component_metrics.json`. GUI persistence, contours, identification and CSV
export are implemented. Production quantification remains blocked until measured isolated
profiles and held-out peak annotations define signal-dependent acceptance limits.
