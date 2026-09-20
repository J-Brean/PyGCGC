# MS-aware overlapping peak detection: full implementation plan

20 September 2026.

## 1. Intended result

PyGCGC should detect a list of chromatographic components. Each component represents
one proposed compound signal and contains:

- a position in first- and second-dimension retention time;
- a soft two-dimensional peak shape;
- a mass spectrum;
- an estimated integrated intensity;
- uncertainty and confidence information;
- a link to the local group in which it was separated.

Component shapes may overlap. A measured pixel can contain contributions from several
components. Analysis rectangles remain an internal way to divide the full chromatogram
into manageable calculations; they are not compounds and should not be the primary
result shown to the user.

The legacy watershed detector remains available as a separate method. It produces
exclusive, hard-edged `BlobSet` footprints and continues to support old projects and
comparison work.

## 2. Plain-language model

At every position in a GC×GC-MS run, the instrument measures a mass spectrum: a list
of intensities at different mass-to-charge values. The new detector explains a small
piece of the run as a sum of compounds:

```text
measured data = peak A shape × peak A spectrum
              + peak B shape × peak B spectrum
              + ...
              + local background
```

Peak A and Peak B may occupy the same pixels. Their different spectra, positions and
chromatographic shapes provide the information needed to estimate their separate
contributions. The software should display those contributions as overlapping contours,
not divide the image into mutually exclusive territories.

## 3. What already exists

The following foundations are implemented and tested:

- `core/candidates.py` finds per-ion peak evidence, groups ions that agree on an apex,
  and makes padded local analysis regions without hard borders.
- `core/overlap.py` builds local ion cubes, fits multiple peak shapes and spectra,
  supports signed constant or bilinear background, chooses component count with BIC,
  predicts fitted data, reports positional uncertainty, and produces overlapping
  component contribution maps.
- A measured-peak template fitter can move copies of an observed alkane peak rather
  than assuming an ideal Gaussian.
- MCR-ALS and repeated-start stability diagnostics provide an independent decomposition.
- Weak fitted components are retained and marked low confidence rather than silently
  discarded.
- Workflow step 4 offers MS-aware detection beside legacy blobs. Candidate centres and
  analysis rectangles are shown in the main view, and a selected region can be opened
  in the overlap dialog.

The missing layer is a durable, automatic component pipeline connecting these pieces.

## 4. Terminology

**Candidate centre** — evidence that one or more ions have a local peak at approximately
the same position. It is a starting hypothesis, not yet a separated compound.

**Analysis region** — a padded local rectangle containing one or more candidate centres.
It limits the amount of data sent to a fit. It is not a peak boundary.

**Component** — one fitted chromatographic shape multiplied by one fitted mass spectrum.
This is the proposed replacement for a blob in the MS-aware route.

**Component group** — all components fitted jointly in one local analysis region. A
group records shared background and fit-quality information.

**Contribution map** — the intensity attributed to one component at every pixel. Maps
from different components can overlap and add together.

**Unresolved group** — a group for which the software can estimate contributions but
cannot support a reliable claim that every component is individually resolved.

## 5. Component data model

Add `core/components.py`. Keep its core classes free of Qt.

### `Component`

The durable result for one proposed peak:

```python
@dataclass
class Component:
    id: str
    group_id: str
    t1: float
    t2: float
    volume: float
    spectrum_mz: np.ndarray
    spectrum_intensity: np.ndarray
    model: str
    shape_parameters: dict[str, float]
    position_se_rows: float
    position_se_columns: float
    relative_to_largest: float
    confidence: str
    reason: str
```

The fitted spectrum should use the run's actual m/z values. Store only nonzero channels
above a small documented relative threshold so project files stay manageable. Preserve
the unthresholded spectrum in memory while a project is open.

`volume` is the integral of that component's contribution map. It is not the sum of a
hard footprint. It may therefore coexist with volume from another component at the same
pixels.

### `ComponentGroup`

The record for one joint fit:

```python
@dataclass
class ComponentGroup:
    id: str
    region: CandidateRegion
    component_ids: tuple[str, ...]
    model: str
    template_source: str | None
    background_order: int
    bic: float
    weighted_residual: float
    explained_fraction: float
    converged: bool
    stable: bool
    unresolved: bool
    warning: str
```

Do not serialize the full contribution arrays. Reconstruct them from the stored model,
shape parameters and spectrum when possible. For free-form MCR results, either store a
compressed local profile or treat MCR as a diagnostic rather than the accepted durable
result in the first implementation.

### `ComponentSet`

One run's complete MS-aware result:

```python
@dataclass
class ComponentSet:
    components: list[Component]
    groups: list[ComponentGroup]
    candidates: CandidateSet
    settings: ComponentDetectionSettings
    source_fingerprint: str
```

It needs helpers equivalent to the useful parts of `BlobSet`: length, nearest component,
sorting by retention position, spectra, and contour generation. It must not expose a
single integer label image, because that would imply exclusive ownership.

## 6. Settings

Introduce `ComponentDetectionSettings`, containing three clearly separated groups:

1. Candidate settings (`CandidateSettings`).
2. Region-building settings: link distance and padding in t1/t2.
3. Separation settings: model family, maximum component count, background order,
   repeated starts, confidence thresholds, and quality gates.

Initial defaults:

- use candidate centres to initialise component positions;
- test one through `min(number of plausible centres, max_components)` components;
- maximum four components per group initially;
- signed bilinear local background;
- five starts per fit;
- BIC for preliminary component count;
- component below 30% of the largest is retained but marked low confidence;
- failed convergence, poor spatial cross-validation or unstable restarts mark the whole
  group unresolved;
- no automatic deletion of a weak component solely because it is weak.

The 30% boundary is provisional and must remain visible in settings and provenance.

## 7. Automatic separation pipeline

Implement a pure-core entry point:

```python
detect_components(run, fold, settings, floor=None, template_library=None,
                  progress=None, cancel=None) -> ComponentSet
```

### Stage 1: candidate detection

Run the existing ion-wise detector. Keep the supporting ion lists and S/N values; they
will initialise spectra and prevent a bright unrelated ion from dominating a weak
candidate.

### Stage 2: construct local groups

Build padded regions, then merge rectangles that overlap substantially or share candidate
centres. Avoid fitting the same physical neighbourhood twice. Keep regions small enough
for local-background assumptions to remain reasonable.

### Stage 3: initialise possible components

For every candidate centre:

- initialise position from its apex;
- initialise its spectrum from a small local, background-subtracted neighbourhood;
- weight its independently supporting ions more strongly;
- estimate starting widths from local t1 and t2 profiles.

Candidates are hypotheses. Closely coincident candidates with nearly identical spectra
should be allowed to merge into one component. A region may also require fewer components
than candidate centres.

### Stage 4: choose and fit the model

The first production-capable route should use an analytic model that is available on
every run. Start with `emg2` or whichever existing family wins held-out real-peak tests.
Measured templates should be preferred only when a suitable local shape library exists.

For each plausible component count:

1. perform repeated-start weighted fitting;
2. fit non-negative component spectra and signed local background;
3. calculate BIC and spatially blocked prediction error;
4. reject numerically duplicated components;
5. compare matched components across starts.

Select the simplest fit that has adequate predictive performance and stable components.
BIC alone must not be the final acceptance rule.

### Stage 5: assess identifiability

For each component record:

- relative amount;
- position uncertainty;
- spectral distinction from neighbouring components;
- amount and position variation across restarts;
- sensitivity to background model;
- whether the component is supported by independent ion apexes.

For the whole group record:

- convergence;
- weighted residual and explained fraction;
- spatial cross-validation error;
- restart stability;
- residual structure;
- whether an alternative component count is nearly as plausible.

The group is unresolved if these checks disagree. Keep its component estimates and show
the reason. Do not silently collapse the group into one peak or silently report uncertain
components as fully resolved.

### Stage 6: reconcile neighbouring groups

Regions at chunk or grouping boundaries may produce duplicate components. Match proposed
components using position, spectral cosine similarity and shape overlap. When two results
represent the same component:

- retain the fit from the group where it lies farther from the region edge and has better
  predictive error; or
- refit the union of both groups if neither result is clearly preferable.

Give component IDs deterministically from run identity plus rounded position and a stable
within-group index. IDs must survive project reloads and small unrelated edits.

### Stage 7: final residual search

Subtract the accepted model prediction locally and run a conservative candidate search on
positive residuals. Add a new candidate only when several ions agree and the residual is
larger than model/background uncertainty. Refit that group once. This catches peaks hidden
under a stronger component without creating an unlimited fit-add-fit loop.

## 8. Visual design

### Main 2D view

The normal MS-aware overlay should show components rather than analysis rectangles:

- one colour per component within a local group;
- a solid contour at 10% of that component's own maximum;
- an optional dotted contour at 50%;
- a small marker at the fitted centre;
- translucent contribution shading when zoomed in;
- overlaps remain visible through colour blending.

Analysis rectangles should be hidden by default. Add a diagnostic toggle such as
**Show analysis regions** for developers and advanced review.

Do not draw a single exterior border around a group. That would visually turn the group
back into a blob.

### Peak table

In MS-aware mode, one row should represent one component, not one region. Suggested columns:

- component number;
- t1 and t2;
- estimated intensity;
- relative size within its group;
- strongest or diagnostic ions;
- confidence;
- group state: resolved or unresolved;
- identification, when available.

Rows belonging to the same group may share subtle background shading or a group label.
Low-confidence components remain visible in orange. Selecting a row highlights only that
component's contours. A separate action selects or opens the whole group.

### Detail view

Rename the current overlap window to **Component group**. It should show:

- observed TIC;
- one contribution map per component;
- the summed fitted signal;
- residual signal;
- component spectra side by side;
- confidence reasons in plain English;
- alternative component-count diagnostics;
- controls for refitting, merging, splitting, adding and suppressing components.

### Manual review operations

All operations must be undoable and recorded as decisions layered over a reproducible
automatic result:

- accept or reject a proposed component;
- retain a component but mark it unresolved;
- add a centre and refit the group;
- merge two fitted components and refit;
- request a split with two seed positions;
- change the local model or background order;
- choose among plausible fits.

Avoid direct painting of exclusive pixel ownership. If a hand-drawn area is supported,
use it as a fitting mask or initialisation hint.

## 9. GUI and state architecture

Rename `BlobState` to `PeakDetectionState` after compatibility shims are in place. During
migration it may continue to expose `by_run` for legacy blobs, but new names should be:

```text
legacy_by_run: RunEntry -> BlobSet
components_by_run: RunEntry -> ComponentSet
method: "ms-aware" | "legacy"
component_decisions: run name -> edit log
```

The state runner should have two explicit entry points internally:

- `run_legacy_all()`;
- `run_components_all()`.

The generic `run_all()` may dispatch from the selected method for the workflow button.
Downstream legacy callers must request legacy data explicitly during migration rather than
changing the user's selected display method.

Separation is much more expensive than candidate finding. Submit one cancellable job per
run and report region-level progress. Cache completed groups as they finish so cancellation
or one failed group does not discard an otherwise useful run.

## 10. Persistence and reproducibility

Project files should save:

- selected method;
- all candidate/component settings;
- source fingerprints for the raw run, fold, background and alignment;
- component and group parameters;
- confidence results;
- manual review decisions;
- implementation version for the component model.

Do not save large dense data cubes. Recompute contribution maps on demand from stored fit
parameters. If the source fingerprint changes, mark the component result stale and offer
to rerun it. Never show a fit from an old fold or background layer as current.

Old projects without component fields load as legacy projects. New projects retain both
legacy results and components when both have been run, so users can compare methods without
recalculation.

## 11. Downstream migration

Migrate one consumer at a time. Keep legacy code paths available until each component path
has validation coverage.

### Identification

This is the most direct first consumer. Search each component's fitted spectrum rather than
an apex spectrum from a mixed blob. Report group confidence beside library-match confidence;
a strong library match does not prove that the chromatographic separation was reliable.

### Non-target feature table and batch comparison

Match components across runs using aligned t1/t2, spectrum similarity and uncertainty.
Use component volume as the feature intensity. A component in an unresolved group remains
available but carries that flag into exports and statistics.

### Siloxane handling

Score fitted component spectra against siloxane references. Suppress or subtract the
siloxane component contribution rather than zeroing every ion in an entire hard blob.
Do not automatically suppress a co-eluting analyte component merely because its group also
contains siloxane.

### Hydrocarbon tools

Use component centres and spectra for homologous-series detection and stencil totals.
The direct-ion alkane ladder can remain independent because it already searches ion images
without requiring blobs.

### Quantification

Target quantification should use a component only when its identity, position and qualifier
ions agree with the target. Calibration and sample runs must use the same measurement rule.
Unresolved targets retain an estimated amount but receive an explicit uncertainty/status
flag. Do not silently replace the existing target-specific integration until known-mixture
recovery demonstrates equal or better bias and precision.

### Reporting and export

Export one row per component with run, group, position, amount, diagnostic ions, confidence,
unresolved reason and identification. Provide a group table and enough provenance to
reproduce the fit.

## 12. Performance plan

The all-ion run is too large for one global tensor decomposition. Keep processing local and
sparse:

- read only ions needed by each region, prioritising candidate-support ions and high local
  variation;
- process independent non-overlapping groups in parallel;
- limit simultaneous dense cubes by a memory budget;
- cache candidate results and completed component groups separately;
- reuse local noise/background estimates;
- defer dense contribution-map construction until display or export;
- record progress as regions completed, with cancellation between starts and groups.

Benchmark candidate detection, fitting time per group, peak memory and complete-run time on
the supplied real acquisition before choosing default worker counts.

## 13. Validation programme

### Unit and synthetic tests

- one isolated component is returned once with correct position, spectrum and volume;
- two spatially overlapping components with distinct spectra are returned separately;
- identical spectra with insufficient chromatographic separation are marked unresolved;
- shared ions do not prevent separation when other diagnostic ions exist;
- negative cleaned-data background is fitted without creating negative component spectra;
- duplicate components across neighbouring regions are reconciled;
- weak estimates are retained and flagged;
- project round-trip preserves IDs, parameters, confidence and manual decisions;
- changing fold/background marks results stale;
- legacy results remain unchanged.

### Software-copy recovery tests

Extend the existing known-area benchmark across:

- t1 and t2 spacing;
- minor/major ratios from equal down through 30%, 20% and 10%;
- total signal strength;
- spectral similarity;
- background type and strength;
- peak-shape mismatch;
- region-edge position.

Report bias, variance, missed components, extra components and unresolved calls. A correct
unresolved call is preferable to a confidently wrong separation.

### Held-out real-data tests

Create a small reviewed set of windows containing:

- isolated clear peaks;
- known co-eluting standards;
- crowded hydrocarbon regions;
- column bleed and other background;
- weak shoulders;
- wrapped peaks near the modulation boundary;
- saturated or clipped peaks.

Reviewers should label presence, approximate centres and any known identities without using
the fitted output as the answer. Keep tuning windows separate from acceptance windows.

### Acceptance gates

Before making MS-aware separation the default:

- false component rate is acceptably low on reviewed background;
- isolated peaks are neither split nor missed more often than legacy detection;
- component count is stable across starts and small setting perturbations;
- reported confidence predicts actual recovery quality;
- typical full runs finish within an agreed time and memory budget.

Before using component volumes for production quantification:

- known-mixture components above the declared operating limit meet the chosen bias and
  precision target across held-out backgrounds;
- results below that limit are reliably marked low confidence or unresolved;
- calibration and sample measurements remain consistent;
- exported uncertainty and unresolved state cannot be mistaken for a fully resolved amount.

## 14. Delivery phases

### Phase A — durable component result

1. Add `Component`, `ComponentGroup`, `ComponentSet` and settings.
2. Convert one existing `OverlapFit` into those objects.
3. Implement serialization, source fingerprints and staleness.
4. Add reconstruction of contribution maps and contours.

**Exit:** one manually fitted region can be accepted, saved, reopened and displayed as
overlapping components with the same amounts and spectra.

### Phase B — automatic run-wide separation

1. Implement `detect_components`.
2. Use candidate centres for initial positions and spectra.
3. Add component-count selection, quality gates and restart stability.
4. Reconcile neighbouring regions and run one residual search.
5. Add progress, cancellation and group-level caching.

**Exit:** pressing Find peaks in MS-aware mode returns a component list for every loaded
run without opening each region manually.

### Phase C — component-first interface

1. Change the main table from regions to components.
2. Draw overlapping component contours by default.
3. Hide analysis rectangles behind a diagnostic toggle.
4. Build the Component group review window and undoable edit log.

**Exit:** users interact with overlapping peaks; regions appear only when inspecting how a
group was calculated.

### Phase D — identification and batch features

1. Search component spectra against libraries.
2. Match components across runs using position plus spectrum.
3. Export component and group tables with confidence.

**Exit:** untargeted identification and feature comparison no longer require legacy blobs.

### Phase E — siloxane, hydrocarbon and quantification migration

1. Remove fitted siloxane contributions component-wise.
2. Read component spectra/amounts in hydrocarbon diagnostics.
3. Add an explicitly enabled component-based target measurement route.
4. Validate it against calibration and held-out known mixtures before changing defaults.

**Exit:** the full workflow can operate without hard blob ownership, while legacy mode
remains available for reproducibility.

## 15. Immediate implementation order

The next coding work should be Phase A, in this order:

1. Create the Qt-free component dataclasses and conversion from `OverlapFit`.
2. Add component/group storage to peak-detection state without changing current legacy
   behaviour.
3. Add a deliberate **Accept fit** action to the component-group window.
4. Replace the selected region's crosses/rectangle with stored overlapping component
   contours after acceptance.
5. Save and reload accepted component results and manual decisions.
6. Only then automate fitting every region, because automatic output needs a stable result
   format, visual review and persistence before it can be scientifically audited.

This sequence makes the desired data product concrete early and prevents temporary analysis
rectangles from becoming an accidental permanent interface.
