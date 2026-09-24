# Roadmap

What is left, in the order it is worth doing. Evidence for every item marked *(audit)* is in
[docs/validation/audit-2026-09-24.md](docs/validation/audit-2026-09-24.md). Delete an item when it
is done, and record it in CHANGELOG.md.

## 1. Data only James can supply

These block validation, not development.

1. **A real calibration sequence**: standards of known composition at several levels, with the
   internal standards spiked, run with the 8 s method. Everything in 🧪 has so far been checked
   only on planted peaks.
2. **Internal-standard concentrations**, and whether the mix always contains the five standards
   found in the current runs (acenaphthene-d10, phenanthrene-d10, n-tetracosane-d50,
   chrysene-d12, perylene-d12). n-Dotriacontane-d66 elutes after the end of these 76-minute
   runs and di-n-octyl phthalate-d4 was not found *(audit)*.
3. **A mass-spectral library** (NIST export via `lib2nist`, or MoNA/MassBank). Every 🏷️ threshold
   is tuned on small synthetic libraries.
4. **A labelled co-eluting standard pair** for the MS-aware peak finder; its amounts must not feed
   quantification until that exists.
5. **A vendor comparison**: one spectrum and one retention time read in ChromSpace/TOF-DS against
   PyGCGC's native reader.

## 2. Recommended changes, with evidence, awaiting a yes

1. **🫧 Stop blocking every ion.** At the default threshold on ATWC, blocking whole legacy blobs
   removed 4.5% of all m/z 57 and 1.9% of m/z 91 signal while catching 6 of 38 siloxane peaks: the
   flagged blobs are mostly hydrocarbon peaks with bleed underneath, and legacy footprints cover
   85% of the image. Recommend `block_all_ions=False` by default, or the per-pixel method, which
   removed 61–76% of m/z 207 with no loss on 57, 91 or 128 *(audit)*.
2. **🫧 Score against both references and lower the threshold.** The bleed reference (207-based)
   misses the cyclic siloxane peaks (73/281-based) and the reverse. Scoring each blob against both
   and taking the larger, at a threshold of 0.15, caught 38 of 38 and 21 of 21 siloxane-like
   blobs on two runs with no clearly non-siloxane blob flagged; blocking only the siloxane ions of
   those blobs removed 69–95% of m/z 207 and 281 and nothing measurable on m/z 57, 91 or 128
   *(audit)*.
3. **⛓️ Re-measure the alkane qualifier ratios for this instrument tune.** The defaults (0.66,
   0.44) came from the earlier acquisition; here n-alkanes give 0.55–0.66 and 0.30–0.43, and the
   gate drops C12 on one run and a late carbon on others. With the gate off the ladders are longer
   and no worse *(audit)*. A per-project "pick from a clicked alkane" already exists; the question
   is the shipped default.
4. **🧪 Pick selective quantifier ions.** On busy ions (m/z 57, 105) planted targets were 15–60%
   off in real matrix; on m/z 136 and 212 within a few per cent. Nothing in the code is wrong; the
   user guide now says so. See feature suggestion 3 below.

## 3. Known gaps

- **MS-aware peaks do not reach every tool.** 🫧 blocking, ⛓️ stencil totals and series, 🗂️ batch
  and the 14 eV path of 🧬 still run the legacy detector. Legacy footprints cover most of the
  image on these busy runs, so anything summed over footprints (stencil totals, blob blocking,
  all-signals ratios) inherits that.
- **The ladder can take a wrong isomer.** One or two carbons per run are off by a modulation or
  more, typically at C11 and C22–C23 where a branched alkane with alkane-like ratios sits close by
  *(audit)*. A spacing check was tried and removed: at an 8 s modulation, quantised positions
  make the spacing too noisy to judge.
- **Results are not saved with a project**: 🔍 blobs, ⛓️ ladders, 🧬 classes, 🧩 clusters and 🧪
  results are recomputed after reopening (settings and decisions are saved).
- **Run type and level can only come from `sequence.csv`**; there is no way to mark a run as a
  calibration level in the app.
- **Blank subtraction** reads the blank through the alignment, but is unchecked on real data, is
  not redone when an alignment changes, and is absent from 🗂️ batch and the CLI.
- **📍 alignment** is unchecked on real drift; clicks do not snap to a peak's centroid.
- **🫧, ⛓️ and 🧬 defaults** were measured on the earlier acquisition; the 🧬 rules have not been
  rechecked on the current runs.
- **🧩 clustering** has never seen a real multi-sample set; its merge distance (0.3) and region
  count (12) are starting points.
- **🏷️ identification**: only Schymanski levels 2a/2b/3 are reachable; MSP parsing is checked only
  on hand-written files.
- **Reports** never carry a real manual-edit log; the mechanism exists but no GUI action passes
  one.
- **Packaging**: `packaging/pygcgc.spec` has never been built.
- **Semi-quantitative amounts** assume an unknown responds like its nearest calibrated target.
- **✂️** PARAFAC across runs has no GUI; PARAFAC2 was never built; MCR-ALS unimodality is an
  approximation.

## 4. Suggested features

In rough order of value for the work these runs are for.

1. **Mark runs in the app.** Edit type, level, calibration set and volume in the Runs table and
   write `sequence.csv`, so a calibration can be set up without a spreadsheet.
2. **Automatic alkane seeding from the 14 eV channel.** On tandem data the n-alkane molecular ion
   names the carbon number, so the ladder could start with no clicks; the audit's ladder check
   already does this for verification.
3. **Quantifier-ion advisor.** For each target, rank its ions by how much of the local signal on
   that ion belongs to the target across the calibration runs, and suggest the most selective one.
4. **Save detection results in the project** (blobs, components, ladders, classes), so reopening
   a large project does not recompute minutes of work.
5. **Faster native reading.** `read_lsc` decodes each scan in Python (about 15 s per 460,000-scan
   file); reading the ion records in one vectorised pass would make it near-instant.
6. **A "report this project" action** that passes the project's review decisions, blob edits and
   stencils into the provenance of an Excel/HTML report.
7. **A settings dialog for 🧬** and a class filter for 🏷️ search.
8. **Split `gui/mainwindow.py`** (about 3,000 lines holding every step's handlers) into one
   controller per step; it is the file every change touches and where parallel edits collide.
