# Roadmap

What is left, in the order it is worth doing. Evidence for every item marked *(audit)* is in
[docs/validation/audit-2026-09-24.md](docs/validation/audit-2026-09-24.md). Delete an item when it
is done, and record it in CHANGELOG.md.

## 1. Data only James can supply

These block validation, not development.

1. **A real calibration sequence**: standards of known composition at several levels, with the
   internal standards spiked, run with the 8 s method. Everything in 🧪 has so far been checked
   only on planted peaks. Types and levels can now be set in the Runs table.
2. **Internal-standard concentrations**, and whether the mix always contains the five standards
   found in the current runs (acenaphthene-d10, phenanthrene-d10, n-tetracosane-d50,
   chrysene-d12, perylene-d12). n-Dotriacontane-d66 elutes after the end of these 76-minute
   runs and di-n-octyl phthalate-d4 was not found *(audit)*.
3. **A mass-spectral library** (NIST export via `lib2nist`, or MoNA/MassBank). Every 🏷️ threshold
   is tuned on small synthetic libraries.
4. **A labelled co-eluting standard pair** (BbF/BkF, benz[a]anthracene/chrysene or
   phenanthrene/anthracene) for the MS-aware peak finder, ✂️ and the retention-index penalty in 🏷️;
   component amounts must not feed quantification until that exists.
5. **Labelled siloxane peaks** from a few runs, to check the 0.15 threshold beyond the two runs it
   was chosen on.
6. **A vendor comparison**: one spectrum and one retention time read in ChromSpace/TOF-DS against
   PyGCGC's native reader.

## 2. Advice from the audit, no code change

- **🧪 Pick selective quantifier ions.** On busy ions (m/z 57, 105) planted targets were 15–60% off
  in real matrix; on m/z 136 and 212 within a few per cent. Qualifier ratios catch an interferent
  with a different spectrum, not one with a similar spectrum *(audit)*.
- **⛓️ Keep the qualifier gate.** Re-measured, n-alkanes here give 71/57 = 0.60 and 85/57 = 0.38,
  but the shipped gate (0.66, 0.44) gives better ladders than one centred on the new values,
  because the branched isomers it must reject have lower ratios still
  (`benchmarks/real_data_audit/alkane_ratios.py`).

## 3. Known gaps

- **MS-aware peaks do not reach every tool.** 🫧 blocking, ⛓️ stencil totals and series, 🗂️ batch
  and the 14 eV path of 🧬 still run the legacy detector. Legacy footprints cover most of the
  image on these busy runs, so anything summed over footprints (stencil totals, all-signals
  ratios) inherits that.
- **The ladder can take a wrong isomer.** One or two carbons per run are off by a modulation or
  more, typically at C11 and C22–C23 where a branched alkane with alkane-like ratios sits close by;
  once one carbon is wrong, the next predictions lean on it *(audit)*. A spacing check and a
  preference for the predicted position were both tried and did not help at an 8 s modulation.
- **Pristane and phytane** share a modulation column with n-C17 and n-C18 at long modulation
  periods, so their offsets are known only to about 0.05 carbon.
- **🧬 classes, 🧩 clusters and 🧪 results** are recomputed after reopening (fast once the cached
  blobs are back; settings and decisions are saved).
- **Blank subtraction** reads the blank through the alignment, but is unchecked on real data, is
  not redone when an alignment changes, is not cached, and is absent from 🗂️ batch and the CLI.
- **📍 alignment** is unchecked on real drift; clicks do not snap to a peak's centroid; batch
  anchor-finding matches each anchor independently rather than by assignment across every peak.
- **🔍 legacy co-elution flag** (single-linkage on each ion's apex) has never been checked against
  a real co-eluting pair.
- **🫧** The per-pixel method uses one reference only. Siloxane ions for 🧩's "leave out siloxane"
  are the bleed reference's.
- **🧬 rules** were checked on the earlier acquisition only.
- **🧩 clustering** has never seen a real multi-sample set; its merge distance (0.3) and region
  count (12) are starting points.
- **🗂️ batch**: the feature-matching radius and gap-fill window are fixed constants; progress is one
  status line.
- **🏷️ identification**: only Schymanski levels 2a/2b/3 are reachable; MSP parsing is checked only
  on hand-written files.
- **✂️**: hand-drawn line and valley splits exist in `core/overlap` with no drawing tool; PARAFAC
  across runs has no GUI; PARAFAC2 was never built; MCR-ALS unimodality is an approximation.
- **Detector saturation** is the brightest value seen three times; estimating it from flat-topped
  apexes, as James asked, awaits real saturated peaks.
- **Reports** never carry a real manual-edit log; the mechanism exists but no GUI action passes
  one.
- **Semi-quantitative amounts** assume an unknown responds like its nearest calibrated target.
- **Packaging**: `packaging/pygcgc.spec` has never been built.

## 4. Suggested features

In rough order of value for the work these runs are for.

1. **Quantifier-ion advisor.** For each target, rank its ions by how much of the local signal on
   that ion belongs to the target across the calibration runs, and suggest the most selective one.
2. **Suggested alkane seeds from the 14 eV channel.** On tandem data the n-alkane molecular ion
   (14n + 2) names the carbon number of a peak on the alkane line. It cannot tell an n-alkane from
   a branched isomer of the same mass, and its signal is weak (1–3% of base at 14 eV), so it should
   propose two seeds for the user to accept, not replace the clicks. In the audit it confirmed 7 to
   11 carbons per run up to C21 and was silent above.
3. **A "report this project" action** that passes the project's review decisions, blob edits and
   stencils into the provenance of an Excel/HTML report.
4. **A settings dialog for 🧬** and a class filter for 🏷️ search.
5. **Split `gui/mainwindow.py`** (about 3,000 lines holding every step's handlers) into one
   controller per step; it is the file every change touches and where parallel edits collide.

Done on 24 September 2026 (see CHANGELOG.md): siloxane scoring against both references at 0.15,
blocking siloxane ions only; run type, level, calibration set and volume in the Runs table; saved
alkane ladders and a results cache for background layers and detections; native reading 6–12 times
faster; the qualifier ratios re-measured (defaults kept, section 2).
