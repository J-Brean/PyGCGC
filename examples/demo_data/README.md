# PyGCGC demonstration data

Seven synthetic GC×GC-EI-ToF-MS runs, made to show every part of PyGCGC: four calibration
standards and three samples, with the supporting files the app reads and the true answer for
every number. The compounds are real; their retention indices and mass spectra are
approximations written for this purpose, not library spectra. Do not use these files to validate
a method.

Rebuild the pack (about 20 s, 100 MB; the fixed seed gives identical files):

```powershell
python examples/demo_data/generate_demo_data.py "C:\path\to\SYNTHETIC_DEMO"
```

## What is in the folder

| File | What it is |
|---|---|
| `data/CAL_L1_1ng.cdf` … `CAL_L4_50ng.cdf` | calibration standards at 1, 5, 20 and 50 ng per compound (decalin, the cycloalkanes and behenic acid at 5×, the oxy- and nitro-PAHs at 0.1×) |
| `data/SAMPLE_1_urban_air.cdf` | traffic and indoor-air mixture: BTEX smeared up t₂ by an overloaded toluene, alkanes, PAHs, aldehydes, terpenes, a large D5 siloxane, phthalates, fragrance, nicotine |
| `data/SAMPLE_2_diesel_exhaust.cdf` | alkanes C9 to C30 on an unresolved hump of 700 minor isomers, pristane and phytane, cycloalkanes, alkylbenzenes, methylnaphthalenes, PAHs, oxy- and nitro-PAHs |
| `data/SAMPLE_3_wood_smoke.cdf` | phenol, cresols, guaiacols, syringol, a tailing levoglucosan, retene, resin and fatty acids, methyl esters, odd-over-even plant-wax alkanes |
| `method.json` | the 8 s modulation period |
| `sequence.csv` | run type, calibration level and sampled volume (m³) of every run |
| `target_table.csv` | 34 targets and five internal standards, with quantifier and qualifier ions, for 🧪 Quantify |
| `calibration_levels.csv` | amount of every target at every level, for 🧪 Quantify |
| `calibration_compounds.csv` | the same compounds in the layout the 📈 Calibrate window reads (one amount column per level) |
| `demo_library.msp` | MSP library of all 150 named compounds with retention indices, for 🏷️ Identify |
| `truth/truth_amounts.csv` | every compound in every run: class, amount (ng), true t₁ and t₂ |
| `truth/truth_runs.json` | each run's phase, retention drift, detector ceiling, noise and threshold |

Every run carries five deuterated internal standards at 10 ng (acenaphthene-d10,
phenanthrene-d10, n-tetracosane-d50, chrysene-d12, perylene-d12) and the same tube background
(D4 to D8 siloxanes, dibutyl phthalate, DEHP, BHT, 2,4-di-tert-butylphenol). The standards also
carry a dichloromethane solvent front.

Problems put in on purpose, for the app to find:

- every run starts at its own phase (1.6 to 2.5 s), so the image is folded wrongly until 🔄 Fold finds it;
- retention drifts up to 3.6 s in t₁ and 0.03 s in t₂ between runs (📍 Align);
- column bleed rises with the oven, a modulator pulse leaves a stripe at the bottom of every modulation, and 80 single-scan spikes are scattered through each run (🧹 Background);
- perylene and perylene-d12 elute after the end of the modulation and wrap round to the bottom of the next one;
- phenanthrene and anthracene sit close together, and benz[a]anthracene and chrysene overlap (🧪 flags them);
- the detector bends over above about 2×10⁵ counts and saturates at 4×10⁵, so the 50 ng points of the strongest compounds fall below the line;
- quantifying alkanes on m/z 57 over the diesel hump over-reads them, as it does on real exhaust.

## A tour of the app

Open the folder: `pygcgc-gui "C:\path\to\SYNTHETIC_DEMO"`. The period comes from
`method.json`, the run types from `sequence.csv`, and the phase of every run is found from the
internal standards in `target_table.csv` as the runs load. Work along the numbered toolbar.

1. **🔄 Fold.** Select the diesel sample. The fold panel says the phase came from five anchors.
   Drag the phase slider to see a wrong phase split the peaks across the top and bottom edges,
   then Ctrl+Z.
2. **📍 Align.** Set the reference to `CAL_L3_20ng`. Select the urban sample, press ➕ Add
   anchor (A) and click the same peak on the reference and on the sample (phenanthrene-d10 at
   about 22 min, 4.0 s, is easy); add four or five spread across the run. 📍 Auto-find on all
   runs finds them elsewhere. The Compare tab (5) blinks, overlays or subtracts the runs.
3. **🧹 Background.** Press ▶ Remove background. Type `73` in the Ion box, select the wood-smoke
   sample, zoom on 33 to 36 min and tick *Before / after split*: behenic acid stands out of the
   bleed once it is removed. Space toggles raw and cleaned; S/N shows signal-to-noise.
4. **🔍 Peaks.** ▶ Find peaks (MS-aware separation). Zoom into the crowded middle of the urban
   sample: overlapping peaks get overlapping contours and a row each in the table.
   Tilted Gaussians with independent widths are fitted by default. **▶ How it works** plays
   a short synthetic explanation with narrow detector-time pulses, an enlarged modulation
   cycle and its matching 2D column; pause, scrub or replay it at any time.
   **Preview fitting…** tries the current region or run without replacing accepted peaks.
   Inspect its fitted sum and residual, then apply the settings to all runs when satisfied.
   - **🫧 Siloxane filterer…** (button under the table): the large D5 in the urban sample (10.5 min,
     0.7 s) and the D4 to D8 series light up in dashed magenta at the default 15% certainty.
   - **⛓️ Hydrocarbon identifier…**: select the diesel sample, Ion `57`, press ➕ Find alkane series
     and click n-C14 (15.0 min) and n-C20 (25.5 min), typing 14 and 20. The ladder grows to
     about C9 to C30; check its ends, where the qualifier gate can take a minor isomer for an
     extra carbon. *Propagate to other runs*, then *Assign peaks to hydrocarbons*.
5. **🏷️ IS.** Select `CAL_L1_1ng`, where the native compounds beside the standards are smallest.
   For each internal standard, select it in the table, press 🏷️ Label peak on this run and click it
   (positions in `truth/truth_amounts.csv`). Go finds all five in every run; the Standards tab
   shows the grid. 📋 Use as target table writes them as `istd` rows.
6. **🧪 Quantify.** Go. The Review tab (3) shows each target across the sequence: pick
   benz[a]anthracene to see it flagged for sharing a peak with chrysene. The Calibration tab (4) shows
   every curve with its relative-error plot; click the 50 ng point of naphthalene and Exclude it.
   The Results tab lists targets, all signals, internal standards and calibrations; compare the
   amounts with `truth/truth_amounts.csv`.
7. **🧬 Classify.** On the wood-smoke sample: phenols, PAHs, acids, esters and alkanes are coloured by
   group; select a row to read why it was called what it was.
8. **🏷️ Identify.** Load `demo_library.msp` and identify the peaks of the urban sample. Hits come with
   a retention-index check (once the alkane ladder exists) and a target-decoy q-value.
9. **🧩 Cluster.** Cluster the diesel sample and fill the peaks with their cluster colour: the
   alkanes, cycloalkanes, alkylbenzenes, naphthalenes and PAHs separate without any rules. Pool every
   run to compare the three samples cluster by cluster.

Elsewhere: ✂️ Resolve overlap opens a peak group in 3D and refits it; 🗂️ Batch runs the whole folder
in parallel and gives a feature table, QC chart, PCA and volcano plot, and exports Excel, HTML or PDF;
📤 Figure saves the 2D view as a 300 dpi figure with a caption file.

### The calibration table (📈 Calibrate)

The main eight-step workflow also offers **Set up calibration and units…** under Assign and
Quantify. It imports the existing target/level inputs into a draft, connects target-to-IS
assignments, checks calibration runs and usable levels, and records actual constant IS amounts
and an output unit. The demo's supplied truth is the source for its amounts; private runs still
need their real spike amounts. Back and Cancel leave active inputs unchanged.

After searching standards, **Review missing or uncertain matches…** walks through exceptions
with candidate spectra and fitted profiles side by side. Reviewed positions feed quantification;
confirmed misses retain missing-standard flags. The workflow explains prerequisites and which
results need rerunning. The separate builder below remains an advanced route for curve editing
and amount-ratio calibration when IS additions vary.

1. In the Runs table select the four `CAL_` runs (Ctrl+click), then press **📈 Calibrate** on the
   toolbar and *Use the 4 runs selected in the main window*. The standards are cleaned, their peaks
   found and their phases lined up; the level amounts (1, 5, 20, 50) are read from the run names.
2. On **2 Compounds**, press *📄 Load a compound table…* and choose `calibration_compounds.csv`.
   Every compound is matched to a peak in every standard and divided by its nearest internal
   standard. (Or press *🔎 Suggest from the standards…* to see what the app finds with no table:
   the five internal standards are recognised by their spectra.) Select a compound to see its
   peak outlined; click another peak to rematch it in one run.
3. On **3 Curves**, each compound has its curve and points. Click the 50 ng point of naphthalene and
   *Leave out*; retype an amount; change weights or model. Everything refits at once and can be undone.
4. On **4 Apply to samples**, tick the three samples and ▶ Quantify. Each compound is found by position
   and spectrum and read off its curve; amber cells say why a result is flagged. 📤 Export results
   writes CSV. Session ▾ saves the whole session to carry on later.

The same steps without the window are in `examples/tutorials/07_calibration_builder.py`.
