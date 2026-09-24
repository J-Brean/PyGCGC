# Real-data checks: compound-class rules

21 September 2026. Untreated 70 eV native file, modulation period 10 s, phase 1.81 s, m/z 38 to 300, 4068 blobs,
n-alkane ladder C12 to C28 (17 carbons, from seeds C13 and C18). Every number below is printed by
`examples/classify_real_data.py`, which also writes `class_map.png` and its caption. No library and no labelled
ground truth exist for this run apart from the alkane ladder, so most checks are of chromatographic consistency,
not identity.

## What was checked, and how far to trust it

| Check | Result | Strength |
|---|---|---|
| Alkane ladder recovered (blob with most m/z 57 at each carbon) | 17 of 17 | In-sample: the alkane thresholds (alkyl share 0.48 to 0.59, olefinicity 0.24 to 0.31) were read from this ladder |
| Spectrum-only alkanes fall on the alkane ridge | 52 blobs, t2 offset median -0.02 s (p10 to p90 -0.06 to 0.01) | Good: no alkane criterion uses t2, and about 35 of the 52 are not ladder members |
| Polarity ordering of classes by t2 offset | alkane -0.02 < alkene 0.20 < indane 1.6 < naphthalene 1.8 < PAH 2.1 < aromatic oxygenate 2.4 s | Consistent with chemistry, but `dt2` is a "prefer" term in most classes, so not independent |
| Deuterated internal standards | acenaphthene-d10 (164), chrysene-d12 (240, one of two blobs), perylene-d12 (264) classify as PAH; fluorene (166) and biphenyl or acenaphthene (154) also PAH | Independent, small |
| 26 synthetic reference spectra typed from familiar EI patterns | every enabled class fires on its own reference, except primary alcohols (below) | Tests the logic, not the thresholds |

## Findings that changed the design

1. **Empty apex spectra.** 319 of 4068 blobs (11.7% of the volume, median S/N 316) have an all-zero apex spectrum
   although their footprint spectra are fine (blob 2824's is an n-alkane, the C21 position). `classify_blobs` falls
   back to the footprint. `identify_blobs` does not, so those blobs are silently unsearched there. The cause is in
   blob detection and was not investigated.
2. **The sulfur class is disabled.** M+2/M of hydrocarbon-like blobs (M the base peak, M 90 to 260) rises with peak
   size: median 0.007 for 3e5 to 3e6 volume, 0.03 for 3e6 to 1e7, 0.05 for 1e7 to 3e7, 0.06 above. Indene,
   methylindenes, naphthalene and methylnaphthalene therefore read as sulfur (12.6% of the volume with the class
   enabled). No per-spectrum count threshold separates them, since a clean acenaphthene at 1.3e5 counts reads 0.006.
3. **"M as a share of ion current" was the wrong dominance test.** It was replaced by M relative to the base peak.
4. **Nitrogen-rule parity separates classes that share ions**: indole (M 117, odd) from indane (M 118, even), and
   an alkene fragment (odd) from a cycloalkane M+ (even).

## Known failures and blind spots

- Saturated peaks distort isotope ratios and can shift the molecular-ion candidate by a mass (chrysene-d12 blob
  3196 picks 241, parity flips, blob left unclassified, not misclassified).
- Primary alcohols look like alkenes at unit mass without m/z 31 (below the range). The 1-octanol reference comes
  out as alkene. Polarity (`dt2`) is the only handle and is a "prefer" term.
- Ketones and aliphatic amines share m/z 58, 72, 86. Only m/z 44 and an odd M separate them, so both may score
  and the tier is lowered.
- Alkenes and cycloalkanes are separated only by an even-mass M+ that is often absent.
- Styrene-type blobs (118, 117, 103) and methylnaphthols are not covered by any class.
- **Unresolved:** 312 of the 428 alkylbenzene blobs (51% of that class's volume) lie at t2 above 4 s and t1 7 to
  19 min, a region holding 36% of the run's volume. Their spectra are alkylbenzene-like, so the class is right
  about the chemistry, but whether this block is real wrapped-around aromatics or solvent-tail streaks cannot be
  told from this run. It makes the alkylbenzene t2-offset box in the figure very wide (median 6.1 s); the
  remaining 116 blobs sit at 1.3 s. The class volumes in the table below are dominated by it.

## Class volumes (all blobs; t1 of 11 min or later in brackets)

alkylbenzene 63.5% (61.8), unclassified 10.3% (11.0), indane, indene and tetralin 10.5% (11.5), aromatic oxygenate
5.0% (5.3), alkene 2.4% (2.3), PAH 2.3% (2.5), alkane 2.2% (1.4), naphthalene 2.1% (2.3), siloxane 1.3% (1.4),
N-heteroaromatic 0.2%, alcohol or ether 0.2%, cycloalkane 0.1%, ketone 0.0%, amine 0.0%. Tiers: 680 confident, 860
probable, 476 possible, 2052 unclassified (mostly weak or bleed-region blobs; 10.3% of the volume).

Nothing in the rules was tuned on the volumes above. Thresholds were set from the alkane ladder, from literature
EI patterns, and by inspecting individual near-miss blobs with `ClassResult.explain`.

## The 14 eV channel (same acquisition)

The 14 eV `.lsc` shares its `acquisition_id` with the 70 eV file and is read on the 70 eV pixel grid (`Fold` sharing).
Soft spectra corroborate and fix M parity, at a 3% M+ threshold.

- The soft M+ recovers 9 of 10 n-alkane molecular masses from the ladder.
- 86% of alkylbenzene calls have M = 8 (mod 14), as an alkylbenzene requires.
- The sulfur-heteroaromatic candidates were not supported by the soft spectra and stay rejected (class disabled).
- About 25 of 4068 calls changed with the soft channel on. It confirms rather than reclassifies.
- It cannot separate alkenes from cycloalkanes (both M = 0 mod 14).
