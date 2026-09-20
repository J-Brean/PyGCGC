# Real data: what these files are, and how to get them into PyGCGC

> Superseded by [NATIVE_LSC_REPORT.md](NATIVE_LSC_REPORT.md), 2026-09-19.
> Further investigation decoded all spectra in all three LSC files. The report
> corrects the SD04 offsets, timestamp interpretation, DAT/DAX sizes and processing
> metadata assumptions below, and proposes native import. This original note is
> retained as investigation history, not the current format specification.

Written after interrogating the six files by hand (hex dumps, struct decoding, cross-checks
against declared header values). No vendor documentation was used because none is public;
everything below is inferred and cross-validated internally, and is marked accordingly.

## 1. What arrived is not what the reader expects

`core/io/andi.py` reads ChromSpace's ANDI/AIA `.cdf` export. These six files are not that.
They are the **raw acquisition output of Markes TOF-DS** (the instrument-control software
that sits underneath ChromSpace), from a BenchTOF unit-mass ToF (serial `TOF-E10140`)
coupled to an Agilent 7890 GC, run in **tandem mode**: every other scanset is acquired at
70 eV and the one in between at 14 eV (`multimode=2`, `mode caption 0=_70eV`,
`mode caption 1=_14eV` in the `.HDR`). HANDOVER.md's assumption ("ChromSpace exports
ANDI") assumed a step in the pipeline (ChromSpace processing) that this data hasn't been
through yet.

| File | Role |
|---|---|
| `.HDR` | Plain text (Windows INI). Full method, GC oven/inlet/valve programme, mass calibration constants, run start/end time, total scanset count. Trivial to parse, and useful regardless of which route below is taken. |
| `.DAT` | Raw combined (both energies interleaved) detector stream for the whole run. Binary, compressed/thresholded, not scan-per-scan waveforms. |
| `.DAX` | A coarse seek index into `.DAT` (one entry per 64 scansets). |
| `..._70eV.lsc`, `..._14eV.lsc` | The two tandem channels, already demultiplexed into separate per-energy files. Self-contained, self-indexed binary containers (see below). |
| `..._70eV_0.1_secs.dbc.lsc` | The 70 eV channel *after* background subtraction. Same container format, same scan grid, much smaller. |

## 2. `.DAT` / `.DAX`: confirmed structure

`.DAT` opens with a 64-byte fixed header. Fields cross-checked against the `.HDR` text:

- bytes 0-3 (`u32`) = **9131** = `meas id` in the `.HDR`, exact match.
- bytes 20-23 (`u32`) = **13215** = `samples per scan` in the `.HDR`, exact match.
- bytes 44-47 (`u32`) = **64** = the `.DAX` block size (below), exact match.

`.DAX` starts with an 8-byte magic `"DAX1    "`, then 12-byte records `[u32 cumulative
scanset count][u64 byte offset in .DAT]`. The count climbs in steps of 64, and the file
holds 6720 records, i.e. 6720 x 64 = 430,080, against the `.HDR`'s `scansets=430034` -
consistent once you allow the last block to be partial. Each block of 64 scansets occupies
roughly 2.5 KB in `.DAT` (not the ~52 KB/scan a raw waveform would need), which fits
`meas data compression=1` in the `.HDR`: this is thresholded/ion-counted data, not raw
digitised waveforms.

This confirms the shape of the format but not enough to extract a usable chromatogram from
`.DAT` alone. It wasn't pursued further because the `.lsc` files (below) are already the
per-channel product and are far more tractable.

## 3. `.lsc`: a self-indexed, self-describing chunk container

All three `.lsc` files share one container format: a flat sequence of chunks, each

```
[4-byte ASCII tag][4-byte little-endian int32 length][length - 8 bytes of payload]
```

where `length` includes its own 8-byte header. Walking the whole of each file this way
parses **cleanly to the last byte, with zero framing errors**, across all three files -
strong evidence the container model is right, not a coincidence.

Chunks found, in order:

| Tag | Count | What it holds |
|---|---|---|
| `FT04` | 1 | 8-byte payload. An int64 that equals (file size - 1,720,216) in every file; turned out to be the byte offset of the `DI02` index (see below) once the whole thing was mapped out. |
| `PD01` | 1 | The method-parameter section of the `.HDR` text (`[TOF-DS Method]`, `[Markes TD]`, `[Agilent 7890 GC]`), verbatim. |
| `SM01` | 1 | The `[Start of measurement]` section of the `.HDR` text, verbatim. |
| `SD04` | 215,015 (in each of the three files) | **One per scan of that channel.** See section 4. |
| `MD01` | 1 | The *entire* `.HDR` file, verbatim, a second time, near the end. |
| `DI02` | 1 | A directory: seek offsets for the `PD01` and `SM01` chunks, then `[u32 SD04 count][one u64 byte-offset per SD04 chunk]`. Checked against the real file offsets of the first several `SD04` chunks - exact match. |

`215,015 (70 eV) + 215,015 (14 eV) = 430,030`, against the `.HDR`'s `scansets=430034` -
matches to within 4, which is what you'd expect from tandem alternation with a couple of
scansets dropped at the very start or end. This is good confirmation that the two `.lsc`
files really are "every other scanset" of the same combined acquisition seen in `.DAT`.

The `DI02` directory means **a `.lsc` file is self-indexed**: you can seek straight to scan
*N* in O(1) without scanning the whole file, and without needing `.DAT`/`.DAX` at all. That
makes `.lsc`, not `.DAT`, the sensible target if a custom reader is ever built.

## 4. `SD04`: one scan, part solved

Every `SD04` payload starts with a fixed 20-byte record:

- bytes 0-7: an **int64 that climbs steadily** through the file, ~7.67 x 10^18 at the
  start. Dividing the total span of this counter across the whole run
  (`18,469,560,000,000` ticks) by the run duration taken from the `.HDR`'s start/end
  timestamps (4300.4 s) gives **~4.29 x 10^9 ticks/second** - a physically sensible clock
  rate for a TOF digitiser. This is very likely the per-scan acquisition time, in raw
  hardware clock ticks, and is internally consistent enough to trust for **retention-time
  reconstruction** once the tick rate is pinned down precisely (it should be derivable
  exactly from the `.HDR`'s own clock-tick/timestamp pair, rather than the rough estimate
  used here).
- bytes 8-15: near-constant across nearby scans (e.g. `1786629098`), and matches the
  `.HDR`'s own `systemtime ticks=1786629098` field almost exactly. A coarse,
  low-resolution wall-clock snapshot, not needed for reconstruction.
- bytes 16-19: varies in a narrow band (4128-4193 in the sample checked), separate from
  the near-identical-looking `4137` seen once in the `.DAT` header. Looks like a
  housekeeping/analog readback, not chromatographic data.
- bytes 20-27 (not fully explored): zero in almost every scan checked, non-zero (a flag?)
  in at least the very last scan of the run.

After that 20-byte header, `SD04` payloads with no signal are exactly 28 bytes long
(header only). Payloads with signal run up to ~4.2 KB and hold the scan's **stick mass
spectrum** - but the exact per-ion encoding (how mass and intensity are packed) was **not**
pinned down. A first attempt (paired float32/uint16 values) produced numbers that don't
look like sensible masses or intensities, so it was not guessed further. Getting this
right needs either vendor documentation or a validated round-trip against a spectrum
exported from TOF-DS's own viewer for a known scan - not something to reverse-engineer
blind on a dataset headed for publication.

## 5. What the background-subtracted file tells us

`..._70eV_0.1_secs.dbc.lsc` has **exactly the same `SD04` count (215,015) and the same
per-scan timestamps** as the plain `..._70eV.lsc` - the first three raw timestamps in both
files are byte-identical. The only difference is size: 60 MB against 446 MB, a ~7.6x
reduction, entirely from smaller `SD04` payloads (less stick data per scan after
background subtraction), not from fewer or re-binned scans.

That means **"0.1_secs" in the filename is very unlikely to be a resampled time axis** -
it looks like a parameter of whatever background-correction step was applied (an
integration or averaging window), not a change to the scan grid. Worth confirming with
whoever ran the background subtraction rather than assuming.

## 6. Recommendation

**Ask first: does anything downstream of TOF-DS (ChromSpace, or another Markes/SepSolve
tool) sit between this raw acquisition and analysis, and can it export ANDI/AIA `.cdf` or
mzML from these files?** If yes, that is by far the safest and cheapest route in - it needs
no new reader, matches the existing `core/io/andi.py` exactly, and carries none of the risk
below. This mirrors HANDOVER.md's original conclusion about `.lsc`, and nothing found here
changes that recommendation - it just explains *why* `.lsc` looked opaque (it's a real,
consistent, but undocumented format, not noise).

If no such export exists (raw TOF-DS output is genuinely all there is), a custom reader is
feasible but not trivial, and should not be built past the point of "we understand the
container" until the `SD04` stick-data encoding is validated against a known-correct
spectrum. Concretely, before writing `core/io/markes_lsc.py`:

1. Get one scan's mass spectrum out of TOF-DS's own viewer (or ChromSpace) for a specific,
   known retention time, and use it as ground truth.
2. Decode that scan's raw `SD04` bytes and confirm which sub-fields correspond to mass and
   intensity, and the exact calibration formula (the `.HDR`'s `Offset T0` / `Mass per Time`
   / `Mass per Time2` fields are very likely part of this, but the functional form -
   probably time as a function of sqrt(mass), the usual TOF relation, rather than a plain
   quadratic in mass - needs confirming, not assuming).
3. Pin down the exact tick-to-second conversion for the bytes 0-7 timestamp (should come
   directly from the `.HDR`'s clock-tick and wall-clock fields, not the rough estimate
   used above).
4. Only then write a reader, using `DI02` for O(1) scan access and reading `.HDR` directly
   for method metadata (steps 5 onwards of HANDOVER.md section 4 - target table, phase,
   background - are unaffected either way, since they operate on whatever `Run` object
   the reader produces).

One thing that's safe to do regardless of which route wins: the `.HDR` file is plain text
and trivial to parse. Even before any spectrum data is readable, it hands over the mass
calibration constants, full GC/valve method, tandem-mode configuration, and exact run
start/end time and scanset count for provenance and cross-checking - worth a small parser
on its own.
