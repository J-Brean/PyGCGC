# real_data

Local instrument files. Everything here except this note is git-ignored, and nothing in it
should be committed without asking James.

Expected, for the tests and benchmarks that use real data:

- the native Markes files (`*_70eV.lsc`, `*_14eV.lsc`, optionally `*.dbc.lsc`, with their
  `.DAT`, `.DAX` and `.HDR` companions);
- `blob_diagnosis/segment.npz`, the cached ten-minute crop of the earlier (10 s modulation)
  acquisition that the MS-aware benchmarks in `benchmarks/` read.

Generated here by the benchmarks: `.cache/` (decoded runs), `quant_validation/`,
`target_table_istd.csv`, `blob_diagnosis/*.json`. Tests that need real files look here or at
`PYGCGC_REAL_DATA`, and skip when nothing is found.
