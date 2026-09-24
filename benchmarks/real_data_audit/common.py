"""Shared set-up for the September 2026 real-data audit scripts.

The four tandem (70 and 14 eV) acquisitions in real_data/ take about 15 s each
to decode, so the decoded runs are cached under real_data/.cache (git-ignored)
the first time a script needs them.
"""

from pathlib import Path

import numpy as np
import scipy.sparse as sp

from pygcgc import Run, read_run

ROOT = Path(__file__).resolve().parents[2]
REAL = ROOT / "real_data"
CACHE = REAL / ".cache"

# The fold the user set in the GUI project for these runs (real_data/project.pygcgc.json).
PERIOD, PHASE = 8.0, 2.77

RUNS = {"ATWC": "2025-10-11_964969_REPSOL_2_ATWC_250807_1",
        "BTWC": "2025-10-13_964767_REPSOL_2_BTWC_250807_1",
        "T191": "2025-10-13_452283_D50T191_250703_1",
        "T60": "2025-10-13_452285_D50T60_250703_1"}


def load(key: str, energy: str = "70eV") -> Run:
    """One channel of one acquisition, from the cache when it is there."""
    stem = f"{RUNS[key]}_{energy}"
    cached = CACHE / f"{stem}.npz"
    if cached.exists():
        z = np.load(cached)
        data = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
        return Run(times=z["times"], mz=z["mz"], data=data, name=f"{key}_{energy}")
    run = read_run(REAL / f"{stem}.lsc")
    CACHE.mkdir(exist_ok=True)
    d = run.data
    np.savez(cached, times=run.times, mz=run.mz, data=d.data, indices=d.indices, indptr=d.indptr,
             shape=np.array(d.shape))
    run.name = f"{key}_{energy}"
    return run
