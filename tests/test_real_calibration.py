"""Opt-in checks against local, private instrument data.

Set PYGCGC_CALIBRATION_DATA to the calibration folder to run this test. Nothing
is written beside the instrument files. A saved manual label is a reference for
recovery and cross-level consistency, not proof of chemical identity or amounts.
"""

import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.real_data
def test_default_detector_recovers_the_labelled_standard_at_each_level(tmp_path):
    value = os.environ.get("PYGCGC_CALIBRATION_DATA")
    if not value:
        pytest.skip("Set PYGCGC_CALIBRATION_DATA to the local calibration folder")
    folder = Path(value)
    project = folder / "project.pygcgc.json"
    assert project.exists(), "Calibration folder needs its saved manually labelled reference"
    output = tmp_path / "recovery.json"
    script = Path(__file__).resolve().parents[1] / "benchmarks" / "validate_calibration_peaks.py"
    subprocess.run([sys.executable, str(script), str(folder), "--out", str(output)],
                   check=True, capture_output=True, text=True, timeout=1800)
    report = json.loads(output.read_text())
    rows = report["results"]
    assert len({row["run"] for row in rows}) == 4
    for row in rows:
        assert row["ms-aware"]["found"], row
        assert row["ms-aware"]["spectral_score"] >= .8, row
        if row["ms-aware"]["peak_unresolved"]:
            assert row["ms-aware"]["needs_review"]
        assert row["component_count"] > 0
        assert row["nearby_components"]
