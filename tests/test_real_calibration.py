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


@pytest.mark.real_data
def test_peak_preview_on_native_calibration_data_preserves_inputs(tmp_path):
    value = os.environ.get("PYGCGC_CALIBRATION_DATA")
    if not value:
        pytest.skip("Set PYGCGC_CALIBRATION_DATA to the local calibration folder")
    pytest.importorskip('PySide6')
    import numpy as np
    from pygcgc.core import fold, read_run, remove_background
    from pygcgc.core.io import crop_run
    from pygcgc.core.components import ComponentDetectionSettings
    from pygcgc.gui.jobs import CancelToken
    from pygcgc.gui.peak_preview import preview_fit
    folder = Path(value)
    saved = (folder/'project.pygcgc.json').read_bytes()
    project = json.loads(saved)
    entry = project['runs'][0]
    path = folder/Path(entry['path']).name
    raw = read_run(path)
    original = fold(raw, project['period'], entry['phase'])
    ref = next(iter(project['istd']['references'].values()))
    crop = crop_run(raw, ref['t1']-60, ref['t1']+60)
    f = fold(crop, project['period'], entry['phase'], t0=original.spec.t0)
    cleaned = remove_background(crop, f).run
    stages = []
    result = preview_fit(cleaned, raw, f, ComponentDetectionSettings(workers=1), path.stem, None,
                         CancelToken(), lambda done, total, message: stages.append(message))
    assert result['count'] > 0 and result['tic'].shape == result['model'].shape
    assert np.isfinite(result['model']).all()
    assert any('Merging' in stage for stage in stages)
    assert any('full spectra' in stage for stage in stages)
    assert (folder/'project.pygcgc.json').read_bytes() == saved
    result.pop('tic')
    result.pop('model')
    result.pop('outlines')
    (tmp_path/'preview_summary.json').write_text(json.dumps(result, default=str, indent=2))
