"""🧬 Compound-class step in the GUI (pytest-qt, headless)."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
pytest.importorskip("pyqtgraph")

from pygcgc.gui.blob_state import MS_AWARE  # noqa: E402
from tests.test_gui import TIMEOUT, window  # noqa: E402,F401  (the fixture)


def _blobs(qtbot, w):
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    return i


def _classified(qtbot, w, i):
    w.start_classify()
    qtbot.waitUntil(lambda: w.classes.result_for(i) is not None and not w.classes.running, timeout=TIMEOUT)
    return w.classes.result_for(i)


def test_classify_is_step_seven_between_hydrocarbons_and_identify(qtbot, window):
    keys = window.workflow.keys()
    assert keys.index("hydro") + 1 == keys.index("classify") == keys.index("identify") - 1
    assert window.step_actions["classify"].text().lstrip("✓▸ ") == "7 🧬 Classify"
    assert window.step_actions["identify"].text().lstrip("✓▸ ") == "8 🏷️ Identify"
    assert window.step_actions["quant"].text().lstrip("✓▸ ") == "9 🧪 Quantify"
    assert window.workflow.step_info("classify")["why"] in window.step_actions["classify"].toolTip()


def test_go_classifies_legacy_blobs_and_shows_why(qtbot, window):
    w = window
    i = _blobs(qtbot, w)
    blobs = w.blobs.blobs_for(i)
    w.step_actions["classify"].trigger()
    w.workflow.goRequested.emit("classify")
    qtbot.waitUntil(lambda: w.classes.result_for(i) is not None, timeout=TIMEOUT)
    result = w.classes.result_for(i)
    assert len(result) == len(blobs) and w.classes.blobs_classified(i) is blobs
    assert w.class_panel.table.rowCount() == len(blobs)
    assert "from 🔍" in w.class_panel.notes.text() and "no 14 eV run loaded" in w.class_panel.notes.text()
    assert w.workflow.is_done("classify")
    w.class_panel.table.setCurrentCell(0, 0)
    assert w.class_panel.explain.toPlainText().startswith("blob ")
    assert w.class_panel.summary.rowCount() >= 1


def test_ms_aware_users_get_classes_without_running_the_legacy_detector(qtbot, window):
    w = window
    w.blobs.set_method(MS_AWARE)
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    assert w.blobs.blobs_for(i) is None
    result = _classified(qtbot, w, i)
    assert len(result) > 0 and w.classes.result_for(i) is result
    assert w.project.entries[i] in w.classes.own_blobs and "found here" in w.class_panel.notes.text()


def test_the_14_eV_run_is_refused_and_found_as_a_partner(qtbot, window):
    w = window
    a, b = w.project.index_of("RUN_A"), w.project.index_of("RUN_B")
    w.project.entries[a].import_meta = {"acquisition_id": "x", "energy_ev": 70.0, "processing_state": "instrument"}
    w.project.entries[b].import_meta = {"acquisition_id": "x", "energy_ev": 14.0, "processing_state": "instrument"}
    assert w.classes.soft_partner(a) == b and w.classes.soft_partner(b) is None
    assert "14 eV" in w.classes.why_not(b) and w.classes.why_not(a) == ""
    w.set_run(b)
    assert not w.workflow._by_key["classify"]["go_enabled"] and "70 eV" in w.class_panel.status.text()
    assert not w.classes.classify(b)
    w.set_run(a)
    assert "RUN_B" in w.class_panel.soft_status.text()


def test_14_eV_partner_is_read_through_the_job_and_feeds_the_soft_features(qtbot, window):
    w = window
    a, b = w.project.index_of("RUN_A"), w.project.index_of("RUN_B")
    w.project.entries[a].import_meta = {"acquisition_id": "x", "energy_ev": 70.0, "processing_state": "instrument"}
    w.project.entries[b].import_meta = {"acquisition_id": "x", "energy_ev": 14.0, "processing_state": "instrument"}
    from pygcgc.core import ClassSettings
    w.classes.set_settings(ClassSettings(min_ions=1, m_min_mass=0))   # synthetic peaks: one or two ions, m/z 50-54
    # RUN_B is a time-shifted copy (for the alignment tests), so its pixels do not sit under RUN_A's blobs.
    # A stand-in 14 eV channel that shares RUN_A's data exercises the plumbing instead.
    w.project.entries[b].run = w.project.entries[a].run
    _blobs(qtbot, w)
    result = _classified(qtbot, w, a)
    assert "14 eV from RUN_B" in w.classes.notes[w.project.entries[a]]
    assert np.isfinite(result.features["soft_M_mass"]).any()
    w.classes.set_use_soft(False)
    assert w.classes.is_stale(a)
    result_off = _classified_again(qtbot, w, a, result)
    assert np.isnan(result_off.features["soft_M_mass"]).all()


def _classified_again(qtbot, w, i, previous):
    w.start_classify()
    qtbot.waitUntil(lambda: w.classes.result_for(i) is not previous and not w.classes.running, timeout=TIMEOUT)
    return w.classes.result_for(i)


def test_results_go_stale_when_rules_or_blobs_change_and_settings_persist(qtbot, window, tmp_path):
    w = window
    i = _blobs(qtbot, w)
    _classified(qtbot, w, i)
    assert not w.classes.is_stale(i)
    w.classes.set_enabled("alkane", False)
    assert w.classes.is_stale(i)
    assert "alkane" not in [c.name for c in w.classes.rules_in_use().classes if c.enabled]
    saved = w.project.to_dict()["classes"]
    assert saved["enabled"] == {"alkane": False}
    w.classes.set_enabled("alkane", True)
    w.classes.from_dict(saved)
    assert w.classes.enabled == {"alkane": False} and not w.classes.results


def test_a_bad_rule_file_is_refused_with_the_reason_and_the_old_rules_kept(qtbot, window, tmp_path):
    w = window
    before = w.classes.rules
    bad = tmp_path / "bad.toml"
    bad.write_text('[ions]\nx = [43]\n[[class]]\nname = "c"\nneed = ["nope: 0.1 -> 0.2"]\n', encoding="utf-8")
    with pytest.raises(ValueError, match="unknown feature 'nope'"):
        w.classes.set_rules_path(bad)
    assert w.classes.rules is before


def test_exported_csv_has_one_row_per_blob(qtbot, window, tmp_path):
    from pygcgc.core import write_classes_csv
    w = window
    i = _blobs(qtbot, w)
    result = _classified(qtbot, w, i)
    path = write_classes_csv(result, w.classes.blobs_classified(i), tmp_path / "c.csv")
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("blob,t1_min") and len(lines) == len(result) + 1
