"""📈 Calibration builder window, driven the way a person would: open the standard
runs, load or click compounds, clean and find blobs, normalise, then edit curves."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
pytest.importorskip("pyqtgraph")

from pygcgc.core.calibration import CalibrationFit  # noqa: E402
from tests.conftest import CAL_LEVELS, GUI_PHASE, IS_AREA, SCALE  # noqa: E402
from tests.test_gui import TIMEOUT, window  # noqa: E402,F401

RUNS = [f"CAL_{lv}" for lv in CAL_LEVELS]


def _compound_csv(path):
    # Amount columns named by level, which the table matches to CAL_L1 ... CAL_L5.
    head = "compound,t1_s,t2_s,role,amount," + ",".join(CAL_LEVELS)
    rows = ["IS,60,1.2,istd,1,,,,,"]
    for name, t1, t2 in (("A", 110, 2.0), ("B", 170, 2.6), ("C", 200, 1.5)):
        rows.append(f"{name},{t1},{t2},target,," + ",".join(f"{v:g}" for v in CAL_LEVELS.values()))
    path.write_text("\n".join([head] + rows) + "\n", encoding="utf-8")
    return path


def _open(qtbot, cal_folder):
    from pygcgc.gui.calibrate_window import CalibrationWindow
    w = CalibrationWindow()
    qtbot.addWidget(w)
    w.show()
    loaded = []
    w.add_runs([cal_folder / "data" / f"{r}.cdf" for r in RUNS], on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and w.jobs.active_count() == 0, timeout=TIMEOUT)
    assert w.project.period == pytest.approx(4.0)          # from the method.json beside the runs
    w.project.set_phases({i: (GUI_PHASE, "manual", "") for i in range(len(w.project.entries))})
    return w


def _process(qtbot, w):
    done = []
    w.run_pipeline(on_ready=lambda: done.append(True))
    qtbot.waitUntil(lambda: bool(done), timeout=TIMEOUT)
    assert all(w.project.blobs.by_run.get(e) is not None for e in w.project.entries)


@pytest.fixture
def cal_window(qtbot, cal_folder, tmp_path):
    w = _open(qtbot, cal_folder)
    w.load_compounds(_compound_csv(tmp_path / "compounds.csv"))
    _process(qtbot, w)
    yield w
    w.jobs.cancel_all()
    w.jobs.wait(5000)


def test_table_to_curves(cal_window):
    w = cal_window
    assert w.session.runs == RUNS
    assert w.session.compound("A").amounts == {f"CAL_{lv}": v for lv, v in CAL_LEVELS.items()}
    assert w.session.compound("IS").is_istd
    assert not w.fits                                     # nothing divided until asked
    w.normalise()
    assert w.tabs.currentWidget() is w.curves
    fit = w.fits["A"]
    assert isinstance(fit, CalibrationFit) and fit.ok
    assert fit.slope == pytest.approx(SCALE * 1.3 / (IS_AREA * 1.5), rel=0.05)
    assert w.curves.compounds.count() == 3
    assert w.curves.table.rowCount() == len(RUNS)


def test_leave_out_put_back_and_undo(cal_window):
    w = cal_window
    w.normalise()
    w._toggle_point("B", "CAL_L5")
    assert w.fits["B"].n == 4 and "CAL_L5" in w.fits["B"].excluded
    w._toggle_point("B", "CAL_L5")
    assert w.fits["B"].n == 5
    w._toggle_point("B", "CAL_L1")
    w.undo()
    assert w.fits["B"].n == 5
    w.redo()
    assert w.fits["B"].n == 4


def test_divide_by_chooses_an_internal_standard(qtbot, cal_window):
    from PySide6 import QtWidgets
    from pygcgc.core.calbuilder import CalCompound
    w = cal_window
    w._edit(lambda s: s.add_compound(CalCompound("IS2", 230.0, 3.5, is_istd=True,
                                                 amounts={r: 1.0 for r in RUNS})))
    w._edit(lambda s: None)                     # refresh twice: combos are reused, never piled up
    row = w.session.names().index("A")
    combo = w.table.cellWidget(row, 5)
    assert combo.itemText(0) == "nearest (IS)"
    combo.setCurrentIndex(combo.findData("IS2"))
    qtbot.waitUntil(lambda: w.session.compound("A").istd == "IS2", timeout=TIMEOUT)
    assert w.points[("A", "CAL_L1")].istd == "IS2"
    combos = [c for c in w.table.findChildren(QtWidgets.QComboBox) if c.isVisible()]
    assert len(combos) == len(w.session.compounds)
    w.undo()
    assert w.session.compound("A").istd == ""


def test_retyping_an_amount_moves_the_point(cal_window):
    from pygcgc.gui.calibrate_window import FIXED_COLUMNS
    w = cal_window
    w.normalise()
    w.edit_cell("C", len(FIXED_COLUMNS) + RUNS.index("CAL_L2"), "")      # blank: no point in that run
    assert "CAL_L2" not in w.session.compound("C").amounts
    assert w.fits["C"].n == 4
    w._set_amount("C", "CAL_L2", 1.0)
    assert w.fits["C"].n == 5


def test_click_rematches_a_compound_and_no_peak_removes_it(cal_window):
    w = cal_window
    w.normalise()
    i = w.project.index_of("CAL_L3")
    w.runs.setCurrentRow(i)
    blobs = w.project.blobs.by_run[w.project.entries[i]]
    c_blob = w.resolved[("C", "CAL_L3")]
    w.select_compound("A")
    w.click_mode.setCurrentIndex(0)                      # match the selected compound to this blob
    w._view_clicked(int(blobs.col[c_blob]), float(blobs.t2[c_blob]))
    assert w.resolved[("A", "CAL_L3")] == c_blob
    assert w.points[("A", "CAL_L3")].pinned
    w._match_current("none")
    assert w.points[("A", "CAL_L3")].blob is None and w.fits["A"].n == 4
    w._set_match("A", "CAL_L3", "auto")
    assert w.resolved[("A", "CAL_L3")] != c_blob and w.fits["A"].n == 5


def test_click_adds_compounds_and_internal_standards(qtbot, cal_folder, monkeypatch):
    from PySide6 import QtWidgets
    w = _open(qtbot, cal_folder)
    _process(qtbot, w)
    i = w.project.index_of("CAL_L4")
    w.runs.setCurrentRow(i)
    blobs = w.project.blobs.by_run[w.project.entries[i]]
    names = iter(["IS", "A"])
    monkeypatch.setattr(QtWidgets.QInputDialog, "getText", lambda *a, **k: (next(names), True))
    is_blob = min(range(len(blobs)), key=lambda k: abs(blobs.t1[k] - 60))
    a_blob = min(range(len(blobs)), key=lambda k: abs(blobs.t1[k] - 110))
    w.click_mode.setCurrentIndex(3)                      # add an internal standard here
    w._view_clicked(int(blobs.col[is_blob]), float(blobs.t2[is_blob]))
    w.click_mode.setCurrentIndex(2)                      # add a compound here
    w._view_clicked(int(blobs.col[a_blob]), float(blobs.t2[a_blob]))
    assert [c.name for c in w.session.compounds] == ["IS", "A"]
    assert w.session.compound("IS").is_istd and not w.session.compound("A").is_istd
    assert w.session.compound("A").t1 == pytest.approx(float(blobs.t1[a_blob]))
    # The new compounds are found in every run, not only the one clicked.
    assert all(w.resolved[("A", r)] is not None for r in RUNS)
    monkeypatch.setattr(QtWidgets.QInputDialog, "getDouble", lambda *a, **k: (1.0, True))
    w.select_compound("IS")
    w._fill_amount()
    assert w.session.compound("IS").amounts == {r: 1.0 for r in RUNS}
    for r, v in zip(RUNS, CAL_LEVELS.values()):
        w._set_amount("A", r, v)
    w.normalise()
    assert isinstance(w.fits["A"], CalibrationFit) and w.fits["A"].ok
    w.jobs.cancel_all()
    w.jobs.wait(5000)


def test_suggested_compounds_give_curves_without_a_table(qtbot, cal_folder):
    from PySide6 import QtCore
    from pygcgc.core.calbuilder import SuggestSettings
    from pygcgc.gui.calibrate_window import SuggestDialog
    w = _open(qtbot, cal_folder)
    _process(qtbot, w)
    dlg = SuggestDialog(w.session.runs, {}, w.blobs_by_run(), SuggestSettings(), w)
    qtbot.addWidget(dlg)
    for r, v in enumerate(CAL_LEVELS.values()):
        dlg.levels.item(r, 1).setText(f"{v:g}")
    dlg.find()
    kinds = sorted(g.kind for g in dlg.found)
    assert kinds.count("calibrant") == 3 and "constant" in kinds
    for k, g in enumerate(dlg.found):                     # keep the constant peak as the IS
        if g.kind == "constant":
            dlg.table.item(k, 0).setCheckState(QtCore.Qt.CheckState.Checked)
    w.add_suggestions(dlg.chosen(), dlg.amounts())
    assert len(w.session.analytes()) == 3 and len(w.session.istds()) >= 1
    assert all(p.blob is not None for p in w.points.values())
    w.normalise()
    assert all(isinstance(f, CalibrationFit) and f.ok for f in w.fits.values())
    w.jobs.cancel_all()
    w.jobs.wait(5000)


def test_export_and_session_round_trip(qtbot, cal_window, tmp_path):
    from pygcgc.gui.calibrate_window import CalibrationWindow
    w = cal_window
    w.normalise()
    i = w.project.index_of("CAL_L2")
    w.runs.setCurrentRow(i)
    w.select_compound("B")
    w._match_current("none")
    w._toggle_point("A", "CAL_L1")
    w._set_override("C", {"model": "linear", "weights": "1/x"})
    names = {p.name for p in w.export(tmp_path)}
    assert names == {"calibration_compounds.csv", "calibration_curves.csv", "calibration_points.csv",
                     "calibration_normalised_blobs.csv"}
    path = w.save_session(tmp_path / "session.json")

    w2 = CalibrationWindow()
    qtbot.addWidget(w2)
    ready = []
    w2.open_session(path, on_ready=lambda: ready.append(True))
    qtbot.waitUntil(lambda: bool(ready), timeout=TIMEOUT)
    assert w2.session.to_dict() == w.session.to_dict()
    assert [e.phase for e in w2.project.entries] == pytest.approx([e.phase for e in w.project.entries])
    assert w2.normalised
    assert w2.points[("B", "CAL_L2")].blob is None
    for name in ("A", "B", "C"):
        assert w2.fits[name].slope == pytest.approx(w.fits[name].slope)
    w2.jobs.cancel_all()
    w2.jobs.wait(5000)


def test_all_outlines_draws_the_same_edges_as_blob_by_blob(cal_window):
    import numpy as np
    from pygcgc.gui.calibrate_window import all_outlines
    blobs = next(iter(cal_window.project.blobs.by_run.values()))

    def segments(x, y):
        s = np.column_stack([np.asarray(x).reshape(-1, 3)[:, :2], np.asarray(y).reshape(-1, 3)[:, :2]])
        return {tuple(np.round(r, 6)) for r in s}

    one_by_one = set()
    for i in range(len(blobs)):
        one_by_one |= segments(*blobs.outline_times(i))
    assert segments(*all_outlines(blobs)) == one_by_one


def test_internal_standards_with_an_ion_anchor_the_phase(qtbot, cal_folder, tmp_path):
    w = _open(qtbot, cal_folder)
    w.project.set_phases({i: (0.0, "manual", "") for i in range(len(w.project.entries))})
    p = tmp_path / "is.csv"
    p.write_text("compound,t1_s,t2_s,role,quant_ion\nIS,60,1.2,istd,50\nA,110,2.0,target,52\n", encoding="utf-8")
    w.load_compounds(p)
    assert [a.name for a in w.project.anchors] == ["IS"]
    done = []
    w.project.estimate_phases(on_done=lambda *_: done.append(True))
    qtbot.waitUntil(lambda: bool(done), timeout=TIMEOUT)
    for e in w.project.entries:
        assert e.phase_method == "anchors" and e.phase == pytest.approx(GUI_PHASE, abs=0.05)
    w.jobs.cancel_all()
    w.jobs.wait(5000)


def test_main_window_opens_the_builder_with_its_period(qtbot, window):
    w = window
    w.act_calibrate.trigger()
    cw = w.calibrate_window
    assert cw.isVisible() and cw.project is not w.project
    assert cw.project.period == pytest.approx(w.project.period)
    assert cw.project.request_autosave() is None             # never writes the main project's file
    cw.close()


def test_a_run_at_the_wrong_phase_is_lined_up_with_the_others(qtbot, cal_folder):
    w = _open(qtbot, cal_folder)
    wrong = w.project.index_of("CAL_L3")
    w.project.set_phase(wrong, GUI_PHASE + 0.4)          # as an independent estimate can be off
    _process(qtbot, w)
    assert w.line_up_phases(min_pairs=3)
    assert w.project.entries[wrong].phase == pytest.approx(GUI_PHASE, abs=0.03)
    assert w.project.entries[wrong].phase_method == "lined up"
    assert all(e.phase == pytest.approx(GUI_PHASE) for i, e in enumerate(w.project.entries) if i != wrong)
    _process(qtbot, w)                                   # clean and find blobs again at the new phase
    assert not w.line_up_phases(min_pairs=3)             # nothing left to move
    w.jobs.cancel_all()
    w.jobs.wait(5000)


def test_standards_get_their_own_period_not_the_main_windows(qtbot, window, cal_folder, tmp_path):
    # The main window's runs are at one period (set by hand here); standards in a folder
    # without a method.json are folded at the period estimated from themselves.
    import shutil
    w = window
    w.project.set_period(6.0)
    folder = tmp_path / "standards"
    folder.mkdir()
    for r in RUNS:
        shutil.copy(cal_folder / "data" / f"{r}.cdf", folder)
    w.show_calibrate_window()
    cw = w.calibrate_window
    loaded = []
    cw.open_folder(folder, on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and cw.jobs.active_count() == 0, timeout=TIMEOUT)
    assert cw.session.runs == RUNS
    assert cw.project.period == pytest.approx(4.0, abs=0.01)
    assert cw.period.value() == pytest.approx(cw.project.period, abs=1e-4)     # the box shows it
    cw.jobs.cancel_all()
    cw.jobs.wait(5000)
    cw.close()
