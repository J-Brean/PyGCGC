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
from pygcgc.gui.blob_state import LEGACY, MS_AWARE  # noqa: E402

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
    w.auto_process = False          # phases are set by hand below, before processing
    w.set_peak_method(LEGACY)       # most tests use legacy blobs: footprints, clicks that add blobs
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
    assert all(w.peaks.get(e) is not None for e in w.project.entries)


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
    assert w.normalised and w.fits                        # divided as soon as there are blobs and an IS
    w.normalise()
    assert w.tabs.currentIndex() == 2
    assert w.tabs.tabText(2).startswith("✓")
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
    w.tabs.setCurrentIndex(1)
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


def test_whole_workflow_from_the_main_window_to_sample_amounts(qtbot, cal_folder, tmp_path):
    """Select the calibration runs in the main window, open 📈, give level amounts, let the
    compounds be suggested, then quantify the main window's samples."""
    import csv
    from PySide6 import QtCore
    from pygcgc.core.calbuilder import SuggestSettings
    from pygcgc.gui.calibrate_window import SuggestDialog
    from pygcgc.gui.mainwindow import MainWindow
    from tests.conftest import SAMPLE_AMOUNTS
    mw = MainWindow()
    qtbot.addWidget(mw)
    mw.show()
    loaded = []
    mw.open_folder(cal_folder, on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and mw.jobs.active_count() == 0, timeout=TIMEOUT)
    qtbot.waitUntil(lambda: all(e.phase_method != "unset" for e in mw.project.entries), timeout=TIMEOUT)
    # 1: highlight the calibration runs in the main window's Runs table
    table = mw.runs_panel.table
    table.clearSelection()
    mode = table.selectionMode()
    table.setSelectionMode(table.SelectionMode.MultiSelection)
    for i, e in enumerate(mw.project.entries):
        if e.name.startswith("CAL_"):
            table.selectRow(i)
    table.setSelectionMode(mode)
    mw.act_calibrate.trigger()
    cw = mw.calibrate_window
    assert len(cw.main_selection()) == len(RUNS)
    assert "5 runs selected" in cw.from_main_btn.text()
    cw.runs_from_main()
    assert cw.peak_method == MS_AWARE                        # the main window's 🔍 default
    # While processing, tab 1 must not send anyone to tab 2, and tab 2 must say it is coming.
    assert cw.processing and "Processing" in cw.standards_banner.text()
    assert "Go to" not in cw.standards_banner.text() and "being processed" in cw.compounds_banner.text()
    qtbot.waitUntil(lambda: bool(cw.peaks) and not cw.processing and cw.jobs.active_count() == 0,
                    timeout=TIMEOUT * 15)
    assert all(type(b).__name__ == "ComponentBlobs" for b in cw.peaks.values())
    assert cw.session.runs == RUNS
    assert all(e.run is mw.project.entries[mw.project.index_of(e.name)].run for e in cw.project.entries)
    assert cw.tabs.currentIndex() == 1                       # moved on to the compounds by itself
    for r, v in zip(RUNS, CAL_LEVELS.values()):              # names say no amount here: type the levels
        cw.set_level(r, v)
    # 2: compounds suggested, the constant one kept as the internal standard
    dlg = SuggestDialog(cw.session.runs, cw.session.levels, cw.blobs_by_run(), SuggestSettings(), cw,
                        cw.usable_masks())
    qtbot.addWidget(dlg)
    dlg.find()
    for k, g in enumerate(dlg.found):
        if g.kind == "constant":
            dlg.table.item(k, 0).setCheckState(QtCore.Qt.CheckState.Checked)
    cw.add_suggestions(dlg.chosen(), dlg.amounts())
    # 3: curves, with nothing to press
    assert sorted(g.kind for g in dlg.found) == ["calibrant"] * 3 + ["constant"], [
        (g.kind, round(g.t1), round(g.slope, 2), len(g.blobs)) for g in dlg.found]
    assert len(cw.fits) == 3 and all(isinstance(f, CalibrationFit) and f.ok for f in cw.fits.values())
    # 4: the main window's samples, ticked by default (the calibration runs are not offered)
    cw.tabs.setCurrentIndex(3)
    assert set(cw.apply.ticked()) == set(SAMPLE_AMOUNTS)
    done = []
    cw.quantify_samples(cw.apply.ticked(), on_done=lambda: done.append(True))
    qtbot.waitUntil(lambda: bool(done), timeout=TIMEOUT * 2)
    assert cw.tabs.tabText(3).startswith("✓")
    truth = {}
    for g in cw.session.analytes():          # which suggested compound is A, B or C: by position
        name = min((("A", 110.0), ("B", 170.0), ("C", 200.0)), key=lambda p: abs(p[1] - g.t1))[0]
        truth[g.name] = name
    for r in cw.sample_results:
        assert r.found and not r.flags, (r.run, r.compound, r.flags)
        assert r.amount == pytest.approx(SAMPLE_AMOUNTS[r.run][truth[r.compound]], rel=0.08)
    path = cw.export_results(tmp_path / "results.csv")
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    assert len(rows) == 6 and all(float(x["concentration"]) == pytest.approx(float(x["amount"]) / 2.0)
                                  for x in rows)
    # A curve change updates the sample results by itself.
    before = {(r.run, r.compound): r.amount for r in cw.sample_results}
    name = cw.session.analytes()[0].name
    cw._toggle_point(name, RUNS[-1])
    after = {(r.run, r.compound): r.amount for r in cw.sample_results}
    assert after[("S1", name)] != before[("S1", name)]
    cw.jobs.cancel_all()
    cw.jobs.wait(5000)
    mw.jobs.cancel_all()
    mw.jobs.wait(5000)
    mw.undo.clear()


def test_the_main_windows_ms_aware_peaks_are_used_as_they_are(qtbot, cal_folder, monkeypatch):
    """Peaks found (and perhaps separated by hand) in the main window's 🔍 step are the ones the
    builder works on: it must not find its own."""
    import pygcgc.gui.calibrate_window as cwmod
    from pygcgc.gui.mainwindow import MainWindow
    mw = MainWindow()
    qtbot.addWidget(mw)
    loaded = []
    mw.open_folder(cal_folder, on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and mw.jobs.active_count() == 0, timeout=TIMEOUT)
    qtbot.waitUntil(lambda: all(e.phase_method != "unset" for e in mw.project.entries), timeout=TIMEOUT)
    mw.blobs.set_method(MS_AWARE)
    mw.run_blobs()
    qtbot.waitUntil(lambda: not mw.blobs.running and len(mw.blobs.component_sets_by_name) == 7,
                    timeout=TIMEOUT * 30)

    def refuse(*_a, **_k):
        raise AssertionError("the builder searched for peaks the main window had already found")
    monkeypatch.setattr(cwmod, "detect_components", refuse)
    mw.show_calibrate_window()
    cw = mw.calibrate_window
    cw.act_line_up.setChecked(False)          # a moved phase would rightly need new peaks
    cw.add_main_entries([e for e in mw.project.entries if e.name.startswith("CAL_")])
    qtbot.waitUntil(lambda: len(cw.peaks) == len(RUNS) and not cw.processing, timeout=TIMEOUT * 3)
    for e in cw.project.entries:
        n_main = sum(len(g.components) for g in mw.blobs.component_sets_by_name[e.name].groups)
        assert len(cw.peaks[e]) == n_main
    cw.jobs.cancel_all()
    cw.jobs.wait(5000)
    mw.jobs.cancel_all()
    mw.jobs.wait(5000)
    mw.undo.clear()


def test_phases_are_left_alone_when_most_runs_disagree(qtbot, cal_folder):
    # On the real standards a bad set of peaks once moved every run's phase. Lining up exists to
    # bring a stray run back to the others; with no majority to bring it back to, nothing moves.
    w = _open(qtbot, cal_folder)
    w.project.set_phases({i: (GUI_PHASE + 0.3 * i, "manual", "") for i in range(len(w.project.entries))})
    _process(qtbot, w)
    before = [e.phase for e in w.project.entries]
    assert not w.line_up_phases(min_pairs=3)
    assert [e.phase for e in w.project.entries] == before
    assert "disagree" in w.statusBar().currentMessage()
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
    cw.auto_process = False                  # only the fold matters here
    loaded = []
    cw.open_folder(folder, on_loaded=lambda: loaded.append(True))
    qtbot.waitUntil(lambda: bool(loaded) and cw.jobs.active_count() == 0, timeout=TIMEOUT)
    assert cw.session.runs == RUNS
    assert cw.project.period == pytest.approx(4.0, abs=0.01)
    assert cw.period.value() == pytest.approx(cw.project.period, abs=1e-4)     # the box shows it
    cw.jobs.cancel_all()
    cw.jobs.wait(5000)
    cw.close()
