"""The matcher's plot explorer, the RI legend, naming by double-click (library or not) and the chemical-space window."""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from PySide6 import QtWidgets  # noqa: E402

from pygcgc.core.chemspace import parse_formula, space_of  # noqa: E402
from pygcgc.gui.chemspace_window import ChemSpaceWindow, space_rows  # noqa: E402
from pygcgc.gui.match_explorer import PLOTS  # noqa: E402
from pygcgc.gui.peak_naming import PeakNameDialog, identity_of  # noqa: E402
from tests.test_gui import TIMEOUT, window  # noqa: E402,F401


def formula_library(path):
    from pygcgc.core.identify import LibraryEntry, write_msp
    entries = [LibraryEntry(name="Compound50", mz=np.array([50]), intensity=np.array([999.0]), formula="C10H22"),
               LibraryEntry(name="Compound51", mz=np.array([51]), intensity=np.array([999.0]), formula="C6H6O")]
    return write_msp(entries, path)


def ready(qtbot, w, tmp_path):
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running, timeout=TIMEOUT)
    w.set_run(w.project.index_of("RUN_A"))
    w.show_library_matcher()
    state = w.project.matcher
    state.load_library(formula_library(tmp_path / "lib.msp"))
    qtbot.waitUntil(lambda: state.library is not None, timeout=TIMEOUT)
    w.library_matcher.match_btn.click()
    qtbot.waitUntil(lambda: state.current_for(w.index) is not None, timeout=TIMEOUT)
    return w.library_matcher, state


def dialog_returning(monkeypatch, typed, accept=True):
    def run(self):
        self.edit.setText(typed)
        return 1 if accept else 0
    monkeypatch.setattr(PeakNameDialog, "exec", run)


# ---------------------------------------------------------------- formulae
def test_formulae_give_the_chemical_space_coordinates():
    decane = space_of("C10H22")
    assert (decane.hc, decane.oc, decane.dbe, decane.c) == (2.2, 0.0, 0.0, 10)
    assert space_of("C6H6").dbe == 4 and space_of("C6H6O").oc == pytest.approx(1 / 6)
    assert space_of("C14D10").hc == pytest.approx(10 / 14) and space_of("C14D10").dbe == 10     # deuterium is hydrogen
    assert space_of("C10H22").family == "CH" and space_of("C6H6O").family == "CHO" and space_of("C2H3Cl").family == "CHX"
    nitrogen = space_of("C5H5N")
    assert nitrogen.dbe == 4 and nitrogen.family == "CHN"
    # homologues share a Kendrick mass defect
    assert space_of("C10H22").kmd == pytest.approx(space_of("C14H30").kmd, abs=1e-6)


def test_unreadable_or_carbonless_formulae_are_refused_not_guessed():
    for text in ("", "   ", "H2O", "C6H5[OH]", "Xx3", "c6h6", "C6H6+", None):
        assert space_of(text) is None
    assert parse_formula("C2H6O") == {"C": 2, "H": 6, "O": 1}


# ---------------------------------------------------------------- the explorer
def test_the_largest_peaks_table_is_gone_and_the_default_plot_is_the_2d_map(qtbot, window, tmp_path):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    assert not hasattr(win, "summary")
    ex = win.explorer
    assert ex.current_kind() == "map" and ex.stack.currentWidget() is ex.map_plot
    assert ex.run_map is not None and ex.image.isVisible()
    found = state.current_for(w.index)
    assert len(ex.peaks) == len(found.matches)
    assert sum(len(item.data["x"]) for item in ex.dots.values()) == len(found.matches)


def test_selecting_a_peak_boxes_it_on_the_map_and_clicking_a_dot_selects_it(qtbot, window, tmp_path):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    ex = win.explorer
    shown = win._shown
    peak = next(p for p in ex.peaks if p.label == shown.query.label)
    rect = ex.box.rect()
    assert ex.box.isVisible() and rect.left() < peak.t1 < rect.right() and rect.top() < peak.t2 < rect.bottom()
    other = next(m for m in state.current_for(w.index).matches if m.query.label != shown.query.label)
    ex.peakClicked.emit(other.query.label)
    assert win._shown.query.label == other.query.label
    moved = next(p for p in ex.peaks if p.label == other.query.label)
    assert ex.ring.data["x"][0] == pytest.approx(moved.t1)


def test_filtered_out_peaks_stay_on_the_map_as_faint_dots(qtbot, window, tmp_path):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    win.name_filter.setText("compound50")
    ex = win.explorer
    assert len(ex.peaks) == len(state.current_for(w.index).matches)
    assert 0 < sum(p.shown for p in ex.peaks) < len(ex.peaks)


def test_every_explorer_plot_draws_and_the_close_up_zooms_to_the_peak(qtbot, window, tmp_path):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    ex = win.explorer
    assert [ex.kind.itemData(i) for i in range(ex.kind.count())] == [k for k, _t in PLOTS]
    for i in range(ex.kind.count()):
        ex.kind.setCurrentIndex(i)
        assert ex.note.text()
    ex.kind.setCurrentIndex(ex.kind.findData("zoom"))
    p = next(p for p in ex.peaks if p.label == win._shown.query.label)
    (x0, x1), (y0, y1) = ex.map_plot.getPlotItem().viewRange()
    assert x0 < p.t1 < x1 and x1 - x0 < 3.0 and y0 < p.t2 < y1
    ex.kind.setCurrentIndex(ex.kind.findData("candidates"))
    assert "candidate" in ex.note.text().lower() or "retention index" in ex.note.text().lower()


def test_the_ri_plot_has_a_legend_and_a_rank_axis(qtbot, window, tmp_path):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    labels = [label.text for _s, label in win.ri_legend.items]
    assert len(labels) == 3 and any("Peak's RI" in t for t in labels) and any("outside the RI gate" in t for t in labels)
    assert "Candidate rank" in win.ri_plot.getPlotItem().getAxis("left").labelText


# ---------------------------------------------------------------- naming
def test_double_clicking_a_name_opens_the_naming_dialog_for_that_peak(qtbot, window, tmp_path, monkeypatch):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    opened = []
    monkeypatch.setattr("pygcgc.gui.peak_naming.rename_peak", lambda parent, project, index, label: opened.append((index, label)))
    win._cell_double_clicked(0, 4)
    assert opened == [(w.index, win.table.item(0, 0).data(256))]
    win._cell_double_clicked(0, 2)
    assert len(opened) == 1                                  # only the name cell does this


def test_a_name_in_no_library_is_displayed_exported_and_never_borrows_a_formula(qtbot, window, tmp_path, monkeypatch):
    import csv
    w = window
    win, state = ready(qtbot, w, tmp_path)
    win.table.selectRow(0)
    label = win._shown.query.label
    dialog_returning(monkeypatch, "my own compound")
    from pygcgc.gui.peak_naming import rename_peak
    assert rename_peak(w, w.project, w.index, label)
    d = state.decisions_for(w.index)[label]
    assert d["name"] == "my own compound" and d["formula"] == "" and d["cas"] == ""
    win._fill_results()
    row = next(r for r in range(win.table.rowCount()) if win.table.item(r, 0).data(256) == label)
    assert win.table.item(row, 4).text() == "my own compound" and win.table.item(row, 5).text() == ""
    assert "named by you" in win.detail_title.text() and "my own compound" in win.detail_title.text()
    out = tmp_path / "ids.csv"
    win.export_to(out)
    exported = {r["peak"]: r for r in csv.DictReader(out.open(encoding="utf-8"))}[label]
    assert exported["decision_name"] == "my own compound" and exported["decision_formula"] == ""


def test_a_library_name_typed_by_hand_keeps_that_compounds_formula(qtbot, window, tmp_path, monkeypatch):
    from pygcgc.gui.peak_naming import rename_peak
    w = window
    win, state = ready(qtbot, w, tmp_path)
    label = win.table.item(0, 0).data(256)
    dialog_returning(monkeypatch, "COMPOUND 51")
    rename_peak(w, w.project, w.index, label)
    d = state.decisions_for(w.index)[label]
    assert d["name"] == "Compound51" and d["formula"] == "C6H6O"
    match = next(m for m in state.current_for(w.index).matches if m.query.label == label)
    assert identity_of(match, d) == ("Compound51", "C6H6O", d["cas"])


def test_a_peak_with_no_candidate_can_still_be_named_cleared_and_promoted(qtbot, window, tmp_path, monkeypatch):
    from pygcgc.core.targets import load_targets
    from pygcgc.gui.peak_naming import rename_peak
    w = window
    win, state = ready(qtbot, w, tmp_path)
    match = state.current_for(w.index).matches[0]
    match.candidates = []
    win._fill_results()
    win._select_label(match.query.label)
    assert win._shown is match and win.rename_btn.isEnabled() and not win.accept_btn.isEnabled()
    assert win._display_name(match, {}) == ""
    dialog_returning(monkeypatch, "hand-named thing")
    assert rename_peak(w, w.project, w.index, match.query.label)
    d = state.decisions_for(w.index)[match.query.label]
    assert d["name"] == "hand-named thing"
    win._select_label(match.query.label)
    assert "hand-named thing" in win.detail_title.text()
    path = tmp_path / "targets.csv"
    match.query.mz = np.array([50, 60, 70])
    match.query.intensity = np.array([10.0, 999.0, 400.0])
    win.promote_to_file(match, None, path, d)
    target = load_targets(path)[0]
    assert target.compound == "hand-named thing" and target.quant_ion == 60 and target.t1 == pytest.approx(match.query.t1)
    win.clear_btn.click()
    assert match.query.label not in state.decisions_for(w.index)


def test_a_rejected_or_unknown_peak_has_no_identity():
    class Match:
        candidates = []
        best = None
    assert identity_of(Match(), {"action": "rejected"}) == ("", "", "")
    assert identity_of(Match(), None) == ("", "", "")


# ---------------------------------------------------------------- chemical space
def test_chemical_space_places_named_peaks_and_counts_the_rest_with_reasons(qtbot, window, tmp_path, monkeypatch):
    from pygcgc.gui.peak_naming import rename_peak
    w = window
    win, state = ready(qtbot, w, tmp_path)
    found = state.current_for(w.index)
    rows, skipped = space_rows(found, state.decisions_for(w.index))
    assert len(rows) + sum(skipped.values()) == len(found.matches)
    assert rows and all(r.space.formula in ("C10H22", "C6H6O") for r in rows)
    label = rows[0].label
    dialog_returning(monkeypatch, "my own compound")
    rename_peak(w, w.project, w.index, label)
    rows2, skipped2 = space_rows(found, state.decisions_for(w.index))
    assert label not in [r.label for r in rows2] and skipped2["no formula (a name no library holds)"] == 1
    state.decide(w.index, label, "rejected")
    _rows3, skipped3 = space_rows(found, state.decisions_for(w.index))
    assert skipped3["no name"] == skipped["no name"] + 1


def test_the_chemical_space_window_draws_every_plot_and_clicking_selects_the_peak(qtbot, window, tmp_path):
    w = window
    win, state = ready(qtbot, w, tmp_path)
    win.space_btn.click()
    cs = win.chemspace
    assert isinstance(cs, ChemSpaceWindow) and cs.isVisible() and cs.rows
    for i in range(cs.kind.count()):
        cs.kind.setCurrentIndex(i)
        for j in range(cs.colour.count()):
            cs.colour.setCurrentIndex(j)
        assert cs.note.text() and cs.about.text()
    cs.kind.setCurrentIndex(cs.kind.findData("vk"))
    other = next(r for r in cs.rows if r.label != win._shown.query.label)
    cs.peakClicked.emit(other.label)
    assert win._shown.query.label == other.label
    assert cs.selected == other.label and cs._ring is not None
    cs.show_box.setCurrentIndex(cs.show_box.findData("named"))
    assert "0 of" in cs.note.text()
    win.space_btn.click()
    assert win.chemspace is cs                                   # one window, brought forward


def test_the_chemical_space_window_follows_a_new_name(qtbot, window, tmp_path, monkeypatch):
    from pygcgc.gui.peak_naming import rename_peak
    w = window
    win, state = ready(qtbot, w, tmp_path)
    win.space_btn.click()
    cs = win.chemspace
    before = len(cs.rows)
    label = cs.rows[0].label
    dialog_returning(monkeypatch, "my own compound")
    rename_peak(w, w.project, w.index, label)
    assert len(cs.rows) == before - 1


def test_the_window_says_so_when_nothing_is_matched_yet(qtbot, window):
    w = window
    w.set_run(w.project.index_of("RUN_A"))
    w.show_library_matcher()
    w.library_matcher.space_btn.click()
    cs = w.library_matcher.chemspace
    assert not cs.rows and "Press Match peaks" in cs.note.text()
    assert isinstance(cs, QtWidgets.QDialog)
