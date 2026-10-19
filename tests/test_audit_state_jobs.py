"""Exercise cancelled and obsolete GUI jobs without relying on thread timing."""

import os
from dataclasses import replace
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
from pygcgc.core.align import AnchorPoint
from pygcgc.core.hydrocarbons import AlkaneSeries
from tests.test_gui import window  # noqa: F401


@pytest.mark.parametrize("running", [False, True])
def test_removing_runs_refits_only_remaining_runs_and_never_refits_a_busy_job(qapp, tmp_path, monkeypatch, running):
    from pygcgc.gui.jobs import JobManager
    from pygcgc.gui.project import Project, RunEntry
    p = Project(JobManager())
    a, b = RunEntry(tmp_path/"a.cdf", "A"), RunEntry(tmp_path/"b.cdf", "B")
    p.entries = [a, b]
    q = p.quant
    q.runquants = {a: object(), b: object()}
    q.running = running
    refitted = []
    monkeypatch.setattr(q, "finish", lambda: refitted.append([r.name for r in q.info()]))
    p.remove([0])
    assert a not in q.runquants
    if running:
        assert not refitted and q.stale
    else:
        assert refitted == [["B"]]


@pytest.mark.parametrize("kind", ["classes", "clusters", "regions", "hydrocarbons"])
@pytest.mark.parametrize("ending", ["cancel", "changed", "removed"])
def test_analysis_jobs_settle_and_reject_obsolete_results(window, monkeypatch, kind, ending):
    w = window
    submitted, completed = [], []
    monkeypatch.setattr(w.jobs, "submit", lambda *a, **kw: submitted.append((a, kw)))
    if kind == "classes":
        state = w.classes
        assert state.classify(0, on_done=completed.append)
    elif kind in ("clusters", "regions"):
        state = w.clusters
        if kind == "regions":
            state.set_unit("regions", 3)
        assert state.cluster(0, on_done=completed.append)
    else:
        state = w.hydrocarbons
        state.find_series(0, [(15, 40., 1.1), (16, 90., 1.9)], on_done=completed.append)
    assert state.running
    args, kwargs = submitted[-1]
    if ending == "cancel":
        kwargs["on_cancel"]()
    else:
        if ending == "changed":
            w.project.entries[0].phase += .1
        else:
            w.project.entries.pop(0)
        args[2](None)  # Must reject before inspecting a result from obsolete inputs.
    assert not state.running and not completed


@pytest.mark.parametrize("change", ["phase", "settings", "method", "removed"])
def test_detection_does_not_publish_obsolete_peaks(window, monkeypatch, change):
    w = window
    submitted, completed = [], []
    monkeypatch.setattr(w.jobs, "submit", lambda *a, **kw: submitted.append((a, kw)))
    assert w.blobs.run_all(on_done=lambda: completed.append(True))
    if change == "phase":
        w.project.entries[0].phase += .1
    elif change == "settings":
        w.blobs.set_settings(replace(w.blobs.settings, snr_min=19))
    elif change == "method":
        w.blobs.set_method("ms-aware")
    else:
        w.project.entries.pop(0)
    for args, _kw in submitted:
        args[2](None)
    assert not w.blobs.running and not completed
    assert not w.blobs.by_run and not w.blobs.component_sets_by_name


@pytest.mark.parametrize("change", ["reference", "point", "phase", "cancel", "failed"])
def test_propagated_anchors_cannot_repopulate_changed_alignment(window, monkeypatch, change):
    w = window
    w.set_reference("RUN_A")
    w.align.add_anchor("RUN_A", "P1", 40, 1.1, 40, 1.1, 50)
    submitted, completed = [], []
    monkeypatch.setattr(w.jobs, "submit", lambda *a, **kw: submitted.append((a, kw)))
    w.align.propagate_to(w.project.loaded(), on_done=lambda: completed.append(True))
    if change == "reference":
        w.align.set_reference("RUN_B")
    elif change == "point":
        w.align.remove_anchor("RUN_A", "P1")
    elif change == "phase":
        w.project.entries[0].phase += .1
    before = w.align.snapshot()
    for number, (args, kwargs) in enumerate(submitted):
        if change == "cancel" and number == 0:
            kwargs["on_cancel"]()
        elif change == "failed" and number == 0:
            args[3]("read failed")
        else:
            args[2]({"P1": AnchorPoint("P1", 40, 1.1, 41, 1.2, 50)})
    assert w.align.snapshot() == before and not completed


@pytest.mark.parametrize("ending", ["cancel", "failed"])
def test_alkane_propagation_failure_does_not_leave_running_flag_set(window, monkeypatch, ending):
    w = window
    # Only its existence matters: the numerical worker is not run by this test.
    w.hydrocarbons.series["RUN_A"] = AlkaneSeries()
    submitted, completed = [], []
    monkeypatch.setattr(w.jobs, "submit", lambda *a, **kw: submitted.append((a, kw)))
    w.hydrocarbons.propagate_series(0, [1, 2], on_done=lambda: completed.append(True))
    assert w.hydrocarbons.running
    for args, kwargs in submitted:
        if ending == "cancel":
            kwargs["on_cancel"]()
        else:
            args[3]("bad data")
    assert not w.hydrocarbons.running and not completed


@pytest.mark.parametrize("kind", ["background", "blank", "siloxane", "phase", "period", "quant"])
def test_preparation_and_measurement_jobs_do_not_apply_after_fold_edits(window, monkeypatch, kind):
    from pygcgc.core.background import BackgroundSettings
    from pygcgc.core.siloxane import SiloxaneSettings
    w = window
    p = w.project
    submitted, completed = [], []
    monkeypatch.setattr(w.jobs, "submit", lambda *a, **kw: submitted.append((a, kw)))
    callback = lambda *args: completed.append(True)
    if kind in ("background", "blank"):
        p.compute_background([0], BackgroundSettings(subtract_blank=kind == "blank"),
                              "RUN_C_BLANK" if kind == "blank" else None, on_done=callback)
    elif kind == "siloxane":
        p.compute_siloxane([0], SiloxaneSettings(method="pixels"), on_done=callback)
    elif kind == "phase":
        p.estimate_phases(on_done=callback)
    elif kind == "period":
        p.auto_fold(on_done=callback)
    else:
        p.quant.requantify([p.entries[0]], on_done=callback)
    p.entries[0].phase += .1
    args, _kwargs = submitted[0]
    args[2]((None, [], p.period, []) if kind == "period" else None)
    assert not completed
    if kind == "quant":
        assert not p.quant.running and not p.quant.runquants


def test_run_metadata_edits_invalidate_quantification_in_flight(window, monkeypatch):
    p = window.project
    submitted, completed = [], []
    monkeypatch.setattr(p.jobs, "submit", lambda *a, **kw: submitted.append((a, kw)))
    p.quant.requantify([p.entries[0]], on_done=lambda: completed.append(True))
    p.entries[0].volume = 12.
    submitted[0][0][2](None)
    assert not p.quant.running and not p.quant.runquants and not completed


def test_decision_edits_during_measurement_do_not_start_overlapping_jobs(window, monkeypatch):
    p = window.project
    submitted = []
    monkeypatch.setattr(p.jobs, "submit", lambda *a, **kw: submitted.append((a, kw)))
    p.quant.requantify([p.entries[0]])
    p.quant.set_decisions(p.quant.decisions.copy(), requantify_runs=[p.entries[0].name])
    assert len(submitted) == 1
    submitted[0][0][2](None)
    assert not p.quant.running and p.quant.stale and not p.quant.runquants
