"""🏷️ Labelling an internal standard on an MS-aware project that also holds legacy blobs."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")
pytest.importorskip("pyqtgraph")

from tests.test_gui import TIMEOUT, window  # noqa: E402,F401


def test_empty_ms_aware_result_does_not_fall_back_to_legacy_peaks(qtbot, window):
    from pygcgc.core.components import ComponentSet
    from pygcgc.gui.blob_state import MS_AWARE
    w = window
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.was_run, timeout=TIMEOUT)
    i = w.project.index_of("RUN_A")
    assert len(w.blobs.blobs_for(i)) > 0
    w.blobs.set_method(MS_AWARE)
    assert w.istd._peaks_for(i) is None
    w.blobs.component_sets_by_name[w.project.entries[i].name] = ComponentSet()
    assert len(w.istd._peaks_for(i)) == 0


def test_label_uses_the_selected_peak_method(qtbot, window):
    # Legacy blobs appear on an MS-aware project once 🫧 or ⛓️ has run. The click picks a
    # component id, so the label must be read from the components, not the legacy blobs
    # (it used to fail with "Couldn't read that peak").
    from pygcgc.gui.blob_state import MS_AWARE

    w = window
    i = w.project.index_of("RUN_A")
    w.set_run(i)
    w.run_blobs()                                     # legacy (the fixture's method)
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.blobs_for(i) is not None, timeout=TIMEOUT)
    w.blobs.set_method(MS_AWARE)
    w.run_blobs()
    qtbot.waitUntil(lambda: not w.blobs.running and w.blobs.components_for(i) is not None
                    and bool(w.blobs.components_for(i).components), timeout=TIMEOUT)
    component = w.blobs.components_for(i).components[0]
    ref = w._nearest_peak_ref(i, component.t1, component.t2)
    assert isinstance(ref, str)
    name = w.istd.library[0].name
    w.istd.label(name, i, ref)
    qtbot.waitUntil(lambda: name in w.istd.references, timeout=TIMEOUT)
    assert w.istd.references[name].t1 == pytest.approx(component.t1, abs=w.project.period)
