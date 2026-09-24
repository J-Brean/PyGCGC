"""Native source selection, persistence, cancellation and identification policy."""
import os
import threading
from pathlib import Path

os.environ.setdefault('QT_QPA_PLATFORM','offscreen')
import numpy as np
import pytest
pytest.importorskip('PySide6')

from pygcgc.gui.jobs import JobManager
from pygcgc.gui.project import Project
from tests.test_markes_lsc import make_lsc


def test_folder_selection_and_saved_source(qtbot,tmp_path):
    raw=make_lsc(tmp_path/'raw.lsc')
    low=make_lsc(tmp_path/'low.lsc',energy=14)
    make_lsc(tmp_path/'dbc.lsc',dbc=True)
    project=Project(JobManager())
    done=[]
    project.open_folder(tmp_path,on_loaded=lambda:done.append(True))
    qtbot.waitUntil(lambda: bool(done),timeout=20000)
    assert [e.path for e in project.entries] == [raw]
    assert project.entries[0].import_meta['energy_ev'] == 70
    project.add_files([low],on_loaded=lambda:done.append(True))
    qtbot.waitUntil(lambda:len(done)==2,timeout=20000)
    saved=tmp_path/'project.json'
    project.save(saved)
    reopened=Project(JobManager())
    reopened.load(saved,on_loaded=lambda:done.append(True))
    qtbot.waitUntil(lambda:len(done)==3,timeout=20000)
    assert [e.import_meta['energy_ev'] for e in reopened.entries] == [70,14]
    assert 'energies' in reopened.quant.ready()[1]
    for a,b in zip(project.entries,reopened.entries):
        np.testing.assert_array_equal(a.run.times,b.run.times)
        np.testing.assert_array_equal(a.run.tic,b.run.tic)


def test_sequence_explicitly_selects_alternative(qtbot,tmp_path):
    make_lsc(tmp_path/'raw.lsc')
    low=make_lsc(tmp_path/'low.lsc',energy=14)
    (tmp_path/'sequence.csv').write_text('file,type,level\nlow.lsc,sample,\n')
    project=Project(JobManager())
    done=[]
    project.open_folder(tmp_path,on_loaded=lambda:done.append(True))
    qtbot.waitUntil(lambda:bool(done),timeout=20000)
    assert [e.path for e in project.entries] == [low]


def test_cancel_import_resets_loading(qtbot,tmp_path,monkeypatch):
    started=threading.Event()
    def wait_for_cancel(path,check_cancel=None,progress=None):
        started.set()
        while True:
            check_cancel()
            threading.Event().wait(.005)
    monkeypatch.setattr('pygcgc.gui.project.read_run',wait_for_cancel)
    jobs=JobManager()
    project=Project(jobs)
    project.add_files([tmp_path/'slow.lsc'])
    qtbot.waitUntil(started.is_set,timeout=20000)
    jobs.cancel_all()
    qtbot.waitUntil(lambda:not project.loading,timeout=20000)
    assert project.entries[0].run is None
    assert 'cancelled' in project.entries[0].error


def test_changed_saved_source_is_reported(qtbot,tmp_path):
    path=make_lsc(tmp_path/'run.lsc')
    project=Project(JobManager())
    done=[]
    project.add_files([path],on_loaded=lambda:done.append(True))
    qtbot.waitUntil(lambda:bool(done),timeout=20000)
    saved=tmp_path/'project.json'
    project.save(saved)
    make_lsc(path,acquisition=300)
    reopened=Project(JobManager())
    reopened.load(saved,on_loaded=lambda:done.append(True))
    qtbot.waitUntil(lambda:len(done)==2,timeout=20000)
    assert reopened.entries[0].run is None
    assert 'changed' in reopened.entries[0].error


def test_real_native_window_renders(qtbot):
    from pygcgc.gui.mainwindow import MainWindow
    folder=Path(os.environ.get('PYGCGC_REAL_DATA',Path(__file__).resolve().parents[1]/'real_data'))
    paths=list(folder.glob('*_70eV.lsc'))
    if not paths: pytest.skip('Local native acquisition unavailable')
    window=MainWindow()
    qtbot.addWidget(window)
    window.show()
    done=[]
    window.project.add_files(paths[:1],on_loaded=lambda:done.append(True))
    qtbot.waitUntil(lambda:bool(done) and window.view2d.result is not None,timeout=60000)
    from pygcgc.core.io.markes_lsc import LSCFile
    with LSCFile(paths[0], validate_index=False) as source:
        assert window.project.entries[0].run.n_scans == len(source.offsets)
    assert not window.project.period_confirmed
    window.close()
