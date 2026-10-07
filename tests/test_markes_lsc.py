"""Structural failures, independent values, and opt-in/local real-file checks."""
import json
import os
import struct
from pathlib import Path

import numpy as np
import pytest

from pygcgc.core.io import discover_runs, read_run
from pygcgc.core.io.markes_lsc import LSCFile, parse_metadata
from pygcgc.core.io.andi import write_andi


def make_lsc(path, *, energy=70, dbc=False, acquisition=100, spectra=None, spectrum_tag=b'SD04'):
    if spectra is None:
        spectra = [(0.98, [(57010, 10.), (57200, 20.), (73000, 5.)]),
                   (1., []), (1.02, [(91000, 7.)])]
    blob = bytearray(b'FT04'+struct.pack('<IQ',16,0))
    index = {}
    def chunk(tag,payload):
        index.setdefault(tag,[]).append(len(blob))
        blob.extend(tag+struct.pack('<I',len(payload)+8)+payload)
    pd = ('[Mass Conversion Settings]\nMass resolution=MASSRESOLUTION_1AMU_CENTROIDED\n'
          'Min mass=38\nMax mass=300\n')
    if dbc: pd += '[DBC]\nModule=DynamicBackgroundCompensation\n'
    chunk(b'PD01',pd.encode())
    sm = f'[Start of measurement]\nsystemtime ticks={acquisition}\nmode caption 0=_70eV\nmode caption 1=_14eV\n'
    chunk(b'SM01',sm.encode())
    for time, ions in spectra:
        second = int(time)
        payload = struct.pack('<5I', round((time-second)*2**32), acquisition+second,
                              0x1020 if energy == 70 else 0x2020,len(ions),0)
        payload += b''.join(struct.pack('<IdI',m,i,0) for m,i in ions)
        chunk(spectrum_tag,payload)
    chunk(b'MD01',sm.encode())
    directory = len(blob)
    payload = b''.join(tag+struct.pack('<I',len(offsets))+struct.pack(f'<{len(offsets)}Q',*offsets)
                       for tag,offsets in index.items())
    chunk(b'DI02',payload)
    struct.pack_into('<Q',blob,8,directory)
    path.write_bytes(blob)
    return path


def test_native_values_rounding_rollover_and_export(tmp_path):
    path = make_lsc(tmp_path/'run.lsc')
    run = read_run(path)
    np.testing.assert_allclose(run.times,[.98,1,1.02],atol=1/2**32)
    np.testing.assert_array_equal(run.tic,[35,0,7])
    assert run.data[0,57-38] == 30
    assert run.meta['merged_points'] == 1
    assert run.meta['energy_ev'] == 70
    assert run.mz.tolist() == list(range(38,301))
    with LSCFile(path) as source:
        np.testing.assert_array_equal(source.spectrum(0)[2], [57.01,57.2,73])
    reread = read_run(write_andi(run,tmp_path/'export.cdf'))
    np.testing.assert_array_equal(reread.data.toarray(),run.data.toarray())
    np.testing.assert_array_equal(reread.times,run.times)


def test_discovery_prefers_raw_but_keeps_ambiguous_acquisitions(tmp_path):
    raw = make_lsc(tmp_path/'arbitrary.lsc')
    low = make_lsc(tmp_path/'another.lsc',energy=14)
    corrected = make_lsc(tmp_path/'processed.lsc',dbc=True)
    assert discover_runs(tmp_path) == [raw]
    assert set(discover_runs(tmp_path,include_alternatives=True)) == {raw,low,corrected}
    other = make_lsc(tmp_path/'other.lsc',acquisition=200)
    assert set(discover_runs(tmp_path)) == {raw,other}
    duplicate = make_lsc(tmp_path/'duplicate.lsc')
    assert set(discover_runs(tmp_path)) == {raw,low,corrected,other,duplicate}


def test_discovery_leaves_out_tof_ds_reexports(tmp_path):
    # "Open TOF Data File" copies (SD01 spectra, often named "<run> 70eV.lsc") sit
    # beside the originals; folder opening skips them instead of failing to load them.
    raw = make_lsc(tmp_path/'run_70eV.lsc')
    make_lsc(tmp_path/'run 70eV.lsc', spectrum_tag=b'SD01')
    assert discover_runs(tmp_path) == [raw]
    with pytest.raises(ValueError, match='no SD04'):
        read_run(tmp_path/'run 70eV.lsc')


@pytest.mark.parametrize('damage', ['truncated','index','length','count','reserved','nan','order','timestamp','version','channel'])
def test_corruption_rejected(tmp_path,damage):
    path=make_lsc(tmp_path/'bad.lsc')
    with LSCFile(path) as source: p=int(source.offsets[0]); p2=int(source.offsets[1])
    b=bytearray(path.read_bytes())
    if damage=='truncated': b=b[:-1]
    elif damage=='index': struct.pack_into('<Q',b,8,len(b)+100)
    elif damage=='length': struct.pack_into('<I',b,p+4,7)
    elif damage=='count': struct.pack_into('<I',b,p+20,2**32-1)
    elif damage=='reserved': struct.pack_into('<I',b,p+40,1)
    elif damage=='nan': struct.pack_into('<d',b,p+32,float('nan'))
    elif damage=='order': struct.pack_into('<I',b,p+44,56000)
    elif damage=='timestamp': b[p2+8:p2+16]=b[p+8:p+16]
    elif damage=='version': b[:4]=b'FT05'
    elif damage=='channel': struct.pack_into('<I',b,p2+16,0x2020)
    path.write_bytes(b)
    with pytest.raises(ValueError,match='byte'): read_run(path)
    # Closed even on failure: Windows permits rename after the reader returns.
    path.rename(tmp_path/'closed.lsc')


def test_cancel_releases_source(tmp_path):
    path=make_lsc(tmp_path/'cancel.lsc')
    class Stop(Exception): pass
    def cancel(): raise Stop()
    with pytest.raises(Stop): read_run(path,check_cancel=cancel)
    path.rename(tmp_path/'closed.lsc')


def test_metadata_duplicate_keys_preserved():
    assert parse_metadata('[A]\nName=first\nName=second') == [['A','Name','first'],['A','Name','second']]


def test_import_identity_in_provenance_and_cache(tmp_path):
    from pygcgc.core.report import build_provenance, write_html_report
    from pygcgc.core.batch import _run_cache_base
    run=read_run(make_lsc(tmp_path/'run.lsc'))
    before=_run_cache_base(run,4,0)
    provenance=build_provenance('test',runs=[run])
    assert provenance.sources[0]['source_sha256'] == run.meta['source_sha256']
    html=write_html_report(tmp_path/'report.html',provenance).read_text()
    assert run.meta['source_sha256'] in html
    run.meta['parser_version'] += 1
    assert _run_cache_base(run,4,0) != before


def test_native_batch_selects_one_channel(tmp_path):
    from pygcgc.core.batch import run_batch, BatchConfig, BatchSettings
    spectra=[(i/50,[(57000,100. if i % 200 in range(50,65) else 2.)]) for i in range(2000)]
    make_lsc(tmp_path/'raw.lsc',spectra=spectra)
    make_lsc(tmp_path/'low.lsc',energy=14,spectra=spectra)
    make_lsc(tmp_path/'dbc.lsc',dbc=True,spectra=spectra)
    result=run_batch(tmp_path,BatchConfig(period=4,settings=BatchSettings(run_quant=False,run_blobs=False)),max_workers=1)
    assert not result.errors
    assert list(result.prepped) == ['raw']


def test_quantification_rejects_mixed_energy(tmp_path):
    from pygcgc.core.quant import quantify_prepared
    raw=read_run(make_lsc(tmp_path/'raw.lsc'))
    low=read_run(make_lsc(tmp_path/'low.lsc',energy=14))
    with pytest.raises(ValueError,match='energies'):
        quantify_prepared([(raw,raw,None,None),(low,low,None,None)],[],[],{})


@pytest.mark.parametrize('suffix', ['_70eV.lsc','_14eV.lsc','_70eV_0.1_secs.dbc.lsc'])
def test_complete_real_files(suffix):
    folder=Path(os.environ.get('PYGCGC_REAL_DATA',Path(__file__).resolve().parents[1]/'real_data'))
    audit_path=folder/'LSC_AUDIT.json'
    if not audit_path.exists(): pytest.skip('Local real-data audit not available')
    references=json.loads(audit_path.read_text())['files']
    reference=next(r for r in references if r['file'].endswith(suffix))
    path=folder/reference['file']
    if not path.exists(): pytest.skip('Local instrument file not available')
    run=read_run(path)
    assert run.n_scans == reference['scans']
    assert run.meta['source_sha256'] == reference['sha256']
    assert run.meta['stored_points'] == reference['points']
    assert run.meta['merged_points'] == reference['duplicate_unit_mass_count']
    assert run.tic.sum() == reference['total_signal']
    assert run.duration == reference['duration_s']
    assert np.count_nonzero(run.tic == 0) == reference['empty_scans']
    # Independent per-scan audit CSV, including every TIC and five unit-ion traces.
    csv_path=folder/(path.stem+'.chromatogram.csv')
    if csv_path.exists():
        expected=np.loadtxt(csv_path,delimiter=',',skiprows=1)
        np.testing.assert_allclose(run.times-run.times[0],expected[:,1],atol=1e-10)
        np.testing.assert_array_equal(run.tic,expected[:,2])
        for j,mass in enumerate([57,73,91,207,281]):
            np.testing.assert_array_equal(run.data[:,mass-38].toarray().ravel(),expected[:,3+j])
    if suffix == '_70eV.lsc':
        from pygcgc.core.io.andi import crop_run
        from pygcgc.core.fold import Fold, FoldSpec
        from pygcgc.core.background import remove_background
        from pygcgc.core.timing import check_timing
        from pygcgc.core.peaks import find_blobs
        report = check_timing(run.times)
        assert report.collisions == 0 and report.missing_scans == 0
        cropped = crop_run(run,650,690)
        original = cropped.data.copy()
        # 4 s is a numerical smoke-test parameter, not the instrument's confirmed period.
        fold = Fold(cropped.times,FoldSpec(period=4))
        np.testing.assert_allclose(fold.tic_image(cropped).sum(),cropped.tic.sum(),rtol=1e-6)
        clean = remove_background(cropped,fold)
        assert (cropped.data != original).nnz == 0
        assert clean.run.data.shape == cropped.data.shape
        assert np.all(np.isfinite(clean.run.data.data))
        # Background residuals intentionally remain signed so noise is unbiased.
        assert np.all(np.isfinite(clean.noise.sigma))
        assert np.all(clean.noise.sigma > 0)
        np.testing.assert_array_equal(clean.run.times,cropped.times)
        find_blobs(clean.run,fold)
