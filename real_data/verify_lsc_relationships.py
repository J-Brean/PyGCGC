"""Independent cross-file checks supplementing audit_lsc.py."""
from pathlib import Path
import json
import mmap
import re
import struct
import numpy as np
from audit_lsc import ION

folder=Path(__file__).parent
paths=[next(folder.glob('*70eV.lsc')),next(folder.glob('*.dbc.lsc'))]
files=[p.open('rb') for p in paths]
maps=[mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ) for f in files]
offsets=[]
for b in maps:
    p=struct.unpack_from('<Q',b,8)[0]
    # Independently use the SD04 index, whose location audit_lsc checked.
    n=struct.unpack_from('<I',b,p+44)[0]
    offsets.append(np.frombuffer(b,dtype='<u8',count=n,offset=p+48).copy())
missing=greater=matched=0
for rp,dp in zip(*offsets):
    spectra=[]
    for b,p in zip(maps,[int(rp),int(dp)]):
        n=struct.unpack_from('<I',b,p+20)[0]
        spectra.append(np.frombuffer(b,dtype=ION,count=n,offset=p+28))
    raw,dbc=spectra
    positions=np.searchsorted(raw['mass'],dbc['mass'])
    valid=positions<len(raw)
    valid[valid] &= raw['mass'][positions[valid]]==dbc['mass'][valid]
    missing+=int(np.count_nonzero(~valid))
    matched+=int(np.count_nonzero(valid))
    greater+=int(np.count_nonzero(dbc['intensity'][valid]>raw['intensity'][positions[valid]]))
del spectra,raw,dbc
for b in maps: b.close()
for f in files: f.close()

dax=next(folder.glob('*.DAX')).read_bytes()
assert dax[:8]==b'DAX1    '
a=np.frombuffer(dax,offset=8,dtype=[('scan','<u4'),('offset','<u8')])
assert a[-1]['offset']==2**64-1
valid=a[:-1]
assert np.all(np.diff(valid['scan'])==64)
dat=next(folder.glob('*.DAT'))
with dat.open('rb') as f:
    ids=[]
    for entry in valid:
        f.seek(int(entry['offset']))
        ids.append(struct.unpack('<I',f.read(4))[0])
assert set(ids)=={9131}
lengths=np.diff(np.append(valid['offset'],np.uint64(dat.stat().st_size)))
header=next(folder.glob('*.HDR')).read_text()
section=header.split('[Mass Calibration Settings 70eV]')[1].split('\n[')[0]
values=dict(re.findall(r'^([^=\r\n]+)=([^\r\n]*)',section,re.M))
residuals=[]
for i in range(1,41):
    t=float(values[f'Time{i}'])
    predicted=float(values['Offset T0'])+float(values['Mass per Time'])*t+float(values['Mass per Time2'])*t*t
    residuals.append(predicted-float(values[f'Mass{i}']))
result=dict(dbc_exact_mass_matches=matched,dbc_masses_not_exactly_in_raw=missing,
    matched_dbc_intensities_greater_than_raw=greater,
    dax_records=len(a),dax_valid_blocks=len(valid),dax_final_scan=int(a[-1]['scan']),
    dat_all_indexed_blocks_measurement_id_9131=True,
    dat_block_bytes_percentiles=np.percentile(lengths,[0,50,100]).tolist(),
    calibration_70eV_quadratic_residual_Da_range=[min(residuals),max(residuals)])
(folder/'LSC_RELATIONSHIPS.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
