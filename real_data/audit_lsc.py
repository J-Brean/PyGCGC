"""Reproducible, read-only investigation of this FT04/SD04 LSC variant.

Run: python real_data/audit_lsc.py
Writes derived JSON/CSV/PNG beside this script; never changes input binaries.
This is an investigative decoder, not a registered PyGCGC reader.
"""
from pathlib import Path
from collections import Counter
import csv
import hashlib
import json
import mmap
import struct
import numpy as np

ION = np.dtype([('mass', '<u4'), ('intensity', '<f8'), ('unknown', '<u4')])

def inspect(path):
    with path.open('rb') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as b:
        chunks, scans, metadata = {}, [], {}
        p = 0
        while p < len(b):
            if p + 8 > len(b): raise ValueError('Truncated chunk header')
            tag, size = struct.unpack_from('<4sI', b, p)
            if size < 8 or p + size > len(b): raise ValueError('Invalid chunk size')
            tag = tag.decode('ascii')
            chunks.setdefault(tag, []).append(p)
            if tag == 'SD04': scans.append((p, size))
            elif tag in ('PD01', 'SM01', 'MD01'):
                metadata[tag] = b[p+8:p+size].decode('cp1252')
            p += size
        assert b[:4] == b'FT04'
        directory = struct.unpack_from('<Q', b, 8)[0]
        assert chunks['DI02'] == [directory]
        p = directory + 8
        end = directory + struct.unpack_from('<I', b, directory+4)[0]
        index = {}
        while p < end:
            tag, count = struct.unpack_from('<4sI', b, p)
            p += 8
            assert p + 8*count <= end
            offsets = struct.unpack_from(f'<{count}Q', b, p)
            p += 8*count
            index[tag.decode()] = list(offsets)
            if tag.decode() in chunks: assert list(offsets) == chunks[tag.decode()]
            else: assert all(x == 0 for x in offsets)
        assert p == end
        n = len(scans)
        ticks = np.zeros(n, dtype=np.uint64)
        counts = np.zeros(n, dtype=np.int64)
        flags = Counter()
        tic = np.zeros(n)
        eics = np.zeros((n, 5))
        masses_of_interest = [57, 73, 91, 207, 281]
        lo, hi, imin, imax = float('inf'), 0, float('inf'), 0
        fractional = off_nominal = duplicate = nonzero_ion_field = 0
        samples = {}
        for i, (p, size) in enumerate(scans):
            tick, flag, count, unknown = struct.unpack_from('<QIII', b, p+8)
            assert unknown == 0
            assert size == 28 + 16*count, (i,size,count)
            ticks[i], counts[i] = tick, count
            flags[flag] += 1
            ions = np.frombuffer(b, dtype=ION, count=count, offset=p+28)
            mass = ions['mass'].astype(float)/1000
            intensity = ions['intensity']
            assert np.all(np.isfinite(intensity)) and np.all(intensity > 0)
            assert np.all(np.diff(mass) > 0)
            nonzero_ion_field += int(np.count_nonzero(ions['unknown']))
            if count:
                lo, hi = min(lo,float(mass.min())), max(hi,float(mass.max()))
                imin, imax = min(imin,float(intensity.min())), max(imax,float(intensity.max()))
                fractional += int(np.count_nonzero(intensity != np.rint(intensity)))
                unit = np.rint(mass).astype(int)
                off_nominal += int(np.count_nonzero(abs(mass-unit) > .3))
                duplicate += int(np.count_nonzero(np.diff(unit) == 0))
                tic[i] = intensity.sum()
                for j,m in enumerate(masses_of_interest): eics[i,j] = intensity[unit == m].sum()
            if i in (15001,30000,100000,200000):
                samples[str(i)] = [[float(mass[k]),float(intensity[k])] for k in np.argsort(intensity)[-10:][::-1]]
            del ions, intensity
        # Subtract integer seconds before floating-point conversion: preserve precision.
        seconds = (ticks >> np.uint64(32)).astype(np.int64)
        fractions = (ticks & np.uint64(0xffffffff)).astype(float)/2**32
        times = seconds-seconds[0] + fractions-fractions[0]
        dt = np.diff(times)
        assert np.all(dt > 0)
        result = dict(file=path.name, bytes=len(b), sha256=hashlib.sha256(b).hexdigest(),
            chunks={k:len(v) for k,v in chunks.items()}, directory_offset=directory,
            index_entries={k:len(v) for k,v in index.items()}, scans=n, points=int(counts.sum()),
            empty_scans=int(np.count_nonzero(counts == 0)), first_nonempty_scan=int(np.flatnonzero(counts)[0]),
            first_nonempty_time_s=float(times[np.flatnonzero(counts)[0]]),
            first_timestamp_words=[int(ticks[0] & np.uint64(0xffffffff)),int(seconds[0])],
            duration_s=float(times[-1]), interval_s=dict(zip(['min','p01','median','p99','max'],map(float,np.percentile(dt,[0,1,50,99,100])))),
            mass_range=[lo,hi], intensity_range=[imin,imax], fractional_intensity_count=fractional,
            off_nominal_point_count=off_nominal, duplicate_unit_mass_count=duplicate,
            nonzero_ion_unknown_count=nonzero_ion_field, flags={hex(k):v for k,v in flags.items()},
            max_points_per_scan=int(counts.max()), total_signal=float(tic.sum()),
            apex_scan=int(tic.argmax()), apex_time_s=float(times[tic.argmax()]), sample_spectra=samples)
        for tag, value in metadata.items():
            (path.parent/(path.stem+'.'+tag+'.txt')).write_text(value,encoding='utf-8')
        return result, times, ticks, tic, eics

def main():
    folder = Path(__file__).parent
    results, arrays = [], {}
    for path in sorted(folder.glob('*.lsc')):
        result, times, ticks, tic, eics = inspect(path)
        results.append(result)
        arrays[path.name] = (times,ticks,tic,eics)
        print(json.dumps(result,indent=2),flush=True)
        with (folder/(path.stem+'.chromatogram.csv')).open('w',newline='') as f:
            writer=csv.writer(f)
            writer.writerow(['scan','seconds_from_first_scan','tic','mz57','mz73','mz91','mz207','mz281'])
            writer.writerows((i,times[i],tic[i],*eics[i]) for i in range(len(times)))
    names = list(arrays)
    raw = next(k for k in names if k.endswith('_70eV.lsc'))
    dbc = next(k for k in names if k.endswith('.dbc.lsc'))
    low = next(k for k in names if k.endswith('_14eV.lsc'))
    delta=(arrays[low][1].astype(np.int64)-arrays[raw][1].astype(np.int64))/2**32
    comparison=dict(all_raw_dbc_timestamps_equal=bool(np.array_equal(arrays[raw][1],arrays[dbc][1])),
                    channel_offset_seconds_percentiles=np.percentile(delta,[0,50,100]).tolist(),
                    dbc_tic_exceeds_raw_scans=int(np.count_nonzero(arrays[dbc][2]>arrays[raw][2])))
    (folder/'LSC_AUDIT.json').write_text(json.dumps(dict(files=results,comparison=comparison),indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(3,1,figsize=(12,9),constrained_layout=True)
    for ax,(name,(times,_,tic,eics)) in zip(axes,arrays.items()):
        ax.plot(times/60,tic,lw=.4)
        ax.set(title=name.split('IS1_')[-1],ylabel='Decoded TIC',xlabel='Minutes from first scan')
    fig.savefig(folder/'LSC_decoded_TIC.png',dpi=150)
    print('Comparison',comparison)

if __name__ == '__main__': main()
