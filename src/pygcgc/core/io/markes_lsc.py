"""Native FT04/SD04 Markes line spectra, validated against TOF-DS 4.2 files.

The indexed centroid reader preserves original mass precision. read_lsc builds
the application's unit-mass CSR representation without applying calibration twice.
Binary semantics are inferred; metadata records that vendor validation is pending.
"""
from __future__ import annotations

import hashlib
import os
import re
import struct
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from ..run import Run

PARSER_VERSION = 1
ION = np.dtype([('mass', '<u4'), ('intensity', '<f8'), ('reserved', '<u4')])


def parse_metadata(text):
    """Preserve duplicate keys and section order as JSON-serializable records."""
    records, section = [], ''
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('[') and line.endswith(']'):
            section = line[1:-1]
        elif '=' in line and not line.startswith(';'):
            key, value = line.split('=', 1)
            records.append([section, key.strip(), value.strip()])
    return records


class LSCFile:
    """Context-managed indexed access; no persistent mappings or borrowed arrays."""
    def __init__(self, path, check_cancel=None, *, validate_index=True):
        self.path = Path(path)
        self.check_cancel = check_cancel or (lambda: None)
        self.validate_index = validate_index

    def fail(self, message, offset=0):
        raise ValueError(f'{self.path.name}: {message} at byte {offset}')

    def read(self, offset, size):
        if offset < 0 or size < 0 or offset + size > self.size:
            self.fail('read outside file bounds', offset)
        self.f.seek(offset)
        data = self.f.read(size)
        if len(data) != size:
            self.fail('file truncated during import', offset)
        return data

    def chunk(self, offset):
        tag, size = struct.unpack('<4sI', self.read(offset, 8))
        if size < 8 or offset + size > self.size:
            self.fail('invalid chunk length', offset)
        return tag, size

    def __enter__(self):
        self.f = self.path.open('rb')
        try:
            stat = os.fstat(self.f.fileno())
            self.size, self.mtime = stat.st_size, stat.st_mtime_ns
            if self.chunk(0) != (b'FT04', 16):
                self.fail('unsupported LSC header (expected FT04)')
            directory = struct.unpack('<Q', self.read(8, 8))[0]
            tag, length = self.chunk(directory)
            if tag != b'DI02' or directory + length != self.size:
                self.fail('invalid DI02 directory', directory)
            p, end, index = directory + 8, directory + length, {}
            while p < end:
                if p + 8 > end:
                    self.fail('truncated directory group', p)
                tag, count = struct.unpack('<4sI', self.read(p, 8))
                if tag in index or p + 8 + count * 8 > end:
                    self.fail('invalid directory count or duplicate group', p)
                index[tag] = np.frombuffer(self.read(p+8, count*8), '<u8').copy()
                p += 8 + count*8
            self.offsets = index.get(b'SD04', np.array([], dtype=np.uint64))
            if not len(self.offsets):
                self.fail('no SD04 spectra')
            # Check the complete framing independently of the directory.
            actual = {}
            p = 0 if self.validate_index else directory
            if not self.validate_index:
                actual = {k:v.tolist() for k,v in index.items() if k in (b'PD01', b'SM01', b'MD01', b'SD04')}
            while p < directory:
                if len(actual.get(b'SD04', [])) % 1024 == 0:
                    self.check_cancel()
                tag, length = self.chunk(p)
                if p + length > directory:
                    self.fail('chunk overlaps directory', p)
                if tag not in (b'FT04', b'PD01', b'SM01', b'MD01', b'SD04'):
                    self.fail(f'unsupported chunk {tag!r}', p)
                actual.setdefault(tag, []).append(p)
                p += length
            for tag, offsets in actual.items():
                if tag != b'FT04' and not np.array_equal(index.get(tag, []), offsets):
                    self.fail(f'directory disagrees with {tag!r} chunks')
            for tag, offsets in index.items():
                if tag not in actual and not (tag in (b'TS01', b'ST01') and np.all(offsets == 0)):
                    self.fail(f'unsupported directory group {tag!r}')
            self.text = {}
            for tag in (b'PD01', b'SM01', b'MD01'):
                entries = actual.get(tag, [])
                if len(entries) != 1:
                    self.fail(f'expected one {tag!r} metadata chunk')
                off = entries[0]
                _, size = self.chunk(off)
                self.text[tag.decode()] = self.read(off+8, size-8).decode('cp1252')
            self.records = {tag: parse_metadata(text) for tag, text in self.text.items()}
            pd = self.text['PD01']
            if 'MASSRESOLUTION_1AMU_CENTROIDED' not in pd:
                self.fail('unsupported spectrum mode; expected unit-mass centroids')
            starts = [v for s,k,v in self.records['SM01'] if k == 'systemtime ticks']
            if not starts:
                self.fail('missing acquisition clock origin')
            self.origin = int(starts[0])
            return self
        except BaseException:
            self.f.close()
            raise

    def __exit__(self, *exc):
        self.f.close()

    def spectrum(self, scan):
        """Return time, flags, original m/z and float64 intensity for one scan."""
        p = int(self.offsets[scan])
        tag, size = self.chunk(p)
        frac, sec, flag, count, reserved = struct.unpack('<5I', self.read(p+8, 20))
        if tag != b'SD04' or size != 28 + 16*count or reserved:
            self.fail(f'invalid or unsupported spectrum {scan}', p)
        ions = np.frombuffer(self.read(p+28, count*16), ION)
        mass = ions['mass'].astype(np.float64)/1000
        intensity = ions['intensity'].copy()
        if np.any(ions['reserved']) or np.any(np.diff(mass) <= 0):
            self.fail(f'unsupported ion layout in scan {scan}', p)
        if np.any(~np.isfinite(intensity)) or np.any(intensity <= 0):
            self.fail(f'invalid intensity in scan {scan}', p)
        return sec-self.origin+frac/2**32, flag, mass, intensity

    def identity(self):
        first = self.spectrum(0)[1]
        mode = (first >> 12) & 0xf
        captions = {k:v for s,k,v in self.records['SM01']}
        caption = captions.get(f'mode caption {mode-1}', '')
        match = re.fullmatch(r'_?(\d+(?:\.\d+)?)eV', caption)
        energy = float(match[1]) if match else None
        processed = 'DynamicBackgroundCompensation' in self.text['PD01']
        # Acquisition identity is independent of filenames and processing PD01.
        acquisition = hashlib.sha256(self.text['SM01'].encode()).hexdigest()
        return dict(acquisition_id=acquisition, energy_ev=energy,
                    processing_state='vendor_dbc' if processed else 'instrument',
                    time_origin_seconds=self.origin,
                    time_origin_policy='whole acquisition clock second; vendor RT offset unvalidated')


def read_lsc(path, *, check_cancel=None, progress=None):
    """Read this supported LSC variant into a Run, preserving every scan."""
    check = check_cancel or (lambda: None)
    with LSCFile(path, check) as source:
        n = len(source.offsets)
        counts = np.empty(n, np.int64)
        for i,p in enumerate(source.offsets):
            if i % 1024 == 0:
                check()
            _, size = source.chunk(int(p))
            count = struct.unpack('<I', source.read(int(p)+20, 4))[0]
            if size != 28+16*count:
                source.fail(f'point count mismatch in scan {i}', int(p))
            counts[i] = count
        indptr = np.concatenate(([0], np.cumsum(counts)))
        values = np.empty(int(indptr[-1]), np.float32)
        indices = np.empty(len(values), np.int32)
        times = np.empty(n, np.float64)
        flag_counts = {}
        lo, hi = 38, 300
        settings = {k.lower():v for s,k,v in source.records['PD01'] if s == 'Mass Conversion Settings'}
        lo, hi = int(float(settings.get('min mass', lo))), int(float(settings.get('max mass', hi)))
        if not 0 <= lo <= hi <= 100000:
            source.fail('invalid declared mass range')
        off_nominal = 0
        channel_bits = source.spectrum(0)[1] & 0xf000
        for i in range(n):
            if i % 1024 == 0:
                check()
                if progress: progress(i, n, f'Decoding {Path(path).name}')
            time, flag, mass, intensity = source.spectrum(i)
            if flag & 0xf000 != channel_bits:
                source.fail(f'mixed energy channels in scan {i}; unsupported LSC variant', int(source.offsets[i]))
            if i and time <= times[i-1]:
                source.fail(f'non-increasing timestamp in scan {i}', int(source.offsets[i]))
            unit = np.rint(mass).astype(np.int64)
            if np.any(unit < lo) or np.any(unit > hi):
                source.fail(f'mass outside declared range in scan {i}', int(source.offsets[i]))
            if np.any(intensity > np.finfo(np.float32).max):
                source.fail(f'intensity exceeds Run precision in scan {i}', int(source.offsets[i]))
            off_nominal += int(np.count_nonzero(abs(mass-unit) > .3))
            a,b = indptr[i:i+2]
            indices[a:b], values[a:b], times[i] = unit-lo, intensity, time
            flag_counts[str(flag)] = flag_counts.get(str(flag), 0)+1
        data = csr_matrix((values, indices, indptr), shape=(n, hi-lo+1))
        original_points = data.nnz
        data.sum_duplicates()
        if not np.all(np.isfinite(data.data)):
            source.fail('unit-mass summation overflow')
        digest = hashlib.sha256()
        source.f.seek(0)
        while block := source.f.read(4*1024*1024):
            check()
            digest.update(block)
        stat = os.fstat(source.f.fileno())
        if (stat.st_size, stat.st_mtime_ns) != (source.size, source.mtime):
            source.fail('file changed during import; retry after acquisition/sync completes')
        meta = dict(format='Markes LSC', parser_version=PARSER_VERSION,
                    source_sha256=digest.hexdigest(), metadata_text=source.text,
                    metadata_records=source.records, flags=flag_counts,
                    stored_points=int(original_points), merged_points=int(original_points-data.nnz),
                    off_nominal_points=off_nominal, vendor_validated=False,
                    modulation_period_confirmed=False, **source.identity())
        if progress: progress(n, n, 'Read')
        return Run(times, np.arange(lo,hi+1), data, Path(path).stem, Path(path), meta)
