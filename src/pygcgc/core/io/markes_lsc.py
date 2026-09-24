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


HEADER = np.dtype([('tag', 'S4'), ('size', '<u4'), ('frac', '<u4'), ('sec', '<u4'),
                   ('flag', '<u4'), ('count', '<u4'), ('reserved', '<u4')])
BLOCK = 16384      # scans decoded together: bounds memory, and sets how often cancel and progress are checked


def read_lsc(path, *, check_cancel=None, progress=None):
    """Read this supported LSC variant into a Run, preserving every scan.

    The file is read into memory once (the SHA-256 is taken of the same bytes)
    and decoded a block of scans at a time with array operations, making the
    checks `LSCFile.spectrum` makes for a single scan. A 460,000-scan file takes
    a second or two; the scan-by-scan loop it replaced took about 25 s."""
    check = check_cancel or (lambda: None)
    name = Path(path).name
    with LSCFile(path, check) as source:
        n = len(source.offsets)
        offsets = source.offsets.astype(np.int64)
        settings = {k.lower():v for s,k,v in source.records['PD01'] if s == 'Mass Conversion Settings'}
        lo, hi = int(float(settings.get('min mass', 38))), int(float(settings.get('max mass', 300)))
        if not 0 <= lo <= hi <= 100000:
            source.fail('invalid declared mass range')
        source.f.seek(0)
        raw = source.f.read()
        if len(raw) != source.size:
            source.fail('file truncated during import')
        buf = np.frombuffer(raw, np.uint8)

        def first(bad, what, scan_of=None):
            """Fail at the first flagged scan, naming it and its byte offset."""
            k = int(np.flatnonzero(bad)[0])
            i = k if scan_of is None else int(scan_of(k))
            source.fail(f'{what} in scan {i}', int(offsets[i]))

        # Spectrum headers (28 bytes each).
        short = offsets + 28 > source.size
        if short.any():
            first(short, 'spectrum header outside the file')
        heads = np.empty(n, HEADER)
        for a in range(0, n, BLOCK):
            check()
            p = offsets[a:a + BLOCK]
            heads[a:a + BLOCK] = buf[p[:, None] + np.arange(28)].view(HEADER).ravel()
        counts = heads['count'].astype(np.int64)
        bad = ((heads['tag'] != b'SD04') | (heads['size'].astype(np.int64) != 28 + 16 * counts)
               | (heads['reserved'] != 0) | (offsets + 28 + 16 * counts > source.size))
        if bad.any():
            first(bad, 'invalid or unsupported spectrum header (tag, length, count or reserved field)')
        flags = heads['flag']
        if np.any(flags & 0xf000 != flags[0] & 0xf000):
            first(flags & 0xf000 != flags[0] & 0xf000, 'mixed energy channels; unsupported LSC variant')
        times = (heads['sec'].astype(np.int64) - source.origin) + heads['frac'] / 2**32
        if n > 1 and np.any(np.diff(times) <= 0):
            first(np.diff(times) <= 0, 'non-increasing timestamp', scan_of=lambda k: k + 1)

        # Ion records (16 bytes each), a block of scans at a time.
        indptr = np.concatenate(([0], np.cumsum(counts)))
        values = np.empty(int(indptr[-1]), np.float32)
        indices = np.empty(len(values), np.int32)
        off_nominal = 0
        f32_max = np.finfo(np.float32).max
        for a in range(0, n, BLOCK):
            check()
            if progress:
                progress(a, n, f'Decoding {name}')
            b = min(a + BLOCK, n)
            A, B = int(indptr[a]), int(indptr[b])
            if A == B:
                continue
            ions = np.concatenate([buf[s:s + 16 * c] for s, c in zip(offsets[a:b] + 28, counts[a:b]) if c]).view(ION)
            local = indptr[a:b + 1] - A
            scan_of = lambda k, local=local, a=a: a + np.searchsorted(local, k, side='right') - 1
            mass = ions['mass'].astype(np.int64)
            steps = np.diff(mass) <= 0
            steps[local[1:-1][(local[1:-1] > 0) & (local[1:-1] < len(mass))] - 1] = False   # a new scan may restart
            if np.any(ions['reserved']) or steps.any():
                first((ions['reserved'] != 0) | np.append(steps, False), 'unsupported ion layout', scan_of)
            intensity = ions['intensity']
            if np.any(~np.isfinite(intensity) | (intensity <= 0)):
                first(~np.isfinite(intensity) | (intensity <= 0), 'invalid intensity', scan_of)
            mz = mass / 1000
            unit = np.rint(mz).astype(np.int64)
            if np.any((unit < lo) | (unit > hi)):
                first((unit < lo) | (unit > hi), 'mass outside declared range', scan_of)
            if np.any(intensity > f32_max):
                first(intensity > f32_max, 'intensity exceeds Run precision', scan_of)
            off_nominal += int(np.count_nonzero(np.abs(mz - unit) > .3))
            indices[A:B] = unit - lo
            values[A:B] = intensity
        data = csr_matrix((values, indices, indptr), shape=(n, hi-lo+1))
        original_points = data.nnz
        data.sum_duplicates()
        if not np.all(np.isfinite(data.data)):
            source.fail('unit-mass summation overflow')
        check()
        digest = hashlib.sha256(raw).hexdigest()
        del buf, raw
        stat = os.fstat(source.f.fileno())
        if (stat.st_size, stat.st_mtime_ns) != (source.size, source.mtime):
            source.fail('file changed during import; retry after acquisition/sync completes')
        flag_values, flag_n = np.unique(flags, return_counts=True)
        meta = dict(format='Markes LSC', parser_version=PARSER_VERSION,
                    source_sha256=digest, metadata_text=source.text,
                    metadata_records=source.records, flags={str(int(f)): int(c) for f, c in zip(flag_values, flag_n)},
                    stored_points=int(original_points), merged_points=int(original_points-data.nnz),
                    off_nominal_points=off_nominal, vendor_validated=False,
                    modulation_period_confirmed=False, **source.identity())
        if progress: progress(n, n, 'Read')
        return Run(times, np.arange(lo,hi+1), data, Path(path).stem, Path(path), meta)
