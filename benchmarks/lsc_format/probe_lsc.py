"""Exploratory binary dumps; run with the project's Python environment."""
from pathlib import Path
import mmap
import struct
import numpy as np

for path in sorted((Path(__file__).resolve().parents[2] / 'real_data').glob('*.lsc')):
    with path.open('rb') as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as b:
        directory = struct.unpack_from('<Q', b, 8)[0]
        print('\nFILE', path.name, 'directory', directory, b[directory:directory+48].hex(' '))
        pos = 0
        scans = []
        while pos < len(b):
            tag, size = struct.unpack_from('<4sI', b, pos)
            assert size >= 8 and pos + size <= len(b)
            if tag == b'SD04': scans.append((pos, size))
            pos += size
        print('scans', len(scans), 'sizes', min(s for _,s in scans),max(s for _,s in scans))
        for i in [0,1,15000,30000,100000,200000,len(scans)-1]:
            p,s = scans[i]
            d = b[p+8:p+s]
            print(i,s, 'u32', struct.unpack_from('<5I',d), 'hex',d[:108].hex(' '))
            print('ions <IdI', [struct.unpack_from('<IdI',d,k) for k in range(20,min(len(d),100),16)])

