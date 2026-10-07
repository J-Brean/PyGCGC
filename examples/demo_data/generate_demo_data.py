"""Synthetic GC×GC-EI-ToF-MS demonstration data for PyGCGC.

    python examples/demo_data/generate_demo_data.py OUTPUT_FOLDER

Writes seven ANDI (.cdf) runs, four calibration standards and three samples, with
every supporting file the app reads, a small spectral library and the ground truth.
Everything is invented: the compounds are real ones, but their retention indices
and mass spectra are approximations written from memory of typical 70 eV spectra,
not copied from NIST. Use it to learn the app and to make screenshots, never to
check a method against a real answer.

What is in a run
----------------
- An 8 s modulation, 50 Hz, m/z 38 to 400, 48 minutes. Every run starts at its
  own phase (1.6 to 2.6 s), as a real acquisition does.
- About 150 named compounds across 17 classes (n-alkanes and their branched and
  cyclic relatives, alkenes, alkylbenzenes, indanes, naphthalenes, PAHs, oxy- and
  nitro-PAHs, N-heteroaromatics, an amine, alcohols, ketones, aldehydes, acids,
  esters, phenols, terpenes) placed by retention index on the n-alkane ladder.
- Five deuterated internal standards at 10 ng in every run, with the spectral
  signatures `core.calbuilder.ISTD_SIGNATURES` recognises. Perylene-d12 and
  perylene elute past the end of the modulation and wrap round.
- Impurities and artefacts: column bleed rising with the oven, cyclic siloxanes
  D4 to D8, phthalates, BHT and 2,4-di-tert-butylphenol from the tubes, a
  dichloromethane solvent front in the standards, an overloaded toluene peak
  smeared up the second dimension, a modulator stripe on m/z 40 and 44, air
  background, single-scan spikes, counting noise, a storage threshold and a
  detector that bends over and then saturates.
- Retention drifts a few seconds in t1 and a few hundredths of a second in t2
  from run to run. The diesel sample sits on an unresolved hump of hundreds of
  minor isomers.

The fixed seed makes every file identical on every machine.
"""

import csv
import json
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.stats import exponnorm, norm

from pygcgc import Run
from pygcgc.core.io import write_andi
from pygcgc.core.siloxane import BLEED_REFERENCE

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "demo_gcxgc")
SEED = 20261007
PERIOD = 8.0               # s
RATE = 50.0                # scans per second
RUN_END = 2880.0           # s (48 min)
MZ = np.arange(38, 401)
NG_TO_AREA = 3.0e5         # summed counts over all ions and scans for 1 ng at response 1
SOFT = 2.0e6               # the detector starts to bend over (v / (1 + v / SOFT))...
CEILING = 4.0e5            # ...and saturates here
NOISE_SD = 10.0            # counts
THRESHOLD = 30.0           # values below this are not stored, as on the instrument
IS_NG = 10.0
rng = np.random.default_rng(SEED)
times = np.arange(int(RUN_END * RATE)) / RATE


# ------------------------------------------------------------------ retention
def alkane_t1(n):
    """First-dimension time (s) of a C`n` n-alkane; the spacing closes up through the oven ramp."""
    n = np.asarray(n, dtype=float)
    return 330.0 + 118.0 * (n - 9) - 0.8 * (n - 9) ** 2


def alkane_t2(t1):
    """Second-dimension time (s) of the n-alkane line at `t1`."""
    return 1.05 + 0.25 * (np.asarray(t1) - 330.0) / 2450.0


# -------------------------------------------------------------------- spectra
def with_isotopes(peaks):
    """Add 13C isotope ions (about 1.1% per carbon, carbons guessed as m/z / 13)."""
    out = dict(peaks)
    for mz, rel in peaks.items():
        c = mz / 13.0
        out[mz + 1] = out.get(mz + 1, 0.0) + rel * 0.011 * c
        out[mz + 2] = out.get(mz + 2, 0.0) + rel * (0.011 * c) ** 2 / 2
    return out


def alkane(n):
    r71 = {9: 0.50, 10: 0.54, 11: 0.58, 12: 0.62}.get(n, 0.66)
    r85 = {9: 0.33, 10: 0.36, 11: 0.39, 12: 0.42}.get(n, 0.44)
    s = {41: 40, 43: 90 if n <= 12 else 75, 55: 25, 56: 8, 57: 100, 69: 12, 70: 6, 71: 100 * r71,
         83: 8, 84: 4, 85: 100 * r85, 97: 5, 111: 3}
    rel, mz = 100 * r85 * 0.62, 99
    while mz <= 14 * n + 2 - 29:
        s[mz] = max(rel, 2.0)
        rel *= 0.72
        mz += 14
    s[14 * n + 2] = 4.0 if n < 20 else 1.5
    return s


def methylalkane(n):
    m = 14 * n + 2
    return {43: 100, 41: 40, 55: 22, 57: 75, 69: 10, 71: 30, 85: 18, 99: 10, 113: 6, m - 43: 12, m - 15: 3}


def cyclohexane(n):
    return {41: 45, 43: 30, 55: 65, 56: 15, 57: 15, 67: 20, 69: 25, 82: 75, 83: 100, 96: 8, 97: 12,
            14 * n: 8}


def alkene(n):
    return {41: 100, 43: 75, 55: 95, 56: 40, 57: 45, 69: 60, 70: 35, 83: 45, 97: 30, 111: 15, 125: 8,
            14 * n: 2}


def alkylbenzene(chain):
    return {51: 5, 65: 10, 77: 6, 91: 100, 92: 55, 105: 7, 78 + 14 * chain: 20}


def trimethylbenzene():
    return {51: 6, 65: 4, 77: 12, 79: 8, 91: 12, 105: 100, 119: 15, 120: 55}


SPECTRA = {
    "dichloromethane": {47: 12, 48: 10, 49: 100, 51: 32, 84: 75, 86: 48, 88: 7},
    "toluene": {39: 10, 51: 6, 62: 4, 63: 8, 65: 12, 91: 100, 92: 70},
    "ethylbenzene": {51: 7, 65: 8, 77: 8, 78: 6, 91: 100, 92: 8, 105: 5, 106: 35},
    "xylene": {51: 10, 65: 6, 77: 12, 79: 8, 91: 100, 105: 25, 106: 65},
    "propylbenzene": {51: 6, 65: 10, 78: 6, 91: 100, 92: 10, 105: 5, 120: 25},
    "alpha-pinene": {39: 15, 41: 20, 77: 30, 79: 25, 91: 45, 92: 35, 93: 100, 105: 8, 121: 12, 136: 10},
    "limonene": {39: 20, 41: 30, 67: 70, 68: 100, 79: 30, 93: 80, 94: 20, 107: 20, 121: 25, 136: 25},
    "phenol": {39: 25, 40: 10, 55: 6, 65: 25, 66: 30, 94: 100},
    "o-cresol": {77: 25, 79: 30, 80: 15, 90: 20, 107: 90, 108: 100},
    "p-cresol": {39: 10, 51: 10, 77: 25, 79: 25, 90: 8, 107: 100, 108: 90},
    "guaiacol": {39: 10, 52: 10, 53: 15, 81: 55, 109: 100, 110: 8, 124: 90},
    "4-methylguaiacol": {67: 15, 77: 15, 95: 30, 123: 90, 138: 100},
    "syringol": {51: 10, 65: 15, 93: 20, 96: 30, 111: 25, 139: 60, 154: 100},
    "indane": {39: 10, 51: 8, 58: 6, 63: 8, 89: 8, 91: 25, 115: 30, 116: 15, 117: 100, 118: 75},
    "tetralin": {51: 8, 78: 10, 91: 45, 104: 100, 115: 15, 117: 25, 131: 10, 132: 70},
    "trans-decalin": {39: 20, 41: 45, 54: 30, 55: 40, 67: 85, 68: 60, 81: 80, 82: 70, 95: 70, 96: 100,
                      109: 15, 138: 45},
    "quinoline": {50: 6, 51: 10, 75: 8, 76: 10, 102: 20, 128: 15, 129: 100},
    "indole": {39: 10, 62: 8, 63: 12, 89: 35, 90: 45, 116: 15, 117: 100},
    "carbazole": {63: 6, 83: 10, 139: 15, 140: 10, 166: 20, 167: 100},
    "nicotine": {42: 25, 84: 100, 92: 8, 105: 10, 119: 10, 133: 30, 161: 25, 162: 20},
    "benzothiazole": {63: 8, 69: 12, 82: 6, 91: 5, 108: 35, 135: 100, 137: 5},
    "N,N-dimethyldodecylamine": {42: 8, 44: 6, 58: 100, 59: 5, 71: 2, 213: 3},
    "1-octanol": {41: 100, 42: 60, 43: 85, 55: 80, 56: 75, 57: 30, 69: 50, 70: 50, 83: 20, 84: 30},
    "1-dodecanol": {41: 85, 43: 90, 55: 100, 56: 50, 57: 35, 69: 70, 70: 45, 83: 55, 97: 30, 111: 10,
                    140: 5},
    "nonanal": {41: 95, 43: 85, 44: 40, 55: 50, 56: 70, 57: 100, 69: 25, 70: 55, 82: 30, 98: 40, 114: 6},
    "decanal": {41: 100, 43: 90, 55: 60, 56: 60, 57: 95, 70: 55, 82: 35, 95: 15, 112: 25, 128: 4},
    "2-undecanone": {41: 25, 43: 95, 57: 20, 58: 100, 59: 55, 71: 30, 85: 10, 112: 8, 127: 6, 170: 4},
    "2-octadecanone": {41: 20, 43: 80, 57: 25, 58: 100, 59: 45, 71: 30, 85: 12, 96: 10, 268: 4},
    "phytone": {43: 95, 57: 25, 58: 100, 59: 35, 71: 30, 85: 15, 109: 10, 124: 8, 250: 2},
    "naphthalene": {51: 6, 63: 4, 64: 5, 77: 3, 102: 7, 126: 6, 127: 12, 128: 100},
    "2-methylnaphthalene": {63: 5, 71: 6, 89: 5, 115: 30, 139: 10, 141: 85, 142: 100},
    "1-methylnaphthalene": {63: 5, 71: 6, 89: 5, 115: 28, 139: 10, 141: 80, 142: 100},
    "dimethylnaphthalene": {115: 25, 128: 8, 141: 75, 153: 10, 155: 30, 156: 100},
    "trimethylnaphthalene": {115: 10, 128: 10, 141: 10, 153: 20, 155: 90, 169: 25, 170: 100},
    "acenaphthylene": {63: 5, 74: 6, 76: 12, 126: 6, 150: 15, 151: 20, 152: 100},
    "acenaphthene": {63: 5, 76: 15, 77: 10, 151: 20, 152: 50, 153: 100, 154: 95},
    "dibenzofuran": {63: 8, 84: 8, 113: 6, 139: 40, 140: 6, 168: 100},
    "fluorene": {63: 5, 82: 10, 83: 8, 139: 8, 163: 12, 164: 10, 165: 95, 166: 100},
    "9-fluorenone": {75: 8, 76: 15, 126: 6, 150: 15, 151: 25, 152: 70, 180: 100},
    "benzophenone": {51: 20, 77: 60, 105: 100, 106: 8, 152: 6, 181: 25, 182: 55},
    "phenanthrene": {76: 8, 89: 12, 151: 6, 152: 10, 176: 20, 177: 10, 178: 100},
    "anthracene": {76: 8, 89: 11, 151: 6, 152: 9, 176: 18, 177: 9, 178: 100},
    "anthraquinone": {76: 30, 126: 8, 150: 15, 151: 25, 152: 70, 180: 85, 208: 100},
    "fluoranthene": {88: 6, 100: 10, 101: 15, 174: 5, 200: 20, 201: 15, 202: 100},
    "pyrene": {88: 6, 100: 12, 101: 18, 174: 5, 200: 22, 201: 14, 202: 100},
    "retene": {101: 6, 189: 10, 202: 15, 203: 25, 204: 35, 219: 100, 234: 70},
    "1-nitropyrene": {100: 15, 174: 10, 189: 35, 200: 40, 201: 65, 202: 20, 217: 40, 247: 100},
    "benz[a]anthracene": {101: 5, 113: 10, 114: 12, 200: 4, 226: 25, 227: 10, 228: 100},
    "chrysene": {101: 5, 113: 9, 114: 13, 200: 4, 226: 28, 227: 9, 228: 100},
    "benzo[b]fluoranthene": {113: 6, 125: 10, 126: 15, 224: 3, 250: 25, 251: 10, 252: 100},
    "benzo[a]pyrene": {113: 6, 125: 11, 126: 16, 224: 3, 250: 24, 251: 9, 252: 100},
    "perylene": {113: 6, 125: 10, 126: 14, 224: 3, 250: 26, 251: 10, 252: 100},
    "palmitic acid": {41: 60, 43: 85, 55: 45, 57: 50, 60: 100, 69: 25, 71: 20, 73: 95, 83: 15, 129: 35,
                      157: 8, 185: 10, 213: 15, 256: 25},
    "stearic acid": {43: 85, 55: 40, 57: 50, 60: 100, 73: 90, 129: 35, 185: 9, 241: 12, 284: 25},
    "behenic acid": {43: 80, 55: 40, 57: 55, 60: 90, 73: 100, 129: 40, 185: 10, 297: 12, 340: 20},
    "dehydroabietic acid": {141: 12, 155: 10, 197: 15, 239: 100, 240: 18, 285: 25, 300: 10},
    "levoglucosan": {43: 30, 57: 45, 60: 100, 69: 15, 70: 20, 71: 15, 73: 50, 85: 8, 98: 15},
    "methyl palmitate": {41: 35, 43: 50, 55: 35, 74: 100, 75: 20, 87: 70, 129: 8, 143: 15, 227: 8,
                         239: 5, 270: 8},
    "methyl stearate": {43: 50, 55: 30, 74: 100, 87: 70, 143: 15, 255: 6, 298: 8},
    "methyl oleate": {41: 70, 55: 100, 69: 70, 74: 45, 83: 50, 87: 35, 97: 40, 222: 10, 264: 25, 296: 5},
    "2,4-di-tert-butylphenol": {41: 10, 57: 30, 91: 5, 163: 5, 191: 100, 192: 15, 206: 30},
    "BHT": {41: 10, 57: 30, 91: 6, 145: 12, 177: 8, 205: 100, 206: 15, 220: 25},
    "galaxolide": {57: 15, 171: 8, 213: 40, 243: 100, 244: 20, 258: 25},
    "diethyl phthalate": {65: 8, 76: 6, 105: 8, 121: 6, 149: 100, 150: 10, 176: 8, 177: 25, 222: 3},
    "dibutyl phthalate": {41: 10, 57: 6, 76: 5, 104: 6, 121: 4, 149: 100, 150: 10, 205: 5, 223: 6},
    "DEHP": {43: 20, 57: 30, 70: 8, 71: 20, 83: 6, 113: 10, 149: 100, 150: 10, 167: 35, 279: 12},
    "pristane": {41: 35, 43: 85, 55: 18, 57: 100, 69: 10, 71: 62, 85: 35, 99: 12, 113: 20, 127: 8,
                 141: 5, 183: 9, 253: 1},
    "phytane": {41: 35, 43: 90, 55: 20, 57: 100, 69: 12, 71: 65, 85: 40, 99: 15, 113: 18, 127: 7, 183: 8},
    "D4": {73: 25, 133: 6, 193: 10, 207: 15, 249: 8, 265: 10, 281: 100},
    "D5": {73: 65, 147: 15, 251: 10, 267: 45, 281: 10, 355: 100},
    "D6": {73: 100, 147: 45, 267: 10, 325: 15, 341: 85},
    "D7": {73: 100, 147: 60, 221: 20, 281: 55, 355: 35},
    "D8": {73: 100, 147: 50, 207: 15, 221: 20, 281: 45, 355: 30},
    # Deuterated internal standards: the quantifier and qualifier ratios are those of
    # core.calbuilder.ISTD_SIGNATURES, measured on the real 2026 runs.
    "acenaphthene-d10": {80: 15, 158: 10, 160: 49.2, 162: 96.3, 164: 100},
    "phenanthrene-d10": {80: 8, 94: 6, 160: 13, 184: 18.6, 188: 100},
    "n-tetracosane-d50": {46: 25, 50: 40, 66: 100, 82: 72.3, 98: 51.2, 114: 32, 130: 20, 146: 12, 162: 8},
    "chrysene-d12": {120: 15, 236: 27.7, 238: 12, 240: 100},
    "perylene-d12": {132: 18.7, 260: 19.9, 262: 10, 264: 100},
}
ALKYL = {2: "ethyl", 3: "propyl", 4: "butyl", 5: "pentyl", 6: "hexyl", 7: "heptyl", 8: "octyl", 9: "nonyl",
         10: "decyl", 11: "undecyl", 12: "dodecyl", 13: "tridecyl", 14: "tetradecyl"}
ALKENE = {12: "1-dodecene", 14: "1-tetradecene", 16: "1-hexadecene", 18: "1-octadecene"}
for n in range(9, 35):
    SPECTRA[f"n-C{n}"] = alkane(n)
for n in range(10, 21):
    SPECTRA[f"n-{ALKYL[n - 6]}cyclohexane"] = cyclohexane(n)
for n, name in ALKENE.items():
    SPECTRA[name] = alkene(n)
for chain in range(4, 11):
    SPECTRA[f"n-{ALKYL[chain]}benzene"] = alkylbenzene(chain)
SPECTRA["1,2,4-trimethylbenzene"] = trimethylbenzene()
SPECTRA["1,3,5-trimethylbenzene"] = trimethylbenzene()

# Compound table: class, retention index, offset above the n-alkane line in t2 (s),
# and response per ng relative to an n-alkane. Isotopes are added when the spectrum is used.
COMPOUNDS = {
    "dichloromethane": ("solvent", 560, 0.9, 1.0),
    "toluene": ("alkylbenzene", 760, 0.9, 1.1),
    "ethylbenzene": ("alkylbenzene", 855, 0.95, 1.1),
    "m/p-xylene": ("alkylbenzene", 865, 0.95, 1.1),
    "o-xylene": ("alkylbenzene", 890, 1.05, 1.1),
    "propylbenzene": ("alkylbenzene", 952, 0.85, 1.1),
    "alpha-pinene": ("terpene", 935, 0.25, 0.9),
    "1,3,5-trimethylbenzene": ("alkylbenzene", 968, 0.95, 1.1),
    "phenol": ("aromatic oxygenate", 980, 2.2, 0.8),
    "1,2,4-trimethylbenzene": ("alkylbenzene", 990, 1.05, 1.1),
    "limonene": ("terpene", 1030, 0.45, 0.9),
    "indane": ("indane", 1040, 1.3, 1.1),
    "trans-decalin": ("cycloalkane", 1057, 0.45, 0.9),
    "quinoline": ("N-heteroaromatic", 1060, 2.6, 1.0),
    "1-octanol": ("alcohol", 1070, 0.75, 0.6),
    "p-cresol": ("aromatic oxygenate", 1075, 2.2, 0.8),
    "o-cresol": ("aromatic oxygenate", 1058, 2.0, 0.8),
    "guaiacol": ("aromatic oxygenate", 1090, 1.8, 0.8),
    "nonanal": ("aldehyde", 1104, 0.55, 0.6),
    "tetralin": ("indane", 1160, 1.45, 1.1),
    "naphthalene": ("naphthalene", 1185, 1.6, 1.3),
    "4-methylguaiacol": ("aromatic oxygenate", 1195, 1.8, 0.8),
    "decanal": ("aldehyde", 1206, 0.5, 0.6),
    "benzothiazole": ("S-heteroaromatic", 1225, 2.5, 0.9),
    "indole": ("N-heteroaromatic", 1290, 2.8, 0.9),
    "2-undecanone": ("ketone", 1293, 0.45, 0.7),
    "2-methylnaphthalene": ("naphthalene", 1300, 1.55, 1.3),
    "1-methylnaphthalene": ("naphthalene", 1315, 1.7, 1.3),
    "syringol": ("aromatic oxygenate", 1350, 2.0, 0.8),
    "2,6-dimethylnaphthalene": ("naphthalene", 1420, 1.5, 1.2),
    "1,3-dimethylnaphthalene": ("naphthalene", 1450, 1.6, 1.2),
    "acenaphthylene": ("PAH", 1460, 2.3, 1.3),
    "acenaphthene-d10": ("internal standard", 1465, 2.08, 1.2),
    "dodecanol": ("alcohol", 1470, 0.7, 0.6),
    "acenaphthene": ("PAH", 1480, 2.2, 1.2),
    "N,N-dimethyldodecylamine": ("amine", 1490, 0.15, 0.7),
    "levoglucosan": ("anhydrosugar", 1495, 3.2, 0.5),
    "BHT": ("aromatic oxygenate", 1515, 0.35, 0.9),
    "2,4-di-tert-butylphenol": ("aromatic oxygenate", 1520, 0.6, 0.9),
    "dibenzofuran": ("aromatic oxygenate", 1520, 2.4, 1.2),
    "2,3,5-trimethylnaphthalene": ("naphthalene", 1545, 1.5, 1.2),
    "fluorene": ("PAH", 1590, 2.5, 1.2),
    "diethyl phthalate": ("phthalate", 1595, 1.7, 0.9),
    "benzophenone": ("aromatic oxygenate", 1630, 2.4, 1.0),
    "pristane": ("branched alkane", 1707, -0.03, 1.0),
    "9-fluorenone": ("oxy-PAH", 1750, 3.0, 1.0),
    "phenanthrene-d10": ("internal standard", 1785, 2.80, 1.3),
    "phenanthrene": ("PAH", 1800, 2.9, 1.3),
    "anthracene": ("PAH", 1818, 3.15, 1.3),
    "phytane": ("branched alkane", 1811, -0.03, 1.0),
    "galaxolide": ("fragrance", 1840, 1.1, 0.9),
    "phytone": ("ketone", 1845, 0.35, 0.7),
    "nicotine": ("N-heteroaromatic", 1360, 1.6, 0.8),
    "methyl palmitate": ("ester", 1925, 0.35, 0.8),
    "palmitic acid": ("carboxylic acid", 1960, 0.7, 0.5),
    "dibutyl phthalate": ("phthalate", 1960, 1.4, 0.9),
    "2-octadecanone": ("ketone", 1990, 0.35, 0.7),
    "anthraquinone": ("oxy-PAH", 1990, 3.4, 0.9),
    "carbazole": ("N-heteroaromatic", 2050, 3.8, 0.9),
    "fluoranthene": ("PAH", 2080, 3.6, 1.3),
    "methyl oleate": ("ester", 2100, 0.45, 0.8),
    "pyrene": ("PAH", 2125, 3.9, 1.3),
    "methyl stearate": ("ester", 2125, 0.33, 0.8),
    "stearic acid": ("carboxylic acid", 2165, 0.65, 0.5),
    "retene": ("PAH", 2230, 2.9, 1.1),
    "dehydroabietic acid": ("carboxylic acid", 2350, 1.6, 0.5),
    "n-tetracosane-d50": ("internal standard", 2360, -0.04, 1.0),
    "chrysene-d12": ("internal standard", 2450, 4.30, 1.3),
    "benz[a]anthracene": ("PAH", 2470, 4.6, 1.3),
    "chrysene": ("PAH", 2480, 4.65, 1.3),
    "1-nitropyrene": ("nitro-PAH", 2500, 5.0, 0.8),
    "DEHP": ("phthalate", 2550, 0.9, 0.9),
    "behenic acid": ("carboxylic acid", 2560, 0.6, 0.5),
    "benzo[b]fluoranthene": ("PAH", 2850, 5.4, 1.3),
    "benzo[a]pyrene": ("PAH", 2920, 6.0, 1.3),
    "perylene-d12": ("internal standard", 2935, 7.56, 1.3),
    "perylene": ("PAH", 2950, 7.69, 1.3),
    "D4": ("siloxane", 990, -0.42, 0.8),
    "D5": ("siloxane", 1160, -0.40, 0.8),
    "D6": ("siloxane", 1330, -0.38, 0.8),
    "D7": ("siloxane", 1500, -0.36, 0.8),
    "D8": ("siloxane", 1660, -0.34, 0.8),
}
for n in range(9, 35):
    COMPOUNDS[f"n-C{n}"] = ("n-alkane", 100 * n, 0.0, 1.0)
for n in range(10, 21):
    COMPOUNDS[f"n-{ALKYL[n - 6]}cyclohexane"] = ("cycloalkane", 100 * n - 12, 0.3, 0.9)
for n, name in ALKENE.items():
    COMPOUNDS[name] = ("alkene", 100 * n - 9, 0.06, 0.9)
for chain in range(4, 11):
    COMPOUNDS[f"n-{ALKYL[chain]}benzene"] = ("alkylbenzene", 100 * (chain + 6) + 55, 0.85 - 0.05 * chain, 1.0)
SPECTRUM_OF = {"m/p-xylene": "xylene", "o-xylene": "xylene", "dodecanol": "1-dodecanol",
               "1,3-dimethylnaphthalene": "dimethylnaphthalene", "2,6-dimethylnaphthalene": "dimethylnaphthalene",
               "2,3,5-trimethylnaphthalene": "trimethylnaphthalene"}


def spectrum_vector(name):
    peaks = with_isotopes(SPECTRA[SPECTRUM_OF.get(name, name)])
    v = np.zeros(MZ.size)
    for mz, rel in peaks.items():
        if MZ[0] <= mz <= MZ[-1]:
            v[mz - MZ[0]] += rel
    return v / v.sum()


def position(name):
    """Nominal (t1, t2) of a compound from its retention index."""
    _, ri, dt2, _ = COMPOUNDS[name]
    t1 = float(alkane_t1(ri / 100.0))
    return t1, float(alkane_t2(t1) + dt2)


# ------------------------------------------------------------ what each run holds
LEVELS = {"L1": 1.0, "L2": 5.0, "L3": 20.0, "L4": 50.0}            # ng on tube
# Stock factors: decalin and cycloalkanes from a 5x more concentrated stock, oxy- and
# nitro-PAHs from a 10x weaker one, mirroring the August 2026 standard scheme.
CAL_MIX = {**{f"n-C{n}": 1.0 for n in range(10, 31, 2)},
           "trans-decalin": 5.0, "n-hexylcyclohexane": 5.0, "n-decylcyclohexane": 5.0,
           "1,2,4-trimethylbenzene": 1.0, "n-butylbenzene": 1.0,
           "naphthalene": 1.0, "2-methylnaphthalene": 1.0, "acenaphthylene": 1.0, "acenaphthene": 1.0,
           "fluorene": 1.0, "phenanthrene": 1.0, "anthracene": 1.0, "fluoranthene": 1.0, "pyrene": 1.0,
           "benz[a]anthracene": 1.0, "chrysene": 1.0, "benzo[b]fluoranthene": 1.0, "benzo[a]pyrene": 1.0,
           "9-fluorenone": 0.1, "anthraquinone": 0.1, "1-nitropyrene": 0.1, "2-octadecanone": 1.0,
           "behenic acid": 5.0}
ISTDS = ["acenaphthene-d10", "phenanthrene-d10", "n-tetracosane-d50", "chrysene-d12", "perylene-d12"]
TUBE_BACKGROUND = {"D4": 0.6, "D5": 1.0, "D6": 0.7, "D7": 0.4, "D8": 0.3, "dibutyl phthalate": 0.8,
                   "DEHP": 1.5, "BHT": 0.4, "2,4-di-tert-butylphenol": 0.5}

urban = {**{f"n-C{n}": round(12.0 * np.exp(-0.5 * ((n - 19) / 5.0) ** 2) + 0.5, 2) for n in range(9, 32)},
         "toluene": 900.0, "ethylbenzene": 25.0, "m/p-xylene": 60.0, "o-xylene": 22.0, "propylbenzene": 6.0,
         "1,3,5-trimethylbenzene": 8.0, "1,2,4-trimethylbenzene": 18.0, "n-butylbenzene": 3.0,
         "alpha-pinene": 6.0, "limonene": 14.0, "nonanal": 9.0, "decanal": 7.0, "naphthalene": 12.0,
         "2-methylnaphthalene": 5.0, "1-methylnaphthalene": 3.0, "acenaphthylene": 1.2, "acenaphthene": 0.8,
         "fluorene": 1.5, "phenanthrene": 3.5, "anthracene": 0.6, "fluoranthene": 1.8, "pyrene": 1.6,
         "benz[a]anthracene": 0.5, "chrysene": 0.7, "benzo[b]fluoranthene": 0.6, "benzo[a]pyrene": 0.4,
         "9-fluorenone": 0.6, "anthraquinone": 0.4, "1-nitropyrene": 0.12, "D4": 4.0, "D5": 40.0,
         "D6": 6.0, "D7": 1.0, "D8": 0.5, "diethyl phthalate": 3.0, "dibutyl phthalate": 4.0, "DEHP": 9.0,
         "BHT": 2.0, "2,4-di-tert-butylphenol": 3.0, "benzophenone": 2.0, "benzothiazole": 1.5,
         "galaxolide": 2.5, "palmitic acid": 8.0, "stearic acid": 4.0, "methyl palmitate": 1.0,
         "2-octadecanone": 1.2, "phytone": 2.0, "trans-decalin": 3.0, "indane": 2.0, "1-octanol": 1.5,
         "N,N-dimethyldodecylamine": 1.0, "nicotine": 2.0, "behenic acid": 2.0}
diesel = {**{f"n-C{n}": round(35.0 * np.exp(-0.5 * ((n - 18) / 4.0) ** 2) + 0.3, 2) for n in range(9, 31)},
          "pristane": 14.0, "phytane": 10.0,
          **{f"n-{ALKYL[n - 6]}cyclohexane": round(8.0 * np.exp(-0.5 * ((n - 16) / 3.5) ** 2), 2) for n in range(10, 21)},
          "trans-decalin": 12.0, "tetralin": 3.0, "indane": 2.5, "toluene": 20.0, "ethylbenzene": 6.0,
          "m/p-xylene": 14.0, "o-xylene": 6.0, "propylbenzene": 4.0, "1,3,5-trimethylbenzene": 6.0,
          "1,2,4-trimethylbenzene": 14.0,
          **{f"n-{ALKYL[c]}benzene": round(4.0 * np.exp(-0.3 * (c - 4)), 2) for c in range(4, 11)},
          "naphthalene": 18.0, "2-methylnaphthalene": 22.0, "1-methylnaphthalene": 12.0,
          "2,6-dimethylnaphthalene": 14.0, "1,3-dimethylnaphthalene": 10.0, "2,3,5-trimethylnaphthalene": 8.0,
          "acenaphthylene": 0.8, "acenaphthene": 2.0, "fluorene": 5.0, "dibenzofuran": 1.0,
          "phenanthrene": 9.0, "anthracene": 1.2, "fluoranthene": 2.5, "pyrene": 3.5,
          "benz[a]anthracene": 0.6, "chrysene": 0.9, "benzo[b]fluoranthene": 0.3, "benzo[a]pyrene": 0.2,
          "9-fluorenone": 1.5, "anthraquinone": 0.8, "1-nitropyrene": 0.35, "DEHP": 1.0, "D5": 0.8,
          "2-octadecanone": 0.5, "carbazole": 0.4, "quinoline": 0.5}
woodsmoke = {**{f"n-C{n}": round((6.0 if n % 2 else 2.0) * np.exp(-0.5 * ((n - 27) / 4.0) ** 2) + 0.3, 2)
                for n in range(15, 34)},
             "phenol": 20.0, "o-cresol": 9.0, "p-cresol": 12.0, "guaiacol": 25.0, "4-methylguaiacol": 15.0,
             "syringol": 30.0, "levoglucosan": 120.0, "retene": 6.0, "dehydroabietic acid": 10.0,
             "palmitic acid": 15.0, "stearic acid": 7.0, "behenic acid": 4.0, "methyl palmitate": 2.0,
             "methyl stearate": 1.5, "methyl oleate": 2.5,
             **{ALKENE[n]: v for n, v in ((12, 2.0), (14, 2.5), (16, 3.0), (18, 2.0))},
             "toluene": 30.0, "ethylbenzene": 4.0, "m/p-xylene": 6.0, "naphthalene": 25.0,
             "2-methylnaphthalene": 6.0, "1-methylnaphthalene": 4.0, "acenaphthylene": 6.0, "acenaphthene": 1.0,
             "fluorene": 3.0, "phenanthrene": 10.0, "anthracene": 2.0, "fluoranthene": 6.0, "pyrene": 5.5,
             "benz[a]anthracene": 2.0, "chrysene": 2.2, "benzo[b]fluoranthene": 1.8, "benzo[a]pyrene": 1.4,
             "perylene": 0.5, "9-fluorenone": 1.2, "anthraquinone": 0.7, "dibenzofuran": 2.0, "indole": 1.5,
             "quinoline": 1.0, "carbazole": 0.8, "2-undecanone": 1.5, "2-octadecanone": 1.0, "1-octanol": 2.0,
             "dodecanol": 1.5, "nonanal": 4.0, "limonene": 3.0, "alpha-pinene": 5.0, "D5": 2.0, "DEHP": 1.5}

RUNS = {
    # name: (type, level, sample volume m3, contents, t1 drift s, t2 drift s)
    "CAL_L1_1ng": ("calibration", "L1", None, None, 0.0, 0.0),
    "CAL_L2_5ng": ("calibration", "L2", None, None, 0.6, 0.004),
    "CAL_L3_20ng": ("calibration", "L3", None, None, -0.5, -0.004),
    "CAL_L4_50ng": ("calibration", "L4", None, None, 0.9, 0.008),
    "SAMPLE_1_urban_air": ("sample", "", 0.005, urban, 2.5, 0.02),
    "SAMPLE_2_diesel_exhaust": ("sample", "", 0.002, diesel, 3.6, 0.03),
    "SAMPLE_3_wood_smoke": ("sample", "", 0.010, woodsmoke, -2.2, -0.015),
}


# --------------------------------------------------------------- building a run
def add_peak(dense, phase, t1c, t2c, area, spectrum, sigma1, sigma2, tau1=0.0, tau2=0.0):
    """Add one compound at scan level: each modulation k releases the share of the
    first-dimension peak that eluted during it, as a second-dimension pulse at t2c."""
    lo, hi = t1c - 5 * sigma1, t1c + 5 * sigma1 + 8 * tau1
    grid = np.arange(lo, hi, 0.05)
    g1 = (exponnorm.pdf(grid, tau1 / sigma1, loc=t1c, scale=sigma1) if tau1 > 0
          else norm.pdf(grid, t1c, sigma1))
    g1 /= g1.sum()
    ions = np.flatnonzero(spectrum > 0)
    for k in range(int(np.floor((lo - phase) / PERIOD)), int(np.ceil((hi - phase) / PERIOD)) + 1):
        start = phase + k * PERIOD
        frac = g1[(grid >= start) & (grid < start + PERIOD)].sum()
        if frac < 1e-4:
            continue
        mu = start + t2c
        i0 = max(int((mu - 6 * sigma2) * RATE), 0)
        i1 = min(int((mu + 6 * sigma2 + 10 * tau2) * RATE) + 1, times.size)
        if i1 <= i0:
            continue
        t = times[i0:i1]
        g2 = (exponnorm.pdf(t, tau2 / sigma2, loc=mu, scale=sigma2) if tau2 > 0
              else norm.pdf(t, mu, sigma2)) / RATE
        dense[i0:i1, ions] += (area * frac) * np.outer(g2, spectrum[ions])


def sigmas(t2):
    return 3.6, 0.042 + 0.011 * t2


truth_rows, amount_rows, run_info = [], [], {}
(OUT / "data").mkdir(parents=True, exist_ok=True)
for run_name, (kind, level, volume, contents, drift1, drift2) in RUNS.items():
    phase = round(float(rng.uniform(1.6, 2.6)), 3)
    if kind == "calibration":
        contents = {c: LEVELS[level] * f for c, f in CAL_MIX.items()}
        contents.update({k: v for k, v in TUBE_BACKGROUND.items() if k not in contents})
    else:
        contents = {**{k: v for k, v in TUBE_BACKGROUND.items()}, **contents}
    contents.update({s: IS_NG for s in ISTDS})
    dense = np.zeros((times.size, MZ.size), dtype=np.float64)

    # Column bleed: from the first column it is modulated into a band at the bottom
    # of each modulation; from the second it is a flat floor. Both rise with the oven.
    mod_start = phase + np.floor((times - phase) / PERIOD) * PERIOD
    t2_scan = times - mod_start
    bleed = np.zeros(MZ.size)
    for mz, rel in BLEED_REFERENCE.items():
        bleed[mz - MZ[0]] = rel / 100.0
    band = 1.0 / (1.0 + np.exp(-(t2_scan - 0.30) / 0.04)) * np.exp(-np.clip(t2_scan - 0.30, 0, None) / 1.4)
    ramp = np.clip(mod_start / RUN_END, 0.0, 1.0)
    first = 150.0 + 1500.0 * ramp ** 3
    second = 30.0 + 900.0 * ramp ** 4
    dense += np.outer(first * band + second, bleed)
    # Air and a pressure pulse from the modulator at the start of every modulation (a t2 stripe).
    stripe = np.exp(-0.5 * ((t2_scan - 0.12) / 0.03) ** 2)
    dense[:, 40 - MZ[0]] += 700.0 + 500.0 * stripe
    dense[:, 44 - MZ[0]] += 350.0 + 300.0 * stripe
    dense[:, 39 - MZ[0]] += 45.0

    # Unresolved hump of the diesel sample, continuous in t1, centred near C18.
    if run_name.startswith("SAMPLE_2"):
        n_scan = np.interp(mod_start, alkane_t1(np.arange(5, 40)), np.arange(5, 40))
        hump = 2.2e3 * np.exp(-0.5 * ((n_scan - 18.5) / 4.0) ** 2) * \
            np.exp(-0.5 * ((t2_scan - alkane_t2(mod_start) - 0.35) / 0.32) ** 2)
        ucm = np.zeros(MZ.size)
        for mz, rel in {41: 40, 43: 80, 55: 50, 57: 100, 67: 20, 69: 40, 71: 55, 81: 18, 83: 30, 85: 30,
                        95: 12, 97: 15, 99: 15, 111: 8, 113: 8}.items():
            ucm[mz - MZ[0]] = rel
        dense += np.outer(hump, ucm / ucm.sum() * 6.0)

    # Named compounds, each shifted by this run's retention drift.
    for name, ng in contents.items():
        cls, ri, dt2, response = COMPOUNDS[name]
        t1c, t2c = position(name)
        t1c += drift1 * t1c / 1500.0
        t2c += drift2
        s1, s2 = sigmas(t2c)
        tau1 = tau2 = 0.0
        if cls in ("carboxylic acid", "anhydrosugar"):
            tau1, tau2 = 3.0, 0.12          # polar compounds tail on both columns
        if name == "toluene" and ng > 100:
            s2, tau2 = 0.12, 1.2            # overloaded: smeared up the second dimension
        add_peak(dense, phase, t1c, t2c, ng * response * NG_TO_AREA, spectrum_vector(name), s1, s2, tau1, tau2)
        amount_rows.append({"run": run_name, "compound": name, "class": cls,
                            "amount_ng": ng, "t1_s": round(t1c, 2), "t2_s": round(t2c, 3)})
    if kind == "calibration":
        # The standards were spiked in dichloromethane, which floods the start of the run.
        add_peak(dense, phase, 175.0, 1.0, 4e3 * NG_TO_AREA, spectrum_vector("dichloromethane"), 9.0, 0.15,
                 tau2=1.5)

    # Minor isomers: branched alkanes everywhere, and hundreds more in the diesel sample.
    n_minor = {"SAMPLE_2": 700, "SAMPLE_1": 220, "SAMPLE_3": 160}.get(run_name[:8], 60)
    for _ in range(n_minor):
        n = float(np.clip(rng.normal(18.0, 4.0), 9.5, 30))
        kind_i = rng.choice(["methylalkane", "cyclohexane", "alkylbenzene", "naphthalene"],
                            p=[0.45, 0.3, 0.17, 0.08])
        t1c = float(alkane_t1(n))
        if kind_i == "methylalkane":
            peaks, dt2 = methylalkane(int(round(n)) + 1), rng.uniform(-0.08, 0.02)
        elif kind_i == "cyclohexane":
            peaks, dt2 = cyclohexane(int(round(n))), rng.uniform(0.15, 0.5)
        elif kind_i == "alkylbenzene":
            peaks, dt2 = alkylbenzene(max(int(round(n)) - 6, 2)), rng.uniform(0.4, 1.0)
        else:
            peaks, dt2 = SPECTRA["trimethylnaphthalene" if n > 15 else "dimethylnaphthalene"], rng.uniform(1.2, 1.8)
        v = np.zeros(MZ.size)
        for mz, rel in with_isotopes(peaks).items():
            if MZ[0] <= mz <= MZ[-1]:
                v[mz - MZ[0]] += rel
        ng = float(rng.lognormal(np.log(0.6 if run_name.startswith("SAMPLE_2") else 0.25), 0.8))
        add_peak(dense, phase, t1c + drift1 * t1c / 1500.0, float(alkane_t2(t1c) + dt2 + drift2),
                 ng * NG_TO_AREA, v / v.sum(), *sigmas(1.2))

    # Detector: bends over, saturates, adds counting and electronic noise, drops what
    # falls below the storage threshold, and records the odd single-scan spike.
    blocks = []
    for i0 in range(0, times.size, 12000):
        chunk = dense[i0:i0 + 12000]
        chunk = chunk / (1.0 + chunk / SOFT)
        chunk = chunk + rng.normal(0.0, 1.0, chunk.shape) * np.sqrt(NOISE_SD ** 2 + 2.0 * chunk)
        chunk = np.minimum(np.rint(chunk), CEILING)
        chunk[chunk < THRESHOLD] = 0.0
        blocks.append(sp.csr_matrix(chunk.astype(np.float32)))
    data = sp.vstack(blocks).tocsr()
    spikes = rng.integers(0, times.size, 80), rng.integers(0, MZ.size, 80)
    data = data + sp.csr_matrix((rng.uniform(2e4, 2e5, 80).astype(np.float32), spikes), shape=data.shape)
    run = Run(times=times.copy(), mz=MZ.copy(), data=sp.csr_matrix(data), name=run_name)
    write_andi(run, OUT / "data" / f"{run_name}.cdf",
               comment="Synthetic PyGCGC demonstration run (examples/demo_data/generate_demo_data.py)")
    run_info[run_name] = {"phase_s": phase, "t1_drift_s_at_25min": drift1, "t2_drift_s": drift2,
                          "type": kind, "level": level, "volume_m3": volume, "stored_points": int(data.nnz),
                          "spikes": 80}
    print(f"{run_name}: phase {phase:.3f} s, {data.nnz:,} stored points")

# ------------------------------------------------------------ supporting files
with open(OUT / "sequence.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["file", "type", "level", "cal_set", "volume_m3"])
    for name, (kind, level, volume, *_rest) in RUNS.items():
        w.writerow([f"{name}.cdf", kind, level, "", "" if volume is None else volume])

(OUT / "method.json").write_text(json.dumps({"modulation_period_s": PERIOD}, indent=2) + "\n")


def quant_ions(name, n=3):
    """The quantifier and two qualifier ions: the largest ions of the spectrum, the
    heavier one preferred when two are within 10% (as an analyst would pick)."""
    v = spectrum_vector(name)
    order = np.argsort(-(v * (1 + 0.001 * MZ)))
    quant = int(MZ[order[0]])
    quals = [int(MZ[i]) for i in order[1:n]]
    return quant, quals, [round(float(v[q - MZ[0]] / v[quant - MZ[0]]), 3) for q in quals]


# The alkanes are quantified on m/z 57 like most laboratories do; the deuterated
# tetracosane on 66, its own base peak.
target_rows = []
for name in ISTDS + list(CAL_MIX):
    t1, t2 = position(name)
    quant, quals, ratios = quant_ions(name)
    if name.startswith("n-C"):
        quant, quals = 57, [71, 85]
        v = spectrum_vector(name)
        ratios = [round(float(v[q - MZ[0]] / v[57 - MZ[0]]), 3) for q in quals]
    target_rows.append([name, "istd" if name in ISTDS else "target", quant,
                        ";".join(map(str, quals)), ";".join(map(str, ratios)), round(t1, 1), 15,
                        round(t2, 2), 0.25, ""])
with open(OUT / "target_table.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["compound", "role", "quant_ion", "qual_ions", "qual_ratios", "t1_expected_s", "t1_window_s",
                "t2_expected_s", "t2_window_s", "istd"])
    w.writerows(target_rows)

with open(OUT / "calibration_levels.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["level"] + [c for c in CAL_MIX] + ISTDS)
    for level, ng in LEVELS.items():
        w.writerow([level] + [f"{ng * f:g}" for f in CAL_MIX.values()] + [f"{IS_NG:g}"] * len(ISTDS))

# The same compounds in the layout the 📈 Calibrate window reads: one amount column per level.
with open(OUT / "calibration_compounds.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["compound", "t1_s", "t2_s", "role", "quant_ion", "istd"] + list(LEVELS))
    quant_of = {row[0]: row[2] for row in target_rows}
    for name in ISTDS + list(CAL_MIX):
        t1, t2 = position(name)
        amounts = [IS_NG] * 4 if name in ISTDS else [ng * CAL_MIX[name] for ng in LEVELS.values()]
        w.writerow([name, round(t1, 1), round(t2, 3), "istd" if name in ISTDS else "target",
                    quant_of[name], ""] + [f"{a:g}" for a in amounts])

# A small library for 🏷️ Identify: every named compound, its retention index, and a
# spectrum perturbed by 10% so that matches are good but not perfect, as with a real library.
with open(OUT / "demo_library.msp", "w", encoding="utf-8") as fh:
    for name, (cls, ri, _dt2, _r) in sorted(COMPOUNDS.items(), key=lambda kv: kv[1][1]):
        if cls == "solvent":
            continue
        v = spectrum_vector(name) * rng.lognormal(0.0, 0.1, MZ.size)
        v = v / v.max() * 999
        keep = np.flatnonzero(v >= 5)
        fh.write(f"Name: {name}\nRI: {ri}\nComment: synthetic; class {cls}\n"
                 f"Num Peaks: {keep.size}\n")
        fh.write("".join(f"{MZ[i]} {v[i]:.0f}\n" for i in keep) + "\n")

(OUT / "truth").mkdir(exist_ok=True)
with open(OUT / "truth" / "truth_amounts.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(amount_rows[0]))
    w.writeheader()
    w.writerows(amount_rows)
(OUT / "truth" / "truth_runs.json").write_text(json.dumps(
    {"period_s": PERIOD, "scan_rate_hz": RATE, "mz_range": [int(MZ[0]), int(MZ[-1])], "detector_ceiling": CEILING,
     "storage_threshold": THRESHOLD, "noise_sd": NOISE_SD, "seed": SEED, "runs": run_info}, indent=2) + "\n")
print(f"wrote {OUT}")
