"""A miniature NIST-style SDF library written on the fly, with invented spectra.

Real NIST data is licensed and never goes in the repository. These files copy
only the *format*: a structure block, labelled fields, CRLF line endings, an
odd `> <CASNO>` header spacing, several peak pairs on one line, and a record
with no CAS number."""

from __future__ import annotations

from pathlib import Path

MOL = ["  NIST    11042915402D 1   1.00000     0.00000", "invented structure block",
       "  1  0  0  0  0  0  0  0  0  0999 V2000",
       "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0", "M  END"]

# name, formula, mw, cas, nistno, peaks (list of (mz, intensity)), estimated RI, error
MAIN = [
    ("Decane", "C10H22", 142, "124185", 1001, [(43, 999), (57, 850), (71, 500), (85, 300), (142, 60)], 1000, 20),
    ("Undecane", "C11H24", 156, "1120214", 1002, [(43, 999), (57, 900), (71, 550), (85, 300), (156, 50)], 1100, 20),
    ("Toluene", "C7H8", 92, "108883", 1003, [(39, 80), (65, 100), (91, 999), (92, 600)], 770, 30),
    ("Propylbenzene", "C9H12", 120, "103651", 1004, [(91, 999), (92, 80), (105, 90), (120, 300)], 950, 25),
    ("Isopropylbenzene look-alike", "C9H12", 120, "98828", 1005, [(91, 999), (92, 80), (105, 95), (120, 290)], 1050, 25),
    ("Compound without CAS", "C5H10O", 86, "", 1006, [(41, 400), (44, 999), (58, 300), (86, 150)], None, None),
]
REP = [("Toluene", "C7H8", 92, "108883", 1003, [(39, 70), (65, 90), (91, 999), (92, 640)], 770, 30)]
# name, formula, cas, nistno, value, column class, data type, program type
RI = [
    ("Decane", "C10H22", "124185", 1001, 1000.0, "Semi-standard non-polar", "Normal alkane RI", "Ramp"),
    ("Decane", "C10H22", "124185", 1001, 1001.0, "Semi-standard non-polar", "Normal alkane RI", "Ramp"),
    ("Decane", "C10H22", "124185", 1001, 999.0, "Semi-standard non-polar", "Linear RI", "Complex"),
    ("Decane", "C10H22", "124185", 1001, 1000.0, "Standard polar", "Normal alkane RI", "Ramp"),
    ("Decane", "C10H22", "124185", 1001, 1040.0, "Semi-standard non-polar", "Kovats RI", "Isothermal"),
    ("Undecane", "C11H24", "1120214", 1002, 1100.0, "All column types",
     "Normal alkane RI value specified by scale definition", ""),               # n-alkanes: fixed by definition
    ("Toluene", "C7H8", "108883", -1, 765.0, "Semi-standard non-polar", "Normal alkane RI", "Ramp"),
    ("Toluene", "C7H8", "108883", -1, 771.0, "Semi-standard non-polar", "Normal alkane RI", "Ramp"),
    ("Propylbenzene", "C9H12", "103651", 1004, 950.0, "Semi-standard non-polar", "Normal alkane RI", "Ramp"),
    ("Isopropylbenzene look-alike", "C9H12", "98828", 1005, 1050.0, "Semi-standard non-polar", "Normal alkane RI", "Ramp"),
    ("COMPOUND without CAS!", "C5H10O", "", -1, 880.0, "Semi-standard non-polar", "Normal alkane RI", "Ramp"),
    ("Nothing in the library", "C2H6", "74840", 9999, 200.0, "Semi-standard non-polar", "Normal alkane RI", "Ramp"),
]


def _peaks_text(peaks, style: str) -> str:
    if style == "pairs":                       # several pairs on a line, as some exports do
        return "\r\n".join(" ".join(f"{m} {i}" for m, i in peaks[k:k + 2]) for k in range(0, len(peaks), 2))
    return "\r\n".join(f"{m} {i}" for m, i in peaks)


def write_ei(path: Path, rows, odd_header: bool = True) -> Path:
    lines = []
    for number, (name, formula, mw, cas, nistno, peaks, est, err) in enumerate(rows):
        lines += [name.upper()] + MOL
        lines += [">  <NAME>", name, "", ">  <SYNONYMS>", "a synonym", "another synonym", "", ">  <FORMULA>", formula, "",
                  ">  <MW>", str(mw), ""]
        if cas:
            lines += [("> <CASNO>" if odd_header and number == 1 else ">  <CASNO>"), cas, ""]
        lines += [">  <NISTNO>", str(nistno), "", ">  <NUM PEAKS>", str(len(peaks)), "", ">  <MASS SPECTRAL PEAKS>",
                  _peaks_text(peaks, "pairs" if number == 2 else "lines"), ""]
        if est is not None:
            lines += [">  <ESTIMATED KOVATS RI>", str(est), "", ">  <RI ESTIMATION ERROR>", str(err), ""]
        lines += ["$$$$"]
    path.write_bytes(("\r\n".join(lines) + "\r\n").encode("latin-1"))
    return path


def write_ri(path: Path, rows) -> Path:
    lines = []
    for name, formula, cas, nistno, value, klass, dtype, prog in rows:
        lines += [f"{name}, ID: C{cas}"] + MOL
        lines += ["", ">  <FORMULA>", formula, "", ">  <NAME>", name, ""]
        if cas:
            lines += ["> <CASNO>", cas, ""]
        if nistno >= 0:
            lines += [">  <NISTNO>", str(nistno), ""]
        lines += [">  <KOVATS INDEX>", f"{value}", "", ">  <COLUMN CLASS>", klass, "", ">  <DATA TYPE>", dtype, "",
                  ">  <PROGRAM TYPE>", prog, "", ">  <AUTHORS>", "Invented et al.", "", "$$$$"]
    path.write_bytes(("\r\n".join(lines) + "\r\n").encode("latin-1"))
    return path


def write_mini_nist(folder: Path) -> dict:
    """`mainlib.sdf`, `replib.sdf` and `nist_ri.sdf` in `folder / "SDF"`."""
    sdf = Path(folder) / "SDF"
    sdf.mkdir(parents=True, exist_ok=True)
    return {"mainlib": write_ei(sdf / "mainlib.sdf", MAIN), "replib": write_ei(sdf / "replib.sdf", REP),
            "nist_ri": write_ri(sdf / "nist_ri.sdf", RI)}
