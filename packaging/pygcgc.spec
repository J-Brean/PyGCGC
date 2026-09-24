# PyInstaller spec for the PyGCGC desktop app.
#
#   uv pip install -e ".[gui]" pyinstaller
#   pyinstaller packaging/pygcgc.spec
#
# Written for milestone 12 but not yet built or run anywhere -- there was no
# environment available to try it in during that milestone. Before trusting
# it: build it, run the frozen PyGCGC.exe (or PyGCGC on Linux/macOS) against
# the dummy data pack, and add whatever `hiddenimports` PyInstaller's own
# warnings (build/PyGCGC/warn-PyGCGC.txt) say are missing -- netCDF4's HDF5
# backend and scipy's compiled extensions are the most likely gaps, since
# they load dynamically rather than through a plain top-level import.
from pathlib import Path

root = Path(SPECPATH).parent  # noqa: F821 -- SPECPATH is injected by PyInstaller

a = Analysis(  # noqa: F821
    [str(root / "packaging" / "pyinstaller_entry.py")],
    pathex=[str(root / "src")],
    binaries=[],
    datas=[(str(root / "src" / "pygcgc" / "presets"), "pygcgc/presets")],   # chemical_classes.toml
    hiddenimports=[
        "netCDF4", "netCDF4.utils", "netCDF4._netCDF4", "cftime",
        "scipy.special._cdflib",
        "pyqtgraph.imageview",
        # scikit-learn's compiled helpers are imported dynamically (HDBSCAN, k-means, SVD in 🧩)
        "sklearn.cluster._hdbscan._tree", "sklearn.cluster._hdbscan._linkage", "sklearn.cluster._hdbscan._reachability",
        "sklearn.utils._typedefs", "sklearn.neighbors._partition_nodes",
    ],
    hookspath=[],
    runtime_hooks=[],
    excludes=[],
)
pyz = PYZ(a.pure, a.zipped_data)  # noqa: F821

exe = EXE(  # noqa: F821
    pyz, a.scripts, [],
    exclude_binaries=True,
    name="PyGCGC",
    debug=False,
    strip=False,
    upx=False,
    console=False,
)

coll = COLLECT(  # noqa: F821
    exe, a.binaries, a.zipfiles, a.datas,
    strip=False,
    upx=False,
    name="PyGCGC",
)
