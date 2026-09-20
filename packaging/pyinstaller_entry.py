"""PyInstaller's entry point: the same thing `pygcgc-gui` runs, as a plain
script PyInstaller can analyse directly (it doesn't follow console-script
entry points in pyproject.toml on its own)."""

import sys

from pygcgc.gui.app import main

if __name__ == "__main__":
    sys.exit(main())
