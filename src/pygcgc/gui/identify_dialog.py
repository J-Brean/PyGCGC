"""🏷️ Identify dialog.

Load an MSP library (NIST via `lib2nist`, MoNA, MassBank, or one built
here from accepted IDs -- "Add to in-house library…"), then search every
blob on the current run against it: a retention index from ⛓️'s own
alkane ladder if one's been found on this run, a weighted dot product
(forward and reverse) plus a peak-overlap ratio for the spectrum, and a
target-decoy q-value and a simplified confidence tier across the whole
batch of best hits. "Promote to target…" writes a confident hit out as a
target-table row, at the position it was actually found.

A one-shot, read-only analysis of whichever run is open when it's shown
-- like ✂️ overlap and 🗂️ batch, there's nothing here to undo or save
(see docs/architecture.md section 5), so it's a dialog, not a docked panel or its
own Project state. Reopening it (the toolbar button) is how it picks up
a run switch; the loaded library stays put across that, since it isn't
tied to any one run.
"""

from __future__ import annotations

from pathlib import Path

from PySide6 import QtCore, QtWidgets

import numpy as np
import pyqtgraph as pg

from ..core.identify import (IdentifySettings, LibraryEntry, build_library,
                             generate_decoys, identify_blobs, identify_components, load_msp, promote_to_target,
                             retention_index, write_msp)
from ..core.targets import write_targets_csv

COLUMNS = ["peak", "t1 (min)", "t2 (s)", "name", "score", "RI (query)", "RI diff", "tier", "q-value"]


LIB_COLUMNS = ["Name", "RI", "Formula", "CAS", "Peaks", "Decoy"]


class LibraryBrowser(QtWidgets.QDialog):
    """List, filter, edit and delete a loaded library's entries, see the
    selected spectrum, and save the real (non-decoy) entries as an MSP file.
    Name, RI, formula and CAS are editable in place. Closing the dialog
    leaves the edited library in `self.library` (rebuilt, so search sees it)."""

    def __init__(self, library, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Library entries")
        self.resize(820, 560)
        self.entries = list(library.entries)
        self.library = library
        self._updating = False
        self.filter = QtWidgets.QLineEdit()
        self.filter.setPlaceholderText("Filter by name, formula or CAS")
        self.filter.textChanged.connect(self._fill)
        self.show_decoys = QtWidgets.QCheckBox("Show decoys")
        self.show_decoys.toggled.connect(self._fill)
        self.table = QtWidgets.QTableWidget(0, len(LIB_COLUMNS))
        self.table.setHorizontalHeaderLabels(LIB_COLUMNS)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.table.horizontalHeader().setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        self.table.itemChanged.connect(self._edited)
        self.table.itemSelectionChanged.connect(self._draw)
        self.plot = pg.PlotWidget()
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.plot.setLabel("bottom", "m/z")
        self.plot.setFixedHeight(150)
        self.sticks = pg.PlotDataItem([], [], connect="pairs", pen=pg.mkPen((209, 149, 224), width=2))
        self.plot.addItem(self.sticks)
        delete = QtWidgets.QPushButton("Delete selected")
        delete.clicked.connect(self._delete)
        save = QtWidgets.QPushButton("Save as MSP…")
        save.clicked.connect(self._save)
        close = QtWidgets.QPushButton("Done")
        close.clicked.connect(self.accept)
        self.status = QtWidgets.QLabel("")
        top = QtWidgets.QHBoxLayout()
        top.addWidget(self.filter, 1)
        top.addWidget(self.show_decoys)
        bottom = QtWidgets.QHBoxLayout()
        bottom.addWidget(delete)
        bottom.addWidget(save)
        bottom.addStretch(1)
        bottom.addWidget(close)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(self.table, 1)
        lay.addWidget(self.plot)
        lay.addWidget(self.status)
        lay.addLayout(bottom)
        self._rows: list[int] = []          # table row -> index into self.entries
        self._fill()

    def _fill(self, *_a) -> None:
        text = self.filter.text().strip().lower()
        self._updating = True
        self._rows = [i for i, e in enumerate(self.entries)
                      if (self.show_decoys.isChecked() or not e.is_decoy)
                      and (not text or text in f"{e.name} {e.formula} {e.cas}".lower())]
        self.table.setRowCount(len(self._rows))
        for r, i in enumerate(self._rows):
            e = self.entries[i]
            values = [e.name, "" if e.ri is None else f"{e.ri:g}", e.formula, e.cas, str(e.mz.size),
                      "yes" if e.is_decoy else ""]
            for c, v in enumerate(values):
                item = QtWidgets.QTableWidgetItem(v)
                if c >= 4:
                    item.setFlags(item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(r, c, item)
        self._updating = False
        self.status.setText(f"{len(self._rows)} of {len(self.entries)} entries shown")

    def _edited(self, item) -> None:
        if self._updating:
            return
        e = self.entries[self._rows[item.row()]]
        text = item.text().strip()
        if item.column() == 0 and text:
            e.name = text
        elif item.column() == 1:
            try:
                e.ri = float(text) if text else None
            except ValueError:
                item.setText("" if e.ri is None else f"{e.ri:g}")
        elif item.column() == 2:
            e.formula = text
        elif item.column() == 3:
            e.cas = text
        self._rebuild()

    def _rebuild(self) -> None:
        self.library = build_library(self.entries)

    def _draw(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if len(rows) != 1:
            self.sticks.setData([], [])
            return
        e = self.entries[self._rows[rows[0]]]
        if e.mz.size == 0:
            self.sticks.setData([], [])
            return
        y = 100.0 * e.intensity / e.intensity.max()
        self.sticks.setData(np.repeat(e.mz.astype(float), 2), np.column_stack([np.zeros_like(y), y]).ravel())

    def _delete(self) -> None:
        drop = {self._rows[r] for r in {i.row() for i in self.table.selectedIndexes()}}
        if not drop:
            return
        self.entries = [e for i, e in enumerate(self.entries) if i not in drop]
        self._rebuild()
        self._fill()

    def _save(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save library", "library.msp", "MSP library (*.msp)")
        if path:
            self.save_to(path)

    def save_to(self, path) -> None:
        real = [e for e in self.entries if not e.is_decoy]
        write_msp(real, path)
        self.status.setText(f"Wrote {len(real)} entries to {path}")


class IdentifyDialog(QtWidgets.QDialog):
    def __init__(self, project, parent=None):
        super().__init__(parent)
        self.setWindowTitle("🏷️ Identify")
        self.resize(920, 640)
        self.project = project
        self._name = None
        self._run = None
        self._fold = None
        self._blobs = None
        self._components = None
        self._series = None
        self.library = None
        self.results: list = []

        top = QtWidgets.QHBoxLayout()
        load_btn = QtWidgets.QPushButton("Load library (MSP)…")
        load_btn.clicked.connect(self._load_library)
        self.library_label = QtWidgets.QLabel("No library loaded.")
        self.library_label.setWordWrap(True)
        self.auto_decoys = QtWidgets.QCheckBox("Generate decoys automatically")
        self.auto_decoys.setChecked(True)
        self.auto_decoys.setToolTip("Add a shuffled decoy per real entry, for the target-decoy q-value -- "
                                    "turn off if the library already carries its own decoys.")
        self.identify_btn = QtWidgets.QPushButton("🏷️ Identify peaks")
        self.identify_btn.setEnabled(False)
        self.identify_btn.clicked.connect(self._run_identify)
        self.settings = IdentifySettings()
        settings_btn = QtWidgets.QPushButton("⚙️")
        settings_btn.setToolTip("Identification settings: RI window and penalty, score weights, thresholds")
        settings_btn.clicked.connect(self._edit_settings)
        browse_btn = QtWidgets.QPushButton("Browse library…")
        browse_btn.setToolTip("List, search, edit, delete and save the loaded library's entries")
        browse_btn.clicked.connect(self._browse_library)
        top.addWidget(load_btn)
        top.addWidget(browse_btn)
        top.addWidget(self.library_label, 1)
        top.addWidget(self.auto_decoys)
        top.addWidget(settings_btn)
        top.addWidget(self.identify_btn)

        self.table = QtWidgets.QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(COLUMNS)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)

        buttons = QtWidgets.QHBoxLayout()
        promote_btn = QtWidgets.QPushButton("Promote to target…")
        promote_btn.clicked.connect(self._promote_selected)
        inhouse_btn = QtWidgets.QPushButton("Add to in-house library…")
        inhouse_btn.clicked.connect(self._add_to_inhouse)
        buttons.addWidget(promote_btn)
        buttons.addWidget(inhouse_btn)
        buttons.addStretch(1)

        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)

        lay = QtWidgets.QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(self.table, 1)
        lay.addLayout(buttons)
        lay.addWidget(self.status)

    # ------------------------------------------------------------------ context
    def set_run(self, name, run, fold, blobs, series, components=None) -> None:
        """Called by the main window with whichever run is open, and its
        own alkane series if ⛓️ has found one (for the retention index;
        `None` is fine, the search just runs on the spectrum alone)."""
        self._name, self._run, self._fold = name, run, fold
        self._blobs, self._series, self._components = blobs, series, components
        self.setWindowTitle(f"🏷️ Identify -- {name}" if name else "🏷️ Identify")
        self.table.setRowCount(0)
        self.results = []
        self.status.setText("")
        self._refresh_identify_enabled()

    def _refresh_identify_enabled(self) -> None:
        energy = self._run.meta.get('energy_ev') if self._run is not None else None
        if energy is not None and energy != 70:
            self.identify_btn.setEnabled(False)
            self.status.setText(f'{energy:g} eV data: the current EI identification workflow requires 70 eV. Select the 70 eV channel.')
            return
        has_components = self._components is not None and len(self._components.components) > 0
        has_blobs = self._blobs is not None and len(self._blobs) > 0
        self.identify_btn.setEnabled(self.library is not None and (has_components or has_blobs))

    # ------------------------------------------------------------------ library
    def _load_library(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Load MSP library", "", "MSP library (*.msp)")
        if path:
            self.load_library_file(path)

    def load_library_file(self, path) -> None:
        """The actual loading logic, separate from the file-picker above
        so it can be driven directly (scripted, or from a test) with a
        path already in hand."""
        try:
            entries = load_msp(path)
        except OSError as exc:
            self.status.setText(f"Could not read {path}: {exc}")
            return
        if not entries:
            self.status.setText(f"{path} parsed to zero entries -- check the file is a plain-text MSP export.")
            return
        if self.auto_decoys.isChecked():
            entries = entries + generate_decoys(entries)
        self.library = build_library(entries)
        n_real = sum(1 for e in self.library.entries if not e.is_decoy)
        n_decoy = len(self.library.entries) - n_real
        self.library_label.setText(f"{Path(path).name}: {n_real} compound(s)"
                                   + (f" + {n_decoy} decoy(s)" if n_decoy else ""))
        self._refresh_identify_enabled()

    def _edit_settings(self) -> None:
        from .dialogs import SettingsDialog
        dlg = SettingsDialog("🏷️ Identification settings", self.settings, parent=self)
        if dlg.exec():
            self.settings = dlg.settings()
            self.status.setText("Settings changed. Identify again to use them.")

    def _browse_library(self) -> None:
        if self.library is None:
            self.status.setText("Load a library first.")
            return
        dlg = LibraryBrowser(self.library, self)
        dlg.exec()
        self.library = dlg.library
        n_real = sum(1 for e in self.library.entries if not e.is_decoy)
        self.library_label.setText(f"{n_real} compound(s) (edited)")
        self._refresh_identify_enabled()

    # ------------------------------------------------------------------ search
    def _run_identify(self) -> None:
        if self.library is None or (self._components is None and self._blobs is None):
            return
        library, blobs, components = self.library, self._blobs, self._components
        series, settings = self._series, self.settings
        self.identify_btn.setEnabled(False)
        self.status.setText("Identifying…")

        def work(token, progress):
            if components is not None and components.components:
                return identify_components(components, library, series, settings)
            return identify_blobs(blobs, library, series, settings)

        def done(results):
            self.results = results
            self._fill_table()
            self._record_identifications()
            self._refresh_identify_enabled()
            n_hit = sum(1 for r in results if r.best is not None)
            source = "component" if components is not None and components.components else "legacy blob"
            self.status.setText(f"{n_hit} of {len(results)} {source}(s) matched.")

        def failed(text):
            self._refresh_identify_enabled()
            self.status.setText(f"Failed: {text}")

        self.project.jobs.submit("🏷️ Identifying", work, done, failed)

    def _fill_table(self) -> None:
        features = self._feature_list()
        self.table.setRowCount(len(self.results))
        for r, ident in enumerate(self.results):
            best = ident.best
            feature = features[ident.blob_index]
            values = [str(ident.blob_index), f"{feature['t1'] / 60:.3f}",
                     f"{feature['t2']:.3f}"]
            if best is None:
                values += ["", "", "", "", "", ""]
            else:
                values += [best.name, f"{best.score:.3f}",
                          "" if best.ri is None else f"{best.ri:.0f}",
                          "" if best.ri_diff is None else f"{best.ri_diff:.0f}",
                          best.tier or "", "" if best.q_value is None else f"{best.q_value:.3f}"]
            for c, v in enumerate(values):
                self.table.setItem(r, c, QtWidgets.QTableWidgetItem(v))

    def _record_identifications(self) -> None:
        """Leave this run's confident hits (tier 2a or 2b, never a decoy) on the
        project, where ⚖️ Compare uses them to pair peaks by identity."""
        if not self._name:
            return
        features = self._feature_list()
        self.project.identifications[self._name] = [
            (float(features[r.blob_index]["t1"]), float(features[r.blob_index]["t2"]), r.best.name)
            for r in self.results
            if r.best is not None and not r.best.is_decoy and r.best.tier in ("2a", "2b")]

    def _feature_list(self) -> list[dict]:
        if self._components is not None and self._components.components:
            return [{"t1": component.t1, "t2": component.t2,
                     "mz": component.mz, "spectrum": component.spectrum}
                    for component in self._components.components]
        blobs = self._blobs
        if blobs is None:
            return []
        use_apex = blobs.apex_spectra.shape[0] == len(blobs) and blobs.apex_spectra.size
        spectra = blobs.apex_spectra if use_apex else blobs.spectra
        return [{"t1": float(blobs.t1[index]), "t2": float(blobs.t2[index]),
                 "mz": blobs.mz, "spectrum": spectra[index]}
                for index in range(len(blobs))]

    def _selected_identification(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if len(rows) != 1 or not self.results:
            return None
        return self.results[rows[0]]

    # ------------------------------------------------------------------ promote / in-house
    def _promote_selected(self) -> None:
        ident = self._selected_identification()
        if ident is None or ident.best is None:
            self.status.setText("Select one identified row first.")
            return
        default = str(self.project.targets_path) if self.project.targets_path else "target_table.csv"
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Append to target table", default, "CSV (*.csv)")
        if path:
            self.promote_to_file(ident, path)

    def promote_to_file(self, ident, path) -> None:
        """The actual write, separate from the save-file dialog above."""
        entry = self.library.entries[ident.best.entry_index]
        feature = self._feature_list()[ident.blob_index]
        target = promote_to_target(entry, ident.best, feature["t1"], feature["t2"])
        append = Path(path).exists()
        write_targets_csv([target], path, append=append)
        self.status.setText(f"{'Appended' if append else 'Wrote'} {target.compound} to {path}. "
                            "Reload the target table to use it.")

    def _add_to_inhouse(self) -> None:
        ident = self._selected_identification()
        if ident is None:
            self.status.setText("Select one row first.")
            return
        default_name = ident.best.name if ident.best is not None else ""
        name, ok = QtWidgets.QInputDialog.getText(self, "Compound name", "Name for the in-house library:",
                                                  text=default_name)
        if not ok or not name.strip():
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Add to in-house library", "inhouse.msp",
                                                        "MSP library (*.msp)")
        if path:
            self.add_to_inhouse_file(ident, name.strip(), path)

    def add_to_inhouse_file(self, ident, name: str, path) -> None:
        """The actual write, separate from the name prompt and save-file
        dialog above."""
        feature = self._feature_list()[ident.blob_index]
        spectrum = np.asarray(feature["spectrum"], dtype=float)
        mz = np.asarray(feature["mz"], dtype=int)
        ri = None
        if self._series is not None and len(self._series) >= 2:
            try:
                ri = retention_index(self._series, float(feature["t1"]))
            except ValueError:
                ri = None
        keep = spectrum > 0
        entry = LibraryEntry(name=name, mz=mz[keep], intensity=spectrum[keep], ri=ri)
        existing = load_msp(path) if Path(path).exists() else []
        write_msp(existing + [entry], path)
        self.status.setText(f"Added {entry.name} to {path}.")
