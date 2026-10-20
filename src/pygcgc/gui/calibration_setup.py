"""A transactional setup wizard connecting named peaks, runs, amounts and units."""

from __future__ import annotations

import math
import re
from dataclasses import replace

from PySide6 import QtCore, QtWidgets
from .theme import fit_to_screen

from ..core.calibration_inputs import CalibrationInputs, MASS_UNITS, filename_amount, calibration_coverage_errors
from ..core.istd import references_to_targets
from ..core.targets import Target, load_targets


class _Page(QtWidgets.QWizardPage):
    def __init__(self, title, initialise=None, validate=None):
        super().__init__()
        self.setTitle(title)
        self.initialise, self.validate = initialise, validate

    def initializePage(self):
        if self.initialise:
            self.initialise()

    def validatePage(self):
        return self.validate() if self.validate else True


def _table(headers):
    widget = QtWidgets.QTableWidget(0, len(headers))
    widget.setHorizontalHeaderLabels(headers)
    widget.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
    widget.horizontalHeader().setStretchLastSection(True)
    widget.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
    return widget


def _text(table, row, column):
    item = table.item(row, column)
    return item.text().strip() if item is not None else ""


class CalibrationSetupWizard(QtWidgets.QWizard):
    applied = QtCore.Signal()

    def __init__(self, project, parent=None, selected_target=None):
        super().__init__(parent)
        self.project = project
        self.setWindowTitle("Set up calibration and units")
        self.setWizardStyle(QtWidgets.QWizard.WizardStyle.ModernStyle)
        fit_to_screen(self, 1000, 680)
        self.setOption(QtWidgets.QWizard.WizardOption.NoBackButtonOnStartPage)
        self.setButtonText(QtWidgets.QWizard.WizardButton.FinishButton, "Apply setup")
        self._original = {t.compound: t for t in project.quant.targets or []}
        self._levels = {lv: dict(amounts) for lv, amounts in (project.quant.levels or {}).items()}
        self._amount_widgets = {}
        self._last_unit = project.quant.amount_unit
        if self._last_unit == "calibration units":
            self._last_unit = "ng"
        self._initial_standards = dict(project.quant.calibration_setup.get("internal_standard_amounts", {}))
        for t in self._original.values():
            values = {a[t.compound] for a in self._levels.values() if t.compound in a}
            if t.is_istd and len(values) == 1:
                self._initial_standards.setdefault(t.compound, values.pop())
        self._draft = None

        compounds = _Page("1. Name compounds and assign internal standards", validate=self._validate_targets)
        compounds.setSubTitle("Import your target table or add named peaks. Positions are seconds; "
                              "leave a target's internal standard blank to use the nearest one.")
        self.targets_table = _table(["Compound", "Role", "Quantifier ion", "t1 (s)", "t2 (s)",
                                     "Internal standard", "Qualifier ions (;)", "Qualifier ratios (;)"])
        import_button = QtWidgets.QPushButton("Import target table…")
        import_button.clicked.connect(self.import_targets)
        add_button = QtWidgets.QPushButton("Add compound")
        add_button.clicked.connect(lambda: self.add_target(Target("")))
        reference_button = QtWidgets.QPushButton("Add labelled internal standards")
        reference_button.clicked.connect(self.add_references)
        peak_button = QtWidgets.QPushButton("Add selected peak")
        self.selected_peak_btn = peak_button
        peak_button.setToolTip("Select exactly one peak in Peaks before opening setup. Its position and "
                              "suggested quantifier are copied; supply and verify its compound name.")
        peak_button.setEnabled(selected_target is not None)
        peak_button.clicked.connect(lambda: self.add_target(replace(selected_target)) if selected_target else None)
        delete_button = QtWidgets.QPushButton("Remove selected")
        delete_button.clicked.connect(self._remove_target)
        row = QtWidgets.QHBoxLayout()
        for button in (import_button, add_button, peak_button, reference_button, delete_button):
            row.addWidget(button)
        self.target_error = QtWidgets.QLabel()
        self.target_error.setWordWrap(True)
        layout = QtWidgets.QVBoxLayout(compounds)
        layout.addLayout(row)
        layout.addWidget(self.targets_table, 1)
        layout.addWidget(self.target_error)
        self.addPage(compounds)
        for target in self._original.values():
            self.add_target(target)
        self.add_references()

        runs = _Page("2. Choose calibration runs and their levels", validate=self._validate_runs)
        runs.setSubTitle("Run types, calibration sets and sample air volumes are saved with the project. "
                         "Filename level names are suggestions; inspect them before continuing.")
        self.runs_table = _table(["Run", "Type", "Level", "Calibration set", "Air volume (m³)"])
        for i in project.loaded():
            e = project.entries[i]
            row = self.runs_table.rowCount()
            self.runs_table.insertRow(row)
            match = re.search(r"(?:^|_)L(\d+)(?:_|$)", e.name, re.I)
            level = e.level or (f"L{match[1]}" if match else "")
            for column, value in ((0, e.name), (2, level), (3, e.cal_set),
                                  (4, "" if e.volume is None else str(e.volume))):
                item = QtWidgets.QTableWidgetItem(value)
                if column == 0:
                    item.setFlags(item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
                self.runs_table.setItem(row, column, item)
            kind = QtWidgets.QComboBox()
            kind.setEditable(True)
            kind.addItems(["sample", "calibration", "tube blank"])
            kind.setCurrentText(e.kind or ("calibration" if match else "sample"))
            self.runs_table.setCellWidget(row, 1, kind)
        self.run_error = QtWidgets.QLabel()
        self.run_error.setWordWrap(True)
        layout = QtWidgets.QVBoxLayout(runs)
        layout.addWidget(self.runs_table, 1)
        layout.addWidget(self.run_error)
        self.addPage(runs)

        amounts = _Page("3. Enter known amounts and their unit", initialise=self._build_amounts,
                        validate=self._validate_amounts)
        amounts.setSubTitle("Use actual amounts on tube, not stock concentrations. Internal-standard "
                            "amounts must be constant. Curve settings require 3 usable levels by default.")
        self.unit_box = QtWidgets.QComboBox()
        self.unit_box.setEditable(True)
        self.unit_box.addItems(["ng", "pg", "ug", "mg"])
        self.unit_box.setCurrentText(self._last_unit)
        self.unit_note = QtWidgets.QLabel("One amount unit is used for calibration inputs and outputs.")
        self.unit_note.setWordWrap(True)
        self.unit_box.currentTextChanged.connect(lambda _text: self.unit_note.setText(
            "Finish editing to convert compatible mass units; custom units retain numbers for checking."))
        self.unit_box.activated.connect(lambda _index: self._commit_unit())
        self.unit_box.lineEdit().editingFinished.connect(self._commit_unit)
        fill = QtWidgets.QPushButton("Copy unambiguous filename suggestions into blank target amounts")
        fill.clicked.connect(self._fill_suggestions)
        self.amounts_table = _table(["Level", "Filename suggestion"])
        self.standards_form = QtWidgets.QFormLayout()
        self.amount_error = QtWidgets.QLabel()
        self.amount_error.setWordWrap(True)
        layout = QtWidgets.QVBoxLayout(amounts)
        unit_row = QtWidgets.QHBoxLayout()
        unit_row.addWidget(QtWidgets.QLabel("Amounts entered and reported in:"))
        unit_row.addWidget(self.unit_box)
        unit_row.addWidget(self.unit_note, 1)
        layout.addLayout(unit_row)
        layout.addWidget(fill)
        layout.addWidget(self.amounts_table, 1)
        layout.addLayout(self.standards_form)
        layout.addWidget(self.amount_error)
        self.addPage(amounts)

        review = _Page("4. Review and apply", initialise=self._build_review,
                       validate=self._validate_review)
        self.summary = QtWidgets.QTextBrowser()
        self.confirm = QtWidgets.QCheckBox(
            "I have checked the compounds, amounts and units; each internal standard was added "
            "at the stated constant amount to every calibration and sample run.")
        self.review_error = QtWidgets.QLabel()
        self.review_error.setWordWrap(True)
        layout = QtWidgets.QVBoxLayout(review)
        layout.addWidget(self.summary, 1)
        layout.addWidget(self.confirm)
        layout.addWidget(self.review_error)
        self.addPage(review)

    def add_target(self, target):
        row = self.targets_table.rowCount()
        self.targets_table.insertRow(row)
        values = [target.compound, target.role, str(target.quant_ion), f"{target.t1:.15g}",
                  f"{target.t2:.15g}", target.istd, ";".join(map(str, target.qual_ions)),
                  ";".join(map(str, target.qual_ratios))]
        for column, value in enumerate(values):
            self.targets_table.setItem(row, column, QtWidgets.QTableWidgetItem(value))

    def add_references(self):
        existing = {_text(self.targets_table, r, 0) for r in range(self.targets_table.rowCount())}
        for target in references_to_targets(list(self.project.istd.references.values())):
            if target.compound not in existing:
                self.add_target(target)

    def import_targets(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Import target table", "", "CSV (*.csv)")
        if not path:
            return
        try:
            imported = load_targets(path)
        except (OSError, ValueError, KeyError) as exc:
            self.target_error.setText(str(exc))
            return
        existing = {_text(self.targets_table, r, 0) for r in range(self.targets_table.rowCount())}
        for target in imported:
            self._original[target.compound] = target
            if target.compound not in existing:
                self.add_target(target)
                existing.add(target.compound)

    def _remove_target(self):
        for row in sorted({i.row() for i in self.targets_table.selectedIndexes()}, reverse=True):
            self.targets_table.removeRow(row)

    def targets(self):
        result = []
        for row in range(self.targets_table.rowCount()):
            v = [_text(self.targets_table, row, c) for c in range(8)]
            parse = lambda text, cast: tuple(cast(x.strip()) for x in text.split(";") if x.strip())
            result.append(replace(self._original.get(v[0], Target(v[0])), compound=v[0], role=v[1].lower(),
                                  quant_ion=int(v[2]), t1=float(v[3]), t2=float(v[4]), istd=v[5],
                                  qual_ions=parse(v[6], int), qual_ratios=parse(v[7], float)))
        return result

    def runs(self):
        return {_text(self.runs_table, r, 0): dict(kind=self.runs_table.cellWidget(r, 1).currentText().strip(),
                  level=_text(self.runs_table, r, 2), cal_set=_text(self.runs_table, r, 3),
                  volume=float(_text(self.runs_table, r, 4)) if _text(self.runs_table, r, 4) else None)
                for r in range(self.runs_table.rowCount())}

    def _validate_targets(self):
        try:
            targets = self.targets()
            probe = CalibrationInputs(targets=targets, confirmed=True)
            errors = [e for e in probe.errors() if not any(s in e for s in
                      ("constant across", "known levels", "distinct non-negative"))]
            if errors:
                raise ValueError("\n".join(errors))
        except (ValueError, TypeError) as exc:
            self.target_error.setText(str(exc))
            return False
        self.target_error.clear()
        return True

    def _validate_runs(self):
        try:
            runs = self.runs()
            calibrations = [r for r in runs.values() if "cal" in r["kind"].lower()]
            if not calibrations or any(not r["level"] for r in calibrations):
                raise ValueError("Mark calibration runs and give each a level name.")
            if any(r['volume'] is not None and (not math.isfinite(r['volume']) or r['volume'] <= 0)
                   for r in runs.values()):
                raise ValueError("Air volume must be positive or blank.")
        except ValueError as exc:
            self.run_error.setText(str(exc))
            return False
        self.run_error.clear()
        return True

    def _save_amounts(self):
        for row in range(self.amounts_table.rowCount()):
            level = _text(self.amounts_table, row, 0)
            values = self._levels.setdefault(level, {})
            for column in range(2, self.amounts_table.columnCount()):
                name = self.amounts_table.horizontalHeaderItem(column).data(QtCore.Qt.ItemDataRole.UserRole)
                text = _text(self.amounts_table, row, column)
                if text:
                    values[name] = float(text)
                else:
                    values.pop(name, None)

    def _suggestion(self, level):
        sources = [name for name, info in self.runs().items() if info['level'] == level
                   and 'cal' in info['kind'].lower()]
        proposed = [filename_amount(n, self.unit_box.currentText()) for n in sources]
        valid = {v for v, _ in proposed if v is not None}
        if proposed and len(valid) == 1 and all(v is not None for v, _ in proposed):
            return proposed[0]
        return None, "; ".join(dict.fromkeys(note for _, note in proposed)) or "Enter known amounts"

    def _build_amounts(self):
        self._save_amounts()
        standard_values = {name: edit.text() for name, edit in self._amount_widgets.items()}
        self._amount_widgets = {}
        while self.standards_form.rowCount():
            self.standards_form.removeRow(0)
        targets = self.targets()
        analytes = [t for t in targets if not t.is_istd]
        self.amounts_table.setColumnCount(2 + len(analytes))
        self.amounts_table.setHorizontalHeaderLabels(["Level", "Filename suggestion"] + [t.compound for t in analytes])
        for col, target in enumerate(analytes, start=2):
            self.amounts_table.horizontalHeaderItem(col).setData(QtCore.Qt.ItemDataRole.UserRole, target.compound)
        levels = list(dict.fromkeys(r['level'] for r in self.runs().values() if 'cal' in r['kind'].lower()))
        self.amounts_table.setRowCount(len(levels))
        for row, level in enumerate(levels):
            for col, value in ((0, level), (1, self._suggestion(level)[1])):
                item = QtWidgets.QTableWidgetItem(value)
                item.setFlags(item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
                self.amounts_table.setItem(row, col, item)
            for col, target in enumerate(analytes, start=2):
                value = self._levels.get(level, {}).get(target.compound)
                self.amounts_table.setItem(row, col, QtWidgets.QTableWidgetItem("" if value is None else f"{value:.12g}"))
        for target in targets:
            if target.is_istd:
                edit = QtWidgets.QLineEdit(standard_values.get(target.compound,
                    str(self._initial_standards.get(target.compound, ""))))
                edit.setPlaceholderText("Enter actual amount; no assumed 1 ng")
                self._amount_widgets[target.compound] = edit
                self.standards_form.addRow(f"{target.compound}: amount in every run", edit)

    def _fill_suggestions(self):
        for row in range(self.amounts_table.rowCount()):
            value, _ = self._suggestion(_text(self.amounts_table, row, 0))
            if value is not None:
                for col in range(2, self.amounts_table.columnCount()):
                    if not _text(self.amounts_table, row, col):
                        self.amounts_table.setItem(row, col, QtWidgets.QTableWidgetItem(f"{value:.12g}"))

    def _unit_changed(self, unit):
        previous = self._last_unit
        self._last_unit = unit
        if unit == previous:
            return
        try:
            self._save_amounts()
            if previous in MASS_UNITS and unit in MASS_UNITS:
                entered = {name: float(edit.text()) for name, edit in self._amount_widgets.items() if edit.text().strip()}
                factor = MASS_UNITS[previous]/MASS_UNITS[unit]
                self._levels = {lv: {name: value*factor for name, value in values.items()}
                                for lv, values in self._levels.items()}
                for name, value in entered.items():
                    self._amount_widgets[name].setText(f"{value*factor:.12g}")
                # Avoid saving the old table over the converted values on rebuild.
                self.amounts_table.setRowCount(0)
                self._build_amounts()
                self.unit_note.setText(f"Converted entered amounts from {previous} to {unit}.")
            else:
                self.unit_note.setText("Custom unit: values are unchanged. Enter the correct numeric amounts "
                                       "in this unit and confirm them on the review page.")
        except ValueError:
            self._last_unit = previous
            self.unit_box.blockSignals(True)
            self.unit_box.setCurrentText(previous)
            self.unit_box.blockSignals(False)
            self.amount_error.setText("Correct non-numeric amounts before changing units.")

    def _commit_unit(self):
        # Editable combo boxes emit incomplete strings for each keystroke. Only
        # convert a completed unit, never an intermediate prefix such as "u".
        unit = self.unit_box.currentText().strip()
        if unit and unit != self._last_unit:
            self._unit_changed(unit)

    def draft(self, confirmed=False):
        self._commit_unit()
        self._save_amounts()
        standards = {name: float(edit.text()) for name, edit in self._amount_widgets.items() if edit.text().strip()}
        levels = {lv: {t.compound: amounts[t.compound] for t in self.targets()
                       if t.compound in amounts and not t.is_istd}
                  for lv, amounts in self._levels.items()}
        for amounts in levels.values():
            amounts.update(standards)
        return CalibrationInputs(self.targets(), levels, self.runs(), self.unit_box.currentText().strip(),
                                 standards, confirmed)

    def _validate_amounts(self):
        try:
            draft = self.draft(confirmed=True)
            errors = draft.errors() + self._coverage_errors(draft)
            if errors:
                raise ValueError("\n".join(errors))
        except ValueError as exc:
            self.amount_error.setText(str(exc))
            return False
        self.amount_error.clear()
        return True

    def _coverage_errors(self, draft):
        q = self.project.quant
        return calibration_coverage_errors(draft.targets, draft.levels, draft.runs, q.settings.calibration, q.decisions)

    def _build_review(self):
        draft = self.draft()
        self.confirm.setChecked(False)
        lines = [f"Amount unit: {draft.amount_unit}", "Output concentration: amount unit / m³ (when air volume is given)",
                 "", "Targets and internal-standard assignment:"]
        lines += [f"  {t.compound}: {t.role}, m/z {t.quant_ion}, t1={t.t1:g} s, t2={t.t2:g} s; "
                  f"IS={t.istd or 'nearest'}" for t in draft.targets]
        lines += ["", "Known amounts by level:"]
        lines += [f"  {lv}: " + ", ".join(f"{name}={amount:g}" for name, amount in values.items())
                  for lv, values in draft.levels.items()]
        lines += ["", "Runs:"] + [f"  {name}: {info['kind']}, {info['level']}, set={info['cal_set'] or '(default)'}"
                                    for name, info in draft.runs.items()]
        lines += ["", "Applying saves inputs in the project and marks existing quantitative results out of date.",
                  "No raw data are changed and no fitting or quantification job starts automatically."]
        self.summary.setPlainText("\n".join(lines))

    def _validate_review(self):
        try:
            if self.project.quant.running:
                raise ValueError("Wait for quantification to finish before applying setup.")
            self._draft = self.draft(self.confirm.isChecked())
            errors = self._draft.errors() + self._coverage_errors(self._draft)
            if errors:
                raise ValueError("\n".join(errors))
            if set(self._draft.runs) != {self.project.entries[i].name for i in self.project.loaded()}:
                raise ValueError("The loaded runs changed while setup was open. Reopen setup for this sequence.")
        except ValueError as exc:
            self.review_error.setText(str(exc))
            return False
        return True

    def accept(self):
        if not self._validate_review():
            return
        draft = self._draft
        q = self.project.quant
        q.targets, q.targets_error = draft.targets, ""
        q._inline_targets = True
        self.project.targets_path = None
        from ..core.phase import Anchor
        self.project.anchors = [Anchor(t.compound, t.quant_ion, t.t1, t.t2, t.t1_window)
                               for t in draft.targets if t.is_istd]
        q.levels, q.levels_path = draft.levels, None
        q.amount_unit = draft.amount_unit
        q.calibration_setup = {"confirmed": True, "internal_standard_amounts": draft.internal_standard_amounts,
                               "constant_internal_standards": True}
        for name, fields in draft.runs.items():
            entry = self.project.entries[self.project.index_of(name)]
            for key, value in fields.items():
                setattr(entry, key, value)
        q._invalidate()
        self.project.anchorsChanged.emit()
        self.project.runsChanged.emit()
        self.project.request_autosave()
        self.applied.emit()
        super().accept()
