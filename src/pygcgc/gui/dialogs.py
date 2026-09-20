"""⚙️ Quantification settings dialog."""

from __future__ import annotations

from dataclasses import fields, replace

from PySide6 import QtWidgets

from ..core import QuantSettings

# field: (label, tooltip, (min, max, step, decimals) or choices)
PEAK_FIELDS = {
    "target_snr_min": ("Target peaks: minimum S/N", "Quantifier-ion peaks weaker than this are ignored", (1, 1000, 1, 1)),
    "snr_min": ("TIC peaks: minimum S/N", "Seeds for 'all signals' peaks", (1, 10000, 1, 1)),
    "floor_snr": ("Footprint edge (S/N)", "Peaks are integrated down to this level", (0.5, 50, 0.5, 1)),
    "wrap_below_s": ("Wrap boundary (s)", "Peaks with t₂ below this belong to the previous modulation. "
                     "Put it in the quiet band just above zero.", (0, 3, 0.05, 2)),
    "min_modulations": ("Minimum modulations", "Peaks must span at least this many modulations", (1, 20, 1, 0)),
    "rim_baseline": ("Subtract local level", "Remove the level around a target (a UCM hump, say)", None),
    "grow_rows": ("Widen footprint (rows)", "Recovers tails of small peaks", (0, 20, 1, 0)),
}
QUANT_FIELDS = {
    "qual_tolerance": ("Qualifier tolerance (%)", "Allowed relative deviation of each qualifier ratio", (1, 200, 1, 0)),
    "t1_tolerance": ("Position tolerance t₁ (s)", "How far from its expected t₁ a target may be", (0.5, 60, 0.5, 1)),
    "t2_tolerance": ("Position tolerance t₂ (s)", "How far from its expected t₂ a target may be", (0.005, 1, 0.005, 3)),
    "ambiguity": ("Ambiguity threshold", "Flag when the runner-up scores at least this fraction of the best",
                  (0.1, 1.0, 0.05, 2)),
    "use_istd_offsets": ("Follow internal standard shifts", "Move each target's expected position with its IS", None),
    "learn_positions": ("Learn positions from calibration runs", "Replace table positions with measured ones", None),
    "t1_scale": ("Distance scale t₁ (s)", "For nearest IS and nearby peaks: this many seconds of t₁…", (1, 600, 5, 0)),
    "t2_scale": ("Distance scale t₂ (s)", "…count the same as this many seconds of t₂", (0.05, 4, 0.05, 2)),
    "surrogate_max_distance": ("Nearby-peak distance", "Largest scaled distance for a semi-quantitative amount",
                               (0.1, 50, 0.5, 1)),
    "is_tolerance": ("IS flag threshold (%)", "Flag internal standards this far from their median area", (5, 500, 5, 0)),
    "saturated_qualifier": ("Saturated target: quantify on a qualifier (optional)",
                            "When the quantifier ion is saturated, use the first unsaturated qualifier ion "
                            "instead: its volume divided by its expected ratio from the target table. Put an M+1 "
                            "ion in the qualifier list, with its ratio, to use one.", None),
    "interpolate_cal_sets": ("Interpolate between calibration sets", "A sample between two calibration blocks "
                             "(a 'cal_set' column in the sequence) is quantified against a blend of both, "
                             "weighted by how far along it sits. Off: use the nearest set before it.", None),
}
CAL_FIELDS = {
    "weights": ("Weights", "Auto picks the scheme with the smallest back-calculation error",
                ["auto", "1", "1/x", "1/x2"]),
    "model": ("Model", "Auto uses Mandel's test", ["auto", "linear", "quadratic"]),
    "force_zero": ("Force through zero", "", None),
    "include_zero": ("Use zero-amount levels in the fit", "", None),
    "accept_re": ("Back-calculation limit (%)", "", (1, 100, 1, 0)),
    "accept_re_lowest": ("…at the lowest level (%)", "", (1, 100, 1, 0)),
    "min_levels": ("Minimum levels", "", (2, 20, 1, 0)),
    "alpha": ("Significance level", "For Mandel's test and the intercept test", (0.001, 0.2, 0.005, 3)),
}


def _editor(value, spec):
    if isinstance(value, bool) or spec is None:
        w = QtWidgets.QCheckBox()
        w.setChecked(bool(value))
        return w, w.isChecked
    if isinstance(spec, list):
        w = QtWidgets.QComboBox()
        w.addItems(spec)
        w.setCurrentText(str(value))
        return w, w.currentText
    lo, hi, step, dec = spec
    if dec == 0:
        w = QtWidgets.QSpinBox()
        w.setRange(int(lo), int(hi))
        w.setSingleStep(int(step))
        w.setValue(int(value))
        return w, w.value
    w = QtWidgets.QDoubleSpinBox()
    w.setRange(lo, hi)
    w.setSingleStep(step)
    w.setDecimals(dec)
    w.setValue(float(value))
    return w, w.value


class SettingsDialog(QtWidgets.QDialog):
    """A form for any settings dataclass: one row per field, typed by the
    field's current value (bool, int, float, string, or a tuple of numbers
    typed as comma-separated text). A field that is itself a dataclass gets
    its own tab. `tips` maps field name to a tooltip; "Restore defaults"
    returns a fresh instance of the same class."""

    def __init__(self, title: str, settings, tips: dict | None = None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(460, 560)
        self._base = settings
        self._restore = False
        self._tips = tips or {}
        self._getters: list = []          # (path, callable)
        tabs = QtWidgets.QTabWidget()
        self._build(tabs, "General", settings, ())
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                             | QtWidgets.QDialogButtonBox.StandardButton.Cancel
                                             | QtWidgets.QDialogButtonBox.StandardButton.RestoreDefaults)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QtWidgets.QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(self._defaults)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(tabs)
        lay.addWidget(buttons)

    def _build(self, tabs, title, obj, path) -> None:
        page = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(page)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(page)
        nested = []
        for f in fields(obj):
            value = getattr(obj, f.name)
            if hasattr(value, "__dataclass_fields__"):
                nested.append((f.name, value))
                continue
            if isinstance(value, dict):
                continue                      # not representable in a form; left as it is
            w, get = _auto_editor(value)
            w.setToolTip(self._tips.get(f.name, ""))
            form.addRow(f.name.replace("_", " "), w)
            self._getters.append((path + (f.name,), get))
        tabs.insertTab(0 if not path else tabs.count(), scroll, title)
        for name, value in nested:
            self._build(tabs, name.replace("_", " "), value, path + (name,))

    def _defaults(self) -> None:
        self._restore = True
        self.accept()

    def settings(self):
        if self._restore:
            return type(self._base)()
        return _rebuild(self._base, (), dict(self._getters))


def _auto_editor(value):
    if isinstance(value, bool):
        w = QtWidgets.QCheckBox()
        w.setChecked(value)
        return w, w.isChecked
    if isinstance(value, int):
        w = QtWidgets.QSpinBox()
        w.setRange(-10**9, 10**9)
        w.setValue(value)
        return w, w.value
    if isinstance(value, float):
        w = QtWidgets.QDoubleSpinBox()
        w.setRange(-1e12, 1e12)
        w.setDecimals(4)
        w.setSingleStep(max(abs(value) / 10, 0.001))
        w.setValue(value)
        return w, w.value
    if isinstance(value, (tuple, list)):
        kind = type(value[0]) if len(value) else float
        w = QtWidgets.QLineEdit(", ".join(str(v) for v in value))
        return w, lambda: tuple(kind(x) for x in w.text().replace(";", ",").split(",") if x.strip())
    if value is None:
        # An optional number (for example "detector ceiling: from the data"): blank means none.
        w = QtWidgets.QLineEdit()
        w.setPlaceholderText("none (automatic)")
        return w, lambda: (float(w.text()) if w.text().strip() else None)
    w = QtWidgets.QLineEdit(str(value))
    return w, w.text


def _rebuild(obj, path, getters):
    changes = {}
    for f in fields(obj):
        value = getattr(obj, f.name)
        if hasattr(value, "__dataclass_fields__"):
            changes[f.name] = _rebuild(value, path + (f.name,), getters)
        elif (path + (f.name,)) in getters:
            try:
                changes[f.name] = getters[path + (f.name,)]()
            except ValueError:
                changes[f.name] = value       # unparseable text keeps the old value
    return replace(obj, **changes)


class QuantSettingsDialog(QtWidgets.QDialog):
    def __init__(self, settings: QuantSettings, parent=None):
        super().__init__(parent)
        self.setWindowTitle("🧪 Quantification settings")
        self._getters = {"peaks": {}, "quant": {}, "cal": {}}
        tabs = QtWidgets.QTabWidget()
        for key, title, obj, spec in (("quant", "Targets", settings, QUANT_FIELDS),
                                      ("peaks", "Peaks", settings.peaks, PEAK_FIELDS),
                                      ("cal", "Calibration", settings.calibration, CAL_FIELDS)):
            page = QtWidgets.QWidget()
            form = QtWidgets.QFormLayout(page)
            for name, (label, tip, rng) in spec.items():
                w, get = _editor(getattr(obj, name), rng)
                w.setToolTip(tip)
                form.addRow(label, w)
                self._getters[key][name] = get
            tabs.addTab(page, title)
        sat = QtWidgets.QWidget()
        sl = QtWidgets.QHBoxLayout(sat)
        sl.setContentsMargins(0, 0, 0, 0)
        self.sat_auto = QtWidgets.QCheckBox("from the data")
        self.sat_value = QtWidgets.QDoubleSpinBox()
        self.sat_value.setRange(1, 1e12)
        self.sat_value.setDecimals(0)
        self.sat_value.setValue(settings.saturation or 1e5)
        self.sat_auto.setChecked(settings.saturation is None)
        self.sat_value.setEnabled(settings.saturation is not None)
        self.sat_auto.toggled.connect(lambda on: self.sat_value.setEnabled(not on))
        sl.addWidget(self.sat_auto)
        sl.addWidget(self.sat_value)
        tabs.widget(0).layout().addRow("Detector ceiling", sat)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Ok
                                             | QtWidgets.QDialogButtonBox.StandardButton.Cancel
                                             | QtWidgets.QDialogButtonBox.StandardButton.RestoreDefaults)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        buttons.button(QtWidgets.QDialogButtonBox.StandardButton.RestoreDefaults).clicked.connect(self._defaults)
        self._restore = False
        note = QtWidgets.QLabel("Changing these marks the results as out of date.")
        note.setStyleSheet("color: palette(mid);")
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(tabs)
        lay.addWidget(note)
        lay.addWidget(buttons)
        self._base = settings

    def _defaults(self) -> None:
        self._restore = True
        self.accept()

    def settings(self) -> QuantSettings:
        if self._restore:
            return QuantSettings()
        g = self._getters
        peaks = replace(self._base.peaks, **{k: f() for k, f in g["peaks"].items()})
        cal = replace(self._base.calibration, **{k: f() for k, f in g["cal"].items()})
        sat = None if self.sat_auto.isChecked() else float(self.sat_value.value())
        return replace(self._base, peaks=peaks, calibration=cal, saturation=sat,
                       **{k: f() for k, f in g["quant"].items()})
