"""Validation and filename suggestions for the guided calibration setup.

Amounts use one declared unit throughout a setup. The table-driven quantifier
uses area ratios and assumes the amount of each internal standard is constant
across the sequence. No amount is inferred as confirmed scientific input.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

from .targets import Target
from .calibration import CalibrationSettings

MASS_UNITS = {"pg": .001, "ng": 1., "ug": 1000., "µg": 1000., "μg": 1000., "mg": 1e6}


def calibration_coverage_errors(targets, levels, runs, settings, decisions=None):
    """Check usable level coverage with the same minimum/zero/exclusion policy as fitting."""
    errors = []
    sets = {info.get('cal_set', '') for info in runs.values() if 'cal' in info.get('kind', '').lower()}
    for cal_set in sets:
        for target in targets:
            if target.is_istd:
                continue
            overrides = decisions.overrides.get(target.compound, {}) if decisions else {}
            policy = CalibrationSettings.from_dict({**settings.to_dict(), **overrides})
            usable = set()
            for name, info in runs.items():
                if 'cal' not in info.get('kind', '').lower() or info.get('cal_set', '') != cal_set:
                    continue
                if decisions and (name in decisions.excluded.get(target.compound, {})
                                  or decisions.status_of(name, target.compound) == 'rejected'):
                    continue
                value = levels.get(info.get('level', ''), {}).get(target.compound)
                if value is not None and math.isfinite(value) and (value > 0 or (policy.include_zero and value == 0)):
                    usable.add(value)
            required = max(2, policy.min_levels)
            if len(usable) < required:
                errors.append(f"{target.compound}, set {cal_set or '(default)'}: {len(usable)} usable known "
                              f"amounts; current curve settings require {required}. "
                              "Zero levels and excluded points count only when the curve settings allow them.")
    return errors


def filename_amount(name: str, unit: str = "ng") -> tuple[float | None, str]:
    """Suggest an explicit amount; an ambiguous leading-zero token stays unknown."""
    found = re.findall(r"(?<![\d.eE+\-])(\d+(?:\.\d+)?)\s*(pg|ng|ug|µg|μg|mg)(?![A-Za-z])", name, re.I)
    if len(found) != 1:
        return None, "No unique filename amount"
    text, source_unit = found[0]
    if len(text) > 1 and text.startswith("0") and "." not in text:
        return None, f"Ambiguous filename amount: {text}{source_unit}; enter the actual amount"
    if unit not in MASS_UNITS:
        return None, f"Filename says {text}{source_unit}; enter in {unit}"
    value = float(text)*MASS_UNITS[source_unit.lower()]/MASS_UNITS[unit]
    if not math.isfinite(value):
        return None, "Filename amount is not finite; enter the actual amount"
    return value, f"Filename suggestion: {value:g} {unit}; confirm before use"


@dataclass
class CalibrationInputs:
    targets: list[Target] = field(default_factory=list)
    levels: dict[str, dict[str, float]] = field(default_factory=dict)
    runs: dict[str, dict] = field(default_factory=dict)
    amount_unit: str = "ng"
    internal_standard_amounts: dict[str, float] = field(default_factory=dict)
    confirmed: bool = False

    def errors(self) -> list[str]:
        errors = []
        names = [t.compound.strip() for t in self.targets]
        if not names or any(not n for n in names):
            errors.append("Give every compound a name and include at least one target.")
        if len({n.casefold() for n in names}) != len(names):
            errors.append("Compound names must be unique.")
        standards = {t.compound for t in self.targets if t.is_istd}
        analytes = [t for t in self.targets if not t.is_istd]
        if not standards or not analytes:
            errors.append("Include at least one target and one internal standard.")
        for t in self.targets:
            if t.role not in ("target", "istd") or t.quant_ion <= 0:
                errors.append(f"{t.compound}: choose a valid role and positive quantifier ion.")
            if not all(math.isfinite(v) and v >= 0 for v in (t.t1, t.t2)):
                errors.append(f"{t.compound}: retention times must be finite and non-negative.")
            if len(t.qual_ions) != len(t.qual_ratios) and t.qual_ratios:
                errors.append(f"{t.compound}: qualifier ions and ratios must have matching lengths.")
            if any(ion <= 0 for ion in t.qual_ions) or any(not math.isfinite(v) or v < 0 for v in t.qual_ratios):
                errors.append(f"{t.compound}: qualifier ions must be positive and ratios finite and non-negative.")
            if t.istd and t.istd not in standards:
                errors.append(f"{t.compound}: choose an internal standard from this table.")
        if not self.amount_unit.strip():
            errors.append("Declare the amount unit.")
        for name in standards:
            amount = self.internal_standard_amounts.get(name)
            if amount is None or not math.isfinite(amount) or amount <= 0:
                errors.append(f"{name}: enter its positive amount, constant across all runs.")
        sets = {}
        for run, info in self.runs.items():
            volume = info.get("volume")
            if volume is not None and (not math.isfinite(volume) or volume <= 0):
                errors.append(f"{run}: sample volume must be positive or blank.")
            if "cal" not in info.get("kind", "").lower():
                continue
            level = info.get("level", "").strip()
            if not level or level not in self.levels:
                errors.append(f"{run}: assign a calibration level with known amounts.")
                continue
            sets.setdefault(info.get("cal_set", ""), set()).add(level)
        if not sets:
            errors.append("Mark at least two known levels as calibration runs.")
        for cal_set, levels in sets.items():
            for t in analytes:
                values = [self.levels.get(lv, {}).get(t.compound) for lv in levels]
                valid = [v for v in values if v is not None and math.isfinite(v) and v >= 0]
                if len(valid) != len(values) or len(set(valid)) < 2:
                    errors.append(f"{t.compound}, set {cal_set or '(default)'}: enter at least two distinct "
                                  "non-negative amounts, covering every used calibration level.")
        if not self.confirmed:
            errors.append("Confirm the amounts, units and constant internal-standard additions.")
        return list(dict.fromkeys(errors))
