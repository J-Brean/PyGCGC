"""📚 Naming a peak: the label shown on the 2D image and the dialog for suggesting another name.

A peak's label is its decided name if the user gave it one (accepted, chosen from the candidates, or typed),
otherwise the matcher's best candidate. The dialog checks a typed name against the library the way
`core/naming.py` describes (capitals, commas, dashes and spaces do not make a new compound) and says plainly
what will be saved. The result is stored as a `named` decision in `MatchState`, with the same content
signature as every other decision, so it only stands while the evidence it was made on is unchanged.
"""

from __future__ import annotations

import numpy as np
from PySide6 import QtCore, QtWidgets

from ..core.naming import check_name, format_cas
from .match_state import NAMING_ACTIONS
from .theme import fit_to_screen

LEVEL_RGB = {"confident": (0, 158, 115), "probable": (86, 180, 233), "tentative": (230, 159, 0), "unknown": (170, 170, 170)}
LEVEL_SYMBOL = {"confident": "●", "probable": "◐", "tentative": "○", "unknown": "✕"}


def identity_of(match, decision) -> tuple[str, str, str]:
    """`(name, formula, CAS)` the peak goes by, any of them possibly empty.

    A decision that names the peak wins. Its formula is the one saved with it, or, for a decision made before
    formulae were saved, that of the candidate with the same name; a name that is in no library has none, and
    the best candidate's formula is never lent to a name the user typed. A rejected or unknown peak has no
    identity. Otherwise it is the best candidate's."""
    action = (decision or {}).get("action")
    if action in ("rejected", "unknown"):
        return "", "", ""
    if action in NAMING_ACTIONS and decision.get("name"):
        name = decision["name"]
        formula, cas = decision.get("formula"), decision.get("cas", "")
        if formula is None:
            same = next((c for c in match.candidates if c.name == name), None)
            formula = same.formula if same is not None else ""
            cas = cas or (same.cas if same is not None else "")
        return name, formula or "", cas or ""
    best = match.best
    return (best.name, best.formula, best.cas) if best is not None else ("", "", "")


def label_for(match, decision) -> tuple[str, tuple]:
    """`(text, rgb)` for a peak's label on the image: its name, how it got it, and how sure the match is."""
    best = match.best
    named = decision is not None and decision.get("action") in NAMING_ACTIONS and decision.get("name")
    if decision is not None and decision.get("action") == "rejected":
        name, source = "(rejected)", ""
    elif named:
        name = decision["name"]
        source = {"named": "named by you", "chosen": "chosen by you", "accepted": "accepted by you"}[decision["action"]]
    elif best is not None:
        name, source = best.name, ""
    else:
        name, source = "(no match)", ""
    second = f"{LEVEL_SYMBOL[match.level]} {match.level}" if not named else source
    if best is not None and not named:
        second += f", MF {best.mf:.0f}"
    if match.flags and not named:
        second += f", {len(match.flags)} flag{'s' if len(match.flags) > 1 else ''}"
    return f"{name}\n{second}", (LEVEL_RGB["confident"] if named else LEVEL_RGB[match.level])


class PeakNameDialog(QtWidgets.QDialog):
    """Suggest another name for one peak, with the checks of `core.naming.check_name` shown as you type."""

    def __init__(self, library, match, run_name: str, share: float, taken, current: dict | None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Name this peak")
        fit_to_screen(self, 560, 520)
        self.library, self.match, self.taken = library, match, list(taken)
        self.check = None
        self.cleared = False
        best = match.best
        head = QtWidgets.QLabel(
            f"<b>{run_name}</b>: peak at t₁ {match.query.t1 / 60:.3f} min, t₂ {match.query.t2:.2f} s"
            + (f", {100 * share:.1f} % of the run's signal" if np.isfinite(share) else "") + "<br>"
            + (f"The library's best match is <b>{best.name}</b> ({match.level}, MF {best.mf:.0f})." if best else
               "The library has no candidate for this peak."))
        head.setWordWrap(True)
        self.edit = QtWidgets.QLineEdit()
        self.edit.setPlaceholderText("Type a name")
        self.edit.setText((current or {}).get("name") or (best.name if best else ""))
        self.edit.textChanged.connect(self._changed)
        self.report = QtWidgets.QLabel("")
        self.report.setWordWrap(True)
        self.report.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.report.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop)
        self.report.setMinimumHeight(120)
        self.suggest = QtWidgets.QListWidget()
        self.suggest.setMaximumHeight(110)
        self.suggest.itemClicked.connect(lambda item: self.edit.setText(item.data(QtCore.Qt.ItemDataRole.UserRole)))
        self.suggest.itemDoubleClicked.connect(lambda item: (self.edit.setText(item.data(QtCore.Qt.ItemDataRole.UserRole)),
                                                             self.accept()))
        for c in match.candidates:
            item = QtWidgets.QListWidgetItem(f"{c.name}   (MF {c.mf:.0f}" + (f", {c.formula}" if c.formula else "") + ")")
            item.setData(QtCore.Qt.ItemDataRole.UserRole, c.name)
            self.suggest.addItem(item)
        self.buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Save
                                                  | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        clear = self.buttons.addButton("Remove my name", QtWidgets.QDialogButtonBox.ButtonRole.DestructiveRole)
        clear.setEnabled(bool(current and current.get("action") in NAMING_ACTIONS))
        clear.clicked.connect(self._clear)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addWidget(head)
        lay.addWidget(QtWidgets.QLabel("Name"))
        lay.addWidget(self.edit)
        lay.addWidget(self.report, 1)
        lay.addWidget(QtWidgets.QLabel("Or click a candidate"))
        lay.addWidget(self.suggest)
        lay.addWidget(self.buttons)
        self._changed()

    def _changed(self, *_a) -> None:
        self.check = check_name(self.edit.text(), self.library, [c.name for c in self.match.candidates], self.taken)
        c = self.check
        save = self.buttons.button(QtWidgets.QDialogButtonBox.StandardButton.Save)
        save.setEnabled(c.status != "empty")
        if c.status == "empty":
            self.report.setText("<span style='color:#D55E00'>Type a name.</span>")
            return
        saved = f"<b>{c.name}</b>" + (f" (CAS {format_cas(c.cas)}" + (f", {c.formula}" if c.formula else "") + ")" if c.cas else "")
        colour = {"exact": "#009E73", "respelled": "#0072B2", "new": "#D55E00"}[c.status]
        words = {"exact": "In the library", "respelled": "Library spelling used", "new": "Not in the library"}[c.status]
        lines = [f"<span style='color:{colour}'><b>{words}</b></span>", f"Will be saved as {saved}."] + c.notes
        if c.suggestions:
            lines.append("Did you mean: " + "; ".join(c.suggestions) + "?")
        self.report.setText("<br>".join(lines))

    def _clear(self) -> None:
        self.cleared = True
        self.accept()

    def result_name(self) -> tuple[str, str, str, str]:
        """`(name to save, CAS, what was typed, formula)`; CAS and formula are empty for a name no library holds."""
        return self.check.name, self.check.cas, self.edit.text().strip(), self.check.formula


def rename_peak(parent, project, index: int, label: str) -> bool:
    """Open the naming dialog for one peak of run `index` and store the outcome. True if anything changed."""
    state = project.matcher
    found = state.current_for(index)
    if found is None:
        QtWidgets.QMessageBox.information(parent, "Name this peak",
                                          "Match this run's peaks to a library first (Assign → Library match). A name belongs "
                                          "to the peak as it was matched, so it can be tied to that evidence.")
        return False
    match = next((m for m in found.matches if m.query.label == label), None)
    if match is None:
        return False
    run = project.entries[index].name
    decisions = state.decisions_for(index)
    taken = []
    for other in found.matches:
        d = decisions.get(other.query.label)
        if other.query.label != label and d and d.get("action") in NAMING_ACTIONS and d.get("name"):
            taken.append((f"at {other.query.t1 / 60:.2f} min", d["name"]))
    total = sum(max(m.query.volume, 0.0) for m in found.matches)
    share = match.query.volume / total if total > 0 else float("nan")
    dialog = PeakNameDialog(state.library, match, run, share, taken, decisions.get(label), parent)
    if not dialog.exec():
        return False
    if dialog.cleared:
        state.clear_decision(index, label)
        return True
    name, cas, typed, formula = dialog.result_name()
    state.decide(index, label, "named", name, cas, typed=typed, formula=formula)
    return True
