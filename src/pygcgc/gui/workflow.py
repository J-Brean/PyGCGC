"""The left-hand workflow panel.

One step at a time: pick a step (from the numbered buttons or the toolbar),
its settings appear below, and a single Go button runs it. When a step
finishes, a banner offers the next one.
"""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets


class WorkflowPanel(QtWidgets.QWidget):
    goRequested = QtCore.Signal(str)        # step key
    stepChanged = QtCore.Signal(str)
    progressChanged = QtCore.Signal()       # a tick or the suggested next step changed

    def __init__(self, parent=None):
        super().__init__(parent)
        self._steps: list[dict] = []
        self._by_key: dict[str, dict] = {}
        self._next_key: str | None = None

        self.title = QtWidgets.QLabel("")
        f = self.title.font()
        f.setPointSizeF(f.pointSizeF() * 1.25)
        f.setBold(True)
        self.title.setFont(f)
        self.title.setFont(QtGui.QFont("Georgia", self.title.font().pointSize(), QtGui.QFont.Weight.Bold))
        self.blurb = QtWidgets.QLabel("")
        self.blurb.setWordWrap(True)
        self.why = QtWidgets.QLabel("")
        self.why.setWordWrap(True)
        self.why.setStyleSheet("color: palette(mid);")
        self.stack = QtWidgets.QStackedWidget()

        self.banner = QtWidgets.QFrame()
        self.banner.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        self.banner_text = QtWidgets.QLabel("")
        self.banner_text.setWordWrap(True)
        self.banner_next = QtWidgets.QPushButton("")
        self.banner_next.clicked.connect(self._go_next)
        bl = QtWidgets.QVBoxLayout(self.banner)
        bl.setContentsMargins(6, 4, 6, 4)
        bl.addWidget(self.banner_text)
        bl.addWidget(self.banner_next)
        self.banner.setVisible(False)

        self.go_btn = QtWidgets.QPushButton("")
        self.go_btn.setMinimumHeight(34)
        f = self.go_btn.font()
        f.setBold(True)
        self.go_btn.setFont(f)
        self.go_btn.clicked.connect(lambda: self.goRequested.emit(self.current()))
        self.go_note = QtWidgets.QLabel("")
        self.go_note.setWordWrap(True)
        self.go_note.setStyleSheet("color: palette(mid);")
        self.go_note.setVisible(False)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(6)
        lay.addWidget(self.title)
        lay.addWidget(self.blurb)
        lay.addWidget(self.why)
        lay.addWidget(self.stack, 1)
        lay.addWidget(self.banner)
        lay.addWidget(self.go_note)
        lay.addWidget(self.go_btn)

    # ---------------------------------------------------------------- setup
    def add_step(self, key: str, label: str, title: str, blurb: str, why: str,
                 widget: QtWidgets.QWidget, go_label: str, scrolled: bool = True) -> None:
        """Steps are chosen from the toolbar (the window owns those buttons); this
        panel shows the open step's settings and its Go button. `label` names the
        step in the next-step banner, `title` is the heading, `blurb` says what the
        step does and `why` says why it sits at this point in the order."""
        n = len(self._steps)
        if scrolled:
            area = QtWidgets.QScrollArea()
            area.setWidgetResizable(True)
            area.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
            area.setWidget(widget)
            host = area
        else:
            host = widget
        self.stack.addWidget(host)
        step = dict(key=key, n=n + 1, label=label, title=title, blurb=blurb, why=why, host=host,
                    go_label=go_label, done=False, go_enabled=True, note="")
        self._steps.append(step)
        self._by_key[key] = step
        if n == 0:
            self.show_step(key, emit=False)
            self._next_key = key

    def next_key(self) -> str | None:
        return self._next_key

    def step_info(self, key: str) -> dict:
        return self._by_key[key]

    def guide_html(self) -> str:
        """Every step, what it does and why it sits where it does (Help → Workflow guide)."""
        rows = "".join(
            f"<p><b>{s['n']}. {s['label']}</b><br>{s['blurb']}<br>"
            f"<i>Why here:</i> {s['why']}<br><i>Go:</i> {s['go_label']}.</p>" for s in self._steps)
        return ("<p>Work along the toolbar from left to right. Each step's settings open in the "
                "<b>Tool</b> tab on the left; press <b>Go</b> to run it on every run. A banner then offers "
                "the next step.</p>" + rows)

    # ---------------------------------------------------------------- state
    def current(self) -> str:
        return self._steps[self.stack.currentIndex()]["key"] if self._steps else ""

    def keys(self) -> list[str]:
        return [s["key"] for s in self._steps]

    def is_done(self, key: str) -> bool:
        return self._by_key[key]["done"]

    def show_step(self, key: str, emit: bool = True) -> None:
        s = self._by_key[key]
        self.stack.setCurrentWidget(s["host"])
        self.title.setText(f"{s['n']}. {s['title']}")
        self.blurb.setText(s["blurb"])
        self.why.setText(f"Why here: {s['why']}")
        self.go_btn.setText(f"▶  {s['go_label']}")
        self._refresh_go(s)
        # The banner belongs to the step that just finished; drop it once the user moves on.
        if self.banner.property("for") not in (None, key):
            self.banner.setVisible(False)
        if emit:
            self.stepChanged.emit(key)

    def mark_done(self, key: str, message: str = "") -> None:
        """A step finished: tick it and offer the next one."""
        s = self._by_key[key]
        s["done"] = True
        nxt = self._following(key)
        self._next_key = nxt
        self.progressChanged.emit()
        text = f"✓ {message}" if message else f"✓ {s['label']} done"
        self.banner_text.setText(text)
        self.banner.setProperty("for", key)
        if nxt is None:
            self.banner_next.setVisible(False)
        else:
            n = self._by_key[nxt]
            self.banner_next.setText(f"Next: {n['label']}  →")
            self.banner_next.setVisible(True)
        self.banner.setVisible(True)

    def clear_done(self, key: str | None = None) -> None:
        for s in self._steps if key is None else [self._by_key[key]]:
            s["done"] = False
        first = next((s["key"] for s in self._steps if not s["done"]), None)
        self._next_key = first
        if key is None:
            self.banner.setVisible(False)
        self.progressChanged.emit()

    def set_go_enabled(self, key: str, enabled: bool, note: str = "") -> None:
        s = self._by_key[key]
        s["go_enabled"], s["note"] = enabled, note
        if key == self.current():
            self._refresh_go(s)

    def set_busy(self, key: str, busy: bool) -> None:
        s = self._by_key[key]
        s["busy"] = busy
        if key == self.current():
            self._refresh_go(s)

    # -------------------------------------------------------------- helpers
    def _following(self, key: str) -> str | None:
        keys = self.keys()
        after = keys[keys.index(key) + 1:] + keys[:keys.index(key)]
        return next((k for k in after if not self._by_key[k]["done"]), None)

    def _go_next(self) -> None:
        if self._next_key:
            self.show_step(self._next_key)

    def _refresh_go(self, s: dict) -> None:
        busy = s.get("busy", False)
        self.go_btn.setEnabled(s["go_enabled"] and not busy)
        self.go_btn.setText("Working…" if busy else f"▶  {s['go_label']}")
        note = "" if s["go_enabled"] else s["note"]
        self.go_note.setText(note)
        self.go_note.setVisible(bool(note))



class QuantStepPanel(QtWidgets.QWidget):
    """Settings side of the 🧪 step: what is loaded, and buttons to load the rest."""
    targetsRequested = QtCore.Signal()
    levelsRequested = QtCore.Signal()
    settingsRequested = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)
        targets = QtWidgets.QPushButton("Load target table…")
        targets.clicked.connect(self.targetsRequested)
        levels = QtWidgets.QPushButton("Load calibration levels…")
        levels.clicked.connect(self.levelsRequested)
        settings = QtWidgets.QPushButton("⚙️ Quantification settings…")
        settings.clicked.connect(self.settingsRequested)
        note = QtWidgets.QLabel("Results, review and calibration curves appear in the tabs and panels "
                                "at the bottom and centre once it has run.")
        note.setWordWrap(True)
        note.setStyleSheet("color: palette(mid);")
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(self.status)
        for b in (targets, levels, settings):
            lay.addWidget(b)
        lay.addWidget(note)
        lay.addStretch(1)

    def set_status(self, targets_ok: bool, why: str, levels_ok: bool, cleaned: bool) -> None:
        rows = [f"{'✓' if targets_ok else '✗'} Target table" + ("" if targets_ok else f": {why}"),
                f"{'✓' if levels_ok else '○'} Calibration levels" + ("" if levels_ok else " (needed for curves)"),
                f"{'✓' if cleaned else '○'} 🧹 Background applied" + ("" if cleaned else " (recommended)")]
        self.status.setText("\n".join(rows))


class IdentifyStepPanel(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        text = QtWidgets.QLabel(
            "Searches every 🔍 blob on the current run against an MSP library, using the retention "
            "index from the ⛓️ ladder if one exists.\n\nGo opens the identifier window. Load a library "
            "there, then press 🏷️ Identify blobs.")
        text.setWordWrap(True)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(text)
        lay.addStretch(1)
