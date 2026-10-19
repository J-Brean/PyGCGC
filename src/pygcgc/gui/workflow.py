"""The left-hand workflow panel.

One step at a time: pick a step (from the numbered buttons or the toolbar),
its settings appear below, and a single Go button runs it. When a step
finishes, a banner offers the next one.
"""

from __future__ import annotations

from PySide6 import QtCore, QtGui, QtWidgets

from .icons import plain
from .theme import fade_in


class WorkflowPanel(QtWidgets.QWidget):
    goRequested = QtCore.Signal(str)        # step key
    stepChanged = QtCore.Signal(str)
    progressChanged = QtCore.Signal()       # a tick or the suggested next step changed
    prerequisiteRequested = QtCore.Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._steps: list[dict] = []
        self._by_key: dict[str, dict] = {}
        self._next_key: str | None = None

        # A sub-step (a tool opened from inside another step) shows a way back to its parent.
        self.back_btn = QtWidgets.QPushButton("")
        self.back_btn.setFlat(True)
        self.back_btn.setStyleSheet("text-align: left; padding: 0;")
        self.back_btn.clicked.connect(self._go_back)
        self.back_btn.setVisible(False)
        self.title = QtWidgets.QLabel("")
        self.title.setWordWrap(True)
        f = self.title.font()
        f.setPointSizeF(f.pointSizeF() * 1.9)
        f.setBold(True)
        self.title.setFont(f)
        self.blurb = QtWidgets.QLabel("")
        self.blurb.setWordWrap(True)
        self.why = QtWidgets.QLabel("")
        self.why.setWordWrap(True)
        self.why.setStyleSheet("color: palette(mid);")
        self.blurb.setStyleSheet("font-size: 10pt;")
        self.stack = QtWidgets.QStackedWidget()
        self.stack.currentChanged.connect(lambda _i: fade_in(self.stack.currentWidget()))

        self.banner = QtWidgets.QFrame()
        self.banner.setObjectName("banner")
        self.banner_text = QtWidgets.QLabel("")
        self.banner_text.setWordWrap(True)
        self.banner_next = QtWidgets.QPushButton("")
        self.banner_next.clicked.connect(self._go_next)
        bl = QtWidgets.QVBoxLayout(self.banner)
        bl.setContentsMargins(12, 10, 12, 10)
        bl.addWidget(self.banner_text)
        bl.addWidget(self.banner_next)
        self.banner.setVisible(False)

        self.go_btn = QtWidgets.QPushButton("")
        self.go_btn.setMinimumHeight(34)
        self.go_btn.setObjectName("go")
        f = self.go_btn.font()
        f.setBold(True)
        self.go_btn.setFont(f)
        self.go_btn.clicked.connect(lambda: self.goRequested.emit(self.current()))
        self.go_note = QtWidgets.QLabel("")
        self.go_note.setWordWrap(True)
        self.go_note.setStyleSheet("color: palette(mid);")
        self.go_note.setVisible(False)
        self.prerequisite_btn = QtWidgets.QPushButton()
        self.prerequisite_btn.setVisible(False)
        self.prerequisite_btn.clicked.connect(lambda: self.prerequisiteRequested.emit(
            self._by_key[self.current()].get("prerequisite", "")))

        self.skip_btn = QtWidgets.QPushButton("Skip alignment")
        self.skip_btn.clicked.connect(self._skip_current)
        self.skip_btn.setVisible(False)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 12, 0, 0)
        lay.setSpacing(8)
        # Some step panels ask for 900 px or more; do not let them widen the workspace.
        lay.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetNoConstraint)
        lay.addWidget(self.back_btn)
        lay.addWidget(self.title)
        lay.addWidget(self.blurb)
        lay.addWidget(self.why)
        lay.addWidget(self.stack, 1)
        lay.addWidget(self.banner)
        lay.addWidget(self.go_note)
        lay.addWidget(self.prerequisite_btn)
        lay.addWidget(self.skip_btn)
        lay.addWidget(self.go_btn)

    def minimumSizeHint(self) -> QtCore.QSize:
        return QtCore.QSize(320, super().minimumSizeHint().height())

    # ---------------------------------------------------------------- setup
    def add_step(self, key: str, label: str, title: str, blurb: str, why: str,
                 widget: QtWidgets.QWidget, go_label: str, scrolled: bool = True,
                 parent: str | None = None, optional: bool = False) -> None:
        """Steps are chosen from the toolbar (the window owns those buttons); this
        panel shows the open step's settings and its Go button. `label` names the
        step in the next-step banner, `title` is the heading, `blurb` says what the
        step does and `why` says why it sits at this point in the order.
        A step with a `parent` is a sub-menu of that step: it has no number or toolbar
        button, is reached from the parent's panel, and is never offered as the next step."""
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
        number = None if parent else sum(1 for s in self._steps if not s["parent"]) + 1
        step = dict(key=key, n=number, label=label, title=title, blurb=blurb, why=why, host=host,
                    go_label=go_label, done=False, skipped=False, optional=optional,
                    go_enabled=True, note="", parent=parent)
        self._steps.append(step)
        self._by_key[key] = step
        if n == 0:
            self.show_step(key, emit=False)
            self._next_key = key

    def next_key(self) -> str | None:
        return self._next_key

    def set_go_label(self, key: str, label: str) -> None:
        """Keep the primary action clear when a step starts or finishes picking."""
        step = self._by_key[key]
        step["go_label"] = label
        if self.current() == key:
            self._refresh_go(step)

    def set_order(self, keys: list[str]) -> None:
        """Order the workflow independently of the settings widgets' stack indices."""
        if len(keys) != len(self._steps) or set(keys) != set(self._by_key):
            raise ValueError("Workflow order must include every step exactly once")
        self._steps = [self._by_key[k] for k in keys]
        number = 0
        for step in self._steps:
            if not step["parent"]:
                number += 1
            step["n"] = None if step["parent"] else number

    def step_info(self, key: str) -> dict:
        return self._by_key[key]

    def guide_html(self) -> str:
        """Every step, what it does and why it sits where it does (Help → Workflow guide)."""
        rows = "".join(
            f"<p><b>{s['n'] if s['n'] else '›'}. {s['label']}</b><br>{s['blurb']}<br>"
            f"{s['why']}<br><i>Go:</i> {s['go_label']}.</p>" for s in self._steps)
        return ("<p>Work along the toolbar from left to right. Each step's settings open in the "
                "<b>Tool</b> tab on the left; press <b>Go</b> to run it on every run. A banner then offers "
                "the next step.</p>" + rows)

    # ---------------------------------------------------------------- state
    def current(self) -> str:
        host = self.stack.currentWidget()
        return next((s["key"] for s in self._steps if s["host"] is host), "")

    def keys(self) -> list[str]:
        return [s["key"] for s in self._steps]

    def is_done(self, key: str) -> bool:
        return self._by_key[key]["done"]

    def is_skipped(self, key: str) -> bool:
        return self._by_key[key]["skipped"]

    def show_step(self, key: str, emit: bool = True) -> None:
        s = self._by_key[key]
        self.stack.setCurrentWidget(s["host"])
        self.title.setText(f"{s['n']}. {s['title']}" if s["n"] else s["title"])
        parent = s["parent"]
        self.back_btn.setVisible(bool(parent))
        if parent:
            self.back_btn.setText(f"‹ Back to {plain(self._by_key[parent]['label'])}")
        self.blurb.setText(s["blurb"])
        self.why.setText(s["why"])
        self.go_btn.setText(f"▶  {s['go_label']}")
        self._refresh_go(s)
        self.skip_btn.setVisible(s["optional"])
        # The banner belongs to the step that just finished; drop it once the user moves on.
        if self.banner.property("for") not in (None, key):
            self.banner.setVisible(False)
        if emit:
            self.stepChanged.emit(key)

    def mark_done(self, key: str, message: str = "") -> None:
        """A step finished: tick it and offer the next one."""
        s = self._by_key[key]
        s["done"] = True
        s["skipped"] = False
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
            self.banner_next.setText(f"Next: {plain(n['label'])}  →")
            self.banner_next.setVisible(True)
        self.banner.setVisible(True)

    def clear_done(self, key: str | None = None) -> None:
        for s in self._steps if key is None else [self._by_key[key]]:
            s["done"] = False
            s["skipped"] = False
        first = next((s["key"] for s in self._steps
                      if not s["done"] and not s["skipped"] and not s["parent"]), None)
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

    def set_guidance(self, key: str, message: str, action: str = "", label: str = "") -> None:
        s = self._by_key[key]
        s['guidance'], s['prerequisite'], s['prerequisite_label'] = message, action, label
        if key == self.current():
            self._refresh_go(s)

    # -------------------------------------------------------------- helpers
    def _following(self, key: str) -> str | None:
        # A sub-step hands on from its parent's position, and sub-steps are never offered.
        key = self._by_key[key]["parent"] or key
        keys = [k for k in self.keys() if not self._by_key[k]["parent"]]
        after = keys[keys.index(key) + 1:]
        return next((k for k in after
                     if not self._by_key[k]["done"] and not self._by_key[k]["skipped"]), None)

    def _skip_current(self) -> None:
        s = self._by_key[self.current()]
        if not s["optional"]:
            return
        s["skipped"] = True
        s["done"] = False
        self._next_key = self._following(s["key"])
        self.progressChanged.emit()
        if self._next_key:
            self.show_step(self._next_key)

    def _go_back(self) -> None:
        parent = self._by_key[self.current()]["parent"]
        if parent:
            self.show_step(parent)

    def _go_next(self) -> None:
        if self._next_key:
            self.show_step(self._next_key)

    def _refresh_go(self, s: dict) -> None:
        busy = s.get("busy", False)
        self.go_btn.setEnabled(s["go_enabled"] and not busy)
        self.go_btn.setText("Working…" if busy else f"▶  {s['go_label']}")
        note = s.get('guidance') or ("" if s["go_enabled"] else s["note"])
        self.go_note.setText(note)
        self.go_note.setVisible(bool(note))
        self.prerequisite_btn.setText(s.get('prerequisite_label', ''))
        self.prerequisite_btn.setVisible(bool(s.get('prerequisite')) and not busy)



class ToolChoicesPanel(QtWidgets.QWidget):
    """A small set of named routes within a workflow step."""
    chosen = QtCore.Signal(str)

    def __init__(self, choices, note="", parent=None):
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        self.buttons = {}
        for key, label in choices:
            button = QtWidgets.QPushButton(label)
            button.clicked.connect(lambda _checked=False, k=key: self.chosen.emit(k))
            layout.addWidget(button)
            self.buttons[key] = button
        self.status = QtWidgets.QLabel(note)
        self.status.setWordWrap(True)
        layout.addWidget(self.status)
        layout.addStretch(1)


class QuantStepPanel(QtWidgets.QWidget):
    """Settings side of the 🧪 step: what is loaded, and buttons to load the rest."""
    targetsRequested = QtCore.Signal()
    levelsRequested = QtCore.Signal()
    settingsRequested = QtCore.Signal()
    calibrateRequested = QtCore.Signal()
    setupRequested = QtCore.Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setup_btn = QtWidgets.QPushButton("Set up calibration and units…")
        self.setup_btn.clicked.connect(self.setupRequested)
        self.unit_label = QtWidgets.QLabel()
        self.unit_label.setWordWrap(True)
        # The guided route: curves built from the standards in their own window, then applied
        # to the other runs. The rest of this panel is the table-driven route.
        calibrate = QtWidgets.QPushButton("📈 Build curves from my standards…")
        calibrate.setStyleSheet("QPushButton { font-weight: bold; padding: 8px; }")
        calibrate.setToolTip("Select your calibration runs in the Runs table first. Compounds are suggested or "
                             "loaded, curves can be edited point by point, then applied to the other runs.")
        calibrate.clicked.connect(self.calibrateRequested)
        guided = QtWidgets.QLabel("Set up compounds, calibration runs, known amounts and units above. "
                                  "Then press Quantify. Existing CSV inputs and the advanced curve builder "
                                  "are available in Advanced calibration tools.")
        guided.setWordWrap(True)
        guided.setStyleSheet("color: palette(mid);")
        self.status = QtWidgets.QLabel("")
        self.status.setWordWrap(True)
        targets = QtWidgets.QPushButton("Load target table…")
        targets.clicked.connect(self.targetsRequested)
        levels = QtWidgets.QPushButton("Load calibration levels…")
        levels.clicked.connect(self.levelsRequested)
        settings = QtWidgets.QPushButton("⚙️ Quantification settings…")
        settings.clicked.connect(self.settingsRequested)
        note = QtWidgets.QLabel("Once it has run, results are in the Results tab on the left, and the "
                                "Review and Calibration tabs above the chromatogram.")
        note.setWordWrap(True)
        note.setStyleSheet("color: palette(mid);")
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addWidget(self.setup_btn)
        lay.addWidget(self.unit_label)
        lay.addWidget(guided)
        lay.addSpacing(10)
        lay.addWidget(self.status)
        advanced = QtWidgets.QToolButton(self)
        advanced.setText("Advanced calibration tools")
        advanced.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QtWidgets.QMenu(advanced)
        for button in (targets, levels, calibrate, settings):
            button.setParent(self)
            button.hide()
            menu.addAction(button.text(), button.click)
        advanced.setMenu(menu)
        lay.addWidget(advanced)
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
