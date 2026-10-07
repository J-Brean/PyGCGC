"""Spectrum panel.

Shows the spectrum under the cursor click (or summed over the region box).
Click an ion to image it; Shift+click adds it to the current sum.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

N_LABELS = 8


class SpectrumView(QtWidgets.QWidget):
    ionClicked = QtCore.Signal(int, bool)     # m/z, add to current expression

    def __init__(self, parent=None):
        super().__init__(parent)
        self.mz = np.array([])
        self.intensity = np.array([])
        self.title = ""
        self.relative = True

        self.header = QtWidgets.QLabel("Click the chromatogram to see a spectrum")
        self.header.setWordWrap(True)
        self.rel_box = QtWidgets.QCheckBox("Relative (%)")
        self.rel_box.setChecked(True)
        self.rel_box.toggled.connect(self._relative_toggled)
        copy = QtWidgets.QPushButton("Copy")
        copy.setToolTip("Copy m/z and intensity as tab-separated text")
        copy.clicked.connect(self.copy)

        self.plot = pg.PlotWidget()
        self.plot.setMenuEnabled(False)
        self.plot.hideButtons()
        self.plot.setLabel("bottom", "m/z")
        self.plot.setMouseEnabled(y=False)
        # Sticks drawn as line pairs always render at least one pixel wide.
        self.sticks = pg.PlotDataItem([], [], connect="pairs",
                                      pen=pg.mkPen((209, 149, 224), width=2))
        self.plot.addItem(self.sticks)
        self.labels: list[pg.TextItem] = []
        self.plot.scene().sigMouseClicked.connect(self._clicked)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.addWidget(self.header)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.rel_box)
        row.addStretch(1)
        row.addWidget(copy)
        lay.addLayout(row)
        lay.addWidget(self.plot, 1)
        hint = QtWidgets.QLabel("Click an ion to image it; Shift+click to add it.")
        hint.setWordWrap(True)
        hint.setStyleSheet("color: palette(mid);")
        lay.addWidget(hint)

    def set_background(self, colour) -> None:
        self.plot.setBackground(colour)

    def set_spectrum(self, mz: np.ndarray, intensity: np.ndarray, title: str) -> None:
        keep = intensity > 0
        self.mz = np.asarray(mz)[keep]
        self.intensity = np.asarray(intensity, dtype=float)[keep]
        self.title = title
        self._draw()

    def clear(self) -> None:
        self.set_spectrum(np.array([]), np.array([]), "")
        self.header.setText("Click the chromatogram to see a spectrum")

    def _relative_toggled(self, on: bool) -> None:
        self.relative = on
        self._draw()

    def _draw(self) -> None:
        for t in self.labels:
            self.plot.removeItem(t)
        self.labels = []
        if self.mz.size == 0:
            self.sticks.setData([], [])
            if self.title:
                self.header.setText(f"{self.title}\nNo signal here")
            return
        y = self.intensity
        top = y.max()
        if self.relative:
            y = 100.0 * y / top
        xs = np.repeat(self.mz.astype(float), 2)
        ys = np.column_stack([np.zeros_like(y), y]).ravel()
        self.sticks.setData(xs, ys)
        self.plot.setLabel("left", "relative intensity (%)" if self.relative else "intensity")
        # Label the biggest peaks, skipping any within 3 m/z of a bigger labelled one.
        chosen: list[int] = []
        for i in np.argsort(y)[::-1]:
            if len(chosen) == N_LABELS or y[i] < 0.02 * y.max():
                break
            if all(abs(self.mz[i] - self.mz[j]) > 3 for j in chosen):
                chosen.append(int(i))
        for i in chosen:
            t = pg.TextItem(str(int(self.mz[i])), anchor=(0.5, 1.0), color=(230, 230, 235))
            t.setPos(float(self.mz[i]), float(y[i]))
            self.plot.addItem(t)
            self.labels.append(t)
        self.plot.setXRange(max(self.mz.min() - 5, 0), self.mz.max() + 5, padding=0)
        self.plot.setYRange(0, float(y.max()) * 1.15, padding=0)
        base = int(self.mz[np.argmax(self.intensity)])
        self.header.setText(f"{self.title}\nBase peak m/z {base}, total {self.intensity.sum():,.0f}")

    def _clicked(self, ev) -> None:
        if self.mz.size == 0 or not self.plot.sceneBoundingRect().contains(ev.scenePos()):
            return
        x = self.plot.plotItem.vb.mapSceneToView(ev.scenePos()).x()
        i = int(np.argmin(np.abs(self.mz - x)))
        if abs(self.mz[i] - x) <= 1.5:
            add = bool(ev.modifiers() & QtCore.Qt.KeyboardModifier.ShiftModifier)
            self.ionClicked.emit(int(self.mz[i]), add)

    def copy(self) -> None:
        lines = ["m/z\tintensity"] + [f"{int(m)}\t{v:.6g}" for m, v in zip(self.mz, self.intensity)]
        QtGui.QGuiApplication.clipboard().setText("\n".join(lines))
