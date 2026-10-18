"""A small, local animation explaining joint multi-component peak fitting.

This is a scripted illustration with synthetic peaks, not an optimiser trace or
a chemical identification. Both views use the same profiles at every frame.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets


class FittingHelpDialog(QtWidgets.QDialog):
    DURATION_MS = 16000
    CAPTIONS = (
        ("1. One signal can hide overlapping peaks",
         "The map shows two synthetic peaks added together. The 1D trace sums the same map "
         "over the second retention time, so their overlap is harder to see."),
        ("2. Propose separate contributions",
         "Ion-by-ion peak positions provide starting guesses. Blue and pink are proposed "
         "components; orange is their sum plus the local background."),
        ("3. Fit positions, widths, tilt and spectra together",
         "Each component has independent widths and a non-negative mass spectrum. Tilt lets "
         "its second retention time drift across modulations. All selected ions constrain the fit."),
        ("4. Keep extra components only when evidence supports them",
         "The automatic search starts with one peak. Extra peaks must improve BIC and held-out "
         "prediction, carry enough signal, and have distinct spectra or clear spatial separation."),
    )

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("How overlapping peaks are fitted")
        self.resize(940, 610)
        self.setMinimumSize(700, 500)
        self.position_ms = 0
        self._x = np.linspace(0, 80, 161)
        self._y = np.linspace(0, 4, 101)
        self._xx, self._yy = np.meshgrid(self._x, self._y)
        # (first time, second time, first width, second width, correlation, height)
        self._truth = np.array([[35., 1.8, 9., .40, -.55, 1.],
                                [45., 2.15, 6., .25, .45, .9]])
        self._guess = np.array([[29., 1.6, 11., .50, 0., .7],
                                [50., 2.4, 9., .40, 0., .6]])
        self._observed = .025 + sum(self._profile(p) for p in self._truth)

        title = QtWidgets.QLabel("From an overlapping signal to separate contributions")
        font = title.font()
        font.setPointSize(font.pointSize() + 3)
        font.setBold(True)
        title.setFont(font)
        self.heading = QtWidgets.QLabel()
        self.heading.setStyleSheet("font-weight: bold;")
        self.caption = QtWidgets.QLabel()
        self.caption.setWordWrap(True)
        self.caption.setMinimumHeight(62)

        self.trace = pg.PlotWidget(title="1D: projection over the second retention time")
        self.map = pg.PlotWidget(title="2D: the same signal and fitted contours")
        self.trace.setLabel("bottom", "First retention time (illustrative seconds)")
        self.trace.setLabel("left", "Summed signal (arbitrary units)")
        self.map.setLabel("bottom", "First retention time (illustrative seconds)")
        self.map.setLabel("left", "Second retention time (s)")
        self.trace.setXRange(0, 80, padding=0)
        self.trace.setYRange(0, self._observed.sum(axis=0).max() * 1.25, padding=0)
        self.map.setRange(xRange=(0, 80), yRange=(0, 4), padding=0)
        for plot in (self.trace, self.map):
            plot.setMenuEnabled(False)
            plot.hideButtons()
            plot.setMouseEnabled(x=False, y=False)
        self.image = pg.ImageItem(axisOrder="row-major")
        self.image.setImage(self._observed, autoLevels=False, levels=(0, self._observed.max()))
        self.image.setRect(QtCore.QRectF(-.25, -.02, 80.5, 4.04))
        self.map.addItem(self.image)
        self.observed_curve = self.trace.plot(self._x, self._observed.sum(axis=0))
        self.sum_curve = self.trace.plot()
        self.component_curves = [self.trace.plot(), self.trace.plot()]
        self.contours = [[self.map.plot(), self.map.plot()] for _ in range(2)]
        self.key = QtWidgets.QLabel()
        self.key.setWordWrap(True)
        self.note = QtWidgets.QLabel(
            "Synthetic illustration; motion is scripted. A fitted component is a proposed contribution, "
            "not a compound identification. Unresolved mixtures still need review.")
        self.note.setWordWrap(True)

        self.play_button = QtWidgets.QPushButton("Pause")
        self.play_button.clicked.connect(self.toggle_playback)
        self.replay_button = QtWidgets.QPushButton("Replay")
        self.replay_button.clicked.connect(self.replay)
        self.slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.slider.setRange(0, self.DURATION_MS)
        self.slider.setAccessibleName("Animation position")
        self.slider.valueChanged.connect(self.set_position)
        self.time_label = QtWidgets.QLabel()
        self.time_label.setMinimumWidth(68)
        close_button = QtWidgets.QPushButton("Close")
        close_button.clicked.connect(self.close)
        controls = QtWidgets.QHBoxLayout()
        for widget in (self.play_button, self.replay_button, self.slider, self.time_label, close_button):
            controls.addWidget(widget, 1 if widget is self.slider else 0)
        plots = QtWidgets.QHBoxLayout()
        plots.addWidget(self.trace, 1)
        plots.addWidget(self.map, 1)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(title)
        layout.addWidget(self.heading)
        layout.addWidget(self.caption)
        layout.addLayout(plots, 1)
        layout.addWidget(self.key)
        layout.addWidget(self.note)
        layout.addLayout(controls)
        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(40)
        self.timer.timeout.connect(self._tick)
        self._restyle()
        self.set_position(0)

    def _profile(self, p):
        x, y, sx, sy, rho, height = p
        dx, dy = (self._xx-x)/sx, (self._yy-y)/sy
        return height * np.exp(-.5*(dx*dx + dy*dy - 2*rho*dx*dy)/(1-rho*rho))

    def _restyle(self):
        # QDialog's stylesheet can map its local Base role onto Window. Use the
        # application palette so these plots match the other themed plots.
        palette = QtWidgets.QApplication.palette()
        dark = palette.color(QtGui.QPalette.ColorRole.Window).lightness() < 128
        colours = ("#28d7e8", "#ff79c6") if dark else ("#007fa3", "#c02680")
        foreground = palette.color(QtGui.QPalette.ColorRole.WindowText)
        base = palette.color(QtGui.QPalette.ColorRole.Base)
        for plot in (self.trace, self.map):
            plot.setBackground(base)
            for axis in ("bottom", "left"):
                plot.getAxis(axis).setPen(foreground)
                plot.getAxis(axis).setTextPen(foreground)
        self.trace.setTitle("1D: projection over the second retention time", color=foreground.name())
        self.map.setTitle("2D: the same signal and fitted contours", color=foreground.name())
        self.image.setColorMap(pg.ColorMap([0, 1], [base, foreground]))
        self.observed_curve.setPen(pg.mkPen(foreground, width=2))
        self.sum_curve.setPen(pg.mkPen("#e69824" if dark else "#a85b00", width=2,
                                      style=QtCore.Qt.PenStyle.DashLine))
        for colour, curve, contours in zip(colours, self.component_curves, self.contours):
            curve.setPen(pg.mkPen(colour, width=2.5))
            for i, contour in enumerate(contours):
                contour.setPen(pg.mkPen(colour, width=2.5,
                                        style=(QtCore.Qt.PenStyle.SolidLine if i == 0
                                               else QtCore.Qt.PenStyle.DotLine)))
        self.key.setText(
            f'<span style="color:{colours[0]}">Blue: component 1</span> &nbsp; '
            f'<span style="color:{colours[1]}">Pink: component 2</span> &nbsp; '
            'Solid neutral: observed signal &nbsp; Orange dashed: fitted sum + background<br>'
            '2D contours: 10% (solid) and 50% (dotted) of each component\'s maximum; overlap is allowed.')

    def set_position(self, milliseconds: int) -> None:
        """Draw matching projections and contours at a chosen illustration time."""
        self.position_ms = int(np.clip(milliseconds, 0, self.DURATION_MS))
        self.slider.blockSignals(True)
        self.slider.setValue(self.position_ms)
        self.slider.blockSignals(False)
        stage = min(self.position_ms // 4000, 3)
        self.heading.setText(self.CAPTIONS[stage][0])
        self.caption.setText(self.CAPTIONS[stage][1])
        self.time_label.setText(f"{self.position_ms/1000:.0f} / 16 s")
        # The guesses become visible first, then converge smoothly to the truth.
        fraction = float(np.clip((self.position_ms-6500)/5500, 0, 1))
        fraction = fraction*fraction*(3-2*fraction)
        parameters = self._guess*(1-fraction) + self._truth*fraction
        profiles = [self._profile(p) for p in parameters]
        self.sum_curve.setData(self._x, (.025 + sum(profiles)).sum(axis=0))
        self.sum_curve.setVisible(stage > 0)
        theta = np.linspace(0, 2*np.pi, 101)
        for p, profile, curve, contours in zip(parameters, profiles, self.component_curves, self.contours):
            curve.setData(self._x, profile.sum(axis=0))
            curve.setVisible(stage > 0)
            x, y, sx, sy, rho, _height = p
            for contour, level in zip(contours, (.10, .50)):
                radius = np.sqrt(-2*np.log(level))
                contour.setData(x + radius*sx*np.cos(theta),
                                y + radius*sy*(rho*np.cos(theta)+np.sqrt(1-rho*rho)*np.sin(theta)))
                contour.setVisible(stage > 0)
        if self.position_ms == self.DURATION_MS:
            self.timer.stop()
            self.play_button.setText("Play")

    def _tick(self):
        self.set_position(self.position_ms + self.timer.interval())

    def toggle_playback(self):
        if self.timer.isActive():
            self.timer.stop()
            self.play_button.setText("Play")
        else:
            if self.position_ms == self.DURATION_MS:
                self.set_position(0)
            self.timer.start()
            self.play_button.setText("Pause")

    def replay(self):
        self.set_position(0)
        if self.isVisible():
            self.timer.start()
            self.play_button.setText("Pause")

    def showEvent(self, event):
        super().showEvent(event)
        self._restyle()
        self.replay()

    def hideEvent(self, event):
        self.timer.stop()
        self.play_button.setText("Play")
        super().hideEvent(event)

    def changeEvent(self, event):
        super().changeEvent(event)
        if event.type() in (QtCore.QEvent.Type.PaletteChange,
                            QtCore.QEvent.Type.ApplicationPaletteChange) and hasattr(self, "key"):
            self._restyle()


def open_fitting_help(owner: QtWidgets.QWidget) -> FittingHelpDialog:
    """Reuse one modeless explanation per host, without starting a fitting job."""
    dialog = getattr(owner, "_fitting_help", None)
    if dialog is None:
        dialog = FittingHelpDialog(owner)
        owner._fitting_help = dialog
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
    return dialog


def fitting_help_button(owner: QtWidgets.QWidget) -> QtWidgets.QToolButton:
    button = QtWidgets.QToolButton(owner)
    button.setText("How it works")
    button.setIcon(button.style().standardIcon(QtWidgets.QStyle.StandardPixmap.SP_MediaPlay))
    button.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
    button.setAccessibleName("Play an explanation of overlapping peak fitting")
    button.setToolTip("Play a 16-second illustrated explanation in 1D and 2D. Uses synthetic peaks.")
    button.clicked.connect(lambda: open_fitting_help(owner))
    return button
