"""A small, local animation explaining joint multi-component peak fitting.

This is a scripted illustration with synthetic peaks, not an optimiser trace or
a chemical identification. The detector trace, cycle enlargement and folded map
use the same sampled profiles at every frame.
"""

from __future__ import annotations

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets
from .theme import fit_to_screen


class FittingHelpDialog(QtWidgets.QDialog):
    DURATION_MS = 16000
    CAPTIONS = (
        ("1. Narrow pulses repeat with each modulation",
         "The detector records narrow pulses: components are 250 ms wide here, with bursts 10 s apart. "
         "Stacking successive 10 s cycles gives the 2D map. The shaded cycle is enlarged below."),
        ("2. Propose separate contributions",
         "Ion-by-ion peak positions provide starting guesses. Blue and pink are proposed "
         "components; orange is their sum plus the local background. Follow the shaded cycle "
         "between the detector trace, enlarged pulse and 2D map."),
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
        fit_to_screen(self, 1040, 780)
        self.setMinimumSize(560, 460)
        self.position_ms = 0
        self._period = 10.
        self._dt = .005
        self._starts = np.arange(0., 80., self._period)
        self._x = self._starts + self._period/2
        self._y = np.arange(0., self._period, self._dt)
        self._detector_time = (self._starts[:, None] + self._y).ravel()
        self._xx, self._yy = np.meshgrid(self._x, self._y)
        # (first time, second time, first width, second width, correlation, height)
        # Conditional pulse FWHM = 2*sqrt(2*log(2))*sy*sqrt(1-rho**2).
        # Set both true components to 250 ms within each modulation.
        pulse_sigma = .250 / (2*np.sqrt(2*np.log(2)))
        self._truth = np.array([[35., 1.8, 14., pulse_sigma/np.sqrt(1-.55**2), -.55, 1.],
                                [45., 1.98, 12., pulse_sigma/np.sqrt(1-.45**2), .45, .9]])
        self._guess = np.array([[29., 1.65, 17., .16, 0., .7],
                                [50., 2.12, 15., .15, 0., .6]])
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

        self.trace = pg.PlotWidget(title="1D: detector trace, narrow pulses every 10 s")
        self.zoom = pg.PlotWidget(title="Selected cycle: resolve the overlapping pulses")
        self.map = pg.PlotWidget(title="2D: the same signal and fitted contours")
        self.trace.setLabel("bottom", "Elapsed detector time (s)")
        self.trace.setLabel("left", "Signal (a.u.)")
        self.zoom.setLabel("bottom", "Second retention time within cycle (s)")
        self.zoom.setLabel("left", "Signal (a.u.)")
        self.map.setLabel("bottom", "First retention time (illustrative seconds)")
        self.map.setLabel("left", "Second retention time (s)")
        self.trace.setXRange(0, 80, padding=0)
        self.trace.setYRange(0, self._observed.max() * 1.25, padding=0)
        self.zoom.setRange(xRange=(1.2, 2.6), yRange=(0, self._observed.max()*1.25), padding=0)
        self.map.setRange(xRange=(0, 80), yRange=(1.2, 2.6), padding=0)
        for plot in (self.trace, self.zoom, self.map):
            plot.setMenuEnabled(False)
            plot.hideButtons()
            plot.setMouseEnabled(x=False, y=False)
        self.image = pg.ImageItem(axisOrder="row-major")
        self.image.setImage(self._observed, autoLevels=False, levels=(0, self._observed.max()))
        self.image.setRect(QtCore.QRectF(0, -self._dt/2, 80, self._period))
        self.map.addItem(self.image)
        self.observed_curve = self.trace.plot(self._detector_time, self._observed.T.ravel())
        self.sum_curve = self.trace.plot()
        self.component_curves = [self.trace.plot(), self.trace.plot()]
        self.zoom_observed = self.zoom.plot()
        self.zoom_sum = self.zoom.plot()
        self.zoom_components = [self.zoom.plot(), self.zoom.plot()]
        self.cycle_bands = []
        for plot in (self.trace, self.map):
            band = pg.LinearRegionItem(values=(30, 40), movable=False)
            band.setZValue(5)
            plot.addItem(band)
            self.cycle_bands.append(band)
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
        plots.addWidget(self.zoom, 1)
        plots.addWidget(self.map, 1)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(title)
        layout.addWidget(self.heading)
        layout.addWidget(self.caption)
        layout.addWidget(self.trace, 1)
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
        for plot in (self.trace, self.zoom, self.map):
            plot.setBackground(base)
            for axis in ("bottom", "left"):
                plot.getAxis(axis).setPen(foreground)
                plot.getAxis(axis).setTextPen(foreground)
        self.trace.setTitle("1D: detector trace, narrow pulses every 10 s", color=foreground.name())
        self.zoom.setTitle(self.zoom.plotItem.titleLabel.text, color=foreground.name())
        self.map.setTitle("2D: the same signal and fitted contours", color=foreground.name())
        self.image.setColorMap(pg.ColorMap([0, 1], [base, foreground]))
        self.observed_curve.setPen(pg.mkPen(foreground, width=2))
        self.zoom_observed.setPen(pg.mkPen(foreground, width=2))
        self.sum_curve.setPen(pg.mkPen("#e69824" if dark else "#a85b00", width=2,
                                      style=QtCore.Qt.PenStyle.DashLine))
        self.zoom_sum.setPen(self.sum_curve.opts["pen"])
        for curve, colour in zip(self.zoom_components, colours):
            curve.setPen(pg.mkPen(colour, width=2.5))
        for band in self.cycle_bands:
            band.setBrush(pg.mkBrush(100, 140, 200, 35))
            band.setHoverBrush(pg.mkBrush(100, 140, 200, 35))
            for line in band.lines:
                line.setPen(pg.mkPen(foreground, width=1, style=QtCore.Qt.PenStyle.DotLine))
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
        """Unfold the same sampled modulation columns into detector time."""
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
        total = .025 + sum(profiles)
        self.sum_curve.setData(self._detector_time, total.T.ravel())
        self.sum_curve.setVisible(stage > 0)
        # Visit the central four cycles, then hold the strongest overlap for review.
        self._selected_cycle = (2 + min(self.position_ms // 3000, 3)
                                if self.position_ms < 12000 else 3)
        selected = self._selected_cycle
        start = self._starts[selected]
        for band in self.cycle_bands:
            band.setRegion((start, start + self._period))
        foreground = QtWidgets.QApplication.palette().color(QtGui.QPalette.ColorRole.WindowText)
        self.zoom.setTitle(f"Cycle {start:.0f}\u2013{start+self._period:.0f} s: enlarged pulses",
                           color=foreground.name())
        self.zoom_observed.setData(self._y, self._observed[:, selected])
        self.zoom_sum.setData(self._y, total[:, selected])
        self.zoom_sum.setVisible(stage > 0)
        for curve, profile in zip(self.zoom_components, profiles):
            curve.setData(self._y, profile[:, selected])
            curve.setVisible(stage > 0)
        theta = np.linspace(0, 2*np.pi, 101)
        for p, profile, curve, contours in zip(parameters, profiles, self.component_curves, self.contours):
            curve.setData(self._detector_time, profile.T.ravel())
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
