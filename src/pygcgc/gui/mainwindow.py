"""PyGCGC main window."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets

from .. import __version__
from ..core import (IonExpressionError, chromatogram_caption, find_homologous_series,
                    parse_ion_expr, plot_chromatogram, qualifier_ratios_at, reference_from_spectrum,
                    save_figure, siloxane_scores, spectrum_at)
from ..core.peaks import make_frame
from ..core.io import file_filter
from .blob_state import LEGACY, MS_AWARE
from .batch_dialog import BatchDialog
from .commands import (ApplyBackground, ApplySiloxane, EditAlignment, EditBlobs, EditComponents, EditDecisions,
                       EditHydrocarbons, SetFold, SetPeriod, SetPhase, SetPhases, SetReference)
from .dialogs import QuantSettingsDialog, SettingsDialog
from .identify_dialog import IdentifyDialog
from .jobs import JobManager, JobsPanel
from .overlap_dialog import OverlapDialog
from .panels import FoldPanel, LayersPanel, RunsPanel
from .panels.align import AlignPanel
from .panels.background import BackgroundPanel
from .panels.blobs import BlobsPanel
from .panels.calibration import CalibrationPanel
from .panels.hydrocarbons import HydrocarbonsPanel
from .panels.results import ResultsPanel
from .panels.siloxane import SiloxanePanel
from .project import Project
from .render import Renderer
from .theme import apply_theme
from .workflow import IdentifyStepPanel, QuantStepPanel, WorkflowPanel
from .views.compare import CompareView
from .views.image2d import SCALE_NAMES, SCALES, ChromatogramView
from .views.review import ReviewView
from .views.spectrum import SpectrumView
from .views.trace1d import OneDView

SHORTCUTS = [
    ("PgUp / PgDn", "Previous / next run (view stays put)"),
    ("M", "Type an ion expression: 178, 178+252, 200:260, 57/71, TIC"),
    ("L", "Cycle colour scale: log, square root, linear"),
    ("E", "Extend the t₂ axis to show wrapped peaks whole"),
    ("R", "Region box: summed spectrum of everything inside"),
    ("D", "Toggle between the 2D image and the raw 1D chromatogram"),
    ("I", "In review: redraw the selected peak's integration footprint (then Ctrl+Enter)"),
    ("Space", "Show raw or 🧹 cleaned data"),
    ("O", "Show or hide 🧪 peak outlines and labels"),
    ("3 / 4 / 5", "🧪 Review gallery / calibration curves / 📍 Compare"),
    ("A", "📍 Add anchor: click the reference, then the same peak here"),
    ("B", "Jump to Compare in blink mode"),
    ("G", "🔍 Merge the selected blobs"),
    ("Y / N / U", "In review or results: accept / reject / clear"),
    ("P", "In review: pick a different peak by clicking the 2D view"),
    ("Esc", "Cancel picking a peak"),
    ("Home", "Reset zoom"),
    ("1 / 2", "2D view / 1D view"),
    ("Ctrl+O", "Open a folder of runs"),
    ("Ctrl+S", "Save the project"),
    ("Ctrl+Z / Ctrl+Shift+Z", "Undo / redo"),
    ("F1", "This list"),
]


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, dark: bool = False):
        super().__init__()
        self.setWindowTitle(f"PyGCGC {__version__}")
        self.resize(1500, 920)
        self.dark = dark
        self.jobs = JobManager(self)
        self.project = Project(self.jobs, self)
        self.undo = QtGui.QUndoStack(self)
        self.renderer = Renderer(self.project, self.jobs, parent=self)
        self.index = -1
        self._pending_phase_undo = False
        self._last_click: float | None = None
        self._picking: tuple[str, str] | None = None
        self._frame_cache: dict = {}
        self.quant = self.project.quant
        self.align = self.project.align
        self.renderer.set_align(self.align)
        self._anchor_stage: str | None = None      # None, "ref", "run"
        self._anchor_pending: dict = {}
        self.blobs = self.project.blobs
        self._blob_add_pending = False
        self.hydrocarbons = self.project.hydrocarbons
        self._alkane_seed_pending: list | None = None
        self._siloxane_pick_pending = False
        self._qual_pick_pending = False
        self._reintegrate: tuple[str, str] | None = None     # (run, compound) whose footprint is being redrawn

        # ---- centre: 2D and 1D tabs
        self.view2d = ChromatogramView()
        self.view1d = OneDView()
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self.view2d, "2D")
        self.tabs.addTab(self.view1d, "1D")
        self.review = ReviewView(self._review_image)
        self.calibration_panel = CalibrationPanel()
        self.tabs.addTab(self.review, "🧪 Review")
        self.tabs.addTab(self.calibration_panel, "Calibration")
        self.compare = CompareView(self.align)
        self.compare.image_source = self._compare_image
        self.tabs.addTab(self.compare, "Compare")
        i3d = self.tabs.addTab(QtWidgets.QLabel("Arrives in milestone 8"), "3D")
        self.tabs.setTabEnabled(i3d, False)
        self.setCentralWidget(self.tabs)

        # ---- docks
        self.act_processed = QtGui.QAction("🧹 Show cleaned", self)
        self.act_processed.setCheckable(True)
        self.act_processed.setChecked(True)
        self.act_processed.setEnabled(False)
        self.act_processed.setShortcut(QtGui.QKeySequence("Space"))
        self.act_processed.setToolTip("Show raw or cleaned data (Space)")
        self.act_processed.toggled.connect(self.project.set_show_processed)
        self.addAction(self.act_processed)
        self.act_aligned = QtGui.QAction("📍 Show aligned", self)
        self.act_aligned.setCheckable(True)
        self.act_aligned.setEnabled(False)
        self.act_aligned.setToolTip("Show this run resampled onto the reference run's grid")
        self.act_aligned.toggled.connect(self.renderer.set_aligned)
        self.addAction(self.act_aligned)
        self.runs_panel = RunsPanel(self.project)
        self.fold_panel = FoldPanel(self.project)
        self.layers_panel = LayersPanel(self.project, self.act_processed, self.act_aligned)
        self.bg_panel = BackgroundPanel()
        self.sil_panel = SiloxanePanel()
        self.sil_panel.set_settings(self.project.siloxane_settings)
        self.align_panel = AlignPanel()
        self.align_panel.set_align(self.align)
        self.blobs_panel = BlobsPanel()
        self.hydro_panel = HydrocarbonsPanel()
        self.spectrum = SpectrumView()
        self.jobs_panel = JobsPanel(self.jobs)
        self.results = ResultsPanel()
        self.history = QtWidgets.QUndoView(self.undo)
        self.history.setEmptyLabel("Opened")
        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)

        # Left: the run list on top, then one workflow panel. Each step's settings
        # replace the last, with a single Go button underneath (see gui/workflow.py).
        # The old panels keep their own primary buttons but hide them; Go does that job.
        for w in (self.bg_panel.apply_all_btn, self.sil_panel.apply_btn, self.blobs_panel.find_btn,
                  self.hydro_panel.find_btn, self.align_panel.propagate_btn, self.fold_panel.fix_btn):
            w.setVisible(False)
        self.quant_step = QuantStepPanel()
        self.identify_step = IdentifyStepPanel()
        self.workflow = WorkflowPanel()
        add = self.workflow.add_step
        add("fold", "🔄 Fold", "Fold and phase",
            "Cuts the detector stream into modulations and stacks them into the 2D image. Set the "
            "modulation period, then estimate the phase (where each modulation starts) of every run.",
            "Every later step works on the folded image. A wrong phase splits each peak across the "
            "top and bottom edges of the image.",
            self.fold_panel, "Fix phase on all runs")
        add("align", "📍 Align", "Alignment",
            "Retention times drift by a few seconds between runs. Pick a reference run, mark the same "
            "peaks (anchors) on it and on one other run, then find them on every run. Skip this with "
            "a single run or no blank.",
            "Aligning now puts every run on the reference run's coordinates, so the background step "
            "can subtract a blank through the alignment and later steps compare like with like.",
            self.align_panel, "Auto-find anchors on all runs", scrolled=False)
        add("background", "🧹 Background", "Background removal",
            "Removes column bleed, the stripe pattern that repeats every modulation, and single-scan "
            "spikes, and can subtract a blank. Preview on this run before applying to every run.",
            "Peaks are found by signal-to-noise, and the noise map comes from this step. It follows "
            "alignment so that a blank can be subtracted through it.",
            self.bg_panel, "Remove background")
        add("blobs", "🔍 Peaks", "Peak detection and separation",
            "Choose MS-aware separation or the older hard-border blob detector, then find peaks on "
            "every run. MS-aware regions can contain several overlapping compounds.",
            "Detection uses the cleaned baseline from the background step. Existing siloxane, ladder "
            "and library-search tools still use legacy blobs while their component-aware paths are built.",
            self.blobs_panel, "Find peaks", scrolled=False)
        add("siloxane", "🫧 Siloxane", "Siloxane removal",
            "Scores each blob against a siloxane reference spectrum and blocks out the ones above the "
            "certainty you set. Choose the reference and the threshold, then press Go.",
            "It works on whole blobs, so it needs step 4 first. It runs before the ladder and "
            "identification so contaminant peaks do not enter them.",
            self.sil_panel, "Remove siloxanes")
        add("hydro", "⛓️ Hydrocarbons", "n-Alkane ladder and stencils",
            "Anchors the n-alkane ladder to give retention indices, and draws group-type stencils "
            "(lumped classes) that transfer between runs. Go starts a two-click flow: click two "
            "n-alkanes on the 2D view and type their carbon numbers.",
            "The ladder is searched for on clean data and gives the retention index that identification "
            "(step 7) uses.",
            self.hydro_panel, "Find alkane series")
        add("identify", "🏷️ Identify", "Identify blobs",
            "Matches each blob's spectrum to a mass spectral library (MSP) and reports a confidence "
            "tier. Go opens the identifier window.",
            "It searches the blobs from step 4 and uses the retention index from step 6 when the "
            "ladder exists, so identification comes after both.",
            self.identify_step, "Open identifier")
        add("quant", "🧪 Quantify", "Quantification",
            "Finds every target and internal standard in every run, fits the calibration curves and "
            "quantifies the samples.",
            "Last: it needs the phase (step 1) and works on the cleaned data (step 3). Everything before "
            "it is there to make its peaks clean and correctly placed.",
            self.quant_step, "Quantify")
        # One left-hand dock with three tabs, so the run list, the layers and the open
        # step's settings share the space instead of competing for it.
        self.left_tabs = QtWidgets.QTabWidget()
        self.left_tabs.addTab(self.runs_panel, "Runs")
        self.tool_tab = self.left_tabs.addTab(self.workflow, "Tool")
        self.left_tabs.addTab(self.layers_panel, "Layers")
        left = self._dock("Workspace", self.left_tabs, QtCore.Qt.DockWidgetArea.LeftDockWidgetArea)
        spec_dock = self._dock("Spectrum", self.spectrum, QtCore.Qt.DockWidgetArea.RightDockWidgetArea)
        self.results_dock = self._dock("🧪 Results", self.results, QtCore.Qt.DockWidgetArea.BottomDockWidgetArea)
        jobs_dock = self._dock("Jobs", self.jobs_panel, QtCore.Qt.DockWidgetArea.BottomDockWidgetArea)
        hist_dock = self._dock("History", self.history, QtCore.Qt.DockWidgetArea.BottomDockWidgetArea)
        log_dock = self._dock("Log", self.log, QtCore.Qt.DockWidgetArea.BottomDockWidgetArea)
        self.tabifyDockWidget(self.results_dock, jobs_dock)
        self.tabifyDockWidget(jobs_dock, hist_dock)
        self.tabifyDockWidget(hist_dock, log_dock)
        jobs_dock.raise_()
        self.resizeDocks([left, spec_dock], [470, 330], QtCore.Qt.Orientation.Horizontal)
        self.resizeDocks([jobs_dock], [190], QtCore.Qt.Orientation.Vertical)
        self.workspace_dock, self.spectrum_dock = left, spec_dock

        self._build_actions()
        self._build_toolbars()
        self._build_menus()
        self._build_status()
        self._connect()
        apply_theme(self, dark)
        self.left_tabs.setTabText(self.tool_tab, "Fold")
        self._default_layout = self.saveState()

    # ================================================================ setup
    def _dock(self, title, widget, area):
        dock = QtWidgets.QDockWidget(title, self)
        dock.setObjectName(f"dock_{title}")
        dock.setWidget(widget)
        self.addDockWidget(area, dock)
        return dock

    @staticmethod
    def _scrolled(widget):
        area = QtWidgets.QScrollArea()
        area.setWidgetResizable(True)
        area.setWidget(widget)
        return area

    def _action(self, text, slot=None, shortcut=None, tip=None, checkable=False):
        a = QtGui.QAction(text, self)
        if shortcut:
            a.setShortcut(QtGui.QKeySequence(shortcut))
            a.setShortcutContext(QtCore.Qt.ShortcutContext.WindowShortcut)
        if tip:
            a.setToolTip(tip)
            a.setStatusTip(tip)
        a.setCheckable(checkable)
        if slot is not None:
            (a.toggled if checkable else a.triggered).connect(slot)
        self.addAction(a)
        return a

    def _build_actions(self):
        A = self._action
        self.act_open = A("📂 Open", self.open_folder_dialog, "Ctrl+O", "Open runs (native LSC defaults to untreated 70 eV)")
        self.act_open_files = A("Open files…", self.open_files_dialog, None, "Choose individual runs, including 14 eV or vendor DBC")
        self.act_open_project = A("Open project…", self.open_project_dialog, "Ctrl+Shift+O")
        self.act_save = A("Save project", self.save_project, "Ctrl+S")
        self.act_save_as = A("Save project as…", self.save_project_as)
        self.act_targets = A("Load target table…", self.load_targets_dialog, None,
                             "Internal standards in the table become phase anchors")
        self.act_sequence = A("Load sequence list…", self.load_sequence_dialog)
        self.act_fix = A("🔄 Fix phase", self.fix_phase, None, "Estimate the phase of every run")
        self.act_background = A("🧹 Background", lambda: self.apply_background(), None,
                                "Remove background and spikes from every run, using the 🧹 panel settings")
        # Toolbar step buttons only open a step's settings on the left; Go (or the
        # act_* actions above, from the menus) does the work.
        self.step_actions = {}
        self._step_names = {}
        step_group = QtGui.QActionGroup(self)
        for key, text, tip in ([
                ("fold", "1 🔄 Fold", "Modulation period and phase"),
                ("align", "2 📍 Align", "Reference run and anchors"),
                ("background", "3 🧹 Background", "Background removal settings"),
                ("blobs", "4 🔍 Peaks", "Peak detection and separation"),
                ("siloxane", "5 🫧 Siloxane", "Siloxane removal settings"),
                ("hydro", "6 ⛓️ Hydrocarbons", "n-Alkane ladder and stencils"),
                ("identify", "7 🏷️ Identify", "Library search"),
                ("quant", "8 🧪 Quantify", "Targets, calibration and quantification")]):
            info = self.workflow.step_info(key)
            act = A(text, lambda _checked=False, k=key: self.workflow.show_step(k), None,
                    f"{tip}. Opens its settings on the left; press Go there to run it.\n"
                    f"Why here: {info['why']}")
            act.setCheckable(True)
            step_group.addAction(act)
            self.step_actions[key] = act
            self._step_names[key] = text
        self.step_actions["fold"].setChecked(True)
        self.act_siloxane = A("🫧 Siloxanes", lambda: self.apply_siloxane(), None,
                              "Fit and subtract the siloxane component on every loaded run, "
                              "on top of 🧹 where that's been applied")
        self.act_align_panel = A("📍 Align", self.show_align_panel, None,
                                 "Reference run, anchors, and batch alignment")
        self.act_add_anchor = A("➕ Add anchor", self.start_add_anchor, "A",
                                "Click the peak on the reference, then the same peak here")
        self.act_compare = A("Compare", lambda: self.tabs.setCurrentWidget(self.compare), "5")
        self.act_blink = A("Blink compare", self._go_blink, "B", "Jump to Compare in blink mode")
        self.act_quant = A("🧪 Quantify", self.run_quant, None,
                           "Find every target and internal standard in every run, fit the calibrations "
                           "and quantify the samples")
        self.act_quant_settings = A("Quantification settings…", self.quant_settings_dialog)
        self.act_levels = A("Load calibration levels…", self.load_levels_dialog, None,
                            "CSV with a level column and one column per compound")
        self.act_export = A("Export 🧪 results…", self.export_results)
        self.act_overlay = A("Peaks", self.view2d.set_overlay_visible, "O",
                             "Show 🧪 peak outlines and labels", checkable=True)
        self.act_overlay.setChecked(True)
        self.act_blobs_panel = A("🔍 Peaks", self.show_blobs_panel, None,
                                 "Choose and inspect either MS-aware peak regions or legacy blobs")
        self.act_find_blobs = A("🔍 Find peaks", self.run_blobs, None,
                                "Run the selected peak-detection method on every loaded run")
        self.act_blob_overlay = A("🔍 Peak detections", self.view2d.set_blob_overlay_visible, None,
                                  "Show detected centres and regions or legacy outlines", checkable=True)
        self.act_add_blob = A("➕ Add blob", self.start_add_blob, None,
                              "Click a peak in the 2D view that automatic detection missed")
        self.act_merge_blobs = A("Merge selected blobs", self.merge_selected_blobs, "G",
                                 "Combine the 🔍 blobs selected in the Blobs panel into one")
        self.act_overlap = A("✂️ Resolve overlap", lambda: self.show_overlap_dialog(), None,
                            "Estimate overlapping component contributions inside an MS-aware region "
                            "or inspect a legacy blob")
        self.act_hydro_panel = A("⛓️ Hydrocarbons", self.show_hydrocarbons_panel, None,
                                 "The n-alkane ladder and group-type stencils")
        self.act_find_alkanes = A("➕ Find alkane series", self.start_find_alkanes, None,
                                  "Click two n-alkanes and type their carbon numbers")
        self.act_stencil_overlay = A("⛓️ Stencil outlines", self.view2d.set_stencil_overlay_visible, None,
                                     "Show every stencil, lifted onto this run's own alkane series",
                                     checkable=True)
        self.act_stencil_overlay.setChecked(True)
        self.act_identify = A("🏷️ Identify", self.show_identify_dialog, None,
                             "Search 🔍 blobs against an MSP library: retention index, weighted dot "
                             "product, target-decoy confidence")
        self.act_batch = A("🗂️ Batch", self.show_batch_dialog, None,
                          "Process a whole folder: cached, parallel pipeline, feature table, QC, PCA and volcano")
        self.act_blob_settings = A("🔍 Peak detection settings…", self.edit_blob_settings, None,
                                   "Settings for the detection method currently selected in step 4")
        self.act_hydro_settings = A("⛓️ Hydrocarbon settings…", self.edit_hydro_settings, None,
                                    "Every alkane-search, qualifier and series setting")
        self.act_export_figure = A("📤 Figure…", self.export_figure, None,
                                   "Save what the 2D view shows as a publication figure (89 or 183 mm wide, "
                                   "300 dpi) with a caption file")
        self.act_reintegrate = A("🧪 Re-integrate peak…", self.reintegrate_current, None,
                                 "Redraw the integration footprint of the peak selected in the Review gallery "
                                 "(key I in the gallery): drag its corners, then Ctrl+Enter")
        self.act_accept_footprint = A("✓ Use this footprint", self.accept_footprint, "Ctrl+Return",
                                      "Integrate the peak inside the orange polygon")
        self.act_accept_footprint.setEnabled(False)
        self.act_reset_footprint = A("↺ Automatic footprint", self.reset_footprint, None,
                                     "Remove the hand-drawn footprint of the peak being edited")
        self.act_reset_footprint.setEnabled(False)
        self.act_review = A("🧪 Review", lambda: self.tabs.setCurrentWidget(self.review), "3")
        self.act_calib = A("Calibration curves", lambda: self.tabs.setCurrentWidget(self.calibration_panel), "4")
        self.act_cancel_pick = A("Cancel pick", self.cancel_pick, "Esc")
        self.act_snr = A("S/N", self.renderer.set_snr, None,
                         "Show signal-to-noise (needs 🧹 on this run)", checkable=True)
        self.act_snr.setEnabled(False)
        self.act_undo = self.undo.createUndoAction(self, "Undo")
        self.act_undo.setShortcut(QtGui.QKeySequence.StandardKey.Undo)
        self.act_redo = self.undo.createRedoAction(self, "Redo")
        self.act_redo.setShortcut(QtGui.QKeySequence("Ctrl+Shift+Z"))
        self.addAction(self.act_undo)
        self.addAction(self.act_redo)
        self.act_prev = A("Previous run", lambda: self.step_run(-1), "PgUp")
        self.act_next = A("Next run", lambda: self.step_run(+1), "PgDown")
        self.act_home = A("Reset zoom", self.view2d.reset_view, "Home")
        self.act_scale = A("Cycle colour scale", self.cycle_scale, "L")
        self.act_ion = A("Choose ion", self.focus_ion_box, "M")
        self.act_extend = A("Extend t₂", self.renderer.set_extended, "E",
                            "Show wrapped peaks whole above the period line", checkable=True)
        self.act_region = A("▭ Region", self.view2d.set_region_tool, "R",
                            "Draggable box: summed spectrum of everything inside", checkable=True)
        self.act_lock = A("🔒 Lock levels", self.set_levels_locked, None,
                          "Keep colour levels when flicking between runs", checkable=True)
        self.act_lock.setChecked(True)
        self.act_auto = A("Auto levels", self.view2d.auto_levels, "Ctrl+L")
        self.act_view2d = A("2D view", lambda: self.tabs.setCurrentIndex(0), "1")
        self.act_view1d = A("1D view", lambda: self.tabs.setCurrentIndex(1), "2")
        self.act_toggle_dim = A("Toggle 2D / 1D", lambda: self.tabs.setCurrentIndex(1 if self.tabs.currentIndex() == 0 else 0),
                                "D", "Switch between the 2D image and the raw 1D chromatogram (reconstructed trace "
                                      "plus the unfolded detector stream)")
        self.act_help = A("Keyboard shortcuts", self.show_shortcuts, "F1")
        self.act_theme = A("Dark theme", self.set_dark, None, checkable=True)
        self.act_theme.setChecked(self.dark)
        self.act_quit = A("Quit", self.close, "Ctrl+Q")

    def _build_toolbars(self):
        tb = self.addToolBar("Main")
        tb.setObjectName("toolbar_main")
        tb.setMovable(False)
        tb.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextOnly)
        # Left to right, the order one works through the data: load, then the eight
        # numbered steps (each opens its settings on the left), then the extras.
        tb.addAction(self.act_open)
        tb.addSeparator()
        for key in ("fold", "align", "background", "blobs", "siloxane", "hydro", "identify", "quant"):
            tb.addAction(self.step_actions[key])
        tb.addSeparator()
        tb.addAction(self.act_overlap)
        tb.addAction(self.act_batch)
        tb.addAction(self.act_export_figure)
        # Closed panels can always be brought back from here (also View → Panels),
        # kept on this bar because the View bar is full.
        tb.addSeparator()
        for dock, tip in ((self.workspace_dock, "Show or hide the Runs / Tool / Layers panel"),
                          (self.spectrum_dock, "Show or hide the spectrum panel")):
            act = dock.toggleViewAction()
            act.setToolTip(tip)
            tb.addAction(act)
        panels = QtWidgets.QToolButton()
        panels.setText("Panels ▾")
        panels.setToolTip("Every panel: tick to show, untick to hide. Reset layout restores the default.")
        panels.setPopupMode(QtWidgets.QToolButton.ToolButtonPopupMode.InstantPopup)
        self.panels_menu = QtWidgets.QMenu(panels)
        panels.setMenu(self.panels_menu)
        tb.addWidget(panels)

        self.addToolBarBreak()
        vb = self.addToolBar("View")
        vb.setObjectName("toolbar_view")
        vb.setMovable(False)
        vb.addWidget(QtWidgets.QLabel(" Ion "))
        self.ion_box = QtWidgets.QComboBox()
        self.ion_box.setEditable(True)
        self.ion_box.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.InsertAtTop)
        self.ion_box.setMinimumWidth(170)
        self.ion_box.addItems(["TIC", "57", "178", "202", "228", "252", "57/71", "178+202+228+252"])
        self.ion_box.lineEdit().setPlaceholderText("TIC, 178, 178+252, 57/71")
        self.ion_box.setToolTip("Ion expression (M). Examples: TIC, 178, 178+202, 200:260, 57/71")
        self.ion_box.lineEdit().returnPressed.connect(self.apply_ion_box)
        self.ion_box.activated.connect(lambda _i: self.apply_ion_box())
        vb.addWidget(self.ion_box)
        vb.addWidget(QtWidgets.QLabel("  Scale "))
        self.scale_box = QtWidgets.QComboBox()
        for s in SCALES:
            self.scale_box.addItem(SCALE_NAMES[s], s)
        self.scale_box.currentIndexChanged.connect(lambda i: self.set_scale(SCALES[i]))
        vb.addWidget(self.scale_box)
        vb.addAction(self.act_snr)
        vb.addAction(self.act_overlay)
        vb.addAction(self.act_blob_overlay)
        vb.addAction(self.act_stencil_overlay)
        vb.addSeparator()
        vb.addAction(self.act_processed)
        vb.addAction(self.act_lock)
        vb.addAction(self.act_auto)
        vb.addAction(self.act_extend)
        vb.addAction(self.act_region)
        vb.addAction(self.act_home)
        vb.addSeparator()
        vb.addAction(self.act_prev)
        vb.addAction(self.act_next)

    def _build_menus(self):
        m = self.menuBar()
        f = m.addMenu("&File")
        for a in (self.act_open, self.act_open_files, None, self.act_open_project, self.act_save,
                  self.act_save_as, None, self.act_sequence, self.act_targets, self.act_levels,
                  None, self.act_export, self.act_export_figure, None, self.act_quit):
            f.addSeparator() if a is None else f.addAction(a)
        e = m.addMenu("&Edit")
        e.addAction(self.act_undo)
        e.addAction(self.act_redo)
        v = m.addMenu("&View")
        for a in (self.act_view2d, self.act_view1d, self.act_toggle_dim, self.act_compare, None,
                  self.act_processed, self.act_aligned, self.act_snr,
                  self.act_blob_overlay, self.act_stencil_overlay, self.act_ion, self.act_scale, self.act_lock,
                  self.act_auto, self.act_extend, self.act_region, self.act_home, None,
                  self.act_prev, self.act_next, None, self.act_theme):
            v.addSeparator() if a is None else v.addAction(a)
        docks = v.addMenu("Panels")
        for menu in (docks, self.panels_menu):
            for d in self.findChildren(QtWidgets.QDockWidget):
                menu.addAction(d.toggleViewAction())
            menu.addSeparator()
            menu.addAction("Reset layout", self.reset_layout)
        # The Process menu follows the workflow order, one headed section per step.
        p = m.addMenu("&Process")
        p.addSection("1  Fold")
        p.addAction(self.act_fix)
        p.addSection("2  Align")
        p.addAction(self.act_align_panel)
        p.addAction(self.act_add_anchor)
        p.addAction(QtGui.QAction("📍 Auto-find anchors on all runs", self, triggered=self.propagate_anchors))
        p.addSection("3  Background")
        p.addAction(self.act_background)
        p.addAction(QtGui.QAction("Remove 🧹 from all runs", self, triggered=self.remove_background_all))
        p.addAction(QtGui.QAction("Background settings…", self, triggered=self.show_background_panel))
        p.addSection("4  Peak detection")
        p.addAction(self.act_blobs_panel)
        p.addAction(self.act_find_blobs)
        p.addAction(self.act_blob_settings)
        p.addAction(self.act_add_blob)
        p.addAction(self.act_merge_blobs)
        p.addAction(self.act_overlap)
        p.addSection("5  Siloxane")
        p.addAction(self.act_siloxane)
        p.addAction(QtGui.QAction("Remove 🫧 from all runs", self, triggered=self.remove_siloxane_all))
        p.addSection("6  Hydrocarbons")
        p.addAction(self.act_hydro_panel)
        p.addAction(self.act_find_alkanes)
        p.addAction(QtGui.QAction("⛓️ Propagate alkanes to all runs", self, triggered=self.propagate_alkanes))
        p.addAction(self.act_hydro_settings)
        p.addSection("7  Identify")
        p.addAction(self.act_identify)
        p.addSection("8  Quantify")
        p.addAction(self.act_quant)
        p.addAction(self.act_quant_settings)
        p.addAction(self.act_review)
        p.addAction(self.act_reintegrate)
        p.addAction(self.act_accept_footprint)
        p.addAction(self.act_reset_footprint)
        p.addAction(self.act_calib)
        p.addAction(self.act_export)
        p.addSection("Whole folder")
        p.addAction(self.act_batch)
        h = m.addMenu("&Help")
        h.addAction("Workflow guide: what each step does and why", self.show_workflow_guide)
        h.addAction(self.act_help)
        h.addAction(QtGui.QAction(f"About PyGCGC {__version__}", self, triggered=self.show_about))

    def _build_status(self):
        sb = self.statusBar()
        self.cursor_label = QtWidgets.QLabel("")
        self.run_label = QtWidgets.QLabel("No runs")
        self.jobs_label = QtWidgets.QLabel("")
        sb.addWidget(self.cursor_label, 1)
        sb.addPermanentWidget(self.jobs_label)
        sb.addPermanentWidget(self.run_label)

    def go(self, key: str):
        """The workflow panel's Go button: run the step that is open."""
        {"fold": self.fix_phase,
         "background": self.apply_background,
         "siloxane": self.apply_siloxane,
         "align": self.propagate_anchors,
         "blobs": self.run_blobs,
         "hydro": self.start_find_alkanes,
         "identify": self._go_identify,
         "quant": self.run_quant}[key]()

    def reset_layout(self):
        self.restoreState(self._default_layout)
        for d in self.findChildren(QtWidgets.QDockWidget):
            d.show()

    def show_workflow_guide(self):
        box = QtWidgets.QMessageBox(self)
        box.setWindowTitle("Workflow guide")
        box.setTextFormat(QtCore.Qt.TextFormat.RichText)
        box.setText(self.workflow.guide_html())
        box.exec()

    def _step_opened(self, key: str):
        """Opening a step brings the Tool tab forward and names it after the step."""
        self.step_actions[key].setChecked(True)
        self.left_tabs.setTabText(self.tool_tab, self._step_names[key].split(" ", 1)[1])
        self.left_tabs.setCurrentIndex(self.tool_tab)
        self.workspace_dock.show()
        self.workspace_dock.raise_()

    def _step_progress(self):
        """Toolbar step buttons carry the progress: ✓ when done, ▸ on the suggested next step."""
        nxt = self.workflow.next_key()
        for key, act in self.step_actions.items():
            mark = "✓ " if self.workflow.is_done(key) else ("▸ " if key == nxt else "")
            act.setText(mark + self._step_names[key])

    def _go_identify(self):
        self.show_identify_dialog()
        self.workflow.mark_done("identify", "Identifier opened. Load a library there and search.")

    def _quant_step_status(self, *_):
        ok, why = self.quant.ready()
        self.quant_step.set_status(ok, why, bool(self.quant.levels), self.project.any_background())
        self.workflow.set_go_enabled("quant", ok and not self.quant.running, why if not ok else "")

    def _connect(self):
        self.project.runsChanged.connect(self._runs_changed)
        self.project.message.connect(self.say)
        self.renderer.ready.connect(self._rendered)
        self.renderer.failed.connect(lambda t: self.say(f"Drawing failed: {t}"))
        self.runs_panel.currentChanged.connect(self.set_run)
        self.runs_panel.removeRequested.connect(self.remove_runs)
        self.fold_panel.phaseEdited.connect(self.edit_phase)
        self.fold_panel.periodEdited.connect(self.edit_period)
        self.fold_panel.periodApplyRequested.connect(self.edit_period)
        self.fold_panel.fixPhaseRequested.connect(self.fix_phase)
        self.fold_panel.estimatePeriodRequested.connect(self.estimate_period)
        self.fold_panel.copyPhaseRequested.connect(self.copy_phase_to_all)
        self.view2d.hovered.connect(self._hovered)
        self.view2d.pixelClicked.connect(self.show_spectrum_at)
        self.view2d.regionChanged.connect(self.show_region_spectrum)
        self.view2d.bandChanged.connect(self._band_from_2d)
        self.view1d.bandChanged.connect(self._band_from_1d)
        self.view1d.timeClicked.connect(self._time_from_1d)
        self.spectrum.ionClicked.connect(self.ion_clicked)
        self.jobs.busyChanged.connect(self._busy)
        self.project.layersChanged.connect(self._layers_changed)
        self.project.runsChanged.connect(lambda: self.bg_panel.set_run_names(self.project.entries))
        self.bg_panel.applyAllRequested.connect(lambda: self.apply_background())
        self.bg_panel.applyThisRequested.connect(
            lambda: self.index >= 0 and self.apply_background([self.index]))
        self.bg_panel.removeAllRequested.connect(self.remove_background_all)
        self.bg_panel.settingsChanged.connect(self._bg_settings_changed)
        self.bg_panel.previewToggled.connect(self._preview_toggled)
        self.bg_panel.splitToggled.connect(self.renderer.set_split)
        self.bg_panel.diagnosticIonChanged.connect(lambda _m: self.update_background_info())
        self.jobs.jobChanged.connect(self._job_logged)
        # 🫧
        self.sil_panel.settingsChanged.connect(self._siloxane_settings_changed)
        self.sil_panel.pickToggled.connect(self._siloxane_pick_toggled)
        self.sil_panel.applyRequested.connect(lambda: self.apply_siloxane())
        self.project.runsChanged.connect(self._sync_siloxane_panel)
        # 🧪
        self.quant.changed.connect(self._quant_changed)
        self.quant.runningChanged.connect(self._quant_running)
        self.view2d.pointClicked.connect(self._point_clicked)
        self.results.gotoRequested.connect(self._goto_row)
        self.results.decisionRequested.connect(self.decide)
        self.results.exportRequested.connect(self.export_results)
        self.results.settingsRequested.connect(self.quant_settings_dialog)
        self.results.runRequested.connect(self.run_quant)
        self.review.gotoRequested.connect(self.goto_peak)
        self.review.decisionRequested.connect(self.decide)
        self.review.pickRequested.connect(self.start_pick)
        self.review.reintegrateRequested.connect(self.start_reintegrate)
        self.calibration_panel.gotoRequested.connect(self.goto_peak)
        self.calibration_panel.overrideRequested.connect(self.override_curve)
        self.calibration_panel.excludeRequested.connect(self.exclude_point)
        # 📍
        self.align.changed.connect(self._align_changed)
        self.project.runsChanged.connect(lambda: self.align_panel.refresh_runs(self.project.entries))
        self.align_panel.referenceChanged.connect(self.set_reference)
        self.align_panel.addAnchorRequested.connect(self.start_add_anchor)
        self.align_panel.removeAnchorRequested.connect(self.remove_anchor)
        self.align_panel.clearRunRequested.connect(self.clear_run_anchors)
        self.align_panel.kindChanged.connect(self.set_align_kind)
        self.align_panel.propagateRequested.connect(self.propagate_anchors)
        self.align_panel.showAlignedToggled.connect(self.act_aligned.setChecked)
        # 🔍
        self.blobs.changed.connect(self._blobs_changed)
        self.blobs.changed.connect(self._hydrocarbons_changed)   # stencil totals depend on 🔍 blobs
        self.blobs.runningChanged.connect(self._blobs_running)
        self.blobs_panel.findRequested.connect(self.run_blobs)
        self.blobs_panel.methodChanged.connect(self._blob_method_changed)
        self.blobs_panel.resolveRequested.connect(self.resolve_region)
        self.blobs_panel.addRequested.connect(self.start_add_blob)
        self.blobs_panel.deleteRequested.connect(self.delete_blobs)
        self.blobs_panel.mergeRequested.connect(self.merge_blobs)
        self.blobs_panel.selectionChanged.connect(self._blob_selection_changed)
        self.blobs_panel.rowActivated.connect(self.goto_blob)
        # ⛓️
        self.hydrocarbons.changed.connect(self._hydrocarbons_changed)
        self.hydrocarbons.runningChanged.connect(self._hydrocarbons_running)
        self.hydro_panel.findRequested.connect(self.start_find_alkanes)
        self.hydro_panel.propagateRequested.connect(self.propagate_alkanes)
        self.hydro_panel.drawRequested.connect(self.view2d.start_stencil_draw)
        self.hydro_panel.cancelDrawRequested.connect(self.view2d.cancel_stencil_draw)
        self.hydro_panel.saveStencilRequested.connect(self.save_stencil)
        self.hydro_panel.removeStencilsRequested.connect(self.remove_stencils)
        self.hydro_panel.responseFactorEdited.connect(self.edit_response_factor)
        self.hydro_panel.settingsRequested.connect(self.edit_hydro_settings)
        self.hydro_panel.findHomologRequested.connect(self.find_homologs)
        self.hydro_panel.homologActivated.connect(self.goto_blob)
        self.hydro_panel.editStencilRequested.connect(self.edit_stencil)
        self.hydro_panel.qualifiersChanged.connect(self._qualifiers_edited)
        self.hydro_panel.pickQualifiersToggled.connect(self._qual_pick_toggled)
        self.hydro_panel.set_qualifiers(self.hydrocarbons.settings)
        # workflow
        self.workflow.goRequested.connect(self.go)
        self.workflow.stepChanged.connect(self._step_opened)
        self.workflow.progressChanged.connect(self._step_progress)
        self.quant_step.targetsRequested.connect(self.load_targets_dialog)
        self.quant_step.levelsRequested.connect(self.load_levels_dialog)
        self.quant_step.settingsRequested.connect(self.quant_settings_dialog)
        self.project.runsChanged.connect(self._quant_step_status)
        self.quant.changed.connect(self._quant_step_status)
        self.project.layersChanged.connect(self._quant_step_status)
        self._quant_step_status()
        self._step_progress()

    # =============================================================== runs
    def open_folder_dialog(self):
        folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Open a folder of runs")
        if folder:
            self.open_folder(folder)

    def open_folder(self, folder, on_loaded=None):
        self.say(f"Opening {folder}")
        self.project.open_folder(folder, on_loaded)

    def open_files_dialog(self):
        files, _ = QtWidgets.QFileDialog.getOpenFileNames(self, "Add runs — select energy or processed files explicitly", "", file_filter())
        if files:
            self.project.add_files(files)

    def load_targets_dialog(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Target table", "", "CSV (*.csv)")
        if path:
            try:
                self.project.load_targets(path, estimate=False)
            except (ValueError, KeyError) as exc:
                self.error("Couldn't read the target table", str(exc))
                return
            self.fix_phase()

    def load_sequence_dialog(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Sequence list", "", "CSV (*.csv)")
        if path:
            self.project.load_sequence(path)

    def remove_runs(self, rows):
        self.project.remove(rows)
        self.undo.clear()      # stored phase snapshots refer to run positions
        self.set_run(min(self.index, len(self.project) - 1))

    def _runs_changed(self):
        n = len(self.project)
        if n and self.index < 0:
            first_loaded = self.project.loaded()
            if first_loaded:
                self.set_run(first_loaded[0])
        self._update_run_label()

    def set_run(self, index: int):
        if not (0 <= index < len(self.project)):
            self.index = -1
            self._update_run_label()
            return
        self.index = index
        self.runs_panel.select(index)
        self.fold_panel.set_index(index)
        self.renderer.set_index(index)
        self._update_run_label()
        self._layers_changed()
        self._update_overlay()
        self._align_changed()
        self._blobs_changed()
        self._hydrocarbons_changed()
        if self.bg_panel.preview.isChecked():
            self.start_preview()
        if self._last_click is not None:
            QtCore.QTimer.singleShot(0, lambda: self.show_spectrum_at(self._last_click))

    def step_run(self, delta: int):
        loaded = self.project.loaded()
        if not loaded:
            return
        if self.index not in loaded:
            self.set_run(loaded[0])
            return
        pos = loaded.index(self.index) + delta
        if 0 <= pos < len(loaded):
            self.set_run(loaded[pos])

    def _update_run_label(self):
        n = len(self.project)
        if n == 0 or self.index < 0:
            self.run_label.setText("No runs" if n == 0 else f"{n} runs")
        else:
            self.run_label.setText(f"{self.project.entries[self.index].name}   run {self.index + 1} of {n}")

    # ================================================================ fold
    def edit_phase(self, value: float):
        if self.index >= 0:
            self.undo.push(SetPhase(self.project, self.index, value))

    def edit_period(self, value: float):
        if abs(value - self.project.period) > 1e-12:
            self.undo.push(SetPeriod(self.project, value))

    def estimate_period(self):
        if not self.project.loaded():
            self.say("Open some runs first")
            return

        def done(before, after):
            self.undo.push(SetFold(self.project, before, after, "🔄 Estimate period and phase"))
            self.workflow.mark_done("fold", "Period and phase estimated. Check the 2D view looks right.")

        self.project.auto_fold(on_done=done, force_period=True)

    def fix_phase(self):
        if not self.project.loaded():
            self.say("Open some runs first")
            return

        def done(before, after):
            self.undo.push(SetPhases(self.project, before, after, "🔄 Fix phase"))
            self.say("Phase estimated for every run")
            self.workflow.mark_done("fold", "Phase estimated for every run. Check the 2D view looks right.")

        self.project.estimate_phases(on_done=done)

    def copy_phase_to_all(self):
        if self.index < 0:
            return
        before = self.project.phase_snapshot()
        p = self.project.entries[self.index].phase
        after = {i: (p, "manual", f"copied from {self.project.entries[self.index].name}")
                 for i in before}
        self.project.set_phases(after)
        self.undo.push(SetPhases(self.project, before, self.project.phase_snapshot(),
                                 f"Phase {p:.3f} s for all runs"))

    # =============================================================== views
    def _rendered(self, result):
        self.view2d.show_result(result)
        self.view1d.show_result(result)
        band = self.view1d.band_range()
        if band:
            self.view2d.set_band(*band)

    def apply_ion_box(self):
        text = self.ion_box.currentText().strip() or "TIC"
        self.set_ion_expression(text)

    def set_ion_expression(self, text: str) -> bool:
        try:
            expr = self.renderer.set_expression(text)
        except IonExpressionError as exc:
            self.say(f"Ion expression: {exc}")
            self.ion_box.lineEdit().setStyleSheet("color: #FB923C;")
            return False
        self.ion_box.lineEdit().setStyleSheet("")
        self.ion_box.setEditText(expr.text())
        if expr.is_ratio and self.view2d.scale == "log":
            self.set_scale("linear")
        return True

    def ion_clicked(self, mz: int, add: bool):
        current = self.renderer.expr
        if add and current.numerator is not None and not current.is_ratio:
            text = "+".join(str(m) for m in (*current.numerator, mz))
        else:
            text = str(mz)
        self.set_ion_expression(text)

    def focus_ion_box(self):
        self.ion_box.setFocus()
        self.ion_box.lineEdit().selectAll()

    def set_scale(self, scale: str):
        self.scale_box.blockSignals(True)
        self.scale_box.setCurrentIndex(SCALES.index(scale))
        self.scale_box.blockSignals(False)
        self.view2d.set_scale(scale)

    def cycle_scale(self):
        self.set_scale(SCALES[(SCALES.index(self.view2d.scale) + 1) % len(SCALES)])
        self.say(f"Colour scale: {SCALE_NAMES[self.view2d.scale]}")

    def set_levels_locked(self, on: bool):
        self.view2d.levels_locked = on

    def _hovered(self, t1, t2, value):
        r = self.view2d.result
        what = "S/N" if (r is not None and r.snr) else self.renderer.expr.label()
        self.cursor_label.setText(f"t₁ {t1 / 60:.3f} min    t₂ {t2:.3f} s    {what}  {value:,.4g}")

    def _band_from_2d(self, a, b):
        self.view1d.set_band(a, b)

    def _band_from_1d(self, a, b):
        self.view2d.set_band(a, b)

    def _time_from_1d(self, t):
        self.show_spectrum_at(t)

    # ============================================================ spectra
    def _current_run(self):
        """The run as currently displayed (cleaned if 🧹 is on and applied)."""
        if self.index < 0 or self.project.entries[self.index].run is None:
            return None
        return self.project.display_run(self.index)

    def show_spectrum_at(self, t: float):
        run = self._current_run()
        if run is None:
            return
        self._last_click = t
        i = run.scan_at(t)
        lo, hi = max(i - 1, 0), min(i + 2, run.n_scans)
        spec = run.spectrum(slice(lo, hi))
        fold = self.project.fold(self.index)
        col, t2 = fold.to_fold(run.times[i])
        name = self.project.entries[self.index].name
        self.spectrum.set_spectrum(run.mz, spec,
                                   f"{name}: t₁ {fold.t1[col] / 60:.3f} min, t₂ {t2:.3f} s "
                                   f"(scans {lo}–{hi - 1})")

    def show_region_spectrum(self, t1a, t1b, t2a, t2b):
        run = self._current_run()
        if run is None:
            return
        fold = self.project.fold(self.index)
        P = self.project.period
        col, t2 = fold.to_fold(run.times)
        start = fold.t1
        c0 = int(np.searchsorted(start, t1a, side="right")) - 1
        c1 = int(np.searchsorted(start, t1b, side="right")) - 1
        in_cols = (col >= c0) & (col <= c1)
        mask = in_cols & (t2 >= t2a) & (t2 <= t2b)
        # Extended rows: time beyond one period belongs to the next modulation.
        prev = ((col - 1) >= c0) & ((col - 1) <= c1) & (t2 + P >= t2a) & (t2 + P <= t2b)
        mask |= prev
        n = int(mask.sum())
        if n == 0:
            self.spectrum.set_spectrum(run.mz, np.zeros(run.n_mz), "Region is empty")
            return
        spec = np.asarray(run.data[np.flatnonzero(mask)].sum(axis=0)).ravel()
        name = self.project.entries[self.index].name
        self.spectrum.set_spectrum(run.mz, spec,
                                   f"{name}: region t₁ {t1a / 60:.2f}–{t1b / 60:.2f} min, "
                                   f"t₂ {t2a:.2f}–{t2b:.2f} s ({n:,} scans)")

    # ========================================================= background
    def show_background_panel(self):
        self.workflow.show_step("background")

    def apply_background(self, indices=None):
        loaded = self.project.loaded()
        if not loaded:
            self.say("Open some runs first")
            return
        indices = loaded if indices is None else indices
        settings, blank = self.bg_panel.settings(), self.bg_panel.blank_name()
        before = self.project.background_snapshot()
        before_settings, before_blank = self.project.bg_settings, self.project.bg_blank
        n = len(indices)

        def done(results):
            after = dict(before)
            after.update(results)
            self.project.set_preview(None, None)
            self.project.set_backgrounds(results, settings, blank)
            label = f"🧹 Background on {n} run{'s' if n > 1 else ''}"
            self.undo.push(ApplyBackground(self.project, before, after, before_settings,
                                           settings, before_blank, blank, label))
            self.act_processed.setChecked(True)
            self.say(label)
            self.workflow.mark_done("background", f"{label}. Space toggles raw and cleaned.")

        self.project.compute_background(indices, settings, blank, on_done=done)

    def remove_background_all(self):
        before = self.project.background_snapshot()
        if not any(v is not None for v in before.values()):
            return
        after = {e: None for e in before}
        self.workflow.clear_done("background")
        self.project.set_backgrounds(after)
        self.undo.push(ApplyBackground(self.project, before, after, self.project.bg_settings,
                                       self.project.bg_settings, self.project.bg_blank,
                                       self.project.bg_blank, "Remove 🧹 from all runs"))

    def _bg_settings_changed(self):
        if self.bg_panel.preview.isChecked():
            self.start_preview()

    # ========================================================= 🫧 siloxanes
    def apply_siloxane(self, indices=None, _retry: bool = False):
        loaded = self.project.loaded()
        if not loaded:
            self.say("Open some runs first")
            return
        indices = loaded if indices is None else indices
        settings = self.project.siloxane_settings
        if settings.method == "blobs":
            # The blob method blocks out 🔍 blobs, so find them first if they are missing or out of date.
            entries = self.project.entries
            if any(self.blobs.blobs_for(i) is None or entries[i] in self.blobs.stale for i in indices):
                if _retry:
                    self.say("🫧 Couldn't find blobs on every run; see the Log tab")
                elif self.blobs.run_all(on_done=lambda: self.apply_siloxane(indices, _retry=True),
                                        method=LEGACY):
                    self.say("🫧 Finding blobs first…")
                else:
                    self.say("Blob detection is already running; press 🫧 again when it finishes")
                return
        before = self.project.siloxane_snapshot()
        n = len(indices)

        def done(results):
            if not results:
                return
            after = dict(before)
            after.update(results)
            self.project.set_siloxanes(results, settings)
            label = f"🫧 Siloxanes on {n} run{'s' if n > 1 else ''}"
            self.undo.push(ApplySiloxane(self.project, before, after, settings, settings, label))
            self.act_processed.setChecked(True)
            if settings.method == "blobs":
                blocked = sum(r.n_blobs for r in results.values())
                msg = f"{label}: blocked out {blocked} siloxane blob{'s' if blocked != 1 else ''}"
            else:
                removed = [r.removed_fraction for r in results.values() if r.ions_used]
                pct = f"{100 * sum(removed) / len(removed):.0f}% average" if removed else "nothing to remove"
                msg = f"{label} ({pct})"
            self.say(msg)
            self.workflow.mark_done("siloxane", msg)

        self.project.compute_siloxane(indices, settings, on_done=done)

    def _siloxane_settings_changed(self):
        """The reference or a threshold was edited in the 🫧 panel. Flags in
        the 🔍 panel follow at once; layers already applied are left as they
        are (and the panel says so) until 🫧 is applied again."""
        settings = self.sil_panel.settings()
        self.project.siloxane_settings = settings
        n_layers = sum(e.siloxane is not None for e in self.project.entries)
        self.sil_panel.set_stale_note(n_layers)
        # Dragging the certainty slider is only useful with the outlines showing.
        if (self.index >= 0 and self.blobs.blobs_for(self.index) is not None
                and not self.act_blob_overlay.isChecked()):
            self.act_blob_overlay.setChecked(True)
        self._blobs_changed()

    def _sync_siloxane_panel(self):
        """After a project opens or a layer is applied, show the settings in force."""
        if self.sil_panel.settings() != self.project.siloxane_settings:
            self.sil_panel.set_settings(self.project.siloxane_settings)

    def _siloxane_pick_toggled(self, on: bool):
        if on:
            if self.index < 0 or self.project.entries[self.index].run is None:
                self.sil_panel.set_picking(False)
                self.say("Select a loaded run first")
                return
            self.cancel_pick()
            self._siloxane_pick_pending = True
            self.sil_panel.set_picking(True)
            self.say("Click a siloxane peak in the 2D view (Esc cancels)")
        else:
            self._siloxane_pick_pending = False

    def _finish_siloxane_pick(self, t1: float, t2: float):
        self._siloxane_pick_pending = False
        self.sil_panel.set_picking(False)
        run = self.project.display_run(self.index)
        fold = self.project.fold(self.index)
        spec = spectrum_at(run, fold, t1, t2, subtract_baseline=self.sil_panel.subtract_baseline())
        ref = reference_from_spectrum(spec, run.mz)
        if len(ref) < self.project.siloxane_settings.min_ions_present:
            self.say("Too few ions stand out at that point. Click nearer the apex, or untick the baseline subtraction.")
            return
        base = max(ref, key=ref.get)
        self.sil_panel.set_reference(
            ref, f"Reference taken at t1 = {t1 / 60:.2f} min, t2 = {t2:.2f} s "
                 f"({len(ref)} ions, base m/z {base}). Edit the table to add or drop ions.")
        self.say(f"🫧 Reference set from the clicked point: {len(ref)} ions, base m/z {base}")

    def remove_siloxane_all(self):
        before = self.project.siloxane_snapshot()
        if not any(v is not None for v in before.values()):
            return
        after = {e: None for e in before}
        settings = self.project.siloxane_settings
        self.workflow.clear_done("siloxane")
        self.project.set_siloxanes(after)
        self.undo.push(ApplySiloxane(self.project, before, after, settings, settings,
                                     "Remove 🫧 from all runs"))

    def _preview_toggled(self, on: bool):
        if on:
            self.start_preview()
        else:
            self.project.set_preview(None, None)

    def start_preview(self):
        if self.index < 0 or self.project.entries[self.index].run is None:
            return
        entry = self.project.entries[self.index]
        settings, blank = self.bg_panel.settings(), self.bg_panel.blank_name()

        def done(results):
            still = (self.bg_panel.preview.isChecked() and 0 <= self.index < len(self.project)
                     and self.project.entries[self.index] is entry)
            if still and entry in results:
                self.project.set_preview(entry, results[entry])

        self.project.compute_background([self.index], settings, blank, on_done=done, preview=True)

    # ======================================================== quantification
    def run_quant(self):
        ok, why = self.quant.ready()
        if not ok:
            self.say(why)
            return
        if not self.project.any_background():
            self.say("🧪 Quantifying raw data. Apply 🧹 first for background-corrected volumes and S/N.")
        if not self.quant.levels:
            self.say("No calibration levels loaded, so curves can't be fitted (File → Load calibration levels)")
        def done():
            self.say("🧪 Quantification finished")
            self.workflow.mark_done("quant", "Quantification finished. Review it in the Review tab (3) "
                                             "and Calibration tab (4), then export from 🧪 Results.")

        self.quant.run(on_done=done)

    def quant_settings_dialog(self):
        dlg = QuantSettingsDialog(self.quant.settings, self)
        if dlg.exec() == QtWidgets.QDialog.DialogCode.Accepted:
            self.quant.set_settings(dlg.settings())

    def load_levels_dialog(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Calibration levels", "", "CSV (*.csv)")
        if path:
            try:
                self.quant.set_levels(path)
            except (ValueError, KeyError, OSError) as exc:
                self.error("Couldn't read the calibration levels", str(exc))

    def export_results(self, folder=None):
        if self.quant.report is None:
            self.say("Nothing to export yet: run 🧪 Quantify first")
            return
        if not folder:
            folder = QtWidgets.QFileDialog.getExistingDirectory(self, "Export 🧪 results to")
        if folder:
            paths = self.quant.report.write(folder)
            if any(e.import_meta for e in self.project.entries):
                from ..core.report import build_provenance
                prov = build_provenance(__version__, period=self.project.period,
                    method={'modulation_period_confirmed': self.project.period_confirmed},
                    settings={'quant': self.quant.settings.to_dict(), 'background': self.project.bg_settings.to_dict()},
                    quant_decisions=self.quant.decisions,
                    runs=[e.run for e in self.project.entries if e.run is not None])
                paths.append(prov.write(Path(folder) / 'provenance.json'))
            self.say(f"Wrote {', '.join(p.name for p in paths)} to {folder}")

    def _quant_running(self, on: bool):
        self.act_quant.setEnabled(not on)
        self.workflow.set_busy("quant", on)
        self._quant_changed()

    def _quant_changed(self):
        q = self.quant
        self.results.set_report(q.report, stale=q.stale, running=q.running)
        self.calibration_panel.set_report(q.report, q.decisions)
        self.review.set_report(q.report)
        self._update_overlay()

    def _update_overlay(self):
        """Outlines of every peak found in this run, and labelled markers for
        targets (lilac) and internal standards (teal); red when rejected."""
        rq = self.quant.run_quant(self.index) if self.index >= 0 else None
        if rq is None or self.quant.report is None:
            self.view2d.set_overlay([], (np.array([]), np.array([])))
            return
        xs, ys = [], []
        blobs = rq.blobs
        for i in range(len(blobs)):
            x, y = blobs.outline_times(i)
            xs.append(x)
            ys.append(y)
        outline = (np.concatenate(xs), np.concatenate(ys)) if xs else (np.array([]), np.array([]))
        status = {(r.run, r.compound): r.status for r in self.quant.report.results}
        markers = []
        for name, hit in rq.hits.items():
            if not hit.found:
                continue
            st = status.get((rq.run, name), "")
            if st == "rejected":
                colour = (194, 65, 12)
            elif hit.role == "istd":
                colour = (92, 198, 188)
            else:
                colour = (209, 149, 224) if st != "check" else (224, 160, 32)
            markers.append({"t1": hit.t1, "t2": hit.t2, "label": name + (" (IS)" if hit.role == "istd" else ""),
                            "colour": colour, "symbol": "s" if hit.role == "istd" else "o",
                            "key": name})
        self.view2d.set_overlay(markers, outline)

    def _hit_time(self, index: int, hit) -> float:
        fold = self.project.fold(index)
        col = int(np.floor((hit.t1 - fold.t1[0]) / fold.spec.period))
        return float(fold.t1[col] + hit.t2)

    def goto_peak(self, run: str, compound: str):
        i = self.project.index_of(run)
        if i is None:
            return
        self.set_run(i)
        self.tabs.setCurrentWidget(self.view2d)
        target = next((t for t in (self.quant.used_targets() or []) if t.compound == compound), None)
        if target is not None:
            self.set_ion_expression(str(target.quant_ion))
        rq = self.quant.run_quant(i)
        hit = rq.hits.get(compound) if rq else None
        if hit is not None and hit.found:
            t1, t2 = hit.t1, hit.t2
        elif hit is not None and hit.expected:
            t1, t2 = hit.expected
        elif target is not None:
            t1, t2 = target.t1, target.t2
        else:
            return

        def centre():
            self.view2d.centre_on(t1, t2)
            if hit is not None and hit.found:
                self.show_spectrum_at(self._hit_time(i, hit))
            self.results.select_target(run, compound)

        QtCore.QTimer.singleShot(60, self, centre)

    def _goto_row(self, row: dict):
        if row.get("compound") and row.get("run"):
            self.goto_peak(row["run"], row["compound"])
        elif row.get("run") is not None and row.get("t1") is not None:
            i = self.project.index_of(row["run"])
            if i is not None:
                self.set_run(i)
                self.tabs.setCurrentWidget(self.view2d)
                QtCore.QTimer.singleShot(60, self, lambda: self.view2d.centre_on(row["t1"], row["t2"]))

    def _edit(self, text: str, change, requantify=()):
        before = self.quant.decisions.copy()
        after = before.copy()
        change(after)
        self.undo.push(EditDecisions(self.quant, before, after, text, requantify))

    def decide(self, pairs: list, status: str):
        if not pairs:
            return

        def change(d):
            for pair in pairs:
                pair = tuple(pair)
                d.accepted.discard(pair)
                d.rejected.discard(pair)
                if status == "accepted":
                    d.accepted.add(pair)
                elif status == "rejected":
                    d.rejected.add(pair)

        what = {"accepted": "Accept", "rejected": "Reject", "clear": "Clear decision on"}[status]
        label = f"{pairs[0][1]} in {pairs[0][0]}" if len(pairs) == 1 else f"{len(pairs)} peaks"
        self._edit(f"🧪 {what} {label}", change)

    def override_curve(self, compound: str, overrides: dict):
        def change(d):
            if overrides:
                d.overrides[compound] = dict(overrides)
            else:
                d.overrides.pop(compound, None)
        self._edit(f"🧪 Curve settings for {compound}", change)

    def exclude_point(self, compound: str, run: str, reason: str):
        def change(d):
            ex = d.excluded.setdefault(compound, {})
            if reason:
                ex[run] = reason
            else:
                ex.pop(run, None)
                if not ex:
                    d.excluded.pop(compound, None)
        verb = f"Exclude {run} from" if reason else f"Include {run} in"
        self._edit(f"🧪 {verb} the {compound} curve", change)

    def start_pick(self, run: str, compound: str):
        self.goto_peak(run, compound)
        self._picking = (run, compound)
        self.say(f"Click the peak for {compound} in {run} (Esc cancels)")

    # ---------------------------------------------- 🧪 hand-drawn footprints
    def reintegrate_current(self):
        tiles = self.review.tiles
        if not tiles:
            self.say("Open the 🧪 Review gallery and select a peak first")
            return
        t = tiles[self.review.current]
        self.start_reintegrate(t.run, t.compound)

    def start_reintegrate(self, run: str, compound: str):
        """Show the peak's integration footprint as an orange polygon on the 2D
        view (its current hand-drawn one, else the bounding box of the automatic
        one) to drag into shape; Ctrl+Enter integrates inside it."""
        i = self.project.index_of(run)
        rq = self.quant.run_quant(i) if i is not None else None
        hit = rq.hits.get(compound) if rq else None
        if hit is None or not hit.found:
            self.say(f"{compound} has no peak in {run} to re-integrate: pick one first (P)")
            return
        self.goto_peak(run, compound)
        fold = self.project.fold(i)
        poly = self.quant.decisions.footprints.get((run, compound))
        if not poly:
            rows, cols = hit.pixels
            t1a, t1b = float(fold.t1[int(cols.min())]), float(fold.t1[int(cols.max())] + fold.spec.period)
            t2a, t2b = (float(rows.min()) - 0.5) * fold.row_pitch, (float(rows.max()) + 0.5) * fold.row_pitch
            poly = [(t1a, t2a), (t1b, t2a), (t1b, t2b), (t1a, t2b)]
        self._reintegrate = (run, compound)
        self.act_accept_footprint.setEnabled(True)
        self.act_reset_footprint.setEnabled(True)
        QtCore.QTimer.singleShot(120, self, lambda: self.view2d.start_footprint_edit(poly))
        self.say(f"Drag the corners of the orange polygon around {compound}, then Ctrl+Enter "
                 "(Esc cancels; Process menu: ↺ Automatic footprint)")

    def _end_reintegrate(self):
        self._reintegrate = None
        self.view2d.cancel_footprint_edit()
        self.act_accept_footprint.setEnabled(False)
        self.act_reset_footprint.setEnabled(False)

    def accept_footprint(self):
        if self._reintegrate is None:
            return
        run, compound = self._reintegrate
        poly = self.view2d.footprint_vertices()
        if len(poly) < 3:
            return

        def change(d):
            d.footprints[(run, compound)] = poly
            d.rejected.discard((run, compound))

        self._end_reintegrate()
        self._edit(f"🧪 Integrate {compound} in {run} by hand", change, requantify=(run,))
        self.say(f"Re-measuring {compound} in {run} inside the drawn footprint")

    def reset_footprint(self):
        if self._reintegrate is None:
            return
        run, compound = self._reintegrate
        self._end_reintegrate()
        if (run, compound) in self.quant.decisions.footprints:
            self._edit(f"🧪 Automatic footprint for {compound} in {run}",
                       lambda d: d.footprints.pop((run, compound), None), requantify=(run,))
        self.say(f"{compound} in {run} back to the automatic footprint")

    def cancel_pick(self):
        if self._anchor_stage is not None:
            self._anchor_stage = None
            self._anchor_pending = {}
            self.say("Anchor picking cancelled")
        elif self._blob_add_pending:
            self._blob_add_pending = False
            self.say("Adding a blob cancelled")
        elif self._alkane_seed_pending is not None:
            self._alkane_seed_pending = None
            self.say("Finding the alkane series cancelled")
        elif self._reintegrate is not None:
            self._end_reintegrate()
            self.say("Re-integration cancelled")
        elif self._siloxane_pick_pending:
            self._siloxane_pick_pending = False
            self.sil_panel.set_picking(False)
            self.say("Picking a siloxane spectrum cancelled")
        elif self._qual_pick_pending:
            self._qual_pick_pending = False
            self.hydro_panel.set_picking_qualifiers(False)
            self.say("Picking the alkane qualifier ratios cancelled")
        elif self._picking is not None:
            self._picking = None
            self.say("Picking cancelled")

    def _point_clicked(self, col: int, t2: float):
        if self._anchor_stage is not None and self.index >= 0:
            fold = self.project.fold(self.index)
            t1 = float(fold.t1[col] + 0.5 * fold.spec.period)
            if self._anchor_stage == "ref":
                self._finish_anchor_ref_click(t1, t2)
            else:
                self._finish_anchor_run_click(t1, t2)
            return
        if self._blob_add_pending and self.index >= 0:
            fold = self.project.fold(self.index)
            t1 = float(fold.t1[col] + 0.5 * fold.spec.period)
            self._finish_add_blob(t1, t2)
            return
        if self._alkane_seed_pending is not None and self.index >= 0:
            fold = self.project.fold(self.index)
            t1 = float(fold.t1[col] + 0.5 * fold.spec.period)
            self._finish_alkane_click(t1, t2)
            return
        if self._siloxane_pick_pending and self.index >= 0:
            fold = self.project.fold(self.index)
            self._finish_siloxane_pick(float(fold.t1[col] + 0.5 * fold.spec.period), t2)
            return
        if self._qual_pick_pending and self.index >= 0:
            fold = self.project.fold(self.index)
            self._finish_qualifier_pick(float(fold.t1[col] + 0.5 * fold.spec.period), t2)
            return
        if self._picking is not None and self.index >= 0:
            run, compound = self._picking
            if self.project.entries[self.index].name != run:
                return
            self._picking = None
            fold = self.project.fold(self.index)
            t1_mid = float(fold.t1[col] + 0.5 * fold.spec.period)

            def change(d):
                d.assigned[(run, compound)] = (t1_mid, float(t2))
                d.rejected.discard((run, compound))

            self._edit(f"🧪 Pick {compound} in {run} by hand", change, requantify=(run,))
            self.say(f"Re-measuring {compound} in {run}")
            return
        self._select_blob_at(col, t2)

    def _select_blob_at(self, col: int, t2: float):
        """Select the nearest centre for the active peak-detection method."""
        if self.index < 0 or not self.view2d.blob_overlay_on:
            return
        if self.blobs.method == MS_AWARE:
            component_set = self.blobs.components_for(self.index)
            if component_set is not None and component_set.components:
                fold = self.project.fold(self.index)
                t1 = float(fold.t1[col] + 0.5 * fold.spec.period)
                components = component_set.components
                distance = np.asarray([
                    np.hypot((component.t1 - t1) / 60.0, (component.t2 - t2) / 0.5)
                    for component in components])
                nearest = int(np.argmin(distance))
                if distance[nearest] < 1.5:
                    self.blobs_panel.select_component(components[nearest].id)
                    self._update_blob_overlay()
                return
            candidates = self.blobs.candidates_for(self.index)
            regions = self.blobs.regions_for(self.index)
            if candidates is None or not regions:
                return
            fold = self.project.fold(self.index)
            t1 = float(fold.t1[col] + 0.5 * fold.spec.period)
            distance = np.hypot((candidates.t1 - t1) / 60.0, (candidates.t2 - t2) / 0.5)
            candidate_index = int(np.argmin(distance))
            if distance[candidate_index] < 1.5:
                region_index = next((i for i, region in enumerate(regions)
                                     if candidate_index in region.candidate_indices), None)
                if region_index is not None:
                    self.blobs_panel.select_region(region_index)
                    self._update_blob_overlay()
            return
        blobs = self.blobs.blobs_for(self.index)
        if blobs is None or len(blobs) == 0:
            return
        fold = self.project.fold(self.index)
        t1 = float(fold.t1[col] + 0.5 * fold.spec.period)
        d = np.hypot((blobs.t1 - t1) / 60.0, (blobs.t2 - t2) / 0.5)
        i = int(np.argmin(d))
        if d[i] < 1.5:      # close enough that this was probably meant for it
            self.blobs_panel.select_row(i)
            self._update_blob_overlay()

    def _review_image(self, run: str, ion: int):
        """Extended quantifier-ion image and frame for the review gallery."""
        i = self.project.index_of(run)
        if i is None or self.project.entries[i].run is None:
            return None
        e = self.project.entries[i]
        clean = e.background.run if e.background is not None else e.run
        fold = self.project.fold(i)
        key = (id(clean), id(fold), ion, self.quant.settings.peaks)
        hit = self._frame_cache.get(key)
        if hit is None:
            frame = make_frame(fold, self.quant.settings.peaks)
            hit = (frame.extend(fold.ion_image(clean, ion)), frame)
            if len(self._frame_cache) > 400:
                self._frame_cache.clear()
            self._frame_cache[key] = hit
        return hit

    # ================================================================ 📍
    def show_align_panel(self):
        self.workflow.show_step("align")

    def _go_blink(self):
        self.tabs.setCurrentWidget(self.compare)
        i = self.compare.mode.findData("blink")
        if i >= 0:
            self.compare.mode.setCurrentIndex(i)

    def _current_single_ion(self) -> int:
        expr = self.renderer.expr
        if expr.numerator and len(expr.numerator) == 1 and not expr.is_ratio:
            return int(expr.numerator[0])
        return 0

    def set_reference(self, name: str):
        name = name or None
        if name == self.align.reference:
            return
        before = self.align.snapshot()
        self.align.set_reference(name)
        label = f"📍 Reference: {name}" if name else "📍 Clear reference"
        self.undo.push(SetReference(self.align, before, name, label))

    def start_add_anchor(self):
        if self.index < 0:
            self.say("Open and select a run first")
            return
        if not self.align.reference:
            self.say("Pick a reference run first (📍 Align)")
            return
        target = self.project.entries[self.index].name
        if target == self.align.reference:
            self.say("Select a run other than the reference to align it")
            return
        self._anchor_pending = {"target": target, "name": self.align.next_anchor_name(),
                                "prev_index": self.index, "ion": self._current_single_ion()}
        self._anchor_stage = "ref"
        self.set_run(self.project.index_of(self.align.reference))
        self.say(f"Click the peak for {self._anchor_pending['name']} on the reference run "
                f"{self.align.reference} (Esc cancels)")

    def _finish_anchor_ref_click(self, t1: float, t2: float):
        p = self._anchor_pending
        p["ref_t1"], p["ref_t2"] = t1, t2
        self._anchor_stage = "run"
        self.set_run(self.project.index_of(p["target"]))
        self.say(f"Now click the same peak on {p['target']} (Esc cancels)")

    def _finish_anchor_run_click(self, t1: float, t2: float):
        p = self._anchor_pending
        before = self.align.snapshot()
        self.align.add_anchor(p["target"], p["name"], p["ref_t1"], p["ref_t2"], t1, t2, p["ion"])
        after = self.align.snapshot()
        self.undo.push(EditAlignment(self.align, before, after,
                                     f"📍 Add anchor {p['name']} on {p['target']}"))
        self.say(f"Anchor {p['name']} added on {p['target']}")
        self._anchor_stage = None
        self._anchor_pending = {}

    def remove_anchor(self, name: str):
        if self.index < 0:
            return
        run = self.project.entries[self.index].name
        before = self.align.snapshot()
        self.align.remove_anchor(run, name)
        self.undo.push(EditAlignment(self.align, before, self.align.snapshot(),
                                     f"📍 Remove anchor {name} from {run}"))

    def clear_run_anchors(self):
        if self.index < 0:
            return
        run = self.project.entries[self.index].name
        pts = list(self.align.points.get(run, []))
        if not pts:
            return
        before = self.align.snapshot()
        for p in pts:
            self.align.remove_anchor(run, p.name)
        self.undo.push(EditAlignment(self.align, before, self.align.snapshot(),
                                     f"📍 Clear anchors on {run}"))

    def set_align_kind(self, kind: str):
        if self.index < 0:
            return
        run = self.project.entries[self.index].name
        if self.align.kind.get(run, "auto") == kind:
            return
        before = self.align.snapshot()
        self.align.set_kind(run, kind)
        self.undo.push(EditAlignment(self.align, before, self.align.snapshot(),
                                     f"📍 {kind} transform for {run}"))

    def propagate_anchors(self):
        loaded = self.project.loaded()
        if not loaded:
            self.say("Open some runs first")
            return
        if not self.align.reference:
            self.say("Pick a reference run first (📍 Align)")
            return
        self.align.propagate_to(loaded, on_done=lambda: self.workflow.mark_done(
            "align", "Anchors searched on every run. Check the residuals, or the Compare tab (5)."))

    def _compare_image(self, index: int):
        if index is None or not (0 <= index < len(self.project.entries)) or not self.align.reference:
            return None
        img = self.align.image_for(index, self.renderer.expr, extended=False)
        if img is None:
            return None
        ref_i = self.project.index_of(self.align.reference)
        ref_fold = self.project.fold(ref_i)
        return img, ref_fold.t1, ref_fold.row_pitch, ref_fold.spec.period, self.project.entries[index].name

    def _align_changed(self):
        self.align_panel.refresh_runs(self.project.entries)
        self.align_panel.set_current_run(self.project.entries[self.index].name if self.index >= 0 else None)
        ref_i = self.project.index_of(self.align.reference) if self.align.reference else -1
        self.compare.set_runs(self.index, ref_i if ref_i is not None else -1)
        has = self.index >= 0 and self.renderer.aligned_available(self.index)
        self.act_aligned.setEnabled(has)
        if not has and self.act_aligned.isChecked():
            self.act_aligned.setChecked(False)

    # ================================================================ 🔍
    def show_blobs_panel(self):
        self.workflow.show_step("blobs")

    def run_blobs(self):
        if not self.project.loaded():
            self.say("Open some runs first")
            return

        def done():
            if self.blobs.method == MS_AWARE:
                self.say("🔍 MS-aware peak detection finished")
                self.workflow.mark_done(
                    "blobs", "MS-aware peak regions found. Select a region to inspect or resolve it.")
            else:
                self.say("🔍 Legacy blob detection finished")
                self.workflow.mark_done("blobs", "Legacy blobs found on every run. Check the table.")

        self.blobs.run_all(on_done=done)

    def _blob_method_changed(self, method: str):
        self.blobs.set_method(method)
        label = "MS-aware separation" if method == MS_AWARE else "legacy blob borders"
        self.say(f"🔍 Peak-detection method: {label}")

    def _blobs_running(self, on: bool):
        self.act_find_blobs.setEnabled(not on)
        self.workflow.set_busy("blobs", on)
        self._blobs_changed()

    def _blobs_changed(self):
        name = self.project.entries[self.index].name if self.index >= 0 else None
        entry = self.project.entries[self.index] if self.index >= 0 else None
        blobs = self.blobs.by_run.get(entry) if entry is not None else None
        candidates = self.blobs.candidates_by_run.get(entry) if entry is not None else None
        regions = self.blobs.regions_by_run.get(entry) if entry is not None else None
        component_set = self.blobs.components_for(self.index) if self.index >= 0 else None
        stale = self.blobs.selected_stale(entry) if entry is not None else False
        sil = self.project.siloxane_settings
        scores = siloxane_scores(blobs, sil) if blobs is not None and len(blobs) else None
        self.blobs_panel.set_current_run(
            name, blobs, stale, self.blobs.running, scores, sil.score_threshold,
            method=self.blobs.method, candidates=candidates, regions=regions,
            component_set=component_set)
        self.sil_panel.set_blob_count(None if scores is None else int((scores >= sil.score_threshold).sum()),
                                      0 if scores is None else len(scores))
        self._update_blob_overlay()

    def _update_blob_overlay(self):
        if self.blobs.method == MS_AWARE:
            candidates = self.blobs.candidates_for(self.index) if self.index >= 0 else None
            regions = self.blobs.regions_for(self.index) if self.index >= 0 else None
            if candidates is None or regions is None:
                self.view2d.set_candidate_overlay(
                    (np.array([]), np.array([])), (np.array([]), np.array([])))
                self.view2d.set_component_overlay([])
                return
            component_set = self.blobs.components_for(self.index)
            accepted_regions = ({group.region_index for group in component_set.groups}
                                if component_set is not None else set())
            outlines = {i: self._candidate_region_outline(candidates, region)
                        for i, region in enumerate(regions) if i not in accepted_regions}
            outline = (np.concatenate([x for x, _ in outlines.values()]),
                       np.concatenate([y for _, y in outlines.values()])) if outlines else (
                           np.array([]), np.array([]))
            selected = [outlines[i] for i in self.blobs_panel.selected_region_indices()
                        if i in outlines]
            highlight = (np.concatenate([x for x, _ in selected]),
                         np.concatenate([y for _, y in selected])) if selected else (
                             np.array([]), np.array([]))
            accepted_candidate_indices = {
                candidate_index for i, region in enumerate(regions) if i in accepted_regions
                for candidate_index in region.candidate_indices}
            keep = np.asarray([i not in accepted_candidate_indices for i in range(len(candidates))])
            self.view2d.set_candidate_overlay(
                outline, (candidates.t1[keep], candidates.t2[keep]), highlight)
            colours = ((0, 220, 255), (255, 105, 180), (255, 205, 60),
                       (120, 235, 120), (190, 130, 255), (255, 140, 70))
            selected_regions = set(self.blobs_panel.selected_region_indices())
            component_overlays = []
            if component_set is not None:
                for group in component_set.groups:
                    for component_index, _component in enumerate(group.components):
                        colour = colours[component_index % len(colours)]
                        for fraction, style in ((0.10, QtCore.Qt.PenStyle.SolidLine),
                                                (0.50, QtCore.Qt.PenStyle.DotLine)):
                            contour = group.contour(component_index, fraction)
                            component_overlays.append({
                                "outline": contour, "colour": colour, "style": style,
                                "width": 3.0 if group.region_index in selected_regions else 1.8})
            self.view2d.set_component_overlay(component_overlays)
            return
        blobs = self.blobs.blobs_for(self.index) if self.index >= 0 else None
        if blobs is None or len(blobs) == 0:
            self.view2d.set_blob_overlay((np.array([]), np.array([])), (np.array([]), np.array([])))
            return
        xs, ys = [], []
        for i in range(len(blobs)):
            x, y = blobs.outline_times(i)
            xs.append(x)
            ys.append(y)
        outline = (np.concatenate(xs), np.concatenate(ys))
        sel = self.blobs_panel.selected_indices()
        if sel:
            xs, ys = [], []
            for i in sel:
                if 0 <= i < len(blobs):
                    x, y = blobs.outline_times(i)
                    xs.append(x)
                    ys.append(y)
            highlight = (np.concatenate(xs), np.concatenate(ys)) if xs else (np.array([]), np.array([]))
        else:
            highlight = (np.array([]), np.array([]))
        flagged = np.flatnonzero(siloxane_scores(blobs, self.project.siloxane_settings)
                                 >= self.project.siloxane_settings.score_threshold)
        xs, ys = [], []
        for i in flagged:
            x, y = blobs.outline_times(int(i))
            xs.append(x)
            ys.append(y)
        siloxane = (np.concatenate(xs), np.concatenate(ys)) if xs else (np.array([]), np.array([]))
        self.view2d.set_blob_overlay(outline, highlight, siloxane)

    @staticmethod
    def _candidate_region_outline(candidates, region):
        x0 = candidates.t1_start + region.col0 * candidates.period
        x1 = candidates.t1_start + region.col1 * candidates.period
        y0, y1 = region.row0 * candidates.pitch, region.row1 * candidates.pitch
        x = np.array([x0, x1, np.nan, x1, x1, np.nan, x1, x0, np.nan, x0, x0, np.nan])
        y = np.array([y0, y0, np.nan, y0, y1, np.nan, y1, y1, np.nan, y1, y0, np.nan])
        return x, y

    def start_add_blob(self):
        if self.index < 0 or self.project.entries[self.index].run is None:
            self.say("Select a loaded run first")
            return
        if self.blobs.blobs_for(self.index) is None:
            self.say("🔍 Find blobs on this run first")
            return
        self.cancel_pick()
        self._blob_add_pending = True
        self.say("Click the peak to add in the 2D view (Esc cancels)")

    def _finish_add_blob(self, t1: float, t2: float):
        self._blob_add_pending = False
        before = self.blobs.snapshot_run(self.index)
        try:
            self.blobs.add(self.index, t1, t2)
        except ValueError as exc:
            self.say(f"Couldn't add a blob there: {exc}")
            return
        after = self.blobs.snapshot_run(self.index)
        self.undo.push(EditBlobs(self.blobs, before, after, "🔍 Add a blob by hand"))

    def delete_blobs(self, indices: list[int]):
        if self.index < 0 or not indices:
            return
        before = self.blobs.snapshot_run(self.index)
        self.blobs.delete(self.index, indices)
        after = self.blobs.snapshot_run(self.index)
        n = len(indices)
        self.undo.push(EditBlobs(self.blobs, before, after, f"🔍 Delete {n} blob{'s' if n != 1 else ''}"))

    def merge_blobs(self, indices: list[int]):
        if self.index < 0 or len(indices) < 2:
            return
        before = self.blobs.snapshot_run(self.index)
        try:
            self.blobs.merge(self.index, indices)
        except ValueError as exc:
            self.say(str(exc))
            return
        after = self.blobs.snapshot_run(self.index)
        self.undo.push(EditBlobs(self.blobs, before, after, f"🔍 Merge {len(indices)} blobs"))

    def merge_selected_blobs(self):
        if self.blobs.method != LEGACY:
            self.say("MS-aware regions may overlap and are not merged as hard borders")
            return
        self.merge_blobs(self.blobs_panel.selected_indices())

    def _blob_selection_changed(self, indices: list[int]):
        self._update_blob_overlay()

    def goto_blob(self, row: int):
        if self.blobs.method == MS_AWARE:
            component = self.blobs_panel.component_for_row(row)
            if component is not None:
                self.tabs.setCurrentIndex(0)
                self.view2d.centre_on(float(component.t1), float(component.t2))
                return
            candidates = self.blobs.candidates_for(self.index) if self.index >= 0 else None
            regions = self.blobs.regions_for(self.index) if self.index >= 0 else None
            region_row = self.blobs_panel.region_for_row(row)
            if (candidates is None or regions is None or region_row is None
                    or not (0 <= region_row < len(regions))):
                return
            indices = np.asarray(regions[region_row].candidate_indices, dtype=int)
            strongest = int(indices[np.argmax(candidates.score[indices])])
            self.tabs.setCurrentIndex(0)
            self.view2d.centre_on(float(candidates.t1[strongest]), float(candidates.t2[strongest]))
            return
        blobs = self.blobs.blobs_for(self.index) if self.index >= 0 else None
        if blobs is None or not (0 <= row < len(blobs)):
            return
        self.tabs.setCurrentIndex(0)
        self.view2d.centre_on(float(blobs.t1[row]), float(blobs.t2[row]))

    def resolve_region(self, row: int):
        self.show_overlap_dialog(candidate_region=row)

    def accept_component_group(self, run_name: str, group) -> None:
        index = self.project.index_of(run_name)
        if index is None:
            return
        before = self.blobs.snapshot_components(index)
        self.blobs.accept_component_group(index, group)
        after = self.blobs.snapshot_components(index)
        self.undo.push(EditComponents(
            self.blobs, before, after,
            f"🔍 Accept {len(group.components)} overlapping component"
            f"{'s' if len(group.components) != 1 else ''}"))
        state = "unresolved" if group.unresolved else "provisionally resolved"
        self.say(f"🔍 Accepted {len(group.components)} overlapping components ({state})")

    # ================================================================ ⛓️
    def show_hydrocarbons_panel(self):
        self.workflow.show_step("hydro")

    def start_find_alkanes(self):
        if self.index < 0 or self.project.entries[self.index].run is None:
            self.say("Select a loaded run first")
            return
        self.cancel_pick()
        self._alkane_seed_pending = []
        self.say("Click the first n-alkane in the 2D view (Esc cancels)")

    def _finish_alkane_click(self, t1: float, t2: float):
        carbon, ok = QtWidgets.QInputDialog.getInt(
            self, "n-Alkane", "Carbon number of the compound you just clicked:", 12, 1, 100)
        if not ok:
            self._alkane_seed_pending = None
            self.say("Finding the alkane series cancelled")
            return
        self._alkane_seed_pending.append((carbon, t1, t2))
        if len(self._alkane_seed_pending) < 2:
            self.say(f"C{carbon} noted. Click a second n-alkane (Esc cancels)")
            return
        seeds = self._alkane_seed_pending
        self._alkane_seed_pending = None
        index = self.index

        def done(series):
            self.say(f"⛓️ {len(series)} alkanes confirmed on {self.project.entries[index].name}")
            self.workflow.mark_done("hydro", f"{len(series)} alkanes confirmed. Propagate them to the "
                                             "other runs if you have more than one, then go on.")

        self.hydrocarbons.find_series(index, seeds, on_done=done)

    def propagate_alkanes(self):
        if self.index < 0:
            self.say("Select the run you found the alkane series on first")
            return
        loaded = self.project.loaded()

        def done():
            self.say("⛓️ Alkane series propagated")

        self.hydrocarbons.propagate_series(self.index, loaded, on_done=done)

    def save_stencil(self, name: str, response_factor: float):
        vertices = self.view2d.stencil_vertices()
        before = self.hydrocarbons.snapshot_stencils()
        if not self.hydrocarbons.add_stencil(self.index, vertices, name, response_factor):
            return
        after = self.hydrocarbons.snapshot_stencils()
        self.undo.push(EditHydrocarbons(self.hydrocarbons, before, after, f"⛓️ Save stencil {name}"))
        self.view2d.cancel_stencil_draw()
        self.hydro_panel.set_drawing(False)

    def remove_stencils(self, names: list[str]):
        if not names:
            return
        before = self.hydrocarbons.snapshot_stencils()
        self.hydrocarbons.remove_stencils(names)
        after = self.hydrocarbons.snapshot_stencils()
        n = len(names)
        self.undo.push(EditHydrocarbons(self.hydrocarbons, before, after,
                                        f"⛓️ Delete {n} stencil{'s' if n != 1 else ''}"))

    def edit_response_factor(self, name: str, value: float):
        before = self.hydrocarbons.snapshot_stencils()
        self.hydrocarbons.set_response_factor(name, value)
        after = self.hydrocarbons.snapshot_stencils()
        self.undo.push(EditHydrocarbons(self.hydrocarbons, before, after, f"⛓️ Response factor for {name}"))

    def _hydrocarbons_running(self, on: bool):
        self.act_find_alkanes.setEnabled(not on)
        self._hydrocarbons_changed()

    def _qualifiers_edited(self):
        ratios, tol, enforce = self.hydro_panel.qualifiers()
        st = self.hydrocarbons
        st.settings = replace(st.settings, qual_ratios=ratios, qual_tol=tol, enforce_qualifiers=enforce)
        self._hydrocarbons_changed()
        self.say("⛓️ Qualifier ratios updated. Find or propagate the alkane series again to use them.")

    def _qual_pick_toggled(self, on: bool):
        if on:
            if self.index < 0 or self.project.entries[self.index].run is None:
                self.hydro_panel.set_picking_qualifiers(False)
                self.say("Select a loaded run first")
                return
            self.cancel_pick()
            self._qual_pick_pending = True
            self.hydro_panel.set_picking_qualifiers(True)
            self.say("Click an n-alkane peak in the 2D view (Esc cancels)")
        else:
            self._qual_pick_pending = False

    def _finish_qualifier_pick(self, t1: float, t2: float):
        self._qual_pick_pending = False
        self.hydro_panel.set_picking_qualifiers(False)
        st = self.hydrocarbons
        run, fold = self.project.display_run(self.index), self.project.fold(self.index)
        try:
            ratios = qualifier_ratios_at(run, fold, t1, t2, st.settings)
        except KeyError:
            self.say("The quantifier or qualifier ions are outside this run's m/z range")
            return
        if not any(r > 0 for r in ratios):
            self.say("No qualifier signal at that point. Click nearer the apex of the alkane.")
            return
        rounded = tuple(round(r, 3) for r in ratios)
        st.settings = replace(st.settings, qual_ratios=rounded)
        self._hydrocarbons_changed()
        txt = ", ".join(f"{i}/{st.settings.quant_ion} = {r:.2f}" for i, r in zip(st.settings.qual_ions, rounded))
        self.say(f"⛓️ Expected ratios set from the clicked peak: {txt}")

    def edit_blob_settings(self):
        ms_aware = self.blobs.method == MS_AWARE
        settings = self.blobs.component_settings if ms_aware else self.blobs.settings
        title = "🔍 MS-aware detection settings" if ms_aware else "🔍 Legacy blob settings"
        dlg = SettingsDialog(title, settings, parent=self)
        if dlg.exec():
            if ms_aware:
                self.blobs.set_component_settings(dlg.settings())
                self.say("🔍 MS-aware settings changed. Find peaks again to use them.")
            else:
                self.blobs.set_settings(dlg.settings())
                self.say("🔍 Legacy settings changed (hand edits cleared). Find peaks again to use them.")

    def edit_hydro_settings(self):
        dlg = SettingsDialog("⛓️ Hydrocarbon settings", self.hydrocarbons.settings, parent=self)
        if dlg.exec():
            self.hydrocarbons.settings = dlg.settings()
            self._hydrocarbons_changed()
            self.say("⛓️ Settings changed. Find or propagate the alkane series again to use them.")

    def find_homologs(self, major: int, minor: int, name: str):
        blobs = self.blobs.blobs_for(self.index) if self.index >= 0 else None
        series = self.hydrocarbons.series_for(self.index) if self.index >= 0 else None
        if blobs is None or len(blobs) == 0:
            self.hydro_panel.set_homologs(None, "Find 🔍 blobs on this run first.")
            return
        if series is None or len(series) < 2:
            self.hydro_panel.set_homologs(None, "Find the alkane series on this run first.")
            return
        hs = find_homologous_series(blobs, (major, minor), series, name, self.hydrocarbons.settings)
        if hs is None:
            self.hydro_panel.set_homologs(None, f"No evenly spaced series carrying m/z {major} and {minor} "
                                                "(one member per carbon, consistent t2). Try a lower diagnostic "
                                                "share in the ⚙️ settings.")
        else:
            self.hydro_panel.set_homologs(hs, f"{len(hs.members)} members. Double-click a row to go to it.")

    def edit_stencil(self, name: str):
        series = self.hydrocarbons.series_for(self.index) if self.index >= 0 else None
        stencil = next((s for s in self.hydrocarbons.stencils.stencils if s.name == name), None)
        if series is None or stencil is None:
            self.say("Find the alkane series on this run first")
            return
        self.tabs.setCurrentIndex(0)
        self.view2d.load_stencil_polygon(stencil.to_absolute(series))
        self.hydro_panel.rf_spin.setValue(stencil.response_factor)
        self.hydro_panel.set_drawing(True)
        self.say(f"Drag the corners of {name}, then 💾 Save to replace it")

    def export_figure(self):
        r = self.view2d.result
        if r is None:
            self.say("Nothing to export: open a run first")
            return
        choice, ok = QtWidgets.QInputDialog.getItem(
            self, "Figure width", "Column width:", ["Single column (89 mm)", "Double column (183 mm)"], 0, False)
        if not ok:
            return
        ion = r.expr.text()
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save figure", f"fig_{r.name}_{ion.replace('/', '_over_')}.png",
            "PNG (*.png);;TIFF (*.tif *.tiff)")
        if path:
            self.save_view_figure(Path(path), 89.0 if choice.startswith("Single") else 183.0)

    def save_view_figure(self, path: Path, width_mm: float) -> None:
        """Write the 2D view's current window, levels and scale as a figure at
        300 dpi with a sibling caption `.txt`. Separate from the file dialogs
        so it can be driven directly."""
        r = self.view2d.result
        xr, yr = self.view2d.plot.vb.viewRange()
        lo, hi = self.view2d.image.getLevels()
        fig = plot_chromatogram(self.view2d._disp_after, r.t1 + 0.5 * r.period, r.t2,
                                t1_range_min=tuple(xr), t2_range_s=(max(yr[0], 0.0), yr[1]), levels=(lo, hi),
                                scale=self.view2d.scale, label=r.expr.text(), width_mm=width_mm, displayed=True)
        save_figure(fig, path, dpi=300)
        state = "cleaned" if r.processed else "raw"
        path.with_suffix(".txt").write_text(chromatogram_caption(
            r.name, f"the ion expression {r.expr.text()}", self.view2d.scale, tuple(xr),
            (max(yr[0], 0.0), yr[1]), r.period, r.phase, f"{state} data, untreated 70 eV"), encoding="utf-8")
        self.say(f"Saved {path.name} and its caption")

    def _hydrocarbons_changed(self):
        self.hydro_panel.set_qualifiers(self.hydrocarbons.settings)
        name = self.project.entries[self.index].name if self.index >= 0 else None
        series = self.hydrocarbons.series_for(self.index) if self.index >= 0 else None
        self.hydro_panel.set_current_run(name, series, self.hydrocarbons.running)
        totals = self.hydrocarbons.totals(self.index) if self.index >= 0 else {}
        self.hydro_panel.set_stencils(self.hydrocarbons.stencils.stencils, totals, series is not None)
        markers = self.hydrocarbons.markers(self.index) if self.index >= 0 else {}
        self.hydro_panel.set_markers(markers)
        self._update_stencil_overlay()

    def _update_stencil_overlay(self):
        series = self.hydrocarbons.series_for(self.index) if self.index >= 0 else None
        if series is None or not self.hydrocarbons.stencils.stencils:
            self.view2d.set_stencil_overlay({})
            return
        try:
            self.view2d.set_stencil_overlay(self.hydrocarbons.stencils.lift(series))
        except ValueError:
            self.view2d.set_stencil_overlay({})

    # ================================================================ ✂️
    def show_overlap_dialog(self, candidate_region: int | None = None):
        if getattr(self, "overlap_dialog", None) is None:
            self.overlap_dialog = OverlapDialog(self.project, parent=self)
            self.overlap_dialog.fitAccepted.connect(self.accept_component_group)
            apply_theme(self, self.dark)
        name = self.project.entries[self.index].name if 0 <= self.index < len(self.project.entries) else None
        run = self.project.display_run(self.index) if self.index >= 0 else None
        fold = self.project.fold(self.index) if self.index >= 0 else None
        blobs = self.blobs.blobs_for(self.index) if self.index >= 0 else None
        use_candidates = self.blobs.method == MS_AWARE or candidate_region is not None
        candidates = self.blobs.candidates_for(self.index) if self.index >= 0 and use_candidates else None
        regions = self.blobs.regions_for(self.index) if self.index >= 0 and use_candidates else None
        series = self.hydrocarbons.series_for(self.index) if self.index >= 0 else None
        self.overlap_dialog.set_run(name, run, fold, blobs, series, candidates=candidates,
                                    regions=regions, candidate_region=candidate_region)
        self.overlap_dialog.show()
        self.overlap_dialog.raise_()
        self.overlap_dialog.activateWindow()

    # ================================================================ 🏷️
    def show_identify_dialog(self):
        if getattr(self, "identify_dialog", None) is None:
            self.identify_dialog = IdentifyDialog(self.project, parent=self)
            apply_theme(self, self.dark)
        name = self.project.entries[self.index].name if 0 <= self.index < len(self.project.entries) else None
        run = self.project.display_run(self.index) if self.index >= 0 else None
        fold = self.project.fold(self.index) if self.index >= 0 else None
        blobs = self.blobs.blobs_for(self.index) if self.index >= 0 else None
        components = (self.blobs.components_for(self.index)
                      if self.index >= 0 and self.blobs.method == MS_AWARE else None)
        series = self.hydrocarbons.series_for(self.index) if self.index >= 0 else None
        self.identify_dialog.set_run(name, run, fold, blobs, series, components=components)
        self.identify_dialog.show()
        self.identify_dialog.raise_()
        self.identify_dialog.activateWindow()

    # ================================================================ 🗂️
    def show_batch_dialog(self):
        if getattr(self, "batch_dialog", None) is None:
            default = None
            if self.project.entries:
                p = self.project.entries[0].path
                default = p.parent.parent if p.parent.name == "data" else p.parent
            self.batch_dialog = BatchDialog(self.project, default_folder=default, parent=self)
            self.batch_dialog.navigateRequested.connect(self._batch_navigate)
            apply_theme(self, self.dark)
        self.batch_dialog.show()
        self.batch_dialog.raise_()
        self.batch_dialog.activateWindow()

    def _batch_navigate(self, run_name: str, t1: float, t2: float):
        i = self.project.index_of(run_name)
        if i is None or self.project.entries[i].run is None:
            self.say(f"Open the folder containing {run_name} to jump to it")
            return
        self.set_run(i)
        self.tabs.setCurrentIndex(0)
        self.view2d.centre_on(t1, t2)

    def _layers_changed(self):
        p = self.project
        self.act_processed.setEnabled(p.any_background() or p.any_siloxane() or p.preview is not None)
        self.act_processed.blockSignals(True)
        self.act_processed.setChecked(p.show_processed)
        self.act_processed.blockSignals(False)
        self.act_processed.changed.emit()
        has = self.index >= 0 and p.entries[self.index].background is not None
        self.act_snr.setEnabled(has)
        if not has and self.act_snr.isChecked():
            self.act_snr.setChecked(False)
        self.update_background_info()

    def update_background_info(self):
        if self.index < 0 or self.project.entries[self.index].run is None:
            self.bg_panel.show_summary("Not applied yet")
            self.bg_panel.show_diagnostics(None, None)
            return
        e = self.project.entries[self.index]
        preview = self.project.preview[1] if (self.project.preview and self.project.preview[0] is e) else None
        res = preview or e.background
        if res is None:
            self.bg_panel.show_summary(f"{e.name}: not cleaned. Tick Preview to try the settings, "
                                       "or 🧹 Apply.")
            self.bg_panel.show_diagnostics(None, None)
            return
        tag = "Preview" if preview is not None else "Applied"
        hidden = "" if (preview is not None or self.project.show_processed) else " (hidden: press Space)"
        self.bg_panel.show_summary(f"{tag} on {e.name}{hidden}. {res.summary()}.")
        mz = int(self.bg_panel.diag_ion.value())
        if not (res.run.mz[0] <= mz <= res.run.mz[-1]):
            self.bg_panel.show_diagnostics(None, None, res.top_ions())
            return
        values, sigma = res.residual_pixels(self.project.fold(self.index), mz)
        dense_enough = mz in res.ions_processed
        self.bg_panel.show_diagnostics(values if dense_enough else np.array([]), sigma, res.top_ions())

    # ============================================================ project
    def save_project(self):
        if self.project.path is None:
            self.save_project_as()
            return
        self._store_view_state()
        self.project.save(self.project.path)
        self.undo.setClean()

    def save_project_as(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save project", "project.pygcgc.json",
                                                        "PyGCGC project (*.json)")
        if path:
            self.project.path = Path(path)
            self.save_project()

    def open_project_dialog(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open project", "", "PyGCGC project (*.json)")
        if path:
            self.open_project(path)

    def open_project(self, path, on_loaded=None):
        try:
            self.index = -1
            self.project.load(path, on_loaded=lambda: (self._restore_view_state(), on_loaded and on_loaded()))
        except (ValueError, KeyError, OSError) as exc:
            self.error("Couldn't open the project", str(exc))
            return
        self.undo.clear()

    def _store_view_state(self):
        self.project.view_state = {
            "ion": self.renderer.expr.text(),
            "scale": self.view2d.scale,
            "extended": self.renderer.extended,
            "snr": self.act_snr.isChecked(),
            "current": self.index,
        }

    def _restore_view_state(self):
        v = self.project.view_state or {}
        if "scale" in v:
            self.set_scale(v["scale"])
        if "extended" in v:
            self.act_extend.setChecked(bool(v["extended"]))
        if "ion" in v:
            self.set_ion_expression(v["ion"])
        if isinstance(v.get("current"), int) and 0 <= v["current"] < len(self.project):
            self.set_run(v["current"])
        self.bg_panel.set_run_names(self.project.entries)
        self.bg_panel.set_values(self.project.bg_settings, self.project.bg_blank)
        self._layers_changed()
        if v.get("snr") and self.act_snr.isEnabled():
            self.act_snr.setChecked(True)

    # ============================================================ misc
    def say(self, text: str):
        self.statusBar().showMessage(text, 6000)
        self.log.appendPlainText(text)

    def error(self, title: str, text: str):
        self.log.appendPlainText(f"{title}: {text}")
        QtWidgets.QMessageBox.warning(self, title, text)

    def _busy(self, n: int):
        self.jobs_label.setText("" if n == 0 else f"⏳ {n} job{'s' if n > 1 else ''}")

    def _job_logged(self, job):
        if job.state == "failed":
            self.log.appendPlainText(f"{job.name} failed: {job.message}\n{job.traceback}")
        elif job.state == "done":
            self.log.appendPlainText(f"{job.name}: done")

    def set_dark(self, on: bool):
        self.dark = on
        apply_theme(self, on)

    def show_shortcuts(self):
        rows = "".join(f"<tr><td><b>{k}</b></td><td>&nbsp;&nbsp;{v}</td></tr>" for k, v in SHORTCUTS)
        QtWidgets.QMessageBox.information(self, "Keyboard shortcuts", f"<table>{rows}</table>")

    def show_about(self):
        QtWidgets.QMessageBox.about(self, "PyGCGC",
                                    f"PyGCGC {__version__}\nGC×GC-EI-ToF-MS processing and viewing.")

    def closeEvent(self, event):
        if not self.undo.isClean() and self.project.entries and self.project.path is not None:
            answer = QtWidgets.QMessageBox.question(
                self, "Unsaved changes", "Save the project before closing?",
                QtWidgets.QMessageBox.StandardButton.Save | QtWidgets.QMessageBox.StandardButton.Discard
                | QtWidgets.QMessageBox.StandardButton.Cancel)
            if answer == QtWidgets.QMessageBox.StandardButton.Cancel:
                event.ignore()
                return
            if answer == QtWidgets.QMessageBox.StandardButton.Save:
                self.save_project()
        self.jobs.cancel_all()
        self.jobs.wait(3000)
        super().closeEvent(event)
