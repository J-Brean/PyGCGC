"""Regenerate the screenshots in docs/images from the synthetic demo data.

    python examples/demo_data/generate_demo_data.py demo_gcxgc
    python docs/images/make_screenshots.py demo_gcxgc [STEP ...]

Drives the real app on the real display (an offscreen platform has no fonts), so
windows appear and disappear while it runs. Every step is what a user would do,
called through the same methods the buttons call. The results cache goes to a
temporary folder and any project file left in the data folder is removed first.
Name steps to make only those screenshots (with everything they depend on).
"""

import csv
import os
import re
import sys
import tempfile
import time
from pathlib import Path

os.environ.setdefault("PYGCGC_CACHE", tempfile.mkdtemp(prefix="pygcgc-shots-cache-"))

from PySide6 import QtCore, QtWidgets

from pygcgc.gui.mainwindow import MainWindow



def main():
    # A worker process of a process pool re-imports this file on Windows; it must not run it.
    DATA = Path(sys.argv[1] if len(sys.argv) > 1 else "demo_gcxgc").resolve()
    ONLY = set(sys.argv[2:])
    OUT = Path(__file__).resolve().parent
    for stale in DATA.glob("*.pygcgc.json"):
        stale.unlink()

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    app.setApplicationName("PyGCGC")
    w = MainWindow()
    w.resize(1680, 1000)
    w.show()


    def pump(seconds=0.2):
        end = time.time() + seconds
        while time.time() < end:
            app.processEvents()
            time.sleep(0.01)


    def wait(cond, timeout=900, what="condition"):
        start = time.time()
        while not cond():
            if time.time() - start > timeout:
                raise TimeoutError(f"timed out waiting for {what}")
            app.processEvents()
            time.sleep(0.02)
        pump(0.3)


    def idle(jobs=None, timeout=1800):
        jobs = jobs or w.jobs
        wait(lambda: jobs.active_count() == 0, timeout, "jobs")
        pump(0.6)


    def shot(name, widget=None):
        if ONLY and name not in ONLY:
            return
        widget = widget or w
        pump(0.8)
        pix = widget.grab()
        if pix.width() > 2000:      # a high-DPI screen grabs at twice the size; keep the files small
            pix = pix.scaledToWidth(2000, QtCore.Qt.TransformationMode.SmoothTransformation)
        pix.save(str(OUT / f"{name}.png"))
        print("saved", name, flush=True)


    def run_index(name):
        return next(i for i, e in enumerate(w.project.entries) if e.name == name)


    def show(run, ion="TIC", t1_min=None, t2_s=None, tab=0):
        """Open `run` in the 2D view on `ion`, zoomed to t1 (minutes) and t2 (seconds) ranges."""
        w.tabs.setCurrentIndex(tab)
        w.set_run(run_index(run))
        w.set_ion_expression(ion)
        idle()
        wait(lambda: w.view2d.result is not None and w.view2d.result.name == run, what=f"{run} to render")
        if t1_min is None:
            w.view2d.reset_view()
        else:
            w.view2d.plot.setRange(xRange=t1_min, yRange=t2_s, padding=0)
        w.view2d.auto_levels()
        pump(0.5)


    def spectrum_of(run, compound):
        """Show the spectrum at a compound's apex, as a click on it would."""
        t1, t2 = truth[(run, compound)]
        fold = w.project.fold(run_index(run))
        col = int((t1 - fold.t1[0]) // fold.spec.period)
        w.show_spectrum_at(float(fold.t1[col] + t2))
        idle()


    # Where each compound really is in each run (the generator's truth), for the clicks a user would make.
    truth = {(r["run"], r["compound"]): (float(r["t1_s"]), float(r["t2_s"]))
             for r in csv.DictReader(open(DATA / "truth" / "truth_amounts.csv", newline=""))}
    ISTDS = {"Acenaphthene-d10": ("acenaphthene-d10", 164), "Phenanthrene-d10": ("phenanthrene-d10", 188),
             "n-Tetracosane-d50": ("n-tetracosane-d50", 66), "Chrysene-d12": ("chrysene-d12", 240),
             "Perylene-d12": ("perylene-d12", 264)}
    CAL_RUNS = ["CAL_L1_1ng", "CAL_L2_5ng", "CAL_L3_20ng", "CAL_L4_50ng"]
    URBAN, DIESEL, SMOKE = "SAMPLE_1_urban_air", "SAMPLE_2_diesel_exhaust", "SAMPLE_3_wood_smoke"
    w.jobs_panel.clear_after_ms = 10

    # ------------------------------------------------------------ 1 open and fold
    loaded = []
    w.open_folder(DATA, on_loaded=lambda: loaded.append(True))
    wait(lambda: loaded and w.view2d.result is not None, what="runs to load")
    idle()
    wait(lambda: all(e.phase_method != "unset" for e in w.project.entries), what="phases")
    w.workflow.show_step("fold")
    show(DIESEL)
    shot("screenshot_fold")
    w.left_tabs.setCurrentWidget(w.runs_panel)
    shot("screenshot_runs")

    # ----------------------------------------------------------------- 2 align
    # The reference is a mid-sequence standard; five internal standards and three alkanes are
    # clicked on the urban sample, then found on every other run.
    w.workflow.show_step("align")
    w.set_reference("CAL_L3_20ng")
    # Perylene-d12 wraps past the end of the modulation, so it is not clicked.
    anchors = [(c, ion) for c, ion in ISTDS.values() if c != "perylene-d12"] + [("n-C14", 57), ("n-C20", 57), ("n-C26", 57)]
    for k, (compound, ion) in enumerate(anchors):
        (r1, r2), (s1, s2) = truth[("CAL_L3_20ng", compound)], truth[(URBAN, compound)]
        w.align.add_anchor(URBAN, w.align.next_anchor_name(), r1, r2, s1, s2, ion)
    w.propagate_anchors()
    idle()
    show(URBAN)
    spectrum_of(URBAN, "phenanthrene-d10")
    shot("screenshot_align")
    w.tabs.setCurrentWidget(w.compare)
    w.set_run(run_index(URBAN))
    w.compare.mode.setCurrentIndex(w.compare.mode.findData("overlay"))
    idle()
    shot("screenshot_align_overlay")

    # ------------------------------------------------------------- 3 background
    w.workflow.show_step("background")
    w.apply_background()
    idle()
    show(SMOKE, "73", (32.6, 36.4), (0.0, 2.6))
    w.bg_panel.split.setChecked(True)
    idle()
    shot("screenshot_background_split")
    w.bg_panel.split.setChecked(False)
    show(DIESEL)
    w.act_snr.setChecked(True)
    spectrum_of(DIESEL, "pristane")
    shot("screenshot_snr")
    w.act_snr.setChecked(False)
    show(DIESEL)
    shot("screenshot_main")

    # ------------------------------------------------------------------ 4 peaks
    # MS-aware separation first (the default), then legacy blobs, which the siloxane filterer,
    # hydrocarbon assignment and batch still read.
    from pygcgc.gui.blob_state import LEGACY, MS_AWARE  # noqa: E402

    w.workflow.show_step("blobs")
    w.blobs.set_method(MS_AWARE)
    w.run_blobs()
    idle()
    w.act_blob_overlay.setChecked(True)
    show(URBAN, "TIC", (17.0, 27.0), (0.6, 3.0))
    spectrum_of(URBAN, "fluorene")
    shot("screenshot_peaks")
    w.blobs.set_method(LEGACY)
    w.run_blobs()
    idle()
    show(URBAN, "TIC", (17.0, 27.0), (0.6, 3.0))
    shot("screenshot_peaks_legacy")

    # ------------------------------------------------------------ 4a siloxane filterer
    w.workflow.show_step("siloxane")
    w.sil_panel.threshold.setValue(15)
    show(URBAN, "73+281+355", (6.0, 22.0), (0.2, 2.2))
    spectrum_of(URBAN, "D5")
    shot("screenshot_siloxane")

    # ----------------------------------------------------------- 4b hydrocarbons
    w.workflow.show_step("hydro")
    i = run_index(DIESEL)
    show(DIESEL, "57")
    seeds = [(n, *truth[(DIESEL, f"n-C{n}")]) for n in (14, 20)]
    w.hydrocarbons.find_series(i, seeds)
    idle()
    w.propagate_alkanes()
    idle()
    w.assign_hydrocarbons()
    idle()
    w.act_blob_overlay.setChecked(False)
    show(DIESEL, "57", (6.0, 42.0), (0.4, 3.0))
    spectrum_of(DIESEL, "n-C17")
    shot("screenshot_hydrocarbons")

    # --------------------------------------------------------- 5 internal standards
    # From here on the MS-aware peaks are used again. The standards are labelled on the lowest
    # standard, where the native compounds next to them are smallest.
    w.blobs.set_method(MS_AWARE)
    idle()
    w.workflow.show_step("istd")
    i = run_index("CAL_L1_1ng")
    show("CAL_L1_1ng")
    for name, (compound, _ion) in ISTDS.items():
        t1, t2 = truth[("CAL_L1_1ng", compound)]
        if t2 >= w.project.period:          # a wrapped standard is labelled where it shows, one modulation on
            t1, t2 = t1 + w.project.period, t2 - w.project.period
        ref = w._nearest_peak_ref(i, t1, t2)
        if ref is not None:
            w.istd.label(name, i, ref)
            idle()
    w.propagate_istd()
    idle()
    spectrum_of("CAL_L1_1ng", "chrysene-d12")
    shot("screenshot_istd")

    # --------------------------------------------------------------- 6 quantify
    w.workflow.show_step("quant")
    w.run_quant()
    idle()
    wait(lambda: not w.quant.running, what="quantification")
    idle()
    show(URBAN, "TIC", (30.0, 38.0), (2.5, 6.5))
    shot("screenshot_quant")
    w.tabs.setCurrentWidget(w.review)
    w.review.show_compound("benz[a]anthracene")
    shot("screenshot_review")
    w.tabs.setCurrentWidget(w.calibration_panel)
    w.calibration_panel.select_compound("fluorene")
    shot("screenshot_calibration")
    w.left_tabs.setCurrentIndex(w.left_tabs.indexOf(w.results.parentWidget().parentWidget()))
    w.tabs.setCurrentIndex(0)
    w.resizeDocks([w.workspace_dock], [1000], QtCore.Qt.Orientation.Horizontal)
    shot("screenshot_results")
    w.resizeDocks([w.workspace_dock, w.spectrum_dock], [470, 380], QtCore.Qt.Orientation.Horizontal)
    w.act_overlay.setChecked(False)

    # ---------------------------------------------------------------- 7 classify
    w.workflow.show_step("classify")
    show(SMOKE)
    w.start_classify()
    idle()
    w.class_panel.overlay_box.setChecked(True)
    w.act_blob_overlay.setChecked(True)          # class colours are drawn on the peak outlines
    show(SMOKE, "TIC", (5.0, 40.0), (0.4, 5.5))
    spectrum_of(SMOKE, "syringol")
    shot("screenshot_classify")
    w.class_panel.overlay_box.setChecked(False)
    w.act_blob_overlay.setChecked(False)

    # ---------------------------------------------------------------- 8 identify
    w.workflow.show_step("identify")
    show(URBAN)
    w.show_identify_dialog()
    d = w.identify_dialog
    d.resize(1300, 760)
    d.load_library_file(DATA / "demo_library.msp")
    idle()
    d._run_identify()
    idle()
    wait(lambda: d.table.rowCount() > 0, what="identification")
    name_col = next(c for c in range(d.table.columnCount()) if d.table.horizontalHeaderItem(c).text() == "name")
    for r in range(d.table.rowCount()):          # unmatched peaks are listed too; show the hits
        item = d.table.item(r, name_col)
        d.table.setRowHidden(r, item is None or not item.text().strip())
    # Scroll past the toluene smear at the start, as a user would.
    t1_col = next(c for c in range(d.table.columnCount()) if d.table.horizontalHeaderItem(c).text() == "t1 (min)")
    first = next(r for r in range(d.table.rowCount()) if not d.table.isRowHidden(r) and float(d.table.item(r, t1_col).text()) > 8)
    d.table.scrollToItem(d.table.item(first, 0), QtWidgets.QAbstractItemView.ScrollHint.PositionAtTop)
    shot("screenshot_identify", d)
    d.close()

    # ----------------------------------------------------------------- 9 cluster
    w.workflow.show_step("cluster")
    show(DIESEL)
    w.start_cluster()
    idle()
    w.cluster_panel.fill_box.setChecked(True)
    show(DIESEL, "TIC", (5.0, 42.0), (0.4, 4.5))
    shot("screenshot_cluster")
    w.cluster_panel.fill_box.setChecked(False)

    # ------------------------------------------------------------ 1D, overlap, batch
    show(DIESEL, "TIC", None, None, tab=1)
    shot("screenshot_1d")
    w.tabs.setCurrentIndex(0)
    w.blobs.set_method(MS_AWARE)
    idle()
    show(URBAN, "TIC", (17.0, 27.0), (0.6, 3.0))
    w.show_overlap_dialog()
    od = w.overlap_dialog
    od.resize(1300, 820)
    idle()
    # The region holding benz[a]anthracene and chrysene, which overlap.
    t1c, t2c = truth[(URBAN, "chrysene")]
    def where(k):
        m = re.search(r"t1=([\d.]+) min\s+t2=([\d.]+) s", od.blob_combo.itemText(k))
        return abs(float(m[1]) * 60 - t1c) / 30 + abs(float(m[2]) - t2c) / 0.3
    od.blob_combo.setCurrentIndex(min(range(od.blob_combo.count()), key=where))
    idle()
    shot("screenshot_overlap", od)
    od.tabs.setCurrentIndex(3)
    od._run_fit()
    idle()
    shot("screenshot_overlap_fit", od)
    od.close()

    w.show_batch_dialog()
    b = w.batch_dialog
    b.resize(1300, 820)
    b.folder_edit.setText(str(DATA))
    b._run_clicked()
    idle(timeout=3600)
    wait(lambda: b.result is not None, timeout=3600, what="batch")
    b.tabs.setCurrentIndex(2)
    shot("screenshot_batch", b)
    b.close()

    # ------------------------------------------------------------ 📈 calibrate window
    w.show_calibrate_window()
    cw = w.calibrate_window
    cw.resize(1680, 1000)
    cw.add_main_entries([e for e in w.project.entries if e.name in CAL_RUNS])
    pump(2)
    wait(lambda: not cw.processing and cw.jobs.active_count() == 0, timeout=3600, what="standards")
    idle(cw.jobs)
    cw.tabs.setCurrentIndex(0)
    shot("screenshot_calibrate_standards", cw)
    cw.load_compounds(DATA / "calibration_compounds.csv")
    idle(cw.jobs)
    cw.tabs.setCurrentIndex(1)
    cw.select_compound("fluorene")
    idle(cw.jobs)
    shot("screenshot_calibrate_compounds", cw)
    cw.tabs.setCurrentIndex(2)
    cw.curves.select("fluorene")
    idle(cw.jobs)
    shot("screenshot_calibrate_curves", cw)
    cw.tabs.setCurrentIndex(3)
    cw.quantify_samples([URBAN, DIESEL, SMOKE])
    idle(cw.jobs)
    wait(lambda: cw.jobs.active_count() == 0 and cw.sample_results, timeout=3600, what="samples")
    idle(cw.jobs)
    shot("screenshot_calibrate_apply", cw)
    cw.close()

    # -------------------------------------------------------------- dark theme
    w.act_theme.setChecked(True)
    w.workflow.show_step("blobs")
    show(DIESEL, "TIC", (5.0, 42.0), (0.4, 4.5))
    spectrum_of(DIESEL, "phenanthrene")
    shot("screenshot_dark")
    w.undo.clear()
    w.jobs.cancel_all()
    print("finished", flush=True)
    os._exit(0)      # skip Qt's teardown, which races the autosave timer


if __name__ == "__main__":
    main()
