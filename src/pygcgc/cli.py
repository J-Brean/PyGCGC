"""Command line for PyGCGC.

    pygcgc info  FILE [FILE ...]
    pygcgc phase FILE [FILE ...] [--period 4] [--targets target_table.csv]
    pygcgc crop  FILE --start 10 --end 14 [--mz-min 35 --mz-max 300] [-o small.cdf]
    pygcgc clean FILE [FILE ...] [--period 4] [--targets CSV] [--outdir cleaned]
    pygcgc quant FOLDER [--outdir results] [--xlsx]
    pygcgc batch FOLDER [--siloxane] [--no-blobs] [--no-quant] [--workers N] [--outdir out]
                 [--report {xlsx,html,pdf,all}]
    pygcgc gui   [FOLDER | PROJECT.json]
    pygcgc fold  FILE [--period 4] [--phase 0.55 | --targets CSV] [--mz 178]
                 [--extend 1.3] [--png out.png] [--npz out.npz]

Without --period, the period is estimated from the run (a few evenly spaced
runs for a sequence) and the estimate is printed; commands stop if none can
be found reliably.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from . import __version__
from .core import (FoldSpec, Fold, auto_phase, check_timing, estimate_period_for_runs, load_anchors,
                   period_from_anchors, phase_for_sequence, read_many, read_run)
from .core.phase import PhaseEstimationError


def _resolve_period(period, runs):
    """The period given, or one estimated from the runs; None if neither."""
    if period is not None:
        return period
    best, found = estimate_period_for_runs(runs)
    if best is None:
        guess = f" (best guess {found[0].period:.4f} s, low confidence)" if found else ""
        print(f"Could not estimate the modulation period from the run{guess}; give --period")
        return None
    print(best.summary())
    return best.period


def _cmd_info(args) -> int:
    for run in read_many(args.files, workers=args.workers):
        rep = check_timing(run.times)
        print(f"== {run.name}")
        print(f"Source: {run.source}")
        if run.meta.get('format') == 'Markes LSC':
            print(f"Native LSC: {run.meta['energy_ev']} eV; {run.meta['processing_state']}")
            print(f"Time reference: {run.meta['time_origin_policy']}")
            print('Validation: vendor spectrum/time comparison pending; the modulation period is not in the file.')
        print(f"Time: {run.times[0]:.2f} to {run.times[-1]:.2f} s  |  m/z {run.mz[0]} to {run.mz[-1]}"
              f"  |  {run.nnz:,} stored points")
        if run.meta.get("comments"):
            print(f"Comment: {run.meta['comments']}")
        print(rep.summary())
        print()
    return 0


def _cmd_phase(args) -> int:
    anchors = load_anchors(args.targets) if args.targets else None
    runs = read_many(args.files, workers=args.workers)
    args.period = _resolve_period(args.period, runs)
    if args.period is None:
        return 2
    estimates = phase_for_sequence(runs, args.period, anchors)
    print(f"{'run':<24}{'phase (s)':>10}  {'method':<9} {'period check'}")
    for run, est in zip(runs, estimates):
        check = ""
        if anchors:
            try:
                p = period_from_anchors(run, anchors, args.period)
                flag = "  <-- check" if p.significant else ""
                check = f"{p.correction * 1e3:+.3f} ± {p.correction_se * 1e3:.3f} ms{flag}"
            except PhaseEstimationError:
                check = "too few anchors"
        print(f"{run.name:<24}{est.phase:>10.3f}  {est.method:<9} {check}")
    return 0


def _cmd_fold(args) -> int:
    run = read_run(args.file)
    args.period = _resolve_period(args.period, [run])
    if args.period is None:
        return 2
    if args.phase is not None:
        phase, how = args.phase, "given"
    else:
        anchors = load_anchors(args.targets) if args.targets else None
        est = auto_phase(run, args.period, anchors)
        phase, how = est.phase, est.summary()
    fold = Fold(run.times, FoldSpec(period=args.period, phase=phase))
    img = fold.ion_image(run, args.mz) if args.mz else fold.tic_image(run)
    t2 = fold.t2
    if args.extend and args.extend > 1:
        img, t2 = fold.extend(img, args.extend)
    print(f"{run.name}: {how}")
    print(f"Image {img.shape[0]} rows × {img.shape[1]} modulations")

    if args.npz:
        np.savez_compressed(args.npz, image=img, t1=fold.t1, t2=t2, period=args.period,
                            phase=phase, mz=np.atleast_1d(args.mz or []))
        print(f"Saved {args.npz}")
    if args.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm

        pos = img[img > 0]
        vmin = np.percentile(pos, 5) if pos.size else 1
        fig, ax = plt.subplots(figsize=(7.0, 3.2), dpi=300)
        ax.imshow(np.clip(img, vmin, None), origin="lower", aspect="auto", cmap="magma",
                  norm=LogNorm(vmin=vmin, vmax=img.max()),
                  extent=[fold.t1[0] / 60, (fold.t1[-1] + args.period) / 60, 0, t2[-1] + fold.row_pitch])
        if args.extend and args.extend > 1:
            # The band above one period repeats the bottom of the next modulation.
            ax.axhspan(args.period, t2[-1] + fold.row_pitch, color="w", alpha=0.12, lw=0)
            ax.axhline(args.period, color="w", lw=0.5, ls="--")
        ax.set_xlabel("$t_1$ (min)")
        ax.set_ylabel("$t_2$ (s)")
        label = "m/z " + "+".join(map(str, args.mz)) if args.mz else "TIC"
        ax.set_title(f"{run.name}, {label}, phase {phase:.3f} s", fontsize=8)
        ax.tick_params(labelsize=7)
        fig.tight_layout()
        fig.savefig(args.png)
        print(f"Saved {args.png}")
    return 0


def _cmd_crop(args) -> int:
    from .core import crop_run, write_andi
    run = read_run(args.file)
    small = crop_run(run, args.start * 60.0, args.end * 60.0, args.mz_min, args.mz_max)
    out = args.output or args.file.with_name(f"{args.file.stem}_crop.cdf")
    write_andi(small, out, comment=f"Cropped from {args.file.name}: {args.start}-{args.end} min")
    size = out.stat().st_size / 1e6
    print(f"Wrote {out} ({small.n_scans:,} scans, m/z {small.mz[0]}-{small.mz[-1]}, {size:.1f} MB)")
    return 0


def _cmd_clean(args) -> int:
    from .core import BackgroundSettings, remove_background, write_andi
    anchors = load_anchors(args.targets) if args.targets else None
    runs = read_many(args.files, workers=args.workers)
    args.period = _resolve_period(args.period, runs)
    if args.period is None:
        return 2
    phases = phase_for_sequence(runs, args.period, anchors) if args.phase is None else None
    args.outdir.mkdir(parents=True, exist_ok=True)
    settings = BackgroundSettings(stiffness=args.stiffness, percentile=args.percentile)
    for k, run in enumerate(runs):
        phase = args.phase if args.phase is not None else phases[k].phase
        fold = Fold(run.times, FoldSpec(period=args.period, phase=phase))
        res = remove_background(run, fold, settings)
        out = write_andi(res.run, args.outdir / f"{run.name}_clean.cdf",
                         comment=f"PyGCGC background removal: {settings.to_dict()}")
        print(f"{run.name}: {res.summary()} -> {out.name}")
    return 0


def find_levels_file(folder: Path) -> Path | None:
    for name in ("calibration_levels.csv", "calibration_levels_ng.csv"):
        if (folder / name).exists():
            return folder / name
    return None


def _cmd_quant(args) -> int:
    import csv
    import json
    from concurrent.futures import ThreadPoolExecutor

    from .core import (BackgroundSettings, QuantSettings, SequenceRun, load_levels, load_targets,
                       quantify_prepared, read_sequence_extras, remove_background)

    folder = args.folder
    period = args.period
    method = json.loads((folder / "method.json").read_text()) if (folder / "method.json").exists() else {}
    if period is None and "modulation_period_s" in method:
        period = float(method["modulation_period_s"])
    targets_path = args.targets or folder / "target_table.csv"
    seq_path = args.sequence or folder / "sequence.csv"
    levels_path = args.levels or find_levels_file(folder)
    for label, pth in (("target table", targets_path), ("sequence", seq_path), ("levels", levels_path)):
        if pth is None or not Path(pth).exists():
            print(f"Missing {label} file")
            return 2
    targets = load_targets(targets_path)
    levels = load_levels(levels_path)
    rows = list(csv.DictReader(open(seq_path, newline="", encoding="utf-8-sig")))
    extras = read_sequence_extras(seq_path)
    files = []
    for row in rows:
        hits = list(folder.rglob(row["file"]))
        if not hits:
            print(f"Can't find {row['file']}")
            return 2
        files.append(hits[0])
    runs = read_many(files, workers=args.workers)
    if period is None:
        method = {**method, "modulation_period_estimated": True}
    period = _resolve_period(period, runs)
    if period is None:
        return 2
    anchors = load_anchors(targets_path)
    phases = phase_for_sequence(runs, period, anchors)

    def prep(k):
        run = runs[k]
        fold = Fold(run.times, FoldSpec(period=period, phase=phases[k].phase))
        res = remove_background(run, fold, BackgroundSettings())
        return run, res.run, fold, res.noise

    with ThreadPoolExecutor(args.workers) as pool:
        prepared = list(pool.map(prep, range(len(runs))))
    info = []
    for run, row in zip(runs, rows):
        ex = extras.get(run.name, {})
        info.append(SequenceRun(run.name, row.get("type", ""), (row.get("level") or "").strip(),
                                ex.get("cal_set", ""), ex.get("volume")))
    settings = QuantSettings(interpolate_cal_sets=not args.no_interpolate)
    report = quantify_prepared(prepared, info, targets, levels, settings,
                               progress=lambda d, t, m: print(f"  {m}"))
    sets = sorted({r.cal_set for r in info if r.cal_set})
    if len(sets) > 1:
        mode = "off" if args.no_interpolate else "on"
        print(f"{len(sets)} calibration sets ({', '.join(sets)}); interpolation {mode}")
    written = report.write(args.outdir)
    n_ok = sum(r.status == "ok" for r in report.results)
    n_check = sum(r.status == "check" for r in report.results)
    fits = [e.get("quant") for e in report.fits.values()]
    good = sum(1 for f in fits if not isinstance(f, str) and f.ok)
    print(f"{len(report.results)} target results: {n_ok} ok, {n_check} to check")
    print(f"{good} of {len(fits)} calibration curves pass")
    for w in written:
        print(f"Wrote {w}")

    from .core.report import build_provenance, write_excel
    prov = build_provenance(__version__, period=period, method=method, settings=settings.to_dict(), runs=runs)
    print(f"Wrote {prov.write(args.outdir / 'provenance.json')}")
    if args.xlsx:
        path = write_excel(args.outdir / "report.xlsx", report=report, provenance=prov)
        print(f"Wrote {path}")
    return 0


def _cmd_batch(args) -> int:
    import json

    from .core import BatchConfig, BatchSettings, pca, run_batch, volcano
    from .core.report import build_provenance, write_excel, write_feature_table_csv, write_html_report, write_pdf_report

    if args.config:
        config = BatchConfig.load(args.config)
    else:
        settings = BatchSettings(run_siloxane=args.siloxane, run_blobs=not args.no_blobs,
                                 run_quant=not args.no_quant)
        config = BatchConfig(period=args.period, settings=settings, reference=args.reference or "")
    if args.save_config:
        config.save(args.save_config)
        print(f"Saved pipeline to {args.save_config}")

    cache_dir = args.cache_dir or (args.folder / ".pygcgc_cache")

    def progress(done, total, message):
        print(f"  [{done}/{total}] {message}")

    result = run_batch(args.folder, config, cache_dir=cache_dir, max_workers=args.workers,
                       progress=progress)
    print(f"{len(result.prepped)} run(s) prepared, {len(result.errors)} failed")
    for name, err in result.errors.items():
        print(f"  {name}: {err}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    if result.report is not None:
        n_ok = sum(r.status == "ok" for r in result.report.results)
        print(f"{len(result.report.results)} target results: {n_ok} ok")
        for w in result.report.write(args.outdir):
            print(f"Wrote {w}")
    run_names = list(result.prepped)
    if result.feature_table:
        path = write_feature_table_csv(result.feature_table, run_names, args.outdir / "feature_table.csv")
        print(f"Wrote {path} ({len(result.feature_table)} features, reference {result.reference})")

    cluster_table = None
    if args.cluster:
        cluster_table = _batch_clusters(result, run_names, config, args)
    if args.regions:
        _batch_regions(result, run_names, args)

    method_path = args.folder / "method.json"
    method = json.loads(method_path.read_text()) if method_path.exists() else {}
    period = result.period
    if result.period_estimated:
        method = {**method, "modulation_period_estimated": True}
    prov = build_provenance(__version__, period=period, method=method, settings=config.settings.to_dict(),
                            runs=[p.run for p in result.prepped.values()])
    print(f"Wrote {prov.write(args.outdir / 'provenance.json')}")

    if args.report != "none":
        qc_values, qc_label = None, ""
        if result.report is not None and result.report.istd:
            first = result.report.istd[0].compound
            rows = [r for r in result.report.istd if r.compound == first]
            qc_values, qc_label = np.array([r.volume for r in rows]), f"{first} (IS)"
        pca_result = pca(result.feature_table, run_names, n_components=2) \
            if result.feature_table and len(run_names) >= 2 else None
        volcano_points = volcano(result.feature_table, list(config.group_a), list(config.group_b)) \
            if result.feature_table and len(config.group_a) >= 2 and len(config.group_b) >= 2 else None
        chart_kwargs = dict(qc_values=qc_values, qc_label=qc_label, pca_result=pca_result,
                            volcano_points=volcano_points)
        if args.report in ("xlsx", "all"):
            path = write_excel(args.outdir / "report.xlsx", report=result.report,
                              feature_table=result.feature_table, run_names=run_names, provenance=prov,
                              cluster_table=cluster_table)
            print(f"Wrote {path}")
        if args.report in ("html", "all"):
            path = write_html_report(args.outdir / "report.html", prov, report=result.report,
                                    feature_table=result.feature_table, run_names=run_names, **chart_kwargs)
            print(f"Wrote {path}")
        if args.report in ("pdf", "all"):
            path = write_pdf_report(args.outdir / "report.pdf", prov, report=result.report,
                                   feature_table=result.feature_table, **chart_kwargs)
            print(f"Wrote {path}")
    return 1 if result.errors and not result.prepped else 0


def _batch_clusters(result, run_names: list, config, args):
    """Pool every run's blobs, cluster them by spectrum, and write the peak and run tables. Returns the cluster
    features (percent of each run's peak volume) for the workbook, or None if there was nothing to cluster."""
    from .core import feature_table_matrix
    from .core.cluster import (ClusterSettings, PeakTable, cluster_features, cluster_peaks, permanova,
                               write_clusters_csv)
    from .core.report import write_feature_table_csv

    tables = [PeakTable.from_blobs(p.blobs, name) for name, p in result.prepped.items() if p.blobs is not None]
    if not any(len(t) for t in tables):
        print("No blobs to cluster (blob detection is off or found nothing)")
        return None
    table = PeakTable.concat(tables)
    clusters = cluster_peaks(table, ClusterSettings(distance_threshold=args.cluster_threshold))
    unassigned = 100.0 * float(table.volume[clusters.labels < 0].sum()) / (float(table.volume.sum()) or 1.0)
    print(f"{len(clusters)} clusters from {len(table)} peaks in {len(tables)} run(s); {unassigned:.1f}% of the volume unassigned")
    print(f"Wrote {write_clusters_csv(table, clusters, args.outdir / 'peak_clusters.csv')}")
    features = cluster_features(table, clusters, run_names, relative=True)
    print(f"Wrote {write_feature_table_csv(features, run_names, args.outdir / 'cluster_run_table.csv')} "
          "(one row per cluster, % of each run's peak volume)")
    if config.group_a and config.group_b and features:
        matrix = feature_table_matrix(features, run_names)
        rows = [run_names.index(n) for n in (*config.group_a, *config.group_b) if n in run_names]
        groups = ["a" if n in config.group_a else "b" for n in (*config.group_a, *config.group_b) if n in run_names]
        try:
            f, p = permanova(matrix[rows], groups)
            print(f"PERMANOVA of the two groups on the cluster shares: pseudo-F = {f:.2f}, p = {p:.3f}")
        except ValueError as exc:
            print(f"PERMANOVA skipped: {exc}")
    return features


def _batch_regions(result, run_names: list, args) -> None:
    """Fit pixel regions on every prepared run at once, write each run's label image (.npy, rows x modulations,
    -1 outside) and the table of each run's share of its signal in each region."""
    from .core import Feature
    from .core.report import write_feature_table_csv
    from .core.segment import SegmentSettings, fit_regions, region_shares

    prepared = [(name, result.prepped[name]) for name in run_names if name in result.prepped]
    model = fit_regions([(p.clean, p.fold) for _, p in prepared], SegmentSettings(n_regions=args.regions))
    labels = {name: model.label_image(p.clean, p.fold) for name, p in prepared}
    shares = {name: region_shares(labels[name], p.clean, p.fold, len(model)) for name, p in prepared}
    first = model.summary(labels[prepared[0][0]], prepared[0][1].clean, prepared[0][1].fold)
    features = [Feature(id=k, ref_t1=first[k].t1_median, ref_t2=first[k].t2_median,
                        volumes={name: float(shares[name][k]) for name, _ in prepared}) for k in range(len(model))]
    for name, lab in labels.items():
        np.save(args.outdir / f"region_labels_{name}.npy", lab)
    path = write_feature_table_csv(features, [n for n, _ in prepared], args.outdir / "region_run_table.csv")
    print(f"{len(model)} pixel regions fitted on {len(prepared)} run(s); wrote {path} "
          "(one row per region, % of each run's signal) and one region_labels_<run>.npy per run")
    for k in range(min(len(model), 8)):
        print(f"  region {k}: top ions {model.top_ions(k, 5)}")


def _cmd_gui(args) -> int:
    from .gui.app import main as gui_main
    return gui_main([str(args.target)] if args.target else [])


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="pygcgc", description="PyGCGC command line")
    p.add_argument("--version", action="version", version=f"PyGCGC {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    pi = sub.add_parser("info", help="summarise files and check scan timing")
    pi.add_argument("files", nargs="+", type=Path)
    pi.add_argument("--workers", type=int, default=4)
    pi.set_defaults(func=_cmd_info)

    pp = sub.add_parser("phase", help="estimate the modulation phase for a sequence")
    pp.add_argument("files", nargs="+", type=Path)
    pp.add_argument("--period", type=float, help="modulation period (s); estimated from the runs if omitted")
    pp.add_argument("--targets", type=Path, help="target table CSV (internal standards used as anchors)")
    pp.add_argument("--workers", type=int, default=4)
    pp.set_defaults(func=_cmd_phase)

    pf = sub.add_parser("fold", help="fold one run into a 2D image")
    pf.add_argument("file", type=Path)
    pf.add_argument("--period", type=float, help="modulation period (s); estimated from the run if omitted")
    pf.add_argument("--phase", type=float, help="phase (s); estimated if omitted")
    pf.add_argument("--targets", type=Path)
    pf.add_argument("--mz", type=int, nargs="*", help="ion(s) to image; TIC if omitted")
    pf.add_argument("--extend", type=float, default=1.0,
                    help="extend the t2 axis to show wrapped peaks whole, e.g. 1.15")
    pf.add_argument("--png", type=Path)
    pf.add_argument("--npz", type=Path)
    pf.set_defaults(func=_cmd_fold)

    pc = sub.add_parser("crop", help="cut a small window out of a run (for test files)")
    pc.add_argument("file", type=Path)
    pc.add_argument("--start", type=float, required=True, help="start (min)")
    pc.add_argument("--end", type=float, required=True, help="end (min)")
    pc.add_argument("--mz-min", type=int)
    pc.add_argument("--mz-max", type=int)
    pc.add_argument("-o", "--output", type=Path)
    pc.set_defaults(func=_cmd_crop)

    pb = sub.add_parser("clean", help="remove background and spikes, write cleaned .cdf files")
    pb.add_argument("files", nargs="+", type=Path)
    pb.add_argument("--period", type=float, help="modulation period (s); estimated from the runs if omitted")
    pb.add_argument("--phase", type=float, help="phase (s); estimated per run if omitted")
    pb.add_argument("--targets", type=Path, help="target table CSV for phase anchors")
    pb.add_argument("--stiffness", type=float, default=1e4)
    pb.add_argument("--percentile", type=float, default=20.0)
    pb.add_argument("--outdir", type=Path, default=Path("cleaned"))
    pb.add_argument("--workers", type=int, default=4)
    pb.set_defaults(func=_cmd_clean)

    pq = sub.add_parser("quant", help="calibrate and quantify a folder (sequence, targets, levels)")
    pq.add_argument("folder", type=Path)
    pq.add_argument("--period", type=float, help="modulation period (s); else method.json, else estimated")
    pq.add_argument("--targets", type=Path)
    pq.add_argument("--sequence", type=Path)
    pq.add_argument("--levels", type=Path, help="calibration levels CSV")
    pq.add_argument("--no-interpolate", action="store_true",
                    help="with more than one calibration set, use the nearest one before a sample "
                         "instead of blending the sets either side of it")
    pq.add_argument("--outdir", type=Path, default=Path("pygcgc_results"))
    pq.add_argument("--workers", type=int, default=4)
    pq.add_argument("--xlsx", action="store_true",
                    help="also write report.xlsx, flagged rows highlighted (needs the report extra)")
    pq.set_defaults(func=_cmd_quant)

    pba = sub.add_parser("batch", help="cached, parallel pipeline over a whole folder; feature table")
    pba.add_argument("folder", type=Path)
    pba.add_argument("--period", type=float, help="modulation period (s); else method.json")
    pba.add_argument("--cache-dir", type=Path, help="step cache; default <folder>/.pygcgc_cache")
    pba.add_argument("--workers", type=int, help="process pool size; default cores - 1")
    pba.add_argument("--outdir", type=Path, default=Path("pygcgc_batch"))
    pba.add_argument("--siloxane", action="store_true", help="also run 🫧 siloxane removal")
    pba.add_argument("--no-blobs", action="store_true", help="skip 🔍 blob detection (and the feature table)")
    pba.add_argument("--no-quant", action="store_true", help="skip 🧪 quantification")
    pba.add_argument("--reference", help="run to build the feature table against; default: "
                                         "the first run with blobs")
    pba.add_argument("--cluster", action="store_true",
                     help="also cluster the peaks of every run by mass spectrum (needs 🔍): writes peak_clusters.csv "
                          "and cluster_run_table.csv, and tests the config's two run groups by PERMANOVA if it has any")
    pba.add_argument("--cluster-threshold", type=float, default=0.3,
                     help="cosine distance at which clusters stop merging (default 0.3; lower gives more clusters)")
    pba.add_argument("--regions", type=int, metavar="N",
                     help="also divide every run into N pixel regions of similar mass spectrum, with no peak "
                          "detection (needs scikit-learn): writes region_run_table.csv and a region_labels_<run>.npy "
                          "per run")
    pba.add_argument("--config", type=Path, help="a saved pipeline JSON; overrides the flags above")
    pba.add_argument("--save-config", type=Path, help="save the settings used as a pipeline JSON")
    pba.add_argument("--report", choices=["none", "xlsx", "html", "pdf", "all"], default="none",
                     help="also write a formatted spreadsheet and/or a document report "
                          "(needs the report/plot extras); provenance.json is always written")
    pba.set_defaults(func=_cmd_batch)

    pg_ = sub.add_parser("gui", help="open the desktop app")
    pg_.add_argument("target", nargs="?", type=Path, help="folder of runs or a project file")
    pg_.set_defaults(func=_cmd_gui)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
