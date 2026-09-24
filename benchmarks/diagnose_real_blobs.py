"""Reproducible local blob audit; generated real-data artifacts stay in real_data/.

Run with the project's Python environment. No production settings are changed.
The cache preserves absolute scan times and shifts the fold origin by whole
modulations only, so the supplied phase keeps its original meaning.
"""
from pathlib import Path
import json
import time
from dataclasses import replace

import numpy as np
from scipy import sparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pygcgc.core.io.markes_lsc import read_lsc
from pygcgc.core.run import Run
from pygcgc.core.fold import fold
from pygcgc.core.background import remove_background, storage_threshold
from pygcgc.core import peaks as p


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[1]
    out = root / "real_data" / "blob_diagnosis"
    out.mkdir(exist_ok=True)
    cache = out / "segment.npz"
    if not cache.exists():
        source = next((root / "real_data").glob("*_70eV.lsc"))
        print("Reading", source.name, flush=True)
        run = read_lsc(source)
        origin = float(run.times[0])
        # Ten-minute context for the 51-modulation background pattern filter.
        use = (run.times >= origin + 600) & (run.times < origin + 1200)
        data = run.data[use]
        np.savez_compressed(cache, times=run.times[use], mz=run.mz,
                            values=data.data, indices=data.indices, indptr=data.indptr,
                            shape=data.shape)
        (out / "segment.json").write_text(json.dumps(dict(
            source=str(source), sha256=run.meta["source_sha256"],
            original_t0=origin, fold_t0=origin + 600, period=10., phase=1.81,
            time_window=[origin + 600, origin + 1200]), indent=2))
    z = np.load(cache)
    meta = json.loads((out / "segment.json").read_text())
    run = Run(z["times"], z["mz"], sparse.csr_matrix(
        (z["values"], z["indices"], z["indptr"]), shape=tuple(z["shape"])), "audit segment")
    f = fold(run, 10., 1.81, t0=meta["fold_t0"])
    print("Segment", run.n_scans, "scans;", f.n_rows, "rows;", f.n_cols, "columns", flush=True)
    bg = remove_background(run, f)
    settings = p.PeakSettings()
    frame = p.make_frame(f, settings)
    xx, yy = np.meshgrid(frame.mid_t1 / 60., np.arange(frame.rows) * frame.pitch)
    focus = (xx >= 14) & (xx <= 16) & (yy >= .6) & (yy <= 2.)
    quiet = (xx >= 14) & (xx <= 16) & (yy >= 5) & (yy <= 9)
    results = dict(segment=meta, shape=[f.n_rows, f.n_cols], scans=run.n_scans, cases={})
    # A small, all-ion cube for subsequent decomposition experiments.
    rr = np.flatnonzero((np.arange(f.n_rows)*f.row_pitch >= .6) & (np.arange(f.n_rows)*f.row_pitch <= 2.))
    cc = np.flatnonzero((frame.mid_t1 >= 840) & (frame.mid_t1 <= 960))
    np.savez_compressed(out / "focus_cube.npz", raw=f.cube(run)[rr][:, cc],
                        clean=f.cube(bg.run)[rr][:, cc], mz=run.mz,
                        t1=frame.mid_t1[cc], t2=rr*f.row_pitch)
    rawfloor = storage_threshold(run)/3.
    cases = [("raw_default", run, settings, None), ("clean_default", bg.run, settings, None),
             ("clean_raw_floor", bg.run, settings, rawfloor),
             ("clean_tic_only", bg.run, replace(settings, ion_snr=False), None),
             ("clean_rel_floor_10pct", bg.run, replace(settings, blob_rel_floor=.1), None),
             ("clean_less_smoothing", bg.run, replace(settings, smooth_modulations=.25, smooth_t2_s=.015), None)]
    plt.rcParams.update({"font.size": 8})
    fig, axes = plt.subplots(len(cases), 2, figsize=(8, 15), constrained_layout=True)
    vmax = float(np.quantile(np.log10(1+np.maximum(frame.extend(f.tic_image(run))[focus | quiet],0)), .999))
    for axrow, (name, work, s, override) in zip(axes, cases):
        start = time.perf_counter()
        raw = frame.extend(f.tic_image(work))
        sm = p._smooth(raw, frame, s)
        floor = storage_threshold(work) / 3. if override is None else override
        tic = sm / p._column_noise(sm, frame, s.noise_block, floor)[None, :]
        ions, cube = p._active_ion_cube(work, f, frame, s)
        smcube, isum = p._ion_snr_cube(cube, frame, s, floor)
        det = np.maximum(tic, isum) if s.ion_snr else tic
        blobs = p.find_blobs(work, f, s, floor=override)
        records = []
        for i in np.unique(blobs.labels[focus]):
            if not i:
                continue
            r, c = np.nonzero(blobs.labels == i)
            records.append(dict(id=int(i), t1=float(blobs.t1[i-1]), t2=float(blobs.t2[i-1]),
                n_pixels=int(blobs.n_pixels[i-1]), n_modulations=int(blobs.n_modulations[i-1]),
                t1_span_s=float((c.max()-c.min()+1)*10), t2_span_s=float((r.max()-r.min()+1)*frame.pitch),
                coelution=bool(blobs.co_elution[i-1])))
        thresholds = np.r_[np.inf, s.blob_rel_floor*blobs.height]
        labelled = blobs.labels > 0
        violations = labelled & (sm < thresholds[blobs.labels])
        info = dict(seconds=time.perf_counter()-start, blobs=len(blobs), floor=floor,
            final_relative_floor_violations=int(violations.sum()),
            active_ions=len(ions), focus_blobs=records,
            focus_mask_fraction=float(np.mean(det[focus] >= s.floor_snr)),
            quiet_mask_fraction=float(np.mean(det[quiet] >= s.floor_snr)),
            quiet_label_fraction=float(np.mean(blobs.labels[quiet] > 0)),
            quiet_score_quantiles=np.quantile(det[quiet], [.1,.5,.9]).tolist(),
            quiet_tic_quantiles=np.quantile(raw[quiet], [.1,.5,.9]).tolist(),
            focus_score_quantiles=np.quantile(det[focus], [.1,.5,.9]).tolist())
        results["cases"][name] = info
        print(name, json.dumps({k:v for k,v in info.items() if k != "focus_blobs"}), flush=True)
        for ax, window in zip(axrow, [(14,16,.6,2.), (14,16,5,9)]):
            mesh = ax.pcolormesh(xx, yy, np.log10(1+np.maximum(raw,0)), cmap="magma", shading="nearest", rasterized=True, vmin=0, vmax=vmax)
            for i in range(len(blobs)):
                tx, ty = blobs.outline_times(i)
                ax.plot(tx/60., ty, color="cyan", lw=.5)
            ax.set(xlim=window[:2], ylim=window[2:], title=name, xlabel="t1 (min)", ylabel="t2 (s)")
    fig.colorbar(mesh, ax=axes, label="log10(1 + TIC); common scale", shrink=.4)
    fig.savefig(out / "comparison.png", dpi=200)
    plt.close(fig)
    # Independent ion evidence, without claiming these traces are pure components.
    fig, axes = plt.subplots(4, 2, figsize=(8, 9), constrained_layout=True)
    ionlist = [43, 57, 71, 85, 91, 105, 128, 142]
    for ax, ion in zip(axes.flat, ionlist):
        img = frame.extend(f.image(bg.run.xic(ion)))
        mesh = ax.pcolormesh(xx, yy, np.log10(1+np.maximum(img,0)), shading="nearest", cmap="magma")
        fig.colorbar(mesh, ax=ax, label="log10(1 + ion signal)")
        ax.set(xlim=(14,16), ylim=(.6,2.), title=f"m/z {ion}", xlabel="t1 (min)", ylabel="t2 (s)")
    fig.savefig(out / "ions.png", dpi=150)
    plt.close(fig)
    (out / "metrics.json").write_text(json.dumps(results, indent=2))
    print("Outputs:", out, flush=True)
