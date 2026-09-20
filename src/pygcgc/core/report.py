"""📤 Reports: provenance, spreadsheet and document exports, and figures.

Every output this module writes says, on its own, how it was made: which
version of PyGCGC, what modulation period and settings were used, and what
a person changed by hand along the way (`Provenance`, built by
`build_provenance`). The manual-edit log is optional and empty by default,
since a scripted `pygcgc batch` run or the batch dialog has no interactive
review to report -- 🗂️ is deliberately decoupled from the per-run Project
state that accumulates 🔍's edit log, 🧪's decisions and ⛓️'s stencils (see
docs/architecture.md section 5); only a caller that has that state can pass it in.

The rest of the module turns a `QuantReport` and/or a feature table into
actual files: a formatted spreadsheet with flagged rows picked out, a
self-contained HTML page, a short print-ready PDF, and individual figures
sized for a single-column (89 mm) or double-column (183 mm) publication
figure. Every chart is matplotlib (`plot` extra), rendered headless
("Agg"); every spreadsheet is openpyxl (`report` extra). Both degrade with
a plain `ImportError` if the extra isn't installed, the same way `pygcgc
fold --png` already does -- nothing here tries to hide that dependency.
"""

from __future__ import annotations

import base64
import csv
import io
import html
import json
import textwrap
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from .batch import PCAResult, VolcanoPoint, westgard_flags
from .quant import QuantReport, calibration_table


# ------------------------------------------------------------- provenance
@dataclass
class Provenance:
    """What made an output: version, method, settings, and anything a
    person changed by hand. Saved as its own small JSON file, and folded
    into every richer report (Excel/HTML/PDF) as well."""
    version: str
    generated: str
    period: float | None = None
    method: dict = field(default_factory=dict)
    settings: dict = field(default_factory=dict)
    edits: dict = field(default_factory=dict)
    sources: list = field(default_factory=list)

    def to_dict(self) -> dict:
        return asdict(self)

    def write(self, path) -> Path:
        path = Path(path)
        path.write_text(json.dumps(self.to_dict(), indent=2, default=str), encoding="utf-8")
        return path


def build_provenance(version: str, *, period: float | None = None, method: dict | None = None,
                     settings: dict | None = None, blob_edits: dict | None = None,
                     quant_decisions=None, stencils=None, runs=()) -> Provenance:
    """Gather what an output needs to say about how it was made.

    `blob_edits` is `BlobState.edits` (per run name, a list of 🔍's
    delete/merge/add operations); `quant_decisions` is a `QuantDecisions`
    (🧪's accept/reject/exclude/override/pick-by-hand record); `stencils`
    is a `StencilSet` (⛓️'s drawn frame). All three are optional and
    contribute nothing when omitted, which is the normal case for a batch
    run: only an interactive Project accumulates them.
    """
    edits: dict = {}
    if blob_edits:
        by_run = {run: len(ops) for run, ops in blob_edits.items() if ops}
        if by_run:
            edits["blob_edits"] = {"total": sum(by_run.values()), "by_run": by_run}
    if quant_decisions is not None:
        d = quant_decisions.to_dict()
        counts = {"accepted": len(d["accepted"]), "rejected": len(d["rejected"]),
                 "picked_by_hand": len(d["assigned"]),
                 "integrated_by_hand": len(d.get("footprints", [])),
                 "excluded_from_curve": sum(len(v) for v in d["excluded"].values()),
                 "overridden_curves": len(d["overrides"])}
        if any(counts.values()):
            edits["quant_decisions"] = counts
    if stencils is not None and getattr(stencils, "stencils", None):
        edits["stencils"] = {"frame": stencils.name, "names": [s.name for s in stencils.stencils]}
    sources = [dict(name=run.name, path=str(run.source), **{k:run.meta[k] for k in (
        'format', 'source_sha256', 'parser_version', 'acquisition_id', 'energy_ev',
        'processing_state', 'time_origin_seconds', 'time_origin_policy', 'vendor_validated') if k in run.meta})
        for run in runs]
    return Provenance(version=version, generated=datetime.now(timezone.utc).isoformat(timespec="seconds"),
                      period=period, method=dict(method or {}), settings=dict(settings or {}), edits=edits,
                      sources=sources)


# ------------------------------------------------------------- plain CSV
def write_feature_table_csv(feature_table: list, run_names: list, path) -> Path:
    """The feature table as CSV, with one flag column: which runs' value
    for this feature was gap-filled (`batch._box_sum`) rather than matched
    to an actual blob, so a spreadsheet reader can tell "not detected,
    estimated as near zero" from "matched a real blob" without guessing
    from the number alone."""
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["feature", "t1_min", "t2_s"] + list(run_names) + ["gap_filled"])
        for f in feature_table:
            w.writerow([f.id, f.ref_t1 / 60.0, f.ref_t2]
                      + [f.volumes.get(n, "") for n in run_names]
                      + ["; ".join(n for n in run_names if f.filled.get(n))])
    return path


def _clean(v):
    if isinstance(v, float) and not np.isfinite(v):
        return None
    if isinstance(v, (list, tuple)):
        return "; ".join(map(str, v))
    return v


# ------------------------------------------------------------- Excel
def write_excel(path, *, report: QuantReport | None = None, feature_table: list | None = None,
                run_names: list | None = None, provenance: Provenance | None = None) -> Path:
    """One workbook, one sheet per table, a light fill on any row that
    needs a look: a target result flagged or not "ok", a calibration that
    didn't pass, an internal standard out of tolerance, or a feature whose
    value for some run was gap-filled rather than measured. Needs the
    `report` extra (openpyxl)."""
    from openpyxl import Workbook
    from openpyxl.styles import Font, PatternFill

    flag_fill = PatternFill(fill_type="solid", start_color="FDE2C8", end_color="FDE2C8")
    filled_fill = PatternFill(fill_type="solid", start_color="DCEBF7", end_color="DCEBF7")
    header_font = Font(bold=True)

    path = Path(path)
    wb = Workbook()
    wb.remove(wb.active)

    def sheet(title, rows, flagged=None, fill=flag_fill):
        ws = wb.create_sheet(title[:31])
        if not rows:
            ws.append(["(none)"])
            return ws
        fields = list(rows[0].keys())
        ws.append(fields)
        for cell in ws[1]:
            cell.font = header_font
        for i, row in enumerate(rows, start=2):
            ws.append([_clean(row[k]) for k in fields])
            if flagged and flagged(row):
                for cell in ws[i]:
                    cell.fill = fill
        for col in ws.columns:
            width = max((len(str(c.value)) for c in col if c.value is not None), default=8)
            ws.column_dimensions[col[0].column_letter].width = min(40, max(8, width + 2))
        return ws

    if report is not None:
        sheet("Targets", [asdict(r) for r in report.results],
             flagged=lambda r: r["status"] in ("check", "rejected") or bool(r["flags"]))
        sheet("Calibrations", calibration_table(report.fits), flagged=lambda r: not r["passes"])
        sheet("All signals", [asdict(r) for r in report.signals])
        sheet("Internal standards", [asdict(r) for r in report.istd], flagged=lambda r: r["flag"])

    if feature_table:
        names = run_names or sorted({n for f in feature_table for n in f.volumes})
        rows = []
        for f in feature_table:
            row = {"feature": f.id, "t1_min": f.ref_t1 / 60.0, "t2_s": f.ref_t2}
            row.update({n: f.volumes.get(n) for n in names})
            row["gap_filled"] = "; ".join(n for n in names if f.filled.get(n))
            rows.append(row)
        sheet("Feature table", rows, flagged=lambda r: bool(r["gap_filled"]), fill=filled_fill)

    if provenance is not None:
        ws = wb.create_sheet("Provenance")
        for key, value in provenance.to_dict().items():
            ws.append([key, json.dumps(value) if isinstance(value, (dict, list)) else value])
        ws.column_dimensions["A"].width = 14
        ws.column_dimensions["B"].width = 80

    if not wb.sheetnames:
        wb.create_sheet("Empty")
    wb.save(path)
    return path


# ------------------------------------------------------------- figures
def _agg_pyplot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_qc_chart(values, mean: float | None = None, sd: float | None = None,
                  label: str = "value", width_mm: float = 89.0):
    """A Westgard-annotated control chart: the mean, +-2SD/+-3SD lines, and
    a ring around every flagged point. Sized for a single-column figure by
    default (89 mm); pass `width_mm=183.0` for a double-column one."""
    plt = _agg_pyplot()
    values = np.asarray(values, dtype=float)
    mean = float(np.nanmean(values)) if mean is None else mean
    sd = ((float(np.nanstd(values, ddof=1)) if len(values) > 1 else 0.0) if sd is None else sd)
    flags = westgard_flags(values, mean=mean, sd=sd)
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, width_mm / 25.4 * 0.62), dpi=300)
    x = np.arange(len(values))
    ax.plot(x, values, "o-", ms=4, lw=0.6, color="#1f6f8b")
    for k, colour in ((2, "#e0a020"), (3, "#c2410c")):
        ax.axhline(mean + k * sd, color=colour, lw=0.6, ls=":")
        ax.axhline(mean - k * sd, color=colour, lw=0.6, ls=":")
    ax.axhline(mean, color="#888888", lw=0.6, ls="--")
    if flags:
        idx = sorted({f.index for f in flags if f.index < len(values)})
        ax.plot(x[idx], values[idx], "o", ms=7, mfc="none", mec="#c2410c", mew=1.2)
    ax.set_xlabel("run order", fontsize=8)
    ax.set_ylabel(label, fontsize=8)
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    return fig


def plot_chromatogram(image, t1_start_s, t2_s, t1_range_min=None, t2_range_s=None, levels=None,
                      scale: str = "log", label: str = "", width_mm: float = 89.0, displayed: bool = False):
    """A GC×GC image as a publication figure: boxed axes with inward ticks,
    Arial 7 to 8 pt, no title, a short horizontal colour bar, magma colours.
    `image` is (rows, modulations); `t1_start_s` gives each modulation's start
    (s) and `t2_s` each row's time (s). `t1_range_min` and `t2_range_s`
    crop to a window; `levels` is (low, high) in the displayed units; `scale`
    is 'log', 'sqrt' or 'linear'. `label` (for example an ion) is written
    inside the panel at the top left. Sized square-ish at `width_mm`. With
    `displayed=True` the image is taken to be in display units already (what
    the app's 2D view shows) and only labelled with `scale`."""
    plt = _agg_pyplot()
    plt.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 7})
    img = np.asarray(image, dtype=float)
    disp = img if displayed else {"log": lambda a: np.log10(1.0 + np.clip(a, 0, None)),
                                  "sqrt": lambda a: np.sqrt(np.clip(a, 0, None)),
                                  "linear": lambda a: a}[scale](img)
    t1 = (np.asarray(t1_start_s, dtype=float) + 0.0) / 60.0
    t2 = np.asarray(t2_s, dtype=float)
    finite = disp[np.isfinite(disp)]
    if levels is None:
        levels = (float(np.quantile(finite, 0.5)), float(np.quantile(finite, 0.9995)))
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, width_mm / 25.4 * 0.9), dpi=300)
    mesh = ax.pcolormesh(t1, t2, disp, cmap="magma", vmin=levels[0], vmax=levels[1], shading="nearest",
                         rasterized=True)
    if t1_range_min is not None:
        ax.set_xlim(*t1_range_min)
    if t2_range_s is not None:
        ax.set_ylim(*t2_range_s)
    ax.set_xlabel("First-dimension retention time (min)", fontsize=8)
    ax.set_ylabel("Second-dimension retention time (s)", fontsize=8)
    ax.tick_params(direction="in", top=True, right=True, length=3, width=0.5, labelsize=7)
    for spine in ax.spines.values():
        spine.set_linewidth(0.6)
    if label:
        ax.text(0.02, 0.97, label, transform=ax.transAxes, va="top", ha="left", fontsize=7, color="white")
    unit = {"log": "log$_{10}$(1 + signal)", "sqrt": "√signal", "linear": "signal"}[scale]
    cax = fig.add_axes([0.62, 0.9, 0.3, 0.025])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=6, length=2, width=0.4)
    cb.set_label(unit, fontsize=6, labelpad=1)
    cb.outline.set_linewidth(0.4)
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.12, top=0.98)
    return fig


def chromatogram_caption(name: str, ion: str, scale: str, t1_range_min, t2_range_s, period: float,
                         phase: float, processing: str) -> str:
    """The sibling caption text for a `plot_chromatogram` figure, in the house
    journal style: what is shown, then the detail needed to reproduce it."""
    win = ""
    if t1_range_min is not None and t2_range_s is not None:
        win = (f" The window spans t₁ = {t1_range_min[0]:.1f} to {t1_range_min[1]:.1f} min and "
               f"t₂ = {t2_range_s[0]:.2f} to {t2_range_s[1]:.2f} s.")
    scale_txt = {"log": "log₁₀(1 + signal)", "sqrt": "square-root", "linear": "linear"}[scale]
    return (f"GC×GC–EI-ToF-MS chromatogram of {name}, showing {ion}.{win} Colours give the signal on a "
            f"{scale_txt} scale. The modulation period was {period:g} s and the phase {phase:.3f} s. "
            f"Data: {processing}. [Add site, period of sampling, instrument settings and the colour-scale "
            f"limits before submission.]\n")


def plot_pca(result: PCAResult, width_mm: float = 89.0):
    """A PC1/PC2 scatter with every run named, sized for publication."""
    plt = _agg_pyplot()
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, width_mm / 25.4), dpi=300)
    scores = result.scores
    ax.scatter(scores[:, 0], scores[:, 1], s=28, color="#1f6f8b", edgecolor="white", linewidth=0.4)
    for (x, y), name in zip(scores[:, :2], result.run_names):
        ax.annotate(name, (x, y), fontsize=6, xytext=(3, 3), textcoords="offset points")
    ev = result.explained_variance
    ax.set_xlabel(f"PC1 ({ev[0]:.0%})" if ev.size else "PC1", fontsize=8)
    ax.set_ylabel(f"PC2 ({ev[1]:.0%})" if ev.size > 1 else "PC2", fontsize=8)
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    return fig


def plot_volcano(points: list[VolcanoPoint], q_threshold: float = 0.05, width_mm: float = 89.0):
    """log2 fold change against -log10(q), significant points picked out."""
    plt = _agg_pyplot()
    fig, ax = plt.subplots(figsize=(width_mm / 25.4, width_mm / 25.4 * 0.8), dpi=300)
    x = np.array([p.log2_fold_change for p in points])
    q = np.clip(np.array([p.q_value for p in points]), 1e-300, None)
    sig = q < q_threshold
    ax.scatter(x[~sig], -np.log10(q[~sig]), s=14, color="#7e8ea6", label="ns")
    ax.scatter(x[sig], -np.log10(q[sig]), s=14, color="#e08420", label="significant")
    ax.axhline(-np.log10(q_threshold), color="#888888", lw=0.6, ls="--")
    ax.set_xlabel("log$_2$ fold change", fontsize=8)
    ax.set_ylabel("$-$log$_{10}(q)$", fontsize=8)
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    return fig


def save_figure(fig, path, dpi: int = 300) -> Path:
    """Save and close a figure. SVG or PNG, chosen from the path's suffix."""
    path = Path(path)
    fig.savefig(path, dpi=dpi)
    _agg_pyplot().close(fig)
    return path


def _fig_to_base64_png(fig, dpi: int = 200) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi)
    _agg_pyplot().close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")


# ------------------------------------------------------------- HTML
_HTML_STYLE = """
body { font-family: -apple-system, "Segoe UI", Arial, sans-serif; margin: 24px; color: #1a1a1a; }
h1 { font-size: 1.3em; }
h2 { font-size: 1.05em; margin-top: 1.6em; border-bottom: 1px solid #ccc; padding-bottom: 2px; }
table { border-collapse: collapse; font-size: 0.82em; margin: 0.5em 0 1.2em; }
th, td { border: 1px solid #ddd; padding: 3px 7px; text-align: right; }
th { background: #f2f2f2; text-align: center; }
td:first-child, th:first-child { text-align: left; }
tr.flag { background: #fde2c8; }
tr.filled { background: #dcebf7; }
.meta { color: #555; font-size: 0.85em; }
img { max-width: 100%; }
"""


def _esc(v) -> str:
    v = _clean(v)
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


def _html_table(rows: list, row_classes: list | None = None) -> str:
    if not rows:
        return "<p><em>(none)</em></p>"
    fields = list(rows[0].keys())
    head = "".join(f"<th>{f}</th>" for f in fields)
    body = []
    for i, row in enumerate(rows):
        cls = f' class="{row_classes[i]}"' if row_classes and row_classes[i] else ""
        cells = "".join(f"<td>{_esc(row[k])}</td>" for k in fields)
        body.append(f"<tr{cls}>{cells}</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table>"


def write_html_report(path, provenance: Provenance, *, report: QuantReport | None = None,
                      feature_table: list | None = None, run_names: list | None = None,
                      qc_values=None, qc_label: str = "", pca_result: PCAResult | None = None,
                      volcano_points: list[VolcanoPoint] | None = None,
                      title: str = "PyGCGC batch report") -> Path:
    """A self-contained HTML page: provenance, whatever tables were given,
    and one chart per statistic that was given data to draw. Every figure
    is a plain embedded PNG (base64), so the file needs nothing else to
    open, now or years from now. Needs the `plot` extra (matplotlib) for
    any of `qc_values`/`pca_result`/`volcano_points`; tables alone need
    neither extra."""
    path = Path(path)
    parts = [f"<!doctype html><meta charset='utf-8'><title>{title}</title><style>{_HTML_STYLE}</style>",
            f"<h1>{title}</h1>",
            f"<p class='meta'>PyGCGC {provenance.version} &middot; generated {provenance.generated}"
            + (f" &middot; modulation period {provenance.period:g} s" if provenance.period else "")
            + "</p>"]
    if provenance.sources:
        parts.append('<h2>Input sources</h2><pre>' + html.escape(json.dumps(provenance.sources, indent=2)) + '</pre>')

    if provenance.edits:
        parts.append("<h2>Manual edits</h2><ul>")
        for key, val in provenance.edits.items():
            parts.append(f"<li><b>{key}</b>: {json.dumps(val)}</li>")
        parts.append("</ul>")
    else:
        parts.append("<p class='meta'>No manual edits recorded for this run.</p>")

    if feature_table:
        names = run_names or sorted({n for f in feature_table for n in f.volumes})
        rows = []
        for f in feature_table:
            row = {"feature": f.id, "t1 (min)": round(f.ref_t1 / 60.0, 3), "t2 (s)": round(f.ref_t2, 3)}
            row.update({n: f.volumes.get(n) for n in names})
            rows.append(row)
        classes = ["filled" if any(f.filled.values()) else "" for f in feature_table]
        parts.append(f"<h2>Feature table</h2>{_html_table(rows, classes)}")

    if report is not None:
        rows = [asdict(r) for r in report.results]
        classes = ["flag" if (r.status in ("check", "rejected") or r.flags) else ""
                  for r in report.results]
        parts.append(f"<h2>Targets</h2>{_html_table(rows, classes)}")
        cal_rows = calibration_table(report.fits)
        parts.append(f"<h2>Calibrations</h2>"
                     f"{_html_table(cal_rows, ['' if r['passes'] else 'flag' for r in cal_rows])}")
        istd_rows = [asdict(r) for r in report.istd]
        parts.append(f"<h2>Internal standards</h2>"
                     f"{_html_table(istd_rows, ['flag' if r.flag else '' for r in report.istd])}")

    figures = []
    if qc_values is not None and len(qc_values):
        figures.append((f"QC chart -- {qc_label}" if qc_label else "QC chart",
                        plot_qc_chart(qc_values, label=qc_label or "value")))
    if pca_result is not None:
        figures.append(("PCA", plot_pca(pca_result)))
    if volcano_points:
        figures.append(("Volcano", plot_volcano(volcano_points)))
    for heading, fig in figures:
        parts.append(f"<h2>{heading}</h2><img src='data:image/png;base64,{_fig_to_base64_png(fig)}'>")

    path.write_text("".join(parts), encoding="utf-8")
    return path


# ------------------------------------------------------------- PDF
def write_pdf_report(path, provenance: Provenance, *, report: QuantReport | None = None,
                     feature_table: list | None = None, qc_values=None, qc_label: str = "",
                     pca_result: PCAResult | None = None, volcano_points: list[VolcanoPoint] | None = None,
                     title: str = "PyGCGC batch report") -> Path:
    """A short, print-ready PDF: one cover page of provenance and counts,
    then one page per chart that was given data to draw. Needs the `plot`
    extra (matplotlib)."""
    plt = _agg_pyplot()
    from matplotlib.backends.backend_pdf import PdfPages

    path = Path(path)
    lines = [title, "", f"PyGCGC {provenance.version}", f"Generated {provenance.generated}"]
    if provenance.period:
        lines.append(f"Modulation period: {provenance.period:g} s")
    if report is not None:
        n_ok = sum(r.status == "ok" for r in report.results)
        lines.append(f"Targets: {len(report.results)} results, {n_ok} ok")
    if feature_table:
        lines.append(f"Feature table: {len(feature_table)} features")
    lines.append("")
    lines.append("Manual edits:" if provenance.edits else "No manual edits recorded.")
    for key, val in provenance.edits.items():
        lines.append(f"  {key}: {json.dumps(val)}")

    with PdfPages(path) as pdf:
        fig, ax = plt.subplots(figsize=(8.27, 5.83), dpi=200)
        ax.axis("off")
        ax.text(0.03, 0.97, "\n".join(lines), va="top", family="monospace", fontsize=9)
        pdf.savefig(fig)
        plt.close(fig)
        for source in provenance.sources:
            fig, ax = plt.subplots(figsize=(8.27, 5.83), dpi=200)
            ax.axis('off')
            text = '\n'.join(textwrap.fill(line, width=95) for line in json.dumps(source, indent=2).splitlines())
            ax.text(.03, .97, 'Input source\n\n'+text, va='top', family='monospace', fontsize=8)
            pdf.savefig(fig)
            plt.close(fig)
        if qc_values is not None and len(qc_values):
            fig = plot_qc_chart(qc_values, label=qc_label or "value", width_mm=170.0)
            pdf.savefig(fig)
            plt.close(fig)
        if pca_result is not None:
            fig = plot_pca(pca_result, width_mm=170.0)
            pdf.savefig(fig)
            plt.close(fig)
        if volcano_points:
            fig = plot_volcano(volcano_points, width_mm=170.0)
            pdf.savefig(fig)
            plt.close(fig)
    return path
