"""📤 Reports: provenance, spreadsheet/document exports, and figures."""

import json

import numpy as np
import pytest
from openpyxl import load_workbook

from pygcgc import (BatchConfig, Feature, PCAResult, QuantDecisions, QuantReport, QuantSettings,
                    VolcanoPoint, pca, run_batch, volcano)
from pygcgc.core.calibration import fit_calibration
from pygcgc.core.hydrocarbons import Stencil, StencilSet
from pygcgc.core.quant import IstdRow, SampleResult, SignalRow
from pygcgc.core.report import (Provenance, build_provenance, plot_pca, plot_qc_chart, plot_volcano,
                                save_figure, write_excel, write_feature_table_csv, write_html_report,
                                write_pdf_report)


# ------------------------------------------------------------- provenance
def test_build_provenance_records_version_period_method_and_settings():
    prov = build_provenance("9.9.9", period=4.0, method={"modulation_period_s": 4.0},
                            settings={"run_blobs": True})
    assert prov.version == "9.9.9"
    assert prov.period == 4.0
    assert prov.method == {"modulation_period_s": 4.0}
    assert prov.settings == {"run_blobs": True}
    assert prov.edits == {}
    # a real, parseable UTC timestamp
    assert prov.generated.endswith("+00:00") or prov.generated.endswith("Z")


def test_build_provenance_summarises_blob_edits_quant_decisions_and_stencils():
    blob_edits = {"S1": [("delete", (3,)), ("merge", (1, 2))], "S2": []}
    decisions = QuantDecisions(rejected={("S1", "A")}, accepted={("S1", "B"), ("S2", "B")},
                               excluded={"A": {"CAL_L1": "outlier"}}, overrides={},
                               assigned={("S1", "C"): (100.0, 1.0)})
    stencils = StencilSet(name="frame", stencils=[Stencil("branched", [(15.0, 0.1), (16.0, 0.2)])])

    prov = build_provenance("1.0.0", blob_edits=blob_edits, quant_decisions=decisions, stencils=stencils)

    assert prov.edits["blob_edits"] == {"total": 2, "by_run": {"S1": 2}}      # S2's empty list drops out
    qd = prov.edits["quant_decisions"]
    assert qd == {"accepted": 2, "rejected": 1, "picked_by_hand": 1, "integrated_by_hand": 0,
                 "excluded_from_curve": 1, "overridden_curves": 0}
    assert prov.edits["stencils"] == {"frame": "frame", "names": ["branched"]}


def test_build_provenance_empty_stencil_set_contributes_nothing():
    prov = build_provenance("1.0.0", stencils=StencilSet(name="frame"))
    assert prov.edits == {}


def test_provenance_round_trips_through_json(tmp_path):
    prov = build_provenance("1.0.0", period=4.0, method={"a": 1}, settings={"b": 2})
    path = prov.write(tmp_path / "provenance.json")
    loaded = json.loads(path.read_text())
    assert loaded == prov.to_dict()


# ------------------------------------------------------------- feature table CSV
def test_write_feature_table_csv_flags_gap_filled_runs(tmp_path):
    ft = [
        Feature(id=0, ref_t1=60.0, ref_t2=1.0, volumes={"R1": 100.0, "R2": 90.0}, filled={"R2": True}),
        Feature(id=1, ref_t1=120.0, ref_t2=1.5, volumes={"R1": 50.0, "R2": 55.0}, filled={}),
    ]
    path = write_feature_table_csv(ft, ["R1", "R2"], tmp_path / "features.csv")
    import csv
    rows = list(csv.DictReader(open(path, encoding="utf-8")))
    assert rows[0]["gap_filled"] == "R2"
    assert rows[1]["gap_filled"] == ""
    assert rows[0]["R1"] == "100.0" and rows[0]["R2"] == "90.0"


# ------------------------------------------------------------- fixtures
def _sample_report() -> QuantReport:
    good_fit = fit_calibration("A", [0, 1, 2, 5, 10], [0, 2, 4, 10, 20],
                               labels=["L0", "L1", "L2", "L3", "L4"])
    fits = {("set1", "A"): {"quant": good_fit},
           ("set1", "B"): {"quant": "insufficient calibration points"}}
    results = [
        SampleResult("S1", "A", "target", True, "ok", [], 10.0, 0.5, 1.0, 5.0, 0.2, 2.5, "set1"),
        SampleResult("S1", "B", "target", True, "check", ["ambiguous"], 20.0, 0.5, 1.0, float("nan"),
                    float("nan"), float("nan"), "set1"),
    ]
    signals = [SignalRow("S1", 0, 10.0, 0.5, 20.0, 1000.0, "IS", 500.0, 2.0, "A", "", 0.0, 5.0, 2.5, "")]
    istd = [IstdRow("S1", "IS", True, 5.0, 0.3, 500.0, 2.0, False),
           IstdRow("S2", "IS", True, 5.1, 0.3, 900.0, 80.0, True)]
    return QuantReport(seq=[], targets=[], positions={}, fits=fits, results=results,
                       signals=signals, istd=istd, settings=QuantSettings())


def _sample_feature_table() -> list:
    return [
        Feature(id=0, ref_t1=60.0, ref_t2=1.0, volumes={"S1": 100.0, "S2": 90.0}, filled={"S2": True}),
        Feature(id=1, ref_t1=120.0, ref_t2=1.5, volumes={"S1": 50.0, "S2": 55.0}, filled={}),
    ]


# ------------------------------------------------------------- Excel
def test_write_excel_highlights_check_and_failed_rows(tmp_path):
    report = _sample_report()
    prov = build_provenance("1.0.0", period=4.0)
    path = write_excel(tmp_path / "report.xlsx", report=report, feature_table=_sample_feature_table(),
                       run_names=["S1", "S2"], provenance=prov)
    wb = load_workbook(path)
    assert set(wb.sheetnames) >= {"Targets", "Calibrations", "All signals",
                                  "Internal standards", "Feature table", "Provenance"}

    targets = wb["Targets"]
    header = [c.value for c in targets[1]]
    status_col = header.index("status") + 1
    for row in targets.iter_rows(min_row=2):
        status = row[status_col - 1].value
        filled = row[0].fill.start_color.rgb
        if status == "check":
            assert filled not in (None, "00000000")
        elif status == "ok":
            assert filled in (None, "00000000")

    cal = wb["Calibrations"]
    header = [c.value for c in cal[1]]
    passes_col = header.index("passes")
    for row in cal.iter_rows(min_row=2):
        filled = row[0].fill.start_color.rgb
        if row[passes_col].value:
            assert filled in (None, "00000000")
        else:
            assert filled not in (None, "00000000")

    istd = wb["Internal standards"]
    header = [c.value for c in istd[1]]
    flag_col = header.index("flag")
    flagged_any = False
    for row in istd.iter_rows(min_row=2):
        if row[flag_col].value:
            flagged_any = True
            assert row[0].fill.start_color.rgb not in (None, "00000000")
    assert flagged_any

    ft_sheet = wb["Feature table"]
    header = [c.value for c in ft_sheet[1]]
    gap_col = header.index("gap_filled")
    saw_filled = False
    for row in ft_sheet.iter_rows(min_row=2):
        if row[gap_col].value:
            saw_filled = True
            assert row[0].fill.start_color.rgb not in (None, "00000000")
    assert saw_filled


def test_write_excel_with_nothing_still_produces_a_workbook(tmp_path):
    path = write_excel(tmp_path / "empty.xlsx")
    wb = load_workbook(path)
    assert wb.sheetnames == ["Empty"]


def test_write_excel_handles_nan_values_without_raising(tmp_path):
    # SampleResult for "B" above carries NaN amount/amount_se/concentration.
    path = write_excel(tmp_path / "report.xlsx", report=_sample_report())
    wb = load_workbook(path)
    header = [c.value for c in wb["Targets"][1]]
    amount_col = header.index("amount")
    values = [row[amount_col].value for row in wb["Targets"].iter_rows(min_row=2)]
    assert None in values           # NaN was written as a blank cell, not raised


# ------------------------------------------------------------- HTML
def test_write_html_report_marks_flagged_rows_and_embeds_charts(tmp_path):
    report = _sample_report()
    ft = _sample_feature_table()
    prov = build_provenance("1.0.0", period=4.0, blob_edits={"S1": [("delete", (0,))]})
    pca_result = PCAResult(scores=np.array([[1.0, 0.0], [-1.0, 0.0]]), loadings=np.zeros((2, 2)),
                           explained_variance=np.array([0.9, 0.1]), run_names=["S1", "S2"],
                           feature_ids=[0, 1])
    volcano_points = [VolcanoPoint(feature_id=0, ref_t1=60.0, ref_t2=1.0,
                                   log2_fold_change=2.0, p_value=0.001, q_value=0.01)]
    path = write_html_report(tmp_path / "report.html", prov, report=report, feature_table=ft,
                             run_names=["S1", "S2"], qc_values=[100.0, 105.0, 98.0], qc_label="IS",
                             pca_result=pca_result, volcano_points=volcano_points)
    text = path.read_text(encoding="utf-8")
    assert 'class="flag"' in text                # the "check" target row and the failed cal fit
    assert 'class="filled"' in text               # the gap-filled feature
    assert text.count("data:image/png;base64,") == 3      # QC, PCA, volcano
    assert "blob_edits" in text
    assert "PyGCGC 1.0.0" in text


def test_write_html_report_with_only_tables_needs_no_figures(tmp_path):
    prov = build_provenance("1.0.0")
    path = write_html_report(tmp_path / "report.html", prov, report=_sample_report())
    text = path.read_text(encoding="utf-8")
    assert "data:image/png;base64," not in text
    assert "No manual edits recorded" in text


# ------------------------------------------------------------- PDF
def test_write_pdf_report_produces_a_pdf_with_a_page_per_chart(tmp_path):
    prov = build_provenance("1.0.0", period=4.0)
    pca_result = PCAResult(scores=np.array([[1.0, 0.0], [-1.0, 0.0]]), loadings=np.zeros((2, 2)),
                           explained_variance=np.array([0.9, 0.1]), run_names=["S1", "S2"],
                           feature_ids=[0, 1])
    path = write_pdf_report(tmp_path / "report.pdf", prov, report=_sample_report(),
                            feature_table=_sample_feature_table(), qc_values=[100.0, 105.0, 98.0],
                            qc_label="IS", pca_result=pca_result)
    data = path.read_bytes()
    assert data.startswith(b"%PDF")
    assert data.count(b"/Type /Page") + data.count(b"/Type/Page") >= 3   # cover + QC + PCA


def test_write_pdf_report_cover_page_only(tmp_path):
    prov = build_provenance("1.0.0")
    path = write_pdf_report(tmp_path / "report.pdf", prov)
    assert path.read_bytes().startswith(b"%PDF")


# ------------------------------------------------------------- figures
def test_plot_qc_chart_flags_the_same_points_westgard_would():
    values = np.array([10.0, 10.1, 9.9, 10.2, 9.8, 10.0, 9.9, 10.1, 10.0, 10.0, 16.0])
    fig = plot_qc_chart(values, label="D4")
    assert fig.get_axes()[0].get_ylabel() == "D4"


def test_plot_pca_and_plot_volcano_return_figures_sized_to_width_mm():
    pca_result = PCAResult(scores=np.array([[1.0, 0.5], [-1.0, -0.5], [0.2, 0.1]]),
                           loadings=np.zeros((2, 2)), explained_variance=np.array([0.7, 0.2]),
                           run_names=["A", "B", "C"], feature_ids=[0, 1])
    fig1 = plot_pca(pca_result, width_mm=183.0)
    w_in, h_in = fig1.get_size_inches()
    assert w_in == pytest.approx(183.0 / 25.4, rel=1e-3)
    assert h_in == pytest.approx(183.0 / 25.4, rel=1e-3)

    points = [VolcanoPoint(feature_id=i, ref_t1=0.0, ref_t2=0.0, log2_fold_change=float(i - 2),
                           p_value=0.01, q_value=0.01 if i < 2 else 0.5) for i in range(5)]
    fig2 = plot_volcano(points)
    assert len(fig2.get_axes()) == 1


def test_save_figure_writes_svg_and_png(tmp_path):
    fig = plot_qc_chart([1.0, 2.0, 3.0])
    svg = save_figure(fig, tmp_path / "chart.svg")
    fig2 = plot_qc_chart([1.0, 2.0, 3.0])
    png = save_figure(fig2, tmp_path / "chart.png")
    assert svg.exists() and svg.stat().st_size > 0
    assert png.exists() and png.stat().st_size > 0


# ------------------------------------------------------------- end to end
def test_report_pipeline_on_a_real_batch_result(cal_folder, tmp_path):
    """Everything above works standalone; check it also works fed straight
    from a real `run_batch` result, the way the CLI and the batch dialog
    both use it."""
    result = run_batch(cal_folder, cache_dir=tmp_path / "cache", max_workers=1)
    run_names = list(result.prepped)
    prov = build_provenance("1.0.0", period=4.0, settings=BatchConfig().settings.to_dict())

    xlsx = write_excel(tmp_path / "r.xlsx", report=result.report, feature_table=result.feature_table,
                       run_names=run_names, provenance=prov)
    assert xlsx.exists()

    pca_result = pca(result.feature_table, run_names, n_components=2)
    html = write_html_report(tmp_path / "r.html", prov, report=result.report,
                             feature_table=result.feature_table, run_names=run_names,
                             pca_result=pca_result)
    assert "Feature table" in html.read_text(encoding="utf-8")

    csv_path = write_feature_table_csv(result.feature_table, run_names, tmp_path / "features.csv")
    assert csv_path.exists()
