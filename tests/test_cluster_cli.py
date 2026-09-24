"""`pygcgc batch --cluster`: peak clusters and the run table from a whole folder."""

import csv

import pytest

from tests.conftest import GUI_PHASE, make_run, write_run_andi

ALKANE_IONS = (57, 71, 85, 99, 113)
AROMATIC_IONS = (65, 77, 91, 105, 119)
WEIGHTS = (1.0, 0.7, 0.4, 0.2, 0.1)
N_AROMATIC = {"A1": 3, "A2": 3, "B1": 8, "B2": 8}     # group B holds far more of the aromatic family


def _peaks(n_alkane, n_aromatic):
    peaks = []
    for j in range(n_alkane + n_aromatic):
        ions = ALKANE_IONS if j < n_alkane else AROMATIC_IONS
        t1, t2 = 61.3 + 55.0 * j, 1.207 + 0.11 * (j % 12)         # off the sampling grid
        peaks += [(t1, t2, ion, 2e5 * w) for ion, w in zip(ions, WEIGHTS)]
    return peaks


@pytest.fixture
def family_folder(tmp_path):
    folder = tmp_path / "family_pack"
    (folder / "data").mkdir(parents=True)
    for k, (name, n_aromatic) in enumerate(N_AROMATIC.items()):
        run = make_run(_peaks(8, n_aromatic), phase=GUI_PHASE, rate=50.0, mz=(50, 130), duration=1100.0, seed=k + 1)
        write_run_andi(run, folder / "data" / f"{name}.cdf")
    (folder / "sequence.csv").write_text("file,type\n" + "".join(f"{n}.cdf,sample\n" for n in N_AROMATIC))
    (folder / "method.json").write_text('{"modulation_period_s": 4.0}')
    return folder


def test_batch_cluster_writes_peak_and_run_tables_and_a_permanova(family_folder, tmp_path, capsys):
    from pygcgc import BatchConfig, BatchSettings
    from pygcgc.cli import main
    config = tmp_path / "pipeline.json"
    BatchConfig(settings=BatchSettings(run_quant=False), group_a=("A1", "A2"), group_b=("B1", "B2")).save(config)
    out = tmp_path / "out"
    assert main(["batch", str(family_folder), "--config", str(config), "--outdir", str(out), "--workers", "1",
                 "--cluster"]) == 0
    text = capsys.readouterr().out
    assert "2 clusters from" in text and "4 run(s)" in text
    assert "PERMANOVA of the two groups" in text
    with open(out / "cluster_run_table.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2 and {"A1", "A2", "B1", "B2"} <= set(rows[0])
    aromatic_row = min(rows, key=lambda r: float(r["A1"]))          # the cluster group A holds least of
    assert float(aromatic_row["B1"]) > float(aromatic_row["A1"])
    assert len((out / "peak_clusters.csv").read_text().splitlines()) > 20


def test_batch_regions_writes_labels_and_a_share_table(family_folder, tmp_path, capsys):
    pytest.importorskip("sklearn")
    import numpy as np
    from pygcgc.cli import main
    out = tmp_path / "out"
    assert main(["batch", str(family_folder), "--outdir", str(out), "--workers", "1", "--no-quant", "--no-blobs",
                 "--regions", "2"]) == 0
    text = capsys.readouterr().out
    assert "2 pixel regions fitted on 4 run(s)" in text
    labels = np.load(out / "region_labels_A1.npy")
    assert labels.ndim == 2 and labels.max() == 1 and labels.min() == -1
    with open(out / "region_run_table.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2 and all(abs(sum(float(r[n]) for r in rows) - 100.0) < 1e-6 for n in ("A1", "B2"))


def test_batch_without_the_flag_writes_no_cluster_files(family_folder, tmp_path):
    from pygcgc.cli import main
    out = tmp_path / "out"
    assert main(["batch", str(family_folder), "--outdir", str(out), "--workers", "1", "--no-quant"]) == 0
    assert not (out / "cluster_run_table.csv").exists()


def test_batch_cluster_sheet_goes_into_the_workbook(family_folder, tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    from pygcgc.cli import main
    out = tmp_path / "out"
    assert main(["batch", str(family_folder), "--outdir", str(out), "--workers", "1", "--no-quant", "--cluster",
                 "--report", "xlsx"]) == 0
    wb = openpyxl.load_workbook(out / "report.xlsx")
    assert "Clusters" in wb.sheetnames
    header = [c.value for c in wb["Clusters"][1]]
    assert header[:3] == ["cluster", "t1_min", "t2_s"] and "B2" in header
