"""Export the available analysis without changing the project's save location."""

from __future__ import annotations

import csv
from dataclasses import asdict, is_dataclass
import json
from pathlib import Path

import numpy as np

from .. import __version__
from ..core.components import write_components_csv
from ..core.chemclass import write_classes_csv
from ..core.cluster import write_clusters_csv
from ..core.istd import write_istd_report_csv
from ..core.libmatch import write_identifications_csv
from ..core.peak_blank import write_rows_csv as write_peak_blank_csv
from ..core.report import build_provenance, write_feature_table_csv
from ..core.targets import write_targets_csv


def _json_value(value):
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot export {type(value).__name__}")


def export_project(project, folder, *, identifications=None) -> list[Path]:
    """Write a project snapshot, all available tables and explicit freshness flags.

    Raw files stay at their original paths. The exported project points to the
    exported target and level tables so those inputs travel with the results.
    Names of per-run output files use row numbers to avoid collisions between
    equal sample names or characters that Windows cannot put in a filename.
    """
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    paths = []

    def write_json(name, value):
        path = folder / name
        path.write_text(json.dumps(value, indent=2, default=_json_value), encoding="utf-8")
        paths.append(path)

    snapshot = project.to_dict()
    for saved, entry in zip(snapshot["runs"], project.entries):
        saved["path"] = str(entry.path.resolve())
    sequence = folder / "sequence.csv"
    with sequence.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["file", "type", "level", "cal_set", "volume_m3"])
        writer.writerows([entry.path.name, entry.kind, entry.level, entry.cal_set,
                          "" if entry.volume is None else entry.volume]
                         for entry in project.entries)
    paths.append(sequence)
    snapshot["sequence"] = "sequence.csv"
    if project.quant.targets:
        paths.append(write_targets_csv(project.quant.targets, folder / "target_table.csv"))
        snapshot["targets"] = "target_table.csv"
    if project.quant.levels:
        levels = project.quant.levels
        compounds = sorted({name for amounts in levels.values() for name in amounts})
        path = folder / "calibration_levels.csv"
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=["level", *compounds])
            writer.writeheader()
            writer.writerows({"level": level, **amounts} for level, amounts in levels.items())
        paths.append(path)
        snapshot["quant"]["levels"] = "calibration_levels.csv"
    if project.quant.report is not None:
        paths.extend(project.quant.report.write(folder))
    matches = project.istd.all_matches()
    if matches:
        paths.append(write_istd_report_csv(matches, folder / "internal_standard_matches.csv", project.istd.settings))

    manifest = {"version": __version__, "raw_files_included": False,
                "amount_unit": project.quant.amount_unit,
                "concentration_unit": f"{project.quant.amount_unit}/m³",
                "quantification_available": project.quant.report is not None,
                "quantification_stale": project.quant.stale, "runs": []}
    write_json("calibration_setup.json", {"amount_unit": project.quant.amount_unit,
              "concentration_unit": f"{project.quant.amount_unit}/m³",
              **project.quant.calibration_setup})
    write_json("standard_review_decisions.json", project.istd.review_decisions)
    manifest['standard_exceptions_remaining'] = len(project.istd.exceptions())
    assignments = {}
    for i, entry in enumerate(project.entries):
        components = project.blobs.components_for(i)
        blobs = project.blobs.blobs_for(i)
        row = {"run": entry.name, "source": str(entry.path), "import_error": entry.error,
               "standard_matches_stale": project.istd.is_stale(entry.name),
               "peaks_stale": project.blobs.selected_stale(entry),
               "classes_stale": entry in project.classes.stale}
        if project.blobs.method == "ms-aware" and components is not None:
            path = folder / f"peaks_{i + 1:03d}.csv"
            paths.append(write_components_csv(components, path, entry.name))
            row["peaks"] = path.name
        elif project.blobs.method == "legacy" and blobs is not None:
            path = folder / f"peaks_{i + 1:03d}.csv"
            with path.open("w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(["run", "peak", "t1_s", "t2_s", "volume", "snr", "co_elution"])
                writer.writerows([entry.name, j + 1, blobs.t1[j], blobs.t2[j], blobs.volume[j],
                                  blobs.snr[j], bool(blobs.co_elution[j])] for j in range(len(blobs)))
            paths.append(path)
            row["peaks"] = path.name
        matched = project.matcher.results_for(i)
        if matched is not None and matched.matches:
            library = project.matcher.library.label if project.matcher.library is not None else ""
            paths.append(write_identifications_csv(folder / f"identifications_{i + 1:03d}.csv", entry.name, matched.matches,
                                                   project.matcher.decisions_for(i), library, project.matcher.settings))
            row["identifications"] = f"identifications_{i + 1:03d}.csv"
            row["identifications_stale"] = project.matcher.current_for(i) is None
        subtraction = project.blank.results_for(i)
        if subtraction is not None and subtraction.rows:
            paths.append(write_peak_blank_csv(folder / f"peak_blank_{i + 1:03d}.csv", entry.name, subtraction.rows,
                                              project.blank.settings))
            row["peak_blank"] = f"peak_blank_{i + 1:03d}.csv"
            row["peak_blank_stale"] = project.blank.current_for(i) is None
        classes = project.classes.result_for(i)
        classified_peaks = project.classes.blobs_classified(i)
        if classes is not None and classified_peaks is not None:
            paths.append(write_classes_csv(classes, classified_peaks, folder / f"classes_{i + 1:03d}.csv"))
        assignments[str(i + 1)] = {"run": entry.name, "classes": classes,
                                   "hydrocarbons": project.hydrocarbons.calls_for(i)}
        manifest["runs"].append(row)
    assignments["clusters"] = {"summary": project.clusters.summary,
                               "peaks": project.clusters.table, "assignments": project.clusters.result,
                               "regions": project.clusters.region_summaries,
                               "region_shares": project.clusters.region_shares,
                               "names": project.clusters.names,
                               "stale": project.clusters.stale}
    if project.clusters.result is not None and project.clusters.table is not None:
        paths.append(write_clusters_csv(project.clusters.table, project.clusters.result,
                                        folder / "peak_clusters.csv", names=project.clusters.name_map()))
    fingerprint = project.clusters.fingerprint()
    if fingerprint is not None:
        paths.append(write_feature_table_csv(fingerprint[0], project.clusters.active_runs(),
                                             folder / "cluster_run_table.csv"))
    if project.clusters.region_labels:
        labels = {f"run_{i + 1:03d}": project.clusters.region_labels[entry.name]
                  for i, entry in enumerate(project.entries) if entry.name in project.clusters.region_labels}
        path = folder / "region_labels.npz"
        np.savez_compressed(path, **labels)
        paths.append(path)
        assignments["clusters"]["region_label_runs"] = {
            f"run_{i + 1:03d}": entry.name for i, entry in enumerate(project.entries)
            if entry.name in project.clusters.region_labels}
    if identifications:
        assignments["library_search"] = identifications
    write_json("assignments.json", assignments)
    write_json("project.pygcgc.json", snapshot)
    provenance = build_provenance(
        __version__, period=project.period,
        method={"modulation_period_confirmed": project.period_confirmed},
        settings={"background": project.bg_settings.to_dict(), "quant": project.quant.settings.to_dict(),
                  "peaks": project.blobs.component_settings.to_dict(),
                  "peak_blank": dict(project.blank.settings.to_dict(), applied_to_quantification=project.blank.apply_to_quant),
                  "library_match": dict(project.matcher.settings.to_dict(), library=project.matcher.source,
                                        alkane_ladder_run=project.matcher.ladder_run or "each run's own")},
        blob_edits=project.blobs.edits, quant_decisions=project.quant.decisions,
        stencils=project.hydrocarbons.stencils,
        runs=[entry.run for entry in project.entries if entry.run is not None])
    paths.append(provenance.write(folder / "provenance.json"))
    manifest["files"] = [path.name for path in paths]
    write_json("manifest.json", manifest)
    return paths
