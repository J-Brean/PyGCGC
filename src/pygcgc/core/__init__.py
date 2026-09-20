"""PyGCGC core: all the science, no GUI."""

from .run import Run
from .io import crop_run, read_many, read_run, register_reader, write_andi
from .timing import TimingReport, check_timing, fill_gaps
from .fold import Fold, FoldSpec, fold
from .phase import (Anchor, PhaseEstimate, PhaseEstimationError, auto_phase,
                    phase_from_anchors, phase_from_void, phase_for_sequence, locate_anchor)
from .period import (PeriodEstimate, estimate_period, estimate_period_for_runs, period_from_anchors,
                     period_from_points)
from .solvent import find_solvent_front
from .targets import Target, load_anchors, load_levels, load_targets, write_targets_csv
from .ionexpr import IonExpr, IonExpressionError, parse_ion_expr
from .background import (BackgroundResult, BackgroundSettings, NoiseMap, arpls, despike,
                         estimate_noise, remove_background)
from .align import (AlignmentError, AlignmentTransform, AnchorPoint, cross_validate,
                    fit_alignment, propagate_anchors, resample_to_reference)

from .peaks import (BlobSet, PeakSettings, add_blob_at, delete_blobs, detect_ceiling, find_blobs,
                    measure_ion_peak, merge_blobs)
from .candidates import (CandidateRegion, CandidateSet, CandidateSettings, candidate_cube,
                         candidate_regions, find_candidates)
from .components import (Component, ComponentDetectionSettings, ComponentGroup, ComponentSet,
                         detect_components, fit_candidate_region, group_from_fit)
from .pixelspec import ratios_to, spectrum_at
from .siloxane import (BLEED_REFERENCE, CYCLIC_REFERENCE, PRESETS as SILOXANE_PRESETS, SILOXANE_IONS,
                       SILOXANE_REFERENCE, SiloxaneResult, SiloxaneSettings, flag_siloxane_blobs,
                       reference_from_spectrum, remove_siloxane_blobs, remove_siloxanes, siloxane_score,
                       siloxane_scores)
from .hydrocarbons import (AlkaneHit, AlkaneSeries, HomologSeries, HydrocarbonSettings, Stencil,
                          StencilSet, find_alkane_series, find_homologous_series, find_markers,
                          propagate_alkane_series, qualifier_ratios_at, ratios_match,
                          stencil_from_polygon, stencil_totals)
from .batch import (BatchConfig, BatchResult, BatchSettings, Feature, PCAResult, PrepResult, QCFlag,
                   StepCache, VolcanoPoint, benjamini_hochberg, build_feature_table,
                   feature_table_matrix, pca, run_batch, volcano, westgard_flags)
from .calibration import CalibrationError, CalibrationFit, CalibrationSettings, fit_calibration
from .quant import (QuantDecisions, QuantReport, QuantSettings, RunQuant, SequenceRun, TargetHit, all_signals,
                    build_calibrations, cal_set_for, bracket_cal_sets, CalBracket, blended_amount,
                    istd_table, learn_positions, quantify_run,
                    quantify_sequence, read_sequence_extras, with_positions,
                    finish_sequence, quantify_prepared)
from .report import (Provenance, build_provenance, chromatogram_caption, plot_chromatogram, plot_pca,
                     plot_qc_chart, plot_volcano,
                     save_figure, write_excel, write_feature_table_csv, write_html_report,
                     write_pdf_report)
from .overlap import (BlobCube, ComponentAssessment, IonApex, MCRResult, MCRStability, OverlapComponent, OverlapCV, OverlapFit, OverlapSettings,
                      PARAFACResult, UniqueIonMap, blob_cube, build_run_tensor, choose_k,
                      assess_components, component_maps, core_consistency, cross_validate_overlap, fit_overlap, fit_template_overlap, mcr_als, mcr_restarts, parafac, predict_overlap, split_by_line,
                      split_by_valley, split_by_watershed, split_half_stability, split_volumes,
                      unique_ion_apexes)
from .identify import (Identification, IdentifySettings, Library, LibraryEntry, SearchHit,
                       assign_qvalues, build_inhouse_entry, build_library, classify,
                       generate_decoys, identify_blobs, identify_components, load_msp, promote_to_target,
                       retention_index, search, write_msp)

__all__ = [
    "Run", "read_run", "read_many", "register_reader", "write_andi", "crop_run",
    "TimingReport", "check_timing", "fill_gaps",
    "Fold", "FoldSpec", "fold",
    "Anchor", "PhaseEstimate", "PhaseEstimationError", "auto_phase",
    "phase_from_anchors", "phase_from_void", "phase_for_sequence", "locate_anchor",
    "PeriodEstimate", "estimate_period", "estimate_period_for_runs", "period_from_anchors", "period_from_points", "find_solvent_front",
    "load_anchors", "Target", "load_targets", "load_levels", "write_targets_csv",
    "IonExpr", "IonExpressionError", "parse_ion_expr",
    "BackgroundResult", "BackgroundSettings", "NoiseMap", "arpls", "despike",
    "estimate_noise", "remove_background",
    "AlignmentError", "AlignmentTransform", "AnchorPoint", "cross_validate",
    "fit_alignment", "propagate_anchors", "resample_to_reference",
    "PeakSettings", "BlobSet", "find_blobs", "measure_ion_peak", "detect_ceiling",
    "delete_blobs", "merge_blobs", "add_blob_at",
    "CandidateSettings", "CandidateSet", "CandidateRegion", "find_candidates",
    "candidate_regions", "candidate_cube",
    "Component", "ComponentDetectionSettings", "ComponentGroup", "ComponentSet",
    "detect_components", "fit_candidate_region", "group_from_fit",
    "SILOXANE_IONS", "SILOXANE_REFERENCE", "BLEED_REFERENCE", "CYCLIC_REFERENCE", "SILOXANE_PRESETS",
    "SiloxaneSettings", "SiloxaneResult", "reference_from_spectrum", "spectrum_at", "ratios_to",
    "siloxane_score", "siloxane_scores", "flag_siloxane_blobs", "remove_siloxanes", "remove_siloxane_blobs",
    "qualifier_ratios_at", "ratios_match",
    "HydrocarbonSettings", "AlkaneHit", "AlkaneSeries", "find_alkane_series",
    "propagate_alkane_series", "HomologSeries", "find_homologous_series", "find_markers",
    "Stencil", "StencilSet", "stencil_from_polygon", "stencil_totals",
    "BatchConfig", "BatchSettings", "BatchResult", "PrepResult", "StepCache", "run_batch",
    "Feature", "build_feature_table", "feature_table_matrix",
    "QCFlag", "westgard_flags", "PCAResult", "pca", "VolcanoPoint", "volcano", "benjamini_hochberg",
    "CalibrationSettings", "CalibrationFit", "CalibrationError", "fit_calibration",
    "QuantSettings", "QuantDecisions", "TargetHit", "RunQuant", "SequenceRun", "quantify_run",
    "learn_positions", "with_positions", "build_calibrations", "quantify_sequence",
    "all_signals", "istd_table", "cal_set_for", "bracket_cal_sets", "CalBracket", "blended_amount",
    "read_sequence_extras",
    "QuantReport", "finish_sequence", "quantify_prepared",
    "plot_chromatogram", "chromatogram_caption", "Provenance", "build_provenance", "write_feature_table_csv", "write_excel",
    "plot_qc_chart", "plot_pca", "plot_volcano", "save_figure", "write_html_report", "write_pdf_report",
    "OverlapSettings", "BlobCube", "blob_cube", "IonApex", "UniqueIonMap", "unique_ion_apexes",
    "split_by_line", "split_by_watershed", "split_by_valley", "split_volumes",
    "OverlapComponent", "OverlapFit", "OverlapCV", "ComponentAssessment",
    "fit_overlap", "fit_template_overlap", "predict_overlap", "component_maps", "assess_components",
    "cross_validate_overlap", "choose_k",
    "MCRResult", "MCRStability", "mcr_als", "mcr_restarts", "PARAFACResult", "parafac", "core_consistency", "split_half_stability",
    "build_run_tensor",
    "retention_index", "IdentifySettings", "LibraryEntry", "Library", "build_library",
    "load_msp", "write_msp", "build_inhouse_entry", "generate_decoys",
    "SearchHit", "search", "assign_qvalues", "classify", "Identification", "identify_blobs",
    "identify_components",
    "promote_to_target",
]
