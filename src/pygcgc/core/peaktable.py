"""Fitted components as the peak list that the downstream modules read.

Classification, clustering and identification were written for legacy blobs
(`peaks.BlobSet`): one row per peak with a t1, t2, volume and a spectrum over the
run's whole m/z axis. A fitted component carries a spectrum only over the ~30 ions
its calculation window used, and a molecular ion is easily not among them. This
module therefore re-solves every component's spectrum on *all* ions with its
fitted shape held fixed (a linear solve, no new peak search), and offers the
result through the same attribute names as a `BlobSet`.

It also accounts for the signal: how much of the run's total ion current the
fitted components explain, where the rest is, and how that depends on intensity.
Both need the cleaned run, so they run on demand, not inside peak finding.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage as ndi

from .candidates import CandidateRegion, _fold_region
from .components import ComponentSet, component_profiles_from_durable
from .peaks import PeakSettings, make_frame


@dataclass
class ComponentPeaks:
    """A `BlobSet`-shaped view of fitted components (only what classify, cluster and identify read)."""
    t1: np.ndarray
    t2: np.ndarray
    volume: np.ndarray
    mz: np.ndarray
    spectra: np.ndarray                 # (n, len(mz)), fitted on every ion
    co_elution: np.ndarray              # True when the component shares its fit with another or is unresolved
    confidence: list = field(default_factory=list)
    component_id: list = field(default_factory=list)
    unresolved: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=bool))
    apex_spectra: np.ndarray = field(default_factory=lambda: np.zeros((0, 0), dtype=np.float32))
    contours: list = field(default_factory=list)   # (t1, t2) outline of each component's 10% contour

    def __len__(self) -> int:
        return int(self.t1.size)

    def outline_times(self, i: int) -> tuple:
        """Outline for drawing, in the same (start, end, NaN) format as BlobSet.outline_times."""
        return self.contours[i]


def respectrum_components(run, fold, component_set: ComponentSet, wrap_below_s: float = 0.5,
                          background_order: int = 1, progress=None, cancel=None
                          ) -> tuple[ComponentPeaks, np.ndarray, np.ndarray]:
    """Full-ion spectra and volumes for every component, with fitted shapes fixed.

    Returns ``(peaks, model_tic, footprint)``. ``model_tic`` is the summed-ion image
    of all components in the extended frame; ``footprint`` marks pixels within three
    fitted widths of a component. Per window: fold every ion into the window, solve
    all ions at once with a shared pixel weight, keep the peak coefficients
    non-negative and let a bilinear background follow the signed residual (the same
    approximation the rapid shape fit uses, so amounts are comparable with it).
    """
    frame = make_frame(fold, PeakSettings(wrap_below_s=wrap_below_s))
    shape = (frame.rows, len(frame.mid_t1))
    model_tic = np.zeros(shape)
    footprint = np.zeros(shape, dtype=bool)
    every_ion = np.arange(run.n_mz)
    rows_out = {"t1": [], "t2": [], "volume": [], "co": [], "conf": [], "id": [], "unres": []}
    spectra, contours = [], []
    groups = component_set.groups
    for done, group in enumerate(groups, start=1):
        if cancel is not None:
            cancel()
        region = CandidateRegion((), group.row0, group.row0 + group.n_rows,
                                 group.col0, group.col0 + group.n_cols)
        data = _fold_region(run, fold, every_ion, region, frame.wrap_rows).astype(np.float64)
        n_rows, n_cols, n_ions = data.shape
        profiles = np.stack([component_profiles_from_durable(c, (n_rows, n_cols))
                             for c in group.components], axis=2) if group.components else None
        if profiles is None:
            continue
        flat = data.reshape(-1, n_ions)
        weight = 1.0 / np.sqrt(np.clip(flat.mean(axis=1), 1.0, None))
        peak_design = profiles.reshape(-1, profiles.shape[2])
        rr, cc = np.indices((n_rows, n_cols), dtype=float)
        rn, cn = rr.ravel() / max(n_rows - 1, 1), cc.ravel() / max(n_cols - 1, 1)
        background = (np.column_stack(((1 - rn) * (1 - cn), (1 - rn) * cn, rn * (1 - cn), rn * cn))
                      if background_order == 1 else np.ones((rn.size, 1)))
        k = peak_design.shape[1]
        coef, *_ = np.linalg.lstsq(np.column_stack((peak_design, background)) * weight[:, None],
                                   flat * weight[:, None], rcond=None)
        peak_coef = np.maximum(coef[:k], 0.0)
        for j, component in enumerate(group.components):
            profile = profiles[:, :, j]
            spectrum = peak_coef[j]
            total = float(spectrum.sum())
            spectra.append(spectrum.astype(np.float32))
            rows_out["t1"].append(component.t1)
            rows_out["t2"].append(component.t2)
            rows_out["volume"].append(float(profile.sum()) * total)
            rows_out["co"].append(len(group.components) > 1 or group.unresolved)
            rows_out["conf"].append(component.confidence)
            rows_out["id"].append(component.id)
            rows_out["unres"].append(group.unresolved)
            contours.append(group.contour(j, 0.10))
            r0, c0 = group.row0, group.col0
            model_tic[r0:r0 + n_rows, c0:c0 + n_cols] += profile * total
            near = (((rr - component.row0) / max(component.sigma_row, 1e-6)) ** 2
                    + ((cc - component.col0) / max(component.sigma_col, 1e-6)) ** 2) <= 9.0
            footprint[r0:r0 + n_rows, c0:c0 + n_cols] |= near
        if progress is not None:
            progress(done, len(groups))
    peaks = ComponentPeaks(
        t1=np.asarray(rows_out["t1"], float), t2=np.asarray(rows_out["t2"], float),
        volume=np.asarray(rows_out["volume"], float), mz=np.asarray(run.mz),
        spectra=np.vstack(spectra) if spectra else np.zeros((0, run.n_mz), np.float32),
        co_elution=np.asarray(rows_out["co"], dtype=bool), confidence=rows_out["conf"],
        component_id=rows_out["id"], unresolved=np.asarray(rows_out["unres"], dtype=bool),
        apex_spectra=np.zeros((0, run.n_mz), np.float32), contours=contours)
    return peaks, model_tic, footprint


def signal_accounting(tic: np.ndarray, model_tic: np.ndarray, footprint: np.ndarray,
                      pitch: float, period: float, t1_start: float, top: int = 10) -> dict:
    """How much of the run's signal the fitted peaks explain, and where the rest is.

    ``tic`` is the cleaned run's folded total-ion image (rows x modulations);
    ``model_tic`` and ``footprint`` are extended-frame images from
    :func:`respectrum_components`. Rows beyond one modulation hold the wrapped
    part of a peak and are folded back onto the following modulation first.
    """
    n_rows, n_cols = tic.shape
    model = model_tic[:n_rows].copy()
    mask = footprint[:n_rows].copy()
    extra = model_tic.shape[0] - n_rows
    if extra > 0:
        model[:extra, 1:] += model_tic[n_rows:n_rows + extra, :-1]
        mask[:extra, 1:] |= footprint[n_rows:n_rows + extra, :-1]
    positive = np.clip(tic, 0.0, None)
    residual = tic - model
    total = float(positive.sum())
    unexplained = np.clip(residual, 0.0, None)
    # Split what is left into a broad, smooth part (an unresolved hump or baseline
    # residue, wider than any peak) and discrete peak-like leftovers. Only the latter
    # are candidates for peaks the finder missed.
    smooth = np.minimum(ndi.gaussian_filter(unexplained, (12.0, 3.0)), unexplained)
    peaklike = unexplained - smooth
    centre = float(np.median(residual))
    noise = 1.4826 * float(np.median(np.abs(residual - centre)))
    patches, n_patch = ndi.label(ndi.gaussian_filter(peaklike, (2.0, 0.5)) > 4.0 * noise / 3.0)
    listed = []
    for label in range(1, n_patch + 1):
        rows, cols = np.nonzero(patches == label)
        weights = np.clip(peaklike[rows, cols], 1e-9, None)
        listed.append((float(peaklike[rows, cols].sum()), float(np.average(cols, weights=weights)),
                       float(np.average(rows, weights=weights)), int(mask[rows, cols].mean() > 0.5)))
    listed.sort(reverse=True)
    decades = []
    if total > 0:
        edges = 10.0 ** np.arange(np.floor(np.log10(max(np.percentile(positive[positive > 0], 1), 1.0))),
                                  np.ceil(np.log10(positive.max())) + 1)
        for lo, hi in zip(edges[:-1], edges[1:]):
            here = (positive >= lo) & (positive < hi)
            share = float(positive[here].sum())
            if share > 0.001 * total:
                decades.append({"pixel_tic_from": float(lo), "pixel_tic_to": float(hi),
                                "percent_of_signal": 100.0 * share / total,
                                "percent_of_that_explained": 100.0 * float(
                                    np.minimum(np.clip(model, 0, None)[here], positive[here]).sum())
                                / max(share, 1e-12)})
    explained = float(np.minimum(np.clip(model, 0.0, None), positive).sum())
    return {
        "total_signal": total,
        "percent_explained": 100.0 * explained / total if total > 0 else float("nan"),
        "percent_unexplained_outside_all_peaks": (
            100.0 * float(unexplained[~mask].sum()) / total if total > 0 else float("nan")),
        "percent_unexplained_inside_peak_footprints": (
            100.0 * float(unexplained[mask].sum()) / total if total > 0 else float("nan")),
        "percent_unexplained_broad_smooth": (
            100.0 * float(smooth.sum()) / total if total > 0 else float("nan")),
        "percent_unexplained_peaklike": (
            100.0 * float(peaklike.sum()) / total if total > 0 else float("nan")),
        "percent_overfitted": (
            100.0 * float(np.clip(model - np.clip(tic, 0, None), 0.0, None).sum()) / total
            if total > 0 else float("nan")),
        "by_pixel_intensity": decades,
        "largest_unexplained_patches": [
            {"t1_s": t1_start + col * period + 0.5 * period, "t2_s": row * pitch,
             "volume": volume, "percent_of_signal": 100.0 * volume / max(total, 1e-12),
             "inside_a_fitted_footprint": bool(inside)}
            for volume, col, row, inside in listed[:top]],
    }


def format_accounting(report: dict) -> str:
    """The few lines worth showing in a status bar."""
    if not np.isfinite(report["percent_explained"]):
        return "No signal to account for."
    line = (f"Fitted peaks explain {report['percent_explained']:.0f}% of the total signal; "
            f"{report['percent_unexplained_outside_all_peaks']:.0f}% lies outside every fitted "
            f"peak and {report['percent_unexplained_inside_peak_footprints']:.0f}% is unexplained "
            "inside a peak footprint")
    line += (f". Of the unexplained signal, {report['percent_unexplained_broad_smooth']:.0f}% of the total is a "
             f"broad smooth hump (not peaks) and {report['percent_unexplained_peaklike']:.0f}% is discrete "
             "peak-like structure the fit missed")
    low = [d for d in report["by_pixel_intensity"] if d["percent_of_signal"] >= 1.0]
    if low:
        lowest = min(low, key=lambda d: d["pixel_tic_from"])
        line += (f". The faintest band ({lowest['pixel_tic_from']:.0e} to {lowest['pixel_tic_to']:.0e} "
                 f"per pixel) holds {lowest['percent_of_signal']:.0f}% of the signal, "
                 f"{lowest['percent_of_that_explained']:.0f}% of it explained")
    missed = [p for p in report["largest_unexplained_patches"]
              if not p["inside_a_fitted_footprint"]][:3]
    if missed:
        line += ". Largest signal outside every fitted peak: " + "; ".join(
            f"t1 {p['t1_s'] / 60.0:.2f} min, t2 {p['t2_s']:.2f} s ({p['percent_of_signal']:.1f}%)"
            for p in missed)
    return line + "."
