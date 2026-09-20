"""Finding a solvent front.

With a thermal modulator the solvent arrives as one huge, broad event early
in the run, tens of seconds long, followed by a long tail. It is far brighter
than anything else, it fills every second-dimension row of the columns it
covers, and it is not tied to the modulation in any useful way. Left in, it
swamps both the period and the phase estimates.

The front is the first stretch, early in the active part of the run, where
the smoothed TIC stays far above its typical level. The tail is followed
until the signal drops back near typical. If nothing stands out that
strongly (most runs without a solvent front) nothing is returned, so a
bright sample region is never mistaken for one.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter1d


def find_solvent_front(times: np.ndarray, tic: np.ndarray, factor: float = 20.0,
                       tail_factor: float = 3.0, smooth_s: float = 10.0,
                       merge_s: float = 30.0, max_start_fraction: float = 0.5,
                       max_span_fraction: float = 0.4) -> tuple[float, float] | None:
    """(start, end) of the solvent front in seconds, or None."""
    times = np.asarray(times, dtype=np.float64)
    tic = np.asarray(tic, dtype=np.float64)
    dt = float(np.median(np.diff(times)))
    width = max(3, int(round(smooth_s / dt)))
    smooth = uniform_filter1d(tic, width, mode="nearest")

    active = np.flatnonzero(smooth > 0)
    if active.size < 10 * width:
        return None
    first, last = active[0], active[-1]
    base = float(np.median(smooth[active]))
    if base <= 0:
        return None

    hot = np.flatnonzero(smooth > factor * base)
    if hot.size == 0 or hot[0] > first + max_start_fraction * (last - first):
        return None

    # Merge hot stretches separated by short dips into one event.
    gap = int(round(merge_s / dt))
    breaks = np.flatnonzero(np.diff(hot) > gap)
    end_hot = hot[breaks[0]] if breaks.size else hot[-1]

    # Follow the tail until the signal is back near typical.
    below = np.flatnonzero(smooth[end_hot:] < tail_factor * base)
    end = end_hot + (below[0] if below.size else last - end_hot)

    start_t, end_t = float(times[hot[0]]), float(times[min(end, times.size - 1)])
    if (end_t - start_t) > max_span_fraction * (times[last] - times[first]):
        return None
    return start_t, end_t
