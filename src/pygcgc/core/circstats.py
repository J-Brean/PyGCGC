"""Circular statistics for second-dimension times.

Second-dimension time wraps round every modulation period, like the hand
of a clock. With a 4 s period, 3.9 s and 0.1 s are only 0.2 s apart, but
their ordinary average (2.0 s) is nowhere near either. The fix is to treat
each time as an angle on a circle one period round, average the angles as
arrows, and turn the direction of the average arrow back into seconds.
The length of that average arrow (between 0 and 1) says how tightly the
values cluster: 1 means identical, near 0 means scattered.
"""

from __future__ import annotations

import numpy as np


def _to_angle(x, period):
    return 2.0 * np.pi * np.asarray(x, dtype=np.float64) / period


def resultant_length(x, period: float, weights=None) -> float:
    a = _to_angle(x, period)
    w = np.ones_like(a) if weights is None else np.asarray(weights, dtype=np.float64)
    c, s = np.sum(w * np.cos(a)), np.sum(w * np.sin(a))
    return float(np.hypot(c, s) / np.sum(w))


def circ_mean(x, period: float, weights=None) -> float:
    """Circular mean in [0, period)."""
    a = _to_angle(x, period)
    w = np.ones_like(a) if weights is None else np.asarray(weights, dtype=np.float64)
    ang = np.arctan2(np.sum(w * np.sin(a)), np.sum(w * np.cos(a)))
    return float((ang / (2.0 * np.pi) * period) % period)


def circ_std(x, period: float, weights=None) -> float:
    """Circular standard deviation, in the same units as x."""
    r = resultant_length(x, period, weights)
    r = min(max(r, 1e-12), 1.0)
    return float(np.sqrt(-2.0 * np.log(r)) * period / (2.0 * np.pi))


def circ_diff(a, b, period: float):
    """Signed shortest difference a - b, in (-period/2, period/2]."""
    d = (np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)) % period
    return np.where(d > period / 2.0, d - period, d)


def circ_median(x, period: float) -> float:
    """The data point with the smallest summed circular distance to the others."""
    x = np.asarray(x, dtype=np.float64) % period
    cost = np.abs(circ_diff(x[:, None], x[None, :], period)).sum(axis=1)
    return float(x[np.argmin(cost)])
