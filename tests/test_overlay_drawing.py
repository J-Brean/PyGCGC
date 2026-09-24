"""Overlay outlines are simplified for speed without changing what is drawn."""

import numpy as np

from pygcgc.gui.views.image2d import _merge_collinear


def _triples(segments):
    x = np.array([[a, b, np.nan] for (a, _c), (b, _d) in segments]).ravel()
    y = np.array([[c, d, np.nan] for (_a, c), (_b, d) in segments]).ravel()
    return x, y


def _covered(x, y, step=0.5):
    """Every half-unit point on every segment, as a set: what the drawing covers."""
    points = set()
    for x0, x1, y0, y1 in zip(x[0::3], x[1::3], y[0::3], y[1::3]):
        n = int(round(max(abs(x1 - x0), abs(y1 - y0)) / step))
        for k in range(n + 1):
            f = k / max(n, 1)
            points.add((round(x0 + f * (x1 - x0), 6), round(y0 + f * (y1 - y0), 6)))
    return points


def test_pixel_edges_are_merged_into_long_segments_without_changing_the_drawing():
    # A 3 x 4 pixel footprint's border, one segment per pixel edge, plus a duplicate
    # shared edge, a diagonal contour piece and a separate blob on the same column line.
    edges = []
    for r in range(4):
        edges += [((0, r), (0, r + 1)), ((3, r), (3, r + 1))]
    for c in range(3):
        edges += [((c, 0), (c + 1, 0)), ((c, 4), (c + 1, 4))]
    edges += [((0, 1), (0, 2)),                  # duplicate of a left edge
              ((5, 5), (6, 7)),                  # diagonal: passes through untouched
              ((0, 6), (0, 7))]                  # same line as the left side, but not touching
    x, y = _triples(edges)
    mx, my = _merge_collinear(x, y)
    assert len(mx) % 3 == 0 and np.isnan(mx[2::3]).all()
    assert _covered(mx, my) == _covered(x, y)
    assert len(mx) // 3 == 6            # left, right, top, bottom, diagonal, separate piece


def test_empty_outline():
    x, y = _merge_collinear(np.array([]), np.array([]))
    assert x.size == 0 and y.size == 0
