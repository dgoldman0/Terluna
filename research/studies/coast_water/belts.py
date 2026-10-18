"""The belts a fresh shore's water table can feed: steady sideways flow from the water's edge with the plants' uptake as
the sink.

Plants draw E (m/s) from the water table wherever it lies within the reach R (m) below the ground. Under a belt of
width W fed only from the water's edge, steady Dupuit flow T h'' = E with the head held at the water's level at the
shore and no flow at the belt's inland edge (beyond it nothing draws and the deep water table stands flat) puts the
water table c(2Wx - x^2) below the water's level at a distance x inland, c = E / (2T), and c W^2 below it at the edge.
On a shore whose ground rises steadily, at a slope s, the belt ends where the ground's height and that fall together
reach the reach:

    s W + c W^2 = R,   W = 2R / (s + sqrt(s^2 + 4cR)).

A steep shore sets W = R/s and a flat one the ground's supply, W = sqrt(R/c) = sqrt(2TR/E), geography's flat-shore
screen. Across a strip of land between two shores, or an island, the same holds from each shore to the divide. Where
the ground does not rise steadily, a point is within reach when its own height and the fall to it under uptake all the
way from the shore stay within the reach; that is exact on a steady rise and errs low behind a ridge, whose plants
cannot reach and so draw nothing.

The arithmetic on the 16 px/deg shores (research.studies.coast_water.shores): a pixel beside the water holds a strip
along each edge facing it, as wide as the belt on the slope from the shoreline, at the edge, to the pixel's own height
at its centre, half its depth inland, and no wider than the pixel; its edges are a staircase, so their length counts
pi/4, the mean for a shore of random direction, as the Jura Island screen does. A pixel further inland holds the share
of its own depth that the belt on the steady slope from the shoreline to its centre covers.
"""
from __future__ import annotations

import numpy as np

STAIRCASE = np.pi / 4.0       # a staircase edge's length against a shore of random direction


def width(slope, c, reach):
    """Belt width (m) on a shore rising at slope (m/m) from the water, with c = E/(2T) (1/m) and the reach (m)."""
    s = np.maximum(np.asarray(slope, float), 0.0)
    c = np.asarray(c, float)
    with np.errstate(divide='ignore'):                      # no uptake on flat ground: the belt has no limit
        return 2.0 * reach / (s + np.sqrt(s * s + 4.0 * c * reach))


def flat_width(transmissivity, uptake, reach):
    """The widest belt the ground can feed on a flat shore (m): sqrt(2TR/E)."""
    return np.sqrt(2.0 * transmissivity * reach / uptake)


def water_table(x, belt, c):
    """The water table's fall below the water's level (m) at x (m) inland of a belt of width belt (m)."""
    x = np.minimum(np.asarray(x, float), belt)
    return c * (2.0 * belt * x - x * x)


def first_ring(edge_km, across_km, rise_m, c, reach):
    """Belt area (km2) of pixels beside the water, by class in the columns: a strip along each edge facing the water,
    as wide as the belt on the slope from the shoreline to the pixel's height at half its depth, at most the depth.
    Returns the strips' areas, their widths (m) and their shore lengths (km) after the staircase factor."""
    edge = np.asarray(edge_km, float)
    across = np.asarray(across_km, float)
    rise = np.nan_to_num(np.asarray(rise_m, float), nan=0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        slope = np.where(across > 0, np.maximum(rise, 0.0) / (0.5 * across * 1000.0), 0.0)
    w = width(slope, c, reach)
    length = STAIRCASE * edge
    area = length * np.minimum(w / 1000.0, across)
    return np.where(edge > 0, area, 0.0), np.where(edge > 0, w, 0.0), length


def inland(distance_km, step_rows, step_cols, height_m, c, reach):
    """Share of each inland pixel's depth (from the shoreline) inside the belt. The shoreline lies half a pixel's
    spacing along the step to the nearest water pixel short of that pixel's centre."""
    d = np.asarray(distance_km, float)
    n = np.maximum(np.maximum(np.abs(step_rows), np.abs(step_cols)), 1).astype(float)
    spacing = d / n
    centre = d - 0.5 * spacing                                # the pixel centre's distance from the shoreline
    slope = np.maximum(np.asarray(height_m, float), 0.0) / (centre * 1000.0)
    w = width(slope, c, reach) / 1000.0
    return np.clip((w - (centre - 0.5 * spacing)) / spacing, 0.0, 1.0)
