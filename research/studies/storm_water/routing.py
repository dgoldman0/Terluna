"""Routing the storm water to the hollows (storm-water study, phase 2a), at any gravity.

A storm's runoff leaves the point as infiltration excess (phase 1 and wetting.py) and meets, in order:

1. The slope. Sheet flow down a hillslope of length L_h is a kinematic wave under Manning's law, q = (sqrt(g / g_E)
   sqrt(S) / n_E) h^(5/3): Manning's n carries sqrt(g), so an Earth n becomes 2.45 n on the Moon, and at the same
   discharge the flow is 1.71 times as deep and 0.58 times as fast. When the excess stops, the water still on the slope
   soaks into it at the soil's capacity as it drains (run-on infiltration); a slower, deeper sheet gives the slope
   longer to take it.
2. The crater bowls below the grid. Their numbers follow the terrain study's saturation equilibrium (1.54 D^-2 per km^2
   with D in km, 0.69-3.9% of geometric saturation): water falling inside a bowl stays there, and a flow path that
   crosses a bowl fills it before going on, so each crossing takes up to the bowl's volume. Young land keeps bowls from
   10 m up to the grid's 1.9 km; mature land keeps those of 100 m and more.
3. The channels. Each 1.9 km cell's flow passes down a reach to its steepest neighbour, losing water to the dry bed by
   Lane's (1983) method (NEH-4 chapter 19), recast for gravity: Lane's K D, the depth a dry bed of effective
   conductivity K takes in during a flow of duration D, becomes the Green and Ampt depth the same bed takes at the
   gravity in question over the Earth's, so that the capillary soaking stays Earth's and the drainage falls to a sixth.
4. The spread of a storm. A storm makes runoff over a footprint A_s; where a reach's catchment is larger, the flow
   reaching it in one event comes from a footprint's worth of it, and the reach sees proportionally more events.

Water ends in a closed depression (any cell below its depression's spill level), the sea or the atlas's lakes, the
humid land beyond the drylands, the sub-grid bowls, or the soil of slopes and channel beds.
"""
from __future__ import annotations

import math

import numpy as np

from research.studies.storm_water import infiltration as IN

ACRE_FT_M3 = 1233.48184
MILE_KM = 1.609344
FOOT_M = 0.3048
INCH_MM = 25.4

# Design choices (each named as such in the README) and Earth values with their sources.
HILLSLOPE = dict(length_m=100.0, slope=0.05, manning_n_earth=0.05)
FOOTPRINT_KM2 = 5.0                # the runoff-making core of a storm: Walnut Gulch's storms mostly cover its 148 km2
                                   # and their cores predict the runoff (Syed et al. 2003, abstract); no core area was
                                   # found, and the basin check passes for 3-9 km2, so it sets this value
WIDTH_M_AT_KM2 = 2.0               # design: channel width w = 2.0 m (A / 1 km2)^0.5
SINUOSITY = 1.2                    # design
WAVE_SPEED_M_S_EARTH = 1.5         # design: a flood wave's speed down Earth's dry sandy channels; on the Moon times the
                                   # Manning factor (0.58 at the same discharge)
# Lane's Table 19-1 (NEH-4 ch. 19, read at the source): class 2, clean sand and gravel in field conditions, 2-5 in/h.
BED = dict(k_eff_in_h=3.0, suction_mm=49.5, theta_e=0.417, se_dry=0.03)   # sand's Green-Ampt suction (Rawls et al. 1983)
# Bowls: the terrain study's equilibrium law and its fresh depth ratio (Pike 1977), equilibrium bowls at half of it.
GEOMETRIC_SATURATION = 1.54
EQUILIBRIUM_SHARE = (0.0069, 0.039)
BOWL_DEPTH_RATIO = 0.5 * 0.196     # design: equilibrium craters are shallower than fresh ones
BOWL_RANGES_M = dict(young=(10.0, 1900.0), mature=(100.0, 1900.0))


def manning_alpha(gravity, slope=HILLSLOPE['slope'], n_earth=HILLSLOPE['manning_n_earth']):
    """alpha in q = alpha h^(5/3) (m^(1/3)/s): Manning's law with its gravity explicit."""
    return math.sqrt(gravity) * math.sqrt(slope) / n_earth


def hillslope_loss(excess_mm, duration_h, capacity_mm_h, gravity, length_m=HILLSLOPE['length_m'], **kw):
    """Share of a storm's point runoff the slope takes back as it drains: the sheet's storage at the end of the excess
    (kinematic equilibrium or what has built up) soaks in at the capacity while the slope drains, over about half the
    equilibrium time. Arrays alike."""
    alpha = manning_alpha(gravity, **kw)
    e = np.maximum(np.asarray(excess_mm, float), 1e-9) / 1000.0 / np.maximum(duration_h, 1e-6) / 3600.0      # m/s
    m = 5.0 / 3.0
    t_eq = (length_m / (alpha * e ** (m - 1.0))) ** (1.0 / m)                                                # s
    stored = np.minimum(e * duration_h * 3600.0 * length_m, length_m * m / (m + 1.0) * (e * length_m / alpha) ** (1.0 / m))
    soaked = np.minimum(stored, 0.5 * np.asarray(capacity_mm_h, float) / 1000.0 / 3600.0 * length_m * t_eq)
    return soaked / (e * duration_h * 3600.0 * length_m), t_eq / 3600.0


def bowl_classes(state, share=None, bins=16):
    """The bowls below the grid: per diameter bin, the number per km2, the diameter (km) and the volume (m3)."""
    lo, hi = BOWL_RANGES_M[state]
    s = math.sqrt(EQUILIBRIUM_SHARE[0] * EQUILIBRIUM_SHARE[1]) if share is None else share
    c = s * GEOMETRIC_SATURATION
    edges = np.geomspace(lo / 1000.0, hi / 1000.0, bins + 1)
    number = c * (edges[:-1] ** -2 - edges[1:] ** -2)
    d = np.sqrt(edges[:-1] * edges[1:])
    volume = math.pi / 8.0 * (d * 1000.0) ** 2 * BOWL_DEPTH_RATIO * d * 1000.0
    return dict(number_km2=number, diameter_km=d, volume_m3=volume)


def bowl_cover(b):
    """Share of the land inside bowls (overlapping as a Poisson field)."""
    return 1.0 - math.exp(-float((b['number_km2'] * math.pi * b['diameter_km'] ** 2 / 4.0).sum()))


def bowl_survival(b, path_km, volume_m3):
    """Share of a flow of volume_m3 (per event) passing path_km of a bowl field: each crossing (a bowl's centre within
    half its diameter of the path) takes up to its volume."""
    rate = b['number_km2'] * b['diameter_km']                                            # crossings per km
    v = np.maximum(np.asarray(volume_m3, float), 1e-9)[..., None]
    return np.exp(-(rate * np.asarray(path_km, float)[..., None] * np.minimum(1.0, b['volume_m3'] / v)).sum(axis=-1))


def ponded_infiltration(hours, k_mm_h, suction_mm, dtheta):
    """Green and Ampt's cumulative infiltration (mm) under ponding after `hours`: t = (F - c ln(1 + F/c)) / K."""
    t = np.asarray(hours, float)
    c = suction_mm * dtheta
    f = np.maximum(np.sqrt(2.0 * c * k_mm_h * t), k_mm_h * t)
    for _ in range(60):
        g = f - c * np.log1p(f / c) - k_mm_h * t
        f = np.maximum(f - g / (f / (c + f)), 1e-9)
    return f


def lane_kd_inches(duration_h, gravity, bed=BED):
    """Lane's K D (inches) recast for gravity: the Earth bed's K D times the ratio of the Green and Ampt depths a dry bed
    takes in at the gravity and at Earth's (conductivity times g, suction over g)."""
    k = bed['k_eff_in_h'] * INCH_MM
    dth = bed['theta_e'] * (1.0 - bed['se_dry'])
    earth = ponded_infiltration(duration_h, k, bed['suction_mm'], dth)
    here = ponded_infiltration(duration_h, k * gravity, bed['suction_mm'] / gravity, dth)
    return bed['k_eff_in_h'] * np.asarray(duration_h, float) * here / earth


def wave_speed_m_s(gravity):
    """A flood wave's speed: Earth's design value times Manning's factor at the same discharge, g^(3/10)."""
    return WAVE_SPEED_M_S_EARTH * gravity ** 0.3


def kd_table(gravity, bed=BED):
    """Lane's recast K D (inches) against the flow's duration (h), for interpolation in log duration."""
    d = np.geomspace(0.01, 2000.0, 200)
    return np.log(d), lane_kd_inches(d, gravity, bed)


def flow_duration_h(d0_h, upstream_path_km, gravity):
    """A flow's duration at a reach: its duration where it formed plus its spread over the travel from the farthest
    source, the path over the wave's speed."""
    return d0_h + np.asarray(upstream_path_km, float) * 1000.0 / wave_speed_m_s(gravity) / 3600.0


def lane_outflow(volume_m3, length_km, width_m, kd_in):
    """Lane's outflow volume (m3) of a reach for an inflow volume (eqs. 19-1, 19-13 to 19-16, in their units)."""
    p = np.maximum(np.asarray(volume_m3, float), 0.0) / ACRE_FT_M3
    x = np.asarray(length_km, float) / MILE_KM
    w = np.asarray(width_m, float) / FOOT_M
    a_unit = -0.00465 * kd_in
    arg = 1.0 - 0.00545 * kd_in / np.maximum(p, 1e-30)
    with np.errstate(divide='ignore', invalid='ignore'):
        k = -1.09 * np.log(np.where(arg > 0.0, arg, 1.0))
        b_unit = np.exp(-k)
        b_xw = np.exp(-k * x * w)
        a_xw = np.where(k > 0.0, a_unit / (1.0 - b_unit) * (1.0 - b_xw), a_unit * x * w)
        q = np.where((arg > 0.0) & (p > -a_xw / np.maximum(b_xw, 1e-300)), a_xw + b_xw * p, 0.0)
    return np.maximum(q, 0.0) * ACRE_FT_M3


def channel_width_m(area_km2):
    return WIDTH_M_AT_KM2 * np.sqrt(np.maximum(area_km2, 1e-3))


def route_events(order, nxt, local, area_m2, length_km, upstream_km2, d0_h, gravity, bowls, footprint_km2=None):
    """Route one event per column through a forest of reaches in topological order.

    order: list of index arrays (fronts; donors before receivers); nxt: receiver index, -1 where the flow leaves the
    forest; local: (cells, columns) m3 of water entering each reach per event; area_m2: each reach's own area;
    length_km, upstream_km2: reach length and catchment; d0_h: each column's flow duration where it forms; gravity in
    Earth's; bowls: per column a bowl field or None. Returns per reach the channel loss, the bowl capture and the outflow (m3 per point-event, the
    events' multiplicity included), the water leaving the forest at each reach, and the per-event inflow volume."""
    n, k = local.shape
    total = local.astype(float).copy()
    area = np.repeat(area_m2[:, None], k, axis=1).astype(float)
    lost, caught, out, event = (np.zeros((n, k)) for _ in range(4))
    width = channel_width_m(upstream_km2)
    fp = (FOOTPRINT_KM2 if footprint_km2 is None else footprint_km2) * 1e6
    logd, kd = kd_table(gravity)
    path = np.zeros(n)
    for f in order:
        t, a = total[f], area[f]
        mult = np.maximum(1.0, a / fp)
        v = t / mult
        dur = flow_duration_h(np.asarray(d0_h, float)[None, :], (path[f] + 0.5 * length_km[f])[:, None], gravity)
        kd_in = np.interp(np.log(dur), logd, kd)
        q = lane_outflow(v, length_km[f][:, None], width[f][:, None], kd_in)
        s = np.ones_like(q)
        for j, b in enumerate(bowls):
            if b is not None:
                s[:, j] = bowl_survival(b, length_km[f], q[:, j])
        lost[f] = mult * (v - q)
        caught[f] = mult * q * (1.0 - s)
        o = mult * q * s
        out[f] = o
        event[f] = v
        r = nxt[f]
        go = r >= 0
        if go.any():
            with np.errstate(divide='ignore', invalid='ignore'):
                live = np.where(t > 0.0, o / t, 0.0)
            np.add.at(total, r[go], o[go])
            np.add.at(area, r[go], (a * live)[go])
            np.maximum.at(path, r[go], (path[f] + length_km[f])[go])
    return lost, caught, out, event


def fronts(nxt):
    """Topological fronts of a forest given each node's receiver (-1 at the forest's edge)."""
    n = nxt.size
    has = nxt >= 0
    indeg = np.bincount(nxt[has], minlength=n).astype(np.int32)
    front = np.flatnonzero(indeg == 0)
    order = []
    while front.size:
        order.append(front)
        r = nxt[front]
        r = r[r >= 0]
        np.subtract.at(indeg, r, 1)
        cand = np.unique(r)
        front = cand[indeg[cand] == 0]
    return order


def event_classes(p, lib, se, classes=4):
    """The library's storms on fresh soil at saturation se, grouped by their runoff into classes of about equal runoff:
    each class's probability per storm, mean runoff (mm), runoff duration (h) and the soil's capacity when it ends."""
    dt = lib['step_min'] / 60.0
    k = len(lib['storms'])
    n = max(len(s) for s in lib['storms'])
    rates = np.zeros((n, k))
    for j, s in enumerate(lib['storms']):
        rates[: len(s), j] = s
    dth = IN.deficit(p, np.full(k, se))
    st = IN.State((k,))
    minutes = np.zeros(k)
    for t in range(n):
        before = st.runoff.copy()
        IN.step(p, st, dth, rates[t], dt)
        minutes += (st.runoff - before) > 1e-6
    IN.end_storm(p, st)
    off = st.runoff
    cap = IN.capacity(p, np.maximum(st.f, 1e-6), dth, st.energy) if p['forming'] else IN.capacity(p, np.maximum(st.f, 1e-6), dth)
    w = np.asarray(lib['weights'], float) / np.sum(lib['weights'])
    order = np.argsort(off)
    run = order[off[order] > 0.05]
    out = []
    if run.size == 0:
        return out
    cum = np.cumsum(off[run] * w[run])
    edges = np.searchsorted(cum / cum[-1], np.linspace(0.0, 1.0, classes + 1)[1:-1], side='right')
    for g in np.split(run, edges):
        if g.size == 0:
            continue
        ww = w[g]
        out.append(dict(probability=float(ww.sum()), runoff_mm=float((off[g] * ww).sum() / ww.sum()),
                        duration_h=float(max((minutes[g] * ww).sum() / ww.sum() * dt, dt)),
                        capacity_mm_h=float((np.minimum(cap[g], 1e3) * ww).sum() / ww.sum()), storms=int(g.size)))
    return out
