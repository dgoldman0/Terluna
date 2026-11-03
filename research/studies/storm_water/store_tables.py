"""The storm-water study's tables for the hollows' store in PlaSim's land scheme (stage (b) of the plan, 3 November 2026;
session-tools/agents/storm_water/stage_b_design.md and stage_b2_report.md on the research drive).

    climate/gcm/.venv/bin/python -m research.studies.storm_water.store_tables     # a few minutes, after run.py

PlaSim's land scheme with the store (climate/gcm/land_store.py on domain/climate, which reads this product as a file)
holds one lumped store of the hollows' water per T21 land cell. For every dryland cell this writes
results/store_tables.json (schema terluna.research.storm-water-store-tables/2):

- per storm library, surface and hollow age, the store's capacity and, at FILLS shares of it, the water area and the
  seepage under each floor design, with the cell's water shared among its hollows by their own yearly inflow (phase
  2b's pools, build_pools), each hollow at most full, and summed with phase 2b's own geometry and seepage laws:
  geography's closed depressions on the 1.9 km grid with their area rising linearly between cells, and the crater
  bowls below the grid as paraboloids at half Pike's fresh depth, wall-fed and channel-fed; the clay seal everywhere,
  the seal below a 1 m stage with sandy banks above, and the seal ten times tighter (phase 2b's sensitivity case);
- for the same cases and each floor, the balanced volume, at which the store's mean inflow meets its mean losses, as
  phase 2b starts its pools (balance_start), which PlaSim's store takes when its restart holds no store;
- per storm library, surface and hollow age, the shares of the slopes' runoff that reach the hollows, the rivers (the
  sea and the humid land) and the ground below the roots (the channel beds' and slopes' recharge), from phase 2a's
  routing; the rest soaks back into the slopes' soil and the channel beds above the roots;
- per storm library and surface, the slopes' runoff per millimetre of convective rain against a 3-day mean of the
  cell's convective rain and the bucket's relative wetness: phase 2a's own wetting tables (the library's storms at
  uniformly random times within a mean, the soil redistributing and drying between them) mixed over the Poisson number
  of storms the mean's rain brings; and per surface the runoff per millimetre of large-scale rain against its 3-day
  mean, the mean's large-scale rain falling as one steady rain of 2 mm/h, as phases 1 and 2a read it;
- per cell, the record's mean convective and large-scale rain, the start of PlaSim's running means of each.

The capacity of a depression with more cells below its pass than phase 2b keeps (14,470 of the 94,858 whose floors lie in
the drylands) is its volume at the pass; phase 2b gives those depressions walls 1 km tall above their sixteenth cell,
which never spill, and holds their area above it at the sixteen cells' (phase2b.depression_curves), which the curve
here keeps.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

from research.studies.storm_water import extract as EX
from research.studies.storm_water import infiltration as IN
from research.studies.storm_water import phase2a as P2
from research.studies.storm_water import phase2b as B
from research.studies.storm_water import routing as RT
from research.studies.storm_water import run as R1
from research.studies.storm_water import storms as ST
from research.studies.storm_water import wetting as W
from shared.constants import JULIAN_YEAR_DAYS, SYNODIC_MONTH_DAYS
from shared.provenance import constants_used

OUT = R1.HERE / 'results' / 'store_tables.json'
SCHEMA = 'terluna.research.storm-water-store-tables/2'
CODE = ('research/studies/storm_water/store_tables.py', 'research/studies/storm_water/phase2b.py',
        'research/studies/storm_water/phase2a.py', 'research/studies/storm_water/routing.py',
        'research/studies/storm_water/infiltration.py', 'research/studies/storm_water/wetting.py',
        'research/studies/storm_water/run.py', 'research/studies/storm_water/extract.py',
        'research/studies/storm_water/storms.py')
FILLS = np.array([0.0, 1e-6, 1e-5, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1.0])     # shares of capacity
RAIN_MM_DAY = np.array([0.0, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0,
                        13.0, 16.0, 20.0, 25.0, 32.0, 40.0, 50.0, 65.0, 80.0, 100.0, 130.0])   # a 3-day mean's rain
FILL_IN_WETNESS = (0.0, 0.01, 0.02, 0.04, 1.0)   # with the soils' own breakpoints (wetness_grid)
FLOORS = ('clay_floored', 'sandy_banked', 'tight_seal')
TIGHT = dict(B.SEAL, k_earth_mm_h=B.SEAL['k_earth_mm_h'] / 10.0)    # phase 2b's sensitivity seal, ten times tighter
LIBRARIES, SURFACES, STATES = R1.BRACKETS, IN.MAIN, B.STATES
DEPTH_TAB = np.linspace(0.0, 30.0, 601)          # phase 2b's seepage table (simulate's depth_tab)
DIGITS = 7                                       # significant digits the product keeps of each table entry


def wetness_grid():
    """The bucket's relative wetness at which the runoff tables are given: FILL_IN_WETNESS and every wetness at which
    a surface's soil reaches one of phase 2a's saturations (SAT_W), so that the convective tables, linear in the
    saturation between those and the saturation linear in the bucket's water, are exact between the grid's points."""
    pts = set(FILL_IN_WETNESS)
    for soil in sorted({IN.SURFACES[n]['soil'] for n in SURFACES}):
        s0, s1 = IN.initial_saturation(IN.SOILS[soil], np.array([0.0, IN.WSMAX_M]))
        for s in P2.SAT_W:
            if s0 < s < s1:
                pts.add(round(float((s - s0) / (s1 - s0)), 12))
    return np.array(sorted(pts))


WETNESS = wetness_grid()


def geometry(v, s, A, V):
    """Stage and area of each hollow at volume v (phase 2b's geometry closure in simulate): the first segment's law
    A = g1 stage (a wedge or a paraboloid, area as the root of the volume), then the piecewise curve (on_curve)."""
    g1 = A[:, 1] / np.maximum(s[:, 1], 1e-9)
    stage = np.sqrt(2.0 * v / np.maximum(g1, 1e-12))
    area = g1 * stage
    big = v > V[:, 1]
    if big.any():
        stage[big], area[big] = B.on_curve(v[big], s[big], A[big], V[big])
    return stage, area


def seepage(dmean, area, clog, floor):
    """Seepage (m3/s) of each hollow with mean depth dmean (m) over its water area (m2) under a floor design, by phase
    2b's laws (simulate's seepage_rate): the seal under the mean depth over the whole area, or with sandy banks the
    seal below the clogged stage's area clog and the banks' sand above it; the tight seal as phase 2b's sensitivity."""
    if floor == 'tight_seal':
        return np.interp(dmean, DEPTH_TAB, B.seal_seepage(DEPTH_TAB, IN.MOON, seal=TIGHT)) * area
    q_seal = np.interp(dmean, DEPTH_TAB, B.seal_seepage(DEPTH_TAB, IN.MOON)) * area
    if floor == 'clay_floored':
        return q_seal
    q_bank = B.bank_seepage(IN.MOON)
    return np.where(area > clog, q_seal / np.maximum(area, 1e-30) * clog + q_bank * (area - clog), q_seal)


def depression_capacity(inp):
    """The capacity (m3) of each dryland depression in phase 2b's order (its inputs' depressions with floors in the
    drylands): its top where at most sixteen cells lie below its pass, else its volume at the pass."""
    keep = inp['floor_t21'] >= 0
    nb = inp['n_below'][keep]
    s, A, V, top = B.depression_curves(inp['hyps_z'][keep].astype(float), inp['hyps_a'][keep].astype(float), nb,
                                       inp['pass_m'][keep].astype(float), inp['vmax_m3'][keep])
    full = nb <= s.shape[1] - 1
    cap = np.where(full, top, np.maximum(inp['vmax_m3'][keep], V[np.arange(nb.size), np.minimum(nb, s.shape[1] - 1)]))
    return cap, ~full


def hollows(inp, b, name, state, cap_d):
    """Phase 2b's pools of one library, surface and hollow age (build_pools): the grid's depressions (their capacity
    cap_d, depression_capacity) and per cell and diameter bin the wall-fed and channel-fed bowls that receive water,
    each with its yearly inflow (m3 per hollow), number, cell and kind (0 depression, 1 wall-fed, 2 channel-fed)."""
    pools = B.build_pools(inp, b, name)
    m = pools['state'] == STATES.index(state)
    h = {k: (v[m] if isinstance(v, np.ndarray) and v.shape[:1] == m.shape else v) for k, v in pools.items()}
    h['top'] = h['top'].copy()
    dep = h['kind'] == 0
    if dep.sum() != cap_d.size:
        raise ValueError('phase 2b\'s depressions are not the drylands\' depressions')
    h['top'][dep] = cap_d
    h['q'] = h['inflow'].sum(axis=1)
    return h


def store_curve(h, ncell, dry_m2, fills=FILLS):
    """Per cell: the capacity (m over the grid's dryland area) and, at each share of it, the water area (share of that
    area) and the seepage (m/s over it) under each floor, the cell's water shared among its hollows by their yearly
    inflow, each hollow at most full (the hollow i holding min(f C q_i / Q, c_i) at the share f of the cell's capacity
    C, with Q the cell's inflow), summed over the hollows with their numbers."""
    w = h['weight']
    cell = h['cell']
    clog = B.area_at_stage(np.full(h['top'].size, B.BANK['clogged_stage_m']), h['s'], h['A'])
    cap_c = np.bincount(cell, weights=w * h['top'], minlength=ncell)
    q_c = np.bincount(cell, weights=w * h['q'], minlength=ncell)
    share = np.where(q_c[cell] > 0, h['q'] / np.maximum(q_c[cell], 1e-30), 0.0)
    area = np.zeros((ncell, fills.size))
    seep = {f: np.zeros((ncell, fills.size)) for f in FLOORS}
    held = np.zeros((ncell, fills.size))
    for k, f in enumerate(fills):
        if f <= 0.0:
            continue
        v = np.minimum(f * cap_c[cell] * share, h['top'])
        _, a = geometry(v, h['s'], h['A'], h['V'])
        dmean = np.where(a > 0.0, v / np.maximum(a, 1e-9), 0.0)
        area[:, k] = np.bincount(cell, weights=w * a, minlength=ncell) / dry_m2
        held[:, k] = np.bincount(cell, weights=w * v, minlength=ncell)
        for floor in FLOORS:
            seep[floor][:, k] = np.bincount(cell, weights=w * seepage(dmean, a, clog, floor), minlength=ncell) / dry_m2
    with np.errstate(invalid='ignore', divide='ignore'):
        held_share = np.where(cap_c[:, None] > 0, held / np.maximum(fills[None, :] * cap_c[:, None], 1e-30), 1.0)
    return cap_c / dry_m2, area, seep, held_share


def read_curve(fill, area, seep, cap, v):
    """The store's water area and seepage at volume v, read as climate/gcm/land_store.curve reads the tables: linear in
    the square root of the share of the capacity, nil at nil volume, held at the full store's above its capacity."""
    cap, v = np.asarray(cap, float), np.asarray(v, float)
    on = (cap > 0.0) & (v > 0.0)
    f = np.where(on, np.minimum(v / np.where(cap > 0.0, cap, 1.0), 1.0), 0.0)
    zs, root = np.sqrt(f), np.sqrt(np.asarray(fill, float))
    k = np.clip(np.searchsorted(root, zs, side='left'), 1, root.size - 1)
    w = (zs - root[k - 1]) / np.maximum(root[k] - root[k - 1], 1e-30)
    rows = np.arange(cap.size)
    fw = area[rows, k - 1] + w * (area[rows, k] - area[rows, k - 1])
    q = seep[rows, k - 1] + w * (seep[rows, k] - seep[rows, k - 1])
    return np.where(on, fw, 0.0), np.where(on, q, 0.0)


def pool_rate(e_day, d):
    """Phase 2b's pools' evaporation over the lunar day (m/s) at the store's mean depth d (m), linear in the logarithm
    of the depth between phase 2b's reference depths and held at their ends (as the offline check reads it)."""
    ln_d = np.log(B.D_REF)
    ld = np.clip(np.log(np.maximum(d, B.D_REF[0])), ln_d[0], ln_d[-1])
    j = np.clip(np.searchsorted(ln_d, ld) - 1, 0, B.D_REF.size - 2)
    w = (ld - ln_d[j]) / (ln_d[j + 1] - ln_d[j])
    rows = np.arange(e_day.shape[0])
    return (1.0 - w) * e_day[rows, j] + w * e_day[rows, j + 1]


def balance_residual(v, fill, area, seep, cap, inflow, rain, e_day):
    """The store's mean gain (m/s over the cell) at volume v: the runoff reaching the hollows on the cell's land
    outside the store's water, inflow (1 - a), and the rain on its water, rain a, less its evaporation at the pools'
    lunar-day rate at its mean depth and its seepage."""
    fw, q = read_curve(fill, area, seep, cap, v)
    d = np.where(fw > 0.0, v / np.maximum(fw, 1e-30), B.D_REF[0])
    return inflow * (1.0 - fw) + fw * (rain - pool_rate(e_day, d)) - q


def balanced_volume(fill, area, seep, cap, inflow, rain, e_day, iterations=60):
    """Each cell's balanced volume, as phase 2b's balance_start finds its pools': the volume at which the store's mean
    inflow and rain meet its mean evaporation and seepage, by halving the range from nil to the capacity, or the
    capacity where the store gains even when full; nil where no runoff reaches it."""
    lo, hi = np.zeros(cap.size), cap.astype(float).copy()
    for _ in range(iterations):
        mid = 0.5 * (lo + hi)
        up = balance_residual(mid, fill, area, seep, cap, inflow, rain, e_day) > 0.0
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    full = balance_residual(cap, fill, area, seep, cap, inflow, rain, e_day) >= 0.0
    return np.where(cap <= 0.0, 0.0, np.where(full, cap, np.where(inflow > 0.0, 0.5 * (lo + hi), 0.0)))


def interpolation_error(fill, area, seep, cap, area_mid, seep_mid, dry_m2):
    """The store's curve read between its fills (read_curve) against the hollows' own sum at the geometric middle of
    each pair of fills: the relative error of the area and of the clay-floored seepage over the cells with a store,
    as the largest over the segments (area-weighted mean and maximum) and per segment (area-weighted mean)."""
    mids = np.sqrt(FILLS[1:-1] * FILLS[2:])
    on = cap > 0
    out = {}
    for name, tab, true in (('area', area, area_mid), ('seepage_clay', seep, seep_mid)):
        errs = []
        for m, f in enumerate(mids, 1):
            est, _ = read_curve(fill, tab, tab, cap, f * cap)
            with np.errstate(invalid='ignore', divide='ignore'):
                errs.append(np.where(true[:, m] > 0.0, np.abs(est / true[:, m] - 1.0), 0.0))
        errs = np.array(errs)
        e = errs.max(axis=0)
        out[name] = dict(max=float(e[on].max()), area_weighted_mean=float((e * dry_m2)[on].sum() / dry_m2[on].sum()),
                         by_segment_area_weighted_mean=[float((v * dry_m2)[on].sum() / dry_m2[on].sum()) for v in errs])
    out['segment_middles'] = [float(v) for v in mids]
    return out


def routing_shares(inp, x):
    """Per library, surface and hollow age, each cell's shares of its slopes' runoff: phase 2a's fates over the runoff
    they came from, both mm a year over the cell (phase2b_inputs' fates and off). Where flows crossing a T21 border
    leave a cell's hollows, rivers and recharge above its own runoff, the three are scaled to it and the soil's share
    is nil; the product records how often and by how much."""
    years = x['pr'].shape[0] / 120.0
    yr = JULIAN_YEAR_DAYS / (12.0 * SYNODIC_MONTH_DAYS)
    out, checks = {}, {}
    for b in LIBRARIES:
        for name in SURFACES:
            rt = R1.by_tenth(inp[f'off__{b}__{name}'].astype(float), x['tenth']).sum(axis=0) / years * yr
            for s in STATES:
                f = {q: inp[f'{q}__{b}__{name}__{s}'].astype(float).sum(axis=0)
                     for q in ('bowls_inside', 'bowls_slope', 'bowls_channel', 'depressions', 'slope_recharge',
                               'bed_recharge', 'sea', 'humid_land')}
                with np.errstate(invalid='ignore', divide='ignore'):
                    share = lambda v: np.where(rt > 0.0, v / np.maximum(rt, 1e-30), 0.0)
                    hol = share(f['bowls_inside'] + f['bowls_slope'] + f['bowls_channel'] + f['depressions'])
                    riv = share(f['sea'] + f['humid_land'])
                    rch = share(f['slope_recharge'] + f['bed_recharge'])
                tot = hol + riv + rch
                over = tot > 1.0
                scale = np.where(over, 1.0 / np.maximum(tot, 1e-30), 1.0)
                hol, riv, rch = hol * scale, riv * scale, rch * scale
                soil = np.where(over, 0.0, np.maximum(1.0 - hol - riv - rch, 0.0))
                out[(b, name, s)] = np.stack([hol, riv, rch, soil], axis=1)
                a = x['area_km2']
                checks[f'{b}__{name}__{s}'] = dict(cells_scaled=int(over.sum()), largest_sum=float(tot.max()),
                                                   area_share_scaled=float(a[over].sum() / a.sum()),
                                                   hollows_share_drylands=float((hol * rt * a).sum() / (rt * a).sum()))
    return out, checks


def mean_depth_mm(lib):
    """The library's mean storm depth (mm), as phase 1's Tables weighs it."""
    w, d = np.asarray(lib['weights'], float), ST.depths(lib)
    return float((w * d).sum() / w.sum())


def per_mm(tot, mu, lam, sat):
    """The runoff per millimetre of a 3-day mean's convective rain: phase 2a's runoff of a mean holding N storms
    (tot, N = 0..n_max on phase 2a's saturations) mixed over a Poisson number of storms of mean lam, over the mean's
    rain lam mu (mu the library's mean storm depth); at lam nil the limit, one storm's runoff over mu."""
    lam, sat = np.asarray(lam, float), np.asarray(sat, float)
    mix = P2.poisson_mix_w(tot, lam.ravel(), sat.ravel()).reshape(lam.shape)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(lam > 0.0, mix / np.maximum(lam * mu, 1e-30), P2.interp_w(tot[1], sat) / mu)


def convective_tables(x, libs, pressure):
    """Per library and surface, the slopes' runoff per millimetre of convective rain (RAIN_MM_DAY, WETNESS): phase
    2a's central wetting tables (wetting.mean_storms, the expected runoff of a 3-day mean holding N storms at random
    times, the soil redistributing and drying between them) mixed over the Poisson number of storms a 3-day mean of
    that rain brings, on soil at the saturation the bucket's relative wetness gives; at nil rain the limit, one storm's
    runoff over its depth. Returns the tables and phase 2a's own rows (for the reproduction check)."""
    dem = P2.demand_mm_h(x)
    wt = x['prc'].astype(float)
    e_rain = float((dem * wt).sum() / wt.sum())
    out, rows = {}, {}
    for b in LIBRARIES:
        mu = mean_depth_mm(libs[b])
        for n in SURFACES:
            p = IN.parameters(n, pressure_pa=pressure)
            c, ps = P2.curve_case(IN.SURFACES[n]['soil'], 'central')
            tot, _, _ = W.mean_storms(p, c, libs[b], P2.SAT_W, P2.N_WET[b], P2.MEAN_HOURS, e_rain, p_shape=ps)
            s = IN.initial_saturation(p['soil'], WETNESS * IN.WSMAX_M)
            lam, sat = np.meshgrid(RAIN_MM_DAY * R1.MEAN_DAYS / mu, s, indexing='ij')
            out[(b, n)] = per_mm(tot, mu, lam, sat)
            rows[(b, n)] = dict(tot=tot, mean_depth_mm=mu, p=p)
            print(f'convective table {b} {n}: {out[(b, n)].min():.4f}-{out[(b, n)].max():.4f}', file=sys.stderr)
    return out, rows, e_rain


def large_scale_tables(pressure):
    """Per surface, the runoff per millimetre of large-scale rain (RAIN_MM_DAY, WETNESS): a 3-day mean of that rain
    falling as one steady rain of phase 1's 2 mm/h on soil at the saturation the bucket's wetness gives (run.runoff's
    reading, which phase 2a keeps); nil where the mean brings 0.02 mm or less or the rate does not exceed the subsoil's
    conductivity, as there."""
    out = {}
    for n in SURFACES:
        p = IN.parameters(n, pressure_pa=pressure)
        s = IN.initial_saturation(p['soil'], WETNESS * IN.WSMAX_M)
        depth, sat = np.meshgrid(RAIN_MM_DAY * R1.MEAN_DAYS, s, indexing='ij')
        g = np.zeros(depth.shape)
        sel = depth > 0.02
        if R1.LS_RATE_MM_H > p['k_sub'] and sel.any():
            g[sel] = IN.steady_rain(p, depth[sel], R1.LS_RATE_MM_H, IN.deficit(p, sat[sel]), np.zeros(int(sel.sum()))) / depth[sel]
        out[n] = g
    return out


def read_table(table, r, w):
    """A runoff table read as climate/gcm/land_store.runoff_share reads it: bilinear in the rain's 3-day mean (mm/day)
    and the bucket's relative wetness, both held within the grids."""
    r = np.clip(np.asarray(r, float), RAIN_MM_DAY[0], RAIN_MM_DAY[-1])
    w = np.clip(np.asarray(w, float), WETNESS[0], WETNESS[-1])
    i = np.clip(np.searchsorted(RAIN_MM_DAY, r, side='left') - 1, 0, RAIN_MM_DAY.size - 2)
    j = np.clip(np.searchsorted(WETNESS, w, side='left') - 1, 0, WETNESS.size - 2)
    a = (r - RAIN_MM_DAY[i]) / (RAIN_MM_DAY[i + 1] - RAIN_MM_DAY[i])
    c = (w - WETNESS[j]) / (WETNESS[j + 1] - WETNESS[j])
    return ((1.0 - a) * (1.0 - c) * table[i, j] + a * (1.0 - c) * table[i + 1, j] + (1.0 - a) * c * table[i, j + 1]
            + a * c * table[i + 1, j + 1])


def reproduction(x, inp, conv_tab, rows, ls_tab):
    """The runoff tables against phase 2a's central runoff (phase2b_inputs' off) on the design years' 3-day means, each
    mean read at its own convective and large-scale rain and bucket: the drylands' totals with phase 2a's own mixing
    recomputed here (which must give phase 2a's runoff) and through the tables, and the tables' error per cell-year."""
    conv_mm = x['prc'].astype(float) * R1.MM_PER_MEAN
    ls_mm = np.maximum(x['pr'].astype(float) * R1.MM_PER_MEAN - conv_mm, 0.0)
    rate = lambda mm: mm / R1.MEAN_DAYS                                         # mm/day
    wet = np.clip(x['mrso'].astype(float) / IN.WSMAX_M, 0.0, 1.0)
    a = x['area_km2'].astype(float)
    out = {}
    for b in LIBRARIES:
        for n in SURFACES:
            r = rows[(b, n)]
            p = r['p']
            s = IN.initial_saturation(p['soil'], x['mrso'].astype(float))
            direct = P2.poisson_mix_w(r['tot'], conv_mm / r['mean_depth_mm'], s)
            ls = np.zeros_like(conv_mm)
            sel = ls_mm > 0.02
            if R1.LS_RATE_MM_H > p['k_sub'] and sel.any():
                ls[sel] = IN.steady_rain(p, ls_mm[sel], R1.LS_RATE_MM_H, IN.deficit(p, s[sel]), np.zeros(int(sel.sum())))
            table = read_table(conv_tab[(b, n)], rate(conv_mm), wet) * conv_mm + read_table(ls_tab[n], rate(ls_mm), wet) * ls_mm
            saved = inp[f'off__{b}__{n}'].astype(float)
            tot = lambda v: float((v.sum(axis=0) * a).sum())
            with np.errstate(invalid='ignore', divide='ignore'):
                err = np.where(saved.sum(axis=0) > 0, np.abs(table.sum(axis=0) / saved.sum(axis=0) - 1.0), 0.0)
            out[f'{b}__{n}'] = dict(direct_over_phase2a=tot(direct + ls) / tot(saved), tables_over_phase2a=tot(table) / tot(saved),
                                    large_scale_share_of_runoff=tot(ls) / max(tot(direct + ls), 1e-30),
                                    cell_year_error_area_weighted_mean=float((err * a).sum() / a.sum()),
                                    cell_year_error_max=float(err.max()))
            print(f'reproduction {b} {n}: {out[f"{b}__{n}"]}', file=sys.stderr)
    return out


def channel_beds(product_2a):
    """The channel beds' share of the slopes' runoff by class, library, surface and hollow age (phase 2a's class sums;
    the per-cell beds' take is not kept in phase 2b's inputs, so the store's tables give it within the soil's share)."""
    fates = product_2a['routing']['fates']
    return {b: {n: {s: {c: fates[b][n][s][c]['channel_bed']['share_of_runoff'] for c in R1.CLASSES + ('drylands',)}
                    for s in STATES} for n in SURFACES} for b in LIBRARIES}


def sig(v, digits=DIGITS):
    """A float kept to the product's significant digits."""
    v = float(v)
    return float(f'{v:.{digits}g}') if v != 0.0 else 0.0


def main():
    x, meta = EX.load()
    with np.load(P2.PHASE2B_INPUTS) as z:
        meta_in = json.loads(str(z['metadata']))
        inp = {k: z[k] for k in z.files if k != 'metadata'}
    ncell = x['cls'].size
    dry_km2 = inp['t21_dry_km2'].astype(float)
    dry_m2 = np.maximum(dry_km2, 1e-9) * 1e6
    libs = R1.libraries(R1.STORMS)
    pressure = R1.mean_pressure_pa(x)
    mean_s = R1.MEAN_DAYS * 86400.0

    # The runoff tables and their reproduction of phase 2a's runoff.
    conv_tab, conv_rows, e_rain = convective_tables(x, libs, pressure)
    ls_tab = large_scale_tables(pressure)
    repro = reproduction(x, inp, conv_tab, conv_rows, ls_tab)

    # The pools' lunar-day evaporation by depth (phase 2b's energy balance), the record's rain and the routing.
    f = B.forcing(x)
    e_tab, _ = B.energy_balance(f, IN.MOON)                                   # (cell, depth, tenth) mm/h
    e_day = e_tab.mean(axis=2) / 1000.0 / 3600.0                             # (cell, depth) m/s
    rain_mean = x['pr'].astype(float).mean(axis=0)                           # m/s
    conv_mean = x['prc'].astype(float).mean(axis=0)
    ls_mean = np.maximum(x['pr'].astype(float) - x['prc'].astype(float), 0.0).mean(axis=0)
    shares, share_checks = routing_shares(inp, x)

    # The store's curves by inflow and their balanced volumes, per library, surface and hollow age.
    cap_d, wide = depression_capacity(inp)
    curves, interp, composition, unseen = {}, {}, {}, {}
    for b in LIBRARIES:
        for n in SURFACES:
            runoff = inp[f'off__{b}__{n}'].astype(float).mean(axis=0) / 1000.0 / mean_s     # m/s over the cell
            for s in STATES:
                h = hollows(inp, b, n, s, cap_d)
                cap, area, seep, held = store_curve(h, ncell, dry_m2)
                _, area_mid, seep_mid, _ = store_curve(h, ncell, dry_m2, fills=np.r_[0.0, np.sqrt(FILLS[1:-1] * FILLS[2:])])
                inflow = shares[(b, n, s)][:, 0] * runoff
                vbal = {fl: balanced_volume(FILLS, area, seep[fl], cap, inflow, rain_mean, e_day) for fl in FLOORS}
                curves[(b, n, s)] = dict(cap=cap, area=area, seep=seep, vbal=vbal, inflow=inflow)
                key = f'{b}__{n}__{s}'
                interp[key] = interpolation_error(FILLS, area, seep['clay_floored'], cap, area_mid,
                                                  seep_mid['clay_floored'], dry_m2)
                w = h['weight'] * h['top']
                wide_h = np.zeros(h['kind'].size, bool)
                wide_h[h['kind'] == 0] = wide
                parts = dict(bowls_wall_fed=h['kind'] == 1, bowls_channel_fed=h['kind'] == 2,
                             depressions_within_sixteen_cells=(h['kind'] == 0) & ~wide_h,
                             depressions_wider=(h['kind'] == 0) & wide_h)
                composition[key] = {k: float(w[m].sum() / w.sum()) for k, m in parts.items()}
                composition[key]['capacity_m_quantiles'] = [float(v) for v in np.quantile(cap[cap > 0], [0.1, 0.5, 0.9])]
                on = cap > 0
                unseen[key] = dict(share_held_by_fill=[float((held[on, k] * dry_m2[on]).sum() / dry_m2[on].sum())
                                                       for k in range(1, FILLS.size)],
                                   balanced_fill_quantiles={fl: [float(v) for v in np.quantile(vbal[fl][on] / cap[on], [0.1, 0.5, 0.9])]
                                                            for fl in FLOORS})
                print(f'curve {key}: capacity median {np.median(cap[on]):.2f} m, area full {np.median(area[on, -1]):.4f}, '
                      f'balanced clay median {np.median(vbal["clay_floored"][on]):.3f} m, {interp[key]}', file=sys.stderr)
                del h

    product_2a = json.loads(R1.OUT.read_text())['phase2a']
    cells = []
    for c in range(ncell):
        store = {b: {n: {s: dict(capacity_m=sig(curves[(b, n, s)]['cap'][c]),
                                 area_share=[sig(v) for v in curves[(b, n, s)]['area'][c]],
                                 seepage_m_s={fl: [sig(v) for v in curves[(b, n, s)]['seep'][fl][c]] for fl in FLOORS},
                                 balanced_m={fl: sig(curves[(b, n, s)]['vbal'][fl][c]) for fl in FLOORS},
                                 inflow_m_s=sig(curves[(b, n, s)]['inflow'][c]))
                         for s in STATES} for n in SURFACES} for b in LIBRARIES}
        row = dict(j=int(x['j'][c]), i=int(x['i'][c]), lat=R1.r3(x['lat'][c]), lon_e=R1.r3(x['lon'][c]),
                   unep=str(x['cls'][c]), area_km2=R1.r3(x['area_km2'][c]), dryland_km2=R1.r3(dry_km2[c]),
                   rain_m_s=dict(convective=sig(conv_mean[c]), large_scale=sig(ls_mean[c]), total=sig(rain_mean[c])),
                   pool_evaporation_m_s=[sig(v) for v in e_day[c]],
                   store=store,
                   shares={b: {n: {s: dict(zip(('hollows', 'rivers', 'recharge', 'soil'),
                                               (sig(v) for v in shares[(b, n, s)][c])))
                                   for s in STATES} for n in SURFACES} for b in LIBRARIES})
        cells.append(row)
    product = dict(
        schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE,
        basis=dict(run=meta['run'], years=meta['years'], cells=ncell, drylands_pressure_pa=pressure,
                   moon_gravity_in_earths=IN.MOON, mean_days=R1.MEAN_DAYS, large_scale_rate_mm_h=R1.LS_RATE_MM_H,
                   storms_per_mean_simulated=P2.N_WET, demand_mm_h_rain_weighted=e_rain,
                   storm_mean_depth_mm={b: mean_depth_mm(libs[b]) for b in LIBRARIES}),
        grid=dict(nlon=int(x['grid_lon'].size), nlat=int(x['grid_lat'].size)),
        grids=dict(fill=[float(v) for v in FILLS], rain_mm_day=[float(v) for v in RAIN_MM_DAY],
                   wetness=[float(v) for v in WETNESS], pool_depth_m=[float(v) for v in B.D_REF]),
        floors=list(FLOORS), surfaces=list(SURFACES), libraries=list(LIBRARIES), hollows=list(STATES),
        floor_designs=dict(clay_floored=dict(seal=B.SEAL), sandy_banked=dict(seal=B.SEAL, bank=B.BANK),
                           tight_seal=dict(seal=TIGHT, note='phase 2b\'s sensitivity case, the loam only there')),
        surface_designs={n: IN.SURFACES[n] for n in SURFACES},
        convective_runoff={b: {n: [[sig(v) for v in row_] for row_ in conv_tab[(b, n)]] for n in SURFACES} for b in LIBRARIES},
        large_scale_runoff={n: [[sig(v) for v in row_] for row_ in ls_tab[n]] for n in SURFACES},
        checks=dict(reproduction_of_phase2a_runoff=repro, capacity_composition=composition,
                    interpolation_between_fills=interp, water_held_by_the_hollows=unseen,
                    routing_shares=share_checks, channel_beds_share_of_runoff=channel_beds(product_2a)),
        cells=cells)
    files = {p: R1.sha256(p) for p in CODE}
    product['producer'] = dict(study='research/studies/storm_water', exporter='store_tables.py', files=files,
                               inputs={str(P2.PHASE2B_INPUTS.relative_to(R1.ROOT)): R1.sha256(P2.PHASE2B_INPUTS.relative_to(R1.ROOT)),
                                       str(EX.OUT.relative_to(R1.ROOT)): meta['npz']['sha256'],
                                       str(R1.OUT.relative_to(R1.ROOT)): R1.sha256(R1.OUT.relative_to(R1.ROOT)),
                                       R1.STORMS: R1.sha256(R1.STORMS)},
                               inputs_schema=meta_in['schema'], drainage_sources=meta_in['drainage_sources'],
                               constants=constants_used(CODE))
    head = json.dumps({k: v for k, v in product.items() if k != 'cells'}, indent=1, default=float)
    rows = ',\n'.join('  ' + json.dumps(c, separators=(',', ':')) for c in cells)
    OUT.write_text(head[:-2] + ',\n "cells": [\n' + rows + '\n ]\n}\n')
    return product


EVIDENCE = ('A model result of the offline storm-water chain on the design years (A28R_301k_design_l15_snow_s934, years '
            '10-19) and geography\'s 1.9 km depression tree, for PlaSim\'s hollows\' store (climate/gcm/land_store.py). '
            'The store\'s curves sum phase 2b\'s hollows with a cell\'s water shared among them by their own yearly '
            'inflow, each at most full, where phase 2b\'s pools fill, spill and dry each on its own; the cost of that '
            'lumping is measured offline (session-tools/agents/storm_water/stage_b2_report.md). The balanced volumes '
            'meet the store\'s mean inflow with its mean losses at phase 2b\'s lunar-day pool evaporation, as phase 2b '
            'starts its pools. The floors are designs (a clay seal scaled from the HAPEX-Sahel plateau pools, sandy banks '
            'above a 1 m stage, a seal ten times tighter), the bowls\' counts the terrain study\'s equilibrium law, their '
            'shape half Pike\'s fresh depth, and the routing shares phase 2a\'s routing under stand-in storms (Walnut '
            'Gulch\'s measured family and CM1\'s 6 km window over wet ground). The runoff tables are phase 2a\'s own '
            'wetting tables for those storms and phase 1\'s steady 2 mm/h large-scale rain (an estimate), Green and Ampt '
            'at lunar gravity on designed surfaces; read at each 3-day mean\'s own rain and bucket they give phase 2a\'s '
            'runoff. No PlaSim run has used these tables.')
READING_RULE = ('Each cell is a T21 point j * nlon + i of PlaSim\'s grid (north to south from 0 E). Depths, areas, rates '
                'and seepage are per unit of the grid\'s dryland area in the cell (dryland_km2), the basis of phase 2a\'s '
                'fates and phase 2b\'s per-cell results, which PlaSim carries over the whole T21 cell as it carries the '
                'bucket. store[library][surface][age].capacity_m is the hollows\' total volume (m of water) and '
                'area_share and seepage_m_s (per floor) the water area\'s share and the floors\' seepage (m/s) with the '
                'store holding the share grids.fill of its capacity, the water shared among the hollows by their yearly '
                'inflow and each hollow at most full; read between fills linearly in the square root of the fill and nil '
                'at nil volume. balanced_m (per floor) is the store\'s balanced volume, its start where a run holds no '
                'store, and inflow_m_s the mean runoff reaching the hollows that the balance takes, before the store\'s '
                'own water area replaces its share of the slopes (the balance: inflow_m_s (1 - a) + a (rain_m_s.total - '
                'e) - seepage = 0, with a the area share and e pool_evaporation_m_s at the store\'s mean depth on '
                'grids.pool_depth_m, linear in the depth\'s logarithm). shares[library][surface][age] split the slopes\' '
                'runoff into the hollows, the rivers (the sea and the humid land), the ground below the roots (recharge) '
                'and the soil (the slopes\' and channel beds\' soak above the roots); they sum to 1. '
                'convective_runoff[library][surface][rain][wetness] is the slopes\' runoff per millimetre of convective '
                'rain and large_scale_runoff[surface][rain][wetness] per millimetre of large-scale rain, bilinear in '
                'grids.rain_mm_day (a 3-day mean of that rain, mm/day; the tables are built on means of basis.mean_days) '
                'and grids.wetness (the bucket\'s water over its 0.5 m capacity), held within both grids. rain_m_s holds '
                'the record\'s mean rain, the start of PlaSim\'s running means.')


if __name__ == '__main__':
    main()
