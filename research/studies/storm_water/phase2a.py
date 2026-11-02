"""The storm-water study, phase 2a: the wetting between storms, and the water the drylands' storms deliver to the hollows.

run.py calls phase2a() after phase 1 and adds its fields to results/storm_water.json. It

1. replaces phase 1's two treatments of the wetting between a mean's storms with wetting.py's redistribution and drying
   between storms at uniformly random times within the mean, a central estimate and a range over the soil curves'
   spread (the bubbling pressure over one standard deviation of the class's suction, Smith et al.'s profile factor 1-2);
2. routes each surface's central runoff on geography's 1.9 km grid (geography.drainage.compute with the runoff
   supplied) into the depression tree, event by event (routing.py): the slopes' run-on infiltration, the crater bowls
   below the grid for young and for mature land, Lane's channel losses recast for gravity, and the storms' footprint;
3. carries the slopes' and channel beds' losses through a lunar day of redistribution in the soil and the bed to find
   what drains below the root zone, the focused recharge;
4. sums, per T21 dryland cell and per surface, by tenth and by year, the water reaching closed depressions and the
   bowls, the water lost in transit and the part of it below the root zone, and per depression its inflow by tenth
   (research/runs/storm_water/hollows_2a.npz, phase 2b's input), and the leverage against the plan's thresholds;
5. runs the routing chain at Earth's gravity on a Walnut Gulch-like basin of 149 km2 against its measured yield.
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np

from research.studies.storm_water import extract as EX
from research.studies.storm_water import infiltration as IN
from research.studies.storm_water import routing as RT
from research.studies.storm_water import run as R1
from research.studies.storm_water import storms as ST
from research.studies.storm_water import wetting as W
from shared.constants import JULIAN_YEAR_DAYS, SYNODIC_MONTH_DAYS

SAT_W = np.array([0.0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.6])
MEAN_HOURS = R1.MEAN_DAYS * 24.0
N_WET = dict(walnut_gulch=24, cm1_6km=12)        # storms per mean simulated exactly; more extrapolate linearly
SUCTION_SD = dict(sandy_loam=45.47 / 11.01, loam=59.38 / 8.89)    # Rawls et al. (1983) Table 2, one SD of log psi_f
# The range: the bubbling pressure over one SD either way, and the profile factor at 2.0, the side on which Richards'
# equation lies (wetting.richards_redistribution agrees with p = 1.7 to 0.01-0.03 in saturation; p = 1 misses by 0.04-0.05).
RANGE_CASES = dict(bubbling_low=dict(hb=-1), bubbling_high=dict(hb=1), profile_factor_2=dict(p=2.0))
ROOT_ZONE_MM = 1000.0            # design: the drylands' plants root to about a metre
LUNAR_DAY_H = SYNODIC_MONTH_DAYS * 24.0
STATES = ('young', 'mature')
HYPSOMETRY_CELLS = 16           # phase 2b: each depression's lowest cells kept for its area-volume curve
PHASE2B_INPUTS = None           # set below once RUNS is known
RUNS = R1.ROOT / 'research' / 'runs' / 'storm_water'
HOLLOWS = RUNS / 'hollows_2a.npz'
PHASE2B_INPUTS = RUNS / 'phase2b_inputs.npz'
WG_BASIN_KM2, WG_BASIN_MM = 149.0, 2.0          # Stone et al. (2008, abstract, via the plan): about 2 mm a year


def curve_case(soil, case, gravity=IN.MOON):
    c = W.curve(soil, gravity)
    hb = RANGE_CASES.get(case, {}).get('hb', 0)
    if hb:
        c['hb'] *= SUCTION_SD[soil] ** hb
    return c, RANGE_CASES.get(case, {}).get('p', W.P_SHAPE)


def interp_w(row, s):
    return np.interp(np.clip(s, SAT_W[0], SAT_W[-1]), SAT_W, row)


def poisson_mix_w(table, lam, s):
    """Phase 1's Poisson mixture on the wetting tables' saturation grid."""
    n_max = table.shape[0] - 1
    top = float(lam.max()) if lam.size else 0.0
    cap = max(n_max, int(top + 10.0 * math.sqrt(top) + 10))
    p = np.exp(-lam)
    acc = p * interp_w(table[0], s)
    slope = table[n_max] - table[n_max - 1]
    for n in range(1, cap + 1):
        p = p * lam / n
        row = table[n] if n <= n_max else table[n_max] + (n - n_max) * slope
        acc = acc + p * interp_w(row, s)
    return acc


def demand_mm_h(x):
    """Priestley and Taylor's demand of each mean (mm/h), as extract.py forms it."""
    return EX.demand(x['tas'].astype(float), x['ps'].astype(float), (x['rss'] + x['rls']).astype(float)) * 3.6e6


def wetting(x, libs, pressure, cache, ls_cache):
    """The central runoff and its range for every surface and library: per (mean, cell) mm, and the class summaries."""
    dem = demand_mm_h(x)
    w = x['prc'].astype(float)
    e_rain = float((dem * w).sum() / w.sum())
    off, summary = {}, {}
    for name in R1.SURFACES:
        p = IN.parameters(name, pressure_pa=pressure)
        soil = IN.SURFACES[name]['soil']
        off[name], summary[name] = {}, {}
        for b in R1.BRACKETS:
            lib = libs[b]
            tab1 = cache(p, lib, 0)
            lam = x['prc'].astype(float) * R1.MM_PER_MEAN / tab1.mean_depth
            s = IN.initial_saturation(p['soil'], x['mrso'].astype(float))
            _, ls = R1.runoff(x, tab1, 'independent', ls_cache=ls_cache)
            cases = {}
            for case in ('central',) + tuple(RANGE_CASES):
                c, ps = curve_case(soil, case)
                tot, layer, evap = W.mean_storms(p, c, lib, SAT_W, N_WET[b], MEAN_HOURS, e_rain, p_shape=ps)
                conv = poisson_mix_w(tot, lam, s)
                cy, ct, cls = R1.summarise(x, conv + ls)
                cases[case] = dict(classes=cls, cell_year=cy, cell_tenth=ct)
                if case == 'central':
                    off[name][b] = conv + ls
                    left = poisson_mix_w(layer, lam, s)
                    dried = poisson_mix_w(evap, lam, s)
            rows = {}
            for k in R1.CLASSES + ('drylands',):
                vals = [cases[q]['classes'][k]['runoff_ratio'] for q in cases]
                rows[k] = dict(central=cases['central']['classes'][k]['runoff_ratio'], low=min(vals), high=max(vals),
                               by_case={q: cases[q]['classes'][k]['runoff_ratio'] for q in cases},
                               share_of_evaporation=cases['central']['classes'][k]['share_of_evaporation'],
                               largest_tenth_share_of_evaporation=cases['central']['classes'][k]['largest_tenth_share_of_evaporation'],
                               runoff_ratio_by_tenth=cases['central']['classes'][k]['runoff_ratio_by_tenth'],
                               runoff_km3_yr=cases['central']['classes'][k]['runoff_km3_yr'])
            rain = x['pr'].astype(float) * R1.MM_PER_MEAN
            summary[name][b] = dict(classes=rows, cell_year=cases['central']['cell_year'],
                                    cell_tenth=cases['central']['cell_tenth'],
                                    layer_left_share_of_rain=float((left.sum(axis=0) @ x['area_km2']) / (rain.sum(axis=0) @ x['area_km2'])),
                                    dried_between_storms_share_of_rain=float((dried.sum(axis=0) @ x['area_km2']) / (rain.sum(axis=0) @ x['area_km2'])))
            print(f"wetting {name} {b}: semi-arid {rows['semi-arid']['central']:.4f} "
                  f"({rows['semi-arid']['low']:.4f}-{rows['semi-arid']['high']:.4f})", file=sys.stderr)
    return off, summary, dict(demand_mm_h_rain_weighted=e_rain, demand_mm_h_mean=float(dem.mean()))


# ---------------------------------------------------------------------------------------------------------------------
# The grid.

def t21_map(x, rows, cols):
    """T21 dryland cell index (0..840, -1 elsewhere) of every 16 px/deg cell."""
    glat, glon = np.asarray(x['grid_lat'], float), np.asarray(x['grid_lon'], float)
    order = np.argsort(-glat)
    lat_sorted = glat[order]
    edges = np.concatenate([[90.0], 0.5 * (lat_sorted[1:] + lat_sorted[:-1]), [-90.0]])
    lat = 90.0 - (np.arange(rows) + 0.5) / 16.0
    band = np.clip(np.searchsorted(-edges, -lat, side='right') - 1, 0, glat.size - 1)
    j_of_row = order[band]
    dlon = 360.0 / glon.size
    lon = (np.arange(cols) + 0.5) / 16.0
    i_of_col = np.round((lon - glon[0]) / dlon).astype(int) % glon.size
    cmap = np.full((glat.size, glon.size), -1, np.int16)
    cmap[x['j'], x['i']] = np.arange(x['j'].size, dtype=np.int16)
    return cmap[j_of_row][:, i_of_col]


class Grid:
    """geography's depression tree on the 16 px/deg grid, with the drylands' cells and the forest of their slopes."""

    def __init__(self, x, cell_mm_day, rows=180 * 16, cols=360 * 16):
        from geography import drainage as dr
        t21 = t21_map(x, rows, cols)
        field = np.where(t21 >= 0, np.asarray(cell_mm_day, np.float32)[np.maximum(t21, 0)], np.float32(0.0))
        atlas, atlas_path, clim_meta, level, h, lat, cell_km, runoff, r = dr.compute(runoff=field)
        del field
        if h.shape != (rows, cols):
            raise RuntimeError(f'expected the {rows} x {cols} grid, got {h.shape}')
        n = rows * cols
        self.rows, self.cols, self.cell_km = rows, cols, cell_km
        self.sources = {str(atlas_path.relative_to(R1.ROOT)): dr.digest(atlas_path)}
        rec = r['receiver'].copy()
        rec[r['sinks']] = r['sinks']
        self.label = r['label']
        self.level = r['level']
        self.outlet = r['outlet']
        self.hollow = ~self.outlet & (h.ravel() < self.level[self.label])
        self.area = r['area'].astype(np.float32)
        self.upstream_km2 = (r['drained_m2'] / 1e6).astype(np.float32)
        self.sink_cell = r['sinks']
        self.floor_m = h.ravel()[r['sinks']].astype(np.float32)
        # The supplied field routed without losses: what reaches the closed depressions and the sea.
        closed = r['closed']
        gen = float((runoff.ravel().astype(float) * dr.MM_DAY * r['area']).sum())
        self.lossless = dict(generated_m3s=gen, to_closed_depressions_m3s=float(r['inflow'][1:][closed[1:]].sum()),
                             to_closed_depressions_share=float(r['inflow'][1:][closed[1:]].sum()) / max(gen, 1e-30),
                             note='the rest joins lakes and depressions that spill with the climatology\'s own water')
        self.tree(r, h.ravel())
        del r, h, runoff
        self.t21 = t21.ravel()
        dry = (self.t21 >= 0) & ~self.outlet
        self.slope_cells = np.flatnonzero(dry & ~self.hollow).astype(np.int32)
        self.hollow_cells = np.flatnonzero(dry & self.hollow).astype(np.int32)
        loc = np.full(n, -1, np.int32)
        loc[self.slope_cells] = np.arange(self.slope_cells.size, dtype=np.int32)
        r_s = rec[self.slope_cells]
        self.receiver = r_s
        self.nxt = np.where(r_s == self.slope_cells, -1, loc[r_s]).astype(np.int32)
        dest = np.full(self.slope_cells.size, 0, np.int8)              # 0 forest, 1 hollow, 2 sea or lake, 3 humid land
        edge = self.nxt < 0
        dest[edge & self.hollow[r_s]] = 1
        dest[edge & (r_s == self.slope_cells)] = 1
        dest[edge & self.outlet[r_s]] = 2
        dest[edge & (dest == 0)] = 3
        self.dest = dest
        del loc
        row = self.slope_cells // cols
        dr_ = r_s // cols - row
        dc = (r_s % cols - self.slope_cells % cols + cols // 2) % cols - cols // 2
        coslat = np.cos(np.radians(90.0 - (row + 0.5) / 16.0))
        self.length_km = (np.hypot(dr_ * cell_km, dc * cell_km * coslat) * RT.SINUOSITY).astype(np.float32)
        self.length_km = np.where(self.length_km > 0, self.length_km, cell_km * RT.SINUOSITY).astype(np.float32)
        self.order = RT.fronts(self.nxt)
        self.floor_t21 = np.where(self.sink_cell >= 0, self.t21[self.sink_cell], -1).astype(np.int32)
        print(f'grid: {self.slope_cells.size:,} dryland slope cells, {self.hollow_cells.size:,} in hollows, '
              f'{len(self.order)} fronts', file=sys.stderr)

    def tree(self, r, hf, k=HYPSOMETRY_CELLS):
        """Phase 2b's depression tree: each depression's own pass (the height of its spill pass, the higher of the
        two cells across it), the depression it spills into (0 the sea), and its hypsometry below the pass: the
        lowest k cells (index, height, area) in height order, their count and the volume at the pass with the area
        rising linearly from the floor cell's height through each next cell's (geography's cells are points)."""
        n_nodes = r['level'].size
        sf, st = r['spill_from'], r['spill_to']
        out = r['outlet']
        self.pass_m = np.full(n_nodes, np.nan, np.float32)
        self.parent = np.zeros(n_nodes, np.int32)
        ok = sf >= 0
        hp = np.where(out[st[ok]], -np.inf, hf[st[ok]])
        self.pass_m[ok] = np.maximum(hf[sf[ok]], hp)
        self.parent[ok] = r['label'][st[ok]]
        cells = np.flatnonzero(~out & (r['label'] > 0))
        lab = r['label'][cells]
        cells = cells[np.lexsort((hf[cells], lab))]
        lab, z = r['label'][cells], hf[cells]
        below = z < self.pass_m[lab]
        cells, lab, z = cells[below], lab[below], z[below]
        a = r['area'][cells]
        start = np.r_[0, np.flatnonzero(lab[1:] != lab[:-1]) + 1]
        size = np.diff(np.r_[start, cells.size])
        rank = np.arange(cells.size) - np.repeat(start, size)
        cum = np.cumsum(a)
        before = cum - a - np.repeat(cum[start] - a[start], size)
        z_next = np.r_[z[1:], 0.0]
        last = np.r_[start[1:] - 1, cells.size - 1]
        z_next[last] = self.pass_m[lab[last]]
        seg = (before + 0.5 * a) * (z_next - z)
        self.n_below = np.bincount(lab, minlength=n_nodes).astype(np.int32)
        self.vmax_m3 = np.bincount(lab, weights=seg, minlength=n_nodes)
        self.area_below_pass = np.bincount(lab, weights=a, minlength=n_nodes)
        keep = rank < k
        self.hyps_cell = np.full((n_nodes, k), -1, np.int32)
        self.hyps_z = np.full((n_nodes, k), np.nan, np.float32)
        self.hyps_a = np.zeros((n_nodes, k), np.float32)
        self.hyps_cell[lab[keep], rank[keep]] = cells[keep]
        self.hyps_z[lab[keep], rank[keep]] = z[keep]
        self.hyps_a[lab[keep], rank[keep]] = a[keep]


def recharge_tables(e_mean):
    """Share of a loss depth below the root zone after a lunar day: channel beds (sand) and slopes (each soil)."""
    depths = np.geomspace(0.1, 3000.0, 60)
    out = dict(depth_mm=depths)
    bed = W.curve('sand')
    se_bed = RT.BED['se_dry']
    below, _, _ = W.below_root_zone(bed, depths, se_bed, LUNAR_DAY_H, e_mean, ROOT_ZONE_MM)
    out['bed'] = below / depths
    for soil in ('sandy_loam', 'loam'):
        c = W.curve(soil)
        se = float(IN.initial_saturation(IN.SOILS[soil], 0.0))
        b, _, _ = W.below_root_zone(c, depths, se, LUNAR_DAY_H, e_mean, ROOT_ZONE_MM)
        out[soil] = b / depths
    return out


def share_below(tab, key, depth_mm):
    return np.interp(np.log(np.maximum(depth_mm, 1e-6)), np.log(tab['depth_mm']), tab[key], left=0.0)


def route_surface(g, x, p, lib, se_ref, rech, gravity=IN.MOON):
    """Per event class and hollow state: the fate of a storm's runoff falling on every dryland cell at once, summed
    per T21 cell (m3 per event), and the inflow per depression (m3 per event)."""
    classes = RT.event_classes(p, lib, se_ref)
    soil = IN.SURFACES[p['name']]['soil'] if p['name'] in IN.SURFACES else 'sandy_loam'
    nk = len(classes) * len(STATES)
    bowls = {s: RT.bowl_classes(s) for s in STATES}
    cols, d0 = [], []
    for c in classes:
        hl, t_eq = RT.hillslope_loss(c['runoff_mm'], c['duration_h'], c['capacity_mm_h'], gravity)
        c['hillslope_loss'] = float(hl)
        c['hillslope_drain_h'] = float(t_eq)
        c['flow_duration_h'] = c['duration_h'] + 0.5 * float(t_eq)
        for s in STATES:
            cols.append((c, s))
            d0.append(c['flow_duration_h'])
    d0 = np.array(d0)
    ncell = x['cls'].size

    def per_t21(cells, values):
        return np.stack([np.bincount(g.t21[cells], weights=values[:, j], minlength=ncell) for j in range(values.shape[1])], axis=1)

    fate = {q: np.zeros((ncell, nk)) for q in ('generated', 'bowls_inside', 'slope_soil', 'bowls_slope', 'channel_bed',
                                               'bowls_channel', 'depressions', 'sea', 'humid_land', 'slope_recharge')}
    local = np.zeros((g.slope_cells.size, nk))
    local_h = np.zeros((g.hollow_cells.size, nk))
    for j, (c, s) in enumerate(cols):
        b = bowls[s]
        inside = RT.bowl_cover(b)
        surv = float(RT.bowl_survival(b, np.array(RT.HILLSLOPE['length_m'] / 2000.0), np.array(1e-6)))
        for cells, arr in ((g.slope_cells, local), (g.hollow_cells, local_h)):
            gen = c['runoff_mm'] / 1000.0 * g.area[cells].astype(float)
            on_slope = gen * (1.0 - inside)
            soak = on_slope * c['hillslope_loss']
            moving = on_slope - soak
            arr[:, j] = moving * surv
            for q, v in (('generated', gen), ('bowls_inside', gen * inside), ('slope_soil', soak),
                         ('bowls_slope', moving * (1.0 - surv))):
                fate[q][:, j] += np.bincount(g.t21[cells], weights=v, minlength=ncell)
            depth = c['runoff_mm'] * c['hillslope_loss']
            fate['slope_recharge'][:, j] += np.bincount(g.t21[cells], weights=soak * share_below(rech, soil, depth),
                                                        minlength=ncell)
    # Water in hollow cells is delivered where it falls.
    fate['depressions'] += per_t21(g.hollow_cells, local_h)
    dep_in = np.stack([np.bincount(g.label[g.hollow_cells], weights=local_h[:, j], minlength=g.level.size)
                       for j in range(nk)], axis=1)
    bowl_list = [bowls[s] for _, s in cols]
    res = route_forest(g, local, d0, gravity, bowl_list, rech)
    for q in ('channel_bed', 'bowls_channel'):
        fate[q] += res[q]
    fate['bowls_channel_by_bin'] = res.pop('bowls_channel_by_bin')
    fate['loss_reach'], fate['bed_m2'] = res.pop('loss_reach'), res.pop('bed_m2')
    out = res['out']
    for code, q in ((1, 'depressions'), (2, 'sea'), (3, 'humid_land')):
        m = g.dest == code
        fate[q] += per_t21(g.slope_cells[m], out[m])
    m = g.dest == 1
    dep_in += np.stack([np.bincount(g.label[g.receiver[m]], weights=out[m, j], minlength=g.level.size) for j in range(nk)], axis=1)
    return classes, cols, fate, dep_in


def route_forest(g, local, d0, gravity, bowl_list, rech):
    """routing.route_events on the drylands' slopes, with each reach's bed losses carried to the recharge."""
    n, k = local.shape
    total = local                                   # the caller's array becomes the running total
    area = np.repeat(g.area[g.slope_cells].astype(float)[:, None], k, axis=1)
    out = np.zeros((n, k), np.float32)
    ncell = int(g.t21.max()) + 1
    acc = {q: np.zeros((ncell, k)) for q in ('channel_bed', 'bowls_channel')}
    nb = bowl_list[0]['number_km2'].size
    by_bin = np.zeros((k, ncell * nb))
    loss_reach = np.zeros((n, k), np.float32)
    width = RT.channel_width_m(g.upstream_km2[g.slope_cells])
    fp = RT.FOOTPRINT_KM2 * 1e6
    t21 = g.t21[g.slope_cells]
    logd, kdt = RT.kd_table(gravity)
    path = np.zeros(n, np.float32)
    chunks = (c for f in g.order for c in np.array_split(f, max(1, f.size // 300000)))
    for f in chunks:
        t, a = total[f], area[f]
        mult = np.maximum(1.0, a / fp)
        v = t / mult
        ln, wd = g.length_km[f].astype(float), width[f].astype(float)
        dur = RT.flow_duration_h(d0[None, :], (path[f] + 0.5 * ln)[:, None], gravity)
        q = RT.lane_outflow(v, ln[:, None], wd[:, None], np.interp(np.log(dur), logd, kdt))
        s = np.ones_like(q)
        for j, b in enumerate(bowl_list):
            s[:, j] = RT.bowl_survival(b, ln, q[:, j])
        loss = mult * (v - q)
        loss_reach[f] = loss
        o = mult * q * s
        for j, b in enumerate(bowl_list):
            # Each crossing takes up to its bowl's volume: the capture splits over the bins as their crossings times
            # the share of the flow each takes (routing.bowl_survival's terms).
            wi = (b['number_km2'] * b['diameter_km'])[None, :] * ln[:, None] * np.minimum(
                1.0, b['volume_m3'][None, :] / np.maximum(q[:, j], 1e-9)[:, None])
            cap = (mult[:, j] * q[:, j] * (1.0 - s[:, j]))[:, None] * wi / np.maximum(wi.sum(axis=1, keepdims=True), 1e-30)
            by_bin[j] += np.bincount((t21[f].astype(np.int64)[:, None] * nb + np.arange(nb)[None, :]).ravel(),
                                     weights=cap.ravel(), minlength=ncell * nb)
        for name, val in (('channel_bed', loss), ('bowls_channel', mult * q * (1.0 - s))):
            for j in range(k):
                acc[name][:, j] += np.bincount(t21[f], weights=val[:, j], minlength=ncell)
        out[f] = o
        r = g.nxt[f]
        go = r >= 0
        if go.any():
            with np.errstate(divide='ignore', invalid='ignore'):
                live = np.where(t > 0.0, o / t, 0.0)
            np.add.at(total, r[go], o[go])
            np.add.at(area, r[go], (a * live)[go])
            np.maximum.at(path, r[go], (path[f] + ln)[go].astype(np.float32))
    acc['bowls_channel_by_bin'] = by_bin.reshape(k, ncell, nb).transpose(1, 2, 0)      # (cell, bin, column)
    acc['out'] = out
    acc['loss_reach'] = loss_reach
    acc['bed_m2'] = (width * g.length_km * 1000.0).astype(np.float32)
    return acc


def delivered(x, off, fate, classes, cols, dep_in, g, rech):
    """Fate per T21 cell by tenth and year (mm a year over the cell) and the depressions' inflow by tenth (m3 a year),
    the classes' event rates set so that each class carries its share of the cell's central runoff."""
    years = x['pr'].shape[0] / 120.0
    yr = JULIAN_YEAR_DAYS / (12.0 * SYNODIC_MONTH_DAYS)
    rt = R1.by_tenth(off, x['tenth']) / years * yr                         # (tenth, cell) mm a year
    prob = np.array([c['probability'] * c['runoff_mm'] for c in classes])
    share = prob / prob.sum()
    out = {}
    # The beds' recharge: each reach's yearly loss spread over its bed, a lunar day's share at a time, and the part of
    # that depth below the root zone after a lunar day of redistribution and drying.
    rate = np.stack([share[i] * rt.sum(axis=0) / c['runoff_mm'] for i, c in enumerate(classes)])     # (class, cell)
    loss_reach, bed_m2 = fate.pop('loss_reach'), fate.pop('bed_m2')
    by_bin = fate.pop('bowls_channel_by_bin')
    for s in STATES:
        jj = [j for j, (_, st) in enumerate(cols) if st == s]
        with np.errstate(divide='ignore', invalid='ignore'):
            w = sum(share[i] * np.where(fate['generated'][:, j, None] > 0, by_bin[:, :, j] / fate['generated'][:, j, None], 0.0)
                    for i, j in enumerate(jj))
        out.setdefault('bowls_channel_by_bin', {})[s] = rt[:, :, None] * w[None, :, :]       # (tenth, cell, bin)
    t21 = g.t21[g.slope_cells]
    per_day = JULIAN_YEAR_DAYS / SYNODIC_MONTH_DAYS
    for si, s in enumerate(STATES):
        jj = [j for j, (_, st) in enumerate(cols) if st == s]
        yearly = sum(rate[i][t21] * loss_reach[:, j] for i, j in enumerate(jj))                    # m3 a year
        depth = yearly / per_day / np.maximum(bed_m2, 1.0) * 1000.0
        rch = np.bincount(t21, weights=yearly * share_below(rech, 'bed', depth), minlength=rt.shape[1])
        gen = sum(rate[i] * fate['generated'][:, j] for i, j in enumerate(jj))
        with np.errstate(divide='ignore', invalid='ignore'):
            frac = np.where(gen > 0, rch / gen, 0.0)
        out.setdefault('bed_recharge', {})[s] = rt * frac[None, :]
    for q, arr in fate.items():
        if q == 'generated':
            continue
        with np.errstate(divide='ignore', invalid='ignore'):
            frac = np.where(fate['generated'] > 0, arr / fate['generated'], 0.0)        # (cell, column)
        out[q] = {}
        for si, s in enumerate(STATES):
            jj = [j for j, (_, st) in enumerate(cols) if st == s]
            w = np.stack([share[i] * frac[:, j] for i, j in enumerate(jj)], axis=1).sum(axis=1)
            out[q][s] = rt * w[None, :]                                                # (tenth, cell) mm a year
    # Depressions: events a year per point in each class, by tenth, at the cell the water comes from.
    dep = {}
    owner = np.full(g.level.size, -1, np.int64)
    owner[g.label[g.hollow_cells]] = g.t21[g.hollow_cells]
    m = g.dest == 1
    owner[g.label[g.receiver[m]]] = g.t21[g.slope_cells[m]]
    have = np.flatnonzero((dep_in.sum(axis=1) > 0) & (owner >= 0))
    for s in STATES:
        jj = [j for j, (_, st) in enumerate(cols) if st == s]
        tot = np.zeros((have.size, 10))
        for i, j in enumerate(jj):
            rate = share[i] * rt[:, owner[have]] / classes[i]['runoff_mm']           # events a year by tenth
            tot += (rate * dep_in[have, j][None, :]).T
        dep[s] = tot
    return out, have, dep


def leverage(x, fates, off):
    """Share of the drylands' rain that leaves the soil's evaporating store, against rain and evaporation."""
    years = x['pr'].shape[0] / 120.0
    yr = JULIAN_YEAR_DAYS / (12.0 * SYNODIC_MONTH_DAYS)
    rain = R1.by_tenth(x['pr'].astype(float) * R1.MM_PER_MEAN, x['tenth']) / years * yr
    evap = R1.by_tenth(-x['evap'].astype(float) * R1.MM_PER_MEAN, x['tenth']) / years * yr
    rows = {}
    for s in STATES:
        leaves = (fates['depressions'][s] + fates['bowls_inside'][s] + fates['bowls_slope'][s] + fates['bowls_channel'][s]
                  + fates['slope_recharge'][s] + fates['bed_recharge'][s] + fates['sea'][s] + fates['humid_land'][s])
        rows[s] = {}
        for name in R1.CLASSES + ('drylands',):
            m = np.ones(x['cls'].size, bool) if name == 'drylands' else (x['cls'] == name)
            a = x['area_km2'][m]
            lt, rt_, et = (leaves[:, m] * a).sum(axis=1), (rain[:, m] * a).sum(axis=1), (evap[:, m] * a).sum(axis=1)
            rows[s][name] = dict(share_of_rain=float(lt.sum() / rt_.sum()), share_of_evaporation=float(lt.sum() / et.sum()),
                                 share_of_evaporation_by_tenth=[float(v) for v in lt / np.maximum(et, 1e-30)],
                                 largest_tenth_share_of_evaporation=float((lt / np.maximum(et, 1e-30)).max()),
                                 passes_yearly=bool(lt.sum() / et.sum() > R1.THRESHOLDS['yearly']),
                                 passes_tenth=bool((lt / np.maximum(et, 1e-30)).max() > R1.THRESHOLDS['tenth']))
    return rows


def class_sums(x, fates, off):
    """Per dryland class and state: each fate as a share of the rain and of the runoff, by year."""
    years = x['pr'].shape[0] / 120.0
    yr = JULIAN_YEAR_DAYS / (12.0 * SYNODIC_MONTH_DAYS)
    rain = R1.by_tenth(x['pr'].astype(float) * R1.MM_PER_MEAN, x['tenth']) / years * yr
    runoff = R1.by_tenth(off, x['tenth']) / years * yr
    out = {}
    for s in STATES:
        out[s] = {}
        for name in R1.CLASSES + ('drylands',):
            m = np.ones(x['cls'].size, bool) if name == 'drylands' else (x['cls'] == name)
            a = x['area_km2'][m]
            r_, o_ = (rain[:, m] * a).sum(), (runoff[:, m] * a).sum()
            row = dict(runoff_share_of_rain=float(o_ / r_))
            for q, v in fates.items():
                tot = (v[s][:, m] * a).sum()
                row[q] = dict(share_of_runoff=float(tot / max(o_, 1e-30)), share_of_rain=float(tot / r_),
                              km3_yr=float(tot * 1e-6))
            out[s][name] = row
    return out


# ---------------------------------------------------------------------------------------------------------------------
# The Earth check: a Walnut Gulch-like basin of 149 km2 on 1.9 km cells.

def synthetic_basin(cell_km, area_km2=WG_BASIN_KM2):
    """A basin of square cells draining to one mouth, a V valley tilted toward it, routed by geography's D8."""
    from geography import drainage as dr
    side = int(round(math.sqrt(area_km2) / cell_km))
    rr, cc = np.meshgrid(np.arange(side), np.arange(side), indexing='ij')
    h = (2.0 * rr + 0.7 * np.abs(cc - side // 2) + 0.01 * cc).astype(np.float32)
    mouth = np.zeros(h.shape, bool)
    mouth[0, side // 2] = True
    rec = dr.d8_receivers(h, mouth, np.zeros(side), cell_km)
    n = side * side
    nxt = np.where(rec == np.arange(n), -1, rec).astype(np.int64)
    area = np.full(n, (cell_km * 1000.0) ** 2)
    up = np.full(n, cell_km ** 2)
    for f in RT.fronts(nxt):
        r = nxt[f]
        go = r >= 0
        np.add.at(up, r[go], up[f][go])
    return nxt, area, up, np.full(n, cell_km * RT.SINUOSITY), int(np.flatnonzero(mouth.ravel())[0])


def walnut_gulch_basin(libs, cell_km):
    """The routing chain at Earth's gravity on Walnut Gulch-like inputs over a 149 km2 basin: its storms on Rawls'
    sandy loam and loam at the wilting point, no bowls, against the basin's measured yield of about 2 mm a year."""
    nxt, area, up, length, outlet = synthetic_basin(cell_km)
    order = RT.fronts(nxt)
    lib = libs['walnut_gulch']
    out = {}
    for soil in ('sandy_loam', 'loam'):
        p = IN.parameters(soil, gravity=1.0)
        s0 = float(IN.initial_saturation(p['soil'], 0.0))
        classes = RT.event_classes(p, lib, s0)
        storms_yr = ST.SEASON_MM / ST.statistics(lib)['mean_depth_mm']
        point = sum(c['probability'] * c['runoff_mm'] for c in classes) * storms_yr          # mm a year at a point
        yield_mm, slope_mm = 0.0, 0.0
        for c in classes:
            hl, t_eq = RT.hillslope_loss(c['runoff_mm'], c['duration_h'], c['capacity_mm_h'], 1.0)
            local = (c['runoff_mm'] / 1000.0 * area * (1.0 - float(hl)))[:, None]
            lost, caught, o, _ = RT.route_events(order, nxt, local, area, length, up,
                                                 np.array([c['duration_h'] + 0.5 * float(t_eq)]), 1.0, [None])
            events = c['probability'] * storms_yr
            yield_mm += events * float(o[outlet, 0]) / area.sum() * 1000.0
            slope_mm += events * c['runoff_mm'] * (1.0 - float(hl))
        out[soil] = dict(point_runoff_mm_yr=point, point_share_of_rain=point / R1.WG_YEAR_MM,
                         off_the_slopes_share_of_rain=slope_mm / R1.WG_YEAR_MM,
                         basin_yield_mm_yr=yield_mm, basin_share_of_rain=yield_mm / R1.WG_YEAR_MM,
                         basin_km2=float(area.sum() / 1e6))
    gm = math.sqrt(out['sandy_loam']['basin_yield_mm_yr'] * out['loam']['basin_yield_mm_yr'])
    out['rule'] = ('Set before the chain was run: it passes when the geometric mean of the two soils\' yields at the '
                   'basin\'s outlet lies within a factor of 3 of the measured 2 mm a year, the factor allowing for the '
                   'basin\'s soils spanning both textures and for the measured yield\'s year-to-year spread.')
    out['geometric_mean_mm_yr'] = gm
    out['passes'] = bool(WG_BASIN_MM / 3.0 <= gm <= 3.0 * WG_BASIN_MM)
    out['measured_mm_yr'] = WG_BASIN_MM
    return out


def ref_saturation(x, soil):
    """The drylands' convective-rain-weighted effective saturation of a soil at storm time."""
    s = IN.initial_saturation(IN.SOILS[soil], x['mrso'].astype(float))
    w = x['prc'].astype(float) * x['area_km2'][None, :]
    return float((s * w).sum() / w.sum())


def save_phase2b_inputs(g, x, off, p2b, union, inflow):
    """Everything phase 2b reads from the routing, in an ignored file run.py rebuilds: the depression tree's fed
    dryland depressions with their dryland ancestors, each one's pass, parent and lowest cells, the bowls' inflows
    per T21 cell by tenth (mm a year) and their channel capture by diameter bin, and each surface's runoff per mean."""
    dry_node = np.r_[False, g.floor_t21 >= 0]                      # by node (0 is the sea)
    nodes = set(int(v) for v in union)
    front = np.array(sorted(nodes), np.int64)
    while front.size:
        par = g.parent[front]
        par = np.unique(par[(par > 0) & dry_node[par]])
        front = np.array([v for v in par.tolist() if v not in nodes], np.int64)
        nodes.update(front.tolist())
    nodes = np.array(sorted(nodes), np.int64)
    pos = np.searchsorted(union, nodes)
    fed = (pos < union.size) & (union[np.minimum(pos, union.size - 1)] == nodes)
    arrays = dict(node=nodes.astype(np.int32), parent=g.parent[nodes], pass_m=g.pass_m[nodes],
                  floor_m=g.floor_m[nodes - 1], sink_cell=g.sink_cell[nodes - 1].astype(np.int64),
                  floor_t21=g.floor_t21[nodes - 1], n_below=g.n_below[nodes], vmax_m3=g.vmax_m3[nodes],
                  area_below_pass_m2=g.area_below_pass[nodes], hyps_cell=g.hyps_cell[nodes], hyps_z=g.hyps_z[nodes],
                  hyps_a=g.hyps_a[nodes], cell_km=np.float64(g.cell_km), cols=np.int64(g.cols))
    for key, a in inflow.items():
        full = np.zeros((nodes.size, 10), np.float32)
        full[fed] = a[pos[fed]]
        arrays['inflow__' + key] = full
    dry = (g.t21 >= 0) & ~g.outlet
    ncell = x['cls'].size
    arrays['t21_dry_km2'] = np.bincount(g.t21[dry], weights=g.area[dry].astype(float), minlength=ncell) / 1e6
    arrays['t21_reach_km'] = np.bincount(g.t21[g.slope_cells], weights=g.length_km.astype(float), minlength=ncell)
    for (b, name), d in p2b.items():
        arrays[f'off__{b}__{name}'] = np.asarray(off[name][b], np.float32)
        for s in STATES:
            arrays[f'bins__{b}__{name}__{s}'] = d['bins'][s].astype(np.float32)
            for q, v in d['fates'].items():
                arrays[f'{q}__{b}__{name}__{s}'] = v[s].astype(np.float32)
    np.savez_compressed(PHASE2B_INPUTS, metadata=json.dumps(dict(
        schema='terluna.research.storm-water-phase2b-inputs/1', drainage_sources=g.sources,
        rule='depressions as geography.drainage nodes; cells as row * cols + column; fates in mm a year over each '
             'T21 cell by tenth; inflow in m3 a year by tenth; off in mm of runoff per 3-day mean')), **arrays)
    print(f'phase 2b inputs: {nodes.size:,} dryland depressions ({int(fed.sum()):,} fed)', file=sys.stderr)


def run_phase2a(x, libs, pressure, cache, ls_cache):
    """Everything phase 2a adds to the product, and the per-cell rows' additions."""
    off, wet, demand = wetting(x, libs, pressure, cache, ls_cache)
    years = x['pr'].shape[0] / 120.0
    yr = JULIAN_YEAR_DAYS / (12.0 * SYNODIC_MONTH_DAYS)
    loam = off['loam']['walnut_gulch'].sum(axis=0) / years * yr / JULIAN_YEAR_DAYS          # mm/day, per cell
    grid = Grid(x, loam)
    rech = recharge_tables(demand['demand_mm_h_mean'])
    routed, fates_cells, leverage_rows, event_rows, hollow_tabs = {}, {}, {}, {}, {}
    p2b = {}
    for b in R1.BRACKETS:
        routed[b], fates_cells[b], leverage_rows[b], event_rows[b] = {}, {}, {}, {}
        for name in IN.MAIN:
            p = IN.parameters(name, pressure_pa=pressure)
            se = ref_saturation(x, IN.SURFACES[name]['soil'])
            classes, cols, fate, dep_in = route_surface(grid, x, p, libs[b], se, rech)
            fates, have, dep = delivered(x, off[name][b], fate, classes, cols, dep_in, grid, rech)
            p2b[(b, name)] = dict(bins=fates.pop('bowls_channel_by_bin'), fates={
                q: fates[q] for q in ('bowls_inside', 'bowls_slope', 'bowls_channel', 'depressions', 'slope_recharge',
                                      'bed_recharge', 'sea', 'humid_land')})
            routed[b][name] = class_sums(x, fates, off[name][b])
            leverage_rows[b][name] = leverage(x, fates, off[name][b])
            event_rows[b][name] = dict(reference_saturation=se, classes=classes)
            fates_cells[b][name] = fates
            hollow_tabs[(b, name)] = (have, dep)
            d = routed[b][name]['young']['drylands']
            print(f"routed {name} {b}: to depressions {d['depressions']['share_of_runoff']:.3f}, bowls "
                  f"{d['bowls_inside']['share_of_runoff'] + d['bowls_slope']['share_of_runoff'] + d['bowls_channel']['share_of_runoff']:.3f}, "
                  f"slopes {d['slope_soil']['share_of_runoff']:.3f}, channels {d['channel_bed']['share_of_runoff']:.3f} (young)",
                  file=sys.stderr)
    # The depressions' inflow by tenth, phase 2b's input (m3 a year arriving in each tenth).
    RUNS.mkdir(parents=True, exist_ok=True)
    arrays, union = {}, np.unique(np.concatenate([v[0] for v in hollow_tabs.values()]))
    for (b, name), (have, dep) in hollow_tabs.items():
        pos = np.searchsorted(union, have)
        for s in STATES:
            a = np.zeros((union.size, 10), np.float32)
            a[pos] = dep[s]
            arrays[f'{b}__{name}__{s}'] = a
    sink = grid.sink_cell[np.maximum(union - 1, 0)]
    meta = dict(schema='terluna.research.storm-water-hollow-inflow/1', node=('geography.drainage depression node '
                'index (1-based; 0 is the sea)'), sink_cell='row * 5760 + column of each depression\'s floor',
                arrays='inflow in m3 a year arriving in each of the ten tenths of the lunar day, per library, surface '
                       'and hollow state', units='m3/yr', drainage_sources=grid.sources)
    below = np.bincount(grid.label[grid.hollow], weights=grid.area[grid.hollow].astype(float), minlength=grid.level.size)
    meta.update(level_m='spill level above the geoid', floor_m='floor height above the geoid',
                area_below_spill_km2='area of the 1.9 km cells below the spill level')
    np.savez_compressed(HOLLOWS, metadata=json.dumps(meta), node=union.astype(np.int32), sink_cell=sink.astype(np.int64),
                        level_m=grid.level[union].astype(np.float32), floor_m=grid.floor_m[np.maximum(union - 1, 0)],
                        area_below_spill_km2=(below[union] / 1e6).astype(np.float32), **arrays)
    hollow_summary = {}
    for (b, name), (have, dep) in hollow_tabs.items():
        for s in STATES:
            tot = dep[s].sum(axis=1)
            q = np.sort(tot)[::-1]
            hollow_summary.setdefault(b, {}).setdefault(name, {})[s] = dict(
                depressions_fed=int((tot > 0).sum()), inflow_km3_yr=float(tot.sum() * 1e-9),
                share_of_inflow_in_largest_hundred=float(q[:100].sum() / max(q.sum(), 1e-30)),
                median_m3_yr=float(np.median(tot[tot > 0])) if (tot > 0).any() else 0.0)
    save_phase2b_inputs(grid, x, off, p2b, union, arrays)
    check = walnut_gulch_basin(libs, grid.cell_km)
    print('walnut gulch basin', json.dumps(check, default=float)[:400], file=sys.stderr)
    bowls = {s: dict(cover=RT.bowl_cover(RT.bowl_classes(s)),
                     cover_range=[RT.bowl_cover(RT.bowl_classes(s, share=v)) for v in RT.EQUILIBRIUM_SHARE],
                     crossings_per_km=float((RT.bowl_classes(s)['number_km2'] * RT.bowl_classes(s)['diameter_km']).sum()),
                     diameters_m=list(RT.BOWL_RANGES_M[s])) for s in STATES}
    product = dict(
        wetting=dict(model='Smith, Corradini and Melone (1993) redistribution with drying between storms at uniformly '
                     'random times within each mean; Brooks and Corey curves (Rawls et al. 1982 via GSSHA/SWMM5 tables; '
                     'Ks = 2 K_GreenAmpt, Rawls et al. 1983 eq. 6)', demand=demand,
                     curves={k: W.CURVES[k] for k in W.CURVES}, profile_factor=W.P_SHAPE, range_cases=RANGE_CASES,
                     suction_sd_factor=SUCTION_SD,
                     surfaces={n: {b: {k: v for k, v in wet[n][b].items() if k not in ('cell_year', 'cell_tenth')}
                                   for b in R1.BRACKETS} for n in R1.SURFACES}),
        routing=dict(drainage_sources=grid.sources, design=dict(hillslope=RT.HILLSLOPE, footprint_km2=RT.FOOTPRINT_KM2,
                                 width_m_at_1_km2=RT.WIDTH_M_AT_KM2, sinuosity=RT.SINUOSITY, bed=RT.BED,
                                 wave_speed_m_s_earth=RT.WAVE_SPEED_M_S_EARTH,
                                 bowl_depth_ratio=RT.BOWL_DEPTH_RATIO, root_zone_mm=ROOT_ZONE_MM),
                     bowls=bowls, lossless_loam_walnut_gulch=grid.lossless, event_classes=event_rows,
                     fates=routed, hollows=hollow_summary,
                     recharge_share_below_root_zone=dict(depth_mm=[float(v) for v in rech['depth_mm'][::6]],
                                                         **{k: [float(v) for v in rech[k][::6]] for k in ('bed', 'sandy_loam', 'loam')})),
        leverage=leverage_rows, walnut_gulch_basin=check,
        hollow_inflow=dict(path=str(HOLLOWS.relative_to(R1.ROOT)), schema=meta['schema']))
    cells = {}
    for c in range(x['cls'].size):
        row = dict(runoff_ratio_2a={n: dict(walnut_gulch=dict(yearly=R1.r3(wet[n]['walnut_gulch']['cell_year'][c]),
                                                              by_tenth=R1.r3(wet[n]['walnut_gulch']['cell_tenth'][:, c])),
                                            cm1_6km=dict(yearly=R1.r3(wet[n]['cm1_6km']['cell_year'][c])))
                                    for n in R1.SURFACES})
        rt = {}
        for n in IN.MAIN:
            f = fates_cells['walnut_gulch'][n]
            rt[n] = {}
            for s in STATES:
                hollows = f['depressions'][s][:, c] + f['bowls_inside'][s][:, c] + f['bowls_slope'][s][:, c] + f['bowls_channel'][s][:, c]
                soil = f['slope_soil'][s][:, c] + f['channel_bed'][s][:, c]
                below = f['slope_recharge'][s][:, c] + f['bed_recharge'][s][:, c]
                rt[n][s] = {k: dict(yearly=R1.r3(v.sum()), by_tenth=R1.r3(v)) for k, v in
                            (('to_hollows_mm_yr', hollows), ('of_which_grid_depressions_mm_yr', f['depressions'][s][:, c]),
                             ('lost_to_soil_mm_yr', soil), ('below_root_zone_mm_yr', below))}
        row['routed_walnut_gulch'] = rt
        cells[c] = row
    return product, cells
