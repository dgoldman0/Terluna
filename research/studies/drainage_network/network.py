"""The drainage network as the Moon's basins fill: which waters fill and spill into which, at what level, over what
area and volume, how fast, and what the seas at the bottom of the chain do.

    ulimit -v 12000000; nice python -m research.studies.drainage_network.network     # results/network.json, about 12 minutes

geography's drainage holds every water below the atlas's common level as one sea at that level and lets the lakes
above it fill by their own balance and spill. The south sea study (research/studies/south_sea_ocean/inflows.py)
let the hollows below the common level, and then the four smaller seas, fill and spill too, under fixed rules.
This study follows each of the six large seas above and below the common level:

- each sea's water at a level is the land joined to its deepest point under that level: a cell joins at the higher
  of its own height and the lowest chain of passes that links its catchment to the sea's (the routing's flood over
  its catchments' passes, run from that sea alone). That gives each sea's area, volume and loss to the air at any
  level, with every basin its water floods on the way up;
- a sea's rim is the lowest pass to another sea's water, and the water its overflow reaches is found by following
  its spill path downhill through a routing;
- the rivers are routed with each sea taking them in at its own level (routing.py: geography's routing with that one
  change), every other hollow and lake filling or drying by its own balance as geography's drainage does;
- a sea that the rivers reaching it at its rim overfill is filled step by step from the common level, its loss
  to the air growing as it spreads, which gives the time it takes; one they cannot fill settles where its loss
  meets them.

The climate is the design years' throughout, with the routing's own budget beside the climate model's: the rain
and the open water's evaporation stay as the climate model made them for the atlas's waters. `--rain` scales the
rain on land and water alike, the range the climate lane carries for the model's rain, with evaporation held. The
heights are LOLA's at 16 pixels per degree above the dry Moon's geoid; the water's own gravity and the crust's sag,
which the atlas round folds in, are not.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ROOT / 'research' / 'runs' / 'drainage_network'
OUT = HERE / 'results' / 'network.json'
SCHEMA = 'terluna.research.drainage-network/1'
ATLAS = 'geography/results/atlas.json'
DRAINAGE = 'geography/turned_corrected/results/drainage.json'
FRAME = 'climate/results/gcm/turned_surface_8d71fd350ca4c6cd.json'
RUN = 'A28R_301k_design_l15_snow_s938_turned_refit_precession'
CLIMATOLOGY = f'climate/gcm/products/climatology_{RUN}_today.npz'
YEARS = [10, 19]
CODE = ('research/studies/drainage_network/network.py', 'research/studies/drainage_network/routing.py',
        'geography/drainage.py', 'geography/topography.py', 'research/studies/south_sea_ocean/inflows.py')
YEAR_S = JULIAN_YEAR_DAYS * JULIAN_DAY
KM3_YR = YEAR_S / 1e9                   # m3/s to km3 a year
CAP_M = 6000.0                          # a sea's tables reach this far above the common level
BELOW_RIM_M = 0.01                      # a full sea stands this far under its pass
FILL_STEPS = 400                        # levels between the common level and a rim for the time to fill
HISTORY_YEARS = 10000                   # how long the seas' joined history is followed
HISTORY_STEP_YEARS = 5.0
RAINS = (1.0, 0.75, 2.0)                # the design years' rain and the ends of the range the climate lane carries
JOINED_BASINS = 4                       # the largest separate basins listed for each sea's water at its rim
RIVERS_A_SEA = 3                        # the largest rivers into each sea whose long profiles are traced
STEEP_M_KM = 20.0                       # a reach falling faster than this (m per km over one cell) counts as steep
PROFILE_KM = 25.0                       # the spacing of the stored long profiles


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------- where each sea's water reaches at a level

def catchment_graph(h, lat, cell_km):
    """Every pit's catchment on the grid and the lowest pass between each pair that touch, with no outlet at all:
    the label of each cell's catchment and the passes (geography.drainage's own pieces)."""
    from geography import drainage as dr
    rows, cols = h.shape
    idx2 = np.arange(rows * cols, dtype=np.int32).reshape(rows, cols)
    receiver = dr.d8_receivers(h, np.zeros(h.shape, bool), lat, cell_km)
    sinks = np.flatnonzero(receiver == np.arange(receiver.size))
    node_of = np.full(receiver.size, -1, np.int32)
    node_of[sinks] = np.arange(sinks.size, dtype=np.int32)
    label = node_of[dr.terminal(receiver)]
    edges = dr.passes(label.reshape(rows, cols), h.astype(np.float32), idx2)
    return label, sinks.size, edges


def joining_level(h_flat, label, n_nodes, edges, seed):
    """The level at which each cell joins the water standing over one seed cell: the higher of its own height and
    the minimax pass height between its catchment and the seed's (the flood run from the seed's catchment alone)."""
    from research.studies.drainage_network.routing import flood_from
    root = int(label[seed])
    level, _, _ = flood_from(n_nodes, edges, [(root, -np.inf)])
    joined = np.maximum(h_flat, level[label].astype(np.float32))
    joined[label == root] = h_flat[label == root]
    return joined


class Table:
    """A water's area, loss to the air and volume at any level: its cells sorted by the level at which they join."""

    def __init__(self, joined, area, loss, h_flat, cap):
        cells = np.flatnonzero(joined < cap)
        order = cells[np.argsort(joined[cells], kind='stable')]
        self.c = joined[order].astype(np.float64)
        self.area = np.cumsum(area[order])
        self.loss = np.cumsum(loss[order])
        self.h_area = np.cumsum(h_flat[order].astype(np.float64) * area[order])

    def at(self, level):
        """Area (m2), loss to the air beyond the rain on it (m3/s) and volume (m3) at a level."""
        i = int(np.searchsorted(self.c, level))
        if i == 0:
            return 0.0, 0.0, 0.0
        return float(self.area[i - 1]), float(self.loss[i - 1]), float(level * self.area[i - 1] - self.h_area[i - 1])

    def level_for(self, supply, below):
        """The level under `below` at which the loss first reaches the supply (m3/s): the level at which the cell
        whose loss takes it there joins, the water standing there over part of that cell. None if it does not."""
        n = int(np.searchsorted(self.c, below))
        hit = np.flatnonzero(self.loss[:n] >= supply)
        return float(self.c[int(hit[0])]) if hit.size else None


# ---------------------------------------------------------------- the network

def spill_targets(h, lat, cell_km, runoff, lake_net, seas, joined, seed, level0):
    """The sea each sea's overflow reaches: its spill path followed downhill through a routing in which the other
    seas take the rivers in at the common level and this one is a hollow like any other. Returns the targets and
    each path's cells from the sea's deepest point."""
    from geography import drainage as dr
    targets, paths = {}, {}
    for k in seas:
        regions = {j: joined[j] < level0 for j in seas if j != k}
        outlet = np.logical_or.reduce(list(regions.values()))
        r = dr.route(h, outlet.reshape(h.shape), lat, cell_km, runoff, lake_net)
        cell, path = seed[k], [seed[k]]
        seen = {cell}
        while not outlet[cell]:
            cell = int(r['receiver'][cell])
            if cell in seen:
                break
            seen.add(cell)
            path.append(cell)
        targets[k] = next((j for j, m in regions.items() if m[cell]), None)
        path = np.array(path, dtype=np.int64)
        # the overflow fills every hollow on its way to its spill level: the water surface along the path
        filled = np.maximum(h.ravel()[path], r['level'][r['label'][path]])
        paths[k] = (path, filled.astype(np.float32))
        del r
    return targets, paths


def route_at(h, lat, cell_km, runoff, lake_net, waters):
    """A routing in which each water, given as (mask over the flattened grid, level), takes the rivers in at its
    own level; returns the routing and the rivers (m3/s) each water receives."""
    from research.studies.drainage_network.routing import route_levels
    level = np.full(h.size, np.nan)
    for mask, lev in waters:
        level[mask] = lev
    outlet = np.isfinite(level)
    r = route_levels(h, outlet.reshape(h.shape), lat, cell_km, runoff, lake_net,
                     np.where(outlet, level, 0.0).reshape(h.shape))
    return r, [float(r['discharge'][mask].sum()) for mask, _ in waters]


def hollows_on_the_way(surface, ground, below):
    """The hollows an overflow fills on its way down a path: the runs of cells whose water surface stands more than
    half a metre over the ground, leaving out the sea below (`below`), where the path ends."""
    flooded = (surface > ground + 0.5) & ~below
    return int(np.count_nonzero(np.diff(np.concatenate([[0], flooded.astype(int)])) == 1))


def fill_history(table, supply, start, stop, steps=FILL_STEPS):
    """A water fed `supply` (m3/s) filled from level `start` towards `stop`, its loss to the air growing as it
    spreads: the years each step takes, and the level where it stops if its loss meets the supply first."""
    levels = np.linspace(start, stop, steps + 1)
    years, settles = 0.0, None
    v_prev = table.at(levels[0])[2]
    for a, b in zip(levels[:-1], levels[1:]):
        mid_loss = table.at(0.5 * (a + b))[1]
        if mid_loss >= supply:
            settles = float(a)
            break
        v = table.at(b)[2]
        years += (v - v_prev) / (supply - mid_loss) / YEAR_S
        v_prev = v
    return years, settles


def place(cell, cols, rotation=None):
    """A cell's latitude and longitude east on today's Moon, and on the turned Moon with a rotation."""
    lat, lon = 90.0 - (cell // cols + 0.5) / 16.0, ((cell % cols + 0.5) / 16.0 + 180.0) % 360.0 - 180.0
    out = dict(today=[round(float(lat), 2), round(float(lon), 2)])
    if rotation is not None:
        from climate.gcm.turned import to_turned
        t_lat, t_lon = to_turned(np.array([lat]), np.array([lon]), rotation)
        out['turned'] = [round(float(t_lat[0]), 2), round(float((t_lon[0] + 180.0) % 360.0 - 180.0), 2)]
    return out


def long_profile(path, surface, discharge, lat_rows, cols, cell_km, level0, lake):
    """A river's long profile along a path of cells, mouth or source first as given: its length, the fall of its
    water surface, its steep reaches (faster than STEEP_M_KM over one cell), its largest single fall between
    neighbouring cells, the lakes it passes through, and the profile every PROFILE_KM."""
    rows_, cols_ = path // cols, path % cols
    dy = np.diff(rows_) * cell_km
    dx = ((np.diff(cols_) + cols // 2) % cols - cols // 2) * cell_km * np.cos(np.radians(lat_rows[rows_[:-1]]))
    step = np.hypot(dx, dy)
    distance = np.concatenate([[0.0], np.cumsum(step)])
    z = surface[path].astype(float) - level0
    fall = -np.diff(z)                                          # positive where the water surface falls along the path
    steep = (np.abs(fall) / np.maximum(step, 1e-9)) > STEEP_M_KM
    marks = np.searchsorted(distance, np.arange(0.0, distance[-1] + PROFILE_KM, PROFILE_KM))
    marks = np.unique(np.clip(marks, 0, path.size - 1))
    lake_runs = np.flatnonzero(np.diff(np.concatenate([[0], lake[path].astype(int), [0]])) == 1)
    return dict(length_km=round(float(distance[-1])), surface_from_m=round(float(z[0])), surface_to_m=round(float(z[-1])),
                steep_km=round(float(step[steep].sum())), steep_reaches=int(np.count_nonzero(np.diff(np.concatenate([[0], steep.astype(int)])) == 1)),
                largest_fall_between_cells_m=round(float(np.abs(fall).max())) if fall.size else 0,
                lakes_passed=int(lake_runs.size),
                discharge_m3s=[round(float(discharge[path[0]])), round(float(discharge[path[-1]]))],
                profile=[[round(float(distance[i])), round(float(z[i]))] for i in marks])


def river_profiles(r, waters, h, lat, cell_km, level0, rotation, count=RIVERS_A_SEA, minimum=0.01):
    """The largest rivers into each water, traced upstream from their mouths along the larger tributary at each
    confluence until the river falls under a hundredth of its mouth's flow."""
    rows, cols = h.shape
    receiver, q = r['receiver'], r['discharge']
    cells = np.arange(receiver.size)
    donors = np.flatnonzero(receiver != cells)
    order = donors[np.argsort(q[donors], kind='stable')]
    upstream = np.full(receiver.size, -1, np.int64)
    upstream[receiver[order]] = order                           # each cell's largest donor (the last written)
    surface = np.where(r['lake'], np.nan_to_num(r['lake_level'], nan=0.0), h.ravel())
    out = {}
    for name, mask in waters.items():
        mouths = donors[mask[receiver[donors]] & ~mask[donors]]
        rivers = []
        for m in mouths[np.argsort(q[mouths])[::-1][:count]]:
            path, cell = [int(m)], int(m)
            while upstream[cell] >= 0 and q[upstream[cell]] >= minimum * q[m] and len(path) < 100000:
                cell = int(upstream[cell])
                path.append(cell)
            rivers.append(dict(mouth=place(int(m), cols, rotation),
                               **long_profile(np.array(path), surface, q, lat, cols, cell_km, level0, r['lake'])))
        out[name] = rivers
    return out


def basins_joined(joined_k, h, cell_area_row, level0, rim, count=JOINED_BASINS):
    """The largest basins a sea's water joins between the common level and its rim that lie below the common level
    themselves: each one's area below the common level, its floor and its centre on today's Moon."""
    from scipy import ndimage
    rows, cols = h.shape
    deep = ((joined_k >= level0 + 1.0) & (joined_k < rim) & (h.ravel() < level0)).reshape(rows, cols)
    lab, n = ndimage.label(deep, structure=np.ones((3, 3), bool))
    if not n:
        return []
    area = np.repeat(cell_area_row[:, None], cols, axis=1)
    sizes = ndimage.sum(area, lab, index=np.arange(1, n + 1))
    out = []
    for i in np.argsort(sizes)[::-1][:count]:
        rr, cc = np.nonzero(lab == i + 1)
        out.append(dict(area_below_the_common_level_km2=round(float(sizes[i]) / 1e6),
                        floor_m=round(float(h[rr, cc].min() - level0)),
                        joins_at_m=round(float(joined_k.reshape(rows, cols)[rr, cc].min() - level0)),
                        centre_today=[round(float(90.0 - (rr.mean() + 0.5) / 16.0), 1),
                                      round(float(((np.angle(np.exp(2j * np.pi * (cc + 0.5) / cols).mean()) / (2 * np.pi)
                                                    * 360.0) + 180.0) % 360.0 - 180.0), 1)]))
    return out


def geometry(log=print):
    """What the rain does not change: the heights, the six seas, where each one's water reaches at each level, the
    passes between them and where each one's overflow goes."""
    from research.studies.south_sea_ocean import inflows as I
    started = time.time()
    atlas, level0, h, lat, cell_km, runoff, lake_net, _, _ = load(1.0)
    rows, cols = h.shape
    h_flat = h.ravel()
    row_area = (cell_km * 1e3) ** 2 * np.cos(np.radians(lat))
    lon = (np.arange(cols) + 0.5) / 16.0
    body = I.label_waters(h < level0)
    named, _ = I.name_seas(body, row_area, lat, ((lon + 180.0) % 360.0) - 180.0, atlas['bodies'])
    seas = list(I.SEAS.values())
    label, n_nodes, edges = catchment_graph(h, lat, cell_km)
    joined, seeds = {}, {}
    for k in seas:
        cells = np.flatnonzero(body.ravel() == named[k]['label'])
        seeds[k] = int(cells[np.argmin(h_flat[cells])])
        joined[k] = joining_level(h_flat, label, n_nodes, edges, seeds[k])
    del label, edges
    log(f'{n_nodes:,} catchments and the seas\' joining levels, {time.time() - started:.0f} s')
    targets, spill_paths = spill_targets(h, lat, cell_km, runoff, lake_net, seas, joined, seeds, level0)
    log(f'spill targets {targets}, {time.time() - started:.0f} s')
    passes = {(a, b): float(joined[a][seeds[b]]) for a in seas for b in seas if a != b}
    return dict(atlas=atlas, level0=level0, h=h, lat=lat, cell_km=cell_km, seas=seas, joined=joined, seeds=seeds,
                targets=targets, spill_paths=spill_paths, passes=passes, catchments=n_nodes)


def network(rain, geo, log=print, evaporation=1.0, record=None):
    """The six seas' balances, rims, fills and chain on one rain, with the open water's evaporation scaled."""
    from geography import drainage as dr
    started = time.time()
    level0, h, lat, cell_km = geo['level0'], geo['h'], geo['lat'], geo['cell_km']
    seas, joined, seeds, targets, passes = geo['seas'], geo['joined'], geo['seeds'], geo['targets'], geo['passes']
    _, _, _, _, _, runoff, lake_net, precip, evap = load(rain, evaporation)
    rows, cols = h.shape
    h_flat = h.ravel()
    row_area = (cell_km * 1e3) ** 2 * np.cos(np.radians(lat))
    area = np.repeat(row_area, cols)
    loss = -lake_net.ravel().astype(np.float64) * dr.MM_DAY * area

    # every sea at the common level: the rivers each receives and the routing's budget
    at_common = [(joined[k] < level0, level0) for k in seas]
    r, rivers = route_at(h, lat, cell_km, runoff, lake_net, at_common)
    sea_cells = np.logical_or.reduce([m for m, _ in at_common])
    lake = r['lake'] & ~sea_cells
    land = ~sea_cells & ~r['lake']
    mm = dr.MM_DAY
    budget = dict(
        land_rain_less_evaporation_m3s=round(float(((precip - evap).ravel() * mm * area)[land].sum())),
        routing_runoff_m3s=round(float((runoff.ravel() * mm * area)[land].sum())),
        lakes_rain_less_evaporation_m3s=round(float((lake_net.ravel() * mm * area)[lake].sum())),
        seas_loss_m3s=round(float(loss[sea_cells].sum())),
        delivered_to_the_seas_m3s=round(sum(rivers)),
        lakes_km2=round(float(area[lake].sum()) / 1e6), closed_lakes_km2=round(float(area[lake & r['closed'][r['label']]].sum()) / 1e6))
    budget['surplus_m3s'] = budget['delivered_to_the_seas_m3s'] - budget['seas_loss_m3s']
    open_water_evaporation = float((evap_lake_of(lake_net, precip).ravel() * mm * area)[sea_cells | lake].sum())
    budget['open_water_evaporation_factor'] = evaporation
    budget['open_water_evaporation_m3s'] = round(open_water_evaporation)
    budget['open_water_evaporation_factor_closing_it'] = round(1.0 + budget['surplus_m3s'] / open_water_evaporation, 3)
    rotation = np.array(json.loads((ROOT / FRAME).read_text())['frame']['rotation_today_to_turned'])
    rivers_profiles = river_profiles(r, {k: m for k, (m, _) in zip(seas, at_common)}, h, lat, cell_km, level0, rotation)
    overflows = {}
    for k in seas:
        path, filled = geo['spill_paths'][k]
        rims_k = min(passes[(k, j)] for j in seas if j != k)
        leave = np.flatnonzero(joined[k][path] >= rims_k - BELOW_RIM_M)
        if leave.size:
            river = path[leave[0]:]
            surface = np.zeros(h_flat.size, np.float32)
            surface[river] = filled[leave[0]:]
            profile = long_profile(river, surface, r['discharge'], lat, cols, cell_km, level0, np.zeros(h_flat.size, bool))
            for key in ('discharge_m3s', 'lakes_passed'):   # the routing at the common level carries no overflow
                profile.pop(key)
            below = np.logical_or.reduce([joined[j][river] < level0 for j in seas if j != k])
            profile['hollows_filled_on_the_way'] = hollows_on_the_way(filled[leave[0]:], h_flat[river], below)
            overflows[k] = dict(leaves=place(int(river[0]), cols, rotation), enters=place(int(river[-1]), cols, rotation),
                                **profile)
    del r
    log(f'budget {budget}, {time.time() - started:.0f} s')

    rims = {k: min(passes[(k, j)] for j in seas if j != k) for k in seas}
    tables = {k: Table(joined[k], area, loss, h_flat, min(level0 + CAP_M, rims[k] + 1.0)) for k in seas}
    out = {}
    for i, k in enumerate(seas):
        rim_to = min((passes[(k, j)], j) for j in seas if j != k)
        rim = rim_to[0]
        tied = [j for j in seas if j != k and passes[(k, j)] <= rim + 1.0]
        target = targets[k] if targets[k] in tied else tied[0]
        a_rim, l_rim, v_rim = tables[k].at(rim - BELOW_RIM_M)
        a0, l0, v0 = tables[k].at(level0)
        waters = [(joined[k] < rim - BELOW_RIM_M, rim - BELOW_RIM_M)] + [
            (joined[j] < level0, level0) for j in seas if j != k and passes[(k, j)] > level0]
        r, got = route_at(h, lat, cell_km, runoff, lake_net, waters)
        supply_rim = got[0]
        del r
        years, settles = fill_history(tables[k], supply_rim, level0, rim - BELOW_RIM_M)
        out[k] = dict(
            atlas_area_km2=round(a0 / 1e6), atlas_volume_km3=round(v0 / 1e9),
            rivers_at_the_common_level_m3s=round(rivers[i]), loss_at_the_common_level_m3s=round(l0),
            rim_m=round(rim - level0, 1), overflow_to=target,
            area_at_the_rim_km2=round(a_rim / 1e6), volume_at_the_rim_km3=round(v_rim / 1e9),
            water_to_fill_km3=round((v_rim - v0) / 1e9),
            rivers_at_the_rim_m3s=round(supply_rim), loss_at_the_rim_m3s=round(l_rim),
            surplus_at_the_rim_m3s=round(supply_rim - l_rim),
            fills=settles is None and supply_rim > l_rim,
            years_to_fill=round(years) if settles is None and supply_rim > l_rim else None,
            settles_on_its_own_rivers_at_m=None if settles is None else round(settles - level0, 1),
            basins_joined_below_the_rim=basins_joined(joined[k], h, row_area, level0, rim - BELOW_RIM_M))
        # the level its own rivers at the common level hold it at, if they cannot fill it
        held = tables[k].level_for(rivers[i], rim - BELOW_RIM_M)
        out[k]['level_on_its_rivers_at_the_common_level_m'] = None if held is None else round(held - level0, 1)
        log(f'{k}: {out[k]}, {time.time() - started:.0f} s')
    # the chain: what each sea receives once the seas above it have filled
    order = sorted(seas, key=lambda k: (-round(rims[k], 1), -out[k]['surplus_at_the_rim_m3s']))
    received = {k: 0.0 for k in seas}
    for k in order:
        o = out[k]
        supply = o['rivers_at_the_rim_m3s'] + received[k]
        o['received_from_seas_above_m3s'] = round(received[k])
        held = tables[k].level_for(supply, rims[k] - BELOW_RIM_M)
        o['with_the_seas_above'] = dict(
            supply_m3s=round(supply), fills=held is None,
            settles_at_m=None if held is None else round(held - level0, 1),
            passes_on_m3s=round(supply - o['loss_at_the_rim_m3s']) if held is None else 0)
        target = o['overflow_to']
        if held is None and target and abs(rims[target] - rims[k]) < 1.0 and \
                out[target].get('with_the_seas_above', {}).get('fills'):
            o['with_the_seas_above']['joins'] = target          # the sea beyond is full at the same pass: they join
            o['with_the_seas_above']['passes_on_m3s'] = 0
        elif held is None and target:
            received[target] += supply - o['loss_at_the_rim_m3s']
    history = joined_history(seas, tables, rims, out, level0, record=record)
    # an overflow river runs only from a sea that fills on this rain, with the seas above it
    overflows = {k: v for k, v in overflows.items()
                 if out[k]['with_the_seas_above']['fills'] and 'joins' not in out[k]['with_the_seas_above']}
    return dict(rain=rain, evaporation=evaporation, level0=level0, budget=budget, seas=out, history=history,
                rivers=rivers_profiles,
                overflow_rivers=overflows, passes_m={f'{a}|{b}': round(v - level0, 1)
                                                                           for (a, b), v in passes.items() if a < b},
                seconds=round(time.time() - started))


def joined_history(seas, tables, rims, out, level0, years=HISTORY_YEARS, step=HISTORY_STEP_YEARS, record=None):
    """The six seas followed together from the atlas's state, all at the common level: each gains its own rivers and
    loses to the air at its level, and once full passes its surplus to the sea its overflow reaches. A sea's own
    rivers are those reaching it at the common level, rising to those reaching it at its rim as it rises. A sea
    that reaches the pass to a sea already standing there joins it, and the history stops following the two.
    Returns the events in order (a sea full, or joined) and each sea's level every 500 years, or at the years in
    `record`, with each sea's state there."""
    grid = {k: np.linspace(min(tables[k].c[0], level0), rims[k] - BELOW_RIM_M, 2000) for k in seas}
    curve = {k: np.array([tables[k].at(L)[2] for L in grid[k]]) for k in seas}
    volume = {k: tables[k].at(level0)[2] for k in seas}
    state = {k: 'filling' for k in seas}
    events, levels, t = [], [], 0.0
    to_level = lambda k, v: float(np.interp(v, curve[k], grid[k]))
    passed = {k: 0.0 for k in seas}
    while t < years:
        arriving, passed = passed, {k: 0.0 for k in seas}          # an overflow arrives a step after it leaves
        for k in seas:
            if state[k] == 'joined':
                continue
            lev = to_level(k, volume[k])
            share = np.clip((lev - level0) / max(rims[k] - level0, 1.0), 0.0, 1.0)
            rivers = (1.0 - share) * out[k]['rivers_at_the_common_level_m3s'] + share * out[k]['rivers_at_the_rim_m3s']
            supply = rivers + arriving[k]
            loss = tables[k].at(lev)[1]
            if state[k] == 'full':
                surplus = supply - tables[k].at(rims[k] - BELOW_RIM_M)[1]
                if surplus >= 0:
                    target = out[k]['overflow_to']
                    if state.get(target) == 'joined' or target is None:
                        continue
                    passed[target] += surplus
                    continue
                state[k] = 'filling'
            volume[k] = min(volume[k] + (supply - loss) * step * YEAR_S, curve[k][-1])
            volume[k] = max(volume[k], curve[k][0])
            if volume[k] >= curve[k][-1] and state[k] == 'filling':
                target = out[k]['overflow_to']
                if target in state and state[target] == 'full' and abs(rims[target] - rims[k]) < 1.0:
                    state[k] = state[target] = 'joined'
                    events.append(dict(year=round(t + step), event=f'{k} reaches its rim and joins {target}',
                                       level_m=round(rims[k] - level0, 1)))
                else:
                    state[k] = 'full'
                    events.append(dict(year=round(t + step), event=f'{k} is full and overflows to {target}',
                                       level_m=round(rims[k] - level0, 1)))
        t += step
        if (record is None and abs(t % 500.0) < 1e-6) or (record is not None and any(abs(t - y) < 1e-6 for y in record)):
            levels.append(dict(year=round(t), **{k: (None if state[k] == 'joined' else round(to_level(k, volume[k]) - level0, 1))
                                                 for k in seas}))
            if record is not None:
                levels[-1]['states'] = dict(state)
    return dict(events=events, levels_m=levels, years=years, step_years=step)


def evap_lake_of(lake_net, precip):
    """The open water's evaporation (mm/day) the routing uses: rain less its lake balance."""
    return precip - lake_net


# ---------------------------------------------------------------- the run

def load(rain, evaporation=1.0):
    from geography import drainage as dr
    from geography import topography as tp
    atlas = json.loads((ROOT / ATLAS).read_text())
    committed = json.loads((ROOT / DRAINAGE).read_text())
    clim = dict(np.load(ROOT / CLIMATOLOGY))
    meta = json.loads(str(clim.pop('metadata')))
    if meta.get('run') != RUN or list(meta.get('years')) != YEARS:
        raise SystemExit(f'{CLIMATOLOGY} is not {RUN} years {YEARS}')
    if committed['climatology']['run'] != RUN or list(committed['climatology']['years']) != YEARS:
        raise SystemExit(f'{DRAINAGE} was not routed on {RUN} years {YEARS}')
    level0 = float(atlas['sea_level_m'])
    h, _, _ = tp.height_above_geoid(dr.PPD, atlas['grid']['geoid_degree'])
    h = h.astype(np.float32)
    rows, cols = h.shape
    lat = 90.0 - (np.arange(rows) + 0.5) / dr.PPD
    cell_km = MOON_RADIUS / 1000.0 * np.radians(1.0 / dr.PPD)
    precip = (dr.climate_on_grid(clim, 'pr_mm_day', rows, cols) * np.float32(rain)).astype(np.float32)
    evap = dr.climate_on_grid(clim, 'evap_mm_day', rows, cols).astype(np.float32)
    runoff = np.maximum(precip - evap, 0.0).astype(np.float32)
    lake_net = (precip - np.float32(evaporation) * dr.lake_evaporation(clim, lat)[:, None].astype(np.float32)).astype(np.float32)
    return atlas, level0, h, lat, cell_km, runoff, lake_net, precip, evap


def build(rains=RAINS, log=print):
    """The design years' rain, then the same with the open water's evaporation raised until the routing's budget at
    the common level closes, then the rain range's ends."""
    started = time.time()
    geo = geometry(log)
    cases = [network(rains[0], geo, log)]
    factor = cases[0]['budget']['open_water_evaporation_factor_closing_it']
    cases.append(network(rains[0], geo, log, evaporation=factor))
    cases += [network(rain, geo, log) for rain in rains[1:]]
    return dict(
        schema=SCHEMA,
        producer=dict(files={p: sha256(ROOT / p) for p in CODE}),
        inputs={ATLAS: sha256(ROOT / ATLAS), DRAINAGE: sha256(ROOT / DRAINAGE), FRAME: sha256(ROOT / FRAME)},
        ignored_inputs={CLIMATOLOGY: sha256(ROOT / CLIMATOLOGY)},
        climatology=dict(run=RUN, years=YEARS),
        evidence=('A calculation on products: geography\'s steady annual routing of the design years\' runoff over '
                  'LOLA\'s 16 pixel/degree heights, with each of the six large seas taking rivers in at its own level '
                  'and filled towards its lowest pass to another sea. The climate is held at the design years\' as '
                  'the atlas\'s waters made it; the routing\'s budget differs from the climate model\'s, and the '
                  'model\'s rain is uncertain by the factor the climate lane carries, so the balances of seas near '
                  'balance are not settled by it. Not validated against any observation.'),
        reading_rule=('Levels in m above the atlas\'s common level; flows in m3/s; areas in km2 and volumes in km3. '
                      'A sea\'s rim is its lowest pass to another sea\'s water; at it the sea covers every basin its '
                      'water joins on the way up, listed as basins_joined_below_the_rim. rivers_at_the_rim is what '
                      'reaches the sea\'s water at its rim with the other seas at the common level; years_to_fill '
                      'fills it from the common level on those rivers with its loss to the air growing as it spreads. '
                      'with_the_seas_above adds what the seas above it pass on once full. budget sets the routing\'s '
                      'land runoff, which leaves out land where evaporation exceeds the rain, beside the climate '
                      'model\'s own rain less evaporation on the same land.'),
        catchments=geo['catchments'], spill_targets=geo['targets'],
        seeds_today={k: [round(float(90.0 - (c // geo['h'].shape[1] + 0.5) / 16.0), 2),
                         round(float(((c % geo['h'].shape[1] + 0.5) / 16.0 + 180.0) % 360.0 - 180.0), 2)]
                     for k, c in geo['seeds'].items()},
        cases=cases, seconds=round(time.time() - started))


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--rain', type=float, nargs='*', default=list(RAINS), help='factors on the design years\' rain')
    args = p.parse_args(argv)
    record = build(tuple(args.rain))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)} in {record["seconds"]} s')
    return 0


if __name__ == '__main__':
    sys.exit(main())
