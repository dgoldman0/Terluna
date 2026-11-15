"""The south sea's circulation with salt: the ocean model's set-up on the sea, a measured trial and the run.

    /media/projectspace/terluna-research/ocean/run.sh -m research.studies.south_sea_ocean.circulation trial --step 1.0 --days 30
    /media/projectspace/terluna-research/ocean/run.sh -m research.studies.south_sea_ocean.circulation run --step 1.0 --years 50

The first two stages reasoned where the sea's brackish water would lie and whether its deeps would stay ventilated.
This runner sets MITgcm (climate/ocean/stratified.py) on the South Pole-Aitken Sea to compute both:

- the grid is spherical in a frame of its own, centred on the sea on the turned Moon, so that its cells keep one
  size from the sea's tropical north to its end at 82 S; the Coriolis parameter is the true latitude's;
- the floor is LOLA's at 16 pixels per degree under the atlas's sea level, averaged into the model's cells; a cell is
  water when the sea covers half of it;
- the surface takes the design run's year 10 as written every 6.5 hours (the record the spray reading used): its
  wind stress, turned into the grid's directions; its net heat flux with the sunlight inside it; and its
  evaporation less rain. The model's sea surface is also drawn towards that record's with a 30-day time scale, since
  the climate model's sea is a slab that carries no heat sideways. The year repeats;
- the rivers enter at their mouths as the six-seas routing of inflows.py places them, and the water the sea loses
  beyond them enters at the mouth of the Moscoviense route, so that the sea's level holds. The first stage found
  that this share of the Moscoviense system's water is what a held sea needs;
- the sea starts at rest at seawater's strength and 27 C throughout, the warmth of its coldest surface.

`trial` runs a few model days from rest, checks that the water stays sane and measures the machine time a model year
takes; `run` continues the experiment in chunks of a model year, each ending in a pickup; `read` reduces what has run.
Neither the monthly tide nor a tracer of the water's age is in this set-up yet.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from climate.gcm.turned import to_turned
from climate.ocean import lunar, mitgcm
from climate.ocean import stratified as ST
from research.studies.south_sea_ocean import inflows as I
from research.studies.south_sea_ocean import run as R
from research.studies.south_sea_ocean import spray as SP

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RUNS = ROOT / 'research' / 'runs' / 'south_sea_ocean'
OUT = HERE / 'results' / 'circulation_trial.json'
SCHEMA = 'terluna.research.south-sea-circulation-trial/1'
MOUTHS = RUNS / 'mouths_six_seas.npz'
PPD = 16
MIN_DEPTH_M = 10.0                      # a cell shallower than two top layers counts as land
THETA0_C = 27.0                         # the sea's coldest surface in the design years (80 S on the turned Moon)
SST_RELAX_DAYS = 30.0
YEAR_S = SP.WRITE_S * 1320              # the forcing's cycle, 357.5 days
MONTH_S = YEAR_S / 12.0                 # the output's interval, 29.8 days: a twelfth of the forcing's year
MOSCOVIENSE_MOUTH_TODAY = (-23.34, 176.09)     # where the Moscoviense route enters (inflows.json, two_seas)
MOSCOVIENSE_ARM_BOX = (-19.0, -9.0, 24.0, 34.0)  # the arm that route enters, in the grid's own degrees
OLD_WATER_DEPTHS_M = (0.0, 100.0, 500.0, 1000.0, 2000.0, 4000.0)
WARMTH_BANDS = ((14, 30), (30, 45), (45, 60), (60, 81))     # degrees south, true latitude on the turned Moon
EEL_WARM_LIMIT_C = 28.0                 # Japanese eel eggs hatch at 19-28 C (the register's row of 6 November)
EEL_DEPTH_M = 1050.0                    # 0.9-1.2 km holds the pressure of the eels' 150-200 m on Earth (seawater lakes)
PATCH_VISC_GRID = 0.0012                # 300 m2/s on a 30 km cell at a 900 s step, falling with the cell's area
SALT_CLASSES = ((0.0, 0.5, 'fresh'), (0.5, 5.0, 'nearly fresh'), (5.0, 18.0, 'brackish'), (18.0, 30.0, 'mixed'),
                (30.0, 34.0, 'diluted sea water'), (34.0, 36.0, 'sea water'), (36.0, 99.0, 'saltier than sea water'))


# ---------------------------------------------------------------- the sea on the model's grid

def frame_of_the_sea(rotation):
    inflows = json.loads((ROOT / R.INFLOWS).read_text())
    lat, lon = inflows['seas'][R.SEA_NAME]['centroid_today']
    t_lat, t_lon = to_turned(np.array([lat]), np.array([lon]), rotation)
    return ST.Frame(float(t_lat[0]), float(t_lon[0]))


def cells(glat, glon, xg, yg):
    """Row and column of grid places in a grid with edges xg, yg, and which fall inside it."""
    j = np.searchsorted(yg, glat) - 1
    i = np.searchsorted(xg, glon) - 1
    return j, i, (j >= 0) & (j < len(yg) - 1) & (i >= 0) & (i < len(xg) - 1)


def floor_on_grid(glat, glon, depth_m, area, step, margin=1, edges=None):
    """A grid of `step` degrees covering the places with one land cell round them (or the grid of the given edges),
    and each cell's depth: the area-weighted mean depth of the sea within it where the sea covers at least half of
    it, else land."""
    if edges is None:
        x0 = (np.floor(glon.min() / step) - margin) * step
        y0 = (np.floor(glat.min() / step) - margin) * step
        nx = int(np.ceil((glon.max() - x0) / step)) + margin
        ny = int(np.ceil((glat.max() - y0) / step)) + margin
        xg, yg = x0 + step * np.arange(nx + 1), y0 + step * np.arange(ny + 1)
    else:
        xg, yg = (np.asarray(e, float) for e in edges)
        nx, ny = len(xg) - 1, len(yg) - 1
    j, i, inside = cells(glat, glon, xg, yg)
    flat = j[inside] * nx + i[inside]
    wet = np.bincount(flat, weights=area[inside], minlength=ny * nx).reshape(ny, nx)
    volume = np.bincount(flat, weights=(area * depth_m)[inside], minlength=ny * nx).reshape(ny, nx)
    grid = mitgcm.Grid(xg=xg, yg=yg, depth=np.zeros((ny, nx)))
    share = wet / (grid.area / 1e6)
    mean = np.where(wet > 0, volume / np.maximum(wet, 1e-30), 0.0)
    grid.depth = np.where((share >= 0.5) & (mean >= MIN_DEPTH_M), mean, 0.0)
    return grid


def sea_floor(step, rotation, frame, edges=None):
    """The South Pole-Aitken Sea's floor on the model's grid, with the fine pixels' area and volume for comparison."""
    from geography import topography as tp
    atlas = json.loads((ROOT / I.ATLAS).read_text())
    level = float(atlas['sea_level_m'])
    h, _, _ = tp.height_above_geoid(PPD, atlas['grid']['geoid_degree'])
    rows, cols = h.shape
    lat = 90.0 - (np.arange(rows) + 0.5) / PPD
    lon = (np.arange(cols) + 0.5) / PPD
    r = lunar.RADIUS / 1000.0
    row_area = r * r * np.radians(1.0 / PPD) ** 2 * np.cos(np.radians(lat))
    body = I.label_waters(h < level)
    named, _ = I.name_seas(body, row_area * 1e6, lat, lon, atlas['bodies'])
    pr, pc = np.nonzero(body == named[R.SEA_NAME]['label'])
    depth = level - h[pr, pc]
    area = row_area[pr]
    t_lat, t_lon = to_turned(lat[pr], lon[pc], rotation)
    glat, glon = frame.to_grid(t_lat, t_lon)
    grid = floor_on_grid(glat, glon, depth, area, step, edges=edges)
    deepest = np.unravel_index(int(np.argmax(grid.depth)), grid.depth.shape)
    grid.depth = mitgcm.keep_connected(grid.depth, [deepest])
    grid = mitgcm.pad_even(grid)
    fine = dict(area_km2=float(area.sum()), volume_km3=float((area * depth).sum() / 1000.0),
                greatest_depth_m=float(depth.max()))
    return grid, fine


def coriolis(grid, frame):
    """2 Omega sin and cos of the true latitude at the cells' centres, and the sine at their south-west corners."""
    yc, xc = np.meshgrid(grid.yc, grid.xc, indexing='ij')
    yg, xg = np.meshgrid(grid.yg[:-1], grid.xg[:-1], indexing='ij')
    centre, _ = frame.to_true(yc, xc)
    corner, _ = frame.to_true(yg, xg)
    f = 2.0 * lunar.OMEGA
    return dict(centre=f * np.sin(np.radians(centre)), corner=f * np.sin(np.radians(corner)),
                cosine=f * np.cos(np.radians(centre)))


# ---------------------------------------------------------------- the record on the model's grid

def nearest_water(lat, lon, water):
    """For every cell of a latitude-longitude grid, the flat index of the nearest water cell (itself where it is
    water)."""
    la, lo = np.meshgrid(lat, lon, indexing='ij')
    pts = ST.unit(la.ravel(), lo.ravel())
    wet = np.flatnonzero(water.ravel())
    _, k = cKDTree(pts[wet]).query(pts)
    return wet[k]


def bilinear(lat, lon, t_lat, t_lon):
    """Indices and weights that carry a field on a grid of latitudes (north first) and evenly spaced longitudes to
    places: four flat indices and four weights for each place. Beyond the last row the nearest row serves."""
    nlat, nlon = len(lat), len(lon)
    south_first = lat[::-1]
    k = np.clip(np.searchsorted(south_first, t_lat) - 1, 0, nlat - 2)
    wy = np.clip((t_lat - south_first[k]) / (south_first[k + 1] - south_first[k]), 0.0, 1.0)
    i_south, i_north = nlat - 1 - k, nlat - 2 - k
    x = np.mod(t_lon - lon[0], 360.0) / (360.0 / nlon)
    j0 = np.floor(x).astype(int) % nlon
    wx = x - np.floor(x)
    j1 = (j0 + 1) % nlon
    index = np.stack((i_south * nlon + j0, i_south * nlon + j1, i_north * nlon + j0, i_north * nlon + j1))
    weight = np.stack(((1 - wy) * (1 - wx), (1 - wy) * wx, wy * (1 - wx), wy * wx))
    return index, weight


def forcing_on_grid(grid, frame, record):
    """The record's fluxes on the model's wet cells through its year (records, ny, nx): the stress turned into the
    grid's directions (Pa), the net heat flux and its sunlight (W/m2, upward), evaporation less rain (m/s, upward)
    and the surface temperature (C). Land cells of the record take their nearest water cell's values first."""
    yc, xc = np.meshgrid(grid.yc, grid.xc, indexing='ij')
    wet = grid.wet
    t_lat, t_lon = frame.to_true(yc[wet], xc[wet])
    fill = nearest_water(record['lat'], record['lon'], record['lsm'] < 0.5)
    index, weight = bilinear(record['lat'], record['lon'], t_lat, t_lon)

    def carry(values):
        v = values.reshape(values.shape[0], -1)[:, fill]
        return sum(v[:, index[k]] * weight[k] for k in range(4))

    east, north = carry(record['tauu']), carry(record['tauv'])
    u, v = frame.vector_to_grid(yc[wet], xc[wet], east, north)
    heat = record['rss'] + record['rls'] + record['hfss'] + record['hfls']        # into the sea
    parts = dict(ustress=u, vstress=v, hflux=-carry(heat), swflux=-carry(record['rss']),
                 sflux=carry(-record['evap'] - record['pr']), climsst=carry(record['ts']) - 273.15)
    out = {}
    for name, values in parts.items():
        full = np.zeros((values.shape[0],) + wet.shape)
        full[:, wet] = values
        out[name] = full
    return out


def read_record(path=SP.RECORD):
    import netCDF4
    names = ('tauu', 'tauv', 'rss', 'rls', 'hfss', 'hfls', 'evap', 'pr', 'ts')
    with netCDF4.Dataset(path) as nc:
        out = {k: np.array(nc[k][:], dtype=np.float64) for k in names}
        out.update(lat=np.array(nc['lat'][:], dtype=float), lon=np.array(nc['lon'][:], dtype=float),
                   lsm=np.array(nc['lsm'][0], dtype=float))
    return out


# ---------------------------------------------------------------- the rivers

def rivers_on_grid(grid, frame, rotation, loss_m3s, spread_km=0.0):
    """The rivers' inflow as a speed into each cell (m/s): the six-seas routing's mouths on the sea, each given to
    its own cell or the nearest wet one, and the rest of the sea's loss at the mouth of the Moscoviense route. With
    `spread_km`, a mouth's water is shared equally among the wet cells within that distance of it, so that a large
    river does not enter through one fine cell."""
    inflows = json.loads((ROOT / R.INFLOWS).read_text())
    want = inflows['rules']['six_seas']['seas'][R.SEA_NAME]['inflow_m3s']
    with np.load(MOUTHS) as z:
        cell, q, water = z['mouth_cell'].astype(np.int64), z['mouth_discharge_m3s'].astype(float), z['mouth_water']
    totals = {int(k): float(q[water == k].sum()) for k in np.unique(water)}
    label = min(totals, key=lambda k: abs(totals[k] - want))
    if abs(totals[label] - want) > 1.0:
        raise SystemExit(f'{MOUTHS} holds no water whose rivers sum to the product\'s {want} m3/s')
    sel = water == label
    pr, pc = np.divmod(cell[sel], 360 * PPD)
    lat = np.concatenate([90.0 - (pr + 0.5) / PPD, [MOSCOVIENSE_MOUTH_TODAY[0]]])
    lon = np.concatenate([(pc + 0.5) / PPD, [MOSCOVIENSE_MOUTH_TODAY[1]]])
    flow = np.concatenate([q[sel], [loss_m3s - totals[label]]])
    if flow[-1] < 0:
        raise SystemExit('The six-seas rivers exceed the sea\'s loss on this grid; no balance to add')
    glat, glon = frame.to_grid(*to_turned(lat, lon, rotation))
    yc, xc = np.meshgrid(grid.yc, grid.xc, indexing='ij')
    wet = np.flatnonzero(grid.wet.ravel())
    tree = cKDTree(ST.unit(yc.ravel()[wet], xc.ravel()[wet]))
    mouths = ST.unit(glat, glon)
    _, k = tree.query(mouths)
    m3s = np.zeros(grid.wet.size)
    near = tree.query_ball_point(mouths, 2.0 * np.sin(0.5 * spread_km * 1000.0 / lunar.RADIUS)) if spread_km > 0 else None
    for n in range(len(flow)):
        share = wet[near[n]] if near is not None and len(near[n]) else wet[k[n]:k[n] + 1]
        m3s[share] += flow[n] / len(share)
    m3s = m3s.reshape(grid.wet.shape)
    info = dict(six_seas_m3s=totals[label], moscoviense_route_m3s=float(flow[-1]), mouths=int(sel.sum()),
                cells_taking_rivers=int((m3s > 0).sum()), largest_cell_m3s=float(m3s.max()),
                largest_cell_m_day=float((m3s / grid.area).max() * 86400.0))
    return m3s / grid.area, info


# ---------------------------------------------------------------- the experiment

def build(step, deltat, tiles=(2, 2), name=None):
    rotation = np.array(json.loads((ROOT / R.FRAME).read_text())['frame']['rotation_today_to_turned'])
    frame = frame_of_the_sea(rotation)
    cache = RUNS / f'circulation_grid_{step:g}deg.npz'
    if cache.exists():
        with np.load(cache) as z:
            grid = mitgcm.Grid(xg=z['xg'], yg=z['yg'], depth=z['depth'])
            fine = json.loads(str(z['fine']))
    else:
        grid, fine = sea_floor(step, rotation, frame)
        RUNS.mkdir(parents=True, exist_ok=True)
        np.savez(cache, xg=grid.xg, yg=grid.yg, depth=grid.depth, fine=json.dumps(fine))
    forcing = forcing_on_grid(grid, frame, read_record())
    loss = float((forcing['sflux'].mean(axis=0) * grid.area).sum())
    runoff, rivers = rivers_on_grid(grid, frame, rotation, loss)
    delr = ST.levels(float(grid.depth.max()))
    salt0 = float(R.practical_from_reference(R.sc.SEAWATER_G_KG))
    rho = float(R.seawater_density(salt0, R.T_REF_C))
    sea = ST.Sea(name=name or f'south_sea_{step:g}deg', grid=grid, deltat=deltat, tiles=tiles, delr=delr,
                 coriolis=coriolis(grid, frame), theta0=THETA0_C, salt0=salt0, rho_const=rho, forcing=forcing,
                 runoff=runoff, forcing_period_s=SP.WRITE_S, sst_relax_s=SST_RELAX_DAYS * 86400.0,
                 viscosity_smag=2.2, output_s=MONTH_S, deep_output_s=YEAR_S)
    wet = grid.wet
    setting = dict(
        step_deg=step, frame_centre_turned=[frame.centre_lat, frame.centre_lon], nx=grid.nx, ny=grid.ny,
        levels=len(delr), top_layer_m=float(delr[0]), wet_columns=int(wet.sum()),
        wet_cells=int((np.cumsum(delr)[None, None, :] - delr[None, None, :] < grid.depth[:, :, None]).sum()),
        cell_km=[float(np.sqrt(grid.area[wet].min()) / 1e3), float(np.sqrt(grid.area[wet].max()) / 1e3)],
        area_km2=float(grid.area[wet].sum() / 1e6), volume_km3=float((grid.area * grid.depth).sum() / 1e9),
        greatest_depth_m=float(grid.depth.max()), lola_pixels=fine,
        true_latitudes=[float(x) for x in np.percentile(
            frame.to_true(*np.meshgrid(grid.yc, grid.xc, indexing='ij'))[0][wet], (0, 100))],
        forcing=dict(records=int(forcing['sflux'].shape[0]), interval_h=SP.WRITE_S / 3600.0,
                     loss_m3s=loss, loss_m_yr=loss / grid.area[wet].sum() * R.YEAR_S,
                     net_heat_w_m2=float(-(forcing['hflux'].mean(axis=0) * grid.area)[wet].sum()
                                         / grid.area[wet].sum()),
                     stress_pa_mean=float(np.hypot(forcing['ustress'], forcing['vstress'])[:, wet].mean()),
                     surface_c=[float(forcing['climsst'][:, wet].min()), float(forcing['climsst'][:, wet].max())]),
        rivers=rivers,
        water_left_over_m3s=float(((forcing['sflux'].mean(axis=0) - runoff) * grid.area).sum()),
        start=dict(theta_c=THETA0_C, practical_salinity=salt0, rho_kg_m3=rho), deltat_s=deltat)
    return sea, setting


def patch_edges(xg, yg, box, fine, ratio=1.15):
    """Edges of a grid that covers the coarse grid's span with cells of `fine` degrees inside box = (x0, x1, y0, y1)
    and cells growing by `ratio` back to the coarse spacing outside it."""
    x0, x1, y0, y1 = box
    return (mitgcm.stretched_edges(float(xg[0]), float(xg[-1]), (x0, x1), fine, float(xg[1] - xg[0]), ratio),
            mitgcm.stretched_edges(float(yg[0]), float(yg[-1]), (y0, y1), fine, float(yg[1] - yg[0]), ratio))


def carry_over(parent_grid, parent_share, parent_fields, grid):
    """A coarser run's fields on another grid's cells, level by level: each cell takes the nearest parent cell that
    holds water at that level; below the parent's deepest water the level above serves."""
    py, px = np.meshgrid(parent_grid.yc, parent_grid.xc, indexing='ij')
    cy, cx = np.meshgrid(grid.yc, grid.xc, indexing='ij')
    target = ST.unit(cy.ravel(), cx.ravel())
    out = {name: np.zeros((parent_share.shape[0],) + cy.shape) for name in parent_fields}
    for k in range(parent_share.shape[0]):
        open_ = parent_share[k] > 0
        if not open_.any():
            for name in out:
                out[name][k] = out[name][k - 1]
            continue
        _, n = cKDTree(ST.unit(py[open_], px[open_])).query(target)
        for name, values in parent_fields.items():
            out[name][k] = values[k][open_][n].reshape(cy.shape)
    return out


def build_patch(box, fine, deltat, parent_step=1.0, spread_km=6.0, name=None, tiles=(2, 2)):
    """The sea with a patch of fine cells round one shore, started from the coarse run's latest yearly means. The
    levels are the coarse run's; rivers are shared among the wet cells within `spread_km` of each mouth."""
    rotation = np.array(json.loads((ROOT / R.FRAME).read_text())['frame']['rotation_today_to_turned'])
    frame = frame_of_the_sea(rotation)
    parent_name = f'south_sea_{parent_step:g}deg'
    parent_folder = mitgcm.RUNS / 'mitgcm' / parent_name
    with np.load(RUNS / f'circulation_grid_{parent_step:g}deg.npz') as z:
        parent_grid = mitgcm.Grid(xg=z['xg'], yg=z['yg'], depth=z['depth'])
    delr = np.array(json.loads((parent_folder / 'experiment.json').read_text())['delr'])
    name = name or f'{parent_name}_patch_{fine:g}deg'
    cache = RUNS / f'circulation_grid_{name}.npz'
    edges = patch_edges(parent_grid.xg, parent_grid.yg, box, fine)
    if cache.exists() and np.array_equal(np.load(cache)['xg'], edges[0]) and np.array_equal(np.load(cache)['yg'], edges[1]):
        with np.load(cache) as z:
            grid = mitgcm.Grid(xg=z['xg'], yg=z['yg'], depth=z['depth'])
    else:
        grid, _ = sea_floor(parent_step, rotation, frame, edges=edges)
        np.savez(cache, xg=grid.xg, yg=grid.yg, depth=grid.depth)
    grid.depth = np.minimum(grid.depth, float(delr.sum()))          # the coarse run's levels hold the floor
    deep = sorted(parent_folder.glob('deep.*.data'))
    if not deep:
        raise SystemExit(f'{parent_folder} holds no yearly means to start from')
    means = mitgcm.read_mds(str(deep[-1])[:-5])
    share = mitgcm.read_mds(str(parent_folder / 'hFacC'))
    start = carry_over(parent_grid, share, dict(theta=means['THETA'], salt=means['SALT'], dye=means['TRAC01']), grid)
    forcing = forcing_on_grid(grid, frame, read_record())
    loss = float((forcing['sflux'].mean(axis=0) * grid.area).sum())
    runoff, rivers = rivers_on_grid(grid, frame, rotation, loss, spread_km=spread_km)
    salt0 = float(R.practical_from_reference(R.sc.SEAWATER_G_KG))
    sea = ST.Sea(name=name, grid=grid, deltat=deltat, tiles=tiles, delr=delr, coriolis=coriolis(grid, frame),
                 theta0=THETA0_C, salt0=salt0, rho_const=float(R.seawater_density(salt0, R.T_REF_C)), forcing=forcing,
                 runoff=runoff, forcing_period_s=SP.WRITE_S, sst_relax_s=SST_RELAX_DAYS * 86400.0,
                 viscosity_smag=2.2, visc_h=0.0, visc_grid=PATCH_VISC_GRID, start=start, output_s=MONTH_S,
                 deep_output_s=YEAR_S)
    wet = grid.wet
    in_box = wet & (grid.xc[None, :] > box[0]) & (grid.xc[None, :] < box[1]) & (grid.yc[:, None] > box[2]) \
        & (grid.yc[:, None] < box[3])
    setting = dict(
        parent=parent_name, started_from_year=int(deep[-1].name.split('.')[1]) * json.loads(
            (parent_folder / 'experiment.json').read_text())['deltat'] / YEAR_S,
        box_grid_degrees=list(box), fine_deg=fine, frame_centre_turned=[frame.centre_lat, frame.centre_lon],
        nx=grid.nx, ny=grid.ny, levels=len(delr), wet_columns=int(wet.sum()), wet_columns_in_the_box=int(in_box.sum()),
        cell_km=[float(np.sqrt(grid.area[wet].min()) / 1e3), float(np.sqrt(grid.area[wet].max()) / 1e3)],
        area_km2=float(grid.area[wet].sum() / 1e6), box_water_km2=float(grid.area[in_box].sum() / 1e6),
        volume_km3=float((grid.area * grid.depth).sum() / 1e9), rivers=rivers, river_spread_km=spread_km,
        forcing=dict(loss_m3s=loss), visc_grid=PATCH_VISC_GRID, deltat_s=deltat,
        water_left_over_m3s=float(((forcing['sflux'].mean(axis=0) - runoff) * grid.area).sum()))
    return sea, setting, in_box


def monitor(folder):
    """The last monitor block's numbers from the model's standard output."""
    import re
    out = {}
    text = (Path(folder) / 'STDOUT.0000').read_text(errors='replace')
    for key in ('time_secondsf', 'dynstat_eta_max', 'dynstat_eta_min', 'dynstat_uvel_max', 'dynstat_uvel_min',
                'dynstat_vvel_max', 'dynstat_vvel_min', 'dynstat_wvel_max', 'dynstat_theta_max', 'dynstat_theta_min',
                'dynstat_salt_max', 'dynstat_salt_min', 'dynstat_salt_mean', 'dynstat_theta_mean', 'advcfl_uvel_max',
                'advcfl_vvel_max', 'advcfl_wvel_max', 'advcfl_W_hf_max', 'ke_mean'):
        found = re.findall(rf'%MON {key}\s+=\s+(\S+)', text)
        if found:
            out[key] = float(found[-1])
    return out


def top_layer(sea):
    """The latest top-layer means: salinity, temperature and flow on the grid (None before any output)."""
    stems = sorted(sea.folder.glob('top.*.data'))
    if not stems:
        return None
    return mitgcm.read_mds(str(stems[-1])[:-5])


def deep_means(sea):
    """The latest whole-depth means (None before any output)."""
    stems = sorted(sea.folder.glob('deep.*.data'))
    return mitgcm.read_mds(str(stems[-1])[:-5]) if stems else None


def cell_volumes(sea):
    """Each model cell's volume (m3): the model's own open share of each cell (hFacC, written when it starts) where
    the run folder holds it, else the floor's. The model rounds thin bottom cells to rock or to a fifth of a layer,
    so its own shares are the ones its output must be weighed with."""
    if (sea.folder / 'hFacC.meta').exists():
        share = mitgcm.read_mds(str(sea.folder / 'hFacC'))
        return share * sea.delr[:, None, None] * sea.grid.area[None]
    top = np.cumsum(sea.delr) - sea.delr
    return np.clip(sea.grid.depth[None] - top[:, None, None], 0.0, sea.delr[:, None, None]) * sea.grid.area[None]


def old_water(sea, dye, volume=None):
    """The share of the water below each depth that has not been at the surface since the start, by the dye."""
    volume = cell_volumes(sea) if volume is None else volume
    top = np.cumsum(sea.delr) - sea.delr
    out = {}
    for depth in OLD_WATER_DEPTHS_M:
        below = top >= depth
        if volume[below].sum() > 0:
            out[f'below {depth:g} m'] = float((dye[below] * volume[below]).sum() / volume[below].sum())
    return out


def by_depth(values, volume):
    """A field's mean at each level over the sea's cells, None where a level holds no water."""
    total = volume.sum(axis=(1, 2))
    mean = (values * volume).sum(axis=(1, 2)) / np.maximum(total, 1e-30)
    return [float(m) if t > 0 else None for m, t in zip(mean, total)]


def warmth(sea, theta, volume, true_lat):
    """The year's temperature: its mean at each level, the top layer by band of true latitude (mean, coldest and
    warmest cell), the depth at which each band's water first stands at or under the eels' warm limit, and the
    water at the depth whose pressure matches the eels' spawning depth on Earth."""
    top = np.cumsum(sea.delr) - sea.delr
    mid = top + 0.5 * sea.delr
    open_ = volume > 0
    wet = open_[0]
    bands, cool = {}, (theta <= EEL_WARM_LIMIT_C) & open_
    first = np.where(cool.any(axis=0), mid[np.argmax(cool, axis=0)], np.nan)
    for a, b in WARMTH_BANDS:
        m = wet & (-true_lat >= a) & (-true_lat < b)
        if not m.any():
            continue
        t, w, d = theta[0][m], sea.grid.area[m], first[m]
        found = np.isfinite(d)
        bands[f'{a}-{b}S'] = dict(
            top_layer_c=dict(mean=float((t * w).sum() / w.sum()), coldest=float(t.min()), warmest=float(t.max())),
            at_or_under_the_limit_from_m=dict(median=float(np.median(d[found])) if found.any() else None,
                                              share_of_area_with_any=float(w[found].sum() / w.sum())))
    k = int(np.abs(mid - EEL_DEPTH_M).argmin())
    t = theta[k][open_[k]]
    return dict(by_depth_c=by_depth(theta, volume), bands=bands, warm_limit_c=EEL_WARM_LIMIT_C,
                at_the_eels_depth=dict(level_mid_m=float(mid[k]), coldest=float(t.min()), median=float(np.median(t)),
                                       warmest=float(t.max())))


def salinity_areas(sea, salt, where=None):
    """The sea's area (km2) by the top layer's salinity class, in g/kg; `where` narrows it to some cells."""
    g_kg = salt * R.sc.SEAWATER_G_KG / float(R.practical_from_reference(R.sc.SEAWATER_G_KG))
    wet, area = sea.grid.wet if where is None else where, sea.grid.area / 1e6
    return {name: float(area[wet & (g_kg >= lo) & (g_kg < hi)].sum()) for lo, hi, name in SALT_CLASSES}


def trial(args):
    sea, setting = build(args.step, args.deltat, name=f'south_sea_{args.step:g}deg_trial')
    exe = sea.binary()
    sea.output_s = sea.deep_output_s = args.days * 86400.0 / 2.0
    sea.write(fresh=True)
    t0 = time.time()
    sea.run(args.days * 86400.0, args.days * 86400.0)
    wall = time.time() - t0
    steps = int(round(args.days * 86400.0 / args.deltat))
    mon = monitor(sea.folder)
    top = top_layer(sea)
    deep = deep_means(sea)
    per_year = wall / steps * (YEAR_S / args.deltat) / 3600.0
    out = dict(schema=SCHEMA, setting=setting, build=exe.parent.name, model_days=args.days, steps=steps,
               wall_s=round(wall, 1), seconds_a_step=wall / steps, ranks=sea.tiles[0] * sea.tiles[1],
               hours_a_model_year=per_year, machine_load_at_the_end=float(open('/proc/loadavg').read().split()[0]),
               monitor=mon, old_water_share=None if deep is None else old_water(sea, deep['TRAC01']),
               top_layer=None if top is None else dict(
                   salinity_areas_km2=salinity_areas(sea, top['SALT']),
                   practical_salinity=[float(top['SALT'][sea.grid.wet].min()), float(top['SALT'][sea.grid.wet].max())],
                   theta_c=[float(top['THETA'][sea.grid.wet].min()), float(top['THETA'][sea.grid.wet].max())],
                   speed_m_s_max=float(np.hypot(top['UVEL'], top['VVEL'])[sea.grid.wet].max())))
    path = OUT.with_name(f'circulation_trial_{args.step:g}deg.json')
    path.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'setting'}, indent=1))
    print(f"{setting['wet_columns']} columns, {setting['wet_cells']} wet cells; {wall / steps:.3f} s a step; "
          f"{per_year:.2f} hours a model year on {out['ranks']} processes")
    return 0


def patch_trial(args):
    box = tuple(args.box)
    sea, setting, in_box = build_patch(box, args.fine, args.deltat, name=f'south_sea_1deg_patch_{args.fine:g}deg_trial')
    exe = sea.binary()
    sea.output_s = sea.deep_output_s = args.days * 86400.0 / 2.0
    sea.write(fresh=True)
    t0 = time.time()
    sea.run(args.days * 86400.0, args.days * 86400.0)
    wall = time.time() - t0
    steps = int(round(args.days * 86400.0 / args.deltat))
    top = top_layer(sea)
    out = dict(schema='terluna.research.south-sea-circulation-patch-trial/1', setting=setting, build=exe.parent.name,
               model_days=args.days, steps=steps, wall_s=round(wall, 1), seconds_a_step=wall / steps,
               ranks=sea.tiles[0] * sea.tiles[1], hours_a_model_year=wall / steps * (YEAR_S / args.deltat) / 3600.0,
               machine_load_at_the_end=float(open('/proc/loadavg').read().split()[0]), monitor=monitor(sea.folder),
               top_layer_in_the_box=dict(salinity_areas_km2=salinity_areas(sea, top['SALT'], in_box),
                                         speed_m_s_max=float(np.hypot(top['UVEL'], top['VVEL'])[in_box].max())))
    path = HERE / 'results' / f'circulation_patch_trial_{args.fine:g}deg.json'
    path.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({k: v for k, v in out.items() if k != 'monitor'}, indent=1))
    print(f"{setting['wet_columns']} columns ({setting['wet_columns_in_the_box']} in the box); {wall / steps:.3f} s a "
          f"step; {out['hours_a_model_year']:.2f} hours a model year on {out['ranks']} processes")
    return 0


def note_year(sea, done, end):
    """After each model year: one line of what the sea looks like, appended to the run's progress file."""
    top, deep = top_layer(sea), deep_means(sea)
    line = dict(year=end * sea.deltat / YEAR_S, wall=time.strftime('%Y-%m-%d %H:%M'), monitor=monitor(sea.folder),
                salinity_areas_km2=salinity_areas(sea, top['SALT']), old_water_share=old_water(sea, deep['TRAC01']))
    with open(sea.folder / 'progress.jsonl', 'a') as f:
        f.write(json.dumps(line) + '\n')


def run(args):
    sea, setting = build(args.step, args.deltat)
    if not sea.folder.exists() or not sea.pickups():
        sea.write(fresh=True)
        (sea.folder / 'setting.json').write_text(json.dumps(setting, indent=1) + '\n')
    sea.run(args.years * YEAR_S, YEAR_S, after_chunk=note_year)
    return 0


def read(args):
    """What has run so far: each model year's note and each month's top layer, as results/circulation_<step>deg.json."""
    sea, setting = build(args.step, args.deltat)
    volume = cell_volumes(sea)
    true_lat, _ = ST.Frame(*setting['frame_centre_turned']).to_true(*np.meshgrid(sea.grid.yc, sea.grid.xc, indexing='ij'))
    years = []
    for path in sorted(sea.folder.glob('deep.*.data')):
        deep = mitgcm.read_mds(str(path)[:-5])
        years.append(dict(year=int(path.name.split('.')[1]) * sea.deltat / YEAR_S,
                          old_water_share=old_water(sea, deep['TRAC01'], volume),
                          old_water_by_depth=by_depth(deep['TRAC01'], volume),
                          salinity_by_depth_g_kg=[None if v is None else v * R.sc.SEAWATER_G_KG
                                                  / float(R.practical_from_reference(R.sc.SEAWATER_G_KG))
                                                  for v in by_depth(deep['SALT'], volume)],
                          warmth=warmth(sea, deep['THETA'], volume, true_lat)))
    factor = R.sc.SEAWATER_G_KG / float(R.practical_from_reference(R.sc.SEAWATER_G_KG))
    wet, months = sea.grid.wet, []
    for path in sorted(sea.folder.glob('top.*.data')):
        top = mitgcm.read_mds(str(path)[:-5])
        iteration = int(path.name.split('.')[1])
        salt = top['SALT'][wet] * factor
        months.append(dict(year=iteration * sea.deltat / YEAR_S, salinity_areas_km2=salinity_areas(sea, top['SALT']),
                           salinity_g_kg=[float(salt.min()), float(salt.max())],
                           theta_c=[float(top['THETA'][wet].min()), float(top['THETA'][wet].max())]))
    out = dict(schema='terluna.research.south-sea-circulation/1',
               evidence=('MITgcm at lunar gravity on one-degree cells under one repeated model year of the fresh-sea '
                         'climate\'s fluxes; a first run from a uniform start, with no tide and with closure constants '
                         'fitted on Earth. Not validated against any observation.'),
               reading_rule=('months are means over a twelfth of the forcing\'s year for the top 5 m; areas in km2 by '
                             'salinity class in g/kg (fresh under 0.5, nearly fresh to 5, brackish to 18, mixed to '
                             '30, diluted sea water to 34, sea water to 36, saltier above). old_water_share is the '
                             'dye: the share of the water below a depth that has not been at the surface since the '
                             'start, as each model year\'s mean; years also give the dye, the salinity and '
                             'the temperature at each level (level_mid_m) over the model\'s own open cells, the '
                             'top layer\'s temperature by band of true latitude, the depth from which the water '
                             'stands at or under the eels\' warm limit, and the water at the eels\' depth. Water '
                             'below about 500 m still holds the run\'s starting temperature in its first decades.'),
               setting=setting, level_mid_m=[float(x) for x in np.cumsum(sea.delr) - 0.5 * sea.delr],
               model_years=len(years), years=years, months=months)
    path = HERE / 'results' / f'circulation_{args.step:g}deg.json'
    path.write_text(json.dumps(out, indent=1) + '\n')
    print(f'{len(years)} model years, {len(months)} months read into {path}')
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest='command', required=True)
    for name, fn in (('trial', trial), ('run', run), ('read', read), ('patch-trial', patch_trial)):
        s = sub.add_parser(name)
        s.add_argument('--box', type=float, nargs=4, default=MOSCOVIENSE_ARM_BOX,
                       help='the fine patch in the grid\'s own degrees: west east south north')
        s.add_argument('--fine', type=float, default=0.1, help='the patch\'s spacing in degrees')
        s.add_argument('--step', type=float, default=1.0, help='the grid\'s spacing in degrees (30 km a degree)')
        s.add_argument('--deltat', type=float, default=1800.0, help='the time step in seconds')
        s.add_argument('--days', type=float, default=30.0)
        s.add_argument('--years', type=float, default=1.0)
        s.set_defaults(fn=fn)
    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == '__main__':
    sys.exit(main())
