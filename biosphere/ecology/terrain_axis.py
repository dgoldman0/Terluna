"""The terrain beside the moisture map on the turned Moon's design climate: the Codex mockups of 7 November tried on the
exact data (8 November 2026; on the corrected design climate from 9 November).

The moisture map classes the land by its water balance alone, and natural ground varies too little in temperature to
divide it further, so the biome analysis looks for a further axis among the conditions the world's builders will not set,
topography among them, and tries the mockups made on 7 November with the exact data where the data exist (register,
7 November). This module computes the terrain those mockups proposed on LOLA's heights above the GRAIL geoid at 16
pixels per degree (1.9 km), and joins each pixel to the moisture class of the climate cell that holds it on the turned
Moon:

- the terrain position of the humidity × terrain mockup: a pixel is a local low or a local high when its height differs
  from the area-weighted mean of the land within 30 km by more than max(100 m, half the land's standard deviation
  there), and intermediate otherwise;
- the six terrain settings of the landscape mockup, which drew them on an invented crater: basin floors (land at least
  10 m below the spill level of the closed depression that holds it, from the drainage routing), steep slopes (over
  10°), crests and rims (local highs), hollows and footslopes (local lows), plains (under 2°) and gentle slopes (the
  rest), each taking the first of those that applies;
- a topographic wetness index, ln(a / tan b), where a is the area drained through a pixel per unit width in the routing
  of the design climate's forty years and b its slope, floored at 0.1°;
- the wetland cues of that mockup's second panel: humid local lows (the humid and wet classes on local lows under 2°)
  and low ground near water (within 15 km of a sea or rain-fed lake, 0–100 m above its surface, under 2° and not a local
  high).

Slopes are centred differences one pixel away on four great-circle bearings, a 3.8 km baseline, with water surfaces in
place of submerged beds. The routing is recomputed from the forty-year climatology on today's ground and must give the
drainage product's lakes cell for cell. The climate is constant within each T21 cell, about 170 km across; land in a
cell the run treats as sea takes the class of the nearest land cell on the turned Moon, and the product counts that land
apart. The thresholds are the mockups' proposals and this module's first choices; the register has adopted none of them
as a biome boundary. The per-pixel fields go to the research drive for the sheets.

    python3 -m biosphere.ecology.terrain_axis [--fields PATH]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

from shared.constants import MOON_RADIUS
from climate.gcm.turned import read_frame, to_turned
from biosphere.ecology.moisture_map import cell_index

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.biosphere.ecology-terrain-axis/1'
MOISTURE = 'biosphere/ecology/results/moisture_map.json'
MOISTURE_SCHEMA = 'terluna.biosphere.ecology-moisture-map/1'
ORIENTATION = 'geology/results/moon_orientation.json'
DRAINAGE = 'geography/turned_corrected_y10-49'              # the lakes routed on the run's forty years
CLIMATOLOGY = 'climate/gcm/products/climatology_A28R_301k_design_l15_snow_s938_turned_refit_precession_y10-49_today.npz'
OUT = HERE / 'results' / 'terrain_axis.json'
FIELDS = Path('/media/projectspace/terluna-research/session-tools/terrain_axis/terrain_axis_16ppd.npz')
PPD = 16
RADIUS_KM, FLOOR_M, STD_FACTOR = 30.0, 100.0, 0.5           # the humidity × terrain mockup's terrain position
PLAIN_DEG, STEEP_DEG = 2.0, 10.0                            # slope classes at the 3.8 km baseline
BASIN_DEPTH_M = 10.0                                        # below the closed depression's spill level
TWI_SLOPE_FLOOR_DEG = 0.1
MARGIN_KM, MARGIN_RISE_M, CUE_SLOPE_DEG = 15.0, 100.0, 2.0  # the wetland mockup's water margin
HUMID = 4                                                   # the moisture map's humid and wet classes, aridity 0.65 and over
POSITIONS = ('local low', 'intermediate', 'local high')
SETTINGS = (dict(key='crest', name='crests and rims'), dict(key='steep', name='steep slopes'),
            dict(key='gentle', name='gentle slopes'), dict(key='plain', name='plains'),
            dict(key='hollow', name='hollows and footslopes'), dict(key='basin', name='basin floors'))
CUES = ('few cues', 'humid local low', 'low ground near water', 'both cues')
EVIDENCE = ('Geometric classes of measured terrain (LOLA above the GRAIL geoid, 16 px/deg) and a wetness index from '
            'the drainage routed on the design climate\'s forty years, joined to the moisture map\'s modelled classes '
            'cell by cell. The thresholds are proposals; no soil, groundwater, flooding or vegetation model enters, and '
            'nothing here is validated against ecology.')
READING_RULE = ('Shares are of the land at 16 px/deg (seas below the atlas\'s sea level and the rain-fed lakes of the '
                'forty-year routing excluded), by exact spherical area. moisture_by_position[m][p] and '
                'settings_by_moisture[m][s] index the moisture map\'s classes (hyper-arid to wet) and the positions or '
                'settings in the order listed. Land the run treats as sea at T21 carries the nearest land cell\'s class '
                'and is counted in land.spread_share; robust follows the moisture map\'s flag for the cell.')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rel(path):
    path = Path(path)
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def row_areas_km2(nr, nc):
    """Exact spherical area of one pixel in each row, north to south."""
    edges = np.radians(90.0 - np.arange(nr + 1) * 180.0 / nr)
    return (MOON_RADIUS / 1e3) ** 2 * (2.0 * math.pi / nc) * (np.sin(edges[:-1]) - np.sin(edges[1:]))


def unit(lat_deg, lon_deg):
    phi, lam = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack([np.cos(phi) * np.cos(lam), np.cos(phi) * np.sin(lam), np.sin(phi)], axis=-1)


def neighbour_mean_std(height, land, lat, radius_km):
    """Area-weighted mean and standard deviation of the land within a great-circle disk around every pixel.

    Row areas weight the pixels, longitude windows wrap, water is left out, and a pixel belongs to the disk when its
    centre does. From the humidity × terrain mockup."""
    nlat, nlon = height.shape
    dlon, dlat = 2 * np.pi / nlon, np.pi / nlat
    phi = np.radians(lat)
    sinp, cosp = np.sin(phi), np.cos(phi)
    area = np.sin(phi + dlat / 2) - np.sin(phi - dlat / 2)
    angle = radius_km * 1000.0 / MOON_RADIUS
    reach = int(np.ceil(angle / dlat)) + 1
    weights = land.astype(float) * area[:, None]
    fields = (weights, weights * height, weights * height ** 2)
    prefix = [np.pad(np.cumsum(np.tile(f, (1, 3)), axis=1), ((0, 0), (1, 0))) for f in fields]
    indices = np.arange(nlon) + nlon
    totals = np.zeros((3, nlat, nlon))
    for i in range(nlat):
        for j in range(max(0, i - reach), min(nlat, i + reach + 1)):
            if abs(phi[i] - phi[j]) > angle:
                continue
            cosine = (np.cos(angle) - sinp[i] * sinp[j]) / (cosp[i] * cosp[j])
            if cosine > 1:
                continue
            half = nlon if cosine <= -1 else int(np.floor(np.arccos(np.clip(cosine, -1, 1)) / dlon + 1e-10))
            if 2 * half + 1 >= nlon:
                for k, f in enumerate(fields):
                    totals[k, i] += f[j].sum()
            else:
                lo, hi = indices - half, indices + half + 1
                for k, pre in enumerate(prefix):
                    totals[k, i] += pre[j, hi] - pre[j, lo]
    del prefix
    mean = np.divide(totals[1], totals[0], out=np.full(height.shape, np.nan), where=totals[0] > 0)
    var = np.divide(totals[2], totals[0], out=np.zeros(height.shape), where=totals[0] > 0) - mean ** 2
    return mean, np.sqrt(np.maximum(var, 0.0))


def slope_degrees(height, lat, lon):
    """Slope from centred differences one grid step away on four great-circle bearings (from the mockup)."""
    nlat, nlon = height.shape
    delta = np.pi / nlat
    phi = np.radians(lat)[:, None]
    lam = np.radians(lon)[None, :]
    padded = np.concatenate([height[:, -1:], height, height[:, :1]], axis=1).astype(np.float64)
    samples = []
    for bearing in (0.0, np.pi, np.pi / 2, 3 * np.pi / 2):
        p2 = np.arcsin(np.clip(np.sin(phi) * np.cos(delta) + np.cos(phi) * np.sin(delta) * np.cos(bearing), -1, 1))
        l2 = lam + np.arctan2(np.sin(bearing) * np.sin(delta) * np.cos(phi), np.cos(delta) - np.sin(phi) * np.sin(p2))
        rows = np.broadcast_to((np.pi / 2 - p2) / delta - 0.5, height.shape)
        cols = (l2 % (2 * np.pi)) / (2 * np.pi / nlon) - 0.5 + 1
        samples.append(map_coordinates(padded, [rows, cols], order=1, mode='nearest', prefilter=False).astype(np.float32))
        del p2, l2, rows, cols
    north = (samples[0] - samples[1]) / (2 * MOON_RADIUS * delta)
    east = (samples[2] - samples[3]) / (2 * MOON_RADIUS * delta)
    return np.degrees(np.arctan(np.hypot(north, east))).astype(np.float32)


def positions(height, land, lat):
    mean, std = neighbour_mean_std(height, land, lat, RADIUS_KM)
    relative = height - mean
    bound = np.maximum(FLOOR_M, STD_FACTOR * std)
    pos = np.ones(height.shape, np.int8)
    pos[relative < -bound] = 0
    pos[relative > bound] = 2
    pos[~land] = -1
    return pos, relative.astype(np.float32), std.astype(np.float32)


def settings(pos, slope, depth, land):
    """The six settings, each pixel taking the first that applies: basin floor, steep, crest, hollow, plain, gentle."""
    s = np.full(pos.shape, 2, np.int8)
    s[slope < PLAIN_DEG] = 3
    s[pos == 0] = 4
    s[pos == 2] = 0
    s[slope > STEEP_DEG] = 1
    s[depth >= BASIN_DEPTH_M] = 5
    s[~land] = -1
    return s


def wetness(drained_m2, slope, px_area_m2):
    tan_b = np.maximum(np.tan(np.radians(slope.astype(np.float64))), math.tan(math.radians(TWI_SLOPE_FLOOR_DEG)))
    with np.errstate(divide='ignore'):                          # the seas drain nothing; they are not land
        return np.log(drained_m2 / np.sqrt(px_area_m2) / tan_b)


def moisture_on_pixels(product, rotation, lat, lon, land):
    """Each land pixel's moisture class and robust flag from the T21 cell that holds it on the turned Moon; land in the
    run's sea cells takes the nearest land cell's, flagged as spread."""
    g, c = product['grid'], product['cells']
    glat, glon = np.array(g['lat']), np.array(g['lon'])
    cls = np.full((glat.size, glon.size), -1, np.int8)
    rob = np.zeros(cls.shape, bool)
    cr, cc = np.array(c['row']), np.array(c['col'])
    cls[cr, cc] = c['moisture']
    rob[cr, cc] = np.array(c['robust'], bool)
    rows, cols = np.nonzero(land)
    tl, tn = to_turned(lat[rows], lon[cols], rotation)
    r, k = cell_index(tl, tn, glat, glon)
    m, rb = cls[r, k], rob[r, k]
    spread = m < 0
    tree = cKDTree(unit(glat[cr], glon[cc]))
    _, j = tree.query(unit(tl[spread], tn[spread]))
    m[spread] = np.array(c['moisture'], np.int8)[j]
    rb[spread] = np.array(c['robust'], bool)[j]
    out = np.full(land.shape, -1, np.int8)
    out[rows, cols] = m
    spread_px = np.zeros(land.shape, bool)
    spread_px[rows, cols] = spread
    robust_px = np.zeros(land.shape, bool)
    robust_px[rows, cols] = rb
    return out, spread_px, robust_px


def water_margin(land, water, height, surface, slope, pos, lat, lon):
    """Low ground near water: within MARGIN_KM of the nearest water pixel, 0 to MARGIN_RISE_M above its surface, under
    CUE_SLOPE_DEG and not a local high. Distances are great-circle, to the water pixels that touch land."""
    touch = np.zeros(water.shape, bool)
    for axis, shift in ((1, 1), (1, -1)):
        touch |= np.roll(land, shift, axis=axis)
    touch[1:] |= land[:-1]
    touch[:-1] |= land[1:]
    wr, wc = np.nonzero(water & touch)
    cand = land & (slope <= CUE_SLOPE_DEG) & (pos != 2)
    lr, lc = np.nonzero(cand)
    tree = cKDTree(unit(lat[wr], lon[wc]))
    chord = 2 * math.sin(MARGIN_KM * 1000.0 / MOON_RADIUS / 2)
    d, j = tree.query(unit(lat[lr], lon[lc]), distance_upper_bound=chord, workers=2)
    near = np.isfinite(d)
    rise = np.full(lr.size, np.nan)
    rise[near] = height[lr[near], lc[near]] - surface[wr[j[near]], wc[j[near]]]
    out = np.zeros(land.shape, bool)
    out[lr, lc] = near & (rise >= 0) & (rise <= MARGIN_RISE_M)
    return out


def weighted_quantiles(values, weights, qs):
    order = np.argsort(values)
    cum = np.cumsum(weights[order])
    return [float(values[order][np.searchsorted(cum, q * cum[-1])]) for q in qs]


def controls():
    """The mockup's independent checks, at a coarse grid: a constant surface, an analytic slope on the sphere, the
    longitude seam."""
    lat = np.linspace(89.375, -89.375, 144)
    lon = np.linspace(0.625, 359.375, 288)
    flat = np.full((144, 288), 432.0)
    every = np.ones(flat.shape, bool)
    m, s = neighbour_mean_std(flat, every, lat, 90.0)
    analytic = np.broadcast_to(1000.0 * np.sin(np.radians(lat))[:, None], flat.shape).copy()
    calc = slope_degrees(analytic, lat, lon).astype(float)
    expected = np.degrees(np.arctan(1000.0 / MOON_RADIUS * np.cos(np.radians(lat))))
    test = np.cos(np.radians(lon))[None, :] * np.cos(np.radians(lat))[:, None] * 1000.0
    m1, _ = neighbour_mean_std(test, every, lat, 120.0)
    m2, _ = neighbour_mean_std(np.roll(test, 21, axis=1), every, lat, 120.0)
    return dict(constant_terrain_mean_max_error_m=float(np.max(np.abs(m - 432.0))),
                constant_terrain_std_max_m=float(np.nanmax(s)),
                constant_terrain_slope_max_deg=float(np.max(slope_degrees(flat, lat, lon))),
                analytic_slope_max_error_deg=float(np.max(np.abs(calc[8:-8] - expected[8:-8, None]))),
                longitude_seam_max_error_m=float(np.nanmax(np.abs(m2 - np.roll(m1, 21, axis=1)))))


def route():
    """The drainage on the forty-year climatology, recomputed for its routing arrays; it must reproduce the drainage
    product's lakes."""
    from geography import drainage as gd
    src = json.loads((ROOT / DRAINAGE / 'results' / 'drainage.json').read_text())['sources']
    clim = ROOT / CLIMATOLOGY
    recorded = {Path(k).name: v for k, v in src.items()}
    if recorded.get(clim.name) != sha(clim):
        raise ValueError(f'{clim.name} differs from the climatology the drainage product routed')
    atlas, _, clim_meta, sea_level, h, lat, cell_km, _, r = gd.compute(climatology=clim)
    saved = np.load(ROOT / DRAINAGE / 'products' / 'drainage_28pct_16ppd.npz')
    lake = np.flatnonzero(r['lake'])
    same = bool(np.array_equal(lake, np.sort(saved['lake_cell'].astype(np.int64))))
    order = np.argsort(saved['lake_cell'])
    level_err = float(np.max(np.abs(r['lake_level'][lake] - saved['lake_level_m'][order]))) if same else math.inf
    if not same or level_err > 0.01:
        raise ValueError('the recomputed routing does not reproduce the drainage product\'s lakes')
    return dict(h=h, lat=lat, sea_level=float(sea_level), outlet=r['outlet'].reshape(h.shape),
                lake=r['lake'].reshape(h.shape), lake_level=r['lake_level'].reshape(h.shape),
                depth=np.where(r['label'] > 0, r['level'][r['label']] - h.ravel(), 0.0).reshape(h.shape),
                drained_m2=r['drained_m2'].reshape(h.shape), clim_meta=clim_meta,
                check=dict(lakes_reproduced=same, lake_cells=int(lake.size), lake_level_max_difference_m=level_err))


def build(fields_path=FIELDS):
    product = json.loads((ROOT / MOISTURE).read_text())
    if product['schema'] != MOISTURE_SCHEMA:
        raise ValueError(f'expected {MOISTURE_SCHEMA} in {MOISTURE}')
    frame = read_frame(ROOT / ORIENTATION)
    run = product['producer']
    rt = route()
    h, lat = rt['h'].astype(np.float32), rt['lat']
    nr, nc = h.shape
    lon = (np.arange(nc) + 0.5) / PPD
    water = rt['outlet'] | rt['lake']
    land = ~water
    surface = h.copy()
    surface[rt['outlet']] = rt['sea_level']
    surface[rt['lake']] = rt['lake_level'][rt['lake']]

    slope = slope_degrees(surface, lat, lon)
    pos, relative, std = positions(h.astype(np.float64), land, lat)
    setting = settings(pos, slope, rt['depth'], land)
    area_km2 = row_areas_km2(nr, nc)
    twi = wetness(rt['drained_m2'], slope, area_km2[:, None] * 1e6).astype(np.float32)
    twi[~land] = np.nan
    moisture, spread, robust = moisture_on_pixels(product, frame['rotation'], lat, lon, land)
    lake_class, _, _ = moisture_on_pixels(product, frame['rotation'], lat, lon, rt['lake'])
    humid_low = land & (moisture >= HUMID) & (pos == 0) & (slope <= CUE_SLOPE_DEG)
    margin = water_margin(land, water, h, surface, slope, pos, lat, lon)
    cue = (humid_low.astype(np.int8) + 2 * margin.astype(np.int8))
    cue[~land] = -1

    w = np.broadcast_to(area_km2[:, None], h.shape)
    wl = w[land]
    land_km2 = float(wl.sum())
    m_l, p_l, s_l, c_l = moisture[land].astype(int), pos[land].astype(int), setting[land].astype(int), cue[land].astype(int)
    joint = np.bincount(m_l * 3 + p_l, weights=wl, minlength=18).reshape(6, 3) / land_km2
    by_setting = np.bincount(m_l * 6 + s_l, weights=wl, minlength=36).reshape(6, 6) / land_km2
    cue_share = np.bincount(c_l, weights=wl, minlength=4) / land_km2
    land_by_m = np.bincount(m_l, weights=wl, minlength=6)
    lake_by_m = np.bincount(lake_class[rt['lake']].astype(int), weights=w[rt['lake']], minlength=6)
    t_l = twi[land].astype(np.float64)
    q10, q50, q90 = weighted_quantiles(t_l, wl, (0.1, 0.5, 0.9))
    twi_by_m = []
    for k in range(6):
        sel = m_l == k
        if not sel.any():
            twi_by_m.append(None)
            continue
        twi_by_m.append(dict(median=weighted_quantiles(t_l[sel], wl[sel], (0.5,))[0],
                             share_above_land_p90=float(wl[sel][t_l[sel] > q90].sum() / wl[sel].sum())))
    moisture_share = joint.sum(axis=1)
    within = (by_setting / np.where(moisture_share[:, None] > 0, moisture_share[:, None], 1.0))

    checks = controls()
    checks.update(rt['check'])
    checks.update(global_area_closure=float(w.sum() / (4 * math.pi * (MOON_RADIUS / 1e3) ** 2) - 1.0),
                  joint_partition_error=float(joint.sum() - 1.0), settings_partition_error=float(by_setting.sum() - 1.0),
                  cue_partition_error=float(cue_share.sum() - 1.0),
                  land_classes_complete=bool((moisture[land] >= 0).all() and (pos[land] >= 0).all() and (setting[land] >= 0).all()),
                  water_unclassed=bool((moisture[water] < 0).all() and (pos[water] < 0).all() and (setting[water] < 0).all()),
                  wetness_finite_on_land=bool(np.isfinite(t_l).all()))

    fields_path = Path(fields_path)
    fields_path.parent.mkdir(parents=True, exist_ok=True)
    meta = dict(schema=SCHEMA + '#fields', grid=dict(pixels_per_degree=PPD, rows=nr, cols=nc,
                reading='rows from 90 N, columns from 0 E on today\'s ground; carry to the turned Moon with '
                        'climate.gcm.turned.to_turned and the orientation product\'s rotation'),
                codes=dict(position=POSITIONS, setting=[s['name'] for s in SETTINGS], cue=CUES,
                           moisture=[x['name'] for x in product['classes']], water=-1))
    np.savez_compressed(fields_path, position=pos, setting=setting, moisture=moisture, spread=spread, robust=robust,
                        cue=cue, wetness=twi.astype(np.float16), slope_deg=slope.astype(np.float16),
                        relative_height_m=relative.astype(np.float16), metadata=json.dumps(meta))

    inputs = [MOISTURE, ORIENTATION, CLIMATOLOGY, f'{DRAINAGE}/results/drainage.json',
              f'{DRAINAGE}/products/drainage_28pct_16ppd.npz', 'geography/results/atlas.json',
              'geography/inputs/ldem_16.img', 'geography/inputs/ldem_4.img', 'geography/inputs/jggrx_0420a_sha.tab',
              'shared/constants.json']
    names = [x['name'] for x in product['classes']]
    return dict(
        schema=SCHEMA,
        producer=dict(domain='biosphere', module='biosphere/ecology/terrain_axis.py',
                      files={rel(p): sha(p)[:16] for p in (Path(__file__), HERE / 'moisture_map.py',
                                                           ROOT / 'geography' / 'drainage.py',
                                                           ROOT / 'geography' / 'topography.py',
                                                           ROOT / 'climate' / 'gcm' / 'turned.py')},
                      run=run['run'], years=run['years'], frame=dict(path=ORIENTATION, sha256=frame['sha256']),
                      inputs={p: sha(ROOT / p) for p in inputs}),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        thresholds=dict(neighbourhood_radius_km=RADIUS_KM, position_floor_m=FLOOR_M, position_std_factor=STD_FACTOR,
                        plain_below_deg=PLAIN_DEG, steep_above_deg=STEEP_DEG, basin_depth_below_spill_m=BASIN_DEPTH_M,
                        wetness_slope_floor_deg=TWI_SLOPE_FLOOR_DEG, margin_km=MARGIN_KM, margin_rise_m=MARGIN_RISE_M,
                        cue_slope_max_deg=CUE_SLOPE_DEG, humid_from_class=names[HUMID],
                        slope_baseline_km=round(2 * MOON_RADIUS * math.pi / nr / 1000, 2),
                        status='proposals from the 7 November mockups and first choices here; none adopted'),
        classes=dict(moisture=names, positions=list(POSITIONS), settings=list(SETTINGS), cues=list(CUES)),
        land=dict(area_km2=round(land_km2), share_of_moon=round(land_km2 / float(w.sum()), 4),
                  spread_share=round(float(wl[spread[land]].sum() / land_km2), 4),
                  robust_share=round(float(wl[robust[land]].sum() / land_km2), 4)),
        shares=dict(moisture=[round(float(x), 4) for x in moisture_share],
                    positions=[round(float(x), 4) for x in joint.sum(axis=0)],
                    moisture_by_position=np.round(joint, 4).tolist(),
                    settings=[round(float(x), 4) for x in by_setting.sum(axis=0)],
                    settings_by_moisture=np.round(by_setting, 4).tolist(),
                    settings_within_moisture=np.round(within, 4).tolist(),
                    cues=[round(float(x), 4) for x in cue_share]),
        lakes=dict(note='the forty years\' rain-fed lakes, each pixel classed by the cell that holds it on the turned '
                        'Moon as the land is; share_of_ground is lake over lake and land',
                   area_km2_by_moisture=[round(float(x)) for x in lake_by_m],
                   share_of_ground_by_moisture=[round(float(a / (a + b)), 4) if a + b > 0 else None
                                                for a, b in zip(lake_by_m, land_by_m)]),
        wetness=dict(index='ln(a / tan b): a the drained area per unit width (m), b the slope floored at 0.1°',
                     land_p10=round(q10, 2), land_median=round(q50, 2), land_p90=round(q90, 2),
                     by_moisture=[None if x is None else dict(median=round(x['median'], 2),
                                  share_above_land_p90=round(x['share_above_land_p90'], 4)) for x in twi_by_m]),
        checks=checks,
        fields=dict(path=str(fields_path), sha256=sha(fields_path)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('--fields', default=str(FIELDS))
    args = ap.parse_args(argv)
    product = build(args.fields)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    s, land = product['shares'], product['land']
    print(f"land {land['area_km2']:,} km2; spread {land['spread_share']:.1%}; positions "
          + ', '.join(f'{n} {x:.1%}' for n, x in zip(POSITIONS, s['positions'])) + '; settings '
          + ', '.join(f"{d['name']} {x:.1%}" for d, x in zip(SETTINGS, s['settings'])) + f' -> {rel(OUT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
