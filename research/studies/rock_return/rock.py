"""The rock the land loses on the design surface: where it first settles, how far below its source, and in what water.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.rock_return.rock   # about four minutes

The terrain study (13 October 2026) worked out, on the earlier drainage, how fast the land's rock leaves the high
ground and what carrying it back would take. This module repeats the part of that reckoning a living return needs, on
the design surface the atlas round settled: the redrawn floor, the seas at their levels and the lakes of the design
routing. Each dry cell's rock is the terrain study's central lowering at its latitude, spread within the study's five
bands of latitude in proportion to local relief, as the study spreads it. It is followed down the design routing to the
first place it can settle (a sea, a lake that spills, a lake that passes nothing on, or a dry hollow where a course
ends), and the way back is measured: the distance over the ground and the height from that place to the cell the rock
left. With the lakes' own relation for particles it is carried on, to find what each lake and sea keeps and how long
each would take to fill. The land-use screen says whose land the rock leaves: wild, hybrid or people's.

The floors' depths are counted too, lake by lake and sea by sea, since what an animal can reach depends on them.
Writes results/rock.json, which run.py reads.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from research.studies.atlas_round import design_routing as DR
from research.studies.atlas_round import sand as SA
from research.studies.sea_nutrients import loads as LO
from research.studies.sea_nutrients import people as PE
from research.studies.terrain_through_time import routing as TR
from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'rock.json'
SCHEMA = 'terluna.research.rock-return-rock/1'
LAND_USE = ROOT / 'research' / 'studies' / 'atlas_round' / 'results' / 'land_use.json'
TERRAIN = ROOT / 'research' / 'studies' / 'terrain_through_time' / 'results' / 'terrain_through_time.json'
ROCK_DENSITY = SA.ROCK_DENSITY               # kg/m3, the terrain study's
PLACES = ('sea', 'lake_that_spills', 'lake_that_passes_nothing_on', 'dry_hollow')
KINDS = {1: 'wild', 2: 'hybrid', 3: 'peoples', 4: 'open_beyond_need'}
DEPTHS_M = (2.0, 20.0, 200.0)                # floor depths counted: marsh and wading, dabbling, a diving bird's reach
SEDIMENT_KG_M3 = (500.0, 1000.0, 1500.0)     # settled sediment, the terrain study's range
SPANS_YR = (1e3, 1e4, 1e5, 1e6, 1e7, 1e8)
FILL_SHARES = (0.001, 0.003, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99)
PERCENTS = (10, 50, 90)


def sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def band_weights(relief, area, eroding, abs_lat):
    """Each eroding cell's relief over the mean relief of the eroding land in its band of latitude (area-weighted), so
    that a band's mean lowering stays the terrain study's while its steep ground gives more."""
    w = np.zeros(relief.size)
    for lo, hi in SA.EROSION_BANDS:
        band = eroding & (abs_lat >= lo) & (abs_lat < hi if hi < 90.0 else abs_lat <= hi)
        if band.any():
            mean = float((relief[band] * area[band]).sum() / area[band].sum())
            w[band] = relief[band] / max(mean, 1e-9)
    return w


def first_settling(receiver, fronts, settles):
    """For every cell, the first cell on its way down where its rock can settle (itself, if it is one)."""
    first = np.arange(receiver.size, dtype=np.int64)
    for front in reversed(fronts):
        f = front[~settles[front]]
        first[f] = first[receiver[f]]
    return first


def arc_m(a, b, cols, ppd=16):
    """Great-circle distance (m) between cells a and b of the grid."""
    def latlon(cell):
        row, col = np.divmod(cell, cols)
        return np.radians(90.0 - (row + 0.5) / ppd), np.radians((col + 0.5) / ppd)
    la1, lo1 = latlon(a)
    la2, lo2 = latlon(b)
    s = np.sin(0.5 * (la2 - la1)) ** 2 + np.cos(la1) * np.cos(la2) * np.sin(0.5 * (lo2 - lo1)) ** 2
    return 2.0 * MOON_RADIUS * np.arcsin(np.sqrt(np.clip(s, 0.0, 1.0)))


def weighted(values, weights):
    """The weighted mean and the 10th, 50th and 90th weighted percentiles."""
    if weights.sum() <= 0:
        return dict(mean=None, percentiles=None)
    order = np.argsort(values)
    v, w = values[order], weights[order]
    cum = np.cumsum(w) / w.sum()
    return dict(mean=float((values * weights).sum() / weights.sum()),
                percentiles={f'{p}%': float(np.interp(p / 100.0, cum, v)) for p in PERCENTS})


def floor_by_depth(depth, area):
    """The area of a water's floor shallower than each depth, as shares of the whole."""
    total = float(area.sum())
    return {f'under {d:g} m': float(area[depth < d].sum() / total) for d in DEPTHS_M} if total > 0 else {}


def still_open(volume, kept_kg, area, ok):
    """The share of the lakes' area whose basins are not yet filled after each span, each lake a box that fills at
    what it keeps and passes nothing more as it fills, for each density of the settled sediment."""
    out = {}
    for rho in SEDIMENT_KG_M3:
        years = np.where(kept_kg > 0, volume * rho / np.maximum(kept_kg, 1e-30), np.inf)
        out[f'{rho:g} kg/m3'] = {f'{t:g}': float(area[ok & (years > t)].sum() / area[ok].sum()) for t in SPANS_YR}
    return out


def fill_quantiles(volume, kept_kg, area, ok, rho=SEDIMENT_KG_M3[1]):
    """The years by which each share of the lakes' area has filled, at the central density: pairs of share and years,
    from the lakes that fill first. A lake that keeps nothing never fills and is left out of the shares it would end."""
    years = np.where(kept_kg > 0, volume * rho / np.maximum(kept_kg, 1e-30), np.inf)[ok]
    a = area[ok]
    order = np.argsort(years)
    cum = np.cumsum(a[order]) / a.sum()
    y = years[order]
    out = {}
    for q in FILL_SHARES:
        k = int(np.searchsorted(cum, q))
        out[f'{q:g}'] = float(y[k]) if k < y.size and np.isfinite(y[k]) else None
    return out


def fields(c, log=print, started=None):
    """The fields the rock's account rests on, from the routing's set-up: each cell's area and ground, the water over
    it, whether it erodes, the rock it loses a year in the central case, and the first cell on its way down where
    that rock can settle. placement.py reads the same fields."""
    started = started or time.time()
    rows, cols, n = c['rows'], c['cols'], c['receiver'].size
    receiver, fronts, donor, body, keep, sea = (c[k] for k in ('receiver', 'fronts', 'donor', 'body', 'keep', 'sea'))
    names, levels = c['names'], c['levels']
    lat_row = 90.0 - (np.arange(rows) + 0.5) / 16.0
    cell_km = MOON_RADIUS / 1000.0 * np.radians(1.0 / 16.0)
    area = np.broadcast_to(c['row_area'][:, None], (rows, cols)).ravel()
    abs_lat = np.abs(np.broadcast_to(lat_row[:, None], (rows, cols)).ravel())
    h = c['h'].ravel()

    # the water's surface over every wet cell, for the relief and for the heights the rock falls through
    with np.load(DR.PRODUCT, allow_pickle=False) as z:
        lake_cell, lake_level = z['lake_cell'].astype(np.int64), z['lake_level_m'].astype(np.float32)
    water = np.full(n, np.nan, np.float32)
    water[lake_cell] = lake_level
    for k, name in enumerate(names):
        water[sea == k] = levels[name]
    surface = np.where(np.isfinite(water), water, h).astype(np.float32).reshape(rows, cols)
    relief = TR.local_relief(surface, lat_row, cell_km).ravel().astype(np.float64)
    del surface
    log(f'relief, {time.time() - started:.0f} s')

    eroding = donor & (body < 0)
    lowering = SA.erosion_m_per_yr(np.broadcast_to(lat_row[:, None], (rows, cols)).ravel()) * band_weights(relief, area, eroding, abs_lat)
    rock = np.where(eroding, lowering * area * ROCK_DENSITY, 0.0)                    # kg a year, the central case
    total = float(rock.sum())
    land_m2 = float(area[eroding].sum())
    log(f'{total / 1e12:.2f} Gt of rock a year off {land_m2 / 1e12:.2f} million km2, {time.time() - started:.0f} s')
    settles = ~donor | (body >= 0) | keep
    first = first_settling(receiver, fronts, settles)
    return dict(rows=rows, cols=cols, n=n, lat_row=lat_row, cell_km=cell_km, area=area, abs_lat=abs_lat, h=h, water=water,
                relief=relief, eroding=eroding, rock=rock, total=total, land_m2=land_m2, first=first,
                lake_cell=lake_cell, lake_level=lake_level)


def build(log=print):
    started = time.time()
    c = LO.prepare(log)
    f = fields(c, log, started)
    rows, cols, n, lat_row, area, h, water, relief, eroding, rock, total, land_m2, first, lake_cell, lake_level = (
        f[k] for k in ('rows', 'cols', 'n', 'lat_row', 'area', 'h', 'water', 'relief', 'eroding', 'rock', 'total', 'land_m2', 'first', 'lake_cell', 'lake_level'))
    receiver, fronts, donor, body, keep, sea = (c[k] for k in ('receiver', 'fronts', 'donor', 'body', 'keep', 'sea'))
    names, levels = c['names'], c['levels']

    # where it first settles, and the way back
    spills_at = np.zeros(n, bool)
    spills_at[body >= 0] = c['spills'][body[body >= 0]]
    place = np.where(~donor[first], 0, np.where(body[first] >= 0, np.where(spills_at[first], 1, 2), 3)).astype(np.int8)
    stop = np.where(np.isfinite(water[first]), water[first], h[first]).astype(np.float64)
    cells = np.flatnonzero(eroding)
    m = rock[cells]
    lift = np.maximum(h[cells].astype(np.float64) - stop[cells], 0.0)
    haul = arc_m(cells, first[cells], cols)
    pl = place[cells]
    del stop
    log(f'the way back, {time.time() - started:.0f} s')

    def ways(mask):
        w = np.where(mask, m, 0.0)
        return dict(rock_share=float(w.sum() / total), haul_m=weighted(haul, w), lift_m=weighted(lift, w),
                    glides_home_share={f'glide {g:g}': float(w[lift * g >= haul].sum() / max(w.sum(), 1e-30)) for g in (15.0, 25.0)})

    everywhere = np.ones(cells.size, bool)
    back = dict(all=ways(everywhere), **{PLACES[k]: ways(pl == k) for k in range(4)})

    # whose land it leaves
    use = json.loads(LAND_USE.read_text())
    grid = Path(use['product']['path'])
    if sha256(grid) != use['product']['sha256']:
        raise SystemExit(f'{grid} does not match the grid {LAND_USE} records')
    kinds, wild_source = {}, {}
    with np.load(grid, allow_pickle=False) as z:
        for key in PE.VARIANTS:
            kind = z['kind_' + key].ravel()
            kc = kind[cells]
            entry = {}
            for code, name in KINDS.items():
                on = kc == code
                entry[name] = dict(land_km2=float(area[cells][on].sum() / 1e6), **ways(on),
                                   by_place={PLACES[k]: float(m[on & (pl == k)].sum() / total) for k in range(4)})
            off = ~np.isin(kc, list(KINDS))
            entry['not_classed'] = dict(land_km2=float(area[cells][off].sum() / 1e6), rock_share=float(m[off].sum() / total))
            kinds[key] = dict(what=PE.VARIANTS[key], **entry)
            wild_source[key] = np.where(eroding & (kind == 1), rock, 0.0)
    log(f'kinds of land, {time.time() - started:.0f} s')

    # what each lake and sea keeps, with the lakes' relation for particles
    spills, lake_area, volume = c['spills'], c['area'], c['volume']
    into, mouth = c['into'], c['mouth']
    cases, kept_central = {}, None
    for case in ('low', 'central', 'high'):
        f = c['factors'][f'particles_{case}']
        p = LO.carry(receiver, fronts, rock, f)
        _, kept = PE.by_lake(c, p, f)
        to_sea = np.bincount(into[mouth], weights=p[mouth] * f[mouth], minlength=len(names))
        closed = float(p[keep].sum())
        cases[case] = dict(kept_by_lakes_that_spill=float(kept[spills].sum() / total), ends_where_water_ends=closed / total,
                           reaches_each_sea=dict(zip(names, (to_sea / total).tolist())),
                           reaches_the_seas=float(to_sea.sum() / total))
        if case == 'central':
            kept_central, sea_central, passing_central = kept, to_sea, p
    order = np.argsort(kept_central)[::-1]
    cum = np.cumsum(kept_central[order]) / kept_central.sum()
    lakes = dict(count=int(spills.sum()), area_km2=float(lake_area[spills].sum() / 1e6),
                 volume_km3=float(volume[spills].sum() / 1e9),
                 keep_half=int(np.searchsorted(cum, 0.5)) + 1, keep_nine_tenths=int(np.searchsorted(cum, 0.9)) + 1,
                 area_still_open_with_no_return=still_open(volume, kept_central, lake_area, spills),
                 years_to_fill_by_share_of_area=fill_quantiles(volume, kept_central, lake_area, spills))
    # the closed waters: what ends at each lowest point that passes nothing on, wet or dry
    ends = passing_central[keep]
    wet_end = body[np.flatnonzero(keep)] >= 0
    closed = dict(lowest_points=int(keep.sum()), under_a_lake=int(wet_end.sum()),
                  rock_share_under_a_lake=float(ends[wet_end].sum() / total), rock_share_on_dry_ground=float(ends[~wet_end].sum() / total))

    # the wild land's part of what each lake keeps
    wild = {}
    for key, source in wild_source.items():
        f = c['factors']['particles_central']
        _, kept_wild = PE.by_lake(c, LO.carry(receiver, fronts, source, f), f)
        share = np.where(kept_central > 0, kept_wild / np.maximum(kept_central, 1e-30), 0.0)
        ok = spills & (kept_central > 0)
        wild[key] = dict(of_what_the_lakes_keep=float(kept_wild[spills].sum() / kept_central[spills].sum()),
                         lake_area_share={'nine tenths wild or more': float(lake_area[ok & (share >= 0.9)].sum() / lake_area[ok].sum()),
                                          'a tenth wild or less': float(lake_area[ok & (share <= 0.1)].sum() / lake_area[ok].sum())},
                         kept_share={'in lakes nine tenths wild or more': float(kept_central[ok & (share >= 0.9)].sum() / kept_central[ok].sum()),
                                     'in lakes a tenth wild or less': float(kept_central[ok & (share <= 0.1)].sum() / kept_central[ok].sum())})
    log(f'lakes, {time.time() - started:.0f} s')

    # the floors' depths
    lake_depth = np.maximum(lake_level.astype(np.float64) - h[lake_cell], 0.0)
    on_spilling = c['spills'][body[lake_cell]]
    floors = dict(lakes_that_spill=floor_by_depth(lake_depth[on_spilling], area[lake_cell][on_spilling]),
                  lakes_that_pass_nothing_on=floor_by_depth(lake_depth[~on_spilling], area[lake_cell][~on_spilling]))
    seas = {}
    for k, name in enumerate(names):
        on = sea == k
        depth = np.maximum(levels[name] - h[on].astype(np.float64), 0.0)
        vol = float((depth * area[on]).sum())
        seas[name] = dict(area_km2=float(area[on].sum() / 1e6), volume_km3=vol / 1e9, floor=floor_by_depth(depth, area[on]),
                          rock_share_central=float(sea_central[k] / total),
                          years_to_fill_on_the_lands_rock={f'{rho:g} kg/m3': (vol * rho / sea_central[k] if sea_central[k] > 0 else None)
                                                           for rho in SEDIMENT_KG_M3})
    all_sea = sea >= 0
    sea_depth = np.maximum(np.array([levels[names[k]] for k in sea[all_sea]]) - h[all_sea].astype(np.float64), 0.0)
    floors['seas'] = floor_by_depth(sea_depth, area[all_sea])

    terrain = json.loads(TERRAIN.read_text())
    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.rock_return.rock', code_sha256=sha256(Path(__file__))),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (LAND_USE, TERRAIN, ROOT / LO.ROUTING_RECORD)},
        ignored_inputs={str(grid): use['product']['sha256']},
        evidence=('The terrain study\'s central lowering by latitude, spread by local relief within its bands, on the '
                  'design surface, and followed down the design routing with Earth\'s reservoir relation for the '
                  'particles a lake keeps at lunar settling. The lowering is a published range for Earth carried to the '
                  'Moon by the conservation study; no erosion has been modelled on this ground. Each lake is a box '
                  'that does not change as it fills. Not validated.'),
        reading_rule=('rock: kg a year off the dry land in the central case; the terrain study\'s slow and fast cases '
                      'are 0.2 and 2 times it. the_way_back: for the rock first settling in each kind of place, its '
                      'share of all the rock and the distance and height from that place to the cell it left, as means '
                      'and percentiles weighted by rock; glides_home_share is the rock for which a flyer released at '
                      'the source could glide back to the place at that glide ratio. kinds: the same by whose land the '
                      'rock leaves, for each wild share. carried_on: with the lakes passing on what Brune\'s curves '
                      'pass, the shares the lakes, the closed waters and the seas keep. lakes: the spilling lakes, how '
                      'few keep half and nine tenths of what the lakes keep, and the share of their area still open '
                      'after each span of years with nothing returned. floors: the share of each water\'s floor '
                      'shallower than each depth.'),
        setting=dict(rock_density_kg_m3=ROCK_DENSITY, central_lowering_m_per_myr=dict(zip(map(str, SA.EROSION_CENTRES_DEG), SA.EROSION_M_PER_MYR)),
                     slow_and_fast_against_central=list(SA.EROSION_RANGE), relief_half_width_km=TR.RELIEF_HALF_WIDTH_KM,
                     terrain_study=dict(mean_lift_m=terrain['upkeep']['mean_lift_m'], mean_distance_km=terrain['upkeep']['mean_distance_km'],
                                        where_the_rock_settles=terrain['upkeep']['where_the_rock_settles_at_the_start'],
                                        rock_gt_yr_central=terrain['upkeep']['rates']['50 m/Myr']['rock_gt_yr'])),
        rock=dict(kg_yr_central=total, land_km2=land_m2 / 1e6, mean_lowering_m_per_myr=total / ROCK_DENSITY / land_m2 * 1e6,
                  mean_relief_m=float((relief[eroding] * area[eroding]).sum() / land_m2)),
        the_way_back=back, kinds=kinds, carried_on=cases, lakes=lakes, closed_waters=closed, wild_part_of_the_lakes=wild,
        floors=floors, seas=seas, seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
