"""Where carriers set the rock down, on the design surface: does a rule of placement return each slope's rock to it?

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.rock_return.placement   # about five minutes

run.py's loop balances the whole: carriers that multiply where sediment stands return what arrives. It does not say
where the rock is set down. This module puts that question on the map. Every cell's rock is gathered where it first
settles (rock.py), and carriers working out from each such place set it down on the land within their range by a
rule:

- worn: in proportion to how fast each piece of ground is wearing, the rule of nesting on bare, freshly worn ground;
- worn_by_district: the same with the wear judged over a district some kilometres across and not cell by cell;
- crest: in proportion to how far the ground stands above its surroundings;
- even: evenly over the land in range.

The share a place's rock sends to a cell is the cell's weight over the summed weights of all the land in range of
that place, so the rock is conserved wherever any land is in range. The result is the rock each cell gets back
against the rock it lost: the share set down on the wrong cell, the same share judged over blocks of the map from 15
to 240 km (rock that is misplaced only within a block is no loss to that block), the share that returns to the
catchment it left, and how fast the ground left unbalanced rises or falls.

It also cuts regions of the design surface for evolve.py, which follows a region's ground through millions of years
under each rule. Writes results/placement.json and the regions to the research drive.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from research.studies.rock_return import rock as RK
from research.studies.sea_nutrients import loads as LO

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'placement.json'
SCHEMA = 'terluna.research.rock-return-placement/1'
REGIONS = Path('/media/projectspace/terluna-research/terrain-return/evolution/regions_16ppd.npz')
REGION_SCHEMA = 'terluna.research.rock-return-regions/1'
RULES = ('worn', 'worn_by_district', 'crest', 'even')
RANGES_KM = (10.0, 20.0, 50.0)       # how far from where the rock lies a carrier sets it down (cases; the mean way back is 9.4 km)
DISTRICT_KM = 10.0                   # the half-width over which the district rule judges wear
SURROUNDINGS_KM = 5.0                # the half-width of the ground a crest is judged against, the terrain study's relief window
BLOCKS = (8, 32, 128)                # cells to a block's side: about 15, 60 and 240 km
REGION_CELLS = 384                   # a region's side in cells, about 730 km from north to south
REGION_CENTRES = dict(korolev=(-4.0, 202.6), southern_highlands=(-45.0, 350.0), mendeleev=(5.7, 140.9), northern_far_side=(40.0, 100.0))


def sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def half_widths(range_km, lat_row, cell_km, cols):
    """A range as cells: the rows to each side, and for each row the columns to each side (more toward the poles)."""
    rows = max(1, int(round(range_km / cell_km)))
    wide = np.minimum(cols // 2 - 1, np.maximum(1, np.round(range_km / (cell_km * np.maximum(np.cos(np.radians(lat_row)), 1e-3))))).astype(np.int64)
    return rows, wide


def lon_sum(field, wide):
    """Each cell's sum of a field over the columns within its row's half-width, wrapping round the Moon."""
    out = np.empty_like(field, dtype=np.float64)
    cols = field.shape[1]
    for i in range(field.shape[0]):
        w = int(wide[i])
        s = np.cumsum(np.concatenate(([0.0], field[i, cols - w:], field[i], field[i, :w])))
        out[i] = s[2 * w + 1:] - s[:cols]
    return out


def lat_sum(field, rows_half):
    """Each cell's sum of a field over the rows within a half-width, stopping at the poles."""
    s = np.cumsum(np.vstack([np.zeros((1, field.shape[1])), field]), axis=0)
    n = field.shape[0]
    top = np.clip(np.arange(n) + rows_half + 1, 0, n)
    bottom = np.clip(np.arange(n) - rows_half, 0, n)
    return s[top] - s[bottom]


def gather(weight_mass, rows_half, wide):
    """For every cell taken as a place where rock lies: the summed weight of the land in range of it."""
    return lat_sum(lon_sum(weight_mass, wide), rows_half)


def spread(per_weight, rows_half, wide):
    """For every cell taken as ground: the sum, over the places it is in range of, of each place's rock per unit weight."""
    return lon_sum(lat_sum(per_weight, rows_half), wide)


def set_down(lying, weight_mass, rows_half, wide):
    """The rock each cell gets back, when rock lying at places is shared over the land in range by weight, and the
    rock that has no land in range. Arrays are grids; lying and the result are masses."""
    reach = gather(weight_mass, rows_half, wide)
    ok = (reach > 0) & (lying > 0)
    per = np.where(ok, lying / np.where(ok, reach, 1.0), 0.0)
    return weight_mass * spread(per, rows_half, wide), float(lying[~ok].sum())


def district_mean(value_mass, area_land, rows_half, wide):
    """A quantity per unit area averaged over the land of each cell's district."""
    a = spread(area_land, rows_half, wide)
    return np.where(a > 0, spread(value_mass, rows_half, wide) / np.where(a > 0, a, 1.0), 0.0)


def weights(rule, f, shape, lat_row, cell_km):
    """Each rule's weight per unit area on the eroding land."""
    eroding = f['eroding'].reshape(shape)
    area = f['area'].reshape(shape)
    wear = np.where(eroding, f['rock'].reshape(shape) / area, 0.0)
    if rule == 'worn':
        return wear
    if rule == 'worn_by_district':
        r, w = half_widths(DISTRICT_KM, lat_row, cell_km, shape[1])
        return np.where(eroding, district_mean(wear * area, np.where(eroding, area, 0.0), r, w), 0.0)
    if rule == 'crest':
        r, w = half_widths(SURROUNDINGS_KM, lat_row, cell_km, shape[1])
        surface = np.where(np.isfinite(f['water']), f['water'], f['h']).astype(np.float64).reshape(shape)
        around = spread(surface * area, r, w) / spread(area, r, w)
        return np.where(eroding, np.maximum(surface - around, 0.0), 0.0)
    if rule == 'even':
        return eroding.astype(np.float64)
    raise ValueError(rule)


def block_sum(field, b):
    rows, cols = field.shape
    return field[:rows - rows % b, :cols - cols % b].reshape(rows // b, b, cols // b, b).sum(axis=(1, 3))


def percentiles(values, weight, percents=(50, 90, 99)):
    order = np.argsort(values)
    cum = np.cumsum(weight[order]) / weight.sum()
    return {f'{p}%': float(np.interp(p / 100.0, cum, values[order])) for p in percents}


def judge(back, f, shape, rows_half, wide, lying_per_weight_at_first):
    """How well the rock set down matches the rock lost."""
    rock = f['rock'].reshape(shape)
    area = f['area'].reshape(shape)
    eroding = f['eroding'].reshape(shape)
    total = f['total']
    net = back - rock
    out = dict(returned_share=float(back.sum() / total),
               set_down_on_the_wrong_cell=float(np.abs(net).sum() / (2.0 * total)),
               wrong_at_blocks={f'{b} cells': float(np.abs(block_sum(net, b)).sum() / (2.0 * total)) for b in BLOCKS})
    rate = np.abs(net[eroding]) / (RK.ROCK_DENSITY * area[eroding]) * 1e6                   # m of rock a million years
    out['ground_left_unbalanced_m_per_myr'] = percentiles(rate, area[eroding])
    out['mean_unbalanced_over_mean_lowering'] = float((np.abs(net[eroding])).sum() / total)
    return out


def own_catchment(f, weight_mass, per, rows_half, wide, shape):
    """The share of the rock that carriers set down in the catchment it left: each cell's weight times the rock per
    unit weight at the place its own rock went to, where that place is in range of it."""
    cols = shape[1]
    cells = np.flatnonzero(f['eroding'])
    first = f['first'][cells]
    r1, c1 = np.divmod(cells, cols)
    r2, c2 = np.divmod(first, cols)
    dc = np.abs(c1 - c2)
    dc = np.minimum(dc, cols - dc)
    near = (np.abs(r1 - r2) <= rows_half) & (dc <= wide[r1])
    return float((weight_mass.ravel()[cells] * per.ravel()[first] * near).sum() / f['total'])


def cut_region(f, c, centre, shape):
    """A square of the design surface about a centre: the ground, the water over it, the rock each cell loses."""
    rows, cols = shape
    lat, lon = centre
    r0 = int(round((90.0 - lat) * 16.0)) - REGION_CELLS // 2
    c0 = int(round(lon * 16.0)) - REGION_CELLS // 2
    rr = np.arange(r0, r0 + REGION_CELLS)
    cc = np.arange(c0, c0 + REGION_CELLS) % cols
    take = lambda a: a.reshape(shape)[np.ix_(rr, cc)]                                   # noqa: E731
    area = take(f['area'])
    eroding = take(f['eroding'])
    lowering = np.where(eroding, take(f['rock']) / (RK.ROCK_DENSITY * area), 0.0)       # m of rock a year
    sea = take(c['sea']) >= 0
    water = take(f['water'])
    lake = np.isfinite(water) & ~sea
    return dict(ground_m=take(f['h']).astype(np.float32), water_m=water.astype(np.float32), sea=sea, eroding=eroding,
                lowering_m_yr=lowering.astype(np.float32), row0=r0, col0=c0, lat_centre=lat, lon_centre=lon,
                north_south_m=1000.0 * f['cell_km'], east_west_m=1000.0 * f['cell_km'] * float(np.cos(np.radians(lat))),
                stats=dict(land_share=float(eroding.mean()), sea_share=float(sea.mean()), lake_share=float(lake.mean()),
                           mean_lowering_m_per_myr=float(lowering[eroding].mean() * 1e6) if eroding.any() else 0.0,
                           ground_range_m=[float(take(f['h']).min()), float(take(f['h']).max())]))


def build(log=print):
    started = time.time()
    c = LO.prepare(log)
    f = RK.fields(c, log, started)
    shape = (f['rows'], f['cols'])
    lat_row, cell_km = f['lat_row'], f['cell_km']
    area = f['area'].reshape(shape)
    cells = np.flatnonzero(f['eroding'])
    lying = np.bincount(f['first'][cells], weights=f['rock'][cells], minlength=f['n']).reshape(shape)
    rules = {}
    for rule in RULES:
        w = weights(rule, f, shape, lat_row, cell_km) * area
        rules[rule] = {}
        for range_km in RANGES_KM:
            rows_half, wide = half_widths(range_km, lat_row, cell_km, shape[1])
            reach = gather(w, rows_half, wide)
            ok = (reach > 0) & (lying > 0)
            per = np.where(ok, lying / np.where(ok, reach, 1.0), 0.0)
            back = w * spread(per, rows_half, wide)
            entry = judge(back, f, shape, rows_half, wide, per)
            entry['rock_with_no_land_in_range'] = float(lying[~ok].sum() / f['total'])
            entry['set_down_in_the_catchment_it_left'] = own_catchment(f, w, per, rows_half, wide, shape)
            rules[rule][f'{range_km:g} km'] = entry
            log(f'{rule} at {range_km:g} km: {entry["set_down_on_the_wrong_cell"]:.3f} on the wrong cell, '
                f'{entry["set_down_in_the_catchment_it_left"]:.3f} in its own catchment, {time.time() - started:.0f} s')

    # the regions evolve.py follows
    regions, arrays = {}, {}
    for name, centre in REGION_CENTRES.items():
        r = cut_region(f, c, centre, shape)
        for key in ('ground_m', 'water_m', 'sea', 'eroding', 'lowering_m_yr'):
            arrays[f'{name}__{key}'] = r[key]
        regions[name] = {k: r[k] for k in ('row0', 'col0', 'lat_centre', 'lon_centre', 'north_south_m', 'east_west_m', 'stats')}
    meta = dict(schema=REGION_SCHEMA, cells=REGION_CELLS, regions=regions,
                reading_rule='For each region, name__ground_m is the redrawn floor (m above the dry Moon\'s geoid), '
                             'name__water_m the level of the lake or sea over a cell (nan on dry ground), name__sea and '
                             'name__eroding flags, and name__lowering_m_yr the rock each dry cell loses in the central '
                             'case (m of rock a year). Rows run from the north, columns to the east; a cell is '
                             'north_south_m by east_west_m at the region\'s centre.')
    REGIONS.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(REGIONS, metadata=np.array(json.dumps(meta)), **arrays)
    log(f'regions to {REGIONS}, {time.time() - started:.0f} s')

    return dict(
        schema=SCHEMA,
        producer=dict(module='research.studies.rock_return.placement',
                      code_sha256={p.name: sha256(p) for p in (Path(__file__), HERE / 'rock.py')}),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (RK.LAND_USE, RK.TERRAIN, ROOT / LO.ROUTING_RECORD)},
        product=dict(path=str(REGIONS), sha256=sha256(REGIONS), schema=REGION_SCHEMA),
        evidence=('Arithmetic on the rock module\'s fields: the rock each cell loses and the place it first settles on '
                  'the design surface. The rules of placement are design guesses for what a carrier might do; no '
                  'animal\'s choice of ground is behind them. A range is a square on the grid, and all the rock that '
                  'settles is taken as within the carriers\' reach. Not validated.'),
        reading_rule=('rules: for each rule and range, returned_share is the rock set down over the rock lost; '
                      'set_down_on_the_wrong_cell is half the summed difference between what cells get back and what '
                      'they lost, over the rock lost (0 when every cell gets back its own, 1 when none does); '
                      'wrong_at_blocks is the same with the map taken in blocks of that many cells to a side; '
                      'set_down_in_the_catchment_it_left is the share of the rock that carriers set down on land '
                      'that drains to the place they took it from; ground_left_unbalanced_m_per_myr gives the rate at '
                      'which the ground rises or falls where the return and the loss differ, by share of the land. '
                      'regions: the squares of the design surface written for evolve.py.'),
        setting=dict(rules=list(RULES), ranges_km=list(RANGES_KM), district_km=DISTRICT_KM, surroundings_km=SURROUNDINGS_KM,
                     blocks_cells=list(BLOCKS), block_km=[b * cell_km for b in BLOCKS], rock_kg_yr=f['total'],
                     mean_lowering_m_per_myr=f['total'] / RK.ROCK_DENSITY / f['land_m2'] * 1e6),
        rules=rules, regions=regions, seconds=round(time.time() - started))


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    out = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=1) + '\n')
    print(f'wrote {OUT} in {out["seconds"]} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
