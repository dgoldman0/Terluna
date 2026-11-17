"""The moisture map on the design climate on the final surface: the run at sunlight 0.930, years 35-49, classed as
biosphere/ecology/moisture_map.py classes the design climate, with its settings taken from the atlas round's design
surface.

    ulimit -v 14000000; OPENBLAS_NUM_THREADS=1 nice -n 10 <climate venv python> -m research.studies.atlas_round.moisture_final

It takes about 10 seconds with the year sums drylands.py cached for the same run and years; without them it reads the
run's fifteen year files, about 3.4 GB.

The moisture map's own code runs unchanged; this module hands it two things the final surface differs in, so that its
product and those that read it stay as they are:

- the frame: the run on the final surface lives in the settled turned frame its surface record names
  (turned_surface_<frame>.json, rotation_today_to_turned and the tilt), where the moisture map reads geology's
  orientation product for the elastic frame;
- the settings: each T21 cell's seas, lakes, closed lakes, coast and tide come from the design surface the globe draws
  (globe_surface.py: every sea at its own level, the near-side sea parted, the governed lakes with their spills),
  where the moisture map reads the atlas's one sea level, its tides and geography's lakes. A 4 px/deg pixel is sea
  where most of it lies under one of the seven seas; the coast is the land-sea pixel edges in each land cell; the tide
  along it is the monthly range at that sea pixel from geography's tide model run on the design surface's own seas
  (tides_final.py, every sea at its own level, the forcing at the seas' turned places); the lakes are the surface's
  lakes off the seas at 16 px/deg, and the closed ones those that do not spill. No lake stands at the common level by
  rule, so the atlas's sea-level hollows are none.

The classes cover the run's T21 land cells, PlaSim's binary mask at its surface's own water share. Lakes smaller than
a cell sit on land cells, and to keep the water's share the mask makes water of some cells that are mostly land, so the
record's land.design_surface gives how the design surface's land and water fall on the run's land and water cells. The
design surface's land in the run's water cells takes the class of the nearest land cell (register, 17 November),
marked as assigned and counted apart (land.design_surface.assigned), and shares.over_the_design_land gives the classes'
shares over all the design surface's land, its land in the run's land cells taking its cell's class.

The window is the run's settled years alone: its earlier years cool toward the settled state, so no decade beside them
is mapped. The polar light's thresholds stay at the elastic frame's tilt (1.613 degrees) where the settled frame's is
1.707. The product is results/moisture_map_final.json, in the moisture map's schema, with its provenance put right.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import moisture_map as bm
from climate.gcm.turned import to_turned
from research.studies.atlas_round import drylands as DL
from research.studies.atlas_round import final_drainage as FD
from research.studies.atlas_round import globe_surface as GS
from research.studies.atlas_round import tides_final as TF
from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'moisture_map_final.json'
BULK = Path('/media/projectspace/terluna-research/atlas-round/moisture_final')
LAKES = BULK / 'design_lakes'                  # the design surface's lakes in geography's drainage product's form
CODE = ('research/studies/atlas_round/moisture_final.py', 'research/studies/atlas_round/globe_surface.py',
        'research/studies/atlas_round/tides_final.py', 'biosphere/ecology/moisture_map.py')
RUN, YEARS = FD.RUN, (35, 49)
SEAS = range(1, 8)
LAKE = 8


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frame_record():
    path, record = FD.frame_record()
    return path, record


def settled_frame(_path=None):
    """The moisture map's frame (climate.gcm.turned.read_frame's form) from the final surface's record."""
    path, record = frame_record()
    return dict(rotation=np.asarray(record['frame']['rotation_today_to_turned'], float),
                tilt_deg=float(record['frame']['tilt_to_ecliptic_deg']), frame=record['frame'], sha256=sha256(path))


def surface_grid():
    """The design surface's water at 16 px/deg on today's grid: each cell's body (0 dry, 1-7 a sea, 8 a lake) and
    whether its water spills."""
    with np.load(GS.PRODUCT, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        rows, cols = meta['grid']['rows'], meta['grid']['cols']
        body = np.zeros(rows * cols, np.int8)
        body[z['water_cell']] = np.where(z['water_body'] > 0, z['water_body'], LAKE)
        spills = np.zeros(rows * cols, bool)
        spills[z['water_cell']] = z['water_spills']
        level = np.full(rows * cols, np.nan, np.float32)
        level[z['water_cell']] = z['water_level_m']
        channel_cell, channel_q = z['channel_cell'], z['channel_discharge_m3s']
    return dict(meta=meta, rows=rows, cols=cols, body=body, spills=spills, level=level,
                channel_cell=channel_cell, channel_q=channel_q)


def map_grid(surface, k=4):
    """The 4 px/deg pixels: sea where most of a pixel lies under one of the seas, and that sea."""
    body = surface['body'].reshape(surface['rows'], surface['cols'])
    r, c = surface['rows'] // k, surface['cols'] // k
    b = body.reshape(r, k, c, k)
    counts = np.stack([(b == s).sum(axis=(1, 3)) for s in SEAS])
    sea_share = counts.sum(axis=0) / (k * k)
    sea_body = np.where(sea_share >= 0.5, np.argmax(counts, axis=0) + 1, 0).astype(np.int8)
    return sea_body


def surface_settings(lat, lon, rotation=None, drainage=None):
    """Seas, lakes, closed lakes, the coast's length and the tide along it for each T21 cell, from the design surface
    (moisture_map.settings's keys and units), each pixel counted in the cell that holds it on the turned Moon."""
    shape = (lat.size, lon.size)
    r_km = MOON_RADIUS / 1e3
    surface = surface_grid()
    sea_body = map_grid(surface)
    sea = sea_body > 0
    landpx = ~sea
    t = np.load(ROOT / bm.GEOGRAPHY['tides_grid'])
    lat4, lon4 = t['lat_deg'], t['lon_deg']
    if sea.shape != (lat4.size, lon4.size):
        raise ValueError('the design surface\'s map grid must be the tide grid')
    step = math.radians(0.25)
    area4 = (r_km * step) ** 2 * np.cos(np.radians(lat4))[:, None] * np.ones((1, lon4.size))
    plat, plon = lat4[:, None] * np.ones((1, lon4.size)), np.ones((lat4.size, 1)) * lon4[None, :]
    if rotation is not None:
        plat, plon = to_turned(plat, plon, rotation)
    rr, cc = bm.cell_index(plat, plon, lat, lon)
    flat = (rr * lon.size + cc).ravel()

    def total(weights):
        return np.bincount(flat, weights=np.asarray(weights, float).ravel(), minlength=lat.size * lon.size).reshape(shape)

    cell_area4 = total(area4)
    out = dict(sea_fraction=total(area4 * sea) / cell_area4, sea_level_lake_fraction=np.zeros(shape),
               atlas_land_fraction=total(area4 * landpx) / cell_area4)
    tide = np.full(sea.shape, np.nan)
    tide[t['rows'], t['cols']] = t['monthly_range_median_m']
    tide[~sea] = np.nan
    coast = np.zeros(sea.shape)
    tide_len = np.zeros(sea.shape)
    tide_w = np.zeros(sea.shape)
    ew = r_km * step
    for shift in (1, -1):
        nb_sea = np.roll(sea, -shift, axis=1)
        nb_tide = np.roll(tide, -shift, axis=1)
        edge = landpx & nb_sea
        coast += edge * ew
        tide_len += np.where(edge & np.isfinite(nb_tide), nb_tide * ew, 0.0)
        tide_w += np.where(edge & np.isfinite(nb_tide), ew, 0.0)
    for shift in (1, -1):
        nb_sea = np.zeros_like(sea)
        nb_tide = np.full(tide.shape, np.nan)
        lat_edge = lat4.copy()
        if shift == 1:
            nb_sea[:-1], nb_tide[:-1] = sea[1:], tide[1:]
            lat_edge[:-1] = 0.5 * (lat4[:-1] + lat4[1:])
        else:
            nb_sea[1:], nb_tide[1:] = sea[:-1], tide[:-1]
            lat_edge[1:] = 0.5 * (lat4[1:] + lat4[:-1])
        length = (r_km * step * np.cos(np.radians(lat_edge)))[:, None] * np.ones((1, lon4.size))
        edge = landpx & nb_sea
        coast += edge * length
        tide_len += np.where(edge & np.isfinite(nb_tide), nb_tide * length, 0.0)
        tide_w += np.where(edge & np.isfinite(nb_tide), length, 0.0)
    out['coast_km'] = total(coast)
    with np.errstate(invalid='ignore', divide='ignore'):
        with_tide = total(tide_w)
        out['tide_m'] = np.where(with_tide > 0, total(tide_len) / with_tide, np.nan)
    cells = np.flatnonzero(surface['body'] == LAKE)
    closed_cell = ~surface['spills'][cells]
    lat16 = 90.0 - ((cells // surface['cols']) + 0.5) / 16.0
    lon16 = ((cells % surface['cols']) + 0.5) / 16.0
    a16 = (r_km * math.radians(1.0 / 16.0)) ** 2 * np.cos(np.radians(lat16))
    if rotation is not None:
        lat16, lon16 = to_turned(lat16, lon16, rotation)
    r16, c16 = bm.cell_index(lat16, lon16, lat, lon)
    f16 = r16 * lon.size + c16
    box = bm.cell_areas_km2(lat, lon)
    lakes = np.bincount(f16, weights=a16, minlength=lat.size * lon.size).reshape(shape)
    closed = np.bincount(f16, weights=a16 * closed_cell, minlength=lat.size * lon.size).reshape(shape)
    out.update(rainfed_lake_fraction=lakes / box, closed_lake_fraction=closed / box, closed_lake_km2=closed,
               sea_level_lake_km2=np.zeros(shape))
    return out


def coverage(lat, lon, geo, rows, cols):
    """How the design surface's land and water fall on the run's land cells (the classes' cells, by row and col) and
    its water cells, in km2."""
    area = bm.cell_areas_km2(lat, lon)
    land = np.zeros(area.shape, bool)
    land[rows, cols] = True
    sea = geo['sea_fraction'] * area
    lake = geo['rainfed_lake_fraction'] * area
    dry = area - np.minimum(sea + lake, area)
    return dict(
        run_land_km2=round(float(area[land].sum())),
        design_land_km2=round(float(dry.sum())),
        design_land_in_run_land_km2=round(float(dry[land].sum())),
        design_lake_in_run_land_km2=round(float(lake[land].sum())),
        design_sea_in_run_land_km2=round(float(sea[land].sum())),
        design_land_in_run_water_km2=round(float(dry[~land].sum())),
        design_land_unclassed_share=round(float(dry[~land].sum() / dry.sum()), 4),
        run_land_cells_mostly_water=int((land & (dry < 0.5 * area)).sum()),
        run_water_cells_mostly_land=int((~land & (dry > 0.5 * area)).sum()),
        note=('The run\'s land cells are PlaSim\'s binary mask at its surface\'s own water share (climate/gcm/boundary.py). '
              'The design surface\'s seas are its 4 px/deg sea pixels and its lakes its 16 px/deg lake cells, each counted '
              'in the T21 cell that holds it on the turned Moon; its land is the rest of each cell. Its land in the run\'s '
              'water cells has no class.'))


def assigned(lat, lon, geo, cells, classes):
    """The design surface's land in the run's water cells, each such cell's land given the class of the nearest of the
    run's land cells (the great circle between the cells' centres on the turned Moon), and the classes' shares over all
    the design surface's land, its land in the run's land cells taking its cell's class; areas in km2."""
    area = bm.cell_areas_km2(lat, lon)
    dry = area - np.minimum((geo['sea_fraction'] + geo['rainfed_lake_fraction']) * area, area)
    rows, cols = np.asarray(cells['row']), np.asarray(cells['col'])
    land = np.zeros(area.shape, bool)
    land[rows, cols] = True
    target = np.argwhere(~land & (dry > 0.0))
    unit = lambda la, lo: np.stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)), axis=-1)  # noqa: E731
    la, lo = np.radians(np.asarray(lat, float)), np.radians(np.asarray(lon, float))
    source, sink = unit(la[rows], lo[cols]), unit(la[target[:, 0]], lo[target[:, 1]])
    nearest = np.argmax(sink @ source.T, axis=1)
    km = np.arccos(np.clip(np.sum(sink * source[nearest], axis=1), -1.0, 1.0)) * MOON_RADIUS / 1e3
    index = np.asarray(cells['moisture'])
    own, given = dry[rows, cols], dry[target[:, 0], target[:, 1]]
    keys = [c['key'] for c in classes]
    by_class = {k: round(float(given[index[nearest] == i].sum())) for i, k in enumerate(keys)}
    in_land = {k: float(own[index == i].sum()) for i, k in enumerate(keys)}
    total = float(own.sum() + given.sum())
    dry_keys = keys[:4]
    return dict(
        assigned=dict(cells=int(target.shape[0]), km2=round(float(given.sum())), by_class_km2=by_class,
                      nearest_land_km=dict(median=round(float(np.median(km))), p90=round(float(np.percentile(km, 90))),
                                           maximum=round(float(km.max()))),
                      rule='the class of the nearest of the run\'s land cells, centre to centre on the turned Moon; '
                           'assigned, not computed for the cell'),
        over_the_design_land=dict(
            land_km2=round(total),
            by_class={k: round((in_land[k] + by_class[k]) / total, 4) for k in keys},
            by_class_in_run_land_cells_only={k: round(in_land[k] / float(own.sum()), 4) for k in keys},
            drylands_share=round(sum(in_land[k] + by_class[k] for k in dry_keys) / total, 4),
            assigned_share=round(float(given.sum()) / total, 4),
            note=('Shares of the design surface\'s land: its land in each of the run\'s land cells takes the cell\'s '
                  'class, its land in the run\'s water cells the class of the nearest land cell (assigned). The shares '
                  'under shares.by_class are of the run\'s land cells, each counted whole.')))


def write_lakes(surface):
    """The design surface's lakes and channels in geography's drainage product's form, which the moisture map names as
    its lakes' source."""
    cells = np.flatnonzero(surface['body'] == LAKE).astype(np.int32)
    (LAKES / 'products').mkdir(parents=True, exist_ok=True)
    (LAKES / 'results').mkdir(parents=True, exist_ok=True)
    meta = dict(schema='terluna.geography.drainage/1', share=0.28, sea_level_m=surface['meta']['common_level_m'],
                grid=surface['meta']['grid'],
                evidence=('The atlas round\'s design surface (globe_surface.py): its lakes off the seas and its channels, '
                          'in the drainage product\'s form; the seas stand at their own levels, so "sea level" is the '
                          'common level only as the heights\' reference.'),
                reading_rule='As geography\'s drainage product.')
    np.savez_compressed(LAKES / 'products' / 'drainage_28pct_16ppd.npz', metadata=json.dumps(meta),
                        channel_cell=surface['channel_cell'], channel_discharge_m3s=surface['channel_q'],
                        lake_cell=cells, lake_level_m=surface['level'][cells], lake_spills=surface['spills'][cells])
    summary = json.loads(GS.OUT.read_text())['summary']
    (LAKES / 'results' / 'drainage.json').write_text(json.dumps(dict(
        schema='terluna.geography.drainage/1', source=str(GS.PRODUCT), source_sha256=sha256(GS.PRODUCT),
        lakes=dict(area_share=summary['lake_share'], bodies=summary['lakes'], volume_km3=summary['lake_volume_km3'],
                   closed_area_share=summary['closed_lake_share']),
        rivers=dict(discharge_to_sea_m3s=summary['rivers_into_the_seas_m3s'])), indent=1) + '\n')
    return LAKES


def build(log=print):
    started = time.time()
    surface = surface_grid()
    lakes = write_lakes(surface)
    frame_path, _ = frame_record()
    bm.RUNS = DL.RUNS
    bm.CACHE = str(DL.CACHE)
    bm.read_frame = settled_frame                    # the run's settled frame from its surface record
    seen = {}

    def settings(lat, lon, rotation=None, drainage=None):
        out = surface_settings(lat, lon, rotation, drainage)
        seen.update(lat=lat, lon=lon, geo=out)
        return out

    bm.settings = settings                          # the design surface's seas, lakes, coasts and tides
    bm.GEOGRAPHY['tides_grid'] = str(TF.PRODUCT)     # the tide of the design surface's own seas
    product = bm.build(RUN, YEARS, (), bm.STORE_YEARS, drainage=str(lakes), orientation=str(frame_path.relative_to(ROOT)))
    cells = product['cells']
    cover = coverage(seen['lat'], seen['lon'], seen['geo'], np.array(cells['row']), np.array(cells['col']))
    given = assigned(seen['lat'], seen['lon'], seen['geo'], cells, product['classes'])
    cover['assigned'] = given['assigned']
    product['land']['design_surface'] = cover
    product['shares']['over_the_design_land'] = given['over_the_design_land']
    coast_km = np.array(cells['coast_km'], float)
    untided = float(coast_km[np.array([v is None for v in cells['tide_m']])].sum())
    unread = [bm.GEOGRAPHY['atlas'], bm.GEOGRAPHY['atlas_grid']]     # the design surface's settings replace the atlas's
    for key in unread:
        product['producer']['inputs'].pop(key, None)
        product['producer']['input_schemas'].pop(key, None)
    product['producer'].update(
        model='biosphere/ecology/moisture_map.py through research/studies/atlas_round/moisture_final.py',
        wrapper=dict(files={f: sha256(ROOT / f) for f in CODE}),
        design_surface=dict(path=str(GS.PRODUCT), sha256=sha256(GS.PRODUCT), record=str(GS.OUT.relative_to(ROOT)),
                            record_sha256=sha256(GS.OUT)),
        tides=dict(path=str(TF.PRODUCT), sha256=sha256(TF.PRODUCT), record=str(TF.OUT.relative_to(ROOT)),
                   record_sha256=sha256(TF.OUT)))
    product['evidence'] = (
        'The moisture map of the design climate on the atlas round\'s final surface (the run at sunlight 0.930, years '
        '35-49, 300.3 K at the ground), classed by biosphere/ecology/moisture_map.py unchanged; the cells are the run\'s, '
        'in the settled turned frame of its surface record, and the settings (seas, lakes, closed lakes, coasts, tides) '
        'come from the round\'s design surface, every sea at its own level and the near-side sea parted. '
        + product['evidence'].replace('geology\'s orientation product', 'the final surface\'s settled frame'))
    product['window'] = dict(years=list(YEARS), note='the run\'s settled years; its earlier years cool toward them, so no '
                                                     'decade beside them is mapped')
    a = given['assigned']
    product['limits'] = [f'The classes cover the run\'s land cells, {cover["run_land_km2"] / 1e6:.1f} million km2: they hold '
                         f'{cover["design_land_in_run_land_km2"] / 1e6:.1f} million km2 of the design surface\'s land, '
                         f'{cover["design_lake_in_run_land_km2"] / 1e6:.1f} million of its lakes and '
                         f'{cover["design_sea_in_run_land_km2"] / 1e6:.1f} million of its seas. Another '
                         f'{cover["design_land_in_run_water_km2"] / 1e6:.1f} million km2 of its land, '
                         f'{100 * cover["design_land_unclassed_share"]:.0f}% of it, lies in cells the run\'s mask makes '
                         'water; it takes the class of the nearest land cell, '
                         f'{a["nearest_land_km"]["median"]} km away at the median and {a["nearest_land_km"]["maximum"]} km '
                         'at most, assigned and counted apart. In lake districts that class errs dry, since the run '
                         'holds no land cell whose air the district\'s own lakes moisten.',
                         'The polar light\'s thresholds are the elastic frame\'s, at a tilt of 1.613 degrees; the settled '
                         'frame\'s is 1.707.',
                         'The tide along the coasts is the monthly range of the design surface\'s seas from geography\'s '
                         'tide model with each sea\'s mean held (tides_final.py), an equilibrium tide checked against '
                         'a shallow-water solution on three seas; '
                         + (f'all {coast_km.sum():,.0f} km of coast on the run\'s land cells have a tide class.' if untided == 0
                            else f'{untided:,.0f} of the {coast_km.sum():,.0f} km of coast on the run\'s land cells lie in '
                                 'cells without a tide class.'),
                         'The gap\'s scale and the night-rain timing are the design run\'s of 1 November, as in the '
                         'moisture map on the corrected climate.']
    product['seconds'] = round(time.time() - started)
    log(f'moisture map on the final surface in {product["seconds"]} s')
    return product


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    product = build(log=lambda *a: print(*a, flush=True))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(product, indent=1) + '\n')
    s = product['shares']
    print(f"{product['land']['cells']} land cells, robust {s['robust']:.3f}")
    for key, row in s['by_class'].items():
        print(f"  {key}: {row['share']:.3f} as run, {row['robust_share']:.3f} robust")
    print('  mix', json.dumps({k: v for k, v in product['mix'].items() if k != 'earth_drylands'})[:1500])
    print('  soda', json.dumps({k: v for k, v in product['settings']['soda_lakes'].items() if k != 'note'}))
    print('  design surface', json.dumps({k: v for k, v in product['land']['design_surface'].items() if k != 'note'}))
    print('  over the design land', json.dumps({k: v for k, v in s['over_the_design_land'].items() if k != 'note'}))
    print(f'wrote {args.out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
