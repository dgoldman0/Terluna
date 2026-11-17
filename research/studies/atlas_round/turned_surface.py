"""The atlas round's surface for the design climate's run on the final surface: the redrawn floor's seas and the governed
lakes, turned into the Moon's settled frame under the final water and averaged onto the GCM grids as
climate/gcm/turned.py averages the atlas's.

    ulimit -v 9000000; OPENBLAS_NUM_THREADS=1 nice -n 10 python -m research.studies.atlas_round.turned_surface   # under a minute

The register (16 November) runs the round's design climate once on the final surface, with the redrawn shores, the
seas at their chosen levels and the new lakes, at the settled turn under the final water, 6.97°. This builds that
surface for the GCM runner:

- the water: the governed routing on the routing's own budget, the design climate's rain (governed_drainage.py), with
  each of the six seas at its level over its cells (Moscoviense, Orientale and Smythii–Marginis full at their rims;
  the near-side ocean, the salt sea and Humboldtianum at the common level) and each lake at its own level, on the
  redrawn floor's 16 px/deg heights (LOLA with the settled case's U - N, redrawn_drainage.load_surface);
- the lakes surface on the atlas's 4 px/deg grid in boundary.lakes_product's form: a cell's water fraction is its share
  of wet 16 px/deg cells, its height the mean of all sixteen, and its water's depth and surface the means over its wet
  cells, so the dry part's height is exactly the mean of its dry cells. Heights and levels stand above the settled
  state's geoid, and the files count them from the redrawn common level;
- the frame: settled_turn.py's frame_today_to_turned and the spin axis's tilt to the ecliptic, the central profile 100
  million years after the water arrives, under the final water;
- the turned surface: turned.turned_product's sampling into the frame, then boundary's GCM grids, written in
  climate/gcm/products as moon_28pct_water_lakes_<governed>_turned_<frame>_<grid>.nc (turned.turned_name), named by
  the governed routing's content and the settled-turn product's, which the runner's experiment A28S reads; and its
  record in turned.py's form (turned.SCHEMA), results/turned_surface_<frame>.json, which the runner's --frame takes.

The dry part's height is averaged as it stands (regrid_signed): the redrawn floor's dry hollows below the common level,
1.75 million km² on the design routing, keep their depth in their cells' mean ground, where every atlas surface held
its hollows below the water level as water and boundary.regrid averages the dry ground from the water level up. The
runner sets PlaSim's ground no lower than the water level (exoplasim_run.surface_files), and the record counts the T21
land cells that this raises and by how much. turned.py and boundary.py are used as they stand, since the products of
other lanes record their content.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from climate.gcm import boundary as bd
from climate.gcm import turned
from research.studies.atlas_round import governed_drainage as GD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.atlas_round import settled_turn as ST

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PRODUCTS = bd.HERE / 'products'
STATE = 'settled'                    # the frame's state: the central profile 100 million years after the water arrives
FACTOR = 4                           # 16 px/deg cells along each axis of a 4 px/deg cell
BELOW_M = 20000.0                    # a level this far below the water level lies under all of the Moon's ground
CODE = ('research/studies/atlas_round/turned_surface.py', 'climate/gcm/turned.py', 'climate/gcm/boundary.py',
        'research/studies/atlas_round/redrawn_drainage.py')


def record_path(frame_sha):
    return HERE / 'results' / f'turned_surface_{frame_sha[:16]}.json'


def coarse(h, sea, lake_cell, lake_level, levels, factor=FACTOR):
    """The lakes surface on the coarse grid from the fine one: h the fine heights (rows north first, columns from 0 E),
    sea each fine cell's index into levels (-1 off the seas), lake_cell and lake_level the lakes' fine cells and
    levels. Returns the coarse water fraction, mean height, mean depth and surface over the wet part (0 and nan where
    dry), the seas' fraction, and counts of the fine cells that need a look: lake cells on a sea, and wet cells whose
    water stands no higher than their ground."""
    rows, cols = h.shape
    level = np.full(h.size, np.nan)
    on_sea = sea.ravel() >= 0
    level[on_sea] = np.asarray(levels, float)[sea.ravel()[on_sea]]
    lakes_on_seas = int(on_sea[lake_cell].sum())
    level[lake_cell] = lake_level
    wet = ~np.isnan(level)
    flat = h.ravel().astype(float)
    dry_water = int((wet & ~(level > flat)).sum())
    depth = np.where(wet, np.maximum(np.nan_to_num(level) - flat, 0.0), 0.0)
    block = lambda a: a.reshape(rows // factor, factor, cols // factor, factor).sum(axis=(1, 3))
    n = float(factor * factor)
    count = block(wet.reshape(rows, cols).astype(float))
    some = count > 0
    return dict(water_fraction=count / n, mean_height_m=block(flat.reshape(rows, cols)) / n,
                mean_water_depth_m=np.where(some, block(depth.reshape(rows, cols)) / np.maximum(count, 1.0), 0.0),
                water_surface_m=np.where(some, block(np.nan_to_num(level).reshape(rows, cols)) / np.maximum(count, 1.0),
                                         np.nan),
                sea_fraction=block(on_sea.reshape(rows, cols).astype(float)) / n,
                lakes_on_seas=lakes_on_seas, water_at_or_below_ground=dry_water)


def frame_from_rotation(rotation, tilt_deg, state, sha256):
    """A frame in turned.read_frame's form from a rotation alone (rows: the turned Moon's axes, toward Earth, the third,
    north, in today's body frame), the spin axis's tilt to the ecliptic in it, the state it belongs to and its
    source's sha256: the north pole and the point under Earth are the rotation's third and first rows in today's
    coordinates, and the pole's move the angle between today's pole and the turned one."""
    r = np.asarray(rotation, float)
    if r.shape != (3, 3) or not np.allclose(r @ r.T, np.eye(3), atol=1e-12) \
            or not np.isclose(np.linalg.det(r), 1.0, atol=1e-12):
        raise ValueError(f'{state}: rotation_today_to_turned is not a proper rotation')
    place = lambda v: dict(lat_deg=float(np.degrees(np.arcsin(np.clip(v[2], -1.0, 1.0)))),
                           lon_deg=float(np.degrees(np.arctan2(v[1], v[0]))))
    return dict(rotation=r, tilt_deg=float(tilt_deg), state=state, sha256=sha256,
                frame=dict(north_pole=place(r[2]), sub_earth_point=place(r[0]),
                           pole_move_deg=round(float(np.degrees(np.arccos(np.clip(r[2, 2], -1.0, 1.0)))), 3)))


def regrid_signed(product, grid, below=BELOW_M):
    """boundary.regrid with the dry ground's height averaged as it stands. boundary.regrid averages the dry part's
    height above the water level from zero up, which clears the rounding of the atlas's surfaces, where no dry ground
    lies below the water level. Counted from a level `below` metres lower, no ground lies under that level, so nothing
    is clipped; the dry ground's height and the water's surface are then counted from the water level again in the
    cells that have a dry or a wet part to average (boundary.regrid gives 0 where a cell has none)."""
    level = product['meta']['level_m']
    wet = np.asarray(product['water_fraction'], float)
    surface = np.asarray(product['water_surface_m'], float) if 'water_surface_m' in product else np.full_like(wet, level)
    bed = surface - np.asarray(product['mean_water_depth_m'], float)
    dry_height = np.where(1.0 - wet > 1e-6, (np.asarray(product['mean_height_m'], float) - wet * bed)
                          / np.maximum(1.0 - wet, 1e-6), np.asarray(product['mean_height_m'], float)) - level
    if np.any(dry_height[wet < 1.0] <= -below) or np.any(surface[wet > 0.0] - level <= -below):
        raise ValueError(f'ground or water lies more than {below:.0f} m below the water level')
    shifted = bd.regrid(dict(product, meta=dict(product['meta'], level_m=level - below)), grid)
    ground = shifted['land_elevation_m']
    out = dict(shifted, land_elevation_m=np.where(ground > 0.0, ground - below, 0.0))
    if 'water_surface_m' in shifted:
        out['water_surface_m'] = np.where(shifted['water_surface_m'] > 0.0, shifted['water_surface_m'] - below, 0.0)
    return out


def t21_summary(fields):
    """A T21 surface's water share and PlaSim's binary mask as the runner makes them, from regridded fields held in
    memory at the precision the files store (as boundary.read_t21 reads a file)."""
    from climate.gcm.exoplasim_run import land_mask
    land = np.asarray(fields['land_fraction'], np.float32).astype(float)
    area = np.asarray(fields['cell_area_fraction'], np.float32).astype(float)
    water = 1.0 - float((land * area).sum() / area.sum())
    mask = land_mask(land, area, water)
    return dict(water_share=water, mask=mask, mask_land_share=float((mask * area).sum() / area.sum()))


def write_surface(prod, sea, frame, *, source, out, inputs, producer, describe, path, samples=turned.SAMPLES):
    """The turned surface of a lakes product (boundary.lakes_product's form; sea its sea cells, whole or in part) in a
    frame (frame_from_rotation's form), written onto boundary's GCM grids in `out` under turned.turned_name with the
    dry ground's height as it stands, and its record at `path` in the form turned.py writes and the runner's --frame
    reads. describe gives the files' producer, source_product and evidence (file_evidence) and the record's evidence.
    Returns the record."""
    frame_sha = frame['sha256'][:16]
    moved = turned.turned_product(prod, sea, frame['rotation'], samples)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    written, grids = [], {}
    for grid in bd.GRIDS:
        fields = regrid_signed(moved, grid)
        file = out / turned.turned_name(prod['drainage_sha256'], frame_sha, grid)
        bd.write_netcdf(file, fields, dict(
            schema=bd.SCHEMA, source_product=describe['source_product'],
            source_sha256=prod['sha256'], orientation_product=source, orientation_sha256=frame['sha256'],
            water_level_m_above_geoid=prod['meta']['level_m'], producer=describe['producer'],
            evidence=describe['file_evidence'],
            reading_rule=('as boundary.py\'s lakes surface: land_fraction is each cell\'s dry share, land_elevation_m '
                          'the mean height of its dry part and water_surface_m of its wet part above the water level, '
                          'ocean_depth_m the mean depth of the wet part, cell_area_fraction its share of the sphere; '
                          'latitude and longitude are the turned Moon\'s; the dry part\'s height is averaged as it '
                          'stands, below the water level where the ground is')))
        written.append(file)
        if grid == 'gaussian_T21':
            grids['today'] = regrid_signed(prod, grid)
            grids['turned'] = fields
    today_t21 = t21_summary(grids['today'])
    new_t21 = bd.read_t21(out / turned.turned_name(prod['drainage_sha256'], frame_sha, 'gaussian_T21'))
    record = dict(
        schema=turned.SCHEMA, producer=producer, inputs=inputs,
        frame=dict(state=frame['state'], rotation_today_to_turned=frame['rotation'].tolist(),
                   north_pole=frame['frame']['north_pole'], sub_earth_point=frame['frame']['sub_earth_point'],
                   pole_move_deg=frame['frame']['pole_move_deg'], tilt_to_ecliptic_deg=frame['tilt_deg']),
        samples_per_axis=samples, routing_run=prod['routing_run'],
        surface=turned.turned_name(prod['drainage_sha256'], frame_sha, 'gaussian_T21'),
        files={file.name: hashlib.sha256(file.read_bytes()).hexdigest() for file in written},
        shares=dict(today=dict(seas=prod['sea_share'], lakes=prod['lake_share'],
                               water=prod['sea_share'] + prod['lake_share']),
                    turned=dict(seas=moved['sea_share'], lakes=moved['lake_share'],
                                water=moved['sea_share'] + moved['lake_share'])),
        t21=dict(turned=dict(water_share=new_t21['water_share'], mask_land_share=new_t21['mask_land_share'],
                             mask_water_cells=int((new_t21['mask'] == 0).sum()), cells=int(new_t21['mask'].size)),
                 today=dict(water_share=today_t21['water_share'], mask_land_share=today_t21['mask_land_share'],
                            mask_water_cells=int((today_t21['mask'] == 0).sum()))),
        land_by_band=dict(today=turned.band_shares(grids['today']), turned=turned.band_shares(grids['turned'])),
        places=[],
        evidence=describe['evidence'],
        reading_rule=('as turned.py\'s records: shares are of the Moon\'s surface on the atlas\'s 4 pixel/degree grid, '
                      'today and turned; t21 gives the T21 files\' water share and the binary mask PlaSim runs; '
                      'land_by_band the land\'s share of each 10° band on the T21 grid, today\'s latitudes and the '
                      'turned Moon\'s; the runner\'s frame reads surface and frame.tilt_to_ecliptic_deg.'))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=1) + '\n')
    return record


def lakes_surface(log=print):
    """The governed routing's water on the redrawn floor in boundary.lakes_product's form on the 4 px/deg grid, with the
    seas' fraction and the record's details."""
    surf = RD.load_surface(log)
    with np.load(GD.GOVERNED_ROUTING, allow_pickle=False) as d:
        meta = json.loads(str(d['metadata']))
        sea, lake_cell, lake_level = d['sea'], d['lake_cell'], d['lake_level_m'].astype(float)
    governed_sha = RD.sha256(GD.GOVERNED_ROUTING)
    h = surf['h']
    level0 = float(surf['level0'])
    if abs(level0 - meta['common_level_m']) > RD.RIM_CHECK_M:
        raise SystemExit(f'the governed routing\'s common level {meta["common_level_m"]:.3f} m is not the floor\'s '
                         f'{level0:.3f} m')
    c = coarse(h, sea.reshape(h.shape), lake_cell, lake_level, [meta['levels_m'][k] for k in meta['seas']])
    rows, cols = c['water_fraction'].shape
    lat = 90.0 - (np.arange(rows) + 0.5) / (RD.dr.PPD / FACTOR)
    lon = (np.arange(cols) + 0.5) / (RD.dr.PPD / FACTOR)
    weights = np.cos(np.radians(lat))[:, None] * np.ones(cols)
    share = lambda a: float((a * weights).sum() / weights.sum())
    wet = c['water_fraction']
    climatology = json.loads(RD.OUT.read_text())['climatology']
    prod = dict(lat_deg=lat, lon_deg=lon, water_fraction=wet, mean_height_m=c['mean_height_m'],
                mean_water_depth_m=c['mean_water_depth_m'],
                water_surface_m=np.where(wet > 0, c['water_surface_m'], level0), meta=dict(level_m=level0),
                sha256=f'{RD.sha256(surf["grids"])[:16]}+{governed_sha[:16]}', drainage_sha256=governed_sha[:16],
                routing_run=climatology['run'], sea_share=share(c['sea_fraction']),
                lake_share=share(wet - c['sea_fraction']))
    details = dict(seas=meta['seas'], levels_above_common_m={k: round(v - level0, 1) for k, v in meta['levels_m'].items()},
                   common_level_m=level0, lakes_on_seas=c['lakes_on_seas'],
                   water_at_or_below_ground=c['water_at_or_below_ground'], climatology=climatology,
                   governed_sha256=governed_sha, grids=str(surf['grids']), grids_sha256=RD.sha256(surf['grids']))
    return prod, c['sea_fraction'], details


def settled_frame():
    """The settled turn under the final water (settled_turn.py) as a frame in turned.read_frame's form."""
    st = json.loads(ST.OUT.read_text())
    return frame_from_rotation(st['frame_today_to_turned'], st['cases']['final_water']['settled']['tilt_to_ecliptic_deg'],
                               STATE, RD.sha256(ST.OUT))


def lifted(path):
    """The T21 land cells whose mean ground lies below the water level, which the runner raises to it: their count,
    their share of the Moon and the area-weighted mean and the deepest of their mean ground."""
    from scipy.io import netcdf_file
    t21 = bd.read_t21(path)
    with netcdf_file(path, 'r', mmap=False) as d:
        ground = np.asarray(d.variables['land_elevation_m'][:], float)
    below = (t21['mask'] > 0) & (ground < 0)
    area = t21['area']
    return dict(cells=int(below.sum()), share_of_the_moon=round(float(area[below].sum() / area.sum()), 5),
                mean_ground_m=round(float((ground[below] * area[below]).sum() / max(area[below].sum(), 1e-300)), 1)
                if below.any() else None,
                deepest_ground_m=round(float(ground[below].min()), 1) if below.any() else None)


def build(log=print):
    started = time.time()
    prod, sea, details = lakes_surface(log)
    frame = settled_frame()
    code = {f: RD.sha256(ROOT / f)[:16] for f in CODE}
    inputs = {str(p.relative_to(ROOT)): RD.sha256(p) for p in (ST.OUT, GD.OUT, RD.OUT)}
    source = str(ST.OUT.relative_to(ROOT))
    path = record_path(frame['sha256'])
    record = write_surface(
        prod, sea, frame, source=source, out=PRODUCTS, inputs=inputs, path=path,
        producer=dict(domain='research', study='atlas_round', files=code),
        describe=dict(
            producer=' '.join(f'{k} {v}' for k, v in code.items()),
            source_product=('the atlas round\'s redrawn floor (LOLA with the settled case\'s U - N) and its governed '
                            'routing (research/studies/atlas_round/governed_drainage.py, governed_routing.npz)'),
            file_evidence=('The atlas round\'s final surface: the redrawn floor at the crust\'s settled state with the '
                           'seas at their chosen levels, the upper seas full at their rims, and the governed lakes of '
                           f'climate run {prod["routing_run"]}\'s rain, carried into the Moon\'s settled frame under '
                           'the final water (research/studies/atlas_round/settled_turn.py); heights above the settled '
                           'state\'s geoid. A storage geometry on the routing, not shorelines shaped by erosion.'),
            evidence=('A calculation: the governed routing\'s seas and lakes on the redrawn floor, averaged onto the '
                      'atlas\'s 4 pixel/degree grid from 16 pixel/degree cells, carried into the settled frame by the '
                      'atlas round\'s rotation and averaged onto the GCM grids with the dry ground\'s height as it '
                      'stands. Heights and levels stand above the settled state\'s geoid and are counted from the '
                      'redrawn common level.')))
    record.update(
        ignored_inputs={str(GD.GOVERNED_ROUTING): details['governed_sha256'], details['grids']: details['grids_sha256']},
        water=dict(common_level_m_above_geoid=round(details['common_level_m'], 3),
                   seas_above_common_level_m=details['levels_above_common_m'], seas=details['seas'],
                   lake_cells_on_seas=details['lakes_on_seas'],
                   wet_cells_at_or_below_their_ground=details['water_at_or_below_ground']),
        ground_below_the_water_level=lifted(PRODUCTS / record['surface']),
        frame_source=dict(product=source, state=STATE,
                          meaning=('the central profile 100 million years after the water arrives, under the final '
                                   'water: the redrawn seas with the upper seas full and the governed lakes')),
        seconds=round(time.time() - started))
    path.write_text(json.dumps(record, indent=1) + '\n')
    log(f'{path.relative_to(ROOT)}: water {record["shares"]["turned"]["water"]:.5f} turned, '
        f'{record["shares"]["today"]["water"]:.5f} today (seas {record["shares"]["today"]["seas"]:.5f}); '
        f'T21 water {record["t21"]["turned"]["water_share"]:.4f}, mask land {record["t21"]["turned"]["mask_land_share"]:.4f}; '
        f'land cells raised to the water level {record["ground_below_the_water_level"]}; {record["seconds"]} s')
    return path, record


def main(argv=None) -> int:
    argparse.ArgumentParser(description=__doc__.splitlines()[0]).parse_args(argv)
    build(log=lambda *a: print(*a, flush=True))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
