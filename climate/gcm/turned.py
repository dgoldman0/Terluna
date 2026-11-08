"""The GCM's surface on the turned Moon: the lakes surface boundary.py builds (the atlas's seas, the drainage estimate's
rain-fed lakes and the atlas's heights) carried into the frame of domain/geology's orientation product, then averaged
onto the GCM grids as boundary.py averages it (5 November 2026).

    python -m climate.gcm.turned <orientation product> --source <branch@commit:path>
        # writes products/moon_28pct_water_lakes_<drainage>_turned_<frame>_<grid>.nc and
        # ../results/gcm/turned_surface_<frame>.json, <frame> the product's sha256[:16]

The seas' and lakes' uneven weight turns the Moon about 8°, and the second set of full simulations runs on the turned
Moon at the elastic state of the seas and lakes (main's register, 4 November; research/studies/moon_orientation). The
frame is the orientation product's `elastic_frame`: the rows of its rotation_today_to_turned are the turned Moon's axes
(toward Earth, the third, north) in today's body frame, so a point's turned unit vector is R n(today) and today's is
R^T n(turned), as geology/orientation.py converts them.

The ground and its water turn with the crust and nothing else changes: the water level, the lakes' levels and every
height stay today's, on the dry Moon's geoid, until the atlas work folds in the water's own gravity and the crust's sag
(register, 4 November). Each cell of the atlas's 4 pixel/degree grid in the turned frame takes the mean of S x S
points spread evenly over it in latitude and longitude, each read from the lakes surface at the point's place today
(the atlas cell holding it): the water fraction and the height as plain means, the water's depth and surface as means
over the wet points. The spin axis's tilt to the ecliptic is the product's for the elastic state, which the runner
takes in place of today's with the surface (exoplasim_run's frame).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from climate.gcm import boundary as bd

SCHEMA = 'terluna.climate.gcm-turned-surface/1'
ORIENTATION_SCHEMA = 'terluna.geology.moon-orientation/1'
SAMPLES = 3                        # points along each axis of an atlas cell
STATE = 'elastic'                  # the orientation product's state the frame belongs to
BANDS = np.arange(-90.0, 90.1, 10.0)   # latitude bands for the record's land shares


def turned_name(drainage_sha, frame_sha, grid):
    """The turned lakes surface built from the drainage product and the orientation product with these content
    hashes (sha256, 16 hex digits)."""
    return f'moon_28pct_water_lakes_{drainage_sha}_turned_{frame_sha}_{grid}.nc'


def unit(lat_deg, lon_deg):
    lat, lon = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat) * np.ones_like(lon)], axis=-1)


def to_turned(lat_deg, lon_deg, rotation):
    """Today's latitude and longitude (deg) to the turned Moon's, longitude in (-180, 180]."""
    p = unit(np.asarray(lat_deg, float), np.asarray(lon_deg, float)) @ np.asarray(rotation, float).T
    return np.degrees(np.arcsin(np.clip(p[..., 2], -1, 1))), np.degrees(np.arctan2(p[..., 1], p[..., 0]))


def to_today(lat_deg, lon_deg, rotation):
    """The turned Moon's latitude and longitude (deg) back to today's, longitude in (-180, 180]."""
    p = unit(np.asarray(lat_deg, float), np.asarray(lon_deg, float)) @ np.asarray(rotation, float)
    return np.degrees(np.arcsin(np.clip(p[..., 2], -1, 1))), np.degrees(np.arctan2(p[..., 1], p[..., 0]))


def read_frame(path):
    """The orientation product's frame: its rotation, the elastic state's tilt to the ecliptic and its places, with the
    product's sha256."""
    path = Path(path)
    product = json.loads(path.read_text())
    if product.get('schema') != ORIENTATION_SCHEMA:
        raise ValueError(f'{path} is not a {ORIENTATION_SCHEMA} product')
    rotation = np.asarray(product['elastic_frame']['rotation_today_to_turned'], float)
    if rotation.shape != (3, 3) or not np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-12) \
            or not np.isclose(np.linalg.det(rotation), 1.0, atol=1e-12):
        raise ValueError(f'{path}: rotation_today_to_turned is not a proper rotation')
    return dict(rotation=rotation, tilt_deg=float(product['report_states'][STATE]['tilt_to_ecliptic_deg']),
                frame=product['elastic_frame'], sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def _cells(lat_deg, lon_deg, grid_lat, grid_lon):
    """The rows and columns of the regular grid (cell centres grid_lat north first, grid_lon from 0 east) that hold
    each point."""
    dlat, dlon = abs(grid_lat[1] - grid_lat[0]), grid_lon[1] - grid_lon[0]
    north, west = grid_lat[0] + dlat / 2, grid_lon[0] - dlon / 2
    rows = np.clip(np.floor((north - lat_deg) / dlat).astype(int), 0, grid_lat.size - 1)
    cols = np.floor(np.mod(lon_deg - west, 360.0) / dlon).astype(int) % grid_lon.size
    return rows, cols


def turn_fields(lat, lon, fields, wet, rotation, samples=SAMPLES):
    """Each field of the regular grid (lat, lon cell centres) carried into the turned frame on the same grid. fields
    maps a name to (array, weighted): a plain mean over each turned cell's points, or, weighted, a mean over its wet
    points (nan where none is wet). wet is the water fraction, carried as a plain mean under the name 'wet'."""
    dlat, dlon = abs(lat[1] - lat[0]), lon[1] - lon[0]
    offsets = (np.arange(samples) + 0.5) / samples - 0.5
    sums = {name: np.zeros((lat.size, lon.size)) for name in fields}
    wet_sum = np.zeros((lat.size, lon.size))
    for a in offsets:
        for b in offsets:
            lat_t = (lat + a * dlat)[:, None] * np.ones(lon.size)
            lon_t = np.ones(lat.size)[:, None] * (lon + b * dlon)
            lat_0, lon_0 = to_today(lat_t, lon_t, rotation)
            rows, cols = _cells(lat_0, lon_0, lat, lon)
            w = wet[rows, cols]
            wet_sum += w
            for name, (field, weighted) in fields.items():
                sums[name] += field[rows, cols] * w if weighted else field[rows, cols]
    out = {'wet': wet_sum / samples ** 2}
    for name, (_, weighted) in fields.items():
        out[name] = np.where(wet_sum > 0, sums[name] / np.maximum(wet_sum, 1e-300), np.nan) if weighted \
            else sums[name] / samples ** 2
    return out


def turned_product(prod, sea, rotation, samples=SAMPLES):
    """The lakes surface (boundary.lakes_product) in the turned frame, in the same form, with its sea and lake
    shares; sea marks the atlas's whole sea cells."""
    level = prod['meta']['level_m']
    fields = dict(height=(prod['mean_height_m'].astype(float), False),
                  depth=(prod['mean_water_depth_m'].astype(float), True),
                  surface=(prod['water_surface_m'].astype(float), True),
                  sea=(sea.astype(float), False))
    t = turn_fields(prod['lat_deg'], prod['lon_deg'], fields, prod['water_fraction'].astype(float), rotation, samples)
    weights = np.cos(np.radians(prod['lat_deg']))[:, None] * np.ones(prod['lon_deg'].size)
    wet = t['wet']
    sea_share = float((t['sea'] * weights).sum() / weights.sum())
    water_share = float((wet * weights).sum() / weights.sum())
    return dict(lat_deg=prod['lat_deg'], lon_deg=prod['lon_deg'], water_fraction=wet, mean_height_m=t['height'],
                mean_water_depth_m=np.where(wet > 0, np.nan_to_num(t['depth']), 0.0),
                water_surface_m=np.where(wet > 0, np.nan_to_num(t['surface'], nan=level), level),
                meta=dict(level_m=level), sha256=prod['sha256'], drainage_sha256=prod['drainage_sha256'],
                routing_run=prod['routing_run'], sea_share=sea_share, lake_share=water_share - sea_share)


def band_shares(fields):
    """The land's share of each latitude band of a regridded surface (BANDS), area-weighted."""
    land, area, lat = fields['land_fraction'], fields['cell_area_fraction'], fields['lat']
    out = []
    for lo, hi in zip(BANDS[:-1], BANDS[1:]):
        rows = (lat >= lo) & (lat < hi)
        out.append(dict(band=f'{lo:+.0f} to {hi:+.0f}',
                        land=round(float((land[rows] * area[rows]).sum() / area[rows].sum()), 4) if rows.any() else None))
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('orientation', type=Path, help="domain/geology's orientation product (moon_orientation.json)")
    parser.add_argument('--source', required=True,
                        help='where the product was read, as branch@commit:path, for the record')
    parser.add_argument('--out', type=Path, default=bd.HERE / 'products')
    parser.add_argument('--samples', type=int, default=SAMPLES)
    args = parser.parse_args(argv)
    if not all((bd.GEOGRAPHY / name).is_file() for name in bd.LAKES):
        print(f'BLOCKED: {" or ".join(bd.LAKES)} missing (python -m geography.atlas; python -m geography.drainage)',
              file=sys.stderr)
        return 1
    frame = read_frame(args.orientation)
    frame_sha = frame['sha256'][:16]
    prod = bd.lakes_product(*bd.LAKES)
    sea = bd.load_product(bd.LAKES[0])['water_label'] > 0
    turned = turned_product(prod, sea, frame['rotation'], args.samples)
    args.out.mkdir(parents=True, exist_ok=True)
    code = {f'gcm/{name}': hashlib.sha256((bd.HERE / name).read_bytes()).hexdigest()[:16]
            for name in ('turned.py', 'boundary.py')}
    written, grids = [], {}
    for grid in bd.GRIDS:
        fields = bd.regrid(turned, grid)
        path = args.out / turned_name(prod['drainage_sha256'], frame_sha, grid)
        bd.write_netcdf(path, fields, dict(
            schema=bd.SCHEMA, source_product=' and '.join(f'geography/products/{name}' for name in bd.LAKES),
            source_sha256=prod['sha256'], orientation_product=args.source, orientation_sha256=frame['sha256'],
            water_level_m_above_geoid=prod['meta']['level_m'],
            producer=' '.join(f'climate/{k} {v}' for k, v in code.items()),
            evidence=('The lakes surface of boundary.py (the atlas\'s seas, the drainage estimate\'s rain-fed lakes of '
                      f'climate run {prod["routing_run"]}\'s runoff, the atlas\'s heights) carried into the turned '
                      'Moon\'s elastic frame by the orientation product\'s rotation, the dry Moon\'s geoid kept; a '
                      'storage geometry, not predicted shorelines.'),
            reading_rule=('as boundary.py\'s lakes surface: land_fraction is each cell\'s dry share, land_elevation_m '
                          'the mean height of its dry part and water_surface_m of its wet part above the water level, '
                          'ocean_depth_m the mean depth of the wet part, cell_area_fraction its share of the sphere; '
                          'latitude and longitude are the turned Moon\'s')))
        written.append(path)
        if grid == 'gaussian_T21':
            grids['today'] = bd.regrid(prod, grid)
            grids['turned'] = fields
    today_t21 = bd.read_t21(args.out / bd.lakes_name(prod['drainage_sha256'], 'gaussian_T21')) \
        if (args.out / bd.lakes_name(prod['drainage_sha256'], 'gaussian_T21')).is_file() else None
    new_t21 = bd.read_t21(args.out / turned_name(prod['drainage_sha256'], frame_sha, 'gaussian_T21'))
    places = []
    for name, p in frame['frame'].get('places', {}).items():
        lat, lon = to_turned(p['today'][0], p['today'][1], frame['rotation'])
        places.append(dict(place=name, today=p['today'], turned=[round(float(lat), 3), round(float(lon), 3)],
                           product_turned=p['turned']))
    record = dict(
        schema=SCHEMA, producer=dict(domain='climate', files=code),
        inputs={**{f'geography/products/{name}': hashlib.sha256((bd.GEOGRAPHY / name).read_bytes()).hexdigest()
                   for name in bd.LAKES}, args.source: frame['sha256']},
        frame=dict(state=STATE, rotation_today_to_turned=frame['rotation'].tolist(),
                   north_pole=frame['frame']['north_pole'], sub_earth_point=frame['frame']['sub_earth_point'],
                   pole_move_deg=frame['frame']['pole_move_deg'], tilt_to_ecliptic_deg=frame['tilt_deg']),
        samples_per_axis=args.samples, routing_run=prod['routing_run'],
        surface=turned_name(prod['drainage_sha256'], frame_sha, 'gaussian_T21'),
        files={path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in written},
        shares=dict(today=dict(seas=prod['sea_share'], lakes=prod['lake_share'],
                               water=prod['sea_share'] + prod['lake_share']),
                    turned=dict(seas=turned['sea_share'], lakes=turned['lake_share'],
                                water=turned['sea_share'] + turned['lake_share'])),
        t21=dict(turned=dict(water_share=new_t21['water_share'], mask_land_share=new_t21['mask_land_share'],
                             mask_water_cells=int((new_t21['mask'] == 0).sum()), cells=int(new_t21['mask'].size)),
                 today=None if today_t21 is None else dict(
                     water_share=today_t21['water_share'], mask_land_share=today_t21['mask_land_share'],
                     mask_water_cells=int((today_t21['mask'] == 0).sum()))),
        land_by_band=dict(today=band_shares(grids['today']), turned=band_shares(grids['turned'])),
        places=places,
        evidence=('A calculation: the lakes surface of boundary.py carried into the turned Moon\'s elastic frame by '
                  'domain/geology\'s rotation and averaged onto the GCM grids; the water level, the lakes\' levels and '
                  'the heights are today\'s, on the dry Moon\'s geoid.'),
        reading_rule=('shares are of the Moon\'s surface on the atlas\'s 4 pixel/degree grid, today and turned; t21 gives '
                      'the T21 files\' water share and the binary mask PlaSim runs; land_by_band the land\'s share of '
                      'each 10° band on the T21 grid, today\'s latitudes and the turned Moon\'s; places compares this '
                      'module\'s rotation with the orientation product\'s own turned coordinates; the runner\'s frame '
                      'reads surface and frame.tilt_to_ecliptic_deg.'))
    folder = bd.RESULTS if Path(args.out).resolve() == (bd.HERE / 'products').resolve() else Path(args.out)
    path = folder / f'turned_surface_{frame_sha}.json'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(record, indent=1) + '\n')
    print(f'{path}: water {record["shares"]["turned"]["water"]:.5f} (today {record["shares"]["today"]["water"]:.5f}); '
          f'T21 mask land {new_t21["mask_land_share"]:.4f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
