"""Boundary conditions and planet parameters for a three-dimensional climate model of the Open Moon.

    python -m climate.gcm.boundary      # writes climate/gcm/products/*.nc and moon_gcm_configuration.json

Nothing here runs a general circulation model. This prepares what any of them
needs, from the repository's own sources:

- planet parameters from shared/constants.json (radius, GM/R^2 gravity, the
  synchronous rotation period, the synodic solar day, the equator's 1.54 deg tilt
  to the ecliptic);
- land-sea fraction, land elevation above the water level and ocean depth,
  conservatively averaged from the geography domain's 1-degree product
  (terluna.geography.land-sea/1) to the grids GCMs commonly use: regular
  4 x 5 and 2 x 2.5 degree, and Gaussian T21 (32 x 64) and T42 (64 x 128);
- the same for the scenario's 28% seas with the rain-fed lakes of the geography
  domain's drainage estimate (terluna.geography.atlas-grid/1 and
  terluna.geography.drainage/1), adding each cell's water surface height:
  sea level for seas, the lake's own level for lakes above it;
- the atmosphere and shield scenarios the 1-D work uses, by reference to their
  data products, so a 3-D run starts from the same composition and sunlight.

The grids are area-conserving block averages of the 1-degree cells; a cell's
land fraction is its dry share, its elevation the mean ground height of its dry
part above the chosen water level, and its ocean depth the mean depth of its wet
part. These are storage geometries (see the geography product's evidence
statement), not shorelines shaped by rainfall, runoff or erosion.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import numpy as np
from scipy.io import netcdf_file

from shared import constants as c

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GEOGRAPHY = ROOT / 'geography' / 'products'
SCHEMA = 'terluna.climate.gcm-boundary/1'
LAKES = ('atlas_28pct_4ppd.npz', 'drainage_28pct_16ppd.npz')          # the scenario's seas, and its rain-fed lakes


def gaussian_latitudes(n):
    """Gaussian latitudes (degrees, north first) and cell edges in sin(latitude)."""
    x, w = np.polynomial.legendre.leggauss(n)
    order = np.argsort(-x)
    x, w = x[order], w[order]
    edges = 1.0 - np.concatenate([[0.0], np.cumsum(w)])
    return np.degrees(np.arcsin(x)), np.clip(edges, -1.0, 1.0)


def regular_latitudes(dlat):
    n = int(round(180.0 / dlat))
    edges_deg = 90.0 - dlat * np.arange(n + 1)
    return 0.5 * (edges_deg[1:] + edges_deg[:-1]), np.sin(np.radians(edges_deg))


GRIDS = {
    'regular_4x5': dict(lat=lambda: regular_latitudes(4.0), nlon=72),
    'regular_2x2.5': dict(lat=lambda: regular_latitudes(2.0), nlon=144),
    'gaussian_T21': dict(lat=lambda: gaussian_latitudes(32), nlon=64),
    'gaussian_T42': dict(lat=lambda: gaussian_latitudes(64), nlon=128),
}


def _overlap(src_edges, dst_edges):
    """Overlap matrix of 1-D intervals (dst x src); edges may run in either direction."""
    s0, s1 = np.minimum(src_edges[:-1], src_edges[1:]), np.maximum(src_edges[:-1], src_edges[1:])
    d0, d1 = np.minimum(dst_edges[:-1], dst_edges[1:]), np.maximum(dst_edges[:-1], dst_edges[1:])
    return np.clip(np.minimum(d1[:, None], s1[None, :]) - np.maximum(d0[:, None], s0[None, :]), 0.0, None)


def combine_lakes(height, sea, lake_cells, lake_levels, fine_cols, ratio, level):
    """Seas and lakes on a coarse grid (the atlas's): water fraction, water depth and water surface height
    (above the geoid) of each cell. `sea` marks whole sea cells; lakes are fine cells (`fine_cols` columns,
    `ratio` fine cells per coarse cell along each axis) on land, whose share of a coarse cell becomes its
    lake fraction and whose mean level its lake surface. The lake bed is taken at the coarse cell's height."""
    rows, cols = (lake_cells // fine_cols) // ratio, (lake_cells % fine_cols) // ratio
    count = np.zeros(height.shape)
    level_sum = np.zeros(height.shape)
    np.add.at(count, (rows, cols), 1.0)
    np.add.at(level_sum, (rows, cols), lake_levels)
    lake = np.where(sea, 0.0, count / ratio ** 2)
    lake_level = np.where(count > 0, level_sum / np.maximum(count, 1.0), level)
    surface = np.where(sea, level, lake_level)
    depth = np.where(sea, level - height, np.where(lake > 0, np.maximum(lake_level - height, 0.0), 0.0))
    return np.where(sea, 1.0, lake), depth, surface


def lakes_product(atlas_name, drainage_name):
    """The atlas's seas with the drainage estimate's rain-fed lakes, in the land-sea product's form."""
    atlas, drainage = load_product(atlas_name), load_product(drainage_name)
    level = atlas['meta']['sea_level_m']
    if abs(drainage['meta']['sea_level_m'] - level) > 0.5:
        raise ValueError('the drainage product was routed at another sea level than the atlas')
    height = atlas['height_m'].astype(float)
    fine = drainage['meta']['grid']
    wet, depth, surface = combine_lakes(height, atlas['water_label'] > 0, drainage['lake_cell'], drainage['lake_level_m'],
                                        fine['cols'], fine['pixels_per_degree'] * 180 // height.shape[0], level)
    return dict(lat_deg=atlas['lat_deg'], lon_deg=atlas['lon_deg'], water_fraction=wet, mean_height_m=height,
                mean_water_depth_m=depth, water_surface_m=surface, meta=dict(level_m=level),
                sha256=f"{atlas['sha256']}+{drainage['sha256']}")


def regrid(product, grid):
    """Area-conserving averages of a land-sea product onto a target grid."""
    lat_src = product['lat_deg']
    dlat = abs(lat_src[1] - lat_src[0])
    src_mu_edges = np.sin(np.radians(np.concatenate([[lat_src[0] + dlat / 2], lat_src - dlat / 2])))
    dlon = product['lon_deg'][1] - product['lon_deg'][0]
    src_lon_edges = np.concatenate([[0.0], product['lon_deg'] + dlon / 2])
    lat, dst_mu_edges = GRIDS[grid]['lat']()
    nlon = GRIDS[grid]['nlon']
    dst_lon_edges = np.linspace(0.0, 360.0, nlon + 1) - 180.0 / nlon          # cells centred on 0, dl, 2dl, ...
    lon = 0.5 * (dst_lon_edges[1:] + dst_lon_edges[:-1])
    a_lat = _overlap(src_mu_edges, dst_mu_edges)                              # sin-lat overlaps (area weights)
    # Longitude overlaps with wrap-around: shift the source cells by -360, 0 and +360 degrees.
    a_lon = sum(_overlap(src_lon_edges + s, dst_lon_edges) for s in (-360.0, 0.0, 360.0))
    wet = product['water_fraction'].astype(float)
    dry = 1.0 - wet
    level = product['meta']['level_m']
    height = product['mean_height_m'].astype(float)
    depth = product['mean_water_depth_m'].astype(float)
    # The wet part's surface: the water level, or each cell's own (lakes above it).
    surface = product['water_surface_m'].astype(float) if 'water_surface_m' in product else np.full_like(wet, level)
    # The dry part's mean height above the water level: cell mean = wet*(surface - depth) + dry*h_dry.
    h_dry = np.where(dry > 1e-6, (height - wet * (surface - depth)) / np.maximum(dry, 1e-6), height) - level

    def average(field, weight):
        num = a_lat @ (field * weight) @ a_lon.T
        den = a_lat @ weight @ a_lon.T
        return num, den

    area = a_lat @ np.ones_like(wet) @ a_lon.T
    water_num, _ = average(wet, np.ones_like(wet))
    land_frac = 1.0 - water_num / area
    elev_num, elev_den = average(np.maximum(h_dry, 0.0), dry)
    depth_num, depth_den = average(depth, wet)
    surf_num, surf_den = average(surface - level, wet)
    out = dict(lat=lat, lon=lon, land_fraction=np.clip(land_frac, 0.0, 1.0),
               land_elevation_m=np.where(elev_den > 0, elev_num / np.maximum(elev_den, 1e-30), 0.0),
               ocean_depth_m=np.where(depth_den > 0, depth_num / np.maximum(depth_den, 1e-30), 0.0),
               cell_area_fraction=area / area.sum())
    if 'water_surface_m' in product:
        out['water_surface_m'] = np.where(surf_den > 0, surf_num / np.maximum(surf_den, 1e-30), 0.0)
    return out


def load_product(name):
    path = GEOGRAPHY / name
    data = np.load(path)
    out = {k: data[k] for k in data.files if k != 'metadata'}
    out['meta'] = json.loads(str(data['metadata']))
    out['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()[:16]
    return out


def planet():
    return dict(
        radius_m=c.MOON_RADIUS, gravity_m_s2=c.MOON_SURFACE_GRAVITY, gm_m3_s2=c.MOON_GM,
        rotation_period_days=c.SIDEREAL_MONTH_DAYS, solar_day_days=c.SYNODIC_MONTH_DAYS,
        rotation_rate_rad_s=2 * math.pi / (c.SIDEREAL_MONTH_DAYS * 86400.0),
        obliquity_deg=c.MOON_EQUATOR_TO_ECLIPTIC_DEG,
        orbit='Earth-Moon barycentre about the Sun: 1 AU, Earth-orbit eccentricity; the Moon rotates synchronously '
              'with its month, so the solar day is the synodic month',
        solar_constant_w_m2=c.SOLAR_CONSTANT,
        notes='Rotation is synchronous; Earth stays fixed in the nearside sky. Lunar eclipses (the Moon in Earth\'s '
              'shadow, a few hours about twice a year) and earthshine (under 0.1 W/m2) are negligible forcings.')


def write_netcdf(path, fields, attrs):
    with netcdf_file(path, 'w') as f:
        f.createDimension('lat', fields['lat'].size)
        f.createDimension('lon', fields['lon'].size)
        for name, dims, units in (('lat', ('lat',), 'degrees_north'), ('lon', ('lon',), 'degrees_east')):
            v = f.createVariable(name, 'f8', dims); v[:] = fields[name]; v.units = units
        for name, units in (('land_fraction', '1'), ('land_elevation_m', 'm above the water level'),
                            ('ocean_depth_m', 'm'), ('cell_area_fraction', '1'),
                            ('water_surface_m', 'm above the water level (mean over the wet part)')):
            if name not in fields:
                continue
            v = f.createVariable(name, 'f4', ('lat', 'lon'))
            v[:] = fields[name].astype(np.float32)
            v.units = units
        for k, v in attrs.items():
            setattr(f, k, v)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'products')
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    products = {}
    for pct in (25, 35):
        name = f'land_sea_{pct}pct_1deg.npz'
        if not (GEOGRAPHY / name).is_file():
            print(f'BLOCKED: {name} missing; run python -m geography.water_inventory --masks 0.25 0.35', file=sys.stderr)
            return 1
        products[pct] = load_product(name)
    sources = {pct: f'geography/products/land_sea_{pct}pct_1deg.npz' for pct in products}
    evidence = {pct: 'Area-conserving averages of hydrostatic LOLA/GRAIL water filling; storage geometry, not '
                     'predicted shorelines.' for pct in products}
    if all((GEOGRAPHY / name).is_file() for name in LAKES):
        products['28_lakes'] = lakes_product(*LAKES)
        sources['28_lakes'] = ' and '.join(f'geography/products/{name}' for name in LAKES)
        evidence['28_lakes'] = ('Area-conserving averages of the atlas\'s seas at the scenario share and the drainage '
                                'estimate\'s rain-fed lakes (runoff from climate run A); lake beds at the atlas\'s '
                                '4 pixel/degree heights. A first estimate, not predicted shorelines.')
    else:
        print(f'no lakes product: {" or ".join(LAKES)} missing (python -m geography.atlas; python -m geography.drainage)',
              file=sys.stderr)
    files = {}
    for pct, prod in products.items():
        for grid in GRIDS:
            fields = regrid(prod, grid)
            path = args.out / (f'moon_{pct}pct_water_{grid}.nc' if isinstance(pct, int) else f'moon_28pct_water_lakes_{grid}.nc')
            attrs = dict(schema=SCHEMA, source_product=sources[pct], source_sha256=prod['sha256'],
                         water_level_m_above_geoid=prod['meta']['level_m'], evidence=evidence[pct])
            write_netcdf(path, fields, attrs)
            files[path.name] = dict(grid=grid, water_percent=pct if isinstance(pct, int) else 28, lakes=pct == '28_lakes',
                                    land_fraction_mean=float(
                np.sum(fields['land_fraction'] * fields['cell_area_fraction'])))
            print(path.name, f"land {files[path.name]['land_fraction_mean']:.3f}")
    code = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]
    config = dict(
        schema='terluna.climate.gcm-configuration/1',
        producer=dict(domain='climate', files={'gcm/boundary.py': code}),
        evidence=('Inputs for a future 3-D climate experiment; no GCM has been run. Planet values derive from '
                  'shared/constants.json; surfaces from the geography product; atmosphere and sunlight from the '
                  'atmosphere and protection products named below.'),
        planet=planet(),
        atmosphere=dict(
            dry_air='Earth O2 partial pressure (21,227 Pa), Earth Ar/N2, CO2 400 ppm; N2 the balance',
            surface_pressures_pa=dict(design=121590.0, lower=101325.0),
            water='prognostic; seas from the land-sea files below',
            ozone=('from atmosphere/middle_atmosphere/results/profiles/<case>.csv for the chosen shield; none behind '
                   'the titania stack'),
            reference_profiles='atmosphere/middle_atmosphere/results/profiles (temperature, humidity, ozone, heating)'),
        sunlight=dict(
            spectrum='TSIS-1 HSRS (atmosphere/radiative_convective inputs)',
            shields='protection/spectra/shield_transmission.json (titania_stack, edge_200nm ... edge_310nm)'),
        surface=dict(land_sea_files=files, ocean='slab mixed layer; 50 m default, 10 and 100 m as sensitivity',
                     land_albedo_note='no vegetation or soil map exists yet; a uniform land albedo is a placeholder'),
    )
    (args.out / 'moon_gcm_configuration.json').write_text(json.dumps(config, indent=1) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
