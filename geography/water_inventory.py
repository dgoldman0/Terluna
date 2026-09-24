"""Water inventories against lunar topography, written as a data product.

    python -m geography.water_inventory            # results/ (tables and JSON)

For a range of common water levels above the geoid: the flooded area, the water
needed (as a global equivalent layer and as mass) and the separate seas. For the
closed depressions: the major joins as water rises, labelled by the named basin
centres they contain. Hydrostatics only; where water actually stays is a climate
and hydrology question.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

from geography import topography as tp, basins as bs, fetch_inputs

HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.geography.water-inventory/1'
EVIDENCE = ('Hydrostatic filling of LOLA topography above the GRAIL geoid (degree 200, plus static Earth tide and '
            'rotation). Level curves assume one common water level; basin capacities are the volumes closed '
            'depressions hold when they first spill. No rainfall, runoff, evaporation, groundwater, ice, erosion or '
            'crustal loading: these are storage geometries for choosing an inventory, not predicted shorelines.')
TARGET_FRACTIONS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.71)


def producer():
    files = ['topography.py', 'basins.py', 'water_inventory.py']
    return dict(domain='geography', files={f: hashlib.sha256((HERE / f).read_bytes()).hexdigest()[:16] for f in files},
                inputs_manifest=hashlib.sha256((HERE / 'inputs.json').read_bytes()).hexdigest()[:16])


def compute(fine_ppd=16, tree_ppd=4, geoid_degree=200):
    height, geoid, grid = tp.height_above_geoid(fine_ppd, geoid_degree)
    area = grid.cell_area_m2
    levels = np.linspace(height.min(), height.max(), 40001)
    frac, vol = bs.level_curve(height, area, levels)
    rows = []
    for target in TARGET_FRACTIONS:
        k = int(np.searchsorted(frac, target))
        sizes, _ = bs.seas(height, area, levels[k])
        rows.append(dict(area_fraction=round(float(frac[k]), 4), level_m=round(float(levels[k]), 1),
                         global_layer_m=round(float(vol[k] / bs.AREA), 1),
                         mass_kg=float(vol[k] * bs.WATER_DENSITY),
                         mean_depth_m=round(float(vol[k] / (frac[k] * bs.AREA)), 0),
                         seas_over_0_1pct=len(sizes), largest_seas=[round(float(x), 4) for x in sizes[:5]]))
    h4, _, g4 = tp.height_above_geoid(tree_ppd, geoid_degree)
    records = bs.depression_tree(h4, g4.cell_area_m2, g4.lat_deg, g4.lon_deg)
    merges = []
    for level, parts in bs.major_merges(records, 2e-3):
        if level > 2000:
            break
        merges.append(dict(level_m=round(level, 0), joining=[
            dict(names=list(d.names), floor_lat=d.floor_lat, floor_lon=d.floor_lon,
                 area_fraction=round(d.area_fraction, 4), global_layer_m=round(d.gel_m, 1)) for d in parts]))
    stats = dict(height_min_m=float(height.min()), height_max_m=float(height.max()),
                 geoid_min_m=float(geoid.min()), geoid_max_m=float(geoid.max()),
                 fine_pixels_per_degree=fine_ppd, tree_pixels_per_degree=tree_ppd, geoid_degree=geoid_degree)
    return dict(levels=rows, merges=merges, stats=stats)


def land_sea(fraction, resolution_deg=1.0, fine_ppd=16, geoid_degree=200):
    """Water fraction, mean ground height and mean water depth on a coarse lat-lon grid for a
    common water level covering `fraction` of the Moon (area-weighted block averages)."""
    height, _, grid = tp.height_above_geoid(fine_ppd, geoid_degree)
    area = grid.cell_area_m2
    levels = np.linspace(height.min(), height.max(), 40001)
    frac, vol = bs.level_curve(height, area, levels)
    level = float(np.interp(fraction, frac, levels))
    k = int(round(resolution_deg * fine_ppd))
    if (180 * fine_ppd) % k or k < 1:
        raise ValueError('Resolution must divide the grid')
    nlat, nlon = height.shape[0] // k, height.shape[1] // k
    blocks = lambda a: a.reshape(nlat, k, nlon, k)
    a = blocks(area)
    wet = blocks(height < level)
    depth = blocks(np.maximum(level - height, 0.0))
    area_sum = a.sum(axis=(1, 3))
    water_fraction = (a * wet).sum(axis=(1, 3)) / area_sum
    mean_height = (a * blocks(height)).sum(axis=(1, 3)) / area_sum
    wet_area = (a * wet).sum(axis=(1, 3))
    mean_depth = np.where(wet_area > 0, (a * depth).sum(axis=(1, 3)) / np.maximum(wet_area, 1e-30), 0.0)
    lat = 90.0 - (np.arange(nlat) + 0.5) * resolution_deg
    lon = (np.arange(nlon) + 0.5) * resolution_deg
    covered = float((water_fraction * area_sum).sum() / area_sum.sum())
    return dict(level_m=level, covered_fraction=covered, global_layer_m=float(np.interp(level, levels, vol) / bs.AREA),
                lat_deg=lat, lon_deg=lon, water_fraction=water_fraction, mean_height_m=mean_height,
                mean_water_depth_m=mean_depth)


def write_masks(fractions, resolution_deg, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    for f in fractions:
        m = land_sea(f, resolution_deg)
        meta = dict(schema='terluna.geography.land-sea/1', producer=producer(), evidence=EVIDENCE,
                    reading_rule=('water_fraction is the share of each cell below one common water level above the '
                                  'geoid; mean_height_m is ground height above the geoid; mean_water_depth_m averages '
                                  'over the wet part only. Longitudes are east, 0-360; latitudes north first.'),
                    level_m=m['level_m'], covered_fraction=m['covered_fraction'], global_layer_m=m['global_layer_m'],
                    resolution_deg=resolution_deg)
        name = out / f'land_sea_{round(100 * f)}pct_{resolution_deg:g}deg.npz'
        np.savez_compressed(name, metadata=json.dumps(meta), lat_deg=m['lat_deg'], lon_deg=m['lon_deg'],
                            water_fraction=m['water_fraction'].astype(np.float32),
                            mean_height_m=m['mean_height_m'].astype(np.float32),
                            mean_water_depth_m=m['mean_water_depth_m'].astype(np.float32))
        print(f"{name.name}: level {m['level_m']/1e3:+.2f} km, covered {m['covered_fraction']:.1%}, "
              f"layer {m['global_layer_m']:.0f} m")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    parser.add_argument('--masks', type=float, nargs='*', help='write land-sea masks for these water fractions')
    parser.add_argument('--resolution', type=float, default=1.0)
    args = parser.parse_args(argv)
    if fetch_inputs.restore()['status'] != 'PASS':
        raise SystemExit('BLOCKED: topography inputs missing; run python -m geography.fetch_inputs --download')
    if args.masks:
        write_masks(args.masks, args.resolution, HERE / 'products')
        return
    result = compute()
    args.out.mkdir(parents=True, exist_ok=True)
    with open(args.out / 'water_levels.csv', 'w', newline='') as handle:
        fields = ['area_fraction', 'level_m', 'global_layer_m', 'mass_kg', 'mean_depth_m', 'seas_over_0_1pct']
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(result['levels'])
    meta = dict(schema=SCHEMA, producer=producer(), evidence=EVIDENCE,
                units=dict(level='m above the geoid (zero-mean equipotential)', global_layer='m of water over the whole Moon',
                           mass='kg at 1000 kg/m3'),
                reading_rule='Levels and layers are hydrostatic storage; do not read them as predicted shorelines.',
                **result)
    (args.out / 'water_inventory.json').write_text(json.dumps(meta, indent=1) + '\n')
    for r in result['levels']:
        print(f"{r['area_fraction']:6.1%} of the Moon: level {r['level_m']/1e3:+6.2f} km, layer {r['global_layer_m']:6.0f} m, "
              f"mean depth {r['mean_depth_m']:5.0f} m, {r['seas_over_0_1pct']} seas > 0.1%")


if __name__ == '__main__':
    main()
