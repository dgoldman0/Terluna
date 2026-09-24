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


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    if fetch_inputs.restore()['status'] != 'PASS':
        raise SystemExit('BLOCKED: topography inputs missing; run python -m geography.fetch_inputs --download')
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
