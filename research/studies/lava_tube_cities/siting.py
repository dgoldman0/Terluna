"""Reproduce city-scale arithmetic and a coarse water screen for literature leads.

Source interpretation is recorded in siting_inputs.json. This calculation does not
infer caves, rank confidence numerically, or approve locations.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from shared.constants import MOON_RADIUS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.lava-tube-city-siting/1'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_bound(spec, key):
    binding = spec['inputs'][key]
    path = ROOT / binding['path']
    if sha256(path) != binding['sha256']:
        raise ValueError(f'Input hash mismatch: {key}')
    if path.suffix == '.npz':
        with np.load(path, allow_pickle=False) as z:
            value = {k: z[k] for k in ('lat_deg', 'lon_deg', 'height_m', 'water_label')}
            value['metadata'] = json.loads(str(z['metadata']))
        actual_schema = value['metadata']['schema']
    else:
        value = json.loads(path.read_text())
        actual_schema = value['schema']
    if actual_schema != binding['schema']:
        raise ValueError(f'Input schema mismatch: {key}')
    return value


def atlas_point(atlas, latitude, longitude):
    """Containing cell; points on boundaries go south/east, longitude wraps."""
    lat = atlas['lat_deg']
    lon = atlas['lon_deg']
    if not -90 <= latitude <= 90:
        raise ValueError('Latitude outside the globe')
    dy = float(lat[0] - lat[1])
    dx = float(lon[1] - lon[0])
    north = float(lat[0]) + dy / 2
    west = float(lon[0]) - dx / 2
    i = int(np.clip(math.floor((north - latitude) / dy), 0, len(lat) - 1))
    j = math.floor(((longitude - west) % 360) / dx) % len(lon)
    sea = atlas['metadata']['sea_level_m']
    neighbours = atlas['height_m'][max(0, i - 1):min(len(lat), i + 2)][:, [(j - 1) % len(lon), j, (j + 1) % len(lon)]]
    level = float(atlas['height_m'][i, j]) - sea
    wet = int(atlas['water_label'][i, j]) > 0
    return {
        'cell_index': [i, j],
        'cell_centre_deg': [float(lat[i]), float(lon[j])],
        'ground_above_common_sea_m': round(level, 3),
        'common_level_water': wet,
        'water_above_cell_ground_m': round(max(-level, 0), 3) if wet else 0,
        'neighbour_ground_above_sea_range_m': [round(float(neighbours.min()) - sea, 3), round(float(neighbours.max()) - sea, 3)],
        'north_south_cell_km': MOON_RADIUS / 1000 * math.radians(dy),
    }


def population_scales(towns, design):
    floor = design['population_cases']['representative']['floor_m2_per_resident']
    rows = []
    for size in towns['population_cases']:
        people = towns['communities'] * size
        rows.append({'residents_per_town': size, 'residents_per_city': people,
                     'residents_four_cities': 4 * people,
                     'residents_cities_plus_towns': 5 * people,
                     'ordinary_floor_m2_per_city': people * floor})
    return rows


def run():
    spec = json.loads((HERE / 'siting_inputs.json').read_text())
    atlas = load_bound(spec, 'atlas')
    towns = load_bound(spec, 'town_construction_inputs')
    design = load_bound(spec, 'town_design')
    old = load_bound(spec, 'land_cities')
    old_rows = {r['key']: r for r in old['tubes']['on_the_atlas']['rows']}
    rows = []
    for candidate in spec['candidates']:
        row = dict(candidate)
        row['coarse_water_screen'] = atlas_point(atlas, candidate['latitude_deg'], candidate['longitude_east_deg'])
        if candidate.get('historical_key'):
            prior = old_rows[candidate['historical_key']]
            row['historical_point_comparison'] = {
                'latitude_deg': prior['lat'], 'longitude_east_deg': prior['lon'],
                'ground_above_sea_m': prior['ground_above_sea_m'],
                'class': prior['class'],
                'same_coordinate': abs(prior['lat'] - candidate['latitude_deg']) < 1e-6 and abs(prior['lon'] - candidate['longitude_east_deg']) < 1e-6,
                'rule': 'Historical point and new coarse cell are different spatial samples; retain both. A different coordinate cannot supply the new candidate elevation.'}
        rows.append(row)
    g = spec['illustrative_geometry']
    gross = g['width_m'] * g['length_m'] * g['levels']
    allocated = gross * g['ordinary_floor_share']
    files = ['research/studies/lava_tube_cities/siting.py',
             'research/studies/lava_tube_cities/siting_inputs.json',
             'research/studies/lava_tube_cities/sources.json']
    return {
        'schema': SCHEMA, 'date': spec['date'], 'evidence': spec['evidence'],
        'inputs': spec['inputs'],
        'producer': {'files': {p: sha256(ROOT / p) for p in files}, 'constants': constants_used(files)},
        'population_scales': population_scales(towns, design),
        'illustrative_geometry': dict(g, gross_floor_m2=gross, ordinary_floor_m2=allocated,
                                      residents_at_inherited_floor_allowance=allocated / design['population_cases']['representative']['floor_m2_per_resident']),
        'candidates': rows,
        'assessment': spec['assessment'],
        'reading_rules': spec['reading_rules'],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'research/runs/lava_tube_city_siting.json')
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    for row in result['candidates']:
        q = row['coarse_water_screen']
        print(row['key'], q['ground_above_common_sea_m'], 'm above common sea;', row['investigation_priority'])


if __name__ == '__main__':
    main()
