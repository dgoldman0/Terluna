"""Inspect the actual exported sites against the source census and estate geometry."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from shared.constants import MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def verify():
    product_path = HERE / 'results/buildout.json'
    product = json.loads(product_path.read_text())
    export = Path(product['selected_sites']['path'])
    assert digest(export) == product['selected_sites']['sha256']
    with gzip.open(export, 'rt', newline='') as f:
        selected = list(csv.DictReader(f))
    ids = {r['crater_id'] for r in selected}
    assert len(ids) == len(selected) == product['programme']['eventual_sites']
    initial = np.array([r['initial_vacuum'] == '1' for r in selected])
    assert initial.sum() == product['programme']['initial_sites']
    data = np.array([[float(r[k]) for k in ('latitude_today_deg', 'longitude_today_deg', 'diameter_m', 'depth_m')]
                     for r in selected])
    by_id = {r['crater_id']: i for i, r in enumerate(selected)}
    matched = 0
    with gzip.open(product['producer']['inputs']['candidates']['path'], 'rt', newline='') as f:
        for row in csv.DictReader(f):
            if row['crater_id'] not in ids:
                continue
            i = by_id[row['crater_id']]
            assert row['confidence_level'] == 'high'
            assert row['adjacent_water_cell'] == '0'
            expected = [float(row[k]) for k in ('latitude_today_deg', 'longitude_today_deg', 'diameter_m', 'depth_m')]
            assert np.array_equal(data[i], expected)
            matched += 1
    assert matched == len(selected)
    design = json.loads((ROOT / 'habitation/crater_gardens/design.json').read_text())
    totals = {}
    for name, p in design['population_cases'].items():
        # Independent scalar contribution per site, then math.fsum for accumulation.
        per_site = [math.pi * (d / 2) ** 2 * .75 * p['non_lake_ground_building_share'] * p['mean_storeys'] / p['floor_m2_per_resident'] for d in data[:, 2]]
        total = math.fsum(per_site)
        initial_total = math.fsum(v for v, keep in zip(per_site, initial) if keep)
        assert math.isclose(total, product['eventual_by_estate_factor']['2.0']['architectural_population'][name], rel_tol=1e-12)
        assert math.isclose(initial_total, product['initial_clustered']['architectural_population'][name], rel_tol=1e-12)
        if name == 'representative':
            assert all(500 <= value <= 1000 for value, keep in zip(per_site, initial) if keep)
        totals[name] = {'eventual': total, 'initial': initial_total}
    lat, lon = np.deg2rad(data[:, 0]), np.deg2rad(data[:, 1])
    unit = np.column_stack((np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)))
    radius = 2 * np.sin(data[:, 2].max() / MOON_RADIUS)
    pairs = cKDTree(unit).query_pairs(radius, output_type='ndarray')
    a, b = unit[pairs[:, 0]], unit[pairs[:, 1]]
    distance = MOON_RADIUS * np.arctan2(np.linalg.norm(np.cross(a, b), axis=1), np.einsum('ij,ij->i', a, b))
    margins = distance - data[pairs[:, 0], 2] - data[pairs[:, 1], 2]
    assert np.all(margins >= -1e-6)
    groups_path = Path(product['initial_groups']['path'])
    assert digest(groups_path) == product['initial_groups']['sha256']
    groups = json.loads(groups_path.read_text())
    for number, group in enumerate(groups):
        members = [i for i, row in enumerate(selected) if row['initial_group'] == str(number)]
        seed = by_id[group['seed_id']]
        assert initial[seed] and len(members) == group['count']
        dots = unit[members] @ unit[seed]
        distances = MOON_RADIUS * np.arctan2(np.linalg.norm(np.cross(unit[members], unit[seed]), axis=1), dots)
        assert np.all(distances <= 25000 + 1e-6)
    record = {'schema': 'terluna.research.vacuum-crater-verification/1',
              'evidence': 'Independent inspection of the exported candidate identities, dimensions, population sums, initial group bounds and spherical circular-estate separation. No physical site validation.',
              'result': 'PASS', 'product_sha256': digest(product_path), 'export_sha256': digest(export),
              'verifier_sha256': digest(Path(__file__)), 'matched_source_rows': matched,
              'initial_rows': int(initial.sum()), 'initial_groups': len(groups),
              'nearby_pairs_checked': len(pairs), 'minimum_estate_clearance_m': float(margins.min()) if len(margins) else None,
              'architectural_population': totals}
    (HERE / 'results/verification.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    verify()
