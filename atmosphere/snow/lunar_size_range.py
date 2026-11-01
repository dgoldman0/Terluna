"""Lunar snow crystal sizes as a range over water supply and crystal number.

Plan step 3.  The steady lunar column of mcsnow_lunar_growth (the 301 K
reference ladder from -15 C to the freezing level, 16.28 km) shares an imposed
water supply Q among a steady incoming crystal number flux F.  Each crystal may
deposit at most Q/F over the full descent, in proportion to the distance it has
fallen, while the frozen improved McSnow component sets how fast it can grow,
how it divides its ice between axes and how fast it falls.  The grid spans light
to heavy snow and three decades of crystal number; every point is also run at
the two combined edges of the crystal model's Earth error band
(earth_error_band.py), so each size carries the model's own uncertainty.  The
results are single crystals arriving at the freezing level; aggregation is plan
step 4.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

from atmosphere.snow import earth_error_band as band

HERE = Path(__file__).resolve().parent
ROOT = band.ROOT
PRODUCT = HERE / 'results/lunar_size_range.json'
WATER_MM_H = (.1, .3, 1., 3., 10.)
FLUX_FACTORS = (1., .3, .1, .03, .01, .003, .001)


def reference():
    entry = band.column_entries()['allowance_F']
    return entry, entry['configuration']['number_flux_m2_s']


def edge_factors():
    product = json.loads(band.PRODUCT.read_text())
    table = product['column']['allowance_F']
    return {name: tuple(table[name]['factors']) for name in ('baseline', 'combined_small', 'combined_large')}


def run():
    entry, f0 = reference()
    edges = edge_factors()
    rows = []
    for water, factor in itertools.product(WATER_MM_H, FLUX_FACTORS):
        case = dict(entry, configuration=dict(entry['configuration'], water_supply_kg_m2_s=water / 3600.,
                                              number_flux_m2_s=f0 * factor))
        sizes = {}
        for name, factors in edges.items():
            label = f'range_q{water:g}_f{factor:g}__{name}'
            r = band.run_column(label, case, factors)
            if not r['complete']:
                raise ValueError('Incomplete column run: ' + label)
            sizes[name] = dict(diameter_m=r['final']['diameter_m'], mass_kg=r['final']['mass_kg'],
                               descent_h=r['final']['elapsed_s'] / 3600,
                               deposited_kg=r['final']['integrated_deposition_kg'],
                               allowed_kg=r['final']['allocated_water_kg'])
        rows.append(dict(water_mm_h=water, flux_factor=factor, number_flux_m2_s=f0 * factor,
                         water_per_crystal_kg=water / 3600. / (f0 * factor) , **sizes))
    unlimited = json.loads(band.PRODUCT.read_text())['column']['unlimited_saturation']
    storm = json.loads(band.PRODUCT.read_text())['storm']
    product = dict(
        schema='terluna.atmosphere.snow-lunar-size-range/1',
        evidence=('Single-crystal sizes at the freezing level from the frozen improved McSnow component in the steady '
                  'prescribed lunar column, with the combined edges of its Earth error band; no aggregation, '
                  'depletion feedback or empirical lunar validation.'),
        column=dict(case=entry['name'], reference_number_flux_m2_s=f0, depth_m=entry['configuration']['cloud_depth_m']),
        edges={k: list(v) for k, v in edges.items()}, grid=rows,
        unlimited_saturation={k: unlimited[k] for k in ('baseline', 'combined_small', 'combined_large')},
        archived_storm={k: storm['archived'][k] for k in ('baseline', 'combined_small', 'combined_large')},
        saturated_storm={k: storm['saturated'][k] for k in ('baseline', 'combined_small', 'combined_large')},
        inputs=dict(earth_error_band_sha256=band.native.sha(band.PRODUCT)),
        producer_sha256=band.native.sha(Path(__file__)), raw_directory=str((band.CACHE / 'column_runs').relative_to(ROOT)),
        limits=['Each crystal receives water in proportion to the distance it falls, from a supply that does not '
                'respond to the crystals; real storms deplete and replenish vapour unevenly.',
                'The column is a straight prescribed profile at water saturation with no vertical air motion.',
                'Single crystals; the flakes that reach the ground are aggregates of them (plan step 4).'])
    PRODUCT.write_text(json.dumps(product, separators=(',', ':'), allow_nan=False) + '\n')
    return product


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    product = run()
    for r in product['grid']:
        b, s, l = (r[k]['diameter_m'] * 1e3 for k in ('baseline', 'combined_small', 'combined_large'))
        print(f"Q={r['water_mm_h']:5g} mm/h  F={r['flux_factor']:6g} F0  m/crystal={r['water_per_crystal_kg']*1e6:9.4f} mg"
              f"  D={b:8.3f} mm  band {s:8.3f}-{l:8.3f}  descent {r['baseline']['descent_h']:6.1f} h")


if __name__ == '__main__':
    main()
