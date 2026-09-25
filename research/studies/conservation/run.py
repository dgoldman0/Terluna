#!/usr/bin/env python3
"""Conservation study: heritage register, scientific targets, gates and community settings.

Run from the repository root after the geography atlas:
    python -m geography.atlas
    python -m research.studies.conservation.run [--out DIR]
The default output is research/studies/conservation/results (small tables kept in Git).

The register joins two layers. The written evaluation in sites.json (values, dependence
on place, proposed treatment and whom to consult) is human judgment. This runner adds
each site's computed situation at the scenario water share: height, depth, pressure,
distances and the water body it lies in.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if __package__ in (None, ''):  # executed as a file path rather than with -m
    sys.path.insert(0, str(ROOT))
from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY
from atmosphere.thermal_column import ColumnConfig
from geography import atlas as ga, nomenclature as nm
from research.studies.conservation import screens as sc

HERE = Path(__file__).resolve().parent
PRODUCT = ROOT / 'geography' / 'results' / 'atlas.json'
SHORE_BAND_M = 150.0     # sites less than this above the sea are shore sites
LAND_MARGIN_M = 50.0     # land counts as land this far above the sea when measuring distances
SURFACE_PA = ColumnConfig().surface_pressure_pa
GCM_BOND_ALBEDO = 0.287  # planetary albedo of the ExoPlaSim design case, run A (climate/gcm/README.md)
BUILD_UP_YEARS = 500.0   # five-century construction scenario of the engineering seed
DESIGN_DEPTH_M = 700.0   # Tranquility community design depth, covering the shoreline's uncertainty


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_product(path=PRODUCT):
    atlas = json.loads(Path(path).read_text())
    if atlas.get('schema') != ga.SCHEMA:
        raise ValueError(f'{path}: expected schema {ga.SCHEMA}')
    grid_path = ROOT / 'geography' / 'products' / f"atlas_{round(100 * atlas['share'])}pct_4ppd.npz"
    if not grid_path.is_file():
        raise SystemExit(f'BLOCKED: {grid_path.name} missing; run python -m geography.atlas')
    arrays = dict(np.load(grid_path))
    if json.loads(str(arrays['metadata']))['schema'] != ga.GRID_SCHEMA:
        raise ValueError(f'{grid_path}: expected schema {ga.GRID_SCHEMA}')
    return atlas, arrays, dict(atlas_json=digest(path), atlas_grid=digest(grid_path))


def classify(height_m, level_m, shore_band_m=SHORE_BAND_M):
    if height_m < level_m:
        return 'submerged'
    return 'shore' if height_m < level_m + shore_band_m else 'dry'


def nearest_km(lat, lon, mask, arrays):
    la, lo = np.meshgrid(arrays['lat_deg'], arrays['lon_deg'], indexing='ij')
    d = ga.great_circle_km(lat, lon, la[mask], lo[mask])
    return float(d.min()) if d.size else float('nan')


def body_names(atlas):
    return {b['label']: (b['water_names'][0] if b['water_names'] else 'crater lake') for b in atlas['bodies']}


def register(sites, atlas, arrays):
    level = atlas['sea_level_m']
    H = arrays['height_m']
    lat = np.array([s['lat'] for s in sites])
    lon = np.array([s['lon'] for s in sites])
    heights, sphere = ga.point_heights(lat, lon)
    names = body_names(atlas)
    grid = ga.tp.Grid(arrays['lat_deg'], arrays['lon_deg'], 4)
    land, water = H >= level + LAND_MARGIN_M, H < level
    rows = []
    for s, h, hs in zip(sites, heights, sphere):
        status = classify(float(h), level)
        i, j = ga.cell_index(grid, s['lat'], s['lon'])
        label = int(arrays['water_label'][i, j])
        row = dict(id=s['id'], name=s['name'], year=s['year'], kind=s['kind'], lat=s['lat'], lon=s['lon'],
                   height_above_geoid_m=round(float(h)), height_above_sphere_m=round(float(hs)), status=status)
        if status == 'submerged':
            depth = level - float(h)
            row.update(depth_m=round(depth), pressure_atm=round(sc.pressure_at_depth(depth, SURFACE_PA) / sc.ATM, 1),
                       earth_dive_equivalent_m=round(sc.earth_dive_equivalent(depth, SURFACE_PA)),
                       water_body=names.get(label, 'crater lake') if label else 'crater lake',
                       nearest_land_km=round(nearest_km(s['lat'], s['lon'], land, arrays)))
        else:
            row.update(above_sea_m=round(float(h) - level), nearest_water_km=round(nearest_km(s['lat'], s['lon'], water, arrays)))
        row.update({k: s[k] for k in ('marks', 'values', 'place_dependence', 'survives', 'matters_to', 'treatment_type',
                                      'proposed_treatment', 'reason', 'status', 'consult', 'open_questions')}
                   | dict(evaluation_status=s['status']))
        row['status'] = status
        rows.append(row)
    return rows


def targets(items, features, atlas):
    by = nm.by_name(features)
    level = atlas['sea_level_m']
    found = [t for t in items if t['name'] in by]
    lat = np.array([by[t['name']]['lat'] for t in found])
    lon = np.array([by[t['name']]['lon'] for t in found])
    heights, _ = ga.point_heights(lat, lon)
    rows = []
    for t, h in zip(found, heights):
        f = by[t['name']]
        under = float(h) < level
        zone = 'B, polar cold trap' if t['cold_trap'] else ('A, future seafloor' if under else 'C, land')
        deadline = ('before the first gas' if t['cold_trap'] else
                    'before flooding' if under else 'before the ground is wetted')
        rows.append(dict(name=t['name'], code=f['code'], lat=round(f['lat'], 2), lon=round(f['lon'], 2),
                         diameter_km=round(f['diameter_km'], 1), record=t['record'],
                         height_relative_to_sea_m=round(float(h) - level),
                         becomes=('lake' if under and t['cold_trap'] else 'seafloor or lake floor' if under else 'land'),
                         zone=zone, sample=deadline))
    missing = [t['name'] for t in items if t['name'] not in by]
    return rows, missing


def gates(reg, atlas):
    tranquility = next(r for r in reg if r['id'] == 'apollo-11')

    def community(depth):
        net = sc.pressure_at_depth(depth, SURFACE_PA) - sc.ATM
        return dict(depth_m=round(depth), net_inward_atm=round(net / sc.ATM, 1),
                    steel_shell_m_at_100m_radius=round(sc.dome_shell_m(net, 100.0, 200e9, 200e6), 2),
                    lunar_concrete_shell_m_at_100m_radius=round(sc.dome_shell_m(net, 100.0, 30e9, 20e6), 2),
                    floor_uplift_basalt_equivalent_m=round(sc.rock_equivalent_m(net)),
                    access_tower_earth_equivalent_m=round(sc.earth_tower_equivalent_m(depth)))
    area = 4 * np.pi * MOON_RADIUS ** 2
    land = atlas['main_land_share'] + atlas['island_share']
    heat_leak_w_m2, carnot_share = 1.0, 0.15
    cop = 40.0 / (295.0 - 40.0) * carnot_share
    return dict(
        about='Screening estimates (screens.py) that date each loss and size the heritage works; domain models supersede them.',
        surface_pressure_pa=SURFACE_PA, air_column_kg_m2=round(SURFACE_PA / MOON_SURFACE_GRAVITY),
        water_depth_per_atm_m=round(sc.ATM / (sc.WATER_DENSITY * MOON_SURFACE_GRAVITY), 1),
        dusty_stage=[dict(sc.saltation_threshold(p * sc.ATM),
                          year_of_linear_build_up=round(BUILD_UP_YEARS * p * sc.ATM / SURFACE_PA))
                     for p in (0.01, 0.03, 0.1, 0.3, 0.6, 1.2)],
        cold_trap_frost_pa=dict(N2_37K=sc.frost_pressure('N2', 37.0), N2_40K=sc.frost_pressure('N2', 40.0),
                                O2_40K=sc.frost_pressure('O2', 40.0)),
        warming_depth_m={f'{kappa:.0e} m2/s': {f'{y} yr': round(sc.warming_depth_m(y, kappa)) for y in (10, 100, 1000)}
                         for kappa in (1e-7, 1e-6)},
        grain_rim_years=dict(laboratory_60_to_200nm=[round(sc.dissolution_years(t, 3e-11), 1) for t in (60e-9, 200e-9)],
                             field_slowdown_factor=[1e2, 1e5]),
        runoff_scaling=sc.runoff_scaling(),
        raindrop_energy_ratio={f'{d} mm': round(sc.raindrop_energy_ratio(d / 1e3, surface_pa=SURFACE_PA), 3)
                               for d in (0.5, 1, 2, 4)},
        tranquility_community=dict(at_scenario_share=community(tranquility['depth_m']),
                                   at_design_depth=community(DESIGN_DEPTH_M)),
        sealed_reference_area=dict(air_load_t_m2=round(SURFACE_PA / MOON_SURFACE_GRAVITY / 1e3, 1),
                                   rock_equivalent_m=round(sc.rock_equivalent_m(SURFACE_PA)),
                                   steel_shell_m_at_150m_radius=round(sc.dome_shell_m(SURFACE_PA, 150.0, 200e9, 200e6), 2),
                                   cold_patch_power_mw_per_km2=round(heat_leak_w_m2 * 1e6 / cop / 1e6, 1),
                                   cold_patch_assumptions='1 W/m2 heat leak through vacuum insulation; 15% of Carnot at 40 K'),
        sediment_return={f'{rate} m/Myr': dict(zip(('ideal_power_gw', 'rock_gt_per_year'),
                                                   (round(p / 1e9, 2), round(m / 1e12, 2))))
                         for rate in (10, 50, 100) for p, m in [sc.sediment_return_watts(rate, land, 3000.0, area)]},
        civilization_power_w=dict(K_1_17=sc.kardashev_watts(1.17)),
        moonlight_on_earth=dict(sc.moonlight_factors(GCM_BOND_ALBEDO), bond_albedo=GCM_BOND_ALBEDO,
                                source='ExoPlaSim design case, run A (climate/gcm/README.md)'),
    )


def community_settings(reg, atlas, target_rows):
    communities = [r for r in reg if r['status'] == 'submerged' and r['treatment_type'].startswith('undersea community')]
    heritage = [dict(site=r['name'], depth_m=r['depth_m'], nearest_land_km=r['nearest_land_km'],
                     treatment=r['treatment_type'], evaluation=r['evaluation_status']) for r in communities]
    coasts = [dict(beside=r['name'], coast_km=r['nearest_land_km']) for r in communities if r['nearest_land_km'] <= 20]
    islands = [dict(summit=i['summit'], area_km2=i['area_km2'], summit_m=i['summit_m'], near=i['features'])
               for i in atlas['islands'] if i['area_km2'] >= 2000]
    far = [dict(body=(b['water_names'] or ['crater lake'])[0], area_km2=b['area_km2'], mean_depth_m=b['mean_depth_m'])
           for b in atlas['bodies'] if abs(b['centroid'][1]) > 90]
    polar = [dict(name=t['name'], height_relative_to_sea_m=t['height_relative_to_sea_m'])
             for t in target_rows if t['zone'].startswith('B') and t['becomes'] == 'land']
    return dict(about=('Candidate settings from geography alone. People choose community sites, with climate, ecology '
                       'and the wishes of those who would live there.'),
                submerged_heritage_communities=heritage, related_coastal_communities=coasts, islands=islands,
                far_side_waters=far, dry_polar_highs=polar)


def write_csv(path, rows, fields):
    with open(path, 'w', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def run(out=HERE / 'results'):
    out = Path(out).resolve()
    if (ROOT / 'immersion') in out.parents or out == ROOT / 'immersion':
        raise ValueError('The conservation study writes outside immersion/')
    atlas, arrays, product_hashes = load_product()
    sites = json.loads((HERE / 'sites.json').read_text())
    items = json.loads((HERE / 'targets.json').read_text())['targets']
    features = nm.features()
    reg = register(sites['sites'], atlas, arrays)
    target_rows, missing = targets(items, features, atlas)
    out.mkdir(parents=True, exist_ok=True)
    header = dict(schema='terluna.research.conservation-register/1', share=atlas['share'],
                  sea_level_m=atlas['sea_level_m'], evaluation=sites['evaluation'],
                  reading_rule=('Computed fields (heights, status, depth, pressure, distances, water body) come from '
                                'the geography atlas at the scenario share; the evaluation fields are written judgment.'))
    (out / 'heritage_register.json').write_text(json.dumps(dict(header, sites=reg), indent=1, ensure_ascii=False) + '\n')
    write_csv(out / 'heritage_register.csv', reg,
              ['id', 'name', 'year', 'kind', 'status', 'depth_m', 'pressure_atm', 'above_sea_m', 'nearest_land_km',
               'nearest_water_km', 'water_body', 'treatment_type', 'evaluation_status'])
    write_csv(out / 'science_targets.csv', target_rows,
              ['name', 'code', 'lat', 'lon', 'diameter_km', 'height_relative_to_sea_m', 'becomes', 'zone', 'sample', 'record'])
    (out / 'gates.json').write_text(json.dumps(gates(reg, atlas), indent=1) + '\n')
    (out / 'community_settings.json').write_text(
        json.dumps(community_settings(reg, atlas, target_rows), indent=1, ensure_ascii=False) + '\n')
    files = ['run.py', 'screens.py', 'sites.json', 'targets.json']
    checkpoint = dict(product=dict(schema=ga.SCHEMA, **product_hashes), share=atlas['share'],
                      sea_level_m=atlas['sea_level_m'], surface_pressure_pa=SURFACE_PA,
                      study_files={f: digest(HERE / f)[:16] for f in files},
                      gazetteer_sha256=digest(ga.fetch_inputs.path(nm.ARCHIVE)), targets_missing=missing)
    (out / 'checkpoint.json').write_text(json.dumps(checkpoint, indent=1) + '\n')
    return reg, target_rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    reg, target_rows = run(args.out)
    for r in reg:
        where = (f"{r['depth_m']:>5} m under, {r['pressure_atm']:>5} atm, land {r['nearest_land_km']:>4} km"
                 if r['status'] == 'submerged' else f"{r['above_sea_m']:>5} m above, water {r['nearest_water_km']:>4} km")
        print(f"{r['name'][:42]:42s} {r['status']:9s} {where}  [{r['evaluation_status']}]")
    print(f"{len(target_rows)} scientific targets placed")


if __name__ == '__main__':
    main()
