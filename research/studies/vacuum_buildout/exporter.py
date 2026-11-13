"""Export a bounded engineering response grid from the pinned provisioning checkout.

Run this recipe in a separate process with --provisioning-root. Consumers read
inputs/engineering.json; they do not import the provisioning implementation.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
COMMIT = '146a558275f486c9fd9a5cf79908b7c6fdf21888'
STRATEGY_SNAPSHOT = ROOT / 'research/studies/complex_returns/inputs/strategy_profiles.json'
STRATEGY_SHA = '68de393762a3411e0fc3fb2956c8d16c936c46019db70995cc3ccd5654449a9c'
SOURCE_STUDY = 'research/studies/initial_industrial_system/'
GRID = (0., .01, .02, .03, .04, .05, .06)
PROFILE_MAP = {
    'raw_resource_price': 'source_refined',
    'commodity_refined': 'commodity_mixed',
    'specialty': 'specialty_mixed',
    'premium_lunar': 'goods_heavy',
    'propellant_logistics': 'logistics_1Mt',
    'compute_digital': 'computing_5GW',
    'services': 'domestic_services',
    'diversified': 'balanced',
    'prior_manufacturing_mix': 'current_large',
    'domestic_services': 'domestic_services',
}
PHYSICAL = (
    'lunar_bulk_kg_year', 'lunar_metal_kg_year', 'lunar_raw_kg_year',
    'lunar_handling_kg_year', 'lunar_process_W', 'resident_mean_W',
    'source_metal_kg_year', 'source_virgin_metal_kg_year', 'source_feed_kg_year',
    'source_water_kg_year', 'source_nitrogen_kg_year', 'source_retained_tails_kg_year',
    'source_process_W', 'source_heat_W', 'depot_W', 'LEO_kg_year',
    'Earth_launches_year', 'local_fabrication_kg_year', 'local_growth_kg_year',
    'retired_components_kg_year', 'recycled_to_source_kg_year',
    'recycled_from_LEO_kg_year', 'imports_by_destination_kg_year',
    'technical_resident_jobs', 'renewal_imports_by_destination_kg_year',
    'growth_imports_by_destination_kg_year', 'source_product_kg_year',
    'commercial_bulk_stock_kg', 'commercial_metal_stock_kg',
    'resident_water_kg_year', 'resident_air_kg_year', 'factory_Earth_feed_kg_year',
)
SCENARIOS = {
    'base': {},
    'lunar_stress': {'lunar': {'metal_yield': .05, 'bulk_yield': .8}},
    'source_stress': {'source': {'water_grade': .05, 'water_recovery': .5,
                                'nitrogen_grade': .0012, 'nitrogen_recovery': .5,
                                'metal_grade': .04, 'metal_recovery': .5}},
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_source(root, snapshot):
    found = subprocess.check_output(['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True).strip()
    if found != COMMIT:
        raise ValueError(f'Expected provisioning commit {COMMIT}; found {found}')
    paths = {row['path'] for row in snapshot['producer']['sources']}
    paths.update(row['path'] for row in snapshot['source_inputs'].values())
    paths.update((SOURCE_STUDY+'results/resource_inventory.json', SOURCE_STUDY+'resource_inventory.py'))
    # Record input manifests as well as the compiled products actually consumed.
    for name in ('industrial_inputs.json', 'transport_inputs.json', 'revenue_inputs.json'):
        manifest = json.loads((root / SOURCE_STUDY / name).read_text())
        for ref in manifest.get('inputs', {}).values():
            if isinstance(ref, dict) and 'path' in ref:
                paths.add(ref['path'])
    identities = {}
    for path in sorted(paths):
        current = sha(root / path)
        pinned = hashlib.sha256(subprocess.check_output(['git', '-C', str(root), 'show', f'{COMMIT}:{path}'])).hexdigest()
        if current != pinned:
            raise ValueError(f'Provisioning source differs from pinned commit: {path}')
        identities[path] = current
    return identities


def expand_architecture(original, source, spec):
    """Add rated equipment and its stock; retain all original route technologies."""
    c, src = deepcopy(original), deepcopy(source)
    factors = {}

    def scale(name, factor):
        factors[name] = factor
        for key in ('mass_kg', 'manufacture_USD2026'):
            c['assets'][name][key] *= factor
        if name in c['engineering']['counts']:
            c['engineering']['counts'][name] = round(c['engineering']['counts'][name] * factor)

    for name in ('excavators', 'haulers', 'bulk_trains', 'metal_trains', 'fabrication_cells',
                 'source_metal_trains', 'source_fabrication_cells'):
        scale(name, 2.)
    for key in ('bulk_kg_year', 'metal_kg_year', 'handling_kg_year', 'local_mean_W', 'source_metal_kg_year'):
        c['peak_flows'][key] *= 2.
    for key in ('source_mean_W', 'source_heat_W'):
        c['peak_flows'][key] *= 1.1
    for name in ('source_chemistry', 'source_solar', 'source_radiators'):
        scale(name, 1.1)
    landers = original['engineering']['counts']['lunar_landers']
    scale('lunar_landers', math.ceil(1.25*landers)/landers)
    scale('construction_and_service_robots', 1.15)
    src['settlement_inputs']['production']['floor_m2_year'] *= 1.15
    d = c['demand']
    firm = c['peak_flows']['local_mean_W']*spec['economics']['capacity_margin']/spec['lunar']['process_availability']
    extra = max(0., firm-d['old_construction_W'])
    total = d['old_total_W']+extra
    for name in ('additional_ground_receivers', 'additional_ground_radiator_controls'):
        scale(name, extra/original['engineering']['added_ground_firm_W'])
    scale('lunar_beam_constellation', total/original['engineering']['total_ground_firm_W'])
    c['engineering'].update(local_industry_firm_W=firm, added_ground_firm_W=extra,
                            total_ground_firm_W=total,
                            added_ground_radiator_m2=c['assets']['additional_ground_radiator_controls']['mass_kg']/d['radiator_controls_kg_m2'])
    for key in ('beam_PV_nameplate_W', 'beam_PV_area_m2', 'beam_converter_radiator_m2'):
        c['engineering'][key] *= factors['lunar_beam_constellation']
    c['engineering']['source_process_radiator_m2'] *= 1.1
    surface = sum(a['mass_kg'] for a in c['assets'].values() if a['destination']=='surface')
    c['peak_flows']['additional_bulk_stock_kg'] = surface*spec['lunar']['infrastructure_bulk_per_hardware_kg']
    c['peak_flows']['additional_metal_stock_kg'] = (
        surface*spec['lunar']['infrastructure_metal_per_hardware_kg']
        + c['assets']['additional_ground_radiator_controls']['mass_kg']/d['radiator_controls_kg_m2']*d['radiator_metal_kg_m2'])
    return c, src, factors


def compact(result, fraction):
    physical = {key: result['physical'][key] for key in PHYSICAL}
    physical['transport'] = {leg: {key: result['physical'][leg][key]
                                  for key in ('trips_year', 'space_supply_kg_year', 'burned_kg_year')}
                             for leg in ('Earth', 'Moon', 'source')}
    return dict(growth_fraction=fraction, physical=physical,
                cash_USD2026=result['cash_USD2026'], capital_formation=result['capital_formation'],
                capacity_use=result['capacity_use'], feasible=result['feasible'],
                binding_capacity=result['binding_capacity'])


def flattened(value, prefix=''):
    for key, item in value.items():
        path = f'{prefix}.{key}' if prefix else key
        if isinstance(item, dict):
            yield from flattened(item, path)
        elif isinstance(item, (float, int)) and not isinstance(item, bool):
            yield path, item


def interpolation_check(low, high, actual):
    a, b = dict(flattened(low)), dict(flattened(high))
    errors = []
    for key, value in flattened(actual):
        if key in ('growth_fraction', 'physical.Earth_launches_year'):
            continue
        estimated = (a[key]+b[key])/2
        relative = abs(estimated-value)/max(1., abs(value))
        errors.append((relative, key, estimated-value))
    worst = max(errors)
    capacities = [(key, abs((low['capacity_use'][key]+high['capacity_use'][key])/2-value))
                  for key, value in actual['capacity_use'].items()]
    return dict(growth_fraction=actual['growth_fraction'], maximum_relative_error=worst[0],
                worst_field=worst[1], worst_absolute_error=worst[2],
                maximum_capacity_ratio_absolute_error=max(value for _, value in capacities),
                actual_feasible=actual['feasible'])


def produce(root):
    if sha(STRATEGY_SNAPSHOT) != STRATEGY_SHA:
        raise ValueError('Frozen strategy-profile snapshot changed')
    snapshot = json.loads(STRATEGY_SNAPSHOT.read_text())
    identities = verify_source(root, snapshot)
    sys.path.insert(0, str(root))
    from research.studies.initial_industrial_system import evolution as model
    q, source, spec, original = model.read_inputs()
    expanded, expanded_source, factors = expand_architecture(original, source, spec)
    inflation = original['cost']['inflation_factor']
    architectures = {'original': (original, source), 'expanded': (expanded, expanded_source)}
    profiles = {}; checks = []; largest = 0.
    for name, source_name in PROFILE_MAP.items():
        args = snapshot['profiles'][source_name]['input_arguments']['physical']
        profiles[name] = dict(source_profile=source_name, input_arguments=args, architectures={})
        for architecture, (stock, src) in architectures.items():
            cases = {}
            for scenario, edits in SCENARIOS.items():
                scenario_spec = deepcopy(spec)
                for section, changes in edits.items():
                    scenario_spec[section].update(changes)
                def annual(g):
                    return model.annual(q, src, scenario_spec, stock, **dict(args,
                        growth_people=250000*g, new_plant_scale=g, commission_fraction=g))
                full = [annual(g) for g in GRID]
                rows = [compact(row, g) for row, g in zip(full, GRID)]
                if architecture=='original' and scenario=='base':
                    old = snapshot['profiles'][source_name]['engineering']
                    for field in ('capacity_use', 'physical', 'cash_USD2026', 'capital_formation'):
                        before, after = dict(flattened(old[field])), dict(flattened(full[0][field]))
                        if before.keys()!=after.keys() or any(not math.isclose(value, after[key], rel_tol=1e-12, abs_tol=1e-7) for key,value in before.items()):
                            raise ValueError(f'Original profile reproduction failed: {name}/{field}')
                for lo, hi in zip(rows[:-1], rows[1:]):
                    middle = (lo['growth_fraction']+hi['growth_fraction'])/2
                    check = interpolation_check(lo, hi, compact(annual(middle), middle))
                    largest = max(largest, check['maximum_relative_error'])
                    checks.append(dict(profile=name, architecture=architecture, scenario=scenario, **check))
                cases[scenario] = rows
            profiles[name]['architectures'][architecture] = dict(scenarios=cases)
    asset_changes = []
    for name, factor in factors.items():
        old, new = original['assets'][name], expanded['assets'][name]
        asset_changes.append(dict(name=name, multiplier=factor, destination=new['destination'],
            additional_mass_kg=new['mass_kg']-old['mass_kg'],
            additional_manufacture_USD2041=(new['manufacture_USD2026']-old['manufacture_USD2026'])*inflation,
            additional_inherited_installed_quote_USD2041=original['cost']['installed_assets_USD2026'][name]*(factor-1)*inflation))
    manufacture = sum(row['additional_manufacture_USD2041'] for row in asset_changes)
    installed = sum(row['additional_inherited_installed_quote_USD2041'] for row in asset_changes)
    build = {}
    for name in profiles:
        rows = profiles[name]['architectures']
        old, new = rows['original']['scenarios']['base'], rows['expanded']['scenarios']['base']
        def cash(row, key): return row['cash_USD2026'][key]
        build[name] = dict(
            added_steady_external_purchases_USD2041_year=(cash(new[0],'total_external_purchases')-cash(old[0],'total_external_purchases'))*inflation,
            added_steady_resource_provision_USD2041_year=(cash(new[0],'resource_provision')-cash(old[0],'resource_provision'))*inflation,
            added_modular_external_purchases_USD2041_per_scale=((cash(new[1],'total_external_purchases')-cash(new[0],'total_external_purchases'))-(cash(old[1],'total_external_purchases')-cash(old[0],'total_external_purchases')))/.01*inflation,
            added_modular_resource_provision_USD2041_per_scale=((cash(new[1],'resource_provision')-cash(new[0],'resource_provision'))-(cash(old[1],'resource_provision')-cash(old[0],'resource_provision')))/.01*inflation,
            added_modular_local_component_value_USD2041_per_scale=((new[1]['capital_formation']['local_growth_value_USD2026']-new[0]['capital_formation']['local_growth_value_USD2026'])-(old[1]['capital_formation']['local_growth_value_USD2026']-old[0]['capital_formation']['local_growth_value_USD2026']))/.01*inflation)
    inventory = json.loads((root / SOURCE_STUDY / 'results/resource_inventory.json').read_text())
    return dict(schema='terluna.research.vacuum-buildout-engineering/1',
        evidence='Pinned producer response with explicit added machinery, power and civil stock. A design sensitivity, not qualified hardware, source reserves, passenger capacity, funding or a settlement schedule.',
        source_commit=COMMIT, source_files=identities,
        producer=dict(path='research/studies/vacuum_buildout/exporter.py', sha256=sha(Path(__file__))),
        profile_definition_source=dict(path=str(STRATEGY_SNAPSHOT.relative_to(ROOT)),sha256=STRATEGY_SHA),
        reproduction='PYTHONDONTWRITEBYTECODE=1 python3 research/studies/vacuum_buildout/exporter.py --provisioning-root /path/to/pinned/terluna-provisioning',
        normalization=dict(population=250000, plant_scale=1., growth_fraction_grid=GRID, inflation_factor=inflation,
                           cash_units='USD2026; multiply by inflation_factor exactly once for constant USD2041'),
        strategy_names=[name for name in PROFILE_MAP if name!='domestic_services'],
        grade_scenarios=SCENARIOS, base_lunar_inputs=spec['lunar'], base_source_inputs=spec['source'],
        architecture={key:dict(assets=stock['assets'], engineering=stock['engineering'], peak_flows=stock['peak_flows'],
            floor_installation_m2_year=src['settlement_inputs']['production']['floor_m2_year'],
            route_installation_km_year=src['settlement_inputs']['production']['route_km_year']) for key,(stock,src) in architectures.items()},
        asset_changes=asset_changes,
        added_initial_cost=dict(manufacture_USD2041=manufacture, inherited_installed_quote_USD2041=installed,
            installed_with_preservation_and_contingency_USD2041=installed*(1+spec['economics']['source_preservation_fraction_of_supplier_capital'])*(1+spec['economics']['programme_contingency']),
            rule='Alternative initial deployment quote; includes inherited Earth-delivery rates, then separate preservation and contingency. Extra seed propellant and programme-specific qualification are unsized. Do not add both this quote and the modular replication bill for the same stock.'),
        incremental_modular_cost=build, profiles=profiles,
        interpolation=dict(method='Piecewise linear between adjacent saved growth fractions; verify interior controls and carry a capacity margin.',
            maximum_checked_relative_error=largest, checks=checks,
            discrete_fields=['physical.Earth_launches_year'],
            rule='Round interpolated launch counts upward. Ferries select payload/tank constraints, so the underlying response can be piecewise; interpolation controls are numerical checks, not a global linearity proof. Never extrapolate beyond the grid.'),
        near_earth_inventory_screen=dict(source=inventory['external_source'], inventory=inventory['inventory']),
        limits=[
            'The profile called raw_resource_price physically uses refined source goods, not shipped raw ore.',
            'All profiles retain the inherited service facilities; goods, computing and mission propellant differ by the recorded arguments. Commercial assets are commissioned at the same fraction as resident/industrial growth.',
            'Added structures, machinery renewal, alloy imports, civil bulk/metal, processing energy and freight enter the annual producer. Added rating alone never substitutes for added machinery.',
            'Floor installation is assumed proportional to the fifteen-percent additional construction robot stock; this is an unqualified productivity relationship.',
            'The producer technical-staff rule does not automatically increase staffing for selective machine additions. Expanded staffing must be screened and funded explicitly by the consumer; recorded jobs are the inherited lower bound.',
            'Humanities cultural/service-kit and premium-processing additions are absent here and must use the same remaining power and staffing headroom.',
            'The original component-lifetime and fixed per-product processing-energy assumptions persist in grade stress. Lower lunar grade therefore adds feed/handling demand without a separately qualified refining-energy law.',
            'Original sites, discrete power geometry, fleet reliability, cargo return and process readiness remain conditional; this response does not establish passenger migration capacity.',
            'No geographic lunar ore reserve or extraction depth is established. Catalogue sphere masses are incomplete hypothetical source stocks, not accessible permitted ore reserves.',
            'Initial deployment quotations, marginal expansion purchases, recurring replacement and local component value are separate accounts. Resource provisions are funding allowances, not realized extraction damage; local values are not automatically additional cash bills.'
        ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--provisioning-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=HERE/'inputs/engineering.json')
    args = parser.parse_args()
    product = produce(args.provisioning_root.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(product, indent=2, allow_nan=False)+'\n')
    print(args.output)
    print('Maximum checked interpolation relative error:', product['interpolation']['maximum_checked_relative_error'])


if __name__ == '__main__':
    main()
