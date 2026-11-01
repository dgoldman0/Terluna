"""Recurring aerophyte products after host upkeep and population replacement.

Reads ecology products, never its model implementations. The result is a carbon
allocation screen for provisioning, not a crop-yield or profitability prediction.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = HERE / 'aerophytes_inputs.json'
SOURCES = HERE / 'aerophytes_sources.json'
OUT = HERE / 'results/aerophytes.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load():
    settings = json.loads(INPUTS.read_text())
    products = {}
    for name, identity in settings['products'].items():
        path = ROOT / identity['path']
        if sha(path) != identity['sha256']:
            raise ValueError(f'Input identity changed: {identity["path"]}')
        data = json.loads(path.read_text())
        if data.get('schema') != identity['schema']:
            raise ValueError(f'Unexpected schema: {identity["path"]}')
        products[name] = data
    return products, settings


def allocate(row, light_fraction, hazard, harvest_fraction):
    """Preserve fixed host expenditure and tenant food before assigning a crop."""
    if not 0 <= light_fraction <= 1 or not 0 <= hazard <= 1 or not 0 <= harvest_fraction <= 1:
        raise ValueError('Fractions must lie in [0,1]')
    gross = row['gpp_c'] * light_fraction
    fixed = (row['photolytic_gas_not_fixed_c'] + row['maintenance_c']
             + row['night_gas_substrate_c'] + row['renewal_biomass_c']
             + row['renewal_respiration_c'] + row['food_to_tenants_c'])
    replacement = hazard * row['construction_c_kg_m2']
    surplus = gross - fixed - replacement
    passes = surplus >= 0
    # Inherited construction yield and tissue carbon fraction, not food chemistry.
    retained_fraction = row['renewal_biomass_c'] / (row['renewal_biomass_c'] + row['renewal_respiration_c'])
    tissue_c_fraction = row['green_tissue_c'] / row['green_tissue_dry_kg']
    allocation = harvest_fraction * surplus if passes else None
    tissue_c = allocation * retained_fraction if passes else None
    respiration = allocation - tissue_c if passes else None
    remainder = surplus - allocation if passes else None
    dry = tissue_c / tissue_c_fraction if passes else None
    return {
        'class': row['cls'], 'diameter_m': row['diameter_m'],
        'light_fraction': light_fraction, 'annual_adult_hazard': hazard,
        'harvest_fraction_of_unallocated_carbon': harvest_fraction,
        'gross_c_kg_m2_year': gross, 'fixed_host_and_tenant_c_kg_m2_year': fixed,
        'population_replacement_c_kg_m2_year': replacement,
        'unallocated_before_harvest_c_kg_m2_year': surplus,
        'carbon_budget_nonnegative': passes,
        'product_assimilate_c_kg_m2_year': allocation,
        'product_retained_c_kg_m2_year': tissue_c,
        'product_construction_respiration_c_kg_m2_year': respiration,
        'remaining_for_other_uses_c_kg_m2_year': remainder,
        'allocation_t_dry_per_projected_ha_year': dry * 10 if passes else None,
        'carbon_closure_residual_kg_m2_year': gross - fixed - replacement - surplus,
        'existing_tenant_sugar_c_kg_m2_year': row['food_to_tenants_c'],
        'existing_green_renewal_t_dry_per_projected_ha_year': row['green_tissue_dry_kg'] * 10,
        'existing_structure_renewal_t_dry_per_projected_ha_year': row['structure_dry_kg'] * 10,
        'inherited_construction_carbon_retained_fraction': retained_fraction,
        'inherited_dry_tissue_carbon_fraction': tissue_c_fraction,
    }


def harvest_phosphorus(dry_t_ha, concentrations_g_kg, return_fraction):
    if dry_t_ha < 0 or not 0 <= return_fraction <= 1 or min(concentrations_g_kg) < 0:
        raise ValueError('Nonnegative mass and concentrations; return fraction in [0,1]')
    exported = [dry_t_ha * concentration for concentration in concentrations_g_kg]
    return {'exported_kg_p_ha_year': exported,
            'returned_upward_kg_p_ha_year': [x * return_fraction for x in exported],
            'unreturned_kg_p_ha_year': [x * (1 - return_fraction) for x in exported]}



def inhabited_clusters(sky, volume, settings):
    """Compare conditional tenant mass with existing whole-district mass budgets.

    This is not a floor-area or residential design. Domatia establishment and
    replacement, extra payload gas and local attachments remain uncosted.
    """
    dm = next(r for r in sky['tenants']['domatia']
              if r['diameter_m'] == settings['diameter_m']
              and r['food_to_tenants_c'] == settings['domatium_food_to_tenants_c'])
    green = next(r for r in sky['food_web']['production']
                 if r['cls'] == 'sky_reef' and r['diameter_m'] == settings['diameter_m']
                 and r['gas_route'] == 'photolysis' and r['gpp_c_kg_m2_year'] == 2)
    budget = allocate(green, 1, settings['green_adult_hazard'], 0)
    if not budget['carbon_budget_nonnegative'] or not dm['mass_gate'] or not dm['solar_gas']:
        raise ValueError('Reference cluster components failed their source or carbon gates')
    rows = []
    for n in settings['domatium_counts']:
        ngreen = n * settings['green_modules_per_domatium']
        available = ngreen * dm['area_m2'] * budget['unallocated_before_harvest_c_kg_m2_year'] / 1000
        support = n * dm['area_m2'] * max(0, -dm['carbon_balance_c_kg_m2_year']) / 1000
        payload = n * dm['steady_tenants_t']
        rows.append({'domatia': n, 'green_modules': ngreen,
                     'sum_projected_area_ha': (n + ngreen) * dm['area_m2'] / 1e4,
                     'conditional_domatia_tenant_allowance_t': payload,
                     'green_residual_after_green_replacement_t_c_year': available,
                     'domatia_routine_support_t_c_year': support,
                     'share_of_green_residual_for_routine_domatia_support': support / available,
                     'green_residual_after_routine_support_t_c_year': available - support,
                     'payload_equivalent_people_by_district_budget': {
                         name: payload / case['total_t_per_person']
                         for name, case in volume['assumptions']['districts'].items()},
                     'supported_floor_area_m2': None,
                     'qualification': 'Mass equivalents, not resident capacity. Projected biological area is not usable floor. District mass already includes generic buildings, floors and services; their actual support design and any excess mass remain open, together with extra gas, domatia replacement, wet load, trim and human use. No green-module payload allowance is counted.'})
    return rows


def result():
    products, settings = load()
    sky, aero = products['sky_ecology'], products['aerophytes']
    assumptions = settings['assumptions']
    phosphorus = sky['phosphorus']['design_guesses']
    rows = []
    for cls, diameter in assumptions['representatives_m'].items():
        body = next(r for r in sky['food_web']['production']
                    if r['cls'] == cls and r['diameter_m'] == diameter
                    and r['gas_route'] == 'photolysis' and r['gpp_c_kg_m2_year'] == 2)
        if not body['all_gates']:
            raise ValueError('Reference body failed its source gates')
        for light, hazard, fraction in itertools.product(assumptions['light_fractions'],
                                                        assumptions['annual_adult_hazards'],
                                                        assumptions['harvest_fractions']):
            row = allocate(body, light, hazard, fraction)
            dry = row['allocation_t_dry_per_projected_ha_year']
            row['leaflike_product_phosphorus'] = (harvest_phosphorus(dry, phosphorus['green_p_g_kg'],
                                                                   phosphorus['food_return'])
                                                   if dry is not None else None)
            rows.append(row)
    membrane = [r['film_m2_year'] for r in aero['resources']['membrane'] if r['body'] == 'colony_2km']
    return {
        'schema': 'terluna.provisioning.aerophyte-products/1',
        'evidence': 'Allocation arithmetic on hypothetical mature aerophytes. Generic new tissue is not a demonstrated edible, saleable or physically supportable crop; no market, farm, nutrient supply or harvest operation is simulated.',
        'reading_rule': 'Carbon flows are kg C per horizontal projected m2 per year; dry allocation is tonnes per projected hectare per year; crop phosphorus is kg P per projected hectare per year. Membrane is m2 per year for one specified 2 km colony. Inhabited clusters report total projected hectares, tonnes of conditional tenant allowance and annual tonnes C; people are payload equivalents under inherited district budgets, not capacity. Projected area is not leaf area, atmospheric volume or surface terrain. The reference uses inherited physiology at a gross fixation of 2 kg C/m2/year; light fractions are prescribed sensitivity cases, not angular transport or the final 301 K weather. Allocate only after fixed host costs, tenant sugar and adult replacement. Crop allocation competes with growth, sailing and extra tenants. New crop maintenance, night storage, decay, wet-load gas and harvest-wound repair are uncomputed deductions from this allocation ceiling. Phosphorus is an additional leaflike harvest loop, excluding existing host demand. Do not add crop potential to an already fully allocated biomass total.',
        'producer': {'files': {'provisioning/aerophytes.py': sha(Path(__file__))},
                     'inputs': {str(INPUTS.relative_to(ROOT)): sha(INPUTS),
                                str(SOURCES.relative_to(ROOT)): sha(SOURCES),
                                **{v['path']: v['sha256'] for v in settings['products'].values()}}},
        'assumptions': assumptions,
        'leaflike_product_p_g_per_kg_dry': phosphorus['green_p_g_kg'],
        'inherited_food_phosphorus_return_fraction': phosphorus['food_return'],
        'cases': rows,
        'inherited_2km_colony_membrane_m2_year': [min(membrane), max(membrane)],
        'inhabited_cluster_comparisons': inhabited_clusters(sky, products['inhabited_volume'], assumptions['inhabited_clusters']),
        'unresolved': ['Edibility, useful composition and processing recovery',
                       'New crop respiration, night storage, decay, wet-load gas and harvest-wound repair',
                       'Wet standing crop, buoyancy, trim and harvest schedules',
                       'Water, nitrogen and phosphorus delivery including new crop demand',
                       'Diffuse light through crowded communities and the final 301 K weather',
                       'Independent fixation by attached crops and their effect on the host',
                       'Access, transport, equipment, labor and comparative costs',
                       'Living-area clusters: domatia replacement and establishment, payload gas, connections, floors and human services',
                       'Competing ecological uses and sustainable wild harvest'],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'research/runs/provisioning/aerophytes.json')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result(), indent=2) + '\n')
    print(args.output)
