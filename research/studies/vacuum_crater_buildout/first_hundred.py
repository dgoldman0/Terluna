"""Conditional quantities, expenditure and production schedule for 100 crater towns.

The monetary and industrial coefficients are explicit planning assumptions.
The original catalogue selection and accepted community architecture are unchanged.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import heapq
import json
import math
from pathlib import Path

from shared.constants import STANDARD_ATMOSPHERE, STEFAN_BOLTZMANN, SYNODIC_MONTH_DAYS, JULIAN_YEAR_DAYS
from shared.provenance import constants_used
from research.studies.vacuum_crater_buildout.run import HERE, ROOT, enclosure, read_product, sha256

SCHEMA = 'terluna.research.first-hundred-crater-communities/1'


def load_inputs():
    spec = json.loads((HERE / 'first_hundred_inputs.json').read_text())
    products = {name: read_product(ROOT / item['path'], item['schema'], item['sha256'])
                for name, item in spec['inputs'].items()}
    if spec['communities'] != 100:
        raise ValueError('This study is scoped to the first 100 communities')
    low, high = products['plan']['initial_community_population']
    if any(not low <= p <= high for p in spec['population_cases']):
        raise ValueError('Population departs from the adopted initial band')
    q = spec['quantity']
    if not math.isclose(sum(q[k] for k in ('ordinary_local_bulk_fraction', 'ordinary_local_metal_fraction', 'ordinary_import_fraction')), 1):
        raise ValueError('Ordinary material fractions must sum to one')
    for key in ('bulk_process_yield', 'metal_product_per_raw_feed'):
        if not 0 < q[key] <= 1:
            raise ValueError(f'Invalid process yield: {key}')
    return spec, products


def production_rates(case):
    return dict(handling_kg=case['effective_handling_systems'] * case['handling_tonnes_hour_system'] * 1000 * case['handling_hours_year'],
                bulk_finished_kg=case['bulk_finished_kg_year'], metal_finished_kg=case['metal_finished_kg_year'],
                pressure_system_kg=case['pressure_system_kg_year'], special_import_kg=case['special_import_kg_year'],
                volatile_import_kg=case['volatile_import_kg_year'], building_floor_m2=case['building_floor_m2_year'],
                road_km=case['road_km_year'])


def quantities(people_per_town, spec, products, production_case='reference'):
    q, service = spec['quantity'], spec['services']
    n, people = spec['communities'], people_per_town * spec['communities']
    design = products['design']
    case = design['population_cases']['representative']
    lake_fraction = design['population_lake_to_settlement_diameter'] ** 2
    area = people * case['floor_m2_per_resident'] / ((1-lake_fraction) * case['non_lake_ground_building_share'] * case['mean_storeys'])
    diameter = math.sqrt(4 * area / n / math.pi)
    conf = deepcopy(products['vacuum_scenarios']['enclosure'])
    conf['material_effective_allowable_Pa'] = q['reference_allowable_Pa']
    one = enclosure(diameter, diameter * conf['reference_depth_over_diameter'],
                    conf['reference_pressure_atm'] * STANDARD_ATMOSPHERE, conf)
    roof = n * one['roof_area_m2']
    garden = area * (1-lake_fraction) * (1-case['non_lake_ground_building_share'])
    ordinary = people * products['economy']['cost_cases'][q['ordinary_stock_case']]['ordinary_stock_kg_per_lunar_resident']
    pressure_ideal = n * one['equivalent_pressure_material_kg']
    pressure_reinforcement = pressure_ideal * q['pressure_detail_mass_factor']
    barriers = roof * q['roof_barrier_kg_m2'] + area * q['ground_barrier_kg_m2']
    pressure = pressure_reinforcement + barriers
    support = roof * q['independent_support_kg_m2_roof']
    foundations = area * q['foundation_kg_m2_plan']
    shield = roof * q['shield_thickness_m'] * conf['shield_density_kg_m3']
    road = q['regional_road_km'] * 1000 * q['road_width_m'] * q['road_layer_m'] * q['road_density_kg_m3']
    op = products['vacuum_scenarios']['operation']
    gas = n * one['gross_gas_mass_kg'] * (1 + q['gas_reserve_fraction'])
    water = dict(lakes_kg=area * lake_fraction * q['lake_mean_depth_m'] * q['water_density_kg_m3'],
                 process_inventory_kg=people * op['water_loop_throughput_kg_per_person_day'] * q['water_inventory_days'],
                 garden_inventory_kg=garden * q['garden_water_kg_m2'])
    basic_W = people * service['basic_kW_resident'] * 1000
    food_W = people * service['food_kW_resident'] * 1000
    garden_peak_W = garden * service['garden_lamp_W_m2_on']
    garden_average_W = garden_peak_W * service['garden_light_hours_day'] / 24
    mean_W, peak_W = basic_W + food_W + garden_average_W, basic_W + food_W + garden_peak_W
    rates = production_rates(spec['production'][production_case])
    thermal_kWh = (rates['bulk_finished_kg'] * q['bulk_manufacturing_MJ_kg'] / 3.6
                   + rates['metal_finished_kg'] * q['metal_manufacturing_kWh_kg'])
    handling_kWh = rates['handling_kg'] * q['handling_MJ_kg'] / 3.6
    industry_W = 1000 * (thermal_kWh / service['industrial_operating_hours_year']
                        + handling_kWh / spec['production'][production_case]['handling_hours_year'])
    firm_W = peak_W + industry_W
    radiator = firm_W / (op['radiator_emissivity'] * STEFAN_BOLTZMANN * service['radiator_temperature_K']**4) * service['radiator_area_factor']
    storage_factor = (1 + op['storage_reserve_fraction']) / (op['storage_discharge_efficiency'] * op['storage_usable_fraction'])
    reserve_Wh = people * service['critical_kW_resident'] * 1000 * service['critical_reserve_hours'] * storage_factor
    night_Wh = mean_W * SYNODIC_MONTH_DAYS * 12 * storage_factor
    bulk = ordinary * q['ordinary_local_bulk_fraction'] + foundations + road
    metal = ordinary * q['ordinary_local_metal_fraction'] + support + radiator * q['radiator_local_metal_kg_m2']
    feed_bulk, feed_metal = bulk / q['bulk_process_yield'], metal / q['metal_product_per_raw_feed']
    feed = shield + feed_bulk + feed_metal
    grading = area * q['grading_mean_depth_m'] * conf['shield_density_kg_m3']
    special = dict(pressure_system_kg=pressure, ordinary_equipment_kg=ordinary*q['ordinary_import_fraction'],
                   construction_plant_kg=q['construction_plant_import_kg'], power_hardware_kg=firm_W*q['power_hardware_import_kg_W'],
                   reserve_hardware_kg=reserve_Wh*q['reserve_import_kg_Wh'],
                   garden_lamps_kg=garden*q['garden_lamp_import_kg_m2'], biosphere_kg=people*q['biosphere_import_kg_resident'],
                   transport_and_docking_kg=n*q['transport_and_docking_import_kg_community'],
                   radiator_controls_kg=radiator*q['radiator_controls_import_kg_m2'])
    establishment_import = special['construction_plant_kg'] + special['power_hardware_kg']
    volatile = gas + sum(water.values())
    demand = dict(handling_kg=feed+grading, bulk_finished_kg=bulk, metal_finished_kg=metal,
                  pressure_system_kg=pressure, special_import_kg=sum(special.values())-establishment_import, volatile_import_kg=volatile,
                  building_floor_m2=people * case['floor_m2_per_resident'], road_km=q['regional_road_km'])
    energy = dict(handling_kWh=(feed+grading)*q['handling_MJ_kg']/3.6,
                  bulk_kWh=bulk*q['bulk_manufacturing_MJ_kg']/3.6, metal_kWh=metal*q['metal_manufacturing_kWh_kg'])
    local_finished = bulk + metal + shield
    return dict(communities=n, population=people, residents_per_community=people_per_town,
                dome_diameter_m=diameter, reference_depth_m=one['depth_m'], total_plan_area_m2=area,
                total_roof_area_m2=roof, garden_area_m2=garden,
                estate_area_km2=area*products['vacuum_scenarios']['selection']['reference_estate_diameter_factor']**2/1e6,
                pressure_ideal_kg=pressure_ideal, pressure_reinforcement_kg=pressure_reinforcement, barriers_kg=barriers,
                pressure_system_kg=pressure, independent_support_kg=support, foundations_kg=foundations,
                shielding_kg=shield, ordinary_stock_kg=ordinary, road_material_kg=road,
                local_finished_kg=local_finished, raw_feed_kg=feed,
                retained_tails_kg=feed_bulk+feed_metal-bulk-metal,
                quarry_area_km2=feed/conf['shield_density_kg_m3']/q['quarry_working_depth_m']/1e6,
                road_corridor_area_km2=q['regional_road_km']*q['road_corridor_width_m']/1000,
                special_imports=special, establishment_import_kg=establishment_import, air_including_reserve_kg=gas, water_inventory=water,
                local_share_of_all_delivered_material=local_finished/(local_finished+sum(special.values())+volatile),
                future_raw_feed_reservation=[dict(usable_return_fraction=f, horizon_years=spec['reserve_diagnostic']['horizon_years'],
                    raw_feed_kg=(feed_bulk+feed_metal)*(1-f)*spec['reserve_diagnostic']['horizon_years']/spec['reserve_diagnostic']['mean_processed_asset_life_years'],
                    annual_raw_feed_kg=(feed_bulk+feed_metal)*(1-f)/spec['reserve_diagnostic']['mean_processed_asset_life_years'],
                    processed_asset_inventory_cycle_years=spec['reserve_diagnostic']['mean_processed_asset_life_years']/(1-f))
                    for f in spec['reserve_diagnostic']['usable_end_to_end_return_fractions']],
                services=dict(basic_W=basic_W, food_W=food_W, garden_average_W=garden_average_W, garden_peak_W=garden_peak_W,
                              mean_occupied_W=mean_W, peak_occupied_W=peak_W, industrial_peak_W=industry_W,
                              firm_capacity_W=firm_W, reserve_Wh=reserve_Wh, full_night_storage_Wh=night_Wh,
                              radiator_emitting_area_m2=radiator),
                production_demand=demand, construction_energy_kWh=energy,
                annual_makeup=dict(air_kg=n*one['gross_gas_mass_kg']*.01, air_inventory_cycle_years=100,
                                   processed_water_kg=people*op['water_loop_throughput_kg_per_person_day']*(1-op['water_recovery_fraction'])*JULIAN_YEAR_DAYS,
                                   boundary='Air: 1% of operating gross inventory; water: only the 100 kg/person/day processed loop at 98% recovery. Lake/crop losses, lock cycles and accidents require additional accounts.'))


def construction_cost(town, spec, products, price_name):
    p, q, s = spec['prices'][price_name], spec['quantity'], town['services']
    ordinary = town['population'] * products['economy']['cost_cases'][p['ordinary_case']]['ordinary_installed_capital_U_per_resident']
    direct = dict(ordinary_installed_capital=ordinary,
                  pressure_system=town['pressure_system_kg']*p['pressure_system_U_kg'],
                  independent_support_and_foundations=(town['independent_support_kg']+town['foundations_kg'])*p['incremental_local_structure_U_kg'],
                  local_shielding=town['shielding_kg']*p['shield_installed_U_kg'],
                  regional_roads=q['regional_road_km']*p['regional_road_U_km'],
                  pressurised_transport_and_docking=town['communities']*p['transport_access_U_community'],
                  dedicated_firm_power=s['firm_capacity_W']*p['firm_power_U_W'],
                  emergency_reserves=s['reserve_Wh']*p['reserve_U_Wh'],
                  exterior_radiators=s['radiator_emitting_area_m2']*p['radiator_U_m2'],
                  civic_garden_lighting=town['garden_area_m2']*p['garden_light_U_m2'],
                  initial_air_and_water=town['production_demand']['volatile_import_kg']*p['volatile_delivered_U_kg'],
                  ordinary_extra_delivery=town['special_imports']['ordinary_equipment_kg']*p['ordinary_import_freight_increment_U_kg'],
                  initial_biosphere_imports=town['special_imports']['biosphere_kg']*p['biosphere_delivered_U_kg'],
                  survey_and_preservation=p['survey_and_preservation_U'],
                  incremental_qualification=p['qualification_U'],
                  regional_industrial_establishment_increment=p['industrial_establishment_increment_U'],
                  construction_utility_operation=sum(town['construction_energy_kWh'].values())/1000*p['construction_utility_variable_U_MWh'],
                  disturbed_ground_restoration=town['raw_feed_kg']*p['restoration_U_kg_raw_feed'])
    subtotal = sum(direct.values())
    contingency = subtotal*p['contingency_fraction']
    removal = subtotal*p['removal_reserve_fraction']
    return dict(direct_U=direct, direct_total_U=subtotal, contingency_U=contingency,
                removal_reserve_U=removal, total_initial_funding_U=subtotal+contingency+removal,
                per_community_U=(subtotal+contingency+removal)/town['communities'],
                per_resident_U=(subtotal+contingency+removal)/town['population'],
                full_night_storage_only_U=s['full_night_storage_Wh']*p['reserve_U_Wh'])


def schedule(demand, communities, case, settings, *, ready_industry=False, establishment_import_kg=0):
    """Monthly conserved inventories plus finite site-team occupancy."""
    rates = production_rates(case)
    if communities < 1 or int(communities) != communities or case['site_teams'] < 1 or case['site_cycle_years'] <= 0:
        raise ValueError('Positive whole community count, teams and site cycle required')
    if set(demand) != set(rates) or any(not math.isfinite(v) or v <= 0 for v in (*demand.values(), *rates.values())):
        raise ValueError('Every production stream needs positive finite demand and capacity')
    step = 1/settings['steps_per_year']
    setup_delivery_years = establishment_import_kg / (rates['special_import_kg']*.25)
    lead = 0 if ready_industry else max(*case['lead_years'].values(),setup_delivery_years)
    inventories = dict.fromkeys(rates, 0.0)
    consumed = dict.fromkeys(rates, 0.0)
    made = dict.fromkeys(rates, 0.0)
    active, completions, starts = [], [], []
    for tick in range(1, round(settings['maximum_years']/step)+1):
        elapsed = tick*step
        year = int((tick-1)/settings['steps_per_year'])
        ramps = settings['ramp_fractions']
        ramp = 1 if ready_industry or year >= len(ramps) else ramps[year]
        for key, rate in rates.items():
            increment = min(rate*ramp*step, max(0.0, demand[key]-made[key]))
            inventories[key] += increment
            made[key] += increment
        while active and active[0] <= elapsed+1e-10:
            completions.append(heapq.heappop(active))
        while len(starts) < communities and len(active) < case['site_teams'] and all(inventories[k]+1e-10*demand[k]/communities >= demand[k]/communities for k in rates):
            starts.append(elapsed)
            for key in rates:
                amount = demand[key]/communities
                inventories[key] = max(0.0, inventories[key]-amount)
                consumed[key] += amount
            heapq.heappush(active, elapsed+case['site_cycle_years'])
        if len(starts) == communities:
            completions.extend(active)
            completions.sort()
            break
    else:
        raise ValueError('Schedule did not finish within the stated horizon')
    equivalent = {k:rates[k]/demand[k]*communities for k in rates}
    equivalent['site_teams'] = case['site_teams']/case['site_cycle_years']
    return dict(lead_years=lead, establishment_delivery_minimum_years=setup_delivery_years, production_and_completion_years=completions[-1], total_years=lead+completions[-1],
                occupation_years={str(k):lead+completions[k-1] for k in (1,10,25,50,100) if k <= communities},
                steady_communities_per_year=equivalent, slowest_steady_stream=min(equivalent,key=equivalent.get),
                last_start_year=lead+starts[-1], consumed=consumed, produced_before_last_start=made,
                remaining_inventory_at_last_start=inventories)


def run():
    spec, products = load_inputs()
    records = {}
    for residents in spec['population_cases']:
        physical = quantities(residents, spec, products)
        prices = {name:construction_cost(physical,spec,products,name) for name in spec['prices']}
        schedules = {}
        for name, production in spec['production'].items():
            own = quantities(residents,spec,products,name)
            schedules[name] = schedule(own['production_demand'],spec['communities'],production,spec['schedule'],establishment_import_kg=own['establishment_import_kg'])
            schedules[name]['power_capacity_W'] = own['services']['firm_capacity_W']
            schedules[name]['reference_price_funding_U'] = construction_cost(own,spec,products,'reference')['total_initial_funding_U']
        ready = schedule(physical['production_demand'],spec['communities'],spec['production']['reference'],spec['schedule'],ready_industry=True,establishment_import_kg=physical['establishment_import_kg'])
        records[str(residents)] = dict(quantities=physical,costs=prices,schedules=schedules,ready_industry_schedule=ready)
    variations = {}
    changes = [('reference',{}),('allowable_400MPa',{'quantity':{'reference_allowable_Pa':400e6}}),
               ('allowable_800MPa',{'quantity':{'reference_allowable_Pa':800e6}}),
               ('shield_7m',{'quantity':{'shield_thickness_m':7}}),
               ('food_53kW_person',{'services':{'food_kW_resident':53}}),
               ('lake_depth_1m',{'quantity':{'lake_mean_depth_m':1}}),
               ('roads_1000km',{'quantity':{'regional_road_km':1000}}),
               ('ordinary_stock_43t',{'quantity':{'ordinary_stock_case':'patient_automation'}}),
               ('ordinary_stock_200t',{'quantity':{'ordinary_stock_case':'expensive_energy_dilute_ore'}})]
    for name, overrides in changes:
        case = deepcopy(spec)
        for section, values in overrides.items():
            case[section].update(values)
        town = quantities(spec['reference_population_per_community'],case,products)
        variations[name] = dict(quantities=town,cost=construction_cost(town,case,products,'reference'),
                               schedule=schedule(town['production_demand'],spec['communities'],spec['production']['reference'],spec['schedule'],establishment_import_kg=town['establishment_import_kg']))
    food = products['food']['mixes']['fully_lamp_lit']['cases']['4']['kept_50']['power_tw']['total']
    files = ['research/studies/vacuum_crater_buildout/first_hundred.py','research/studies/vacuum_crater_buildout/run.py','research/studies/vacuum_crater_buildout/first_hundred_inputs.json']
    return dict(schema=SCHEMA, evidence=spec['evidence'], numeraire=spec['numeraire'], inputs=spec['inputs'],
                producer=dict(files={p:sha256(ROOT/p) for p in files},constants=constants_used(files)),
                food_process_bracket_kW_resident=[x/4 for x in food], cases=records, sensitivities=variations,
                construction_cost_boundary=spec['price_rules'],
                unresolved=['Site selection, heritage judgement, scientific sampling, geological reserves and future wild/dark-land corridors.',
                            'Qualified whole pressure boundary, anchors, independent shielding supports and transition/removal design.',
                            'Radiation transport and civic light quality; 3 m and 7 m remain mass comparisons.',
                            'Industrial yields, fleet productivity, supplier capacity, prices and qualification duration.',
                            'Regional roads, port/relay location, traffic capacity and a food/heat/water system with measured performance.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'research/runs/first_hundred.json')
    args = parser.parse_args()
    result = run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    for size, case in result['cases'].items():
        print(size,'residents/town:',round(case['costs']['reference']['total_initial_funding_U']/1e9,2),'billion U;',round(case['schedules']['reference']['total_years'],2),'years')


if __name__ == '__main__':
    main()
