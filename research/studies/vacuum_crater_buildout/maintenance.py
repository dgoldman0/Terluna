"""Annual physical upkeep, variable services and renewal of one crater town."""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path

from shared.constants import JULIAN_YEAR_DAYS
from shared.provenance import constants_used
from research.studies.vacuum_crater_buildout.run import HERE, ROOT, read_product, sha256

SCHEMA='terluna.research.crater-town-maintenance/1'


def load_inputs():
    spec=json.loads((HERE/'maintenance_inputs.json').read_text())
    source={k:read_product(ROOT/v['path'],v['schema'],v['sha256']) for k,v in spec['inputs'].items()}
    ref=spec['regimes']['reference']['life_years']
    annual=source['annual_economy']['incremental_asset_lives_years']
    if ref['firm_power'] != annual['orbital_power_and_beaming'] or ref['emergency_storage'] != annual['critical_energy_reserves']:
        raise ValueError('Reference utility lives differ from the stated economic transfer')
    return spec,source


def assets(town, construction, price_name, spec):
    q=town['quantities']; p=construction['prices'][price_name]
    direct=town['costs'][price_name]['direct_U']; n=q['communities']
    ordinary=direct['ordinary_installed_capital']/n
    share=spec['ordinary_durable_capital_fraction']
    operation_share=q['services']['peak_occupied_W']/q['services']['firm_capacity_W']
    return dict(ordinary_durable=ordinary*share,ordinary_equipment=ordinary*(1-share),
                pressure_reinforcement=q['pressure_reinforcement_kg']*p['pressure_system_U_kg']/n,
                pressure_barriers=q['barriers_kg']*p['pressure_system_U_kg']/n,
                independent_support=q['independent_support_kg']*p['incremental_local_structure_U_kg']/n,
                foundations=q['foundations_kg']*p['incremental_local_structure_U_kg']/n,
                roads=direct['regional_roads']/n,transport_docking=direct['pressurised_transport_and_docking']/n,
                firm_power=direct['dedicated_firm_power']*operation_share/n,
                emergency_storage=direct['emergency_reserves']/n,
                radiators=direct['exterior_radiators']*operation_share/n,
                garden_lighting=direct['civic_garden_lighting']/n)


def replacement_count(horizon,life):
    if not math.isfinite(horizon) or not math.isfinite(life) or horizon <= 0 or life <= 0:
        raise ValueError('Positive finite horizon and life required')
    return max(0,math.ceil(horizon/life)-1)


def account(town, construction, spec, price_name='reference',regime_name='reference'):
    q=town['quantities']; n=q['communities']; people=q['population']/n
    regime=spec['regimes'][regime_name]; prices=construction['prices'][price_name]
    capital=assets(town,construction,price_name,spec)
    life=regime['life_years']
    if set(capital) != set(life) or any(v <= 0 for v in life.values()):
        raise ValueError('Every maintained asset needs a positive service life')
    if not 0 <= regime['processed_water_recovery_fraction'] <= 1:
        raise ValueError('Recovery must lie in [0,1]')
    upkeep={k:v*regime['routine_upkeep_fraction_year'] for k,v in capital.items()}
    replacement_value=dict(capital)
    # Renewed ordinary imported equipment again incurs its extra vacuum freight.
    replacement_value['ordinary_equipment'] += town['costs'][price_name]['direct_U']['ordinary_extra_delivery']/n
    operation_share=q['services']['peak_occupied_W']/q['services']['firm_capacity_W']
    recipe=construction['quantity']
    local_masses=dict(independent_support=q['independent_support_kg']/n,
                      foundations=q['foundations_kg']/n,roads=q['road_material_kg']/n,
                      radiators=q['services']['radiator_emitting_area_m2']*operation_share*recipe['radiator_local_metal_kg_m2']/n)
    replacement_utility={}
    for key,mass in local_masses.items():
        metal=key in ('independent_support','radiators')
        energy_kWh=mass*(recipe['metal_manufacturing_kWh_kg'] if metal else recipe['bulk_manufacturing_MJ_kg']/3.6)
        # Conservative raw-equivalent handling; no unproved recycling energy credit.
        yield_fraction=recipe['metal_product_per_raw_feed'] if metal else recipe['bulk_process_yield']
        energy_kWh += mass/yield_fraction*recipe['handling_MJ_kg']/3.6
        replacement_utility[key]=energy_kWh/1000*spec['replacement_utility_U_MWh'][price_name]
        replacement_value[key] += replacement_utility[key]
    renewal={k:v/life[k] for k,v in replacement_value.items()}
    annual_energy_MWh=q['services']['mean_occupied_W']/n/1e6*JULIAN_YEAR_DAYS*24
    operating_air=q['air_including_reserve_kg']/n/(1+construction['quantity']['gas_reserve_fraction'])
    air=operating_air*regime['annual_air_loss_fraction']
    # Read the independently pinned water-loop throughput before applying
    # this scenario's recovery.
    reference_recovery=construction['inputs']['vacuum_scenarios']
    vacuum=read_product(ROOT/reference_recovery['path'],reference_recovery['schema'],reference_recovery['sha256'])
    throughput=vacuum['operation']['water_loop_throughput_kg_per_person_day']
    water=people*throughput*(1-regime['processed_water_recovery_fraction'])*JULIAN_YEAR_DAYS
    extra_water=people*regime['additional_water_makeup_kg_person_day']*JULIAN_YEAR_DAYS
    biosphere=q['special_imports']['biosphere_kg']/n*regime['annual_biosphere_import_fraction']
    shield_rehandled=q['shielding_kg']/n*spec['shield_rehandling_fraction_year']
    variable=dict(utility_operation=annual_energy_MWh*prices['construction_utility_variable_U_MWh'],
                  air_makeup=air*prices['volatile_delivered_U_kg'],
                  processed_water_makeup=water*prices['volatile_delivered_U_kg'],
                  additional_water_makeup=extra_water*prices['volatile_delivered_U_kg'],
                  biosphere_makeup=biosphere*prices['biosphere_delivered_U_kg'],
                  retained_shield_rehandling=shield_rehandled*prices['shield_installed_U_kg'],
                  shield_rehandling_utility=shield_rehandled*recipe['handling_MJ_kg']/3.6/1000*spec['replacement_utility_U_MWh'][price_name])
    routine=sum(upkeep.values()); services=sum(variable.values()); fund=sum(renewal.values())
    margin=(routine+services)*spec['routine_operating_margin_fraction']
    total=routine+services+margin+fund
    horizons={}
    for years in spec['horizons_years']:
        counts={k:replacement_count(years,life[k]) for k in capital}
        replacement={k:replacement_value[k]*count for k,count in counts.items()}
        recurring=(routine+services+margin)*years
        horizons[str(years)]=dict(component_replacements=counts,full_replacement_U=replacement,
                                  recurring_upkeep_services_margin_U=recurring,
                                  total_vacuum_expenditure_U=recurring+sum(replacement.values()),
                                  boundary='Excludes initial build and final removal. Full replacements are strictly before endpoint; annual reserve contributions are not also charged.')
    return dict(population=people,price_case=price_name,maintenance_regime=regime_name,
                initial_build_share_U=town['costs'][price_name]['per_community_U'],
                maintained_installed_value_U=capital,full_replacement_value_U=replacement_value,
                local_replacement_utility_U_per_event=replacement_utility,asset_lives_years=life,
                routine_upkeep_U_year=upkeep,variable_services_U_year=variable,renewal_reserve_U_year=renewal,
                totals_U_year=dict(routine_upkeep=routine,variable_services=services,operating_margin=margin,
                                   recurring_cash_before_renewal=routine+services+margin,renewal_fund=fund,
                                   total_annual_funding=total),
                annual_funding_U_resident=total/people,annual_funding_fraction_of_initial_build=total/town['costs'][price_name]['per_community_U'],
                operating_quantities=dict(electricity_MWh_year=annual_energy_MWh,
                     air_makeup_kg_year=air,air_inventory_cycle_years=1/regime['annual_air_loss_fraction'] if regime['annual_air_loss_fraction'] else None,
                     processed_water_makeup_kg_year=water,additional_water_makeup_kg_year=extra_water,
                     biosphere_makeup_kg_year=biosphere,retained_shield_rehandled_kg_year=shield_rehandled,
                     pressure_replacement_kg_year=(q['pressure_reinforcement_kg']/life['pressure_reinforcement']+q['barriers_kg']/life['pressure_barriers'])/n),
                regional_100_town_requirements=dict(total_annual_funding_U=total*n,
                     steady_air_and_water_import_kg_year=(air+water+extra_water)*n,
                     construction_dedicated_volatile_import_kg_year=construction['production']['reference']['volatile_import_kg_year'],
                     rule='Operating deliveries add to the construction-only initial-fill capacity; assume separately contracted capacity rather than borrowing it without changing the build schedule.'),
                horizons=horizons)


def run():
    spec,source=load_inputs()
    records={}
    for population,town in source['construction_result']['cases'].items():
        records[population]={price:{regime:account(town,source['construction_inputs'],spec,price,regime)
                                  for regime in spec['regimes']} for price in source['construction_inputs']['prices']}
    files=['research/studies/vacuum_crater_buildout/maintenance.py','research/studies/vacuum_crater_buildout/maintenance_inputs.json','research/studies/vacuum_crater_buildout/run.py']
    return dict(schema=SCHEMA,evidence=spec['evidence'],numeraire=source['construction_result']['numeraire'],
                inputs=spec['inputs'],producer=dict(files={p:sha256(ROOT/p) for p in files},constants=constants_used(files)),
                rules=spec['rules'],cases=records)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'research/runs/crater_town_maintenance.json')
    args=parser.parse_args(); result=run()
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(result,indent=2)+'\n')
    for size,cases in result['cases'].items():
        c=cases['reference']['reference']
        print(size,'residents:',{k:round(v/1e6,3) for k,v in c['totals_U_year'].items()},'million U/year')


if __name__=='__main__':
    main()
