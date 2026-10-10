"""Conditional sky-reef symbiosis requirements, not a population or defence model.

All area densities are per horizontal reef FOOTPRINT unless explicitly named
exterior. Carbon is kg C; phosphorus is kg P; biomass is kg live mass. Annual
replacement is a population turnover budget, not an individual's growth rate.
Run as a module; default output is ignored research/runs/sky_ecology/symbioses.json.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

from research.studies.sky_ecology import stoichiometry as st

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent.parent
SCHEMA = 'terluna.research.sky-reef-symbioses/1'


def _nonnegative(**values):
    for name, value in values.items():
        if not math.isfinite(value) or value < 0:
            raise ValueError(f'{name} must be finite and nonnegative')


def _fraction(name, value, positive=False):
    if not math.isfinite(value) or not 0 <= value <= 1 or (positive and value == 0):
        raise ValueError(f'{name} must be in {"(0, 1]" if positive else "[0, 1]"}')


def population(live_density, individual_live_kg, footprint_m2, dry_fraction,
               active_fraction, exterior_per_footprint):
    """Collective biomass does not specify individual size. Counts can be fractional
    expectations over area, not evidence of enough individuals at every point."""
    _nonnegative(live_density=live_density, individual_live_kg=individual_live_kg,
                 footprint_m2=footprint_m2, exterior_per_footprint=exterior_per_footprint)
    if min(individual_live_kg, footprint_m2, exterior_per_footprint) <= 0:
        raise ValueError('individual mass, footprint and area ratio must be positive')
    _fraction('dry_fraction', dry_fraction, positive=True)
    _fraction('active_fraction', active_fraction)
    density = live_density / individual_live_kg
    return dict(live_kg_per_footprint_m2=live_density, individual_live_kg=individual_live_kg,
                live_kg=live_density*footprint_m2, dry_kg=live_density*footprint_m2*dry_fraction,
                individuals=density*footprint_m2, active_individuals=density*footprint_m2*active_fraction,
                individuals_per_footprint_m2=density,
                active_per_footprint_m2=density*active_fraction,
                individuals_per_exterior_m2=density/exterior_per_footprint,
                active_per_exterior_m2=density*active_fraction/exterior_per_footprint)


def response_seconds(spacing_m, speed_m_s, activation_delay_s, path_factor):
    """Worst point in a square nest grid lies spacing/sqrt(2) from nearest nest.
    Delay plus prescribed travel time only: arrival says nothing about efficacy."""
    _nonnegative(spacing=spacing_m, speed=speed_m_s, delay=activation_delay_s, path=path_factor)
    if speed_m_s <= 0 or path_factor < 1:
        raise ValueError('positive speed and path factor at least one required')
    return activation_delay_s + path_factor*spacing_m/math.sqrt(2)/speed_m_s


def spacing_for_deadline(deadline_s, speed_m_s, activation_delay_s, path_factor):
    """Largest allowed spacing under this travel rule; None means no positive
    spacing can meet the deadline, even before fighting or recruitment."""
    response_seconds(0., speed_m_s, activation_delay_s, path_factor)
    _nonnegative(deadline=deadline_s)
    if deadline_s <= activation_delay_s:
        return None
    return (deadline_s-activation_delay_s)*speed_m_s*math.sqrt(2)/path_factor


def annual_requirements(live_density, dry_fraction, turnover_per_year, carbon_per_dry,
                        phosphorus_per_dry, carbon_assimilation,
                        population_production_efficiency, phosphorus_assimilation):
    """P/B replacement at prescribed turnover. Production efficiency is new animal
    C / assimilated C at the POPULATION level, already including respiratory costs.
    It cannot be reused for additional metabolism or assumed valid for every taxon.
    Holding that efficiency fixed while changing turnover also changes implied
    maintenance. Near-zero turnover requires independent metabolic calibration;
    this screen cannot diagnose an individual's feeding or maintenance requirement.
    Positive turnover is required: zero turnover here would spuriously erase upkeep.
    P requirement is gross replacement incorporation, not automatically fresh import;
    any internal recycling must be supplied explicitly without counting it as new P.
    """
    _nonnegative(live_density=live_density, turnover=turnover_per_year)
    if turnover_per_year == 0:
        raise ValueError('positive turnover required for this replacement-budget screen')
    for name, value in [('dry', dry_fraction), ('C/dry', carbon_per_dry),
                        ('P/dry', phosphorus_per_dry), ('C assimilation', carbon_assimilation),
                        ('production efficiency', population_production_efficiency),
                        ('P assimilation', phosphorus_assimilation)]:
        _fraction(name, value, positive=True)
    if carbon_per_dry + phosphorus_per_dry > 1:
        raise ValueError('C and P cannot exceed total dry tissue')
    replacement_dry = live_density*dry_fraction*turnover_per_year
    new_c, new_p = replacement_dry*carbon_per_dry, replacement_dry*phosphorus_per_dry
    food_c = new_c/carbon_assimilation/population_production_efficiency
    return dict(live_kg_per_footprint_m2=live_density, turnover_per_year=turnover_per_year,
                replacement_dry_kg_m2_year=replacement_dry,
                replacement_c_kg_m2_year=new_c, replacement_p_kg_m2_year=new_p,
                food_c_required_kg_m2_year=food_c,
                unassimilated_c_kg_m2_year=food_c*(1-carbon_assimilation),
                respiratory_c_kg_m2_year=food_c*carbon_assimilation-new_c,
                dietary_p_required_kg_m2_year=new_p/phosphorus_assimilation)


def allocate_food(guilds, host_food_c, traits):
    """One host sugar budget shared by all guilds. External dietary supplies are
    exclusive allocations crossing the boundary, never unremoved resident prey.
    Report deficits, without converting a budget closure into viable populations.
    No extra respiration charge; no P credit for sugar or for unmodelled recycling.
    """
    _nonnegative(host_food_c=host_food_c)
    shares = [g['host_budget_share'] for g in guilds]
    for share in shares:
        _fraction('host_budget_share', share)
    if sum(shares) > 1 + 1e-12:
        raise ValueError('guilds over-allocate the shared host food budget')
    rows = []
    for g in guilds:
        ext_c = g['external_food_c_supply_kg_m2_year']
        diet_p = g['external_diet_p_supply_kg_m2_year']
        _nonnegative(external_c=ext_c, dietary_p=diet_p)
        if ext_c == 0 and diet_p > 0:
            raise ValueError('external food P requires paired external food C')
        r = annual_requirements(g['live_kg_per_footprint_m2'],
                                turnover_per_year=g['turnover_per_year'], **traits)
        need_c = r['food_c_required_kg_m2_year']
        allocation = host_food_c*g['host_budget_share']
        external_used = min(ext_c, need_c)
        host_used = min(allocation, need_c-external_used)
        consumed_p = diet_p*external_used/ext_c if ext_c else 0.
        usable_p = consumed_p*traits['phosphorus_assimilation']
        external_c_for_p = (r['replacement_p_kg_m2_year']*ext_c/
                           (diet_p*traits['phosphorus_assimilation'])) if diet_p else None
        r.update(name=g['name'], host_c_allocation_kg_m2_year=allocation,
                 host_c_used_kg_m2_year=host_used,
                 host_c_unused_allocation_kg_m2_year=allocation-host_used,
                 external_c_supply_kg_m2_year=ext_c,
                 external_c_used_kg_m2_year=external_used,
                 external_c_unused_supply_kg_m2_year=ext_c-external_used,
                 external_c_required_for_p_kg_m2_year=external_c_for_p,
                 food_c_shortfall_kg_m2_year=max(0., need_c-host_used-external_used),
                 external_diet_p_supply_kg_m2_year=diet_p,
                 external_diet_p_consumed_kg_m2_year=consumed_p,
                 external_diet_p_unconsumed_kg_m2_year=diet_p-consumed_p,
                 assimilable_p_supply_kg_m2_year=usable_p,
                 unassimilated_consumed_p_kg_m2_year=consumed_p-usable_p,
                 incorporation_p_shortfall_kg_m2_year=max(0., r['replacement_p_kg_m2_year']-usable_p))
        rows.append(r)
    return dict(guilds=rows, host_food_c_kg_m2_year=host_food_c,
                host_c_allocated_kg_m2_year=host_food_c*sum(shares),
                host_c_unallocated_kg_m2_year=host_food_c*(1-sum(shares)),
                host_c_used_kg_m2_year=sum(r['host_c_used_kg_m2_year'] for r in rows),
                external_c_used_kg_m2_year=sum(r['external_c_used_kg_m2_year'] for r in rows))


def nutrient_import(gross_external_food_p_kg_m2_year, deposition_share_on_reef,
                    retention_share, initially_bioavailable_share, exports_from_reef_p_kg_m2_year):
    """Boundary P account. Recycling is not new import; retained P is not necessarily
    immediately available. Exports exclude deposition losses already subtracted here.
    Negative net import is allowed. Do not add this account to guild food supplies.
    """
    _nonnegative(gross_external_food_p=gross_external_food_p_kg_m2_year,
                 exports=exports_from_reef_p_kg_m2_year)
    for name, value in [('deposition', deposition_share_on_reef), ('retention', retention_share),
                        ('availability', initially_bioavailable_share)]:
        _fraction(name, value)
    deposited = gross_external_food_p_kg_m2_year*deposition_share_on_reef
    retained = deposited*retention_share
    return dict(gross_external_food_p_kg_m2_year=gross_external_food_p_kg_m2_year,
                deposited_p_kg_m2_year=deposited, deposition_loss_p_kg_m2_year=deposited-retained,
                retained_import_p_kg_m2_year=retained,
                initially_bioavailable_import_p_kg_m2_year=retained*initially_bioavailable_share,
                exports_p_kg_m2_year=exports_from_reef_p_kg_m2_year,
                net_total_p_import_kg_m2_year=retained-exports_from_reef_p_kg_m2_year)


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def results(root=ROOT, scenarios_path=HERE/'symbioses_scenarios.json',
            sources_path=HERE/'symbioses_sources.json'):
    """Small algebraic screen; source-register changes invalidate the product hash."""
    scenario = json.loads(scenarios_path.read_text())
    if scenario['schema'] != 'terluna.research.sky-reef-symbioses-scenarios/1':
        raise ValueError('unexpected symbiosis scenario schema')
    host_budget = scenario['host_food_c_kg_per_footprint_m2_year']
    _nonnegative(host_budget=host_budget)
    if host_budget == 0:
        raise ValueError('positive host budget required for displayed budget fractions')
    product_path = root/'research/studies/sky_ecology/results/sky_ecology.json'
    product = json.loads(product_path.read_text())
    if product['schema'] != 'terluna.research.sky-ecology/1':
        raise ValueError('unexpected sky ecology schema')
    reef = next(r for r in product['tenants']['real_estate']['rows'] if r['body'] == 'raft_module')
    surface_ratio = sum(reef['zone_m2_per_m2'].values())
    traits = dict(dry_fraction=scenario['dry_fraction'], carbon_per_dry=st.CONSUMER_C_PER_DRY,
                  phosphorus_per_dry=sum(st.TENANT_P_G_KG)/len(st.TENANT_P_G_KG)/1000,
                  phosphorus_assimilation=sum(st.CONSUMER_P_ASSIMILATION)/len(st.CONSUMER_P_ASSIMILATION),
                  carbon_assimilation=scenario['carbon_assimilation'],
                  population_production_efficiency=scenario['population_production_efficiency'])
    densities = scenario['density_live_kg_per_footprint_m2']
    abundance = [population(b, m, scenario['footprint_m2'], scenario['dry_fraction'],
                            scenario['active_fraction'], surface_ratio)
                 for b, m in itertools.product(densities, scenario['individual_live_kg'])]
    demand = []
    for b, t in itertools.product(densities, scenario['turnover_per_year']):
        r = annual_requirements(b, turnover_per_year=t, **traits)
        r['fraction_of_whole_host_food_c_budget'] = r['food_c_required_kg_m2_year']/scenario['host_food_c_kg_per_footprint_m2_year']
        r['sugar_only_p_incorporation_shortfall_kg_m2_year'] = r['replacement_p_kg_m2_year']
        demand.append(r)
    response = scenario['response']
    arrivals = [dict(spacing_m=s, speed_m_s=v, response_s=response_seconds(s, v,
                    response['activation_delay_s'], response['path_factor']))
                for s, v in itertools.product(response['nest_spacing_m'], response['speed_m_s'])]
    spacing = [dict(deadline_s=d, speed_m_s=v, maximum_spacing_m=spacing_for_deadline(d, v,
                   response['activation_delay_s'], response['path_factor']))
               for d, v in itertools.product(response['deadline_s'], response['speed_m_s'])]
    paths = {'symbioses_scenarios.json': scenarios_path, 'symbioses_sources.json': sources_path,
             'results/sky_ecology.json': product_path}
    return dict(schema=SCHEMA,
                producer=dict(model='research/studies/sky_ecology/symbioses.py',
                              model_sha256=_digest(Path(__file__)),
                              dependencies={'stoichiometry.py': _digest(Path(st.__file__))},
                              inputs={name: _digest(path) for name, path in paths.items()}),
                evidence='Algebraic requirements for illustrative scenarios; no defended population, coexistence, persistence or empirically calibrated rates.',
                reading_rule='All m2 denominators are reef footprint except explicitly named exterior densities. Live and dry mass differ. Response is only delay plus travel. Food allocations share ONE host budget. Paired external food is used first, then host sugar; P is credited only in the consumed external food. Fixed production efficiency with varying turnover changes implied maintenance, not independently calibrated metabolism; no individual feeding inference or near-zero-turnover maintenance claim. P incorporation is not fresh P import. Nutrient boundary account is separate; do not add it to food supplies. No efficacy or adequacy flag.',
                scenarios=scenario, traits=traits, exterior_per_footprint=surface_ratio,
                abundance=abundance, annual_requirements=demand, response=arrivals, deadline_spacing=spacing,
                mixed_guild_budget=allocate_food(scenario['mixed_guild_example'],
                    scenario['host_food_c_kg_per_footprint_m2_year'], traits),
                nutrient_import=nutrient_import(**{k: v for k, v in scenario['nutrient_import_example'].items() if k != 'boundary'}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'research/runs/sky_ecology/symbioses.json')
    args = parser.parse_args()
    data = results()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
