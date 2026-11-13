"""Conditional process-gate cost components, with unpriced obligations explicit."""
from __future__ import annotations
import math
from .lunar import weighted_quantiles,linear_estimate

def crf(rate,life):
    if not math.isfinite(rate) or rate<0 or life<=0:raise ValueError('Invalid capital recovery inputs')
    return 1/life if rate==0 else rate/(1-(1+rate)**(-life))

def unit_feed_cost(coeff,case,rate):
    l,s,p,e=(coeff[k] for k in ('lunar','source','prices','economics'))
    year_hours=coeff['constants']['JULIAN_YEAR_DAYS']*coeff['constants']['JULIAN_DAY']/3600
    dollars=coeff['normalization']['USD2026_to_USD2041'];scale=case['equipment_cost_multiplier']
    upkeep=e['supplier_upkeep_fraction_year']; margin=e['capacity_margin']
    def equipment(mass,price,life,annual_capacity):
        return mass*price*dollars*scale*margin*(crf(rate,life)+upkeep)/annual_capacity
    availability=min(1,l['mine_availability']*case['availability_multiplier'])
    lunar_exc=equipment(l['excavator_kg'],p['lunar_hardware_USD2026_kg'],l['robot_life_years'],l['excavator_kg_hr']*year_hours*availability)
    lunar_haul=equipment(l['hauler_kg'],p['lunar_hardware_USD2026_kg'],l['robot_life_years'],l['hauler_payload_kg']/l['haul_roundtrip_hr']*year_hours*availability)
    lunar_energy=l['handling_kWh_kg']*case['electricity_USD2041_kWh']
    lunar_restore=e['local_restoration_USD2026_kg_feed']*dollars
    lunar_handling=lunar_exc+lunar_haul+lunar_energy+lunar_restore
    avail=min(1,l['process_availability']*case['availability_multiplier'])
    bulk_plant=equipment(l['bulk_train_hardware_kg'],p['lunar_hardware_USD2026_kg'],l['plant_life_years'],l['bulk_train_kg_hr']*year_hours*avail)
    bulk_energy=l['bulk_kWh_kg']*case['electricity_USD2041_kWh']
    source_avail=min(1,s['mine_availability']*case['availability_multiplier'])
    source_feed=s['mine_module_kg_hr']*year_hours*source_avail
    source_mine=equipment(s['mine_module_hardware_kg'],p['source_mine_USD2026_kg'],coeff['source_mine_life_years'],source_feed)
    energy=s['ore_kWh_kg']*case['ore_energy_multiplier']
    # An explicitly provisional chemistry-capital proxy per active thermal W.
    # Feed energy is not also added to inherited aggregate metal processing.
    thermal_W=s['mine_module_kg_hr']*energy*1000
    processor_capital=thermal_W*s['chemistry_USD2026_W']*dollars*scale*margin*(crf(rate,s['chemistry_life_years'])+upkeep)/source_feed
    source_energy=energy*case['electricity_USD2041_kWh']
    transition=coeff['source_transition_USD2026_per_kg_feed']*dollars
    feed=source_mine+processor_capital+source_energy+transition
    yield_water=case['water_grade']*case['water_recovery']
    if yield_water<=0 or yield_water>1:raise ValueError('Invalid recovered water yield')
    return {'lunar_handling_USD2041_kg_feed':lunar_handling,'lunar_handling_components':{'excavation':lunar_exc,'haul':lunar_haul,'energy':lunar_energy,'restoration_allowance':lunar_restore},
        'lunar_bulk_USD2041_kg_product':lunar_handling/l['bulk_yield']+bulk_plant+bulk_energy,
        'lunar_bulk_components':{'handling':lunar_handling/l['bulk_yield'],'plant':bulk_plant,'energy':bulk_energy},
        'source_USD2041_kg_feed':feed,'source_feed_components':{'mine':source_mine,'processor_proxy':processor_capital,'energy':source_energy,'transition_allowance':transition},
        'source_water_USD2041_kg':feed/yield_water,'source_water_kWh_kg':energy/yield_water,'source_water_feed_kg_kg':1/yield_water}

def cost_product(coeff,config,pixels):
    trials={name:unit_feed_cost(coeff,v,config['real_capital_rate']) for name,v in config['cases'].items()}
    weights=[p['solid_angle_sr'] for p in pixels]; total=sum(weights); oxide={}
    for key in ('FeO','Al2O3','TiO2'):
        grades=[p['abundances'][key] for p in pixels]
        gq=weighted_quantiles(grades,weights,(.05,.5,.95))
        # Quantiles of grades carry a spatial range, not a distribution of future mines.
        burden=[None if g<=0 else 1/(g*.7) for g in reversed(gq)]
        formal_nonpositive=sum(w for p,w in zip(pixels,weights) if linear_estimate(p,{key:1},1.96)['formal_interval'][0]<=0)/total
        oxide[key]={'grade_spatial_p05_p50_p95':gq,'feed_kg_per_kg_oxide_equivalent_at_70percent_recovery_high_median_low_grade':burden,'handling_only_USD2041_kg_equivalent_reference':[None if b is None else b*trials['reference']['lunar_handling_USD2041_kg_feed'] for b in burden], 'grade_recovery_handling_scenario_envelope':[trials['favorable']['lunar_handling_USD2041_kg_feed']/(gq[2]*max(config['oxide_equivalent_recovery_trials'])),None if gq[0]<=0 else trials['adverse']['lunar_handling_USD2041_kg_feed']/(gq[0]*min(config['oxide_equivalent_recovery_trials']))], 'zero_grade_surface_fraction':sum(w for g,w in zip(grades,weights) if g==0)/total,'formal_1_96sigma_lower_nonpositive_surface_fraction':formal_nonpositive,'full_cost_upper_bound':None,'interpretation':'Oxide-equivalent composition accounting only. Not separated oxides or metals; reduction, beneficiation, finishing and their yield/cost remain unpriced. Formal nonpositive tails have no finite inferred cost upper bound.'}
    # Pure yield sensitivity holds the process costs fixed, and does not redo ferry feedback.
    grade_yield=[{'grade':g,'recovery':r,'USD2041_kg_water':trials['reference']['source_USD2041_kg_feed']/(g*r),'feed_vs_old_10percent_grade_70percent_recovery':.07/(g*r)} for g in (.03,.0684,.1) for r in (.5,.7,.8)]
    return {'cases':trials,'water_grade_recovery_sensitivity':grade_yield,'lunar_oxide_feed_burden':oxide,'unpriced':['source prospecting and sampling programme','native-condition archive construction and perpetual care','site-specific exclusion buffers and future-accessible stocks','delivered plant freight and site construction','labour and common supporting infrastructure','water purification and product certification','dust and exhaust capture specific hardware','waste/residue containment and safe delivery qualification','source-to-depot transport and propellant feedback','beneficiation, reduction and finishing of named metal products'], 'rule':'Partial process-gate screening costs. Each product is considered separately; common feed costs are not additive coproduct charges. Frozen fully purchased-equipment valuations are not external cash for a self-fabricating system. No omitted item has an implicit zero cost.'}
