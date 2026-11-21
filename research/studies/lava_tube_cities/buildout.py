"""Conditional construction and sustaining budget for four cities and 100 towns.

Products are read through schema/hash bindings. No catalogue lead is selected or
promoted to a measured cavern. Every added design/price input lives in the JSON.
"""
from __future__ import annotations
from copy import deepcopy
import argparse
import hashlib
import json
import math
from pathlib import Path

from shared.constants import (GAS_CONSTANT, JULIAN_DAY, JULIAN_YEAR_DAYS,
                              MOON_SURFACE_GRAVITY, STANDARD_ATMOSPHERE,
                              STEFAN_BOLTZMANN, SYNODIC_MONTH_DAYS)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.lava-tube-buildout/1'
YEAR_HOURS = JULIAN_YEAR_DAYS * JULIAN_DAY / 3600


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_inputs(spec=None):
    spec = spec if spec is not None else json.loads((HERE / 'buildout_inputs.json').read_text())
    source = {}
    for key, binding in spec['inputs'].items():
        path = ROOT / binding['path']
        if sha(path) != binding['sha256']:
            raise ValueError(f'Input hash mismatch: {key}')
        value = json.loads(path.read_text())
        if 'schema' in binding and value.get('schema') != binding['schema']:
            raise ValueError(f'Input schema mismatch: {key}')
        source[key] = value
    return spec, source


def food_system(s, source, population):
    f = s['food']; basis = source['food']
    fish_edible = population * f['fish_edible_g_person_day'] / 1000 * JULIAN_YEAR_DAYS
    fish_live = fish_edible / f['fillet_yield']
    feed = fish_live * f['feed_conversion_ratio']
    microbial_protein = population * f['microbial_protein_g_person_day'] / 1000 * JULIAN_YEAR_DAYS
    microbial_biomass = microbial_protein / f['microbial_protein_mass_fraction']
    fish_kcal_day = basis['literature']['eat_lancet_kcal_day']['value']['fish'] * f['fish_edible_g_person_day'] / basis['literature']['eat_lancet_g_day']['value']['fish']
    plant_kcal_person_day = f['intake_kcal_person_day'] - fish_kcal_day - f['microbial_protein_g_person_day'] * f['protein_kcal_g']
    plant_kcal_year = plant_kcal_person_day * population * JULIAN_YEAR_DAYS
    gross_human_crop_kcal = plant_kcal_year * f['crop_loss_factor'] * f['variety_factor']
    feed_crop_kcal = feed * f['feed_kcal_kg'] * (1-f['feed_non_crop_fraction']) * f['crop_loss_factor']
    crop_kcal = gross_human_crop_kcal + feed_crop_kcal
    active_bed = crop_kcal / JULIAN_YEAR_DAYS / f['crop_productivity_kcal_m2_day']
    installed_bed = active_bed * f['crop_capacity_factor']
    floor = installed_bed / f['crop_layers'] * f['crop_floor_service_factor']
    light_W = active_bed * f['crop_light_mol_m2_day'] / (f['crop_LED_umol_J'] * 1e-6) / JULIAN_DAY
    crop_W = light_W * (1+f['crop_auxiliary_fraction_of_light'])
    fish_W = fish_live * basis['derived']['fish']['kwh_per_kg_fish'] * 1000 / YEAR_HOURS
    microbial_W = microbial_protein * f['microbial_kWh_kg_protein'] * 1000 / YEAR_HOURS
    fresh_crop_kg = gross_human_crop_kcal / f['crop_edible_dry_kcal_kg'] / (1-f['fresh_crop_water_mass_fraction'])
    processing_W = fresh_crop_kg * f['food_handling_kWh_kg'] * 1000 / YEAR_HOURS
    tanks = fish_live * f['fish_capacity_factor'] / f['fish_production_kg_m3_year']
    ras_floor = tanks / f['fish_tank_depth_m'] * f['fish_floor_service_factor']
    reserve = population*f['intake_kcal_person_day']*f['stored_food_days']/f['stored_food_kcal_kg']*(1+f['stored_food_extra_fraction'])
    crop_water = active_bed * f['crop_water_flux_kg_m2_day'] * JULIAN_YEAR_DAYS
    ras_purge = fish_live * basis['derived']['fish']['water_l_per_kg_fish']
    nutrient = {'N':(population*f['food_N_kg_person_year']+feed*f['feed_N_fraction'])*(1-f['nutrient_recovery']),
                'P':(population*f['food_P_kg_person_year']+feed*f['feed_P_fraction'])*(1-f['nutrient_recovery']),
                'K':population*f['food_K_kg_person_year']*(1-f['nutrient_recovery']),
                'trace':population*f['trace_nutrients_kg_person_year'],
                'carbon':population*f['carbon_makeup_kg_person_year']}
    nutrient_delivered = nutrient['carbon']/f['carbon_feed_carbon_mass_fraction'] + sum(v for k,v in nutrient.items() if k != 'carbon')*f['nutrient_salt_mass_factor']
    working_feed = population*(f['working_carbon_kg_person']/f['carbon_feed_carbon_mass_fraction']+f['working_NPK_kg_person']*f['nutrient_salt_mass_factor'])
    return dict(population=population,edible_fish_kg_year=fish_edible,live_fish_kg_year=fish_live,feed_kg_year=feed,
                microbial_protein_kg_year=microbial_protein,microbial_biomass_kg_year=microbial_biomass,
                plant_intake_kcal_person_day=plant_kcal_person_day,fish_intake_kcal_person_day=fish_kcal_day,
                microbial_counted_kcal_person_day=f['microbial_protein_g_person_day']*f['protein_kcal_g'],
                diet_protein_g_person_day=plant_kcal_person_day/f['crop_edible_dry_kcal_kg']*f['plant_edible_protein_fraction']*1000+f['fish_edible_g_person_day']*basis['derived']['fish']['fillet_protein_share']+f['microbial_protein_g_person_day'],
                crop_kcal_year=crop_kcal,feed_crop_kcal_year=feed_crop_kcal,active_crop_bed_m2=active_bed,
                installed_crop_bed_m2=installed_bed,farm_floor_m2=floor,fish_tanks_m3=tanks,fish_floor_m2=ras_floor,
                fermentation_floor_m2=population*f['microbial_floor_m2_person'],food_handling_floor_m2=population*f['food_handling_floor_m2_person'],
                fresh_crop_kg_year=fresh_crop_kg,stored_food_kg=reserve,
                lighting_W=light_W,crop_W=crop_W,fish_W=fish_W,microbial_W=microbial_W,processing_W=processing_W,
                mean_food_W=crop_W+fish_W+microbial_W+processing_W,
                capacity_food_W=(crop_W+microbial_W+processing_W)*f['crop_capacity_factor']+fish_W*f['fish_capacity_factor'],
                crop_water_throughput_kg_year=crop_water,crop_water_makeup_kg_year=crop_water*(1-f['crop_water_recovery']),
                ras_purge_kg_year=ras_purge,ras_water_makeup_kg_year=ras_purge*(1-f['ras_purge_recovery']),
                nutrient_element_makeup_kg_year=nutrient,nutrient_feed_delivered_kg_year=nutrient_delivered,initial_working_feed_kg=working_feed,one_city_out_crop_capacity_fraction=(s['cities']-1)/s['cities']*f['crop_capacity_factor'],
                stock_production_electricity_kWh=reserve*f['stored_food_kcal_kg']/crop_kcal*crop_W/1000*YEAR_HOURS)


def geometry(s, food):
    g=s['geometry']; cities=s['cities']; people=cities*s['city_population']; total=food['population']
    floor={'ordinary':people*g['ordinary_floor_m2_person'],'parks':people*g['park_floor_m2_city_resident'],
           'municipal':people*g['municipal_floor_m2_city_resident'],'agriculture':food['farm_floor_m2'],
           'aquaculture':food['fish_floor_m2'],'fermentation':food['fermentation_floor_m2'],
           'food_handling':food['food_handling_floor_m2'],'light_industry':total*s['industry']['production_floor_m2_person']}
    radius=g['pressure_hall_diameter_m']/2
    # A triangle bound in the elliptical norm proves cross-section containment
    # for the assumed aligned circular halls, including their lateral clearance.
    outer_center=(g['lanes']-1)*(2*radius+g['pressure_lane_clearance_m'])/2
    a=g['cave_width_m']/2;b=g['cave_height_m']
    packing_bound=math.hypot(outer_center/a,radius/b)+radius/min(a,b)
    if packing_bound>1 or any(abs(z)>=1 for z in g['deck_height_fractions']):
        raise ValueError('Pressure halls or decks do not fit the assumed cave section')
    chord=sum(2*radius*math.sqrt(1-z*z) for z in g['deck_height_fractions'])
    unit_floor=chord*g['isolation_length_m']*g['usable_floor_fraction']
    modules_lane=math.ceil(sum(floor.values())/cities/g['lanes']/unit_floor)
    modules=modules_lane*g['lanes']*cities
    cylinder_length=modules*g['isolation_length_m']
    # Each module closes its loads through two hemispherical ends. Their axial
    # length is included in the required natural corridor, without counting head floor.
    cave_length=modules_lane*(g['isolation_length_m']+2*radius)
    cylinder_area=2*math.pi*radius*cylinder_length
    head_area=modules*4*math.pi*radius**2
    multiplier=1+g['connector_surface_fraction']
    area=(cylinder_area+head_area)*multiplier
    # Orthogonal composite membrane directions: hoop plus axial in a cylinder,
    # two equal tangential directions in a sphere. Detail factor is additional.
    cylinder_areal=1.5*STANDARD_ATMOSPHERE*radius/g['effective_allowable_Pa']*g['reinforcement_density_kg_m3']
    head_areal=STANDARD_ATMOSPHERE*radius/g['effective_allowable_Pa']*g['reinforcement_density_kg_m3']
    reinforcement=(cylinder_area*cylinder_areal+head_area*head_areal)*g['pressure_detail_factor']*multiplier
    a=g['cave_width_m']/2;b=g['cave_height_m'];h=(a-b)**2/(a+b)**2
    perimeter=math.pi*(a+b)*(1+3*h/(10+math.sqrt(4-3*h)))/2+2*a
    cave_area=perimeter*cave_length*cities
    volume=(math.pi*radius**2*cylinder_length+modules*4/3*math.pi*radius**3)*multiplier
    air=STANDARD_ATMOSPHERE*g['air_molar_mass_kg_mol']/(GAS_CONSTANT*g['air_temperature_K'])*volume
    return dict(floor_demand_m2=floor,required_floor_m2=sum(floor.values()),available_floor_m2=modules*unit_floor,
                pressure_modules=modules,modules_per_city=modules//cities,modules_per_lane=modules_lane,
                equivalent_cave_length_m_per_city=cave_length,cave_width_m=g['cave_width_m'],cave_height_m=g['cave_height_m'],
                pressure_hall_diameter_m=2*radius,cross_section_containment_bound=packing_bound,total_envelope_area_m2=area,cave_surface_treated_m2=cave_area,
                pressurized_volume_m3=volume,air_operating_kg=air,air_with_reserve_kg=air*(1+g['gas_reserve_fraction']),
                pressure_reinforcement_kg=reinforcement,pressure_liner_kg=area*g['liner_kg_m2'],
                independent_support_kg=area*g['support_kg_m2_envelope'],rock_stabilization_kg=cave_area*g['rock_stabilization_kg_m2'],
                grading_kg=g['cave_width_m']*cave_length*cities*g['grading_depth_m']*g['rock_density_kg_m3'],
                module_end_force_N=STANDARD_ATMOSPHERE*math.pi*radius**2,
                overburden_weight_pressure_Pa=g['roof_cover_m']*g['rock_density_kg_m3']*MOON_SURFACE_GRAVITY,
                overburden_ratio_to_atmosphere=g['roof_cover_m']*g['rock_density_kg_m3']*MOON_SURFACE_GRAVITY/STANDARD_ATMOSPHERE,
                rule='Overburden ratio is a weight diagnostic only, not a stability criterion. Independent pressure shells do not claim natural rock as the pressure vessel.')


def replacement_count(years, life):
    if years<=0 or life<=0:raise ValueError('Positive horizon and asset life required')
    return max(0,math.ceil(years/life)-1)


def scenario(s, source, price_name='reference', wear='reference'):
    old=source['towns']['cases'][s['town_case']];q=old['quantities'];n=q['communities'];tp=source['town_inputs']['prices'][price_name]
    p=s['prices'][price_name];g=s['geometry'];u=s['utilities'];r=s['resources'];ind=s['industry'];m=s['maintenance']
    city_people=s['cities']*s['city_population'];population=city_people+q['population']
    f=food_system(s,source,population);geo=geometry(s,f)
    allocation=ind['ordinary_capital_food_share']+ind['ordinary_capital_industry_share']
    if not 0<=allocation<1:raise ValueError('Production allocation must lie in [0,1)')
    ordinary_price=source['economy']['cost_cases'][tp['ordinary_case']]['ordinary_installed_capital_U_per_resident']
    price_scale=ordinary_price/source['economy']['cost_cases']['advanced_industrial']['ordinary_installed_capital_U_per_resident']
    direct_old=old['costs'][price_name]['direct_U']
    credits={'embedded_production_capital':direct_old['ordinary_installed_capital']*allocation,
             'embedded_production_extra_freight':direct_old['ordinary_extra_delivery']*allocation,
             'dedicated_power':direct_old['dedicated_firm_power'],'storage':direct_old['emergency_reserves'],
             'radiators':direct_old['exterior_radiators'],'construction_electricity':direct_old['construction_utility_operation']}
    town_direct=sum(direct_old.values())-sum(credits.values())
    goods=population*ind['goods_kg_person_year'];virgin_goods=goods*(1-ind['goods_material_recovery'])
    industry_W=goods*ind['goods_specific_kWh_kg']*1000/YEAR_HOURS
    recovery_W=goods*ind['goods_material_recovery']*ind['recovery_service_kWh_kg_processed']*1000/YEAR_HOURS
    park_mean=geo['floor_demand_m2']['parks']*u['garden_average_W_m2']
    town_mean=q['services']['mean_occupied_W']-q['services']['food_W']
    town_peak=q['services']['peak_occupied_W']-q['services']['food_W']
    city_basic=city_people*u['basic_kW_city_resident']*1000
    mean_end_use=town_mean+city_basic+park_mean+f['mean_food_W']+industry_W+recovery_W
    peak_end_use=town_peak+city_basic+park_mean*u['garden_peak_to_mean']+f['capacity_food_W']+industry_W+recovery_W
    mean_bus=mean_end_use*(1+u['cooling_heat_pump_electric_fraction'])
    occupied_firm=peak_end_use*(1+u['cooling_heat_pump_electric_fraction'])*(1+u['firm_redundancy_fraction'])
    total_firm=occupied_firm+u['construction_power_W']
    heat=mean_bus/u['receiver_efficiency']
    emissivity=source['vacuum']['operation']['radiator_emissivity']
    flux=emissivity*STEFAN_BOLTZMANN*u['radiator_K']**4
    rad_area=total_firm/u['receiver_efficiency']/flux*u['radiator_capacity_factor']
    retained_share=occupied_firm/total_firm
    storage_factor=(1+source['vacuum']['operation']['storage_reserve_fraction'])/(source['vacuum']['operation']['storage_discharge_efficiency']*source['vacuum']['operation']['storage_usable_fraction'])
    storage=population*(u['critical_kW_person']+u['food_preservation_critical_kW_person'])*1000*u['critical_reserve_hours']*storage_factor
    route_mass=s['network']['additional_route_km']*1000*s['network']['route_width_m']*s['network']['route_layer_m']*s['network']['route_density_kg_m3']
    ordinary_city_stock=city_people*r['ordinary_stock_kg_person']*(1-allocation)
    shared_stock=f['farm_floor_m2']*s['food']['crop_structural_kg_m2_floor']+f['installed_crop_bed_m2']*s['food']['crop_equipment_kg_m2_bed']+(f['fish_floor_m2']+f['fermentation_floor_m2']+f['food_handling_floor_m2'])*ind['food_service_structure_kg_m2']+geo['floor_demand_m2']['light_industry']*ind['plant_and_structure_kg_m2']
    new_local_bulk=(ordinary_city_stock+shared_stock)*r['local_bulk_share']+route_mass
    new_local_metal=(ordinary_city_stock+shared_stock)*r['local_metal_share']+geo['independent_support_kg']+geo['rock_stabilization_kg']+s['cities']*s['network']['access_local_metal_kg_city']+rad_area*source['town_inputs']['quantity']['radiator_local_metal_kg_m2']
    new_raw=new_local_bulk/r['bulk_yield']+new_local_metal/r['metal_yield']
    new_tails=new_raw-new_local_bulk-new_local_metal
    old_removed_stock=q['ordinary_stock_kg']*allocation
    old_rad_metal=q['services']['radiator_emitting_area_m2']*source['town_inputs']['quantity']['radiator_local_metal_kg_m2']
    removed_bulk=old_removed_stock*r['local_bulk_share'];removed_metal=old_removed_stock*r['local_metal_share']+old_rad_metal
    removed_raw=removed_bulk/r['bulk_yield']+removed_metal/r['metal_yield']
    town_raw=q['raw_feed_kg']-removed_raw
    conversion_kWh=(new_local_bulk-removed_bulk)*r['bulk_MJ_kg']/3.6+(new_local_metal-removed_metal)*r['metal_kWh_kg']+(new_raw-removed_raw+geo['grading_kg'])*r['handling_MJ_kg']/3.6
    construction_kWh=sum(q['construction_energy_kWh'].values())+conversion_kWh+f['stock_production_electricity_kWh']
    city_water=city_people*u['city_water_inventory_kg_person']+f['installed_crop_bed_m2']*s['food']['crop_water_inventory_kg_m2_bed']+f['fish_tanks_m3']*source['town_inputs']['quantity']['water_density_kg_m3']
    city_biosphere=city_people*u['city_biosphere_kg_person']+population*s['food']['seed_nursery_kg_person']
    new_pressure=geo['pressure_reinforcement_kg']+geo['pressure_liner_kg']
    additions={
      'city_ordinary_nonproduction_capital':city_people*ordinary_price*(1-allocation),
      'city_ordinary_extra_delivery':ordinary_city_stock*r['import_share']*tp['ordinary_import_freight_increment_U_kg'],
      'pressure_reinforcement':geo['pressure_reinforcement_kg']*tp['pressure_system_U_kg'],
      'pressure_liners':geo['pressure_liner_kg']*tp['pressure_system_U_kg'],
      'independent_support':geo['independent_support_kg']*tp['incremental_local_structure_U_kg'],
      'rock_stabilization':geo['rock_stabilization_kg']*tp['incremental_local_structure_U_kg'],
      'grading':geo['grading_kg']*tp['shield_installed_U_kg'],
      'escape_and_access':s['cities']*s['network']['city_escape_access_sets']*s['network']['access_set_U_reference']*price_scale,
      'regional_routes':s['network']['additional_route_km']*tp['regional_road_U_km'],
      'transport_and_docking':population*s['network']['city_transport_allowance_U_person_reference']*price_scale,
      'farm_structure':f['farm_floor_m2']*p['food_floor_U_m2'],
      'farm_equipment_and_lighting':f['installed_crop_bed_m2']*p['food_bed_U_m2'],
      'aquaculture':f['live_fish_kg_year']*s['food']['fish_capacity_factor']*p['ras_U_kg_year_capacity'],
      'fermentation':f['microbial_biomass_kg_year']*s['food']['crop_capacity_factor']*p['microbial_U_kg_year_capacity'],
      'food_processing':f['food_handling_floor_m2']*p['food_handling_U_m2'],
      'light_industry':geo['floor_demand_m2']['light_industry']*p['industry_U_m2'],
      'integrated_firm_power':total_firm*tp['firm_power_U_W'],
      'integrated_storage':storage*tp['reserve_U_Wh'],
      'integrated_radiators':rad_area*tp['radiator_U_m2'],
      'city_garden_lighting':geo['floor_demand_m2']['parks']*tp['garden_light_U_m2'],
      'initial_air_water':(geo['air_with_reserve_kg']+city_water)*tp['volatile_delivered_U_kg'],
      'initial_biosphere':city_biosphere*tp['biosphere_delivered_U_kg'],
      'initial_local_food_stock':f['stored_food_kg']*p['initial_food_variable_U_kg'],
      'initial_carbon_and_nutrient_stocks':f['initial_working_feed_kg']*p['nutrient_U_kg'],
      'survey_and_qualification':p['survey_qualification_U'],
      'regional_startup':p['regional_startup_U'],
      'construction_utility_operation':construction_kWh/1000*tp['construction_utility_variable_U_MWh'],
      'disturbed_ground_restoration':(new_raw+geo['grading_kg']-removed_raw)*tp['restoration_U_kg_raw_feed']}
    total_direct=town_direct+sum(additions.values())
    initial=total_direct*(1+tp['contingency_fraction']+tp['removal_reserve_fraction'])
    tm=source['town_maintenance']['cases'][s['town_case']][price_name][wear]
    asset={};replacement={};life={}
    for name,value in tm['maintained_installed_value_U'].items():
        if name in ('firm_power','emergency_storage','radiators'):continue
        factor=(1-allocation) if name.startswith('ordinary_') else 1
        key='town_'+name;asset[key]=value*n*factor;replacement[key]=tm['full_replacement_value_U'][name]*n*factor;life[key]=tm['asset_lives_years'][name]
    life_factor={'longer_life':1.25,'reference':1,'higher_wear':0.625}[wear]
    routine_fraction=source['maintenance_inputs']['regimes'][wear]['routine_upkeep_fraction_year']
    def add(name,value,extra=0):
        key='new_'+name;asset[key]=value;replacement[key]=value+extra;life[key]=m['life_years'][name]*life_factor
    ordinary=additions['city_ordinary_nonproduction_capital']
    add('ordinary_durable',ordinary*source['maintenance_inputs']['ordinary_durable_capital_fraction'])
    add('ordinary_equipment',ordinary*(1-source['maintenance_inputs']['ordinary_durable_capital_fraction']),additions['city_ordinary_extra_delivery'])
    for name in ('pressure_reinforcement','pressure_liners','independent_support','rock_stabilization','farm_structure','aquaculture','fermentation','food_processing'):
        add(name,additions[name])
    add('access',additions['escape_and_access']);add('roads',additions['regional_routes']);add('transport',additions['transport_and_docking'])
    lamp_share=m['farm_lighting_cost_share'];ind_share=m['industry_structure_cost_share']
    add('farm_equipment',additions['farm_equipment_and_lighting']*(1-lamp_share));add('farm_lighting',additions['farm_equipment_and_lighting']*lamp_share)
    add('light_industry_structure',additions['light_industry']*ind_share);add('light_industry_equipment',additions['light_industry']*(1-ind_share))
    add('power',additions['integrated_firm_power']*retained_share);add('storage',additions['integrated_storage']);add('radiators',additions['integrated_radiators']*retained_share);add('garden_lighting',additions['city_garden_lighting'])
    # On-site structural renewal purchases conversion electricity as a separate
    # supplier service. Food/industry installed replacements already include it.
    local_replacement_mass={'independent_support':geo['independent_support_kg'],'rock_stabilization':geo['rock_stabilization_kg'],'roads':route_mass,'access':s['cities']*s['network']['access_local_metal_kg_city'],'radiators':rad_area*retained_share*source['town_inputs']['quantity']['radiator_local_metal_kg_m2']}
    for key,mass in local_replacement_mass.items():
        metal=key!='roads';yield_fraction=r['metal_yield'] if metal else r['bulk_yield']
        ekWh=mass*(r['metal_kWh_kg'] if metal else r['bulk_MJ_kg']/3.6)+mass/yield_fraction*r['handling_MJ_kg']/3.6
        replacement['new_'+key]+=ekWh/1000*source['maintenance_inputs']['replacement_utility_U_MWh'][price_name]
    routine=sum(asset.values())*routine_fraction
    renewal={k:replacement[k]/life[k] for k in asset}
    variable={'town_'+k:v*n for k,v in tm['variable_services_U_year'].items() if k!='utility_operation'}
    domestic_water=city_people*source['vacuum']['operation']['water_loop_throughput_kg_per_person_day']*(1-source['maintenance_inputs']['regimes'][wear]['processed_water_recovery_fraction'])*JULIAN_YEAR_DAYS
    city_garden_water=city_people*source['maintenance_inputs']['regimes'][wear]['additional_water_makeup_kg_person_day']*JULIAN_YEAR_DAYS
    industry_water=population*ind['process_water_kg_person_day']*(1-ind['process_water_recovery'])*JULIAN_YEAR_DAYS
    new_water=domestic_water+city_garden_water+industry_water+f['crop_water_makeup_kg_year']+f['ras_water_makeup_kg_year']
    air_loss=geo['air_operating_kg']*source['maintenance_inputs']['regimes'][wear]['annual_air_loss_fraction']
    food_staff=population/1000*s['food']['farm_staff_per_1000_people'];industry_staff=population/1000*ind['staff_per_1000_people']
    variable.update(electricity=mean_bus/1e6*YEAR_HOURS*tp['construction_utility_variable_U_MWh'],
       new_air_makeup=air_loss*tp['volatile_delivered_U_kg'],new_water_makeup=new_water*tp['volatile_delivered_U_kg'],
       nutrient_carbon_makeup=f['nutrient_feed_delivered_kg_year']*p['nutrient_U_kg'],
       city_biosphere_makeup=city_biosphere*source['maintenance_inputs']['regimes'][wear]['annual_biosphere_import_fraction']*tp['biosphere_delivered_U_kg'],
       everyday_local_material_feed=virgin_goods*(1-ind['carbon_polymer_share_of_virgin_goods'])*tp['incremental_local_structure_U_kg'],
       everyday_imported_carbon_feed=virgin_goods*ind['carbon_polymer_share_of_virgin_goods']*p['nutrient_U_kg'],
       everyday_process_consumables_and_services=goods*p['goods_processing_U_kg'],
       specialty_consumable_imports=population*ind['special_import_kg_person_year']*p['special_import_U_kg'],
       production_operators=(food_staff+industry_staff)*p['staff_U_person_year'],
       basic_civic_services=population*p['civic_U_person_year'])
    recurring_before_margin=routine+sum(variable.values());margin=recurring_before_margin*m['operating_margin_fraction'];recurring=recurring_before_margin+margin
    fund=sum(renewal.values());annual=recurring+fund
    horizons={str(y):{'recurring_U':recurring*y,'actual_replacements_U':sum(replacement[k]*replacement_count(y,life[k]) for k in asset),
       'initial_plus_vacuum_cash_U':initial+recurring*y+sum(replacement[k]*replacement_count(y,life[k]) for k in asset)} for y in m['horizons_years']}
    farm_mass=shared_stock;new_special=ordinary_city_stock*r['import_share']+farm_mass*r['import_share']+new_pressure+total_firm*source['town_inputs']['quantity']['power_hardware_import_kg_W']+storage*source['town_inputs']['quantity']['reserve_import_kg_Wh']+rad_area*source['town_inputs']['quantity']['radiator_controls_import_kg_m2']+r['additional_construction_plant_import_kg']+population*s['network']['transport_import_kg_person']+s['cities']*s['network']['access_import_kg_city']+geo['floor_demand_m2']['parks']*source['town_inputs']['quantity']['garden_lamp_import_kg_m2']
    retained_town_special=sum(q['special_imports'].values())-old_removed_stock*r['import_share']-q['special_imports']['power_hardware_kg']-q['special_imports']['reserve_hardware_kg']-q['special_imports']['radiator_controls_kg']
    # Preserve the inherited town local-stock identity, removing only the credited
    # stock and old radiator metal, then add each new physical material once.
    local_total=q['local_finished_kg']-removed_bulk-removed_metal+new_local_bulk+new_local_metal
    raw_total=town_raw+new_raw;initial_volatiles=q['air_including_reserve_kg']+sum(q['water_inventory'].values())+geo['air_with_reserve_kg']+city_water
    special_total=retained_town_special+new_special+city_biosphere+f['initial_working_feed_kg']
    town_water=(tm['variable_services_U_year']['processed_water_makeup']+tm['variable_services_U_year']['additional_water_makeup'])*n/tp['volatile_delivered_U_kg']
    town_air=tm['variable_services_U_year']['air_makeup']*n/tp['volatile_delivered_U_kg']
    bio_makeup=(variable['town_biosphere_makeup']+variable['city_biosphere_makeup'])/tp['biosphere_delivered_U_kg']
    operating_import={'water':new_water+town_water,'air':air_loss+town_air,'biosphere':bio_makeup,
       'food_nutrients_and_carbon':f['nutrient_feed_delivered_kg_year'],
       'everyday_carbon_feed':virgin_goods*ind['carbon_polymer_share_of_virgin_goods'],
       'specialist_consumables':population*ind['special_import_kg_person_year']}
    # Physical replacement masses use the same expected service lives as the U
    # account. Commissioning-only construction plant is not renewed here.
    import_mass={
       'town_ordinary_equipment':q['special_imports']['ordinary_equipment_kg']*(1-allocation),
       'town_pressure_reinforcement':q['pressure_reinforcement_kg'],
       'town_pressure_barriers':q['barriers_kg'],
       'town_transport_docking':q['special_imports']['transport_and_docking_kg'],
       'town_garden_lighting':q['special_imports']['garden_lamps_kg'],
       'new_ordinary_equipment':ordinary_city_stock*r['import_share'],
       'new_pressure_reinforcement':geo['pressure_reinforcement_kg'],
       'new_pressure_liners':geo['pressure_liner_kg'],
       'new_farm_structure':f['farm_floor_m2']*s['food']['crop_structural_kg_m2_floor']*r['import_share'],
       'new_farm_equipment':f['installed_crop_bed_m2']*s['food']['crop_equipment_kg_m2_bed']*r['import_share']*(1-lamp_share),
       'new_farm_lighting':f['installed_crop_bed_m2']*s['food']['crop_equipment_kg_m2_bed']*r['import_share']*lamp_share,
       'new_aquaculture':f['fish_floor_m2']*ind['food_service_structure_kg_m2']*r['import_share'],
       'new_fermentation':f['fermentation_floor_m2']*ind['food_service_structure_kg_m2']*r['import_share'],
       'new_food_processing':f['food_handling_floor_m2']*ind['food_service_structure_kg_m2']*r['import_share'],
       'new_light_industry_structure':geo['floor_demand_m2']['light_industry']*ind['plant_and_structure_kg_m2']*r['import_share']*ind_share,
       'new_light_industry_equipment':geo['floor_demand_m2']['light_industry']*ind['plant_and_structure_kg_m2']*r['import_share']*(1-ind_share),
       'new_power':occupied_firm*source['town_inputs']['quantity']['power_hardware_import_kg_W'],
       'new_storage':storage*source['town_inputs']['quantity']['reserve_import_kg_Wh'],
       'new_radiators':rad_area*retained_share*source['town_inputs']['quantity']['radiator_controls_import_kg_m2'],
       'new_garden_lighting':geo['floor_demand_m2']['parks']*source['town_inputs']['quantity']['garden_lamp_import_kg_m2'],
       'new_transport':population*s['network']['transport_import_kg_person'],
       'new_access':s['cities']*s['network']['access_import_kg_city']}
    renewal_import={k:v/life[k] for k,v in import_mass.items()}
    total_floor=geo['required_floor_m2']+q['population']*g['ordinary_floor_m2_person']
    demand={'handling':q['production_demand']['handling_kg']-removed_raw+new_raw+geo['grading_kg'],
            'local_bulk':q['production_demand']['bulk_finished_kg']-removed_bulk+new_local_bulk,
            'local_metal':q['production_demand']['metal_finished_kg']-removed_metal+new_local_metal,
            'pressure_components':new_pressure+q['pressure_system_kg'],'other_imports':special_total-new_pressure-q['pressure_system_kg'],
            'volatile_imports':initial_volatiles,'floor':total_floor,'route':s['network']['additional_route_km']+source['town_inputs']['quantity']['regional_road_km']}
    years={k:v/s['production'][k+('_m2_year' if k=='floor' else '_km_year' if k=='route' else '_kg_year')] for k,v in demand.items()}
    bottleneck=max(years,key=years.get)
    elapsed=s['production']['lead_years']+s['production']['ramp_allowance_years']+max(years.values())+s['production']['final_commissioning_years']
    return dict(price_case=price_name,wear_case=wear,population=population,city_population=city_people,town_population=q['population'],
       food=f,geometry=geo,utilities={'mean_load_components_W':{'food':f['mean_food_W'],'basic_services':city_basic+q['services']['basic_W'],'civic_gardens':park_mean+q['services']['garden_average_W'],'everyday_industry':industry_W+recovery_W,'cooling':mean_bus-mean_end_use},'mean_end_use_W':mean_end_use,'mean_electric_bus_W':mean_bus,'mean_heat_rejection_W':heat,'occupied_firm_W':occupied_firm,'construction_firm_W':u['construction_power_W'],'total_firm_W':total_firm,'reserve_Wh':storage,'radiator_emitting_m2':rad_area,'occupied_radiator_emitting_m2':rad_area*retained_share,'full_dark_half_cycle_storage_Wh':mean_bus*SYNODIC_MONTH_DAYS*24/2*storage_factor},
       construction={'original_town_initial_U':old['costs'][price_name]['total_initial_funding_U'],'town_direct_credits_U':credits,'retained_town_direct_U':town_direct,'additions_direct_U':additions,'direct_U':total_direct,'contingency_U':total_direct*tp['contingency_fraction'],'transition_reserve_U':total_direct*tp['removal_reserve_fraction'],'initial_U':initial,'net_addition_to_original_towns_U':initial-old['costs'][price_name]['total_initial_funding_U'],'initial_U_resident':initial/population},
       maintenance={'assets_U':asset,'replacement_value_U':replacement,'life_years':life,'routine_upkeep_U_year':routine,'variable_operations_U_year':variable,'operating_margin_U_year':margin,'recurring_U_year':recurring,'renewal_U_year':renewal,'renewal_total_U_year':fund,'total_U_year':annual,'U_resident_year':annual/population,'horizons':horizons},
       resources={'local_finished_kg':local_total,'raw_feed_kg':raw_total,'additional_retained_tails_kg':new_tails,'additional_grading_kg':geo['grading_kg'],'special_imports_kg':special_total,'initial_air_water_kg':initial_volatiles,'local_fraction_delivered':local_total/(local_total+special_total+initial_volatiles),'new_processed_raw_feed_200year_reservation_kg':new_raw*(1+r['horizon_years']/r['replacement_cycle_years']*(1-r['end_to_end_return_fraction'])),'quarry_area_m2':raw_total/g['rock_density_kg_m3']/r['quarry_working_depth_m'],'new_route_corridor_m2':s['network']['additional_route_km']*1000*s['network']['route_corridor_width_m'],'construction_electricity_kWh':construction_kWh,'everyday_goods_kg_year':goods,'virgin_goods_kg_year':virgin_goods,'new_water_makeup_kg_year':new_water,'new_air_makeup_kg_year':air_loss,'air_inventory_cycle_years':geo['air_operating_kg']/air_loss if air_loss else None,'operating_import_kg_year':operating_import,'renewal_import_kg_year':renewal_import,'total_annual_import_kg':sum(operating_import.values())+sum(renewal_import.values()),'initial_local_food_kg':f['stored_food_kg'],'local_goods_makeup_kg_year':virgin_goods*(1-ind['carbon_polymer_share_of_virgin_goods']),'food_operators':food_staff,'industry_operators':industry_staff},
       schedule={'production_demand':demand,'capacity_years':years,'bottleneck':bottleneck,'completion_screen_years':elapsed,'rule':'Throughput diagnostic with assumed lead, ramp and final commissioning; simultaneous supplier readiness and routes are conditions, not a construction forecast.'})


def run():
    spec,source=read_inputs()
    cases={name:scenario(spec,source,name,wear) for name,wear in [('lower','longer_life'),('reference','reference'),('upper','higher_wear')]}
    sensitivity={}
    changes={'crop_yield_half':('food','crop_productivity_kcal_m2_day',25),
             'LED_2_5':('food','crop_LED_umol_J',2.5),'allowable_400MPa':('geometry','effective_allowable_Pa',400e6),
             'routes_6000km':('network','additional_route_km',6000),'radiators_300K':('utilities','radiator_K',300),
             'crop_water_98percent':('food','crop_water_recovery',0.98)}
    for name,(section,key,value) in changes.items():
        changed=deepcopy(spec);changed[section][key]=value;c=scenario(changed,source)
        sensitivity[name]={'initial_U':c['construction']['initial_U'],'annual_U':c['maintenance']['total_U_year'],'mean_electric_bus_W':c['utilities']['mean_electric_bus_W'],'cave_length_m_per_city':c['geometry']['equivalent_cave_length_m_per_city'],'annual_import_kg':c['resources']['total_annual_import_kg']}
    changed=deepcopy(spec);changed['industry']['ordinary_capital_food_share']=0;changed['industry']['ordinary_capital_industry_share']=0
    c=scenario(changed,source);sensitivity['no_generic_production_credit']={'initial_U':c['construction']['initial_U'],'annual_U':c['maintenance']['total_U_year']}
    ref=cases['reference'];town_prices=source['town_inputs']['prices']['reference']
    sensitivity['town_capital_already_spent_no_recovery']={
        'initial_U':ref['construction']['initial_U']+sum(ref['construction']['town_direct_credits_U'].values())*(1+town_prices['contingency_fraction']+town_prices['removal_reserve_fraction']),
        'rule':'Prior-spending diagnostic: retain the complete original town bill, then buy the integrated end state without salvage recovery. Temporary bridge operation is additional; this is not a construction-sequence simulation.'}
    files=['research/studies/lava_tube_cities/'+p for p in ['buildout.py','buildout_inputs.json','buildout_sources.json']]
    return {'schema':SCHEMA,'evidence':spec['evidence'],'date':spec['date'],'numeraire':source['towns']['numeraire'],'inputs':spec['inputs'],
            'producer':{'files':{p:sha(ROOT/p) for p in files},'constants':constants_used(files)},'cases':cases,'sensitivities':sensitivity,'rules':spec['rules']}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=ROOT/'research/runs/lava_tube_buildout.json');args=parser.parse_args()
    result=run();args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    for key,c in result['cases'].items():
        print(key,round(c['construction']['initial_U']/1e9,3),'billion U initial;',round(c['maintenance']['total_U_year']/1e9,3),'billion U/year;',round(c['utilities']['mean_electric_bus_W']/1e9,3),'GW')


if __name__=='__main__':main()
