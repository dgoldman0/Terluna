"""Owned industrial supply and procurement screen; no validated mine or trajectory.

All payloads, including industrial replacement hardware, drive the propellant,
source, power and fleet accounts. The fixed point closes those requirements.
Internal transfer tariffs are diagnostics and never added to the owner budget.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
from shared.constants import (JULIAN_DAY, JULIAN_YEAR_DAYS, STANDARD_GRAVITY,
    MOON_RADIUS, MOON_GM, SOLAR_CONSTANT, SPEED_OF_LIGHT, STEFAN_BOLTZMANN)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.initial-industrial-system/1'
HOURS = JULIAN_DAY * JULIAN_YEAR_DAYS / 3600


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_inputs(spec=None):
    s = spec if spec is not None else json.loads((HERE/'industrial_inputs.json').read_text())
    source = {}
    for k, b in s['inputs'].items():
        p = ROOT/b['path']
        if sha(p) != b['sha256']:
            raise ValueError(f'Input hash mismatch: {k}')
        d = json.loads(p.read_text())
        if 'schema' in b and d.get('schema') != b['schema']:
            raise ValueError(f'Input schema mismatch: {k}')
        source[k] = d
    return s, source


def with_case(s, name):
    q = deepcopy(s)
    for path, value in s['case_overrides'][name].items():
        a, b = path.split('.')
        q[a][b] = value
    return q


def round_trip(dv, isp, dry_ratio, return_payload_ratio=0):
    """Fuel for outbound unit cargo and an empty or hardware-loaded return.

All return fuel travels outbound. Return payload is picked up at destination.
Separate loss allowance is added by the caller, not used as free thrust.
"""
    if min(dv, dry_ratio, return_payload_ratio) < 0 or isp <= 0:
        raise ValueError('Invalid rocket input')
    r = math.exp(dv/(isp*STANDARD_GRAVITY))
    back = (dry_ratio+return_payload_ratio)*(r-1)
    outward = (1+dry_ratio+back)*(r-1)
    return dict(outbound=outward, returning=back, total=outward+back, mass_ratio=r)


def settlement_demand(s, src):
    c = src['settlement']['cases']['reference']
    q = src['towns']['cases'][src['settlement_inputs']['town_case']]['quantities']
    p = c['schedule']['production_demand']; r = c['resources']; l = s['lunar']
    roads = q['road_material_kg'] + src['settlement_inputs']['network']['additional_route_km']*1000*src['settlement_inputs']['network']['route_width_m']*src['settlement_inputs']['network']['route_layer_m']*src['settlement_inputs']['network']['route_density_kg_m3']
    air = q['air_including_reserve_kg']+c['geometry']['air_with_reserve_kg']
    operating = r['operating_import_kg_year']
    renewal = r['renewal_import_kg_year']
    pressure_renewal = sum(v for k,v in renewal.items() if 'pressure' in k)
    return dict(population=c['population'],bulk=p['local_bulk'],metal=p['local_metal'],
        loose_and_grading=p['handling']-p['local_bulk']/l['bulk_yield']-p['local_metal']/l['metal_yield'],
        pressure=p['pressure_components'],earth_special=r['special_imports_kg']-l['old_embedded_plant_kg'],
        water=r['initial_air_water_kg']-air,air=air,roads=roads,
        renewal_bulk=roads/l['road_life_years']+(p['local_bulk']-roads)/l['bulk_life_years'],
        renewal_metal=p['local_metal']/l['metal_life_years'],
        yearly_water=operating['water'],yearly_air=operating['air'],
        yearly_earth=sum(renewal.values())+sum(v for k,v in operating.items() if k not in ('water','air')),
        yearly_pressure=pressure_renewal,goods_virgin=r['local_goods_makeup_kg_year'],
        occupied_W=c['utilities']['occupied_firm_W'],occupied_mean_W=c['utilities']['mean_electric_bus_W'],
        old_construction_W=c['utilities']['construction_firm_W'],
        old_total_W=c['utilities']['total_firm_W'],old_radiator_m2=c['utilities']['radiator_emitting_m2'],
        reserve_Wh=c['utilities']['reserve_Wh'],
        power_hardware_kg_W=src['town_inputs']['quantity']['power_hardware_import_kg_W'],
        radiator_controls_kg_m2=src['town_inputs']['quantity']['radiator_controls_import_kg_m2'],
        radiator_metal_kg_m2=src['town_inputs']['quantity']['radiator_local_metal_kg_m2'],
        U_pressure_kg=src['town_inputs']['prices']['reference']['pressure_system_U_kg'],
        U_volatile_kg=src['town_inputs']['prices']['reference']['volatile_delivered_U_kg'],
        U_other_freight_kg=src['town_inputs']['prices']['reference']['ordinary_import_freight_increment_U_kg'])


def flows(s, d, assets, construction=True):
    """Annual rated work; existing fleet also needs renewal after occupation."""
    H = s['schedule']['delivery_equivalent_years']; f = 1/H if construction else 0
    l=s['lunar']; t=s['transport']; z=s['source']
    destination = {k:sum(a['mass_kg']*a.get('import_fraction',1)*(f+1/a['life_years']) for a in assets.values() if a['destination']==k) for k in ('surface','depot','source','leo')}
    surf_mass = sum(a['mass_kg'] for a in assets.values() if a['destination']=='surface')
    added_bulk=surf_mass*l['infrastructure_bulk_per_hardware_kg']
    added_metal=surf_mass*l['infrastructure_metal_per_hardware_kg']
    added_radiator_metal=assets.get('additional_ground_radiator_controls',{}).get('mass_kg',0)/d['radiator_controls_kg_m2']*d['radiator_metal_kg_m2']
    # Independent radiator panel metal is additional to factory frames.
    added_metal+=added_radiator_metal
    bulk=(d['bulk']+added_bulk)*f+d['renewal_bulk']+added_bulk/l['bulk_life_years']+d['goods_virgin']
    metal=(d['metal']+added_metal)*f+d['renewal_metal']+added_metal/l['metal_life_years']
    fresh_bulk=(d['bulk']+added_bulk)*f+(d['renewal_bulk']+added_bulk/l['bulk_life_years'])*(1-l['recovery'])+d['goods_virgin']
    fresh_metal=(d['metal']+added_metal)*f+(d['renewal_metal']+added_metal/l['metal_life_years'])*(1-l['recovery'])
    raw=fresh_bulk/l['bulk_yield']+fresh_metal/l['metal_yield']
    handling=raw+d['loose_and_grading']*f
    water=d['water']*f+d['yearly_water']; air=d['air']*f+d['yearly_air']
    earth_surface=d['earth_special']*f+d['yearly_earth']+destination['surface']
    moon=earth_surface+water+air
    land=round_trip(t['lander_dv_oneway_m_s'],t['Isp_s'],t['lander_dry_to_payload'])
    land_prop=moon*land['total']*t['propellant_margin']
    water_per_prop=z['water_to_hydrogen_mass_ratio']/(1+t['mixture_O_to_H'])
    if t['mixture_O_to_H'] != z['oxygen_to_hydrogen_mass_ratio']:
        raise ValueError('This water-balanced case requires nominal stoichiometric propellant')
    oxygen=air*(1-z['air_nitrogen_mass_fraction']); nitrogen=air*z['air_nitrogen_mass_fraction']
    # Electrolyze propellant mass plus air oxygen. All hydrogen enters the
    # engines, giving a small shift below the nominal mixture of eight.
    depot_water=land_prop+oxygen
    depot_hydrogen=depot_water/z['water_to_hydrogen_mass_ratio']
    depot_oxygen_to_fuel=land_prop-depot_hydrogen
    if depot_oxygen_to_fuel < 0:
        raise ValueError('Air oxygen exceeds this coupled propellant process')
    depot_oxygen_excess=oxygen
    nea=water+nitrogen+depot_water
    return_hardware=destination['source']
    nr=round_trip(t['nea_dv_oneway_m_s'],t['Isp_s'],t['nea_dry_to_payload'],return_hardware/nea)
    nea_prop=nea*nr['total']*t['propellant_margin']
    source_water=water+depot_water+nea_prop*water_per_prop
    ore_water=source_water/(z['water_grade']*z['water_recovery'])
    ore_nitrogen=nitrogen/(z['nitrogen_grade']*z['nitrogen_recovery'])
    tug=assets.get('asteroid_tugs',{'mass_kg':0,'life_years':t['nea_life_years']})
    local_tug_stock=tug['mass_kg']*z['tug_local_mass_fraction']
    source_metal=local_tug_stock*(f+1/tug['life_years'])
    source_virgin_metal=local_tug_stock*(f+(1-l['recovery'])/tug['life_years'])
    ore_metal=source_virgin_metal/(z['metal_grade']*z['metal_recovery'])
    ore=max(ore_water,ore_nitrogen,ore_metal)
    local_kWh=bulk*l['bulk_kWh_kg']+metal*l['metal_kWh_kg']+handling*l['handling_kWh_kg']
    source_process_kWh=ore*z['ore_kWh_kg']+nitrogen*z['nitrogen_kWh_kg']+source_metal*z['metal_refining_and_fabrication_kWh_kg']
    source_chem_kWh=nea_prop/(1+t['mixture_O_to_H'])*z['electrolysis_liquefaction_kWh_kg_H2']
    depot_chem_kWh=depot_hydrogen*z['electrolysis_liquefaction_kWh_kg_H2']
    chem_heat_fraction=1-z['hydrogen_HHV_kWh_kg']/z['electrolysis_liquefaction_kWh_kg_H2']
    er=round_trip(t['earth_dv_oneway_m_s'],t['Isp_s'],t['earth_dry_to_payload'])
    earth_deep=earth_surface+destination['depot']+destination['source']
    earth_prop=earth_deep*er['total']*t['propellant_margin']
    return dict(bulk_kg_year=bulk,metal_kg_year=metal,raw_lunar_kg_year=raw,handling_kg_year=handling,
        additional_bulk_stock_kg=added_bulk,additional_metal_stock_kg=added_metal,
        earth_surface_kg_year=earth_surface,earth_deep_kg_year=earth_deep,
        moon_payload_kg_year=moon,moon_water_kg_year=water,moon_air_kg_year=air,
        lunar_propellant_kg_year=land_prop,depot_water_for_propellant_kg_year=depot_water,
        depot_excess_oxygen_kg_year=depot_oxygen_excess-oxygen,depot_hydrogen_kg_year=depot_hydrogen,depot_oxygen_to_fuel_kg_year=depot_oxygen_to_fuel,depot_effective_O_to_H=depot_oxygen_to_fuel/depot_hydrogen,
        nea_payload_kg_year=nea,nea_return_hardware_kg_year=return_hardware,
        nea_propellant_kg_year=nea_prop,source_water_kg_year=source_water,source_nitrogen_kg_year=nitrogen,
        source_ore_kg_year=ore,source_ore_water_bound_kg_year=ore_water,source_ore_nitrogen_bound_kg_year=ore_nitrogen,
        source_retained_tails_kg_year=ore-source_water-nitrogen-source_virgin_metal,source_metal_kg_year=source_metal,source_virgin_metal_kg_year=source_virgin_metal,source_ore_metal_bound_kg_year=ore_metal,
        local_mean_W=local_kWh*1000/HOURS,source_mean_W=(source_process_kWh+source_chem_kWh)*1000/HOURS,
        depot_mean_W=depot_chem_kWh*1000/HOURS,
        source_heat_W=(source_process_kWh+source_chem_kWh*chem_heat_fraction)*1000/HOURS,
        depot_heat_W=depot_chem_kWh*chem_heat_fraction*1000/HOURS,
        earth_transfer_propellant_kg_year=earth_prop,LEO_kg_year=earth_deep+earth_prop+destination['leo'],
        payload_destination_kg_year=destination,earth_round_trip=er,lander_round_trip=land,nea_round_trip=nr)


def design_assets(s,d,f):
    l=s['lunar']; z=s['source']; t=s['transport']; p=s['power']; prices=s['prices']
    m=s['economics']['capacity_margin']; assets={}; counts={}
    def add(name,kg,cost,life,destination,group):
        assets[name]=dict(mass_kg=kg,manufacture_USD2026=cost,life_years=life,destination=destination,group=group)
    def machine(name,rate,unit_rate,kg,availability,life,dest,price,group):
        n=math.ceil(rate*m/(unit_rate*HOURS*availability)); counts[name]=n
        add(name,n*kg,n*kg*price,life,dest,group)
    machine('excavators',f['handling_kg_year'],l['excavator_kg_hr'],l['excavator_kg'],l['mine_availability'],l['robot_life_years'],'surface',prices['lunar_hardware_USD2026_kg'],'lunar')
    machine('haulers',f['handling_kg_year'],l['hauler_payload_kg']/l['haul_roundtrip_hr'],l['hauler_kg'],l['mine_availability'],l['robot_life_years'],'surface',prices['lunar_hardware_USD2026_kg'],'lunar')
    for name,key,rate,kg in [('bulk_trains','bulk_kg_year',l['bulk_train_kg_hr'],l['bulk_train_hardware_kg']),('metal_trains','metal_kg_year',l['metal_train_kg_hr'],l['metal_train_hardware_kg']),('fabrication_cells','metal_kg_year',l['fabrication_cell_kg_hr'],l['fabrication_cell_hardware_kg'])]:
        machine(name,f[key],rate,kg,l['process_availability'],l['plant_life_years'],'surface',prices['lunar_hardware_USD2026_kg'],'lunar')
    add('construction_and_service_robots',l['construction_plant_kg'],l['construction_plant_kg']*prices['lunar_hardware_USD2026_kg'],l['robot_life_years'],'surface','lunar')
    machine('source_mine_trains',f['source_ore_kg_year'],z['mine_module_kg_hr'],z['mine_module_hardware_kg'],z['mine_availability'],l['plant_life_years'],'source',prices['source_mine_USD2026_kg'],'source')
    machine('source_metal_trains',f['source_metal_kg_year'],l['metal_train_kg_hr'],l['metal_train_hardware_kg'],l['process_availability'],l['plant_life_years'],'source',prices['source_mine_USD2026_kg'],'source')
    machine('source_fabrication_cells',f['source_metal_kg_year'],l['fabrication_cell_kg_hr'],l['fabrication_cell_hardware_kg'],l['process_availability'],l['plant_life_years'],'source',prices['source_mine_USD2026_kg'],'source')
    for name,key,dest,group in [('source_chemistry','source_mean_W','source','source'),('depot_chemistry','depot_mean_W','depot','landing')]:
        # Includes capture, cleanup and thermal-process auxiliaries, not only electrolyzer stacks.
        cap=f[key]*m
        add(name,cap*z['chemistry_hardware_kg_W'],cap*z['chemistry_USD2026_W'],z['chemistry_life_years'],dest,group)
    local_firm=f['local_mean_W']*m/l['process_availability']
    extra_ground=max(0,local_firm-d['old_construction_W'])
    extra_rad=extra_ground*d['old_radiator_m2']/d['old_total_W']
    add('additional_ground_receivers',extra_ground*d['power_hardware_kg_W'],extra_ground*prices['ground_receiver_USD2026_W'],20,'surface','lunar')
    add('additional_ground_radiator_controls',extra_rad*d['radiator_controls_kg_m2'],extra_rad*prices['ground_radiator_assembly_USD2026_m2'],30,'surface','lunar')
    # The existing 2 GW construction receivers/radiators are credited in extra_ground.
    total_firm=d['old_total_W']+extra_ground
    source_DC=total_firm/(p['RF_efficiency']*p['beam_capture']*p['receiver_efficiency'])*p['stations']/2/p['end_of_life_factor']
    add('lunar_beam_constellation',source_DC/p['space_system_W_kg'],source_DC*p['space_USD2026_W'],p['space_life_years'],'depot','power')
    for name,key,dest,group in [('source_solar','source_mean_W','source','source'),('depot_solar','depot_mean_W','depot','landing')]:
        W=f[key]*m/p['end_of_life_factor']/(p['source_solar_flux_fraction'] if dest=='source' else 1)
        add(name,W/p['space_system_W_kg'],W*p['space_USD2026_W'],p['space_life_years'],dest,group)
    flux=p['emissivity']*STEFAN_BOLTZMANN*p['radiator_K']**4
    for name,key,dest,group in [('source_radiators','source_heat_W','source','source'),('depot_radiators','depot_heat_W','depot','landing')]:
        area=f[key]*m/flux
        add(name,area*p['radiator_kg_m2'],area*p['radiator_kg_m2']*prices['lunar_hardware_USD2026_kg'],20,dest,group)
    for name,key,payload,dry,cycles,life,dest,group in [('lunar_landers','moon_payload_kg_year',t['lander_payload_kg'],t['lander_dry_to_payload'],t['lander_cycles_year'],t['lander_life_years'],'depot','landing'),('asteroid_tugs','nea_payload_kg_year',t['nea_payload_kg'],t['nea_dry_to_payload'],1/t['nea_cycle_years'],t['nea_life_years'],'source','source'),('Earth_transfer_tugs','earth_deep_kg_year',t['earth_payload_kg'],t['earth_dry_to_payload'],t['earth_cycles_year'],t['earth_tug_life_years'],'leo','earth_transfer')]:
        n=math.ceil(f[key]*m/(payload*cycles*t['fleet_availability']));counts[name]=n
        fraction=1-z['tug_local_mass_fraction'] if name=='asteroid_tugs' else 1
        add(name,n*payload*dry,n*payload*dry*prices['spacecraft_USD2026_kg']*fraction,life,dest,group)
        assets[name]['import_fraction']=fraction
    tank_mass=(f['nea_payload_kg_year']+f['lunar_propellant_kg_year'])*t['depot_inventory_years']*t['tank_dry_to_contents']
    add('depots_and_orbital_workshops',tank_mass,tank_mass*prices['spacecraft_USD2026_kg'],t['tank_life_years'],'depot','landing')
    return assets,dict(counts=counts,local_industry_firm_W=local_firm,total_ground_firm_W=total_firm,
        added_ground_firm_W=extra_ground,added_ground_radiator_m2=extra_rad,beam_PV_nameplate_W=source_DC,
        beam_PV_area_m2=source_DC/(SOLAR_CONSTANT*p['PV_efficiency']),process_radiator_flux_W_m2=flux,
        beam_converter_radiator_m2=source_DC*p['end_of_life_factor']*(1-p['RF_efficiency'])/flux,
        source_process_radiator_m2=f['source_heat_W']*m/flux,depot_process_radiator_m2=f['depot_heat_W']*m/flux)


def solve(s,d):
    assets={}
    for iteration in range(160):
        f=flows(s,d,assets); new,engineering=design_assets(s,d,f)
        error=max(abs(new[k]['mass_kg']-assets.get(k,{}).get('mass_kg',0))/max(new[k]['mass_kg'],1) for k in new)
        assets=new
        if error<1e-10:
            return assets,flows(s,d,assets),engineering,iteration+1,error
    raise ValueError('Industrial replacement/deployment mass balance did not converge')


def budget(s,d,assets,f):
    t=s['transport']; e=s['economics']; prices=s['prices']; H=s['schedule']['delivery_equivalent_years']
    inflation=(1+e['inflation'])**(e['automation_year']-e['base_year'])
    prop=f['earth_round_trip']['total']*t['propellant_margin']
    def earth_freight(kg,deep=True):
        return kg*((1+prop)*t['LEO_USD2026_kg']+prop*t['earth_propellant_USD2026_kg']) if deep else kg*t['LEO_USD2026_kg']
    installed={k:a['manufacture_USD2026']+earth_freight(a['mass_kg']*a.get('import_fraction',1),a['destination']!='leo') for k,a in assets.items()}
    by_group={k:sum(installed[n] for n,a in assets.items() if a['group']==k) for k in ('lunar','source','landing','earth_transfer','power')}
    hardware=sum(installed.values())
    preservation=hardware*e['source_preservation_fraction_of_supplier_capital']
    research=e['research_qualification_USD2026']
    # Earth seed water commissions a small surface/source fraction before source delivery.
    water_per_prop=s['source']['water_to_hydrogen_mass_ratio']/(1+t['mixture_O_to_H'])
    surface=sum(a['mass_kg'] for a in assets.values() if a['destination']=='surface')
    source=sum(a['mass_kg'] for a in assets.values() if a['destination']=='source')
    seed_water=water_per_prop*s['schedule']['seed_reserve_factor']*(surface*s['schedule']['seed_surface_fraction']*f['lander_round_trip']['total']+source*s['schedule']['seed_source_fraction']*(f['nea_round_trip']['mass_ratio']-1))*t['propellant_margin']
    seed_cost=earth_freight(seed_water)+seed_water*t['earth_propellant_USD2026_kg']
    capital_direct=hardware+preservation+research+seed_cost
    capital=capital_direct*(1+e['programme_contingency'])
    fixed_operation=hardware*e['supplier_upkeep_fraction_year']+e['staff']*e['loaded_staff_USD2026_year']
    renewal=sum(installed[k]/a['life_years'] for k,a in assets.items())
    initial_goods=d['pressure']*prices['pressure_USD2026_kg']+(d['earth_special']-d['pressure'])*prices['other_special_USD2026_kg']
    initial_freight=earth_freight(d['earth_special'])
    initial_restoration=(f['raw_lunar_kg_year']-flows(s,d,assets,False)['raw_lunar_kg_year'])*H*e['local_restoration_USD2026_kg_feed']
    initial_stock=initial_goods+initial_freight+initial_restoration
    # Full-equivalent operating allowance pays for owned processing, internal transport
    # and replacement during production. No internal per-kg tariff is added again.
    campaign_supplier=(H+s['schedule']['commission_years'])*(fixed_operation+renewal)
    steady=flows(s,d,assets,False)
    yearly_goods=d['yearly_pressure']*prices['pressure_USD2026_kg']+(d['yearly_earth']-d['yearly_pressure'])*prices['other_special_USD2026_kg']
    yearly_freight=earth_freight(d['yearly_earth'])
    restoration=steady['raw_lunar_kg_year']*e['local_restoration_USD2026_kg_feed']
    yearly_direct=fixed_operation+renewal+yearly_goods+yearly_freight+restoration
    yearly=yearly_direct*(1+e['programme_contingency'])
    occupied_equivalent_years=(H+s['schedule']['ramp_years']/2)/2+s['schedule']['commission_years']
    campaign_settlement_makeup=occupied_equivalent_years*(yearly_goods+yearly_freight+restoration)
    campaign=(capital_direct+initial_stock+campaign_supplier+campaign_settlement_makeup)*(1+e['programme_contingency'])
    # Level annual supplier tariffs: capital recovery at zero real interest, not a
    # second cash bill. Infrastructure imports/renewals already enter rated flows.
    shared_annual=(research+preservation+seed_cost)/H+e['staff']*e['loaded_staff_USD2026_year']
    fixed_group={g:(sum(installed[k]*(1/a['life_years']+e['supplier_upkeep_fraction_year']) for k,a in assets.items() if a['group']==g)+shared_annual*by_group[g]/hardware)*(1+e['programme_contingency']) for g in by_group}
    source_unit=fixed_group['source']/f['nea_payload_kg_year']
    landing_unit=fixed_group['landing']/f['moon_payload_kg_year']+f['depot_water_for_propellant_kg_year']/f['moon_payload_kg_year']*source_unit
    earth_unit=earth_freight(1)+fixed_group['earth_transfer']/f['earth_deep_kg_year']+landing_unit
    # All air/water cargo uses the water equivalent rate for a bounded pooled tariff.
    volatile_unit=source_unit+landing_unit
    U_pressure=d['U_pressure_kg'];U_volatile=d['U_volatile_kg'];U_freight=d['U_other_freight_kg']
    basket_U=d['pressure']*U_pressure+(d['water']+d['air'])*U_volatile+(d['earth_special']-d['pressure'])*U_freight
    basket_USD=d['pressure']*(prices['pressure_USD2026_kg']+earth_unit)+(d['water']+d['air'])*volatile_unit+(d['earth_special']-d['pressure'])*earth_unit
    out=dict(inflation_factor=inflation,installed_assets_USD2026=installed,groups_USD2041={k:v*inflation for k,v in by_group.items()},
        preservation_USD2041=preservation*inflation,research_USD2041=research*inflation,seed_water_kg=seed_water,seed_USD2041=seed_cost*inflation,
        capital_direct_USD2041=capital_direct*inflation,capital_with_contingency_USD2041=capital*inflation,
        initial_stock_manufacture_USD2041=initial_goods*inflation,initial_stock_Earth_freight_USD2041=initial_freight*inflation,
        initial_stock_restoration_USD2041=initial_restoration*inflation,
        campaign_supplier_operation_and_renewal_USD2041=campaign_supplier*inflation,
        construction_period_settlement_makeup_USD2041=campaign_settlement_makeup*inflation,
        construction_occupied_equivalent_years=occupied_equivalent_years,
        industrial_and_initial_material_campaign_USD2041=campaign*inflation,
        annual_supplier_operation_USD2041=fixed_operation*inflation,annual_supplier_renewal_USD2041=renewal*inflation,
        annual_settlement_import_manufacture_USD2041=yearly_goods*inflation,annual_settlement_import_freight_USD2041=yearly_freight*inflation,
        annual_local_restoration_USD2041=restoration*inflation,annual_supply_direct_USD2041=yearly_direct*inflation,
        annual_supply_with_contingency_USD2041=yearly*inflation,
        tariffs_USD2041_kg={'Earth_to_Moon_freight':earth_unit*inflation,'pooled_volatile_delivery':volatile_unit*inflation,'pressure_manufactured_and_delivered':(prices['pressure_USD2026_kg']+earth_unit)*inflation},
        U_calibration={'basket_U':basket_U,'basket_USD2041':basket_USD*inflation,'USD2041_per_U':basket_USD/basket_U*inflation,
          'component_USD2041_per_U':{'pressure':(prices['pressure_USD2026_kg']+earth_unit)*inflation/U_pressure,'volatiles':volatile_unit*inflation/U_volatile,'freight':earth_unit*inflation/U_freight},
          'scope':'Proposed procurement-basket calibration only. Pressure systems, air/water, and incremental freight for other site imports at inherited 60, 1 and 10 U/kg; no general U exchange rate or blanket conversion of prior settlement capital.'})
    return out,steady


def case(s,src):
    d=settlement_demand(s,src); assets,f,eng,iterations,error=solve(s,d); cost,steady=budget(s,d,assets,f)
    p=s['power'];t=s['transport'];e=s['economics'];sched=s['schedule']
    period=2*math.pi*math.sqrt((MOON_RADIUS+p['orbital_altitude_m'])**3/MOON_GM)/3600
    visible_arc=2*math.acos(MOON_RADIUS/((MOON_RADIUS+p['orbital_altitude_m'])*math.cos(math.radians(p['latitude_limit_deg']))))
    visible_min=math.floor(visible_arc/(2*math.pi/p['stations']))
    wavelength=SPEED_OF_LIGHT/p['microwave_Hz'];distance=2*MOON_RADIUS+p['orbital_altitude_m']
    spot=1.22*wavelength*distance/p['aperture_m']
    usable_storage_fraction=src['settlement_inputs']['utilities']['critical_reserve_hours']*(src['settlement_inputs']['utilities']['critical_kW_person']+src['settlement_inputs']['utilities']['food_preservation_critical_kW_person'])*1000*d['population']/d['reserve_Wh']
    local_heat=max(0,eng['local_industry_firm_W']-d['old_construction_W'])
    source_initial=f['source_ore_kg_year']*sched['delivery_equivalent_years']
    source_reservation=source_initial+steady['source_ore_kg_year']*sched['reserve_horizon_years']
    elapsed=sched['bootstrap_years']+sched['delivery_equivalent_years']+sched['ramp_years']/2+sched['commission_years']
    eng.update(orbital_period_hours=period,minimum_geometric_sources_visible=visible_min,
        sources_after_one_failure_and_one_lunar_eclipse=visible_min-2,
        far_range_m=distance,first_null_radius_m=spot,
        receiver_active_area_m2=eng['total_ground_firm_W']/p['receiver_mean_flux_W_m2'],
        full_occupied_load_reserve_hours=d['reserve_Wh']*usable_storage_fraction/d['occupied_mean_W'],
        launch_flights_year=math.ceil(f['LEO_kg_year']/t['launch_payload_kg']),
        launch_pads=math.ceil(f['LEO_kg_year']/t['launch_payload_kg']/t['launches_per_pad_year']),
        chemical_exhaust_kg_year=f['lunar_propellant_kg_year']+f['nea_propellant_kg_year']+f['earth_transfer_propellant_kg_year'],
        Earth_launcher_propellant_kg_year=f['LEO_kg_year']*t['launcher_propellant_kg_per_kg_LEO'],
        Earth_launcher_methane_kg_year=f['LEO_kg_year']*t['launcher_propellant_kg_per_kg_LEO']/(1+t['launcher_O_to_CH4']),
        Earth_gross_launcher_CO2_kg_year=f['LEO_kg_year']*t['launcher_propellant_kg_per_kg_LEO']/(1+t['launcher_O_to_CH4'])*t['methane_CO2_mass_ratio'],
        Earth_synthetic_fuel_mean_W=f['LEO_kg_year']*t['launcher_propellant_kg_per_kg_LEO']/(1+t['launcher_O_to_CH4'])*t['synthetic_methane_kWh_kg']*1000/HOURS)
    return dict(demand=d,assets=assets,peak_flows=f,steady_flows=steady,engineering=eng,cost=cost,
        convergence={'iterations':iterations,'relative_mass_residual':error},
        schedule={'start_year':e['automation_year'],'first_source_capacity_year':e['automation_year']+sched['bootstrap_years'],
          'full_rate_year':e['automation_year']+sched['bootstrap_years']+sched['ramp_years'],
          'conditional_completion_year':e['automation_year']+elapsed,'elapsed_years':elapsed},
        resource_reservation={'source_campaign_upper_kg':source_initial,'source_campaign_plus_200year_kg':source_reservation,
          'equivalent_spherical_body_diameter_m':(6*source_reservation/(math.pi*s['source']['source_density_kg_m3']))**(1/3),
          'rule':'An equivalent volume, not a target body, mineable fraction or permission. Campaign envelope includes full steady replenishment throughout delivery; no permitted extractable fraction is assumed.'})


def electric_tug(s, reference):
    """Close a single water-electric tug's power/heat/dry-mass requirement."""
    q=s['electric_tug_comparison']; p=s['power']; cargo=s['transport']['nea_payload_kg']
    back=reference['peak_flows']['nea_return_hardware_kg_year']/reference['peak_flows']['nea_payload_kg_year']
    ve=q['Isp_s']*STANDARD_GRAVITY; results=[]
    for eta in q['efficiencies']:
        dry=q['structural_dry_to_payload']*cargo
        for iteration in range(100):
            prop=cargo*round_trip(q['dv_oneway_m_s'],q['Isp_s'],dry/cargo,back)['total']
            W=prop*ve**2/(2*eta*q['powered_days_per_roundtrip']*JULIAN_DAY)
            array_W=W/p['end_of_life_factor']
            rad=W*(1-eta)/(p['emissivity']*STEFAN_BOLTZMANN*p['radiator_K']**4)
            new=q['structural_dry_to_payload']*cargo+array_W/p['space_system_W_kg']+rad*p['radiator_kg_m2']+W*q['power_processing_kg_W']
            if abs(new-dry)<1e-5:break
            dry=new
        else:raise ValueError('Electric tug mass did not converge')
        results.append(dict(efficiency=eta,cargo_kg=cargo,dry_kg=new,propellant_water_kg=prop,power_W=W,PV_nameplate_W=array_W,radiator_m2=rad,average_thrust_N=prop*ve/(q['powered_days_per_roundtrip']*JULIAN_DAY),powered_days=q['powered_days_per_roundtrip']))
    return results


def run(spec=None):
    s,src=read_inputs(spec)
    sensitivities={}
    for name,launch,water,years in [('launch_50',50,1,20),('water_makeup_quarter',200,.25,20),('stock_delivery_10_years',200,1,10),('stock_delivery_40_years',200,1,40)]:
        q=deepcopy(s); q['transport']['LEO_USD2026_kg']=launch; q['schedule']['delivery_equivalent_years']=years
        data=deepcopy(src); data['settlement']['cases']['reference']['resources']['operating_import_kg_year']['water']*=water
        c=case(q,data);sensitivities[name]={'water_makeup_factor':water,'launch_USD2026_kg':launch,'delivery_equivalent_years':years,'capital_USD2041':c['cost']['capital_with_contingency_USD2041'],'campaign_USD2041':c['cost']['industrial_and_initial_material_campaign_USD2041'],'annual_USD2041':c['cost']['annual_supply_with_contingency_USD2041'],'source_ore_kg_year':c['peak_flows']['source_ore_kg_year'],'completion_year':c['schedule']['conditional_completion_year']}
    cases={name:case(with_case(s,name),src) for name in s['case_overrides']}
    return {'schema':SCHEMA,'date':s['date'],'evidence':'Scenario engineering and procurement arithmetic; not an empirical forecast, flight design, reserve certification or settlement site approval.',
        'inputs':s['inputs'],'producer':{'path':'research/studies/initial_industrial_system/run.py','sha256':sha(__file__),'configuration':{'path':'research/studies/initial_industrial_system/industrial_inputs.json','sha256':sha(HERE/'industrial_inputs.json')},'sources':{'path':'research/studies/initial_industrial_system/industrial_sources.json','sha256':sha(HERE/'industrial_sources.json')}},
        'cases':cases,'sensitivities':sensitivities,'electric_tug_comparison':electric_tug(s,cases['reference']),
        'rules':['No previous settlement U or physical product is modified.','Credit the 150 million kg embedded construction plant once; retain the existing ground utility inventory and add only industrial capacity beyond its 2 GW allocation.','Prices are 2026-real scenario inputs; one 15-year inflation multiplication gives 2041 dollars. All later spending is at flat 2041 prices.','Supplier capital, operating costs and renewal are paid in the owner budget; diagnostic internal tariffs are never added to that budget.','The material campaign excludes civil installation and commissioning labour within settlements, passenger transport, settlement civic, food-growing and household operations during phased occupation, finance, planetary terraforming, and unqualified source/site losses.']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'research/runs/initial_industrial_system.json');a=p.parse_args()
    result=run();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    for name,c in result['cases'].items():
        q=c['cost']; print(name,'capital T$',round(q['capital_with_contingency_USD2041']/1e12,3),'campaign T$',round(q['industrial_and_initial_material_campaign_USD2041']/1e12,3),'annual B$',round(q['annual_supply_with_contingency_USD2041']/1e9,2),'iterations',c['convergence']['iterations'])
if __name__=='__main__':main()
