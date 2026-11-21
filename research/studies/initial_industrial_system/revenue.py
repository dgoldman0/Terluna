"""Consolidated external receipts, Earth returns and resource-horizon scenarios.

Internal purchases and resident income are not additional external receipts.
The annual net measures are whole-system external cash balances, not GDP or
corporate profit. Prospective markets are assumptions, never booked contracts.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
from . import run as industry
from . import transport
from shared.constants import EARTH_GM, EARTH_RADIUS, JULIAN_DAY, JULIAN_YEAR_DAYS, STANDARD_GRAVITY, STEFAN_BOLTZMANN
HERE,ROOT=industry.HERE,industry.ROOT
SCHEMA='terluna.research.initial-system-revenue/1'


def read_inputs(spec=None):
    q=deepcopy(spec) if spec is not None else json.loads((HERE/'revenue_inputs.json').read_text())
    src={}
    for k,b in q['inputs'].items():
        path=ROOT/b['path']
        if industry.sha(path)!=b['sha256']:raise ValueError(f'Input hash mismatch: {k}')
        src[k]=json.loads(path.read_text())
    return q,src


def supply_cases(src):
    s=src['industrial_configuration'];t=src['transport'];tc=src['transport_configuration'];out={}
    for name,c in t['cases'].items():out[name]=(deepcopy(s),c,transport.scenario(tc,name))
    for name in ['split_chemical','space_chemical','electric_earth_leg']:
        out['local_structures_'+name]=(industry.with_case(s,'local_fleet_fabrication'),t['source_structure_localization'][name],transport.scenario(tc,name))
    ss=deepcopy(s);ss['power']['space_system_W_kg']=200
    qq=transport.scenario(tc,'electric_both_legs');qq['electric']['total_efficiency']=.5
    out['electric_development']=(ss,t['sensitivities']['electric_50pct_200Wkg']['case'],qq)
    out['favourable_earth_chemical']=(industry.with_case(s,'lower'),src['industrial']['cases']['lower'],transport.scenario(tc,'earth_chemical'))
    return out


def trip_burns(cargo,back,dry,dv,isp,margin,mode,retention=1):
    """Burns per trip; tank constraints are checked against constructed vehicles."""
    k=math.expm1(dv/(isp*STANDARD_GRAVITY))
    if mode=='earth_chemical':
        b=(dry+back)*k;a=(cargo+dry+b)*k
        a*=margin;b*=margin
        return {'out':a,'back':b,'Earth_fuel':a+b,'space_supply':0,'tank':a+b,'burned':a+b}
    if mode=='split_chemical':
        a=(cargo+dry)*k*margin;b=(back+dry)*k*margin
        return {'out':a,'back':b,'Earth_fuel':a,'space_supply':b/retention,'tank':max(a,b),'burned':a+b}
    a=(cargo+dry)*k*margin
    carried=a/retention;b=(dry+carried+back)*k*margin
    return {'out':a,'back':b,'Earth_fuel':0,'space_supply':(carried+b)/retention,'tank':carried+b,'burned':a+b}


def ferry(cargo_year,back_year,*,payload,dry,dv,isp,margin,mode,retention,tank_limit,burn_limit=None):
    """Reduce loading, if needed, to keep the original tanks and vehicle.

Annual rates may be fractional averages. Extra trips are tested against the
fleet's rated capacity separately, preserving its 25 percent spare margin.
"""
    n=max(cargo_year/payload,back_year/payload,1e-12)
    def at(trips):return trip_burns(cargo_year/trips,back_year/trips,dry,dv,isp,margin,mode,retention)
    def fits(trips):
        b=at(trips)
        return b['tank']<=tank_limit*(1+1e-12) and (burn_limit is None or b['burned']<=burn_limit*(1+1e-12))
    low=n
    for _ in range(100):
        if fits(n):break
        n*=1.15
    else:raise ValueError('Return cargo cannot fit this vehicle')
    if n>low:
        high=n
        for _ in range(70):
            mid=(low+high)/2
            if fits(mid):high=mid
            else:low=mid
        n=high
    b=at(n)
    return {'trips_year':n,'outbound_kg_trip':cargo_year/n,'return_goods_kg_trip':back_year/n,
            'tank_kg_trip':b['tank'],'tank_limit_kg':tank_limit,
            **{k+'_kg_year':n*b[k] for k in ['Earth_fuel','space_supply','burned']}}


def commerce(q,s,c,m,disabled=()):
    """Explicit extra facilities; existing resident facilities retain their scope."""
    t,p=s['transport'],s['power'];pop=q['boundary']['population'];assets={};receipts={};operating={}
    civil=0;initial_bulk=0;initial_metal=0;initial_water=0;initial_air=0
    op_surface=0;op_leo=0;mean_surface_W=0;staff=0
    prices=s['prices'];f=q['manufacturing'];v=q['visitors'];comp=q['computing'];sv=q['services'];rr=q['return']
    def asset(name,mass,price,life,destination,sortie_freight=False):
        assets[name]={'mass_kg':mass,'manufacture_USD2026':price,'life_years':life,'destination':destination,'freight_in_sorties':sortie_freight}
    def equipment(name,mass,life,price=500,destination='surface'):
        asset(name,mass,mass*price,life,destination)
    nameplate=pop*f['nameplate_finished_kg_per_resident_year'] if 'manufacturing' not in disabled else 0
    delivered_fraction=(1-rr['loss_probability_flight'])*f['delivered_quality_yield']
    made=min(nameplate*f['saleable_yield']*m['factory_utilization'],m['factory_customer_budget_USD2026_year']/m['factory_USD2026_kg']/delivered_fraction)
    delivered=made*delivered_fraction;packed=made*(1+f['transport_packaging_mass_fraction'])
    material=made/f['saleable_yield'];factory_import=material*f['Earth_input_mass_fraction']
    local_export_metal=material-factory_import+made*f['transport_packaging_mass_fraction']
    receipts['manufactured_goods_at_Earth']=delivered*m['factory_USD2026_kg']
    equipment('export_factory',nameplate*f['equipment_kg_per_nameplate_kg_year'],f['equipment_life_years'],prices['lunar_hardware_USD2026_kg'])
    civil+=nameplate*f['civil_and_commissioning_USD2026_per_nameplate_kg_year']
    factory_floor=nameplate*f['floor_m2_per_nameplate_kg_year'];initial_bulk+=factory_floor*f['bulk_kg_m2'];initial_metal+=factory_floor*f['metal_kg_m2'];initial_air+=factory_floor*f['air_kg_m2']
    op_surface+=factory_import;operating['factory_imports']=factory_import*prices['other_special_USD2026_kg']
    staff+=f['staff'] if nameplate else 0
    mean_surface_W+=made*f['finish_kWh_per_finished_kg']*1000/industry.HOURS
    # Capsule shields are single-use operating inputs; durable capsule renewal
    # is separate. Relaunch includes both reused and newly replaced buses once.
    flights=packed/rr['capsule_payload_kg'];capsules=math.ceil(flights*s['economics']['capacity_margin']/(rr['cycles_year']*t['fleet_availability']))
    bus=rr['capsule_payload_kg']*rr['capsule_dry_to_payload']*(1-rr['shield_fraction_of_dry'])
    shield=rr['capsule_payload_kg']*rr['capsule_dry_to_payload']*rr['shield_fraction_of_dry']
    loss=rr['loss_probability_flight'];life_flights=(1-(1-loss)**rr['durable_life_flights'])/loss if loss else rr['durable_life_flights']
    life=min(rr['durable_calendar_life_years'],life_flights*capsules/flights) if flights else rr['durable_calendar_life_years']
    equipment('return_capsule_buses',capsules*bus,life,prices['spacecraft_USD2026_kg'],'leo')
    assets['return_capsule_buses']['freight_in_sorties']=True
    prop_per_flight=(rr['capsule_payload_kg']+bus+shield)*math.expm1(rr['deorbit_dv_m_s']/(rr['deorbit_Isp_s']*STANDARD_GRAVITY))
    op_leo+=flights*(bus+shield+prop_per_flight)
    operating['entry_shields']=flights*shield*prices['spacecraft_USD2026_kg']
    operating['entry_propellant']=flights*prop_per_flight*t['earth_propellant_USD2026_kg']
    operating['Earth_recovery_and_delivery']=packed*rr['Earth_ground_service_USD2026_per_payload_kg']
    initial_shields=capsules*shield*prices['spacecraft_USD2026_kg']
    beds=math.ceil(m['visitors_year']*v['stay_nights']/(JULIAN_YEAR_DAYS*v['occupancy']))
    guest_years=m['visitors_year']*v['stay_nights']/JULIAN_YEAR_DAYS
    receipts['destination_services']=m['visitors_year']*v['stay_nights']*m['visitor_USD2026_night']
    operating['visitor_Earth_services']=receipts['destination_services']*v['Earth_service_share_of_receipts']
    equipment('visitor_life_support_and_food',beds*v['surface_hardware_kg_bed'],v['equipment_life_years'],prices['lunar_hardware_USD2026_kg'])
    civil+=beds*v['civil_USD2026_bed'];initial_water+=beds*v['initial_water_kg_bed'];initial_air+=beds*v['initial_air_kg_bed']
    initial_bulk+=beds*v['local_bulk_kg_bed'];initial_metal+=beds*v['local_metal_kg_bed'];mean_surface_W+=beds*v['mean_W_bed'];staff+=beds*v['staff_per_bed']
    water=c['demand']['yearly_water']*guest_years/pop;air=c['demand']['yearly_air']*guest_years/pop
    settlement_op_earth=q['_operating_settlement_Earth_kg']*guest_years/pop
    op_surface+=settlement_op_earth;operating['visitor_consumables']=settlement_op_earth*prices['other_special_USD2026_kg']
    W=m['compute_W'] if 'computing' not in disabled else 0
    receipts['batch_computing']=W*comp['billable_duty']*m['compute_USD2026_per_W_billable_year']
    asset('batch_IT',W*comp['IT_kg_W']*comp['space_packaging_factor'],W*comp['IT_USD2026_W'],comp['IT_life_years'],'leo')
    asset('batch_communications',W*comp['communications_kg_W'],W*comp['communications_USD2026_W'],comp['IT_life_years'],'leo')
    array=W/p['end_of_life_factor'];rad=W/(p['emissivity']*STEFAN_BOLTZMANN*p['radiator_K']**4)
    asset('batch_solar',array/p['space_system_W_kg'],array*p['space_USD2026_W'],p['space_life_years'],'leo')
    equipment('batch_radiators',rad*p['radiator_kg_m2'],20,prices['lunar_hardware_USD2026_kg'],'leo')
    operating['computing_Earth_services']=receipts['batch_computing']*comp['Earth_service_share_of_receipts'];staff+=W/1e6*comp['staff_per_MW']
    receipts['research_services']=m['research_receipts_USD2026_year']
    receipts['media_and_licensing']=m['media_receipts_USD2026_year']
    receipts['remote_professional_services']=m['remote_staff']*m['remote_USD2026_staff_year']
    for sector,count,kg,life,price,share,receipt in [
      ('research',m['research_staff'],sv['research_hardware_kg_staff'],sv['research_hardware_life_years'],prices['lunar_hardware_USD2026_kg'],sv['research_Earth_service_share'],'research_services'),
      ('media',m['media_staff'],sv['media_hardware_kg_staff'],sv['media_hardware_life_years'],prices['spacecraft_USD2026_kg'],sv['media_Earth_service_share'],'media_and_licensing'),
      ('remote',m['remote_staff'],sv['remote_hardware_kg_staff'],sv['remote_hardware_life_years'],prices['spacecraft_USD2026_kg'],sv['remote_Earth_service_share'],'remote_professional_services')]:
        equipment(sector+'_equipment',count*kg,life,price)
        operating[sector+'_Earth_services']=receipts[receipt]*share;staff+=count
    civil+=m['research_staff']*sv['research_civil_USD2026_staff']
    lab_floor=m['research_staff']*sv['research_floor_m2_staff'];initial_bulk+=lab_floor*sv['bulk_kg_m2'];initial_metal+=lab_floor*sv['metal_kg_m2'];initial_air+=lab_floor*sv['air_kg_m2']
    reagents=m['research_staff']*sv['research_consumables_kg_staff_year'];op_surface+=reagents;operating['research_consumables']=reagents*prices['other_special_USD2026_kg']
    mean_surface_W+=m['research_staff']*sv['research_W_staff']+(m['media_staff']+m['remote_staff'])*sv['office_W_staff']
    # Independent Earth mission procurement competes with Earth-provided fuel.
    original=industry.round_trip(t['earth_dv_oneway_m_s'],t['Isp_s'],t['earth_dry_to_payload'])['total']*t['propellant_margin']
    Earth_fuel_alternative=(1+original)*(t['LEO_USD2026_kg']+t['earth_propellant_USD2026_kg'])
    fuel_price=min(m['propellant_USD2026_kg'],q['propellant_customers']['Earth_alternative_price_fraction_ceiling']*Earth_fuel_alternative)
    fuel_qty=min(m['propellant_customer_kg_year'],max(0,m['mission_customer_budget_USD2026_year']-receipts['research_services'])/fuel_price)
    receipts['Earth_mission_propellant']=fuel_qty*fuel_price
    renew_mass={dest:sum(a['mass_kg']/a['life_years'] for a in assets.values() if a['destination']==dest and not a['freight_in_sorties']) for dest in ['surface','leo']}
    initial_mass={dest:sum(a['mass_kg'] for a in assets.values() if a['destination']==dest) for dest in ['surface','leo']}
    manufacture=sum(a['manufacture_USD2026'] for a in assets.values())
    renewal_manufacture=sum(a['manufacture_USD2026']/a['life_years'] for a in assets.values())
    return dict(receipts_USD2026=receipts,operating_purchases_USD2026=operating,assets=assets,
                civil_and_commissioning_USD2026=civil,initial_shields_USD2026=initial_shields,initial_shields_kg=capsules*shield,
                manufacture_USD2026=manufacture,renewal_manufacture_USD2026=renewal_manufacture,
                operating_surface_kg=op_surface,operating_LEO_kg=op_leo,renewal_mass_kg=renew_mass,initial_mass_kg=initial_mass,
                goods_made_kg=made,goods_delivered_saleable_kg=delivered,return_packed_kg=packed,
                factory_input_kg=material,factory_rejects_kg=material-made,local_export_metal_kg=local_export_metal,
                initial_bulk_kg=initial_bulk,initial_metal_kg=initial_metal,initial_water_kg=initial_water,initial_air_kg=initial_air,
                guest_years=guest_years,visitor_beds=beds,guest_water_kg_year=water,guest_air_kg_year=air,
                added_surface_mean_W=mean_surface_W,compute_peak_W=W,compute_billable_duty=comp['billable_duty'],
                propellant_customer_kg=fuel_qty,propellant_sale_USD2026_kg=fuel_price,Earth_fuel_alternative_USD2026_kg=Earth_fuel_alternative,
                return_capsules=capsules,return_flights_year=flights,expected_capsule_life_flights=life_flights,
                return_entry_kinetic_W=flights*(rr['capsule_payload_kg']+bus+shield)*EARTH_GM/(EARTH_RADIUS+rr['LEO_altitude_m'])/2/(JULIAN_YEAR_DAYS*JULIAN_DAY),
                additional_resident_jobs=staff)


def logistics(s,c,tc,b,renewal=True):
    """Operate retained fleets with explicit return goods; do not sell reserve."""
    old=c['steady_flows'];peak=c['peak_flows'];t=s['transport'];z=s['source'];l=s['lunar'];e=s['economics'];p=s['power']
    extra_surface=b['operating_surface_kg']+(b['renewal_mass_kg']['surface'] if renewal else 0)
    extra_leo=b['operating_LEO_kg']+(b['renewal_mass_kg']['leo'] if renewal else 0)
    cargo=old['earth_deep_kg_year']+extra_surface
    mode=tc['earth_route'];back=b['return_packed_kg']
    if mode in ('earth_chemical','split_chemical'):
        dry=t['earth_payload_kg']*t['earth_dry_to_payload'];dv=t['earth_dv_oneway_m_s'];isp=t['Isp_s']
        tank=t['earth_payload_kg']*industry.round_trip(dv,isp,t['earth_dry_to_payload'])['total']*t['propellant_margin']
        burn_limit=None;cycle=JULIAN_YEAR_DAYS/t['earth_cycles_year'];count=c['engineering']['counts']['Earth_transfer_tugs'];ret=tc['chemical']['handling_retention']
    else:
        v=peak['Earth_vehicle'];dry=v['dry_kg'];dv=v['delta_v_oneway_m_s'];isp=v['Isp_s'];tank=v['burns']['peak_tank_contents_kg']
        burn_limit=v['burns']['ascent_kg']+v['burns']['descent_kg'] if mode=='electric' else None
        cycle=v['cycle_days'];count=c['engineering']['counts']['Earth_route_ships'];ret=tc['electric' if mode=='electric' else 'chemical']['handling_retention']
    earth=ferry(cargo,back,payload=t['earth_payload_kg'],dry=dry,dv=dv,isp=isp,margin=t['propellant_margin'],mode=mode,retention=ret,tank_limit=tank,burn_limit=burn_limit)
    moon_cargo=old['moon_payload_kg_year']+extra_surface+b['guest_water_kg_year']+b['guest_air_kg_year']
    P=t['lander_payload_kg'];D=P*t['lander_dry_to_payload'];dv=t['lander_dv_oneway_m_s']
    tank=P*industry.round_trip(dv,t['Isp_s'],D/P)['total']*t['propellant_margin']
    moon=ferry(moon_cargo,back,payload=P,dry=D,dv=dv,isp=t['Isp_s'],margin=t['propellant_margin'],mode='earth_chemical',retention=1,tank_limit=tank)
    # Record the actual origin after using the full-return-fuel burn formula.
    moon['Earth_fuel_kg_year']=0;moon['space_supply_kg_year']=moon['burned_kg_year']
    land_prop=moon['burned_kg_year'];oxygen_extra=b['guest_air_kg_year']*(1-z['air_nitrogen_mass_fraction'])
    nitrogen=old['source_nitrogen_kg_year']+b['guest_air_kg_year']*z['air_nitrogen_mass_fraction']
    old_earth_space=old.get('Earth_route_source_water_kg_year',0)
    delivered=old['nea_payload_kg_year']+b['guest_water_kg_year']+(nitrogen-old['source_nitrogen_kg_year'])+(land_prop-old['lunar_propellant_kg_year'])+oxygen_extra+(earth['space_supply_kg_year']-old_earth_space)+b['propellant_customer_kg']
    if tc['source_route']=='electric':
        source_v=peak['Source_vehicle'];source_dry=source_v['dry_kg'];source_dv=source_v['delta_v_oneway_m_s'];source_isp=source_v['Isp_s'];source_cycle=source_v['cycle_days']
        source_count=c['engineering']['counts']['Source_route_ships'];source_retention=tc['electric']['handling_retention']
        source_limit=source_v['burns']['peak_tank_contents_kg']
    else:
        source_dry=t['nea_payload_kg']*t['nea_dry_to_payload'];source_dv=t['nea_dv_oneway_m_s'];source_isp=t['Isp_s'];source_cycle=t['nea_cycle_years']*JULIAN_YEAR_DAYS
        source_count=c['engineering']['counts']['asteroid_tugs'];source_retention=1
        source_limit=t['nea_payload_kg']*peak['nea_round_trip']['total']*t['propellant_margin']
    source=ferry(delivered,old['nea_return_hardware_kg_year'],payload=t['nea_payload_kg'],dry=source_dry,dv=source_dv,isp=source_isp,margin=t['propellant_margin'],mode='earth_chemical',retention=1,tank_limit=source_limit)
    source['Earth_fuel_kg_year']=0;source['space_supply_kg_year']=source['burned_kg_year']/source_retention
    source_water=(delivered-nitrogen)+source['space_supply_kg_year']
    ore=max(source_water/(z['water_grade']*z['water_recovery']),nitrogen/(z['nitrogen_grade']*z['nitrogen_recovery']),old['source_ore_metal_bound_kg_year'])
    process_kWh=ore*z['ore_kWh_kg']+nitrogen*z['nitrogen_kWh_kg']+old['source_metal_kg_year']*z['metal_refining_and_fabrication_kWh_kg']
    source_chemical_kWh=0 if tc['source_route']=='electric' else source['burned_kg_year']/z['water_to_hydrogen_mass_ratio']*z['electrolysis_liquefaction_kWh_kg_H2']
    source_W=(process_kWh+source_chemical_kWh)*1000/industry.HOURS
    heat_fraction=1-z['hydrogen_HHV_kWh_kg']/z['electrolysis_liquefaction_kWh_kg_H2']
    source_heat=(process_kWh+source_chemical_kWh*heat_fraction)*1000/industry.HOURS
    depot_extra=land_prop-old['lunar_propellant_kg_year']+oxygen_extra+b['propellant_customer_kg']
    if mode in ('space_chemical','split_chemical'):depot_extra+=earth['space_supply_kg_year']-old_earth_space
    depot_extra_W=depot_extra/z['water_to_hydrogen_mass_ratio']*z['electrolysis_liquefaction_kWh_kg_H2']*1000/industry.HOURS
    depot_W=old['depot_mean_W']+depot_extra_W;depot_heat=old['depot_heat_W']+depot_extra_W*heat_fraction
    metal_renew=b['initial_metal_kg']/l['metal_life_years'];bulk_renew=b['initial_bulk_kg']/l['bulk_life_years']
    virgin_metal=b.get('commissioning_metal_kg_year',0);virgin_bulk=b.get('commissioning_bulk_kg_year',0)
    metal=b['local_export_metal_kg']+metal_renew+virgin_metal;bulk=bulk_renew+virgin_bulk
    raw_extra=(b['local_export_metal_kg']+virgin_metal+metal_renew*(1-l['recovery']))/l['metal_yield']+(virgin_bulk+bulk_renew*(1-l['recovery']))/l['bulk_yield']
    local_W=old['local_mean_W']+(metal*l['metal_kWh_kg']+bulk*l['bulk_kWh_kg']+raw_extra*l['handling_kWh_kg'])*1000/industry.HOURS+b['added_surface_mean_W']
    def trip_capacity(n,days):return n*JULIAN_YEAR_DAYS/days*t['fleet_availability']/e['capacity_margin']
    ratios={
      'Earth_ferry_trips':earth['trips_year']/trip_capacity(count,cycle),
      'lunar_ferry_trips':moon['trips_year']/trip_capacity(c['engineering']['counts']['lunar_landers'],JULIAN_YEAR_DAYS/t['lander_cycles_year']),
      'source_ferry_trips':source['trips_year']/trip_capacity(source_count,source_cycle),
      'source_feed':ore/peak['source_ore_kg_year'],'source_process_power':source_W/peak['source_mean_W'],
      'source_heat_rejection':source_heat/peak['source_heat_W'],'depot_power':depot_W/peak['depot_mean_W'],
      'depot_heat_rejection':depot_heat/peak['depot_heat_W'],'lunar_industrial_power':local_W/peak['local_mean_W'],
      'lunar_bulk_production':(old['bulk_kg_year']+bulk)/peak['bulk_kg_year'],
      'lunar_metal_production':(old['metal_kg_year']+metal)/peak['metal_kg_year'],
      'lunar_handling':(old['handling_kg_year']+raw_extra)/peak['handling_kg_year']}
    if max(ratios.values())>1+1e-8:raise ValueError('Commercial activity exceeds retained rated capacity: '+str(ratios))
    old_Earth_prop=old.get('earth_launch_propellant_kg_year',old['earth_transfer_propellant_kg_year'])
    LEO=old['LEO_kg_year']+extra_surface+extra_leo+earth['Earth_fuel_kg_year']-old_Earth_prop
    return dict(Earth=earth,Moon=moon,source=source,LEO_kg_year=LEO,
                Earth_launches_year=math.ceil(LEO/t['launch_payload_kg']),
                source_water_kg_year=source_water,source_nitrogen_kg_year=nitrogen,source_ore_kg_year=ore,
                source_retained_tails_kg_year=ore-source_water-nitrogen-old['source_virgin_metal_kg_year'],
                raw_lunar_kg_year=old['raw_lunar_kg_year']+raw_extra,local_industry_mean_W=local_W,
                source_process_mean_W=source_W,depot_mean_W=depot_W,capacity_use_fraction_of_rated=ratios,
                source_virgin_metal_kg_year=old['source_virgin_metal_kg_year'],
                Earth_cargo_kg_year=cargo,lunar_cargo_kg_year=moon_cargo,source_cargo_kg_year=delivered,
                local_added_bulk_kg_year=bulk,local_added_metal_kg_year=metal,
                rule='All capacity ratios preserve the original separate spare margin; availability and tank limits are also applied.')


def accounts(q,s,c,tc,m,disabled=()):
    b=commerce(q,s,c,m,disabled);full=logistics(s,c,tc,b);op=logistics(s,c,tc,b,False)
    I=c['cost']['inflation_factor'];t=s['transport'];e=s['economics'];p=s['prices'];fin=q['finance'];old=c['steady_flows'];cost=c['cost']
    old_Earth_prop=old.get('earth_launch_propellant_kg_year',old['earth_transfer_propellant_kg_year'])
    prop_per_deep=old_Earth_prop/old['earth_deep_kg_year']
    freight=(1+prop_per_deep)*t['LEO_USD2026_kg']+prop_per_deep*t['earth_propellant_USD2026_kg']
    settlement_opex=q['_operating_settlement_Earth_kg']*(p['other_special_USD2026_kg']+freight)
    settlement_renewal=(q['_renewal_settlement_pressure_kg']*p['pressure_USD2026_kg']+(q['_renewal_settlement_Earth_kg']-q['_renewal_settlement_pressure_kg'])*p['other_special_USD2026_kg']+q['_renewal_settlement_Earth_kg']*freight)
    staff_total=e['staff']*e['loaded_staff_USD2026_year']
    upkeep_total=cost['annual_supplier_operation_USD2041']/I-staff_total
    baseline_external_services=upkeep_total*fin['maintenance_services_Earth_fraction']+staff_total*fin['staff_Earth_fraction']
    supplier_renewal=cost['annual_supplier_renewal_USD2041']/I
    ground=c['assets']['additional_ground_radiator_controls']
    supplier_renewal-=ground['manufacture_USD2026']/ground['life_years']*(1-fin['ground_radiator_assembly_Earth_fraction'])
    def added_freight(f):return (f['LEO_kg_year']-old['LEO_kg_year'])*t['LEO_USD2026_kg']+(f['Earth']['Earth_fuel_kg_year']-old_Earth_prop)*t['earth_propellant_USD2026_kg']
    extra_freight_op=added_freight(op);extra_freight_renew=added_freight(full)-extra_freight_op
    # Add the commercial stock over five years while allowing the full steady
    # portfolio to operate. This is a conservative capacity screen, not a ramp
    # dispatch schedule. Existing annual plant costs continue throughout.
    years=q['boundary']['commercial_ramp_years'];ramp=deepcopy(b)
    ramp['operating_surface_kg']+=b['initial_mass_kg']['surface']/years
    ramp['operating_LEO_kg']+=(b['initial_mass_kg']['leo']+b['initial_shields_kg'])/years
    ramp['guest_water_kg_year']+=b['initial_water_kg']/years
    ramp['guest_air_kg_year']+=b['initial_air_kg']/years
    ramp['commissioning_bulk_kg_year']=b['initial_bulk_kg']/years
    ramp['commissioning_metal_kg_year']=b['initial_metal_kg']/years
    commissioning=logistics(s,c,tc,ramp)
    initial_freight=(added_freight(commissioning)-added_freight(full))*years
    initial_resource_provision=((commissioning['source_ore_kg_year']-full['source_ore_kg_year'])*fin['source_transition_USD2026_per_kg_feed']+(commissioning['raw_lunar_kg_year']-full['raw_lunar_kg_year'])*e['local_restoration_USD2026_kg_feed'])*years
    commercial_capital_direct=b['manufacture_USD2026']+initial_freight+b['civil_and_commissioning_USD2026']+b['initial_shields_USD2026']
    maintenance=commercial_capital_direct*fin['annual_maintenance_fraction_commerce_capital']*fin['commerce_maintenance_Earth_fraction']
    op_components={'baseline_Earth_services_and_staff':baseline_external_services,'baseline_settlement_recurring_imports':settlement_opex,
                   'commercial_recurring_imports_and_Earth_services':sum(b['operating_purchases_USD2026'].values()),'additional_freight_for_operations':extra_freight_op,'commercial_Earth_maintenance_services':maintenance}
    renewal_components={'baseline_supplier_replacement':supplier_renewal,'baseline_settlement_replacement_imports':settlement_renewal,'commercial_equipment_manufacture':b['renewal_manufacture_USD2026'],'additional_freight_for_replacement':extra_freight_renew}
    gross=sum(b['receipts_USD2026'].values());opex=sum(op_components.values());renewal=sum(renewal_components.values())
    resource=full['source_ore_kg_year']*fin['source_transition_USD2026_per_kg_feed']+full['raw_lunar_kg_year']*e['local_restoration_USD2026_kg_feed']
    operating=gross-opex;after_renewal=operating-renewal;net=after_renewal-resource
    extra_capital=commercial_capital_direct*(1+e['programme_contingency'])*I
    capital=cost['capital_with_contingency_USD2041']+extra_capital
    campaign=cost['industrial_and_initial_material_campaign_USD2041']+extra_capital
    return {'gross_receipts_USD2041':gross*I,'receipts_by_activity_USD2041':{k:v*I for k,v in b['receipts_USD2026'].items()},
            'external_operating_cost_USD2041':opex*I,'operating_cost_components_USD2041':{k:v*I for k,v in op_components.items()},
            'operating_cash_balance_USD2041':operating*I,'external_renewal_USD2041':renewal*I,'renewal_components_USD2041':{k:v*I for k,v in renewal_components.items()},
            'cash_after_external_renewal_USD2041':after_renewal*I,'resource_transition_provision_USD2041':resource*I,
            'balance_after_renewal_and_resource_provision_USD2041':net*I,'debt_service_USD2041':None,
            'annual_procurement_contingency_stress_USD2041':(opex+renewal)*fin['import_service_reserve_fraction']*I,
            'commercial_capital_USD2041':extra_capital,'industrial_and_commercial_establishment_USD2041':capital,
            'industrial_material_campaign_plus_commercial_capital_USD2041':campaign,
            'establishment_over_gross_receipts_years':capital/(gross*I) if gross else None,
            'campaign_cash_payback_years_if_constant_positive_net':campaign/(net*I) if net>0 else None,
            'renewal_affordability':{'available_after_operation_and_resource_USD2041':(operating-resource)*I,
               'fraction_of_external_renewal_affordable':max(0,operating-resource)/renewal if renewal else None,
               'renewal_reduction_alone_can_close_gap':operating>=resource,
               'rule':'Static diagnostic with sales and other costs fixed. It is not an engineered localization, retirement or refurbishment scenario.'},
            'baseline_upkeep_USD2026':upkeep_total,'baseline_staff_USD2026':staff_total,
            'commercial':b,'physical':full,
            'commercial_commissioning':{'years':years,'annual_combined_flows':commissioning,
               'initial_freight_USD2041':initial_freight*I,
               'initial_resource_provision_USD2041':initial_resource_provision*I,
               'capital_and_initial_resource_provision_USD2041':extra_capital+initial_resource_provision*I,
               'additional_source_feed_kg':(commissioning['source_ore_kg_year']-full['source_ore_kg_year'])*years,
               'additional_source_water_kg':(commissioning['source_water_kg_year']-full['source_water_kg_year'])*years,
               'rule':'Initial solids, air, water, equipment and shield inventory are additional. Full steady commercial operations plus five-year stock placement fit retained rated capacity. Baseline operating and renewal bills continue during these years; they are not another capital invoice. This is not a timed dispatch or transition-finance solution.'},
            'resident_export_and_industrial_jobs':b['additional_resident_jobs']+e['staff']*(1-fin['staff_Earth_fraction']),
            'net_scope':'Whole-system external operating cash balance, then external replacement and an internally reserved resource provision. No GDP, full corporate profit, tax/dividend policy or debt-service claim.'}


def portfolio(q,s,c,tc,m):
    disabled=[];rejected=[]
    potential=accounts(q,s,c,tc,m)
    current=potential
    for _ in range(3):
        changed=False
        for activity in ['manufacturing','computing']:
            if activity in disabled:continue
            alternative=accounts(q,s,c,tc,m,disabled+[activity])
            contribution=current['balance_after_renewal_and_resource_provision_USD2041']-alternative['balance_after_renewal_and_resource_provision_USD2041']
            if contribution < 0:
                disabled.append(activity);rejected.append({'activity':activity,'incremental_balance_after_renewal_and_resource_provision_USD2041':contribution,'rule':'Suppressed in recurring portfolio; hypothetical demand and its negative contribution remain visible.'})
                current=alternative;changed=True
        if not changed:break
    current['rejected_activities']=rejected
    current['unscreened_potential_receipts_USD2041']=potential['gross_receipts_USD2041']
    return current


def calibration(q,src):
    """Scoped price comparisons; an arbitrary unit cannot repair a cash deficit."""
    ref=q['markets']['reference'];pop=q['boundary']['population'];f=q['manufacturing'];v=q['visitors'];r=q['return']
    goods=pop*f['nameplate_finished_kg_per_resident_year']*f['saleable_yield']*ref['factory_utilization']*(1-r['loss_probability_flight'])*f['delivered_quality_yield']
    nights=ref['visitors_year']*v['stay_nights'];U=goods*f['retained_U_price_per_kg']+nights*v['retained_U_per_night']
    I=src['industrial']['cases']['reference']['cost']['inflation_factor']
    return {'goods_at_Earth_kg':goods,'destination_guest_nights':nights,'fixed_basket_U':U,
            'destination_service_USD2041_per_U_by_market_price_case':{k:m['visitor_USD2026_night']*I/v['retained_U_per_night'] for k,m in q['markets'].items()},
            'hypothetical_manufacturing_delivery_basket_USD2041_per_U':{k:(goods*m['factory_USD2026_kg']+nights*m['visitor_USD2026_night'])*I/U for k,m in q['markets'].items()},
            'proposed_working_USD2041_per_U':q['normalization']['proposed_USD2041_per_U'],
            'original_pressure_component_USD2041_per_U':src['industrial']['cases']['reference']['cost']['U_calibration']['component_USD2041_per_U']['pressure'],
            'reference_component_USD2041_per_U':{'manufactured_goods_at_Earth':ref['factory_USD2026_kg']*I/f['retained_U_price_per_kg'],'destination_guest_night':ref['visitor_USD2026_night']*I/v['retained_U_per_night']},
            'original_procurement_basket_USD2041_per_U':src['industrial']['cases']['reference']['cost']['U_calibration']['USD2041_per_U'],
            'rule':'Proposed normalization from destination-service prices, with a scoped original pressure-procurement cross-check. The manufacturing basket assumes a common Earth-delivery scope that the previous U quote did not establish, and is not an independent calibration. No market exchange rate, peg or universal relative-price consistency is established; rescaling units cannot finance a deficit.'}


def run(spec=None):
    q,src=read_inputs(spec);settle=src['settlement']['cases']['reference']['resources']
    q['_operating_settlement_Earth_kg']=sum(v for k,v in settle['operating_import_kg_year'].items() if k not in ('water','air'))
    q['_renewal_settlement_Earth_kg']=sum(settle['renewal_import_kg_year'].values())
    q['_renewal_settlement_pressure_kg']=sum(v for k,v in settle['renewal_import_kg_year'].items() if 'pressure' in k)
    cases={}
    for name,(s,c,tc) in supply_cases(src).items():
        cases[name]={market:portfolio(q,s,c,tc,m) for market,m in q['markets'].items()}
    rate=calibration(q,src)
    for group in cases.values():
        for market,c in group.items():
            # Use one reference normalization across scenarios, so changing the
            # unit cannot disguise economic differences between the cases.
            x=rate['proposed_working_USD2041_per_U']
            c['proposed_normalization_U']={k.removesuffix('_USD2041')+'_U':c[k]/x for k in ['gross_receipts_USD2041','operating_cash_balance_USD2041','balance_after_renewal_and_resource_provision_USD2041']}
            feed=c['physical']['source_ore_kg_year'];water=c['physical']['source_water_kg_year']
            c['resource_horizons']={'years_for_specified_usable_feed_pools':{str(int(pool)):pool/feed for pool in q['resources']['usable_feed_pools_kg']},
               'requirements':{str(years):{'raw_feed_kg':feed*years,'water_kg':water*years,'raw_lunar_kg':c['physical']['raw_lunar_kg_year']*years,'equivalent_source_diameter_m':(6*feed*years/(math.pi*q['resources']['equivalent_body_density_kg_m3']))**(1/3)} for years in q['resources']['horizons_years']},
               'catalogue_gross_volume_duration_at_assumed_density_years':{k:g['mass_at_assumed_density_kg']/feed for k,g in src['inventory']['inventory']['groups'].items()},
               'rule':'Gross catalogue-volume ratios and specified usable-pool lifetimes are conditional arithmetic, not certified resource endurance or permission to extract entire bodies.'}
    focus=cases[q['boundary']['comparison_focus']]['reference'];I=src['industrial']['cases']['reference']['cost']['inflation_factor'];fin=q['finance']
    sensitivities={'maintenance_services_Earth_fraction':{str(f):focus['balance_after_renewal_and_resource_provision_USD2041']+(fin['maintenance_services_Earth_fraction']-f)*focus['baseline_upkeep_USD2026']*I for f in fin['maintenance_services_Earth_fraction_sensitivity']},
       'source_transition_USD2026_kg_feed':{str(f):focus['balance_after_renewal_and_resource_provision_USD2041']+(fin['source_transition_USD2026_per_kg_feed']-f)*focus['physical']['source_ore_kg_year']*I for f in fin['source_transition_sensitivity_USD2026_per_kg_feed']}}
    return {'schema':SCHEMA,'date':q['date'],'evidence':q['evidence'],'inputs':q['inputs'],
       'producer':{'path':'research/studies/initial_industrial_system/revenue.py','sha256':industry.sha(__file__),
         'configuration':{'path':'research/studies/initial_industrial_system/revenue_inputs.json','sha256':industry.sha(HERE/'revenue_inputs.json')},
         'sources':{'path':'research/studies/initial_industrial_system/revenue_sources.json','sha256':industry.sha(HERE/'revenue_sources.json')},
         'model_dependencies':[{'path':'research/studies/initial_industrial_system/'+name,'sha256':industry.sha(HERE/name)} for name in ['run.py','transport.py']]},
       'boundary':q['boundary'],'cases':cases,'U_calibration':rate,'focus_sensitivities':sensitivities,'inventory_summary':src['inventory']['inventory'],'rules':q['rules']}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'research/runs/initial_revenue.json');a=p.parse_args()
    r=run();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(r,indent=2)+'\n')
    for name,group in r['cases'].items():
        c=group['reference'];print(name,'gross/op/net B$',*[round(c[k]/1e9,2) for k in ['gross_receipts_USD2041','operating_cash_balance_USD2041','balance_after_renewal_and_resource_provision_USD2041']],'rejected',[v['activity'] for v in c['rejected_activities']])
    print('U normalization',r['U_calibration']['proposed_working_USD2041_per_U'])
if __name__=='__main__':main()
