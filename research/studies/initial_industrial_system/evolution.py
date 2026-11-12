"""Component renewal and competing use of the installed industrial capacity.

The original fleet geometry and rated capacity remain fixed at each evaluation.
New machines enter a later step; component partitions are proposed designs.
"""
from __future__ import annotations
from copy import deepcopy
import json
import math
from pathlib import Path
from . import run as base
from . import revenue
from shared.constants import JULIAN_YEAR_DAYS
HERE,ROOT=base.HERE,base.ROOT


def read_inputs(spec=None):
    q=deepcopy(spec) if spec is not None else json.loads((HERE/'economy_inputs.json').read_text())
    src={}
    for key,ref in q['inputs'].items():
        p=ROOT/ref['path']
        if base.sha(p)!=ref['sha256']:raise ValueError('Input hash mismatch: '+key)
        src[key]=json.loads(p.read_text())
    s=base.with_case(src['industrial_inputs'],'local_fleet_fabrication')
    c=src['transport']['source_structure_localization']['space_chemical']
    return q,src,s,c


def component_class(name):
    if 'pressure' in name:return 'pressure'
    if 'computer' in name or 'communications' in name:return 'computer'
    if 'solar' in name or 'beam_constellation' in name or name.endswith('_power') or 'receivers' in name:return 'solar'
    if 'radiator' in name:return 'radiator'
    if 'lander' in name:return 'lander'
    if 'tug' in name or 'hull' in name:return 'vehicle'
    if 'chemistry' in name or 'cryocooling' in name:return 'chemistry'
    if 'tank' in name or 'depot' in name or name.endswith('_storage'):return 'tank'
    if 'structure' in name:return 'building'
    if any(x in name for x in ['trains','cells','robot','excavator','hauler']):return 'machine'
    return 'ordinary'


def stocks(q,src,c,goods=0,compute=0,service=None,refinery_capacity=0):
    """Existing asset inventories plus capacity-scaled commercial additions."""
    items=[]
    for name,a in c['assets'].items():
        frac=a.get('import_fraction',1)
        price=a['manufacture_USD2026']/a['mass_kg']/frac if a['mass_kg'] and frac else 0
        items.append(dict(name=name,mass=a['mass_kg'],price=price,life=a['life_years'],destination=a['destination'],owner='industry',inherited_local=1-frac))
    city=src['settlement']['cases']['reference']
    for name,annual in city['resources']['renewal_import_kg_year'].items():
        items.append(dict(name=name,mass=annual*city['maintenance']['life_years'][name],price=src['industrial_inputs']['prices']['pressure_USD2026_kg' if 'pressure' in name else 'other_special_USD2026_kg'],life=city['maintenance']['life_years'][name],destination='surface',owner='settlement',inherited_local=0))
    e=q['exports'];plate=goods/(e['factory_yield']*e['factory_utilization'])
    if goods:
        items.append(dict(name='export_factory',mass=plate*e['factory_equipment_kg_per_nameplate_kg_year'],price=src['industrial_inputs']['prices']['lunar_hardware_USD2026_kg'],life=e['factory_equipment_life_years'],destination='surface',owner='commerce',inherited_local=0))
    if compute:
        items.append(dict(name='source_computer',mass=compute*e['compute_IT_kg_W'],price=e['compute_IT_USD2026_W']/e['compute_IT_kg_W'],life=e['compute_life_years'],destination='source',owner='commerce',inherited_local=0))
        items.append(dict(name='source_communications',mass=compute*e['compute_communications_kg_W'],price=e['compute_communications_USD2026_W']/e['compute_communications_kg_W'],life=e['compute_life_years'],destination='source',owner='commerce',inherited_local=0))
    if refinery_capacity:
        mass=sum(c['assets'][n]['mass_kg'] for n in ['source_metal_trains','source_fabrication_cells'])*refinery_capacity/c['peak_flows']['source_metal_kg_year']
        items.append(dict(name='source_extra_metal_trains',mass=mass,price=src['industrial_inputs']['prices']['source_mine_USD2026_kg'],life=c['assets']['source_metal_trains']['life_years'],destination='source',owner='commerce',inherited_local=0))
    if service:
        for name,a in service['assets'].items():
            if a['mass_kg']:
                items.append(dict(name=name,mass=a['mass_kg'],price=a['manufacture_USD2026']/a['mass_kg'],life=a['life_years'],destination=a['destination'],owner='commerce',inherited_local=0))
    return items


def services(src,s,c,enabled=False):
    if not enabled:return None
    spec=deepcopy(src['revenue_inputs']);city=src['settlement']['cases']['reference']
    spec['_operating_settlement_Earth_kg']=sum(v for k,v in city['resources']['operating_import_kg_year'].items() if k not in ('air','water'))
    market=deepcopy(spec['markets']['reference']);market['propellant_customer_kg_year']=0
    return revenue.commerce(spec,s,c,market,disabled=('manufacturing','computing'))


def components(q,items,profile,pop_scale,plant_scale,growth_scale=0,new_plant_scale=0,commerce_growth=0,renewal_factor=1):
    cfg=q['profiles'][profile];l=q['production'];result=[]
    for a in items:
        kind=component_class(a['name']);part=q['components'][kind]
        scale=plant_scale if a['owner']=='industry' else pop_scale if a['owner']=='settlement' else 1
        grow=new_plant_scale if a['owner']=='industry' else growth_scale if a['owner']=='settlement' else commerce_growth
        if cfg['modular']:
            struct=max(a['inherited_local'],part['structural_mass_fraction'])
            local=struct if cfg['localize'] else a['inherited_local']
            life=max(a['life'],part['structural_life_years'])
        else:struct=a['inherited_local'];local=struct;life=a['life']
        mass=a['mass'];renew_struct=mass*struct*scale/life*renewal_factor;renew_active=mass*(1-struct)*scale/a['life']*renewal_factor
        local_renew=mass*local*scale/life*renewal_factor;local_new=mass*local*grow
        imported_renew=renew_struct+renew_active-local_renew
        imported_new=mass*(1-local)*grow
        multiplier=l['special_component_price_multiplier'] if cfg['modular'] else 1
        # Imported active components receive the explicit price-density sensitivity;
        # imported passive structures retain the inherited unit price.
        renew_cost=((renew_struct-local_renew)+renew_active*multiplier)*a['price']
        growth_cost=(mass*(struct-local)*grow+mass*(1-struct)*grow*multiplier)*a['price']
        local_price=l['local_fabrication_value_USD2026_kg']
        result.append({**a,'class':kind,'stock_scale':scale,'structural_share':struct,'local_share':local,'structural_life':life,
            'local_renew_kg':local_renew,'local_growth_kg':local_new,'import_renew_kg':imported_renew,'import_growth_kg':imported_new,
            'import_renew_USD2026':renew_cost,'import_growth_USD2026':growth_cost,
            'local_renew_value_USD2026':local_renew*local_price,'local_growth_value_USD2026':local_new*local_price,
            'annual_CFC_USD2026':renew_cost+local_renew*local_price})
    return result


def annual(q,src,s,c,profile='inherited',population=250000,plant_scale=1,growth_people=0,new_plant_scale=0,goods=0,compute=0,propellant_sales=0,source_goods=0,extra_metal=0,service_portfolio=False,commission_fraction=0,qualification=False,renewal_factor=1,refinery_capacity=0):
    """One year's material/power demand against the capacity already installed.

Growth uses the founding system's quantities per supported resident as an
intensity comparison. It does not approve repeated four-city/100-town siting.
"""
    pop0=c['demand']['population'];ps=population/pop0;gs=growth_people/pop0
    d=c['demand'];t=s['transport'];z=s['source'];l=s['lunar'];e=s['economics'];peak=c['peak_flows'];old=c['steady_flows'];ex=q['exports'];prod=q['production'];K=plant_scale
    cfg=q['profiles'][profile]
    b=services(src,s,c,service_portfolio)
    parts=components(q,stocks(q,src,c,goods,compute,b,refinery_capacity),profile,ps,K,gs,new_plant_scale,commission_fraction,renewal_factor)
    R=renewal_factor
    imports={v:0. for v in ['surface','depot','source','leo']}
    renew_imports=imports.copy();growth_imports=imports.copy();local={v:0. for v in imports};retired=local.copy();local_growth=local.copy()
    for a in parts:
        dest=a['destination'];renew_imports[dest]+=a['import_renew_kg'];growth_imports[dest]+=a['import_growth_kg']
        local[dest]+=a['local_renew_kg']+a['local_growth_kg'];retired[dest]+=a['local_renew_kg'];local_growth[dest]+=a['local_growth_kg']
    city=src['settlement']['cases']['reference'];op=city['resources']['operating_import_kg_year']
    resident_consumables=sum(v for k,v in op.items() if k not in ('air','water'))*ps
    # Ordinary process supplies supplement module renewal; the inherited full
    # percentage-of-capital upkeep remains a financial comparison, not a mass.
    normal_parts=parts if renewal_factor else components(q,stocks(q,src,c,goods,compute,b,refinery_capacity),profile,ps,K)
    industrial_consumables=prod['process_consumables_fraction_of_imported_renewal']*sum(a['import_renew_kg'] for a in normal_parts if a['owner']=='industry') if cfg['modular'] else 0
    for dest in imports:imports[dest]=renew_imports[dest]+growth_imports[dest]
    imports['surface']+=resident_consumables+(b['operating_surface_kg'] if b else 0)
    imports['source']+=industrial_consumables
    imports['surface']+=d['earth_special']*gs-sum(a['mass'] for a in parts if a['owner']=='settlement')*gs
    if imports['surface']<0:raise ValueError('Growth equipment exceeds inherited settlement import stock')
    recovered=prod['source_recovery']
    surface_virgin=local_growth['surface']+retired['surface']*(1-l['recovery'])
    source_product=source_goods/ex['factory_yield']+source_goods*.05+extra_metal*(1+q['coproducts']['stock_containment_kg_per_kg'])
    source_virgin=sum(local_growth[k]+retired[k]*(1-recovered) for k in ['source','depot','leo'])+source_product
    surface_alloy=surface_virgin*prod['alloy_import_fraction_of_virgin_local_metal'] if cfg['localize'] else 0
    source_alloy=source_virgin*prod['alloy_import_fraction_of_virgin_local_metal'] if cfg['localize'] else 0
    imports['surface']+=surface_alloy;imports['source']+=source_alloy
    # All source-produced structures used away from the source travel outward;
    # recovered off-site metal returns through the same explicit network.
    source_to_depot_structures=local['depot']+local['leo']
    depot_to_source_recycled=(retired['depot']+retired['leo'])*recovered
    leo_to_depot_recycled=retired['leo']*recovered
    depot_to_leo_structures=local['leo']
    gross_feed=goods/ex['factory_yield'];earth_factory_feed=gross_feed*ex['factory_Earth_input_fraction']
    factory_local=gross_feed-earth_factory_feed+goods*.05
    imports['surface']+=earth_factory_feed
    # Cargo-entry hardware retains the prior explicit construction and turnover
    # assumptions; no passenger use of the cargo capsule is inferred.
    rv=src['revenue_inputs'];rr=rv['return'];packed=goods*(1+rv['manufacturing']['transport_packaging_mass_fraction']);source_packed=source_goods*(1+rv['manufacturing']['transport_packaging_mass_fraction']);packed+=source_packed;flights=packed/rr['capsule_payload_kg']
    bus=rr['capsule_payload_kg']*rr['capsule_dry_to_payload']*(1-rr['shield_fraction_of_dry']);shield=rr['capsule_payload_kg']*rr['capsule_dry_to_payload']*rr['shield_fraction_of_dry']
    expected=sum((1-rr['loss_probability_flight'])**n for n in range(rr['durable_life_flights']))
    deorbit=(rr['capsule_payload_kg']+bus+shield)*math.expm1(rr['deorbit_dv_m_s']/(rr['deorbit_Isp_s']*base.STANDARD_GRAVITY))
    imports['leo']+=flights*(bus+shield+deorbit)
    water=d['yearly_water']*ps+d['water']*gs;air=d['yearly_air']*ps+d['air']*gs
    if b:water+=b['guest_water_kg_year']+b['initial_water_kg']*commission_fraction;air+=b['guest_air_kg_year']+b['initial_air_kg']*commission_fraction
    lunar_bulk=d['renewal_bulk']*ps*R+d['bulk']*gs+d['goods_virgin']*ps+peak['additional_bulk_stock_kg']*(K*R/l['bulk_life_years']+new_plant_scale)
    lunar_metal=d['renewal_metal']*ps*R+d['metal']*gs+peak['additional_metal_stock_kg']*(K*R/l['metal_life_years']+new_plant_scale)+local['surface']+factory_local
    factory_plate=goods/(ex['factory_yield']*ex['factory_utilization']);factory_floor=factory_plate*ex['factory_floor_m2_per_nameplate_kg_year']
    lunar_bulk+=factory_floor*ex['factory_bulk_kg_m2']*R/l['bulk_life_years'];lunar_metal+=factory_floor*ex['factory_metal_kg_m2']*R/l['metal_life_years']
    virgin_bulk=(d['renewal_bulk']*ps*R+peak['additional_bulk_stock_kg']*K*R/l['bulk_life_years']+factory_floor*ex['factory_bulk_kg_m2']*R/l['bulk_life_years'])*(1-l['recovery'])+d['goods_virgin']*ps+d['bulk']*gs+peak['additional_bulk_stock_kg']*new_plant_scale
    virgin_metal=(d['renewal_metal']*ps*R+peak['additional_metal_stock_kg']*K*R/l['metal_life_years']+factory_floor*ex['factory_metal_kg_m2']*R/l['metal_life_years'])*(1-l['recovery'])+d['metal']*gs+peak['additional_metal_stock_kg']*new_plant_scale+surface_virgin-surface_alloy+factory_local
    commercial_bulk=factory_floor*ex['factory_bulk_kg_m2']+(b['initial_bulk_kg'] if b else 0)
    commercial_metal=factory_floor*ex['factory_metal_kg_m2']+(b['initial_metal_kg'] if b else 0)
    # The factory floor replacement is already present above; add service civil stock.
    service_bulk=(b['initial_bulk_kg']/l['bulk_life_years']*R if b else 0)
    service_metal=(b['initial_metal_kg']/l['metal_life_years']*R if b else 0)
    lunar_bulk+=commercial_bulk*commission_fraction+service_bulk
    lunar_metal+=commercial_metal*commission_fraction+service_metal
    virgin_bulk+=commercial_bulk*commission_fraction+service_bulk*(1-l['recovery'])
    virgin_metal+=commercial_metal*commission_fraction+service_metal*(1-l['recovery'])
    # Factory air is additional; the sealed industrial floor uses the earlier leak allowance.
    factory_air=factory_floor*rv['manufacturing']['air_kg_m2']
    air+=factory_air*(commission_fraction+d['yearly_air']/d['air'])
    qualification_hardware=prod['qualification_import_hardware_kg']/prod['qualification_years'] if qualification else 0
    imports['source']+=qualification_hardware
    raw_lunar=virgin_bulk/l['bulk_yield']+virgin_metal/l['metal_yield'];handling=raw_lunar+d['loose_and_grading']*gs
    local_W=(lunar_bulk*l['bulk_kWh_kg']+lunar_metal*l['metal_kWh_kg']+handling*l['handling_kWh_kg']+local['surface']*prod['local_fabrication_finish_kWh_kg']+goods*ex['factory_finish_kWh_kg'])*1000/base.HOURS+(b['added_surface_mean_W'] if b else 0)
    ecargo=imports['surface']+imports['depot']+imports['source']+leo_to_depot_recycled
    v=peak['Earth_vehicle'];ret=src['transport_inputs']['chemical']['handling_retention']
    earth=revenue.ferry(ecargo,packed+depot_to_leo_structures,payload=t['earth_payload_kg'],dry=v['dry_kg'],dv=v['delta_v_oneway_m_s'],isp=v['Isp_s'],margin=t['propellant_margin'],mode='space_chemical',retention=ret,tank_limit=v['burns']['peak_tank_contents_kg'])
    mooncargo=imports['surface']+water+air;P=t['lander_payload_kg'];D=P*t['lander_dry_to_payload']
    landlimit=P*base.round_trip(t['lander_dv_oneway_m_s'],t['Isp_s'],D/P)['total']*t['propellant_margin']
    moon=revenue.ferry(mooncargo,packed-source_packed,payload=P,dry=D,dv=t['lander_dv_oneway_m_s'],isp=t['Isp_s'],margin=t['propellant_margin'],mode='earth_chemical',retention=1,tank_limit=landlimit)
    N=air*z['air_nitrogen_mass_fraction'];O=air-N
    depot_water=moon['burned_kg_year']+O+earth['space_supply_kg_year']+propellant_sales
    sourcecargo=water+N+depot_water+source_to_depot_structures+source_packed
    sourceback=imports['source']+depot_to_source_recycled
    sdry=t['nea_payload_kg']*t['nea_dry_to_payload'];source_limit=t['nea_payload_kg']*peak['nea_round_trip']['total']*t['propellant_margin']
    source=revenue.ferry(sourcecargo,sourceback,payload=t['nea_payload_kg'],dry=sdry,dv=t['nea_dv_oneway_m_s'],isp=t['Isp_s'],margin=t['propellant_margin'],mode='earth_chemical',retention=1,tank_limit=source_limit)
    for route in [moon,source]:route['Earth_fuel_kg_year']=0;route['space_supply_kg_year']=route['burned_kg_year']
    source_water=water+depot_water+source['burned_kg_year']
    source_metal=sum(local[k] for k in ['source','depot','leo'])+source_product
    ore=max(source_water/(z['water_grade']*z['water_recovery']),N/(z['nitrogen_grade']*z['nitrogen_recovery']),(source_virgin-source_alloy)/(z['metal_grade']*z['metal_recovery']))
    process=ore*z['ore_kWh_kg']+N*z['nitrogen_kWh_kg']+source_metal*z['metal_refining_and_fabrication_kWh_kg']
    source_chem=source['burned_kg_year']/z['water_to_hydrogen_mass_ratio']*z['electrolysis_liquefaction_kWh_kg_H2']
    source_W=(process+source_chem)*1000/base.HOURS+compute
    hf=1-z['hydrogen_HHV_kWh_kg']/z['electrolysis_liquefaction_kWh_kg_H2']
    source_heat=(process+source_chem*hf)*1000/base.HOURS+compute
    depot_W=depot_water/z['water_to_hydrogen_mass_ratio']*z['electrolysis_liquefaction_kWh_kg_H2']*1000/base.HOURS
    def trips(n,days):return n*K*JULIAN_YEAR_DAYS/days*t['fleet_availability']/e['capacity_margin']
    ratios={'lunar_bulk':lunar_bulk/(peak['bulk_kg_year']*K),'lunar_metal':lunar_metal/(peak['metal_kg_year']*K),'lunar_handling':handling/(peak['handling_kg_year']*K),'lunar_process_power':local_W/(peak['local_mean_W']*K),
      'ground_firm_power':(d['occupied_W']*ps+local_W*e['capacity_margin']/l['process_availability'])/(c['engineering']['total_ground_firm_W']*K),
      'source_metal':source_metal/(peak['source_metal_kg_year']*K+refinery_capacity),'source_feed':ore/(peak['source_ore_kg_year']*K),'source_power':source_W/(peak['source_mean_W']*K),'source_heat':source_heat/(peak['source_heat_W']*K),'depot_power':depot_W/(peak['depot_mean_W']*K),
      'Earth_ferry':earth['trips_year']/trips(c['engineering']['counts']['Earth_route_ships'],v['cycle_days']),
      'lunar_ferry':moon['trips_year']/trips(c['engineering']['counts']['lunar_landers'],JULIAN_YEAR_DAYS/t['lander_cycles_year']),
      'source_ferry':source['trips_year']/trips(c['engineering']['counts']['asteroid_tugs'],t['nea_cycle_years']*JULIAN_YEAR_DAYS)}
    build=src['settlement_inputs']['production'];demand=city['schedule']['production_demand']
    ratios['floor_installation']=(demand['floor']*gs+factory_floor*commission_fraction+(b['initial_bulk_kg']/rv['services']['bulk_kg_m2']*commission_fraction if b else 0))/(build['floor_m2_year']*K)
    ratios['route_installation']=demand['route']*gs/(build['route_km_year']*K)
    saleable=(goods+source_goods)*(1-rr['loss_probability_flight'])*rv['manufacturing']['delivered_quality_yield']
    import_make=sum(a['import_renew_USD2026']+a['import_growth_USD2026'] for a in parts)
    residual_growth_import=d['earth_special']*gs-sum(a['mass'] for a in parts if a['owner']=='settlement')*gs
    import_make+=qualification_hardware*s['prices']['source_mine_USD2026_kg']
    import_make+=max(0,residual_growth_import)*s['prices']['other_special_USD2026_kg']
    operating_import=(resident_consumables+industrial_consumables+surface_alloy+source_alloy+earth_factory_feed)*s['prices']['other_special_USD2026_kg']
    capsule_renew=flights*bus/expected*s['prices']['spacecraft_USD2026_kg']*R
    capsule_op=flights*(shield*s['prices']['spacecraft_USD2026_kg']+deorbit*t['earth_propellant_USD2026_kg'])+packed*rr['Earth_ground_service_USD2026_per_payload_kg']
    leo_mass=ecargo-leo_to_depot_recycled+imports['leo']
    launch_cost=leo_mass*t['LEO_USD2026_kg']
    remote_staff=e['staff']*K*.5*e['loaded_staff_USD2026_year']
    if cfg['modular']:external_services=remote_staff+prod.get('external_engineering_contract_USD2026_year',2e9)*K
    else:
        old_staff=e['staff']*e['loaded_staff_USD2026_year'];external_services=(c['cost']['annual_supplier_operation_USD2041']/c['cost']['inflation_factor']-old_staff)*src['revenue_inputs']['finance']['maintenance_services_Earth_fraction']*K+remote_staff
    service_operating=sum(b['operating_purchases_USD2026'].values()) if b else 0
    # Allocate the qualification programme over its explicit five-year trial;
    # hardware is included in that budget, with freight charged above.
    qualification_services=max(0,prod['qualification_capital_USD2026']-prod['qualification_import_hardware_kg']*s['prices']['source_mine_USD2026_kg'])/prod['qualification_years'] if qualification else 0
    capsule_stock=math.ceil(flights*e['capacity_margin']/(rr['cycles_year']*t['fleet_availability']))*(bus+shield)*s['prices']['spacecraft_USD2026_kg']
    commissioning_services=(factory_plate*rv['manufacturing']['civil_and_commissioning_USD2026_per_nameplate_kg_year']+(b['civil_and_commissioning_USD2026'] if b else 0))*commission_fraction
    commissioning_services+=sum(a['import_growth_USD2026']+a['local_growth_value_USD2026'] for a in parts if a['owner']=='commerce' and a['destination']=='source')*q['valuation']['stock_installation_fraction_of_hardware_value']
    import_make+=capsule_stock*commission_fraction
    external_services+=service_operating+qualification_services
    stage_imports=import_make+operating_import+capsule_renew+capsule_op+launch_cost+external_services
    resource=ore*src['revenue_inputs']['finance']['source_transition_USD2026_per_kg_feed']+raw_lunar*e['local_restoration_USD2026_kg_feed']
    return {'profile':profile,'population':population,'plant_scale':K,'growth_people_year':growth_people,'new_plant_scale_year':new_plant_scale,'finished_export_kg_year':goods,'source_export_kg_year':source_goods,'stock_metal_kg_year':extra_metal,'refinery_capacity_kg_year':refinery_capacity,'commission_fraction':commission_fraction,'qualification':qualification,'service_portfolio':service_portfolio,'service':b,'saleable_export_kg_year':saleable,'compute_W':compute,'propellant_sales_kg_year':propellant_sales,'components':parts,
      'capacity_use':ratios,'binding_capacity':max(ratios,key=ratios.get),'feasible':max(ratios.values())<=1+1e-8,
      'physical':{'lunar_bulk_kg_year':lunar_bulk,'lunar_metal_kg_year':lunar_metal,'lunar_raw_kg_year':raw_lunar,'lunar_handling_kg_year':handling,'lunar_process_W':local_W,'resident_mean_W':d['occupied_mean_W']*ps,'source_metal_kg_year':source_metal,'source_virgin_metal_kg_year':source_virgin-source_alloy,'source_feed_kg_year':ore,'source_water_kg_year':source_water,'source_nitrogen_kg_year':N,'source_retained_tails_kg_year':ore-source_water-N-source_virgin+source_alloy,'source_process_W':source_W,'source_heat_W':source_heat,'depot_W':depot_W,'Earth':earth,'Moon':moon,'source':source,'LEO_kg_year':leo_mass,'Earth_launches_year':math.ceil(leo_mass/t['launch_payload_kg']),'local_fabrication_kg_year':local,'local_growth_kg_year':local_growth,'retired_components_kg_year':retired,'recycled_to_source_kg_year':depot_to_source_recycled,'recycled_from_LEO_kg_year':leo_to_depot_recycled,'imports_by_destination_kg_year':imports,'technical_resident_jobs':e['staff']*K*.5+(b['additional_resident_jobs'] if b else 0)+factory_plate/(c['demand']['population']*rv['manufacturing']['nameplate_finished_kg_per_resident_year'])*rv['manufacturing']['staff']+compute/1e6*rv['computing']['staff_per_MW']+city['resources']['food_operators']*ps+city['resources']['industry_operators']*ps,'renewal_imports_by_destination_kg_year':renew_imports,'growth_imports_by_destination_kg_year':growth_imports,'source_product_kg_year':source_product,'commercial_bulk_stock_kg':commercial_bulk,'commercial_metal_stock_kg':commercial_metal,'resident_water_kg_year':water,'resident_air_kg_year':air,'factory_Earth_feed_kg_year':earth_factory_feed},
      'cash_USD2026':{'equipment_imports':import_make,'commercial_civil_installation_allowance':commissioning_services,'operating_material_imports':operating_import,'capsule_replacement':capsule_renew,'capsule_operations':capsule_op,'Earth_launch_service':launch_cost,'industrial_external_services':external_services,'total_external_purchases':stage_imports,'resource_provision':resource},
      'capital_formation':{'local_renewal_value_USD2026':sum(a['local_renew_value_USD2026'] for a in parts),'local_growth_value_USD2026':sum(a['local_growth_value_USD2026'] for a in parts),'imported_renewal_USD2026':sum(a['import_renew_USD2026'] for a in parts),'imported_growth_USD2026':sum(a['import_growth_USD2026'] for a in parts)+max(0,residual_growth_import)*s['prices']['other_special_USD2026_kg'],'component_CFC_USD2026':sum(a['annual_CFC_USD2026'] for a in parts)}}


def frontier(q,src,s,c,variable,profile='local_structures',**fixed):
    balanced=fixed.pop('balanced_growth',False)
    start=annual(q,src,s,c,profile=profile,**fixed)
    if not start['feasible']:return {'status':'base_infeasible','maximum':0,'at_limit':start}
    lo=0.;hi={'goods':1e8,'source_goods':1e8,'extra_metal':1e8,'compute':1e9,'growth_people':10000,'new_plant_scale':.02,'propellant_sales':1e8}[variable]
    def trial(v):
        args={**fixed,variable:v}
        if variable=='growth_people' and balanced:args['new_plant_scale']=v/c['demand']['population']
        args.pop('balanced_growth',None)
        return annual(q,src,s,c,profile=profile,**args)
    for _ in range(40):
        if not trial(hi)['feasible']:break
        hi*=2
    for _ in range(55):
        mid=(lo+hi)/2
        if trial(mid)['feasible']:lo=mid
        else:hi=mid
    return {'status':'capacity_only','maximum':lo,'at_limit':trial(lo),'rule':'One-at-a-time physical frontier with the other specified uses fixed; maxima cannot be added, and do not establish demand or financing.'}


def coproducts(q,s,c,a):
    """Joint-feed screen; raw retained value is not earned annual cash."""
    f=a['physical'];z=s['source'];feed=f['source_feed_kg_year'];actual={'water':f['source_water_kg_year'],'nitrogen':f['source_nitrogen_kg_year'],'aggregate_recoverable_metals':f['source_virgin_metal_kg_year']}
    grade={'water':z['water_grade'],'nitrogen':z['nitrogen_grade'],'aggregate_recoverable_metals':z['metal_grade']}
    recovery={'water':z['water_recovery'],'nitrogen':z['nitrogen_recovery'],'aggregate_recoverable_metals':z['metal_recovery']}
    rows={k:{'feed_element_or_water_equivalent_kg_year':feed*grade[k],'recovery_scenario_potential_kg_year':feed*grade[k]*recovery[k],'actually_required_virgin_kg_year':actual[k],'unused_recovery_potential_kg_year':max(0,feed*grade[k]*recovery[k]-actual[k])} for k in actual}
    metal_spare=max(0,c['peak_flows']['source_metal_kg_year']*a['plant_scale']-f['source_metal_kg_year'])
    available=min(rows['aggregate_recoverable_metals']['unused_recovery_potential_kg_year'],metal_spare)
    available=min(available,max(0,c['peak_flows']['source_mean_W']*a['plant_scale']-f['source_process_W'])*base.HOURS/1000/z['metal_refining_and_fabrication_kWh_kg'])
    available=min(available,max(0,c['peak_flows']['source_heat_W']*a['plant_scale']-f['source_heat_W'])*base.HOURS/1000/z['metal_refining_and_fabrication_kWh_kg'])
    target=q['coproducts']['metal_refinery_expansion_kg_year'];ratio=target/c['peak_flows']['source_metal_kg_year']
    hardware=sum(c['assets'][k]['mass_kg'] for k in ['source_metal_trains','source_fabrication_cells'])*ratio
    part=q['components']['machine']['structural_mass_fraction'];local=hardware*part;imported=hardware-local
    # This is a design requirement for additional trains, not an installed plant
    # or credited revenue; fabrication and power must compete with other uses.
    return {'joint_feed_kg_year':feed,'products':rows,'current_extra_metal_refining_capacity_kg_year':available,
            'controlled_residue_kg_year':f['source_retained_tails_kg_year'],'unassayed_materials':q['coproducts']['unassayed'],
            'refinery_expansion_requirement':{'extra_capacity_kg_year':target,'extra_process_W':target*z['metal_refining_and_fabrication_kWh_kg']*1000/base.HOURS,'new_hardware_kg':hardware,'local_structures_kg':local,'Earth_functional_hardware_kg':imported,'hardware_manufacture_USD2026':local*q['production']['local_fabrication_value_USD2026_kg']+imported*s['prices']['source_mine_USD2026_kg'],'stage':'Unbuilt capacity comparison; delivery, installation, qualification and financing are additional.'},
            'raw_residue_revenue_USD2026':0,'rule':'Recoverable-equivalent content is not separated product, proved reserve, sale or added annual revenue. Joint extraction costs are charged once. Actual internally used metal enters fabrication; additional refinement requires capacity and buyers or funded stock accumulation.'}
