"""Stocks, dated procurement and physical requirements for an enclosed economy.

The model consumes a frozen engineering product. It never imports the former
research implementation or turns population into a fleet multiplier.
"""
from __future__ import annotations
from copy import deepcopy
from functools import lru_cache
import math
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, STANDARD_GRAVITY, STEFAN_BOLTZMANN

HOURS=JULIAN_DAY*JULIAN_YEAR_DAYS/3600
DESTS=('surface','source','depot','leo')


def burns(cargo,back,dry,dv,isp,margin,mode,retention=1.):
    k=math.expm1(dv/(isp*STANDARD_GRAVITY))
    if mode=='earth_chemical':
        b=(dry+back)*k;a=(cargo+dry+b)*k
        a*=margin;b*=margin
        return dict(out=a,back=b,tank=a+b,burned=a+b,space_supply=a+b)
    a=(cargo+dry)*k*margin;carried=a/retention
    b=(dry+carried+back)*k*margin
    return dict(out=a,back=b,tank=carried+b,burned=a+b,space_supply=(carried+b)/retention)


def ferry(cargo,back,route):
    """Algebraic equivalent of the inherited affine tank-limited trip solver."""
    args={k:route[k] for k in ('dry','dv','isp','margin','mode','retention')}
    empty=burns(0,0,**args);unit_args=dict(args,dry=0.);one=burns(1,0,**unit_args);ret=burns(0,1,**unit_args)
    room=route['tank_limit']-empty['tank']
    if room<=0:raise ValueError('Vehicle cannot carry its own return propellant')
    n=max(cargo/route['payload'],back/route['payload'],
          (one['tank']*cargo+ret['tank']*back)/room)
    if n<=0:return dict(trips=0.,out=0.,back=0.,tank=0.,burned=0.,space_supply=0.)
    b=burns(cargo/n,back/n,**args)
    return dict(trips=n,**{k:v*n for k,v in b.items()})


def routes(data):
    s=data['industrial'];t=s['transport'];v=data['original_architecture']['peak_flows']['Earth_vehicle']
    out={'Earth':dict(payload=t['earth_payload_kg'],dry=v['dry_kg'],dv=v['delta_v_oneway_m_s'],isp=v['Isp_s'],
        margin=t['propellant_margin'],mode='space_chemical',retention=data['transport']['chemical']['handling_retention'],
        tank_limit=v['burns']['peak_tank_contents_kg'],cycles=JULIAN_YEAR_DAYS/v['cycle_days'])}
    for name,prefix in [('Moon','lander'),('source','nea')]:
        payload=t[prefix+'_payload_kg'];dry=payload*t[prefix+'_dry_to_payload'];dv=t[prefix+'_dv_oneway_m_s']
        args=dict(dry=dry,dv=dv,isp=t['Isp_s'],margin=t['propellant_margin'],mode='earth_chemical',retention=1.)
        tank=burns(payload,0,**args)['tank']
        if name=='source':tank=payload*data['original_architecture']['peak_flows']['nea_round_trip']['total']*t['propellant_margin']
        out[name]=dict(payload=payload,**args,tank_limit=tank,cycles=t['lander_cycles_year'] if name=='Moon' else 1/t['nea_cycle_years'])
    return out


def catalogue(data):
    """One mass unit per asset; active and structural lifetimes stay distinct."""
    out=deepcopy(data['templates']);q=data['economy'];s=data['industrial'];ex=q['exports'];rv=data['revenue']
    for a in out.values():a.update(material='metal',civil=False)
    def add(name,owner,price,life,share=0.,structlife=None,dest='surface',material='metal',civil=False):
        out[name]=dict(name=name,owner=owner,price=price,life=life,structural_share=share,
            structural_life=structlife or life,destination=dest,material=material,civil=civil,mass=0.)
    add('export_factory','commerce',s['prices']['lunar_hardware_USD2026_kg'],ex['factory_equipment_life_years'],.6,30)
    add('source_computer','commerce',ex['compute_IT_USD2026_W']/ex['compute_IT_kg_W'],ex['compute_life_years'])
    out['source_computer']['destination']='source'
    add('source_communications','commerce',ex['compute_communications_USD2026_W']/ex['compute_communications_kg_W'],ex['compute_life_years'],dest='source')
    sv=rv['services']
    add('service_equipment','service',s['prices']['lunar_hardware_USD2026_kg'],sv['research_hardware_life_years'],.6,30)
    add('visitor_equipment','hospitality',s['prices']['lunar_hardware_USD2026_kg'],rv['visitors']['equipment_life_years'],.6,30)
    rr=rv['return']
    add('return_capsule_buses','commerce',s['prices']['spacecraft_USD2026_kg'],rr['durable_calendar_life_years'],dest='leo')
    for name,owner,mat,life in [('residential_bulk','settlement','bulk',s['lunar']['bulk_life_years']),
        ('residential_roads','settlement','bulk',s['lunar']['road_life_years']),('residential_metal','settlement','metal',s['lunar']['metal_life_years']),
        ('industrial_civil_bulk','industry','bulk',s['lunar']['bulk_life_years']),('industrial_civil_metal','industry','metal',s['lunar']['metal_life_years']),
        ('commercial_civil_bulk','commerce','bulk',s['lunar']['bulk_life_years']),('commercial_civil_metal','commerce','metal',s['lunar']['metal_life_years']),
        ('service_civil_bulk','service','bulk',s['lunar']['bulk_life_years']),('service_civil_metal','service','metal',s['lunar']['metal_life_years']),
        ('visitor_civil_bulk','hospitality','bulk',s['lunar']['bulk_life_years']),('visitor_civil_metal','hospitality','metal',s['lunar']['metal_life_years'])]:
        add(name,owner,q['production']['local_bulk_value_USD2026_kg'] if mat=='bulk' else q['production']['structural_metal_value_USD2026_kg'],life,1.,life,material=mat,civil=True)
    return out


def commercial_targets(data,row):
    ex=data['economy']['exports'];rr=data['revenue']['return'];sv=data['revenue']['services'];s=data['industrial']
    goods=row.get('goods_capacity',row.get('goods',0.));compute=row.get('compute_capacity',row.get('compute',0.))
    plate=goods/(ex['factory_yield']*ex['factory_utilization']);jobs=row.get('service_jobs',0.)
    floor=plate*ex['factory_floor_m2_per_nameplate_kg_year'];service_floor=jobs*sv['research_floor_m2_staff'];beds=row.get('visitor_beds',0.);v=data['revenue']['visitors']
    flights=(row.get('goods',0.)+row.get('source_goods',0.))*(1+data['revenue']['manufacturing']['transport_packaging_mass_fraction'])/rr['capsule_payload_kg']
    capsules=math.ceil(flights*s['economics']['capacity_margin']/(rr['cycles_year']*s['transport']['fleet_availability']))
    bus=rr['capsule_payload_kg']*rr['capsule_dry_to_payload']*(1-rr['shield_fraction_of_dry'])
    return dict(export_factory=plate*ex['factory_equipment_kg_per_nameplate_kg_year'],
        source_computer=compute*ex['compute_IT_kg_W'],source_communications=compute*ex['compute_communications_kg_W'],
        service_equipment=jobs*sv['research_hardware_kg_staff'],return_capsule_buses=capsules*bus,
        commercial_civil_bulk=floor*ex['factory_bulk_kg_m2'],commercial_civil_metal=floor*ex['factory_metal_kg_m2'],
        service_civil_bulk=service_floor*sv['bulk_kg_m2'],service_civil_metal=service_floor*sv['metal_kg_m2'],
        visitor_equipment=beds*v['surface_hardware_kg_bed'],visitor_civil_bulk=beds*v['local_bulk_kg_bed'],visitor_civil_metal=beds*v['local_metal_kg_bed'])


def residential_targets(data,population):
    f=population/data['normalization']['population'];d=data['demand']
    out={n:a['mass']*f for n,a in data['templates'].items() if a['owner']=='settlement'}
    out.update(residential_bulk=(d['bulk']-d['roads'])*f,residential_roads=d['roads']*f,residential_metal=d['metal']*f)
    return out


def blank_flow():
    return dict(imports={k:0. for k in DESTS},local={k:0. for k in DESTS},recycled={k:0. for k in DESTS},
        bulk=0.,bulk_recycled=0.,external_capital=0.,external_replacement=0.,local_capital=0.,local_replacement=0.,
        service_capital=0.,service_replacement=0.,retired_mass=0.,new_floor=0.,capsule_purchases_kg=0.,new_visitor_beds=0.,industry_replacement_import_kg=0.,
        service_external_capital=0.,service_external_replacement=0.,hospitality_external_capital=0.,hospitality_external_replacement=0.,capital_imports={k:0. for k in DESTS},capital_local={k:0. for k in DESTS},capital_bulk=0.,hospitality_capital=0.,hospitality_replacement=0.)


@lru_cache(maxsize=None)
def lifetime_scale(life,shape):
    # Annual cohorts use beginning-of-year surviving stock. Calibrate the
    # discrete survival sum to the inherited mean; otherwise sampling a
    # continuous Weibull adds approximately half a service year for free.
    lo=0.;hi=2*life/math.gamma(1+1/shape)
    for _ in range(60):
        mid=(lo+hi)/2
        mean=sum(math.exp(-(age/mid)**shape) for age in range(math.ceil(15*life)+1))
        if mean<life:lo=mid
        else:hi=mid
    return (lo+hi)/2


@lru_cache(maxsize=None)
def survival(age,life,shape=3.):
    if age<0:return 1.
    if shape is None:return float(age<life)
    return math.exp(-(age/lifetime_scale(life,shape))**shape)


def cohort_ledger(data,rows,targets,seed):
    """Annual components expire at their recorded age; surplus is not replaced.

    Orders and their manufacture/freight occur two years before service. Seed
    hardware has eight deployment cohorts (short-lived units use their life),
    corresponding to the explicitly assumed preproduction deployment window.
    """
    cats=catalogue(data);start=rows[0]['year'];end=rows[-1]['year'];lead=2
    flows={y:blank_flow() for y in range(start-lead,end+1)};stocks={y:{} for y in range(start,end+1)};events=[]
    shape=data.get('lifetime_shape',3.)
    for name,a in cats.items():
        for part,share,life in [('structure',a['structural_share'],a['structural_life']),('active',1-a['structural_share'],a['life'])]:
            if share<=0:continue
            phases=max(1,min(8,math.ceil(life)))
            denominator=sum(survival(j+1,life,shape) for j in range(phases))
            if denominator<=0:raise ValueError('Seed deployment ages exceed every component lifetime')
            cohorts=[(start-j-1,seed.get(name,0.)*share/denominator) for j in range(phases)]
            curve=[survival(age,life,shape) for age in range(end-start+phases+2)]
            for row in rows:
                y=row['year'];expired=alive=0.;surviving=[]
                for birth,mass in cohorts:
                    age=y-birth;fraction=curve[age]
                    if y>start:expired+=mass*(curve[max(0,age-1)]-fraction)
                    if fraction>1e-12:
                        surviving.append((birth,mass));alive+=mass*fraction
                cohorts=surviving;required=targets[y].get(name,0.)*share
                added=max(0.,required-alive)
                if added>1e-6:
                    cohorts.append((y,added));order=y-lead;replacement=min(added,expired)
                    growth=added-replacement
                    for workyear in (order,order+1):
                        f=flows[workyear];amount=added/2;growth_work=growth/2;replace_work=replacement/2
                        unit=a['price'] if part=='active' or a['civil'] else data['economy']['production']['local_fabrication_value_USD2026_kg']
                        if part=='active':
                            f['imports'][a['destination']]+=amount
                            f['external_capital']+=growth_work*unit;f['external_replacement']+=replace_work*unit
                            f['capital_imports'][a['destination']]+=growth_work
                            if a['owner']=='industry':f['industry_replacement_import_kg']+=replace_work
                            if a['owner'] in ('service','hospitality'):
                                f[a['owner']+'_external_capital']+=growth_work*unit
                                f[a['owner']+'_external_replacement']+=replace_work*unit
                            if name=='return_capsule_buses':f['capsule_purchases_kg']+=amount
                        elif a['material']=='bulk':
                            f['bulk']+=amount;f['capital_bulk']+=growth_work
                        else:
                            f['local'][a['destination']]+=amount;f['capital_local'][a['destination']]+=growth_work
                        if part=='structure':
                            f['local_capital']+=growth_work*unit;f['local_replacement']+=replace_work*unit
                        if a['owner']=='service':
                            f['service_capital']+=growth_work*unit;f['service_replacement']+=replace_work*unit
                        if a['owner']=='hospitality':
                            f['hospitality_capital']+=growth_work*unit;f['hospitality_replacement']+=replace_work*unit
                        if name in ('commercial_civil_bulk','service_civil_bulk'):f['new_floor']+=growth_work/data['economy']['exports']['factory_bulk_kg_m2']
                    events.append(dict(asset=name,component=part,order_year=order,commission_year=y,mass_kg=added,replacement_kg=replacement))
                if expired:
                    f=flows[y];f['retired_mass']+=expired
                    if part=='structure':
                        if a['material']=='bulk':f['bulk_recycled']+=expired*data['industrial']['lunar']['recovery']
                        else:f['recycled'][a['destination']]+=expired*(data['industrial']['lunar']['recovery'] if a['destination']=='surface' else data['economy']['production']['source_recovery'])
                available=alive+added
                old=stocks[y].get(name,math.inf);stocks[y][name]=min(old,available/share)
    previous_beds=rows[0].get('visitor_beds',0.)
    for row in rows:
        beds=row.get('visitor_beds',0.);added=max(0.,beds-previous_beds);previous_beds=beds
        for y in (row['year']-2,row['year']-1):flows[y]['new_visitor_beds']+=added/2
    return flows,stocks,events


def physical(data,row,f,bank=None):
    """Mass balance followed by the inherited chemical freight and process laws."""
    s=data['industrial'];d=data['demand'];z=s['source'];l=s['lunar'];q=data['economy'];ex=q['exports'];rv=data['revenue'];rr=rv['return'];t=s['transport'];p=q['production']
    ps=row['population']/data['normalization']['population'];gs=row.get('new_residential_people',0.)/data['normalization']['population']
    bank=deepcopy(bank or {**{k:0. for k in DESTS},'bulk':0.})
    imports=dict(f['imports']);local=dict(f['local']);virgin={}
    for dest in DESTS:
        bank[dest]+=f['recycled'][dest];used=min(local[dest],bank[dest]);bank[dest]-=used;virgin[dest]=local[dest]-used
    bank['bulk']+=f['bulk_recycled'];rebulk=min(f['bulk'],bank['bulk']);bank['bulk']-=rebulk
    goods=row.get('goods',0.);source_goods=row.get('source_goods',0.);fuel=row.get('propellant_sales',0.);compute=row.get('compute',0.)
    resident_consumables=sum(v for k,v in data['settlement']['operating_import_kg_year'].items() if k not in ('water','air'))*ps
    industrial_consumables=p['process_consumables_fraction_of_imported_renewal']*f.get('industry_replacement_import_kg',0.)
    guest_years=row.get('visitor_nights_year',0.)/JULIAN_YEAR_DAYS
    visitor_consumables=resident_consumables/row['population']*guest_years if row['population'] else 0.
    research_consumables=row.get('service_jobs',0.)*rv['services']['research_consumables_kg_staff_year']
    imports['surface']+=resident_consumables+visitor_consumables+research_consumables
    imports['source']+=industrial_consumables
    factory_feed=goods/ex['factory_yield'];earthfeed=factory_feed*ex['factory_Earth_input_fraction'];imports['surface']+=earthfeed
    factory_local=factory_feed-earthfeed+goods*rv['manufacturing']['transport_packaging_mass_fraction']
    source_product=source_goods/ex['factory_yield']+source_goods*rv['manufacturing']['transport_packaging_mass_fraction']
    surface_alloy=virgin['surface']*p['alloy_import_fraction_of_virgin_local_metal']
    source_alloy=(sum(virgin[k] for k in ('source','depot','leo'))+source_product)*p['alloy_import_fraction_of_virgin_local_metal']
    imports['surface']+=surface_alloy;imports['source']+=source_alloy
    freshbulk=f['bulk']-rebulk+d['goods_virgin']*ps
    freshmetal=virgin['surface']-surface_alloy+factory_local
    bulk=f['bulk']+d['goods_virgin']*ps;metal=local['surface']+factory_local
    raw=freshbulk/l['bulk_yield']+freshmetal/l['metal_yield'];handling=raw+d['loose_and_grading']*gs
    lunar_W=(bulk*l['bulk_kWh_kg']+metal*l['metal_kWh_kg']+handling*l['handling_kWh_kg']+local['surface']*p['local_fabrication_finish_kWh_kg']+goods*ex['factory_finish_kWh_kg'])*1000/HOURS+row.get('service_power_W',0.)+row.get('premium_power_W',0.)+row.get('visitor_beds',0.)*rv['visitors']['mean_W_bed']
    water=d['yearly_water']*(ps+guest_years/data['normalization']['population'])+d['water']*gs+f.get('new_visitor_beds',0.)*rv['visitors']['initial_water_kg_bed']
    air=d['yearly_air']*(ps+guest_years/data['normalization']['population'])+d['air']*gs+f.get('new_visitor_beds',0.)*rv['visitors']['initial_air_kg_bed']
    # New commercial enclosure air follows the explicitly priced floor proxy.
    newfloor=f['new_floor'];air+=newfloor*rv['manufacturing']['air_kg_m2']
    packed=(goods+source_goods)*(1+rv['manufacturing']['transport_packaging_mass_fraction']);source_packed=source_goods*(1+rv['manufacturing']['transport_packaging_mass_fraction'])
    flights=packed/rr['capsule_payload_kg'];bus=rr['capsule_payload_kg']*rr['capsule_dry_to_payload']*(1-rr['shield_fraction_of_dry']);shield=rr['capsule_payload_kg']*rr['capsule_dry_to_payload']*rr['shield_fraction_of_dry']
    deorbit=(rr['capsule_payload_kg']+bus+shield)*math.expm1(rr['deorbit_dv_m_s']/(rr['deorbit_Isp_s']*STANDARD_GRAVITY))
    # Relaunch used buses once per sortie; newly bought buses are already within
    # this stream and must not receive a second launch charge in their purchase year.
    imports['leo']+=flights*(bus+shield+deorbit)-f.get('capsule_purchases_kg',0.)
    recycle_back=f['recycled']['depot']+f['recycled']['leo'];leo_recycle=f['recycled']['leo']
    r=routes(data)
    earth=ferry(imports['surface']+imports['depot']+imports['source']+leo_recycle,packed+local['leo'],r['Earth'])
    moon=ferry(imports['surface']+water+air,packed-source_packed,r['Moon'])
    nitrogen=air*z['air_nitrogen_mass_fraction'];oxygen=air-nitrogen
    depot_water=moon['burned']+oxygen+earth['space_supply']+fuel
    source_cargo=water+nitrogen+depot_water+local['depot']+local['leo']+source_packed
    source=ferry(source_cargo,imports['source']+recycle_back,r['source'])
    source_water=water+depot_water+source['burned'];source_metal=sum(local[k] for k in ('source','depot','leo'))+source_product
    source_virgin=sum(virgin[k] for k in ('source','depot','leo'))+source_product-source_alloy
    feed=max(source_water/(z['water_grade']*z['water_recovery']),nitrogen/(z['nitrogen_grade']*z['nitrogen_recovery']),source_virgin/(z['metal_grade']*z['metal_recovery']))
    process_kWh=feed*z['ore_kWh_kg']+nitrogen*z['nitrogen_kWh_kg']+source_metal*z['metal_refining_and_fabrication_kWh_kg']
    chem=source['burned']/z['water_to_hydrogen_mass_ratio']*z['electrolysis_liquefaction_kWh_kg_H2']
    source_W=(process_kWh+chem)*1000/HOURS+compute
    heat=(process_kWh+chem*(1-z['hydrogen_HHV_kWh_kg']/z['electrolysis_liquefaction_kWh_kg_H2']))*1000/HOURS+compute
    depot_W=depot_water/z['water_to_hydrogen_mass_ratio']*z['electrolysis_liquefaction_kWh_kg_H2']*1000/HOURS
    ground=d['occupied_W']*ps+lunar_W*s['economics']['capacity_margin']/l['process_availability']
    operating_material=(resident_consumables+visitor_consumables+research_consumables+industrial_consumables+surface_alloy+source_alloy+earthfeed)*s['prices']['other_special_USD2026_kg']
    capsule_op=flights*(shield*s['prices']['spacecraft_USD2026_kg']+deorbit*t['earth_propellant_USD2026_kg'])+packed*rr['Earth_ground_service_USD2026_per_payload_kg']
    LEO=sum(imports.values())
    return dict(lunar_raw_kg_year=raw,lunar_bulk_kg_year=bulk,lunar_metal_kg_year=metal,lunar_recycled_bulk_kg_year=rebulk,
        lunar_virgin_metal_kg_year=freshmetal,lunar_handling_kg_year=handling,lunar_process_W=lunar_W,
        source_feed_kg_year=feed,source_water_kg_year=source_water,source_nitrogen_kg_year=nitrogen,source_virgin_metal_kg_year=source_virgin,
        source_metal_kg_year=source_metal,source_process_W=source_W,source_heat_W=heat,source_feed_process_W=feed*z['ore_kWh_kg']*1000/HOURS,
        source_product_refining_W=(nitrogen*z['nitrogen_kWh_kg']+source_metal*z['metal_refining_and_fabrication_kWh_kg'])*1000/HOURS,
        depot_W=depot_W,depot_water_kg_year=depot_water,ground_firm_W=ground,Earth=earth,Moon=moon,source=source,
        source_cargo_kg_year=source_cargo,
        source_cargo_breakdown_kg_year=dict(resident_and_guest_water=water,nitrogen=nitrogen,depot_water=depot_water,local_offsite_structures=local['depot']+local['leo'],source_exports=source_packed),
        depot_water_breakdown_kg_year=dict(lunar_ferry_fuel=moon['burned'],oxygen_makeup=oxygen,Earth_ferry_fuel=earth['space_supply'],external_fuel_sales=fuel),
        source_water_breakdown_kg_year=dict(resident_and_guest_water=water,depot_water=depot_water,source_ferry_fuel=source['burned']),
        source_retained_tails_kg_year=feed-source_water-nitrogen-source_virgin,LEO_kg_year=LEO,imports_by_destination_kg_year=imports,
        local_fabrication_by_destination_kg_year=local,recycled_inventory_kg=bank,
        factory_floor_new_m2=newfloor,construction_floor_m2=data['settlement']['production_demand']['floor']*gs+newfloor,
        operating_material_USD2026=operating_material,capsule_operations_USD2026=capsule_op,Earth_launch_USD2026=LEO*t['LEO_USD2026_kg'],
        resource_provision_USD2026=raw*s['economics']['local_restoration_USD2026_kg_feed']+feed*rv['finance']['source_transition_USD2026_per_kg_feed'])


def requirements(data,row,f):
    """Rated machinery required by actual flows, with inherited spare margins."""
    s=data['industrial'];l=s['lunar'];z=s['source'];p=s['power'];t=s['transport'];d=data['demand'];margin=s['economics']['capacity_margin'];original=data['original_architecture']
    out=residential_targets(data,row['population']);out.update(commercial_targets(data,row))
    def machine(name,flow,rate,unit,availability):out[name]=math.ceil(flow*margin/(rate*HOURS*availability))*unit
    machine('excavators',f['lunar_handling_kg_year'],l['excavator_kg_hr'],l['excavator_kg'],l['mine_availability'])
    machine('haulers',f['lunar_handling_kg_year'],l['hauler_payload_kg']/l['haul_roundtrip_hr'],l['hauler_kg'],l['mine_availability'])
    for prefix,key in [('', 'lunar_metal_kg_year'),('source_', 'source_metal_kg_year')]:
        machine(prefix+'metal_trains',f[key],l['metal_train_kg_hr'],l['metal_train_hardware_kg'],l['process_availability'])
        machine(prefix+'fabrication_cells',f[key],l['fabrication_cell_kg_hr'],l['fabrication_cell_hardware_kg'],l['process_availability'])
    machine('bulk_trains',f['lunar_bulk_kg_year'],l['bulk_train_kg_hr'],l['bulk_train_hardware_kg'],l['process_availability'])
    machine('source_mine_trains',f['source_feed_kg_year'],z['mine_module_kg_hr'],z['mine_module_hardware_kg'],z['mine_availability'])
    out['construction_and_service_robots']=l['construction_plant_kg']*f['construction_floor_m2']/data['settlement']['production']['floor_m2_year']
    ground=f['ground_firm_W'];industrialfirm=f['lunar_process_W']*margin/l['process_availability']
    out['additional_ground_receivers']=industrialfirm*d['power_hardware_kg_W']
    area=industrialfirm*d['old_radiator_m2']/d['old_total_W'];out['additional_ground_radiator_controls']=area*d['radiator_controls_kg_m2']
    sourceDC=ground/(p['RF_efficiency']*p['beam_capture']*p['receiver_efficiency'])*p['stations']/2/p['end_of_life_factor']
    out['lunar_beam_constellation']=sourceDC/p['space_system_W_kg']
    radflux=p['emissivity']*STEFAN_BOLTZMANN*p['radiator_K']**4
    for pre,W,heat in [('source',f['source_process_W'],f['source_heat_W']),('depot',f['depot_W'],f['depot_W']*(1-z['hydrogen_HHV_kWh_kg']/z['electrolysis_liquefaction_kWh_kg_H2']))]:
        chemical_W=max(0.,W-row.get('compute',0.)) if pre=='source' else W
        out[pre+'_chemistry']=chemical_W*margin*z['chemistry_hardware_kg_W']
        out[pre+'_solar']=W*margin/p['end_of_life_factor']/(p['source_solar_flux_fraction'] if pre=='source' else 1)/p['space_system_W_kg']
        out[pre+'_radiators']=heat*margin/radflux*p['radiator_kg_m2']
    r=routes(data)
    counts={name:math.ceil(f[name]['trips']*margin/(a['cycles']*t['fleet_availability'])) for name,a in r.items()}
    out['lunar_landers']=counts['Moon']*r['Moon']['dry'];out['asteroid_tugs']=counts['source']*r['source']['dry']
    for name,a in original['assets'].items():
        if name.startswith('Earth_'):out[name]=counts['Earth']*a['mass_kg']/original['engineering']['counts']['Earth_route_ships']
    out['depots_and_orbital_workshops']=(f['source_cargo_kg_year']+f['Moon']['burned'])*t['depot_inventory_years']*t['tank_dry_to_contents']
    for name,driver in [('LEO_route_tanks','out'),('LEO_route_storage_controls','out'),('depot_route_storage_controls','space_supply')]:
        peak=original['peak_flows'];oldstock=peak['Earth_route_LEO_inventory_kg'] if driver=='out' else peak['Earth_route_depot_inventory_kg']
        out[name]=original['assets'][name]['mass_kg']*f['Earth'][driver]*data['transport']['depots']['inventory_years']/oldstock
    surface=sum(v for n,v in out.items() if n in data['templates'] and data['templates'][n]['owner']=='industry' and data['templates'][n]['destination']=='surface')
    out['industrial_civil_bulk']=surface*l['infrastructure_bulk_per_hardware_kg']
    out['industrial_civil_metal']=surface*l['infrastructure_metal_per_hardware_kg']+area*d['radiator_metal_kg_m2']
    return out


def seed_stock(data,first,bootstrap_years=8):
    """Size one seed against construction throughput, then retain its assets.

    Complete seed industrial equipment is procured from outside. This avoids
    assuming an unbuilt local industry manufactures its own first active plant.
    The eight-year average is a throughput screen, not a qualified deployment.
    """
    cats=catalogue(data);base=residential_targets(data,first['population']);base.update(commercial_targets(data,first));seed=dict(base)
    row=dict(first,new_residential_people=first['population']/bootstrap_years,population=first['population']/2)
    def demands(stock):
        f=blank_flow()
        for n,mass in stock.items():
            a=cats[n];share=a['structural_share'];amount=mass/bootstrap_years
            def gross_factor(life):
                phases=max(1,min(8,math.ceil(life)));average=sum(survival(j+1,life,data.get('lifetime_shape',3.)) for j in range(phases))/phases
                return 1/average
            structural_factor=gross_factor(a['structural_life']);active_factor=gross_factor(a['life'])
            if a['owner']=='industry' and not a['civil']:
                amount*=share*structural_factor+(1-share)*active_factor
                f['imports'][a['destination']]+=amount
                f['external_capital']+=amount*a['price']
            else:
                imported=amount*(1-share)*active_factor;made=amount*share*structural_factor
                f['imports'][a['destination']]+=imported;f['external_capital']+=imported*a['price']
                if a['material']=='bulk':f['bulk']+=made
                else:f['local'][a['destination']]+=made
                unit=a['price'] if a['civil'] else data['economy']['production']['local_fabrication_value_USD2026_kg']
                f['local_capital']+=made*unit
                if a['owner']=='service':
                    f['service_capital']+=imported*a['price']+made*unit;f['service_external_capital']+=imported*a['price']
        ex=data['economy']['exports'];sv=data['revenue']['services']
        f['new_floor']=(first.get('goods_capacity',first.get('goods',0.))/(ex['factory_yield']*ex['factory_utilization'])*ex['factory_floor_m2_per_nameplate_kg_year']+first.get('service_jobs',0.)*sv['research_floor_m2_staff'])/bootstrap_years
        f['new_visitor_beds']=first.get('visitor_beds',0.)/bootstrap_years
        return f
    history=[]
    for iteration in range(300):
        flow=demands(seed);phys=physical(data,row,flow);needed=requirements(data,row,phys);needed.update(base)
        # Seed also supports opening operation and the first construction order.
        opening=requirements(data,first,physical(data,first,blank_flow()))
        needed={k:max(v,opening.get(k,0.)) for k,v in needed.items()}
        error=max((max(0.,v-seed.get(k,0.))/max(1.,v) for k,v in needed.items()),default=0.)
        history.append(error)
        seed={k:max(seed.get(k,0.),v) for k,v in needed.items()}
        if error<1e-9:break
    else:raise ValueError('Seed stock/flow requirements did not converge')
    flow=demands(seed);phys=physical(data,row,flow);I=data['normalization']['USD2026_to_USD2041']
    original_value=sum(a['mass']*a['price'] for a in cats.values() if a['owner']=='industry' and not a['civil'])
    seed_value=sum(seed.get(n,0.)*a['price'] for n,a in cats.items() if a['owner']=='industry' and not a['civil'])
    ratio=seed_value/original_value;econ=data['industrial']['economics'];staff=econ['staff']*ratio
    resident_jobs=staff/4+(data['settlement']['food_operators']+data['settlement']['industry_operators'])*first['population']/data['normalization']['population']/2
    resident_payroll=resident_jobs*econ['loaded_staff_USD2026_year']*bootstrap_years*I
    external_payroll=staff/4*econ['loaded_staff_USD2026_year']*bootstrap_years*I
    cost=dict(external_hardware_USD2041=flow['external_capital']*bootstrap_years*I,
        resident_payroll_USD2041=resident_payroll,external_payroll_USD2041=external_payroll,
        external_engineering_support_USD2041=data['economy']['production']['external_engineering_contract_USD2026_year']*ratio/2*bootstrap_years*I,
        delivery_and_operating_material_USD2041=(phys['Earth_launch_USD2026']+phys['operating_material_USD2026']+phys['capsule_operations_USD2026'])*bootstrap_years*I,
        resource_provision_USD2041=phys['resource_provision_USD2026']*bootstrap_years*I,
        common_programme_USD2041=data['common_programme_USD2041'],
        identified_service_hardware_USD2041=flow['service_external_capital']*bootstrap_years*I)
    return seed,dict(stock_mass_kg=seed,bootstrap_years=bootstrap_years,cost=cost,resident_payroll_USD2041=resident_payroll,local_fabrication_account_value_USD2041=flow['local_capital']*bootstrap_years*I,
        lunar_raw_kg=phys['lunar_raw_kg_year']*bootstrap_years,source_feed_kg=phys['source_feed_kg_year']*bootstrap_years,
        physical_annual_average=phys,iterations=iteration+1,maximum_relative_requirement_residual=history[-1],
        commission_years=list(range(first['year']-8,first['year'])),
        evidence='Seed cohorts commissioned in the eight years before opening; gross procurement includes survival losses to opening. Complete industrial equipment imported once; local enclosure/civil construction and average half-occupied seed operations included. Detailed preproduction staging, staff accommodation and source qualification remain unsized.')


def solve(data,annual,*,max_iterations=100,tolerance=1e-7):
    """Close dated component procurements and demand-sized productive capacity.

    Capacity plans increase monotonically within the fixed-point solution.
    Each year's target is independent of the population multiplier and can
    decline after construction demand subsides. Surviving surplus stays on hand.
    """
    rows=deepcopy(annual)
    if not rows or any(b['year']!=a['year']+1 for a,b in zip(rows,rows[1:])):raise ValueError('Consecutive annual request rows required')
    if any(row.get(k,0)<0 for row in rows for k in ('population','goods','source_goods','compute','propellant_sales','new_residential_people')):raise ValueError('Negative physical demand')
    # Two forward years prevent final observations receiving free expired plant.
    end=rows[-1]['year']
    for y in range(end+1,end+41):rows.append(dict(rows[-1],year=y,new_residential_people=0.,commissioned_people=0.))
    seed,seed_record=seed_stock(data,rows[0]);targets={r['year']:requirements(data,r,physical(data,r,blank_flow())) for r in rows}
    residuals=[]
    for iteration in range(max_iterations):
        flows,stocks,events=cohort_ledger(data,rows,targets,seed)
        bank=None;maximum=0.;next_targets={}
        for index,row in enumerate(rows):
            workrow=dict(row,new_residential_people=(row.get('new_residential_people',0.)+(rows[index-1].get('new_residential_people',0.) if index else 0.))/2)
            y=row['year'];f=physical(data,workrow,flows[y],bank);bank=f['recycled_inventory_kg'];need=requirements(data,row,f)
            previous=targets[y]
            error=max((max(0.,v-previous.get(k,0.))/max(1.,v) for k,v in need.items()),default=0.)
            maximum=max(maximum,error)
            next_targets[y]={k:max(v,previous.get(k,0.)) for k,v in need.items()}
        residuals.append(maximum);targets=next_targets
        if maximum<tolerance:break
    converged=residuals[-1]<tolerance
    flows,stocks,events=cohort_ledger(data,rows,targets,seed);bank=None;result=[];I=data['normalization']['USD2026_to_USD2041'];cats=catalogue(data)
    original_value=sum(a['mass']*a['price'] for a in cats.values() if a['owner']=='industry' and not a['civil'])
    cumulative_lunar=cumulative_source=0.
    for index,row in enumerate(rows):
        y=row['year'];f=flows[y];opening_bank=deepcopy(bank)
        workrow=dict(row,new_residential_people=(row.get('new_residential_people',0.)+(rows[index-1].get('new_residential_people',0.) if index else 0.))/2)
        physical_row=physical(data,workrow,f,bank);bank=physical_row['recycled_inventory_kg'];need=requirements(data,row,physical_row)
        ratios={n:need[n]/max(1e-30,stocks[y].get(n,0.)) for n in need if need[n]>1e-8}
        # Capacity comparison uses required machinery including the inherited25%
        # margin. A ratio>1 is an unsatisfied requirement, never added capacity.
        maximum=max(ratios.values(),default=0.);feasible=converged and maximum<=1+1e-6
        industrial_value=sum(min(stocks[y].get(n,0.),need.get(n,0.))*a['price'] for n,a in cats.items() if a['owner']=='industry' and not a['civil'])
        staffing=data['industrial']['economics']['staff']*industrial_value/original_value
        ex=data['economy']['exports'];rv=data['revenue'];P0=data['normalization']['population']
        factory_jobs=row.get('goods_capacity',row.get('goods',0.))/(ex['factory_yield']*ex['factory_utilization'])/(P0*rv['manufacturing']['nameplate_finished_kg_per_resident_year'])*rv['manufacturing']['staff']
        compute_jobs=row.get('compute_capacity',row.get('compute',0.))/1e6*rv['computing']['staff_per_MW']
        base_jobs=(data['settlement']['food_operators']+data['settlement']['industry_operators'])*row['population']/P0
        core_resident_jobs=staffing/2+base_jobs+factory_jobs+compute_jobs
        commercial_jobs=row.get('service_jobs',0.)+row.get('premium_jobs',0.)+row.get('logistics_jobs',0.)+row.get('visitor_beds',0.)*rv['visitors']['staff_per_bed']
        total_resident_jobs=core_resident_jobs+commercial_jobs
        ratios['resident_workforce']=total_resident_jobs/max(1.,row['population']*.4)
        maximum=max(ratios.values(),default=0.);feasible=converged and maximum<=1+1e-6
        payroll=core_resident_jobs*data['industrial']['economics']['loaded_staff_USD2026_year']
        external_staff=staffing/2*data['industrial']['economics']['loaded_staff_USD2026_year']+data['economy']['production']['external_engineering_contract_USD2026_year']*industrial_value/original_value
        # Evaluate the same physical system without new capital manufacture.
        # Replacement, occupied residents, sold outputs and recycling stay in it.
        opflow=deepcopy(f)
        for dest in DESTS:
            opflow['imports'][dest]-=f['capital_imports'][dest]
            opflow['local'][dest]-=f['capital_local'][dest]
        opflow['bulk']-=f['capital_bulk'];opflow['new_floor']=0.;opflow['new_visitor_beds']=0.
        # Capsule purchase masses travel in the sortie stream; subtract only
        # replacement purchases when forming the operating counterfactual.
        capbus=sum(ev['mass_kg']-ev['replacement_kg'] for ev in events if ev['asset']=='return_capsule_buses' and ev['order_year'] in (y,y-1))/2
        opflow['capsule_purchases_kg']=max(0.,f['capsule_purchases_kg']-capbus)
        op_physical=physical(data,dict(workrow,new_residential_people=0.),opflow,opening_bank)
        variable=('Earth_launch_USD2026','operating_material_USD2026','capsule_operations_USD2026')
        construction_logistics=sum(physical_row[k]-op_physical[k] for k in variable)
        cap_resource=physical_row['resource_provision_USD2026']-op_physical['resource_provision_USD2026']
        production_value=(row.get('goods',0.)+row.get('source_goods',0.))*data['economy']['production']['structural_metal_value_USD2026_kg']
        construction_fraction=f['local_capital']/max(1.,f['local_capital']+f['local_replacement']+production_value)
        construction_payroll=staffing/2*data['industrial']['economics']['loaded_staff_USD2026_year']*construction_fraction
        costs=dict(capital_external_USD2041=(f['external_capital']+construction_logistics)*I,
            capital_local_USD2041=construction_payroll*I,
            local_capital_account_value_USD2041=f['local_capital']*I,local_replacement_account_value_USD2041=f['local_replacement']*I,
            capital_resource_provision_USD2041=cap_resource*I,
            operating_external_replacement_USD2041=f['external_replacement']*I,
            operating_external_logistics_materials_USD2041=sum(op_physical[k] for k in variable)*I,
            operating_external_services_USD2041=external_staff*I,
            operating_local_replacement_USD2041=0.,
            operating_local_payroll_USD2041=(payroll-construction_payroll)*I,
            resident_payroll_USD2041=payroll*I,
            operating_resource_provision_USD2041=op_physical['resource_provision_USD2026']*I,
            resource_provision_USD2041=physical_row['resource_provision_USD2026']*I,
            identified_service_capital_USD2041=(f['service_external_capital']+construction_payroll*max(0.,f['service_capital']-f['service_external_capital'])/max(1.,f['local_capital']))*I,
            identified_service_replacement_USD2041=f['service_external_replacement']*I,
            identified_hospitality_capital_USD2041=f['hospitality_external_capital']*I,
            identified_hospitality_replacement_USD2041=f['hospitality_external_replacement']*I,
            Earth_launch_USD2041=physical_row['Earth_launch_USD2026']*I,
            operating_material_USD2041=physical_row['operating_material_USD2026']*I,
            capsule_operations_USD2041=physical_row['capsule_operations_USD2026']*I)
        service_stock=service_mean_renewal=0.
        for n,mass in commercial_targets(data,row).items():
            a=cats[n]
            if a['owner']!='service':continue
            share=a['structural_share'];localprice=a['price'] if a['civil'] else data['economy']['production']['local_fabrication_value_USD2026_kg']
            service_stock+=mass*((1-share)*a['price']+share*localprice)
            service_mean_renewal+=mass*((1-share)*a['price']/a['life']+share*localprice/a['structural_life'])
        costs['identified_service_stock_value_USD2041']=service_stock*I
        costs['identified_service_mean_renewal_allowance_USD2041']=service_mean_renewal*I
        costs['operating_total_USD2041']=sum(costs[k] for k in ('operating_external_replacement_USD2041','operating_external_logistics_materials_USD2041','operating_external_services_USD2041','operating_local_replacement_USD2041','operating_local_payroll_USD2041','operating_resource_provision_USD2041'))
        costs['capital_total_USD2041']=costs['capital_external_USD2041']+costs['capital_local_USD2041']+costs['capital_resource_provision_USD2041']
        direct_total=(f['external_capital']+f['external_replacement']+sum(physical_row[k] for k in variable)+external_staff+payroll+physical_row['resource_provision_USD2026'])*I
        costs['accounting_residual_USD2041']=costs['capital_total_USD2041']+costs['operating_total_USD2041']-direct_total
        yearly_events=[ev for ev in events if ev['order_year'] in (y,y-1)]
        drivers={}
        for ev in yearly_events:
            n=ev['asset'];a=cats[n];unit=a['price'] if ev['component']=='active' or a['civil'] else data['economy']['production']['local_fabrication_value_USD2026_kg']
            if ev['component']=='active':drivers[n]=drivers.get(n,0.)+ev['mass_kg']/2*unit*I
        tug_import=sum(ev['mass_kg']/2 for ev in yearly_events if ev['asset']=='asteroid_tugs' and ev['component']=='active')
        physical_row['source_tug_component_import_kg_year']=tug_import
        physical_row['source_tug_count']=stocks[y].get('asteroid_tugs',0.)/routes(data)['source']['dry']
        cumulative_lunar+=physical_row['lunar_raw_kg_year'];cumulative_source+=physical_row['source_feed_kg_year']
        if y<=end:
            result.append(dict(year=y,population=row['population'],request=row,physical=physical_row,cash=costs,
                industrial_staff_total=staffing,technical_resident_jobs=total_resident_jobs,core_resident_jobs=core_resident_jobs,commercial_resident_jobs=commercial_jobs,
                top_imported_component_purchases_USD2041=dict(sorted(drivers.items(),key=lambda kv:kv[1],reverse=True)[:8]),
                required_asset_mass_kg=need,available_asset_mass_kg=stocks[y],capacity_use=ratios,
                maximum_capacity_use=maximum,feasible=feasible,
                actual_sales_capacity_fraction=1. if feasible else 0.,
                shortfalls={k:v-1 for k,v in ratios.items() if v>1+1e-6},
                cumulative_lunar_raw_kg=cumulative_lunar,cumulative_source_feed_kg=cumulative_source))
    # Procurements before opening belong to the seed once, not year8 operations.
    late=blank_flow()
    for y in range(rows[0]['year']-2,rows[0]['year']):
        f=flows[y]
        for dest in DESTS:
            late['imports'][dest]+=f['imports'][dest];late['local'][dest]+=f['local'][dest]
        for k in ('bulk','external_capital','external_replacement','local_capital','local_replacement','service_capital','service_replacement','new_floor','new_visitor_beds','capsule_purchases_kg'):late[k]+=f[k]
    late_phys=physical(data,dict(rows[0],population=0.,new_residential_people=0.,goods=0.,source_goods=0.,propellant_sales=0.,compute=0.,service_jobs=0.,visitor_beds=0.,visitor_nights_year=0.,service_power_W=0.,premium_power_W=0.),late)
    seed_record['additional_preproduction_cost_USD2041']=(late['external_capital']+late['external_replacement']+late_phys['Earth_launch_USD2026']+late_phys['operating_material_USD2026']+late_phys['resource_provision_USD2026'])*I
    seed_record['lunar_raw_kg']+=late_phys['lunar_raw_kg_year'];seed_record['source_feed_kg']+=late_phys['source_feed_kg_year']
    seed_record['cost']['additional_preproduction_cost_USD2041']=seed_record['additional_preproduction_cost_USD2041']
    seed_record['total_resource_provision_USD2041']=seed_record['cost']['resource_provision_USD2041']+late_phys['resource_provision_USD2026']*I
    seed_record['cost']['total_USD2041']=sum(v for k,v in seed_record['cost'].items() if k!='identified_service_hardware_USD2041')
    return dict(seed=seed_record,annual=result,component_orders=[r for r in events if r['order_year']<=end],
        convergence=dict(converged=converged,iterations=iteration+1,relative_requirement_residual=residuals[-1],history=residuals),
        first_capacity_shortfall=next((r['year'] for r in result if not r['feasible']),None))
