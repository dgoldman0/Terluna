"""Compare freight architectures without changing the initial industrial baseline.

The iteration includes each alternative's imported equipment and replacement,
source extraction, power and fuel. Delta-v and trip times are design envelopes,
not integrated trajectories. The cases do not grant extraction permission.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import json
import math
from pathlib import Path
from . import run as base
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, STANDARD_GRAVITY, STEFAN_BOLTZMANN

HERE, ROOT = base.HERE, base.ROOT
SCHEMA = 'terluna.research.industrial-transport/1'


def read_inputs(spec=None):
    q = deepcopy(spec) if spec is not None else json.loads((HERE/'transport_inputs.json').read_text())
    for name, record in q['inputs'].items():
        if base.sha(ROOT/record['path']) != record['sha256']:
            raise ValueError(f'Input hash mismatch: {name}')
    s, src = base.read_inputs()
    return q, s, src


def depot_fuelled_burns(cargo, dry, dv, isp, margin=1, retention=1, inventory_fraction=0):
    """LEO -> depot cargo; depot -> LEO fuel for the following cargo trip.

A is ascent fuel at LEO. A/retention is the fuel cargo loaded at the depot.
B brakes that fuel cargo and the tug into LEO. No aerocapture is credited.
"""
    if min(cargo, dry, dv, inventory_fraction) < 0 or isp <= 0 or margin < 1 or not 0 < retention <= 1:
        raise ValueError('Invalid depot-fuelled route')
    k = math.expm1(dv/(isp*STANDARD_GRAVITY))*margin
    ascent = (cargo+dry)*k
    backhaul = ascent*(1+inventory_fraction)/retention
    descent = (dry+backhaul)*k
    return dict(ascent_kg=ascent,backhaul_kg=backhaul,descent_kg=descent,
                supplied_kg=(backhaul+descent)/retention,
                peak_tank_contents_kg=backhaul+descent,
                inventory_delivery_kg=ascent*inventory_fraction,
                handling_loss_kg=(ascent*(1+inventory_fraction)+backhaul+descent)*(1/retention-1))


def ship(s, q, earth=True, source_return_ratio=0, installed_dry_kg=None, construction=True):
    """Close vehicle mass against tanks, electricity, heat rejection and drives."""
    t, p = s['transport'], s['power']
    mode = q['earth_route'] if earth else q['source_route']
    cargo = t['earth_payload_kg'] if earth else t['nea_payload_kg']
    electric = mode == 'electric'
    c = q['electric'] if electric else q['chemical']
    dv = (c['earth_dv_oneway_m_s'] if earth else c['source_dv_oneway_m_s']) if electric else t['earth_dv_oneway_m_s']
    isp = c['Isp_s'] if electric else t['Isp_s']
    days = c['earth_powered_days'] if earth and electric else c.get('source_powered_days',0)
    cycle = (c['earth_cycle_days'] if earth else c['source_cycle_days']) if electric else JULIAN_YEAR_DAYS/t['earth_cycles_year']
    if electric and (not 0 < c['total_efficiency'] < 1 or cycle < days or days <= 0):
        raise ValueError('Invalid electric efficiency or mission duration')
    flux = p['emissivity']*STEFAN_BOLTZMANN*p['radiator_K']**4
    dry = cargo*(c['structural_dry_to_payload'] if electric else c['fixed_dry_to_payload'])
    if installed_dry_kg is not None:dry=installed_dry_kg
    for iteration in range(300):
        if earth:
            inventory_fraction=q['depots']['inventory_years']/s['schedule']['delivery_equivalent_years'] if construction else 0
            burns = depot_fuelled_burns(cargo,dry,dv,isp,t['propellant_margin'],c['handling_retention'],inventory_fraction)
            prop = burns['ascent_kg']+burns['descent_kg']
        else:
            r = base.round_trip(dv,isp,dry/cargo,source_return_ratio)
            prop = cargo*r['total']*t['propellant_margin']
            burns = dict(ascent_kg=cargo*r['outbound']*t['propellant_margin'],
                         descent_kg=cargo*r['returning']*t['propellant_margin'],
                         supplied_kg=prop/c['handling_retention'],peak_tank_contents_kg=prop,
                         handling_loss_kg=prop/c['handling_retention']-prop,
                         return_hardware_kg=cargo*source_return_ratio)
        parts = {}
        def part(name,mass,price,life):
            parts[name]=dict(mass_kg=mass,manufacture_USD2026=price,life_years=life)
        hull = cargo*(c['structural_dry_to_payload'] if electric else c['fixed_dry_to_payload'])
        tanks = burns['peak_tank_contents_kg']*(c['tank_dry_to_propellant'] if electric else c['tank_dry_to_peak_propellant'])
        life = c['hull_life_years'] if electric else t['earth_tug_life_years']
        part('hull_tanks',hull+tanks,(hull+tanks)*s['prices']['spacecraft_USD2026_kg'],life)
        if electric:
            ve = isp*STANDARD_GRAVITY
            W = prop*ve**2/(2*c['total_efficiency']*days*JULIAN_DAY)
            heat = W*(1-c['total_efficiency'])
            drive_life = min(life,c['drive_powered_life_hours']/(days*24*JULIAN_YEAR_DAYS/cycle*t['fleet_availability']))
            part('drive',W*c['drive_kg_W'],W*c['drive_USD2026_W'],drive_life)
            array_life = c['earth_array_life_years'] if earth else c['source_array_life_years']
            rad_life = c['radiator_life_years']
        else:
            W = burns['peak_tank_contents_kg']*c['cooling_electric_W_per_stored_kg']
            heat = W
            part('cryocooling',W*c['cooling_hardware_kg_W'],W*c['cooling_hardware_USD2026_W'],10)
            array_life,rad_life = p['space_life_years'],20
        array_W = W/p['end_of_life_factor']/(1 if earth else p['source_solar_flux_fraction'])
        part('solar',array_W/p['space_system_W_kg'],array_W*p['space_USD2026_W'],array_life)
        area = heat/flux
        part('radiators',area*p['radiator_kg_m2'],area*p['radiator_kg_m2']*s['prices']['lunar_hardware_USD2026_kg'],rad_life)
        new = sum(a['mass_kg'] for a in parts.values())
        if not math.isfinite(new) or new > 100*cargo:
            raise ValueError('Vehicle mass/power iteration diverges within this trip-time envelope')
        if installed_dry_kg is not None or abs(new-dry) < 1e-10*max(new,1):
            return dict(cargo_kg=cargo,dry_kg=dry if installed_dry_kg is not None else new,burns=burns,parts=parts,
                        electric_bus_W=W,solar_nameplate_W=array_W,radiator_emitting_m2=area,
                        cycle_days=cycle,powered_days=days,delta_v_oneway_m_s=dv,Isp_s=isp,
                        thrust_N=prop*isp*STANDARD_GRAVITY/(days*JULIAN_DAY) if electric else None,
                        source_return_ratio=source_return_ratio,iterations=iteration+1)
        dry = new
    raise ValueError('Vehicle mass did not converge')


def split_ship(s,q):
    """Existing chemical hull; each end supplies its departing burn."""
    t=s['transport'];cargo=t['earth_payload_kg'];dry=cargo*t['earth_dry_to_payload']
    k=math.expm1(t['earth_dv_oneway_m_s']/(t['Isp_s']*STANDARD_GRAVITY))*t['propellant_margin']
    ascent=(cargo+dry)*k;descent=dry*k;ret=q['chemical']['handling_retention']
    return dict(cargo_kg=cargo,dry_kg=dry,burns={'ascent_kg':ascent,'descent_kg':descent,'backhaul_kg':0,
        'supplied_kg':descent/ret,'handling_loss_kg':descent*(1/ret-1),'inventory_delivery_kg':0},
        parts={'hull_tanks':{'mass_kg':dry,'manufacture_USD2026':dry*s['prices']['spacecraft_USD2026_kg'],'life_years':t['earth_tug_life_years']}},
        electric_bus_W=0,solar_nameplate_W=0,radiator_emitting_m2=0,cycle_days=JULIAN_YEAR_DAYS/t['earth_cycles_year'],
        powered_days=0,delta_v_oneway_m_s=t['earth_dv_oneway_m_s'],Isp_s=t['Isp_s'],source_return_ratio=0)


def flows(s,d,assets,q,construction=True):
    f = base.flows(s,d,assets,construction)
    if q['earth_route']=='earth_chemical':
        return f
    p,z,t = s['power'],s['source'],s['transport']
    split=q['earth_route']=='split_chemical'
    v = split_ship(s,q) if split else ship(s,q)
    if not construction and not split:v=ship(s,q,installed_dry_kg=v['dry_kg'],construction=False)
    earth_trips = f['earth_deep_kg_year']/v['cargo_kg']
    normal=depot_fuelled_burns(v['cargo_kg'],v['dry_kg'],v['delta_v_oneway_m_s'],v['Isp_s'],t['propellant_margin'],(q['chemical'] if q['earth_route']=='space_chemical' else q['electric'])['handling_retention'])
    if split:normal=v['burns']
    depot_stock=earth_trips*normal['supplied_kg']*q['depots']['inventory_years']
    depot_fill=depot_stock/s['schedule']['delivery_equivalent_years'] if construction else 0
    earth_feed = earth_trips*v['burns']['supplied_kg']+depot_fill
    earth_loss = earth_trips*v['burns']['handling_loss_kg']
    f['earth_launch_propellant_kg_year']=earth_trips*v['burns']['ascent_kg'] if split else 0.0
    f['earth_transfer_propellant_kg_year']=earth_trips*(v['burns']['ascent_kg']+v['burns']['descent_kg'])
    f['Earth_route_source_water_kg_year']=earth_feed
    f['Earth_route_space_fuel_consumed_kg_year']=f['earth_transfer_propellant_kg_year']-f['earth_launch_propellant_kg_year']
    f['Earth_route_handling_loss_kg_year']=earth_loss
    f['Earth_route_inventory_fill_kg_year']=earth_trips*v['burns']['inventory_delivery_kg']+depot_fill
    f['Earth_route_LEO_inventory_kg']=0 if split else earth_trips*v['burns']['ascent_kg']*q['depots']['inventory_years']
    f['Earth_route_depot_inventory_kg']=depot_stock
    f['Earth_route_fuel_backhaul_kg_year']=earth_trips*v['burns']['backhaul_kg']
    f['LEO_kg_year']=f['earth_deep_kg_year']+f['payload_destination_kg_year']['leo']+f['earth_launch_propellant_kg_year']
    f['nea_payload_kg_year']+=earth_feed
    f['Earth_vehicle']=v
    if q['source_route']=='electric':
        installed_dry=None
        if not construction:
            count=assets['Source_hull_tanks']['vehicle_count']
            installed_dry=sum(a['mass_kg'] for k,a in assets.items() if k.startswith('Source_'))/count
        source_v=ship(s,q,False,f['nea_return_hardware_kg_year']/f['nea_payload_kg_year'],installed_dry)
        trips=f['nea_payload_kg_year']/source_v['cargo_kg']
        f['nea_propellant_kg_year']=trips*(source_v['burns']['ascent_kg']+source_v['burns']['descent_kg'])
        source_water=trips*source_v['burns']['supplied_kg']
        f['Source_vehicle']=source_v
        source_chem_kWh=0
        f['nea_round_trip']['mass_ratio']=math.exp(q['electric']['source_dv_oneway_m_s']/(q['electric']['Isp_s']*STANDARD_GRAVITY))
    else:
        nr=base.round_trip(t['nea_dv_oneway_m_s'],t['Isp_s'],t['nea_dry_to_payload'],f['nea_return_hardware_kg_year']/f['nea_payload_kg_year'])
        f['nea_round_trip']=nr
        f['nea_propellant_kg_year']=f['nea_payload_kg_year']*nr['total']*t['propellant_margin']
        source_water=f['nea_propellant_kg_year']
        source_chem_kWh=source_water/z['water_to_hydrogen_mass_ratio']*z['electrolysis_liquefaction_kWh_kg_H2']
    f['Source_route_water_kg_year']=source_water
    f['source_water_kg_year']=f['moon_water_kg_year']+f['depot_water_for_propellant_kg_year']+earth_feed+source_water
    f['source_ore_water_bound_kg_year']=f['source_water_kg_year']/(z['water_grade']*z['water_recovery'])
    f['source_ore_kg_year']=max(f['source_ore_water_bound_kg_year'],f['source_ore_nitrogen_bound_kg_year'],f['source_ore_metal_bound_kg_year'])
    f['source_retained_tails_kg_year']=f['source_ore_kg_year']-f['source_water_kg_year']-f['source_nitrogen_kg_year']-f['source_virgin_metal_kg_year']
    process_kWh=f['source_ore_kg_year']*z['ore_kWh_kg']+f['source_nitrogen_kg_year']*z['nitrogen_kWh_kg']+f['source_metal_kg_year']*z['metal_refining_and_fabrication_kWh_kg']
    heat_fraction=1-z['hydrogen_HHV_kWh_kg']/z['electrolysis_liquefaction_kWh_kg_H2']
    f['source_mean_W']=(process_kWh+source_chem_kWh)*1000/base.HOURS
    f['source_heat_W']=(process_kWh+source_chem_kWh*heat_fraction)*1000/base.HOURS
    if q['earth_route'] in ('space_chemical','split_chemical'):
        extra_kWh=earth_feed/z['water_to_hydrogen_mass_ratio']*z['electrolysis_liquefaction_kWh_kg_H2']
        f['depot_mean_W']+=extra_kWh*1000/base.HOURS
        f['depot_heat_W']+=extra_kWh*heat_fraction*1000/base.HOURS
    # New inventories are annualized over construction alongside other initial
    # stocks. Actual commissioning must fill reserves before releasing service.
    # Added tanks and cooling are sized separately; the existing depot already
    # gains tank capacity from the greater incoming source cargo.
    return f


def design_assets(s,d,f,q):
    assets,eng=base.design_assets(s,d,f)
    if q['earth_route']=='earth_chemical':return assets,eng
    split=q['earth_route']=='split_chemical'
    if not split:
        del assets['Earth_transfer_tugs'];del eng['counts']['Earth_transfer_tugs']
    t,p,e=s['transport'],s['power'],s['economics']
    def add(name,mass,cost,life,dest,group):
        assets[name]=dict(mass_kg=mass,manufacture_USD2026=cost,life_years=life,destination=dest,group=group)
    for prefix,key,flow,dest,group in [('Earth','Earth_vehicle','earth_deep_kg_year','leo','earth_transfer'),('Source','Source_vehicle','nea_payload_kg_year','source','source')]:
        if key not in f or (split and prefix=='Earth'):continue
        if prefix=='Source':del assets['asteroid_tugs'];del eng['counts']['asteroid_tugs']
        v=f[key];n=math.ceil(f[flow]*e['capacity_margin']/(v['cargo_kg']*JULIAN_YEAR_DAYS/v['cycle_days']*t['fleet_availability']))
        eng['counts'][prefix+'_route_ships']=n
        for part,a in v['parts'].items():
            add(prefix+'_'+part,a['mass_kg']*n,a['manufacture_USD2026']*n,a['life_years'],dest,group)
            assets[prefix+'_'+part]['vehicle_count']=n
    chemical=q['earth_route'] in ('space_chemical','split_chemical');c=q['chemical'] if chemical else q['electric']
    leo_stock=f['Earth_route_LEO_inventory_kg']
    depot_stock=f['Earth_route_depot_inventory_kg']
    for name,stock,dest in [('LEO',leo_stock,'leo'),('depot',depot_stock,'depot')]:
        if name=='LEO':
            kg=stock*q['depots']['tank_mass_per_stored_kg']
            add(name+'_route_tanks',kg,kg*s['prices']['spacecraft_USD2026_kg'],q['depots']['life_years'],dest,'earth_transfer')
        W=stock*(c['cooling_electric_W_per_stored_kg'] if chemical else c['thermal_control_W_per_stored_kg'])
        array=W/p['end_of_life_factor'];rad=W/(p['emissivity']*STEFAN_BOLTZMANN*p['radiator_K']**4)*p['radiator_kg_m2']
        equipment=W*(c['cooling_hardware_kg_W'] if chemical else c['thermal_control_kg_W'])
        equipment_cost=W*(c['cooling_hardware_USD2026_W'] if chemical else c['thermal_control_USD2026_W'])
        add(name+'_route_storage_controls',equipment+rad+array/p['space_system_W_kg'],equipment_cost+rad*s['prices']['lunar_hardware_USD2026_kg']+array*p['space_USD2026_W'],10,dest,'earth_transfer')
    eng['additional_route_storage_kg']=leo_stock+depot_stock
    eng['LEO_route_storage_kg']=leo_stock
    eng['depot_route_storage_kg']=depot_stock
    return assets,eng


def launch_split(f,t):
    prop=f.get('earth_launch_propellant_kg_year',f['earth_transfer_propellant_kg_year'])
    cargo=f['LEO_kg_year']-prop
    return dict(LEO_kg_year=f['LEO_kg_year'],transfer_propellant_from_Earth_kg_year=prop,
                manufactured_cargo_and_hardware_kg_year=cargo,propellant_fraction=prop/f['LEO_kg_year'],
                launches_at_100t=math.ceil(f['LEO_kg_year']/t['launch_payload_kg']))


def case(s,src,q):
    d=base.settlement_demand(s,src)
    if q['earth_route']=='earth_chemical':
        c=base.case(s,src)
        c['launch_split']={k:launch_split(c[k],s['transport']) for k in ('peak_flows','steady_flows')}
        return c
    assets={}
    for iteration in range(250):
        f=flows(s,d,assets,q)
        new,eng=design_assets(s,d,f,q)
        error=max(abs(a['mass_kg']-assets.get(k,{}).get('mass_kg',0))/max(a['mass_kg'],1) for k,a in new.items())
        assets=new
        if error<1e-10:break
    else:raise ValueError('Industrial replacement and transport iteration did not converge')
    f=flows(s,d,assets,q);steady=flows(s,d,assets,q,False)
    # Owner cash purchases exclude internally supplied propellant. Source mines,
    # fuel plants, storage, fleet electricity and their renewal remain in assets.
    accounting_f=deepcopy(f)
    accounting_f['earth_round_trip']['total']=f['earth_launch_propellant_kg_year']/f['earth_deep_kg_year']/s['transport']['propellant_margin']
    cost,_=base.budget(s,d,assets,accounting_f)
    cost.pop('U_calibration');cost.pop('tariffs_USD2041_kg')
    e,t=s['economics'],s['transport'];I=cost['inflation_factor'];margin=1+e['programme_contingency']
    v=f['Earth_vehicle']
    seed_LEO=0 if q['earth_route']=='split_chemical' else eng['counts']['Earth_route_ships']*q['depots']['seed_earth_fleet_fraction']*v['burns']['ascent_kg']
    seed_cost=seed_LEO*(t['LEO_USD2026_kg']+t['earth_propellant_USD2026_kg'])*I
    cost['additional_route_seed_LEO_kg']=seed_LEO
    cost['additional_route_seed_USD2041']=seed_cost
    cost['capital_direct_USD2041']+=seed_cost
    cost['capital_with_contingency_USD2041']+=seed_cost*margin
    cost['industrial_and_initial_material_campaign_USD2041']+=seed_cost*margin
    return dict(demand=d,assets=assets,peak_flows=f,steady_flows=steady,engineering=eng,cost=cost,
                convergence={'iterations':iteration+1,'relative_mass_residual':error},
                launch_split={k:launch_split(vv,t) for k,vv in [('peak_flows',f),('steady_flows',steady)]})


def scenario(q,name):
    cfg=deepcopy(q);cfg.update(q['architecture_cases'][name]);return cfg


def run(spec=None):
    q,s,src=read_inputs(spec)
    cases={name:case(s,src,scenario(q,name)) for name in q['architecture_cases']}
    comparisons={}
    for name,overrides in q['sensitivities'].items():
        v=scenario(q,'electric_both_legs')
        for path,value in overrides.items():a,b=path.split('.');v[a][b]=value
        try:comparisons[name]={'status':'closed','case':case(s,src,v)}
        except ValueError as exc:comparisons[name]={'status':'infeasible_in_this_sizing_envelope','reason':str(exc)}
    localization={}
    local_s=base.with_case(s,'local_fleet_fabrication')
    for name in ['earth_chemical','split_chemical','space_chemical','electric_earth_leg']:
        localization[name]=case(local_s,src,scenario(q,name))
    # Additional power-specific-mass and efficiency targets are coupled here;
    # they remain a development case, distinct from the central comparison.
    advanced_s=deepcopy(s);advanced_s['power']['space_system_W_kg']=200
    advanced_q=scenario(q,'electric_both_legs');advanced_q['electric']['total_efficiency']=0.5
    comparisons['electric_50pct_200Wkg']={'status':'closed','case':case(advanced_s,src,advanced_q)}
    launch_prices={}
    for price in [50,200,1000]:
        ss=deepcopy(s);ss['transport']['LEO_USD2026_kg']=price
        launch_prices[str(price)]={}
        for name in q['architecture_cases']:
            c=cases[name] if price==200 else case(ss,src,scenario(q,name))
            launch_prices[str(price)][name]={k:c['cost'][k] for k in ('capital_with_contingency_USD2041','industrial_and_initial_material_campaign_USD2041','annual_supply_with_contingency_USD2041')}
    return dict(schema=SCHEMA,date=q['date'],evidence=q['evidence'],inputs=q['inputs'],
                producer={'path':'research/studies/initial_industrial_system/transport.py','sha256':base.sha(__file__),
                  'configuration':{'path':'research/studies/initial_industrial_system/transport_inputs.json','sha256':base.sha(HERE/'transport_inputs.json')},
                  'sources':{'path':'research/studies/initial_industrial_system/transport_sources.json','sha256':base.sha(HERE/'transport_sources.json')}},
                cases=cases,sensitivities=comparisons,source_structure_localization=localization,launch_price_sensitivity=launch_prices,rules=q['rules'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'research/runs/industrial_transport.json');a=p.parse_args()
    result=run();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    for name,c in result['cases'].items():
        b=c['cost'];print(name,'launches',c['launch_split']['peak_flows']['launches_at_100t'],c['launch_split']['steady_flows']['launches_at_100t'],'capital/campaign T$',round(b['capital_with_contingency_USD2041']/1e12,3),round(b['industrial_and_initial_material_campaign_USD2041']/1e12,3),'annual B$',round(b['annual_supply_with_contingency_USD2041']/1e9,3))
    for name,c in result['sensitivities'].items():print(name,c['status'])
if __name__=='__main__':main()
