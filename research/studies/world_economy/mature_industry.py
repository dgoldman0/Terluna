"""Mature 1.2-atmosphere industry: physical ratios, costs and recurring orders."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from shared.constants import MOON_SURFACE_GRAVITY, MOON_RADIUS, STANDARD_GRAVITY
from shared.provenance import ROOT, constants_used
from .industrial_physics import ideal_pump_work
from .mature_economics import component_cost, compare_locations, renewal_demand

HERE=Path(__file__).resolve().parent
INPUT_PATH=HERE/'mature_industry_inputs.json'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(path=INPUT_PATH, root=ROOT):
    s=json.loads(Path(path).read_text())
    if s['schema']!='terluna.research.mature-industry-inputs/1':
        raise ValueError('Unexpected mature industry inputs')
    products={}
    for name,spec in s['products'].items():
        p=Path(root)/name
        if digest(p)!=spec['sha256']:raise ValueError(f'Input changed; review and rebind: {name}')
        d=json.loads(p.read_text())
        if d['schema']!=spec['schema']:raise ValueError(f'Unexpected product schema: {name}')
        products[name]=d
    return s,products


def environment_ratios(earth_pressure, moon_pressure, target, volume, speed, process_load):
    if not 0 < target < min(earth_pressure,moon_pressure) or min(volume,speed)<=0 or process_load<0:
        raise ValueError('Positive pressures, volume and speed; target below surroundings')
    g=MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    density=moon_pressure/earth_pressure
    work_e=ideal_pump_work(volume,earth_pressure,target)['ideal_reversible_work_J']
    work_m=ideal_pump_work(volume,moon_pressure,target)['ideal_reversible_work_J']
    time_e=volume/speed*math.log(earth_pressure/target)
    time_m=volume/speed*math.log(moon_pressure/target)
    pressure=(moon_pressure-target)/(earth_pressure-target)
    return dict(earth_pressure_Pa=earth_pressure,moon_pressure_Pa=moon_pressure,
        gravity_ratio=g,pressure_ratio=moon_pressure/earth_pressure,
        pressure_load_ratio=pressure,ideal_buckling_thickness_ratio=math.sqrt(pressure),
        initial_gas_inventory_ratio=density,ideal_pump_work_ratio=work_m/work_e,
        earth_ideal_pump_work_J=work_e,moon_ideal_pump_work_J=work_m,
        constant_speed_pumpdown_ratio=time_m/time_e,
        added_ideal_pumpdown_seconds=time_m-time_e,
        internal_process_load_minimum_inlet_speed_m3_s=process_load/target,
        added_ideal_fixed_throughput_compression_power_W=process_load*math.log(density),
        liquid_hydrostatic_stress_ratio=g,
        same_stress_self_weight_height_ratio=1/g,
        capillary_length_and_same_Bond_size_ratio=g**-.5,
        same_liquid_Rayleigh_number_length_ratio=g**(-1/3),
        same_length_liquid_Bond_and_Rayleigh_ratio=g,
        gas_Rayleigh_ratio_same_T_properties=g*density*density,
        self_weight_beam_same_absolute_sag_span_ratio=g**-.25,
        self_weight_beam_same_stress_span_ratio=g**-.5,
        same_mass_wing_area_lift_coefficient_speed_ratio=math.sqrt(g/density),
        ideal_hover_induced_power_same_mass_disk_area_ratio=g**1.5/math.sqrt(density),
        buoyant_supported_mass_per_volume_ratio=density,
        wind_pressure_same_speed_ratio=density,
        interpretation='Initial pumpdown is a no-added-gas diagnostic, separate from the fixed internal-load speed requirement. Same material/process geometry unless a changed length is explicitly identified. Gas comparisons hold composition, temperature, dynamic viscosity and conductivity fixed. Hover is ideal induced power only; buoyant mass per volume does not increase from reduced gravity itself. No factory-cost or certified size multiplier.')


def centrifuge_screen(depth, gradient_fraction, speed, coriolis_fraction):
    if min(depth,gradient_fraction,coriolis_fraction)<=0 or speed<0:
        raise ValueError('Positive dimensions/tolerances and nonnegative process speed')
    radius=depth/gradient_fraction
    omega=math.sqrt(MOON_SURFACE_GRAVITY/radius)
    coriolis=2*omega*speed/MOON_SURFACE_GRAVITY
    return dict(radial_work_depth_m=depth, allowed_radial_g_fraction=gradient_fraction,
        radius_for_radial_gradient_m=radius, rotation_rpm=omega*60/(2*math.pi),
        transverse_speed_m_s=speed, maximum_Coriolis_to_gravity_ratio=coriolis,
        radius_for_selected_maximum_Coriolis_ratio_m=4*speed**2/(coriolis_fraction**2*MOON_SURFACE_GRAVITY),
        lunar_central_field_radial_change_over_same_height_fraction=1-(MOON_RADIUS/(MOON_RADIUS+depth))**2)


def cost_cases(s, physical):
    c=s['cost_screen'];rows=[]
    for recipe,shares in c['recipes'].items():
        for saving in c['gravity_step_savings']:
            for site in c['array_site_U_per_qualified_unit']:
                for destination in c['delivery_U_per_qualified_unit']['moon']:
                    locations={}
                    for world in ('earth','moon','array'):
                        factors={'gravity_step':1-(c['earth_step_saving'] if world=='earth' else saving)}
                        if world=='moon':factors.update(pressure_gas=physical['ideal_pump_work_ratio'],
                            pressure_shell=physical['ideal_buckling_thickness_ratio'],pump_cycle=physical['constant_speed_pumpdown_ratio'])
                        if world=='array':factors.update(pressure_gas=0,pressure_shell=0,pump_cycle=0)
                        locations[world]=component_cost(c['base_attempt_cost_U'],shares,factors,c['accepted_yield'],
                            c['annual_qualified_output'],site*c['annual_qualified_output'] if world=='array' else 0,
                            c['delivery_U_per_qualified_unit'][world][destination])
                    rows.append(dict(recipe=recipe,assumed_gravity_step_saving=saving,array_site_U_per_good=site,
                        destination=destination,locations=locations,comparison=compare_locations(locations,'gravity_step')))
    yield_cases=[]
    for row in rows:
        if (row['assumed_gravity_step_saving']==c['reference_step_saving'] and
            row['array_site_U_per_good']==c['reference_array_site_U_per_qualified_unit']):
            for accepted in c['yield_sensitivity']:
                locations=dict(row['locations']);old=locations['moon']
                locations['moon']=component_cost(old['base_attempt_cost_U'],old['shares'],old['factors'],accepted,
                    c['annual_qualified_output'],delivery_cost_per_good=old['delivery_cost_per_qualified_unit_U'])
                yield_cases.append(dict(recipe=row['recipe'],destination=row['destination'],moon_accepted_yield=accepted,
                    comparison=compare_locations(locations,'gravity_step')))
    return rows,yield_cases


def results():
    s,p=load_inputs();old=p['research/studies/world_economy/results/industrial_physics.json']
    pressure=sorted({r['against_vacuum']['internal_pressure_Pa'] for r in old['enclosure_cases']})
    e=s['environment'];phys=environment_ratios(pressure[0],pressure[-1],e['target_vacuum_Pa'],e['chamber_m3'],e['effective_speed_m3_s'],e['fixed_process_gas_load_Pa_m3_s'])
    c=s['centrifuges'];cent=[centrifuge_screen(d,f,v,c['coriolis_tolerance_fraction']) for d in c['radial_work_depths_m'] for f in c['allowed_radial_g_fraction'] for v in c['transverse_process_speeds_m_s']]
    costs,yields=cost_cases(s,phys)
    construction=p['research/studies/world_economy/results/construction.json']
    con=next(r for r in construction['cases'] if r['cost_case']=='advanced_industrial' and r['water_case']=='surface_only')
    aerial=next(r for r in con['installed_assets'] if r['element']=='aerial_carriers_and_control')
    comparison=p['research/studies/world_economy/results/comparison.json']
    cfc=next(r for r in comparison['mature_accounts']['close_parity']['incremental_cfc_by_asset'] if r['element']=='aerial_carriers_and_control')
    r=s['renewal'];array=con['included_array_structure_comparison']
    demand={'aerial':renewal_demand(aerial['quantity'],cfc['life_years'],r['aerial_eligible_fraction'],r['aerial_award_share'],r['aerial_capacity_kg_year'],r['aerial_budget_U_year'],aerial['unit_cost_U']),
        'array_structure':[renewal_demand(array['mass_kg'],life,r['array_eligible_fraction'],r['array_award_share'],r['array_capacity_kg_year'],r['array_budget_U_year'],array['separately_valued_U']/array['mass_kg']) for life in r['array_structure_life_years']],
        'shield_gross_renewal_Gt_year':p['resources/results/industry.json']['plant_traffic']['upkeep_gt_yr']}
    names=('mature_industry.py','mature_industry_inputs.json','mature_economics.py','industrial_physics.py')
    return dict(schema='terluna.research.mature-industry/1',producer=dict(files={str((HERE/n).relative_to(ROOT)):digest(HERE/n) for n in names},constants=constants_used([__file__]),input_products=s['products']),
        evidence=s['evidence'],environment=phys,centrifuges=cent,cost_cases=costs,yield_sensitivities=yields,recurring_demand=demand,
        assumptions={k:s[k] for k in ('environment','centrifuges','cost_screen','renewal')},boundaries=s['boundaries'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'research/runs/world_economy/mature_industry.json')
    a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n');print(a.output)


if __name__=='__main__':main()
