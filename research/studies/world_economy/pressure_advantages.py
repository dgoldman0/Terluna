"""Useful pressure effects at 1.2 atm, separated from lunar gravity.

These are matched-condition engineering screens, not product qualification or
new economic output. The atmosphere domain supplies a reusable water-property
law; existing domain results are neither imported nor regenerated.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from atmosphere.radiative_convective import thermodynamics as water
from shared.constants import GAS_CONSTANT, MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY
from shared.provenance import ROOT, constants_used
from .industrial_physics import ideal_pump_work

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE / 'pressure_advantages_inputs.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def positive(*values):
    if any(isinstance(v, bool) or not math.isfinite(v) or v <= 0 for v in values):
        raise ValueError('Finite positive quantities required')


def boiling_temperature(pressure):
    """Invert the existing liquid-water saturation law within 273.16--400 K."""
    positive(pressure)
    lo, hi = water.TRIPLE_POINT_K, 400.
    if not float(water.saturation_pressure(lo)) <= pressure <= float(water.saturation_pressure(hi)):
        raise ValueError('Pressure outside this boiling-temperature bracket')
    for _ in range(70):
        mid = (lo + hi) / 2
        if float(water.saturation_pressure(mid)) < pressure:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def pressure_scalings(ratio):
    positive(ratio)
    return dict(density_and_gross_buoyant_mass_per_volume=ratio,
                volume_for_same_buoyant_mass=1/ratio,
                similar_envelope_area_for_same_buoyant_mass=ratio**(-2/3),
                wing_speed_at_same_lift=ratio**-.5,
                ideal_induced_hover_power_at_same_weight_and_disk=ratio**-.5,
                wind_power_and_dynamic_load_at_same_speed=ratio,
                natural_convection_Rayleigh=ratio**2,
                heat_transfer_if_Nusselt_scales_as_Rayleigh_quarter=ratio**.5,
                heat_transfer_if_Nusselt_scales_as_Rayleigh_third=ratio**(2/3))


def air_circuit(density, mass_flow, area, loss_coefficient, fixed_loss, efficiency):
    """Fixed-geometry low-Mach air circuit; K and fan efficiency are prescribed."""
    positive(density, mass_flow, area, efficiency)
    if efficiency > 1 or min(loss_coefficient, fixed_loss) < 0 or not all(
            math.isfinite(v) for v in (loss_coefficient, fixed_loss)):
        raise ValueError('Finite nonnegative losses and efficiency at most one')
    volume = mass_flow/density
    speed = volume/area
    dynamic_loss = loss_coefficient*density*speed**2/2
    return dict(density_kg_m3=density, mass_flow_kg_s=mass_flow,
                volume_flow_m3_s=volume, speed_m_s=speed,
                dynamic_pressure_loss_Pa=dynamic_loss, fixed_loss_Pa=fixed_loss,
                total_pressure_loss_Pa=dynamic_loss+fixed_loss,
                electric_power_W=volume*(dynamic_loss+fixed_loss)/efficiency)


def clamp_case(ambient, pressure_drop, volume, speed):
    positive(ambient, pressure_drop, volume, speed)
    if pressure_drop >= ambient:
        raise ValueError('Finite pumpdown requires positive final absolute pressure')
    final = ambient-pressure_drop
    return dict(ambient_Pa=ambient, clamping_pressure_Pa=pressure_drop,
                final_absolute_Pa=final,
                constant_speed_no_load_pumpdown_s=volume/speed*math.log(ambient/final),
                ideal_work_J=ideal_pump_work(volume, ambient, final)['ideal_reversible_work_J'])


def compressor(inlet, outlet, temperature, gamma):
    positive(inlet, outlet, temperature, gamma)
    if outlet < inlet or gamma <= 1:
        raise ValueError('Outlet at least inlet and gamma greater than one')
    ratio = outlet/inlet
    return dict(inlet_Pa=inlet, outlet_Pa=outlet, pressure_ratio=ratio,
                reversible_isothermal_J_mol=GAS_CONSTANT*temperature*math.log(ratio),
                ideal_isentropic_J_mol=GAS_CONSTANT*gamma/(gamma-1)*temperature*(ratio**((gamma-1)/gamma)-1))


def liquid_margin(ambient, temperature, density, gravity):
    positive(ambient, temperature, density, gravity)
    vapour = float(water.saturation_pressure(temperature))
    margin = ambient-vapour
    return dict(ambient_Pa=ambient, temperature_K=temperature, density_kg_m3=density,
                gravity_m_s2=gravity, vapour_pressure_Pa=vapour,
                pressure_above_vapour_Pa=margin, unboiled_surface_possible=margin>=0,
                ideal_suction_column_limit_m=margin/(density*gravity) if margin>=0 else None)


def pressure_level_height(surface_pressure, target_pressure, temperature, molar_mass):
    """Spherical, isothermal ideal-gas column; not a climate sounding."""
    positive(surface_pressure, target_pressure, temperature, molar_mass)
    if target_pressure > surface_pressure:
        raise ValueError('Target pressure above surface pressure')
    scale = GAS_CONSTANT*temperature/(molar_mass*MOON_SURFACE_GRAVITY)
    x = scale*math.log(surface_pressure/target_pressure)
    if x >= MOON_RADIUS:
        raise ValueError('Target outside the finite spherical-column potential')
    return MOON_RADIUS*x/(MOON_RADIUS-x)


def load_inputs(path=INPUT_PATH, root=ROOT):
    s=json.loads(Path(path).read_text())
    if s['schema']!='terluna.research.pressure-advantages-inputs/1':
        raise ValueError('Unexpected pressure advantages input schema')
    products={}
    for name,spec in s['products'].items():
        p=Path(root)/name
        if digest(p)!=spec['sha256']:
            raise ValueError(f'Input changed; review and rebind: {name}')
        data=json.loads(p.read_text())
        if data['schema']!=spec['schema']:
            raise ValueError(f'Unexpected input schema: {name}')
        products[name]=data
    return s,products


def results():
    s,products=load_inputs()
    env=products['research/studies/world_economy/results/mature_industry.json']['environment']
    p0,p1=env['earth_pressure_Pa'],env['moon_pressure_Pa']
    e=s['environment']; t=e['temperature_K']; ratio=p1/p0
    y=e['moon_oxygen_mole_fraction']; ye=e['reference_oxygen_mole_fraction']
    mass=(1-y)*water.MOLAR_MASS['N2']+y*water.MOLAR_MASS['O2']
    mass_e=(1-ye)*water.MOLAR_MASS['N2']+ye*water.MOLAR_MASS['O2']
    cp=((1-y)*water.MOLAR_CP['N2']+y*water.MOLAR_CP['O2'])/mass
    rho0=p0*mass/(GAS_CONSTANT*t); rho1=p1*mass/(GAS_CONSTANT*t)
    c=s['cooling']; flow=c['heat_transport_W']/(cp*c['air_temperature_rise_K'])
    fans=[]
    for fixed in c['fixed_pressure_losses_Pa']:
        args=(c['duct_area_m2'],c['loss_coefficient'],fixed,c['fan_efficiency'])
        base=air_circuit(rho0,flow,*args)
        same_mass=air_circuit(rho1,flow,*args)
        same_volume=air_circuit(rho1,flow*ratio,*args)
        fans.append(dict(fixed_loss_Pa=fixed,one_atm=base,locked_same_heat_transport=same_mass,
                         locked_same_volume_flow=same_volume,
                         same_heat_transport_power_ratio=same_mass['electric_power_W']/base['electric_power_W'],
                         same_volume_flow_power_ratio=same_volume['electric_power_W']/base['electric_power_W']))
    v=s['vacuum_assisted']; drop=v['required_clamping_pressure_Pa']
    clamps=[clamp_case(p,drop,v['rigid_volume_m3'],v['constant_pump_speed_m3_s']) for p in (p0,p1)]
    compression=[]
    for outlet in s['compression']['outlet_absolute_pressures_Pa']:
        rows=[compressor(p,outlet,t,s['compression']['gamma']) for p in (p0,p1)]
        compression.append(dict(outlet_Pa=outlet,one_atm=rows[0],locked=rows[1],
            isothermal_work_saving_fraction=1-rows[1]['reversible_isothermal_J_mol']/rows[0]['reversible_isothermal_J_mol'],
            isentropic_work_saving_fraction=1-rows[1]['ideal_isentropic_J_mol']/rows[0]['ideal_isentropic_J_mol']))
    sat=float(water.saturation_pressure(t))
    names=('pressure_advantages.py','pressure_advantages_inputs.json','pressure_advantages_sources.json','industrial_physics.py')
    return dict(schema='terluna.research.pressure-advantages/1',
        producer=dict(files={str((HERE/n).relative_to(ROOT)):digest(HERE/n) for n in names},
                      thermodynamic_law=dict(path=str(Path(water.__file__).relative_to(ROOT)),sha256=digest(water.__file__)),
                      constants=constants_used([__file__]),input_products=s['products']),
        evidence=s['evidence'],reading_rule=s['reading_rule'],assumptions=s,
        environment=dict(reference_Pa=p0,locked_Pa=p1,temperature_K=t,pressure_ratio=ratio,
                         same_composition_reference_density_kg_m3=rho0,locked_density_kg_m3=rho1,
                         binary_composition_corrected_density_ratio=ratio*mass/mass_e,
                         molar_mass_kg_mol=mass,mass_cp_J_kg_K=cp),
        pressure_only_scalings=pressure_scalings(ratio),
        ideal_isothermal_height_of_one_atm_m=pressure_level_height(p1,p0,t,mass),
        gas_partial_pressures=dict(reference_oxygen_Pa=p0*ye,locked_oxygen_Pa=p1*y,
            dry_oxygen_ratio=p1*y/(p0*ye),
            saturated_water_interface_oxygen_ratio=(p1-sat)*y/((p0-sat)*ye),
            reference_trace_CO2_Pa=p0*e['trace_CO2_mole_fraction'],locked_trace_CO2_Pa=p1*e['trace_CO2_mole_fraction']),
        cooling= fans,
        vacuum_assisted=dict(same_absolute_bag_pressure_Pa=v['bag_absolute_pressure_Pa'],
            one_atm_clamping_Pa=p0-v['bag_absolute_pressure_Pa'],locked_clamping_Pa=p1-v['bag_absolute_pressure_Pa'],
            added_clamping_force_N_per_m2=p1-p0,fixed_force_cases=clamps,
            fixed_force_pumpdown_ratio=clamps[1]['constant_speed_no_load_pumpdown_s']/clamps[0]['constant_speed_no_load_pumpdown_s'],
            fixed_force_ideal_work_ratio=clamps[1]['ideal_work_J']/clamps[0]['ideal_work_J']),
        liquid=dict(boiling_temperature_K=[boiling_temperature(p) for p in (p0,p1)],
            added_pressure_margin_Pa=p1-p0,
            added_water_column_m=(p1-p0)/(s['liquid']['density_kg_m3']*MOON_SURFACE_GRAVITY),
            cases=[liquid_margin(p,temp,s['liquid']['density_kg_m3'],g)
                   for p,g in ((p0,STANDARD_GRAVITY),(p0,MOON_SURFACE_GRAVITY),(p1,MOON_SURFACE_GRAVITY))
                   for temp in s['liquid']['temperatures_K']]),
        compression=compression,
        conditional_whole_cost_savings=[dict(fan_electricity_cost_share=share,
            whole_cost_saving_fraction=share*(1-fans[0]['same_heat_transport_power_ratio']))
            for share in c['fan_electricity_cost_shares']])


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,default=ROOT/'research/runs/world_economy/pressure_advantages.json')
    args=p.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__=='__main__':
    main()
