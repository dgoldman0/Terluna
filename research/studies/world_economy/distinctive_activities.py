"""Matched-condition references for activities in the Open Moon environment.

Natural gravity, rotating analogues, flight and liquid geometry are compared.
These calculations establish neither human tolerances nor market advantage.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from shared.constants import (JULIAN_DAY, MOON_RADIUS, MOON_SURFACE_GRAVITY,
                              SIDEREAL_MONTH_DAYS, STANDARD_GRAVITY)
from shared.provenance import ROOT, constants_used

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE / 'distinctive_activities_inputs.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def positive(*values):
    if any(isinstance(v, bool) or not math.isfinite(v) or v <= 0 for v in values):
        raise ValueError('Finite positive quantities required')


def nonnegative(*values):
    if any(isinstance(v, bool) or not math.isfinite(v) or v < 0 for v in values):
        raise ValueError('Finite nonnegative quantities required')


def flight_scalings(gravity_ratio, density_ratio):
    """Same mass, wing/disk area and aerodynamic coefficients; still air."""
    positive(gravity_ratio, density_ratio)
    return dict(speed_at_fixed_lift_coefficient=math.sqrt(gravity_ratio/density_ratio),
                aerodynamic_power_at_fixed_lift_and_drag_coefficients=gravity_ratio**1.5/math.sqrt(density_ratio),
                ideal_induced_hover_power=gravity_ratio**1.5/math.sqrt(density_ratio),
                work_per_distance_at_fixed_lift_to_drag=gravity_ratio,
                parasite_drag_at_fixed_speed_and_coefficient=density_ratio)


def ballistic(gravity, vertical_speed):
    """Prescribed takeoff speed; no aerodynamic drag or biomechanics."""
    positive(gravity)
    nonnegative(vertical_speed)
    return dict(height_m=vertical_speed**2/(2*gravity),
                return_time_s=2*vertical_speed/gravity)


def rotating_venue(gravity, radius, height, arc_width, transverse_speed):
    """Uniform rotation; height points inward, width is circumferential arc."""
    positive(gravity, radius)
    nonnegative(height, arc_width, transverse_speed)
    if height >= radius or arc_width > 2*math.pi*radius:
        raise ValueError('Venue must lie inside the rim and span at most one turn')
    omega = math.sqrt(gravity/radius)
    return dict(radius_m=radius, angular_speed_rad_s=omega, rpm=omega*60/(2*math.pi),
                fractional_gravity_change=height/radius,
                down_direction_change_degrees=math.degrees(arc_width/radius),
                coriolis_to_rim_gravity=2*omega*transverse_speed/gravity)


def radius_requirements(gravity, height, arc_width, transverse_speed,
                        gravity_tolerance, direction_tolerance_degrees, coriolis_tolerance):
    """Chosen fidelity constraints, not comfort or habitat feasibility limits."""
    positive(gravity, height, arc_width, transverse_speed, gravity_tolerance,
             direction_tolerance_degrees, coriolis_tolerance)
    if gravity_tolerance >= 1 or direction_tolerance_degrees >= 180:
        raise ValueError('Use a local venue with fractional gravity tolerance below one')
    requirements = dict(gravity_gradient_m=height/gravity_tolerance,
                        direction_change_m=arc_width/math.radians(direction_tolerance_degrees),
                        coriolis_m=4*transverse_speed**2/(gravity*coriolis_tolerance**2))
    radius = max(requirements.values())
    return dict(constraint_radii=requirements, combined_radius_m=radius,
                venue=rotating_venue(gravity,radius,height,arc_width,transverse_speed))


def hydrostatic_depth(surface_pressure, absolute_pressure, density, gravity):
    """Constant-density reference; an apparatus pressure budget, not diving advice."""
    positive(surface_pressure, absolute_pressure, density, gravity)
    if absolute_pressure < surface_pressure:
        raise ValueError('Target absolute pressure below the water surface pressure')
    return (absolute_pressure-surface_pressure)/(density*gravity)


def gravity_wave(gravity, wavelength, depth):
    """Linear, inviscid gravity wave over a flat bottom; no surface tension."""
    positive(gravity, wavelength, depth)
    k = 2*math.pi/wavelength
    omega = math.sqrt(gravity*k*math.tanh(k*depth))
    return dict(angular_frequency_rad_s=omega, phase_speed_m_s=omega/k,
                period_s=2*math.pi/omega)


def load_inputs(path=INPUT_PATH, root=ROOT):
    spec=json.loads(Path(path).read_text())
    if spec['schema'] != 'terluna.research.distinctive-activities-inputs/1':
        raise ValueError('Unexpected distinctive activities input schema')
    products={}
    for name, identity in spec['products'].items():
        path=Path(root)/name
        if digest(path) != identity['sha256']:
            raise ValueError(f'Input changed; review and rebind: {name}')
        data=json.loads(path.read_text())
        if data['schema'] != identity['schema']:
            raise ValueError(f'Unexpected input schema: {name}')
        products[name]=data
    return spec, products


def results():
    spec,products=load_inputs()
    env=next(iter(products.values()))['environment']
    ratio=MOON_SURFACE_GRAVITY/STANDARD_GRAVITY
    venue=spec['venue']
    requirements=radius_requirements(MOON_SURFACE_GRAVITY,
        venue['height_m'], venue['circumferential_arc_width_m'], venue['transverse_speed_m_s'],
        venue['gravity_tolerance_fraction'],venue['direction_tolerance_degrees'],
        venue['coriolis_tolerance_fraction'])
    omega_moon=2*math.pi/(SIDEREAL_MONTH_DAYS*JULIAN_DAY)
    water=spec['water']; waves=spec['waves']
    depths={name:hydrostatic_depth(pressure,water['target_absolute_pressure_Pa'],
        water['density_kg_m3'],gravity) for name,pressure,gravity in
        [('earth',env['reference_Pa'],STANDARD_GRAVITY),('moon',env['locked_Pa'],MOON_SURFACE_GRAVITY)]}
    files=['research/studies/world_economy/distinctive_activities.py',
           'research/studies/world_economy/distinctive_activities_inputs.json',
           'research/studies/world_economy/distinctive_activities_sources.json']
    return dict(schema='terluna.research.distinctive-activities/1',evidence=spec['evidence'],
        reading_rule=spec['reading_rule'],assumptions=spec,
        producer=dict(files={f:digest(ROOT/f) for f in files},constants=constants_used(files)),
        environment=dict(earth_gravity_m_s2=STANDARD_GRAVITY,moon_gravity_m_s2=MOON_SURFACE_GRAVITY,
            gravity_ratio=ratio,matched_temperature_K=env['temperature_K'],
            earth_pressure_Pa=env['reference_Pa'],moon_pressure_Pa=env['locked_Pa'],
            binary_density_ratio=env['binary_composition_corrected_density_ratio']),
        flight=flight_scalings(ratio,env['binary_composition_corrected_density_ratio']),
        ballistic={name:ballistic(gravity,spec['vertical_takeoff_speed_m_s']) for name,gravity in
                   [('earth',STANDARD_GRAVITY),('moon',MOON_SURFACE_GRAVITY)]},
        rotating_analogue=requirements,
        natural_moon_venue=dict(
            fractional_gravity_change=1-(MOON_RADIUS/(MOON_RADIUS+venue['height_m']))**2,
            down_direction_change_degrees=math.degrees(venue['circumferential_arc_width_m']/MOON_RADIUS),
            maximum_coriolis_to_gravity=2*omega_moon*venue['transverse_speed_m_s']/MOON_SURFACE_GRAVITY),
        water=dict(depth_m=depths,depth_ratio=depths['moon']/depths['earth'],
                   same_depth_gauge_pressure_ratio=ratio,
                   same_liquid_capillary_length_ratio=ratio**-.5,
                   gravity_waves={name:gravity_wave(gravity,waves['wavelength_m'],waves['depth_m'])
                     for name,gravity in [('earth',STANDARD_GRAVITY),('moon',MOON_SURFACE_GRAVITY)]}))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'research/runs/world_economy/distinctive_activities.json')
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
