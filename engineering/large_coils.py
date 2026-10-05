"""Screening hardware and renewal of large closed current paths.

Engineering parameters are scenarios, not demonstrated complete systems.
Thin circular loop self-inductance and its energy derivative give the
magnetic tensile lower estimate. Terrain, gravity, bending and faults remain.
"""
import numpy as np
from shared import constants as K
from .electromagnetic import radiator, collector


def winding(radius, current, field_limit=10., bundle_radius=None, heat_leak=.05,
            cop=.005, current_density=1e8, specific_strength=1e6,
            density=6000., terminal_current=50000., joint_spacing=10000., joint_ohm=1e-9):
    current = abs(current); perimeter = 2*np.pi*radius
    packed = np.sqrt(current/(np.pi*current_density*.5))
    field_radius = K.VACUUM_PERMEABILITY*current/(2*np.pi*field_limit)
    bundle = max(.1, packed, field_radius) if bundle_radius is None else bundle_radius
    if not 0 < bundle < radius/10:
        return dict(admissible=False, reason='Bundle is outside the thin-loop model',
                    radius_m=radius, ampere_turns=current, bundle_radius_m=bundle)
    logarithm = np.log(8*radius/bundle)
    inductance = K.VACUUM_PERMEABILITY*radius*(logarithm-2)
    energy = .5*inductance*current**2
    tension = K.VACUUM_PERMEABILITY*current**2/(4*np.pi)*(logarithm-1)
    conductor = perimeter*current/current_density*density
    support = perimeter*tension/specific_strength
    # 10 kg/m2 cold envelope is an additional assumed warm-shell allowance.
    envelope = perimeter*2*np.pi*(bundle+.25)
    # Total turn length / splice interval, with terminal current through each splice.
    turn_length = perimeter*current/terminal_current
    joints = turn_length/joint_spacing*terminal_current**2*joint_ohm
    cold = envelope*heat_leak+joints
    electrical = cold/cop; rad = radiator(electrical+cold)
    cryostat = envelope*10.
    refrigerator = electrical/300.
    total = conductor+support+cryostat+refrigerator+rad['mass_kg']
    generation = collector(electrical)
    return dict(admissible=True, radius_m=radius, ampere_turns=current,
        moment_A_m2=current*np.pi*radius**2, bundle_radius_m=bundle,
        conductor_pack_radius_m=packed, bundle_self_field_T=K.VACUUM_PERMEABILITY*current/(2*np.pi*bundle),
        engineering_current_density_A_m2=current_density, specific_strength_J_kg=specific_strength,
        terminal_current_A=terminal_current, aggregate_turn_length_m=turn_length,
        inductance_H=inductance, self_energy_J=energy, hoop_tension_N=tension,
        self_line_load_N_m=tension/radius, conductor_mass_kg=conductor,
        ideal_tensile_mass_kg=support, cryostat_mass_kg=cryostat,
        cold_envelope_m2=envelope, heat_leak_W_m2=heat_leak, cop=cop,
        joints_cold_W=joints, cold_W=cold, refrigerator_W=electrical,
        refrigerator_mass_kg=refrigerator, radiator=rad, generation=generation,
        screened_mass_without_generation_kg=total,
        screened_mass_with_generation_kg=total+generation['installed_mass_kg'],
        uncertainty='No foundations, terrain, atmospheric convection, gravity/bending, deployment, fault dump or night storage. The shell and refrigerator mass are engineering scenarios; joints extrapolate metre-scale measurements.')


def renewal(mass, energy, life_years=100., recovery=.999, manufacture_J_kg=1e8):
    seconds = life_years*K.JULIAN_YEAR_DAYS*K.JULIAN_DAY
    return dict(life_years=life_years, gross_kg_s=mass/seconds,
        fresh_kg_s=mass/seconds*(1-recovery), manufacture_W=mass/seconds*manufacture_J_kg,
        full_recharge_W=energy/seconds, unrecovered_field_energy_per_replacement_J=energy,
        scope='Full field loss at replacement is a sensitivity; recovery and useful life are unqualified.')


def upstream_width(distance, target_radius, speed=400e3, lateral_speed=60e3, angle_deg=5.):
    """Kinematic erosion/offset sensitivity, not a magnetotail solution."""
    flight = distance/speed
    erosion = lateral_speed*flight
    offset = distance*np.tan(np.deg2rad(angle_deg))
    return dict(distance_m=distance, target_radius_m=target_radius, wind_speed_m_s=speed,
        transverse_speed_m_s=lateral_speed, direction_error_deg=angle_deg,
        transit_s=flight, lateral_erosion_m=erosion, centre_offset_m=offset,
        required_obstacle_radius_m=target_radius+erosion+offset,
        evidence='Parametric kinematic margin; magnetic guidance can change refilling in either direction.')
