"""Preliminary square sail hardware accounting, SI units throughout.

This is a conditional sizing screen, not a qualified spacecraft design. Every
material strength, modulus and subsystem specific mass below is an explicitly
assumed admissible design target, not a measured property or flight heritage.
The square membrane has fixed, straight edges and uniform biaxial tension. A
mass-accounted perimeter space truss supplies those edge reactions. Its local
member screens do not prove global stiffness, deployment or fatigue life.

Momentum capacity is an input envelope. This module does not turn an inherited
rotating-frame torque integral into a validated wheel demand. Use the assembled
locked inertia with a conservative attitude dynamics calculation separately.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from math import ceil, isfinite, pi, sqrt


def _positive(**values: float) -> None:
    if any(not isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f"Expected finite positive values: {values}")


def _nonnegative(**values: float) -> None:
    if any(not isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f"Expected finite nonnegative values: {values}")


@dataclass(frozen=True)
class Material:
    """Assumed density, modulus and admissible stress, including desired margins."""
    density_kg_m3: float = 1600.0
    youngs_modulus_Pa: float = 100e9
    admissible_stress_Pa: float = 500e6

    def __post_init__(self):
        _positive(**asdict(self))


def membrane_coefficient(odd_terms: int = 31) -> float:
    """Centre sag coefficient for -T laplacian(w)=q on a fixed square.

    w(centre)=C q L^2/T. Double odd sine series, with `odd_terms` terms
    in each direction. C converges to approximately 0.07367.
    """
    if not isinstance(odd_terms, int) or odd_terms < 1:
        raise ValueError("odd_terms must be a positive integer")
    return 16 / pi**4 * sum(
        (-1)**(i + j) / ((2*i+1)*(2*j+1)*((2*i+1)**2+(2*j+1)**2))
        for i in range(odd_terms) for j in range(odd_terms))


def square_membrane(side_m: float, pressure_Pa: float, tension_N_m: float,
                    areal_mass_kg_m2: float, *, odd_terms: int = 31) -> dict:
    """Small-deflection static sag and first frequency of an ideal membrane.

    f11=sqrt(T/sigma)/(sqrt(2)*L). Static slope is evaluated along the
    middle of one edge; its convergence is slower than centre displacement.
    Uniform |pressure| also bounds static deflection for spatially signed
    pressure whose magnitude is no larger, by the Poisson comparison principle.
    """
    _positive(side_m=side_m, tension_N_m=tension_N_m,
              areal_mass_kg_m2=areal_mass_kg_m2)
    _nonnegative(pressure_Pa=pressure_Pa)
    coefficient = membrane_coefficient(odd_terms)
    edge_slope_coefficient = 16 / pi**3 * sum(
        (-1)**j / ((2*j+1)*((2*i+1)**2+(2*j+1)**2))
        for i in range(odd_terms) for j in range(odd_terms))
    frequency = sqrt(tension_N_m / areal_mass_kg_m2) / (sqrt(2)*side_m)
    return dict(centre_sag_m=coefficient*pressure_Pa*side_m**2/tension_N_m,
                edge_midpoint_slope_rad=edge_slope_coefficient*pressure_Pa*side_m/tension_N_m,
                fundamental_frequency_Hz=frequency, fundamental_period_s=1/frequency,
                edge_in_plane_force_N=tension_N_m*side_m,
                centre_coefficient=coefficient, odd_terms=odd_terms)


def _tube_area(force_N: float, length_m: float, radius_m: float,
               material: Material) -> tuple[float, float]:
    """Thin circular tube; pinned Euler and admissible stress, no local buckling."""
    stress_limit = min(material.admissible_stress_Pa,
                       pi**2*material.youngs_modulus_Pa*radius_m**2/(2*length_m**2))
    return force_N/stress_limit, stress_limit


def perimeter_truss(side_m: float, tension_N_m: float, depth_m: float = 100.,
                    bay_m: float = 50., chord_radius_m: float = 2.,
                    material: Material = Material(), joint_mass_fraction: float = .2,
                    allowable_edge_bow_m: float = 5.) -> dict:
    """Four square-section space trusses close the membrane's in-plane loads.

    Each side is a simply supported beam under line load T, with Mmax=T L^2/8
    and compression from the adjacent side C=T L/2. Four chords share C; the
    two compression chords each take M/(2d). All chords use that maximum size.
    Every bay has tension-only X braces on all four faces and four transverse
    posts. Two faces carry the in-plane shear; the other faces have the same
    bracing as a declared reserve. Post force is conservatively T*bay/2.

    Half the allowed edge bow is assigned to Euler-Bernoulli bending, half to
    ideal pin-jointed web shear. For four corner chords I=A_chord*d^2;
    bending bow=5*T*L^4/(384*E*I). The two active diagonal faces give effective
    shear rigidity 2*E*A_brace*bay*d^2/diagonal^3 and shear bow=T*L^2/(8*GA).
    All four faces have the same brace size. Finite post compliance is omitted.

    This accounts for finite anchors and edge reactions, but leaves joint
    stiffness, tube local buckling and complete-frame eigenmodes unresolved.
    """
    _positive(side_m=side_m, tension_N_m=tension_N_m, depth_m=depth_m,
              bay_m=bay_m, chord_radius_m=chord_radius_m, allowable_edge_bow_m=allowable_edge_bow_m)
    _nonnegative(joint_mass_fraction=joint_mass_fraction)
    if depth_m >= side_m or 4*chord_radius_m >= min(bay_m, depth_m):
        raise ValueError("Truss proportions are inconsistent")
    n = ceil(side_m/bay_m)
    bay = side_m/n
    bending = tension_N_m*side_m**2/8
    axial = tension_N_m*side_m/2
    chord_load = bending/(2*depth_m)+axial/4
    chord_area, chord_limit = _tube_area(chord_load, bay, chord_radius_m, material)
    chord_strength_buckling_area = chord_area
    chord_area = max(chord_area, 5*tension_N_m*side_m**4/
                     (384*material.youngs_modulus_Pa*depth_m**2*(allowable_edge_bow_m/2)))
    diagonal = sqrt(bay**2+depth_m**2)
    brace_load = tension_N_m*side_m/2 * diagonal/(2*depth_m)
    brace_area = brace_load/material.admissible_stress_Pa
    brace_area = max(brace_area, tension_N_m*side_m**2*diagonal**3/
                     (16*material.youngs_modulus_Pa*bay*depth_m**2*(allowable_edge_bow_m/2)))
    post_load = tension_N_m*bay/2
    post_area, post_limit = _tube_area(post_load, depth_m, chord_radius_m, material)
    rho = material.density_kg_m3
    chords = 16*side_m*chord_area*rho
    braces = 32*n*diagonal*brace_area*rho
    posts = 16*(n+1)*depth_m*post_area*rho
    joints = joint_mass_fraction*(chords+braces+posts)
    bow_bending = 5*tension_N_m*side_m**4/(384*material.youngs_modulus_Pa*chord_area*depth_m**2)
    shear_rigidity = 2*material.youngs_modulus_Pa*brace_area*bay*depth_m**2/diagonal**3
    bow_shear = tension_N_m*side_m**2/(8*shear_rigidity)
    return dict(mass_kg=chords+braces+posts+joints,
                components_kg=dict(chords=chords, x_braces=braces, posts=posts, joints=joints),
                bays_per_side=n, actual_bay_m=bay, max_bending_moment_N_m=bending,
                side_axial_compression_N=axial, chord_compression_N=chord_load,
                chord_area_m2=chord_area, chord_effective_stress_limit_Pa=chord_limit,
                chord_strength_and_buckling_area_m2=chord_strength_buckling_area,
                edge_bow_bending_m=bow_bending, edge_bow_shear_m=bow_shear,
                ideal_edge_bow_m=bow_bending+bow_shear,
                chord_wall_thickness_m=chord_area/(2*pi*chord_radius_m),
                chord_thin_wall_ratio=chord_area/(2*pi*chord_radius_m**2),
                post_effective_stress_limit_Pa=post_limit,
                outside_span_m=side_m+depth_m, half_depth_m=depth_m/2,
                unresolved=["tube local buckling", "joint and global-frame stiffness",
                            "deployment", "fatigue and micrometeoroid damage"])


def hoop_store(required_momentum_N_m_s: float, radius_m: float,
               rim_speed_m_s: float, *, density_kg_m3: float = 1600.,
               admissible_stress_Pa: float = 2e9, axes: int = 3,
               support_mass_fraction: float = .5) -> dict:
    """Reversible thin hoops, with full input capacity on each independent axis.

    H=m*r*v, I=m*r^2, E=m*v^2/2, hoop stress=rho*v^2. Two separated
    hoops per axis share capacity equally. Housing, bearings, spokes and motor
    mechanical mass are an explicit fractional allowance, not a solved design;
    no burst containment is established. In-plane axes protrude by radius_m.
    """
    _nonnegative(required_momentum_N_m_s=required_momentum_N_m_s,
                 support_mass_fraction=support_mass_fraction)
    _positive(radius_m=radius_m, rim_speed_m_s=rim_speed_m_s,
              density_kg_m3=density_kg_m3, admissible_stress_Pa=admissible_stress_Pa)
    if not isinstance(axes, int) or axes < 1:
        raise ValueError("axes must be a positive integer")
    mass = required_momentum_N_m_s/(radius_m*rim_speed_m_s)
    stress = density_kg_m3*rim_speed_m_s**2
    return dict(axes=axes, hoops=2*axes, rim_mass_per_axis_kg=mass,
                total_rim_mass_kg=axes*mass, support_mass_kg=axes*mass*support_mass_fraction,
                total_mass_kg=axes*mass*(1+support_mass_fraction),
                axial_inertia_per_axis_kg_m2=mass*radius_m**2,
                angular_speed_rad_s=rim_speed_m_s/radius_m,
                stored_energy_per_axis_J=mass*rim_speed_m_s**2/2,
                total_maximum_stored_energy_J=axes*mass*rim_speed_m_s**2/2,
                hoop_stress_Pa=stress, stress_utilization=stress/admissible_stress_Pa,
                stress_screen_passes=stress <= admissible_stress_Pa,
                out_of_plane_half_extent_m=radius_m,
                capacity_rule="full vector-envelope capacity independently on every axis",
                unresolved=["spokes and hoop support at zero spin", "bearings and motor detail",
                            "rotor burst containment", "gyroscopic coupling"])


@dataclass(frozen=True)
class TileDesign:
    side_m: float
    film_areal_mass_kg_m2: float
    target_total_areal_mass_kg_m2: float
    pressure_Pa: float
    angular_acceleration_rad_s2: float
    required_momentum_N_m_s: float
    required_torque_N_m: float
    orbit_period_s: float
    allowable_sag_m: float = 5.
    membrane_dynamic_factor: float = 2.
    substrate_areal_mass_kg_m2: float = .001
    substrate_density_kg_m3: float = 1400.
    substrate_admissible_stress_Pa: float = 100e6
    angular_speed_rad_s: float = 0.
    frame_depth_m: float = 100.
    frame_bay_m: float = 50.
    chord_radius_m: float = 2.
    frame_allowable_edge_bow_m: float = 5.
    structure: Material = field(default_factory=Material)
    rotor_radius_m: float = 100.
    rotor_rim_speed_m_s: float = 1000.
    rotor_density_kg_m3: float = 1600.
    rotor_admissible_stress_Pa: float = 2e9
    rotor_support_mass_fraction: float = .5
    trim_area_fraction: float = .25
    trim_hardware_areal_mass_kg_m2: float = .002
    trim_switch_J_m2: float = 5.
    trim_switches_per_orbit: int = 20
    trim_switch_duration_s: float = 20000.
    trim_hold_W_m2: float = 0.
    housekeeping_W: float = 1e4
    eclipse_s: float = 14400.
    slew_duration_s: float = 20000.
    recharge_s: float = 200000.
    bus_voltage_V: float = 10000.
    wire_loss_fraction: float = .02
    wire_resistivity_ohm_m: float = 2.8e-8
    wire_density_kg_m3: float = 2700.
    battery_specific_energy_J_kg: float = 720000.
    battery_depth_of_discharge: float = .8
    converter_efficiency: float = .9
    electronics_specific_power_W_kg: float = 2000.
    generator_specific_power_W_kg: float = 200.
    radiator_specific_heat_rejection_W_kg: float = 100.
    avionics_and_sensors_kg: float = 5000.
    mass_growth_fraction: float = .2

    def __post_init__(self):
        for key, value in asdict(self).items():
            if isinstance(value, (int, float)):
                _nonnegative(**{key: value})
        _positive(side_m=self.side_m, film_areal_mass_kg_m2=self.film_areal_mass_kg_m2,
                  target_total_areal_mass_kg_m2=self.target_total_areal_mass_kg_m2,
                  allowable_sag_m=self.allowable_sag_m,
                  substrate_areal_mass_kg_m2=self.substrate_areal_mass_kg_m2,
                  substrate_density_kg_m3=self.substrate_density_kg_m3,
                  substrate_admissible_stress_Pa=self.substrate_admissible_stress_Pa,
                  bus_voltage_V=self.bus_voltage_V, wire_loss_fraction=self.wire_loss_fraction,
                  battery_specific_energy_J_kg=self.battery_specific_energy_J_kg,
                  electronics_specific_power_W_kg=self.electronics_specific_power_W_kg,
                  generator_specific_power_W_kg=self.generator_specific_power_W_kg,
                  radiator_specific_heat_rejection_W_kg=self.radiator_specific_heat_rejection_W_kg,
                  recharge_s=self.recharge_s, orbit_period_s=self.orbit_period_s,
                  trim_switch_duration_s=self.trim_switch_duration_s)
        if not 0 < self.converter_efficiency <= 1 or not 0 < self.battery_depth_of_discharge <= 1:
            raise ValueError("Efficiencies and depth of discharge must lie in (0, 1]")
        if not 0 <= self.trim_area_fraction <= 1 or self.wire_loss_fraction >= 1:
            raise ValueError("Invalid area or wiring-loss fraction")
        if self.membrane_dynamic_factor < 1 or self.rotor_radius_m > self.side_m/10:
            raise ValueError("Invalid dynamic factor or overlapping rotor station geometry")
        if not isinstance(self.trim_switches_per_orbit, int):
            raise ValueError("trim_switches_per_orbit must be a nonnegative integer")


def tile_budget(design: TileDesign) -> dict:
    """Close the mass ledger for explicit idealizations; report unmet targets.

    The load bound is qsolar + sigma*(|alpha|+|omega|^2)*L/sqrt(2), a
    conservative norm bound including combined-axis rotation. Sigma includes
    the full local trim coating wherever it is present. Full solar pressure is
    conservative in the translating frame: uniform membrane acceleration can
    cancel some of it. A factor of two is an undamped single-mode step screen,
    not a bound for arbitrary resonant forcing. Frame edges are assumed fixed.

    Power scenario: simultaneous peak torque on all three axes, six distributed
    local motor buses, and battery capacity for one full zero-to-maximum energy
    swing on all axes plus every orbit's trim switches in eclipse. Switch energy
    is specified at trim-device terminals; zero hold power is an optimistic
    bistable-device target. Peak trim power is at least both the nominal event
    rate and the orbit/eclipse envelope averages. These independent conservative
    envelopes do not establish a realizable switching chronology. Local motor
    buses avoid sending peak
    regenerative motor power through a ten-kilometre trunk; the trunk is still
    conservatively sized for total motor peak power.
    """
    d = design
    area = d.side_m**2
    sigma = (d.film_areal_mass_kg_m2+d.substrate_areal_mass_kg_m2+
             (d.trim_hardware_areal_mass_kg_m2 if d.trim_area_fraction else 0.))
    inertial_pressure = sigma*(d.angular_acceleration_rad_s2+d.angular_speed_rad_s**2)*d.side_m/sqrt(2)
    pressure = d.pressure_Pa+inertial_pressure
    # A quiescent sheet still needs prescribed prestress; zero load is not a
    # meaningful automatic tension-sizing case.
    _positive(total_pressure_bound_Pa=pressure)
    tension = membrane_coefficient()*pressure*area*d.membrane_dynamic_factor/d.allowable_sag_m
    membrane = square_membrane(d.side_m, pressure, tension, sigma)
    substrate_thickness = d.substrate_areal_mass_kg_m2/d.substrate_density_kg_m3
    stress = tension/substrate_thickness
    membrane.update(pressure_bound_Pa=pressure, solar_pressure_Pa=d.pressure_Pa,
                    slew_pressure_bound_Pa=inertial_pressure, tension_N_m=tension,
                    dynamic_sag_screen_m=membrane['centre_sag_m']*d.membrane_dynamic_factor,
                    substrate_stress_Pa=stress,
                    substrate_stress_utilization=stress/d.substrate_admissible_stress_Pa,
                    load_path_assumption='added substrate carries prestress; optical-film strain accommodation unresolved',
                    slew_duration_over_first_period=d.slew_duration_s/membrane['fundamental_period_s'])
    frame = perimeter_truss(d.side_m, tension, d.frame_depth_m, d.frame_bay_m,
                            d.chord_radius_m, d.structure,
                            allowable_edge_bow_m=d.frame_allowable_edge_bow_m)
    rotors = hoop_store(d.required_momentum_N_m_s, d.rotor_radius_m, d.rotor_rim_speed_m_s,
                       density_kg_m3=d.rotor_density_kg_m3,
                       admissible_stress_Pa=d.rotor_admissible_stress_Pa,
                       support_mass_fraction=d.rotor_support_mass_fraction)
    peak_mechanical = 3*d.required_torque_N_m*rotors['angular_speed_rad_s']
    trim_area = d.trim_area_fraction*area
    switch_energy_J = trim_area*d.trim_switch_J_m2
    switch_cycle_energy_J = switch_energy_J*d.trim_switches_per_orbit
    switch_average_W = switch_cycle_energy_J/d.orbit_period_s
    switch_nominal_peak_W = switch_energy_J/d.trim_switch_duration_s if d.trim_switches_per_orbit else 0.
    switch_eclipse_peak_W = switch_cycle_energy_J/d.eclipse_s if d.eclipse_s else 0.
    switch_peak_W = max(switch_nominal_peak_W, switch_average_W, switch_eclipse_peak_W)
    hold_W = trim_area*d.trim_hold_W_m2
    trim_peak_W = switch_peak_W+hold_W
    peak_electrical = ((peak_mechanical+trim_peak_W)/d.converter_efficiency+d.housekeeping_W)/(1-d.wire_loss_fraction)
    energy_swing = rotors['total_maximum_stored_energy_J']
    trim_eclipse_energy_J = switch_cycle_energy_J+hold_W*d.eclipse_s
    reserve_J = ((energy_swing+trim_eclipse_energy_J)/d.converter_efficiency+
                 d.housekeeping_W*d.eclipse_s)/(1-d.wire_loss_fraction)
    battery_nominal_J = reserve_J/d.battery_depth_of_discharge
    # Full-cycle charging includes the whole eclipse reserve, not just rotor
    # losses; regenerative recovery could lower this but is not credited.
    # The reserve already contains the complete cycle's switch energy. Do not
    # add the average switch energy to recharge a second time. Take the larger
    # of recharge-plus-continuous demand and the steady average trim demand.
    continuous_bus_W = (d.housekeeping_W+hold_W/d.converter_efficiency)/(1-d.wire_loss_fraction)
    recharge_generation_W = continuous_bus_W+reserve_J/(d.recharge_s*d.converter_efficiency*(1-d.wire_loss_fraction))
    mean_trim_bus_W = (d.housekeeping_W+(hold_W+switch_average_W)/d.converter_efficiency)/(1-d.wire_loss_fraction)
    generation_W = max(recharge_generation_W, mean_trim_bus_W)
    heat_W = peak_electrical-peak_mechanical
    round_trip_m = d.side_m
    branch_power_W = peak_electrical/6
    wire_area = d.wire_resistivity_ohm_m*round_trip_m*branch_power_W/(d.bus_voltage_V**2*d.wire_loss_fraction)
    wire_mass = 6*round_trip_m*wire_area*d.wire_density_kg_m3
    components = dict(film=d.film_areal_mass_kg_m2*area,
                      load_bearing_substrate=d.substrate_areal_mass_kg_m2*area,
                      perimeter_frame=frame['mass_kg'],
                      rotor_rims=rotors['total_rim_mass_kg'],
                      rotor_mechanical_support=rotors['support_mass_kg'],
                      optical_trim=d.trim_area_fraction*area*d.trim_hardware_areal_mass_kg_m2,
                      motor_power_electronics=peak_electrical/d.electronics_specific_power_W_kg,
                      power_generation=generation_W/d.generator_specific_power_W_kg,
                      eclipse_energy_storage=battery_nominal_J/d.battery_specific_energy_J_kg,
                      motor_heat_rejection=heat_W/d.radiator_specific_heat_rejection_W_kg,
                      power_conductors=wire_mass,
                      avionics_sensors=d.avionics_and_sensors_kg)
    nominal_mass = sum(components.values())
    components['mass_growth_allowance'] = nominal_mass*d.mass_growth_fraction
    total_mass = sum(components.values())
    target_mass = d.target_total_areal_mass_kg_m2*area

    # Locked principal inertia: film/substrate uniform square; frame uniformly
    # around perimeter. Three rotor axis-pairs at ±(a,0), ±(0,a), ±(a,a),
    # a=L/2, are attached directly to the finite frame and include hoop
    # self-inertia. Housings at the same points. Remaining
    # electronics, supplies, trim, harness and reserve mass spread uniformly
    # along the perimeter as an explicit conservative placement assumption.
    sheet_mass = components['film']+components['load_bearing_substrate']+components['optical_trim']
    per_axis_mass = rotors['rim_mass_per_axis_kg']
    a = d.side_m/2
    rotor_self = 2*per_axis_mass*d.rotor_radius_m**2  # sum over x/y/z axes
    station_mass = per_axis_mass*(1+d.rotor_support_mass_fraction)
    rotor_diag = [2*station_mass*a*a+rotor_self,
                  2*station_mass*a*a+rotor_self,
                  4*station_mass*a*a+rotor_self]
    perimeter_mass = total_mass-sheet_mass-rotors['total_mass_kg']
    inertia = [sheet_mass*area/12+perimeter_mass*area/6+rotor_diag[0],
               sheet_mass*area/12+perimeter_mass*area/6+rotor_diag[1],
               sheet_mass*area/6+perimeter_mass*area/3+rotor_diag[2]]
    # The z-axis pair at ±(a,a) has product of inertia -m*x*y; preserve it.
    off_diagonal = -station_mass*a*a
    tensor = [[inertia[0], off_diagonal, 0.], [off_diagonal, inertia[1], 0.], [0., 0., inertia[2]]]
    uniform_reference = target_mass*area/12
    flags = []
    if total_mass > target_mass:
        flags.append('mass target exceeded')
    if stress > d.substrate_admissible_stress_Pa:
        flags.append('membrane substrate stress exceeds assumed admissible value')
    if not rotors['stress_screen_passes']:
        flags.append('rotor hoop stress exceeds assumed admissible value')
    if frame['chord_thin_wall_ratio'] > .1:
        flags.append('thin-wall chord approximation invalid')
    flags.extend(['assembled inertia requires a new attitude and orbit calculation',
                  'frame flexibility and membrane resonances unresolved',
                  'motor attachment strength and complete out-of-plane frame slew loads unresolved',
                  'optical-film/substrate strain compatibility and modified spectrum unresolved',
                  'rotor containment, detailed hardware and deployment unresolved'])
    return dict(design=asdict(d), evidence='conditional analytic sizing screen; no validated tile design',
                membrane=membrane, frame=frame, rotors=rotors,
                power=dict(peak_mechanical_W=peak_mechanical, peak_electrical_W=peak_electrical,
                           battery_nominal_J=battery_nominal_J, generation_W=generation_W,
                           peak_waste_heat_W=heat_W, branch_wire_area_m2=wire_area,
                           wire_round_trip_length_m=round_trip_m, independent_branches=6,
                           electronics_scope='motor and trim converters',
                           recharge_generation_W=recharge_generation_W,
                           mean_trim_bus_W=mean_trim_bus_W),
                optical_trim=dict(one_axis_torque_upper_bound_N_m=d.pressure_Pa*area*d.side_m/4*d.trim_area_fraction,
                                  rule='all trim area on one axis at quarter-side lever, full pressure modulation; axes share area',
                                  qualified_authority=False,
                                  switch_energy_per_event_J=switch_energy_J,
                                  switch_energy_per_orbit_J=switch_cycle_energy_J,
                                  average_switch_power_W=switch_average_W,
                                  nominal_event_switch_power_W=switch_nominal_peak_W,
                                  eclipse_envelope_switch_power_W=switch_eclipse_peak_W,
                                  peak_switch_power_W=switch_peak_W, continuous_hold_power_W=hold_W,
                                  eclipse_device_energy_J=trim_eclipse_energy_J,
                                  switch_window_sum_over_orbit=d.trim_switches_per_orbit*d.trim_switch_duration_s/d.orbit_period_s,
                                  electrical_evidence='bistable device proposal, no qualified device; all cycle switches reserved in eclipse',
                                  qualified_electrical_device=False),
                motor_attachment=dict(force_couple_per_station_N=d.required_torque_N_m/(2*d.frame_depth_m),
                                      support_rule='two stations per axis attach to perimeter; torque shared equally across frame-depth couple',
                                      detailed_joint_design=False),
                components_kg=components, total_mass_kg=total_mass,
                total_areal_mass_kg_m2=total_mass/area, target_mass_kg=target_mass,
                mass_margin_kg=target_mass-total_mass, mass_ratio_to_target=total_mass/target_mass,
                locked_inertia_kg_m2=tensor,
                inertia_diagonal_ratio_to_target_uniform_square=[inertia[0]/uniform_reference,
                                                                  inertia[1]/uniform_reference,
                                                                  inertia[2]/(2*uniform_reference)],
                geometry=dict(outside_frame_span_m=frame['outside_span_m'],
                              frame_out_of_plane_half_extent_m=d.frame_depth_m/2,
                              rotor_out_of_plane_half_extent_m=d.rotor_radius_m,
                              conservative_pair_extra_envelope_m=2*max(d.rotor_radius_m,
                                  sqrt(2)*d.frame_depth_m/2, d.allowable_sag_m),
                              envelope_rule='sufficient whole-tile bound only; compact rotor stations require swept-component collision analysis',
                              rotor_station_centres_rule='axis pairs attached to frame at ±(L/2,0), ±(0,L/2), ±(L/2,L/2)'),
                findings=flags, validated_design=False)
