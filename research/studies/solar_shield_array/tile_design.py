"""Tile hardware, local faults and fleet service for research/array-industry.

Run: OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.tile_design

This is a reproducible preliminary design assessment. It can reject a specified
hardware concept or service policy; it cannot certify the complete ring fleet.
It reads the branch's rolled-law products and keeps all earlier results intact.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import time

import numpy as np
import scipy

from engineering.array_tile import TileDesign, tile_budget
from engineering.array_momentum import angular_velocity_from_frames, integrate_rotor_momentum
from engineering.array_service import ServiceCase, coverage_multiplicity, missing_ray_fraction, simulate_service
from protection.dynamics.tile_failure import FailureConfig, simulate_failure, rotation_bound_window_s
from shared import constants as K
from shared.provenance import constants_used
from .attitude_schemes import body_frames, rolled

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/tile_design.json'
FILES = ['research/studies/solar_shield_array/tile_design.py',
         'engineering/array_tile.py', 'engineering/array_momentum.py', 'engineering/array_service.py',
         'protection/dynamics/tile_failure.py',
         'research/studies/solar_shield_array/attitude_schemes.py', 'protection/dynamics/ring_bundle.py',
         'protection/dynamics/optical.py', 'shared/provenance.py']
INPUTS = {
    'attitude': ('research/studies/solar_shield_array/results/attitude_schemes.json', 'terluna.research.attitude-schemes/1'),
    'layout': ('research/studies/solar_shield_array/results/ring_layout_rolled.json', 'terluna.research.ring-layout/1'),
    'heat': ('research/studies/solar_shield_array/results/array_heat.json', 'terluna.research.array-heat/1'),
    'glow': ('research/studies/solar_shield_array/results/glow.json', 'terluna.research.fleet-glow/1'),
    'films': ('protection/spectra/annulus_film.json', 'terluna.protection.annulus-film/1'),
    'sources': ('research/studies/solar_shield_array/tile_sources.json', 'terluna.research.tile-hardware-sources/1'),
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs():
    data = {}
    for key, (name, schema) in INPUTS.items():
        data[key] = json.loads((ROOT/name).read_text())
        if schema and data[key].get('schema') != schema:
            raise ValueError(f'Unexpected schema for {name}')
    if data['layout']['variant']['law'] != 'rolled_80_deg' or data['layout']['smoke']:
        raise ValueError('Expected the full rolled-law layout')
    return data


def ideal_attitude(radius_m, samples=1441):
    """One circular lunar orbit with the branch's 80-degree roll command.

    This isolates inertia and correct momentum transport. It is not the stored
    eccentric/perturbed orbit, and does not replay its shadow torques.
    """
    n = math.sqrt(K.MOON_GM/radius_m**3)
    t = np.linspace(0, 2*math.pi/n, samples)
    phase = n*t
    states = np.column_stack([radius_m*np.cos(phase), radius_m*np.sin(phase), np.zeros(samples),
                              -radius_m*n*np.sin(phase), radius_m*n*np.cos(phase), np.zeros(samples)])
    roll = rolled(phase, math.radians(80), 25, 30, -1)
    frames = body_frames(states, math.radians(1.5), roll)
    omega = angular_velocity_from_frames(t, frames)
    alpha = np.gradient(omega, t, axis=0, edge_order=2)
    return t, states, frames, omega, alpha


def momentum_screen(inertia, radius_m, samples=1441):
    """Required rotor motion with gravity-gradient torque and no optical unload.

    Zero initial rotor momentum. This is a declared degraded-mode comparison,
    not a guaranteed upper bound on all radiation/shadow torques or an optimum.
    """
    t, states, frames, omega, _ = ideal_attitude(radius_m, samples)
    inertia = np.asarray(inertia)
    radial_body = np.einsum('tji,tj->ti', frames, states[:, :3]/radius_m)
    gravity = 3*K.MOON_GM/radius_m**3*np.cross(radial_body, radial_body@inertia.T)
    p = integrate_rotor_momentum(t, frames, inertia, gravity, body_angular_velocity_rad_s=omega)
    return dict(samples=samples, radius_m=radius_m, duration_s=float(t[-1]),
                initial_relative_rotor_momentum_N_m_s=[0, 0, 0],
                peak_rotor_momentum_N_m_s=p['peak_rotor_momentum_N_m_s'],
                peak_axis_momentum_N_m_s=np.max(np.abs(p['rotor_momentum_body_N_m_s']), axis=0).tolist(),
                peak_axis_motor_torque_N_m=np.max(np.abs(p['rotor_momentum_rate_body_N_m']), axis=0).tolist(),
                peak_bearing_gyroscopic_torque_N_m=p['peak_rotor_gyroscopic_torque_N_m'],
                peak_body_reaction_torque_N_m=p['peak_rotor_reaction_torque_N_m'],
                derivative_balance_residual_N_m=p['max_derivative_balance_residual_N_m'],
                final_rotor_momentum_N_m_s=p['rotor_momentum_body_N_m_s'][-1].tolist(),
                external_model='lunar point-mass gravity gradient, zero radiation and unloading torque',
                evidence='rigid assembled-inertia requirement on an ideal circular roll; no control allocation or flexible dynamics')


def size_momentum_hardware(design, radius_m, margin=1.2, maximum_iterations=40):
    """Resize the rotor/power mass and re-evaluate its inertia until consistent.

    The fixed-point only closes this rigid, gravity-only requirement screen.
    It gives no credit for optical unloading and is not an operating controller.
    Never shrink the initial allocation: its historical radiation load is not
    present in the new ideal history. Report non-convergence without accepting it.
    """
    trace = []
    d = design
    for _ in range(maximum_iterations):
        budget = tile_budget(d)
        m = momentum_screen(budget['locked_inertia_kg_m2'], radius_m)
        required_h = max(m['peak_axis_momentum_N_m_s'])*margin
        required_torque = max(m['peak_axis_motor_torque_N_m'])*margin
        trace.append(dict(total_mass_kg=budget['total_mass_kg'],
                          capacity_per_axis_N_m_s=d.required_momentum_N_m_s,
                          motor_torque_per_axis_N_m=d.required_torque_N_m,
                          required_axis_H_with_margin_N_m_s=required_h,
                          required_axis_torque_with_margin_N_m=required_torque))
        converged = (required_h <= d.required_momentum_N_m_s*1.001 and
                     required_torque <= d.required_torque_N_m*1.001)
        if converged:
            break
        d = replace(d, required_momentum_N_m_s=max(d.required_momentum_N_m_s, required_h),
                    required_torque_N_m=max(d.required_torque_N_m, required_torque))
    budget['assembled_momentum_screen'] = m
    budget['rigid_sizing_iteration'] = dict(converged=converged, tolerance_fraction=.001,
                                           requirement_margin=margin, trace=trace,
                                           periodic_controller_demonstrated=False,
                                           scope='one-orbit rigid requirement with nonzero terminal rotor momentum; external unloading and flexible control unresolved')
    budget['momentum_capacity_utilization_in_screen'] = max(m['peak_axis_momentum_N_m_s'])/budget['design']['required_momentum_N_m_s']
    budget['motor_torque_utilization_in_screen'] = max(m['peak_axis_motor_torque_N_m'])/budget['design']['required_torque_N_m']
    return budget


def hardware_cases(data):
    attitude = data['attitude']
    inherited = attitude['schemes']['balanced_rolled_80_deg']['by_areal_mass']['26_g_m2']['store_N_m_s']['trim_0.25']
    radius = data['layout']['design']['middle_radius_km']*1000
    command_t, _, _, omega, alpha = ideal_attitude(radius)
    alpha_max = float(np.linalg.norm(alpha, axis=1).max())
    omega_max = float(np.linalg.norm(omega, axis=1).max())
    sail = data['glow']['results']['recommended']['costs']['sail_loading']
    films = {'annulus': next(x for x in data['films']['films']
                            if x['titania_um'] == .1 and x['silica_um'] == 4 and x['coated']),
             'window': data['films']['window_stack']}
    cases = {}
    for name, film in films.items():
        target = sail[f'{name}_tile_g_m2']/1000
        pressure_coefficient = sail['annulus_coefficient' if name == 'annulus' else 'window_dynamic_coefficient']
        # The inherited H is only an initial allocation. New inertia is checked
        # separately; no linear rescaling of its old uniform-square load is accepted.
        base = TileDesign(side_m=10000, film_areal_mass_kg_m2=film['areal_mass_g_m2']/1000,
                          target_total_areal_mass_kg_m2=target,
                          pressure_Pa=pressure_coefficient*K.SOLAR_CONSTANT/K.SPEED_OF_LIGHT,
                          angular_acceleration_rad_s2=alpha_max,
                          angular_speed_rad_s=omega_max,
                          orbit_period_s=float(command_t[-1]),
                          required_momentum_N_m_s=inherited, required_torque_N_m=1e6)
        for suffix, settings in [('reference', {}), ('shallow_frame', {'frame_depth_m': 50.}),
                                 ('looser_flatness', {'allowable_sag_m': 20.}),
                                 ('bow_50m', {'frame_allowable_edge_bow_m': 50.}),
                                 ('wide_rotors', {'rotor_radius_m': 400.})]:
            d = replace(base, **settings)
            result = size_momentum_hardware(d, radius)
            result['sail_acceleration_over_inherited_target'] = target/result['total_areal_mass_kg_m2']
            result['nominal_shingle_parallel_frame_clearance_m'] = 8500*math.sin(math.radians(1.5))-d.frame_depth_m
            result['nominal_shingle_frame_gate_passes'] = result['nominal_shingle_parallel_frame_clearance_m'] >= 150
            result['optical_qualification'] = 'Added substrate, trim, frame and power hardware not included in inherited spectra/glow'
            result['pressure_convention'] = dict(inherited_dynamic_coefficient=pressure_coefficient,
                                                  bare_film_coefficient=film['pressure_coefficient'],
                                                  rule='use the glow/dynamics redirector convention; bare-film spectra alone are not the flight force')
            cases[f'{name}_{suffix}'] = result
        # A geometry/structure scaling comparison, not a resized ring layout.
        # H scales as L^3 here to retain a radiation-torque-like allocation;
        # pure inertia would scale as L^4 at fixed areal mass and angular history.
        for side in (1000., 3000.):
            ratio = side/base.side_m
            d = replace(base, side_m=side, frame_depth_m=side/100, frame_bay_m=side/200,
                        chord_radius_m=side/5000, frame_allowable_edge_bow_m=side/2000,
                        allowable_sag_m=side/2000, rotor_radius_m=side/100,
                        required_momentum_N_m_s=inherited*ratio**3,
                        required_torque_N_m=1e6*ratio**3)
            result = size_momentum_hardware(d, radius)
            result['interpretation'] = 'Size sensitivity only: frame/sag geometry scales with side, H and torque initially as L^3; avionics/housekeeping fixed per tile. Requires a new formation, population and optical model.'
            cases[f'{name}_side_{side:g}m'] = result
    # Numerical convergence is distinct from the physical adequacy of the model.
    ref = cases['annulus_reference']['locked_inertia_kg_m2']
    coarse = momentum_screen(ref, radius, 721)
    fine = momentum_screen(ref, radius, 2881)
    baseline = cases['annulus_reference']['assembled_momentum_screen']
    return dict(inherited_momentum_allocation_N_m_s=inherited,
                inherited_limit='Old store() integrates body-axis torque components without transport; this is an allocation, not a validated capacity requirement.',
                initial_torque_allocation_per_axis_N_m=1e6,
                ideal_command_peak_angular_acceleration_rad_s2=alpha_max,
                ideal_command_peak_angular_speed_rad_s=float(np.linalg.norm(omega, axis=1).max()),
                current_areal_mass_allocations_g_m2={name: sail[f'{name}_tile_g_m2'] for name in films},
                window_bare_film_excess_over_complete_allocation_g_m2=films['window']['areal_mass_g_m2']-sail['window_tile_g_m2'],
                variants=cases,
                momentum_convergence=dict(coarse=coarse, baseline=baseline, fine=fine,
                    relative_peak_H_change=abs(fine['peak_rotor_momentum_N_m_s']/baseline['peak_rotor_momentum_N_m_s']-1)),
                evidence='One explicit perimeter-truss architecture and declared hardware targets; exceeding a target rejects this architecture under these assumptions, not every sail architecture.')


def local_failures(data):
    radius = data['layout']['design']['middle_radius_km']*1000
    cases = {}
    for thickness in (1., 50., 100.):
        for acceleration in (1e-6, 1e-5, 1e-4):
            c = FailureConfig(radius_m=radius, envelope_thickness_m=thickness,
                              residual_acceleration_m_s2=(acceleration, 0, 0), duration_s=21600)
            cases[f'depth_{thickness:g}_a_{acceleration:g}'] = simulate_failure(c)
    for delay in (600., 1800., 3600.):
        c = FailureConfig(radius_m=radius, envelope_thickness_m=50.,
                          residual_acceleration_m_s2=(1e-5, 0, 0), response_delay_s=delay,
                          backup_acceleration_cap_m_s2=1e-4, duration_s=21600)
        cases[f'backup_after_{delay:g}s'] = simulate_failure(c)
    ref = FailureConfig(radius_m=radius, envelope_thickness_m=50., duration_s=21600,
                        residual_acceleration_m_s2=(1e-5, 0, 0), max_step_s=10.)
    fine = simulate_failure(ref)
    coarse = cases['depth_50_a_1e-05']
    initial = coarse['initial_minimum_surface_clearance_m']
    return dict(cases=cases, timestep_check=dict(fine=fine,
                safety_time_difference_s=abs(coarse['first_safety_crossing_s']-fine['first_safety_crossing_s']),
                contact_time_difference_s=abs(coarse['first_envelope_contact_s']-fine['first_envelope_contact_s'])),
                attitude_guarantee_loss=[dict(bounded_error_rate_deg_s=rate,
                    envelope_thickness_m=50.,
                    bound_lost_s=rotation_bound_window_s(initial, 150, 10000, math.radians(rate), 50.))
                    for rate in (.0001, .001, .01)],
                evidence='Finite local envelopes and continuous drift/recovery, then conditional population accounting. No spacecraft removal or replacement insertion trajectory is demonstrated.',
                hardware_coupling='1 m is an optimistic thin-envelope sensitivity; 50/100 m are explicit frame-depth cases. Rotor placements and moving parts require separate swept geometry.')


def fleet_service(data):
    # Retain the prior inventory; the new structural mass does not justify quietly
    # changing ring geometry, radii, optical area or the global tile count.
    population = round(max(data['heat']['fleet']['tiles']))
    base = ServiceCase(population=population)
    scenarios = {
        'local_spares_1day': base,
        'plant_delivery_27days': replace(base, failed_slot_response_days=27.),
        'capacity_below_renewal': replace(base, capacity_over_planned=.8),
        'dispatch_stops_14days': replace(base, outage_start_days=120., outage_duration_days=14.),
        'correlated_10000_slots': replace(base, correlated_failure_count=10000),
        'hazard_0.001_per_year': replace(base, failure_rate_per_year=.001),
        'hazard_0.1_per_year': replace(base, failure_rate_per_year=.1),
    }
    height_pitch = data['layout']['design']['height_step_m']
    multiplicity = coverage_multiplicity(8900., 9890., 8500., height_pitch)
    results = {}
    for name, config in scenarios.items():
        p = simulate_service(config)
        # Average the nonlinear optical comparison over the actual daily trace,
        # not f(mean(q)); retain the conservative all-single-cover bound too.
        q = np.array([row['missing']/population for row in p['daily_trace']])
        p['ideal_independent_peak_missing_ray_fraction'] = missing_ray_fraction(p['peak_missing_fraction'], multiplicity)
        p['ideal_independent_daily_sample_mean_missing_ray_fraction'] = float(np.mean([missing_ray_fraction(x, multiplicity) for x in q]))
        if config.correlated_failure_count:
            p['ideal_independent_peak_missing_ray_fraction'] = None
            p['ideal_independent_daily_sample_mean_missing_ray_fraction'] = None
        results[name] = p
    return dict(population=population, scenarios=results,
                ideal_service_lattice=dict(optical_width_m=8900., optical_height_m=9890., pitch_m=[8500.,height_pitch],
                                          coverage_multiplicity=multiplicity),
                interpretation='Uniform independent failure optical comparison only, baseline seams/holes excluded. Correlated-case optical estimates are withheld. Every queue assumes failed hardware is contained; local fault response must establish that separately.',
                rate_evidence='0.001, 0.01 and 0.1 failures per tile-year are sensitivities, not measured reliability; 15-year renewal is an assumed policy.')


def industry_comparison(data, hardware, service):
    """Cost of copying the rejected reference design, not a new fleet solution."""
    share = data['heat']['fleet']['window_ring_share']
    annulus = hardware['variants']['annulus_reference']
    window = hardware['variants']['window_reference']
    mean_mass = (1-share)*annulus['total_mass_kg']+share*window['total_mass_kg']
    n = service['population']
    base = service['scenarios']['local_spares_1day']
    jobs_per_day = base['planned_per_day']+base['expected_failure_per_day_at_full_population']
    return dict(window_share=share, mean_reference_tile_mass_kg=mean_mass,
                inventory_if_copied_Gt=n*mean_mass/1e12,
                whole_tile_exchange_Gt_per_year=jobs_per_day*K.JULIAN_YEAR_DAYS*mean_mass/1e12,
                whole_tile_exchange_kg_s=jobs_per_day*mean_mass/K.JULIAN_DAY,
                interpretation='Counterfactual hardware mass at the old population and class share. Not a re-flown shield or a film-only factory load. Whole-tile exchange is a conservative service policy; component reuse and manufacturing processes are not designed.')


def run():
    started = time.process_time()
    data = load_inputs()
    hardware = hardware_cases(data)
    failures = local_failures(data)
    service = fleet_service(data)
    return dict(schema='terluna.research.array-tile-design/1',
                producer=dict(files={name: digest(ROOT/name) for name in FILES}, constants=constants_used(FILES)),
                inputs={name: dict(sha256=digest(ROOT/name), schema=data[key].get('schema'))
                        for key, (name, _) in INPUTS.items()},
                evidence='Preliminary analytic hardware design, conserved rigid-body momentum on an ideal orbit, local circular-orbit fault dynamics, and aggregate whole-population service queues. Not a validated tile, integrated photon controller, qualified reliability estimate, or full-fleet trajectory simulation.',
                reading_rule='All suffixes carry SI units except explicitly named days, years, degrees and g/m2. Hardware variants are alternatives; never add their masses. Local cases stop physical claims at contact. Service counts conditional on fault containment and finite transit; their independent flat-lattice optical fractions are not total shield transmission.',
                baseline=dict(branch='research/array-industry', attitude='balanced_rolled_80_deg',
                              orbit_model='circular 20,000 km for new local/dynamics screens; rolled product remains an inherited design input',
                              previous_products_modified=False),
                hardware=hardware, local_failures=failures, fleet_service=service,
                industry_comparison=industry_comparison(data, hardware, service),
                checks=dict(momentum_relative_peak_refinement=hardware['momentum_convergence']['relative_peak_H_change'],
                            fault_safety_event_refinement_s=failures['timestep_check']['safety_time_difference_s'],
                            fault_contact_event_refinement_s=failures['timestep_check']['contact_time_difference_s'],
                            queue_balance=all(all(x==0 for x in s['balance_residual'].values()) for s in service['scenarios'].values())),
                resources=dict(cpu_s=time.process_time()-started,
                               max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                               thread_environment={k: os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS')},
                               versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args()
    product = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(product, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(args.output), checks=product['checks'], resources=product['resources']), indent=2))


if __name__ == '__main__':
    main()
