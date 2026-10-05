"""Bounded primary-magnet comparison: vacuum fields and engineering screens.

Reads saved atmospheric and tile products; does not propagate a new fleet or
solve a magnetosphere. Run from the repository root with numerical threads=1.
"""
import hashlib
import json
from pathlib import Path
import resource
import signal
import time
import numpy as np
from scipy.spatial.transform import Rotation

from shared import constants as K
from shared.provenance import constants_used
from protection.magnetic_fields import field, loop_field, sphere, regional_loops, mutual_load, mutual_inductance, loop_samples
from protection.dynamics.ephemeris import Ephemeris
from protection.dynamics.model import geometry, relative_gravity
from engineering.large_coils import winding, renewal, upstream_width
from engineering.electromagnetic import collector, radiator, plasma, storage

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
RAW = ROOT/'research/runs/solar_shield_array/magnetic_architecture'
REFERENCE_MOMENT = 1.5e21
SMALL_MOMENT = 10**19.5
SOURCES = ['protection/magnetic_fields.py', 'engineering/large_coils.py',
    'engineering/electromagnetic.py', 'protection/dynamics/ephemeris.py',
    'protection/dynamics/model.py', 'protection/dynamics/optical.py',
    'research/studies/solar_shield_array/magnetic_architecture_run.py',
    'research/studies/solar_shield_array/magnetic_architecture_sources.json']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plain(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    raise TypeError(type(value))


def summarize(loops, include_mutual=True):
    hardware = [winding(l['radius'], l['current']) for l in loops]
    mass = sum(h['screened_mass_with_generation_kg'] for h in hardware)
    self_energy = sum(h['self_energy_J'] for h in hardware)
    interactions = []; mutual_energy = 0.
    for i, source in enumerate(loops):
        for j, target in enumerate(loops[i+1:], i+1):
            if not include_mutual: continue
            v = mutual_load(source, target, 512)
            reverse = mutual_load(target, source, 512)
            inductance = mutual_inductance(source, target, 256)
            coarse = mutual_inductance(source, target, 128)
            mutual_energy += inductance*source['current']*target['current']
            angular = v['torque_N_m']+reverse['torque_N_m']+np.cross(target['centre']-source['centre'], v['force_N'])
            interactions.append(dict(source=i, target=j, **v, mutual_H=inductance,
                mutual_refinement_relative=abs(inductance-coarse)/max(abs(inductance), 1e-30),
                force_reciprocity_relative=float(np.linalg.norm(v['force_N']+reverse['force_N'])/max(np.linalg.norm(v['force_N']), 1.)),
                angular_reciprocity_relative=float(np.linalg.norm(angular)/max(np.linalg.norm(v['torque_N_m']), np.linalg.norm(reverse['torque_N_m']), np.linalg.norm(target['centre']-source['centre'])*np.linalg.norm(v['force_N']), 1.))))
    power = sum(h['refrigerator_W'] for h in hardware)
    dark_seconds = .5*K.SYNODIC_MONTH_DAYS*K.JULIAN_DAY
    energy = self_energy+mutual_energy
    dark = power*dark_seconds
    return dict(loop_count=len(loops), hardware=hardware,
        conductor_kg=sum(h['conductor_mass_kg'] for h in hardware),
        ideal_tensile_kg=sum(h['ideal_tensile_mass_kg'] for h in hardware),
        screened_mass_kg=mass,
        refrigeration_W=power,
        collector_area_m2=sum(h['generation']['area_m2'] for h in hardware),
        self_energy_J=self_energy, mutual_energy_J=mutual_energy,
        total_magnetic_energy_J=self_energy+mutual_energy,
        mutual_loads=interactions,
        half_lunar_night=dict(duration_s=dark_seconds, delivered_energy_J=dark,
            battery=storage(dark, power),
            field_reserve_initial_current_multiplier=float(np.sqrt(1+dark/(.95*energy))),
            scope='Unshadowed solar generation in main mass; night requires a grid, other firm source, or additional storage. Field reserve is one-pass extraction above nominal current, before its extra mass/cooling, converter and fault design.'),
        replacement=[renewal(mass, energy, life) for life in [20., 100., 1000.]])


def rescale(loops, ratio):
    return [{**l, 'current': l['current']*ratio} for l in loops]


def obstacle_field(radius, obstacle, samples=361):
    theta = np.linspace(0, np.pi/2, samples)
    p = obstacle*np.c_[np.sin(theta), np.zeros(samples), np.cos(theta)]
    b = np.linalg.norm(loop_field(p, radius), axis=1)
    return float(b.min())


def regional_case(name, radius, count=4, historical=False):
    loops = regional_loops(radius, REFERENCE_MOMENT, count, historical)
    low = np.linalg.norm(field(sphere(4*K.MOON_RADIUS), loops), axis=1)
    fine_points = sphere(4*K.MOON_RADIUS, 121, 240)
    fine = np.linalg.norm(field(fine_points, loops), axis=1)
    at3 = np.linalg.norm(field(sphere(3*K.MOON_RADIUS), loops), axis=1)
    outage = [float(np.min(np.linalg.norm(field(fine_points, loops[:i]+loops[i+1:]), axis=1))) for i in range(count)]
    routes = np.concatenate([loop_samples(l)[0] for l in loops])
    extents = np.linalg.norm(routes, axis=1)-K.MOON_RADIUS
    row = dict(name=name, radius_m=radius, loop_count=count, loops=loops,
        reference_moment_A_m2=REFERENCE_MOMENT,
        min_B_at4R_T=float(fine.min()), max_B_at4R_T=float(fine.max()),
        min_B_at3R_T=float(at3.min()), angular_grid_relative_change=abs(float(low.min()/fine.min()-1)),
        loss_of_one_min_B_at4R_T=outage,
        loss_of_one_current_multiplier=float(fine.min()/min(outage)),
        route_height_range_m=[float(extents.min()), float(extents.max())],
        reference=summarize(loops), reduced_moment=summarize(rescale(loops, SMALL_MOMENT/REFERENCE_MOMENT)),
        pressure_sizing=[])
    for target_R, bmin in [(3., at3.min()), (4., fine.min())]:
        for pressure in [2e-9, 20e-9, 100e-9]:
            factor = np.sqrt(2*K.VACUUM_PERMEABILITY*pressure)/bmin
            s = summarize(rescale(loops, factor), include_mutual=False)
            row['pressure_sizing'].append(dict(target_radius_R=target_R, pressure_Pa=pressure,
                field_amplification_assumed=1., required_moment_A_m2=REFERENCE_MOMENT*factor,
                current_scale=factor, screened_mass_kg=s['screened_mass_kg'],
                conductor_kg=s['conductor_kg'], ideal_tensile_kg=s['ideal_tensile_kg'],
                refrigeration_W=s['refrigeration_W'], collector_area_m2=s['collector_area_m2'],
                self_energy_J=s['self_energy_J']))
    reserve_factor = np.sqrt(2*K.VACUUM_PERMEABILITY*100e-9)/min(outage)
    reserve = summarize(rescale(loops, reserve_factor), include_mutual=False)
    row['four_radius_100nPa_one_out_reserve'] = dict(all_installed_current_multiplier=reserve_factor,
        screened_installed_mass_kg=reserve['screened_mass_kg'],
        all_energized_refrigeration_W=reserve['refrigeration_W'],
        evidence='Every circuit installed for the sampled worst single outage; no transient or fault propagation model.')
    return row


def held_acceleration(samples, geo, distance, radius, orientation):
    phi = np.arange(64)*2*np.pi/64
    i, j = (1, 2) if orientation == 'axial' else (0, 1)
    q = distance*geo['frame'][:, 0, None, :]+radius*(np.cos(phi)[None, :, None]*geo['frame'][:, i, None, :]+np.sin(phi)[None, :, None]*geo['frame'][:, j, None, :])
    qa = distance*geo['frame_dd'][:, 0, None, :]+radius*(np.cos(phi)[None, :, None]*geo['frame_dd'][:, i, None, :]+np.sin(phi)[None, :, None]*geo['frame_dd'][:, j, None, :])
    g = relative_gravity(q, {k:v[:, None, :] for k, v in samples['positions'].items()}, samples['moon_a'][:, None, :])
    required = qa-g
    net = required.mean(axis=1)
    offset = q-distance*geo['frame'][:, 0, None, :]
    torque = np.cross(offset, required).mean(axis=1)
    equivalent = np.linalg.norm(net, axis=1)+np.linalg.norm(torque, axis=1)/radius
    return dict(net_vector_m_s2=net, norm_m_s2=np.linalg.norm(net, axis=1),
        torque_per_kg_N_m=torque, control_equivalent_m_s2=equivalent,
        distributed_norm_m_s2=np.linalg.norm(required, axis=2).mean(axis=1),
        max_differential_m_s2=float(np.max(np.linalg.norm(required-net[:, None], axis=2))),
        sunward_m_s2=np.sum(net*geo['frame'][:, 0], axis=1))


def held_budget(w, acceleration, obstacle_radius, pressure, specific_power=300.):
    """Conservative scalar mass/power closure including seven-day fuel buffer.

    Requires controlling a prescribed Sun-following path. The triangle bound
    treats wind and photons as opposing useful thrust at every epoch. Failure
    of this sufficient closure alone does not prove all trajectories impossible.
    """
    ve = 30000.; efficiency = .7; cant = np.cos(np.pi/4); margin = 1.25
    # Histories have uniform time spacing; endpoints receive half weight.
    average = lambda x: float(np.trapezoid(x)/(len(x)-1))
    w_per_n = ve/(2*efficiency*cant)
    gen = collector(1., specific_power=specific_power)
    # Array plus PMAD, separate 1 kW/kg thrusters and a 600 K thruster radiator.
    installed_per_W = margin*(gen['installed_mass_kg']+.001+radiator(1-efficiency, 600.)['mass_kg'])
    buffer_per_W = 7*K.JULIAN_DAY*2*efficiency/ve**2
    total_per_W = installed_per_W+buffer_per_W
    base_mass = w['screened_mass_without_generation_kg']
    cryo = w['refrigerator_W']
    photon_per_W = margin*gen['absorbed_photon_force_N']
    # Orientation-independent upper projected envelope; other structures remain.
    silhouette = 4*np.pi*w['radius_m']*(w['bundle_radius_m']+.25)
    fixed_photon = .9*K.SOLAR_CONSTANT*silhouette/K.SPEED_OF_LIGHT
    wind = 2*pressure*np.pi*obstacle_radius**2
    background_torque = w['moment_A_m2']*5e-9
    background_couple_force = background_torque/w['radius_m']
    a = acceleration['control_equivalent_m_s2']; peak = float(a.max())
    translation = acceleration['norm_m_s2']
    feedback = w_per_n*(total_per_W*peak+photon_per_W)
    numerator = cryo+w_per_n*(base_mass*peak+fixed_photon+wind+background_couple_force)
    out = dict(generator_specific_power_W_kg=specific_power, peak_net_acceleration_m_s2=float(translation.max()),
        peak_control_equivalent_m_s2=peak, mean_control_equivalent_m_s2=average(a),
        max_attitude_torque_per_mass_N_m_kg=float(np.linalg.norm(acceleration['torque_per_kg_N_m'], axis=1).max()),
        background_field_T=5e-9, background_torque_envelope_N_m=background_torque,
        mean_net_acceleration_m_s2=average(translation), minimum_sunward_acceleration_m_s2=float(acceleration['sunward_m_s2'].min()),
        distributed_to_net_mean_ratio=average(acceleration['distributed_norm_m_s2'])/average(translation),
        differential_stress_scale_J_kg=acceleration['max_differential_m_s2']*w['radius_m'],
        feedback_gain=feedback, sufficient_scalar_closure=feedback < 1,
        generator_only_required_specific_power_W_kg=margin*w_per_n*peak,
        perfect_reversal_wind_force_N=wind, envelope_absorbed_photon_force_N=fixed_photon,
        winding_shadow_area_upper_m2=silhouette, dry_winding_mass_kg=base_mass,
        evidence='Monthly DE440 finite-loop net translation and attitude torque with ideal rim couples, plus a 5 nT uniform-field torque envelope; conservative propulsion/power/fuel closure. No free orbit, validated wake or plume compatibility.')
    if feedback >= 1:
        out['failure'] = 'No finite positive solution to this installed-hardware scalar closure'
        return out
    peak_W = numerator/(1-feedback)
    total_mass = base_mass+total_per_W*peak_W
    force = total_mass*a+fixed_photon+wind+background_couple_force+photon_per_W*peak_W
    propulsion = force*w_per_n
    out.update(total_carried_mass_kg=total_mass, peak_total_bus_W=peak_W,
        mean_propulsion_W=average(propulsion), mean_total_bus_W=average(propulsion)+cryo,
        propellant_kg_s=average(propulsion)*2*efficiency/ve**2,
        seven_day_upper_fuel_buffer_kg=buffer_per_W*peak_W,
        installed_generation_W=margin*peak_W,
        collector_area_m2=margin*peak_W*gen['area_m2'],
        first_intercepted_collector_W=margin*peak_W*gen['incident_W'],
        conversion_loss_W=margin*peak_W*(gen['incident_W']-1.),
        mean_propulsion_over_23_835TW=average(propulsion)/23.835e12)
    return out


def main():
    before = time.monotonic(); signal.alarm(300)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    RAW.mkdir(parents=True, exist_ok=True)
    inputs = ['research/studies/protection_architecture/results/design_point.json',
        'research/studies/solar_shield_array/results/electromagnetic.json',
        'research/studies/solar_shield_array/results/electromagnetic_audit.json',
        'research/studies/solar_shield_array/scenario.json', 'protection/dynamics/ephemeris.json']
    parents = {p:json.loads((ROOT/p).read_text()) for p in inputs}
    em = parents[inputs[1]]
    out = dict(schema='terluna.research.primary-magnetic-architecture/1',
        producer=dict(source_hashes={p:digest(ROOT/p) for p in SOURCES}, constants=constants_used(SOURCES),
                      input_hashes={p:digest(ROOT/p) for p in inputs}, parent_commit='479816d288bdcecf6b1c36f0a8b4e2f328e3cd07'),
        evidence='Executed finite-loop vacuum fields, saved atmospheric loss and geometry reads, monthly prescribed-path DE440 force and engineering sensitivities. No plasma or new fleet trajectory.',
        accepted_magnetosphere=False, accepted_trajectory=False, accepted_return=False, accepted_cycle=False,
        original_target_W=23.835e12, target_changed=False, optical_target_radius_m=4*K.MOON_RADIUS,
        superconducting_planetary_ring_selected=False, atmosphere_cases=[], regional=[], upstream=[],
        upstream_margins=[], holding=[], tile_interactions=[], sensitivities=[], raw_inputs={},
        reading_rule='Read fields, loss screens, hardware assumptions and prescribed-path budgets at their individual scope; pressure equivalence grants no retention/dose or coverage.')
    for p in parents[inputs[0]]['points']:
        if p['shield'] != 'titania_stack' or p.get('level') != 'standard' or p.get('hole_factor') != 1.: continue
        series = p['moment_sweep']; idx = int(np.argmin(abs(np.array(series['moment_A_m2'])-SMALL_MOMENT)))
        out['atmosphere_cases'].append({k:p[k] for k in ['treatment', 'activity', 'exobase_radius_R', 'exobase_temperature_k', 'uv_driven_loss_kg_s', 'total_at_protected_radius_kg_s', 'moment_for_exosphere_loss_A_m2']} |
            dict(reduced_moment_A_m2=series['moment_A_m2'][idx],
                 reduced_moment_exosphere_loss_kg_s=series['loss_kg_s'][idx],
                 reduced_moment_total_loss_kg_s=series['loss_kg_s'][idx]+p['uv_driven_loss_kg_s'],
                 assumed_boundary_field_gain=2., wind_pressure_scope='Stored typical-wind case; not a storm calculation'))
    cases = [('historical_four_100', 100e3, 4, True), ('regional_four_250', 250e3, 4, False),
             ('regional_four_500', 500e3, 4, False), ('regional_four_750', 750e3, 4, False),
             ('regional_four_1000', 1000e3, 4, False), ('polar_two_500', 500e3, 2, False),
             ('polar_two_1000', 1000e3, 2, False)]
    for name, radius, count, historical in cases:
        out['regional'].append(regional_case(name, radius, count, historical))
    out['large_loop_benchmark'] = dict(radius_m=2e6, moment_A_m2=REFERENCE_MOMENT,
        hardware=winding(2e6, REFERENCE_MOMENT/(np.pi*(2e6)**2)),
        scope='Free-loop scaling only. Lunar-encircling deployment excluded by retained decision.')
    for radius in [100e3, 500e3, 1000e3, 2000e3, 4000e3]:
        for obstacle in [4*K.MOON_RADIUS, upstream_width(15e6, 4*K.MOON_RADIUS)['required_obstacle_radius_m']]:
            unit = obstacle_field(radius, obstacle)
            for pressure in [2e-9, 20e-9, 100e-9]:
                current = np.sqrt(2*K.VACUUM_PERMEABILITY*pressure)/unit
                w = winding(radius, current)
                out['upstream'].append(dict(radius_m=radius, obstacle_radius_m=obstacle,
                    pressure_Pa=pressure, hardware=w,
                    field_sizing_scope='Minimum finite-loop vacuum field over sphere, 361 polar samples; no downstream width guarantee',
                    angular_refinement_relative=abs(obstacle_field(radius, obstacle, 181)/unit-1)))
    for distance in [15e6, 78e6, 150e6, 1500e6]:
        for lateral in [0., 30e3, 60e3, 100e3]:
            for angle in [0., 5., 10.]:
                out['upstream_margins'].append(upstream_width(distance, 4*K.MOON_RADIUS, lateral_speed=lateral, angle_deg=angle))
    eph = Ephemeris(epoch=parents[inputs[3]]['epoch_tdb'])
    t = np.linspace(0., 30*K.JULIAN_DAY, 241)
    samples = eph.sample(t); geo = geometry(samples)
    out['ephemeris'] = dict(epoch_tdb=eph.epoch, kernel=eph.manifest, days=30., samples=len(t))
    eph.close()
    histories = {'t':t}
    for distance in [15e6, 78e6, 150e6, 1500e6]:
        for radius in [500e3, 2000e3, 4000e3]:
            for orientation in ['axial', 'transverse']:
                a = held_acceleration(samples, geo, distance, radius, orientation)
                key = f'd{int(distance/1000)}_r{int(radius/1000)}_{orientation}'
                histories[key] = a['net_vector_m_s2']
                histories[key+'_torque_per_kg'] = a['torque_per_kg_N_m']
                for width in ['ideal_zero_erosion', 'kinematic_60km_s_5deg']:
                    obstacle = 4*K.MOON_RADIUS if width.startswith('ideal') else upstream_width(distance, 4*K.MOON_RADIUS)['required_obstacle_radius_m']
                    current = np.sqrt(2*K.VACUUM_PERMEABILITY*2e-9)/obstacle_field(radius, obstacle)
                    w = winding(radius, current)
                    if not w['admissible']:
                        out['holding'].append(dict(distance_m=distance, radius_m=radius, orientation=orientation, width_case=width, hardware=w, failed=True))
                        continue
                    out['holding'].append(dict(distance_m=distance, radius_m=radius, orientation=orientation, width_case=width,
                        obstacle_radius_m=obstacle, hardware=w,
                        budgets=[held_budget(w, a, obstacle, 2e-9, sp) for sp in [300., 1000.]],
                        monthly_grid_peak_change=float(abs(a['control_equivalent_m_s2'][::2].max()/a['control_equivalent_m_s2'].max()-1))))
    np.savez_compressed(RAW/'holding.npz', **histories)
    out['raw_output'] = dict(path=str((RAW/'holding.npz').relative_to(ROOT)), sha256=digest(RAW/'holding.npz'))
    # Actual tile centres and attitudes, sampled at 0, 6, 12 and final-return hours.
    paths = ['research/runs/solar_shield_array/closure/closure_fine_0.npz',
             'research/runs/solar_shield_array/closure/closure_fine_1.npz',
             'research/runs/solar_shield_array/return_repair/repair_trial_local_refined_fine.npz']
    frames = []
    for path in paths:
        if digest(ROOT/path) != em['raw_inputs'][path]: raise ValueError('Changed saved trajectory: '+path)
        out['raw_inputs'][path] = digest(ROOT/path)
        z = np.load(ROOT/path)
        for idx in ([0, -1] if path == paths[0] else [-1]):
            frames.append((float(z['t'][idx]), z['state'][idx, :, :3].copy(), z['frames'][idx].copy()))
    for row in [out['regional'][0], out['regional'][2], out['regional'][4]]:
        for moment in [SMALL_MOMENT, REFERENCE_MOMENT]:
            records = []
            for ti, q, frame in frames:
                for phase in np.arange(4)*np.pi/2:
                    # Explicit orientation scenario; no lunar body ephemeris is claimed.
                    rot = Rotation.from_euler('x', K.EARTH_OBLIQUITY_DEG, degrees=True).as_matrix() @ Rotation.from_rotvec(np.array([0., 0., phase+2*np.pi*ti/(K.SIDEREAL_MONTH_DAYS*K.JULIAN_DAY)])).as_matrix()
                    loops = rescale(row['loops'], moment/REFERENCE_MOMENT)
                    oriented = [{**l, 'centre':rot@l['centre'], 'normal':rot@l['normal']} for l in loops]
                    b = field(q, oriented)
                    torque = np.linalg.norm(np.cross(frame[:, 2]*1e13, b), axis=1)
                    delta = 100.; force = np.column_stack([((field(q+np.eye(3)[j]*delta, oriented)-field(q-np.eye(3)[j]*delta, oriented)) @ (frame[:, 2]*1e13))/(2*delta) for j in range(3)])
                    records.append(dict(time_s=ti, lunar_spin_phase_rad=phase, max_tile_B_T=float(np.linalg.norm(b, axis=1).max()),
                        max_tile_torque_N_m=float(torque.max()), max_tile_force_N=float(np.linalg.norm(force, axis=1).max())))
            out['tile_interactions'].append(dict(architecture=row['name'], moment_A_m2=moment, samples=records,
                orientation_scope='Spin axis set to ecliptic north; four initial lunar longitudes; true lunar pole/libration and plasma currents omitted. Actual centres and tile attitudes retained.'))
    for radius in [100e3, 500e3, 1000e3]:
        loops = regional_loops(radius, REFERENCE_MOMENT, 4, historical=radius==100e3)
        for limit in [2.4, 10., 20.]:
            for leak, cop in [(.01, .01), (.05, .005), (.1, .002)]:
                h = [winding(radius, l['current'], field_limit=limit, heat_leak=leak, cop=cop) for l in loops]
                out['sensitivities'].append(dict(radius_m=radius, field_limit_T=limit, heat_leak_W_m2=leak, cop=cop,
                    screened_mass_kg=sum(w['screened_mass_with_generation_kg'] for w in h),
                    refrigeration_W=sum(w['refrigerator_W'] for w in h)))
    out['joint_sensitivity'] = []
    for spacing in [1000., 10000., 100000.]:
        for resistance in [.5e-9, 1e-9, 2e-9]:
            h = winding(500e3, out['regional'][2]['loops'][0]['current'], joint_spacing=spacing, joint_ohm=resistance)
            out['joint_sensitivity'].append(dict(spacing_m=spacing, resistance_ohm=resistance,
                cold_joint_W_four_loops=4*h['joints_cold_W'], refrigeration_W_four_loops=4*h['refrigerator_W']))
    out['material_sensitivity'] = []
    for strength in [1e5, 3e5, 1e6]:
        for density in [3e7, 1e8]:
            for row in [out['regional'][0], out['regional'][2], out['regional'][4]]:
                hw = [winding(l['radius'], l['current'], specific_strength=strength, current_density=density) for l in row['loops']]
                out['material_sensitivity'].append(dict(architecture=row['name'], specific_strength_J_kg=strength,
                    current_density_A_m2=density, screened_mass_kg=sum(h['screened_mass_with_generation_kg'] for h in hw)))
    # Independent parametric plasma-current ring. Guide infrastructure is unknown.
    me = K.VACUUM_PERMEABILITY*K.ELEMENTARY_CHARGE**2/(4*np.pi*K.CLASSICAL_ELECTRON_RADIUS)
    radius = 1e7; current = REFERENCE_MOMENT/(np.pi*radius**2); minor = 1e5
    inductance = K.VACUUM_PERMEABILITY*radius*(np.log(8*radius/minor)-2)
    magnetic_energy = .5*inductance*current**2
    carriers = []
    for electron_eV in [100., 1000., 100000.]:
        energy = electron_eV*K.ELEMENTARY_CHARGE; gamma = 1+energy/(me*K.SPEED_OF_LIGHT**2)
        v = K.SPEED_OF_LIGHT*np.sqrt(1-gamma**-2)
        number = current*2*np.pi*radius/(K.ELEMENTARY_CHARGE*v)
        carriers.append(dict(electron_energy_eV=electron_eV, drift_m_s=v, count=number,
            electron_kinetic_energy_J=number*energy, neutralizing_proton_mass_kg=number*K.PROTON_MASS,
            density_m3=number/(2*np.pi*radius*np.pi*minor**2)))
    out['plasma_current'] = dict(radius_m=radius, minor_radius_m=minor, current_A=current,
        moment_A_m2=REFERENCE_MOMENT, magnetic_energy_J=magnetic_energy, carriers=carriers,
        thin_current_path_hoop_tension_N=K.VACUUM_PERMEABILITY*current**2/(4*np.pi)*(np.log(8*radius/minor)-1),
        field_loss_power=[dict(energy_decay_s=tau, sustain_W=magnetic_energy/tau,
                              sustain_at_10percent_drive_W=10*magnetic_energy/tau) for tau in [1., 1000., K.JULIAN_DAY, 30*K.JULIAN_DAY]],
        accepted=False, evidence='Exploratory current inventory and energy-decay sensitivity. No equilibrium, guide-field mass, current closure, efficiency, stability or loss time established. It is not a superconducting planetary ring.')
    out['plasma_scales'] = [plasma(), plasma(20., 800., 20., 15.), plasma(100., 800., 30., 30.)]
    out['reduced_moment_pressure_sensitivity'] = [dict(pressure_Pa=p, boundary_field_gain=gain,
        central_dipole_standoff_R=(gain*K.VACUUM_PERMEABILITY/(4*np.pi)*SMALL_MOMENT/np.sqrt(2*K.VACUUM_PERMEABILITY*p))**(1/3)/K.MOON_RADIUS)
        for p in [out['plasma_scales'][0]['ram_pressure_Pa'], 2e-9, 20e-9, 100e-9] for gain in [1., 2.]]
    out['resources'] = dict(elapsed_s=time.monotonic()-before, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        budget_wall_s=900, processes=1, numerical_threads=1, address_space_GiB=2,
        raw_output_MiB=sum(p.stat().st_size for p in RAW.iterdir() if p.is_file())/1024**2)
    out['completed'] = True
    (HERE/'results/magnetic_architecture.json').write_text(json.dumps(out, indent=2, default=plain)+'\n')
    print(json.dumps(out['resources']), flush=True)


if __name__ == '__main__':
    main()
