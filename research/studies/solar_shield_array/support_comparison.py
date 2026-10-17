"""Material substitutions, finite sail supports and spin necessary conditions.

The preserved tile design supplies the whole-spacecraft ledger. New support
architectures are component screens, not replacements for that ledger.
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
from pathlib import Path
import resource
import time

import numpy as np

from engineering.array_tile import Material, TileDesign, membrane_coefficient, perimeter_truss
from engineering.array_support import (SupportInfeasible, scalloped_support,
                                       normal_acceleration_coefficients, box_normal_acceleration_envelope)
from engineering.array_spin import ballasted_disk_screen, disk_grid_coverage
from shared import constants as K
from shared.provenance import constants_used
from . import tile_design as td

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/support_comparison.json'
FILES = sorted(set(td.FILES + [
    'research/studies/solar_shield_array/support_comparison.py',
    'engineering/array_support.py', 'engineering/array_spin.py']))
INPUTS = {
    'baseline': ('research/studies/solar_shield_array/results/tile_design.json', 'terluna.research.array-tile-design/1'),
    'materials': ('research/studies/solar_shield_array/support_materials.json', 'terluna.research.array-support-materials/1'),
    'sources': ('research/studies/solar_shield_array/support_sources.json', 'terluna.research.array-support-sources/1'),
}


def material(record):
    return Material(**{k: record[k] for k in ('density_kg_m3', 'youngs_modulus_Pa', 'admissible_stress_Pa')})


def inherited_design(baseline, kind='annulus', side=10000.):
    key = f'{kind}_reference' if side == 10000 else f'{kind}_side_{side:g}m'
    values = dict(baseline['hardware']['variants'][key]['design'])
    values['structure'] = Material(**values['structure'])
    values['required_momentum_N_m_s'] = baseline['hardware']['inherited_momentum_allocation_N_m_s']*(side/10000)**3
    values['required_torque_N_m'] = 1e6*(side/10000)**3
    return TileDesign(**values)


def coverage_inventory(side_m, retreat_m, population):
    """Sufficient planar rectangle containment, not a fleet ray trace.

    Place a uniformly scaled 8900 x 9890 m active rectangle inside the loaded
    sail's central square. Scale both old pitches by the same factor. Shapes
    are prescribed after prestress; manufacture must compensate shortening.
    """
    core = side_m-2*retreat_m
    scale = min(side_m/10000, core/8900, core/9890)
    if not 0 < scale <= 1:
        raise ValueError('Expected positive geometry no larger than the reference')
    multiplier = 1/scale**2
    return dict(guaranteed_loaded_core_side_m=core, pitch_scale=scale,
                pitch_x_m=8500*scale, pitch_y_m=9286.456*scale,
                inventory_multiplier=multiplier, equivalent_population=population*multiplier,
                reference_gross_area_per_unit_m2=1e8/multiplier,
                rule='sufficient local planar active-rectangle containment; preserves old overlap ratios; projection, scallop seams, shadow forces and full fleet coverage require replay')


def material_substitutions(baseline, catalogue, radius):
    rows = []
    for kind in ('annulus', 'window'):
        base = inherited_design(baseline, kind)
        for entry in catalogue['compression_candidates']:
            budget = td.size_momentum_hardware(replace(base, structure=material(entry)), radius)
            rows.append(dict(kind=kind, material_id=entry['id'],
                             specific_stiffness_relative_to_baseline=(entry['youngs_modulus_Pa']/entry['density_kg_m3'])/(100e9/1600),
                             budget=budget))
    return rows


def support_case(design, compression, cable, retreat_fraction, population, *,
                 tension, pressure=0., acceleration=0., depth=None, case_id='',
                 equipment_tip_force=0., equipment_tip_moment=0.):
    side = design.side_m
    depth = depth or side/100
    coverage = coverage_inventory(side, side*retreat_fraction, population)
    case = dict(id=case_id, side_m=side, edge_retreat_fraction=retreat_fraction,
                boom_depth_m=depth, compression_material_id=compression['id'],
                cable_material_id=cable['id'], tension_N_m=tension,
                pressure_for_support_Pa=pressure, structural_acceleration_m_s2=acceleration,
                equipment_tip_force_N=equipment_tip_force, equipment_tip_moment_N_m=equipment_tip_moment,
                coverage=coverage, complete_spacecraft_budget=False, validated_design=False)
    try:
        support = scalloped_support(side, tension, side*retreat_fraction,
                    material(compression), material(cable), boom_depth_m=depth,
                    bay_m=side/200, chord_radius_m=side/5000, hub_side_m=side/10,
                    allowable_anchor_retreat_m=side/2000,
                    allowable_lateral_deflection_m=side/2000,
                    membrane_normal_pressure_Pa=pressure,
                    structural_acceleration_m_s2=acceleration,
                    tip_transverse_force_N=equipment_tip_force, tip_moment_N_m=equipment_tip_moment)
    except SupportInfeasible as exc:
        return dict(case, converged=False, failure=str(exc))
    area = support['geometry']['retained_area_m2']
    sheet_sigma = design.film_areal_mass_kg_m2+design.substrate_areal_mass_kg_m2+design.trim_area_fraction*design.trim_hardware_areal_mass_kg_m2
    sheet_mass = area*sheet_sigma
    subtotal = sheet_mass+support['mass_kg']
    normalization = coverage['reference_gross_area_per_unit_m2']
    normalized = subtotal/normalization
    # Whole-sheet electrical/control hardware is deliberately not borrowed
    # from a differently shaped spacecraft. Expose this missing closure.
    inertia = np.asarray(support['inertia_kg_m2'])+sheet_mass*np.asarray(support['geometry']['membrane_area_inertia_per_mass_m2'])
    base_population = coverage['equivalent_population']
    required_gap = 150.
    gap = coverage['pitch_x_m']*math.sin(math.radians(1.5))-2*support['half_depth_m']
    stress = tension/(design.substrate_areal_mass_kg_m2/design.substrate_density_kg_m3)
    return dict(case, converged=True, support=support, retained_sheet_mass_kg=sheet_mass,
                support_and_sheet_subtotal_kg=subtotal,
                subtotal_kg_per_reference_gross_m2=normalized,
                subtotal_ratio_to_inherited_total_allowance=normalized/design.target_total_areal_mass_kg_m2,
                support_and_sheet_locked_inertia_kg_m2=inertia.tolist(),
                substrate_stress_Pa=stress,
                substrate_stress_utilization=stress/design.substrate_admissible_stress_Pa,
                nominal_parallel_support_clearance_m=gap,
                inherited_150m_support_gate_passes=gap >= required_gap,
                fixed_avionics_inventory_kg=base_population*design.avionics_and_sensors_kg,
                fixed_housekeeping_inventory_W=base_population*design.housekeeping_W,
                assumed_annual_replacement_jobs=base_population*(1/15+.01),
                omitted=['rotors and motor supports', 'power generation, storage and wiring',
                         'avionics and mass-growth reserve', 'closed equipment-force and flexible attitude iteration',
                         'cable normal shape, optical strain and global spacecraft equilibrium'],
                interpretation='conditional component subtotal under prescribed loads; never an accepted total mass, optical performance or trajectory')


def numpy_scalar(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Unsupported JSON value: {type(value)}")


def run():
    start = time.perf_counter()
    data = {}
    for key, (name, schema) in INPUTS.items():
        data[key] = json.loads((ROOT/name).read_text())
        if data[key]['schema'] != schema:
            raise ValueError(f'Unexpected schema: {name}')
    baseline, catalogue = data['baseline'], data['materials']
    comp = {x['id']: x for x in catalogue['compression_candidates']}
    cable = {x['id']: x for x in catalogue['tension_candidates']}
    radius = baseline['hardware']['variants']['annulus_reference']['assembled_momentum_screen']['radius_m']
    population = baseline['fleet_service']['population']
    d = inherited_design(baseline)
    ref = baseline['hardware']['variants']['annulus_reference']
    t, states, frames, omega, alpha = td.ideal_attitude(radius, 2881)
    radial = np.einsum('tji,tj->ti', frames, states[:, :3]/radius)
    coeff = normal_acceleration_coefficients(omega, alpha, radial, K.MOON_GM/radius**3)
    normal = float(d.side_m/2*np.max(np.sum(np.abs(coeff), axis=1)))
    old_normal = (d.angular_acceleration_rad_s2+d.angular_speed_rad_s**2)*d.side_m/math.sqrt(2)
    refinement = []
    for samples in (1441, 2881, 5761):
        tt, ss, ff, ww, aa = td.ideal_attitude(radius, samples)
        rr = np.einsum('tji,tj->ti', ff, ss[:, :3]/radius)
        cc = normal_acceleration_coefficients(ww, aa, rr, K.MOON_GM/radius**3)
        refinement.append(dict(samples=samples, maximum_normal_acceleration_m_s2=float(d.side_m/2*np.max(np.sum(np.abs(cc), axis=1)))))
    tension = ref['membrane']['tension_N_m']
    fixed = []
    for fraction in (.005, .05, .1, .2):
        fixed.append(support_case(d, comp['generic_composite_baseline'], cable['generic_cable_reference'], fraction,
                                 population, tension=tension, case_id=f'geometry_{fraction:g}'))
    for entry in catalogue['tension_candidates']:
        fixed.append(support_case(d, comp['m55j_2500_ud'], entry, .05, population,
                                 tension=tension, case_id=f'cable_{entry["id"]}'))
    loaded = []
    for side in (1000., 3000., 10000.):
        ds = inherited_design(baseline, side=side)
        accel = normal*side/d.side_m
        sigma_local = ds.film_areal_mass_kg_m2+ds.substrate_areal_mass_kg_m2+ds.trim_hardware_areal_mass_kg_m2
        q = ds.pressure_Pa+sigma_local*accel
        load_tension = membrane_coefficient()*q*side**2*ds.membrane_dynamic_factor/ds.allowable_sag_m
        for fraction in (.005, .05, .1):
            for depth_ratio in (.005, .01):
                loaded.append(support_case(ds, comp['m55j_2500_ud'], cable['kevlar49_yarn'], fraction,
                    population, tension=load_tension, pressure=ds.membrane_dynamic_factor*q,
                    acceleration=ds.membrane_dynamic_factor*box_normal_acceleration_envelope(
                        [side/2+math.sqrt(2)*side*depth_ratio/2]*2+[side*depth_ratio/2+ds.chord_radius_m+ds.allowable_sag_m],
                        omega, alpha, radial, K.MOON_GM/radius**3), depth=side*depth_ratio,
                    case_id=f'loaded_{side:g}_{fraction:g}_{depth_ratio:g}'))
    # Explicit equipment-reaction sensitivity at fixed assumed loads, never
    # pass it off as an iterated complete spacecraft or as hardware mass.
    equipment = []
    for force, moment in ((0., 0.), (100., 1e6), (1000., 1e7)):
        q = d.pressure_Pa+(d.film_areal_mass_kg_m2+d.substrate_areal_mass_kg_m2+d.trim_hardware_areal_mass_kg_m2)*normal
        equipment.append(support_case(d, comp['m55j_2500_ud'], cable['kevlar49_yarn'], .05,
            population, tension=membrane_coefficient()*q*d.side_m**2*2/d.allowable_sag_m,
            pressure=2*q, acceleration=2*box_normal_acceleration_envelope(
                [d.side_m/2+math.sqrt(2)*d.frame_depth_m/2]*2+[d.frame_depth_m/2+d.chord_radius_m+d.allowable_sag_m],
                omega, alpha, radial, K.MOON_GM/radius**3), case_id=f'equipment_{force:g}_{moment:g}',
            equipment_tip_force=force, equipment_tip_moment=moment))
    # Same complete perimeter ledger, new explicit normal-pressure envelope.
    refined = td.size_momentum_hardware(replace(d, angular_speed_rad_s=0.,
                      angular_acceleration_rad_s2=float(normal*math.sqrt(2)/d.side_m)), radius)
    refined_high_modulus = td.size_momentum_hardware(replace(d, structure=material(comp['m55j_2500_ud']),
                      angular_speed_rad_s=0., angular_acceleration_rad_s2=float(normal*math.sqrt(2)/d.side_m)), radius)
    refined_high_modulus['equivalent_load_input_rule'] = 'same resolved normal-pressure adapter as refined_load_perimeter; M55J composite proxy replaces baseline material'
    refined['equivalent_load_input_rule'] = 'alpha input is an adapter to the old pressure ledger, not the physical angular acceleration; original commanded attitude still drives momentum sizing'
    spin = []
    disk_grid = disk_grid_coverage(5000., 8500., 9286.456)
    for law in ('fixed_normal_pressure_benchmark', 'rolled_command'):
        for speed in (.0001, .0003, .001, .003, .01):
            fixed_law = law.startswith('fixed')
            result = ballasted_disk_screen(radius_m=5000.,
                areal_mass_kg_m2=d.film_areal_mass_kg_m2+d.substrate_areal_mass_kg_m2+d.trim_area_fraction*d.trim_hardware_areal_mass_kg_m2,
                spin_rad_s=speed, solar_pressure_Pa=d.pressure_Pa, allowable_sag_m=5.,
                load_bearing_thickness_m=d.substrate_areal_mass_kg_m2/d.substrate_density_kg_m3,
                admissible_stress_Pa=d.substrate_admissible_stress_Pa,
                times_s=t, omega_rad_s=np.zeros_like(omega) if fixed_law else omega,
                alpha_rad_s2=np.zeros_like(alpha) if fixed_law else alpha,
                radial_body=radial, mean_motion_squared_s2=0. if fixed_law else K.MOON_GM/radius**3)
            photon_bound = 2*K.SOLAR_CONSTANT/K.SPEED_OF_LIGHT*result['film_area_m2']*result['radius_m']
            result['direct_solar_photon_torque_upper_bound_N_m'] = photon_bound
            result['rigid_required_torque_over_direct_photon_bound'] = result['peak_full_rigid_nongravitational_torque_N_m']/photon_bound
            result['photon_bound_rule'] = 'all incident photons transfer twice their momentum at the full radius; ignores incidence and filtering losses; no assisted torque or internal counter-momentum hardware'
            multiplier = disk_grid['count_multiplier_at_fixed_protected_area']
            spin.append(dict(law=law, spin_rad_s=speed, result=result,
                subtotal_kg_per_reference_gross_m2=result['ideal_film_plus_ballast_mass_kg']*multiplier/1e8,
                equivalent_population=population*multiplier, complete_spacecraft_budget=False))
    records = {name: dict(schema=schema, sha256=td.digest(ROOT/name)) for name, schema in INPUTS.values()}
    product = dict(schema='terluna.research.array-support-comparison/1',
        producer=dict(files={name: td.digest(ROOT/name) for name in FILES}, constants=constants_used(FILES)),
        inputs=records, baseline_commit='e1196f5',
        evidence='Executed conditional material/component calculations, inherited whole-tile ledger substitutions and source review; no accepted new architecture or qualified material life.',
        reading_rule='Compare like load and area denominators. Only material_substitutions and the refined_load perimeter cases retain the prior whole-tile ledger. Support cases are incomplete subtotals. Spin cases are necessary-condition diagnostics. Nonconvergence and failed geometric screens remain failures.',
        material_substitutions=material_substitutions(baseline, catalogue, radius),
        loading=dict(original_norm_acceleration_m_s2=old_normal,
                     resolved_normal_acceleration_m_s2=float(normal), refinement=refinement,
                     assumptions='flat rigid sheet, ideal circular command, lunar linear gravity gradient; full unrelieved solar pressure; factor-two load screen is not a general resonance bound'),
        refined_load_perimeter=refined, refined_load_high_modulus_perimeter=refined_high_modulus,
        fixed_prestress=dict(tension_N_m=tension, perimeter=perimeter_truss(10000., tension), cases=fixed,
                            scope='in-plane prestress; no normal or equipment loads; all structures unqualified and member constraints differ from inherited perimeter'),
        command_loaded_supports=loaded, equipment_load_sensitivity=equipment,
        spinning_membranes=dict(cases=spin, planar_grid=disk_grid,
                               scope='uniform disk, independently ballasted edge, no edge hoop tension; no complete actuator, attachment or normal edge dynamics'),
        runtime_s=time.perf_counter()-start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    OUT.write_text(json.dumps(product, indent=2, allow_nan=False, default=numpy_scalar)+'\n')
    print(json.dumps(dict(output=str(OUT), runtime_s=product['runtime_s'],
                         material_cases=len(product['material_substitutions']), support_cases=len(fixed)+len(loaded)+len(equipment), spin_cases=len(spin))))
    return product


if __name__ == '__main__':
    run()
