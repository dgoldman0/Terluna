"""Bounded algebraic geometry screen; no coupled attitude or trajectory search.

Run with one BLAS thread. Every design and failure goes to ignored JSONL;
only distributions, strata and representative records enter the small product.
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
import math
from pathlib import Path
import resource
import time

import numpy as np

from engineering.array_support import (SupportInfeasible, scalloped_support,
    normal_acceleration_coefficients, box_normal_acceleration_envelope)
from engineering.array_tile import membrane_coefficient
from shared import constants as K
from shared.provenance import constants_used
from . import support_comparison as previous
from . import tile_design as td

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/geometry_screen.json'
ARCHIVE = ROOT/'research/runs/solar_shield_array_geometry_screen'
FILES = sorted(set(previous.FILES+[
    'research/studies/solar_shield_array/geometry_screen.py']))
INPUTS = dict(previous.INPUTS, comparison=(
    'research/studies/solar_shield_array/results/support_comparison.json',
    'terluna.research.array-support-comparison/1'))
GRID = dict(side_m=[1000., 2000., 3000., 5000., 7500., 10000.],
    retreat_fraction=[.005, .025, .05, .075, .1],
    depth_fraction=[.003, .005, .0075, .01, .015],
    compression_material_id=['generic_composite_baseline', 'as4_8552_ud', 'm55j_2500_ud'],
    cable_material_id=['kevlar49_yarn', 'spectra1000_100'],
    normal_load_factor=[.25, .5, 1.])
SHINGLES = [1.5, 2., 3.]
PITCH_FACTORS = [.9, 1.]
EVIDENCE = ('Conditional algebraic sizing and sufficient local planar geometry; '
    'no achieved lower-load attitude, complete spacecraft, validated material life, '
    'coupled flexible control, full-fleet coverage or collision certificate.')


def mechanical_configurations(grid=None):
    """Stable enumeration; spacing never triggers a new mechanical sizing call."""
    grid = GRID if grid is None else grid
    names = list(GRID)
    if set(grid) != set(names) or any(not grid[name] for name in names):
        raise ValueError('Expected exactly the nonempty mechanical grid axes')
    for index, values in enumerate(product(*(grid[name] for name in names))):
        yield dict(id=f'm{index:05d}', **dict(zip(names, values)))


def load_inputs():
    values = {}
    for key, (path, schema) in INPUTS.items():
        values[key] = json.loads((ROOT/path).read_text())
        if values[key].get('schema') != schema:
            raise ValueError(f'Unexpected input schema: {path}')
    return values


def context_from_inputs(data):
    baseline, comparison = data['baseline'], data['comparison']
    ref = baseline['hardware']['variants']['annulus_reference']
    design = ref['design']
    lattice = baseline['fleet_service']['ideal_service_lattice']
    safety_values = {row['config']['safety_clearance_m']
                     for row in baseline['local_failures']['cases'].values()}
    if len(safety_values) != 1:
        raise ValueError('The inherited clearance gate is ambiguous')
    radius = ref['assembled_momentum_screen']['radius_m']
    # Use the same sampled command as the previous normal-load product. The
    # new load factors scale its acceleration coefficients, not a real attitude.
    samples = next(row['samples'] for row in comparison['loading']['refinement']
                   if math.isclose(row['maximum_normal_acceleration_m_s2'],
                       comparison['loading']['resolved_normal_acceleration_m_s2'], rel_tol=1e-12))
    _, states, frames, omega, alpha = td.ideal_attitude(radius, samples)
    radial = np.einsum('tji,tj->ti', frames, states[:, :3]/radius)
    n2 = K.MOON_GM/radius**3
    coeff = normal_acceleration_coefficients(omega, alpha, radial, n2)
    reproduced = float(design['side_m']/2*np.max(np.sum(np.abs(coeff), axis=1)))
    normal = comparison['loading']['resolved_normal_acceleration_m_s2']
    if not math.isclose(reproduced, normal, rel_tol=1e-12):
        raise ValueError('Inherited command and normal-load product disagree')
    old_support = next(row['support'] for row in comparison['command_loaded_supports']
        if row['side_m'] == design['side_m'] and row['boom_depth_m'] == design['frame_depth_m'])
    context = dict(design=design, lattice=lattice,
        population=baseline['fleet_service']['population'],
        safety_clearance_m=safety_values.pop(), reference_side_m=design['side_m'],
        normal_acceleration_m_s2=normal, samples=samples, omega=omega,
        alpha=alpha, radial=radial, mean_motion_squared_s2=n2,
        hub_fraction=old_support['hub_member']['length_m']/design['side_m'],
        compression={row['id']: row for row in data['materials']['compression_candidates']},
        cables={row['id']: row for row in data['materials']['tension_candidates']})
    context['inherited_values'] = dict(
        reference_side_m=design['side_m'], reference_gross_area_m2=design['side_m']**2,
        reference_population=context['population'], service_lattice=lattice,
        safety_clearance_m=context['safety_clearance_m'], design=design,
        normal_acceleration_m_s2=normal, command_samples=samples,
        hub_fraction=context['hub_fraction'])
    return context


def mechanical_case(config, context):
    """One support subtotal, including its finite-envelope acceleration bound."""
    row = dict(config, converged=False, complete_spacecraft_budget=False,
               validated_design=False, lower_load_attitude_demonstrated=False)
    d = context['design']
    side = config['side_m']
    scale = side/context['reference_side_m']
    depth = side*config['depth_fraction']
    chord_radius = d['chord_radius_m']*scale
    factor = config['normal_load_factor']
    normal = context['normal_acceleration_m_s2']*scale*factor
    # Full trim-hardware density over the sheet is the inherited local-load
    # envelope; only its actual area fraction contributes to integrated mass.
    local_sigma = (d['film_areal_mass_kg_m2']+d['substrate_areal_mass_kg_m2']+
                   d['trim_hardware_areal_mass_kg_m2'])
    sigma = (d['film_areal_mass_kg_m2']+d['substrate_areal_mass_kg_m2']+
             d['trim_area_fraction']*d['trim_hardware_areal_mass_kg_m2'])
    pressure = d['membrane_dynamic_factor']*(d['pressure_Pa']+local_sigma*normal)
    tension = membrane_coefficient()*pressure*side**2/d['allowable_sag_m']
    bend_allowance = d['frame_allowable_edge_bow_m']*scale
    deformation_allowance = d['allowable_sag_m']+bend_allowance
    half_spans = np.array([side/2+depth/(2*math.sqrt(2))+chord_radius]*2+
                         [depth/2+chord_radius])+deformation_allowance
    def acceleration(spans):
        return factor*d['membrane_dynamic_factor']*box_normal_acceleration_envelope(
            spans, context['omega'], context['alpha'], context['radial'],
            context['mean_motion_squared_s2'])
    assumed_acceleration = acceleration(half_spans)
    try:
        # This is algebraic envelope/self-load closure, not a coupled solve.
        for iteration in range(8):
            support = scalloped_support(side, tension, side*config['retreat_fraction'],
                previous.material(context['compression'][config['compression_material_id']]),
                previous.material(context['cables'][config['cable_material_id']]),
                boom_depth_m=depth, bay_m=d['frame_bay_m']*scale,
                chord_radius_m=chord_radius, hub_side_m=context['hub_fraction']*side,
                allowable_anchor_retreat_m=bend_allowance,
                allowable_lateral_deflection_m=bend_allowance,
                membrane_normal_pressure_Pa=pressure,
                structural_acceleration_m_s2=assumed_acceleration)
            required_spans = np.array([support['outside_span_m']/2]*2+
                [support['half_depth_m']])+deformation_allowance
            required_acceleration = acceleration(required_spans)
            if required_acceleration <= assumed_acceleration*(1+1e-9):
                break
            half_spans = np.maximum(half_spans, required_spans)*(1+1e-6)
            assumed_acceleration = acceleration(half_spans)
        else:
            raise SupportInfeasible('Finite material-envelope self-load did not close')
    except SupportInfeasible as exc:
        return dict(row, failure=str(exc), tension_N_m=tension,
                    membrane_normal_pressure_Pa=pressure)
    area = support['geometry']['retained_area_m2']
    sheet_mass = sigma*area
    stress = tension/(d['substrate_areal_mass_kg_m2']/d['substrate_density_kg_m3'])
    return dict(row, converged=True, retained_area_m2=area,
        retained_sheet_mass_kg=sheet_mass, support_mass_kg=support['mass_kg'],
        components_kg=support['components_kg'],
        support_and_sheet_subtotal_kg=sheet_mass+support['mass_kg'],
        tension_N_m=tension, membrane_normal_pressure_Pa=pressure,
        normal_membrane_acceleration_m_s2=normal,
        structural_acceleration_m_s2=assumed_acceleration,
        required_finite_envelope_acceleration_m_s2=required_acceleration,
        acceleration_envelope_half_spans_m=half_spans.tolist(),
        finite_envelope_iterations=iteration+1,
        support_outside_span_m=support['outside_span_m'],
        support_half_depth_m=support['half_depth_m'],
        deformation_allowance_m=deformation_allowance,
        protected_core_side_m=side-2*side*config['retreat_fraction'],
        support_shadow_upper_bound_m2=support['projected_support_obscuration_upper_bound_m2'],
        support_shadow_fraction_upper_bound=support['projected_support_obscuration_upper_bound_m2']/area,
        substrate_stress_utilization=stress/d['substrate_admissible_stress_Pa'],
        thin_wall_screen_passes=support['thin_wall_screen_passes'],
        membrane_normal_force_N=support['membrane_normal_resultant_N'],
        membrane_normal_force_per_boom_N=support['membrane_normal_force_per_boom_N'],
        cable_axial_strain=support['cable_axial_strain'])


def spacing_case(mechanical, context, shingle_deg, pitch_factor, index):
    """Pose-projected inscribed rectangle and parallel x-neighbour envelope.

    Projection is about the local y axis. The shingle angle affects only this
    service-geometry scenario: mechanical loads still use the inherited 1.5°
    command or an explicit scaled-load hypothetical. No new command is solved.
    """
    if not 0 <= shingle_deg < 90 or not 0 < pitch_factor <= 1:
        raise ValueError('Expected an acute service angle and pitch factor in (0,1]')
    row = dict(id=f"{mechanical['id']}_s{index}", mechanical_id=mechanical['id'],
        side_m=mechanical['side_m'], retreat_fraction=mechanical['retreat_fraction'],
        depth_fraction=mechanical['depth_fraction'],
        compression_material_id=mechanical['compression_material_id'],
        cable_material_id=mechanical['cable_material_id'],
        normal_load_factor=mechanical['normal_load_factor'],
        shingle_deg=shingle_deg, pitch_factor=pitch_factor,
        converged=mechanical['converged'], complete_spacecraft_budget=False,
        validated_design=False, lower_load_attitude_demonstrated=False)
    if not mechanical['converged']:
        return dict(row, failure=mechanical['failure'], filters={'mechanical_converged': False})
    d, lattice = context['design'], context['lattice']
    cosine = math.cos(math.radians(shingle_deg))
    core = mechanical['protected_core_side_m']
    scale = min(mechanical['side_m']/context['reference_side_m'],
        core*cosine/lattice['optical_width_m'], core/lattice['optical_height_m'])*pitch_factor
    pitches = [pitch*scale for pitch in lattice['pitch_m']]
    projected_width = lattice['optical_width_m']*scale
    projected_height = lattice['optical_height_m']*scale
    multiplier = 1/scale**2
    normalization = context['reference_side_m']**2/multiplier
    population = context['population']*multiplier
    subtotal = mechanical['support_and_sheet_subtotal_kg']/normalization
    half_depth = mechanical['support_half_depth_m']+mechanical['deformation_allowance_m']
    gap = pitches[0]*math.sin(math.radians(shingle_deg))-2*half_depth
    member_ok = mechanical['thin_wall_screen_passes'] and mechanical['substrate_stress_utilization'] <= 1
    bus = {}
    for name, fraction in [('fixed_inherited_equipment', 1.),
                           ('area_scaled_equipment_hypothesis', (mechanical['side_m']/context['reference_side_m'])**2)]:
        mass = d['avionics_and_sensors_kg']*fraction
        power = d['housekeeping_W']*fraction
        bus[name] = dict(assumed_mass_per_unit_kg=mass, assumed_housekeeping_per_unit_W=power,
            subtotal_with_bus_kg_per_reference_gross_m2=subtotal+mass/normalization,
            fleet_bus_inventory_kg=mass*population,
            fleet_housekeeping_W=power*population,
            bus_requirements_demonstrated=False, power_system_mass_included=False)
    return dict(row, pitch_scale=scale, pitch_m=pitches,
        projected_protected_rectangle_m=[projected_width, projected_height],
        physical_protected_rectangle_m=[projected_width/cosine, projected_height],
        protected_core_side_m=core, reference_gross_area_per_unit_m2=normalization,
        equivalent_population=population, inventory_multiplier=multiplier,
        support_and_sheet_subtotal_kg=mechanical['support_and_sheet_subtotal_kg'],
        subtotal_kg_per_reference_gross_m2=subtotal,
        nominal_parallel_envelope_clearance_m=gap,
        required_parallel_clearance_m=context['safety_clearance_m'],
        support_shadow_fraction_upper_bound=mechanical['support_shadow_fraction_upper_bound'],
        bus_scenarios=bus,
        filters=dict(mechanical_converged=True,
            member_model_domain=member_ok,
            thin_wall_model_domain=mechanical['thin_wall_screen_passes'],
            substrate_stress=mechanical['substrate_stress_utilization'] <= 1,
            parallel_150m_clearance=gap >= context['safety_clearance_m'],
            component_mass_within_inherited_total_allowance=subtotal <= d['target_total_areal_mass_kg_m2'],
            geometry_and_member_screens=member_ok and gap >= context['safety_clearance_m']))


def pareto_ids(rows, *, bus=None):
    """Minimize subtotal and population; maximize the conditional clearance.

    Compare only a common load stratum supplied by the caller. Equal objective
    vectors survive together, making ties and duplicate material outcomes clear.
    """
    def coordinates(row):
        mass = (row['subtotal_kg_per_reference_gross_m2'] if bus is None else
                row['bus_scenarios'][bus]['subtotal_with_bus_kg_per_reference_gross_m2'])
        return (mass, row['equivalent_population'], -row['nominal_parallel_envelope_clearance_m'])
    ordered = sorted(rows, key=lambda row: (*coordinates(row), row['id']))
    front, points = [], []
    for row in ordered:
        point = np.asarray(coordinates(row))
        if points:
            old = np.asarray(points)
            if np.any(np.all(old <= point, axis=1) & np.any(old < point, axis=1)):
                continue
        front.append(row['id'])
        points.append(point)
    return front


def distributions(rows):
    result = {}
    for field in ('subtotal_kg_per_reference_gross_m2', 'equivalent_population',
                  'nominal_parallel_envelope_clearance_m', 'support_shadow_fraction_upper_bound'):
        values = [row[field] for row in rows if row['converged']]
        result[field] = (dict(zip(('minimum', 'p10', 'median', 'p90', 'maximum'),
            map(float, np.quantile(values, [0, .1, .5, .9, 1])))) if values else None)
    return result


def summarize(mechanical, rows, context):
    strata = []
    eligible = [row for row in rows if row['filters'].get('geometry_and_member_screens')]
    for factor in sorted({row['normal_load_factor'] for row in rows}):
        group = [row for row in eligible if row['normal_load_factor'] == factor]
        for bus in (None, 'fixed_inherited_equipment', 'area_scaled_equipment_hypothesis'):
            front = pareto_ids(group, bus=bus)
            strata.append(dict(normal_load_factor=factor, bus_scenario=bus or 'support_and_sheet_only',
                eligible_count=len(group), pareto_count=len(front), dominated_count=len(group)-len(front),
                pareto_ids=front))
    representatives = []
    for factor, side, material in product(sorted({r['normal_load_factor'] for r in rows}),
        sorted({r['side_m'] for r in rows}), sorted({r['compression_material_id'] for r in rows})):
        group = [r for r in rows if r['converged'] and r['normal_load_factor'] == factor
                 and r['side_m'] == side and r['compression_material_id'] == material]
        if not group:
            continue
        passing = [r for r in group if r['filters']['geometry_and_member_screens']]
        chosen = min(passing or group, key=lambda r: (r['subtotal_kg_per_reference_gross_m2'], r['id']))
        representatives.append(dict(selection='least_component_mass_in_geometry_member_stratum' if passing
            else 'least_component_mass_failed_geometry_or_member_screen', record=chosen))
    all_filters = sorted({key for row in rows for key in row['filters']})
    return dict(mechanical_count=len(mechanical), configuration_count=len(rows),
        converged_mechanical_count=sum(row['converged'] for row in mechanical),
        mechanical_failure_reasons=dict(Counter(row['failure'] for row in mechanical if not row['converged'])),
        filter_counts={name: sum(bool(row['filters'].get(name, False)) for row in rows) for name in all_filters},
        all_converged_distributions=distributions(rows),
        geometry_member_pass_distributions=distributions(eligible),
        pareto_objectives=['minimize normalized component or component-plus-assumed-bus mass',
                           'minimize equivalent population', 'maximize parallel-envelope clearance'],
        pareto_rule='within common normal-load and bus scenarios, after member and parallel-clearance filters; dominance is conditional on this model and is never physical impossibility',
        pareto_strata=strata, representatives=representatives)


def write_jsonl(path, rows):
    with path.open('w') as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n')
    return dict(path=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path),
        sha256=td.digest(path), rows=len(rows), bytes=path.stat().st_size,
        reading_rule='Every enumerated configuration is retained, including failures; ignored reproducible local run artifact.')


def run(*, archive_dir=ARCHIVE, summary_out=OUT, grid=None, shingle_angles=None,
        pitch_factors=None):
    start = time.perf_counter()
    archive_dir, summary_out = Path(archive_dir), Path(summary_out)
    grid = GRID if grid is None else grid
    angles = SHINGLES if shingle_angles is None else shingle_angles
    pitches = PITCH_FACTORS if pitch_factors is None else pitch_factors
    configurations = list(mechanical_configurations(grid))
    if len(configurations)*len(angles)*len(pitches) > 20000:
        raise ValueError('The cheap screen is capped at 20,000 configurations')
    data = load_inputs()
    context = context_from_inputs(data)
    mechanical = []
    for index, config in enumerate(configurations):
        mechanical.append(mechanical_case(config, context))
        if (index+1) % 270 == 0:
            print(json.dumps(dict(mechanical_completed=index+1, total=len(configurations),
                elapsed_s=time.perf_counter()-start)), flush=True)
    rows = [spacing_case(case, context, angle, pitch, i)
        for case in mechanical for i, (angle, pitch) in enumerate(product(angles, pitches))]
    summary = summarize(mechanical, rows, context)
    archive_dir.mkdir(parents=True, exist_ok=True)
    archives = dict(mechanical=write_jsonl(archive_dir/'mechanical.jsonl', mechanical),
                    configurations=write_jsonl(archive_dir/'configurations.jsonl', rows))
    result = dict(schema='terluna.research.array-geometry-screen/1',
        producer=dict(files={name: td.digest(ROOT/name) for name in FILES}, constants=constants_used(FILES)),
        inputs={path: dict(schema=schema, sha256=td.digest(ROOT/path)) for path, schema in INPUTS.values()},
        baseline_commit='669ddd3', evidence=EVIDENCE,
        reading_rule='Representative candidates and model failures only. No universal optimum or architecture-wide infeasibility claim. Load-factor and service-pose choices do not demonstrate feasible commands. All old products remain unchanged.',
        inherited_fields=context['inherited_values'],
        field_provenance=dict(design='baseline.hardware.variants.annulus_reference.design',
            lattice='baseline.fleet_service.ideal_service_lattice', population='baseline.fleet_service.population',
            safety_clearance='baseline.local_failures.cases.*.config.safety_clearance_m',
            normal_load='comparison.loading.resolved_normal_acceleration_m_s2',
            hub_fraction='comparison.command_loaded_supports[reference side/depth].support.hub_member.length_m / reference side'),
        grid=dict(mechanical=grid, shingle_deg=angles, pitch_factor=pitches),
        assumptions=dict(load_factors='scale all normal inertial/gravity coefficients of the inherited sampled command; solar pressure unchanged; lower factors are hypothetical load envelopes, not solved attitudes',
            load_density='full trim hardware density applied everywhere for the normal-load envelope; actual trim area fraction used for mass',
            geometry='scalloped support uses inherited dimension ratios; membrane sag stays at the inherited absolute tolerance',
            envelope='finite outer tube/brace/cable envelope plus prescribed membrane sag and member deflection allowance; joint protrusions and global deformation remain unbounded',
            coverage='sufficient inscribed rectangle projected at the scenario shingle angle; inherited pitch-to-active-rectangle ratios preserved; pitch factor below one increases overlap and population',
            clearance='parallel along-pitch neighbour plane separation minus two finite normal half-envelopes; 150 m inherited gate unchanged; no sliding adjacent-ring, fault or fleet collision replay',
            bus='fixed inherited avionics/housekeeping and hypothetical area-scaled equipment; power generation/storage/wiring, actuators, rotor supports and equipment-load closure omitted',
            classification='member and geometric screens organize follow-up candidates; conservative failed screens cannot reject all supports or attitudes'),
        archives=archives, summary=summary,
        runtime_s=time.perf_counter()-start, peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    summary_out.parent.mkdir(parents=True, exist_ok=True)
    summary_out.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(summary_out), mechanical_cases=len(mechanical),
        configurations=len(rows), runtime_s=result['runtime_s'], archive_bytes=sum(v['bytes'] for v in archives.values()))))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive-dir', type=Path, default=ARCHIVE)
    parser.add_argument('--summary-out', type=Path, default=OUT)
    args = parser.parse_args()
    run(archive_dir=args.archive_dir, summary_out=args.summary_out)
