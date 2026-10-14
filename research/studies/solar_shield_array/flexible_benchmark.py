"""Fixed-case verification and cost benchmark for the reduced flexible fleet.

This runner does not optimize a design or controller. All cases start from one
preserved support scenario and compare declared discretization changes. Their
short duration and sampled rays cannot establish orbital service or fleet safety.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import time

import numpy as np

from engineering.array_flexible import model_from_support
from shared import constants as K
from shared.provenance import constants_used
from . import tile_design as td

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/flexible_benchmark.json'
MODELS_OUT = HERE/'results/flexible_models.json'
ARCHIVE = ROOT/'research/runs/solar_shield_array_flexible_benchmark'
FILES = [
    'research/studies/solar_shield_array/flexible_benchmark.py',
    'engineering/array_flexible.py',
    'protection/dynamics/flexible_fleet.py',
    'protection/dynamics/flexible_optics.py',
    'protection/dynamics/flexible_environment.py',
    'shared/provenance.py']
MODEL_FILES = ['research/studies/solar_shield_array/flexible_benchmark.py',
               'engineering/array_flexible.py', 'shared/provenance.py']
INPUTS = {
    'baseline': ('research/studies/solar_shield_array/results/tile_design.json',
                 'terluna.research.array-tile-design/1'),
    'support': ('research/studies/solar_shield_array/results/support_comparison.json',
                'terluna.research.array-support-comparison/1'),
    'materials': ('research/studies/solar_shield_array/support_materials.json',
                  'terluna.research.array-support-materials/1'),
    'support_sources': ('research/studies/solar_shield_array/support_sources.json',
                       'terluna.research.array-support-sources/1'),
    'flexible_sources': ('engineering/flexible_sources.json',
                        'terluna.research.flexible-array-source-notes/1'),
    'optics_sources': ('research/studies/solar_shield_array/flexible_optics_sources.json',
                      'terluna.research.flexible-optics-source-notes/1'),
    'film': ('protection/spectra/annulus_film.json', 'terluna.protection.annulus-film/1')}
EVIDENCE = ('Executed fixed-case reduced-model numerical comparisons; no optimized '
    'candidate, empirically validated structure, full control/hardware closure, '
    'long-duration orbital service, continuous coverage or collision certificate.')
MODEL_READING_RULE = ('SI units. Nodes in metres about the nodal centre of mass, nodal masses in kg, normal-only '
    'modes mass-normalized in kg^(-1/2) with modal coordinates in m sqrt(kg), omega_squared in s^-2. Compare '
    'physical probes across meshes and bases, not modal coordinates. Support masses are lumped at member '
    'endpoints (metadata.nodal_inertia_relative_frobenius_error against the distributed reference) and the hub '
    'restraint is an effective stiffness; read metadata.mass_rule, support_rule and geometry_rule before using the '
    'model for attitude or control.')
ABSOLUTE_TOLERANCES = dict(position_m=1e-5, velocity_m_s=1e-8,
    quaternion=1e-11, angular_velocity_rad_s=1e-13,
    modal_peak_displacement_m=1e-7, modal_peak_velocity_m_s=1e-9,
    wheel_momentum_N_m_s=1e-4, work_J=1e-3, angular_impulse_N_m_s=.01)


class BenchmarkBudgetExceeded(RuntimeError):
    """A resource stop carries an attempted RHS stage, never an accepted state."""
    def __init__(self, elapsed_s, rhs_calls, stage_time_s, stage_state):
        super().__init__('Per-case RHS wall-time budget exceeded; no accepted integration result')
        self.context = dict(elapsed_wall_s=elapsed_s, rhs_evaluations=rhs_calls,
            last_evaluated_stage_time_s=float(stage_time_s),
            last_evaluated_stage_state=np.asarray(stage_state).tolist(),
            state_interpretation='Attempted integrator stage, not an accepted physical trajectory or completed prefix.')


def canonical_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def load_inputs():
    data = {}
    for key, (path, schema) in INPUTS.items():
        data[key] = json.loads((ROOT/path).read_text())
        if data[key].get('schema') != schema:
            raise ValueError(f'Unexpected input schema: {path}')
    return data


def select_reference(data):
    """A fixed prior case, independent of the new screen's candidate ordering."""
    reference = data['baseline']['hardware']['variants']['annulus_reference']
    design = reference['design']
    candidates = [row for row in data['support']['command_loaded_supports']
        if row['side_m'] == design['side_m'] and row['edge_retreat_fraction'] == .05
        and row['boom_depth_m'] == design['frame_depth_m']
        and row['compression_material_id'] == 'm55j_2500_ud'
        and row['cable_material_id'] == 'kevlar49_yarn']
    if len(candidates) != 1 or not candidates[0]['converged']:
        raise ValueError('Expected exactly one converged fixed reference support')
    films = [row for row in data['film']['films'] if math.isclose(
        row['areal_mass_g_m2']/1000, design['film_areal_mass_kg_m2'], rel_tol=1e-12)]
    if len(films) != 1:
        raise ValueError('Film mass must identify exactly one inherited optical row')
    optics = films[0]['bolometric']
    if not math.isclose(sum(optics[name] for name in ('A', 'R', 'T')), 1., abs_tol=1e-12):
        raise ValueError('Inherited optical fractions do not close')
    sigma = (design['film_areal_mass_kg_m2']+design['substrate_areal_mass_kg_m2']+
             design['trim_area_fraction']*design['trim_hardware_areal_mass_kg_m2'])
    return dict(support=candidates[0]['support'], support_case_id=candidates[0]['id'],
        design=design, sheet_areal_mass_kg_m2=sigma, optical_fractions=optics,
        optical_row_description={key: films[0][key] for key in ('titania_um', 'silica_um', 'coated', 'areal_mass_g_m2')},
        radius_m=reference['assembled_momentum_screen']['radius_m'],
        hub_equipment_mass_kg=design['avionics_and_sensors_kg'],
        pitch_m=data['baseline']['fleet_service']['ideal_service_lattice']['pitch_m'])


def fixed_cases():
    """Declared cases only; refinements alter one discretization at a time."""
    base = dict(mesh_divisions=8, retained_modes=6, source_rings=0,
        max_step_s=10., output_step_s=5., duration_s=60., tiles=1,
        gravity=True, photons=True, controller='none',
        rtol=2e-8, absolute_tolerance_profile='dimensioned_v1', method='DOP853',
        integration_wall_budget_s=300.)
    specifications = [
        ('reference', {}),
        ('mesh12', dict(mesh_divisions=12)),
        ('mesh16', dict(mesh_divisions=16)),
        ('modes12', dict(mesh_divisions=12, retained_modes=12)),
        ('modes20', dict(mesh_divisions=12, retained_modes=20)),
        ('modes40', dict(mesh_divisions=12, retained_modes=40)),
        ('modes80', dict(mesh_divisions=12, retained_modes=80)),
        ('modes170', dict(mesh_divisions=12, retained_modes=170)),
        ('full_mesh8', dict(mesh_divisions=8, retained_modes=82)),
        ('full_mesh16', dict(mesh_divisions=16, retained_modes=290)),
        ('full_mesh16_step', dict(mesh_divisions=16, retained_modes=290, max_step_s=.5)),
        ('step5', dict(max_step_s=5.)),
        ('step2_5', dict(max_step_s=2.5)),
        ('source7', dict(source_rings=1)),
        ('source19', dict(source_rings=2)),
        ('cost1', dict(duration_s=30., max_step_s=5.)),
        ('cost3', dict(duration_s=30., max_step_s=5., tiles=3)),
        ('solver_rk45', dict(duration_s=30., max_step_s=5., method='RK45')),
        ('solver_radau', dict(duration_s=30., max_step_s=5., method='Radau')),
        ('tolerance_scalar', dict(duration_s=2., max_step_s=1., output_step_s=1.,
                                  absolute_tolerance_profile='legacy_scalar_1e-9')),
        ('tolerance_dimensioned', dict(duration_s=2., max_step_s=1., output_step_s=1.)),
        ('unforced_energy', dict(gravity=False, photons=False)),
        ('bounded_photon_trim', dict(controller='photon_trim_pd')),
        ('lowest_mode_period', dict(period_probe=True))]
    return [dict(base, **changes, id=name) for name, changes in specifications]


def export_model(model, reference, producer, inputs):
    """Explicit geometry/mass/mode bytes crossing the engineering boundary."""
    masses = np.asarray(model.nodal_masses_kg)
    centre = np.average(model.nodes_m, axis=0, weights=masses)
    payload = dict(nodes_m=(np.asarray(model.nodes_m)-centre).tolist(),
        triangles=np.asarray(model.triangles).tolist(), nodal_masses_kg=masses.tolist(),
        modes=np.asarray(model.modes).tolist(),
        omega_squared_s2=np.asarray(model.omega_squared_s2).tolist(),
        absorptivity=reference['optical_fractions']['A'],
        reflectivity=reference['optical_fractions']['R'], damping_ratio=0.,
        envelope_margin_m=math.sqrt(reference['support']['half_depth_m']**2+
            2*((reference['support']['outside_span_m']-reference['design']['side_m'])/2)**2),
        metadata=dict(model.metadata, centre_subtracted_m=centre.tolist(),
            corner_indices=np.asarray(model.corner_indices).tolist(),
            membrane_node_count=model.membrane_node_count,
            source_support_case_id=reference['support_case_id'],
            material_catalogue_qualification='none; preserved conditional material assumptions',
            whole_spacecraft_budget_closed=False))
    return dict(schema='terluna.engineering.array-flexible/1',
        model=payload, model_sha256=canonical_hash(payload),
        producer_hashes=producer['files'], constants=producer['constants'],
        input_hashes={path: value['sha256'] for path, value in inputs.items()},
        evidence='Explicit numerical reduced-model input, not empirical material or hardware validation.',
        reading_rule=MODEL_READING_RULE)


def model_key(case):
    return f"mesh{case['mesh_divisions']}_modes{case['retained_modes']}"


def build_models(reference, cases, *, output=MODELS_OUT):
    start = time.perf_counter()
    provenance = dict(files={name: td.digest(ROOT/name) for name in MODEL_FILES},
                      constants=constants_used(MODEL_FILES))
    inputs = {path: dict(schema=schema, sha256=td.digest(ROOT/path))
              for path, schema in INPUTS.values()}
    models, resources = {}, {}
    for case in cases:
        key = model_key(case)
        if key in models:
            continue
        before = time.perf_counter()
        model = model_from_support(reference['support'], reference['sheet_areal_mass_kg_m2'],
            mesh_divisions=case['mesh_divisions'], retained_modes=case['retained_modes'],
            hub_equipment_mass_kg=reference['hub_equipment_mass_kg'])
        models[key] = export_model(model, reference, provenance, inputs)
        resources[key] = dict(build_wall_s=time.perf_counter()-before,
            nodes=len(model.nodes_m), triangles=len(model.triangles), retained_modes=len(model.omega_squared_s2))
    collection = dict(schema='terluna.research.array-flexible-model-inputs/1',
        producer=provenance, inputs=inputs, models=models,
        reference_case=dict(support_case_id=reference['support_case_id'],
            sheet_areal_mass_kg_m2=reference['sheet_areal_mass_kg_m2'],
            hub_equipment_mass_kg=reference['hub_equipment_mass_kg'],
            optical_fractions=reference['optical_fractions'],
            optical_row_description=reference['optical_row_description']),
        resources=dict(models=resources, total_build_wall_s=time.perf_counter()-start),
        evidence='Explicit fixed-case reduced input products for implementation benchmarks; no design selection.',
        reading_rule=('One reduced structural input per mesh and retained-mode count of the fixed benchmark cases. '
                      'Each model carries its own schema, hash, provenance, evidence and reading rule; check its '
                      'model_sha256 before use.'))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(collection, separators=(',', ':'), allow_nan=False)+'\n')
    return collection


def initial_modal_bump(payload, side_m, amplitude_m=.01):
    """Project one fixed smooth asymmetric physical displacement into each basis."""
    nodes, modes = np.asarray(payload['nodes_m']), np.asarray(payload['modes'])
    masses = np.asarray(payload['nodal_masses_kg'])
    distance = ((nodes[:, 0]-.15*side_m)**2+(nodes[:, 1]+.1*side_m)**2)/(side_m/5)**2
    bump = amplitude_m*np.exp(-distance)
    return modes[:, 2, :].T@(masses*bump)


def finite_json(value, path='', special=None):
    """Retain diagnostic nonfinite locations explicitly instead of emitting NaN."""
    if special is None:
        special = []
    if isinstance(value, np.ndarray):
        return finite_json(value.tolist(), path, special)
    if isinstance(value, np.generic):
        return finite_json(value.item(), path, special)
    if isinstance(value, dict):
        return {str(key): finite_json(item, f'{path}.{key}' if path else str(key), special)
                for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [finite_json(item, f'{path}[{index}]', special) for index, item in enumerate(value)]
    if isinstance(value, float) and not math.isfinite(value):
        special.append(dict(path=path, value=str(value)))
        return None
    return value


def difference_metrics(first, second):
    """Absolute differences remain visible when a reference quantity is tiny."""
    result = {}
    for key in ('last_position_m', 'last_velocity_m_s', 'last_corner_positions_m',
                'last_shape_probe_displacements_m', 'last_total_photon_force_N'):
        a, b = np.asarray(first[key]), np.asarray(second[key])
        if a.shape != b.shape:
            raise ValueError(f'Incomparable refinement observables: {key}')
        difference = float(np.max(np.linalg.norm((a-b).reshape(-1, a.shape[-1]), axis=1)))
        scale = float(np.max(np.linalg.norm(b.reshape(-1, b.shape[-1]), axis=1)))
        result[key] = dict(maximum_absolute_norm_difference=difference,
            relative_to_second_maximum_norm=difference/scale if scale else None,
            second_maximum_norm=scale)
    return result


def aggregate_refinements(cases):
    by_id = {row['id']: row for row in cases}
    sequences = dict(mesh=['reference', 'mesh12', 'mesh16'],
        retained_modes=['mesh12', 'modes12', 'modes20'],
        extended_retained_modes=['modes20', 'modes40', 'modes80', 'modes170'],
        complete_basis_mesh=['full_mesh8', 'modes170', 'full_mesh16'],
        complete_basis_time=['full_mesh16', 'full_mesh16_step'],
        maximum_step=['reference', 'step5', 'step2_5'],
        solar_source=['reference', 'source7', 'source19'])
    result = {}
    for name, sequence in sequences.items():
        pairs = []
        for a, b in zip(sequence, sequence[1:]):
            first, second = by_id[a], by_id[b]
            comparable = (first['integration_reached_requested_end'] and second['integration_reached_requested_end'])
            pairs.append(dict(first=a, second=b, comparable=comparable,
                differences=difference_metrics(first['observables'], second['observables']) if comparable else None,
                initial_projection_differences=difference_metrics(first['initial_observables'], second['initial_observables']) if comparable else None))
        result[name] = dict(cases=sequence, adjacent_pairs=pairs,
            reading_rule='Observed fixed-case discretization sensitivity over this short prefix, not a convergence order, error bound or general accuracy guarantee.')
    return result


def resolution_assessment(cases):
    """Declare probe tolerances and show modal error against a complete basis.

    This is a chosen numerical comparison criterion, not a material, optical
    or flight requirement. Full basis removes modal truncation only for its
    particular discrete normal finite-element model.
    """
    by_id = {row['id']: row for row in cases}
    absolute_tolerance_m, relative_tolerance = 1e-4, .01
    def compare(first, second):
        a, b = by_id[first], by_id[second]
        if not (a['integration_reached_requested_end'] and b['integration_reached_requested_end']):
            return dict(first=first, second=second, comparable=False,
                        meets_declared_probe_tolerance=False)
        metric = difference_metrics(a['observables'], b['observables'])['last_shape_probe_displacements_m']
        tolerance = absolute_tolerance_m+relative_tolerance*metric['second_maximum_norm']
        return dict(first=first, second=second, comparable=True,
            maximum_shape_probe_difference_m=metric['maximum_absolute_norm_difference'],
            reference_maximum_shape_probe_displacement_m=metric['second_maximum_norm'],
            declared_probe_tolerance_m=tolerance,
            meets_declared_probe_tolerance=metric['maximum_absolute_norm_difference'] <= tolerance)
    modal = [compare(name, 'modes170') for name in
             ('mesh12', 'modes12', 'modes20', 'modes40', 'modes80', 'modes170')]
    mesh = [compare('full_mesh8', 'modes170'), compare('modes170', 'full_mesh16')]
    time = compare('full_mesh16', 'full_mesh16_step')
    return dict(absolute_probe_tolerance_m=absolute_tolerance_m,
        relative_probe_tolerance=relative_tolerance,
        criterion='maximum common physical shape-probe difference <= 0.1 mm + 1% of reference peak probe displacement; a declared benchmark criterion, not a physical requirement or certified error bound',
        fixed_mesh12_modal_errors_against_complete_basis=modal,
        complete_basis_mesh_comparisons=mesh,
        complete_basis_time_comparison=time,
        complete_basis_definition='all N-3 elastic normal degrees of freedom after three affine rigid normal modes are removed: mesh8=82, mesh12=170, mesh16=290',
        final_complete_basis_mesh_pair_meets_declared_probe_tolerance=mesh[-1]['meets_declared_probe_tolerance'],
        complete_basis_time_pair_meets_declared_probe_tolerance=time['meets_declared_probe_tolerance'],
        reading_rule='Do not infer mesh convergence from a six-mode mesh sequence. Full discrete bases remove modal truncation for each mesh, but leave mesh geometry, assumed stiffness, endpoint inertia and all physical-model limitations.')


def solver_comparison(cases):
    by_id = {row['id']: row for row in cases}
    reference = by_id['cost1']
    rows = []
    for name in ('cost1', 'solver_rk45', 'solver_radau'):
        row = by_id[name]
        comparable = row['integration_reached_requested_end'] and reference['integration_reached_requested_end']
        rows.append(dict(id=name, method=row['configuration']['method'], resources=row['resources'],
            completed=row['integration_reached_requested_end'],
            terminal_difference_from_DOP853=difference_metrics(row['observables'], reference['observables']) if comparable else None))
    return dict(cases=rows, physical_duration_s=30.,
        rule='Same geometry, state, loads, requested tolerances and maximum step. Agreement between methods is a numerical check, not a validated error bound. Cost is measured for this fixed prefix only.')


def tolerance_comparison(cases):
    by_id = {row['id']: row for row in cases}
    scalar, dimensioned = by_id['tolerance_scalar'], by_id['tolerance_dimensioned']
    comparable = scalar['integration_reached_requested_end'] and dimensioned['integration_reached_requested_end']
    return dict(scalar_case=scalar['id'], dimensioned_case=dimensioned['id'], comparable=comparable,
        physical_observable_differences=difference_metrics(scalar['observables'], dimensioned['observables']) if comparable else None,
        scalar_over_dimensioned_wall_ratio=scalar['resources']['integration_wall_s']/dimensioned['resources']['integration_wall_s'],
        scalar_rhs_evaluations=scalar['resources']['rhs_evaluations'],
        dimensioned_rhs_evaluations=dimensioned['resources']['rhs_evaluations'],
        final_work_difference_J=scalar['diagnostic_end']['work_J']-dimensioned['diagnostic_end']['work_J'] if comparable else None,
        final_external_angular_impulse_difference_N_m_s=(np.asarray(scalar['diagnostic_end']['angular_impulse_N_m_s'])-
            np.asarray(dimensioned['diagnostic_end']['angular_impulse_N_m_s'])).tolist() if comparable else None,
        reading_rule='A reproducible two-second tolerance/startup check. Dimensioned absolute tolerances express distinct physical units; they do not relax rtol or establish general physical accuracy.')


def resource_comparison(cases):
    by_id = {row['id']: row for row in cases}
    one, three = by_id['cost1'], by_id['cost3']
    return dict(one_tile=one['resources'], three_tiles=three['resources'],
        wall_ratio_three_over_one=three['resources']['integration_wall_s']/one['resources']['integration_wall_s'] if one['resources']['integration_wall_s'] else None,
        rhs_ratio_three_over_one=three['resources']['rhs_evaluations']/one['resources']['rhs_evaluations'] if one['resources']['rhs_evaluations'] else None,
        both_completed=one['integration_reached_requested_end'] and three['integration_reached_requested_end'],
        comparable_requested_duration_s=one['configuration']['duration_s'],
        rule='Measured one/three-tile fixed-prefix cost only; no linear large-fleet extrapolation, optimizer budget or accepted service duration.')


def period_configuration(case, model_product, measured_reference):
    """One fixed physical-period probe, capped only by an explicit cost estimate."""
    frequencies = np.sqrt(model_product['model']['omega_squared_s2'])
    period = 2*math.pi/float(frequencies.min())
    maximum_step = min(period/100, 2*math.pi/float(frequencies.max())/12)
    resource = measured_reference['resources']
    cost_per_rhs = resource['integration_wall_s']/max(resource['rhs_evaluations'], 1)
    estimated_wall = (12*math.ceil(period/maximum_step)+50)*cost_per_rhs
    limited = estimated_wall > 300
    duration = min(period, 1000.) if limited else period
    result = dict(case, duration_s=duration, max_step_s=min(maximum_step, duration/100),
        output_step_s=duration/100, lowest_mode_period_s=period,
        requested_period_fraction=duration/period,
        full_period_wall_estimate_s=estimated_wall, period_cost_cap_applied=limited,
        period_cost_rule='12 DOP853 stages per bounded step plus 50 startup calls, priced at the measured reference cost per RHS; if projected wall time exceeds 300 s, requested duration is capped at 1000 s. This is a rough resource estimate, not solver accuracy.')
    return result


def probe_interpolation(payload, side_m):
    """P1 interpolation at common body-plane locations across different meshes."""
    probes = np.array([[x, y] for x in (-.25*side_m, 0., .25*side_m)
                       for y in (-.25*side_m, 0., .25*side_m)])
    nodes = np.asarray(payload['nodes_m'])
    triangles = np.asarray(payload['triangles'])
    vertices = nodes[triangles, :2]
    first, second = vertices[:, 1]-vertices[:, 0], vertices[:, 2]-vertices[:, 0]
    determinant = first[:, 0]*second[:, 1]-first[:, 1]*second[:, 0]
    indices, weights = [], []
    for probe in probes:
        delta = probe-vertices[:, 0]
        u = (delta[:, 0]*second[:, 1]-delta[:, 1]*second[:, 0])/determinant
        v = (first[:, 0]*delta[:, 1]-first[:, 1]*delta[:, 0])/determinant
        found = np.flatnonzero((u >= -1e-10) & (v >= -1e-10) & (u+v <= 1+1e-10))
        if not len(found):
            raise ValueError('A physical comparison probe lies outside the membrane mesh')
        index = found[0]
        indices.append(triangles[index]); weights.append([1-u[index]-v[index], u[index], v[index]])
    return probes, np.asarray(indices), np.asarray(weights)


def target_centre(t, index, case, reference):
    offset = (index-(case['tiles']-1)/2)*reference['pitch_m'][0]
    radius = reference['radius_m']+offset
    if case['gravity']:
        angular = math.sqrt(K.MOON_GM/radius**3)
        phase = angular*t
        return (np.array([radius*math.cos(phase), radius*math.sin(phase), 0.]),
                np.array([-radius*angular*math.sin(phase), radius*angular*math.cos(phase), 0.]))
    return np.array([radius, 0., 0.]), np.zeros(3)


def state_absolute_tolerances(bodies):
    values = ABSOLUTE_TOLERANCES
    vectors = []
    for body in bodies:
        mode_max = np.max(np.linalg.norm(body.modes, axis=1), axis=0)
        vectors.append(np.r_[np.full(3, values['position_m']),
            np.full(3, values['velocity_m_s']), np.full(4, values['quaternion']),
            np.full(3, values['angular_velocity_rad_s']),
            values['modal_peak_displacement_m']/mode_max,
            values['modal_peak_velocity_m_s']/mode_max,
            np.full(3, values['wheel_momentum_N_m_s']), values['work_J'],
            np.full(3, values['angular_impulse_N_m_s'])])
    return np.concatenate(vectors)


def execute_case(case, model_product, reference, archive_dir):
    from protection.dynamics.flexible_fleet import Body, Fleet, rotation
    from protection.dynamics.flexible_optics import solar_disk_quadrature, receiver_rays

    started = time.perf_counter()
    cpu_start = time.process_time()
    payload = model_product['model']
    side = reference['design']['side_m']
    optical = reference['optical_fractions']
    trimmed = case['controller'] == 'photon_trim_pd'
    options = dict(absorptivity=optical['A'], reflectivity=optical['R'],
        envelope_margin_m=payload['envelope_margin_m'],
        trim_max_reflectivity=min(1-optical['A'], optical['R']+.1) if trimmed else None,
        trim_fraction=reference['design']['trim_area_fraction'])
    bodies = [Body.from_product(model_product, **options) for _ in range(case['tiles'])]
    directions, weights = solar_disk_quadrature([0., 0., 1.],
        math.asin(K.SUN_RADIUS/K.AU), rings=case['source_rings'])
    pressure = (K.SOLAR_CONSTANT/K.SPEED_OF_LIGHT/float(weights@directions[:, 2])) if case['photons'] else 0.
    # This reference shingle is the prior attitude's initial value. The present
    # free-body integration does not impose the prior large roll command.
    shingle = math.radians(1.5)
    quaternion = np.array([math.cos(shingle/2), 0., math.sin(shingle/2), 0.])
    target_frame = rotation(quaternion)
    q = initial_modal_bump(payload, side)
    initial = []
    for index in range(case['tiles']):
        position, velocity = target_centre(0., index, case, reference)
        initial.append(dict(position_m=position, velocity_m_s=velocity,
            quaternion_wxyz=quaternion, omega_rad_s=[1e-5, -.5e-5, .25e-5], q=q))
    initial_inertia = np.asarray(payload['metadata']['nodal_reference_inertia_kg_m2'])
    nominal_force = np.array([0., 0., -pressure*reference['support']['geometry']['retained_area_m2']*(optical['A']+2*optical['R'])])
    def controller(t, index, state):
        target_position, target_velocity = target_centre(t, index, case, reference)
        # A fixed 5% larger photon acceleration is a tracking test trajectory,
        # not an orbital design or a guaranteed feasible six-axis demand.
        desired_acceleration = 1.05*nominal_force/bodies[index].mass_kg
        target_position += .5*desired_acceleration*t*t
        target_velocity += desired_acceleration*t
        force = bodies[index].mass_kg*(desired_acceleration+
            1e-4*(target_position-state['position_m'])+.02*(target_velocity-state['velocity_m_s']))
        error = target_frame.T@rotation(state['quaternion_wxyz'])
        angle_error = .5*np.array([error[2, 1]-error[1, 2], error[0, 2]-error[2, 0], error[1, 0]-error[0, 1]])
        torque = -initial_inertia@(1e-7*angle_error+1e-4*state['omega_rad_s'])
        return dict(force_world_N=force, torque_body_N_m=torque)
    fleet = Fleet(bodies, source_directions=directions, source_weights=weights,
        pressure_Pa=pressure, gravitational_parameter_m3_s2=K.MOON_GM if case['gravity'] else 0.,
        controller=controller if trimmed else None, self_shadow=True)
    y0 = fleet.initial_state(initial)
    integrate_start = time.perf_counter()
    absolute_tolerances = (np.full(len(y0), 1e-9) if case['absolute_tolerance_profile'] == 'legacy_scalar_1e-9'
                           else state_absolute_tolerances(bodies))
    original_rhs = fleet.rhs
    def bounded_rhs(t, y):
        elapsed = time.perf_counter()-integrate_start
        if elapsed > case['integration_wall_budget_s']:
            raise BenchmarkBudgetExceeded(elapsed, fleet.rhs_calls, t, y)
        return original_rhs(t, y)
    fleet.rhs = bounded_rhs
    solution = fleet.integrate(y0, case['duration_s'], output_step_s=case['output_step_s'],
        max_step_s=case['max_step_s'], rtol=case['rtol'], atol=absolute_tolerances, method=case['method'])
    integration_wall = time.perf_counter()-integrate_start
    probes, probe_nodes, probe_weights = probe_interpolation(payload, side)
    corners = np.asarray(payload['metadata']['corner_indices'])
    diagnostics, sampled, observable_history = [], [], []
    xextent = .4*side+(case['tiles']-1)/2*reference['pitch_m'][0]
    receiver_offsets = np.array([[x, y, -side] for x in np.linspace(-xextent, xextent, 5)
                                for y in np.linspace(-.4*side, .4*side, 5)])
    for t, y in zip(solution.t, solution.y.T):
        evaluation = fleet.evaluate(t, y)
        diagnostic = fleet.diagnostics(t, y)
        diagnostics.append(diagnostic)
        centre, _ = target_centre(t, (case['tiles']-1)/2, case, reference)
        receiver = centre+receiver_offsets
        coverage = receiver_rays(evaluation['meshes'],
            np.repeat(receiver, len(directions), axis=0), np.tile(directions, (len(receiver), 1)))
        covered = coverage['geometric_covered'].reshape(len(receiver), len(directions))
        transmission = coverage['energy_transmittance'].reshape(len(receiver), len(directions))
        controls = evaluation['controls']
        sampled.append(dict(time_s=float(t),
            minimum_receiver_source_coverage_fraction=float(np.min(covered@weights)),
            mean_receiver_source_coverage_fraction=float(np.mean(covered@weights)),
            maximum_receiver_bolometric_transmission=float(np.max(transmission@weights)),
            maximum_force_request_residual_N=max(float(np.linalg.norm(c['force_residual_N'])) for c in controls),
            maximum_torque_request_residual_N_m=max(float(np.linalg.norm(c['torque_residual_N_m'])) for c in controls),
            maximum_trim_fraction_of_available_range=max(float(np.max(c['trim_command'])) for c in controls),
            allocation_success_by_tile=[c['allocation_success'] for c in controls],
            allocation_status_by_tile=[c['allocation_status'] for c in controls],
            allocation_scaled_residual_N_by_tile=[c['allocation_scaled_residual_N'] for c in controls],
            photon_force_balance_residual_N=float(np.linalg.norm(evaluation['photon_result']['force_balance_residual_N'])),
            photon_torque_balance_residual_N_m=float(np.linalg.norm(evaluation['photon_result']['torque_balance_residual_N_m'])),
            max_slope=diagnostic['max_slope'], max_deformation_m=diagnostic['max_deformation_m'],
            clearance_lower_bound_m=diagnostic['clearance_lower_bound_m']))
        shape = [np.einsum('pi,pij->pj', probe_weights,
            (local-body.nodes_m)[probe_nodes]) for body, local in zip(bodies, evaluation['local'])]
        observable_history.append(dict(last_position_m=[state['position_m'].tolist() for state in evaluation['states']],
            last_velocity_m_s=[state['velocity_m_s'].tolist() for state in evaluation['states']],
            last_corner_positions_m=[world[corners].tolist() for world in evaluation['world']],
            last_shape_probe_displacements_m=np.asarray(shape).tolist(),
            last_total_photon_force_N=evaluation['photon_result']['total_force'].tolist()))
    energy = np.array([d['energy_minus_work_J'] for d in diagnostics])
    momentum = np.array([d['angular_minus_impulse_N_m_s'] for d in diagnostics])
    delta_energy = float(np.max(np.abs(energy-energy[0])))
    delta_momentum = float(np.max(np.linalg.norm(momentum-momentum[0], axis=1)))
    archive_path = Path(archive_dir)/(case['id']+'.npz')
    np.savez_compressed(archive_path, times_s=solution.t, states=solution.y,
                        body_offsets=fleet.offsets)
    special = []
    result = dict(id=case['id'], configuration=case, model_key=model_key(case),
        model_sha256=model_product['model_sha256'], integration_success=bool(solution.success),
        integration_reached_requested_end=bool(solution['completed_duration']),
        validity_stop=bool(solution['validity_stop']), solver_message=solution.message,
        stop_reason=solution['stop_reason'],
        last_sampled_time_s=float(solution.t[-1]),
        completed_lowest_mode_period_fraction=float(solution.t[-1])/(2*math.pi/math.sqrt(min(payload['omega_squared_s2']))),
        validity_event_times_s=[np.asarray(event).tolist() for event in solution.t_events],
        accepted_operating_design=False, continuous_coverage_demonstrated=False,
        complete_spacecraft_budget=False, optimized=False,
        pressure_Pa=pressure, source_count=len(directions), initial_bump_amplitude_m=.01,
        absolute_tolerances=dict(physical_settings=ABSOLUTE_TOLERANCES if case['absolute_tolerance_profile'] == 'dimensioned_v1' else None,
            effective_solver_vector=absolute_tolerances.tolist(),
            rule='Per-mode absolute tolerance bounds that mode alone by the stated maximum physical displacement/velocity; combined truncation and solution error are not bounded by this setting.'),
        shape_probe_body_xy_m=probes.tolist(), receiver_count=len(receiver_offsets),
        model_approximation=dict(area_relative_error=payload['metadata']['area_relative_error'],
            nodal_inertia_relative_frobenius_error=payload['metadata']['nodal_inertia_relative_frobenius_error'],
            elastic_normal_degrees_of_freedom=len(payload['nodes_m'])-3,
            retained_modes=len(payload['omega_squared_s2']),
            complete_discrete_normal_basis=len(payload['omega_squared_s2']) == len(payload['nodes_m'])-3,
            highest_retained_angular_frequency_rad_s=math.sqrt(max(payload['omega_squared_s2'])),
            note='Polygon boundary may overstate coverage; endpoint-lumped support inertia is a model approximation that mesh/mode refinement does not remove.'),
        control=dict(kind=case['controller'], wheel_hardware_included=False,
            trim_max_reflectivity=options['trim_max_reflectivity'], trim_area_fraction=options['trim_fraction'],
            control_law_rule='diagnostic grey-reflection PD requests; actual finite photon allocation and residuals reported; no UV, stability or hardware qualification'),
        conservation=dict(maximum_energy_minus_work_drift_J=delta_energy,
            drift_fraction_of_initial_total_energy=delta_energy/max(abs(energy[0]), 1e-300),
            maximum_angular_minus_impulse_drift_N_m_s=delta_momentum,
            drift_fraction_of_initial_total_angular_momentum=delta_momentum/max(float(np.linalg.norm(momentum[0])), 1e-300),
            rule='Point-mass/model energy and angular balance; orbital totals can hide small elastic errors. The unforced case removes orbital translation/gravity but still tests only the retained model.'),
        diagnostic_start=diagnostics[0], diagnostic_end=diagnostics[-1],
        maximum_sampled_slope=max(d['max_slope'] for d in diagnostics),
        maximum_sampled_deformation_m=max(d['max_deformation_m'] for d in diagnostics),
        maximum_quaternion_norm_error=max(d['max_quaternion_norm_error'] for d in diagnostics),
        sampled_history=sampled, initial_observables=observable_history[0], observables=observable_history[-1],
        state_archive=dict(path=str(archive_path.relative_to(ROOT)) if archive_path.is_relative_to(ROOT) else str(archive_path),
            sha256=td.digest(archive_path), bytes=archive_path.stat().st_size),
        resources=dict(integration_wall_s=integration_wall,
            diagnostics_and_archive_wall_s=time.perf_counter()-integrate_start-integration_wall,
            total_case_wall_s=time.perf_counter()-started, process_cpu_s=time.process_time()-cpu_start,
            rhs_evaluations=int(solution.nfev), returned_samples=len(solution.t),
            counted_rhs_evaluations=int(solution['rhs_calls']),
            jacobian_evaluations=int(solution.njev), lu_decompositions=int(solution.nlu),
            state_variables=len(y0), triangles=case['tiles']*len(bodies[0].triangles),
            nodes=case['tiles']*len(bodies[0].nodes_m), source_directions=len(directions),
            peak_process_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
    result = finite_json(result, special=special)
    result['nonfinite_diagnostic_fields'] = special
    return result


def run(*, archive_dir=ARCHIVE, output=OUT, models_output=MODELS_OUT):
    started = time.perf_counter()
    cpu_started = time.process_time()
    archive_dir, output = Path(archive_dir), Path(output)
    archive_dir.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    reference = select_reference(data)
    plan = fixed_cases()
    model_collection = build_models(reference, plan, output=models_output)
    cases = []
    for configuration in plan:
        if configuration.get('period_probe'):
            configuration.update(period_configuration(configuration,
                model_collection['models'][model_key(configuration)],
                next(row for row in cases if row['id'] == 'reference')))
        start = time.perf_counter()
        print(json.dumps(dict(starting_case=configuration['id'], elapsed_s=time.perf_counter()-started)), flush=True)
        try:
            row = execute_case(configuration, model_collection['models'][model_key(configuration)], reference, archive_dir)
        except BenchmarkBudgetExceeded as exc:
            row = dict(id=configuration['id'], configuration=configuration,
                integration_reached_requested_end=False, integration_success=False,
                failure_type=type(exc).__name__, failure=str(exc),
                accepted_operating_design=False, continuous_coverage_demonstrated=False,
                complete_spacecraft_budget=False, optimized=False, stop_reason='wall_time_budget',
                last_evaluated_stage=exc.context,
                resources=dict(integration_wall_s=exc.context['elapsed_wall_s'],
                    rhs_evaluations=exc.context['rhs_evaluations'], cost_is_failed_attempt=True))
        except (ValueError, RuntimeError, np.linalg.LinAlgError, FloatingPointError) as exc:
            row = dict(id=configuration['id'], configuration=configuration,
                integration_reached_requested_end=False, integration_success=False,
                failure_type=type(exc).__name__, failure=str(exc),
                accepted_operating_design=False, continuous_coverage_demonstrated=False,
                complete_spacecraft_budget=False, optimized=False,
                resources=dict(integration_wall_s=time.perf_counter()-start, rhs_evaluations=0,
                               cost_is_failed_attempt=True))
        row.setdefault('model_key', model_key(configuration))
        row.setdefault('model_sha256', model_collection['models'][model_key(configuration)]['model_sha256'])
        cases.append(row)
        print(json.dumps(dict(finished_case=row['id'], completed=row['integration_reached_requested_end'],
                              wall_s=time.perf_counter()-start)), flush=True)
    inputs = {path: dict(schema=schema, sha256=td.digest(ROOT/path)) for path, schema in INPUTS.values()}
    model_path = Path(models_output)
    inputs[str(model_path.relative_to(ROOT)) if model_path.is_relative_to(ROOT) else str(model_path)] = dict(
        schema=model_collection['schema'], sha256=td.digest(model_path))
    result = dict(schema='terluna.research.array-flexible-benchmark/1',
        producer=dict(files={path: td.digest(ROOT/path) for path in FILES}, constants=constants_used(FILES)),
        inputs=inputs, baseline_commit='669ddd3', evidence=EVIDENCE,
        reading_rule='Fixed short numerical tests, not a coupled search. Preserve failures and invalid prefixes. Compare physical probes across mesh/mode bases, not modal coordinates. Positive sampled coverage is not continuous coverage.',
        fixed_reference=model_collection['reference_case'],
        case_plan=plan, cases=cases, refinements=aggregate_refinements(cases),
        solver_comparison=solver_comparison(cases),
        absolute_tolerance_pilot=tolerance_comparison(cases),
        resolution_assessment=resolution_assessment(cases),
        resources=dict(model_build=model_collection['resources'], one_three_tile=resource_comparison(cases),
            total_wall_s=time.perf_counter()-started, total_process_cpu_s=time.process_time()-cpu_started,
            peak_process_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024),
        limits=['normal-only linear prestress model with assumed finite support compliances',
            'endpoint-lumped support inertia and polygonal approximation to circular scallops',
            'short circular lunar point-mass benchmark, fixed source direction, no Earth/Sun gravity or eclipses',
            'zero assumed damping; no in-plane motion, wrinkling, thermal distortion or nonlinear constitutive response',
            'triangle-centroid photons with sampled solar directions and receiver rays; hardware shadows unmeshed',
            'instantaneous AABB separation lower bounds; no continuous clearance/contact certificate',
            'no solved large rolled service command, propulsion, storage, bus or equipment reaction budget',
            'no coupled parameter search or optimizer executed'])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(output), completed=sum(r['integration_reached_requested_end'] for r in cases),
        cases=len(cases), runtime_s=result['resources']['total_wall_s'])))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive-dir', type=Path, default=ARCHIVE)
    parser.add_argument('--output', type=Path, default=OUT)
    parser.add_argument('--models-output', type=Path, default=MODELS_OUT)
    args = parser.parse_args()
    run(archive_dir=args.archive_dir, output=args.output, models_output=args.models_output)
