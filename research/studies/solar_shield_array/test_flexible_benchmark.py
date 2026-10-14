"""Integration, refinement comparability and provenance of the fixed benchmark."""
from copy import deepcopy
import json
import math

import numpy as np
import pytest

from protection.dynamics.flexible_fleet import Body
from shared.provenance import constants_changed
from . import flexible_benchmark as fb


@pytest.fixture(scope='module')
def reference():
    return fb.select_reference(fb.load_inputs())


@pytest.fixture(scope='module')
def model_product(reference):
    model = fb.model_from_support(reference['support'], reference['sheet_areal_mass_kg_m2'],
        mesh_divisions=4, retained_modes=3, hub_equipment_mass_kg=reference['hub_equipment_mass_kg'])
    producer = dict(files={path: fb.td.digest(fb.ROOT/path) for path in fb.MODEL_FILES}, constants={})
    inputs = {path: dict(sha256=fb.td.digest(fb.ROOT/path)) for path, _ in fb.INPUTS.values()}
    return fb.export_model(model, reference, producer, inputs)


def test_fixed_case_plan_has_matched_refinements_and_no_design_search():
    cases = fb.fixed_cases()
    by_id = {row['id']: row for row in cases}
    assert len(cases) == len(by_id) == 24
    for key, names in [('mesh_divisions', ['reference', 'mesh12', 'mesh16']),
                       ('retained_modes', ['mesh12', 'modes12', 'modes20']),
                       ('max_step_s', ['reference', 'step5', 'step2_5']),
                       ('source_rings', ['reference', 'source7', 'source19'])]:
        records = [{k: v for k, v in by_id[name].items() if k not in ('id', key)} for name in names]
        assert records[0] == records[1] == records[2]
    one, three = by_id['cost1'], by_id['cost3']
    assert {k: v for k, v in one.items() if k not in ('id', 'tiles')} == {k: v for k, v in three.items() if k not in ('id', 'tiles')}
    for name in ('solver_rk45', 'solver_radau'):
        assert {k: v for k, v in by_id[name].items() if k not in ('id', 'method')} == {k: v for k, v in one.items() if k not in ('id', 'method')}
    full, refined = by_id['full_mesh16'], by_id['full_mesh16_step']
    assert {k: v for k, v in full.items() if k not in ('id', 'max_step_s')} == {k: v for k, v in refined.items() if k not in ('id', 'max_step_s')}


def test_complete_basis_cases_remove_exactly_the_three_rigid_normal_modes():
    by_id = {row['id']: row for row in fb.fixed_cases()}
    for name in ('full_mesh8', 'modes170', 'full_mesh16', 'full_mesh16_step'):
        row = by_id[name]
        node_count = (row['mesh_divisions']+1)**2+4
        assert row['retained_modes'] == node_count-3
    assert by_id['modes40']['mesh_divisions'] == by_id['modes80']['mesh_divisions'] == 12
    assert by_id['full_mesh16_step']['max_step_s'] == .5


def test_resolution_assessment_distinguishes_mode_mesh_and_time_errors():
    def record(identifier, displacement):
        observable = dict(last_position_m=[[0., 0., 0.]], last_velocity_m_s=[[0., 0., 0.]],
            last_corner_positions_m=[[[0., 0., 0.]]],
            last_shape_probe_displacements_m=[[[0., 0., displacement]]],
            last_total_photon_force_N=[0., 0., 0.])
        return dict(id=identifier, integration_reached_requested_end=True, observables=observable)
    values = {'mesh12': .01, 'modes12': .02, 'modes20': .05, 'modes40': .0995,
        'modes80': .1001, 'modes170': .1, 'full_mesh8': .08,
        'full_mesh16': .103, 'full_mesh16_step': .103000001}
    assessment = fb.resolution_assessment([record(key, value) for key, value in values.items()])
    modal = {row['first']: row for row in assessment['fixed_mesh12_modal_errors_against_complete_basis']}
    assert not modal['modes20']['meets_declared_probe_tolerance']
    assert modal['modes40']['meets_declared_probe_tolerance']
    assert modal['modes170']['maximum_shape_probe_difference_m'] == 0
    assert not assessment['final_complete_basis_mesh_pair_meets_declared_probe_tolerance']
    assert assessment['complete_basis_time_pair_meets_declared_probe_tolerance']


def test_reference_is_fixed_prior_support_and_mass_matched_optical_row(reference):
    inputs = fb.load_inputs()
    prior = next(row for row in inputs['support']['command_loaded_supports']
                 if row['id'] == reference['support_case_id'])
    assert prior['support'] == reference['support']
    assert prior['compression_material_id'] == 'm55j_2500_ud'
    assert prior['cable_material_id'] == 'kevlar49_yarn'
    assert reference['optical_row_description']['areal_mass_g_m2']/1000 == pytest.approx(reference['design']['film_areal_mass_kg_m2'])
    assert sum(reference['optical_fractions'].values()) == pytest.approx(1.)


def test_export_closes_mass_centres_geometry_and_rejects_tampering(model_product, reference):
    body = Body.from_product(model_product)
    expected = reference['support']['mass_kg']+reference['support']['geometry']['retained_area_m2']*reference['sheet_areal_mass_kg_m2']+reference['hub_equipment_mass_kg']
    assert body.mass_kg == pytest.approx(expected)
    assert np.linalg.norm(np.average(body.nodes_m, axis=0, weights=body.nodal_masses_kg)) < 1e-10
    assert model_product['model_sha256'] == fb.canonical_hash(model_product['model'])
    damaged = deepcopy(model_product)
    damaged['model']['nodal_masses_kg'][0] *= 1.01
    with pytest.raises(ValueError, match='hash mismatch'):
        Body.from_product(damaged)
    metadata = model_product['model']['metadata']
    assert metadata['nodal_inertia_relative_frobenius_error'] > 0
    assert metadata['area_relative_error'] > 0
    assert not metadata['whole_spacecraft_budget_closed']


def test_common_physical_probes_reproduce_affine_fields(model_product, reference):
    probes, indices, weights = fb.probe_interpolation(model_product['model'], reference['design']['side_m'])
    nodes = np.asarray(model_product['model']['nodes_m'])
    field = 2+3*nodes[:, 0]-.2*nodes[:, 1]
    reconstructed = np.einsum('pi,pi->p', weights, field[indices])
    assert reconstructed == pytest.approx(2+3*probes[:, 0]-.2*probes[:, 1], abs=1e-10)
    assert np.allclose(weights.sum(axis=1), 1.)
    assert (weights >= -1e-10).all()


def test_initial_bump_is_elastic_projection_with_no_rigid_mass_motion(model_product, reference):
    body = Body.from_product(model_product)
    q = fb.initial_modal_bump(model_product['model'], reference['design']['side_m'])
    displacement = np.einsum('nik,k->ni', body.modes, q)
    assert np.linalg.norm(q) > 0
    assert np.linalg.norm(np.sum(body.nodal_masses_kg[:, None]*displacement, axis=0)) < 1e-7
    assert np.linalg.norm(np.sum(np.cross(body.nodes_m, body.nodal_masses_kg[:, None]*displacement), axis=0)) < 1e-3


def test_period_case_uses_physical_frequencies_and_explicit_cost_cap(model_product):
    config = next(row for row in fb.fixed_cases() if row['id'] == 'lowest_mode_period')
    fast = dict(resources=dict(integration_wall_s=1., rhs_evaluations=1000))
    full = fb.period_configuration(config, model_product, fast)
    frequencies = np.sqrt(model_product['model']['omega_squared_s2'])
    assert full['duration_s'] == pytest.approx(2*math.pi/frequencies.min())
    assert full['max_step_s'] <= full['duration_s']/100
    assert full['max_step_s'] <= 2*math.pi/frequencies.max()/12
    slow = dict(resources=dict(integration_wall_s=1000., rhs_evaluations=1))
    limited = fb.period_configuration(config, model_product, slow)
    assert limited['period_cost_cap_applied']
    assert limited['duration_s'] == min(1000., limited['lowest_mode_period_s'])
    assert limited['requested_period_fraction'] <= 1.


def test_dimensioned_absolute_tolerances_follow_physical_mode_amplitudes(model_product):
    body = Body.from_product(model_product)
    values = fb.state_absolute_tolerances([body])
    assert len(values) == body.state_size
    assert np.isfinite(values).all() and (values > 0).all()
    mode_max = np.max(np.linalg.norm(body.modes, axis=1), axis=0)
    assert values[13:13+body.k]*mode_max == pytest.approx(np.full(body.k, 1e-7))
    assert values[13+body.k:13+2*body.k]*mode_max == pytest.approx(np.full(body.k, 1e-9))
    assert values[16+2*body.k] == 1e-3
    assert values[17+2*body.k:] == pytest.approx(np.full(3, .01))


def test_wall_budget_reports_only_an_attempted_stage(model_product, reference, tmp_path):
    case = dict(fb.fixed_cases()[0], integration_wall_budget_s=1e-30)
    with pytest.raises(fb.BenchmarkBudgetExceeded) as error:
        fb.execute_case(case, model_product, reference, tmp_path)
    assert 'not an accepted' in error.value.context['state_interpretation']
    assert error.value.context['elapsed_wall_s'] > case['integration_wall_budget_s']
    assert not list(tmp_path.glob('*.npz'))


def test_nonfinite_diagnostics_are_disclosed_instead_of_written_as_nan():
    special = []
    value = fb.finite_json(dict(a=np.array([1., np.inf]), b=float('nan')), special=special)
    assert value == dict(a=[1., None], b=None)
    assert {record['path'] for record in special} == {'a[1]', 'b'}
    json.dumps(value, allow_nan=False)


def test_short_driver_records_all_output_samples_controls_coverage_and_states(model_product, reference, tmp_path):
    case = dict(fb.fixed_cases()[0], id='test_short', mesh_divisions=4, retained_modes=3,
                duration_s=.1, output_step_s=.05, max_step_s=.05)
    result = fb.execute_case(case, model_product, reference, tmp_path)
    assert result['integration_reached_requested_end']
    assert result['stop_reason'] == 'completed'
    assert [row['time_s'] for row in result['sampled_history']] == pytest.approx([0., .05, .1])
    assert not result['accepted_operating_design']
    assert not result['continuous_coverage_demonstrated']
    assert not result['complete_spacecraft_budget']
    assert not result['optimized']
    assert result['resources']['rhs_evaluations'] > 0
    for row in result['sampled_history']:
        assert 0 <= row['minimum_receiver_source_coverage_fraction'] <= 1
        assert 0 <= row['maximum_receiver_bolometric_transmission'] <= 1
        assert row['maximum_trim_fraction_of_available_range'] == 0
        assert row['allocation_success_by_tile'] == [None]
    archive = result['state_archive']
    path = tmp_path/'test_short.npz'
    assert fb.td.digest(path) == archive['sha256']
    with np.load(path) as data:
        assert data['states'].shape[1] == len(result['sampled_history'])
        assert np.isfinite(data['states']).all()
        assert data['times_s'][-1] == .1


def test_preserved_product_hashes_and_all_fixed_cases():
    if not fb.OUT.exists():
        pytest.skip('Run the fixed benchmark to generate its product')
    result = json.loads(fb.OUT.read_text())
    assert result['schema'] == 'terluna.research.array-flexible-benchmark/1'
    assert set(result['producer']['files']) == set(fb.FILES)
    for path, digest in result['producer']['files'].items():
        assert fb.td.digest(fb.ROOT/path) == digest
    for path, record in result['inputs'].items():
        assert fb.td.digest(fb.ROOT/path) == record['sha256']
        assert json.loads((fb.ROOT/path).read_text())['schema'] == record['schema']
    assert not constants_changed(result['producer']['constants'])
    assert {row['id'] for row in result['cases']} == {row['id'] for row in fb.fixed_cases()}
    for case in result['cases']:
        assert not case['accepted_operating_design']
        if not case['integration_success']:
            assert case.get('failure') or case.get('solver_message')
    models = json.loads(fb.MODELS_OUT.read_text())
    for product in models['models'].values():
        assert fb.canonical_hash(product['model']) == product['model_sha256']
        for path, digest in product['producer_hashes'].items():
            assert fb.td.digest(fb.ROOT/path) == digest
        for path, digest in product['input_hashes'].items():
            assert fb.td.digest(fb.ROOT/path) == digest
    for case in result['cases']:
        assert case.get('model_sha256', models['models'][fb.model_key(case['configuration'])]['model_sha256']) == models['models'][fb.model_key(case['configuration'])]['model_sha256']
    assert {row['method'] for row in result['solver_comparison']['cases']} == {'DOP853', 'RK45', 'Radau'}
