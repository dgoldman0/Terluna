"""Check analytic-screen bookkeeping, boundaries and preserved provenance."""
from itertools import product
import json
import math

import numpy as np
import pytest

from shared.provenance import constants_changed
from . import geometry_screen as gs


@pytest.fixture(scope='module')
def context():
    return gs.context_from_inputs(gs.load_inputs())


@pytest.fixture(scope='module')
def mechanical(context):
    config = dict(id='check', side_m=10000., retreat_fraction=.05,
        depth_fraction=.01, compression_material_id='m55j_2500_ud',
        cable_material_id='kevlar49_yarn', normal_load_factor=1.)
    return gs.mechanical_case(config, context)


def test_deterministic_bounded_enumeration_separates_mechanics_and_spacing():
    a, b = list(gs.mechanical_configurations()), list(gs.mechanical_configurations())
    assert a == b
    assert len(a) == len({row['id'] for row in a}) == 2700
    assert len(a)*len(gs.SHINGLES)*len(gs.PITCH_FACTORS) == 16200
    actual = {tuple(row[name] for name in gs.GRID) for row in a}
    expected = set(product(*(gs.GRID[name] for name in gs.GRID)))
    assert actual == expected
    with pytest.raises(ValueError, match='nonempty'):
        list(gs.mechanical_configurations(dict(gs.GRID, side_m=[])))


def test_input_context_reads_physical_fields_without_rounded_reference_values(context):
    data = gs.load_inputs()
    old = data['baseline']['hardware']['variants']['annulus_reference']['design']
    assert context['design'] == old
    assert context['reference_side_m'] == old['side_m']
    assert context['lattice'] == data['baseline']['fleet_service']['ideal_service_lattice']
    assert context['population'] == data['baseline']['fleet_service']['population']
    assert context['normal_acceleration_m_s2'] == data['comparison']['loading']['resolved_normal_acceleration_m_s2']
    assert context['safety_clearance_m'] == 150


def test_mass_force_and_finite_material_envelope_close(mechanical, context):
    row = mechanical
    assert row['converged']
    assert row['support_mass_kg'] == pytest.approx(sum(row['components_kg'].values()))
    d = context['design']
    sigma = d['film_areal_mass_kg_m2']+d['substrate_areal_mass_kg_m2']+d['trim_area_fraction']*d['trim_hardware_areal_mass_kg_m2']
    assert row['retained_sheet_mass_kg'] == pytest.approx(sigma*row['retained_area_m2'])
    assert row['support_and_sheet_subtotal_kg'] == pytest.approx(row['support_mass_kg']+row['retained_sheet_mass_kg'])
    assert row['membrane_normal_force_N'] == pytest.approx(row['retained_area_m2']*row['membrane_normal_pressure_Pa'])
    assert row['membrane_normal_force_N'] == pytest.approx(4*row['membrane_normal_force_per_boom_N'])
    required = np.array([row['support_outside_span_m']/2]*2+[row['support_half_depth_m']])+row['deformation_allowance_m']
    assert np.all(np.asarray(row['acceleration_envelope_half_spans_m']) >= required*(1-1e-9))
    actual_acceleration = gs.box_normal_acceleration_envelope(required, context['omega'], context['alpha'], context['radial'], context['mean_motion_squared_s2'])
    actual_acceleration *= row['normal_load_factor']*d['membrane_dynamic_factor']
    assert actual_acceleration == pytest.approx(row['required_finite_envelope_acceleration_m_s2'])
    assert row['structural_acceleration_m_s2'] >= actual_acceleration*(1-1e-9)
    assert not row['complete_spacecraft_budget']
    assert not row['lower_load_attitude_demonstrated']


def test_load_factors_leave_solar_pressure_unchanged_and_do_not_prove_attitudes(context, mechanical):
    config = {name: mechanical[name] for name in gs.GRID}
    low = gs.mechanical_case(dict(config, id='low', normal_load_factor=.25), context)
    d = context['design']
    solar = d['membrane_dynamic_factor']*d['pressure_Pa']
    assert low['membrane_normal_pressure_Pa']-solar == pytest.approx(.25*(mechanical['membrane_normal_pressure_Pa']-solar))
    assert not low['lower_load_attitude_demonstrated']
    assert not low['validated_design']


def test_projected_inscribed_rectangle_and_inventory_normalization(mechanical, context):
    for angle, pitch_factor in product(gs.SHINGLES, gs.PITCH_FACTORS):
        row = gs.spacing_case(mechanical, context, angle, pitch_factor, 0)
        physical = row['physical_protected_rectangle_m']
        projected = row['projected_protected_rectangle_m']
        assert max(physical) <= row['protected_core_side_m']*(1+1e-12)
        assert physical[0]*math.cos(math.radians(angle)) == pytest.approx(projected[0])
        assert physical[1] == pytest.approx(projected[1])
        for axis in range(2):
            reference_width = context['lattice'][('optical_width_m', 'optical_height_m')[axis]]
            assert row['pitch_m'][axis]/projected[axis] == pytest.approx(context['lattice']['pitch_m'][axis]/reference_width)
        assert row['equivalent_population']*row['reference_gross_area_per_unit_m2'] == pytest.approx(context['population']*context['reference_side_m']**2)
        assert row['subtotal_kg_per_reference_gross_m2']*context['population']*context['reference_side_m']**2 == pytest.approx(row['support_and_sheet_subtotal_kg']*row['equivalent_population'])


def test_additional_overlap_costs_population_and_reduces_parallel_clearance(mechanical, context):
    nominal = gs.spacing_case(mechanical, context, 3., 1., 0)
    dense = gs.spacing_case(mechanical, context, 3., .9, 1)
    assert dense['equivalent_population']/nominal['equivalent_population'] == pytest.approx(1/.9**2)
    assert dense['subtotal_kg_per_reference_gross_m2']/nominal['subtotal_kg_per_reference_gross_m2'] == pytest.approx(1/.9**2)
    assert dense['nominal_parallel_envelope_clearance_m'] < nominal['nominal_parallel_envelope_clearance_m']
    for row in (nominal, dense):
        assert row['filters']['parallel_150m_clearance'] == (row['nominal_parallel_envelope_clearance_m'] >= 150)
    assert nominal['nominal_parallel_envelope_clearance_m'] == pytest.approx(
        nominal['pitch_m'][0]*math.sin(math.radians(3.))-2*(mechanical['support_half_depth_m']+mechanical['deformation_allowance_m']))


def test_bus_scenarios_charge_inventory_but_do_not_invent_power_hardware(mechanical, context):
    row = gs.spacing_case(mechanical, context, 3., 1., 0)
    for bus in row['bus_scenarios'].values():
        assert bus['fleet_bus_inventory_kg'] == pytest.approx(bus['assumed_mass_per_unit_kg']*row['equivalent_population'])
        assert bus['fleet_housekeeping_W'] == pytest.approx(bus['assumed_housekeeping_per_unit_W']*row['equivalent_population'])
        assert bus['subtotal_with_bus_kg_per_reference_gross_m2'] == pytest.approx(row['subtotal_kg_per_reference_gross_m2']+bus['assumed_mass_per_unit_kg']/row['reference_gross_area_per_unit_m2'])
        assert not bus['bus_requirements_demonstrated']
        assert not bus['power_system_mass_included']


def test_failures_survive_all_spacing_choices_and_archive(monkeypatch, context, tmp_path):
    def failed(*args, **kwargs):
        raise gs.SupportInfeasible('injected failure')
    monkeypatch.setattr(gs, 'scalloped_support', failed)
    mechanical = gs.mechanical_case(next(gs.mechanical_configurations()), context)
    rows = [gs.spacing_case(mechanical, context, angle, pitch, i)
        for i, (angle, pitch) in enumerate(product(gs.SHINGLES, gs.PITCH_FACTORS))]
    record = gs.write_jsonl(tmp_path/'failed.jsonl', rows)
    reloaded = [json.loads(line) for line in (tmp_path/'failed.jsonl').read_text().splitlines()]
    assert len(reloaded) == record['rows'] == 6
    assert record['sha256'] == gs.td.digest(tmp_path/'failed.jsonl')
    for row in reloaded:
        assert not row['converged']
        assert row['failure'] == 'injected failure'
        assert 'subtotal_kg_per_reference_gross_m2' not in row
        assert not row['complete_spacecraft_budget']


def test_pareto_matches_exhaustive_dominance_and_preserves_equal_ties():
    coordinates = [(1, 3, 2), (2, 2, 2), (3, 1, 2), (2, 3, 3),
                   (1, 3, 2), (4, 4, 4), (.5, 5, 1)]
    rows = [dict(id=str(i), subtotal_kg_per_reference_gross_m2=p[0],
        equivalent_population=p[1], nominal_parallel_envelope_clearance_m=-p[2])
        for i, p in enumerate(coordinates)]
    expected = {str(i) for i, p in enumerate(coordinates) if not any(
        all(a <= b for a, b in zip(q, p)) and any(a < b for a, b in zip(q, p))
        for q in coordinates)}
    assert set(gs.pareto_ids(rows)) == expected
    assert set(gs.pareto_ids(list(reversed(rows)))) == expected
    assert {'0', '4'} <= expected


def test_summary_retains_failed_strata_and_separates_load_hypotheses(mechanical, context):
    rows = [gs.spacing_case(mechanical, context, angle, pitch, i)
        for i, (angle, pitch) in enumerate(product(gs.SHINGLES, gs.PITCH_FACTORS))]
    lower = [dict(row, id='lower_'+row['id'], normal_load_factor=.25) for row in rows]
    summary = gs.summarize([mechanical], rows+lower, context)
    assert summary['configuration_count'] == 12
    assert {stratum['normal_load_factor'] for stratum in summary['pareto_strata']} == {.25, 1.}
    for stratum in summary['pareto_strata']:
        assert stratum['eligible_count'] == stratum['pareto_count']+stratum['dominated_count']
        assert all(identifier.startswith('lower_') == (stratum['normal_load_factor'] == .25)
                   for identifier in stratum['pareto_ids'])


def test_committed_summary_binds_exact_producer_inputs_and_configuration_counts():
    if not gs.OUT.exists():
        pytest.skip('Run the bounded screen to produce its summary')
    result = json.loads(gs.OUT.read_text())
    assert result['schema'] == 'terluna.research.array-geometry-screen/1'
    assert set(result['producer']['files']) == set(gs.FILES)
    assert set(result['inputs']) == {path for path, _ in gs.INPUTS.values()}
    for path, digest in result['producer']['files'].items():
        assert gs.td.digest(gs.ROOT/path) == digest
    for path, schema in gs.INPUTS.values():
        assert result['inputs'][path] == dict(schema=schema, sha256=gs.td.digest(gs.ROOT/path))
    assert not constants_changed(result['producer']['constants'])
    assert result['summary']['mechanical_count'] == 2700
    assert result['summary']['configuration_count'] == 16200
    assert result['archives']['mechanical']['rows'] == 2700
    assert result['archives']['configurations']['rows'] == 16200
    assert sum(row['bytes'] for row in result['archives'].values()) < 100_000_000
    for archive in result['archives'].values():
        path = gs.ROOT/archive['path']
        if path.exists():
            assert gs.td.digest(path) == archive['sha256']
