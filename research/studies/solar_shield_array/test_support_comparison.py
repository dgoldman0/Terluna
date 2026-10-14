"""Integration and evidence checks for the support/material comparison.

Domain tests cover the analytic formulas. These checks catch stale inputs,
lost failure cases, changed area denominators and accidental promotion of a
component calculation into a complete spacecraft claim.
"""
from itertools import product as cartesian_product
import json
import math

import numpy as np
import pytest

from shared.provenance import constants_changed
from . import support_comparison as sc


@pytest.fixture(scope='module')
def records():
    output = json.loads(sc.OUT.read_text())
    inputs = {name: json.loads((sc.ROOT/path).read_text())
              for name, (path, _) in sc.INPUTS.items()}
    return output, inputs['baseline'], inputs['materials']


def support_rows(output):
    return (output['fixed_prestress']['cases']+
            output['command_loaded_supports']+
            output['equipment_load_sensitivity'])


def test_product_binds_exact_sources_inputs_and_constants(records):
    output, _, _ = records
    assert output['schema'] == 'terluna.research.array-support-comparison/1'
    assert set(output['producer']['files']) == set(sc.FILES)
    assert set(output['inputs']) == {path for path, _ in sc.INPUTS.values()}
    for path, digest in output['producer']['files'].items():
        assert sc.td.digest(sc.ROOT/path) == digest, path
    for path, schema in sc.INPUTS.values():
        assert output['inputs'][path]['schema'] == schema
        assert json.loads((sc.ROOT/path).read_text())['schema'] == schema
        assert output['inputs'][path]['sha256'] == sc.td.digest(sc.ROOT/path)
    assert not constants_changed(output['producer']['constants'])


def test_material_roles_do_not_promote_fibers_to_compression(records):
    _, _, catalogue = records
    sources = {source['id']: source for source in catalogue['sources']}
    seen = set()
    for name, role in [('compression_candidates', 'compression'),
                       ('tension_candidates', 'tension')]:
        for entry in catalogue[name]:
            assert entry['id'] not in seen
            seen.add(entry['id'])
            assert entry['role'] == role
            assert entry['evidence'] and entry['properties_basis']
            assert not entry['qualified_for_array_service']
            for field in ('density_kg_m3', 'youngs_modulus_Pa', 'admissible_stress_Pa'):
                assert math.isfinite(entry[field]) and entry[field] > 0
            for identifier in entry['source_ids']:
                assert sources[identifier]['url'].startswith('https://')
                assert sources[identifier]['read_extent']
            if entry['experimental']:
                assert role == 'tension'
    by_id = {entry['id']: entry for entry in catalogue['materials']}
    assert by_id['silicon_carbide_grade_f']['reported_tensile_strength_Pa'] is None
    assert by_id['silicon_carbide_grade_f']['reported_compressive_strength_relation'] == '>'
    assert by_id['m55j_fiber']['youngs_modulus_Pa'] > by_id['m55j_2500_ud']['youngs_modulus_Pa']


def test_generic_material_reproduces_preserved_complete_ledger(records):
    output, baseline, catalogue = records
    rows = {(row['kind'], row['material_id']): row
            for row in output['material_substitutions']}
    expected = set(cartesian_product(('annulus', 'window'),
                   (entry['id'] for entry in catalogue['compression_candidates'])))
    assert set(rows) == expected
    for kind in ('annulus', 'window'):
        old = baseline['hardware']['variants'][kind+'_reference']
        new = rows[kind, 'generic_composite_baseline']['budget']
        assert new['components_kg'] == pytest.approx(old['components_kg'], rel=1e-10)
        assert new['locked_inertia_kg_m2'] == pytest.approx(np.array(old['locked_inertia_kg_m2']), rel=1e-10)
        assert new['target_mass_kg'] == pytest.approx(old['target_mass_kg'])


def test_complete_ledgers_close_mass_and_rigid_requirements_without_certification(records):
    output, _, _ = records
    ledgers = [row['budget'] for row in output['material_substitutions']]
    ledgers.extend(output[name] for name in (
        'refined_load_perimeter', 'refined_load_high_modulus_perimeter'))
    for budget in ledgers:
        d = budget['design']
        assert sum(budget['components_kg'].values()) == pytest.approx(budget['total_mass_kg'])
        assert budget['total_mass_kg']/d['side_m']**2 == pytest.approx(budget['total_areal_mass_kg_m2'])
        assert budget['target_mass_kg'] == pytest.approx(d['target_total_areal_mass_kg_m2']*d['side_m']**2)
        assert np.linalg.eigvalsh(budget['locked_inertia_kg_m2']).min() > 0
        iteration = budget['rigid_sizing_iteration']
        if iteration['converged']:
            limit = (1+iteration['tolerance_fraction'])/iteration['requirement_margin']
            assert budget['momentum_capacity_utilization_in_screen'] <= limit
            assert budget['motor_torque_utilization_in_screen'] <= limit
        assert not iteration['periodic_controller_demonstrated']
        assert not budget['validated_design']
        assert not budget['optical_trim']['qualified_authority']


def test_every_planned_support_case_is_preserved_including_failures(records):
    output, _, catalogue = records
    fixed_ids = {row['id'] for row in output['fixed_prestress']['cases']}
    expected_fixed = {f'geometry_{fraction:g}' for fraction in (.005, .05, .1, .2)}
    expected_fixed |= {'cable_'+entry['id'] for entry in catalogue['tension_candidates']}
    assert fixed_ids == expected_fixed
    loaded = output['command_loaded_supports']
    assert len(loaded) == 18
    assert {(row['side_m'], row['edge_retreat_fraction'], row['boom_depth_m']/row['side_m'])
            for row in loaded} == set(cartesian_product((1000., 3000., 10000.),
                                                       (.005, .05, .1), (.005, .01)))
    for row in support_rows(output):
        assert not row['complete_spacecraft_budget']
        assert not row['validated_design']
        if not row['converged']:
            assert row['failure']
            assert 'support_and_sheet_subtotal_kg' not in row
        else:
            assert row['omitted']


def test_nonconvergence_is_a_recorded_failure_not_a_synthetic_mass(monkeypatch, records):
    _, baseline, catalogue = records
    def reject(*args, **kwargs):
        raise sc.SupportInfeasible('injected self-load nonconvergence')
    monkeypatch.setattr(sc, 'scalloped_support', reject)
    result = sc.support_case(sc.inherited_design(baseline),
        catalogue['compression_candidates'][0], catalogue['tension_candidates'][0],
        .05, baseline['fleet_service']['population'], tension=1., case_id='nonconverged')
    assert not result['converged']
    assert result['id'] == 'nonconverged'
    assert result['failure'] == 'injected self-load nonconvergence'
    assert 'support_and_sheet_subtotal_kg' not in result
    assert not result['validated_design']


def test_equal_protected_area_normalization_charges_population_and_fixed_equipment(records):
    output, baseline, _ = records
    population = baseline['fleet_service']['population']
    lattice = baseline['fleet_service']['ideal_service_lattice']
    d = baseline['hardware']['variants']['annulus_reference']['design']
    for row in support_rows(output):
        coverage = row['coverage']
        scale = coverage['pitch_scale']
        assert scale > 0
        # Both pitches preserve the same overlap ratio, while the guaranteed
        # optical rectangle fits inside the loaded membrane's central square.
        assert coverage['pitch_x_m']/lattice['pitch_m'][0] == pytest.approx(scale)
        assert coverage['pitch_y_m']/lattice['pitch_m'][1] == pytest.approx(scale)
        assert max(lattice['optical_width_m'], lattice['optical_height_m'])*scale <= coverage['guaranteed_loaded_core_side_m']*(1+1e-12)
        assert coverage['equivalent_population']*coverage['reference_gross_area_per_unit_m2'] == pytest.approx(population*1e8)
        if row['converged']:
            assert row['subtotal_kg_per_reference_gross_m2']*population*1e8 == pytest.approx(
                row['support_and_sheet_subtotal_kg']*coverage['equivalent_population'])
            assert row['fixed_avionics_inventory_kg'] == pytest.approx(coverage['equivalent_population']*d['avionics_and_sensors_kg'])
            assert row['fixed_housekeeping_inventory_W'] == pytest.approx(coverage['equivalent_population']*d['housekeeping_W'])


def test_support_subtotals_use_retained_area_and_transfer_applied_normal_load(records):
    output, baseline, _ = records
    d = baseline['hardware']['variants']['annulus_reference']['design']
    sigma = d['film_areal_mass_kg_m2']+d['substrate_areal_mass_kg_m2']+d['trim_area_fraction']*d['trim_hardware_areal_mass_kg_m2']
    for row in support_rows(output):
        if not row['converged']:
            continue
        support = row['support']
        retained = support['geometry']['retained_area_m2']
        assert 0 < retained < row['side_m']**2
        assert row['retained_sheet_mass_kg'] == pytest.approx(retained*sigma)
        assert support['mass_kg'] == pytest.approx(sum(support['components_kg'].values()))
        assert row['support_and_sheet_subtotal_kg'] == pytest.approx(support['mass_kg']+retained*sigma)
        membrane_force = retained*row['pressure_for_support_Pa']
        assert support['membrane_normal_resultant_N'] == pytest.approx(membrane_force)
        assert 4*support['membrane_normal_force_per_boom_N'] == pytest.approx(membrane_force)
        boom = support['radial_boom']
        assert boom['root_transverse_force_N'] == pytest.approx(
            membrane_force/4+row['equipment_tip_force_N']+
            row['structural_acceleration_m_s2']*boom['mass_kg'], rel=2e-9)
        tensor = np.asarray(row['support_and_sheet_locked_inertia_kg_m2'])
        assert tensor == pytest.approx(tensor.T)
        assert np.linalg.eigvalsh(tensor).min() > 0
        assert tensor == pytest.approx(np.asarray(support['inertia_kg_m2'])+
            retained*sigma*np.asarray(support['geometry']['membrane_area_inertia_per_mass_m2']))


def test_spin_subtotals_charge_coverage_and_keep_diagnostic_evidence(records):
    output, baseline, _ = records
    spins = output['spinning_membranes']
    grid = spins['planar_grid']
    count = grid['count_multiplier_at_fixed_protected_area']
    assert not grid['original_grid_fully_covered']
    assert count > 1
    by_law = {}
    for case in spins['cases']:
        result = case['result']
        by_law[case['law'], case['spin_rad_s']] = result
        mass = result['ideal_film_plus_ballast_mass_kg']
        assert mass == pytest.approx(result['film_mass_kg']+result['ideal_rim_ballast_mass_kg'])
        assert case['subtotal_kg_per_reference_gross_m2'] == pytest.approx(mass*count/1e8)
        assert not case['complete_spacecraft_budget']
        assert case['equivalent_population'] == pytest.approx(
            baseline['fleet_service']['population']*count)
        assert not result['accepted_complete_design']
        assert result['unresolved']
        assert result['exact_axisymmetric_static_sag_m'] <= result['constant_tension_static_sag_bound_m']*(1+1e-12)
        axial_inertia = result['locked_axisymmetric_inertia_kg_m2'][2][2]
        min_spin = result['minimum_actual_axial_spin_rad_s']
        max_spin = result['maximum_actual_axial_spin_rad_s']
        assert 0 < min_spin <= max_spin
        assert result['spin_angular_momentum_N_m_s'] == pytest.approx(axial_inertia*max_spin)
        assert result['minimum_spin_angular_momentum_N_m_s'] == pytest.approx(axial_inertia*min_spin)
        assert result['spin_kinetic_energy_J'] == pytest.approx(.5*axial_inertia*max_spin**2)
    assert len(by_law) == 10
    for speed in (.0001, .0003, .001, .003, .01):
        fixed = by_law['fixed_normal_pressure_benchmark', speed]
        rolled = by_law['rolled_command', speed]
        assert fixed['peak_spin_precession_torque_N_m'] == 0
        assert rolled['peak_spin_precession_torque_N_m'] > 0
        assert rolled['peak_effective_pressure_Pa'] > fixed['peak_effective_pressure_Pa']


def test_time_sample_refinement_supports_the_reported_normal_load(records):
    output, _, _ = records
    refinement = output['loading']['refinement']
    assert [r['samples'] for r in refinement] == sorted(r['samples'] for r in refinement)
    values = [r['maximum_normal_acceleration_m_s2'] for r in refinement]
    assert min(values) > 0
    assert abs(values[-1]/values[-2]-1) < .005


def test_finite_tube_radius_removes_apparent_marginal_clearance_pass(records):
    output, _, _ = records
    case = next(row for row in output['command_loaded_supports']
                if row['id'] == 'loaded_10000_0.05_0.005')
    separation = case['coverage']['pitch_x_m']*math.sin(math.radians(1.5))
    # Centreline-depth accounting would incorrectly pass this case. The
    # finite tubes occupy four additional metres across the two surfaces.
    assert separation-case['boom_depth_m'] > 150
    assert case['nominal_parallel_support_clearance_m'] < 150
    assert not case['inherited_150m_support_gate_passes']
    assert case['nominal_parallel_support_clearance_m'] == pytest.approx(
        separation-case['boom_depth_m']-2*case['support']['maximum_member_outer_radius_m'])
