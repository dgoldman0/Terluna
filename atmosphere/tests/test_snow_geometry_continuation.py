"""The geometry control preserves mass pivots, seed forcing and frozen numerics."""
import json

import numpy as np
import pytest

from atmosphere.snow import flakes as sf, geometry_continuation as gc
from atmosphere.snow.column_experiments import solve_column
from atmosphere.snow.sectional import seed_distribution


def test_control_preserves_baseline_mass_pivots_and_seed_distribution():
    masses = gc.baseline_mass_grid()
    baseline = np.geomspace(float(gc.BASE.mass(2e-5)), float(gc.BASE.mass(10.)), 96)
    np.testing.assert_array_equal(masses, baseline)
    np.testing.assert_array_equal(seed_distribution(masses, gc.GEOMETRY),
                                  seed_distribution(masses, gc.BASE))
    assert gc.GEOMETRY.size(masses[-1]) == pytest.approx(7.079457843841377)
    assert gc.GEOMETRY.size(masses[-1]) < gc.BASE.size(masses[-1])


def test_continuation_is_continuous_and_retains_small_particle_physics():
    join = gc.GEOMETRY.anchor_m
    below = np.geomspace(2e-5, join, 40)
    np.testing.assert_array_equal(gc.GEOMETRY.mass(below), gc.BASE.mass(below))
    np.testing.assert_array_equal(gc.GEOMETRY.area(below), gc.BASE.area(below))
    masses = gc.GEOMETRY.mass(np.array([join, 2*join, 4*join]))
    np.testing.assert_allclose(masses / masses[0], [1., 4., 16.], rtol=1e-15)
    np.testing.assert_allclose(gc.GEOMETRY.size(masses), [join, 2*join, 4*join], rtol=1e-15)


def test_source_velocity_and_duration_remain_the_baseline_values():
    world, _ = sf.moon_column()
    args = gc.column_arguments(world)
    from atmosphere.snow import aggregation as ag
    seed_d = np.sqrt(ag.START_M[0]*ag.START_M[1])
    state = world.state(world.height_of(sf.T0 + args['top_c']))
    expected = float(ag.speeds(np.array([seed_d]), np.array([gc.BASE.mass(seed_d)]),
                              world.gravity, state, gc.BASE.area_ratio)[0][0])
    assert args['injection_velocity_m_s'] == expected
    duration = (args['fill'] + args['gather'])*(world.height_of(sf.T0+args['top_c'])-world.z_freeze)/expected
    assert duration / 3600. == pytest.approx(132.5892799336426, rel=1e-9)
    assert args['maximum_diameter_m'] == float(gc.GEOMETRY.size(args['mass_grid'][-1]))


def test_small_column_is_identical_when_mass_domain_is_below_join():
    masses = np.geomspace(float(gc.BASE.mass(2e-5)), float(gc.BASE.mass(.003)), 24)
    args = dict(layers=2, mass_grid=masses, mass_bins=24, maximum_diameter_m=.003,
                fill=.01, gather=.01, max_step_s=13., growth_scheme='limited')
    base = solve_column(sf.earth_column(), -2., 1e4, 1., geometry=gc.BASE, **args)
    control = solve_column(sf.earth_column(), -2., 1e4, 1., geometry=gc.GEOMETRY, **args)
    for key in ('arriving_numbers_m2', 'final_number_m3', 'time_mean_number_m3'):
        np.testing.assert_array_equal(base[key], control[key])
    assert base['mass_ledger_kg_m2'] == control['mass_ledger_kg_m2']
    assert base['diagnostics'] == control['diagnostics']
    exported = gc.export(control, {}, .1)
    assert not exported['optical_prediction_eligible']
    np.testing.assert_array_equal(exported['native_arrays']['final_number_m3'], control['final_number_m3'])
    json.dumps(exported, allow_nan=False)


def test_identity_change_during_calculation_rejects_product(monkeypatch):
    identities = iter(({'files': {'driver': 'before'}}, {'files': {'driver': 'after'}}))
    monkeypatch.setattr(gc, 'identity', lambda binding: next(identities))
    monkeypatch.setattr(gc, 'solve_column', lambda *args, **kwargs: {})
    with pytest.raises(RuntimeError, match='Producer changed during run'):
        gc.run()


def test_completed_product_is_bound_and_native_arrays_close_budget():
    if not gc.OUT.exists():
        pytest.skip('The full geometry sensitivity has not completed')
    result = json.loads(gc.OUT.read_text())
    _, binding = sf.moon_column()
    assert result['producer'] == gc.identity(binding)
    assert result['specification'] == gc.SPECIFICATION
    native = result['native_arrays']
    masses = np.asarray(native['masses_kg'])
    np.testing.assert_array_equal(masses, gc.baseline_mass_grid())
    np.testing.assert_array_equal(native['diameters_m'], gc.GEOMETRY.size(masses))
    numbers = np.asarray(native['final_number_m3'])
    dz = np.diff(native['layer_edges_m'])
    remaining = np.sum(numbers*masses*dz[:, None])
    budget = result['mass_ledger_kg_m2']
    assert remaining == pytest.approx(budget['remaining'], rel=1e-12)
    arrivals = np.asarray(native['arriving_numbers_m2'])
    flux = float(np.dot(masses, arrivals))*3600./result['gathering_seconds']
    assert flux == pytest.approx(result['summary']['snowfall_mm_h'], rel=1e-12)
    assert abs(budget['relative_residual']) < 1e-10
    assert result['domain_adequacy']['passed'] == (
        result['domain_adequacy']['mass_fraction_exited'] <= result['domain_adequacy']['tolerance'])
    assert not result['optical_prediction_eligible']
