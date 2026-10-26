"""Analytical transport solutions separate mass-grid spreading from time error."""
import json
import math

import numpy as np
import pytest

from atmosphere.snow import growth_transport as gt
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_changed
from shared.baselines import recorded_file


def test_constant_growth_characteristics_and_grid_poisson_distribution():
    m = np.arange(1., 81.)
    n = np.zeros_like(m); n[0] = 1.
    times = np.array([0., 2.])
    continuous = gt.characteristics([1.], [1.], times, lambda x: np.ones_like(x), 100.)
    assert continuous['masses'][-1, 0] == pytest.approx(3., rel=1e-12)
    discrete = gt.fixed_grid(m, n, times, np.ones_like(m), 100., max_step_s=.002)
    k = m-1
    exact = np.exp(-2. + k*np.log(2.) - np.array([math.lgamma(x+1) for x in k]))
    np.testing.assert_allclose(discrete['numbers'][-1], exact, rtol=0, atol=3e-7)
    assert np.dot(discrete['numbers'][-1], (m-3.)**2) == pytest.approx(2., rel=1e-12)
    assert discrete['upper_exit_mass'][-1] < 1e-50


def test_linear_growth_geometric_grid_has_known_spatial_variance():
    ratio, rate, duration = 1.1, .4, 1.
    m = ratio**np.arange(100)
    n = np.zeros_like(m); n[0] = 1.
    result = gt.fixed_grid(m, n, [0., duration], rate*m, 100., max_step_s=.001)
    final = result['numbers'][-1]
    mean = np.dot(m, final)
    variance = np.dot(final, (m-mean)**2)
    assert mean == pytest.approx(math.exp(rate*duration), rel=1e-7)
    assert variance/mean**2 == pytest.approx(math.exp(rate*(ratio-1)*duration)-1., rel=5e-6)
    reference = gt.characteristics([1.], [1.], [0., duration], lambda x: rate*x, 100.)
    assert reference['masses'][-1, 0] == pytest.approx(math.exp(rate*duration), rel=1e-10)


def test_fixed_grid_time_refinement_converges_without_removing_spatial_spread():
    m = np.arange(1., 65.)
    n = np.zeros_like(m); n[0] = 1.
    errors = []
    for step in (.04, .02, .01):
        result = gt.fixed_grid(m, n, [0., 1.], np.ones_like(m), 100., max_step_s=step)
        final = result['numbers'][-1]
        errors.append(abs(final[0]-math.exp(-1.)))
        assert np.dot(final, (m-2.)**2) == pytest.approx(1., rel=1e-12)
    assert errors[1] < .27*errors[0]
    assert errors[2] < .27*errors[1]


def test_supply_cap_preserves_mass_and_zero_weight_support_tracers():
    result = gt.characteristics([1., 2., 3.], [10., 0., 0.], [0., 2.], lambda x: np.ones_like(x), 2.)
    np.testing.assert_allclose(result['masses'][-1], [1.4, 2.4, 3.4], rtol=1e-12)
    assert result['deposited_mass'][-1] == pytest.approx(4., rel=1e-12)
    np.testing.assert_array_equal(result['numbers'], [10., 0., 0.])


def test_projection_conserves_number_and_mass_without_repeated_transport():
    pivots = np.array([1., 2., 4., 8.])
    masses, numbers = np.array([1.25, 3., 6.]), np.array([2., .1, .001])
    projected = gt.project_once(pivots, masses, numbers)
    assert projected.sum() == pytest.approx(numbers.sum(), rel=1e-14)
    assert np.dot(projected, pivots) == pytest.approx(np.dot(masses, numbers), rel=1e-14)
    with pytest.raises(ValueError, match='leaves'):
        gt.project_once(pivots, [9.], [1.])


def test_physical_seed_quadrature_preserves_known_number_and_mass():
    geometry = GEOMETRIES['bf95_dmax_dendrites_area']
    low, high = 1e-4, 4e-4
    mean = (geometry.mass_coefficient*(high**geometry.mass_exponent-low**geometry.mass_exponent)
            / (geometry.mass_exponent*math.log(high/low)))
    for size in (32, 64, 128):
        mass, number = gt.physical_seeds(geometry, size, 1e5)
        assert number.sum() == pytest.approx(1e5, rel=1e-14)
        assert np.dot(number, mass) == pytest.approx(1e5*mean, rel=1e-14)
        assert number[-2:].sum() == 0
        np.testing.assert_allclose(geometry.size(mass[-2:]), [low, high], rtol=1e-14)


def test_growth_product_binds_inputs_and_quantifies_time_reference_controls():
    if not gt.OUT.exists():
        pytest.skip('growth comparison has not been produced')
    product = json.loads(gt.OUT.read_text())
    assert product['schema'] == gt.SCHEMA
    for category in ('files', 'inputs'):
        for path, expected in product['producer'][category].items():
            assert gt.sha(recorded_file(path, expected)) == expected, path
    assert not constants_changed(product['producer']['constants'])
    assert all(row['passed'] for row in product['controls']['reference'])
    assert all(row['passed'] for row in product['controls']['time'].values())
    for run in product['runs']:
        n = np.asarray(run['final_numbers_m3'])
        assert np.isfinite(n).all() and (n >= 0).all()
        for row in run['snapshots']:
            ledger = row['water_ledger']
            source = ledger['initial']+ledger['deposited']
            residual = source-ledger['remaining']-ledger['upper_exit']
            assert abs(residual)/source < 1e-10
            assert abs(row['number_residual'])/row['number_m3'] < 1e-10
            assert 0 <= row['mass_share_beyond_once_projected_support'] <= 1
