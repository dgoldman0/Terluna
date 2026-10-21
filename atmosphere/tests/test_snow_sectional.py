"""Independent equation benchmarks and conservation controls for the repair."""
import math

import numpy as np
import pytest

from atmosphere.snow import flakes as sf
from atmosphere.snow.habits import GEOMETRIES
from atmosphere.snow.sectional import CollisionOperator, column, integrate_box, seed_distribution

GEOMETRY = GEOMETRIES['bf95_dmax_dendrites_area']


def test_constant_kernel_matches_exact_discrete_distribution():
    # Monodisperse unit mass and number, K=1: n_k=(t/2)^(k-1)/(1+t/2)^(k+1).
    mass = np.arange(1., 81.)
    initial = np.zeros_like(mass); initial[0] = 1.
    operator = CollisionOperator(mass, np.ones((mass.size, mass.size)))
    time = 2.
    result = integrate_box(operator, initial, time, max_step_s=.004)
    exact = (time / 2)**(mass - 1) / (1 + time / 2)**(mass + 1)
    np.testing.assert_allclose(result['numbers'], exact, rtol=0, atol=7e-7)
    assert result['numbers'].sum() == pytest.approx(1 / (1 + time / 2), abs=3e-7)
    assert np.dot(mass, result['numbers']) == pytest.approx(1., abs=1e-13)
    assert result['upper_exit_mass'] < 1e-15


def test_additive_kernel_matches_exact_borel_distribution_and_moments():
    # For K(i,j)=i+j, n_k=e^-t exp(-kx)(kx)^(k-1)/k!, x=1-e^-t.
    mass = np.arange(1., 129.)
    initial = np.zeros_like(mass); initial[0] = 1.
    operator = CollisionOperator(mass, mass[:, None] + mass[None, :])
    time = .5
    result = integrate_box(operator, initial, time, max_step_s=.0005)
    x = 1 - math.exp(-time)
    exact = np.exp(-time - mass*x + (mass-1)*np.log(mass*x)
                   - np.array([math.lgamma(k + 1) for k in mass]))
    np.testing.assert_allclose(result['numbers'], exact, rtol=0, atol=3e-7)
    assert result['numbers'].sum() == pytest.approx(math.exp(-time), abs=5e-8)
    assert np.dot(mass, result['numbers']) == pytest.approx(1., abs=2e-13)
    assert np.dot(mass**2, result['numbers']) == pytest.approx(math.exp(2*time), rel=2e-6)


def test_rk2_time_refinement_has_second_order_error_reduction():
    mass = np.arange(1., 65.)
    initial = np.zeros_like(mass); initial[0] = 1.
    operator = CollisionOperator(mass, np.ones((mass.size, mass.size)))
    errors = []
    for step in (.04, .02, .01):
        result = integrate_box(operator, initial, 1., max_step_s=step)
        errors.append(abs(result['numbers'][0] - 1 / 1.5**2))
    assert errors[1] < .27 * errors[0]
    assert errors[2] < .27 * errors[1]


def test_geometric_mass_refinement_reduces_known_second_moment_error():
    # Constant-kernel exact M2(t)=1+t; fixed pivots conserve M0/M1 but diffuse M2.
    errors = []
    for bins in (25, 49, 97):
        mass = np.geomspace(1., 256., bins)
        initial = np.zeros_like(mass); initial[0] = 1.
        operator = CollisionOperator(mass, np.ones((bins, bins)))
        result = integrate_box(operator, initial, 1., max_step_s=.005)
        errors.append(abs(np.dot(result['numbers'], mass**2) - 2.))
        assert np.dot(result['numbers'], mass) + result['upper_exit_mass'] == pytest.approx(1., abs=1e-13)
    assert errors[1] < .6 * errors[0]
    assert errors[2] < .6 * errors[1]


def test_binary_breakup_channel_preserves_mass_and_retains_fractional_fragments():
    mass = np.arange(1., 33.)
    kernel = np.zeros((32, 32)); kernel[7, 15] = kernel[15, 7] = 2.
    losses = np.zeros_like(kernel); losses[7, 15] = losses[15, 7] = 3.2
    initial = np.zeros(32); initial[7], initial[15] = .002, .003
    operator = CollisionOperator(mass, kernel, sticking=.25, loss_per_breakup=losses, fragment_mass=1.)
    derivative, exit_number, exit_mass = operator.evaluate(initial)
    events = 2 * initial[7] * initial[15]
    assert derivative.sum() == pytest.approx(events * (-.25 + .75 * 3.2), abs=1e-19)
    assert derivative[0] == pytest.approx(events * .75 * 3.2)
    assert np.dot(derivative, mass) + exit_mass.sum() == pytest.approx(0., abs=1e-18)
    assert exit_number.sum() == 0.


def test_domain_exits_are_conserved_without_returning_mass_to_grid():
    mass = np.array([1., 2., 4.])
    kernel = np.ones((3, 3))
    operator = CollisionOperator(mass, kernel)
    initial = np.array([0., 0., .125])
    derivative, exit_number, exit_mass = operator.evaluate(initial)
    assert derivative[-1] < 0
    assert exit_number[1] == pytest.approx(.5 * .125**2)
    assert exit_mass[1] == pytest.approx(8 * exit_number[1])
    assert np.dot(derivative, mass) + exit_mass.sum() == pytest.approx(0., abs=1e-15)
    losses = np.full((3, 3), .5)
    operator = CollisionOperator(mass, kernel, sticking=0., loss_per_breakup=losses, fragment_mass=.25)
    derivative, exit_number, exit_mass = operator.evaluate(np.array([1., 0., 0.]))
    assert exit_number[0] == pytest.approx(1.5)
    assert exit_mass[0] == pytest.approx(.5)
    assert np.dot(derivative, mass) + exit_mass.sum() == pytest.approx(0., abs=1e-15)


def test_stiff_box_remains_positive_with_mass_accounting():
    mass = np.arange(1., 33.)
    initial = np.zeros(32); initial[0] = 1.; initial[-1] = 1e-200
    operator = CollisionOperator(mass, np.full((32, 32), 1e5))
    result = integrate_box(operator, initial, .01, max_step_s=1., cfl=.3)
    assert (result['numbers'] >= 0).all()
    assert np.dot(result['numbers'], mass) + result['upper_exit_mass'] == pytest.approx(1., abs=1e-12)
    assert result['numbers'][-1] > 0


def test_seed_projection_preserves_number_and_physical_mass_across_grids():
    low, high = 1e-4, 4e-4
    expected_mass = (GEOMETRY.mass_coefficient * (high**GEOMETRY.mass_exponent - low**GEOMETRY.mass_exponent)
                     / (GEOMETRY.mass_exponent * math.log(high / low)))
    for bins in (20, 47, 128):
        mass = np.geomspace(GEOMETRY.mass(2e-5), GEOMETRY.mass(1.), bins)
        seed = seed_distribution(mass, GEOMETRY)
        assert seed.sum() == pytest.approx(1., rel=2e-14)
        assert np.dot(seed, mass) == pytest.approx(expected_mass, rel=2e-14)


def test_column_area_independence_exact_gathering_and_mass_ledger():
    args = dict(geometry=GEOMETRY, layers=3, mass_bins=32, maximum_diameter_m=.03,
                fill=.1, gather=.2, max_step_s=13.)
    one = column(sf.earth_column(), -2., 1e4, 1., **args)
    many = column(sf.earth_column(), -2., 1e4, 1., horizontal_area_m2=256., **args)
    for key in ('arriving_numbers_m2', 'final_number_m3', 'time_mean_number_m3'):
        np.testing.assert_array_equal(one[key], many[key])
    assert one['mass_ledger_kg_m2'] == many['mass_ledger_kg_m2']
    for key, value in one['extensive_mass_ledger_kg'].items():
        assert many['extensive_mass_ledger_kg'][key] == 256 * value
    assert abs(one['mass_ledger_kg_m2']['relative_residual']) < 2e-13
    assert one['diagnostics']['final_time_s'] == pytest.approx(.3 * one['crossing_hours'] * 3600., rel=1e-14)
    assert one['seconds'] == pytest.approx(.2 * one['crossing_hours'] * 3600., rel=1e-14)
    assert np.dot(one['arriving_numbers_m2'], one['masses']) > 0
    assert one['domain_adequacy']['passed']


def test_column_reports_inadequate_mass_domain_and_keeps_exits_out_of_snowfall():
    result = column(sf.earth_column(), -2., 1e4, 3., geometry=GEOMETRY, layers=2, mass_bins=16,
                    minimum_diameter_m=5e-5, maximum_diameter_m=4e-4, fill=.5, gather=.5)
    ledger = result['mass_ledger_kg_m2']
    assert ledger['upper_domain_exit'] > 0
    assert not result['domain_adequacy']['passed']
    assert abs(ledger['relative_residual']) < 1e-12
    assert np.dot(result['arriving_numbers_m2'], result['masses']) <= ledger['freezing_exit']
    assert not result['domain_adequacy']['physical_size_cutoff']


def test_rejects_invalid_kernel_population_and_missing_seed_range():
    with pytest.raises(ValueError, match='symmetric'):
        CollisionOperator([1., 2.], [[1., 2.], [3., 1.]])
    with pytest.raises(ValueError, match='injected'):
        seed_distribution(GEOMETRY.mass(np.array([2e-4, 1e-3])), GEOMETRY)
    operator = CollisionOperator([1., 2.], np.ones((2, 2)))
    with pytest.raises(ValueError, match='nonnegative'):
        operator.evaluate([-1., 1.])


def test_top_ghost_allocation_limits_fractional_growth_without_losing_water():
    mass = np.array([1., 2., 4.])
    kernel = np.zeros((3, 3)); kernel[0, 2] = kernel[2, 0] = 1.
    operator = CollisionOperator(mass, kernel)
    derivative, exit_number, exit_mass = operator.evaluate(np.array([1., 0., 1e-100]))
    # Product mass 5 is 3/4 of the resolved mass-4 pivot and 1/4 of ghost mass 8.
    assert derivative[0] == pytest.approx(-1e-100, rel=1e-14, abs=0.)
    assert derivative[2] == pytest.approx(-.25e-100, rel=1e-14, abs=0.)
    assert exit_number[1] == pytest.approx(.25e-100, rel=1e-14, abs=0.)
    assert exit_mass[1] == pytest.approx(2e-100, rel=1e-14, abs=0.)
    assert np.dot(derivative, mass) + exit_mass.sum() == 0.


def test_vertical_refinement_preserves_boundary_forcing_and_exact_duration():
    args = dict(geometry=GEOMETRY, mass_bins=24, maximum_diameter_m=.03, fill=.01, gather=.02)
    coarse = column(sf.earth_column(), -2., 1e4, 1., layers=2, **args)
    fine = column(sf.earth_column(), -2., 1e4, 1., layers=7, **args)
    assert coarse['crossing_hours'] == fine['crossing_hours']
    assert coarse['seconds'] == fine['seconds']
    for key in ('injection_reference_velocity_m_s', 'injected_number_flux_m2_s', 'injected_mass_flux_kg_m2_s'):
        assert coarse['settings'][key] == fine['settings'][key]
    assert coarse['mass_ledger_kg_m2']['injected'] == pytest.approx(fine['mass_ledger_kg_m2']['injected'], rel=1e-14)
