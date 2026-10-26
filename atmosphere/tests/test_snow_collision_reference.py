"""Independent exact-event, capture and bookkeeping controls."""
import math

import numpy as np
import pytest

from atmosphere.snow import collision_reference as c, flakes as sf
from shared.baselines import recorded_file


def test_two_particles_have_exact_waiting_time_and_mass_sum():
    seed = 42
    expected = np.random.default_rng(seed).exponential(3./2.)
    result = c.coagulate([1., 4.], 3., [0., expected*.99, expected*1.01],
                         lambda a, b: 2., seed=seed)
    assert [len(s['masses_kg']) for s in result['snapshots']] == [2, 2, 1]
    np.testing.assert_array_equal(result['snapshots'][-1]['masses_kg'], [5.])
    np.testing.assert_array_equal(result['snapshots'][-1]['monomer_counts'], [2])


def test_three_particle_constant_kernel_matches_analytic_death_process():
    # N: 3 -> 2 at rate3, then 2 -> 1 at rate1 for K=V=1.
    t, repetitions = .4, 3000
    counts = np.zeros(3)
    for seed in range(repetitions):
        result = c.coagulate([1., 1., 1.], 1., [t], lambda a, b: 1., seed=seed)
        counts[len(result['snapshots'][0]['masses_kg'])-1] += 1
    p3 = math.exp(-3*t)
    p2 = 1.5*(math.exp(-t)-p3)
    probability = np.array([1-p2-p3, p2, p3])
    error = np.abs(counts/repetitions-probability)
    assert np.all(error < 6*np.sqrt(probability*(1-probability)/repetitions))


def test_pure_birth_constant_rate_matches_poisson_mean_variance_and_tail():
    result = c.pure_birth(1., .01, 2., [0., 1., 3.], lambda a, b: 1.,
                          samples=8192, maximum_events=100, seed=3)
    assert np.all(result['counts'][0] == 0)
    for index, expected in enumerate((2., 6.), start=1):
        sample = result['counts'][index]
        assert abs(sample.mean()-expected) < 6*math.sqrt(expected/sample.size)
        assert abs(sample.var()-expected) < 6*math.sqrt((expected+2*expected**2)/sample.size)
    # P[X=0] = exp(-2), checked separately from the first two moments.
    observed = np.mean(result['counts'][1] == 0)
    expected = math.exp(-2)
    assert abs(observed-expected) < 6*math.sqrt(expected*(1-expected)/8192)
    assert result['log_event_cap_probability_bound'] < math.log(1e-12)


def test_physical_kernel_keeps_identical_mass_encounters_zero_and_mass_conserved():
    world, _ = sf.moon_column()
    kernel = c.SnowKernel(world)
    masses = c.ordinary_initial()
    np.testing.assert_array_equal(kernel(masses, masses), np.zeros_like(masses))
    result = c.coagulate(np.tile(masses, 8), len(masses)*8/1e5,
                         [0., 7200., 21600.], kernel, seed=51)
    for snapshot in result['snapshots']:
        assert snapshot['masses_kg'].sum() == pytest.approx(np.tile(masses, 8).sum(), rel=1e-13)
        assert snapshot['monomer_counts'].sum() == len(masses)*8
        assert (snapshot['masses_kg'] > 0).all()


def test_nested_grids_retain_identical_initial_particles_and_capture_parents():
    base = c.nested_grid()
    initial = c.ordinary_initial()
    for factor in (1, 2, 4):
        masses = c.nested_grid(factor)
        np.testing.assert_array_equal(masses[::factor], base)
        numbers = c._project(masses, initial, len(initial)/1e5)
        np.testing.assert_array_equal(masses[numbers > 0], initial)
        assert numbers.sum() == pytest.approx(1e5)
        assert np.dot(numbers, masses) == pytest.approx(initial.mean()*1e5, rel=1e-15)


def test_prescribed_bath_has_separate_number_water_and_domain_ledgers():
    masses = np.array([.01, 1., 1.2, 1.44, 1.728, 2.1, 4., 8.])
    result = c.bath_fixed_pivot(masses, 1., .01, 2., [0., 1., 3.], lambda a, b: 1.)
    for row in result:
        assert abs(row['number_residual']) < 1e-12
        assert abs(row['relative_mass_residual']) < 1e-12
        assert row['exponential_composition_error'] < 1e-12
        assert (row['numbers'] >= 0).all()
    assert result[-1]['supplied_capture_mass_kg'] == pytest.approx(.06, rel=1e-8)


def test_insufficient_capture_event_domain_is_rejected():
    with pytest.raises(ValueError, match='insufficiently bounded'):
        c.pure_birth(1., .01, 2., [3.], lambda a, b: 1., samples=10, maximum_events=5)


def test_unequal_pair_rates_select_the_exact_first_collision_distribution():
    def kernel(a, b):
        small, large = np.minimum(a, b), np.maximum(a, b)
        return (1.*((small == 1.) & (large == 2.))
                + 2.*((small == 1.) & (large == 4.))
                + 3.*((small == 2.) & (large == 4.)))
    repetitions = 3000
    counts = np.zeros(3)
    for seed in range(repetitions):
        result = c.coagulate([1., 2., 4.], 1., [100.], kernel, seed=seed)
        final = tuple(sorted(result['snapshots'][0]['masses_kg']))
        counts[{(3., 4.): 0, (2., 5.): 1, (1., 6.): 2}[final]] += 1
    probability = np.array([1., 2., 3.])/6.
    assert np.all(np.abs(counts/repetitions-probability)
                  < 6*np.sqrt(probability*(1-probability)/repetitions))


def test_saved_product_identities_populations_and_ledgers():
    import hashlib
    import json
    product = json.loads(c.OUT.read_text())
    assert product['schema'] == c.SCHEMA
    for key in ('files', 'inputs'):
        for path, expected in product['producer'][key].items():
            assert hashlib.sha256(recorded_file(path, expected).read_bytes()).hexdigest() == expected
    for grid in product['capture']['fixed_pivot']:
        mass = np.array(grid['mass_pivots_kg'])
        for row in grid['snapshots']:
            number = np.array(row['numbers_per_initial_collector'])
            assert np.all(number >= 0)
            assert number.sum()+row['upper_exit_number'] == pytest.approx(1., abs=1e-12)
            expected = product['capture']['initial_collector_mass_kg']+row['supplied_capture_mass_kg']
            assert number@mass+row['upper_exit_mass_kg'] == pytest.approx(expected, rel=1e-12)
    for ensemble in product['ordinary_box']['exact_event_ensembles']:
        assert ensemble['maximum_relative_mass_residual'] < 1e-12
        assert ensemble['maximum_rate_sum_relative_discrepancy'] < 1e-11
        for row in ensemble['snapshots']:
            assert np.all(np.diff(row['mass_cdf_mean']) >= -1e-13)
            assert len(row['maximum_particle_mass_kg_per_realization']) == ensemble['independent_realizations']
    for grid in product['ordinary_box']['fixed_pivot']+[product['ordinary_box']['timestep_control']]:
        for row in grid['snapshots']:
            assert min(row['numbers_m3']) >= 0
            assert abs(row['relative_mass_residual']) < 1e-12
            assert np.all(np.diff(row['mass_cdf']) >= -1e-13)
