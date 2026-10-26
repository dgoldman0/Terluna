"""Physical limits and independent controls for the snow optical screen."""
import math
import json

import numpy as np
import pytest

from atmosphere.snow import optics
from shared.provenance import constants_changed
from shared.baselines import recorded_file


def test_extinction_keeps_area_under_population_rebinning():
    # Two distinct populations combined into one bin must retain their total
    # number-weighted area, even though their representative diameter changes.
    resolved = optics.summarize_distribution([.001, .2], [100., 1.], [1e-7, .01])
    merged = optics.summarize_distribution([.003], [101.], [(100e-7+.01)/101.])
    assert resolved['extinction_m_inverse'] == pytest.approx(merged['extinction_m_inverse'])
    beta = resolved['extinction_m_inverse']
    assert math.exp(-beta*resolved['beam_5_percent_distance_m']) == pytest.approx(.05)
    assert resolved['q_ext_one_sensitivity']['beam_5_percent_distance_m'] == pytest.approx(2*resolved['beam_5_percent_distance_m'])
    assert sum(x['extinction_share'] for x in resolved['size_contributions']) == pytest.approx(1.)
    assert next(x['extinction_share'] for x in resolved['size_contributions'] if x['diameter_min_m'] == .1) > .99


def test_submetre_component_does_not_redistribute_excluded_mass_or_area():
    result = optics.summarize_distribution([.2, 2., 20.], [10., .1, .001], [.01, 1., 100.])
    assert result['extinction_m_inverse'] == pytest.approx(.6)
    assert result['below_one_metre_component']['extinction_m_inverse'] == pytest.approx(.2)
    assert result['below_one_metre_component']['fraction_of_total_extinction'] == pytest.approx(1/3)
    assert result['size_contributions'][-1]['extinction_share'] == pytest.approx(1/3)


def test_empty_snow_and_invalid_distributions():
    empty = optics.summarize_distribution([], [], [])
    assert empty['extinction_m_inverse'] == 0.
    assert empty['beam_5_percent_distance_m'] is None
    assert all(x['fraction'] == 1 for x in empty['transmission'])
    with pytest.raises(ValueError):
        optics.summarize_distribution([.1], [-1.], [.01])
    with pytest.raises(ValueError):
        optics.summarize_distribution([.1], [1., 2.], [.01])


@pytest.mark.parametrize('g', [0., .7, .8, .9])
def test_phase_probability_and_mean_cosine(g):
    mu, weights = np.polynomial.legendre.leggauss(200)
    probability = 2*np.pi*weights*optics.phase_hg(mu, g)
    assert probability.sum() == pytest.approx(1., abs=1e-11)
    assert (mu*probability).sum() == pytest.approx(g, abs=1e-11)


def test_isotropic_light_preserves_colour_and_zero_path():
    source_xyz = np.array([2., 3., 4.])
    scalar = source_xyz*4*np.pi
    for g in (.0, .9):
        result = optics.thin_radiance([0., 0., 0.], scalar, 20., tau=.05, g=g)
        assert result == pytest.approx(source_xyz*-math.expm1(-.05))
    assert optics.thin_radiance(source_xyz, scalar, 90., tau=0.) == pytest.approx(np.zeros(3))
    with pytest.raises(ValueError):
        optics.thin_radiance(source_xyz, scalar, 90., tau=1.)


def test_vertical_tangent_and_occulted_shell_geometry():
    edges = np.array([10., 12., 15.])
    assert optics.shell_path_lengths(11., 1., edges) == pytest.approx([1., 3.])
    tangent = optics.shell_path_lengths(11., 0., edges)
    assert tangent == pytest.approx([math.sqrt(12**2-11**2), math.sqrt(15**2-11**2)-math.sqrt(12**2-11**2)])
    assert optics.shell_path_lengths(11., -1., edges) is None


def test_finite_sun_in_vacuum_and_below_local_horizon():
    edges = np.array([1e6, 1.1e6])
    zero = np.zeros((1, 2))
    assert optics.direct_transmission(0., 90., edges, zero) == pytest.approx([1., 1.])
    assert optics.direct_transmission(0., 0., edges, zero) == pytest.approx([.5, .5])
    assert optics.direct_transmission(0., -10., edges, zero) == pytest.approx([0., 0.])
    assert optics.direct_transmission(10000., -5., edges, zero) == pytest.approx([1., 1.])
    elevated_horizon = -math.degrees(math.acos(edges[0]/(edges[0]+10000)))
    assert optics.direct_transmission(10000., elevated_horizon, edges, zero) == pytest.approx([.5, .5], abs=1e-12)
    assert optics.direct_transmission(10000., elevated_horizon-1., edges, zero) == pytest.approx([0., 0.])


def test_conservative_multiple_scattering_energy_and_beer_lambert_control():
    result = optics.slab_flux(2., g=.8, photons=30000, seed=201)
    assert result['energy_residual'] == 0.
    assert result['transmitted']+result['reflected'] == 1.
    assert abs(result['unscattered']-math.exp(-2)) < 5*result['unscattered_standard_error']
    assert result['transmitted'] > result['unscattered']
    assert result['mean_scatter_events'] > 1.


def test_absorbing_slab_and_vacuum_controls():
    empty = optics.slab_flux(0., photons=1000)
    assert empty['transmitted'] == empty['unscattered'] == 1.
    assert empty['mean_scatter_events'] == 0.
    opaque = optics.slab_flux(1., mu0=.5, albedo=0., photons=30000, seed=892)
    assert opaque['reflected'] == 0.
    assert opaque['mean_scatter_events'] == 0.
    assert opaque['energy_residual'] == 0.
    assert opaque['transmitted'] == opaque['unscattered']
    assert abs(opaque['transmitted']-math.exp(-2)) < 5*opaque['transmitted_standard_error']
    with pytest.raises(RuntimeError):
        optics.slab_flux(10., photons=1000, max_events=1)


def test_cloud_droplets_area_mass_scaling():
    beta = optics.droplet_extinction(1e-4, 10e-6)
    assert optics.droplet_extinction(2e-4, 20e-6) == pytest.approx(beta)
    assert .0149 < beta < .0151


def test_height_interpolation_keeps_scalar_irradiance():
    r = np.array([10., 20.])
    a = np.array([-1., 1.])
    f = np.array([[1., 3.], [5., 7.]])
    assert optics._bilinear(r, a, f, 15., 0.) == 4.
    with pytest.raises(ValueError):
        optics._bilinear(r, a, f, 21., 0.)


def test_optics_rejects_populations_without_the_supplied_audit_assessment(tmp_path, monkeypatch):
    monkeypatch.setattr(optics, 'ROOT', tmp_path)
    population, depth, audit = [tmp_path / (name+'.json') for name in ('population', 'depth', 'audit')]
    population.write_text(json.dumps({'schema': 'terluna.atmosphere.snow-population/1'}))
    depth.write_text(json.dumps({'schema': 'terluna.atmosphere.snow-growth-depth/1'}))
    audit.write_text(json.dumps({'schema': 'terluna.atmosphere.snow-numerical-audit/1'}))
    with pytest.raises(ValueError, match='requires a current assessment'):
        optics.build(population, growth_depth_path=depth, numerical_audit_path=audit)
    # A warning tied to a different audit is also insufficient evidence.
    population.write_text(json.dumps({
        'schema': 'terluna.atmosphere.snow-population/1',
        'numerical_assessment': {'schema': 'terluna.atmosphere.snow-assessment/1',
                                 'audit': {'path': 'audit.json', 'sha256': 'stale'}}}))
    with pytest.raises(ValueError, match='requires a current assessment'):
        optics.build(population, growth_depth_path=depth, numerical_audit_path=audit)


def test_committed_optics_product_records_conserved_inputs_and_light_height():
    if not optics.OUT.exists():
        pytest.skip('optical product has not yet been generated')
    product = json.loads(optics.OUT.read_text())
    assert product['schema'] == optics.SCHEMA
    for path, expected in product['producer']['files'].items():
        assert optics.digest(recorded_file(path, expected)) == expected
    assert not constants_changed(product['producer']['constants'])
    for path, expected in product['inputs'].items():
        assert optics.digest(recorded_file(path, expected)) == expected
    population = json.loads(optics.POPULATION.read_text())
    growth_depth = json.loads(optics.GROWTH_DEPTH.read_text())
    audit = json.loads(optics.NUMERICAL_AUDIT.read_text())
    for source in (population, growth_depth):
        assert source['numerical_assessment']['audit']['sha256'] == optics.digest(optics.NUMERICAL_AUDIT)
        for path in ('atmosphere/snow/aggregation.py', 'atmosphere/snow/habits.py'):
            assert source['producer']['files'][path] == audit['producer']['files'][path] == optics.digest(optics.ROOT/path)
    for run in population['runs'] + growth_depth['runs']:
        row = next(r for r in product['distributions'] if r['run_id'] == run['id'] and r['level'] == 'freezing_level')
        d = run['freezing_level']
        expected = 2*np.dot(d['number_m3'], d['projected_area_m2'])
        assert row['extinction_m_inverse'] == pytest.approx(expected)
    assert {r['altitude_m'] for r in product['lighting']['rows']} == {29200., 35000., 40000., 45500.}
    assert max(abs(r['relative_difference']) for r in product['lighting']['surface_direct_controls']) < .005
    assert max(abs(r['diffuse_scalar_relative_difference']) for r in product['lighting']['coarse_moment_controls']) < .03
    assert max(abs(r['direct_normal_relative_difference']) for r in product['lighting']['finite_disk_controls']) < .001
    for slab in product['multiple_scattering_slabs']:
        assert slab['energy_residual'] == 0.
        assert abs(slab['unscattered']-slab['unscattered_exact']) < 6*slab['unscattered_standard_error']+4/slab['photons']
