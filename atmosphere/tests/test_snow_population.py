"""Conservation and compatibility checks for the snow distribution export."""
import numpy as np
import pytest

from atmosphere.snow import aggregation as ag, flakes as sf
from atmosphere.snow.habits import GEOMETRIES
from atmosphere.snow.population import histogram, distribution, arriving_distribution, CloudSamples, ArrivalHistory


def test_binning_preserves_number_mass_area_and_flux_without_median_shortcut():
    d = np.array([1e-4, 0.001, 0.01001, 0.01002, 0.7])
    n = np.array([100., 20., 3., 7., 0.01])
    m, a, v = 0.0121*d**1.9, 0.25*np.pi/4*d*d, np.array([0.1, 0.2, 0.3, 0.5, 0.6])
    out = distribution(histogram(d, n, m, a, v))
    got = out['totals']
    for key, value in [('number_m3', n.sum()), ('mass_kg_m3', (n*m).sum()),
                       ('projected_area_m_inv', (n*a).sum()), ('number_flux_m2_s', (n*v).sum()),
                       ('mass_flux_kg_m2_s', (n*m*v).sum())]:
        assert got[key] == pytest.approx(value, rel=1e-13)
    assert np.dot(out['number_m3'], out['projected_area_m2']) == pytest.approx((n*a).sum())


def test_bins_refuse_to_silently_drop_outside_sizes():
    with pytest.raises(ValueError, match='exceeds'):
        histogram(np.array([1e-6]), *[np.ones(1) for _ in range(4)])


def test_legacy_geometry_preserves_the_original_stochastic_trajectory():
    earth = sf.earth_column()
    args = dict(top_c=-2., number_per_m3=1e4, snowfall_mm_h=1., bins=4, inject=4, fill=1., gather=0.5)
    old = ag.column(earth, **args)
    new = ag.column(earth, geometry=GEOMETRIES['legacy_bf95'], **args)
    np.testing.assert_array_equal(new['masses'], old['masses'])
    np.testing.assert_array_equal(new['numbers'], old['numbers'])


def test_population_ledger_and_flux_concentration_conversion():
    earth, geometry = sf.earth_column(), GEOMETRIES['bf95_dmax']
    observer = CloudSamples(earth, geometry, -2., layers=3, sample_s=60.)
    result = ag.column(earth, -2., 1e4, 1., bins=4, inject=4, fill=1., gather=0.5,
                       geometry=geometry, observer=observer)
    ledger = result['mass_ledger_kg_m2']
    assert abs(ledger['residual']) < 1e-10 * (ledger['injected'] + ledger['deposited'])
    out = arriving_distribution(result, earth, geometry)
    assert out['totals']['mass_flux_kg_m2_s'] == pytest.approx(
        np.dot(result['numbers'], result['masses'])/result['seconds'], rel=1e-13)
    assert observer.count > 0
    assert sum(x['totals']['mass_kg_m3'] for x in observer.export()) > 0


def test_arrival_history_conserves_flux_and_keeps_empty_windows():
    earth, geometry = sf.earth_column(), GEOMETRIES['bf95_dmax']
    history = ArrivalHistory(earth, geometry, window_s=100.)
    mass = np.array([1e-6, 2e-6])
    number = np.array([20., 30.])
    history(mass, number, 10., 1.)
    history(mass, number, 90., 1.)
    history(mass, number, 250., 1.)
    rows = history.export()
    assert len(rows) == 2  # incomplete third window deliberately excluded
    assert rows[0]['totals']['mass_flux_kg_m2_s'] == pytest.approx(2*np.dot(mass, number)/100.)
    assert rows[0]['totals']['number_flux_m2_s'] == pytest.approx(2*number.sum()/100.)
    assert rows[1]['totals']['mass_flux_kg_m2_s'] == 0.


def test_committed_population_preserves_input_identity_and_extensive_moments():
    import json
    from atmosphere.snow import population as pop
    from shared.provenance import constants_changed
    if not pop.OUT.exists():
        pytest.skip('Population product not generated yet')
    product = json.loads(pop.OUT.read_text())
    assert product['schema'] == pop.SCHEMA
    assert not constants_changed(product['producer']['constants'])
    for path, expected in product['producer']['files'].items():
        assert pop.sha(pop.ROOT/path) == expected, path
    for path, expected in product['producer']['inputs'].items():
        assert pop.sha(pop.ROOT/path) == expected, path
    assert {r['id'] for r in product['runs']} == {r['id'] for r in pop.cases()}
    for run in product['runs']:
        assert abs(run['mass_ledger_kg_m2']['relative_residual']) < 1e-10
        boundary = run['freezing_level']
        assert boundary['totals']['mass_flux_kg_m2_s'] * 3600 == pytest.approx(run['summary']['snowfall_mm_h'])
        for d in [boundary] + run['layers']:
            n, a, mass, flux, mass_flux = [np.asarray(d[k]) for k in
                ('number_m3', 'projected_area_m2', 'mass_kg_m3', 'number_flux_m2_s', 'mass_flux_kg_m2_s')]
            assert all(np.isfinite(x).all() and (x >= 0).all() for x in (n, a, mass, flux, mass_flux))
            for key, value in [('number_m3', n.sum()), ('mass_kg_m3', mass.sum()),
                               ('projected_area_m_inv', (n*a).sum()), ('number_flux_m2_s', flux.sum()),
                               ('mass_flux_kg_m2_s', mass_flux.sum())]:
                assert d['totals'][key] == pytest.approx(value, rel=1e-12)


def test_rare_large_model_tail_keeps_its_area_and_mass():
    from atmosphere.snow.population import large_tail
    geometry = GEOMETRIES['bf95_dmax_dendrites_area']
    d = np.array([.1, 1120.])
    m, n = geometry.mass(d), np.array([100., .01])
    out = distribution(histogram(d, n, m, geometry.area(d), np.full(2, .6)))
    assert out['totals']['projected_area_m_inv'] == pytest.approx(np.dot(n, geometry.area(d)))
    assert out['totals']['mass_kg_m3'] == pytest.approx(np.dot(n, m))
    tail = large_tail(dict(masses=m, numbers=n, seconds=3600.), geometry)
    assert tail['maximum_sampled_diameter_m'] == pytest.approx(1120.)
    assert tail['thresholds'][-1]['arrivals_per_m2_hour'] == pytest.approx(.01)


def test_open_tail_has_no_invented_upper_edge_or_lost_moments():
    geometry = GEOMETRIES['bf95_dmax_dendrites_area']
    d = np.array([.1, 11., 1e4])
    n = np.array([100., .01, 1e-12])
    m, a = geometry.mass(d), geometry.area(d)
    out = distribution(histogram(d, n, m, a, np.ones(3)))
    assert out['open_ended_last_bin']
    assert out['size_edges_m'][-1] is None
    assert out['size_edges_m'][-2] == 10.
    assert out['totals']['mass_kg_m3'] == pytest.approx(np.dot(n, m))
    assert np.dot(out['number_m3'], out['projected_area_m2']) == pytest.approx(np.dot(n, a))
    assert out['number_m3'][-1] == pytest.approx(n[1:].sum())
