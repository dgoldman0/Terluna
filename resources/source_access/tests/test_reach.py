"""The velocity change each near-Earth source asks of freight, and the survey's water within each reach."""
import gzip
import json
import math
import numpy as np
import pytest
from resources.source_access import reach as r


def test_formulas_reproduce_shoemaker_and_helin_table_1():
    for name, a, e, i, F, dv, at in r.TABLE_1:
        assert abs(float(r.figure_of_merit(a, e, i, at)) - F) < 0.007, name
        assert abs(float(r.leo_rendezvous_km_s(a, e, i, at)) - dv) < 0.1, name


def test_an_orbit_like_the_earths_needs_only_the_correction_from_the_depot():
    assert float(r.depot_leg_m_s(1.0, 0.0, 0.0)) == pytest.approx(500.0, abs=1e-6)
    # From low Earth orbit the same body asks for escape from low orbit as well.
    escape = math.sqrt(2 * r.EARTH_GM / (r.EARTH_RADIUS + 1000 * r.LEO_ALTITUDE_KM)) / 1000
    circular = escape / math.sqrt(2)
    assert float(r.leo_rendezvous_km_s(1.0, 0.0, 0.0)) == pytest.approx(escape - circular + 0.5, abs=1e-9)


def test_the_depot_leg_saves_the_escape_from_low_orbit_and_grows_with_inclination_and_aphelion():
    a = np.array([0.9, 1.1, 1.4, 2.0, 2.6])
    e = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    for i in (0.0, 10.0, 25.0):
        leg = r.depot_leg_m_s(a, e, i) / 1000
        leo = r.leo_rendezvous_km_s(a, e, i)
        assert np.all(leg < leo)
    assert np.all(np.diff([float(r.depot_leg_m_s(1.4, 0.2, i)) for i in (0, 5, 10, 20, 30)]) > 0)
    # With the perihelion held near the Earth's orbit, a farther aphelion asks more; at a fixed semimajor axis a larger
    # eccentricity can ask less, as the perihelion comes nearer the Earth's orbit.
    q = 1.05
    aphelia = (1.2, 1.6, 2.0, 2.6, 3.2)
    assert np.all(np.diff([float(r.depot_leg_m_s((q + Q) / 2, (Q - q) / (Q + q), 5.0)) for Q in aphelia]) > 0)


def test_bins_take_the_first_upper_edge_not_below_and_hold_the_rest_beyond():
    edges = np.asarray(r.EDGES_M_S)
    assert list(r.bin_of([0.0, 100.0, 100.5, edges[-1], edges[-1] + 1, np.inf])) == [0, 0, 1, len(edges) - 1, len(edges), len(edges)]


@pytest.mark.skipif(not (r.HERE / 'results/reach.json.gz').exists(), reason='The product has not been built')
def test_product_reproduces_the_survey_and_is_current():
    product = json.loads(gzip.decompress((r.HERE / 'results/reach.json.gz').read_bytes()))
    for path, value in product['dependencies'].items():
        assert r.digest(r.ROOT / path) == value, path
    assert product['producer']['sha256'] == r.digest(r.HERE / 'reach.py')
    assert all(product['catalogued_pools_reproduced'].values())
    survey = json.loads(gzip.decompress(r.POPULATION.read_bytes()))
    stated = survey['scenarios']['central']['pools']['all_nea']['total_quantiles_kg']['water_equivalent_kg']
    within = product['contained_water_within_kg']
    for key in ('p05', 'p50', 'p95'):
        values = np.asarray(within[key])
        assert np.all(np.diff(values) >= -1e-6 * values[-1])
        assert values[-1] == pytest.approx(stated[key], rel=1e-9)
    assert len(within['p50']) == len(product['depot_leg_upper_edges_m_s']) + 1
    for row in product['table_1_check']:
        assert abs(row['leo_rendezvous_km_s'] - row['published_km_s']) < 0.1
