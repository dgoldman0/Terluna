import json
import math

import pytest

from shared import constants as K
from shared.provenance import constants_changed
from . import array_heat as ah


def test_a_film_at_the_sunward_crossing_radiates_what_it_absorbs_from_both_faces():
    row = ah.tile(.03, .1, 4.)
    emitted = (row['emissivity_sunward']+row['emissivity_shaded'])*K.STEFAN_BOLTZMANN*row['temperature_K']**4
    assert emitted == pytest.approx(row['absorbed_W_m2'], rel=1e-6)
    # A few micrometres of silica carry the shaded face's share with well under a millikelvin across them.
    assert row['across_film_K'] < 1e-3 and row['carnot_across_film'] < 1e-5


def test_the_films_emit_in_the_thermal_infrared_as_silica_does():
    (front, back), thickness = ah.emissivity(.1, 4., 180.)
    assert .2 < front < .7 and .2 < back < .7
    (window_front, window_back), _ = ah.emissivity(1., 10., 180.)
    assert window_back > back                    # ten micrometres of silica emit more than four
    assert thickness > 4.1


def test_a_rotor_stores_momentum_cheapest_wide_and_slow():
    h, r, v = 2e9, 100., 300.
    mass, energy = h/(v*r), .5*h*v/r
    assert energy == pytest.approx(.5*mass*v**2)
    # Four times the radius at the same speed: a quarter of the mass and of the energy.
    assert .5*h*v/(4*r) == pytest.approx(energy/4)


def test_the_moon_intercepts_a_small_share_of_heat_released_at_the_ring_radius():
    share = ah.placement(dict(requirements=dict(climate_top_of_air_W_m2=1173.25)))['moon_share_of_isotropic_heat']
    ratio = K.MOON_RADIUS/(ah.DESIGN['orbit_radius_km']*1e3)
    assert share == pytest.approx(ratio**2/4, rel=1e-2)


@pytest.mark.skipif(not ah.OUT.exists(), reason='array_heat.json not written')
def test_product_is_current():
    product = json.loads(ah.OUT.read_text())
    assert product['schema'] == ah.SCHEMA
    for name, value in product['producer']['files'].items():
        assert ah.digest(ah.ROOT/name) == value, name
    for name, value in product['inputs'].items():
        assert ah.digest(ah.ROOT/name) == value, name
    assert not constants_changed(product['producer']['constants'])
    assert product['design'] == json.loads(json.dumps(ah.DESIGN))
    assert all(product['checks'].values())
    # The fleet the comparison carries: about 28 million tiles and 46 Gt.
    assert 2.7e7 < product['fleet']['tiles'][0] < product['fleet']['tiles'][1] < 2.9e7
    assert all(math.isclose(m, 46., rel_tol=.02) for m in product['fleet']['optical_mass_Gt'])
