"""Platelet barriers, dissolved loss, shadows, gas cells, food and falls."""
import json
import math

import pytest

from research.studies.floater_viability import gas_biology as gb
from research.studies.floater_viability import resources as r


def test_platelet_barrier_limits():
    assert gb.platelet_barrier(0., 1000.) == 1.
    assert gb.platelet_barrier(.02, 5000.) < gb.platelet_barrier(.02, 500.) < gb.platelet_barrier(.02, 100.)
    # A hundredfold cut, as stripping the guanine layer reverses, needs alpha*phi near 200.
    assert gb.platelet_barrier(.04, 5000.) == pytest.approx(.96/101.)


def test_dissolved_loss_is_linear_and_the_exchanger_returns_its_share():
    a = gb.dissolved_gas_loss(1., 1e5)
    b = gb.dissolved_gas_loss(2., 1e5)
    assert b['mol_m2_s'] == pytest.approx(2*a['mol_m2_s'])
    c = gb.dissolved_gas_loss(1., 1e5, exchanger_efficiency=.9)
    assert c['mol_m2_s'] == pytest.approx(.1*a['mol_m2_s'])
    out = gb.evaluate()
    assert json.loads(json.dumps(out, allow_nan=False)) == out


def test_shadow_geometry():
    theta = r.SUN_ANGULAR_DIAMETER_RAD
    s = r.shadow(200., 10000.)
    assert s['umbra_diameter_m'] == pytest.approx(200-10000*theta)
    assert s['umbra_length_m'] == pytest.approx(200/theta)
    small = r.shadow(40., 10000.)
    assert not small['full_shadow'] and small['peak_light_removed'] < .49
    assert math.degrees(theta) == pytest.approx(.533, abs=2e-3)


def test_gas_cells_food_and_falls():
    one = r.gas_cell_area(100., 20., 1)
    assert one == pytest.approx(2*math.pi*100+2*math.pi*10*100)
    food = r.food_people(1000., 1., .5, .1)
    assert food['people'] == pytest.approx(50*4000/(2500*r.JULIAN_YEAR_DAYS))
    fall = r.draped_fall(10., 1.4)
    g = r.MOON_GM/r.MOON_RADIUS**2
    assert .5*1.4*1.3*fall['sink_m_s']**2 == pytest.approx(10*g)
    out = r.evaluate()
    assert json.loads(json.dumps(out, allow_nan=False)) == out
    assert out['sky_ship_cells'][0]['length_m'] == 245.


def test_same_loss_thickness_is_the_permeability_ratio():
    out = gb.evaluate()
    for row in out['platelets']:
        assert row['thickness_for_same_loss_ratio'] == pytest.approx(row['permeability_ratio'])


def test_fall_through_a_uniform_profile_matches_the_one_density_fall():
    one = r.draped_fall(10., 1.2, height_m=10000.)
    integrated = r.fall_minutes(10., 10000., [0., 20.], [1.2, 1.2])
    # draped_fall takes gravity at the start height; the integral takes it at each height (1% over 10 km)
    assert integrated == pytest.approx(one['minutes_from_height'], rel=5e-3)
    assert integrated < one['minutes_from_height']
    assert r.fall_minutes(10., 10000., [0., 20.], [1.4, .9]) < r.fall_minutes(10., 10000., [0., 20.], [1.4, 1.4])
