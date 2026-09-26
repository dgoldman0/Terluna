"""Checks on Earth's tidal effect on molecular escape."""
import json
import math

import numpy as np
import pytest

from shared.constants import EARTH_MOON_DISTANCE, MOON_RADIUS
from atmosphere.loss_response import tides as td


def test_l1_matches_the_earth_moon_values():
    c_l1, x_l1 = td.l1_jacobi()
    assert c_l1 == pytest.approx(3.1883, abs=2e-4)
    assert (1 - td.MU - x_l1) * EARTH_MOON_DISTANCE / 1e3 == pytest.approx(58000, rel=0.01)
    assert td.HILL * EARTH_MOON_DISTANCE / MOON_RADIUS == pytest.approx(35.3, abs=0.1)


def test_the_integration_keeps_the_jacobi_constant():
    rng = np.random.default_rng(3)
    state, _, _, _, r_c = td.launch(3.0, 200, rng)
    outcome, _, drift = td.follow(state, r_c, t_max=10 * td.DAY_UNITS)
    assert drift.max() < 1e-6
    assert set(np.unique(outcome)) <= {1, 2, 3}


def test_launch_speeds_lie_between_the_l1_threshold_and_two_body_escape():
    rng = np.random.default_rng(4)
    state, v2, v_min2, v_esc2, r_c = td.launch(3.0, 500, rng)
    assert np.all(v2 >= v_min2) and np.all(v2 <= v_esc2)
    assert np.allclose(td.moon_distance(*state[:3]), r_c)
    # outward only
    up = (state[:3] - np.array([1 - td.MU, 0, 0])[:, None]) / r_c
    assert np.all(np.sum(up * state[3:], axis=0) > 0)


def test_multiplier_limits():
    v_esc2 = 1.0
    v_min2 = np.full(4000, 0.87)
    v2 = np.linspace(0.87, 1.0, 4000)
    none_escape = np.ones(4000, dtype=int)
    all_escape = np.full(4000, 2)
    m_none, _, bound = td.multipliers(v2, v_min2, v_esc2, none_escape, lambdas=(11.0,))
    m_all, _, _ = td.multipliers(v2, v_min2, v_esc2, all_escape, lambdas=(11.0,))
    assert m_none[0] == 1.0
    assert m_all[0] == pytest.approx(bound[0], rel=1e-3)


def test_erkaev_factor_is_the_feasibility_baselines():
    xi = td.HILL * EARTH_MOON_DISTANCE / (3.0 * MOON_RADIUS)
    k = 1 - 1.5 / xi + 0.5 / xi ** 3
    assert td.erkaev(3.0, 11.0) == pytest.approx((1 + k * 11) / 12 * math.exp(11 * (1 - k)))


def test_stored_table_lies_between_jeans_and_the_energy_bound():
    table = json.loads((td.HERE / 'results' / 'tidal_escape.json').read_text())
    assert table['schema'] == td.SCHEMA
    for row in table['radii']:
        m, bound = np.array(row['multiplier']), np.array(row['energy_bound'])
        assert np.all(m >= 1.0) and np.all(m <= bound * 1.02)
        assert np.all(np.diff(m) > 0)
        assert row['escape_fraction_above_two_body_escape'] > 0.99
    at3 = td.multiplier(table, 3.0, 11.0)
    assert 1.8 < at3 < 2.8
    assert td.multiplier(table, 4.0, 11.0) > at3 > td.multiplier(table, 2.0, 11.0)
