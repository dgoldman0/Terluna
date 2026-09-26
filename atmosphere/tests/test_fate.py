"""Checks on where the escaping air goes."""
import json
import math

import numpy as np
import pytest

from atmosphere.loss_response import fate
from atmosphere.loss_response import tides as td


def test_launch_follows_the_outgoing_flux_above_the_l1_energy():
    rng = np.random.default_rng(5)
    state, eps, eps_min, r_c = fate.launch_escaping(3.0, 11.0, 40000, rng)
    assert np.all(eps >= eps_min)
    # p(eps) ~ eps exp(-eps) above eps_min has mean excess (eps_min + 2) / (eps_min + 1)
    expected = np.mean((eps_min + 2) / (eps_min + 1))
    assert np.mean(eps - eps_min) == pytest.approx(expected, rel=0.02)
    assert np.allclose(td.moon_distance(*state[:3]), r_c)


def test_the_sun_turns_once_a_synodic_month():
    assert fate.SUN_RATE * td.DAY_UNITS * 29.53059 == pytest.approx(2 * math.pi)


def test_the_tail_lies_behind_earth():
    earth = np.array([-td.MU, 0.0, 0.0])
    points = np.array([earth + [-1.0, 0, 0], earth + [1.0, 0, 0], earth + [0, 1.0, 0]]).T
    inside = fate.in_tail(points, np.zeros(3), np.zeros(3))     # Sun along +x at time zero
    assert list(inside) == [True, False, False]


def test_following_records_only_beyond_the_hill_sphere():
    rng = np.random.default_rng(6)
    state, _, _, r_c = fate.launch_escaping(3.0, 11.0, 60, rng)
    outcome, t, left_hill, records = fate.follow(state, r_c, t_max=20 * td.DAY_UNITS)
    assert set(np.unique(outcome)) <= {1, 2, 3, 4}
    for _, _, pos in records:
        assert np.all(td.moon_distance(*pos) > td.HILL)


def test_stored_cloud_grows_with_lifetime_and_its_fates_add_up():
    out = json.loads((fate.HERE / 'results' / 'escape_fate.json').read_text())
    assert out['schema'] == fate.SCHEMA
    masses = [c['cloud_kg_per_kg_s'] for c in out['clouds']]
    assert masses == sorted(masses)
    for c in out['clouds']:
        assert sum(c['fates'].values()) == pytest.approx(1.0)
        assert c['fates']['destroyed_by_sunlight'] > 0.9
        # the cloud cannot hold more than the lifetime's worth of escape
        assert c['cloud_kg_per_kg_s'] < c['lifetime_days'] * 86400
