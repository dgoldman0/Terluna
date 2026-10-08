import json
from types import SimpleNamespace

import numpy as np

from shared import constants as K
from shared.provenance import constants_changed
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from . import plane_motion as pm

PRODUCT = json.loads(pm.OUT.read_text())


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.plane-motion/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(pm.FILES)
    for name, value in producer['files'].items():
        assert pm.digest(pm.ROOT/name) == value, name
    assert producer['inputs']['results/frozen_rings.json'] == pm.digest(pm.FROZEN_OUT)
    if DEFAULT_KERNEL.exists():
        assert producer['inputs']['de440s.bsp'] == pm.digest(DEFAULT_KERNEL)
    assert not constants_changed(producer['constants'])


def test_the_tilt_puts_the_sunward_crossing_at_the_height_asked_for():
    sun = np.array([1.5e11, 0., 0.])
    pole = np.array([0., 0., 1.])
    assert np.isclose(pm.tilt_for(sun, pole, pole, 15000., 3000.), np.degrees(np.arcsin(.2)))
    tilted = np.array([0., -np.sin(.09), np.cos(.09)])
    beta = pm.tilt_for(sun, tilted, pole, 15000., 3000.)
    up = np.cos(np.radians(beta))*np.array([1., 0., 0.])+np.sin(np.radians(beta))*tilted
    assert np.isclose(15000.*up@pole, 3000.)


def test_a_ring_in_a_fixed_plane_crosses_the_three_planes_where_its_circle_does():
    radius, beta = 15e6, np.radians(10.)
    times = np.linspace(0., 3*86400., 4000)
    n = np.sqrt(K.MOON_GM/radius**3)
    up = np.array([np.cos(beta), 0., np.sin(beta)])
    side = np.array([0., -1., 0.])
    p = radius*(np.cos(n*times)[:, None]*up+np.sin(n*times)[:, None]*side)
    traj = np.concatenate([p, np.zeros_like(p)], axis=1)[:, None, :]
    env = SimpleNamespace(at=lambda t: dict(positions=dict(sun=np.array([1.5e11, 0., 0.]))))
    rows = pm.passages(env, times, traj, np.array([0., 0., 1.]), [4e6])[0]
    for row in rows:
        assert len(row['t']) >= 1
        x = row['x']
        expected = np.sin(beta)*np.sqrt(radius**2-x**2)
        assert np.allclose(row['y'], expected, rtol=1e-4)
        assert np.allclose(row['r'], radius, rtol=1e-6)


def test_a_common_shift_and_rotation_leave_no_differential_and_a_stretch_opens_every_overlap():
    heights = np.linspace(-6e6, 6e6, 9)
    xs = np.c_[-4e6*np.ones(9), np.zeros(9), 4e6*np.ones(9)]
    days = np.arange(5.)
    rigid = 1e3*days[:, None, None]+1e-4*days[:, None, None]*xs[None]
    _, _, residual, change = pm.split(heights, xs, rigid, 9286.)
    assert np.abs(residual).max() < 1e-6 and np.abs(change).max() < 1e-9
    stretch = 1e-3*days[:, None, None]*heights[None, :, None]*np.ones((1, 1, 3))
    _, _, _, change = pm.split(heights, xs, stretch, 9286.)
    assert np.allclose(change[-1], 4e-3*9286.)


def test_photon_steering_falls_short_of_the_stretch_of_the_strip_pattern():
    for row in PRODUCT['cases']:
        assert row['common']['oversizing_km_per_edge'] >= 0.
        assert row['differential']['overlap_change_per_step_m']['opening_max'] >= 0.
        if 'photon_steering' in row:
            assert row['photon_steering']['crossing_height_km_per_day'] > 0.


def test_the_common_shift_follows_the_suns_excursion_from_the_moons_orbit_plane():
    for row in PRODUCT['cases']:
        if row['plane'] == 'moon_orbit':
            amplitude = row['radius_km']*np.sin(np.radians(5.145))
            assert .7*amplitude < row['common']['centred_oversizing_km_per_edge'] < 1.2*amplitude


def test_the_strips_tear_apart_at_15000_km_and_hold_best_near_20000_km():
    by = {(r['plane'], r['radius_km']): r for r in PRODUCT['cases']}
    assert max(by[('moon_orbit', 15000.)]['photon_steering']['needed_over_available_by_ring']) > 10.
    edges = [by[('moon_orbit', r)]['differential']['residual_km_max'] for r in (15000., 19000., 20000.)]
    assert edges[0] > edges[1] > edges[2]
    for r in (21000., 22000.):
        assert by[('moon_orbit', r)]['differential']['overlap_change_per_step_m']['opening_max'] > 2000.
    middle = by[('moon_orbit', 20000.)]['photon_steering']['needed_over_available_by_ring'][2:7]
    assert max(middle) < 1.
