"""The aimed-projectiles study: its product regenerates and binds its inputs, its flights reduce to the closed forms
where they exist and meet their sources, and the shooter's and the dodging target's figures follow their mechanics."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from humanities.sport import envelope as env
from humanities.sport import state
from humanities.sport.projectiles import aimed as A
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(A.OUT.read_text())
FT = 0.3048


def test_the_projectiles_product_regenerates_and_binds_its_code_constants_and_inputs():
    assert PRODUCT == json.loads(json.dumps(A.results()))
    assert PRODUCT['schema'] == A.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16] == h, path
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_flights_reduce_to_the_vacuum_and_fall_at_their_terminal_speed():
    for p in A.PROJECTILES.values():
        for _, g, rho in A.AIRS:
            # without air, a level shot from h0 falls g t^2 / 2 at x = v0 t
            t, x, y, z, v = A.fly(p, g, 0.0, [0.0], x_max=5.0)
            i = int(np.nanargmax(x[0]))
            assert A.MUZZLE_M - z[0, i] == pytest.approx(0.5 * g * t[i] ** 2, rel=1e-6, abs=1e-9)
            assert x[0, i] == pytest.approx(p['m_s'] * t[i], rel=1e-9)
            # dropped from rest, the speed follows v_t tanh(g t / v_t), tending to v_t = sqrt(g / k), drag equal weight
            q = dict(p, m_s=1e-9)
            t, x, y, z, v = A.fly(q, g, rho, [-90.0], h0=1e5, dt=0.01, t_max=60.0)
            vt = math.sqrt(g / A.drag_k(p, rho))
            for i in (100, 1000, len(t) - 1):
                assert v[0, i] == pytest.approx(vt * math.tanh(g * t[i] / vt), rel=1e-5)


def test_the_small_bore_bullet_keeps_the_makers_table_and_the_clay_its_measured_throw():
    p = A.PROJECTILES['small-bore bullet']
    t, x, y, z, v = A.fly(p, state.G_E, state.RHO_EARTH, [0.0], x_max=100.0)
    for yards, ft_s in ((50, 977), (100, 908)):
        speed = A.at_distance(t, x, y, z, v, yards * 0.9144)[0, 3]
        assert speed / FT == pytest.approx(ft_s, rel=0.02), yards
    high = A.SKEET['high house, long']
    earth = A.clay_flight(state.G_E, state.RHO_EARTH, high['height_m'], high['pitch_rad'])
    assert min(high['measured_m']) * 0.97 < earth['range_m'] < max(high['measured_m']) * 1.03
    # a lunar clay thrown by the same machine climbs on its lift and goes more than twice as far
    moon = A.clay_flight(state.G_M, state.RHO_MOON, high['height_m'], high['pitch_rad'])
    assert moon['range_m'] > 2 * earth['range_m'] and moon['apex_m'] > 2 * earth['apex_m']


def test_a_lunar_shot_drops_a_sixth_as_far_over_the_same_short_flight():
    for name in ('arrow', 'small-bore bullet'):
        moon, earth = (PRODUCT['projectiles'][name][w]['table'] for w in ('moon', 'earth'))
        for m, e in zip(moon, earth):
            if m['flat']['drop_m'] and e['flat']['drop_m'] and m['flat']['flight_s'] < 0.6:
                # nearly the same flight time, so the drop scales with g
                assert m['flat']['drop_m'] / e['flat']['drop_m'] == pytest.approx(state.G_M / state.G_E, rel=0.12)


def test_recoil_keeps_momentum_and_rocks_the_shooter_six_times_as_far_on_the_moon():
    dancer = state.BODIES['dancer']
    moon, earth = A.recoil(dancer, state.G_M), A.recoil(dancer, state.G_E)
    for name, (kg, v, gas, gun) in A.SHOOTERS.items():
        J = kg * v + gas
        assert moon[name]['impulse_N_s'] == pytest.approx(J, abs=1e-3)
        assert moon[name]['speed_change_m_s'] == pytest.approx(J / (dancer.mass_kg + gun), abs=1e-4)
    shot = 'shotgun, 24 g trap load'
    assert moon[shot]['rock_back_deg'] / earth[shot]['rock_back_deg'] == pytest.approx(state.G_E / state.G_M, rel=0.01)


def test_earths_ankle_stiffness_stands_a_shooter_still_on_the_moon_but_not_on_earth():
    for body in state.BODIES.values():
        assert not A.sway(body, state.G_E)['stiffness_0.91_of_earth_mgh']['passively_stable']
        moon = A.sway(body, state.G_M)
        assert all(row['passively_stable'] for row in moon.values())
        assert moon['stiffness_0.91_of_earth_mgh']['share_of_toppling'] == pytest.approx(0.91 * state.G_E / state.G_M,
                                                                                         rel=1e-3)


def test_a_seen_paintball_is_dodged_only_from_far_and_never_by_a_player_in_the_air():
    for world in ('moon', 'earth'):
        dodge = PRODUCT['paintball_dodge'][world]['by_body']['dancer']
        for posture, row in dodge['postures'].items():
            esc = [v for v in row['escape_m'].values()]
            assert esc == sorted(esc), (world, posture)          # a longer flight never shortens the escape
            if posture == 'aloft':
                assert row['least_dodge_m'] is None and max(esc) == 0.0
    moon = PRODUCT['paintball_dodge']['moon']['by_body']['dancer']['postures']
    earth = PRODUCT['paintball_dodge']['earth']['by_body']['dancer']['postures']
    for posture in ('ready', 'upright', 'moving, foot down'):
        assert moon[posture]['least_dodge_m'] > earth[posture]['least_dodge_m']
    # ducking half a metre runs on gravity's clock
    duck = PRODUCT['paintball_dodge']
    assert duck['moon']['duck_by_gravity_s']['0.5_m'] == pytest.approx(env.dip_s(0.5, state.G_M), abs=0.01)
