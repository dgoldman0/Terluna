"""The parkour moves: their product regenerates and binds its inputs, the units are the surface moves' and join the
athletic vocabulary in one graph, and what they change against the floor follows the surface moves' mechanics."""
import hashlib
import json
import math
from pathlib import Path

import pytest

from humanities.sport import envelope as env
from humanities.sport import parkour as P
from humanities.sport import state, surfaces
from humanities.sport import vocabulary as V
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(P.OUT.read_text())
UNITS = {u['id']: u for u in P.load()['units']}


def test_the_parkour_product_binds_its_inputs_and_regenerates_with_one_route_rebuilt():
    assert PRODUCT['schema'] == P.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16] == h, path
    assert not constants_changed(PRODUCT['producer']['constants'])
    again = json.loads(json.dumps(dict(units=[dict(id=u['id'], by_body=P.figures(u)) for u in P.load()['units']],
                                       wall=P.wall_against_floor(), curved=P.curved_wall())))
    assert again['units'] == [dict(id=u['id'], by_body=u['by_body']) for u in PRODUCT['units']]
    assert again['wall'] == PRODUCT['wall_against_floor'] and again['curved'] == PRODUCT['curved_wall_least_m_s']
    case = ('a 1 m barrier across the way', 'moon', state.MU, 'dancer')
    stored = next(r for r in PRODUCT['routes'] if (r['arena'], r['world'], r['grip'], r['body']) == case)
    assert json.loads(json.dumps(P._route(case))) == stored


def test_the_units_are_the_surface_moves_and_evaluate_everywhere():
    assert tuple(UNITS) == P.UNITS == surfaces.PARKOUR_UNITS + ('rise_from_low',)
    for u in UNITS.values():
        assert u['lunar'].split(':')[0] in ('native', 'adapted', 'native in its hazard'), u['id']
        for body in state.BODIES.values():
            for _, g in state.WORLDS:
                for mu in state.GRIPS:
                    t = P.times(u, body, g, mu)
                    values = [t['duration_s'], t['tell_s'], t['commits_s'], t['recover_s'], *t['figures'].values()]
                    assert all(math.isfinite(x) and x >= 0 for x in values), (u['id'], body.key, g, mu)
                    assert t['figures'].get('possible', 1) == 1, (u['id'], body.key, g, mu)


def test_the_parkour_moves_join_the_athletic_vocabulary_in_one_graph():
    units, graph = P.joined_graph()
    names = {u['id'] for u in units}
    assert V.reachable(graph, 'ready_crouch') == names
    back = {k: [a for a in graph if k in graph[a]] for k in graph}
    assert V.reachable(back, 'ready_crouch') == names
    # a jump leads to a catch on a ledge, a hang to a climb or a swing, and low cover back to the crouch
    assert 'cat_leap' in graph['crouch_jump'] and {'climb_up', 'lache'} <= set(graph['cat_leap'])
    assert 'rise_from_low' in graph['slide'] and 'ready_crouch' not in graph['slide']


def test_on_the_moon_a_wall_turns_a_runner_where_the_floor_cannot():
    wall = PRODUCT['wall_against_floor']
    for mu in state.GRIPS:
        moon, earth = wall['moon'][str(mu)], wall['earth'][str(mu)]
        # the wall's push does not depend on the floor's grip; the floor's stance does, and on Earth it does far more
        assert moon['wall_change_m_s'] == earth['wall_change_m_s'] == wall['moon']['0.8']['wall_change_m_s']
        assert moon['floor_one_stance_m_s'] == pytest.approx(mu * state.G_M * env.STRIDE_S, abs=0.01)
        # the same push against a stance that scales with g: the Moon's advantage over Earth's is g_E / g_M
        assert moon['wall_over_floor'] / earth['wall_over_floor'] == pytest.approx(state.G_E / state.G_M, rel=0.03)
        assert moon['wall_over_floor'] > 8
    # the floor's cut through the wall's turn takes far longer than the wall's 0.38 s contact on the Moon
    assert wall['moon']['0.8']['floor_cut_total_s'] > 10 * wall['moon']['0.8']['contact_s']


def test_pulling_down_behind_cover_beats_gravity_and_the_rise_from_low_stops_on_gravitys_clock():
    dancer = state.BODIES['dancer']
    moon = P.times(UNITS['pull_down'], dancer, state.G_M, state.MU)['figures']
    earth = P.times(UNITS['pull_down'], dancer, state.G_E, state.MU)['figures']
    assert moon['saves_s'] > 0.4 > 0.1 > earth['saves_s'] > 0
    assert moon['gravity_alone_s'] == pytest.approx(env.dip_s(0.5, state.G_M), abs=1e-3)
    for _, g in state.WORLDS:
        ns = P.physics(dancer, g, state.MU)
        rise = P.times(UNITS['rise_from_low'], dancer, g, state.MU)['duration_s']
        assert rise == pytest.approx(ns['rise_s'](dancer.ready_com_m - ns['slide_com']))
