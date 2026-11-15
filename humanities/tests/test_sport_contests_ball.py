"""The ball contests: a pass against an interceptor and a contested lob. Their screens hold the physics they rest on (the
ball's flight, the stances' holds, reach that grows with time), and the contests' product binds them and regenerates
one case of each."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from humanities.sport import contest_lob as LOB
from humanities.sport import contest_pass as PASS
from humanities.sport import contests, state, windows

ROOT = Path(__file__).resolve().parents[2]
MOON = (state.G_M, state.RHO_MOON)
EARTH = (state.G_E, state.RHO_EARTH)
SMALL_GRID = [(x, y) for x in (1.0, 5.0, 9.0) for y in (0.25, 1.25, 3.25)]


def test_the_pass_fan_stays_within_the_stance_and_more_head_start_never_shrinks_the_shadow():
    for g, rho in (MOON, EARTH):
        dancer = state.BODIES['dancer']
        for stance, (kind, lever) in PASS.STANCES.items():
            v0, p0, limit, meta = PASS._passes(dancer, kind, lever, g, rho, state.MU, 10.0)
            held = windows.release_limit(dancer, kind, PASS.BALL, g, state.MU,
                                         lever=lever if lever else dancer.release_above_com_m)['m_s']
            assert limit == pytest.approx(held)
            assert np.linalg.norm(v0, axis=1).max() == pytest.approx(max(held, PASS.MIN_SPEED_M_S))
        c = PASS.contest(g, rho, state.MU, 10.0, 'chest pass, feet under', grid=SMALL_GRID)
        names = ['react'] + [f'head_start_{h:g}_s' for h in PASS.HEAD_STARTS_S]
        for cell in c['cells']:
            shadowed = [not cell[n]['safe'] for n in names]
            # once a place is controlled, starting earlier keeps it controlled
            assert shadowed == sorted(shadowed), cell


def test_a_reacting_interceptor_beside_the_line_is_beaten_and_one_on_the_passer_with_a_head_start_is_not():
    for g, rho in (MOON, EARTH):
        c = PASS.contest(g, rho, state.MU, 10.0, 'chest pass, feet under', grid=[(5.0, 0.25), (0.5, 0.0)])
        mid, at_passer = c['cells']
        assert mid['react']['safe']
        # standing at the passer's hands with a second's start, every pass is touched before the catch
        assert not at_passer['head_start_1_s']['safe']


def test_a_lob_lands_where_it_is_aimed_and_hangs_longer_on_the_moon():
    for el in LOB.ELEVATIONS_DEG:
        moon, earth = LOB.lob(*MOON, el), LOB.lob(*EARTH, el)
        for shot in (moon, earth):
            end = shot['path'][-1]
            # the last sample before the floor sits no higher than the ball fell over the step before it
            assert end[0] == pytest.approx(LOB.LOB_M, abs=0.15)
            assert 0 <= end[2] <= shot['path'][-2][2] - end[2] + 1e-9
        assert moon['t'][-1] > 1.8 * earth['t'][-1]


def test_higher_lobs_gather_more_players_and_the_moon_gathers_more():
    for g, rho in (MOON, EARTH):
        areas = [LOB.contest(g, rho, state.MU, el, '3 m short, on the line')['crowd_area_m2']
                 for el in (60.0, 70.0)]
        assert areas[1] > areas[0] > 0
    moon = LOB.contest(*MOON, state.MU, 70.0, '3 m short, on the line')
    earth = LOB.contest(*EARTH, state.MU, 70.0, '3 m short, on the line')
    assert moon['crowd_area_m2'] > 2 * earth['crowd_area_m2']
    # the first touch on a high lob waits for the ball to come down into the jumpers' reach
    for c in (moon, earth):
        assert c['b']['touch_s'] == pytest.approx(c['earliest_touch_s'], abs=0.05)
        assert c['b']['slack_s'] >= 0.5


def test_the_contests_product_binds_its_code_and_regenerates_one_case_of_each_ball_contest():
    if not contests.OUT.exists():
        pytest.skip('results/contests.json is built by python -m humanities.sport.contests')
    product = json.loads(contests.OUT.read_text())
    assert product['schema'] == contests.SCHEMA
    for path, h in {**product['producer']['files'], **product['producer']['inputs']}.items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest()[:16] == h, path
    # the quickest case to rebuild: Earth's 6 m chest pass (a lunar case takes a minute or more)
    stored = next(c for c in product['pass_against_interceptor']['cases']
                  if c['world'] == 'earth' and c['grip'] == state.MU and c['distance_m'] == 6.0
                  and c['stance'] == 'chest pass, feet under' and 'interceptor' not in c)
    again = PASS._case(('earth', *EARTH, state.MU, 6.0, 'chest pass, feet under', 'dancer'))
    assert json.loads(json.dumps(again)) == stored
    stored = [c for c in product['contested_lob']['cases'] if c['world'] == 'earth' and c['grip'] == state.MU
              and c['b_place'] == '3 m beside']
    again, _ = LOB._case(('earth', *EARTH, state.MU, '3 m beside'))
    assert json.loads(json.dumps(again)) == stored
