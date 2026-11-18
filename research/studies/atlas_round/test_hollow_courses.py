"""Checks of the courses joined through the corrected routing's lakes (hollow_courses.py): on the hand-worked hollow
of three lowest points and on random ground every course reaches a sea, the flow on dry ground is the routing's own,
and a lake's exit carries what the lake passes on."""

import numpy as np
import pytest

from research.studies.atlas_round import hollow_courses as HC
from research.studies.atlas_round import lake_rule as LR
from research.studies.atlas_round.test_lake_rule import U, random_case


def routed(h, outlet, lat, runoff, net, level, sources=()):
    s = LR.hollow_structure(h, outlet, lat, 1.0, level)
    terms = LR.hollow_terms(s, runoff, net)
    state = LR.hollow_state(s, terms, sources)
    r = LR.hollow_routing(s, terms, state)
    return s, terms, state, r, HC.connect(s, terms, state, r, h.shape)


def ends(receiver, outlet):
    """Where each cell's course ends, following the receivers."""
    at = receiver.copy()
    for _ in range(64):
        at = at[at]
    return at, outlet[at]


def profile(q):
    """test_lake_rule's row: a sea, a rim at 10, lowest points A at 2, B at 1 and C at 3 behind sills at 6 and 9,
    and a rim at 12; every cell under water loses 1 unit and q units are put in at A."""
    h = np.array([[-5.0, 10.0, 2.0, 6.0, 1.0, 9.0, 3.0, 12.0]], np.float32)
    outlet = np.zeros(h.shape, bool)
    outlet[0, 0] = True
    net = np.full(h.shape, -1.0, np.float32)
    return routed(h, outlet, np.zeros(1), np.zeros(h.shape, np.float32), net, np.zeros(h.shape), [(2, q * U)] if q else ())


def test_a_full_hollow_of_three_lowest_points_is_crossed_and_left_at_its_pass():
    s, _, _, r, c = profile(6.0)                                    # the whole hollow full, 1 unit over the outer pass
    at, sea = ends(c['receiver'], s['outlet'])
    assert sea.all() and (at == 0).all()                            # every course reaches the sea
    assert c['receiver'][2] == 1                                    # A lies by the outer pass: its water goes over it
    assert c['discharge'][2] / U == pytest.approx(1.0)              # and it carries what the lake passes on
    assert c['discharge'][1] / U == pytest.approx(1.0) and c['discharge'][0] / U == pytest.approx(1.0)
    assert c['receiver'][4] in (3, 2) and c['receiver'][6] in (5, 4)    # B's and C's water goes on towards A
    assert c['passes'][2] and not c['passes'][[4, 6]].any()         # B and C took in nothing of their own to pass


def test_a_lake_that_spills_into_a_hollow_not_full_ends_there():
    s, _, _, r, c = profile(3.5)                                    # A and B one lake to the sill at 9, 0.5 into C
    assert c['discharge'][6] / U == pytest.approx(0.5)              # C's lowest point receives the spill
    assert not c['passes'][6]                                       # and passes nothing on
    at, sea = ends(c['receiver'], s['outlet'])
    assert sea.all()                                                # the courses still lead to the sea
    assert c['discharge'][0] == 0.0                                 # which receives nothing


@pytest.mark.parametrize('seed', [1, 2, 3, 4, 5, 6])
def test_on_random_ground_courses_reach_a_sea_and_carry_the_routings_flow(seed):
    h, outlet, lat, runoff, net, level, sources = random_case(seed)
    s, terms, state, r, c = routed(h, outlet, lat, runoff, net, level, sources)
    flat_out = s['outlet']
    _, sea = ends(c['receiver'], flat_out)
    assert sea.all()
    dry = ~r['lake'] & ~flat_out & ~s['sink'] & ~c['crossing']
    assert c['discharge'][dry] == pytest.approx(r['discharge'][dry], abs=1e-9 * max(r['discharge'].max(), 1.0))
    assert c['discharge'][flat_out].sum() == pytest.approx(r['discharge'][flat_out].sum())
    # a lake that passes nothing on stops the drained area
    area_to_seas = c['drained_m2'][flat_out].sum()
    assert 0.0 < area_to_seas < s['area'][~flat_out].sum()          # some ground drains to lakes that keep their water
