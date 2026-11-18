"""Tests of the waterways (waterways.py): the main course upstream, where water is put in, and lakes crossed straight."""

import json

import numpy as np
import pytest

from research.studies.atlas_round import waterways as W


def test_main_donor_is_the_one_carrying_most():
    # cells 1, 2 and 3 drain into 0; 4 drains into 2
    receiver = np.array([-1, 0, 0, 0, 2])
    q = np.array([60.0, 10.0, 30.0, 20.0, 25.0])
    donor = W.main_donors(receiver, q)
    assert donor.tolist() == [2, -1, 4, -1, -1]


def test_follow_goes_up_the_main_flow_and_stops_where_water_is_put_in():
    receiver = np.array([-1, 0, 1, 2, 3, 1])
    q = np.array([5000.0, 5000.0, 4800.0, 4700.0, 150.0, 200.0])
    donor = W.main_donors(receiver, q)
    added = q - np.bincount(receiver[receiver >= 0], weights=q[receiver >= 0], minlength=q.size)
    cells, put_in = W.follow(0, donor, q, added, head=100.0, put_in=1000.0)
    assert cells == [0, 1, 2, 3] and put_in                  # cell 3 adds 4,550 beyond its donor's 150: water put in
    cells, put_in = W.follow(0, donor, q, np.zeros(6), head=100.0, put_in=1000.0)
    assert cells == [0, 1, 2, 3, 4] and not put_in           # with nothing put in it runs to the headwaters
    cells, _ = W.follow(0, donor, q, np.zeros(6), head=180.0, put_in=1000.0)
    assert cells == [0, 1, 2, 3]                             # and stops where the flow falls under the head's


def test_lakes_are_crossed_straight():
    cells = [10, 11, 12, 13, 14, 15, 16]
    lakes = np.array([False, True, True, True, False, True, False])
    kept, across, n = W.straighten(cells, lakes)
    assert kept == [10, 11, 13, 14, 15, 16]
    assert across == [False, False, True, False, False, False]
    assert n == 2


@pytest.mark.skipif(not W.OUT.exists(), reason='the waterways have not been traced')
def test_record_holds_named_waterways_largest_first():
    r = json.loads(W.OUT.read_text())
    flows = [w['flow_mouth_m3s'] for w in r['waterways']]
    assert flows == sorted(flows, reverse=True) and flows[-1] >= W.MOUTH_M3S
    assert len({w['anchor'] for w in r['waterways']}) == len(r['waterways'])
    names = [w['name'] for w in r['waterways']]
    assert all(names) and len(set(names)) == len(names)
    for w in r['waterways']:
        assert w['lake_km'] + w['channel_km'] in (w['length_km'] - 1, w['length_km'], w['length_km'] + 1)
        assert w['fed_by'] == 'headwaters' or w['flow_head_m3s'] >= W.PUT_IN_M3S
