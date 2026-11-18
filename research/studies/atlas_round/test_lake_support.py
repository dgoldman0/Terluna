"""Tests of the check of the design surface's lakes against its heights (lake_support.py)."""

import json

import numpy as np
import pytest

from research.studies.atlas_round import lake_support as LS


def hollow():
    """A hollow of two floors side by side under one rim: the left floor holds a lake at 100 m, the right is dry."""
    h = np.full((7, 12), 200.0, np.float32)
    h[2:5, 2:5] = 20.0                      # the left floor
    h[2:5, 5:8] = 30.0                      # the right floor, joined to the left with no sill between
    h[3, 9] = 50.0                          # a pit beyond the rim, not joined
    level = np.full(h.shape, np.nan, np.float32)
    level[2:5, 2:5] = 100.0
    return h, level


def test_a_lake_beside_lower_dry_ground_is_not_held():
    h, level = hollow()
    lake = np.isfinite(level)
    loose = LS.unheld(h, level, lake)
    assert loose[2:5, 4].all()                               # the lake's cells along the dry floor
    assert not loose[2:5, 2:4].any()                         # its other shores stand against the rim
    reach = LS.ground_below(h, level, lake)
    assert np.isfinite(reach[2:5, 5:8]).all() and (reach[2:5, 5:8] == 100.0).all()
    assert not np.isfinite(reach[3, 9])                      # the pit beyond the rim is not joined
    assert np.isfinite(reach).sum() == 9


def test_a_lake_held_by_its_rim_counts_nothing():
    h, level = hollow()
    h[2:5, 5:8] = 200.0                                      # the right floor filled in: rim all round
    lake = np.isfinite(level)
    assert not LS.unheld(h, level, lake).any()
    assert not np.isfinite(LS.ground_below(h, level, lake)).any()


def test_the_highest_level_reaches_where_two_lakes_do():
    h = np.full((5, 11), 200.0, np.float32)
    h[2, 1:10] = 10.0                                        # one trough
    level = np.full(h.shape, np.nan, np.float32)
    level[2, 1], level[2, 9] = 50.0, 80.0                    # a lake at each end
    reach = LS.ground_below(h, level, np.isfinite(level))
    assert (reach[2, 2:9] == 80.0).all()


def test_measure_counts_shores_bodies_and_ground():
    h, level = hollow()
    area = np.ones(h.shape)
    m, _ = LS.measure(h, level, np.isfinite(level), area)
    assert m['cells'] == 9 and m['shore_cells'] == 8 and m['shore_cells_not_held'] == 3
    assert m['bodies'] == 1 and m['bodies_with_such_a_shore'] == 1 and m['area_share_in_those_bodies'] == 1.0
    assert m['ground_below_and_joined_km2'] == 9 and m['ground_below_depth_m']['largest'] == 70


@pytest.mark.skipif(not LS.OUT.exists(), reason='the check has not been run')
def test_record_is_consistent():
    r = json.loads(LS.OUT.read_text())
    for kind in ('lakes', 'seas'):
        k = r[kind]
        assert 0 <= k['shore_cells_not_held'] <= k['shore_cells']
        assert k['bodies_with_such_a_shore'] <= k['bodies']
    p = r['catchment_pairs']
    assert p['in_one_catchment'] + p['in_two_catchments'] == p['pairs']
