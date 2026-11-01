"""Exact-moment controls for the bounded collision projection diagnostic."""
import json

import numpy as np
import pytest

from atmosphere.snow import collision_projection as cp


def test_one_sticking_event_conserves_first_moment_and_adds_known_second_moment():
    # The physical product 11 lies between pivots 10 and 14: 3/4 and 1/4.
    row = cp.probe(np.array([1., 3., 10., 14., 30.]), 0, 2)
    assert row['number_change_per_event'] == -1.
    assert row['relative_first_moment_residual'] == 0.
    assert row['exact_second_moment_change_kg2'] == 20.
    assert row['projected_second_moment_change_kg2'] == 23.
    assert row['excess_second_moment_change_kg2'] == 3.
    assert row['upper_pivot_allocation_number_fraction'] == .25
    assert row['collector_jump_second_moment_ratio'] == 4.


def test_probe_declares_its_interior_small_increment_scope():
    with pytest.raises(ValueError, match='below the next'):
        cp.probe(np.array([3., 10., 12., 20.]), 0, 1)
    with pytest.raises(ValueError, match='interior'):
        cp.probe(np.array([1., 10., 20.]), 0, 2)


def test_saved_product_matches_source_and_moment_inequalities():
    result = json.loads(cp.OUT.read_text())
    assert result['producer'] == cp.identity()
    assert len(result['rows']) == 15
    for diameter in (.01, .1, .15, .6, 1.):
        rows = [r for r in result['rows'] if r['requested_collector_diameter_m'] == diameter]
        excess = [r['excess_second_moment_fraction'] for r in rows]
        assert excess[0] > excess[1] > excess[2] > 0.
        for row in rows:
            assert row['excess_second_moment_fraction'] < (row['geometric_mass_ratio']-1)/2
            assert row['number_change_per_event'] == pytest.approx(-1., abs=1e-14)
    assert result['checks']['maximum_excess_formula_relative_residual'] < 1e-11
    assert result['settings']['time_integration'] is False
