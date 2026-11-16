import numpy as np

from research.studies.drainage_network import stops as S


def test_runs_are_found_from_end_to_end():
    flags = np.array([True, True, False, True, False, False, True])
    assert S.runs_of(flags) == [(0, 2), (3, 4), (6, 7)]
    assert S.runs_of(np.zeros(3, bool)) == []
    # two hollows that meet on the path, each at its own level, stay apart
    levels = np.array([5.0, 5.0, 3.0, 3.0, 0.0, 1.0])
    assert S.runs_of(np.array([True, True, True, True, False, True]), levels) == [(0, 2), (2, 4), (5, 6)]


def test_a_hollow_is_the_low_ground_joined_to_its_seed_below_its_level():
    h = np.array([[9, 9, 9, 9, 9, 9],
                  [9, 2, 3, 9, 1, 9],
                  [9, 3, 4, 9, 1, 9],
                  [9, 9, 9, 9, 9, 9]], float)
    cells = S.component_below(h, 5.0, np.array([7]))           # the left hollow, from its deepest cell
    assert cells.reshape(h.shape)[1:3, 1:3].all() and cells.sum() == 4   # the right hollow stays apart


def test_each_cell_keeps_its_first_and_second_body_and_when_the_second_is_its_surface():
    ident = np.array([[1, 0, 2, 0],      # year 0: a sea over the first cell, another over the third
                      [1, 7, 2, 0],      # a hollow fills on an overflow
                      [3, 7, 2, 5],      # a rising sea covers the first cell; a new water over the last
                      [3, 0, 4, 5]])     # a fourth body over the third cell: more than two
    first, second, second_bits, extra = S.two_bodies(ident)
    assert first.tolist() == [1, 7, 2, 5]
    assert second.tolist() == [3, 0, 4, 0]
    assert second_bits.tolist() == [0b1100, 0, 0b1000, 0]
    assert extra == 0
    ident[1, 2] = 6
    assert S.two_bodies(ident)[3] == 1
