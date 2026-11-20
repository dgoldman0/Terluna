"""Tests of look_tiles.py's parts that need no surface: the blocks a pack is worked in, the far average and the decoding."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import appearance  # noqa: E402
import look_tiles as LT  # noqa: E402


def pack(valid, r=1, c=2):
    return dict(r=r, c=c, valid=list(np.asarray(valid, bool).ravel()))


def test_a_pack_is_worked_in_blocks_of_two_by_two_tiles_that_hold_a_valid_tile():
    valid = np.zeros((4, 5), bool)
    valid[0, 1:4] = True          # the first pair of rows: columns 1 to 3, so the blocks at columns 0 and 2
    valid[1, 3] = True
    valid[3, 4] = True            # the second pair: the last column alone, a block one tile wide
    assert LT.blocks(pack(valid)) == [(4, 10, 2, 2), (4, 12, 2, 2), (6, 14, 2, 1)]


def test_a_block_with_no_valid_tile_is_left_out_and_every_valid_tile_lies_in_one_block():
    valid = np.zeros((4, 5), bool)
    valid[2:, :] = True
    got = LT.blocks(pack(valid, r=0, c=0))
    assert got == [(2, 0, 2, 2), (2, 2, 2, 2), (2, 4, 2, 1)]
    assert LT.blocks(pack(np.zeros((4, 5), bool))) == []
    covered = np.zeros((4, 5), int)
    for tr, tc, nr, nc in got:
        covered[tr:tr + nr, tc:tc + nc] += 1
    assert (covered[valid] == 1).all()


def test_the_far_mean_is_the_mean_of_each_four_by_four():
    a = np.arange(8 * 8 * 3, dtype=np.float32).reshape(8, 8, 3)
    m = LT.far_mean(a)
    assert m.shape == (2, 2, 3)
    assert np.allclose(m[0, 0], a[:4, :4].mean(axis=(0, 1)))
    assert np.allclose(m[1, 1], a[4:, 4:].mean(axis=(0, 1)))


def test_decoding_undoes_the_pages_encoding():
    linear = np.linspace(0.0, 1.0, 256, dtype=np.float32)
    encoded = np.round(appearance.srgb_encode(linear) * 255).astype(np.uint8)
    assert np.abs(LT.decode(encoded) - linear).max() < 0.012
