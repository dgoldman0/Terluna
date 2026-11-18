"""Tests of the design surface at 64 px/deg (fine_surface.py): the heights keep each coarse cell's mean and stay smooth
across cell edges, and the water follows the fine heights at the design's levels."""

import json

import numpy as np
import pytest

from research.studies.atlas_round import fine_surface as FS


def test_lift_keeps_every_block_mean_without_steps():
    rng = np.random.default_rng(1)
    coarse = np.cumsum(np.cumsum(rng.normal(0, 30, (24, 48)), axis=0), axis=1).astype(np.float32)
    fine, residual = FS.lift_to_means(coarse)
    assert np.abs(FS.block_mean(fine) - coarse).max() < 1e-3
    assert residual < 0.05 * np.abs(coarse).max()
    # no step at block edges larger than the field's own slope would give: compare jumps across and inside blocks
    jumps = np.abs(np.diff(fine, axis=1))
    across = jumps[:, FS.K - 1::FS.K][:, :-1].mean()
    inside = np.delete(jumps, np.s_[FS.K - 1::FS.K], axis=1).mean()
    assert across < 1.5 * inside


def test_fine_heights_keep_the_fine_shape_and_the_design_means():
    rng = np.random.default_rng(2)
    detail = rng.normal(0, 50, (24, 48)).astype(np.float32) + 5000.0     # another datum entirely
    # a gentle offset, as the geoid and the water's load are: a few metres per fine cell at most
    offset = (np.linspace(-30, 30, 12)[None, :] + np.linspace(0, 20, 6)[:, None]).astype(np.float32)
    design = FS.block_mean(detail) - 5000.0 + offset                      # the design: the same ground, a smooth offset
    fine, _ = FS.fine_heights(design, detail.copy())
    assert np.abs(FS.block_mean(fine) - design).max() < 1e-3
    d = detail - FS.expand(FS.block_mean(detail))
    s = fine - FS.expand(FS.block_mean(fine))
    assert np.corrcoef(d.ravel(), s.ravel())[0, 1] > 0.98


def test_a_ramp_stays_a_ramp():
    # the shape is kept whole: a smooth slope gets no sawtooth at the coarse cells' edges
    x = np.arange(48, dtype=np.float32)
    ramp = np.repeat((3.0 * x)[None, :], 24, axis=0)
    design = FS.block_mean(ramp) + 100.0
    fine, _ = FS.fine_heights(design, ramp.copy())
    assert np.allclose(fine, ramp + 100.0, atol=1e-3)


def test_fine_water_follows_the_heights_at_the_level():
    # one water cell (level 0) beside land cells; fine heights: a cove below the level joined to the water, a pond
    # below the level cut off by a ridge, and an islet above the level inside the water cell
    water16 = np.zeros((3, 3), bool)
    water16[1, 0] = True
    level16 = np.full((3, 3), np.nan, np.float32)
    level16[1, 0] = 0.0
    body16 = np.zeros((3, 3), np.int16)
    body16[1, 0] = 3
    h = np.full((12, 12), 50.0, np.float32)
    h[4:8, 0:4] = -10.0                    # the water cell's floor
    h[5, 2] = 5.0                          # an islet in it
    h[5:7, 4:6] = -3.0                     # a cove in the land cell beside it, joined to the water
    h[1:3, 9:11] = -3.0                    # a pond far from the water: a land cell not beside it
    h[9:11, 5:7] = -3.0                    # a pond beside the water's cell but cut off by dry ground
    code, level = FS.fine_water(h, water16, level16, body16)
    assert code[4, 0] == 3 and code[5, 2] == 0            # water, and the islet dry
    assert code[5, 4] == 3 and code[6, 5] == 3            # the cove joins the water
    assert code[1, 9] == 0                                # the far pond stays dry
    assert code[9, 5] == 0                                # the cut-off pond stays dry
    assert np.isnan(level[0, 0]) and level[4, 0] == 0.0


@pytest.mark.skipif(not FS.OUT.exists(), reason='the fine surface has not been built')
def test_fine_surface_record_checks():
    r = json.loads(FS.OUT.read_text())
    c = r['checks']
    assert c['block_mean_largest_off_m'] <= 0.5             # int16 half-metres round each fine cell by 0.25 m at most
    assert abs(c['water_share_fine'] - c['water_share_16ppd']) < 0.01
