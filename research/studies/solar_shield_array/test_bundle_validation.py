import json

import numpy as np

from shared.provenance import constants_changed
from protection.dynamics import ring_bundle as rb
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from . import bundle_validation as bv

PRODUCT = json.loads(bv.OUT.read_text())


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.bundle-validation/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(bv.FILES)
    for name, value in producer['files'].items():
        assert bv.digest(bv.ROOT/name) == value, name
    assert producer['inputs']['results/photon_control.json'] == bv.digest(bv.PHOTON)
    if DEFAULT_KERNEL.exists():
        assert producer['inputs']['de440s.bsp'] == bv.digest(DEFAULT_KERNEL)
    assert not constants_changed(producer['constants'])


def test_a_ray_through_a_tile_is_covered_by_half_its_side_and_one_past_its_edge_is_not():
    tile = np.array([[15e6, 0., 0., 0., 1800., 0.]])
    frame = rb.tile_frames(tile, np.radians(1.5))[0]
    origin = np.array([[0., 0., 0.], [0., 0., 0.]])
    target = np.array([tile[0, :3], tile[0, :3]+5100.*frame[:, 1]])
    d = target-origin
    d /= np.linalg.norm(d, axis=1)[:, None]
    hits, margin = bv.coverage(origin, d, tile, np.radians(1.5))
    assert hits.tolist() == [1, 0]
    assert np.isclose(margin[0], rb.CLEAR/2, rtol=1e-6)


def test_the_patch_stays_clear_and_covered_with_every_tile_propagated():
    assert PRODUCT['interior_clearance']['min_m'] > 150.
    assert PRODUCT['service']['along_ring_overlap_min_m'] > 0. and PRODUCT['service']['adjacent_ring_overlap_min_m'] > 0.
    rays = PRODUCT['rays']
    assert rays['tested'] > 100000 and rays['uncovered'] == 0
    assert rays['margin_min_m'] > 0.


def test_holding_interior_tiles_to_their_slots_asks_little_of_photon_keeping():
    keeping = PRODUCT['keeping']
    assert keeping['interior_tiles_m_s2']['p95'] < .25*keeping['sail_acceleration_normal_incidence_m_s2']
    assert keeping['slot_offset_m']['interior_max'] < 200.
