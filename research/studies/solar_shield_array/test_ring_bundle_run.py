import json

from shared.provenance import constants_changed
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from .ring_bundle_run import ANALYSIS, FILES, OUT, ROOT, SEED, digest

PRODUCT = json.loads(OUT.read_text())


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.ring-bundle/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(FILES+ANALYSIS)
    for name, value in producer['files'].items():
        assert digest(ROOT/name) == value, name
    assert producer['inputs']['results/joint_seed.json'] == digest(SEED)
    if DEFAULT_KERNEL.exists():
        assert producer['inputs']['de440s.bsp'] == digest(DEFAULT_KERNEL)
    assert not constants_changed(producer['constants'])


def test_free_rings_need_keeping_and_a_lone_ring_keeps_its_shingle():
    assert PRODUCT['segment']['clearance']['min_m'] < 150.
    assert PRODUCT['continuous']['clearance']['min_m'] < 150.
    lone = PRODUCT['lone_ring_gravity']
    assert lone['hours'] >= 120.
    assert 200. < lone['shingle_normal_range_m'][0] <= lone['shingle_normal_range_m'][1] < 250.
    loss = PRODUCT['starting_shadow_loss']
    assert loss['top_ring'] < loss['interior_rings']


def test_free_continuous_rings_fail_first_at_an_edge_and_keep_overlapping_along_each_ring():
    part = PRODUCT['continuous']
    edges = {0, PRODUCT['design']['rings']-1}
    pair = part['clearance']['first_below_150_m_pair']
    assert {pair[0], pair[2]} & edges
    assert part['interior_clearance']['first_below_150_m_h'] > part['clearance']['first_below_150_m_h']
    for name in ('segment', 'continuous'):
        assert PRODUCT[name]['service']['along_ring_overlap_min_m'] > 0


def test_the_tide_alone_separates_rings_slowly_in_shape_and_light_separates_them_fast():
    gravity = PRODUCT['continuous_gravity_only']['element_drift']
    light = PRODUCT['continuous']['element_drift']
    assert PRODUCT['continuous_gravity_only']['orbits'] >= 15
    assert gravity['interior_e_change_median_m'] < 50.
    assert gravity['interior_e_change_max_m'] < 150. and max(gravity['edge_e_change_max_m']) < 150.
    assert light['interior_e_change_median_m'] > 10*gravity['interior_e_change_median_m']


def test_the_tide_turns_the_stack_together_and_only_breathes_it():
    gravity = PRODUCT['continuous_gravity_only']
    fixed = gravity['element_drift']['plane_departure_final_m']
    node = gravity['node_frame']
    assert node['node_turn_deg'][-1] > 15.
    assert min(fixed) > 3000. and max(fixed)-min(fixed) < 50.
    assert max(node['departure_final_m']) < 200. < 500. < min(node['departure_max_m'])


def test_along_ring_overlap_compares_only_consecutive_tiles():
    import numpy as np
    from shared import constants as K
    from protection.dynamics import ring_bundle as rb
    from .ring_bundle_run import overlaps
    # One ring carrying two separate groups of tiles, as continuous rings carry copies around each representative.
    a, pitch, tilt = 15.0e6, 8.5e3, np.radians(1.5)
    slot = np.r_[np.arange(-3, 4), np.arange(10, 14)]
    phase = slot*pitch/a
    speed = np.sqrt(K.MOON_GM/a)
    states = np.stack([np.r_[a*np.cos(p), a*np.sin(p), 0., -speed*np.sin(p), speed*np.cos(p), 0.] for p in phase])
    along, across = overlaps(states, np.zeros(len(slot), int), slot, np.array([1., 0., 0.]), tilt,
                             np.ones(len(slot), bool))
    assert np.isclose(along, rb.CLEAR*np.cos(tilt)-pitch, atol=5.)
    assert across == np.inf
