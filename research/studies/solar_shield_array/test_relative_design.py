import json

from shared.provenance import constants_changed
from .relative_design import CLEARANCE, FILES, OUT
from .relative_diagnosis import PACKAGE, ROOT, SEED, digest

PRODUCT = json.loads(OUT.read_text())


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.relative-design/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(FILES)
    for name, value in producer['files'].items():
        assert digest(ROOT/name) == value, name
    assert producer['inputs']['results/joint_seed.json'] == digest(SEED)
    assert producer['inputs']['local_continuation/ephemeris.npz'] == digest(PACKAGE/'ephemeris.npz')
    assert not constants_changed(producer['constants'])


def test_every_finite_patch_meets_a_neighbour_somewhere_in_the_orbit():
    patches = [c for c in PRODUCT['cases'] if c['kind'] != 'ring_bundle']
    assert {c['kind'] for c in patches} == {'equal_energy_patch', 'drifting_patch'}
    assert {c['attitude'] for c in patches} == {'radial', 'sun'}
    assert all(c['min_surface_distance_m'] < 10. for c in patches)


def test_ring_bundles_meeting_both_conditions_stay_clear_and_overlap_through_service():
    chosen = [c for c in PRODUCT['cases'] if c['kind'] == 'ring_bundle' and c['orbits'] == 3]
    assert chosen
    for c in chosen:
        assert c['samples_below_clearance'] == 0 and c['min_surface_distance_m'] >= CLEARANCE
        assert c['adjacent_row_overlap_min_m'] > 0 and c['along_row_overlap_min_m'] > 0


def test_rings_need_time_shifted_starts_and_own_attitudes_in_ephemeris_gravity():
    replay = PRODUCT['ephemeris_ring_bundle']
    circles = [r for r in replay if not r['time_shifted_rings']]
    common = [r for r in replay if r['time_shifted_rings'] and r['attitude'] == 'common']
    own = [r for r in replay if r['time_shifted_rings'] and r['attitude'] == 'own']
    assert circles and common and own
    assert all(r['min_surface_distance_m'] < CLEARANCE for r in circles+common)
    assert all(r['samples_below_clearance'] == 0 and r['min_surface_distance_m'] >= CLEARANCE for r in own)
    assert all(r['hours'] >= 90. for r in own)
    assert all(r['centre_slot_offset_m'] == 0. for r in common+own)
