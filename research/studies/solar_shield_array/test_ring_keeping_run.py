import json

from shared.provenance import constants_changed
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from .ring_keeping_run import FILES, OUT, ROOT, SEED, digest

PRODUCT = json.loads(OUT.read_text())


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.ring-keeping/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(FILES)
    for name, value in producer['files'].items():
        assert digest(ROOT/name) == value, name
    assert producer['inputs']['results/joint_seed.json'] == digest(SEED)
    if DEFAULT_KERNEL.exists():
        assert producer['inputs']['de440s.bsp'] == digest(DEFAULT_KERNEL)
    assert not constants_changed(producer['constants'])


def test_kept_bundle_record_is_complete():
    keeping = PRODUCT['keeping']
    assert len(keeping['per_orbit_m_s']) == PRODUCT['design']['orbits']
    assert keeping['peak_acceleration_m_s2'] <= PRODUCT['design']['limit_m_s2']*(1+1e-9)
    assert len(keeping['mean_m_s_per_orbit_by_ring']) == PRODUCT['bundle']['rings']
    checks = PRODUCT['checks']
    assert checks['clearance']['samples'] > 0 and checks['service']['samples'] > 0


def test_interior_rings_stay_clear_and_every_contact_involves_an_edge_ring():
    checks = PRODUCT['checks']
    assert checks['interior_clearance']['samples_below_150_m'] == 0
    edges = {0, PRODUCT['bundle']['rings']-1}
    pair = checks['clearance']['first_below_150_m_pair']
    assert {pair[0], pair[2]} & edges
    assert checks['service']['along_ring_overlap_min_m'] > 0 and checks['service']['adjacent_ring_overlap_min_m'] > 0


def test_interior_keeping_cost_is_the_tides_and_grows_away_from_the_free_ring():
    import numpy as np
    gravity = PRODUCT['gravity_only']
    orbits = gravity['orbits']
    light = np.array(PRODUCT['keeping']['per_orbit_m_s'])[:orbits]
    tide = np.array(gravity['per_orbit_m_s'])
    rings = PRODUCT['bundle']['rings']
    mid = rings//2
    interior = [k for k in range(1, rings-1) if k != mid]
    assert np.allclose(light[:, interior], tide[:, interior], rtol=.05)
    assert abs(light[:, interior].mean()/tide[:, interior].mean()-1) < .02
    assert np.all(light[:, [0, rings-1]] > 1.2*tide[:, [0, rings-1]])
    below = tide[:, 1:mid]
    assert np.all(np.diff(below, axis=1) < 0)
    radial, along, normal = np.array(gravity['radial_along_normal_m_s_per_orbit_by_ring']).T
    assert np.all(normal[interior] > radial[interior]+along[interior])
