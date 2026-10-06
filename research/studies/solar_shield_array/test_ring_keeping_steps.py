import json

from shared.provenance import constants_changed
from .ring_keeping_steps import FILES, KEPT, OUT, ROOT, digest

PRODUCT = json.loads(OUT.read_text())


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.ring-keeping-steps/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(FILES)
    for name, value in producer['files'].items():
        assert digest(ROOT/name) == value, name
    assert producer['inputs']['results/ring_keeping.json'] == digest(KEPT)
    assert not constants_changed(producer['constants'])


def test_every_kept_orbit_has_matched_phase_steps_for_all_pairs():
    rows = PRODUCT['orbits']
    kept = json.loads(KEPT.read_text())
    assert len(rows) == kept['design']['orbits']
    assert all(len(r['steps_m']) == kept['bundle']['rings']-1 for r in rows)


def test_interior_steps_stay_near_the_design_and_edge_steps_shrink():
    rows = PRODUCT['orbits']
    step = PRODUCT['designed_step_m']
    assert all(.85*step < r['interior_min_m'] <= r['interior_max_m'] < 1.15*step for r in rows)
    assert min(r['top_edge_m'] for r in rows) < .25*step
