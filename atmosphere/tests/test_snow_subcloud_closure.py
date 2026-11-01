"""The adopted sub-cloud closure: its tables reproduce the fifth test, its families join smoothly, its product is bound
and carries the ladder's lapse."""
import hashlib
import json
import math
from pathlib import Path

import pytest

from atmosphere.snow import subcloud_closure as cl

PRODUCT = json.loads(cl.PRODUCT.read_text())
FAMILIES = cl.load(PRODUCT)


def test_product_was_written_by_this_module_from_the_capped_references():
    assert PRODUCT['schema'] == cl.SCHEMA
    assert PRODUCT['producer_sha256'] == hashlib.sha256(Path(cl.__file__).read_bytes()).hexdigest()
    for name, digest in PRODUCT['inputs']['files_sha256'].items():
        assert hashlib.sha256((cl.HERE / name).read_bytes()).hexdigest() == digest, name
    for fam, (names, top) in cl.FAMILIES.items():
        assert PRODUCT['families'][fam] == dict(references=list(names), top_m=top)
        for name in names:
            assert PRODUCT['runs'][name]['top_m'] == top
            assert all(r['k_per_m'] > 0.0 for r in PRODUCT['tables'][name])


def test_tables_reproduce_the_fifth_tests_closure():
    for v in cl.test_vectors(PRODUCT):
        stored = next(x for x in PRODUCT['test_vectors'] if x['name'] == v['name'])
        assert v['D'] == pytest.approx(stored['D'], rel=1e-12)
        assert v['D'] == pytest.approx(v['D_fifth_test'], rel=0.025), v['name']


def test_matching_is_linear_between_neighbours_and_held_at_the_ends():
    nodes = [math.log(x) for x in (1.0, 3.0, 10.0)]
    assert list(cl.match(math.log(0.5), nodes)) == [1.0, 0.0, 0.0]
    assert list(cl.match(math.log(20.0), nodes)) == [0.0, 0.0, 1.0]
    w = cl.match(0.5 * (nodes[1] + nodes[2]), nodes)
    assert w.sum() == pytest.approx(1.0) and w[1] == pytest.approx(0.5) and w[2] == pytest.approx(0.5)


@pytest.mark.parametrize('edge', [cl.TOP_S, cl.TOP_D])
def test_the_families_join_without_a_step_at_their_tops(edge):
    sub = lambda z: 0.05
    # a cloud base off the tables' layer edges, where the flux line's layers change by one
    d = [cl.closure(FAMILIES, 4100.0, 0.3 / 3600.0, top, sub, step=50.0)['D'] for top in (edge - 1.0, edge + 1.0)]
    assert d[0] == pytest.approx(d[1], rel=1e-3)


def test_tops_beyond_the_families_take_the_nearest_family_by_depth_below_the_top():
    sub = lambda z: 0.05
    high = cl.closure(FAMILIES, 4000.0, 0.3 / 3600.0, 45000.0, sub, step=50.0)
    low = cl.closure(FAMILIES, 2000.0, 0.3 / 3600.0, 9000.0, sub, step=50.0)
    assert high['weight_on_D'] == 1.0 and low['weight_on_D'] == 0.0
    assert math.isfinite(high['D']) and math.isfinite(low['D']) and high['D'] > 0.0 and low['D'] > 0.0


def test_moister_air_sublimates_less_and_saturated_air_none():
    f = 0.3 / 3600.0
    dry = cl.closure(FAMILIES, 4000.0, f, 24000.0, lambda z: 0.08, step=50.0)['D']
    moist = cl.closure(FAMILIES, 4000.0, f, 24000.0, lambda z: 0.02, step=50.0)['D']
    assert dry > moist > 0.0
    assert cl.closure(FAMILIES, 4000.0, f, 24000.0, lambda z: 0.0, step=50.0)['D'] == 0.0


def test_the_product_carries_the_ladders_lapse_through_both_families_tops():
    lapse = PRODUCT['inputs']['column']['lapse_k_per_m']
    for fam, top_c in (('S', -15.0), ('D', -30.0)):
        assert PRODUCT['families'][fam]['top_m'] * lapse == pytest.approx(-top_c, abs=0.01)
