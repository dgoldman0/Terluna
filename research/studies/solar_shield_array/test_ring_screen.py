import json

from shared.provenance import constants_changed
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from .ring_screen import FILES, OUT, RADII_KM, ROOT, TILTS_DEG, digest

PRODUCT = json.loads(OUT.read_text())


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.ring-screen/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(FILES)
    for name, value in producer['files'].items():
        assert digest(ROOT/name) == value, name
    if DEFAULT_KERNEL.exists():
        assert producer['inputs']['de440s.bsp'] == digest(DEFAULT_KERNEL)
    assert not constants_changed(producer['constants'])


def test_every_ring_case_is_screened_against_both_planes():
    cases = PRODUCT['cases']
    assert len(cases) == 2*len(RADII_KM)*len(TILTS_DEG)
    assert {c['plane'] for c in cases} == {'ecliptic', 'moon_orbit'}
    assert all(c['crossings'] > 100 for c in cases)
    rows = PRODUCT['differential']['between_tilts']+PRODUCT['differential']['between_radii']
    assert {r['plane'] for r in rows} == {'ecliptic', 'moon_orbit'}
