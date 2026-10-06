import json

from shared.provenance import constants_changed
from .zoned_aperture import FILES, HERE, HOLDING, OUT, ROOT, digest

PRODUCT = json.loads(OUT.read_text())
CASES = {c['label']: c for c in PRODUCT['cases']}


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.zoned-aperture/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(FILES)
    for name, value in producer['files'].items():
        assert digest(ROOT/name) == value, name
    assert producer['inputs']['results/holding.json'] == digest(HOLDING)
    assert producer['inputs']['scenario.json'] == digest(HERE/'scenario.json')
    assert not constants_changed(producer['constants'])


def test_the_published_quadrature_is_reproduced():
    published = json.loads(HOLDING.read_text())['aperture_quadrature']['variable_distance']
    check = CASES['published_check']
    assert abs(check['mean_power_TW']/published['mean_power_TW']-1) < 1e-3
    assert abs(check['propellant_kg_s']/published['propellant_kg_s']-1) < 1e-3


def test_lighter_annulus_lowers_holding_power_and_the_core_share_stays_fixed():
    order = ['visible_passing_annulus_50g', 'annulus_20g', 'annulus_10g', 'annulus_5g']
    power = [CASES[k]['mean_power_TW'] for k in order]
    assert all(a > b for a, b in zip(power, power[1:]))
    cores = {round(CASES[k]['mean_power_core_TW'], 6) for k in order}
    assert len(cores) == 1
    assert CASES['core_26g_annulus_5g']['mean_power_TW'] < CASES['annulus_5g']['mean_power_TW']
