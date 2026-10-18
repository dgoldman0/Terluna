"""The waters' oxygen: the products' binding, the saturation, the lakes' selection and what the runs must keep."""
import hashlib
import json

import numpy as np
import pytest

from biosphere.ecology import lake as L
from research.studies.lake_oxygen import lakes as K
from research.studies.lake_oxygen import run as R


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


@pytest.fixture(scope='module')
def chosen():
    return json.loads(K.OUT.read_text())


def present(path):
    return (R.ROOT / path).is_file()


def digest(path):
    return hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest()


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == R.SCHEMA
    assert not product['setting']['quick']
    for path, h in product['producer']['files'].items():
        assert digest(path) == h, path
    for path, h in product['producer']['inputs'].items():
        if path.startswith('research/runs/') and not present(path):
            continue                                     # the GCM's forcing lives on the research drive
        assert digest(path) == h, path


def test_lakes_product_binds_its_code_and_inputs(chosen):
    assert chosen['schema'] == K.SCHEMA
    for path, h in chosen['producer']['files'].items():
        assert digest(path) == h, path
    for path, h in chosen['producer']['inputs'].items():
        if path.startswith('geography/inputs/') and not present(path):
            continue                                     # LOLA and the gazetteer are fetched, not committed
        assert digest(path) == h, path


def test_the_lakes_reproduce_the_drainage_summary_and_lie_near_its_median(chosen):
    summary = json.loads(K.DRAINAGE_SUMMARY.read_text())['lakes']
    assert chosen['lakes']['bodies'] == summary['bodies']
    assert chosen['lakes']['median_local_depth_m'] == pytest.approx(summary['median_depth_m'], abs=1.0)
    for lake in chosen['chosen']:
        lo, hi = (float(x) for x in lake['band'].split('-'))
        assert lo <= abs(lake['lat']) < hi
        assert K.DEPTH_WINDOW[0] <= lake['deepest_m'] / chosen['lakes']['median_local_depth_m'] <= K.DEPTH_WINDOW[1]
        a = lake['hypsography']['area_fraction']
        assert a[0] == 1.0 and a[-1] == 0.0 and all(x >= y for x, y in zip(a, a[1:]))
        assert lake['crater'] is not None


@pytest.mark.skipif(not present('research/runs/lake_oxygen/forcing.json'), reason='the GCM forcing is not here')
def test_the_forcing_comes_from_the_current_extract(product):
    meta = json.loads((R.ROOT / 'research/runs/lake_oxygen/forcing.json').read_text())
    assert meta['schema'] == 'terluna.research.lake-oxygen-forcing/1'
    for path, h in meta['producer']['files'].items():
        assert digest(path) == h, path
    assert meta['producer']['read']['research/studies/lake_oxygen/results/lakes.json'] == digest(K.OUT)


def test_saturation_follows_the_fit_at_each_target(product):
    sat = product['saturation']
    checks = sat['checks']
    assert checks['earth_fresh_20c'] == pytest.approx(checks['expected']['earth_fresh_20c'], abs=0.005)
    assert checks['earth_fresh_35c'] == pytest.approx(checks['expected']['earth_fresh_35c'], abs=0.005)
    p = R.DESIGN_PRESSURE_PA / 101325.0
    for name, entry in sat['targets'].items():
        expected = float(L.oxygen_saturation(entry['sea_median_c'], p, sat['salinity_g_kg'], sat['o2_fraction']))
        assert entry['mg_l'] == pytest.approx(expected, abs=0.002)
        assert entry['seawater_mg_l'] < entry['mg_l']
    hot, resort, lock = (sat['targets'][k]['mg_l'] for k in ('hot', 'resort', 'lock'))
    assert hot < resort < lock
    assert sat['o2_fraction'] == pytest.approx(21227.0 / 121590.0, rel=1e-9)


def test_the_targets_temperatures_are_climate_targets_own(product):
    ct = json.loads(R.TARGETS.read_text())
    for name, ct_name in R.TARGET_NAMES.items():
        tgt = next(t for t in ct['targets'] if t['name'] == ct_name)
        assert product['saturation']['targets'][name]['sea_median_c'] == pytest.approx(tgt['climate']['sea_median_c'],
                                                                                       abs=0.005)
        assert product['setting']['share_of_the_way'][name] == tgt['interpolated_share_of_the_way']


def test_every_run_conserves_heat_and_oxygen(product):
    for group, cases in product['runs'].items():
        for key, case in cases.items():
            assert abs(case['budget']['heat_relative_error']) < 1e-8, (group, key)
            assert case['budget']['oxygen_relative_error'] < 1e-8, (group, key)


def test_production_meets_each_scenario(product):
    names = [s.name for s in R.scenarios()]
    for key, case in product['runs']['shallow'].items():
        for name, made in zip(names, case['budget']['production_g_m2_day']):
            assert made == pytest.approx(R.TROPHIC[name.split('/')[0]]['gpp'], rel=0.02), (key, name)


def test_no_lake_freezes_and_the_deep_runs_settle(product):
    for group, cases in product['runs'].items():
        for key, case in cases.items():
            assert case['ice_max_m'] == 0.0, (group, key)
    for key, case in product['runs']['deep'].items():
        spin = case['settling']['spin_up']
        drift = abs(case['settling']['bottom_t_last_two_years'] - case['settling']['bottom_t_first_two_years'])
        assert spin['settled'] or drift < 0.05, key


def test_the_slab_follows_the_gcms_own_water(product):
    for check in product['checks']['slab']:
        for name in ('lunar_free_convection', 'earth_free_convection'):
            assert abs(check[name]['bias_k']) < 1.5, (check['site'], check['run'], name)
            assert check[name]['correlation'] > 0.85, (check['site'], check['run'], name)
        # free convection carried to lunar gravity leaves the water warmer than Earth's constant does
        assert check['lunar_free_convection']['model_mean_c'] > check['earth_free_convection']['model_mean_c']


def test_the_answers_are_the_central_runs(product):
    for target, a in product['answers'].items():
        for site, o in a['overturn'].items():
            case = product['runs']['deep'][f'{target}/{site}']
            assert o['complete_overturns'] == case['mixing']['complete_overturns']
        for site, entry in a['shallows_late_night'].items():
            case = product['runs']['shallow'][f'{target}/{site}']
            for trophic, v in entry.items():
                n = case['oxygen'][f'{trophic}/renewal/supply']['nights']['surface_minimum']
                assert v['median'] == n[1]
                assert v['bracket_of_medians'][0] <= v['median'] <= v['bracket_of_medians'][1]


def test_scenarios_carry_their_sources():
    s = {x.name: x for x in R.scenarios()}
    assert len(s) == len(R.TROPHIC) * len(R.TRANSFER) * len(R.KINETICS)
    eu = s['eutrophic/renewal/supply']
    assert eu.gpp == pytest.approx(2.0 * 32.0 / 12.0) and eu.attenuation == pytest.approx(1.7)
    assert eu.sediment + eu.deep_water == pytest.approx(0.90)                      # Muller et al.'s mean
    reservoir = s['new_reservoir/renewal/supply']
    assert reservoir.sediment + reservoir.deep_water == pytest.approx(1.31)        # and their highest
    assert s['oligotrophic/renewal/warm'].activation_ev == pytest.approx(0.62)
    assert R.SOLUTE['co2'] == pytest.approx(44.01 / 32.0 * 0.284e-3)
    assert np.isclose(R.kd(1.0), 1.7)


def test_each_windows_humidity_spans_its_opening_and_closing_records():
    # PlaSim writes the humidity of the step that opens each window, held through it: a window's mean is the mean of
    # its own record and the next, and the last window keeps its own
    rng = np.random.default_rng(3)
    steps = rng.random(1001)                              # humidity at each model step of 10 windows of 100 steps
    records = steps[:-1:100]                              # what the output holds: each window's opening step
    windows = R.window_humidity(records)
    assert windows[:-1] == pytest.approx(0.5 * (steps[0:-100:100] + steps[100::100])[:9])
    assert windows[-1] == records[-1]
    linear = np.linspace(0.010, 0.020, 1001)              # a steady trend: the mean of the ends is the window mean
    assert R.window_humidity(linear[:-1:100])[:-1] == pytest.approx(
        [linear[i * 100:(i + 1) * 100 + 1].mean() for i in range(9)], rel=1e-12)
