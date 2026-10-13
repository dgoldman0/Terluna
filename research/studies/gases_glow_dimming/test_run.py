"""The gases, the glow and the dimming in the GCM: the products' binding, their regeneration and the arithmetic."""
import hashlib
import json
import math

import pytest

from research.studies.gases_glow_dimming import glow_absorption as A
from research.studies.gases_glow_dimming import run as R


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


@pytest.fixture(scope='module')
def absorption():
    return json.loads(A.OUT.read_text())


def close(a, b, rel=1e-2, abs_tol=1e-9):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_tol)


def gcm_inputs_present():
    paths = [R.CLIMATOLOGY.format(r) for r, _ in [R.CONTROL, *R.RUNS.values()]]
    paths += [R.PROGRESS.format(r) for r, _ in [R.CONTROL, *R.RUNS.values()]] + [R.EXTRACT]
    return all((R.ROOT / p).is_file() for p in paths)


def test_products_bind_their_producers_and_inputs(product, absorption):
    assert product['schema'] == R.SCHEMA and absorption['schema'] == A.SCHEMA
    for p in (product, absorption):
        for path, h in {**p['producer']['files'], **p['producer']['inputs']}.items():
            assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == h, path
    assert set(product['producer']['inputs']) == set(R.INPUTS)
    assert set(absorption['producer']['inputs']) == {str(p.relative_to(A.ROOT)) for p in A.input_files()}


def test_a_fresh_run_matches_the_product(product):
    if not gcm_inputs_present():
        pytest.skip('the GCM climatologies, run records or extract have not been generated')
    assert json.loads(R.dumps(R.results())) == product


def test_a_fresh_absorption_check_matches_its_product(absorption):
    if not (A.ROOT / 'atmosphere' / 'radiative_convective' / 'cache' / 'ck').is_dir():
        pytest.skip('the correlated-k tables have not been built')
    assert json.loads(A.dumps(A.results())) == absorption


def test_the_shares_of_infrared_from_space_add_to_one(absorption):
    """Each source's share reaching the ground and its shares absorbed in PlaSim's ten layers make up all of it; the
    colder the films, the less of their glow reaches the ground, and the returned infrared reaches it more than the
    films' own glow does."""
    for col in [*absorption['columns'].values(), absorption['at_design']]:
        sources = [col['returned_infrared'], *col['glow'].values()]
        for s in sources:
            assert len(s['by_plasim_layer']) == 10
            assert close(s['ground'] + sum(s['by_plasim_layer']), 1.0, abs_tol=1e-3)
        grounds = [col['glow'][f'{t:.0f}_k']['ground'] for t in A.GLOW_K]
        assert all(a < b for a, b in zip(grounds, grounds[1:]))
        assert col['returned_infrared']['ground'] > col['glow']['190_k']['ground']


def test_the_runs_are_the_runs_described(product):
    s = product['setup']
    assert s['co2_equivalent_ppm'] == 755.0
    assert close(s['calibration_check_first_ten_means_w_m2'], R.CALIBRATION['forcing_w_m2_first_ten_means'])
    # 755 ppm warms PlaSim's air as much as the middle pair's all-sky forcing, at PlaSim's own forcing per e-fold
    per_efold = R.CALIBRATION['forcing_w_m2_first_mean'] / math.log(R.CALIBRATION['co2_ppm'] / 400.0)
    assert close(per_efold * math.log(755.0 / 400.0), s['gases_forcing_allsky_w_m2'], rel=0.01)
    runs = product['runs']
    assert close(runs['gases']['added_dimming_percent'], s['gases_dimming_percent'], rel=0.005)
    assert close(runs['gases_and_glow']['added_dimming_percent'],
                 s['gases_dimming_percent'] + s['glow_watt_for_watt_dimming_percent'], rel=0.005)
    assert close(runs['gases_and_glow_held']['added_dimming_percent'],
                 s['gases_dimming_percent'] + s['glow_held_dimming_percent'], rel=0.005)
    day, night = R.GLOW['top_infrared_w_m2'].values()
    assert close((day + night) / 2, s['glow_infrared_w_m2'], rel=0.001)
    for name, r in runs.items():
        cfg = r['series']['configuration']
        assert cfg['co2_ppm'] == 755.0
        assert (cfg['top_infrared_w_m2'] is None) == (name == 'gases')
        assert cfg['initial_state']['path'].endswith('A28_dim5_moon/model/MOST_REST.00029')


def test_changes_are_run_less_control_and_growth_follows_evaporation(product):
    light = json.loads((R.ROOT / 'research/studies/light_limitation/results/light_limitation.json').read_text())
    for r in product['runs'].values():
        m = r['means']
        for v in (m['surface_k'], m['rain_mm_day']['land'], m['land_bands']['0-15']['evaporation_mm_day']):
            assert close(v['change'], v['run'] - v['control'], rel=1e-2, abs_tol=2e-3)
        g = r['rain_fed_growth']
        for e in (0, 1):
            total = sum(light['land']['bands'][k]['rain_fed_gt_c'][e]
                        * m['land_bands'][k]['evaporation_mm_day']['percent'] / 100.0 for k in light['land']['bands'])
            assert close(g['change_gt_c'][e], total, rel=1e-2)
            assert close(g['share_of_living_moon_percent'][e], 100.0 * total / light['ledger']['total_gt_c'][e],
                         rel=1e-2)
        f = r['first_means']
        assert close(f['air_radiation'], f['toa_net'] - f['surface_radiation'], rel=1e-2, abs_tol=2e-3)


def test_the_gases_dimming_holds_the_climate_and_the_watt_for_watt_glow_cools_it(product):
    runs = product['runs']
    assert abs(runs['gases']['means']['surface_k']['change']) < 0.3
    assert runs['gases_and_glow']['means']['surface_k']['change'] < -1.0
    assert abs(runs['gases_and_glow_held']['means']['surface_k']['change']) < 0.2
    for r in runs.values():
        assert r['means']['rain_mm_day']['all']['percent'] < 0.0
    # held at the ground, the glow warms the upper air (sigma below 0.25) and leaves the lowest air as it was
    held = runs['gases_and_glow_held']['layers']
    upper = [t for s, t in zip(held['sigma'], held['temperature_change_k']) if s < 0.25]
    lower = [t for s, t in zip(held['sigma'], held['temperature_change_k']) if s > 0.8]
    assert min(upper) > 0.5 and max(abs(t) for t in lower) < 0.2
    for name in ('gases_and_glow', 'gases_and_glow_held'):
        assert runs[name]['first_means']['glow_share_reaching_ground'] < 0.1


def test_the_held_point_lies_at_the_design_temperature(product):
    g, runs = product['glow'], product['runs']
    ww, held = runs['gases_and_glow'], runs['gases_and_glow_held']
    t = lambda r: r['means']['surface_k']['run']
    assert close(g['response_k_per_percent'],
                 (t(held) - t(ww)) / (g['watt_for_watt_dimming_percent'] - g['held_run_dimming_percent']), rel=2e-3)
    h = g['held']
    assert abs(h['means']['surface_k']['change']) < 0.01
    assert close(h['efficacy'], h['glow_dimming_percent'] / g['watt_for_watt_dimming_percent'], rel=2e-3)
    assert close(h['plan_dimming_percent'], runs['gases']['added_dimming_percent'] + h['glow_dimming_percent'],
                 rel=2e-3)
    assert 0.5 < h['efficacy'] < 0.9
    assert close(h['glow_rain_percent']['all'],
                 h['means']['rain_mm_day']['all']['percent'] - runs['gases']['means']['rain_mm_day']['all']['percent'],
                 rel=1e-2)
    assert close(h['glow_rain_percent_per_w_m2'], h['glow_rain_percent']['all'] / g['infrared_w_m2'], rel=1e-2)


def test_the_pairs_scale_with_their_warming(product):
    ledgers = json.loads((R.ROOT / 'research/studies/domain_synthesis/results/ledgers.json').read_text())
    weighed = ledgers['gases']['pairs_weighed']
    middle = product['pairs']['middle']
    for key, p in product['pairs'].items():
        k = weighed[key]['allsky_w_m2'] / weighed['middle']['allsky_w_m2']
        assert close(p['rain_percent']['all'], k * middle['rain_percent']['all'], rel=2e-3)
    assert product['pairs']['low']['rain_percent']['all'] > middle['rain_percent']['all'] \
        > product['pairs']['earth_2019']['rain_percent']['all']
