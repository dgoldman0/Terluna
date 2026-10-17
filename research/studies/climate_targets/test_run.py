"""What each target climate takes: the product's binding and regeneration, and the arithmetic behind it."""
import hashlib
import json
import math

import pytest

from research.studies.climate_targets import run as R
from research.studies.domain_synthesis import ledgers as L


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


def extract_present():
    return (R.ROOT / R.EXTRACT).is_file()


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == R.SCHEMA
    for path, h in product['producer']['files'].items():
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == h, path
    assert set(product['producer']['inputs']) == set(R.INPUTS)
    for path, h in product['producer']['inputs'].items():
        if path == R.EXTRACT and not extract_present():
            continue
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == h, path


@pytest.mark.skipif(not extract_present(), reason='the GCM extract (research/runs, ignored) is not here')
def test_a_fresh_run_matches_the_product(product):
    assert json.loads(json.dumps(R.results())) == product


def test_the_held_pair_is_the_domain_synthesis(product):
    assert R.HELD_PAIR == f'{L.HELD_CH4_PPM:g}+{L.HELD_N2O_PPM:g}'
    assert product['heat_beyond_sunlight']['gases']['held_pair_ppm'] == R.HELD_PAIR


def test_the_test_settles_near_where_its_last_years_stand(product):
    s = product['settled']
    assert 0.0 < s['feedback_w_m2_k'][0] <= s['feedback_w_m2_k'][1]
    test = product['runs']['test']
    assert test['ground_k'] < s['balance_k'][0] <= s['balance_k'][1] < test['ground_k'] + 1.0
    # its clear sky already sends out more than a saturated column can, so its balance rests on dry air and cloud
    assert s['clear_sky_outgoing_over_the_saturated_column_w_m2'] > 0.0


def test_the_heat_beyond_sunlight_adds_up(product):
    h = product['heat_beyond_sunlight']
    absorbed = product['runs']['test']['top_of_air_w_m2']['absorbed']
    for e in (0, 1):
        total = sum(h[k]['absorbed_sunlight_w_m2'][e] for k in ('gases', 'glow', 'beamed_power'))
        assert math.isclose(h['answered_in_every_target']['absorbed_sunlight_w_m2'][e], total, rel_tol=1e-12)
        assert math.isclose(h['answered_in_every_target']['dimming_percent'][e], 100.0 * total / absorbed, rel_tol=1e-12)
    g = h['glow']
    assert math.isclose(g['absorbed_sunlight_w_m2'][0], g['efficacy'] * g['infrared_w_m2']
                        + (1.0 - 0.28) * g['sunlight_back_w_m2'], rel_tol=1e-3)
    w = h['water_top_up']
    year = R.K.JULIAN_YEAR_DAYS * R.K.JULIAN_DAY
    area = 4.0 * math.pi * R.K.MOON_RADIUS ** 2
    assert math.isclose(w['first_millennium_w_m2'][1],
                        w['first_millennium_central_km3_yr'][1] * 1e12 / year * w['fall_mj_per_kg'][1] * 1e6 / area,
                        rel_tol=1e-12)


def test_cooler_targets_take_more_dimming_and_leave_more_margin(product):
    t = {x['name']: x for x in product['targets']}
    warm, mid, cool = t['no further dimming'], t['300 K'], t['295 K']
    assert warm['dimming_percent'] == [0.0, 0.0]
    for e in (0, 1):
        assert 0.0 < mid['dimming_percent'][e] < cool['dimming_percent'][e]
        assert warm['margin_to_the_saturated_column_w_m2'][e] < mid['margin_to_the_saturated_column_w_m2'][e]
        assert mid['margin_to_the_saturated_column_w_m2'][e] < cool['margin_to_the_saturated_column_w_m2'][e]
    assert cool['sunlight_scale'][1] < mid['sunlight_scale'][1] < warm['sunlight_scale'][1] < R.DESIGN_SCALE
    # the interpolated climate lies between the two runs, and threshold shares stay with the runs
    a, b = product['runs']['design'], product['runs']['test']
    for k in ('land_desert', 'wet_bulb_median_c', 'sea_median_c', 'land_rain_mm_day'):
        lo, hi = sorted((a[k], b[k]))
        assert lo <= mid['climate'][k] <= hi, k
    assert not any(k in mid['climate'] for k in R.RUNS_ONLY)


def test_the_extract_agrees_with_the_climate_readme(product):
    """The figures the GCM README states for the settled test (Where the corrected Moon settles)."""
    test, design = product['runs']['test'], product['runs']['design']
    assert round(test['ground_k'], 1) == 306.5 and round(design['ground_k'], 1) == 295.0
    assert round(100 * test['land_desert'], 1) == 18.9 and round(100 * design['land_desert'], 1) == 37.6
    assert round(test['wet_bulb_median_c'], 1) == 25.9 and round(test['wet_bulb_p99_c'], 1) == 31.0
    assert round(100 * test['wet_bulb_above_28'], 1) == 17.6
    assert round(test['sea_median_c'], 1) == 34.0
