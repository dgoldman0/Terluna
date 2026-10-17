"""What each target climate takes: the product's binding and regeneration, and the arithmetic behind it."""
import hashlib
import json
import math

import numpy as np
import pytest

from research.studies.climate_targets import extract as X
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


@pytest.mark.skipif(not extract_present(), reason='the GCM extract (research/runs, ignored) is not here')
def test_the_extract_comes_from_the_current_extract_code():
    """The extract records the code that made it and the files it read; the small ones are checked here, the output
    years by a fresh extract."""
    ex = json.loads((R.ROOT / R.EXTRACT).read_text())
    assert ex['schema'] == X.SCHEMA
    for path, h in ex['producer']['files'].items():
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == h, path
    for path, h in ex['producer']['read'].items():
        p = R.ROOT / path
        if p.is_file() and p.stat().st_size < 50_000_000:
            assert hashlib.sha256(p.read_bytes()).hexdigest() == h, path


def test_the_extract_weights_quantiles_by_area():
    rng = np.random.default_rng(0)
    x = rng.normal(size=20001)
    assert X.wquantile(x, np.ones_like(x), 0.5) == pytest.approx(np.median(x), abs=1e-12)
    assert X.wquantile(x, np.ones_like(x), 0.99) == pytest.approx(np.percentile(x, 99), abs=1e-3)
    # a cell twice another's area counts as two of it
    a, b = np.array([1.0, 2.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0])
    assert X.wquantile(b, np.array([1.0, 2.0, 1.0]), 0.5) == pytest.approx(X.wquantile(a, np.ones(4), 0.5))
    assert X.share(np.array([27.0, 29.0]), np.array([3.0, 1.0]), 28.0) == pytest.approx(0.25)


def test_the_coasts_join_the_seam_and_never_a_pole():
    sea = np.zeros((6, 8), bool)
    sea[-1, 3] = True                                    # a sea cell in the row nearest one pole
    near = X.beside(sea)
    assert not near[0].any() and near[-2:, 2:5].all()
    sea = np.zeros((6, 8), bool)
    sea[2, 0] = True                                     # and one on the longitude seam
    assert X.beside(sea)[1:4, -1].all()


def test_the_freezing_height_is_the_first_crossing_from_the_ground():
    ta = np.array([200.0, 230.0, 260.0, 270.0, 280.0, 290.0]).reshape(1, 6, 1, 1)     # levels from the top down
    z = np.array([30e3, 20e3, 10e3, 6e3, 3e3, 1e3]).reshape(1, 6, 1, 1)
    assert X.freezing_height(ta, z).item() == pytest.approx(3e3 + (280.0 - 273.15) / 10.0 * 3e3)


def test_the_held_pair_is_the_domain_synthesis(product):
    assert R.HELD_PAIR == f'{L.HELD_CH4_PPM:g}+{L.HELD_N2O_PPM:g}'
    assert product['heat_beyond_sunlight']['gases']['held_pair_ppm'] == R.HELD_PAIR


def test_the_test_settles_near_where_its_last_years_stand(product):
    s = product['settled']
    assert 0.0 < s['feedback_w_m2_k'][0] <= s['feedback_w_m2_k'][1]
    test = product['runs']['test']
    assert test['ground_k'] < s['balance_k'][0] <= s['balance_k'][1] < test['ground_k'] + 1.0
    # its unsaturated clear sky already sends out more than a saturated column can, and its cloud warms it on net
    assert s['clear_sky_outgoing_over_the_saturated_column_w_m2'] > 0.0
    assert s['cloud_effect_w_m2']['sunlight'] < 0.0 < s['cloud_effect_w_m2']['infrared']
    assert s['cloud_effect_w_m2']['sunlight'] + s['cloud_effect_w_m2']['infrared'] > 0.0


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
        for gap in ('gap_to_the_saturated_ceiling_w_m2', 'gap_with_the_added_heat_answered_w_m2',
                    'gap_with_the_added_heat_unanswered_w_m2'):
            assert warm[gap][e] < mid[gap][e] < cool[gap][e], gap
            assert warm[gap][0] <= warm[gap][1] and cool[gap][0] <= cool[gap][1], gap
    assert cool['sunlight_scale'][1] < mid['sunlight_scale'][1] < warm['sunlight_scale'][1] < R.DESIGN_SCALE
    for x in (warm, mid, cool):
        assert x['sunlight_scale'][1] < x['sunlight_scale_of_its_design_case'][0] or x is warm
    assert warm['sunlight_scale_of_its_design_case'] == [R.DESIGN_SCALE, R.DESIGN_SCALE]
    # each fit's own feedback with its own balance
    fits = product['settled']['fits']
    for x, target in ((mid, 300.0), (cool, 295.0)):
        assert x['absorbed_sunlight_to_shed_w_m2'] == sorted(f['feedback_w_m2_k'] * (f['balance_k'] - target) for f in fits)
    # the warm Moon's clear sky with the added heat unanswered would have to send out more than a saturated column can
    assert warm['gap_with_the_added_heat_unanswered_w_m2'][1] < 0.0 < warm['gap_with_the_added_heat_answered_w_m2'][0]
    # the interpolated climate lies between the two runs, and threshold shares stay with the runs
    a, b = product['runs']['design'], product['runs']['test']
    for k in ('land_desert', 'wet_bulb_median_c', 'sea_median_c', 'land_rain_mm_day'):
        lo, hi = sorted((a[k], b[k]))
        assert lo <= mid['climate'][k] <= hi, k
    assert not any(k in mid['climate'] for k in R.RUNS_ONLY)


def test_the_product_agrees_with_the_climate_readme(product):
    """The figures the GCM README states for the settled test (Where the corrected Moon settles), with the design
    case's, which it gives in brackets."""
    test, design = product['runs']['test'], product['runs']['design']
    assert round(test['ground_k'], 1) == 306.5 and round(design['ground_k'], 1) == 295.0
    assert round(100 * test['land_desert'], 1) == 18.9 and round(100 * design['land_desert'], 1) == 37.6
    # wet-bulb at 2 m (median, 99th percentile), air at 2 m, the same at 0.9 km, and the seas' median
    for run, figures in ((test, (26.9, 31.9, 33.5, 26.3, 31.3, 32.1, 34.6)),
                         (design, (16.7, 21.9, 20.4, 16.0, 21.3, 19.0, 24.2))):
        land = run['heat']['land']
        got = (land['at_2_m']['wet_bulb']['median'], land['at_2_m']['wet_bulb']['p99'], land['at_2_m']['air']['median'],
               land['wet_bulb']['median'], land['wet_bulb']['p99'], land['air']['median'], run['seas']['all']['median'])
        assert tuple(round(x, 1) for x in got) == figures
    land = test['heat']['land']
    assert round(100 * land['at_2_m']['above_28']) == 32 and round(100 * land['at_2_m']['above_31'], 1) == 3.3
    assert round(100 * land['above_28'], 1) == 23.6 and round(100 * land['above_31'], 1) == 1.6
    assert round(100 * land['counted_by_cell']['above_28'], 1) == 17.6
    assert round(test['freezing_level_over_land_km']) == 36 and round(design['freezing_level_over_land_km']) == 23
    # the place-by-place table: wet-bulb at 2 m (median, 99th percentile) and time above 28 C (at 2 m, 0.9 km up)
    table = {('bands', '0-15'): (27.8, 32.4, 46, 37), ('bands', '15-30'): (27.7, 32.0, 44, 34),
             ('bands', '30-45'): (27.2, 31.5, 35, 24), ('bands', '45-60'): (26.4, 30.8, 19, 12),
             ('bands', '60-90'): (25.6, 29.3, 6.2, 2.9), ('below_1_km',): (28.4, 32.6, 57, 46),
             ('from_4_km',): (25.4, 29.1, 5.7, 2.5), ('coasts',): (27.3, 32.1, 39, 30)}
    for path, (median, p99, at_2_m, at_level) in table.items():
        h = test['heat']
        for k in path:
            h = h[k]
        assert (round(h['at_2_m']['wet_bulb']['median'], 1), round(h['at_2_m']['wet_bulb']['p99'], 1)) == (median, p99)
        for share, figure in ((h['at_2_m']['above_28'], at_2_m), (h['above_28'], at_level)):
            assert round(100 * share, 1 if figure < 10 else 0) == figure, path
