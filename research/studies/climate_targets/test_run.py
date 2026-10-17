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
    for e in (0, 1):
        assert math.isclose(g['absorbed_sunlight_w_m2'][e], g['efficacy'][e] * g['infrared_w_m2']
                            + (1.0 - 0.28) * g['sunlight_back_w_m2'], rel_tol=1e-3)
    # the glow's two efficacies: the corrected radiation's, which the product finds, and the held branch of 13 October
    assert g['efficacy'] == sorted(g['efficacy_sources'].values())
    assert g['efficacy_sources']['corrected_radiation_17_october'] == product['glow_efficacy']['efficacy']
    held = json.loads((R.ROOT / R.GASES_GLOW).read_text())['glow']['held']['efficacy']
    assert g['efficacy_sources']['held_branch_13_october'] == held
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


def test_the_gap_with_the_heat_answered_spans_every_corner(product):
    """The gap's ends are the least and most load over the ends of the feedback, the glow's efficacy and the power."""
    h = product['heat_beyond_sunlight']
    clear = product['settled']['clear_sky_absorbed_w_m2']
    saturated = product['settled']['saturated_column_max_outgoing_w_m2']
    absorbed = product['runs']['test']['top_of_air_w_m2']['absorbed']
    gas = h['gases']['absorbed_sunlight_w_m2'][0]
    for t in product['targets']:
        gaps = [saturated - clear * (1.0 - (d + 100.0 * (gas + g + p) / absorbed) / 100.0) - out
                for d in t['dimming_percent'] for g in h['glow']['absorbed_sunlight_w_m2']
                for p, out in zip(h['beamed_power']['absorbed_sunlight_w_m2'], h['to_send_out_w_m2'])]
        assert t['gap_with_the_added_heat_answered_w_m2'] == pytest.approx([min(gaps), max(gaps)], abs=1e-9)


def test_the_glow_efficacy_follows_from_the_glow_pair(product):
    """The corrected radiation's efficacy, found as the gases, glow and dimming study found its 0.71."""
    g = product['glow_efficacy']
    d, t = g['added_dimming_percent'], g['ground_k']
    # the pair carries 13 October's added dimming on the design scale of 0.970
    assert round(d['watt_for_watt'], 2) == 4.22 and round(d['held'], 2) == 3.52
    response = (t['held'] - t['watt_for_watt']) / (d['watt_for_watt'] - d['held'])
    assert g['response_k_per_percent'] == pytest.approx(response, rel=1e-12)
    hold = d['held'] - (t['design'] - t['held']) / response
    assert g['dimming_that_holds_percent'] == pytest.approx(hold, rel=1e-12)
    gases = json.loads((R.ROOT / R.GASES_GLOW).read_text())['runs']['gases']['added_dimming_percent']
    assert g['gases_dimming_percent'] == gases
    glow_w_m2 = (hold - gases) / 100.0 * g['design_absorbed_sunlight_w_m2']
    assert g['efficacy'] == pytest.approx(glow_w_m2 / g['infrared_w_m2'], rel=1e-12)
    assert g['design_absorbed_sunlight_w_m2'] == product['runs']['design']['top_of_air_w_m2']['absorbed']
    # the figures the climate README and the register give: 0.61 K below, about 0.59 in place of 0.71
    assert round(t['design'] - t['held'], 2) == 0.61 and round(g['efficacy'], 2) == 0.59
    assert round(g['in_points_of_dimming'], 2) == 0.55
    assert round(g['with_the_gases_at_their_forcing']['efficacy'], 2) == 0.65
    # at 0.59 the glow's term falls from 4.24 to 3.59 W/m2, 0.29 points off the shared dimming
    glow = product['heat_beyond_sunlight']['glow']
    assert [round(x, 2) for x in glow['absorbed_sunlight_w_m2']] == [3.59, 4.24]
    absorbed = product['runs']['test']['top_of_air_w_m2']['absorbed']
    drop = 100.0 * (glow['absorbed_sunlight_w_m2'][1] - glow['absorbed_sunlight_w_m2'][0]) / absorbed
    assert round(drop, 2) == 0.29
    shared = product['shared_dimming']['dimming_percent']
    assert shared[1] - shared[0] > drop


@pytest.mark.skipif(not extract_present(), reason='the GCM extract (research/runs, ignored) is not here')
def test_the_glow_pair_is_bound_to_its_runs_records(product):
    ex = json.loads((R.ROOT / R.EXTRACT).read_text())
    for run in (R.DESIGN, R.WW, R.HELD):
        assert f'climate/gcm/runs/{run}/progress.json' in ex['producer']['read'], run
    pair, g = ex['glow_pair'], product['glow_efficacy']
    assert (pair[R.DESIGN]['ground_k'], pair[R.WW]['ground_k'], pair[R.HELD]['ground_k']) == (
        g['ground_k']['design'], g['ground_k']['watt_for_watt'], g['ground_k']['held'])
    for run in (R.WW, R.HELD):
        assert pair[run]['co2_ppm'] == 755.0 and pair[run]['top_infrared_w_m2'] == dict(day=10.1718, night=0.5165)
    assert pair[R.DESIGN]['sunlight_scale'] == R.DESIGN_SCALE


def test_the_dimming_curve_runs_from_the_lock_to_where_the_warm_moon_settles(product):
    c = product['dimming_curve']
    absorbed = product['runs']['test']['top_of_air_w_m2']['absorbed']
    fits = product['settled']['fits']
    assert c['per_kelvin_percent'] == pytest.approx(sorted(100.0 * f['feedback_w_m2_k'] / absorbed for f in fits))
    assert [round(x, 2) for x in c['per_kelvin_percent']] == [0.73, 0.96]
    rows = {r['ground_k']: r for r in c['targets']}
    assert min(rows) == 295.0 and max(rows) == 307.0 and len(rows) == 25
    t = {x['name']: x for x in product['targets']}
    for k, name in ((295.0, '295 K'), (300.0, '300 K')):
        for field in ('dimming_percent', 'with_the_heat_beyond_sunlight_percent', 'sunlight_scale',
                      'sunlight_scale_of_its_design_case'):
            assert rows[k][field] == pytest.approx(t[name][field], rel=1e-12), (k, field)
    temps = sorted(rows)
    for lo, hi in zip(temps, temps[1:]):
        for e in (0, 1):
            assert rows[hi]['dimming_percent'][e] < rows[lo]['dimming_percent'][e]
            assert rows[hi]['sunlight_scale'][e] > rows[lo]['sunlight_scale'][e]
    # 307 K lies between the two fits' balances: one fit would dim a little, the other brighten
    assert rows[307.0]['dimming_percent'][0] < 0.0 < rows[307.0]['dimming_percent'][1]
    shared = product['shared_dimming']['dimming_percent']
    for r in c['targets']:
        assert r['with_the_heat_beyond_sunlight_percent'] == pytest.approx(
            [r['dimming_percent'][e] + shared[e] for e in (0, 1)], rel=1e-12)


def test_half_sunlit_marks_a_mean_the_sun_is_up_for_half_of():
    w = 2.0 * np.pi / 10.0                                     # a tenth of the lunar day

    def mean(start, end):
        """The mean of the Sun's clipped cosine over hour angles start to end, at the equator."""
        x = np.linspace(start, end, 200001)
        return np.trapezoid(np.maximum(np.cos(x), 0.0), x) / (end - start)

    assert X.half_sunlit(w) == pytest.approx(mean(np.pi / 2.0 - w / 2.0, np.pi / 2.0 + w / 2.0), rel=1e-6)
    assert mean(np.pi / 2.0 - 0.6 * w, np.pi / 2.0 + 0.4 * w) > X.half_sunlit(w)       # the Sun up for 60% of it
    assert mean(np.pi / 2.0 - 0.4 * w, np.pi / 2.0 + 0.6 * w) < X.half_sunlit(w)       # and for 40%
    assert mean(np.pi / 2.0 - w, np.pi / 2.0) > X.half_sunlit(w)                       # wholly in the day
    assert round(X.half_sunlit(w), 3) == 0.078


def test_heat_in_weighs_the_selected_land_by_area():
    tw = np.array([[[20.0, 29.0]], [[22.0, 31.5]]])            # two means of two cells: wet-bulb
    air = tw + 5.0
    land = np.ones(tw.shape, bool)
    wl = np.array([[1.0, 3.0]])
    sel = np.array([[[True, True]], [[False, True]]])
    out = X.heat_in(dict(at_2_m=(tw, air)), wl, land, sel)
    assert out['share_of_land_time'] == pytest.approx((1.0 + 3.0 + 3.0) / 8.0)
    assert out['at_2_m']['above_28'] == pytest.approx(6.0 / 7.0)
    assert out['at_2_m']['above_31'] == pytest.approx(3.0 / 7.0)
    assert out['at_2_m']['air']['top'] == 36.5


def test_the_day_and_the_night_split_the_lands_time(product):
    for run in ('design', 'test'):
        land = product['runs'][run]['heat']['land']
        day, night = land['day'], land['night']
        assert day['share_of_land_time'] + night['share_of_land_time'] == pytest.approx(1.0)
        assert abs(day['share_of_land_time'] - 0.5) < 0.01
        for height in ('at_2_m', 'lowest_level'):
            for var in ('wet_bulb', 'air'):
                for q in ('median', 'p99'):
                    assert day[height][var][q] > night[height][var][q], (run, height, var, q)
                whole = (land['at_2_m'] if height == 'at_2_m' else land)[var]['median']
                assert night[height][var]['median'] < whole < day[height][var]['median'], (run, height, var)
        assert day['at_2_m']['above_28'] >= night['at_2_m']['above_28']
    # the targets carry them, the resort's interpolated between the runs
    t = {x['name']: x for x in product['targets']}
    a, b = product['runs']['design'], product['runs']['test']
    for k in R.DAY_NIGHT:
        lo, hi = sorted((a[k], b[k]))
        assert lo <= t['300 K']['climate'][k] <= hi, k
        assert t['no further dimming']['climate'][k] == b[k]
        assert t['295 K']['climate'][k] == pytest.approx(a[k], abs=0.05), k      # the design case is at 295.02 K


def test_the_short_weather_stretches_the_tails_a_little(product):
    s = product['smoothing']
    assert s['run'] == 'A28_corrected_970_eddies' and s['years'] == [0, 1] and s['hours_each'] == 6.5
    assert s['means_to_a_3_day_mean'] == 11
    design = product['runs']['design']['heat']['land']
    for height, mine in (('at_2_m', design['at_2_m']), ('lowest_level', design)):
        for var in ('wet_bulb', 'air'):
            x = s[height][var]
            assert abs(x['median']['at_6_5_hours'] - x['median']['in_3_day_means']) < 0.05
            for q in ('p99', 'top'):
                assert 0.0 < x[q]['excess_k'] < 1.5 and x[q]['ratio_above_the_median'] > 1.0, (height, var, q)
                assert x[q]['excess_k'] == pytest.approx(x[q]['at_6_5_hours'] - x[q]['in_3_day_means'])
            # the 3-day means the two years make match the design case's own 3-day output over its years 10-19
            assert abs(x['median']['in_3_day_means'] - mine[var]['median']) < 0.1
            assert abs(x['p99']['in_3_day_means'] - mine[var]['p99']) < 0.1
    for t in product['targets']:
        c = t['climate']
        for k, (m, p, height, var) in R.SMOOTHED.items():
            r = s[height][var]['p99']['ratio_above_the_median']
            assert t['at_6_5_hours'][f'{k}_p99_c'] == pytest.approx(c[m] + r * (c[p] - c[m]), rel=1e-12)
            assert t['at_6_5_hours'][f'{k}_p99_c'] > c[p]


def test_the_product_agrees_with_the_study_readme(product):
    """The figures the study README gives for the additions of 17 October on domain/climate."""
    t = {x['name']: x for x in product['targets']}
    warm, resort, lock = t['no further dimming'], t['300 K'], t['295 K']
    one = lambda xs: [round(x, 1) for x in xs]
    assert one(product['shared_dimming']['dimming_percent']) == [3.6, 5.3]
    assert one(product['heat_beyond_sunlight']['answered_in_every_target']['absorbed_sunlight_w_m2']) == [8.0, 12.0]
    assert [round(x, 2) for x in product['shared_dimming']['growth_cost_through_light_percent']] == [0.25, 1.56]
    assert one(product['shared_dimming']['land_rain_cut_percent']) == [2.0, 7.3]
    for x, total, scale, answered in ((warm, [3.6, 5.3], [0.919, 0.935], [1.0, 1.8]),
                                      (resort, [8.8, 11.8], [0.855, 0.885], [13.3, 17.3]),
                                      (lock, [12.5, 16.6], [0.809, 0.849], [22.0, 28.6])):
        assert one(x['with_the_heat_beyond_sunlight_percent']) == total, x['name']
        assert [round(v, 3) for v in x['sunlight_scale']] == scale, x['name']
        assert one(x['gap_with_the_added_heat_answered_w_m2']) == answered, x['name']
    rows = {r['ground_k']: r for r in product['dimming_curve']['targets']}
    table = {296.0: ([8.2, 10.4], [0.870, 0.891], [11.7, 15.7], [0.818, 0.856]),
             302.0: ([3.8, 4.6], [0.925, 0.934], [7.3, 9.9], [0.874, 0.899]),
             306.0: ([0.8, 0.8], [0.962, 0.962], [4.4, 6.1], [0.911, 0.928])}
    for k, (dim, own, total, scale) in table.items():
        r = rows[k]
        assert one(r['dimming_percent']) == dim and one(r['with_the_heat_beyond_sunlight_percent']) == total, k
        assert [round(v, 3) for v in r['sunlight_scale_of_its_design_case']] == own, k
        assert [round(v, 3) for v in r['sunlight_scale']] == scale, k
    # day and night: the wet-bulb at 2 m (median, 99th percentile) and the air at 2 m (median), day then night
    for x, figures in ((warm, (27.1, 26.8, 32.1, 31.6, 34.6, 32.9)), (resort, (21.5, 20.8, 26.5, 25.6, 27.6, 25.3)),
                       (lock, (17.3, 16.2, 22.3, 21.0, 22.2, 19.5))):
        c = x['climate']
        got = (c['day_wet_bulb_2m_median_c'], c['night_wet_bulb_2m_median_c'], c['day_wet_bulb_2m_p99_c'],
               c['night_wet_bulb_2m_p99_c'], c['day_land_air_2m_median_c'], c['night_land_air_2m_median_c'])
        assert tuple(round(v, 1) for v in got) == figures, x['name']
    land = product['runs']['test']['heat']['land']
    assert round(100 * land['day']['at_2_m']['above_28']) == 35
    assert round(100 * land['night']['at_2_m']['above_28']) == 29
    # the weather shorter than three days, and the 99th percentile of the wet-bulb at 2 m in 6.5-hour means
    s = product['smoothing']['at_2_m']['wet_bulb']
    assert (round(s['p99']['excess_k'], 2), round(s['top']['excess_k'], 2)) == (0.32, 0.94)
    assert [round(x['at_6_5_hours']['wet_bulb_2m_p99_c'], 1) for x in (warm, resort, lock)] == [32.2, 26.5, 22.2]
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


def test_the_glow_pair_agrees_with_the_climate_readme(product):
    """The figures the GCM README gives for the glow pair on the corrected radiation (The gases, the glow and the
    dimming), which the product binds."""
    g = product['glow_efficacy']
    d, t = g['added_dimming_percent'], g['ground_k']
    assert (round(t['design'], 2), round(t['watt_for_watt'], 2), round(t['held'], 2)) == (295.02, 293.45, 294.41)
    assert round(d['watt_for_watt'] - d['held'], 2) == 0.70 and round(t['held'] - t['watt_for_watt'], 2) == 0.96
    assert round(g['response_k_per_percent'], 2) == 1.37
    assert round(d['held'] - g['dimming_that_holds_percent'], 2) == 0.45
    assert round(g['dimming_that_holds_percent'], 2) == 3.07
    assert (round(g['gases_dimming_percent'], 2), round(g['glow_dimming_percent'], 2)) == (1.69, 1.38)
    assert (round(g['glow_w_m2'], 2), round(g['infrared_w_m2'], 2), round(g['efficacy'], 2)) == (3.15, 5.34, 0.59)
    assert round(g['design_absorbed_sunlight_w_m2'], 1) == 228.5
