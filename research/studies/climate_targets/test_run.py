"""What each target climate takes: the product's binding and regeneration, the arithmetic behind it, and the figures
the study README and the GCM README state, read from the READMEs themselves."""
import hashlib
import json
import math
import re

import numpy as np
import pytest

from research.studies.climate_targets import extract as X
from research.studies.climate_targets import run as R
from research.studies.domain_synthesis import ledgers as L

STUDY_README = R.HERE / 'README.md'
GCM_README = R.ROOT / 'climate' / 'gcm' / 'README.md'
NUMBER = re.compile(r'[−+-]?\d+(?:,\d{3})*(?:\.\d+)?')
SMALL = 50_000_000                   # files the provenance test hashes again; the output years are checked by a re-run


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


def extract_present():
    return (R.ROOT / R.EXTRACT).is_file()


# ---------------------------------------------------------------------------------------------- reading the READMEs
def section(path, heading):
    """A README section's text, from its heading to the next heading of the same or a higher level."""
    lines = path.read_text().splitlines()
    level = len(heading) - len(heading.lstrip('#'))
    start = lines.index(heading)
    end = next((i for i in range(start + 1, len(lines))
                if lines[i].startswith('#') and len(lines[i]) - len(lines[i].lstrip('#')) <= level), len(lines))
    return '\n'.join(lines[start + 1:end])


def flat(text):
    return re.sub(r'\s+', ' ', text)


def table(text, header):
    """The rows of the markdown table whose header row starts with a cell, keyed by their first cell."""
    rows, found = {}, False
    for line in text.splitlines():
        if not line.startswith('|'):
            if found:
                break
            continue
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        if not found:
            found = cells[0] == header
        elif not set(cells[0]) <= set('-'):
            rows[cells[0]] = cells[1:]
    assert found and rows, header
    return rows


def numbers(text):
    """The numbers a cell or a sentence states, in order, as written."""
    return [m.replace('−', '-').replace(',', '').lstrip('+') for m in NUMBER.findall(text)]


def agree(written, values, what=''):
    """Each number as written against a product value, at the decimals the README gives it."""
    assert len(written) == len(values), (what, written, values)
    for w, v in zip(written, values):
        decimals = len(w.split('.')[1]) if '.' in w else 0
        assert f'{v:.{decimals}f}' in (w, w.replace('-', '−')) or f'{v:.{decimals}f}' == f'{float(w):.{decimals}f}', \
            (what, w, v)


def found(pattern, text):
    m = re.search(pattern, flat(text))
    assert m, pattern
    return [g.replace('−', '-') for g in m.groups()]


# ------------------------------------------------------------------------------------------ binding and regeneration
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
    """The extract records the code that made it and the files it read, with their sizes; the small ones are hashed
    again here, the output years by a fresh extract."""
    ex = json.loads((R.ROOT / R.EXTRACT).read_text())
    assert ex['schema'] == X.SCHEMA
    for path, h in ex['producer']['files'].items():
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == h, path
    assert set(ex['producer']['sizes']) == set(ex['producer']['read'])
    small = [p for p, size in ex['producer']['sizes'].items() if size < SMALL]
    missing = [p for p in small if not (R.ROOT / p).is_file()]
    if missing:
        pytest.skip(f'recorded files are not here: {", ".join(missing)}')
    for path in small:
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == ex['producer']['read'][path], path


# ------------------------------------------------------------------------------------------------ the extract's parts
def test_the_extract_weights_quantiles_by_area():
    rng = np.random.default_rng(0)
    x = rng.normal(size=20001)
    assert X.wquantile(x, np.ones_like(x), 0.5) == pytest.approx(np.median(x), abs=1e-12)
    assert X.wquantile(x, np.ones_like(x), 0.99) == pytest.approx(np.percentile(x, 99), abs=1e-3)
    # weights scale out, and a cell twice another's area counts as two of it at the median
    w = rng.uniform(0.5, 2.0, size=x.size)
    assert X.wquantile(x, 3.0 * w, 0.9) == pytest.approx(X.wquantile(x, w, 0.9), rel=1e-12)
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
    z = np.array([30e3, 20e3, 10e3, 6e3, 3e3, 1e3])                                     # levels from the top down
    columns = np.array([[200.0, 230.0, 260.0, 270.0, 280.0, 290.0],                     # an ordinary column
                        [200.0, 230.0, 280.0, 270.0, 278.0, 290.0],                     # a warm layer aloft
                        [200.0, 230.0, 260.0, 270.0, 280.0, 271.0],                     # its lowest level frozen
                        [280.0, 285.0, 290.0, 295.0, 300.0, 305.0]])                    # never freezing
    h = X.freezing_height(columns.T.reshape(1, 6, 1, 4), np.repeat(z[:, None], 4, 1).reshape(1, 6, 1, 4))[0, 0]
    assert h[0] == pytest.approx(3e3 + (280.0 - 273.15) / 10.0 * 3e3)
    assert h[1] == pytest.approx(3e3 + (278.0 - 273.15) / 8.0 * 3e3)        # the first crossing from the ground
    assert h[2] == pytest.approx(1e3)
    assert np.isnan(h[3])


def test_half_sunlit_marks_a_mean_the_sun_is_up_for_half_of():
    w = 2.0 * np.pi / 10.0

    def mean(start, end):
        """The mean of the Sun's clipped cosine over hour angles start to end, at the equator."""
        x = np.linspace(start, end, 200001)
        return np.trapezoid(np.maximum(np.cos(x), 0.0), x) / (end - start)

    assert X.half_sunlit(w) == pytest.approx(mean(np.pi / 2.0 - w / 2.0, np.pi / 2.0 + w / 2.0), rel=1e-6)
    assert mean(np.pi / 2.0 - 0.6 * w, np.pi / 2.0 + 0.4 * w) > X.half_sunlit(w)       # the Sun up for 60% of it
    assert mean(np.pi / 2.0 - 0.4 * w, np.pi / 2.0 + 0.6 * w) < X.half_sunlit(w)       # and for 40%
    assert mean(np.pi / 2.0 - w, np.pi / 2.0) > X.half_sunlit(w)                       # wholly in the day


def test_a_mean_takes_the_sun_through_its_synodic_share_of_the_day():
    """143 half-hour steps at the Sun's synodic rate, from the run's sidereal rotation and year."""
    config = dict(model=dict(steps_per_lunar_day=1430, writes_per_lunar_day=10, timestep_min=30.0),
                  planet=dict(rotationperiod=27.321661, year=365.25636))
    window = X.sun_window(config)
    assert np.degrees(window) == pytest.approx(36.32, abs=0.01)
    assert round(X.half_sunlit(window), 4) == 0.0786


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


def test_ranked_tenths_keep_the_same_means_together():
    x = np.arange(20.0).reshape(20, 1)                         # one cell, twenty means
    y = -x
    t = X.ranked_tenths(x, dict(x=x, y=y), 2)
    assert t['warmest']['x'][0] == pytest.approx(18.5) and t['warmest']['y'][0] == pytest.approx(-18.5)
    assert t['coolest']['x'][0] == pytest.approx(0.5) and t['coolest']['y'][0] == pytest.approx(-0.5)


def test_the_wet_bulb_takes_the_air_where_the_air_is_supersaturated():
    t, p = np.array([290.0]), np.array([1.2e5])
    assert X.wet_bulb(t, 2.0 * X.qsat(t, p), p)[0] == pytest.approx(290.0, abs=1e-6)
    assert X.wet_bulb(t, X.qsat(t, p), p)[0] == pytest.approx(290.0, abs=1e-6)


# ------------------------------------------------------------------------------------------------ the product's parts
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


def test_the_land_is_the_models_own(product):
    for run in ('design', 'test'):
        lw = product['runs'][run]['land_and_water']
        assert lw['land_cells'] == 1289 and lw['water_cells_with_land_fraction_above_half'] == 166
        assert round(100 * lw['water_share_of_area'], 2) == 40.05
        assert product['runs'][run]['land_columns_that_never_freeze'] == 0.0
        land = product['runs'][run]['land']
        bands = [b['desert'] for b in land['bands'].values()]
        assert min(bands) >= 0.0 and max(bands) <= 1.0 and 0.0 < land['all']['desert'] < 1.0


def test_the_2_m_air_is_the_lowest_level_carried_down_the_dry_adiabat(product):
    for run in ('design', 'test'):
        land = product['runs'][run]['heat']['land']
        a = land['at_2_m']['air_down_the_dry_adiabat']
        assert abs(a['mean_k']) < 0.05 and a['sd_k'] < 0.1 and 1.0 < a['above_the_lowest_level_k'] < 2.0
        assert 0.0 <= land['saturated']['capped_at_2_m'] < 0.01
        p = land['surface_pressure_bar']
        assert 1.0 < p['p05'] < p['p50'] < p['p95'] < 1.25


def test_the_tenths_pair_their_means(product):
    for run in ('design', 'test'):
        for height in ('at_2_m', None):
            heat = product['runs'][run]['heat']['land']
            t = (heat['at_2_m'] if height else heat)['tenths']
            wb, air = t['by_wet_bulb'], t['by_air']
            # each ranking holds the most of its own quantity, at both ends
            assert wb['warmest']['wet_bulb'] >= air['warmest']['wet_bulb']
            assert air['warmest']['air'] >= wb['warmest']['air']
            assert wb['coolest']['wet_bulb'] <= air['coolest']['wet_bulb']
            assert air['coolest']['air'] <= wb['coolest']['air']
            for s in (wb['warmest'], wb['coolest'], air['warmest'], air['coolest']):
                assert 0.0 <= s['afternoon'] <= s['day'] <= 1.0
            # the air's warmest tenth is mostly the day's afternoons, its coolest mostly nights
            assert air['warmest']['day'] > 0.9 and air['warmest']['afternoon'] > 0.6 and air['coolest']['day'] < 0.25


def test_the_heat_beyond_sunlight_adds_up(product):
    h = product['heat_beyond_sunlight']
    absorbed = product['runs']['test']['top_of_air_w_m2']['absorbed']
    for e in (0, 1):
        total = sum(h[k]['absorbed_sunlight_w_m2'][e] for k in ('gases', 'glow', 'beamed_power'))
        assert math.isclose(h['answered_in_every_target']['absorbed_sunlight_w_m2'][e], total, rel_tol=1e-12)
        assert math.isclose(h['answered_in_every_target']['dimming_percent'][e], 100.0 * total / absorbed,
                            rel_tol=1e-12)
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
        assert x['absorbed_sunlight_to_shed_w_m2'] == sorted(f['feedback_w_m2_k'] * (f['balance_k'] - target)
                                                             for f in fits)
    # the hot world's clear sky with the added heat unanswered would have to send out more than a saturated column can
    assert warm['gap_with_the_added_heat_unanswered_w_m2'][1] < 0.0 < warm['gap_with_the_added_heat_answered_w_m2'][0]
    # the interpolated climate lies between the two runs, threshold shares stay with the runs, and the lock, a little
    # below the design case's 295.02 K, takes the design case's own climate
    a, b = product['runs']['design'], product['runs']['test']
    for k in ('land_desert', 'wet_bulb_median_c', 'sea_median_c', 'land_rain_mm_day'):
        lo, hi = sorted((a[k], b[k]))
        assert lo <= mid['climate'][k] <= hi, k
    assert not any(k in mid['climate'] for k in R.RUNS_ONLY)
    assert cool['interpolated_share_of_the_way'] == 0.0
    for k, v in cool['climate'].items():
        if k != 'ground_k':
            assert v == a[k], k


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
    """The pair fixes the total the gases and the glow take; with the gases at their forcing the glow takes the rest,
    and the other split gives the same total."""
    g = product['glow_efficacy']
    d, t = g['added_dimming_percent'], g['ground_k']
    assert round(d[R.WW], 2) == 4.22 and round(d[R.HELD], 2) == 3.52     # 13 October's added dimming on 0.970
    response = (t[R.HELD] - t[R.WW]) / (d[R.WW] - d[R.HELD])
    assert g['response_k_per_percent'] == pytest.approx(response, rel=1e-12)
    hold = d[R.HELD] - (t[R.DESIGN] - t[R.HELD]) / response
    assert g['dimming_that_holds_percent'] == pytest.approx(hold, rel=1e-12)
    absorbed = product['runs']['design']['top_of_air_w_m2']['absorbed']
    assert g['design_absorbed_sunlight_w_m2'] == absorbed
    total = hold / 100.0 * absorbed
    assert g['gases_and_glow_w_m2'] == pytest.approx(total, rel=1e-12)
    gas = product['heat_beyond_sunlight']['gases']['forcing_w_m2']
    assert g['gases_w_m2'] == gas and g['efficacy'] == pytest.approx((total - gas) / g['infrared_w_m2'], rel=1e-12)
    other = g['with_the_gases_at_their_points']
    points = json.loads((R.ROOT / R.GASES_GLOW).read_text())['runs']['gases']['added_dimming_percent']
    assert other['gases_w_m2'] == pytest.approx(points / 100.0 * absorbed, rel=1e-12)
    for split in (g, other):
        assert split['gases_w_m2'] + split['glow_w_m2'] == pytest.approx(total, rel=1e-12)
    assert round(g['efficacy'], 3) == 0.645 and round(other['efficacy'], 2) == 0.59
    # at 0.645 the glow's term falls from 4.24 to 3.89 W/m2, 0.16 points off the shared dimming
    glow = product['heat_beyond_sunlight']['glow']
    assert [round(x, 2) for x in glow['absorbed_sunlight_w_m2']] == [3.89, 4.24]
    drop = 100.0 * (glow['absorbed_sunlight_w_m2'][1] - glow['absorbed_sunlight_w_m2'][0]) / \
        product['runs']['test']['top_of_air_w_m2']['absorbed']
    assert round(drop, 2) == 0.16


@pytest.mark.skipif(not extract_present(), reason='the GCM extract (research/runs, ignored) is not here')
def test_the_glow_pair_is_bound_to_its_runs_records(product):
    ex = json.loads((R.ROOT / R.EXTRACT).read_text())
    for run in (R.DESIGN, R.WW, R.HELD):
        assert f'climate/gcm/runs/{run}/progress.json' in ex['producer']['read'], run
    pair, g = ex['glow_pair'], product['glow_efficacy']
    for run in (R.DESIGN, R.WW, R.HELD):
        assert pair[run]['ground_k'] == g['ground_k'][run]
    for run in (R.WW, R.HELD):
        assert pair[run]['co2_ppm'] == 755.0 and pair[run]['top_infrared_w_m2'] == dict(day=10.1718, night=0.5165)
    assert pair[R.DESIGN]['sunlight_scale'] == R.DESIGN_SCALE
    assert ex['producer']['read'][X.START].startswith('cc83ad0f')


def test_the_dimming_curve_runs_from_the_lock_to_where_the_hot_world_settles(product):
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


def test_the_day_and_the_night_split_the_lands_time(product):
    for run in ('design', 'test'):
        land = product['runs'][run]['heat']['land']
        day, night = land['day'], land['night']
        assert day['share_of_land_time'] + night['share_of_land_time'] == pytest.approx(1.0)
        assert abs(day['share_of_land_time'] - 0.5) < 0.01
        assert round(day['czen_above']['times_cos_latitude'], 4) == 0.0786
        assert round(day['czen_above']['window_degrees'], 2) == 36.32
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
        assert t['no further dimming']['climate'][k] == b[k] and t['295 K']['climate'][k] == a[k]


def test_the_wet_bulb_per_kelvin_of_air_is_bound(product):
    for run in ('design', 'test'):
        s = product['wet_bulb_per_kelvin_of_air'][run]
        land = product['runs'][run]['heat']['land']
        assert s['air_c'] == pytest.approx(land['day']['at_2_m']['air']['median'], abs=1e-9)
        assert s['wet_bulb_c'] == pytest.approx(land['day']['at_2_m']['wet_bulb']['median'], abs=1e-9)
        for at in s['at']:
            p, q, t = np.array([at['pressure_bar'] * 1e5]), np.array([at['specific_humidity']]), s['air_c'] + 273.15
            assert X.wet_bulb(np.array([t]), q, p)[0] - 273.15 == pytest.approx(s['wet_bulb_c'], abs=1e-6)
            change = X.wet_bulb(np.array([t + 0.5]), q, p)[0] - X.wet_bulb(np.array([t - 0.5]), q, p)[0]
            assert change == pytest.approx(at['k_per_k'], abs=1e-9)
    assert [round(x, 2) for x in product['wet_bulb_per_kelvin_of_air']['test']['k_per_k']] == [0.26, 0.28]
    assert [round(x, 2) for x in product['wet_bulb_per_kelvin_of_air']['design']['k_per_k']] == [0.37, 0.39]


def test_the_short_weather_stretches_the_tails_a_little(product):
    s = product['smoothing']
    assert s['run'] == 'A28_corrected_970_eddies' and s['years'] == [0, 1] and s['hours_each'] == 6.5
    assert s['means_to_a_3_day_mean'] == 11
    assert s['start']['path'] == X.START and s['start']['sha256'].startswith('cc83ad0f')
    assert 0.0 < s['supersaturated_lowest_level']['in_3_day_means'] < s['supersaturated_lowest_level']['at_6_5_hours']
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


# ---------------------------------------------------------------------------------------- the study README's figures
def targets(product):
    t = {x['name']: x for x in product['targets']}
    return t['no further dimming'], t['300 K'], t['295 K']


def test_the_study_readme_states_the_heat_beyond_the_sunlight(product):
    text = section(STUDY_README, '## The heat beyond the sunlight')
    rows = table(text, 'What the plan adds')
    h = product['heat_beyond_sunlight']
    glow = rows["The fleet's glow, both halves rolled"]
    agree(numbers(glow[0]), h['glow']['absorbed_sunlight_w_m2'], 'glow W/m2')
    agree(numbers(glow[1]), h['glow']['dimming_percent'], 'glow %')
    agree(numbers(glow[2])[1:3], h['glow']['efficacy'], 'efficacy')
    agree(numbers(rows['**Together**'][0]), h['answered_in_every_target']['absorbed_sunlight_w_m2'], 'together')
    agree(numbers(rows['**Together**'][1]), h['answered_in_every_target']['dimming_percent'], 'together %')
    g = product['glow_efficacy']
    lo, hi = found(r'The glow\'s efficacy runs from ([\d.]+) to ([\d.]+)\.', text)
    agree([lo, hi], h['glow']['efficacy'])
    agree(found(r'take ([\d.]+)% of its sunlight, ([\d.]+) W/m² of its ([\d.]+) W/m² absorbed', text),
          [g['dimming_that_holds_percent'], g['gases_and_glow_w_m2'], g['design_absorbed_sunlight_w_m2']])
    agree(found(r'the glow takes the rest, ([\d.]+) W/m², ([\d.]+) of its ([\d.]+) W/m² of infrared', text),
          [g['glow_w_m2'], g['efficacy'], g['infrared_w_m2']])
    other = g['with_the_gases_at_their_points']
    agree(found(r'larger absorbed sunlight, ([\d.]+) W/m², and leaves the glow ([\d.]+),', text),
          [other['gases_w_m2'], other['efficacy']])
    drop = 100.0 * (h['glow']['absorbed_sunlight_w_m2'][1] - h['glow']['absorbed_sunlight_w_m2'][0]) / \
        product['runs']['test']['top_of_air_w_m2']['absorbed']
    agree(found(r'falls from ([\d.]+) to ([\d.]+) W/m², which takes ([\d.]+) points off', text),
          h['glow']['absorbed_sunlight_w_m2'][::-1] + [drop])


def test_the_study_readme_states_the_three_targets(product):
    rows = table(section(STUDY_README, '## The three targets'), '')
    warm, resort, lock = targets(product)
    for key, field in (('With the heat beyond the sunlight', 'with_the_heat_beyond_sunlight_percent'),
                       ('Sunlight scale with the added heat answered too', 'sunlight_scale'),
                       ('That gap with the added heat at full watts: answered; unanswered', None)):
        for cell, x in zip(rows[key], (warm, resort, lock)):
            if field:
                agree(numbers(cell), x[field], key)
            else:
                agree(numbers(cell), x['gap_with_the_added_heat_answered_w_m2']
                      + x['gap_with_the_added_heat_unanswered_w_m2'], key)
    rain = rows['Rain on land, a day']
    agree(numbers(rain[0]), [warm['land_rain_mm_day']])
    for cell, x in zip(rain[1:], (resort, lock)):
        agree(numbers(cell), [x['land_rain_mm_day']] + x['land_rain_with_the_dimmings_own_cut_mm_day'], 'rain')
    for cell, x in zip(rows['Land under 0.5 mm a day'], (warm, resort, lock)):
        agree(numbers(cell), [100 * x['climate']['land_desert']], 'desert')
    for key, pair in (('Land wet-bulb at 2 m: median; 99th percentile', ('wet_bulb_2m_median_c', 'wet_bulb_2m_p99_c')),
                      ('Land wet-bulb 0.9 km up: median; 99th percentile', ('wet_bulb_median_c', 'wet_bulb_p99_c')),
                      ('Land air, median: at 2 m; 0.9 km up', ('land_air_2m_median_c', 'land_air_median_c'))):
        for cell, x in zip(rows[key], (warm, resort, lock)):
            agree(numbers(cell)[:2], [x['climate'][k] for k in pair], key)
    for cell, x in zip(rows["Seas' surface, median"], (warm, resort, lock)):
        agree(numbers(cell), [x['climate']['sea_median_c']])
    twenty = rows["Where the land's mean profile reaches 20 °C"]
    agree(numbers(twenty[0]) + numbers(twenty[1]), [warm['climate']['height_of_20c_km'],
                                                    resort['climate']['height_of_20c_km']])
    assert lock['climate']['height_of_20c_km'] is None
    for cell, x in zip(rows["Where the land's mean profile freezes"], (warm, resort, lock)):
        agree(numbers(cell), [x['climate']['freezing_level_km']])


def test_the_study_readme_states_the_dimming_curve(product):
    rows = table(section(STUDY_README, '## The dimming at every temperature'), 'Ground')
    curve = {r['ground_k']: r for r in product['dimming_curve']['targets']}
    for k in range(295, 307):
        r, cells = curve[float(k)], rows[f'{k} K']
        for cell, field in zip(cells, ('dimming_percent', 'sunlight_scale_of_its_design_case',
                                       'with_the_heat_beyond_sunlight_percent', 'sunlight_scale')):
            written, values = numbers(cell), r[field]
            agree(written, values if len(written) == 2 else values[:1], (k, field))
    warm = targets(product)[0]
    agree(numbers(rows['306.8–307.1 K, where it settles'][2]), warm['with_the_heat_beyond_sunlight_percent'])
    agree(numbers(rows['306.8–307.1 K, where it settles'][3]), warm['sunlight_scale'])


def test_the_study_readme_states_the_day_and_the_night(product):
    text = section(STUDY_README, '## Day and night')
    rows = table(text, "Land: the day's means; the night's")
    keys = {'Wet-bulb at 2 m: median': 'wet_bulb_2m_median_c', 'Wet-bulb at 2 m: 99th percentile': 'wet_bulb_2m_p99_c',
            'Air at 2 m: median': 'land_air_2m_median_c', 'Air at 2 m: 99th percentile': 'land_air_2m_p99_c',
            'Wet-bulb 0.9 km up: median': 'wet_bulb_median_c', 'Wet-bulb 0.9 km up: 99th percentile': 'wet_bulb_p99_c',
            'Air 0.9 km up: median': 'land_air_median_c', 'Air 0.9 km up: 99th percentile': 'land_air_p99_c'}
    for label, k in keys.items():
        for cell, x in zip(rows[label], targets(product)):
            agree(numbers(cell), [x['climate'][f'day_{k}'], x['climate'][f'night_{k}']], label)
    day = product['runs']['test']['heat']['land']['day']
    agree(found(r'through ([\d.]+)° of its hour angle', text), [day['czen_above']['window_degrees']])
    agree(found(r'has a czen of ([\d.]+) times the cosine', text), [day['czen_above']['times_cos_latitude']])
    s = product['wet_bulb_per_kelvin_of_air']
    agree(found(r'moves the wet-bulb temperature ([\d.]+)–([\d.]+) K at the hot world\'s day medians and ([\d.]+)–'
                r'([\d.]+) K at the lock\'s', text), s['test']['k_per_k'] + s['design']['k_per_k'])
    land = product['runs']['test']['heat']['land']
    agree(found(r'passes 28 °C for ([\d.]+)% of the day\'s land time and ([\d.]+)% of the night\'s, and 31 °C for '
                r'([\d.]+)% and ([\d.]+)%; at 0\.9 km, ([\d.]+)% and ([\d.]+)%, and ([\d.]+)% and ([\d.]+)%', text),
          [100 * land[p][h][k] for h in ('at_2_m', 'lowest_level') for k in ('above_28', 'above_31')
           for p in ('day', 'night')])
    for run, figure in (('test', 'stands ([\\d.]+) K above the ground\'s own temperature'), ('design', None)):
        coolest = product['runs'][run]['heat']['land']['at_2_m']['tenths']['by_air']['coolest']
        if figure:
            agree(found(figure, text), [coolest['air'] - coolest['ground']])
        else:
            agree(found(r'temperature, and ([\d.]+) K in the lock', text), [coolest['air'] - coolest['ground']])


def test_the_study_readme_states_the_short_weather(product):
    text = section(STUDY_README, '## The weather shorter than three days')
    s = product['smoothing']
    rows = table(text, 'Design case, land, its two 6.5-hour years')
    for label, height, var in (('Wet-bulb at 2 m', 'at_2_m', 'wet_bulb'), ('Air at 2 m', 'at_2_m', 'air'),
                               ('Wet-bulb 0.9 km up', 'lowest_level', 'wet_bulb'),
                               ('Air 0.9 km up', 'lowest_level', 'air')):
        x = s[height][var]
        agree(numbers(rows[label][0]), [x['median']['at_6_5_hours'], x['median']['in_3_day_means']], label)
        for cell, q in zip(rows[label][1:], ('p99', 'top')):
            agree(numbers(cell), [x[q]['at_6_5_hours'], x[q]['in_3_day_means'], x[q]['excess_k'],
                                  x[q]['ratio_above_the_median']], (label, q))
    agree(found(r'In the 6\.5-hour means ([\d.]+)% of the land\'s lowest-level means are supersaturated, against '
                r'([\d.]+)% of the 3-day means', text),
          [100 * s['supersaturated_lowest_level'][k] for k in ('at_6_5_hours', 'in_3_day_means')])
    rows = table(text, 'Land, 99th percentile')
    for label, k in (('Wet-bulb at 2 m', 'wet_bulb_2m'), ('Air at 2 m', 'land_air_2m'),
                     ('Wet-bulb 0.9 km up', 'wet_bulb'), ('Air 0.9 km up', 'land_air')):
        for cell, x in zip(rows[label], targets(product)):
            agree(numbers(cell), [x['at_6_5_hours'][f'{k}_p99_c'], x['climate'][f'{k}_p99_c']], label)


# ------------------------------------------------------------------------------------------ the GCM README's figures
def settled():
    return section(GCM_README, '### Where the corrected Moon settles')


def test_the_climate_readme_states_the_settled_heat(product):
    """The GCM README's prose figures for the settled test (Where the corrected Moon settles), with the design case's
    in brackets."""
    text = settled()
    test, design = product['runs']['test'], product['runs']['design']
    lt, ld = test['heat']['land'], design['heat']['land']
    agree(found(r'from ([\d.]+)% to ([\d.]+)%\. The rain is still heaviest', text),
          [100 * design['land_desert'], 100 * test['land_desert']])
    poles = [test['land']['bands'][b]['ground_c'] for b in ('60-75', '75-90')]
    tropics = [test['land']['bands'][b]['ground_c'] for b in ('0-15', '15-30')]
    agree(found(r'poleward of 60° is the warmest, ([\d.]+)–([\d.]+) °C against ([\d.]+)–([\d.]+) °C within 30°', text),
          sorted(poles) + sorted(tropics))
    agree(found(r'the land\'s air at 2 m has a median of ([\d.]+) °C \(([\d.]+) in the design case\), ([\d.]+) °C at '
                r'the 99th percentile and ([\d.]+) °C at most', text),
          [lt['at_2_m']['air']['median'], ld['at_2_m']['air']['median'], lt['at_2_m']['air']['p99'],
           lt['at_2_m']['air']['top']])
    agree(found(r'Its wet-bulb temperature at 2 m has a median of ([\d.]+) °C \(([\d.]+)\), ([\d.]+) °C at the 99th '
                r'percentile \(([\d.]+)\) and ([\d.]+) °C at most \(([\d.]+)\), and passes 28 °C for ([\d.]+)% of the '
                r'land\'s time and 31 °C for ([\d.]+)%', text),
          [lt['at_2_m']['wet_bulb']['median'], ld['at_2_m']['wet_bulb']['median'], lt['at_2_m']['wet_bulb']['p99'],
           ld['at_2_m']['wet_bulb']['p99'], lt['at_2_m']['wet_bulb']['top'], ld['at_2_m']['wet_bulb']['top'],
           100 * lt['at_2_m']['above_28'], 100 * lt['at_2_m']['above_31']])
    agree(found(r'the air\'s median is ([\d.]+) °C \(([\d.]+)\) and the wet-bulb temperature\'s ([\d.]+) °C '
                r'\(([\d.]+)\), with ([\d.]+) °C at the 99th percentile \(([\d.]+)\), passing 28 °C for ([\d.]+)% of '
                r'the land\'s time and 31 °C for ([\d.]+)%', text),
          [lt['air']['median'], ld['air']['median'], lt['wet_bulb']['median'], ld['wet_bulb']['median'],
           lt['wet_bulb']['p99'], ld['wet_bulb']['p99'], 100 * lt['above_28'], 100 * lt['above_31']])
    agree(found(r'The seas\' surface has a median of ([\d.]+) °C \(([\d.]+)\) and ([\d.]+) °C at the 99th percentile',
                text), [test['seas']['all']['median'], design['seas']['all']['median'], test['seas']['all']['p99']])
    agree(found(r'columns crosses freezing (\d+) km above the ground \((\d+) km in the design case\)', text),
          [test['freezing_level_over_land_km'], design['freezing_level_over_land_km']])
    lw = test['land_and_water']
    agree(found(r'so (\d+) cells with a land fraction of 0\.50–0\.59 run as water and water covers ([\d.]+)% of '
                r'the Moon', text), [lw['water_cells_with_land_fraction_above_half'], 100 * lw['water_share_of_area']])
    a = lt['at_2_m']['air_down_the_dry_adiabat']
    agree(found(r'it stands ([\d.]+) K above the lowest level\'s air and within ([\d.]+) K of that air times '
                r'σ\^−0\.286 \(standard deviation ([\d.]+) K\)', text),
          [a['above_the_lowest_level_k'], 0.02, a['sd_k']])
    assert abs(a['mean_k']) < 0.02 and abs(ld['at_2_m']['air_down_the_dry_adiabat']['mean_k']) < 0.05
    agree(found(r'binds in ([\d.]+)% of the land\'s means here and ([\d.]+)% in the design case', text),
          [100 * lt['saturated']['capped_at_2_m'], 100 * ld['saturated']['capped_at_2_m']])
    agree(found(r'put the lowest level\'s time above 28 °C at ([\d.]+)%', text),
          [100 * lt['readout_of_17_october']['above_28']])


def test_the_climate_readme_states_the_deserts(product):
    rows = table(settled(), 'Land, years 10–19')
    design, test = product['runs']['design']['land'], product['runs']['test']['land']
    for label, band in (('0–15°', '0-15'), ('15–30°', '15-30'), ('30–45°', '30-45'), ('45–60°', '45-60'),
                        ('60–75°', '60-75'), ('75–90°', '75-90')):
        for cell, land in zip(rows[label], (design, test)):
            agree(numbers(cell), [100 * land['bands'][band]['desert'], land['bands'][band]['rain_mm_day']], label)
    for cell, land in zip(rows['Below 2 km; 2–4 km; 4 km and up'], (design, test)):
        agree(numbers(cell), [100 * land['heights'][k]['desert'] for k in ('below_2_km', '2_to_4_km', 'from_4_km')])
    for cell, land in zip(rows['All land'], (design, test)):
        agree(numbers(cell), [100 * land['all']['desert'], land['all']['rain_mm_day']], 'all land')


def test_the_climate_readme_states_the_places(product):
    """The place-by-place table, the design case's figures in brackets."""
    rows = table(settled(), 'Land')
    test, design = product['runs']['test']['heat'], product['runs']['design']['heat']
    places = {'0–15°': ('bands', '0-15'), '15–30°': ('bands', '15-30'), '30–45°': ('bands', '30-45'),
              '45–60°': ('bands', '45-60'), '60–90°': ('bands', '60-90'), 'Below 1 km': ('below_1_km',),
              '4 km and up': ('from_4_km',), 'Coasts, the land beside the sea': ('coasts',)}
    assert set(rows) == set(places)
    for label, path in places.items():
        t, d = test, design
        for k in path:
            t, d = t[k], d[k]
        wet, time, air = rows[label]
        agree(numbers(wet), [t['at_2_m']['wet_bulb']['median'], t['at_2_m']['wet_bulb']['p99'],
                             d['at_2_m']['wet_bulb']['median'], d['at_2_m']['wet_bulb']['p99']], label)
        assert d['at_2_m']['above_28'] == d['above_28'] == 0.0
        agree(numbers(time), [100 * t['at_2_m']['above_28'], 100 * t['above_28'], 0.0], label)
        agree(numbers(air), [t['at_2_m']['air']['median'], d['at_2_m']['air']['median']], label)


def test_the_climate_readme_states_the_tenths_seas_and_ladder(product):
    text = settled()
    rows = table(text, 'Land, at 2 m')
    test, design = (product['runs'][r]['heat']['land']['at_2_m']['tenths'] for r in ('test', 'design'))
    for label, by, end in (('Warmest tenth by the wet-bulb temperature', 'by_wet_bulb', 'warmest'),
                           ('Coolest tenth by the wet-bulb temperature', 'by_wet_bulb', 'coolest'),
                           ('Warmest tenth by the air', 'by_air', 'warmest'),
                           ('Coolest tenth by the air', 'by_air', 'coolest')):
        t, d = test[by][end], design[by][end]
        wet, air, when = rows[label]
        agree(numbers(wet), [t['wet_bulb'], d['wet_bulb']], label)
        agree(numbers(air), [t['air'], d['air']], label)
        agree(numbers(when), [100 * t['day'], 100 * t['afternoon'], 100 * d['day'], 100 * d['afternoon']], label)
    agree(found(r'the 2 m air stands ([\d.]+) K above the ground\'s own temperature \(([\d.]+) K in the design case\)',
                text), [test['by_air']['coolest']['air'] - test['by_air']['coolest']['ground'],
                        design['by_air']['coolest']['air'] - design['by_air']['coolest']['ground']])
    seas = {r: product['runs'][r]['seas']['bands'] for r in ('test', 'design')}
    agree(found(r'The seas\' surface has a median of ([\d.]+) °C within 15° of the equator and ([\d.]+) °C poleward of '
                r'60° \(([\d.]+) and ([\d.]+)\)', text),
          [seas['test']['0-15']['median'], seas['test']['60-90']['median'], seas['design']['0-15']['median'],
           seas['design']['60-90']['median']])
    lt = {x['height_km']: x for x in product['runs']['test']['ladder']}
    ld = {x['height_km']: x for x in product['runs']['design']['ladder']}
    agree(found(r'the air 10 km up, ([\d.]+) °C at ([\d.]+) atm, is near the design case\'s at the lowest level, '
                r'([\d.]+) °C at ([\d.]+) atm; 20 km up it is ([\d.]+) °C at ([\d.]+) atm and 30 km up ([\d.]+) °C at '
                r'([\d.]+) atm \(([\d.]+) and (−?[\d.]+) °C', text),
          [lt[10]['air_c'], lt[10]['pressure_atm'], ld[0.9]['air_c'], ld[0.9]['pressure_atm'], lt[20]['air_c'],
           lt[20]['pressure_atm'], lt[30]['air_c'], lt[30]['pressure_atm'], ld[20]['air_c'], ld[30]['air_c']])


def test_the_climate_readme_states_the_glow_pair(product):
    """The GCM README's figures for the glow pair on the corrected radiation (The gases, the glow and the dimming)."""
    text = section(GCM_README, '### The gases, the glow and the dimming')
    g = product['glow_efficacy']
    d, t = g['added_dimming_percent'], g['ground_k']
    agree(found(r'The held branch ends ([\d.]+) K below the design case', text), [t[R.DESIGN] - t[R.HELD]])
    agree(found(r'The two branches differ by ([\d.]+) points of dimming and ([\d.]+) K, ([\d.]+) K a point, so holding '
                r'the temperature takes ([\d.]+) points less than the held branch: ([\d.]+)% of the design case\'s '
                r'sunlight, ([\d.]+) W/m² of its ([\d.]+) W/m² absorbed', text),
          [d[R.WW] - d[R.HELD], t[R.HELD] - t[R.WW], g['response_k_per_percent'],
           d[R.HELD] - g['dimming_that_holds_percent'], g['dimming_that_holds_percent'], g['gases_and_glow_w_m2'],
           g['design_absorbed_sunlight_w_m2']])
    agree(found(r'With the gases at their forcing, ([\d.]+) W/m², as on 13 October, the glow takes ([\d.]+) W/m², '
                r'([\d.]+) of its ([\d.]+) W/m² of infrared', text),
          [g['gases_w_m2'], g['glow_w_m2'], g['efficacy'], g['infrared_w_m2']])
    other = g['with_the_gases_at_their_points']
    agree(found(r'Kept at their ([\d.]+)% of this radiation\'s sunlight, ([\d.]+) W/m², the gases would leave the glow '
                r'([\d.]+), with the same total', text),
          [other['gases_dimming_percent'], other['gases_w_m2'], other['efficacy']])
