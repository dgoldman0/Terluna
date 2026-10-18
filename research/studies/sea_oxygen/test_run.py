"""The seas' oxygen: the products' binding, the seas' geometry, the 301 K construction and what the runs must keep."""
import hashlib
import json

import numpy as np
import pytest

from biosphere.ecology import lake as L
from research.studies.lake_oxygen import run as R
from research.studies.sea_oxygen import columns as C
from research.studies.sea_oxygen import extract as X
from research.studies.sea_oxygen import run as RUN
from research.studies.sea_oxygen import seas as K
from research.studies.sea_oxygen import setting as T


@pytest.fixture(scope='module')
def product():
    return json.loads(RUN.OUT.read_text())


@pytest.fixture(scope='module')
def seas():
    return json.loads(K.OUT.read_text())


@pytest.fixture(scope='module')
def columns(product):
    if product['columns'] is None:
        pytest.skip('the sea columns wait for the revised lake model')
    return product


def present(path):
    return (RUN.ROOT / path).is_file()


def digest(path):
    return hashlib.sha256((RUN.ROOT / path).read_bytes()).hexdigest()


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == RUN.SCHEMA
    for path, h in product['producer']['files'].items():
        assert digest(path) == h, path
    for path, h in product['producer']['inputs'].items():
        if path.startswith('research/runs/') and not present(path):
            continue                                     # the GCM's weather and the columns live on the research drive
        assert digest(path) == h, path


def test_seas_product_binds_its_code_and_inputs(seas):
    assert seas['schema'] == K.SCHEMA
    for path, h in seas['producer']['files'].items():
        assert digest(path) == h, path
    for path, h in seas['producer']['inputs'].items():
        if path.startswith('geography/inputs/') and not present(path):
            continue                                     # LOLA, the geoid and the gazetteer are fetched, not committed
        assert digest(path) == h, path


def test_the_seas_are_the_atlas_seas_on_the_fine_grid(seas):
    s = seas['seas']
    assert s['sea_level_m'] == seas['checks']['atlas_sea_level_m'] == seas['checks']['drainage_sea_level_m']
    assert s['area_share_of_moon'] == pytest.approx(seas['checks']['atlas_share'], abs=1e-3)
    for b in seas['bodies']:
        if 'atlas_area_km2' in b:
            assert b['area_km2'] == pytest.approx(b['atlas_area_km2'], rel=0.1), b['body']
    # every square kilometre of sea is stood for by exactly one run basin
    total = sum(sum(v for k, v in b['stands_for'].items() if k.endswith('km2')) for b in seas['run_basins'])
    assert total == pytest.approx(s['area_km2'], rel=1e-4)
    assert sum(b['stands_for']['share_of_the_seas'] for b in seas['run_basins']) == pytest.approx(1.0, abs=1e-3)


def test_each_run_basin_is_described_whole(seas):
    for b in seas['run_basins']:
        a = b['hypsography']['area_fraction']
        assert a[0] == 1.0 and a[-1] == 0.0 and all(x >= y for x, y in zip(a, a[1:])), b['key']
        assert 0.0 < b['mean_depth_m'] < b['deepest_m']
        assert sum(c['area_km2'] for c in b['t21_cells']) == pytest.approx(b['area_km2'], rel=1e-3)
        shelves = list(b['shelf_area_share'].values())
        assert all(x <= y for x, y in zip(shelves, shelves[1:]))
        for other, sill in b['sills_m'].items():
            o = next(x for x in seas['run_basins'] if x['key'] == other)
            assert sill is not None and sill <= min(b['deepest_m'], o['deepest_m'])
            if b['basin_of'].startswith('a region'):
                assert sill < K.PARTITION_M                  # basins apart at the partition depth join above it
    keys = {b['key']: b for b in seas['run_basins']}
    assert max(keys, key=lambda k: keys[k]['area_km2']) in keys
    assert max(b['deepest_m'] for b in seas['run_basins']) == max(x['deepest_m'] for x in seas['bodies'])


def test_the_t21_grid_is_gaussian():
    lat = K.t21_lat()
    assert lat[0] == pytest.approx(85.7606, abs=1e-3) and lat[15] == pytest.approx(2.7689, abs=1e-3)
    assert np.allclose(lat, -lat[::-1])


@pytest.mark.skipif(not present('research/runs/sea_oxygen/forcing.json'), reason='the GCM weather is not here')
def test_the_forcing_comes_from_the_current_extract(seas):
    meta = json.loads((RUN.ROOT / 'research/runs/sea_oxygen/forcing.json').read_text())
    assert meta['schema'] == X.SCHEMA
    for path, h in meta['producer']['files'].items():
        assert digest(path) == h, path
    assert meta['producer']['read']['research/studies/sea_oxygen/results/seas.json'] == digest(K.OUT)
    for s in meta['sites']:
        assert s['lsm'] < 0.5 or s['land_cell']
        assert s['water_surface_m'] < X.SEA_LEVEL_WITHIN_M or s['land_cell']


@pytest.mark.skipif(not present('research/runs/sea_oxygen/forcing.json'), reason='the GCM weather is not here')
def test_a_sea_cell_is_read_at_its_own_level():
    # forcing_for() carries the air from a site's elevation to its lake's level; at a sea cell both are the seas'
    # level, so the forcing is the cell's own weather, as lake_oxygen's windows() reads a water cell unmoved
    st = RUN.setting()
    i = next(i for i, s in enumerate(st.meta['sites']) if not s['land_cell'])
    site = st.site(i)
    moved, unmoved = R.windows(st, 'lock', site, target='resort'), R.windows(st, 'lock', site, at_lake=False)
    for k in ('sw', 'lw', 'air', 'q', 'u10', 'p'):
        assert np.array_equal(moved[k], unmoved[k]), k


def test_the_300k_construction_reproduces_lake_oxygen(product):
    r = product['reproduction']
    assert r['proved']
    assert r['construction']['share_equal'] and r['construction']['ladder_equal'] and r['construction']['lapse_equal']
    assert all(r['waters_identical'].values()) and all(r['sun_phase_identical'].values())
    assert r['saturation_identical']
    assert r['resort_300k']['three_day_means_median_mg_l'] == r['resort_300k']['lake_oxygen'][
        'three_day_means_median_mg_l']


def test_the_301k_share_is_climate_targets_basis(product):
    ct = json.loads(R.TARGETS.read_text())
    a, b = ct['runs']['design']['ground_k'], ct['runs']['test']['ground_k']
    assert product['setting']['ground_k'] == 301.0
    assert product['setting']['share_of_the_way'] == pytest.approx((301.0 - a) / (b - a), rel=1e-12)
    assert 0.5 < product['setting']['share_of_the_way'] < 0.54


def test_saturation_at_301k_lies_between_the_targets(product):
    lo = json.loads(R.OUT.read_text())['saturation']['targets']
    s = product['saturation']['at_301k']
    assert lo['hot']['mg_l'] < s['mg_l'] < lo['resort']['mg_l']
    assert lo['resort']['sea_median_c'] < s['sea_median_c'] < lo['hot']['sea_median_c']
    p = R.DESIGN_PRESSURE_PA / 101325.0
    expected = float(L.oxygen_saturation(s['sea_median_c'], p, product['setting']['salinity_g_kg'],
                                         product['setting']['o2_fraction']))
    assert s['mg_l'] == pytest.approx(expected, abs=2e-3)


def test_every_run_conserves_heat_and_oxygen_and_none_freezes(columns):
    product = columns
    # lake_oxygen's measures and its tests' limits (budget_errors()): the heat error against the column's own heat and
    # as a flux over its surface, rebuilt from the saved decade (run.heat_budget()), and the oxygen error against each
    # oxygen column's content
    for group in ('deep', 'variants', 'pans'):
        for name, case in product[group].items():
            b = case['budget']
            assert abs(b['heat_relative_error']) < 1e-8 and abs(b['heat_error_w_m2']) < 1e-6, (group, name)
            assert b['oxygen_relative_error'] < 1e-6, (group, name)
            assert case['ice_max_m'] == 0.0, (group, name)
    # a deep column takes in next to no net heat over its decade: as a drift of its temperature (its heat is counted
    # from 0 C), less than the settle test allows a pass where it settled, and a few hundredths of a kelvin where it
    # still ages; the saved figure, the error against that net heat (or a millionth of the column's heat), is the
    # rebuilt error's own scale
    for name, case in {**product['deep'], **product['variants']}.items():
        b = case['budget']
        drift_k = abs(b['net_heat_share']) * case['profiles']['deepest_tenth_t_c']['mean']
        settled = case['spin_up']['at_the_analysis_step']['settled']
        assert drift_k < (R.SETTLED_K if settled else 0.05), name
        scale = max(abs(b['net_heat_share']), 1e-6)
        assert b['heat_relative_error'] == pytest.approx(b['saved_heat_relative_error'] * scale, rel=0.02, abs=1e-13)


def test_every_deep_column_runs_the_full_spin_up(columns):
    product = columns
    for name, case in {**product['deep'], **product['variants']}.items():
        spin = case['spin_up']
        expected = C.FROM_PASSES if spin['started_from'] else C.SPIN_PASSES
        assert spin['passes'] == expected and spin['years'] == 10 * expected, name
        assert spin['read_from_pass'] == expected // 2, name
        fine = spin['at_the_analysis_step']
        assert fine['settled'] or fine['passes'] == fine['limit'] == C.FINE_LIMIT, name
        assert 0.0 <= spin['passes_with_a_complete_overturn'] <= 1.0
        for trophic, o in spin['oxygen'].items():
            assert 0.0 <= o['anoxic_share'] <= o['hypoxic_share'] + 1e-9 <= 1.0 + 1e-9, (name, trophic)


def test_the_coldest_column_makes_the_densest_water(columns):
    product = columns
    for key, b in product['bounds'].items():
        assert b['deep_t_c']['coldest'] <= b['deep_t_c']['typical'] + 0.05, key
        # the screen reads its interval from the saved decade's overturns, as the oxygen is read, the coarse beside it
        coldest = product['deep'][f'{key}_coldest']
        n = len(coldest['profiles']['overturns_begin_years'])
        assert b['overturns']['coldest']['saved_decade'] == n > 0, key
        assert n <= coldest['complete_overturns'], key          # an overturn touches one lunar cycle or two
        assert b['screen']['renewal_interval_coldest_years'] == pytest.approx(coldest['years'] / n, abs=0.006), key
        assert b['screen']['coarse']['renewal_interval_coldest_years'] == pytest.approx(
            100.0 / coldest['spin_up']['overturns_a_century'], abs=0.006), key
        for s in (b['screen'], b['screen']['coarse']):
            shares = [m['cold_end_share'] for m in s['margins'].values()]
            assert all(x <= y for x, y in zip(shares, shares[1:])), key       # a wider margin takes in more
            for margin, m in s['margins'].items():
                if m['renewal_interval_with_exchange_years'] is None:
                    continue
                assert m['renewal_interval_with_exchange_years'] >= s['renewal_interval_coldest_years'] - 0.05
                for trophic, o in m['deepest_tenth_mean_with_exchange'].items():
                    assert o is None or o <= b['oxygen'][trophic]['coldest']['mean'] + 1e-9, (key, margin, trophic)
        # the reference demand is the central columns', which the warm kinetics would raise at the deep water's warmth
        assert b['warm_kinetics_factor']['coldest'] > 1.0, key


def test_the_product_labels_its_unsettled_columns(columns):
    product = columns
    for key, a in product['answers'].items():
        for label, row in a.items():
            if not isinstance(row, dict) or 'analysis_step' not in row:
                continue
            step = row['analysis_step']
            assert step['settled'] == row['settled_at_the_analysis_step'], (key, label)
            assert (step['reading'] == 'settled') == step['settled'], (key, label)
            if not step['settled']:
                assert step['oxygen_is'] or step['age_is'], (key, label)
    # Procellarum's column under the other form overturns once at the analysis step, in its last pass, after a run of
    # passes without one: its decade is the column just renewed
    shih = product['answers']['procellarum']['typical_shih_bouffard']['analysis_step']
    assert not shih['settled'] and shih['last_pass_overturned'] and shih['passes_overturning'] == 1
    assert shih['oxygen_is'] == 'an upper value' and shih['age_is'] == 'a lower value'
    assert shih['bottom_age_years']['oldest'] > 30.0
    # every view of the seas names the unsettled columns it takes
    for view, e in product['seas_wide'].items():
        for key, a in product['answers'].items():
            pick = a.get(view) or a['typical']
            assert (key in e['unsettled_columns']) == (not pick['analysis_step']['settled']), (view, key)


def test_the_basins_turn_with_the_moon(columns):
    # the internal Rossby radius c/|f| at middle latitudes, the equatorial radius (c/beta)^(1/2) at the equator
    from research.studies.sea_oxygen import analysis as A
    omega = 2.0 * np.pi / (A.SIDEREAL_MONTH_DAYS * 86400.0)
    assert A.rossby_radius_km(0.1, 45.0) == pytest.approx(0.1 / (2 * omega * np.sin(np.pi / 4)) / 1000.0, rel=1e-12)
    assert A.rossby_radius_km(0.1, 0.0) == pytest.approx(np.sqrt(0.1 * A.MOON_RADIUS / (2 * omega)) / 1000.0,
                                                         rel=1e-12)
    for key, r in columns['rotation'].items():
        lo, mid, hi = r['wave']['c_m_s']
        assert 0.0 < lo <= mid <= hi < 1.0, key
        assert r['wave']['days'] > 0.9 * r['wave']['of_days'], key
    # the large basins span many radii, so the Moon's rotation shapes their answer to the wind
    for key in ('procellarum', 'ingenii', 'nubium', 'humboldtianum', 'crisium'):
        assert columns['rotation'][key]['scale_in_radii'] > 5.0, key


def test_oxygen_falls_with_the_trophic_state(columns):
    product = columns
    for name, case in product['deep'].items():
        means = [case['oxygen'][t]['deepest_tenth']['mean'] for t in RUN.TROPHIC]
        assert all(a >= b - 0.05 for a, b in zip(means, means[1:])), name
    for name, case in product['pans'].items():
        lows = [case['nights'][t]['surface_minimum'][1] for t in RUN.TROPHIC]
        assert all(a >= b - 0.05 for a, b in zip(lows, lows[1:])), name


def test_the_slab_follows_the_gcms_own_water(columns):
    product = columns
    assert product['checks']['slab']
    for check in product['checks']['slab']:
        # the pairing the columns take, at the GCM's lowest level as PlaSim pairs it, follows the GCM's own water
        # closely, and more closely than the 10 m pairing at every cell
        own, ten = check['lowest_level'], check['ten_metres']
        for name in ('lunar_free_convection', 'earth_free_convection'):
            assert own[name]['correlation'] > 0.9, (check['cell'], check['run'], name)
            assert own[name]['correlation'] >= ten[name]['correlation'], (check['cell'], check['run'], name)
        # PlaSim's own pairing and free-convection constant follow its slab ocean closely; carried to lunar gravity
        # the free convection leaves the water warmer, as the review found at the lakes' water cells
        assert abs(own['earth_free_convection']['bias_k']) < 0.6, (check['cell'], check['run'])
        assert own['lunar_free_convection']['model_mean_c'] > own['earth_free_convection']['model_mean_c']


def test_the_physics_is_lake_oxygens(columns):
    product = columns
    for name, case in {**product['deep'], **product['variants']}.items():
        ph = case['physics']
        if ph['deep_mixing'] == 'dissipation':
            assert ph['cap'] == R.CAP['central'], name
        else:
            assert ph['cap'] == R.CAP_WIND, name
        if ph['deep_mixing'] == 'wind':
            assert ph['wind_share'] in R.WIND_SHARE.values(), name
        assert case['spin_up']['passes'] >= C.FROM_PASSES, name
    for key in product['cells']:
        assert product['deep'][f'{key}_typical']['physics']['wind_share'] == R.WIND_SHARE['central']
        for eta in R.ETA_BRACKET:
            v = product['variants'][f'{key}_typical_wind_share_{eta}']
            assert v['physics']['wind_share'] == R.WIND_SHARE[eta]
    # Osborn's rule everywhere but the Lake number's rows and the other form's, which runs the central wind share
    forms = {n: c['diffusivity'] for n, c in {**product['deep'], **product['variants']}.items()}
    assert sorted(n for n, f in forms.items() if f == RUN.FORM) == sorted(f'{k}_typical_{RUN.FORM}'
                                                                          for k in RUN.FORM_BASINS)
    for name, form in forms.items():
        ph = {**product['deep'], **product['variants']}[name]['physics']
        assert (form == 'lake_number') == (ph['deep_mixing'] == 'lake_number'), name
        if form == RUN.FORM:
            assert ph['deep_mixing'] == 'wind' and ph['wind_share'] == R.WIND_SHARE['central'], name
        else:
            assert form in ('osborn', 'lake_number'), name


def test_the_other_form_mixes_the_deep_water_less(columns):
    # Shih et al.'s regimes give 2 nu Re_b^(1/2) where Osborn's rule gives 0.2 nu Re_b, less above Re_b = 100, so the
    # same wind leaves the deep water older, overturning no more often and holding no more oxygen
    assert sorted(columns['forms']) == sorted(RUN.FORM_BASINS)
    for key, f in columns['forms'].items():
        osborn, shih = f['osborn'], f['shih_bouffard']
        assert shih['deep_age_years'] >= osborn['deep_age_years'], key
        assert shih['overturns_a_century'] <= osborn['overturns_a_century'], key
        for t in RUN.TROPHIC:
            assert shih['oxygen'][t]['deepest_tenth_mean'] <= osborn['oxygen'][t]['deepest_tenth_mean'] + 0.05, (key, t)
            assert f['range']['oxygen'][t] == sorted([shih['oxygen'][t]['deepest_tenth_mean'],
                                                      osborn['oxygen'][t]['deepest_tenth_mean']]), (key, t)


def test_the_pans_bracket_lessens_the_floors_ventilation(columns):
    # the fixed dissipation and its cap mix the pans' floors far less than the wind's share spread over the thin water
    # below their mixed layers
    assert columns['pans_bracket']
    for name, b in columns['pans_bracket'].items():
        for t, row in b['nights'].items():
            assert row['fixed_dissipation']['bottom_time_hypoxic_share'] >= row['wind']['bottom_time_hypoxic_share'] \
                - 1e-3, (name, t)


@pytest.mark.skipif(not present('research/runs/sea_oxygen/columns.json'), reason='the columns are not here')
def test_the_index_binds_every_saved_column(columns):
    product = columns
    index = json.loads((RUN.ROOT / 'research/runs/sea_oxygen/columns.json').read_text())
    assert index['schema'] == RUN.INDEX_SCHEMA
    assert product['producer']['inputs']['research/runs/sea_oxygen/columns.json'] == digest(
        'research/runs/sea_oxygen/columns.json')
    for name, entry in index['cases'].items():
        assert digest(f'research/runs/sea_oxygen/cases/{entry["summary"]}') == entry['summary_sha256'], name
        for f, h in entry['files'].items():
            path = f'research/runs/sea_oxygen/cases/{f}'
            assert present(path), path
        if entry['diffusivity'] == RUN.FORM:                 # the other form's code is the code that ran
            assert entry['case']['form_code'] == RUN.FORM_CODE == index['diffusivity_forms']['form_code'], name
    assert len(index['cases']) == product['index']['cases']


def test_the_cells_follow_their_rule(product):
    for key, ch in product['cells'].items():
        t, c = ch['typical'], ch['coldest']
        assert c['area_km2'] >= T.COLD_MIN_KM2 or ch['cells'] == 1, key
        lo, hi = ch['coldest_spell_c_range']                  # the ranges are rounded to 0.001 K
        assert lo - 1e-3 <= c['min'] <= hi + 1e-3, key
        assert ch['mean_c_range'][0] - 1e-3 <= t['mean'] <= ch['mean_c_range'][1] + 1e-3, key
        same = t['index'] == c['index'] or abs(c['min'] - t['min']) < T.SAME_K
        assert ch['one_column'] == same, key
        assert 0.0 < ch['area_share_within_same_k_of_coldest'] <= 1.0


def test_the_product_says_what_it_holds(product):
    if product['columns'] is None:
        assert 'wait' in product['status']
        assert product['cases_planned']
        for name, case in product['cases_planned'].items():
            assert case['basin'] in product['seas']['run_basins'], name
            assert case['cell'] in {product['cells'][case['basin']][c]['index'] for c in ('typical', 'coldest')}
    else:
        assert product['columns']['deep'] == len(product['deep'])


@pytest.mark.skipif(not present('research/runs/sea_oxygen/forcing.json'), reason='the GCM weather is not here')
def test_the_planned_cases_are_the_runners(product):
    if product['columns'] is None:
        planned, _ = RUN.cases(RUN.setting())
        assert set(planned) == set(product['cases_planned'])
