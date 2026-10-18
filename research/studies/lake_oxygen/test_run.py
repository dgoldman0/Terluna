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
            continue                                     # LOLA and the gazetteer are fetched
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
    # against the column's own heat and oxygen, and as a mean flux over the lake: the implicit diffusion's round-off,
    # which the wind's large diffusivities across near-neutral water raise
    for group, cases in product['runs'].items():
        for key, case in cases.items():
            b = case['budget']
            assert abs(b['heat_relative_error']) < 1e-8 and abs(b['heat_error_w_m2']) < 1e-6, (group, key)
            assert b['oxygen_relative_error'] < 1e-6 and b['oxygen_error_g_m2_day'] < 1e-6, (group, key)


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
        key = (check['site'], check['run'])
        level, ten = check['lowest_level'], check['ten_metres']
        # paired as PlaSim pairs it and with PlaSim's own Earth constant, the slab follows the GCM's slab ocean
        assert abs(level['earth_free_convection']['bias_k']) < 0.5, key
        assert level['earth_free_convection']['correlation'] > 0.9, key
        for pairing in (level, ten):
            assert pairing['lunar_free_convection']['correlation'] > 0.85, key
            # free convection carried to lunar gravity leaves the water warmer than Earth's constant does
            assert pairing['lunar_free_convection']['model_mean_c'] > pairing['earth_free_convection']['model_mean_c']
        # the 10 m pairing exchanges more heat for the same air and runs the water cooler
        assert ten['earth_free_convection']['model_mean_c'] < level['earth_free_convection']['model_mean_c'], key


def test_the_answers_are_the_central_runs(product):
    for label, group in (('central', 'deep'), ('eta_low', 'deep_eta_low'), ('eta_high', 'deep_eta_high')):
        for target, a in product['answers'][label].items():
            for site, o in a['overturn'].items():
                case = product['runs'][group][f'{target}/{site}']
                assert o['complete_overturns'] == case['mixing']['complete_overturns']
                assert o['lake_number'] == case['lake_number']
            for site, entry in a['shallows_late_night'].items():
                case = product['runs']['shallow'][f'{target}/{site}']
                for trophic, v in entry.items():
                    n = case['oxygen'][f'{trophic}/renewal/supply']['nights']['surface_minimum']
                    assert v['median'] == n[1]
                    assert v['bracket_of_medians'][0] <= v['median'] <= v['bracket_of_medians'][1]
                    assert set(v['labile_brackets']) == {'labile_2d', 'labile_7d'}
                    # the warm kinetics are a transient beside the bracket, which spans the transfers and the pool
                    lows = [case['oxygen'][f'{trophic}/{tr}/supply']['nights']['surface_minimum'][1]
                            for tr in R.TRANSFER] + list(v['labile_brackets'].values())
                    assert v['bracket_of_medians'] == [min(lows), max(lows)]
                    assert len(v['warm_transient']) == 2


def test_the_wind_sets_the_deep_mixing_with_its_bracket_and_the_sensitivity_rows(product):
    runs = product['runs']
    for eta, group in (('central', 'deep'), ('low', 'deep_eta_low'), ('high', 'deep_eta_high')):
        assert set(runs[group]) == set(runs['deep'])
        for key, case in runs[group].items():
            mix = case['deep_mixing']
            assert case['case']['deep_mixing'] == mix['treatment'] == 'wind'
            assert mix['wind_share'] == R.WIND_SHARE[eta] and mix['cap_m2_s'] == R.CAP_WIND
            # the raised cap holds no interface back; the old cap held many under the same wind
            assert mix['share_held_by_the_cap'] < 1e-3, (group, key)
            assert mix['dissipation_w_kg']['mean'] > 0.0
    for key, case in runs['wind_old_cap'].items():
        assert case['deep_mixing']['cap_m2_s'] == R.CAP['central']
        assert case['deep_mixing']['share_held_by_the_cap'] > runs['deep'][f'{key}/0-15']['deep_mixing'][
            'share_held_by_the_cap']
    # a cap ten times the wind's changes nothing: it held nothing back
    for key, case in runs['wind_cap_raised'].items():
        central = runs['deep'][f'{key}/0-15']
        mix = case['deep_mixing']
        assert mix['cap_m2_s'] == R.CAP_WIND_RAISED and mix['share_held_by_the_cap'] == 0
        n = central['mixing']['complete_overturns']
        assert abs(case['mixing']['complete_overturns'] - n) <= max(1, 0.1 * n), key
        for trophic in R.TROPHIC:
            a = case['oxygen'][f'{trophic}/renewal/supply']['deepest_tenth']['mean']
            b = central['oxygen'][f'{trophic}/renewal/supply']['deepest_tenth']['mean']
            assert abs(a - b) <= 0.1, (key, trophic, a, b)
    # a larger share of the wind's energy dissipates more below the mixed layer, in proportion at the median (the mean
    # takes the spikes of the hours when a deep mixed layer leaves little water below it)
    for key, case in runs['deep'].items():
        mid = case['deep_mixing']['dissipation_w_kg']['median']
        for group, eta in (('deep_eta_low', 'low'), ('deep_eta_high', 'high')):
            other = runs[group][key]['deep_mixing']['dissipation_w_kg']['median']
            assert other / mid == pytest.approx(R.WIND_SHARE[eta] / R.WIND_SHARE['central'], rel=0.05), (group, key)
    for key, case in runs['deep_fixed'].items():
        target, band, diss = key.split('/')
        assert case['case']['deep_mixing'] == 'dissipation' and case['case']['dissipation'] == R.DISSIPATION[diss]
        assert case['case']['cap'] == R.CAP['central']
    assert {v['case']['deep_mixing'] for v in runs['deep_lake_number'].values()} == {'lake_number'}
    for group in ('deep', 'deep_eta_low', 'deep_eta_high'):
        for key, case in runs[group].items():
            ln = case['lake_number']
            assert set(ln) == set(R.LAKE_NUMBER_BELOW_M)
            for rule in ln.values():
                assert 0.0 < rule['p10'] <= rule['median'] <= rule['p90'], (group, key)
            # the month's thermocline lies deeper than the day's warm layer, and the Lake number grows with z_T
            assert ln['month']['median'] >= ln['day']['median'], (group, key)


def test_the_deep_lakes_settle_at_the_analysis_step(product):
    for group in ('deep', 'deep_eta_low', 'deep_eta_high', 'deep_lake_number', 'deep_fixed', 'wind_old_cap',
                  'wind_cap_raised', 'sensitivity'):
        for key, case in product['runs'][group].items():
            spun = case['settling']['spin_up']
            fine = spun['at_the_analysis_step']
            assert fine['step_s'] == case['case']['step_s'], (group, key)
            assert 2 <= fine['passes'] <= fine['limit'], (group, key)
            assert len(fine['bottom_age_days_by_pass']) == fine['passes']


def test_iron_reduction_nets_off_the_co2_it_consumes():
    salts = (0.125 * 55.845 + 0.25 * 61.017) * 0.8e-3          # the Fe2+ and HCO3- of a gram of O2 equivalents
    co2 = 0.25 * 44.01 * 0.284e-3                             # 7/32 consumed and 1/32 not made, at Kivu's coefficient
    assert R.IRON['per_debt'] == pytest.approx((1.0 - R.METHANE_SHARE) * (salts - co2))
    assert R.IRON['per_debt'] / ((1.0 - R.METHANE_SHARE) * salts) == pytest.approx(0.824, abs=0.005)
    assert R.IRON['cap'] == pytest.approx((0.120 * 55.845 + 0.240 * 61.017) * 0.8e-3 - 0.240 * 44.01 * 0.284e-3)


def test_scenarios_carry_their_sources():
    s = {x.name: x for x in R.scenarios()}
    assert len(s) == len(R.TROPHIC) * (len(R.TRANSFER) * len(R.KINETICS) + 2)
    assert s['oligotrophic/renewal/supply'].labile_days == R.LABILE_DAYS['central']
    assert s['oligotrophic/renewal/supply/labile_2d'].labile_days == 2.0
    assert s['oligotrophic/renewal/supply/labile_7d'].labile_days == 7.0
    assert s['eutrophic/cole_caraco/warm'].reference_c == R.SOURCE_C == 6.0
    assert s['eutrophic/renewal/supply'].reference_c == 20.0 and s['eutrophic/renewal/supply'].activation_ev == 0.0
    eu = s['eutrophic/renewal/supply']
    assert eu.gpp == pytest.approx(2.0 * 32.0 / 12.0) and eu.attenuation == pytest.approx(1.7)
    assert eu.sediment + eu.deep_water == pytest.approx(0.90)                      # Muller et al.'s mean
    reservoir = s['new_reservoir/renewal/supply']
    assert reservoir.sediment + reservoir.deep_water == pytest.approx(1.31)        # and their highest
    assert s['oligotrophic/renewal/warm'].activation_ev == pytest.approx(0.62)
    assert R.SOLUTE['co2'] == pytest.approx(44.01 / 32.0 * 0.284e-3)
    assert np.isclose(R.kd(1.0), 1.7)
    # the iron bracket: Matano's 120 uM of Fe(II) with its bicarbonate at the model's salinity term, net of the CO2
    # iron reduction consumes (test_iron_reduction_nets_off_the_co2_it_consumes)
    assert R.IRON['cap'] < (0.120 * 55.845 + 0.240 * 61.017) * 0.8e-3


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


# --- the mixing's form below the mixed layer (mixing_forms.py) -----------------------------------------------------

def test_the_energetic_regime_follows_shih_and_bouffard():
    from research.studies.lake_oxygen import mixing_forms as MF
    ph = L.Physics(deep_mixing='wind', max_diffusivity=1.0)
    eps = 1e-9
    reb = np.array([1.0, 50.0, 100.0, 1e4])
    n2 = eps / (MF.NU * reb)
    k = MF.shih_bouffard(n2, ph, eps)
    assert k[0] == 0.0                                                  # the diffusive regime: molecular alone
    assert k[1] == pytest.approx(0.2 * eps / n2[1])                     # the intermediate: Osborn's rule
    assert k[3] == pytest.approx(2.0 * MF.NU * np.sqrt(1e4))            # the energetic: 2 nu Re_b^(1/2)
    above = MF.shih_bouffard(np.array([eps / (MF.NU * 100.0 * (1 + 1e-9))]), ph, eps)[0]
    assert above == pytest.approx(k[2], rel=1e-6)                       # continuous at Re_b = 100
    # Osborn's rule over the energetic regime's is 0.1 Re_b^(1/2): ten times at Re_b = 10^4
    assert 0.2 * eps / n2[3] / k[3] == pytest.approx(0.1 * np.sqrt(1e4))
    # unstable interfaces and the cap as lake.py's
    assert MF.shih_bouffard(np.array([-1e-6]), ph, eps)[0] == L.diffusivity(np.array([-1e-6]), ph, eps)[0]
    # the wrapper restores lake.py's own functions
    with MF.mixing_form('shih_bouffard', 1e5):
        assert L.diffusivity is MF.shih_bouffard and L.wind_dissipation is not MF._WIND
    assert L.diffusivity is MF._OSBORN and L.wind_dissipation is MF._WIND


@pytest.fixture(scope='module')
def forms():
    from research.studies.lake_oxygen import mixing_forms as MF
    if not MF.OUT.is_file():
        pytest.skip('the mixing forms run has not written its product')
    return json.loads(MF.OUT.read_text())


def test_the_mixing_forms_product_binds_its_code_and_inputs(forms):
    from research.studies.lake_oxygen import mixing_forms as MF
    assert forms['schema'] == MF.SCHEMA
    for path, h in forms['producer']['files'].items():
        assert digest(path) == h, path
    for path, h in forms['producer']['inputs'].items():
        if path.startswith('research/runs/') and not present(path):
            continue
        assert digest(path) == h, path


def test_the_energetic_regime_mixes_the_deep_lakes_less(forms):
    for target, rows in forms['comparison']['lake_oxygen_convention'].items():
        o, s = rows['osborn'], rows['shih_bouffard']
        assert s['complete_overturns'] <= o['complete_overturns'], target
        assert s['bottom_age_years'] >= o['bottom_age_years'], target
        for trophic in R.TROPHIC:
            assert s['deepest_tenth'][trophic]['mean'] <= o['deepest_tenth'][trophic]['mean'] + 0.05, (target, trophic)
    for name, run in forms['runs'].items():
        b = run['budget']
        assert abs(b['heat_relative_error']) < 1e-8 and b['oxygen_relative_error'] < 1e-6, name
        assert run['deep_mixing']['share_held_by_the_cap'] < 1e-3, name
