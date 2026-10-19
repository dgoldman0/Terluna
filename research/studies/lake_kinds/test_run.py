"""The lake kinds: the product's binding, the 301 K construction, the runs' physics and the reuse of the other
studies' sizes."""
import hashlib
import json
from types import SimpleNamespace

import numpy as np
import pytest

from research.studies.lake_anoxia import run as A
from research.studies.lake_kinds import census as CE
from research.studies.lake_kinds import ladder as LD
from research.studies.lake_kinds import run as K
from research.studies.lake_kinds import setting as S
from research.studies.lake_oxygen import run as R


@pytest.fixture(scope='module')
def product():
    if not K.OUT.is_file():
        pytest.skip('the lake kinds product follows the lake model\'s corrections of 18 October')
    return json.loads(K.OUT.read_text())


def present(path):
    return (K.ROOT / path).is_file()


def digest(path):
    return hashlib.sha256((K.ROOT / path).read_bytes()).hexdigest()


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == K.SCHEMA
    for path, h in product['producer']['files'].items():
        assert digest(path) == h, path
    for path, h in product['producer']['inputs'].items():
        if (path.startswith(('research/runs/', 'geography/inputs/', 'geography/products/', 'climate/gcm/products/'))
                and not present(path)):
            continue                                     # the runs, LOLA and the grids are fetched or regenerated
        assert digest(path) == h, path


def test_the_share_of_the_way_is_climate_targets_own_rule():
    targets = json.loads(R.TARGETS.read_text())
    a, b = targets['runs']['design']['ground_k'], targets['runs']['test']['ground_k']
    i = S.interpolation(targets, 301.0)
    assert i['share'] == pytest.approx((301.0 - a) / (b - a), rel=1e-12)
    assert i['share'] == pytest.approx(0.5197, abs=1e-4)
    # at 300 K the construction is climate_targets' own share and ladder, and lake_oxygen's lapse rate
    t300 = next(t for t in targets['targets'] if t['name'] == '300 K')
    j = S.interpolation(targets, 300.0)
    assert j['share'] == t300['interpolated_share_of_the_way']
    assert [x['air_c'] for x in j['ladder']] == [x['air_c'] for x in t300['climate']['ladder']]
    by = {x['height_km']: x['air_c'] for x in t300['climate']['ladder']}
    assert j['lapse_k_per_m'] == (by[0.9] - by[5]) / 4100.0
    assert a < 300.0 < 301.0 < b and j['lapse_k_per_m'] < i['lapse_k_per_m']


def test_the_construction_reproduces_the_waters_oxygen_resort(product):
    rep = product['setting']['reproduction_at_300k']
    assert all(rep['construction'][k] for k in ('share_equal', 'ladder_equal', 'lapse_equal'))
    assert set(rep['runs']) == {'deep/0-15', 'deep/15-30', 'deep/30-45', 'shallow/0-15'}
    assert all(v['identical'] for v in rep['runs'].values())
    assert rep['runner']['identical']
    assert product['setting']['resort']['share_of_the_way'] == pytest.approx(0.5197, abs=1e-4)


def test_every_run_conserves_heat_and_oxygen_and_none_freezes(product):
    for band, depths in product['ladder'].items():
        for depth, rows in depths.items():
            for trophic, v in rows.items():
                b = v['budget']                 # against the column's own heat and oxygen, and as a mean flux
                assert abs(b['heat_relative_error']) < 1e-8 and abs(b['heat_error_w_m2']) < 1e-6, (band, depth, trophic)
                assert b['oxygen_relative_error'] < 1e-6 and b['oxygen_error_g_m2_day'] < 1e-6, (band, depth, trophic)
                assert v['ice_max_m'] == 0.0
    for group in ('giants', 'meromictic'):
        for name, v in product[group].items():
            b = v['budget']
            assert abs(b['heat_relative_error']) < 1e-8 and abs(b['heat_error_w_m2']) < 1e-6, name
            assert b['oxygen_relative_error'] < 1e-6 and b['oxygen_error_g_m2_day'] < 1e-6, name
            assert v['ice_max_m'] == 0.0, name


def test_the_ladder_holds_every_case(product):
    assert set(product['ladder']) == set(LD.BAND_NAMES)
    for band, depths in product['ladder'].items():
        assert set(depths) == {f'{d:g}' for d in LD.LADDER_M}
        for rows in depths.values():
            assert set(rows) == set(LD.TROPHIC)
    assert set(product['ladder_eta']) == set(LD.ETAS)
    for eta, ladder in product['ladder_eta'].items():
        for band, depths in ladder.items():
            for depth, rows in depths.items():
                assert set(rows) == (set(LD.ETA_TROPHIC) if float(depth) in LD.ETA_LADDER_M else set()), depth
    giants = [n for n in LD.cases() if n.startswith('giant_')]
    assert set(product['giants']) == set(giants)
    assert len(product['meromictic']) == len(LD.meromictic_cases())
    # each giant passed over the record to its regime, in coarse steps and at the analysis step, before the fifty
    # years it saves, and each ladder lake settled at the analysis step before its analysed decade
    for name, g in product['giants'].items():
        assert g['spin_up']['passes'] >= 2 and g['years'] >= 48.0, name
        assert g['spin_up']['at_the_analysis_step']['passes'] >= 2, name
    ladders = [product['ladder']] + list(product['ladder_eta'].values())
    for ladder in ladders:
        for band, depths in ladder.items():
            for depth, rows in depths.items():
                for trophic, v in rows.items():
                    fine = v['passes_at_the_analysis_step']['passes']
                    if float(depth) >= LD.FINE_MIN_M:
                        assert fine >= 2, (band, depth, trophic)
                    else:                              # the shallower rungs take the quicker path
                        assert fine is None, (band, depth, trophic)


def test_the_wind_sets_the_deep_mixing_and_the_cap_holds_none_back(product):
    for eta, ladder in [('central', product['ladder'])] + list(product['ladder_eta'].items()):
        for band, depths in ladder.items():
            for depth, rows in depths.items():
                for trophic, v in rows.items():
                    mix = v['deep_mixing']
                    assert mix['treatment'] == 'wind' and mix['wind_share'] == R.WIND_SHARE[eta]
                    assert mix['share_held_by_the_cap'] < 1e-3, (eta, band, depth, trophic)
    for name, g in product['giants'].items():
        assert g['deep_mixing']['wind_share'] == R.WIND_SHARE[g['eta']]
        assert g['deep_mixing']['share_held_by_the_cap'] < 1e-3, name
    # the 352 m lakes dissipate more below the mixed layer as the wind's share grows, in proportion at the median (the
    # mean takes the spikes of the hours when a deep mixed layer leaves little water below it)
    for band in LD.BAND_NAMES:
        for trophic in LD.ETA_TROPHIC:
            mid = product['ladder'][band]['352'][trophic]['deep_mixing']['dissipation_w_kg']['median']
            for eta in LD.ETAS:
                other = product['ladder_eta'][eta][band]['352'][trophic]['deep_mixing']['dissipation_w_kg']['median']
                assert other / mid == pytest.approx(R.WIND_SHARE[eta] / R.WIND_SHARE['central'], rel=0.05), \
                    (band, trophic, eta)


def test_deeper_floors_are_mixed_through_on_fewer_nights(product):
    for band, depths in product['ladder'].items():
        for trophic in LD.TROPHIC:
            shares = [depths[f'{d:g}'][trophic]['nights_mixing_to_the_floor'] for d in LD.LADDER_M]
            assert all(a >= b - 0.02 for a, b in zip(shares, shares[1:])), (band, trophic, shares)
            floor = [depths[f'{d:g}'][trophic]['floor_ventilated_in_half_the_nights'] for d in LD.LADDER_M]
            assert floor[0] >= floor[-1], (band, trophic)


def test_the_ventilated_depth_is_where_the_nights_stop_reaching_the_floor(product):
    thr = product['thresholds']
    for band in LD.BAND_NAMES:
        for trophic in LD.TROPHIC:
            shares = [product['ladder'][band][f'{d:g}'][trophic]['nights_mixing_to_the_floor'] for d in LD.LADDER_M]
            z = thr['ventilated_depth_m'][band][trophic]
            for key, level in K.SHARE_OF_NIGHTS.items():
                assert z[key] == pytest.approx(K.floor_crossing(LD.LADDER_M, shares, level)[0], rel=1e-3)
            # more of the nights reach the shallower floors
            assert z['three_quarters'] <= z['half'] * (1 + 1e-3) and z['half'] <= z['quarter'] * (1 + 1e-3)
            # the 352 m lakes' ordinary nights stop short of their floor, and every night's median takes the overturns
            row = product['ladder'][band]['352'][trophic]
            ordinary = thr['ordinary_night_mixing_352m_m'][band][trophic]
            assert ordinary['median'] == (row['night_mixing_short_of_the_floor_m'] or row['night_mixing_m'])[2]
            assert ordinary['median'] <= row['night_mixing_m'][2] + 1e-9
        # clearer water lets the ordinary nights mix deeper
        o = thr['ordinary_night_mixing_352m_m'][band]
        assert o['oligotrophic']['median'] >= o['eutrophic']['median']
    census_z = product['census']['ventilated_depth_m']
    for band in LD.BAND_NAMES:
        assert census_z[band] == pytest.approx(thr['ventilated_depth_m'][band][K.CENTRAL_TROPHIC]['half'], rel=1e-3)
    # the wind's share below the mixed layer moves the overturns, while the mixed layer's own energy balance sets the
    # ordinary night's depth
    low, high = product['ladder_eta']['low'], product['ladder_eta']['high']
    for band in LD.BAND_NAMES:
        for trophic in LD.ETA_TROPHIC:
            share = [x[band]['352'][trophic]['complete_overturn_share'] for x in (low, product['ladder'], high)]
            assert share[0] <= share[1] + 0.02 and share[1] <= share[2] + 0.02, (band, trophic, share)
            mid = thr['ordinary_night_mixing_352m_m'][band][trophic]['median']
            for eta in LD.ETAS:
                other = product['thresholds_eta'][eta]['ordinary_night_mixing_352m_m'][band][trophic]['median']
                assert abs(other - mid) < 3.0, (band, trophic, eta)


def test_the_census_counts_every_body_once(product):
    cen = product['census']
    summary = json.loads(CE.DRAINAGE_SUMMARY.read_text())['lakes']
    assert cen['total']['bodies'] == summary['bodies']
    assert cen['total']['closed_bodies'] <= summary['closed_bodies']        # the mostly closed of the summary's
    assert sum(v['bodies'] for v in cen['kinds'].values()) == summary['bodies']
    assert sum(v['area_km2'] for v in cen['kinds'].values()) == pytest.approx(cen['total']['area_km2'], rel=1e-3)
    assert sum(v['closed_bodies'] for v in cen['kinds'].values()) == cen['total']['closed_bodies']
    for kind, v in cen['kinds'].items():
        assert sum(x['bodies'] for x in cen['by_band'][kind].values()) == v['bodies']
        assert sum(x['bodies'] for x in cen['by_size'][kind].values()) == v['bodies']
        assert sum(cen['floor_by_depth_km2'][kind]) == pytest.approx(v['area_km2'], rel=2e-3, abs=1.0)
        assert 0.0 <= v['ventilated_floor_share'] <= 1.0 + 1e-6
    assert cen['kinds']['self_ventilated']['ventilated_floor_share'] == pytest.approx(1.0, abs=1e-3)
    # the deeper the kind, the less of its floor lies within the ventilated depth
    assert cen['kinds']['giant']['ventilated_floor_share'] <= cen['kinds']['deep']['ventilated_floor_share'] + 1e-9


def test_circulation_reuses_lake_anoxias_sizes():
    anoxia = json.loads(A.OUT.read_text())
    for depth in ('352', '1925'):
        theirs = anoxia['remedies']['targets']['resort']['circulation'][depth]['kw_per_km2']
        assert K.circulation_kw_per_km2(float(depth)) == pytest.approx(theirs, rel=1e-3)


def test_isolated_hypolimnion_matches_lake_anoxias_hazard_map():
    from atmosphere.radiative_convective.thermodynamics import earthlike_air
    anoxia = json.loads(A.OUT.read_text())
    oxygen = json.loads(R.OUT.read_text())['setting']
    chem = json.loads(K.SEA_CHEMISTRY.read_text())
    held = next(r for r in chem['runs'] if all(r[k] == v for k, v in R.HELD_FRESH.items()))
    dry = earthlike_air(R.DESIGN_PRESSURE_PA, A.DESIGN_CO2_PPM)
    pco2 = A.DESIGN_CO2_PPM * 1e-6 * A.SURFACE_ATM
    gas = SimpleNamespace(gravity=oxygen['gravity'], salinity=oxygen['salinity_g_kg'],
                          o2_fraction=oxygen['o2_fraction'], n2_frac=float(dry.fractions['N2']),
                          ar_frac=float(dry.fractions['Ar']), pco2=pco2,
                          ta=A.G.alkalinity_at(float(held['rows'][-1]['lake_ph']), pco2, 21.15, A.LAKE_IONIC)[0])
    deep = anoxia['remedies']['targets']['resort']['deep_t_c']
    for depth in ('30', '50', '352'):
        for trophic in ('oligotrophic', 'eutrophic'):
            theirs = anoxia['hazard_map'][depth]['traits'][trophic]['moon']['resort']
            assert K.bubble_years(gas, float(depth), deep, trophic) == pytest.approx(theirs, rel=2e-3)


def test_the_moorings_add_up(product):
    m = product['needs']['monitoring']
    assert sum(m['by_kind'].values()) == m['moorings']
    assert sum(m['by_size'].values()) == m['moorings']
    assert m['by_kind']['self_ventilated'] == 0


def test_overturn_events_are_counted_once(product):
    for band, depths in product['overturn_events'].items():
        for depth, v in depths.items():
            assert v['events_below_2'] <= v['events']
            assert len(v['surface_min']) == v['events_below_2'] == len(v['hours_below_2'])
            assert all(x < 2.0 for x in v['surface_min'])
            ladder = product['ladder'][band][depth]['oligotrophic']['after_overturns']
            assert ladder['events'] == v['events'] and ladder['events_below_2'] == v['events_below_2']


def test_the_permissible_load_is_lake_anoxias_at_its_depth():
    anoxia = json.loads(A.OUT.read_text())
    p = anoxia['phosphorus']
    z = float(np.mean([v['mean_depth_m'] for v in anoxia['deep_lakes'].values()]))      # the 352 m lakes' mean depth
    for key, f in A.SETTLING.items():
        assert K.permissible_load(z, p['hydraulic_load_m_yr'], f) == pytest.approx(
            p['permissible_load_g_m2_yr'][key]['oligotrophic'], rel=2e-3)
    assert K.permissible_load(3.0, p['hydraulic_load_m_yr']) < K.permissible_load(100.0, p['hydraulic_load_m_yr'])
    # anoxic floors can at most cancel the retention at steady state: lake_anoxia's ceiling, at any depth
    for c, v in p['without_retention']['permissible_load_g_m2_yr'].items():
        if c == 'oligotrophic':
            for depth in (3.0, 112.0, 1000.0):
                assert K.permissible_load(depth, p['hydraulic_load_m_yr'], K.RETENTION['anoxic_floors']) == \
                    pytest.approx(v, rel=1e-9)

def test_each_kind_carries_the_phosphorus_rule_with_its_retention_brackets(product):
    anoxia = json.loads(A.OUT.read_text())
    q_s = anoxia['phosphorus']['hydraulic_load_m_yr']
    for kind in K.KINDS:
        rule = product['needs'][kind]['oligotrophic_load_g_p_m2_yr']
        assert set(rule) == set(K.RETENTION)
        # with the retention cancelled the load that holds 10 mg/m3 is the throughflow's, whatever the depth
        assert rule['anoxic_floors'] == pytest.approx(A.OECD_TP['oligotrophic'] * q_s / 1000.0, rel=1e-3)
        if not product['needs'][kind]['bodies']:
            continue                                      # an empty kind has no depth for its retention
        assert rule['anoxic_floors'] < rule['stokes'] < rule['aggregates'] < rule['earth']


def test_the_rungs_bracket_the_self_ventilated_kind(product):
    cen = product['census']
    rungs = cen['within_each_rung_by_band']
    by_band = cen['by_band']['self_ventilated']
    for band, rows in rungs.items():
        if band == 'all':
            continue
        # the ventilated depth lies between the 100 m and 352 m rungs, so the census there brackets the kind
        assert rows['100']['area_km2'] <= by_band[band]['area_km2'] <= rows['352']['area_km2'], band
        assert all(rows[f'{a:g}']['bodies'] <= rows[f'{b:g}']['bodies'] for a, b in zip(LD.LADDER_M, LD.LADDER_M[1:]))
    for rung in rungs['all']:
        assert rungs['all'][rung]['bodies'] == sum(v[rung]['bodies'] for k, v in rungs.items() if k != 'all')
    assert rungs['all']['352']['area_km2'] >= cen['kinds']['self_ventilated']['area_km2'] >= rungs['all']['100'][
        'area_km2']


def test_the_springs_and_the_oxygenation_read_what_they_say(product):
    rows = {k for v in product['spring_water'].values() for k in v['at']}
    assert len(rows) == 1                                   # one temperature, the 352 m lakes' deep water's
    t_deep = float(rows.pop())
    lad = product['ladder']
    assert t_deep == pytest.approx(np.mean([lad[b]['352'][K.CENTRAL_TROPHIC]['bottom_t_c']['mean']
                                            for b in LD.BAND_NAMES]), abs=0.01)     # to the product's rounding
    excess = [v['at'][f'{t_deep:g}']['density_excess_kg_m3'] for v in product['spring_water'].values()]
    assert product['spring_fed']['density_excess_kg_m3'] == [min(e[0] for e in excess), max(e[1] for e in excess)]
    oxy = product['needs']['oxygenation_for_deep_water_life']
    kinds = product['census']['kinds']
    assert oxy['deep_and_giant_lakes_area_km2'] == pytest.approx(kinds['deep']['area_km2'] + kinds['giant']['area_km2'],
                                                                 rel=1e-3)
    built = product['needs']['built_lakes']
    assert built['first_decade_oxygen_delivered_t_km2_day'][0] == pytest.approx(
        built['first_decade_oxygen_consumed_t_km2_day'])
