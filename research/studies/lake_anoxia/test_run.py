"""What anoxic depths do: the product's binding, the gas physics it rests on and what its numbers must keep."""
import hashlib
import json

import numpy as np
import pytest

from biosphere.ecology import lake as L
from research.studies.lake_anoxia import gases as G
from research.studies.lake_anoxia import run as A


@pytest.fixture(scope='module')
def product():
    return json.loads(A.OUT.read_text())


def digest(path):
    return hashlib.sha256((A.ROOT / path).read_bytes()).hexdigest()


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == A.SCHEMA
    for path, h in product['producer']['files'].items():
        assert digest(path) == h, path
    for path, h in product['producer']['inputs'].items():
        if path.startswith('research/runs/') and not (A.ROOT / path).is_file():
            continue                                     # the saved columns live on the research drive
        assert digest(path) == h, path


def test_solubilities_are_the_standard_values():
    assert float(G.solubility('CH4', 20.0)) == pytest.approx(1.55, abs=0.02)        # mmol/L/atm
    assert float(G.solubility('N2', 20.0)) == pytest.approx(0.69, abs=0.02)
    assert float(G.solubility('Ar', 25.0)) == pytest.approx(1.40, abs=0.03)
    assert float(G.solubility('CO2', 25.0)) == pytest.approx(33.9, abs=0.5)
    o2 = float(G.solubility('O2', 20.0)) * 0.20946 * (1.0 - float(L.vapour_pressure_atm(20.0))) * 32.0 / 1000.0
    assert o2 * 1000.0 == pytest.approx(float(L.oxygen_saturation(20.0)), rel=0.01)  # Weiss against Benson & Krause


def test_ebullition_share_reproduces_lake_soppen():
    # Langenegger et al. (2019): at 10 m and 5 C, production near 8 mmol m-2 d-1 leaves 10-17% as bubbles
    share = G.ebullition_share(8.0 / 86400e3, 10.0, 5.0, 0.944 / 1.01325, 9.81, 0.78 * 0.944 / 1.01325)
    assert 0.08 < share < 0.2
    deeper = G.ebullition_share(8.0 / 86400e3, 20.0, 5.0, 0.944 / 1.01325, 9.81, 0.78 * 0.944 / 1.01325)
    assert deeper < share
    assert G.ebullition_share(0.5 / 86400e3, 10.0, 5.0, 0.944 / 1.01325, 9.81, 0.74) == 0.0
    # under the Moon's own air, 1.2 atm with its nitrogen at 0.95 atm against Earth's 1 atm and 0.76, a floor bubbles
    # as Earth's floor at 10 m does at about six times that depth: the extra surface pressure and the extra nitrogen
    # nearly cancel, and the weaker gravity spreads the same pressure over six times the depth
    from scipy.optimize import brentq
    earth = G.ebullition_share(8.0 / 86400e3, 10.0, 20.0, 1.0, 9.80665, 0.76)
    moon_at = lambda z: G.ebullition_share(8.0 / 86400e3, z, 20.0, 1.2, 1.6242, 0.95) - earth
    depth = brentq(moon_at, 20.0, 200.0)
    assert 5.5 < depth / 10.0 < 7.0


def test_bubbles_lose_more_methane_the_deeper_they_start(product):
    s = product['survival']
    for target, rows in s['moon'].items():
        for key, row in rows.items():
            assert all(a >= b - 1e-9 for a, b in zip(row, row[1:])), (target, key)
        assert all(c >= r - 1e-9 for c, r in zip(rows['large/rigid'], rows['small/rigid']))
    assert s['earth_20c']['5mm/rigid'][2] > 0.85                 # most of a 5 mm bubble from 10 m reaches the air
    assert G.rise_velocity(0.005, 20.0, 1.6242) < G.rise_velocity(0.005, 20.0, 9.81)


def test_carbonate_round_trip():
    ta, dic = G.alkalinity_at(7.72, 400e-6 * 1.2, 21.15, A.LAKE_IONIC)
    ph, co2 = G.speciate(dic, ta, 21.15, A.LAKE_IONIC)
    assert float(ph) == pytest.approx(7.72, abs=1e-4)
    assert float(co2) == pytest.approx(float(G.co2_kh(21.15)) * 400e-6 * 1.2, rel=1e-4)
    assert 0.3e-3 < ta < 0.6e-3                                   # a soft lake


def test_an_isolated_hypolimnion_reaches_its_bubble_point_sooner_on_the_moon(product):
    rows = product['hazard_map']
    depths = sorted(rows, key=float)
    for trophic in A.TRAITS:
        for target in A.TARGET_ORDER:
            years = [rows[d]['traits'][trophic]['moon'][target] for d in depths]
            assert all(a < b for a, b in zip(years, years[1:]))
        for d in depths:
            r = rows[d]['traits'][trophic]
            assert max(r['moon'].values()) < r['earth_at_25c']


def test_the_deep_lakes_stay_below_their_bubble_point(product):
    for name, lake in product['deep_lakes'].items():
        for trophic, v in lake['traits'].items():
            assert v['saturation_ratio_bottom_max'] < 0.5, (name, trophic)
            assert v['saturation_ratio_max'] < 1.0, (name, trophic)


def test_methane_rises_with_the_trophic_state_and_scales_with_the_area(product):
    area = product['moon_lakes']['area_km2'] * 1e6
    for target, v in product['methane'].items():
        for rule in A.RELEASE_RULES:
            central = [v['traits'][t]['by_release_rule'][rule]['moon_tg_yr']['central'] for t in A.TRAITS]
            assert all(a < b for a, b in zip(central, central[1:])), (target, rule)
    # the brackets hold the central flux at each wind share, the bubbles' as well as the water's, and the stock
    # rule releases at least what the net rise does
    for group in (product, product['eta_low'], product['eta_high']):
        for target, v in group['methane'].items():
            for t in A.TRAITS:
                ebul = v['traits'][t]['ebullition_mol_m2_yr']
                assert ebul['low'] <= ebul['central'] <= ebul['high'], (target, t)
                rules = v['traits'][t]['by_release_rule']
                for rule, row in rules.items():
                    assert row['moon_tg_yr']['central'] == pytest.approx(
                        row['flux_g_m2_yr']['central'] * area / 1e12, rel=1e-3)
                    for key in ('moon_tg_yr', 'pelagic_mol_m2_yr'):
                        assert row[key]['low'] <= row[key]['central'] <= row[key]['high'], (target, t, rule, key)
                assert rules['stock']['released_at_overturns_mol_m2_yr'] >= \
                    rules['net_rise']['released_at_overturns_mol_m2_yr'] - 1e-12, (target, t)
                assert rules['stock']['pelagic_mol_m2_yr']['central'] >= \
                    rules['net_rise']['pelagic_mol_m2_yr']['central'] - 1e-12, (target, t)
    shares = sum(c['area_share'] for c in product['moon_lakes']['classes'].values())
    assert shares == pytest.approx(1.0, abs=2e-4)                   # four figures each


def test_uncut_lakes_hold_more_methane_than_cut_ones(product):
    assert set(product['air']) == set(A.RELEASE_RULES)
    for rule, air in product['air'].items():
        for target, rows in air['targets'].items():
            for trophic, row in rows.items():
                c = row['central']
                assert all(u > k for u, k in zip(c['lakes_uncut']['steady_ppm_by_soils'], c['steady_ppm_by_soils']))
                other = product['air']['stock' if rule == 'net_rise' else 'net_rise']['targets'][target][trophic]
                if rule == 'stock':                    # the stock rule's air holds at least the net rise's methane
                    assert all(a >= b - 1e-9 for a, b in zip(c['lakes_uncut']['steady_ppm_by_soils'],
                                                             other['central']['lakes_uncut']['steady_ppm_by_soils']))


def test_the_oxygenation_bounds_are_fixed_by_supply(product):
    for group in (product, product['eta_low'], product['eta_high']):
        for target, rem in group['remedies']['targets'].items():
            for trophic, row in rem['traits'].items():
                p = A.R.TROPHIC[trophic]
                demand = p['sediment'] + p['deep_water']
                lo, hi = row['oxygenation']['o2_t_km2_day']
                assert lo == pytest.approx(demand * A.INDUCED_DEMAND[0]) and hi == pytest.approx(
                    demand * A.INDUCED_DEMAND[1])
                # the warm kinetics' transient stands beside the bounds, above them in water warmer than 6 C
                assert row['warm_transient']['o2_t_km2_day'] == pytest.approx(hi * rem['warm_factor'], rel=1e-3)
                assert rem['warm_factor'] > 1.0


def test_the_dissolved_loads_phosphorus_under_each_settling(product):
    p = product['phosphorus']
    tp = p['tp_from_dissolved_load_mg_m3']
    # with Earth's settling the dissolved load keeps the lakes ultra-oligotrophic; slower settling at the Moon's
    # gravity raises it towards the oligotrophic limit
    assert max(tp['earth']['lakes_352m']) < A.OECD_TP['ultraoligotrophic']
    assert min(tp['stokes']['lakes_352m']) < A.OECD_TP['oligotrophic']
    for i in range(2):
        assert tp['earth']['lakes_352m'][i] < tp['aggregates']['lakes_352m'][i] < tp['stokes']['lakes_352m'][i]
    permissible = p['permissible_load_g_m2_yr']
    assert permissible['earth']['oligotrophic'] > permissible['aggregates']['oligotrophic'] > \
        permissible['stokes']['oligotrophic']
    for target, rows in p['targets'].items():
        # Nurnberg's release over the floor's anoxic days, a tenth of it to all of it for the low sulphate
        for trophic, row in rows.items():
            expected = [A.P_RELEASE[trophic] * row['anoxic_factor_days_yr'] * f / 1000.0
                        for f in A.LOW_SULPHATE_FACTOR]
            assert row['internal_load_g_m2_yr'] == pytest.approx(expected, rel=2e-3, abs=1e-6), (target, trophic)
        # the internal load is added outside the retention (Nurnberg 1984)
        for key, f in A.SETTLING.items():
            q = p['hydraulic_load_m_yr']
            ext = tp[key]['lakes_352m'][0]
            internal = rows['eutrophic']['internal_load_g_m2_yr'][0]
            assert rows['eutrophic']['tp_with_internal_mg_m3'][key][0] == pytest.approx(ext + 1000.0 * internal / q,
                                                                                       rel=2e-3)
    # the warmer targets leave the richer floors anoxic for longer
    for group in (product, product['eta_low'], product['eta_high']):
        rows = group['phosphorus']['targets']
        for trophic in ('eutrophic', 'new_reservoir'):
            af = [rows[t][trophic]['anoxic_factor_days_yr'] for t in ('lock', 'resort', 'hot')]
            assert af[0] <= af[1] + 1e-6 <= af[2] + 2e-6, (trophic, af)
    # at steady state a lake keeps at most what flows in, whatever its floors release: the retention cancelled
    ceiling = p['without_retention']
    q = p['hydraulic_load_m_yr']
    for i, load in enumerate(p['loads_g_m2_yr']['dissolved']):
        assert ceiling['tp_from_dissolved_load_mg_m3'][i] == pytest.approx(1000.0 * load / q, rel=1e-3)
        for key in A.SETTLING:                                   # to the product's four figures
            ratio = ceiling['tp_from_dissolved_load_mg_m3'][i] / tp[key]['lakes_352m'][i]
            assert ratio == pytest.approx(ceiling['over_the_oxic_lake'][key], rel=1e-3)
    assert ceiling['over_the_oxic_lake']['earth'] > ceiling['over_the_oxic_lake']['aggregates'] > \
        ceiling['over_the_oxic_lake']['stokes'] > 1.0
    for target, rows in p['targets'].items():
        for trophic, row in rows.items():
            for key in A.SETTLING:
                steady = row['tp_with_internal_at_steady_state_mg_m3'][key]
                assert all(a <= b + 1e-9 for a, b in zip(steady, row['tp_with_internal_mg_m3'][key]))
                assert max(steady) <= max(ceiling['tp_from_dissolved_load_mg_m3']) + 1e-9


def test_each_overturn_releases_within_the_two_rules(product):
    for group in (product, product['eta_low'], product['eta_high']):
        for name, lake in group['deep_lakes'].items():
            for trophic, v in lake['traits'].items():
                before, gaps = v['inventory_before_overturns_mol_m2'], v['days_since_the_previous_overturn']
                made = v['methane_made_mol_m2_yr']['central']
                rise, stock = v['released_at_overturns_mol_m2']['net_rise'], v['released_at_overturns_mol_m2']['stock']
                assert len(rise) == len(stock) == len(before) == len(gaps) == lake['overturn_events']
                for r, k, b, g in zip(rise, stock, before, gaps):
                    # neither rule releases more than the stock there is or the methane made since the last overturn,
                    # to the product's four figures
                    tol = 2e-3
                    assert r <= k * (1 + tol) + 1e-12 and k <= b * (1 + tol) + 1e-9, (name, trophic)
                    assert k <= made * g / 365.25 * (1 + tol) + 1e-9, (name, trophic)
                assert sum(stock) <= made * lake['years'] * 1.001 + 1e-9, (name, trophic)


def test_the_oxycline_share_is_kivus_and_the_giants_lower_the_flux(product):
    assert A.OXYCLINE_ESCAPE['central'] == pytest.approx(0.2 / 35.0)
    for group in (product, product['eta_low'], product['eta_high']):
        for target, v in group['methane'].items():
            for trophic, row in v['traits'].items():
                for rule, r in row['by_release_rule'].items():
                    giants, plain = r['moon_tg_yr_giants_on_the_oxycline'], r['moon_tg_yr']
                    assert giants['central'] <= plain['central'] * 1.0001, (target, trophic, rule)


def test_the_wind_shares_stand_side_by_side(product):
    for eta in ('low', 'high'):
        group = product[f'eta_{eta}']
        for key in ('deep_lakes', 'deepest', 'methane', 'air', 'phosphorus', 'remedies', 'degassing', 'hazard_map'):
            assert key in group and key in product
        assert {n.rsplit(f'_eta_{eta}', 1)[0] for n in group['deep_lakes']} == set(product['deep_lakes'])
    assert product['setting']['wind_share'] == {e: A.R.WIND_SHARE[e] for e in A.ETAS}
