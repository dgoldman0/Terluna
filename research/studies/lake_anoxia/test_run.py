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
    # at lunar gravity the same floor bubbles to about six times the depth
    earth = G.ebullition_share(8.0 / 86400e3, 10.0, 20.0, 1.0, 9.81, 0.76)
    moon = G.ebullition_share(8.0 / 86400e3, 10.0 * 9.81 / 1.6242, 20.0, 1.0, 1.6242, 0.76)
    assert moon == pytest.approx(earth, rel=1e-6)


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
        central = [v['traits'][t]['moon_tg_yr']['central'] for t in A.TRAITS]
        assert all(a < b for a, b in zip(central, central[1:])), target
        for t in A.TRAITS:
            row = v['traits'][t]
            assert row['moon_tg_yr']['central'] == pytest.approx(row['flux_g_m2_yr']['central'] * area / 1e12, rel=1e-3)
            assert row['moon_tg_yr']['low'] <= row['moon_tg_yr']['central'] <= row['moon_tg_yr']['high']
    shares = sum(c['area_share'] for c in product['moon_lakes']['classes'].values())
    assert shares == pytest.approx(1.0, abs=2e-4)                   # four figures each


def test_uncut_lakes_hold_more_methane_than_cut_ones(product):
    for target, rows in product['air']['targets'].items():
        for trophic, row in rows.items():
            c = row['central']
            assert all(u > k for u, k in zip(c['lakes_uncut']['steady_ppm_by_soils'], c['steady_ppm_by_soils']))


def test_the_dissolved_load_alone_keeps_the_lakes_oligotrophic(product):
    p = product['phosphorus']
    assert max(p['tp_from_dissolved_load_mg_m3']['lakes_352m']) < A.OECD_TP['ultraoligotrophic']
    assert max(p['loads_g_m2_yr']['dissolved']) < p['permissible_load_g_m2_yr']['oligotrophic']
    for target, rows in p['targets'].items():
        assert rows['eutrophic']['internal_load_g_m2_yr'][0] > max(p['loads_g_m2_yr']['dissolved'])
