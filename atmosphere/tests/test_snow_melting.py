"""Melting of snowflakes below the freezing level: the onset, the heat and vapour bookkeeping, the melting flake's
geometry and the budgets' conservation, on synthetic soundings."""
import math

import numpy as np
import pytest

from atmosphere.snow import melting as ml
from shared.constants import STANDARD_GRAVITY


def earth(rh, lapse=6e-3):
    return ml.straight_sounding(STANDARD_GRAVITY, ml.EARTH_AIR, 1500., 85000., lapse, rh)


def test_saturated_air_melts_from_the_freezing_level():
    col = earth(1.0)
    speeds = ml.DropSpeeds(col, np.linspace(0., 1500., 4))
    r = ml.follow(0.01, col, 1500., speeds)
    assert r['onset_depth_m'] == pytest.approx(0., abs=5.)
    assert 0 < r['half_melted_depth_m'] < r['melted_depth_m'] < 1500.


@pytest.mark.parametrize('rh,depth', sorted(ml.MATSUO_SASYO_ONSET_M.items()))
def test_onset_in_subsaturated_air_matches_matsuo_and_sasyo(rh, depth):
    col = earth(rh)
    speeds = ml.DropSpeeds(col, np.linspace(0., 1500., 4))
    r = ml.follow(0.02, col, 1500., speeds)
    assert r['onset_depth_m'] == pytest.approx(depth, rel=0.06)
    assert r['onset_temperature_c'] == pytest.approx(ml.heymsfield_onset_c(rh, 850.), abs=0.2)


def test_heat_and_vapour_bookkeeping_close():
    col = earth(0.9)
    speeds = ml.DropSpeeds(col, np.linspace(0., 1500., 4))
    m = ml.aggregate_mass(0.02)
    z = 1500. - 400.                                     # below the onset, while melting
    for water_share in (0.05, 0.4, 0.9):
        r = ml.rates(m * (1 - water_share), m * water_share, z, col, speeds, dict(capacitance=1., speed=1.))
        assert ml.LF * r['melt'] + ml.LS0 * r['sub'] + ml.LV0 * r['evap'] == pytest.approx(r['heat'], rel=1e-9)
        assert r['sub'] + r['evap'] == pytest.approx(r['vapour'], rel=1e-9)
        assert r['sub'] / r['vapour'] == pytest.approx(1 - water_share, rel=1e-9)


def test_dry_flake_takes_mitra_geometry_and_mcsnow_speed():
    col = earth(1.0)
    speeds = ml.DropSpeeds(col, np.linspace(0., 1500., 4))
    d = 0.05
    m = ml.aggregate_mass(d)
    st = ml.flake_state(m, 0., 1500., col, speeds, dict(capacitance=1., speed=1.))
    assert st['d_m'] == pytest.approx(d, rel=1e-9)
    e = math.sqrt(1 - 0.09)
    assert st['c'] == pytest.approx(0.8 * d / 2 * e / math.asin(e), rel=1e-9)
    assert st['v'] == pytest.approx(ml.dry_speed(m, col.state(1500.), STANDARD_GRAVITY)[0], rel=1e-9)
    wet = ml.flake_state(m * 1e-6, m * (1 - 1e-6), 1500., col, speeds, dict(capacitance=1., speed=1.))
    assert wet['d_m'] == pytest.approx(wet['d_eq'], rel=1e-3)
    assert wet['v'] == pytest.approx(wet['v_r'], rel=1e-3)


def test_budgets_conserve_mass_through_melting_and_rain():
    col = earth(0.85)
    table = ml.response_table(col, 1500., spans=(0.002, 0.01, 0.03))
    s = table['stack']
    total = s['ice'] + s['water'] + s['sub'] + s['evap']
    assert np.allclose(total, 1., atol=2e-4)
    assert np.all(np.diff(s['melt'], axis=1) >= -1e-9)
    pop = ml.population(table, table['masses'], np.ones(3))
    assert np.allclose(pop['ice'] + pop['water'] + pop['sub'] + pop['evap'], 1., atol=2e-4)


def test_heymsfield_fits_are_reproduced_at_their_pressures():
    assert ml.heymsfield_onset_c(0.8734, 1000.) == pytest.approx(1.0, rel=1e-3)
    assert ml.heymsfield_onset_c(0.907, 500.) == pytest.approx(1.0, rel=1e-3)


def _product():
    import json
    if not ml.PRODUCT.exists():
        pytest.skip('results/snow_melting.json not built')
    return json.loads(ml.PRODUCT.read_text())


def test_product_budgets_close_and_depths_are_ordered():
    product = _product()
    assert product['schema'] == ml.SCHEMA
    for box in product['moon']['boxes'].values():
        for comp in box['composites'].values():
            for band in comp['bands'].values():
                for pop in band['populations'].values():
                    t = pop['totals']
                    assert t['sublimated'] + t['evaporated'] + t['ground_water'] + t['ground_ice'] == pytest.approx(1., abs=2e-3)
                    marks = [pop['depths_ice_gone_m'][k] for k in ('10', '50', '90', '99')]
                    if None not in marks:
                        assert marks == sorted(marks)
                for row in band['flakes']:
                    if row['onset_depth_m'] is not None and row['melted_depth_m'] is not None:
                        assert row['onset_depth_m'] <= row['half_melted_depth_m'] <= row['melted_depth_m']


def test_product_earth_check_holds():
    product = _product()
    for row in product['earth_check']['onset_comparison']:
        assert row['model_onset_depth_m'] == pytest.approx(row['matsuo_sasyo_m'], rel=0.06)
        assert row['model_onset_c'] == pytest.approx(row['heymsfield_fit_c'], abs=0.2)
