"""Rain, vapour, droplet and trim balances for the water component."""
import json
import math

import numpy as np
import pytest

from research.studies.aerophytes import environment as env
from research.studies.aerophytes import water as w


def test_vapour_density_near_room_temperature():
    assert w.vapour_density(20.) == pytest.approx(.01729, rel=5e-3)
    assert w.vapour_density(20., .5) == pytest.approx(.5*w.vapour_density(20.))


def test_need_falls_with_concentration_and_humidity():
    p = 100108.
    ref = w.daily_need(1.83, 14.2, 5., .6, p)
    ccm = w.daily_need(1.83, 14.2, 5., .6, p, ci_ca=.3)
    assert ccm['water_kg_m2_year'] == pytest.approx(ref['water_kg_m2_year']*(1-.7)/(1-.3))
    humid = w.daily_need(1.83, 14.2, 0., .9, p)
    assert humid['water_kg_m2_year'] < ref['water_kg_m2_year']/5
    assert ref['water_kg_m2_day'] == pytest.approx(ref['water_kg_m2_year']/w.JULIAN_YEAR_DAYS)


def test_zonal_rain_reproduces_the_gcm_moon_mean():
    clim = np.load(w.ROOT/w.INPUT_FILES[0])
    zonal = w.zonal_rain(clim['pr_mm_day'], clim['lat'])
    conv = json.loads((w.ROOT/'climate/results/gcm/convection_A28_dim5_moon.json').read_text())
    assert zonal['moon_mean_mm_day'] == pytest.approx(conv['moon_mean']['pr'], rel=2e-3)
    assert w.rain_share_meeting(clim['pr_mm_day'], clim['lat'], 0.) == 1.
    assert w.rain_share_meeting(clim['pr_mm_day'], clim['lat'], 1e3) == 0.


def test_rain_height_factor_is_one_at_the_ground():
    z = np.array([0., 5., 10.])
    out = w.rain_height_factor(dict(storm_height_km=z, storm_rain_water_g_kg=np.array([[1., 0.], [1., 1.], [0., 0.]]),
                                    air_density_kg_m3=np.array([1.4, 1.3, 1.2])))
    assert out['flux_over_surface'][0] == 1.
    assert out['flux_over_surface'][2] == 0.


def test_trim_rule_scales_with_gas_column():
    a = w.trim_pressure_per_kg(133.3, 1.204, 100108.)
    assert a == pytest.approx(100108/(1.204*133.3))
    assert w.allowed_swing(a, 133.3, 1.204, 100108.) == pytest.approx(1.)
    assert w.trim_pressure_per_kg(1333., 1.204, 100108.) == pytest.approx(a/10, rel=1e-3)


def test_sorbent_heat_balance_closes():
    row = w.hygroscopic_uptake(14.2, .64, .3, 100108., 5., 5.4, 1.204)
    assert row['latent_w_m2'] == pytest.approx((5+5.4)*row['warming_k'], rel=1e-6)
    assert row['uptake_kg_m2_day'] > 0
    assert w.hygroscopic_uptake(14.2, .64, .7, 100108.)['uptake_kg_m2_day'] == 0.


def test_vapour_work_and_dew_threshold_limits():
    assert w.vapour_work(.9, .9, 15.) == 0.
    assert w.vapour_work(.64, .99, 15.) == pytest.approx(8.314462618*288.15/.01801528*math.log(.99/.64), rel=1e-6)
    assert w.dew_threshold(14.2, 0.) == pytest.approx(1.)
    assert w.dew_threshold(14.2, 2.) < w.dew_threshold(14.2, 1.) < 1


def test_impaction_efficiency_limits():
    assert w.impaction_efficiency(0.) == 0.
    assert w.impaction_efficiency(1e6) == pytest.approx(1., rel=1e-5)
    st = w.stokes_number(10e-6, 1e-4, 1.)
    assert st == pytest.approx(2*1000*(5e-6)**2/(9*1.8e-5*5e-5))
    assert w.impaction_efficiency(math.pi/2) == pytest.approx(.5)
    assert w.droplet_capture(.1, 1., .3) == pytest.approx(.1e-3*.3*86400)


def test_evaluate_needs_and_closure_widen_with_thrift():
    out = w.evaluate()
    assert json.loads(json.dumps(out, allow_nan=False)) == out
    needs = {n['name']: n['water_kg_m2_day'] for n in out['needs']}
    assert needs['reference'] == pytest.approx(1.254, rel=2e-3)
    closure = {c['name']: c for c in out['rain_closure'] if c['catch_share'] == 1. and c['height_factor'] == 1.}
    assert closure['reference']['moon_share'] < closure['concentrating_cool_leaves']['moon_share'] < closure['concentrating_humid_layer']['moon_share']
    assert 0 < out['ring_rain_by_hour']['sector_fraction'] < .5


def test_closure_latitude_interpolates_to_the_need():
    zonal = dict(latitude_deg=[-45., -15., 15., 45.], rain_mm_day=[1., 3., 3., 0.])
    out = w.rain_closing_latitude(zonal, 1.5)
    assert out['north_deg'] == pytest.approx(15+30*(3-1.5)/3)
    assert out['south_deg'] == pytest.approx(-15-30*(3-1.5)/2)
    assert w.rain_closing_latitude(zonal, 4.) == dict(north_deg=None, south_deg=None)
    assert w.rain_closing_latitude(zonal, .5)['south_deg'] == -45.


def test_daylight_storage_counts_the_lit_hours_of_the_dry_gap():
    assert w.daylight_in_gap(5., 115.) == pytest.approx(95.)
    assert w.daylight_in_gap(300., 20.) == pytest.approx(100.)
    row = w.storage_daylight(1.25, 21., 5., 115.)
    assert row['daylight_dry_days'] == pytest.approx(21*95/360)
    assert row['storage_kg_m2'] == pytest.approx(2*1.25*21*95/360)
