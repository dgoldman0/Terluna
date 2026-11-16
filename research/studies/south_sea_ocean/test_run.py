import json

import numpy as np
import pytest

from research.studies.sea_chemistry import run as sc
from research.studies.south_sea_ocean import run as R

NEEDS = [R.ROOT / R.ATLAS_GRID, R.ROOT / R.CLIMATOLOGY]
needs_products = pytest.mark.skipif(not all(p.exists() for p in NEEDS),
                                    reason='the atlas grid or the design years\' climatology is not on this machine')


def test_the_ions_sum_to_the_reference_salinity():
    grams = R.ion_grams_per_kg()
    assert sum(grams.values()) == pytest.approx(sc.SEAWATER_G_KG)
    assert 0.0 < grams['minor'] < 0.2                      # bromide, strontium, boron and fluoride
    assert grams['cl'] / sc.SEAWATER_G_KG == pytest.approx(0.55, abs=0.005)
    assert grams['na'] / sc.SEAWATER_G_KG == pytest.approx(0.307, abs=0.005)


def test_seawaters_charges_balance_within_its_minor_ions():
    s = sc.SEAWATER
    cations = s['na'] + s['k'] + 2 * s['mg'] + 2 * s['ca']
    anions = s['cl'] + 2 * s['so4'] + s['alkalinity']
    assert cations == pytest.approx(anions, rel=2e-3)


def test_the_salt_agrees_with_the_sea_chemistry_and_the_seawater_lakes():
    """The sea chemistry gives 174,000 Gt for this sea and the seawater lakes 35.94 Mt of salt per km3."""
    ocean = R.salt(4.84e6, 1.0)
    assert ocean['salt_kg'] == pytest.approx(1.74e17, rel=5e-3)
    assert ocean['salt_kg'] / 4.84e6 == pytest.approx(35.94e9, rel=1e-3)
    assert sum(ocean['ions_kg'].values()) == pytest.approx(ocean['salt_kg'])


def test_a_weaker_sea_holds_its_share_of_the_salt_at_its_own_density():
    ocean, half = R.salt(1000.0, 1.0), R.salt(1000.0, 0.5)
    assert half['salinity_g_kg'] == pytest.approx(0.5 * ocean['salinity_g_kg'])
    assert half['density_kg_m3'] < ocean['density_kg_m3']
    assert half['ions_kg']['cl'] / ocean['ions_kg']['cl'] == pytest.approx(0.5 * half['density_kg_m3']
                                                                             / ocean['density_kg_m3'])


def grids():
    lat = 90.0 - (np.arange(720) + 0.5) / 4.0
    lon = (np.arange(1440) + 0.5) / 4.0 - 180.0
    glat = np.linspace(85.0, -85.0, 32)
    glon = np.arange(64) * 5.625
    return lat, lon, glat, glon


def test_the_sea_keeps_its_area_on_the_models_grid_and_lands_in_the_nearest_cells():
    lat, lon, glat, glon = grids()
    mask = np.zeros((720, 1440), bool)
    rows = np.flatnonzero((lat < -40.0) & (lat > -50.0))
    cols = np.flatnonzero((lon > 10.0) & (lon < 30.0))
    mask[np.ix_(rows, cols)] = True
    area = R.on_model_grid(mask, lat, lon, np.eye(3), glat, glon)
    expected = (R.cell_areas_km2(lat, lon)[rows] * cols.size).sum()
    assert area.sum() == pytest.approx(expected)
    filled_rows = np.flatnonzero(area.sum(axis=1) > 0)
    assert glat[filled_rows].max() < -35.0 and glat[filled_rows].min() > -55.0
    filled_cols = np.flatnonzero(area.sum(axis=0) > 0)
    assert glon[filled_cols].min() >= 5.0 and glon[filled_cols].max() <= 35.0


def test_a_quarter_turn_moves_the_sea_with_the_frame():
    lat, lon, glat, glon = grids()
    mask = np.zeros((720, 1440), bool)
    mask[np.ix_(np.flatnonzero(np.abs(lat) < 5.0), np.flatnonzero(np.abs(lon) < 5.0))] = True
    quarter = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])     # 90 degrees east about the pole
    area = R.on_model_grid(mask, lat, lon, quarter, glat, glon)
    filled_cols = np.flatnonzero(area.sum(axis=0) > 0)
    assert np.all(np.abs(glon[filled_cols] - 90.0) <= 10.0)


def test_the_water_balance_reads_only_the_models_water_cells():
    glat, glon = np.linspace(85.0, -85.0, 4), np.arange(4) * 90.0
    area = np.zeros((4, 4))
    area[1, :2] = (300.0, 100.0)
    clim = dict(lat=glat, lon=glon, lsm=np.zeros((4, 4)), pr_mm_day=np.full((4, 4), 1.0),
                evap_mm_day=np.full((4, 4), 3.0), ts_k=np.full((4, 4), 300.0))
    clim['lsm'][1, 1] = 1.0                                  # a land cell of the model under part of the sea
    clim['evap_mm_day'][1, 1] = 0.2
    b = R.water_balance(area, clim, 400.0)
    assert b['share_of_area_in_model_water_cells'] == pytest.approx(0.75)
    assert b['evaporation_mm_day'] == pytest.approx(3.0)
    assert b['net_evaporation_m_yr'] == pytest.approx(2.0 * 365.25 / 1000.0)
    assert b['net_evaporation_km3_yr'] == pytest.approx(2.0 * 365.25 / 1000.0 * 400.0 * 1e-3)
    low, high = b['net_evaporation_m_yr_over_the_rain_range']
    assert low == pytest.approx((3.0 - 2.0 * 1.0) * 0.36525) and high == pytest.approx((3.0 - 0.75) * 0.36525)


def test_rivers_that_match_the_loss_leave_the_level_alone():
    r = R.rivers_against_need(700.0, 700.0, 0.7, 1.0e6)
    assert r['surplus_km3_yr'] == 0.0 and r['level_change_m_yr'] == 0.0
    assert r['balanced_area_share_of_sea'] == pytest.approx(1.0)
    short = R.rivers_against_need(350.0, 700.0, 0.7, 1.0e6)
    assert short['share_of_need'] == pytest.approx(0.5) and short['level_change_m_yr'] == pytest.approx(-0.35)
    assert short['surplus_m3s'] == pytest.approx(-350.0e9 / R.YEAR_S)


def test_density_slopes_have_seawaters_size_and_sign():
    slopes = R.density_slopes(35.165, 29.5)
    assert 0.30 < slopes['kg_m3_per_k_of_cooling'] < 0.36      # thermal expansion near 3.3e-4 per K at 30 C
    assert 0.72 < slopes['kg_m3_per_g_kg'] < 0.78
    rate = R.salting(0.7, 35.0, slopes)['100_m']
    assert rate['g_kg_a_year'] == pytest.approx(35.0 * 0.7 / 100.0)
    assert rate['equal_to_cooling_k_a_year'] == pytest.approx(rate['kg_m3_a_year'] / slopes['kg_m3_per_k_of_cooling'])


def test_the_floors_take_grows_with_time_and_carries_its_salt():
    uptake = json.loads((R.ROOT / R.UPTAKE).read_text())
    ocean = R.salt(1000.0, 1.0)
    bands = {'45-60': 1000.0, '60-75': 500.0}
    out = R.floor_uptake(bands, uptake, ocean, 1000.0)
    for case in ('slow', 'central', 'fast'):
        water = [out[case][y]['water_km3'] for y in R.UPTAKE_YEARS]
        assert water == sorted(water) and water[0] > 0.0
        for y in R.UPTAKE_YEARS:
            row = out[case][y]
            assert row['salt_kg'] == pytest.approx(row['water_km3'] * 1e9 * ocean['density_kg_m3'] * 35.165e-3)
    assert out['slow']['10000']['water_km3'] < out['central']['10000']['water_km3'] < out['fast']['10000']['water_km3']
    held = uptake['cases']['central']['bands']['45-60']['columns']['sea_floor']['by_year']['100']['held_m']
    held2 = uptake['cases']['central']['bands']['60-75']['columns']['sea_floor']['by_year']['100']['held_m']
    assert out['central']['100']['water_km3'] == pytest.approx((1000.0 * held + 500.0 * held2) * 1e-3)


@needs_products
def test_the_product_regenerates():
    assert json.loads(R.OUT.read_text()) == json.loads(json.dumps(R.build()))


def test_the_products_headline_relations_hold():
    d = json.loads(R.OUT.read_text())
    sea, own = d['sea'], d['water'][R.SEA_NAME]
    assert sum(d['floor']['sea_area_by_band_km2'].values()) == pytest.approx(sea['area_km2'], rel=5e-3)
    assert own['share_of_area_in_model_water_cells'] > 0.9
    assert own['net_evaporation_km3_yr'] == pytest.approx(own['net_evaporation_m_yr'] * sea['area_km2'] * 1e-3)
    assert d['held']['fresh_inflow_km3_yr'] == own['net_evaporation_km3_yr']
    assert d['held']['years_to_evaporate_its_volume'] == pytest.approx(sea['volume_km3'] / own['net_evaporation_km3_yr'])
    shares = [d['rivers'][rule][R.SEA_NAME]['share_of_need'] for rule in ('every_water', 'six_seas', 'two_seas')]
    assert shares == sorted(shares)                            # each rule lets more of the land's water through
    ocean = d['salt']['strengths']['ocean']
    assert ocean['ions_kg']['cl'] / ocean['salt_kg'] == pytest.approx(0.55, abs=0.005)
    for hollow in d['filling']['hollows']:
        assert hollow['years_to_fill'] == pytest.approx(hollow['volume_above_the_common_level_km3']
                                                        / hollow['inflow_km3_yr'])
    for name, s in d['filling']['smaller_seas'].items():
        assert (s['years_to_fill'] is None) == (s['surplus_km3_yr'] <= 0), name
