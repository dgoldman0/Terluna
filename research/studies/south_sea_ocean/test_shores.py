import json
import os

import numpy as np
import pytest

from research.studies.south_sea_ocean import run as R
from research.studies.south_sea_ocean import shores as SH


def test_a_pixel_faces_the_sea_beside_it_and_across_the_seam():
    sea = np.zeros((4, 8), bool)
    sea[1, 0] = True                           # touches (1, 7) across the seam, (0, 0), (2, 0) and (1, 1)
    sea[3, 4] = True                           # on the last row: nothing lies beyond the pole
    cols = 8
    index = np.array([1 * cols + 7, 0 * cols + 0, 2 * cols + 0, 1 * cols + 1, 2 * cols + 7, 3 * cols + 3,
                      2 * cols + 4, 0 * cols + 4])
    assert SH.facing(index, sea, cols).tolist() == [True, True, True, True, False, True, True, False]


def test_rain_classes_break_at_half_one_and_two_millimetres():
    assert SH.rain_class([0.0, 0.49, 0.5, 0.99, 1.0, 1.99, 2.0, 9.0]).tolist() == [0, 0, 1, 1, 2, 2, 3, 3]
    assert [c[0] for c in SH.RAIN_CLASSES] == ['desert', 'half_dry', 'moist', 'wet']


def test_the_rain_is_read_over_a_place_and_at_the_nearest_land():
    lat = np.array([60.0, 0.0, -60.0])
    lon = np.arange(0.0, 360.0, 90.0)
    pr = np.arange(12.0).reshape(3, 4)
    lsm = np.zeros((3, 4))
    lsm[0, 1] = 1.0                            # the only land cell, at 60 N, 90 E
    over, land, is_land = SH.rain_at(np.array([-50.0, 55.0]), np.array([265.0, 100.0]),
                                     dict(lat=lat, lon=lon, pr_mm_day=pr, lsm=lsm))
    assert over.tolist() == [pr[2, 3], pr[0, 1]]
    assert land.tolist() == [pr[0, 1], pr[0, 1]]
    assert is_land.tolist() == [False, True]


def test_the_high_water_spreads_its_amplitude_over_the_slope_and_stops_at_the_pixel():
    # a pixel 2 km deep that stands 20 m above the sea at half its depth: a slope of 1 in 50
    assert SH.high_water_strip(1.6, 20.0, 2.0) == pytest.approx(80.0)
    assert SH.high_water_strip(1.6, 0.0, 2.0) == pytest.approx(2000.0)         # flat ground floods to its edge
    assert SH.high_water_strip(100.0, 20.0, 2.0) == pytest.approx(2000.0)
    assert SH.high_water_strip(1.6, np.nan, 2.0) == pytest.approx(2000.0)


def test_a_belt_ends_where_a_steep_shore_outruns_the_reach_and_never_exceeds_its_pixel():
    coast = dict(edge=np.array([2.0, 2.0]), across=np.array([2.0, 2.0]), rise=np.array([500.0, 0.5]),
                 area=np.array([4.0, 4.0]))
    inland = dict(distance=np.array([3.0]), step_rows=np.array([2]), step_cols=np.array([0]),     # two pixels inland
                  height=np.array([900.0]), area=np.array([4.0]))
    uptake, t, reach = 1.94e-3 / 86400.0, 5.2e-4, 15.0
    area, width, length, far = SH.belts_for(coast, inland, uptake, t, reach)
    assert width[0] == pytest.approx(reach / 0.5, rel=0.01)                    # a slope of 1 in 2: W = R / s
    assert width[1] == pytest.approx(np.sqrt(2.0 * t * reach / uptake), rel=0.05)   # nearly flat: the ground's supply
    assert length == pytest.approx(np.pi / 4.0 * coast['edge'])
    assert area[0] == pytest.approx(length[0] * width[0] / 1000.0)
    assert (area <= coast['area']).all() and far[0] == 0.0


def test_percentiles_weigh_by_length():
    assert SH.percentiles([1.0, 2.0, 3.0], [1.0, 1.0, 8.0]) == [1.0, 3.0, 3.0]
    assert SH.percentiles([1.0], [0.0]) == [None, None, None]


def test_the_product_is_current_and_its_parts_add_up():
    out = json.loads(SH.OUT.read_text())
    assert out['schema'] == SH.SCHEMA
    for path, digest in out['producer']['files'].items():
        assert SH.sha256(SH.ROOT / path) == digest, path
    for path in (SH.I.ATLAS, SH.COAST_WATER, SH.I.FRAME):
        assert SH.sha256(SH.ROOT / path) == out['inputs'][path], path
    s = out['shore']
    assert s['mainland_km'] + s['island_km'] == pytest.approx(s['length_km'])
    assert sum(s['by_rain_km'].values()) == pytest.approx(s['length_km'])
    assert sum(s['by_rain_at_the_cell_over_it_km'].values()) == pytest.approx(s['length_km'])
    assert sum(b['all'] for b in s['by_turned_latitude_km'].values()) == pytest.approx(s['length_km'])
    assert sum(b['desert'] for b in s['by_turned_latitude_km'].values()) == pytest.approx(s['by_rain_km']['desert'])
    assert s['pixels_beside_two_seas'] == 0
    assert s['length_km'] < s['moon_km']['seas'] < s['moon_km']['every_water']
    assert out['sea']['area_km2'] == pytest.approx(out['sea']['atlas_area_km2'], rel=0.02)
    # every shore poleward of 60 S on the turned Moon lies in desert
    for band in ('60-75S', '75-90S'):
        assert s['by_turned_latitude_km'][band]['desert'] == pytest.approx(s['by_turned_latitude_km'][band]['all'])


def test_the_belts_draw_sea_water_at_seawaters_strength():
    out = json.loads(SH.OUT.read_text())
    b = out['fresh_sea_belts']
    assert b['desert_shore_km'] == pytest.approx(out['shore']['by_rain_km']['desert'])
    assert b['salt_kg_per_m3_of_sea_water'] == pytest.approx(R.salt(1.0, 1.0)['salt_kg'] / 1e9)
    assert b['salt_kg_per_m3_of_sea_water'] == pytest.approx(35.94, abs=0.01)
    for case in b['cases'].values():
        assert case['belt_km2'] == pytest.approx(case['beside_the_water_km2'] + case['further_inland_km2'])
        assert case['salt_in_that_draw_mt_yr'] == pytest.approx(case['draw_km3_yr'] * b['salt_kg_per_m3_of_sea_water'])
        assert case['draw_km3_yr'] * 1e9 == pytest.approx(
            case['belt_km2'] * 1e6 * b['uptake_mm_day'] / 1000.0 * 365.25)
        assert 100.0 < case['belt_per_km_of_shore_m'] < 500.0
    cases = b['cases']
    assert cases['slow_15m']['belt_km2'] < cases['central_15m']['belt_km2'] < cases['sandy_15m']['belt_km2']
    assert cases['central_15m']['belt_km2'] < cases['central_25m']['belt_km2']
    t = out['tide']
    assert sum(t['by_rain_km2'].values()) == pytest.approx(t['covered_at_high_water_km2'])
    assert 1.0 < t['typical_monthly_range_m_p10_p50_p90'][1] < 6.0


@pytest.mark.skipif(os.environ.get('TERLUNA_SOUTH_SEA_REGENERATE') != '1',
                    reason='the 16 pixel/degree shores, a minute and 2 GB; set TERLUNA_SOUTH_SEA_REGENERATE=1')
def test_the_shores_regenerate():
    assert json.loads(SH.OUT.read_text()) == json.loads(json.dumps(SH.build()))
