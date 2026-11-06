import json
import math

import numpy as np
import pytest

from research.studies.seawater_lakes import sites as S


def test_the_pass_is_the_lowest_notch_in_the_ridge():
    h = np.zeros((5, 7))
    h[:, 3] = [30.0, 12.0, 25.0, 8.0, 40.0]                 # a ridge with notches of 12 and 8
    start = np.zeros(h.shape, bool)
    start[:, 0] = True
    goal = np.zeros(h.shape, bool)
    goal[:, 6] = True
    assert S.minimax_pass(h, start, goal) == pytest.approx(8.0)


def test_lake_geometry_below_a_level():
    heights = np.array([-10.0, -4.0, 2.0])
    areas = np.array([1.0, 2.0, 3.0])
    g = S.lake_geometry(heights, areas, 0.0)
    assert g['area_km2'] == pytest.approx(3.0)
    assert g['volume_km3'] == pytest.approx((10.0 * 1.0 + 4.0 * 2.0) * 1e-3)
    assert g['greatest_depth_m'] == pytest.approx(10.0)
    assert g['mean_depth_m'] == pytest.approx(6.0)


def test_larval_spread_is_sqrt_2kt():
    assert S.larval_spread_km(1.0, 120.0) == pytest.approx(math.sqrt(2.0 * 120.0 * 86400.0) / 1e3)


@pytest.mark.skipif(not S.FRAME.exists(), reason='the frame record sits in the domain/climate worktree')
def test_the_turn_reproduces_the_frame_records_places():
    rec = json.loads(S.FRAME.read_text())
    rot = rec['frame']['rotation_today_to_turned']
    for p in rec['places']:
        lat, lon = S.to_turned(*p['today'], rot)
        assert float(lat) == pytest.approx(p['turned'][0], abs=2e-3)
        assert (float(lon) - p['turned'][1] + 180.0) % 360.0 - 180.0 == pytest.approx(0.0, abs=2e-3)


@pytest.mark.skipif(not (S.CACHE.exists() and S.LAND_READING.exists()), reason='the run cache sits on the research drive')
def test_the_reader_matches_the_turned_years_land_reading():
    """Rain, demand and aridity at the places turned_years_land.py read on the same run."""
    c, _ = S.read_run()
    ref = json.loads(S.LAND_READING.read_text())['places']
    for name, p in ref.items():
        row = p['turned_years']
        j, i, _ = S.nearest_cell(c, row['lat'], row['lon_e'], row['land'])
        assert c['lat'][j] == pytest.approx(row['lat']) and c['lon'][i] == pytest.approx(row['lon_e'])
        assert round(float(c['rain'][j, i] * 365.0)) == pytest.approx(row['rain_mm_yr'], abs=1), name
        assert round(float(c['pet'][j, i] * 365.0)) == pytest.approx(row['demand_mm_yr'], abs=1), name
        assert float(c['ai'][j, i]) == pytest.approx(row['aridity'], abs=1.5e-3), name


@pytest.mark.skipif(not (S.CACHE.exists() and S.FRAME.exists() and S.ATLAS_GRID.exists()),
                    reason='the screen reads the research drive, the climate worktree and the atlas grid')
def test_the_sites_product_regenerates():
    assert json.loads(S.OUT.read_text()) == json.loads(json.dumps(S.build()))
