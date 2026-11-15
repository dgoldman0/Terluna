"""The terrain beside the moisture map binds its inputs, its controls hold, its classes follow their rules on synthetic
ground, its shares close, and its per-pixel fields (when the research drive holds them) give its shares and carry each
pixel's moisture class from the cell that holds it on the turned Moon."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from biosphere.ecology import terrain_axis as ta
from biosphere.ecology.moisture_map import cell_index
from climate.gcm.turned import read_frame, to_turned

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(ta.OUT.read_text())
IGNORED = ('geography/inputs/', 'climate/gcm/products/', '/products/')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_product_binds_its_files_and_inputs():
    assert PRODUCT['schema'] == ta.SCHEMA
    for path, h in PRODUCT['producer']['files'].items():
        assert digest(ROOT / path)[:16] == h, path
    moisture = json.loads((ROOT / ta.MOISTURE).read_text())
    assert (PRODUCT['producer']['run'], PRODUCT['producer']['years']) == (moisture['producer']['run'],
                                                                         moisture['producer']['years'])
    for path, h in PRODUCT['producer']['inputs'].items():
        if (ROOT / path).exists():
            assert digest(ROOT / path) == h, path
        else:                                          # ignored inputs a fresh checkout may lack
            assert any(k in path for k in IGNORED), path
    assert PRODUCT['producer']['frame']['sha256'] == read_frame(ROOT / ta.ORIENTATION)['sha256']


def test_the_controls_hold():
    c = PRODUCT['checks']
    assert c['constant_terrain_mean_max_error_m'] < 1e-6 and c['constant_terrain_std_max_m'] < 1e-3
    assert c['constant_terrain_slope_max_deg'] < 1e-10
    assert c['analytic_slope_max_error_deg'] < 2e-5
    assert c['longitude_seam_max_error_m'] < 1e-6
    assert c['lakes_reproduced'] and c['lake_level_max_difference_m'] <= 0.01
    assert abs(c['global_area_closure']) < 1e-12
    for key in ('joint_partition_error', 'settings_partition_error', 'cue_partition_error'):
        assert abs(c[key]) < 1e-9, key
    assert c['land_classes_complete'] and c['water_unclassed'] and c['wetness_finite_on_land']


def test_the_settings_take_the_first_rule_that_applies():
    # one pixel per case: basin beats steep, steep beats crest, crest beats plain, hollow beats plain, plain, gentle
    pos = np.array([[0, 2, 2, 0, 1, 1, 1]], np.int8)
    slope = np.array([[15, 15, 1, 1, 1, 5, 5]], np.float32)
    depth = np.array([[20, 0, 0, 0, 0, 0, 9.9]])
    land = np.ones(pos.shape, bool)
    s = ta.settings(pos, slope, depth, land)
    assert s.tolist() == [[5, 1, 0, 4, 3, 2, 2]]
    land[0, 3] = False
    assert ta.settings(pos, slope, depth, land)[0, 3] == -1


def test_a_crater_has_a_low_floor_and_a_high_rim():
    lat = np.linspace(89.875, -89.875, 720)
    lon = np.linspace(0.125, 359.875, 1440)
    la, lo = np.meshgrid(lat[300:420], lon[600:840], indexing='ij')
    r_km = np.hypot((la - lat[360]) * 30.3, (lo - lon[720]) * 30.3 * np.cos(np.radians(la)))
    bowl = np.where(r_km < 20, -1500 * (1 - (r_km / 20) ** 2), 0) + 600 * np.exp(-((r_km - 22) / 4) ** 2)
    h = np.zeros((720, 1440))
    h[300:420, 600:840] = bowl
    land = np.ones(h.shape, bool)
    pos, rel, std = ta.positions(h, land, lat)                   # the whole sphere: the grid's spacing follows its rows
    centre = pos[360, 720]
    rim = pos[360, 720 + int(round(22 / (30.3 * 0.25 * math.cos(math.radians(lat[360])))))]
    assert centre == 0 and rim == 2 and pos[100, 100] == 1


def test_slope_of_a_tilted_plane_and_the_wetness_index():
    lat = np.linspace(89.875, -89.875, 720)
    lon = np.linspace(0.125, 359.875, 1440)
    plane = np.broadcast_to(1000.0 * np.sin(np.radians(lat))[:, None], (720, 1440)).copy()
    s = ta.slope_degrees(plane, lat, lon)
    expected = np.degrees(np.arctan(1000.0 / ta.MOON_RADIUS * np.cos(np.radians(lat))))
    assert np.max(np.abs(s[20:-20].astype(float) - expected[20:-20, None])) < 1e-4
    a = np.array([1e6, 1e8, 1e8, 1e8])
    sl = np.array([5.0, 5.0, 0.01, 0.0])
    w = ta.wetness(a, sl, 3.6e6)
    assert w[1] > w[0] and w[2] == pytest.approx(w[3])          # more drained area wetter; flat ground floored
    assert w[2] == pytest.approx(math.log(1e8 / math.sqrt(3.6e6) / math.tan(math.radians(ta.TWI_SLOPE_FLOOR_DEG))))


def test_the_shares_close():
    s = PRODUCT['shares']
    joint = np.array(s['moisture_by_position'])
    assert joint.sum() == pytest.approx(1.0, abs=2e-4)
    assert np.allclose(joint.sum(axis=1), s['moisture'], atol=2e-4)
    assert np.allclose(joint.sum(axis=0), s['positions'], atol=2e-4)
    by = np.array(s['settings_by_moisture'])
    assert np.allclose(by.sum(axis=0), s['settings'], atol=2e-4)
    assert np.allclose(by.sum(axis=1), s['moisture'], atol=2e-4)
    within = np.array(s['settings_within_moisture'])
    assert np.allclose(within.sum(axis=1)[np.array(s['moisture']) > 0], 1.0, atol=1e-3)
    assert sum(s['cues']) == pytest.approx(1.0, abs=2e-4)
    land = PRODUCT['land']
    assert 0 < land['spread_share'] < 0.5 and 0 < land['robust_share'] < 1


def test_the_lakes_by_class_are_the_drainage_products_lakes():
    path = ROOT / ta.DRAINAGE / 'products' / 'drainage_28pct_16ppd.npz'
    if not path.exists():
        pytest.skip('the forty-year drainage grid is an ignored product this checkout lacks')
    with np.load(path) as z:
        rows = z['lake_cell'].astype(np.int64) // 5760
    total = ta.row_areas_km2(2880, 5760)[rows].sum()
    lakes = PRODUCT['lakes']
    assert sum(lakes['area_km2_by_moisture']) == pytest.approx(total, abs=6)
    assert all(s is None or 0 <= s < 1 for s in lakes['share_of_ground_by_moisture'])


@pytest.fixture(scope='module')
def fields():
    path = Path(PRODUCT['fields']['path'])
    if not path.exists():
        pytest.skip('the per-pixel fields are on the research drive, which is not mounted')
    assert digest(path) == PRODUCT['fields']['sha256']
    with np.load(path) as z:
        return {k: z[k] for k in ('moisture', 'position', 'setting', 'spread', 'cue')}


def test_the_fields_give_the_shares(fields):
    m, pos, s, cue = fields['moisture'], fields['position'], fields['setting'], fields['cue']
    land = m >= 0
    assert np.array_equal(land, pos >= 0) and np.array_equal(land, s >= 0) and np.array_equal(land, cue >= 0)
    area = np.broadcast_to(ta.row_areas_km2(*m.shape)[:, None], m.shape)[land]
    total = area.sum()
    assert total == pytest.approx(PRODUCT['land']['area_km2'], abs=1)
    joint = np.bincount(m[land].astype(int) * 3 + pos[land], weights=area, minlength=18).reshape(6, 3) / total
    assert np.allclose(joint, PRODUCT['shares']['moisture_by_position'], atol=6e-5)
    by = np.bincount(m[land].astype(int) * 6 + s[land], weights=area, minlength=36).reshape(6, 6) / total
    assert np.allclose(by, PRODUCT['shares']['settings_by_moisture'], atol=6e-5)
    assert float(area[fields['spread'][land]].sum() / total) == pytest.approx(PRODUCT['land']['spread_share'], abs=6e-5)


def test_each_pixel_carries_the_class_of_the_cell_that_holds_it(fields):
    moisture = json.loads((ROOT / ta.MOISTURE).read_text())
    g, c = moisture['grid'], moisture['cells']
    cls = {(r, k): v for r, k, v in zip(c['row'], c['col'], c['moisture'])}
    rot = read_frame(ROOT / ta.ORIENTATION)['rotation']
    m, spread = fields['moisture'], fields['spread']
    rows, cols = np.nonzero((m >= 0) & ~spread)
    pick = np.random.default_rng(8).choice(rows.size, 2000, replace=False)
    lat = 90.0 - (rows[pick] + 0.5) / ta.PPD
    lon = (cols[pick] + 0.5) / ta.PPD
    tl, tn = to_turned(lat, lon, rot)
    r, k = cell_index(tl, tn, np.array(g['lat']), np.array(g['lon']))
    assert [cls[(int(a), int(b))] for a, b in zip(r, k)] == m[rows[pick], cols[pick]].tolist()
    rows_s, cols_s = np.nonzero(spread)
    tl, tn = to_turned(90.0 - (rows_s[:500] + 0.5) / ta.PPD, (cols_s[:500] + 0.5) / ta.PPD, rot)
    r, k = cell_index(tl, tn, np.array(g['lat']), np.array(g['lon']))
    assert all((int(a), int(b)) not in cls for a, b in zip(r, k))      # spread land lies in the run's sea cells
