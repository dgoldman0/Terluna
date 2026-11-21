"""Checks of the made sand seas' siting: the following of drainage to its ends, an outline's rule, the taper, the
haulage's arithmetic, what a margin sheds, the dunes' rule, and the written product."""
import hashlib
import json
import math

import numpy as np
import pytest

from research.studies.atlas_round import sand as S
from research.studies.atlas_round import sand_seas as M


def test_drainage_is_followed_to_its_first_end():
    # ten cells in a row draining eastward: cell 4 is an end, cell 9 its own receiver
    receiver = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 9], np.int32)
    stays = np.zeros(10, bool)
    stays[4] = True
    assert M.follow(receiver, stays).tolist() == [4, 4, 4, 4, 4, 9, 9, 9, 9, 9]
    long = np.minimum(np.arange(5000, dtype=np.int32) + 1, 4999)          # a way of five thousand cells ends at its last
    assert (M.follow(long, np.zeros(5000, bool)) == 4999).all()


def test_an_outline_takes_the_best_ground_round_its_seed_up_to_its_size():
    rows, cols = np.mgrid[0:120, 0:160]
    score = -np.hypot(rows - 60.0, cols - 70.0).astype(np.float32)        # best at the seed, falling away evenly
    allowed = np.ones(score.shape, bool)
    allowed[:, 90:] = False                                               # ground the sea may not cover, east of column 90
    area_row = np.full(120, 4.0)                                          # km2 a cell
    mask, t = M.outline((60, 70), allowed, score, area_row, 8000.0)
    assert mask[60, 70] and not mask[:, 90:].any()
    assert abs(mask.sum() * 4.0 - 8000.0) < 0.05 * 8000.0
    from scipy import ndimage
    assert ndimage.label(mask)[1] == 1 and not (ndimage.binary_fill_holes(mask) & ~mask).any()    # one piece, no holes
    assert mask[60, 88] and not mask[5, 70]                               # it runs up to the barred ground, within the rounding's cell, and stays round the seed
    small, _ = M.outline((60, 70), allowed, score, area_row, 400.0)
    assert small.sum() < mask.sum() and (mask | ~small).all()             # a smaller sea lies inside a larger one


def test_a_sea_thins_to_its_edge():
    mask = np.zeros((40, 40), bool)
    mask[5:35, 5:35] = True
    share = M.taper(mask)
    assert share[20, 20] == 1.0 and share[0, 0] == 0.0
    assert share[5, 20] == pytest.approx(M.TAPER_EDGE) and M.TAPER_EDGE < share[7, 20] < 1.0
    assert share[5 + int(M.TAPER_CELLS) + 1, 20] == 1.0
    blown = M.taper(mask, towards_deg=90.0)                               # a drift towards the east: a clean western edge, a feathered eastern one
    assert blown[20, 7] == 1.0 and blown[20, 34 - 3] < 0.6 and blown[20, 5] == blown[20, 34] == pytest.approx(M.TAPER_EDGE)
    assert np.allclose(blown[7, 12:28], blown[32, 12:28][::1])            # the two sides along the drift are alike
    narrow = M.taper(mask, aspect=0.5)                                    # columns half as wide as rows: the margin as wide on the ground
    assert narrow[20, 4 + 2 * 4] == pytest.approx(narrow[4 + 4, 20]) and narrow[4 + 4, 20] < 1.0   # eight columns in from the ground off the sea, as four rows


def test_haulage_is_tonnes_times_distance_at_the_loops_energy():
    # two ends on the equator, 100 km and 400 km east of the sea, a million cubic metres a year each; a sea's shore is not hauled from
    deg = lambda km: math.degrees(km * 1e3 / S.MOON_RADIUS)               # noqa: E731
    src = dict(lat=np.zeros(3), lon=np.array([deg(100.0), deg(400.0), deg(50.0)]), m3_yr=np.array([1e6, 1e6, 5e6]),
               kind=np.array([0, 1, 2]), basalt=np.array([1.0, 0.0, 1.0]))
    h = M.haul(src, 0.0, 0.0, need_km3=2.0, years=(1000.0, 2500.0), radii_km=(250.0, 500.0))
    near, both = h['by_radius']
    assert near['sand_km3_yr'] == pytest.approx(1e-3) and near['years'] == 2000 and near['mean_haul_km'] == 100 and near['basalt_share'] == 1.0
    watts = 1e6 * 1.6 * 100.0 * M.HAUL_KJ_PER_T_KM * 1e3 / S.SECONDS_PER_YEAR
    assert near['power_gw'] == pytest.approx(watts / 1e9, abs=5e-4)
    assert both['sand_km3_yr'] == pytest.approx(2e-3) and both['years'] == 1000 and both['mean_haul_km'] == 250 and both['basalt_share'] == 0.5
    assert both['energy_j'] == pytest.approx(2e9 * 1.6 * 250.0 * M.HAUL_KJ_PER_T_KM * 1e3, rel=1e-2)
    assert h['within_years']['2500'] == near and h['within_years']['1000'] == both
    assert M.haul(src, 0.0, 0.0, need_km3=2.0, years=(500.0,), radii_km=())['within_years']['500'] is None   # more than the ends give


def test_a_margin_sheds_downwind_what_the_drift_can_carry():
    shape = (S.ROWS, S.COLS)
    row, col = S.ROWS // 2, 2000
    mask = np.zeros(shape, bool)
    mask[row - 5:row + 5, col:col + 20] = True
    thick = np.where(mask, 10.0, 0.0).astype(np.float32)
    q = np.full(shape, 5.0, np.float32)                                   # m2 a year
    az = np.full(shape, 90.0, np.float32)                                 # eastward
    water = np.zeros(shape, bool)
    lat, area = S.fine_lat(), S.cell_area_km2()
    years = 1000.0
    thin, left, shed_km3, wet_km3 = M.margin_loss(mask, thick, q, az, water, lat, area, years)
    width = math.sqrt(area[row] * 1e6)
    assert shed_km3 == pytest.approx(10 * 5.0 * width * years / 1e9, rel=0.02) and wet_km3 == 0.0   # the ten cells of its eastern edge
    assert (thin[row - 5:row + 5, col + 20] > 0).all() and thin[:, :col + 20].max() == 0.0          # set down just east of it
    whole = lambda a: float((a * area[:, None]).sum())                    # noqa: E731
    assert whole(thin) + whole(left) == pytest.approx(whole(thick), rel=1e-5)                       # the sand is kept
    assert left[row, col + 19] < 10.0 and left[row, col + 5] == 10.0
    water[row - 5:row + 5, col + 20] = True                               # water at its eastern edge takes what it sheds
    thin, _left, shed_km3, wet_km3 = M.margin_loss(mask, thick, q, az, water, lat, area, years)
    assert wet_km3 == pytest.approx(shed_km3) and thin.max() == 0.0
    huge = M.margin_loss(mask, thick, q * 1e6, az, np.zeros(shape, bool), lat, area, years)         # a drift past all its sand takes no more than the cell holds
    assert huge[1].min() >= 0.0 and huge[2] == pytest.approx(10 * 10.0 * area[row] * 1e6 / 1e9, rel=0.02)


def test_dunes_follow_the_sand_studys_rules():
    star = M.dunes_of(0.27, 21.0, 0.42, 25.0)
    assert star['form'] == 'star dunes' and star['crest_spacing_m'] == 1500 and star['height_m'] == 75 and star['limited_by'] == 'the sand at hand'
    assert M.dunes_of(0.89, 29.0, 0.48, 30.0)['form'] == 'crescent dunes and transverse ridges'
    assert M.dunes_of(0.5, 20.0, 0.4, 30.0)['form'] == 'linear dunes'
    young = M.dunes_of(0.9, 20.0, 0.5, 30.0, years=100.0)                 # a hundred years: the drift has gathered little
    assert young['limited_by'] == 'the drift and the years' and young['crest_spacing_m'] == round(math.sqrt(2 * S.CREST_ASPECT * 10.0 * 100.0))


def test_product_and_record():
    if not (M.PRODUCT.exists() and M.OUT.exists()):
        pytest.skip('the product is on the research drive')
    record = json.loads(M.OUT.read_text())
    assert record['schema'] == M.SCHEMA
    for code in M.CODE:
        if record['producer']['files'][code] != hashlib.sha256((M.ROOT / code).read_bytes()).hexdigest():
            pytest.skip(f'{code} has changed since the record was written')
    assert record['product']['sha256'] == S.sha256(M.PRODUCT)
    z = np.load(M.PRODUCT, allow_pickle=False)
    meta = json.loads(str(z['metadata']))
    assert meta['schema'] == M.SCHEMA and [s['name'] for s in meta['seas']] == list(record['seas'])
    made, share, thick, part = z['made_sea'], z['sand_share'], z['sand_thickness_m'], z['made_sand_part']
    area = S.cell_area_km2()[:, None]
    supply = record['supply']
    assert sum(supply['km3_yr'].values()) == pytest.approx(supply['land_sheds_km3_yr'], rel=2e-3) and supply['unreached_km3_yr'] == 0.0
    for sea in meta['seas']:
        g = record['seas'][sea['name']]
        inside = made == sea['id']
        assert float((area * inside).sum()) == pytest.approx(g['area_km2'], abs=1.0)
        assert float((share * area)[inside].sum() / (area * inside).sum()) > 0.97 and part[inside].min() == 1.0   # the sea is sand, thin only on its swells,
        #                                                                                                        and the sand is its own
        assert float((thick * area)[inside].sum() / 1e3) == pytest.approx(g['sand_km3'], rel=0.12)   # its sand lies inside it, but for the fringe
        plan = g['plan']
        assert plan['years'] <= plan['deliver_years'] and 0.0 <= plan['basalt_share'] <= 1.0
        assert plan['power_gw'] < M.LOOP_RANGE_GW[1]
    assert int(made.max()) == len(meta['seas']) and float(part[share == 0].max()) == 0.0
    gain = record['map']['seas']['sand_km2'] - record['map']['open']['sand_km2']
    assert gain == pytest.approx(sum(g['area_km2'] for g in record['seas'].values()), rel=0.15)   # the seas and their fringes
