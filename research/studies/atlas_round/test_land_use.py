"""Tests of the first sorting of the design surface's land into wild, hybrid and open land (land_use.py)."""

import json

import numpy as np
import pytest

from research.studies.atlas_round import land_use as LU


def test_the_climate_cell_is_found_by_turned_latitude_and_longitude():
    edges = [90.0, 60.0, 0.0, -60.0, -90.0]
    rows, cols = LU.t21_cell([75.0, 30.0, -30.0, -89.0], [0.0, 22.4, 22.6, 359.0], edges, 8)
    assert rows.tolist() == [0, 1, 2, 3]
    assert cols.tolist() == [0, 0, 1, 0]                     # column 0 is centred on longitude 0 and 45 degrees wide


def test_turning_keeps_a_place_under_no_rotation_and_moves_the_pole_under_one():
    lat, lon = LU.turned([10.0, -45.0], [20.0, 300.0], np.eye(3))
    assert np.allclose(lat, [10.0, -45.0]) and np.allclose(lon, [20.0, 300.0])
    tilt = np.radians(7.0)
    about_y = np.array([[np.cos(tilt), 0, -np.sin(tilt)], [0, 1, 0], [np.sin(tilt), 0, np.cos(tilt)]])
    lat, _ = LU.turned([90.0], [0.0], about_y)
    assert lat[0] == pytest.approx(83.0, abs=1e-6)


def test_distances_follow_the_great_circle_within_a_few_percent():
    rows, cols = 360, 720
    lat, lon, cell_km, _ = LU.grid(rows, cols)
    source = np.zeros((rows, cols), bool)
    r0, c0 = 99, 200                                         # 40.25 N, 100.25 E
    source[r0, c0] = True
    d = LU.distance_km(source, lat, cell_km, band_rows=10, reach_km=300.0)
    rr, cc = np.meshgrid(np.arange(r0 - 12, r0 + 13), np.arange(c0 - 12, c0 + 13), indexing='ij')
    true = LU.great_circle_km(lat[rr], lon[cc], lat[r0], lon[c0])
    near = true < 250.0
    assert np.all(np.abs(d[rr, cc][near] - true[near]) <= 0.08 * true[near] + 0.5 * cell_km)
    assert d[0, 0] == 300.0                                  # capped where nothing lies within reach


def test_distances_wrap_round_the_seam():
    rows, cols = 360, 720
    lat, lon, cell_km, _ = LU.grid(rows, cols)
    source = np.zeros((rows, cols), bool)
    source[180, 0] = True
    d = LU.distance_km(source, lat, cell_km, band_rows=10)
    assert d[180, cols - 1] == pytest.approx(d[180, 1], rel=1e-6) and d[180, cols - 1] < 20.0


def test_a_lake_across_the_seam_is_one_lake():
    lake = np.zeros((5, 10), bool)
    lake[2, 0:2] = lake[2, 9] = True                         # three cells joined across the seam
    lake[0, 5] = True                                        # one cell alone
    big = LU.great_lakes(lake, np.ones(5), min_km2=3.0)
    assert big[2, 0] and big[2, 1] and big[2, 9] and not big[0, 5]


def test_outlets_follow_the_receivers_to_the_sea():
    receiver = np.array([1, 2, 3, 3, 4])                     # 0 -> 1 -> 2 -> 3, a sea cell; 4 drains nowhere
    sea = np.array([-1, -1, -1, 0, -1], np.int8)
    assert LU.outlets([0, 2, 4], receiver, sea).tolist() == [3, 3, 4]


def test_ranks_follow_height_latitude_and_moisture():
    moisture = np.array([1, 1, 5, 4, 1, 2, 5], np.int8)
    turned_lat = np.array([50.0, -60.0, 10.0, 40.0, 10.0, 20.0, 40.0], np.float32)
    up = np.array([5000.0, 100.0, 100.0, 100.0, 5000.0, 100.0, 6000.0], np.float32)
    assert LU.habitation_rank(moisture, turned_lat, up).tolist() == [1, 2, 4, 3, 5, 6, 6]


def land(n=1200, seed=3):
    """A made land of equal cells: half dry and half wet, a quarter of the wet cells candidate valleys."""
    rng = np.random.default_rng(seed)
    moisture = np.where(np.arange(n) % 2 == 0, 2, 5).astype(np.int8)
    candidate = (moisture == 5) & (np.arange(n) % 8 == 1)
    shore = rng.uniform(0, 100, n).astype(np.float32)
    place = rng.permutation(n)
    return np.ones(n), moisture, candidate, shore, place


def test_the_sorting_keeps_the_wild_share_in_every_group_and_the_hybrid_share_of_the_land():
    area, moisture, candidate, shore, place = land()
    fixed = np.zeros(area.size, np.int8)
    kind, reach = LU.sort_land(area, moisture, moisture.astype(np.int64), fixed, candidate, shore, place, 0.30, 0.20)
    for c in (2, 5):
        assert np.mean(kind[moisture == c] == LU.WILD) == pytest.approx(0.30, abs=0.002)
    assert np.mean(kind == LU.HYBRID) == pytest.approx(0.20, abs=0.002)
    assert (kind[candidate] == LU.HYBRID).all()              # the candidate valleys first
    gardens = (kind == LU.HYBRID) & ~candidate
    assert (moisture[gardens] == 2).all()                    # gardens only in the dry class
    assert shore[gardens].max() == pytest.approx(reach)
    assert shore[gardens].max() <= shore[(moisture == 2) & ~gardens].min()   # and the nearest water first
    wet_free = (moisture == 5) & ~candidate
    assert place[wet_free & (kind == LU.WILD)].min() > place[wet_free & (kind == LU.OPEN)].max()   # the wild takes the last-ranked


def test_fixed_places_keep_their_kind_and_count_towards_the_shares():
    area, moisture, candidate, shore, place = land()
    fixed = np.zeros(area.size, np.int8)
    fixed[:10] = LU.WILD
    fixed[10:20] = LU.OPEN
    kind, _ = LU.sort_land(area, moisture, moisture.astype(np.int64), fixed, candidate, shore, place, 0.50, 0.20)
    assert (kind[:10] == LU.WILD).all() and (kind[10:20] == LU.OPEN).all()
    assert np.mean(kind == LU.WILD) == pytest.approx(0.50, abs=0.003)


def test_the_mosaics_go_where_the_variant_puts_them():
    area, moisture, candidate, shore, place = land()
    fixed = np.zeros(area.size, np.int8)
    strata = moisture.astype(np.int64)
    kind, _ = LU.sort_land(area, moisture, strata, fixed, candidate, shore, place, 0.30, 0.20, mosaics='habitation')
    assert (kind[candidate] == LU.OPEN).all() and np.mean(kind == LU.HYBRID) == pytest.approx(0.20, abs=0.002)
    kind, _ = LU.sort_land(area, moisture, strata, fixed, candidate, shore, place, 0.30, 0.20, mosaics='beside')
    assert np.mean(kind == LU.HYBRID) == pytest.approx(0.20 + candidate.mean(), abs=0.002)
    with pytest.raises(ValueError):
        LU.sort_land(area, moisture, strata, fixed, candidate, shore, place, 0.30, 0.20, mosaics='elsewhere')


def test_people_take_the_best_open_land_of_their_band_and_a_shortfall_is_reported():
    area = np.ones(10)
    kind = np.array([3, 3, 3, 3, 1, 2, 3, 3, 3, 3], np.int8)
    band = np.array([1, 1, 1, 1, 1, 1, 4, 4, 0, 0], np.int8)
    place = np.array([3, 0, 2, 1, 4, 5, 6, 7, 8, 9])
    used, table = LU.place_people(area, kind, band, place, {1: 2.0, 4: 5.0})
    assert used.tolist() == [False, True, False, True, False, False, True, True, False, False]
    assert table[1] == dict(open_km2=4, need_km2=2, placed_km2=2, short_km2=0)
    assert table[4] == dict(open_km2=2, need_km2=5, placed_km2=2, short_km2=3)


def test_the_needs_by_band_are_the_studies_figures():
    low, high = LU.band_needs(0.30, 0), LU.band_needs(0.30, 1)
    assert low[4] == pytest.approx(680000 + 0.475 * 100000) and high[4] == pytest.approx(1860000 + 0.475 * 571000)
    assert LU.band_needs(0.50, 1)[4] == pytest.approx(2170000 + 0.475 * 571000)   # the fields under half wild
    assert high[1] == 16900 + 193000 and high[3] == 595000
    assert sum(LU.band_needs(0.30, 1).values()) + 2464 < 0.15 * 22.13e6


@pytest.mark.skipif(not LU.OUT.exists(), reason='the sorting has not been run')
def test_record_is_consistent_and_bound_to_its_code_and_inputs():
    r = json.loads(LU.OUT.read_text())
    assert r['schema'] == LU.SCHEMA
    assert r['producer']['files'] == {f: LU.sha256(LU.ROOT / f) for f in LU.CODE}       # made by this code
    for path, digest in r['inputs'].items():
        assert LU.sha256(LU.ROOT / path) == digest                                      # from these inputs
    inv = r['inventory']
    for k, share in inv['moisture_map_shares_over_the_design_land'].items():
        assert inv['by_moisture_class'][k]['share'] == pytest.approx(share, abs=0.005)  # the moisture map's own count
    assert set(r['variants']) == {v['key'] for v in LU.VARIANTS}
    for v in r['variants'].values():
        m = v['moon']
        assert m['wild'] + m['hybrid'] + m['open'] == pytest.approx(1.0, abs=1e-3)
        assert m['wild'] == pytest.approx(v['wild_share_set'], abs=0.002)
        assert m['filled_at_the_high_end'] <= m['open']
        for group in v['by'].values():
            for row in group.values():
                assert row['wild'] + row['hybrid'] + row['open'] == pytest.approx(1.0, abs=2e-3)
                assert row['filled'] <= row['open'] + 1e-6
        for row in v['by']['moisture_class'].values():
            assert row['wild'] == pytest.approx(v['wild_share_set'], abs=0.003)
        assert v['by']['night_regime']['dark']['wild'] == pytest.approx(v['wild_share_set'], abs=0.003)
        for end in ('low', 'high'):
            for row in v['room']['by_band'][end].values():
                assert row['placed_km2'] <= row['open_km2']
                assert abs(row['short_km2'] - max(row['need_km2'] - row['open_km2'], 0)) <= 1
    assert r['variants']['A']['moon']['hybrid'] == pytest.approx(0.20, abs=0.001)
    assert r['variants']['C']['hybrid']['candidate_valleys_km2'] == 0
    assert r['jura']['hybrid_share_of_the_region'] + r['jura']['wild_share_of_the_region'] == pytest.approx(1.0, abs=1e-3)
    assert len(r['heritage']) == 27
