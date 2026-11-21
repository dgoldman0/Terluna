"""Checks of the rivers' losses in the drylands: a channel's width and water surface, a cell's flow length, a valley's
floodplain by height above its river and distance from it, the nearest river found across the 0/360 seam, a wetland's
area capped by its inflow, and the product once it has been computed."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import river_losses as RL

CELL_KM = 1.9


def test_a_channel_is_widened_for_lunar_gravity():
    q = np.array([100.0, 10000.0])
    assert RL.LUNAR_WIDENING == pytest.approx(1.435, abs=0.005)
    np.testing.assert_allclose(RL.channel_width_m(q, 3.0), 3.0 * RL.LUNAR_WIDENING * np.sqrt(q))
    assert RL.channel_width_m(np.array([-5.0]), 3.0)[0] == 0.0


def test_flow_length_follows_the_receiver():
    rows, cols = 4, 8
    lat = np.array([60.0, 20.0, -20.0, -60.0])
    idx = np.arange(rows * cols)
    receiver = idx.copy()
    receiver[1 * cols + 2] = 1 * cols + 3            # east, at 20 N
    receiver[1 * cols + 4] = 2 * cols + 4            # south
    receiver[1 * cols + 5] = 2 * cols + 6            # south-east
    receiver[0 * cols + 7] = 0 * cols + 0            # east across the seam, at 60 N
    length = RL.flow_length_km(receiver, rows, cols, lat, CELL_KM)
    assert length[1 * cols + 2] == pytest.approx(CELL_KM * np.cos(np.radians(20.0)))
    assert length[1 * cols + 4] == pytest.approx(CELL_KM)
    assert length[1 * cols + 5] == pytest.approx(CELL_KM * np.hypot(1.0, np.cos(np.radians(20.0))))
    assert length[0 * cols + 7] == pytest.approx(CELL_KM * 0.5)
    assert length[0] == 0.0                         # a pit


def valley(rows=21, cols=40):
    """A river along row 10 at height 0, the ground rising 1 m a cell on both sides, near the equator."""
    lat = np.linspace(0.5, -0.5, rows)
    h = np.abs(np.arange(rows) - 10)[:, None] * np.ones(cols)
    river = np.zeros((rows, cols), bool)
    river[10] = True
    return h.astype(np.float32), river, lat


def test_a_floodplain_is_the_flat_land_near_its_river():
    h, river, lat = valley()
    near = RL.nearest_river(river)
    dry = np.ones(h.size, bool)
    land = ~river.ravel()
    for hand, rows_each_side in ((2.0, 2), (5.0, 5), (10.0, 5)):     # 10 km reaches five cells of 1.9 km
        plain = RL.floodplain(h, near, dry, land, lat, CELL_KM, hand).reshape(h.shape)
        assert plain.sum() == 2 * rows_each_side * h.shape[1], hand
        assert not plain[10].any()
    wet = dry.copy()
    wet[: 10 * h.shape[1]] = False                                  # the north bank is not dry land
    plain = RL.floodplain(h, near, wet, land, lat, CELL_KM, 5.0).reshape(h.shape)
    assert not plain[:10].any() and plain[11:16].all()


def test_the_nearest_river_is_found_across_the_seam():
    river = np.zeros((5, 200), bool)
    river[2, 0] = True
    near = RL.nearest_river(river, pad=8)
    assert near[2 * 200 + 199] == 2 * 200 + 0
    assert near[2 * 200 + 5] == 2 * 200 + 0


def test_a_wetland_settles_where_its_loss_meets_its_inflow():
    loss = np.array([1.0 / RL.SECONDS_PER_YEAR] * 2)                # 1 m a year
    area = np.array([100.0, 100.0])                                 # km2, needing 3.17 m3/s
    kept = RL.capped_area_km2(area, loss, np.array([10.0, 1.0]))
    assert kept[0] == 100.0
    assert kept[1] == pytest.approx(RL.SECONDS_PER_YEAR / 1e6)      # 1 m3/s keeps 31.6 km2 wet


@pytest.fixture(scope='module')
def product():
    if not RL.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.river_losses first')
    return json.loads(RL.OUT.read_text())


def test_the_product_holds_together(product):
    for budget in product['budgets'].values():
        lengths = [budget['rivers'][f'{q:.0f}']['length_km'] for q in RL.Q_CLASSES]
        assert lengths == sorted(lengths, reverse=True) and lengths[0] > 0
        narrow, wide = (budget['channel'][f'{c:g}']['evaporation_m3s'] for c in RL.WIDTH)
        assert 0 < narrow < wide
        plains = [budget['floodplains'][f'{h:g}']['area_km2'] for h in RL.HAND_M]
        assert plains == sorted(plains)
    closed, routing = product['budgets']['closed_budget'], product['budgets']['routing_budget']
    assert closed['dry_land_share_of_land'] >= routing['dry_land_share_of_land']
    at = product['floodplains_at']
    assert at['kept_wet_km2'] <= at['area_km2'] + 1
    for plain in at['largest']:
        assert plain['kept_wet_km2'] <= plain['area_km2'] + 1
    assert 0 < product['gcm']['mask_land_fraction_at_most'] < 1
