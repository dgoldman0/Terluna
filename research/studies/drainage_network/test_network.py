import json

import numpy as np
import pytest

from geography import drainage as dr
from research.studies.drainage_network import network as N
from research.studies.drainage_network import routing as RL

ROWS, COLS = 20, 40
LAT = np.linspace(1.9, -1.9, ROWS)              # near the equator: cells about square
CELL_KM = 1.0


def two_basins(ridge=30.0):
    """A deep basin in columns 3-12 and a shallow one in 25-34, a ridge between with one pass of height `ridge`."""
    col = np.arange(COLS)
    h = np.full((ROWS, COLS), 100.0)
    h[2:18, 3:13] = -200.0 + 5.0 * np.abs(col[3:13] - 7.5)[None, :]
    h[2:18, 25:35] = -50.0 + 5.0 * np.abs(col[25:35] - 29.5)[None, :]
    h[9:11, 13:25] = ridge                          # a channel between them, as high as the pass
    h += 0.001 * np.arange(ROWS)[:, None]
    return h


def slope_with_a_bowl(rows=24, cols=48):
    """Land rising both ways from a sea in columns 0-5 (longitude wraps), with a bowl on the slope."""
    col = np.arange(cols)
    distance = np.minimum(np.abs(col - 2.5), cols - np.abs(col - 2.5)) - 3.0
    h = np.repeat((np.maximum(distance, 0.0) * 10.0)[None, :], rows, axis=0) + 0.01 * np.arange(rows)[:, None]
    h[:, :6] = -50.0
    r, c = np.meshgrid(np.arange(rows), col, indexing='ij')
    h -= 80.0 * np.exp(-((r - 12) ** 2 + (c - 20) ** 2) / 8.0)
    return h, h < -1.0, np.linspace(2.0, -2.0, rows)


@pytest.mark.parametrize('lake_net', [0.5, -20.0])
def test_with_one_level_the_routing_is_geographys(lake_net):
    h, outlet, lat = slope_with_a_bowl()
    runoff, net = np.ones_like(h), np.full_like(h, lake_net)        # a bowl that spills, then one too dry to
    theirs = dr.route(h, outlet, lat, 1.0, runoff, net)
    ours = RL.route_levels(h, outlet, lat, 1.0, runoff, net, np.where(outlet, float(h[outlet].max()), 0.0))
    assert ours['outlet_nodes'] == 1
    for key in ('receiver', 'label', 'level', 'discharge', 'drained_m2', 'lake', 'closed', 'inflow', 'passed_on'):
        assert np.array_equal(ours[key], theirs[key]), key
    assert np.array_equal(np.isnan(ours['lake_level']), np.isnan(theirs['lake_level']))


def test_seas_at_different_levels_each_set_the_spill_of_their_own_hollows():
    # land rising both ways from a low sea in columns 0-5 (level -40 m), a higher sea sunk into it in columns 28-33
    # (level 22 m, its shore above 30 m), and a hollow on the low sea's slope whose pass lies between the two levels
    rows, cols = 24, 48
    col = np.arange(cols)
    away = np.minimum(np.abs(col - 2.5), cols - np.abs(col - 2.5))
    h = np.repeat((-40.0 + 4.0 * away)[None, :], rows, axis=0) + 0.01 * np.arange(rows)[:, None]
    low, high = away < 3.0, (col >= 28) & (col <= 33)
    h[:, low], h[:, high] = -50.0, 20.0
    r, c = np.meshgrid(np.arange(rows), col, indexing='ij')
    h -= 40.0 * np.exp(-((r - 12) ** 2 + (c - 12) ** 2) / 4.0)
    outlet = np.repeat((low | high)[None, :], rows, axis=0)
    level = np.where(low[None, :], -40.0, 22.0) * np.ones((rows, 1))
    lat = np.linspace(2.0, -2.0, rows)
    runoff, lake_net = np.ones_like(h), np.full_like(h, 0.5)
    one = dr.route(h, outlet, lat, 1.0, runoff, lake_net)
    two = RL.route_levels(h, outlet, lat, 1.0, runoff, lake_net, level)
    hollow = 12 * cols + 12
    assert two['outlet_nodes'] == 2
    assert one['level'][one['label'][hollow]] == float(h[outlet].max())        # one water: the higher sea's 20 m
    assert -40.0 < two['level'][two['label'][hollow]] < 0.0                    # its own pass toward the low sea
    assert two['level'][0] == -40.0 and two['level'][1] == 22.0                # the seas keep their levels
    land = ~two['outlet'] & ~two['lake']
    area = two['area']
    expected = float((np.where(land, runoff.ravel(), 0.0) * dr.MM_DAY * area).sum()
                     + (np.where(two['lake'], lake_net.ravel(), 0.0) * dr.MM_DAY * area).sum())
    assert np.isclose(two['discharge'][two['outlet']].sum(), expected, rtol=1e-9)   # all the water reaches the seas


def test_a_cell_joins_a_water_at_its_height_or_the_pass_on_the_way():
    h = two_basins(ridge=30.0)
    label, n, edges = N.catchment_graph(h, LAT, CELL_KM)
    seed = int(np.argmin(h))
    joined = N.joining_level(h.ravel(), label, n, edges, seed).reshape(h.shape)
    assert joined[10, 7] == pytest.approx(h[10, 7])                         # in the seed's basin: its own height
    assert joined[10, 29] == pytest.approx(30.0, abs=0.01)                  # beyond the pass: the pass's height
    assert joined[0, 0] == pytest.approx(100.0, abs=0.1)


def test_a_table_gives_area_loss_and_volume_at_a_level():
    joined = np.array([-10.0, -5.0, 0.0, 5.0])
    h = joined.astype(np.float32)
    area = np.full(4, 2.0)
    loss = np.array([1.0, 1.0, -0.5, 1.0])
    t = N.Table(joined, area, loss, h, cap=100.0)
    a, l_, v = t.at(1.0)                           # three cells under 1 m: depths 11, 6 and 1 m
    assert (a, l_, v) == pytest.approx((6.0, 1.5, 2.0 * (11.0 + 6.0 + 1.0)))
    assert t.level_for(1.6, below=100.0) == pytest.approx(-5.0)           # the second cell's loss takes it past 1.6
    assert t.level_for(1.6, below=-5.0) is None                             # not under -5 m
    assert t.level_for(5.0, below=100.0) is None


def test_a_basin_of_one_area_fills_in_its_volume_over_its_surplus():
    joined = np.linspace(-100.0, 0.0, 101)          # 101 cells of 1 km2 at every metre: about a uniform area
    area = np.full(101, 1e6)
    loss = np.full(101, 1e-3)                      # m3/s a cell
    t = N.Table(joined, area, loss, joined.astype(np.float32), cap=10.0)
    supply = 1.0
    years, settles = N.fill_history(t, supply, -50.0, -1.0)
    v = t.at(-1.0)[2] - t.at(-50.0)[2]
    assert settles is None
    assert years == pytest.approx(v / (supply - 0.075) / N.YEAR_S, rel=0.05)
    years, settles = N.fill_history(t, 0.04, -100.0, -1.0)                # a supply it meets near 40 cells
    assert settles == pytest.approx(-61.0, abs=2.0)


def test_a_long_profile_measures_the_rivers_fall():
    cols = 40
    lat_rows = np.zeros(10)
    path = np.array([5 * cols + c for c in range(10)])                    # ten cells along a row, mouth first
    surface = np.zeros(10 * cols, np.float32)
    surface[path] = [0.0, 1.0, 2.0, 60.0, 61.0, 62.0, 62.0, 62.0, 63.0, 64.0]
    lake = np.zeros(10 * cols, bool)
    lake[path[5:8]] = True
    q = np.zeros(10 * cols)
    q[path] = np.linspace(100.0, 10.0, 10)
    p = N.long_profile(path, surface, q, lat_rows, cols, 1.0, 0.0, lake)
    assert p['length_km'] == 9 and p['surface_from_m'] == 0 and p['surface_to_m'] == 64
    assert p['steep_reaches'] == 1 and p['steep_km'] == 1 and p['largest_fall_between_cells_m'] == 58
    assert p['lakes_passed'] == 1 and p['discharge_m3s'] == [100, 10]


def test_an_upper_sea_fills_and_its_overflow_raises_the_one_below():
    # the upper sea's rivers exceed its loss; the lower one's do not until the upper overflows into it
    def table(depth, area_m2, loss_m3s):
        joined = np.linspace(-depth, 0.0, 201)
        return N.Table(joined, np.full(201, area_m2 / 201), np.full(201, loss_m3s / 201), joined.astype(np.float32), 50.0)
    tables = dict(upper=table(100.0, 1e9, 1.0), lower=table(100.0, 1e10, 50.0))
    rims = dict(upper=0.0, lower=0.0)
    out = dict(upper=dict(rivers_at_the_common_level_m3s=20.0, rivers_at_the_rim_m3s=20.0, overflow_to='lower'),
               lower=dict(rivers_at_the_common_level_m3s=30.0, rivers_at_the_rim_m3s=30.0, overflow_to='upper'))
    h = N.joined_history(['upper', 'lower'], tables, rims, out, -100.0, years=2000.0, step=1.0)
    first = h['events'][0]
    assert first['event'].startswith('upper is full')
    v_up = tables['upper'].at(-N.BELOW_RIM_M)[2] - tables['upper'].at(-100.0)[2]
    assert first['year'] == pytest.approx(v_up / (20.0 - 0.5) / N.YEAR_S, rel=0.1)


def test_the_product_holds_its_chain_and_budget():
    path = N.OUT
    if not path.exists():
        pytest.skip('results/network.json has not been written')
    out = json.loads(path.read_text())
    assert out['schema'] == N.SCHEMA
    assert [c['rain'] for c in out['cases']] == [N.RAINS[0], N.RAINS[0]] + list(N.RAINS[1:])
    closed = out['cases'][1]['budget']
    assert out['cases'][1]['evaporation'] == out['cases'][0]['budget']['open_water_evaporation_factor_closing_it']
    assert abs(closed['surplus_m3s']) < 0.1 * closed['seas_loss_m3s']        # the raised evaporation all but closes it
    for c in out['cases']:
        b = c['budget']
        assert b['surplus_m3s'] == b['delivered_to_the_seas_m3s'] - b['seas_loss_m3s']
        assert b['routing_runoff_m3s'] >= b['land_rain_less_evaporation_m3s']
        for name, s in c['seas'].items():
            assert s['area_at_the_rim_km2'] >= s['atlas_area_km2'] and s['overflow_to'] in c['seas']
            assert s['fills'] == (s['years_to_fill'] is not None)
        years = [e['year'] for e in c['history']['events']]
        assert years == sorted(years)
