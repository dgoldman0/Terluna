"""Checks of the atlas round's drainage on the redrawn floor: the design state's overflows, the single cuts out of a full
sea, the least-cut routes, and the product once it has been computed."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.drainage_network import network as NW

ROWS, COLS = 24, 48
LAT = np.linspace(2.0, -2.0, ROWS)                # near the equator: cells about square
CELL_KM = 1.0


def terrain(east_pass=None):
    """A low sea in columns 0-5 (wrapping round), land rising east of it, and a walled basin in columns 24-30 whose
    lowest pass, 35 m at row 12, column 23, opens west onto the slope down to the sea, with a hollow on the way.
    `east_pass` opens a second, higher pass in the basin's east wall onto the land that runs round to the sea's
    other shore."""
    col = np.arange(COLS)
    h = np.empty((ROWS, COLS))
    h[:, :6] = -50.0
    h[:, 6:23] = 2.0 * (col[6:23] - 5)[None, :]                           # 2 m at column 6 to 34 m at column 22
    h[:, 23:32] = 60.0                                                     # the basin's walls
    h[3:21, 24:31] = 10.0                                                  # its floor
    h[12, 23] = 35.0                                                       # its lowest pass
    h[:, 32:] = (44.0 - 2.5 * (col[32:] - 32))[None, :]                    # 44 m at column 32 to 6.5 m at column 47
    if east_pass is not None:
        h[12, 31] = east_pass
    r, c = np.meshgrid(np.arange(ROWS), col, indexing='ij')
    h -= 6.0 * ((r == 12) & (c >= 11) & (c <= 13))                         # a hollow on the western slope
    h += 0.001 * np.arange(ROWS)[:, None] + 0.0001 * col[None, :]
    return h.astype(np.float32)


def small_geo(h, level0, runoff, lake_net, seeds, upper='high', lower='low'):
    """The pieces of redrawn_drainage.geometry on a small grid, for named seeds: the upper sea's spill path found as
    network.spill_targets finds it, through a routing in which the lower sea takes the rivers in."""
    label, sinks, edges = RD.catchment_graph(h, LAT, CELL_KM)
    flat = h.ravel()
    seas = list(seeds)
    joined = {k: NW.joining_level(flat, label, sinks.size, edges, seeds[k]) for k in seas}
    outlet = joined[lower] < level0
    r = RD.dr.route(h, outlet.reshape(h.shape), LAT, CELL_KM, runoff, lake_net)
    path = [seeds[upper]]
    while not outlet[path[-1]] and len(path) < h.size:
        path.append(int(r['receiver'][path[-1]]))
    passes = {(a, b): float(joined[a][seeds[b]]) for a in seas for b in seas if a != b}
    return dict(level0=level0, h=h, lat=LAT, cell_km=CELL_KM, seas=seas, joined=joined, seeds=seeds,
                targets={upper: lower}, spill_paths={upper: (np.array(path), None)}, passes=passes,
                catchments=sinks.size, graph=(label, sinks, edges))


SEEDS = dict(high=12 * COLS + 27, low=12 * COLS + 2)


def test_a_full_seas_surplus_runs_from_past_its_pass_through_the_hollow_to_the_sea_below():
    h = terrain()
    runoff = np.ones_like(h)
    lake_net = np.full_like(h, -0.1)
    geo = small_geo(h, -1.0, runoff, lake_net, SEEDS)
    state = RD.design_state(geo, runoff, lake_net, upper=('high',), lower=('low',))
    assert abs(state['levels']['high'] - (35.0 + 0.0001 * 23 + 0.012 - NW.BELOW_RIM_M)) < 0.01
    assert state['entry']['high'] == 12 * COLS + 22                        # the first cell past the pass, on the slope
    surplus = state['own']['high'] - state['loss']['high']
    assert surplus > 0 and state['put_in']['high'] == pytest.approx(surplus)
    reaches, path = state['paths']['high']
    assert reaches == 'low'
    assert {12 * COLS + 12} <= set(path.tolist())                          # it runs through the hollow
    # the sea below receives the surplus, less what the hollow's water loses to the air
    hollow = np.flatnonzero(state['r']['lake'])
    hollow_loss = float((0.1 * RD.dr.MM_DAY * state['r']['area'])[hollow].sum())
    assert state['supply']['low'] == pytest.approx(state['own']['low'] + surplus, abs=hollow_loss + 1e-9)
    assert state['supply']['low'] > state['own']['low'] + 0.9 * surplus


def test_a_higher_pass_out_of_a_full_sea_is_listed_with_the_cut_it_needs():
    h = terrain(east_pass=45.0)
    runoff = np.ones_like(h)
    lake_net = np.full_like(h, -0.1)
    geo = small_geo(h, -1.0, runoff, lake_net, SEEDS)
    state = RD.design_state(geo, runoff, lake_net, upper=('high',), lower=('low',))
    ways = RD.single_cuts(geo, state, None, upper=('high',), deepest=20.0, same_gap_km=12.0)['high']
    assert len(ways) == 2 and all(w['reaches'] == 'low' for w in ways)
    assert ways[0]['cut_m'] == pytest.approx(0.0, abs=0.05)                # the natural overflow
    assert ways[1]['cut_m'] == pytest.approx(10.0, abs=0.1)                # the east wall's pass, 10 m above the water
    assert ways[0]['mouth'] != ways[1]['mouth']


def test_the_least_cut_route_goes_round_by_the_shallower_cut():
    # five catchments: a full sea at 100 m (0), hollows filling to 60 m (1) and 50 m (2), and two seas at 0 m (3, 4).
    # To sea 3 the water overflows naturally through hollow 1; to sea 4 it could cut 30 m straight into hollow 2, or
    # overflow into hollow 1 and cut 20 m from there into hollow 2.
    rows, cols = 4, 8
    cell = lambda r, c: r * cols + c
    u = np.array([0, 1, 0, 2, 1])
    v = np.array([1, 3, 2, 4, 2])
    w = np.array([100.005, 60.0, 130.0, 50.0, 80.0], np.float32)
    cu = np.array([cell(0, 1), cell(1, 1), cell(0, 2), cell(2, 2), cell(1, 3)])
    cv = np.array([cell(1, 0), cell(2, 0), cell(2, 3), cell(3, 3), cell(2, 4)])
    sinks = np.array([cell(0, 0), cell(1, 0), cell(2, 2), cell(3, 0), cell(3, 6)])
    geo = dict(graph=(None, sinks, (u, v, w, cu, cv)), h=np.zeros((rows, cols)), lat=np.linspace(1.5, -1.5, rows),
               level0=0.0)
    state = dict(names=['full', 'near', 'far'])
    levels = np.array([100.0, 60.0, 50.0, 0.0, 0.0])
    in_sea = np.array([0, -1, -1, 1, 2], np.int8)
    routes = RD.least_cut_routes(geo, state, levels, in_sea, 'full', None)
    assert routes['near']['deepest_cut_m'] == pytest.approx(0.0, abs=0.01) and routes['near']['cuts'] == []
    assert routes['far']['deepest_cut_m'] == pytest.approx(20.0)
    assert [c['depth_m'] for c in routes['far']['cuts']] == [20.0]
    assert routes['far']['cuts'][0]['pass_above_the_common_level_m'] == 80.0
    assert routes['far']['passes'] == 3


@pytest.fixture(scope='module')
def product():
    if not RD.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.redrawn_drainage first')
    return json.loads(RD.OUT.read_text())


def test_the_redrawn_floor_reproduces_the_first_steps_rims(product):
    assert all(v['agree'] for v in product['rims_against_the_first_step'].values())
    assert product['common_level']['against_the_atlas_m'] == pytest.approx(39.17, abs=0.05)


def test_the_full_seas_spill_and_their_overflows_reach_the_seas_below(product):
    for budget in ('routing_budget', 'closed_budget'):
        seas = product['design'][budget]['seas']
        arrived = sum(seas[k]['from_the_full_seas_m3s'] for k in RD.LOWER)
        sent = sum(seas[k]['passes_on_m3s'] for k in RD.UPPER)
        for k in RD.UPPER:
            assert seas[k]['balance_m3s'] > 0 and seas[k]['reaches'] in RD.LOWER + RD.UPPER
        # what the full seas send arrives below, less what the hollows on the way lose to the air
        assert 0.95 * sent <= arrived <= sent + 1


def test_the_network_cases_hold_together(product):
    assert [c['rain'] for c in product['cases']] == [1.0, 1.0, 0.75, 2.0]
    closed = product['cases'][1]['budget']
    assert abs(closed['surplus_m3s']) < 0.05 * closed['seas_loss_m3s']
