"""Checks of the atlas round's governed drainage: a branch's put-ins and take-outs along a route of cuts, a take-out
acting on the lake the water passes through, the weirs setting each held sea at its need, and the product once it has
been computed."""
import json
from types import SimpleNamespace

import numpy as np
import pytest

from research.studies.atlas_round import governed_drainage as GD
from research.studies.atlas_round import redrawn_drainage as RD
from research.studies.drainage_network import network as NW

ROWS, COLS = 24, 48
LAT = np.linspace(2.0, -2.0, ROWS)                # near the equator: cells about square
CELL_KM = 1.0


def chain_search():
    """Six catchments: a full sea at 100 m (0), slope ground whose hollow would stand at 130 m (1), hollows standing
    at 80 m (2), 60 m (3) and 40 m (4), and a sea at 0 m (5); edges 0-1 (pass 120), 1-2 (110), 2-3 (80), 3-4 (75),
    4-5 (40) and 0-5 (130)."""
    lev = [100.0, 130.0, 80.0, 60.0, 40.0, 0.0]
    src = np.array([0, 1, 2, 3, 4, 0])
    dst = np.array([1, 2, 3, 4, 5, 5])
    wt = [120.0, 110.0, 80.0, 75.0, 40.0, 130.0]
    c_to = np.array([101, 102, 103, 104, 105, 106])
    return SimpleNamespace(lev_l=lev, wt_l=wt, src=src, dst=dst, c_to=c_to)


def test_a_branch_is_put_in_where_it_runs_on_by_itself_and_taken_out_where_a_cut_leaves_a_lake():
    search = chain_search()
    sinks = np.array([10, 11, 12, 13, 14, 15])
    in_sea = np.array([0, -1, -1, -1, -1, 1], np.int8)
    names = ['full', 'below']
    # through the slope's ground at 100 m (cuts of 20 and 10 m), into the 80 m hollow, on by itself into the 60 m
    # hollow, out of it by a 15 m cut into the 40 m hollow, and on by itself into the sea
    water = np.zeros(110, bool)
    water[106] = True                                                      # the 0-5 pass opens onto the sea's water
    chain = [(0, 100.0), (1, 100.0), (2, 80.0), (3, 60.0), (4, 40.0)]
    ops, onto, cuts, free = GD.branch_operations(chain, search, sinks, in_sea, names, water)
    assert ops == [(102, +1, None), (13, -1, 60.0), (104, +1, None)]
    assert onto is None and free == [102, 104]
    assert [(e, round(d, 1)) for e, d in cuts] == [(0, 20.0), (1, 10.0), (3, 15.0)]
    # a cut that opens straight onto the sea puts nothing into the rivers
    ops, onto, cuts, free = GD.branch_operations([(5, 100.0)], search, sinks, in_sea, names, water)
    assert ops == [] and onto == 'below' and free == []


def three_seas():
    """A low sea in columns 0-5, land rising east of it with a hollow on the slope, a walled basin in columns 24-30 whose
    lowest pass (35 m, row 12, column 23) opens west onto the slope and a second pass (45 m, column 31) east onto land
    falling to another sea in columns 40-44, and a ridge in columns 45-47 between that sea and the low one."""
    col = np.arange(COLS)
    h = np.empty((ROWS, COLS))
    h[:, :6] = -50.0
    h[:, 6:23] = 2.0 * (col[6:23] - 5)[None, :]
    h[:, 23:32] = 60.0
    h[3:21, 24:31] = 10.0
    h[12, 23] = 35.0
    h[12, 31] = 45.0
    h[:, 32:40] = (44.0 - 2.5 * (col[32:40] - 32))[None, :]
    h[:, 40:45] = -50.0
    h[:, 45:] = 80.0
    r, c = np.meshgrid(np.arange(ROWS), col, indexing='ij')
    h -= 6.0 * ((r == 12) & (c >= 11) & (c <= 13))
    h += 0.001 * np.arange(ROWS)[:, None] + 0.0001 * col[None, :]
    runoff = np.ones_like(h)
    runoff[:, 23:32] = 30.0                                                # a wet basin
    lake_net = np.full_like(h, -5.0)                                       # dry hollows
    lake_net[:, :6] = lake_net[:, 40:45] = -15.0                           # and seas drier still
    lake_net[:, 23:32] = -0.1
    return h.astype(np.float32), runoff, lake_net


def small_geo(h, level0, runoff, lake_net, seeds, upper):
    label, sinks, edges = RD.catchment_graph(h, LAT, CELL_KM)
    flat = h.ravel()
    seas = list(seeds)
    joined = {k: NW.joining_level(flat, label, sinks.size, edges, seeds[k]) for k in seas}
    outlet = np.logical_or.reduce([joined[k] < level0 for k in seas if k != upper])
    r = RD.dr.route(h, outlet.reshape(h.shape), LAT, CELL_KM, runoff, lake_net)
    path = [seeds[upper]]
    while not outlet[path[-1]] and len(path) < h.size:
        path.append(int(r['receiver'][path[-1]]))
    passes = {(a, b): float(joined[a][seeds[b]]) for a in seas for b in seas if a != b}
    return dict(level0=level0, h=h, lat=LAT, cell_km=CELL_KM, seas=seas, joined=joined, seeds=seeds,
                spill_paths={upper: (np.array(path), None)}, passes=passes, catchments=sinks.size,
                graph=(label, sinks, edges))


SEEDS = dict(high=12 * COLS + 27, low=12 * COLS + 2, other=12 * COLS + 42)


def test_the_weirs_send_the_held_sea_its_need_and_the_rest_through_the_cut():
    h, runoff, lake_net = three_seas()
    geo = small_geo(h, -1.0, runoff, lake_net, SEEDS, 'high')
    state = RD.design_state(geo, runoff, lake_net, upper=('high',), lower=('low', 'other'))
    levels, in_sea = RD.node_levels(geo, state)
    search = RD.CutSearch(geo, levels, in_sea, state['names'], 'high')
    allowed = search.shallowest(only=('other',))['other']
    assert allowed == pytest.approx(10.0, abs=RD.CUT_PRECISION_M)          # the east pass, 10 m above the water
    chain, tidy = search.route('other', allowed)
    _, sinks, _ = geo['graph']
    ops, onto, cuts, free = GD.branch_operations(chain, search, sinks, in_sea, state['names'], state['sea_of'] >= 0)
    assert tidy and onto is None and len(free) == 1 and all(sign > 0 for _, sign, _ in ops)
    branches = {'other': ([(c, s) for c, s, _ in ops], onto)}
    gov = GD.govern(geo, state, runoff, lake_net, branches, source='high', salt_sea='low', held=None, rest='other',
                    full_seas=(), settle=1e-9, log=lambda *a: None)
    loss, own = state['loss'], state['own']
    assert loss['low'] > own['low']                                        # the low sea needs the basin's water
    assert gov['supply']['low'] == pytest.approx(loss['low'], abs=1e-6)    # and gets just that over the natural pass
    rest = gov['surplus'] - gov['share']['low']
    assert gov['share']['other'] == pytest.approx(rest)
    assert gov['supply']['other'] == pytest.approx(own['other'] + rest, rel=1e-6)


def test_a_take_out_at_a_lake_the_water_passes_leaves_the_sea_below_as_it_was():
    h, runoff, lake_net = three_seas()
    geo = small_geo(h, -1.0, runoff, lake_net, SEEDS, 'high')
    state = RD.design_state(geo, runoff, lake_net, upper=('high',), lower=('low', 'other'))
    r = state['r']
    put = 12 * COLS + 20                                                   # on the slope above the hollow
    lake = GD.lake_on_the_way(r, h.ravel(), put, float(r['level'][r['label'][12 * COLS + 11]]))
    assert lake is not None and lake // COLS == 12 and 11 <= lake % COLS <= 13
    area = np.repeat((CELL_KM * 1e3) ** 2 * np.cos(np.radians(LAT)), COLS)
    lake_full = ~r['outlet'] & (h.ravel() < r['level'][r['label']])
    waters = [(state['masks'][k], state['levels'][k]) for k in state['names']]
    plain, base = GD.route_with(geo, runoff, lake_net, waters, state['names'], [], area, lake_full)
    routed, with_ops = GD.route_with(geo, runoff, lake_net, waters, state['names'], [(put, 2.0), (lake, -2.0)], area, lake_full)
    assert with_ops['low'] == pytest.approx(base['low'], abs=1e-9)
    assert routed['discharge'][12 * COLS + 15] == pytest.approx(plain['discharge'][12 * COLS + 15] + 2.0, abs=1e-6)


@pytest.fixture(scope='module')
def product():
    if not GD.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.governed_drainage first')
    return json.loads(GD.OUT.read_text())


def test_the_held_seas_get_their_need_and_the_rest_reaches_smythii_marginis(product):
    for case in product['cases'].values():
        seas = case['seas']
        assert abs(seas[GD.SALT_SEA]['balance_m3s']) <= 2 and abs(seas[GD.HELD_BRANCH]['balance_m3s']) <= 2
        released, arrived = case['released_m3s'], case['arrived_m3s']
        assert sum(released.values()) == pytest.approx(case['moscoviense_surplus_m3s'], abs=2)
        for k in released:
            assert arrived[k] <= released[k] + 2                            # the lakes on the way only take
        assert seas[GD.REST_BRANCH]['balance_m3s'] == pytest.approx(case['full_seas_pass_on_m3s'][GD.REST_BRANCH], abs=2)
