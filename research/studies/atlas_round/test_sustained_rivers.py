"""Checks of the sustained rivers: on small grids, the replay of the routing reproduces route_levels exactly, the
catchments' cascade gives each water what the replay gives it, a river's loss is taken where it falls and reaches the
lakes below, and the storage holds still on a steady year, keeps its water and carries a lake's outflow through a dry
year; Moscoviense's put-ins follow its fixed shares; a river's stems follow the larger branch, pass small lakes and end
at large ones; and, once computed, the product and its record."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import sustained_rivers as SR
from research.studies.drainage_network.routing import route_levels

CELL_KM = 1.9
YEAR_S = SR.YEAR_S


def basin(rows=24, cols=48):
    """A slope rising to the east, a bowl on it that fills to a lake and spills west, a sea along the west edge and the
    six other waters as small blocks on high ground, each at its own level; with runoff (mm/day) on the land and the
    open water losing 3 mm/day beyond its rain."""
    lat = np.linspace(2.0, -2.0, rows)
    c = np.arange(cols)[None, :] * np.ones((rows, 1))
    r = np.arange(rows)[:, None] * np.ones((1, cols))
    h = 500.0 + 12.0 * c + 0.37 * r
    h -= 160.0 * np.exp(-((r - 12.0) ** 2 / 18.0 + (c - 24.0) ** 2 / 30.0))
    water = np.full(rows * cols, -1, np.int8)
    levels = [0.0] * 7
    blocks = {3: (slice(0, rows), slice(0, 2))}                 # the salt sea along the west edge
    for i, k in enumerate((0, 1, 2, 4, 5, 6)):
        blocks[k] = (slice(0, 2), slice(34 + 2 * i, 36 + 2 * i))
    for k, (rs, cs) in blocks.items():
        cell = np.zeros((rows, cols), bool)
        cell[rs, cs] = True
        water[cell.ravel()] = k
        levels[k] = float(h[cell].max()) + 1.0 + k
    runoff = np.full((rows, cols), 1.0)
    lake_net = np.full((rows, cols), -3.0)
    return h.astype(np.float32), water, levels, lat, runoff, lake_net


def structure(net, neck=0.0):
    """A design structure with no put-ins or take-outs, for the small grid."""
    meta = dict(fractions={'south_pole_aitken': 0.5, 'humboldtianum': 0.2, 'smythii_marginis': 0.3},
                neck_flow=neck, level0=0.0)
    empty = np.zeros(0)
    return dict(meta=meta, op_cell=empty.astype(np.int64), op_kind=empty.astype(np.int8), op_sign=empty,
                op_lake=empty.astype(bool), op_design=empty)


@pytest.fixture(scope='module')
def small():
    h, water, levels, lat, runoff, lake_net = basin()
    net = SR.Network(h, water, levels, lat, CELL_KM, log=lambda *a: None)
    return net, h, water, levels, lat, runoff, lake_net


def test_the_replay_is_route_levels(small):
    net, h, water, levels, lat, runoff, lake_net = small
    outlet = (water >= 0).reshape(h.shape)
    level = np.where(water >= 0, np.asarray(levels)[np.maximum(water, 0)], 0.0).reshape(h.shape)
    ref = route_levels(h, outlet, lat, CELL_KM, runoff, lake_net, level)
    q, out = net.replay(*net.sources(runoff, lake_net))
    assert np.array_equal(q, ref['discharge'])
    assert np.array_equal(net.receiver, ref['receiver'])
    assert np.array_equal(out[net.n_out:], ref['passed_on'][net.n_out:])
    lake = ref['lake'] & ~ref['outlet']
    assert lake.sum() > 20                                       # the bowl holds a lake
    assert ref['passed_on'][net.n_out:].max() > 0.0              # and spills


def test_the_cascade_gives_each_water_what_the_replay_does(small):
    net, h, water, levels, lat, runoff, lake_net = small
    runoff_c, lake = net.sources(runoff, lake_net)
    q, out = net.replay(runoff_c, lake)
    _, out_c, supply = net.cascade(runoff_c, lake)
    np.testing.assert_allclose(out_c, out, rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(supply, net.supply_of(q), rtol=1e-12, atol=1e-9)
    assert supply[3] > 0.0
    q2, _ = net.replay(runoff_c, lake, passed_on=out)           # each pit passing on what it passes on anyway
    np.testing.assert_allclose(q2, q, rtol=1e-12, atol=1e-9)


def test_a_loss_is_taken_where_it_falls_and_never_more_than_the_river_carries(small):
    net, h, water, levels, lat, runoff, lake_net = small
    runoff_c, lake = net.sources(runoff, lake_net)
    q, out = net.replay(runoff_c, lake)
    k = int(np.argmax(out))                                     # the bowl's lake
    above = int(np.flatnonzero(net.receiver == net.sinks[k - net.n_out])[0])
    while net.lake_full[above]:                                 # a cell on the ground upstream of the lake
        above = int(np.flatnonzero(net.receiver == above)[0])
    loss = np.zeros(net.n)
    loss[above] = 0.25 * q[above]
    q_l, out_l = net.replay(runoff_c, lake, loss=loss)
    assert out_l[k] == pytest.approx(out[k] - loss[above], rel=1e-9)
    assert net.supply_of(q_l)[3] == pytest.approx(net.supply_of(q)[3] - loss[above], rel=1e-9)
    loss[above] = 10.0 * q[above]                               # more than it carries: it runs dry there
    q_l, _ = net.replay(runoff_c, lake, loss=loss)
    assert q_l[net.receiver[above]] == pytest.approx(q[net.receiver[above]] - q[above], rel=1e-9, abs=1e-9)


def test_the_head_over_an_outlet_balances_the_step():
    area, c, s = np.array([1e9, 1e6, 5e8]), np.array([1e7, 1e9, 0.0]), np.array([3e9, 2e7, 1e9])
    h = SR.solve_head(area, c, s)
    np.testing.assert_allclose(area * h + c * h ** 1.5, s, rtol=1e-10)
    assert h[2] == pytest.approx(2.0)                           # no outlet: the water only rises


def test_a_row_table_interpolates_as_numpy_does():
    xs = np.array([[0.0, -1.0, -3.0, -6.0], [0.0, -2.0, -2.0, -5.0]])
    ys = np.array([[-10.0, -8.0, -4.0, 0.0], [-5.0, -4.0, -4.0, -1.0]])
    x = np.array([-2.0, -4.0])
    got = SR.row_interp(x, xs, ys)
    assert got[0] == pytest.approx(np.interp(2.0, -xs[0], ys[0]))
    assert got[1] == pytest.approx(np.interp(4.0, [2.0, 5.0], [-4.0, -1.0]))
    assert SR.row_interp(np.array([-9.0]), xs[:1], ys[:1])[0] == 0.0       # below the table: empty


def test_moscovienses_put_ins_follow_its_fixed_shares():
    kinds = np.array([SR.OPS.index(k) for k in ('salt_sea_weir', 'humboldtianum_branch', 'smythii_marginis_branch',
                                                'smythii_marginis_branch', 'orientale', 'neck')], np.int8)
    s = dict(meta=dict(fractions={'south_pole_aitken': 0.6, 'humboldtianum': 0.15, 'smythii_marginis': 0.25}),
             op_kind=kinds, op_sign=np.array([1.0, 1.0, 1.0, -1.0, 1.0, 1.0]))
    q = SR.op_flows(s, 1000.0, 300.0, 400.0, 50.0)
    np.testing.assert_allclose(q, [600.0, 150.0, 250.0, -250.0, 300.0, 50.0])


@pytest.fixture(scope='module')
def stored(small):
    net, h, water, levels, lat, runoff, lake_net = small
    s = structure(net)
    design = SR.without_storage(net, s, runoff, lake_net)
    runoff_c, lake = net.sources(runoff, lake_net)
    design['inflow'] = net.cascade(runoff_c, lake)[0]
    return net, s, design, runoff, lake_net


def test_the_storage_holds_still_on_the_design_year(stored):
    net, s, design, runoff, lake_net = stored
    st = SR.Storage(net, s, design, SR.RL.WIDTH[0], log=lambda *a: None)
    st.initial(design, lake_net, design['inflow'])
    res = st.year(runoff, lake_net)
    lakes = st.lake_nodes
    np.testing.assert_allclose(res['out'][lakes], design['out'][lakes], rtol=1e-6, atol=1e-6)
    q = st.field(runoff, lake_net, res)
    np.testing.assert_allclose(q, design['q'], rtol=1e-6, atol=1e-6)


def test_the_storage_keeps_its_water_and_carries_a_lake_through_a_dry_year(stored):
    net, s, design, runoff, lake_net = stored
    st = SR.Storage(net, s, design, SR.RL.WIDTH[0], log=lambda *a: None)
    st.initial(design, lake_net, design['inflow'])
    k = int(np.argmax(np.where(np.arange(net.n_nodes) >= net.n_out, design['out'], 0.0)))
    i = st.slot[k]
    wet = 1.5 * runoff
    runoff_c, lake = net.sources(wet, lake_net)
    inflow = net.cascade(runoff_c, lake)[0][k]
    v0 = st.vol[i]
    res = st.year(wet, lake_net)
    assert st.vol[i] > 0.0                                       # above its spill level all year
    assert st.vol[i] - v0 == pytest.approx(YEAR_S * (inflow + lake[k] - res['out'][k]), rel=1e-6)
    dry = 0.3 * runoff
    res = st.year(dry, lake_net)
    plain = SR.without_storage(net, s, dry, lake_net)
    assert res['out'][k] > plain['out'][k] + 0.5                 # the lake's head carries water into the dry year


def twin_hollows(rows=24, cols=48):
    """A plateau with two hollows joined by a saddle under the notch through which the western one spills to a slope
    down to the sea, and the six other waters as blocks on the plateau."""
    lat = np.linspace(2.0, -2.0, rows)
    c = np.arange(cols)[None, :] * np.ones((rows, 1))
    r = np.arange(rows)[:, None] * np.ones((1, cols))
    h = 1000.0 + 0.01 * c + 0.003 * r
    h[:, :8] = 100.0 + 110.0 * c[:, :8] + 0.003 * r[:, :8]
    for centre in (19.0, 27.0):
        h -= 200.0 * np.exp(-((r - 12.0) ** 2 + (c - centre) ** 2) / 8.0)
    h[12, 8:15] = np.linspace(930.0, 975.0, 7)
    water = np.full(rows * cols, -1, np.int8)
    levels = [0.0] * 7
    blocks = {3: (slice(0, rows), slice(0, 1))}
    for i, k in enumerate((0, 1, 2, 4, 5, 6)):
        blocks[k] = (slice(0, 2), slice(34 + 2 * i, 36 + 2 * i))
    for k, (rs, cs) in blocks.items():
        cell = np.zeros((rows, cols), bool)
        cell[rs, cs] = True
        water[cell.ravel()] = k
        levels[k] = float(h[cell].max()) + 1.0 + k
    return h.astype(np.float32), water, levels, lat, np.full((rows, cols), 2.0), np.full((rows, cols), -1.0)


def test_a_lake_of_two_hollows_stores_as_one_body():
    h, water, levels, lat, runoff, lake_net = twin_hollows()
    net = SR.Network(h, water, levels, lat, CELL_KM, log=lambda *a: None)
    inner = np.flatnonzero(net.inner)
    assert inner.size >= 1                                       # the eastern hollow lies under the western's level
    k = int(net.body[inner[0]])
    assert not net.inner[k] and net.level[inner[0]] == net.level[k]
    s = structure(net)
    design = SR.without_storage(net, s, runoff, lake_net)
    runoff_c, lake = net.sources(runoff, lake_net)
    design['inflow'] = net.cascade(runoff_c, lake)[0]
    st = SR.Storage(net, s, design, SR.RL.WIDTH[0], log=lambda *a: None)
    node_area = np.bincount(net.label[net.lake_full], weights=net.area[net.lake_full], minlength=net.n_nodes)
    assert st.slot[inner].max() < 0                              # the inner catchments hold no volume of their own
    assert st.area_head[st.slot[k]] == pytest.approx(node_area[net.body == k].sum())
    assert st.area_head[st.slot[k]] > node_area[k]
    st.initial(design, lake_net, design['inflow'])
    res = st.year(runoff, lake_net)
    np.testing.assert_allclose(res['out'], design['out'], rtol=1e-6, atol=1e-6)


def path_grid():
    """A river of ten cells into a sea through a small lake, with a tributary: cells 0-2 on land, 3-4 the lake (4 its
    pit, which spills to 5), 5-6 on land, 7 the sea; 8 feeds 2."""
    n = 10
    receiver = np.array([1, 2, 3, 4, 5, 6, 7, 7, 2, 9])
    q = np.array([150.0, 160.0, 290.0, 0.0, 0.0, 280.0, 290.0, 0.0, 120.0, 0.0])
    wet = np.zeros(n, bool)
    wet[[3, 4]] = True
    outlet = np.zeros(n, bool)
    outlet[7] = True
    water = np.full(n, -1, np.int8)
    water[7] = 3
    sink = np.zeros(n, bool)
    sink[4] = True
    node_of = np.full(n, -1)
    node_of[4] = 5
    body = np.zeros(n, np.int64)
    body[[3, 4]] = 1
    return receiver, q, wet, outlet, water, sink, node_of, body


def test_a_river_runs_through_a_small_lake_and_ends_at_a_large_one():
    receiver, q, wet, outlet, water, sink, node_of, body = path_grid()
    cells = np.array([0, 1, 2, 5, 6, 8])
    passes = np.zeros(6)
    passes[5] = 300.0
    down = SR.next_down(cells, q, wet, outlet, water, receiver, sink, node_of, passes, body, np.array([0.0, 40.0]), 100.0)
    assert down[2] == ('river', 5) and down[6] == ('sea', 3) and down[8] == ('river', 2)
    stems = SR.stems(cells, q, down)
    by_head = {s['cells'][0]: s for s in stems}
    assert by_head[0]['cells'] == [0, 1, 2, 5, 6] and by_head[0]['ends'] == ('sea', 3)
    assert by_head[8]['cells'] == [8] and by_head[8]['ends'] == ('joins', 2)
    big = SR.next_down(cells, q, wet, outlet, water, receiver, sink, node_of, passes, body, np.array([0.0, 400.0]), 100.0)
    assert big[2] == ('lake', 1)
    passes[5] = 0.0
    shut = SR.next_down(cells, q, wet, outlet, water, receiver, sink, node_of, passes, body, np.array([0.0, 40.0]), 100.0)
    assert shut[2] == ('closed', 1)


def test_a_year_under_the_threshold_makes_a_river_intermittent():
    flows = np.array([[150.0, 90.0, 0.5], [130.0, 120.0, 40.0], [170.0, 110.0, 60.0]])
    c = SR.classify(flows, 100.0)
    np.testing.assert_allclose(c['lowest'], [130.0, 90.0, 0.5])
    np.testing.assert_allclose(c['second'], [150.0, 110.0, 40.0])
    assert c['years'].tolist() == [3, 2, 0]
    assert c['sustained'].tolist() == [True, False, False]
    assert c['stops'].tolist() == [False, False, True]


def test_a_lake_across_the_seam_is_one_body():
    rows, cols = 4, 8
    wet = np.zeros(rows * cols, bool)
    wet[[1 * cols + 0, 1 * cols + 7, 2 * cols + 7, 3 * cols + 3]] = True
    body, km2 = SR.lake_bodies(wet, rows, cols, np.full(rows * cols, 2e6))
    assert body[1 * cols + 0] == body[1 * cols + 7] == body[2 * cols + 7] != body[3 * cols + 3]
    assert km2[body[1 * cols + 0]] == pytest.approx(6.0)


# ---------------------------------------------------------------- the product

@pytest.fixture(scope='module')
def record():
    if not SR.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.sustained_rivers first')
    return json.loads(SR.OUT.read_text())


def test_the_years_procedure_reproduces_the_design_routing(record):
    check = record['design_check']
    assert check['replay_largest_difference_m3s'] == 0.0
    assert check['years_procedure_largest_difference_m3s'] < 0.01
    main = check['main_part']
    assert main['routing_balance_m3s'] - main['design_spill_m3s'] == pytest.approx(main['offset_held_m3s'], abs=1.0)
    spread = record['against_the_design']['years_mean_over_design']
    assert 0.9 < spread['median'] < 1.1
    for width in ('narrow', 'wide'):
        steady = record['storage'][width]
        assert steady['largest_change_of_a_lake_outflow_m3s'] < 0.01
        for design, held in steady['outflows_m3s'].values():
            assert held == pytest.approx(design, abs=0.1)


def test_the_classes_add_up(record):
    c = record['classes']
    for k, v in c['by_storage'].items():
        assert v['rivers_sustained_to_their_end'] + v['rivers_intermittent_at_their_end'] == c['rivers']
        assert v['km_sustained'] + v['km_intermittent'] == pytest.approx(c['major_km'], abs=2)
    for row in c['largest_sustained'] + c['intermittent_largest']:
        for k in ('none', 'narrow', 'wide'):
            assert row[k]['lowest_m3s'] <= row[k]['second_lowest_m3s'] <= max(row[k]['mean_m3s'], row[k]['second_lowest_m3s'])
        assert row['after_losses_m3s']['high'] <= row['after_losses_m3s']['low'] + 0.1 <= row['design_m3s'] + 1.0
    for row in c['largest_sustained']:
        assert row['wide']['sustained'] and row['wide']['lowest_m3s'] >= c['threshold_m3s']


def test_the_losses_reach_the_seas(record):
    routed = record['losses']['routed']
    assert routed['rivers_into_the_seas_m3s']['high'] <= routed['rivers_into_the_seas_m3s']['low'] <= \
        routed['rivers_into_the_seas_m3s']['none']
    for budget in record['losses']['budgets'].values():
        narrow, wide = (budget['channel'][f'{c:g}']['evaporation_m3s'] for c in SR.RL.WIDTH)
        assert 0 < narrow < wide


def test_the_product_matches_its_record(record):
    if not SR.PRODUCT.exists():
        pytest.skip('the product is on the research drive')
    assert SR.sha256(SR.PRODUCT) == record['product']['sha256']
    assert SR.sha256(SR.REACHES) == record['reaches']['sha256']
    with np.load(SR.PRODUCT, allow_pickle=False) as z:
        meta = json.loads(str(z['metadata']))
        flows, low, major, inter, reach = z['flow_m3s'], z['sustained_m3s'], z['major'], z['intermittent'], z['reach']
    assert meta['schema'] == SR.PRODUCT_SCHEMA and flows.shape[0] == len(meta['years'])
    np.testing.assert_allclose(low, flows.min(axis=0))
    assert np.array_equal(inter, major & (low < meta['threshold_m3s']))
    assert (reach[major] >= 0).all() and (reach[~major] < 0).all()
