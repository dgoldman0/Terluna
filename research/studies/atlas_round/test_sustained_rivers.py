"""Checks of the sustained rivers: on small grids, the network on the corrected rule's lakes gives the joined courses'
flows back exactly, the catchments' cascade gives each water what the replay gives it, a flow put in is carried down
and a river's loss is taken where it falls, a lake the design leaves below its pass passes nothing, and the storage
holds still on a steady year, keeps its water, carries a lake's outflow through a dry year and raises a lake below its
pass in wet years; Moscoviense's put-ins follow its fixed shares and each cut its part of them; a river's stems follow
the larger branch, pass small lakes and end at large ones; and, once computed, the product and its record."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import hollow_courses as HC
from research.studies.atlas_round import lake_rule as LR
from research.studies.atlas_round import sustained_rivers as SR

CELL_KM = 1.9
YEAR_S = SR.YEAR_S
QUIET = lambda *a: None                                                         # noqa: E731


def basin(rows=24, cols=48, lake_loss=3.0):
    """A slope rising to the east, a bowl on it that fills to a lake and spills west, a sea along the west edge and the
    six other waters as small blocks on high ground, each at its own level; with runoff (mm/day) on the land and the
    open water losing `lake_loss` mm/day beyond its rain."""
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
    runoff = np.full((rows, cols), 1.0, np.float32)
    lake_net = np.full((rows, cols), -lake_loss, np.float32)
    return h.astype(np.float32), water, levels, lat, runoff, lake_net


def network(h, water, levels, lat, runoff, lake_net, sources=()):
    """The corrected rule's routing of a small floor with its courses joined, and the network on it."""
    outlet = (water >= 0).reshape(h.shape)
    level = np.where(water >= 0, np.asarray(levels)[np.maximum(water, 0)], 0.0).reshape(h.shape)
    hollows = LR.hollow_structure(h, outlet, lat, CELL_KM, level)
    terms = LR.hollow_terms(hollows, runoff, lake_net)
    state = LR.hollow_state(hollows, terms, sources)
    r = LR.hollow_routing(hollows, terms, state)
    c = HC.connect(hollows, terms, state, r, h.shape)
    return SR.Network(hollows, state, r, c, water, len(SR.WATERS), h.shape, log=QUIET), r, c


def structure(ops=()):
    """A design structure for the small grid, the near-side sea one sea across its neck: `ops` as (cell, kind, sign,
    part of the branch's share)."""
    meta = dict(fractions={'south_pole_aitken': 0.5, 'humboldtianum': 0.2, 'smythii_marginis': 0.3},
                neck_flow=0.0, joined=True, level0=0.0)
    return dict(meta=meta, op_cell=np.array([c for c, _, _, _ in ops], np.int64),
                op_kind=np.array([SR.OPS.index(k) for _, k, _, _ in ops], np.int8),
                op_sign=np.array([x for _, _, x, _ in ops], float), op_part=np.array([p for _, _, _, p in ops], float))


@pytest.fixture(scope='module')
def small():
    h, water, levels, lat, runoff, lake_net = basin()
    net, r, c = network(h, water, levels, lat, runoff, lake_net)
    return net, r, c, runoff, lake_net


def test_the_replay_is_the_joined_courses(small):
    net, r, c, runoff, lake_net = small
    src = net.sources(runoff, lake_net)
    q, out = net.replay(src, passed_on=net.design_out)          # the design's own lake outflows: the routing itself
    assert np.array_equal(q, c['discharge'])
    assert np.array_equal(net.receiver, c['receiver'])
    q_y, out_y = net.replay(src)                                # the years' procedure: each lake's balance taken again
    np.testing.assert_allclose(q_y, c['discharge'], rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(out_y, net.design_out, rtol=1e-9, atol=1e-9)
    assert (r['lake'] & ~r['outlet']).sum() > 20                # the bowl holds a lake
    assert net.design_out.max() > 0.0 and not net.closed.any()  # and spills


def test_the_cascade_gives_each_water_what_the_replay_does(small):
    net, r, c, runoff, lake_net = small
    src = net.sources(runoff, lake_net)
    q, out = net.replay(src)
    _, out_c, supply = net.cascade(src)
    np.testing.assert_allclose(out_c, out, rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(supply, net.supply_of(q), rtol=1e-12, atol=1e-9)
    assert supply[3] > 0.0
    q2, _ = net.replay(src, passed_on=out)                      # each lake passing on what it passes on anyway
    np.testing.assert_allclose(q2, q, rtol=1e-12, atol=1e-9)


def upstream_of_the_lake(net, k):
    """A cell on dry ground whose water runs into the lake that leaves by catchment k."""
    above = int(np.flatnonzero(net.receiver == net.sinks[k - net.n_out])[0])
    while net.lake_full[above]:
        above = int(np.flatnonzero(net.receiver == above)[0])
    return above


def test_a_flow_put_in_is_carried_down_and_one_taken_out_of_a_lake_leaves_it(small):
    net, r, c, runoff, lake_net = small
    src = net.sources(runoff, lake_net)
    q, out = net.replay(src)
    _, _, supply = net.cascade(src)
    k = int(np.argmax(out))                                     # the bowl's lake
    above = upstream_of_the_lake(net, k)
    lake_cell = int(net.lake_cells[0])
    ops = [(above, 0.4 * out[k]), (lake_cell, -0.1 * out[k])]
    q_o, out_o = net.replay(src, ops)
    _, out_c, supply_o = net.cascade(src, ops)
    assert out_o[k] == pytest.approx(1.3 * out[k], rel=1e-9)
    np.testing.assert_allclose(out_c, out_o, rtol=1e-12, atol=1e-9)
    np.testing.assert_allclose(supply_o, net.supply_of(q_o), rtol=1e-12, atol=1e-9)
    assert supply_o[3] - supply[3] == pytest.approx(0.3 * out[k], rel=1e-9)
    assert q_o[above] == pytest.approx(q[above] + 0.4 * out[k], rel=1e-9)


def test_a_loss_is_taken_where_it_falls_and_never_more_than_the_river_carries(small):
    net, r, c, runoff, lake_net = small
    src = net.sources(runoff, lake_net)
    q, out = net.replay(src)
    k = int(np.argmax(out))                                     # the bowl's lake
    above = upstream_of_the_lake(net, k)
    loss = np.zeros(net.n)
    loss[above] = 0.25 * q[above]
    q_l, out_l = net.replay(src, loss=loss)
    assert out_l[k] == pytest.approx(out[k] - loss[above], rel=1e-9)
    assert net.supply_of(q_l)[3] == pytest.approx(net.supply_of(q)[3] - loss[above], rel=1e-9)
    loss[above] = 10.0 * q[above]                               # more than it carries: it runs dry there
    q_l, _ = net.replay(src, loss=loss)
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


def test_moscovienses_put_ins_follow_its_fixed_shares_and_each_cut_its_part():
    s = structure([(0, 'salt_sea_weir', 1.0, 1.0), (1, 'humboldtianum_branch', 1.0, 1.0),
                   (2, 'smythii_marginis_branch', 1.0, 1.0), (3, 'smythii_marginis_branch', -1.0, 0.9),
                   (4, 'smythii_marginis_branch', 1.0, 0.9), (5, 'orientale', 1.0, 1.0), (6, 'neck', 1.0, 1.0)])
    s['meta']['fractions'] = {'south_pole_aitken': 0.6, 'humboldtianum': 0.15, 'smythii_marginis': 0.25}
    q = SR.op_flows(s, 1000.0, 300.0, 400.0, 50.0)
    np.testing.assert_allclose(q, [600.0, 150.0, 250.0, -225.0, 225.0, 300.0, 50.0])
    assert SR.ops_list(s, q)[3] == (3, -225.0)


@pytest.fixture(scope='module')
def stored(small):
    net, r, c, runoff, lake_net = small
    s = structure()
    design = SR.without_storage(net, s, runoff, lake_net)
    design['inflow'] = net.cascade(design['sources'])[0]
    return net, s, design, runoff, lake_net


def test_the_storage_holds_still_on_the_design_year(stored):
    net, s, design, runoff, lake_net = stored
    st = SR.Storage(net, s, design, SR.RL.WIDTH[0], log=QUIET)
    st.initial(design, runoff, lake_net, design['inflow'])
    res = st.year(runoff, lake_net)
    lakes = st.lake_nodes
    np.testing.assert_allclose(res['out'][lakes], design['out'][lakes], rtol=1e-6, atol=1e-6)
    q = st.field(runoff, lake_net, res)
    np.testing.assert_allclose(q, design['q'], rtol=1e-6, atol=1e-6)


def test_the_storage_keeps_its_water_and_carries_a_lake_through_a_dry_year(stored):
    net, s, design, runoff, lake_net = stored
    st = SR.Storage(net, s, design, SR.RL.WIDTH[0], log=QUIET)
    st.initial(design, runoff, lake_net, design['inflow'])
    k = int(np.argmax(np.where(np.arange(net.n_nodes) >= net.n_out, design['out'], 0.0)))
    i = st.slot[k]
    wet = 1.5 * runoff
    src = net.sources(wet, lake_net)
    inflow = net.cascade(src)[0][k]
    v0 = st.vol[i]
    res = st.year(wet, lake_net)
    assert st.vol[i] > 0.0                                       # above its pass all year
    assert st.vol[i] - v0 == pytest.approx(YEAR_S * (inflow + src['lake'][k] - res['out'][k]), rel=1e-6)
    dry = 0.3 * runoff
    res = st.year(dry, lake_net)
    plain = SR.without_storage(net, s, dry, lake_net)
    assert res['out'][k] > plain['out'][k] + 0.5                 # the lake's head carries water into the dry year


def test_a_lake_below_its_pass_passes_nothing_and_rises_in_wet_years():
    h, water, levels, lat, runoff, lake_net = basin(lake_loss=60.0)
    net, r, c = network(h, water, levels, lat, runoff, lake_net)
    k = int(np.flatnonzero(net.closed)[0])                       # the bowl's lake, which its rivers cannot fill
    assert net.closed.sum() == 1 and 0 < net.lake_cells.size < net.cover_cells.size
    s = structure()
    design = SR.without_storage(net, s, runoff, lake_net)
    assert design['out'][k] == 0.0
    np.testing.assert_allclose(design['q'], c['discharge'], rtol=1e-9, atol=1e-9)
    wet = SR.without_storage(net, s, 20.0 * runoff, lake_net)    # without storage a wet year cannot raise it to its pass
    assert wet['out'][k] == 0.0
    design['inflow'] = net.cascade(design['sources'])[0]
    st = SR.Storage(net, s, design, SR.RL.WIDTH[0], log=QUIET)
    st.initial(design, runoff, lake_net, design['inflow'])
    i = st.slot[k]
    v0 = st.vol[i]
    assert -st.vol_full[i] < v0 < 0.0                            # it stands below its pass, not empty
    # full to its pass it would cover ground that now yields runoff, and lose more than reaches it
    cover = net.cover_cells[net.cover_node == k]
    dry = cover[~net.lake_full[cover]]
    full = (lake_net.ravel()[cover] * SR.MM_DAY * net.area[cover]).sum() - (runoff.ravel()[dry] * SR.MM_DAY * net.area[dry]).sum()
    table = st.table(runoff, lake_net)
    assert table[i, 0] == pytest.approx(full, rel=1e-9) and design['inflow'][k] + table[i, 0] < 0.0
    res = st.year(runoff, lake_net)
    assert res['out'][k] == 0.0 and st.vol[i] == pytest.approx(v0, rel=1e-6)
    res = st.year(3.0 * runoff, lake_net)
    assert st.vol[i] > v0 and res['out'][k] == 0.0               # a wet year raises it, still under its pass
    for _ in range(200):
        res = st.year(40.0 * runoff, lake_net)
    assert st.vol[i] > 0.0 and res['out'][k] > 0.0               # and wet years enough bring it over


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
    return (h.astype(np.float32), water, levels, lat, np.full((rows, cols), 2.0, np.float32),
            np.full((rows, cols), -1.0, np.float32))


def test_a_lake_of_two_hollows_stores_as_one_body():
    h, water, levels, lat, runoff, lake_net = twin_hollows()
    net, r, c = network(h, water, levels, lat, runoff, lake_net)
    inner = np.flatnonzero(net.through)
    assert inner.size >= 1                                       # the two hollows stand as one water
    k = int(net.body[inner[0]])
    assert not net.through[k] and net.level[inner[0]] == net.level[k]
    s = structure()
    design = SR.without_storage(net, s, runoff, lake_net)
    np.testing.assert_allclose(design['q'], c['discharge'], rtol=1e-9, atol=1e-9)
    design['inflow'] = net.cascade(design['sources'])[0]
    st = SR.Storage(net, s, design, SR.RL.WIDTH[0], log=QUIET)
    own = net.area[net.lake_cells[net.label[net.lake_cells] == k]].sum()
    whole = net.area[net.lake_cells[net.lake_node == k]].sum()
    assert st.slot[inner].max() < 0                              # the lowest points inside hold no volume of their own
    assert st.area_head[st.slot[k]] == pytest.approx(whole) and whole > own > 0.0
    st.initial(design, runoff, lake_net, design['inflow'])
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
    assert check['routing_made_again_largest_difference_m3s'] == 0.0
    assert check['replay_largest_difference_m3s'] == 0.0
    settled = SR.GD.SETTLE_LIMIT_M                              # the design's weirs are set within this of each balance
    assert check['years_procedure_largest_difference_m3s'] < settled
    for designed, routed in check['years_procedure_outflows_m3s'].values():
        assert routed == pytest.approx(designed, abs=settled)
    near = check['near_side']                                   # each part settled within 0.2% of what it receives
    assert abs(near['main']['left_m3s']) < 0.003 * near['main']['rivers_m3s']
    assert abs(near['south']['left_m3s']) < 0.003 * (near['south']['rivers_m3s'] + near['neck_flow_m3s'])
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
