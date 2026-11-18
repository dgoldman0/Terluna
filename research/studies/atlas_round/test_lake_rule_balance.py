"""Checks of the settings made again on the corrected lake rule (lake_rule_balance.py): on a row of seas worked by
hand the weirs give each held sea its need after the lakes on the way, a cut takes no more than its lake can spare,
and the near-side sea's parts stand where their losses meet their rivers, over the crest and under it."""

import json

import numpy as np
import pytest

from research.studies.atlas_round import lake_rule as LR
from research.studies.atlas_round import lake_rule_balance as LB
from research.studies.atlas_round import near_side_neck as N
from research.studies.drainage_network import network as NW

U = LR.MM_DAY * 1e6                  # m3/s of 1 mm/day over a cell of 1 km2: the unit of the hand-worked cases
SALT, HELD, REST, MAIN_SEA, SOURCE = 0, 9, 14, 17, 3


def row_of_seas(own=12.0):
    """One row of ground round the Moon. A source sea stands at 40 with a river of `own` units; the salt sea, the held
    sea, the rest's sea and the main sea stand at 0 and lose 4, 3, 0.5 and 1 unit. The source's pass opens on cell 1,
    which runs to the salt sea. The held sea's branch is put in at cell 6, a hollow that loses 1 unit before it
    overflows to the held sea. The rest's branch is put in at cell 11, a hollow that loses 2 units before it overflows
    to the rest's sea, whose own pass opens on cell 16, which runs to the main sea."""
    h = np.array([[-5.0, 5.0, 50.0, 35.0, 45.0, 48.0, 10.0, 15.0, 7.0, -5.0, 30.0, 8.0, 12.0, 6.0, -5.0, 50.0, 10.0, -5.0,
                   60.0, 20.0]], np.float32)
    water = np.full(h.size, -1, np.int8)
    names = list(LB.WATERS)
    for cell, name in ((SALT, LB.SALT_SEA), (HELD, LB.HELD), (REST, LB.REST), (MAIN_SEA, 'near_side_main'), (SOURCE, LB.SOURCE)):
        water[cell] = names.index(name)
    levels = [0.0] * len(names)
    levels[names.index(LB.SOURCE)] = 40.0
    runoff = np.zeros(h.shape, np.float32)
    runoff[0, 4] = own
    net = np.zeros(h.shape, np.float32)
    for cell, loss in ((SALT, 4.0), (HELD, 3.0), (REST, 0.5), (MAIN_SEA, 1.0), (SOURCE, 1.0), (6, 1.0), (11, 2.0)):
        net[0, cell] = -loss
    d = dict(h=h, lat=np.zeros(1), cell_km=1.0)
    s = LB.structure_of(d, water, levels)
    return s, LR.hollow_terms(s, runoff, net), water, names


def works(take_outs=False):
    rest = [(11, 1.0), (11, -1.0), (13, 1.0)] if take_outs else [(11, 1.0)]
    return dict(source_pass=1, full_pass={'orientale': 19, 'smythii_marginis': 16}, held_branch=[(6, 1.0)], rest_branch=rest)


@pytest.mark.parametrize('take_outs', [False, True])
def test_the_weirs_give_each_held_sea_its_need_after_the_lakes_on_the_way(take_outs):
    s, terms, water, names = row_of_seas()
    gov = LB.set_weirs(s, terms, water, names, works(take_outs), within=1e-9, log=lambda *a: None)
    at = names.index
    assert not gov['unmet']
    assert gov['surplus'] / U == pytest.approx(11.0)                # the source's river less its own loss
    share = {k: v / U for k, v in gov['share'].items()}
    assert share == pytest.approx({LB.SALT_SEA: 4.0, LB.HELD: 4.0, LB.REST: 3.0})    # the held branch's lake takes 1
    supply = gov['supply'] / U
    assert supply[at(LB.SALT_SEA)] == pytest.approx(4.0) and supply[at(LB.HELD)] == pytest.approx(3.0)
    assert supply[at(LB.REST)] == pytest.approx(1.0)                # the rest's branch after its lake has taken 2
    assert gov['over'][LB.REST] / U == pytest.approx(0.5)           # the rest's sea passes on what it has to spare
    assert supply[at('near_side_main')] == pytest.approx(0.5)
    r = LR.hollow_routing(s, terms, gov['state'])
    gathered = np.array([r['discharge'][water == i].sum() for i in range(len(names))])
    assert gathered == pytest.approx(gov['supply'])                 # the rivers gathered bring the seas the same
    assert r['lake'][[6, 11]].all() and r['lake_spills'][6] and float(r['left'].min()) >= 0.0
    if take_outs:
        assert [q / U for q in gov['flows']['rest_branch']] == pytest.approx([3.0, 1.0, 1.0], abs=1e-6)
        assert r['full'][LB.roots_of(r['tree'])[r['label'][11]]]    # the lake the cut leaves from stands full


def test_a_cut_takes_no_more_than_its_lake_can_spare():
    s, terms, _, _ = row_of_seas()
    root = LB.roots_of(s['tree'])
    branch = [(11, 1.0), (11, -1.0), (13, 1.0)]
    for share, beyond in ((5.0, 3.0), (2.5, 0.5), (1.0, 0.0)):       # the lake needs 2 units to stand full
        state, flows, sources = LB.branch_state(s, terms, [], [(branch, share * U)], root)
        assert [q / U for q in flows[0]] == pytest.approx([share, beyond, beyond], abs=1e-6)
        assert bool(state['full'][root[s['label'][11]]]) == (share >= 2.0)
        assert sum(q for _, q in sources) / U == pytest.approx(share, abs=1e-6)


def test_a_held_sea_that_cannot_be_met_takes_all_the_source_has_left():
    s, terms, water, names = row_of_seas(own=6.0)                   # a surplus of 5: the salt sea's 4, then 1 left
    gov = LB.set_weirs(s, terms, water, names, works(), within=1e-9, log=lambda *a: None)
    share = {k: v / U for k, v in gov['share'].items()}
    assert share == pytest.approx({LB.SALT_SEA: 4.0, LB.HELD: 1.0, LB.REST: 0.0})
    assert gov['unmet'][LB.HELD] / U == pytest.approx(3.0)          # its lake takes the 1 and the sea gets nothing
    assert gov['supply'][names.index(LB.HELD)] / U == pytest.approx(0.0, abs=1e-6)


def parts(main=100, south=50):
    """Two parts whose cells each join at their own height, a metre apart, and lose 1 m3/s each."""
    tab = {}
    for key, count in (('main', main), ('south', south)):
        joined = np.arange(count, dtype=float)
        tab[key] = NW.Table(joined, np.full(count, 1e6), np.ones(count), joined, 1000.0)
    return tab, N.Neck(60.0, [60.0, 61.0, 70.0], [0.0, 100.0, 1000.0], [0.0, 100.0, 100.0], 1000.0)


def test_under_the_crest_each_part_stands_on_its_own_rivers():
    tab, neck = parts()
    z = LB.part_levels(tab, dict(main=40.0, south=10.0), neck)
    assert not z['main_over_the_crest'] and z['neck_flow'] == 0.0
    assert (z['main'], z['south']) == (39.0, 9.0)                   # the cell whose loss brings each to its rivers


def test_over_the_crest_the_main_part_sends_the_south_what_the_neck_passes():
    tab, neck = parts()
    z = LB.part_levels(tab, dict(main=70.0, south=10.0), neck)
    assert z['main_over_the_crest'] and 0.0 < z['neck_flow'] <= 10.0 and z['main'] >= 60.0
    assert abs(tab['main'].at(z['main'])[1] - (70.0 - z['neck_flow'])) <= 1.0
    assert abs(tab['south'].at(z['south'])[1] - (10.0 + z['neck_flow'])) <= 1.0


def test_a_basin_that_joins_at_one_level_comes_under_the_water_whole():
    """Ten cells of a main part: three join at 0, 1 and 2 m and lose 1 each; four join together at 5 m, a basin beyond
    a sill, losing 5, gaining 4, losing 5 and gaining 5, so 1 in all; three more join at 6, 7 and 8 m and lose 1."""
    joined = np.array([0.0, 1.0, 2.0, 5.0, 5.0, 5.0, 5.0, 6.0, 7.0, 8.0, 3.0, 4.0], np.float32)
    part = np.array([1] * 10 + [2] * 2, np.int8)
    loss = np.array([1.0, 1.0, 1.0, 5.0, -4.0, 5.0, -5.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    d = dict(h=np.zeros((1, 12), np.float32), level0=0.0)
    st = dict(part=part, joined_main=joined, joined_south=joined)
    tab = LB.part_tables(d, st, loss, np.ones(12))
    assert N.level_for(tab['main'], 3.5, 100.0) == 5.0              # between the loss below the sill and with the basin
    assert N.level_for(tab['main'], 5.0, 100.0) == 6.0              # the basin whole loses 1: the water stands past it
    cell_by_cell = NW.Table(np.where(part == 1, joined, np.inf).astype(float), np.ones(12), loss, np.zeros(12), 1000.0)
    assert N.level_for(cell_by_cell, 5.0, 100.0) == 5.0             # the search cell by cell stops inside the basin
    assert tab['main'].at(6.0)[1] == 4.0 and tab['main'].at(5.0)[1] == 3.0


def test_the_record_is_the_modules():
    if not LB.OUT.exists():
        pytest.skip('the settings have not been made again here')
    record = json.loads(LB.OUT.read_text())
    assert record['schema'] == LB.SCHEMA
    assert record['producer']['files'] == {f: LB.sha256(LB.ROOT / f) for f in LB.CODE}
    seas = record['seas']['rebalanced']
    for held in (LB.SALT_SEA, LB.HELD):                             # the weirs hold them at the common level
        assert seas[held]['level_above_common_m'] == 0.0 and abs(seas[held]['balance_m3s']) <= 2
    for part in ('main', 'south'):                                  # each part stands where its loss meets its supply
        assert abs(record['near_side']['left_out_of_balance_m3s'][part]) < LB.PART_SETTLED_M3S
        assert abs(seas[f'near_side_{part}']['balance_m3s']) < LB.PART_SETTLED_M3S + 1
    assert record['lakes']['rebalanced']['shore_cells_not_held'] == 0
    assert abs(record['budget_gap_m3s']) < 1.0
    for cut in record['works']['take_outs']:
        assert cut['stands_full'] or cut['takes_m3s'] == 0.0
