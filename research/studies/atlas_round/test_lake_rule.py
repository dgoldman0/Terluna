"""Checks of the corrected lake rule (lake_rule.py) on small grids: a hollow of three lowest points fills step by step as
worked by hand; on random ground no lake stands above dry ground beside it, the budget closes, and each lake's own
balance agrees with what the hollows' tree gave it; where every hollow has one lowest point the rule is route_levels'."""

import json

import numpy as np
import pytest
from scipy import ndimage

from research.studies.atlas_round import lake_rule as LR
from research.studies.atlas_round import lake_support as LS
from research.studies.drainage_network.routing import route_levels

U = LR.MM_DAY * 1e6                  # m3/s of 1 mm/day over a cell of 1 km2: the unit of the hand-worked cases


def profile(q, runoff_at=None, floor_runoff=True):
    """One row of ground round the Moon: a sea, then a hollow with three lowest points A, B and C behind sills at 6
    and 9 and an outer pass at 10 on A's side; every cell under water loses 1 unit, and q units are put in at A."""
    h = np.array([[-5.0, 10.0, 2.0, 6.0, 1.0, 9.0, 3.0, 12.0]], np.float32)
    outlet = np.zeros(h.shape, bool)
    outlet[0, 0] = True
    runoff = np.zeros(h.shape, np.float32)
    if runoff_at is not None:
        runoff[0, runoff_at[0]] = runoff_at[1]
    net = np.full(h.shape, -1.0, np.float32)
    return LR.route_hollows(h, outlet, np.zeros(1), 1.0, runoff, net, np.zeros(h.shape), sources=[(2, q * U)] if q else (),
                            floor_runoff=floor_runoff)


def state(r):
    wet = np.flatnonzero(r['lake']).tolist()
    levels = sorted(set(np.round(r['lake_level'][r['lake']], 3).tolist()))
    return wet, levels, round(float(r['discharge'][0]) / U, 6)


def test_a_hollow_of_three_lowest_points_fills_step_by_step():
    assert state(profile(0.5)) == ([], [], 0.0)                     # less than A's lowest cell loses: no lake
    assert state(profile(1.5)) == ([2], [6.0], 0.0)                 # A stands at its sill; the rest is lost in B
    assert state(profile(2.5)) == ([2, 4], [6.0], 0.0)              # A and B both at the sill between them
    assert state(profile(3.5)) == ([2, 3, 4], [9.0], 0.0)           # one lake over the sill, full to the next
    assert state(profile(4.5)) == ([2, 3, 4, 6], [9.0], 0.0)        # C filled from it to the sill between them
    assert state(profile(6.0)) == ([2, 3, 4, 5, 6], [10.0], 1.0)    # the whole hollow full, 1 unit over the outer pass
    r = profile(3.5)
    assert r['lake_spills'][[2, 3, 4]].all() and r['discharge'][6] / U == pytest.approx(0.5)   # A and B's lake spills into C
    r = profile(2.5)
    assert not r['lake_spills'][[2, 4]].any()                       # joined at the sill, a lake that passes nothing on


def test_water_is_followed_to_where_it_ends():
    r = profile(3.5)
    lakes, end = LR.trace(r, 2)
    assert end == ('dry', 6)                                        # over the sill into C, which holds no lake
    assert [t for t, _ in lakes] == [r['body'][2]] and lakes[0][1] / U == pytest.approx(0.5)
    r = profile(4.5)
    lakes, end = LR.trace(r, 3)                                     # one lake of all three, passing nothing on
    assert end == ('lake', r['body'][6]) and r['body'][6] == r['body'][2] and lakes[-1][1] == 0.0
    r = profile(6.0)
    lakes, end = LR.trace(r, 4)                                     # over the outer pass and down to the sea
    assert end == ('sea', 0) and lakes[0][1] / U == pytest.approx(1.0)


def test_dry_ground_below_a_pass_yields_runoff_unless_asked_not_to():
    wet, levels, sea = state(profile(0.0, runoff_at=(5, 2.0)))      # the sill between B and C, in B's catchment
    assert (wet, levels, sea) == ([2, 4], [6.0], 0.0)               # B fills and spills into A, which fills to the sill
    assert state(profile(0.0, runoff_at=(5, 2.0), floor_runoff=False)) == ([], [], 0.0)


def random_case(seed, rows=60, cols=96):
    """Broad basins pocked with smaller hollows, a sea over the lowest 15%, a climate in which open water mostly loses
    (here and there it gains), and four flows put in on the ground."""
    rng = np.random.default_rng(seed)
    broad = ndimage.gaussian_filter(rng.normal(size=(rows, cols)), 6.0, mode=('nearest', 'wrap'))
    fine = ndimage.gaussian_filter(rng.normal(size=(rows, cols)), 1.0, mode=('nearest', 'wrap'))
    h = (1000.0 * broad / broad.std() + 400.0 * fine / fine.std()).astype(np.float32)
    level = float(np.percentile(h, 15))
    outlet = h < level
    lat = np.linspace(55.0, -55.0, rows)
    runoff = (rng.gamma(0.6, 1.5, size=(rows, cols)) * (rng.random((rows, cols)) < 0.7)).astype(np.float32)
    net = (rng.normal(-3.0, 2.0, size=(rows, cols))).astype(np.float32)
    ground = np.flatnonzero(~outlet.ravel())
    sources = [(int(c), q * U) for c, q in zip(rng.choice(ground, 4, replace=False), (300.0, 120.0, 60.0, 25.0))]
    return h, outlet, lat, runoff, net, np.where(outlet, level, 0.0), sources


@pytest.mark.parametrize('seed', [1, 2, 3, 4, 5, 6])
@pytest.mark.parametrize('floor_runoff', [True, False])
def test_on_random_ground_shores_hold_and_the_budget_closes(seed, floor_runoff):
    h, outlet, lat, runoff, net, level, sources = random_case(seed)
    rows, cols = h.shape
    r = LR.route_hollows(h, outlet, lat, 1.0, runoff, net, level, sources=sources, floor_runoff=floor_runoff)
    lake = r['lake'].reshape(rows, cols)
    tree, full = r['tree'], r['full']
    n_all = full.size
    sills = full & (tree['parent'] >= 0) & ~full[np.maximum(tree['sibling'], 0)]
    assert lake.sum() > 100 and r['lake_spills'].any() and (lake & ~r['lake_spills'].reshape(rows, cols)).any()
    if floor_runoff:                                                # the ground exercises the rule's every part
        assert full[tree['n_nodes']:].any()                         # hollows of several catchments standing full
        assert (r['over'][sills] > 0).any()                         # a full hollow spilling into one that is not
    # no lake stands above dry ground beside it
    water_level = np.where(lake, r['lake_level'].reshape(rows, cols), np.nan).astype(np.float64)
    water_level[outlet] = level[outlet]
    assert not LS.unheld(h.astype(np.float64), water_level, lake, drop=1e-6).any()
    # a body of touching lake cells stands at one level
    bodies, count = ndimage.label(lake, structure=np.ones((3, 3), bool))
    wrap = bodies[:, 0][(bodies[:, 0] > 0) & (bodies[:, -1] > 0)], bodies[:, -1][(bodies[:, 0] > 0) & (bodies[:, -1] > 0)]
    top = ndimage.maximum(water_level, bodies, np.arange(1, count + 1))
    low = ndimage.minimum(water_level, bodies, np.arange(1, count + 1))
    assert np.abs(top - low).max() < 1e-6
    assert np.abs(top[wrap[0] - 1] - top[wrap[1] - 1]).max(initial=0.0) < 1e-6
    # the budget: dry ground's runoff, the put-in and the lakes' balance against the seas and what lakes not full keep
    b = r['budget']
    total = b['runoff_of_dry_ground'] + b['put_in'] + b['lakes_rain_less_evaporation']
    assert total == pytest.approx(b['reaches_the_seas'] + b['kept_by_lakes_not_full'], abs=1e-9 * max(abs(total), 1.0))
    assert (r['left'] >= -1e-12).all()
    # each lake's own balance, from the cells: what flows in over its shore and starts on it, less its loss, is what
    # the tree says it passes on or keeps
    flat = lake.ravel()
    donors = np.flatnonzero(~flat & ~r['outlet'] & (r['receiver'] != np.arange(flat.size)))
    into = donors[flat[r['receiver'][donors]]]
    gain = np.zeros(flat.size)
    np.add.at(gain, r['receiver'][into], r['discharge'][into])
    gain[flat] += r['own_flow'][flat] + r['net'][flat]
    dry_pits = r['sinks'][~flat[r['sinks']]]                        # a lowest point with no lake keeps what reaches it
    spare = gain[flat].sum() + r['discharge'][dry_pits].sum()
    leaves_out = full & (tree['parent'] < 0) & (np.arange(n_all) >= tree['n_out'])
    passed = r['out'][leaves_out].sum() + r['over'][sills].sum()
    assert spare == pytest.approx(passed + r['left'][r['made'] & ~full].sum(), abs=1e-9 * max(abs(spare), 1.0))
    # no lake loses more than it is given; two full hollows joined at their sill are one lake, though the sill's own
    # cell stays dry between them
    ids = r['body'][flat]
    by_body = np.bincount(np.unique(ids, return_inverse=True)[1], weights=gain[flat])
    assert by_body.min() > -1e-9 * max(abs(total), 1.0)
    first = np.unique(ids, return_index=True)[1]
    assert np.allclose(r['lake_level'][flat], r['lake_level'][flat][first][np.unique(ids, return_inverse=True)[1]])


def test_where_every_hollow_has_one_lowest_point_the_rule_is_route_levels():
    rows, cols = 30, 48
    y, x = np.mgrid[0:rows, 0:cols].astype(np.float64)
    h = 20.0 * x + 3.0 * np.abs(y - 15.0)                           # a slope down to the first columns
    for cy, cx, depth in ((8, 14, 90.0), (21, 14, 60.0), (8, 28, 120.0), (21, 30, 150.0), (14, 40, 100.0)):
        d2 = (y - cy) ** 2 + (x - cx) ** 2
        h -= depth * np.exp(-d2 / 6.0)                              # bowls well apart, each with one lowest point
    h = h.astype(np.float32)
    outlet = np.zeros(h.shape, bool)
    outlet[:, :2] = True                                            # the slope's top drains round the wrap to it
    h[outlet] = -50.0                                               # a level sea below all the ground
    lat = np.zeros(rows)
    rng = np.random.default_rng(11)
    runoff = rng.gamma(1.0, 2.0, size=h.shape).astype(np.float32)
    net = np.full(h.shape, -1.0, np.float32)
    net[:, 26:] = -30.0                                             # the upper bowls dry, the lower ones overflowing
    level = np.where(outlet, float(h[outlet].max()), 0.0)
    ref = route_levels(h, outlet, lat, 1.0, runoff, net, level)
    r = LR.route_hollows(h, outlet, lat, 1.0, runoff, net, level, floor_runoff=False)
    assert (r['tree']['parent'][r['tree']['n_out']:] < 0).all()     # no hollow lies inside another
    assert ref['lake'].sum() > 10 and ref['closed'].any() and (~ref['closed'][ref['outlet_nodes']:]).any()
    assert (r['lake'] == ref['lake']).all()
    assert np.allclose(r['lake_level'][r['lake']], ref['lake_level'][ref['lake']])
    assert r['discharge'][outlet.ravel()].sum() == pytest.approx(ref['discharge'][outlet.ravel()].sum())
    dry = ~ref['lake'] & ~outlet.ravel()
    assert np.allclose(r['discharge'][dry], ref['discharge'][dry])


def test_the_record_is_the_modules():
    if not LR.OUT.exists():
        pytest.skip('the lake rule has not been routed here')
    record = json.loads(LR.OUT.read_text())
    assert record['schema'] == LR.SCHEMA
    assert record['producer']['files'] == {f: LR.sha256(LR.ROOT / f) for f in LR.CODE}
    assert record['control_largest_difference_m3s'] < 1e-3
    for key in ('corrected', 'corrected_without_floor_runoff'):
        case = record['cases'][key]
        assert case['lakes']['shore_cells_not_held'] == 0
        assert abs(case['budget_gap_m3s']) < 1.0
    assert record['cases']['present']['lakes']['shore_cells_not_held'] > 50_000


@pytest.mark.parametrize('seed', [1, 4])
def test_what_reaches_each_sea_is_what_the_rivers_gather(seed):
    h, outlet, lat, runoff, net, level, sources = random_case(seed)
    water = np.where(outlet.ravel(), (np.arange(h.size) % h.shape[1] > h.shape[1] // 2).astype(np.int8), -1).astype(np.int8)
    s = LR.hollow_structure(h, outlet, lat, 1.0, level)
    terms = LR.hollow_terms(s, runoff, net)
    state = LR.hollow_state(s, terms, sources)
    r = LR.hollow_routing(s, terms, state)
    gathered = [r['discharge'][water == i].sum() for i in (0, 1)]
    assert min(gathered) > 0.0                                      # both waters receive rivers
    assert LR.reaches(s, terms, state, water) == pytest.approx(gathered)
    whole = LR.route_hollows(h, outlet, lat, 1.0, runoff, net, level, sources=sources)
    assert np.array_equal(whole['discharge'], r['discharge']) and np.array_equal(whole['lake'], r['lake'])
