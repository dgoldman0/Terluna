"""Tests of the rock's placement on the map and of the ground's evolution: the sums, the kernels and the products."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from research.studies.rock_return import evolve as EV
from research.studies.rock_return import placement as PL

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---- placement's sums --------------------------------------------------------------------------------------------
def test_a_range_is_more_columns_toward_the_poles():
    lat = np.array([0.0, 60.0, 89.9])
    rows, wide = PL.half_widths(20.0, lat, 1.9, 5760)
    assert rows == 11 and wide[0] == 11 and wide[1] == 21 and wide[2] > 1000


def test_sums_along_a_row_wrap_round_the_moon_and_stop_at_the_poles():
    field = np.zeros((3, 10))
    field[1, 0] = 1.0
    wide = np.array([1, 2, 1])
    s = PL.lon_sum(field, wide)
    assert s[1].tolist() == [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0] and s[0].sum() == 0.0
    ones = np.ones((4, 3))
    assert PL.lat_sum(ones, 1)[:, 0].tolist() == [2.0, 3.0, 3.0, 2.0]


def test_gathering_and_spreading_pair_the_same_cells():
    rng = np.random.default_rng(3)
    w, per = rng.random((12, 20)), rng.random((12, 20))
    wide = np.array([2, 2, 3, 3, 4, 4, 4, 4, 3, 3, 2, 2])
    assert math.isclose((per * PL.gather(w, 2, wide)).sum(), (w * PL.spread(per, 2, wide)).sum(), rel_tol=1e-12)


def test_rock_set_down_is_conserved_and_follows_the_weights():
    lying = np.zeros((9, 16))
    lying[4, 5] = 10.0
    w = np.zeros((9, 16))
    w[4, 4], w[4, 6], w[0, 0] = 1.0, 3.0, 5.0                       # two cells in range, one far off
    back, lost = PL.set_down(lying, w, 1, np.ones(9, np.int64))
    assert lost == 0.0 and math.isclose(back[4, 4], 2.5) and math.isclose(back[4, 6], 7.5) and back[0, 0] == 0.0
    back, lost = PL.set_down(lying, np.zeros((9, 16)), 1, np.ones(9, np.int64))
    assert lost == 10.0 and back.sum() == 0.0


def test_blocks_and_shares_of_land():
    assert PL.block_sum(np.ones((5, 9)), 2).tolist() == [[4.0] * 4, [4.0] * 4]
    p = PL.percentiles(np.array([1.0, 2.0, 3.0, 4.0]), np.array([1.0, 1.0, 1.0, 1.0]), (50,))
    assert 1.9 < p['50%'] < 2.1


# ---- the kernels of the ground's evolution -----------------------------------------------------------------------
def bowl(n=11):
    y, x = np.mgrid[0:n, 0:n]
    r = np.hypot(y - n // 2, x - n // 2)
    ground = 100.0 + 0.5 * r                                         # a plain rising gently away from the middle
    ground[r < 3.5] = 60.0 + 5.0 * r[r < 3.5]                        # with a hollow cut in it
    outlet = np.zeros((n, n), bool)
    outlet[0, :] = outlet[-1, :] = outlet[:, 0] = outlet[:, -1] = True
    return ground.ravel().astype(np.float64), outlet.ravel(), n


def test_the_flood_fills_a_hollow_to_its_lip_and_every_cell_drains_out():
    h, outlet, n = bowl()
    filled, receiver, order = EV.flood(h, outlet, n, n, 1e-6)
    assert order.size == h.size and (filled >= h).all()
    inside = filled - h > 1e-3
    assert inside.sum() >= 20 and np.ptp(filled[inside]) < 1e-3     # one level, the lip's
    assert (np.diff(filled[order]) >= 0).all()
    for i in np.flatnonzero(~outlet):                                # following receivers reaches an outlet
        j, steps = i, 0
        while not outlet[j]:
            j, steps = receiver[j], steps + 1
            assert steps < h.size


def test_the_steepest_way_and_the_land_drained():
    h = np.array([[9.0, 9.0, 9.0], [9.0, 5.0, 9.0], [9.0, 1.0, 9.0]]).ravel()
    outlet = np.zeros(9, bool)
    outlet[7] = True
    receiver = EV.steepest(h, outlet, 3, 3, 1.0, 1.0)
    assert receiver[4] == 7 and receiver[1] == 4 and receiver[7] == 7
    order = np.argsort(h, kind='stable')
    drained = EV.accumulate(receiver, order, 2.0)
    assert drained[7] == 18.0 and drained[4] >= 4.0
    assert math.isclose(EV.lengths(receiver, 3, 1.0, 2.0)[0], math.hypot(2.0, 1.0))      # cell 0 falls to the centre


def test_cutting_is_implicit_in_the_cell_below_and_stops_at_water():
    h = np.array([10.0, 20.0, 40.0])
    receiver, order = np.array([0, 0, 1]), np.array([0, 1, 2])
    drained, length = np.array([3.0, 2.0, 1.0]), np.ones(3)
    dry = np.zeros(3, bool)
    new = EV.erode(h, receiver, order, drained, length, dry, h, 0.5, 1.0, 1.0)
    assert new[0] == 10.0 and math.isclose(new[1], (20.0 + 1.0 * 10.0) / 2.0) and math.isclose(new[2], (40.0 + 0.5 * new[1]) / 1.5)
    wet = np.array([False, True, False])
    level = np.array([10.0, 30.0, 40.0])
    new = EV.erode(h, receiver, order, drained, length, wet, level, 0.5, 1.0, 1.0)
    assert new[1] == 20.0 and math.isclose(new[2], (40.0 + 0.5 * 30.0) / 1.5)             # cut to the water's surface


def test_cut_rock_comes_to_rest_in_the_first_hollow_or_at_the_outlet():
    receiver, order = np.array([0, 0, 1, 2, 4]), np.array([0, 1, 2, 3, 4])
    in_hollow = np.array([False, False, True, False, False])
    outlet = np.array([True, False, False, False, False])
    assert EV.entries(receiver, order, in_hollow, outlet).tolist() == [0, 0, 2, 2, -1]


def test_a_hollow_fills_from_the_bottom_to_one_level_and_no_higher_than_its_lip():
    h = np.array([0.0, 1.0, 2.0, 9.0])
    filled = np.array([3.0, 3.0, 3.0, 9.0])
    cells, starts = np.array([0, 1, 2]), np.array([0, 3])
    over = EV.fill_hollows(h, filled, cells, starts, np.array([2.0]), 1.0)
    assert over == 0.0 and np.allclose(h[:3], [1.5, 1.5, 2.0])
    h = np.array([0.0, 1.0, 2.0, 9.0])
    over = EV.fill_hollows(h, filled, cells, starts, np.array([10.0]), 1.0)
    assert math.isclose(over, 4.0) and np.allclose(h[:3], 3.0)


def test_carriers_share_rock_over_the_land_round_it():
    g = dict(shape=(21, 21), north_south=1900.0, east_west=1900.0)
    lying = np.zeros(441)
    lying[220] = 8.0
    weight = np.zeros(441)
    weight[219], weight[221] = 1.0, 3.0                              # two cells as near as each other
    back, unplaced = EV.set_down(lying, weight, g, range_m=4000.0)
    assert math.isclose(back[219], 2.0) and math.isclose(back[221], 6.0) and unplaced.sum() == 0.0
    far = np.zeros(441)
    far[0] = 5.0                                                     # the only weighted land lies out of range
    back, unplaced = EV.set_down(lying, far, g, range_m=4000.0)
    assert unplaced[220] == 8.0 and back.sum() == 0.0
    spread = EV.near(lying, g, 4000.0)
    assert math.isclose(spread.sum(), 8.0, rel_tol=1e-6) and spread[220] == spread.max()


def small_region(n=14, seed=5):
    rng = np.random.default_rng(seed)
    y, x = np.mgrid[0:n, 0:n]
    r = np.hypot(y - n / 2, x - n / 2)
    ground = 500.0 - 20.0 * r + 150.0 * np.exp(-((r - 4.0) / 1.2) ** 2) - 200.0 * np.exp(-(r / 2.5) ** 2) + rng.normal(0.0, 1.0, (n, n))
    return dict(ground=ground.astype(np.float32), water=np.full((n, n), np.nan, np.float32), sea=np.zeros((n, n), bool),
                eroding=np.ones((n, n), bool), lowering=np.full((n, n), 5.0e-5, np.float32), east_west_m=1900.0, north_south_m=1900.0)


@pytest.mark.parametrize('water', EV.WATERS)
@pytest.mark.parametrize('rule', EV.RULES)
def test_a_step_keeps_the_rock(water, rule):
    region = small_region()
    g, h = EV.grid_of(region), EV.start_of(region)
    total = dict(cut=0.0, laid=0.0, set_down=0.0, gone=0.0, slid=0.0)
    new, pending = h, None
    for _ in range(4):
        new, moved, pending = EV.step(new, g, water, rule, 2.0e-6, 25000.0, pending=pending)
        for k in total:
            total[k] += moved[k]
    assert total['cut'] > 0
    # the ground's change, the rock held where it slid onto water and the rock gone sum to nothing
    assert abs((new - h).sum() * g['cell_area'] + pending.sum() + total['gone']) < 1e-9 * total['cut']


def test_roundness_tells_a_disc_from_an_ellipse():
    y, x = np.mgrid[0:61, 0:61]
    disc = np.hypot(y - 30, x - 30) < 20
    ellipse = ((y - 30) / 10.0) ** 2 + ((x - 30) / 25.0) ** 2 < 1.0
    assert EV.roundness(disc, 1.0, 1.0) < 0.03 < 0.2 < EV.roundness(ellipse, 1.0, 1.0)
    assert EV.roundness(disc[:3, :3], 1.0, 1.0) is None
    assert EV.roundness(disc, 2.0, 1.0) > 0.2                         # cells twice as wide make it an ellipse on the ground


def test_the_rate_is_scaled_to_the_first_lowering():
    region = small_region(n=2 * EV.MARGIN + 8)                       # large enough to leave ground inside the margin
    for water in EV.WATERS:
        rate, target = EV.calibrate(region, water)
        assert rate > 0 and math.isclose(target, 5.0e-5, rel_tol=1e-6)


# ---- the products ------------------------------------------------------------------------------------------------
@pytest.fixture(scope='module')
def placement():
    return json.loads(PL.OUT.read_text())


@pytest.fixture(scope='module')
def evolution():
    return json.loads(EV.OUT.read_text())


def test_the_placement_product_is_whole(placement):
    assert placement['schema'] == PL.SCHEMA
    assert placement['producer']['code_sha256'] == {'placement.py': sha256(HERE / 'placement.py'), 'rock.py': sha256(HERE / 'rock.py')}
    for rel, digest in placement['inputs'].items():
        assert sha256(ROOT / rel) == digest, rel
    for rule in PL.RULES:
        for row in placement['rules'][rule].values():
            assert 0.999 < row['returned_share'] + row['rock_with_no_land_in_range'] < 1.001
            blocks = [row['set_down_on_the_wrong_cell']] + [row['wrong_at_blocks'][f'{b} cells'] for b in PL.BLOCKS]
            assert all(a >= b for a, b in zip(blocks, blocks[1:]))                   # misplacing shrinks as the map is taken in larger blocks


def test_what_the_write_up_says_of_placement(placement):
    rules = placement['rules']
    worn = [rules['worn'][f'{r:g} km']['set_down_on_the_wrong_cell'] for r in PL.RANGES_KM]
    assert worn[0] > worn[1] > worn[2] and worn[2] < 0.06 < 0.25 < worn[0]
    for r in PL.RANGES_KM:
        key = f'{r:g} km'
        assert rules['worn'][key]['set_down_on_the_wrong_cell'] < rules['even'][key]['set_down_on_the_wrong_cell'] < rules['crest'][key]['set_down_on_the_wrong_cell']
        assert rules['crest'][key]['set_down_on_the_wrong_cell'] > 0.5
    own = [rules['worn'][f'{r:g} km']['set_down_in_the_catchment_it_left'] for r in PL.RANGES_KM]
    assert own[0] > own[1] > own[2] and own[2] < 0.05
    assert all(rules[rule]['50 km']['wrong_at_blocks']['128 cells'] < 0.02 for rule in PL.RULES)


def test_the_evolution_product_is_whole(evolution, placement):
    assert evolution['schema'] == EV.SCHEMA
    assert evolution['producer']['code_sha256'] == sha256(HERE / 'evolve.py') and evolution['producer']['compiled']
    assert evolution['inputs']['research/studies/rock_return/results/placement.json'] == sha256(PL.OUT)
    assert evolution['ignored_inputs'][str(EV.REGIONS)] == placement['product']['sha256']
    for name, by_water in evolution['regions'].items():
        for water in EV.WATERS:
            for case in EV.CASES:
                row = by_water[water][case]
                assert set(row['years']) == {f'{y:g}' for y in EV.YEARS}
                for at in row['years'].values():
                    rock = at['rock']
                    assert abs(rock['ground'] + rock['held'] + rock['gone']) < 1e-4 * rock['cut']       # the floors' roughness is the remainder
                    assert -1.0 <= at['forms_kept'] <= 1.0


# ---- the return as designed --------------------------------------------------------------------------------------
def test_healing_gives_no_ground_more_than_it_lost_and_keeps_the_rock():
    g = dict(shape=(21, 21), north_south=1900.0, east_west=1900.0)
    rng = np.random.default_rng(2)
    lost = rng.random(441) * 5.0
    lying = np.zeros(441)
    lying[[50, 220, 400]] = [300.0, 40.0, 5.0]
    back, left = EV.heal(lying, lost, g, range_m=8000.0)
    assert math.isclose(back.sum() + left.sum(), lying.sum(), rel_tol=1e-12)
    assert (back <= lost + 1e-12).all() and left[50] > 0.0                      # the large heap has rock to spare
    plenty = np.full(441, 1000.0)
    back, left = EV.heal(plenty, lost, g, range_m=8000.0)
    assert np.allclose(back, lost)                                               # with rock to spare everywhere, all is made good
    back, left = EV.heal(np.zeros(441), lost, g)
    assert back.sum() == 0.0 and left.sum() == 0.0


def test_ground_too_steep_fails_until_it_is_not_and_keeps_its_rock():
    n = 9
    h = np.zeros((n, n))
    h[4, 4] = 5000.0                                                  # a needle on a plain of 1.9 km cells
    h = h.ravel()
    fixed = np.zeros(n * n, bool)
    hold = np.zeros(n * n, bool)
    pending = np.zeros(n * n)
    tan = math.tan(math.radians(35.0))
    before = h.sum()
    moved = EV.slide(h, fixed, hold, pending, tan, n, n, 1900.0, 1900.0, 1.0, 200)
    assert moved > 0 and math.isclose(h.sum(), before, rel_tol=1e-12) and pending.sum() == 0.0
    z = h.reshape(n, n)
    assert max(np.abs(np.diff(z, axis=0)).max(), np.abs(np.diff(z, axis=1)).max()) <= tan * 1900.0 * (1 + 1e-6)
    assert z[4, 4] < 5000.0 and z[4, 4] == z.max()
    # onto water or an outlet the rock does not pile: it is held as come to rest there
    h = np.zeros(n * n)
    h[40] = 5000.0
    hold[41] = True
    EV.slide(h, fixed, hold, pending, tan, n, n, 1900.0, 1900.0, 1.0, 200)
    assert h[41] == 0.0 and pending[41] > 0 and math.isclose(h.sum() + pending.sum(), 5000.0, rel_tol=1e-12)


def test_creep_smooths_and_keeps_the_rock():
    g = dict(shape=(5, 5), north_south=100.0, east_west=100.0, kappa=0.01)
    h = np.zeros(25)
    h[12] = 10.0
    land = np.ones(25, bool)
    new = EV.creep(h, land, g, 1000.0)
    assert math.isclose(new.sum(), 10.0) and new[12] < 10.0 and new[11] > 0.0
    land[11] = False
    new = EV.creep(h, land, g, 1000.0)
    assert new[11] == 0.0 and math.isclose(new.sum(), 10.0)           # nothing creeps onto water
    assert (EV.creep(h, land, dict(g, kappa=0.0), 1000.0) == h).all()


@pytest.fixture(scope='module')
def process():
    from research.studies.rock_return import process as PR
    return PR, json.loads(PR.OUT.read_text())


def test_the_process_product_is_whole(process, placement):
    PR, product = process
    assert product['schema'] == PR.SCHEMA and product['producer']['compiled']
    assert product['producer']['code_sha256'] == {'process.py': sha256(HERE / 'process.py'), 'evolve.py': sha256(HERE / 'evolve.py')}
    assert product['ignored_inputs'][str(EV.REGIONS)] == placement['product']['sha256']
    for name, design in PR.DESIGN.items():
        for reading in PR.READINGS:
            for label in design['cases']:
                row = product['regions'][name][reading][label]
                assert set(row['years']) == {f'{y:g}' for y in PR.YEARS}
                limit = row['settings'].get('critical_deg', EV.CRITICAL_DEG)
                for at in row['years'].values():
                    rock = at['rock']
                    assert abs(rock['ground'] + rock['held'] + rock['gone']) < 1e-4 * rock['cut']
                    assert at['mean_slope_deg'] < limit
    for reading in PR.READINGS:
        assert set(product['checks']['grain'][reading]) == {'1', '2'} and set(product['checks']['step'][reading]) == {'half_the_step'}


def test_what_the_write_up_says_of_the_return_as_designed(process):
    PR, product = process
    regions = product['regions']
    for name in PR.DESIGN:
        shared, capped = regions[name]['shared_by_wear']['all_of_it'], regions[name]['no_more_than_it_lost']['all_of_it']
        # both readings keep the forms through thirty million years, and neither raises the map's own summits
        assert shared['years']['3e+07']['forms_kept'] > 0.9 and capped['years']['3e+07']['forms_kept'] > 0.9
        for row in (shared, capped):
            assert all(at['summits']['median_m'] < 0 for at in row['years'].values())
        # with ground that fails, nothing grows without limit
        assert all(at['heights']['highest'] < shared['start']['heights']['highest'] + 1500.0 for at in shared['years'].values())
        assert all(at['heights']['highest'] <= capped['start']['heights']['highest'] for at in capped['years'].values())
        # loose ground changes the forms kept by little
        for reading in PR.READINGS:
            a, b = regions[name][reading]['all_of_it'], regions[name][reading]['all_of_it_on_loose_ground']
            assert all(abs(a['years'][y]['forms_kept'] - b['years'][y]['forms_kept']) < 0.03 for y in a['years'])
    south, korolev = regions['southern_highlands'], regions['korolev']
    assert south['no_more_than_it_lost']['all_of_it']['years']['1e+08']['forms_kept'] > 0.9 > 0.65 > south['shared_by_wear']['all_of_it']['years']['1e+08']['forms_kept']
    assert south['shared_by_wear']['all_of_it']['years']['1e+08']['share_risen_over_300_m'] > 0.15
    assert korolev['shared_by_wear']['all_of_it']['years']['1e+08']['share_of_land_at_the_failing_slope'] > 0.05
    far = 'all_of_it_within_50_km'
    for reading in PR.READINGS:                                       # the longer range serves Korolev under either reading
        assert korolev[reading][far]['years']['1e+08']['forms_kept'] > korolev[reading]['all_of_it']['years']['1e+08']['forms_kept'] + 0.05
        assert korolev[reading]['flyers_and_waders_alone']['years']['3e+07']['forms_kept'] < 0.4
    assert korolev['no_more_than_it_lost'][far]['years']['1e+08']['forms_kept'] > 0.94
    # the model's own grain and step do not shape the answer
    for reading in PR.READINGS:
        grain = product['checks']['grain'][reading]
        assert abs(grain['1']['years']['3e+07']['forms_kept'] - grain['2']['years']['3e+07']['forms_kept']) < 0.02
        half = product['checks']['step'][reading]['half_the_step']['years']['3e+07']['forms_kept']
        assert abs(half - south[reading]['all_of_it']['years']['3e+07']['forms_kept']) < 0.002


# ---- hollows that fill and hollows that form ---------------------------------------------------------------------
def test_the_count_tells_old_hollows_from_new():
    from research.studies.rock_return import hollows as HO
    n = 30
    ground = np.full((n, n), 100.0)
    ground[5:10, 5:10] = 50.0                                        # one hollow of 25 cells at the start
    region = dict(ground=ground.astype(np.float32), water=np.full((n, n), np.nan, np.float32), sea=np.zeros((n, n), bool),
                  eroding=np.ones((n, n), bool), lowering=np.zeros((n, n), np.float32), east_west_m=1000.0, north_south_m=1000.0, margin=2)
    g, h0 = EV.grid_of(region), EV.start_of(region)
    inside = EV.interior(region).reshape(g['shape'])
    start, old = HO.census(h0, g, inside)
    assert start['number']['10-100 km2'] == 1 and math.isclose(start['volume_km3'], 25 * 50.0 / 1e3, rel_tol=1e-4) and old.sum() == 25
    later = h0.reshape(n, n).copy()
    later[5:10, 7] = 100.0                                           # rock set down across it cuts it in two
    later[20:24, 20:24] = 80.0                                       # and a new hollow of 16 cells stands elsewhere
    row, _ = HO.census(later.ravel(), g, inside, old)
    kinds = row['hollows_of_10_km2_or_more']
    assert kinds['on_old_ground']['number'] == 2 and kinds['on_new_ground']['number'] == 1 and kinds['on_both']['number'] == 0
    assert math.isclose(row['ground_in_a_hollow_now_and_not_at_the_start'], 16 / inside.sum())
    assert math.isclose(row['ground_in_a_hollow_at_the_start_and_not_now'], 5 / inside.sum())
    assert row['volume_km3'] < start['volume_km3'] + 16 * 20.0 / 1e3 + 1e-9


def test_the_hollows_product_and_what_the_write_up_says_of_it():
    from research.studies.rock_return import hollows as HO
    product = json.loads(HO.OUT.read_text())
    assert product['schema'] == HO.SCHEMA
    assert product['producer']['code_sha256'] == {name: sha256(HERE / name) for name in ('hollows.py', 'evolve.py', 'process.py')}
    assert product['inputs']['research/studies/rock_return/results/process.json'] == sha256(HERE / 'results' / 'process.json')
    for name, region in product['regions'].items():
        start = region['start']
        shared, capped = region['shared_by_wear']['1e+08'], region['no_more_than_it_lost']['1e+08']
        counted = lambda row: sum(v['number'] for v in row['hollows_of_10_km2_or_more'].values())      # noqa: E731
        # shared by wear the hollows grow in number and shrink in volume; new ones form, and they are small
        assert counted(shared) > start['hollows_of_10_km2_or_more']
        assert 0.7 < shared['volume_km3'] / start['volume_km3'] < 0.85
        new = shared['hollows_of_10_km2_or_more']['on_new_ground']
        assert new['number'] >= 10 and new['area_km2'] < 0.05 * sum(v['area_km2'] for v in shared['hollows_of_10_km2_or_more'].values())
        # taken no more than it lost, almost none form
        assert capped['hollows_of_10_km2_or_more']['on_new_ground']['number'] <= 10
        assert capped['hollows_of_10_km2_or_more']['on_new_ground']['area_km2'] < 0.1 * new['area_km2']
    korolev = product['regions']['korolev']
    assert korolev['no_more_than_it_lost']['1e+08']['volume_km3'] < 0.2 * korolev['start']['volume_km3']
    assert abs(korolev['shared_by_wear']['1e+08']['share_of_ground'] - korolev['start']['share_of_ground']) < 0.01
