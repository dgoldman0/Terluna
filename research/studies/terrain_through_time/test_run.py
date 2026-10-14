"""The terrain study: the routing summary's binding and balances, and the second stage's regeneration and physics."""
import hashlib
import json
import math
import os

import pytest

from research.studies.terrain_through_time import routing as rt
from research.studies.terrain_through_time import run as T
from shared.provenance import constants_changed


@pytest.fixture(scope='module')
def routing():
    return json.loads(rt.OUT.read_text())


@pytest.fixture(scope='module')
def product():
    return json.loads(T.OUT.read_text())


def digest(path):
    return hashlib.sha256((rt.ROOT / path).read_bytes()).hexdigest()


# ---- The routing summary ---------------------------------------------------------------------------------------------

def test_routing_binds_its_code_and_inputs(routing):
    assert routing['schema'] == rt.SCHEMA
    for path, h in {**routing['producer']['files'], **routing['producer']['inputs']}.items():
        assert digest(path) == h, path
    assert set(routing['producer']['files']) == set(rt.CODE)
    assert set(routing['producer']['inputs']) == set(rt.INPUTS)
    assert set(routing['producer']['ignored_inputs']) == set(rt.IGNORED)
    assert not constants_changed(routing['producer']['constants'])
    for path, h in routing['producer']['ignored_inputs'].items():
        if (rt.ROOT / path).exists():
            assert digest(path) == h, path


def test_routing_reproduced_the_drainage_product(routing):
    for key, (mine, theirs) in routing['reproduced_drainage'].items():
        assert mine == theirs, key


def test_rock_lands_somewhere_and_lakes_add_up(routing):
    for name, shares in routing['carrying_back']['first_settles_in'].items():
        assert math.isclose(sum(shares.values()), 1.0, rel_tol=1e-5), name
    assert math.isclose(sum(routing["upkeep_at_start"]["kept_share"].values()), 1.0, rel_tol=1e-5)
    assert math.isclose(sum(b['rock_share'] for b in routing['land']['by_band'].values()), 1.0, rel_tol=1e-6)
    u = routing['units']
    assert sum(c['lakes'] for c in u['by_class'].values()) == u['spilling_lakes'] + u['closed_lakes']
    assert math.isclose(sum(c['area_km2'] for c in u['by_class'].values()), u['lake_area_km2'], rel_tol=1e-5)
    assert math.isclose(u['lake_water_km3'], routing['reproduced_drainage']['lake_volume_km3'][0], rel_tol=1e-4)
    c = routing['carrying_back']
    assert c['mean_lift_to_heads_m'] > c['mean_lift_m'] > 0
    assert 0.3 < c['heads_share_of_land'] < 0.36


def test_lakes_only_close_and_fill_in_order(routing):
    years = [f'{y:g}' for y in rt.SURVIVAL_YEARS]
    for name, f in routing['fill'].items():
        surv = f['survival']
        for cls in surv.values():
            counts = [cls['open_count'][y] for y in years]
            assert all(a >= b for a, b in zip(counts, counts[1:])), name
            assert counts[0] <= cls['lakes']
        assert sum(v['lakes'] for k, v in surv.items() if k != 'all') == surv['all']['lakes']
        left = list(f['lake_capacity_left'].values())
        assert all(a >= b for a, b in zip(left, left[1:])), name
        cum = list(f['to_common_level_cumulative_gt'].values())
        assert all(a <= b for a, b in zip(cum, cum[1:])), name
        # Once everything upstream is full the common level takes the whole land's rock.
        assert math.isclose(f['to_common_level_gt_yr']['1e+09'], f['case']['rock_gt_yr'], rel_tol=1e-3), name
    med = {n: routing['fill'][n]['survival']['all']['median_life_yr'] for n in ('earliest', 'central', 'latest')}
    assert med['earliest'] < med['central'] < med['latest']
    rock_only = routing['fill']['central_rock_only']['survival']['all']['median_life_yr']
    assert med['central'] < rock_only


def test_channels_add_up(routing):
    ch = routing['channels']
    assert math.isclose(sum(b["length_km"] for b in ch["by_band"].values()), ch["total_km"], rel_tol=1e-5)
    for b in ch['by_band'].values():
        assert math.isclose(sum(v for *_, v in b['by_slope_and_log_qs']), b['length_km'], rel_tol=1e-4)
        assert math.isclose(sum(v for _, v in b['by_log_hardness']), b['length_km'], rel_tol=1e-4)


def test_deltas_grow_and_stay_within_their_waters(routing):
    d = routing['deltas']
    for name in rt.DELTA_SHARE:
        rows = [d[name]['by_year'][f'{y:g}'] for y in rt.DELTA_YEARS]
        for key in ('lake_delta_km2', 'common_level_delta_km2', 'mouths_over_1_km2'):
            assert all(a[key] <= b[key] for a, b in zip(rows, rows[1:])), (name, key)
        assert all(0.0 <= r['lake_delta_share_of_lake_area'] <= 1.0 for r in rows)
    ten = [d[n]['by_year']['10000']['lake_delta_km2'] for n in ('latest', 'central', 'earliest')]
    assert ten[0] < ten[1] < ten[2]
    big = [x['area_km2']['10000'] for x in d['largest_central']]
    assert big == sorted(big, reverse=True)
    for x in d['largest_central']:
        assert x['area_km2']['100'] <= x['area_km2']['1000'] <= x['area_km2']['10000']
        assert math.isclose(x['advance_km']['10000'], math.sqrt(2.0 * x['area_km2']['10000'] / math.pi), rel_tol=1e-4)
        assert x['area_km2']['10000'] <= x['receiving_water_km2'] * (1.0 + 1e-6)


@pytest.mark.skipif(not os.environ.get('TERLUNA_TERRAIN_REGENERATE')
                    or not all((rt.ROOT / p).exists() for p in rt.IGNORED),
                    reason='set TERLUNA_TERRAIN_REGENERATE=1 with the grid inputs restored (1-2 minutes, 5 GB)')
def test_routing_regenerates(routing):
    product, _ = rt.compute()
    fresh = json.loads(json.dumps(rt.rnd(product)))
    assert fresh == {k: v for k, v in routing.items() if k != 'producer'}


# ---- The second stage ------------------------------------------------------------------------------------------------

def test_product_binds_its_producer_and_inputs(product):
    assert product['schema'] == T.SCHEMA
    for path, h in {**product['producer']['files'], **product['producer']['inputs']}.items():
        assert digest(path) == h, path
    assert set(product['producer']['inputs']) == set(T.INPUTS)


def test_a_fresh_run_matches_the_product(product):
    assert json.loads(T.dumps(T.results())) == product


def test_crater_diffusion_is_scale_free_and_ordered(product):
    nd = product['small_craters']['dimensionless_kappa_t_over_d2']
    assert nd['0.5'] < nd['0.25'] < nd['0.1']
    assert 0.015 < nd['0.5'] < 0.04
    rows = product['small_craters']['craters']
    for a, b in zip(rows, rows[1:]):
        ratio = (b['diameter_m'] / a['diameter_m']) ** 2
        for k in a['half_depth_by_creep_yr']:
            assert math.isclose(b['half_depth_by_creep_yr'][k] / a['half_depth_by_creep_yr'][k], ratio, rel_tol=1e-6)
        assert all(b['half_depth_together_yr'][k] > a['half_depth_together_yr'][k] for k in a['half_depth_together_yr'])


def test_saturation_rises_with_load_and_the_rivers_sit_below_it(product):
    t = product['travertine']
    omegas = [v['omega'] for v in t['at_load'].values()]
    assert all(a < b for a, b in zip(omegas, omegas[1:]))
    need = t['load_needed']
    assert t['at_load']['1']['omega'] < 1.0 < need['designed'] < need['plitvice'] < need['tufa']
    for corner in t['corners'].values():
        shares = [r['rivers']['tufa'] for r in corner['by_year'].values()]
        assert all(a >= b - 1e-12 for a, b in zip(shares, shares[1:]))
    assert t['corners']['high']['by_year']['1000']['rivers']['tufa'] > 0.9
    assert t['corners']['low']['by_year']['0']['rivers']['tufa'] < 0.01


def test_upkeep_scales_with_the_rock(product):
    rates = product['upkeep']['rates']
    base = rates['10 m/Myr']
    for row in rates.values():
        k = row['rock_gt_yr'] / base['rock_gt_yr']
        for key in ('lift_gw_in_place', 'lift_gw_to_heads', 'haulage_tonne_km_yr', 'spread_in_place_mm_per_century'):
            assert math.isclose(row[key], k * base[key], rel_tol=1e-4), key
        assert row['lift_gw_to_heads'] > row['lift_gw_in_place']
    assert math.isclose(product['upkeep']['dredging_fuel_mj_per_kg'], 43.0, rel_tol=0.02)


def test_dams_hold_on_less_under_lunar_floods(product):
    cases = product['dam_builders']['cases']
    for limit in ('dams_persist', 'with_blowouts'):
        for regime in ('steady', 'flashy'):
            lunar = cases[f'{limit}/{regime}/equal_shields_number']['dams_hold_km']
            earth = cases[f'{limit}/{regime}/earth_discharge_times_slope']['dams_hold_km']
            assert 0 < lunar < earth < product['dam_builders']['channel_km']
        assert cases[f'{limit}/flashy/equal_shields_number']['dams_hold_km'] <= \
            cases[f'{limit}/steady/equal_shields_number']['dams_hold_km']
    assert cases['dams_persist/steady/equal_shields_number']['dams_hold_km'] <= \
        cases['with_blowouts/steady/equal_shields_number']['dams_hold_km']


def test_the_seas_fill_within_the_sea_chemistry_s_range(product):
    s = product['seas']
    lo, hi = s['sea_chemistry_box_seas_full_yr']
    assert s['earliest']['seas_full_yr'] < s['central']['seas_full_yr'] < s['latest']['seas_full_yr']
    assert 0.5 * lo < s['earliest']['seas_full_yr'] and s['latest']['seas_full_yr'] < 2.0 * hi


def test_cases_match_the_sea_chemistry_and_conservation():
    from research.studies.sea_chemistry import run as S
    dens = sorted({c['dry_density_kg_m3'] for c in rt.CASES.values()})
    assert dens[0] == S.DRY_DENSITY[0] and dens[-1] == S.DRY_DENSITY[1]
    gates = json.loads((rt.ROOT / rt.GATES).read_text())
    assert sorted(float(k.split()[0]) for k in gates['sediment_return']) == \
        sorted({c['erosion_m_myr'] for c in rt.CASES.values()})
