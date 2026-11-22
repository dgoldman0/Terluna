"""Checks of the seas' phosphorus budget: the tracer on a made-up drainage, the budget's rule by hand, and the
products against themselves and their inputs."""
import json
import os
from pathlib import Path

import numpy as np
import pytest

from research.studies.canopy_colony import column as CL
from research.studies.canopy_colony import run as CC
from research.studies.sea_nutrients import dust as DU
from research.studies.sea_nutrients import loads as LO
from research.studies.sea_nutrients import run as R

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


@pytest.fixture(scope='module')
def loads():
    return json.loads(R.LOADS.read_text())


def cylinder(n=24, thick=50.0, k=1e-5, saturation=6.0):
    mid = (np.arange(n) + 0.5) * thick
    return CL.from_level_volumes(mid, np.full(n, thick), np.full(n, k), saturation)


# --- the tracer

def test_water_is_carried_down_in_order_and_a_lake_keeps_its_share():
    # cells 0 and 1 drain to 2, 2 to 3, 3 to the outlet 4; cell 5 drains to 3 as well
    receiver = np.array([2, 2, 3, 4, 4, 3])
    donor = np.array([True, True, True, True, False, True])
    fronts = LO.kahn_fronts(receiver, donor)
    order = np.concatenate(fronts)
    assert sorted(order.tolist()) == [0, 1, 2, 3, 5]
    assert order.tolist().index(2) > max(order.tolist().index(0), order.tolist().index(1))
    assert order.tolist().index(3) > max(order.tolist().index(2), order.tolist().index(5))
    source = np.array([1.0, 2.0, 0.0, 0.0, 0.0, 4.0])
    passing = LO.carry(receiver, fronts, source, np.ones(6))
    assert passing.tolist() == [1.0, 2.0, 3.0, 7.0, 7.0, 4.0]                      # everything reaches the outlet
    factor = np.array([1.0, 1.0, 0.25, 1.0, 1.0, 0.0])                             # a lake at 2 keeps three quarters; 5 keeps all
    kept = LO.carry(receiver, fronts, source, factor)
    assert kept[3] == pytest.approx(0.75) and kept[4] == pytest.approx(0.75)


def test_a_lake_is_the_cells_at_one_level():
    cols = 4
    h = np.array([[10.0, 2.0, 3.0, 10.0], [10.0, 10.0, 10.0, 10.0]], np.float32)
    rt = dict(lake_cell=np.array([1, 2], np.int32), lake_level_m=np.array([5.0, 5.0], np.float32),
              receiver=np.array([1, 2, 3, 3, 4, 5, 6, 7], np.int32), discharge_m3s=np.array([1.0, 2.0, 3.0, 3.0, 0, 0, 0, 0]))
    body, volume, area, outflow, leaving = LO.lake_bodies(rt, h, np.array([100.0, 100.0]), cols)
    assert body.tolist() == [-1, 0, 0, -1, -1, -1, -1, -1]
    assert volume.tolist() == [100.0 * 3.0 + 100.0 * 2.0] and area.tolist() == [200.0]
    assert leaving.tolist() == [2] and outflow.tolist() == [3.0]                  # the water leaves by cell 2 alone


# --- the budget's rule, by hand

def test_the_supply_sets_the_export_and_the_oxygen_spent_below():
    assert R.P_MG_PER_G_C == pytest.approx(1e3 * 30.974 / (106.0 * 12.011))        # 24.3 mg of phosphorus a gram of carbon
    assert R.O2_MMOL_PER_MG_P == pytest.approx(106.0 / 30.974)                    # 3.42 mmol of oxygen a mg
    shares = dict(below_250_m=0.2, below_the_trunks_reach=0.02)
    safe = {'300 years': dict(export_g_c_m2_yr=20.0), '1000 years': dict(export_g_c_m2_yr=9.0), '3000 years': dict(export_g_c_m2_yr=5.0)}
    out = R.what_a_supply_does(10.0, shares, safe)
    assert out['export_g_c_m2_yr'] * 0.2 * R.P_MG_PER_G_C == pytest.approx(10.0)   # what passes 250 m is what came in
    assert out['export_with_the_colony_g_c_m2_yr'] == pytest.approx(10.0 * out['export_g_c_m2_yr'])
    assert out['oxygen_spent_below_250_m_mol_m2_yr'] == pytest.approx(10.0 * 106.0 / 30.974 / 1e3)
    assert out['unaided_sea_holds_through_years'] == 3000.0                        # 2.05 g C is under every limit
    assert R.what_a_supply_does(60.0, shares, safe)['unaided_sea_holds_through_years'] == 300.0
    assert R.what_a_supply_does(200.0, shares, safe)['unaided_sea_holds_through_years'] == 0.0


def test_the_colony_shares_what_passes_250_m():
    col = cylinder()
    stretch = CC.sinking_forms()['moon']
    s = R.passing_shares(col, stretch)
    assert s['below_250_m'] == pytest.approx(float(CL.flux_share(250.0, stretch)))    # a cylinder loses none to slopes
    assert s['taken_by_the_feeders'] + s['returned_by_the_trunk'] + s['below_the_trunks_reach'] == pytest.approx(s['below_250_m'])
    assert s['below_the_trunks_reach'] == pytest.approx(0.5 * float(CL.flux_share(900.0, stretch)))


def test_the_safe_export_is_the_edge_of_holding():
    col = cylinder(k=2e-6)
    stretch = CC.sinking_forms()['moon']
    edge = {years: R.safe_export(col, stretch, years) for years in (300.0, 1000.0)}
    assert edge[300.0] > edge[1000.0] > 0.0                                         # a longer hold takes a poorer sea
    assert R.lowest_oxygen(col, 0.98 * edge[1000.0], stretch, 1000.0) >= CC.HYPOXIC_MG_L
    assert R.lowest_oxygen(col, 1.05 * edge[1000.0], stretch, 1000.0) < CC.HYPOXIC_MG_L


def test_salt_water_frees_part_of_the_particles_and_fresh_water_none():
    load = dict(area_m2=1e12, share_of_the_lands_dissolved_loss=0.01, dissolved_share_from_other_seas=0.001,
                share_of_the_lands_particles=dict(low=0.002, central=0.003, high=0.004))
    land = dict(dissolved_kg_yr=[1e8, 2e8], particulate_kg_yr=[1e9, 4e9])
    fresh = R.supply_cases(load, land, [0.001, 0.02], salt=False)
    salt = R.supply_cases(load, land, [0.001, 0.02], salt=True)
    assert fresh['rivers_dissolved'] == pytest.approx([1.1, 2.2]) and fresh['rivers_particles'] == pytest.approx([2.0, 16.0])
    assert fresh['freed_from_particles'] == [0.0, 0.0] and fresh['rivers_alone'] == pytest.approx([1.1, 2.2])
    assert salt['freed_from_particles'] == pytest.approx([0.2, 4.8]) and salt['rivers_alone'] == pytest.approx([1.3, 7.0])
    assert fresh['rivers_and_dust'] == pytest.approx([2.1, 22.2])


def test_dust_is_a_share_of_the_sand_carried_and_falls_over_the_whole_moon():
    assert DU.dust_raised(36.0, 1.65e-6) == pytest.approx(1.65e-6 * 36.0 * 1600.0)        # kg per m2 of moving ground a year
    assert DU.phosphorus_fall(1.0e9, 5.0e-4, 4.0e13) == pytest.approx(1e6 * 1.0e9 * 5.0e-4 / 4.0e13)   # mg per m2 a year
    sand = dict(winds=dict(over_sand=dict(central=dict(by_group={
        'hyper-arid': dict(drift_potential_m2_yr=dict(mean=40.0)), 'arid': dict(drift_potential_m2_yr=dict(mean=30.0)),
        'semi-arid': dict(drift_potential_m2_yr=dict(mean=20.0))}))),
        land=dict(by_class_km2={'hyper-arid': 1.0, 'arid': 1.0, 'semi-arid': 2.0}))
    assert DU.drift_potential(sand, 'central') == pytest.approx((40.0 + 30.0 + 2 * 20.0) / 4.0)
    out = json.loads(DU.OUT.read_text())
    for rel, digest in out['inputs'].items():
        assert DU.sha256(ROOT / rel) == digest, rel
    central, outer = out['natural']['phosphorus_mg_m2_yr']['central'], out['natural']['phosphorus_mg_m2_yr']['outer']
    assert outer[0] < central[0] < central[1] < outer[1] < 5.0
    assert central[1] == pytest.approx(10.0 * central[0] * out['setting']['soil_phosphorus_grade'][1] / out['setting']['soil_phosphorus_grade'][0])
    assert out['inherited_mg_m2_yr'] == pytest.approx([0.7, 28.0]) and central[1] < out['inherited_mg_m2_yr'][0]
    assert out['with_earths_alpha_factor'] == pytest.approx(1.0 / 0.165)
    raised = [c['dust_kg_per_m2_of_moving_ground_yr'] for c in out['natural']['central']]
    assert 0.05 < raised[0] < raised[1] < 2.0                                      # the rate of Earth's active deserts
    assert json.loads(json.dumps(DU.build())) == out


def test_nitrogen_is_taken_up_with_the_phosphorus():
    assert R.N_MG_PER_MG_P == pytest.approx(16.0 * 14.007 / 30.974)
    load = dict(share_of_the_lands_runoff=0.2, of_it_crossed_no_lake=0.05)
    sea = dict(area_km2=1.0e6)
    out = R.nitrogen(load, sea, dict(rivers_and_dust=[1.0, 4.0]), land_km2=2.0e7, lightning=[0.1, 5.0])
    assert out['taken_up_with_the_phosphorus_mg_m2_yr'] == pytest.approx([R.N_MG_PER_MG_P, 4.0 * R.N_MG_PER_MG_P])
    assert out['rivers_mg_m2_yr'] == pytest.approx([1e3 * 0.1 * 2.0e13 * 0.01 / 1e12, 1e3 * 0.5 * 2.0e13 * 0.2 / 1e12])
    assert out['arrives_mg_m2_yr'] == pytest.approx([20.1, 2005.0])
    assert out['fixing_needed_in_the_sea_mg_m2_yr'] == pytest.approx([0.0, 4.0 * R.N_MG_PER_MG_P - 20.1])


# --- the products

def test_the_lakes_keep_most_of_the_lands_phosphorus(loads):
    seas = loads['seas']
    dissolved = sum(s['share_of_the_lands_dissolved_loss'] for s in seas.values())
    assert 0.03 < dissolved < 0.10                                                 # the phosphorus study found 6.5% on the atlas
    for case in ('low', 'central', 'high'):
        assert 0.02 < sum(s['share_of_the_lands_particles'][case] for s in seas.values()) < 0.06
    for s in seas.values():
        p = s['share_of_the_lands_particles']
        assert p['low'] <= p['central'] <= p['high'] and s['share_of_the_lands_dissolved_loss'] <= s['share_of_the_lands_runoff']
    runoff = sum(s['share_of_the_lands_runoff'] for s in seas.values())
    assert runoff + loads['census']['runoff_ending_in_closed_lakes'] == pytest.approx(1.0, abs=2e-3)
    record = json.loads((ROOT / LO.ROUTING_RECORD).read_text())['summary']['weirs']
    assert loads['transfers']['moscoviense']['outflow_m3s'] == pytest.approx(record['moscoviense_surplus_m3s'], abs=2.0)
    assert loads['water_entering_from_upper_seas_m3s'] > sum(t['outflow_m3s'] for k, t in loads['transfers'].items() if k != 'near_side_main')
    for t in loads['transfers'].values():
        assert t['passes_on'] == pytest.approx(1.0 / (1.0 + np.sqrt(t['residence_yr']))) and t['passes_on'] < 0.05


def test_the_product_binds_its_inputs(product):
    for rel, digest in product['inputs'].items():
        assert R.sha256(ROOT / rel) == digest, rel


def test_the_rivers_bring_far_less_than_the_nutrient_poor_state_takes(product, loads):
    setting = product['setting']
    for name, sea in product['seas'].items():
        supply, safe = sea['supply_mg_m2_yr'], sea['supply_the_unaided_sea_can_take']
        needed = setting['supply_that_would_sustain_the_reference_export_mg_m2_yr'][name]
        assert supply['rivers_alone'][1] < 0.05 * needed and supply['rivers_and_dust'][1] < 0.3 * needed
        assert supply['rivers_alone'][1] < safe['3000 years']['supply_mg_m2_yr']   # on its rivers a sea holds for millennia
        limits = [safe[h]['supply_mg_m2_yr'] for h in ('300 years', '1000 years', '3000 years')]
        assert limits == sorted(limits, reverse=True)
        for case in sea['what_the_supply_does'].values():
            for end in case.values():
                share = sea['shares_of_the_exports_phosphorus']['below_250_m']
                assert end['export_g_c_m2_yr'] * share * R.P_MG_PER_G_C == pytest.approx(end['supply_mg_m2_yr'])
                assert end['oxygen_spent_below_250_m_mol_m2_yr'] == pytest.approx(end['supply_mg_m2_yr'] * R.O2_MMOL_PER_MG_P / 1e3)
        assert sea['what_the_supply_does']['rivers_alone']['high']['unaided_sea_holds_through_years'] == 3000.0
        assert sea['what_the_supply_does']['rivers_and_dust']['high']['unaided_sea_holds_through_years'] == 3000.0
        assert supply['dust'] == setting['dust_mg_m2_yr']['used'] and supply['dust_inherited'] == pytest.approx([0.7, 28.0])
        n = sea['nitrogen']
        assert n['arrives_mg_m2_yr'][0] <= n['arrives_mg_m2_yr'][1] and 0.0 <= n['fixing_needed_in_the_sea_mg_m2_yr'][1] < 25.0
        assert n['taken_up_with_the_phosphorus_mg_m2_yr'][1] == pytest.approx(supply['rivers_and_dust'][1] * R.N_MG_PER_MG_P)
        d = sea['delivered_water']
        assert d['70 mg/m3']['oxygen_it_takes_below_mol_m2'] == pytest.approx(70.0 * d['1 mg/m3']['oxygen_it_takes_below_mol_m2'])
        assert 200.0 < sea['water_above_250_m_per_m2'] <= 250.0
    total = product['all_seas']['kt_yr']['rivers_dissolved']
    share = product['all_seas']['share_of_the_lands_loss_reaching_the_seas']['dissolved']
    assert total[0] >= share * setting['land_loss_kg_yr']['dissolved'][0] / 1e6 * 0.999   # with what the seas pass on


def test_the_product_regenerates(product):
    again = json.loads(json.dumps(R.build(log=lambda *a: None)))
    for key in ('setting', 'all_seas', 'seas'):
        assert again[key] == product[key], key


@pytest.mark.skipif(not os.environ.get('TERLUNA_SEA_NUTRIENTS_REGENERATE'),
                    reason='reads the design routing on the research drive; set TERLUNA_SEA_NUTRIENTS_REGENERATE=1')
def test_the_loads_regenerate(loads):
    again = json.loads(json.dumps(LO.build(log=lambda *a: None)))
    for key in ('census', 'transfers', 'seas'):
        assert again[key] == loads[key], key
