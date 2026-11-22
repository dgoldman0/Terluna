"""Checks of the canopy colony study: the column's physics by hand, the parts' sizing rules, and the products against
themselves and their inputs."""
import json
import os
from pathlib import Path

import numpy as np
import pytest

from biosphere.ecology import lake as L
from research.studies.canopy_colony import column as CL
from research.studies.canopy_colony import currents as CU
from research.studies.canopy_colony import run as R
from research.studies.canopy_colony import seas as SE
from shared.constants import MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


def cylinder(n=20, thick=50.0, k=1e-5, saturation=6.0, taper=0.0):
    mid = (np.arange(n) + 0.5) * thick
    volume = thick * (1.0 - taper * np.arange(n) / n)
    return CL.from_level_volumes(mid, volume, np.full(n, k), saturation)


# --- the column, by hand

def test_the_moons_sinking_is_gravity_over_viscosity():
    assert R.sinking_stretch(R.EARTH_FLUX_WATER_C, gravity=STANDARD_GRAVITY) == pytest.approx(1.0)
    moon = R.sinking_stretch()
    assert moon == pytest.approx(STANDARD_GRAVITY / MOON_SURFACE_GRAVITY * L.viscosity(28.0) / L.viscosity(7.5))
    assert moon == pytest.approx(3.60, abs=0.02)
    assert R.sinking_stretch(decay_q10=2.0) == pytest.approx(moon * 2.0 ** 2.05)
    assert CL.flux_share(400.0) == pytest.approx(4.0 ** -0.86)                     # Earth: 30% of the export passes 400 m
    assert CL.flux_share(400.0, moon) == pytest.approx((1.0 + 3.0 * moon) ** -0.86)
    assert CL.flux_share(50.0, moon) == 1.0                                        # nothing is lost above the reference


def test_a_column_follows_its_seas_floor():
    col = cylinder(taper=0.5)
    assert col.area[0] == 1.0 and col.area[-1] == 0.0 and (np.diff(col.area) <= 0).all()
    assert col.volume.sum() == pytest.approx(50.0 * (1.0 - 0.5 * np.arange(20) / 20).sum())
    by_area = CL.from_areas((np.arange(8) + 0.5) * 100.0, [0.0, 400.0, 600.0], [10.0, 10.0, 0.0], np.ones(8), 7.0)
    assert by_area.n == 6 and by_area.edge[-1] == 600.0                            # levels under the deepest water go
    assert by_area.volume.sum() == pytest.approx(400.0 + 0.5 * 200.0)              # a cylinder to 400 m, a cone below


def test_the_demand_is_what_the_flux_loses_and_the_feeders_take_their_share():
    col = cylinder(taper=0.5)
    at = lambda z: float(np.interp(z, col.edge, col.area))                         # noqa: E731
    for stretch in (1.0, 3.6):
        d, caught = CL.demand(col, 100.0, stretch)
        assert caught == 0.0 and (d[col.edge[1:] <= 100.0] == 0.0).all()
        assert (d * col.volume).sum() == pytest.approx(100.0 * at(100.0))          # all that passes 100 m is spent below
        half, caught = CL.demand(col, 100.0, stretch, 0.5, 250.0)
        assert caught == pytest.approx(0.5 * 100.0 * CL.flux_share(250.0, stretch) * at(250.0))
        assert (half * col.volume).sum() + caught == pytest.approx(100.0 * at(100.0))
        deep = col.edge[:-1] >= 250.0
        assert half[deep] == pytest.approx(0.5 * d[deep]) and half[col.edge[1:] <= 250.0] == pytest.approx(d[col.edge[1:] <= 250.0])
    slow, _ = CL.demand(col, 100.0, 3.6)
    fast, _ = CL.demand(col, 100.0, 1.0)
    assert (slow * col.volume)[col.mid > 400.0].sum() < 0.5 * (fast * col.volume)[col.mid > 400.0].sum()


def test_the_steady_column_carries_each_levels_demand_up_through_the_water_above():
    col = cylinder(n=6, thick=100.0, k=1e-5)
    d = np.zeros(6)
    d[4] = 0.02                                                                    # g/m3 a year in one level
    spent = CL.steady(col, d, CL.operator(col))
    step = 0.02 * 100.0 * 100.0 / (1e-5 * CL.YEAR_S)                               # flux times spacing over diffusivity
    assert spent[0] == 0.0 and np.diff(spent)[:4] == pytest.approx(step) and spent[5] == pytest.approx(spent[4])
    end, at, means = CL.integrate(col, d, 4000.0, marks=(10.0,), dt=0.5, yearly=2)
    assert end == pytest.approx(spent, rel=1e-3) and (at[10.0] < spent + 1e-12).all() and len(means) == 2
    held, _, _ = CL.integrate(col, 100.0 * d, 4000.0, dt=0.5)
    assert held.max() == col.saturation                                            # water cannot spend more than it held


def test_the_trunk_replaces_its_foots_water_and_lifts_the_water_above():
    col = cylinder(n=6, thick=100.0, k=1e-10)                                      # next to no mixing: the trunk alone
    d = np.zeros(6)
    d[3] = 0.05
    op = CL.operator(col, 2.0, (300.0, 400.0))                                     # 2 m a year let out in level 3
    spent = CL.steady(col, d, op)
    assert spent[3] == pytest.approx(0.05 * 100.0 / 2.0, rel=1e-3)                 # demand over the flushing rate
    assert spent[1:3] == pytest.approx(spent[3], rel=1e-3) and spent[0] == 0.0     # the rising water keeps what it spent
    assert spent[4:] == pytest.approx(spent[3], rel=1e-3)                          # still water below comes level with it
    spread = CL.operator(col, 2.0, (300.0, 500.0))                                 # let out over two levels
    assert -spread[4, 4] * col.volume[4] == pytest.approx(1.0, rel=1e-3) and spread[3, 4] * col.volume[3] == pytest.approx(1.0, rel=1e-3)
    assert -spread[3, 3] * col.volume[3] == pytest.approx(2.0, rel=1e-3)           # its own metre and the one from below


def test_the_fit_finds_the_share_a_run_was_made_with():
    col = cylinder(n=12, thick=40.0, k=3e-4)
    d, _ = CL.demand(col, 90.0, 1.0)
    k = col.k.copy()
    k[col.mid < 300.0] *= 0.2
    _, _, used = CL.integrate(col.with_k(k), d, 3, yearly=3)
    share, err = CL.calibrated_share(col, d, [list(u) for u in used])
    assert share == pytest.approx(0.2, rel=0.06) and err < 0.02


# --- the parts

def test_the_trunks_channels_carry_what_their_walls_can_cool():
    col = cylinder(n=20, thick=50.0)
    warmth = np.linspace(30.0, 27.0, 20)
    water = {'250-500 m': dict(share_of_area_denser=0.2, excess_density_kg_m3=dict(median=-0.2, p5=-1.0))}
    out = R.trunk_engineering(col, 2.0, (250.0, 500.0), warmth, dict(surface_water_cooled_to_the_foot=water))
    gradient = (out['surface_c'] - out['foot_c']) / out['mean_length_m']
    lag = out['lags']['0.3 K']
    assert lag['channel_m3_a_day'] == pytest.approx(R.NUSSELT * np.pi * R.HEAT_DIFFUSIVITY_M2_S * 0.3 / gradient * 86400.0)
    assert lag['channels_per_km2'] * lag['channel_m3_a_day'] * 365.25 == pytest.approx(2.0e6)
    assert out['lags']['0.03 K']['channels_per_km2'] == pytest.approx(10.0 * lag['channels_per_km2'])
    assert out['pumping']['head_pa']['median'] == pytest.approx(0.2 * MOON_SURFACE_GRAVITY * out['mean_length_m'])
    assert out['uncooled_pipe_warms_the_foot_k_a_century'] == pytest.approx(
        lag['foot_warms_k_a_century_if_nothing_carries_it_off'] * (out['surface_c'] - out['foot_c']) / 0.3)


def test_the_feeders_surface_takes_the_share_asked():
    currents = dict(central_run='r', runs=dict(r=dict(water_past_a_drifting_canopy_m_s=dict(feeders=dict(median=0.04, p10=0.01)))))
    net = R.feeders_net(0.5, 3.6, currents)
    for label, need in net['surface_times_efficiency_m2_per_m2'].items():
        u = 0.04 if label.startswith('median') else 0.01
        w = float(label.split(', ')[1].split(' ')[0]) / 86400.0
        assert 1.0 - np.exp(-need * u / w) == pytest.approx(0.5)
    assert net['sinking_m_day'] == pytest.approx([20.0 / 3.6, 100.0 / 3.6])
    assert R.feeders_net(0.75, 3.6, currents)['facing_surface_share_of_the_sea'][0] == pytest.approx(
        2.0 * net['facing_surface_share_of_the_sea'][0])


def test_the_capsules_calcium_is_counted_against_the_rivers_and_the_stock():
    sea = dict(rivers_m_yr=0.5, mean_depth_m=1000.0)
    out = R.calcium(0.07, sea, all_rivers_m3s=250000.0, stock_mol_m3=10.0)
    assert out['capsules_made_mol_o2_m2_yr'] == pytest.approx(0.1) and out['calcium_mol_m2_yr'] == pytest.approx(0.2)
    per_m3 = 150e12 / R.CALCITE_G_PER_MOL / (250000.0 * CL.YEAR_S)
    assert out['rivers_calcium_mol_m2_yr'] == pytest.approx(0.5 * per_m3) and 0.15 < per_m3 < 0.25
    assert out['stock_drawn_percent_a_millennium'] == pytest.approx(100.0 * (0.2 - 0.5 * per_m3) * 1000.0 / 1e4)
    assert out['light_carbon_forgone_g_m2_yr'] == pytest.approx(0.1 * CL.G_C_PER_MOL)


def test_a_sea_is_measured_by_depth_and_latitude():
    out = SE.sea_table(np.array([5.0, 150.0, 900.0]), np.array([1e6, 1e6, 2e6]), np.array([10.0, 50.0, 60.0]))
    assert out['area_km2'] == 4.0 and out['mean_depth_m'] == pytest.approx((5.0 + 150.0 + 1800.0) / 4.0)
    below = {b['depth_m']: b['area_km2'] for b in out['below']}
    assert below[0] == 4.0 and below[100] == 3.0 and below[850] == 2.0 and 1000 not in below
    assert out['share_poleward_of_45'] == 0.75 and out['latitudes'] == [10.0, 60.0]


def test_the_water_past_a_canopy_is_the_difference_of_two_layers():
    mid = np.array([30.0, 200.0, 400.0])
    volume = np.ones((3, 2, 2))
    volume[2, 1, 1] = 0.0                                                          # one column ends above the foot
    u = np.stack([np.full((2, 2), 0.05), np.full((2, 2), 0.01), np.zeros((2, 2))])
    v = np.zeros_like(u)
    assert CU.layer_mean(u, volume, mid, (150.0, 250.0)) == pytest.approx(np.full((2, 2), 0.01))
    assert np.isnan(CU.layer_mean(u, volume, mid, (300.0, 500.0))[1, 1])
    run = dict(u=u, v=v, volume=volume, mid=mid, area=np.ones((2, 2)), year=1.0, file='x',
               salt=np.stack([np.full((2, 2), 35.5), np.full((2, 2), 35.0), np.full((2, 2), 35.0)]),
               theta=np.stack([np.full((2, 2), 30.0), np.full((2, 2), 28.0), np.full((2, 2), 27.0)]))
    assert CU.weighted_percentiles(np.array([1.0, 2.0, 3.0, np.nan]), np.ones(4), (50,))[0] == pytest.approx(1.5, abs=0.51)
    water = CU.trunk_water(run, lambda s, t: 1000.0 + 0.75 * s - 0.3 * t, (300.0, 500.0))
    assert water['share_of_area_denser'] == 1.0 and water['excess_density_kg_m3']['median'] == pytest.approx(0.375)
    assert water['temperature_step_k'] == pytest.approx(3.0)


# --- the products

def test_the_product_binds_its_inputs(product):
    for rel, digest in product['inputs'].items():
        assert R.sha256(ROOT / rel) == digest, rel
    source = ROOT / R.MIXING_SOURCE                    # present once domain/hydrosphere's reading has joined this branch
    if source.exists():
        assert source.read_bytes() == R.MIXING.read_bytes()


def test_the_column_follows_the_paused_run(product):
    col = product['south_sea']['column']
    fit = col['fit']
    assert 0.05 < fit['share_of_the_mean_diffusivity_above_300_m'] < 0.3 and fit['log_error'] < 0.2
    mid = np.array(col['level_mid_m'])
    run, mine = np.array(fit['run_used_mg_l'], float), np.array(fit['column_used_mg_l'])
    deep = (mid > 320.0) & (mid < 900.0)                                          # the still layer, outside the fit
    assert mine[:, deep] == pytest.approx(run[:, deep], rel=0.25)
    k = np.array(col['diffusivity_m2_s'])
    assert (k[(mid > 390.0) & (mid < 900.0)] < 1e-6).all() and k[mid > 1400.0].min() > 1e-3
    assert col['open_area_at_each_edge'][0] == 1.0 and 1900.0 < col['mean_depth_m'] < 2200.0


def test_slower_sinking_keeps_more_oxygen_and_more_export_less(product):
    unaided = product['south_sea']['unaided']
    for name in product['setting']['exports_g_o2_m2_yr']:
        for mark in ('100 years', '300 years', '1000 years'):
            eels = [unaided[form][name][mark]['eels_layer_mg_l'] for form in ('earth', 'moon', 'moon_warm')]
            assert eels == sorted(eels)
    for form in unaided:
        by_export = [unaided[form][name]['300 years']['eels_layer_mg_l'] for name in ('poorest', 'nutrient-poor', 'middling')]
        assert by_export == sorted(by_export, reverse=True)
    placement = product['south_sea']['placement']
    for form, p in placement.items():
        bands = {k: v for k, v in p['share_of_demand'].items() if k != 'of it on the floor'}
        assert sum(bands.values()) == pytest.approx(1.0)
    assert placement['moon']['share_of_demand']['100-250 m'] > 0.75 > placement['earth']['share_of_demand']['100-250 m']
    assert placement['earth']['still_sinking_at_m']['400'] == pytest.approx(4.0 ** -0.86)


def test_each_sized_system_holds_its_sea(product):
    setting, systems = product['setting'], product['south_sea']['systems']
    for form, by_export in systems.items():
        for name, cases in by_export.items():
            export = setting['exports_g_o2_m2_yr'][name]
            for label, s in cases.items():
                w2, w4 = s['trunk_m_yr_for_2_mg_l_everywhere'], s['trunk_m_yr_with_the_eels_layer_at_4_mg_l']
                assert w2 is not None and w4 is not None and 0.0 <= w2 <= w4 + 1e-9
                steady = s['steady_with_that_trunk']
                assert steady['lowest_mg_l'] >= setting['hypoxic_mg_l'] - 1e-6
                assert steady['eels_layer_mg_l'] >= setting['eels_mg_l'] - 1e-6
                assert s['capsules_needed_mol_o2_m2_yr'] < export / CL.G_O2_PER_MOL * 0.3
            for trunk in ('short', 'long'):
                need = [cases[f'catch {c:g}, {trunk} trunk']['capsules_needed_mol_o2_m2_yr'] for c in (0.25, 0.5, 0.75)]
                flow = [cases[f'catch {c:g}, {trunk} trunk']['trunk_m_yr_with_the_eels_layer_at_4_mg_l'] for c in (0.25, 0.5, 0.75)]
                assert need == sorted(need, reverse=True) and flow == sorted(flow, reverse=True)
                assert need[0] == pytest.approx(3.0 * need[2])                   # what is left is what was not caught
            assert cases['catch 0.5, long trunk']['capsules_needed_mol_o2_m2_yr'] < cases['catch 0.5, short trunk']['capsules_needed_mol_o2_m2_yr']
    poor, rich = systems['moon']['poorest']['catch 0.5, long trunk'], systems['moon']['middling']['catch 0.5, long trunk']
    ratio = setting['exports_g_o2_m2_yr']['middling'] / setting['exports_g_o2_m2_yr']['poorest']
    assert rich['capsules_needed_mol_o2_m2_yr'] == pytest.approx(ratio * poor['capsules_needed_mol_o2_m2_yr'])
    assert rich['caught_mol_o2_m2_yr'] == pytest.approx(ratio * poor['caught_mol_o2_m2_yr'])


def test_the_central_system_needs_every_part(product):
    central = product['south_sea']['central']
    years = central['thousand_years']
    whole = years['whole']['1000 years']
    assert whole['under_2_mg_l_between_m'] is None and whole['eels_layer_mg_l'] > 4.0
    for missing in ('without the feeders', 'without the trunk', 'feeders alone'):
        assert years[missing]['1000 years']['under_2_mg_l_between_m'] is not None
        assert years[missing]['1000 years']['eels_layer_mg_l'] < whole['eels_layer_mg_l']
    late = years['without the capsules']['first_millennium_with_under_2_mg_l_below_the_foot']
    assert late is not None and late > 1000.0                                      # the stock below lasts, then fails
    loop = central['loop']
    shares = loop['nutrients']['regenerated_share_of_the_export']
    assert sum(shares['without_the_colony'].values()) == pytest.approx(1.0, abs=1e-6)
    assert sum(shares['with_the_feeders'].values()) + shares['taken_by_the_feeders'] == pytest.approx(1.0, abs=1e-6)
    assert loop['catch']['g_c_m2_yr'] == pytest.approx(loop['export']['g_c_m2_yr'] * loop['catch']['share_of_the_export'], rel=1e-3)
    assert loop['trunk']['oxygen_given_up_at_the_foot_mol_m2_yr'] < loop['trunk']['oxygen_carried_mol_m2_yr']
    assert max(loop['trunk']['pumping_share_of_the_catchs_energy'].values()) < 0.05
    p = loop['nutrients']['phosphorus']
    assert p['lost_below_the_trunks_reach_mg_m2_yr'] < 0.2 * p['without_the_colony_lost_below_250_m_mg_m2_yr']
    assert p['river_phosphorus_that_replaces_the_loss_mg_m3'] * p['rivers_m_yr'] == pytest.approx(p['lost_below_the_trunks_reach_mg_m2_yr'])


def test_the_fresh_seas_stand_on_their_own_floors(product):
    fresh = product['fresh_seas']
    assert fresh['saturation_mg_l'] == pytest.approx(7.66, abs=0.01)               # the sea oxygen study's fresh water
    seas = json.loads(R.SEAS.read_text())['seas']
    assert set(fresh['seas']) == {k for k, v in seas.items() if not v['salt']} and len(fresh['seas']) == 6
    for name, s in fresh['seas'].items():
        assert s['mean_depth_m'] == pytest.approx(seas[name]['mean_depth_m'])
        for form, case in s['cases'].items():
            sized = case['system']
            assert sized['trunk_m_yr_for_2_mg_l_everywhere'] is not None
            assert case['calcium']['rivers_calcium_mol_m2_yr'] == pytest.approx(
                case['calcium']['rivers_calcium_mol_m3'] * seas[name]['rivers_m_yr'])
        assert s['cases']['moon']['system']['capsules_needed_mol_o2_m2_yr'] > s['cases']['moon_warm']['system']['capsules_needed_mol_o2_m2_yr']
    shallow, deep = fresh['seas']['near_side_south'], fresh['seas']['moscoviense']
    assert shallow['share_of_floor_above_250_m'] > 3.0 * deep['share_of_floor_above_250_m']
    assert shallow['cases']['moon']['system']['capsules_needed_mol_o2_m2_yr'] < deep['cases']['moon']['system']['capsules_needed_mol_o2_m2_yr']


def test_the_product_regenerates(product):
    again = json.loads(json.dumps(R.build(log=lambda *a: None)))
    for key in ('setting', 'south_sea', 'fresh_seas'):
        assert again[key] == product[key], key


@pytest.mark.skipif(not os.environ.get('TERLUNA_CANOPY_COLONY_REGENERATE'),
                    reason='reads the design routing and the south sea runs on the research drive; set '
                           'TERLUNA_CANOPY_COLONY_REGENERATE=1 (the currents in the ocean environment)')
def test_the_seas_and_the_currents_regenerate():
    seas = json.loads(R.SEAS.read_text())
    assert json.loads(json.dumps(SE.build(log=lambda *a: None)['seas'])) == seas['seas']
    currents = json.loads(R.CURRENTS.read_text())
    again = json.loads(json.dumps(CU.build()))
    assert again['runs'] == currents['runs'] and again['surface_water_cooled_to_the_foot'] == currents['surface_water_cooled_to_the_foot']
