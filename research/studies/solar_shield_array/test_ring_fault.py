import json
import math

import numpy as np
import pytest

from shared import constants as K
from shared.provenance import constants_changed
from . import ring_fault as rf


def product():
    return json.loads(rf.OUT.read_text())


needs_product = pytest.mark.skipif(not rf.OUT.exists(), reason='ring_fault.json not written')


@needs_product
def test_product_is_current():
    p = product()
    assert p['schema'] == rf.SCHEMA
    assert set(p['producer']['files']) == set(rf.FILES)
    for name, value in p['producer']['files'].items():
        assert rf.digest(rf.ROOT/name) == value, name
    for name, value in p['inputs'].items():
        assert rf.digest(rf.ROOT/name) == value, name
    assert not constants_changed(p['producer']['constants'])
    assert p['design'] == json.loads(json.dumps(rf.DESIGN))
    assert all(p['checks'].values())


def test_a_uniform_plate_reproduces_the_attitude_studys_turn_times():
    rows = rf.turn_check(rf.load(rf.ATTITUDE))
    assert {'5_g_m2', '26_g_m2', '62.7_g_m2'} <= set(rows)
    for row in rows.values():
        assert row['over_plate'] == pytest.approx([1.]*len(row['over_plate']), abs=1e-3)


def test_the_shingle_gate_reproduces_the_tile_studys_lost_guarantee():
    clock = rf.attitude_clock(rf.load(rf.ATTITUDE), rf.load(rf.TILE))
    for row in clock['guarantee_loss_check']:
        assert row['screen_s'] == pytest.approx(row['product_s'], rel=1e-2)
    # The 1.5 degree shingle at an 8.5 km pitch stands 222.5 m off its neighbour, as the tile study has it.
    assert clock['angles']['envelope_0_m']['shingle_gap_m'] == pytest.approx(222.50, abs=.01)


def test_a_light_tile_reaches_its_neighbour_within_hours_and_a_dumped_store_within_minutes():
    clock = rf.attitude_clock(rf.load(rf.ATTITUDE), rf.load(rf.TILE))
    light = clock['by_tile']['annulus_allocation']['envelope_0_m']
    assert max(light['drift_hours_median_p95_torque']['contact_rad']) < 2.
    assert light['store_dump_seconds']['chosen']['contact_rad'] < 600.
    assert light['wrong_command_minutes']['contact_rad'] < 60.
    # Drift goes as the square root of inertia over torque.
    heavy = clock['by_tile']['kept_run']
    ratio = math.sqrt(heavy['inertia_kg_m2']/clock['by_tile']['annulus_allocation']['inertia_kg_m2']
                      * clock['by_tile']['annulus_allocation']['torque_N_m'][0]/heavy['torque_N_m'][0])
    assert heavy['envelope_0_m']['drift_hours_median_p95_torque']['gate_rad'][0] == pytest.approx(
        light['drift_hours_median_p95_torque']['gate_rad'][0]*ratio, rel=1e-9)


def test_free_rings_scale_to_the_design_loading_as_the_inverse_square_root_of_acceleration():
    keep = rf.keeping_clock(rf.load(rf.PHOTON), rf.load(rf.BUNDLE), rf.load(rf.TILE), rf.load(rf.KEPT))
    assert keep['free_rings_product_h']['edge_h'] < keep['free_rings_product_h']['interior_h']
    scale = math.sqrt(rf.DESIGN['areal_masses_kg_m2']['annulus_allocation']/keep['run']['areal_mass_kg_m2'])
    assert keep['free_rings_design_h']['edge_h'] == pytest.approx(keep['free_rings_product_h']['edge_h']*scale)
    # The keeping demand at the design loading lies inside the tile study's span of residual accelerations, so the
    # failed-tile times are interpolated, not extrapolated.
    for a in keep['keeping_demand_m_s2']['design']:
        assert 1e-6 <= a <= 1e-4
    single = keep['single_tile_design']
    assert single['residual_peak']['contact_min'] < single['residual_median']['contact_min']
    # The kept run's cost per orbit, spread over the orbit, agrees with the photon study's median demand.
    assert .5 < keep['kept_cost']['over_photon_median_demand'][0] < 3.


def test_a_rolled_tile_departs_faster_than_a_radial_facing_one():
    clock = rf.attitude_clock(rf.load(rf.ATTITUDE), rf.load(rf.TILE))
    for name, row in clock['rolled_law'].items():
        radial = clock['by_tile'][name]['envelope_0_m']['drift_hours_median_p95_torque']
        assert row['hours']['contact_rad'][0] < radial['contact_rad'][0]
        assert min(row['missing_torque_N_m']) > max(clock['by_tile'][name]['torque_N_m'])


def test_threshold_counts_follow_the_margin_the_cover_and_the_dispatch_rate():
    tile, uv = rf.load(rf.TILE), rf.load(rf.UV)
    power = rf.load(rf.POWER)
    shade = rf.shade_threshold(tile, uv, power, dict(whole_tile=power['ring']['transmitting_tiles']))
    n, single, rate = shade['tiles'], shade['single_cover_share'], shade['dispatch_per_day']
    for row in shade['latch'].values():
        assert row['tiles_contiguous'] == pytest.approx([m*n for m in row['margin_over_standard']])
        assert row['tiles_spread'] == pytest.approx([c/single for c in row['tiles_contiguous']])
        assert row['dispatch_days'][:2] == pytest.approx([c/rate for c in row['tiles_contiguous']])
    # The behind-the-edge latch is the tighter one; the standard level already takes up to a fifth of its transmission.
    assert max(shade['latch']['edge_200_nm']['tiles_contiguous']) < min(shade['latch']['titania']['tiles_spread'])
    assert .1 < shade['standard_level'][1]/rf.DESIGN['latch_transmission']['edge_200_nm'][0] < .2


def test_whole_tile_transmitters_alone_cannot_open_the_shade_past_the_latch():
    p = product() if rf.OUT.exists() else None
    if p is None:
        pytest.skip('ring_fault.json not written')
    lowest = min(min(r['margin_over_standard']) for r in p['shade']['latch'].values())
    assert max(p['shade']['transmitters_all_lost_uncovered']['whole_tile_collateral_0']) < lowest
    # Spread as tenth-tile collectors over ten times the tiles, they can.
    assert max(p['shade']['transmitters_all_lost_uncovered']['tenth_tile_collateral_0']) > lowest


def test_radius_neighbours_at_the_nodes_come_from_opposite_sides_of_the_stack():
    layout = rf.load(rf.LAYOUT)
    nodes = rf.nodes(layout)
    speeds = nodes['radius_neighbours']['speed_m_s']
    assert speeds['median'] <= speeds['p90'] <= speeds['max'] <= nodes['widest_pair_m_s']
    assert .5 < nodes['radius_neighbours']['opposite_side_share'] < .7
    # The widest pair is twice the orbital speed times the sine of half the stack's tilt span.
    v = math.sqrt(K.MOON_GM/(1e3*rf.DESIGN['design_radius_km']))
    span = math.radians(nodes['tilt_range_deg'][1]-nodes['tilt_range_deg'][0])
    assert nodes['widest_pair_m_s'] == pytest.approx(2*v*math.sin(span/2), rel=.06)
    # Strip neighbours differ by one 9.29 km height step and pass at a fraction of a metre a second.
    assert nodes['strip_neighbours_speed_m_s']['max'] < 1.
    # The matched layout's flown rings fall in the same radius order as their heights' ideal radii.
    flown = layout['flown']['matched']['rings']
    table = layout['stack']['ideal_radius_km']['by_height']
    ideal = np.interp([r['height_km'] for r in flown], [t['height_km'] for t in table],
                      [t['radius_km'] for t in table])
    assert list(np.argsort(ideal)) == list(np.argsort([r['radius_km'] for r in flown]))


def test_the_recoil_is_the_radio_power_over_the_speed_of_light():
    power = rf.load(rf.POWER)
    recoil = rf.recoil(power, rf.load(rf.ATTITUDE), rf.load(rf.BUNDLE))
    ring = power['ring']
    for e in (0, 1):
        assert recoil['in_use_N'][e] == pytest.approx(ring['rf_in_use_tw'][e]*1e12/ring['transmitting_tiles'][e]
                                                      / K.SPEED_OF_LIGHT)
        assert recoil['rated_N'][e] == pytest.approx(ring['tile']['whole_tile_rf_gw'][e]*1e9/K.SPEED_OF_LIGHT)
    # Sunlight absorbed over a whole tile pushes with about 444 N; the beams' step is a few to a fifth of it.
    assert recoil['whole_tile_absorbed_N'] == pytest.approx(444., rel=1e-2)
    assert .04 < min(recoil['step_share_of_absorbed_push']['in_use'])
    assert max(recoil['step_share_of_absorbed_push']['rated']) < .25


def test_stores_follow_the_restart_times_and_the_grids_black_start():
    power = rf.load(rf.POWER)
    beams = rf.beams(power, rf.FAULTS)
    week = power['reserves']['ring_fault']
    total = week['critical_total_tw']
    # A week of the critical loads is the power screen's 425-1,270 TWh.
    assert [t*24*week['days'] for t in total] == pytest.approx(week['energy_twh'], rel=2e-3)
    grid = rf.DESIGN['grid_restoration_days'][1]
    for f in rf.FAULTS:
        for side, key in (('without', 'without'), ('with_', 'with')):
            s = beams['stores'][f['id']][key]
            stopped = 'beams' in f['reach'][side] and f['restart_days'][side][1] > 0
            assert s['days'] == pytest.approx(f['restart_days'][side][1]+grid if stopped else 0.)
            assert s['energy_twh'] == pytest.approx([t*24*s['days'] for t in total])
            assert s['within_a_week'] == (s['days'] <= week['days'])
        assert beams['stores'][f['id']]['with']['days'] <= beams['stores'][f['id']]['without']['days']
    assert beams['towns_weeks_on_reserve_lift'] == pytest.approx([1/s for s in week['towns']['share_of_reserve_lift']])
    assert beams['local_share_of_ground']['flows'][0] < beams['local_share_of_ground']['flows'][1] < 1.


def test_the_latch_leads_to_losses_with_cycle_times_of_millions_to_billions_of_years():
    swollen = rf.swollen_loss(rf.load(rf.LOSS))
    year = K.JULIAN_YEAR_DAYS*K.JULIAN_DAY
    for row in swollen.values():
        assert row['cycle_years'] == pytest.approx([rf.DESIGN['atmosphere_kg']/r/year for r in row['kg_s']])
    # The design's budgets of 1-100 kg/s keep cycles of 1-100 billion years; these are far shorter.
    assert max(max(r['cycle_years']) for r in swollen.values()) < 1e9


def test_a_whole_tile_collector_shares_the_rings_loading_only_near_two_tenths_of_a_kilogram():
    recoil = rf.recoil(rf.load(rf.POWER), rf.load(rf.ATTITUDE), rf.load(rf.BUNDLE))
    assert .15 < recoil['matched_sandwich_kg_m2'] < .3
    # The 1978 reference's 10 kg/kW transmitter is more than ten times heavier than that.
    assert min(recoil['reference_over_matched']) > 10.


def test_the_halloween_storms_lost_satellites_at_about_the_loosest_budget():
    tile, uv, power = rf.load(rf.TILE), rf.load(rf.UV), rf.load(rf.POWER)
    shade = rf.shade_threshold(tile, uv, power, dict(whole_tile=power['ring']['transmitting_tiles']))
    t = rf.targets(shade)
    assert t['halloween_2003']['loss_share'] == pytest.approx([2/1000, 2/800])
    assert min(t['halloween_2003']['loss_share_over_loosest_budget']) > 1.
    assert shade['fleet_rebuild_years'] == pytest.approx(shade['tiles']/shade['dispatch_per_day']/K.JULIAN_YEAR_DAYS)


def test_fault_classes_are_complete_and_the_separation_never_widens_a_reach():
    assert len(rf.FAULTS) >= 5
    for f in rf.FAULTS:
        assert f['kind'] in rf.LABELS
        assert set(f['modes']) == set(rf.FUNCTIONS)
        assert set(f['reach']['with_']) <= set(f['reach']['without']) <= set(rf.FUNCTIONS)
        lo, hi = f['restart_days']['without']
        assert 0 <= lo <= hi
        assert f['restart_days']['with_'][1] <= hi
        for function in rf.FUNCTIONS:
            assert (f['modes'][function] != 'unaffected') == (function in f['reach']['without']), (f['id'], function)
    # Only physical faults reach the film itself.
    assert {f['id'] for f in rf.FAULTS if 'film' in f['reach']['without']} == {'F8', 'F9'}
    # With the separation the faults left that can open the shade are space weather, shared navigation and timing,
    # and the physical ones.
    table = {row['id']: row for row in rf.fault_table()}
    assert {k for k, row in table.items() if row['opens_shade_with']} == {'F6', 'F7', 'F8', 'F9'}
    assert all(row['opens_shade_without'] for k, row in table.items() if k not in ('F1', 'F2'))


def test_a_meteoroid_storm_holes_far_less_than_the_latch_allows():
    tile, uv, power = rf.load(rf.TILE), rf.load(rf.UV), rf.load(rf.POWER)
    shade = rf.shade_threshold(tile, uv, power, dict(whole_tile=power['ring']['transmitting_tiles']))
    storm = rf.meteoroid_storm(uv, shade)
    assert storm['holing_per_year_at_largest_holes'] == pytest.approx(16*uv['holing_rate_per_year_hole_factor_1'])
    assert max(storm['over_smallest_margin']) < .01


def test_the_upper_airs_bracket_is_heat_capacity_times_warming_over_heat():
    tile, uv, power = rf.load(rf.TILE), rf.load(rf.UV), rf.load(rf.POWER)
    shade = rf.shade_threshold(tile, uv, power, dict(whole_tile=power['ring']['transmitting_tiles']))
    air = rf.upper_air(rf.load(rf.LOSS), shade)
    assert air['column_kg_m2'] == pytest.approx(.3/K.MOON_SURFACE_GRAVITY)
    fastest = air['warmed_capacity_J_m2_K'][0]*rf.DESIGN['warming_K'][0]/air['threshold_heat_W_m2'][1]
    assert air['days'][0] == pytest.approx(fastest/K.JULIAN_DAY)
    # On these assumptions the air takes weeks at the least, longer than the shade's mechanical clocks of hours.
    assert air['days'][0] > 7.


def test_targets_and_the_trim_at_stake():
    t = rf.targets(dict(tiles=28233094, latch=dict(x=dict(tiles_contiguous=[30000., 60000.]))))
    assert t['per_year']['10000_years'] == pytest.approx(rf.DESIGN['accepted_probability']/1e4)
    assert t['aviation_per_year_of_continuous_operation'] == pytest.approx(1e-9*365.25*24)
    assert t['uncontrolled_at_once_budget_tiles']['x'] == pytest.approx([10000., 6000.])
    trim = rf.trim_at_stake()
    d = rf.DESIGN['trim']
    assert trim['sunlight_share_if_transmitted'] == pytest.approx(
        [d['range']*d['reflected'][0]*d['cover'][0], d['range']*d['reflected'][1]*d['cover'][1]])
    assert trim['uv_transmission_at_stake'] == 0.


def test_every_precedent_is_in_the_sources_file_with_its_reading():
    sources = {s['id']: s for s in rf.load(rf.SOURCES)['sources']}
    for f in rf.FAULTS:
        for ref in f['precedents']:
            assert ref in sources, (f['id'], ref)
    for s in sources.values():
        assert s['read'] and s['used_for'] and s['citation']
