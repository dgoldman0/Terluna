"""The energy-options screen: its binding and regeneration, the heat each option leaves per delivered watt, the
coupling between the beams' heat and the dimming that powers the tiles, the air's absorption by height and frequency,
the night's store, the mixes' sums, and the design run's extract it reads."""
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np
import pytest

from provisioning import energy_extract as X
from provisioning import energy_options as E
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(E.OUT.read_text())
CLIMATE = json.loads(X.OUT.read_text())
POWER = json.loads((ROOT / E.POWER).read_text())


def close(a, b, rel=2e-3):
    return math.isclose(a, b, rel_tol=rel, abs_tol=1e-9)


def test_product_binds_its_producer_inputs_sources_and_constants():
    assert PRODUCT['schema'] == E.SCHEMA
    prod = PRODUCT['producer']
    assert set(prod['files']) == set(E.CODE) and set(prod['inputs']) == set(E.INPUTS)
    for path, h in {**prod['files'], **prod['inputs'], **prod['sources']}.items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == h, path
    assert not constants_changed(prod['constants'])
    for key in ('evidence', 'reading_rule', 'assumptions'):
        assert PRODUCT[key]


def test_a_fresh_run_matches_the_product():
    assert json.loads(E.dumps(E.results())) == PRODUCT


def test_every_outside_value_is_labelled_with_its_kind_and_basis():
    kinds = {'project', 'published', 'recalled', 'estimate', 'design guess', 'published/estimate'}
    for key, entry in PRODUCT['assumptions'].items():
        assert entry['kind'] in kinds, key
        assert len(entry['basis']) > 10, key
    sources = json.loads(E.SOURCES.read_text())['sources']
    assert PRODUCT['sources_count'] == len(sources) and all(s['access'] for s in sources)


def test_the_loads_are_the_synthesis_array_case():
    loads = PRODUCT['loads']
    assert close(loads['total_tw'][0], 22.1, 5e-3) and close(loads['total_tw'][1], 105.5, 5e-3)
    for e in (0, 1):
        assert close(loads['sky_tw'][e] + loads['ground_tw'][e], loads['total_tw'][e], 1e-3)
    # The power screen's delivered total.
    for e in (0, 1):
        assert close(loads['total_tw'][e], POWER['reserves']['ring_fault']['delivered_total_tw'][e], 5e-3)


def test_each_beams_heat_parts_add_up_and_the_air_is_a_small_part():
    chain = POWER['absorption']['link']['chain']
    for key, pair in PRODUCT['beams'].items():
        for o in pair:
            assert close(sum(o['heat_parts_w_per_w'].values()), o['heat_w_per_w'], 2e-3), key
            assert close(o['tile_electricity_w_per_w'] * chain['tx'], o['heat_w_per_w'], 2e-3), key
    floor = 1.0 / (chain['capture'] * chain['rx'] * chain['bus'])
    assert close(PRODUCT['beams']['ground_no_air'][0]['heat_w_per_w'], floor, 1e-3)
    for o in PRODUCT['beams']['chosen_system']:
        assert o['heat_parts_w_per_w']['air'] < 0.13
        assert o['heat_parts_w_per_w']['air'] / o['heat_w_per_w'] < 0.09
        assert o['heat_parts_w_per_w']['air'] < o['heat_parts_w_per_w']['rectenna']


def test_the_chosen_beams_are_the_power_screens():
    link = POWER['absorption']['link']
    chosen = PRODUCT['beams']['chosen_system']
    assert close(chosen[0]['heat_w_per_w'], link['with_the_best_tiles']['heat_per_w_delivered'][0], 2e-3)
    assert close(chosen[1]['heat_w_per_w'], link['heat_per_w_delivered_with_air'][0], 2e-3)
    assert close(PRODUCT['beams']['towns_2_45_ghz'][1]['heat_w_per_w'], link['towns_heat_per_w'], 2e-3)
    assert close(PRODUCT['beams']['ground_2_45_ghz'][1]['heat_w_per_w'], link['ground_heat_per_w'], 2e-3)


def test_the_dimming_that_answers_the_beams_cannot_power_them_but_the_resorts_dimming_can():
    b = PRODUCT['basis']
    assert close(b['point']['heat_tw'], 0.01 * b['point']['absorbed_sunlight_w_m2'] * b['moon_area_m2'] / 1e12, 1e-3)
    for o in PRODUCT['beams']['chosen_system']:
        assert close(o['points'], o['heat_w_per_w'] / b['point']['heat_tw'], 2e-3)
        assert close(o['dimming_light_tw_per_tw'], o['points'] * b['point']['light_at_the_shield_tw'], 2e-3)
        for share, conv in zip(o['self_supply_share'], (0.2, 0.3)):
            assert close(share, o['dimming_light_tw_per_tw'] * conv / o['tile_electricity_w_per_w'], 2e-3)
            assert share < 0.3
    # The light the resort's whole dimming takes covers the tiles' electricity at both ends, at 20% conversion alone.
    low = b['resort_dimming_electricity_at_the_lower_conversion_tw']
    assert low[0] > b['tile_electricity_tw'][0] and low[1] > b['tile_electricity_tw'][1]


def test_the_air_takes_less_higher_up_and_frequency_barely_matters_below_10_ghz():
    for name, col in PRODUCT['air']['columns'].items():
        rows = col['rows']
        for f in E.FREQS_GHZ:
            for e in (90.0, 30.0):
                losses = [rows[f'{f:g}_ghz_{h:g}_km_{e:g}_deg']['db'] for h in E.HEIGHTS_KM]
                assert all(a > b for a, b in zip(losses, losses[1:])), (name, f, e)
        z = {f: rows[f'{f:g}_ghz_0_km_90_deg']['share_lost'] for f in E.FREQS_GHZ}
        assert 0.75 < z[1.0] / z[2.45] < 1.0, name                 # oxygen's flat relaxation spectrum
        assert rows['2.45_ghz_0_km_90_deg']['vapour_share'] < 0.06, name
    base = PRODUCT['air']['columns']['1.2 atm, 295 K']['rows']['2.45_ghz_0_km_90_deg']
    assert close(base['db'], POWER['absorption']['by_elevation']['2.45_ghz']['lowland']['90_deg']['db'], 1e-3)
    thin = PRODUCT['air']['columns']['1.0 atm, 301 K']['rows']['2.45_ghz_0_km_90_deg']['db']
    assert 0.75 < thin / base['db'] < 0.9


def test_the_store_solver_balances_the_cycle():
    hours = np.linspace(-175.0, 175.0, 36)
    shape = np.maximum(np.cos(np.radians(hours)), 0.0)
    for eta in (0.4, 0.8, 1.0):
        out = E.solve_store(shape, eta)
        s = shape / shape.mean()
        c = out['generation_over_load']
        direct = np.mean(np.minimum(c * s, 1.0))
        stored = np.mean(np.maximum(c * s - 1.0, 0.0))
        assert close(direct + eta * stored, 1.0, 1e-6)
        assert 0.45 < out['store_cycle_fraction'] < 0.7                # the night and the low sun
    flat = E.solve_store(np.ones(36), 0.5)
    assert close(flat['generation_over_load'], 1.0, 1e-6) and flat['store_share'] < 1e-9
    for band in PRODUCT['local_solar']['by_band'].values():
        for store in band['stores'].values():
            assert all(0.0 <= x < E.MONTH_H for x in store['store_twh_per_tw'])


def test_panels_take_more_of_the_grounds_sunlight_than_the_beams_dimming():
    light = PRODUCT['local_solar']['sunlight_taken_from_the_ground_w_per_w']
    assert min(light['panels_global_grid']) > 4 * max(light['beams_dimming'])
    routes = PRODUCT['local_solar']['routes']
    for name, r in routes.items():
        assert r['heat_w_per_w'][0] < r['heat_w_per_w'][1] < PRODUCT['beams']['chosen_system'][0]['heat_w_per_w'], name
        assert r['panel_km2_per_tw'][0] > 5 * PRODUCT['beams']['ground_2_45_ghz'][0]['ground_closed_km2_per_tw'], name


def test_the_mixes_add_up_and_order_by_heat():
    mixes = PRODUCT['mixes']
    loads = PRODUCT['loads']
    link = POWER['absorption']['link']
    for e in (0, 1):
        heat = loads['sky_tw'][e] * link['towns_heat_per_w'] + loads['ground_tw'][e] * link['ground_heat_per_w']
        assert close(mixes['A']['heat_tw'][e], heat, 2e-3)
        assert (mixes['B']['heat_tw'][e] < mixes['A2']['heat_tw'][e] < mixes['A']['heat_tw'][e]
                < mixes['E']['heat_tw'][e] < mixes['F']['heat_tw'][e])
        moved = loads['delivered_tw']['close_computing'][e] * PRODUCT['orbit']['computing']['moved_share']
        assert close(mixes['A2']['heat_tw'][e], heat - moved * link['ground_heat_per_w'], 2e-3)
        for m in mixes.values():
            assert close(m['points'][e], m['heat_tw'][e] / PRODUCT['basis']['point']['heat_tw'], 2e-3)
    assert close(mixes['A']['points'][0], PRODUCT['basis']['beams_now_points'][0], 2e-2)
    assert close(mixes['A']['points'][1], PRODUCT['basis']['beams_now_points'][1], 2e-2)
    assert mixes['C']['ground_km2'][0] > 5 * mixes['A']['ground_km2'][0]


def test_the_winds_are_low_in_the_gcm_and_in_the_cloud_resolving_rings():
    wind = PRODUCT['other_flows']['wind']
    earth = PRODUCT['assumptions']['earth_land_wind_power_density_w_m2']['value'][0]
    assert wind['share_below_cut_in'] > 0.8 and wind['power_density_w_m2'] < 0.05 * earth
    for name, ring in wind['cm1_rings'].items():
        for where in ('land_10_m', 'water_10_m', '1_km'):
            row = ring[where]
            assert row['power_density_bound_w_m2'] < 0.4 * earth, (name, where)
            if row['power_density_w_m2'] is not None:
                assert row['power_density_w_m2'] <= row['power_density_bound_w_m2'], (name, where)


def test_the_extract_reads_the_design_run():
    assert CLIMATE['schema'] == X.SCHEMA and CLIMATE['run']['name'] == X.RUN
    assert CLIMATE['run']['years'] == list(X.YEARS) and CLIMATE['run']['records'] == 120 * len(X.YEARS)
    assert CLIMATE['producer']['files']['provisioning/energy_extract.py'] == hashlib.sha256(
        (ROOT / 'provisioning/energy_extract.py').read_bytes()).hexdigest()
    whole = CLIMATE['whole']
    assert whole['ghi_w_m2'] < whole['ghi_clear_w_m2'] < whole['toa_w_m2']
    for band in CLIMATE['bands'].values():
        assert band['land_ghi_w_m2'] < band['land_toa_w_m2']
        assert abs(np.mean(band['day_profile_ghi_w_m2']) - band['land_ghi_w_m2']) < 0.15 * band['land_ghi_w_m2']
    # The run's files, where the research drive is mounted, are the ones the extract read (sizes only; the full
    # hashes are rechecked by rerunning the extract).
    for path, rec in CLIMATE['run']['files'].items():
        full = ROOT / path
        if not full.exists():
            pytest.skip('the design run is not on this machine')
        assert full.stat().st_size == rec['bytes'], path


@pytest.mark.skipif(not os.environ.get('TERLUNA_RERUN_EXTRACT'), reason='reads 2.3 GB of GCM output')
def test_a_fresh_extract_matches_the_product():
    pytest.importorskip('netCDF4')               # the GCM's environment: climate/gcm/.venv
    assert json.loads(X.dumps(X.results())) == CLIMATE
