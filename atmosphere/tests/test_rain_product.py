"""The rain product binds its code and inputs, reproduces from them, and carries the CM1 check consistently."""
import hashlib
import json

import numpy as np
import pytest

from shared.constants import MOON_SURFACE_GRAVITY
from shared.provenance import constants_changed
from atmosphere.rain import drops, product as rain, targets as climate_targets

PRODUCT = json.loads(rain.OUT.read_text())
DATA, BINDING = climate_targets.load()
READ = climate_targets.read(DATA)
WORLDS = ('moon', 'earth', 'moon_air_earth_gravity')


def test_header_and_bindings():
    """The code and the in-tree inputs by hash, the climate targets among them as targets.py binds them."""
    assert PRODUCT['schema'] == rain.SCHEMA
    assert PRODUCT['evidence'] and PRODUCT['reading_rule'] and PRODUCT['assumptions']
    producer = PRODUCT['producer']
    assert set(producer['files']) == set(rain.CODE)
    for kind, entries in BINDING.items():
        assert producer[kind] == entries, kind
    for path, digest in {**producer['files'], **producer.get('inputs', {})}.items():
        assert hashlib.sha256((rain.ROOT / path).read_bytes()).hexdigest() == digest, path
    assert not constants_changed(producer['constants'])
    assert DATA['schema'] == climate_targets.SCHEMA


def test_targets_air_reproduces_the_climate_targets():
    for key in rain.TARGETS:
        air = PRODUCT['targets'][key]
        climate = READ[key]['climate']
        assert air['air_2m_k'] == pytest.approx(climate['land_air_2m_median_c'] + 273.15, abs=1e-9)
        assert rain.wet_bulb(air['air_2m_k'], air['specific_humidity'], rain.MOON_PRESSURE_PA) == pytest.approx(
            climate['wet_bulb_2m_median_c'] + 273.15, abs=1e-6)
        assert 0.5 < air['relative_humidity'] < 0.75
    assert [PRODUCT['targets'][k]['air_2m_k'] for k in rain.TARGETS] == sorted(
        (PRODUCT['targets'][k]['air_2m_k'] for k in rain.TARGETS), reverse=True)
    # The resort is the resort band's candidate at 301 K; the lock is the design case and the hot world the test.
    row = next(r for r in DATA['resort_band']['rows'] if r['ground_k'] == 301.0)
    resort = PRODUCT['targets']['resort']
    assert (resort['held_in'], resort['ground_k'], resort['share_of_the_way']) == (
        'resort_band.rows', 301.0, row['share_of_the_way'])
    assert resort['air_2m_k'] == row['comfort']['land_air_2m_median_c'] + 273.15
    assert resort['wet_bulb_2m_k'] == row['comfort']['wet_bulb_2m_median_c'] + 273.15
    assert PRODUCT['targets']['lock']['share_of_the_way'] == 0.0
    assert PRODUCT['targets']['hot_world']['share_of_the_way'] == 1.0


def test_values_reproduce_from_the_code():
    air = rain.target_air(READ)
    assert air == PRODUCT['targets']
    columns = rain.columns(air['resort'])
    speeds = PRODUCT['fall_speed']['resort']
    i = PRODUCT['diameters_mm'].index(2.0)
    for world in WORLDS:
        assert drops.terminal_speed(2e-3, columns[world]) == pytest.approx(speeds[world]['terminal_speed_m_s'][i],
                                                                           rel=1e-12)
    moon = columns['moon']
    assert moon.cloud_base_m == pytest.approx(PRODUCT['survival']['resort']['moon']['cloud_base_m'], rel=1e-12)
    row = next(r for r in PRODUCT['survival']['resort']['moon']['drops'] if r['initial_diameter_mm'] == 3.0)
    again = drops.fall(3e-3, moon, sample_heights_m=[moon.cloud_base_m, 0.0])
    assert again['surface_diameter_m'] * 1e3 == pytest.approx(row['surface_diameter_mm'], rel=1e-9)


def test_fall_speeds_and_survival_hold_together():
    for key in rain.TARGETS:
        speeds = PRODUCT['fall_speed'][key]
        for world in WORLDS:
            v = [x for x in speeds[world]['terminal_speed_m_s'] if x is not None]
            assert np.all(np.diff(v) > 0)
        assert all(m < e for m, e in zip(speeds['moon']['terminal_speed_m_s'][:12],
                                         speeds['earth']['terminal_speed_m_s'][:12]))
        assert speeds['moon']['stable_diameter_mm'] > 2.4 * speeds['earth']['stable_diameter_mm']
        moon, earth = PRODUCT['survival'][key]['moon'], PRODUCT['survival'][key]['earth']
        assert moon['cloud_base_m'] > 5 * earth['cloud_base_m']
        assert moon['smallest_surviving_mm'] > 2 * earth['smallest_surviving_mm']
        for side in (moon, earth):
            survived = [r['initial_diameter_mm'] for r in side['drops'] if r['status'] == 'surface']
            assert survived and min(survived) >= side['smallest_surviving_mm']
            assert all(r['status'] == 'surface' for r in side['drops'] if r['initial_diameter_mm'] > side[
                'smallest_surviving_mm'])
            shares = {(s['dm_mm'], s['shape_mu']): s['water_survival_fraction'] for s in side['rain']['spectra']}
            assert shares[(0.75, 0.0)] < shares[(1.5, 0.0)] < shares[(2.5, 0.0)] <= 1
        assert (moon['rain']['central_encounters']['tau_by_surface_water']
                > 10 * earth['rain']['central_encounters']['tau_by_surface_water'])


def test_spectra_hold_on_a_finer_grid():
    for world in ('moon', 'earth'):
        coarse = PRODUCT['survival']['resort'][world]['rain']['spectra']
        fine = PRODUCT['grid_check']['resort'][world]['spectra']
        assert [(c['dm_mm'], c['shape_mu']) for c in coarse] == [(f['dm_mm'], f['shape_mu']) for f in fine]
        for c, f in zip(coarse, fine):
            assert abs(c['water_survival_fraction'] - f['water_survival_fraction']) < 5e-4


def test_column_sensitivities():
    """The column's dry adiabat takes the repository's heat capacity; the parcel started at the GCM's lowest level and
    the transport model's heat capacity, the recovered solver's, move the cloud base by tens of metres. The GCM's mean
    air over land is more stable than the column and implies about 1.158 atm at the ground."""
    for key in rain.TARGETS:
        sens = PRODUCT['column_sensitivity'][key]
        moon = PRODUCT['survival'][key]['moon']
        base = moon['cloud_base_m']
        lowest = sens['lowest_level']
        assert 850 < lowest['height_m'] < 950 and -30 < lowest['cloud_base_change_m'] < 0
        assert lowest['cloud_base_m'] == pytest.approx(base + lowest['cloud_base_change_m'], abs=1e-6)
        cp = sens['transport_heat_capacity']
        assert MOON_SURFACE_GRAVITY / cp['column_heat_capacity_j_kg_k'] * 1e3 == pytest.approx(moon['lapse_k_km'])
        assert cp['heat_capacity_j_kg_k'] / cp['column_heat_capacity_j_kg_k'] - 1 == pytest.approx(0.0083, abs=5e-4)
        assert 0.009 < cp['cloud_base_change_m'] / base < 0.012
        assert 0 < cp['smallest_surviving_change'] < 0.01
        air = sens['gcm_land_air']
        assert air['height_km'] == [0.9, 2.0, 5.0]
        assert 0.7 < air['gcm_lapse_k_km'] < 1.0 < moon['lapse_k_km']
        assert 2.5 < air['gcm_air_k'][2] - air['column_air_k'][2] < 4.5
        assert air['ground_pressure_atm'] == pytest.approx(1.158, abs=0.002)


def test_cm1_check():
    factor = (9.81 / MOON_SURFACE_GRAVITY)**0.6
    for key in rain.TARGETS:
        check = PRODUCT['cm1_check'][key]
        assert check['cm1_slowing_at_every_size'] == pytest.approx(factor, rel=1e-12)
        assert 1.0 < check['crossing_diameter_mm'] < 1.5
        slowing = check['slowing_by_gravity']
        assert np.all(np.diff(slowing) < 0)
        mm = check['diameter_mm']
        assert slowing[mm.index(0.2)] > 4.3 and slowing[mm.index(6.0)] < 2.0
        assert np.all(np.array(check['rigid_sphere_slowing_by_gravity']) >= np.array(slowing) - 1e-9)
        ratio = np.array(check['cm1_over_solver'])
        assert ratio[mm.index(0.2)] > 1.5 and ratio[mm.index(2.0)] < 0.9
        assert check['cm1_mean_fall_speed_cap_m_s'] == pytest.approx(
            np.interp(np.log(check['solver_diameter_at_the_cap_mm']), np.log(mm), check['solver_speed_m_s']),
            rel=0.02)
        loss = check['water_lost_per_metre_cm1_over_solver']
        assert loss[mm.index(0.2)] < 0.7 and loss[mm.index(2.0)] > 1.05
