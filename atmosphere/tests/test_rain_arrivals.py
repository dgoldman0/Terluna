"""The arrivals product binds its code and inputs, agrees with the rain product, and closes its water ledgers."""
import hashlib
import json
import math

import numpy as np
import pytest

from shared.provenance import constants_changed
from atmosphere.rain import arrivals, drops, product as rain, targets as climate_targets

PRODUCT = json.loads(arrivals.OUT.read_text())
RAIN = json.loads(rain.OUT.read_text())


def test_header_and_bindings():
    assert PRODUCT['schema'] == arrivals.SCHEMA
    assert PRODUCT['evidence'] and PRODUCT['reading_rule'] and PRODUCT['assumptions']
    producer = PRODUCT['producer']
    assert set(producer['files']) == set(arrivals.CODE)
    for kind, entries in climate_targets.load()[1].items():
        assert producer[kind] == entries, kind
    for path, digest in {**producer['files'], **producer.get('inputs', {})}.items():
        assert hashlib.sha256((arrivals.ROOT / path).read_bytes()).hexdigest() == digest, path
    assert not constants_changed(producer['constants'])
    assert [(s['dm_mm'], s['shape_mu']) for s in PRODUCT['spectra']] == list(rain.SPECTRA)


def test_survival_shares_are_the_rain_products():
    for key in rain.TARGETS:
        for world in arrivals.WORLDS:
            ours = PRODUCT['targets'][key][world]
            theirs = RAIN['survival'][key][world]
            assert ours['air']['cloud_base_m'] == pytest.approx(theirs['cloud_base_m'], rel=1e-12)
            shares = {(s['dm_mm'], s['shape_mu']): s['water_survival_fraction'] for s in theirs['rain']['spectra']}
            for row in ours['spectra']:
                assert row['water_survival_fraction'] == pytest.approx(shares[(row['dm_mm'], row['shape_mu'])],
                                                                       rel=1e-9, abs=1e-15)


def test_ground_fluxes_close_the_water_ledger():
    for key in rain.TARGETS:
        for world in arrivals.WORLDS:
            side = PRODUCT['targets'][key][world]
            size = np.array([0.0 if x is None else x for x in side['classes']['surface_diameter_mm']]) * 1e-3
            reached = np.array([s == 'surface' for s in side['classes']['status']])
            for row in side['spectra']:
                if row['ground_flux_per_mm_h'] is None:
                    assert not row['water_survival_fraction']
                    continue
                flux = np.array(row['ground_flux_per_mm_h'])
                assert np.all(flux[~reached] == 0) and np.all(flux >= 0)
                assert np.sum(math.pi / 6 * size**3 * flux) * arrivals.MM_H_PER_M_S == pytest.approx(1.0, rel=1e-12)
                assert row['cloud_base_rain_per_ground_rain'] * row['water_survival_fraction'] == pytest.approx(1.0)
                assert row['arrivals_m2_s_per_mm_h'] == pytest.approx(flux.sum(), rel=1e-12)


def test_the_moons_rain_arrives_slower_and_sparser():
    central = (PRODUCT['central_spectrum']['dm_mm'], PRODUCT['central_spectrum']['shape_mu'])
    for key in rain.TARGETS:
        rows = {w: next(s for s in PRODUCT['targets'][key][w]['spectra'] if (s['dm_mm'], s['shape_mu']) == central)
                for w in arrivals.WORLDS}
        assert rows['moon']['arrivals_m2_s_per_mm_h'] < rows['earth']['arrivals_m2_s_per_mm_h']
        assert rows['moon']['arrival_speed_m_s_by_number']['p50'] < 0.5 * rows['earth']['arrival_speed_m_s_by_number'][
            'p50']
        assert rows['moon']['water_survival_fraction'] < rows['earth']['water_survival_fraction']


def test_a_class_reproduces_from_the_code():
    air = rain.target_air(climate_targets.read(climate_targets.load()[0]))
    moon = rain.columns(air['lock'])['moon']
    side = PRODUCT['targets']['lock']['moon']['classes']
    # The largest class, far from the survival threshold: near it the adaptive integration turns the rounding of a
    # stored size into a change near the solver's tolerance.
    again = drops.fall(side['cloud_base_diameter_mm'][-1] * 1e-3, moon, sample_heights_m=[moon.cloud_base_m, 0.0])
    assert again['surface_diameter_m'] * 1e3 == pytest.approx(side['surface_diameter_mm'][-1], rel=1e-6)
    assert again['surface_speed_m_s'] == pytest.approx(side['impact_speed_m_s'][-1], rel=1e-6)


def test_drips_fall_below_both_bounds_and_faster_from_higher():
    for key in rain.TARGETS:
        for world in arrivals.WORLDS:
            side = PRODUCT['targets'][key][world]
            g = side['air']['gravity_m_s2']
            for mm in arrivals.DRIP_DIAMETERS_MM:
                rows = [r for r in side['drips'] if r['release_diameter_mm'] == mm]
                speeds = [r['impact_speed_m_s'] for r in rows]
                assert np.all(np.diff(speeds) > 0)
                for r in rows:
                    assert r['impact_speed_m_s'] < math.sqrt(2 * g * r['release_height_m'])
                    # A small drip that reaches its terminal speed and then shrinks as it evaporates lags its own
                    # falling terminal speed, by 0.03% for a 1 mm drip from 4 m on the Moon.
                    assert r['impact_speed_m_s'] < 1.001 * r['terminal_speed_m_s']
                    assert 0.98 < (r['surface_diameter_mm'] / mm)**3 <= 1.0
        moon, earth = (PRODUCT['targets'][key][w]['drips'] for w in arrivals.WORLDS)
        assert all(m['impact_speed_m_s'] < e['impact_speed_m_s'] for m, e in zip(moon, earth))


def test_the_tail_reaches_the_moons_ground_and_breaks_up_on_earth():
    for key in rain.TARGETS:
        moon, earth = (PRODUCT['targets'][key][w]['tail'] for w in arrivals.WORLDS)
        assert all(r['status'] == 'surface' and r['surviving_water_fraction'] > 0.85 for r in moon)
        assert all(r['status'] == 'above the stability screen' for r in earth)
        assert np.all(np.diff([r['impact_speed_m_s'] for r in moon]) > 0)
