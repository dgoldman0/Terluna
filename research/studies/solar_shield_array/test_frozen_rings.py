import json

import numpy as np

from shared import constants as K
from shared.provenance import constants_changed
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from . import frozen_rings as fr

PRODUCT = json.loads(fr.OUT.read_text())
FIRST, LAST = PRODUCT['passes'][0], PRODUCT['passes'][-1]


def key(row):
    return row['plane'], row['radius_km'], row['tilt_deg'], row['sigma_kg_m2']


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.frozen-rings/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(fr.FILES)
    for name, value in producer['files'].items():
        assert fr.digest(fr.ROOT/name) == value, name
    if DEFAULT_KERNEL.exists():
        assert producer['inputs']['de440s.bsp'] == fr.digest(DEFAULT_KERNEL)
    assert not constants_changed(producer['constants'])


def test_a_ring_starts_with_the_eccentricity_asked_for():
    sun = np.array([1.5e11, 0., 0.])
    pole = np.array([0., 0., 1.])
    state, up, along = fr.start(sun, pole, 15000., 20., -.05, .01)
    e = fr.eccentricity(state[None])[0]
    assert np.isclose(e@up, -.05) and np.isclose(e@along, .01)
    assert np.isclose(np.linalg.norm(state[:3]), 15e6)
    # Retrograde about the pole.
    assert np.cross(state[:3], state[3:])@pole < 0


def test_the_forced_eccentricity_points_away_from_the_sun_and_grows_as_tiles_get_lighter():
    rows = {key(r): r for r in LAST}
    for (plane, radius, tilt, sigma), row in rows.items():
        cx, cy = row['centre']
        assert cx < 0 and abs(np.degrees(np.arctan2(cy, cx))) > 150.
    for plane in fr.DESIGN['planes']:
        for radius in fr.DESIGN['radii_km']:
            for tilt in fr.DESIGN['tilts_deg']:
                sizes = [np.hypot(*rows[(plane, radius, tilt, s)]['centre'])
                         for s in sorted(fr.DESIGN['areal_masses_kg_m2'], reverse=True) if (plane, radius, tilt, s) in rows]
                assert all(b > a for a, b in zip(sizes, sizes[1:]))


def test_starting_on_the_forced_eccentricity_shrinks_the_swing_of_heavy_tiles():
    first = {key(r): r for r in FIRST}
    for row in LAST:
        if row['sigma_kg_m2'] == max(fr.DESIGN['areal_masses_kg_m2']):
            assert row['circle_radius'] < .7*first[key(row)]['circle_radius']


def test_light_tiles_dip_inside_four_lunar_radii_and_tilted_rings_at_15000_km_lag_the_sun():
    light = [r for r in LAST if r['sigma_kg_m2'] == min(fr.DESIGN['areal_masses_kg_m2'])]
    assert light and min(r['periapsis_min_km'] for r in light) < 4*K.MOON_RADIUS/1e3
    lagging = [r for r in LAST if r['radius_km'] == 15000. and r['tilt_deg'] > 0]
    assert lagging and all(abs(r['height_drift_km_per_day']) > 10. for r in lagging)
