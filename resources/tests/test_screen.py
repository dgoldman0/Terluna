"""The resources screen binds its inputs and its arithmetic stays checkable."""
import hashlib
import json
from pathlib import Path

import numpy as np

from resources import screen
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(screen.OUT.read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def test_product_binds_producer_inputs_and_constants():
    assert PRODUCT['schema'] == screen.SCHEMA
    for path, h in {**PRODUCT['producer']['files'], **PRODUCT['producer']['inputs']}.items():
        assert digest(ROOT / path) == h, path
    assert not constants_changed(PRODUCT['producer']['constants'])


def test_the_oxygens_hydrogen_would_burn_all_of_it():
    h = PRODUCT['hydrogen']
    assert abs(h['share_of_air_oxygen_it_would_burn'] - 1.0) < 1e-3   # O2 + 2 H2 -> 2 H2O
    assert 0.25 < h['mole_fraction_if_released'] < 0.27


def test_fall_energy_is_the_gravity_well():
    gm_r = screen.MOON_GM / screen.MOON_RADIUS / 1e6
    fall = PRODUCT['delivery']['fall_mj_per_kg']
    assert np.isclose(fall['far_terminal'], gm_r * (1 - screen.MOON_RADIUS / screen.TERMINAL_RADIUS_M), atol=0.01)
    assert all(1.0 < x < gm_r for x in fall['exobase'])
    assert 1.5 < PRODUCT['delivery']['exobase_radius_r'][0] < PRODUCT['delivery']['exobase_radius_r'][1] < 4.0


def test_film_titania_and_silica_fit_inside_the_fleets_mass():
    film = PRODUCT['film']
    for tio2, sio2, total in zip(film['tio2_gt'], film['sio2_gt'], film['optical_mass_gt']):
        assert tio2 + sio2 < total
        assert 0.05 < tio2 / (tio2 + sio2) < 0.15
    for recovery, feeds in film['fresh_feed_kg_s'].items():
        for feed, upkeep in zip(feeds, film['upkeep_kg_s']):
            assert np.isclose(feed, upkeep * (1 - float(recovery)), rtol=0.01)


def test_rock_takes_oxygen_far_faster_than_space():
    sink = PRODUCT['oxygen_sink']
    assert all(t > 10 for t in sink['times_all_escape'])
    assert all(1 < c < 1000 for c in sink['o2_cycle_myr'])


def test_seas_freshen_slowly():
    seas = PRODUCT['fresh_seas']
    for ions, everything in zip(seas['years_to_1_g_kg']['cations'], seas['years_to_1_g_kg']['all_dissolved']):
        assert everything < ions
    assert min(seas['years_to_1_g_kg']['all_dissolved']) > 500
