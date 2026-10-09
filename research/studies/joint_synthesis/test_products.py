"""The joint synthesis's products bind the files and inputs that made them."""
import hashlib
import json
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PRODUCTS = ('night_sky', 'ledgers', 'magnets_map', 'places', 'upper_air_no', 'air_chemistry')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


@pytest.mark.parametrize('name', PRODUCTS)
def test_product_binds_producer_and_inputs(name):
    d = json.loads((HERE / 'results' / f'{name}.json').read_text())
    assert d['schema'].startswith('terluna.research.joint-synthesis')
    for path, h in d['producer']['files'].items():
        assert digest(ROOT / path) == h, path
    for path, h in d['producer'].get('inputs', {}).items():
        assert digest(ROOT / path) == h, path


def test_loss_rates_carry_cycle_times():
    air = json.loads((HERE / 'results' / 'ledgers.json').read_text())['air']
    for row in air['losses'].values():
        assert row['cycle_years']


def test_ledger_powers_are_in_watts():
    energy = json.loads((HERE / 'results' / 'ledgers.json').read_text())['energy']
    assert 1e11 < min(energy['film_plant_w']) and max(energy['film_plant_w']) < 1e13      # 0.15-2.9 TW
    assert 6e16 < min(energy['tiles_absorb_w']) and max(energy['tiles_absorb_w']) < 9e16  # 65-83 PW


def test_earthlight_follows_the_earths_height():
    night = json.loads((HERE / 'results' / 'night_sky.json').read_text())['places']
    places = json.loads((HERE / 'results' / 'places.json').read_text())['places']
    centre, far = night['nearside, 0 deg']['diffuse_rho_0.001'], night['far side, 0 deg']['diffuse_rho_0.001']
    russell = places['Oceanus Procellarum by Russell']['night']
    smythii = places['eastern Smythii headland']['night']
    # Under the full Earth overhead the night stays above practical dusk; toward the limb the Earth sinks.
    assert centre['natural_darkest_lux'] > 2.98 > russell['natural_midnight_lux'] > smythii['natural_midnight_lux']
    assert smythii['natural_midnight_lux'] > 10 * far['natural_midnight_lux']
    assert places['eastern Smythii headland']['side'] == 'limb'


def test_earthlight_evaluator_matches_the_dated_calendar():
    import numpy as np
    from illumination.calendar.validate import Response
    from shared.constants import EARTH_MOON_DISTANCE
    transfer = json.loads((ROOT / 'illumination' / 'calendar' / 'results' / 'transfer.json').read_text())
    response = Response(transfer)
    calendar = json.loads((ROOT / 'research' / 'studies' / 'sea_appearance' / 'results' / 'lighting_calendar.json').read_text())
    for site in calendar['sites']:
        s = site['series']
        for i in range(0, len(s['earth_elevation_deg']), 30):
            if s['earth_elevation_deg'][i] < 5 or s['ground_lux_from_earth'][i] < 0.1:
                continue
            phase = float(np.degrees(np.arccos(2 * s['earth_lit_fraction'][i] - 1)))
            g = dict(sun_elevation_deg=s['sun_elevation_deg'][i], sun_radius_deg=0.27, sun_distance_factor=1.0,
                     earth_phase_deg=phase, earth_elevation_deg=s['earth_elevation_deg'][i], earth_radius_deg=0.95,
                     earth_bright_limb_rad=0.0, earth_distance_m=EARTH_MOON_DISTANCE, earth_sunlight_factor=1.0,
                     source_separation_deg=90.0)
            # The dated month's Earth distance and phase law set its light above the air; compare the path through it.
            above = np.interp(phase, transfer['phase_deg'], transfer['earth']['above_air_lux'])
            ground = response.light(g)['earth'] * s['earth_above_air_lux'][i] / above
            assert ground == pytest.approx(s['ground_lux_from_earth'][i], rel=0.08), site['short']
