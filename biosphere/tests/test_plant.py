"""Checks on the whole-plant carbon cycle: calibration, strategies, cost accounting and stored results."""
import json
import math

import numpy as np
import pytest

from biosphere.canopy import plant


class SimpleTables:
    """Photosynthesis proportional to leaf area and to the sine of the Sun, with a little twilight."""

    def photosynthesis(self, world, lai, elevation, temp):
        light = np.sin(np.radians(np.clip(elevation, 0.0, 90.0))) + 0.02 * np.clip(1 + elevation / 20.0, 0.0, 1.0)
        return 60.0 * light * np.broadcast_to(np.asarray(lai, dtype=float), elevation.shape) / 5.0

    def par(self, world, elevation):
        return 2000.0 * np.sin(np.radians(np.clip(elevation, 0.0, 90.0))) + 40.0 * np.clip(1 + elevation / 20.0, 0.0, 1.0)


def trace(day_c=23.0, night_c=20.0):
    hour = (np.arange(plant.STEPS) + 0.5) / plant.STEPS * 360.0 - 180.0
    return np.where(np.abs(hour) < 90, day_c, night_c)


def test_calibration_gives_the_earth_stand_its_carbon_use_efficiency():
    tables, temp = SimpleTables(), trace()
    for cue in (0.4, 0.5):
        r25 = plant.calibrate(tables, temp, cue)
        earth = plant.cycle(tables, 'earth', 0.0, temp, 86400.0, r25)
        assert earth['carbon_use_efficiency'] == pytest.approx(cue, rel=1e-9)


def test_an_efficiency_the_leaves_alone_rule_out_is_refused_with_their_share():
    # At 0.79 the stand keeps 1 - 1.25 * 0.79, about 1%, of its photosynthesis for upkeep: less than its leaves use.
    with pytest.raises(ValueError, match='of photosynthesis, exceeds the budget'):
        plant.calibrate(SimpleTables(), trace(), 0.79)


def test_idling_lowers_the_store_and_raises_growth():
    tables, temp = SimpleTables(), trace()
    r25 = plant.calibrate(tables, temp)
    runs = [plant.cycle(tables, 'moon', 0.0, temp, plant.MOON_DAY_S, r25, night=n) for n in (1.0, 0.5, 0.25, 0.1)]
    assert np.all(np.diff([r['store'] for r in runs]) < 0)
    assert np.all(np.diff([r['growth'] for r in runs]) > 0)
    # the store covers at least the respiration in the dark once photosynthesis stops
    assert runs[0]['store'] > 0.9 * runs[0]['respiration_in_the_dark'] - runs[0]['twilight_gpp']


def test_regrowing_the_canopy_pays_for_a_whole_canopy():
    tables, temp = SimpleTables(), trace()
    r25 = plant.calibrate(tables, temp)
    regrow = plant.cycle(tables, 'moon', 0.0, temp, plant.MOON_DAY_S, r25, regrow=True)
    keep = plant.cycle(tables, 'moon', 0.0, temp, plant.MOON_DAY_S, r25)
    cost = 5.0 * plant.LEAF_DRY_G_M2 * plant.CARBON_FRACTION * (1 + plant.GROWTH_RESPIRATION)
    assert regrow['new_canopy'] == pytest.approx(cost, rel=1e-6)
    assert regrow['growth'] < keep['growth'] and regrow['gpp'] < keep['gpp']


def test_along_cycle_wraps_around_midnight():
    centres = (np.arange(36) + 0.5) * 10.0 - 180.0
    values = np.cos(np.radians(centres))
    out = plant.along_cycle(centres, values, steps=720)
    hour = (np.arange(720) + 0.5) / 720 * 360.0 - 180.0
    assert out == pytest.approx(np.cos(np.radians(hour)), abs=0.02)


def test_leaf_respiration_follows_leaf_area():
    assert plant.leaf_respiration(2.5, 22.0) == pytest.approx(0.5 * plant.leaf_respiration(5.0, 22.0))


@pytest.fixture(scope='module')
def product():
    path = plant.RESULTS / 'plant.json'
    if not path.is_file():
        pytest.skip('results/plant.json has not been generated')
    return json.loads(path.read_text())


def test_stored_product(product):
    assert product['schema'] == plant.SCHEMA
    assert product['earth']['carbon_use_efficiency'] == pytest.approx(plant.CUE_EARTH, rel=1e-4)
    for side in plant.SIDES:
        for band in plant.BANDS:
            r = product['results'][side][band]
            ladder = [r[k] for k in ('earth_like', 'idle_half', 'idle_quarter', 'idle_tenth')]
            assert np.all(np.diff([x['store'] for x in ladder]) < 0)
            assert np.all(np.diff([x['growth'] for x in ladder]) > 0)
            assert r['regrow_leaves']['growth'] < r['earth_like']['growth']
    hours = product['twilight']['moon']['hours_after_sunset_above']
    assert hours['50.0'] < hours['10.0'] < hours['1.0']
    near = product['results']['near']
    assert near['60']['earth_like']['dark_hours'] < near['0']['earth_like']['dark_hours'] < 354.0
    no_twilight = product['sensitivity']['without_twilight']['earth_like']
    assert no_twilight['store'] > near['0']['earth_like']['store']
