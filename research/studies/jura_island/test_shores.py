import json
import math

import numpy as np
import pytest

from research.studies.jura_island import shores
from shared.constants import MOON_RADIUS
from shared.provenance import constants_changed, constants_used


def cone(height_m=600.0, radius_km=30.0):
    """A cone island at the frame's centre on the 16 px/deg grid, the sea round it."""
    lat = 41.0 + (np.arange(8 * shores.PPD) + 0.5) / shores.PPD
    lon = 318.0 + (np.arange(12 * shores.PPD) + 0.5) / shores.PPD
    km = 2 * math.pi * MOON_RADIUS / 1000 / 360
    y = (lat[:, None] - shores.CENTRE[0]) * km
    x = (lon[None, :] - shores.CENTRE[1]) * km * np.cos(np.radians(lat))[:, None]
    r = np.hypot(x, y)
    return lat, lon, height_m * (1 - r / radius_km)


def test_a_cone_island_measures_as_a_cone():
    lat, lon, h = cone()
    m = shores.measure(lat, lon, h)
    assert m['island']['area_km2'] == pytest.approx(math.pi * 30.0**2, rel=0.03)
    assert m['island']['summit_m'] == pytest.approx(600.0, rel=0.05)    # the apex lies between pixel centres
    below_250 = sum(v for k, v in m['heights_share'].items() if int(k.split('-')[1].split()[0]) <= 250)
    assert below_250 == pytest.approx(1 - (1 - 250 / 600) ** 2, abs=0.03)
    assert m['coast_km'] == pytest.approx(2 * math.pi * 30.0, rel=0.1)
    # a 2% shore: the low case's 15 m reach makes a belt 0.56 km wide (its supply limit) round 188 km of coast
    low, high = m['belts']['low']['area_km2'], m['belts']['high']['area_km2']
    assert 0 < low < high
    assert low == pytest.approx(2 * math.pi * 30.0 * min(0.015 / 0.02, 0.56), rel=0.5)


def test_saved_product_pins_current_dependencies_and_reproduces(tmp_path, monkeypatch):
    saved = json.loads(shores.OUT.read_text())
    assert saved['schema'] == shores.SCHEMA
    provenance = saved['producer']
    for path, recorded in provenance['files'].items():
        assert shores.digest(shores.ROOT / path) == recorded, path
    assert not constants_changed(provenance['constants'])
    assert provenance['constants'] == constants_used([shores.ROOT / p for p in provenance['files']])
    if not all((shores.ROOT / p).exists() for p in shores.INPUTS):
        pytest.skip('LOLA and GRAIL inputs are not fetched (geography/fetch_inputs.py)')
    for path, recorded in provenance['inputs'].items():
        assert shores.digest(shores.ROOT / path) == recorded, path
    monkeypatch.setattr(shores, 'OUT', tmp_path / 'shores.json')
    shores.run()
    assert json.loads(shores.OUT.read_text()) == saved
