"""Checks of the drylands' reading: the shares weigh land cells by their area in UNEP's classes, and the product holds
together once computed: the classes add up, the two windows on the final surface share their land, and the rivers'
green land is a part of it."""
import json

import numpy as np
import pytest

from biosphere.ecology import moisture_map as bm
from research.studies.atlas_round import drylands as DL


def test_the_shares_weigh_land_by_area_in_unep_classes():
    lat = np.array([60.0, 0.0, -60.0])
    lon = np.array([0.0, 120.0, 240.0])
    land = np.array([[True, True, False], [True, True, True], [False, True, True]])
    ai = np.array([[0.01, 0.1, 9.0], [0.3, 0.6, 2.0], [9.0, 0.3, 2.0]])
    c = dict(lat=lat, lon=lon, land=land, ai=ai, rain=np.ones((3, 3)), pet=np.full((3, 3), 2.0))
    summary, cls, area = DL.shares(c)
    total = area[land].sum()
    expect = {name: area[land & (cls == k)].sum() / total for k, name in enumerate(DL.CLASSES)}
    for name in DL.CLASSES:
        assert summary['classes'][name] == pytest.approx(expect[name], abs=5e-5)
    assert sum(summary['classes'].values()) == pytest.approx(1.0, abs=5e-4)
    assert summary['arid_and_semi_arid'] == pytest.approx(expect['arid'] + expect['semi-arid'], abs=1e-4)
    assert summary['land_rain_mm_day'] == 1.0 and summary['land_demand_mm_day'] == 2.0
    assert list(bm.unep_class(np.array([0.01, 0.1, 0.3, 0.6, 2.0]))) == [0, 1, 2, 3, 4]


@pytest.fixture(scope='module')
def product():
    if not DL.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.drylands first')
    return json.loads(DL.OUT.read_text())


def test_the_product_holds_together(product):
    windows = product['windows']
    assert set(windows) == {label for label, _, _, _ in DL.WINDOWS}
    for w in windows.values():
        assert sum(w['classes'].values()) == pytest.approx(1.0, abs=5e-4)
        assert w['drylands'] == pytest.approx(sum(w['classes'][n] for n in DL.CLASSES[:4]), abs=5e-4)
    final = [w for label, w in windows.items() if label.startswith('final_surface')]
    assert len({w['land_km2'] for w in final}) == 1                 # the same surface
    for w in final:
        for g in w['rivers_green_land'].values():
            assert 0 <= g['share_of_the_land'] < 0.01
            assert g['arid_and_semi_arid_left'] <= w['arid_and_semi_arid']
