import json
import os

import numpy as np
import pytest

from research.studies.south_sea_ocean import inflows as I


def test_waters_join_across_the_seam_and_stay_apart_elsewhere():
    below = np.zeros((6, 12), bool)
    below[1, 0] = below[2, 11] = True          # touch diagonally across the seam
    below[4, 3:5] = True                       # one water in the middle
    below[4, 7] = True                         # another, two cells away
    body = I.label_waters(below)
    assert body[1, 0] == body[2, 11] != 0
    assert body[4, 3] == body[4, 4] != 0
    assert len({body[1, 0], body[4, 3], body[4, 7]}) == 3
    assert body[0, 0] == 0


def test_centroid_and_separation_on_the_sphere():
    lat = np.array([10.0, 10.0])
    lon = np.array([179.0, -179.0])
    c = I.centroid(lat, lon, np.ones(2))
    assert c[0] == pytest.approx(10.0, abs=0.01) and abs(c[1]) == pytest.approx(180.0, abs=1e-6)
    assert I.separation_deg((0.0, 0.0), (0.0, 90.0)) == pytest.approx(90.0)
    assert I.separation_deg((90.0, 0.0), (-90.0, 50.0)) == pytest.approx(180.0)


def test_the_six_seas_are_named_by_the_nearest_atlas_centroid():
    lat = 90.0 - (np.arange(18) + 0.5) * 10.0
    lon = (np.arange(36) + 0.5) * 10.0 - 180.0
    body = np.zeros((18, 36), int)
    centres = {433: (26.0, -12.0), 1417: (-56.0, -172.0), 874: (5.0, 88.0), 213: (58.0, 85.0), 1296: (-19.0, -94.0),
               714: (27.0, 148.0)}
    sizes = {433: 4, 1417: 3, 874: 1, 213: 1, 1296: 1, 714: 1}
    for k, (label, (la, lo)) in enumerate(centres.items(), start=1):
        r, c = int(np.abs(lat - la).argmin()), int(np.abs(lon - lo).argmin())
        body[r, c:c + sizes[label]] = k
    area = np.cos(np.radians(lat)) * 1e10
    atlas_bodies = [dict(label=label, centroid=list(c), area_km2=1.0) for label, c in centres.items()]
    named, _ = I.name_seas(body, area, lat, lon, atlas_bodies)
    assert set(named) == set(I.SEAS.values())
    assert named['south_pole_aitken']['label'] == 2 and named['humboldtianum']['label'] == 4


def test_the_product_is_current_and_reproduces_geographys_drainage():
    d = json.loads(I.OUT.read_text())
    assert d['schema'] == I.SCHEMA
    for path, digest in {**d['producer']['files'], **d['inputs']}.items():
        assert I.sha256(I.ROOT / path) == digest, f'{path} changed since inflows.json was written'
    committed = json.loads((I.ROOT / I.DRAINAGE).read_text())
    assert round(d['rules']['every_water']['to_the_outlets_m3s']) == committed['rivers']['discharge_to_sea_m3s']
    for name, sea in d['seas'].items():
        assert sea['area_km2'] == pytest.approx(sea['atlas_area_km2'], rel=0.07), name   # 16 against 4 pixels/degree


def test_each_rule_lets_more_water_through_to_the_two_great_seas():
    rules = json.loads(I.OUT.read_text())['rules']
    for name in I.TWO:
        flows = [rules[rule]['seas'][name]['inflow_m3s'] for rule in ('every_water', 'six_seas', 'two_seas')]
        assert flows == sorted(flows)
    assert rules['every_water']['to_other_waters_m3s'] > 0.8 * rules['every_water']['to_the_outlets_m3s']
    assert rules['two_seas']['to_other_waters_m3s'] == pytest.approx(0.0, abs=1.0)
    for rule in rules.values():
        for sea in rule['seas'].values():
            assert sea['inflow_m3s'] >= sum(m['discharge_m3s'] for m in sea['largest']) - 1.0
            assert sea['inflow_km3_yr'] == pytest.approx(sea['inflow_m3s'] * I.KM3_YR, abs=0.06)


@pytest.mark.skipif(os.environ.get('TERLUNA_SOUTH_SEA_REGENERATE') != '1',
                    reason='three routings, 3.5 minutes and 5 GB; set TERLUNA_SOUTH_SEA_REGENERATE=1')
def test_the_inflows_regenerate():
    assert json.loads(I.OUT.read_text()) == json.loads(json.dumps(I.build()))
