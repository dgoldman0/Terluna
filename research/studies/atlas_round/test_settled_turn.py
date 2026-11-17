"""Checks of the atlas round's settled turn: it reproduces the geology lane's own settled turn from the lane's load, its
masses add up, and its frame is a rotation."""
import json

import numpy as np
import pytest

from research.studies.atlas_round import settled_turn as ST


@pytest.fixture(scope='module')
def product():
    if not ST.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.settled_turn first')
    return json.loads(ST.OUT.read_text())


def test_the_lanes_own_settled_turn_is_reproduced(product):
    listed, reproduced = product['series_settled']['listed'], product['series_settled']['reproduced']
    for key in ('pole_deg', 'pole_toward_lon_deg', 'sub_earth_lat_deg', 'sub_earth_lon_deg'):
        assert reproduced[key] == pytest.approx(listed[key], abs=0.01), key


def test_the_final_waters_masses_add_up_and_its_frame_is_a_rotation(product):
    m = product['final_water_mass_kg']
    assert m['all'] == pytest.approx(m['seas'] + m['lakes'], rel=1e-9)
    frame = np.array(product['frame_today_to_turned'])
    assert np.allclose(frame @ frame.T, np.eye(3), atol=1e-9) and np.linalg.det(frame) == pytest.approx(1.0)
    settled = product['cases']['final_water']['settled']
    assert np.degrees(np.arccos(frame[2, 2])) == pytest.approx(settled['pole_deg'], abs=0.01)
