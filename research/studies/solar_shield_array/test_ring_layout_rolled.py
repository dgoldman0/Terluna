"""The ring stack laid out under the rolled law's push: its binding, and what the smaller push changes."""
import json

import numpy as np
import pytest

from . import ring_layout as rl
from . import ring_layout_rolled as rr


@pytest.fixture(scope='module')
def product():
    if not rr.OUT.exists():
        pytest.skip('results/ring_layout_rolled.json has not been generated')
    return json.loads(rr.OUT.read_text())


@pytest.fixture(scope='module')
def built():
    return json.loads(rl.OUT.read_text())


def test_the_rolled_layout_binds_its_code_and_flies_the_laws_drive(product):
    assert product['schema'] == 'terluna.research.ring-layout/1' and not product['smoke']
    files = product['producer']['files']
    assert sorted(files) == sorted(rl.FILES+['research/studies/solar_shield_array/ring_layout_rolled.py'])
    for name, value in files.items():
        assert rl.digest(rl.ROOT/name) == value, name
    assert product['producer']['inputs']['results/frozen_rings.json'] == rl.digest(rl.FROZEN_OUT)
    assert product['producer']['inputs']['results/plane_motion.json'] == rl.digest(rl.PLANES_OUT)
    v = product['variant']
    scheme = json.loads(rl.ATTITUDE.read_text())['schemes'][rr.LAW]
    assert v['eccentricity_drive_over_as_built'] == pytest.approx(scheme['eccentricity_drive']['over_as_built'])
    assert v['flown_areal_mass_kg_m2'] == pytest.approx(rr.SIGMA/v['eccentricity_drive_over_as_built'])
    assert product['design']['sigma_kg_m2'] == rr.SIGMA                 # the tiles' own, the steering reach's
    assert v['plane_effects_not_flown']['crossing_height_rate_km_per_day'] == scheme['crossing_height_rate_km_per_day']


def test_the_smaller_push_shrinks_the_forced_eccentricity_and_raises_the_perilunes(product, built):
    """Every calibrated ring's forced eccentricity falls, by about the drive near the Moon's orbit plane and by less
    toward the stack's edges, where the tide's share of the balance grows; the flown rings' night halves ride
    farther out."""
    ratio = np.array(product['calibration']['forced_eccentricity'])/np.array(built['calibration']['forced_eccentricity'])
    drive = product['variant']['eccentricity_drive_over_as_built']
    assert np.all(ratio < 1.) and np.all(ratio > .9*drive)
    tilts = np.abs(np.array(product['calibration']['tilts_deg']))
    assert np.all(ratio[:, tilts == tilts.max()].mean() > ratio[:, tilts == tilts.min()].mean())
    perilune = lambda rings: np.array([r['radius_km']*(1-r['forced_eccentricity'])/(1+r['forced_eccentricity'])
                                       for r in rings])
    flown, ref = product['flown']['matched']['rings'], built['flown']['matched']['rings']
    assert perilune(flown).min() > perilune(ref).min()
    assert max(r['forced_eccentricity'] for r in flown) < max(r['forced_eccentricity'] for r in ref)
