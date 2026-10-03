"""Analytic cross-shore geometry independent of the lunar elevation archive."""
import json

import numpy as np
import pytest

from geography.shore_profiles import extract, numerical_profile
from shared.constants import MOON_RADIUS


def test_oblique_plane_recovers_its_normal_depth_and_shore(tmp_path):
    longitude=93.+np.arange(-40,41)/256
    latitude=2.+np.arange(-40,41)/256
    metres=MOON_RADIUS*np.pi/180*np.array([np.cos(np.deg2rad(2.)),1.])
    east,north=np.meshgrid((longitude-93)*metres[0],(latitude-2)*metres[1])
    depth=30.-(.6*east+.8*north)/20
    terrain=tmp_path/'terrain.npz'
    np.savez(terrain,longitude_deg=longitude,latitude_deg=latitude,height_m=-1000-depth,
             metadata=json.dumps(dict(schema='terluna.geography.coastal-grid/1',sea_level_m=-1000)))
    names=('west_face','east_face','southern_shore')
    stations=tmp_path/'stations.json'
    stations.write_text(json.dumps(dict(schema='terluna.climate.shore-exposure/1',cases=[dict(
        name='h852_s1_d36',sites={name:dict(longitude_deg=93.,latitude_deg=2.) for name in names})])))
    product=extract(terrain,stations)
    for profile in product['profiles'].values():
        np.testing.assert_allclose(profile['normal_east_north'],[.6,.8],atol=1e-12)
        np.testing.assert_allclose(profile['depth_m'],30-np.array(profile['distance_m'])/20,atol=1e-10)
        assert profile['shore_distance_m']==pytest.approx(600)
        assert profile['mean_submerged_slope']==pytest.approx(.05)
    for kind in ('rock','slope'):
        x,bed,settings=numerical_profile(product,'west_face',kind=kind,dx=.5,slope_denominator=20)
        np.testing.assert_allclose(bed[x<=100],30)
        # Both constructions share this exact analytic shoreline.
        assert settings['still_water_shore_m']==pytest.approx(700)
        assert np.interp(700,x,bed)==pytest.approx(0,abs=1e-10)
        assert bed[-1]<=-15+1e-10
