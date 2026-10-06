import json

import numpy as np

from shared import constants as K
from shared.provenance import constants_changed
from protection.dynamics import ring_bundle as rb
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from . import photon_control as pc

PRODUCT = json.loads(pc.OUT.read_text())


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.photon-control/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(pc.FILES)
    for name, value in producer['files'].items():
        assert pc.digest(pc.ROOT/name) == value, name
    assert producer['inputs']['results/ring_keeping.json'] == pc.digest(pc.KEPT_OUT)
    if DEFAULT_KERNEL.exists():
        assert producer['inputs']['de440s.bsp'] == pc.digest(DEFAULT_KERNEL)
    assert not constants_changed(producer['constants'])


def test_a_tile_facing_the_sun_gets_the_specular_force_and_tilting_turns_it():
    normal = np.array([[0., 0., 1.]])
    sun = np.array([[0., 0., -1.]])
    force = pc.sail(normal, sun, np.array([K.SOLAR_CONSTANT]), np.array([1.]), .138)
    assert np.allclose(force[0], [0., 0., 2*.138*K.SOLAR_CONSTANT/K.SPEED_OF_LIGHT*rb.CLEAR**2])
    frames = np.eye(3)[None]
    mass = np.array([5e5])
    demand = np.array([[1., 0., 0.]])
    small = pc.reach(frames, sun, np.array([K.SOLAR_CONSTANT]), np.array([1.]), .138, mass, demand, 1., 0.)
    large = pc.reach(frames, sun, np.array([K.SOLAR_CONSTANT]), np.array([1.]), .138, mass, demand, 5., 0.)
    assert 0. < small[0] < large[0]
    # A tilt of delta turns the force by about delta at normal incidence.
    assert np.isclose(small[0], force[0, 2]/mass[0]*np.sin(np.radians(1.)), rtol=.01)


def test_more_tilt_and_trim_reach_more_of_the_keeping_demand():
    levels = list(PRODUCT['translation'])
    for part in ('interior', 'edges'):
        reach = [PRODUCT['translation'][level][part]['reach_over_demand'] for level in levels]
        assert all(b > a for a, b in zip(reach, reach[1:]))
    assert PRODUCT['translation']['2_deg_0.25']['interior']['reach_over_demand'] > 1.
    assert abs(sum(PRODUCT['first_tilt_level_that_reaches'].values())-1.) < 1e-9


def test_shingle_shadows_put_a_steady_pitch_torque_on_each_tile_beyond_reflectivity_trim():
    attitude = PRODUCT['attitude']
    pitch = attitude['by_axis']['pitch_about_orbit_normal']
    assert pitch['radiation_orbit_mean_N_m_median'] > 5*abs(attitude['by_axis']['roll_about_along_track']
                                                             ['radiation_orbit_mean_N_m_median'])
    assert attitude['trim_needed_lit_median'] > 1.
    assert attitude['pressure_centre_offset_m']['median'] > 500.
    # The 1.5 degree shingle tilt holds the tile off gravity gradient's equilibrium in pitch.
    assert pitch['gravity_mean_N_m'] < 0. and abs(pitch['gravity_mean_N_m']) > 1000.


def test_gravity_torque_scales_with_areal_mass():
    attitude = PRODUCT['attitude']
    by_mass = attitude['by_areal_mass']
    sigma = PRODUCT['tile_areal_mass_kg_m2']
    for name, row in by_mass.items():
        areal = float(name.split('_')[0])/1e3
        assert np.isclose(row['gravity_torque_N_m_median'], attitude['gravity_torque_N_m_median']*areal/sigma, rtol=.02)
        assert np.isclose(row['sail_reach_factor'], sigma/areal)
