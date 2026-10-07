import json

import numpy as np

from shared.provenance import constants_changed
from protection.dynamics import ring_bundle as rb
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from . import attitude_schemes as at

PRODUCT = json.loads(at.OUT.read_text())
GEOM = (15e6, 1000., 9286.456, np.radians(1.5))


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.attitude-schemes/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(at.FILES)
    for name, value in producer['files'].items():
        assert at.digest(at.ROOT/name) == value, name
    assert producer['inputs']['results/ring_keeping.json'] == at.digest(at.KEPT_OUT)
    assert producer['inputs']['results/frozen_rings.json'] == at.digest(at.FROZEN_OUT)
    if DEFAULT_KERNEL.exists():
        assert producer['inputs']['de440s.bsp'] == at.digest(DEFAULT_KERNEL)
    assert not constants_changed(producer['constants'])


def test_turned_frames_stay_orthonormal_and_roll_turns_about_the_line_of_centres():
    state = np.array([15e6, 0., 0., 0., 1800., 0.])
    for roll, pitch in ((0., 0.), (.7, 0.), (0., .3), (-1.2, .1)):
        f = at.body_frames(state, np.radians(1.5), roll, pitch)
        assert np.allclose(f.T@f, np.eye(3), atol=1e-12)
    # The velocity is the roll axis: it keeps its components in the turned frame.
    plain = at.body_frames(state, np.radians(1.5))
    rolled = at.body_frames(state, np.radians(1.5), .9)
    t = state[3:]/np.linalg.norm(state[3:])
    assert np.allclose(t@plain, t@rolled, atol=1e-12)


def test_parallel_squares_are_as_far_apart_as_their_planes_or_their_edges():
    axes = np.eye(3)
    assert np.isclose(at.square_gap(axes, np.array([[300., 200., 50.]]), 1e4)[0], 50.)
    assert np.isclose(at.square_gap(axes, np.array([[10300., 0., 0.]]), 1e4)[0], 300.)
    assert np.isclose(at.square_gap(axes, np.array([[10300., 10400., 0.]]), 1e4)[0], 500.)


def test_pitching_far_from_radial_meets_the_next_ring_but_rolling_to_the_free_side_does_not():
    # The shingle alone keeps 222 m between a ring's own tiles.
    assert np.isclose(at.worst_gap(0., 0., 0., *GEOM[:3], 8500., GEOM[3], 3, 10.), 8500*np.sin(GEOM[3]), rtol=1e-3)
    assert at.worst_gap(0., 0., np.radians(10.), *GEOM[:3], 8500., GEOM[3], 3, 10.) < 150.
    assert at.worst_gap(np.radians(45.), np.radians(-60.), 0., *GEOM[:3], 8500., GEOM[3], 3, 10.) > 150.
    assert at.worst_gap(np.radians(45.), np.radians(9.), 0., *GEOM[:3], 8500., GEOM[3], 3, 10.) < 150.
    assert at.worst_gap(np.radians(135.), np.radians(60.), 0., *GEOM[:3], 8500., GEOM[3], 3, 10.) > 150.


def test_rolled_laws_rest_radial_in_service_and_at_the_nodes_and_steering_ends_before_the_night():
    u = np.radians([0., 20., 90., -90.])
    assert np.allclose(at.rolled(u, 1., 25., 30.), 0.)
    assert at.rolled(np.radians(57.), 1., 25., 30.) < 0 < at.rolled(np.radians(150.), 1., 25., 30.)
    phases = np.radians(np.linspace(-180., 180., 3601))
    for law in (at.rolled(phases, 1., 25., 30.), at.steering(phases, 1., 25., 30.)):
        assert np.max(np.abs(np.diff(law))) < .01
    assert np.allclose(at.steering(np.radians([-150., -60., 179.9]), 1., 25., 30.), 0., atol=1e-4)


def test_edge_on_turns_collide_with_the_next_ring_and_rolled_laws_keep_the_shingle_gap():
    laws = PRODUCT['clearance']['laws']
    assert laws['edge_on_pitch']['min_gap_m'] < 150. and laws['edge_on_pitch']['share_of_orbit_below_clearance'] > .5
    for name in ('as_built', 'rolled_30_deg', 'rolled_60_deg', 'rolled_80_deg', 'steering'):
        assert laws[name]['min_gap_m'] > 200.
    for row in PRODUCT['clearance']['envelope']:
        spans = row['pitch_allowed_deg']
        assert min(s[0] for s in spans) >= -5. and max(s[1] for s in spans) <= 8.


def test_the_as_built_scheme_reproduces_photon_control():
    photon = json.loads((at.HERE/'results/photon_control.json').read_text())
    built = PRODUCT['schemes']['as_built']
    assert np.isclose(built['radiation_orbit_mean_N_m']['pitch'],
                      photon['attitude']['by_axis']['pitch_about_orbit_normal']['radiation_orbit_mean_N_m_median'],
                      rtol=.05)
    assert np.isclose(built['pressure_centre_offset_m']['median'],
                      photon['attitude']['pressure_centre_offset_m']['median'], rtol=.05)


def test_a_centred_filter_removes_most_of_the_steady_pitch_torque_without_changing_the_push():
    built = PRODUCT['schemes']['as_built']
    for name in ('balanced_200_m', 'balanced_400_m', 'balanced_560_m'):
        row = PRODUCT['schemes'][name]
        assert row['radiation_orbit_mean_N_m']['pitch'] < .5*built['radiation_orbit_mean_N_m']['pitch']
        assert row['store_N_m_s']['trim_0.25']['median'] < .2*built['store_N_m_s']['trim_0.25']['median']
        assert abs(row['eccentricity_drive']['over_as_built']-1.) < .01
    torque = [PRODUCT['schemes'][f'balanced_{d:g}_m']['radiation_orbit_mean_N_m']['pitch'] for d in (200., 400., 800.)]
    assert torque[0] < torque[1] < torque[2]


def test_rolling_halves_the_push_at_most_and_needs_large_attitude_torques():
    built = PRODUCT['schemes']['as_built']
    drives = [PRODUCT['schemes'][f'rolled_{c:g}_deg']['eccentricity_drive']['over_as_built'] for c in (30., 60., 80.)]
    assert 1. > drives[0] > drives[1] > drives[2] > .4
    for c in (30., 60., 80.):
        row = PRODUCT['schemes'][f'rolled_{c:g}_deg']
        assert row['attitude_torque_N_m']['p95'] > 10*built['attitude_torque_N_m']['p95']
    assert PRODUCT['edge_on_pitch']['eccentricity_drive_over_as_built'] < drives[2]
    # A roll on one side of the orbit only turns the plane: the crossing height moves by under a kilometre a day.
    rate = PRODUCT['schemes']['steering']['crossing_height_rate_km_per_day']['median']
    assert .1 < abs(rate) < 2.


def test_a_moving_mass_trades_ballast_for_store():
    rows = PRODUCT['moving_mass']['62.7_g_m2']['by_reach']
    reach = sorted(rows, key=lambda k: float(k[:-2]))
    store = [rows[k]['store_N_m_s']['trim_0.25'] for k in reach]
    ballast = [rows[k]['ballast_over_tile_mass'] for k in reach]
    assert all(b > a for a, b in zip(ballast, ballast[1:]))
    assert store[-1] < store[0]
    assert PRODUCT['moving_mass']['62.7_g_m2']['target_offset_m']['median'] > 500.
