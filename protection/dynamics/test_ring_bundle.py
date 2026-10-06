import numpy as np

from shared import constants as K
from . import ring_bundle as rb
from .fast_parallel import load as parallel_load


class Env:
    """A minimal environment: Moon point gravity and a fixed Sun far along +x."""
    sigma, central_fraction = .05, .137950

    def at(self, t):
        z = np.zeros(3)
        sun = np.array([K.AU, 0., 0.])
        return dict(positions=dict(sun=sun, earth=np.array([-3.844e8, 0., 0.])),
                    sun_v=z, earth_v=z, sun_a=z, earth_a=z, moon_a=z, moon_v=z)

    def gravity(self, t, q):
        return -K.MOON_GM*q/np.linalg.norm(q, axis=-1)[..., None]**3


def circular(radius, phase, normal=np.array([0., 0., 1.])):
    r = radius*np.array([np.cos(phase), np.sin(phase), 0.])
    v = np.sqrt(K.MOON_GM/radius)*np.cross(normal, r/radius)
    return np.r_[r, v]


def test_a_time_track_starts_on_its_centre_and_follows_its_trajectory_both_ways():
    centre = circular(15e6, .3)
    lags = np.array([-223., -15., 0., 15., 223.])
    states = rb.time_track(lambda t, s: Env().gravity(t, s[:, :3]), centre, lags)
    assert np.array_equal(states[2], centre)
    expected = rb.kepler_shift(np.repeat(centre[None], len(lags), axis=0), lags)
    assert np.allclose(states[:, :3], expected[:, :3], rtol=0., atol=1e-2)
    one_sided = rb.time_track(lambda t, s: Env().gravity(t, s[:, :3]), centre, lags[lags >= 0])
    assert np.array_equal(one_sided[0], centre)


def test_built_rings_put_their_middle_tile_on_the_designed_centre():
    centres = []

    def track(centre, lags):
        centres.append(centre)
        return rb.time_track(lambda t, s: Env().gravity(t, s[:, :3]), centre, lags)
    states, ring, slot = rb.build(15e6, track, 3, 5, 8.5e3, 600., 9.3e3, np.array([1., 0., 0.]),
                                  np.array([0., 0., 1.]), -.2)
    for k, centre in enumerate(centres):
        assert np.array_equal(states[(ring == k) & (slot == 2)][0], centre)


def test_tile_frames_are_right_handed_radial_facing_and_tilted_about_the_normal():
    s = np.array([circular(15e6, .3), circular(15.2e6, -1.)])
    tilt = np.radians(1.5)
    f = rb.tile_frames(s, tilt)
    for k in range(2):
        assert np.allclose(f[k].T@f[k], np.eye(3), atol=1e-12)
        assert np.isclose(np.linalg.det(f[k]), 1.)
        radial = s[k, :3]/np.linalg.norm(s[k, :3])
        assert np.isclose(np.degrees(np.arccos(f[k][:, 2]@radial)), 1.5)
        assert np.allclose(f[k][:, 1], [0., 0., 1.], atol=1e-12)


def test_a_lone_tile_gets_the_parallel_backend_force_in_its_own_frame():
    s = circular(15e6, .2)[None]
    tilt = np.radians(1.5)
    forces = rb.RingBundleForces(Env(), tilt, [1.25])
    mine, _, _ = forces.light(0., s)
    theirs = parallel_load(Env(), 0., s, rb.tile_frames(s, tilt)[0])['reflected_force_N']
    assert np.allclose(mine, theirs, rtol=1e-12, atol=0.)


def test_rings_are_time_shifted_strings_with_the_designed_steps():
    def track(centre, lags):
        n = np.linalg.norm(centre[3:])/np.linalg.norm(centre[:3])
        phase = np.arctan2(centre[1], centre[0])
        radius = np.linalg.norm(centre[:3])
        return [circular(radius, phase+n*lag) if abs(centre[2]) < 1 else
                np.r_[centre[:3], centre[3:]] for lag in lags]
    states, ring, slot = rb.build(15e6, track, 3, 5, 8.5e3, 600., 0., np.array([1., 0., 0.]),
                                  np.array([0., 0., 1.]), -.2)
    radii = np.linalg.norm(states[:, :3], axis=1)
    for k in range(3):
        assert np.allclose(radii[ring == k], 15e6+(k-1)*600.)
        q = states[ring == k][np.argsort(slot[ring == k]), :3]
        assert np.allclose(np.linalg.norm(np.diff(q, axis=0), axis=1), 8.5e3, rtol=1e-6)


def test_own_attitude_clearance_reads_shingles_and_side_by_side_gaps():
    tilt = np.radians(1.5)
    a = circular(15e6, 0.)
    b = circular(15e6, 8.5e3/15e6)
    d, _ = rb.own_attitude_clearance(np.array([a, b]), tilt)
    assert abs(d-8.5e3*np.sin(tilt)) < 5.
    far = circular(15e6, 10.2e3/15e6)
    d, _ = rb.own_attitude_clearance(np.array([a, far]), 0.)
    assert 150. < d < 250.


def test_kepler_shift_follows_circular_and_eccentric_two_body_motion():
    from scipy.integrate import solve_ivp
    radius = 15e6
    n = np.sqrt(K.MOON_GM/radius**3)
    shifted = rb.kepler_shift(circular(radius, .4)[None], np.array([.01/n]))[0]
    assert np.allclose(shifted, circular(radius, .41), rtol=0, atol=1e-6*radius)
    eccentric = circular(radius, .4)
    eccentric[3:] *= 1.02
    def rhs(t, y):
        return np.r_[y[3:], -K.MOON_GM*y[:3]/np.linalg.norm(y[:3])**3]
    for dt in (450., -450., 3*3600.):
        reference = solve_ivp(rhs, (0., dt), eccentric, method='DOP853', rtol=1e-12, atol=1e-6).y[:, -1]
        mine = rb.kepler_shift(eccentric[None], np.array([dt]))[0]
        assert np.linalg.norm(mine[:3]-reference[:3]) < 1e-3
        assert np.linalg.norm(mine[3:]-reference[3:]) < 1e-6
    assert np.allclose(rb.kepler_shift(eccentric[None], np.array([0.]))[0], eccentric, rtol=0, atol=1e-9)


def test_copies_place_ring_neighbours_without_duplicating_representatives():
    reps = np.array([circular(15e6, 0.), circular(15.0006e6, 0.002)])
    lags = 8.5e3/np.linalg.norm(reps[:, 3:], axis=1)
    model = rb.ContinuousRings(Env(), np.radians(1.5), [1.25, 1.25], lags, reach=2)
    states, ring, index = model.copies(reps)
    assert np.all(index[:2] == 0) and np.count_nonzero(index == 0) == 2
    pitch = np.linalg.norm(states[(ring == 0) & (index == 1)][0, :3]-reps[0, :3])
    assert abs(pitch-8.5e3) < 1.
    # Ring 1's representative sits 0.002 rad ahead, about 3.5 lags; copies of ring 1 reach back to it.
    behind = index[ring == 1].min()
    assert behind <= -4
    assert len({(r, j) for r, j in zip(ring, index)}) == len(ring)
