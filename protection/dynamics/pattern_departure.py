"""Temporary layer expansion and bounded shared-attitude departure proposals.

The reference and gravity response maps are proposal generators. Coupled
finite-source radiation and shadow forces must verify every accepted prefix.
No local plane or centre trajectory is held by translational feedback.
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.spatial.transform import Rotation, RotationSpline

from shared import constants as K
from .active_formation import solar_frame, lattice
from .cycling import ray_geometry, visible_sun
from .pattern_return import smooth_arc
from .parallel_shadow import minimum_clearance
from .optical import unit, length


def layer_labels(rows=19):
    offsets = lattice(rows, 8500., 4000.)[0]
    labels = offsets[:, 0]/4000.
    return labels-labels.mean()


def sun_frame(env, t):
    u, b, c = solar_frame(env.at(t)).T
    return np.column_stack([b, -c, -u])


def departure_command(env, start, end, delay, duration, axis, sign):
    """Quintic zero-end-rate/acceleration turn from face-on to feathered.

    axis is an initial solar-frame square edge (0 or 1), and sign is +/-1.
    The shared frame tracks the slow actual solar ephemeris throughout.
    """
    times = np.unique(np.r_[np.arange(start, end, 60.), end,
        start+delay, start+delay+duration])
    f = np.clip((times-start-delay)/duration, 0., 1.)
    angle = sign*np.pi/2*(10*f**3-15*f**4+6*f**5)
    vector = np.zeros((len(times), 3)); vector[:, axis] = angle
    frames = np.array([sun_frame(env, t) for t in times])@Rotation.from_rotvec(vector).as_matrix()
    return RotationSpline(times, Rotation.from_matrix(frames))


def expansion_maps(env, initial, start, end, burn):
    """Gravity variational response to one layer-dependent smooth impulse.

    The centre follows an unshadowed, Sun-facing sail reference. Its common
    force is only a screen approximation once the squares change attitude.
    """
    window = np.array([[start, start+burn]])
    def rhs(t, y):
        state = y[:6]; phi = y[6:42].reshape(6, 6); response = y[42:].reshape(6, 3)
        sample = env.at(t); u = solar_frame(sample)[:, 0]
        rays = ray_geometry(state[:3], sample)
        light = visible_sun(rays['sun'], rays['earth'], rays['moon'])
        cosine = unit(rays['sun'])@u
        scale = 2*env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*env.sigma)
        sail = -u*scale*light*cosine**2*(K.AU/length(rays['sun']))**2*(9890/10000)**2
        eps = 10.
        gravity = env.gravity(t, state[:3]+np.r_[np.eye(3), -np.eye(3)]*eps)
        gradient = ((gravity[:3]-gravity[3:])/(2*eps)).T
        jacobian = np.block([[np.zeros((3, 3)), np.eye(3)], [gradient, np.zeros((3, 3))]])
        drive = np.r_[np.zeros((3, 3)), np.eye(3)*smooth_arc(t, window)[0]]
        return np.r_[state[3:], env.gravity(t, state[:3])+sail,
            (jacobian@phi).ravel(), (jacobian@response+drive).ravel()]
    sol = solve_ivp(rhs, [start, end], np.r_[initial.mean(axis=0), np.eye(6).ravel(), np.zeros(18)],
        method='DOP853', rtol=2e-11, atol=1e-9, max_step=120., dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def expansion_state(maps, initial, impulses, times):
    y = maps.sol(np.atleast_1d(times)).T
    phi = y[:, 6:42].reshape(-1, 6, 6); response = y[:, 42:].reshape(-1, 6, 3)
    return y[:, None, :6]+np.einsum('tij,nj->tni', phi, initial-initial.mean(axis=0))+np.einsum('tij,nj->tni', response, impulses)


def frozen_paths(initial, frame, labels):
    """Static reachability diagnostic, without assigning maneuver cost."""
    rows = []
    for extra in [0., 2000., 4000., 6000., 8000.]:
        state = initial.copy()
        state[:, :3] -= labels[:, None]*extra*frame[:, 2]
        for axis in [0, 1]:
            for sign in [-1, 1]:
                values = []
                for angle in np.linspace(0., np.pi/2, 91):
                    vector = np.zeros(3); vector[axis] = sign*angle
                    values.append(minimum_clearance(state, frame@Rotation.from_rotvec(vector).as_matrix())[0])
                rows.append(dict(extra_layer_spacing_m=extra, axis=axis, sign=sign,
                    minimum_sampled_clearance_m=float(min(values)),
                    worst_angle_deg=int(np.argmin(values)), sampled_angles=91))
    return rows
