"""Rigid-tile gravity and reflected-band torques on the parallel-array model.

The spatial first moment of a shadow union supplies photon torque. The
illumination and per-aperture centre-ray approximation match the existing
force model; flexible deformation, absorption and thermal torque remain open.
"""
import numpy as np

from shared import constants as K
from .cycling import ray_geometry
from .fleet import sun_points
from .natural_pattern import retarded_sun
from .parallel_shadow import source_neighbours
from .collection import require_uneclipsed_tiles
from .optical import unit, length


def rectangle_moments(lo, hi, active):
    """Exact area and first x/y moments of batched axis-aligned unions."""
    active = active & np.all(hi > lo, axis=-1)
    lo, hi = np.where(active[..., None], lo, 0.), np.where(active[..., None], hi, 0.)
    edges = np.sort(np.concatenate([lo[..., 0], hi[..., 0]], axis=1), axis=1)
    mid = (edges[:, 1:]+edges[:, :-1])/2
    present = active[:, None, :] & (mid[..., None] > lo[:, None, :, 0]) & (mid[..., None] < hi[:, None, :, 0])
    order = np.argsort(lo[..., 1], axis=1)
    low = np.take_along_axis(lo[..., 1], order, axis=1)[:, None, :]
    high = np.take_along_axis(hi[..., 1], order, axis=1)[:, None, :]
    present = np.take_along_axis(present, order[:, None, :], axis=2)
    top = np.maximum.accumulate(np.where(present, high, -np.inf), axis=2)
    previous = np.concatenate([np.full(top.shape[:2]+(1,), -np.inf), top[..., :-1]], axis=2)
    bottom = np.maximum(low, previous)
    used = present & (high > bottom)
    height = np.where(used, high-bottom, 0.).sum(axis=2)
    moment_y = np.where(used, (high**2-bottom**2)/2, 0.).sum(axis=2)
    width = np.diff(edges, axis=1)
    return np.c_[np.sum(width*height, axis=1), np.sum(width*mid*height, axis=1),
                 np.sum(width*moment_y, axis=1)]


def illuminated_moments(q, frame, source, clear_side=9890.):
    """Visible clear area and its first moments in the two square edges."""
    neighbour, valid, _ = source_neighbours(q, source)
    centre = q.mean(axis=0); basis = frame[:, [2, 0, 1]]
    xyz, s = (q-centre)@basis, (source-centre)@basis
    qi, qj = xyz[:, None, :], xyz[neighbour]
    to_source = s[0]-qi[..., 0]
    with np.errstate(divide='ignore', invalid='ignore'):
        fraction = (qj[..., 0]-qi[..., 0])/to_source
        mag = to_source/(s[0]-qj[..., 0])
    active = valid & (fraction > 0) & (fraction < 1)
    mag = np.where(active, mag, 1.)
    shift = s[1:]+(qj[..., 1:]-s[1:])*mag[..., None]-qi[..., 1:]
    h = clear_side/2
    lo, hi = np.maximum(-h, shift-h*mag[..., None]), np.minimum(h, shift+h*mag[..., None])
    blocked = rectangle_moments(lo, hi, active)
    return np.c_[np.clip(clear_side**2-blocked[:, 0], 0., clear_side**2), -blocked[:, 1:]]


def gravity_torque(env, t, q, frame, side=10000., order=3):
    """Specific body-frame torque from uniform square mass, Nm/kg."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    x, y = np.meshgrid(nodes*side/2, nodes*side/2)
    offsets = x.ravel()[:, None]*frame[:, 0]+y.ravel()[:, None]*frame[:, 1]
    points = q[:, None, :]+offsets[None, :, :]
    residual = env.gravity(t, points)-env.gravity(t, q)[:, None, :]
    torque = np.einsum('k,nkj->nj', np.outer(weights, weights).ravel()/4,
                      np.cross(offsets[None, :, :], residual))
    return torque@frame


def parallel_load(env, t, state, frame, suns=16, rotation=.317, side=10000., gravity_order=3):
    """Actual-model gravity/radiation torque plus force and optical inventory."""
    sample = env.at(t); q = state[:, :3]
    require_uneclipsed_tiles(ray_geometry(q, sample), side)
    sources = sun_points(retarded_sun(sample), suns, rotation)
    mass = env.sigma*side**2; clear = side-110.; count = len(q)
    unshadowed = np.zeros(count); first = np.zeros(count)
    front = np.zeros(count); back = np.zeros(count)
    radiation = np.zeros((count, 3)); force = np.zeros((count, 3))
    for source in sources:
        area, mx, my = illuminated_moments(q, frame, source, clear).T
        ray = source-q; cosine = unit(ray)@frame[:, 2]
        irradiance = K.SOLAR_CONSTANT*(K.AU/length(ray))**2/suns
        power_density = irradiance*abs(cosine)
        unique = power_density*area
        first += unique; unshadowed += power_density*clear**2
        front += np.where(cosine >= 0, unique, 0.)
        back += np.where(cosine < 0, unique, 0.)
        pressure = -2*env.central_fraction/K.SPEED_OF_LIGHT*power_density*cosine
        force += (pressure*area)[:, None]*frame[:, 2]
        radiation += pressure[:, None]*np.c_[my, -mx, np.zeros(count)]/mass
    optical = dict(unshadowed_aperture_W=unshadowed,
        first_intercept_bolometric_equivalent_W=first,
        front_first_intercept_W=front, back_first_intercept_W=back,
        redirected_band_optical_W=first*env.central_fraction,
        mutual_shadow_equivalent_loss_W=unshadowed-first)
    return dict(gravity_torque_per_mass=gravity_torque(env, t, q, frame, side, gravity_order),
        radiation_torque_per_mass=radiation, reflected_force_N=force, optical=optical)
