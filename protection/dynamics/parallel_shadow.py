"""Two-sided parallel finite squares with dynamic, finite-source blockers.

Used for uneclipsed validation prefixes with arbitrary shared attitude. Source
directions may illuminate either face. Source-dependent candidate lists replace
the fixed eight-neighbour stencil when a pattern rotates toward feathering.
Partial body eclipses explicitly reject this local model.
"""
import numpy as np
from scipy.spatial import cKDTree

from shared import constants as K
from .cycling import sail_basis, ray_geometry, visible_sun
from .fleet import sun_points
from .natural_pattern import rectangle_union, finite_gravity, retarded_sun
from .optical import unit, length


def source_neighbours(q, source, side=10000., every=False):
    """Conservative perspective bounding disks; no lattice topology assumed."""
    centre = q.mean(axis=0)
    relative = q-centre; source = source-centre
    u = unit(source); _, b, c = sail_basis(source)
    depth = relative@u
    distance = length(source)
    magnification = distance/(distance-depth)
    xy = np.c_[relative@b, relative@c]*magnification[:, None]
    # A finite square lies in a 3-D circumsphere. Perspective magnification
    # and the centre's off-axis angle bound its projected displacement.
    radius = side/np.sqrt(2)
    maximum = distance/(distance-depth.max()-radius)
    transverse = float(length(relative-depth[:, None]*u).max())
    bound = 2*radius*maximum*(1+(transverse+radius)/(distance-depth.max()-radius))+1.
    if every:
        neighbour = np.tile(np.arange(len(q)), (len(q), 1))
        return neighbour, ~np.eye(len(q), dtype=bool), bound
    pairs = cKDTree(xy).query_pairs(bound, output_type='ndarray')
    counts = np.bincount(pairs.ravel(), minlength=len(q)) if len(pairs) else np.zeros(len(q), int)
    width = max(1, int(counts.max()))
    neighbour = np.zeros((len(q), width), int); valid = np.zeros_like(neighbour, bool)
    used = np.zeros(len(q), int)
    for i, j in pairs:
        neighbour[i, used[i]] = j; valid[i, used[i]] = True; used[i] += 1
        neighbour[j, used[j]] = i; valid[j, used[j]] = True; used[j] += 1
    return neighbour, valid, bound


def source_illumination(q, frame, source, clear_side=9890., every=False):
    """Exact rectangle unions on either side of a shared square plane."""
    neighbour, valid, _ = source_neighbours(q, source, every=every)
    centre = q.mean(axis=0)
    # Coordinates are normal, edge a, edge b.
    basis = frame[:, [2, 0, 1]]
    xyz = (q-centre)@basis; s = (source-centre)@basis
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
    return np.clip(1-rectangle_union(lo, hi, active)/clear_side**2, 0., 1.)


class ParallelPattern:
    def __init__(self, env, command, suns=8, side=10000.):
        self.env, self.command, self.suns, self.side = env, command, suns, side

    def acceleration(self, t, state, diagnostics=False):
        sample = self.env.at(t); frame = self.command(t).as_matrix()
        q = state[:, :3]
        rays = ray_geometry(q, sample)
        if np.min(visible_sun(rays['sun'], rays['earth'], rays['moon'])) < 1-1e-12:
            raise ValueError('Partial body eclipse outside the uneclipsed validation-prefix model')
        sources = sun_points(retarded_sun(sample), self.suns, .317)
        coefficient = np.zeros(len(q)); light = []
        n = frame[:, 2]
        for source in sources:
            visible = source_illumination(q, frame, source, self.side-110.)
            incident = source-q
            cosine = -unit(incident)@n
            # The sign changes with the illuminated face; opposite solar
            # directions at exact feathering need not have identical blockers.
            coefficient += visible*cosine*abs(cosine)*(K.AU/length(incident))**2/self.suns
            light.append(visible)
        scale = 2*self.env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*self.env.sigma)
        sail = (scale*((self.side-110.)/self.side)**2*coefficient)[:, None]*n
        gravity = finite_gravity(self.env, t, q, frame[:, 0], frame[:, 1], self.side)
        if diagnostics:
            return gravity+sail, np.array(light), sail
        return gravity+sail


def minimum_clearance(state, frame, side=10000., margin=100.):
    q = state[:, :3]
    pairs = cKDTree(q).query_pairs(np.sqrt(2)*side+margin, output_type='ndarray')
    if not len(pairs):
        return np.inf, None, None
    coordinates = (q[pairs[:, 0]]-q[pairs[:, 1]])@frame
    separation = np.max(abs(coordinates)-np.array([side, side, 0.]), axis=1)
    index = int(np.argmin(separation))
    return float(separation[index]), pairs[index], coordinates[index]
