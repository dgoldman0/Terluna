"""Finite near-natural local patterns, with explicit roll and coupled light.

The affine shooting model proposes initial velocities only. Individual tiles
are subsequently propagated with finite-area gravity and mutual interception.
Shared Sun-facing planes permit exact rectangle unions for each sampled source;
the neighbour reduction has a separately evaluated all-pairs exclusion guard.
"""
import numpy as np
from scipy.spatial.transform import Rotation, RotationSpline

from shared import constants as K
from .active_retime import UV, WEIGHTS
from .active_formation import solar_frame
from .cycling import ray_geometry, visible_sun, sail_basis
from .fleet import sun_points, projected_polygons
from .optical import unit, length, arriving_ray_vectors


def transported_roll(times, orientations):
    """Minimum normal-to-normal transport, then a smooth SO(3) command.

    Normal signs describe the same two-sided plane. The returned spline still
    requires rate, acceleration and beam-path checks on the propagated state.
    """
    frames = [orientations[0].copy()]
    for proposed in orientations[1:]:
        before = frames[-1]
        n = proposed[:, 2].copy()
        if n@before[:, 2] < 0:
            n *= -1
        axis = np.cross(before[:, 2], n)
        sine = length(axis)
        rotation = Rotation.from_rotvec(axis/max(sine, 1e-100)*
            np.arctan2(sine, n@before[:, 2])).as_matrix()
        frames.append(rotation@before)
    return RotationSpline(times, Rotation.from_matrix(frames))


def finite_gravity(env, t, q, a, b, side=10000.):
    a, b = np.broadcast_to(a, q.shape), np.broadcast_to(b, q.shape)
    points = q[:, None, :]+side/2*(UV[None, :, :1]*a[:, None, :]+
                                 UV[None, :, 1:]*b[:, None, :])
    return np.einsum('j,ijk->ik', WEIGHTS, env.gravity(t, points))


def rectangle_union(lo, hi, active):
    """Union area of batched, axis-aligned rectangles by exact x strips.

    Inputs have shape (batch, rectangles, 2). Degenerate or inactive rectangles
    contribute zero. No assumptions about pair/triple overlap multiplicity.
    """
    active = active & np.all(hi > lo, axis=-1)
    lo = np.where(active[..., None], lo, 0.)
    hi = np.where(active[..., None], hi, 0.)
    edges = np.sort(np.concatenate([lo[..., 0], hi[..., 0]], axis=1), axis=1)
    mid = (edges[:, 1:]+edges[:, :-1])/2
    present = active[:, None, :] & (mid[..., None] > lo[:, None, :, 0]) & (mid[..., None] < hi[:, None, :, 0])
    order = np.argsort(lo[..., 1], axis=1)
    low = np.take_along_axis(lo[..., 1], order, axis=1)[:, None, :]
    high = np.take_along_axis(hi[..., 1], order, axis=1)[:, None, :]
    present = np.take_along_axis(present, order[:, None, :], axis=2)
    top = np.maximum.accumulate(np.where(present, high, -np.inf), axis=2)
    previous = np.concatenate([np.full(top.shape[:2]+(1,), -np.inf), top[..., :-1]], axis=2)
    vertical = np.where(present, np.maximum(high-np.maximum(low, previous), 0.), 0.).sum(axis=2)
    return np.sum(np.diff(edges, axis=1)*vertical, axis=1)


def retarded_sun(sample):
    return arriving_ray_vectors(sample['positions']['sun'], sample['positions']['earth'],
        sample['sun_v'], sample['earth_v'], sample['sun_a'], sample['earth_a'])[0]


def mutual_visibility(state, frame, sample, neighbour, valid, suns=8, clear_side=9890.):
    """Visible area for each finite-Sun quadrature direction and each tile."""
    xyz = state[:, :3]@frame
    sources = sun_points(retarded_sun(sample), suns, .317)
    qi, qj = xyz[:, None, :], xyz[neighbour]
    upstream = valid & (qj[..., 0] > qi[..., 0])
    fractions = []
    h = clear_side/2
    for s in sources@frame:
        mag = (s[0]-qi[..., 0])/(s[0]-qj[..., 0])
        centre = s[1:]+(qj[..., 1:]-s[1:])*mag[..., None]-qi[..., 1:]
        lo = np.maximum(-h, centre-h*mag[..., None])
        hi = np.minimum(h, centre+h*mag[..., None])
        fractions.append(np.clip(1-rectangle_union(lo, hi, upstream)/clear_side**2, 0, 1))
    return np.array(fractions), sources


def neighbour_exclusion(state, frame, sample, neighbour, valid, clear_side=9890., allowance=50.):
    """All non-stencil pairs must be separated in projected x OR y.

    Solar-cone widening includes the off-axis source direction, finite distance,
    and a numerical/moving-interceptor allowance. Positive values exclude all
    other tiles at this instant. Time intervals need a relative-motion bound.
    """
    xyz = state[:, :3]@frame
    s = retarded_sun(sample)@frame
    cone = (K.SUN_RADIUS+np.max(np.abs(s[1:]))+np.max(np.abs(xyz[:, 1:])))/(
        s[0]-np.max(xyz[:, 0]))
    delta = np.abs(xyz[:, None, :]-xyz[None, :, :])
    margin = np.max(delta[..., 1:], axis=-1)-clear_side-cone*(delta[..., 0]+clear_side)-allowance
    margin[np.diag_indices(len(state))] = np.inf
    i = np.broadcast_to(np.arange(len(state))[:, None], neighbour.shape)
    margin[i[valid], neighbour[valid]] = np.inf
    return float(margin.min())


class NaturalPattern:
    def __init__(self, env, neighbour, valid, suns=8, side=10000.):
        self.env, self.neighbour, self.valid, self.suns, self.side = env, neighbour, valid, suns, side

    def forces(self, t, state, diagnostics=False):
        sample = self.env.at(t)
        frame = solar_frame(sample)
        normal = -frame[:, 0]
        light, sources = mutual_visibility(state, frame, sample, self.neighbour,
                                           self.valid, self.suns, self.side-110.)
        rays = ray_geometry(state[:, :3], sample)
        natural_light = visible_sun(rays['sun'], rays['earth'], rays['moon'])
        # The local pilot accepts only completely uneclipsed passages. Combined
        # body/array partial eclipse must use a joint source/area integral.
        if np.min(natural_light) < 1-1e-12:
            raise ValueError('Partial body eclipse invalidates this local illumination model')
        incident = sources[:, None, :]-state[None, :, :3]
        cosine = np.sum(-unit(incident)*normal, axis=-1)
        strength = np.mean(light*cosine**2*(K.AU/length(incident))**2, axis=0)
        scale = 2*self.env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*self.env.sigma)
        sail = (scale*((self.side-110.)/self.side)**2*strength)[:, None]*normal
        acceleration = finite_gravity(self.env, t, state[:, :3], frame[:, 1], frame[:, 2], self.side)+sail
        if diagnostics:
            return acceleration, light.mean(axis=0), sail
        return acceleration

    def derivative(self, t, flat):
        state = flat.reshape(-1, 6)
        return np.c_[state[:, 3:], self.forces(t, state)].ravel()


def local_coverage(state, sample, centre, radius=10000., limb_count=8, rotation=0., interior=0):
    """Polygon-union coverage on the entire local disk, sampled Sun sources.

    Shared Sun-facing physical squares have 9890 m clear width. Each source
    checks a continuum of receiver points, with a polygonal circle boundary.
    """
    import shapely
    from shapely.geometry import Point
    sun = retarded_sun(sample)
    u = unit(sun); _, b, c = sail_basis(sun)
    q = state[:, :3]-(state[:, 3:]+sample['moon_v'])*np.maximum(state[:, :3]@u, 0)[:, None]/K.SPEED_OF_LIGHT
    qc = centre[:3]-(centre[3:]+sample['moon_v'])*max(centre[:3]@u, 0)/K.SPEED_OF_LIGHT
    intercept = sun+(qc-sun)*(sun@u)/(sun@u-qc@u)
    disk = Point(intercept@b, intercept@c).buffer(radius, quad_segs=128)
    signs = np.array([[-1,-1],[1,-1],[1,1],[-1,1]])*9890/2
    vertices = q[:, None, :]+b*signs[None, :, :1]+c*signs[None, :, 1:]
    angle = rotation+np.arange(limb_count)*2*np.pi/limb_count
    sources = np.r_[sun[None, :], sun+K.SUN_RADIUS*(np.cos(angle)[:, None]*b+np.sin(angle)[:, None]*c)]
    if interior:
        sources = np.r_[sources, sun_points(sun, interior, rotation)]
    full = disk; worst = 0.; witness = None
    fractions = []
    for i, source in enumerate(sources):
        union = shapely.union_all(projected_polygons(vertices, source, u, b, c))
        missing = disk.difference(union)
        fractions.append(1-missing.area/disk.area)
        full = full.intersection(union)
        if missing.area > worst:
            worst = missing.area; point = missing.representative_point()
            witness = dict(source_index=i, receiver_xy_m=[point.x, point.y], missed_area_m2=worst)
    return dict(all_sampled_sun_coverage=float(full.area/disk.area),
        minimum_source_coverage=float(min(fractions)), worst_missed_ray=witness,
        window_centre_distance_m=float(np.hypot(intercept@b, intercept@c)), sources=len(sources))
