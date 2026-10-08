"""Simultaneous finite membranes: beam guards, ray unions and swept exclusions.

All tile centres are separate initial-value problems at a common ephemeris
epoch. Finite squares are ideal rigid optical surfaces. Geometry is evaluated
at reception time with first-order moving-tile retardation; no orbit is reused
at another time. The force law remains the cycling model's central-ray ideal
spectral actuator. A population with mutual shadows needs coupled forces or
exclusion of the conflicting complete trajectories.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from shared import constants as K
from .cycling import ray_geometry, sail_basis, service_coast_angles, visible_sun
from .optical import length, unit, arriving_ray_vectors


def plane_normal(sun, angles):
    e, b, c = sail_basis(sun)
    a, h = angles[..., 0], angles[..., 1]
    return e*np.cos(a)[..., None]+np.sin(a)[..., None]*(
        b*np.cos(h)[..., None]+c*np.sin(h)[..., None])


def beam_margins(q, sun, normal, sample, protected_radius, side, pointing=.001):
    """Clear the full solar cone, tile corners, pointing error and target motion.

    Earth's exclusion radius includes 100 km above its surface. The motion
    allowance uses barycentric target speed through the entire flight time.
    Pointing is a normal error, hence twice its value in the reflected beam.
    """
    e = -unit(sun)
    outgoing = e-2*np.sum(e*normal, axis=-1)[..., None]*normal
    solar = np.arcsin(K.SUN_RADIUS/length(sun))
    result = []
    for vector, radius, velocity in (
        (-q, protected_radius, sample["moon_v"]),
        (sample["positions"]["earth"]-q, K.EARTH_RADIUS+100e3, sample["earth_v"]),
    ):
        d = length(vector)
        inflation = side/np.sqrt(2)+length(velocity)*d/K.SPEED_OF_LIGHT+100.
        angle = np.arctan2(length(np.cross(outgoing, vector)), np.sum(outgoing*vector, axis=-1))
        result.append(angle-np.arcsin(np.clip((radius+inflation)/d, 0, 1))-solar-2*pointing)
    return np.stack(result, axis=-1)


def safe_normals(state, sample, protected_radius, side, pointing=.001):
    """Service/coast normal, nudged away from Earth and the protected Moon."""
    q = state[..., :3]
    rays = ray_geometry(q, sample)
    # Expand service to include every corner of the finite membrane.
    angles = service_coast_angles(state, rays, sample, protected_radius+side/np.sqrt(2))
    n = plane_normal(rays["sun"], angles)
    e, b, _ = sail_basis(rays["sun"])
    changed = np.zeros(q.shape[:-1], dtype=bool)
    for _ in range(3):
        margins = beam_margins(q, rays["sun"], n, sample, protected_radius, side, pointing)
        for j, target in enumerate((-q, sample["positions"]["earth"]-q)):
            bad = margins[..., j] < .001
            if not np.any(bad):
                continue
            changed |= bad
            outgoing = e-2*np.sum(e*n, axis=-1)[..., None]*n
            target = unit(target)
            tangent = outgoing*np.sum(outgoing*target, axis=-1)[..., None]-target
            tangent = np.where((length(tangent) < 1e-10)[..., None], b, tangent)
            tangent = unit(tangent)
            delta = np.maximum(.001-margins[..., j], 0)
            new_out = outgoing*np.cos(delta)[..., None]+tangent*np.sin(delta)[..., None]
            difference = e-new_out
            candidate = difference/np.maximum(length(difference)[..., None], 1e-100)
            n = np.where(bad[..., None], candidate, n)
    return n, rays, changed


def fleet_acceleration(env, t, state, protected_radius, side, pointing=.001, diagnostics=False):
    sample = env.at(t)
    n, rays, changed = safe_normals(state, sample, protected_radius, side, pointing)
    cosine = np.maximum(np.sum(n*(-unit(rays["sun"])), axis=-1), 0)
    visibility = visible_sun(rays["sun"], rays["earth"], rays["moon"])
    scale = 2*env.central_fraction*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*env.sigma)
    sail = (scale*(K.AU/length(rays["sun"]))**2*visibility*cosine**2)[..., None]*n
    acceleration = env.gravity(t, state[..., :3])+sail
    if diagnostics:
        return dict(normal=n, acceleration=acceleration, sail=sail, visibility=visibility,
                    projected=cosine, beam_guard_changed=changed, **rays,
                    margins=beam_margins(state[..., :3], rays["sun"], n, sample,
                                         protected_radius, side, pointing))
    return acceleration


def initial_fleet(env, radii, planes, phases, arrangement="isotropic", phase_offset=0.):
    """Distinct osculating circular lunar seeds, all initialized at t=0.

    Plane normals cover the retrograde hemisphere in the instantaneous lunar
    frame. These are initial guesses, not presumed circular ephemeris orbits.
    """
    states, labels = [], []
    frame = env.frames(0.)
    golden = (np.sqrt(5)-1)/2
    for family, radius in enumerate(radii):
        for j in range(planes):
            z = -1. if arrangement == "planar" else -(j+.5)/planes
            phi = 2*np.pi*(j*golden+family*.173)
            h = np.array([np.sqrt(1-z*z)*np.cos(phi), np.sqrt(1-z*z)*np.sin(phi), z])
            axis = np.array([1., 0., 0.]) if abs(h[0]) < .9 else np.array([0., 1., 0.])
            b = unit(np.cross(h, axis)); c = np.cross(h, b)
            angle = 2*np.pi*(np.arange(phases)/phases+(j*golden+family*.37+phase_offset)/phases)
            q = radius*(np.cos(angle)[:, None]*b+np.sin(angle)[:, None]*c)
            v = np.sqrt(K.MOON_GM/radius)*(-np.sin(angle)[:, None]*b+np.cos(angle)[:, None]*c)
            states.append(np.column_stack([q@frame.T, v@frame.T]))
            labels += [[family, j, k] for k in range(phases)]
    return np.concatenate(states), np.asarray(labels)


def square_vertices(q, normal, sun, side):
    """Square with a reproducible roll angle; side is the clear optical width."""
    _, b, c = sail_basis(sun)
    a = b-normal*np.sum(b*normal, axis=-1)[..., None]
    a = np.where((length(a) < 1e-8)[..., None], c, a)
    a = unit(a); d = np.cross(normal, a)
    signs = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])*.5*side
    return q[:, None, :]+a[:, None, :]*signs[None, :, :1]+d[:, None, :]*signs[None, :, 1:]


def sun_points(sun, count, rotation=0.):
    """Equal-area disk quadrature, with an independent rotation for refinement."""
    _, b, c = sail_basis(sun)
    i = np.arange(count)
    r = K.SUN_RADIUS*np.sqrt((i+.5)/count)
    angle = i*np.pi*(3-np.sqrt(5))+rotation
    return sun+r[:, None]*(np.cos(angle)[:, None]*b+np.sin(angle)[:, None]*c)


def projected_polygons(vertices, source, u, b, c):
    """Perspective ray intersections with the target's central normal plane."""
    import shapely
    delta = vertices-source
    factor = -np.dot(source, u)/np.einsum("...i,i->...", delta, u)
    p = source+factor[..., None]*delta
    xy = np.stack([p@b, p@c], axis=-1)
    return shapely.polygons(xy)


def coverage_snapshot(state, normal, sample, side, protected_radius, sun_count=32,
                      seam=10., placement=50., rotation=0., circle_resolution=256):
    """Union of finite retarded square shadows for each finite-Sun direction.

    Seams remove total `seam` width, plus `placement` at each edge. Overlaps
    count once. The intersection across solar directions measures areas fully
    covered on that quadrature. It is not a certificate for unsampled rays.
    Earth eclipses are excluded from the incident-ray denominator.
    """
    import shapely
    from shapely.geometry import Point, GeometryCollection
    q, v = state[:, :3], state[:, 3:]
    sun, _ = arriving_ray_vectors(sample["positions"]["sun"], sample["positions"]["earth"],
        sample["sun_v"], sample["earth_v"], sample["sun_a"], sample["earth_a"])
    u = unit(sun); _, b, c = sail_basis(sun)
    tau = np.maximum(q@u, 0)/K.SPEED_OF_LIGHT
    # Inertial position at interception relative to the reception-time Moon.
    retarded = q-(v+sample["moon_v"])*tau[:, None]
    clear_side = side-seam-2*placement
    if clear_side <= 0:
        raise ValueError("Seam and placement allowance consume the membrane")
    radial = length(retarded-u*(retarded@u)[:, None])
    keep = ((retarded@u > 0) &
            (radial < protected_radius+side+np.maximum(retarded@u, 0)*K.SUN_RADIUS/length(sun)))
    vertices = square_vertices(retarded[keep], normal[keep], np.broadcast_to(sun, retarded[keep].shape), clear_side)
    disk = Point(0, 0).buffer(protected_radius, quad_segs=circle_resolution)
    core = Point(0, 0).buffer(K.MOON_RADIUS, quad_segs=circle_resolution)
    complete = disk; natural_complete = disk
    covered = incident = core_covered = core_incident = multiple = 0.
    sources = sun_points(sun, sun_count, rotation)
    # Body motion during light flight is included at first order. The projected
    # circular terrestrial silhouette approximation is recorded in the product.
    earth = sample["positions"]["earth"]
    earth = earth-sample["earth_v"]*max(np.dot(earth, u), 0)/K.SPEED_OF_LIGHT
    for source in sources:
        if len(vertices):
            polygons = projected_polygons(vertices, source, u, b, c)
            merged = shapely.intersection(shapely.union_all(polygons), disk)
            summed = float(np.sum(shapely.area(shapely.intersection(polygons, disk))))
        else:
            merged = GeometryCollection(); summed = 0.
        eclipse = GeometryCollection()
        ed = np.dot(earth, u)
        if 0 < ed < np.dot(source, u):
            magnification = np.dot(source, u)/(np.dot(source, u)-ed)
            centre = source+(earth-source)*magnification
            eclipse = Point(np.dot(centre, b), np.dot(centre, c)).buffer(
                K.EARTH_RADIUS*magnification, quad_segs=circle_resolution)
            eclipse = eclipse.intersection(disk)
        available = disk.difference(eclipse)
        covered += merged.intersection(available).area
        incident += available.area
        core_available = core.difference(eclipse)
        core_covered += merged.intersection(core_available).area
        core_incident += core_available.area
        multiple += max(summed-merged.area, 0.)
        complete = complete.intersection(merged.union(eclipse))
        natural_complete = natural_complete.intersection(eclipse)
    fully_available_area = disk.area-natural_complete.area
    return dict(ray_coverage=covered/incident if incident else None,
        core_ray_coverage=core_covered/core_incident if core_incident else None,
        all_sampled_sun_coverage=(complete.area-natural_complete.area)/fully_available_area
            if fully_available_area else None,
        natural_eclipse_fraction=1-incident/(sun_count*disk.area),
        excess_intersection_area_fraction=multiple/(sun_count*disk.area),
        active_finite_tiles=int(np.count_nonzero(keep)),
        circle_area_relative_error=1-disk.area/(np.pi*protected_radius**2))


def swept_conflicts(times, states, side, acceleration_bound=.15, clearance=100.,
                    solar_angular_radius=.0048, sun_directions=None):
    """Conservative pair exclusions over every integration-output interval.

    Linear closest approaches are widened by a*dt²/4 for two accelerated
    centres. Circumscribed squares make the collision guard orientation
    independent. A separate projected guard identifies possible mutual solar
    shadowing, including the whole solar disk. It deliberately over-rejects.
    """
    collision, shadow = {}, {}
    threshold = np.sqrt(2)*side+clearance
    radius = length(states[..., :3])
    max_radius = float(radius.max())
    for k, dt in enumerate(np.diff(times)):
        q0, q1 = states[k, :, :3], states[k+1, :, :3]
        mid = (q0+q1)/2; move = q1-q0
        sag = acceleration_bound*dt*dt/4
        travel = float(length(move).max())
        def collect(coords, delta, limit, record, axial=None):
            pairs = cKDTree(coords).query_pairs(limit+travel+sag, output_type="ndarray")
            if not len(pairs):
                return
            i, j = pairs.T
            start = coords[i]-coords[j]-(delta[i]-delta[j])/2
            change = delta[i]-delta[j]
            f = np.clip(-np.sum(start*change, axis=-1)/np.maximum(np.sum(change*change, axis=-1), 1e-100), 0, 1)
            lower = length(start+f[:, None]*change)-sag
            bound = np.full(len(pairs), limit)
            if axial is not None:
                # Plane rotation over a 300-s interval is also bounded.
                depth = np.abs(axial[i]-axial[j])+travel
                bound = threshold+solar_angular_radius*depth+3e-7*dt*max_radius
            for ii, jj, distance in zip(i[lower < bound], j[lower < bound], lower[lower < bound]):
                key = (int(ii), int(jj))
                record[key] = min(record.get(key, float("inf")), float(distance))
        collect(mid, move, threshold, collision)
        if sun_directions is not None:
            u = unit(sun_directions[k]+sun_directions[k+1]); _, b, c = sail_basis(u)
            # Include off-duty and nightside members: shading their nonzero
            # beam-avoidance force would also change their subsequent orbits.
            ids = np.arange(len(mid))
            projected = np.column_stack([mid[ids]@b, mid[ids]@c])
            delta = np.column_stack([move[ids]@b, move[ids]@c])
            local = {}
            collect(projected, delta, threshold+2*solar_angular_radius*max_radius, local, mid[ids]@u)
            for (i, j), distance in local.items():
                key = (int(ids[i]), int(ids[j]))
                shadow[key] = min(shadow.get(key, float("inf")), distance)
    return collision, shadow


def conflict_free_subset(count, edges):
    """Deterministic greedy independent set; no claim of optimal inventory."""
    neighbours = [set() for _ in range(count)]
    for i, j in edges:
        neighbours[i].add(j); neighbours[j].add(i)
    remaining = set(range(count)); selected = []
    while remaining:
        i = min(remaining, key=lambda j: (len(neighbours[j] & remaining), j))
        selected.append(i)
        remaining.difference_update(neighbours[i] | {i})
    return np.sort(selected)
