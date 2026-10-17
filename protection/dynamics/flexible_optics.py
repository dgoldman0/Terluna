"""Direct photon loading and sampled coverage on actual deformed triangles.

Each call uses the supplied world vertices. Mesh normals, interception area,
mutual transmission and receiver coverage therefore change with deformation.
Area integration is one centroid sample per triangle. It is exact for a flat,
uniformly illuminated face, but shadow boundaries require mesh refinement.
Neither mesh/source sampling nor a positive sampled coverage result is a
continuous-space or continuous-time coverage guarantee.

Ray intersections use the Moller--Trumbore barycentric construction (1997,
doi:10.1080/10867651.1997.10487468; paper mirrored at
https://my.eng.utah.edu/~cs6965/papers/MT97.pdf). The photon law is direct
momentum balance for absorption, specular reflection and straight transmission.
NASA's richer reference model also includes diffuse reflection and thermal
emission (https://ntrs.nasa.gov/citations/20150021404); those effects, secondary
reflected rays, spectral variation and hardware not explicitly meshed are
outside this initial model. Optical coefficients are supplied assumptions.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import cos, isfinite, pi, sin, sqrt

import numpy as np


def _vectors(value, name):
    a = np.asarray(value, dtype=float)
    if a.ndim == 1:
        a = a[None, :]
    if a.ndim != 2 or a.shape[1] != 3 or not len(a) or not np.isfinite(a).all():
        raise ValueError(f"{name} must be a finite nonempty (N,3) array")
    return a


def _unit_vectors(value, name):
    a = _vectors(value, name)
    if not np.allclose(np.linalg.norm(a, axis=1), 1., rtol=1e-10, atol=1e-12):
        raise ValueError(f"{name} must contain unit vectors")
    return a


@dataclass(frozen=True)
class MeshTile:
    """Two-sided, zero-thickness triangle filter with scalar/per-face fractions.

    The same supplied optical properties apply to both sides. Face winding
    changes no physical force: its normal is oriented toward each source ray.
    A mesh may fold, so a ray can cross one tile more than once at distinct
    distances. Coincident triangles of one tile count as one optical crossing.
    """
    vertices: object
    triangles: object
    absorptivity: object
    reflectivity: object
    transmissivity: object

    def __post_init__(self):
        vertices = _vectors(self.vertices, 'vertices')
        triangles = np.asarray(self.triangles)
        if (triangles.ndim != 2 or triangles.shape[1] != 3 or not len(triangles)
                or not np.issubdtype(triangles.dtype, np.integer)
                or (triangles < 0).any() or (triangles >= len(vertices)).any()):
            raise ValueError("triangles must be a nonempty (F,3) array of valid integer vertex indices")
        object.__setattr__(self, 'vertices', vertices)
        object.__setattr__(self, 'triangles', triangles.astype(np.intp, copy=False))
        coefficients = []
        for name in ('absorptivity', 'reflectivity', 'transmissivity'):
            value = np.asarray(getattr(self, name), dtype=float)
            if value.ndim == 0:
                value = np.full(len(triangles), float(value))
            if (value.shape != (len(triangles),) or not np.isfinite(value).all()
                    or (value < 0).any() or (value > 1).any()):
                raise ValueError(f"{name} must be a finite scalar or per-face fraction in [0,1]")
            object.__setattr__(self, name, value)
            coefficients.append(value)
        if not np.allclose(sum(coefficients), 1., rtol=0., atol=1e-12):
            raise ValueError("Absorption, reflection and transmission must sum to one on each face")


def solar_disk_quadrature(direction, angular_radius_rad: float, rings: int = 1):
    """Positive central/ring quadrature for a uniform-radiance spherical cap.

    Returns (directions, weights). Counts for rings=0,1,2 are 1,7,19; general
    count is 1+3*rings*(rings+1). Each ring j has 6j equal-azimuth samples.
    Its equal-solid-angle annulus has weight 6j/count, sampled at mid area.
    The central cap has weight 1/count, sampled at the centre. This simple
    convergent construction is not a high-order or limb-darkened solar rule.

    Weights sum to one, describing relative directional intensity. A supplied
    pressure S/c then means the sum of ray normal-intensities/c. For exact
    irradiance S on a plane normal to the cap centre, use pressure
    (S/c)/sum(weights*(directions@centre)). This is a tiny solar-angle correction.
    """
    centre = _unit_vectors(direction, 'direction')
    if len(centre) != 1:
        raise ValueError("direction must be one unit vector")
    centre = centre[0]
    if not isfinite(angular_radius_rad) or not 0 <= angular_radius_rad < pi/2:
        raise ValueError("angular_radius_rad must be finite and in [0,pi/2)")
    if isinstance(rings, bool) or not isinstance(rings, int) or rings < 0:
        raise ValueError("rings must be a nonnegative integer")
    if rings == 0 or angular_radius_rad == 0:
        return centre[None, :].copy(), np.ones(1)
    # Choose the least-aligned Cartesian axis for a well-conditioned basis.
    reference = np.eye(3)[np.argmin(np.abs(centre))]
    first = np.cross(centre, reference)
    first /= np.linalg.norm(first)
    second = np.cross(centre, first)
    count = 1+3*rings*(rings+1)
    directions = [centre]
    weights = [1/count]
    cumulative = 1/count
    cap_height = 2*sin(angular_radius_rad/2)**2
    for ring in range(1, rings+1):
        points = 6*ring
        fraction = points/count
        mid = cumulative+fraction/2
        z = 1-mid*cap_height
        radial = sqrt(max(0., (1-z)*(1+z)))
        # Stagger successive rings to avoid long aligned rows of source nodes.
        phase = (ring % 2)*pi/points
        for index in range(points):
            phi = 2*pi*index/points+phase
            directions.append(z*centre+radial*(cos(phi)*first+sin(phi)*second))
            weights.append(1/count)
        cumulative += fraction
    return np.asarray(directions), np.asarray(weights)


def _scene(tiles):
    tiles = list(tiles)
    if any(not isinstance(tile, MeshTile) for tile in tiles):
        raise ValueError("tiles must contain MeshTile instances")
    if not tiles:
        return dict(tiles=[], points=np.empty((0, 3, 3)), tile=np.empty(0, dtype=int),
                    face=np.empty(0, dtype=int), area=np.empty(0), normal=np.empty((0, 3)),
                    centre=np.empty((0, 3)), absorption=np.empty(0),
                    reflection=np.empty(0), transmission=np.empty(0),
                    bounds_min=np.empty((0, 3)), bounds_max=np.empty((0, 3)))
    points = np.concatenate([tile.vertices[tile.triangles] for tile in tiles])
    if not np.isfinite(points).all():
        raise ValueError("Updated mesh positions must remain finite")
    first, second = points[:, 1]-points[:, 0], points[:, 2]-points[:, 0]
    cross = np.cross(first, second)
    magnitude = np.linalg.norm(cross, axis=1)
    scale = np.linalg.norm(first, axis=1)*np.linalg.norm(second, axis=1)
    if (magnitude <= 1e-14*scale).any() or (scale == 0).any():
        raise ValueError("A triangle is collapsed or numerically degenerate")
    return dict(tiles=tiles, points=points,
                tile=np.concatenate([np.full(len(tile.triangles), i, dtype=int) for i, tile in enumerate(tiles)]),
                face=np.concatenate([np.arange(len(tile.triangles)) for tile in tiles]),
                area=magnitude/2, normal=cross/magnitude[:, None], centre=points.mean(axis=1),
                absorption=np.concatenate([tile.absorptivity for tile in tiles]),
                reflection=np.concatenate([tile.reflectivity for tile in tiles]),
                transmission=np.concatenate([tile.transmissivity for tile in tiles]),
                bounds_min=np.asarray([tile.vertices.min(axis=0) for tile in tiles]),
                bounds_max=np.asarray([tile.vertices.max(axis=0) for tile in tiles]))


def _ray_aabb_candidates(origins, directions, minimum, maximum, epsilon):
    """Ray/slab intersection for each ray and tile; includes boundary tolerance."""
    delta0 = minimum[None, :, :]-epsilon-origins[:, None, :]
    delta1 = maximum[None, :, :]+epsilon-origins[:, None, :]
    direction = directions[:, None, :]
    nonzero = np.abs(direction) > 1e-15
    lower = np.divide(delta0, direction, out=np.zeros_like(delta0), where=nonzero)
    upper = np.divide(delta1, direction, out=np.zeros_like(delta1), where=nonzero)
    entry, leave = np.minimum(lower, upper), np.maximum(lower, upper)
    inside = (delta0 <= 0) & (delta1 >= 0)
    entry = np.where(nonzero, entry, np.where(inside, -np.inf, np.inf))
    leave = np.where(nonzero, leave, np.where(inside, np.inf, -np.inf))
    return (leave.min(axis=2) >= np.maximum(entry.max(axis=2), 0.))


def _trace(scene, origins, directions, *, ray_epsilon_m, merge_tolerance_m,
           chunk_size, triangle_chunk_size, exclude_faces=None, exclude_tiles=None):
    """Chunked ray/triangle tests, then distinct crossing transmission products."""
    if (not isfinite(ray_epsilon_m) or ray_epsilon_m < 0 or
            not isfinite(merge_tolerance_m) or merge_tolerance_m < 0):
        raise ValueError("Ray and crossing-merge tolerances must be finite and nonnegative")
    if any(isinstance(x, bool) or not isinstance(x, int) or x < 1
           for x in (chunk_size, triangle_chunk_size)):
        raise ValueError("Chunk sizes must be positive integers")
    origins = _vectors(origins, 'origins')
    directions = _unit_vectors(directions, 'directions')
    try:
        origins, directions = np.broadcast_arrays(origins, directions)
    except ValueError as exc:
        raise ValueError("Ray origins and directions must have matching or singleton counts") from exc
    count = len(origins)
    first_distance = np.full(count, np.inf)
    first_tile = np.full(count, -1, dtype=int)
    first_face = np.full(count, -1, dtype=int)
    transmittance = np.ones(count)
    crossings = np.zeros(count, dtype=int)
    filtered = np.zeros(count, dtype=bool)
    faces = scene['points']
    for start in range(0, count, chunk_size):
        stop = min(count, start+chunk_size)
        hits = [[] for _ in range(stop-start)]
        o, d = origins[start:stop], directions[start:stop]
        candidates = _ray_aabb_candidates(o, d, scene['bounds_min'], scene['bounds_max'],
                                          max(ray_epsilon_m, merge_tolerance_m))
        for fstart in range(0, len(faces), triangle_chunk_size):
            fstop = min(len(faces), fstart+triangle_chunk_size)
            proposed = np.arange(fstart, fstop)
            candidate_pairs = candidates[:, scene['tile'][proposed]]
            kept_faces = np.any(candidate_pairs, axis=0)
            if not kept_faces.any():
                continue
            indices = proposed[kept_faces]
            rows_kept = np.flatnonzero(np.any(candidate_pairs[:, kept_faces], axis=1))
            block = faces[indices]
            selected_o, selected_d = o[rows_kept], d[rows_kept]
            e1, e2 = block[:, 1]-block[:, 0], block[:, 2]-block[:, 0]
            pvec = np.cross(selected_d[:, None, :], e2[None, :, :])
            determinant = np.einsum('rfi,fi->rf', pvec, e1)
            det_scale = np.linalg.norm(e1, axis=1)*np.linalg.norm(e2, axis=1)
            valid = np.abs(determinant) > 1e-13*det_scale[None, :]
            inverse = np.divide(1., determinant, out=np.zeros_like(determinant), where=valid)
            tvec = selected_o[:, None, :]-block[None, :, 0, :]
            u = np.einsum('rfi,rfi->rf', tvec, pvec)*inverse
            qvec = np.cross(tvec, e1[None, :, :])
            v = np.einsum('ri,rfi->rf', selected_d, qvec)*inverse
            distance = np.einsum('fi,rfi->rf', e2, qvec)*inverse
            barycentric_tolerance = 1e-12
            valid &= ((u >= -barycentric_tolerance) & (v >= -barycentric_tolerance)
                      & (u+v <= 1+barycentric_tolerance) & (distance > ray_epsilon_m))
            valid &= candidates[rows_kept[:, None], scene['tile'][indices][None, :]]
            if exclude_faces is not None:
                valid &= indices[None, :] != exclude_faces[start+rows_kept, None]
            if exclude_tiles is not None:
                valid &= scene['tile'][indices][None, :] != exclude_tiles[start+rows_kept, None]
            rows, columns = np.nonzero(valid)
            for r, c in zip(rows, columns):
                hits[rows_kept[r]].append((float(distance[r, c]), int(indices[c])))
        for local, ray_hits in enumerate(hits):
            index = start+local
            if not ray_hits:
                continue
            ray_hits.sort()
            first_distance[index] = ray_hits[0][0]
            first_tile[index] = scene['tile'][ray_hits[0][1]]
            first_face[index] = scene['face'][ray_hits[0][1]]
            groups = []
            latest = {}
            for distance, face in ray_hits:
                tile = int(scene['tile'][face])
                fraction = float(scene['transmission'][face])
                previous = latest.get(tile)
                if previous is not None and abs(distance-groups[previous][0]) <= merge_tolerance_m:
                    # Shared-edge coefficients may differ. Use the lower
                    # transmission on the measure-zero seam, never multiply twice.
                    groups[previous][1] = min(groups[previous][1], fraction)
                else:
                    latest[tile] = len(groups)
                    groups.append([distance, fraction])
            crossings[index] = len(groups)
            transmittance[index] = float(np.prod([entry[1] for entry in groups]))
            filtered[index] = any(entry[1] < 1 for entry in groups)
    return dict(geometric_covered=crossings > 0, filtered_covered=filtered,
                first_hit_distance_m=first_distance, first_hit_tile=first_tile,
                first_hit_face=first_face, crossing_count=crossings,
                energy_transmittance=transmittance,
                reading_rule="Direct-ray scalar transmission; coverage means sampled intersections, not a UV guarantee. Coincident faces within each tile merge using minimum transmission.")


def receiver_rays(tiles, origins, directions, *, ray_epsilon_m=1e-7,
                  merge_tolerance_m=1e-7, chunk_size=128, triangle_chunk_size=256):
    """Trace receiver-to-source rays through all distinct mesh crossings.

    The same optical band must be used across tiles. An opaque UV model may
    use zero transmissivity; broadband coefficients cannot certify UV removal.
    Geometric coverage, any attenuating crossing and scalar transmitted energy
    are deliberately separate outputs. Reflections into a receiver are omitted.
    """
    return _trace(_scene(tiles), origins, directions, ray_epsilon_m=ray_epsilon_m,
                  merge_tolerance_m=merge_tolerance_m, chunk_size=chunk_size,
                  triangle_chunk_size=triangle_chunk_size)


def photon_forces(tiles, source_directions, source_weights, pressure_Pa: float, *,
                  torque_origin=None, self_shadow=True, ray_epsilon_m=1e-7,
                  merge_tolerance_m=1e-7, chunk_size=128, triangle_chunk_size=256,
                  external_transmission=None):
    """Centroid photon forces with upstream multilayer transmission feedback.

    Source directions point *toward* the source. For n oriented toward that
    ray, c=n.s and incoming k=-s, F=(S/c_light)*A*c*(a*k-2*r*c*n).
    Weights and upstream straight transmission multiply this force. Face force
    is equally distributed to its three vertices: total force, torque and work
    against affine triangle velocities are conserved by this discretization.
    It does not integrate partly shadowed faces exactly.

    Optional external_transmission(centres, direction) supplies a finite (F,)
    array in [0,1] for attenuation outside this mesh scene, such as eclipsing
    bodies. It is evaluated on the actual world centroids for each positive-
    weight source sample, and multiplies the mesh crossing transmission.
    """
    if not isfinite(pressure_Pa) or pressure_Pa < 0:
        raise ValueError("pressure_Pa must be finite and nonnegative")
    source = _unit_vectors(source_directions, 'source_directions')
    weights = np.asarray(source_weights, dtype=float)
    if (weights.shape != (len(source),) or not np.isfinite(weights).all()
            or (weights < 0).any() or not np.isclose(weights.sum(), 1., rtol=0., atol=1e-12)):
        raise ValueError("Source weights must be finite nonnegative and sum to one")
    if not isinstance(self_shadow, (bool, np.bool_)):
        raise ValueError("self_shadow must be boolean")
    if external_transmission is not None and not callable(external_transmission):
        raise ValueError("external_transmission must be callable or None")
    origin = np.zeros(3) if torque_origin is None else np.asarray(torque_origin, dtype=float)
    if origin.shape != (3,) or not np.isfinite(origin).all():
        raise ValueError("torque_origin must be a finite three-vector")
    scene = _scene(tiles)
    count = len(scene['points'])
    force = np.zeros((count, 3))
    reflection_derivative = np.zeros((count, 3))
    incidence = np.zeros(count)
    attenuated_incidence = np.zeros(count)
    unoccluded = np.zeros(count)
    direct_weight = np.zeros(count)
    external_weight = np.zeros(count)
    if count:
        for direction, weight in zip(source, weights):
            if weight == 0:
                continue
            trace = _trace(scene, scene['centre'], direction,
                           ray_epsilon_m=ray_epsilon_m, merge_tolerance_m=merge_tolerance_m,
                           chunk_size=chunk_size, triangle_chunk_size=triangle_chunk_size,
                           exclude_faces=np.arange(count),
                           exclude_tiles=None if self_shadow else scene['tile'])
            external = (np.ones(count) if external_transmission is None else
                        np.asarray(external_transmission(scene['centre'], direction), dtype=float))
            if (external.shape != (count,) or not np.isfinite(external).all()
                    or (external < 0).any() or (external > 1).any()):
                raise ValueError("External transmission must return a finite (F,) array in [0,1]")
            transmission = trace['energy_transmittance']*external
            dot = scene['normal']@direction
            cosine = np.abs(dot)
            oriented = scene['normal']*np.where(dot >= 0, 1., -1.)[:, None]
            impulse = (-scene['absorption'][:, None]*direction
                       -2*scene['reflection'][:, None]*cosine[:, None]*oriented)
            force += (pressure_Pa*weight*transmission*scene['area']*cosine)[:, None]*impulse
            reflection_derivative += (-2*pressure_Pa*weight*transmission*scene['area']*cosine**2)[:, None]*oriented
            incidence += weight*cosine
            attenuated_incidence += weight*transmission*cosine
            unoccluded += weight*(~trace['geometric_covered'])*(external == 1.)
            direct_weight += weight*transmission
            external_weight += weight*external
    nodal, face_forces, face_incidence, source_transmission, unoccluded_weights = [], [], [], [], []
    attenuated_face_incidence, tile_forces, tile_torques = [], [], []
    derivatives, face_normals, face_areas = [], [], []
    external_transmissions = []
    nodal_total = np.zeros(3)
    nodal_torque = np.zeros(3)
    for index, tile in enumerate(scene['tiles']):
        selected = scene['tile'] == index
        local_force = force[selected]
        nodes = np.zeros_like(tile.vertices)
        for vertex_index in range(3):
            np.add.at(nodes, tile.triangles[:, vertex_index], local_force/3)
        local_sum = nodes.sum(axis=0)
        local_torque = np.cross(tile.vertices-origin, nodes).sum(axis=0)
        nodal.append(nodes)
        face_forces.append(local_force)
        face_incidence.append(incidence[selected])
        attenuated_face_incidence.append(attenuated_incidence[selected])
        source_transmission.append(direct_weight[selected])
        external_transmissions.append(external_weight[selected])
        unoccluded_weights.append(unoccluded[selected])
        derivatives.append(reflection_derivative[selected])
        face_normals.append(scene['normal'][selected])
        face_areas.append(scene['area'][selected])
        tile_forces.append(local_sum)
        tile_torques.append(local_torque)
        nodal_total += local_sum
        nodal_torque += local_torque
    face_total = force.sum(axis=0)
    face_torque = np.cross(scene['centre']-origin, force).sum(axis=0)
    return dict(nodal_forces=nodal, face_forces=face_forces,
                face_incidence=face_incidence,
                face_attenuated_incidence=attenuated_face_incidence,
                face_direct_source_transmission=source_transmission,
                face_external_source_transmission=external_transmissions,
                external_transmission_applied=external_transmission is not None,
                transmission_rule="Source-weighted direct transmission combines distinct mesh crossings and the external callback. External transmission is also reported separately; unoccluded weights require no mesh crossing and external transmission equal to one.",
                face_unoccluded_source_weights=unoccluded_weights,
                reflection_force_derivative=derivatives,
                reflection_derivative_rule="Local derivative at fixed absorption and fixed upstream transmission; changing reflection also changes downstream transmission, requiring full recomputation.",
                face_normals=face_normals, face_areas=face_areas,
                tile_forces=tile_forces, tile_torques=tile_torques,
                total_force=nodal_total, total_torque=nodal_torque,
                torque_origin=origin, force_balance_residual_N=nodal_total-face_total,
                torque_balance_residual_N_m=nodal_torque-face_torque,
                source_samples=len(source), face_samples=count,
                evidence="direct momentum balance on deformed triangle centroids with sampled mutual transmission",
                unresolved=["mesh/source refinement at shadow boundaries and seams",
                            "secondary reflected rays, diffuse scattering and thermal recoil",
                            "spectral/angle-dependent optical properties and unmeshed hardware",
                            "continuous-time coverage and self-contact"])
