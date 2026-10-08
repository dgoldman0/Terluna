"""Near-natural orbit libraries and a finite-square area-capacity relaxation.

Fractional replication preserves the phase measure of each column. Its summed
areas are necessary cell-capacity constraints, not union coverage. No density
interpolation is used to stand in for independently propagated phases.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linprog

from shared import constants as K
from .cycling import sail_basis, ray_geometry
from .fleet import square_vertices, projected_polygons
from .optical import arriving_ray_vectors, length, unit


def seed_library(env, radii, azimuths, tilts, phases):
    sample = env.at(0.)
    sun, _ = arriving_ray_vectors(sample['positions']['sun'], sample['positions']['earth'],
        sample['sun_v'], sample['earth_v'], sample['sun_a'], sample['earth_a'])
    u = unit(sun); _, b, c = sail_basis(sun)
    initial, descriptions = [], []
    golden = (np.sqrt(5)-1)/2
    for radius in radii:
        for azimuth in azimuths:
            for tilt in tilts:
                h = np.cos(tilt)*(np.cos(azimuth)*b+np.sin(azimuth)*c)+np.sin(tilt)*u
                if np.dot(h, env.frames(0.)[:, 2]) > 0:
                    h = -h
                a = unit(u-h*np.dot(h, u)); d = np.cross(h, a)
                offset = (len(descriptions)*golden) % 1
                angle = 2*np.pi*(np.arange(phases)/phases+offset)
                q = radius*(np.cos(angle)[:, None]*a+np.sin(angle)[:, None]*d)
                v = np.sqrt(K.MOON_GM/radius)*(-np.sin(angle)[:, None]*a+np.cos(angle)[:, None]*d)
                initial.append(np.c_[q, v])
                descriptions.append(dict(radius_m=float(radius), azimuth_rad=float(azimuth),
                    tilt_rad=float(tilt), phase_offset_turns=offset))
    return np.concatenate(initial), descriptions


def receiver_cells(radius, rings=3, sectors=8):
    """Inscribed polygon partition, with exact analytic area error reported."""
    from shapely.geometry import Polygon
    cells = []
    for ring in range(rings):
        inner, outer = radius*np.sqrt(np.array([ring, ring+1])/rings)
        for sector in range(sectors):
            angle = np.linspace(sector*2*np.pi/sectors, (sector+1)*2*np.pi/sectors, 65)
            outer_xy = outer*np.c_[np.cos(angle), np.sin(angle)]
            inner_xy = inner*np.c_[np.cos(angle[::-1]), np.sin(angle[::-1])]
            cells.append(Polygon(np.r_[outer_xy, inner_xy]))
    return np.array(cells, dtype=object)


def sources_with_limb(sun, limb_count=4, rotation=0.):
    _, b, c = sail_basis(sun)
    angles = rotation+np.arange(limb_count)*2*np.pi/limb_count
    return np.r_[sun[None, :], sun+K.SUN_RADIUS*(np.cos(angles)[:, None]*b+np.sin(angles)[:, None]*c)]


def capacity_rows(state, normals, sample, cells, families, phases, side=10000., limb_count=4):
    """Cell capacity per million tiles; physical duplicates are a relaxation.

    Each shadow is a true 9.89 km clear square, with first-order interception
    motion. Overlap is summed and available cell area removes Earth eclipse.
    Central plus limb sources are constraints, not an irradiance quadrature.
    """
    import shapely
    from shapely.geometry import Point
    sun, _ = arriving_ray_vectors(sample['positions']['sun'], sample['positions']['earth'],
        sample['sun_v'], sample['earth_v'], sample['sun_a'], sample['earth_a'])
    u = unit(sun); _, b, c = sail_basis(sun)
    q = state[:, :3]-(state[:, 3:]+sample['moon_v'])*np.maximum(state[:, :3]@u, 0)[:, None]/K.SPEED_OF_LIGHT
    keep = q@u > 0
    ids = np.flatnonzero(keep)
    # Use the same tile-local retarded incident vector to set roll as the
    # finite-extent force quadrature. The target-plane frame stays Moon-centred.
    tile_sun = ray_geometry(state[:, :3], sample)['sun']
    corners = square_vertices(q[keep], normals[keep], tile_sun[keep], side-110.)
    rows, labels = [], []
    earth = sample['positions']['earth']
    earth = earth-sample['earth_v']*max(np.dot(earth, u), 0)/K.SPEED_OF_LIGHT
    for source_id, source in enumerate(sources_with_limb(sun, limb_count)):
        polys = projected_polygons(corners, source, u, b, c)
        available = cells
        depth = earth@u
        if 0 < depth < source@u:
            mag = (source@u)/(source@u-depth)
            centre = source+(earth-source)*mag
            eclipse = Point(centre@b, centre@c).buffer(K.EARTH_RADIUS*mag, quad_segs=256)
            available = shapely.difference(cells, eclipse)
        for cell_id, cell in enumerate(available):
            if cell.area < 1.:
                continue
            area = shapely.area(shapely.intersection(polys, cell))
            summed = np.bincount(ids//phases, weights=area, minlength=families)
            rows.append(1e6*summed/phases/cell.area)
            labels.append([source_id, cell_id])
    return np.asarray(rows), labels


def cover_capacity(matrix, admissible=None):
    """Minimize millions of tiles with A x >= 1; preserve uncovered rows."""
    matrix = np.asarray(matrix, float)
    if admissible is None:
        admissible = np.ones(matrix.shape[1], bool)
    ids = np.flatnonzero(admissible)
    empty = np.flatnonzero(np.max(matrix[:, ids], axis=1) <= 1e-15) if len(ids) else np.arange(len(matrix))
    if len(empty):
        return dict(success=False, reason='empty_capacity_rows', empty_rows=empty.tolist())
    result = linprog(np.ones(len(ids)), A_ub=-matrix[:, ids], b_ub=-np.ones(len(matrix)),
                     bounds=(0, None), method='highs', options={'time_limit': 30.})
    if not result.success:
        return dict(success=False, reason=result.message)
    x = np.zeros(matrix.shape[1]); x[ids] = result.x
    achieved = matrix@x
    dual = -result.ineqlin.marginals
    bottlenecks = np.argsort(dual)[-12:][::-1]
    return dict(success=True, inventory_tiles=float(x.sum()*1e6), weights_million_tiles=x.tolist(),
        minimum_capacity=float(achieved.min()), dual_lower_bound_million_tiles=float(dual.sum()),
        dual_max_column_price=float(np.max(matrix[:, ids].T@dual)),
        bottleneck_rows=bottlenecks.tolist(), bottleneck_duals=dual[bottlenecks].tolist())


def conditional_cost(tiles, dv_per_day, prop, hardware_multiplier=1., duty=1.):
    """Explicit scenario, including seven-day fuel and installed-power feedback.

    The multiplier supplies fixed extra vehicle mass relative to optical mass.
    Constant recurring impulse and a selected peak/mean duty are assumptions.
    This scalar closure does not repropagate the changed sail acceleration.
    """
    mass0 = tiles*10000.**2*.05
    cant = np.cos(np.deg2rad(prop['cant_deg']))
    a = dv_per_day/K.JULIAN_DAY
    per_mass_power = a*prop['exhaust_velocity_m_s']/(2*prop['efficiency']*cant)
    power_fraction = prop['peak_margin_factor']*per_mass_power/(duty*prop['specific_power_W_kg'])
    fuel_fraction = a*prop['propellant_buffer_days']*K.JULIAN_DAY/(prop['exhaust_velocity_m_s']*cant)
    denominator = 1-power_fraction-fuel_fraction
    if denominator <= 0:
        return dict(closed=False, denominator=float(denominator))
    mass = mass0*hardware_multiplier/denominator
    return dict(closed=True, optical_mass_kg=float(mass0), total_scenario_mass_kg=float(mass),
        recurring_dv_m_s_per_day=float(dv_per_day), peak_duty=float(duty),
        extra_fixed_mass_multiplier=float(hardware_multiplier),
        propulsion_TW=float(mass*per_mass_power/1e12),
        propellant_kg_s=float(mass*a/(prop['exhaust_velocity_m_s']*cant)),
        installed_propulsion_power_TW=float(prop['peak_margin_factor']*mass*per_mass_power/duty/1e12))


def orientation_matrices(state, normal, sun, side=10000.):
    corners = square_vertices(state[:, :3], normal, sun, side)
    a = (corners[:, 1]-corners[:, 0])/side
    b = (corners[:, 3]-corners[:, 0])/side
    return np.stack([a, b, normal], axis=-1)


def rotation_step_angle(before, after):
    """Minimum rotation of an unlabelled, two-sided ideal square.

    Four rolls and the two normal signs describe the same square. Minimizing
    over these symmetries avoids counting a basis branch flip as a slew.
    This is a necessary rate bound; labelled hardware or asymmetric coatings
    can require a longer motion. Safe intermediate reflected beams are separate.
    """
    minimum = np.full(np.asarray(before).shape[:-2], np.inf)
    for sign in [-1, 1]:
        for k in range(4):
            angle = k*np.pi/2
            c, s = np.cos(angle), np.sin(angle)
            symmetry = np.array([[c,-s,0], [s,c,0], [0,0,1.]])@np.diag([1,sign,sign])
            cosine = np.clip((np.sum(before*(after@symmetry), axis=(-2,-1))-1)/2, -1, 1)
            minimum = np.minimum(minimum, np.arccos(cosine))
    return minimum
