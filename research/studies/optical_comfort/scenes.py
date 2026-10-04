"""Diffuse surface transport in explicitly prescribed metre-scale scenes.

Cosine-weighted Sobol paths integrate illuminance. The collimated Sun is
evaluated separately at each surface; the angular sky is sampled on escaping
paths. Photopic neutral Lambertian reflectance makes this transport scalar.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import math

import numpy as np
from scipy.stats import qmc

EPS = 1e-6  # ray-origin offset, metres; far below all scene features
SOBOL_BOUNCES = 18  # fixed dimension permits nested depth/sample comparisons


@dataclass(frozen=True)
class Rectangle:
    name: str
    axis: int
    coordinate: float
    low: tuple[float, float]
    high: tuple[float, float]
    reflectance: float


@dataclass(frozen=True)
class Scene:
    name: str
    ground_reflectance: float = 0.1
    patch_half_width: float = 4.0
    patch_reflectance: float = 0.1
    rectangles: tuple[Rectangle, ...] = ()

    def __post_init__(self):
        values = [self.ground_reflectance, self.patch_reflectance]
        for r in self.rectangles:
            if r.axis not in (0, 1, 2) or np.any(np.asarray(r.low) >= np.asarray(r.high)):
                raise ValueError('Invalid rectangle')
            values.append(r.reflectance)
        if not np.all(np.isfinite(values)) or min(values) < 0 or max(values) >= 1:
            raise ValueError('Reflectances must be finite and lie in [0, 1)')

    @property
    def maximum_reflectance(self):
        return max([self.ground_reflectance, self.patch_reflectance]
                   + [r.reflectance for r in self.rectangles])

    def record(self):
        return asdict(self)


def intersect(scene, origins, rays):
    """Nearest opaque surface; ID -1 is sky, 0 is ground, >0 is a rectangle.

    Rectangles are opaque and reflect from both faces with the same rho.
    The ground's finite pale patch changes reflectance without a coplanar mesh.
    """
    origins, rays = np.broadcast_arrays(np.asarray(origins, float), np.asarray(rays, float))
    count = len(rays)
    distance = np.full(count, np.inf)
    normal = np.zeros((count, 3))
    rho = np.zeros(count)
    ids = np.full(count, -1, dtype=int)
    t = np.divide(-origins[:, 2], rays[:, 2], out=np.full(count, np.inf), where=rays[:, 2] < -1e-14)
    hit = (t > EPS / 4) & np.isfinite(t)
    distance[hit] = t[hit]
    normal[hit, 2] = 1
    ids[hit] = 0
    xy = origins[hit, :2] + rays[hit, :2] * t[hit, None]
    patch = np.all(np.abs(xy) <= scene.patch_half_width, axis=1)
    rho[hit] = np.where(patch, scene.patch_reflectance, scene.ground_reflectance)
    for index, r in enumerate(scene.rectangles, 1):
        ax = r.axis
        other = [k for k in range(3) if k != ax]
        t = np.divide(r.coordinate - origins[:, ax], rays[:, ax],
                      out=np.full(count, np.inf), where=np.abs(rays[:, ax]) > 1e-14)
        candidates = np.flatnonzero((t > EPS / 4) & (t < distance))
        p = origins[candidates][:, other] + rays[candidates][:, other] * t[candidates, None]
        take = candidates[np.all((p >= r.low) & (p <= r.high), axis=1)]
        distance[take] = t[take]
        normal[take] = 0
        normal[take, ax] = -np.sign(rays[take, ax])
        rho[take] = r.reflectance
        ids[take] = index
    return distance, normal, rho, ids


def cosine_directions(samples, normals):
    normals = np.broadcast_to(normals, (len(samples), 3))
    helper = np.zeros_like(normals)
    helper[:, 2] = 1
    helper[np.abs(normals[:, 2]) > .9] = [1, 0, 0]
    tangent = np.cross(helper, normals)
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    bitangent = np.cross(normals, tangent)
    radius = np.sqrt(samples[:, 0])
    phi = 2 * math.pi * samples[:, 1]
    return (radius[:, None] * np.cos(phi)[:, None] * tangent
            + radius[:, None] * np.sin(phi)[:, None] * bitangent
            + np.sqrt(1 - samples[:, 0])[:, None] * normals)


def probe(scene, light, position, normal, power=17, seed=17, max_bounces=18):
    """Incident lux with separate sky, reflected-Sun and direct-Sun budgets.

    max_bounces counts reflected surface events, including ground. The omitted
    tail is estimated conservatively along surviving sampled paths using a
    radiance supersolution. Sampling error is assessed separately with a nested
    quarter-size sequence and a second scrambled sequence in the study runner.
    """
    if not 1 <= max_bounces <= SOBOL_BOUNCES or power < 4:
        raise ValueError('Unsupported integration settings')
    normal, position = np.asarray(normal, float), np.asarray(position, float)
    if not np.isclose(np.linalg.norm(normal), 1) or not np.all(np.isfinite(position)):
        raise ValueError('Finite position and unit normal required')
    count = 2 ** power
    samples = qmc.Sobol(d=2 * (SOBOL_BOUNCES + 1), scramble=True, seed=seed).random_base2(power)
    rays = cosine_directions(samples[:, :2], normal)
    origins = np.broadcast_to(position + EPS * normal, (count, 3)).copy()
    live = np.arange(count)
    throughput = np.full(count, math.pi)
    sky_values, sun_values = np.zeros(count), np.zeros(count)
    first = np.full(count, -2, dtype=int)
    tail = np.zeros(count)
    rho_max = scene.maximum_reflectance
    # If all incident radiance <= M, outgoing surface radiance <= rho*(M+DNI/pi).
    bound = max(light.sky.maximum * light.diffuse_h,
                rho_max * light.dni / (math.pi * (1 - rho_max)))
    for bounce in range(max_bounces + 1):
        t, normals, rho, ids = intersect(scene, origins, rays)
        if bounce == 0:
            first[:] = ids
        miss = ids == -1
        if np.any(miss):
            sky_values[live[miss]] += throughput[miss] * light.luminance(rays[miss])
        hit = ~miss
        if bounce == max_bounces:
            tail[live[hit]] = throughput[hit] * bound
            break
        live, throughput = live[hit], throughput[hit]
        if not len(live):
            break
        points = origins[hit] + rays[hit] * t[hit, None]
        normals, rho = normals[hit], rho[hit]
        origins = points + EPS * normals
        solar_cosine = np.maximum(normals @ light.sun, 0)
        shadow_t = intersect(scene, origins, np.broadcast_to(light.sun, origins.shape))[0]
        visible = ~np.isfinite(shadow_t)
        sun_values[live] += throughput * rho / math.pi * light.dni * solar_cosine * visible
        throughput *= rho
        rays = cosine_directions(samples[live, 2 * (bounce + 1):2 * (bounce + 2)], normals)
        # Black paths contribute exactly zero at every later bounce.
        keep = throughput > 0
        live, throughput, rays, origins = live[keep], throughput[keep], rays[keep], origins[keep]
        if not len(live):
            break
    sun_visible = not np.isfinite(intersect(scene, (position + EPS * normal)[None], light.sun[None])[0][0])
    direct = light.dni * max(float(normal @ light.sun), 0) * sun_visible
    all_values = sky_values + sun_values
    result = dict(
        total_lux=float(np.mean(all_values)) + direct,
        ambient_lux=float(np.mean(all_values)), sun_direct_lux=direct,
        sky_unobstructed_lux=float(np.mean(np.where(first == -1, sky_values, 0))),
        sky_reflected_lux=float(np.mean(np.where(first != -1, sky_values, 0))),
        sun_reflected_lux=float(np.mean(sun_values)),
        ground_in_view_lux=float(np.mean(np.where(first == 0, all_values, 0))),
        structures_in_view_lux=float(np.mean(np.where(first > 0, all_values, 0))),
        coarse_total_lux=float(np.mean(all_values[:count // 4])) + direct,
        sampled_tail_upper_lux=float(np.mean(tail)),
    )
    return result


def scene_set():
    roof = Rectangle('roof', 2, 3., (-3., -3.), (3., 3.), .3)
    def walls(half, rho):
        return (Rectangle('left_wall', 0, -half, (-half, 0.), (half, 3.), rho),
                Rectangle('right_wall', 0, half, (-half, 0.), (half, 3.), rho),
                Rectangle('back_wall', 1, half, (-half, 0.), (half, 3.), rho))
    return [Scene('open_dark'), Scene('open_pale_patch', patch_reflectance=.8),
            Scene('pale_courtyard', patch_reflectance=.8, rectangles=walls(4., .8)),
            Scene('canopy_dark', rectangles=(roof,)),
            Scene('canopy_pale', patch_reflectance=.8, rectangles=(roof,)),
            Scene('screened_canopy_pale', patch_reflectance=.8, rectangles=(roof,) + walls(3., .3))]


def with_step(scene):
    step = (Rectangle('step_tread', 2, .15, (-1., -2.), (1., -1.), .2),
            Rectangle('step_riser', 1, -1., (-1., 0.), (1., .15), .2),
            Rectangle('step_back', 1, -2., (-1., 0.), (1., .15), .2),
            Rectangle('step_left', 0, -1., (-2., 0.), (-1., .15), .2),
            Rectangle('step_right', 0, 1., (-2., 0.), (-1., .15), .2))
    return replace(scene, name=scene.name + '_with_step', rectangles=scene.rectangles + step)


PROBES = {
    'eye_level': ((0., 0., 1.6), (0., -1., 0.)),
    'eye_to_workplane': ((0., 0., 1.6), tuple(np.array([0., -1., -.8]) / np.sqrt(1.64))),
    'ground_centre': ((0., 0., 0.), (0., 0., 1.)),
    'workplane': ((0., -1., .8), (0., 0., 1.)),
    'step_tread': ((0., -1.5, .15), (0., 0., 1.)),
    'step_riser': ((0., -1., .075), (0., 1., 0.)),
}
