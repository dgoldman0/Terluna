"""Reflection by a wind-roughened water surface: Fresnel facets with Gaussian slopes and Smith's shadowing.

The surface is a field of facets whose gradients follow a two-dimensional Gaussian with a given east-north
covariance; Cox and Munk (1954) measured the sea's slopes close to Gaussian, and their small skewness and
peakedness are left out here. Each facet reflects by Fresnel's law for unpolarized light. Smith's (1967)
shadowing function for a Gaussian surface gives the share of facets hidden from the observer and of reflected
rays that strike another wave.

Two results follow. A small bright source, the Sun or the Earth, appears as glitter:
L = E_n rho(omega) p(g) G / (4 cos(theta_v) cos^4(beta)) (Cox and Munk 1954), with E_n its irradiance normal to
the beam, g the facet gradient that mirrors it to the observer, beta that facet's tilt and G the shadowing. The
sky reflects as the same kernel integrated over the facets the observer sees. A reflected ray that points into
the water or meets another wave reflects once more off a level surface; Mobley (2015, Appl. Opt. 54, 4828)
finds multiple reflection a few percent of the reflected light outside grazing views.

Directions are unit vectors (east, north, up) pointing away from the surface: v toward the observer, s toward
the source.
"""
from __future__ import annotations

import numpy as np
from scipy.special import erfc

SQRT_PI = np.sqrt(np.pi)


def refractive_index(wavelength_nm, temperature_c=20.0, salinity=35.0):
    """Seawater relative to air, Quan and Fry (1995, Appl. Opt. 34, 3477): 400-700 nm, 0-30 C, salinity 0-35."""
    w = np.asarray(wavelength_nm, float)
    t, s = temperature_c, salinity
    return (1.31405 + (1.779e-4 - 1.05e-6 * t + 1.6e-8 * t * t) * s - 2.02e-6 * t * t
            + (15.868 + 0.01155 * s - 0.00423 * t) / w - 4382.0 / w ** 2 + 1.1455e6 / w ** 3)


def fresnel(cos_incidence, n):
    """Unpolarized Fresnel reflectance of a facet from air into water of index n."""
    c = np.clip(np.asarray(cos_incidence, float), 0.0, 1.0)
    root = np.sqrt(np.maximum(n * n - 1 + c * c, 0.0))
    rs = ((c - root) / (c + root)) ** 2
    rp = ((n * n * c - root) / (n * n * c + root)) ** 2
    return 0.5 * (rs + rp)


def smith_lambda(direction, covariance):
    """Smith's Lambda for rays leaving the surface in the given directions (shape (..., 3))."""
    d = np.asarray(direction, float)
    horizontal = np.hypot(d[..., 0], d[..., 1])
    with np.errstate(divide="ignore", invalid="ignore"):
        unit = np.stack([d[..., 0], d[..., 1]], -1) / horizontal[..., None]
        sigma = np.sqrt(np.einsum("...i,ij,...j->...", unit, covariance, unit))
        a = d[..., 2] / horizontal / (np.sqrt(2) * sigma)
        value = (np.exp(-a * a) - a * SQRT_PI * erfc(a)) / (2 * a * SQRT_PI)
    return np.where(horizontal < 1e-12, 0.0, np.where(d[..., 2] > 0, value, np.inf))


def slope_density(gradient, covariance):
    g = np.asarray(gradient, float)
    inverse = np.linalg.inv(covariance)
    return np.exp(-0.5 * np.einsum("...i,ij,...j->...", g, inverse, g)) / (2 * np.pi * np.sqrt(np.linalg.det(covariance)))


def glint(source, view, covariance, n):
    """Glitter factor K: a source of irradiance E_n normal to its beam appears with radiance E_n K."""
    s, v = np.asarray(source, float), np.asarray(view, float)
    half = s + v
    half = half / np.linalg.norm(half, axis=-1, keepdims=True)
    up = (s[..., 2] > 0) & (v[..., 2] > 0) & (half[..., 2] > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        gradient = -half[..., :2] / half[..., 2:3]
        k = (fresnel(np.sum(half * v, -1), n) * slope_density(gradient, covariance)
             / (4 * v[..., 2] * half[..., 2] ** 4)
             / (1 + smith_lambda(v, covariance) + smith_lambda(s, covariance)))
    return np.where(up, k, 0.0)


def disk_glint(centre, radius_rad, view, covariance, n, points=37):
    """Glitter factor averaged over a uniform disk source, from a sunflower pattern of points of equal area."""
    centre = np.asarray(centre, float) / np.linalg.norm(centre)
    a = np.cross(centre, [0.0, 0.0, 1.0] if abs(centre[2]) < 0.9 else [1.0, 0.0, 0.0])
    a /= np.linalg.norm(a)
    b = np.cross(centre, a)
    k = np.arange(points)
    r = radius_rad * np.sqrt((k + 0.5) / points)
    phi = k * np.pi * (3 - np.sqrt(5))
    directions = (np.cos(r)[:, None] * centre
                  + np.sin(r)[:, None] * (np.cos(phi)[:, None] * a + np.sin(phi)[:, None] * b))
    return np.mean([glint(d, view, covariance, n) for d in directions], axis=0)


def facet_nodes(covariance, order=40):
    """Gauss-Hermite nodes and weights for gradients distributed with the given covariance."""
    x, w = np.polynomial.hermite.hermgauss(order)
    values, vectors = np.linalg.eigh(covariance)
    xx, yy = np.meshgrid(x, x, indexing="ij")
    unit = np.stack([xx.ravel(), yy.ravel()], -1) * np.sqrt(2 * np.maximum(values, 0.0))
    return unit @ vectors.T, np.outer(w, w).ravel() / np.pi


def sky_reflection(view, covariance, n, sky, order=40):
    """Radiance reflected toward the observer from a sky given as sky(directions) -> (N, channels).

    The visible facets are weighted by their projected area toward the observer and Smith's masking, which
    makes the weights sum to one. Rays that escape reach the sky; rays sent into the water or onto another
    wave reflect once more off a level surface.
    """
    v = np.asarray(view, float)
    gradient, weight = facet_nodes(covariance, order)
    normal = np.c_[-gradient, np.ones(len(gradient))]
    normal /= np.linalg.norm(normal, axis=1, keepdims=True)
    cos_v = normal @ v
    lam_v = smith_lambda(v, covariance)
    w = weight * np.maximum(cos_v, 0.0) / (v[2] * normal[:, 2]) / (1 + lam_v)
    reflected = 2 * cos_v[:, None] * normal - v
    rho = fresnel(cos_v, n)
    lam_i = smith_lambda(reflected, covariance)
    escape = np.where(reflected[:, 2] > 0, (1 + lam_v) / (1 + lam_v + lam_i), 0.0)
    level = reflected * [1.0, 1.0, -1.0]
    second = np.where((reflected[:, 2] < 0)[:, None], level, reflected)
    second[:, 2] = np.maximum(second[:, 2], 1e-6)
    second /= np.linalg.norm(second, axis=1, keepdims=True)
    direct = sky(np.where((reflected[:, 2] > 0)[:, None], reflected, second))
    again = sky(second) * fresnel(second[:, 2], n)[:, None]
    radiance = (w * rho * escape) @ direct + (w * rho * (1 - escape)) @ again
    return radiance, dict(weight_sum=float(w.sum()), escaping=float((w * escape).sum()))


def sky_reflection_batch(views, covariance, n, sky, order=24, block=2000):
    """sky_reflection for many view directions at once (views shape (V, 3)).

    Returns the mean reflected radiance (V, channels) and the mean of its square, so the facet-to-facet
    spread of radiance that shows the waves follows as their difference.
    """
    views = np.atleast_2d(np.asarray(views, float))
    gradient, weight = facet_nodes(covariance, order)
    normal = np.c_[-gradient, np.ones(len(gradient))]
    normal /= np.linalg.norm(normal, axis=1, keepdims=True)
    mean, square = [], []
    for start in range(0, len(views), block):
        v = views[start:start + block]
        cos_v = v @ normal.T                                           # (V, G)
        lam_v = smith_lambda(v, covariance)[:, None]
        w = weight[None, :] * np.maximum(cos_v, 0.0) / (v[:, 2:3] * normal[None, :, 2]) / (1 + lam_v)
        reflected = 2 * cos_v[..., None] * normal[None] - v[:, None, :]  # (V, G, 3)
        rho = fresnel(cos_v, n)
        lam_i = smith_lambda(reflected, covariance)
        up = reflected[..., 2] > 0
        escape = np.where(up, (1 + lam_v) / (1 + lam_v + lam_i), 0.0)
        second = np.where(up[..., None], reflected, reflected * [1.0, 1.0, -1.0])
        second[..., 2] = np.maximum(second[..., 2], 1e-6)
        second /= np.linalg.norm(second, axis=-1, keepdims=True)
        flat = second.reshape(-1, 3)
        seen = sky(flat).reshape(*second.shape[:2], -1)                 # (V, G, channels)
        again = fresnel(second[..., 2], n)
        value = rho[..., None] * seen * (escape + (1 - escape) * again)[..., None]
        mean.append(np.einsum("vg,vgc->vc", w, value))
        square.append(np.einsum("vg,vgc->vc", w, value ** 2))
    return np.concatenate(mean), np.concatenate(square)
