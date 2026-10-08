"""The solved sky seen from a height, and the air between an observer and what it sees.

The solved spherical sky (solved_transport.py) gives the radiance reaching a ground observer from the
scattering sources its successive orders computed. Here the same formal integration runs from an observer at
a small height above the ground, and stops at chosen distances: the light scattered into the line of sight
in front of a surface (the aerial perspective, or path radiance) and the transmittance to it. At the end of a
ray that escapes the air the in-scattered light is the sky's radiance; a ray that meets the ground ends there,
and the ground's own light is the caller's to add.

The sources are those of solved_transport.radiance: the stored diffuse moments (all orders but the last, which
only measured convergence) and the direct beam, at each path segment's mid radius and its local source
elevation. A distance inside a segment takes that segment's mid-point source over the partial length, with
the segment's exact homogeneous-layer weights.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

from illumination.sky.solver import bracket, ray_path
from illumination.sky.solved_transport import sample_at


@njit(cache=True)
def _segment_sources(path, data, beam, rg, ag, sr, sa, ss, ca, mu, cp, nl):
    """Scattering source (per unit volume scattering) at each segment's mid-point, per channel."""
    nu = ss * mu + ca * math.sqrt(max(0., 1 - mu * mu)) * cp
    out = np.zeros((len(path), nl))
    for k in range(len(path)):
        r, qr = path[k, 0], path[k, 1]
        ms = max(-1., min(1., ss * path[k, 2] + ca * path[k, 3] * cp))
        aa = math.asin(ms)
        qt = (nu - ms * qr) / math.sqrt(max(1e-16, 1 - ms * ms))
        qp2 = max(0., 1 - qr * qr - qt * qt)
        ri, u = bracket(rg, r)
        ai, v = bracket(ag, aa)
        bi, bu = bracket(sr, r)
        bj, bv = bracket(sa, aa)
        for l in range(nl):
            source = (sample_at(data, ri, ai, u, v, l, 0) + qr * qr * sample_at(data, ri, ai, u, v, l, 1) +
                      qt * qt * sample_at(data, ri, ai, u, v, l, 2) + qp2 * sample_at(data, ri, ai, u, v, l, 3) +
                      2 * qr * qt * sample_at(data, ri, ai, u, v, l, 4) +
                      (1 + nu * nu) * sample_at(beam, bi, bj, bu, bv, l, 0))
            out[k, l] = 3 / (16 * math.pi) * max(0., source)
    return out


@njit(parallel=True, cache=True)
def along_rays(data, beam, rg, ag, sr, sa, edges, sca, absorb, shells, suns, elevations, azimuths, height,
               distances, weights, xyz):
    """In-scattered XYZ along view rays up to each distance, and the transmittance per channel.

    suns: source elevations (rad); elevations, azimuths: view directions (rad, azimuth from the source);
    height: the observer's height above the lowest shell (m); distances (m, increasing; use a large value for
    the whole ray); weights: per-channel source strength; xyz: per-channel XYZ weights.
    Returns inscatter (suns, elevations, azimuths, distances, 3), transmittance (elevations, distances,
    channels) and the length of each ray through the air (elevations).
    """
    ns, ne, na, nd, nl = len(suns), len(elevations), len(azimuths), len(distances), sca.shape[1]
    inscatter = np.zeros((ns, ne, na, nd, 3))
    transmittance = np.ones((ne, nd, nl))
    length = np.zeros(ne)
    radius = edges[0] + height
    for iv in prange(ne):
        mu = math.sin(elevations[iv])
        path, ground, cg, sg = ray_path(radius, mu, edges[0], edges[-1], shells, 1., 1.)
        n = len(path)
        start = np.zeros(n)
        for k in range(1, n):
            start[k] = start[k - 1] + path[k - 1, 6]
        length[iv] = start[n - 1] + path[n - 1, 6]
        # Per segment: optical depth rate, scattering share and the transmittance at its start.
        chi = np.zeros((n, nl))
        share = np.zeros((n, nl))
        t0 = np.ones((n, nl))
        for k in range(n):
            layer = min(len(edges) - 2, max(0, np.searchsorted(edges, path[k, 0]) - 1))
            for l in range(nl):
                bs = sca[layer, l]
                chi[k, l] = bs + absorb[layer, l]
                share[k, l] = bs / chi[k, l] if chi[k, l] > 0 else 0.
                if k + 1 < n:
                    t0[k + 1, l] = t0[k, l] * math.exp(-chi[k, l] * path[k, 6])
        # Where each requested distance ends: the segment and the length into it.
        segment = np.zeros(nd, np.int64)
        into = np.zeros(nd)
        for m in range(nd):
            d = min(distances[m], length[iv])
            k = 0
            while k + 1 < n and start[k + 1] <= d:
                k += 1
            segment[m] = k
            into[m] = min(d - start[k], path[k, 6])
            for l in range(nl):
                transmittance[iv, m, l] = t0[k, l] * math.exp(-chi[k, l] * into[m])
        for i in range(ns):
            ca, ss = math.cos(suns[i]), math.sin(suns[i])
            for j in range(na):
                src = _segment_sources(path, data, beam, rg, ag, sr, sa, ss, ca, mu, math.cos(azimuths[j]), nl)
                total = np.zeros(3)
                m = 0
                for k in range(n):
                    # Requested distances that end inside this segment take a partial segment.
                    while m < nd and segment[m] == k:
                        partial = total.copy()
                        for l in range(nl):
                            dt = chi[k, l] * into[m]
                            value = t0[k, l] * (-math.expm1(-dt)) * share[k, l] * src[k, l] * weights[l]
                            partial[0] += value * xyz[l, 0]
                            partial[1] += value * xyz[l, 1]
                            partial[2] += value * xyz[l, 2]
                        inscatter[i, iv, j, m] = partial
                        m += 1
                    for l in range(nl):
                        dt = chi[k, l] * path[k, 6]
                        value = t0[k, l] * (-math.expm1(-dt)) * share[k, l] * src[k, l] * weights[l]
                        total[0] += value * xyz[l, 0]
                        total[1] += value * xyz[l, 1]
                        total[2] += value * xyz[l, 2]
    return inscatter, transmittance, length


def sky_arguments(solution):
    """The solved solution's arrays in along_rays' order (moments without the last order, as evaluate uses)."""
    return (solution["moments"] - solution["last_order"], solution["beam"], solution["r"], solution["a"],
            solution["sr"], solution["sa"], solution["edges"], solution["scattering"], solution["absorption"],
            solution["shells"])


def view(solution, source_elevation_deg, elevations_deg, azimuths_deg, distances_m, height_m=0.0, weights=None):
    """In-scattered XYZ (elevations, azimuths, distances, 3) and transmittance (elevations, distances, channels)
    for one source; azimuths are measured from the source's."""
    nl = solution["scattering"].shape[1]
    w = np.ones(nl) if weights is None else np.asarray(weights, float)
    inscatter, transmittance, length = along_rays(
        *sky_arguments(solution), np.radians([float(source_elevation_deg)]), np.radians(elevations_deg),
        np.radians(azimuths_deg), float(height_m), np.asarray(distances_m, float), w,
        np.ascontiguousarray(solution["xyz"], dtype=np.float64))
    return inscatter[0], transmittance, length
