"""Numba kernels of the sea-scene path tracer.

Frame: metres east (x), north (y) and up (z) from the sea-level point below the camera, the Moon's centre at
(0, 0, -R). The sea and the terrain are heightfields over the tangent plane lowered by the parabolic drop
(x^2 + y^2) / 2R, within 0.3 m of the sphere at 75 km. Each sea or terrain point is shaded in its own local
frame, whose vertical points away from the Moon's centre, so that grazing views near the horizon keep their
true angles.

The sea: the largest wave tile is explicit geometry, traced with a two-level grid of maximum heights; every
tile's slopes are filtered over the pixel's footprint by their mipmapped first and second moments (LEAN
mapping, Olano and Baker 2010), and the slopes left inside the footprint, with the waves no tile resolves,
form a Gaussian distribution of facets with that mean and covariance (Bruneton, Neyret and Holzschuch 2010).
Facets are drawn from the distribution of visible normals (Heitz and d'Eon 2014; Beckmann slopes stretched to
the covariance and sheared to the mean). Reflection is Fresnel's for unpolarized light; Smith's height-
correlated shadowing sends reflected rays that meet another wave, or point into the water, to one more
reflection off a level surface, as illumination/water_surface/reflection.py does. The Sun's glitter is the
Cox-Munk formula over the facet distribution widened by the disk; the Earth's glitter combines sampling its
disk and sampling the facets by multiple importance (Veach 1997), so its phase shows in calm water.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

SQRT_PI = math.sqrt(math.pi)
INV_SQRT_PI = 1.0 / SQRT_PI


# ---------------------------------------------------------------------------------------------- small helpers
@njit(cache=True, inline="always")
def dot3(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


@njit(cache=True, inline="always")
def norm3(a):
    n = math.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])
    return np.array([a[0] / n, a[1] / n, a[2] / n])


@njit(cache=True)
def erfinv(x):
    """Giles (2010) single-precision approximation refined by two Newton steps on math.erf."""
    x = min(max(x, -1 + 1e-15), 1 - 1e-15)
    w = -math.log((1.0 - x) * (1.0 + x))
    if w < 5.0:
        w -= 2.5
        p = 2.81022636e-08
        p = 3.43273939e-07 + p * w
        p = -3.5233877e-06 + p * w
        p = -4.39150654e-06 + p * w
        p = 0.00021858087 + p * w
        p = -0.00125372503 + p * w
        p = -0.00417768164 + p * w
        p = 0.246640727 + p * w
        p = 1.50140941 + p * w
    else:
        w = math.sqrt(w) - 3.0
        p = -0.000200214257
        p = 0.000100950558 + p * w
        p = 0.00134934322 + p * w
        p = -0.00367342844 + p * w
        p = 0.00573950773 + p * w
        p = -0.0076224613 + p * w
        p = 0.00943887047 + p * w
        p = 1.00167406 + p * w
        p = 2.83297682 + p * w
    y = p * x
    for _ in range(2):
        y -= (math.erf(y) - x) / (2 * INV_SQRT_PI * math.exp(-y * y))
    return y


@njit(cache=True)
def fresnel(c, n):
    c = min(max(c, 0.0), 1.0)
    root = math.sqrt(max(n * n - 1 + c * c, 0.0))
    rs = ((c - root) / (c + root)) ** 2
    rp = ((n * n * c - root) / (n * n * c + root)) ** 2
    return 0.5 * (rs + rp)


@njit(cache=True)
def local_frame(p, R):
    """Rotation taking camera-frame vectors to the local frame at p (vertical away from the Moon's centre)."""
    up = norm3(np.array([p[0], p[1], p[2] + R]))
    # Rodrigues rotation taking up to +z.
    ax, ay = up[1], -up[0]                 # up x z
    s = math.sqrt(ax * ax + ay * ay)
    c = up[2]
    m = np.eye(3)
    if s > 1e-15:
        kx, ky = ax / s, ay / s
        k = np.array([[0.0, 0.0, ky], [0.0, 0.0, -kx], [-ky, kx, 0.0]])
        kk = np.array([[kx * kx, kx * ky, 0.0], [kx * ky, ky * ky, 0.0], [0.0, 0.0, 0.0]])
        m = c * np.eye(3) + s * k + (1 - c) * kk
    return m


@njit(cache=True, inline="always")
def apply(m, v):
    return np.array([m[0, 0] * v[0] + m[0, 1] * v[1] + m[0, 2] * v[2],
                     m[1, 0] * v[0] + m[1, 1] * v[1] + m[1, 2] * v[2],
                     m[2, 0] * v[0] + m[2, 1] * v[1] + m[2, 2] * v[2]])


@njit(cache=True, inline="always")
def apply_t(m, v):
    return np.array([m[0, 0] * v[0] + m[1, 0] * v[1] + m[2, 0] * v[2],
                     m[0, 1] * v[0] + m[1, 1] * v[1] + m[2, 1] * v[2],
                     m[0, 2] * v[0] + m[1, 2] * v[1] + m[2, 2] * v[2]])


# ------------------------------------------------------------------------------------- heightfield tracing
@njit(cache=True, inline="always")
def ray_height(o, d, t, inv2r):
    x = o[0] + t * d[0]
    y = o[1] + t * d[1]
    return o[2] + t * d[2] + (x * x + y * y) * inv2r


@njit(cache=True, inline="always")
def min_ray_height(o, d, t0, t1, inv2r):
    """Minimum over [t0, t1] of the ray's height above the curved datum (a convex quadratic in t)."""
    a = (d[0] * d[0] + d[1] * d[1]) * inv2r
    b = d[2] + 2 * (o[0] * d[0] + o[1] * d[1]) * inv2r
    h0, h1 = ray_height(o, d, t0, inv2r), ray_height(o, d, t1, inv2r)
    low = min(h0, h1)
    if a > 0:
        tv = -b / (2 * a)
        if t0 < tv < t1:
            low = min(low, ray_height(o, d, tv, inv2r))
    return low


@njit(cache=True, inline="always")
def node(h, i, j, periodic):
    n0, n1 = h.shape[0], h.shape[1]
    if periodic:
        return h[j % n0, i % n1]
    return h[min(max(j, 0), n0 - 1), min(max(i, 0), n1 - 1)]


@njit(cache=True)
def surface_height(h, x, y, x0, y0, dx, periodic):
    """Bilinear height at (x, y) from nodes at (x0 + i dx, y0 + j dx) (rows j northward, columns i eastward)."""
    u = (x - x0) / dx
    v = (y - y0) / dx
    i = int(math.floor(u))
    j = int(math.floor(v))
    fu, fv = u - i, v - j
    return ((1 - fu) * (1 - fv) * node(h, i, j, periodic) + fu * (1 - fv) * node(h, i + 1, j, periodic) +
            (1 - fu) * fv * node(h, i, j + 1, periodic) + fu * fv * node(h, i + 1, j + 1, periodic))


@njit(cache=True)
def gap(o, d, t, h, x0, y0, dx, periodic, offset, inv2r):
    x = o[0] + t * d[0]
    y = o[1] + t * d[1]
    return ray_height(o, d, t, inv2r) - offset - surface_height(h, x, y, x0, y0, dx, periodic)


@njit(cache=True)
def refine(o, d, ta, tb, h, x0, y0, dx, periodic, offset, inv2r):
    """Bisection for the crossing between ta (above) and tb (below)."""
    for _ in range(40):
        tm = 0.5 * (ta + tb)
        if gap(o, d, tm, h, x0, y0, dx, periodic, offset, inv2r) > 0:
            ta = tm
        else:
            tb = tm
        if tb - ta < 1e-4 * max(1.0, ta) * 1e-3:
            break
    return 0.5 * (ta + tb)


@njit(cache=True)
def cell_hit(o, d, t0, t1, h, i, j, x0, y0, dx, periodic, offset, inv2r):
    """First crossing inside one cell's stretch of the ray, or -1."""
    top = max(node(h, i, j, periodic), node(h, i + 1, j, periodic), node(h, i, j + 1, periodic),
              node(h, i + 1, j + 1, periodic)) + offset
    if min_ray_height(o, d, t0, t1, inv2r) > top:
        return -1.0
    steps = 4
    ta = t0
    ga = gap(o, d, ta, h, x0, y0, dx, periodic, offset, inv2r)
    if ga <= 0:
        return t0
    for s in range(1, steps + 1):
        tb = t0 + (t1 - t0) * s / steps
        gb = gap(o, d, tb, h, x0, y0, dx, periodic, offset, inv2r)
        if gb <= 0:
            return refine(o, d, ta, tb, h, x0, y0, dx, periodic, offset, inv2r)
        ta, ga = tb, gb
    return -1.0


@njit(cache=True)
def trace_heightfield(o, d, t_start, t_end, h, blocks, block, x0, y0, dx, periodic, offset, inv2r):
    """First crossing of the ray with the heightfield between t_start and t_end, or -1.

    blocks holds the maximum node height of each block of block x block cells (periodic or not, as h).
    """
    if t_end <= t_start:
        return -1.0
    ncells_x = h.shape[1] if periodic else h.shape[1] - 1
    ncells_y = h.shape[0] if periodic else h.shape[0] - 1
    bw = block * dx
    nbx, nby = blocks.shape[1], blocks.shape[0]
    t = t_start
    # Block-level DDA.
    px, py = o[0] + t * d[0], o[1] + t * d[1]
    bi = int(math.floor((px - x0) / bw))
    bj = int(math.floor((py - y0) / bw))
    sx = 1 if d[0] > 0 else -1
    sy = 1 if d[1] > 0 else -1
    big = 1e30
    tdx = bw / abs(d[0]) if abs(d[0]) > 1e-15 else big
    tdy = bw / abs(d[1]) if abs(d[1]) > 1e-15 else big
    tx = (x0 + (bi + (1 if sx > 0 else 0)) * bw - o[0]) / d[0] if abs(d[0]) > 1e-15 else big
    ty = (y0 + (bj + (1 if sy > 0 else 0)) * bw - o[1]) / d[1] if abs(d[1]) > 1e-15 else big
    for _ in range(1000000):
        if t >= t_end:
            return -1.0
        t_next = min(tx, ty, t_end)
        inside = periodic or (0 <= bi < nbx and 0 <= bj < nby)
        if inside:
            top = blocks[bj % nby, bi % nbx] + offset
            if min_ray_height(o, d, t, t_next, inv2r) <= top:
                # Cell-level DDA within this block's stretch.
                ta = t
                qx, qy = o[0] + ta * d[0], o[1] + ta * d[1]
                ci = int(math.floor((qx - x0) / dx))
                cj = int(math.floor((qy - y0) / dx))
                cdx = dx / abs(d[0]) if abs(d[0]) > 1e-15 else big
                cdy = dx / abs(d[1]) if abs(d[1]) > 1e-15 else big
                cx = (x0 + (ci + (1 if sx > 0 else 0)) * dx - o[0]) / d[0] if abs(d[0]) > 1e-15 else big
                cy = (y0 + (cj + (1 if sy > 0 else 0)) * dx - o[1]) / d[1] if abs(d[1]) > 1e-15 else big
                for _c in range(4 * block + 8):
                    tb = min(cx, cy, t_next)
                    valid = periodic or (0 <= ci < ncells_x and 0 <= cj < ncells_y)
                    if valid and tb > ta:
                        hit = cell_hit(o, d, ta, tb, h, ci, cj, x0, y0, dx, periodic, offset, inv2r)
                        if hit >= 0:
                            return hit
                    if tb >= t_next:
                        break
                    ta = tb
                    if cx < cy:
                        ci += sx
                        cx += cdx
                    else:
                        cj += sy
                        cy += cdy
        elif not periodic:
            # Leaving the grid for good?
            if (bi < 0 and sx < 0) or (bi >= nbx and sx > 0) or (bj < 0 and sy < 0) or (bj >= nby and sy > 0):
                return -1.0
        t = t_next
        if tx < ty:
            bi += sx
            tx += tdx
        else:
            bj += sy
            ty += tdy
    return -1.0


@njit(cache=True)
def slab(o, d, level, inv2r):
    """Ray parameters where the ray's height above the curved datum crosses a level: (first, second) or (-1,-1)."""
    a = (d[0] * d[0] + d[1] * d[1]) * inv2r
    b = d[2] + 2 * (o[0] * d[0] + o[1] * d[1]) * inv2r
    c = o[2] + (o[0] * o[0] + o[1] * o[1]) * inv2r - level
    if a < 1e-300:
        if abs(b) < 1e-300:
            return -1.0, -1.0
        t = -c / b
        return (t, 1e30) if b < 0 else (-1e30, t)
    disc = b * b - 4 * a * c
    if disc < 0:
        return -1.0, -1.0
    r = math.sqrt(disc)
    return (-b - r) / (2 * a), (-b + r) / (2 * a)


@njit(cache=True)
def trace_sea(o, d, sea_h, sea_blocks, sea_block, sea_dx, hmin, hmax, inv2r, t_limit):
    """First hit of the explicit wave surface (periodic tile over the curved datum), or -1."""
    t0, t1 = slab(o, d, hmax, inv2r)
    if t1 < 0 and t0 < 0:
        return -1.0
    c = o[2] + (o[0] * o[0] + o[1] * o[1]) * inv2r
    t_in = 0.0 if c <= hmax else t0
    if t_in < 0:
        return -1.0
    u0, u1 = slab(o, d, hmin, inv2r)
    t_out = t1
    if u0 >= 0 and u0 > t_in:
        t_out = min(t_out, u0)          # below every trough: a crossing came before
    t_out = min(t_out, t_limit)
    if t_out <= t_in:
        return -1.0
    return trace_heightfield(o, d, t_in, t_out + 1e-6, sea_h, sea_blocks, sea_block, 0.0, 0.0, sea_dx, True,
                             0.0, inv2r)


# ------------------------------------------------------------------------------------------ slope statistics
@njit(cache=True)
def lean_fetch(data, offsets, sizes, spacing, levels, x, y, footprint, out):
    """Add one tile's filtered slope moments at (x, y) over a footprint (m) to out[0:5]."""
    lam = math.log2(max(footprint / spacing, 1e-9))
    top = levels - 1
    if lam <= 0:
        l0, f = 0, 0.0
    elif lam >= top:
        l0, f = top, 0.0
    else:
        l0 = int(math.floor(lam))
        f = lam - l0
    for which in range(2):
        level = l0 + which
        weight = (1 - f) if which == 0 else f
        if weight <= 0 or level > top:
            continue
        n = sizes[level]
        scale = 2.0 ** level
        u = x / (spacing * scale) - (scale - 1) / (2 * scale)
        v = y / (spacing * scale) - (scale - 1) / (2 * scale)
        i = int(math.floor(u))
        j = int(math.floor(v))
        fu, fv = u - i, v - j
        base = offsets[level]
        for k in range(5):
            plane = base + k * n * n
            a = data[plane + (j % n) * n + (i % n)]
            b = data[plane + (j % n) * n + ((i + 1) % n)]
            c = data[plane + ((j + 1) % n) * n + (i % n)]
            e = data[plane + ((j + 1) % n) * n + ((i + 1) % n)]
            out[k] += weight * ((1 - fu) * (1 - fv) * a + fu * (1 - fv) * b + (1 - fu) * fv * c + fu * fv * e)


@njit(cache=True)
def slope_statistics(lean_data, lean_offsets, lean_sizes, lean_spacing, lean_levels, lean_tiles, x, y, along,
                     minor, major, residual):
    """Mean slope (2) and covariance (3: xx, yy, xy) over an anisotropic footprint, summed over the tiles."""
    taps = int(min(8, max(1, math.ceil(major / max(minor, 1e-9)))))
    size = max(minor, major / taps)
    mean = np.zeros(2)
    cov = np.array([residual[0], residual[1], residual[2]])
    for tile in range(lean_tiles):
        acc = np.zeros(5)
        for k in range(taps):
            s = (k + 0.5) / taps - 0.5
            lean_fetch(lean_data, lean_offsets[tile], lean_sizes[tile], lean_spacing[tile], lean_levels[tile],
                       x + s * major * along[0], y + s * major * along[1], size, acc)
        for k in range(5):
            acc[k] /= taps
        mean[0] += acc[0]
        mean[1] += acc[1]
        cov[0] += max(acc[2] - acc[0] * acc[0], 0.0)
        cov[1] += max(acc[3] - acc[1] * acc[1], 0.0)
        cov[2] += acc[4] - acc[0] * acc[1]
    return mean, cov


@njit(cache=True)
def smith_lambda(w, mean, cov):
    """Smith's Lambda for a ray leaving a Gaussian-slope surface of the given mean and covariance (local frame)."""
    hx, hy = w[0], w[1]
    hor = math.sqrt(hx * hx + hy * hy)
    if hor < 1e-12:
        return 0.0
    ux, uy = hx / hor, hy / hor
    zs = w[2] - mean[0] * hx - mean[1] * hy           # sheared to the mean plane
    if zs <= 0:
        return 1e30
    sigma = math.sqrt(max(ux * ux * cov[0] + uy * uy * cov[1] + 2 * ux * uy * cov[2], 1e-20))
    a = zs / hor / (math.sqrt(2.0) * sigma)
    if a > 25:
        return 0.0
    return (math.exp(-a * a) - a * SQRT_PI * math.erfc(a)) / (2 * a * SQRT_PI)


@njit(cache=True)
def slope_pdf(sx, sy, mean, cov):
    dx_, dy_ = sx - mean[0], sy - mean[1]
    det = cov[0] * cov[1] - cov[2] * cov[2]
    if det <= 1e-30:
        det = 1e-30
    q = (cov[1] * dx_ * dx_ - 2 * cov[2] * dx_ * dy_ + cov[0] * dy_ * dy_) / det
    return math.exp(-0.5 * q) / (2 * math.pi * math.sqrt(det))


@njit(cache=True)
def sample_visible_11(theta, u1, u2):
    """Visible slopes of the unit Beckmann distribution for an incidence angle theta in the x-z plane
    (Heitz and d'Eon 2014, with Jakob's numerical inversion as in Mitsuba)."""
    if theta < 1e-4:
        r = math.sqrt(-math.log(max(1.0 - u1, 1e-300)))
        phi = 2 * math.pi * u2
        return r * math.cos(phi), r * math.sin(phi)
    tan_t = math.tan(theta)
    cot_t = 1.0 / tan_t
    a, c = -1.0, math.erf(cot_t)
    sample_x = max(u1, 1e-6)
    fit = 1 + theta * (-0.876 + theta * (0.4265 - 0.0594 * theta))
    b = c - (1 + c) * (1 - sample_x) ** fit
    normalization = 1.0 / (1 + c + INV_SQRT_PI * tan_t * math.exp(-cot_t * cot_t))
    for _ in range(12):
        if not (a <= b <= c):
            b = 0.5 * (a + c)
        inv = erfinv(b)
        value = normalization * (1 + b + INV_SQRT_PI * tan_t * math.exp(-inv * inv)) - sample_x
        derivative = normalization * (1 - inv * tan_t)
        if abs(value) < 1e-7:
            break
        if value > 0:
            c = b
        else:
            a = b
        if derivative != 0:
            b -= value / derivative
    return erfinv(b), erfinv(2.0 * max(u2, 1e-6) - 1.0)


@njit(cache=True)
def sample_visible_slope(v, mean, cov, u1, u2):
    """A facet slope drawn from the visible-normal distribution of Gaussian slopes N(mean, cov) seen from v."""
    vz = v[2] - mean[0] * v[0] - mean[1] * v[1]
    # Principal axes of the covariance.
    tr = cov[0] + cov[1]
    df = cov[0] - cov[1]
    root = math.sqrt(max(0.25 * df * df + cov[2] * cov[2], 0.0))
    l1, l2 = max(0.5 * tr + root, 1e-14), max(0.5 * tr - root, 1e-14)
    angle = 0.5 * math.atan2(2 * cov[2], df)
    ca, sa = math.cos(angle), math.sin(angle)
    # View in the principal frame, stretched to unit roughness (Beckmann alpha^2 = 2 sigma^2).
    a1, a2 = math.sqrt(2 * l1), math.sqrt(2 * l2)
    px = ca * v[0] + sa * v[1]
    py = -sa * v[0] + ca * v[1]
    w = norm3(np.array([a1 * px, a2 * py, max(vz, 1e-12)]))
    theta = math.acos(min(1.0, w[2]))
    phi = math.atan2(w[1], w[0])
    sx, sy = sample_visible_11(theta, u1, u2)
    cp, sp = math.cos(phi), math.sin(phi)
    rx, ry = cp * sx - sp * sy, sp * sx + cp * sy
    rx, ry = a1 * rx, a2 * ry
    # Back to east-north, and to the mean plane. Mitsuba's slope is the height gradient (normal (-g, 1)).
    gx = (ca * rx - sa * ry) + mean[0]
    gy = (sa * rx + ca * ry) + mean[1]
    return gx, gy


# ---------------------------------------------------------------------------------------------- sky lookups
@njit(cache=True)
def interp_index(grid, x):
    n = len(grid)
    if x <= grid[0]:
        return 0, 0.0
    if x >= grid[n - 1]:
        return n - 2, 1.0
    lo, hi = 0, n - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if grid[mid] <= x:
            lo = mid
        else:
            hi = mid
    return lo, (x - grid[lo]) / (grid[lo + 1] - grid[lo])


@njit(cache=True)
def sky_xyz(sky, sky_el, az_step, d):
    """Sky radiance (XYZ) in direction d (local frame of the viewer), from a map over elevation and azimuth."""
    el = math.degrees(math.asin(min(1.0, max(-1.0, d[2]))))
    az = math.degrees(math.atan2(d[0], d[1])) % 360.0
    i, u = interp_index(sky_el, el)
    na = sky.shape[1]
    q = az / az_step
    j = int(math.floor(q)) % na
    v = q - math.floor(q)
    j1 = (j + 1) % na
    out = np.zeros(3)
    for k in range(3):
        out[k] = ((1 - u) * ((1 - v) * sky[i, j, k] + v * sky[i, j1, k]) +
                  u * ((1 - v) * sky[i + 1, j, k] + v * sky[i + 1, j1, k]))
    return out


@njit(cache=True)
def earth_disk(d, earth):
    """Radiance (XYZ) of the Earth's lit disk in direction d, or zeros.

    earth: [dir(3), angular radius, sun direction (3), radiance XYZ (3)].
    """
    out = np.zeros(3)
    c = d[0] * earth[0] + d[1] * earth[1] + d[2] * earth[2]
    rho = earth[3]
    if c < math.cos(rho):
        return out
    # The point of the Earth's sphere seen in direction d: its normal tilts from -centre by the offset.
    off = np.array([d[0] - c * earth[0], d[1] - c * earth[1], d[2] - c * earth[2]])
    s = math.sqrt(dot3(off, off))
    frac = min(1.0, math.sin(math.acos(min(1.0, c))) / math.sin(rho))
    normal = np.array([-earth[0], -earth[1], -earth[2]]) * math.sqrt(max(0.0, 1 - frac * frac))
    if s > 0:
        normal = normal + off / s * frac
    if normal[0] * earth[4] + normal[1] * earth[5] + normal[2] * earth[6] <= 0:
        return out
    out[0], out[1], out[2] = earth[7], earth[8], earth[9]
    return out


@njit(cache=True)
def aerial(lut_s, lut_t, lut_el, lut_az_step, lut_d, d, distance, radiance):
    """Radiance at the eye from a surface of radiance (XYZ) at a distance along d: T L + S."""
    el = math.degrees(math.asin(min(1.0, max(-1.0, d[2]))))
    az = math.degrees(math.atan2(d[0], d[1])) % 360.0
    i, u = interp_index(lut_el, el)
    k, w = interp_index(lut_d, distance)
    na = lut_s.shape[1]
    q = az / lut_az_step
    j = int(math.floor(q)) % na
    v = q - math.floor(q)
    j1 = (j + 1) % na
    out = np.zeros(3)
    for c in range(3):
        s = ((1 - u) * (1 - v) * ((1 - w) * lut_s[i, j, k, c] + w * lut_s[i, j, k + 1, c]) +
             (1 - u) * v * ((1 - w) * lut_s[i, j1, k, c] + w * lut_s[i, j1, k + 1, c]) +
             u * (1 - v) * ((1 - w) * lut_s[i + 1, j, k, c] + w * lut_s[i + 1, j, k + 1, c]) +
             u * v * ((1 - w) * lut_s[i + 1, j1, k, c] + w * lut_s[i + 1, j1, k + 1, c]))
        t = ((1 - u) * ((1 - w) * lut_t[i, k, c] + w * lut_t[i, k + 1, c]) +
             u * ((1 - w) * lut_t[i + 1, k, c] + w * lut_t[i + 1, k + 1, c]))
        out[c] = t * radiance[c] + s
    return out


@njit(cache=True)
def aerial_bands(lut_s, lut_tb, lut_el, lut_az_step, lut_d, d, distance, bands):
    """As aerial, for a surface radiance given as XYZ per spectral band (B, 3), each band transmitted alone."""
    el = math.degrees(math.asin(min(1.0, max(-1.0, d[2]))))
    az = math.degrees(math.atan2(d[0], d[1])) % 360.0
    i, u = interp_index(lut_el, el)
    k, w = interp_index(lut_d, distance)
    na = lut_s.shape[1]
    q = az / lut_az_step
    j = int(math.floor(q)) % na
    v = q - math.floor(q)
    j1 = (j + 1) % na
    out = np.zeros(3)
    for c in range(3):
        s = ((1 - u) * (1 - v) * ((1 - w) * lut_s[i, j, k, c] + w * lut_s[i, j, k + 1, c]) +
             (1 - u) * v * ((1 - w) * lut_s[i, j1, k, c] + w * lut_s[i, j1, k + 1, c]) +
             u * (1 - v) * ((1 - w) * lut_s[i + 1, j, k, c] + w * lut_s[i + 1, j, k + 1, c]) +
             u * v * ((1 - w) * lut_s[i + 1, j1, k, c] + w * lut_s[i + 1, j1, k + 1, c]))
        out[c] = s
    for band in range(bands.shape[0]):
        t = ((1 - u) * ((1 - w) * lut_tb[i, k, band] + w * lut_tb[i, k + 1, band]) +
             u * ((1 - w) * lut_tb[i + 1, k, band] + w * lut_tb[i + 1, k + 1, band]))
        for c in range(3):
            out[c] += t * bands[band, c]
    return out


@njit(cache=True)
def terrain_radiance(tex, x, y, x0, y0, dx):
    """Lambertian terrain radiance (B, 3) at (x, y), bilinear in the radiance texture."""
    u = (x - x0) / dx
    v = (y - y0) / dx
    ny, nx = tex.shape[0], tex.shape[1]
    i = min(max(int(math.floor(u)), 0), nx - 2)
    j = min(max(int(math.floor(v)), 0), ny - 2)
    fu = min(max(u - i, 0.0), 1.0)
    fv = min(max(v - j, 0.0), 1.0)
    return ((1 - fu) * (1 - fv) * tex[j, i] + fu * (1 - fv) * tex[j, i + 1] +
            (1 - fu) * fv * tex[j + 1, i] + fu * fv * tex[j + 1, i + 1])


# ------------------------------------------------------------------------------------------------- shading
# params (float64): see scene.PARAMS for the layout.
P_R, P_SEA_DX, P_SEA_HMIN, P_SEA_HMAX, P_N, P_F0 = 0, 1, 2, 3, 4, 5
P_TER_X0, P_TER_Y0, P_TER_DX, P_TER_OFFSET, P_TER_HMAX, P_HAS_TERRAIN = 6, 7, 8, 9, 10, 11
P_TEX_X0, P_TEX_Y0, P_TEX_DX = 12, 13, 14
P_SKY_AZ_STEP, P_LUT_AZ_STEP, P_PIXEL, P_INTERREFLECT, P_SEA_BLOCK, P_TER_BLOCK = 15, 16, 17, 18, 19, 20
P_TAN_HALF, P_ASPECT, P_SUN_TEST = 21, 22, 23
# sun: dir(3), radius, normal irradiance XYZ(3), disk radiance XYZ(3)
# earth: dir(3), radius, sun dir(3), radiance XYZ(3), solid angle, up (0/1)


@njit(cache=True)
def terrain_hit(o, d, t_max, ter_h, ter_blocks, params):
    if params[P_HAS_TERRAIN] == 0:
        return -1.0
    inv2r = 0.5 / params[P_R]
    lo, hi = slab(o, d, params[P_TER_HMAX], inv2r)
    t_end = t_max
    if hi > 0:
        t_end = min(t_end, hi)              # beyond this the ray stays above every summit
    elif lo < 0:
        return -1.0
    return trace_heightfield(o, d, 0.0, t_end, ter_h, ter_blocks, int(params[P_TER_BLOCK]), params[P_TER_X0],
                             params[P_TER_Y0], params[P_TER_DX], False, params[P_TER_OFFSET], inv2r)


@njit(cache=True)
def source_blocked(p, s, ter_h, ter_blocks, params):
    """Whether terrain hides a source in direction s (camera frame) from point p."""
    if params[P_HAS_TERRAIN] == 0 or s[2] > params[P_SUN_TEST]:
        return False
    q = np.array([p[0] + 0.05 * s[0], p[1] + 0.05 * s[1], p[2] + 0.05])
    return terrain_hit(q, s, 2e5, ter_h, ter_blocks, params) >= 0


@njit(cache=True)
def sky_local(sky, sky_el, params, d_local, earth, include_earth):
    out = sky_xyz(sky, sky_el, params[P_SKY_AZ_STEP], d_local)
    if include_earth and earth[11] > 0:
        e = earth_disk(d_local, earth)
        out += e
    return out


@njit(cache=True)
def shade_sea(p, d_in, minor, major, depth, u1, u2, u3, u4, params, sea_h, sea_blocks, lean_data, lean_offsets,
              lean_sizes, lean_spacing, lean_levels, residual, ter_h, ter_blocks, tex, sky, sky_el, lut_s, lut_tb,
              lut_el, lut_d, sun, earth, water):
    """Radiance (XYZ) leaving the sea at p toward the incoming ray's origin."""
    R = params[P_R]
    inv2r = 0.5 / R
    n = params[P_N]
    m_loc = local_frame(p, R)
    v = apply(m_loc, np.array([-d_in[0], -d_in[1], -d_in[2]]))
    if v[2] < 1e-6:
        v[2] = 1e-6
        v = norm3(v)
    hor = math.sqrt(d_in[0] * d_in[0] + d_in[1] * d_in[1])
    along = np.array([d_in[0] / hor, d_in[1] / hor]) if hor > 1e-12 else np.array([1.0, 0.0])
    mean, cov = slope_statistics(lean_data, lean_offsets, lean_sizes, lean_spacing, lean_levels,
                                 lean_offsets.shape[0], p[0], p[1], along, minor, major, residual)
    gx, gy = sample_visible_slope(v, mean, cov, u1, u2)
    m = norm3(np.array([-gx, -gy, 1.0]))
    cvm = max(dot3(v, m), 1e-9)
    F = fresnel(cvm, n)
    r = np.array([2 * cvm * m[0] - v[0], 2 * cvm * m[1] - v[1], 2 * cvm * m[2] - v[2]])
    lam_v = smith_lambda(v, mean, cov)
    above = r[2] - mean[0] * r[0] - mean[1] * r[1] > 0
    escape = 0.0
    if above:
        lam_r = smith_lambda(r, mean, cov)
        escape = (1 + lam_v) / (1 + lam_v + lam_r)
    out = np.zeros(3)
    earth_up = earth[11] > 0
    omega_e = earth[10]
    g1v = 1.0 / (1.0 + lam_v)
    # The Earth as seen in this point's local frame, for lookups with local directions.
    earth_l = earth.copy()
    e_loc = apply(m_loc, np.array([earth[0], earth[1], earth[2]]))
    s_loc = apply(m_loc, np.array([earth[4], earth[5], earth[6]]))
    earth_l[0], earth_l[1], earth_l[2] = e_loc[0], e_loc[1], e_loc[2]
    earth_l[4], earth_l[5], earth_l[6] = s_loc[0], s_loc[1], s_loc[2]
    # The reflected ray that escapes: other waves, the coast, or the sky (the Earth weighted against its sampling).
    if escape > 0:
        r_cam = apply_t(m_loc, r)
        done = False
        if depth == 0 and minor < params[P_INTERREFLECT] and r_cam[2] < 0.2:
            q = np.array([p[0] + 0.01 * r_cam[0], p[1] + 0.01 * r_cam[1], p[2] + 0.01 * r_cam[2] + 0.002])
            t_sea = trace_sea(q, r_cam, sea_h, sea_blocks, int(params[P_SEA_BLOCK]), params[P_SEA_DX],
                              params[P_SEA_HMIN], params[P_SEA_HMAX], inv2r, 3000.0)
            t_ter = terrain_hit(q, r_cam, t_sea if t_sea > 0 else 2e5, ter_h, ter_blocks, params)
            if t_ter >= 0:
                x = q + t_ter * r_cam
                bands = terrain_radiance(tex, x[0], x[1], params[P_TEX_X0], params[P_TEX_Y0], params[P_TEX_DX])
                out += F * escape * aerial_bands(lut_s, lut_tb, lut_el, params[P_LUT_AZ_STEP], lut_d, r_cam,
                                                 t_ter, bands)
                done = True
            elif t_sea >= 0:
                x = q + t_sea * r_cam
                second = shade_sea_level(x, r_cam, t_sea * params[P_PIXEL] + minor, params, lean_data, lean_offsets,
                                         lean_sizes, lean_spacing, lean_levels, residual, sky, sky_el, earth, water,
                                         u3, u4)
                out += F * escape * second
                done = True
        elif r_cam[2] < 0.2 and params[P_HAS_TERRAIN] > 0:
            t_ter = terrain_hit(np.array([p[0], p[1], p[2] + 0.02]), r_cam, 2e5, ter_h, ter_blocks, params)
            if t_ter >= 0:
                x = p + t_ter * r_cam
                bands = terrain_radiance(tex, x[0], x[1], params[P_TEX_X0], params[P_TEX_Y0], params[P_TEX_DX])
                out += F * escape * aerial_bands(lut_s, lut_tb, lut_el, params[P_LUT_AZ_STEP], lut_d, r_cam,
                                                 t_ter, bands)
                done = True
        if not done:
            out += F * escape * sky_xyz(sky, sky_el, params[P_SKY_AZ_STEP], r)
            if earth_up:
                e = earth_disk(r, earth_l)
                if e[1] > 0:
                    m_slope = slope_pdf(gx, gy, mean, cov) / m[2] ** 4
                    pb = g1v * m_slope / (4 * v[2])
                    pl = 1.0 / omega_e
                    out += F * escape * e * pb * pb / (pb * pb + pl * pl)
    # The rest reflects once more off a level surface.
    if escape < 1:
        r2 = norm3(np.array([r[0], r[1], max(abs(r[2]), 1e-6)]))
        out += F * (1 - escape) * fresnel(r2[2], n) * sky_local(sky, sky_el, params, r2, earth_l, True)
    # The Earth sampled over its disk.
    if earth_up:
        cos_max = math.cos(earth[3])
        z = 1 - u3 * (1 - cos_max)
        phi = 2 * math.pi * u4
        sz = math.sqrt(max(0.0, 1 - z * z))
        e_dir = np.array([earth[0], earth[1], earth[2]])
        a = norm3(np.array([-e_dir[1], e_dir[0], 0.0])) if abs(e_dir[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
        b = np.array([e_dir[1] * a[2] - e_dir[2] * a[1], e_dir[2] * a[0] - e_dir[0] * a[2],
                      e_dir[0] * a[1] - e_dir[1] * a[0]])
        l_cam = z * e_dir + sz * math.cos(phi) * a + sz * math.sin(phi) * b
        l = apply(m_loc, l_cam)
        if l[2] > 0:
            e = earth_disk(l_cam, earth)
            if e[1] > 0 and not source_blocked(p, l_cam, ter_h, ter_blocks, params):
                h = norm3(np.array([v[0] + l[0], v[1] + l[1], v[2] + l[2]]))
                if h[2] > 0:
                    d_h = slope_pdf(-h[0] / h[2], -h[1] / h[2], mean, cov) / h[2] ** 4
                    lam_l = smith_lambda(l, mean, cov)
                    g2 = 1.0 / (1 + lam_v + lam_l)
                    pl = 1.0 / omega_e
                    pb = g1v * d_h / (4 * v[2])
                    value = fresnel(dot3(v, h), n) * d_h * g2 / (4 * v[2] * pl)
                    out += e * value * pl * pl / (pl * pl + pb * pb)
    # The Sun's glitter: the Cox-Munk formula with the disk's spread added to the facets'.
    s_cam = np.array([sun[0], sun[1], sun[2]])
    s = apply(m_loc, s_cam)
    if s[2] > -sun[3] and sun[5] > 0 and not source_blocked(p, s_cam, ter_h, ter_blocks, params):
        h = norm3(np.array([v[0] + s[0], v[1] + s[1], v[2] + s[2]]))
        if h[2] > 0:
            spread = sun[3] * sun[3] / 16
            wide = np.array([cov[0] + spread, cov[1] + spread, cov[2]])
            pdf = slope_pdf(-h[0] / h[2], -h[1] / h[2], mean, wide)
            lam_s = smith_lambda(s, mean, wide) if s[2] > 0 else 0.0
            k = fresnel(dot3(v, h), n) * pdf / (4 * v[2] * h[2] ** 4) / (1 + smith_lambda(v, mean, wide) + lam_s)
            out[0] += k * sun[4]
            out[1] += k * sun[5]
            out[2] += k * sun[6]
    # The water's own light, through the facet.
    t = (1 - F) / (1 - params[P_F0])
    out[0] += t * water[0]
    out[1] += t * water[1]
    out[2] += t * water[2]
    return out


@njit(cache=True)
def shade_sea_level(p, d_in, footprint, params, lean_data, lean_offsets, lean_sizes, lean_spacing, lean_levels,
                    residual, sky, sky_el, earth, water, u1, u2):
    """A second, statistical look at the sea for rays reflected onto another wave: sky reflection only."""
    R = params[P_R]
    n = params[P_N]
    m_loc = local_frame(p, R)
    v = apply(m_loc, np.array([-d_in[0], -d_in[1], -d_in[2]]))
    if v[2] < 1e-6:
        v[2] = 1e-6
        v = norm3(v)
    mean, cov = slope_statistics(lean_data, lean_offsets, lean_sizes, lean_spacing, lean_levels,
                                 lean_offsets.shape[0], p[0], p[1], np.array([1.0, 0.0]), footprint, footprint,
                                 residual)
    gx, gy = sample_visible_slope(v, mean, cov, u1, u2)
    m = norm3(np.array([-gx, -gy, 1.0]))
    cvm = max(dot3(v, m), 1e-9)
    F = fresnel(cvm, n)
    r = np.array([2 * cvm * m[0] - v[0], 2 * cvm * m[1] - v[1], 2 * cvm * m[2] - v[2]])
    r = norm3(np.array([r[0], r[1], max(abs(r[2]), 1e-6)]))
    earth_l = earth.copy()
    e_loc = apply(m_loc, np.array([earth[0], earth[1], earth[2]]))
    s_loc = apply(m_loc, np.array([earth[4], earth[5], earth[6]]))
    earth_l[0], earth_l[1], earth_l[2] = e_loc[0], e_loc[1], e_loc[2]
    earth_l[4], earth_l[5], earth_l[6] = s_loc[0], s_loc[1], s_loc[2]
    out = F * sky_local(sky, sky_el, params, r, earth_l, True)
    t = (1 - F) / (1 - params[P_F0])
    out[0] += t * water[0]
    out[1] += t * water[1]
    out[2] += t * water[2]
    return out


@njit(cache=True)
def camera_direction(u, v, camera, params):
    th = params[P_TAN_HALF]
    a = params[P_ASPECT]
    d = np.array([camera[3] + (2 * u - 1) * th * camera[6] + (1 - 2 * v) * th / a * camera[9],
                  camera[4] + (2 * u - 1) * th * camera[7] + (1 - 2 * v) * th / a * camera[10],
                  camera[5] + (2 * u - 1) * th * camera[8] + (1 - 2 * v) * th / a * camera[11]])
    return norm3(d)


@njit(cache=True)
def trace_camera(o, d, rnd, params, sea_h, sea_blocks, lean_data, lean_offsets, lean_sizes, lean_spacing,
                 lean_levels, residual, ter_h, ter_blocks, tex, sky, sky_el, lut_s, lut_t, lut_tb, lut_el, lut_d, sun,
                 earth, water):
    """Radiance (XYZ) reaching the eye along d, and what it met (0 sky, 1 sea, 2 land)."""
    inv2r = 0.5 / params[P_R]
    t_sea = trace_sea(o, d, sea_h, sea_blocks, int(params[P_SEA_BLOCK]), params[P_SEA_DX], params[P_SEA_HMIN],
                      params[P_SEA_HMAX], inv2r, 1e6)
    t_ter = terrain_hit(o, d, t_sea if t_sea > 0 else 3e5, ter_h, ter_blocks, params)
    if t_ter >= 0:
        x = o + t_ter * d
        bands = terrain_radiance(tex, x[0], x[1], params[P_TEX_X0], params[P_TEX_Y0], params[P_TEX_DX])
        return aerial_bands(lut_s, lut_tb, lut_el, params[P_LUT_AZ_STEP], lut_d, d, t_ter, bands), 2
    if t_sea >= 0:
        p = o + t_sea * d
        m_loc = local_frame(p, params[P_R])
        vz = abs(apply(m_loc, d)[2])
        minor = t_sea * params[P_PIXEL]
        major = minor / max(vz, 1e-3)
        sea = shade_sea(p, d, minor, major, 0, rnd[0], rnd[1], rnd[2], rnd[3], params, sea_h, sea_blocks,
                        lean_data, lean_offsets, lean_sizes, lean_spacing, lean_levels, residual, ter_h, ter_blocks,
                        tex, sky, sky_el, lut_s, lut_tb, lut_el, lut_d, sun, earth, water)
        return aerial(lut_s, lut_t, lut_el, params[P_LUT_AZ_STEP], lut_d, d, t_sea, sea), 1
    out = sky_xyz(sky, sky_el, params[P_SKY_AZ_STEP], d)
    if d[0] * sun[0] + d[1] * sun[1] + d[2] * sun[2] >= math.cos(sun[3]):
        out[0] += sun[7]
        out[1] += sun[8]
        out[2] += sun[9]
    if earth[11] > 0:
        out += earth_disk(d, earth)
    return out, 0


@njit(parallel=True, cache=True)
def render(width, height, row0, row1, spp_side, seed, camera, params, sea_h, sea_blocks, lean_data, lean_offsets,
           lean_sizes, lean_spacing, lean_levels, residual, ter_h, ter_blocks, tex, sky, sky_el, lut_s, lut_t,
           lut_tb, lut_el, lut_d, sun, earth, water):
    """Mean XYZ radiance per pixel for rows row0..row1, the standard error of its luminance and what each pixel
    met (fraction of samples on sea and on land). spp_side^2 stratified samples per pixel."""
    rows = row1 - row0
    mean = np.zeros((rows, width, 3), np.float32)
    error = np.zeros((rows, width), np.float32)
    kind = np.zeros((rows, width, 2), np.float32)
    o = np.array([camera[0], camera[1], camera[2]])
    n = spp_side * spp_side
    for jj in prange(rows):
        j = row0 + jj
        for i in range(width):
            np.random.seed((seed * 1000003 + j * 7919 + i * 104729) % 2147483647)
            acc = np.zeros(3)
            acc2 = 0.0
            sea = 0.0
            land = 0.0
            for s in range(n):
                su = (s % spp_side + np.random.random()) / spp_side
                sv = (s // spp_side + np.random.random()) / spp_side
                d = camera_direction((i + su) / width, (j + sv) / height, camera, params)
                rnd = np.random.random(4)
                L, what = trace_camera(o, d, rnd, params, sea_h, sea_blocks, lean_data, lean_offsets, lean_sizes,
                                       lean_spacing, lean_levels, residual, ter_h, ter_blocks, tex, sky, sky_el, lut_s,
                                       lut_t, lut_tb, lut_el, lut_d, sun, earth, water)
                acc += L
                acc2 += L[1] * L[1]
                if what == 1:
                    sea += 1
                elif what == 2:
                    land += 1
            m = acc / n
            mean[jj, i, 0], mean[jj, i, 1], mean[jj, i, 2] = m[0], m[1], m[2]
            error[jj, i] = math.sqrt(max(acc2 / n - m[1] * m[1], 0.0) / max(n - 1, 1))
            kind[jj, i, 0] = sea / n
            kind[jj, i, 1] = land / n
    return mean, error, kind
