"""A sea surface drawn from its directional spectrum at any gravity, for rendering explicit waves.

The surface is a linear random sea: a sum of Fourier modes with independent Gaussian amplitudes whose
variance is the spectrum's, on periodic square tiles. Waves from hundreds of metres down to the capillary
roll-off do not fit one grid, so the wavenumbers are shared among tiles of decreasing size (cascades, as in
Tessendorf's 2001 ocean): each wavenumber belongs to exactly one tile, which samples it at least four times
per wavelength. The slopes of the waves shorter than the last tile's limit form a residual covariance, for a
reflection model to take statistically.

The waves the wave model resolves come from its directional spectrum (SWAN's variance density in m2/Hz/degree,
Cartesian propagation directions); the shorter ones from the short-wave regime of the unified spectrum
(short_waves.py), above the same cutoff as the sea-slope study uses. Both carry the gravity they were computed
at, so lunar waves keep their lunar lengths.

Gravity waves carry the horizontal displacement of first-order Lagrangian (Gerstner) motion, the choppy wave
model of Nouguier, Guérin and Chapron (2009, J. Geophys. Res. 114, C09012), which sharpens crests and widens
troughs while keeping the elevation spectrum to first order. Tiles whose waves are shorter than a few times
the capillary length stay linear, since capillary waves do not steepen that way. The displaced surface is
resampled onto a regular grid (Eulerian heights and slopes) by fixed-point inversion of the displacement.

Coordinates are metres east and north; slopes are the gradient of the height, east and north.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.ndimage import map_coordinates

from illumination.water_surface.model import wavenumber
from illumination.water_surface.short_waves import capillary_scales, phase_speed, short_wave_level

# Tile sizes and grids: 2048 m with 1 m cells down to 1.43 m with 1.4 mm cells. The sizes are mutually
# incommensurate, so the tiles' repeats do not line up.
CASCADES = ((2048.0, 2048), (131.0, 1024), (13.7, 1024), (1.43, 1024))
SAMPLES_PER_WAVELENGTH = 4.0


@dataclass
class Tile:
    """One periodic tile: its size, grid and the wavenumbers it holds."""

    size_m: float
    cells: int
    k_low: float
    k_high: float
    choppy: bool
    height: np.ndarray | None = None          # Eulerian heights (m), only for tiles kept as geometry
    slope: np.ndarray | None = None           # Eulerian slopes (2, cells, cells)
    expected_covariance: np.ndarray = field(default_factory=lambda: np.zeros((2, 2)))
    expected_variance: float = 0.0

    @property
    def spacing_m(self):
        return self.size_m / self.cells


def swan_density(frequency_hz, direction_deg, variance_m2_hz_deg, depth_m, gravity):
    """Wavenumber-vector density F(kx, ky) (m^4) of a SWAN directional spectrum, as a function of (kx, ky).

    F dkx dky = E df dtheta: F = E (180/pi) (df/dk) / k, with df/dk the group speed over 2 pi from
    omega^2 = g k tanh(k depth). E is interpolated linearly in frequency and periodically in direction,
    and is zero outside the printed frequencies.
    """
    f = np.asarray(frequency_hz, float)
    theta = np.asarray(direction_deg, float)
    e = np.asarray(variance_m2_hz_deg, float)
    order = np.argsort(theta % 360)
    theta, e = theta[order] % 360, e[:, order]
    k_grid = wavenumber(f, depth_m, gravity)
    step = 360.0 / len(theta)

    def density(kx, ky):
        k = np.hypot(kx, ky)
        out = np.zeros_like(k)
        inside = (k >= k_grid[0]) & (k <= k_grid[-1])
        kk = k[inside]
        # Frequency of each wavenumber and the group speed from the finite-depth dispersion relation.
        tanh = np.tanh(kk * depth_m)
        omega = np.sqrt(gravity * kk * tanh)
        group = 0.5 * omega / kk * (1 + 2 * kk * depth_m / np.sinh(np.minimum(2 * kk * depth_m, 700)))
        freq = omega / (2 * np.pi)
        i = np.clip(np.searchsorted(f, freq) - 1, 0, len(f) - 2)
        u = np.clip((freq - f[i]) / (f[i + 1] - f[i]), 0, 1)
        angle = np.degrees(np.arctan2(ky[inside], kx[inside])) % 360
        j = np.floor((angle - theta[0]) / step).astype(int) % len(theta)
        v = ((angle - theta[0]) / step) % 1.0
        j1 = (j + 1) % len(theta)
        value = ((1 - u) * ((1 - v) * e[i, j] + v * e[i, j1]) + u * ((1 - v) * e[i + 1, j] + v * e[i + 1, j1]))
        out[inside] = value * (180 / np.pi) * group / (2 * np.pi) / kk
        return out

    return density


def short_density(u_star, toward_deg, gravity, k_from, peak_speed=np.inf):
    """Wavenumber-vector density of the unified spectrum's short-wave regime above k_from (Elfouhaily 1997, 40, 49).

    F(k, phi) = k^-4 B_h(k) [1 + Delta(k) cos 2(phi - wind)] / 2 pi, with B_h the short waves' curvature
    spectrum; its slope integral is short_waves.short_wave_slopes for the same arguments.
    """
    k_m, c_m = capillary_scales(gravity)
    level = float(short_wave_level(u_star, gravity))
    wind = np.radians(toward_deg)

    def density(kx, ky):
        k = np.hypot(kx, ky)
        out = np.zeros_like(k)
        inside = k >= k_from
        if level <= 0 or not inside.any():
            return out
        kk = k[inside]
        c = phase_speed(kk, gravity)
        b = 0.5 * level * c_m / c * np.exp(-0.25 * (kk / k_m - 1) ** 2)
        delta = np.tanh(np.log(2) / 4 + 4 * (c / peak_speed) ** 2.5 + 0.13 * u_star / c_m * (c_m / c) ** 2.5)
        phi = np.arctan2(ky[inside], kx[inside])
        out[inside] = kk ** -4 * b * (1 + delta * np.cos(2 * (phi - wind))) / (2 * np.pi)
        return out

    return density


def residual_covariance(u_star, toward_deg, gravity, k_from, peak_speed=np.inf, points=4000):
    """East-north slope covariance of the short waves above k_from, which no tile resolves."""
    k_m, c_m = capillary_scales(gravity)
    level = float(short_wave_level(u_star, gravity))
    if level <= 0 or k_from >= 15 * k_m:
        return np.zeros((2, 2))
    k = np.geomspace(k_from, 15 * k_m, points)
    c = phase_speed(k, gravity)
    b = 0.5 * level * c_m / c * np.exp(-0.25 * (k / k_m - 1) ** 2)
    delta = np.tanh(np.log(2) / 4 + 4 * (c / peak_speed) ** 2.5 + 0.13 * u_star / c_m * (c_m / c) ** 2.5)
    total = float(np.trapezoid(b, np.log(k)))
    difference = float(0.5 * np.trapezoid(delta * b, np.log(k)))
    along, across = (total + difference) / 2, (total - difference) / 2
    a = np.radians(toward_deg)
    c_, s_ = np.cos(a), np.sin(a)
    return np.array([[along * c_ * c_ + across * s_ * s_, (along - across) * c_ * s_],
                     [(along - across) * c_ * s_, along * s_ * s_ + across * c_ * c_]])


def wavevectors(size_m, cells):
    k = 2 * np.pi * np.fft.fftfreq(cells, d=size_m / cells)
    return np.meshgrid(k, k, indexing="xy")          # kx varies along columns (east), ky along rows (north)


def plan(cascades=CASCADES, gravity=None, choppy_limit=None):
    """Tiles with their wavenumber shares: each holds k up to pi / (2 cells) of its own grid, 4 per wavelength."""
    tiles, k_low = [], 0.0
    for size, cells in cascades:
        k_high = 2 * np.pi / (SAMPLES_PER_WAVELENGTH * size / cells)
        choppy = True if choppy_limit is None else k_high <= choppy_limit
        tiles.append(Tile(size, cells, k_low, k_high, choppy))
        k_low = k_high
    return tiles


def _periodic_sample(field_, rows, cols):
    return map_coordinates(field_, [rows, cols], order=1, mode="grid-wrap")


def realize(densities, depth_m, gravity, seed, cascades=CASCADES, geometry_tiles=1, chop=1.0,
            inversion_steps=6):
    """Draw the sea surface on every tile.

    densities: functions F(kx, ky) whose sum is the surface's wavenumber-vector density (m^4).
    geometry_tiles: how many of the largest tiles keep their heights (the rest keep slopes only).
    Gravity-wave tiles (all wavenumbers below a quarter of the capillary wavenumber) are displaced by the
    choppy wave model with the factor chop (1 is first-order Lagrangian motion).
    Returns the tiles with Eulerian slopes, heights where kept, and each tile's expected slope covariance and
    elevation variance (the spectrum's integral over its wavenumbers).
    """
    k_m, _ = capillary_scales(gravity)
    tiles = plan(cascades, choppy_limit=k_m / 4)
    rng = np.random.default_rng(seed)
    for index, tile in enumerate(tiles):
        kx, ky = wavevectors(tile.size_m, tile.cells)
        k = np.hypot(kx, ky)
        share = (k >= tile.k_low) & (k < tile.k_high) & (k > 0)
        dk2 = (2 * np.pi / tile.size_m) ** 2
        density = np.zeros_like(k)
        for d in densities:
            density[share] += d(kx[share], ky[share])
        power = density * dk2
        tile.expected_variance = float(power.sum())
        tile.expected_covariance = np.array([[np.sum(kx * kx * power), np.sum(kx * ky * power)],
                                             [np.sum(kx * ky * power), np.sum(ky * ky * power)]])
        # Complex Gaussian amplitudes with E|a|^2 = 2 F dk^2; the real part of their sum has variance sum F dk^2.
        amplitude = (rng.standard_normal(k.shape) + 1j * rng.standard_normal(k.shape)) * np.sqrt(power)
        n2 = tile.cells ** 2

        def synth(spectrum):
            return np.real(np.fft.ifft2(spectrum)) * n2

        gx, gy = synth(1j * kx * amplitude), synth(1j * ky * amplitude)
        keep_height = index < geometry_tiles
        height = synth(amplitude) if keep_height else None
        if not tile.choppy or chop == 0:
            tile.slope = np.stack([gx, gy]).astype(np.float32)
            tile.height = None if height is None else height.astype(np.float32)
            continue
        # Lagrangian displacement D = Re sum i k^ coth(k h) a e^{ik.X} and its Jacobian J = I + dD/dX.
        with np.errstate(divide="ignore", invalid="ignore"):
            factor = np.where(k > 0, chop / np.tanh(np.minimum(k * depth_m, 50.0)) / k, 0.0)
        dx_, dy_ = synth(1j * kx * factor * amplitude), synth(1j * ky * factor * amplitude)
        jxx = 1 - synth(kx * kx * factor * amplitude)
        jyy = 1 - synth(ky * ky * factor * amplitude)
        jxy = -synth(kx * ky * factor * amplitude)
        det = jxx * jyy - jxy * jxy
        det = np.where(np.abs(det) < 1e-3, np.copysign(1e-3, det), det)
        # Eulerian slope at each Lagrangian point: grad_x eta = J^-T grad_X eta.
        sx = (jyy * gx - jxy * gy) / det
        sy = (-jxy * gx + jxx * gy) / det
        # Invert x = X + D(X) on the regular grid by fixed-point iteration (grid units, periodic).
        cell = tile.spacing_m
        rows, cols = np.mgrid[0:tile.cells, 0:tile.cells].astype(float)
        X, Y = cols.copy(), rows.copy()
        for _ in range(inversion_steps):
            X = cols - _periodic_sample(dx_, Y, X) / cell
            Y = rows - _periodic_sample(dy_, Y, X) / cell
        tile.slope = np.stack([_periodic_sample(sx, Y, X), _periodic_sample(sy, Y, X)]).astype(np.float32)
        tile.height = None if height is None else _periodic_sample(height, Y, X).astype(np.float32)
    return tiles


def lean_pyramid(slope):
    """Mipmapped first and second slope moments for filtering (LEAN mapping, Olano and Baker 2010).

    Level 0 holds (sx, sy, sx^2, sy^2, sx sy) per cell; each level averages 2 by 2 cells of the one below,
    periodically, so that a level's covariance (second moments minus squared means) is the slope spread
    inside its cells.
    """
    sx, sy = slope.astype(np.float64)
    level = np.stack([sx, sy, sx * sx, sy * sy, sx * sy])
    levels = [level.astype(np.float32)]
    while level.shape[1] > 1:
        level = 0.25 * (level[:, 0::2, 0::2] + level[:, 1::2, 0::2] + level[:, 0::2, 1::2] + level[:, 1::2, 1::2])
        levels.append(level.astype(np.float32))
    return levels


def sample_covariance(tile):
    """Slope covariance of a drawn tile (the realization's own, to compare with its expectation)."""
    sx, sy = tile.slope.astype(np.float64)
    sx, sy = sx - sx.mean(), sy - sy.mean()
    return np.array([[np.mean(sx * sx), np.mean(sx * sy)], [np.mean(sx * sy), np.mean(sy * sy)]])
