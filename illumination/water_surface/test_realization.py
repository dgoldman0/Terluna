"""Checks of the drawn sea surface against the spectra it is drawn from."""
import numpy as np
import pytest

from illumination.water_surface import realization as r
from illumination.water_surface.model import slope_moments, wavenumber
from illumination.water_surface.short_waves import short_wave_slopes
from shared.constants import MOON_SURFACE_GRAVITY as G

SMALL = ((512.0, 512), (33.0, 256), (3.1, 256), (0.29, 256))


def swan_like(depth=80.0):
    """A wind sea at lunar gravity on SWAN's grid: JONSWAP in frequency, cos^2s spread about 30 degrees."""
    f = 0.03 * 1.1 ** np.arange(30)
    f = f[f <= 0.497]
    f = np.r_[f, 0.497]
    direction = np.arange(5.0, 360.0, 10.0)
    fp = 0.12
    sigma = np.where(f <= fp, 0.07, 0.09)
    jonswap = f ** -5 * np.exp(-1.25 * (fp / f) ** 4) * 3.3 ** np.exp(-(f - fp) ** 2 / (2 * sigma ** 2 * fp ** 2))
    jonswap *= 0.4 / np.trapezoid(jonswap, f)          # elevation variance 0.4 m2 (Hs 2.5 m)
    spread = np.cos(np.radians(direction - 30.0) / 2) ** 8
    spread /= spread.sum() * 10.0
    return f, direction, jonswap[:, None] * spread[None, :], depth


def test_tiles_share_the_spectrum_and_reproduce_the_study_slopes():
    f, direction, e, depth = swan_like()
    k_from = float(wavenumber(np.array([f[-1]]), depth, G)[0])
    u_star, toward = 0.15, 40.0
    densities = [r.swan_density(f, direction, e, depth, G), r.short_density(u_star, toward, G, k_from)]
    tiles = r.plan(SMALL)
    assert tiles[0].k_low == 0 and all(a.k_high == b.k_low for a, b in zip(tiles, tiles[1:]))
    drawn = r.realize(densities, depth, G, seed=1, cascades=SMALL)
    total = sum(t.expected_covariance for t in drawn) + r.residual_covariance(u_star, toward, G, drawn[-1].k_high)
    # The sea-slope study's integrals: the resolved spectrum and the short waves above its last frequency.
    ff, ee = f, e
    for _ in range(3):
        x = np.empty(2 * len(ff) - 1); x[::2] = ff; x[1::2] = (ff[1:] + ff[:-1]) / 2
        y = np.empty((len(x), ee.shape[1])); y[::2] = ee; y[1::2] = (ee[1:] + ee[:-1]) / 2
        ff, ee = x, y
    resolved = slope_moments(ff, direction, ee, depth, G)["covariance"]
    s = short_wave_slopes(u_star, G, k_from)
    a = np.radians(toward)
    unit = np.array([np.cos(a), np.sin(a)])
    across_unit = np.array([-np.sin(a), np.cos(a)])
    short = s["along"] * np.outer(unit, unit) + s["across"] * np.outer(across_unit, across_unit)
    study = resolved + short
    assert np.trace(total) == pytest.approx(np.trace(study), rel=0.03)
    assert np.linalg.eigvalsh(total) == pytest.approx(np.linalg.eigvalsh(study), rel=0.06)


def test_drawn_tiles_match_their_expected_statistics():
    f, direction, e, depth = swan_like()
    k_from = float(wavenumber(np.array([f[-1]]), depth, G)[0])
    densities = [r.swan_density(f, direction, e, depth, G), r.short_density(0.15, 40.0, G, k_from)]
    drawn = r.realize(densities, depth, G, seed=2, cascades=SMALL, chop=0.0)
    for tile in drawn:
        assert np.trace(r.sample_covariance(tile)) == pytest.approx(np.trace(tile.expected_covariance), rel=0.06)
    assert np.var(drawn[0].height) == pytest.approx(drawn[0].expected_variance, rel=0.15)


def test_lagrangian_displacement_sharpens_crests():
    # One wave of steepness kA = 0.15: a trochoid, whose Eulerian mean falls by kA^2/2 and whose heights skew
    # by 3/8 kA^4 / sigma^3 = 1.06 kA (crests narrow, troughs wide).
    size, cells = 256.0, 256
    k0, amplitude = 2 * np.pi / size * 16, 0.15 / (2 * np.pi / size * 16)

    def single(kx, ky):
        on = (np.abs(kx - k0) < 1e-9) & (np.abs(ky) < 1e-9)
        return np.where(on, amplitude ** 2 / 2 / (2 * np.pi / size) ** 2, 0.0)

    tile = r.realize([single], 1e4, G, seed=0, cascades=((size, cells),), chop=1.0)[0]
    h = tile.height.astype(float)
    rms = np.sqrt(2 * np.var(h))                          # the drawn amplitude (Rayleigh-distributed)
    assert h.mean() == pytest.approx(-k0 * rms ** 2 / 2, rel=0.1)
    assert np.mean((h - h.mean()) ** 3) / np.var(h) ** 1.5 == pytest.approx(1.06 * k0 * rms, rel=0.1)
    # A random sea: the mean level falls by the sum of k A^2 / 2 over the modes, as the choppy wave model has it.
    f, direction, e, depth = swan_like()
    density = r.swan_density(f, direction, e * 0.25, depth, G)
    choppy = r.realize([density], depth, G, seed=3, cascades=SMALL[:1], chop=1.0)[0]
    kx, ky = r.wavevectors(*SMALL[0])
    k = np.hypot(kx, ky)
    share = (k > 0) & (k < choppy.k_high)
    shift = -np.sum(k[share] * density(kx[share], ky[share])) * (2 * np.pi / SMALL[0][0]) ** 2
    assert choppy.height.mean() == pytest.approx(shift, rel=0.15)
    linear = r.realize([density], depth, G, seed=3, cascades=SMALL[:1], chop=0.0)[0]
    assert np.trace(r.sample_covariance(choppy)) == pytest.approx(np.trace(r.sample_covariance(linear)), rel=0.12)


def test_lean_levels_hold_the_slope_spread():
    f, direction, e, depth = swan_like()
    tile = r.realize([r.swan_density(f, direction, e, depth, G)], depth, G, seed=4, cascades=SMALL[:1])[0]
    levels = r.lean_pyramid(tile.slope)
    top = levels[-1][:, 0, 0].astype(float)
    assert top[2] - top[0] ** 2 + top[3] - top[1] ** 2 == pytest.approx(np.trace(r.sample_covariance(tile)), rel=1e-4)
    assert levels[0].shape == (5, 512, 512) and levels[1].shape == (5, 256, 256)
