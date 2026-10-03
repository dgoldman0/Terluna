"""Resolved gravity-wave slopes for optical studies of a water surface.

The covariance follows from differentiating a linear random wave field. Its
frequency limits belong to the supplied spectrum. A glitter image additionally
needs a slope distribution, shorter waves, masking and an observer geometry.
"""
from __future__ import annotations

from datetime import datetime
from pathlib import Path

import numpy as np


def wavenumber(frequency_hz, depth_m, gravity_m_s2):
    """Solve omega**2 = g k tanh(k depth), with a bounded positive bracket."""
    f = np.asarray(frequency_hz, float)
    if (not np.isfinite(f).all() or np.any(f <= 0) or
            not np.isfinite([depth_m, gravity_m_s2]).all() or min(depth_m, gravity_m_s2) <= 0):
        raise ValueError('Positive finite frequencies, depth and gravity required')
    target = (2 * np.pi * f) ** 2 / gravity_m_s2
    low = np.zeros_like(f)
    high = target + np.sqrt(target / depth_m)
    while np.any(high * np.tanh(high * depth_m) < target):
        high *= 2
    for _ in range(64):
        middle = (low + high) / 2
        below = middle * np.tanh(middle * depth_m) < target
        low = np.where(below, middle, low)
        high = np.where(below, high, middle)
    return (low + high) / 2


def slope_moments(frequency_hz, direction_deg, variance_m2_hz_deg, depth_m, gravity_m_s2):
    """Integrate printed variance and k_i k_j E over a periodic directional grid.

    Directions are Cartesian propagation angles: zero east, 90 north. The
    frequency interpolant is linear in each moment density; integration uses
    the trapezoidal rule, with both end frequencies included and no tail.
    """
    f, angle, e = map(lambda x: np.asarray(x, float),
                      (frequency_hz, direction_deg, variance_m2_hz_deg))
    if (f.ndim != 1 or angle.ndim != 1 or len(f) < 2 or len(angle) < 3 or
            e.shape != (len(f), len(angle)) or not np.isfinite(e).all() or np.any(e < 0) or
            not np.isfinite(angle).all() or np.any(np.diff(f) <= 0) or
            not np.allclose(np.diff(angle), 360 / len(angle), rtol=0, atol=1e-5)):
        raise ValueError('Invalid directional variance grid')
    k = wavenumber(f, depth_m, gravity_m_s2)
    unit = np.column_stack((np.cos(np.radians(angle)), np.sin(np.radians(angle))))
    per_direction = np.trapezoid(e * k[:, None] ** 2, f, axis=0) * (360 / len(angle))
    covariance = (unit.T * per_direction) @ unit
    elevation_density = e.sum(axis=1) * (360 / len(angle))
    slope_density = elevation_density * k ** 2
    return dict(variance_m2=float(np.trapezoid(elevation_density, f)), covariance=covariance,
                frequency_hz=f, elevation_density=elevation_density, slope_density=slope_density,
                dispersion_relative_residual=float(np.max(abs(
                    gravity_m_s2 * k * np.tanh(k * depth_m) / (2 * np.pi * f) ** 2 - 1))))


def fraction_above(frequency, density, cutoff):
    """Fraction above an exact cutoff in a nonnegative piecewise-linear density."""
    f, y = np.asarray(frequency), np.asarray(density)
    total = np.trapezoid(y, f)
    if total == 0:
        return None
    if cutoff <= f[0]:
        return 1.0
    if cutoff >= f[-1]:
        return 0.0
    keep = f > cutoff
    x = np.r_[cutoff, f[keep]]
    v = np.r_[np.interp(cutoff, f, y), y[keep]]
    return float(np.trapezoid(v, x) / total)


def iter_swan(path):
    """Stream time-dependent SWAN v1 Cartesian-direction variance spectra.

    This reader accepts the format used by the lunar-seas coastal product. Dry
    records carry an explicit availability mask. AFREQ and RFREQ are returned
    to the caller, which must establish the current convention before use.
    """
    with Path(path).open() as stream:
        lines = (s.strip() for s in stream if s.strip() and not s.lstrip().startswith('$'))

        def line():
            try:
                return next(lines)
            except StopIteration as exc:
                raise ValueError('Truncated SWAN spectrum') from exc

        def expect(word):
            if line().split()[0] != word:
                raise ValueError(f'Expected SWAN field {word}')

        if line().split()[:2] != ['SWAN', '1']:
            raise ValueError('SWAN v1 required')
        expect('TIME')
        if int(line().split()[0]) != 1:
            raise ValueError('Unsupported SWAN time coding')
        expect('LONLAT')
        nloc = int(line().split()[0])
        locations = np.array([[float(v) for v in line().split()] for _ in range(nloc)])
        reference = line().split()[0]
        if reference not in ('AFREQ', 'RFREQ'):
            raise ValueError('Unsupported frequency reference')
        nf = int(line().split()[0])
        frequency = np.array([float(line().split()[0]) for _ in range(nf)])
        expect('CDIR')
        nd = int(line().split()[0])
        direction = np.array([float(line().split()[0]) for _ in range(nd)])
        expect('QUANT')
        if int(line().split()[0]) != 1:
            raise ValueError('One variance quantity required')
        expect('VaDens')
        expect('m2/Hz/degr')
        exception = float(line().split()[0])
        if (nloc < 1 or locations.shape != (nloc, 2) or not np.isfinite(locations).all() or
                nf < 2 or nd < 3 or not np.isfinite(frequency).all() or
                np.any(frequency <= 0) or np.any(np.diff(frequency) <= 0) or
                not np.isfinite(direction).all() or
                not np.allclose(np.diff(direction), 360 / nd, rtol=0, atol=1e-5) or
                not np.isfinite(exception) or exception >= 0):
            raise ValueError('Invalid SWAN axes or exception value')
        previous = None
        for stamp in lines:
            time = datetime.strptime(stamp.split()[0], '%Y%m%d.%H%M%S')
            if previous is not None and time <= previous:
                raise ValueError('SWAN times must increase')
            previous = time
            values = np.zeros((nloc, nf, nd))
            available = np.ones(nloc, bool)
            for i in range(nloc):
                mode = line()
                if mode == 'NODATA':
                    available[i] = False
                elif mode == 'FACTOR':
                    scale = float(line())
                    table = np.array([[float(v) for v in line().split()] for _ in range(nf)])
                    if (not np.isfinite(scale) or scale <= 0 or table.shape != (nf, nd) or
                            not np.isfinite(table).all() or np.any(table < 0)):
                        raise ValueError('Invalid SWAN variance table')
                    values[i] = table * scale
                elif mode != 'ZERO':
                    raise ValueError(f'Unsupported SWAN record {mode}')
            yield dict(time=time, locations=locations, frequency=frequency, direction=direction,
                       frequency_reference=reference, variance=values, available=available)
        if previous is None:
            raise ValueError('Empty SWAN record sequence')
