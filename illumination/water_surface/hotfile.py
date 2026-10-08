"""Directional wave spectra from SWAN's restart files.

A nonstationary SWAN run that writes a restart file at the end of each segment
keeps the full directional spectrum of every grid point at those times, where
its printed spectra cover only chosen points.
"""
from __future__ import annotations

import struct
from datetime import datetime
from pathlib import Path

import numpy as np


def read_hotfile(path, nodes=None):
    """Directional variance spectra from a SWAN 41.51 unformatted hotfile of one time.

    The Fortran records hold the version and project, the grid's point count and
    its two dimensions, every point's longitude and latitude, the frequencies (Hz)
    and Cartesian directions (degrees), the time, then each point's action density
    N(sigma, theta) in SWAN's units (m2 s per rad/s per rad), frequency by
    direction. Variance density in m2/Hz/degree is sigma N times 2 pi and pi/180.
    nodes, an array of (longitude, latitude), selects the nearest grid points.
    """
    data = Path(path).read_bytes()
    offsets, position = [], 0
    while position < len(data):
        (length,) = struct.unpack('<i', data[position:position + 4])
        end = position + 4 + length
        if length < 0 or end + 4 > len(data) or struct.unpack('<i', data[end:end + 4])[0] != length:
            raise ValueError('Broken Fortran record in SWAN hotfile')
        offsets.append((position + 4, length))
        position = end + 4

    def record(i):
        start, length = offsets[i]
        return data[start:start + length]

    count, nx, ny = struct.unpack('<3i', record(5))
    if not record(0).decode().strip().startswith('41.51') or count != nx * ny:
        raise ValueError('SWAN 41.51 hotfile of a regular grid required')
    lonlat = np.array([struct.unpack('<2f', record(6 + i)) for i in range(count)], float)
    i = 6 + count
    nf = struct.unpack('<i', record(i))[0]
    frequency = np.array([struct.unpack('<f', record(i + 1 + j))[0] for j in range(nf)], float)
    i += 1 + nf
    nd = struct.unpack('<i', record(i))[0]
    # Single precision blurs SWAN's exactly spaced bin centres by up to 2e-5 degrees.
    direction = np.round([struct.unpack('<f', record(i + 1 + j))[0] for j in range(nd)], 4)
    i += 1 + nd
    time = datetime.strptime(record(i).decode().split()[0], '%Y%m%d.%H%M%S')
    i += 1
    if (len(offsets) != i + 2 * count or np.any(np.diff(frequency) <= 0) or
            not np.allclose(np.diff(direction), 360 / nd, rtol=0, atol=1e-3) or
            any(offsets[i + 2 * p + 1][1] != 4 * nf * nd for p in range(count))):
        raise ValueError('Unexpected SWAN hotfile layout')
    if nodes is None:
        chosen = np.arange(count)
    else:
        nodes = np.atleast_2d(np.asarray(nodes, float))
        chosen = np.array([int(np.argmin(np.hypot(lonlat[:, 0] - x, lonlat[:, 1] - y))) for x, y in nodes])
    action = np.stack([np.frombuffer(record(i + 2 * p + 1), '<f4').reshape(nf, nd) for p in chosen]).astype(float)
    variance = action * (2 * np.pi * frequency)[None, :, None] * 2 * np.pi * np.pi / 180
    return dict(time=time, lonlat=lonlat[chosen], frequency=frequency, direction=direction,
                variance=variance, index=chosen)
