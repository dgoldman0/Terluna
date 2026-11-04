"""The Moon's tide as forcing for the ocean models.

Geography's tide model (geography/tides.py) writes Earth's and the Sun's
degree-2 and degree-3 pull as 20 monomials of the unit position, with
coefficients through time from the lunar ephemeris (degree 2 carries the
Moon's Love numbers, 1 + k2 - h2). The equilibrium height of the forcing at a
point is coefficients(t) @ monomials(point).

The water's own gravity adds to it. Each monomial's self-attraction is
computed once from the seas' equilibrium response to that monomial (each sea's
volume held, spherical harmonics to degree 60, fresh water at the seas'
warmth) and added to the monomial, so the forcing basis is B_k = m_k + S_k.
The dynamic models then feel the equilibrium tide's self-attraction; their
own departures from equilibrium (within 7% on the main seas) change it by a
few per cent of 9%.

Model time zero is 2026-01-01 00:00 TT, the tide record's start.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from climate.ocean import lunar

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SAL_PRODUCT = HERE / "products" / "tide_self_attraction_4ppd.npz"
SAL_SCHEMA = "terluna.climate.ocean-tide-self-attraction/1"
SAL_BOX = dict(lat=(-40.0, 80.0), lon=(-95.0, 85.0))


def unit_points(lat_deg, lon_deg):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)), axis=-1)


def monomials(lat_deg, lon_deg):
    """The 20 monomials at points (..., 20)."""
    from geography import tides
    return tides.monomial_values(unit_points(lat_deg, lon_deg))


def coefficients(seconds):
    """Monomial coefficients (time, 20), metres of equilibrium height, at model times in seconds."""
    from geography import tides
    return tides.forcing_coefficients(tides.JD_START + np.asarray(seconds, dtype=float) / 86400.0)


def lines(count=12):
    """The forcing's strongest lines over 2026-2045: [(line, complex monomial amplitude (20,))]."""
    from geography import tides
    days = np.arange(0.0, tides.YEARS * 365.25, tides.STEP_DAYS)
    found, _, amplitude, _ = tides.constituents(days, tides.forcing_coefficients(tides.JD_START + days), count=count)
    return list(zip(found, amplitude))


def line(argument):
    for entry, amplitude in lines():
        if entry["argument"].replace(" ", "") == argument.replace(" ", ""):
            return entry, amplitude
    raise KeyError(argument)


# ------------------------------------------------- self-attraction

def build_self_attraction(path=SAL_PRODUCT):
    """Each monomial's self-attraction (m of equilibrium height) on the atlas's quarter-degree nodes in a box
    round the near side's seas; written once (about two minutes)."""
    from geography import tides
    grid = tides.load_seas()
    bodies = tides.sea_bodies()
    response = tides.Response(grid, bodies.keys(), water_density=lunar.WATER_DENSITY)
    lat, lon = grid["lat"], grid["lon"]
    la, lo = np.meshgrid(lat, lon, indexing="ij")
    values = tides.monomial_values(unit_points(la, lo))
    fields = []
    for k in range(values.shape[-1]):
        height, _ = response.solve(values[..., k])
        fields.append(response.self_attraction(height))
    fields = np.array(fields)                                   # (20, lat, lon), lat north first, lon 0..360
    lon180 = np.where(lon > 180, lon - 360, lon)
    order = np.argsort(lon180)
    rows = np.flatnonzero((lat >= SAL_BOX["lat"][0]) & (lat <= SAL_BOX["lat"][1]))[::-1]
    cols = order[(lon180[order] >= SAL_BOX["lon"][0]) & (lon180[order] <= SAL_BOX["lon"][1])]
    box = fields[:, rows][:, :, cols].astype(np.float32)
    metadata = dict(schema=SAL_SCHEMA, monomials=list(tides.MONOMIALS), water_density_kg_m3=lunar.WATER_DENSITY,
                    harmonic_degree=tides.DEGREE, atlas_level_m=grid["level_m"], atlas_sha256=grid["source_sha256"],
                    producer={str(p.relative_to(ROOT)): tides.sha256(p) for p in (Path(__file__), ROOT / "geography/tides.py")},
                    reading_rule="sal[k, j, i]: equilibrium height (m) added by the seas' own attraction when the forcing "
                                 "is monomial k, at lat_deg[j] (south first) and lon_deg[i] (-180..180); add to the monomial.")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, sal=box, lat_deg=lat[rows], lon_deg=lon180[cols], metadata=np.array(json.dumps(metadata)))
    return Path(path)


def self_attraction(lat_deg, lon_deg, path=SAL_PRODUCT):
    """Each monomial's self-attraction at points (..., 20), bilinear from the product."""
    from scipy.interpolate import RegularGridInterpolator
    if not Path(path).exists():
        build_self_attraction(path)
    with np.load(path, allow_pickle=False) as source:
        if json.loads(source["metadata"].item())["schema"] != SAL_SCHEMA:
            raise ValueError("unsupported self-attraction product")
        sal, lat, lon = source["sal"].astype(float), source["lat_deg"], source["lon_deg"]
    lat_deg, lon_deg = np.broadcast_arrays(np.asarray(lat_deg, float), np.asarray(lon_deg, float))
    points = np.stack((lat_deg.ravel(), lon_deg.ravel()), axis=-1)
    out = np.empty(lat_deg.shape + (sal.shape[0],))
    for k in range(sal.shape[0]):
        f = RegularGridInterpolator((lat, lon), sal[k], bounds_error=True)
        out[..., k] = f(points).reshape(lat_deg.shape)
    return out


def basis(lat_deg, lon_deg, attraction=True):
    """The forcing basis at points (..., 20): monomials, plus their self-attraction."""
    b = monomials(lat_deg, lon_deg)
    return b + self_attraction(lat_deg, lon_deg) if attraction else b


# ------------------------------------------------- time series

def ramp(seconds, ramp_days):
    """tanh ramp from rest, as SCHISM's dramp: tanh(2 t / ramp)."""
    if ramp_days <= 0:
        return np.ones_like(np.asarray(seconds, float))
    return np.tanh(2.0 * np.asarray(seconds, float) / (ramp_days * 86400.0))


def ephemeris_series(seconds):
    """Coefficients through time from the ephemeris, minus their 2026-2045 mean (the mean is the static tide the
    geoid already holds)."""
    from geography import tides
    days = np.arange(0.0, tides.YEARS * 365.25, tides.STEP_DAYS)
    mean = tides.forcing_coefficients(tides.JD_START + days).mean(axis=0)
    return coefficients(seconds) - mean


def line_series(seconds, amplitude, frequency_cpd):
    """Coefficients through time for one line: Re[A exp(i w t)]."""
    w = 2.0 * math.pi * frequency_cpd / 86400.0
    return np.real(np.asarray(amplitude)[None, :] * np.exp(1j * w * np.asarray(seconds, float))[:, None])


def write_mitgcm(path, grid_spec, seconds, coefs, attraction=True):
    """MITgcm's tidal forcing file: equilibrium height F (m) on a regular grid, one record per time
    (big-endian float64; the namelist scales it by -g into the geopotential)."""
    lat = grid_spec["lat0"] + grid_spec["step"] * np.arange(grid_spec["nlat"])
    lon = grid_spec["lon0"] + grid_spec["step"] * np.arange(grid_spec["nlon"])
    la, lo = np.meshgrid(lat, lon, indexing="ij")
    b = basis(la, lo, attraction)                                          # (nlat, nlon, 20)
    with open(path, "wb") as out:
        for c in coefs:
            np.asarray(b @ c, dtype=">f8").tofile(out)
    return Path(path)


def write_schism(folder, lat_nodes, lon_nodes, seconds, coefs, attraction=True):
    """SCHISM's lunar tide files: lunar_tide.th (coefficients through time) and lunar_tide_basis.in
    (each node's basis), for ilunartide = 1 of the lunar patch."""
    folder = Path(folder)
    step = float(seconds[1] - seconds[0])
    if not np.allclose(np.diff(seconds), step) or abs(seconds[0]) > 1e-9:
        raise ValueError("SCHISM's lunar tide records start at zero and are evenly spaced")
    with open(folder / "lunar_tide.th", "w") as out:
        out.write(f"{len(coefs)} {coefs.shape[1]} {step:.6f}\n")
        np.savetxt(out, coefs, fmt="%.9e")
    b = basis(np.asarray(lat_nodes), np.asarray(lon_nodes), attraction)
    with open(folder / "lunar_tide_basis.in", "w") as out:
        out.write(f"{b.shape[0]} {b.shape[1]}\n")
        np.savetxt(out, b, fmt="%.9e")
