"""Fresh water's optics and the light's passage through the surface and down the column.

The illumination domain's seawater optics (illumination/water_column/model.py) carry two Earth terms that the Open
Moon's fresh seas change; the domain's fresh_water.py corrects both, and this module reads its pure water:

- pure-water scattering: the domain's table is seawater's at 20 C and 38.4 g/kg (Zhang et al. 2009). Morel's
  measurements, which Zhang et al.'s fluctuation model reproduces within 1%, put seawater at 38.4 g/kg 30% above pure
  water, with the salt's share close to linear in salinity; at the seas' 0.14 g/kg the scattering is pure water's to
  0.1%, 0.770 of the table's;
- pure-water absorption below 420 nm: the table's values below 380 nm are Smith and Baker's (1981) bounds for the
  clearest natural waters (0.141 per metre at 300 nm), ten times pure water's. The IOCCG absorption protocol's
  composite (2018, Table 1.1) replaces them: Jonasz and Fournier's model from 250 to 340 nm and Morel et al.'s
  (2007) values from 340 to 415 nm, which join Pope and Fry at 420 nm.

Coloured dissolved organic matter keeps the domain's exponential slope, 0.0176 per nm, carried into the ultraviolet;
phytoplankton absorption stays at its 400 nm value below 400 nm, as the domain holds it; developed fines keep their
exponential. A steeper ultraviolet slope for dissolved organic matter is offered as a sensitivity.

At the surface a flat water passes the direct beam by Fresnel's equations for unpolarised light and refracts it; the
sky's diffuse light is treated as isotropic radiance split into zenith streams by Gauss-Legendre quadrature in the
cosine of the zenith angle, each passed by its own Fresnel transmission. Under the surface each stream fades with
Lee et al.'s (2005) diffuse attenuation evaluated at that stream's zenith above the water, so the slant streams of the
diffuse light fade faster than an overhead beam.
"""
from __future__ import annotations

import numpy as np

from illumination.water_column import fresh_water as fresh, model as domain

SEAWATER_TABLE_SALINITY = fresh.SEAWATER_TABLE_SALINITY
SALT_SCATTERING_AT_TABLE = fresh.SALT_SCATTERING_AT_TABLE
SEA_SALINITY = fresh.SEAS_SALINITY
PURE_WATER_UV = fresh.PURE_WATER_UV
JOIN_NM = fresh.JOIN_NM
WATER_INDEX = 1.333                 # pure water near 20-30 C in the visible
scattering_factor = fresh.scattering_factor
pure_water = fresh.pure_water


def dissolved_organic(wavelength_nm, absorption_440, uv_slope=None, uv_from_nm=350.0):
    """The domain's exponential; uv_slope (per nm) steepens it below uv_from_nm, a sensitivity."""
    w = np.asarray(wavelength_nm, float)
    a = domain.dissolved_organic(w, absorption_440)
    if uv_slope is None:
        return a
    knee = domain.dissolved_organic(uv_from_nm, absorption_440)
    return np.where(w < uv_from_nm, knee * np.exp(-uv_slope * (w - uv_from_nm)), a)


def water(wavelength_nm, chlorophyll_mg_m3=0.0, cdom_440=0.0, fines_g_m3=0.0, soil="developed",
          salinity_g_kg=SEA_SALINITY, uv_absorption=True, cdom_uv_slope=None):
    """Absorption and backscattering (per metre) of a water: the domain's constituents over fresh pure water."""
    w = np.asarray(wavelength_nm, float)
    a_w, b_w = pure_water(w, salinity_g_kg, uv_absorption)
    a_p, _, bb_p = domain.phytoplankton(w, chlorophyll_mg_m3)
    a_g = dissolved_organic(w, cdom_440, cdom_uv_slope)
    a_f = bb_f = np.zeros_like(w)
    if fines_g_m3 > 0:
        a_f, _, bb_f = domain.fines(w, fines_g_m3, soil)
    return dict(a=a_w + a_p + a_g + a_f, bb=0.5 * b_w + bb_p + bb_f)


def fresnel_reflectance(zenith_rad, index=WATER_INDEX):
    """Reflectance of a flat water surface for unpolarised light arriving at a zenith angle."""
    i = np.clip(np.asarray(zenith_rad, float), 1e-9, np.pi / 2 - 1e-9)
    t = np.arcsin(np.sin(i) / index)
    s = (np.sin(i - t) / np.sin(i + t)) ** 2
    p = (np.tan(i - t) / np.tan(i + t)) ** 2
    return 0.5 * (s + p)


def sky_streams(n=8):
    """Zenith angles (degrees above the water) and irradiance shares of an isotropic sky's streams."""
    x, wq = np.polynomial.legendre.leggauss(n)
    mu, wq = 0.5 * (x + 1), 0.5 * wq
    share = 2 * mu * wq                      # an isotropic sky's horizontal irradiance goes as 2 mu d(mu)
    return np.degrees(np.arccos(mu)), share / share.sum()


def diffuse_transmission(n=8, index=WATER_INDEX):
    """The share of an isotropic sky's horizontal irradiance a flat surface passes."""
    zenith, share = sky_streams(n)
    return float(np.sum(share * (1 - fresnel_reflectance(np.radians(zenith), index))))


def kd(iop, zenith_deg):
    """Lee et al.'s (2005) diffuse attenuation (per metre) for light arriving at a zenith angle above the water."""
    return domain.diffuse_attenuation(iop["a"], iop["bb"], zenith_deg)


class Column:
    """Light under a water for one set of spectral channels: direct and diffuse kernels against depth.

    below(e_direct, e_diffuse, sun_zenith_deg) gives the downwelling irradiance in each channel at each depth from the
    horizontal direct and diffuse irradiance on the water (any linear unit per channel)."""

    def __init__(self, wavelength_nm, depth_m, n_streams=8, **water_kw):
        self.depth = np.asarray(depth_m, float)
        self.iop = water(wavelength_nm, **water_kw)
        zenith, share = sky_streams(n_streams)
        transmit = 1 - fresnel_reflectance(np.radians(zenith))
        # Diffuse kernel (channels x depths): the streams' shares, transmissions and attenuations summed.
        self.diffuse = np.zeros((len(self.iop["a"]), len(self.depth)))
        for z, s, t in zip(zenith, share, transmit):
            self.diffuse += s * t * np.exp(-np.outer(kd(self.iop, z), self.depth))

    def direct(self, sun_zenith_deg):
        t = 1 - fresnel_reflectance(np.radians(sun_zenith_deg))
        return t * np.exp(-np.outer(kd(self.iop, sun_zenith_deg), self.depth))

    def below(self, e_direct, e_diffuse, sun_zenith_deg):
        out = np.asarray(e_diffuse, float)[:, None] * self.diffuse
        if np.any(np.asarray(e_direct) > 0):
            out = out + np.asarray(e_direct, float)[:, None] * self.direct(sun_zenith_deg)
        return out


def dna_damage(wavelength_nm):
    """Setlow's (1974) generalised DNA-damage action spectrum as the analytic form TEMIS uses (Bernhard and Seckmeyer
    1997), normalised to 1 at 300 nm. Its ultraviolet-A tail levels off at 3.1e-5 of the 300 nm value."""
    w = np.asarray(wavelength_nm, float)
    f = lambda x: np.exp(13.82 * (1 / (1 + np.exp((x - 310.0) / 9.0)) - 1))
    return f(w) / f(300.0)


def sun_elevation(lat_deg, hour_angle_deg, declination_deg=0.0):
    """The Sun's elevation (degrees) at a latitude and hour angle from local noon."""
    la, d = np.radians(lat_deg), np.radians(declination_deg)
    h = np.radians(np.asarray(hour_angle_deg, float))
    return np.degrees(np.arcsin(np.sin(la) * np.sin(d) + np.cos(la) * np.cos(d) * np.cos(h)))
