"""The Earth's light at the Moon by wavelength: Glenar et al.'s spectrum at Robinson et al.'s visual brightness.

Glenar et al. (2019, Icarus 321, 841) condensed NASA Virtual Planetary Laboratory models of the whole Earth,
March and April 2008, into the Earth's diurnally averaged radiance over its visible disk, seen from the Moon,
as a function of wavelength (0.30-2.49 micrometres), phase angle g and the sub-observer latitude theta:

    I(g, theta) = I00 PHI(g, b) [1 + A_L (fL(theta) - 0.190)]

with PHI a symmetric double Henyey-Greenstein phase curve of width b, normalized to one at full phase, and
fL(theta) = 0.181 + 7.578e-4 theta + 1.1517e-5 theta^2 the mean land fraction of the visible disk. Their fit
holds for phase angles below 60 degrees and sub-observer latitudes within 30 degrees of the equator, within
10-12% at 0.3 micrometres and more in the near infrared; beyond 60 degrees the curve is extrapolated. The
day-to-day turn of the continents and clouds (their diurnal term) is left out. The Earth's irradiance at the
Moon is I times the solid angle of its disk.

The model's Earth is fainter than the observed one. Robinson et al. (2025, Planet. Sci. J., arXiv:2507.22258)
combined ground and spacecraft photometry from 5 to 144 degrees of phase into the Earth's visual (0.4-0.7
micrometre) phase curve: a geometric albedo of 0.242, against the long-quoted 0.367 that rests on Danjon's
extrapolated earthshine, and an analytic fit I/F = f P_HG(180 - alpha; g) / (pi P_HG(180; g)) with f = 0.23 and
g = -0.33. calibrated_irradiance keeps Glenar et al.'s spectral shape and scales it, at each phase angle, to
that visual brightness, weighting by the solar spectrum over 400-700 nm. Beyond the 144 degrees their data
reach, the scale stays at its 144-degree value and Glenar et al.'s curve carries the Earth to dark at new Earth.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np

from illumination.earthlight.fetch_inputs import path
from shared.constants import EARTH_RADIUS

TABLE = "Earthshine_Spring_Approximation_Revised_Sept2018.txt"
VALID_PHASE_DEG = 60.0
VISUAL_NORMALISATION, VISUAL_ASYMMETRY = 0.23, -0.33     # Robinson et al. (2025), equation 14
VISUAL_BAND_NM = (400.0, 700.0)
OBSERVED_PHASE_DEG = 144.0       # the largest phase angle in Robinson et al.'s data


@lru_cache(maxsize=1)
def coefficients():
    """Wavelength (nm), zero-phase radiance I00 (W m-2 sr-1 nm-1), width b and land coefficient A_L (read-only)."""
    rows = []
    started = False
    for line in path(TABLE).read_text(encoding="latin-1").splitlines():
        if line.strip().startswith("Wvl"):
            started = True
            continue
        if started and line.strip():
            rows.append([float(v) for v in line.split()])
    t = np.array(rows)
    if t.shape[1] != 11:
        raise ValueError("Unexpected layout of the earthshine table")
    return dict(wavelength_nm=t[:, 0] * 1000, radiance_00=t[:, 1] / 1000, width=t[:, 3], land=t[:, 5],
                uncertainty=t[:, 10])


def _lobes(cos_g, b):
    return 0.5 * (1 - b * b) / (1 - 2 * b * cos_g + b * b) ** 1.5 + 0.5 * (1 - b * b) / (1 + 2 * b * cos_g + b * b) ** 1.5


def phase_curve(phase_deg, width):
    """PHI(g, b): the double Henyey-Greenstein light curve, one at full phase."""
    c = np.cos(np.radians(phase_deg))
    return 0.5 * _lobes(c, width) / _lobes(1.0, width) * (1 + c)


def land_fraction(latitude_deg):
    return 0.181 + 7.578e-4 * latitude_deg + 1.1517e-5 * latitude_deg ** 2


def radiance(wavelength_nm, phase_deg, sub_latitude_deg=0.0):
    """The Earth's disk-averaged radiance (W m-2 sr-1 nm-1) at the given wavelengths."""
    c = coefficients()
    i = (c["radiance_00"] * phase_curve(phase_deg, c["width"])
         * (1 + c["land"] * (land_fraction(sub_latitude_deg) - 0.190)))
    return np.interp(np.asarray(wavelength_nm, float), c["wavelength_nm"], i)


def band_radiance(edges_nm, phase_deg, sub_latitude_deg=0.0):
    """Mean radiance over each band [low, high] (nm), from the table's points within it."""
    c = coefficients()
    i = (c["radiance_00"] * phase_curve(phase_deg, c["width"])
         * (1 + c["land"] * (land_fraction(sub_latitude_deg) - 0.190)))
    out = []
    for low, high in np.atleast_2d(edges_nm):
        grid = np.linspace(low, high, 41)
        out.append(np.trapezoid(np.interp(grid, c["wavelength_nm"], i), grid) / (high - low))
    return np.array(out)


def solid_angle(distance_m):
    """Solid angle of the Earth's disk from a given distance (sr)."""
    return 2 * np.pi * (1 - np.sqrt(1 - (EARTH_RADIUS / np.asarray(distance_m, float)) ** 2))


def visual_phase_curve(phase_deg, asymmetry=VISUAL_ASYMMETRY):
    """The Earth's visual phase function, one at full phase: a Henyey-Greenstein lobe at scattering angle 180 - alpha."""
    def lobe(cos_theta):
        return (1 - asymmetry ** 2) / (1 + asymmetry ** 2 - 2 * asymmetry * cos_theta) ** 1.5
    return lobe(-np.cos(np.radians(phase_deg))) / lobe(-1.0)


def visual_intensity(phase_deg):
    """The Earth's visual disk-averaged intensity over the solar flux, I/F (per steradian)."""
    return VISUAL_NORMALISATION * visual_phase_curve(phase_deg) / np.pi


def calibrated_irradiance(bands_nm, solar_band_irradiance, phase_deg, distance_m, sub_latitude_deg=0.0):
    """The Earth's irradiance at the Moon in each band (W/m2): Glenar's spectrum scaled to Robinson's visual level.

    solar_band_irradiance is the unshielded sunlight at the Earth's distance from the Sun in the same bands (W/m2).
    """
    bands = np.atleast_2d(np.asarray(bands_nm, float))
    width = bands[:, 1] - bands[:, 0]
    visual = (bands[:, 0] >= VISUAL_BAND_NM[0]) & (bands[:, 1] <= VISUAL_BAND_NM[1])
    anchor = min(float(phase_deg), OBSERVED_PHASE_DEG)
    modelled = ((band_radiance(bands, anchor, sub_latitude_deg) * width)[visual].sum()
                / np.asarray(solar_band_irradiance, float)[visual].sum())
    intensity = band_radiance(bands, phase_deg, sub_latitude_deg) * width          # W m-2 sr-1 per band
    return intensity * visual_intensity(anchor) / modelled * solid_angle(distance_m)
