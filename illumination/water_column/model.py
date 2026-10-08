"""The colour of seawater: what its constituents absorb and scatter, and the reflectance that follows.

Absorption a, scattering b and backscattering b_b (per metre) add over the water's constituents:

- pure seawater: absorption from Pope and Fry (1997) and Kou et al. (1993), scattering from Zhang et al. (2009)
  at 20 C and 38.4 g/kg, as NASA's ocean-colour software tabulates them; half of it is backward;
- phytoplankton and the particles that vary with them (Earth's Case 1 water): absorption A Chl^E of Bricaud et al.
  (1998), 400-700 nm, held at its 400 nm value below and left out above, where water absorbs a hundred times
  more; scattering 0.416 Chl^0.766 at 550 nm (Loisel and Morel 1998) with Morel and Maritorena's (2001)
  spectral slope and backscattering ratio;
- coloured dissolved organic matter: a_g(440) exp(-0.0176 (lambda - 440)), the mean slope of Babin et al. (2003);
- suspended regolith fines: the single-scattering albedo of lunar soil's finest fraction (under 10 micrometres)
  from the Lunar Soil Characterization Consortium's reflectance through Hapke's isotropic model, carried into
  water through Hapke's equivalent slab (the grains' internal transmission stays, their surface reflection
  follows the index relative to water), for equant grains of a stated size and density. Diffraction adds
  forward scattering equal to the geometric cross-section; a share of the geometrically scattered light,
  0.06, goes backward, which gives the backscattering ratios of 0.015-0.03 measured for mineral particles.

The remote-sensing reflectance is that of Lee et al. (2002): u = b_b/(a + b_b), r_rs = u (0.089 + 0.125 u) below
the surface and R_rs = 0.52 r_rs / (1 - 1.7 r_rs) above it, per steradian. The downwelling light's diffuse
attenuation follows Lee et al. (2005): K_d = (1 + 0.005 theta_s) a + 4.18 (1 - 0.52 exp(-10.8 a)) b_b.
"""
from __future__ import annotations

import re

import numpy as np
from scipy.optimize import brentq

from illumination.water_column.fetch_inputs import path

RELAB_INCIDENCE_DEG, RELAB_EMISSION_DEG = 30.0, 0.0
CDOM_SLOPE = 0.0176
BACKWARD_SHARE = 0.06
# Real index and grain density: mare soils of pyroxene, glass and ilmenite near 1.7 and 3.1 g/cm3; highland
# soils of plagioclase near 1.6 and 2.9 g/cm3 (Carrier, Olhoeft and Mendell 1991, Lunar Sourcebook, ch. 9).
SOILS = {"12001": dict(name="Apollo 12 mare soil, Oceanus Procellarum", index=1.70, density_g_m3=3.1e6),
         "10084": dict(name="Apollo 11 high-titanium mare soil, Mare Tranquillitatis", index=1.70, density_g_m3=3.1e6),
         "62231": dict(name="Apollo 16 highland soil, Descartes", index=1.60, density_g_m3=2.9e6)}


def _table(name, header_lines=0, delimiter=None):
    rows = []
    for line in path(name).read_text(encoding="latin-1").splitlines()[header_lines:]:
        line = line.strip()
        if line and not line.startswith(("/", "!")) and line[0].isdigit():
            rows.append([float(v) for v in (line.split(delimiter) if delimiter else line.split())])
    return np.array(rows)


def seawater(wavelength_nm):
    """Absorption and scattering of pure seawater (per metre)."""
    t = _table("water_spectra.dat")
    w = np.asarray(wavelength_nm, float)
    return np.interp(w, t[:, 0], t[:, 1]), np.interp(w, t[:, 0], t[:, 2])


def phytoplankton(wavelength_nm, chlorophyll_mg_m3):
    """Absorption, scattering and backscattering of phytoplankton and the particles that vary with them."""
    t = _table("aph_bricaud_1998.txt", delimiter=",")
    w = np.asarray(wavelength_nm, float)
    chl = float(chlorophyll_mg_m3)
    if chl <= 0:
        zero = np.zeros_like(w)
        return zero, zero, zero
    absorption = np.interp(w, t[:, 0], t[:, 1]) * chl ** np.interp(w, t[:, 0], t[:, 2])
    absorption = np.where(w > t[-1, 0], 0.0, absorption)
    slope = 0.5 * (np.log10(chl) - 0.3) if chl < 2 else 0.0
    scattering = 0.416 * chl ** 0.766 * (w / 550.0) ** slope
    ratio = 0.002 + 0.01 * (0.5 - 0.25 * np.log10(chl)) * (w / 550.0) ** slope
    return absorption, scattering, ratio * scattering


def dissolved_organic(wavelength_nm, absorption_440):
    return absorption_440 * np.exp(-CDOM_SLOPE * (np.asarray(wavelength_nm, float) - 440.0))


def soil_reflectance(sample, fraction=1):
    """Wavelength (nm) and RELAB reflectance of an LSCC soil: fraction 1 is <10 um, 3 is 10-20, 5 is 20-45, 7 is <45."""
    text = path(f"{sample}Kdata.txt").read_text(encoding="latin-1")
    rows = []
    for line in re.split(r"\r\n|\r|\n", text)[1:]:
        cells = line.split("\t")
        if len(cells) > fraction and cells[0].strip() and cells[fraction].strip():
            rows.append((float(cells[0]), float(cells[fraction])))
    return np.array(rows).T


def _hapke_h(x, albedo):
    """Hapke's (2002) approximation to Chandrasekhar's H function for isotropic scatterers."""
    gamma = np.sqrt(1 - albedo)
    r0 = (1 - gamma) / (1 + gamma)
    return 1 / (1 - albedo * x * (r0 + (1 - 2 * r0 * x) / 2 * np.log((1 + x) / x)))


def reflectance_factor(albedo, incidence_deg=RELAB_INCIDENCE_DEG, emission_deg=RELAB_EMISSION_DEG):
    """Bidirectional reflectance factor of a half-space of isotropic scatterers, without the opposition effect."""
    mu0, mu = np.cos(np.radians(incidence_deg)), np.cos(np.radians(emission_deg))
    return albedo / 4 * _hapke_h(mu0, albedo) * _hapke_h(mu, albedo) / (mu0 + mu)


def single_scattering_albedo(reflectance):
    """Invert the isotropic Hapke reflectance factor at RELAB's geometry."""
    return np.array([brentq(lambda w: reflectance_factor(w) - r, 1e-9, 1 - 1e-12) for r in np.atleast_1d(reflectance)])


def _slab(index):
    external = ((index - 1) / (index + 1)) ** 2 + 0.05
    internal = 1 - 4 / (index * (index + 1) ** 2)
    path_length = 2 / 3 * (index ** 2 - (index ** 2 - 1) ** 1.5 / index)
    return external, internal, path_length


def albedo_in_water(albedo_in_air, grain_index, water_index):
    """Carry a grain's single-scattering albedo from air into water through Hapke's equivalent slab.

    The internal transmission Theta = exp(-alpha <D>) follows from the albedo in air; in water the same grain
    keeps its absorption coefficient, its mean internal path scales with the relative index, and its surface
    reflects as that index gives. Grains darker than their own surface reflection count as opaque.
    """
    se, si, d_air = _slab(grain_index)
    theta = np.clip((albedo_in_air - se) / ((1 - se) * (1 - si) + si * (albedo_in_air - se)), 0.0, 1.0)
    relative = grain_index / np.asarray(water_index, float)
    se_w, si_w, d_water = _slab(relative)
    theta_w = theta ** (d_water / d_air)
    return se_w + (1 - se_w) * (1 - si_w) * theta_w / (1 - si_w * theta_w)


def fines(wavelength_nm, concentration_g_m3, sample="12001", diameter_um=4.0, water_index=1.34):
    """Absorption, scattering and backscattering of suspended regolith fines (per metre)."""
    soil = SOILS[sample]
    w = np.asarray(wavelength_nm, float)
    lam, refl = soil_reflectance(sample)
    albedo_air = np.interp(w, lam, single_scattering_albedo(refl))
    albedo = albedo_in_water(albedo_air, soil["index"], water_index)
    cross_section = concentration_g_m3 * 3 / (2 * soil["density_g_m3"] * diameter_um * 1e-6)
    return cross_section * (1 - albedo), cross_section * (1 + albedo), cross_section * albedo * BACKWARD_SHARE


def water(wavelength_nm, chlorophyll_mg_m3=0.0, cdom_440=0.0, fines_g_m3=0.0, soil="12001",
          fines_diameter_um=4.0):
    """Total and per-constituent inherent optical properties of a water."""
    w = np.asarray(wavelength_nm, float)
    a_w, b_w = seawater(w)
    a_p, b_p, bb_p = phytoplankton(w, chlorophyll_mg_m3)
    a_g = dissolved_organic(w, cdom_440)
    if fines_g_m3 > 0:
        a_f, b_f, bb_f = fines(w, fines_g_m3, soil, fines_diameter_um)
    else:
        a_f = b_f = bb_f = np.zeros_like(w)
    parts = dict(seawater=(a_w, b_w, 0.5 * b_w), phytoplankton=(a_p, b_p, bb_p),
                 dissolved_organic=(a_g, np.zeros_like(w), np.zeros_like(w)), fines=(a_f, b_f, bb_f))
    a = sum(p[0] for p in parts.values())
    b = sum(p[1] for p in parts.values())
    bb = sum(p[2] for p in parts.values())
    return dict(a=a, b=b, bb=bb, parts=parts)


def remote_sensing_reflectance(a, bb):
    """Above-surface remote-sensing reflectance (per steradian) of an optically deep water, Lee et al. (2002)."""
    u = np.asarray(bb, float) / (np.asarray(a, float) + bb)
    below = u * (0.089 + 0.125 * u)
    return 0.52 * below / (1 - 1.7 * below)


def diffuse_attenuation(a, bb, sun_zenith_deg):
    """Diffuse attenuation of downwelling light (per metre), Lee et al. (2005)."""
    a = np.asarray(a, float)
    return (1 + 0.005 * sun_zenith_deg) * a + 4.18 * (1 - 0.52 * np.exp(-10.8 * a)) * bb
