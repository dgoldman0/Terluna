"""Pure water as the Open Moon's fresh seas hold it: model.py's constituents over fresh water's own optics.

model.py's pure water is seawater as Earth's ocean-colour tables give it, and two of its terms change in the Open
Moon's fresh seas (about 0.14 g/kg):

- scattering: the table's is seawater's at 20 C and 38.4 g/kg (Zhang et al. 2009). Morel's measurements, which
  Zhang et al.'s fluctuation model reproduces within 1%, put seawater at 38.4 g/kg 30% above pure water, with the
  salt's share close to linear in salinity; at the seas' 0.14 g/kg the scattering is pure water's to 0.1%, 0.770 of
  the table's;
- absorption below 420 nm: the table's values below 380 nm are Smith and Baker's (1981) bounds for the clearest
  natural waters (0.141 per metre at 300 nm), ten times pure water's. The IOCCG absorption protocol's composite
  (2018, Table 1.1) replaces them: Jonasz and Fournier's model from 250 to 340 nm and Morel et al.'s (2007) values
  from 340 to 415 nm, which join Pope and Fry at 420 nm.

water() is model.water() over this pure water. Callers that keep model.py's seawater (the sea-appearance study's
waters, regimes and scenes) move here when they are next rerun.
"""
from __future__ import annotations

import numpy as np

from illumination.water_column import model

SEAWATER_TABLE_SALINITY = 38.4      # g/kg, model.py's table scattering (Zhang et al. 2009 at 20 C)
SALT_SCATTERING_AT_TABLE = 0.30     # seawater at 38.4 g/kg over pure water (Morel 1966, 1968; Zhang et al. 2009)
SEAS_SALINITY = 0.14                # g/kg, the Open Moon's held fresh seas (biosphere and sea-chemistry design)
# Pure water's absorption (per metre), IOCCG Protocol Series vol. 1 (2018), Table 1.1: Jonasz and Fournier (2007)
# 250-340 nm, Morel et al. (2007) 340-415 nm; Pope and Fry (1997) continues from 420 nm in model.py's table.
PURE_WATER_UV = np.array([
    (250, 0.045), (255, 0.0392), (260, 0.0344), (265, 0.0303), (270, 0.0269), (275, 0.0240), (280, 0.0216),
    (285, 0.0194), (290, 0.0176), (295, 0.0160), (300, 0.0147), (305, 0.0134), (310, 0.0124), (315, 0.0114),
    (320, 0.0106), (325, 0.0098), (330, 0.0092), (335, 0.0085), (340, 0.0080), (345, 0.0075), (350, 0.0071),
    (355, 0.0068), (360, 0.0066), (365, 0.0063), (370, 0.0060), (375, 0.0056), (380, 0.0052), (385, 0.0050),
    (390, 0.0048), (395, 0.0047), (400, 0.0046), (405, 0.0046), (410, 0.0046), (415, 0.0046)])
JOIN_NM = 420.0


def scattering_factor(salinity_g_kg):
    """Pure water's molecular scattering at a salinity over model.py's seawater at 38.4 g/kg."""
    return (1 + SALT_SCATTERING_AT_TABLE * salinity_g_kg / SEAWATER_TABLE_SALINITY) / (1 + SALT_SCATTERING_AT_TABLE)


def pure_water(wavelength_nm, salinity_g_kg=SEAS_SALINITY, uv_absorption=True):
    """Absorption and scattering (per metre) of pure water at a salinity. uv_absorption=False sets the absorption
    below 420 nm to zero, the bound beyond every measurement."""
    w = np.asarray(wavelength_nm, float)
    a_table, b_table = model.seawater(w)
    if uv_absorption:
        a_uv = np.interp(w, PURE_WATER_UV[:, 0], PURE_WATER_UV[:, 1])
        # From 415 to 420 nm the composite runs straight to Pope and Fry's 420 nm value, 0.00006 per metre below it.
        top = float(model.seawater(JOIN_NM)[0])
        a_uv = np.where(w > PURE_WATER_UV[-1, 0], PURE_WATER_UV[-1, 1] + (top - PURE_WATER_UV[-1, 1])
                        * (w - PURE_WATER_UV[-1, 0]) / (JOIN_NM - PURE_WATER_UV[-1, 0]), a_uv)
    else:
        a_uv = np.zeros_like(w)
    return np.where(w < JOIN_NM, a_uv, a_table), b_table * scattering_factor(salinity_g_kg)


def water(wavelength_nm, chlorophyll_mg_m3=0.0, cdom_440=0.0, fines_g_m3=0.0, soil="12001",
          fines_diameter_um=4.0, salinity_g_kg=SEAS_SALINITY):
    """model.water() with fresh pure water: total and per-constituent absorption, scattering and backscattering."""
    w = np.asarray(wavelength_nm, float)
    a_w, b_w = pure_water(w, salinity_g_kg)
    a_p, b_p, bb_p = model.phytoplankton(w, chlorophyll_mg_m3)
    a_g = model.dissolved_organic(w, cdom_440)
    if fines_g_m3 > 0:
        a_f, b_f, bb_f = model.fines(w, fines_g_m3, soil, fines_diameter_um)
    else:
        a_f = b_f = bb_f = np.zeros_like(w)
    parts = dict(pure_water=(a_w, b_w, 0.5 * b_w), phytoplankton=(a_p, b_p, bb_p),
                 dissolved_organic=(a_g, np.zeros_like(w), np.zeros_like(w)), fines=(a_f, b_f, bb_f))
    return dict(a=sum(p[0] for p in parts.values()), b=sum(p[1] for p in parts.values()),
                bb=sum(p[2] for p in parts.values()), parts=parts)
