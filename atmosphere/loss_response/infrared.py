"""What the thermal column needs to radiate in the infrared: its coolants above the base and CO2's band escape.

The thermal column (atmosphere/thermal_column.py) radiates CO2's 15-um band, NO's 5.3-um band and atomic oxygen's
fine-structure lines in excess of what air at the base temperature emits, given each coolant's mole fraction against
log pressure above the 0.3 Pa base and the chance that a photon of CO2's band leaves through a column of CO2. This
module supplies both for one middle-atmosphere case:

- CO2 holds the middle atmosphere's 400 ppm at the base and above it follows the steady balance of eddy and molecular
  diffusion with no net flux, d ln x / d ln p = (m_CO2 / m - 1) D / (D + K): well mixed while the middle atmosphere's
  eddy mixing (scaled for the Moon) outruns molecular diffusion, separating toward its own scale height above the
  homopause. Its binary diffusion in air is Massman's (1998), 0.1381 cm^2/s at 273.15 K and 1 atm, rising as T^1.81
  at fixed pressure, at the base temperature throughout. Photolysis by the light through gaps and beyond the aperture
  is left out.
- Atomic oxygen and NO hold their mole fractions at the base, carried well mixed, as the limb tables carry the base's
  atomic oxygen; the atoms the oxygen step finds above the base are left out, so the cooling is a lower estimate.
- CO2's band escape is the strength-weighted mean over HITRAN's lines of the band, 500-850 cm^-1 (the bending
  fundamental, its hot bands and isotopologues, as the radiative model reads them), of a Doppler line's escape
  through the column at the base temperature, for photons emitted into the hemisphere toward it. The thermal column
  takes a level's escape as the mean of the ways up to space and down to the base, each column scaled for the
  Doppler width at its temperature.
"""
from __future__ import annotations
import csv
import functools
import json
import math
from pathlib import Path

import numpy as np

from shared.constants import AVOGADRO, BOLTZMANN, SPEED_OF_LIGHT
from atmosphere.thermal_column import doppler_escape_table
from atmosphere.radiative_convective import spectroscopy as sp, ck

ROOT = Path(__file__).resolve().parents[2]
MIDDLE = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results'
CO2_PPM = 400.0                          # the middle atmosphere's CO2 (its README; limb_heat.DESIGN)
BAND_CM = (500.0, 850.0)                 # CO2's 15-um band system
LOG_PRESSURE = np.round(np.arange(0.0, 40.0001, 0.25), 4)     # log pressures above the base for the coolants
ESCAPE_COLUMNS_CM2 = np.geomspace(1e10, 1e24, 85)
CO2_DIFFUSION = (0.1381, 273.15, 1.81)  # Massman (1998): D cm^2/s at T K and 1 atm, and its temperature exponent
CO2_MOLAR = 0.0440095
STANDARD_ATMOSPHERE_PA = 101325.0


@functools.lru_cache(maxsize=None)
def _band_lines():
    lines = sp.load_lines('CO2', ck.S_MIN['CO2']).select(*BAND_CM)
    return lines, sp.PartitionSums(lines.gid)


@functools.lru_cache(maxsize=None)
def band_escape(t_k):
    """CO2's band escape at temperature t_k: columns (cm^-2) and the share of a photon's chances to leave through
    each, for photons emitted into the hemisphere toward it."""
    lines, sums = _band_lines()
    c2 = sp.C2
    strength = (lines.strength * sums.ratio(lines.gid, t_k) * np.exp(-c2 * lines.elower * (1 / t_k - 1 / sp.T_REF))
                * -np.expm1(-c2 * lines.nu / t_k) / -np.expm1(-c2 * lines.nu / sp.T_REF))
    mass_kg = lines.mass * 1e-3 / AVOGADRO
    width = lines.nu * np.sqrt(2 * BOLTZMANN * t_k / mass_kg) / SPEED_OF_LIGHT          # 1/e half width, cm^-1
    sigma = strength / (math.sqrt(math.pi) * width)                                  # line centre, cm^2
    weight = strength / strength.sum()
    escape = np.array([weight @ doppler_escape_table(n * sigma) for n in ESCAPE_COLUMNS_CM2])
    return ESCAPE_COLUMNS_CM2.copy(), np.minimum(np.minimum.accumulate(escape), 1.0)


def base_mixing(case_name, species, base_pa):
    """A species' mole fraction at the base in the middle atmosphere's stored profile (log-interpolated)."""
    rows = list(csv.DictReader(open(MIDDLE / 'profiles' / f'{case_name}.csv', newline='')))
    p = np.array([float(r['p_pa']) for r in rows])
    x = np.array([float(r[species]) for r in rows])
    order = np.argsort(p)
    return float(np.exp(np.interp(math.log(base_pa), np.log(p[order]), np.log(np.maximum(x[order], 1e-30)))))


def co2_profile(case_name, base_pa, base_temperature_k, mean_molar, x=LOG_PRESSURE):
    """CO2's mole fraction at log pressures x above the base: well mixed below the homopause, separating above it."""
    from atmosphere.middle_atmosphere import run as middle_run
    case = middle_run.cases()[case_name]
    tropopause = json.loads((MIDDLE / 'cases' / f'{case_name}.json').read_text())['tropopause_pa']
    p = base_pa * np.exp(-np.asarray(x, float))
    eddy = case.mixing().profile(p, tropopause)                                       # cm^2/s
    d0, t0, power = CO2_DIFFUSION
    molecular = d0 * (base_temperature_k / t0) ** power * STANDARD_ATMOSPHERE_PA / p      # cm^2/s
    share = molecular / (molecular + eddy)
    rate = (CO2_MOLAR / mean_molar - 1.0) * share
    integral = np.r_[0.0, np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(x))]
    return CO2_PPM * 1e-6 * np.exp(-integral)


def column_fields(case_name, base_pa, base_temperature_k, mean_molar):
    """The thermal column's infrared-cooling fields for a middle-atmosphere case (ColumnConfig keywords)."""
    columns, escape = band_escape(round(float(base_temperature_k), 3))
    return dict(coolant_log_pressure=tuple(float(v) for v in LOG_PRESSURE),
                co2_mole_fraction=tuple(float(v) for v in co2_profile(case_name, base_pa, base_temperature_k,
                                                                      mean_molar)),
                o_mole_fraction=tuple([base_mixing(case_name, 'O', base_pa)] * LOG_PRESSURE.size),
                no_mole_fraction=tuple([base_mixing(case_name, 'NO', base_pa)] * LOG_PRESSURE.size),
                co2_escape_column_cm2=tuple(float(v) for v in columns), co2_escape=tuple(float(v) for v in escape))
