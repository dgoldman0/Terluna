"""The recovered rain solver of 11 October 2026, as far as the port now departs from it.

Since 25 October the port takes the vapour's diffusivity from Massman's (1998) measured fit (drops.vapour_diffusivity).
The recovered solver took rainprops' Chapman-Enskog value: the vapour's binary diffusion through each gas, with
Lennard-Jones parameters, combined by Blanc's law with the vapour's own self-diffusion term (Reid, Prausnitz and
Sherwood 1977), which stands 17.6% below Massman's at every temperature and pressure the columns reach. This module
keeps that value, so that the port can still be checked against the recovered solver's saved outputs, and loads the
solver as it stood on 24 October from its byte-for-byte copy, against which the products built before the fix are
checked (shared/baselines.py). No product reads it.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

from atmosphere.rain import drops

BASELINE = Path(__file__).resolve().parent / 'baseline' / 'drops_2026-10-24.py'


def chapman_enskog_diffusivity(temperature_k, pressure_pa, vapour_fraction, dry_fractions):
    """The vapour's diffusivity (m2/s) as rainprops and the recovered solver computed it."""
    fractions = dict(dry_fractions)
    fractions['N2'] = fractions.get('N2', 0.0) + fractions.pop('Ar', 0.0)
    names = [n for n in fractions if fractions[n] > 0] + ['H2O']
    x = np.array([fractions[n] * (1 - vapour_fraction) for n in names[:-1]] + [vapour_fraction])
    molar = np.array([drops.GASES[n][0] for n in names])
    sigma = np.array([drops.GASES[n][1] for n in names])
    well = np.array([drops.GASES[n][2] for n in names])
    pair_sigma = 0.5 * (sigma + drops.GASES['H2O'][1])
    t_star = temperature_k / np.sqrt(well * drops.GASES['H2O'][2])
    omega_d = (1.06036 / t_star**0.15610 + 0.193 / np.exp(0.47635 * t_star)
               + 1.03587 / np.exp(1.52996 * t_star) + 1.76474 / np.exp(3.89411 * t_star))
    binary = (1.858e-3 * np.sqrt(1 / (molar * 1e3) + 1 / (drops.WATER_MOLAR_MASS * 1e3)) / pair_sigma**2
              * temperature_k**1.5 / (pressure_pa / 101325.0) / omega_d * 1e-4)
    return float(1.0 / np.sum(x / binary))


class Column:
    """A drops.Column whose air carries the recovered solver's vapour diffusivity; every other property is the
    column's own."""

    def __init__(self, column):
        self.column = column

    def __getattr__(self, name):
        return getattr(self.column, name)

    def state(self, altitude_m=0.0):
        s = self.column.state(altitude_m)
        s['vapour_diffusivity_m2_s'] = chapman_enskog_diffusivity(
            s['temperature_k'], s['pressure_pa'], s['vapour_pressure_pa'] / s['pressure_pa'], self.column.air.fractions)
        return s


def baseline():
    """The rain solver as it stood on 24 October 2026, loaded from its byte-for-byte copy (baseline/, listed in
    research/provenance.json): the products built before the diffusivity fix are checked against it."""
    name = 'atmosphere.rain.baseline_drops_2026_10_24'
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, BASELINE)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]
