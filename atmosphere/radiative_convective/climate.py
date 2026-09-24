"""Inverse climate calculation: outgoing longwave and absorbed sunlight against surface temperature.

For a prescribed column (thermodynamics.py) at surface temperature Ts, compute
the clear-sky outgoing longwave radiation (OLR) and the global-mean absorbed
solar radiation (ASR). Where ASR(Ts) = OLR(Ts) with OLR rising faster, the
column is a possible equilibrium; if ASR exceeds the largest OLR any Ts can
reach, the water-rich column has no equilibrium (the runaway-greenhouse limit).
This is the standard one-dimensional habitable-zone method (Kasting 1988;
Kopparapu et al. 2013). Clouds are not included here; their effect enters as
a separate, explicit assumption.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import hashlib
import json
import math
from pathlib import Path
import numpy as np

from shared.constants import SOLAR_CONSTANT
from atmosphere.radiative_convective import thermodynamics as th, optics, longwave as lw, shortwave as sw
from atmosphere.radiative_convective import spectroscopy as sp

UV_GRID = sp.Grid.span(20000.0, 49500.0, 2.0)   # 202-500 nm: Rayleigh and O2-O2 only
# Spectral-shield transmission, a data product of the protection domain.
SHIELD_PRODUCT = Path(__file__).resolve().parents[2] / 'protection' / 'spectra' / 'shield_transmission.json'
SHIELDS = ('none', 'titania_stack', 'edge_200nm', 'edge_310nm')


def load_shield(name):
    """Transmission function of wavelength (nm) for a named shield, or None for unfiltered sunlight."""
    if name in (None, 'none'):
        return None
    raw = SHIELD_PRODUCT.read_bytes()
    data = json.loads(raw)
    if data['schema'] != 'terluna.protection.shield-transmission/1':
        raise ValueError('Unexpected shield product schema')
    w = np.asarray(data['wavelength_nm'], dtype=float)
    t = np.asarray(data['transmission'][name], dtype=float)

    def transmission(lam_nm):
        return np.interp(lam_nm, w, t, left=t[0], right=t[-1])
    transmission.label = name
    transmission.product_sha256 = hashlib.sha256(raw).hexdigest()[:16]
    return transmission


@dataclass(frozen=True)
class Scenario:
    name: str
    planet: th.Planet
    dry_pressure_pa: float
    co2_ppm: float
    humidity: str = 'manabe_wetherald'      # or 'saturated', or 'uniform:<value>'
    stratosphere_k: float = 200.0
    surface_albedo: float = 0.13
    longwave_levels: tuple = (70, 20)       # troposphere (to 1% of ps) and upper levels
    shortwave_levels: tuple = (40, 8)
    n_mu0: int = 6

    def air(self):
        return th.earthlike_air(self.dry_pressure_pa, self.co2_ppm)

    def rh(self):
        if self.humidity == 'manabe_wetherald':
            return th.manabe_wetherald()
        if self.humidity == 'saturated':
            return th.uniform_rh(1.0)
        if self.humidity.startswith('uniform:'):
            return th.uniform_rh(float(self.humidity.split(':')[1]))
        raise ValueError(f'Unknown humidity profile {self.humidity}')

    def column(self, ts, levels):
        spec = th.ColumnSpec(self.planet, self.air(), ts, self.stratosphere_k, self.rh(),
                             troposphere_levels=levels[0], upper_levels=levels[1])
        return th.build_column(spec)


def longwave(col, cfg=optics.LONGWAVE, n_angles=4, diffusivity=None, bins=10.0):
    nu, tau = optics.optical_depth(col, cfg)
    up, down = lw.fluxes(nu, tau, col['t_k'], col['t_k'][0], n_angles=n_angles, diffusivity=diffusivity)
    f_up, f_down = lw.integrate(nu, up), lw.integrate(nu, down)
    # Spectral OLR in bins for diagnostics.
    edges = np.arange(nu[0], nu[-1] + bins, bins)
    idx = np.clip(((nu - nu[0]) // bins).astype(int), 0, len(edges) - 2)
    spectral = np.bincount(idx, weights=up[-1] * (nu[1] - nu[0]), minlength=len(edges) - 1)
    return dict(olr=float(f_up[-1]), surface_down=float(f_down[0]), surface_up=float(f_up[0]),
                net_up=f_up - f_down, olr_bins=(edges[:-1], spectral))


def _uv_optics(col, nu):
    """Absorption (O2-O2 and O2-N2 collision-induced) on the ultraviolet-visible grid, surface-first."""
    rows = []
    n_tot = th.layer_number_density(col)
    tables = [(a, b, sp.CollisionInduced(f)) for a, b, f in sp.CIA_PAIRS if 'O2' in (a, b) and 'H2O' not in (a, b)]
    for k in range(len(col['layer_p_pa'])):
        x = col['layer_x_h2o'][k]
        frac = {s: v * (1 - x) for s, v in col['dry_fractions'].items()}
        tau = np.zeros(nu.size)
        for a, b, table in tables:
            tau += table.coefficient(nu, col['layer_t_k'][k]) * frac[a] * frac[b] * n_tot[k] ** 2 * col['layer_dz_cm'][k]
        rows.append(tau)
    return np.vstack(rows)


def shortwave_response(col, surface_albedo, cfg=optics.SHORTWAVE, n_mu0=6, spherical=True):
    """Global-mean clear-sky reflection and surface absorption per unit solar spectral irradiance.

    Returned per spectral grid, so any sunlight spectrum (for example a shielded
    one) applies as a dot product without repeating the radiative transfer.
    """
    radius = (col['planet'].radius_m + col['z_m'])[::-1] if spherical else None
    nu, tau_abs = optics.optical_depth(col, cfg)
    parts = []
    for nu_k, tau_k in ((nu, tau_abs), (UV_GRID.nu, _uv_optics(col, UV_GRID.nu))):
        tau_r = optics.rayleigh_optical_depth(col, nu_k)
        res = sw.planetary_mean(nu_k, np.ones(nu_k.size), tau_k[::-1], tau_r[::-1], np.zeros_like(tau_r),
                                surface_albedo, radius, n_mu0=n_mu0)
        parts.append((nu_k, res['weights'], res['spectral_reflected'], res['spectral_surface']))
    return parts


def apply_sunlight(parts, shield=None):
    """Global-mean solar budget (W/m^2) of a response for sunlight filtered by `shield`."""
    reflected = surface = incident = 0.0
    longward = None
    for nu, weights, reflect_unit, surface_unit in parts:
        spectrum, lw_part = sw.solar_spectrum(nu, shield=shield)
        if longward is None:
            longward = lw_part
        reflected += float((reflect_unit * spectrum) @ weights)
        surface += float((surface_unit * spectrum) @ weights)
        incident += 0.25 * float(spectrum @ weights)
    # Sunlight longward of 5 um (below the grid) is taken as absorbed in the atmosphere.
    incident_total = incident + 0.25 * longward
    return dict(incident=incident_total, reflected=reflected, albedo=reflected / incident_total,
                absorbed=incident_total - reflected, surface_absorbed=surface,
                atmosphere_absorbed=incident_total - reflected - surface,
                longward_of_grid=0.25 * longward, grid_incident=incident)


def shortwave(col, surface_albedo, cfg=optics.SHORTWAVE, n_mu0=6, spherical=True, shield=None):
    """Global-mean clear-sky solar budget of a column (W/m^2)."""
    return apply_sunlight(shortwave_response(col, surface_albedo, cfg, n_mu0, spherical), shield)


def evaluate(scenario: Scenario, ts: float, solar=True, longwave_too=True, shields=('none',)):
    """Clear-sky energy budget of a scenario at surface temperature ts, one row per shield."""
    col_lw = scenario.column(ts, scenario.longwave_levels)
    base = dict(scenario=scenario.name, planet=scenario.planet.name, dry_pressure_pa=scenario.dry_pressure_pa,
                co2_ppm=scenario.co2_ppm, humidity=scenario.humidity, stratosphere_k=scenario.stratosphere_k,
                surface_albedo=scenario.surface_albedo, ts_k=ts,
                surface_pressure_pa=col_lw['surface_pressure_pa'], tropopause_pa=col_lw['tropopause_pa'],
                tropopause_height_km=float(np.interp(math.log(col_lw['tropopause_pa']), np.log(col_lw['p_pa'][::-1]),
                                                     col_lw['z_m'][::-1]) / 1e3),
                stratospheric_h2o=col_lw['tropopause_h2o'],
                water_column_kg_m2=float(np.sum(col_lw['layer_mass_kg_m2'] * col_lw['layer_x_h2o'] * th.MOLAR_MASS['H2O']
                                                / (th.MOLAR_MASS['H2O'] * col_lw['layer_x_h2o']
                                                   + scenario.air().molar_mass * (1 - col_lw['layer_x_h2o'])))),
                air_column_kg_m2=float(col_lw['layer_mass_kg_m2'].sum()))
    if longwave_too:
        lwr = longwave(col_lw)
        base.update(olr=lwr['olr'], surface_down_lw=lwr['surface_down'], surface_up_lw=lwr['surface_up'])
    rows = []
    parts = shortwave_response(scenario.column(ts, scenario.shortwave_levels), scenario.surface_albedo,
                               n_mu0=scenario.n_mu0) if solar else None
    for name in shields:
        row = dict(base, shield=name)
        if solar:
            shield = load_shield(name)
            swr = apply_sunlight(parts, shield)
            row.update(asr=swr['absorbed'], albedo=swr['albedo'], incident_sw=swr['incident'],
                       surface_sw=swr['surface_absorbed'], atmosphere_sw=swr['atmosphere_absorbed'],
                       sw_longward_of_grid=swr['longward_of_grid'],
                       shield_product=getattr(shield, 'product_sha256', ''))
            if 'olr' in row:
                row['net_toa'] = row['asr'] - row['olr']
        rows.append(row)
    return rows
