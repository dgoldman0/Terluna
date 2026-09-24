"""Ultraviolet and visible sunlight in a column: photolysis rates and solar heating.

Sunlight from 202 to 850 nm in 1-nm bins (TSIS-1 HSRS through the radiative-
convective model's solar spectrum, optionally filtered by a spectral shield) is
followed through the column with the delta-two-stream solver of
atmosphere/radiative_convective/shortwave.py, including its pseudo-spherical
direct beam. Absorbers are O2 (Herzberg continuum), O3, NO2, NO3, N2O5, N2O,
HNO3, H2O2 and the O2-O2 and O2-N2 collision-induced bands; air scatters by
Rayleigh's law (Bodhaine et al. 1999). Cross-sections and quantum yields are the
JPL-recommended laboratory values listed in inputs.json.

Photolysis rates use the actinic flux: the direct beam averaged exactly over each
layer's exponential attenuation, plus twice the diffuse hemispheric fluxes
(isotropic diffuse light, the usual two-stream convention). Heating is the
divergence of the net flux, all absorbed energy taken as local heat. Global
means average the sunlit hemisphere over Gauss nodes in the solar zenith cosine
and halve the result for the night side.

Not included: sunlight below 202 nm (the TSIS-1 spectrum starts there, and every
shield considered blocks it), H2O and HO2 photolysis, temperature dependence of
the HNO3, H2O2, N2O and N2O5 cross-sections, O(1D) from O3 beyond 340 nm, and
pressure-dependent O2 Herzberg absorption.
"""
from __future__ import annotations
import math
import re
import numpy as np

from shared.constants import PLANCK, SPEED_OF_LIGHT
from atmosphere.radiative_convective import shortwave as sw, spectroscopy as sp
from atmosphere.radiative_convective.thermodynamics import layer_number_density
from atmosphere.middle_atmosphere import fetch_inputs

EDGES_NM = np.arange(202.0, 851.0, 1.0)
CENTRES_NM = 0.5 * (EDGES_NM[1:] + EDGES_NM[:-1])
HEATING_LIMIT_NM = 500.0      # heating is taken from this module below 500 nm, from the line-by-line model above
ABSORBERS = ('O2', 'O3', 'NO2', 'NO3', 'N2O5', 'N2O', 'HNO3', 'H2O2')


def _rows(name):
    """Numeric rows of an atlas file: (lambda, value) points or (lo, hi, value) intervals."""
    out = []
    for line in fetch_inputs.path(name).read_text(encoding='latin-1').splitlines():
        if '(' in line or not line.strip():
            continue
        m = re.match(r'^\s*([\d.]+)\s*-\s*([\d.]+)\s+([\d.eE+-]+)\s*$', line)
        if m:
            out.append((float(m.group(1)), float(m.group(2)), float(m.group(3))))
            continue
        parts = line.split()
        out.append((float(parts[0]), float(parts[1])))
    return out


def _on_bins(name, fill_below=0.0, fill_above=0.0):
    """Tabulated data evaluated at the bin centres: intervals piecewise constant, points linear."""
    rows = _rows(name)
    values = np.full(CENTRES_NM.size, np.nan)
    if len(rows[0]) == 3:
        for lo, hi, v in rows:
            values[(CENTRES_NM >= lo) & (CENTRES_NM < hi)] = v
        lo, hi = rows[0][0], rows[-1][1]
    else:
        x = np.array([r[0] for r in rows]); y = np.array([r[1] for r in rows])
        inside = (CENTRES_NM >= x[0]) & (CENTRES_NM <= x[-1])
        values[inside] = np.interp(CENTRES_NM[inside], x, y)
        lo, hi = x[0], x[-1]
    values = np.where(np.isnan(values) & (CENTRES_NM < lo), fill_below, values)
    values = np.where(np.isnan(values) & (CENTRES_NM > hi), fill_above, values)
    return np.nan_to_num(values, nan=0.0)


class CrossSections:
    """Absorption cross-sections (cm^2) and photolysis channels on the 1-nm bins."""

    O1D_TEMPS = (203, 223, 253, 273, 298, 321)

    def __init__(self):
        room = _on_bins('o3_jpl2010_293-298K_121.6-827.5nm.txt')
        cold = _on_bins('o3_jpl2010_218K_196-343nm.txt')
        self.o3_room, self.o3_cold = room, np.where(cold > 0, cold, room)
        self.no2_cold = _on_bins('no2_jpl2010_220K_240-662.5nm.txt')
        self.no2_warm = _on_bins('no2_jpl2010_294K_240-662.5nm.txt')
        # The Herzberg continuum table starts at 205 nm; 202-205 nm holds its first value.
        o2 = _rows('o2_jpl2010_298K_205-245nm.txt')
        self.fixed = dict(
            O2=_on_bins('o2_jpl2010_298K_205-245nm.txt', fill_below=o2[0][1]),
            NO3=_on_bins('no3_jpl2010_298K_403-691nm.txt'),
            N2O5=_on_bins('n2o5_jpl2010_195-300K_200-420nm.txt'),
            N2O=_on_bins('n2o_jpl2010_298K_160-240nm.txt'),
            HNO3=_on_bins('hno3_jpl2010_298K_192-350nm.txt'),
            H2O2=_on_bins('h2o2_jpl2010_298K_190-350nm.txt'))
        self.o1d = np.array([_on_bins(f'o1d_yield_matsumi2002_{t}K.txt', fill_below=0.9) for t in self.O1D_TEMPS])
        self.no2_yield = np.array([_on_bins(f'no2_yield_troe2000_{t}K.txt', fill_below=1.0) for t in (248, 298)])
        self.no3_no2 = np.array([_on_bins(f'no3_no2+o_yield_johnston1996_{t}K.txt', fill_below=1.0) for t in (190, 230, 298)])
        self.no3_no = np.array([_on_bins(f'no3_no+o2_yield_johnston1996_{t}K.txt') for t in (190, 230, 298)])

    @staticmethod
    def _between(t, temps, tables):
        """Linear interpolation in temperature between tabulated sets, clamped at the ends."""
        t = np.asarray(t, dtype=float)[..., None]
        temps = np.asarray(temps, dtype=float)
        tc = np.clip(t, temps[0], temps[-1])
        i = np.clip(np.searchsorted(temps, tc[..., 0], side='right') - 1, 0, len(temps) - 2)
        w = (tc - temps[i][..., None]) / (temps[i + 1] - temps[i])[..., None]
        return (1 - w) * tables[i] + w * tables[i + 1]

    def sigma(self, species, t):
        """Cross-sections (layers, bins) for layer temperatures t."""
        t = np.asarray(t, dtype=float)
        if species == 'O3':
            return self._between(t, (218.0, 295.0), np.array([self.o3_cold, self.o3_room]))
        if species == 'NO2':
            return self._between(t, (220.0, 294.0), np.array([self.no2_cold, self.no2_warm]))
        return np.broadcast_to(self.fixed[species], t.shape + (CENTRES_NM.size,))

    def channels(self, t):
        """Photolysis channels: name -> (absorber, cross-section times quantum yield (layers, bins))."""
        o3 = self.sigma('O3', t)
        phi = self._between(t, self.O1D_TEMPS, self.o1d)
        no3 = self.sigma('NO3', t)
        no2 = self.sigma('NO2', t)
        return {
            'O2 -> O + O': ('O2', self.sigma('O2', t)),
            'O3 -> O1D + O2': ('O3', o3 * phi),
            'O3 -> O + O2': ('O3', o3 * (1 - phi)),
            'NO2 -> NO + O': ('NO2', no2 * self._between(t, (248.0, 298.0), self.no2_yield)),
            'NO3 -> NO2 + O': ('NO3', no3 * self._between(t, (190.0, 230.0, 298.0), self.no3_no2)),
            'NO3 -> NO + O2': ('NO3', no3 * self._between(t, (190.0, 230.0, 298.0), self.no3_no)),
            'N2O5 -> NO2 + NO3': ('N2O5', self.sigma('N2O5', t)),
            'N2O -> N2 + O1D': ('N2O', self.sigma('N2O', t)),
            'HNO3 -> OH + NO2': ('HNO3', self.sigma('HNO3', t)),
            'H2O2 -> OH + OH': ('H2O2', self.sigma('H2O2', t)),
        }


def solar_bins(shield=None, step_nm=0.02):
    """Top-of-atmosphere irradiance per bin normal to the beam: energy (W/m^2) and photons (cm^-2 s^-1)."""
    lam = np.arange(EDGES_NM[0], EDGES_NM[-1] + step_nm / 2, step_nm)
    nu = 1e7 / lam
    per_wn, _ = sw.solar_spectrum(nu[::-1], shield=shield)
    per_nm = per_wn[::-1] * nu**2 / 1e7
    photons = per_nm * lam * 1e-9 / (PLANCK * SPEED_OF_LIGHT) * 1e-4
    energy = np.empty(CENTRES_NM.size); count = np.empty(CENTRES_NM.size)
    for i in range(CENTRES_NM.size):
        sl = slice(int(round(i / step_nm)), int(round((i + 1) / step_nm)) + 1)
        energy[i] = np.trapezoid(per_nm[sl], lam[sl])
        count[i] = np.trapezoid(photons[sl], lam[sl])
    return energy, count


class Sunlight:
    """Photolysis rates and ultraviolet-visible heating of a column."""

    def __init__(self, shield=None, n_mu0=6):
        self.xs = CrossSections()
        self.energy, self.photons = solar_bins(shield)
        self.mu0, self.w0 = sw.disk_nodes(n_mu0)
        nu = 1e7 / CENTRES_NM
        self.nu = nu
        self.cia = [(a, b, sp.CollisionInduced(f)) for a, b, f in sp.CIA_PAIRS
                    if 'O2' in (a, b) and 'H2O' not in (a, b)]
        self._fixed = (None, None)

    def _fixed_optics(self, col):
        """O2, collision-induced and Rayleigh optical depths, which depend only on the column (cached)."""
        if self._fixed[0] is col:
            return self._fixed[1]
        t = col['layer_t_k']
        dz = col['layer_dz_cm']
        n_tot = layer_number_density(col)
        o2 = col['dry_fractions']['O2'] * (1 - col['layer_x_h2o']) * n_tot * dz
        tau_abs = self.xs.sigma('O2', t) * o2[:, None]
        for k in range(t.size):
            x = col['layer_x_h2o'][k]
            frac = {s: v * (1 - x) for s, v in col['dry_fractions'].items()}
            for a, b, table in self.cia:
                tau_abs[k] += table.coefficient(self.nu, float(t[k])) * frac[a] * frac[b] * n_tot[k] ** 2 * dz[k]
        co2 = col['dry_fractions'].get('CO2', 0.0)
        weight = col['layer_column_cm2'] * (1 - 0.25 * col['layer_x_h2o'])
        tau_sca = weight[:, None] * sp.rayleigh_cross_section(self.nu, co2)[None, :]
        self._fixed = (col, (tau_abs, tau_sca, o2))
        return self._fixed[1]

    def optics(self, col, densities):
        """Absorption and Rayleigh optical depths (layers surface-first, bins) and per-absorber columns."""
        fixed_abs, tau_sca, o2 = self._fixed_optics(col)
        t = col['layer_t_k']
        dz = col['layer_dz_cm']
        tau_abs = fixed_abs.copy()
        columns = {'O2': o2}
        for s in ABSORBERS[1:]:
            n = densities.get(s)
            if n is None:
                continue
            columns[s] = n * dz
            tau_abs += self.xs.sigma(s, t) * columns[s][:, None]
        return tau_abs, tau_sca, columns

    def evaluate(self, col, densities, surface_albedo, spherical=True):
        """Global-mean photolysis rates (s^-1 per layer) and heating (W/m^2 per layer, below HEATING_LIMIT_NM).

        densities: species -> number density (cm^-3) at layers, surface first.
        """
        tau_abs, tau_sca, _ = self.optics(col, densities)
        top_down = (tau_abs[::-1], tau_sca[::-1], np.zeros_like(tau_sca))
        radius = (col['planet'].radius_m + col['z_m'])[::-1] if spherical else None
        nlay = tau_abs.shape[0]
        actinic = np.zeros((nlay, CENTRES_NM.size))      # photons cm^-2 s^-1 per unit TOA photon flux
        absorbed = np.zeros((nlay, CENTRES_NM.size))     # per unit TOA flux normal to the beam
        surface = np.zeros(CENTRES_NM.size)
        toa_net = np.zeros(CENTRES_NM.size)
        surface_net = np.zeros(CENTRES_NM.size)
        below = CENTRES_NM < HEATING_LIMIT_NM
        by_node = []                                      # local heating (W/m^2 per layer) at each mu0 node
        for m, w in zip(self.mu0, self.w0):
            direct, diffuse, up = sw.column_fluxes(*top_down, surface_albedo, m, radius)
            by_node.append((((direct + diffuse - up)[:-1] - (direct + diffuse - up)[1:])[:, below]
                            @ self.energy[below])[::-1])
            surface += 0.5 * w * (direct[-1] + diffuse[-1]) * self.energy
            toa_net += 0.5 * w * (direct[0] - up[0]) * self.energy
            surface_net += 0.5 * w * (direct[-1] + diffuse[-1] - up[-1]) * self.energy
            slant = -np.log(np.maximum(direct[1:] / np.maximum(direct[:-1], 1e-300), 1e-300))
            beam_top = direct[:-1] / m
            mean_beam = beam_top * np.where(slant > 1e-6, -np.expm1(-slant) / np.maximum(slant, 1e-300), 1 - slant / 2)
            diffuse_mean = (diffuse[:-1] + diffuse[1:] + up[:-1] + up[1:])
            layer_actinic = mean_beam + diffuse_mean            # 2 x mean of (down + up) diffuse
            net = direct + diffuse - up
            actinic += 0.5 * w * layer_actinic[::-1]
            absorbed += 0.5 * w * (net[:-1] - net[1:])[::-1]
        rates = {}
        for name, (absorber, sigma_phi) in self.xs.channels(col['layer_t_k']).items():
            rates[name] = (sigma_phi * actinic * self.photons[None, :]).sum(axis=1)
        heating_bins = absorbed * self.energy[None, :]
        return dict(rates=rates, heating=heating_bins[:, below].sum(axis=1), heating_nodes=np.array(by_node).T,
                    heating_all=heating_bins.sum(axis=1), actinic=actinic, surface_spectrum=surface,
                    toa_net=toa_net, surface_net=surface_net)

    def noon_surface(self, col, densities, albedo, mu0, spherical=True):
        """Downward surface irradiance per bin (W/m^2) for one solar zenith cosine."""
        tau_abs, tau_sca, _ = self.optics(col, densities)
        radius = (col['planet'].radius_m + col['z_m'])[::-1] if spherical else None
        direct, diffuse, _ = sw.column_fluxes(tau_abs[::-1], tau_sca[::-1], np.zeros_like(tau_sca), albedo, mu0, radius)
        return (direct[-1] + diffuse[-1]) * self.energy


def erythemal_weight(lam_nm):
    """CIE erythemal action spectrum (McKinlay and Diffey 1987; ISO 17166)."""
    lam = np.asarray(lam_nm, dtype=float)
    return np.where(lam <= 298, 1.0, np.where(lam <= 328, 10 ** (0.094 * (298 - lam)),
                    np.where(lam <= 400, 10 ** (0.015 * (140 - lam)), 0.0)))


def uv_index(surface_per_bin):
    """UV index of a surface irradiance per 1-nm bin: 40 m^2/W times the erythemally weighted irradiance."""
    return 40.0 * float(np.sum(surface_per_bin * erythemal_weight(CENTRES_NM)))
