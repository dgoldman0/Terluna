"""Steady-state photochemistry of oxygen, hydrogen and nitrogen species in a column.

Fourteen species are solved at layer mid-points: O, O(1D), O3, H, OH, HO2,
H2O2, H2, NO, NO2, NO3, N2O5, HNO3 and N2O. N2, O2, water vapour and the air
density are fixed by the column. Rate coefficients are the NASA/JPL evaluation
19-5 (Burkholder et al. 2019) except where a reaction lists another source.
Species move by eddy diffusion of their mixing ratios through spherical shells
(the Moon's middle atmosphere reaches several hundred kilometres, where shell
area matters), with deposition and fixed-mixing-ratio conditions at the ground,
first-order washout of soluble species in the troposphere, and no flux through
the top. The steady state is found by implicit (backward-Euler) time steps of
growing length with Newton iterations, the photolysis rates being refreshed from
the current composition as the solution evolves.

Left out, and stated with every result: chlorine and bromine chemistry,
methane, CO and hydrocarbon chemistry (so tropospheric ozone production),
HO2NO2, heterogeneous reactions on aerosols, ions, molecular diffusion and
escape through the top, and any photolysis below 202 nm.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from shared.constants import BOLTZMANN
from atmosphere.radiative_convective.thermodynamics import layer_number_density

SPECIES = ('O', 'O1D', 'O3', 'H', 'OH', 'HO2', 'H2O2', 'H2', 'NO', 'NO2', 'NO3', 'N2O5', 'HNO3', 'N2O')
FIXED = ('O2', 'N2', 'H2O')
INDEX = {s: i for i, s in enumerate(SPECIES)}
DU = 2.6867e16   # molecules cm^-2 per Dobson unit


def arrhenius(a, e_over_r):
    return lambda t, m, h2o: a * np.exp(-e_over_r / t)


def termolecular(k0, n, kinf=None, minf=0.0):
    """JPL falloff form with Fc = 0.6; k0(T) = k0 (T/300)^-n, kinf(T) = kinf (T/300)^-m."""
    def k(t, m, h2o):
        low = k0 * (t / 300.0) ** (-n) * m
        if kinf is None:
            return low
        high = kinf * (t / 300.0) ** (-minf)
        ratio = low / high
        return low / (1 + ratio) * 0.6 ** (1.0 / (1 + np.log10(ratio) ** 2))
    return k


def _ho2_ho2(t, m, h2o):
    return (3.0e-13 * np.exp(460.0 / t) + 2.1e-33 * m * np.exp(920.0 / t)) * (1 + 1.4e-21 * h2o * np.exp(2200.0 / t))


def _oh_hno3(t, m, h2o):
    k0 = 2.4e-14 * np.exp(460.0 / t)
    k2 = 2.7e-17 * np.exp(2199.0 / t)
    k3 = 6.5e-34 * np.exp(1335.0 / t) * m
    return k0 + k3 / (1 + k3 / k2)


_NO2_NO3 = termolecular(2.4e-30, 3.0, 1.6e-12, -0.1)


def _n2o5_decomposition(t, m, h2o):
    return _NO2_NO3(t, m, h2o) / (2.7e-27 * np.exp(11000.0 / t))


# (reactants, products, rate(T, [M], [H2O]) in cm^3 s^-1 per the listed reactants, source)
REACTIONS = [
    (('O', 'O2'), ('O3',), termolecular(6.1e-34, 2.4), 'JPL 19-5'),
    (('O', 'O3'), ('O2', 'O2'), arrhenius(8.0e-12, 2060.0), 'JPL 19-5'),
    (('O', 'O'), ('O2',), lambda t, m, h2o: 4.7e-33 * (300.0 / t) ** 2 * m, 'Campbell and Gray 1973'),
    (('O1D', 'N2'), ('O', 'N2'), arrhenius(2.15e-11, -110.0), 'JPL 19-5'),
    (('O1D', 'O2'), ('O', 'O2'), arrhenius(3.3e-11, -55.0), 'JPL 19-5'),
    (('O1D', 'O3'), ('O2', 'O2'), arrhenius(1.2e-10, 0.0), 'JPL 19-5'),
    (('O1D', 'O3'), ('O2', 'O', 'O'), arrhenius(1.2e-10, 0.0), 'JPL 19-5'),
    (('O1D', 'H2O'), ('OH', 'OH'), arrhenius(1.63e-10, -60.0), 'JPL 19-5'),
    (('O1D', 'H2'), ('OH', 'H'), arrhenius(1.2e-10, 0.0), 'JPL 19-5'),
    (('O1D', 'N2O'), ('NO', 'NO'), arrhenius(7.25e-11, -20.0), 'JPL 19-5'),
    (('O1D', 'N2O'), ('N2', 'O2'), arrhenius(4.63e-11, -20.0), 'JPL 19-5'),
    (('H', 'O2'), ('HO2',), termolecular(5.3e-32, 1.8, 9.5e-11, -0.4), 'JPL 19-5'),
    (('H', 'O3'), ('OH', 'O2'), arrhenius(1.4e-10, 470.0), 'JPL 19-5'),
    (('H', 'HO2'), ('OH', 'OH'), arrhenius(7.2e-11, 0.0), 'JPL 19-5'),
    (('H', 'HO2'), ('H2', 'O2'), arrhenius(6.9e-12, 0.0), 'JPL 19-5'),
    (('H', 'HO2'), ('H2O', 'O'), arrhenius(1.6e-12, 0.0), 'JPL 19-5'),
    (('OH', 'O'), ('O2', 'H'), arrhenius(1.8e-11, -180.0), 'JPL 19-5'),
    (('HO2', 'O'), ('OH', 'O2'), arrhenius(3.0e-11, -200.0), 'JPL 19-5'),
    (('OH', 'O3'), ('HO2', 'O2'), arrhenius(1.7e-12, 940.0), 'JPL 19-5'),
    (('HO2', 'O3'), ('OH', 'O2', 'O2'), arrhenius(1.0e-14, 490.0), 'JPL 19-5'),
    (('OH', 'HO2'), ('H2O', 'O2'), arrhenius(4.8e-11, -250.0), 'JPL 19-5'),
    (('HO2', 'HO2'), ('H2O2', 'O2'), _ho2_ho2, 'JPL 19-5, with the water-vapour enhancement'),
    (('OH', 'OH'), ('H2O', 'O'), arrhenius(1.8e-12, 0.0), 'JPL 19-5'),
    (('OH', 'OH'), ('H2O2',), termolecular(6.9e-31, 1.0, 2.6e-11, 0.0), 'JPL 19-5'),
    (('OH', 'H2O2'), ('H2O', 'HO2'), arrhenius(1.8e-12, 0.0), 'JPL 19-5'),
    (('O', 'H2O2'), ('OH', 'HO2'), arrhenius(1.4e-12, 2000.0), 'JPL 19-5'),
    (('OH', 'H2'), ('H2O', 'H'), arrhenius(2.8e-12, 1800.0), 'JPL 19-5'),
    (('NO', 'O3'), ('NO2', 'O2'), arrhenius(3.0e-12, 1500.0), 'JPL 19-5'),
    (('NO2', 'O'), ('NO', 'O2'), arrhenius(5.1e-12, -210.0), 'JPL 19-5'),
    (('NO', 'HO2'), ('NO2', 'OH'), arrhenius(3.44e-12, -260.0), 'JPL 19-5'),
    (('NO', 'O'), ('NO2',), termolecular(9.1e-32, 1.5, 3.0e-11, 0.0), 'JPL 19-5'),
    (('NO2', 'O'), ('NO3',), termolecular(3.4e-31, 1.6, 2.3e-11, -0.2), 'JPL 19-5'),
    (('NO2', 'OH'), ('HNO3',), termolecular(1.8e-30, 3.0, 2.8e-11, 0.0), 'JPL 19-5'),
    (('HNO3', 'OH'), ('H2O', 'NO3'), _oh_hno3, 'JPL 19-5'),
    (('NO2', 'O3'), ('NO3', 'O2'), arrhenius(1.2e-13, 2450.0), 'JPL 19-5'),
    (('NO3', 'NO'), ('NO2', 'NO2'), arrhenius(1.5e-11, -170.0), 'JPL 19-5'),
    (('NO3', 'O'), ('NO2', 'O2'), arrhenius(1.0e-11, 0.0), 'JPL 19-5'),
    (('NO3', 'OH'), ('NO2', 'HO2'), arrhenius(2.2e-11, 0.0), 'JPL 19-5'),
    (('NO3', 'HO2'), ('OH', 'NO2', 'O2'), arrhenius(3.5e-12, 0.0), 'JPL 19-5'),
    (('NO2', 'NO3'), ('N2O5',), _NO2_NO3, 'JPL 19-5'),
    (('N2O5',), ('NO2', 'NO3'), _n2o5_decomposition, 'JPL 19-5 forward rate and equilibrium constant'),
    (('H', 'NO2'), ('OH', 'NO'), arrhenius(4.0e-10, 340.0), 'JPL 19-5'),
]

# Photolysis channels (names as in photolysis.Sunlight) -> (absorber, products)
PHOTOLYSIS = {
    'O2 -> O + O': ('O2', ('O', 'O')),
    'O3 -> O1D + O2': ('O3', ('O1D', 'O2')),
    'O3 -> O + O2': ('O3', ('O', 'O2')),
    'NO2 -> NO + O': ('NO2', ('NO', 'O')),
    'NO3 -> NO2 + O': ('NO3', ('NO2', 'O')),
    'NO3 -> NO + O2': ('NO3', ('NO', 'O2')),
    'N2O5 -> NO2 + NO3': ('N2O5', ('NO2', 'NO3')),
    'N2O -> N2 + O1D': ('N2O', ('N2', 'O1D')),
    'HNO3 -> OH + NO2': ('HNO3', ('OH', 'NO2')),
    'H2O2 -> OH + OH': ('H2O2', ('OH', 'OH')),
}


@dataclass(frozen=True)
class Boundary:
    """Lower boundary: fixed surface mixing ratios, deposition velocities (cm/s) and washout.

    Washout removes soluble species below the tropopause at the rate at which
    precipitation replaces the column's water vapour (precipitation / water column).
    """
    fixed: tuple = (('N2O', 3.3e-7), ('H2', 5.3e-7))
    deposition: tuple = (('O3', 0.1), ('HNO3', 1.0), ('H2O2', 0.5), ('N2O5', 1.0), ('NO2', 0.1))
    soluble: tuple = ('HNO3', 'H2O2', 'N2O5')
    precipitation_m_per_year: float = 1.0


@dataclass(frozen=True)
class Mixing:
    """Eddy diffusion (cm^2/s): K_troposphere below the tropopause; above it
    K_min (p_tropopause / p)^0.5, the n^-1/2 growth of Massie and Hunten (1981),
    capped at K_max. `scale` multiplies the whole profile (for the Moon, 36 keeps
    the mixing time per scale height of the Earth values)."""
    troposphere: float = 1.0e5
    minimum: float = 3.0e3
    maximum: float = 1.0e7
    scale: float = 1.0

    def profile(self, p_levels, p_tropopause):
        k = np.where(p_levels >= p_tropopause, self.troposphere,
                     np.minimum(self.minimum * np.sqrt(p_tropopause / p_levels), self.maximum))
        return self.scale * k


def fixed_densities(col):
    n = layer_number_density(col)
    x = col['layer_x_h2o']
    dry = col['dry_fractions']
    return dict(M=n, O2=dry['O2'] * (1 - x) * n, N2=dry['N2'] * (1 - x) * n, H2O=x * n)


class Column:
    """Transport geometry and fixed composition of a column for the chemistry."""

    def __init__(self, col, p_tropopause, mixing: Mixing, boundary: Boundary):
        self.col = col
        self.nl = col['layer_p_pa'].size
        self.fixed = fixed_densities(col)
        self.t = col['layer_t_k']
        r = col['planet'].radius_m * 100.0 + col['z_m'] * 100.0          # cm, levels
        zc = col['layer_z_m'] * 100.0
        self.area = r**2                                                 # at levels (relative)
        self.volume = (r[1:]**3 - r[:-1]**3) / 3.0                        # per layer (same units)
        k_levels = mixing.profile(col['p_pa'], p_tropopause)
        n_levels = col['p_pa'] / (BOLTZMANN * col['t_k']) * 1e-6
        # Conductance between layer centres k and k+1 through level k+1.
        self.conductance = k_levels[1:-1] * n_levels[1:-1] / (zc[1:] - zc[:-1])
        self.boundary = boundary
        self.tropo = col['layer_p_pa'] >= p_tropopause
        water = np.sum(col['layer_mass_kg_m2'] * col['layer_x_h2o'] * 0.018015
                       / (0.018015 * col['layer_x_h2o'] + 0.02897 * (1 - col['layer_x_h2o'])))
        self.washout = (boundary.precipitation_m_per_year * 1000.0 / 3.156e7) / max(water, 1e-9)

    def transport(self):
        """Sparse coefficients: d n_k/dt contributions a[k, j] * n_j, the same for every species."""
        m = self.fixed['M']
        rows, cols, vals = [], [], []
        for k in range(self.nl - 1):
            # Flux upward through level k+1: -C (f_{k+1} - f_k), f = n / M.
            c = self.conductance[k] * self.area[k + 1]
            for (row, sign) in ((k, 1.0), (k + 1, -1.0)):
                vol = self.volume[row]
                rows += [row, row]
                cols += [k + 1, k]
                vals += [sign * c / (m[k + 1] * vol), -sign * c / (m[k] * vol)]
        return np.array(rows), np.array(cols), np.array(vals)


def _rates_and_jacobian(n, fixed, t, j_rates, include_jac=True):
    """Net production (species, layers) and its Jacobian blocks (layers, species, species)."""
    ns, nl = len(SPECIES), n.shape[1]
    conc = {s: n[i] for i, s in enumerate(SPECIES)}
    conc.update(fixed)
    dndt = np.zeros((ns, nl))
    jac = np.zeros((nl, ns, ns)) if include_jac else None
    m, h2o = fixed['M'], fixed['H2O']

    def apply(k, reactants, products):
        rate = k.copy()
        for s in reactants:
            rate = rate * conc[s]
        for s in reactants:
            if s in INDEX:
                dndt[INDEX[s]] -= rate
        for s in products:
            if s in INDEX:
                dndt[INDEX[s]] += rate
        if include_jac:
            for a, s in enumerate(reactants):
                if s not in INDEX:
                    continue
                d = k.copy()
                for b, s2 in enumerate(reactants):
                    if b != a:
                        d = d * conc[s2]
                col = INDEX[s]
                for s3 in reactants:
                    if s3 in INDEX:
                        jac[:, INDEX[s3], col] -= d
                for s3 in products:
                    if s3 in INDEX:
                        jac[:, INDEX[s3], col] += d

    for reactants, products, kfun, _ in REACTIONS:
        apply(np.broadcast_to(kfun(t, m, h2o), (nl,)).astype(float), reactants, products)
    for name, (absorber, products) in PHOTOLYSIS.items():
        apply(np.asarray(j_rates.get(name, np.zeros(nl)), dtype=float), (absorber,), products)
    return dndt, jac


@dataclass
class State:
    n: np.ndarray                     # (species, layers) cm^-3
    j_rates: dict = field(default_factory=dict)
    heating: np.ndarray = None        # W/m^2 absorbed per layer below the heating limit
    info: dict = field(default_factory=dict)


def initial_state(col, boundary: Boundary, ozone_ppm=1.0):
    m = layer_number_density(col)
    n = np.full((len(SPECIES), m.size), 1e-30) * m
    n[INDEX['O3']] = ozone_ppm * 1e-6 * m
    for s, f in boundary.fixed:
        n[INDEX[s]] = f * m
    return State(n=n)


def solve(col, sunlight, surface_albedo, p_tropopause, mixing=Mixing(), boundary=Boundary(), state=None,
          lightning_no_per_cm2_s=0.0, max_steps=4000, tol=1e-7, verbose=False):
    """Steady-state composition of the column. Returns a State."""
    geo = Column(col, p_tropopause, mixing, boundary)
    nl, ns = geo.nl, len(SPECIES)
    fixed = geo.fixed
    start = state or initial_state(col, boundary)
    n = np.maximum(start.n.copy(), 0.0)
    trow, tcol, tval = geo.transport()
    fixed_surface = {INDEX[s]: f for s, f in boundary.fixed}
    deposition = {INDEX[s]: v for s, v in boundary.deposition}
    soluble = [INDEX[s] for s in boundary.soluble]
    # Lightning NO, spread through the troposphere in proportion to air mass.
    lightning = np.zeros(nl)
    if lightning_no_per_cm2_s > 0:
        mass = np.where(geo.tropo, col['layer_mass_kg_m2'], 0.0)
        lightning = lightning_no_per_cm2_s * mass / mass.sum() * geo.area[0] / geo.volume

    def refresh(nn):
        dens = {s: nn[INDEX[s]] for s in ('O3', 'NO2', 'NO3', 'N2O5', 'N2O', 'HNO3', 'H2O2')}
        out = sunlight.evaluate(col, dens, surface_albedo)
        return out['rates'], out

    j_rates, light = refresh(n)

    def residual_and_matrix(nn, dt, old):
        dndt, jac = _rates_and_jacobian(nn, fixed, geo.t, j_rates)
        # Transport, deposition, washout and sources.
        for i in range(ns):
            dndt[i] += np.bincount(trow, weights=tval * nn[i][tcol], minlength=nl)
            if i in deposition:
                dndt[i, 0] -= deposition[i] * nn[i, 0] * geo.area[0] / geo.volume[0]
            if i in soluble:
                dndt[i] -= np.where(geo.tropo, geo.washout, 0.0) * nn[i]
        dndt[INDEX['NO']] += lightning
        res = (nn - old) / dt - dndt
        # Matrix of the Newton system, rows (layer, species).
        rows, cols, vals = [], [], []
        base = np.arange(nl) * ns
        for i in range(ns):
            for j in range(ns):
                v = -jac[:, i, j]
                if i == j:
                    v = v + 1.0 / dt
                    if i in deposition:
                        v = v.copy(); v[0] += deposition[i] * geo.area[0] / geo.volume[0]
                    if i in soluble:
                        v = v + np.where(geo.tropo, geo.washout, 0.0)
                nz = v != 0
                rows.append(base[nz] + i); cols.append(base[nz] + j); vals.append(v[nz])
            rows.append(trow * ns + i); cols.append(tcol * ns + i); vals.append(-tval)
        for i, f in fixed_surface.items():
            res[i, 0] = nn[i, 0] - f * fixed['M'][0]
        rows = np.concatenate(rows); cols = np.concatenate(cols); vals = np.concatenate(vals)
        keep = np.ones(rows.size, dtype=bool)
        for i in fixed_surface:
            keep &= rows != i
        rows = np.concatenate([rows[keep], [i for i in fixed_surface]])
        cols = np.concatenate([cols[keep], [i for i in fixed_surface]])
        vals = np.concatenate([vals[keep], np.ones(len(fixed_surface))])
        mat = coo_matrix((vals, (rows, cols)), shape=(nl * ns, nl * ns)).tocsc()
        return res.T.ravel(), mat

    dt, time_total, steps = (1e-4 if state is None else 1e2), 0.0, 0
    floor = 1e-30 * fixed['M']
    history = []
    while steps < max_steps:
        steps += 1
        old = n.copy()
        trial = n.copy()
        ok = False
        for _ in range(8):
            res, mat = residual_and_matrix(trial, dt, old)
            try:
                delta = splu(mat).solve(-res).reshape(nl, ns).T
            except RuntimeError:
                break
            trial = trial + delta
            if not np.all(np.isfinite(trial)):
                break
            trial = np.maximum(trial, floor)
            change = np.max(np.abs(delta) / np.maximum(np.abs(trial), 1e-6 * fixed['M'] * 1e-12 + 1.0))
            if change < 1e-3:
                ok = True
                break
        if not ok:
            dt /= 4.0
            if dt < 1e-10:
                raise RuntimeError('Chemistry step size collapsed')
            continue
        rel = np.max(np.abs(trial - old) / np.maximum(trial, 1e-10 * fixed['M']))
        n = trial
        time_total += dt
        history.append((dt, rel))
        if dt < 1e9 or steps % 5 == 0:
            j_rates, light = refresh(n)
        if verbose and steps % 20 == 0:
            print(f'step {steps} dt={dt:.2e} rel={rel:.2e} O3 column={column_du(col, n):.1f} DU', flush=True)
        if dt > 1e15 and rel < tol:
            j_rates, light = refresh(n)
            break
        dt = min(dt * (2.0 if rel < 0.1 else 1.2), 1e17)
    return State(n=n, j_rates=j_rates, heating=light['heating'],
                 info=dict(steps=steps, time_s=time_total, converged=steps < max_steps, light=light))


def column_du(col, n):
    return float(np.sum(n[INDEX['O3']] * col['layer_dz_cm']) / DU)
