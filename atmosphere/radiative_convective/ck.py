"""Correlated-k thermal radiation built from the line-by-line spectroscopy.

Band k-distributions of H2O (lines and MT_CKD continuum), CO2, O3 and the
collision-induced absorption are computed once from line-by-line spectra on a
node grid in pressure, temperature and (for water) vapour fraction, and cached.
Each spectrum is resolved to a quarter of the Voigt half-width of a CO2 line at
that pressure (0.0002 to 0.01 cm^-1), so the low-pressure line cores that set
stratospheric and mesospheric cooling are sampled. At run time the gases in a
layer are combined by random overlap with resorting and rebinning (Amundsen et
al. 2017, A&A 598, A97) and the fluxes follow from the linear-in-optical-depth
source function of longwave.py, per band and g-point, with band-integrated
Planck radiances. Accuracy against the line-by-line model is recorded in the
tests and the validation record.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import numpy as np

from shared.constants import AVOGADRO, BOLTZMANN, SPEED_OF_LIGHT
from atmosphere.radiative_convective import spectroscopy as sp, longwave as lw
from atmosphere.radiative_convective.thermodynamics import layer_number_density, species_column

BANDS = ((1.0, 350.0), (350.0, 500.0), (500.0, 630.0), (630.0, 700.0), (700.0, 820.0), (820.0, 980.0),
         (980.0, 1080.0), (1080.0, 1180.0), (1180.0, 1390.0), (1390.0, 1480.0), (1480.0, 1800.0),
         (1800.0, 2080.0), (2080.0, 2250.0), (2250.0, 2380.0), (2380.0, 2600.0), (2600.0, 3000.0))
P_NODES = 10.0 ** (np.arange(-4, 23) / 4.0)        # 0.1 Pa to 3.2e5 Pa
T_NODES = np.arange(100.0, 400.0, 25.0)            # 100 to 375 K
X_NODES = np.array([0.0, 0.01, 0.05, 0.2, 0.5])    # water vapour mole fraction
S_MIN = {'H2O': 1e-28, 'CO2': 1e-28, 'O3': 1e-24}
TABLE_VERSION = 1
CACHE = sp.CACHE / 'ck'


def g_quadrature():
    """16 g-points: 8 Gauss points on [0, 0.9] and 8 on [0.9, 1]; weights sum to 1."""
    x, w = np.polynomial.legendre.leggauss(8)
    g1, w1 = 0.45 * (x + 1), 0.45 * w
    g2, w2 = 0.9 + 0.05 * (x + 1), 0.05 * w
    return np.concatenate([g1, g2]), np.concatenate([w1, w2])


G, W = g_quadrature()
NG = G.size
NB = len(BANDS)


def resolution(p_pa, t_k=200.0):
    """Spectral step for a node: a quarter of the Voigt HWHM of a CO2 line at 667 cm^-1."""
    alpha = 667.0 / SPEED_OF_LIGHT * math.sqrt(2 * math.log(2) * BOLTZMANN * t_k * AVOGADRO * 1e3 / 44.0)
    gamma = 0.07 * p_pa / 101325.0
    voigt = 0.5346 * gamma + math.sqrt(0.2166 * gamma**2 + alpha**2)
    return float(np.clip(0.25 * voigt, 2e-4, 0.01)), voigt


def _kdist(sigma):
    """k at the g-points of a band from its sampled cross-sections (inverse CDF)."""
    s = np.sort(sigma)
    pos = G * (s.size - 1)
    return np.interp(pos, np.arange(s.size), s)


_STATE = {}


def _load_state():
    if _STATE:
        return
    lines = {m: sp.load_lines(m, S_MIN[m]) for m in S_MIN}
    gids = np.concatenate([l.gid for l in lines.values()])
    _STATE.update(lines=lines, sums=sp.PartitionSums(gids), continuum=sp.WaterContinuum(),
                  cia={f'{a}-{b}': sp.CollisionInduced(f) for a, b, f in sp.CIA_PAIRS})


def table_key():
    _load_state()
    if 'key' in _STATE:
        return _STATE['key']
    digest = hashlib.sha256()
    for m in sorted(S_MIN):
        l = _STATE['lines'][m]
        for arr in (l.nu, l.strength, l.gamma_air, l.gamma_self, l.elower, l.n_air, l.delta_air):
            digest.update(np.ascontiguousarray(arr).tobytes())
    digest.update(repr((BANDS, P_NODES.tolist(), T_NODES.tolist(), X_NODES.tolist(), G.tolist(), S_MIN,
                        TABLE_VERSION, sp.KERNEL_VERSION, sp.VOIGT_CORE)).encode())
    _STATE['key'] = digest.hexdigest()[:16]
    return _STATE['key']


def _node(args):
    molecule, ip, it, ix = args
    path = CACHE / table_key() / f'{molecule}_p{ip:02d}_t{it:02d}_x{ix:d}.npy'
    if path.is_file():
        return path
    p, t = float(P_NODES[ip]), float(T_NODES[it])
    x = float(X_NODES[ix]) if molecule == 'H2O' else 0.0
    h, voigt = resolution(p, t)
    coarse = max(1, int(round(0.01 / h)))
    fine = min(1.0, max(40.0 * voigt, 0.03))
    out = np.zeros((NB, NG))
    for b, (lo, hi) in enumerate(BANDS):
        grid = sp.Grid.span(lo, hi, h)
        sigma = sp.line_absorption(grid, _STATE['lines'][molecule], _STATE['sums'], t, p, x * p,
                                   pedestal=(molecule == 'H2O'), coarse_factor=coarse, fine_halfwidth=fine)
        if molecule == 'H2O':
            s_c, f_c = _STATE['continuum'].cross_sections(grid.nu, t, p, x)
            sigma = sigma + s_c + f_c
        out[b] = _kdist(sigma)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix('.tmp.npy')
    np.save(tmp, out)
    tmp.replace(path)
    return path


def _path(m, ip, it, ix):
    return CACHE / table_key() / f'{m}_p{ip:02d}_t{it:02d}_x{ix:d}.npy'


def routine(m, ip, it, ix):
    """Nodes built in advance. Water fractions above 0.05 below 1e4 Pa, and temperatures of 150 K
    or less above 3e4 Pa, occur only in states far from these climates; such nodes are computed
    when a column first needs them (CKLongwave.ensure), never filled in."""
    if m == 'H2O' and X_NODES[ix] > 0.05 and P_NODES[ip] < 1e4:
        return False
    return not (T_NODES[it] <= 150.0 and P_NODES[ip] > 3e4)


def build(processes=0, everything=False, reverse=False):
    """Compute the missing routine nodes (or every node) in parallel; returns the table arrays."""
    _load_state()
    jobs = []
    for m in S_MIN:
        for ip in range(P_NODES.size):
            for it in range(T_NODES.size):
                for ix in (range(X_NODES.size) if m == 'H2O' else (0,)):
                    if (everything or routine(m, ip, it, ix)) and not _path(m, ip, it, ix).is_file():
                        jobs.append((m, ip, it, ix))
    procs = processes or max(1, (os.cpu_count() or 2) // 2)
    if jobs:
        jobs.sort(key=lambda j: j[1], reverse=reverse)
        with mp.get_context('fork').Pool(procs) as pool:
            for n, _ in enumerate(pool.imap_unordered(_node, jobs, chunksize=1), 1):
                if n % 50 == 0 or n == len(jobs):
                    print(f'ck table: {n}/{len(jobs)} nodes', flush=True)
    return load()


def load():
    """The table as built so far; nodes not yet computed are NaN (see CKLongwave.ensure)."""
    key = table_key()
    tables = {}
    for m in S_MIN:
        nx = X_NODES.size if m == 'H2O' else 1
        arr = np.full((P_NODES.size, T_NODES.size, nx, NB, NG), np.nan)
        for ip in range(P_NODES.size):
            for it in range(T_NODES.size):
                for ix in range(nx):
                    path = _path(m, ip, it, ix)
                    if path.is_file():
                        arr[ip, it, ix] = np.log(np.maximum(np.load(path), 1e-60))
                    elif routine(m, ip, it, ix):
                        raise FileNotFoundError(f'Correlated-k node {path.name} missing; run ck.build()')
        tables[m] = arr
    # Collision-induced absorption: temperature-only tables (cm^5/molecule^2), 1 cm^-1 sampling.
    cia = {}
    for name, table in _STATE['cia'].items():
        arr = np.empty((T_NODES.size, NB, NG))
        for it, t in enumerate(T_NODES):
            for b, (lo, hi) in enumerate(BANDS):
                nu = np.arange(lo, hi, 0.5)
                arr[it, b] = _kdist(table.coefficient(nu, float(t)))
        cia[name] = arr
    return dict(key=key, gas=tables, cia=cia)


def _weights(values, nodes, log=False):
    """Bracketing indices and linear weights for values on a node axis (clamped)."""
    v = np.log(values) if log else np.asarray(values, dtype=float)
    n = np.log(nodes) if log else nodes
    v = np.clip(v, n[0], n[-1])
    i = np.clip(np.searchsorted(n, v, side='right') - 1, 0, n.size - 2)
    w = (v - n[i]) / (n[i + 1] - n[i])
    return i, w


def gas_k(table, p, t, x=None):
    """Interpolated k (layers, bands, g) from a gas table, log-linear in p, T and x."""
    ip, wp = _weights(p, P_NODES, log=True)
    it, wt = _weights(t, T_NODES)
    if table.shape[2] > 1:
        ix, wx = _weights(x, X_NODES)
    else:
        ix, wx = np.zeros_like(ip), np.zeros_like(wp)
    out = 0.0
    for dp_, ap in ((0, 1 - wp), (1, wp)):
        for dt_, at in ((0, 1 - wt), (1, wt)):
            for dx_, ax in ((0, 1 - wx), (1, wx)):
                if table.shape[2] == 1 and dx_ == 1:
                    continue
                weight = (ap * at * (ax if table.shape[2] > 1 else 1.0))[:, None, None]
                out = out + weight * table[ip + dp_, it + dt_, np.minimum(ix + dx_, table.shape[2] - 1)]
    return np.exp(out)


def cia_k(table, t):
    it, wt = _weights(t, T_NODES)
    return (1 - wt)[:, None, None] * table[it] + wt[:, None, None] * table[it + 1]


def overlap(tau_a, tau_b):
    """Random overlap with resorting and rebinning of two (layers, bands, g) optical depths."""
    nl, nb, ng = tau_a.shape
    comb = (tau_a[..., :, None] + tau_b[..., None, :]).reshape(nl, nb, ng * ng)
    wts = np.broadcast_to((W[:, None] * W[None, :]).ravel(), comb.shape)
    order = np.argsort(comb, axis=-1)
    tau_s = np.take_along_axis(comb, order, axis=-1)
    w_s = np.take_along_axis(wts, order, axis=-1)
    cw = np.cumsum(w_s, axis=-1)
    cwt = np.cumsum(w_s * tau_s, axis=-1)
    edges = np.concatenate([[0.0], np.cumsum(W)])
    out = np.empty_like(tau_a)
    flat_cw = cw.reshape(-1, ng * ng)
    flat_cwt = cwt.reshape(-1, ng * ng)
    res = np.empty((flat_cw.shape[0], ng))
    for r in range(flat_cw.shape[0]):
        xs = np.concatenate([[0.0], flat_cw[r]])
        ys = np.concatenate([[0.0], flat_cwt[r]])
        cum = np.interp(edges, xs, ys)
        res[r] = np.diff(cum) / W
    return res.reshape(nl, nb, ng)


def band_planck(t):
    """Band-integrated Planck radiance (W m^-2 sr^-1) for temperatures t (any shape) -> (..., bands)."""
    t = np.asarray(t, dtype=float)
    out = np.empty(t.shape + (NB,))
    for b, (lo, hi) in enumerate(BANDS):
        nu = np.linspace(lo, hi, 400)
        vals = lw.planck(nu[None, :], t.reshape(-1, 1))
        out[..., b] = np.trapezoid(vals, nu, axis=-1).reshape(t.shape)
    return out


class CKLongwave:
    def __init__(self, tables=None):
        _load_state()
        self.tables = tables or load()

    def ensure(self, p, t, x):
        """Compute and load any table node the layers' interpolation brackets still lack."""
        ip, _ = _weights(p, P_NODES, log=True)
        it, _ = _weights(t, T_NODES)
        ix, _ = _weights(x, X_NODES)
        for m, table in self.tables['gas'].items():
            needed = set()
            for a, b, c in zip(ip, it, ix):
                for dp_ in (0, 1):
                    for dt_ in (0, 1):
                        for dx_ in ((0, 1) if m == 'H2O' else (None,)):
                            k = (a + dp_, b + dt_, 0 if dx_ is None else c + dx_)
                            if np.isnan(table[k][0, 0]):
                                needed.add(k)
            for k in sorted(needed):
                _node((m, *k))
                table[k] = np.log(np.maximum(np.load(_path(m, *k)), 1e-60))

    def optical_depth(self, col, extra=None):
        """Combined layer optical depth (layers, bands, g). `extra` maps molecule -> layer columns (cm^-2)."""
        p, t, x = col['layer_p_pa'], col['layer_t_k'], col['layer_x_h2o']
        self.ensure(p, t, x)
        columns = {'H2O': species_column(col, 'H2O'), 'CO2': species_column(col, 'CO2')}
        if extra:
            columns.update(extra)
        tau = gas_k(self.tables['gas']['H2O'], p, t, x) * columns['H2O'][:, None, None]
        for m in ('CO2', 'O3'):
            n = columns.get(m)
            if n is not None and np.any(n > 0):
                tau = overlap(tau, gas_k(self.tables['gas'][m], p, t) * n[:, None, None])
        n_tot = layer_number_density(col)
        dry = col['dry_fractions']
        cia = 0.0
        for a, b, _ in sp.CIA_PAIRS:
            fa = dry.get(a, 0.0) * (1 - x) if a != 'H2O' else x
            fb = dry.get(b, 0.0) * (1 - x) if b != 'H2O' else x
            cia = cia + cia_k(self.tables['cia'][f'{a}-{b}'], t) * (fa * fb * n_tot**2 * col['layer_dz_cm'])[:, None, None]
        return overlap(tau, cia)

    def fluxes(self, col, t_surface=None, extra=None, n_angles=4):
        """Upward and downward fluxes (W/m^2) at levels, surface first, and per band at the top."""
        tau = self.optical_depth(col, extra)
        nl = tau.shape[0]
        ts = col['t_k'][0] if t_surface is None else t_surface
        b = band_planck(col['t_k'])                  # (levels, bands)
        bs = band_planck(np.array([ts]))[0]          # (bands,)
        mu, wmu = lw.angles(n_angles)
        up = np.zeros((nl + 1, NB, NG)); down = np.zeros((nl + 1, NB, NG))
        for m, wt in zip(mu, wmu):
            weight = 2 * math.pi * wt * m
            i_down = np.zeros((NB, NG)); dl = [i_down]; terms = []
            for k in range(nl - 1, -1, -1):
                e, a, c = lw._linear_source_terms(tau[k] / m)
                terms.append((e, a, c))
                i_down = i_down * e + b[k][:, None] * a - (b[k] - b[k + 1])[:, None] * c
                dl.append(i_down)
            terms.reverse(); dl.reverse()
            i_up = np.broadcast_to(bs[:, None], (NB, NG)).copy()
            up[0] += weight * i_up
            for k in range(nl):
                e, a, c = terms[k]
                i_up = i_up * e + b[k + 1][:, None] * a - (b[k + 1] - b[k])[:, None] * c
                up[k + 1] += weight * i_up
            for k in range(nl + 1):
                down[k] += weight * dl[k]
        f_up = (up * W).sum(axis=(1, 2)); f_down = (down * W).sum(axis=(1, 2))
        return dict(up=f_up, down=f_down, olr=float(f_up[-1]), surface_down=float(f_down[0]),
                    olr_bands=(up[-1] * W).sum(axis=1))

    def heating_rate(self, col, cp, t_surface=None, extra=None):
        """Longwave heating rate of each layer (K/day)."""
        f = self.fluxes(col, t_surface, extra)
        net = f['up'] - f['down']
        absorbed = net[:-1] - net[1:]
        return absorbed / (col['layer_mass_kg_m2'] * cp) * 86400.0, f

    def fluxes_nonlte(self, col, t_surface=None, extra=None, o_mixing=None, nonlte=None):
        """Fluxes as fluxes(), with CO2 15-um emission out of local thermodynamic equilibrium aloft.

        In the bands of `nonlte` (500-820 cm^-1) and the layers above its start pressure, the source
        function is the two-level one, S = (J + eps B)/(1 + eps): J the absorption-weighted mean
        intensity of the layer, B the Planck function, eps the ratio of collisional deactivation of
        CO2(01101) by N2, O2 and O to its spontaneous emission. S follows from one linear solve per
        band (the band cores stay optically thick up to the model top, where iterating on J would
        converge very slowly), using isothermal-layer transfer. The difference that non-LTE makes
        in that transfer is added to the linear-in-optical-depth LTE fluxes, so where eps is large
        the result is exactly fluxes().
        """
        nonlte = nonlte or NonLTE()
        f = self.fluxes(col, t_surface, extra)
        tau = self.optical_depth(col, extra)
        ts = col['t_k'][0] if t_surface is None else t_surface
        eps = nonlte.epsilon(col, o_mixing)
        free = np.where(col['layer_p_pa'] < nonlte.start_pa)[0]
        if free.size == 0:
            return dict(f, epsilon=eps)
        b_layer = band_planck(col['layer_t_k'])
        b_surface = band_planck(np.array([ts]))[0]
        d_up = np.zeros(col['p_pa'].size); d_down = np.zeros(col['p_pa'].size)
        source = {}
        for band in nonlte.bands:
            tb = tau[:, band, :]
            lte = b_layer[:, band].copy()
            known = lte.copy(); known[free] = 0.0
            j_known, _, _ = _isothermal_transfer(tb, known[None, :], np.array([b_surface[band]]))
            units = np.zeros((free.size, lte.size)); units[np.arange(free.size), free] = 1.0
            j_unit, _, _ = _isothermal_transfer(tb, units, np.zeros(free.size))
            lam = j_unit[:, free].T                              # response of layer i to unit source in j
            a = np.diag(1.0 + eps[free]) - lam
            rhs = eps[free] * lte[free] + j_known[0, free]
            s = lte.copy()
            s[free] = np.linalg.solve(a, rhs)
            source[band] = s
            both = np.vstack([s, lte])
            _, up, down = _isothermal_transfer(tb, both, np.full(2, b_surface[band]))
            d_up += up[0] - up[1]
            d_down += down[0] - down[1]
        up, down = f['up'] + d_up, f['down'] + d_down
        return dict(up=up, down=down, olr=float(up[-1]), surface_down=float(down[0]),
                    olr_bands=f['olr_bands'], epsilon=eps, source=source)


@dataclass(frozen=True)
class NonLTE:
    """Two-level non-LTE for the CO2 bending mode (15 um).

    Einstein A of the v2 fundamental 1.5 s^-1; V-T deactivation of CO2(01101) by N2 (and O2, taken
    equal) 7e-17 T^0.5 + 6.7e-10 exp(-83.8 T^-1/3) cm^3 s^-1, and by atomic oxygen 3e-12 (T/300)^0.5
    cm^3 s^-1, as compiled by Lopez-Puertas and Taylor (2001), Non-LTE Radiative Transfer in the
    Atmosphere; the oxygen rate is uncertain by about a factor of two, which `o_rate_factor` varies.
    """
    bands: tuple = (2, 3, 4)
    start_pa: float = 50.0
    einstein_a: float = 1.5
    o_rate_factor: float = 1.0
    n2_rate_factor: float = 1.0

    def epsilon(self, col, o_mixing=None):
        t = col['layer_t_k']
        n = layer_number_density(col)
        x = col['layer_x_h2o']
        dry = col['dry_fractions']
        k_m = self.n2_rate_factor * (7e-17 * np.sqrt(t) + 6.7e-10 * np.exp(-83.8 * t ** (-1.0 / 3.0)))
        k_o = self.o_rate_factor * 3e-12 * np.sqrt(t / 300.0)
        n_m = (dry.get('N2', 0.0) + dry.get('O2', 0.0)) * (1 - x) * n
        n_o = (np.zeros_like(n) if o_mixing is None else np.asarray(o_mixing) * n)
        return (k_m * n_m + k_o * n_o) / self.einstein_a


def _isothermal_transfer(tau, source, surface, n_angles=4):
    """Isothermal-layer transfer of one band for a batch of layer sources.

    tau: (layers, g); source: (batch, layers); surface: (batch,) upward radiance at the ground.
    Returns the absorption-weighted mean intensity of each layer (batch, layers) and the upward and
    downward fluxes at levels (batch, levels), all band-integrated (W m^-2 sr^-1 and W m^-2).
    """
    nb, nl = source.shape
    mu, wmu = lw.angles(n_angles)
    jbar = np.zeros((nb, nl, NG))
    up = np.zeros((nb, nl + 1)); down = np.zeros((nb, nl + 1))
    for m, wt in zip(mu, wmu):
        x = tau / m
        e = np.exp(-x)
        mean = np.where(x > 1e-6, -np.expm1(-x) / np.maximum(x, 1e-300), 1 - x / 2)
        i = np.zeros((nb, NG))
        levels_down = [i]
        for k in range(nl - 1, -1, -1):
            s = source[:, k][:, None]
            jbar[:, k] += 0.5 * wt * (s + (i - s) * mean[k])
            i = i * e[k] + s * (1 - e[k])
            levels_down.append(i)
        levels_down.reverse()
        i = np.broadcast_to(surface[:, None], (nb, NG)).copy()
        levels_up = [i]
        for k in range(nl):
            s = source[:, k][:, None]
            jbar[:, k] += 0.5 * wt * (s + (i - s) * mean[k])
            i = i * e[k] + s * (1 - e[k])
            levels_up.append(i)
        weight = 2 * math.pi * wt * m
        up += weight * (np.stack(levels_up, axis=1) * W).sum(-1)
        down += weight * (np.stack(levels_down, axis=1) * W).sum(-1)
    absorb = tau * W
    j = (jbar * absorb[None]).sum(-1) / np.maximum(absorb.sum(-1), 1e-300)[None]
    return j, up, down
