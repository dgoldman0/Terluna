"""Layer optical depths of a column for the thermal and solar calculations.

Absorption from HITRAN lines (H2O, CO2 and, in the solar range, O2), the MT_CKD
water continuum and HITRAN collision-induced absorption, plus Rayleigh scattering
for the solar range. Layers are independent, so they are computed in parallel
worker processes that share the loaded spectroscopy.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import multiprocessing as mp
import os
import numpy as np

from atmosphere.radiative_convective import spectroscopy as sp
from atmosphere.radiative_convective.thermodynamics import layer_number_density, species_column


@dataclass(frozen=True)
class OpticsConfig:
    """Spectral grid and line selection. Intensity floors are HITRAN 296-K values."""
    nu_min: float = 1.0
    nu_max: float = 3000.0
    step: float = 0.01
    molecules: tuple = ('H2O', 'CO2')
    s_min: tuple = (('H2O', 1e-28), ('CO2', 1e-28), ('O2', 0.0))
    continuum: bool = True
    cia: bool = True
    coarse_factor: int = 10
    tabulated: tuple = ('CO2', 'O2')  # air-broadened molecules served from a (ln p, T) node table
    processes: int = 0              # 0: half the available cores

    @property
    def grid(self):
        return sp.Grid.span(self.nu_min, self.nu_max, self.step)

    def floor(self, molecule):
        return dict(self.s_min)[molecule]


LONGWAVE = OpticsConfig()
SHORTWAVE = OpticsConfig(nu_min=2000.0, nu_max=20000.0, step=0.05, molecules=('H2O', 'CO2', 'O2'), coarse_factor=4)

_STATE: dict = {}

# Node table for air-broadened molecules: cross-sections at p = 10^(i/10) Pa and
# T = 10 j K, computed on demand, cached on disk, interpolated log-linearly in p
# and T. Its error against direct calculation is bounded by the tests.
TABLE_P_PER_DECADE = 10
TABLE_T_STEP = 10.0
_FLOOR = 1e-40


def _load(cfg: OpticsConfig):
    key = (cfg.molecules, cfg.s_min)
    if _STATE.get('key') == key:
        return
    lines = {m: sp.load_lines(m, cfg.floor(m)) for m in cfg.molecules}
    gids = np.concatenate([l.gid for l in lines.values()])
    _STATE.clear()
    _STATE.update(key=key, lines=lines, sums=sp.PartitionSums(gids), continuum=sp.WaterContinuum(),
                  cia=[(a, b, sp.CollisionInduced(f)) for a, b, f in sp.CIA_PAIRS])


def _table_key(cfg, molecule):
    keys = _STATE.setdefault('table_keys', {})
    grid_key = (molecule, cfg.nu_min, cfg.nu_max, cfg.step, cfg.coarse_factor)
    if grid_key not in keys:
        lines = _STATE['lines'][molecule]
        digest = hashlib.sha256()
        for arr in (lines.nu, lines.strength, lines.gamma_air, lines.elower, lines.n_air, lines.delta_air):
            digest.update(np.ascontiguousarray(arr).tobytes())
        digest.update(repr((grid_key, sp.KERNEL_VERSION, sp.VOIGT_CORE)).encode())
        keys[grid_key] = f'{molecule}_{digest.hexdigest()[:16]}'
    return keys[grid_key]


def _node_path(cfg, molecule, i, j):
    return sp.CACHE / 'nodes' / _table_key(cfg, molecule) / f'p{i:+04d}_t{j:03d}.npy'


def _node(args):
    grid, cfg, molecule, i, j = args
    path = _node_path(cfg, molecule, i, j)
    if not path.is_file():
        sigma = sp.line_absorption(grid, _STATE['lines'][molecule], _STATE['sums'], TABLE_T_STEP * j,
                                   10.0 ** (i / TABLE_P_PER_DECADE), 0.0, coarse_factor=cfg.coarse_factor)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix('.tmp.npy')
        np.save(tmp, sigma.astype(np.float32))
        tmp.replace(path)
    return path


def _brackets(p, t):
    x = np.log10(p) * TABLE_P_PER_DECADE
    i0 = int(np.floor(x)); wi = x - i0
    y = t / TABLE_T_STEP
    j0 = int(np.floor(y)); wj = y - j0
    return i0, wi, j0, wj


def _tabulated(cfg, molecule, p, t):
    i0, wi, j0, wj = _brackets(p, t)
    out = 0.0
    for di, a in ((0, 1 - wi), (1, wi)):
        for dj, b in ((0, 1 - wj), (1, wj)):
            if a * b == 0:
                continue
            out = out + a * b * np.log(np.maximum(np.load(_node_path(cfg, molecule, i0 + di, j0 + dj)), _FLOOR))
    return np.exp(out)


def _layer(args):
    grid, cfg, layer = args
    p, t, x_h2o = layer['p'], layer['t'], layer['x_h2o']
    nu = grid.nu
    tau = np.zeros(grid.size)
    for m in cfg.molecules:
        n = layer['columns'][m]
        if n <= 0:
            continue
        if m in cfg.tabulated:
            sigma = _tabulated(cfg, m, p, t)
        else:
            self_p = x_h2o * p if m == 'H2O' else 0.0
            sigma = sp.line_absorption(grid, _STATE['lines'][m], _STATE['sums'], t, p, self_p,
                                       pedestal=(m == 'H2O'), coarse_factor=cfg.coarse_factor)
        tau += sigma * n
    if cfg.continuum and layer['columns']['H2O'] > 0:
        s_c, f_c = _STATE['continuum'].cross_sections(nu, t, p, x_h2o)
        tau += (s_c + f_c) * layer['columns']['H2O']
    if cfg.cia:
        n_tot = layer['n_total']
        for a, b, table in _STATE['cia']:
            xa, xb = layer['fractions'].get(a, 0.0), layer['fractions'].get(b, 0.0)
            if xa > 0 and xb > 0:
                tau += table.coefficient(nu, t) * xa * xb * n_tot**2 * layer['dz_cm']
    return tau


def layer_inputs(col):
    """Per-layer state needed by the optical-depth workers."""
    n_tot = layer_number_density(col)
    dry = col['dry_fractions']
    out = []
    for k in range(len(col['layer_p_pa'])):
        x = col['layer_x_h2o'][k]
        fractions = {s: v * (1 - x) for s, v in dry.items()}
        fractions['H2O'] = x
        out.append(dict(p=float(col['layer_p_pa'][k]), t=float(col['layer_t_k'][k]), x_h2o=float(x),
                        n_total=float(n_tot[k]), dz_cm=float(col['layer_dz_cm'][k]), fractions=fractions,
                        columns={m: float(species_column(col, m)[k]) for m in ('H2O', 'CO2', 'O2')}))
    return out


def _map(func, jobs, cfg):
    procs = cfg.processes or max(1, (os.cpu_count() or 2) // 2)
    if procs == 1 or len(jobs) <= 1:
        return [func(j) for j in jobs]
    with mp.get_context('fork').Pool(min(procs, len(jobs))) as pool:
        return pool.map(func, jobs, chunksize=1)


def optical_depth(col, cfg: OpticsConfig = LONGWAVE):
    """Absorption optical depth (layers x wavenumbers) of a column."""
    _load(cfg)
    grid = cfg.grid
    layers = layer_inputs(col)
    needed = set()
    for layer in layers:
        i0, wi, j0, wj = _brackets(layer['p'], layer['t'])
        for m in cfg.tabulated:
            if m in cfg.molecules and layer['columns'][m] > 0:
                needed.update((m, i0 + di, j0 + dj) for di in (0, 1) for dj in (0, 1))
    missing = [(grid, cfg, m, i, j) for m, i, j in sorted(needed) if not _node_path(cfg, m, i, j).is_file()]
    _map(_node, missing, cfg)
    rows = _map(_layer, [(grid, cfg, layer) for layer in layers], cfg)
    return grid.nu, np.vstack(rows)


def rayleigh_optical_depth(col, nu):
    """Rayleigh scattering optical depth (layers x wavenumbers)."""
    co2 = col['dry_fractions'].get('CO2', 0.0)
    sigma = sp.rayleigh_cross_section(nu, co2)
    # Water vapour scatters about 0.75 times as strongly per molecule as air.
    weight = col['layer_column_cm2'] * (1 - 0.25 * col['layer_x_h2o'])
    return weight[:, None] * sigma[None, :]
