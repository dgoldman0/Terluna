"""Line-by-line reference fluxes on the equilibrium profiles, for checking a 3-D model's radiation code.

    python -m atmosphere.middle_atmosphere.benchmarks     # writes results/radiation_benchmarks.json

A general circulation model's fast radiation scheme has to reproduce, on the
same column, the fluxes this repository computes line by line. For each listed
case this rebuilds the final profile (levels, layer temperature, water vapour,
ozone) and records:

- thermal fluxes up and down at every level from the 0.01 cm^-1 line-by-line
  model (H2O, CO2 and O3 lines, MT_CKD 4.3, collision-induced absorption),
  beside the correlated-k fluxes the equilibrium used;
- solar fluxes (direct, diffuse down, up) at every level for sunlight filtered
  by the case's shield, at solar zenith cosines 1.0 and 0.5 and as the global
  mean over the sunlit disk: 202-500 nm from photolysis.py (O2, O3, NO2 and
  other trace absorbers, collision-induced bands, Rayleigh scattering) and
  500 nm-5 um from the line-by-line solar model with ozone's Chappuis band.

The ozone lines use the 0.01 cm^-1 grid, which does not resolve Doppler cores at
low pressure; heating rates above ~10 Pa are therefore less exact than fluxes.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import multiprocessing as mp
import numpy as np

from shared.constants import SOLAR_CONSTANT
from atmosphere.radiative_convective import optics, longwave as lw, shortwave as sw, spectroscopy as sp, ck, climate as cl
from atmosphere.middle_atmosphere import equilibrium as eq, photolysis as ph, run as rn, lbl_check

HERE = Path(__file__).resolve().parent
CASES = ('earth_control', 'moon_1.2atm_titania_stack', 'moon_1.0atm_titania_stack', 'moon_1.2atm_edge_200nm')
MU0 = (1.0, 0.5)
_O3 = {}


def _o3_sigma(args):
    p, t = args
    if not _O3:
        _O3['lines'] = sp.load_lines('O3', ck.S_MIN['O3'])
        _O3['sums'] = sp.PartitionSums(_O3['lines'].gid)
    grid = optics.LONGWAVE.grid
    return sp.line_absorption(grid, _O3['lines'], _O3['sums'], t, p, 0.0, coarse_factor=optics.LONGWAVE.coarse_factor)


def thermal(col, o3_columns, processes):
    nu, tau = optics.optical_depth(col, optics.LONGWAVE)
    jobs = [(float(p), float(t)) for p, t in zip(col['layer_p_pa'], col['layer_t_k'])]
    with mp.get_context('fork').Pool(processes) as pool:
        for k, sigma in enumerate(pool.imap(_o3_sigma, jobs, chunksize=1)):
            tau[k] += sigma * o3_columns[k]
    up, down = lw.fluxes(nu, tau, col['t_k'], col['t_k'][0])
    return lw.integrate(nu, up), lw.integrate(nu, down)


def solar(col, air, o3_mixing, shield, albedo, state_densities):
    """Level fluxes (surface first) of direct, diffuse-down and up sunlight, W/m^2, per zenith case."""
    o3_columns = o3_mixing * col['layer_column_cm2']
    sun = ph.Sunlight(shield)
    tau_abs_uv, tau_sca_uv, _ = sun.optics(col, state_densities)
    below = ph.CENTRES_NM < ph.HEATING_LIMIT_NM
    nu, tau_abs = optics.optical_depth(col, eq.SW_CONFIG)
    lam = 1e7 / nu
    tau_abs = tau_abs + o3_columns[:, None] * np.interp(lam, ph.CENTRES_NM, ph.CrossSections().o3_room,
                                                        left=0.0, right=0.0)[None, :]
    tau_r = optics.rayleigh_optical_depth(col, nu)
    spectrum, longward = sw.solar_spectrum(nu, shield=shield)
    weights = np.full(nu.size, nu[1] - nu[0]); weights[0] *= 0.5; weights[-1] *= 0.5
    radius = (col['planet'].radius_m + col['z_m'])[::-1]
    nodes, node_w = sw.disk_nodes(6)

    def at(mu0):
        """Fluxes for sunlight at zenith cosine mu0, per unit area of ground (incident S * mu0)."""
        out = np.zeros((3, col['p_pa'].size))
        d, f, u = sw.column_fluxes(tau_abs_uv[::-1][:, below], tau_sca_uv[::-1][:, below],
                                   np.zeros_like(tau_sca_uv[::-1][:, below]), albedo, mu0, radius)
        e = sun.energy[below]
        out += np.array([d @ e, f @ e, u @ e])[:, ::-1]
        chunk = 40000
        for a in range(0, nu.size, chunk):
            sl = slice(a, a + chunk)
            d, f, u = sw.column_fluxes(tau_abs[::-1, sl], tau_r[::-1, sl], np.zeros_like(tau_r[::-1, sl]),
                                       albedo, mu0, radius)
            e = spectrum[sl] * weights[sl]
            out += np.array([d @ e, f @ e, u @ e])[:, ::-1]
        return out

    result = {f'mu0={m:g}': at(m) for m in MU0}
    mean = sum(0.5 * w * at(m) for m, w in zip(nodes, node_w))
    result['global_mean'] = mean
    return result, dict(beyond_5um_w_m2=0.25 * float(longward or 0.0),
                        below_202nm='not included (every shield considered blocks it)')


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--results', type=Path, default=HERE / 'results')
    parser.add_argument('--processes', type=int, default=6)
    parser.add_argument('--only', nargs='*')
    args = parser.parse_args(argv)
    plan = rn.cases()
    out_path = args.results / 'radiation_benchmarks.json'
    record = json.loads(out_path.read_text()) if out_path.is_file() else dict(
        schema='terluna.atmosphere.radiation-benchmark/1',
        evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
        reading_rule=('Arrays are surface first. Level fluxes in W/m^2: thermal up/down; solar direct (on a '
                      'horizontal surface), diffuse down and up, for sunlight at the stated zenith cosine '
                      '(incident S mu0, S the solar constant filtered by the shield) and for the global mean '
                      '(incident S/4). Mixing ratios are mole fractions of the layer. The profiles are the '
                      'equilibrium results in profiles/<case>.csv.'),
        solar_constant_w_m2=SOLAR_CONSTANT, cases={})
    solver_lw = ck.CKLongwave()
    for name in (args.only or CASES):
        case = plan[name]
        prof, t, trop, o3, col, _ = lbl_check.rebuild(case, args.results)
        rows = [dict(zip(r.keys(), r.values())) for r in __import__('csv').DictReader(
            open(args.results / 'profiles' / f'{name}.csv', newline=''))]
        dens = {s: np.array([float(r[s]) for r in rows]) * eq.ch.fixed_densities(col)['M']
                for s in ('O3', 'NO2', 'NO3', 'N2O5', 'N2O', 'HNO3', 'H2O2')}
        up_lbl, down_lbl = thermal(col, o3 * col['layer_column_cm2'], args.processes)
        f_ck = solver_lw.fluxes(col, extra={'O3': o3 * col['layer_column_cm2']})
        shield = cl.load_shield(case.shield)
        sol, notes = solar(col, prof.air, o3, shield, case.surface_albedo, dens)
        record['cases'][name] = dict(
            planet=case.planet.name, gravity_m_s2=case.planet.surface_gravity, radius_m=case.planet.radius_m,
            dry_air_fractions=col['dry_fractions'], shield=case.shield, surface_albedo=case.surface_albedo,
            surface_temperature_k=case.surface_temperature_k,
            level_p_pa=col['p_pa'].tolist(), level_t_k=col['t_k'].tolist(), level_z_m=col['z_m'].tolist(),
            layer_p_pa=col['layer_p_pa'].tolist(), layer_t_k=col['layer_t_k'].tolist(),
            layer_h2o=col['layer_x_h2o'].tolist(), layer_o3=o3.tolist(),
            thermal_lbl=dict(up=up_lbl.tolist(), down=down_lbl.tolist()),
            thermal_ck=dict(up=f_ck['up'].tolist(), down=f_ck['down'].tolist()),
            solar={k: dict(direct=v[0].tolist(), diffuse_down=v[1].tolist(), up=v[2].tolist()) for k, v in sol.items()},
            solar_notes=notes)
        g = sol['global_mean']
        print(f"{name:32s} OLR lbl {up_lbl[-1]:7.2f} ck {f_ck['olr']:7.2f}  DLR lbl {down_lbl[0]:7.2f} ck "
              f"{f_ck['surface_down']:7.2f}  global-mean ASR {g[0][-1] + g[1][-1] - g[2][-1]:7.2f} W/m2", flush=True)
        out_path.write_text(json.dumps(record) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
