"""Cosmic-ray muons through a deep air column, and the ionization they make, computed with MCEq (matrix cascade
equations; Fedynitch et al.) in its own environment:

    /media/projectspace/terluna-research/venvs/mceq/bin/python -m atmosphere.electricity.muons earth
    /media/projectspace/terluna-research/venvs/mceq/bin/python -m atmosphere.electricity.muons moon

writes results/muon_ionization_<case>.json.

MCEq follows the hadronic cascade from the primary cosmic rays (the Hillas-Gaisser H3a spectrum at the top of the
air, unmodulated, with no geomagnetic cutoff) through pion and kaon production and decay to muons, with the muons'
energy loss and decay, along each zenith angle of a curved column. 'earth' is the CORSIKA US standard atmosphere;
'moon' is the Open Moon's: the mean profile of the CM1 equatorial box below its top (150 km, 22 hPa) and above it
the middle-atmosphere model's design-shield column (titania stack, 1.2 atm), joined in pressure and integrated
hydrostatically with gravity falling off with height, on the Moon's radius.

At each vertical depth the muon flux over eight zenith angles (cos 0.15 to 1, Gauss-Legendre) gives the scalar
flux, and the ionization per gram of air is that flux times the muons' total stopping power in air (PDG; it includes
the knock-on electrons and radiative losses, deposited nearby) over 35 eV per ion pair, as in CRAC:CRII, plus the
electrons of muon decays: in flight, 7/20 of the muon's energy; at rest, 35 MeV for each muon that stops.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
SCHEMA = 'terluna.atmosphere.muon-ionization/1'
W_EV = 35.0                                    # eV per ion pair, as CRAC:CRII uses
MU_MASS_GEV = 0.1056584
MU_LIFETIME_S = 2.1969811e-6
C_CM_S = 2.99792458e10
DECAY_IN_FLIGHT_SHARE = 7.0 / 20.0             # mean electron energy over muon energy, relativistic decay
DECAY_AT_REST_GEV = 0.035
RD = 287.04                                    # J/(kg K), as CM1
COS_LOW = 0.15
N_ANGLES = 8
INTERACTION_MODEL = 'SIBYLL23E'
CM1_PROFILE = ROOT / 'climate' / 'results' / 'crm' / 'mixed_phase_box_0e.json'
UPPER_PROFILE = HERE.parent / 'middle_atmosphere' / 'results' / 'profiles' / 'moon_1.2atm_titania_stack.csv'
DEPTHS = dict(earth=(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1030),
              moon=(100, 200, 300, 500, 700, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500,
                    7000, 7400))


def stopping_power(path: Path):
    """The PDG muon table for dry air: kinetic energy (GeV) and total stopping power (GeV cm2/g), from 10 MeV up."""
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 9:
            try:
                t_mev, dedx = float(parts[0]), float(parts[7])
            except ValueError:
                continue
            if t_mev >= 10.0:
                rows.append((t_mev / 1000.0, dedx / 1000.0))
    t, s = np.array(rows).T
    return t, s


def moon_profile(gravity: float, radius_m: float):
    """Heights (m), pressures (Pa), temperatures (K) and densities (kg/m3) of the Open Moon's column: CM1's mean
    profile below its top, the middle-atmosphere model's temperatures above, by pressure, integrated hydrostatically."""
    cm1 = json.loads(CM1_PROFILE.read_text())['mean_profile']
    z = np.array(cm1['z_km']) * 1000.0
    p = np.array(cm1['pressure_hpa']) * 100.0
    t = np.array(cm1['temperature_c']) + 273.15
    upper = list(csv.DictReader(UPPER_PROFILE.open()))
    p_up = np.array([float(r['p_pa']) for r in upper])
    t_up = np.array([float(r['t_k']) for r in upper])
    order = np.argsort(p_up)
    p_up, t_up = p_up[order], t_up[order]
    zs, ps, ts = list(z), list(p), list(t)
    for pp in np.geomspace(p[-1], p_up[0], 400)[1:]:
        tt = float(np.interp(np.log(pp), np.log(p_up), t_up))
        g = gravity * (radius_m / (radius_m + zs[-1])) ** 2
        t_mean = 0.5 * (tt + ts[-1])
        zs.append(zs[-1] + RD * t_mean / g * np.log(ps[-1] / pp))
        ps.append(pp)
        ts.append(tt)
    zs, ps, ts = (np.array(a) for a in (zs, ps, ts))
    return zs, ps, ts, ps / (RD * ts)


def build(case: str):
    """An MCEq run for the case and its density model, with the target heights (cm) of each vertical depth."""
    import MCEq.config as config
    import crflux.models as pm
    from MCEq.core import MCEqRun
    from MCEq.geometry.geometry import EarthGeometry
    import MCEq.geometry.density_profiles as dp

    if case == 'moon':
        sys.path.insert(0, str(ROOT))
        from shared.constants import MOON_RADIUS, MOON_SURFACE_GRAVITY
        z, p, t, rho = moon_profile(MOON_SURFACE_GRAVITY, MOON_RADIUS)
        config.r_E = MOON_RADIUS
        config.h_obs = 0.0
        config.h_atm = float(z[-1])
        config.max_density = float(rho[0] / 1000.0)
        log_rho = np.log(rho / 1000.0)                                          # g/cm3
        h_cm = z * 100.0

        class OpenMoonAtmosphere(dp.EarthsAtmosphere):
            def get_density(self, h):
                h = np.asarray(h, float)
                return np.exp(np.interp(h, h_cm, log_rho, right=log_rho[-1] - (h - h_cm[-1]) / 3.0e6))

        density = OpenMoonAtmosphere(geometry=EarthGeometry())
        column = dict(z_m=z, p_pa=p, t_k=t, rho=rho, gravity=MOON_SURFACE_GRAVITY)
    else:
        density = ('CORSIKA', ('BK_USStd', None))
        column = None
    run = MCEqRun(interaction_model=INTERACTION_MODEL, primary_model=(pm.HillasGaisser2012, 'H3a'), theta_deg=0.0,
                  density_model=density)
    return run, column


def depth_heights(run, case: str, column) -> dict:
    """The height (cm) of each target vertical depth (g/cm2), from the vertical path."""
    run.set_theta_deg(0.0)
    dm = run.density_model
    heights = {}
    for x in DEPTHS[case]:
        if x >= dm.max_X:
            continue
        heights[x] = float(dm.X2h(x))
    return heights


def ionization(case: str) -> dict:
    run, column = build(case)
    t_tab, s_tab = stopping_power(HERE / 'inputs' / 'muE_air_dry_1_atm.txt')
    e = run.e_grid                                                              # kinetic energy (GeV)
    de = run.e_widths
    s = np.interp(e, t_tab, s_tab)                                              # GeV cm2/g
    gamma_ = 1.0 + e / MU_MASS_GEV
    beta = np.sqrt(1.0 - 1.0 / gamma_ ** 2)
    heights = depth_heights(run, case, column)
    nodes, weights = np.polynomial.legendre.leggauss(N_ANGLES)
    cos_t = COS_LOW + (1.0 - COS_LOW) * (nodes + 1.0) / 2.0
    w_cos = weights * (1.0 - COS_LOW) / 2.0
    scalar = {x: np.zeros_like(e) for x in heights}
    vertical = {}
    for c, w in zip(cos_t, w_cos):
        run.set_theta_deg(float(np.degrees(np.arccos(c))))
        dm = run.density_model
        slant = np.array([float(dm.h2X(heights[x])) for x in heights])
        order = np.argsort(slant)
        run.solve(int_grid=slant[order])
        for i, k in enumerate(order):
            x = list(heights)[k]
            flux = run.get_solution('mu+', grid_idx=i) + run.get_solution('mu-', grid_idx=i)  # cm-2 s-1 sr-1 GeV-1
            scalar[x] += 2.0 * np.pi * w * flux
    run.set_theta_deg(0.0)
    dm = run.density_model
    slant = np.array([float(dm.h2X(heights[x])) for x in heights])
    run.solve(int_grid=np.sort(slant))
    for i, x in enumerate(heights):
        flux = run.get_solution('mu+', grid_idx=i) + run.get_solution('mu-', grid_idx=i)
        p = np.sqrt(e ** 2 + 2.0 * e * MU_MASS_GEV)
        vertical[x] = dict(above_1_gev_c_m2_s_sr=float((flux * de)[p >= 1.0].sum() * 1e4),
                           total_m2_s_sr=float((flux * de).sum() * 1e4))
    rows = []
    for x, h in heights.items():
        rho = float(run.density_model.get_density(h))                         # g/cm3 at that height
        phi = scalar[x]
        direct = float((phi * s * de).sum()) * 1e9 / W_EV                      # ion pairs per g per s
        in_flight = float((phi / (beta * C_CM_S) / (gamma_ * MU_LIFETIME_S) / rho * de
                           * DECAY_IN_FLIGHT_SHARE * (e + MU_MASS_GEV)).sum()) * 1e9 / W_EV
        at_rest = float(phi[0] * s[0]) * DECAY_AT_REST_GEV * 1e9 / W_EV          # muons slowing past the grid's foot
        rows.append(dict(depth_g_cm2=float(x), height_km=h / 1e5, air_density_g_cm3=rho,
                         scalar_flux_cm2_s=float((phi * de).sum()), vertical=vertical[x],
                         ion_pairs_per_g_s=dict(muons=direct, decay_in_flight=in_flight, decay_at_rest=at_rest,
                                                total=direct + in_flight + at_rest)))
    out = dict(schema=SCHEMA, case=case, model=dict(mceq='MCEq 1.4.2', interaction=INTERACTION_MODEL,
                                                   primaries='Hillas-Gaisser 2012 H3a, no cutoff, unmodulated',
                                                   zenith='%d Gauss-Legendre angles, cos %.2f to 1' % (N_ANGLES, COS_LOW),
                                                   w_ev=W_EV), rows=rows)
    if column is not None:
        out['column'] = dict(surface_depth_g_cm2=float(run.density_model.max_X), top_km=float(column['z_m'][-1] / 1000.0),
                             joined_at_hpa=float(column['p_pa'][110] / 100.0),
                             upper_air='atmosphere/middle_atmosphere/results/profiles/moon_1.2atm_titania_stack.csv')
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('case', choices=('earth', 'moon'))
    args = parser.parse_args(argv)
    result = ionization(args.case)
    result['producer'] = dict(domain='atmosphere', files={f'electricity/{Path(__file__).name}':
                                                          hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f'muon_ionization_{args.case}.json'
    path.write_text(json.dumps(result, indent=1) + '\n')
    for r in result['rows']:
        print(f"{r['depth_g_cm2']:7.0f} g/cm2  {r['height_km']:7.1f} km  muons {r['ion_pairs_per_g_s']['muons']:9.3g}"
              f"  decays {r['ion_pairs_per_g_s']['decay_in_flight'] + r['ion_pairs_per_g_s']['decay_at_rest']:9.3g}"
              f"  I_v(>1 GeV/c) {r['vertical']['above_1_gev_c_m2_s_sr']:8.3g} m-2 s-1 sr-1")
    return 0


if __name__ == '__main__':
    sys.exit(main())
