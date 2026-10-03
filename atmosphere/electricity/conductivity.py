"""The air's small ions and conductivity: ion pairs made at rate q are lost by recombination with each other and by
attachment to aerosol particles and cloud droplets, so that in steady state

    q = alpha n^2 + S n,   S = sum_j beta_j Z_j,

with n the density of either sign, alpha the ion-ion recombination coefficient and S the attachment rate to particles
of concentration Z_j and attachment coefficient beta_j. The conductivity is sigma = e n (mu+ + mu-), each mobility
scaling inversely with the air's number density, and charge in the air relaxes with time constant eps0 / sigma.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ELEMENTARY_CHARGE = 1.602176634e-19          # C
EPSILON_0 = 8.8541878128e-12                 # F/m


def ion_density(q, alpha, sink=0.0):
    """Steady small-ion density (m-3) of each sign for production q (m-3 s-1), recombination alpha (m3 s-1) and
    attachment rate sink (s-1)."""
    q, alpha, sink = (np.asarray(a, float) for a in (q, alpha, sink))
    return 2.0 * q / (sink + np.sqrt(sink ** 2 + 4.0 * alpha * q))      # the positive root, stable when sink >> alpha n


def mobility(mu_standard, number_density, standard_number_density):
    """A small ion's mobility (m2 V-1 s-1) at the air's number density, from its reduced mobility at standard
    conditions."""
    return mu_standard * standard_number_density / np.asarray(number_density, float)


def conductivity(n, mu_positive, mu_negative):
    """Conductivity (S/m) from equal densities n (m-3) of positive and negative small ions."""
    return ELEMENTARY_CHARGE * np.asarray(n, float) * (np.asarray(mu_positive, float) + np.asarray(mu_negative, float))


def relaxation_time(sigma):
    """The time (s) over which charge in air of conductivity sigma (S/m) leaks away: eps0 / sigma."""
    return EPSILON_0 / np.asarray(sigma, float)


# --- The ion physics of air (literature note on the data drive, ionization_conductivity.md, sections 4.1-4.6) ------

BOLTZMANN = 1.380649e-23                     # J/K
W_EV = 35.0                                  # eV per ion pair, as CRAC:CRII
MOBILITY_STANDARD = (1.36e-4, 1.53e-4)       # m2/(V s), positive and negative small ions (Hirsikko et al. 2011, Horrak)
P0, T0 = 101325.0, 273.15
DETACHMENT_FREE = 0.0
DUST_RADIUS_UM = 0.05                        # aerosol radius for the attachment bracket
DROPLETS_M3 = 1.0e8                          # cloud droplets, as the CM1 runs set them (100 per cm3)
NONMUON_ATTENUATION_G_CM2 = 140.0            # the nucleonic and electromagnetic cascade below Earth's sea-level depth


def recombination(t_k, number_density_m3):
    """Ion-ion recombination coefficient (m3/s) of Brasseur and Chatel (1983): 6e-8 (300/T)^0.5 + 6e-26 [M] (300/T)^4
    cm3/s with [M] in cm-3, as recommended from the ground to 20 km (Zauner-Wieczorek et al. 2022)."""
    t = np.asarray(t_k, float)
    m_cm3 = np.asarray(number_density_m3, float) * 1e-6
    return (6.0e-8 * (300.0 / t) ** 0.5 + 6.0e-26 * m_cm3 * (300.0 / t) ** 4) * 1e-6


def mobility_air(mu_standard, p_pa, t_k):
    """A small ion's mobility (m2/(V s)) at pressure p and temperature T: inversely with pressure and as T^0.6 at
    fixed pressure (Tammet 1997, 1998)."""
    return mu_standard * (P0 / np.asarray(p_pa, float)) * (np.asarray(t_k, float) / T0) ** 0.6


def aerosol_attachment(radius_um, number_m3):
    """Attachment rate (1/s) of small ions to aerosol of one radius (um), with the coefficient of Tinsley and Zhou
    (2006) as printed by Baumgaertner et al. (2013): 4.36e-5 r - 9.2e-8 cm3/s above 0.01 um."""
    beta_cm3 = 4.36e-5 * radius_um - 9.2e-8
    return beta_cm3 * 1e-6 * np.asarray(number_m3, float)


def droplet_attachment(t_k, mobility, number_m3, radius_m):
    """Attachment rate (1/s) of small ions to cloud droplets by diffusion, 4 pi D N r with D = mu k T / e
    (Baumgaertner et al. 2013 Eqs. 7-8)."""
    diffusivity = mobility * BOLTZMANN * np.asarray(t_k, float) / ELEMENTARY_CHARGE
    return 4.0 * np.pi * diffusivity * np.asarray(number_m3, float) * np.asarray(radius_m, float)


def conductivity_of(q_m3_s, p_pa, t_k, sink_s=0.0):
    """Total conductivity (S/m) of air making q ion pairs per m3 per s at p and T, with an attachment sink (1/s)."""
    n_air = np.asarray(p_pa, float) / (BOLTZMANN * np.asarray(t_k, float))
    n = ion_density(q_m3_s, recombination(t_k, n_air), sink_s)
    mu_pos, mu_neg = (mobility_air(m, p_pa, t_k) for m in MOBILITY_STANDARD)
    return conductivity(n, mu_pos, mu_neg), n


def crac_crii(zip_path, pc_gv: float, phi_mv: float):
    """Depths (g/cm2) and cosmic-ray ion production (ion pairs per g per s) from the CRAC:CRII v2 tables at cutoff
    rigidity pc_gv and modulation potential phi_mv, interpolated linearly in both."""
    import zipfile
    depths, values = [], []
    with zipfile.ZipFile(zip_path) as z:
        for name in sorted(n for n in z.namelist() if n.startswith('S_') and n.endswith('.RES')):
            lines = z.read(name).decode('latin-1').split('\n')
            phi = np.array(lines[0].split()[1:], float)
            rows = np.array([l.split() for l in lines[1:] if l.strip()], float)
            pc, table = rows[:, 0], rows[:, 1:]
            by_phi = np.array([np.interp(phi_mv, phi, r) for r in table])
            depths.append(int(name[2:8]) / 100.0)
            values.append(float(np.interp(pc_gv, pc, by_phi)))
    order = np.argsort(depths)
    return np.array(depths)[order], np.array(values)[order]


def us_standard_atmosphere(z_km):
    """Pressure (Pa) and temperature (K) of the 1976 US standard atmosphere up to 32 km."""
    z = np.asarray(z_km, float)
    t = np.where(z <= 11.0, 288.15 - 6.5 * z, np.where(z <= 20.0, 216.65, 216.65 + (z - 20.0)))
    p = np.where(z <= 11.0, 101325.0 * (t / 288.15) ** 5.2559,
                 np.where(z <= 20.0, 22632.1 * np.exp(-(z - 11.0) * 1000.0 * 9.80665 / (287.053 * 216.65)),
                          5474.89 * (216.65 / t) ** 34.163))
    return p, t


def earth_column(z_km, pc_gv: float = 3.0, phi_mv: float = 400.0, sink_s=0.0):
    """Earth's cosmic-ray ion production (per m3 per s) and total and positive conductivity (S/m) at heights z_km of
    the US standard atmosphere, at cutoff pc_gv and modulation phi_mv, with an attachment sink (1/s)."""
    from atmosphere.electricity import fetch_inputs
    p, t = us_standard_atmosphere(z_km)
    depth = p / 9.80665 / 10.0
    dep, y = crac_crii(fetch_inputs.path('CRII_tables.zip'), pc_gv, phi_mv)
    q = np.exp(np.interp(depth, dep[dep > 0], np.log(y[dep > 0]))) * p / (287.053 * t) * 1e-3 * 1e6
    sigma, n = conductivity_of(q, p, t, sink_s)
    return q, sigma, ELEMENTARY_CHARGE * n * mobility_air(MOBILITY_STANDARD[0], p, t)


def gringel_positive(z_km, solar='minimum'):
    """Gringel's (1978) fit of the measured fair-weather positive conductivity over Germany (S/m), as given by
    Nicoll (2012) Eq. 9."""
    a, b, c, d = ((0.636, 0.36008, -0.008605, 0.00010331) if solar == 'minimum'
                  else (0.66837, 0.35653, -0.0095435, 0.00012313))
    z = np.asarray(z_km, float)
    return 1e-14 * np.exp(a + b * z + c * z ** 2 + d * z ** 3)


# --- The Open Moon's column -----------------------------------------------------------------------------------------

HERE = Path(__file__).resolve().parent
RESULTS = HERE / 'results'
SCHEMA = 'terluna.atmosphere.conductivity-column/1'
CM1_PROFILE = HERE.parents[1] / 'climate' / 'results' / 'crm' / 'mixed_phase_box_0e.json'
AEROSOL_M3 = (0.0, 1.0e8, 1.0e9)             # clear air: none, 100 and 1000 per cm3 of 0.05-um particles
CLOUD_WATER_G_M3 = (0.05, 0.1, 0.2)          # the charging zone's cloud water spans these (mixed_phase results)
EVIDENCE = ('Ion production from the CRAC:CRII v2 tables at zero cutoff above Earth\'s sea-level depth, with their '
            'muon part replaced by MCEq\'s muons for the Open Moon\'s column, and MCEq\'s muons and their decay electrons '
            'alone below it; recombination (Brasseur and Chatel 1983), small-ion mobility scaled with pressure and '
            'temperature (Tammet), attachment to aerosol (Tinsley and Zhou 2006) and to cloud droplets by diffusion. The '
            'same chain reproduces Earth\'s measured fair-weather conductivity (Gringel 1978) within 2-8 % from 5 to 30 km. '
            'Ground radioactivity, absent over the sea and confined to the lowest kilometres over land, is left out.')


def muon_yield(path: Path):
    """Depths (g/cm2) and the ionization per gram (ion pairs g-1 s-1) of muons and their decay electrons, from an
    MCEq product of muons.py."""
    rows = json.loads(Path(path).read_text())['rows']
    return np.array([r['depth_g_cm2'] for r in rows]), np.array([r['ion_pairs_per_g_s']['total'] for r in rows])


def moon_ion_production(depth_g_cm2, phi_mv: float, muons_moon: Path, muons_earth: Path):
    """Ion pairs per g per s at each depth of the Open Moon's column: Earth's polar CRAC:CRII production less Earth's
    MCEq muons (the nucleonic and electromagnetic cascade, attenuated beyond 1025 g/cm2) plus the Moon's MCEq muons."""
    from atmosphere.electricity import fetch_inputs
    x = np.asarray(depth_g_cm2, float)
    crac_x, crac_y = crac_crii(fetch_inputs.path('CRII_tables.zip'), 0.1, phi_mv)
    ex, ey = muon_yield(muons_earth)
    mx, my = muon_yield(muons_moon)
    log_interp = lambda xx, xs, ys: np.exp(np.interp(xx, xs, np.log(np.maximum(ys, 1e-30))))
    deepest = crac_x[-1]
    cascade_at = lambda xx: np.maximum(log_interp(np.minimum(xx, deepest), crac_x[crac_x > 0], crac_y[crac_x > 0])
                                       - log_interp(np.minimum(xx, deepest), ex, ey), 0.0)
    cascade = cascade_at(x) * np.exp(-np.maximum(x - deepest, 0.0) / NONMUON_ATTENUATION_G_CM2)
    return cascade + log_interp(x, mx, my), cascade


def moon_column(phi_mv: float = 400.0) -> dict:
    """The conductivity of the Open Moon's air from the ground to the CM1 top, in clear air and in cloud."""
    from shared.constants import MOON_SURFACE_GRAVITY
    prof = json.loads(CM1_PROFILE.read_text())['mean_profile']
    z = np.array(prof['z_km'])
    p = np.array(prof['pressure_hpa']) * 100.0
    t = np.array(prof['temperature_c']) + 273.15
    rho = np.array(prof['air_density_kg_m3'])
    depth = p / MOON_SURFACE_GRAVITY / 10.0                                  # g/cm2
    y, cascade = moon_ion_production(depth, phi_mv, RESULTS / 'muon_ionization_moon.json',
                                     RESULTS / 'muon_ionization_earth.json')
    q = y * rho * 1e-3                                                       # ion pairs per cm3 per s
    q_m3 = q * 1e6
    mu_pos = mobility_air(MOBILITY_STANDARD[0], p, t)
    clear = {f'{z_cm3 / 1e6:g}_per_cm3': conductivity_of(q_m3, p, t, aerosol_attachment(DUST_RADIUS_UM, z_cm3))[0]
             for z_cm3 in AEROSOL_M3}
    cloud = {}
    for lwc in CLOUD_WATER_G_M3:
        radius = (3.0 * lwc * 1e-3 / (4.0 * np.pi * 1000.0 * DROPLETS_M3)) ** (1.0 / 3.0)
        sink = droplet_attachment(t, mu_pos, DROPLETS_M3, radius)
        cloud[f'{lwc:g}_g_m3'] = conductivity_of(q_m3, p, t, sink)[0]
    return dict(z_km=z.tolist(), pressure_hpa=(p / 100.0).tolist(), temperature_k=t.tolist(), depth_g_cm2=depth.tolist(),
                ion_pairs_per_g_s=y.tolist(), cascade_share=(cascade / np.maximum(y, 1e-30)).tolist(),
                ion_pairs_per_cm3_s=q.tolist(),
                clear_air_s_m={k: v.tolist() for k, v in clear.items()},
                cloud_s_m={k: v.tolist() for k, v in cloud.items()})


def main(argv=None) -> int:
    import argparse, hashlib
    parser = argparse.ArgumentParser(description='The Open Moon\'s ionization and conductivity column')
    parser.parse_args(argv)
    result = dict(schema=SCHEMA, evidence=EVIDENCE,
                  reading_rule=('Columns are the CM1 equatorial box\'s mean profile (second lunar day). Conductivity is '
                                'the total of both signs. Clear air is given without aerosol and with 100 and 1000 per cm3 '
                                'of 0.05-um particles; cloud with 100 droplets per cm3 holding the stated cloud water. '
                                'Solar minimum and maximum are modulation potentials of 400 and 1000 MV.'),
                  solar_minimum=moon_column(400.0), solar_maximum=moon_column(1000.0))
    result['producer'] = dict(domain='atmosphere', files={f'electricity/{Path(__file__).name}':
                                                          hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / 'conductivity_moon.json'
    path.write_text(json.dumps(result, indent=1) + '\n')
    col = result['solar_minimum']
    for k in range(0, len(col['z_km']), 10):
        print(f"{col['z_km'][k]:6.1f} km {col['depth_g_cm2'][k]:6.0f} g/cm2  q {col['ion_pairs_per_cm3_s'][k]:8.3g} /cm3/s  "
              f"clear {col['clear_air_s_m']['0_per_cm3'][k]:8.2e}  cloud(0.1) {col['cloud_s_m']['0.1_g_m3'][k]:8.2e} S/m")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
