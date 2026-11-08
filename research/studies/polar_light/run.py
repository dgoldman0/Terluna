"""Polar light on the turned Moon (8 November 2026).

What light the polar regions of the Open Moon receive after the Moon's turn, from 60 degrees of latitude to the poles in
both hemispheres, for plants and for people, and how much of it the spin axis's tilt and its precession set.

- The Sun's path. A mean-geometry ephemeris in the turned Moon's own frame. The Moon stays in Cassini state 2, so its
  spin axis keeps its tilt eps to the ecliptic (1.613 deg in the design state, up to 2.5 deg while the pores fill) and
  precesses with the orbit's node every 18.6 years. The Sun's declination is asin(sin eps sin u), u advancing through
  the eclipse year; the Sun's distance follows Earth's orbit with the run's year and eccentricity. The precession's phase
  psi is the value of u at perihelion; the run holds psi at 270 deg, the south's summer at perihelion every year. The
  Earth stands at sub-Earth latitude -asin(sin(i + eps) sin F), F = D + u, with D the synodic phase (0 at full Earth),
  and at sub-Earth longitude -L sin M, the optical libration in longitude over the anomalistic month (shared constants).
- Refraction. Bouguer's invariant n r cos(h) integrated through the design column (the middle atmosphere's 1.2-atm
  profiles at 288 and 298 K, blended in log density to the solved column's 294.9 K surface), with the refractivity of
  standard air at 550 nm (Edlen 1966) scaled by number density. The Earth control's column is the check against Earth's
  measured horizontal refraction of about 34'.
- The light at the ground. The light calendar's evaluator (illumination/calendar, Response) gives direct and diffuse
  sunlight and the Earth's light by elevation and phase on the solved 1.2-atm column. Photosynthetic photons per lux
  come from the solved standard solution's 114 channels: the diffuse ground flux of each channel at every solar
  elevation from -90 to 90 deg, and the direct beam by slant optical depth through its 100 shells. Clouds are the run's
  cover by latitude and hour angle (climatology y10-49) in a two-state sky, clear or overcast, overcast passing 0.25 of
  the clear sky's light (Kasten & Czeplak 1980 at eight oktas).
- Terrain. LOLA's heights above the GRAIL geoid at 16 px/deg, with the seas and the forty-year routing's rain-fed lakes
  as flat water, carried to the turned frame and resampled on polar stereographic grids of 1.9 km about each turned
  pole; horizons in 72 azimuths out to 265 km with the Moon's curvature exact; the Sun's centre against the horizon
  through a lunar day at each declination; and the clear sky's diffuse light and twilight that each horizon still
  admits, from the solved atlas's angular radiance.

    nice -n 10 ionice -c 3 python3 -m research.studies.polar_light.run [--figures] [--workers 4]
    climate/gcm/.venv/bin/python -m research.studies.polar_light.run --extract-run   # the run's sunlight (diagnostic)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import time
from pathlib import Path

import numpy as np

from shared.constants import (ANOMALISTIC_MONTH_DAYS, AU, AVOGADRO, BOLTZMANN, EARTH_MOON_DISTANCE, EARTH_RADIUS,
                              LIBRATION_LONGITUDE_DEG, MOON_ORBIT_INCLINATION_DEG, MOON_RADIUS, PLANCK, SPEED_OF_LIGHT,
                              SUN_RADIUS, SYNODIC_MONTH_DAYS)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.polar-light/1'
OUT = HERE / 'results' / 'polar_light.json'
DRIVE = Path('/media/projectspace/terluna-research/session-tools/agents/polar_light')
SPHERICAL = Path('/media/projectspace/terluna-research/optical-comfort-runs/spherical')
SOLUTION = SPHERICAL / 'moon_1.2atm_standard.npz'       # the solved sky's standard solution (its atlas's source)
ATLAS = SPHERICAL / 'moon_1.2atm_atlas.npz'

RUN = 'A28R_301k_design_l15_snow_s924_turned_refit'
RUN_YEARS = (10, 49)
RUN_DIR = ROOT / 'climate' / 'gcm' / 'runs' / RUN / 'model'
RUN_EXTRACT = DRIVE / 'run_sunlight_y10-49.npz'
HEIGHTS_CACHE = DRIVE / 'heights_16ppd.npy'
HORIZON_CACHE = DRIVE / 'horizons_{cap}.npz'
TENTHS = np.array([-180.0 + 36.0 * (k + 0.5) for k in range(10)])     # hour-angle bins of the run's 3-day means

INPUTS = {
    'transfer': ('illumination/calendar/results/transfer.json', 'terluna.illumination.calendar-transfer/1'),
    'solved_sky': ('illumination/sky/results/solved_sky.json', 'terluna.illumination.solved-spherical-sky/1'),
    'orientation': ('geology/results/moon_orientation.json', 'terluna.geology.moon-orientation/1'),
    'nodal_period': ('geology/interior_sources.json', None),
    'climatology': (f'climate/gcm/products/climatology_{RUN}_y10-49.npz', 'terluna.climate.gcm-climatology/1'),
    'plant': ('biosphere/canopy/results/plant.json', None),
    'living_light': ('research/studies/living_light/results/living_light.json', 'terluna.research.living-light/1'),
    'light_limitation': ('research/studies/light_limitation/results/light_limitation.json',
                         'terluna.research.light-limitation/1'),
    'atlas': ('geography/results/atlas.json', 'terluna.geography.atlas/1'),
    'drainage': ('geography/turned_y10-49/results/drainage.json', 'terluna.geography.drainage/1'),
    'lakes': ('geography/turned_y10-49/products/drainage_28pct_16ppd.npz', None),
    'profile_288': ('atmosphere/middle_atmosphere/results/profiles/moon_1.2atm_titania_stack.csv', None),
    'profile_298': ('atmosphere/middle_atmosphere/results/profiles/moon_1.2atm_titania_stack_ts298.csv', None),
    'profile_earth': ('atmosphere/middle_atmosphere/results/profiles/earth_control.csv', None),
    'case_288': ('atmosphere/middle_atmosphere/results/cases/moon_1.2atm_titania_stack.json', None),
    'case_298': ('atmosphere/middle_atmosphere/results/cases/moon_1.2atm_titania_stack_ts298.json', None),
    'case_earth': ('atmosphere/middle_atmosphere/results/cases/earth_control.json', None),
    'ldem_16': ('geography/inputs/ldem_16.img', None),
    'ldem_4': ('geography/inputs/ldem_4.img', None),
    'gravity': ('geography/inputs/jggrx_0420a_sha.tab', None),
}
IMPORTED = ('illumination/calendar/validate.py', 'illumination/calendar/disk.py', 'illumination/calendar/geometry.py',
            'climate/gcm/turned.py', 'geography/topography.py')

# ---- study settings -------------------------------------------------------------------------------------------------
TILTS = (1.613, 2.5)                    # the design state with seas and lakes; the elastic crust's top while pores fill
PSI_RUN, PSI_NORTH = 270.0, 90.0        # precession phase: the south's summer at perihelion (the run), the north's
LATS = np.concatenate([np.arange(60.0, 85.0, 1.0), np.arange(85.0, 90.01, 0.5)])
TABLE_LATS = (60.0, 62.0, 64.0, 65.0, 66.0, 68.0, 70.0, 72.0, 75.0, 78.0, 80.0, 82.0, 84.0, 86.0, 88.0, 89.0, 90.0)
SEASONS = np.arange(0.0, 360.0, 5.0)    # eclipse-year phase u at the lunar day's local noon
SEASON_TABLE = np.arange(0.0, 360.0, 30.0)
STEP_H = 0.25                           # time step through a lunar day, hours
SECTORS = {'near side': 0.0, 'limb': 90.0, 'far side': 180.0}   # turned longitude; the sub-Earth point at 0 N 0 E
DARK_UMOL = 1.0                         # the plant model's dark: under 1 umol/m2/s of photosynthetic light
OVERCAST = 0.25                         # overcast global light / clear (Kasten & Czeplak 1980, eight oktas)
OVERCAST_RANGE = (0.1, 0.5)             # sensitivity, thick low overcast to thin high cloud
EDLEN_550 = 2.7784e-4                   # (n - 1) of standard dry air at 550 nm, 15 C and 101,325 Pa (Edlen 1966)
STANDARD_AIR = (101325.0, 288.15)
PEOPLE_LUX = {'walking': 1.0, 'faces at 4 m': 1.6, 'civil twilight on Earth': 3.0, 'reading': 8.0,
              'close work outdoors': 100.0, 'indoor reading and writing': 500.0}
EK_UMOL = (100.0, 300.0)                # the light-limitation study's saturation onset, Ek (umol/m2/s)
VGPM_HALF_MOL = 4.1                     # the light-limitation study's VGPM half-saturation daily light (mol/m2/d)
DLI_NEEDS = {'seedlings and cuttings': (5.0, 10.0), 'impatiens, a shade crop': (8.0, 12.0), 'lettuce': (12.0, 17.0),
             'spinach': (14.0, 20.0), 'tomato': (20.0, 30.0)}   # mol/m2/d recommended (Stallknecht 2025, table 3)
FAUST_RANGE = (5.0, 19.0, 43.0)         # bedding plants' dry mass rose with DLI from 5 to 43, begonia's and impatiens' to 19
# ---- terrain --------------------------------------------------------------------------------------------------------
DX_KM = MOON_RADIUS / 1000.0 * math.radians(1.0 / 16.0)          # 1.895 km, 16 px/deg at the pole
CAP_DEG = 30.0                                                    # analysis: 60 deg to the pole
REACH_PX = np.unique(np.round(np.geomspace(1.0, 140.0, 48))).astype(int)   # ray samples, grid pixels, to 265 km
N_AZ = 72
TERRAIN_DECS = 360                      # Sun positions per lunar day for direct sunlight
TERRAIN_SKY_STEPS = 48                  # Sun positions per lunar day for the diffuse light at terrain pixels
TERRAIN_SUBSAMPLE = 2                   # every second pixel in each direction for the diffuse light

SOURCES = [
    dict(key='edlen1966', citation='Edlen B 1966, The refractive index of air, Metrologia 2, 71-80',
         doi_or_url='https://doi.org/10.1088/0026-1394/2/2/002',
         use='EDLEN_550: (n - 1) of standard dry air at 550 nm from the dispersion formula, 2.7784e-4 at 15 C and '
             '101,325 Pa, scaled by number density through the column',
         read='the dispersion formula for standard air as widely reproduced; the paper itself was not opened'),
    dict(key='kasten_czeplak1980', citation='Kasten F, Czeplak G 1980, Solar and terrestrial radiation dependent on '
                                            'the amount and type of cloud, Solar Energy 24, 177-189',
         doi_or_url='https://doi.org/10.1016/0038-092X(80)90391-6',
         use='OVERCAST: the clear-sky index under N oktas of cloud, 1 - 0.75 (N/8)^3.4, 0.25 at eight oktas; fitted to '
             'ten years of hourly data at Hamburg with the Sun above 5 deg, and applied here beyond that range, at low '
             'Sun and in twilight',
         read='the relation and its fitting range as review literature reproduces them (search of 8 November 2026); '
              'the paper itself was not opened'),
    dict(key='faust2005', citation='Faust JE, Holcombe V, Rajapakse NC, Layne DR 2005, The effect of daily light '
                                   'integral on bedding plant growth and flowering, HortScience 40, 645-649',
         doi_or_url='https://doi.org/10.21273/HORTSCI.40.3.645',
         use='FAUST_RANGE: bedding plants grown at 5-43 mol/m2/d; total dry mass rose with DLI over the whole range '
             'except in begonia and impatiens, which rose to 19 mol/m2/d',
         read='the abstract'),
    dict(key='stallknecht2025', citation='Stallknecht E 2025, Calculating and using daily light integral (DLI): an '
                                         'introductory guide, Virginia Cooperative Extension SPES-720',
         doi_or_url='https://www.pubs.ext.vt.edu/SPES/spes-720/spes-720.html',
         use='DLI_NEEDS: recommended DLI for seedlings and cuttings 5-10, impatiens 8-12, lettuce 12-17, spinach '
             '14-20 and tomato 20-30 mol/m2/d (table 3, citing Faust et al. 2005, Pramuk & Runkle 2005, Dou et al. 2018)',
         read='the guide\'s table 3 and its text on the categories'),
    dict(key='usno_twilight', citation='U.S. Naval Observatory, Astronomical Applications Department, Rise, set and '
                                       'twilight definitions',
         doi_or_url='https://aa.usno.navy.mil/faq/RST_defs',
         use='civil twilight on Earth ends with the Sun 6 deg down, the limit for ordinary outdoor activity without '
             'lamps; its light on the ground is the solved sky\'s Earth control at -6 deg, 3.0 lux',
         read='the definitions page'),
]


# ==== identity ======================================================================================================
def sha(path, n=64):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()[:n]


def read_json(key):
    path, schema = INPUTS[key]
    data = json.loads((ROOT / path).read_text())
    if schema and data.get('schema') != schema:
        raise ValueError(f'{path}: expected {schema}')
    return data


def plain(x, digits=4):
    """JSON-ready numbers rounded to significant digits."""
    if isinstance(x, dict):
        return {str(k): plain(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v, digits) for v in x]
    if isinstance(x, np.ndarray):
        return plain(x.tolist(), digits)
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        x = float(x)
        if not math.isfinite(x):
            return None
        return float(f'{x:.{digits}g}')
    return x


# ==== astronomy =====================================================================================================
def periods():
    """The cycles the polar light follows: the synodic month (shared constants), the run's year and Earth's orbital
    eccentricity (climatology metadata), the node's regression (geology's source) and the eclipse year they give."""
    meta = json.loads(str(np.load(ROOT / INPUTS['climatology'][0])['metadata']))
    planet = meta['configuration']['planet']
    node = next(s['values_used']['node_days'] for s in read_json('nodal_period')['sources'] if s['key'] == 'nodal_period')
    year = float(planet['year'])
    return dict(synodic_days=SYNODIC_MONTH_DAYS, synodic_hours=SYNODIC_MONTH_DAYS * 24.0, year_days=year,
                eccentricity=float(planet['eccentricity']), node_days=float(node),
                eclipse_year_days=1.0 / (1.0 / year + 1.0 / node), run_tilt_deg=float(planet['obliquity']))


def declination(u_deg, eps_deg):
    return np.degrees(np.arcsin(np.sin(np.radians(eps_deg)) * np.sin(np.radians(u_deg))))


def distance_factor(u_deg, psi_deg, ecc):
    """(1 AU / r)^2 with the true anomaly nu = u - psi."""
    nu = np.radians(np.asarray(u_deg) - psi_deg)
    return ((1.0 + ecc * np.cos(nu)) / (1.0 - ecc ** 2)) ** 2


def season_weights(u_deg, psi_deg, per):
    """The time spent per unit of u: 1 / (d nu/dt + d psi/dt), the Sun faster near perihelion, normalised."""
    ecc = per['eccentricity']
    nu = np.radians(np.asarray(u_deg) - psi_deg)
    rate = (360.0 / per['year_days']) * (1 + ecc * np.cos(nu)) ** 2 / (1 - ecc ** 2) ** 1.5 + 360.0 / per['node_days']
    w = 1.0 / rate
    return w / w.sum()


def lunar_days(u0, eps_deg, per, lon=0.0, m0=None, step_h=STEP_H):
    """One lunar day per eclipse-year phase u0 (deg, at local noon), local midnight to midnight: hour angle H, the
    Sun's declination and subsolar longitude, the sub-Earth latitude and longitude (deg), time from noon (h)."""
    u0 = np.atleast_1d(np.asarray(u0, float))
    n = int(round(per['synodic_hours'] / step_h))
    t = (np.arange(n) + 0.5) / n * per['synodic_hours'] - per['synodic_hours'] / 2
    H = t / per['synodic_hours'] * 360.0
    u = u0[:, None] + 360.0 * (t / 24.0) / per['eclipse_year_days']
    lam_ss = lon - H[None, :] * np.ones_like(u)
    D = 180.0 - lam_ss                                   # 0 at full Earth, when the Sun stands over the far side's centre
    F = D + u                                            # the Moon's argument of latitude
    inc = math.radians(MOON_ORBIT_INCLINATION_DEG + eps_deg)
    beta_e = -np.degrees(np.arcsin(math.sin(inc) * np.sin(np.radians(F))))
    if m0 is None:
        m0 = (np.arange(u0.size) * 137.50776) % 360.0    # spread the anomalistic phase over the days (golden angle)
    M = np.asarray(m0)[:, None] + 360.0 * (t / 24.0) / ANOMALISTIC_MONTH_DAYS
    lam_e = -LIBRATION_LONGITUDE_DEG * np.sin(np.radians(M))
    return dict(t_h=t, H=H, u=u, dec=declination(u, eps_deg), lam_ss=lam_ss, beta_e=beta_e, lam_e=lam_e)


def sun_elevation(lat, lon, dec, lam_ss):
    phi, d = math.radians(lat), np.radians(dec)
    return np.degrees(np.arcsin(np.sin(phi) * np.sin(d) + np.cos(phi) * np.cos(d) * np.cos(np.radians(lon - lam_ss))))


def unit(lat_deg, lon_deg):
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la) * np.ones_like(lo)], axis=-1)


# ==== refraction ====================================================================================================
def read_profile(path):
    rows = np.genfromtxt(ROOT / path, delimiter=',', names=True)
    return rows['z_km'] * 1e3, rows['p_pa'], rows['t_k']


def column_density(case='moon'):
    """Number density (m^-3) against height (m): the 1.2-atm profiles blended in log density to the solved column's
    surface temperature, with the solved surface added; or the Earth control."""
    if case == 'earth':
        z, p, t = read_profile(INPUTS['profile_earth'][0])
        c = read_json('case_earth')
        z0, p0, t0 = 0.0, c['surface_pressure_pa'], c['surface_temperature_k']
        return np.concatenate([[z0], z]), np.concatenate([[p0 / (BOLTZMANN * t0)], p / (BOLTZMANN * t)])
    solved = read_json('solved_sky')['worlds']['moon_1.2atm']['optical_column']['column']
    ts, ps = solved['surface_temperature_k'], solved['surface_pressure_pa']
    c1, c2 = read_json('case_288'), read_json('case_298')
    w = (ts - c1['surface_temperature_k']) / (c2['surface_temperature_k'] - c1['surface_temperature_k'])
    z1, p1, t1 = read_profile(INPUTS['profile_288'][0])
    z2, p2, t2 = read_profile(INPUTS['profile_298'][0])
    z = np.unique(np.concatenate([z1, z2]))
    ln = (1 - w) * np.interp(z, z1, np.log(p1 / (BOLTZMANN * t1))) + w * np.interp(z, z2, np.log(p2 / (BOLTZMANN * t2)))
    return np.concatenate([[0.0], z]), np.concatenate([[ps / (BOLTZMANN * ts)], np.exp(ln)])


def refraction(z_m, number_density, radius_m, apparent_deg, n_s=40000):
    """Astronomical refraction (deg) for an observer on the ground at the given apparent elevations: the bending
    integral -int k n'/(n sqrt(n^2 r^2 - k^2)) dr with k = n0 r0 cos(h), in s = sqrt(r - r0), which removes the
    horizon's singularity; log density linear in height between the profile's levels."""
    ns = STANDARD_AIR[0] / (BOLTZMANN * STANDARD_AIR[1])
    lnN = np.log(number_density)
    top = z_m[-1]
    s = np.linspace(0.0, math.sqrt(top), n_s)
    z = s ** 2
    nn = np.exp(np.interp(z, z_m, lnN))
    dlnN = np.gradient(lnN, z_m)
    slope = np.interp(z, z_m, dlnN)                       # d ln N / dz
    refr = EDLEN_550 * nn / ns                            # n - 1
    n = 1.0 + refr
    dndr = refr * slope
    r = radius_m + z
    out = []
    for h in np.atleast_1d(apparent_deg):
        k = n[0] * radius_m * math.cos(math.radians(h))
        root = np.sqrt(np.maximum((n * r) ** 2 - k ** 2, 0.0))
        f = np.zeros_like(s)
        ok = s > 0
        f[ok] = -k * dndr[ok] / (n[ok] * root[ok]) * 2 * s[ok]
        if h == 0.0:                                      # the limit at s -> 0
            d = 2 * n[0] * radius_m * (n[0] + radius_m * dndr[0])
            f[0] = -k * dndr[0] / (n[0] * math.sqrt(d)) * 2
        out.append(math.degrees(np.trapezoid(f, s)))
    return np.array(out), EDLEN_550 * nn[0] / ns


def refraction_summary():
    out = {}
    elev = np.array([0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0])
    for case, radius in (('moon', MOON_RADIUS), ('earth', EARTH_RADIUS)):
        z, N = column_density(case)
        r, n0 = refraction(z, N, radius, elev)
        # the density scale height over the lowest 20 km, for the exponential-atmosphere check
        sel = z <= 20e3
        H = -1.0 / np.polyfit(z[sel], np.log(N[sel]), 1)[0]
        out[case] = dict(apparent_elevation_deg=elev, refraction_deg=r, surface_refractivity=n0,
                         density_scale_height_km=H / 1e3,
                         exponential_estimate_deg=math.degrees(n0 * math.sqrt(math.pi * radius / (2 * H))))
    return out


# ==== light at the ground ===========================================================================================
class Light:
    """Clear-sky light on level ground: the calendar evaluator's sunlight and Earthlight (lux), and photons per lux."""

    def __init__(self, transfer, solution_path=SOLUTION):
        from illumination.calendar.validate import Response
        self.response = Response(transfer)
        self.sun_radius = math.degrees(math.asin(SUN_RADIUS / AU))
        self.elev = np.round(np.arange(-90.0, 90.0001, 0.02), 4)
        base = dict(sun_radius_deg=self.sun_radius, sun_distance_factor=1.0, earth_phase_deg=180.0,
                    source_separation_deg=180.0, earth_radius_deg=0.0)
        parts = [self.response.light({**base, 'sun_elevation_deg': float(e)}, include_earth=False) for e in self.elev]
        self.direct = np.array([p['solar_direct'] for p in parts])
        self.diffuse = np.array([p['solar_diffuse'] for p in parts])
        self.photons = photon_ratios(solution_path)
        self.q_dir = np.interp(self.elev, self.photons['direct_elevation_deg'], self.photons['direct_umol_per_lux'])
        self.q_dif = np.interp(self.elev, self.photons['diffuse_elevation_deg'], self.photons['diffuse_umol_per_lux'])
        self.earth_elev = np.arange(-3.0, 50.001, 0.25)
        self.earth_phase = self.response.phase_angles
        er = math.degrees(math.asin(EARTH_RADIUS / EARTH_MOON_DISTANCE))
        table = np.zeros((self.earth_phase.size, self.earth_elev.size))
        for i, ph in enumerate(self.earth_phase):
            for j, e in enumerate(self.earth_elev):
                g = dict(sun_elevation_deg=-90.0, sun_radius_deg=self.sun_radius, sun_distance_factor=1.0,
                         earth_phase_deg=float(ph), earth_elevation_deg=float(e), earth_radius_deg=er,
                         earth_bright_limb_rad=0.0, earth_distance_m=EARTH_MOON_DISTANCE, earth_sunlight_factor=1.0,
                         source_separation_deg=180.0)
                table[i, j] = self.response.light(g)['earth']
        self.earth_table = table

    def sun(self, h):
        """Direct and diffuse sunlight (lux) and photons (umol/m2/s) at the mean distance."""
        d, f = np.interp(h, self.elev, self.direct), np.interp(h, self.elev, self.diffuse)
        q = d * np.interp(h, self.elev, self.q_dir) + f * np.interp(h, self.elev, self.q_dif)
        return d, f, q

    def earth(self, h, phase):
        """The Earth's light (lux) from the (phase, elevation) table, bilinear; zero below the table."""
        i = np.clip((phase - self.earth_phase[0]) / (self.earth_phase[1] - self.earth_phase[0]), 0,
                    self.earth_phase.size - 1.000001)
        j = np.clip((h - self.earth_elev[0]) / (self.earth_elev[1] - self.earth_elev[0]), 0,
                    self.earth_elev.size - 1.000001)
        i0, j0 = i.astype(int), j.astype(int)
        fi, fj = i - i0, j - j0
        t = self.earth_table
        v = ((1 - fi) * (1 - fj) * t[i0, j0] + fi * (1 - fj) * t[i0 + 1, j0] + (1 - fi) * fj * t[i0, j0 + 1] +
             fi * fj * t[i0 + 1, j0 + 1])
        return np.where(h < self.earth_elev[0], 0.0, v)


def photon_ratios(path=SOLUTION):
    """Photosynthetic photons (400-700 nm, umol/m2/s) per lux for the diffuse sky light at each solar elevation of the
    solved standard solution (-90 to 90 deg) and for the direct beam (slant optical depth through the 100 shells)."""
    z = np.load(path)
    energy, xyz, wl, band = z['energy'], z['xyz'], z['wavelength'], z['band']
    par = (band[:, 0] >= 400.0) & (band[:, 1] <= 700.0)
    per_flux = energy / (PLANCK * SPEED_OF_LIGHT / (wl * 1e-9)) / AVOGADRO * 1e6
    moments = z['moments'][0, :, :, 5]                                   # diffuse ground flux per unit incident
    a = np.degrees(z['a'])
    dif_lux, dif_par = moments @ xyz[:, 1], moments[:, par] @ per_flux[par]
    edges, ext = z['edges'], z['scattering'] + z['absorption']
    elev = np.concatenate([np.arange(0.0, 3.0, 0.05), np.arange(3.0, 90.001, 0.5)])
    b = edges[0] * np.cos(np.radians(elev))[:, None]
    path_m = np.sqrt(np.maximum(edges[None, 1:] ** 2 - b ** 2, 0)) - np.sqrt(np.maximum(edges[None, :-1] ** 2 - b ** 2, 0))
    T = np.exp(-(path_m @ ext))
    dir_lux, dir_par = T @ xyz[:, 1], T[:, par] @ per_flux[par]
    ok = dif_lux > 0
    return dict(diffuse_elevation_deg=a[ok], diffuse_umol_per_lux=dif_par[ok] / dif_lux[ok],
                direct_elevation_deg=elev, direct_umol_per_lux=dir_par / dir_lux,
                direct_normal_lux=dir_lux, diffuse_lux=dif_lux[ok], channels=int(wl.size), par_channels=int(par.sum()))


# ==== clouds ========================================================================================================
def cloud_cover():
    """The run's cloud cover by latitude row and hour angle (climatology y10-49, clt_sun), and its run diagnostics."""
    z = np.load(ROOT / INPUTS['climatology'][0])
    return dict(lat=z['lat'].astype(float), hour=z['hour_angle_deg'].astype(float), clt=z['clt_sun'].astype(float),
                clt_mean=z['clt_mean'].astype(float), lsm=z['lsm'].astype(float))


def cover_at(cc, lat, H):
    """Cloud cover at a latitude (same hemisphere's rows, nearest pole row held beyond it) and hour angles H (deg)."""
    rows = np.flatnonzero(np.sign(cc['lat']) == np.sign(lat))
    la = cc['lat'][rows]
    order = np.argsort(la)
    by_row = cc['clt'][rows][order]
    col = np.array([np.interp(lat, la[order], by_row[:, k]) for k in range(cc['hour'].size)])
    hh = np.concatenate([cc['hour'] - 360.0, cc['hour'], cc['hour'] + 360.0])
    return np.interp((np.asarray(H) + 180.0) % 360.0 - 180.0, hh, np.tile(col, 3))


def cover_pixels(cc, lat_pix, H):
    """Cloud cover at many points: bilinear in latitude (rows of the points' hemisphere) and hour angle."""
    out = np.zeros(lat_pix.shape)
    hh = np.concatenate([cc['hour'] - 360.0, cc['hour'], cc['hour'] + 360.0])
    for hemi in (1.0, -1.0):
        sel = np.sign(lat_pix) == hemi
        if not sel.any():
            continue
        rows = np.flatnonzero(np.sign(cc['lat']) == hemi)
        order = np.argsort(cc['lat'][rows])
        la, tab = cc['lat'][rows][order], np.tile(cc['clt'][rows][order], (1, 3))
        x = np.interp(lat_pix[sel], la, np.arange(la.size))
        r0 = np.clip(np.floor(x).astype(int), 0, la.size - 2)
        fr = x - r0
        hcol = np.interp((H[sel] + 180.0) % 360.0 - 180.0, hh, np.arange(hh.size))
        c0 = np.clip(np.floor(hcol).astype(int), 0, hh.size - 2)
        fc = hcol - c0
        out[sel] = ((1 - fr) * ((1 - fc) * tab[r0, c0] + fc * tab[r0, c0 + 1]) +
                    fr * ((1 - fc) * tab[r0 + 1, c0] + fc * tab[r0 + 1, c0 + 1]))
    return out


def run_sunlight_diagnostic(path=RUN_EXTRACT):
    """The run's all-sky over clear-sky surface sunlight by polar row and tenth of the lunar day (from the extract)."""
    if not path.exists():
        return None
    z = np.load(path)
    lat = z['lat']
    rows = np.flatnonzero(np.abs(lat) > 55)
    out = dict(path=str(path), sha256=sha(path), run=str(z['run']), years=z['years'].tolist(), tenths_deg=TENTHS,
               rows={})
    for r in rows:
        n = np.maximum(z['n_land'][:, r], 1)
        out['rows'][f'{lat[r]:.2f}'] = dict(
            land_all_over_clear=z['rss_land'][:, r] / np.maximum(z['rsscs_land'][:, r], 1e-12),
            land_all_over_clear_day=z['rss_land'][:, r].sum() / max(z['rsscs_land'][:, r].sum(), 1e-12),
            clear_sky_w_m2=z['rsscs_land'][:, r] / n, mean_cos_zenith=z['czen_land'][:, r] / n,
            cloud_cover=z['clt_land'][:, r] / n)
    return out


# ==== flat ground through the lunar day =============================================================================
def site_days(light, per, lat, eps, sector_lon=None, psi=None, cc=None, overcast=OVERCAST):
    """One site's 72 seasonal lunar days: the Sun's path, and the light at the ground clear and under the run's clouds
    (two-state); sector_lon adds the Earth's light at that turned longitude."""
    from illumination.calendar.geometry import geometry_from_vectors
    lon = 0.0 if sector_lon is None else sector_lon
    s = lunar_days(SEASONS, eps, per, lon=lon)
    h = sun_elevation(lat, lon, s['dec'], s['lam_ss'])
    d, f, q = light.sun(h)
    if psi is not None:
        fd = distance_factor(SEASONS, psi, per['eccentricity'])[:, None]
        d, f, q = d * fd, f * fd, q * fd
    earth = np.zeros_like(h)
    q_earth = np.zeros_like(h)
    e_elev = np.full_like(h, np.nan)
    if sector_lon is not None:
        n = h.size
        g = dict(earth=unit(s['beta_e'].ravel(), s['lam_e'].ravel()), earth_distance_m=np.full(n, EARTH_MOON_DISTANCE),
                 sun=unit(s['dec'].ravel(), s['lam_ss'].ravel()), sun_distance_m=np.full(n, AU))
        geo = geometry_from_vectors(g, lon, lat)
        e_elev = geo['earth_elevation_deg'].reshape(h.shape)
        earth = light.earth(e_elev, geo['earth_phase_deg'].reshape(h.shape))
        q_earth = earth * np.interp(e_elev, light.elev, light.q_dif)
    lux = d + f + earth
    ppfd = q + q_earth
    cover = cover_at(cc, lat, s['H'])[None, :] * np.ones_like(h) if cc is not None else np.zeros_like(h)
    return dict(s=s, h=h, direct=d, diffuse=f, earth=earth, earth_elev=e_elev, lux=lux, ppfd=ppfd, ppfd_earth=q_earth,
                cover=cover, overcast=overcast)


def summarise_days(x, dt_h, sun_radius, refr0, levels=PEOPLE_LUX, dark=DARK_UMOL):
    """Per lunar day: hours of the Sun above the horizon, extreme elevations, photons and the hours under the plant
    model's dark and under (and over) people's light levels; clear, under the run's clouds and overcast."""
    h, lux, ppfd, c, ov = x['h'], x['lux'], x['ppfd'], x['cover'], x['overcast']
    mean_f = 1.0 - (1.0 - ov) * c
    out = dict(
        hours_centre=(h > 0).sum(1) * dt_h, hours_limb=(h > -sun_radius).sum(1) * dt_h,
        hours_limb_refracted=(h > -(sun_radius + refr0)).sum(1) * dt_h,
        max_elevation=h.max(1), min_elevation=h.min(1),
        mol_clear=ppfd.sum(1) * dt_h * 3600e-6, mol_cloud=(ppfd * mean_f).sum(1) * dt_h * 3600e-6,
        mol_overcast=ppfd.sum(1) * ov * dt_h * 3600e-6,
        mol_earth=x['ppfd_earth'].sum(1) * dt_h * 3600e-6,
        dark_clear=(ppfd < dark).sum(1) * dt_h,
        dark_cloud=(((1 - c) * (ppfd < dark)) + c * (ppfd * ov < dark)).sum(1) * dt_h,
        dark_overcast=(ppfd * ov < dark).sum(1) * dt_h,
        darkest_lux=lux.min(1), darkest_ppfd=ppfd.min(1), brightest_lux=lux.max(1),
        mean_lux_clear=lux.mean(1), mean_lux_cloud=(lux * mean_f).mean(1),
        earth_max_lux=x['earth'].max(1),
        earth_share_at_darkest=np.take_along_axis(x['earth'] / np.maximum(lux, 1e-30), lux.argmin(1)[:, None], 1)[:, 0],
        mean_cover=c.mean(1),
        above_ek_low=(ppfd > EK_UMOL[0]).sum(1) * dt_h, above_ek_high=(ppfd > EK_UMOL[1]).sum(1) * dt_h)
    n24 = int(round(24.0 / dt_h))
    ext = np.concatenate([ppfd, ppfd[:, :n24]], axis=1)
    run24 = np.cumsum(ext, axis=1)
    win = (run24[:, n24:] - run24[:, :-n24]) * dt_h * 3600e-6
    out['dli24_min'], out['dli24_max'] = win.min(1), win.max(1)
    for name, level in levels.items():
        out[f'below {name}'] = (lux < level).sum(1) * dt_h
        out[f'below {name} cloud'] = (((1 - c) * (lux < level)) + c * (lux * ov < level)).sum(1) * dt_h
    return out


# ==== terrain =======================================================================================================
def heights():
    """LOLA heights above the GRAIL geoid at 16 px/deg (m), cached on the drive with the inputs' identity."""
    atlas = read_json('atlas')
    meta = DRIVE / 'heights_16ppd.json'
    ident = {k: sha(ROOT / INPUTS[k][0]) for k in ('ldem_16', 'ldem_4', 'gravity')}
    ident['geoid_degree'] = atlas['grid']['geoid_degree']
    if HEIGHTS_CACHE.exists() and meta.exists() and json.loads(meta.read_text()) == ident:
        return np.load(HEIGHTS_CACHE), ident
    from geography import topography as tp
    h, _, _ = tp.height_above_geoid(16, atlas['grid']['geoid_degree'])
    h = h.astype(np.float32)
    np.save(HEIGHTS_CACHE, h)
    meta.write_text(json.dumps(ident))
    return h, ident


def water_surface(h):
    """The ground and water surface (m): seas at the atlas's level, the forty-year routing's lakes at their levels."""
    drainage = read_json('drainage')
    level = float(drainage['sea_level_m'])
    lakes = np.load(ROOT / INPUTS['lakes'][0])
    surf = h.copy()
    sea = h < level
    surf[sea] = level
    lake = np.zeros(h.shape, bool)
    lake.flat[lakes['lake_cell']] = True
    surf.flat[lakes['lake_cell']] = lakes['lake_level_m']
    return surf, sea, lake, level


def stereo_grid(cap, reach_deg):
    """A polar stereographic grid about the turned pole (km): unit vectors, turned latitude, true cell area."""
    R = MOON_RADIUS / 1000.0
    c_max = math.radians(CAP_DEG + reach_deg)
    rho_max = 2 * R * math.tan(c_max / 2)
    n = int(math.ceil(rho_max / DX_KM))
    x = (np.arange(2 * n + 1) - n) * DX_KM
    X, Y = np.meshgrid(x, x)
    rho2 = X ** 2 + Y ** 2
    sgn = 1.0 if cap == 'north' else -1.0
    nz = sgn * (4 * R * R - rho2) / (4 * R * R + rho2)
    nx, ny = 4 * R * X / (4 * R * R + rho2), 4 * R * Y / (4 * R * R + rho2)
    lat = np.degrees(np.arcsin(np.clip(nz, -1, 1)))
    lon = np.degrees(np.arctan2(ny, nx))
    colat = np.radians(90.0 - sgn * lat)
    area = DX_KM ** 2 * ((1 + np.cos(colat)) / 2) ** 2                   # grid area / k^2, k = 2 / (1 + cos c)
    return dict(n=n, x=x, nvec=np.stack([nx, ny, nz], -1), lat=lat, lon=lon, area=area, sgn=sgn)


def sample_today(grid, field, rotation, order=1):
    """A 16 px/deg field of today's ground read at the grid's turned points."""
    from scipy.ndimage import map_coordinates
    from climate.gcm.turned import to_today
    lat0, lon0 = to_today(grid['lat'], grid['lon'], rotation)
    rows = (90.0 - lat0) * 16.0 - 0.5
    cols = np.mod(lon0, 360.0) * 16.0 - 0.5
    pad = np.concatenate([field[:, -1:], field, field[:, :1]], axis=1)
    return map_coordinates(pad, [rows, cols + 1.0], order=order, mode='nearest')


_H = {}                                 # the horizon workers' arrays, inherited through fork


def _horizon_chunk(az):
    """Horizon tangents for one block of azimuths (worker; arrays in _H)."""
    n0, nvec_n, nvec_e, h0, surf, n, sgn = (_H[k] for k in ('n0', 'north', 'east', 'h0', 'surf', 'n', 'sgn'))
    R = MOON_RADIUS / 1000.0
    out = np.full((n0.shape[0], len(az)), -np.inf)
    for k, A in enumerate(az):
        a = math.radians(A)
        dvec = math.cos(a) * nvec_n + math.sin(a) * nvec_e
        best = out[:, k]
        for px in REACH_PX:
            d = px * DX_KM
            ang = d / R
            ca, sa = math.cos(ang), math.sin(ang)
            p = ca * n0 + sa * dvec
            den = 1.0 + sgn * p[:, 2]
            ix = np.rint(2 * R * p[:, 0] / den / DX_KM).astype(np.int64) + n
            iy = np.rint(2 * R * p[:, 1] / den / DX_KM).astype(np.int64) + n
            hj = surf[iy, ix]
            tanv = ((R + hj) * ca - (R + h0)) / ((R + hj) * sa)
            np.maximum(best, tanv, out=best)
    return out


def local_axes(n0):
    """Local north and east unit vectors at each point of the turned frame."""
    east = np.cross(np.array([0.0, 0.0, 1.0]), n0)
    norm = np.linalg.norm(east, axis=1)
    pole = norm < 1e-12                                   # at the pole itself azimuths take +y as east
    east[pole] = np.array([0.0, 1.0, 0.0])
    norm[pole] = 1.0
    east /= norm[:, None]
    return np.cross(n0, east), east


def horizons(surf_km, grid, mask, workers=4):
    """Horizon elevation (deg) in N_AZ azimuths (clockwise from north) for each analysis pixel of a cap: the highest
    elevation of the ground or water surface along each great circle out to REACH_PX, curvature exact."""
    idx = np.flatnonzero(mask.ravel())
    n0 = grid['nvec'].reshape(-1, 3)[idx]
    north, east = local_axes(n0)
    _H.update(n0=n0, north=north, east=east, h0=surf_km.ravel()[idx], surf=surf_km, n=grid['n'], sgn=grid['sgn'])
    az = np.arange(N_AZ) * 360.0 / N_AZ
    blocks = np.array_split(az, max(1, workers) * 3)
    if workers > 1:
        import multiprocessing as mp
        with mp.get_context('fork').Pool(workers) as pool:
            parts = pool.map(_horizon_chunk, blocks)
    else:
        parts = [_horizon_chunk(b) for b in blocks]
    _H.clear()
    return np.degrees(np.arctan(np.concatenate(parts, axis=1))).astype(np.float32), n0, north, east


def interp_horizon(hz, az_deg):
    """Horizon elevation at each pixel's azimuth (deg), linear between the N_AZ directions."""
    x = (az_deg % 360.0) / (360.0 / N_AZ)
    k0 = np.floor(x).astype(np.int64) % N_AZ
    f = x - np.floor(x)
    rows = np.arange(hz.shape[0])
    return hz[rows, k0] * (1 - f) + hz[rows, (k0 + 1) % N_AZ] * f


def sunlit_hours(hz, n0, north, east, decs, refr=None, steps=TERRAIN_DECS):
    """Hours a lunar day with the Sun's centre above each pixel's horizon (and above flat ground), per declination."""
    lit = np.zeros((len(decs), n0.shape[0]))
    flat = np.zeros((len(decs), n0.shape[0]))
    for i, dec in enumerate(decs):
        for lam in np.arange(steps) * 360.0 / steps:
            s = unit(dec, lam)
            elev = np.degrees(np.arcsin(np.clip(n0 @ s, -1, 1)))
            if refr is not None:
                elev = elev + np.interp(elev, refr[0], refr[1], left=refr[1][0], right=0.0)
            az = np.degrees(np.arctan2(east @ s, north @ s))
            lit[i] += elev > interp_horizon(hz, az)
            flat[i] += elev > 0
    scale = SYNODIC_MONTH_DAYS * 24.0 / steps
    return lit * scale, flat * scale


def sky_table(atlas_path=ATLAS):
    """The solved atlas's radiance integrated from the horizon up, G[sun, azimuth, mu] = int_0^mu L mu' dmu' (lux per
    radian of azimuth, mu = sin elevation, the atlas's interpolation rule), and its open-sky total by Sun height."""
    z = np.load(atlas_path)
    L = z['xyz'][..., 1]                                   # (suns, elevations, azimuths 0..180), cd/m2
    mu_nodes = np.sin(np.radians(z['elevations']))
    mu = np.linspace(0.0, 0.5, 501)                       # up to 30 deg
    Li = np.stack([np.stack([np.interp(mu, mu_nodes, L[i, :, j]) for j in range(L.shape[2])], 0)
                   for i in range(L.shape[0])], 0)        # (suns, az, mu)
    integrand = Li * mu
    G = np.concatenate([np.zeros(Li.shape[:2] + (1,)),
                        np.cumsum(0.5 * (integrand[..., 1:] + integrand[..., :-1]) * np.diff(mu), axis=-1)], axis=-1)
    total = z['diffuse_horizontal_xyz'][:, 1]
    return dict(suns=z['suns'], az=z['azimuths'], mu=mu, G=G, open=total)


def blocked_fraction(tab, hz, sun_elev, sun_az):
    """The share of the clear sky's diffuse light on horizontal ground that each pixel's horizon blocks, for the Sun
    at the given elevation (scalar or per pixel) and azimuth (per pixel)."""
    suns = tab['suns']
    se = np.clip(np.asarray(sun_elev, float), suns[0], suns[-1]) * np.ones(hz.shape[0])
    i = np.clip(np.searchsorted(suns, se) - 1, 0, suns.size - 2)
    w = (se - suns[i]) / (suns[i + 1] - suns[i])
    dmu = tab['mu'][1] - tab['mu'][0]
    daz = 2 * math.pi / N_AZ
    mu_k = np.clip(np.sin(np.radians(np.maximum(hz, 0.0))), 0, tab['mu'][-1] - 1e-9) / dmu
    blocked = np.zeros(hz.shape[0])
    open_ = (1 - w) * tab['open'][i] + w * tab['open'][i + 1]
    step = tab['az'][1] - tab['az'][0]
    for k in range(N_AZ):
        rel = np.abs(((k * 360.0 / N_AZ) - sun_az + 180.0) % 360.0 - 180.0)
        j = np.rint(rel / step).astype(np.int64)
        m0 = mu_k[:, k].astype(np.int64)
        fm = mu_k[:, k] - m0
        g0 = tab['G'][i, j, m0] * (1 - fm) + tab['G'][i, j, m0 + 1] * fm
        g1 = tab['G'][i + 1, j, m0] * (1 - fm) + tab['G'][i + 1, j, m0 + 1] * fm
        blocked += ((1 - w) * g0 + w * g1) * daz
    return np.clip(blocked / np.maximum(open_, 1e-30), 0.0, 1.0)


def terrain_light(light, tab, hz, n0, north, east, decs, cc, lat_pix, steps=TERRAIN_SKY_STEPS, overcast=OVERCAST):
    """Light at terrain pixels through a lunar day per declination: direct sunlight where the Sun clears the horizon,
    the diffuse light the horizon admits; mean photons, the plant model's dark hours clear and under the run's clouds,
    and hours under people's levels."""
    out = []
    dt = SYNODIC_MONTH_DAYS * 24.0 / steps
    for dec in decs:
        ppfd = np.zeros((steps, n0.shape[0]))
        lux = np.zeros((steps, n0.shape[0]))
        cover = np.zeros((steps, n0.shape[0]))
        blk = np.zeros((steps, n0.shape[0]))
        for k, lam in enumerate(np.arange(steps) * 360.0 / steps):
            s = unit(dec, lam)
            elev = np.degrees(np.arcsin(np.clip(n0 @ s, -1, 1)))
            az = np.degrees(np.arctan2(east @ s, north @ s))
            d, f, q = light.sun(elev)
            vis = elev > interp_horizon(hz, az)
            b = blocked_fraction(tab, hz, elev, az)
            qd = d * np.interp(elev, light.elev, light.q_dir)
            qf = f * np.interp(elev, light.elev, light.q_dif)
            ppfd[k] = qd * vis + qf * (1 - b)
            lux[k] = d * vis + f * (1 - b)
            blk[k] = b
            lon = np.degrees(np.arctan2(n0[:, 1], n0[:, 0]))
            cover[k] = cover_pixels(cc, lat_pix, (lon - lam + 180.0) % 360.0 - 180.0)
        c = cover
        res = dict(dec=dec,
                   mol_day_clear=ppfd.sum(0) * dt * 3600e-6 / SYNODIC_MONTH_DAYS,
                   mol_day_cloud=(ppfd * (1 - (1 - overcast) * c)).sum(0) * dt * 3600e-6 / SYNODIC_MONTH_DAYS,
                   dark_clear=hours_below(ppfd, DARK_UMOL, dt),
                   dark_cloud=hours_below_cloud(ppfd, DARK_UMOL, dt, c, overcast),
                   below_reading=hours_below(lux, PEOPLE_LUX['reading'], dt),
                   below_work=hours_below(lux, PEOPLE_LUX['close work outdoors'], dt),
                   darkest_lux=lux.min(0), blocked_at_darkest=np.take_along_axis(blk, lux.argmin(0)[None, :], 0)[0])
        out.append(res)
    return out


def below_fraction(series, level):
    """The share of each interval between neighbouring samples (steps, pixels; periodic) spent below a level, by
    log-linear crossing (daily_life's rule)."""
    y = np.log(np.maximum(series, 1e-30))
    line = math.log(level)
    y1 = np.roll(y, -1, axis=0)
    both = (y < line) & (y1 < line)
    fall = (y >= line) & (y1 < line)
    rise = (y < line) & (y1 >= line)
    frac = np.where(fall | rise, (y - line) / np.where(y != y1, y - y1, 1.0), 0.0)
    return both + fall * (1 - frac) + rise * frac


def hours_below(series, level, dt):
    return (below_fraction(series, level) * dt).sum(0)


def hours_below_cloud(series, level, dt, cover, overcast):
    """Expected hours below a level in the two-state sky: clear with probability 1 - cover, overcast with cover."""
    c = 0.5 * (cover + np.roll(cover, -1, axis=0))
    return (((1 - c) * below_fraction(series, level) + c * below_fraction(series * overcast, level)) * dt).sum(0)


BANDS = ((60.0, 70.0), (70.0, 80.0), (80.0, 85.0), (85.0, 90.01))


def cap_surface(cap, workers=4, log=print):
    """One cap's grid, surface, water and horizons; horizons are cached on the drive under the surface's identity."""
    from climate.gcm.turned import read_frame
    h16, ident = heights()
    surf, sea, lake, level = water_surface(h16)
    frame = read_frame(ROOT / INPUTS['orientation'][0])
    reach = REACH_PX[-1] * DX_KM / (MOON_RADIUS / 1000.0) * 180.0 / math.pi + 0.5
    grid = stereo_grid(cap, reach)
    surf_km = sample_today(grid, surf, frame['rotation'], order=1) / 1000.0
    sea_g = sample_today(grid, sea.astype(np.float32), frame['rotation'], order=0) > 0.5
    lake_g = sample_today(grid, lake.astype(np.float32), frame['rotation'], order=0) > 0.5
    colat = 90.0 - grid['sgn'] * grid['lat']
    mask = colat <= CAP_DEG
    key = hashlib.sha256(np.ascontiguousarray(surf_km.astype(np.float32)).tobytes() +
                         json.dumps([REACH_PX.tolist(), N_AZ, DX_KM, CAP_DEG]).encode()).hexdigest()[:16]
    cache = Path(str(HORIZON_CACHE).format(cap=cap))
    idx = np.flatnonzero(mask.ravel())
    if cache.exists() and str(np.load(cache)['key']) == key:
        hz = np.load(cache)['horizon'].astype(np.float32)
        n0 = grid['nvec'].reshape(-1, 3)[idx]
        north, east = local_axes(n0)
        log(f'{cap}: horizons from the cache {cache.name}')
    else:
        t = time.time()
        hz, n0, north, east = horizons(surf_km, grid, mask, workers=workers)
        np.savez_compressed(cache, horizon=hz.astype(np.float16), key=key)
        log(f'{cap}: horizons for {idx.size} pixels in {time.time() - t:.0f} s')
    return dict(grid=grid, surf_km=surf_km, sea=sea_g.ravel()[idx], lake=lake_g.ravel()[idx], idx=idx, hz=hz, n0=n0,
                north=north, east=east, lat=grid['lat'].ravel()[idx], area=grid['area'].ravel()[idx],
                height_km=surf_km.ravel()[idx], key=key, cache=str(cache), sea_level_m=level, heights=ident,
                frame_sha256=frame['sha256'])


def band_share(values, weight, lat, cond, bands=BANDS):
    """Area share of each |latitude| band (and the whole cap) where cond holds."""
    out = {}
    alat = np.abs(lat)
    for lo, hi in bands:
        sel = (alat >= lo) & (alat < hi)
        out[f'{lo:g}-{min(hi, 90):g}'] = float(weight[sel & cond].sum() / weight[sel].sum()) if sel.any() else None
    out['60-90'] = float(weight[cond].sum() / weight.sum())
    return out


def band_quantiles(values, weight, lat, cond, qs=(0.1, 0.5, 0.9), bands=BANDS):
    out = {}
    alat = np.abs(lat)
    for lo, hi in list(bands) + [(60.0, 90.01)]:
        sel = (alat >= lo) & (alat < hi) & cond
        if sel.sum() < 10:
            out[f'{lo:g}-{min(hi, 90):g}'] = None
            continue
        v, w = values[sel], weight[sel]
        o = np.argsort(v)
        cw = np.cumsum(w[o]) / w.sum()
        out[f'{lo:g}-{min(hi, 90):g}'] = [float(np.interp(q, cw, v[o])) for q in qs]
    return out


def terrain(light, cc, refr, workers=4, log=print):
    """Direct sunlight, permanent shade and the light the horizon admits on each cap's ground at 16 px/deg."""
    tab = sky_table()
    out = {}
    for cap in ('north', 'south'):
        c = cap_surface(cap, workers=workers, log=log)
        sgn = 1.0 if cap == 'north' else -1.0
        land = ~(c['sea'] | c['lake'])
        w = c['area']
        res = dict(pixels=int(c['idx'].size), area_km2=float(w.sum()), land_share=float(w[land].sum() / w.sum()),
                   sea_share=float(w[c['sea']].sum() / w.sum()), lake_share=float(w[c['lake']].sum() / w.sum()),
                   height_km_quantiles=band_quantiles(c['height_km'], w, c['lat'], land),
                   horizon_cache=dict(path=c['cache'], key=c['key']))
        hmax = c['hz'].max(axis=1)
        res['highest_horizon_deg_quantiles'] = band_quantiles(hmax, w, c['lat'], np.ones_like(land))
        svf = 1.0 - (np.sin(np.radians(np.maximum(c['hz'], 0.0))) ** 2).mean(axis=1)
        res['sky_view_factor_quantiles'] = band_quantiles(svf, w, c['lat'], np.ones_like(land))
        t = time.time()
        decs = [sgn * d for d in (-2.5, -TILTS[0], 0.0, TILTS[0], 2.5)]            # winter to summer for this cap
        names = ('winter, tilt 2.5', 'winter, tilt 1.613', 'equinox', 'summer, tilt 1.613', 'summer, tilt 2.5')
        lit, flat = sunlit_hours(c['hz'], c['n0'], c['north'], c['east'], decs)
        lit_r, _ = sunlit_hours(c['hz'], c['n0'], c['north'], c['east'], [decs[3], decs[4]], refr=refr)
        log(f'{cap}: direct sunlight in {time.time() - t:.0f} s')
        sun = {}
        for k, name in enumerate(names):
            lost = 1.0 - lit[k] / np.maximum(flat[k], 1e-9)
            sun[name] = dict(
                declination_deg=decs[k],
                no_direct_sun_share=band_share(None, w, c['lat'], lit[k] == 0),
                no_direct_sun_share_land=band_share(None, w[land], c['lat'][land], lit[k][land] == 0),
                flat_ground_without_sun_share=band_share(None, w, c['lat'], flat[k] == 0),
                hours_lit_quantiles=band_quantiles(lit[k], w, c['lat'], np.ones_like(land)),
                share_of_flat_hours_lost_mean={b: v for b, v in zip(
                    [f'{lo:g}-{min(hi, 90):g}' for lo, hi in BANDS] + ['60-90'],
                    [float((lost * w)[(np.abs(c['lat']) >= lo) & (np.abs(c['lat']) < hi) & (flat[k] > 0)].sum() /
                           max(w[(np.abs(c['lat']) >= lo) & (np.abs(c['lat']) < hi) & (flat[k] > 0)].sum(), 1e-9))
                     for lo, hi in list(BANDS) + [(60.0, 90.01)]])})
        res['direct_sun'] = sun
        res['no_direct_sun_all_year'] = dict(
            tilt_1613=band_share(None, w, c['lat'], lit[3] == 0), tilt_2500=band_share(None, w, c['lat'], lit[4] == 0),
            tilt_1613_refracted=band_share(None, w, c['lat'], lit_r[0] == 0),
            tilt_2500_refracted=band_share(None, w, c['lat'], lit_r[1] == 0),
            tilt_1613_area_km2=float(w[lit[3] == 0].sum()), tilt_2500_area_km2=float(w[lit[4] == 0].sum()),
            tilt_1613_land_area_km2=float(w[(lit[3] == 0) & land].sum()))
        # the light the horizon admits, on every second pixel in each direction
        gi = c['idx'] // c['grid']['nvec'].shape[1]
        gj = c['idx'] % c['grid']['nvec'].shape[1]
        sub = (gi % TERRAIN_SUBSAMPLE == 0) & (gj % TERRAIN_SUBSAMPLE == 0)
        t = time.time()
        tl = terrain_light(light, tab, c['hz'][sub], c['n0'][sub], c['north'][sub], c['east'][sub],
                           [decs[1], decs[2], decs[3]], cc, c['lat'][sub])
        log(f'{cap}: terrain light on {int(sub.sum())} pixels in {time.time() - t:.0f} s')
        ws, ls, lands = w[sub], c['lat'][sub], land[sub]
        shaded_eq = (lit[2] == 0)[sub]
        light_out = {}
        for name, r in zip(('winter', 'equinox', 'summer'), tl):
            light_out[name] = dict(
                declination_deg=r['dec'],
                dli_clear_quantiles=band_quantiles(r['mol_day_clear'], ws, ls, lands),
                dli_cloud_quantiles=band_quantiles(r['mol_day_cloud'], ws, ls, lands),
                dark_hours_clear_quantiles=band_quantiles(r['dark_clear'], ws, ls, lands),
                dark_hours_cloud_quantiles=band_quantiles(r['dark_cloud'], ws, ls, lands),
                share_with_dark_hours_clear=band_share(None, ws[lands], ls[lands], r['dark_clear'][lands] > 0.5),
                share_with_dark_hours_cloud=band_share(None, ws[lands], ls[lands], r['dark_cloud'][lands] > 0.5),
                hours_below_reading_quantiles=band_quantiles(r['below_reading'], ws, ls, lands),
                hours_below_work_quantiles=band_quantiles(r['below_work'], ws, ls, lands),
                darkest_lux_quantiles=band_quantiles(r['darkest_lux'], ws, ls, lands),
                blocked_share_at_darkest_quantiles=band_quantiles(r['blocked_at_darkest'], ws, ls, lands),
                shaded_ground=dict(
                    note='pixels without direct sunlight through the equinox lunar day',
                    dli_clear_quantiles=band_quantiles(r['mol_day_clear'], ws, ls, shaded_eq),
                    dark_hours_clear_quantiles=band_quantiles(r['dark_clear'], ws, ls, shaded_eq),
                    darkest_lux_quantiles=band_quantiles(r['darkest_lux'], ws, ls, shaded_eq)))
        res['light'] = light_out
        res['light_pixels'] = int(sub.sum())
        res['maps'] = dict(lit_equinox=lit[2], flat_equinox=flat[2], lit_summer=lit[3], hz_max=hmax, svf=svf,
                           sub=sub, dark_cloud_winter=tl[0]['dark_cloud'], dli_clear_equinox=tl[1]['mol_day_clear'],
                           grid_n=c['grid']['n'], idx=c['idx'], land=land)
        res['frame_sha256'] = c['frame_sha256']
        res['heights_identity'] = c['heights']
        res['sea_level_m'] = c['sea_level_m']
        out[cap] = res
    return out


# ==== run-sunlight extraction (needs netCDF4) =======================================================================
def extract_run_sunlight(path=RUN_EXTRACT, years=RUN_YEARS):
    """The run's surface sunlight, all-sky (rss) and clear-sky (rsscs), its cloud cover (clt) and mean solar height
    (czen), summed by latitude row and tenth of the lunar day over every longitude and over land, for years 10-49.
    Needs netCDF4 (climate/gcm/.venv); reads the run's output and writes nothing but the extract."""
    import netCDF4
    from climate.crm import gcm_site as gs
    names = ('rss', 'rsscs', 'clt', 'czen')
    sums, digests = None, {}
    for y in range(years[0], years[1] + 1):
        f = RUN_DIR / f'MOST.{y:05d}.nc'
        with netCDF4.Dataset(f) as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            land = np.asarray(d['lsm'][:], float)[0] > 0.5
            fields = {k: np.asarray(d[k][:], float) for k in names}
        h = hashlib.sha256()
        for a in (lat, lon, *(fields[k] for k in names)):
            h.update(np.ascontiguousarray(a).tobytes())
        digests[f.name] = dict(fields_read_sha256=h.hexdigest()[:32], bytes=f.stat().st_size)
        hour = (lon[None, :] - gs.subsolar_longitude(fields['czen'], lat, lon)[:, None] + 180.0) % 360.0 - 180.0
        tenth = np.floor((hour + 180.0) / 36.0).astype(int) % 10
        if sums is None:
            sums = {f'{k}_{part}': np.zeros((10, lat.size)) for k in names + ('n',) for part in ('all', 'land')}
        for k in range(10):
            sel = (tenth == k)[:, None, :]
            for part, mask in (('all', sel), ('land', sel & land[None, :, :])):
                m = np.broadcast_to(mask, fields['rss'].shape)
                for name in names:
                    sums[f'{name}_{part}'][k] += np.where(m, fields[name], 0.0).sum(axis=(0, 2))
                sums[f'n_{part}'][k] += m.sum(axis=(0, 2))
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, lat=lat, tenths=TENTHS, **sums, digests=json.dumps(digests), run=RUN, years=np.array(years))
    return path


# ==== assembly ======================================================================================================
def flat_ground(light, per, cc, refr0, log=print):
    """Every latitude's 72 seasonal lunar days for each tilt: the Sun alone at the run's precession phase, at the
    north's and at the mean distance, and with the Earth in each sector at the run's phase."""
    out = {}
    for eps in TILTS:
        rows = {}
        for hemi in (1.0, -1.0):
            for lat0 in LATS:
                lat = hemi * float(lat0)
                r = {}
                for name, psi in (('run', PSI_RUN), ('north_perihelion', PSI_NORTH), ('mean_distance', None)):
                    r[name] = summarise_days(site_days(light, per, lat, eps, psi=psi, cc=cc), STEP_H,
                                             light.sun_radius, refr0)
                for sector, lon in SECTORS.items():
                    r[sector] = summarise_days(site_days(light, per, lat, eps, sector_lon=lon, psi=PSI_RUN, cc=cc),
                                               STEP_H, light.sun_radius, refr0)
                rows[lat] = r
        out[eps] = rows
        log(f'flat ground, tilt {eps}: {len(rows)} latitudes')
    return out


def annual(values, weights):
    return float(np.sum(values * weights))


def season_index(phases=SEASON_TABLE):
    return [int(np.argmin(np.abs(SEASONS - p))) for p in phases]


def derive(flat, per, refr, light):
    """The product's flat-ground sections from the per-day summaries."""
    P = per['synodic_days']
    iseason = season_index()
    w_run = season_weights(SEASONS, PSI_RUN, per)
    w_north = season_weights(SEASONS, PSI_NORTH, per)
    sun_path, light_tab, annual_tab, earth_tab = {}, {}, {}, {}
    for eps, rows in flat.items():
        key = f'tilt {eps:g}'
        sp, lt, an, et = {}, {}, {}, {}
        for lat, r in rows.items():
            m = r['run']
            sp[f'{lat:g}'] = dict(
                hours_centre=m['hours_centre'][iseason], hours_limb=m['hours_limb'][iseason],
                hours_limb_refracted=m['hours_limb_refracted'][iseason],
                refraction_adds_hours=(m['hours_limb_refracted'] - m['hours_limb'])[iseason],
                max_elevation_deg=m['max_elevation'][iseason], min_elevation_deg=m['min_elevation'][iseason],
                hours_centre_range=[m['hours_centre'].min(), m['hours_centre'].max()],
                max_elevation_range=[m['max_elevation'].min(), m['max_elevation'].max()])
            near = r['near side']
            an[f'{lat:g}'] = dict(
                dli_clear=annual(m['mol_clear'], w_run) / P, dli_cloud=annual(m['mol_cloud'], w_run) / P,
                dli_overcast=annual(m['mol_overcast'], w_run) / P,
                dli_clear_range=[m['mol_clear'].min() / P, m['mol_clear'].max() / P],
                dli_cloud_range=[m['mol_cloud'].min() / P, m['mol_cloud'].max() / P],
                dli24_lowest=m['dli24_min'].min(), dli24_highest=m['dli24_max'].max(),
                dark_hours_clear_mean=annual(m['dark_clear'], w_run), dark_hours_clear_max=m['dark_clear'].max(),
                dark_hours_cloud_mean=annual(m['dark_cloud'], w_run), dark_hours_cloud_max=m['dark_cloud'].max(),
                dark_hours_overcast_max=m['dark_overcast'].max(),
                lunar_days_with_dark_share_clear=float(np.sum(w_run * (m['dark_clear'] > 0))),
                lunar_days_with_dark_share_cloud=float(np.sum(w_run * (m['dark_cloud'] > 0.5))),
                darkest_lux_clear=m['darkest_lux'].min(), darkest_ppfd_clear=m['darkest_ppfd'].min(),
                darkest_lux_near_side=near['darkest_lux'].min(),
                hours_below={name: dict(clear_mean=annual(m[f'below {name}'], w_run),
                                        clear_max=m[f'below {name}'].max(),
                                        cloud_mean=annual(m[f'below {name} cloud'], w_run),
                                        cloud_max=m[f'below {name} cloud'].max()) for name in PEOPLE_LUX},
                hours_above_ek=dict(low=annual(m['above_ek_low'], w_run), high=annual(m['above_ek_high'], w_run)),
                dli_clear_north_perihelion=annual(r['north_perihelion']['mol_clear'], w_north) / P,
                dli_clear_mean_distance=float(np.mean(r['mean_distance']['mol_clear'])) / P)
            if abs(lat) in TABLE_LATS:
                lt[f'{lat:g}'] = dict(
                    dli_clear=m['mol_clear'][iseason] / P, dli_cloud=m['mol_cloud'][iseason] / P,
                    dli24_min=m['dli24_min'][iseason], dli24_max=m['dli24_max'][iseason],
                    dark_hours_clear=m['dark_clear'][iseason], dark_hours_cloud=m['dark_cloud'][iseason],
                    dark_hours_overcast=m['dark_overcast'][iseason],
                    darkest_lux=m['darkest_lux'][iseason], darkest_ppfd=m['darkest_ppfd'][iseason],
                    brightest_lux=m['brightest_lux'][iseason],
                    mean_lux_clear=m['mean_lux_clear'][iseason], mean_lux_cloud=m['mean_lux_cloud'][iseason],
                    hours_below_reading_cloud=m['below reading cloud'][iseason],
                    hours_below_close_work_cloud=m['below close work outdoors cloud'][iseason],
                    hours_below_indoor_reading_cloud=m['below indoor reading and writing cloud'][iseason],
                    hours_above_ek_low=m['above_ek_low'][iseason], hours_above_ek_high=m['above_ek_high'][iseason],
                    dli_clear_north_perihelion=r['north_perihelion']['mol_clear'][iseason] / P,
                    mean_cover=m['mean_cover'][iseason])
                et[f'{lat:g}'] = {sector: dict(
                    earth_max_lux=r[sector]['earth_max_lux'][iseason],
                    earth_share_at_darkest=r[sector]['earth_share_at_darkest'][iseason],
                    darkest_lux=r[sector]['darkest_lux'][iseason],
                    earth_mol_per_lunar_day=r[sector]['mol_earth'][iseason]) for sector in SECTORS}
        sun_path[key], light_tab[key], annual_tab[key], earth_tab[key] = sp, lt, an, et
    return sun_path, light_tab, annual_tab, earth_tab


def first_in(rows, lats, hemi, case, cond):
    """The lowest |latitude| from which cond holds for the given case's lunar days."""
    ok = [l for l in lats if cond(rows[hemi * l][case])]
    return min(ok) if ok else None


def thresholds(flat, per):
    """Latitudes where the polar day and night begin and where the plant model's dark hours end, per tilt."""
    out = {}
    full = per['synodic_hours'] - STEP_H
    for eps, rows in flat.items():
        o = {}
        for hemi, name in ((1.0, 'north'), (-1.0, 'south')):
            lats = sorted(abs(l) for l in rows if np.sign(l) == hemi)
            def first(cond):
                ok = [l for l in lats if cond(rows[hemi * l]['run'])]
                return min(ok) if ok else None
            o[name] = dict(
                polar_day_centre_from=first(lambda m: m['hours_centre'].max() >= full),
                polar_day_limb_refracted_from=first(lambda m: m['hours_limb_refracted'].max() >= full),
                polar_night_centre_from=first(lambda m: m['hours_centre'].min() <= 0),
                polar_night_limb_refracted_from=first(lambda m: m['hours_limb_refracted'].min() <= 0),
                no_dark_hours_clear_from=first(lambda m: m['dark_clear'].max() <= 0),
                no_dark_hours_cloud_from=first(lambda m: m['dark_cloud'].max() <= 0.5),
                no_dark_hours_overcast_from=first(lambda m: m['dark_overcast'].max() <= 0),
                walking_light_always_from=first(lambda m: m['below walking cloud'].max() <= 0.5),
                reading_light_always_clear_from=first(lambda m: m['below reading'].max() <= 0),
                reading_light_always_cloud_from=first(lambda m: m['below reading cloud'].max() <= 0.5),
                no_dark_hours_clear_from_north_perihelion=first_in(rows, lats, hemi, 'north_perihelion',
                                                                   lambda m: m['dark_clear'].max() <= 0),
                no_dark_hours_cloud_from_north_perihelion=first_in(rows, lats, hemi, 'north_perihelion',
                                                                   lambda m: m['dark_cloud'].max() <= 0.5))
        out[f'tilt {eps:g}'] = o
    return out


def toa_swing(per, eps, psi, band=(75.0, 90.0)):
    """The seasonal swing of the lunar-day mean sunlight above the air over a polar band (W/m2), each hemisphere."""
    from shared.constants import SOLAR_CONSTANT
    u = np.arange(0.0, 360.0, 1.0)
    dec = np.radians(declination(u, eps))
    f = distance_factor(u, psi, per['eccentricity'])
    out = {}
    for hemi, name in ((1.0, 'north'), (-1.0, 'south')):
        lats = np.radians(np.linspace(band[0], band[1], 61) * hemi)
        wts = np.cos(lats)
        H0 = np.arccos(np.clip(-np.tan(lats)[:, None] * np.tan(dec)[None, :], -1, 1))
        q = SOLAR_CONSTANT * f / math.pi * (H0 * np.sin(lats)[:, None] * np.sin(dec) +
                                             np.cos(lats)[:, None] * np.cos(dec) * np.sin(H0))
        qb = (q * wts[:, None]).sum(0) / wts.sum()
        out[name] = dict(mean=float(qb.mean()), swing=float(qb.max() - qb.min()))
    return out


def tilt_and_precession(annual_tab, light_tab, sun_path, thr, per):
    """What moves between the tilts and with the precession's phase, as differences of the computed numbers."""
    out = dict(by_latitude={})
    a1, a2 = annual_tab['tilt 1.613'], annual_tab['tilt 2.5']
    l1, l2 = light_tab['tilt 1.613'], light_tab['tilt 2.5']
    s1, s2 = sun_path['tilt 1.613'], sun_path['tilt 2.5']
    iu90, iu270 = 3, 9                                     # u = 90 (north's midsummer) and 270 in SEASON_TABLE
    for lat in l1:
        hemi = 1.0 if float(lat) > 0 else -1.0
        summer, winter = (iu90, iu270) if hemi > 0 else (iu270, iu90)
        out['by_latitude'][lat] = dict(
            summer_max_elevation=[s1[lat]['max_elevation_deg'][summer], s2[lat]['max_elevation_deg'][summer]],
            summer_hours_centre=[s1[lat]['hours_centre'][summer], s2[lat]['hours_centre'][summer]],
            winter_hours_centre=[s1[lat]['hours_centre'][winter], s2[lat]['hours_centre'][winter]],
            summer_dli_clear=[l1[lat]['dli_clear'][summer], l2[lat]['dli_clear'][summer]],
            winter_dli_clear=[l1[lat]['dli_clear'][winter], l2[lat]['dli_clear'][winter]],
            winter_dark_hours_cloud=[l1[lat]['dark_hours_cloud'][winter], l2[lat]['dark_hours_cloud'][winter]],
            annual_dli_clear=[a1[lat]['dli_clear'], a2[lat]['dli_clear']],
            summer_dli_clear_precession=[l1[lat]['dli_clear'][summer], l1[lat]['dli_clear_north_perihelion'][summer]],
            winter_dli_clear_precession=[l1[lat]['dli_clear'][winter], l1[lat]['dli_clear_north_perihelion'][winter]])
    out['toa_polar_band_75_90'] = {
        'tilt 1.613, run phase (south summer at perihelion)': toa_swing(per, 1.613, PSI_RUN),
        'tilt 1.613, north summer at perihelion': toa_swing(per, 1.613, PSI_NORTH),
        'tilt 2.5, run phase': toa_swing(per, 2.5, PSI_RUN),
        'tilt 1.613, circular orbit': toa_swing(dict(per, eccentricity=0.0), 1.613, PSI_RUN)}
    out['thresholds'] = thr
    return out


def overcast_sensitivity(light, per, cc, refr0):
    """The flat ground's dark hours and photons under the run's clouds with an overcast passing 0.1 or 0.5 of the light
    (tilt 1.613, the run's phase), beside the central 0.25."""
    P = per['synodic_days']
    w = season_weights(SEASONS, PSI_RUN, per)
    out = {}
    for ov in OVERCAST_RANGE:
        first, rows = {}, {}
        for hemi, name in ((1.0, 'north'), (-1.0, 'south')):
            ok = []
            for lat0 in LATS:
                lat = hemi * float(lat0)
                m = summarise_days(site_days(light, per, lat, TILTS[0], psi=PSI_RUN, cc=cc, overcast=ov), STEP_H,
                                   light.sun_radius, refr0)
                if m['dark_cloud'].max() <= 0.5:
                    ok.append(float(lat0))
                if float(lat0) in TABLE_LATS:
                    rows[f'{lat:g}'] = dict(dli_cloud=annual(m['mol_cloud'], w) / P,
                                            dark_hours_cloud_mean=annual(m['dark_cloud'], w),
                                            dark_hours_cloud_max=m['dark_cloud'].max())
            first[name] = min(ok) if ok else None
        out[f'overcast {ov:g}'] = dict(no_dark_hours_cloud_from=first, by_latitude=rows)
    return out


def compact_terrain(t):
    """The terrain section without its maps."""
    return {cap: {k: v for k, v in r.items() if k != 'maps'} for cap, r in t.items()}


def figures(product, flat, terr, light, per):
    """Figures on the research drive (out of Git)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    DRIVE.mkdir(parents=True, exist_ok=True)
    written = []
    # 1. hours of Sun a lunar day and its highest elevation against latitude, by season and tilt
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    seasons = ((90.0, 'north midsummer'), (0.0, 'equinox'), (270.0, 'north midwinter'))
    styles = {TILTS[0]: '-', TILTS[1]: '--'}
    colours = ('#c0392b', '#7f8c8d', '#2c3e80')
    for eps in TILTS:
        rows = flat[eps]
        lats = np.array(sorted(l for l in rows if l > 0))
        for (u0, name), col in zip(seasons, colours):
            k = int(np.argmin(np.abs(SEASONS - u0)))
            hrs = [rows[l]['run']['hours_limb_refracted'][k] for l in lats]
            top = [rows[l]['run']['max_elevation'][k] for l in lats]
            axes[0].plot(lats, hrs, styles[eps], color=col, label=f'{name}, tilt {eps:g} deg')
            axes[1].plot(lats, top, styles[eps], color=col, label=f'{name}, tilt {eps:g} deg')
    axes[0].set_ylabel('hours of Sun a lunar day (upper limb, refracted; of 708.7)')
    axes[1].set_ylabel('highest Sun elevation in the lunar day (deg)')
    for ax in axes:
        ax.set_xlabel('latitude (deg N; the south mirrors it half an eclipse year later)')
        ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.suptitle('The Sun\'s path over the turned Moon\'s north polar region')
    p = DRIVE / 'sun_path_north.png'
    fig.savefig(p, dpi=110)
    plt.close(fig)
    written.append(p)
    # 2. light through a lunar day at four latitudes, clear, north midwinter and midsummer, tilt 1.613
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True, sharey=True)
    for ax, (u0, name) in zip(axes, ((270.0, 'north midwinter (u = 270)'), (90.0, 'north midsummer (u = 90)'))):
        for lat in (60.0, 65.0, 70.0, 80.0, 90.0):
            x = site_days(light, per, lat, 1.613, psi=PSI_RUN, cc=None)
            k = int(np.argmin(np.abs(SEASONS - u0)))
            ax.semilogy(x['s']['t_h'], x['ppfd'][k], label=f'{lat:g} deg')
        ax.axhline(DARK_UMOL, color='k', lw=0.8, ls='--')
        ax.text(-350, DARK_UMOL * 1.2, 'the plant model\'s dark, 1 umol/m2/s', fontsize=8)
        ax.set_xlabel('hours from local noon')
        ax.set_title(name)
        ax.set_ylim(0.05, 2000)
    axes[0].set_ylabel('photosynthetic photons on level ground, clear sky (umol/m2/s)')
    axes[1].legend(title='latitude N')
    p = DRIVE / 'lunar_day_light.png'
    fig.savefig(p, dpi=110)
    plt.close(fig)
    written.append(p)
    # 3. terrain maps per cap: equinox sunlit share, sky view, winter dark hours under clouds
    for cap, r in terr.items():
        mp = r['maps']
        n = mp['grid_n']
        size = 2 * n + 1
        def canvas(values, idx, step=1):
            img = np.full(size * size, np.nan)
            img[idx] = values
            return img.reshape(size, size)[::step, ::step]
        ext = [-n - 0.5, n + 0.5, -n - 0.5, n + 0.5]
        fig, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
        share = mp['lit_equinox'] / np.maximum(mp['flat_equinox'], 1e-9)
        im = axes[0].imshow(canvas(share, mp['idx']), origin='lower', cmap='cividis', vmin=0, vmax=1, extent=ext)
        fig.colorbar(im, ax=axes[0], label='share of flat ground\'s Sun hours, equinox lunar day')
        im = axes[1].imshow(canvas(mp['hz_max'], mp['idx']), origin='lower', cmap='magma', vmin=0, vmax=20, extent=ext)
        fig.colorbar(im, ax=axes[1], label='highest horizon (deg)')
        sub_idx = mp['idx'][mp['sub']]
        im = axes[2].imshow(canvas(mp['dark_cloud_winter'], sub_idx, TERRAIN_SUBSAMPLE), origin='lower', cmap='viridis',
                             vmin=0, vmax=150, extent=ext)
        fig.colorbar(im, ax=axes[2], label='hours under 1 umol/m2/s, midwinter lunar day, the run\'s clouds')
        for ax in axes:
            ax.set_xticks([])
            ax.set_yticks([])
            for c_deg in (10, 20, 30):
                rho = 2 * MOON_RADIUS / 1000 * math.tan(math.radians(c_deg) / 2) / DX_KM
                ax.add_patch(plt.Circle((0, 0), rho, fill=False, color='w', lw=0.6, ls=':'))
        fig.suptitle(f'Turned {cap} polar cap, 60 deg to the pole (dotted: 80, 70, 60 deg), polar stereographic, 1.9 km')
        p = DRIVE / f'terrain_{cap}.png'
        fig.savefig(p, dpi=100)
        plt.close(fig)
        written.append(p)
    return [dict(path=str(p), sha256=sha(p)) for p in written]


EVIDENCE = ('A calculation on committed products: a mean-geometry ephemeris of the turned Moon in Cassini state 2; the '
            'light calendar\'s clear-sky evaluator on the solved 1.2-atm column; photons per lux from the solved '
            'standard solution\'s channels; the run\'s cloud cover in a two-state sky with an Earth overcast '
            'transmission; LOLA terrain with the atlas\'s seas and the routed lakes. Not observations, not a climate '
            'run, not validated against measured light on any Moon; the twilight and clear sky are the solved '
            'molecular column\'s (no aerosol, no refraction in the scattering).')
READING_RULE = ('Latitudes are the turned Moon\'s (north positive). The eclipse-year phase u is 0 at the Sun\'s northward '
                'crossing of the lunar equator, 90 at the north\'s midsummer and 270 at the south\'s; season tables hold '
                'u = 0, 30, ..., 330 at the lunar day\'s local noon. A lunar day is 708.7 hours, midnight to midnight. '
                'dli values are mol/m2 per 24 hours averaged over the lunar day (the lunar day\'s photons / 29.53); '
                'dli24 the lowest and highest 24-hour window inside it. Hours are per lunar day. "cloud" is the '
                'expected value in the two-state sky with the run\'s cover; "overcast" the whole day under cloud. '
                'Light is on level, open ground unless a terrain section says otherwise; the run\'s precession phase '
                '(south summer at perihelion) unless a key names another. Annual values weight each season by the '
                'time the Sun spends there. Pairs in tilt_and_precession are [tilt 1.613, tilt 2.5] or [run phase, '
                'north summer at perihelion]. Terrain quantiles are [10th, 50th, 90th] percentiles by area.')
ASSUMPTIONS = [
    'Mean geometry: uniform rotation, the Sun\'s hour angle at the mean synodic rate, the Sun\'s declination from the '
    'eclipse-year phase, the Earth\'s orbit with the run\'s year and eccentricity; a mean Earth distance and sinusoidal '
    'libration in longitude; good to a few tenths of a degree, not an ephemeris for dates',
    'Clear-sky light is the calendar evaluator\'s: one horizontally uniform 1.2-atm molecular column, Lambertian ground '
    'of albedo 0.1, no aerosol; refraction is computed apart and enters only the hours the Sun is up',
    'Clouds: the run\'s cover by latitude and hour angle in a two-state sky; overcast passes 0.25 of the clear sky\'s '
    'light at every Sun height and in twilight (an Earth relation); the run\'s own surface sunlight cannot serve as a '
    'cloud transmission in the polar rows (see clouds.run_sunlight_diagnostic)',
    'Terrain: horizons on a 1.9 km grid to 265 km, nearest-neighbour samples; light on horizontal ground, no light '
    'reflected from the terrain itself, no sloped surfaces; seas at the atlas\'s level and the routed lakes as flat water',
    'Earthlight carries sunlight\'s photons per lux for the sky at the Earth\'s height; at most a few lux at the poles',
]


def build(workers=4, with_figures=False, log=print):
    t0 = time.time()
    per = periods()
    transfer = read_json('transfer')
    light = Light(transfer)
    cc = cloud_cover()
    ref = refraction_summary()
    refr0 = float(ref['moon']['refraction_deg'][0])
    refr_curve = (np.asarray(ref['moon']['apparent_elevation_deg']), np.asarray(ref['moon']['refraction_deg']))
    log(f'tables in {time.time() - t0:.0f} s')
    flat = flat_ground(light, per, cc, refr0, log=log)
    sun_path, light_tab, annual_tab, earth_tab = derive(flat, per, ref, light)
    thr = thresholds(flat, per)
    tp = tilt_and_precession(annual_tab, light_tab, sun_path, thr, per)
    sens = overcast_sensitivity(light, per, cc, refr0)
    log(f'flat ground done at {time.time() - t0:.0f} s')
    terr = terrain(light, cc, refr_curve, workers=workers, log=log)
    log(f'terrain done at {time.time() - t0:.0f} s')
    plant = read_json('plant')
    ll = read_json('light_limitation')
    living = read_json('living_light')
    solved = read_json('solved_sky')
    atlas_head = json.loads(str(np.load(ATLAS)['meta']))
    checks = dict(
        solution_is_atlas_source=sha(SOLUTION) == atlas_head['source_solution_sha256'],
        atlas_is_solved_sky_atlas=sha(ATLAS) == solved['worlds']['moon_1.2atm']['atlas']['sha256'],
        direct_normal_lux_overhead=dict(
            channels=float(np.interp(90.0, light.photons['direct_elevation_deg'], light.photons['direct_normal_lux'])),
            solved_sky=next(s['direct_normal_lux'] for s in solved['worlds']['moon_1.2atm']['samples']
                            if s['sun_deg'] == 90.0)),
        diffuse_lux_at_minus_10=dict(
            channels=float(np.interp(-10.0, light.photons['diffuse_elevation_deg'], light.photons['diffuse_lux'])),
            calendar=float(np.interp(-10.0, light.elev, light.diffuse))),
        equinox_hours_half_lunar_day=float(flat[TILTS[0]][70.0]['mean_distance']['hours_centre'][0]),
        refraction_earth_control_deg=float(ref['earth']['refraction_deg'][0]),
        plant_model_twilight_par_at_minus_25=dict(
            this_study=float(light.sun(np.array([-25.0]))[2][0]),
            plant_model=float(np.interp(-25.0, plant['twilight']['moon']['elevation_deg'],
                                        plant['twilight']['moon']['par_umol_m2_s']))))
    files = [Path(__file__)] + [ROOT / f for f in IMPORTED]
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='research', study='research/studies/polar_light', runner='research/studies/polar_light/run.py',
                      files={str(f.relative_to(ROOT)): sha(f) for f in files},
                      inputs={INPUTS[k][0]: sha(ROOT / INPUTS[k][0]) for k in INPUTS},
                      drive_inputs={str(SOLUTION): sha(SOLUTION), str(ATLAS): sha(ATLAS)},
                      constants=constants_record(), software=dict(python=platform.python_version(),
                                                                  numpy=np.__version__),
                      run=RUN, climatology_years=list(RUN_YEARS)),
        evidence=EVIDENCE, reading_rule=READING_RULE, assumptions=ASSUMPTIONS, sources=SOURCES,
        settings=dict(tilts_deg=TILTS, precession_phase_run_deg=PSI_RUN, precession_phase_north_deg=PSI_NORTH,
                      latitudes=LATS, season_phases_deg=SEASON_TABLE, step_hours=STEP_H, dark_umol=DARK_UMOL,
                      overcast_transmission=OVERCAST, people_lux=PEOPLE_LUX, ek_umol=EK_UMOL,
                      vgpm_half_mol=VGPM_HALF_MOL, dli_needs=DLI_NEEDS, faust_range=FAUST_RANGE,
                      sectors_turned_longitude=SECTORS, periods=per,
                      sun_semidiameter_deg=light.sun_radius,
                      earth_latitude_libration_deg={f'tilt {e:g}': MOON_ORBIT_INCLINATION_DEG + e for e in TILTS},
                      earth_longitude_libration_deg=LIBRATION_LONGITUDE_DEG),
        refraction=ref,
        photons_per_lux=dict(
            elevation_deg=[-30, -25, -20, -10, -5, 0, 2, 5, 10, 20, 30, 60, 90],
            diffuse=np.interp([-30, -25, -20, -10, -5, 0, 2, 5, 10, 20, 30, 60, 90],
                              light.photons['diffuse_elevation_deg'], light.photons['diffuse_umol_per_lux']),
            direct=np.where(np.array([-30, -25, -20, -10, -5, 0, 2, 5, 10, 20, 30, 60, 90]) >= 0, np.interp(
                [-30, -25, -20, -10, -5, 0, 2, 5, 10, 20, 30, 60, 90], light.photons['direct_elevation_deg'],
                light.photons['direct_umol_per_lux']), np.nan),
            channels=light.photons['channels'], par_channels=light.photons['par_channels']),
        clear_sky_by_sun_elevation=clear_sky_table(light),
        sun_path=sun_path, light=light_tab, annual=annual_tab, earth=earth_tab, thresholds=thr,
        tilt_and_precession=tp,
        clouds=dict(cover_by_row={f'{la:.2f}': dict(mean=float(cc['clt'][i].mean()), by_hour_angle=cc['clt'][i][::3])
                                  for i, la in enumerate(cc['lat']) if abs(la) > 55},
                    hour_angle_deg=cc['hour'][::3], overcast_sensitivity=sens,
                    run_sunlight_diagnostic=run_sunlight_diagnostic()),
        terrain=compact_terrain(terr),
        comparisons=dict(
            plant_model_dark_hours=plant['results']['near']['60']['earth_like']['dark_hours'],
            plant_model_bands=plant['settings']['bands_abs_latitude'],
            moisture_map_rule='plant dark hours fall to none at 65 deg (biosphere/ecology/moisture_map.py, dark_hours)',
            light_limitation_ek_umol=ll['guesses']['ek_umol_m2_s'],
            living_light_polar_plain_75N_darkest_lux=living['needs']['places'].get('polar_plain', {}).get('darkest_lux')),
        checks=checks, seconds=round(time.time() - t0, 1))
    if with_figures:
        product['figures'] = figures(product, flat, terr, light, per)
    out = plain(product)
    out['producer']['constants'] = constants_record()          # exact, for the constants test
    return out


def clear_sky_table(light, elev=(-35, -30, -25, -20, -15, -10, -5, -2, 0, 1, 2, 3, 5, 7, 10, 15, 20, 30)):
    """Clear-sky light on open level ground by the Sun's elevation at the mean distance: direct, diffuse, photons."""
    e = np.array(elev, float)
    d, f, q = light.sun(e)
    return dict(elevation_deg=e, direct_lux=d, diffuse_lux=f, total_lux=d + f, diffuse_share=f / np.maximum(d + f, 1e-30),
                umol=q, direct_normal_lux=np.where(e > 0, d / np.maximum(np.sin(np.radians(e)), 1e-9), 0.0))


def constants_record():
    from shared.provenance import constants_used
    return constants_used([Path(__file__)])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--extract-run', action='store_true', help='read the run\'s sunlight (needs netCDF4)')
    ap.add_argument('--figures', action='store_true')
    ap.add_argument('--workers', type=int, default=4)
    ap.add_argument('--output', default=str(OUT))
    args = ap.parse_args(argv)
    if args.extract_run:
        print(extract_run_sunlight())
        return 0
    product = build(workers=args.workers, with_figures=args.figures,
                    log=lambda m: print(time.strftime('%H:%M:%S'), m, flush=True))
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    print(f'wrote {out}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
