"""The GCM's forcing at one site, for a CM1 box no ring passes through: the GCM's own air there, its mean vertical wind,
and the heat and vapour its circulation brings by local time, from the GCM's own budget.

    climate/gcm/.venv/bin/python -m climate.crm.gcm_site A28R_301k_cold_moist 0 5 -5 -45 sea      # about a minute

CM1's boxes have taken their forcing from the rings they extend: the ring's upper air, the GCM's mean vertical wind
along it, and the heating and moistening of the rings' own day-night circulation. A box where no ring runs takes all
three from the GCM's column instead, over its cells within half_lat_deg and half_lon_deg of the site whose surface is
the site's (sea at sea level, or land).

The heat and vapour come from the GCM's own budget at each of its levels, the storage less the physics. The Moon's
column at 301 K holds about 220 kg/m2 of vapour, so its circulation moves water through a column at about a hundred
times the rate the column keeps of it, and the net at a site, a millimetre a day, is a residual that differences of
the winds and humidity across the GCM's cells cannot resolve: at the Near-side Sea centred and spectral differences of
the 3-day means bring the column +0.1 to +2.5 mm/day where the GCM's own budget takes away 0.7. PlaSim writes its
physics' heating and moistening of each layer as it applies them (with --heating-diagnostics and
--moisture-diagnostics): the heating at each grid point (code 21), here projected onto the model's T21 harmonics as its
spectral step takes it, and the moistening as that step applies it (code 40). The storage less these is what the
GCM's circulation brings, sideways and up or down, horizontal diffusion included. These diagnostics are snapshots of
a single time step at the end of each 3-day mean, so their mean over the ten times of day they fall at carries the
weather's noise and the rain's bursts; the vapour over the column is therefore held to the GCM's exact means, its
precipitable water's change less the evaporation less the rain, the difference spread over the levels in proportion
to the vapour each holds.

CM1 applies the GCM's mean vertical wind to its own air through the column. PlaSim's postprocessor writes omega (wap)
about a hundred times too weak at the equator against the meridional mass flux of the same winds, so the wind here is
the GCM's continuity on its own spectral fields: the divergence and the surface pressure's gradient from the written
winds and pressure projected onto the harmonics, integrated down each column. The table CM1 takes as prescribed (its
var11 forcing) holds the circulation's heating and moistening with the mean wind's own advection of the GCM's air
added back, so that in the GCM's state the two give the GCM's budget.

The means (ta, hus, ua, va, ps, czen, pr, evap, prw) cover a tenth of the lunar day each, 36 degrees of hour angle,
centred on the site's hour angle from the subsolar longitude of the mean solar height; the snapshots fall at each mean's
end, 18 degrees later. The year files' hus records the humidity at each interval's start, so its means come from
gcm/year_files.py, each the average of a record and the next, with each year file's first and last records left out
(from 27 October; the boxes forced before then took each record as its interval's mean). Each level's series is fitted
by least squares with the day's mean and first four harmonics, those of the means divided by the damping a 36-degree
mean gives them.

The table stands on CM1's own levels, the mean wind's advection added back on the box's starting air, so that the two
cancel at the start to the accuracy of the two differences (centred here, CM1's third-order upstream). The GCM holds
its tendencies as layer values between half levels, so each layer's is laid on CM1's levels by pressure, linear
within the layer about its value, its slope the monotonized-central limit of its neighbours' (flat in the top and
bottom layers), which keeps every layer's total; the part of the top layer above CM1's top, about 2% of the column's
mass, is left out. The first boxes (2026-10-21) took the table and the mean wind to 16 km only, as the rings'
day-night forcing did, interpolated the layers' values in height with the lowest held to the ground, and nudged their
air to the GCM's from 8 km, fully from 16 km. Over the Near-side Sea that left out the 35.5 W/m2 the circulation takes
from the 76% of the column's mass above 16 km, where the GCM's convection forms most of its rain; the nudging put back
83.5 W/m2 of heat there in its place, the box's convection stayed below 10 km and it rained 0.06 mm/day against the
GCM's 2.7; the interpolation made the column 0.55 mm/day drier than the GCM's. The forcing now reaches the top, and the
box holds its air to the GCM's only above the top of the GCM's convection (convection_top, the lowest level above
which its convection forms less than 1% of its convective rain): from the GCM's first level above that top, fully
from the next. The GCM's temperature falls to its top level at both sites (to 166 K at 129 km), so it resolves no
tropopause to hold from. A box cannot make the planet's day-night circulation, whose wind turns with the Sun (over the
Near-side Sea the wind 0.9 km up swings between about 3.7 m/s toward the east around noon and 3.8 m/s toward the west
in the evening), so its mean wind is nudged toward the GCM's by local time, a window of the day at a time.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RD, CP, LV = 287.04, 1005.7, 2.501e6
KAPPA = RD / CP
WATER_DENSITY = 1000.0                     # PlaSim's, in its diagnostics' water fluxes
TRUNCATION = 21                            # the GCM's triangular truncation
HARMONICS = 4
BINS = 36                                  # the table CM1 reads: bin centres from -175 degrees, 10 degrees apart
TOP_M = None                               # the table's top: the whole column (the first boxes stopped at 16 km)
CONVECTION_TOP_SHARE = 0.01                # the levels above the GCM's convection form less than this of its rain
MEANS = ('ta', 'hus', 'ua', 'va', 'ps', 'ts', 'tas', 'czen', 'pr', 'evap', 'prw', 'mrso')
SNAPSHOTS = ('dt_diab', 'dq_dyn')          # code 40, written under that name: the physics' moistening as applied
SCHEMES = ('dt_conv', 'dt_lsp')            # the convective and large-scale heating: where the GCM forms its rain


class Spectral:
    """Triangular spherical harmonics on a Gaussian grid (rows north to south): analysis of grid fields, and synthesis
    of a field, its derivative in mu = sin(latitude) or in longitude."""

    def __init__(self, nlat: int, nlon: int, truncation: int = TRUNCATION):
        from scipy.special import assoc_legendre_p_all
        nodes, weights = np.polynomial.legendre.leggauss(nlat)
        order = np.argsort(-nodes)
        self.mu, self.w = nodes[order], weights[order]
        p = np.asarray(assoc_legendre_p_all(truncation, truncation, self.mu, norm=True, diff_n=1))
        self.p, self.dp = p[0, :, :truncation + 1, :], p[1, :, :truncation + 1, :]       # (n, m, lat)
        self.t, self.nlon = truncation, nlon
        self.m = np.arange(truncation + 1)
        self.cos = np.sqrt(1.0 - self.mu**2)

    def analyse(self, f):
        """Coefficients (..., n, m) of a grid field (..., lat, lon)."""
        fourier = np.fft.rfft(np.asarray(f, float), axis=-1)[..., :self.t + 1] / self.nlon
        return np.einsum('...jm,nmj,j->...nm', fourier, self.p, self.w)

    def synthesise(self, s, along=None):
        """The field (..., lat, lon), or its derivative along 'mu' or 'lambda' (radians), from coefficients."""
        f = np.einsum('...nm,nmj->...jm', s, self.dp if along == 'mu' else self.p)
        if along == 'lambda':
            f = f * 1j * self.m
        full = np.zeros(f.shape[:-1] + (self.nlon // 2 + 1,), complex)
        full[..., :self.t + 1] = f * self.nlon
        return np.fft.irfft(full, n=self.nlon, axis=-1)

    def project(self, f):
        return self.synthesise(self.analyse(f))


def omega(u, v, ps, sigma, radius, sp: Spectral, wind: Spectral):
    """PlaSim's continuity on its spectral fields: omega (Pa/s) at the full levels (time, lev, lat, lon) from the
    winds and surface pressure (Pa), omega = sigma ps v.grad(ln ps) less ps times the sum from the top of
    (div v + v.grad(ln ps)) over the layers above and the part of its own layer above the level, each layer's value
    its full level's, as PlaSim's vertical sums take them; the divergence from (u cos, v cos) on harmonics one degree
    higher (wind), as a T21 flow's components carry them, and the gradient of ln ps from its T21 projection (sp)."""
    cos = wind.cos[:, None]
    div = (wind.synthesise(wind.analyse(u * cos), 'lambda') / (radius * cos**2)
           + wind.synthesise(wind.analyse(v * cos), 'mu') / radius)
    lnps = sp.analyse(np.log(ps))
    gx = sp.synthesise(lnps, 'lambda') / (radius * cos)
    gy = sp.synthesise(lnps, 'mu') * cos / radius
    vgrad = u * gx[:, None] + v * gy[:, None]
    edges = plasim_edges(sigma)
    thick = np.diff(edges)[None, :, None, None]
    total = div + vgrad
    above = np.cumsum(total * thick, axis=1) - total * thick + total * (sigma - edges[:-1])[None, :, None, None]
    return ps[:, None] * (sigma[None, :, None, None] * vgrad - above)


def plasim_edges(sigma):
    from climate.crm.cm1_run import plasim_half_levels
    return np.asarray(plasim_half_levels(sigma), float)


def site_cells(lat, lon, lsm, height, site):
    """The GCM cells of the site: within half_lat_deg and half_lon_deg of it, sea at sea level or land."""
    rows = np.abs(lat - site['lat_deg']) <= site['half_lat_deg']
    cols = np.abs((lon - site['lon_deg'] + 180.0) % 360.0 - 180.0) <= site['half_lon_deg']
    cells = rows[:, None] & cols[None, :]
    if site['surface'] == 'sea':
        cells &= ~lsm & (height < 1.0)
    elif site['surface'] == 'land':
        cells &= lsm
    else:
        raise ValueError(f"surface is sea or land, not {site['surface']!r}")
    if not cells.any():
        raise ValueError(f'no {site["surface"]} cells at the site {site}')
    return cells


def subsolar_longitude(czen, lat, lon):
    """The longitude (degrees) each output's mean solar height centres on, near the equator: the circular mean of
    longitude weighted by the mean cosine of the solar zenith, which a mean over a window of the day keeps symmetric
    about the window's centre."""
    rows = np.abs(lat) < 20.0
    weight = czen[:, rows, :].mean(axis=1)                                  # (time, lon)
    a = np.radians(lon)
    return np.degrees(np.arctan2((weight * np.sin(a)).sum(axis=1), (weight * np.cos(a)).sum(axis=1))) % 360.0


def level_heights(ta, hus, sigma, gravity):
    """Heights (m) of the full levels above the ground (..., lev), top first, from the virtual temperatures."""
    tv = ta * (1.0 + 0.608 * hus)
    z = np.zeros_like(tv)
    z[..., -1] = RD * tv[..., -1] / gravity * np.log(1.0 / sigma[-1])
    for k in range(sigma.size - 2, -1, -1):
        z[..., k] = z[..., k + 1] + RD * 0.5 * (tv[..., k] + tv[..., k + 1]) / gravity * np.log(sigma[k + 1] / sigma[k])
    return z


def basis(hour_deg, harmonics=HARMONICS):
    a = np.radians(np.asarray(hour_deg, float))
    return np.column_stack([np.ones_like(a)] + [f(n * a) for n in range(1, harmonics + 1) for f in (np.cos, np.sin)])


def harmonic_coefficients(hour_deg, values, harmonics=HARMONICS):
    """Least-squares coefficients of a series values (time, ...) against the hour angle: the mean, then the cosine and
    sine of each harmonic."""
    b = basis(hour_deg, harmonics)
    values = np.asarray(values, float)
    coef, *_ = np.linalg.lstsq(b, values.reshape(len(b), -1), rcond=None)
    return coef.reshape((b.shape[1],) + values.shape[1:])


def undamped(coef, window_deg, harmonics=HARMONICS):
    """The coefficients with each harmonic divided by the damping a mean over window_deg of hour angle gives it."""
    out = np.array(coef, float)
    for n in range(1, harmonics + 1):
        x = math.radians(n * window_deg / 2.0)
        out[2 * n - 1:2 * n + 1] /= math.sin(x) / x
    return out


def evaluate(coef, hour_deg, harmonics=HARMONICS):
    return np.tensordot(basis(hour_deg, harmonics), coef, axes=(1, 0))


def rate(coef, hour_deg, day_s, harmonics=HARMONICS):
    """The time derivative (per second) of a fitted series at the hour angles, which grow 360 degrees a solar day."""
    a = np.radians(np.asarray(hour_deg, float))
    d = np.column_stack([np.zeros_like(a)] + [g for n in range(1, harmonics + 1)
                                              for g in (-n * np.sin(n * a), n * np.cos(n * a))])
    return np.tensordot(d, coef, axes=(1, 0)) * 2.0 * math.pi / day_s


def remap(values, edges_pa, targets_pa):
    """A field the GCM holds as layer values (..., layer) between the pressures edges_pa (top first, increasing) as
    means over the target layers between targets_pa (in either order): linear in pressure within each source layer
    about its value, its slope the monotonized-central limit of the neighbouring layers' (flat in the top and bottom
    layers), integrated over the part of each target layer the source layer covers. Each source layer's integral in
    pressure is kept wherever the targets cover it."""
    v = np.asarray(values, float)
    e = np.asarray(edges_pa, float)
    mid = 0.5 * (e[1:] + e[:-1])
    slope = np.zeros_like(v)
    if mid.size > 2:
        up = (v[..., 1:-1] - v[..., :-2]) / (mid[1:-1] - mid[:-2])
        down = (v[..., 2:] - v[..., 1:-1]) / (mid[2:] - mid[1:-1])
        centre = (v[..., 2:] - v[..., :-2]) / (mid[2:] - mid[:-2])
        limit = 2.0 * np.minimum(np.abs(up), np.abs(down))
        slope[..., 1:-1] = np.where(up * down > 0.0, np.sign(centre) * np.minimum(np.abs(centre), limit), 0.0)
    t = np.asarray(targets_pa, float)
    lo, hi = np.minimum(t[1:], t[:-1]), np.maximum(t[1:], t[:-1])
    out = np.zeros(v.shape[:-1] + lo.shape)
    for k in range(mid.size):
        a, b = np.clip(lo, e[k], e[k + 1]), np.clip(hi, e[k], e[k + 1])
        out += v[..., k, None] * (b - a) + slope[..., k, None] * 0.5 * ((b - mid[k])**2 - (a - mid[k])**2)
    return out / (hi - lo)


def column_pressures(zw, theta_at, qv_at, surface_pa, gravity, steps=20):
    """Pressure (Pa) at the heights zw (the ground up) in hydrostatic balance with air of potential temperature
    theta_at(z) and vapour qv_at(z), from the surface pressure: the Exner function falls by g / (cp theta_v) per metre,
    integrated by the trapezoid rule on steps sub-layers of each layer."""
    zw = np.asarray(zw, float)
    z = np.concatenate([np.linspace(a, b, steps, endpoint=False) for a, b in zip(zw[:-1], zw[1:])] + [zw[-1:]])
    inv = 1.0 / (np.asarray(theta_at(z), float) * (1.0 + 0.608 * np.asarray(qv_at(z), float)))
    exner = (surface_pa / 1e5)**KAPPA - gravity / CP * np.r_[0.0, np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(z))]
    return 1e5 * exner[::steps]**(1.0 / KAPPA)


def convection_top(formed, z, share=CONVECTION_TOP_SHARE):
    """The top of the GCM's convection from the rain it forms in each layer (top first) at heights z (m): the lowest
    level above which less than share of it forms; and the heights from which a box holds its air to the GCM's, the
    next level up and, fully, the one above that (None past the GCM's top)."""
    formed = np.maximum(np.asarray(formed, float), 0.0)
    if formed.sum() <= 0.0:
        raise ValueError('the GCM forms no convective rain here')
    above = np.r_[0.0, np.cumsum(formed)[:-1]] / formed.sum()            # formed above each level
    k = int(np.flatnonzero(above < share).max())
    level = lambda i: float(z[i]) if i >= 0 else None
    return dict(level=k, top_m=float(z[k]), share=share, held_from_m=level(k - 1), held_fully_from_m=level(k - 2))


def site_forcing(run, years, site, window_deg=36.0, top_m=TOP_M, series=False):
    """The GCM's column at a site over the years: its mean air and surface, its mean vertical wind, the table CM1
    prescribes by local time, and the budgets behind it. Arrays by level run top first, as PlaSim writes them. With
    series, also the cells' mean of every record behind the fits, its index in the run's records, and the fits'
    coefficients (site_sequence reads them)."""
    import netCDF4
    from climate.crm.cm1_run import gcm_soil, planet, solar_day_s
    from climate.gcm import year_files as yf
    p = planet()
    gravity, radius = p['gravity_m_s2'], p['radius_m']
    folder = HERE.parent / 'gcm' / 'runs' / run / 'model'
    keep = {k: [] for k in ('hour', 'theta', 'q', 'u', 'v', 'z', 'rho', 'omega', 'heat', 'moist', 'ps', 'ts', 'tas',
                            'rain', 'evaporation', 'water', 'wetness', 'conv', 'lsp', 'index')}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            missing = [k for k in MEANS + SNAPSHOTS + SCHEMES if k not in d.variables]
            if missing:
                raise ValueError(f'{run} year {year} lacks {missing}: run it with --heating-diagnostics '
                                 '--moisture-diagnostics')
            f = yf.fields(d, MEANS + SNAPSHOTS + SCHEMES)            # the humidity's means: gcm/year_files.py
            count = len(d['time'])
            keep['index'].append(year * count + np.arange(count)[yf.records(d)])   # each record's place in the run
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            lsm = np.asarray(d['lsm'][:], float)[0] > 0.5
        if year == years[0]:
            height = (np.loadtxt(folder.parent / 'inputs' / 'moon_topography.sra', skiprows=1).ravel()
                      .reshape(lsm.shape) / gravity)
            cells = site_cells(lat, lon, lsm, height, site)
            jj, ii = np.nonzero(cells)
            weight = np.cos(np.radians(lat[jj]))
            weight = weight / weight.sum()
            thick = np.diff(plasim_edges(sigma))
            soil = gcm_soil(folder)
            sp, wind = Spectral(lat.size, lon.size), Spectral(lat.size, lon.size, TRUNCATION + 1)
        at = lambda a: a[..., jj, ii]                                        # (..., cells)
        mean = lambda a: np.tensordot(a, weight, axes=(-1, 0))               # the cells' mean
        grid_ps = f['ps'] * 100.0                                            # hPa in the output
        ps = at(grid_ps)                                                     # (time, cells)
        ta, hus = at(f['ta']), at(f['hus'])                                  # (time, lev, cells)
        exner = (sigma[None, :, None] * ps[:, None, :] / 1e5)**KAPPA
        rho = sigma[None, :, None] * ps[:, None, :] / (RD * ta * (1.0 + 0.608 * hus))
        water = WATER_DENSITY * gravity / (ps[:, None, :] * thick[None, :, None])     # water flux to kg/kg/s
        keep['theta'].append(mean(ta / exner))
        keep['q'].append(mean(hus))
        keep['u'].append(mean(at(f['ua'])))
        keep['v'].append(mean(at(f['va'])))
        keep['rho'].append(mean(rho))
        keep['z'].append(mean(np.moveaxis(level_heights(np.moveaxis(ta, 1, -1), np.moveaxis(hus, 1, -1), sigma,
                                                        gravity), -1, 1)))
        keep['omega'].append(mean(at(omega(f['ua'], f['va'], grid_ps, sigma, radius, sp, wind))))
        keep['heat'].append(mean(at(sp.project(f['dt_diab'])) / exner))     # the physics' heating, as theta, at T21
        keep['moist'].append(mean(at(f['dq_dyn']) * water))                 # the physics' moistening as applied
        formed = CP * ps[:, None, :] * thick[None, :, None] / (gravity * LV)   # K/s of latent heating to kg/m2/s
        keep['conv'].append(mean(at(f['dt_conv']) * formed))
        keep['lsp'].append(mean(at(f['dt_lsp']) * formed))
        keep['ps'].append(mean(ps))
        keep['ts'].append(mean(at(f['ts'])))
        keep['tas'].append(mean(at(f['tas'])))
        keep['rain'].append(mean(at(f['pr'])) * WATER_DENSITY)               # kg/m2/s
        keep['evaporation'].append(-mean(at(f['evap'])) * WATER_DENSITY)    # upward positive
        keep['water'].append(mean(at(f['prw'])))                             # kg/m2
        keep['wetness'].append(mean(at(np.minimum(1.0, f['mrso'] / (soil['drhsfull'] * soil['wsmax'])))))
        keep['hour'].append((site['lon_deg'] - subsolar_longitude(f['czen'], lat, lon) + 180.0) % 360.0 - 180.0)
    s = {k: np.concatenate(v) for k, v in keep.items()}
    tm = lambda k: s[k].mean(axis=0)
    day_s = solar_day_s()
    hours = (np.arange(BINS) + 0.5) * 360.0 / BINS - 180.0
    centre, end = s['hour'], (s['hour'] + window_deg / 2.0 + 180.0) % 360.0 - 180.0
    fit_mean = lambda k: undamped(harmonic_coefficients(centre, s[k]), window_deg)
    fit_snap = lambda k: harmonic_coefficients(end, s[k])
    theta_c, q_c, water_c, rain_c, evap_c = (fit_mean(k) for k in ('theta', 'q', 'water', 'rain', 'evaporation'))
    u_c, v_c = fit_mean('u'), fit_mean('v')
    heat_c, moist_c = fit_snap('heat'), fit_snap('moist')
    z, rho = tm('z'), tm('rho')
    w_gcm = -tm('omega') / (rho * gravity)
    below = z <= (np.inf if top_m is None else top_m)
    w_mean = np.where(below, w_gcm, 0.0)                                     # none above a top, where one is set
    mass = tm('ps') * thick / gravity                                        # kg/m2 in each layer
    dz = lambda a: np.gradient(a, z, axis=-1)
    theta_h, q_h = evaluate(theta_c, hours), evaluate(q_c, hours)            # (bins, lev)
    heating = rate(theta_c, hours, day_s) - evaluate(heat_c, hours)          # the circulation's, by local time
    moistening = rate(q_c, hours, day_s) - evaluate(moist_c, hours)
    # The column's vapour held to the exact means: the precipitable water's change less the evaporation, plus the rain.
    exact = rate(water_c, hours, day_s) - evaluate(evap_c, hours) + evaluate(rain_c, hours)          # kg/m2/s
    snapped = moistening @ mass
    share = q_h * mass[None, :] / (q_h @ mass)[:, None]                      # each level's share of the vapour
    correction = (exact - snapped)[:, None] * share / mass[None, :]
    moistening = moistening + correction
    table_th = heating + w_mean[None, :] * dz(theta_h)                       # the mean wind's own part added back
    table_q = moistening + w_mean[None, :] * dz(q_h)
    table_th, table_q = np.where(below, table_th, 0.0), np.where(below, table_q, 0.0)
    circ_th, circ_q = np.where(below, heating, 0.0), np.where(below, moistening, 0.0)
    per_day = 86400.0
    extra = dict(s, coefficients=dict(theta=theta_c, q=q_c, u=u_c, v=v_c), mass=mass) if series else None
    phase = np.round(((centre + 180.0) % 360.0) / window_deg - 0.5).astype(int) % int(round(360.0 / window_deg))
    out = dict(
        run=run, years=list(years), site=dict(site), window_deg=window_deg, harmonics=HARMONICS, truncation=TRUNCATION,
        cells=[[float(lat[j]), float(lon[i])] for j, i in zip(jj, ii)],
        cell_heights_m=[float(height[j, i]) for j, i in zip(jj, ii)],
        outputs=int(centre.size), phases=int(np.unique(phase).size), sigma=sigma.tolist(),
        profile=dict(z_m=z.tolist(), theta_k=tm('theta').tolist(), qv_kg_kg=tm('q').tolist(), u_m_s=tm('u').tolist(),
                     v_m_s=tm('v').tolist(), density_kg_m3=rho.tolist()),
        surface=dict(pressure_pa=float(tm('ps')), ground_k=float(tm('ts')), air_2m_k=float(tm('tas')),
                     wetness=float(tm('wetness')),
                     height_m=float(np.mean([height[j, i] for j, i in zip(jj, ii)]))),
        vertical_wind_m_s=w_mean.tolist(), gcm_vertical_wind_m_s=w_gcm.tolist(),
        day_night=dict(hour_angle_deg=hours.tolist(), z_m=z.tolist(), theta_k_s=table_th.tolist(),
                       qv_kg_kg_s=table_q.tolist(), top_m=top_m),
        circulation=dict(hour_angle_deg=hours.tolist(), z_m=z.tolist(), theta_k_s=circ_th.tolist(),
                         qv_kg_kg_s=circ_q.tolist(), top_m=top_m),
        convection=dict(convective_rain_formed_mm_day=(tm('conv') * per_day).tolist(),
                        large_scale_rain_formed_mm_day=(tm('lsp') * per_day).tolist(),
                        **convection_top(tm('conv'), z)),
        wind_by_hour=dict(hour_angle_deg=hours.tolist(), z_m=z.tolist(), u_m_s=evaluate(u_c, hours).tolist(),
                          v_m_s=evaluate(v_c, hours).tolist()),
        budget=dict(
            circulation_heating_theta_k_day=(heating.mean(axis=0) * per_day).tolist(),
            circulation_moistening_g_kg_day=(moistening.mean(axis=0) * per_day * 1e3).tolist(),
            physics_heating_theta_k_day=(heat_c[0] * per_day).tolist(),            # the fits' means over the day
            physics_moistening_g_kg_day=(moist_c[0] * per_day * 1e3).tolist(),
            mean_wind_advection_theta_k_day=(-w_mean * dz(tm('theta')) * per_day).tolist(),
            mean_wind_advection_qv_g_kg_day=(-w_mean * dz(tm('q')) * per_day * 1e3).tolist(),
            table_mean_theta_k_day=(table_th.mean(axis=0) * per_day).tolist(),
            table_mean_qv_g_kg_day=(table_q.mean(axis=0) * per_day * 1e3).tolist(),
            column_vapour_mm_day=dict(exact=float(exact.mean() * per_day), from_snapshots=float(snapped.mean() * per_day),
                                      exact_by_hour=[float(x) for x in exact * per_day],
                                      from_snapshots_by_hour=[float(x) for x in snapped * per_day],
                                      held_by_hour=[float(x) for x in (moistening @ mass) * per_day]),
            column_heating_w_m2=dict(physics=float(CP * (heat_c[0] * mass * (tm('ps') * sigma / 1e5)**KAPPA).sum()),
                                     rain_latent=float(LV * tm('rain'))),
            rain_mm_day=float(tm('rain') * per_day), evaporation_mm_day=float(tm('evaporation') * per_day)))
    return (out, extra) if series else out


def to_cm1(forcing, zw, theta_at=None, qv_at=None, gravity=None):
    """The forcing as CM1's box reads it: the mean vertical wind on its w levels zw (zero at the ground, falling to
    zero at the top past the GCM's top level), and the table of prescribed heating (potential temperature, K/s) and
    moistening (kg/kg/s) by local time. Given the box's starting air (theta_at and qv_at, functions of height), the
    table stands on CM1's own scalar levels: each GCM layer's heating and moistening by the circulation laid on them
    by pressure (remap), the levels' pressures those of the starting air in balance from the GCM's surface pressure,
    with the mean wind's advection of that starting air added back, so that at the start CM1's own advection by the
    mean wind takes it away again; nil above the forcing's top, where it has one. Without them it stands on the GCM's
    levels with the mean wind's advection of the GCM's air added back."""
    z = np.asarray(forcing['profile']['z_m'], float)[::-1]                # the ground up
    w = np.asarray(forcing['vertical_wind_m_s'], float)[::-1]
    zw = np.asarray(zw, float)
    zz, ww = np.r_[0.0, z, zw[-1]], np.r_[0.0, w, 0.0]
    w_zw = np.interp(zw, zz, ww)
    hours = np.asarray(forcing['day_night']['hour_angle_deg'], float)
    if theta_at is None:
        dn = forcing['day_night']
        th = np.asarray(dn['theta_k_s'], float)[:, ::-1]
        qv = np.asarray(dn['qv_kg_kg_s'], float)[:, ::-1]
        return dict(w_m_s=[float(x) for x in w_zw],
                    day_night=dict(z_m=np.r_[0.0, z], theta_k_s=np.c_[th[:, :1], th], qv_kg_kg_s=np.c_[qv[:, :1], qv],
                                   hour_angle_deg=hours))
    if gravity is None:
        from climate.crm.cm1_run import planet
        gravity = planet()['gravity_m_s2']
    zh = 0.5 * (zw[1:] + zw[:-1])
    w_zh = 0.5 * (w_zw[1:] + w_zw[:-1])
    circ = forcing['circulation']
    ps = forcing['surface']['pressure_pa']
    edges = plasim_edges(forcing['sigma']) * ps                            # the GCM's half levels, top first
    p_w = column_pressures(zw, theta_at, qv_at, ps, gravity)
    on_box = lambda a: remap(a, edges, p_w)                                # (bins, CM1's levels from the ground)
    th0, q0 = np.asarray(theta_at(zh), float), np.asarray(qv_at(zh), float)
    top = circ.get('top_m')
    below = zh <= (np.inf if top is None else top)
    th = np.where(below, on_box(circ['theta_k_s']) + (w_zh * np.gradient(th0, zh))[None, :], 0.0)
    qv = np.where(below, on_box(circ['qv_kg_kg_s']) + (w_zh * np.gradient(q0, zh))[None, :], 0.0)
    return dict(w_m_s=[float(x) for x in w_zw], pressure_pa=p_w,
                day_night=dict(z_m=zh, theta_k_s=th, qv_kg_kg_s=qv, hour_angle_deg=hours))


def wind_windows(forcing, zh, start_hour_deg, day_s, length_s, bins=BINS):
    """The GCM's wind by local time as a series of windows for CM1's nudging: each window a bin of the day, from the
    box's start through its length, with the wind (m/s) on the box's scalar levels zh at the bin's middle hour angle,
    held at the GCM's lowest level down to the ground and at its top level above."""
    wind = forcing['wind_by_hour']
    z = np.asarray(wind['z_m'], float)[::-1]
    hours = np.asarray(wind['hour_angle_deg'], float)
    u, v = (np.asarray(wind[k], float)[:, ::-1] for k in ('u_m_s', 'v_m_s'))
    step = day_s / bins
    out, t = [], 0.0
    while t < length_s:
        middle = (start_hour_deg + 360.0 * (t + 0.5 * step) / day_s + 180.0) % 360.0 - 180.0
        b = int(np.argmin(np.abs((hours - middle + 180.0) % 360.0 - 180.0)))
        out.append(dict(time_s=(t, t + step), hour_angle_deg=float(middle),
                        u_m_s=np.interp(zh, z, u[b]), v_m_s=np.interp(zh, z, v[b])))
        t += step
    return out



# The GCM's own weather at a site (2026-11-01). A box forced by the fits above follows the mean circulation by local
# time and carries none of the weather the GCM's evaporation rests on (gap_boxes/findings_sea.md). site_sequence keeps
# the fits as the forcing's mean part and adds each record's departure from them, so that a box can follow the
# design years' own sequence: the states (air, wind) by their departures at the records' centres, the circulation's
# heating and moistening by their departures at the boundaries between consecutive records (where the physics'
# snapshots fall), each the storage between the two records less the physics' snapshot, less the same operator on the
# fits, with the column's vapour departure held to the exact means' budget across the boundary. The records are the
# year files' 3-day means (the weather between the means) or, from a rerun writing its means more often, means of a
# few hours (the weather within them); window_deg is then the rerun's interval in hour angle.
def record_seconds(run, year=False) -> float:
    """The length (s) of one output interval of a GCM run, or with year of one model year, from its progress record."""
    import json
    m = json.loads((HERE.parent / 'gcm' / 'runs' / run / 'progress.json').read_text())['configuration']['model']
    step = m['timestep_min'] * 60.0
    return step * m['steps_per_lunar_day'] * (m['lunar_days_per_year'] if year else 1.0 / m['writes_per_lunar_day'])


def site_sequence(run, years, site, window_deg=36.0, top_m=TOP_M, fits=None):
    """The GCM's sequence at a site: the mean forcing (site_forcing, from fits = (run, years) when given, else from
    these records) and the departures of these records from its fits. Times are seconds from the run's start; arrays
    by level run top first. Departures above the forcing's top, where it has one, are nil."""
    f, s = site_forcing(run, years, site, window_deg, top_m, series=True)
    if fits is not None:                                                   # the mean part from other years' records
        f, s_fit = site_forcing(*fits, site, 36.0, top_m, series=True)       # the year files' 3-day means
    else:
        s_fit = s
    step = record_seconds(run)
    centre = s['hour']
    end = (centre + window_deg / 2.0 + 180.0) % 360.0 - 180.0
    fc, fe = s_fit['hour'], (s_fit['hour'] + f['window_deg'] / 2.0 + 180.0) % 360.0 - 180.0
    raw = {k: harmonic_coefficients(fc, s_fit[k]) for k in ('theta', 'q', 'u', 'v', 'water', 'rain', 'evaporation')}
    phys = dict(theta=harmonic_coefficients(fe, s_fit['heat']), q=harmonic_coefficients(fe, s_fit['moist']))
    # The fits' window means at a record's phase: the means' own fit, damped from window_deg of the fits' records to
    # this sequence's window; for the year files' own records the fit itself.
    damp = lambda coef: coef if fits is None else redamped(coef, f['window_deg'], window_deg)
    mean_at = lambda k, h: evaluate(damp(raw[k]), h)
    states = {k: s[k] - mean_at(k, centre) for k in ('theta', 'q', 'u', 'v')}
    i = np.nonzero(np.diff(s['index']) == 1)[0]                            # consecutive records
    depart = {}
    for k, physics in (('theta', 'heat'), ('q', 'moist')):
        own = (s[k][i + 1] - s[k][i]) / step - s[physics][i]
        fit = (mean_at(k, centre[i + 1]) - mean_at(k, centre[i])) / step - evaluate(phys[k], end[i])
        depart[k] = own - fit
    budget = lambda w, e, r: (w[1] - w[0]) / step - 0.5 * (e[0] + e[1]) + 0.5 * (r[0] + r[1])
    exact = budget(*([s[k][i], s[k][i + 1]] for k in ('water', 'evaporation', 'rain')))
    exact_fit = budget(*([mean_at(k, centre[i]), mean_at(k, centre[i + 1])] for k in ('water', 'evaporation', 'rain')))
    mass = s_fit['mass']
    q_b = 0.5 * (s['q'][i] + s['q'][i + 1])
    share = q_b * mass[None, :] / (q_b @ mass)[:, None]
    depart['q'] = depart['q'] + ((exact - exact_fit) - depart['q'] @ mass)[:, None] * share / mass[None, :]
    z = np.asarray(f['profile']['z_m'], float)
    below = z <= (np.inf if top_m is None else top_m)
    depart = {k: np.where(below, v, 0.0) for k, v in depart.items()}
    t = (s['index'] + 0.5) * step
    order = np.argsort(t)
    unwrapped = np.degrees(np.unwrap(np.radians(centre[order])))
    slope, h0 = np.polyfit(t[order], unwrapped, 1)                          # the site's hour angle against time
    return dict(forcing=f, run=run, years=list(years), record_s=step, year_s=record_seconds(run, year=True),
                window_deg=window_deg,
                sun=dict(period_s=float(360.0 / slope), hour_at_zero_deg=float(h0),
                         residual_deg=float(np.abs(unwrapped - (h0 + slope * t[order])).max())),
                coefficients=s_fit['coefficients'], mean_profile=dict(theta_k=s_fit['theta'].mean(axis=0),
                                                                      qv_kg_kg=s_fit['q'].mean(axis=0)),
                states=dict(time_s=t, hour_angle_deg=centre, **{f'{k}_departure': v for k, v in states.items()}),
                boundaries=dict(time_s=(s['index'][i] + 1.0) * step, theta_k_s=depart['theta'], qv_kg_kg_s=depart['q'],
                                column_vapour_kg_m2_s=exact, column_vapour_fit_kg_m2_s=exact_fit),
                sea_k=dict(time_s=t, value=s['ts']), rain_kg_m2_s=s['rain'], evaporation_kg_m2_s=s['evaporation'])


def redamped(coef, from_deg, to_deg, harmonics=HARMONICS):
    """Coefficients of a fit to means over from_deg of hour angle, as means over to_deg would carry them."""
    return undamped(coef, from_deg, harmonics) / undamped(np.ones_like(coef), to_deg, harmonics)


def gcm_time(seq, hour_deg, after_s):
    """The first time (s from the GCM run's start) at or after after_s when the site's hour angle is hour_deg."""
    period, h0 = seq['sun']['period_s'], seq['sun']['hour_at_zero_deg']
    base = (hour_deg - h0) / 360.0 * period
    return base + period * math.ceil((after_s - base) / period - 1e-9)


def sequence_state(seq, key, t_s, hour_deg):
    """The GCM's air or wind (theta, q, u or v) at its levels at GCM times t_s (s) and the site's hour angles there:
    the fit by local time (undamped) plus the records' departures interpolated linearly in time between their
    centres; (times, levels), top first."""
    t_s, hour_deg = np.atleast_1d(np.asarray(t_s, float)), np.atleast_1d(np.asarray(hour_deg, float))
    st = seq['states']
    order = np.argsort(st['time_s'])
    tt, dep = st['time_s'][order], st[f'{key}_departure'][order]
    if t_s.min() < tt[0] or t_s.max() > tt[-1]:
        raise ValueError(f'times {t_s.min():.0f}-{t_s.max():.0f} s lie outside the records ({tt[0]:.0f}-{tt[-1]:.0f} s)')
    depart = np.stack([np.interp(t_s, tt, dep[:, k]) for k in range(dep.shape[1])], axis=-1)
    return evaluate(seq['coefficients'][key], hour_deg) + depart


def sequence_tendency(seq, key, t_s):
    """The departure of the GCM's circulation heating (theta, K/s) or moistening (q, kg/kg/s) from its fit by local
    time at GCM times t_s, linear in time between the boundaries; (times, levels), top first."""
    t_s = np.atleast_1d(np.asarray(t_s, float))
    b = seq['boundaries']
    order = np.argsort(b['time_s'])
    tt, dep = b['time_s'][order], b[f'{key}_k_s' if key == 'theta' else 'qv_kg_kg_s'][order]
    if t_s.min() < tt[0] or t_s.max() > tt[-1]:
        raise ValueError(f'times {t_s.min():.0f}-{t_s.max():.0f} s lie outside the boundaries ({tt[0]:.0f}-{tt[-1]:.0f} s)')
    return np.stack([np.interp(t_s, tt, dep[:, k]) for k in range(dep.shape[1])], axis=-1)


# The damped gravity wave's boxes (5 November 2026; session-tools/agents/gap_boxes/yardstick/rules_wtg.txt and the offline
# step's report). A box whose large-scale vertical wind answers its own heating (cm1_run.DGW_PATCHES) takes three things
# from the GCM beyond the forcing above. The air it is held near and draws in is the site's column by local time: in the
# design years the column's band stands within 0.17 K of the tropical mean at every tenth of the lunar day at the three
# sites, and converging air comes from the region around the site. Its table holds the circulation's horizontal part
# alone (the circulation's tendencies less the vertical advection by the GCM's own vertical wind at that hour, its omega
# from continuity fitted by local time like the air), since the vertical advection is the box's to make. The first
# baroclinic gravity wave's speed comes from the tropical mean, the large-scale air the wave crosses.
def dgw_forcing(run, years, site, window_deg=36.0):
    """The GCM's inputs to a damped-gravity-wave box at a site: site_forcing's result (forcing); the column's air by
    local time (air_by_hour: potential temperature, K, and vapour, kg/kg, at the BINS hour angles, the fits' values);
    its vertical wind by local time (m/s); the circulation's horizontal part by local time (horizontal: theta K/s and
    vapour kg/kg/s); and the column sums (mm/day of vapour by hour) of the circulation, its vertical advection by the
    GCM's own vertical wind and the remainder. Arrays by level run top first."""
    from climate.crm.cm1_run import planet
    out, extra = site_forcing(run, years, site, window_deg=window_deg, series=True)
    g = planet()['gravity_m_s2']
    hours = np.asarray(out['circulation']['hour_angle_deg'], float)
    coef = extra['coefficients']
    theta_h, q_h = evaluate(coef['theta'], hours), evaluate(coef['q'], hours)
    rho, z, mass = extra['rho'].mean(axis=0), extra['z'].mean(axis=0), extra['mass']
    w_h = -evaluate(undamped(harmonic_coefficients(extra['hour'], extra['omega']), window_deg), hours) / (rho[None, :] * g)
    dz = lambda a: np.gradient(a, z, axis=-1)
    circ_th = np.asarray(out['circulation']['theta_k_s'], float)
    circ_q = np.asarray(out['circulation']['qv_kg_kg_s'], float)
    horiz_th, horiz_q = circ_th + w_h * dz(theta_h), circ_q + w_h * dz(q_h)   # less the vertical advection -w d/dz
    column = lambda a: (a @ mass) * 86400.0
    return dict(run=run, years=list(years), site=dict(site), window_deg=window_deg, forcing=out,
                air_by_hour=dict(hour_angle_deg=hours.tolist(), z_m=z.tolist(), theta_k=theta_h.tolist(),
                                 qv_kg_kg=q_h.tolist()),
                vertical_wind_by_hour=dict(hour_angle_deg=hours.tolist(), z_m=z.tolist(), w_m_s=w_h.tolist()),
                horizontal=dict(hour_angle_deg=hours.tolist(), z_m=z.tolist(), theta_k_s=horiz_th.tolist(),
                                qv_kg_kg_s=horiz_q.tolist(), top_m=out['circulation']['top_m']),
                column_vapour_mm_day=dict(circulation=column(circ_q).tolist(), horizontal=column(horiz_q).tolist(),
                                          vertical_advection=column(circ_q - horiz_q).tolist()))


def first_mode(run, years, lat_limit_deg=30.0):
    """The first baroclinic gravity wave's speed over the tropics (m/s): the WKB estimate c = (1/pi) * integral of N dz
    on the lunar-day mean of the GCM's tropical-mean potential temperature (cells within lat_limit_deg, cos-latitude
    weights), from the ground to the top of the tropical convection (convection_top's rule on the tropics' convective
    rain by level), and to the model's top level beside it."""
    import netCDF4
    from climate.crm.cm1_run import planet
    from climate.gcm import year_files as yf
    g = planet()['gravity_m_s2']
    folder = HERE.parent / 'gcm' / 'runs' / run / 'model'
    sums, count = {}, 0
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            f = yf.fields(d, ('ta', 'hus', 'ps', 'dt_conv'))
            lat = np.asarray(d['lat'][:], float)
            sigma = np.asarray(d['lev'][:], float)
        ps = f['ps'] * 100.0
        w = np.where(np.abs(lat) <= lat_limit_deg, np.cos(np.radians(lat)), 0.0)[:, None] * np.ones(ps.shape[-1])[None, :]
        w /= w.sum()
        tm = lambda a: np.tensordot(a, w, axes=([-2, -1], [0, 1])).sum(axis=0)          # summed over records
        thick = np.diff(plasim_edges(sigma))
        for key, value in (('theta', f['ta'] / (sigma[None, :, None, None] * ps[:, None] / 1e5) ** KAPPA),
                           ('ta', f['ta']), ('q', f['hus']),
                           ('conv', f['dt_conv'] * CP * ps[:, None] * thick[None, :, None, None] / (g * LV))):
            sums[key] = sums.get(key, 0.0) + tm(value)
        count += ps.shape[0]
    mean = {k: v / count for k, v in sums.items()}
    z = level_heights(mean['ta'], mean['q'], sigma, g)                                 # top first
    top = convection_top(mean['conv'], z)
    zb, thb = z[::-1], mean['theta'][::-1]                                              # the ground up
    n = np.sqrt(np.maximum(g * np.diff(np.log(thb)) / np.diff(zb), 0.0))
    cum = np.r_[n[0] * zb[0], n[0] * zb[0] + np.cumsum(n * np.diff(zb))]              # from the ground to each level
    kt = int(np.argmin(np.abs(zb - top['top_m'])))
    return dict(run=run, years=list(years), lat_limit_deg=lat_limit_deg, records=count, convection_top_m=top['top_m'],
                c_m_s=float(cum[kt] / math.pi), c_to_model_top_m_s=float(cum[-1] / math.pi),
                z_m=zb.tolist(), theta_k=thb.tolist())


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('run')
    parser.add_argument('first', type=int)
    parser.add_argument('last', type=int)
    parser.add_argument('lat', type=float)
    parser.add_argument('lon', type=float)
    parser.add_argument('surface', choices=('sea', 'land'))
    parser.add_argument('--half-lat', type=float, default=5.6)
    parser.add_argument('--half-lon', type=float, default=8.5)
    args = parser.parse_args(argv)
    site = dict(lat_deg=args.lat, lon_deg=args.lon, half_lat_deg=args.half_lat, half_lon_deg=args.half_lon,
                surface=args.surface)
    f = site_forcing(args.run, (args.first, args.last), site)
    print(f"{len(f['cells'])} cells {f['cells']}, {f['outputs']} means at {f['phases']} phases of the day")
    print(f"surface {f['surface']}")
    b = f['budget']
    print('sigma  z km  w mm/s | heat K/d: circ  phys  mean-w  table | vapour g/kg/d: circ  phys  mean-w  table')
    for k in range(len(f['sigma']) - 1, -1, -1):
        print(f"{f['sigma'][k]:6.3f} {f['profile']['z_m'][k] / 1000:5.1f} {f['vertical_wind_m_s'][k] * 1000:6.2f} | "
              f"{b['circulation_heating_theta_k_day'][k]:7.3f} {b['physics_heating_theta_k_day'][k]:6.3f} "
              f"{b['mean_wind_advection_theta_k_day'][k]:7.3f} {b['table_mean_theta_k_day'][k]:6.3f} | "
              f"{b['circulation_moistening_g_kg_day'][k]:7.3f} {b['physics_moistening_g_kg_day'][k]:6.3f} "
              f"{b['mean_wind_advection_qv_g_kg_day'][k]:7.3f} {b['table_mean_qv_g_kg_day'][k]:6.3f}")
    c = b['column_vapour_mm_day']
    print(f"column vapour brought by the circulation: exact {c['exact']:.2f} mm/day, from the snapshots "
          f"{c['from_snapshots']:.2f}; rain {b['rain_mm_day']:.2f}, evaporation {b['evaporation_mm_day']:.2f}")
    print(f"column physics heating {b['column_heating_w_m2']['physics']:.1f} W/m2, of it the rain's latent heat "
          f"{b['column_heating_w_m2']['rain_latent']:.1f}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
