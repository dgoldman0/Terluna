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


def site_forcing(run, years, site, window_deg=36.0, top_m=TOP_M):
    """The GCM's column at a site over the years: its mean air and surface, its mean vertical wind, the table CM1
    prescribes by local time, and the budgets behind it. Arrays by level run top first, as PlaSim writes them."""
    import netCDF4
    from climate.crm.cm1_run import gcm_soil, planet, solar_day_s
    from climate.gcm import year_files as yf
    p = planet()
    gravity, radius = p['gravity_m_s2'], p['radius_m']
    folder = HERE.parent / 'gcm' / 'runs' / run / 'model'
    keep = {k: [] for k in ('hour', 'theta', 'q', 'u', 'v', 'z', 'rho', 'omega', 'heat', 'moist', 'ps', 'ts', 'tas',
                            'rain', 'evaporation', 'water', 'wetness', 'conv', 'lsp')}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            missing = [k for k in MEANS + SNAPSHOTS + SCHEMES if k not in d.variables]
            if missing:
                raise ValueError(f'{run} year {year} lacks {missing}: run it with --heating-diagnostics '
                                 '--moisture-diagnostics')
            f = yf.fields(d, MEANS + SNAPSHOTS + SCHEMES)            # the humidity's means: gcm/year_files.py
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
    phase = np.round(((centre + 180.0) % 360.0) / window_deg - 0.5).astype(int) % int(round(360.0 / window_deg))
    return dict(
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
