"""The storms a dryland point receives: storm libraries, each a set of point hyetographs with weights, and the
disaggregation of a GCM cell's 3-day mean of convective rain into an expected number of them.

A library is the input that sets the storms' intensities, so a later product can replace it. Two exist now and bracket
the intensities:

- walnut_gulch(): Earth's measured point storms at the USDA's Walnut Gulch watershed in south-eastern Arizona, as a
  parametric family fitted to measured point statistics. Each storm rises linearly to its peak over RISE_MIN and then
  decays exponentially with the time constant tau; storm depths are exponential with mean mu. mu and tau are fitted to
  Osborn's (1983) seasonal counts and durations of point rain above 50, 100, 150 and 200 mm/h (his Table 1) for a
  season of 203 mm of thunderstorm rain, two thirds of the 305 mm year (Mendez et al. 2003). The annual maxima of the
  5-, 15-, 30- and 60-minute intensities at the 2- and 10-year return periods (Mendez et al. 2003, Tables 1-4, the
  means of their gauge groups) are held out of the fit and test it.
- from_storm_lives(path): any CM1 storm-tracker product (schema terluna.climate.crm-storm-lives/1). Each whole track's
  rain beneath its core, the largest column rate at each output, is one storm. The 6 km window of box_0e_elec_corrected
  is the floor: its columns average 36 km2 and its outputs ten minutes, and it stands on wet ground at the equator. The
  dry box's windows replace it through the same reader when they run.

A cell's 3-day mean brings D mm of convective rain to every point on average; with a library whose storms average mu mm
a point receives a Poisson number of storms with mean D/mu, drawn from the library by weight.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

STEP_MIN = 1.0                          # every library's hyetographs are resampled to one-minute steps
RISE_MIN = 5.0                          # design choice: the Walnut Gulch family's rise to its peak
SEASON_MM = 2.0 / 3.0 * 305.0           # Walnut Gulch's thunderstorm rain a year (Mendez et al. 2003)
OSBORN_1983 = dict(threshold_mm_h=(50.0, 100.0, 150.0, 200.0), events_per_season=(4.6, 0.94, 0.29, 0.06),
                   minutes_per_event=(7.0, 4.3, 3.0, 2.4), min_minutes=2.0)
MENDEZ_2003 = dict(duration_min=(5.0, 15.0, 30.0, 60.0), two_year_mm_h=(91.9, 62.6, 41.7, 23.9),
                   ten_year_mm_h=(146.3, 95.3, 66.0, 39.0))
BINS = 48                               # equal-probability depth bins of the fitted family
SCHEMA = 'terluna.research.storm-water-storm-library/1'


def unit_shape(tau_min, step=0.05):
    """The family's hyetograph for 1 mm of rain on a fine grid: times (min) and rates (mm/min)."""
    t = np.arange(0.0, RISE_MIN + 12.0 * tau_min, step)
    r = np.where(t < RISE_MIN, t / RISE_MIN, np.exp(-(t - RISE_MIN) / tau_min))
    return t, r / (r.sum() * step)


def window_fraction(tau_min, minutes):
    """The largest share of a storm's depth falling within any window of the given length."""
    step = 0.05
    t, r = unit_shape(tau_min, step)
    c = np.concatenate(([0.0], np.cumsum(r) * step))
    n = int(round(minutes / step))
    return float(np.max(c[n:] - c[:-n])) if n < c.size else 1.0


def exceedance(mu, tau, threshold, min_minutes=OSBORN_1983['min_minutes']):
    """Per season of SEASON_MM: the expected number of storms whose rate exceeds the threshold (mm/h) for at least
    min_minutes, and their mean minutes above it, for exponential depths of mean mu (mm) and decay tau (min)."""
    d = np.linspace(0.0, 40.0 * mu, 20001)[1:]
    pdf = np.exp(-d / mu) / mu
    peak = d / ((RISE_MIN / 2.0 + tau) / 60.0)                                   # mm/h
    with np.errstate(divide='ignore', invalid='ignore'):
        dur = np.where(peak > threshold, RISE_MIN * (1.0 - threshold / peak) + tau * np.log(peak / threshold), 0.0)
    hit = dur >= min_minutes
    dd = d[1] - d[0]
    p = float((pdf * hit).sum() * dd)
    n = SEASON_MM / mu * p
    mean = float((pdf * hit * dur).sum() * dd / max(p, 1e-300))
    return n, mean


def annual_maximum(mu, tau, minutes, years):
    """The rate (mm/h) of the annual maximum over the given window at a return period (years): the largest storm of a
    Poisson season of SEASON_MM / mu exponential storms, its quantile 1 - 1/years, times the window's share."""
    q = 1.0 - 1.0 / years
    largest = mu * (np.log(SEASON_MM / mu) - np.log(-np.log(q)))
    return window_fraction(tau, minutes) * largest / (minutes / 60.0)


def fit_misfit(mu, tau):
    o = OSBORN_1983
    res = []
    for x, n_obs, t_obs in zip(o['threshold_mm_h'], o['events_per_season'], o['minutes_per_event']):
        n, t = exceedance(mu, tau, x)
        res += [np.log(max(n, 1e-9) / n_obs), np.log(max(t, 1e-9) / t_obs)]
    return np.array(res)


def fit_walnut_gulch():
    """mu (mm) and tau (min) minimising the squared log misfit to Osborn's eight numbers: a grid, then refined."""
    best = None
    for mu in np.linspace(3.0, 16.0, 27):
        for tau in np.linspace(4.0, 40.0, 37):
            e = float((fit_misfit(mu, tau) ** 2).sum())
            if best is None or e < best[0]:
                best = (e, mu, tau)
    from scipy.optimize import minimize
    r = minimize(lambda v: float((fit_misfit(*v) ** 2).sum()), [best[1], best[2]], method='Nelder-Mead',
                 options=dict(xatol=1e-4, fatol=1e-8))
    return float(r.x[0]), float(r.x[1])


def walnut_gulch():
    """The fitted family as a library: BINS storms at the conditional mean depths of equal-probability bins."""
    mu, tau = fit_walnut_gulch()
    q = np.linspace(0.0, 1.0, BINS + 1)
    a, b = 1.0 - q[:-1], 1.0 - q[1:]
    with np.errstate(divide='ignore', invalid='ignore'):
        part = np.where(b > 0, b * (1.0 - np.log(np.where(b > 0, b, 1.0))), 0.0)
    depth = mu * (a * (1.0 - np.log(a)) - part) / (q[1:] - q[:-1])
    t, r = unit_shape(tau, 0.05)
    per_min = np.add.reduceat(r * 0.05, np.arange(0, t.size, int(round(STEP_MIN / 0.05))))
    per_min = per_min[: np.nonzero(np.cumsum(per_min) < 0.9995)[0][-1] + 2]
    per_min /= per_min.sum()
    storms = [d * per_min * 60.0 / STEP_MIN for d in depth]                     # mm/h each minute
    fit = dict(osborn=[dict(threshold_mm_h=x, events=exceedance(mu, tau, x)[0], minutes=exceedance(mu, tau, x)[1],
                            measured_events=n, measured_minutes=m)
                       for x, n, m in zip(OSBORN_1983['threshold_mm_h'], OSBORN_1983['events_per_season'],
                                          OSBORN_1983['minutes_per_event'])],
               mendez=[dict(minutes=k, two_year=annual_maximum(mu, tau, k, 2.0), ten_year=annual_maximum(mu, tau, k, 10.0),
                            measured_two_year=a2, measured_ten_year=a10)
                       for k, a2, a10 in zip(MENDEZ_2003['duration_min'], MENDEZ_2003['two_year_mm_h'],
                                             MENDEZ_2003['ten_year_mm_h'])])
    return dict(schema=SCHEMA, name='walnut_gulch', step_min=STEP_MIN, storms=storms, weights=np.full(BINS, 1.0 / BINS),
                parameters=dict(mean_depth_mm=mu, decay_min=tau, rise_min=RISE_MIN, season_mm=SEASON_MM),
                fit=fit, source='Osborn (1983), WRR 19(4):1036-1042, Table 1 (fitted); Mendez, Goodrich and Osborn '
                                '(2003), JAWRA 39(3):611-621, Tables 1-4 (held out)',
                evidence='A parametric family fitted to measured point statistics of Earth\'s thunderstorms at Walnut '
                         'Gulch; Earth\'s storms, not the Moon\'s')


def from_storm_lives(path, root=None):
    """Each whole track of a CM1 storm-tracker product as a storm: its rain beneath the core at each output."""
    p = Path(path) if root is None else Path(root) / path
    raw = p.read_bytes()
    d = json.loads(raw)
    if d.get('schema') != 'terluna.climate.crm-storm-lives/1':
        raise ValueError(f'{path}: schema {d.get("schema")} is not a CM1 storm-lives product')
    rep = int(round(float(d['output_s']) / 60.0 / STEP_MIN))
    storms = [np.repeat(np.asarray(t['series']['rain_mm_h'], float), rep) for t in d['tracks'] if t['whole']]
    storms = [s for s in storms if s.sum() > 0.0]
    return dict(schema=SCHEMA, name=d['case'], step_min=STEP_MIN, storms=storms,
                weights=np.full(len(storms), 1.0 / len(storms)), parameters=dict(output_s=float(d['output_s'])),
                source=f'{path} (sha256 {hashlib.sha256(raw).hexdigest()})',
                evidence='CM1 model result: the rain beneath each whole tracked core, its largest 6 km column at each '
                         'ten-minute output; a floor on point intensities')


def depths(lib):
    return np.array([s.sum() * lib['step_min'] / 60.0 for s in lib['storms']])


def statistics(lib):
    """The library's storms: mean depth, duration with rain, peak and their weighted percentiles."""
    w = np.asarray(lib['weights'], float)
    d = depths(lib)
    dur = np.array([(s > 0.1).sum() * lib['step_min'] for s in lib['storms']])
    peak = np.array([s.max() for s in lib['storms']])

    def pct(x, q):
        o = np.argsort(x)
        c = np.cumsum(w[o]) / w.sum()
        return float(np.interp(q, c, x[o]))
    return dict(storms=len(d), mean_depth_mm=float((w * d).sum() / w.sum()),
                depth_mm_p50_p90=[pct(d, 0.5), pct(d, 0.9)], minutes_raining_mean=float((w * dur).sum() / w.sum()),
                peak_mm_h_p50_p90_max=[pct(peak, 0.5), pct(peak, 0.9), float(peak.max())],
                rain_weighted_peak_mm_h=float((w * d * peak).sum() / (w * d).sum()))


def storm_counts(convective_mm, lib):
    """The expected number of storms at a point for each mean's convective rain (mm)."""
    return np.asarray(convective_mm, float) / statistics(lib)['mean_depth_mm']
