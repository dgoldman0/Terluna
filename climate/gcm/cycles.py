"""Weather cycles beyond the lunar day in an ExoPlaSim run's 3-day records, tested against red noise with the
family-wise error controlled (8 November 2026).

    climate/gcm/.venv/bin/python -m climate.gcm.cycles A28R_301k_design_l15_snow_s924_turned_refit:10-49 \
        --compare 10-28 --surrogates 10000 --work <folder for the series>   # results/gcm/cycles_<run>_y10-49.json

The screen reruns the bounded screen of the turned design run's years 0-28 made on 7 November (`plasim-cycle-screen` in
`session-tools/codex_terrain_axis_2026-11-07` on the research drive) on the chosen years, with the significance test that
screen lacked. Its series are the same: area means (weighted by the cosine of latitude) of six fields over 39 regions in
the run's own frame, the globe, its land, its sea and 36 boxes of 30 degrees of latitude by 60 of longitude, at the year
files' 3-day means (143 steps of 30 minutes, 2.979 days, 120 to a file year of 357.5 days). The fields are the 2 m air
temperature, the rain (mm/day), the surface pressure, the 10 m winds and the sunlight arriving at the top of the air
(rst - rsut, since PlaSim writes the outgoing flux as negative). That sunlight is the orbit's and the rotation's, with
the twilight scheme's sunlight beyond the terminator, which follows the cloud above, so its series show the lines the
setup itself makes over a little weather noise. The preprocessing is the first screen's: each series loses its mean, a
linear trend and the first four harmonics of the solar day; for the residual screen also the first four harmonics of the
365.25-day forcing year and the solar harmonics' sidebands one and two forcing years away.

The test. Each residual's spectrum is the first screen's Welch estimate (segments of 768 records, 2,288 days, half
overlapping, each detrended), divided by the spectrum of the AR(1) process with the residual's lag-one correlation and
variance. The statistic of a series is the largest of these ratios over periods of 8-2,200 days, away from the solar day
and its second and third harmonics by 4%. Its distribution under red noise comes from AR(1) surrogates at the series'
correlation and variance, each put through the same regression, spectrum and fit, so a peak's p-value, the share of
surrogate maxima at or above its ratio, allows for the search over frequency (the single-step max-statistic of Westfall
and Young, 1993). Holm's step-down (1979) then controls the family-wise error across the 234 series: it holds under any
dependence between them, and the series here depend strongly on each other, since the globe, land and sea contain the
boxes and the six fields share the weather. Benjamini and Hochberg's false discovery rate is given beside it as the more
permissive reading. A peak within a series that Holm keeps is significant on its own when its p-value times the number
of series is at most the level: with a continuum of frequencies no series ever has all its hypotheses rejected, so the
step-down gains nothing over Bonferroni at the level of single peaks.

Lines at known frequencies (the forcing year's harmonics and the solar day's yearly sidebands) are tested by their
fitted amplitudes against the same surrogates' amplitudes, with Holm across the series for each line. The halves of the
window are screened alike, and for each candidate the amplitude and phase (in the run's time, from the start of year 0)
at its frequency are compared between them: a forced line keeps its phase, and the difference of the two halves'
complex amplitudes against the red noise's spread at that frequency gives a chi-square test of one line in both (for a
peak the search found, whose frequency the whole window fixes, that test leans toward agreement, and the halves' own
screens judge persistence). The
file year's comb, which restarts at each year file could leave, is tested as the variance of the residual's 120-record
composite against the surrogates'.

A second null answers what AR(1) cannot. The series remember far longer than their lag-one correlation implies (the
seas' 50 m mixed layers integrate the weather), so against AR(1) the longest periods stand out in many series at once.
Each residual's autocorrelation over a year (lags of 1-120 records) is fitted as the sum of a fast and a slow AR(1), and
surrogates of that sum go through the same regression, spectrum and ratio, so a peak that passes this null as well stands
above the series' own red continuum.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy import signal

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-cycles/1'
FIELDS = ('tas', 'pr', 'ps', 'uas', 'vas', 'incoming')
UNITS = dict(tas='K', pr='mm/day', ps='hPa', uas='m/s', vas='m/s', incoming='W/m2')
RECORDS = 120                       # records in a year file
STEPS = 143                         # model steps of 30 minutes in a record (NSTPW)
DT = STEPS * 30 / 1440              # days in a record, 2.979
FILE_YEAR = RECORDS * DT            # 357.5 days
NPERSEG, NOVERLAP = 768, 384        # the first screen's Welch segments of 2,288 days, half overlapping
BAND = (8.0, 2200.0)                # periods searched, days
GUARD = 0.04                        # relative distance kept from the solar day and its second and third harmonics
RHO_CLIP = (-0.98, 0.995)
ALPHA = 0.05
CANDIDATE_BANDS = {'460-520 d': (460.0, 520.0), '700-1100 d': (700.0, 1100.0)}
MEMORY_LAGS = 120                   # records of autocorrelation the two-timescale null fits, a file year
EVIDENCE = ('Spectral tests of one model realization: ExoPlaSim T21 3-day regional means over the stated years, against '
            'AR(1) red noise fitted to each series, with the family-wise error controlled across frequencies by surrogate '
            'maxima and across the series by Holm. A significant peak is a departure from AR(1) noise in this run. The '
            'attribution tags are rules applied to computed facts (the forcing catalogue, the halves\' phase test, the '
            'sunlight series); the written attribution lives with the climate README. The seas are 50 m mixed layers '
            'without currents, so the run\'s internal variability on long periods is a slab ocean\'s; nothing here is an '
            'observation of the Moon or a forecast.')
READING_RULE = ('Periods in days, from records 2.979 days apart. Amplitudes are half the peak-to-trough range of a fitted '
                'sinusoid in the field\'s units (tas K, pr mm/day, ps hPa, uas and vas m/s, incoming W/m2); phases are in '
                'degrees of y = A cos(2 pi t / period - phase), t in days from the start of model year 0 (file years of '
                '357.5 days). Regions are latitude a:b N, longitude c:d E in the turned frame. p_series is the within-series '
                'max-statistic p-value against AR(1) surrogates (floor 1/(surrogates + 1)); series_holm is the series\' '
                'Holm-adjusted p-value over all series; p_bonferroni is p_series times the number of series; significant '
                'means both at most 0.05. windows: "main" is the stated years, "compare" the first screen\'s years 10-28, '
                'the halves split the main window. forced_lines give each known line\'s significant series (Holm per line) '
                'and, for the strongest, the halves\' amplitudes and phases with p_same_line, large when one line fits both '
                'halves; for a peak the search found, the frequency is fitted to the whole window, which favours one where the '
                'halves agree, so there p_same_line leans toward agreement and the halves\' own ratios and p-values are '
                'the test of persistence. candidates settle the first screen\'s named features. persistence_of_main_peaks carries each '
                'significant peak of the main window into the halves with its attribution tag: forced (rule) when a line of '
                'the catalogue lies within half the record\'s resolution, p_same_line >= 0.01 and the sunlight of some '
                'region has a significant peak in the same bin; otherwise continuum (rule) when the two-timescale null '
                'explains it (p_memory_bonferroni > 0.05), internal (rule) when no line lies near and p_same_line < 0.01, '
                'undetermined otherwise. memory gives each series\' fitted fast share, fast and '
                'slow lag-one correlations and misfit, and the slow component\'s e-folding in days.')


def solar_day(rotation_days, year_days):
    """The solar (synodic) day from the sidereal rotation and the year."""
    return 1 / (1 / rotation_days - 1 / year_days)


def forcing_year_days(year_days, steps_per_day=48):
    """The year as the model steps it: a whole number of 30-minute steps, 365.25 days for the configured 365.25636."""
    return round(year_days * steps_per_day) / steps_per_day


def alias(f):
    """The frequency (per day) at which a line of frequency f appears in records DT apart, from 0 to the Nyquist."""
    fs = 1 / DT
    return np.abs(np.asarray(f, float) - np.round(np.asarray(f, float) / fs) * fs)


def regions(lat, lon, land):
    """The first screen's 39 regions as area weights (cells, regions), with their names: the globe, land, sea and boxes of
    30 degrees of latitude by 60 of longitude, longitudes east from -180."""
    w = np.cos(np.radians(lat))[:, None] * np.ones((1, len(lon)))
    w = w / w.sum()
    out = {'global': w, 'land': w * land, 'sea': w * ~land}
    lons = (lon + 180) % 360 - 180
    for a, b in [(-90, -60), (-60, -30), (-30, 0), (0, 30), (30, 60), (60, 90)]:
        for c in range(-180, 180, 60):
            mask = (lat[:, None] >= a) & (lat[:, None] < b) & (lons[None, :] >= c) & (lons[None, :] < c + 60)
            out[f'{a}:{b}N,{c}:{c + 60}E'] = w * mask
    out = {k: v / v.sum() for k, v in out.items() if v.sum() > 0}
    names = list(out)
    return names, np.stack([out[n].ravel() for n in names], axis=1)


def read(folder, first, last):
    """Regional series of the six fields over the year files first..last of a run folder (time first), the equatorial
    Sun's first longitude harmonic for the clock, and each year's digest as the first screen took it (SHA-256 of the six
    fields as float32, rain in mm/day, then czen)."""
    import netCDF4
    series = {k: [] for k in FIELDS}
    sun, digests = [], {}
    names = weights = None
    for year in range(first, last + 1):
        path = Path(folder) / 'model' / f'MOST.{year:05d}.nc'
        h = hashlib.sha256()
        with netCDF4.Dataset(path) as d:
            steps = np.asarray(d['time'][:])
            if len(steps) != RECORDS or np.any(np.diff(steps) != STEPS):
                raise ValueError(f'{path} does not hold {RECORDS} records {STEPS} steps apart')
            lat, lon = np.asarray(d['lat'][:]), np.asarray(d['lon'][:])
            land = np.asarray(d['lsm'][0]) > 0.5
            if names is None:
                names, weights = regions(lat, lon, land)
            for k in FIELDS:
                if k == 'incoming':
                    a = np.asarray(d['rst'][:] - d['rsut'][:], dtype='float32')
                else:
                    a = np.asarray(d[k][:], dtype='float32')
                if k == 'pr':
                    a *= 86400e3
                h.update(a.tobytes())
                series[k].append(a.reshape(RECORDS, -1) @ weights)
            cz = np.asarray(d['czen'][:], dtype='float32')
            h.update(cz.tobytes())
            sun.append(cz[:, np.abs(lat) < 20, :].mean(axis=1) @ np.exp(1j * np.radians(lon)))
        st = path.stat()
        digests[str(year)] = dict(fields_sha256=h.hexdigest(), bytes=st.st_size, mtime_ns=st.st_mtime_ns)
        print(f'read year {year}', flush=True)
    t = (np.arange(first * RECORDS, (last + 1) * RECORDS) + 0.5) * DT
    return dict(t=t, names=names, weights=weights, lat=lat, lon=lon, land=land, sun=np.concatenate(sun),
                digests=digests, **{k: np.concatenate(v) for k, v in series.items()})


def clock(t, sun):
    """The solar day the run keeps: the period of the equatorial Sun's longitude (the phase of czen's first longitude
    harmonic within 20 degrees of the equator) fitted over the records, with the fit's rms phase error in degrees."""
    phase = np.unwrap(np.angle(sun))
    fit = np.polyfit(t, phase, 1)
    err = phase - np.polyval(fit, t)
    return float(abs(2 * np.pi / fit[0])), float(np.degrees(np.std(err)))


def year_fit(t, insolation, solar, trial=None):
    """The year that fits the global sunlight's yearly swing best: for each trial year the rms residual of the mean, a
    trend, four solar harmonics and four harmonics of the trial year."""
    trial = np.arange(355.0, 370.0001, 0.01) if trial is None else np.asarray(trial, float)
    rms = []
    for year in trial:
        X = basis(t, [k / solar for k in range(1, 5)] + [k / year for k in range(1, 5)])
        coef = np.linalg.lstsq(X, insolation, rcond=None)[0]
        rms.append(float(np.std(insolation - X @ coef)))
    rms = np.asarray(rms)
    return float(trial[np.argmin(rms)]), rms


def frequencies(solar, annual, kind):
    """Frequencies (per day) the regression removes: 'solar', the solar day's first four harmonics; 'forced', with the
    forcing year's first four and the solar harmonics' sidebands one and two years away (the first screen's residual);
    'extended', with the file year's first four harmonics and the solar harmonics 5-40, which the records alias."""
    fs = [k / solar for k in range(1, 5)]
    if kind != 'solar':
        fs += [k / annual for k in range(1, 5)]
        fs += [k / solar + s * j / annual for k in range(1, 5) for j in (1, 2) for s in (-1, 1)]
    if kind == 'extended':
        fs += [k / FILE_YEAR for k in range(1, 5)]
        fs += [k / solar for k in range(5, 41)]
    return fs


def line_names(solar, annual):
    """Names and frequencies of the 'forced' basis's lines beyond the solar harmonics, in its column order."""
    names = [f'annual {k}' for k in range(1, 5)]
    names += [f'solar {k} {"+" if s > 0 else "-"} annual {j}' for k in range(1, 5) for j in (1, 2) for s in (-1, 1)]
    return names, frequencies(solar, annual, 'forced')[4:]


def basis(t, freqs):
    """Regression columns: a constant, a linear trend, and a sine and cosine at each frequency (in the run's time)."""
    cols = [np.ones(len(t)), (t - t.mean()) / np.ptp(t)]
    for f in freqs:
        cols += [np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t)]
    return np.stack(cols, axis=1)


def ratio_spectrum(r, dt=DT):
    """The first screen's spectrum of residuals r (time last): Welch's estimate, the lag-one correlation, the AR(1)
    spectrum of that correlation scaled to the residual's variance, and their ratio."""
    r = np.atleast_2d(r)
    f, p = signal.welch(r, fs=1 / dt, axis=-1, nperseg=NPERSEG, noverlap=NOVERLAP, detrend='linear')
    a = r[:, :-1] - r[:, :-1].mean(axis=1, keepdims=True)
    b = r[:, 1:] - r[:, 1:].mean(axis=1, keepdims=True)
    rho = np.clip((a * b).sum(axis=1) / np.sqrt((a * a).sum(axis=1) * (b * b).sum(axis=1)), *RHO_CLIP)
    red = (1 - rho[:, None] ** 2) / (1 + rho[:, None] ** 2 - 2 * rho[:, None] * np.cos(2 * np.pi * f[None] * dt))
    red *= (np.var(r, axis=1) / np.trapezoid(red, f, axis=1))[:, None]
    return f, p, rho, red, p / np.maximum(red, 1e-30)


def searched(f, solar):
    """The frequencies the screen searches: periods of 8-2,200 days, 4% away from the solar day and its harmonics 2 and 3."""
    with np.errstate(divide='ignore'):
        period = np.where(f > 0, 1 / np.where(f > 0, f, 1), np.inf)
    keep = (period > BAND[0]) & (period < BAND[1])
    for m in (1, 2, 3):
        keep &= np.abs(period - solar / m) / (solar / m) > GUARD
    return keep


def band_bins(f, lo, hi):
    """Welch bins whose period lies within lo-hi days, the band widened by half a bin on each side."""
    df = f[1] - f[0]
    return (f >= 1 / hi - df / 2) & (f <= 1 / lo + df / 2)


def ar1(rho, n, size, rng, sd=1.0):
    """size stationary AR(1) series of length n with lag-one correlation rho and standard deviation sd (time last)."""
    e = rng.standard_normal((size, n)) * (sd * np.sqrt(1 - rho * rho))
    e[:, 0] /= np.sqrt(1 - rho * rho)
    return signal.lfilter([1.0], [1.0, -rho], e, axis=1)


def file_year_stat(r):
    """Variance of a residual's composite over the 120 records of the file year, as a share of its variance (time last;
    the window starts at a year file's start)."""
    r = np.atleast_2d(r)
    k = r.shape[1] // RECORDS
    comp = r[:, :k * RECORDS].reshape(r.shape[0], k, RECORDS).mean(axis=1)
    comp = comp - comp.mean(axis=1, keepdims=True)
    return (comp ** 2).mean(axis=1) / r.var(axis=1)


def acf(r, lags):
    """Autocorrelation of a series at the given lags (records)."""
    r = np.asarray(r, float) - np.mean(r)
    d = float(np.dot(r, r))
    return np.array([np.dot(r[:-L], r[L:]) / d for L in lags])


def fit_memory(r, kmax=MEMORY_LAGS):
    """The sum of a fast and a slow AR(1) whose autocorrelation best fits the series' over lags 1-kmax records: the fast
    share of the variance, the fast and slow lag-one correlations, and the rms misfit of the autocorrelation."""
    from scipy import optimize
    lags = np.arange(1, kmax + 1)
    target = acf(r, lags)

    def misfit(q):
        w, fast, slow = q
        return w * fast ** lags + (1 - w) * slow ** lags - target
    best = None
    for start in ([0.7, 0.5, 0.97], [0.9, 0.2, 0.99], [0.5, -0.3, 0.95]):
        sol = optimize.least_squares(misfit, start, bounds=([0, -0.99, 0], [1, 0.99, 0.9995]))
        if best is None or sol.cost < best.cost:
            best = sol
    w, fast, slow = (float(v) for v in best.x)
    return w, fast, slow, float(np.sqrt(np.mean(best.fun ** 2)))


def two_ar1(w, fast, slow, n, size, rng, sd=1.0):
    """size series of a fast and a slow AR(1) summed (time last), the fast holding the share w of the variance sd**2."""
    return ar1(fast, n, size, rng, sd * np.sqrt(w)) + ar1(slow, n, size, rng, sd * np.sqrt(1 - w))


def memory_screen(obs, t, solar, annual, surrogates, rng, batch=1000, bands=CANDIDATE_BANDS, progress=None):
    """The screen's statistic under the two-timescale null: each series' fitted memory and the surrogate maxima of the
    same ratio over the searched band and each named band, through the same regression and spectrum (obs from observe)."""
    X = basis(t, frequencies(solar, annual, 'forced'))
    P = np.linalg.pinv(X)
    R, f, mask = obs['R'], obs['f'], obs['mask']
    bins = {k: band_bins(f, *v) & mask for k, v in bands.items()}
    S, n = R.shape
    out = dict(memory=np.empty((S, 4)), smax=np.empty((S, surrogates)),
               band_null={k: np.empty((S, surrogates)) for k in bands})
    start = time.time()
    for s in range(S):
        w, fast, slow, err = fit_memory(R[s])
        out['memory'][s] = (w, fast, slow, err)
        sd = float(R[s].std())
        for lo in range(0, surrogates, batch):
            b = min(batch, surrogates - lo)
            x = two_ar1(w, fast, slow, n, b, rng, sd)
            rat = ratio_spectrum(x - (x @ P.T) @ X.T)[4]
            out['smax'][s, lo:lo + b] = rat[:, mask].max(axis=1)
            for k, bb in bins.items():
                out['band_null'][k][s, lo:lo + b] = rat[:, bb].max(axis=1)
        if progress and (s + 1) % 39 == 0:
            print(f'{progress}: {s + 1}/{S} series, {time.time() - start:.0f} s', flush=True)
    return out


def memory_path(work, name, a, b):
    return Path(work) / f'memory_{name}_y{a}-{b}.npz'


def holm(p):
    """Holm's step-down adjusted p-values."""
    p = np.asarray(p, float)
    m = len(p)
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(np.argsort(p, kind='stable')):
        running = max(running, min(1.0, (m - rank) * p[i]))
        adj[i] = running
    return adj


def bh(p):
    """Benjamini and Hochberg's adjusted p-values (q-values)."""
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p, kind='stable')
    q = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


def sinusoid(r, t, f0, width):
    """The sinusoid that fits r best within f0 +- width/2: its frequency, amplitude and phase (y = A cos(2 pi f t - phase),
    in the run's time) and its complex amplitude cos + i sin."""
    grid = np.linspace(f0 - width / 2, f0 + width / 2, 81)
    grid = grid[grid > 0]
    power = np.abs(np.exp(-2j * np.pi * np.outer(grid, t)) @ (r - r.mean()))
    f = float(grid[np.argmax(power)])
    X = np.stack([np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t), np.ones(len(t))], axis=1)
    a, b, _ = np.linalg.lstsq(X, r, rcond=None)[0]
    return f, float(np.hypot(a, b)), float(np.degrees(np.arctan2(a, b))), complex(b, a)


def fit_line(r, t, f):
    """Complex amplitude (cos + i sin) of a sinusoid at frequency f fitted to r with a constant, in the run's time."""
    X = np.stack([np.sin(2 * np.pi * f * t), np.cos(2 * np.pi * f * t), np.ones(len(t))], axis=1)
    a, b, _ = np.linalg.lstsq(X, r, rcond=None)[0]
    return complex(b, a)


def screen(Y, t, solar, annual, surrogates, rng, batch=1000, bands=CANDIDATE_BANDS, progress=None):
    """Screen the series Y (time, series) on the times t against AR(1) red noise.

    Returns the residuals (series, time), the frequencies, the observed spectra, ratios, red spectra and correlations,
    the searched mask, each series' surrogate maxima (series, surrogates), the surrogate band maxima for each named band,
    the observed and surrogate amplitudes of the 'forced' lines, and the observed and surrogate file-year statistics."""
    out = observe(Y, t, solar, annual, bands)
    X, P, R, rho, mask, bins = out.pop('X'), out.pop('P'), out['R'], out['rho'], out['mask'], out.pop('bins')
    S, n = R.shape
    nlines = (X.shape[1] - 10) // 2
    out.update(smax=np.empty((S, surrogates)), amp_null=np.empty((S, nlines, surrogates)),
               fy_null=np.empty((S, surrogates)), band_null={k: np.empty((S, surrogates)) for k in bands})
    start = time.time()
    for s in range(S):
        sd = float(R[s].std())
        for lo in range(0, surrogates, batch):
            b = min(batch, surrogates - lo)
            x = ar1(float(rho[s]), n, b, rng, sd)
            c = x @ P.T
            res = x - c @ X.T
            rat = ratio_spectrum(res)[4]
            out['smax'][s, lo:lo + b] = rat[:, mask].max(axis=1)
            out['amp_null'][s, :, lo:lo + b] = np.hypot(c[:, 10::2], c[:, 11::2]).T
            out['fy_null'][s, lo:lo + b] = file_year_stat(res)
            for k, bb in bins.items():
                out['band_null'][k][s, lo:lo + b] = rat[:, bb].max(axis=1) if bb.any() else np.nan
        if progress and (s + 1) % 39 == 0:
            print(f'{progress}: {s + 1}/{S} series, {time.time() - start:.0f} s', flush=True)
    return out


def observe(Y, t, solar, annual, bands=CANDIDATE_BANDS):
    """The observed side of the screen: residuals of the 'forced' regression (series, time), spectra, ratios, the
    searched mask, the lines' fitted amplitudes, the file-year statistic and the band maxima."""
    X = basis(t, frequencies(solar, annual, 'forced'))
    P = np.linalg.pinv(X)
    coef = Y.T @ P.T
    R = Y.T - coef @ X.T
    f, p, rho, red, ratio = ratio_spectrum(R)
    mask = searched(f, solar)
    bins = {k: band_bins(f, *v) & mask for k, v in bands.items()}
    return dict(X=X, P=P, bins=bins, R=R, f=f, p=p, rho=rho, red=red, ratio=ratio, mask=mask,
                amp_obs=np.hypot(coef[:, 10::2], coef[:, 11::2]), fy_obs=file_year_stat(R),
                band_obs={k: ratio[:, b].max(axis=1) if b.any() else np.full(len(R), np.nan) for k, b in bins.items()})


def save_screen(path, sc, p_lines, labels, made_by, seed):
    """Keep a window's surrogate results on disk, so that the summary can be remade without recomputing them."""
    np.savez_compressed(path, smax=sc['smax'], fy_null=sc['fy_null'], p_lines=p_lines, made_by=made_by, seed=seed,
                        labels=np.array([f'{k} {r}' for k, r in labels]), f=sc['f'], ratio=sc['ratio'], red=sc['red'],
                        rho=sc['rho'], mask=sc['mask'], amp_obs=sc['amp_obs'],
                        **{'band_null ' + k: v for k, v in sc['band_null'].items()})


def load_screen(path, Y, t, solar, annual, bands=CANDIDATE_BANDS):
    """A window's screen from its kept surrogate results, the observed side recomputed and checked against them."""
    z = np.load(path, allow_pickle=False)
    sc = observe(Y, t, solar, annual, bands)
    for k in ('X', 'P', 'bins'):
        sc.pop(k)
    if not np.allclose(sc['ratio'], z['ratio'], rtol=1e-6, atol=1e-8 * float(np.abs(z['ratio']).max())):
        raise ValueError(f'{path} was screened from other series')
    sc.update(smax=z['smax'], fy_null=z['fy_null'], band_null={k: z['band_null ' + k] for k in bands})
    return sc, z['p_lines'], str(z['made_by'])


def exceed(null, obs):
    """Monte Carlo p-value: (1 + surrogates at or above the observed value) / (1 + surrogates), along the last axis."""
    null = np.asarray(null)
    return (1 + (null >= np.asarray(obs)[..., None]).sum(axis=-1)) / (1 + null.shape[-1])


def peaks(sc, labels, t, keep_p=0.2):
    """Each series' spectral peaks in the searched band with their within-series max-statistic p-values, keeping those
    at or below keep_p, and the series' own p-value (its largest searched ratio)."""
    f, mask, ratio = sc['f'], sc['mask'], sc['ratio']
    df = f[1] - f[0]
    S = ratio.shape[0]
    series_max = np.where(mask[None], ratio, -np.inf).max(axis=1)
    p_series = exceed(sc['smax'], series_max)
    rows = []
    for s in range(S):
        found, _ = signal.find_peaks(ratio[s])
        for i in found[mask[found]]:
            pv = float(exceed(sc['smax'][s], ratio[s, i]))
            if pv > keep_p:
                continue
            fr, amp, ph, _ = sinusoid(sc['R'][s], t, f[i], df)
            rows.append(dict(series=s, field=labels[s][0], region=labels[s][1], bin_period_days=float(1 / f[i]),
                             period_days=1 / fr, ratio=float(ratio[s, i]), p_series=pv, amplitude=amp, phase_deg=ph,
                             variance_share=amp ** 2 / 2 / float(sc['R'][s].var()), rho1=float(sc['rho'][s])))
    return p_series, series_max, rows


def catalogue(solar, annual, kmax=40, fmax=1 / BAND[0]):
    """Lines the setup makes, as they appear in the records: the solar harmonics 1-kmax and their sidebands one and two
    forcing years away, aliased; the forcing year's harmonics; the file year's harmonics."""
    lines = []
    for k in range(1, kmax + 1):
        for j in (-2, -1, 0, 1, 2):
            name = f'solar {k}' + ('' if j == 0 else f' {"+" if j > 0 else "-"} annual {abs(j)}')
            fa = float(alias(k / solar + j / annual))
            lines.append((name + (' (aliased)' if abs(k / solar + j / annual - fa) > 1e-12 else ''), fa))
    lines += [(f'annual {k}', k / annual) for k in range(1, 9)]
    lines += [(f'file year {k}', k / FILE_YEAR) for k in range(1, int(FILE_YEAR * fmax) + 1)]
    return [(n, f) for n, f in lines if 0 < f <= fmax * 1.05]


def nearest(f, lines):
    """The catalogue line nearest frequency f: its name, frequency and the offset f - line (per day)."""
    i = int(np.argmin([abs(f - fl) for _, fl in lines]))
    return lines[i][0], lines[i][1], f - lines[i][1]


def red_level(sc, s, f):
    """The fitted red spectrum (one-sided density) of series s at frequency f."""
    return float(np.interp(f, sc['f'], sc['red'][s]))


def same_line(c1, c2, var1, var2):
    """p-value that two windows' complex amplitudes c1 and c2 differ only by red noise of coefficient variances var1 and
    var2 (chi-square with two degrees of freedom)."""
    return float(np.exp(-abs(c1 - c2) ** 2 / (2 * (var1 + var2))))


def window(data, first, last):
    """The records of years first..last from a series dictionary that holds them."""
    t = data['t']
    keep = (t >= first * FILE_YEAR) & (t < (last + 1) * FILE_YEAR)
    Y = np.concatenate([data[k][keep] for k in FIELDS], axis=1)
    return t[keep], Y


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(args):
    from climate.gcm.exoplasim_run import RUNS
    name, years = args.run.split(':')
    first, last = (int(v) for v in years.split('-'))
    folder = (Path(args.runs) if args.runs else RUNS) / name
    progress_bytes = (folder / 'progress.json').read_bytes()
    planet = json.loads(progress_bytes)['configuration']['planet']
    solar = solar_day(planet['rotationperiod'], planet['year'])
    annual = forcing_year_days(planet['year'])
    work = Path(args.work)
    work.mkdir(parents=True, exist_ok=True)
    cache = work / f'series_{name}_y{first}-{last}.npz'
    if cache.exists() and not args.reread:
        z = np.load(cache, allow_pickle=False)
        data = {k: z[k] for k in z.files if k not in ('names', 'digests')}
        data['names'] = list(z['names'])
        data['digests'] = json.loads(str(z['digests']))
    else:
        data = read(folder, first, last)
        np.savez(cache, **{k: v for k, v in data.items() if k not in ('names', 'digests')},
                 names=np.array(data['names']), digests=json.dumps(data['digests']))
    names = data['names']
    labels = [(k, r) for k in FIELDS for r in names]
    S = len(labels)
    observed_solar, phase_rms = clock(data['t'], data['sun'])
    best_year, _ = year_fit(data['t'], data['incoming'][:, names.index('global')], solar)
    lines = catalogue(solar, annual)
    lnames, lfreqs = line_names(solar, annual)
    rng = np.random.default_rng(args.seed)

    windows = {'main': (first, last)}
    mid = (first + last + 1) // 2
    memory_windows = {'main': (first, last), 'second half': (mid, last)}
    if args.memory_only:
        mrng = np.random.default_rng([args.seed, 1])
        for key, (a, b) in memory_windows.items():
            t, Y = window(data, a, b)
            nsur = args.memory_surrogates if key == 'main' else args.memory_half_surrogates
            ms = memory_screen(observe(Y, t, solar, annual), t, solar, annual, nsur, mrng, progress=f'memory {key} {a}-{b}')
            np.savez_compressed(memory_path(work, name, a, b), memory=ms['memory'], smax=ms['smax'],
                                made_by=digest(Path(__file__))[:16], seed=args.seed,
                                **{'band_null ' + k: v for k, v in ms['band_null'].items()})
        return None
    if args.compare:
        a, b = (int(v) for v in args.compare.split('-'))
        windows['compare'] = (a, b)
    windows['first half'] = (first, mid - 1)
    windows['second half'] = (mid, last)
    screens, made_by, surrogates = {}, {}, {}
    for key, (a, b) in windows.items():
        t, Y = window(data, a, b)
        nsur = args.half_surrogates if key.endswith('half') and args.half_surrogates else args.surrogates
        path = work / f'screen_{name}_y{a}-{b}.npz'
        if args.from_screens:
            sc, p_lines, made_by[key] = load_screen(path, Y, t, solar, annual)
            nsur = sc['smax'].shape[1]
        else:
            sc = screen(Y, t, solar, annual, nsur, rng, progress=f'{key} {a}-{b}')
            p_lines = exceed(sc['amp_null'], sc['amp_obs'])             # (series, lines)
            made_by[key] = digest(Path(__file__))[:16]
            save_screen(path, sc, p_lines, labels, made_by[key], args.seed)
        surrogates[key] = int(nsur)
        p_series, series_max, rows = peaks(sc, labels, t)
        h, q = holm(p_series), bh(p_series)
        holm_lines = np.stack([holm(p_lines[:, j]) for j in range(p_lines.shape[1])], axis=1)
        p_fy = exceed(sc['fy_null'], sc['fy_obs'])
        p_band = {k: exceed(sc['band_null'][k], sc['band_obs'][k]) for k in CANDIDATE_BANDS}
        for row in rows:
            row['p_bonferroni'] = min(1.0, S * row['p_series'])
            row['series_holm'] = float(h[row['series']])
            row['series_bh'] = float(q[row['series']])
            row['significant'] = bool(row['series_holm'] <= ALPHA and row['p_bonferroni'] <= ALPHA)
            ln, lf, off = nearest(1 / row['period_days'], lines)
            row['nearest_forced_line'] = dict(name=ln, period_days=1 / lf, offset_in_resolution=off * (t[-1] - t[0] + DT))
        screens[key] = dict(sc=sc, t=t, Y=Y, rows=rows, p_series=p_series, holm=h, bh=q, series_max=series_max,
                            p_lines=p_lines, holm_lines=holm_lines, p_fy=p_fy, holm_fy=holm(p_fy), p_band=p_band)
        (work / f'peaks_{name}_y{a}-{b}.json').write_text(json.dumps(dict(surrogates=nsur, rows=rows), indent=0) + '\n')
    memory = {}
    mrng = np.random.default_rng([args.seed, 1])
    for key, (a, b) in memory_windows.items():
        w, path = screens[key], memory_path(work, name, a, b)
        if not path.exists():
            nsur = args.memory_surrogates if key == 'main' else args.memory_half_surrogates
            ms = memory_screen(w['sc'], w['t'], solar, annual, nsur, mrng, progress=f'memory {key} {a}-{b}')
            np.savez_compressed(path, memory=ms['memory'], smax=ms['smax'], made_by=digest(Path(__file__))[:16],
                                seed=args.seed, **{'band_null ' + k: v for k, v in ms['band_null'].items()})
        z = np.load(path, allow_pickle=False)
        ms = dict(memory=z['memory'], smax=z['smax'], band_null={k: z['band_null ' + k] for k in CANDIDATE_BANDS})
        made_by['memory ' + key] = str(z['made_by'])
        surrogates['memory ' + key] = int(ms['smax'].shape[1])
        ms['p_series'] = exceed(ms['smax'], w['series_max'])
        ms['holm'] = holm(ms['p_series'])
        ms['p_band'] = {k: exceed(ms['band_null'][k], w['sc']['band_obs'][k]) for k in CANDIDATE_BANDS}
        memory[key] = ms
    product = summarise(args, name, first, last, folder, progress_bytes, planet, solar, annual, observed_solar,
                        phase_rms, best_year, data, labels, lnames, lfreqs, lines, screens, windows, memory)
    product['method']['surrogates'] = (f'Surrogates per series and window ({surrogates}): AR(1) at the residual\'s '
                                       'correlation and variance for the screens, the fitted two-timescale sum for '
                                       '"memory", each through the same regression, spectrum and fit; seed '
                                       f'{args.seed} (the memory null on its own stream).')
    product['producer']['screens_made_by'] = {k: {'gcm/cycles.py': v} for k, v in made_by.items()}
    out = Path(args.out) if args.out else RESULTS / f'cycles_{name}_y{first}-{last}.json'
    out.write_text(json.dumps(product, indent=1) + '\n')
    print(f'wrote {out}')
    return product


def summarise(args, name, first, last, folder, progress_bytes, planet, solar, annual, observed_solar, phase_rms,
              best_year, data, labels, lnames, lfreqs, lines, screens, windows, memory):
    """The data product: the clock, each window's significant peaks and lines, the candidates of the first screen settled,
    the file year's comb, and each finding's attribution."""
    S = len(labels)
    main, halves = screens['main'], (screens['first half'], screens['second half'])
    T = {k: float(v['t'][-1] - v['t'][0] + DT) for k, v in screens.items()}

    def window_summary(key):
        w = screens[key]
        sig = [r for r in w['rows'] if r['significant']]
        sig.sort(key=lambda r: r['p_series'])
        groups = {}
        for r in sig:
            groups.setdefault(round(r['bin_period_days'], 1), []).append(r)
        lines_out = {}
        for j, ln in enumerate(lnames):
            keep = np.where(w['holm_lines'][:, j] <= ALPHA)[0]
            order = keep[np.argsort(-w['sc']['amp_obs'][keep, j])]
            lines_out[ln] = dict(period_days=1 / lfreqs[j], series_significant=int(len(keep)),
                                 by_field={k: int(sum(labels[s][0] == k for s in keep)) for k in FIELDS},
                                 largest=[dict(field=labels[s][0], region=labels[s][1],
                                               amplitude=float(w['sc']['amp_obs'][s, j]),
                                               p_series=float(w['p_lines'][s, j]), holm=float(w['holm_lines'][s, j]))
                                          for s in order[:6]])
        mem = {}
        if key in memory:
            m = memory[key]
            mem = dict(memory_null=dict(
                series_significant_holm=int((m['holm'] <= ALPHA).sum()),
                by_field={k: int(sum((m['holm'][s] <= ALPHA) for s in range(S) if labels[s][0] == k)) for k in FIELDS},
                significant=[dict(field=labels[s][0], region=labels[s][1], p_series=float(m['p_series'][s]),
                                  holm=float(m['holm'][s]), ratio=float(w['series_max'][s]),
                                  bin_period_days=float(1 / w['sc']['f'][np.where(w['sc']['mask'])[0][
                                      np.argmax(w['sc']['ratio'][s][w['sc']['mask']])]]))
                             for s in np.argsort(m['p_series']) if m['holm'][s] <= ALPHA][:40]))
        return dict(
            years=list(windows[key]), records=int(len(w['t'])), span_days=T[key], **mem,
            series_significant_holm=int((w['holm'] <= ALPHA).sum()),
            series_significant_holm_by_field={k: int(sum((w['holm'][s] <= ALPHA) for s in range(S) if labels[s][0] == k))
                                              for k in FIELDS},
            series_significant_bh=int((w['bh'] <= ALPHA).sum()),
            significant_peaks_by_period=[dict(
                bin_period_days=per, series=len(g), fields=sorted({r['field'] for r in g}),
                strongest=[{k: r[k] for k in ('field', 'region', 'period_days', 'ratio', 'p_series', 'p_bonferroni',
                                               'series_holm', 'amplitude', 'variance_share', 'nearest_forced_line')}
                           for r in sorted(g, key=lambda r: -r['ratio'])[:4]]) for per, g in sorted(groups.items())],
            lines=lines_out,
            file_year=dict(series_significant_holm=int((w['holm_fy'] <= ALPHA).sum()),
                           by_field={k: int(sum((w['holm_fy'][s] <= ALPHA) for s in range(S) if labels[s][0] == k))
                                     for k in FIELDS},
                           smallest=[dict(field=labels[s][0], region=labels[s][1], share=float(w['sc']['fy_obs'][s]),
                                          p_series=float(w['p_fy'][s]), holm=float(w['holm_fy'][s]))
                                     for s in np.argsort(w['p_fy'])[:5]]),
            surrogate_maximum_quantiles=dict(
                median_over_series={q: float(np.median(np.quantile(w['sc']['smax'], q, axis=1)))
                                    for q in (0.5, 0.95, 0.99, 0.999)}))

    half_coef = [np.linalg.lstsq(basis(h['t'], frequencies(solar, annual, 'forced')), h['Y'], rcond=None)[0]
                 for h in halves]       # each half's own fit of every line, in the run's time

    def line_halves(j, s):
        """One forced line in one series: amplitude and phase in each half, and the test of one line in both."""
        f = lfreqs[j]
        cs = [complex(hc[11 + 2 * j, s], hc[10 + 2 * j, s]) for hc in half_coef]
        var = [red_level(h['sc'], s, f) / T[k] for h, k in zip(halves, ('first half', 'second half'))]
        return dict(amplitude=[abs(c) for c in cs], phase_deg=[float(np.degrees(np.angle(c))) for c in cs],
                    noise_amplitude=[float(np.sqrt(v)) for v in var], p_same_line=same_line(cs[0], cs[1], *var))

    # The forced lines: in the main window, each line's strongest series and how the halves agree.
    forced_lines = {}
    for j, ln in enumerate(lnames):
        keep = np.where(main['holm_lines'][:, j] <= ALPHA)[0]
        same = np.array([line_halves(j, s)['p_same_line'] for s in keep])
        strongest = {}
        for k in FIELDS:
            mine = [s for s in keep if labels[s][0] == k]
            top = sorted(mine, key=lambda s: -main['sc']['amp_obs'][s, j])[:2]
            strongest[k] = [dict(region=labels[s][1], amplitude=float(main['sc']['amp_obs'][s, j]),
                                 holm=float(main['holm_lines'][s, j]), halves=line_halves(j, s)) for s in top]
        forced_lines[ln] = dict(period_days=1 / lfreqs[j], series_significant=int(len(keep)),
                                by_field={k: int(sum(labels[s][0] == k for s in keep)) for k in FIELDS},
                                halves_one_line=int((same >= 0.01).sum()) if len(keep) else 0,
                                strongest_by_field=strongest)

    # The first screen's named candidates.
    idx = {lab: i for i, lab in enumerate(labels)}

    def band_case(field, region, band):
        s = idx[(field, region)]
        lo, hi = CANDIDATE_BANDS[band]
        out = {}
        for key, w in screens.items():
            f, ratio = w['sc']['f'], w['sc']['ratio'][s]
            bb = band_bins(f, lo, hi) & w['sc']['mask']
            i = np.where(bb)[0][np.argmax(ratio[bb])]
            pf, pp = signal.periodogram(w['sc']['R'][s], fs=1 / DT, window='hann', detrend='linear')
            sel = (pf >= 1 / hi * 0.75) & (pf <= 1 / lo * 1.25)
            half_peak = float(1 / pf[np.where(sel)[0][np.argmax(pp[sel])]])
            fr, amp, ph, _ = sinusoid(w['sc']['R'][s], w['t'], f[i], f[1] - f[0])
            extra = {}
            if key in memory:
                extra = dict(p_band_memory_null=float(memory[key]['p_band'][band][s]),
                             holm_band_memory_null_across_series=float(holm(memory[key]['p_band'][band])[s]))
            out[key] = dict(years=list(windows[key]), bin_period_days=float(1 / f[i]), ratio=float(ratio[i]), **extra,
                            p_band=float(w['p_band'][band][s]), p_series_search=float(exceed(w['sc']['smax'][s], ratio[i])),
                            holm_band_across_series=float(holm(w['p_band'][band])[s]),
                            periodogram_peak_days=half_peak, fitted_period_days=1 / fr, amplitude=amp, phase_deg=ph,
                            variance_share=amp ** 2 / 2 / float(w['sc']['R'][s].var()))
        allp = {key: dict(series_holm_significant=int((holm(w['p_band'][band]) <= ALPHA).sum()),
                          smallest=[dict(field=labels[s2][0], region=labels[s2][1], p_band=float(w['p_band'][band][s2]),
                                         holm=float(holm(w['p_band'][band])[s2]))
                                    for s2 in np.argsort(w['p_band'][band])[:5]]) for key, w in screens.items()}
        return dict(field=field, region=region, band_days=[lo, hi], windows=out, across_series=allp)

    def named(ln, cases):
        j = lnames.index(ln)
        return [dict(field=k, region=r, amplitude=float(main['sc']['amp_obs'][idx[(k, r)], j]),
                     holm=float(main['holm_lines'][idx[(k, r)], j]), halves=line_halves(j, idx[(k, r)]))
                for k, r in cases]

    first_screen_annual = [('tas', '-90:-60N,60:120E'), ('pr', '0:30N,-180:-120E'), ('tas', '30:60N,-120:-60E')]
    candidates = {
        'annual cycle and harmonics': {ln: dict(forced_lines[ln], first_screen_regions=named(ln, first_screen_annual))
                                       for ln in lnames[:4]},
        '27.32 d (solar + annual)': dict(forced_lines['solar 1 + annual 1'], first_screen_regions=named(
            'solar 1 + annual 1', [('tas', '-90:-60N,60:120E'), ('pr', '-30:0N,0:60E'), ('incoming', '-60:-30N,-60:0E')])),
        '32.13 d (solar - annual)': dict(forced_lines['solar 1 - annual 1'], first_screen_regions=named(
            'solar 1 - annual 1', [('tas', '30:60N,-60:0E'), ('pr', '-30:0N,-180:-120E'), ('incoming', '-60:-30N,-60:0E')])),
        '460-520 d rain': band_case('pr', '0:30N,-180:-120E', '460-520 d'),
        '700-1100 d temperature': band_case('tas', '30:60N,-120:-60E', '700-1100 d'),
    }

    # The extended regression: how the known lines' amplitudes move when the aliased harmonics and the file year join.
    t, Y = main['t'], main['Y']
    Xf, Xe = basis(t, frequencies(solar, annual, 'forced')), basis(t, frequencies(solar, annual, 'extended'))
    cf, ce = np.linalg.lstsq(Xf, Y, rcond=None)[0], np.linalg.lstsq(Xe, Y, rcond=None)[0]
    extended = {}
    for j, ln in enumerate(lnames[:8]):
        af, ae = np.hypot(cf[10 + 2 * j], cf[11 + 2 * j]), np.hypot(ce[10 + 2 * j], ce[11 + 2 * j])
        s = int(np.argmax(np.where(main['holm_lines'][:, j] <= ALPHA, af, -1)))
        extended[ln] = dict(field=labels[s][0], region=labels[s][1], amplitude_forced_basis=float(af[s]),
                            amplitude_extended_basis=float(ae[s]))
    resid_f = Y - Xf @ cf
    resid_e = Y - Xe @ ce
    inc = [i for i, lab in enumerate(labels) if lab[0] == 'incoming']
    extended['incoming residual rms, forced to extended basis (W/m2), median over regions'] = [
        float(np.median(resid_f[:, inc].std(axis=0))), float(np.median(resid_e[:, inc].std(axis=0)))]

    # Persistence of each significant peak of the main window in the halves.
    persistence = []
    significant = sorted([r for r in main['rows'] if r['significant']], key=lambda r: (r['p_series'], -r['ratio']))
    chosen = [r for r in significant if r['field'] != 'incoming'][:120] + [r for r in significant
                                                                          if r['field'] == 'incoming'][:12]
    for r in chosen:
        s, f = r['series'], 1 / r['period_days']
        hs = []
        cs, var = [], []
        for h, key in zip(halves, ('first half', 'second half')):
            fh, ratio = h['sc']['f'], h['sc']['ratio'][s]
            i = int(np.argmin(abs(fh - f)))
            c = fit_line(h['sc']['R'][s], h['t'], f)
            cs.append(c)
            var.append(red_level(h['sc'], s, f) / T[key])
            hs.append(dict(ratio=float(ratio[i]), p_series=float(exceed(h['sc']['smax'][s], ratio[i])),
                           amplitude=abs(c), phase_deg=float(np.degrees(np.angle(c)))))
        p_same = same_line(cs[0], cs[1], *var)
        match = abs(r['nearest_forced_line']['offset_in_resolution']) <= 0.5
        p_mem = float(exceed(memory['main']['smax'][s], r['ratio']))
        sun = [o for o in main['rows'] if o['significant'] and o['field'] == 'incoming'
               and abs(o['bin_period_days'] - r['bin_period_days']) < 1e-6]
        in_sun = any(o['region'] == r['region'] for o in sun)
        tag = ('forced (rule)' if match and p_same >= 0.01 and sun else
               'continuum (rule)' if min(1.0, S * p_mem) > ALPHA else
               'internal (rule)' if not match and p_same < 0.01 else 'undetermined (rule)')
        persistence.append(dict(field=r['field'], region=r['region'], period_days=r['period_days'],
                                bin_period_days=r['bin_period_days'], ratio=r['ratio'], p_bonferroni=r['p_bonferroni'],
                                amplitude=r['amplitude'], variance_share=r['variance_share'],
                                nearest_forced_line=r['nearest_forced_line'], halves=hs, p_same_line=p_same,
                                p_memory_null=p_mem, p_memory_bonferroni=min(1.0, S * p_mem),
                                in_sunlight_same_region=bool(in_sun), in_sunlight_regions=len(sun), attribution=tag))

    inputs = {'run': name, 'folder': str(folder), 'years': [first, last],
              'progress_json_sha256': hashlib.sha256(progress_bytes).hexdigest(),
              'year_files': data['digests'],
              'digest_rule': 'SHA-256 of tas, pr (mm/day), ps, uas, vas, rst - rsut and czen as float32, in that order, '
                             'as the first screen recorded them'}
    return dict(
        schema=SCHEMA, run=name, frame='turned (the run\'s own grid)', years=[first, last],
        clock=dict(record_days=DT, file_year_days=FILE_YEAR, forcing_year_days=annual,
                   forcing_year_best_fit_days=best_year, solar_days_config=solar, solar_days_observed=observed_solar,
                   solar_phase_rms_degrees=phase_rms, sampling_over_solar=solar / DT,
                   rotation_days=planet['rotationperiod'], year_days_config=planet['year'],
                   eccentricity=planet['eccentricity'], obliquity_deg=planet['obliquity']),
        series=dict(fields=list(FIELDS), units=UNITS, regions=data['names'], count=len(labels)),
        method=dict(
            preprocessing='Mean, linear trend and solar harmonics 1-4 removed by least squares; the residual screen also '
                          'removes the forcing year\'s harmonics 1-4 and the solar harmonics\' sidebands one and two '
                          'forcing years away (the first screen\'s basis).',
            spectrum=f'Welch, {NPERSEG}-record segments ({NPERSEG * DT:.0f} d), {NOVERLAP} overlap, linear detrend per '
                     'segment; ratio to the AR(1) spectrum of the residual\'s lag-one correlation (clipped to '
                     f'{RHO_CLIP}) scaled to its variance.',
            search=f'Periods {BAND[0]:g}-{BAND[1]:g} d, {GUARD:.0%} away from the solar day and its harmonics 2 and 3.',
            surrogates='',
            within_series='Single-step max-statistic p-value: share of surrogate maxima at or above a peak\'s ratio '
                          '(Westfall and Young 1993).',
            across_series=f'Holm (1979) step-down over the {S} series at {ALPHA}, valid under any dependence (the regions '
                          'nest and the fields share the weather); a peak is significant when its series passes Holm '
                          'and its p-value times the number of series is at most the level. Benjamini and Hochberg '
                          '(1995) given as the permissive reading.',
            memory_null='Two-timescale red noise: each residual\'s autocorrelation over lags 1-120 records fitted as '
                        'the sum of a fast and a slow AR(1); surrogates of the sum through the same regression, spectrum '
                        'and ratio give p-values under a null that holds the series\' long memory.',
            lines='Known lines tested by fitted amplitude against the surrogates\' amplitudes, Holm across series per line.',
            halves='Amplitude and phase at the line\'s frequency in each half (run time from the start of year 0); '
                   'p_same_line is the chi-square (2 dof) test that the halves\' complex amplitudes differ only by the '
                   'fitted red noise at that frequency.',
            red_noise_reference='AR(1) red-noise spectrum: Gilman, Fuglister and Mitchell (1963).'),
        windows={k: window_summary(k) for k in screens},
        forced_lines=forced_lines, candidates=candidates, extended_basis=extended,
        persistence_of_main_peaks=persistence,
        memory={k: dict(window=list(memory_windows_of(windows)[k]), columns=['fast share', 'fast rho1', 'slow rho1',
                                                                          'acf rms misfit'],
                        slow_efolding_days={f'{labels[s][0]} {labels[s][1]}': (float(-DT / np.log(m['memory'][s, 2]))
                                                                             if 0 < m['memory'][s, 2] < 1 else None)
                                            for s in range(S)},
                        fits={f'{labels[s][0]} {labels[s][1]}': [round(float(v), 5) for v in m['memory'][s]]
                              for s in range(S)}) for k, m in memory.items()},
        inputs=inputs, evidence=EVIDENCE, reading_rule=READING_RULE,
        producer=dict(domain='climate', files={'gcm/cycles.py': digest(Path(__file__))[:16]},
                      command=' '.join(sys.argv), python=sys.version.split()[0], numpy=np.__version__),
    )


def memory_windows_of(windows):
    first, last = windows['main']
    return {'main': (first, last), 'second half': windows['second half']}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('run', help='run:first-last, e.g. A28R_301k_design_l15_snow_s924_turned_refit:10-49')
    parser.add_argument('--compare', help='a further window first-last screened alike, e.g. 10-28')
    parser.add_argument('--surrogates', type=int, default=10000)
    parser.add_argument('--half-surrogates', type=int, help='surrogates for the halves (default: --surrogates)')
    parser.add_argument('--seed', type=int, default=20261108)
    parser.add_argument('--work', required=True, help='folder for the series and spectra (on the research drive)')
    parser.add_argument('--runs', help='folder holding the run folders (default: the runner\'s)')
    parser.add_argument('--out', help='product path (default: results/gcm/cycles_<run>_y<first>-<last>.json)')
    parser.add_argument('--reread', action='store_true', help='read the year files even if the series are cached')
    parser.add_argument('--memory-surrogates', type=int, default=5000)
    parser.add_argument('--memory-half-surrogates', type=int, default=2000)
    parser.add_argument('--memory-only', action='store_true', help='compute and keep only the two-timescale screens')
    parser.add_argument('--from-screens', action='store_true',
                        help='remake the summary from the surrogate results kept in --work, without recomputing them')
    run(parser.parse_args(argv))
    return 0


if __name__ == '__main__':
    sys.exit(main())
