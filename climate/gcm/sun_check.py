"""Where a GCM run's Sun shines: its sunlight by longitude and by latitude, against PlaSim's solar-position code and
against the orbit the run asked for.

    climate/gcm/.venv/bin/python -m climate.gcm.sun_check A28_dim5 --years 15 24

writes ../results/gcm/sun_check_<run>.json in seconds, from the run's output, namelists and log.

Averaged over the years, the sunlight arriving along a row of latitude should be the same at every longitude, and
the rows' means should follow the Moon's 1.54 degree tilt. PlaSim's `solang` (plasim/src/radmod.f90) turns the Sun
at the planet's rotation rate, ROTSPD turns per 1440 minutes, and restarts its clock every model day of NTSPD steps,
which PlaSim makes a solar day: `solang_mean_mu` gives the map that makes. PlaSim uses the tilt it is given only with
a fixed orbit (NFIXORB=1). Otherwise it computes Earth's orbit for its start year from Berger's (1978) series and
reports that orbit in its log.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-sun-check/1'
EVIDENCE = ('Time means of the run\'s own output over the years (rst - rsut, the sunlight arriving at the top of the '
            'air, as means over each output interval), with its namelists and the orbit its log reports. The solang '
            'map is PlaSim\'s arithmetic for the Sun\'s hour angle (NGENKEPLERIAN=0, no dawn threshold) with the Sun '
            'on the equator; the annual means integrate the daily-mean sunlight over the orbit, uniformly in true '
            'longitude.')
READING_RULE = ('along_row.relative is each longitude\'s sunlight over its row\'s mean, 1 everywhere for a Sun that '
                'sweeps evenly; solang_relative is the same from PlaSim\'s solar-position arithmetic alone. by_latitude '
                'compares the rows\' means with the annual means for the orbit PlaSim\'s log reports and for the one '
                'the namelist asked for. A run with a sound Sun has relative within 0.01 of 1 and rows that follow '
                'the orbit asked for.')


def plasim_day_steps(rotspd: float, days_per_year: int, mpstep_min: float) -> int:
    """Steps in PlaSim's model day (NTSPD), the solar day it derives in plasim.f90 from the rotation and the
    rotations it counts in a year, N_DAYS_PER_YEAR."""
    if rotspd == 1.0:
        solar = 86400.0
    else:
        sidereal = 86400.0 / rotspd
        solar = sidereal * days_per_year / (days_per_year - 1) if days_per_year != 1 else sidereal
    return int(round(solar)) // int(round(mpstep_min * 60))


def solang_mean_mu(lat_deg: float, nlon: int, day_steps: int, mpstep_min: float, rotspd: float,
                   declination_deg: float = 0.0) -> np.ndarray:
    """Mean cosine of the Sun's zenith angle by longitude over one model day, as solang has it: the clock counts
    whole minutes from the start of the model day and turns the Sun ROTSPD times per 1440 of them."""
    minutes = np.floor(np.arange(day_steps) * mpstep_min)
    hour = minutes[:, None] * rotspd * 2 * np.pi / 1440.0 + np.arange(nlon)[None, :] * 2 * np.pi / nlon - np.pi
    lat, dec = np.radians(lat_deg), np.radians(declination_deg)
    return np.maximum(0.0, np.sin(dec) * np.sin(lat) + np.cos(lat) * np.cos(dec) * np.cos(hour)).mean(axis=0)


def annual_mean_insolation(lat_deg, obliquity_deg: float, eccentricity: float, sunlight_w_m2: float) -> np.ndarray:
    """Annual-mean daily sunlight at the top of the air by latitude. Time weighted by the inverse square of the
    distance is uniform in true longitude, so the mean over the orbit is a plain mean over it."""
    along = (np.arange(3600) + 0.5) * 2 * np.pi / 3600
    dec = np.arcsin(np.sin(np.radians(obliquity_deg)) * np.sin(along))[:, None]
    phi = np.radians(np.atleast_1d(np.asarray(lat_deg, float)))[None, :]
    h0 = np.arccos(np.clip(-np.tan(phi) * np.tan(dec), -1.0, 1.0))
    mu = (h0 * np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.sin(h0)) / np.pi
    return sunlight_w_m2 / np.sqrt(1.0 - eccentricity ** 2) * mu.mean(axis=0)


def namelist(path: Path) -> dict:
    """A PlaSim namelist's numeric KEY = value pairs."""
    out = {}
    for line in path.read_text().splitlines():
        m = re.match(r'\s*(\w+)\s*=\s*([-+.\w]+)', line)
        if m:
            try:
                out[m.group(1).upper()] = float(m.group(2).lower().replace('d', 'e'))
            except ValueError:
                pass
    return out


def logged_orbit(log: Path) -> dict:
    """The eccentricity and obliquity PlaSim reports using, from a year's diagnostic log."""
    text = log.read_text(errors='replace')
    grab = lambda label: float(re.search(label + r'\s*=\s*([-+.\dEe]+)', text).group(1))
    return dict(obliquity_deg=grab(r'Obliquity \(deg\)'), eccentricity=grab(r'\* Eccentricity'))


def check(run: str, years) -> dict:
    import netCDF4
    folder = HERE / 'runs' / run / 'model'
    nl = {**namelist(folder / 'plasim_namelist'), **namelist(folder / 'planet_namelist')}
    arriving, count = 0.0, 0
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(folder / f'MOST.{year:05d}.nc') as d:
            arriving = arriving + (np.asarray(d['rst'][:], float) - np.asarray(d['rsut'][:], float)).sum(axis=0)
            count += d['time'].size
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
    arriving = arriving / count
    steps = plasim_day_steps(nl['ROTSPD'], int(nl['N_DAYS_PER_YEAR']), nl['MPSTEP'])
    sweep = steps * nl['MPSTEP'] * nl['ROTSPD'] * 360.0 / 1440.0
    j = int(np.argmin(np.abs(lat)))
    row = arriving[j] / arriving[j].mean()
    model = solang_mean_mu(lat[j], lon.size, steps, nl['MPSTEP'], nl['ROTSPD'])
    model = model / model.mean()
    used = logged_orbit(folder / f'MOST_DIAG.{years[0]:05d}')
    asked = dict(obliquity_deg=nl['OBLIQ'], eccentricity=nl['ECCEN'])
    zonal = arriving.mean(axis=1)
    expect_used = annual_mean_insolation(lat, used['obliquity_deg'], used['eccentricity'], nl['GSOL0'])
    expect_asked = annual_mean_insolation(lat, asked['obliquity_deg'], asked['eccentricity'], nl['GSOL0'])
    rms = lambda a: float(np.sqrt(np.mean((zonal - a) ** 2)))
    weight = np.cos(np.radians(lat))
    return dict(schema=SCHEMA, run=run, years=list(years), evidence=EVIDENCE, reading_rule=READING_RULE,
                model_day=dict(steps=steps, days=steps * nl['MPSTEP'] / 1440.0, sun_sweep_deg=sweep,
                               skip_back_deg=sweep - 360.0),
                along_row=dict(lat_deg=float(lat[j]), lon_deg=lon.tolist(), arriving_w_m2=arriving[j].tolist(),
                               relative=row.tolist(), solang_relative=model.tolist(),
                               relative_range=[float(row.min()), float(row.max())],
                               correlation_with_solang=float(np.corrcoef(row, model)[0, 1])),
                orbit=dict(asked=asked, used=used),
                by_latitude=dict(lat_deg=lat.tolist(), arriving_w_m2=zonal.tolist(), orbit_used_w_m2=expect_used.tolist(),
                                 orbit_asked_w_m2=expect_asked.tolist(), rms_against_used_w_m2=rms(expect_used),
                                 rms_against_asked_w_m2=rms(expect_asked)),
                global_w_m2=float(np.sum(zonal * weight) / weight.sum()))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('run')
    parser.add_argument('--years', type=int, nargs=2, default=(15, 24))
    args = parser.parse_args(argv)
    result = check(args.run, tuple(args.years))
    result['producer'] = dict(domain='climate', files={'gcm/sun_check.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
    from climate.crm.ring_crossings import rounded
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f'sun_check_{args.run}.json'
    path.write_text(json.dumps(rounded(result, 4), indent=1) + '\n')
    a, b = result['along_row'], result['by_latitude']
    print(f"model day {result['model_day']['steps']} steps ({result['model_day']['days']:.2f} d): the Sun sweeps "
          f"{result['model_day']['sun_sweep_deg']:.1f} deg and skips back {result['model_day']['skip_back_deg']:.1f}")
    print(f"row at {a['lat_deg']:.1f}: sunlight {a['relative_range'][0]:.3f}-{a['relative_range'][1]:.3f} of its mean; "
          f"correlation with solang {a['correlation_with_solang']:.4f}")
    print(f"tilt asked {result['orbit']['asked']['obliquity_deg']:.3f}, used {result['orbit']['used']['obliquity_deg']:.3f}; "
          f"rows' means rms {b['rms_against_used_w_m2']:.1f} W/m2 from the orbit used, {b['rms_against_asked_w_m2']:.1f} from the one asked")
    print(path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
