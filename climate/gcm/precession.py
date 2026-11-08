"""The Moon's spin axis precessing with its orbit's node, in PlaSim (8 November 2026).

The Moon stays in Cassini state 2: its spin axis precesses with the node of its orbit, so the equinox, where the lunar
equator crosses the ecliptic, regresses once a nodal period while the perihelion of the orbit the Moon shares with Earth
around the Sun stays put. The Sun's tilt seasons then follow the eclipse year, 346.6 days, and drift through the
distance's year once a nodal period, so each pole has the summers nearer the Sun for half of every 18.6 years. PlaSim
as the runs had it holds the equinox fixed: the turned design run's southern summers fall at perihelion every year.

PlaSim's orbit (radmod.f90, orb_decl, from CCSM) places the Sun by the mean longitude at the equinox, lambm0, and the
longitude of perihelion from the equinox plus pi, mvelpp. EDITS add to both the angle the equinox has regressed by since
step PRECSTEP0, a whole turn in PRECYEARS orbital years. The mean anomaly, lambm - mvelpp, and with it the Sun's
distance keep the orbit's year, while the Sun's longitude from the equinox, and so its declination, gains a turn each
nodal period. The general Keplerian orbit (gen_orb_decl) takes the angle on its perihelion's longitude. PRECYEARS 0, the
default, adds nothing, so a run without it computes what it did before.

The nodal period, 6,798.38 days, is geology's source (interior_sources.json, nodal_period), as the orientation study
and the polar light study take it. PlaSim's orbit turns once in n_steps_per_year steps, the 365.256-day sidereal year
in its 30-minute steps (17,532; plasim.f90), whatever the length of the run's output years (12 lunar days, 17,160
steps), so PRECYEARS is the nodal period in those orbital years, 18.61, and the declination repeats every eclipse year,
346.6 days. A run counts the regression from its own first step, the step its initial restart holds, so a branch from a
run without precession continues that run's orbit unbroken.

    <climate venv>/bin/python -m climate.gcm.precession check FOLDER YEAR [--out JSON]
        # a run's zonal-mean cosine of the Sun's zenith angle in one year against the orbit's, with and without the
        # run's precession
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import struct
import sys

import numpy as np

from climate.gcm import levels as lv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCES = ROOT / 'geology' / 'interior_sources.json'
RUNS = HERE / 'runs'
VERNAL_EQUINOX = 80.5 / 365.0          # orb_decl's place of the equinox in the orbit's year

MARKER = 'Terluna: the equinox regresses with the node'
EDITS = [
    ("      real :: rasc=0.  ! Solar right ascension in radians\n",
     "      real :: rasc=0.  ! Solar right ascension in radians\n"
     "      real :: precyears = 0.0 ! Terluna: the equinox regresses with the node, once in precyears years\n"
     "      real :: precstep0 = 0.0 ! Terluna: the step from which the regression is counted\n"
     "      real :: precang = 0.0   ! Terluna: the angle the equinox has regressed by, radians\n"),
    ("      namelist/radmod_nl/ndcycle,ncstsol,solclat,solcdec,no3,co2        &\n",
     "      namelist/radmod_nl/precyears,precstep0 ! Terluna: the equinox's regression\n"
     "      namelist/radmod_nl/ndcycle,ncstsol,solclat,solcdec,no3,co2        &\n"),
    ("      call mpbcr(rayscale)\n",
     "      call mpbcr(rayscale)\n      call mpbcr(precyears)\n      call mpbcr(precstep0)\n"),
    ("      if (ngenkeplerian == 0) then\n"
     "          call orb_decl(zcday, eccen, mvelpp, lambm0, obliqr, orbnu, lambm, rasc, zdecl, eccf)\n",
     "      precang = 0.0 ! Terluna: the angle the equinox has regressed by since step precstep0\n"
     "      if (precyears > 0.0) precang = TWOPI * (real(nstep) - precstep0) / (precyears * real(n_steps_per_year))\n"
     "      if (ngenkeplerian == 0) then\n"
     "          call orb_decl(zcday, eccen, mvelpp + precang, lambm0 + precang, obliqr, orbnu, lambm, rasc, "
     "zdecl, eccf)\n"),
    ("          call gen_orb_decl(zcday, eccen, obliqr, mvelpp, orbnu, lambm, rasc, zdecl, eccf)\n",
     "          call gen_orb_decl(zcday, eccen, obliqr, mvelpp + precang, orbnu, lambm, rasc, zdecl, eccf)\n"),
]
PATCHES = [('radmod.f90', MARKER, EDITS)]


# ==== the setting =====================================================================================================
def node_days(path=SOURCES) -> float:
    """The regression period of the lunar orbit's node, days, from geology's sources."""
    sources = json.loads(Path(path).read_text())['sources']
    return float(next(s['values_used']['node_days'] for s in sources if s['key'] == 'nodal_period'))


def setting(year_days: float) -> dict:
    """The precession a run records: the nodal period, the orbit's year it is counted in and the period in those years.
    The step it counts from, step0, joins it once the run's initial state is known (exoplasim_run.run)."""
    node = node_days()
    return dict(node_days=node, year_days=float(year_days), period_years=node / float(year_days),
                source=f'{SOURCES.relative_to(ROOT)} (nodal_period)')


def restart_step(path) -> int:
    """The step a PlaSim restart stands at, its nstep record."""
    for name, data in lv.read_records(path):
        if name.strip() == 'nstep':
            return int(struct.unpack('<q' if len(data) == 8 else '<i', data)[0])
    raise ValueError(f'{path} holds no nstep record')


def namelist(entry: dict) -> list:
    """radmod_namelist's entries for a run's precession entry, as (name, text)."""
    return [('PRECYEARS', repr(float(entry['period_years']))), ('PRECSTEP0', f"{int(entry['step0'])}.0")]


# ==== PlaSim's orbit ==================================================================================================
def orbit_constants(eccen, mvelp_deg):
    """orb_params's mvelpp and lambm0 (radians) for a fixed orbit of eccentricity eccen and longitude of perihelion
    mvelp_deg."""
    mvelpp = np.radians(mvelp_deg + 180.0)
    beta = np.sqrt(1.0 - eccen ** 2)
    lambm0 = 2.0 * ((0.5 * eccen + 0.125 * eccen ** 3) * (1.0 + beta) * np.sin(mvelpp)
                    - 0.25 * eccen ** 2 * (0.5 + beta) * np.sin(2.0 * mvelpp)
                    + 0.125 * eccen ** 3 * (1.0 / 3.0 + beta) * np.sin(3.0 * mvelpp))
    return float(mvelpp), float(lambm0)


def orb_decl(calday, eccen, mvelpp, lambm0, obliqr):
    """PlaSim's orb_decl: the Sun's declination (rad) and distance factor (1 AU / r)^2 at the fraction calday of the
    orbit's year."""
    lambm = lambm0 + (np.asarray(calday, float) - VERNAL_EQUINOX) * 2.0 * np.pi
    lmm = lambm - mvelpp
    sinl = np.sin(lmm)
    lamb = lambm + eccen * (2.0 * sinl + eccen * (1.25 * np.sin(2.0 * lmm)
                                                  + eccen * ((13.0 / 12.0) * np.sin(3.0 * lmm) - 0.25 * sinl)))
    invrho = (1.0 + eccen * np.cos(lamb - mvelpp)) / (1.0 - eccen ** 2)
    return np.arcsin(np.sin(obliqr) * np.sin(lamb)), invrho ** 2


def angle(nstep, step0, period_years, steps_per_year):
    """The patch's precang: the angle (rad) the equinox has regressed by at step nstep, a turn in period_years orbital
    years of steps_per_year steps."""
    nstep = np.asarray(nstep, float)
    if period_years <= 0.0:
        return np.zeros_like(nstep)
    return 2.0 * np.pi * (nstep - step0) / (period_years * steps_per_year)


def sun(nstep, steps_per_year, eccen, mvelp_deg, obliq_deg, period_years=0.0, step0=0):
    """The Sun's declination (rad) and distance factor at step nstep as solang finds them, with the patch's precession;
    steps_per_year is PlaSim's n_steps_per_year."""
    nstep = np.asarray(nstep)
    mvelpp, lambm0 = orbit_constants(eccen, mvelp_deg)
    a = angle(nstep, step0, period_years, steps_per_year)
    return orb_decl(np.mod(nstep, steps_per_year) / steps_per_year, eccen, mvelpp + a, lambm0 + a,
                    np.radians(obliq_deg))


def orbit_steps(model_dir: Path) -> int:
    """PlaSim's n_steps_per_year for a run: the sidereal year in its steps, nint(SIDEREAL_YEAR / (MPSTEP * 60))."""
    year = namelist_value(model_dir / 'planet_namelist', 'SIDEREAL_YEAR')
    step = namelist_value(model_dir / 'plasim_namelist', 'MPSTEP')
    return int(np.rint(year / (step * 60.0)))


def mean_cos_zenith(lat_deg, delta):
    """The mean over every hour angle of the cosine of the Sun's zenith angle while the Sun is up, 0 below: the zonal
    mean of PlaSim's gmu0 at one step (dawn 0), latitudes (deg) by declinations (rad)."""
    phi = np.radians(np.asarray(lat_deg, float))[:, None]
    d = np.asarray(delta, float)[None, :]
    h = np.arccos(np.clip(-np.tan(phi) * np.tan(d), -1.0, 1.0))
    return (h * np.sin(phi) * np.sin(d) + np.cos(phi) * np.cos(d) * np.sin(h)) / np.pi


# ==== a run against the orbit =========================================================================================
def namelist_value(path: Path, key: str) -> float:
    """A number from a PlaSim namelist file."""
    m = re.search(rf'^\s*{key}\s*=\s*([-+0-9.eEdD]+)', Path(path).read_text(), re.M)
    if m is None:
        raise ValueError(f'{path} sets no {key}')
    return float(m.group(1).replace('d', 'e').replace('D', 'e'))


def planet_namelist(model_dir: Path) -> dict:
    """The orbit PlaSim read: ECCEN, MVELP and OBLIQ from the run's planet_namelist."""
    return {key: namelist_value(model_dir / 'planet_namelist', key) for key in ('ECCEN', 'MVELP', 'OBLIQ')}


def check(folder: str, year: int) -> dict:
    """A run's czen (PlaSim's code 318, the mean of gmu0 over each output interval) in one year, as zonal means, against
    the orbit's with the run's precession and with the equinox held: the root mean square and the largest difference
    over every latitude and interval. The year starts at the step its previous year's restart holds (the initial state's
    for a first year); an interval of NSTPW steps takes the steps from its start."""
    import netCDF4
    rundir = RUNS / folder
    model = rundir / 'model'
    progress = json.loads((rundir / 'progress.json').read_text())
    cfg = progress['configuration']
    m = cfg['model']
    per_year = orbit_steps(model)
    per_record = m['steps_per_lunar_day'] // m['writes_per_lunar_day']
    if year > 0:
        start = restart_step(model / f'MOST_REST.{year - 1:05d}')
    elif isinstance(cfg.get('initial_state'), dict):
        start = restart_step(cfg['initial_state']['path'])
    else:
        start = 0
    orbit = planet_namelist(model)
    prec = cfg.get('precession') or {}
    period, step0 = prec.get('period_years', 0.0), prec.get('step0', start)
    with netCDF4.Dataset(model / f'MOST.{year:05d}.nc') as d:
        lat = np.asarray(d.variables['lat'][:], float)
        czen = np.asarray(d.variables['czen'][:], float).mean(axis=-1)             # (records, lat): zonal means
    steps = start + np.arange(czen.shape[0] * per_record)
    out = {}
    for case, years in (('with_run_precession', period), ('equinox_held', 0.0)):
        delta, _ = sun(steps, per_year, orbit['ECCEN'], orbit['MVELP'], orbit['OBLIQ'], years, step0)
        expected = mean_cos_zenith(lat, delta).reshape(lat.size, czen.shape[0], per_record).mean(axis=2).T
        diff = czen - expected
        out[case] = dict(rms=float(np.sqrt(np.mean(diff ** 2))), max_abs=float(np.abs(diff).max()),
                         max_abs_poleward_of_60=float(np.abs(diff[:, np.abs(lat) >= 60.0]).max()))
    delta, _ = sun(steps, per_year, orbit['ECCEN'], orbit['MVELP'], orbit['OBLIQ'], period, step0)
    return dict(run=folder, year=year, start_step=int(start), orbit_steps_per_year=int(per_year),
                steps_per_record=int(per_record), orbit=orbit, precession=prec or None,
                precession_angle_deg=[float(np.degrees(angle(s, step0, period, per_year)))
                                      for s in (steps[0], steps[-1])],
                declination_deg_range=[float(np.degrees(delta.min())), float(np.degrees(delta.max()))], **out)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest='command', required=True)
    c = sub.add_parser('check', help="a run's sun in one year against the orbit")
    c.add_argument('folder')
    c.add_argument('year', type=int)
    c.add_argument('--out', type=Path, help='write the result as JSON there')
    args = parser.parse_args(argv)
    result = check(args.folder, args.year)
    text = json.dumps(result, indent=1)
    if args.out:
        args.out.write_text(text + '\n')
    print(text)
    return 0


if __name__ == '__main__':
    sys.exit(main())
