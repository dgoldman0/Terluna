"""Run the Open Moon's ExoPlaSim experiments one simulated year at a time, resumable after any stop.

    climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --years 40          # run or resume A
    climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --status            # progress so far
    climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --stop              # stop it cleanly
    climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --years 15 --folder A_full_column_clouds \
        --set cloud_water=full_column --start-from climate/gcm/runs/A/model/MOST_REST.00039   # a branch

Run from the repository root with the GCM environment (climate/gcm/.venv, which holds ExoPlaSim and
MPI-built executables). Each experiment lives in climate/gcm/runs/<name>/ (ignored by Git; on the
author's machine climate/gcm/runs links to the external drive, /media/projectspace/terluna-research/
gcm-runs). Every year's output and restart is kept, about 80 MB a year; --prune keeps only the last
five years' output and three restarts, plus every tenth year, for small disks:

- model/        ExoPlaSim's working folder: executable, namelists, restart files MOST_REST.<year>,
                output MOST.<year>.nc (ten means per lunar day);
- progress.json the completed years with a summary of each, written atomically after every year;
- run.lock      the process (and process group) that owns the run, so two runners never share a
                folder.

Resuming: ExoPlaSim writes each year's end state to MOST_REST.<year> and starts the next year from
it. After a stop, the runner rebuilds the same model, hands it the last completed year's restart
and continues, so at most the year in progress is lost. Anything the interrupted year left behind
is removed first, including model processes still working in the folder after their runner was
killed. The runner leads its own process group, and a stop signal (--stop, or SIGTERM, SIGINT or
SIGHUP to the runner) takes ExoPlaSim's MPI ranks down with it. A run refuses to resume if anything
that determines the model changed: its settings, the number of MPI ranks, ExoPlaSim's version or
the input files the runner writes (land mask, topography, spectrum). Tested by killing a run with
SIGKILL a minute into its second year: the resumed year's restart and output were byte-identical
to an uninterrupted run's.

The planet (radius, gravity, rotation, orbit), the land-sea mask and topography come from
climate/gcm/products (python -m climate.gcm.boundary). The mask is binary per cell, as PlaSim
needs: cells are made sea in order of least land fraction until the sea covers the experiment's
water share. The shield enters as sunlight: the solar spectrum (TSIS-1 HSRS) times the shield's
transmission (protection/spectra/shield_transmission.json) sets both the total flux and the split
between ExoPlaSim's two solar bands. The air is the 1-D models' dry composition (O2 at 21,227 Pa,
Earth's Ar/N2, 400 ppm CO2, N2 the balance), with the pressure set at the water level as on Earth:
land stands above it, so the global mean surface pressure is a few per cent lower. The land albedo
is a uniform placeholder. PlaSim fits 12 lunar days into its year, which makes its solar day 29.8
Earth days (0.9% longer than the real 29.53; the rotation keeps the true 27.32 days) and its year
357.5 days. Output is ten means per lunar day. The random seed is fixed, so runs repeat exactly.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import sys
import time
import numpy as np
from shared.constants import STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = HERE / 'runs'
PRODUCTS = HERE / 'products'
SHIELDS = ROOT / 'protection' / 'spectra' / 'shield_transmission.json'
EARTH_RADIUS_M = 6371000.0          # ExoPlaSim takes the radius in units of this
O2_PA = 21227.0                     # Earth's O2 partial pressure, as in the 1-D models
AR_PER_N2 = 0.00934 / 0.78084       # Earth's ratio
CO2_PPM = 400.0
MIN_FREE_GB = 3.0                   # pause rather than fill the disk
KEEP_RESTARTS = 3                   # with --prune: most recent restarts kept, and every tenth year's
KEEP_OUTPUT_YEARS = 5               # with --prune: most recent outputs kept, and every tenth year's

EXPERIMENTS = {
    'A': dict(pressure_pa=121590.0, water=25, shield='titania_stack', mldepth=50.0, purpose='design case'),
    'B': dict(pressure_pa=121590.0, water=35, shield='titania_stack', mldepth=50.0, purpose='water sensitivity'),
    'C': dict(pressure_pa=101325.0, water=25, shield='titania_stack', mldepth=50.0, purpose='pressure sensitivity'),
    # The scenario's 28% seas (shared/scenarios/water.json) with the geography domain's rain-fed lakes as
    # water cells at their own levels: about 40% of the surface under water. A restart carries PlaSim's
    # land-sea mask and orography, so a change of coastline needs a cold start.
    'A28': dict(pressure_pa=121590.0, water=28, lakes=True, shield='titania_stack', mldepth=50.0,
                purpose='design case with the atlas seas and rain-fed lakes'),
}
# PlaSim fits a whole number of solar days into its year: 12 lunar days of 1430 half-hour steps
# (29.8 Earth days each, 0.9% longer than the real 29.53, with the true 27.32-day rotation), 357.5
# Earth days in all. Output: 10 means per lunar day (about 3 days each), 120 a year.
MODEL = dict(resolution='T21', layers=10, timestep_min=30.0, land_albedo=0.2, steps_per_lunar_day=1430,
             writes_per_lunar_day=10, lunar_days_per_year=12, seed=1,
             # Radiation calibrated against the line-by-line model on cloud-free lunar columns at 280 and
             # 310 K (see the README): a multiplier on ExoPlaSim's spectrum-derived Rayleigh coefficient,
             # and PlaSim's water-vapour continuum coefficient. 1.0 and 0.024 reproduce the unmodified
             # model, which absorbed about 14 W/m2 too much sunlight and emitted 12 W/m2 too little when warm.
             rayleigh_scale=1.8, h2o_continuum=0.004,
             # How the radiation sees cloud water (CLOUD_WATER below). Runs A and B used 'plasim', ExoPlaSim
             # unmodified; 'earth_path' corrects its Earth-fitted formula for the Moon's gravity.
             cloud_water='earth_path',
             # Seconds per day in PlaSim's convective cloud cover, 0.245 + 0.125 ln(convective rain in mm per
             # day), capped at 0.8 (Slingo's fit to Earth). PlaSim divides the rain by its own solar day, 29.8
             # Earth days on the Moon, so it reads the rain about 30 times too heavy and raises the cover by
             # 0.42, up to the cap. 86400 counts it per Earth day, as fitted; 0 keeps PlaSim's behaviour (runs
             # A and B and the cloud-water bracket).
             convective_day_s=86400.0,
             # PlaSim's clear-sky diagnostic (1 on, 0 off): the radiation computed a second time without
             # clouds, for the clouds' radiative effect. It changes no prognostic field.
             clear_sky=1,
             # A multiplier on the sunlight the shield passes, the same at every wavelength: 0.92 is a shield
             # that passes 8% less. 1.0 keeps the shield product's transmission.
             sunlight_scale=1.0)
# Cloud water. PlaSim uses it only for the clouds' radiation: rain comes from its condensation and
# convection schemes, which rain out any excess vapour at once. It diagnoses the water with CCM3's
# formula, 0.21 g/m3 at the ground falling off with geometric height on a scale of
# 700 ln(1 + precipitable water in kg/m2) m, fitted to Earth. At a sixth of Earth's gravity a pressure
# level stands six times higher, so on the Moon the formula leaves clouds above about 20 km, most of
# the model's cloud cover, almost without water. The Terluna patch to rainmod.f90 allows two corrections:
#   plasim       the formula as PlaSim has it;
#   earth_path   heights and precipitable water counted as on Earth (the formula in pressure terms):
#                each cloud holds the water path it would hold at the same pressures on Earth;
#   full_column  the same mixing ratio kept over the Moon's column, which holds six times Earth's mass
#                per pascal of pressure: six times Earth's water path, the upper bound.
CLOUD_WATER = ('plasim', 'earth_path', 'full_column')
# Kept in each year's 3-day means (ExoPlaSim writes about 100 fields by default, 100 MB a year).
OUTPUT = ['ta', 'ua', 'va', 'hus', 'cl',
          'ts', 'tas', 'maxt', 'mint', 'ps', 'psl', 'pr', 'prc', 'evap', 'prw', 'clt', 'sic', 'sit', 'snd', 'mrso',
          'rst', 'rsut', 'rlut', 'ntr', 'nbr', 'rss', 'rls', 'hfss', 'hfls', 'alb', 'czen', 'lsm']
# PlaSim's clear-sky fluxes, written when the diagnostic is on but unknown to ExoPlaSim's postprocessor,
# which the runner teaches them. Unlike the other fluxes they are snapshots at each output time, not
# means over it: sound for global and lunar-day averages, not for single times.
CLEAR_SKY = {101: ('rsscs', 'surface_net_shortwave_flux_clear_sky'),
             102: ('rlscs', 'surface_net_longwave_flux_clear_sky'),
             103: ('rstcs', 'toa_net_shortwave_flux_clear_sky'),
             104: ('rltcs', 'toa_net_longwave_flux_clear_sky')}


def output_variables(model: dict) -> list:
    return OUTPUT + ([name for name, _ in CLEAR_SKY.values()] if model['clear_sky'] else [])


def cloud_water_namelist(choice: str, gravity: float) -> tuple:
    """CWGREF and CWSCALE for the patched cloud scheme: the gravity at which the cloud-water formula
    counts heights and precipitable water (0: the planet's own), and a multiplier on the water."""
    if choice == 'plasim':
        return 0.0, 1.0
    if choice == 'earth_path':
        return STANDARD_GRAVITY, gravity / STANDARD_GRAVITY
    if choice == 'full_column':
        return STANDARD_GRAVITY, 1.0
    raise ValueError(f'cloud_water must be one of {CLOUD_WATER}, not {choice!r}')


# The Terluna patches to ExoPlaSim's source, applied once, idempotently. Radiation (plasim/src/radmod.f90):
# a namelist multiplier, RAYSCALE, on the Rayleigh coefficient the model derives from the stellar spectrum.
RADMOD_EDITS = [
    ("      real :: rcoeff = 1.0       ! Rayleigh scattering coefficient for cross section dependence\n",
     "      real :: rcoeff = 1.0       ! Rayleigh scattering coefficient for cross section dependence\n"
     "      real :: rayscale = 1.0     ! Terluna: multiplier on rcoeff, calibrated against line-by-line\n"),
    ("        rcoeff = (zcross1 + zcross2) * zsolar1 / z1 / zchi !Using default zsolar=0.517 here\n",
     "        rcoeff = (zcross1 + zcross2) * zsolar1 / z1 / zchi !Using default zsolar=0.517 here\n"
     "        rcoeff = rcoeff * rayscale ! Terluna: calibration multiplier\n"),
    ("     &               ,nsimplealbedo,nstarfile,starfile,starfilehr,minwavel\n",
     "     &               ,nsimplealbedo,nstarfile,starfile,starfilehr,minwavel,rayscale\n"),
    ("      call mpbcr(minwavel)\n",
     "      call mpbcr(minwavel)\n      call mpbcr(rayscale)\n"),
]
# Clouds (plasim/src/rainmod.f90, mkclouds): CWGREF, the gravity at which the cloud-water formula counts
# heights and precipitable water (0, the default, keeps the planet's own), and CWSCALE, a multiplier on
# the diagnosed water. The defaults leave the model as it was.
RAINMOD_EDITS = [
    ("      real :: clwfac    = 0. ! smothing for cloud suppression\n",
     "      real :: clwfac    = 0. ! smothing for cloud suppression\n"
     "      real :: cwgref    = 0. ! Terluna: gravity for the cloud-water heights (0: the planet's)\n"
     "      real :: cwscale   = 1. ! Terluna: multiplier on the diagnosed cloud water\n"),
    ("rbeta,rcritmod,rcritslope", "rbeta,rcritmod,rcritslope,cwgref,cwscale"),
    ("      call mpbcr(gamma)\n", "      call mpbcr(gamma)\n      call mpbcr(cwgref)\n      call mpbcr(cwscale)\n"),
    ("      real zzf(NHOR,NLEV),zzh(NHOR),zdh(NHOR)\n",
     "      real zzf(NHOR,NLEV),zzh(NHOR),zdh(NHOR)\n      real zgc ! Terluna: gravity for the heights below\n"),
    ("      dcc(:,1:NLEV)=0.\n", "      dcc(:,1:NLEV)=0.\n      zgc=ga\n      if(cwgref > 0.) zgc=cwgref\n"),
    ("gascon/ga*ALOG(sigmah(jlev-1)/sigmah(jlev))", "gascon/zgc*ALOG(sigmah(jlev-1)/sigmah(jlev))"),
    ("gascon/ga*ALOG(sigma(1)/sigmah(1))*0.5", "gascon/zgc*ALOG(sigma(1)/sigmah(1))*0.5"),
    ("      zzh(:)=700.*ALOG(1.+dqvi(:))\n", "      zzh(:)=700.*ALOG(1.+dqvi(:)*(ga/zgc))\n"),
    ("        dql(:,jlev)=MAX(dql(:,jlev),1.E-9)\n", "        dql(:,jlev)=MAX(cwscale*dql(:,jlev),1.E-9)\n"),
]
# Convective cloud cover (rainmod.f90, mkclouds): CONVDAY, the seconds per day in which the cover formula
# counts rain (0, the default, keeps the model's solar day).
CONVECTIVE_EDITS = [
    ("      real :: cwscale   = 1. ! Terluna: multiplier on the diagnosed cloud water\n",
     "      real :: cwscale   = 1. ! Terluna: multiplier on the diagnosed cloud water\n"
     "      real :: convday   = 0. ! Terluna: seconds per day in the convective cloud formula (0: solar day)\n"),
    ("rbeta,rcritmod,rcritslope,cwgref,cwscale", "rbeta,rcritmod,rcritslope,cwgref,cwscale,convday"),
    ("      call mpbcr(cwscale)\n", "      call mpbcr(cwscale)\n      call mpbcr(convday)\n"),
    ("      zrfac = solar_day * 1000.0 ! convert m/s into mm/day\n",
     "      zrfac = solar_day * 1000.0 ! convert m/s into mm/day\n"
     "      if(convday > 0.) zrfac = convday * 1000.0 ! Terluna: per convday seconds\n"),
]
# (file, marker, edits) in the order applied; a marker appears only in the text its patch adds.
PATCHES = [('radmod.f90', 'Terluna: calibration multiplier', RADMOD_EDITS),
           ('rainmod.f90', 'Terluna: gravity for the heights below', RAINMOD_EDITS),
           ('rainmod.f90', 'Terluna: seconds per day in the convective cloud formula', CONVECTIVE_EDITS)]


def ensure_patched() -> dict:
    """Apply the Terluna patches to ExoPlaSim's source where missing (removing the compiled executables
    so ExoPlaSim rebuilds them), and return each patched file's hash for the run key."""
    import exoplasim
    src = Path(exoplasim.__file__).parent / 'plasim' / 'src'
    hashes = {}
    for name, marker, edits in PATCHES:
        path = src / name
        text = path.read_text()
        if marker not in text:
            for old, new in edits:
                if text.count(old) != 1:
                    raise RuntimeError(f'ExoPlaSim {name} is not the version the Terluna patch expects near {old.strip()!r}')
                text = text.replace(old, new)
            path.write_text(text)
            for exe in (src.parent / 'run').glob('most_plasim_*.x'):
                exe.unlink()
        hashes[path.stem] = hashlib.sha256(path.read_bytes()).hexdigest()[:16]
    return hashes


def _atomic_write(path: Path, text: str):
    tmp = path.with_name(path.name + '.part')
    tmp.write_text(text)
    tmp.replace(path)


def configuration(name, ncpus, overrides=None):
    """Everything that determines the run apart from the input files the runner writes."""
    from importlib.metadata import version
    exp = EXPERIMENTS[name]
    planet = json.loads((PRODUCTS / 'moon_gcm_configuration.json').read_text())['planet']
    dry_bar = exp['pressure_pa'] / 1e5
    co2 = CO2_PPM * 1e-6 * dry_bar
    n2_ar = dry_bar - O2_PA / 1e5 - co2
    model = {**MODEL, **(overrides or {})}
    cloud_water_namelist(model['cloud_water'], planet['gravity_m_s2'])        # refuses an unknown choice
    if model['convective_day_s'] < 0:
        raise ValueError('convective_day_s must be 0 (the solar day) or a positive number of seconds')
    if not 0.0 < model['sunlight_scale'] <= 1.5:
        raise ValueError('sunlight_scale must lie between 0 and 1.5')
    return dict(experiment=name, **exp, model=model, output=output_variables(model), ncpus=ncpus,
                exoplasim=version('exoplasim'), source=ensure_patched(),
                planet=dict(radius=planet['radius_m'] / EARTH_RADIUS_M, gravity=planet['gravity_m_s2'],
                            rotationperiod=planet['rotation_period_days'], year=365.25636,
                            obliquity=planet['obliquity_deg'], eccentricity=0.016715),
                gases_bar=dict(pN2=n2_ar / (1 + AR_PER_N2), pO2=O2_PA / 1e5, pAr=n2_ar * AR_PER_N2 / (1 + AR_PER_N2),
                               pCO2=co2))


def run_key(cfg, inputs) -> str:
    """Changes if the settings or any written input file change, so a resume never mixes models."""
    digest = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode())
    for path in (inputs['landmap'], inputs['topomap'], inputs['starspec'], inputs['starspec'][:-4] + '_hr.dat'):
        digest.update(Path(path).read_bytes())
    return digest.hexdigest()[:16]


def _write_sra(path: Path, code: int, field: np.ndarray):
    """A PlaSim surface file: an 8-integer header, then the field (north to south, from 0 E), 8 per line."""
    nlat, nlon = field.shape
    lines = [f'{code:10d}{0:10d}{20090101:10d}{0:10d}{nlon:10d}{nlat:10d}{0:10d}{0:10d}']
    flat = field.ravel()
    for i in range(0, flat.size, 8):
        lines.append(''.join(f'{v:10.3f}' for v in flat[i:i + 8]))
    path.write_text('\n'.join(lines) + '\n')


def land_mask(land_fraction, area, water_share):
    """PlaSim's binary mask: cells become sea in order of least land fraction until the sea covers
    `water_share` of the area (the last cell may overshoot by less than its own area)."""
    order = np.argsort(land_fraction, axis=None, kind='stable')
    covered = np.cumsum(area.ravel()[order]) / area.sum()
    sea = np.zeros(land_fraction.size, dtype=bool)
    sea[order[:int(np.searchsorted(covered, water_share)) + 1]] = True
    return (~sea).reshape(land_fraction.shape).astype(float)


def surface_files(cfg, folder: Path):
    """Land mask (code 172) and surface geopotential (code 129) at T21 from the geography product. With
    lakes, the water share is the product's (seas and lakes) and water cells stand at their water surface."""
    import netCDF4
    lakes = cfg.get('lakes', False)
    name = f"moon_{cfg['water']}pct_water_lakes_gaussian_T21.nc" if lakes else f"moon_{cfg['water']}pct_water_gaussian_T21.nc"
    with netCDF4.Dataset(PRODUCTS / name) as d:
        land = np.asarray(d['land_fraction'][:], dtype=float)
        height = np.asarray(d['land_elevation_m'][:], dtype=float)
        area = np.asarray(d['cell_area_fraction'][:], dtype=float)
        surface = np.asarray(d['water_surface_m'][:], dtype=float) if lakes else np.zeros_like(land)
    share = 1.0 - float((land * area).sum() / area.sum()) if lakes else cfg['water'] / 100.0
    mask = land_mask(land, area, share)
    geopotential = np.where(mask > 0, np.maximum(height, 0.0), np.maximum(surface, 0.0)) * cfg['planet']['gravity']
    folder.mkdir(parents=True, exist_ok=True)
    _write_sra(folder / 'moon_landmask.sra', 172, mask)
    _write_sra(folder / 'moon_topography.sra', 129, geopotential)
    return folder / 'moon_landmask.sra', folder / 'moon_topography.sra', float((mask * area).sum() / area.sum())


def filtered_spectrum(grid_um, w_nm, f_nm, shield_nm, shield_t):
    """Solar F_lambda (W m^-2 um^-1) through the shield on `grid_um`. Past the measured spectrum's end
    it falls off as the Sun's 5772 K blackbody, matched there: holding the last value instead would
    invent an infrared tail that ExoPlaSim counts in its near-infrared band."""
    h, c, k = 6.62607015e-34, 299792458.0, 1.380649e-23
    planck = lambda um: 1.0 / ((um * 1e-6) ** 5 * np.expm1(h * c / (um * 1e-6 * k * 5772.0)))
    w_um, f_um = w_nm / 1000.0, f_nm * 1000.0
    grid = np.asarray(grid_um, dtype=float)
    sun = np.where(grid <= w_um[-1], np.interp(grid, w_um, f_um), f_um[-1] * planck(grid) / planck(w_um[-1]))
    return sun * np.interp(grid * 1000.0, shield_nm, shield_t)


def sunlight(cfg, folder: Path):
    """Solar spectrum times the shield's transmission, in ExoPlaSim's two spectrum files; returns the flux."""
    from atmosphere.radiative_convective import fetch_inputs
    rows = []
    with open(fetch_inputs.path('TSIS1_HSRS_stride100.csv'), newline='') as handle:
        for row in csv.reader(handle):
            try:
                rows.append([float(row[0]), float(row[1])])
            except (ValueError, IndexError):
                pass
    w_nm, f = np.array(rows).T                                        # nm, W m^-2 nm^-1
    product = json.loads(SHIELDS.read_text())
    t = np.interp(w_nm, product['wavelength_nm'], product['transmission'][cfg['shield']])
    passed = float(np.trapezoid(f * t, w_nm))
    covered = float(np.trapezoid(f, w_nm))
    # Beyond the spectrum's range (mostly past 2.7 um) the solar constant's remainder passes at the
    # transmission of the longest measured wavelengths.
    tail = max(1361.0 - covered, 0.0) * float(np.interp(w_nm[-1], product['wavelength_nm'], product['transmission'][cfg['shield']]))
    flux = passed + tail
    import exoplasim
    source = Path(exoplasim.__file__).parent
    hi = np.concatenate([np.geomspace(0.2, 0.75, 1025)[:-1], np.geomspace(0.75, 100.0, 1024)])
    lo = np.loadtxt(source / 'wvref.txt')
    speed = 299792458.0                                               # ExoPlaSim's own scaling of these files
    folder.mkdir(parents=True, exist_ok=True)
    stem = folder / f"sun_{cfg['shield']}"
    for grid, name in ((hi, stem.name + '_hr'), (lo, stem.name)):
        values = filtered_spectrum(grid, w_nm, f, product['wavelength_nm'], product['transmission'][cfg['shield']]) / speed
        text = ' Wavelength    Flux  \n' + ''.join(f'{x} {y}\n' for x, y in zip(grid, values))
        (folder / f'{name}.dat').write_text(text)
    return f'{stem}.dat', flux


def prepare_inputs(cfg, rundir: Path) -> dict:
    folder = rundir / 'inputs'
    landmap, topomap, land_share = surface_files(cfg, folder)
    starspec, flux = sunlight(cfg, folder)
    return dict(landmap=str(landmap), topomap=str(topomap), starspec=starspec,
                flux_w_m2=flux * cfg['model']['sunlight_scale'], shield_flux_w_m2=flux, land_share=land_share)


def build_model(cfg, rundir: Path, ncpus: int, restart: Path | None, inputs: dict):
    import exoplasim as exo
    landmap, topomap, starspec, flux = inputs['landmap'], inputs['topomap'], inputs['starspec'], inputs['flux_w_m2']
    m = cfg['model']
    model = exo.Model(resolution=m['resolution'], layers=m['layers'], ncpus=ncpus, precision=8,
                      workdir=str(rundir / 'model'), modelname=f"moon_{cfg['experiment']}", outputtype='.nc',
                      crashtolerant=False)
    p = cfg['planet']
    model.configure(flux=flux, starspec=starspec, **cfg['gases_bar'], gravity=p['gravity'], radius=p['radius'],
                    rotationperiod=p['rotationperiod'], synchronous=False, year=p['year'],
                    eccentricity=p['eccentricity'], obliquity=p['obliquity'], landmap=str(landmap),
                    topomap=str(topomap), mldepth=cfg['mldepth'], seaice=True, ozone=False,
                    soilalbedo=m['land_albedo'], timestep=m['timestep_min'], snapshots=None,
                    restartfile=str(restart) if restart else None)
    # The namelist holds file names in 80 characters, too few for an absolute path here: the model
    # reads the spectrum from short names in its own working folder instead.
    workdir = Path(model.workdir)
    shutil.copy2(starspec, workdir / 'sun.dat')
    shutil.copy2(starspec[:-4] + '_hr.dat', workdir / 'sun_hr.dat')
    model._edit_namelist('radmod_namelist', 'STARFILE', "'sun.dat'")
    model._edit_namelist('radmod_namelist', 'STARFILEHR', "'sun_hr.dat'")
    # A fixed seed makes runs repeatable; the generator's state then travels in every restart file.
    model._edit_namelist('plasim_namelist', 'SEED', str(m['seed']))
    model._edit_namelist('plasim_namelist', 'NSTPW', str(m['steps_per_lunar_day'] // m['writes_per_lunar_day']))
    model._edit_namelist('radmod_namelist', 'RAYSCALE', str(m['rayleigh_scale']))
    model._edit_namelist('radmod_namelist', 'TH2OC', str(m['h2o_continuum']))
    gref, scale = cloud_water_namelist(m['cloud_water'], p['gravity'])
    model._edit_namelist('rainmod_namelist', 'CWGREF', str(gref))
    model._edit_namelist('rainmod_namelist', 'CWSCALE', str(scale))
    model._edit_namelist('plasim_namelist', 'NDIAGCF', str(m['clear_sky']))
    model._edit_namelist('rainmod_namelist', 'CONVDAY', str(float(m['convective_day_s'])))
    from exoplasim import pyburn
    for code, (name, long_name) in CLEAR_SKY.items():
        pyburn.ilibrary[str(code)] = [name, long_name, 'W m-2']
        pyburn.slibrary[name] = [code, long_name, 'W m-2']
    model.cfgpostprocessor(ftype='regular', extension='.nc', variables=cfg['output'],
                           times=m['writes_per_lunar_day'] * m['lunar_days_per_year'], timeaverage=True,
                           interpolatetimes=False)
    return model


def _free_gb(path: Path) -> float:
    return shutil.disk_usage(path).free / 1e9


def _clean_partial(workdir: Path, year: int):
    """Remove whatever an interrupted year left in the working folder."""
    for name in ('plasim_output', 'plasim_status', 'plasim_snapshot', 'plasim_diag', 'plasim_hcadence', 'Abort_Message'):
        (workdir / name).unlink(missing_ok=True)
    for p in workdir.glob(f'MOST*.{year:05d}*'):
        p.unlink()


def _prune(workdir: Path, year: int):
    for p in workdir.glob('MOST_REST.*'):
        y = int(p.name.split('.')[1])
        if y <= year - KEEP_RESTARTS and y % 10 != 9:
            p.unlink()
    for p in workdir.glob('MOST.?????.nc'):
        y = int(p.name.split('.')[1])
        if y <= year - KEEP_OUTPUT_YEARS and y % 10 != 9:
            p.unlink()


def model_reported_sunlight(workdir: Path) -> dict:
    """The band split and Rayleigh coefficient ExoPlaSim derived from the spectrum (its diagnostic file)."""
    out = {}
    for path in sorted(workdir.glob('MOST_DIAG.*'))[:1]:
        for line in path.read_text(errors='replace').splitlines():
            for key, label in (('energy_below_0.75um', 'Energy fraction below 0.75 microns:'),
                               ('rayleigh_coefficient', 'Rayleigh scattering coefficient:')):
                if label in line:
                    out[key] = float(line.split(':')[1])
    return out


def _members(pgid: int) -> list:
    """Processes in a process group, read from /proc."""
    out = []
    for stat in Path('/proc').glob('[0-9]*/stat'):
        try:
            fields = stat.read_text().rsplit(')', 1)[1].split()       # state, ppid, pgrp, ...
            if int(fields[2]) == pgid:
                out.append(int(stat.parent.name))
        except (OSError, IndexError, ValueError):
            pass
    return out


def _using(workdir: Path) -> list:
    """Processes whose working folder is this run's model folder: ExoPlaSim's MPI ranks and mpiexec."""
    out = []
    for cwd in Path('/proc').glob('[0-9]*/cwd'):
        try:
            if Path(os.readlink(cwd)) == workdir.resolve():
                out.append(int(cwd.parent.name))
        except OSError:
            pass
    return out


def _terminate(pids, grace=30.0):
    """SIGTERM, then SIGKILL whatever is left after `grace` seconds."""
    pids = [p for p in set(pids) if p != os.getpid()]
    for pid in pids:
        try:
            os.kill(pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    deadline = time.time() + grace
    while time.time() < deadline and any(Path(f'/proc/{p}').exists() for p in pids):
        time.sleep(0.5)
    for pid in pids:
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False


def stop(name: str, folder: str | None = None) -> int:
    rundir = RUNS / (folder or name)
    lock = rundir / 'run.lock'
    if not lock.exists():
        print(f'{rundir.name}: not running')
        return 0
    pid, pgid = (int(v) for v in lock.read_text().split()[:2])
    if _alive(pid):
        os.kill(pid, signal.SIGTERM)                               # the runner stops its own group
        deadline = time.time() + 60.0
        while time.time() < deadline and _alive(pid):
            time.sleep(0.5)
    _terminate(_members(pgid) + _using(rundir / 'model'))
    lock.unlink(missing_ok=True)
    print(f'{rundir.name}: stopped')
    return 0


def _on_signal(signum, frame):
    """Take the MPI ranks down with the runner, then let the finally blocks clean up."""
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, signal.SIG_IGN)
    if os.getpgrp() == os.getpid():
        try:
            os.killpg(os.getpid(), signal.SIGTERM)
        except OSError:
            pass
    sys.exit(128 + signum)

def summarise_year(path: Path) -> dict:
    """Headline numbers from one year of output (3-day means), for watching the run settle at a glance."""
    import netCDF4
    with netCDF4.Dataset(path) as d:
        clear = 'rstcs' in d.variables
        keys = ('ts', 'tas', 'ntr', 'clt', 'pr', 'sic', 'lsm', 'rst', 'rsut', 'rlut') + (('rstcs', 'rltcs') if clear else ())
        v = {k: np.asarray(d[k][:], dtype=float) for k in keys}
        lat = np.asarray(d['lat'][:], dtype=float)
        times = d['time'].size
    expected = MODEL['writes_per_lunar_day'] * MODEL['lunar_days_per_year']
    if times != expected:
        raise RuntimeError(f'{path.name}: {times} output times, expected {expected}; the calendar is not as assumed')
    nodes, weights = np.polynomial.legendre.leggauss(lat.size)       # Gaussian weights of the T21 rows
    w = weights[np.argsort(np.argsort(np.sin(np.radians(lat))))][:, None] * np.ones(v['ts'].shape[-1])
    mean = lambda x, m=None: float((x * w * (1 if m is None else m)).sum() / (w * (1 if m is None else m)).sum())
    year = lambda k: v[k].mean(axis=0)
    land = v['lsm'][0] > 0.5
    swing = v['tas'].max(axis=0) - v['tas'].min(axis=0)
    equator = (np.abs(lat) < 12.0)[:, None] & land
    pole = (np.abs(lat) > 70.0)[:, None] * np.ones_like(land)
    # rsut is the reflected sunlight, negative (upward); rst the net, so the incoming is rst - rsut.
    albedo = -mean(year('rsut')) / mean(year('rst') - year('rsut'))
    clouds = {}
    if clear:                                   # the clouds' effect on the top-of-atmosphere net flux
        sw, lw = mean(year('rst') - year('rstcs')), mean(year('rlut') - year('rltcs'))
        clouds = dict(cloud_effect_sw_w_m2=round(sw, 2), cloud_effect_lw_w_m2=round(lw, 2),
                      cloud_effect_net_w_m2=round(sw + lw, 2))
    return dict(
        surface_k=round(mean(year('ts')), 2), air_2m_k=round(mean(year('tas')), 2),
        toa_net_w_m2=round(mean(year('ntr')), 2), planetary_albedo=round(albedo, 4), **clouds,
        cloud_cover=round(mean(year('clt')), 3),
        precipitation_mm_day=round(mean(year('pr')) * 86400e3, 3),
        sea_ice_share_of_sea=round(mean(year('sic'), ~land), 3),
        warmest_air_k=round(float(v['tas'].max()), 1), coldest_air_k=round(float(v['tas'].min()), 1),
        equator_land_day_night_swing_k=round(mean(swing, equator), 1),
        equator_land_air_k=round(mean(year('tas'), equator), 1),
        polar_air_k=round(mean(year('tas'), pole), 1),
        output_mb=round(path.stat().st_size / 1e6, 1))


def run(name: str, years: int, ncpus: int, folder: str | None = None, prune: bool = False,
        overrides: dict | None = None, start_from: str | None = None) -> int:
    cfg = configuration(name, ncpus, overrides)
    rundir = RUNS / (folder or name)
    rundir.mkdir(parents=True, exist_ok=True)
    lock = rundir / 'run.lock'
    if lock.exists():
        pid, pgid = (int(v) for v in lock.read_text().split()[:2])
        if _alive(pid):
            print(f'{rundir.name} is already running (process {pid}); not starting a second runner.', file=sys.stderr)
            return 1
        _terminate(_members(pgid) + _using(rundir / 'model'))     # orphans of a killed runner
        lock.unlink()
    leftover = _using(rundir / 'model')
    if leftover:
        _terminate(leftover)
    try:
        os.setpgid(0, 0)                                          # lead a group the ranks will share
    except PermissionError:
        pass                                                      # already a session leader
    _atomic_write(lock, f'{os.getpid()} {os.getpgrp()} {time.strftime("%Y-%m-%d %H:%M:%S")}\n')
    progress_path = rundir / 'progress.json'
    try:
        inputs = prepare_inputs(cfg, rundir)
        stored = json.loads(progress_path.read_text()) if progress_path.exists() else None
        # A branch starts from another run's restart; its identity is part of the run.
        if start_from:
            cfg['initial_state'] = dict(path=str(start_from), sha256=hashlib.sha256(Path(start_from).read_bytes()).hexdigest()[:16])
        elif stored:
            cfg['initial_state'] = stored['configuration'].get('initial_state', 'exoplasim default')
        else:
            cfg['initial_state'] = 'exoplasim default'
        cfg['key'] = run_key(cfg, inputs)
        progress = stored or dict(configuration=cfg, inputs=inputs, years=[])
        if progress['configuration']['key'] != cfg['key']:
            print(f'{rundir.name}: the model changed since this run began (key '
                  f"{progress['configuration']['key']} vs {cfg['key']}); start it in a new folder.", file=sys.stderr)
            return 1
        done = len(progress['years'])
        workdir = rundir / 'model'
        restart = workdir / f'MOST_REST.{done - 1:05d}' if done else None
        if not done and isinstance(cfg['initial_state'], dict):
            restart = Path(cfg['initial_state']['path'])
        if restart is not None and not restart.is_file():
            print(f'{name}: restart for year {done - 1} is missing; cannot resume.', file=sys.stderr)
            return 1
        if workdir.exists():
            _clean_partial(workdir, done)
        if restart is not None:                                       # keep it safe from configure()'s copy
            kept = rundir / 'resume_from'
            shutil.copy2(restart, kept)
            restart = kept
        model = build_model(cfg, rundir, ncpus, restart, inputs)
        model.currentyear = done
        for year in range(done, years):
            if _free_gb(rundir) < MIN_FREE_GB:
                print(f'{name}: less than {MIN_FREE_GB} GB free; pausing after year {year - 1}.', file=sys.stderr)
                return 2
            start = time.time()
            model.run(years=1, crashifbroken=True)
            out, rest = workdir / f'MOST.{year:05d}.nc', workdir / f'MOST_REST.{year:05d}'
            if not (out.is_file() and rest.is_file()):
                print(f'{name}: year {year} produced no output or restart; stopping.', file=sys.stderr)
                return 1
            entry = dict(year=year, seconds=round(time.time() - start, 1), **summarise_year(out))
            if 'model_sunlight' not in progress:
                # The Sun's light splits about evenly across 0.75 um; anything far off means a broken spectrum.
                progress['model_sunlight'] = model_reported_sunlight(workdir)
                split = progress['model_sunlight'].get('energy_below_0.75um', float('nan'))
                if not 0.3 < split < 0.7:
                    print(f'{rundir.name}: the model puts {split} of the sunlight below 0.75 um; the spectrum '
                          'is wrong. Stopping.', file=sys.stderr)
                    return 1
            progress['years'].append(entry)
            _atomic_write(progress_path, json.dumps(progress, indent=1) + '\n')
            if prune:
                _prune(workdir, year)
            print(f'{name} year {year}: {entry}', flush=True)
        return 0
    finally:
        lock.unlink(missing_ok=True)


def status(folder: str) -> int:
    rundir = RUNS / folder
    lock = rundir / 'run.lock'
    running = lock.exists() and _alive(int(lock.read_text().split()[0]))
    path = rundir / 'progress.json'
    years = json.loads(path.read_text())['years'] if path.exists() else []
    state = 'running' if running else ('stopped' if years else 'not started')
    print(f"{folder}: {len(years)} years done; {state}" + (f', year {len(years)} in progress' if running else ''))
    for entry in years[-5:]:
        print(' ', entry)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('experiment', choices=sorted(EXPERIMENTS))
    parser.add_argument('--years', type=int, default=40)
    parser.add_argument('--ncpus', type=int, default=8)
    parser.add_argument('--status', action='store_true')
    parser.add_argument('--stop', action='store_true')
    parser.add_argument('--folder', help='run folder under climate/gcm/runs (default: the experiment name)')
    parser.add_argument('--prune', action='store_true', help='keep only recent and every tenth year')
    parser.add_argument('--set', action='append', default=[], metavar='KEY=VALUE',
                        help='override a model setting, e.g. rayleigh_scale=1.5 (part of the run key)')
    parser.add_argument('--start-from', help='restart file to branch a new run from')
    args = parser.parse_args(argv)
    if args.status:
        return status(args.folder or args.experiment)
    if args.stop:
        return stop(args.experiment, args.folder)
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(sig, _on_signal)
    overrides = {}
    for item in args.set:
        key, value = item.split('=', 1)
        if key not in MODEL:
            parser.error(f'unknown model setting {key}')
        overrides[key] = type(MODEL[key])(value)
    return run(args.experiment, args.years, args.ncpus, args.folder, args.prune, overrides, args.start_from)


if __name__ == '__main__':
    sys.exit(main())
