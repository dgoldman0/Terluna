"""Run the Open Moon's ExoPlaSim experiments one simulated year at a time, resumable after any stop.

    climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --years 40          # run or resume A
    climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --status            # progress so far
    climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --stop              # stop it cleanly
    climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --years 15 --folder A_full_column_clouds \
        --set cloud_water=full_column --start-from climate/gcm/runs/A/model/MOST_REST.00039   # a branch
    climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A28 --years 20 --folder A28_gases_glow \
        --co2-ppm 755 --top-infrared 10.1718 0.5165 --set sunlight_scale=0.909867 \
        --start-from climate/gcm/runs/A28_dim5_moon/model/MOST_REST.00029        # CO2 and infrared from space

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
is a uniform placeholder. --co2-ppm sets the CO2 the air carries instead of 400 ppm, with N2 the balance; PlaSim's
radiation knows only water vapour, CO2 and ozone, so a run can carry another gas's warming as the CO2 that warms as
much (its CO2 equivalent), and the run records the choice. --top-infrared adds infrared from space at the top of the
air, by day and by night (TOP_INFRARED_EDITS). PlaSim counts 12 of the Moon's rotations in a 360-day year, which makes
its model day 29.8 Earth days (the rotation keeps the true 27.32 days); a model year is 12 of them, and the
output ten means per model day. The random seed is fixed, so runs repeat exactly. A patch to PlaSim's
solang turns the Sun evenly at the synodic rate, one turn per 29.53 days, and the orbit is held to the
Moon's elements with its 1.54 degree tilt. Runs before 2026-09-29 had PlaSim's own clock, which sweeps the
Sun unevenly, and Earth's orbit for 1 AD with a 23.7 degree tilt: see the README ("How the Moon is set
up"). Since 2026-10-14 the corrected radiation is the standard, as the author decided: new runs take the curved air's
light and infrared (--twilight full) and CO2 and O2 in PlaSim's sunlight (--gas-bands), with the water and Rayleigh
multipliers refitted with them, in the twilight copy of ExoPlaSim; --twilight none and --no-gas-bands leave them out.
Every earlier run had neither, so none of them resumes under the new defaults. Since 2026-10-19 the curved air's terms
come from the ledger on the design climate's 301 K column (--twilight-ledger 301k; 294.9k gives the earlier one) and
each layer's are taken in the clear share above it (--twilight-clear-share, TWILIGHT_CLEAR_EDITS), and experiment A28R
reads the lakes surface built from the geography domain's current routing (boundary.lakes_name), which needs a cold
start. Since 2026-10-26 --levels l14 runs PlaSim on fourteen levels with more of them near the ground (levels.py), its
restart carried there by levels.py's conservative remap and its twilight terms, longwave correction and cold cloud's
factors on those levels, and --sbm-shallow-cap holds a shallow cloud's top within the depth test's threshold of its own
base in place of PlaSim's 700 hPa (SBM_CAP_PATCHES); without either a run is as before. --levels l15 adds a thin bottom
layer, the ground to 0.6 km, to the fourteen, and --longwave-fit l15_provisional gives a run on it the longwave
correction refitted on its own levels' columns. --convection massflux runs the mass-flux scheme of massflux.py in
Kuo's place (rainmod's nprc 3; massflux_fortran.py's patch, applied only for a run that sets it); a run without it is as
before. --free-convection lunar takes the free-convection constant in the heat and vapour transfer over water (Miller
and others 1992; FREE_CONVECTION_EDITS) from Earth's 0.0016 to the Moon's gravity and surface pressure, and
--massflux-start boundary_layer starts the mass-flux scheme's surface parcel at a multiple of the boundary layer's
convective velocity (massflux_fortran.START_PATCHES). --mixing-length lunar stretches PlaSim's Blackadar length in its
vertical diffusion (fluxmod's vdiff_lamm, 160 m) by the ratio of the gravities, and --boundary-layer nonlocal runs the
nonlocal boundary layer of pbl_scheme.py in vdiff (pbl_fortran.PBL_PATCHES, fluxmod's nlpbl); a run with either also
writes the vapour each process gives each layer and the heating of falling precipitation's evaporation (BUDGET), with
PlaSim's convection slot completed for the Betts-Miller and mass-flux schemes (pbl_fortran.DIAGNOSTIC_PATCHES).
--surface-gust beljaars adds Beljaars's gust, the convective velocity of the boundary layer's eddies, to the wind in
PlaSim's surface transfer over every surface, in place of the free-convection term over water and the least wind of
1 m/s (surface_gust.py; GUST_PATCHES over the nonlocal scheme's patch, whose depth rule it takes); a run with it and
the moisture diagnostics writes the vapour by process too.
--sbm-diagnostic-night, a diagnostic and never a configuration, holds the Betts-Miller scheme's deep convection off
wherever the Sun is down (SBM_NIGHT_PATCHES); a run with it and the moisture diagnostics writes the vapour by process
too. --large-scale-condensate cloud holds PlaSim's large-scale condensate as cloud water in its layer above freezing,
evaporated where the air is unsaturated and turned to rain by autoconversion and accretion (condensate.py), and
--humidity-mean writes the humidity's true 3-day mean as hus_mean (HUS_MEAN_PATCHES); without them a run is as
before.
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
from shared.constants import EXOPLASIM_FREE_CONVECTION, STANDARD_ATMOSPHERE, STANDARD_GRAVITY
from climate.gcm import plasim_ice_optics as pio
from climate.gcm import levels as lv
from climate.gcm import massflux_fortran as mff
from climate.gcm import pbl_fortran as pbf
from climate.gcm import pbl_scheme as pbs
from climate.gcm import surface_gust as sgf
from climate.gcm import condensate as cnd

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
    # The same seas with the lakes of the routing the geography domain carries since 2026-09-29 (the A28_dim5_moon
    # climate), on the surface named by that drainage product's content (boundary.lakes_name): 10.34% of the Moon in
    # lakes against run A's 12.02%. A28's surface, run A's lakes, stays as every A28 run read it. A new coastline
    # needs a cold start.
    'A28R': dict(pressure_pa=121590.0, water=28, lakes='98842fc2a648f889', shield='titania_stack', mldepth=50.0,
                 purpose='design case with the atlas seas and the rain-fed lakes of the current routing'),
}
# PlaSim's model day is 1430 half-hour steps (29.8 Earth days; with the synodic Sun clock the Sun's own day
# is the real 29.53). A model year is 12 of them, 17,160 steps (whole_lunar_days; ExoPlaSim's own year is
# 17,280). Output: 10 means per model day (about 3 days each), 120 a year.
MODEL = dict(resolution='T21', layers=10, timestep_min=30.0, land_albedo=0.2, steps_per_lunar_day=1430,
             writes_per_lunar_day=10, lunar_days_per_year=12, seed=1,
             # Radiation calibrated against the line-by-line model (see the README): a multiplier on ExoPlaSim's
             # spectrum-derived Rayleigh coefficient and PlaSim's water-vapour continuum coefficient, first
             # fitted on cloud-free lunar columns at 280 and 310 K (1.8 and 0.004, runs before 2026-09-29),
             # the Rayleigh multiplier then refitted beside the water absorption below. 1.0 and 0.024
             # reproduce the unmodified model, which absorbed about 14 W/m2 too much sunlight and emitted
             # 12 W/m2 too little when warm.
             rayleigh_scale=2.237, h2o_continuum=0.004,
             # A multiplier on water vapour's absorptivity in PlaSim's near-infrared band (SHORTWAVE_EDITS below;
             # 1.0 is PlaSim's own). With the Rayleigh multiplier refitted beside it, it brings PlaSim's clear-sky
             # split of the sunlight to line-by-line's on 11 columns of the design case (plasim_shortwave.py),
             # with the Moon's slant factors and the CO2 and O2 bands since 2026-10-14
             # (../results/gcm/shortwave_fit_A28_dim5_lunar.json); runs before had 1.356 and 2.215, fitted without.
             sw_water_absorption=1.292,
             # How the radiation sees cloud water (CLOUD_WATER below). Runs A and B used 'plasim', ExoPlaSim
             # unmodified; 'earth_path' corrects its Earth-fitted formula for the Moon's gravity.
             cloud_water='earth_path',
             # Seconds per day in PlaSim's convective cloud cover, 0.245 + 0.125 ln(convective rain in mm per
             # day), capped at 0.8 (Slingo's fit to Earth). PlaSim divides the rain by its own solar day, 29.8
             # Earth days on the Moon, so it reads the rain about 30 times too heavy and raises the cover by
             # 0.42, up to the cap. 86400 counts it per Earth day, as fitted; 0 keeps PlaSim's behaviour (runs
             # A and B and the cloud-water bracket).
             convective_day_s=86400.0,
             # PlaSim's critical relative humidity for cloud in its lowest layer (sigma 0.983, the ground to about
             # 1.8 km), above which cover grows as ((rh - rcrit)/(1 - rcrit))^2. 0 keeps PlaSim's own,
             # max(0.85, sigma) = 0.983; a lower one gives low cloud and fog more often (the night fog test of
             # 2026-10-01, against CM1's foggy nights).
             surface_cloud_rcrit=0.0,
             # PlaSim's clear-sky diagnostic (1 on, 0 off): the radiation computed a second time without
             # clouds, for the clouds' radiative effect. It changes no prognostic field.
             clear_sky=1,
             # A multiplier on the sunlight the shield passes, the same at every wavelength: 0.92 is a shield
             # that passes 8% less. 1.0 keeps the shield product's transmission.
             sunlight_scale=1.0,
             # The Sun's clock (the Terluna patch to solang, SUN_CLOCKS below): 'synodic' turns the Sun
             # evenly, once per 29.53-day lunar day. 'plasim' keeps PlaSim's own clock, which sweeps the Sun
             # 392.5 degrees per 29.8-day model day and skips it back 32.5 (README, "How the Moon is set
             # up"); every run before 2026-09-29 had it.
             sun_clock='synodic',
             # 1 holds the orbit to the planet's own elements (ExoPlaSim's fixedorbit): the Moon's 1.54 degree
             # tilt. With 0, ExoPlaSim's default, PlaSim replaces them with Earth's orbit for its start year,
             # 1 AD, and its 23.7 degree tilt; every run before 2026-09-29 had that.
             fixed_orbit=1,
             # 1 makes each model year exactly lunar_days_per_year model days (17,160 steps), all of it in
             # the output's 120 means. 0 keeps ExoPlaSim's 360-day year (17,280 steps), whose last 120 steps
             # fall outside them.
             whole_lunar_days=1)
# The Sun's clock, PlaSim's NSUNCLOCK (the index): its own; one turn per model day of steps_per_lunar_day
# steps (29.8 days); or the synodic rate the rotation and the orbit give, continuous in time.
SUN_CLOCKS = ('plasim', 'model_day', 'synodic')
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
# Kept in each year's 3-day means (ExoPlaSim writes about 100 fields by default, 100 MB a year). The surface stress,
# tauu and tauv (PlaSim's codes 180 and 181, means over each output interval), joined on 2026-10-19 so that the waters'
# studies can set the 10 m wind from the model's own stress beside the neutral profile from the lowest level. ExoPlaSim's
# postprocessor also names uas and vas, the 10 m wind (codes 165 and 166), but PlaSim 3.4.2 writes neither (outmod.f90),
# and the postprocessor then leaves a requested name out of the file without an error, so they are not asked for;
# test_exoplasim_run checks every name here against the codes PlaSim writes.
OUTPUT = ['ta', 'ua', 'va', 'hus', 'cl',
          'ts', 'tas', 'maxt', 'mint', 'ps', 'psl', 'pr', 'prc', 'evap', 'prw', 'clt', 'sic', 'sit', 'snd', 'mrso',
          'rst', 'rsut', 'rlut', 'ntr', 'nbr', 'rss', 'rls', 'hfss', 'hfls', 'alb', 'czen', 'lsm', 'tauu', 'tauv']
# The codes PlaSim 3.4.2 writes to its raw output (outmod.f90's writegp calls on unit 40), and the names ExoPlaSim's
# postprocessor derives from them (pyburn.dataset): a name outside both comes out of the postprocessor as nothing.
# The diagnostics' own codes (the heating's 21-29 among them) come on top when their switches are on, as does 233,
# the humidity's true mean, which HUS_MEAN_PATCHES writes through a variable only with nhusmean 1.
PLASIM_CODES = (101, 102, 103, 104, 105, 106, 107, 110, 129, 130, 133, 138, 139, 140, 141, 142, 143, 144, 146, 147, 152,
                155, 159, 160, 161, 162, 164, 165, 166, 167, 169, 170, 172, 173, 174, 175, 176, 177, 178, 179, 180,
                181, 182, 183, 184, 201, 202, 203, 204, 205, 207, 208, 209, 210, 211, 218, 221, 229, 230, 232, 265,
                266, 267, 268, 269, 298, 318, 319, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 410, 411)
DERIVED = ('ua', 'va', 'spd', 'ps', 'psl', 'pr', 'ntr', 'nbr', 'hfns', 'wfn', 'zg', 'hur', 'wap', 'wa', 'stf', 'psi')
# PlaSim's clear-sky fluxes, written when the diagnostic is on but unknown to ExoPlaSim's postprocessor,
# which the runner teaches them. Unlike the other fluxes they are snapshots at each output time, not
# means over it: sound for global and lunar-day averages, not for single times.
CLEAR_SKY = {101: ('rsscs', 'surface_net_shortwave_flux_clear_sky'),
             102: ('rlscs', 'surface_net_longwave_flux_clear_sky'),
             103: ('rstcs', 'toa_net_shortwave_flux_clear_sky'),
             104: ('rltcs', 'toa_net_longwave_flux_clear_sky')}
# PlaSim's gridpoint diagnostics (NDIAGGP = 1) write the temperature tendency of each process at every level as
# codes 21-29, K/s: the sum of all of them that the gridpoint diabatic step (gridpointd) hands the dynamics, the
# nudging of the top level (mknudge), the surface's sensible heat into the lowest layer (mkshfl), vertical diffusion
# (vdiff), longwave and shortwave radiation (radmod), large-scale condensation (mklsp), Kuo convection (kuo) and dry
# convective adjustment (mkdca). The dynamics' own tendency is not among them.
HEATING = {21: ('dt_diab', 'temperature_tendency_summed_over_the_physics'),
           22: ('dt_nudge', 'temperature_tendency_due_to_top_level_nudging'),
           23: ('dt_sfc', 'temperature_tendency_due_to_surface_sensible_heat'),
           24: ('dt_vdif', 'temperature_tendency_due_to_vertical_diffusion'),
           25: ('dt_lw', 'temperature_tendency_due_to_longwave_radiation'),
           26: ('dt_sw', 'temperature_tendency_due_to_shortwave_radiation'),
           27: ('dt_lsp', 'temperature_tendency_due_to_large_scale_condensation'),
           28: ('dt_conv', 'temperature_tendency_due_to_convection'),
           29: ('dt_dca', 'temperature_tendency_due_to_dry_convective_adjustment')}
# With the same diagnostics PlaSim also writes the water vapour each layer gains, as a water flux (m/s: the specific
# humidity tendency times the layer's mass over the density of water): code 31 the sum over the physics, code 32 the fix
# that fills negative humidity (miscmod), and codes 40 and 41 the two parts of the spectral step that changes the
# humidity after the physics (spectrald, mkdqtgp): the physics' moistening as that step applies it, carried to the
# model's harmonics and back, then the horizontal diffusion, which sums to nothing over the globe. Code 40 keeps the
# name dq_dyn under which the runs have written it since 21 October, but no code holds the dynamics' own transport: the
# circulation's moistening is the storage less code 40 (climate/crm/gcm_site.py, gcm/night_water.py). With them the
# runner writes omega (wap), which the postprocessor derives from the spectral fields. A column model forced by the GCM
# takes these as the GCM's own budget, which closes by construction (climate/crm/gcm_site.py).
MOISTENING = {31: ('dq_diab', 'water_vapour_gained_in_the_layer_summed_over_the_physics'),
              32: ('dq_fill', 'water_vapour_gained_in_the_layer_by_the_negative_humidity_fix'),
              40: ('dq_dyn', 'water_vapour_gained_in_the_layer_by_the_physics_as_the_spectral_step_applies_it'),
              41: ('dq_dyn_late', 'water_vapour_gained_in_the_layer_by_horizontal_diffusion')}
# A run with the boundary-layer options or the night switch (a diagnostic) and the moisture diagnostics also writes the
# vapour each process gives each layer (2026-10-27), from slots PlaSim fills whenever its gridpoint diagnostics are on:
# the surface's evaporation into the lowest layer (mkevap), vertical diffusion (vdiff), large-scale condensation
# (mklsp), convection (Kuo's routine, and with pbl_fortran.DIAGNOSTIC_PATCHES the Betts-Miller and mass-flux schemes)
# and dry adjustment (mkdca), each a water flux as code 31 is, and the evaporation of falling rain and snow (mkrain),
# whose vapour PlaSim writes as a specific humidity tendency and whose heating is code 30.
BUDGET = {30: ('dt_evap', 'temperature_tendency_due_to_evaporation_of_falling_precipitation', 'K s-1'),
          33: ('dq_sfc', 'water_vapour_gained_in_the_lowest_layer_from_the_surface', 'm s-1'),
          34: ('dq_vdif', 'water_vapour_gained_in_the_layer_by_vertical_diffusion', 'm s-1'),
          35: ('dq_lsp', 'water_vapour_gained_in_the_layer_by_large_scale_condensation', 'm s-1'),
          36: ('dq_conv', 'water_vapour_gained_in_the_layer_by_convection', 'm s-1'),
          37: ('dq_dca', 'water_vapour_gained_in_the_layer_by_dry_convective_adjustment', 'm s-1'),
          39: ('dq_evap', 'specific_humidity_tendency_due_to_evaporation_of_falling_precipitation', 'kg kg-1 s-1')}


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
# The Sun's clock (radmod.f90, solang): NSUNCLOCK 0, the default, keeps PlaSim's hour angle, which turns at
# the rotation rate and restarts every model day. 1 turns the Sun once per model day. 2 turns it at the
# synodic rate, 1/sidereal_day - 1/sidereal_year, continuously from step 0 (for runs that are not
# synchronous; slowdown does not apply). Reported upstream as a bug in solang (README).
SUN_EDITS = [
    ("      real :: rayscale = 1.0     ! Terluna: multiplier on rcoeff, calibrated against line-by-line\n",
     "      real :: rayscale = 1.0     ! Terluna: multiplier on rcoeff, calibrated against line-by-line\n"
     "      integer :: nsunclock = 0   ! Terluna: the Sun's clock (0 PlaSim's, 1 a turn per model day, 2 synodic)\n"),
    (",nsimplealbedo,nstarfile,starfile,starfilehr,minwavel,rayscale\n",
     ",nsimplealbedo,nstarfile,starfile,starfilehr,minwavel,rayscale,nsunclock\n"),
    ("      call mpbcr(rayscale)\n", "      call mpbcr(rayscale)\n      call mpbci(nsunclock)\n"),
    ("      zrtim = rotspd * TWOPI / 1440.0         ! scale time   to radians\n",
     "      zrtim = rotspd * TWOPI / 1440.0         ! scale time   to radians\n"
     "      if (nsunclock == 1) zrtim = TWOPI / (ntspd * mpstep) ! Terluna: one turn per model day\n"
     "      if (nsunclock == 2) zturns = mod(real(nstep) * mpstep * 60.0                              &\n"
     "     &    * (1.0 / sidereal_day - 1.0 / sidereal_year), 1.0) ! Terluna: synodic turns since step 0\n"),
    ("         if (ngenkeplerian==1) zhangle = zhangle - rasc\n",
     "         if (ngenkeplerian==1) zhangle = zhangle - rasc\n"
     "         if (nsunclock == 2) zhangle = zturns * TWOPI + jlon * zrlon - PI ! Terluna: synodic clock\n"),
]
# Water vapour's sunlight (radmod.f90, swr): SWH2OSCALE multiplies the absorptivity of PlaSim's near-infrared
# band, Lacis and Hansen's formula fitted to Earth's water paths, which absorbs too little over the Moon's much
# longer ones; the band's transmission is kept above 0.001. 1, the default, keeps PlaSim's. The form was chosen
# with plasim_shortwave.py, which holds the fit.
SHORTWAVE_EDITS = [
    ("      integer :: nsunclock = 0   ! Terluna: the Sun's clock (0 PlaSim's, 1 a turn per model day, 2 synodic)\n",
     "      integer :: nsunclock = 0   ! Terluna: the Sun's clock (0 PlaSim's, 1 a turn per model day, 2 synodic)\n"
     "      real :: swh2oscale = 1.0   ! Terluna: multiplier on water vapour's shortwave absorptivity\n"),
    (",minwavel,rayscale,nsunclock\n", ",minwavel,rayscale,nsunclock,swh2oscale\n"),
    ("      call mpbci(nsunclock)\n", "      call mpbci(nsunclock)\n      call mpbcr(swh2oscale)\n"),
    ("       ztwvtu(:)=1.-2.9*zwv(:)/((1.+141.5*zwv(:))**0.635+5.925*zwv(:))  &\n"
     "     &            /zsolar2\n",
     "       ztwvtu(:)=max(0.001,1.-swh2oscale*2.9*zwv(:)/((1.+141.5*zwv(:))**0.635+5.925*zwv(:))  &\n"
     "     &            /zsolar2) ! Terluna: calibrated absorptivity\n"),
    ("       ztwv(:)=(1.-2.9*zwv(:)/((1.+141.5*zwv(:))**0.635+5.925*zwv(:))   &\n",
     "       ztwv(:)=max(0.001,1.-swh2oscale*2.9*zwv(:)/((1.+141.5*zwv(:))**0.635+5.925*zwv(:))   &\n"),
    ("     &         /(1.-2.9*zwv(:)/((1.+141.5*zwv(:))**0.635+5.925*zwv(:))  &\n",
     "     &         /max(0.001,1.-swh2oscale*2.9*zwv(:)/((1.+141.5*zwv(:))**0.635+5.925*zwv(:))  &\n"),
]
# Infrared arriving from space at the top of the air, W/m2 over the day and the night hemispheres: the ring fleet's
# glow and the Moon's own infrared its films return. PlaSim's longwave scheme takes the downward flux at the model's
# top as zbd(:,0), the emission of the top layer (its ExoPlaSim note marks the place for an outside source), and
# carries it down with the column's broadband transmissivities, so what enters there is absorbed in the air as the
# air's own infrared is and the rest reaches the ground. 0 by day and night, the default, leaves the model as it was.
TOP_INFRARED_EDITS = [
    ("      real :: swh2oscale = 1.0   ! Terluna: multiplier on water vapour's shortwave absorptivity\n",
     "      real :: swh2oscale = 1.0   ! Terluna: multiplier on water vapour's shortwave absorptivity\n"
     "      real :: irtopday = 0.0     ! Terluna: infrared from space at the top of the air by day, W/m2\n"
     "      real :: irtopnight = 0.0   ! Terluna: and by night\n"),
    (",minwavel,rayscale,nsunclock,swh2oscale\n", ",minwavel,rayscale,nsunclock,swh2oscale,irtopday,irtopnight\n"),
    ("      call mpbcr(swh2oscale)\n", "      call mpbcr(swh2oscale)\n      call mpbcr(irtopday)\n      call mpbcr(irtopnight)\n"),
    ("      zbd(:,0)=SBK*zttop**4 !We could add IR flux from M dwarf host star here? --AYP\n",
     "      zbd(:,0)=SBK*zttop**4 !We could add IR flux from M dwarf host star here? --AYP\n"
     "      where (gmu0(:) > 0.)    ! Terluna: infrared from space at the top of the air\n"
     "       zbd(:,0)=zbd(:,0)+irtopday\n"
     "      elsewhere\n"
     "       zbd(:,0)=zbd(:,0)+irtopnight\n"
     "      endwhere\n"),
]
# The curved air's sunlight (radmod.f90; climate/gcm/twilight.py). PlaSim lights a column only while its Sun is above
# the horizon and stretches the direct beam's water path with Earth's horizon magnification, 35 / sqrt(1 + 1224 mu^2).
# NTWILIGHT 1 adds, at the end of swr, what the Moon's tall curved air takes beyond a flat column of the same air, by
# the solar elevation at each point (solang keeps the cosine below the horizon too): to each layer, and to the ground
# in the column's clear share, from the table twilight.dat the runner writes from ../results/gcm/twilight.json. The
# extra flux enters at the top, so the energy budget carries it. NTWILW 1 adds, at the end of lwr, the infrared the same
# geometry lets out beyond a flat column, each layer's and the ground's at the design column's temperatures scaled with
# their own to the fourth power. ZMLUNAR (a, b, c) replaces Earth's magnification with the Moon's,
# (mu + a) / (mu^2 + b mu + c). The defaults, 0, leave the model as it was.
TWILIGHT_EDITS = [
    ("      real :: irtopnight = 0.0   ! Terluna: and by night\n",
     "      real :: irtopnight = 0.0   ! Terluna: and by night\n"
     "      integer :: ntwilight = 0   ! Terluna: 1 adds the curved air's twilight and side light (twilight.dat)\n"
     "      real :: zmlunar(3) = 0.0   ! Terluna: the Moon's slant water path (mu+a)/(mu^2+b mu+c); 0 keeps Earth's\n"
     "      integer, parameter :: NTWLMAX = 400 ! Terluna: elevations the twilight table may hold\n"
     "      integer :: ntwl = 0        ! Terluna: elevations it holds\n"
     "      real :: twlel(NTWLMAX) = 0.0       ! Terluna: solar elevation, degrees\n"
     "      real :: twlg(NTWLMAX) = 0.0        ! Terluna: added to the ground, per unit sunlight at the beam\n"
     "      real :: twlq(NLEV,NTWLMAX) = 0.0   ! Terluna: added to each layer, the same\n"
     "      real :: gmuz(NHOR) = 0.0   ! Terluna: cosine of the solar zenith angle, below the horizon too\n"
     "      integer :: ntwilw = 0      ! Terluna: 1 adds the curved air's extra infrared escape (twilight.dat)\n"
     "      real :: twlw(NLEV+1) = 0.0 ! Terluna: extra infrared to space from each layer and the ground\n"
     "      real :: twlt(NLEV+1) = 0.0 ! Terluna: the temperatures it is given at\n"),
    (",minwavel,rayscale,nsunclock,swh2oscale,irtopday,irtopnight\n",
     ",minwavel,rayscale,nsunclock,swh2oscale,irtopday,irtopnight  &\n     &               ,ntwilight,zmlunar,ntwilw\n"),
    ("      call mpbcr(irtopnight)\n",
     "      call mpbcr(irtopnight)\n"
     "      call mpbci(ntwilight)\n"
     "      call mpbcrn(zmlunar,3)\n"
     "      call mpbci(ntwilw)\n"
     "      if (ntwilight == 1 .or. ntwilw == 1) then   ! Terluna: the curved air's twilight table\n"
     "       if (mypid == NROOT) then\n"
     "        open(47,file='twilight.dat',form='formatted',status='old')\n"
     "        read(47,*) ntwl\n"
     "        if (ntwl < 2 .or. ntwl > NTWLMAX) stop 'twilight.dat: too few or too many elevations'\n"
     "        do jtwl = 1, ntwl\n"
     "         read(47,*) twlel(jtwl),twlg(jtwl),(twlq(jlev,jtwl),jlev=1,NLEV)\n"
     "        enddo\n"
     "        read(47,*) (twlw(jlev),jlev=1,NLEV+1)\n"
     "        read(47,*) (twlt(jlev),jlev=1,NLEV+1)\n"
     "        close(47)\n"
     "       endif\n"
     "       call mpbci(ntwl)\n"
     "       call mpbcrn(twlel,NTWLMAX)\n"
     "       call mpbcrn(twlg,NTWLMAX)\n"
     "       call mpbcrn(twlq,NLEV*NTWLMAX)\n"
     "       call mpbcrn(twlw,NLEV+1)\n"
     "       call mpbcrn(twlt,NLEV+1)\n"
     "      endif\n"),
    ("      gmu0(:) = 0.0\n", "      gmu0(:) = 0.0\n      gmuz(:) = 0.0 ! Terluna: the cosine below the horizon too\n"),
    ("         zmuz=sin(zdecl)*sid(jlat)+cola(jlat)*cos(zdecl)*cos(zhangle)\n"
     "         if (zmuz > zdawn) gmu0(jhor) = zmuz\n",
     "         zmuz=sin(zdecl)*sid(jlat)+cola(jlat)*cos(zdecl)*cos(zhangle)\n"
     "         if (zmuz > zdawn) gmu0(jhor) = zmuz\n"
     "         gmuz(jhor) = zmuz\n"),
    ("          zmuz=solslatsdec+solclatcdec/PI\n         endif\n         if (zmuz > zdawn) gmu0(jhor) = zmuz\n",
     "          zmuz=solslatsdec+solclatcdec/PI\n         endif\n         if (zmuz > zdawn) gmu0(jhor) = zmuz\n"
     "         gmuz(jhor) = zmuz\n"),
    ("       zm(:)=35./SQRT(1.+1224.*zmu0(:)*zmu0(:))\n",
     "       zm(:)=35./SQRT(1.+1224.*zmu0(:)*zmu0(:))\n"
     "       zm(:)=merge((zmu0(:)+zmlunar(1))                                       &\n"
     "     &            /max(zmu0(:)*zmu0(:)+zmlunar(2)*zmu0(:)+zmlunar(3),1.e-12)    &\n"
     "     &            ,zm(:),zmlunar(3) > 0.) ! Terluna: the Moon's slant path\n"),
    ("      return\n      end subroutine swr\n",
     "      if (ntwilight == 1) call twilight_term ! Terluna: the curved air's twilight and side light\n"
     "      return\n      end subroutine swr\n"
     "\n"
     "      subroutine twilight_term\n"
     "      use radmod\n"
     "      implicit none\n"
     "!\n"
     "!     Terluna: the sunlight the Moon's curved air adds to the flat column (climate/gcm/twilight.py), by the\n"
     "!     solar elevation at each point, below the horizon too: to each layer and, in the column's clear share,\n"
     "!     to the ground, carried down from the top as net downward flux\n"
     "!\n"
     "      integer :: jhor, jlev, jlo, jhi, jmid\n"
     "      real :: zel, zu, zsum, zclear, zsun\n"
     "      zsun = gsol0 * gdist2\n"
     "      do jhor = 1, NHOR\n"
     "       zel = asin(max(-1.0, min(1.0, gmuz(jhor)))) * 180.0 / PI\n"
     "       if (zel < twlel(1) .or. zel > twlel(ntwl)) cycle\n"
     "       jlo = 1\n"
     "       jhi = ntwl\n"
     "       do while (jhi - jlo > 1)\n"
     "        jmid = (jlo + jhi) / 2\n"
     "        if (twlel(jmid) > zel) then\n"
     "         jhi = jmid\n"
     "        else\n"
     "         jlo = jmid\n"
     "        endif\n"
     "       enddo\n"
     "       zu = (zel - twlel(jlo)) / (twlel(jhi) - twlel(jlo))\n"
     "       zclear = 1.0\n"
     "       do jlev = 1, NLEV\n"
     "        zclear = zclear * (1.0 - dcc(jhor,jlev) * nclouds)\n"
     "       enddo\n"
     "       zsum = zsun * ((1.0 - zu) * twlg(jlo) + zu * twlg(jhi)) * zclear\n"
     "       dswfl(jhor,NLEP) = dswfl(jhor,NLEP) + zsum\n"
     "       dfd(jhor,NLEP) = dfd(jhor,NLEP) + zsum\n"
     "       do jlev = NLEV, 1, -1\n"
     "        zsum = zsum + zsun * ((1.0 - zu) * twlq(jlev,jlo) + zu * twlq(jlev,jhi))\n"
     "        dswfl(jhor,jlev) = dswfl(jhor,jlev) + zsum\n"
     "        dfd(jhor,jlev) = dfd(jhor,jlev) + zsum\n"
     "       enddo\n"
     "      enddo\n"
     "      return\n"
     "      end subroutine twilight_term\n"),
    ("      dlwfl(:,:)=dftu(:,:)+dftd(:,:)\n!\n      return\n      end subroutine lwr\n",
     "      dlwfl(:,:)=dftu(:,:)+dftd(:,:)\n"
     "      if (ntwilw == 1) call twilight_lw ! Terluna: the curved air's extra infrared escape\n"
     "!\n      return\n      end subroutine lwr\n"
     "\n"
     "      subroutine twilight_lw\n"
     "      use radmod\n"
     "      implicit none\n"
     "!\n"
     "!     Terluna: the infrared the Moon's curved air loses beyond a flat column (climate/gcm/twilight.py): each\n"
     "!     layer's and the ground's extra escape to space at the design column's temperatures, scaled with their\n"
     "!     own to the fourth power, carried up to the top as upward flux (negative in PlaSim)\n"
     "!\n"
     "      integer :: jhor, jlev\n"
     "      real :: zsum\n"
     "      do jhor = 1, NHOR\n"
     "       zsum = twlw(NLEP) * (dt(jhor,NLEP) / twlt(NLEP))**4\n"
     "       dftu(jhor,NLEP) = dftu(jhor,NLEP) - zsum\n"
     "       dlwfl(jhor,NLEP) = dlwfl(jhor,NLEP) - zsum\n"
     "       do jlev = NLEV, 1, -1\n"
     "        zsum = zsum + twlw(jlev) * (dt(jhor,jlev) / twlt(jlev))**4\n"
     "        dftu(jhor,jlev) = dftu(jhor,jlev) - zsum\n"
     "        dlwfl(jhor,jlev) = dlwfl(jhor,jlev) - zsum\n"
     "       enddo\n"
     "      enddo\n"
     "      return\n"
     "      end subroutine twilight_lw\n"),
]
# A simplified Betts-Miller convection scheme (rainmod.f90, mksbm; the author's direction of 2026-10-21, the climate
# staying on ExoPlaSim), chosen by PlaSim's own switch: nprc 1 keeps Kuo's scheme exactly, 2 calls mksbm in its place.
# Kuo's turns the column's moisture supply into rain wherever an undiluted surface parcel is buoyant, so over the
# Near-side Sea it rains 2.9 mm/day after midnight while the GCM's circulation sinks there and CM1's resolved
# convection makes almost none (climate/crm/README.md). mksbm follows Frierson (2007): a parcel from the lowest level
# gives the reference temperature up to the top of its buoyancy, and the reference vapour is rh_bm of saturation there;
# where the parcel is buoyant the temperature and vapour relax toward them over tau_bm, the rain is the vapour removed,
# and the heating is shifted by one amount through the convecting layers to release the rain's latent heat, keeping
# the column's moist enthalpy. Two lunar options, off by default: eps_bm mixes surrounding air into the rising parcel,
# a share 1 - exp(-eps_bm ln(p_below/p)) per layer, and cin_bm (J/kg) holds back deep convection below a parcel whose
# inhibition exceeds it. A buoyant column too dry to rain, held back, or topped below pdeep takes PlaSim's shallow
# scheme, as under Kuo. The rain formed at each level follows the vapour removed there, for mkrain's evaporation.
SBM_MARKER = 'Terluna: the simplified Betts-Miller scheme'
SBM_CODE = """
!     ================
!     SUBROUTINE MKSBM
!     ================

      subroutine mksbm
      use rainmod
!
!     Terluna: the simplified Betts-Miller scheme (Frierson 2007, J. Atmos. Sci. 64, 1959-1976), in place of Kuo's
!     with nprc = 2. A parcel lifted from the lowest level, dry and then moist adiabatically with PlaSim's saturation,
!     and mixing in a share 1 - exp(-eps_bm*ln(p_below/p)) of the surrounding air per layer, gives the reference
!     temperature from the lowest level to the top of its buoyancy; the reference vapour is rh_bm of saturation there.
!     Where the parcel is buoyant (and its inhibition below cin_bm J/kg, when that is set) the temperature and vapour
!     relax toward the references over tau_bm s; the rain is the vapour removed, and the heating is shifted by one
!     amount through the convecting layers so that it releases the rain's latent heat, keeping the moist enthalpy.
!
      integer :: ishallow(NHOR), ilift(NHOR), itop(NHOR)
      integer :: jhor, jlev, jtop, jlfc
      real :: zte(NHOR,NLEV), zqe(NHOR,NLEV), zdtdts(NHOR,NLEV), zdqdts(NHOR,NLEV), zdtc(NHOR,NLEV)
      real :: ztp(NLEV), zdtr(NLEV), zdqr(NLEV), zm(NLEV), zcp(NLEV), zp(NLEV)
      real :: ztnew, zqnew, zmix, zcin, zpq, zpt, zcpm, zcorr, zcond, ztau
!
      ztau=max(tau_bm,deltsec2)
      ishallow(:)=0
      ilift(:)=-1
      itop(:)=-1
      icctot(:)=0
      icclev(:,:)=0
      dprcl(:,:)=0.
      dprscl(:,:)=0.
      zdtdts(:,:)=0.
      zdqdts(:,:)=0.
      zdtc(:,:)=0.
      if(ndiaggp==1) then
       do jlev=1,NLEV
        dgp3d(:,jlev,16)=-dqdt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
        dgp3d(:,jlev,18)=dqt(:,jlev)*dp(:)*dsigma(jlev)/ga/1000.
       enddo
      end if
      zte(:,:)=dt(:,1:NLEV)+dtdt(:,1:NLEV)*deltsec2
      zqe(:,:)=dq(:,1:NLEV)+dqdt(:,1:NLEV)*deltsec2
!
      do jhor=1,NHOR
       zp(:)=dp(jhor)*sigma(:)
       zm(:)=dp(jhor)*dsigma(:)/ga
       zcp(:)=acpd*(1.+ADV*zqe(jhor,:))
!      the parcel, from the lowest level up
       ztnew=zte(jhor,NLEV)
       zqnew=zqe(jhor,NLEV)
       call sbmsat(ztnew,zqnew,zp(NLEV))
       ztp(NLEV)=ztnew
       do jlev=NLEM,1,-1
        ztnew=ztnew*(sigma(jlev)/sigma(jlev+1))**akap
        if(eps_bm > 0.) then
         zmix=1.-exp(-eps_bm*log(sigma(jlev+1)/sigma(jlev)))
         ztnew=(1.-zmix)*ztnew+zmix*zte(jhor,jlev)
         zqnew=(1.-zmix)*zqnew+zmix*zqe(jhor,jlev)
        endif
        call sbmsat(ztnew,zqnew,zp(jlev))
        ztp(jlev)=ztnew
       enddo
!      its level of free convection, the lowest above the start where it is buoyant, and its top, the last level
!      above that it stays buoyant without a break
       jlfc=0
       do jlev=NLEM,1,-1
        if(ztp(jlev) > zte(jhor,jlev)) then
         jlfc=jlev
         exit
        endif
       enddo
       if(jlfc == 0) cycle
       jtop=jlfc
       do jlev=jlfc-1,1,-1
        if(ztp(jlev) <= zte(jhor,jlev)) exit
        jtop=jlev
       enddo
       zcin=0.
       do jlev=jlfc+1,NLEV
        zcin=zcin+gascon*max(0.,zte(jhor,jlev)-ztp(jlev))*dsigma(jlev)/sigma(jlev)
       enddo
!      the relaxation toward the references, from the lowest level to the top
       zpq=0.
       zpt=0.
       zcpm=0.
       do jlev=jtop,NLEV
        zdtr(jlev)=(ztp(jlev)-zte(jhor,jlev))/ztau
        zdqr(jlev)=(rh_bm*sbmqs(ztp(jlev),zp(jlev))-zqe(jhor,jlev))/ztau
        zpq=zpq-zdqr(jlev)*zm(jlev)
        zpt=zpt+zcp(jlev)*zdtr(jlev)*zm(jlev)
        zcpm=zcpm+zcp(jlev)*zm(jlev)
       enddo
       if(zpq <= 0. .or. zp(jtop) > pdeep .or. (cin_bm > 0. .and. zcin > cin_bm)) then
        if(nshallow == 1) then
         ishallow(jhor)=1
         ilift(jhor)=NLEP
         itop(jhor)=jtop
        endif
        cycle
       endif
!      deep convection: the heating shifted to release the rain's latent heat, the rain where the vapour is removed
       zcorr=(ALV*zpq-zpt)/zcpm
       zcond=0.
       do jlev=jtop,NLEV
        zdtr(jlev)=zdtr(jlev)+zcorr
        zcond=zcond+max(0.,-zdqr(jlev))*zm(jlev)
       enddo
       do jlev=jtop,NLEV
        dtdt(jhor,jlev)=dtdt(jhor,jlev)+zdtr(jlev)
        dqdt(jhor,jlev)=dqdt(jhor,jlev)+zdqr(jlev)
        zdtc(jhor,jlev)=zdtr(jlev)
        if(zdqr(jlev) < 0.) then
         dprcl(jhor,jlev)=-zdqr(jlev)*zpq/zcond*dsigma(jlev)
         icclev(jhor,jlev)=1
         icctot(jhor)=icctot(jhor)+1
         if(sigma(jlev) < 0.4 .and. jlev == jtop) icclev(jhor,jlev)=2
        endif
       enddo
      enddo
!
      if(nshallow == 1) then
       call mkshallow(zte,zqe,ishallow,ilift,itop,zdtdts,zdqdts)
       dtdt(:,1:NLEV)=dtdt(:,1:NLEV)+zdtdts(:,1:NLEV)
       dqdt(:,1:NLEV)=dqdt(:,1:NLEV)+zdqdts(:,1:NLEV)
      endif
      if(ndiaggp==1) then
       do jlev=1,NLEV
        dgp3d(:,jlev,8)=zdtc(:,jlev)+zdtdts(:,jlev)
       enddo
      end if
      return
!
      contains
!
      real function sbmqs(pt,pp)
!     saturation specific humidity, as Kuo's parcel takes it
      real, intent(in) :: pt,pp
      real :: zqs
      zqs=rdbrv*ra1*exp(ra2*(pt-TMELT)/(pt-ra4))/pp
      zqs=min(zqs,rdbrv)
      sbmqs=zqs/(1.-(1./rdbrv-1.)*zqs)
      end function sbmqs
!
      subroutine sbmsat(pt,pq,pp)
!     a supersaturated parcel condensed to saturation in two steps, as Kuo's, the ice's latent heat below freezing
      real, intent(inout) :: pt,pq
      real, intent(in) :: pp
      real :: zqs,zcor,zlcp,zdqsdt,zrain
      integer :: jit
      do jit=1,2
       zqs=rdbrv*ra1*exp(ra2*(pt-TMELT)/(pt-ra4))/pp
       zqs=min(zqs,rdbrv)
       zcor=1./(1.-(1./rdbrv-1.)*zqs)
       zqs=zqs*zcor
       if(zqs >= pq) return
       if(pt < TMELT) then
        zlcp=ALS/(acpd*(1.+ADV*pq))
       else
        zlcp=ALV/(acpd*(1.+ADV*pq))
       endif
       zdqsdt=ra2*(TMELT-ra4)*zqs*zcor/((pt-ra4)*(pt-ra4))
       zrain=(pq-zqs)/(1.+zlcp*zdqsdt)
       pt=pt+zlcp*zrain
       pq=pq-zrain
      enddo
      end subroutine sbmsat
!
      end subroutine mksbm
"""
SBM_EDITS = [
    ("      real :: gamma   = 0.01  ! tuning parameter for evaporation of precip.\n",
     "      real :: gamma   = 0.01  ! tuning parameter for evaporation of precip.\n"
     "      real :: tau_bm  = 7200. ! Terluna: the simplified Betts-Miller scheme's relaxation time (s), nprc = 2\n"
     "      real :: rh_bm   = 0.7   ! its reference relative humidity\n"
     "      real :: eps_bm  = 0.    ! its parcel's mixing per e-folding of pressure (0: undiluted)\n"
     "      real :: cin_bm  = 0.    ! its limit on the parcel's inhibition (J/kg; 0: none)\n"),
    ("      namelist/rainmod_nl/kbeta,nprl,nprc,ndca,ncsurf,nmoment,nshallow  &\n",
     "      namelist/rainmod_nl/tau_bm,rh_bm,eps_bm,cin_bm   ! Terluna: the simplified Betts-Miller scheme\n"
     "      namelist/rainmod_nl/kbeta,nprl,nprc,ndca,ncsurf,nmoment,nshallow  &\n"),
    ("      call mpbcr(gamma)\n",
     "      call mpbcr(gamma)\n      call mpbcr(tau_bm)\n      call mpbcr(rh_bm)\n      call mpbcr(eps_bm)\n"
     "      call mpbcr(cin_bm)\n"),
    ("      if(nprc == 1) call mkdqtgp\n", "      if(nprc >= 1) call mkdqtgp\n"),
    ("      if(nprc == 1) call kuo\n",
     "      if(nprc == 1) call kuo\n      if(nprc == 2) call mksbm   ! Terluna: the simplified Betts-Miller scheme\n"),
    ("      end subroutine kuo\n", "      end subroutine kuo\n" + SBM_CODE),
]
# (file, marker, edits) in the order applied; a marker appears only in the text its patch adds.
PATCHES = [('radmod.f90', 'Terluna: calibration multiplier', RADMOD_EDITS),
           ('rainmod.f90', 'Terluna: gravity for the heights below', RAINMOD_EDITS),
           ('rainmod.f90', 'Terluna: seconds per day in the convective cloud formula', CONVECTIVE_EDITS),
           ('radmod.f90', "Terluna: the Sun's clock", SUN_EDITS),
           ('radmod.f90', 'Terluna: calibrated absorptivity', SHORTWAVE_EDITS),
           ('radmod.f90', 'Terluna: infrared from space at the top of the air', TOP_INFRARED_EDITS),
           ('rainmod.f90', SBM_MARKER, SBM_EDITS)]
# A CAPE trigger for the simplified Betts-Miller scheme (24 October), off by default: with cape_bm > 0 a column convects
# deep only where a second parcel, lifted as the reference's but mixing epst_bm per e-folding of pressure, has at least
# cape_bm J/kg of CAPE from its own level of free convection up through the layers it stays buoyant, counted as the
# scheme counts its inhibition; a column held back takes the shallow scheme. It is a patch of its own on the scheme's
# text, applied only for a run that sets it, so a source the scheme's patch has already reached takes it in place and
# runs without it keep their source; with cape_bm 0 the scheme's arithmetic is the same.
SBM_TRIGGER_MARKER = 'Terluna: the CAPE trigger'
SBM_TRIGGER_BLOCK = """\
!      Terluna: the CAPE trigger: deep convection only where a parcel mixing epst_bm per e-folding of pressure has
!      at least cape_bm J/kg of CAPE, from its own level of free convection up through the layers it stays buoyant
       zcapet=0.
       if(cape_bm > 0.) then
        ztnew=zte(jhor,NLEV)
        zqnew=zqe(jhor,NLEV)
        call sbmsat(ztnew,zqnew,zp(NLEV))
        ztt(NLEV)=ztnew
        do jlev=NLEM,1,-1
         ztnew=ztnew*(sigma(jlev)/sigma(jlev+1))**akap
         if(epst_bm > 0.) then
          zmix=1.-exp(-epst_bm*log(sigma(jlev+1)/sigma(jlev)))
          ztnew=(1.-zmix)*ztnew+zmix*zte(jhor,jlev)
          zqnew=(1.-zmix)*zqnew+zmix*zqe(jhor,jlev)
         endif
         call sbmsat(ztnew,zqnew,zp(jlev))
         ztt(jlev)=ztnew
        enddo
        jlfct=0
        do jlev=NLEM,1,-1
         if(ztt(jlev) > zte(jhor,jlev)) then
          jlfct=jlev
          exit
         endif
        enddo
        do jlev=jlfct,1,-1
         if(ztt(jlev) <= zte(jhor,jlev)) exit
         zcapet=zcapet+gascon*(ztt(jlev)-zte(jhor,jlev))*dsigma(jlev)/sigma(jlev)
        enddo
       endif
"""
SBM_TRIGGER_EDITS = [
    ("      real :: cin_bm  = 0.    ! its limit on the parcel's inhibition (J/kg; 0: none)\n",
     "      real :: cin_bm  = 0.    ! its limit on the parcel's inhibition (J/kg; 0: none)\n"
     "      real :: cape_bm = 0.    ! Terluna: the CAPE trigger: deep only above this CAPE (J/kg; 0: none)\n"
     "      real :: epst_bm = 0.    ! its parcel's mixing per e-folding of pressure\n"),
    ("      namelist/rainmod_nl/tau_bm,rh_bm,eps_bm,cin_bm   ! Terluna: the simplified Betts-Miller scheme\n",
     "      namelist/rainmod_nl/tau_bm,rh_bm,eps_bm,cin_bm   ! Terluna: the simplified Betts-Miller scheme\n"
     "      namelist/rainmod_nl/cape_bm,epst_bm   ! Terluna: the CAPE trigger\n"),
    ("      call mpbcr(cin_bm)\n", "      call mpbcr(cin_bm)\n      call mpbcr(cape_bm)\n      call mpbcr(epst_bm)\n"),
    ("      real :: ztnew, zqnew, zmix, zcin, zpq, zpt, zcpm, zcorr, zcond, ztau\n",
     "      real :: ztnew, zqnew, zmix, zcin, zpq, zpt, zcpm, zcorr, zcond, ztau\n"
     "      real :: ztt(NLEV), zcapet   ! Terluna: the CAPE trigger\n      integer :: jlfct\n"),
    ("       if(zpq <= 0. .or. zp(jtop) > pdeep .or. (cin_bm > 0. .and. zcin > cin_bm)) then\n",
     SBM_TRIGGER_BLOCK + "       if(zpq <= 0. .or. zp(jtop) > pdeep .or. (cin_bm > 0. .and. zcin > cin_bm) .or. "
     "(cape_bm > 0. .and. zcapet < cape_bm)) then\n"),
]
SBM_TRIGGER_PATCHES = [('rainmod.f90', SBM_TRIGGER_MARKER, SBM_TRIGGER_EDITS)]
# Two more options for the simplified Betts-Miller scheme (25 October), each off by default. The cloud-depth test
# (dpdeep_bm > 0, Pa) lets a column convect deep only where the parcel's cloud, from the condensation level of the
# lowest level's air (found by bisection in log pressure with PlaSim's saturation) up to the scheme's top level, is
# more than dpdeep_bm thick; a thinner cloud takes PlaSim's shallow scheme, which does not rain. It is the IFS's rule
# (documentation Cy36r1, Part IV, section 6.4: a cloud from its LCL up thicker than 200 hPa is deep), where PlaSim's
# own pdeep tests the top's pressure alone. The mixing within layers (nsub_bm > 1, with eps_bm > 0) carries the
# parcel across each layer in nsub_bm equal steps of log pressure, each lifting it dry adiabatically, mixing in its
# step's share of the surrounding air interpolated there in log pressure, and condensing the mixture, so the mixing
# acts through the thick layers rather than once at their tops; one step is the scheme as before, and eps_bm 0 never
# takes the steps. Like the CAPE trigger it is a patch of its own on the scheme's text, applied only for a run that
# sets either option, on a source with or without the trigger; with both off the routine's arithmetic is unchanged.
SBM_DEPTH_MARKER = 'Terluna: the cloud-depth test and the mixing within layers'
SBM_DEPTH_LOOP_OLD = """\
       ztp(NLEV)=ztnew
       do jlev=NLEM,1,-1
        ztnew=ztnew*(sigma(jlev)/sigma(jlev+1))**akap
        if(eps_bm > 0.) then
         zmix=1.-exp(-eps_bm*log(sigma(jlev+1)/sigma(jlev)))
         ztnew=(1.-zmix)*ztnew+zmix*zte(jhor,jlev)
         zqnew=(1.-zmix)*zqnew+zmix*zqe(jhor,jlev)
        endif
        call sbmsat(ztnew,zqnew,zp(jlev))
        ztp(jlev)=ztnew
       enddo
"""
SBM_DEPTH_LOOP_NEW = """\
       ztp(NLEV)=ztnew
       do jlev=NLEM,1,-1
        if(eps_bm > 0. .and. nsub_bm > 1) then
!        Terluna: the mixing within layers: nsub_bm equal steps of log pressure, each lifting the parcel dry
!        adiabatically, mixing in its share of the air interpolated there in log pressure and condensing the mixture
         zsub=log(sigma(jlev+1)/sigma(jlev))/real(nsub_bm)
         zmix=1.-exp(-eps_bm*zsub)
         zstep=exp(-zsub)
         zlift=exp(-akap*zsub)
         zsig0=sigma(jlev+1)
         do jsub=1,nsub_bm
          if(jsub < nsub_bm) then
           zsig1=zsig0*zstep
           ztnew=ztnew*zlift
          else
           zsig1=sigma(jlev)
           ztnew=ztnew*(zsig1/zsig0)**akap
          endif
          zw=real(jsub)/real(nsub_bm)
          ztnew=(1.-zmix)*ztnew+zmix*((1.-zw)*zte(jhor,jlev+1)+zw*zte(jhor,jlev))
          zqnew=(1.-zmix)*zqnew+zmix*((1.-zw)*zqe(jhor,jlev+1)+zw*zqe(jhor,jlev))
          call sbmsat(ztnew,zqnew,dp(jhor)*zsig1)
          zsig0=zsig1
         enddo
        else
         ztnew=ztnew*(sigma(jlev)/sigma(jlev+1))**akap
         if(eps_bm > 0.) then
          zmix=1.-exp(-eps_bm*log(sigma(jlev+1)/sigma(jlev)))
          ztnew=(1.-zmix)*ztnew+zmix*zte(jhor,jlev)
          zqnew=(1.-zmix)*zqnew+zmix*zqe(jhor,jlev)
         endif
         call sbmsat(ztnew,zqnew,zp(jlev))
        endif
        ztp(jlev)=ztnew
       enddo
"""
SBM_DEPTH_TEST = """\
!      Terluna: the cloud-depth test: deep convection only where the parcel's cloud, from the condensation level of
!      the lowest level's air up to the top, is more than dpdeep_bm thick; a thinner cloud takes the shallow scheme
       if(dpdeep_bm > 0.) then
        zpbase=sbmlcl(zte(jhor,NLEV),zqe(jhor,NLEV),zp(NLEV),zp(jtop))
        if(zpbase-zp(jtop) <= dpdeep_bm) then
         if(nshallow == 1) then
          ishallow(jhor)=1
          ilift(jhor)=NLEP
          itop(jhor)=jtop
         endif
         cycle
        endif
       endif
"""
SBM_DEPTH_LCL = """\
!
      real function sbmlcl(pt,pq,pp,ptop)
!     Terluna: the pressure at which air of temperature pt and vapour pq at pp, lifted dry adiabatically with its
!     vapour kept, reaches saturation: pp if it is saturated there, ptop if it is still unsaturated at ptop, and
!     otherwise found by bisection in log pressure between the two
      real, intent(in) :: pt,pq,pp,ptop
      real :: zlo,zhi,zmid,zpm
      integer :: jit
      if(pq >= sbmqs(pt,pp)) then
       sbmlcl=pp
       return
      endif
      if(pq < sbmqs(pt*(ptop/pp)**akap,ptop)) then
       sbmlcl=ptop
       return
      endif
      zlo=log(ptop)
      zhi=log(pp)
      do jit=1,24
       zmid=0.5*(zlo+zhi)
       zpm=exp(zmid)
       if(pq >= sbmqs(pt*(zpm/pp)**akap,zpm)) then
        zlo=zmid
       else
        zhi=zmid
       endif
      enddo
      sbmlcl=exp(0.5*(zlo+zhi))
      end function sbmlcl
"""
SBM_DEPTH_EDITS = [
    ("      real :: cin_bm  = 0.    ! its limit on the parcel's inhibition (J/kg; 0: none)\n",
     "      real :: cin_bm  = 0.    ! its limit on the parcel's inhibition (J/kg; 0: none)\n"
     "      real :: dpdeep_bm = 0.  ! Terluna: the cloud-depth test: deep only above this depth (Pa; 0: no test)\n"
     "      integer :: nsub_bm = 1  ! the parcel's mixing steps within each layer (1: one a layer)\n"),
    ("      namelist/rainmod_nl/tau_bm,rh_bm,eps_bm,cin_bm   ! Terluna: the simplified Betts-Miller scheme\n",
     "      namelist/rainmod_nl/tau_bm,rh_bm,eps_bm,cin_bm   ! Terluna: the simplified Betts-Miller scheme\n"
     f"      namelist/rainmod_nl/dpdeep_bm,nsub_bm   ! {SBM_DEPTH_MARKER}\n"),
    ("      call mpbcr(cin_bm)\n", "      call mpbcr(cin_bm)\n      call mpbcr(dpdeep_bm)\n      call mpbci(nsub_bm)\n"),
    ("      real :: ztnew, zqnew, zmix, zcin, zpq, zpt, zcpm, zcorr, zcond, ztau\n",
     "      real :: ztnew, zqnew, zmix, zcin, zpq, zpt, zcpm, zcorr, zcond, ztau\n"
     "      real :: zsub, zstep, zlift, zsig0, zsig1, zw, zpbase   ! Terluna: the cloud-depth test and the mixing within "
     "layers\n"
     "      integer :: jsub\n"),
    (SBM_DEPTH_LOOP_OLD, SBM_DEPTH_LOOP_NEW),
    ("!      the relaxation toward the references, from the lowest level to the top\n",
     SBM_DEPTH_TEST + "!      the relaxation toward the references, from the lowest level to the top\n"),
    ("      end function sbmqs\n", "      end function sbmqs\n" + SBM_DEPTH_LCL),
]
SBM_DEPTH_PATCHES = [('rainmod.f90', SBM_DEPTH_MARKER, SBM_DEPTH_EDITS)]
# The shallow cap (26 October), off by default and applied only for a run that sets it, after the depth patch it needs.
# PlaSim's shallow scheme (mkshallow) holds every shallow cloud's top at or below the lowest level above pdeepth, 700 hPa,
# Earth's lower troposphere and about 28 km up in the Moon's air. With nshcap_bm 1 and the depth test set, a column the
# scheme hands to the shallow scheme after the depth test (too dry to rain, or held back by its inhibition or the CAPE
# trigger) gets a top no higher than dpdeep_bm above its cloud's base, the highest level at or below its parcel's top
# within that depth of the condensation level, and mkshallow's absolute cap gives way to it; a cloud the depth test
# itself hands over is already that thin. With nshcap_bm 0 the routine's arithmetic and mkshallow's are as before.
SBM_CAP_MARKER = 'Terluna: the shallow cap'
SBM_CAP_EDITS = [
    ("      integer :: nsub_bm = 1  ! the parcel's mixing steps within each layer (1: one a layer)\n",
     "      integer :: nsub_bm = 1  ! the parcel's mixing steps within each layer (1: one a layer)\n"
     f"      integer :: nshcap_bm = 0 ! {SBM_CAP_MARKER}: 1 holds a shallow top within dpdeep_bm of its base\n"),
    (f"      namelist/rainmod_nl/dpdeep_bm,nsub_bm   ! {SBM_DEPTH_MARKER}\n",
     f"      namelist/rainmod_nl/dpdeep_bm,nsub_bm   ! {SBM_DEPTH_MARKER}\n"
     f"      namelist/rainmod_nl/nshcap_bm   ! {SBM_CAP_MARKER}\n"),
    ("      call mpbci(nsub_bm)\n", "      call mpbci(nsub_bm)\n      call mpbci(nshcap_bm)\n"),
    ("        if(nshallow == 1) then\n         ishallow(jhor)=1\n         ilift(jhor)=NLEP\n         itop(jhor)=jtop\n"
     "        endif\n        cycle\n       endif\n"
     "!      deep convection: the heating shifted to release the rain's latent heat, the rain where the vapour is removed\n",
     "        if(nshallow == 1) then\n         ishallow(jhor)=1\n         ilift(jhor)=NLEP\n         itop(jhor)=jtop\n"
     f"!        {SBM_CAP_MARKER}: the top no higher than dpdeep_bm above the cloud's base\n"
     "         if(nshcap_bm == 1 .and. dpdeep_bm > 0.) then\n"
     "          do jlev=jtop,NLEV\n"
     "           if(zp(jlev) >= zpbase-dpdeep_bm) exit\n"
     "          enddo\n"
     "          itop(jhor)=min(jlev,NLEV)\n"
     "         endif\n"
     "        endif\n        cycle\n       endif\n"
     "!      deep convection: the heating shifted to release the rain's latent heat, the rain where the vapour is removed\n"),
    ("      ktop(:)=MAX(ktop(:),kktop(:))\n",
     f"      if(nprc /= 2 .or. nshcap_bm /= 1) ktop(:)=MAX(ktop(:),kktop(:)) ! {SBM_CAP_MARKER} replaces it\n"),
]
SBM_CAP_PATCHES = [('rainmod.f90', SBM_CAP_MARKER, SBM_CAP_EDITS)]
# The night switch (27 October), a diagnostic and never a configuration, off by default and applied only for a run that
# sets it. With ndiagnight_bm 1 the simplified Betts-Miller scheme convects deep nowhere the Sun is down: where the
# geometric cosine of the Sun's zenith at the grid point is at most 0 (radmod's gmuz, the twilight patch's cosine that
# solang sets each step before the convection runs, below the horizon too), a column that would convect deep takes
# PlaSim's shallow scheme as a column held back by the inhibition limit or the CAPE trigger takes it, its top held by
# the shallow cap where that is set. The scheme reads the cosine through sbmnight, a routine added to radmod.f90 and
# called without a module dependency, so rainmod.f90 still compiles on PlaSim's module alone. The patch edits the text
# of the CAPE trigger's, the depth test's and the cap's patches and goes over them; with ndiagnight_bm 0 the routine's
# arithmetic is as before.
SBM_NIGHT_MARKER = 'Terluna: the night switch, a diagnostic'
SBM_NIGHT_GETTER = f"""
!     ===================
!     SUBROUTINE SBMNIGHT
!     ===================

      subroutine sbmnight(pmuz)
      use radmod
!
!     {SBM_NIGHT_MARKER} (rainmod's ndiagnight_bm): the geometric cosine of the Sun's zenith at each grid
!     point, below the horizon too, as solang set it this step
!
      real, intent(out) :: pmuz(NHOR)
      pmuz(:)=gmuz(:)
      return
      end subroutine sbmnight
"""
SBM_NIGHT_EDITS = [
    (f"      integer :: nshcap_bm = 0 ! {SBM_CAP_MARKER}: 1 holds a shallow top within dpdeep_bm of its base\n",
     f"      integer :: nshcap_bm = 0 ! {SBM_CAP_MARKER}: 1 holds a shallow top within dpdeep_bm of its base\n"
     f"      integer :: ndiagnight_bm = 0 ! {SBM_NIGHT_MARKER}: 1 holds deep convection off where the Sun is down\n"),
    (f"      namelist/rainmod_nl/nshcap_bm   ! {SBM_CAP_MARKER}\n",
     f"      namelist/rainmod_nl/nshcap_bm   ! {SBM_CAP_MARKER}\n"
     f"      namelist/rainmod_nl/ndiagnight_bm   ! {SBM_NIGHT_MARKER}\n"),
    ("      call mpbci(nshcap_bm)\n", "      call mpbci(nshcap_bm)\n      call mpbci(ndiagnight_bm)\n"),
    ("      real :: ztt(NLEV), zcapet   ! Terluna: the CAPE trigger\n",
     "      real :: ztt(NLEV), zcapet   ! Terluna: the CAPE trigger\n"
     f"      real :: zmuz(NHOR)   ! {SBM_NIGHT_MARKER}: the geometric cosine of the Sun's zenith\n"),
    ("      ztau=max(tau_bm,deltsec2)\n",
     "      ztau=max(tau_bm,deltsec2)\n      zmuz(:)=1.\n      if(ndiagnight_bm == 1) call sbmnight(zmuz)\n"),
    ("       if(zpq <= 0. .or. zp(jtop) > pdeep .or. (cin_bm > 0. .and. zcin > cin_bm) .or. "
     "(cape_bm > 0. .and. zcapet < cape_bm)) then\n",
     f"!      {SBM_NIGHT_MARKER}: no deep convection where the Sun is down, the shallow scheme instead\n"
     "       if(zpq <= 0. .or. zp(jtop) > pdeep .or. (cin_bm > 0. .and. zcin > cin_bm) .or. "
     "(cape_bm > 0. .and. zcapet < cape_bm) &\n"
     "     &    .or. (ndiagnight_bm == 1 .and. zmuz(jhor) <= 0.)) then\n"),
]
SBM_NIGHT_RADMOD_EDITS = [("      end subroutine solang\n", "      end subroutine solang\n" + SBM_NIGHT_GETTER)]
SBM_NIGHT_PATCHES = [('rainmod.f90', SBM_NIGHT_MARKER, SBM_NIGHT_EDITS),
                     ('radmod.f90', SBM_NIGHT_MARKER, SBM_NIGHT_RADMOD_EDITS)]


# Two switches in PlaSim's own namelist, read before the restart, so that the restart can see them (both 0 by default,
# when the model and its output are as before): ncwat 1 holds large-scale condensate as cloud water
# (condensate.CLOUD_WATER_PATCHES), the restart carrying it as dqcw and a restart without it giving none
# (restart_has_array); nhusmean 1 writes the humidity's true 3-day mean as code 233 (HUS_MEAN_PATCHES), the gridpoint
# humidity PlaSim adds every step (aadq), where code 133 holds the humidity at each interval's start (year_files.py).
SWITCH_MARKER = 'Terluna: switches'
SWITCH_PLASIMMOD_EDITS = [
    ("      integer :: nqspec   = 1   ! 1: spectral q   0: gridpoint q (semi-Langrangian)\n",
     "      integer :: nqspec   = 1   ! 1: spectral q   0: gridpoint q (semi-Langrangian)\n"
     f"      integer :: ncwat    = 0   ! {SWITCH_MARKER}: 1 holds large-scale condensate as cloud water\n"
     "      integer :: nhusmean = 0   ! 1 writes the humidity's true mean as code 233 (low-I/O mode)\n"),
    ("      real :: dql(NHOR,NLEP)  = 0.     ! Liquid water content\n",
     "      real :: dql(NHOR,NLEP)  = 0.     ! Liquid water content\n"
     f"      real :: dqcw(NHOR,NLEV) = 0.     ! {SWITCH_MARKER}: the cloud water mklsp holds (kg/kg)\n"),
]
SWITCH_PLASIM_EDITS = [
    ("                   , sigh    , nenergy , nener3d , nsponge , dampsp     &\n                   , l_aero\n",
     "                   , sigh    , nenergy , nener3d , nsponge , dampsp     &\n"
     "                   , l_aero  , ncwat   , nhusmean\n"),
    ("      call mpbci(ndiaggp3d ) ! no of 3d gp diagnostic arrays\n",
     "      call mpbci(ndiaggp3d ) ! no of 3d gp diagnostic arrays\n"
     f"      call mpbci(ncwat    ) ! {SWITCH_MARKER}: cloud water\n"
     "      call mpbci(nhusmean ) ! the humidity's true mean\n"),
    ("      call mpputgp('dqsat'  ,dqsat  ,NHOR,NLEV)\n",
     "      call mpputgp('dqsat'  ,dqsat  ,NHOR,NLEV)\n"
     f"      if (ncwat == 1) call mpputgp('dqcw'   ,dqcw   ,NHOR,NLEV)   ! {SWITCH_MARKER}\n"),
    ("      call mpgetgp('dqsat'  ,dqsat  ,NHOR,NLEV)\n",
     "      call mpgetgp('dqsat'  ,dqsat  ,NHOR,NLEV)\n"
     f"      if (ncwat == 1) then   ! {SWITCH_MARKER}: the cloud water, none in a restart without it\n"
     "         icwread = 0\n"
     "         if (mypid == NROOT) call restart_has_array('dqcw',icwread)\n"
     "         call mpbci(icwread)\n"
     "         if (icwread == 1) call mpgetgp('dqcw'   ,dqcw   ,NHOR,NLEV)\n"
     "      endif\n"),
]
SWITCH_RESTART_EDITS = [
    ("      end subroutine get_restart_array\n",
     "      end subroutine get_restart_array\n"
     f"\n!     {SWITCH_MARKER}: whether the restart holds the array yn (kfound 1) or not (0)\n"
     "      subroutine restart_has_array(yn,kfound)\n"
     "      use restartmod\n"
     "      character (len=*) :: yn\n"
     "      integer :: kfound\n"
     "      kfound = 0\n"
     "      do j = 1 , nresnum\n"
     "         if (trim(yn) == trim(yresnam(j))) kfound = 1\n"
     "      enddo\n"
     "      return\n"
     "      end subroutine restart_has_array\n"),
]
PLASIM_SWITCH_PATCHES = [('plasimmod.f90', SWITCH_MARKER, SWITCH_PLASIMMOD_EDITS),
                         ('plasim.f90', SWITCH_MARKER, SWITCH_PLASIM_EDITS),
                         ('restartmod.f90', SWITCH_MARKER, SWITCH_RESTART_EDITS)]
HUS_MEAN_MARKER = "Terluna: the humidity's true mean"
HUS_MEAN_CODE = 233
HUS_MEAN = ('hus_mean', 'specific_humidity_mean_over_the_output_interval', '1')
_HUS_LOWIO = ("        if (nqspec == 0) then ! Semi Langrangian advection active\n"
              "           do jlev = 1 , NLEV\n"
              "              aadq(:,jlev) = aadq(:,jlev)/real(naccuout)\n"
              "              call writegp(40,aadq(1,jlev),133,jlev)\n"
              "           enddo\n"
              "        endif\n")
HUS_MEAN_EDITS = [
    (_HUS_LOWIO,
     _HUS_LOWIO
     + f"        if (nqspec == 1 .and. nhusmean == 1) then ! {HUS_MEAN_MARKER}, code 233\n"
     "           ihuscode = 233   ! through a variable: PlaSim writes it only with nhusmean 1\n"
     "           do jlev = 1 , NLEV\n"
     "              aadq(:,jlev) = aadq(:,jlev)/real(naccuout)\n"
     "              call writegp(40,aadq(1,jlev),ihuscode,jlev)\n"
     "           enddo\n"
     "        endif\n"),
]
HUS_MEAN_PATCHES = [('outmod.f90', HUS_MEAN_MARKER, HUS_MEAN_EDITS)]
# CO2 and O2 in PlaSim's sunlight (radmod.f90, swr; climate/gcm/plasim_shortwave.py, GAS_BANDS). PlaSim's two bands hold
# no CO2 or O2, whose strongest lines the Moon's deep column takes high up, so it puts that sunlight where its water
# is. NGASBANDS 1 adds CO2 and O2 in the near-infrared band and O2 in the visible, each an absorptivity of its own
# pressure-scaled path (ZGASQ, fitted in ../results/gcm/shortwave_gases.json) and the gases' mass fractions (ZGASW),
# along the Moon's slant path for well-mixed air (ZMAIRLUNAR), multiplied into the layers' transmissions as water's
# is. With the fractions at 0, the default, every transmission is multiplied by 1 and the model is as it was.
def _gas_absorptivity(g, u):
    return f"zgasq(1,{g})*{u}/((1.+zgasq(2,{g})*{u})**zgasq(3,{g})+zgasq(4,{g})*{u})"


def _gas_layer(g, down, up):
    a = lambda u: _gas_absorptivity(g, u)
    return (f"        zgp(:)=zgyl(:,jlev,{g})\n"
            f"        ztg(:)=(1.-{a('zgp(:)')})/ztgt(:,{g})\n"
            f"        ztgt(:,{g})=ztgt(:,{g})*ztg(:)\n"
            f"        zgp(:)=zgyt(:,{g})+zbetta*(zgut(:,{g})-zgul(:,jlev,{g}))\n"
            f"        ztgu(:)=ztgtu(:,{g})/(1.-{a('zgp(:)')})\n"
            f"        ztgtu(:,{g})=ztgtu(:,{g})/ztgu(:)\n"
            f"        {down}(:)={down}(:)*ztg(:)\n"
            f"        {up}(:)={up}(:)*ztgu(:)\n")


GAS_EDITS = [
    ("      real :: twlt(NLEV+1) = 0.0 ! Terluna: the temperatures it is given at\n",
     "      real :: twlt(NLEV+1) = 0.0 ! Terluna: the temperatures it is given at\n"
     "      integer :: ngasbands = 0   ! Terluna: 1 adds CO2 and O2 to the sunlight's absorbers\n"
     "      real :: zgasq(5,3) = 0.0   ! Terluna: their fits (a, b, c, d, n): CO2 and O2 near infrared, O2 visible\n"
     "      real :: zgasw(2) = 0.0     ! Terluna: the mass fractions of CO2 and O2\n"
     "      real :: zmairlunar(3) = 0.0 ! Terluna: the Moon's slant factor for well-mixed air (0: the water's)\n"),
    ("     &               ,ntwilight,zmlunar,ntwilw\n",
     "     &               ,ntwilight,zmlunar,ntwilw,ngasbands,zgasq,zgasw,zmairlunar\n"),
    ("      call mpbci(ntwilw)\n",
     "      call mpbci(ntwilw)\n"
     "      call mpbci(ngasbands)   ! Terluna: CO2 and O2 in the sunlight\n"
     "      call mpbcrn(zgasq,15)\n"
     "      call mpbcrn(zgasw,2)\n"
     "      call mpbcrn(zmairlunar,3)\n"),
    ("      real ztwvt(NHOR),ztwvtu(NHOR)! total water vapor transmissivities (d/u)\n",
     "      real ztwvt(NHOR),ztwvtu(NHOR)! total water vapor transmissivities (d/u)\n"
     "      real zgu(NHOR,3),zgut(NHOR,3),zgyt(NHOR,3)  ! Terluna: CO2 and O2 paths: layer, total, effective total\n"
     "      real zgul(NHOR,NLEV,3),zgyl(NHOR,NLEV,3)   ! Terluna: their amounts above each layer's bottom\n"
     "      real ztgt(NHOR,3),ztgtu(NHOR,3),ztg(NHOR),ztgu(NHOR),zgp(NHOR),zmg(NHOR)\n"),
    ("       zwvt(:)=0.\n       zywvt(:)=0.\n",
     "       zwvt(:)=0.\n       zywvt(:)=0.\n"
     "       zmg(:)=merge((zmu0(:)+zmairlunar(1))                                   &\n"
     "     &             /max(zmu0(:)*zmu0(:)+zmairlunar(2)*zmu0(:)+zmairlunar(3),1.e-12) &\n"
     "     &             ,zm(:),zmairlunar(3) > 0.) ! Terluna: the slant path for CO2 and O2\n"
     + ''.join(f"       zgut(:,{g})=0.\n       zgyt(:,{g})=0.\n" for g in (1, 2, 3))),
    ("        zywvl(:,jlev)=zywvt(:)\n",
     "        zywvl(:,jlev)=zywvt(:)\n"
     + ''.join(f"        zgu(:,{g})=dsigma(jlev)*dp(:)/ga*zgasw({1 if g == 1 else 2})*(sigma(jlev)*dp(:)/1.e5)**zgasq(5,{g})\n"
               f"        zgut(:,{g})=zgut(:,{g})+zgu(:,{g})\n"
               f"        zgyt(:,{g})=zcs(:)*(zgyt(:,{g})+zmg(:)*zgu(:,{g}))+(1.-zcs(:))*(zgyt(:,{g})+zbetta*zgu(:,{g}))\n"
               f"        zgul(:,jlev,{g})=zgut(:,{g})\n"
               f"        zgyl(:,jlev,{g})=zgyt(:,{g})\n" for g in (1, 2, 3))),
    ("     &            /zsolar2) ! Terluna: calibrated absorptivity\n",
     "     &            /zsolar2) ! Terluna: calibrated absorptivity\n"
     + ''.join(f"       ztgt(:,{g})=1.\n"
               f"       zgp(:)=zgyt(:,{g})+zbetta*zgut(:,{g})\n"
               f"       ztgtu(:,{g})=1.-{_gas_absorptivity(g, 'zgp(:)')}\n" for g in (1, 2, 3))),
    ("        zto3tu(:)=zto3tu(:)/zto3u(:)\n",
     "        zto3tu(:)=zto3tu(:)/zto3u(:)\n" + _gas_layer(3, 'zto3', 'zto3u')),
    ("       ztwvtu(:)=ztwvtu(:)/ztwvu(:)\n",
     "       ztwvtu(:)=ztwvtu(:)/ztwvu(:)\n" + _gas_layer(1, 'ztwv', 'ztwvu') + _gas_layer(2, 'ztwv', 'ztwvu')),
]
# The longwave by height (2026-10-17): Sasamori's broadband absorptivities put PlaSim's longwave cooling too low in the
# Moon's moist, massive air, by up to 12 W/m2 a layer against the correlated-k model, which its own continuum and path
# exponents cannot move (plasim_longwave.py). NLWCORR 1 adds, after the curved air's escape, each layer's extra escape
# to space fitted over land and over sea, scaled with its own temperature to the fourth power, from lwcorr.dat.
LW_CORRECTION_EDITS = [
    ("      real :: twlt(NLEV+1) = 0.0 ! Terluna: the temperatures it is given at\n",
     "      real :: twlt(NLEV+1) = 0.0 ! Terluna: the temperatures it is given at\n"
     "      integer :: nlwcorr = 0     ! Terluna: 1 adds the longwave by height (lwcorr.dat)\n"
     "      real :: lwcl(NLEV+1) = 0.0 ! Terluna: the longwave by height, extra escape of each layer over land\n"
     "      real :: lwcs(NLEV+1) = 0.0 ! Terluna: the same over sea\n"
     "      real :: lwct(NLEV+1) = 0.0 ! Terluna: the temperatures they are given at\n"),
    ("     &               ,ntwilight,zmlunar,ntwilw,ngasbands,zgasq,zgasw,zmairlunar\n",
     "     &               ,ntwilight,zmlunar,ntwilw,ngasbands,zgasq,zgasw,zmairlunar,nlwcorr\n"),
    ("      call mpbci(ntwilw)\n",
     "      call mpbci(ntwilw)\n"
     "      call mpbci(nlwcorr)\n"
     "      if (nlwcorr == 1) then   ! Terluna: the longwave by height\n"
     "       if (mypid == NROOT) then\n"
     "        open(48,file='lwcorr.dat',form='formatted',status='old')\n"
     "        read(48,*) (lwcl(jlev),jlev=1,NLEV+1)\n"
     "        read(48,*) (lwcs(jlev),jlev=1,NLEV+1)\n"
     "        read(48,*) (lwct(jlev),jlev=1,NLEV+1)\n"
     "        close(48)\n"
     "       endif\n"
     "       call mpbcrn(lwcl,NLEV+1)\n"
     "       call mpbcrn(lwcs,NLEV+1)\n"
     "       call mpbcrn(lwct,NLEV+1)\n"
     "      endif\n"),
    ("      if (ntwilw == 1) call twilight_lw ! Terluna: the curved air's extra infrared escape\n",
     "      if (ntwilw == 1) call twilight_lw ! Terluna: the curved air's extra infrared escape\n"
     "      if (nlwcorr == 1) call lw_correction ! Terluna: the longwave by height\n"),
    ("      end subroutine twilight_lw\n",
     "      end subroutine twilight_lw\n"
     "\n"
     "      subroutine lw_correction\n"
     "      use radmod\n"
     "      implicit none\n"
     "!\n"
     "!     Terluna: the longwave by height (climate/gcm/plasim_longwave.py): each layer's and the ground's extra\n"
     "!     escape to space, fitted over land and over sea to the correlated-k model, scaled with its own\n"
     "!     temperature to the fourth power and carried up to the top as upward flux (negative in PlaSim)\n"
     "!\n"
     "      integer :: jhor, jlev\n"
     "      real :: zsum, zland\n"
     "      do jhor = 1, NHOR\n"
     "       zland = dls(jhor)\n"
     "       zsum = (zland*lwcl(NLEP) + (1.-zland)*lwcs(NLEP)) * (dt(jhor,NLEP) / lwct(NLEP))**4\n"
     "       dftu(jhor,NLEP) = dftu(jhor,NLEP) - zsum\n"
     "       dlwfl(jhor,NLEP) = dlwfl(jhor,NLEP) - zsum\n"
     "       do jlev = NLEV, 1, -1\n"
     "        zsum = zsum + (zland*lwcl(jlev) + (1.-zland)*lwcs(jlev)) * (dt(jhor,jlev) / lwct(jlev))**4\n"
     "        dftu(jhor,jlev) = dftu(jhor,jlev) - zsum\n"
     "        dlwfl(jhor,jlev) = dlwfl(jhor,jlev) - zsum\n"
     "       enddo\n"
     "      enddo\n"
     "      return\n"
     "      end subroutine lw_correction\n"),
]
# Ice optics for cold cloud (2026-10-17): PlaSim gives every cloud liquid water's optics, which for the Moon's thin,
# cold upper cloud reflect a seventh to a half of the sunlight ice crystals with the same water would, and at 63-81 km
# trap twice their infrared. NICEOPT 1 gives cloud colder than TICEOPT (-38 °C) CAM3's ice optics (plasim_ice_optics.py):
# Ebert and Curry's in sunlight through the atmosphere domain's PIFM two-stream, and their grey absorption in the
# infrared, at the effective radius Kristjánsson and Mitchell's table gives for the layer's temperature.
def _fortran_list(values, fmt, per_line=8):
    rows = [', '.join(fmt.format(v) for v in values[i:i + per_line]) for i in range(0, len(values), per_line)]
    return ' &\n     &   '.join(row + ',' if k < len(rows) - 1 else row for k, row in enumerate(rows))


ICE_ROUTINES = (
    "\n"
    "      subroutine ice_cloud_sw(jlev,zmu0,zcs,losun,zr1,zr1s,zr2,zt2,zr2s,zt2s)\n"
    "      use radmod\n"
    "      implicit none\n"
    "!\n"
    "!     Terluna: ice optics for cold cloud (climate/gcm/plasim_ice_optics.py). Where a sunlit layer's cloud is colder\n"
    "!     than ticeopt, its reflection and transmission come from Ebert and Curry's (1992) ice optics as CAM3 and RRTMG\n"
    "!     carry them, at the radius ice_radius gives, through the atmosphere domain's delta-scaled PIFM two-stream:\n"
    "!     band 1 from their visible band, band 2 as the sunlight-weighted sum of their four near-infrared bands. As in\n"
    "!     PlaSim, light under cloud higher up (1 - zcs) counts as diffuse.\n"
    "!\n"
    "      integer, intent(in) :: jlev\n"
    "      real, intent(in) :: zmu0(NHOR), zcs(NHOR)\n"
    "      logical, intent(in) :: losun(NHOR)\n"
    "      real, intent(inout) :: zr1(NHOR), zr1s(NHOR), zr2(NHOR), zt2(NHOR), zr2s(NHOR), zt2s(NHOR)\n"
    f"      real, parameter :: zea = {pio.EC_A:.3e}, zeb = {pio.EC_B:.3f}   ! extinction per g/m2: zea + zeb/r, r in microns\n"
    f"      real, parameter :: zec(5) = (/{_fortran_list(pio.EC_C, '{:.3e}')}/)   ! 1 - omega: zec + zed r\n"
    f"      real, parameter :: zed(5) = (/{_fortran_list(pio.EC_D, '{:.3e}')}/)\n"
    f"      real, parameter :: zee(5) = (/{_fortran_list(pio.EC_E, '{:.3e}')}/)   ! asymmetry: zee + zef r\n"
    f"      real, parameter :: zef(5) = (/{_fortran_list(pio.EC_F, '{:.3e}')}/)\n"
    f"      real, parameter :: zew(4) = (/{_fortran_list(pio.NIR_WEIGHTS, '{:.4f}')}/)   ! bands 2-5 in band 2\n"
    "      real, external :: ice_radius\n"
    "      integer :: jhor, jb\n"
    "      real :: zre, zwp, ztau, zmu, zrd, ztd, zrf, ztf, zr, zt, zrs, zts\n"
    "      do jhor = 1, NHOR\n"
    "       if (.not. losun(jhor)) cycle\n"
    "       if (dcc(jhor,jlev) <= 0. .or. dt(jhor,jlev) >= ticeopt) cycle\n"
    "       zre = ice_radius(dt(jhor,jlev))\n"
    "       zwp = min(1000.,1000.*dql(jhor,jlev)*dp(jhor)/ga*dsigma(jlev))\n"
    "       ztau = zwp*(zea+zeb/zre)\n"
    f"       zmu = max(zmu0(jhor),{pio.MU_FLOOR:.0e})\n"
    "       call ice_pifm(ztau,1.-zec(1)-zed(1)*zre,min(zee(1)+zef(1)*zre,0.9999),zmu,zrd,ztd,zrf,ztf)\n"
    "       zr1s(jhor) = zrf\n"
    "       zr1(jhor) = zcs(jhor)*zrd+(1.-zcs(jhor))*zrf\n"
    "       zr = 0.\n"
    "       zt = 0.\n"
    "       zrs = 0.\n"
    "       zts = 0.\n"
    "       do jb = 2, 5\n"
    "        call ice_pifm(ztau,1.-zec(jb)-zed(jb)*zre,min(zee(jb)+zef(jb)*zre,0.9999),zmu,zrd,ztd,zrf,ztf)\n"
    "        zr = zr+zew(jb-1)*zrd\n"
    "        zt = zt+zew(jb-1)*ztd\n"
    "        zrs = zrs+zew(jb-1)*zrf\n"
    "        zts = zts+zew(jb-1)*ztf\n"
    "       enddo\n"
    "       zr2s(jhor) = zrs\n"
    "       zt2s(jhor) = zts\n"
    "       zr2(jhor) = zcs(jhor)*zr+(1.-zcs(jhor))*zrs\n"
    "       zt2(jhor) = zcs(jhor)*zt+(1.-zcs(jhor))*zts\n"
    "      enddo\n"
    "      return\n"
    "      end subroutine ice_cloud_sw\n"
    "\n"
    "      subroutine ice_cloud_lw(jlev,ztaucc)\n"
    "      use radmod\n"
    "      implicit none\n"
    "!\n"
    "!     Terluna: ice optics for cold cloud in the infrared: Ebert and Curry's grey absorption, 0.005 + 1/r m2/g, in\n"
    "!     place of acllwr where the layer is colder than ticeopt, with PlaSim's diffusivity factor\n"
    "!\n"
    "      integer, intent(in) :: jlev\n"
    "      real, intent(inout) :: ztaucc(NHOR)\n"
    "      real, external :: ice_radius\n"
    "      integer :: jhor\n"
    "      real :: zwp\n"
    "      do jhor = 1, NHOR\n"
    "       if (dt(jhor,jlev) >= ticeopt) cycle\n"
    "       zwp = 1000.*dql(jhor,jlev)*dp(jhor)/ga*dsigma(jlev)\n"
    "       ztaucc(jhor) = 1.-dcc(jhor,jlev)*(1.-exp(-1.66*(0.005+1./ice_radius(dt(jhor,jlev)))*zwp))\n"
    "      enddo\n"
    "      return\n"
    "      end subroutine ice_cloud_lw\n"
    "\n"
    "      real function ice_radius(zt)\n"
    "      implicit none\n"
    "!\n"
    "!     Terluna: the ice effective radius (microns) at temperature zt (K): Kristjansson and Mitchell's table for\n"
    "!     hexagonal columns at 180-274 K as CAM3 and WRF's RRTMG carry it (reicalc), and below 180 K its own trend\n"
    "!\n"
    "      real, intent(in) :: zt\n"
    f"      real, parameter :: zretab(95) = (/{_fortran_list(pio.RETAB, '{:.6g}')}/)\n"
    "      integer :: ji\n"
    "      real :: zc\n"
    "      if (zt < 180.) then\n"
    "       ice_radius = zretab(1)*exp(log(zretab(2)/zretab(1))*(zt-180.))\n"
    "      else\n"
    "       ji = min(int(zt-179.),94)\n"
    "       zc = zt-real(179+ji)\n"
    "       ice_radius = zretab(ji)*(1.-zc)+zretab(ji+1)*zc\n"
    "      endif\n"
    "      return\n"
    "      end function ice_radius\n"
    "\n"
    "      subroutine ice_pifm(ptau,pom,pg,pmu,prdir,ptdir,prdif,ptdif)\n"
    "      implicit none\n"
    "!\n"
    "!     Terluna: a homogeneous layer's delta-scaled PIFM two-stream (Zdunkowski et al. 1980, Joseph et al. 1976) as\n"
    "!     in the atmosphere domain's shortwave.layer_properties: the direct beam's reflection and total transmission per\n"
    "!     unit flux on a horizontal surface, then the diffuse reflection and transmission\n"
    "!\n"
    "      real, intent(in) :: ptau, pom, pg, pmu\n"
    "      real, intent(out) :: prdir, ptdir, prdif, ptdif\n"
    "      real :: zf, zts, zoms, zgs, zg1, zg2, zg3, zg4, zlam, zmu, zgam, ze, zeb, zden, zdet, za, zb, zup, zdn\n"
    "      real :: za1, za2\n"
    "      zf = pg*pg\n"
    "      zts = (1.-pom*zf)*ptau\n"
    "      zoms = min((1.-zf)*pom/(1.-pom*zf),1.-1.e-9)\n"
    "      zgs = (pg-zf)/(1.-zf)\n"
    "      zg1 = (8.-zoms*(5.+3.*zgs))/4.\n"
    "      zg2 = 3.*zoms*(1.-zgs)/4.\n"
    "      zg3 = (2.-3.*zgs*pmu)/4.\n"
    "      zg4 = 1.-zg3\n"
    "      zlam = sqrt(max(zg1*zg1-zg2*zg2,1.e-30))\n"
    "      zmu = pmu\n"
    "      if (abs(1.-(zlam*zmu)**2) < 1.e-6) zmu = zmu*(1.+1.e-4)\n"
    "      zgam = zg2/(zg1+zlam)\n"
    "      ze = exp(-zlam*zts)\n"
    "      zeb = exp(-zts/zmu)\n"
    "      zden = 1.-zgam*zgam*ze*ze\n"
    "      prdif = zgam*(1.-ze*ze)/zden\n"
    "      ptdif = ze*(1.-zgam*zgam)/zden\n"
    "      zdet = 1./zmu**2-zlam*zlam\n"
    "      za = zg1+1./zmu\n"
    "      zb = zg1-1./zmu\n"
    "      zup = zoms*(-zb*zg3-zg2*zg4)/zdet\n"
    "      zdn = -zoms*(za*zg4+zg2*zg3)/zdet\n"
    "      za1 = (zgam*ze*zup*zeb-zdn)/zden\n"
    "      za2 = -zgam*ze*za1-zup*zeb\n"
    "      prdir = max((zgam*za1+za2*ze+zup)/zmu,0.)\n"
    "      ptdir = max((za1*ze+zgam*za2+zdn*zeb)/zmu,0.)+zeb\n"
    "      return\n"
    "      end subroutine ice_pifm\n")
ICE_OPTICS_EDITS = [
    ("      real :: lwct(NLEV+1) = 0.0 ! Terluna: the temperatures they are given at\n",
     "      real :: lwct(NLEV+1) = 0.0 ! Terluna: the temperatures they are given at\n"
     "      integer :: niceopt = 0     ! Terluna: ice optics for cold cloud, 1 for cloud colder than ticeopt\n"
     f"      real :: ticeopt = {pio.ICE_BELOW_K:.2f}   ! Terluna: the temperature below which cloud is ice\n"),
    ("     &               ,ntwilight,zmlunar,ntwilw,ngasbands,zgasq,zgasw,zmairlunar,nlwcorr\n",
     "     &               ,ntwilight,zmlunar,ntwilw,ngasbands,zgasq,zgasw,zmairlunar,nlwcorr,niceopt,ticeopt\n"),
    ("      call mpbci(nlwcorr)\n",
     "      call mpbci(nlwcorr)\n"
     "      call mpbci(niceopt)\n"
     "      call mpbcr(ticeopt)\n"),
    ("         zrcl2(:,jlev)=zcs(:)*zrcl2(:,jlev)+(1.-zcs(:))*zrcl2s(:,jlev)\n"
     "        endwhere\n",
     "         zrcl2(:,jlev)=zcs(:)*zrcl2(:,jlev)+(1.-zcs(:))*zrcl2s(:,jlev)\n"
     "        endwhere\n"
     "        if (niceopt == 1) call ice_cloud_sw(jlev,zmu0,zcs,losun,zrcl1(:,jlev),zrcl1s(:,jlev) & ! Terluna: ice\n"
     "     &                       ,zrcl2(:,jlev),ztcl2(:,jlev),zrcl2s(:,jlev),ztcl2s(:,jlev))\n"),
    ("        ztaucc0(:,jlev)=1.-dcc(:,jlev)*(1.-exp(zzf3*dql(:,jlev)*dp(:)))\n"
     "       endif\n",
     "        ztaucc0(:,jlev)=1.-dcc(:,jlev)*(1.-exp(zzf3*dql(:,jlev)*dp(:)))\n"
     "       endif\n"
     "       if (niceopt == 1 .and. clgray <= 0.) call ice_cloud_lw(jlev,ztaucc0(:,jlev)) ! Terluna: ice\n"),
    ("      end subroutine lw_correction\n",
     "      end subroutine lw_correction\n" + ICE_ROUTINES),
]
# The curved air's terms in the clear share (2026-10-19). The twilight patch puts its sunlight into each layer of every
# column and only the ground's into the column's clear share, and lets out its infrared in every column, cloudy or
# not, while the terms are clear-sky values (twilight.json, twilight_301k.json). NTWLCLEAR 1 adds each layer's sunlight
# in the share of the column clear of cloud above it, by PlaSim's own random overlap of cover (its zcs in swr), and lets
# each layer's and the ground's extra infrared through the clouds above it at PlaSim's own longwave transmissivity of
# each cloudy layer (ztaucc0 in lwr, ice optics included). 0, the default, leaves the model as it was.
TWILIGHT_CLEAR_EDITS = [
    ("      real :: twlt(NLEV+1) = 0.0 ! Terluna: the temperatures it is given at\n",
     "      real :: twlt(NLEV+1) = 0.0 ! Terluna: the temperatures it is given at\n"
     "      integer :: ntwlclear = 0   ! Terluna: 1 puts the twilight's terms in the clear share above each layer\n"
     "      real :: twltau(NHOR,NLEV) = 1.0 ! Terluna: each layer's longwave cloud transmissivity, from lwr\n"),
    ("     &               ,ntwilight,zmlunar,ntwilw,ngasbands,zgasq,zgasw,zmairlunar,nlwcorr,niceopt,ticeopt\n",
     "     &               ,ntwilight,zmlunar,ntwilw,ngasbands,zgasq,zgasw,zmairlunar,nlwcorr,niceopt,ticeopt  &\n"
     "     &               ,ntwlclear\n"),
    ("      call mpbcr(ticeopt)\n", "      call mpbcr(ticeopt)\n      call mpbci(ntwlclear)\n"),
    ("      real :: zel, zu, zsum, zclear, zsun\n", "      real :: zel, zu, zsum, zclear, zsun, zabove(NLEV)\n"),
    ("       zclear = 1.0\n"
     "       do jlev = 1, NLEV\n"
     "        zclear = zclear * (1.0 - dcc(jhor,jlev) * nclouds)\n"
     "       enddo\n",
     "       zclear = 1.0\n"
     "       do jlev = 1, NLEV\n"
     "        zabove(jlev) = zclear      ! Terluna: the clear share above each layer\n"
     "        zclear = zclear * (1.0 - dcc(jhor,jlev) * nclouds)\n"
     "       enddo\n"
     "       if (ntwlclear /= 1) zabove(:) = 1.0\n"),
    ("        zsum = zsum + zsun * ((1.0 - zu) * twlq(jlev,jlo) + zu * twlq(jlev,jhi))\n",
     "        zsum = zsum + zsun * ((1.0 - zu) * twlq(jlev,jlo) + zu * twlq(jlev,jhi)) * zabove(jlev)\n"),
    ("      integer :: jhor, jlev\n"
     "      real :: zsum\n"
     "      do jhor = 1, NHOR\n"
     "       zsum = twlw(NLEP) * (dt(jhor,NLEP) / twlt(NLEP))**4\n",
     "      integer :: jhor, jlev\n"
     "      real :: zsum, zclear, zabove(NLEV)\n"
     "      do jhor = 1, NHOR\n"
     "       zclear = 1.0\n"
     "       do jlev = 1, NLEV\n"
     "        zabove(jlev) = zclear      ! Terluna: the longwave's clear share above each layer\n"
     "        if (ntwlclear == 1) zclear = zclear * twltau(jhor,jlev)\n"
     "       enddo\n"
     "       zsum = twlw(NLEP) * (dt(jhor,NLEP) / twlt(NLEP))**4 * zclear\n"),
    ("        zsum = zsum + twlw(jlev) * (dt(jhor,jlev) / twlt(jlev))**4\n",
     "        zsum = zsum + twlw(jlev) * (dt(jhor,jlev) / twlt(jlev))**4 * zabove(jlev)\n"),
    ("       if (niceopt == 1 .and. clgray <= 0.) call ice_cloud_lw(jlev,ztaucc0(:,jlev)) ! Terluna: ice\n",
     "       if (niceopt == 1 .and. clgray <= 0.) call ice_cloud_lw(jlev,ztaucc0(:,jlev)) ! Terluna: ice\n"
     "       if (ntwlclear == 1) twltau(:,jlev) = ztaucc0(:,jlev) ! Terluna: for the twilight's clear share\n"),
]
# The ice optics' weights from the namelist (2026-10-19): world-crafting's patch 11 (its lab branch at 79ea224), which
# turns the four near-infrared weights of the ice patch into ZICEW with the ice patch's constants as its default, so a
# run that writes nothing is the model as it was. Runs with ice optics write the shares of the design shield's
# sunlight (plasim_ice_optics.SHIELD_NIR_WEIGHTS) unless --ice-weights stock keeps the unfiltered Sun's.
ICE_WEIGHTS_EDITS = [
    (f"      real :: ticeopt = {pio.ICE_BELOW_K:.2f}   ! Terluna: the temperature below which cloud is ice\n",
     f"      real :: ticeopt = {pio.ICE_BELOW_K:.2f}   ! Terluna: the temperature below which cloud is ice\n"
     f"      real :: zicew(4) = (/{_fortran_list(pio.NIR_WEIGHTS, '{:.4f}')}/) "
     "! world-crafting: ice bands 2-5 in band 2, patch 10's\n"),
    ("     &               ,ntwlclear\n", "     &               ,ntwlclear,zicew\n"),
    ("      call mpbci(ntwlclear)\n", "      call mpbci(ntwlclear)\n      call mpbcrn(zicew,4)\n"),
    (f"      real, parameter :: zew(4) = (/{_fortran_list(pio.NIR_WEIGHTS, '{:.4f}')}/)   ! bands 2-5 in band 2\n", ""),
    ("        zr = zr+zew(jb-1)*zrd\n", "        zr = zr+zicew(jb-1)*zrd\n"),
    ("        zt = zt+zew(jb-1)*ztd\n", "        zt = zt+zicew(jb-1)*ztd\n"),
    ("        zrs = zrs+zew(jb-1)*zrf\n", "        zrs = zrs+zicew(jb-1)*zrf\n"),
    ("        zts = zts+zew(jb-1)*ztf\n", "        zts = zts+zicew(jb-1)*ztf\n"),
]
# The cold cloud's lunar factors (2026-10-20): ZICEWF and ZICERF scale, level by level, the water path and the crystal
# radius the ice optics give cloud colder than TICEOPT, in sunlight and in the infrared alike (cold_cloud.py: the
# factors that carry PlaSim's Earth water and Kristjansson and Mitchell's Earth crystals to lunar gravity). Their
# defaults of 1 leave the model as it was; the factors touch only how the radiation sees the cloud.
ICE_FACTOR_EDITS = [
    ("      real :: zicew(4) = (/" + _fortran_list(pio.NIR_WEIGHTS, '{:.4f}') + "/) "
     "! world-crafting: ice bands 2-5 in band 2, patch 10's\n",
     "      real :: zicew(4) = (/" + _fortran_list(pio.NIR_WEIGHTS, '{:.4f}') + "/) "
     "! world-crafting: ice bands 2-5 in band 2, patch 10's\n"
     "      real :: zicewf(NLEV) = 1.0 ! Terluna: the cold cloud's water factor by level\n"
     "      real :: zicerf(NLEV) = 1.0 ! Terluna: the cold cloud's crystal radius factor by level\n"),
    ("     &               ,ntwlclear,zicew\n", "     &               ,ntwlclear,zicew,zicewf,zicerf\n"),
    ("      call mpbcrn(zicew,4)\n", "      call mpbcrn(zicew,4)\n      call mpbcrn(zicewf,NLEV)\n      call mpbcrn(zicerf,NLEV)\n"),
    ("       zre = ice_radius(dt(jhor,jlev))\n", "       zre = ice_radius(dt(jhor,jlev))*zicerf(jlev)\n"),
    ("       zwp = min(1000.,1000.*dql(jhor,jlev)*dp(jhor)/ga*dsigma(jlev))\n",
     "       zwp = min(1000.,1000.*dql(jhor,jlev)*dp(jhor)/ga*dsigma(jlev)*zicewf(jlev))\n"),
    ("       zwp = 1000.*dql(jhor,jlev)*dp(jhor)/ga*dsigma(jlev)\n"
     "       ztaucc(jhor) = 1.-dcc(jhor,jlev)*(1.-exp(-1.66*(0.005+1./ice_radius(dt(jhor,jlev)))*zwp))\n",
     "       zwp = 1000.*dql(jhor,jlev)*dp(jhor)/ga*dsigma(jlev)*zicewf(jlev)\n"
     "       ztaucc(jhor) = 1.-dcc(jhor,jlev)*(1.-exp(-1.66*(0.005+1./(ice_radius(dt(jhor,jlev))*zicerf(jlev)))*zwp))\n"),
]
# The 10 m wind (2026-10-19). PlaSim 3.4.2 writes no wind below its lowest level, about 0.9 km up in the lunar air, and
# the lake model brought that wind down on a neutral log profile. NWIND10 1 has the surface layer (surflx) take the
# 10 m wind on Geleyn's (1988) profile through the lowest layer, u(z) = u_L (ln(1 + z/z0) - (z/z_L)(B_N - B_D)) / B_D,
# from the drag the model itself uses: B_N = ln(1 + z_L/z0) = k/sqrt(C_DN) and B_D = k/sqrt(C_DN f_m), f_m the
# stability factor of PlaSim's momentum transfer (zrifm). The neutral limit is the log profile; stable air brings less
# of the wind down, unstable air more. The output writes its mean over the steps since the last output, codes 165 and
# 166 (uas, vas), counting its own steps, since PlaSim reads its restart before this namelist. 0, the default, leaves
# the model as it was.
WIND10_MOD_EDITS = [
    ("      real :: atsa(NHOR)  = 0. ! accumulated surface air temperature\n",
     "      real :: atsa(NHOR)  = 0. ! accumulated surface air temperature\n"
     "      integer :: nwind10 = 0      ! Terluna: 1 writes the 10 m wind (codes 165 and 166)\n"
     "      integer :: nu10acc = 0      ! Terluna: the 10 m wind's steps since its last output\n"
     "      real :: du10(NHOR)  = 0.    ! Terluna: the 10 m wind (m/s)\n"
     "      real :: dv10(NHOR)  = 0.\n"
     "      real :: au10(NHOR)  = 0.    ! Terluna: its sum since its last output\n"
     "      real :: av10(NHOR)  = 0.\n"),
]
WIND10_FLUX_EDITS = [
    ("      namelist/fluxmod_nl/nvdiff,nshfl,nevap,nstress,ntsa",
     "      namelist/fluxmod_nl/nvdiff,nshfl,nevap,nstress,ntsa,nwind10"),
    ("      call mpbci(ntsa)\n", "      call mpbci(ntsa)\n      call mpbci(nwind10)\n"),
    ("      real znl(NHOR)       ! z of lowermost level (m)\n",
     "      real znl(NHOR)       ! z of lowermost level (m)\n"
     "      real zbn,zbd,zf10    ! Terluna: the 10 m wind, Geleyn's profile\n"),
    ("       dtransm(jhor)=ztrans*zrifm(jhor)\n",
     "       dtransm(jhor)=ztrans*zrifm(jhor)\n"
     "       if (nwind10 == 1) then   ! Terluna: the 10 m wind from the lowest layer's drag\n"
     "        zbn=LOG(zbz0(jhor)+1.)\n"
     "        zbd=zbn/SQRT(zrifm(jhor))\n"
     "        zf10=(LOG(1.+10./dz0(jhor))-10./znl(jhor)*(zbn-zbd))/zbd\n"
     "        zf10=MAX(0.,MIN(1.,zf10))\n"
     "        du10(jhor)=zf10*du(jhor,NLEV)\n"
     "        dv10(jhor)=zf10*dv(jhor,NLEV)\n"
     "       endif\n"),
]
WIND10_OUT_EDITS = [
    ("      call writegp(40,atsa,167,0)\n",
     "      call writegp(40,atsa,167,0)\n"
     "      if (nwind10 == 1) then   ! Terluna: the 10 m wind's mean since its last output\n"
     "       au10(:)=au10(:)/real(max(nu10acc,1))\n"
     "       call writegp(40,au10,165,0)\n"
     "       av10(:)=av10(:)/real(max(nu10acc,1))\n"
     "       call writegp(40,av10,166,0)\n"
     "      endif\n"),
    ("      atsa(:)=atsa(:)+dtsa(:)\n",
     "      atsa(:)=atsa(:)+dtsa(:)\n"
     "      if (nwind10 == 1) then\n"
     "       au10(:)=au10(:)+du10(:)\n"
     "       av10(:)=av10(:)+dv10(:)\n"
     "       nu10acc=nu10acc+1\n"
     "      endif\n"),
    ("      atsa(:)=0.\n", "      atsa(:)=0.\n      au10(:)=0.\n      av10(:)=0.\n      nu10acc=0\n"),
]
WIND10_PATCHES = [('plasimmod.f90', 'Terluna: 1 writes the 10 m wind', WIND10_MOD_EDITS),
                  ('fluxmod.f90', "Terluna: the 10 m wind, Geleyn's profile", WIND10_FLUX_EDITS),
                  ('outmod.f90', "Terluna: the 10 m wind's mean since its last output", WIND10_OUT_EDITS)]
# The free-convection constant over water (2026-10-26). In unstable air over water surflx raises the transfer of heat and
# vapour by the interpolation of Miller, Beljaars and Palmer (1992, J. Climate 5, 418-434, their eq. 9),
# C_H = C_HN (1 + C_R^1.25)^0.8 with C_R = 0.0016 (theta_vs - theta_v1)^(1/3) / (C_HN |U|), whose still-air limit,
# C_H |U| = 0.0016 (theta_vs - theta_v1)^(1/3), stands for the free convection of tank experiments (their eq. 8):
# E / rho = 0.17 (g kappa^2 / (theta_v nu))^(1/3) (theta_vs - theta_v1)^(1/3) (q_s - q_1), with kappa the air's thermal
# diffusivity and nu its kinematic viscosity. At 300 K and 1013 hPa that law's constant is 0.0017; they took 0.0016, to
# sit near Monin-Obukhov similarity at low wind speeds. With kappa = k / (rho c_p), nu = mu / rho and rho theta_v = p / R
# at the ground, the constant goes as (g / p_s)^(1/3) at a given temperature, the conductivity and viscosity adding
# about T^0.3 (free_convection_setting). FREE_CONVECTION_EDITS make it fluxmod's namelist variable CFREE, PlaSim's 0.0016
# unless a run sets it, so a run without the setting computes as before.
FREE_CONVECTION_MARKER = 'Terluna: the free-convection constant over water'
FREE_CONVECTION_EDITS = [
    ('      real :: vdiff_d    = 5.    !        "\n',
     '      real :: vdiff_d    = 5.    !        "\n'
     f'      real :: cfree      = 0.0016 ! {FREE_CONVECTION_MARKER} (m s-1 K-1/3; Miller et al. 1992)\n'),
    ("     &                ,zumin,vdiff_lamm,vdiff_b,vdiff_c,vdiff_d\n",
     "     &                ,zumin,vdiff_lamm,vdiff_b,vdiff_c,vdiff_d,cfree\n"),
    ("      call mpbcr(vdiff_d)\n", "      call mpbcr(vdiff_d)\n      call mpbcr(cfree)\n"),
    ("         zrifh(jhor)=(1.+(0.0016*zdth**(1./3.)/SQRT(zabsu2(jhor))       &\n",
     "         zrifh(jhor)=(1.+(cfree*zdth**(1./3.)/SQRT(zabsu2(jhor))        &\n"),
]
FREE_CONVECTION_PATCHES = [('fluxmod.f90', FREE_CONVECTION_MARKER, FREE_CONVECTION_EDITS)]
# The twilight, gas, longwave and ice patches go only into a copy of ExoPlaSim of its own, which every run imports
# first since the corrected radiation became the standard (2026-10-14); the shared installation stays as the earlier
# runs had it. The twilight's clear share follows them, from 2026-10-19, and the ice weights and the 10 m wind after it,
# and the cold cloud's lunar factors from 2026-10-20.
TWILIGHT_PATCHES = [('radmod.f90', "Terluna: the curved air's twilight table", TWILIGHT_EDITS),
                    ('radmod.f90', 'Terluna: CO2 and O2 in the sunlight', GAS_EDITS),
                    ('radmod.f90', 'Terluna: the longwave by height', LW_CORRECTION_EDITS),
                    ('radmod.f90', 'Terluna: ice optics for cold cloud', ICE_OPTICS_EDITS),
                    ('radmod.f90', "Terluna: 1 puts the twilight's terms in the clear share", TWILIGHT_CLEAR_EDITS),
                    ('radmod.f90', 'world-crafting: ice bands 2-5 in band 2', ICE_WEIGHTS_EDITS)] + WIND10_PATCHES
TWILIGHT_PATCHES += [('radmod.f90', "Terluna: the cold cloud's water factor by level", ICE_FACTOR_EDITS)]
TWILIGHT_PACKAGE = Path('/media/projectspace/terluna-research/exoplasim-twilight')
# The curved air's ledgers (twilight.py): on the design climate's 301 K column, the standard since 2026-10-19, and on
# the 294.9 K column every run from 14 to 19 October took. The Moon's slant factors (ZMLUNAR, ZMAIRLUNAR) stay the 294.9 K
# ledger's, with which the shortwave calibration was refitted (shortwave_fit_A28_dim5_lunar.json); the 301 K ledger
# records its own column's beside them.
TWILIGHT_LEDGERS = {'301k': ROOT / 'climate' / 'results' / 'gcm' / 'twilight_301k.json',
                    '294.9k': ROOT / 'climate' / 'results' / 'gcm' / 'twilight.json'}
TWILIGHT_PRODUCT = TWILIGHT_LEDGERS['301k']
SLANT_PRODUCT = TWILIGHT_LEDGERS['294.9k']
# The longwave by height (plasim_longwave.py): fitted on the design case's columns near 295 K, and refitted on 19 October
# on the columns of the 20-year test with the correction as it passed 301.2 K (its year 2), since the first does not hold
# there (longwave_check_301k.json); new runs with --longwave-correction take the refit. 'l15_provisional' (2026-10-26) is
# the correction refitted on the 15-level month trial's own cloud-free columns (plasim_longwave.py refit), for tests on
# those levels until the configuration chosen for a long run is refitted; a run on another set refuses it.
LONGWAVE_FITS = {'301k': ROOT / 'climate' / 'results' / 'gcm' / 'longwave_A28_corrected_970_lw_ice.json',
                 '295k': ROOT / 'climate' / 'results' / 'gcm' / 'longwave_A28_corrected_970_heating.json',
                 'l15_provisional': ROOT / 'climate' / 'results' / 'gcm'
                 / 'longwave_A28R_301k_cold_sbm_rh07_deep200_l15cap_month_refit.json'}
LONGWAVE_PRODUCT = LONGWAVE_FITS['301k']
GASES_PRODUCT = ROOT / 'climate' / 'results' / 'gcm' / 'shortwave_gases.json'
GAS_ORDER = ('co2_nir', 'o2_nir', 'o2_vis')                     # PlaSim's ZGASQ columns
MOLAR = dict(pN2=0.0280134, pO2=0.0319988, pAr=0.039948, pCO2=0.0440095)   # kg/mol
# 'curved' adds the curved light and the curved escape of infrared, 'full' the shells' widening on both sides as well
# (twilight.json's cone_w_m2 and longwave.cone), and both put in the Moon's slant path; 'path' puts in that alone.
TWILIGHT_MODES = ('curved', 'full', 'path')


def use_twilight_package():
    """Put the twilight copy of ExoPlaSim first on the import path, copying the installed package there first if
    it is missing; before anything imports exoplasim."""
    import importlib.util
    if 'exoplasim' in sys.modules:
        raise RuntimeError('exoplasim is already imported; the twilight copy must come first')
    target = TWILIGHT_PACKAGE / 'exoplasim'
    if not target.exists():
        source = Path(importlib.util.find_spec('exoplasim').origin).parent
        tmp = TWILIGHT_PACKAGE / 'exoplasim.part'
        shutil.copytree(source, tmp, ignore=shutil.ignore_patterns('__pycache__'))
        tmp.rename(target)
    sys.path.insert(0, str(TWILIGHT_PACKAGE))


def ensure_patched(twilight: bool = False, trigger: bool = False, depth: bool = False, cap: bool = False,
                   massflux: bool = False, free_convection: bool = False, massflux_start: bool = False,
                   boundary_layer: bool = False, night: bool = False, humidity_mean: bool = False,
                   cloud_water: bool = False, gust: bool = False) -> dict:
    """Apply the Terluna patches to ExoPlaSim's source where missing (removing the compiled executables
    so ExoPlaSim rebuilds them), with the twilight patch in the twilight copy, and the scheme's CAPE trigger
    (SBM_TRIGGER_PATCHES), its cloud-depth test and mixing within layers (SBM_DEPTH_PATCHES) and its shallow cap
    (SBM_CAP_PATCHES, after the depth patch), the mass-flux scheme (massflux_fortran.MASSFLUX_PATCHES, over the
    Betts-Miller patch) with its surface parcel's starting speed (massflux_fortran.START_PATCHES, over it), and the
    free-convection constant over water (FREE_CONVECTION_PATCHES, after the 10 m wind's), and the nonlocal boundary
    layer with the convection's vapour diagnostic (pbl_fortran.PBL_PATCHES over the free-convection patch and
    DIAGNOSTIC_PATCHES over the Betts-Miller and mass-flux patches), and the night switch, a diagnostic
    (SBM_NIGHT_PATCHES, over the CAPE trigger's, the depth test's and the cap's patches, with the convection's vapour
    diagnostic), and the switches in PlaSim's namelist (PLASIM_SWITCH_PATCHES) with the humidity's true mean
    (HUS_MEAN_PATCHES) and the cloud water (condensate.CLOUD_WATER_PATCHES, with the convection's vapour diagnostic),
    and Beljaars's gust in the surface transfer (surface_gust.GUST_PATCHES, over the nonlocal boundary layer's patch,
    with the convection's vapour diagnostic), for a run that sets them, and return each patched file's hash for the
    run key."""
    import exoplasim
    src = Path(exoplasim.__file__).parent / 'plasim' / 'src'
    if twilight != str(src).startswith(str(TWILIGHT_PACKAGE)):
        raise RuntimeError(f'the twilight patch belongs only in {TWILIGHT_PACKAGE}, not in {src}')
    if free_convection and not twilight:
        raise RuntimeError('the free-convection constant goes into the twilight copy, after the 10 m wind patch')
    if massflux_start and not massflux:
        raise RuntimeError("the starting speed's patch goes over the mass-flux scheme's")
    if boundary_layer and not twilight:
        raise RuntimeError('the boundary-layer patch goes into the twilight copy, over the free-convection patch')
    if gust and not twilight:
        raise RuntimeError("the gust's patch goes into the twilight copy, over the boundary-layer patch")
    if night and not twilight:
        raise RuntimeError("the night switch reads the twilight patch's cosine below the horizon; it goes into the "
                           'twilight copy')
    hashes = {}
    for name, marker, edits in (PATCHES + (TWILIGHT_PATCHES if twilight else [])
                                + (SBM_TRIGGER_PATCHES if trigger or night else [])
                                + (SBM_DEPTH_PATCHES if depth or cap or night else [])
                                + (SBM_CAP_PATCHES if cap or night else [])
                                + (mff.MASSFLUX_PATCHES if massflux or boundary_layer or night or cloud_water
                                   else [])
                                + (mff.START_PATCHES if massflux_start else [])
                                + (FREE_CONVECTION_PATCHES if free_convection or boundary_layer or gust else [])
                                + (pbf.PBL_PATCHES if boundary_layer or gust else [])
                                + (sgf.GUST_PATCHES if gust else [])
                                + (pbf.DIAGNOSTIC_PATCHES if boundary_layer or night or cloud_water or gust else [])
                                + (SBM_NIGHT_PATCHES if night else [])
                                + (PLASIM_SWITCH_PATCHES if humidity_mean or cloud_water else [])
                                + (HUS_MEAN_PATCHES if humidity_mean else [])
                                + (cnd.CLOUD_WATER_PATCHES if cloud_water else [])):
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


PLASIM_RAIN_GAMMA = 0.01   # rainmod's gamma at T21: the share of a layer's saturation deficit its rain evaporates a step


def rain_evaporation_setting(gamma) -> dict:
    """PlaSim's evaporation of falling rain and snow as the configuration records it: rainmod's gamma, the share of
    each unsaturated layer's saturation deficit the precipitation passing through it evaporates in a time step (with
    the latent cooling's correction, at most the precipitation there), in place of PlaSim's own 0.01."""
    gamma = float(gamma)
    if not 0.0 < gamma <= 1.0:
        raise ValueError('the rain evaporation constant must lie above 0 and at most 1')
    return dict(gamma=gamma, plasim_default=PLASIM_RAIN_GAMMA)


FREE_CONVECTION_CHOICES = ('earth', 'lunar')


def free_convection_setting(choice: str, gravity: float, pressure_pa: float) -> dict:
    """The free-convection constant over water (FREE_CONVECTION_EDITS) as the configuration records it. 'lunar' takes
    PlaSim's 0.0016 m s-1 K-1/3, set at Earth's gravity and sea-level pressure, to the planet's gravity and the
    experiment's surface pressure: the tank experiments' law behind it, 0.17 (g kappa^2 / (theta_v nu))^(1/3), goes as
    (g / p_s)^(1/3) at a given temperature, since kappa^2 / nu = k^2 / (rho c_p^2 mu) and rho theta_v = p_s / R at the
    ground, with the conductivity k and the viscosity mu set by the temperature alone. The sea's temperature, near 300 K on
    both worlds, and the air's composition are taken as Earth's; their factor, about (T / T_E)^0.3, is left out.
    'earth' is PlaSim's own constant, which the configuration leaves unrecorded."""
    if choice not in FREE_CONVECTION_CHOICES[1:]:
        raise ValueError(f"the free-convection constant is set only for 'lunar' (PlaSim's Earth value is the default), "
                         f"not {choice!r}")
    if not (gravity > 0.0 and pressure_pa > 0.0):
        raise ValueError('the free-convection constant needs a positive gravity and surface pressure')
    factor = (gravity / STANDARD_GRAVITY * STANDARD_ATMOSPHERE / pressure_pa) ** (1.0 / 3.0)
    return dict(choice=choice, constant=EXOPLASIM_FREE_CONVECTION * factor, earth_constant=EXOPLASIM_FREE_CONVECTION,
                factor=factor, gravity_m_s2=float(gravity), surface_pressure_pa=float(pressure_pa),
                earth_reference=dict(gravity_m_s2=STANDARD_GRAVITY, surface_pressure_pa=STANDARD_ATMOSPHERE),
                scaling='(g / p_s)^(1/3), from 0.17 (g kappa^2 / (theta_v nu))^(1/3) with rho theta_v = p_s / R',
                source='Miller, Beljaars and Palmer (1992), J. Climate 5, 418-434, eqs. (8) and (9)')


MIXING_LENGTH_CHOICES = ('earth', 'lunar')
BOUNDARY_LAYER_CHOICES = ('local', 'nonlocal')


def mixing_length_setting(choice: str, gravity: float) -> dict:
    """PlaSim's Blackadar length in vdiff as the configuration records it. 'lunar' stretches vdiff_lamm, 160 m for
    momentum (sqrt(7.5) times that for heat and vapour), by the ratio of the gravities, the rule climate/crm/cm1_run.py
    applies to CM1's asymptotic length l_inf: a column stretched by g_E / g keeps its velocities and pressures and
    lengthens its lengths by that ratio. It is a scaling assumption the CM1 boxes share, with no independent support
    yet, so runs take it as a sensitivity. 'earth' is PlaSim's own length, which the configuration leaves unrecorded."""
    if choice not in MIXING_LENGTH_CHOICES[1:]:
        raise ValueError(f"the mixing length is set only for 'lunar' (PlaSim's Earth value is the default), not {choice!r}")
    if not gravity > 0.0:
        raise ValueError('the lunar mixing length needs a positive gravity')
    momentum = pbs.lunar_length(gravity)
    return dict(choice=choice, momentum_m=momentum, heat_m=momentum * (3.0 * pbs.STABILITY * 0.5) ** 0.5,
                earth_momentum_m=pbs.LAMBDA_EARTH_M, stretch=STANDARD_GRAVITY / gravity, gravity_m_s2=float(gravity),
                rule="vdiff_lamm times g_E / g, as climate/crm/cm1_run.py stretches CM1's l_inf",
                status='a scaling assumption the CM1 boxes share, taken as a sensitivity')


def boundary_layer_setting(choice: str) -> dict:
    """The nonlocal boundary layer in vdiff (pbl_scheme.py, pbl_fortran.PBL_PATCHES) as the configuration records it,
    with every constant fluxmod's namelist receives; 'local' is PlaSim's own diffusion, left unrecorded."""
    if choice not in BOUNDARY_LAYER_CHOICES[1:]:
        raise ValueError(f"the boundary layer is set only for 'nonlocal' (PlaSim's local diffusion is the default), "
                         f'not {choice!r}')
    return dict(scheme=choice, parameters=pbf.fortran_values(),
                sources=('Holtslag and Boville (1993), J. Climate 6, 1825-1842, as CAM3 implements it (Collins and '
                         'others 2004, NCAR/TN-464+STR, section 4.11.2; hb_diff.F90); the entrainment flux of Hong, Noh '
                         'and Dudhia (2006), Mon. Wea. Rev. 134, 2318-2341'),
                code='climate/gcm/pbl_scheme.py (nonlocal_column), as Fortran in climate/gcm/pbl_fortran.py (pblcol)')


SURFACE_GUST_CHOICES = ('none', 'beljaars')


def surface_gust_setting(choice: str) -> dict:
    """Beljaars's gust in PlaSim's surface transfer (surface_gust.py, GUST_PATCHES) as the configuration records it:
    beta w* joins the lowest level's wind in the bulk formulas over every surface, the wind there at least 0.1 m/s, and
    over water the free-convection term goes; w* comes from the surface fluxes of the step before and the boundary
    layer's depth by the nonlocal scheme's rule and constants. 'none' is PlaSim's own transfer, left unrecorded."""
    if choice not in SURFACE_GUST_CHOICES[1:]:
        raise ValueError(f"the surface gust is set only for 'beljaars' (PlaSim's own transfer is the default), "
                         f'not {choice!r}')
    return dict(form=choice, beta=sgf.GUST_BETA, least_wind_m_s=sgf.GUST_UMIN,
                replaces=dict(free_convection_over_water=sgf.EXOPLASIM_CFREE, least_wind_m_s=sgf.EXOPLASIM_ZUMIN),
                depth=dict(rule="the nonlocal scheme's bulk Richardson number with the thermals' excess (pblcol)",
                           parameters=pbf.fortran_values()),
                fluxes='the surface fluxes of the step before (dshfl, devap, dust3)',
                sources=('Beljaars (1995), Q. J. R. Meteorol. Soc. 121, 255-270; CM1 sfclay (vconvc 1, least wind '
                         '0.1 m/s); Fairall and others (2003), J. Climate 16, 571-591'),
                code='climate/gcm/surface_gust.py (surface_transfer, gust_column; gustnl and gustcol in Fortran)')


SBM_DEFAULTS = dict(tau_s=7200.0, rh=0.7, eps=0.0, cin_j_kg=0.0)   # Frierson's 2 h and 70%; the lunar options off
SBM_TRIGGER_SETTINGS = ('cape_j_kg', 'trigger_eps')                   # the CAPE trigger, recorded only when set
SBM_DEPTH_SETTINGS = ('deep_hpa', 'mixing_steps')                     # SBM_DEPTH_PATCHES' options, likewise
SBM_CAP_SETTINGS = ('shallow_cap',)                                   # SBM_CAP_PATCHES' option, likewise
SBM_NIGHT_SETTINGS = ('diagnostic_night',)                            # SBM_NIGHT_PATCHES' diagnostic, likewise
SBM_NIGHT_RECORD = dict(rule='no deep convection where the geometric cosine of the Sun\'s zenith at the grid point is '
                             'at most 0; such a column takes PlaSim\'s shallow scheme, its top held by the shallow cap '
                             'where that is set',
                        status='a diagnostic, never a configuration (27 October 2026)',
                        patch='climate/gcm/exoplasim_run.py SBM_NIGHT_PATCHES (rainmod ndiagnight_bm)')


def convection_setting(scheme='sbm', **values) -> dict:
    """PlaSim's deep convection as the configuration records it: 'massflux' for the mass-flux scheme (massflux.py, as
    massflux_fortran.py patches it in, rainmod's nprc 3) with every one of its parameters as the routine's defaults give
    them, which takes no setting of the Betts-Miller scheme's; or 'sbm' for the simplified Betts-Miller scheme of the
    Terluna patch (SBM_EDITS, rainmod's nprc 2) with its relaxation time tau_s, reference relative humidity rh, the
    parcel's mixing per e-folding of pressure eps and its inhibition limit cin_j_kg (0: none), each defaulting to
    SBM_DEFAULTS; left out, the configuration keeps Kuo's scheme (nprc 1). The CAPE trigger (SBM_TRIGGER_PATCHES)
    joins the record only when cape_j_kg is set above 0: the CAPE a parcel mixing trigger_eps per e-folding of
    pressure (default 0) must have for the column to convect deep. The options of SBM_DEPTH_PATCHES join it only
    when set: deep_hpa above 0, the depth (hPa) the parcel's cloud must exceed for deep convection, and mixing_steps
    above 1, the steps of the parcel's mixing within each layer, which need eps above 0. The shallow cap of
    SBM_CAP_PATCHES joins it only when shallow_cap is true, and it needs the depth test. The night switch of
    SBM_NIGHT_PATCHES, a diagnostic, joins it as SBM_NIGHT_RECORD only when diagnostic_night is true."""
    if scheme == 'massflux':
        values = dict(values)
        start = values.pop('start', None)
        given = sorted(k for k, v in values.items() if v not in (None, False))
        if given:
            raise ValueError(f'the mass-flux scheme takes none of the Betts-Miller settings: {given}')
        out = dict(scheme='massflux', nprc=3, parameters=mff.fortran_values())
        if start not in (None, 'ifs'):                    # the surface parcel's starting speed (massflux_fortran)
            out['start'] = dict(rule=start, parameters=mff.start_values(dict(mff.mf.PARAMETERS, start=start)))
        return out
    if scheme != 'sbm':
        raise ValueError(f"the convection scheme is 'sbm' or 'massflux' (Kuo's is the default, left out), "
                         f"not {scheme!r}")
    unknown = (set(values) - set(SBM_DEFAULTS) - set(SBM_TRIGGER_SETTINGS) - set(SBM_DEPTH_SETTINGS)
               - set(SBM_CAP_SETTINGS) - set(SBM_NIGHT_SETTINGS))
    if unknown:
        raise ValueError(f'unknown settings of the simplified Betts-Miller scheme: {sorted(unknown)}')
    s = {k: float(values.get(k, v) if values.get(k) is not None else v) for k, v in SBM_DEFAULTS.items()}
    if s['tau_s'] <= 0.0 or not 0.0 < s['rh'] <= 1.0 or s['eps'] < 0.0 or s['cin_j_kg'] < 0.0:
        raise ValueError(f'the simplified Betts-Miller settings are out of range: {s}')
    cape, mix = values.get('cape_j_kg'), values.get('trigger_eps')
    if cape is not None and float(cape) > 0.0:
        s.update(cape_j_kg=float(cape), trigger_eps=float(mix or 0.0))
        if s['trigger_eps'] < 0.0:
            raise ValueError(f"the CAPE trigger's mixing must not be negative: {s['trigger_eps']}")
    elif (cape is not None and float(cape) < 0.0) or (mix is not None and float(mix) != 0.0):
        raise ValueError('the CAPE trigger takes a positive cape_j_kg, and its mixing only with it')
    deep, steps = values.get('deep_hpa'), values.get('mixing_steps')
    if deep is not None:
        if not 0.0 <= float(deep) < 1000.0:
            raise ValueError(f'the cloud-depth test takes a depth from 0 to 1000 hPa, not {deep}')
        if float(deep) > 0.0:
            s['deep_hpa'] = float(deep)
    if steps is not None:
        if int(steps) != steps or int(steps) < 1:
            raise ValueError(f'the mixing within layers takes a whole number of steps from 1, not {steps}')
        if int(steps) > 1:
            if s['eps'] <= 0.0:
                raise ValueError('the mixing within layers needs the parcel to mix (eps above 0)')
            s['mixing_steps'] = int(steps)
    if values.get('shallow_cap'):
        if 'deep_hpa' not in s:
            raise ValueError("the shallow cap holds a top within the depth test's threshold, which needs the depth test")
        s['shallow_cap'] = True
    if values.get('diagnostic_night'):
        s['diagnostic_night'] = dict(SBM_NIGHT_RECORD)
    return dict(scheme='sbm', nprc=2, **s)


def convection_namelist(c: dict) -> list:
    """The rainmod namelist entries (name, text) the scheme's record sets: its settings, the CAPE trigger's and the
    cloud-depth patch's only when the record holds them, the depth in Pa as PlaSim reads it; for the mass-flux scheme
    every parameter of its record."""
    entries = [('NPRC', str(c['nprc']))]
    if c.get('scheme') == 'massflux':
        start = c['start']['parameters'] if 'start' in c else {}
        return entries + [(name.upper(), mff.fortran_literal(value))
                          for name, value in list(c['parameters'].items()) + list(start.items())]
    for key, name in (('tau_s', 'TAU_BM'), ('rh', 'RH_BM'), ('eps', 'EPS_BM'), ('cin_j_kg', 'CIN_BM'),
                      ('cape_j_kg', 'CAPE_BM'), ('trigger_eps', 'EPST_BM')):
        if key in c:                                                  # the CAPE trigger only when set
            entries.append((name, repr(float(c[key]))))
    if 'deep_hpa' in c:
        entries.append(('DPDEEP_BM', repr(float(c['deep_hpa']) * 100.0)))
    if 'mixing_steps' in c:
        entries.append(('NSUB_BM', str(int(c['mixing_steps']))))
    if c.get('shallow_cap'):
        entries.append(('NSHCAP_BM', '1'))
    if c.get('diagnostic_night'):
        entries.append(('NDIAGNIGHT_BM', '1'))                       # the night switch, a diagnostic
    return entries


def uses_depth_patch(convection: dict | None) -> bool:
    """Whether a run's convection settings (as given to configuration) set an option of SBM_DEPTH_PATCHES."""
    return bool(convection) and ((convection.get('deep_hpa') or 0.0) > 0.0
                                 or (convection.get('mixing_steps') or 1) > 1)


def uses_massflux(convection: dict | None) -> bool:
    """Whether a run's convection settings (as given to configuration) choose the mass-flux scheme."""
    return bool(convection) and convection.get('scheme') == 'massflux'


def uses_massflux_start(convection: dict | None) -> bool:
    """Whether a run's convection settings (as given to configuration) set the mass-flux scheme's starting speed."""
    return uses_massflux(convection) and convection.get('start') not in (None, 'ifs')


def uses_cap_patch(convection: dict | None) -> bool:
    """Whether a run's convection settings (as given to configuration) set the shallow cap of SBM_CAP_PATCHES."""
    return bool(convection) and bool(convection.get('shallow_cap'))


def uses_night_switch(convection: dict | None) -> bool:
    """Whether a run's convection settings (as given to configuration) set the night switch of SBM_NIGHT_PATCHES, a
    diagnostic."""
    return bool(convection) and bool(convection.get('diagnostic_night'))


def cold_cloud_factors(path, half=None) -> dict:
    """The central factors of a cold_cloud.py product, its top layer first, as the configuration records them: the
    ice optics' water and crystal radius on PlaSim's top levels (ZICEWF, ZICERF), the levels below left at 1. With a
    run's half levels (half) it refuses a product whose layers are not the run's top layers."""
    path = Path(path)
    product = json.loads(path.read_text())
    if product.get('schema') != 'terluna.climate.gcm-cold-cloud-factors/1':
        raise ValueError(f'{path} is not a cold_cloud.py product')
    layers = sorted(product['layers'], key=lambda r: r['layer'])
    if [r['layer'] for r in layers] != list(range(1, len(layers) + 1)):
        raise ValueError(f'{path} does not hold the top layers from the first')
    factors = [(float(r['factors']['water']), float(r['factors']['radius'])) for r in layers]
    if half is not None:
        edges = np.r_[0.0, np.asarray(half, float)]
        for r in layers:
            k = r['layer']
            if abs(r['sigma_top'] - edges[k - 1]) > 1e-5 or abs(r['sigma_bottom'] - edges[k]) > 1e-5:
                raise ValueError(f"{path}: its layer {k} spans sigma {r['sigma_top']:.5f}-{r['sigma_bottom']:.5f}, the "
                                 f"run's {edges[k - 1]:.5f}-{edges[k]:.5f}; its factors belong to other levels")
    if not all(0.0 < w < 1e3 and 0.0 < r < 1e2 for w, r in factors):
        raise ValueError(f'{path} holds factors outside 0-1000 for the water or 0-100 for the radius')
    resolved = path.resolve()
    return dict(product=str(resolved.relative_to(ROOT)) if resolved.is_relative_to(ROOT) else str(resolved),
                sha256=hashlib.sha256(path.read_bytes()).hexdigest()[:16], run=product['run'],
                water_factors=[round(w, 6) for w, _ in factors], radius_factors=[round(r, 6) for _, r in factors])


def configuration(name, ncpus, overrides=None, co2_ppm=None, top_infrared=None, twilight='full', gas_bands=True,
                  heating=False, longwave_correction=False, ice_optics=False, twilight_ledger='301k',
                  twilight_clear_share=True, longwave_fit='301k', wind10=False, ice_weights='shield',
                  cold_cloud=None, rain_evaporation=None, moisture=False, convection=None, levels=None,
                  free_convection=None, mixing_length=None, boundary_layer=None, large_scale_condensate=None,
                  humidity_mean=False, surface_gust=None):
    """Everything that determines the run apart from the input files the runner writes. co2_ppm, when given, replaces the
    air's 400 ppm of CO2, top_infrared, (day, night) W/m2, adds infrared from space at the top of the air, and twilight
    (TWILIGHT_MODES) adds the curved air's sunlight, its term, the Moon's slant path or both, from the ledger
    twilight_ledger names (TWILIGHT_LEDGERS), its terms in the clear share above each layer with twilight_clear_share
    (TWILIGHT_CLEAR_EDITS), and gas_bands CO2 and O2 in PlaSim's sunlight (GAS_EDITS); heating writes each process's
    temperature tendency (HEATING), longwave_correction the longwave by height (LW_CORRECTION_EDITS) from the fit
    longwave_fit names (LONGWAVE_FITS), and ice_optics ice optics for cloud colder than -38 °C (ICE_OPTICS_EDITS) with
    the near-infrared weights ice_weights names (ICE_WEIGHTS_EDITS), and wind10 the 10 m wind (WIND10_PATCHES);
    cold_cloud, a cold_cloud.py product, scales the ice optics' water and radius on its layers by its central factors
    (ICE_FACTOR_EDITS); rain_evaporation replaces PlaSim's constant for the evaporation of falling rain and snow
    (rainmod's gamma); levels, a set of levels.LEVEL_SETS other than PlaSim's own, runs PlaSim on it, with the twilight
    terms and the longwave correction from their products on that set (levels_product) and the cold cloud's factors
    checked against its top layers; free_convection 'lunar' sets the free-convection constant over water for the
    planet's gravity and surface pressure (free_convection_setting, FREE_CONVECTION_EDITS); mixing_length 'lunar'
    stretches vdiff's Blackadar length by the gravities (mixing_length_setting) and boundary_layer 'nonlocal' runs the
    nonlocal boundary layer (boundary_layer_setting); with either of them, or with the Betts-Miller scheme's night
    switch (a diagnostic, SBM_NIGHT_PATCHES), the moisture diagnostics also write the vapour by process (BUDGET);
    large_scale_condensate 'cloud' holds large-scale condensate as cloud water (condensate.py; with the moisture
    diagnostics its own diagnostics and the vapour by process too), and humidity_mean writes the humidity's true mean
    (HUS_MEAN_PATCHES), and surface_gust 'beljaars' Beljaars's gust in the surface transfer in place of the
    free-convection term and the least wind (surface_gust_setting; with the moisture diagnostics the vapour by process
    too). Each is recorded when given, and left out the configuration is the one every earlier run had;
    the 294.9 K ledger without the clear share records the twilight as the runs of 14 to 19 October did, and the 295 K
    fit the longwave correction as those of 17 October did."""
    from importlib.metadata import version
    exp = EXPERIMENTS[name]
    planet = json.loads((PRODUCTS / 'moon_gcm_configuration.json').read_text())['planet']
    dry_bar = exp['pressure_pa'] / 1e5
    if co2_ppm is not None and not 0.0 < co2_ppm < 1e5:
        raise ValueError('co2_ppm must lie between 0 and 100,000')
    if top_infrared is not None and (len(top_infrared) != 2 or not all(0.0 <= x < 1000.0 for x in top_infrared)):
        raise ValueError('top_infrared is (day, night) W/m2, each from 0 to 1000')
    co2 = (CO2_PPM if co2_ppm is None else co2_ppm) * 1e-6 * dry_bar
    n2_ar = dry_bar - O2_PA / 1e5 - co2
    model = {**MODEL, **(overrides or {})}
    level_set, half = level_settings(model, overrides, levels)
    cloud_water_namelist(model['cloud_water'], planet['gravity_m_s2'])        # refuses an unknown choice
    if model['convective_day_s'] < 0:
        raise ValueError('convective_day_s must be 0 (the solar day) or a positive number of seconds')
    if not 0.0 <= model['surface_cloud_rcrit'] < 1.0:
        raise ValueError('surface_cloud_rcrit must be 0 (PlaSim\'s own) or lie below 1')
    if not 0.0 < model['sunlight_scale'] <= 1.5:
        raise ValueError('sunlight_scale must lie between 0 and 1.5')
    if model['sun_clock'] not in SUN_CLOCKS:
        raise ValueError(f'sun_clock must be one of {SUN_CLOCKS}, not {model["sun_clock"]!r}')
    if model['fixed_orbit'] not in (0, 1) or model['whole_lunar_days'] not in (0, 1):
        raise ValueError('fixed_orbit and whole_lunar_days are 0 or 1')
    if not 0.0 < model['sw_water_absorption'] <= 3.0:
        raise ValueError('sw_water_absorption must lie between 0 and 3')
    if twilight is not None and twilight not in TWILIGHT_MODES:
        raise ValueError(f'twilight must be one of {TWILIGHT_MODES}')
    if twilight_ledger not in TWILIGHT_LEDGERS:
        raise ValueError(f'twilight_ledger must be one of {tuple(TWILIGHT_LEDGERS)}')
    if longwave_fit not in LONGWAVE_FITS:
        raise ValueError(f'longwave_fit must be one of {tuple(LONGWAVE_FITS)}')
    if large_scale_condensate not in (None, 'cloud'):
        raise ValueError("large_scale_condensate is 'cloud' or None (PlaSim's own: rain at once)")
    cloud = large_scale_condensate == 'cloud'
    if surface_gust is not None and free_convection is not None:
        raise ValueError("the gust drops the free-convection term over water; set one of them")
    cfg = dict(experiment=name, **exp, model=model, output=output_variables(model), ncpus=ncpus,
               exoplasim=version('exoplasim'),
               source=ensure_patched(True, trigger=bool(convection and (convection.get('cape_j_kg') or 0.0) > 0.0),
                                     depth=uses_depth_patch(convection), cap=uses_cap_patch(convection),
                                     massflux=uses_massflux(convection),
                                     free_convection=free_convection is not None,
                                     massflux_start=uses_massflux_start(convection),
                                     boundary_layer=mixing_length is not None or boundary_layer is not None,
                                     night=uses_night_switch(convection), humidity_mean=bool(humidity_mean),
                                     cloud_water=cloud, gust=surface_gust is not None),
               planet=dict(radius=planet['radius_m'] / EARTH_RADIUS_M, gravity=planet['gravity_m_s2'],
                           rotationperiod=planet['rotation_period_days'], year=365.25636,
                           obliquity=planet['obliquity_deg'], eccentricity=0.016715),
               gases_bar=dict(pN2=n2_ar / (1 + AR_PER_N2), pO2=O2_PA / 1e5, pAr=n2_ar * AR_PER_N2 / (1 + AR_PER_N2),
                              pCO2=co2))
    if level_set != lv.DEFAULT:
        cfg['levels'] = dict(name=level_set, half_levels=[float(x) for x in half],
                             full_levels=[float(x) for x in lv.full_levels(half)],
                             plasim_namelist=dict(lv.namelist_entries(level_set)),
                             set_by='climate/gcm/levels.py (LEVEL_SETS)')
    if co2_ppm is not None:
        cfg['co2_ppm'] = co2_ppm
    if top_infrared is not None:
        cfg['top_infrared_w_m2'] = dict(day=float(top_infrared[0]), night=float(top_infrared[1]))
    if longwave_correction:
        fit = longwave_product(longwave_fit, level_set, half.size)
        product = json.loads(fit.read_text())
        cfg['longwave_correction'] = dict(product=str(fit.relative_to(ROOT)),
                                          sha256=hashlib.sha256(fit.read_bytes()).hexdigest()[:16],
                                          **{k: product['correction'][k] for k in ('land', 'sea', 'tref')})
    if ice_optics:
        cfg['ice_optics'] = dict(below_k=pio.ICE_BELOW_K, near_ir_weights=list(pio.NIR_WEIGHTS),
                                 radius='Kristjansson and Mitchell, as CAM3 and RRTMG carry it; its trend below 180 K',
                                 optics='Ebert and Curry (1992), as RRTMG carries them (iceflag 1 and 0)')
        if ice_weights == 'shield':
            cfg['ice_optics'].update(near_ir_weights=list(pio.SHIELD_NIR_WEIGHTS),
                                     weights_source="the design shield's sunlight, written as ZICEW (world-crafting's "
                                                    'patch 11, lab branch 79ea224)')
    if rain_evaporation is not None:
        cfg['rain_evaporation'] = rain_evaporation_setting(rain_evaporation)
    if free_convection is not None:
        cfg['free_convection'] = free_convection_setting(free_convection, planet['gravity_m_s2'], exp['pressure_pa'])
    if mixing_length is not None:
        cfg['mixing_length'] = mixing_length_setting(mixing_length, planet['gravity_m_s2'])
    if boundary_layer is not None:
        cfg['boundary_layer'] = boundary_layer_setting(boundary_layer)
    if surface_gust is not None:
        cfg['surface_gust'] = surface_gust_setting(surface_gust)
    if cloud:
        cfg['large_scale_condensate'] = cnd.setting()
    if humidity_mean:
        cfg['humidity_mean'] = dict(code=HUS_MEAN_CODE, name=HUS_MEAN[0],
                                    field='the gridpoint humidity PlaSim adds every step (outmod aadq)',
                                    patch='climate/gcm/exoplasim_run.py HUS_MEAN_PATCHES (pumamod nhusmean)')
        cfg['output'] = cfg['output'] + [HUS_MEAN[0]]
    if convection is not None:                       # the Betts-Miller or the mass-flux scheme in Kuo's place
        cfg['convection'] = convection_setting(**convection)
    if cold_cloud is not None:
        if not ice_optics:
            raise ValueError("the cold cloud's factors act on the ice optics, which need --ice-optics")
        cfg['ice_optics']['cold_cloud'] = cold_cloud_factors(cold_cloud, half=half)
    if wind10:
        cfg['wind10'] = dict(codes=[165, 166],
                             profile="Geleyn's (1988) profile through the lowest layer, from the drag PlaSim uses with "
                                     'its stability factor (WIND10_PATCHES)')
        cfg['output'] = cfg['output'] + ['uas', 'vas']
    if heating:
        cfg['heating_diagnostics'] = True
        cfg['output'] = cfg['output'] + [name for name, _ in HEATING.values()]
    if moisture:
        if not heating:
            raise ValueError('the moisture diagnostics come with the heating diagnostics, which need --heating-diagnostics')
        cfg['moisture_diagnostics'] = True
        cfg['output'] = cfg['output'] + [name for name, _ in MOISTENING.values()] + ['wap']
        if (mixing_length is not None or boundary_layer is not None or uses_night_switch(convection) or cloud
                or surface_gust is not None):
            cfg['budget_diagnostics'] = sorted(BUDGET)
            cfg['output'] = cfg['output'] + [name for name, _, _ in BUDGET.values()]
        if cloud:
            cfg['cloud_water_diagnostics'] = sorted(cnd.DIAGNOSTICS)
            cfg['output'] = cfg['output'] + [name for name, _, _ in cnd.DIAGNOSTICS.values()]
    if twilight is not None:
        ledger = TWILIGHT_LEDGERS[twilight_ledger]
        terms = ledger if level_set == lv.DEFAULT else levels_product(ledger, level_set)
        lunar = json.loads(SLANT_PRODUCT.read_text())['magnification']
        cfg['twilight'] = dict(mode=twilight, product=str(terms.relative_to(ROOT)),
                               sha256=hashlib.sha256(terms.read_bytes()).hexdigest()[:16],
                               magnification=[lunar['a'], lunar['b'], lunar['c']])
        if terms != ledger:
            cfg['twilight']['ledger'] = str(ledger.relative_to(ROOT))
        if ledger != SLANT_PRODUCT:
            cfg['twilight'].update(slant_product=str(SLANT_PRODUCT.relative_to(ROOT)),
                                   slant_sha256=hashlib.sha256(SLANT_PRODUCT.read_bytes()).hexdigest()[:16])
        if twilight_clear_share and twilight in ('curved', 'full'):
            cfg['twilight']['clear_share'] = ('each layer\'s sunlight in the cover clear above it, its infrared and the '
                                              'ground\'s through the clouds above at their longwave transmissivity')
    if gas_bands:
        fits = json.loads(GASES_PRODUCT.read_text())
        air = json.loads(SLANT_PRODUCT.read_text())['magnification_air']
        cfg['gas_bands'] = dict(product=str(GASES_PRODUCT.relative_to(ROOT)),
                                sha256=hashlib.sha256(GASES_PRODUCT.read_bytes()).hexdigest()[:16],
                                fits=[fits['fits'][k] for k in GAS_ORDER],
                                mass_fractions=gas_mass_fractions(cfg['gases_bar']),
                                magnification_air=[air['a'], air['b'], air['c']])
    return cfg


def level_settings(model: dict, overrides: dict | None, levels: str | None) -> tuple:
    """The level set a run takes (levels.LEVEL_SETS, PlaSim's own when None) and its half levels, setting the model's
    layer count to the set's; a layer count other than the set's is refused, so per-layer inputs always match."""
    level_set = levels or lv.DEFAULT
    half = lv.half_levels(level_set)                                           # refuses an unknown set
    if level_set == lv.DEFAULT and model['layers'] != half.size:
        raise ValueError(f"{model['layers']} layers need a level set of as many (levels); PlaSim's own has {half.size}")
    if level_set != lv.DEFAULT:
        if (overrides or {}).get('layers', half.size) != half.size:
            raise ValueError(f"level set {level_set} has {half.size} layers, not {overrides['layers']}")
        model['layers'] = half.size
    return level_set, half


def levels_product(product: Path, level_set: str) -> Path:
    """A 10-level input product's counterpart on another level set, written beside it as <name>_<set>.json by the
    script that wrote it: the twilight ledger's GCM terms (twilight.py levels) and the longwave correction
    (plasim_longwave.py levels)."""
    path = product.with_name(f'{product.stem}_{level_set}.json')
    if not path.exists():
        raise FileNotFoundError(f'{path.relative_to(ROOT)} is missing: the script that wrote {product.name} writes it for '
                                f'level set {level_set} (twilight.py levels, plasim_longwave.py levels)')
    return path


def longwave_product(longwave_fit: str, level_set: str, layers: int) -> Path:
    """The longwave correction a run on level_set reads (LONGWAVE_FITS): a fit made on that set's own columns as it
    is, a fit on PlaSim's ten levels carried onto the set (levels_product); a fit made on another set is refused, and
    so is a correction without a value for each of the run's layers and the ground."""
    fit = LONGWAVE_FITS[longwave_fit]
    product = json.loads(fit.read_text())
    own = product.get('levels', {}).get('name', lv.DEFAULT)                  # the set the fit was made on
    if own != level_set and own != lv.DEFAULT:
        raise ValueError(f'{fit.name} was fitted on level set {own}, and the run takes level set {level_set}')
    path = fit if own == level_set else levels_product(fit, level_set)
    correction = json.loads(path.read_text())['correction']
    if any(len(correction[k]) != layers + 1 for k in ('land', 'sea', 'tref')):
        raise ValueError(f'{path.name} does not hold a value for each of {layers} layers and the ground')
    return path


def gas_mass_fractions(gases: dict) -> dict:
    """The mass fractions of CO2 and O2 (GAS_EDITS' ZGASW) in air of these partial pressures (bar)."""
    total = sum(gases.values())
    molar = sum(gases[k] / total * MOLAR[k] for k in gases)
    return dict(CO2=gases['pCO2'] / total * MOLAR['pCO2'] / molar, O2=gases['pO2'] / total * MOLAR['pO2'] / molar)


def twilight_table(path: Path, mode: str, ledger: Path | None = None):
    """A twilight ledger's GCM terms (TWILIGHT_PRODUCT without one) as PlaSim's patch reads them: the number of
    elevations, then for each the solar elevation (degrees), what the curved air adds to the ground and to each layer,
    top first, per unit sunlight at the beam; then the extra infrared each layer and the ground lose to space (W/m2)
    and the temperatures it is given at. 'full' adds the shells' widening to both."""
    product = json.loads((ledger or TWILIGHT_PRODUCT).read_text())
    term, cone, lw = product['gcm_term'], product['gcm_term_cone'], product['longwave']
    full = mode == 'full'
    rows = [f"{len(term['elevation_deg'])}"]
    for i, el in enumerate(term['elevation_deg']):
        layers = [q[i] + (c[i] if full else 0.0) for q, c in zip(term['layers'], cone['layers'])]
        rows.append(' '.join(f'{x:.6e}' for x in [el, term['ground'][i]] + layers))
    escape = [a + (b if full else 0.0) for a, b in zip(lw['gcm_layers']['curved'], lw['gcm_layers']['cone'])]
    rows.append(' '.join(f'{x:.6e}' for x in escape + [lw['ground_curved_w_m2']]))
    rows.append(' '.join(f'{x:.6e}' for x in lw['gcm_layers']['temperature_k'] + [lw['ground_temperature_k']]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(rows) + '\n')
    return path


def run_key(cfg, inputs) -> str:
    """Changes if the settings or any written input file change, so a resume never mixes models."""
    digest = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode())
    for path in (inputs['landmap'], inputs['topomap'], inputs['starspec'], inputs['starspec'][:-4] + '_hr.dat',
                 *([inputs['twilight']] if 'twilight' in inputs else [])):
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


def surface_product(water, lakes):
    """The T21 surface an experiment reads: without lakes the seas alone; lakes True, the surface of 2026-09-25 with run
    A's lakes, kept as it is (boundary.SUPERSEDED); lakes a drainage product's content hash, the surface built from
    that drainage (boundary.lakes_name)."""
    if lakes is True:
        return f'moon_{water}pct_water_lakes_gaussian_T21.nc'
    if lakes:
        return f'moon_{water}pct_water_lakes_{lakes}_gaussian_T21.nc'
    return f'moon_{water}pct_water_gaussian_T21.nc'


def surface_files(cfg, folder: Path):
    """Land mask (code 172) and surface geopotential (code 129) at T21 from the geography product. With
    lakes, the water share is the product's (seas and lakes) and water cells stand at their water surface."""
    import netCDF4
    lakes = cfg.get('lakes', False)
    name = surface_product(cfg['water'], lakes)
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
    out = dict(landmap=str(landmap), topomap=str(topomap), starspec=starspec,
               flux_w_m2=flux * cfg['model']['sunlight_scale'], shield_flux_w_m2=flux, land_share=land_share)
    if cfg.get('twilight', {}).get('mode') in ('curved', 'full'):
        out['twilight'] = str(twilight_table(folder / 'twilight.dat', cfg['twilight']['mode'],
                                             ROOT / cfg['twilight']['product']))
    return out


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
                    eccentricity=p['eccentricity'], obliquity=p['obliquity'], fixedorbit=bool(m['fixed_orbit']),
                    landmap=str(landmap), topomap=str(topomap), mldepth=cfg['mldepth'], seaice=True, ozone=False,
                    soilalbedo=m['land_albedo'], timestep=m['timestep_min'], snapshots=None,
                    runsteps=m['steps_per_lunar_day'] * m['lunar_days_per_year'] if m['whole_lunar_days'] else None,
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
    for name, text in lv.namelist_entries(cfg.get('levels', {}).get('name', lv.DEFAULT)):
        model._edit_namelist('plasim_namelist', name, text)               # a level set's half levels and top friction
    model._edit_namelist('radmod_namelist', 'RAYSCALE', str(m['rayleigh_scale']))
    model._edit_namelist('radmod_namelist', 'NSUNCLOCK', str(SUN_CLOCKS.index(m['sun_clock'])))
    model._edit_namelist('radmod_namelist', 'SWH2OSCALE', str(m['sw_water_absorption']))
    model._edit_namelist('radmod_namelist', 'TH2OC', str(m['h2o_continuum']))
    gref, scale = cloud_water_namelist(m['cloud_water'], p['gravity'])
    model._edit_namelist('rainmod_namelist', 'CWGREF', str(gref))
    model._edit_namelist('rainmod_namelist', 'CWSCALE', str(scale))
    model._edit_namelist('plasim_namelist', 'NDIAGCF', str(m['clear_sky']))
    model._edit_namelist('rainmod_namelist', 'CONVDAY', str(float(m['convective_day_s'])))
    if 'rain_evaporation' in cfg:
        model._edit_namelist('rainmod_namelist', 'GAMMA', repr(float(cfg['rain_evaporation']['gamma'])))
    if 'convection' in cfg:
        for name, text in convection_namelist(cfg['convection']):
            model._edit_namelist('rainmod_namelist', name, text)
    if 'large_scale_condensate' in cfg:                               # condensate.py's cloud water
        model._edit_namelist('plasim_namelist', 'NCWAT', '1')
        for name, text in cnd.namelist_values(cfg['large_scale_condensate']['parameters']).items():
            model._edit_namelist('rainmod_namelist', name, text)
    if 'humidity_mean' in cfg:
        model._edit_namelist('plasim_namelist', 'NHUSMEAN', '1')
    if m['surface_cloud_rcrit'] > 0.0:                                 # the lowest layer is the last of RCRIT
        model._edit_namelist('rainmod_namelist', f"RCRIT({m['layers']})", str(float(m['surface_cloud_rcrit'])))
    if 'top_infrared_w_m2' in cfg:
        model._edit_namelist('radmod_namelist', 'IRTOPDAY', str(cfg['top_infrared_w_m2']['day']))
        model._edit_namelist('radmod_namelist', 'IRTOPNIGHT', str(cfg['top_infrared_w_m2']['night']))
    if 'twilight' in inputs:
        shutil.copy2(inputs['twilight'], workdir / 'twilight.dat')
        model._edit_namelist('radmod_namelist', 'NTWILIGHT', '1')
        model._edit_namelist('radmod_namelist', 'NTWILW', '1')
        if cfg['twilight'].get('clear_share'):
            model._edit_namelist('radmod_namelist', 'NTWLCLEAR', '1')
    if 'longwave_correction' in cfg:
        c = cfg['longwave_correction']
        (workdir / 'lwcorr.dat').write_text(''.join(' '.join(f'{x:.6e}' for x in c[k]) + '\n' for k in ('land', 'sea', 'tref')))
        model._edit_namelist('radmod_namelist', 'NLWCORR', '1')
    if 'ice_optics' in cfg:
        model._edit_namelist('radmod_namelist', 'NICEOPT', '1')
        model._edit_namelist('radmod_namelist', 'TICEOPT', str(cfg['ice_optics']['below_k']))
        if cfg['ice_optics'].get('weights_source'):
            for i, w in enumerate(cfg['ice_optics']['near_ir_weights'], 1):
                model._edit_namelist('radmod_namelist', f'ZICEW({i})', repr(float(w)))
        factors = cfg['ice_optics'].get('cold_cloud')
        if factors:
            for i, (w, r) in enumerate(zip(factors['water_factors'], factors['radius_factors']), 1):
                model._edit_namelist('radmod_namelist', f'ZICEWF({i})', repr(float(w)))
                model._edit_namelist('radmod_namelist', f'ZICERF({i})', repr(float(r)))
    if 'wind10' in cfg:
        model._edit_namelist('fluxmod_namelist', 'NWIND10', '1')
    if 'free_convection' in cfg:
        model._edit_namelist('fluxmod_namelist', 'CFREE', repr(float(cfg['free_convection']['constant'])))
    if 'mixing_length' in cfg:
        model._edit_namelist('fluxmod_namelist', 'VDIFF_LAMM', repr(float(cfg['mixing_length']['momentum_m'])))
    if 'boundary_layer' in cfg:
        model._edit_namelist('fluxmod_namelist', 'NLPBL', '1')
        for name, value in cfg['boundary_layer']['parameters'].items():
            model._edit_namelist('fluxmod_namelist', name.upper(), repr(float(value)))
    if 'surface_gust' in cfg:
        model._edit_namelist('fluxmod_namelist', 'NGUST', '1')
        model._edit_namelist('fluxmod_namelist', 'GUST_BETA', repr(float(cfg['surface_gust']['beta'])))
        model._edit_namelist('fluxmod_namelist', 'GUST_UMIN', repr(float(cfg['surface_gust']['least_wind_m_s'])))
    if cfg.get('twilight', {}).get('magnification'):
        for i, x in enumerate(cfg['twilight']['magnification'], 1):
            model._edit_namelist('radmod_namelist', f'ZMLUNAR({i})', repr(float(x)))
    if 'gas_bands' in cfg:
        g = cfg['gas_bands']
        model._edit_namelist('radmod_namelist', 'NGASBANDS', '1')
        for j, q in enumerate(g['fits'], 1):
            for i, x in enumerate(q, 1):
                model._edit_namelist('radmod_namelist', f'ZGASQ({i},{j})', repr(float(x)))
        model._edit_namelist('radmod_namelist', 'ZGASW(1)', repr(float(g['mass_fractions']['CO2'])))
        model._edit_namelist('radmod_namelist', 'ZGASW(2)', repr(float(g['mass_fractions']['O2'])))
        for i, x in enumerate(g['magnification_air'], 1):
            model._edit_namelist('radmod_namelist', f'ZMAIRLUNAR({i})', repr(float(x)))
    from exoplasim import pyburn
    for code, (name, long_name) in CLEAR_SKY.items():
        pyburn.ilibrary[str(code)] = [name, long_name, 'W m-2']
        pyburn.slibrary[name] = [code, long_name, 'W m-2']
    if cfg.get('heating_diagnostics'):
        model._edit_namelist('plasim_namelist', 'NDIAGGP', '1')
        for code, (name, long_name) in HEATING.items():
            pyburn.ilibrary[str(code)] = [name, long_name, 'K s-1']
            pyburn.slibrary[name] = [code, long_name, 'K s-1']
    if cfg.get('moisture_diagnostics'):
        for code, (name, long_name) in MOISTENING.items():
            pyburn.ilibrary[str(code)] = [name, long_name, 'm s-1']
            pyburn.slibrary[name] = [code, long_name, 'm s-1']
    if cfg.get('budget_diagnostics'):
        for code, (name, long_name, units) in BUDGET.items():
            pyburn.ilibrary[str(code)] = [name, long_name, units]
            pyburn.slibrary[name] = [code, long_name, units]
    if cfg.get('cloud_water_diagnostics'):                            # four slots beyond PlaSim's 21
        model._edit_namelist('plasim_namelist', 'NDIAGGP3D', str(cnd.EXTRA_SLOTS))
        for code, (name, long_name, units) in cnd.DIAGNOSTICS.items():
            pyburn.ilibrary[str(code)] = [name, long_name, units]
            pyburn.slibrary[name] = [code, long_name, units]
    if 'humidity_mean' in cfg:
        pyburn.ilibrary[str(HUS_MEAN_CODE)] = list(HUS_MEAN)
        pyburn.slibrary[HUS_MEAN[0]] = [HUS_MEAN_CODE, HUS_MEAN[1], HUS_MEAN[2]]
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

def summarise_year(path: Path, times_per_year: int | None = None) -> dict:
    """Headline numbers from one year of output, for watching the run settle at a glance. The year holds
    times_per_year means, the runner's 120 (about 3 days each) unless the run writes more often."""
    import netCDF4
    with netCDF4.Dataset(path) as d:
        clear = 'rstcs' in d.variables
        keys = ('ts', 'tas', 'ntr', 'clt', 'pr', 'sic', 'lsm', 'rst', 'rsut', 'rlut') + (('rstcs', 'rltcs') if clear else ())
        v = {k: np.asarray(d[k][:], dtype=float) for k in keys}
        lat = np.asarray(d['lat'][:], dtype=float)
        times = d['time'].size
    expected = times_per_year or MODEL['writes_per_lunar_day'] * MODEL['lunar_days_per_year']
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
        overrides: dict | None = None, start_from: str | None = None, co2_ppm: float | None = None,
        top_infrared: tuple | None = None, twilight: str | None = 'full', gas_bands: bool = True,
        heating: bool = False, longwave_correction: bool = False, ice_optics: bool = False,
        twilight_ledger: str = '301k', twilight_clear_share: bool = True, longwave_fit: str = '301k',
        wind10: bool = False, ice_weights: str = 'shield', cold_cloud: str | None = None,
        rain_evaporation: float | None = None, moisture: bool = False, convection: dict | None = None,
        levels: str | None = None, free_convection: str | None = None, mixing_length: str | None = None,
        boundary_layer: str | None = None, large_scale_condensate: str | None = None,
        humidity_mean: bool = False, surface_gust: str | None = None) -> int:
    use_twilight_package()
    cfg = configuration(name, ncpus, overrides, co2_ppm, top_infrared, twilight, gas_bands, heating, longwave_correction,
                        ice_optics, twilight_ledger, twilight_clear_share, longwave_fit, wind10=wind10,
                        ice_weights=ice_weights, cold_cloud=cold_cloud, rain_evaporation=rain_evaporation,
                        moisture=moisture, convection=convection, levels=levels, free_convection=free_convection,
                        mixing_length=mixing_length, boundary_layer=boundary_layer,
                        large_scale_condensate=large_scale_condensate, humidity_mean=humidity_mean,
                        surface_gust=surface_gust)
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
            writes = cfg['model']['writes_per_lunar_day'] * cfg['model']['lunar_days_per_year']
            entry = dict(year=year, seconds=round(time.time() - start, 1), **summarise_year(out, writes))
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
    parser.add_argument('--co2-ppm', type=float, help='CO2 the air carries, ppm, in place of 400 (part of the run key)')
    parser.add_argument('--top-infrared', type=float, nargs=2, metavar=('DAY', 'NIGHT'),
                        help='infrared from space at the top of the air, W/m2 by day and by night '
                             '(part of the run key)')
    parser.add_argument('--gas-bands', action=argparse.BooleanOptionalAction, default=True,
                        help="CO2 and O2 in PlaSim's sunlight from ../results/gcm/shortwave_gases.json, the standard "
                             'since 2026-10-14 (part of the run key)')
    parser.add_argument('--longwave-correction', action='store_true',
                        help='the longwave by height, fitted to the correlated-k model (LONGWAVE_FITS; part of the '
                             'run key)')
    parser.add_argument('--longwave-fit', choices=tuple(LONGWAVE_FITS), default='301k',
                        help='which fit --longwave-correction adds: the refit on the 301 K columns (the standard since '
                             '2026-10-19), the fit near 295 K the runs of 17 October took, or l15_provisional, refitted '
                             'on the 15-level month trial for runs on l15 (part of the run key)')
    parser.add_argument('--ice-optics', action='store_true',
                        help='ice optics for cloud colder than -38 C (plasim_ice_optics.py; part of the run key)')
    parser.add_argument('--ice-weights', choices=('shield', 'stock'), default='shield',
                        help="the ice optics' near-infrared weights: the shares of the design shield's sunlight "
                             "(ZICEW, the standard since 2026-10-19) or the unfiltered Sun's the earlier runs took "
                             '(part of the run key)')
    parser.add_argument('--cold-cloud', metavar='PRODUCT',
                        help="a cold_cloud.py product whose central factors scale the ice optics' water and crystal "
                             'radius on its layers, carrying them to lunar gravity (ICE_FACTOR_EDITS; needs '
                             '--ice-optics; part of the run key)')
    parser.add_argument('--convection', choices=('kuo', 'sbm', 'massflux'), default='kuo',
                        help="PlaSim's deep convection: Kuo's scheme, its own, the simplified Betts-Miller scheme of "
                             'the Terluna patch (SBM_EDITS), or the mass-flux scheme (massflux_fortran.py; part of the '
                             'run key)')
    parser.add_argument('--massflux-start', choices=tuple(mff.START_CHOICES), default='ifs',
                        help="the mass-flux scheme's surface parcel's starting speed: the IFS's (the larger of 1 m/s "
                             "and the excess's velocity scale) or c_start times the boundary layer's convective velocity "
                             "w* (massflux_fortran.START_PATCHES; part of the run key)")
    parser.add_argument('--sbm-tau', type=float, metavar='SECONDS', help='its relaxation time (7200)')
    parser.add_argument('--sbm-rh', type=float, metavar='SHARE', help='its reference relative humidity (0.7)')
    parser.add_argument('--sbm-eps', type=float, help="its parcel's mixing per e-folding of pressure (0)")
    parser.add_argument('--sbm-cin', type=float, metavar='J_KG', help="its limit on the parcel's inhibition (0: none)")
    parser.add_argument('--sbm-cape', type=float, metavar='J_KG',
                        help='its CAPE trigger: deep convection only where a parcel has this CAPE (0: none; '
                             'SBM_TRIGGER_PATCHES)')
    parser.add_argument('--sbm-trigger-eps', type=float, metavar='EPS',
                        help="the CAPE trigger's parcel mixing per e-folding of pressure (0)")
    parser.add_argument('--sbm-deep-hpa', type=float, metavar='HPA',
                        help="its cloud-depth test: deep convection only where the parcel's cloud, from its "
                             'condensation level to its top, is thicker than this (0: none; SBM_DEPTH_PATCHES)')
    parser.add_argument('--sbm-mixing-steps', type=int, metavar='N',
                        help="the parcel's mixing within each layer in N steps of log pressure (1: one a layer; "
                             'needs --sbm-eps; SBM_DEPTH_PATCHES)')
    parser.add_argument('--sbm-diagnostic-night', action='store_true',
                        help='a diagnostic, never a configuration: no deep convection in the simplified Betts-Miller '
                             "scheme wherever the Sun is down (the geometric cosine of its zenith at the grid point at "
                             "most 0), such columns taking PlaSim's shallow scheme (SBM_NIGHT_PATCHES; with "
                             '--moisture-diagnostics it writes the vapour by process; part of the run key)')
    parser.add_argument('--sbm-shallow-cap', action='store_true',
                        help="with --sbm-deep-hpa, a shallow cloud's top no higher than that depth above its own base, "
                             "in place of PlaSim's 700 hPa (SBM_CAP_PATCHES; part of the run key)")
    parser.add_argument('--levels', choices=tuple(lv.LEVEL_SETS), default=lv.DEFAULT,
                        help="PlaSim's levels: its own ten, or a set of levels.py (l14: fourteen, more of them near "
                             'the ground; l15: those and a bottom layer 0.6 km thick), which needs a restart on the '
                             'same levels (levels.py restart; part of the run key)')
    parser.add_argument('--rain-evaporation', type=float, metavar='GAMMA',
                        help="PlaSim's constant for the evaporation of falling rain and snow (rainmod's gamma, "
                             f'{PLASIM_RAIN_GAMMA} its own; part of the run key)')
    parser.add_argument('--wind10', action='store_true',
                        help="the 10 m wind on Geleyn's profile from the lowest layer's drag, codes 165 and 166 "
                             '(WIND10_PATCHES; part of the run key)')
    parser.add_argument('--free-convection', choices=FREE_CONVECTION_CHOICES, default='earth',
                        help="the free-convection constant in the heat and vapour transfer over water (Miller and "
                             "others 1992): PlaSim's Earth value, 0.0016 m s-1 K-1/3, or the lunar one, scaled as "
                             "(g / p_s)^(1/3) to the Moon's gravity and the surface pressure (FREE_CONVECTION_EDITS; "
                             'part of the run key)')
    parser.add_argument('--mixing-length', choices=MIXING_LENGTH_CHOICES, default='earth',
                        help="Blackadar's asymptotic mixing length in PlaSim's vertical diffusion: its Earth value, "
                             '160 m for momentum, or that stretched by the ratio of the gravities, a sensitivity '
                             '(mixing_length_setting; part of the run key)')
    parser.add_argument('--boundary-layer', choices=BOUNDARY_LAYER_CHOICES, default='local',
                        help="PlaSim's local vertical diffusion alone, or the nonlocal boundary layer beneath it "
                             '(pbl_scheme.py; pbl_fortran.PBL_PATCHES; part of the run key)')
    parser.add_argument('--surface-gust', choices=SURFACE_GUST_CHOICES, default='none',
                        help="Beljaars's gust, the boundary layer's convective velocity, added to the wind in PlaSim's "
                             'surface transfer over every surface in place of the free-convection term over water and '
                             'the least wind of 1 m/s (surface_gust.py; GUST_PATCHES; part of the run key)')
    parser.add_argument('--large-scale-condensate', choices=('rain', 'cloud'), default='rain',
                        help="what PlaSim's large-scale condensation makes of its condensate: rain at once, as "
                             'PlaSim does, or cloud water held in its layer above freezing, evaporated where the air '
                             'is unsaturated and turned to rain by autoconversion and accretion (condensate.py; '
                             'CLOUD_WATER_PATCHES; part of the run key)')
    parser.add_argument('--humidity-mean', action='store_true',
                        help="write the humidity's true 3-day mean as hus_mean, code 233, where a year file's hus "
                             "holds the humidity at each interval's start (HUS_MEAN_PATCHES; year_files.py; part of "
                             'the run key)')
    parser.add_argument('--heating-diagnostics', action='store_true',
                        help="PlaSim's temperature tendency of each process at every level, codes 21-29 (part of "
                             'the run key)')
    parser.add_argument('--moisture-diagnostics', action='store_true',
                        help="with --heating-diagnostics, the water vapour each layer gains from the physics, the "
                             "negative-humidity fix and the dynamics, codes 31, 32, 40 and 41, and omega (part of the "
                             'run key)')
    parser.add_argument('--twilight', choices=TWILIGHT_MODES + ('none',), default='full',
                        help="the curved air's sunlight and infrared from the twilight ledger: the curved "
                             "parts, all of it (full, the standard since 2026-10-14), the Moon's slant path alone, or "
                             'none (part of the run key)')
    parser.add_argument('--twilight-ledger', choices=tuple(TWILIGHT_LEDGERS), default='301k',
                        help="the ledger the terms come from: the design climate's 301 K column (twilight_301k.json, "
                             'the standard since 2026-10-19) or the 294.9 K column (twilight.json) the runs of 14 to 19 '
                             'October took (part of the run key)')
    parser.add_argument('--twilight-clear-share', action=argparse.BooleanOptionalAction, default=True,
                        help="the twilight's terms in the clear share above each layer (TWILIGHT_CLEAR_EDITS), the "
                             'standard since 2026-10-19; --no-twilight-clear-share adds them in every column as the '
                             'earlier runs did (part of the run key)')
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
    if args.sbm_shallow_cap and args.convection != 'sbm':
        parser.error('--sbm-shallow-cap belongs to the simplified Betts-Miller scheme (--convection sbm)')
    if args.sbm_diagnostic_night and args.convection != 'sbm':
        parser.error('--sbm-diagnostic-night, a diagnostic, belongs to the simplified Betts-Miller scheme '
                     '(--convection sbm)')
    if args.massflux_start != 'ifs' and args.convection != 'massflux':
        parser.error('--massflux-start belongs to the mass-flux scheme (--convection massflux)')
    sbm_given = [k for k in ('sbm_tau', 'sbm_rh', 'sbm_eps', 'sbm_cin', 'sbm_cape', 'sbm_trigger_eps', 'sbm_deep_hpa',
                             'sbm_mixing_steps') if getattr(args, k) is not None]
    if sbm_given and args.convection != 'sbm':
        parser.error(f'{sbm_given} belong to the simplified Betts-Miller scheme (--convection sbm)')
    convection = None
    if args.convection == 'sbm':
        convection = dict(tau_s=args.sbm_tau, rh=args.sbm_rh, eps=args.sbm_eps, cin_j_kg=args.sbm_cin,
                          cape_j_kg=args.sbm_cape, trigger_eps=args.sbm_trigger_eps, deep_hpa=args.sbm_deep_hpa,
                          mixing_steps=args.sbm_mixing_steps, shallow_cap=args.sbm_shallow_cap or None,
                          diagnostic_night=args.sbm_diagnostic_night or None)
    elif args.convection == 'massflux':
        convection = dict(scheme='massflux', start=None if args.massflux_start == 'ifs' else args.massflux_start)
    return run(args.experiment, args.years, args.ncpus, args.folder, args.prune, overrides, args.start_from,
               args.co2_ppm, tuple(args.top_infrared) if args.top_infrared else None,
               None if args.twilight == 'none' else args.twilight, args.gas_bands, args.heating_diagnostics,
               args.longwave_correction, args.ice_optics, args.twilight_ledger, args.twilight_clear_share,
               args.longwave_fit, args.wind10, args.ice_weights, args.cold_cloud, args.rain_evaporation,
               args.moisture_diagnostics, convection, args.levels,
               None if args.free_convection == 'earth' else args.free_convection,
               None if args.mixing_length == 'earth' else args.mixing_length,
               None if args.boundary_layer == 'local' else args.boundary_layer,
               None if args.large_scale_condensate == 'rain' else args.large_scale_condensate, args.humidity_mean,
               None if args.surface_gust == 'none' else args.surface_gust)


if __name__ == '__main__':
    sys.exit(main())
