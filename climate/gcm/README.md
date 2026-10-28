# Three-dimensional climate: preparation

This folder holds what the general circulation model (GCM) experiments need from
the repository, the questions only such a model can answer, and the runner for
ExoPlaSim, the first model installed (2026-09-25).

| File | Holds |
|---|---|
| [boundary.py](boundary.py) | Builds the planet parameters and the surface files below from `shared/constants.json` and the geography product; the lakes surface is named by its drainage product's content, and its committed record is `../results/gcm/boundary_lakes_<drainage>.json` (`terluna.climate.gcm-boundary-lakes/1`) |
| [wave_coverage.py](wave_coverage.py) | Recovers years 20–29 of the corrected design case over Smythii–Marginis: 71.5-hour means at the lowest model level, with source hashes and explicit surface-forcing requirements ([coastal-wave study](../waves/coastal.md)) |
| [wave_snapshots.py](wave_snapshots.py) | Runs an isolated copy of the corrected year-29 restart, capturing global three-hourly winds and surface stress: 241 snapshots over 30 Earth days and 481 over 60 days. Matching controls preserve the final climate state and ordinary output exactly ([weather episode](../waves/weather.md), [full-cycle waves](../waves/cycle.md)) |
| `products/moon_{25,35}pct_water_{grid}.nc` | Land fraction, land elevation above the water level and ocean depth on regular 4°×5° and 2°×2.5° grids and on Gaussian T21 and T42 grids (NetCDF 3). These are area-conserving averages of the LOLA/GRAIL hydrostatic filling |
| `products/moon_28pct_water_lakes_{grid}.nc` | The 28% seas with run A's rain-fed lakes (25 September), which every A28 run read. Its drainage product was overwritten on 29 September, so it is kept as it is and never rebuilt ([the lakes for a cold start](#the-lakes-for-a-cold-start)) |
| `products/moon_28pct_water_lakes_98842fc2a648f889_{grid}.nc` | The same seas with the lakes of the current routing, from the drainage product with that sha256 prefix; experiment `A28R` reads it |
| `products/moon_gcm_configuration.json` | Radius, gravity, rotation (synchronous, 27.32 d) and solar day (29.53 d), obliquity 1.54°, solar constant, and pointers to the atmosphere, ozone and shield products |
| [exoplasim_run.py](exoplasim_run.py) | Runs the ExoPlaSim experiments one simulated year at a time, resumable after any stop |
| [levels.py](levels.py) | PlaSim's vertical levels: its own ten, the fourteen of `l14`, with more of them near the ground, and the fifteen of `l15`, which add a bottom layer 0.6 km thick; the namelist a run on them writes, and the conservative remap that carries a restart onto them (`levels.py restart`) |
| [cloud_water.py](cloud_water.py) | Recomputes from a year of output the cloud water PlaSim's radiation sees, layer by layer, as PlaSim has it and with the runner's corrections |
| [plasim_longwave.py](plasim_longwave.py) | PlaSim's clear-sky longwave column by column, checked against its own heating, and the correction by height fitted to the correlated-k model that the runner's `--longwave-correction` adds (`../results/gcm/longwave_<run>.json`); its `check` sets the correction against correlated-k on states from 295 to 306.6 K (`../results/gcm/longwave_check_301k.json`), its `levels` carries a ten-level fit onto another level set, and its `refit` refits the correction alone on a run's own levels (`../results/gcm/longwave_<run>_refit.json`) |
| [plasim_ice_optics.py](plasim_ice_optics.py) | The ice optics the runner's `--ice-optics` gives cloud colder than −38 °C, beside PlaSim's liquid-water optics (`../results/gcm/ice_optics.json`) |
| [high_cloud.py](high_cloud.py) | How far the top three layers' water and crystal radius could move a run's climate: first-order changes at the top of the air and the ground (`../results/gcm/high_cloud_<run>.json`) |
| [rain_evaporation.py](rain_evaporation.py) | How much of its falling rain a run evaporates, measured from the rain its schemes form against the rain at the ground, against CM1's lunar rain, and the constant that matches it (`../results/gcm/rain_evaporation.json`) |
| [cold_cloud.py](cold_cloud.py) | The cold cloud at lunar gravity: factors on the top layers' water and crystal radius from slabs of rising air on both worlds ([atmosphere/snow/cirrus.py](../../atmosphere/snow/cirrus.py)), for the runner's `--cold-cloud` (`../results/gcm/cold_cloud_<run>.json`) |
| [compare.py](compare.py) | The numbers behind the run tables, and where a run still warming or cooling would settle |
| [climatology.py](climatology.py) | A run's surface climatology over chosen years, for display and other lanes: cloud in a low deck (the lowest 34 km) and a high deck, precipitation, evaporation, soil moisture, temperature, water vapour, sea ice and snow, with cloud composited against the Sun by latitude and local time, and land air temperature the same way for all land and for the near and far sides. It writes `products/climatology_<run>.npz` (`terluna.climate.gcm-climatology/1`) |
| [sun_check.py](sun_check.py) | Where a run's Sun shines: its sunlight along the equator against PlaSim's solar-position arithmetic, and by latitude against the orbit the run asked for and the one PlaSim's log reports. It writes `../results/gcm/sun_check_<run>.json` in seconds |
| [radiation_check.py](radiation_check.py) | How the GCM divides the clear-sky sunlight arriving over the lunar day (reflected, absorbed in the air, absorbed by the ground) against the line-by-line model, on 11 columns of a run from dry land to seas at 60°. It writes `../results/gcm/radiation_check_<run>.json` in about half an hour on 8 cores |
| [site_winds.py](site_winds.py) | A run's winds at one site on every model layer, for other lanes: each layer's height above the model ground (hypsometric, with the design air), that ground's height above the model's sea level, the speed's percentiles, mean cube and yearly highest, the mean vector, and the speed by local hour angle. It writes `../results/gcm/site_winds_<run>_<site>.json` (`terluna.climate.gcm-site-winds/1`), which is committed; the summit of the far-side highlands is the first site |
| [global_winds.py](global_winds.py) | A run's winds and air over the whole Moon at fixed heights above sea level, for other lanes: the air's density (mean and 1st–99th percentiles), pressure and temperature, the speed's percentiles, highest, yearly highest and 50-year return level, and the speed by band of latitude, all weighted by area over the ground below each height. It writes `../results/gcm/global_winds_<run>.json` (`terluna.climate.gcm-global-winds/1`), which is committed; the sky-ship study reads the flight band from it |
| [zonal_winds.py](zonal_winds.py) | A run's wind as a vector over the whole Moon, for other lanes and for the drifting towns and organisms of [zonal_winds.md](zonal_winds.md): the eastward and northward wind by latitude and height and, in the Sun's frame, by local solar hour; how often it holds the time of day or turns the Sun back; the solar day of a body drifting with it; the shear over 2.5–15 km; and passive routes at fixed heights through the steady Sun-frame flow and through the 3-day means in sequence, with the points where the steady flow stops. It writes `../results/gcm/zonal_winds_<run>.json` (`terluna.climate.gcm-zonal-winds/1`) and a supporting `.npz` of Sun-frame percentiles, shear by hour and routes, both committed |
| [convection.py](convection.py) | A run's convective rain over the whole Moon, for other lanes: maps of the convection scheme's rain and all rain, the means of the convective rain's 1.5 power and square for nonlinear consumers, the share of outputs raining, and the land share, with the source files' hashes. It writes `../results/gcm/convection_<run>.json` in seconds; atmosphere/electricity/occurrence.py reads it |
| [twilight.py](twilight.py) | What the Moon's tall curved air does to sunlight and infrared against flat columns of the same air: the energy ledger by ground, GCM layer, shell and solar elevation, and the terms the GCM's and CM1's twilight patches take. It writes `../results/gcm/twilight.json` on the 294.9 K column (`terluna.climate.gcm-twilight/1`) and, with `--ts 301`, `../results/gcm/twilight_301k.json` on the design climate's (`/2`), both committed ([the curved air's light](#the-curved-airs-light)) |

`products/` is regenerated by `python -m climate.gcm.boundary` and, for the
climatology, `climate/gcm/.venv/bin/python -m climate.gcm.climatology A28_corrected_970:10-19`. It
is not committed. The site winds come from
`climate/gcm/.venv/bin/python -m climate.gcm.site_winds A28_dim5:15-24 summit 5.375 201.375`, and the winds over the whole Moon from
`climate/gcm/.venv/bin/python -m climate.gcm.global_winds A28_dim5:15-24`. Main corrected the GCM's Sun, tilt and
sunlight split on 2026-09-29, and its design run is now `A28_dim5_moon`, years 20–29. The same two products for that
run, `site_winds_A28_dim5_moon_summit.json` and `global_winds_A28_dim5_moon.json`, were made on 2026-10-08 for the
[infrastructure review](../../research/studies/infrastructure_review/README.md) with the same commands, the run given
as `A28_dim5_moon:20-29`.

## What only 3-D can answer

The one-dimensional work (`atmosphere/radiative_convective`,
`atmosphere/middle_atmosphere`) fixes the clear-sky radiation. The following
remain open, and each needs a circulation:

1. **Clouds.** The balance temperature moves about 7 K for every 10 W/m² of
   cloud effect on the Moon, against about 4 K on Earth. With Earth-like
   clouds (−20 W/m²) the 1.2-atm Moon balances near 284 K behind the titania
   shield and near 290 K behind a 200-nm edge that lets ozone form
   (`atmosphere/middle_atmosphere`, with its self-consistent stratosphere). The
   cloud effect of a slowly rotating world with a month-long day may differ
   strongly from Earth's. Slow-rotator GCM studies find thick dayside clouds,
   but at higher sunlight or with tidal locking (see "Cloud amount" below).
2. **The month-long day and night.** At 1.2 atm the lunar atmosphere holds 7.3
   times Earth's mass per square metre (six times per pascal of pressure). Its heat capacity equals an 18 m layer of
   water at 1.2 atm, which buffers the 15-day night. How much land and sea
   temperatures swing, and whether dayside convection and nightside inversions
   dominate, depends on the day–night circulation.
3. **The poles.** With a 1.54° tilt, the poles receive grazing sunlight all
   year. Without heat transport they would be far below freezing. Their
   temperature, sea ice and any permanent ice caps are set by transport.
4. **Water.** Where it rains, how the nearside and South Pole–Aitken seas
   exchange moisture with land, and whether moisture is cold-trapped at the
   poles.
5. **The middle atmosphere's day–night swing.** Radiation alone would swing the
   ozone-heated upper air by up to about 60 K between the lunar day and night
   (`atmosphere/middle_atmosphere/results`), or 1–10 K behind the titania
   shield. The tides and overturning this drives set the temperature at the
   base of the thermosphere, and through it the exobase. A GCM for this
   question needs a model top near 0.1 Pa, 550–820 km up on the Moon.

## Candidate models

| Model | Radiation | Suits | Effort |
|---|---|---|---|
| ExoPlaSim (PlaSim, Python-driven) | Band scheme designed for Earth-like air | Fast first look at slow rotation, the day–night cycle, land/sea and poles at T21 | `pip install exoplasim` plus a Fortran compiler and MPI |
| ROCKE-3D (NASA GISS ModelE) | SOCRATES with correlated-k, validated for N2–O2–H2O atmospheres at other pressures | The most complete physics for this problem (clouds, sea ice, dynamic ocean option); used for slowly rotating planets | ModelE source from NASA GISS, Fortran build with NetCDF and MPI; days of setup |
| LMD Generic PCM | Correlated-k tables built for the atmosphere | Flexible composition and radiation; a strong exoplanet record | Fortran build, k-table generation |
| Isca (GFDL-based) | Grey, RRTM or SOCRATES | Controlled experiments, easy configuration | Python and Fortran build |

A practical sequence:

1. ExoPlaSim first, to map the regime (cloud cover, day–night swing, polar
   temperatures, circulation pattern), at 1.0 and 1.2 atm with 25% and 35%
   water.
2. Then ROCKE-3D or the LMD model for the cases that matter, with radiation
   checked against this repository's line-by-line fluxes on the same profiles:
   `atmosphere/middle_atmosphere/results/radiation_benchmarks.json` holds four
   columns (the Earth control, the Moon behind the titania stack at 1.2 and
   1.0 atm, and behind a 200-nm edge). For each it gives thermal fluxes line by
   line and by correlated-k, and solar fluxes at zenith cosines 1.0 and 0.5
   and as a global mean.

GPUs do not help these codes; all run on CPU cores with MPI. On this machine a
T21 ExoPlaSim year takes about 3.3 minutes on eight cores. ROCKE-3D fits it too,
with a slab ocean (its Q-flux ocean with no ocean heat transport, as in the
PlaSim runs) at its standard 4°×5° grid, which GISS used for its own lunar runs
(Aleinov et al. 2019, GRL 46, 5107). The Moon's small radius makes those cells
3.7 times narrower than on Earth and shortens the model's time step, so a lunar
year should cost two to four times the 11 CPU-hours GISS published for an Earth
year (Way et al. 2017): about 2–5 hours on eight cores, by estimate. A 20–30
year run started from PlaSim's settled state then takes a few weeks of sessions
with rests, and a one-year benchmark would pin the cost down. Rented compute is
the faster route.

## ExoPlaSim runs (2026-09-25)

**The sunlight of runs before 2026-09-29 has two faults (found and fixed that
day).** They have Earth's 23.7° tilt in place of the Moon's 1.54°, and a Sun
that sweeps unevenly across longitude (both under "How the Moon is set up"
below). Their global total of sunlight is right. Where it falls is not, so
their results that vary with latitude (above all near the poles), with the
seasons or with longitude are not the Moon's, and their global means can differ
from the corrected runs'. The runner now fixes both, and the radiation is
recalibrated with the Sun corrected.

ExoPlaSim 3.4.2 is installed in `climate/gcm/.venv` (ignored by Git), built with
gfortran and Open MPI. Its postprocessor needs `meson` and `ninja` in the same
environment to compile its FFT module.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --years 40   # run or resume
climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --status     # progress
climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A --stop       # stop cleanly
climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A28 --years 20 --folder A28_gases_glow --co2-ppm 755 \
    --top-infrared 10.1718 0.5165 --set sunlight_scale=0.909867 \
    --start-from climate/gcm/runs/A28_dim5_moon/model/MOST_REST.00029   # a branch with CO2 and infrared from space
```

Runs live in `climate/gcm/runs/<experiment>/`. On the author's machine that folder
links to the external drive, `/media/projectspace/terluna-research/gcm-runs`.
Each run keeps:
- every year's restart, about 7 MB;
- every year's output, about 75 MB: 120 means per model year, ten per lunar day,
  of the fields that answer the questions above.

`--prune` keeps less, for small disks. A simulated year takes about 2½ minutes on
8 MPI ranks at T21 (32×64 cells, 170 km on the Moon), about 3½ with the clear-sky
diagnostic the cloud runs switched on, and all ranks together use about 0.3 GB.

**Restartable.** The runner advances one year at a time. After each year it writes
`progress.json` atomically, with a summary of the year:
- global mean temperature, top-of-atmosphere balance and planetary albedo;
- with the clear-sky diagnostic, the clouds' effect in sunlight and infrared;
- cloud cover, precipitation and sea ice;
- the range of equatorial land's 3-day means over the year (day and night with
  the weather) and the polar temperature.

A stopped run resumes from the last completed year's restart, so at most one year
(a few minutes) is lost. It first removes the interrupted year's files and any
model processes still running in the folder. Two runners cannot share a folder.
A run refuses to resume if anything that determines the model has changed: its
settings, the number of ranks, ExoPlaSim's version or the input files the runner
writes. The seed is fixed, so runs repeat exactly. The test: a run killed with
`kill -9` a minute into its second year left ten model processes running; the
resume stopped them, and its second year's restart and output were
byte-identical to an uninterrupted run's.

**How the Moon is set up.**
- **Planet:** radius, gravity, the true 27.32-day rotation and the Moon's 1.54°
  tilt from `products/`, the orbit held to them (`fixed_orbit`; runs before
  2026-09-29 have Earth's, see "Its tilt was Earth's" below).
- **Air:** the 1-D models' dry composition, with 1.2 atm at the water level. Land
  stands above it, so the global mean surface pressure is 1.16 atm and ranges
  from 1.03 to 1.22 atm.
- **Land and sea:** the T21 land-sea file, made binary by filling cells with sea in
  order of least land until the sea covers the water share exactly.
- **Ocean:** a 50 m slab with sea ice.
- **Land albedo:** a uniform 0.2, a labelled placeholder.
- **Clouds:** PlaSim's diagnostic scheme, with the cloud water its radiation sees
  set by the runner's `cloud_water` setting (see "Clouds" below).
- **Sunlight:** the TSIS-1 solar spectrum times the shield's transmission. For the
  titania stack that gives 1235 W/m², and it sets the split between ExoPlaSim's two
  solar bands.
- **Calendar:** ExoPlaSim counts the Moon's rotations in a 360-day year, 12,
  and PlaSim makes its model day 12/11 of the 27.32-day rotation: 29.8 days,
  1,430 half-hour steps. A model year is 12 of them, 17,160 steps, and the
  output holds ten means per model day (`whole_lunar_days`; before 2026-09-29
  each year ran ExoPlaSim's 17,280 steps, and the output numbers slipped 120
  steps against the Sun from one year to the next). The orbit follows the real
  year. The Sun keeps its own day of 29.53 (below), 3.2° of its path behind the
  model day each day, so day–night analyses sort by the Sun's position
  (`czen`).
- **Its Sun swept 393° a day and skipped back 33° (found and fixed 2026-09-29).** PlaSim's
  `solang` (in `radmod.f90`) turns the Sun at the planet's rotation rate, one
  turn in 27.32 days, and restarts its clock at the start of every model day of
  29.8. Each model day the Sun therefore sweeps 392.5° and then skips back
  32.5°, the same every day, so longitudes from 148° E to 180° E get noon twice
  a lunar day. Each row of latitude still gets its right mean. Along it,
  longitudes from 270° E through 0° to 57° E get 8.3% less and the far side
  around 164° E up to 17% more. [sun_check.py](sun_check.py) reproduces the
  pattern of `A28_dim5` from `solang`'s arithmetic alone (correlation 1.0000;
  [sun_check_A28_dim5.json](../results/gcm/sun_check_A28_dim5.json)). The skip
  is 1/(N − 1) of a turn for a planet whose rotations ExoPlaSim counts as N in a
  360-day year: about a degree for N near 360, and 33° for the Moon's 12. It is
  in ExoPlaSim 3.4.2, the latest release, and in its current source. The
  runner's patch to `solang` (`sun_clock`) now turns the Sun evenly at the
  synodic rate that the rotation and the orbit give, one turn per 29.53 days. In
  the two-year test `A28_dim5_sun_test` the equator's sunlight varies along the
  row only by what 12.1 lunar days in a year leave, 0.991–1.017 of its mean, as
  `solang`'s new arithmetic predicts
  ([sun_check_A28_dim5_sun_test.json](../results/gcm/sun_check_A28_dim5_sun_test.json)).
- **Its tilt was Earth's (found and fixed 2026-09-29).** ExoPlaSim uses the obliquity it
  is given only when the orbit is fixed (`fixedorbit=True`, PlaSim's
  `NFIXORB=1`). The runner kept the default, with which PlaSim replaces the
  namelist's orbit with Earth's for its start year, 1 AD, from Berger's (1978)
  series: obliquity 23.695° and eccentricity 0.0175. Every run's log reports
  that orbit, the calibration runs' included, and the runs have Earth's
  seasons. In `A28_dim5` the rows' mean sunlight follows that orbit within
  1.8 W/m² rms and misses the Moon's by 44: 4% short on the equator, 6% over at
  58°, 29% over at 69°, and 149 W/m² at 86° against the Moon's 29, where the
  Sun stays down for about five of the year's twelve lunar days. The runner now
  holds the orbit (`fixed_orbit`): the test's log reports 1.5424° and 0.016715,
  and its rows follow the Moon's orbit within 0.1 W/m² rms.
- **Together** the two faults left 328 W/m² of the 373 due at 0° E over a
  lunar day: the sweep took about 30 and the tilt about 15. Every result below
  from runs before 2026-09-29 that varies with latitude, season or longitude
  carries them, as did the products other lanes built from those runs. Those
  products were regenerated from the corrected design case ("The corrected
  design case" below), except the forcing the CM1 cases were set up with and the
  atlas's clouds, which came from run A until the globe took the design case's on
  16 October.

**Radiation, checked against line-by-line and calibrated.** The line-by-line model
(`atmosphere/radiative_convective`) was run on the GCM's own cloud-free columns:
the same temperature and humidity profiles, the same surface albedo, and
sunlight through the titania stack. The unmodified model was fine when cool and
dry, with outgoing longwave within +3 ± 4 W/m². When warm and moist (310 K,
450 kg/m² of water vapour) it went wrong in two ways:
- it reflected too little sunlight: clear-sky albedo 0.203 against 0.248;
- it emitted 12 ± 2 W/m² too little longwave.

Together that is about 26 W/m² of spurious heating, which drove the uncalibrated
design case past 310 K and still warming (`runs/A_uncalibrated`). Two
adjustments bring it within about ±2 W/m² of line-by-line at both 280 and 310 K:
- **Rayleigh multiplier (1.8).** ExoPlaSim derives a Rayleigh coefficient from
  the spectrum, normalised to a 5772 K blackbody's ultraviolet. For the
  titania-filtered Sun that gives 0.21, too low. A one-line patch to PlaSim's
  `radmod.f90` (the `RAYSCALE` namelist multiplier) multiplies it by 1.8. The
  runner applies the patch itself, and it reproduces the unmodified model
  exactly when set to 1.
- **Water-vapour continuum (0.004).** PlaSim's water-vapour absorptivity is fitted
  to Earth's range of water paths. At lunar paths its continuum term over-absorbs,
  so its coefficient `TH2OC` is lowered from 0.024 to 0.004.

The trials, branched for a year from the uncalibrated run's cold (year 1) and
warm (year 9) states, are in `runs/calibration/`:

| Rayleigh multiplier, continuum | Albedo, cold: GCM / line-by-line | Albedo, warm | Longwave GCM − line-by-line, cold / warm |
|---|---|---|---|
| 1.0, 0.024 (unmodified) | – | 0.203 / 0.248 | +3.4 / −11.8 W/m² |
| 1.5, 0.012 | 0.250 / 0.275 | 0.235 / 0.247 | −0.5 / −8.7 W/m² |
| 1.9, 0.004 | 0.276 / 0.276 | 0.261 / 0.247 | +2.1 / −2.1 W/m² |
| **1.8, 0.004 (used until 2026-09-29)** | about −0.007 | about +0.007 | about ±2 W/m² |

What the calibration does not reach is the clouds. PlaSim diagnoses them from
humidity, with Earth-tuned optical properties, and the cloud effect moves the
balance about 7 K per 10 W/m². On the Moon one of those Earth fits, the clouds'
water, fails outright (see "Clouds" below). ExoPlaSim's global mean temperature
therefore reads as "given PlaSim's clouds", and a GCM with cloud physics
(ROCKE-3D) remains the way to settle it.

**How it splits the sunlight (checked 2026-09-29).** The calibration compared
clear-sky albedo and outgoing longwave on the global-mean column.
[radiation_check.py](radiation_check.py) compares how the GCM divides the
sunlight arriving over the lunar day in clear sky with line-by-line, on 11
columns of `A28_dim5` (years 15–24): land from its driest to its most humid, a
tropical highland, and seas from the equator to 60°
([radiation_check_A28_dim5.json](../results/gcm/radiation_check_A28_dim5.json)).
In every column the GCM's air absorbs 7.7–8.2 points of the arriving sunlight
less than line-by-line's, and its ground takes 7–12 points more. Its reflection
matches over land (within 0.5 points, 1.8 low in the driest column) and runs
1.6–4.1 points low over the seas, lowest at 60°. At 0° E that is 26 W/m² too
little absorbed in the air and 25 W/m² too much at the ground. The line-by-line
side keeps the Sun over the equator all year, so beyond the tropics, where the
runs' Sun follows Earth's seasons, the comparison also carries a difference in
Sun heights; the tropical columns show the same 8 points. Neither calibrated
setting (Rayleigh scattering, the longwave continuum) moves sunlight into the
air. PlaSim's near-infrared water-vapour absorption, fitted to Earth's water
paths, is the likely place to recalibrate.

**Recalibrated (2026-09-29).** PlaSim's near-infrared band absorbs through
Lacis and Hansen's (1974) formula, fitted to Earth's water paths. It levels off
toward half the sunlight, so over the Moon's paths, several times Earth's, it
absorbs too little. [plasim_shortwave.py](plasim_shortwave.py) repeats PlaSim's
clear-sky shortwave for one column, and on the corrected Sun's two-year test it
reproduces the GCM's own clear-sky split within 0.4 points on all 11 columns
([shortwave_validate_A28_dim5_sun_test.json](../results/gcm/shortwave_validate_A28_dim5_sun_test.json)).
Fitted against line-by-line on the radiation check's columns
([shortwave_fit_A28_dim5.json](../results/gcm/shortwave_fit_A28_dim5.json)), a
multiplier of 1.356 on the absorptivity (the runner's `sw_water_absorption`, a
patch to `swr` that keeps the band's transmission above 0.001), with the
Rayleigh multiplier raised from 1.8 to 2.215 to keep the reflection, leaves the
air 0.3 points and the reflection 0.7 points rms from line-by-line; before, the
air was 8.1 points short. Of three ways tried to raise the absorption (a longer
path, a larger absorptivity, a term growing with the path), the larger
absorptivity carried over best to two columns the fit did not see: a warmer,
moister column from the full shield's run (303 K, 338 kg/m² of water vapour),
where the air stays 0.7 points short against 9.7 before, and a cooler one from
the 8%-dimmer run, 0.2 points short. In the corrected design case's own output
(years 20–29), on 11 new columns, the GCM's clear air takes line-by-line's share
within 0.6 points rms, its ground within 0.6 and its reflection within 0.8
([radiation_check_A28_dim5_moon.json](../results/gcm/radiation_check_A28_dim5_moon.json));
at 0° E its air takes 31.4%, as line-by-line's does.

### Run A, the design case

Run A is 1.2 atm with 25% water, behind the titania stack, with calibrated
radiation, for 40 years. It started from ExoPlaSim's cold state: 250 K, with the
seas 87% frozen. The ice was gone by year 3, and the climate had settled by
year 27, with the top-of-atmosphere balance within ±0.3 W/m². Years 30–39:

| Global | |
|---|---|
| Surface temperature | 295.5 K (22 °C); sea surface 299.2 K |
| Top-of-atmosphere balance | +0.2 W/m² (220.4 absorbed, 220.2 emitted) |
| Planetary albedo | 0.287 |
| Cloud cover | 28% |
| Precipitation | 3.3 mm/day |
| Water vapour | 198 kg/m² |
| Sea ice | none |

- **A warm, even world.** Air near the ground is 295 ± 1 K from the equator to
  the poles, over land and sea. The run's poles had Earth's tilt, though, and
  with it Earth's polar summers: at 86° five times the sunlight of the Moon's
  grazing Sun (see "Its tilt is Earth's"). A slowly turning Moon's circulation
  spans the globe, and its heavy air carries heat efficiently, but how warm the
  Moon's own poles stay waits on the corrected runs. The single column put 80°
  land at 276–283 K. It could exchange heat only through its free atmosphere,
  not through this deep mixing.
- **The month-long day is mild on land.** Over equatorial land the air near the
  ground averages 297.3 K in sunlit 3-day means and 293.9 K in dark ones. Over a
  year, 3-day means at a place range about 9 K at the equator and 11 K at 60°
  (with the runs' Earth seasons), day and night together with the weather. The
  extremes anywhere are 288–307 K
  (15–34 °C). The sea runs a few kelvin warmer than the land.
- **Rain follows the equator.** It is 5.8 mm/day there and 0.1–0.2 mm/day at
  the poles; cloud cover falls from 36% to 10%.
- **Clouds are the weak point.** PlaSim's clear-sky diagnostic, switched on for a
  rerun of year 39 (`runs/A_patch_check`, identical to the original), gives the
  clouds' effect on the top-of-atmosphere balance: −11.0 W/m² in sunlight and
  +0.3 W/m² in the infrared, −10.7 W/m² net. (An earlier estimate of −4 W/m²,
  read off the 1-D balance temperatures, was wrong.) Earth's clouds give about
  −47, +28 and −20 W/m². The infrared effect is nearly nil because the model's
  high clouds hold no water; see "Clouds" below.

### Run B, more water

Run B is the same as A but with 35% of the surface under water. Years 30–39
against run A (the run tables come from
`climate/gcm/.venv/bin/python -m climate.gcm.compare A:30-39 B:30-39`):

| | A, 25% water | B, 35% water |
|---|---|---|
| Surface temperature | 295.5 K | 297.2 K (still rising about 0.05 K a year) |
| Sea surface | 299.2 K | 300.4 K |
| Air over land | 295.4 K | 296.6 K |
| Planetary albedo | 0.287 | 0.282 |
| Precipitation | 3.30 mm/day | 3.64 mm/day |
| Water vapour | 198 kg/m² | 220 kg/m² |
| Equatorial land air, sunlit / dark 3-day means | 297.3 / 293.9 K | 297.9 / 295.2 K |
| Equatorial land, range of 3-day means at a place in a year | 8.7 K | 7.3 K |
| Extremes of 3-day means | 288–307 K | 290–308 K |

More sea makes the Moon about 2 K warmer and wetter, and it damps the land's
day–night swing. The sea is darker than the placeholder land albedo of 0.2, and
it keeps more vapour in the air. Neither run has sea ice.

### Clouds: an Earth-fitted cloud-water formula, and a bracket

PlaSim carries no cloud water. Its condensation and convection schemes rain out
any excess vapour at once, so the rain (3.3 mm/day in run A) does not depend on
it. For the clouds' radiation, how much sunlight they reflect and how much
infrared they trap, PlaSim diagnoses their water separately. It uses a formula
from NCAR's CCM3 (Kiehl et al. 1996) fitted to Earth: 0.21 g/m³ at the ground,
falling off with geometric height on a scale of 700 ln(1 + PW) m, where PW is the
precipitable water in kg/m². At a sixth of Earth's gravity a pressure level
stands six times higher, so on the Moon the formula empties the upper clouds.
Run A's final year, layer by layer (`python -m climate.gcm.cloud_water`; in-cloud
water paths weighted by cloud cover, from 3-day means):

| Level (σ) | Height on the Moon | Cloud cover | Water path, PlaSim as it is | Earth's water path | Full column |
|---|---|---|---|---|---|
| 0.12 | 95 km | 10% | 0 g/m² | 1.5 g/m² | 9 g/m² |
| 0.21 | 73 km | 9% | 0 | 5 | 32 |
| 0.32 | 56 km | 4% | 0 | 14 | 86 |
| 0.44 and 0.57 | 41 and 29 km | 1.4% each | 0.1 and 1.3 | 32 and 61 | 190 and 370 |
| 0.70 | 18 km | 2% | 16 | 98 | 590 |
| 0.82 | 10 km | 3% | 105 | 129 | 780 |
| 0.92 | 4 km | 3% | 333 | 126 | 760 |
| 0.98 | 0.9 km | 2% | 291 | 58 | 348 |

The upper layers hold most of the model's cloud cover, together about a quarter
of the sky, and in PlaSim as it is no water: they neither reflect sunlight nor
trap infrared. The lowest layers, in a column with six times Earth's mass per
pascal of pressure, come out wetter than Earth's clouds.

A patch to PlaSim's `rainmod.f90`, applied by the runner like the radiation
patch, adds two settings to the formula: the gravity at which it counts heights
and precipitable water, and a multiplier on the water. The runner's
`cloud_water` setting uses them three ways:
- `plasim`: the formula as PlaSim has it (runs A and B);
- `earth_path`: heights and precipitable water counted as on Earth, which puts
  the formula in pressure terms, and each cloud holding the water path it would
  hold at the same pressures on Earth;
- `full_column`: the same mixing ratio over the Moon's column, which holds six
  times Earth's mass per pascal: six times Earth's water path. This is the upper bound: the deeper clouds rain none of it out.

With `plasim` the patch changes nothing: year 39 of run A, rerun with the
patched build and the clear-sky diagnostic, came out byte-identical
(`runs/A_patch_check`). New runs default to `earth_path` and to the
convective fix below. Runs A and B predate these settings, so the runner now
refuses to extend them in place; a continuation branches from their last
restart with `--start-from`.

Both corrections were branched from run A's year 39 and run for 15 years with
the clear-sky diagnostic on (`runs/A_earth_path_clouds`,
`runs/A_full_column_clouds`). Each settled within about a decade. Years 10–14
against run A (A's cloud effect is from its rerun year 39; "settles near" is
`compare --settle` over years 1–14):

| | A: PlaSim's cloud water | Earth's water path | Full column |
|---|---|---|---|
| Surface temperature | 295.5 K | 297.8 K, settling near 298.0 | 292.3 K, settling near 292.0 |
| Sea surface | 299.2 K | 301.0 K | 296.6 K |
| Air, equator / poleward of 70° | 295.8 / 295.2 K | 298.0 / 297.6 K | 292.9 / 292.2 K |
| Equatorial land air, sunlit / dark 3-day means | 297.3 / 293.9 K | 299.2 / 296.2 K | 294.3 / 291.0 K |
| Extremes of 3-day means | 288–307 K | 290–309 K | 284–303 K |
| Planetary albedo | 0.287 | 0.302 | 0.347 |
| Cloud effect: sunlight, infrared, net | −11.0, +0.3, −10.7 W/m² | −16.2, +8.4, −7.8 W/m² | −29.2, +14.1, −15.2 W/m² |
| Cloud cover | 28% | 22% | 23% |
| Precipitation | 3.30 mm/day | 3.33 mm/day | 2.62 mm/day |
| Water vapour | 198 kg/m² | 234 kg/m² | 168 kg/m² |
| Sea ice | none | none | none |

- **Watery high clouds warm the Moon.** Given Earth's water path, the high clouds
  reflect some sunlight but trap more infrared, day and night. The net cloud
  effect weakens from −10.7 to −7.8 W/m², and the Moon settles near 298 K,
  2.5 K warmer than run A. Only with six times Earth's water does reflection win
  clearly: −15.2 W/m² net, and the Moon settles near 292 K.
- **Given PlaSim's cloud amount, the design case lies between 292 and 298 K.**
  Earth's water path is the best single estimate. Cloud water is limited by how
  fast clouds rain out, and rain here falls at Earth's rate. The full column is
  the upper bound on cloud water. Run A's 295.5 K is not a bound: PlaSim's
  formula takes water out of the high clouds and adds it to the low ones, and
  the two errors partly cancel.
- **The pattern holds.** In every case the air is within 0.7 K from the equator
  to the poles, sunlit and dark 3-day means over equatorial land differ by about
  3 K, and there is no sea ice. Across the bracket the global temperature moves
  about 0.8 K per W/m² of cloud effect.
- **Cloud amount is not bracketed.** Cover stays at 22–28%, against Earth's
  roughly two-thirds. PlaSim diagnoses it from humidity and convective rain with
  Earth-fitted formulas, and it builds no thick deck under the Sun. See "Cloud
  amount" below for what published models imply.

### Convective cloud: a second Earth-unit slip

PlaSim's convective cloud cover is Slingo's fit to Earth, 0.245 + 0.125 ln P
with P the convective rain in mm per day, capped at 0.8. PlaSim divides the rain
by its own solar day, which on the Moon is 29.8 Earth days (2,575,191 s), so the
formula reads the rain about 30 times too heavy and gives 0.42 more cover, up to
the cap. A second patch to `rainmod.f90` adds `CONVDAY`, the day length in which
the formula counts rain. The runner's `convective_day_s` sets it: 86400 by
default (an Earth day, as fitted), 0 for PlaSim's behaviour. With 0 the patched
build reproduces year 14 of `A_earth_path_clouds` byte for byte, and with 86400
it reproduces the corrected run's first year (`runs/regression_convday_off` and
`_on`).

`runs/A_corrected_clouds` branches from `A_earth_path_clouds` at year 14 with
both corrections and runs 15 years. Years 10–14:

| | Earth's water path | + rain counted per Earth day |
|---|---|---|
| Surface temperature | 297.8 K, settling near 298.0 | 299.6 K, settling near 299.6 |
| Sea surface | 301.0 K | 302.4 K |
| Air, equator / poleward of 70° | 298.0 / 297.6 K | 299.7 / 299.3 K |
| Equatorial land air, sunlit / dark 3-day means | 299.2 / 296.2 K | 300.8 / 298.0 K |
| Planetary albedo | 0.302 | 0.293 |
| Cloud effect: sunlight, infrared, net | −16.2, +8.4, −7.8 W/m² | −13.4, +7.6, −5.8 W/m² |
| Cloud cover: global (under a high Sun / at night) | 22% (37% / 15%) | 20% (33% / 14%) |
| Precipitation | 3.33 mm/day | 3.50 mm/day |
| Water vapour | 234 kg/m² | 263 kg/m² |
| Sea ice | none | none |

- **With both slips corrected and PlaSim's cloud amount, the design case is
  about 299.6 K (26.5 °C).** Fewer convective clouds reflect less sunlight, and
  the Moon settles 1.6 K warmer than with the cloud-water fix alone.
- **The Moon is winterless.** Over years 10–14 the ground never falls below
  286.5 K (13 °C) in any half-hour, on a 4 km plateau at 75°S late in the
  night, and its hottest half-hour is 315.9 K (43 °C), over equatorial land near
  noon. The air near the ground stays within 292.0–310.6 K (19–38 °C) in 3-day
  means. No snow or sea ice forms. Even the full-column clouds' 292 K case
  bottoms out at 277 K (4 °C) on the ground. The freezing level stands 32–34 km
  up (23–25 km in the 292 K case), and the model's highest land is 8.2 km: at a
  sixth of Earth's gravity the air cools only about 1 K per km, so mountains
  barely reach cold air. The coldest ground sits 13–15 K below the global mean,
  so frost on the highest ground would need a Moon about as cool as Earth on
  average.
- **An output artefact to screen.** In 0.8% of 3-day windows the output's lowest
  model level reads up to about 20 K colder than the air near the ground, which
  no real atmosphere sustains. It most likely comes from converting the model's
  spectral temperatures to the grid over steep terrain. The surface fields are
  unaffected; analyses of the lowest level should screen it.
- **What the output's temperatures are.** `tas` is the lowest model layer (the
  bottom 1.8 km, centred 0.9 km up) carried to the ground along a dry adiabat
  and averaged over each 3-day window. `maxt` and `mint` are the extremes of the
  ground surface temperature within the window: PlaSim's `outmod.f90` tracks
  `dt(:,NLEP)`, the surface, for them. Heat, cold and day–night figures in this
  file use `tas` for the air and `maxt`/`mint` for the ground; on 2026-09-26
  the tables below were corrected where they had listed ground extremes as air.

### A28: the 28% seas with their lakes

The scenario's water share is 28% (`shared/scenarios/water.json`), and the
geography domain's drainage estimate adds rain-fed lakes above sea level. Run
`A28` takes both: the boundary builder combines the atlas's seas with the
drainage estimate's lakes (`products/moon_28pct_water_lakes_*.nc`, with each
cell's water surface height), and PlaSim's binary mask then makes 759 of the
2,048 T21 cells water, 40% of the surface. 201 of them are lakes, standing at
their own levels up to 8.1 km above sea level. A restart carries PlaSim's
land-sea mask and orography, so a new coastline needs a cold start: `A28` ran
30 years from ExoPlaSim's cold state with both cloud corrections, paused
overnight after year 7 and resumed. Years 25–29, against the corrected 25% run:

| | 25% water (`A_corrected_clouds`) | 28% seas and lakes, 40% water (`A28`) |
|---|---|---|
| Surface temperature | 299.6 K | 302.8 K, settling near 303.2 |
| Sea surface | 302.4 K | 304.3 K |
| Air over land; within 12° of the equator / poleward of 70° | 299.3; 299.7 / 299.3 K | 302.4; 302.4 / 302.0 K |
| Land air, annual means: 5th / 50th / 95th percentile | 296.6 / 299.3 / 301.5 K | 300.3 / 302.3 / 304.4 K |
| Coldest land cell, annual mean | 295.5 K, 7.2 km up | 297.1 K, 8.4 km up |
| Air near the ground, 3-day means, coldest–warmest | 292.0–310.6 K | 294.5–312.8 K |
| Ground surface, half-hourly extremes | 286.5–315.9 K | 289.4–317.8 K |
| Planetary albedo | 0.293 | 0.278 |
| Cloud effect: sunlight, infrared, net | −13.4, +7.6, −5.8 W/m² | −13.0, +6.7, −6.3 W/m² |
| Cloud cover | 20% | 17% |
| Precipitation; land under 0.5 / over 5 mm/day | 3.50 mm/day; 3% / 28% | 3.70 mm/day; 2% / 20% |
| Water vapour | 263 kg/m² | 307 kg/m² |

- **More water makes the Moon warmer and no more varied.** The water is darker
  than the placeholder land and fills the air with vapour: the Moon settles
  3.6 K warmer, near 303 K (30 °C), and annual means over land still span only
  about 4 K. The lakes, a mean 4.5 km above sea level, run 0.8 K cooler than the
  seas. Nowhere is cool: the coldest land, 8.4 km up, averages 24 °C.
- **Lakes and climate agree.** The lakes came from run A's runoff; routed again
  with `A28`'s own climatology (`python -m climate.gcm.climatology A28:20-29`,
  then `python -m geography.drainage --climatology
  climate/gcm/products/climatology_A28.npz`), they cover 12.9% instead of 12.0%
  and would flip 84 T21 cells, raising the water from 40.1% to 40.9%, worth
  about 0.1 K. The geography domain now carries that routing.

### A dimmer shield: cooling the 28% Moon

At 303 K the 28% Moon is hot and humid. Over land, 2% of 3-day means of the air
reach 35 °C, and the humid heat peaks near 32 °C wet-bulb, close to where healthy
people fail. The runner's `sunlight_scale` dims the shield evenly at every wavelength;
branches start from `A28`'s year 29, whose coastline they keep. After 8 years
each, with the settling fit (`compare --settle`) and heat measures from years
6–7 (their extremes are still falling):

| | Full shield | 3% dimmer | 5% dimmer | 8% dimmer |
|---|---|---|---|---|
| Sunlight at the top of the air | 1,235 W/m² | 1,198 W/m² | 1,173 W/m² | 1,136 W/m² |
| Settles near | 303 K | 297 K | 294 K | 288 K |
| Air over land, 3-day means, 90th / 99th percentile | 33 / 36 °C | 29 / 32 °C | 27 / 30 °C | 23 / 26 °C |
| Wet-bulb of those means, 99th percentile | 29 °C | 26 °C | 24 °C | 21 °C |
| Coldest 3-day mean of the air; coldest half-hour of the ground | 22; 18 °C | 18; 13 °C | 15; 8 °C | 10; 3 °C |

The response is linear, about 1.9 K per 1% of sunlight, as the energy balance
alone predicts. Cloud cover rises as the Moon cools (17% at 303 K, 25% at
293 K), but the clouds' net effect holds at −6 to −7 W/m²: the extra reflection
and the extra infrared they trap cancel. Whether the cover itself would rise
like this is beyond PlaSim's diagnostic clouds; a cloud-resolving model run at
two temperatures could test it.

`A28_dim5` (5% dimmer) was run on to 25 years and has settled. Years 20–24,
against the full shield:

| | Full shield (`A28`, years 25–29) | 5% dimmer (`A28_dim5`, years 20–24) |
|---|---|---|
| Surface temperature | 302.8 K | 294.2 K, settling near 293.8 |
| Sea surface | 304.3 K | 297.2 K |
| Land air, annual means: 5th / 50th / 95th percentile | 300.5 / 302.5 / 304.5 K | 292.0 / 293.5 / 295.4 K |
| Equatorial land, air at the warmest / coldest time of the lunar day | 32 / 28 °C | 24 / 19 °C |
| Equatorial land, ground at its typical extremes | 35 / 26 °C | 27 / 15 °C |
| Land above 5 km, air at the warmest / coldest time | 27 / 26 °C | 20 / 18 °C |
| Air over land, 3-day means: 90th / 99th percentile / warmest | 33 / 36 / 40 °C | 24 / 27 / 31 °C |
| Wet-bulb of those means: median / 99th percentile / highest | 24 / 29 / 32 °C | 17 / 22 / 25 °C |
| Coldest 3-day mean of the air; coldest half-hour of the ground | 21; 16 °C | 12; 3 °C, no frost |
| Cloud cover; cloud effect net | 17%; −6.3 W/m² | 25%; −6.2 W/m² |
| Precipitation; land under 0.5 mm/day | 3.70 mm/day; 2% | 2.89 mm/day; 9% |
| Water vapour | 307 kg/m² | 173 kg/m² |

- **The 5% Moon is mild and safe to live in.** The air over equatorial land
  runs about 24 °C in the lunar afternoon and 19 °C late in the night, and the
  ground swings from 27 to 15 °C, a wider swing than at 303 K because the air
  holds half the vapour. Over land the warmest 3-day mean of the air is 31 °C,
  and humid heat stays well under the danger range.
- **A little more variety, mostly in water.** Dry land (under 0.5 mm/day)
  triples to 9%, while annual means over land still span only about 3 K; land
  above 5 km runs about 2.5 K cooler than equatorial land on average, 4–5 K by
  day and 1–2 K at night.
- **Its lakes stay close to the full shield's.** Routed with its own climate
  (`python -m climate.gcm.climatology A28_dim5:15-24`, then
  `python -m geography.drainage --climatology climate/gcm/products/climatology_A28_dim5.npz --out <folder>`),
  the lakes cover 11.6% instead of 12.9% and the rivers carry 25% less; at T21
  that would flip 105 cells and lower the water from 40.1% to 39.6%, worth
  about 0.1 K. The geography domain now carries this routing.
- **Humid heat is estimated** from each 3-day mean of the air with the humidity
  of the model's lowest layer (about 0.9 km up); the air at the ground could be
  about a degree more humid, and the 3-day means smooth the peaks of single
  afternoons.
- **Winter lives aloft.** The air cools about 0.6 K per km near the ground and
  0.8 K per km on average, so the 0 °C level stands 25.4 km above sea level
  (24.5–26.3 km over places and times, the same at the equator and poles, by
  day and by night), −10 °C at 36 km, −20 °C at 45 km and −40 °C at 61 km. The
  highest land, 8.4 km, never reaches it. The air thins on the way: at 0 °C it
  is 0.72 atm with 12.8 kPa of oxygen, what Earth's air gives at about 4,200 m;
  at −10 °C 0.58 atm and 10.3 kPa (about 6,100 m); at −20 °C 0.48 atm and
  8.5 kPa (about 7,700 m). Snow forms in the clouds above 25 km and melts on
  the way down. The corrected design case (below) keeps the same upper air: over
  its seas the freezing level stands where `A28_dim5`'s does, and the flight
  band's air differs by at most 0.2 °C.

The author chose the 5% dimmer shield (2026-09-26). Other lanes now read its
corrected run, below.

### The corrected design case

`A28_dim5_moon` is the chosen design, the 28% seas with lakes behind the 5%
dimmer shield, with the GCM's Sun, tilt and sunlight split corrected
(2026-09-29). It branched from `A28_dim5`'s year 24 and ran 30 years. It settled
within about five years: over years 20–29 its surface warmed 0.08 K per decade
and the top of the atmosphere gained 0.07 W/m². Along the equator its sunlight
is even to 0.1%, every row gets the Moon's sunlight within 0.1 W/m² rms
([sun_check_A28_dim5_moon.json](../results/gcm/sun_check_A28_dim5_moon.json)),
and its clear-sky split matches line-by-line on 11 columns (above,
"Recalibrated"). Against the uncorrected run
(`python -m climate.gcm.compare A28_dim5:15-24 A28_dim5_moon:20-29`):

| | `A28_dim5`, years 15–24 | `A28_dim5_moon`, years 20–29 |
|---|---|---|
| Surface temperature | 294.3 K | 294.9 K |
| Sea surface | 297.3 K | 297.0 K |
| Air over land | 293.7 K | 294.3 K |
| Air poleward of 70° | 293.2 K | 292.6 K |
| Coldest 3-day mean anywhere | 284.5 K | 287.4 K |
| Year's range of 3-day means, land at 55–65° | 10.8 K | 7.3 K |
| Planetary albedo | 0.293 | 0.280 |
| Cloud cover; cloud effect in sunlight / net | 25%; −15.7 / −6.3 W/m² | 23.5%; −8.5 / −1.6 W/m² |
| Precipitation | 2.91 mm/day | 2.24 mm/day |
| Land under 0.5 mm/day | 9% | 42% |
| Water vapour | 174 kg/m² | 176 kg/m² |

- **The Moon stays warm and even.** Its surface is 0.6 K warmer. The poles now
  get only the grazing sunlight of a 1.54° tilt, 29 W/m² at 86° against 149,
  and the air poleward of 70° stays at 292.6 K, 0.6 K cooler: the heavy air's
  circulation carries the heat, as the uncorrected runs found. Without seasons
  the year's range at a place shrinks, and no 3-day mean anywhere falls below
  14 °C.
- **Its rain is weaker and keeps to the tropics.** Tropical land keeps its rain,
  4.8 mm/day within 15° of the equator, but land between 30° and 45° gets
  1.3 mm/day instead of 2.8, and between 45° and 60° 0.2 instead of 1.9. Two
  changes do it. With only the Sun and tilt corrected (`A28_dim5_sun_test`, two
  years), the high latitudes already dry out and 36% of the land gets under
  0.5 mm/day: Earth's seasons had carried summer rain to 60°. The recalibrated
  sunlight then cuts the total rain by a further 18%, because the air now
  absorbs sunlight that used to evaporate water at the ground.
- **Clouds matter less.** Cover falls from 25% to 23.5%, and their net effect
  from −6.3 to −1.6 W/m², mostly in sunlight. About a third of that comes with
  the Sun and tilt corrected (−13 W/m² in sunlight in the two-year test), the
  rest with the air absorbing more of the sunlight around them.
- **The chosen shield still gives about 294 K**: 294.9 K, where the uncorrected
  run gave 294.3.
- **What reads it was regenerated from it:** the climatology
  (`python -m climate.gcm.climatology A28_dim5_moon:20-29`), geography's lakes
  and groundwater, the surface light at 294.9 K, the canopy models and the CM1
  comparisons. The lakes cover 10.3% of the Moon instead of 11.6%: those within
  15° of the equator grow a little, and those between 30° and 60° shrink to a
  fraction. Built into the GCM's T21 coastline they would change 111 of its
  2,048 cells and lower the water from 40.1% to 38.3%. The coastline needs a cold
  start; its surface was built for the 301 K design case on 19 October ([the lakes
  for a cold start](#the-lakes-for-a-cold-start)). The atlas's clouds then still
  came from run A.

### Its winds as vectors: the days of a drifting town

[zonal_winds.py](zonal_winds.py) keeps the direction of the design case's wind, which
[global_winds.py](global_winds.py) drops, and reads it in the frame of the slowly moving Sun, for the roaming sky
towns the author adopted on 9 October ([register](../../research/decisions.md#life-and-people)) and for aerophytes. The subsolar point moves west at v☉ = 2π(R + z) cos φ / 29.53 d, 4.30 m/s on the equator 10 km up, and a body
drifting east at u lives a solar day of 29.53 v☉ / |u + v☉| Earth days.
`climate/gcm/.venv/bin/python -m climate.gcm.zonal_winds A28_dim5_moon:20-29` reads years 20–29 in about 20 minutes
on one core and writes [zonal_winds_A28_dim5_moon.json](../results/gcm/zonal_winds_A28_dim5_moon.json) with its
`.npz`; [zonal_winds.md](zonal_winds.md) reads them for towns and aerophytes, by height and latitude.

- **Above about 8 km the air superrotates, and a drifting town lives a short day.** The mean eastward wind is
  0.85 m/s at 7.5 km, 3.5 m/s at 20 km, 10.7 m/s at 40 km and 21.9 m/s at 75 km. The median drifting day on the
  equator is 21.8 days at 10 km, 9.3 at 40 km and 4.8 at 75 km, and at 60° 17.3, 4.3 and 2.5 days.
- **Near the ground the flow follows the Sun.** On the equator at 2.5 km the air parts 15° of hour angle before dawn
  and meets 64° past noon (the cloud-resolving ring put these at 45° before dawn and 20° past noon), so a drifting
  town lingers through the evening: its day averages 34 days and it has the Sun up 44% of the time. There the hour at
  each place explains 85% of the eastward wind's variation about its zonal mean and the Sun's frame alone 50%; from
  10 km up they explain 19–39% and 7–14% of a variation of 1.0–1.6 m/s rms.
- **The wind holds the hour within 1 m/s in 18–31% of 3-day means in the lowest 5 km at 40–75°,** up to 66% in the
  late afternoon, and turns the Sun back there in 9–17%. From 20 km up no place equatorward of 76° holds the hour in
  30% of 3-day means. The steady Sun-frame flow stops only poleward of 65°, and in the 3-day means no route
  stays within a 30° span of hours longer than 45 days.
- **Drifting routes gather.** At 2.5 km the median distance between neighbouring routes falls from 300 km to 12 km in
  30 days, and within 180 days every route at that height shares one track at the model's resolution; at 40–60 km a
  poleward drift of 0.03–0.08 m/s brings routes into the polar caps within a few years.
- **The shear over 15 km** is about 4 m/s in the lowest 20 km and 6 m/s in and below the 35–45 km flight band,
  where it is the eastward wind's growth with height and changes by about 10% through the solar day.

The product's statistics run over heights 0–90 km above the model's sea level every 2.5 km, the model's latitude
rows, and hour angles in degrees east of the subsolar point in 36 bins, as [climatology.py](climatology.py) counts
them; its reading rule names every array. They describe the large-scale flow of the 3-day means at T21, and the
routes are carried at fixed heights by that flow alone.

### The gases, the glow and the dimming

On 13 October the author asked for the GCM to run again with the held gas pair, the ring fleet's glow and the dimming
that answers them. The runner gained two options for it. Each is part of the run key and recorded only when given, so
earlier runs keep their keys:

- `--co2-ppm` sets the air's CO2 in place of 400 ppm, N2 the balance. PlaSim's longwave radiation has water vapour,
  CO2 and ozone only, so another gas's warming enters as the CO2 that warms as much. A calibration pair from
  `A28_dim5_moon`'s year-29 restart (`A28_cal_control`, `A28_cal_co2_1022`, one year each) gives 5.27 W/m² at the top
  of the air for 400 to 1,022 ppm over the first 3-day mean, 5.63 W/m² for each e-fold.
- `--top-infrared DAY NIGHT` adds infrared from space to the downward flux at the top of PlaSim's air, over the
  hemisphere where the Sun is up and the one where it is down, and PlaSim carries it down with its broadband
  transmissivities (`TOP_INFRARED_EDITS` in [exoplasim_run.py](exoplasim_run.py)). At zero the patched build repeats
  the unpatched one exactly (`A28_cal_control_patched`).

Three branches from the same restart flew 20 years each, about 3⅓ minutes a model year on 4 ranks. Their years 10–19
against the design case's years 20–29 are in
[research/studies/gases_glow_dimming](../../research/studies/gases_glow_dimming/README.md):

| | `A28_gases` | `A28_gases_glow` | `A28_gases_glow_held` |
|---|---|---|---|
| CO2; infrared from space by day / night | 755 ppm; none | 755 ppm; 10.17 / 0.52 W/m² | the same |
| Sunlight scale (the design case's 0.95) | 0.933907 | 0.909867 | 0.916522 |
| Surface temperature (design case 294.90 K) | 295.03 K | 293.64 K | 294.84 K |
| Precipitation (2.24 mm/day) | 2.22 mm/day | 2.07 mm/day | 2.16 mm/day |
| Cloud cover (23.5%) | 23.5% | 20.1% | 19.6% |

- **The gases' dimming holds the design temperature** and takes 0.9% of the rain.
- **The glow's infrared heats the upper air.** About 4% of it reaches the ground. Matched watt for watt by dimming, it
  cools the Moon 1.3 K; the dimming that holds the temperature is 0.71 of that. Held there, the air above about 55 km
  warms 0.3–2.1 K, the convective rain falls 3.9% and the high cloud's cover falls by a fifth.

**The glow pair on the corrected radiation**, which the author asked on 16 October to run soon, flew on 17 October
from `A28_corrected_970` with 13 October's added dimming (`A28_corrected_gases_glow`, `A28_corrected_gases_glow_held`).
Their years 10–19 against the design case's, from the runs' yearly means:

| | `A28_corrected_970` | `A28_corrected_gases_glow` | `A28_corrected_gases_glow_held` |
|---|---|---|---|
| Sunlight scale | 0.970 | 0.929066 | 0.935856 |
| Surface temperature | 295.02 K | 293.45 K, cooling 0.1 K a decade | 294.41 K |
| Precipitation | 2.31 mm/day | 2.11 mm/day | 2.19 mm/day |
| Cloud cover | 19.8% | 17.0% | 16.6% |

On this radiation the same dimming cools more. The held branch ends 0.61 K below the design case, where on the old
radiation it ended within 0.06 K. The two branches differ by 0.70 points of dimming and 0.96 K, 1.37 K a point, so
holding the temperature takes 0.45 points less than the held branch: 3.07% of the design case's sunlight, 7.03 W/m² of
its 228.5 W/m² absorbed, for the gases and the glow together. The pair fixes only that total, since no gases-only
branch flew on this radiation. With the gases at their forcing, 3.58 W/m², as on 13 October, the glow takes
3.45 W/m², 0.645 of its 5.34 W/m² of infrared; on the old radiation it took 0.71. Kept at their 1.69% of this
radiation's sunlight, 3.87 W/m², the gases would leave the glow 0.59, with the same total.
[The target-climate study](../../research/studies/climate_targets/README.md#the-heat-beyond-the-sunlight) binds both
to the runs' records. The held branch rains 5.5% less, against 3.6% on the old radiation.
The same morning the longwave was found to cool the lower air too much (below), and the author asked for the pair to
run again with its correction in place, from a design case that lands on the design temperature under it; the lakes
and the climate's other readers follow whatever that changes in the design case's rain.

### The longwave by height

On 17 October the 45° N ring's heat budget turned to the GCM's own heating by process (`A28_corrected_970_heating`,
the runner's `--heating-diagnostics`). On its cloud-free columns PlaSim's longwave cooled the lowest 1.8 km 10.6 W/m²
more than the atmosphere domain's correlated-k model, which follows the line-by-line spectroscopy, and the next layer
11.4 W/m² more ([the CRM's account](../crm/README.md#what-else-divides-cm1-and-the-gcm-at-the-crossing)). The author
asked for it to be refitted by height. [plasim_longwave.py](plasim_longwave.py) copies PlaSim's longwave, Sasamori's
broadband absorptivities with the curved air's escape, and reproduces its own heating on the run's cloud-free 3-day
means within 0.01 K a day on average. Against correlated-k on about 2,100 of those columns, from every latitude:

| Layer from the ground | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | Outgoing |
|---|---|---|---|---|---|---|---|---|---|---|
| PlaSim's error, W/m² | −9.7 | −11.9 | −6.7 | −5.1 | −3.6 | +0.4 | +4.6 | +9.7 | +12.0 | +3.7 |

It cools the lower air too much and the upper air too little, alike at every latitude, and the lowest layer by
16 W/m² over the sea against 6 over land. Its continuum was fitted to the outgoing longwave alone, and its own settings
cannot move the cooling: refitting the continuum, with or without a term on the vapour's own pressure, lowers the
fit's score from 738 to 697, and a power of pressure on the paths to 644. So each layer gets an extra escape to space,
fitted over land and over sea and scaled with its own temperature to the fourth power, as the curved air's escape is,
the top layer's carrying the outgoing longwave to correlated-k's
([longwave_A28_corrected_970_heating.json](../results/gcm/longwave_A28_corrected_970_heating.json)). It runs from
+12 W/m² high up to −6 W/m² in the lowest layer over land and −16 over the sea, and removes the mean error of every
layer over land and over sea, leaving 1–5 W/m² from column to column and within 2.5 W/m² in every band of latitude.
The runner's `--longwave-correction` writes it to `lwcorr.dat` for a patch that adds it after the curved air's escape
(`lw_correction` in radmod.f90). It does not hold near 301 K, and since 19 October the runner takes its refit there
([the longwave correction near 301 K](#the-longwave-correction-near-301-k)).

In the GCM the patch reproduces the copy: a year from the design case with it (`A28_corrected_970_lwcheck`, 2 ranks)
gives PlaSim's own longwave heating within 0.02 K a day of the copy's on its cloud-free columns, the lowest layer now
cooling 0.24 K a day where it cooled 0.45. That year also shows what the correction does to the climate. The top two
levels cool 11 and 7 K, high cloud forms in them, from 20% to 46% of the sky, and its longwave effect rises from 6.6 to
19.7 W/m²; the lower air warms 2–3.6 K, the ground 3.0 K over the year and 4.8 K by its end, and the top of the air
gains 19.6 W/m², still warming. The upper levels' extra cooling is what the correlated-k model gives, and PlaSim's
cloud, fitted to Earth by relative humidity, answers it.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave validate A28_corrected_970_heating --year 0   # a second
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave fit A28_corrected_970_heating --year 0 --every 20   # 2 min
```

### The high cloud and its ice

On 17 October, before the 20-year test of the correction, the author asked for a look at PlaSim's high cloud: its
humidity thresholds, how coarse its top levels are, and CM1 as a check that may itself lean on Earth.

- **Its thresholds sit near where ice appears.** PlaSim reckons saturation over liquid water at every temperature
  (Magnus-Tetens in `rainmod.f90`), for its cloud and for rain-out alike. Cloud grows as ((RH − rcrit)/(1 − rcrit))²
  above rcrit = max(0.85, σ, 1 − σ). In the cold upper air ice saturation is 1.6–1.7 times lower, so the high cloud
  starts at 136–149% of ice saturation and vapour stays in the air up to 160–169%. Clean air makes ice by itself at
  154–167% at those temperatures (homogeneous freezing, Koop et al. 2000, as Kärcher and Lohmann 2002 fit it);
  ice-nucleating particles start cirrus at 110–130%. After the correction the 81–107 km layer sits at 105% of ice
  saturation in the 3-day means and 29% of them pass the cloud threshold: the state of Earth's tropical tropopause
  layer, at the same pressures and temperatures, where thin cirrus is common.
- **Its top layers are Earth's.** The top three are 0–9 kPa (above 107 km, reaching the top of the air), 9–20 kPa
  (81–107 km) and 20–32 kPa (63–81 km), PlaSim's σ levels for Earth at the same pressures. Through the 81–107 km layer
  the temperature runs from about 180 to 204 K and ice saturation changes 55-fold, while one humidity and one cloud
  fraction stand for the whole layer. On Earth these levels lie at the tropopause with warmer air above. Behind the
  titania stack the Moon's air keeps cooling upward, 150–177 K above an 86-km tropopause in the
  [middle-atmosphere model](../../atmosphere/middle_atmosphere/README.md), so the coldest, most saturated air falls in
  the two layers the grid resolves worst.
- **CM1 checks less than it seemed.** The 45° N ring holds ice in every column above 107 km, but that ice comes from
  its start, and below it the ring's air is held to the GCM's vapour
  ([the CRM's account](../crm/README.md#what-else-divides-cm1-and-the-gcm-at-the-crossing)). Where it makes ice of
  its own, at 63–107 km in 10–33% of columns, it holds 1.5–3.5 times the water PlaSim gives its clouds there.
- **The correction's upper part stands.** On the ring's cloud-free columns RRTMG run offline cools 81–107 km by
  0.15 K a day in clear sky where the correlated-k model cools it 0.17, and CM1's own RRTMG, which also carries the
  ring's ice above and its own layers above the top, by 0.13 (its stored tendency read back to RRTMG's rate;
  [the CRM's account](../crm/README.md#what-else-divides-cm1-and-the-gcm-at-the-crossing)). On the GCM's own columns
  PlaSim cools that layer 0.04 K a day where correlated-k cools it 0.19. Its sunlight there was fitted by layer on
  14 October (above).
- **The cloud's optics were liquid water's.** PlaSim gives every cloud Stephens' (1978) fit of reflection to water
  path and 0.1 m²/g of infrared absorption. Ice with the same water reflects more sunlight and, from crystals of about
  15 µm up, traps less infrared.

[plasim_ice_optics.py](plasim_ice_optics.py) gives cloud colder than −38 °C, where cloud water freezes by itself, CAM3's
ice optics as RRTMG carries them: Ebert and Curry's (1992) extinction, single-scattering albedo and asymmetry in
sunlight, through the atmosphere domain's delta-scaled PIFM two-stream (PlaSim's near-infrared band as the sum of
their four bands, weighted by the TSIS-1 sunlight), and their grey infrared absorption, 0.005 + 1/r m²/g. The
effective radius r comes from the layer's temperature through Kristjánsson and Mitchell's table for hexagonal columns,
180–274 K, as CAM3 and WRF's RRTMG carry it. A fully cloudy layer alone, as means over the sunlit disk
([ice_optics.json](../results/gcm/ice_optics.json)):

| Cloud after the correction | Water | PlaSim's reflection; emissivity | Ice: radius; reflection; emissivity |
|---|---|---|---|
| 81–107 km, 189 K | 1.5 g/m² | 1.1%; 0.22 | 9.7 µm; 7.9%; 0.24 |
| 63–81 km, 217 K | 4.7 g/m² | 7.6%; 0.54 | 32.8 µm; 7.3%; 0.24 |

Two stretches lie past what the sources measured. Below 180 K, where the table ends, the radius follows the table's
own trend, 5.5% smaller per kelvin, to 2.6 µm at 165 K, the top layer's temperature after the correction (5% of its
sky clouded in the check year). Ebert and Curry fitted crystals of 10 µm and more, which the table gives from 190 K, so
the cloud of the top two layers rests on their formulas carried to smaller crystals. Small ice absorbs strongly in the
infrared per gram, as the formula's 1/r does. The lunar side may also differ: in a sixth of Earth's gravity crystals
can grow larger before they fall, which would reflect less and trap less. The range for the high cloud, from its water
and its crystals' size, follows the 20-year test.

The runner's `--ice-optics` puts it into the twilight copy (`ICE_OPTICS_EDITS`: `ice_cloud_sw` after PlaSim sets each
layer's liquid optics in `swr`, `ice_cloud_lw` after its cloud transmissivities in `lwr`), switched by `NICEOPT`. The
patch's own Fortran, compiled as PlaSim compiles, reproduces the Python copy to nine figures, and with `NICEOPT` 0 the
rebuilt model repeats the longwave check year exactly (`A28_corrected_970_lwcheck_rebuilt`). A year from the design
case with both the correction and the ice optics (`A28_corrected_970_icecheck`, 4 ranks), against the uncorrected year
and the correction alone, global means:

| | Uncorrected (`A28_corrected_970_heating`) | Correction (`A28_corrected_970_lwcheck`) | Correction and ice |
|---|---|---|---|
| Ground: the year; its last 36 days | 294.99; 295.20 K | 298.04; 299.93 K | 297.46; 298.77 K |
| Top of the air: the year; its last 36 days | −0.2; −0.9 W/m² | +19.5; +22.3 W/m² | +13.1; +15.5 W/m² |
| Cloud effect: sunlight; infrared | −6.4; +6.6 W/m² | −7.1; +19.7 W/m² | −10.4; +16.4 W/m² |
| Cloud cover: all; 81–107 km; 63–81 km | 19.6%; 5.6%; 6.1% | 45.7%; 29.9%; 13.6% | 43.1%; 27.0%; 14.1% |
| Top two levels | 175.6, 196.2 K | 164.5, 189.3 K | 165.6, 189.2 K |

Ice keeps nearly all the cloud, 27% of the sky at 81–107 km against 30%. It reflects 3.3 W/m² more sunlight
and traps 3.4 W/m² less infrared, so the high cloud's net warming over the uncorrected year falls from 12.4 to
5.7 W/m². The top of the air still gains 13 W/m² over the year; the 20-year test with both settles near 307 K
([below](#where-the-corrected-moon-settles)).

```sh
climate/gcm/.venv/bin/python -m climate.gcm.plasim_ice_optics write          # seconds
climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A28 --folder A28_corrected_970_icecheck --years 1 --ncpus 4 \
    --set sunlight_scale=0.970 --longwave-correction --ice-optics --heating-diagnostics \
    --start-from climate/gcm/runs/A28_corrected_970/model/MOST_REST.00019
```

### Where the corrected Moon settles

The 20-year test with the longwave correction and the ice optics (`A28_corrected_970_lw_ice`, 17 October, from
`A28_corrected_970`'s year 19 at the design sunlight scale of 0.970, 4 ranks, about 4.7 minutes a year) settles near
307 K. Over years 10–19 the ground averages 306.5 K, 11.5 K above the design case, and still warms 0.08 K a year, while
the top of the air gains 0.24 W/m² over the last five years. Fitted against the ground's temperature, the gain falls
2.2 W/m² for each kelvin over years 5–19, which puts the balance at 306.8 K, and 1.65 W/m² over years 0–19, at 307.1 K.
Global means over years 10–19 (`progress.json`; the air is at 2 m, its warmest and coldest the extremes of the year's
3-day means, the equator's land within 12°, the poles beyond 70°):

| Years 10–19 | Design case (`A28_corrected_970`) | Correction and ice (`A28_corrected_970_lw_ice`) |
|---|---|---|
| Ground; air | 295.0; 294.3 K | 306.5; 307.0 K |
| Top of the air; planetary albedo | +0.1 W/m²; 0.270 | +0.6 W/m²; 0.278 |
| Cloud effect: sunlight; infrared; net | −6.4; +6.8; +0.3 W/m² | −10.8; +16.8; +6.0 W/m² |
| Cloud cover; rain | 19.8%; 2.31 mm/day | 38.3%; 3.32 mm/day |
| Air: warmest; coldest | 303.3; 287.1 K | 314.9; 297.7 K |
| Equator's land: air; its swing through the year's 3-day means | 294.4; 6.7 K | 306.4; 5.5 K |
| Poles' air | 292.7 K | 307.8 K |

The warmth reaches every latitude. The poles' air warms 15 K and the equator's land 12 K, and the land's surface
poleward of 60° is the warmest, 34.7–34.9 °C against 31.9–32.5 °C within 30° of the equator. The rain rises 43% and
reaches further poleward, and the land under 0.5 mm a day, the desert measure of 14 October, halves, from 37.6% to
18.9% counting every cell's land and from 42.3% to 20.1% on the model's own land cells. The rain is still heaviest
near the equator:

| Land, years 10–19: every cell's land (the model's land cells) | Design case: desert; rain | Correction and ice: desert; rain |
|---|---|---|
| 0–15° | 0%; 4.88 mm/day (0%; 4.62) | 0%; 6.09 mm/day (0%; 5.45) |
| 15–30° | 1%; 2.87 mm/day (1%; 2.78) | 1%; 4.22 mm/day (1%; 3.98) |
| 30–45° | 34%; 1.52 mm/day (38%; 1.42) | 3%; 3.01 mm/day (2%; 2.93) |
| 45–60° | 85%; 0.24 mm/day (86%; 0.23) | 22%; 1.22 mm/day (18%; 1.27) |
| 60–75° | 100%; 0.06 mm/day (100%; 0.06) | 93%; 0.17 mm/day (93%; 0.17) |
| 75–90° | 100%; 0.04 mm/day (100%; 0.04) | 100%; 0.04 mm/day (100%; 0.04) |
| Below 2 km; 2–4 km; 4 km and up | 46%; 36%; 4% (47%; 39%; 9%) | 24%; 16%; 1% (23%; 18%; 2%) |
| All land | 37.6%; 2.25 mm/day (42.3%; 1.98) | 18.9%; 3.34 mm/day (20.1%; 2.94) |

The rain falls on a cell's land as well as on its water, so the land's rain and deserts count every cell's land
fraction and area with that cell's own rain, and the model's own land cells alone give the dry bound. The 18% of the
boundary file's land that lies in cells the model runs as water is mostly coast and lakeshore, and those cells rain
5.15 mm a day here. The T21 grid cannot resolve the coasts, and coast and island hydrology, geography's groundwater
work and the shore study to come, decides where between the two the real coasts lie.

The heat at the ground is humid. In the 3-day means, each about a tenth of the lunar day, the land's air at 2 m has a
median of 33.6 °C (20.4 in the design case), 39.3 °C at the 99th percentile and 42.3 °C at most. Its wet-bulb
temperature at 2 m has a median of 26.9 °C (16.6), 31.5 °C at the 99th percentile (21.8) and 34.1 °C at most (24.0),
and passes 28 °C for 31% of the land's time and 31 °C for 2.1%. At the lowest level, about 0.9 km above the ground,
the air's median is 32.1 °C (19.0) and the wet-bulb temperature's 26.3 °C (15.9), with 30.9 °C at the 99th percentile
(21.1), passing 28 °C for 22.3% of the land's time and 31 °C for 0.8%. The seas' surface has a median of 34.5 °C
(24.1) and 36.8 °C at the 99th percentile. The air cools about 1 K a kilometre upward, and the median of the land's
columns crosses freezing 36 km above the ground (23 km in the design case).

The heat's land is the model's own, the cells it runs as land. The runner turns cells to water, least land first,
until the planned water share is met (`land_mask` in [exoplasim_run.py](exoplasim_run.py)), so 166 cells with a land
fraction of 0.50–0.59 run as water and water covers 40.05% of the Moon against the planned 40.02%. The land's heat is
weighted by land area and the seas' by water area, their quantiles included. The output's air at 2 m is the lowest
level's carried down the dry adiabat to the ground's pressure (`mkshfl` in PlaSim's fluxmod.f90): over land it stands
1.5 K above the lowest level's air and within 0.02 K of that air times σ^−0.286 (standard deviation 0.07 K). The
wet-bulb temperature at 2 m carries the lowest level's specific humidity down with it, capped at saturation, which
binds in none of the land's means here and 0.002% in the design case. Together they read a well-mixed lowest
kilometre, as by day. The output's humidity is a snapshot: ExoPlaSim 3.4.2 computes it only at the steps where
mod(nstep, nafter) is zero (plasim.f90, line 3218) and adds that held value at every step of the mean (outmod.f90,
line 2556), so each record holds the humidity at its mean's first step, while the air, the pressure and the Sun's
cosine are true means. Each mean's humidity is therefore the mean of its start and its end, the next record's, and
the last mean of years 10–19, with no next record, leaves the heat's figures. The wet-bulb temperature is the
psychrometric one over water with Bolton's saturation, capped at the air's own temperature where the air is
supersaturated: with the humidity as written, 1.7% of the design case's land means at the lowest level, and with
each mean's own, 0.24%. The figures come from the run's output files and the T21 boundary file through the
target-climate study's [extract](../../research/studies/climate_targets/extract.py), and its product carries them.
The first readout of 17 October took the humidity as written and counted every cell with a land fraction above one
half alike, and put the lowest level's time above 28 °C at 17.6%.

Place by place, in the same 3-day means of years 10–19, with the design case's figures in brackets:

| Land | Wet-bulb at 2 m: median; 99th percentile | Time above 28 °C wet-bulb: at 2 m; 0.9 km up | Air at 2 m: median |
|---|---|---|---|
| 0–15° | 27.9; 32.0 °C (18.3; 22.4) | 49%; 38% (0) | 33.5 °C (21.0) |
| 15–30° | 27.8; 31.6 °C (17.9; 21.8) | 46%; 34% (0) | 33.5 °C (20.9) |
| 30–45° | 27.2; 31.1 °C (16.8; 21.2) | 34%; 23% (0) | 33.3 °C (20.6) |
| 45–60° | 26.4; 30.3 °C (15.9; 18.8) | 17%; 9.6% (0) | 33.3 °C (20.2) |
| 60–90° | 25.6; 28.8 °C (15.1; 17.7) | 3.6%; 1.4% (0) | 34.5 °C (19.6) |
| Below 1 km | 28.3; 32.1 °C (17.8; 22.3) | 56%; 44% (0) | 34.6 °C (21.6) |
| 4 km and up | 25.4; 28.8 °C (16.0; 19.6) | 4.6%; 1.7% (0) | 30.7 °C (18.9) |
| Coasts, the land beside the sea | 27.4; 31.6 °C (17.1; 21.8) | 38%; 28% (0) | 33.8 °C (20.7) |

Each cell's warmest and coolest tenth of its 3-day means, ranked by one quantity with the other over the same means,
averaged over the land at 2 m (the design case in brackets):

| Land, at 2 m | Wet-bulb | Air | In the day; its afternoon |
|---|---|---|---|
| Warmest tenth by the wet-bulb temperature | 29.5 °C (19.2) | 34.1 °C (23.5) | 68%; 56% (85%; 67%) |
| Coolest tenth by the wet-bulb temperature | 24.0 °C (14.1) | 33.2 °C (19.5) | 34%; 7% (27%; 5%) |
| Warmest tenth by the air | 27.5 °C (18.4) | 36.9 °C (24.8) | 98%; 70% (98%; 80%) |
| Coolest tenth by the air | 26.4 °C (15.4) | 31.7 °C (18.5) | 19%; 3% (18%; 1%) |

The air's warmest tenth falls mostly in its afternoons, and the wet-bulb temperature's mostly in the day, 56% of it in
the afternoons. The wet-bulb temperature's coolest tenth comes in air near the median's warmth: in the hot world the
wet-bulb temperature follows the humidity more than the air. In the air's coolest tenth, mostly nights, the 2 m air
stands 1.1 K above the ground's own temperature (1.8 K in the design case); in calm night air the air near the ground
cools below the lowest level's air carried down, as the cloud-resolving ring's calm, saturated land nights show ([the
CRM's account](../crm/README.md#what-the-ring-shows)), so the 2 m figures read warm at night. The seas' surface has a
median of 35.3 °C within 15° of the equator and 33.2 °C poleward of 60° (25.8 and 20.4). Over the land, in land-area
means at each height above the ground, the air 10 km up, 22.3 °C at 0.96 atm, is near the design case's at the lowest
level, 19.6 °C at 1.14 atm; 20 km up it is 12.9 °C at 0.79 atm and 30 km up 5.1 °C at 0.65 atm (2.6 and −6.9 °C at the
same heights in the design case).

The base of the thermosphere barely follows the ground: the middle-atmosphere cases at ground temperatures of 288 and
308 K ([results/cases](../../atmosphere/middle_atmosphere/results/cases/), titania stack at 1.2 atm) put the 0.3 Pa
base at 163.3 and 163.4 K. Through the escape model's thermal column at solar maximum, behind the film alone, the
exobase warms from 163.6 to 163.7 K and the escape rises from 4.4 × 10⁻⁸ to 5.8 × 10⁻⁸ kg/s, an atmospheric cycle
time of about 1.5 × 10¹⁸ years; with 0.1% of the far ultraviolet transmitted, from 3.2 × 10⁻³ to 3.8 × 10⁻³ kg/s,
about 2.3 × 10¹³ years (`exobase_warm.py` in `session-tools/climate_settled_2026-10-17` on the research drive). Which
target the Moon takes, the hot world it settles on here, the resort between or the design temperature, the lock, is
the author's open question ([register](../../research/decisions.md#the-target-climate-17-october-2026)); new runs wait
on it.

### The curved air's light

On 14 October the author asked whether the GCM and CM1 take enough account of the Moon's scattering. Both light an
air column only while its Sun is above that column's horizon, and both treat each column as flat layers; PlaSim also
stretches the low Sun's water path with Earth's horizon magnification, 35 at the horizon. The Moon's air is tall
against its radius: its mass sits on average 48.6 km up, and a shell 50 km up has 5.8% more area than the ground
under it. [twilight.py](twilight.py) sets what that geometry does to the sunlight and the infrared, against flat
columns of the same air, and writes [../results/gcm/twilight.json](../results/gcm/twilight.json):

```sh
PY=/media/projectspace/terluna-research/venvs/optics/bin/python     # the sky model's solver needs numba
$PY -m climate.gcm.twilight optics      # the near infrared's line optics, cached on the drive
$PY -m climate.gcm.twilight longwave    # the longwave line optics, cached
$PY -m climate.gcm.twilight solve       # the sky's multiple scattering in 12 batches, 40 minutes on 6 threads
$PY -m climate.gcm.twilight ledger      # the product, 2.5 minutes
# each with --ts 301 for the design climate's column: twilight_301k.json, about 50 minutes on 4 threads in all
# (gcm-runs/twilight/ledger_301k.sh on the research drive)
```

**How the ledger is built.** The design column's line-by-line optics from 360 nm to 5 µm are grouped into 10,697
channels that keep the absorbed energy within 1% on slant columns; the 10,415 above 1,250 nm keep, within the same 1%,
what the beam loses above the bottom of every layer as well. The 282 channels below 1,250 nm that are not opaque
take the spherical sky model's multiple scattering at its standard quality; the rest take the direct beam through the
same 100 spherical shells with the diffuse light of a pseudo-spherical two-stream column. The flat and
pseudo-spherical columns are the radiative-convective model's two-stream on the same layers, and PlaSim's is
[plasim_shortwave.py](plasim_shortwave.py) on the solved column at the design run's levels and calibration. The
infrared leaving the top is traced through the same layers, flat and as shells, whose rays above the ground's limb
cross the air twice; its flat side matches the radiative-convective model's own solver within 0.7%, and a thin shell
at 46 km lets out 1.2250 times its flat escape against the exact 1.2253. All of it is clear sky over ground of albedo
0.1, without refraction.

**The sky model's channels hold its colours and overstate its absorbed energy.** Its 114 production channels average
optical depth over wavenumbers that mix windows with opaque water and O₂ lines, and an averaged depth absorbs more
than its members: with the Sun overhead they absorb 104.5 W/m² above 500 nm against the fine grid's 66.2, 58% too
much. The ledger's channels agree with the fine grid within 0.4–0.8% in total and, with the Sun overhead, within
0.9 W/m² in every GCM layer.

**Channels have to keep the heights as well.** The ledger's first channels were bisected on the whole column alone.
Narrow line cores that go dark high up and wings that absorb low down can take the same share of the whole column, and
a channel that mixes them puts its sunlight too high: above 1,250 nm, with the Sun overhead, those channels put
19.0 W/m² at 66–85 km against the fine grid's 10.7 and 27.7 at 35–49 km against 35.1, with the column's total right
within 0.2%. The test of the beam's loss above every layer corrects it with ten times the channels. The figures here
are the corrected ledger's; the two 20-year tests below ran with the first one's terms.

**What the geometry does**, global means at the design column's temperatures, W/m², against the pseudo-spherical
column of the same air, whose totals PlaSim's calibrated arithmetic matches within 1.1 W/m²:

| | Sunlight to the ground | Sunlight to the air | Infrared to space | Net |
|---|---|---|---|---|
| Cone: the shells widen above each square metre | 0 | +4.0 | +10.5 more | −6.5 |
| Curved light and curved escape | +1.5 | +4.8 | +5.6 more | +0.7 |
| All | +1.5 | +8.8 | +16.2 more | −5.8 |

- **The twilight is long and bright.** The curved light gives 2.6 W/m² with the Sun below the horizon (1.4 to the
  ground) and 3.0 W/m² with it 0–10° up, where the beam crosses each column from the side and the flat column counts
  only what enters its top: at 1° the air takes 33.9 W/m² per square metre of ground against 7.5. Overhead the
  sphere's light on the flat column's layers matches the flat column within 0.1% in the air; its ground takes 1.3%
  less, the exact scattering reflecting more than the two-stream.
- **The infrared leaves by the same geometry.** Most of it leaves from 15–100 km, where the shells are widest, and the
  upper air also sees space below its horizontal down to the ground's limb. At the same temperatures the shells lose
  16.2 W/m² more than flat layers, which outweighs the 10.3 W/m² more sunlight they take.
- **Almost none of it reaches the thermosphere.** Of the curved light's 4.8 W/m² in the air, 1.0 falls below 15 km,
  2.3 at 15–50 km and 0.14 above 150 km; 0.0007 W/m² falls above the thermosphere's 0.3 Pa base, 519 km up. The
  exobase tables run on millionths of a W/m², and that light comes in CO₂ and O₂ line cores whose molecules also cool
  and, at those pressures, re-emit much of what they absorb; the column here holds them in local equilibrium and stops
  at 0.1 Pa, so the middle-atmosphere model is where to judge it.
- **PlaSim places its absorbed sunlight too low.** Its column's total matches; by layer it puts nothing in the top two
  (above 85 km), where the line-by-line column takes 4.2 W/m², 2.3 W/m² less in the third (66–85 km) and 2.1 W/m²
  more in each of the fourth to sixth (24–66 km).
- **The Moon's horizon slant factor is 13.5 for its water and 7.1 for well-mixed air**, against PlaSim's 35. The
  rational form (μ + 0.18401) / (μ² + 0.17006 μ + 0.013593) follows the water's within 0.30% at every elevation and is
  exactly 1 overhead. With it in PlaSim's arithmetic the calibration's water multiplier goes from 1.3564 to 1.3619 and
  its Rayleigh multiplier from 2.2146 to 2.2160 to keep the daily-mean split it was fitted to.

**Into the models.** `twilight.json` holds the terms by solar elevation for the GCM's ten layers and the ground
(`gcm_term`, the curved light, and `gcm_term_cone`) and for the design column's 100 shells (`cm1_term`,
`cm1_term_cone`), and the infrared by layer and shell with the temperatures it is given at (`longwave`). The runner's
`--twilight` puts them into PlaSim through `TWILIGHT_EDITS`: at the end of swr the sunlight by each column's solar
elevation, solang keeping the cosine below the horizon, the ground's share in the column's clear share; at the end of
lwr each layer's and the ground's extra infrared, scaled with its own temperature to the fourth power; and the Moon's
slant factor in place of Earth's. `curved` takes the curved light and escape, `full` the cones as well, `path` the
slant factor alone. The patch goes only into a copy of ExoPlaSim of its own on the drive, which those runs import
first, so the shared installation and every other run's key stay as they were. CM1 takes the same terms through the
`moon_omp_twilight` build
([the CRM's account](../crm/README.md#the-45-rings-on-the-design-ground-with-the-curved-airs-light)).

**The 20-year test.** `A28_dim5_moon_geometry` branches from the design case's year-29 restart with `--twilight full`,
the refitted water and Rayleigh multipliers and the design case's settings otherwise (sunlight scale 0.95), and flew
20 years on 4 ranks, about 3¾ minutes a model year. Its years 10–19 against the design case's years 20–29:

| | Design case | With the curved air | Change (year-to-year spread) |
|---|---|---|---|
| Surface temperature | 294.90 K | 294.24 K | −0.65 K (0.04) |
| Air at 2 m | 294.37 K | 293.58 K | −0.79 K (0.03) |
| Cloud cover | 23.5% | 26.2% | +2.8 points, nearly all in the high deck |
| Clouds' net effect at the top | −1.6 W/m² | +2.2 W/m² | +3.8 W/m² |
| Precipitation | 2.24 mm/day | 2.20 mm/day | −1.6% |

- **The clouds give back most of what the geometry takes.** In the first year the clear sky lost 5.4 W/m² at the top
  of the air, as the ledger has it, and the clouds returned 3.5 W/m², reflecting less sunlight and holding more
  infrared; the climate settled within about eight years.
- **The evenings keep some warmth.** Land air cools 0.8–0.95 K at most local times, but at dusk at 30–60° only
  0.3–0.6 K, the twilight holding the evening air up; dawn gains a little. The 3-day output means resolve local time
  to about a tenth of the lunar day.
- **The rain moves poleward**: 0.19 mm/day less over land within 15° of the equator, 0.02–0.05 mm/day more over land
  at 30–60°.
- **The sunlight side alone would mislead.** A first branch with only the curved sunlight
  (`A28_dim5_moon_twilight`, stopped after four years) warmed 1.1, 2.1, 2.9 and 3.7 K in its first four years.

**PlaSim's CO₂ and O₂.** PlaSim's sunlight has two bands: below 0.75 µm Rayleigh scattering and ozone, above it
water vapour alone, Lacis and Hansen's absorptivity of a pressure-weighted water path. CO₂ and O₂ are not in either,
and the calibration's water multiplier absorbs their share where the water is. The Moon's column holds 7.2 times
Earth's air over each square metre, so their strongest lines go dark high up: with the Sun overhead CO₂ alone takes
6.1 W/m² in the GCM's top layer (112–575 km), where PlaSim takes nothing. `plasim_shortwave.py gases` fits each band
an absorptivity a u / ((1 + b u)^c + d u) of its own path u, each layer's mass of the gas weighted by (p / 1 bar)^n,
to the line-by-line direct beam's loss above each of the column's 100 levels at slant factors 1 to 13
([../results/gcm/shortwave_gases.json](../results/gcm/shortwave_gases.json)):

| Band | Overhead Sun takes | n | Largest error, share of the band |
|---|---|---|---|
| CO₂ near infrared (2.0, 2.7 and 4.3 µm) | 22.2 W/m² | 0.83 | 0.06% |
| O₂ near infrared (A band, 1.27 µm) | 8.9 W/m² | 0.95 | 0.13% |
| O₂ visible (B and γ bands) | 3.4 W/m² | 0.88 | 0.04% |

The runner's `--gas-bands` puts them into the twilight copy through `GAS_EDITS`, each multiplied into its PlaSim band's
layer transmissions as water's are, along the Moon's slant path for well-mixed air (7.1 at the horizon), with the gases'
mass fractions from the run's air (CO₂ 6.1×10⁻⁴, O₂ 0.194). With the bands and the Moon's slant factors the
calibration refits to a water multiplier of 1.292 and a Rayleigh multiplier of 2.237, keeping the air's daily-mean
share within 0.014 points on average and 0.24 rms
([../results/gcm/shortwave_fit_A28_dim5_lunar.json](../results/gcm/shortwave_fit_A28_dim5_lunar.json)). The design column's absorbed sunlight by layer, global means, W/m², the line-by-line pseudo-spherical column against
PlaSim's arithmetic on the same column:

| GCM layer | Line-by-line | PlaSim as run | PlaSim with the bands |
|---|---|---|---|
| 112–575 km | 2.52 | 0.00 | 2.19 |
| 85–112 km | 1.69 | 0.01 | 1.25 |
| 66–85 km | 4.49 | 2.24 | 3.12 |
| 49–66 km | 11.95 | 14.06 | 13.87 |
| 35–49 km | 16.53 | 18.60 | 17.81 |
| 24–35 km | 16.56 | 18.61 | 17.60 |
| 14–24 km | 14.71 | 15.88 | 14.92 |
| 7–14 km | 11.79 | 12.36 | 11.56 |
| 2–7 km | 8.11 | 8.16 | 7.61 |
| 0–2 km | 3.19 | 2.73 | 2.54 |
| Whole column | 91.54 | 92.65 | 92.47 |

With the bands PlaSim's layers come within 0.97 W/m² rms of the line-by-line column's, against 1.71 W/m² as run. The
rest of the gap is water's, about 2 W/m² that PlaSim puts at 35–66 km in place of 66–85 km and the lowest 7 km: with
the Sun overhead the line-by-line water takes 9.2 W/m² at 66–85 km, just under the design column's 86 km tropopause,
and PlaSim's 5.2. A pressure exponent on PlaSim's water path does not move it (`plasim_shortwave.py water`,
[../results/gcm/shortwave_water_A28_dim5.json](../results/gcm/shortwave_water_A28_dim5.json)): a lower exponent counts
more of the upper water, but most of what it adds falls at 49–66 km, where PlaSim already takes too much, and every
exponent below 1 fits the layers worse. Lacis and Hansen's curve of growth takes sunlight from the first water it
meets, while the line-by-line water's strongest lines go dark higher up and leave less for the layers below; and
PlaSim's water multiplier also stands in for what its bands leave out, water's visible bands and the
collision-induced bands low in the column. In the middle-atmosphere column the misplaced sunlight is worth
0.88 K of ground as run and 0.19 K with the bands
([its account](../../atmosphere/middle_atmosphere/README.md#the-curved-airs-light-2026-10-14)).

**The 20-year test with the bands.** `A28_dim5_moon_geometry_gases` branches from the same restart with
`--twilight full --gas-bands` and the multipliers refitted with them. Its years 10–19 against the design case's 20–29
and the curved air's test without the bands:

| | Design case | Curved air | Curved air with CO₂ and O₂ |
|---|---|---|---|
| Surface temperature | 294.90 K | 294.24 K | 292.02 K, still cooling 0.3 K a decade |
| Air at 2 m | 294.37 K | 293.58 K | 291.57 K |
| Cloud cover; high deck; low deck | 23.5%; 14.3%; 5.3% | 26.2%; 16.2%; 5.3% | 20.0%; 11.6%; 5.3% |
| Clouds' net effect at the top | −1.6 W/m² | +2.2 W/m² | +0.5 W/m² |
| Precipitation; over land | 2.24; 2.04 mm/day | 2.20; 2.01 mm/day | 2.07; 1.90 mm/day |
| Land under 0.5 mm/day; over 2 mm/day | 42.2%; 39.2% | 41.1%; 39.8% | 43.4%; 37.6% |
| Water vapour | 176 kg/m² | 162 kg/m² | 138 kg/m² |
| Air in the top layer, 112–575 km | 170.0 K | 169.3 K | 176.6 K |

- **The bands' heating clears the high cirrus, and the Moon cools 2.2 K more.** The top layer warms 7 K and its cloud
  falls from 2.0% to 0.05%, and the next layer's, at 85–112 km round the 86 km tropopause, from 9.4% to 3.7%. The
  cirrus held infrared: the clouds' longwave effect falls from 8.9 to 6.8 W/m². The column without clouds moved its
  ground 0.57 K with the bands, so most of the cooling is this response of PlaSim's diagnostic, Earth-tuned clouds.
- **The middle air cools about twice as much as the ground**: 5.0–5.7 K at 35–85 km against the design case.
- **The ground cools evenly.** Land air is 2.6–2.95 K cooler at every latitude and time of day, least at noon in the
  tropics and at dusk, most in the small hours; dusk at 30–60° keeps the twilight's warmth, 2.68 K cooler against
  2.9 K at night.
- **The rain falls 7.6%**, over land 7.5% in the tropics and 5.5% at 30–60°: 2.6% for each kelvin against the
  dimming runs' 3.1%, the ground taking 5.6 W/m² more sunlight through thinner cloud and drier air.
- **Holding the design temperature would take about 1.5–1.7% less dimming**, at the dimming runs' 1.72–1.9 K for
  each 1% of sunlight: a sunlight scale near 0.965 in place of 0.95, run below.

**The new standard and the dimming that holds the design temperature.** With the curved air, the bands and the corrected
ledger's terms as the runner's defaults, `A28_corrected_9655` branches from the design case's year-29 restart and
`A28_corrected_9685` from its year 19, and the new design case `A28_corrected_970` from that run's year 19. Their years
10–19 against the old design case's 20–29:

| | Old design case | Scale 0.9655 | Scale 0.9685 | Scale 0.970, the design case |
|---|---|---|---|---|
| Surface temperature | 294.90 K | 294.40 K | 294.76 K, drifting 0.004 K a decade | 295.02 K, drifting 0.01 K a decade |
| Air at 2 m | 294.37 K | 293.70 K | 294.02 K | 294.27 K |
| Cloud cover; in the top two layers | 23.5%; 4.2% and 8.5% | 20.1%; 0.2% and 5.3% | 19.9%; 0.2% and 5.5% | 19.8%; 0.2% and 5.6% |
| Precipitation | 2.24 mm/day | 2.25 mm/day | 2.29 mm/day | 2.31 mm/day |
| Sunlight taken by the ground; by the air | 120.6; 90.6 W/m² | 126.1; 101.1 W/m² | 126.1; 101.7 W/m² | 126.3; 102.2 W/m² |
| Water vapour | 176 kg/m² | 163 kg/m² | 167 kg/m² | 170 kg/m² |

- **The design temperature needs a sunlight scale near 0.970**, a dimming of about 3.0% below the shield product in
  place of the design case's 5%. The two runs move 1.13 K for each 1% of sunlight, against 1.7–1.9 K in the dimming
  runs of the old radiation, and their cloud hardly changes between them.
- **The rain comes back ahead of the temperature.** At 0.9685 the Moon rains 2.5% more than the design case while
  0.14 K cooler at the ground: the ground takes 5.6 W/m² more sunlight under the thinner high cloud.
- **`A28_corrected_970` is the design case since 14 October**, the author having accepted it: 295.02 K at the ground
  and 294.27 K at 2 m, 3.0% dimmer than the shield product, with 3.3% more rain than the old design case. It replaces
  `A28_dim5_moon` as the source of the climatology and of the CM1 rings' forcing
  ([the CRM's account](../crm/README.md#the-45-rings-on-the-corrected-radiations-design-case)).
- **Its land is wetter.** Over land it rains 2.10 mm/day against 2.04; 40.3% of the land gets under 0.5 mm/day against
  42.2%, and 41.7% over 2 mm/day against 39.2%. The gain lies poleward of 15°: 4% more rain at 15–30°, 12% more at
  30–45° and 21% more at 45–60°, while the tropics keep theirs (`climatology A28_corrected_970:10-19`). The products
  that read the old design case's climatology (geography's lakes and groundwater, the surface light, the canopy models)
  still read `A28_dim5_moon`'s until they are regenerated.

### The design case at 301 K: its preparation (19 October)

The resort at 301 K is the design climate since 19 October
([register](../../research/decisions.md#the-design-climate-and-the-waters-19-october-2026)). Its design case waits for
the author's go-ahead with its machine time. On the same day the author approved, as its preparation: a T21 surface with
today's lakes for a cold start, the cold start's spin-up sized from the runs' records, the longwave correction checked
near 301 K, the curved air's ledger rebuilt on a 301 K column with its terms in the clear share, and the surface stress
in the output. A check of the world-crafting lab's copy of this code (`session-tools/world_crafting_check_2026-10-19.md`
on the research drive, items 1–4 and 7) raised these points. Nothing here runs the GCM.

#### The lakes for a cold start

The T21 surface every A28 run read, `moon_28pct_water_lakes_gaussian_T21.nc`, was built on 25 September from the lakes
of run A's runoff. It records the drainage product `060219b12f480619`. That product was overwritten under its name on 29
September by the `A28_dim5_moon` routing (`98842fc2a648f889`), and no copy survives, so the surface cannot be
regenerated from what it records. It stays as it is, the record of the A28 runs. `python -m climate.gcm.boundary
--lakes-only` builds the surface from today's drainage and atlas, as main's geography gives them, under a name taken
from the drainage product's content: `moon_28pct_water_lakes_98842fc2a648f889_<grid>.nc`. A later routing writes a
surface of its own beside it. Its attributes name the surface it supersedes and say that that surface's source is lost.
The products are not committed, so
[boundary_lakes_98842fc2a648f889.json](../results/gcm/boundary_lakes_98842fc2a648f889.json) records the inputs and files
by sha256, the shares and every T21 cell that changes. Copies of the drainage product and of both surfaces, under names
that carry their content, are on the research drive (`gcm-runs/boundary`, with `SHA256SUMS`).

| | Run A's lakes (25 September) | The current routing's lakes (`A28_dim5_moon`'s runoff) |
|---|---|---|
| Lakes; seas | 12.02%; 28.00% | 10.34%; 28.00% |
| Water on the T21 grid | 40.02% | 38.34% |
| PlaSim's binary mask: land; water cells of 2,048 | 59.95%; 759 | 61.62%; 734 |

The mask changes in 111 cells, 7.1% of the Moon's area: 68 turn from water to land and 43 from land to water. Within 15°
of the equator 22 cells become water and 8 land, as the tropical lakes grow; at 15–60° 60 become land and 15 water, as
the lakes there shrink; poleward of 60° 6 become water. The runner's experiment `A28R` reads the new surface, and a cold
start is needed, since a restart carries PlaSim's mask and orography.

#### How long a cold start takes

Three runs started from ExoPlaSim's own initial state, 250 K at the ground with the seas 87% frozen: A and B for 40
years and A28 for 30. The warm branch with the longwave correction and ice optics, `A28_corrected_970_lw_ice`, shows how
the present physics relaxes after a step of 11.5 K. A run counts as settled when the linear trend of its global ground
temperature over the last ten years is under 0.1 K a decade and the mean imbalance at the top of the air over the last
five is under 0.2 W/m². Exponential fits to each run's approach, from two years after its sea ice has gone, give the
year each criterion passes (`gcm-runs/spinup_2026-10-19/spinup.py` and its `spinup.json` on the research drive):

| | A, 25% water | B, 35% water | A28, 28% seas and lakes | Present physics (`_lw_ice`, warm branch) |
|---|---|---|---|---|
| Sea ice gone by | year 3 | year 4 | year 4 | none at the start |
| E-folding of the ground; of the imbalance | 6.1; 5.3 years | 7.1; 6.1 years | 6.3; 5.7 years | 4.0; 4.6 years |
| Its sea's surface; its land's ground | 6.0; 6.2 years | 6.9; 7.2 years | 6.2; 6.4 years | 4.1; 3.9 years |
| At the record's end: ten-year trend; five-year imbalance | 0.16 K a decade; 0.14 W/m² (year 39) | 0.47; 0.33 (39) | 1.5; 0.79 (29) | 0.75; 0.24 (19) |
| The fit's year: trend under 0.01 K a year; imbalance within 0.2 W/m² of its settling value | 36; 26 | 43; 31 | 40; 32 | 22; 20 |
| The imbalance's settling value | +0.11 W/m² | +0.27 W/m² | +0.30 W/m² | −0.25 W/m² |

- **None of the records reaches both criteria.** Run A comes closest, its imbalance passing from year 34 and its trend
  still 0.16 K a decade at year 39. A ten-year trend lags the rate at its middle by about five years and a five-year
  mean by about two, so the criteria pass that much after the fit's years.
- **The imbalance settles a few tenths away from zero**, +0.11 to +0.30 W/m² in the cold starts and −0.25 in the warm
  branch, the model's own residual (the settled design case `A28_corrected_970` holds +0.1 W/m² with a spread of 0.2
  from year to year). The criterion of 0.2 W/m² is therefore read against each run's own settling value as well as
  against zero.
- **No part of the model is slower than the whole.** The 50 m slab sea and the land's ground relax together with the
  global mean, and so does the soil's water (6.9 years in A, 2.3 in A28). PlaSim's soil is 12.4 m deep in five layers
  with a diffusivity of 7.5 × 10⁻⁷ m²/s (`landmod.f90`), whose slowest mode decays in about 2.6 years, and the output
  holds no deep soil temperature to follow. The sea ice is gone within four years. What sets the pace is the heat of the
  slab, the air and its vapour against a weak feedback.
- **The present physics relaxes faster**, 4.0–4.6 years against 5.3–7.1, as its stronger feedback (1.65–2.2 W/m² a
  kelvin, [where the corrected Moon settles](#where-the-corrected-moon-settles)) would make it. Taken from A28's state
  at the end of its ice, 282 K at year 4 with an imbalance of 25 W/m², a cold start on it passes the imbalance criterion
  near year 28 and the trend criterion near year 33 (the fits' years 26 and 29 and the criteria's lags); with the old
  physics' time scale, near years 34 and 45.

Machine time follows from the known rate, 20 years in about 1.5 hours on four ranks: 4.5 minutes a year for the 20-year
test with the correction, the ice optics and the heating diagnostics, 3.6 for the design case without them. A cold start
that settles by year 33 and then gives the ten years the design case's statistics take needs about 45 years: 3.4 hours
on four ranks, or up to 4 hours if its cold first years run slower, as run A's did.

#### The run's plan

Proposed for the author's go-ahead; nothing has run.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.boundary --lakes-only    # the surface, in the checkout the run starts from
climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A28R --folder A28R_301k --years 15 --ncpus 4 \
    --set sunlight_scale=0.925 --longwave-correction --ice-optics --heating-diagnostics
# then the same with --years 30 and --years 45, each block resuming where the last stopped
```

- **A cold start on the new surface**, from ExoPlaSim's own initial state, with the physics of the 20-year test and the
  standards of 19 October: the curved air's terms from the 301 K ledger in the clear share, the CO₂ and O₂ bands, the
  longwave refit at 301 K and the ice optics, all of them the runner's defaults or the flags above. The runner applies
  the clear-share patch to the twilight copy of ExoPlaSim at the run's start, and ExoPlaSim rebuilds PlaSim, a few
  minutes.
- **45 years on four ranks**, about 3.4–4 hours, in three blocks of 15 years (70–80 minutes each) with half-hour rests
  between them. Every model year writes its restart and the progress record, so a crash or a stop costs at most the year
  in progress, about 5 minutes; the runner refuses to resume if anything that sets the model has changed.
- **The sunlight scale, 0.925.** The target-climate study puts 301 K at 0.916–0.927 for the 20-year test's physics on
  run A's lakes and the 294.9 K terms ([the dimming at every
  temperature](../../research/studies/climate_targets/README.md#the-dimming-at-every-temperature)). The design run
  departs from that physics four ways: 1.7 points less water, about −0.2 K at the routings' 0.1 K a point; the 301 K
  ledger, 1.5 W/m² more sunlight in the air and 0.4 W/m² more infrared at the test's state, about +0.5 K; the clear
  share, −1.9 to +2.5 W/m², −0.9 to +1.2 K, its lower end for thin ice aloft; and the longwave refit, whose clear-sky
  outgoing longwave stands about 0.7 W/m² above the 295 K fit's at the test's state, about −0.3 K before its clouds
  answer. 0.925 lies between the study's middle and the thin-ice case; the landing is uncertain by about a kelvin.
- **A check at years 20–25.** By then a cold start on this physics lies within about 0.1–0.5 K of where it settles, and
  `compare --settle` gives that value. If it lies more than 0.3 K from 301 K, the remaining years branch from that
  year's restart on the same coastline with the scale moved 0.73–0.96% for each kelvin (the target-climate study's
  rate), which settles within about 15 years.
- **Output**: the runner's 3-day means with the surface stress, the clear-sky fluxes and the heating by process, about
  86 MB a year, 3.9 GB in all on the research drive. The design years are the ten after the criteria pass.
- **Its alternatives, the author's to weigh**: a branch from the 20-year test on run A's lakes, 20 years in 1.5 hours,
  keeping 1.7 points more water as a stated departure; or a cold start from a warmer initial state, saving perhaps ten
  years. PlaSim's initial profile takes its ground temperature from its namelist's `TGR`, 288 K by default; how far the
  slab sea and its ice follow it is untried here.

#### The longwave correction near 301 K

The correction of 17 October ([the longwave by height](#the-longwave-by-height)) was fitted on the design case's
cloud-free columns near 295 K, and the GCM scales each layer's term with its temperature to the fourth power alone. No
run has settled at 301 K, so `plasim_longwave.py check` takes the 20-year test with the correction and ice optics
(`A28_corrected_970_lw_ice`) as it passes 299.6, 301.2 and 302.5 K (its years 1–3, warming about 1.5 K a year) and where
it settles near 306.6 K (year 15), with the fit's own year at 295 K beside them. Those years carry the correction's own
upper air, high cloud and vapour, which an interpolation between the two settled runs would not. On each it sets
PlaSim's copy, with and without the correction, against the correlated-k model: up to 1,500 columns clear of cloud, as
the fit took them, and every 193rd column of the year whatever its cloud, both codes clear-sky
([longwave_check_301k.json](../results/gcm/longwave_check_301k.json), with the samples cached on the research drive
under `gcm-runs/longwave_check`). The largest of the nine layers' mean errors, over land and over sea, W/m²:

| State | PlaSim alone | The fit near 295 K | Its refit at 301.2 K |
|---|---|---|---|
| 295.0 K, the fit's own year: clear columns; all | 11.9–16.0; 12.7–17.8 | 0.2–0.4; 1.0–1.7 | 5.4–5.8; 4.6–6.0 |
| 299.6 K: clear columns; all | 17.5–18.9; 18.2–18.4 | 4.8–6.3; 5.5–5.9 | 1.8–2.7; 1.8–4.2 |
| 301.2 K: clear columns; all | 18.9–19.6; 19.0–19.6 | 5.9–6.8; 6.2–6.7 | 0.2–0.3 (fitted there); 0.6–1.4 |
| 302.5 K: clear columns; all | 18.8–20.1; 17.9–19.6 | 6.9–7.6; 6.4–6.6 | 3.3–3.6; 3.3–4.3 |
| 306.6 K: clear columns; all | 21.7–22.3; 20.5–21.4 | 13.7–14.9; 12.6–14.0 | 13.3–14.1; 12.4–13.9 |

- **The fit near 295 K does not hold at 301 K.** Its layers stray up to 5.9 W/m² over land and 6.8 over the sea at 301.2
  K, 2.4–2.7 times the 2.5 W/m² it left in any band of latitude on its own columns. The second layer from the ground
  cools 5.9–6.8 W/m² too much, the fourth and the ninth 3.4–4.4 too little and the seventh 4.0 too much. At 302.5 K the
  four layers that stray most stray 4.2–6.9 W/m² over land and 5.0–7.6 over the sea, the largest the lowest layer over
  the sea: this is the figure `research/plan.md`'s "5–7 W/m² at 302 K" stands for, which had no product behind it before
  this check.
- **PlaSim's own error grows much faster than the temperature's fourth power.** Between 295 and 301.2 K its second layer
  from the ground goes from 11.9 to 19.1 W/m² of extra cooling, its third from 6.7 to 9.3, and its outgoing longwave
  from 3.7 to 7.3 W/m² too high, while T⁴ rises 8%. The likely cause is the vapour: the Moon's moister air carries water
  paths beyond those Sasamori's broadband absorptivities were fitted to on Earth.
- **Refitted on the 301.2 K columns** (`plasim_longwave fit A28_corrected_970_lw_ice --year 2 --every 1`, 1,553
  cloud-free columns, [longwave_A28_corrected_970_lw_ice.json](../results/gcm/longwave_A28_corrected_970_lw_ice.json)),
  the correction leaves every layer within 0.3 W/m² on those columns, within 1.4 on all of that year's columns and
  within 2.1 in every band of latitude. Its extra escape grows at the top, 3.1 to 8.1 W/m² over land and 13.9 to 16.4
  over the sea, and in the second layer from the ground, 12.0 to 19.0 W/m², and changes sign in the middle air. It holds
  within 1.8–2.7 W/m² on the clear columns 1.6 K cooler and 3.3–3.6 W/m² 1.3 K warmer, and not at all at 306.6 K.
- **The correction therefore belongs to its state.** The runner's `--longwave-correction` takes the refit since 19
  October (`--longwave-fit 301k`; `295k` gives the fit the runs of 17 October took). A design case that settles more
  than about a kelvin from 301.2 K, or whose own columns differ from the test's, needs its correction checked on its own
  years, which the heating diagnostics make possible, and refitted if its layers stray more than 2.5 W/m². Correlated-k
  is the target throughout: its own error against the line-by-line model, up to about 5 W/m² in net flux at some level
  by the world-crafting check (item 2), is not located here.

#### The curved air's light on a 301 K column

The ledger of 14 October ([the curved air's light](#the-curved-airs-light)) was built on one clear column,
illumination's design column at 294.9 K (`illumination/surface_light/model.py`, `DESIGN_TS`), with 205 kg/m² of water.
`twilight.py --ts 301` builds the column the same way at 301 K: the moist adiabat from a 301 K ground with
Manabe–Wetherald humidity up to the tropopause, and above it the solved upper air and trace gases of the middle
atmosphere's two stored cases that bracket it, `moon_1.2atm_titania_stack_ts298` and `_ts308`, weighted 0.7 and 0.3, as
the design column weights its cases at 288 and 298 K. It is the atmosphere domain's own clear column at that ground,
self-consistent in temperature, vapour and tropopause. It holds 331 kg/m² of water against 231–241 kg/m² in the GCM's
3-day means as the 20-year test passes 301.2 K (its clear columns, then all), as the 294.9 K column held 205 against
158–170 in the design case: the radiative-convective humidity runs a fifth to a third wetter than PlaSim's. Its
tropopause stands at 11.5 kPa, 111.5 km up, against 15.1 kPa and 95.8 km.
[twilight_301k.json](../results/gcm/twilight_301k.json) (`terluna.climate.gcm-twilight/2`) records the column's source.
Global means, W/m², each ledger at its own column's temperatures:

| | 294.9 K column (`twilight.json`) | 301 K column (`twilight_301k.json`) |
|---|---|---|
| Sunlight the geometry adds: ground; air (cone + curved light) | +1.53; +8.79 (3.99 + 4.80) | +1.53; +10.31 (4.95 + 5.36) |
| Sunlight in all | +10.32 | +11.85 |
| Infrared the geometry lets out: cone + curved | 16.15 (10.52 + 5.63) | 19.85 (13.69 + 6.16) |
| Net | −5.83 | −8.00 |
| The same infrared scaled to the test's 3-day means at 301.2 K, as the GCM scales it | 17.61 | 17.97 |
| Water's slant factor at the horizon; well-mixed air's | 13.54; 7.14 | 12.71; 7.04 |
| Channels; solved by the sky's multiple scattering | 10,697; 282 | 11,559; 267 |

- **More of the sunlight is taken in the air**, 1.5 W/m² more, most of it in the third to fifth layers from the top,
  where the warmer, wetter, taller troposphere puts more vapour. The ground's twilight is unchanged.
- **The infrared escapes from higher, warmer layers.** The column's third layer from the top is 15 K warmer, so the
  shells' wider area lets out 3.7 W/m² more at its own temperatures. Scaled with the fourth power to the GCM's own
  layers at 301.2 K, whose upper air the longwave correction keeps colder, the two ledgers' escapes differ by 0.36 W/m²:
  the scaling carries most of the temperature, and the sunlight's 1.5 W/m² is the main change.
- **The GCM keeps the calibration's slant factors.** The shortwave calibration was refitted with the 294.9 K ledger's
  slant factors, so the runner takes the Moon's slant path from `twilight.json` and the terms from `twilight_301k.json`
  (`SLANT_PRODUCT` and `TWILIGHT_LEDGERS`); the 301 K column's own factors, 6% lower for water at the horizon, are
  recorded beside them.
- **The dynamics' thin shell.** The 301 K column holds 81,028 kg/m² of air over each square metre of ground with gravity
  falling with height and 85,876 kg/m² counting the shells' widening, against the GCM's p_s/g_s of 74,861 kg/m² (+8.2%
  and +14.7%; at 294.9 K +7.2% and +13.4%). PlaSim's and CM1's dynamics hold one gravity and a thin shell, so their heat
  capacity and vapour mass stay at p_s/g_s; the radiation's calibrations were fitted against the g(z) column, which
  absorbs part of the difference. No setting changes it.

**The terms in the clear share.** The twilight patch put each layer's sunlight into every column and only the ground's
into the column's clear share, and let out the infrared in every column, while the terms are clear-sky values. A new
patch, `TWILIGHT_CLEAR_EDITS` (`NTWLCLEAR`, the runner's default since 19 October; `--no-twilight-clear-share` gives the
earlier behaviour), adds each layer's sunlight in the share of the column clear of cloud above it, PlaSim's own random
overlap of cover, and lets each layer's and the ground's infrared through the clouds above it at PlaSim's own longwave
transmissivity of each cloudy layer, ice optics included. With it off the patched subroutines repeat the originals bit
for bit (a test compiles both on a stand-in module), and the patched `radmod.f90` compiles with PlaSim's own flags
against the deployed build's modules. A first size from the 20-year test's 3-day means at 301.2 K, whose cover in the
top three layers is 6%, 30% and 11% and leaves 55% of the sky clear above the ground: the clear share takes 3.7 W/m² of
the terms' sunlight out of the layers and keeps in the column 1.9 W/m² of their infrared if the high cloud passes three
quarters of it, as thin ice of 1.5–4.7 g/m² does (emissivity 0.24), or 6.2 W/m² if every cloud were opaque: −1.9 to +2.5
W/m² at the top of the air, about −0.9 to +1.2 K (`gcm-runs/spinup_2026-10-19/clear_share_sizing.py` on the research
drive). The test's own cloud is the guide; the design case's will be its own.

**The optics caches' keys.** The near infrared's and the longwave's line optics are cached on the drive under keys that
held the column's heights and the line settings alone, so new HITRAN files or a change to the line code would have kept
an old cache. Their keys now hold the column's levels (pressure, temperature, water and heights), the line settings, the
producing function's source and the hashes of the line model's code, the atmosphere products the column is built from
and the restored spectroscopy files, as the sky model's own optics record them (`solved_optics.inputs`, which refuses a
spectroscopy file whose checksum differs from its manifest). The sky batches' key adds the solver's own code
(`illumination/sky/solver.py`) to the transport's; the 294.9 K ledger's batches keep their names, so rebuilding that
ledger would solve its sky again.

**The scattering orders.** The sky model's own rule refuses a solution whose scattering orders stopped at the solver's
cap (`illumination/sky/solved_sky.py`, which sends it to `solved_checks.finish_orders`). Both ledgers' first batch
(366–513 nm) left one channel above the tolerance of 2 × 10⁻⁵ at the cap of 100 orders. `solve` now completes such
channels, and `ledger` refuses batches that are not converged. In the 301 K ledger the channel ends at order 106 and
`orders_converged` is true; the six orders add 8.8 × 10⁻⁶ W/m² to the ground and 3 × 10⁻⁷ W/m² to the air. Completed on
a copy, the 294.9 K ledger's open channel ends at order 104 and adds 5.8 × 10⁻⁶ and 2 × 10⁻⁷ W/m², so its
`orders_converged: false` stands for a leftover a million times smaller than its terms; `twilight.json` itself is left
as every run from 14 to 19 October read it (`gcm-runs/twilight/analysis/leftover_orders.py` and its two
`leftover_*.json`).

#### The surface stress and the 10 m wind in the output

The lake studies take their 10 m wind from a neutral log profile from the lowest model level, about 0.9 km up. To let
the next study bracket it with the model's own surface layer, the runner's output now holds the surface stress, `tauu`
and `tauv` (PlaSim's codes 180 and 181, means over each 3-day output interval). ExoPlaSim's postprocessor also names the
10 m wind, `uas` and `vas` (codes 165 and 166, `pyburn.py`), but PlaSim 3.4.2 never writes them: its output routine
(`outmod.f90`) has no 10 m wind, and the postprocessor then leaves a requested name out of the file without an error.
The author asked for a 10 m wind of PlaSim's own, and since 19 October `--wind10` gives one (`WIND10_PATCHES`, in the
runner's copy of ExoPlaSim only). PlaSim's surface layer takes it on Geleyn's (1988) profile through the lowest layer,
u(z) = u_L (ln(1 + z/z0) − (z/z_L)(B_N − B_D))/B_D, from the drag it uses itself: B_N = ln(1 + z_L/z0) = k/√C_DN,
and B_D = k/√(C_DN f_m) with f_m the stability factor of its momentum transfer. Neutral air gives the log profile,
stable air brings less of the wind down and unstable air more; over the sea at the lowest level's 885 m the neutral
factor is 0.72, and f_m of 0.5 and 1.5 give 0.51 and 0.88. The output writes its mean over the steps since the last
output as codes 165 and 166, which the runner then asks for as `uas` and `vas`; the mean counts its own steps, since
PlaSim reads its restart before this namelist. Without the switch the model is as it was. A test compiles the surface
layer's lines and checks them against the formula, and another checks every output name against the codes PlaSim's
output routine writes.
The stress also gives the friction velocity, u* = (|τ|/ρ)^½, as a check on the profile's.

#### The ice optics' weights

The ice optics take PlaSim's near-infrared band as four of Ebert and Curry's bands weighted by their shares of the
unfiltered TSIS-1 Sun, 0.6164, 0.2283, 0.0811 and 0.0742 (`plasim_ice_optics.NIR_WEIGHTS`). Through the titania stack
the same rule gives 0.6258, 0.2262, 0.0748 and 0.0733 (checked on 19 October), which lowers one ice layer's absorption
in that band by 2–3%, a few hundredths of a W/m² over the Moon, and only with `--ice-optics`. The world-crafting lab's
patch 11 (github.com/dgoldman0/world-crafting, branch `lab`, commit 79ea224,
`models/exoplasim/patches/11-ice-weights.patch`, sha256 `ad31c4a23479e683`) reads the weights from the namelist as
`ZICEW`, its default the constants, so that a runner can write the shield's. Its patches 01–10 give this copy's
`radmod.f90` and `rainmod.f90` byte for byte from ExoPlaSim 3.4.2's own sources (`a49791572fa8f5b2` and
`44c0668de3e8f7db`; the sources recovered both by undoing Terluna's edits and from an unpatched install agree), and
patch 11 applies after them, giving `radmod.f90` `50e3df63c3d9872e`. At the author's word it is ported unchanged in
substance as `ICE_WEIGHTS_EDITS`, after the clear share's patch, in the runner's copy of ExoPlaSim only; runs with ice
optics now write the shield's weights as `ZICEW` (`plasim_ice_optics.SHIELD_NIR_WEIGHTS`, tested against
`shielded_nir_weights()`), and `--ice-weights stock` writes nothing and keeps the earlier runs' constants.

#### The design years (20 October)

The cold start (`A28R_301k`, sunlight 0.925, on the rebuilt lakes) ran 25 years heading near 305 K, about 4 K warmer
than the target study's runs, so a branch from its year-24 restart (`A28R_301k_s0894`) ran 15 years at 0.894. A twin of
that branch without the clear share (`A28R_301k_s0894_noclear`) measured the clear share: it cools the climate by
about 1 K (their year 14 at 301.05 and 301.93 K, settling near 301.0 and 301.9 K by `compare --settle`). The
design years (`A28R_301k_s0897`) ran on from the branch's year-14 restart at 0.897, 15 years on 4 ranks
(`run_A28R_301k_design.sh` on the research drive).

- **They settled near 301.3 K.** Years 5–14 average 301.25 K at the ground, every year within 0.13 K of the others,
  with the top of the air in balance (+0.02 W/m²) and a drift of +0.015 K a year; years 10–14 average 301.29 K. On
  20 October the author took years 5–14 as the design years, at sunlight 0.897, 7.5% dimmer than the 0.970 of the
  295 K design case.
- **`compare --settle` read them low until fixed.** It ended a run where its imbalance reached the leftover
  +0.62 W/m² of the 20-year test (`A28_corrected_970_lw_ice`, years 10–19), which these runs do not keep, so it read
  299.4 K for the design years and 300.3 K for the branch at 0.894. It now takes each run's own leftover from its
  drift: each year-to-year warming against the two years' mean imbalance, N = N_left + C dT/dt, at no warming, with
  its standard error ([compare.py](compare.py)). The design years read 301.2–301.3 K, the branch at 0.894 301.0 K and
  its twin 301.9–302.0 K, each within 0.1 K. So the branch at 0.894 was on target already, and the design years at
  0.897 sit 0.2–0.3 K above it.
- **The longwave refit holds on them, at its edge**
  ([longwave_check_A28R_301k_s0897.json](../results/gcm/longwave_check_A28R_301k_s0897.json), `plasim_longwave check`
  with `--case`). On years 5, 9 and 14 its largest layer mean over land or sea is 1.9–2.65 W/m², against 6.2–6.9 for
  the 295 K fit and about 20 for PlaSim alone. One sample of six, year 14's cloud-free columns, passes the 2.5 W/m²
  the rule allows, by 0.15. A refit would need a new design run, so it waits for the next change to the GCM's physics.
- **The high cloud's water and crystals are worth several kelvin** ([high_cloud.py](high_cloud.py),
  [high_cloud_A28R_301k_s0897.json](../results/gcm/high_cloud_A28R_301k_s0897.json)). The top three layers over the
  design years, with the ice optics the GCM runs with:

  | Layer | Cover | In-cloud water | Temperature | Crystal radius | Reflection; emissivity |
  |---|---|---|---|---|---|
  | Above 107 km | 9.5% | 1.16 g/m² | 162 K | 2.2 µm | 19.0%; 0.58 |
  | 81–107 km | 33.6% | 2.03 g/m² | 192 K | 11.5 µm | 8.3%; 0.27 |
  | 63–81 km | 11.5% | 6.31 g/m² | 222 K | 35.2 µm | 8.1%; 0.30 |

  Scaling every layer's water or crystal radius, first order (the cloud's amount, the vapour and the circulation
  held, the layers counted apart, the sunlight averaged over day and night), at the top of the air and at the ground
  through the branches' feedback of 1.20–1.56 W/m² a kelvin:

  | Change | At the top of the air | At the ground |
  |---|---|---|
  | Water halved | −8.6 W/m² | −5.5 to −7.1 K |
  | Water 1.5 times | +6.5 W/m² | +4.2 to +5.4 K |
  | Water 2.5 times | +15.2 W/m² | +9.8 to +12.7 K |
  | Water 3.5 times | +20.3 W/m² | +13.0 to +16.9 K |
  | Crystals half the size | +10.6 W/m² | +6.8 to +8.8 K |
  | Crystals twice the size | −7.8 W/m² | −5.0 to −6.5 K |

  Most of each change is the infrared: more water or smaller crystals trap more of it than they reflect sunlight. A
  factor of two in either moves the design climate by 5–7 K, more than any other uncertainty measured in it so far,
  and the ice CM1 made of its own at these heights, 1.5–3.5 times PlaSim's water, would be worth +4 to +17 K. The
  coldest layer's crystals rest on the table's trend below 180 K. It is a scale for the high cloud, not a run: a GCM
  answer needs the cloud's water and crystals from a lunar treatment of the cold cloud.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave check --case A28R_301k_s0897:5 --case A28R_301k_s0897:9 \
    --case A28R_301k_s0897:14 --output climate/results/gcm/longwave_check_A28R_301k_s0897.json      # minutes, cached
climate/gcm/.venv/bin/python -m climate.gcm.high_cloud A28R_301k_s0897 --years 5-14                 # seconds
```

#### The cold cloud at lunar gravity (20 October)

On 20 October the author asked for the high cloud's water and crystals at lunar gravity, as factors on the GCM's Earth
values for cloud colder than −38 °C ([register](../../research/decisions.md), 20 October). [cold_cloud.py](cold_cloud.py)
runs slabs of rising air on both worlds through the design years' top three layers, with homogeneous freezing,
kinetically limited growth and fall-out ([the snow's account](../../atmosphere/snow/README.md#the-high-clouds-crystals)),
at the same pressure velocity: the one at which Earth's slab holds PlaSim's own water path there, 0.053–0.058 Pa/s.
Where the design years condense vapour in those layers (a quarter of the 3-day means at 81–107 km, 7% at 63–81 km),
their large-scale condensation against the humidity of the layer below gives an ascent of 0.03–0.15 Pa/s between
quartiles, 0.07 at the median, about the same
([cold_cloud_A28R_301k_s0897.json](../results/gcm/cold_cloud_A28R_301k_s0897.json)):

| Layer | PlaSim's water; radius | Factors: water; radius | Their range |
|---|---|---|---|
| Above 107 km | 1.16 g/m²; 2.2 µm | 6.1; 1.00 | 6.1–6.5; 1.00 |
| 81–107 km | 2.03 g/m²; 11.5 µm | 36.6; 1.68 | 12.8–52.8; 0.98–3.05 |
| 63–81 km | 6.31 g/m²; 35.2 µm | 42.1; 1.60 | 36.2–55.3; 1.44–1.79 |

The range takes every slab of 2–8 kPa at ascents within a factor of three of the central one, and the deposition
coefficient and the haze's number at it. At first order (high_cloud.effect, the run's own cloud and cover held) the
factors make each top layer nearly black in the infrared, its emissivity 0.995–1.000 against 0.27–0.58, and let it
reflect 53–57% of the sunlight a fully cloudy layer receives against 8–19%: +17.8 W/m² at the top of the air, nearly
all from the 81–107 km layer (+15.2), whose third of the sky now radiates to space at 192 K; at 63–81 km the added
reflection and trapping cancel (−0.3). Through the 301 K branches' feedback that is +11.4 to +14.8 K at the ground.
Across the corners of the ranges it runs from +9.0 W/m² (the most water on the smallest crystals) to +22.0 (the least
on the largest), +5.8 to +18.3 K: the infrared saturates at nearly every corner, the reflection does not, and more
water on smaller crystals sends back more sunlight.

The runner's `--cold-cloud PRODUCT` writes the central factors into the namelist level by level (ZICEWF and ZICERF,
`ICE_FACTOR_EDITS`), which scale the ice optics' water path and radius in sunlight and in the infrared; the run's
configuration records them with the product's hash. PlaSim's cover stays its own: it draws high cloud where the air
passes 136–149% of ice saturation, while a lunar ice cloud, whose crystals stay, would outlast the supersaturation
that formed it, so a run gives the thicker cloud on PlaSim's cover and no more cover. In sunlight the patch keeps
PlaSim's cap of 1,000 g/m² on a layer's water path, above the means (258 g/m² at 63–81 km).

At their defaults of 1 the patched model repeats the design run's year 14 from its year-13 restart in every one of its
56 output variables (`A28R_301k_cold_regression`). The run `A28R_301k_cold` carried the factors from the design years'
end for 20 years at the same sunlight, 0.897 (`run_A28R_301k_cold.sh` on the research drive), and settles near
298.5 K (±0.1 K by `compare --settle` over its years 5–19), 2.8 K below the design years. Its years 10–19 against them:

| | Design years (5–14) | With the lunar cold cloud (10–19) |
|---|---|---|
| Ground; sea surface | 301.25; 302.6 K | 298.66; 300.1 K |
| Planetary albedo | 0.298 | 0.370 |
| Cloud cover: all; 81–107 km | 49.6%; 33.6% | 34.8%; 14.5% |
| Cloud effect: sunlight; infrared; net | −18.2; +29.2; +11.1 W/m² | −37.5; +48.6; +11.1 W/m² |
| Sunlight absorbed at the ground | 107.7 W/m² | 93.8 W/m² |
| Evaporation's heat; precipitation; water vapour | 78.1 W/m²; 2.70 mm/day; 248 kg/m² | 63.5 W/m²; 2.20 mm/day; 223 kg/m² |
| Warmest; coldest 3-day mean anywhere | 309.3; 292.3 K | 306.5; 289.8 K |

- **The cover retreats where the ice thickens.** From the first year the thick ice warms its own layer at 81–107 km,
  absorbing more sunlight in the near infrared (0.08 against 0.04 K a day) and more of the infrared from below, and
  the layer's relative humidity falls; PlaSim's cover there, drawn from the humidity, falls from a third of the sky to
  a seventh and stays there.
- **The ground loses sunlight.** The cover that remains reflects about half the sunlight it receives: the planetary
  albedo rises from 0.30 to 0.37, and the ground absorbs 14 W/m² less. The infrared the thicker cloud traps warms the
  air above it, the top two levels by 2–4 K, while every level below cools by 2–3 K.
- **The cloud's net effect at the top of the air stays +11.1 W/m²**, its sunlight and infrared terms each larger by
  about 19 W/m²: the change moves sunlight from the ground to the upper air and to space.
- **The water cycle slows.** The ground evaporates 15 W/m² less: 19% less rain and 10% less vapour.

With the cover held, the factors gave +17.8 W/m² at first order and a warmer Moon; the run's own cover turns that
into a cooler one. The cover is PlaSim's diagnosis from humidity over liquid water, which a lunar ice cloud whose
crystals stay would outlast, so its retreat is the least certain part of the run. Holding 301 K with this cloud would
take about 2% more sunlight (1.9% by the run's own feedback of 1.22 W/m² a kelvin on its 180 W/m² of absorbed
sunlight, 2.3% at the target study's 0.845% a kelvin), a sunlight scale near 0.915–0.918 in place of 0.897; the
design years stand until the author chooses.

The longwave refit of 19 October holds on the run while it stays near 301 K and not at its end
([longwave_check_A28R_301k_cold.json](../results/gcm/longwave_check_A28R_301k_cold.json), `plasim_longwave check`,
which now also judges every case it is given): its largest layer mean over land or sea is 1.9 W/m² in year 5
(299.4 K) and 3.1–3.6 W/m² in years 12 and 19 (298.8 and 298.5 K), against the rule's 2.5 and PlaSim's own 18–19,
the excess in the lowest layer over the sea, which it cools too much. A design case that takes the cold cloud's
factors takes its own refit with them.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.cold_cloud A28R_301k_s0897 --years 5-14                 # about 16 minutes
climate/gcm/.venv/bin/python -m climate.gcm.exoplasim_run A28R --folder A28R_301k_cold --years 20 --ncpus 4 \
    --start-from climate/gcm/runs/A28R_301k_s0897/model/MOST_REST.00014 --set sunlight_scale=0.897 \
    --longwave-correction --ice-optics --heating-diagnostics --wind10 \
    --cold-cloud climate/results/gcm/cold_cloud_A28R_301k_s0897.json
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave check --case A28R_301k_cold:5 --case A28R_301k_cold:12 \
    --case A28R_301k_cold:19 --output climate/results/gcm/longwave_check_A28R_301k_cold.json    # about 2 minutes
```

#### The evaporation of falling rain (21 October)

PlaSim rains out what its large-scale and convective schemes condense at once, and each unsaturated layer the rain
falls through then evaporates a share gamma of its saturation deficit each time step (rainmod's mkrain), corrected for
the latent cooling and at most the rain there: one constant, 0.01 at T21, set for Earth, which knows nothing of the
drops' size or fall. [rain_evaporation.py](rain_evaporation.py) measures a run's evaporated share two ways from its
output, the rain the schemes form (their heating, written before mkrain) against the rain at the ground, and mkrain's
cooling left in the summed physics tendency; the two agree within 1.1 points
([rain_evaporation.json](../results/gcm/rain_evaporation.json)). The heating diagnostics are snapshots and the rain a
mean, so only long means compare.

The cold cloud's run evaporates 29% of its falling rain, the same over land and sea, and 32% on the equatorial land
about 0° E where CM1's small box stands. On Earth that would be an ordinary share. CM1's box, with Morrison's rain
limits carried to lunar gravity, evaporates 48% of the water leaving its cloud, and 86% with CM1's own limits: lunar
drops of a size lose far more water falling slowly through six times the air, and the lunar limits let them grow
large enough to bring the loss back to 48% ([the CRM's account](../crm/README.md)). Two years at each of two
constants from the cold cloud's run's year-19 restart (`run_rain_gamma.sh` on the research drive):

| gamma | Evaporated: global; at the box | Rain formed; at the ground | Cloud; albedo | Ground after two years |
|---|---|---|---|---|
| 0.01, PlaSim's own (years 18–19) | 28.6%; 31.8% | 3.05; 2.17 mm/day | 34.8%; 0.370 | 298.50 K |
| 0.02 | 44.5%; 41.0% | 3.73; 2.07 mm/day | 40.1%; 0.394 | 298.11 K, top −0.3 W/m² |
| 0.04 | 60.1%; 59.4% | 4.50; 1.80 mm/day | 53.2%; 0.437 | 296.73 K, top −4.9 W/m² |

- **CM1's 48% at the box comes at gamma 0.026**, 2.6 times PlaSim's own, interpolated in log gamma.
- **More evaporation moistens the lower air, and PlaSim's low cloud grows.** The cloud's effect on the sunlight
  strengthens by 7 W/m² at 0.02 and 20 W/m² at 0.04, and the Moon cools; at 0.026 a settled run would likely sit
  1–2 K below the cold cloud's 298.5 K. Less of the rain formed reaches the ground.
- **The match rests on one box.** CM1's small box over one lunar day, which rains a third as much as the GCM at its
  place, is the only lunar measure of the share.

**More years, and the Near-side Sea** (21 October). At the author's word the constant is tested before the design run
takes it, by a pair of CM1 boxes on the design climate over the Near-side Sea at 5° S, 45° W and over the land at 0° E
([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)). The measure now pools
the runs at each constant (`run_column_budget.sh` on the research drive carried the two branches to five years each,
and the cold-cloud run's branch `A28R_301k_cold_moist` adds six years at 0.01), and sets each box beside the GCM at its
own place, the equatorial land about 0° E or the sea at sea level within 11° of latitude and 20° of longitude of the
sea box's site:

| gamma | Years | Evaporated: global | Land about 0° E | The Near-side Sea |
|---|---|---|---|---|
| 0.01, PlaSim's own | 14 | 28.6% | 29.8% | 34.5% |
| 0.02 | 5 | 44.3% | 43.9% | 49.7% |
| 0.04 | 5 | 60.2% | 58.2% | 75.8% |

Over the Near-side Sea PlaSim evaporates more of its rain than at the land box, and the gap widens with gamma: its
rain there is lighter, and at 0.04 a quarter of it reaches the ground. With the pooled years the small box's 48% comes
at gamma 0.024 at its place; at that constant the Near-side Sea would evaporate about 57%. At 0.02 the ground settles
near 298.0 K over the branch's five years, and at 0.04 it is still falling, 295.5 K at year 4.

**Where the GCM forms its rain, and why the boxes' shares do not yet test gamma** (21 October). At the Near-side Sea
PlaSim's convection scheme forms 3.5 of the 4.0 mm/day of rain formed. It forms 3.2 of it between 18 and 72 km, and
the rain falls that far through the air below. The first sea box held its air to the GCM's above 8–16 km and took the
GCM's forcing to 16 km only. Its nudging supplied the upper air's heat in place of convection, its clouds stayed below
10 km, and it rained 0.06 mm/day against the GCM's 2.72. Its share, 53% of the water leaving its cloud, comes from
rain that falls a fraction of the GCM's way ([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)).
The small box's 48% carries the same doubt: its nudging supplied 50 of the 90 W/m² its radiation removed. CM1's
boundary layer also put its density times the surface's heat and water into the air. The pair runs again with the
forcing through the column, the air held only above the GCM's convection, and that fixed (domain/climate 13e7123).
- **The sea box** now convects to about 70 km. It rains 31% of the GCM's rain and evaporates 61% of its falling water,
  gamma about 0.027 by the pooled shares. That figure is provisional: the box's upper air was still cooling, and its
  rain falls at other hours than PlaSim's (below).
- **The land box** rains 92% of the GCM's rain at the same hours of the day, and also evaporates 61% of its falling
  water. Over that land PlaSim evaporates only 58.2% even at 0.04, so the land wants about 0.046, beyond the runs.

No single constant fits both places. At one gamma PlaSim evaporates far more of its rain over the sea than on land,
75.8% against 58.2% at 0.04, while CM1 evaporates about the same at both. The design run's constant waits on that,
and on PlaSim's convection (below).

The runner's `--rain-evaporation GAMMA` writes the constant into PlaSim's rain namelist, recorded in the run's
configuration only when given.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.rain_evaporation A28R_301k_cold:12-19 A28R_301k_cold_moist:0-5 \
    A28R_301k_cold_g02:0-4 A28R_301k_cold_g04:0-4                                          # a minute
```

#### PlaSim's convection at the Moon's long day (21 October)

**What the sea box showed.** With the forcing through the column, CM1's sea box convects as deep as PlaSim's
convection forms its rain, up to about 70 km
([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)). Over its second lunar
day it rained 31% of the GCM's rain, but unevenly through the day:
- **Noon to midnight.** It rained 34–58% of the GCM's rain.
- **Midnight to noon.** It rained almost nothing. At those hours the GCM's own circulation sinks over the site, warming
  the upper air by about 95 W/m² and drying the column by 2–4.6 mm/day, while PlaSim rains 2.9 mm/day after midnight
  and 0.9 mm/day in the morning.

**Why PlaSim rains through the night.** Its deep convection is a Kuo scheme (`rainmod.f90`, `nprc` 1). Wherever a
parcel lifted undiluted from the surface is buoyant, it turns the column's moisture supply into rain, the surface's
evaporation included, so the sinking air's warming and drying do not stop it. After midnight the site's sea evaporates
3.4 mm/day and the scheme rains 2.9 of it.

**Why PlaSim still carries it.** The leading forecast and climate models dropped Kuo schemes decades ago; ECMWF moved
to Tiedtke's mass-flux scheme in 1989. The closure has been criticized for reversing the causality between convection
and the convergence it largely makes (Emanuel, Neelin and Bretherton 1994, *QJRMS* 120; Arakawa 2004, *J. Climate* 17).
PlaSim, built for speed and portability as a model of intermediate complexity, kept the Kuo scheme with Tiedtke's
1983 shallow convection. ExoPlaSim kept PlaSim's physics (Paradise and others 2022, *MNRAS* 511), and its namelist only
switches the scheme on or off.

**The replacement** (the author's direction of 21 October; the climate stays on ExoPlaSim). A simplified Betts–Miller
scheme (Frierson 2007, *J. Atmos. Sci.* 64) is written in Fortran as a Terluna patch to `rainmod.f90`, behind PlaSim's
own switch: `nprc` 1 keeps Kuo exactly, and 2 selects the new scheme.
- **How it works.** It relaxes an unstable column's temperature toward a moist adiabat lifted from the lowest level, and
  its vapour toward a fixed relative humidity of that profile, over a relaxation time. It rains the vapour removed, with
  Frierson's correction keeping the column's moist energy. So it rains only while the column is both unstable and
  moister than its reference, and the night's sinking air shuts it off as it shuts off CM1's convection.
- **Its parameters.** Earth's values are 2 hours and 70%; both are to be fitted to CM1.
- **What it costs.** Per column it does what Kuo does, one parcel lifted through the 10 levels with PlaSim's own
  saturation formula and a few column sums, so it should run at about Kuo's speed.
- **The grid.** Each of PlaSim's 10 layers is 5–30 km thick on the Moon, too coarse for a plume scheme's mixing to
  show, so a relaxation scheme fits the grid.
- **Two lunar options, off by default.** One is a reference lifted with some mixing of surrounding air, since the sea
  box's air at 40–60 km drifted 8–11 K below the GCM's, which Kuo heats toward an undiluted parcel. The other is an
  inhibition limit: the box did not rain in the mornings at 55 J/kg.
- **Shallow convection.** PlaSim's shallow scheme stays for columns that are unstable but too dry to rain.
- **What it leaves out.** It has no plumes, downdrafts, cold pools or organization, and its rain follows the
  instability after a fixed delay.

**Its tests.** The scheme is written (`SBM_CODE` and `SBM_EDITS` in [exoplasim_run.py](exoplasim_run.py), the
last of `PATCHES`; `--convection sbm` with `--sbm-tau`, `--sbm-rh`, `--sbm-eps` and `--sbm-cin` selects it), with a
Python copy ([sbm.py](sbm.py)). The first test is done (21 October,
[test_sbm.py](../tests/test_sbm.py)): the Fortran routine, compiled alone against a stub of PlaSim's module, gives
the Python copy's tendencies to 10⁻¹⁰ on five columns, and the copy keeps the column's moist enthalpy, rains the
vapour it removes in proportion to the relaxation rate, leaves dry and stable columns alone, sends a buoyant column
too dry to rain to the shallow scheme, and with mixing cools the parcel through its convecting layers. The patch
applies once to ExoPlaSim's own `rainmod.f90` with Kuo's path kept. On the Near-side Sea's mean column from the
GCM the scheme rains only with a reference humidity of 60% or less, so the reference humidity is the first thing
the branches have to fit.

**In the GCM** (21 October; `run_after_pair.sh` on the research drive). One year through the patched model with Kuo's
scheme from the design climate's year-19 restart, `A28R_301k_cold_kuo_check`, matches `A28R_301k_cold_moist`'s year 0
bit for bit, in 237 s against 233 s. One year of the new scheme at a reference humidity of 0.6,
`A28R_301k_cold_sbm_smoke`, ran without fault in 236.5 s, Kuo's speed. Against Kuo's year it ends:

| After one year | New scheme, rh 0.6 | Kuo |
|---|---|---|
| Ground | 298.66 K | 298.41 K |
| Air at 2 m | 299.03 K | 298.22 K |
| Rain | 2.39 mm/day | 2.19 mm/day |
| Cloud cover | 36% | 35% |
| Longwave cloud warming | 52.0 W/m² | 48.7 W/m² |
| Imbalance at the top of the air | +3.76 W/m² | +0.07 W/m² |

So it is still warming after its first year, and its settled climate needs longer runs.

**The new scheme settled** (the author's go-ahead of 22 October; `run_sbm_branches.sh` and `sbm_settled_analysis.py`
on the research drive). Two 15-year branches from the same restart, at the design climate's sunlight (0.897) with
every design flag, each at about Kuo's speed (240–306 s a year on a shared machine): `A28R_301k_cold_sbm_rh07` at
Frierson's settings (a reference humidity of 0.7, 2 hours) and `A28R_301k_cold_sbm_rh06` at 0.6. Both reach balance at
the top of the air by years 8–10. `compare --settle` puts them at 299.6 and 300.1 K, each within 0.1 K, against Kuo's
298.4 K, with leftover imbalances of −0.03 ± 0.08 and +0.12 ± 0.09 W/m² against Kuo's −0.18 ± 0.12. Over their last
five years against Kuo's settled years:

| Settled | Kuo, `A28R_301k_cold` 12–19 | New scheme, rh 0.7 (10–14) | New scheme, rh 0.6 (10–14) |
|---|---|---|---|
| Ground | 298.61 K | 299.56 K | 300.13 K |
| Sea surface | 300.1 K | 300.5 K | 300.9 K |
| Air over land | 298.2 K | 299.7 K | 300.5 K |
| Air poleward of 70° | 297.6 K | 299.6 K | 300.3 K |
| Equatorial land by day; by night | 299.1; 297.5 K | 300.2; 298.9 K | 301.1; 299.6 K |
| Coldest; warmest 3-day mean anywhere | 289.8; 306.5 K | 290.9; 308.0 K | 291.5; 308.7 K |
| Rain | 2.19 mm/day | 2.32 mm/day | 2.45 mm/day |
| Water vapour | 222 kg/m² | 245 kg/m² | 236 kg/m² |
| Cloud cover; planetary albedo | 34.8%; 0.370 | 35.2%; 0.378 | 34.9%; 0.370 |
| Cloud effect: sunlight; infrared | −37.5; 48.6 W/m² | −40.0; 50.1 W/m² | −37.6; 51.0 W/m² |
| Land under 0.5 mm/day of rain | 28.0% | 26.8% | 26.8% |

At the same sunlight the new scheme warms the Moon by 1.0–1.5 K at the ground, 1.5–2.3 K in the air over land and
2.0–2.7 K poleward of 70°. The column holds 6–10% more water and rains 6–12% more. Rain at the two CM1 box sites by
local time (mm/day, the GCM's years as in the table; the CM1 boxes' second lunar day, forced by Kuo's GCM):

| Near-side Sea | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise | Mean |
|---|---|---|---|---|---|
| Kuo | 0.92 | 3.50 | 3.54 | 2.93 | 2.72 |
| New scheme, rh 0.7 | 0.31 | 2.88 | 3.33 | 1.45 | 1.99 |
| New scheme, rh 0.6 | 0.74 | 2.93 | 3.38 | 2.68 | 2.43 |
| CM1 box | 0.00 | 1.19 | 2.04 | 0.17 | 0.85 |

| Land at 0° E | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise | Mean |
|---|---|---|---|---|---|
| Kuo | 1.91 | 11.06 | 0.96 | 0.21 | 3.57 |
| New scheme, rh 0.7 | 1.05 | 15.10 | 0.47 | 0.08 | 4.20 |
| New scheme, rh 0.6 | 3.27 | 12.67 | 0.47 | 0.11 | 4.15 |
| CM1 box | 2.87 | 8.59 | 1.53 | 0.07 | 3.28 |

- **Over the sea at rh 0.7** the rain after midnight halves and the morning's falls by two-thirds, toward CM1's
  pattern. The nights still rain about eight times CM1's. At rh 0.6 the nights keep nearly Kuo's rain.
- **Over land** both settings rain more and sharpen the afternoon peak, 15.1 and 12.7 mm/day against Kuo's 11.1 and
  CM1's 8.6, while the evening's rain halves, away from CM1's spread through the evening.
- **Against CM1 near the ground.** The corrected CM1 boxes ran 0.23 K warmer than Kuo's GCM at 2 m over the sea and
  1.15 K over land ([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)), so
  the new scheme's warming has the right sign and runs larger. Those boxes were forced by Kuo's GCM, and CM1 boxes
  forced by the new GCM's own budget are the comparison that tests it.

**The checks.**
- **The longwave.** The 301 K refit, on each branch's last year, leaves a largest layer error of 3.38 W/m², against
  the rule's 2.5 ([longwave_check_sbm.json](../results/gcm/longwave_check_sbm.json); `plasim_longwave check
  --case`). The new scheme's climate therefore takes a refit of its own, as the cold cloud's design case already does.
- **The exobase.** It stays as the 14 October check left it. The middle atmosphere's exobase moves by under 0.4 K
  across grounds of 268–308 K, and the water at PlaSim's top level falls from Kuo's 0.034 to 0.020–0.022 ppmv, with the
  air there 1.2–1.6 K cooler.
- **The energy.** The leftover imbalances above are no larger than Kuo's.

The new scheme's climate lands 1.0–1.5 K nearer the design climate's 301 K than Kuo's at the same sunlight, so a
design case with it would need less of the 2% more sunlight the cold cloud's factors called for. Its adoption for the
design climate, and its reference humidity, are the author's call; the CM1 boxes forced by it come next, at the
author's word. It would change the convective rain maps the electricity study reads ([convection.py](convection.py)).

**The scheme on the boxes' air and on its own** (24 October; [sbm_offline.py](sbm_offline.py),
[sbm_offline.json](../results/gcm/sbm_offline.json), 15 seconds). With the new scheme at rh 0.7 the GCM still rains
1.45 mm/day over the Near-side Sea after midnight, where the CM1 sea box rains 0.17. Run on through a third lunar day
as its upper air settles, the box rains 0.14 mm/day after midnight and 2.23 and 2.97 mm/day in the afternoon and
evening ([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)). The boxes
were forced by Kuo's GCM, and the new GCM's nights are harsher: after midnight the rh 0.7 branch's circulation warms the column above 16 km
by 131 W/m², against Kuo's 57, and dries it by 5.2 mm/day, against 3.3 (`sbm_forcing_compare.py` on the research
drive). The scheme's parcel rises undiluted, and cloud-resolving models respond far more to dry air aloft than schemes
with unmixed parcels (Derbyshire and others 2004, *QJRMS* 130). So before a branch with the parcel's mixing
(`eps_bm`), the scheme ran offline at the two box sites on two sets of columns, at its 2-hour relaxation time:
- **The settings.** A reference humidity of 0.6, 0.7 or 0.8 with mixing of 0, 0.5, 1, 2, 4, 8 or 16 per e-folding of
  pressure, and at rh 0.7 inhibition limits of 10 and 30 J/kg.
- **The constants.** The runs' own: the gas constant 288.38 J/kg/K from the planet namelist and ExoPlaSim's kappa of
  0.286.
- **The rain.** The scheme's rate at once, the rain formed before mkrain evaporates any of it.
- **CM1's columns.** Each box's domain mean every 3 hours over its second lunar day, laid on PlaSim's ten levels at
  the box's surface pressure. The 192-km domain stands for a GCM cell, the boxes' air was shaped by Kuo's forcing,
  and the scheme's instantaneous rate is set against CM1's realised rain.
- **The GCM's own columns.** The 3-day means in each site cell, six over the sea and four on land, over years 10–14
  of each settled branch, with the scheme at that branch's reference humidity.

**On CM1's air the scheme rains against CM1's hours.** In mm/day, with the share of the columns on which the scheme
convects and rains:

| Near-side Sea, CM1's columns | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise |
|---|---|---|---|---|
| CM1's own rain | 0.00 | 1.19 | 2.04 | 0.17 |
| rh 0.7, no mixing | 20 (39%) | 34 (81%) | 31 (47%) | 119 (100%) |
| rh 0.7, eps 1 | 138 (100%) | 154 (100%) | 226 (100%) | 287 (100%) |
| rh 0.7, eps 4 | 273 (100%) | 284 (100%) | 479 (100%) | 484 (100%) |
| rh 0.7, eps 8 | 0 (0%) | 148 (46%) | 140 (27%) | 537 (98%) |
| rh 0.7, eps 16 | 0 | 0 | 0 | 0 |
| rh 0.6, no mixing | 446 (100%) | 493 (100%) | 478 (100%) | 571 (100%) |
| rh 0.8, no mixing | 0 | 0 | 0 | 0 |
| rh 0.7, no mixing, inhibition limit 10 J/kg | 0 (0%) | 0 (0%) | 30 (41%) | 57 (52%) |

| Land at 0° E, CM1's columns | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise |
|---|---|---|---|---|
| CM1's own rain | 2.91 | 8.59 | 1.53 | 0.08 |
| rh 0.7, no mixing | 294 (56%) | 0 (0%) | 305 (76%) | 449 (100%) |
| rh 0.7, eps 1 | 379 (100%) | 97 (100%) | 505 (100%) | 600 (100%) |
| rh 0.7, eps 8 | 489 (75%) | 196 (27%) | 582 (73%) | 12 (2%) |
| rh 0.7, eps 16 | 0 | 0 | 151 (24%) | 0 |
| rh 0.8, no mixing | 101 (37%) | 0 (0%) | 7 (34%) | 7 (22%) |

- **The rain at once** is the vapour above the scheme's reference over the relaxation time. At rh 0.7 without mixing
  the vapour CM1's sea columns hold above 70% of the parcel's saturation, averaged with nil where they hold less, comes
  to 2.6–2.8 kg/m² in the afternoon and evening and 9.9 kg/m² after midnight, when CM1's own convection has stopped
  spending what the sea evaporates. Then the scheme's parcel from the lowest level is buoyant up to σ 0.32, about
  55 km up, against an inhibition of 14 J/kg, and CM1's own columns hold 1,420 J/kg of CAPE against 3 J/kg of
  inhibition.
- **Mixing** cools the reference, and the cooler reference holds less vapour, so the rain rises with the mixing until
  the diluted parcel loses its buoyancy. At eps 8 that happens in the sea's mornings and in half its afternoons, while
  98% of its post-midnight columns still convect; at 16 the scheme stops at every hour at both sites.
- **The other options** act in the wrong order. A reference humidity of 0.8 holds all the sea's vapour, and the scheme
  never rains there. An inhibition limit holds back the mornings and afternoons first, where the parcel meets
  70–88 J/kg, against 10–14 J/kg after sunset.
- **The mapping and the constants.** On layer means by mass, or on CM1's own 111 levels, the numbers move, but the
  scheme still convects after 78–100% of the sea's midnights at every mixing up to 16, at 46–577 mm/day. The
  constants of [test_sbm.py](../tests/test_sbm.py) in place of the runs' change the rain by up to 11%.

No setting passes the test, which asks for CM1's rain within half a millimetre a day or half its value at the sea in
every quarter and over land in the afternoon and evening. The settings that bring the sea's nights and mornings near
CM1's do it by stopping the scheme at the sea at every hour (rh 0.8 with mixing up to 2, or a mixing of 16); every
other rains at least 24 mm/day there after midnight, 140 times CM1's rain. Over land the scheme without mixing
convects on none of CM1's afternoons and on all its nights, yet in the GCM it rains at CM1's hours, 15.1 mm/day in the
afternoon and 0.08 after midnight. A model running the scheme keeps its column near the reference, and its rain
follows what the circulation, the sea and the Sun supply, so the scheme's rate on a cloud-resolving model's air does
not show when it would rain in the GCM.

**On its own columns, mixing moves the sea's rain out of the night** (mm/day; the rh 0.7 branch, the scheme at rh 0.7):

| Near-side Sea, the GCM's columns | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise |
|---|---|---|---|---|
| The branch's convective rain at the ground | 0.31 | 2.80 | 3.30 | 1.44 |
| The rain its scheme formed (snapshots) | 0.49 | 2.74 | 5.49 | 1.67 |
| Offline, no mixing | 0.01 | 1.99 | 6.69 | 3.15 |
| Offline, eps 0.5 | 0.23 | 5.12 | 5.77 | 2.01 |
| Offline, eps 1 | 0.07 | 7.72 | 5.59 | 0.60 |
| Offline, eps 2 | 0.00 | 6.89 | 3.56 | 0.24 |

| Land at 0° E, the GCM's columns | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise |
|---|---|---|---|---|
| The branch's convective rain at the ground | 0.92 | 13.28 | 0.46 | 0.08 |
| Offline, no mixing | 0.14 | 15.17 | 1.50 | 0.32 |
| Offline, eps 1 | 0.09 | 26.99 | 0.36 | 0.00 |
| Offline, eps 2 | 0.68 | 52.39 | 0.00 | 0.00 |

- **Without mixing** the offline scheme gives the branch's own pattern at the sea, the evening wettest and the morning
  nearly dry, from noon to sunrise within a factor of two of the rain the branch's scheme formed. It convects on 3–8% of
  the sea's cell means then, and the branch's own scheme formed rain in 4–11% of its snapshots there.
- **Mixing** takes the parcel's buoyancy first where the branch's air sinks. At eps 1 the parcel is buoyant in 5% of
  the sea's cell means after midnight and 28% in the afternoon, against 40% and 62% without mixing. The sea's rain after
  midnight falls to a fifth and its afternoon's rises fourfold; over land the afternoon nearly doubles and the evening
  falls to a quarter, 0.36 mm/day where CM1 rains 1.53.
- **The rh 0.6 branch** answers differently. On its columns a mixing of 0.5 cuts the sea's rain at every hour, after
  midnight from 1.41 to 0.44 mm/day and in the afternoon from 2.42 to 0.38, and the land's afternoon from 19.4 to 8.3.
- **The evidence is thin.** At eps 1 only 6 of the 900 sea cell means after midnight convect, and the five years give
  0.0–2.6 mm/day there.

These numbers are the scheme's first response to columns its branch built without mixing. A branch with mixing would
cool its upper air toward the diluted parcel until convection resumed, and its rain would again follow the supply, so
only such a branch shows where that rain falls. Among the settings tried, the screen on the GCM's own columns points to
eps 1 at rh 0.7; on CM1's air that setting fails the test with the rest.

**A CAPE trigger, tested on the GCM's own columns** (24 October, at the author's request for one more attempt;
`cape_limit` and `eps_trigger` in [sbm.py](sbm.py), `SBM_TRIGGER_PATCHES` in [exoplasim_run.py](exoplasim_run.py),
`trigger_test` in [sbm_offline.json](../results/gcm/sbm_offline.json)). The diagnosis of why the CM1 sea box stays
quiet after midnight ([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e))
weighed the three trigger families against the boxes:
- **Boundary-layer air reaching its LFC.** This would favour the sea's nights: the box's 6-km boundary layer then holds
  the LFC in 59% of its columns, against an inhibition of 0.4 J/kg.
- **The forcing's dCAPE** (Xie and Zhang 2000). It treats the rainy evening and the quiet night alike, at −0.7 to
  −2.4 J/kg/h through both.
- **A CAPE threshold.** This follows the box: its storms start once the afternoon's circulation raises the CAPE above
  about 1,800 J/kg, and they stop as they spend it.

So the trigger is the CAPE family, with a mixing of its own for the dilute form (Neale and others 2008). With
`cape_bm` above 0, a column convects deep only where a second parcel, lifted as the reference's but mixing `epst_bm` per
e-folding of pressure, has at least `cape_bm` J/kg of CAPE. That CAPE is counted as the scheme counts its inhibition,
from the parcel's own LFC up through the layers it stays buoyant, and a column held back takes the shallow scheme.
Off by default, it is a patch of its own on the scheme's text, applied only for a run that sets it (`--sbm-cape`,
`--sbm-trigger-eps`), so runs without it keep their source. The checks ([test_sbm.py](../tests/test_sbm.py)):
- **Fortran against Python.** The patched Fortran routine gives the Python copy's tendencies on four columns with the
  trigger on, to 10⁻¹⁰.
- **With the trigger off,** the routine prints exactly what the scheme printed before on five columns. The refactored
  Python copy matches its earlier self bit for bit on 2,880 of the GCM's columns.
- **The patch** applies once to both ExoPlaSim's own source and the twilight copy's patched one. That was checked in
  memory; neither source was changed.

Five settings, on the rh 0.7 branch's cell columns of years 10–14, beside the scheme without the trigger and with its
reference mixed at eps 1. The offline rain of each quarter is turned into the branch's own rain (pr) times the
setting's offline rain over the scheme's without the trigger; mm/day:

| | Sea, sunrise–noon | Sea, noon–sunset | Sea, sunset–midnight | Sea, midnight–sunrise | Land, noon–sunset | Land, sunset–midnight |
|---|---|---|---|---|---|---|
| CM1 | 0.00 | 2.23 | 2.97 | 0.14 | 8.59 | 1.53 |
| The branch, without the trigger | 0.31 | 2.88 | 3.33 | 1.45 | 15.10 | 0.47 |
| The reference mixed at eps 1 | 1.39 | 11.19 | 2.79 | 0.28 | 26.87 | 0.11 |
| CAPE above 50 J/kg | 0.00 | 0.70 | 1.69 | 0.54 | 9.04 | 0.04 |
| CAPE above 100 J/kg | 0.00 | 0.68 | 1.49 | 0.21 | 8.23 | 0.04 |
| CAPE above 200 J/kg | 0.00 | 0.27 | 0.87 | 0.09 | 6.43 | 0.00 |
| CAPE above 5 J/kg, mixing 1 | 0.00 | 0.06 | 0.13 | 0.00 | 4.05 | 0.04 |
| CAPE above 5 J/kg, mixing 4 | 0.00 | 0.00 | 0.06 | 0.00 | 3.64 | 0.00 |

No setting passes. The test asks for the sea's rain after midnight at most half the branch's and within 0.5 mm/day of
CM1's, with the afternoon and evening at both sites within a factor of two of CM1's. The offline rain itself fails the
same way.
- **The undiluted thresholds cut the nights most.** At 100 J/kg the sea's night falls to 0.21 mm/day and the land's
  afternoon to 8.23, next to CM1's 0.14 and 8.59. But the sea's afternoon falls to 0.68, and the land's evening, already
  a third of CM1's, to 0.04. The branch's afternoon convection fires on columns holding little CAPE, a median of
  32 J/kg on the sea's raining columns against 75 in the evening, because its circulation's ascent destabilizes them
  continuously and the scheme rains what it supplies.
- **The dilute forms find no CAPE.** On every quarter's raining columns the median CAPE of a parcel mixing at 1 is nil:
  the branch's convection relaxes its troposphere toward the undiluted parcel's profile, against which a mixed parcel
  stays cooler. So a dilute threshold stops convection at almost every hour.
- **What the test cannot show.** It takes the scheme's first response to columns a branch built without the trigger.
  A running model would let CAPE build until the threshold passed, fastest where the circulation destabilizes, as it
  does in the afternoon; how much of either quarter's loss would last, only a branch could show.

So the screen ran no branch, and it left open whether the box's quiet nights belonged to its 6 km grid. Boxes at 3 and
2 km then kept the sea's nights quiet
([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)), and the threshold of
about 100 J/kg that the screen pointed to went to a two-year branch, set aside below.

**The trigger's branch, set aside** (25 October; `A28R_301k_cold_sbm_rh07_cape100`, two years from the design climate's
year-19 restart with the trigger at 100 J/kg and no mixing; `sbm_cape_test_analysis.txt` in `session-tools` on the
research drive). Over years 0–1, against the rh 0.7 branch at the same stage, it halved the sea's rain after midnight,
from 2.06 to 1.04 mm/day, where CM1 rains 0.14. It raised the sea's afternoon to 4.51 mm/day against 3.14 (CM1 2.23),
and it left the land's afternoon near 14 mm/day (14.23 against 14.55; CM1 8.6). In year 1 its ground ran 0.6 K cooler,
298.4 against 299.0 K in the global mean. The trigger stays off by default and is not adopted, because nothing
explained why the GCM would need it.

**Why the scheme rains after midnight** (25 October, at the author's go-ahead; no new runs;
[box_entrainment.py](../crm/box_entrainment.py),
[box_entrainment_near_side_sea.json](../results/crm/box_entrainment_near_side_sea.json)). The CM1 box's nights stay
shallow because its cloud mixes in its surroundings at 0.2–1.3 × 10⁻³ m⁻¹ and loses the little buoyancy a parcel has
between 6 and 10 km, beneath the deep layer that holds nearly all its CAPE
([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)). The rh 0.7 branch's six
sea cells over years 10–14 show the scheme's side in two records: the convection's heating at the last step of each
3-day mean (one step every 71.5 hours in each cell, 900 after midnight), and the means themselves, on which
[sbm.py](sbm.py) traces the decision at once.
- **The running scheme rains mostly from convection topped at 4–10 km.** After midnight 36 of the 900 steps rain,
  forming 1.67 mm/day on average against the branch's 1.44 at the ground. Their heating shows the top the scheme chose:
  24 of them heat only the two lowest layers, topped at σ 0.924 (about 4.2 km up), and 4 more reach σ 0.823 (10.2 km).
  So 64% of the rain formed after midnight (45–89% over 300 resamplings of the steps) comes from convection topped at or
  below 10 km, against 62% (49–76%) in the evening and 22% (9–38%) in the afternoon.
- **Its split between deep and shallow weighs moisture alone.** The scheme convects deep wherever its undiluted
  parcel from the lowest level is buoyant at some level above it and the layers from the parcel's top down hold, in
  sum, more vapour than 70% of the parcel's saturation. It has no test of depth: PlaSim's `pdeep`, which sends
  convection topped below a set pressure to the shallow scheme, stands at 999999 Pa, off. When the parcel finds
  buoyancy only at σ 0.924, just above the moist lowest layer, the scheme relaxes the lowest 6.6 km toward 70% of the
  parcel's saturation and rains the vapour above it, 27 mm/day on average while it lasts.
- **The means smooth the steps at which it fires.** Run at once on the 888 cell means after midnight, the scheme finds
  its parcel buoyant on 358, mostly at σ 0.21–0.32 (56–73 km) above a median inhibition of 83 J/kg, and the vapour
  above its reference on only 30, which it rains on. On those 30 the vapour above the reference sums to a median of
  5.5 kg/m². The excess lies at σ 0.924, 12.7 kg/m² there on average, against deficits higher up and in the lowest
  layer.
- **Against CM1's night**, on PlaSim's levels (CM1's domain mean interpolated in log pressure), the GCM is warmer at
  every level: by 1.9 K at σ 0.983, 3.7–3.8 K at σ 0.92–0.82 and 5.5–11.1 K at σ 0.70–0.32, 19–56 km up. Its relative
  humidity is 0.60–0.66 below 30 km against CM1's 0.67–0.95, 0.86 at σ 0.924, where CM1's 6 km boundary layer and its
  cloud stand. The Kuo GCM that forced the box ran 1.2–2.5 K cooler than this branch at every level to σ 0.32. So the
  box holds a cooler boundary layer, moister near its top, under an upper troposphere 8–11 K colder than the GCM's,
  and its undiluted parcel has 1,415 J/kg of CAPE with nil inhibition.
- **A plume at CM1's rates finds no buoyancy on the GCM's columns.** Lifted from PlaSim's lowest level through each mean
  column laid on CM1's levels, dropping its condensate as the scheme's parcel does, a plume mixing at any of CM1's
  measured rates (0.18–1.34 × 10⁻³ m⁻¹, 9–64 per e-folding of pressure as `eps_bm` counts it) is buoyant nowhere on the
  30 night means where the scheme rains at once, nor on 99–100% of the afternoon's 27 and the evening's 74. Undiluted,
  that plume, which leaves out the ice's heat that the scheme's parcel takes, is buoyant nowhere on 80% of the night's
  30 and 41% of the afternoon's 27. The scheme keeps its troposphere near its undiluted parcel's temperature, and in
  such a column CM1's entrainment would stop convection below the scheme's top at every hour, the afternoon's included.

So the scheme rains after midnight in two ways, and the box shows the physics each leaves out. Most of its night rain
comes from convection topped at 4–10 km, which it rains out because the warm sea keeps the boundary layer above 70% of
the parcel's saturation; CM1's boundary layer holds such vapour too, 9.9 kg/m² above the reference on its night columns,
and spreads it into a cloud layer 33–45 hPa deep that barely rains (0.14 mm/day). The rest comes from steps whose
undiluted parcel runs to 19–56 km; in the box's night air a cloud mixing at the measured rates ends at 7–8 km. The
evidence is thin where it counts most: 36 raining steps after midnight in five years, single steps whose columns the
output does not hold, and the means' columns in their place.

**What it points to** (proposals of 25 October; the next part builds the first two and tests them).
- **Decide deep or shallow by the cloud's depth.** The IFS (documentation Cy36r1, 2010, Part IV, section 6.4, read)
  calls a column deep when its test parcel's cloud, from condensation to where its updraught stops, is thicker than
  200 hPa, and shallow below that. If lunar convection is Earth's stretched by the ratio of the gravities, depths in
  pressure carry over, and the threshold falls between CM1's night cloud (33–45 hPa) and its afternoon cloud (215–270
  hPa from its base to its median top). On the running scheme's steps a 200 hPa test would class those topped at
  σ 0.924 as shallow, about 70 hPa above the lowest level and 43% of the rain formed after midnight, and would set the
  four at σ 0.823, about 195 hPa and 21% more, at its edge. PlaSim's `pdeep` already makes that split by the top's
  pressure.
- **Lift the scheme's parcel as CM1's cloud mixes.** In CM1's night air a parcel mixing at the measured rates stops
  by 8.3 km, where the box's cloud ends. On the GCM's present columns it finds no buoyancy at any hour, so the
  scheme would at first stop everywhere, and its troposphere would cool toward the box's until diluted parcels rose
  again. On PlaSim's ten levels a parcel mixing at these rates exchanges 42–98% of its air within the lowest step
  alone, so it would have to be integrated within each layer, as the plume here is on CM1's levels.
- **A mass-flux scheme with entraining plumes** (Tiedtke 1989's family) is the larger change: its shallow branch would
  carry the boundary layer's vapour into a cloud layer and detrain it there, as CM1's night does, and its deep branch
  could start from the afternoon's measured 0.16–0.33 × 10⁻³ m⁻¹.

Each of these changes the column a running model builds, so the first response measured here does not show where the
rain would then fall; a branch would.

**The cloud-depth test and the mixing within layers** (25 October, at the author's approval; `SBM_DEPTH_PATCHES` in
[exoplasim_run.py](exoplasim_run.py), `deep_pa` and `steps` in [sbm.py](sbm.py),
[sbm_depth_mixing.json](../results/gcm/sbm_depth_mixing.json); `run_depth_tests.sh` and `sbm_depth_test_analysis.py`
in `session-tools/agents/climate-scheme` on the research drive). The first two proposals are built as two options of
the scheme, each off by default and recorded in a run's configuration only when set:
- **The cloud-depth test** (`--sbm-deep-hpa`; `dpdeep_bm` in Pa). A column convects deep only where the parcel's cloud,
  from its base to the scheme's top level, is thicker than the threshold; a thinner cloud goes to PlaSim's shallow
  scheme, which mixes temperature and vapour through the layers it spans and forms no rain. The base is the lifting
  condensation level, the IFS's cloud base: the pressure at which the lowest level's air, lifted dry adiabatically
  with its vapour kept, reaches PlaSim's saturation, found by bisection in log pressure. The top is the scheme's own,
  the last level at which its parcel stays buoyant. The IFS takes the level where its updraught's speed vanishes, at
  or above the end of the buoyancy, and its test parcel mixes and may start above the lowest level. PlaSim's `pdeep`
  tests the top's pressure alone, so over high ground, deep basins or a high cloud base one level would mark clouds of
  different depths; the test is a parameter of its own and `pdeep` stays off.
- **The mixing within layers** (`--sbm-mixing-steps`; `nsub_bm`, with `--sbm-eps`). The parcel crosses each layer in
  equal steps of log pressure; each step lifts it dry adiabatically, mixes in its share of the surrounding air,
  interpolated linearly in log pressure, and condenses the mixture. Over a layer the parcel still exchanges a share
  1 − exp(−ε ln(p_below/p)) of its air, one step a layer is the scheme as before, and a mixing of 0 never takes the
  steps.

**The checks** ([test_sbm.py](../tests/test_sbm.py)).
- **Fortran against Python.** Compiled with ExoPlaSim's own options (`-O3 -fcheck=all -ffpe-trap=invalid,zero,overflow
  -fdefault-real-8`), the routine gives the Python copy's tendencies and rain to 10⁻¹⁰ on 16 columns with the options
  on, alone, together and with the CAPE trigger, with the two patches applied in either order, and hands the same
  columns to the shallow scheme at the same top.
- **Both options off.** The routine prints every digit the scheme printed before on 25 columns, the CAPE trigger's
  cases among them, with or without the trigger in the source; 17 significant digits fix every bit.
- **The patch** applies once to both installed sources, the twilight copy's, which carries the trigger, and
  ExoPlaSim's own, which has neither patch, with the trigger before it, after it or absent. It goes into a source only
  for a run that sets either option, as the trigger's does; the first test put it into the twilight copy.
- **In the GCM.** One year with neither option, `A28R_301k_cold_sbm_rh07_offcheck`, compiled from the twilight copy
  with the scheme, the trigger and this patch in its source, reproduces the rh 0.7 branch's year 0, whose source held
  the scheme alone: all 61 output variables and the year's restart file match bit for bit. Runs that set neither
  option compute as before. The source's hash in their run key is new (`rainmod` aa2a624377cd8b58), so a run begun
  before the patch cannot resume, as after the trigger's.

**The steps converge at first order.** Against 1,024 steps a layer, the reference's largest error below σ 0.2 halves
with each doubling of the steps:

| Largest error below σ 0.2 (K), steps a layer | 1 | 2 | 4 | 8 | 16 | 32 | 64 |
|---|---|---|---|---|---|---|---|
| The rh 0.7 branch's sea column (years 10–14), mixing 2.4 | 0.98 | 0.55 | 0.30 | 0.15 | 0.08 | 0.04 | 0.02 |
| The same, mixing 4.8 | 0.77 | 0.48 | 0.27 | 0.14 | 0.07 | 0.04 | 0.02 |
| CM1's sea box, day 3, on PlaSim's levels, mixing 2.4 | 2.04 | 1.14 | 0.60 | 0.31 | 0.15 | 0.08 | 0.04 |
| The same, mixing 4.8 | 2.25 | 1.37 | 0.75 | 0.39 | 0.20 | 0.10 | 0.05 |

The tests take 16 steps a layer, which hold the reference within 0.08 K of the limit on the GCM's column and 0.2 K on
the box's colder one; the top two levels, above the convecting troposphere, keep errors of 0.7 K. In the GCM the 16
steps cost a quarter more time a year (below).

**On CM1's mean column the integrated reference reproduces the box's air** (the sea box's third lunar day, laid on
PlaSim's levels at its surface pressure, at box_upper_air's heights; the reference minus the box's air, K):

| At | 18.9 km | 29.8 km | 41.8 km | 55.8 km |
|---|---|---|---|---|
| The box's air | 277.18 | 265.39 | 251.49 | 231.76 |
| Undiluted reference | +3.24 | +5.13 | +5.90 | +6.78 |
| Mixing 2.4, 16 steps | +1.30 | +1.32 | +0.50 | +0.39 |
| Mixing 4.8, 16 steps | +0.55 | +0.28 | −0.41 | +0.06 |
| Mixing 2.4 and 4.8, one step | +1.21, +0.13 | +1.46, +0.16 | +0.62, −0.40 | +0.64, +0.24 |
| Zero-buoyancy plumes at 0.5 and 1 × 10⁻⁴ m⁻¹ ([box_upper_air.py](../crm/box_upper_air.py)) | +1.72, +0.45 | +2.56, +0.58 | +2.21, −0.66 | +1.97, −1.73 |

At 2.4 the reference lies within 1.3 K of the box's air at 19–56 km and at 4.8 within 0.6 K, where the plumes at the
matching rates lie within 2.6 and 1.7 K. On this column the steps move the reference by at most 0.4 K at these
heights; on the GCM's warmer, moister column one step errs by up to 1.0 K below σ 0.2.

**Three two-year tests** (from the design climate's year-19 restart with every flag of the rh 0.7 branch and the depth
test at 200 hPa: alone, `A28R_301k_cold_sbm_rh07_deep200`; with mixing 2.4, `…_deep200_eps24`; with mixing 4.8,
`…_deep200_eps48`; the mixing at 16 steps a layer). Each is set beside the rh 0.7 branch over the same years from the
same restart, and beside the settled CM1 boxes (the sea box's third lunar day, the land box's second). Rain by local
time, mm/day, over years 0–1 and, after the semicolon, year 1:

| Near-side Sea | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise | Mean |
|---|---|---|---|---|---|
| rh 0.7 branch | 0.51; 0.66 | 3.14; 2.68 | 2.61; 1.95 | 2.06; 2.04 | 2.06; 1.81 |
| Depth test | 0.25; 0.19 | 3.66; 3.15 | 3.07; 2.66 | 1.55; 1.46 | 2.11; 1.85 |
| Depth test, mixing 2.4 | 0.20; 0.06 | 1.96; 1.44 | 0.94; 0.52 | 0.54; 0.44 | 0.90; 0.60 |
| Depth test, mixing 4.8 | 0.18; 0.26 | 1.05; 2.04 | 0.58; 0.73 | 0.66; 0.70 | 0.61; 0.92 |
| CM1 box | 0.00 | 2.23 | 2.97 | 0.14 | 1.34 |

| Land at 0° E | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise | Mean |
|---|---|---|---|---|---|
| rh 0.7 branch | 0.95; 0.99 | 14.55; 15.42 | 0.48; 0.65 | 0.14; 0.11 | 4.09; 4.42 |
| Depth test | 0.80; 0.83 | 15.79; 17.95 | 0.45; 0.47 | 0.06; 0.05 | 4.34; 4.97 |
| Depth test, mixing 2.4 | 0.60; 0.92 | 16.87; 17.16 | 0.71; 1.06 | 0.02; 0.02 | 4.62; 4.93 |
| Depth test, mixing 4.8 | 0.35; 0.47 | 27.12; 19.72 | 0.05; 0.06 | 0.06; 0.06 | 7.01; 5.24 |
| CM1 box | 2.87 | 8.59 | 1.53 | 0.07 | 3.27 |

| Years 0–1 | rh 0.7 branch | Depth test | Mixing 2.4 | Mixing 4.8 | CM1 box |
|---|---|---|---|---|---|
| Sea site at 41.8 km; 55.8 km | 260.06; 241.91 K | 259.90; 241.68 | 258.11; 239.58 | 257.85; 239.46 | 251.49; 231.76 |
| The same in year 1 | 260.30; 242.16 K | 260.06; 241.85 | 258.22; 239.58 | 257.80; 239.35 | |
| Ground; air at 2 m | 298.85; 299.09 K | 298.46; 298.50 | 300.57; 300.62 | 301.26; 301.43 | |
| Rain; its convective share | 2.23 mm/day; 94% | 2.18; 94% | 1.95; 95% | 1.87; 91% | |
| Imbalance at the top of the air | +2.35 W/m² | +1.45 | +2.63 | +2.84 | |
| Planetary albedo; cloud cover | 0.376; 35.8% | 0.380; 35.4% | 0.360; 33.8% | 0.346; 32.4% | |
| Cloud effect: sunlight; infrared | −39.4; 50.7 W/m² | −40.5; 50.0 | −34.5; 49.3 | −30.2; 48.4 | |
| Top two levels: air | 164.30, 193.46 K | 164.29, 192.85 | 163.24, 192.15 | 163.18, 192.46 | |
| Top two levels: water | 0.021, 3.28 ppmv | 0.022, 2.97 | 0.016, 2.52 | 0.015, 2.61 | |

The global numbers for year 1 alone, and every level of the sea site's column, are in the results file.

**What each change does.**
- **The depth test** hands the thin clouds to the shallow scheme. The sea's rain after midnight falls by a quarter, to
  1.55 mm/day (1.46 in year 1), and its morning's by half; all that is left after midnight is convective, so it comes
  from clouds deeper than 200 hPa. The vapour those clouds held back rains later in the day, about 0.5 mm/day more in
  each of the sea's afternoon and evening, the evening reaching CM1's 2.97 and the afternoon moving away from CM1's
  2.23. Over land the afternoon rises to 15.8 mm/day and the night falls to CM1's 0.06–0.07. The upper air moves by
  0.2 K, the ground runs 0.4 K cooler, and the water at σ 0.119 falls by a tenth.
- **The mixing** cools the reference aloft. Within two years the sea site's air at 42–56 km cools by 2.0–2.5 K
  toward the box's, a quarter of the 8.6–10.2 K gap, while the air below 19 km warms by 1.0–3.8 K (σ 0.70–0.98, most
  at 10 km), so at 18.9 km the gap to the box grows from 5.0 to 5.9 and 6.6 K. The lowest level's vapour at the site
  rises from 12.1 g/kg to 14.1 and 14.5 (as mixing ratios), against 11.8 in the box's lowest 0.5 km. Convection weakens
  at the sea at every hour. At 2.4 the sea rains 0.54 mm/day after midnight, 0.20 of it convective, and 0.94 in the
  evening against CM1's 2.97; at 4.8 its convection stops after midnight and in the morning, and the night's
  0.66 mm/day comes from large-scale condensation. Over land the afternoon storms grow, to 16.9 and 27.1 mm/day (17.2
  and 19.7 in year 1) against CM1's 8.6. The evening rises toward CM1's 1.53 at 2.4 (0.71; 1.06 in year 1) and all but
  stops at 4.8.
- **The ground warms as the convective cloud thins.** PlaSim draws its convective cloud cover from the convective
  rain (`mkclouds`: 0.245 + 0.125 ln of the rate in mm/day, spread over the convecting layers). With mixing the cloud
  cover falls by 2.0 and 3.4 points and the cloud's shading of sunlight by 4.9 and 9.3 W/m², while the ground warms by
  1.7 and 2.4 K; the thinner cloud as the warming's cause is an inference from these numbers.
- **The top two levels**, for the exobase check, cool by about 1 K with mixing and hold a quarter less water, 0.015–
  0.016 ppmv at σ 0.038 against 0.021; the depth test alone leaves σ 0.038 as it was and cools σ 0.119 by 0.6 K.

No test reproduces CM1's day at both sites. The depth test removes the thin clouds' rain and leaves the upper air as it
was; the mixing moves the upper air toward the box's but takes the sea's evening rain, sharpens the land's afternoon
and warms the air near the ground. Each test is one realization over the first two years of a spin-up that the branch
needed 8–10 years to settle, with the imbalance at the top still 1.5–2.8 W/m².

**What the runs cost.** The depth test took 258 and 245 s a year and the run without either option 250 s, beside the
branch's 270 and 284. At 16 steps the mixing took 326 and 316 s at 2.4 and 314 and 305 at 4.8, a quarter more. The
machine was shared, so the times are approximate. The four runs hold 1.55 GB on the research drive.

**PlaSim's cap on shallow tops.** The shallow scheme holds its tops at or below the lowest level above 700 hPa
(`pdeepth`, PlaSim's Earth value). A cloud the depth test hands it could meet that cap only with its base above
900 hPa. The lowest air of 1.0–1.4% of the cell means is still unsaturated at 900 hPa without the mixing, 4.6% at 2.4
and 12.8% at 4.8; how many of those columns the depth test sends to the shallow scheme the means do not show.

**What remains.**
- **The longwave refit.** In the mixing tests the sea site's column cools by up to 2.9 K at 29–72 km and warms by up
  to 3.8 K at 1–19 km, away from the columns the 301 K refit was fitted on, and on CM1's colder air RRTMG cools
  35–63 km 17–18 W/m² more than PlaSim's correction. A mixing that is adopted takes a refit of the longwave
  correction.
- **The air near the ground.** The GCM's lowest level holds more vapour than the box's, worth 3–4 K aloft, and no
  mixing setting reaches it; the mixing tests raise it further.
- **The 15-year run** of any setting is the author's call after this stage.

**Sources.** The IFS documentation (Cy36r1, 2010, Part IV, section 6.4) was read from the text on the research drive
(`session-tools/agents/climate-explain/ifs_cy36r1_part4.txt`): the LCL as the cloud base, the 200 hPa depth from it to
where the updraught's speed vanishes, the deep test parcel's mixing and departure levels. The zero-buoyancy plume
(Singh and O'Gorman 2013) is box_upper_air's calculation and was not read here. PlaSim's own `rainmod.f90` (`kuo`,
`mkshallow`, `mkclouds`) and ExoPlaSim's compiler options (`most_compiler`) were read in the installed copy.

```sh
python3 -m pytest climate/tests/test_sbm.py               # with gfortran: the Fortran against Python, the options off
S=/media/projectspace/terluna-research/session-tools/agents/climate-scheme
$S/smoke.sh deep200 150                                    # a smoke test of a variant: start, check, stop
$S/run_depth_tests.sh deep200                              # about 8 minutes on 4 workers; rerun to resume
$S/run_depth_tests.sh deep200_eps24                        # and deep200_eps48, about 11 minutes each
$S/run_depth_tests.sh offcheck                             # one year without either option, for the bit check
climate/gcm/.venv/bin/python $S/sbm_depth_test_analysis.py > $S/sbm_depth_test_analysis.txt   # the results file
```

#### Fourteen levels near the ground (26 October)

**Why** (the author's approval of 26 October). PlaSim's ten levels stand near 0.9, 4.0, 9.9, 18.2, 29, 42, 59, 79, 108
and 166 km at the design air's 51 km scale height. CM1's settled sea box holds a boundary layer 6 km deep and, at night,
cloud from its level of free convection at 6.2 km to tops at 7–8 km
([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)). Both lie inside the one
layer between 6.4 and 13.7 km, so no convection scheme on the ten levels can represent the box's night or the boundary
layer's ventilation.

**The level set** ([levels.py](levels.py), `l14`). PlaSim puts each full level midway in sigma between the half levels
that bound its layer, so the half levels fix the full levels. The set keeps PlaSim's own six upper layers, from 23.3 km
up, and its lowest layer, from the ground to 1.73 km. Six new half levels at 15.5, 11.5, 9.0, 7.0, 5.0 and 3.1 km (sigma
0.7379, 0.7981, 0.8382, 0.8717, 0.9066 and 0.9410) put the full levels between them at 19.2, 13.5, 10.2, 8.0, 6.0, 4.0
and 2.4 km:

| Full levels, km at the 51 km scale height, the top first | |
|---|---|
| PlaSim's ten | 166, 108, 79, 59, 42, 29, 18.2, 9.9, 4.0, 0.86 |
| The fourteen | 166, 108, 79, 59, 42, 29, 19.2, 13.5, 10.2, 8.0, 6.0, 4.0, 2.4, 0.86 |
| The fourteen on CM1's sea box (its day-3 mean, at equal shares of its surface pressure) | 127, 91, 70, 54, 40, 28, 19.3, 13.7, 10.5, 8.2, 6.2, 4.2, 2.5, 0.90 |

- **What stays.** The upper air, the top two levels the exobase check reads, the three layers the cold cloud's factors
  describe, and the lowest level, where the surface fluxes and the 10 m wind are taken.
- **What moves.** The level at 18.2 km goes to 19.2 km and the one at 9.9 km to 10.2 km. With the half level at
  23.3 km fixed, a full level at 18.2 km would put the next half level at 13.7 km, so the next full level would fall
  near 12 km under a layer three times as thick as its own. Moving the 18.2 km level lets the layers thin smoothly,
  never to less than half of the one above.
- **On CM1's sea box** the new levels near 6 and 8 km stand at 6.2 and 8.2 km, at the box's level of free convection and
  just above its night cloud's tops at 7.3–7.9 km, with a half level between them at 7.2 km.

**The options**, each off by default and recorded in a run's configuration only when set. `--levels l14` sets the layer
count and writes the half levels (NEQSIG −1 and SIGH) and the Rayleigh friction PlaSim presets on its top two levels
only on its ten-level grid (20 and 100 days, TFRC), since those two levels are the same. `--sbm-shallow-cap` is the cap
fix below. The runner refuses a layer count other than its level set's, so the per-layer inputs always match the
levels. With neither option, the configuration the runner builds for the flags of `A28R_301k_cold_sbm_rh07_deep200`
matched that run's record in every entry, checked before the cap's patch entered the source and changed the source hash
every later run records.

**The inputs on fourteen levels**, each made by the script that made the ten-level one.
- **The twilight** ([twilight_301k_l14.json](../results/gcm/twilight_301k_l14.json); `twilight levels`). The 301 K
  ledger keeps its terms by shell as well as by layer, so its curved light, cone and infrared escape are shared among
  the new layers by their overlap in sigma, as the ledger shares them among the ten. The same arithmetic on the ten
  levels gives the ledger's own terms back within 6 × 10⁻⁸ per unit of sunlight at the beam, against terms up to
  4.6 × 10⁻³, its infrared escape within 2 × 10⁻⁵ W/m² and its layer temperatures within 0.001 K. On the fourteen
  levels the shared layers keep their terms, the new layers sum to what the ten-level layers they replace held, and the
  column's totals stay at 5.36 W/m² of curved light, 4.95 of cone and 19.82 of infrared escape. The ledger's shells are
  7.5 km thick near the ground, so the layers below 7.5 km take one shell's terms per unit mass and its temperature,
  298.5 K, as the ten-level layers there did.
- **The longwave correction** ([longwave_A28_corrected_970_lw_ice_l14.json](../results/gcm/longwave_A28_corrected_970_lw_ice_l14.json);
  `plasim_longwave levels`). The 301 K refit's extra escape per unit mass is carried by the same remap as the restart
  (below): the shared layers keep their coefficients, the new layers' sum to the old ones' (−29.46 W/m² over land and
  −31.65 over sea), and the reference temperatures are remapped as a profile. It is not fitted on fourteen-level
  columns.
- **The cold cloud's factors** ([cold_cloud_A28R_301k_s0897.json](../results/gcm/cold_cloud_A28R_301k_s0897.json))
  scale the top three layers, which both grids share, so they carry over unchanged. The runner now checks a product's
  layer bounds against the run's half levels and refuses a level set that changes them.

**The restart** (`levels.py restart`; `gcm-runs/A28R_301k_cold_l14_restart/MOST_REST.00019` on the research drive).
Sigma counts mass, each layer holding ps Δσ/g, so the year-19 restart is carried by a conservative remap in sigma. Each
old layer's profile is rebuilt linearly about its full level with the centred slope of its neighbours, and each new
layer takes the mean of that profile over its sigma range. The surface pressure, and with it every column's mass, stays
as it was, and every column keeps its integrals of water, temperature times mass and momentum to within 10⁻¹³ of
their size. The remap is one linear operator in the vertical, so it acts on the spectral vorticity, divergence,
temperature and humidity of both time levels directly, and the shared layers are copied bit for bit. The restart's
accumulated means are zero at a year's end, and the first step's cloud cover and cloud water are clipped to their
ranges. On year 19's last twelve 3-day means the remapped temperature at the new levels lies within 0.4 K of a plain
interpolation in ln sigma on average (2.2 K at most, at 2.4 km, where the old layer's slope spans the boundary layer),
its vapour up to 3% below the interpolation's, which would add water, and the hydrostatic thickness of the regridded
layers changes by 4–5 m.

**The shallow cap** (`SBM_CAP_PATCHES`, `--sbm-shallow-cap`, with `--sbm-deep-hpa`). PlaSim's shallow scheme holds every
shallow cloud's top at or below the lowest level above 700 hPa (`pdeepth`), Earth's lower troposphere and about 28 km
in the Moon's air. With the option, a column the new scheme hands to the shallow scheme after the depth test (too dry
to rain, or held back by its inhibition or the CAPE trigger) gets a top no higher than the depth test's threshold above
its own cloud base, the condensation level, and PlaSim's absolute cap gives way to it. A cloud the depth test itself
sends over is already that thin. It is a patch of its own, applied only for a run that sets it, after the depth patch.
- **Its tests** ([test_shallow_cap.py](../tests/test_shallow_cap.py)). With the option off the routine prints every
  digit it printed before on ten columns, with and without the CAPE trigger; with it on, the Fortran hands the shallow
  scheme the tops the Python copy ([sbm.py](sbm.py), `shallow_cap`) gives, six of them capped. The patch applies once
  to both installed sources after the depth patch.
- **In the GCM.** One ten-level year with the flags of `A28R_301k_cold_sbm_rh07_deep200` and neither new option,
  `A28R_301k_cold_sbm_rh07_deep200_offcheck`, built after the fourteen-level executable from the source holding the cap
  (`rainmod` 9db0dc616690cb6f), reproduces that run's year 0 bit for bit: all 61 output variables and the year's restart.
  So ten-level runs build and run as before. A run begun before the patch reached the twilight copy cannot resume, its
  source hash being new, as after the depth patch.

**Cost and stability.** One model month from the year-19 restart took 34.9 s at fourteen levels and 22.5 s at ten, and
the two-year test took 414 and 441 s a year, against 278 s for a ten-level year the same morning (the machine was
shared, so the times are approximate): about one and a half times the cost. A year's output grows from 130 to 172 MB.
Neither the month nor the two years shows a fault: every field is finite, the humidity is nowhere negative, the
strongest wind at any level stays under 26 m/s, and the temperature falls smoothly through the new levels.

**Two years at fourteen levels** (`A28R_301k_cold_sbm_rh07_deep200_l14cap`, from the carried restart with every flag of
`A28R_301k_cold_sbm_rh07_deep200` and the shallow cap; `levels_test_analysis.py` in `session-tools/agents/climate-levels`
on the research drive, [levels_l14_test.json](../results/gcm/levels_l14_test.json)). It is set beside the ten-level run
with the same flags, the rh 0.7 branch over the same years from the same restart, and the settled CM1 boxes. A second
fourteen-level run without the cap, `…_l14`, tells the levels' part of the change from the cap's. Years 0–1, and after
the semicolon year 1:

| Near-side Sea, mm/day | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise | Mean |
|---|---|---|---|---|---|
| 14 levels with the cap | 0.05; 0.02 | 2.87; 2.57 | 0.88; 0.52 | 0.09; 0.02 | 0.95; 0.77 |
| 14 levels, no cap | 0.04; 0.02 | 2.70; 2.70 | 1.08; 1.15 | 0.13; 0.19 | 0.97; 1.00 |
| 10 levels (deep200) | 0.25; 0.19 | 3.66; 3.15 | 3.07; 2.66 | 1.55; 1.46 | 2.11; 1.85 |
| rh 0.7 branch | 0.51; 0.66 | 3.14; 2.68 | 2.61; 1.95 | 2.06; 2.04 | 2.06; 1.81 |
| CM1 box | 0.00 | 2.23 | 2.97 | 0.14 | 1.34 |

| Land at 0° E, mm/day | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise | Mean |
|---|---|---|---|---|---|
| 14 levels with the cap | 0.01; 0.00 | 13.33; 13.05 | 0.54; 0.37 | 0.01; 0.00 | 3.52; 3.46 |
| 14 levels, no cap | 0.09; 0.17 | 14.74; 16.88 | 0.25; 0.08 | 0.00; 0.00 | 3.83; 4.43 |
| 10 levels (deep200) | 0.80; 0.83 | 15.79; 17.95 | 0.45; 0.47 | 0.06; 0.05 | 4.34; 4.97 |
| rh 0.7 branch | 0.95; 0.99 | 14.55; 15.42 | 0.48; 0.65 | 0.14; 0.11 | 4.09; 4.42 |
| CM1 box | 2.87 | 8.59 | 1.53 | 0.07 | 3.27 |

| Years 0–1 | 14 levels, cap | 14 levels, no cap | 10 levels (deep200) | rh 0.7 branch | CM1 box |
|---|---|---|---|---|---|
| Sea site at 41.8; 55.8 km | 259.32; 241.20 K | 259.30; 241.22 | 259.90; 241.68 | 260.06; 241.91 | 251.49; 231.76 |
| Sea site's lowest level: vapour; air | 12.88 g/kg; 300.70 K | 12.92; 300.80 | 12.10; 299.39 | 12.14; 299.73 | 11.8 g/kg (lowest 0.5 km) |
| Sea site: air at 2 m; ground | 302.15; 302.51 K | 302.25; 302.57 | 300.83; 301.83 | 301.16; 301.71 | |
| Ground; air at 2 m | 299.51; 299.97 K | 299.54; 300.01 | 298.46; 298.50 | 298.85; 299.09 | |
| Rain; its convective share | 2.17 mm/day; 88% | 2.17; 88% | 2.18; 94% | 2.23; 94% | |
| Imbalance at the top of the air | +1.40 W/m² | +1.63 | +1.45 | +2.35 | |
| Planetary albedo; cloud cover | 0.378; 36.3% | 0.377; 36.0% | 0.380; 35.4% | 0.376; 35.8% | |
| Cloud effect: sunlight; infrared | −39.9; 49.7 W/m² | −39.6; 49.7 | −40.5; 50.0 | −39.4; 50.7 | |
| Top two levels: air | 164.04, 192.66 K | 164.04, 192.69 | 164.29, 192.85 | 164.30, 193.46 | |
| Top two levels: water | 0.0204, 2.93 ppmv | 0.0203, 2.94 | 0.0217, 2.97 | 0.0213, 3.28 | |

The vapour is a mixing ratio, the sea site's numbers are lunar-day means over its cells, and the global numbers for
year 1 alone and the sea site's column level by level are in the results file.

**What the levels change.**
- **The scheme convects deep far less.** On the sea site's 3-day means the parcel from the lowest level is buoyant in
  78–92% of them at fourteen levels in each quarter of the day, against 23–53% at ten, its top mostly at 60–80 km. In
  73–92% of all the means the column then holds less vapour than 70% of the parcel's saturation, against 13–48% at
  ten, so the scheme hands it to the shallow scheme, which forms no rain; the depth test's thin clouds make up 0–5%
  (`night_tops.py` beside the analysis). The means smooth the steps at which the scheme fires, so these shares describe the columns' usual state.
  Globally the convective share of the rain falls from 94% to 88%, the rain itself staying at 2.17 mm/day.
- **The sea's nights go quiet by that route.** After midnight the sea rains 0.09–0.13 mm/day beside CM1's 0.14, where
  the ten levels rained 1.55, and its mornings fall to CM1's near-nothing. In 88% of the means after midnight the parcel
  is buoyant to tens of kilometres and the column too dry to rain; CM1's night cloud stops at 7–8 km as it mixes in its
  surroundings, so the GCM reaches the box's rain here through its reference humidity, a different mechanism.
- **The sea's evenings lose most of their rain**, 0.88–1.08 mm/day against CM1's 2.97, which the ten levels matched
  (3.07), and the afternoon comes down to 2.70–2.87, nearer CM1's 2.23. So the sea's mean rain halves, to 0.95 mm/day
  against CM1's 1.34.
- **Over land** the afternoon falls from 15.8 to 13.3–14.7 mm/day, toward CM1's 8.6, while the mornings' 0.8 mm/day
  vanishes against CM1's 2.9 and the evening stays at 0.25–0.54 against CM1's 1.5.
- **The air near the ground warms and moistens.** The sea site's lowest level runs 1.3–1.4 K warmer and holds 0.8 g/kg
  more vapour, away from CM1's 11.8 g/kg, and the global ground and air at 2 m warm by 1.05 and 1.5 K. The physics now
  leaves vapour worth 44 W/m² of latent heat in the air below 5 km, against 23 at ten levels, and takes it out at
  5–12 km, where the deep convection rained it out of the ten levels' layer at 6–14 km (`low_budget.py` beside the
  analysis). The first month already showed the lowest level 1.0 K warmer, so the change sets in within weeks; whether
  the moister boundary layer leaves the columns too dry for the reference or follows from it, the runs do not separate.
- **Aloft** the sea site's air cools by 0.5–1.2 K above 13 km, 0.6 K at 41.8 km and 0.5 K at 55.8 km, a twentieth of
  the way to CM1's box, as the deep convection's heating weakens. The top two levels cool by 0.2–0.25 K and the
  topmost holds 6% less water.
- **The cap moves little of this.** The two fourteen-level runs differ by under 0.05 K in the global ground and under
  0.15 K in the upper air. In year 1 the cap takes the sea's evening further from CM1 (0.52 against 1.15 mm/day) and
  brings the land's afternoon toward it (13.1 against 16.9 mm/day).

Each run is one realization over the first two years of a spin-up that the rh 0.7 branch needed 8–10 years to settle,
with the imbalance at the top still 1.1–1.6 W/m² in year 1, so these show the direction and early size of each change.

**The longwave correction on fourteen-level columns** ([longwave_check_l14.json](../results/gcm/longwave_check_l14.json);
`plasim_longwave check --case … --correction …`). On the capped test's year 1 the carried 301 K refit leaves a largest
layer mean of 6.71 W/m² over sea and 6.51 over land on 1,400 cloud-free columns, and 5.73–5.94 on a sample of 1,274 of
every column, against the rule's 2.5. PlaSim alone leaves 12.7–13.4. The errors lie in the layers at 19 and 29 km, which
cool 4.1–6.7 W/m² too much, the second of them a layer both grids share: PlaSim's cooling of a layer depends on the
layers around it, so the correction fitted on ten levels does not carry over. The same check of the ten-level refit on
the ten-level run's year 1 leaves 3.31 W/m² (`longwave_check_deep200.json` beside the analysis on the research drive). A
refit on the test's own cloud-free columns, made only to see what one would leave and not adopted, holds the sample of
every column within 1.3 W/m²; with the correlated-k targets cached it takes about a minute
(`longwave_refit_probe.py`, beside the analysis).

**What the runs cost.** The two fourteen-level runs took 414 and 441 s and 514 and 454 s a year, the second while the
longwave checks shared the machine; the month trials 35 and 23 s; the ten-level check year 278 s. All of this work
holds 1.4 GB on the research drive.

**What remains before a 15-year run.**
- **The scheme's reference on these levels.** At rh 0.7 the reference now leaves most deep, buoyant columns too dry to
  rain. That quiets the sea's nights, which CM1's cloud quiets by mixing, and takes the evening's rain with it. Whether the
  levels then need another reference humidity, the parcel's mixing within layers or the mass-flux scheme under
  construction (which can take the set from `levels.py`) is the next question, and the boundary layer's warming and
  moistening needs tracing with it.
- **The longwave correction** needs a refit on fourteen-level columns of the configuration chosen for that run; the
  carried one fails the rule by 4 W/m² in the layers at 19–29 km.
- **The cap** moves little at fourteen levels; adopting it, and the levels themselves, is the author's call.
- **The cold cloud's factors** describe the top three layers of a ten-level Kuo run, as for the other branches; they
  carry over unchanged, and a settled fourteen-level design climate would let them be refitted with the rest.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.levels show l14                     # the set, its heights and layers
climate/gcm/.venv/bin/python -m climate.gcm.levels restart \
    /media/projectspace/terluna-research/gcm-runs/A28R_301k_cold/model/MOST_REST.00019 \
    /media/projectspace/terluna-research/gcm-runs/A28R_301k_cold_l14_restart/MOST_REST.00019 --to l14   # seconds
climate/gcm/.venv/bin/python -m climate.gcm.twilight levels --ts 301 --levels l14             # twilight_301k_l14.json
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave levels --levels l14    # the refit carried onto the levels
L=/media/projectspace/terluna-research/session-tools/agents/climate-levels
$L/run_levels_tests.sh month14                   # and month10: one model month each, for the cost
$L/run_levels_tests.sh offcheck10                # a ten-level year with the cap's patch in the source, for the bit check
$L/run_levels_tests.sh test14                    # two years, about 15 minutes on 4 workers; test14nocap without the cap
climate/gcm/.venv/bin/python $L/levels_test_analysis.py > $L/levels_test_analysis.txt   # the results file
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave check --case A28R_301k_cold_sbm_rh07_deep200_l14cap:1 \
    --correction refit_301k_l14=climate/results/gcm/longwave_A28_corrected_970_lw_ice_l14.json \
    --output climate/results/gcm/longwave_check_l14.json                          # about a minute
```

#### Fifteen levels: a thin bottom layer (26 October)

**Why** (the author's approval of 26 October). The offline mass-flux scheme ([massflux.md](massflux.md)) lifts its
updraughts from PlaSim's lowest layer, which spans the ground to 1.73 km on both the ten- and the fourteen-level grids.
That layer averages away the near-ground air CM1's updraughts carry: it holds 1.4–2.3 kJ/kg less moist static energy
than CM1's air at 50 m. A thin bottom layer hands a scheme that air directly.

**The level set** ([levels.py](levels.py), `l15`). The set keeps every half level of `l14` and adds one at σ 0.9883,
0.60 km at the 51 km scale height, splitting the lowest layer in two. PlaSim puts each full level midway in σ between
its half levels, so the lowest level stands at σ 0.99415, 0.30 km, and the next at σ 0.97745, 1.16 km; the thirteen
levels above them are the fourteen-level set's own.

| Full levels, km at the 51 km scale height, the top first | |
|---|---|
| The fourteen | 166, 108, 79, 59, 42, 29, 19.2, 13.5, 10.2, 8.0, 6.0, 4.0, 2.4, 0.86 |
| The fifteen | 166, 108, 79, 59, 42, 29, 19.2, 13.5, 10.2, 8.0, 6.0, 4.0, 2.4, 1.16, 0.30 |
| The fifteen on CM1's sea box (its day-3 mean, at equal shares of its surface pressure) | 127, 91, 70, 54, 40, 28, 19.3, 13.7, 10.5, 8.2, 6.2, 4.2, 2.5, 1.21, 0.31 |

With one half level added, 0.3 km is close to the lowest the bottom level can stand while the bottom layer keeps at
least half the mass of the layer above it, the limit the fourteen-level set kept between neighbours: the limit falls
at a half level of σ 0.9889 and a lowest level of 0.28 km, and this set's ratio is 0.54. The bottom layer holds about
850 kg/m² of air, a third of the fourteen-level bottom layer's 2,420 kg/m² and two and a half times what PlaSim's
lowest layer holds on Earth.

**How PlaSim's surface layer takes the thin layer** (`fluxmod.f90` read in the twilight copy, beside the runs below;
`surface_layer_check.py` beside the analysis).
- **Height and drag.** Each step `surflx` works out the lowest level's height from that level's temperature and the
  ground's pressure, z_L = −R ½(T_L + θ_v) ln σ_L / g, with the Moon's gravity and the air's gas constant: 0.31 km at
  fifteen levels against 0.90 km on the coarser grids. The neutral transfer coefficient is the log law to that height,
  (κ / ln(z_L/z0 + 1))², and the stability factors are ECHAM's Louis (1979) functions of the bulk Richardson number
  g z_L Δθ_v / (θ_v |U|²) and of z_L/z0. Every term takes the level's own height and the Moon's gravity, so the
  coefficients follow the level down: the neutral coefficient rises from 0.0043 to 0.0063 over land (on PlaSim's
  preset z0 of 2 m) and from 0.00059 to 0.00068 over sea (on Charnock's z0), while the wind it multiplies falls. For a
  given wind and temperature difference a level 0.3 km up on the Moon has the Richardson number of a level 50 m up on
  Earth, near the first levels of the Earth models those functions were made for; the 0.9 km level had that of a level
  150 m up.
- **The time step.** The stress, heat and vapour enter the lowest layer implicitly, x′ = (x + K x_s)/(1 + K), which
  stays stable for a layer of any mass.
- **The fluxes.** Over the first month from the same restart the fifteen-level run gives 13.2 W/m² of sensible heat,
  64.0 W/m² of latent heat, a stress of 0.038 Pa and a 10 m wind of 1.24 m/s, against 12.1, 62.7, 0.039 and 1.17 at
  fourteen levels, the 10 m wind on Geleyn's profile from a level a third as high. Over two years with the same
  longwave correction they are 15.6, 62.3, 0.043 and 1.28 against 13.4, 62.9, 0.044 and 1.19.
- **The mixing above.** `vdiff` is ECHAM's implicit K-diffusion with Blackadar's mixing length κz λ/(λ + κz),
  λ = 160 m and heights counted at a fixed 250 K; the new interface at 0.6 km takes its own Richardson number and a
  mixing length of about 90 m. Over the two years the thin layer takes 15.4 W/m² of sensible heat from the ground and
  passes 7.4 W/m² of it up by diffusion and 4.4 W/m² by the dry adjustment, where the fourteen-level bottom layer passed
  2.1 and 3.1 W/m² of its 13.3.

Nothing in the surface layer or the diffusion needs changing for the thin layer. The reading turned up these points,
each the same on every grid:
- **An Earth constant over the sea.** In unstable air over the sea PlaSim raises the heat and vapour transfer by the
  free-convection term of Miller et al. (1992), whose constant, 0.0016 m s⁻¹ K^−1/3, multiplies Δθ_v^1/3. A
  free-convection transfer scales as (g/ρ)^1/3, so at the Moon's gravity and the design air's density the constant
  would be about 0.00085. The ground is warmer than the air in 88–96% of the sea's 3-day means, and on those means
  PlaSim's constant multiplies the transfer by 2.2–2.5 where the lunar one would give 1.6–1.7: it may add up to about a
  third to the sea's evaporation and sensible heat in light winds. The means understate the wind and so overstate the
  term. Scaling it would be a one-constant patch to `surflx` (`zrifh` over the sea) and is the author's decision; this
  work leaves it as it is.
- **Placeholders and floors.** The land's roughness is PlaSim's preset 2 m everywhere. The wind PlaSim feeds the
  transfer has a floor of 1 m/s, under which the lowest level's 3-day mean wind falls in 22% of the land's means at
  0.3 km against 18% at 0.9 km. The diffusion's 160 m mixing length and 250 K height scale are Earth's.
- **Two settings outside the surface layer that follow the grid.** PlaSim's critical humidity for large-scale cloud is
  max(0.85, σ) at each level, so the thin layer's is 0.994, against 0.983 for the fourteen-level bottom layer, and the
  layer above takes 0.977; the bottom layer's cloud covers 0.1% of the sky on both grids, and the runner's
  `surface_cloud_rcrit` sets its threshold. PlaSim's air at 2 m is the lowest level's temperature taken
  dry-adiabatically to the ground's pressure, so it now stands for the air 0.3 km up; with the same longwave correction
  it runs 0.9 K cooler than at fourteen levels over the two years.

**The inputs on fifteen levels**, each made by the script that made the fourteen-level one.
- **The twilight** ([twilight_301k_l15.json](../results/gcm/twilight_301k_l15.json); `twilight levels`). The thirteen
  upper layers keep their fourteen-level terms, the two bottom layers' terms sum to the fourteen-level bottom layer's,
  and the column's totals stay as they were.
- **The longwave correction carried from ten levels**
  ([longwave_A28_corrected_970_lw_ice_l15.json](../results/gcm/longwave_A28_corrected_970_lw_ice_l15.json);
  `plasim_longwave levels`) keeps each split group's escape. On the month trial's columns it leaves 6.4–7.3 W/m² in
  its worst layer, as it did on fourteen levels.
- **The provisional refit**
  ([longwave_A28R_301k_cold_sbm_rh07_deep200_l15cap_month_refit.json](../results/gcm/longwave_A28R_301k_cold_sbm_rh07_deep200_l15cap_month_refit.json);
  `plasim_longwave refit`, the runner's `--longwave-fit l15_provisional`). The correction alone, refitted in 18 s on the
  month trial's 883 cloud-free columns, leaves 1.9–2.1 W/m² on a sample of 107 of the month's columns of every kind. A
  month after the remap is a short sample, so the fit is provisional; the final refit belongs to the configuration
  chosen for a long run. The runner reads a fit made on a run's own level set as it is and refuses it on any other set.
- **The cold cloud's factors** scale the top three layers, which all three grids share.
- **The restart** (`levels.py restart --to l15`; `gcm-runs/A28R_301k_cold_l15_restart/MOST_REST.00019`). The year-19
  restart goes from ten levels to fifteen in one conservative remap. Every column keeps its integrals of water,
  temperature times mass and momentum to within 10⁻¹³ of their size, and the six shared upper layers are copied bit for
  bit.

**The shared source.** Building the fifteen-level executable left the ten- and fourteen-level executables and every
source file byte for byte as they were. The ten- and fourteen-level month trials, run again after the build, repeat
the earlier ones bit for bit in all 61 output variables and their restarts (`bit_check.py`).

**Cost and stability.** Run within three minutes of each other, the month trials took 31.9 s at fifteen levels, 28.7 s
at fourteen and 20.4 s at ten, so a fifteen-level year costs about 1.1 times a fourteen-level one and 1.6 times a
ten-level one. The two-year runs took 448–464 s a year while checks shared the machine part of the time, and wrote
182–183 MB of output a year against 172 at fourteen levels and 130 at ten. Neither the month nor the four years shows a
fault: every field is finite, the humidity is nowhere negative, the strongest wind at any level stays under 27 m/s, and
the temperature falls smoothly through the new levels. All of this work holds 1.3 GB on the research drive.

**Two years at fifteen levels** ([levels_l15_test.json](../results/gcm/levels_l15_test.json);
`levels15_test_analysis.py` in `session-tools/agents/climate-levels15` on the research drive). The baseline,
`A28R_301k_cold_sbm_rh07_deep200_l15cap`, takes every flag of the fourteen-level test and the provisional refit. A
second run, `…_l15cap_carried`, takes the fourteen-level test's own correction carried onto the levels: beside the
fourteen-level test it shows what the thin layer alone changes, and beside the baseline what the refit changes. Both
start from the carried restart. Years 0–1, and after the semicolon year 1:

| Near-side Sea, mm/day | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise | Mean |
|---|---|---|---|---|---|
| 15 levels, provisional refit | 0.12; 0.11 | 1.45; 0.90 | 2.68; 2.85 | 1.98; 1.74 | 1.55; 1.40 |
| 15 levels, carried correction | 0.30; 0.21 | 2.69; 3.03 | 0.87; 1.16 | 1.48; 1.03 | 1.31; 1.34 |
| 14 levels | 0.05; 0.02 | 2.87; 2.57 | 0.88; 0.52 | 0.09; 0.02 | 0.95; 0.77 |
| 10 levels (deep200) | 0.25; 0.19 | 3.66; 3.15 | 3.07; 2.66 | 1.55; 1.46 | 2.11; 1.85 |
| CM1 box | 0.00 | 2.23 | 2.97 | 0.14 | 1.34 |

| Land at 0° E, mm/day | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise | Mean |
|---|---|---|---|---|---|
| 15 levels, provisional refit | 0.50; 0.03 | 17.84; 19.72 | 0.41; 0.54 | 0.03; 0.04 | 4.77; 5.24 |
| 15 levels, carried correction | 1.19; 1.47 | 21.14; 20.04 | 0.99; 1.72 | 0.02; 0.00 | 5.92; 5.96 |
| 14 levels | 0.01; 0.00 | 13.33; 13.05 | 0.54; 0.37 | 0.01; 0.00 | 3.52; 3.46 |
| 10 levels (deep200) | 0.80; 0.83 | 15.79; 17.95 | 0.45; 0.47 | 0.06; 0.05 | 4.34; 4.97 |
| CM1 box | 2.87 | 8.59 | 1.53 | 0.07 | 3.27 |

| Years 0–1 | 15, provisional refit | 15, carried | 14 levels | 10 levels | CM1 box |
|---|---|---|---|---|---|
| Sea site at 41.8; 55.8 km | 258.53; 240.12 K | 259.33; 241.22 | 259.32; 241.20 | 259.90; 241.68 | 251.49; 231.76 |
| Sea site's lowest level: vapour; air | 12.34 g/kg; 299.82 K | 13.43; 301.07 | 12.88; 300.70 | 12.10; 299.39 | |
| Sea site's lowest 0.6 km: vapour; air | 12.34 g/kg; 299.82 K | 13.43; 301.07 | 13.34; 301.40 | 12.41; 300.07 | 11.8 g/kg (lowest 0.5 km) |
| Sea site: air at 2 m; ground | 300.33; 301.39 K | 301.57; 302.29 | 302.15; 302.51 | 300.83; 301.83 | |
| Ground; air at 2 m | 297.82; 297.60 K | 299.03; 299.08 | 299.51; 299.97 | 298.46; 298.50 | |
| Rain; its convective share | 2.18 mm/day; 86% | 2.15; 87% | 2.17; 88% | 2.18; 94% | |
| Imbalance at the top of the air | −3.70 W/m² | +0.63 | +1.40 | +1.45 | |
| Planetary albedo; cloud cover | 0.378; 34.9% | 0.380; 36.4% | 0.378; 36.3% | 0.380; 35.4% | |
| Cloud effect: sunlight; infrared | −39.6; 42.6 W/m² | −40.4; 49.6 | −39.9; 49.7 | −40.5; 50.0 | |
| Top two levels: air | 165.82, 191.56 K | 164.07, 192.64 | 164.04, 192.66 | 164.29, 192.85 | |
| Top two levels: water | 0.0349, 2.33 ppmv | 0.0204, 2.93 | 0.0204, 2.93 | 0.0217, 2.97 | |

The lowest level stands at 0.30 km on fifteen levels and 0.86 km on the others. The lowest 0.6 km on the coarser grids
is the mass-weighted mean over σ 0.9883–1 of the conservative remap's reconstruction of their bottom layer. The vapour
is a mixing ratio, the sea site's numbers are lunar-day means over its cells, and year 1's global numbers and the sea
site's column level by level are in the results file.

**What the thin layer changes** (the carried-correction run against the fourteen-level test).
- **The sea rains after midnight again**, 1.48 mm/day against 0.09 at fourteen levels and CM1's 0.14, while its evening
  stays at 0.87 against CM1's 2.97; its mean, 1.31 mm/day, reaches CM1's 1.34 by that route. On the sea site's 3-day
  means the scheme's parcel, now lifted from 0.3 km, meets thin cloud topping out near 1.2 km in 11–20% of the evening
  and night means, where the fourteen-level grid found it in 2–3%, and convects deep in 3–5% of them against 0–1%
  (`night_tops15.py`).
- **Over land** the afternoon rises from 13.3 to 21.1 mm/day, further from CM1's 8.6, while the morning (1.19 against
  CM1's 2.87) and the evening (0.99 against 1.53) move toward CM1.
- **The near-ground air.** The thin layer over the sea holds 13.4 g/kg at 301.1 K, what the fourteen-level grid's own
  profile implied for its lowest 0.6 km (13.3 g/kg), and 1.6 g/kg more than CM1's 11.8. The coarser layer averaged that
  air faithfully; the excess over CM1 belongs to the GCM's near-ground air itself.
- **Little else moves.** The upper air and the top two levels stay within 0.1 K and 1% of the fourteen-level run's.
  The sea site's air below 10 km runs 0.1–0.35 K cooler, the global ground 0.5 K cooler, and the imbalance at the top of
  the air is 0.8 W/m² smaller.

**What the provisional refit changes** (the baseline against the carried-correction run).
- **The longwave by layer comes close to correlated-k.** On the baseline's own year 1
  ([longwave_check_l15.json](../results/gcm/longwave_check_l15.json)) the refit holds the clear-sky longwave within
  3.0 W/m² of correlated-k by layer on 1,424 cloud-free columns and within 2.6 W/m² on a sample of 1,274 of every
  column, its worst layer at 42 km, against 8.2 and 7.7 W/m² for the carried correction, worst at 19 km, and
  13.3–13.6 W/m² for PlaSim alone. It meets the rule's 2.5 W/m² at every layer but the one at 42 km.
- **The top levels and the balance at the top of the air move.** The refit's top-layer term, which carries the
  column's outgoing longwave in the fit, is 3–4 W/m² smaller than the ten-level fit's. The top level warms 1.7 K and
  holds 70% more water, the second cools 1.1 K and holds 20% less, and the cloud in the top two levels falls from 7.3%
  and 13.6% of the sky to 5.5% and 10.0%. The cloud's infrared effect falls 7 W/m² and the run loses 3.3–4.1 W/m² at the
  top of the air: its ground ends year 1 1.6 K below the carried run's, still cooling about 0.6 K a year.
- **The sea's rain and upper air move.** The sea site's air cools 0.2–1.9 K at every level below the top (up to 2.3 K
  in year 1), its afternoon rain moves into the evening (2.68 mm/day, near CM1's 2.97) and the night (1.98 against
  CM1's 0.14), and its air at 41.8 and 55.8 km cools 0.8–1.1 K, about a tenth of the way to CM1's box.

Each run is one realization of the first two years after the remap, so these show the direction and early size of the
changes.

**What remains before the mass-flux tests.**
- **Their baseline.** Either fifteen-level run can serve. The baseline with the provisional refit has the clear-sky
  longwave within about 3 W/m² by layer and a climate cooling about 0.6 K a year as its top levels' cloud adjusts; the
  run with the carried correction keeps the fourteen-level energy balance and carries layer errors of up to 8 W/m² at
  19 km.
- **The final longwave refit** belongs to the configuration chosen for a long run. Its top-layer term moves the top two
  levels the exobase check reads, so that check belongs beside the refit.
- **The free-convection constant over the sea** is an Earth value on every grid; scaling it to the Moon is the author's
  call.
- **Adopting the fifteen levels** is the author's call.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.levels show l15                     # the set, its heights and layers
climate/gcm/.venv/bin/python -m climate.gcm.levels restart \
    /media/projectspace/terluna-research/gcm-runs/A28R_301k_cold/model/MOST_REST.00019 \
    /media/projectspace/terluna-research/gcm-runs/A28R_301k_cold_l15_restart/MOST_REST.00019 --to l15   # seconds
climate/gcm/.venv/bin/python -m climate.gcm.twilight levels --ts 301 --levels l15             # twilight_301k_l15.json
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave levels --levels l15    # the 301 K refit carried onto them
L=/media/projectspace/terluna-research/session-tools/agents/climate-levels15
$L/run_levels15_tests.sh month15                 # one model month; it builds the fifteen-level executable
$L/run_levels15_tests.sh month14b                # and month10b: the earlier month trials again, for the bit check
climate/gcm/.venv/bin/python $L/bit_check.py A28R_301k_cold_sbm_rh07_deep200_l14cap_month_recheck \
    A28R_301k_cold_sbm_rh07_deep200_l14cap_month A28R_301k_cold_sbm_rh07_deep200_month_recheck \
    A28R_301k_cold_sbm_rh07_deep200_month                                          # every variable and the restarts
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave refit A28R_301k_cold_sbm_rh07_deep200_l15cap_month  # 18 s
$L/run_levels15_tests.sh test15                  # two years with the provisional refit, about 15 minutes on 4 workers
$L/run_levels15_carried.sh                       # two years with the carried correction
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave check --case A28R_301k_cold_sbm_rh07_deep200_l15cap:1 \
    --correction provisional_l15=climate/results/gcm/longwave_A28R_301k_cold_sbm_rh07_deep200_l15cap_month_refit.json \
    --correction carried_301k_l15=climate/results/gcm/longwave_A28_corrected_970_lw_ice_l15.json \
    --output climate/results/gcm/longwave_check_l15.json                          # about a minute
climate/gcm/.venv/bin/python $L/surface_layer_check.py \
    l15_refit_y01=A28R_301k_cold_sbm_rh07_deep200_l15cap:0-1 \
    l15_carried_y01=A28R_301k_cold_sbm_rh07_deep200_l15cap_carried:0-1 \
    l14_y01=A28R_301k_cold_sbm_rh07_deep200_l14cap:0-1 l10_y01=A28R_301k_cold_sbm_rh07_deep200:0-1   # and the months
climate/gcm/.venv/bin/python $L/levels15_test_analysis.py > $L/levels15_test_analysis.txt   # the results file
```

#### The free-convection constant over water and the mass-flux parcel's starting speed (26 October)

The author approved two changes on 26 October, each tested in one two-year run on fifteen levels with the provisional
longwave refit from the fifteen-level restart: PlaSim's free-convection constant over water at its lunar value, and a
starting speed for the mass-flux scheme's surface parcel taken from the boundary layer, which [massflux.md](massflux.md)
records as the third round. Both runs are in [trigger_l15_test.json](../results/gcm/trigger_l15_test.json).

**The constant and its source.** In unstable air over water `surflx` raises the transfer coefficient for heat and
vapour by the interpolation of Miller, Beljaars and Palmer (1992, *J. Climate* 5, 418–434, their eq. 9):
C_H = C_HN (1 + C_R^1.25)^0.8 with C_R = 0.0016 (θ_vs − θ_v1)^1/3 / (C_HN |U|), C_HN the neutral coefficient. In still
air the transfer tends to C_H |U| = 0.0016 Δθ_v^1/3, the free convection of tank experiments, which the paper gives as
its eq. 8: E/ρ = 0.17 (g κ²/(θ_v ν))^1/3 Δθ_v^1/3 Δq, with κ the air's thermal diffusivity and ν its kinematic
viscosity. At 300 K and 1013 hPa that law's coefficient is 0.0017 m s⁻¹ K^−1/3; the authors took 0.0016 so that the
interpolation stays near Monin–Obukhov similarity at low wind speeds. PlaSim applies the same factor to evaporation and
sensible heat over its water cells, the seas and lakes.

**How it scales.** With κ = k/(ρ c_p) and ν = μ/ρ, where the conductivity k and the viscosity μ depend on temperature
alone, κ²/ν = k²/(ρ c_p² μ), and at the ground ρ θ_v = p_s/R. The coefficient therefore goes as (g/p_s)^1/3 at a given
temperature, k²/μ adding about T^0.3. The Moon's gravity alone takes 0.0016 to 0.00088; the design air's surface
pressure, 1215.9 hPa at the water level and 1.2 times Earth's sea-level pressure, lowers κ²/ν by that factor and the
constant by its cube root, and the lunar constant is 0.0016 × 0.517 = 0.00083 m s⁻¹ K^−1/3. The sea's temperature, near
300 K on both worlds, and the air's composition are taken as Earth's. (The fifteen-level section above estimated about
0.00085 from (g/ρ)^1/3.)

**The option.** `--free-convection lunar` sets the constant, now fluxmod's namelist variable `cfree`, to 0.0016 (g p_E /
(g_E p_s))^1/3 from the planet's gravity and the experiment's pressure (`FREE_CONVECTION_EDITS` and
`free_convection_setting` in [exoplasim_run.py](exoplasim_run.py)); without it `cfree` keeps 0.0016, and the run's
configuration records nothing new. The patch went into the shared twilight source on 26 October, changing
`fluxmod.f90`'s hash from acdbd360… to 2f106bc4…. Month trials without the option, run on the patched source, repeat the
earlier months bit for bit in all 61 output variables and their restarts: on ten and fifteen levels with the
Betts–Miller scheme and on fifteen with the mass-flux scheme. Compiled with PlaSim's flags, the namelist variable at its
default gives every transfer factor bit for bit as the literal did
([test_exoplasim_run.py](../tests/test_exoplasim_run.py)). A smoke month with the option keeps every field finite.

**Two years with the lunar constant** (`A28R_301k_cold_sbm_rh07_deep200_l15cap_lunarfc`, the baseline's exact flags with
the option; years 0–1, quarters from sunrise):

| | Baseline | Lunar constant | CM1 sea box |
|---|---|---|---|
| Sea site's evaporation by quarter; mean, mm/day | 3.73, 4.65, 3.53, 3.39; 3.82 | 3.42, 4.74, 3.23, 3.13; 3.62 | 3.10, 2.79, 2.95, 3.08; 2.98 |
| Sea site's latent; sensible heat, W/m² | 110.5; 6.5 | 104.7; 8.1 | 86.2; 6.7 |
| Sea site's lowest 0.6 km: air; vapour | 299.82 K; 12.34 g/kg | 299.81; 12.00 | 299.75; 11.73 |
| Sea site's level at 1.16 km: air; vapour | 298.54 K; 11.68 g/kg | 298.47; 11.57 | |
| Sea site's ground | 301.39 K | 301.96 | |
| Sea site's rain by quarter, mm/day | 0.12, 1.45, 2.68, 1.98 | 0.24, 1.52, 2.55, 1.65 | 0.00, 2.23, 2.97, 0.14 |
| Land site's rain by quarter, mm/day | 0.50, 17.84, 0.41, 0.03 | 0.07, 16.03, 0.43, 0.02 | 2.87, 8.59, 1.53, 0.07 |
| Global evaporation, mm/day; latent and sensible heat, W/m² | 2.172; 62.9, 16.0 | 2.101; 60.8, 15.8 | |
| Global rain, mm/day; its convective share; cloud cover | 2.183; 86.0%; 34.9% | 2.118; 85.0%; 35.1% | |
| Global ground; imbalance at the top of the air | 297.81 K; −3.70 W/m² | 297.91; −3.31 | |
| Top two levels: air; water | 165.82, 191.56 K; 0.0349, 2.34 ppmv | 165.65, 191.36; 0.0335, 2.25 | |

The vapour is a mixing ratio, the sea site's lowest 0.6 km is the fifteen-level bottom layer, and CM1's column is
averaged by mass over the same 0–0.6 km (its lowest 0.5 km holds 11.77 g/kg at 299.83 K). CM1's evaporation is its
surface vapour flux over the sea box's third lunar day, as its own water budget gives it.

**What the constant changes.** It lowers the sea site's evaporation by 0.20 mm/day, 5%, which removes a quarter of the
baseline's 0.84 mm/day excess over CM1's box; three-quarters of the excess remains. The site's lowest 0.6 km dries by
0.34 g/kg, from 0.61 above CM1's to 0.27 above, at the same temperature, and the sea's ground warms 0.57 K as it loses
4 W/m² less to the air in latent and sensible heat. Globally the evaporation and the
rain fall 3%, convection's share of the rain stays at 85–86%, and the imbalance at the top of the air moves 0.4 W/m²
toward balance. The sea's night rain falls from 1.98 to 1.65 mm/day and the land's afternoon from 17.8 to 16.0, both
toward CM1's, while the land's morning falls from 0.50 to 0.07 (CM1 2.87). The top two levels cool 0.2 K and hold 4%
less water. Year 1 alone gives the same evaporation change, 3.74 to 3.54 mm/day at the sea site. The term acts in
light winds: the lowest level's wind over the site averages 3.1–3.2 m/s in both runs.

**The mass-flux parcel's starting speed** ([massflux.md](massflux.md), the third round). With `--massflux-start
boundary_layer` the scheme's surface test parcel starts at 2.2 w*, the boundary layer's convective velocity from
PlaSim's own surface fluxes and the depth of its dry convective layer, the speed of a grid cell's strongest thermals as
CM1's 2 km box of the sea's night measures it (the 99.9th percentile of its vertical wind, 2.1–2.3 w*, with Lenschow and
others' 1980 spread of 0.68 w*). On CM1's own w* that rule gives 2.5–2.6 m/s at the land's onset by day, where CM1
measures 2.6–2.9 m/s, and 1.55 m/s in the sea's night, where CM1 measures 1.1–1.5 m/s with grid spacing from 6 to 2 km.
In PlaSim the speed is 0.7–1.3 m/s over the sea and nil to 0.8 m/s over the wet land, and on the runs' 3-day mean
columns the share of parcels reaching their cloud moves by at most two points. Two years with it
(`A28R_301k_cold_massflux_wstart_l15`) repeat the earlier mass-flux run's global state, convection's share of the rain
42.4% against 86.0% in the baseline, which tripped the stop rule the author set at about 75%. Of the brief's criteria
the run meets the sea's evening alone, 2.99 mm/day against CM1's 2.97. The starting speed hardly decides the trigger:
PlaSim's dry convective layer over the sea is 0.5–4.2 km deep under a condensation level 4–6 km up, and the parcel's
sub-cloud mixing spends any start across the stable air between. A candidate cause of that shallow layer, untested, is
PlaSim's vertical diffusion, which keeps Earth's 160 m mixing length while its coefficient in σ goes as g³. The option
stays in the runner, off by default, and runs without it are bit for bit as before.

```sh
T=/media/projectspace/terluna-research/session-tools/agents/climate-trigger
$T/run_trigger_tests.sh fc_smoke15               # one model month with the lunar constant; it applies the patch
$T/run_trigger_tests.sh fc_month15               # and fc_month10, fc_mfmonth15: the months without it, for the bit check
climate/gcm/.venv/bin/python $T/bit_check.py A28R_301k_cold_sbm_rh07_deep200_l15cap_month_fcsource \
    A28R_301k_cold_sbm_rh07_deep200_l15cap_month_mfsource A28R_301k_cold_massflux_l15_month_fcsource \
    A28R_301k_cold_massflux_l15_month A28R_301k_cold_sbm_rh07_deep200_month_fcsource \
    A28R_301k_cold_sbm_rh07_deep200_month_mfsource                                 # every variable and the restarts
$T/run_trigger_tests.sh fc15                     # two years with the lunar constant, about 16 minutes on 4 workers
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave check --case A28R_301k_cold_massflux_wstart_l15:1 \
    --case A28R_301k_cold_sbm_rh07_deep200_l15cap_lunarfc:1 \
    --correction provisional_l15=climate/results/gcm/longwave_A28R_301k_cold_sbm_rh07_deep200_l15cap_month_refit.json \
    --output climate/results/gcm/longwave_check_trigger_l15.json                  # about a minute
climate/gcm/.venv/bin/python $T/trigger_analysis.py > $T/trigger_analysis.txt   # trigger_l15_test.json
```

#### PlaSim's boundary-layer mixing at lunar gravity (27 October)

The section above named PlaSim's vertical diffusion as a candidate cause of the shallow dry convective layer over the
sea. On 27 October the author agreed that PlaSim's mixing near the ground carries Earth's scales and asked for the
first step of the check proposed with that result: read how PlaSim mixes its lowest air, then compare what the mixing
carries with what the CM1 boxes' boundary layers carry. [pbl_mixing.py](pbl_mixing.py) makes the comparison by quarter
of the lunar day at the two CM1 sites, from the fifteen-level baseline's years 0–1
(`A28R_301k_cold_sbm_rh07_deep200_l15cap`) and the boxes' settled lunar days, and writes
[pbl_mixing_l15.json](../results/gcm/pbl_mixing_l15.json). It reads existing output only.

**How PlaSim mixes its lowest air.** Each time step PlaSim's vertical diffusion (`vdiff`, fluxmod.f90) gives heat,
vapour and momentum the diffusivity K = l² |∂V/∂z| F(Ri). The mixing length is Blackadar's, l = λκz/(λ + κz), with
λ = 160 m for momentum and 160 √7.5 = 438 m for heat and vapour and the heights z counted at 250 K with the planet's
gravity; F is Louis's stability function in ECHAM's form with b = c = d = 5, and the shear takes the difference in wind
speed between levels, at least 1 m/s. The g³ in the coefficient's σ form converts K from heights to σ: ∂σ/∂z =
−gσ/(RT) enters twice for the diffusion and once more for the shear, so gravity enters as the physics requires. Later
in the step, after the convection scheme, the dry adjustment (`mkdca`, rainmod.f90) mixes potential temperature and
vapour evenly through any layers where the dry potential temperature falls with height. Of these settings the two
mixing lengths and the wind floor were set for Earth. The CM1 boxes' boundary-layer scheme (`ipbl` 2) mixes with
K = l² √max(S² − N²/Pr, small), its shear taken column by column from the box's resolved winds, and with an asymptotic
length that [cm1_run.py](../crm/cm1_run.py) stretches by the gravities from 75 m to 453 m.

**What the length does at lunar heights.** PlaSim's diffusivity grows as l² at every stability these columns reach.
ECHAM's unstable function stops depending on l only where its term 75 l² √(−Ri) ζ exceeds 1, ζ being a geometric
factor near 1/(5z²); at lunar heights l/z is small, and that needs instability far beyond any of the 3-day means. Near
the ground κz sets l whatever λ is. With the lengths stretched by the gravities, to 966 m for momentum and 2.65 km for
heat and vapour, the length for heat at the half levels near 0.6, 1.8, 3.2 and 5.2 km above the sea grows from 141,
254, 312 and 350 m to 194, 491, 767 and 1,050 m, and the diffusivity on the site's means grows 1.5–3 times at 0.6–1.8
km, 4–5 times at 3.2 km and 5–9 times at 5.2 km.

**The sea site.** The table reads the mixed layer's depth by one rule in both models, the first height above the GCM's
lowest full level (0.31 km) where the mean virtual potential temperature has risen 0.5 K, and gives vapour as a mixing
ratio at the GCM's full levels near 0.3 and 6.2 km. The heat is what each model carries upward through its half levels:
the GCM's from its own heating by vertical diffusion and dry adjustment, which PlaSim writes as single time steps at the
end of each 3-day mean and which is binned here by the time of those steps, half a mean (18° of the lunar day) after its
midpoint; CM1's from the three parts in which it writes each flux: the resolved part, the advection's implicit
diffusion, which belongs to the resolved advection, and its boundary-layer scheme's subgrid part. CM1's effective
diffusivity is its vapour flux over its mean gradient, where its vapour falls by at least 0.05 g/kg per km.

| Sea site | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise |
|---|---|---|---|---|
| Mixed layer's depth, km: GCM; CM1 | 2.3; 3.2 | 2.9; 2.6 | 5.8; 5.6 | 4.3; 6.5 |
| Vapour at 0.3 km, g/kg: GCM; CM1 | 12.16; 11.18 | 12.52; 11.26 | 12.11; 12.64 | 12.57; 11.44 |
| Vapour at 6.2 km, g/kg: GCM; CM1 | 7.46; 10.54 | 8.47; 10.18 | 9.45; 10.50 | 9.24; 11.09 |
| Surface sensible; latent heat, W/m²: GCM | 8.3; 108 | 4.2; 135 | 4.9; 102 | 8.5; 98 |
| The same, CM1 | 7.1; 90 | 2.0; 81 | 9.2; 85 | 8.3; 89 |
| Heat carried up through 0.6, 3.2, 5.2, 7.2 km, W/m²: GCM | 0.9, 0.1, −0.6, −0.6 | −0.9, 1.5, 0.9, −0.7 | 6.4, 7.8, −0.1, −0.3 | 8.6, 6.8, −0.9, −1.0 |
| The same, CM1 | 4.8, −1.2, −0.9, −1.7 | 0.5, −2.8, −1.0, −1.2 | 8.9, −6.2, −12.2, −7.3 | 10.9, 4.1, −5.8, −9.3 |
| Vapour carried up through 5.2 km by CM1, W/m² (resolved advection) | 0.1 (0.3) | 28.2 (25.4) | 112.3 (93.3) | 28.2 (26.3) |
| K at 5.2 km, m²/s: GCM's median on its means; CM1's scheme; CM1's effective | 0.1; 4; n/a | 0.3; 20; 82 | 0.3; 17; 62 | 0.1; 22; 103 |

Near the ground PlaSim carries heat at about CM1's rate, more than half of it by its dry adjustment: after sunset and
after midnight 6–9 W/m² rise through 0.6 km in the GCM against 9–11 W/m² in CM1. Higher up the two part. After midnight
CM1's mixed layer reaches 6.5 km and the GCM's 4.3 km, and in the morning 3.2 km against 2.3. At every hour CM1 keeps
its vapour nearly uniform up to 6.2 km, losing 0.35–2.1 g/kg over that depth, while the GCM's falls 2.7–4.7 g/kg; the
GCM's lowest level holds 1.0–1.3 g/kg more than CM1's in three quarters and 0.5 g/kg less after sunset, and its air at
6.2 km holds 1.1–3.1 g/kg less in all four. CM1 carries vapour upward through 5.2 km at 28 W/m² in the afternoon and
after midnight and at 112 W/m² after sunset, 83–93% of it by resolved advection, the advection's implicit diffusion
included, and at the top of its mixed layer it carries heat down, 6–12 W/m² at 5–7 km after sunset and after midnight,
drawing the warmer air above into the layer. Its own scheme mixes at 17–22 m²/s at 5.2 km in those quarters. PlaSim's
air there has a median Richardson number of 3–8 and a median diffusivity of 0.1–0.3 m²/s, so it carries almost no heat
down (0.1–1.0 W/m² at 5–7 km), and its dry adjustment acts only where the air is dry-unstable. With the lengths
stretched, the diffusivity at 5.2 km after midnight would be 5 m²/s on the means. (Qualified later on 27 October: the
CM1 box here was forced by the Kuo run; forced by the fifteen-level GCM's own circulation, its mixed layer after
midnight is 3.6–4.5 km deep, [below](#the-fifteen-level-run-over-four-years-and-the-sea-box-it-forces-27-october).)

**The other two runs.** With the mass-flux scheme, which convects at half the baseline's share, the GCM's mixed layer
over the sea is shallower at every hour, 1.7–4.0 km, and its lowest level moister, 14.7–16.5 g/kg: above its mixed
layer the GCM moves vapour mostly by its convection scheme. With the lunar free-convection constant the depths stay
within a third of a kilometre of the baseline's and the lowest level is 0.2–0.6 g/kg drier.

**The land site.** There the comparison measures the two models' different air more than their mixing. CM1's land box
is 1.8–4.0 K warmer in virtual potential temperature and 1.2–5.5 g/kg moister through its lowest 6 km than the GCM's
land cells, which put −5 to +6 W/m² into sensible heat by day against CM1's 27–37 W/m², and 196–229 W/m² into latent
heat against 146–171. By the common rule the GCM's daytime mixed layer over land, 2.4–3.0 km, is as deep as CM1's,
1.6–2.8 km.

**What it suggests.** PlaSim's mixing lengths are Earth's. Stretched by the gravities, the scaling CM1's setup assumes
for its own scheme and one no experiment here has tested, they strengthen PlaSim's diffusion mostly inside its mixed
layer and leave it at a few m²/s in the stable air above, where the two models differ. CM1's boundary layer carries
vapour and heat through its whole depth by eddies that rise from the ground and, at the layer's top, mix in the stable
air above; PlaSim's diffusion is local, and its dry adjustment acts only on dry-unstable air. Transport through the
whole boundary layer, with entrainment at its top, is therefore a candidate for part of the difference, and the
comparison does not size that part. PlaSim's convection scheme also carries vapour above the mixed layer, the run writes
no vapour by process, and the two models differ as well in their convection, circulation, radiation and land surface. In
the afternoon and evening the mixed layers are as deep in both models while their vapour profiles still differ, so depth
alone does not settle the moisture. A scheme that mixes the sub-cloud layer more deeply would also have to keep the
stable layer above 6 km that caps CM1's night clouds, or it could mend the humidity and bring the wrong storms. Two
families of boundary-layer schemes set such transport from the boundary layer's own depth and convective velocity, so
their scales follow the Moon's: nonlocal K-profile schemes with a countergradient term and a depth from a bulk
Richardson number (Holtslag and Boville 1993, as in CCM3 and CAM3), and eddy-diffusivity mass-flux schemes, whose dry
updrafts become the shallow clouds (Siebesma, Soares and Teixeira 2007). Neither paper has been read here.

**Evidence and limits.** The GCM's profiles, winds and surface fluxes are 3-day means of one two-year run in the site
cells, and its heating by process is single time steps at the end of those means, ten a lunar day. The heat comparison
holds PlaSim's diffusion and dry adjustment against CM1's whole transport, its resolved motions included, and leaves out
PlaSim's convection scheme, whose heating mixes transport with condensation. Recomputed on the 3-day means, the
diffusion's heat flux disagrees with the run's own in sign and size within the mixed layer, and the mean diffusivity
exceeds the median severalfold where a minority of the means are unstable, so the diffusivities and Richardson numbers
on the means place each half level's regime and leave its transport to the run's own heating. The run writes no vapour
tendency by process: PlaSim computes the diffusion's and the dry adjustment's (its gridpoint diagnostics 14 and 17), and
the runner leaves them out. CM1's values come from one realization per box on a 6 km grid, whose boundary-layer eddies
the resolution test found about 30% weaker than at 3 and 2 km.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.pbl_mixing --output climate/results/gcm/pbl_mixing_l15.json   # seconds
climate/gcm/.venv/bin/python -m climate.gcm.pbl_mixing --run A28R_301k_cold_massflux_l15 \
    --output /media/projectspace/terluna-research/session-tools/agents/climate-pbl/pbl_mixing_A28R_301k_cold_massflux_l15.json
```

#### A nonlocal boundary layer for the Moon: the single-column test (27 October)

On 27 October the author approved the next step after the comparison above: an option for lunar mixing lengths, a
nonlocal boundary-layer scheme for PlaSim with an explicit entrainment flux, a Python copy tested in a single column
against the CM1 boxes, the scheme as Fortran in `vdiff`, and two-year runs if the column passed a stop rule set in
advance: the new scheme had to close at least half of the column's gap to CM1 in the sea night's mixed depth and in its
0.3–6.2 km vapour drop. It closed neither, so no GCM run was made. The scheme is in [pbl_scheme.py](pbl_scheme.py), the
column test in [pbl_column.py](pbl_column.py) with its results in
[pbl_column_l15.json](../results/gcm/pbl_column_l15.json), and the Fortran in [pbl_fortran.py](pbl_fortran.py).

**The scheme.** It is Holtslag and Boville's nonlocal K-profile (1993) as CAM3 implements it, with the explicit
entrainment flux of Hong, Noh and Dudhia's YSU scheme (2006):
- **The depth h** comes from the bulk Richardson number between the lowest level and the levels above, in Vogelezang
  and Holtslag's form with 100 u*² added to the shear. In unstable air the lowest level's virtual potential temperature
  gains the thermals' excess 8.5 (w′θv′)₀/w_m, and h is the height where the number reaches zero, the inversion's base,
  as YSU sets it when the entrainment flux is explicit; in stable air the critical number is CAM3's 0.3.
- **The velocity scales** come from u* and w*: w_m = (u*³ + 0.6 w*³)^1/3 for momentum and w_m/Pr for heat and vapour,
  the Prandtl number matched to surface-layer similarity at 0.1 h; in stable air u*/φ_h(z/L), Pr = 1.
- **The profile** K = κ w z (1 − z/h)² holds inside h, where the diffusivity is the larger of the profile's and PlaSim's
  local one. Above h PlaSim's local diffusion runs alone, with whichever mixing length the run takes.
- **The countergradient flux** K γ, γ = 7.2 w* (w′c′)₀/(w_m² h), carries heat and vapour.
- **The entrainment flux** at h is −0.15 (θv/g)(w*³ + 5u*³)/h for θv, −0.15 times the surface buoyancy flux in free
  convection. Heat and vapour cross the inversion at the rate it sets, each in proportion to its own jump between the
  levels bracketing h, and the flux falls as (z/h)³ below.
- **The numerics** follow Holtslag and Boville's appendix B: the nonlocal fluxes enter as an explicit step before
  `vdiff`'s implicit solve, so the column keeps the sums of σ-thickness times temperature, vapour and wind exactly as
  `vdiff` keeps them. A vapour profile the explicit step would leave negative keeps its own.

Two choices are the scheme's own, each an interpretation. PlaSim's lowest interface stands 0.6 km up, at or above the
lunar surface layer's top (0.1 h) for any h below 6 km; there CAM3's surface-layer forms would be evaluated near z/L =
−20, beyond the range the Dyer–Hicks forms were fitted on, and their heat velocity u*/φ_h grows without bound as u*
falls. So the scheme takes the outer layer's velocity scales at every interface inside h. And the dry adjustment
(`mkdca`) stays as it is, the brief's default; how often it still acts is reported below. Holtslag and Boville ran CCM2
without dry adjustment below its top three layers.

**Its surface fluxes.** `surflx` runs before `vdiff` in each step and leaves its fluxes in PlaSim's module: the sensible
heat `dshfl` and the evaporation `devap` (downward) and u*³ (`dust3`, from the stress `mkstress` computes). The scheme
divides them by the surface air's density p_s/(R T_s) and, for the sensible heat, by the moist heat capacity `mkshfl`
uses, which recovers the kinematic fluxes PlaSim's own transfer formula gave.

**The constants.** Every constant of the scheme is dimensionless:

| Constant | Value | Source |
|---|---|---|
| Thermals' excess b | 8.5 | HB93 (3.12), (A14): b = 0.85 d, d ≈ 10 from Troen and Mahrt; CAM `fak` |
| Countergradient a | 7.2 | HB93 (3.10), (A14): a = 0.85² d; CAM `fakn` |
| Surface layer's share of h | 0.1 | HB93 (A13); CAM `sffrac` |
| Unstable gradient constants | 15, 15 | Dyer (1974) via HB93 (A6), (A8); CAM `betam`, `betah` |
| Stable gradient constant | 5 | HB93 (A3), (A5); CAM `betas` |
| Critical bulk number, stable; unstable | 0.3; 0 | CAM3 (4.471); YSU (A12) |
| u*² factor in the bulk number | 100 | Vogelezang and Holtslag (1996) via CAM3 (4.471); CAM `fac` |
| Entrainment ratio; its u*³ factor | 0.15; 5 | YSU (A8), (A9), after Noh and others (2003) |

The IFS and CAM's 'HBR' variant take 0.2 for the entrainment ratio. The implementations hold Earth limits that the
scheme leaves out: CAM's floor on h of 700 u* (0.07/f with f at Earth's middle latitudes; the Moon's f is 5 × 10⁻⁶ s⁻¹
at most), its 50 m rule for cloud in the lowest layer, its 400 hPa ceiling on h, its 30 m free-air length, its least
diffusivity of 0.01 m²/s and least u* of 0.01 m/s; WRF's YSU (v4.0) caps K at 1000 m²/s and the entrainment zone at
100 m, takes the inversion's jump as 0.001 K/m times h, caps the countergradient terms at 3 K and 2 g/kg, and sets its
stable critical number over water from Earth's Coriolis parameter, 10⁻⁴ s⁻¹; Hong and Pan (1996) capped the excess at
3 K. The one cap kept is YSU's on the entrainment velocity, at most w_m, the layer's own scale.

**PlaSim's own Earth settings.** `--mixing-length lunar` sets `vdiff_lamm` to 160 m × g_E/g = 966 m (2.65 km for heat
and vapour), the rule [cm1_run.py](../crm/cm1_run.py) applies to CM1's `l_inf`. That rule is a scaling assumption the
CM1 boxes share, with no independent support yet, so the runs here take it as a sensitivity. The wind floor `zumin`, 1
m/s, is a velocity, which a column stretched by the gravities keeps, and it matters: on the CM1 boxes' winds laid on
fifteen levels the difference in wind speed between neighbouring levels from 0.6 to 7.2 km is under 1 m/s at every
output (medians 0.03–0.41 m/s), so the floor sets the shear in `vdiff`'s diffusivity throughout the lunar boundary
layer, and the lowest level's speed falls below it on 23% of the sea box's outputs and 67% of the land box's, where
`surflx` holds the squared wind at 1 m²/s². It is left as it is. `vdiff` also counts its heights at a fixed 250 K.

**The column test.** Each column takes a CM1 box's air on PlaSim's fifteen levels a quarter of the lunar day before the
quarter compared and runs two quarters with the box's radiation, day-night table, mean ascent, surface fluxes,
friction velocity and wind, in PlaSim's order: the surface fluxes and the vertical diffusion, radiation, the mass-flux
convection scheme with PlaSim's evaporation of its rain, the dry adjustment, then the large-scale terms and a saturation
adjustment. The control runs Python copies of PlaSim's `vdiff` and `mkdca`, which the Fortran tests below hold to the
model's own routines. All four boundary layers start from the same air under the same forcing: Earth's or the lunar
length, each with PlaSim's local diffusion alone or the nonlocal scheme beneath it. Means over the compared quarter
(mixing ratio; the mixed depth by the comparison's rule above, the first height where θv has risen 0.5 K above its
value at the lowest level):

| Column (quarter compared) | Measure | CM1 | PlaSim (160 m) | 966 m | 160 m, nonlocal | 966 m, nonlocal |
|---|---|---|---|---|---|---|
| Sea night (from sunset) | Mixed depth, km | 6.48 | 4.87 | 4.57 | 4.46 | 4.41 |
| | Vapour drop 0.3–6.2 km, g/kg | 0.34 | 3.34 | 3.08 | 2.37 | 2.29 |
| | θv rise 6.2–8.2 km, K | 1.38 | 1.96 | 1.86 | 1.59 | 1.58 |
| | Vapour carried up through 3.2, 5.1, 7.1 km, W/m² | 40, 28, 22 | 12, 2, 0 | 14, 5, 0 | 23, 6, 0 | 25, 8, 1 |
| | Rain at the ground, mm/day | 0.16 | 1.06 | 0.97 | 0.77 | 0.73 |
| Sea evening (from noon) | Mixed depth, km | 5.58 | 4.72 | 4.90 | 5.14 | 4.79 |
| | Vapour drop, g/kg | 2.11 | 3.61 | 3.05 | 2.28 | 2.81 |
| Land afternoon (from sunrise) | Mixed depth, km | 2.75 | 6.86 | 6.51 | 3.54 | 3.55 |
| | Vapour drop, g/kg | 3.03 | 5.42 | 5.19 | 4.20 | 4.11 |
| | Vapour carried up through 3.3, 5.2, 7.3 km, W/m² | 181, 231, 242 | 60, 32, 1 | 68, 46, 11 | 116, 75, 47 | 119, 80, 53 |

The column's fluxes are those of the vertical diffusion and the dry adjustment; CM1's are its whole transport,
resolved and parameterized. In the sea night the nonlocal scheme with the lunar length leaves the mixed layer 0.46 km
shallower than PlaSim's own mixing, 29% further from CM1, and closes 35% of the vapour drop's gap; the lunar length
alone leaves it 0.30 km shallower and closes 8% of the drop's gap. The stop rule tripped. The stable layer that caps
CM1's night cloud stays in every column (P4: its rise is 1.6–2.0 K against CM1's 1.4). Over land in the afternoon the
scheme closes 81% of the mixed depth's gap and 55% of the vapour drop's, and the dry adjustment, which acts on
two-thirds of PlaSim's steps there, stops acting.

**What the columns show.** In all four sea-night columns the convection scheme goes deep on every step, its tops 57–58
km up, and on 94–96% of the evening's steps before them, where CM1's night cloud stops at 7–8 km. Its heating of 4–6
km, 5–10 W/m², holds an inversion 3.1–3.2 km up, the nonlocal scheme's median h (CM1's own scheme gives 6.1 km), and
the scheme's entrainment, 3.9 mm/s, does not lift it through the quarter. Over the night quarter the sea evaporates 3.1
mm/day into each column, the large-scale terms remove 3.5–3.9, the column rains 0.7–1.1 and its store falls by 1.4–1.5;
radiation cools it by 185 W/m², the large-scale terms warm it by 75, and it loses 70–81 W/m² of heat. Each budget
balances to within 0.13 mm/day and a few W/m².

Started instead from CM1's air at midnight, a sensitivity, all four columns keep CM1's night layer through the quarter
(6.53–6.58 km; vapour drop 0.00–0.53 g/kg, the nonlocal scheme mixing it to nearly nothing against CM1's 0.34), and the
convection goes deep on 43–57% of steps, as on CM1's own columns (42%, [massflux.md](massflux.md)). So the column loses
CM1's night layer while it runs through the evening under deep convection, whichever boundary layer it has.

The sea's night also needs mixing the scheme leaves out. CM1's night layer is topped by cloud: radiation cools its
7–10 km by 12 W/m² after midnight, and 10–12 W/m² of heat goes down through 6–7 km, at least two-thirds of it by CM1's
own subgrid mixing in the cloud layer and the rest resolved. That flux is two-thirds to four-fifths of the surface's
buoyancy flux, 15 W/m² in heat units, where a layer driven from the ground entrains 0.15–0.2 of it, which points to
entrainment driven mostly from the cloud's top (an interpretation). The IFS carries such a term: a K-profile driven by
the radiative cooling at cloud top (Lock and others 2000) and 0.24 ΔR/(ρc_p) in its entrainment (Cy40r1, section
3.3).

**The Fortran** (`PBL_PATCHES`, over the free-convection patch): fluxmod gains the switch `nlpbl` and the scheme's
constants in its namelist; `vdiff`, with the switch on, keeps the state it diffuses, calls `pblnl` after its local
diffusivities and takes its tendencies against the kept state; `pblcol` is the Python copy step for step. A second
patch (`DIAGNOSTIC_PATCHES`) completes PlaSim's gridpoint diagnostic 16, the vapour each layer gains by convection,
which the Betts–Miller and mass-flux routines preset to minus the tendency so far and never complete; with either option
a run would also write codes 30, 33–37 and 39, the vapour by process and the heating of falling rain's evaporation
(`BUDGET` in [exoplasim_run.py](exoplasim_run.py)). The runner's `--mixing-length lunar` and `--boundary-layer nonlocal`
record their settings and write fluxmod's namelist. Checks ([test_pbl_fortran.py](../tests/test_pbl_fortran.py),
[test_pbl_scheme.py](../tests/test_pbl_scheme.py)), compiled with ExoPlaSim's flags against a stub of PlaSim's module:
- the patched `vdiff` gives the Python copy's tendencies of temperature, vapour and wind to 1e-10 on 24 lunar and
  Earth-like columns (unstable, calm, windy, stable, neutral and nearly dry; fifteen levels, ten, and 40 equal layers),
  with either length and the switch on or off; they agree to 2.3e-12 of each tendency's scale;
- it keeps each column's sums of σ-thickness times temperature, vapour and wind to 1e-12 of their content;
- the Python dry adjustment gives PlaSim's `mkdca` tendencies to 1e-10 on 18 columns in two sets, each set's stable
  column left as it was;
- with the switch off it prints every digit the unpatched `vdiff` prints;
- the patched fluxmod.f90 and rainmod.f90 compile with PlaSim's modules on fifteen levels.

No run has set the options, so the shared twilight source does not carry the boundary-layer patch. The convection's
vapour diagnostic (`DIAGNOSTIC_PATCHES`) entered it on 27 October with the night switch of the diagnostic run below,
where months without either repeat earlier months bit for bit and the run's vapour budget closes; for the boundary-layer
patch the runner-level checks of runs without the options and the smoke months remain.

**Evidence and limits.** The columns take CM1's radiation and surface fluxes as given, so neither answers the column's
own state, and by the compared quarter they run 0.1–2.1 K colder and 0.3–2.8 g/kg moister near the ground than CM1, the
sea night most (1.8–2.1 K and 1.3–2.1 g/kg). CM1's values come from one realization per box on a 6 km grid. The
convection in the columns is the mass-flux scheme, whose deep nights the earlier tests found
([massflux.md](massflux.md)); the GCM's baseline runs the Betts–Miller scheme with PlaSim's shallow scheme, which have
no Python copy here, so the columns test the boundary layer under one convection scheme. The implicit diffusion is taken
over the column's 1800 s step, where PlaSim takes it over twice that.

**Sources.** The PDFs and code are in `session-tools/agents/climate-pbl-scheme/literature` on the research drive.
- **Holtslag and Boville (1993), *J. Climate* 6, 1825–1842.** Read: sections 1–3 and appendices A–B, pp. 1825–1831 and
  1840–1841 (journals.ametsoc.org, sha256 981cef3c…); the comparisons of sections 4–5 seen in part.
- **Collins and others (2004), NCAR/TN-464+STR, the CAM 3.0 description.** Read: section 4.11, pp. 149–155
  (opensky.ucar.edu, sha256 a836d66e…). CAM's `hb_diff.F90` and `pbl_utils.F90` (ESCOMP/CAM, cam_cesm2_2_rel) read for
  the implementation's constants and limits.
- **Hong, Noh and Dudhia (2006), *Mon. Wea. Rev.* 134, 2318–2341.** Read: sections 1–3 and appendices A–B
  (journals.ametsoc.org, sha256 7d0fde5b…); section 4 in part. WRF v4.0's `module_bl_ysu.F` read for its limits.
- **IFS Documentation Cy40r1, Part IV, section 3.3** (the K-profile, top-driven mixing and entrainment). Read, from the
  copy the mass-flux work keeps.
- **Not read:** Troen and Mahrt (1986), Dyer (1974), Vogelezang and Holtslag (1996), Noh and others (2003) and Lock and
  others (2000), taken as the read sources give them; the WRF technical note; Siebesma, Soares and Teixeira (2007), kept
  for the eddy-diffusivity mass-flux scheme the author may take up later.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.pbl_column     # three columns, four boundary layers each; 2 minutes
python3 -m pytest climate/tests/test_pbl_scheme.py climate/tests/test_pbl_fortran.py   # with gfortran, about 20 s
```

#### The single column with CM1's own convection (27 October)

The column test above gave two reasons why the nonlocal scheme left the sea's night layer shallow: the mass-flux scheme
convects deep on every night step and heats 4–6 km, which holds the boundary layer near 3 km, and CM1's night boundary
layer is topped by cloud and mixed from above. On 27 October the author approved a test of both. The same columns run
under two more convection setups ([pbl_column.py](pbl_column.py) `--test convection`, results in
[pbl_column_convection_l15.json](../results/gcm/pbl_column_convection_l15.json)), and
[box_lid.py](../crm/box_lid.py) splits what CM1 carries through the top of its night boundary layer by the way it
carries it ([box_lid_near_side_sea_night.json](../results/crm/box_lid_near_side_sea_night.json)). No GCM run was made
and no Fortran changed. A second test, of CM1's own boundary-layer length, is in
[the CM1 account](../crm/README.md#the-sea-box-with-earths-boundary-layer-length-27-october).

**The setups.** Every column keeps the box's radiation, day-night table, mean ascent and surface fluxes, the starts a
quarter early and the 2×2 of boundary layers of the first test; only its convection changes.
- **CM1's own convection.** In place of the mass-flux scheme and PlaSim's evaporation of its rain, the column takes
  CM1's tendencies of its microphysics and of its resolved transport (box_budget's groups: the box's advection with its
  implicit diffusion, the imposed ascent left out), for potential temperature as temperature at constant pressure and
  for vapour. They are laid on the fifteen levels by layer means in pressure and interpolated linearly between the
  box's 3-hourly outputs, each output's budget being one time step's mean over the box's columns. The column carries
  no condensate, so the microphysics' tendencies stand for the condensate's formation and loss, and the column's
  saturation adjustment rains out whatever its own air still condenses. PlaSim's mixing then stands for CM1's
  boundary-layer scheme alone. CM1's resolved transport includes the box's resolved boundary-layer eddies, so below the
  boundary layer's top the nonlocal scheme counts part of them a second time; the like-for-like comparison is PlaSim's
  local mixing against CM1's subgrid scheme at the lid.
- **The control: shallow convection at night.** The mass-flux scheme runs as in the first test, its deep and mid-level
  branches off while the Sun is down (`SHALLOW_ONLY`, a setting of the Python column alone: the deep test parcels
  depart from no level and the mid-level search stops at once). Its shallow branch, a cloud rooted at the ground and at
  most 200 hPa thick, about 9 km near the lunar ground, can still convect. The land's afternoon runs by day and repeats
  the first test.

Below 16 km the groups' tendencies equal the convergence of CM1's own diagnosed resolved fluxes to within 0.6 W/m² of
heat and 0.9 W/m² of vapour in each band and quarter of the sea box's third lunar day. Above 16 km they hold a further
3–11 W/m² of heating in every quarter, 11 W/m² after midnight, that no flux's convergence accounts for; CM1 takes its
advection at fixed heights, where the mean vertical motion of a column that expands and contracts through the day and
the compressible divergence term both appear (an interpretation). The column takes the groups as they are and reckons
the transport's flux at its lid from the ground up.

**CM1's lid, split.** CM1 writes each vertical flux at its w levels in three parts (domaindiag.F, `scalar_flux`): the
resolved part, w times the scalar's departure from its level's mean at sixth order (`thfr`, `qvfr`); the advection's
implicit diffusion, w times its fifth-order interpolation less the sixth-order one (`thfd`, `qvfd`), which belongs to
the resolved advection; and the subgrid part, −khv times the gradient (`thfs`, `qvfs`), which with `ipbl` 2 and no LES
model is the boundary-layer scheme's. Over the held nights' second day, at the top of their boundary layer (upward
positive; the scheme's part in brackets):

| Held night, second day | 6 km grid | 3 km | 2 km |
|---|---|---|---|
| Heat through 6.45 km, W/m² | −12.3 (−10.1) | −10.9 (−8.7) | −10.5 (−8.3) |
| Heat through 6.99 km | −10.2 (−9.1) | −8.9 (−7.6) | −8.5 (−7.0) |
| Vapour through 6.45 km | 15.5 (11.4) | 14.5 (10.1) | 14.5 (9.8) |
| Vapour through 6.99 km | 15.2 (11.7) | 14.2 (10.1) | 14.3 (9.7) |
| The scheme's khv at 6.45; 6.99 km, m²/s | 10.5; 9.8 | 9.1; 8.2 | 8.7; 7.7 |

The scheme carries 78–89% of the heat and 68–77% of the vapour on every grid, and with the implicit diffusion 87–96%
and 75–83%. The first split of these runs counted the diffusion with the scheme and gave 87–95% and 74–83%. The
resolved part carries 0.4–1.5 W/m² of the heat and 2.6–3.7 W/m² of the vapour. Over the sea box's third lunar day after
midnight the scheme carries 10.1 of 10.7 W/m² of heat down and 14.1 of 22.6 W/m² of vapour up through 6.45–6.99 km,
while inside the layer, at 1.8–5.5 km, resolved motion carries 90–98% of the vapour.

**The columns.** Means over the compared quarter, by convection setup, for the four boundary layers in the order
Earth's length (160 m), the lunar length (966 m), each with the nonlocal scheme:

| Column (CM1: depth; drop; rise) | Convection | Mixed depth, km | Vapour drop 0.3–6.2 km, g/kg | θv rise 6.2–8.2 km, K |
|---|---|---|---|---|
| Sea evening (5.58; 2.11; 1.25) | Mass-flux scheme | 4.72, 4.90, 5.14, 4.79 | 3.61, 3.04, 2.28, 2.81 | 1.52, 1.23, 1.37, 1.39 |
| | CM1's own | 3.30, 3.16, 4.59, 4.61 | 4.75, 4.83, 1.75, 1.86 | 0.94, 0.74, 1.12, 0.93 |
| | Shallow at night | 5.01, 6.58, 7.05, 6.80 | 2.44, 1.42, 0.35, 0.63 | 1.21, 1.24, 1.04, 1.11 |
| Sea night (6.48; 0.34; 1.38) | Mass-flux scheme | 4.87, 4.57, 4.46, 4.41 | 3.34, 3.08, 2.37, 2.29 | 1.96, 1.86, 1.59, 1.57 |
| | CM1's own | 5.28, 5.24, 5.06, 5.15 | 1.06, 0.77, −0.65, −0.65 | 2.89, 2.80, 2.54, 2.56 |
| | Shallow at night | 6.69, 6.47, 6.15, 6.34 | 1.76, 1.62, 0.96, 0.82 | 1.42, 1.36, 0.82, 0.85 |
| Land afternoon (2.75; 3.03; 0.60) | Mass-flux scheme | 6.86, 6.51, 3.54, 3.55 | 5.42, 5.19, 4.20, 4.11 | 0.57, 0.53, 0.32, 0.33 |
| | CM1's own | 1.41, 1.40, 6.18, 6.27 | 2.69, 2.72, −0.32, −0.31 | 1.20, 1.14, 0.58, 0.55 |

What each column carries through its interface nearest the lid, 7.1–7.3 km, beside CM1's own through 6.45–6.99 km
(W/m², upward; heat; vapour):

| Through the lid | Sea evening | Sea night | Land afternoon |
|---|---|---|---|
| CM1: its boundary-layer scheme | −7.2; 19.0 | −10.1; 14.1 | −4.6; 22.5 |
| CM1: all its parts | −9.5; 88.0 | −10.7; 22.6 | −14.4; 244.1 |
| Column with CM1's convection: the prescribed transport | −2.2; 68.8 | 0.2; 6.8 | −5.0; 224.5 |
| PlaSim's local mixing beneath it, 160 m; 966 m | −0.04; 0.08 and −0.52; 0.58 | −0.02; 0.04 and −0.27; 0.57 | −0.05; 0.25 and −0.58; 3.17 |
| The nonlocal scheme beneath it, 160 m; 966 m | −0.04; 0.08 and −0.49; 0.74 | −0.03; 0.06 and −0.30; 0.72 | −10.0; 25.2 and −10.0; 24.4 |
| Column with shallow convection at night, PlaSim's local mixing, 160 m; 966 m | −0.03; 0.03 and −0.37; 0.13 | −0.04; 0.16 and −0.49; 2.20 | −0.13; 1.05 and −1.35; 10.9, by day as the first test |

By the stop rule's threshold, half of the gap to CM1 in both the mixed depth and the vapour drop, measured from the
testbed with PlaSim's own mixing (the mass-flux scheme, 160 m):
- **The sea night.** CM1's own convection closes 25% of the depth's gap and 76% of the drop's; with the lunar length
  23% and 86%. The control closes 87% and 52%, and with the lunar length 99% and 57%: it alone passes, in both
  measures. In it no convection heats the layers at 4.2 and 6.1 km, which took 16.6 W/m² from the scheme in the first
  test. Within each setup the nonlocal scheme fails the threshold: beneath CM1's convection it mixes the vapour until it
  rises with height (−0.65 g/kg) and leaves the layer at 5.1 km, and under the control it closes 56–66% of the drop's
  remaining gap while the depth's remaining 0.21 km grows to 0.33 km with Earth's length and shrinks to 0.14 km with
  the lunar one.
- **The sea evening.** CM1's own convection widens both gaps, to a layer 3.30 km deep with a drop of 4.75 g/kg; the
  control closes 35% and 78%.
- **The land afternoon.** CM1's own convection closes 67% of the depth's gap and 86% of the drop's; the control repeats
  the testbed.

**What the columns show.** Over the sea at night the column's shallow layer comes mostly from the mass-flux scheme's
deep night convection. With it off, PlaSim's own local mixing deepens the layer to 6.5–6.7 km and keeps CM1's lid,
1.36–1.42 K against 1.38, while its vapour still drops 1.6–1.8 g/kg against CM1's 0.34. Most of that mixing is the dry
adjustment's: with Earth's length it carries 89–93% of the heat that goes up through 0.6 and 1.8 km.

The exchange through the lid belongs to CM1's boundary-layer scheme, and PlaSim's local mixing carries almost none of
it. Like for like, beneath CM1's convection after midnight, PlaSim's diffusion and dry adjustment carry 0.02–0.27 W/m²
of heat and 0.04–0.57 W/m² of vapour through 7.1 km, 0.2–4% of what CM1's scheme carries through 6.45–6.99 km; in the
control with the lunar length they carry 0.49 and 2.2 W/m². With CM1's convection prescribed and the scheme's exchange
missing, the lid doubles to 2.8–2.9 K and the layer stops at 5.2–5.3 km, while the drop closes 76–86% of its gap. In
those columns CM1's resolved transport and microphysics heat the layer near 8.2 km, at the top of CM1's night cloud, by
3.9 W/m² and cool the layer near 6.1 km by 1.6 W/m². In CM1 the scheme carries 10 W/m² of heat down through the lid,
cooling the air above it and warming the air below. The reading, an interpretation, is that CM1's resolved motion keeps
its layer's vapour nearly uniform and its scheme's exchange keeps the lid weak. In the evening the column with CM1's
convection also lacks the scheme's 19 W/m² of vapour through the lid and the third of the vapour it carries at 0.6 km
(36.5 of 109 W/m²), and its layer stays at 3.2–3.3 km with the vapour held low.

So the column's night gap splits into two parts, and they rest on different evidence. The convection scheme's deep night
convection accounts for most of the depth's gap, which a GCM run can test. The rest, the lid's exchange, is CM1's own
boundary-layer scheme at a length stretched by the gravities. Run at Earth's length
([the CM1 test](../crm/README.md#the-sea-box-with-earths-boundary-layer-length-27-october)), the sea box keeps its night
layer inside the margins set before the run, 6.15 km deep against 6.48 with a drop of 0.53 g/kg against 0.35, while the
scheme's exchange through the lid falls by more than 90%. CM1's night layer therefore does not rest on the stretched
length, and the column with CM1's tendencies prescribed, which cannot adjust, overstates what the exchange does (an
interpretation).

**Evidence and limits.** The columns take CM1's radiation and surface fluxes as given, as in the first test, and CM1's
convective tendencies are single time steps every 3 hours, interpolated, from one realization per box on a 6 km grid.
CM1's resolved transport carries its resolved boundary-layer eddies, so beneath it the nonlocal scheme's results count
them twice. The water budgets close to round-off once the floor on vapour is counted (0.02 mm/day with CM1's convection
and 0.13 mm/day in the control after midnight), and the heat budgets leave at most 5.4 W/m². The rain of the columns
with CM1's convection, 0.08–1.05 mm/day, is the saturation adjustment's alone. The first test's product regenerates
unchanged from the extended module, apart from that module's hash.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.pbl_column --test convection   # 36 columns; a minute
climate/gcm/.venv/bin/python -m climate.crm.box_lid                         # the held nights' lid, a second
python3 -m pytest climate/tests/test_pbl_column.py climate/tests/test_box_lid.py
```

#### Deep convection held off at night: a diagnostic run (27 October)

The single column above found the mass-flux scheme's deep convection on every night step holding the sea's night mixed
layer at 4.9 km against CM1's 6.5, and PlaSim's own local mixing deepening it to 6.5–6.7 km once that convection was
off. On 27 October the author approved one GCM run to test that in the model itself, as a diagnostic and never a
configuration: the fifteen-level baseline (`A28R_301k_cold_sbm_rh07_deep200_l15cap`) with the Betts–Miller scheme's deep
convection held off wherever the Sun is down. The question was whether the night boundary layer then deepens toward
CM1's and its lowest level dries. It does neither. After midnight the sea site still rains 1.98 mm/day, now as
large-scale condensation that heats the air at 3–16 km by 42 W/m², and its mixed layer stands at 3.17 km, 1.1 km
shallower than the baseline's and half of CM1's 6.48 km, with 0.4 g/kg more vapour in its lowest level. Results are in
[night_switch_diagnostic_l15.json](../results/gcm/night_switch_diagnostic_l15.json).

**The switch** (`--sbm-diagnostic-night`, SBM_NIGHT_PATCHES in [exoplasim_run.py](exoplasim_run.py), off by default).
Wherever the geometric cosine of the Sun's zenith at a grid point is at most 0, a column the Betts–Miller scheme would
convect deep takes PlaSim's shallow scheme instead, with the shallow cap, as a column the inhibition limit or the CAPE
trigger holds back takes it. The cosine is the twilight patch's `gmuz`, which `solang` sets each step before the
convection runs; `mksbm` reads it through `sbmnight`, a routine added to radmod.f90 and called without a module
dependency, so rainmod.f90 still compiles on PlaSim's module alone. The runner records the switch in the run's
configuration as a diagnostic, writes `ndiagnight_bm` into rainmod's namelist, and with the moisture diagnostics writes
the vapour by process (BUDGET), the convection's slot completed by pbl_fortran's `DIAGNOSTIC_PATCHES`.

**The shared source.** The smoke month put both patches into the twilight copy: rainmod.f90's hash went from 42122b3c…
to 58d9cf92… and radmod.f90's from 1eb3996c… to ec58dd8f…. Runs without the switch repeat earlier runs bit for bit: the
fifteen-level Betts–Miller month, the fifteen-level mass-flux month and the ten-level Betts–Miller month, rerun on the
patched source, match the months run before the patch in their restarts and, where those months keep their output, in
all 61 output variables (`bit_check_night.json` in `session-tools/agents/climate-entrainment` on the research drive).
The tests ([test_sbm_night.py](../tests/test_sbm_night.py)), compiled with ExoPlaSim's flags against a stub of PlaSim's
module:
- with the switch off, or on where the Sun is up, the patched routine prints every digit the unpatched one prints on 17
  columns;
- with the Sun down it prints every digit the unpatched routine prints under a CAPE limit no column meets, so every
  column takes the held-back branch, and the columns that convect deep by day go to the shallow scheme;
- the patch applies once over the CAPE trigger's, the depth test's and the cap's patches, the patched rainmod.f90 and
  radmod.f90 compile with PlaSim's modules on fifteen levels, and `sbmnight` hands on the cosine below the horizon too;
- the runner records the switch only when it is set, labelled a diagnostic, and refuses it without the Betts–Miller
  scheme.

**The smoke month** (the run's flags over one model month) keeps every field finite and the humidity nowhere negative,
its strongest wind 26.2 m/s. Its convective rain is nil in the 3-day means that lie wholly in the night and 2.81 mm/day
in the others. Its vapour budget closes: the vapour each process gives each layer sums to the physics' total to 6 × 10⁻⁵
mm/day over the globe and 0.03 mm/day in any column.

**Two years** (`A28R_301k_cold_sbm_rh07_deep200_l15cap_diagnostic_night`) from the fifteen-level restart on four
processes took 437 and 453 s a model year against the baseline's 464 and 448, and wrote 206 MB a year against 182, the
difference the vapour by process, 0.58 GB in all. Years 0–1 at the two CM1 sites, each cell the baseline; the diagnostic
run; CM1 (the mixed depth by pbl_mixing's rule, vapour as a mixing ratio at the full levels near 0.3 and 6.2 km, the lid
the virtual potential temperature's rise from 6.2 to 8.2 km; rain with its convective part in brackets):

| By quarter | Sunrise–noon | Noon–sunset | Sunset–midnight | Midnight–sunrise |
|---|---|---|---|---|
| Sea: mixed depth, km | 2.26; 2.60; 3.25 | 2.86; 3.72; 2.64 | 5.80; 4.91; 5.58 | 4.27; 3.17; 6.48 |
| Sea: vapour at 0.3 km, g/kg | 12.16; 12.34; 11.18 | 12.52; 12.63; 11.26 | 12.11; 12.27; 12.64 | 12.57; 12.96; 11.44 |
| Sea: vapour at 6.2 km, g/kg | 7.46; 8.44; 10.54 | 8.47; 8.94; 10.18 | 9.45; 9.77; 10.49 | 9.24; 10.67; 11.09 |
| Sea: vapour drop, g/kg | 4.70; 3.90; 0.64 | 4.05; 3.69; 1.08 | 2.66; 2.50; 2.15 | 3.33; 2.29; 0.35 |
| Sea: lid, K | 1.12; 1.00; 1.30 | 0.71; 0.80; 0.22 | 1.44; 1.17; 1.25 | 1.62; 1.27; 1.44 |
| Sea: rain, mm/day | 0.12 (0.06); 0.31 (0.22); 0.00 | 1.45 (1.14); 2.51 (2.13); 2.23 | 2.68 (2.51); 2.08 (0.22); 2.97 | 1.98 (1.42); 1.98 (0.02); 0.14 |
| Land: rain, mm/day | 0.50 (0.36); 0.96 (0.76); 2.87 | 17.84 (15.69); 15.18 (12.96); 8.59 | 0.41 (0.40); 0.06 (0.04); 1.53 | 0.03 (0.03); 0.00 (0.00); 0.07 |

| Over the run | Baseline | Diagnostic run | CM1 |
|---|---|---|---|
| Sea site's air at 41.8; 55.8 km, K | 258.53; 240.12 | 258.44; 239.98 | 251.49; 231.76 |
| Ground, global, K | 297.81 | 297.92 | |
| Rain, global, mm/day; convection's share | 2.183; 0.860 | 2.190; 0.799 | |
| Imbalance at the top, W/m² | −3.70 | −3.41 | |
| Top two levels: air, K; water, ppmv | 165.8, 191.6; 0.035, 2.33 | 165.7, 191.7; 0.033, 2.36 | |

**What replaces the night's convection.** Over the sea the night still rains at the baseline's rate, 1.98 mm/day after
midnight, but the convection scheme gives 0.02 of it against 1.42, and large-scale condensation forms the rest where the
grid's layers now saturate. That condensation heats the air at 3–16 km by 42 W/m² after midnight and by 75 W/m² after
sunset, against 15 and 5 W/m² in the baseline, most of it at 5–9 km, with cloud covering 0.17–0.19 of the cells at 3–7
km; the convection's own heating, 34 W/m² after midnight in the baseline and mostly at 16–64 km, falls to nil. With the
latent heat released at 5–9 km the mixed layer stops lower, 3.17 km after midnight against 4.27 and 4.91 km after sunset
against 5.80, and the heat carried up through 1.8 and 3.2 km after midnight falls from 5.5 and 6.8 to 3.8 and 4.8 W/m².
The air at 6.2 km gains 1.4 g/kg after midnight and the lowest level 0.4 g/kg, so the vapour drop shrinks from 3.33 to
2.29 g/kg (CM1 0.35) while the lowest level moves further from CM1's, 12.96 against 11.44 g/kg; the lid weakens from
1.62 to 1.27 K against CM1's 1.44.

By day the sea's afternoon convection rains 2.13 mm/day against 1.14, near CM1's 2.23 in all, and the land's afternoon
15.2 against 17.8 (CM1 8.6). Aloft the sea site's air at 41.8 and 55.8 km cools 0.1 K, and the top two levels change by
−0.2 and +0.1 K. Over the globe the convection's share of the rain falls from 0.86 to 0.80, the ground warms 0.11 K and
the imbalance at the top moves from −3.70 to −3.41 W/m².

So holding off the night's deep convection widens the GCM's gap to CM1's night layer. The column's shallow-night control
deepened its layer to 6.7 km under CM1's large-scale terms, prescribed and unchanged whatever its convection did. In the
GCM the sea keeps evaporating 3.6 mm/day after midnight, and the vapour the night convection rained out now condenses on
the grid at 5–9 km (an interpretation of the budgets above). The GCM's shallow night layer therefore does not come from
its night convection alone: with that convection off, the grid's own condensation holds the layer lower still.

**Evidence and limits.** A diagnostic, in one realization of the first two years of a spin-up from the fifteen-level
restart; 3-day means in the six sea cells and four land cells, and CM1 one realization per box on a 6 km grid. Year 1
alone gives the same answers: after midnight the mixed layer is 3.04 km against 4.26, the lowest level 12.78 against
12.36 g/kg, and the rain 1.80 mm/day (0.00 convective) against 1.74 (1.36). Only the diagnostic run writes the vapour by
process; the heating by process comes from both runs. The longwave check (`plasim_longwave check` on year 1,
[longwave_check_night_l15.json](../results/gcm/longwave_check_night_l15.json)) holds the provisional refit within 2.24
W/m² of correlated-k in its worst layer on the run's columns.

```sh
E=/media/projectspace/terluna-research/session-tools/agents/climate-entrainment
$E/run_night_tests.sh night_smoke15        # one model month with the switch; it applies the patch
$E/run_night_tests.sh night_month15        # and night_mfmonth15, night_month10: the months without it
climate/gcm/.venv/bin/python $E/bit_check.py A28R_301k_cold_sbm_rh07_deep200_l15cap_month_nightsource \
    A28R_301k_cold_sbm_rh07_deep200_l15cap_month_mfsource …                  # every variable and the restarts
$E/run_night_tests.sh night15              # two years, about 15 minutes on 4 workers
climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave check \
    --case A28R_301k_cold_sbm_rh07_deep200_l15cap_diagnostic_night:1 \
    --correction provisional_l15=climate/results/gcm/longwave_A28R_301k_cold_sbm_rh07_deep200_l15cap_month_refit.json \
    --output climate/results/gcm/longwave_check_night_l15.json               # under a minute
climate/gcm/.venv/bin/python $E/night_analysis.py > $E/night_analysis.txt    # night_switch_diagnostic_l15.json
python3 -m pytest climate/tests/test_sbm_night.py                             # with gfortran, about 8 s
```

#### An entrainment that follows the cloud's environment: the column test (27 October)

The mass-flux scheme mixes its updraughts by one fixed profile, 0.8 × 10⁻⁴ m⁻¹ aloft and 2 × 10⁻⁴ m⁻¹ more decaying over
2 km above the cloud base, shallow clouds twice that, whatever air surrounds them. CM1's sea box mixes its cloudy
updraughts at 1.1–1.3 × 10⁻³ m⁻¹ at night by h and qt (1.2–1.5 over qt's resampled range) and 0.2–0.3 × 10⁻³ m⁻¹ in the
afternoon ([the CM1 account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)). On 27 October the
author approved a form in which the scheme's entrainment follows the cloud's environment, tested in the single column
against criteria set in advance, with a Fortran routine and a two-year run only if the column passed. On CM1's own
soundings the form gives CM1's night and afternoon rates by qt to within an eighth and most of their contrast, but in
the running column the sea's night still convects deep and its afternoon and evening lose their rain, so the column
failed and no Fortran or GCM run was made ([massflux_entrainment.py](massflux_entrainment.py),
[massflux_entrainment_l15.json](../results/gcm/massflux_entrainment_l15.json)).

**The form** (massflux.py, `entrainment_form` 'environment'; the fixed profile stays the default). Above the cloud base
the updraught and the test parcels mix ε = ε_bg + ε_env (1.3 − RH) (q_s/q_s,base)³, twice that in shallow clouds, RH
being the surrounding air's relative humidity and q_s, q_s,base its saturation specific humidity at the level and at the
cloud base. The second term is the IFS's entrainment (Cy40r1 (6.7)), the form Bechtold and others (2008) introduced to
make the IFS's convection follow the free troposphere's humidity: drier surroundings mix in faster, and the cube of the
saturation ratio, which mimics a cloud ensemble whose larger members reach higher, makes the rate fall with height. Its
constants are the IFS's: ε_env = 1.75 × 10⁻³ m⁻¹ (1.8 × 10⁻³ in Bechtold and others 2014 and de Rooy and others 2013),
the offset 1.3 and the shallow clouds' factor 2. Four choices are the form's own here, each an interpretation:
- **ε_env per metre is Earth's.** Where CM1 measures no rate the scheme stretches Earth's per-metre rates by the
  gravities, which would give 2.9 × 10⁻⁴ m⁻¹ and a night rate of 0.33 × 10⁻³ m⁻¹, a quarter of CM1's. CM1's lunar clouds
  mix per metre at Earth's rates (its night cloud as BOMEX's 0.8–1.7 × 10⁻³ m⁻¹), so the choice between the two values
  rests on CM1's measured rates, a calibration against the reference process.
- **ε_bg stays** the CM1 boxes' 0.8 × 10⁻⁴ m⁻¹ aloft, which the IFS's vanishing saturation ratio leaves out and which
  holds CM1's diluted upper air.
- **The cloud mixes where it is heavier than its surroundings.** The IFS mixes only where its updraught is buoyant;
  CM1's night cloud mixes 1.1–1.3 × 10⁻³ m⁻¹ while 0.05–0.10 K heavier than its level.
- **The rest of the scheme stays**: the test parcels' sub-cloud mixing (1.12 × 10⁻⁴ m⁻¹), the detrainment ε (1.6 − RH)
  and the mid-level clouds' deep rate (the IFS doubles theirs).

**On CM1's own air.** Each box's 3-D output over its settled lunar day, sampled as box_entrainment samples it: the form
on the quarter's mean sounding at the levels where CM1's cloudy updraughts are sampled, its base their lowest sampled
level, weighted by their mass flux below 16 km as CM1's bulk rate is; deep or shallow by the scheme's 200 hPa rule from
that base to the columns' median cloud top. The plumes rise from the lowest 0.5 km through the same sounding, dropping
their condensate (rates per metre, × 10⁻³):

| Quarter | CM1 by qt; by h | The form (kind) | Plume caps: the form as deep; as shallow; the fixed profile as deep, km | CM1's tops: median, 90th percentile, km |
|---|---|---|---|---|
| Sea night | 1.34; 1.10 | 1.18 (shallow) | 7.0; none; 8.5 | 7.3, 7.9 |
| Sea morning | 0.50; 0.32 | 1.61 (shallow) | 26.5; none; 33.3 | 6.7, 7.3 |
| Sea afternoon | 0.33; 0.22 | 0.37 (deep) | 25.3; none; 36.8 | 22.1, 27.8 |
| Sea evening | 0.31; 0.24 | 1.00 (shallow); 0.50 as deep | none; none; 7.6 | 7.3, 25.8 |
| Land morning | 0.20; 0.12 | 0.53 (deep) | 26.9; none; 31.0 | 25.8, 33.8 |
| Land afternoon | 0.20; 0.13 | 0.44 (deep) | none; none; 33.8 | 25.8, 33.8 |

At the sea the form gives the night 1.18 × 10⁻³ m⁻¹ and the afternoon 0.37 × 10⁻³, a contrast of 3.2 against CM1's 4.1
by qt and 4.9 by h, and its plumes stop where CM1's clouds do: at 7.0 km at night and 25.3 km in the afternoon, where
the fixed profile's reach 8.5 and 36.8 km. The contrast comes from the shallow clouds' factor and the saturation ratio,
which stands at 0.2–0.4 where the afternoon's mass flux runs, 4–8 km above its base; the humidity works the other way,
the night cloud's surroundings at relative humidity 0.96–0.98 against the afternoon's 0.77–0.89. In the other quarters
the form mixes 1.6–3.2 times CM1's rates by qt, and its plume never becomes buoyant in the sea's evening or the land's
afternoon.

**The column test** (pbl_column's testbed on fifteen levels, PlaSim's local mixing at Earth's length, each column
started a quarter early; criteria `CRITERIA`, set before the run from the brief; the sea night's from the testbed's own
gaps):

| Column | Criterion | Fixed profile | The form | CM1 | Met |
|---|---|---|---|---|---|
| Sea night | Convective tops below 10 km on at least 75% of steps | 0% | 6% | tops 7.3 km | no |
| | Mixed depth at least 5.68 km | 4.87 | 6.59 | 6.48 | yes |
| | Vapour drop 0.3–6.2 km at most 1.84 g/kg | 3.34 | 0.84 | 0.34 | yes |
| | Lid at least 0.69 K | 1.96 | 1.13 | 1.38 | yes |
| | Rain at the ground at most 0.5 mm/day | 1.06 | 0.18 | 0.16 | yes |
| Sea afternoon | Rain 1.14–3.41 mm/day | 2.26 | 0.99 | 2.27 | no |
| Sea evening | Rain at least 1.97 mm/day | 3.26 | 0.69 | 2.96 | no |
| Land afternoon | Rain 4.36–13.07 mm/day | 9.37 | 7.55 | 8.71 | yes |
| Both afternoons | Tops above 20 km on half the convecting steps | 99%, 100% | 100%, 100% | | yes |
| Day columns | Air at 41.8; 55.8 km no further from CM1's than the fixed profile's (0.1 K) | 0.0–1.2; 1.8–2.9 K | 0.2–2.7; 2.9–4.5 K | | no, in all six |

At night the form leaves the scheme convecting deep on 94% of the steps, from a base at 5.2 km, the top of the mixed
layer, where the fixed profile's clouds start at 1.8 km, and with a fifth of their mass flux (1.7 against 7.8 g m⁻² s⁻¹
at the base, median), and none of its rain reaches the ground. With the convection's heating lifted above the layer, the
column's night layer comes near CM1's in depth, lid and rain. By day over the sea the stronger mixing near the base
keeps the test parcels from deep clouds: the afternoon convects mid-level from bases near 34 km on 95% of its steps and
its convective rain at the ground falls to 0.03 mm/day against 1.97, the evening's to 0.19 against 2.38, and large-scale
condensation gives most of the little rain that falls. With less deep heating every day column drifts colder aloft,
0.2–4.5 K below CM1's at 41.8–55.8 km against 0.0–2.9 K with the fixed profile. The running column from the GCM's air
(massflux.scm on fifteen levels at the sea, three lunar days from the rh 0.7 branch's site mean, 9.7 and 11.3 K warmer
than CM1's at 41.8 and 55.8 km) cools there to +1.7 and −0.3 K with the form against +4.0 and +3.4 K with the fixed
profile. These are the values regenerated with the humidity's means ([year_files.py](year_files.py)); the start's
change, under 1% of its vapour and 0.002 K at those heights, moved the third day's air by up to 0.7 K and its
quarters' rain by up to 1.2 mm/day, so this single column gives the direction of a change, not its size.

So on CM1's mean soundings the form mixes at CM1's night and afternoon rates and caps its plumes at CM1's cloud tops,
while in the running column, whose lowest level stays moister than CM1's at night (12.4 against 11.4 g/kg) under a
weaker lid and whose surface parcel carries the surface layer's excess, the night's deep convection survives the night's
faster mixing and the sea's daytime convection does not survive the afternoon's (an interpretation of the steps'
records). A scratch sensitivity with the IFS's own rule, mixing only where the updraught is buoyant
(`sensitivity_buoyant_only.py` beside the run scripts), convected deep on 97% of the night's steps and rained 0.99 and
0.89 mm/day in the sea's afternoon and evening, so the failure does not rest on that choice.

**Sources.** The PDFs and their text are in `session-tools/agents/climate-massflux/literature` and
`session-tools/agents/climate-entrainment/literature` on the research drive.
- **IFS Documentation Cy40r1, Part IV, sections 6.3–6.4.** Read: (6.7) the entrainment, (6.8)–(6.9) the detrainment,
  (6.21) the deep test parcel's 0.4 times the rate.
- **Bechtold and others (2008), *QJRMS* 134, 1337–1351.** Read as ECMWF Technical Memorandum 556, section 2.1: the
  turbulent and organised entrainment, ε = c₀F_ε + c₁(q_s − q)/q_s F_ε with F_ε = (q_s/q_s,b)², c₀ and c₁ of order 10⁻⁴
  m⁻¹.
- **Bechtold and others (2014), *J. Atmos. Sci.* 71, 734–753.** Read as Technical Memorandum 705, appendix A: 1.8 × 10⁻³
  m⁻¹, shallow twice deep, "as also supported by CRM data".
- **De Rooy and others (2013), *QJRMS* 139, 1–19.** Read: sections 2, 4.3–4.5 and 5 (KNMI's copy, sha256 d3b6d4fc…): the
  IFS's form as (26) with its 200 hPa factor of two; LES rates that fall as clouds deepen; Gregory's (2001) βB/w² (β
  0.03 deep, 0.06 shallow) and β/z, which gave the ECMWF single column too little sensitivity to humidity; and the
  finding that detrainment varies more than entrainment and sets the mass flux's shape.
- **Not read:** Derbyshire and others (2004), *QJRMS* 130, 3055–3079 (Wiley refused the download), taken as de Rooy and
  others and Bechtold and others (2008) describe it; Gregory (2001), Neggers and others (2002) and Simpson (1971), taken
  as de Rooy and others give them. A buoyancy-based form such as Gregory's gives CM1's night cloud, heavier than its
  surroundings, no mixing, and the radius a cloud-size form (0.2/R) needs is not measured here.

**Evidence and limits.** The CM1 check sets mean soundings against bulk rates from one realization of each box, which
the bulk-plume method gives at about half the directly measured mixing (Romps 2010). The columns take CM1's radiation
and surface fluxes as given and run 1.8–2.1 K colder and 1.3–2.1 g/kg moister near the ground than CM1 by the sea's
night with the fixed profile (the first test's account). The running column from the GCM's air has a mixed boundary
layer in place of PlaSim's diffusion.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.massflux_entrainment   # the CM1 check, the columns, the GCM's air; a minute
python3 -m pytest climate/tests/test_massflux_entrainment.py         # under a second
```

#### The sea's night water by process and height (27 October)

After the diagnostic run above, the author approved an analysis of the sea's night water. The question was what carries
vapour up to the 5–12 km layer where the GCM condenses after dark, where that vapour condenses and how much of it falls,
set against CM1's sea box. No run was made. [night_water.py](night_water.py) reads three GCM runs and the box quarter by
quarter on the fifteen-level grid's layers at the site, with results in
[night_water_l15.json](../results/gcm/night_water_l15.json). After midnight both models lift about the same vapour into
the 5–12 km layer, and they part over what happens to it there. CM1's cloud re-evaporates most of what it condenses.
PlaSim rains it, and the GCM's circulation over the layer changes with where its latent heat is released.

**The budgets.** The GCM runs are the Kuo run on ten levels whose budget forced the box (`A28R_301k_cold_moist`, years
0–5), the fifteen-level Betts–Miller baseline (years 0–1) and the same run with deep convection held off at night (years
0–1), the only one that writes the vapour of every process.
- **Storage and physics.** Each layer's storage over a quarter comes from the 3-day means, fitted with the day's mean
  and four harmonics undamped for their windows, the humidity's means as [year_files.py](year_files.py) reads them
  (below), and the column's from the precipitable water. Its physics comes from the single time steps written at the end
  of each mean, about 60 a quarter in two years and 180 in six. The circulation, sideways and vertical transport with
  the horizontal diffusion, is the storage less the physics as applied.
- **Code 40.** The runner names it `dq_dyn`, and until now described it as the dynamics' transport. It is the physics'
  moistening as the spectral step applies it: spectrald adds the diabatic tendencies, and mkdqtgp writes them back to
  the grid. Its global means over the held-off run's year 1 match the physics' sum (code 31) within 0.15 mm/day in every
  layer and to 0.001 mm/day over the column. [crm/gcm_site.py](../crm/gcm_site.py) reads it rightly.
- **The column and the noise.** A quarter's column physics is exactly its evaporation less its rain, from the means. The
  snapshots' column sum misses that by −0.6 to +1.9 mm/day, and the product keeps each year's or two-year block's values
  as a measure of that noise.
- **CM1.** Its own terms close its vapour budget to 0.3 mm/day in each height range and 0.36 mm/day over the column in
  every quarter, and to 0.03 mm/day after midnight. The groups are the prescribed circulation (the day-night table and
  the imposed mean ascent; the nudging is nil there), the resolved transport, the boundary-layer scheme with the sea's
  evaporation, and the microphysics. The microphysics splits into condensation, evaporation from cloud and from rain,
  deposition and sublimation, whose sum matches CM1's net within a few percent.

**After midnight, the columns** (mm/day):

| Column | Evaporation | Rain (convective) | Circulation | Storage |
|---|---|---|---|---|
| Betts–Miller, fifteen levels | 3.40 | 2.08 (1.53) | −4.47 | −3.15 |
| Deep convection held off at night | 3.60 | 2.20 (0.04) | −5.67 | −4.26 |
| Kuo, ten levels (the box's forcing) | 3.44 | 2.99 (2.19) | −3.31 | −2.85 |
| CM1's sea box | 3.08 | 0.14 | −4.90 | −1.90 |

**After midnight, the layer at 5.2–11.8 km**, with each run's range over its years or blocks where its snapshots enter:

| 5.2–11.8 km | Vapour lifted into it, mm/day | Net condensation, mm/day | Radiation, W/m² | Circulation's heating, W/m² |
|---|---|---|---|---|
| Betts–Miller | not separable from its convection | 0.16 large-scale (0.08–0.24) | −24 | +34 |
| Held off | 0.77 (0.67–0.87): diffusion 0.43, shallow scheme 0.34 | 1.54 (1.36–1.73) | −28 | −9 (−14 to −3) |
| Kuo | not separable from its convection | 0.27 large-scale (0.21–0.33) | −20 | +13 (11–15) |
| CM1 | 0.90 resolved, −0.01 its scheme | 0.21: 1.14 condensed, 0.85 and 0.08 re-evaporated | −22 | +18 |

**The supply is alike.** In the held-off run PlaSim's diffusion and shallow scheme lift 0.77 mm/day into the layer. In
CM1 the resolved motion lifts 0.90 mm/day, and its boundary-layer scheme nothing. The circulation's net vapour there is
−0.91 mm/day (−1.42 to −0.39) in the held-off run and −0.49 in CM1.

**CM1's cloud recycles its water.** It condenses 1.14 mm/day in the layer and re-evaporates 0.85 mm/day from cloud and
0.08 from falling rain at the same heights, so only 0.21 mm/day leaves the layer as falling water. The domain-mean
relative humidity stands at 0.97 and 0.92 at 5.2–9.3 km, with cloud on 5–7% of the domain holding 0.02 kg/m² of cloud
water. Radiation cools the layer by 22 W/m² and the resolved eddies by 5. The prescribed circulation warms it by
18 W/m², 7 of them from the imposed mean sinking above 7.5 km, and the cloud's net latent heat by 6, so the layer
cools by only 4 W/m² over the quarter. In the Kuo run itself the circulation warms the same heights by 11–15 W/m².

**The GCM's circulation answers its latent heat.** In the baseline, deep convection releases its heat above 12 km. The
circulation then warms 5.2–11.8 km by 34 W/m², and large-scale condensation forms only 0.16 mm/day there; the night's
2.08 mm/day of rain is mostly deep convection rooted in the boundary layer.

In the held-off run the same heights condense 1.54 mm/day, which heats them by 45 W/m². Against it stand 28 W/m² of
longwave cooling and 9 W/m² of cooling by the circulation, which now lifts the air there, and the layer warms by 5 W/m²
over the quarter. Only 0.03 mm/day re-evaporates in the layer: mklsp passes all its condensate to mkrain in the same
step, and PlaSim holds no cloud water between steps. Its cloud covers 17–19% of the cells at 3–7 km, against CM1's 5–7%
at 5–9 km. Below 5.2 km every GCM run's circulation cools the air by 24–29 W/m² after midnight, against the box's
14 W/m².

**The box's own column.** Over the lunar day the box's forcing takes away 2.05 mm/day, against 0.79 mm/day for the Kuo
run's own circulation at the same cells. The gap holds in every quarter, 0.7–1.6 mm/day. The box's
lowest air holds 11.3–11.5 g/kg after midnight against the GCM runs' 12.0–13.2.

The forcing applies the GCM's mean vertical wind, rising below 7 km and sinking above 7.5 km, to the box's own air. That
circulation brings in vapour in proportion to how much moister the low air is than the air above, and the box, mixing
its vapour up through 6.5 km, keeps that difference small; on its well-mixed air the ascent also cools less (an
interpretation of the table-and-wind construction in [crm/gcm_site.py](../crm/gcm_site.py)). The box's quiet night
therefore rests partly on its own deeper, drier boundary layer, acting through the forcing. Before midnight the box
rains 2.97 mm/day from deep convection, against the baseline's 2.66, so the gap between the models is the second half of
the night.

**What the budgets point to.** The two models agree on how much vapour reaches the night's cloud layer and differ in
what that layer does with it.
- **CM1** re-evaporates three-quarters of its condensate in place. Its net latent heat stays small enough for the
  prescribed circulation's warming to carry the layer through the night, and its 6.5 km mixed layer keeps the low air
  drier.
- **PlaSim** either convects deep from the boundary layer or, with that held off, condenses on the grid and rains all of
  it. Its circulation then lifts the layer further.

The fix therefore belongs in PlaSim's moist physics of that layer: condensate that persists and can re-evaporate, which
PlaSim lacks, and convection that stays shallow at night, with the deep night layer the two would allow. Neither the
circulation alone nor the vapour that the diffusion and the shallow scheme lift is the cause. These are interpretations
of the budgets; whether PlaSim with cloud water of its own would keep CM1's night is untested. (Qualified later on 27
October: the box here was forced by the Kuo run, and forced by the fifteen-level GCM's own circulation CM1's sea night
rains at 71–82% of the GCM's rate, so the budgets measure PlaSim against a night that came with the Kuo forcing;
[below](#the-fifteen-level-run-over-four-years-and-the-sea-box-it-forces-27-october).)

**Two defects found on the way, corrected the same day.**
- **The runner's description of code 40** (`MOISTENING` in [exoplasim_run.py](exoplasim_run.py)) called it the
  dynamics' transport. It now names the physics' moistening as the spectral step applies it.
- **The mass-flux scheme's GCM columns** (`gcm_site_columns` in [massflux.py](massflux.py)) took code 40, plus the
  surface's evaporation, as the boundary layer's non-convective moistening, and lacked the dynamics' heating. They now
  take the circulation's heating and moistening at each mean's hour from the run's own budget at the site.
  - **Affected:** the closure computed offline on GCM columns, in `on_gcm_columns` and the agents' GCM-column probes.
    It had used the physics' total, convection included, with the evaporation counted twice.
  - **Not affected:** the trigger shares, which test parcels alone, and the GCM runs, whose Fortran takes the model's
    own tendencies.

**The humidity in the year files.** The runs' `hus` is not a mean. PlaSim in low-I/O mode adds each field to its
3-day mean every step. With spectral humidity, though, the field it adds for `hus` is recomputed only at output steps
(plasim.f90, the gridpoint step's `mod(nstep,nafter)` test), so each record holds the humidity at the end of the
interval before. That field is not in the restart either, which leaves each year file's first record with about 2% of
the vapour. The other fields are true means.
- **The evidence.** The column's vapour from a record matches the average of that record's precipitable-water mean and
  the one before (7.6 kg/m² rms over cells and records) better than its own (10.4).
- **The rule for existing runs.** [year_files.py](year_files.py) takes each interval's mean as the average of the record
  and the next, and leaves out each year file's first and last records. That brings the column's vapour within 5.7 kg/m²
  of the precipitable water, with no bias (8.2 against 14.2 in the Kuo run). Here it brings the layers' summed storage
  after midnight to within 0.3 mm/day of the precipitable water's, against 0.5–1.0 with the records taken as means.
- **The readers that take it.** gcm_site (the CM1 boxes' forcing), pbl_mixing, massflux, sbm_offline, box_upper_air,
  plasim_longwave's fitting columns and this module. Their products are regenerated, except the longwave fits and
  checks, which keep the columns they were made on until the next refit.
- **What moved.** The night quarters barely did: the sea's mixed layer after midnight moves from 4.27 to 4.28 km. The
  day quarters moved most. In the morning the mixed layer goes from 2.26 to 2.49 km and the lowest level gains
  0.45 g/kg; in the afternoon the layer goes from 2.86 to 2.52 km and the air at 6.2 km gains 0.44 g/kg.
- **The box's forcing.** Tested on the Kuo run, the record fix alone moved the forcing table's moistening by
  0.004 g/kg/day on average, against values of about 0.15. The held profile's lowest level came out 0.1 g/kg (0.8%)
  drier, and the mean ascent moved by at most 2%.
- **Not yet changed.** The climatology, the cross-domain studies' extracts and the rings' readers still take `hus` as
  it stands. Their long means are low by about 0.8% of the vapour, and their profiles by time of day lag by half an
  interval.

**Evidence and limits.** These are diagnostics of existing runs. The GCM's physics by process are single time steps,
averaged over the six sea cells, so a layer's value carries their noise. The years and blocks give its size: the
circulation's vapour at 5.2–11.8 km after midnight differs by up to 1.03 mm/day between the held-off run's two years,
while its heat terms differ by at most 11 W/m². The GCM's circulation of heat includes the difference between the
grid-point physics and the physics as applied, which gcm_site projects onto T21. The fifteen-level runs are the first
two years of a spin-up. CM1 is one realization of its third lunar day on a 6 km grid; over that day its column dries by
0.37 mm/day, more than half of it above 9 km. The box was forced by the Kuo run, and the fifteen-level runs' own
circulation, which differs from it, has met no box.

```sh
climate/gcm/.venv/bin/python -m climate.gcm.night_water     # three runs and the box; under a minute
python3 -m pytest climate/tests/test_night_water.py          # under a second
```

#### Cloud water for PlaSim's large-scale condensation (27 October)

The night water budget above found PlaSim raining all it condenses in the sea's night cloud layer, where CM1 holds
most of the same water as cloud and re-evaporates it. On 27 October the author approved building cloud water into
PlaSim, tested first with the night switch as a diagnostic against a stop rule set in advance. The scheme is in
[condensate.py](condensate.py), with its Python copy, its Fortran and the fall speed it takes from the lunar fall law,
and its tests in [test_condensate.py](../tests/test_condensate.py).

**The scheme** (`--large-scale-condensate cloud`; PlaSim's `ncwat` 1; condensate.CLOUD_WATER_PATCHES).
- **Held.** mklsp hands what condenses in a layer at or above freezing, by the layer's temperature after condensing,
  to that layer's cloud water. Below freezing the condensate falls at once, as before.
- **Evaporated.** mkcwat, after mklsp and before mkrain, goes down each column. Where the layer's air is unsaturated,
  the cloud water evaporates toward saturation, with mklsp's correction for the latent heat, at most the water there.
  Cloud water in a layer now below freezing goes to the precipitation.
- **Turned to rain.** What remains turns to rain by Khairoutdinov and Kogan's (2000) autoconversion,
  1350 qc^2.47 Nc^-1.79, and accretion, 67 (qc qr)^1.15 (kg/kg/s), as CM1's Morrison scheme takes them, with the CM1
  boxes' 100 droplets per cm³. The rain's mixing ratio qr is the rain formed above over the air's density and its
  fall speed. From the lunar fall law with Marshall and Palmer's intercept, which is Earth's, that speed is
  0.99 m/s × (F / 3 mm/day)^0.20, within 1.4% from 0.3 to 30 mm/day.
- **Left out.** The water is neither advected nor mixed, and radiation keeps PlaSim's diagnosed clouds. The restart
  carries the water (`dqcw`); a run starting from a restart without it starts with none.

**The humidity's true mean** (`--humidity-mean`; `nhusmean` 1; HUS_MEAN_PATCHES) writes as `hus_mean` (code 233) the
gridpoint humidity PlaSim adds every step. Both switches sit in PlaSim's own namelist, read before the restart
(PLASIM_SWITCH_PATCHES in [exoplasim_run.py](exoplasim_run.py), with a check of whether the restart holds the water).

**The checks.**
- **The routines.** mklsp and mkcwat, compiled with ExoPlaSim's flags against a stub of PlaSim's module, give the
  Python copy's tendencies, cloud water and precipitation to 1e-10 on lunar columns. With the switch off, mklsp
  prints every digit the unpatched routine prints. The copy keeps each layer's water and energy to round-off. The
  patches apply once, undo exactly, and compile with PlaSim's modules.
- **The model with the options off.** With the patches in the shared source, four months repeat their references bit
  for bit, in every output variable and in their restarts. They are the fifteen-level Betts–Miller month, the
  fifteen-level mass-flux month, the ten-level Betts–Miller month and the night switch's month with the vapour by
  process, all made on the night patch's source; the configurations differ only in the source's hashes.
- **The smoke month** (the night switch's flags with both options) keeps every field finite, ending at 298.41 K
  against the night switch's month's 298.42 K, with 2.17 mm/day of rain against 2.19. Its `hus_mean` gives the
  column's vapour within 0.001 kg/m² of the precipitable water, its first record included.
  - **Cloud water worldwide.** Condensation puts 0.54 mm/day into cloud; 0.11 re-evaporates and 0.47 turns to rain.
    The cloud holds 0.10 kg/m² on average, 13 kg/m² in the thickest column.
  - **Cost.** The month took 37 s against 32 s.
- **The shared source.** rainmod.f90, plasimmod.f90, plasim.f90, restartmod.f90 and outmod.f90 carry the patches from
  this month on. Their new hashes stop any run begun on the older source from resuming; such a run continues as a new
  one from its restart. The humidity mean's code 233 reaches the output through a variable, so that the other branches'
  runners, whose tests list the codes PlaSim always writes, still pass against the shared source; the smoke month rerun
  that way repeats the first bit for bit, in all 73 variables and its restart.

**The stop rule, set before the two-year run.** The run takes the night switch's two years with the cloud water and
the humidity's mean (`A28R_301k_cold_sbm_rh07_deep200_l15cap_diagnostic_night_cloud`), from the fifteen-level restart.
It is judged against the night switch's run without cloud water and against CM1, over the sea after midnight.
- **Pass needs both of these.**
  - **The rain.** The site's rain after midnight falls to at most half of 2.20 mm/day, that is 1.10 mm/day or less.
  - **The cloud layer.** The water 5.2–11.8 km loses to the precipitation after midnight falls to at most half of
    1.54 mm/day, that is 0.77 or less. It is counted as the water that layer's cloud hands to the precipitation,
    plus what mklsp sends straight down.
- **And all of these.**
  - **The day's rain.** The sea's afternoon rain stays within half to one and a half times CM1's 2.23 mm/day. The
    land's afternoon stays no further from CM1's 8.59 than the night switch's 15.18.
  - **The globe.** The imbalance at the top stays within 1 W/m² of −3.41, the global rain within 5% of
    2.19 mm/day, and the ground within 0.5 K of 297.92 K.
  - **The cloud's size.** The sea's column after midnight holds at most 1.0 kg/m² of cloud water, about eight times
    the condensate CM1's night column holds.
- **Reported beside them.** The circulation's heating of 5.2–11.8 km after midnight, the night's mixed layer and its
  lowest level's vapour, the share of the cloud's water that re-evaporates, and the cloud water worldwide.

**The two years fail the rule.** The run took 450 and 437 s a model year against the night switch's 437 and 453.
Results are in [cloud_water_night_l15.json](../results/gcm/cloud_water_night_l15.json) and, by layer,
[night_water_l15.json](../results/gcm/night_water_l15.json) (`sbm_l15_night_off_cloud`).
- **The rain after midnight** is 2.51 mm/day over the sea against 2.20 without the cloud water; the rule asked for 1.10
  or less. Its two years give 2.24 and 2.77 against 2.63 and 1.76, so the rise lies within the years' spread.
- **The cloud layer** at 5.2–11.8 km hands 2.07 mm/day to the precipitation after midnight against 1.55; the rule asked
  for 0.77. It condenses more than before, 2.31 mm/day into its cloud (2.59 and 2.01 in the two years) against the
  1.54 the run without cloud water rained at once (1.73 and 1.36), and re-evaporates 0.49 of it from cloud, 21%,
  where CM1's night cloud re-evaporates 75%. Before midnight it condenses 2.22 mm/day against 1.42 and re-evaporates
  12%.
- **The circulation** over that layer cools it by 17 W/m² after midnight against 9, and by 32 against 17 before
  midnight, while the condensation's latent heat there rises from 45 to 67 W/m² and from 41 to 64.
- **The day.** The sea's afternoon rains 3.02 mm/day against 2.47 (CM1 2.23), inside the rule's band. The land's
  afternoon rains 17.4 against 14.9, further from CM1's 8.6, which fails the rule too.
- **The globe** stays inside the rule. The imbalance at the top is −3.91 against −3.41 W/m², the rain 2.172 against
  2.191 mm/day, the ground 297.85 against 297.93 K, and cloud cover 0.370 against 0.362.
- **The cloud's size.** The sea's column holds 0.85 kg/m² of cloud water after midnight, 0.50 of it at 5.2–11.8 km.
  That is inside the rule's 1.0 and about six and a half times CM1's night condensate, 0.13 kg/m².

**What it shows.** CM1's night cloud recycles its water at the scale of its eddies, condensing in updraughts and
evaporating where they mix with drier air, while the layer's mean stays just below saturation. A store of cloud water
on the grid's mean evaporates only where the whole cell falls below saturation. The night's cooling, by radiation and
by the ascent the condensation drives, rarely allows that, so the water stays, the layer stays saturated, and its
condensation and latent heat grow with the ascent over it. That is an interpretation of the budgets, which show the
growth but not its cause. Cloud water without a subgrid spread of humidity, which a statistical cloud scheme carries,
does not give the sea's night CM1's balance. The switch stays off. The humidity's true mean is available to every new
run.

**Evidence and limits.** These are two years each, one realization, from the fifteen-level restart; the layer terms
are the night budget's single steps, whose spread the two years give. The constants are CM1's microphysics, not fitted
to the GCM, and the fall speed's intercept is Earth's. The water's lack of transport and its absence from the radiation
are the scheme's choices, so they are untested.

```sh
S=/media/projectspace/terluna-research/session-tools/night_water
$S/cw_runs.sh cw_smoke15                   # one month with both options; it applies the patches
$S/cw_runs.sh cw_month15                   # and cw_mfmonth15, cw_month10, cw_nightmonth15: the months without them
climate/gcm/.venv/bin/python $S/bit_check.py A28R_301k_cold_sbm_rh07_deep200_l15cap_month_cwsource \
    A28R_301k_cold_sbm_rh07_deep200_l15cap_month_nightsource …    # every variable and the restarts
$S/cw_runs.sh cw_night15                   # two years, about 15 minutes on 4 workers
climate/gcm/.venv/bin/python -m climate.gcm.night_water
climate/gcm/.venv/bin/python $S/stop_rule.py climate/results/gcm/night_water_l15.json \
    climate/results/gcm/cloud_water_night_l15.json
python3 -m pytest climate/tests/test_condensate.py     # with gfortran, about 10 s
```

#### The fifteen-level run over four years, and the sea box it forces (27 October)

The comparisons above set PlaSim's sea nights against the CM1 box forced by the Kuo run, whose circulation differs from
the fifteen-level runs'. On 27 October the author chose a box forced by the fifteen-level GCM itself, after running
that GCM two years further to see whether it settles ([the CRM's account](../crm/README.md#the-sea-box-forced-by-the-fifteen-level-gcm-27-october)).

**Four years of the fifteen-level baseline.** The run `A28R_301k_cold_sbm_rh07_deep200_l15cap` was made again from its
start on the shared source as it now stands (`A28R_301k_cold_sbm_rh07_deep200_l15cap_cwsource`, 31 minutes on four
cores). Its first two years repeat the original's in all 61 output fields and both restarts.

| Year | Ground, K | Imbalance at the top, W/m² | Rain, mm/day | Sea site's rain after midnight, mm/day |
|---|---|---|---|---|
| 0 | 298.10 | −4.12 | 2.19 | 2.31 |
| 1 | 297.53 | −3.27 | 2.17 | 1.85 |
| 2 | 297.15 | −2.63 | 2.11 | 2.33 |
| 3 | 296.87 | −1.90 | 2.07 | 2.62 |

- **It settles,** its cooling slowing from 0.57 to 0.38 and 0.28 K a year. A line through the imbalance against the
  ground reaches zero at 295.7 K, and an exponential approach fitted to the ground levels off at 296.3 K with an
  e-folding time of 2.7 years. At the design sunlight that is about 3.5 K below the ten-level Betts–Miller run's
  299.6 K. The provisional longwave fit cools this configuration by 3.7 W/m² at the top, so its climate waits on its
  own refit.
- **The sea's night does not follow the drift.** Years 2 and 3 fall inside 1.4–2.8 mm/day, the first two years' range
  widened by their difference, the rule set before the run, so the box took all four years. The evening's rain rose
  each year, from 2.53 to 3.52 mm/day, and the circulation's heating of 5.2–11.8 km after midnight swung between 13
  and 34 W/m² from year to year, the single steps' noise.

**The box.** Forced by those four years, CM1's sea box rains after midnight 0.50 mm/day over its second lunar day,
1.62 over its third and 1.88 over its fourth, against 0.14 under the Kuo forcing and the GCM's 2.28. By the rule set
before the run, CM1's quiet nights came with the Kuo forcing: under the fifteen-level GCM's own circulation the sea's
night rains at 71–82% of the GCM's rate, its boundary layer 3.6–4.5 km deep against the GCM's 4.3. The two now part in
the evening, where the box rains 1.06–1.20 mm/day against the GCM's 2.96 under 18–26 J/kg of inhibition. On the
box's own air PlaSim's scheme convects deep at a quarter to a third of the afternoon and evening outputs and at
none after midnight. Split by 15° of hour angle, those outputs fall where the box's own storms rain, at 60–120°, and
none in the box's lull from 135° to midnight, where the GCM rains 2.3–3.6 mm/day: through the evening the scheme keeps
the box's timing on the box's air, and the models part where their air differs. The quarters alone were first read as
the scheme convecting through the inhibition that holds the box's storms down, which the split corrects.

**What it qualifies.** The night water budget above, the boundary-layer comparison and the cloud water's test set
PlaSim's nights against the Kuo-forced box. Its 6.5 km night layer belongs to that box's forcing, and the reading that
the gap came from PlaSim keeping no cloud water to re-evaporate compared PlaSim with that box; CM1's night cloud still
re-evaporates 83% of its condensate under the new forcing. Their measurements of PlaSim stand. For the Kuo GCM itself
the comparison was like for like, and its night rain stays unsupported by CM1.

**The late evening (28 October).** The box's lull holds on 3 and 2 km grids as on 6 km, by a rule set before the runs,
so on the box's own air CM1 does not support the GCM's rain from 135° to midnight
([the CRM's account](../crm/README.md#the-late-evening-on-finer-grids-and-where-the-models-part-28-october)). In those hours the GCM's
air is 1.0–1.7 K warmer below 8 km and up to 1 g/kg moister below 4 km, and on its mean column the scheme convects deep
from 5.4 to 55 km after 20 J/kg of inhibition, where on the box's it stays shallow. Through the lunar day the GCM's
sea, within 0.3 K of the box's, evaporates 9–50% more under a wind at 10 m 1.5–2.2 times the box's. Two sources of the
low air's difference were put to the test next: PlaSim's Earth free-convection constant over the sea, a quarter of its
extra evaporation on 26 October, which carries none of the late-evening rain (below); and the box's own wind, which its
3-day nudging holds to 44–67% of the GCM's speed through most of the lunar day (the CRM's account).

```sh
/media/projectspace/terluna-research/session-tools/option3/runs.sh l15_four_years      # 31 minutes on four cores
climate/gcm/.venv/bin/python /media/projectspace/terluna-research/session-tools/option3/settling.py   # the rule
climate/gcm/.venv/bin/python /media/projectspace/terluna-research/session-tools/option3/forcing_audit.py
```

#### The free-convection constant over four years, and the two models' light-wind exchange (28 October)

The late evening's account put the two models' parting in their low air, which the sea's evaporation feeds, and
PlaSim's free-convection constant over water was one candidate for the GCM's extra evaporation. Its two-year test of
26 October (above) lowered the sea site's evaporation by 5%. Two questions were taken up: what the lunar constant does
over the four years the box is forced by, read by a rule set before the run; and how each model treats light winds
over water.

**The run.** `A28R_301k_cold_sbm_rh07_deep200_l15cap_lunarfc_cwsource` takes the 26 October run's flags
(`--free-convection lunar` on the fifteen-level baseline) for four years from the same restart on the shared source as
it now stands, beside the baseline's four (`_cwsource`). Its years 0–1 repeat the 26 October run bit for bit in all 61
output fields and both restarts; the configurations differ only in the hashes of modules patched since, their options
off. It ran in 29 minutes on four cores.

**The rule, set before the run was read.** The sea site's rain over 135° to midnight, the late evening where CM1's box
stays quiet, over the four years: at most 1.95 mm/day, a third below the baseline's 2.93, would mean the Earth constant
carries a material part of that rain; within 0.5 mm/day of 2.93, the spread of the years (2.74–3.08 in the baseline,
2.18–3.12 in the two-year lunar run), would mean it carries none the years can show.

| Four years, from the same restart | Baseline (Earth constant) | Lunar constant |
|---|---|---|
| Globe: ground by year, K | 298.10, 297.53, 297.15, 296.87 | 298.17, 297.65, 297.30, 297.02 |
| Globe: imbalance at the top; rain; cloud cover, four-year means | −2.98 W/m²; 2.137 mm/day; 35.2% | −2.73; 2.087; 35.4% |
| Sea site by quarter from sunrise: the sea's surface, K | 300.77, 301.32, 301.16, 300.84 | 301.43, 301.97, 301.82, 301.48 |
| Evaporation by quarter; daily mean, mm/day | 3.61, 4.54, 3.50, 3.37; 3.75 | 3.34, 4.67, 3.20, 3.12; 3.58 |
| Air at 0.3 km by quarter, K | 298.83, 300.24, 299.66, 298.72 | 298.90, 300.35, 299.67, 298.69 |
| Vapour at 0.3 km by quarter, g/kg | 12.36, 12.11, 11.93, 12.29 | 11.99, 11.62, 11.77, 12.08 |
| Rain by quarter; daily mean, mm/day | 0.15, 1.71, 3.02, 2.18; 1.78 | 0.20, 1.18, 2.40, 1.72; 1.38 |
| Rain at 60–75°, 75–90°, 90–105°, 105–120°, 120–135°, mm/day | 2.47, 4.90, 4.67, 2.73, 1.89 | 1.38, 4.13, 2.28, 1.84, 2.03 |
| Rain at 135–150°, 150–165°, 165–180°, mm/day | 2.82, 2.29, 3.64 | 2.94, 2.86, 2.47 |
| Late evening, 135–180°: four years; by year, mm/day | 2.93; 2.74, 2.90, 3.08, 3.00 | 2.76; 2.18, 3.12, 2.13, 3.54 |

The vapour is a mixing ratio. The site's values are the 3-day means at the six cells the box's forcing comes from,
weighted by the cosine of latitude and binned by the site's hour angle; the air at 0.3 km is laid on the box's heights.

**By the rule the Earth constant carries none of the late-evening rain the years can show:** 2.76 mm/day against 2.93,
with the lunar run's years ranging from 2.13 to 3.54. What the constant does change, the sea's energy balance shapes.
The free-convection term's factor on the transfer over the site's unstable water falls from 1.71–2.14 to 1.35–1.53 on
the 3-day means; the sea warms 0.65 K at every hour, and the air–sea contrast in virtual temperature grows from 2.0–2.9
to 2.8–3.8 K, so the evaporation falls only 5% and rises in the afternoon. The site's rain falls 22%, most through
sunset (at 90–105° from 4.67 to 2.28 mm/day, against CM1's 2.68–3.19) and after midnight (from 2.18 to 1.72), while
the late evening holds. The low air dries by 0.2–0.5 g/kg at the same temperature, and in the late evening it stays
about 1.3 K warmer and 0.9 g/kg moister than the box's. Year by year the globe runs 0.07–0.15 K warmer and rains 2%
less, and it stands 0.25 W/m² nearer balance at the top; the two runs settle alike, and a line through the imbalance
against the ground reaches zero at 295.6 K against 295.7.

**Light winds over water in the two models** (`light_wind.py` beside the runs). In light winds over a warm sea both
models raise the surface's transfer of heat and vapour beyond what the mean wind carries, each by a formula calibrated
on Earth:
- **PlaSim** (`surflx`) multiplies the transfer over unstable water by Miller, Beljaars and Palmer's (1992)
  interpolation toward the free convection of tank experiments, whose constant 0.0016 m s⁻¹ K^−1/3 carries Earth's
  gravity. On the baseline's 3-day means at the site it raises the transfer 1.71–2.14 times, where PlaSim's own
  unstable factor without the term gives 1.10–1.11 and the lunar constant, the tank law scaled by (g/p_s)^1/3, gives
  1.33–1.53.
- **CM1** (`sfclay`, `sfcmodel` 2, unpatched in every Terluna build until today) adds a gust to the wind in its bulk
  formulas: over land Beljaars's (1995) convective velocity of the boundary layer's eddies, (g/T_s · h · B)^1/3 from
  the build's gravity, the boundary layer's depth h and its buoyancy flux B; over water the old MM5 formula, a gust of
  √Δθ_v m/s for Δθ_v in K, which has no gravity in it. In the sea box forced by the fifteen-level GCM that gust is
  1.6–2.1 m/s, against a resolved wind of 1.3–1.7 m/s at 50 m, and raises the wind in the bulk formulas 1.6–1.8 times.
  Scaled to lunar gravity as a buoyancy velocity it would give 1.14–1.18; Beljaars's gust from the box's own fluxes
  and depth, at 1.2 w*, gives 1.10–1.17.

At lunar gravity a kelvin of contrast between sea and air carries a sixth of the buoyancy it carries on Earth, so both
Earth formulas likely overstate the light-wind exchange over the lunar sea, by similar amounts, and neither explains
the gap between the models. Beljaars (1995) showed that in free convection the wind in the surface transfer is the
velocity of the boundary layer's large eddies, so that similarity theory needs no one-third law once that velocity
enters the bulk formula; the COARE ocean-flux algorithm carries the same gust, 1.25 w*, with a boundary layer of 600 m
over Earth's tropical ocean (Fairall and others 2003). On that physics the lunar boost falls below even the tank law's
lunar constant (an interpretation). In CM1's box under the GCM's wind, Beljaars's gust over water lowers the
evaporation from 3.42 to 3.21 mm/day and keeps the late evening quiet, while the night still rains 2.3 mm/day
([the CRM's account](../crm/README.md#the-sea-box-with-its-wind-following-the-gcms-and-the-gust-over-water-at-lunar-gravity-28-october)).

**Evidence and limits.** One four-year run of each, from one restart, both still settling by 0.3–0.5 K a year. The
transfer factors are computed on 3-day means, which smooth the wind the model works on and so overstate the term; the
GCM's 10 m wind from its components' means understates its speed. The lunar constant scales the tank law alone; the
gust's scalings use CM1's domain means. Whether the lunar constant, a Beljaars gust or neither suits PlaSim is the
author's call.

```sh
W=/media/projectspace/terluna-research/session-tools/wind_fc
$W/runs.sh lunarfc_four_years                                                   # 29 minutes on four cores
climate/gcm/.venv/bin/python $W/bit_check.py A28R_301k_cold_sbm_rh07_deep200_l15cap_lunarfc_cwsource \
    A28R_301k_cold_sbm_rh07_deep200_l15cap_lunarfc                               # years 0-1, every field and restart
climate/gcm/.venv/bin/python $W/fc_compare.py                                   # the globe, the site and the rule
climate/gcm/.venv/bin/python $W/light_wind.py A28R_301k_cold_sbm_rh07_deep200_l15cap_cwsource:0-3 \
    A28R_301k_cold_sbm_rh07_deep200_l15cap_lunarfc_cwsource:0-3                  # both models' light-wind factors
```

### Cloud amount: what published models say

From a literature search, with Yang et al. (2014) and Way et al. (2018) checked
against the papers:

- **At this sunlight, slow rotation changes the clouds only a little.** Slowly
  rotating planets that are not tidally locked, near Earth's sunlight, move
  their albedo by −0.03 to +0.04 against the same model's 1-day run in ROCKE-3D
  with Earth's continents (Way et al. 2018; Del Genio et al. 2019) and in ECHAM6
  (Salameh et al. 2018). CAM3 ocean worlds are the exception, +0.11 and 20 K
  cooler at 256 days (Yang et al. 2014). The thick day-side deck needs 1.5 times
  Earth's sunlight or tidal locking.
- **The Moon should be on the no-deck side.** Way et al. place the switch to a
  circulation driven by day and night where the solar day exceeds the air's
  radiative time scale, p c_p T/(g F), 1–2 months for Earth. The Moon's column,
  7.3 times Earth's, gives it roughly 7–15 months against a 29.5-day day, so the
  day–night contrast stays weak, as the runs show.
- **Absolute amounts cannot be borrowed.** Five GCMs give the same Earth-rotation
  ocean world 43–86% cover (Yang et al. 2019), and a cloud-resolving model gave
  27–41% where CAM3 and ExoCAM gave 82–97% on tidally locked planets (Yang et al.
  2023). PlaSim's 20–28% is not by itself a bias.
- **What transfers is the change with rotation:** a shortwave cloud effect 1.0–1.3
  times the fast run's, a longwave one 0.6–0.85 times. Testing ExoPlaSim against
  that needs a 1-day-rotation run of the same Moon. The absolute amount needs a
  model that makes its own clouds: a cloud-resolving model at lunar gravity, or
  ROCKE-3D.

| Run | State |
|---|---|
| A, design case | Done (40 years, above). The uncalibrated run and one with a broken spectrum tail are kept for reference |
| B, 35% water | Done (40 years, above) |
| Cloud-water bracket | Done (15 years each, above): `A_earth_path_clouds`, `A_full_column_clouds`, and `A_patch_check` for the patch's regression test |
| Corrected clouds | Done (15 years, above): `A_corrected_clouds`, both corrections; the design case with PlaSim's cloud amount |
| A28, 28% seas with lakes | Done (30 years from cold, above): the design case with the scenario's water, 303 K |
| Dimmer shield (`A28_dim3`, `A28_dim5`, `A28_dim8`) | 3% and 8% run 8 years (paused); 5% run 25 years and settled near 294 K (above) |
| Sun test (`A28_dim5_sun_test`) | Done (2 years from `A28_dim5`'s year 24): the corrected Sun and tilt with the old sunlight split |
| Corrected design case (`A28_dim5_moon`) | Done (30 years from `A28_dim5`'s year 24, above): the 5% dimmer shield with the Sun, tilt and sunlight split corrected, 294.9 K |
| Thin convective cloud (`A28_dim5_moon_thin_convective`) | Done (15 years from `A28_dim5_moon`'s year 29, 2026-10-01): `convective_day_s` 2,400 s, the convective cloud near its floor; 296.6 K over years 5–14, settling near 297.1 K, warmer alike from the ground to 18 km. A test of whether the GCM's daytime cloud separates it from CM1 ([../crm](../crm/README.md), "The GCM's daytime cloud") |
| More low cloud (`A28_dim5_moon_low_cloud`) | Done (15 years from `A28_dim5_moon`'s year 29, 2026-10-01): `surface_cloud_rcrit` 0.6, the lowest layer's cloud threshold lowered from 0.983; cloud cover 43%, 289.9 K over years 5–14, settling near 289.0 K, cooler about alike from the ground to 18 km. A test of whether CM1's foggy nights explain the gap ([../crm](../crm/README.md), "The GCM's low cloud at night") |
| C, 1.0 atm | Not run: no current decision needs it. The 1-D balance puts 1.0 atm 0.2–1.3 K cooler than 1.2 atm, depending on the clouds, and 1.2 atm is the design pressure. The runner keeps it defined |
| Gases, glow and dimming (`A28_gases`, `A28_gases_glow`, `A28_gases_glow_held`) | Done (20 years each from `A28_dim5_moon`'s year 29, 2026-10-13, above), with the calibration branches `A28_cal_control`, `A28_cal_co2_1022` and `A28_cal_control_patched` (one year each) |
| Eddies every 6.5 hours (`A28_corrected_970_eddies`) | Done (2 years from `A28_corrected_970`'s year 19, 2026-10-17): `writes_per_lunar_day` 110, eleven means to each standard 3-day mean, for the eddies' transport the 3-day means miss at 45° N, 2% of it ([the CRM's account](../crm/README.md#what-else-divides-cm1-and-the-gcm-at-the-crossing)) |
| Heating by process (`A28_corrected_970_heating`) | Done (1 year from `A28_corrected_970`'s year 19 on 2 ranks, 2026-10-17): `--heating-diagnostics`, PlaSim's temperature tendency of each process at every level (codes 21-29). At 47° N its longwave cools the lowest 1.8 km 10.6 W/m² more than the correlated-k model and the next layer 11.4 W/m² more: its water continuum was fitted to the outgoing longwave, which leaves its cooling by height free ([the CRM's account](../crm/README.md#what-else-divides-cm1-and-the-gcm-at-the-crossing)) |
| Moistening by process (`A28R_301k_cold_moist`) | Done (6 years from `A28R_301k_cold`'s year 19 on 4 ranks, 2026-10-21): `--moisture-diagnostics` with `--heating-diagnostics`, PlaSim's physics moistening of each layer as a water flux (code 31), its negative-humidity fix (32) and the two stages of its diabatic spectral step, the physics' moistening as it applies it (40) and the horizontal diffusion after it (41), with omega (`wap`). The postprocessor's omega runs about a hundredth of the continuity of the written winds over the equatorial band, which the band's meridional mass flux matches, and is not used. The CM1 boxes on the design climate take their forcing from this budget ([the CRM's account](../crm/README.md#boxes-on-the-design-climate-the-near-side-sea-and-0-e)) |
| The glow pair on the corrected radiation (`A28_corrected_gases_glow`, `A28_corrected_gases_glow_held`) | Done (20 years each from `A28_corrected_970`'s year 19, 2026-10-17), with 13 October's added dimming on the design scale of 0.970: the held branch 0.61 K below the design case (above). It runs again on the radiation with both corrections |
| The longwave by height, a year (`A28_corrected_970_lwcheck`) | Done (1 year from `A28_corrected_970`'s year 19 on 2 ranks, 2026-10-17): `--longwave-correction`, the patch reproducing its Python copy; the high cloud forms (above) |
| Ice optics, a year (`A28_corrected_970_icecheck`, `A28_corrected_970_lwcheck_rebuilt`) | Done (1 year each from `A28_corrected_970`'s year 19, 2026-10-17): the correction with `--ice-optics` on 4 ranks, and the longwave year again on the rebuilt model without them on 2 ranks, which repeats it exactly (above) |
| The 20-year test of both (`A28_corrected_970_lw_ice`) | Done (20 years from `A28_corrected_970`'s year 19 at 0.970 with `--longwave-correction --ice-optics` and the heating diagnostics, on 4 ranks, 2026-10-17): 306.5 K over years 10–19, settling at 306.8–307.1 K ([above](#where-the-corrected-moon-settles)) |
| The design case at 301 K (`A28R`) | Not run: prepared on 19 October, a cold start on the current lakes for about 45 years on 4 ranks, waiting for the author's go-ahead ([its plan](#the-runs-plan)) |
| D, E | D needs ExoPlaSim's ozone profile from the 1-D results; E's trace gases now enter as their CO2 equivalent (`--co2-ppm`, above) |

## First experiments

| Run | Pressure | Water | Shield (sunlight) | Purpose |
|---|---|---|---|---|
| A | 1.2 atm | 25% | titania stack | Design case |
| B | 1.2 atm | 35% | titania stack | Water sensitivity |
| C | 1.0 atm | 25% | titania stack | Pressure sensitivity |
| D | 1.2 atm | 25% | edge_200nm with 1-D ozone | UV-passing shield |
| E | 1.2 atm | 25% | titania stack, trace-gas forcing | Warming lever |

Each needs a slab ocean (50 m, with 10 and 100 m as bounds), a uniform land
albedo as a labelled placeholder until a surface-cover map exists, and CO2 at
400 ppm. Useful outputs:

- global and regional surface temperature, with day–night and pole–equator
  ranges;
- cloud fraction and cloud radiative effect;
- precipitation and sea ice;
- the circulation;
- stratospheric temperature and its day–night swing.
