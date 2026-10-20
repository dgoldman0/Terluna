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
| [cloud_water.py](cloud_water.py) | Recomputes from a year of output the cloud water PlaSim's radiation sees, layer by layer, as PlaSim has it and with the runner's corrections |
| [plasim_longwave.py](plasim_longwave.py) | PlaSim's clear-sky longwave column by column, checked against its own heating, and the correction by height fitted to the correlated-k model that the runner's `--longwave-correction` adds (`../results/gcm/longwave_<run>.json`); its `check` sets the correction against correlated-k on states from 295 to 306.6 K (`../results/gcm/longwave_check_301k.json`) |
| [plasim_ice_optics.py](plasim_ice_optics.py) | The ice optics the runner's `--ice-optics` gives cloud colder than −38 °C, beside PlaSim's liquid-water optics (`../results/gcm/ice_optics.json`) |
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
- **The correction's upper part stands.** On the ring's cloud-free columns RRTMG and the correlated-k model agree on
  the clear-sky cooling at 81–107 km, 0.18 and 0.17 K a day; on the GCM's own PlaSim cools that layer 0.04 K a day where
  correlated-k cools it 0.19. Its sunlight there was fitted by layer on 14 October (above).
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
