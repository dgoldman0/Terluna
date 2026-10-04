# Waves and tides of the Open Moon's seas

This folder holds the wave climate of the Open Moon's seas and the coasts the
waves reach. It covers how the light lunar winds build waves at one sixth of
Earth's gravity, where the waves run, which shores take their energy, and how
they break and run up. The seas' monthly tide comes from the
[geography tide product](../../geography/README.md#the-monthly-tide). SWAN
(Simulating WAves Nearshore) computes the sea states and SWASH the individual
waves at the shore, both at lunar gravity, on the 28% atlas and the corrected
design climate A28_dim5_moon.

## What the seas are like

**Light winds raise big, slow waves.** Over the seas the winds are light. The
GCM's surface stress corresponds to 10 m winds of 2.6 m/s at the median and
5.8 m/s at the 99th percentile over the equatorial seas, close to CM1's 2.15
and 5.84 m/s. At one sixth of Earth's gravity the same wind, fetch and duration
build waves 4–5 times as high as on Earth, with periods five times as long
([first experiments](first_experiments.md)).

**Sea states through a lunar cycle.** Driven by 60 days of the GCM's weather,
with the first lunar cycle spinning the seas up and the second measured:

| | Nearside sea | Smythii–Marginis |
|---|---:|---:|
| Area | 5.29 million km² | 285,000 km² |
| Mean Hs over area and time | 1.39 m | 0.88 m |
| Area mean through the cycle | 1.0–1.8 m | 0.4–1.5 m |
| Share of area and time with Hs ≥ 1 m | 68% | 34% |
| Share of area and time with Hs ≥ 2 m | 17% | 1.3% |
| Largest Hs (mean period) | 4.9 m (19.3 s) | 3.4 m (15.8 s) |

The nearside sea, eighteen times larger, stays rougher than Smythii–Marginis
through the whole cycle. Its mean periods run 8–15 s, and metre-scale seas
peak near 17 s. North of 10° N its waves run steadily toward the south-west;
in the south their direction varies with a northward tendency, so the seas
travel toward the equator from both sides. The height exceeded 5% of the time
is 2.4 m at the median node, 1.7 times its mean, and 4.3 m at most. Swell that
has outrun the local wind carries a third to two thirds of the wave energy
([nearside sea](nearside.md), [Smythii–Marginis](cycle.md)).

**How the waves look.** Lunar waves keep Earth's shapes and steepness (height
4–5% of wavelength) and move at 0.41 times the speed of an Earth wave of the
same length. A typical nearside wave, with an 11.5 s mean period, is 34 m long
and travels at 3.0 m/s; an Earth wave that long has a 4.7 s period and
travels at 7.3 m/s. The largest sea, with a 22 s peak period, has crests
127 m apart moving at 5.7 m/s. Spray takes 2.5 times as long to fall from the
same height.

**Where the energy reaches the shore.** Over the nearside sea the shores the
equatorward seas run onto take the most energy. On average through the cycle
these are western Oceanus Procellarum (595–731 W/m), the southern shore of Mare
Imbrium by Delisle and at the eastern end of Montes Carpatus (583–584 W/m) and
eastern Mare Fecunditatis (575 W/m), against 87 W/m along the median coast
([nearside sea](nearside.md)). At the eastern Smythii headland, computed on
474 m to 118 m grids over native terrain, the western face takes 12.5 times the
eastern face's wave energy through the cycle, and the exposed side switches
when the sea turns ([shore](shore.md), [coastal history](coastal_history.md)).
Followed wave by wave on today's steep slopes, the 98th percentile of run-up is
2.1 m at the exposed western face and 1.3 m on the southern shore (largest
2.4 and 1.7 m in 30-minute records). On 1:20 to 1:100 beaches it is 0.9–1.1 m
([run-up](runup.md)).

**The water level moves with a monthly tide.** The orbit's eccentricity and
the librations move the Earth's tide through each month. The typical monthly
range is 3.7 m over the nearside sea (5.5–6.0 m at Fecunditatis, Nubium and
Humorum), 2.1 m in the South Pole–Aitken sea and 0.55 m in Smythii–Marginis,
with 0.78 m at the headland. Mare Fecunditatis joins the nearside sea through a
narrow strait whose geometry and friction set its fortnightly tide: twice the
equilibrium with bed friction for slow flows, 0.8–1.05 times it with the
friction of the metre-per-second flows the strait would carry. Air pressure
tilts the seas by at most 0.18 m, and wind setup is about 1 cm
([geography](../../geography/README.md#the-monthly-tide), [review](review.md)).

**At a shore the waves and the tide meet.** Followed through the month at four
shores of the nearside sea, with the tide of a real month matched to the GCM's
Sun, the water rises and falls about 3 m on the northern shores and 7.6–7.7 m on
the southern ones, and storms arrive on their own days: southern Imbrium and
southern Nubium peak together at 3.5–4.0 m, while western Procellarum stays
mostly below 1 m for ten days and then holds 2.7–3.4 m for six, with swell up
to 100 m from crest to crest ([a month at four shores](nearside.md#a-month-at-four-shores)).

**The shores are regolith.** The seas flood today's regolith-mantled ground. At
lunar gravity the same waves move grains about six times larger, and fine grains
settle about six times slower. Newly flooded slopes rework toward gentle,
fine-grained shores, and coastal water stays turbid.

## How the picture is built

The calculations form one chain; each study's record holds its checks:

| Study | What it computes | Record |
|---|---|---|
| Source audit | How SWAN 41.51 carries gravity; the lunar source changes | [audit.md](audit.md) |
| First experiments | Earth against Moon over a strip; early growth; the Smythii–Marginis pilot and assumed slopes | [first_experiments.md](first_experiments.md) |
| Recovered winds and coastal resolution | CM1 ring and GCM wind records; eastern Smythii at 7.6–1.9 km | [coastal.md](coastal.md) |
| A weather week | The GCM's surface stress coupled to SWAN | [weather.md](weather.md) |
| Smythii–Marginis through a lunar cycle | Two cycles of waves, the second measured | [cycle.md](cycle.md) |
| The nearside sea through a lunar cycle | The same for half of the Moon's water, with swell, coastal power and a month at four shores with their tide | [nearside.md](nearside.md) |
| Shore exposure | Directional energy at the headland on 118 m terrain | [shore.md](shore.md) |
| Coastal history | Arrivals and durations at the headland through the cycle | [coastal_history.md](coastal_history.md) |
| Breaking and run-up | Individual waves on rock and beach profiles (SWASH) | [runup.md](runup.md) |
| Review | What the chain establishes and where its uncertainty lies | [review.md](review.md) |
| Monthly tide | Every sea's tide, with a dynamic check | [geography](../../geography/README.md#the-monthly-tide) |

The atmosphere enters as the GCM's surface stress, every three hours over 60
days at T21. The coupling inverts that stress through SWAN's own drag law, so
SWAN's friction velocity equals the GCM's, and the GCM's sea roughness follows
Charnock's law at lunar gravity. The first lunar cycle spins the seas up and
the second is measured. Regional and coastal nests pass hourly directional
spectra down to 1.9 km, 474 m and 118 m grids over native LOLA terrain, and
SWASH takes the shoreward spectrum onto one-dimensional profiles. The figures
are in [visualization/waves](../../visualization/waves/).

## The physics at lunar gravity

Under Froude similarity the dimensionless parts of the wave physics carry over
to lunar gravity unchanged: the breaking indices, the whitecapping and the
four-wave interactions. Three parts change with gravity:

- the air–sea momentum transfer, which the GCM's Charnock roughness carries into the forced runs;
- the capillary physics of centimetre waves, below the scales computed here;
- sediment, which has its own calculation to come.

Two source changes make SWAN gravity-consistent. AGROW's 1 Hz growth knee
scales with gravity, and the air density fixed in the code, 1.28 kg/m³, becomes
CM1's 1.404 kg/m³ over the seas. With depth breaking switched off, SWAN 41.51
reads an uninitialised variable, so every run keeps breaking on.

## Evidence

- **The code.** The lunar handling of SWAN and SWASH matches their sources. All
  52 wave tests pass with the executables, including the lunar-radius,
  energy-flux, nesting, group-delay and SWASH gravity controls.
- **Basin grid.** The seas run on the atlas's 1-degree nodes. In the pilot,
  0.5-degree nodes changed local heights by up to 41% and the open sea by under
  0.002%.
- **Timestep.** Between 150 and 75 s, Smythii–Marginis's metre-scale seas change
  by under 1% and seas of a few decimetres by up to 16%. The nearside sea's
  300 s run matches a 150 s rerun of its stormiest three days within 3.6% for
  metre-scale seas (0.12% at the 95th percentile); seas of a few decimetres
  differ by up to 86% at single nodes, and the time above 1 m by at most
  2.6 hours.
- **Spectral band.** Seas under about 0.75 m put more than 1% of their variance
  in the top two frequency intervals in 37% of the coastal spectra; the
  nearside reference spectra pass in 1,864 of 1,888 cases.
- **Forcing.** The T21 GCM's lowest level sits about 885 m above the sea, and
  one lunar cycle supplies each sea's statistics.
- **Coasts.** Between 237 m and 118 m grids, shoreline heights change by 2.8% at
  the 95th percentile and by up to 39% at single nodes. SWASH reproduces an
  irregular input spectrum within 5–6%.
- **Tide.** The ephemeris matches JPL Horizons within 0.05° and 11 km.
  Smythii–Marginis and the South Pole–Aitken sea follow their equilibrium tide
  within 1.4% and 6%, and the nearside sea's monthly lines within 7%. The
  GCM's Sun clock, which places the shore month's tide, reproduces the GCM's own
  sunlight to 0.0001°.

## Next calculations

- Further lunar cycles, for the variability between months. The atmosphere
  costs seconds per month, and the nearside sea about 3.3 hours per two cycles.
- The South Pole–Aitken sea's waves.
- The monthly tide's water level in the coastal and run-up calculations, and
  Mare Fecunditatis's strait at 16 px/deg with nonlinear friction.
- A 0.5-degree basin for Smythii–Marginis's cycle and coastal histories.
- One run with SWAN's Janssen physics, which derives Charnock's roughness from
  the configured gravity, to test the growth and drag physics together.
- A wider frequency band for weak seas.
- Sediment supply and beach evolution on regolith coasts.
- The seas' composition and density, a choice for the author.
- Vertical heat and mixing columns and spectral underwater light, from the
  atmosphere's surface fluxes and the illumination products.

## Builds and storage

SWAN 41.51 and SWASH 12.01 come from their official archives, pinned by SHA-256
and licensed GPL-3.0-or-later by Delft University of Technology. The
repository keeps the patches and hashes, and the builds live on the research
drive:

- `python -m climate.waves.build --download --build-root <dir>` makes the stock
  and AGROW builds. `build_coupled(<dir>, density, optimize=True, openmp=True)`
  in [build.py](build.py) adds the coupled-air OpenMP build used since the
  lunar-cycle studies. The recorded builds are in
  `/media/projectspace/terluna-research/wave-runs/swan_builds/lunar-seas-swan/`,
  and the coupled OpenMP executable's SHA-256 begins `8883fc9e`.
- [swash_build.py](swash_build.py) builds SWASH under
  `research/runs/waves/swash_source/`.
- `research/runs/waves` links to `/media/projectspace/terluna-research/wave-runs`,
  which holds the inputs, atmospheric archives and every run. Completed runs are
  reused once their input, executable and output hashes match.

## Reproduction

Each study's record ends with its commands. Run them in this order: the builds,
then [first experiments](first_experiments.md#reproduction),
[coastal](coastal.md#reproduction), [weather](weather.md#reproduction),
[cycle](cycle.md#reproduction), [nearside](nearside.md#reproduction),
[shore](shore.md#reproduction-and-storage),
[coastal history](coastal_history.md#reproduction-and-storage) and
[run-up](runup.md#reproduction-and-storage). The review's numbers come from
`python -m climate.waves.review_checks`, the tide from `python -m geography.tides`
and the shore month from `python -m climate.waves.shore_month`.
