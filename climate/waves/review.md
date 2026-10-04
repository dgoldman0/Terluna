# Review of the wave studies, 3 October 2026

Codex wrote the wave studies in this folder on 2 and 3 October 2026. This
review checked them against the SWAN and SWASH sources, the 60-day atmospheric
archive and the 28% atlas. [review_checks.py](review_checks.py) recomputes every
number quoted here into [results/review_checks.json](results/review_checks.json).
On 3 October the author asked for the review to be locked in and chose the next
work: the nearside sea's waves, a tide product in geography and a synthesis of
the wave studies ([decisions register](../../research/decisions.md), Water and
climate).

## What the wave chain establishes

The lunar handling of SWAN 41.51 matches its source code. Gravity reaches the
dispersion relation, Komen wind input and whitecapping, and the four-wave
interactions. `COORDINATES SPHERICAL` reads the lunar radius. A source change
replaces the air density fixed in the code (1.28 kg/m³) with CM1's 1.404 kg/m³
over the seas, and another scales the AGROW frequency knee with gravity. With
depth breaking switched off, SWAN 41.51 reads an uninitialised variable
(`KTETA`). The studies avoid that defect by running with breaking on.

The runs forced by the GCM carry lunar air–sea momentum transfer. The coupling
inverts the GCM's surface stress through SWAN's own Wu drag law, so SWAN's
friction velocity equals the GCM's. ExoPlaSim sets open-sea roughness by
Charnock's law with the planet's gravity, `z0 = max(0.018 u*²/g, 1.5×10⁻⁵ m)`.
Across 203,389 open-sea samples of the archive, the stored roughness matches
that law at lunar gravity: the ratio is 1.000000 at the median and between
0.999999 and 1.000001 from the 5th to the 95th percentile. At Earth's gravity
the same law gives roughness 6.04 times smaller. Wu's Earth calibration
therefore enters only the first strip experiment.

Over the equatorial seas the two atmospheric models agree to about 20%. The
GCM's stress corresponds to neutral 10 m winds of 2.6, 4.4 and 5.8 m/s at the
50th, 90th and 99th percentiles; CM1's equatorial ring has 2.15, 4.00 and
5.84 m/s over its water. Their disagreement about the air near the ground
concerns the land.

With the SWAN and SWASH executables, all 52 wave tests pass. They include the
lunar-radius, energy-flux, nesting, group-delay and SWASH gravity controls. The
continuous cycle's headline recomputes from its product: at the reference point
Hs is at least 1 m for 309.6 hours, with a longest stretch of 98.1 hours.

## What the studies leave out

**The seas beyond Smythii–Marginis.** Smythii–Marginis covers 285,000 km², 2.7%
of the Moon's standing water. The nearside sea holds half of it: 5.29 million km²,
with the Earth-facing coasts from Procellarum to Fecunditatis. The South
Pole–Aitken sea holds 22%, 2.32 million km². The 60-day global stress archive
already covers both. Smythii–Marginis comes close to full development under its
strongest weather: its peak basin-mean stress of 0.036 Pa corresponds to a
4.47 m/s wind, whose fully developed sea (`g Hs/U² = 0.24`) has Hs 2.95 m, and
the basin reached 3.38 m. The larger seas' waves therefore depend on what
Smythii–Marginis cannot show: their own winds, longer events, and swell crossing
thousands of kilometres.

**The monthly tide.** The atlas's equipotential already holds Earth's static
tide. The orbit's eccentricity and the optical librations move that tide
through each month, and every coastal and run-up calculation holds the water
level fixed. The [geography tide product](../../geography/README.md#the-monthly-tide)
computes it for every sea, with each sea's volume held, the water's
self-attraction and a dynamic check. Typical monthly ranges:

| Sea | Median over the sea's area | 95th percentile | Largest in 19 years |
|---|---:|---:|---:|
| Nearside | 3.7 m | 6.0 m | 9.8 m |
| South Pole–Aitken | 2.1 m | 3.5 m | 6.7 m |
| Smythii–Marginis | 0.55 m | 1.4 m | 2.8 m |

At the eastern Smythii headland the typical monthly range is 0.78 m. That
moves the waterline 40–80 m on the 1:50 and 1:100 beaches, against run-up of
1–2.4 m. The review's first estimate (Keplerian distance and the leading
librations, [review_checks.py](review_checks.py)) gave 3.4, 2.1 and 0.50 m. Air
pressure tilts the seas by at most 0.18 m (nearside) and 0.09 m
(Smythii–Marginis), 6.0 times the Earth response per hectopascal, small because
the GCM's pressure differences across a sea are small. Wind setup across
Smythii–Marginis is about 1 cm at its peak stress.

**Coasts of regolith.** The shore profiles interpolate 118 m LOLA pixels over
today's regolith-covered surface, with mean slopes of 1:3.8 to 1:6.5. The
stress that starts a grain moving scales with gravity, so the same waves move
grains about six times larger. Fine grains settle about six times slower, so
coastal water stays turbid. Newly flooded slopes will be reworked toward
gentle, fine-grained shores, which makes the beach profiles the central case
and the steep rock profiles an end-member. Sediment also sets underwater light
for the biosphere.

**What the waves look like.** The strongest arrival at the western face has Hs
1.8 m and a 12.2 s mean period. Those waves are 39 m long and travel at
3.2 m/s; an Earth wave of the same length has a 5.0 s period and travels at
7.8 m/s. The basin's largest sea, Hs 3.4 m with an 18.3 s peak period, has
86 m crests moving at 4.7 m/s. Lunar seas keep Earth's shapes and steepness
(Hs/wavelength 0.04–0.05) and move at 0.41 times the speed of an Earth sea of
the same wavelength. Spray takes 2.5 times as long to fall from the same height.

## Where the uncertainty lies

The largest uncertainties sit upstream of the coast:

- **Basin grid.** In the pilot, refining the basin from 1° to 0.5° changed
  local heights by up to 41%. The continuous cycle and every coastal history use
  the 1° grid.
- **Forcing.** The forcing is a T21 GCM whose lowest level sits about 885 m above
  the sea (σ 0.9833, scale height 52.5 km). A 10–20% difference in wind between
  the models becomes roughly 20–45% in the height of a developed sea.
- **One cycle.** The statistics come from one lunar cycle. Extending the
  atmospheric record costs seconds per month.
- **The strip experiment's drag.** Charnock's law at lunar gravity gives drag
  0.85, 1.11 and 1.34 times Wu's at the strip's three winds (2.15, 4.00 and
  5.84 m/s), and friction velocity 0.92, 1.05 and 1.16 times Wu's. Wave height
  grows between the first and second power of friction velocity, from
  fetch-limited growth to full development. The strip's lunar heights are
  therefore 8–15% high at the median wind, 5–11% low at the 90th-percentile
  wind and 16–34% low at the 99th.
- **Weak seas.** These need a wider frequency band: 37% of the coastal spectra
  and 129 of 708 cycle spectra put more than 1% of their variance in the top
  two frequency intervals.

The numerical refinements listed as next steps at the headland (5–6% SWASH
input errors, 5% beach-grid sensitivity, iterations per step) are second-order
beside these.

Under Froude similarity, the dimensionless coefficients transfer to lunar
gravity unchanged: the breaking indices, whitecapping and the four-wave
coefficients. Three parts of the physics change with gravity:

- air–sea momentum transfer, through Charnock's parameter and the GCM's surface layer;
- the capillary physics of centimetre waves;
- sediment.

## Writing and figures

The write-ups state each result next to its numerical checks, and close many
paragraphs on evidence disclaimers ("remains open" or "further work" 21 times).
The README follows the studies in the order they ran. The synthesis README
replaces it with what the seas are like, followed by the evidence.

The figures are clean, and their renderers check the displayed arrays against
the products. They need three things:

- **Geography in the basin maps.** These are bare 1° cells, with no coastline at
  the atlas's resolution, no relief, no names, no wave direction and no locator.
- **A scale that shows the third phase.** The three-phase shore map shares one
  linear 0–700 W/m scale, which leaves the third phase (the eastern face at
  52.9 W/m) nearly black.
- **The two missing figures.** `waves.png` and `wind_coverage.png` are linked
  and were never rendered.

## Smaller items

- **Seawater density.** `shared/scenarios/waves.json` makes 1,025 kg/m³ seawater
  a shared scenario. The seas' composition is undecided: their water is
  imported, and the Moon is poor in chlorine. The waves change by 2.5% between
  fresh and sea water.
- **SWAN builds in `/tmp`.** The builds behind the cycle and coastal runs live
  in `/tmp/terluna-lunar-seas-swan`, which a reboot clears. SWASH lives on the
  research drive.
- **Full maps in JSON.** The committed results hold 3.6 MB of JSON, mostly full
  maps (`pilot.json` 1.7 MB).

## Open, with recommendations

The author has chosen none of these yet:

- **Janssen sensitivity.** SWAN's `GEN3 JANSSEN` derives the same Charnock
  roughness from the configured gravity (`SWIND4` in `swancom3.ftn`). One
  Janssen run would test drag and growth physics together.
- **A 0.5° basin parent.** This would apply to the continuous cycle and the
  coastal histories, before further refinement at the headland.
- **Further lunar cycles.** These give variability between months.
- **A wider frequency band.** Needed for weak seas.
- **The seas' composition and density.** For the author to set.
- **Sediment.** Sediment supply and beach evolution on regolith coasts.
- **Builds on the research drive.** Move the SWAN builds there, beside SWASH.
- **Basin figures.** Draw them over relief and coastline, with names and wave
  direction.
