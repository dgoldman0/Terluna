# The Open Moon's aerosol

5 October 2026. What the Open Moon's seas, forests, plains and polar lands would put into
the air near the ground, built from the project's own designs and climate runs and from the
literature, and what those particles do to the air's electrical conductivity and to cloud
droplets. It serves the [atmospheric-electricity study](../atmospheric_electricity/README.md),
whose conductivity column carried a single assumed aerosol (100 particles per cm³ of
0.05 µm radius, κ 0.3), and the storm runs, which assume 100 cloud droplets per cm³.

`python -m research.studies.open_moon_aerosol.run` (a few seconds) writes
[results/open_moon_aerosol.json](results/open_moon_aerosol.json) (schema
`terluna.research.open-moon-aerosol/1`). It reads the GCM climatology
(`python -m climate.gcm.climatology A28_dim5_moon:20-29`, a product kept out of Git),
geography's atlas and drainage, CM1's ring summaries and the conductivity column.
[sources.json](sources.json) lists the sources. The four literature notes behind every
parameter, with their quotes, derivations and gaps, are on the data drive
(`atmospheric-electricity/aerosol_marine.md`, `aerosol_biogenic.md`, `aerosol_dust_fog.md`,
`aerosol_budget.md`).

This is a structured estimate. Each region's particles are a measured Earth analogue changed
by lunar factors, in clean, central and loaded cases that run from the fewest particles to
the most. Nothing in it is a measurement of the Open Moon.

## The regions and their weather

The regions follow the climate runs and the ecology register's landscapes
([lunar_cycle_ecology](../lunar_cycle_ecology/README.md), F1–F3). Areas come from the GCM
(rain) and geography (seas and lakes); the weather near the ground comes from CM1's rings,
which resolve the air the GCM's lowest layer, 2.5 km deep in the lunar air, averages over.

| Region | Share of the Moon | Day: humidity, rain, wind, mixed layer | Night: humidity, fog, wind, mixed layer |
|---|---|---|---|
| Seas | 28.0 % | 82 %, 0.9 mm/day, 2.2 m/s, 1.7 km | 89 %, none, 1.7 m/s, 1.6 km |
| Lakes | 10.3 % | 83 %, 1.5 mm/day, 2.5 m/s, 2.2 km | 89 %, none, 1.9 m/s, 2.2 km |
| Wet land (forests and wet margins) | 35.6 % | 76 %, 3.0 mm/day, 2.1 m/s, 5.0 km | 97 %, 9 %, 0.6 m/s, 400 m |
| Fog desert (plains, 30–60°) | 16.4 % | 68 %, 0.45 mm/day, 2.0 m/s, 5.1 km | 93 %, 6 %, 0.55 m/s, 280 m |
| Polar dry land | 9.6 % | 74 %, none, calm, 660 m | 95 %, 8 %, calm, 70 m |

The land poleward of 30° gets under 0.5 mm of rain a day under humid air, with dew most
nights: a fog desert like the Namib, with no dry-air desert on low ground.

## What the designs put into the air

**The seas.** The ecology register's living seas, plankton and seaweeds chosen to release
sulfur and iodine, give marine air by day like the remote ocean's: particles formed aloft
from their sulfur and mixed down (aerosol_marine.md Sect. 7). How much sulfur depends on the
seawater's dimethyl sulfide, a design choice (1, 3 or 10 nM). The calm lunar winds raise
little spray: 0.7–5 particles per cm³, with whitecaps set by lunar gravity through Titan's
roughness relation. Seas rich in magnesium and carbonate make spray whose κ depends on the
sea's composition, which hydrology has not yet set. Seaweed coasts make iodine particle
bursts in daylight at low tide, 10⁴–10⁶ per cm³ near the shore (O'Dowd et al. 2002), local
and left out of the area means. Through the night particle formation stops and the marine air
thins toward the Antarctic winter's.

**The forests and wet margins.** Over clean forest the number of particles comes from
particles formed high in the troposphere and brought down by storms (Wang et al. 2016); what
the forest emits sets their mass, growth and κ (aerosol_biogenic.md Sect. 7.0). The lunar
forests' cooler air (21–27 °C by day) cuts their isoprene to about 0.4 of an Amazon day's,
and they make about half its secondary organic aerosol; their light is about the Amazon's.
Isoprene protects leaves from heat, which the lunar climate barely asks of them, so how much
the engineered trees emit is a design choice. The Open Moon's low ionization takes away most
of the particle formation ions drive, a third of cloud nuclei on Earth (Kirkby et al. 2016;
Gordon et al. 2017). Particles live 30 days or more against Earth's week, but mix through a
layer three to six times deeper. The ten-day warm, humid nights of the night food web, with
their decomposers, add fungal nanoparticle bursts (Lawler et al. 2020) and spores; the dusk
fruit fall's fermentation vapours make no particles. Seabird-like colonies and grazers
release ammonia, which makes particles only with acids the clean air lacks.

**The plains of the fog desert.** The seas' sulfur reaches the plains at 1–8 parts per
trillion of sulfur dioxide (derived here), against the 0.7 ppb under which savanna plains
form particles almost daily (Vakkari et al. 2013): the plains are clean, like the forest,
with less vapour to grow their particles. Their dust comes from developed soils (clays,
carbonates, some salt; κ about 0.03), lifted in storm outflows by day: lunar gravity lowers
the wind needed to start soil moving to 0.6 of Earth's and dust settles six times slower,
but dew and fog wet the ground most nights and crusts raise the threshold about fivefold
(aerosol_dust_fog.md Sect. 8). The closed lakes' magnesium and calcium chloride bitterns stay
damp above 28–33 % humidity, which the Open Moon's air never falls below, so the salt flats
raise no dust.

**The polar dry land.** Calm and rainless in CM1's polar rings, with an Arctic-summer
aerosol and a night layer only 70 m deep.

**Fire.** Ignited by lightning at the Open Moon's 0.003–0.006 ground strikes per km² a year,
fire keeps 0.1–39 particles per cm³ in the air, about 2 in the central case, unless lunar
fuels burn as readily as boreal forest (aerosol_budget.md Sect. 7.8); it is left out of the
area means.

## The particles, the conductivity and the cloud nuclei

At the ground, central case with the clean-to-loaded range; conductivity in units of
10⁻¹⁵ S/m, weighted with the hours in fog; the relaxation time is ε₀/σ.

| Region | Day: particles per cm³ | Day: conductivity | Night: particles | Night: conductivity | Cloud nuclei at 0.3 %, day / night |
|---|---|---|---|---|---|
| Seas | 318 (100–658) | 1.5 (0.9–3.6) | 94 (24–214) | 3.0 (1.7–5.1) | 162 / 59 |
| Lakes | 195 (28–784) | 2.4 (0.4–5.2) | 177 (16–769) | 2.4 (0.4–5.4) | 92 / 83 |
| Wet land | 195 (28–784) | 2.4 (0.4–5.3) | 224 (4–933) | 2.1 (0.3–5.2) | 92 / 56 |
| Fog desert | 197 (28–921) | 2.7 (0.3–5.4) | 146 (3–883) | 3.1 (0.3–5.4) | 79 / 39 |
| Polar dry land | 155 (62–441) | 2.4 (0.9–3.9) | 52 (9–274) | 3.8 (1.0–4.7) | 94 / 16 |
| Area mean | 226 (51–738) | 2.2 (0.6–4.7) | 153 (11–643) | 2.7 (0.8–5.2) | 110 / 53 |

- The ground's conductivity is about 2.2–2.7 × 10⁻¹⁵ S/m over the Moon, 0.6–0.7 of the
  column's single assumed aerosol (3.7 × 10⁻¹⁵), and charge near the ground relaxes in about
  an hour (55–67 minutes) against the column's 40. The loaded case reaches 5.5–10 hours over
  the forests, the lakes and the fog desert, the clean case half an hour.
- In fog the small ions go to the droplets: the conductivity falls to 8–9 × 10⁻¹⁷ S/m, 2–4 %
  of clear air's, for the 6–9 % of night hours fog covers the land; the literature's derived
  range is 1–30 % (aerosol_dust_fog.md Sect. 8.8).
- Cloud nuclei at 0.3 % supersaturation run 80–160 per cm³ by day, about the 100 droplets per
  cm³ the storm runs assume. Through the night the area mean halves (110 to 53 per cm³): by
  half or more over the seas, the fog desert and the polar land, by 39 % over the wet land and
  10 % over the lakes, whose deep night air loses little to the ground.
- Each design lever moves the area-mean conductivity by 15 % or less: no isoprene-emitting
  trees ×1.15 by day, every tree emitting ×0.87; seas at 1 or 10 nM of dimethyl sulfide ×1.10
  or ×0.93; fully crusted or bare, grazed plains ×1.01 or ×0.96; no night fungal bursts or
  twice as many ×1.05 or ×0.95 by night. The case spread, how many particles the air above
  supplies, is the larger uncertainty.

**Checks.** The same machinery with the clean Amazon's measured particles gives 770, 1,250
and 1,950 small ions per cm³ at an ionization of 2.5, 5 and 10 per cm³ per second, against
549 and 856 measured inside the rainforest at 97 % humidity (Wimmer et al. 2018, who give no
ionization rate; air over land usually takes 5–10). With the remote ocean's particles it gives
2.6–3.2 × 10⁻¹⁴ S/m at 1.5–2 ion pairs per cm³ per second, against 1.0–2.3 × 10⁻¹⁴ measured
over the Indian Ocean and the Arabian Sea (Kamra et al. 1997; Siingh et al. 2005). Both lean
high, by up to two to three times: the machinery may take up too few ions (inside a forest the
leaves take up ions it leaves out), so the conductivities above are more likely high than low.

## Open questions

- **Aloft.** The estimate is for the air near the ground. The particles at the storms'
  charging levels (30–60 km), which storms loft and few rains remove, are not estimated.
- **From the domains:** the seas' composition (hydrology), the plains' crust cover and the
  fuels' readiness to burn (biosphere), and the extent of seaweed coasts.
- **Design choices for the author:** how much isoprene the engineered trees emit, the seas'
  dimethyl sulfide, crusts and grazing on the plains, and the night decomposers.
- **The conductivity column** still carries its single aerosol; the regional states here
  could replace it near the ground.
