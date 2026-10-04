# How the Open Moon's seas look

This study works out how the seas look by day, through the long evening, under
earthlight and in the dark far-side night. It couples the
[wave studies](../../../climate/waves/README.md), the
[monthly tide](../../../geography/README.md#the-monthly-tide), the solved sky and
cloud light of [illumination](../../../illumination/), the
[optical comfort](../optical_comfort/README.md) and
[cloud twilight](../cloud_twilight/README.md) studies, and the living water of the
[biosphere](../../../biosphere/). The plan below was adopted on 2026-10-04.

## What frames it

- **The evening lasts days.** The Sun sinks 0.51° an hour: its disk takes an
  hour to set, and reaching 18° below the horizon takes 35 hours. The tall
  atmosphere keeps the sky bright far below the horizon. The
  [solved sky](../../../illumination/sky/results/solved_sky.json) gives the ground
  2,210 lux with the Sun 6° down, 267 lux at 18° and 24 lux at 30°, where Earth's
  sky would be in full night. The register's practical dawn and dusk of about
  five hours still awaits a photometric definition for this atmosphere; the
  lighting calendar supplies one.
- **The nearside night is earthlit.** The Earth stays nearly fixed in the
  nearside sky and is near full through the night. Its 1.9° disk gives about
  13 lux above the air at full Earth, a ten-thousandth of direct sunlight
  ([illumination/ephemeris.py](../../../illumination/ephemeris.py): a Lambert-phase
  sphere of geometric albedo 0.367).
- **The far side's nights are dark,** lit by the stars and by an airglow not yet
  estimated for this atmosphere.
- **Around the limb the Earth rises and sets with the libration.** Over
  Smythii–Marginis, where the wave studies' headland stands, it spends part of
  each month above the horizon and part below.
- **The horizon is close:** 2.6 km for an eye 2 m above the water, 18.6 km from
  100 m.

## Decisions

The author's, on 2026-10-04:

- **The seas are living.** The optical work uses best guesses with stated bases
  and ranges, labelled as design guesses: plankton pigments by kind of sea
  (nutrient-poor open water, productive coasts and river mouths) and through the
  15-day day and night, dissolved organic matter from the rain-fed lakes and
  rivers, regolith fines near the coasts, and bioluminescence in the dark-night
  seas where the [lunar-cycle ecology](../lunar_cycle_ecology/README.md) supports it.
  Earth's ocean optics are the reference and engineered life a design variable.
- **Cloud ice keeps the cloud scheme's own sizes.** The cloud-light model
  ([microphysics.py](../../../illumination/cloud_light/microphysics.py)) takes each
  cell's ice size from the scheme's distributions, which CM1's radiation cap does
  not touch; 140 µm serves as a check. The atmospheric-electricity work found the
  box's ice reaching 157 µm while CM1's radiation reads ice only up to 140 µm, so
  the saved cloud fields grew under anvils of large crystals that their radiation
  took as too thick. The study reports the ice sizes in the clouds it uses.
- **Four coasts for the renderings,** each showing what the others do not:
  western Oceanus Procellarum by Russell, the most exposed coast, where a full
  Earth stands about 14° up in the east-southeast over the open water at night;
  the eastern Smythii headland, where the Earth rises and sets over the western
  sea; southern Mare Nubium, with the largest tide; and a South Pole–Aitken coast
  under the dark far-side night.

## Building blocks

1. **The surface's reflection** ([illumination/water_surface](../../../illumination/water_surface/)):
   Fresnel reflection of the sky and the glitter of the Sun's and Earth's disks
   over a distribution of slopes. The resolved slopes come from the SWAN spectra,
   as in the optical comfort study's water slopes; the shorter waves, which carry
   much of the slope, come from the unified spectrum of Elfouhaily et al. (1997),
   whose gravity–capillary range holds gravity and surface tension explicitly,
   evaluated at lunar gravity.
2. **The water's colour:** absorption and scattering of seawater from published
   tables, the living constituents above and suspended regolith fines, mare basalt
   and highland, from Apollo soil spectra. The same optics give the light under
   water for the ecology work.
3. **The light sources:** the solved sky, extended below 30° of solar depression
   to find where night begins; the earthlight, with a published Earth reflectance
   spectrum, and the earthlit sky from the sky solver; clouds from the cloud-light
   model; the stars, with airglow a labelled placeholder until the atmosphere
   domain estimates it.
4. **Whitecaps and foam** from the waves' breaking, at lunar gravity, where
   bubbles rise more slowly.

## Work order

1. **A lighting calendar.** Hour by hour through the month at six coasts (the four
   shores of the [shore month](../../../climate/waves/nearside.md#a-month-at-four-shores),
   the Smythii headland and a South Pole–Aitken coast): the Sun's and the Earth's
   positions, the Earth's phase, the light reaching the ground from the sky and
   from the Earth, which source dominates and which mode of vision applies. With
   it, a map of where the Earth stands over every sea and where it rises and sets.
2. **The reflection and water-colour models.**
3. **Results by regime** (day, the long evening, earthlit night, dark night):
   radiance and colour by direction of view, the glitter paths and the contrast of
   the waves.
4. **Scientific renderings of the four coasts** in [visualization](../../../visualization/),
   extending the [B1 spectral path tracer](../../../visualization/reference-renderer/),
   judged against photographs of seascapes.

Models belong to their domains (optics in illumination, sea states in climate,
the tide and the seas in geography, the living water's guesses in the biosphere);
this folder holds the coupled runners, their results and the write-up.
