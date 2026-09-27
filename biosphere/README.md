# Biosphere

The [accepted companion seed](../ensemble/papers/74a946e2-8f98/seed-02.md) owns biological viability, human biological requirements, ecology and biogeochemical persistence across the whole Moon.

## Executable requirements models

[long_night.py](long_night.py) now provides two bounded calculations:

- Periodic carbon reserves: cycle-integrated surplus and the maximum cumulative deficit determine the minimum reservoir capacity for a prescribed production/demand cycle. The necessary-and-sufficient condition is proved under its model assumptions in [findings](../research/findings.md), with an O(N) implementation tested against exhaustive searches.
- Mixed aquatic oxygen: production, respiration, exchange and a finite concentration ceiling are balanced exactly for each interval. Outgassing and unmet aerobic demand are reported explicitly. A mathematical periodic state can still fail the selected concentration or demand criteria.

Trait rates, suppression, Q10, exchange coefficients, dissolved-oxygen ceiling and thresholds are hypothetical inputs. Units and every parameter are recorded. The models identify requirements; their outputs are not survival, health, reproduction, nutrient-cycle or ecosystem evidence. Excess carbon is unallocated surplus, not proven growth. The temperature/irradiance interface demonstrates how a climate trace alters these requirements; tissue injury remains unmodelled there, and light saturation is in the canopy model below.

[Results](../research/studies/environment_screens/results/) include 180 carbon cases, 20 oxygen cases and 18 climate-to-carbon demonstrations. Species-level calibration and long-night experiments are the next evidence step.

## Canopy photosynthesis

[canopy/](canopy/README.md) follows the illumination domain's clear-sky light at the
ground into a plant canopy by wavelength: PROSPECT-D leaves, sunlit and shaded leaves
layer by layer, and the C3 leaf model of the CO₂ study with light-limited electron
transport. At the equator, with a leaf area index of 5, the Moon's canopy absorbs 35%
fewer photons over the cycle than the same canopy under Earth's clear sky and fixes 18%
more carbon (11% more at equal CO₂): its diffuse sky spreads the light over leaves that
can use it. Through the 354-hour night the leaves alone need a store of 33 g C per m² of
ground, thirty times an Earth night's. [plant.py](canopy/plant.py) follows the whole plant
through the cycle, with the chosen climate's temperatures and the long twilight of the
Moon's tall sky, which leaves only 238 of the 354 night hours dark at the equator and 34
at 60°. Every evergreen way through the night closes its carbon budget; idling at a
quarter of the daytime upkeep cuts the store from 46 to 10 g C per m² of ground, while
regrowing the canopy every cycle leaves a tenth of the growth or none. Stomata, water,
acclimation to continuous light, full sugar stores and whether leaves survive weeks of
darkness are not yet in it. [fruit.py](canopy/fruit.py) designs fruit for the lunar
day on that plant, with the fruit's developmental program as the design variable. A
7.5 kg fruit whose program fits one sunlit half (about 200 degree-days at the equator,
a third of today's watermelon's) is set at sunrise and ripe at sunset. It fills at up
to 0.9 kg a day on the daylight surplus and never draws on the night store. With 60%
of the plant's new tissue going to fruit, the stand yields 3.7 kg per m² every lunar
cycle, 51% more than the same stand on Earth. The plant keeps enough growth to rebuild
its canopy every 2.4 cycles, and its store stays at 46 g C per m². Today's watermelon
takes 45 days, lives through a night, and needs a store of 116 g C per m² to reach
full size.

## Ecology of the lunar cycle

The [ecology register](../research/studies/lunar_cycle_ecology/README.md) collects the
design ideas built on these models: plants for the long night and the long day, the
day fruit, the warm night's food web, plains of storage organs and forests of fruit,
and the nutrient and mineral cycles, including a sea-to-forest migrant cycle. Each idea
is set against the models and the literature, and the register recommends which to
carry forward.

## Remaining biological work

The wider portfolio retains soils, aquatic communities, detrital/subsurface habitats, varied plant architectures, aerial exchange and human developmental requirements. Nutrient compartments, ecological interactions, plant hydraulics, structural support and complete life cycles still need separate models and empirical tests. Megaforests remain one candidate within this scope.

Run `python -m pytest biosphere` and `python -m research.studies.environment_screens.run` from the root. Provisional research sources and their actual access status are in [environment_sources.json](../research/studies/environment_screens/sources.json); older sources and archives remain linked through [archive_status.md](../research/archive_status.md).
