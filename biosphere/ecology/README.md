# Ecology of the Open Moon

The engineered organisms of the Open Moon and the communities they form: biomes and food webs, soils, the cycles
of nitrogen and the trace gases, the light and darkness life lives by, the waters, fire, dispersal and containment.
The work opened on 9 October 2026 on the branch `domain/ecology`, inside the biosphere domain, which owns biology
and ecology together. It starts from three things:
- what the project already holds ([state.md](state.md));
- the literature across ecology, plant and animal physiology, microbiology, photobiology and atmospheric chemistry
  ([literature.md](literature.md), with its [sources](sources.json));
- a first screen of the committed products ([screen.py](screen.py), [results/first_screen.json](results/first_screen.json)).


The author's standing directions frame it.
- Organisms are engineered, with Earth species' traits as reference points and design variables.
- Surface materials become developed, weathered and biologically worked soils.
- The Moon should hold a mix of climates and biomes, with pleasantly warm days at the hottest.
- The ring fleet keeps its light off the night zones.
- Imports of carbon and other elements are part of the plan, and the design problem is closing each cycle with
  local return paths.

## What the project already says

Each finding below follows from results the repository already holds, read together, with the literature that
bears on it. Numbers marked *screen* come from [results/first_screen.json](results/first_screen.json); the rest name
their source. The evidence ranges from solid physics to Earth analogues carried to the Moon, and each finding says
which.

### The light under the shield

1. **No vertebrate on the Open Moon makes vitamin D from sunlight.** The titania stack passes nothing below 320 nm.
   With the Sun overhead the ground gets no UV-B photons, against Earth's 6.8 µmol m⁻² s⁻¹, and the band that makes
   previtamin D carries 2.5×10⁻⁷ of Earth's light (*screen*). The previtamin D3 action spectrum lies below 315 nm
   (MacLaughlin et al. 1982; CIE 174:2006), and the revised in-vivo spectrum lies 5 nm shorter still (Young et al.
   2021). Earth's food webs make their vitamin D with UV-B somewhere in the chain. Here every vertebrate needs another
   route: engineered biosynthesis that runs without UV-B, calcium handling that does without the vitamin, or vitamin
   D supplied by people. Cats and dogs show that mammals living on dietary vitamin D are viable; basking reptiles are
   the hardest case. People need supplements or UV-B lamps. *Solid for the light and the chemistry.*
2. **Bees and birds lose most of their ultraviolet.** A honeybee's UV receptor (peak 344 nm, Peitsch et al. 1992)
   catches 6.9% of Earth's light with the Sun overhead and 6.4% at 30°, while its blue and green receptors catch 49%
   and 66%. Its UV signal against green is therefore a tenth of Earth's. A bird's UV-sensitive cone (370 nm) catches
   22% (*screen*). What UV remains lies at 350–400 nm. Flower UV guides and the UV skylight insects use fade.
   Engineered flowers and pollinators can signal at 350–400 nm or in the violet and visible; UV patterning in flowers
   follows the UV a habitat receives (Koski & Ashman 2015). Bees read the polarized skylight at moonlight levels
   (Dacke et al. 2003), so polarization navigation may survive on the light that is left. *Solid for the light; the
   behaviour is to test.*
3. **Sunlight barely disinfects.** Inactivation spectra falling a decade per 30–50 nm give the ground 1–6% of
   Earth's disinfecting dose (*screen*). Solar disinfection of water needs more than 108 kJ m⁻² at 295–385 nm
   (Ubomba-Jaswa et al. 2009). Pathogens therefore last longer in lakes, on wet leaves and in the air. CM1's land
   nights are 97–99% humid at 19–21 °C with fog and dew, leaf wetness of a hundred hours and more, and Earth's leaf
   fungi infect within 6–24 hours of wetness. Disease resistance becomes a design requirement of the evergreen
   strategy. *Solid in direction; the size depends on each pathogen's action spectrum.*
4. **Plants lose their UV-B cue, though only partly.** UVR8, the UV-B receptor, responds up to about 335–350 nm (Rai
   et al. 2020), where the stack passes 0.2–5%. Crops under UV-blocking greenhouse films grow as well or better, carry
   15–35% less flavonoid and anthocyanin, and see less pest pressure (Katsoulas et al. 2020). UV-B also primes the
   jasmonate defences against herbivores. The Moon's plants would grow well, colour paler and defend less. *Earth
   evidence; the size of the remaining UVR8 signal is uncomputed.*
5. **Nearside nights are brighter than any nocturnal threshold Earth's life knows.** Earthlight keeps the nearside at
   2.3–2.9 lux at its darkest (decisions register), 10–60 times typical full moonlight on Earth, 0.05–0.3 lux (Kyba et
   al. 2017), and about a hundred times the 0.01–0.03 lux at which light at night suppresses melatonin in fishes and
   rodents (Grubisic et al. 2019). Earth-calibrated nocturnal life fits the far side's low latitudes, about
   0.001 lux at midnight, which is also where the summit metropolis stands. The fleet decision keeps the night
   dark; nearside organisms still need raised circadian thresholds. The night's requirement can be set by zone.
   Earth's lunar-cued behaviour, such as coral spawning in the dark hour after sunset (Lin et al. 2021), has to be
   rebuilt on the Moon's own month, where the 29.53-day cycle is the day itself. *Solid for the light.*

### An air without its cleanser

6. **Methane lasts centuries, and a wet Moon makes plenty of it.** Under the stack the air holds no OH (the joint
   synthesis), so methane leaves only into soils. Soils taking it up at Earth's rate per area over the 62% of the
   Moon that is dry land remove 1.1–4.8 Tg a year for each ppm in the air, a lifetime of 370–1,640 years (*screen*).
   Lakes on 10.3% of the Moon emitting at Earth's 2–30 g m⁻² a year, with 0.5–2 million km² of wet margins at
   10–50 g m⁻², give 13–218 Tg a year. Methane then settles near 20 ppm at central values (59 Tg a year), and
   anywhere from 3 to 200 ppm across the ranges, approached over several centuries (*screen*). Deep crater lakes
   (median 352 m) emit less per area than Earth's shallow ponds. Bubbles leave sediment more easily at a sixth of
   Earth's water pressure, though, and skip oxidation in the water column. Methane oxidizers in soils and waters are
   therefore a requirement from the start. Wetland extent, lake depth and anoxic volume become variables of the air's
   composition (Bastviken et al. 2011). The greenhouse effect of tens of ppm in this deep column is uncomputed.
   *Solid in direction; the rates are Earth analogues, good to a factor of three to ten.*
7. **N2O has no sink at all.** N2O is destroyed by light below about 240 nm and by excited oxygen atoms made by
   ultraviolet, and the stack removes both. Soils take up at most 2% of Earth's N2O sources (Schlesinger 2013).
   Earth-like soils and seas would emit 1.0–1.3 Tg N a year, raising N2O by 0.3–0.4 ppb a year: 0.2 ppm over the
   500-year construction and 3–4 ppm in 10,000 years. The light passing the swarm's gaps gives a lifetime of 45,000
   years to 4 million years (*screen*). The biomes' nitrogen cycle has to end in N2, so complete denitrification,
   with the N2O reductase NosZ, is a design requirement. *Solid chemistry; the sources are Earth analogues.*
8. **Ethylene may decide whether cereals work.** Crop stands give off ethylene; closed chambers reach 40–120 ppb
   (Wheeler et al. 1996), and 50 ppb held continuously costs wheat 36% and rice 63% of their yield (Klassen & Bugbee
   2002). Without OH and ozone its sinks are soils and rain. The same holds for isoprene from trees and dimethyl
   sulfide from seas: seawater at 1–10 nM of dimethyl sulfide would hold the air near 2–17 ppb at equilibrium, around
   the human odour threshold (digest arithmetic, [state.md](state.md)). The aerosol study's sulfate and organic
   particles assumed oxidation that the joint synthesis later found absent. *Uncomputed; the hazard is Earth
   evidence.*
9. **"No oxidant" holds for the bare column.** Soils emitting NO at Earth's natural rate would put 15 to several
   thousand times lightning's NO into the air. NO2 still splits at a third of Earth's rate, because 350–420 nm
   passes, which can make some ozone, and ozone with alkenes or nitrous acid in light can make some OH. A living,
   peopled column needs its own chemistry run, with biological and urban emissions and their sinks (plan item 3).
   *The ratio is solid; the chemistry is open.*

### A world of long reach

10. **Pollen and spores stay aloft for weeks.** Particles fall at a sixth of Earth's speed. A 25 µm pollen grain
    falls at 3.1 mm/s (Earth 19), so a 5 km daytime mixed layer holds it for 19 days and CM1's 12 km noon layer for 45,
    against under a day on Earth. In a light 3 m/s wind 19 days carries it about 4,800 km. A 10 µm spore stays
    100–250 days, unless rain removes it first, and rain falls mostly in the tropics. Plumed seeds fall at 0.37 of
    Earth's speed and reach 2.7 times as far (*screen*). Founding biomes leak into one another; isolation distances
    for engineered plants fail; pathogens and allergens travel Moon-wide. Transgenic pollen already reaches 21 km on
    Earth (Watrud et al. 2004). Synthetic auxotrophy resists both mutation and gene transfer (Mandell et al. 2015),
    but it needs a supplied amino acid, which suits crops and bioreactors. For wild engineered organisms no
    demonstrated firewall exists, so what may live wild is an ethical and governance decision. *Solid physics.*
11. **Flight is cheap for animals too.** Hovering takes about a sixteenth of Earth's power per unit mass, and a
    kilometre on a wing a sixth of the energy (the ecology register and the sky-fleet study). Large flying grazers and
    long-range pollinators are natural on this world. Insects gain the aerodynamics, and diffusion through their
    tracheae at 1.2 atm runs at 0.83 of Earth's at the same oxygen partial pressure, which works against giant
    insects. *Solid.*

### Nutrients and soils

12. **Molybdenum may limit nitrogen fixation before nitrogen does.** Biology must fix almost all of the biomes' new
    nitrogen, since lightning gives at most 0.2% (the joint synthesis, confirmed from Schumann & Huntrieser 2007).
    Nitrogenase needs molybdenum, or vanadium or iron in its alternative forms. Molybdenum limits free-living fixation
    in weathered soils (Barron et al. 2009), and the bulk silicate Moon holds about 0.019 ppm (Sossi et al. 2024).
    Molybdenum and vanadium belong on the import list beside carbon and potassium, and the alternative nitrogenases
    belong in the engineered fixers. Fixation itself is affordable: Earth's cost of 6–12 g C per g N fixed is 0.2–2%
    of the stand's growth over the long day (digest arithmetic). *Literature cross-link; lunar molybdenum data are
    sparse.*
13. **Soil is the slow variable, and fresh basalt runs short of nitrogen first and phosphorus later.** Natural soil
    forms at 0.017–0.036 mm a year (Montgomery 2007), about 8,000–18,000 years for 30 cm, and lasting humus builds at
    0.7% of plant production (Schlesinger 1990). On Hawaii's basalt, young soils are short of nitrogen and old ones of
    phosphorus, and after 4 million years the forests live on phosphorus in dust from 6,000 km away (Chadwick et al.
    1999). The Moon has no continental dust, so its long-term phosphorus supply is a designed return. Pioneers on
    fresh basalt are nitrogen fixers and chemolithotrophs (Kelly et al. 2014). Biological crusts fix nearly half of
    Earth's land nitrogen (Elbert et al. 2012), and a moss mat held Hekla's succession for 170–700 years
    (Vilmundardóttir et al. 2018), so the order of introduction decides which state a region reaches. Microbes
    weather basalt strongly at any gravity (Olsson-Francis et al. 2012; Cockell et al. 2020). Soils are built by
    organic loading and fine material. *Earth evidence.*
14. **Long daylight breaks litter down without the soil.** Light breaks plant litter down through lignin's absorption
    in the UV and visible (Austin & Ballaré 2010). The stack removes the UV-B share, while 354 hours of daylight
    lengthen the exposure. Dry biomes, 42% of the land under 0.5 mm of rain a day, could lose litter carbon to the air
    without it entering the soil. *Direction from Earth evidence.*

### Climate, biomes and fire

15. **The biome mix comes from water and light.** The GCM's coldest 3-day mean anywhere is 14 °C, and the highest
    land is about 10 K cooler than the lowlands. No frost, alpine, boreal or seasonal biome exists on natural ground.
    The variety the author asked for on 2026-09-25 ("definitely more of a mix of climates and biomes") lies in:
    - rain: wet tropics at 4.8 mm a day, a fog desert at 30–60°, rainless poles;
    - light: the Earthlit nearside, the dark far-side equator, and poles that never fall below about
      1 µmol m⁻² s⁻¹ of light poleward of about 61°;
    - lakes, tides and soils.

    Cold places are the built ones, the winter parks at 13.5–15 km. The only year-like signals are the Earth's orbit,
    ±3.4% in sunlight, and the tides' beat, 347–412 days. *Solid on the models' climate.*
16. **Fire: ignition is rare and spread is uncertain.** Natural ignition is about a thousandth of Earth's
    (0.0016 ground strikes per km² a year from the joint synthesis). At 17.5% oxygen, tests at Earth's gravity put
    fire in the suppressed band, greatly suppressed below 18.5% and off below 16% (Belcher et al. 2010). Lunar gravity
    lowers materials' oxygen limit by 2–6 points (Ferkul & Olson 2011, in the port fire study), so natural fuels may
    burn here. No one has measured them at 17.5%, 1.2 atm and 0.16 g. Earth's fire-kept grassy biomes, over 70% of
    its burnt area (Bond & Zaloumis 2016), close over into forest where climate allows, unless people graze or burn
    them. *Ignition is solid; spread is open.*

### Waters

17. **The first seas are fresh, mixed each month, and eutrophic at first.** Weathering alone takes 1,100–89,000 years
    to bring the seas to 1 g/kg (the resources screen), so unless the delivered water carries salts the living seas
    begin as fresh to hard waters. Ice delivered fresh would leave them near pH 5.6 and undersaturated in carbonate,
    so calcifiers cannot build until weathering supplies alkalinity. New reservoirs run eutrophic for a decade
    (Paterson et al. 2019). The 354-hour heating and cooling give monthly stratification and overturn, which keeps
    waters oxygenated. Deep crater lakes, with a median depth of 352 m, risk anoxia, and productive shallows run short
    of oxygen late in the night at the night winds' gas exchange. The living-water guesses take Earth seawater at
    38.4 g/kg. *Solid if the water arrives fresh; its composition is open.*
18. **The monthly tide makes flats exposed for days.** Ranges of 2.5–3.7 m on shores sloping 1:50–1:100 expose about
    20,000–30,000 km² for days each month, the most productive kind of coast on Earth and a methane and N2O source.
    Low water drifts through the lunar day on the tidal lines' beat of about a year (digest arithmetic). *Rough.*

### Gravity and the long cycle

19. **Plants sense lunar gravity badly.** The threshold for plant gravity sensing lies at 0.1–0.3 g (Kiss 2014), and
    simulated Moon gravity disturbed root meristems more than simulated weightlessness (Manzano et al. 2018). Root
    growth regulation is an engineering target. Mice at lunar gravity keep their soleus muscle mass (2023).
    *Earth-based simulation.*
20. **Retuning crop clocks to the 708.7-hour cycle is probably the highest-value change.** Continuous-light injury in
    tomato comes from circadian asynchrony (Velez-Ramirez et al. 2017), and one locus gives tolerance with up to 20%
    more yield (Velez-Ramirez et al. 2014). Plants run their starch out about 24 hours after the last dawn whatever
    the real dawn (Graf et al. 2010). The canopy model's idling strategies already assume clocks set to the lunar
    day. *Earth evidence.*

## Statements elsewhere that need correcting

The evaluation ([state.md](state.md#statements-elsewhere-that-need-correcting)) lists stale or conflicting
statements in other folders. These belong to main and are not changed on this branch. The largest are:
- the ecology register's numbers, which predate the corrected design case;
- the aerosol study's oxidation chemistry, which predates the no-oxidant finding;
- the forests' leaf area of 22, against the canopy model's optimum near 6;
- the joint synthesis's flat 2.3 lux of Earthlight at every nearside place, against 0.48 lux at Procellarum by Russell.

## Files

| File | Holds |
|---|---|
| [screen.py](screen.py) | The first screen: ultraviolet by action spectrum, settling at lunar gravity, methane and N2O with soils as the only sink |
| [results/first_screen.json](results/first_screen.json) | Its product (schema `terluna.biosphere.ecology-first-screen/1`) |
| [state.md](state.md) | What the project already holds for ecology, what others assume of it, and statements to correct |
| [literature.md](literature.md) | The literature by topic, with what it means for the Open Moon |
| [sources.json](sources.json) | Every source cited, its use and how far it was read |

`python -m biosphere.ecology.screen` rewrites the product; `python -m pytest biosphere/tests/test_ecology_screen.py`
checks it.
