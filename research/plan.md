# Core-first research plan

The next full manuscript is **Constructing and Sustaining an Open Moon**, paper `cef466e6-d9f8`. Develop its supporting studies alongside the writing; the four companions subsequently deepen their distinct arguments. Use several planetary/regional/vertical possibilities rather than treating the opening forest-and-lake scene as the whole Moon.

## Where the work stands (8 October 2026)

The scenario family is set: an open 1.2-atm atmosphere, standing water over 28% of the surface, the titania-stack
spectrum dimmed a further 5%, and water, nitrogen and oxygen delivered by the Solar-System-wide resource operation
([decisions.md](decisions.md)). The corrected GCM and the CM1 rings and boxes give the climate as a range: they
disagree on the warmth and humidity of the air over land, and the climate programme is paused. Two seas' waves and
every sea's monthly tide are computed. The light is computed from the solved spherical sky, dated for 2000–2500, with
the Earth's measured brightness and a photometric dusk 82 hours after an equatorial sunset. The storms are
electrified and their lightning measured at one site. The solar shield's lead design is the ring fleet, with its
retention screen, heat and industry. The summit port and its metropolis are sized as concepts with their flyers.

On 8 October these lines were merged in a joint integration ([integration/](integration/README.md)), and a
[synthesis](studies/joint_synthesis/README.md) set each line's assumptions about the others against their answers.
Its main findings: the ring fleet as designed lights the Moon's night at 25–66 lux; the air under the shield holds
no OH or ozone, so released gases leave only through soils, rain and escape; lightning over the Moon is a third of
the electrified box's rate and supplies at most 0.2% of the biosphere's nitrogen; the summit port's crown would start
its own lightning; haze barely touches the light; the regional magnets reach the south pole's heritage.

## The domain branches

Since 9 October three domain branches carry their own work until a joint integration: `domain/resources` (where
every material comes from and how its cycle closes), `domain/provisioning` (how resources and energy become what
people and the array need) and `domain/ecology` (engineered organisms and their communities, in biosphere/ecology).
The author's decisions in that work sit in each branch's copy of the register. The key points so far:

- **The planned climate** (provisioning). The investigated climate is the planned climate: what the shield and the
  array add to the Moon's heat is answered so that the Moon keeps it.
- **The fleet's own infrared glow** (provisioning). The films absorb 3–4.5% of the sunlight they pass, 65–83 PW over
  the fleet, and their glow puts 5.7–9.9 W/m² on the Moon, 1.9–3.4% of the design sunlight. It is answered with added
  solar protection, by redirecting it away from the Moon, or both. Films that absorb less cut it however the stack
  lies; a heat mirror on Moon-facing faces is weakened by the stack's two overlapping layers.
- **Computing and industry** (provisioning). Against the glow they do little, since keeping a watt of it off the
  Moon takes 221–303 W out of the films. On the light the fleet turns away or collects they are significant, with
  demand the limit. Every industrial use is weighed by what it puts on the films, how it dims and where its heat
  goes; a terawatt used on the Moon adds 0.026 W/m², one released in orbit 5×10⁻⁵ W/m².
- **Nitrogen** (resources). Ordinary icy bodies hold the air's nitrogen bound, and the Kuiper belt's small bodies
  alone hold 80–800 times the need, so the conservation principles hold; the author asked for the extraction to be
  planned in detail, feasible and sustainable and protecting the major bodies.
- **Opening findings** in each branch's README: the imports' redox account, delivery heat and the film plant's flows
  (resources); the night's energy, food land per person and transit at lunar gravity (provisioning); no vitamin D
  from sunlight, methane and N2O building up in an air with no OH, and pollen aloft for weeks (ecology).

On 9 October the author set near horizons beside the billion-year frame (the first thousand years with the build,
then several thousand), a working scale of 20–30 billion affluent people across the Solar System, computing as the
leading use of the array's collected energy, and the aerial biomes almost first. The studies that followed:

- **Computing** (provisioning). The night half's 851 PW outruns any computing the array can build. Mass binds first,
  7–34 Mt in orbit per terawatt, so 1,000 TW weighs 15–74% of the fleet; heat placement binds next, 0.1 K of lunar
  warming at 460–1,100 TW from the ring radius and 5,600 times less from the Earth–Sun hubs. Control and today's
  science take megawatts, a metre-scale weather twin of the Moon fills the array, AI for 20–30 billion people takes
  14–2,100 TW today, and brain emulation at scale and Matrioshka brains pass it.
- **Population** (provisioning). The scale needs growth at the pre-industrial pace and an affluent band of 2–6 kW of
  primary energy a person. People's heat on the Moon is answered by added dimming, 0.009–0.014% a terawatt used
  there. The array's habitats need 990 t of shield a person, and weekly rocket commuting costs 18–46 kW. The array
  serves Earth with computing and power and prototypes Earth's own arrays. The study's placements of people left the
  sky empty and read heat and runoff as ceilings; the author's shared comparison cases (below) replace them.
- **Industry and flows** (resources). Concentrated industry flies on platforms of its own, and heat counts by where it
  is released: per terawatt, 0.026 W/m² on the Moon, 5×10⁻⁵ from the ring radius, 10⁻⁸ from the hubs. On the near
  horizons the air needs only its CO2 returned; the array's habitats are the largest stream after the build, and the
  operation then shrinks 2,400–380,000 times.
- **The sky, land and people** (ecology). The sky is the Moon's largest habitat, wet and lightly ionized, with room for
  flyers of 60 m and food webs fed from below. The vegetated Moon's own ethylene settles at 94–233 ppb with soils its
  only sink, above the level that costs wheat a third of its yield. Each billion people living by Earth's practice
  grazes 1.2 times the fog desert and lights a third to three-fifths of the dry land; by the lunar design it crops
  2.5–16% of the rain-fed land and lights 3–6%. People's light reaches the sky's night before their crops fill the land.
- **The placement of people** follows a comparison of what each way of living holds and costs, in the shared cases
  below.

The author approved the order of the population and lifestyle work on 9 October: the checks and these corrections;
the winds by latitude and height from the saved climate run, which set the sky towns' days; the ways of living, built
from the project's designs; what each holds and costs; how people live; then distributions as outcomes for the author
to choose among. Alongside, on ecology: the sky giants' real limits (fittings, seams, gusts, damage, water gathered
from the air, storms), growth over centuries, and a taxonomic system for the biosphere, from which the floaters take
their new name.

Steps one to five and the ecology work landed on 9 October:

- **Winds aloft** (provisioning, `climate/gcm/zonal_winds.md`).
  - Above about 8 km the design climate's air superrotates, so a town drifting with it lives a day of 22 days at
    10 km and 9 at 40 km over the equator.
  - The wind holds the hour only in the lowest 5 km at 40–75°, from the late afternoon into the night.
  - Passive routes gather onto shared tracks: near the ground within months, over the poles at 40–60 km within years.
- **Ways of living** (provisioning, `studies/ways_of_living`). Eighteen settings built from the project's designs,
  each costed per billion residents.
  - Every mix carries every shared case, and towns take 0.6–10% of the dry land.
  - Food grows fastest with the Moon's people, the high end irrigated from the fresh seas.
  - Heat costs 0.08–1.8% more dimming, and light fills the dark night before towns fill the land.
  - Five billion aloft live in about 870,000 towns, with the lifting gas's upkeep their largest cost.
- **Daily life** (provisioning, `studies/daily_life`). Each setting through the 709-hour cycle.
  - People keep a 24-hour day indoors for 11–45 kWh a person a night, and the dusk lasts 82–137 hours as shared time.
  - Height buys sun, and no human has lived in partial gravity beyond 75 hours.
  - The array answers inside conversation's gap and Earth does not.
- **The floater giants** (ecology, `studies/floater_viability`, schema 2, which supersedes the version summarized in
  [integration/floater_viability.md](integration/floater_viability.md)).
  - A round body lives from about 50 m to 0.7–3 km across, by fibre and light, and a colony of 60–150 m modules has no
    structural size limit.
  - A 2-km colony grows in 36–544 years and a 1-km round giant in 73–1,186.
  - Rain meets the water need equatorward of 35–38°, and vapour and humid layers serve the higher latitudes.
  - Megaforest crowns raise juveniles to 50–60 m.
- **The taxonomy** (ecology, `biosphere/ecology/taxonomy.md`). Designed organisms keep the Earth lineage of their
  chassis, founders are registered, and descendants are named as clades from the founders. A first register holds
  36 designed groups.

Waiting on the author:
- the floaters' name, from six candidates, and six choices about the taxonomic system;
- how much of the Moon is kept wild and dark (the studies carry 30% and half);
- the distribution of people, step six, from the costs above.

## Next work, in order

1. **The Moon's night with the ring fleet.** The author has decided that the fleet avoids lighting the night zones,
   its light redirected away from the Moon and used where it can be, for computing for instance. The biosphere and
   human work set the requirement's value; the shield's night-side attitude (rolls within the bundle's clearance),
   its films' scatter and switched states are tested against it in the shield's dynamics. Flights of the kept orbits
   take about an hour of CPU each set. This joins the planned branch research/array-industry, with the trim and store
   hardware and power, heat reuse, computing and industry
   ([array_industry.md](studies/solar_shield_array/array_industry.md)).
2. **The crown as a lightning conductor.** The tower's induced charge and upward leaders in the saved storm fields,
   a protection concept for its triggered flashes and the hydrogen berths beside it. A terrain-aware electrified
   storm needs CM1's field solver over terrain, which touches the paused climate programme and waits for the author.
3. **The air's chemistry under the shield.** NO's deposition and the lunar night's chemistry in the middle
   atmosphere, then the biomes' and cities' emissions in an air with no OH; about an hour of runs once added.
4. **The biomes' nitrogen and the dark night** as biological requirements: biological fixation, and the light that
   dark-night organisms tolerate (biosphere, lunar-cycle ecology).
5. **The chosen port form's programme**, so the crown's traffic and berths follow the building the author chose.
6. **The regional magnets' design and turn**, weighed with the heritage register's written evaluations.
7. **Atmospheric electricity's stage 3**, the global circuit and the conducting upper air, at the author's go-ahead.
8. **Sea appearance's open items**: the full renderings of the four coasts on a quieter machine.

The climate programme stays paused, and no rented compute is commissioned.

## Conservation through the transformation

The [conservation study](studies/conservation/README.md) runs two programmes.

- **Scientific preservation** means sampling, excavation, storage in original condition, and research before and during the transformation.
- **Heritage preservation** means designated places, kept by engineering and grounded in ethics, with a written evaluation for each site.
- **Replacements** stand in for the Moon's scientific functions: its role as a calibration standard, far-side radio quiet, laser-ranging targets and low orbit.

Each irreversible step waits for the record it destroys, so the gates in the study enter the construction sequence beside atmosphere, water and biology.

The atlas at the selected 28% water share places everything. Its next refinements are:

- groundwater storage in the porous crust and crustal loading, which fix the delivered inventory and the shoreline;
- lake levels from each catchment's rainfall and evaporation in a climate run at 28%;
- LOLA polar stereographic grids for the cold traps.

The heritage register's written evaluation goes to the author for review. The appearance of the Open Moon from Earth, and the effect of its moonlight on Earth's life and astronomy, belong to illumination.

Water, nitrogen and oxygen come from the Solar-System-wide resource operation (study decision D4). The September baseline's local oxygen production is superseded, and the engineering companion's next seed carries the correction.

## Research threshold for the core

Seek several coherent environmental possibilities, a credible account of light, functional biological requirements for a varied biome portfolio, quantified representative spatial opportunities, and compatible construction/renewal cases. Some features can remain explicit design hypotheses. The root should derive expressive detail from this range while retaining its scientifically grounded dream purpose. The [synthesis](studies/joint_synthesis/README.md#against-the-research-threshold-for-the-core) sets the work of 8 October against it.

Current sessions support ordinary numerical models, parameter sweeps, moderate rasters, geometry, optimization and evidence synthesis. Specialist model installation and runtime need direct tests. Full plasma/dose performance, material durability, complete factories and multigenerational biological outcomes are additional validation programs. No long-running service or outside compute has been commissioned by this plan.

Use the existing editorial charter and full-source reading requirements when these results enter the manuscripts. Research notes and numerical checks are not scholarly clearance.


## Aerial biosphere and inhabited volume follow-up, 9 October 2026

The [cross-domain evaluation](integration/aerial_biosphere_population.md) supersedes interpreting the domain
branches' independent illustrative population placements as capacities. Use [one scenario file](../shared/scenarios/population.json)
for new inventories. Much floating habitation may roam, with collision/storm avoidance and maneuver reserves.

Next discriminating tasks are parcel-resolved aerial reproduction and removal; wet floater mass, lift gas and night
survival; complete phosphorus return including harvest and sediment; and roaming district control, services and
shared ecological occupancy. Restricted UV-B is an open candidate under retained short-wave rejection; adoption
needs living-air chemistry, vertical biological doses and the current non-LTE exobase/film thermal calculation.
The bounded screens do not reopen the paused climate programme or select a population distribution.

## Floater requirements follow-up, 9 October 2026

The [main summary](integration/floater_viability.md) records the ecology-owned coupled evaluation. Persistent
twilight and wind-assisted Sun-following are included. Next: one assembled wet gas-barrier test and a coupled
water/gas/pressure/flight route, then development and recruitment. The climate programme remains paused.
