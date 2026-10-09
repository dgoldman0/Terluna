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
  primary energy a person. The Moon's heat budget holds 1.4–2.2 billion at 0.1 K with industry in orbit, and far more
  with added dimming; its runoff binds near 4 billion, which fresh seas may lift. The array's habitats need 990 t of
  shield a person, and weekly rocket commuting costs 18–46 kW. The array serves Earth with computing and power and
  prototypes Earth's own arrays.
- **Industry and flows** (resources). Concentrated industry flies on platforms of its own, and heat counts by where it
  is released: per terawatt, 0.026 W/m² on the Moon, 5×10⁻⁵ from the ring radius, 10⁻⁸ from the hubs. On the near
  horizons the air needs only its CO2 returned; the array's habitats are the largest stream after the build, and the
  operation then shrinks 2,400–380,000 times.
- **The sky, land and people** (ecology). The sky is the Moon's largest habitat, wet and lightly ionized, with room for
  flyers of 60 m and food webs fed from below. The vegetated Moon's own ethylene settles at 94–233 ppb with soils its
  only sink, above the level that costs wheat a third of its yield. With Earth's practice the Moon holds about half a
  billion affluent people; with the lunar design, about 3–15 billion; people's light darkens the sky's night before any
  diet reaches its crop limit.
- **Open across the branches:** the placement of people. The resources screen reads the scale with Earth at 10.3
  billion and 2–14 billion in the array, the population study with 16–20 billion on Earth and 0.1–1 billion in the
  array; the author's placement settles both.

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
