# Aerial production, collectors and nutrient retention

The author asked on 9 October 2026 for true aeroplankton, large photosynthetic organisms and their consumers to
be evaluated together, with a functioning phosphorus cycle and food harvested from the aerial biomes themselves.
These are adopted research directions. This study quantifies their requirements. It establishes neither a
sustainable harvest nor the biomass or species richness of a completed biosphere.

Run `python -m research.studies.aerial_ecology.run`. The [product](results/aerial_ecology.json) records inputs,
producer hashes, constants, assumptions and units. The [source register](sources.json) distinguishes observations,
laboratory work and hypotheses. [Ecology's earlier screen](../../../biosphere/ecology/people_and_land.md) supplies
the atmospheric habitat and phosphorus-fallout analogue. Its once-through phosphorus arithmetic remains useful;
recycled nutrients can support further production.

## What the evidence supports

Microbes can grow in sampled cloud water, and recent fog observations support growth associated with droplets.
Photosynthetic algae and cyanobacteria have been recovered from clouds. These observations support investigating
resident aeroplankton, while leaving its net reproduction over complete airborne life cycles unresolved. The
2026 fog study concerns mainly C1-consuming heterotrophs; its activity is not a measurement of aerial carbon
fixation. Its phototrophic contribution was small. Terrestrial pollen-feeding spiders establish that biological
collectors can obtain airborne food, with a substrate taking their webs' aerodynamic loads. They do not establish
a free-flying consumer's energy balance.

Large buoyant photosynthetic organisms remain a design hypothesis. The Sagan–Salpeter Jovian ecology paper is a
theoretical precedent, with a different atmosphere and lifting-gas problem. Terluna needs its own wet-mass,
gas-production, membrane-permeability, fire, storm, reproduction and long-night solutions. Engineering organisms
is an available design direction, rather than evidence that any particular set of traits is jointly attainable.

## Reproduction must outrun removal

For an initially sparse population, the screen uses

\[
r=f\ln(2)/t_d-1/t_{\rm removal}-m.
\]

Here `f` is the fraction of the entire lunar cycle with usable water, energy and nutrients, `t_d` is doubling time
under those conditions, and `m` is additional mortality. A positive `r` is an invasion condition, not a carrying
capacity. Continuous rain removal, settling, grazing and resource limitation require a spatial life-cycle model.
An Eulerian cloud fraction cannot substitute for an organism's wet residence time. A photosynthetic organism's
favourable fraction also includes the light cycle; heterotrophy consumes organic carbon produced elsewhere.

At an assumed additional mortality of 0.01/day, a 3.6-day doubling time needs favourable conditions for **22.5%**
of the cycle if removal takes 30 days, or **10.4%** if it takes 100 days. At 19.5-day doubling, the requirements are
122% and 56.3%. The former is impossible in this equation. Three-day removal defeats even continuously growing
3.6-day organisms. The 1-day case is an optimistic engineered assumption. The 2.5-day sensitivity is informed by
a field-derived heterotroph estimate; the 3.6–19.5-day values originate in cloud-water incubations. None measures
lunar phytoplankton.

The earlier 919 million km³ warm-air volume gives 9,190–91,900 tonnes of wet cells at uniform concentrations of
10,000–100,000 cells/m³ and one picogram per cell. This is a dimensional illustration; the air is neither
uniformly wet nor uniformly populated. Low gravity slows sedimentation but does not fix precipitation removal.
Dry-particle fall times are reported separately from rain-removal sensitivities.

## Filter-feeding has a speed and drag budget

For unit frontal area, captured mass is `C v eta_capture`; chemical intake is that mass times food energy.
Assimilation reduces the usable fraction. Powered collection must also supply
`rho Cd v³ / (2 eta_propulsive)` and maintenance. All coefficients are explicit scenario inputs.

At 10 m/s, 50% capture, 60% assimilation, 25% propulsion efficiency and 1 W/m² maintenance, the previous gross
energy estimates remain 0.9, 9 and 90 W/m² for 0.01, 0.1 and 1 mg dry biomass/m³. The more informative results are:

| Effective drag coefficient | Break-even dry food concentration | Net power at 1 mg/m³ |
|---|---:|---:|
| 0.02 | 0.907 mg/m³ | +5 W/m² |
| 0.1 | 4.46 mg/m³ | −187 W/m² |
| 1 | 44.5 mg/m³ | −2,347 W/m² |

The lowest drag and 50% capture are an optimistic combination, not a demonstrated biological filter. Actual
capture efficiency and drag vary together with particle size and speed. Filter renewal, clogging, moisture loss,
the carrier's metabolism and time spent outside feeding patches can remove the small positive margin. The
product includes 1 and 3 m/s comparisons: slower collection reduces drag strongly but sweeps less food.

A tethered collector can draw flow energy from the wind. Wind-shear harvesting needs a quantified reaction path
and control scheme. An organism drifting with a uniform wind has no wind-relative feeding flow. A buoyant body
can feed by settling, active motion or local currents, whose energy costs and velocities must be supplied.

**Working direction:** favour small collectors, retained liquid microhabitats, concentrated colonies and food
around floaters; retain large mobile filter feeders as a demanding hypothesis. No whale-sized animal is inferred
from a gross intake table.

## A dense airborne food field also changes the light

For dry spheres of density 1,000 kg/m³ and extinction efficiency 2, the mass extinction is `3 Q/(2 rho d)`.
A 1-km layer holding 0.1 mg/m³ has optical depth **0.3** for 1-µm particles or **0.03** for 10-µm particles.
At 1 mg/m³ these become **3** and **0.3**. These are geometric optical screens, not Mie calculations, cloud
condensation predictions or climate forcing. Hygroscopic growth changes the effective size and scattering.

Larger aggregates reduce extinction per dry mass but generally settle faster. Low-density forms change both
budgets and must carry their own structure. Patchiness can spare other regions while preserving dense local food
fields; it reduces the space and time over which a consumer can collect. Neither a worldwide haze nor a safe
coverage fraction is adopted.

## Large photosynthetic organisms and their wet mass

The [coupled floater evaluation](../floater_viability/README.md) supersedes these lift-only values as a design
assessment. It includes pressure-bearing structure, gas barriers, active and inert mass, shared photons, water
and trim, reproduction, persistent twilight and wind-assisted Sun-following. The table below remains a historical
mass inequality; it does not select a sustainable size.

For a spherical lifting envelope, gross supported mass per projected square metre is `4 R delta_rho/3`. The skin
has four times the projected area. Allocating half the gross lift to skin and biomass, with a 1 kg/m² skin and
90% water in living biomass, yields these minimum radii:

| Altitude above sea level | 1 kg dry biomass/m² projected | 10 kg dry biomass/m² projected |
|---|---:|---:|
| 10 km | 19.1 m | 142 m |
| 20 km | 22.5 m | 167 m |
| 40 km | 31.4 m | 233 m |

The lift comes from the committed sky-ship product. These are mass inequalities, not organisms or structural
designs. The remaining half of gross lift is a reserve for unmodelled organs, bracing, control and cargo; its
sufficiency is unknown. A partly loaded organism must trim to neutral buoyancy, for example by varying gas volume
with an air-filled compartment. An unused mass allowance cannot remain as an upward force at the assumed altitude.
Biological construction of that trim system is unvalidated. The product also varies water fraction from 80% to
95%. Retained water is a substantial
lift burden. Long-lived wet structures can improve nutrient retention and provide stable habitats, but do not
eliminate gas leakage, drought, freezing or mechanical failure.

A dark half-cycle lasts about 14.77 days. Carbon loss at constant fractional respiration of 0.1%, 1% or 3% daily
would consume about 1.5%, 13.7% or 35.8% of initial carbon. That exponential screen makes respiration proportional
to the remaining pool; it omits a maintenance floor, tissue failure and storage chemistry. These are carbon-loss
sensitivities, separate from the annual NPP assumption and from the existing whole-plant survival model.
Twilight and Earthlight can shorten darkness without necessarily supplying enough
photosynthetic energy. Low-metabolism states and stored carbon must be evaluated with actual leaf traits.

## Production, standing biomass and harvest

The central arithmetic uses **1% of lunar surface area** as the sum of organism projected areas, 1,000 g C/m²/year
of net primary production, 45% carbon in dry matter and a 5% harvest. This is 379,323 km² of projected biological
area. It does not use the Moon's smaller sun-facing disk or multiply sunlight by atmospheric depth.

| Quantity | Conditional central result |
|---|---:|
| Net carbon production | 379 Mt C/year |
| New dry biomass | 843 Mt/year |
| Dry standing biomass at 1–10 kg/m² | 0.379–3.79 Gt |
| Phosphorus in that stock at Redfield C:P | 4.15–41.5 Mt P |
| Phosphorus assimilated into production | 9.23 Mt P/year |
| Dry material removed at 5% harvest | 42.1 Mt/year |
| Food-energy equivalents, all edible and recovered | 184.6 million people |
| Food-energy equivalents, half edible and 80% recovered | 73.9 million people |
| Food-energy equivalents, quarter edible and 70% recovered | 32.3 million people |

Food energy uses 4,000 kcal/kg dry matter and 2,500 kcal/person/day. It establishes neither nutritional completeness
nor a safe harvest fraction. The gross figure is less than 2% of a ten-billion lunar population. The more moderate
edible/recovery case is 0.74%. No dependable wild-aerial yield is credited toward residential carrying capacity.

The NPP assumption corresponds to 1.27 W of net biomass chemical energy per square metre of projected organism
area at 18 MJ/kg dry matter. It is a requirement for the light/physiology model. Photosynthetic floaters intercept
sunlight that would otherwise illuminate air, clouds, land or water; the net change in whole-Moon production must
subtract any displaced production beneath them. A grazer eating surface-grown food relocates that production.
Its carbon is counted once. Human harvest is removed before the illustrative 3–15% trophic-transfer sensitivities
are applied; animal production and human harvest cannot both spend the full primary-production budget.

## The cycle and the community

Phosphorus fallout at the inherited Earth analogue is 0.0266–1.062 Mt/year over the Moon. At 10% accessibility,
0.00266–0.106 Mt/year reaches the aerial organisms. This is an internal transfer from other reservoirs, with
unknown lunar rates. Capturing 10% of global transfer with 1% projected coverage requires extended collectors,
concentration, motion or repeated encounters; capture is a separate unproven requirement. Even with no harvest it
would require roughly 98.85–99.97% retention of the central annual
P throughput. Harvesting 5% exports about **0.461 Mt P/year** at the same stoichiometry: a separate recovery route
is essential. Species can depart from Redfield proportions, so the resources evaluation varies C:P rather than
treating a marine average as an engineered-floater law.

The [resources branch](https://github.com/dgoldman0/Terluna/tree/domain/resources) owns the conservative P cycle,
including sky, land, water, sediment and human recovery. Recycling fractions must describe complete routes,
including capture and delivery, rather than a laboratory extraction step alone. Nitrogen, iron, molybdenum and
other trace elements require corresponding accounts. Nitrogen fixation has its own carbon and cofactor demand.

There is room to investigate cloud phototrophs and heterotrophs, droplet grazers, long-lived photosynthetic
floaters, collectors, tissue grazers, predators and decomposers. Persistent floater microhabitats can create niche
structure even in a well-mixed atmosphere. Energy supply, isolation, life cycles, disturbance and coexistence
determine realized diversity. **No species count follows from the available atmospheric volume.**

The next discriminating calculation couples parcel hydration and light histories to reproduction, rain removal
and nutrient uptake, then closes wet floater mass, gas retention and the night budget for explicit traits.
Measured size-dependent filter capture/drag and biological aerosol optics are needed before expanding the food
or occupancy scenarios. The existing climate programme remains paused; these analytical requirements do not
authorize a new GCM campaign.


Executed numerical and repository checks, including unsuccessful full-suite gates and their causes, are in the
[check record](checks.json). The new aerial and UV checks pass; the whole repository is not reported green.
