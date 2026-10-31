# Additional production in the deep living sky

31 October 2026. **A richly occupied sky could add roughly 1–5 Gt C of primary organic
production per year to the whole Moon**, under the reference engineered-organism traits and
10–20% total projected biological area studied here. That is approximately a 5–25% addition
to the earlier working estimate near 20 Gt C/year. The occupancy is a scenario, rather than
an inferred natural equilibrium or an adopted allocation. The present models do not establish
a likely coverage or a maximum carrying capacity; they do support a substantial aerial
contribution if water delivery, nutrient return and the organisms' combined traits work.

**Scattered light is only partly represented in this estimate.** The inherited surface
production includes the Moon's direct and diffuse light and its effect within plant canopies.
The new calculation of mutual shading among aerial organisms uses a vertical-beam column;
its homogenized surface response spreads a prescribed light reduction, but does not solve
angular scattering through the living atmosphere. The 1–5 Gt C/year result is therefore a
provisional scenario scale, not a scattering-resolved prediction or a likely interval.

The existing [surface-light calculation](../../../illumination/surface_light/README.md)
gives 37% diffuse photosynthetic light with the Sun overhead, 56% at 30° elevation and 80%
at 10°, against 7%, 12% and 30% in its Earth control. These are the older clear-sky columns,
not final 301 K fields. The [tenant-light calculation](../sky_ecology/tenants.md) includes
light scattered upward by the air below an organism, isotropic diffuse illumination from
both hemispheres and obstruction by neighbouring reef modules. At 10 km, a lone body's
underside receives 19–30% of the light on a horizontal surface at its top over ground albedos
0.05–0.2; reef undersides receive 11–22%. This supplies a physical basis for photosynthetic
flanks and undersides. Their independent carbon fixation has not been added to the new
whole-sky account.

Strong scattering can illuminate surfaces hidden from the direct Sun and distribute light
among leaves with less saturation. Consequently, complete local loss of ground production
beneath an intercepted footprint is a deliberately severe sensitivity, rather than the
expected lunar shadow. However, oblique rays also cross more organisms, and absorbing or
reflecting bodies change the diffuse field itself. Neither the direction nor the magnitude
of the whole-world correction follows from the diffuse fraction alone. An angular calculation
with the organisms, atmosphere and ground sharing photons is needed before using the table
to infer attainable cover. Applying the canopy model's diffuse-light gain a second time to
the inherited surface production would double count it.

The former 0.1–2% cover was an illustrative choice in the
[sky-ecology study](../sky_ecology/food_web.md), not a physical limit. This calculation extends
that study with a shared photon column, explicit surface displacement and resource requirements.
[run.py](run.py) reads the committed products and [parameters](parameters.json), producing
[results/aerial_productivity.json](results/aerial_productivity.json). It is a conditional
screen, with source identities and evidence states in [sources.json](sources.json).

| Sum of projected organism areas, as share of lunar surface | Primary organic production aloft, Gt C/year | Additional whole-world production after surface shading, Gt C/year |
|---|---:|---:|
| 1% | 0.26 | 0.07–0.23 |
| 5% | 1.28 | 0.36–1.15 |
| 10% | 2.56 | 0.77–2.31 |
| 20% | 5.12 | 1.71–4.63 |
| 30%, a further occupancy sensitivity | 7.68 | 2.81–6.96 |

The table uses an inherited mixture of 50% sky reefs, 30% round giants and 20% lesser
aerophytes by projected area: 60-m reef modules, 1-km giants and 100-m lesser bodies. The
reference gross fixation is 2 kg C per projected m² per year, with photolytic lifting gas,
low maintenance, night suppression and an annual adult hazard of 0.003. These are engineered
trait assumptions, whose joint attainability remains open. No new productivity is credited
merely for higher altitude, reflection from clouds, Sun-following or the area of an underside.

The quantity counted is **new host tissue and structure, including replacement, plus the
host's sugar supplied to its food web**. The old mature-host range counted only specified
tissue production. Here the representative mixture makes 0.575 kg C/m²/year of that tissue
and exports 0.100 kg C/m²/year of sugar. Consumer and decomposer production spends this
carbon; adding it again would count the same primary production twice. The calculation is
therefore an accounted part of aerial NPP, rather than a complete carbon account of every
aerial organism and metabolite.

The mature bodies' renewal and replacement programmes stay fixed while enough photosynthesis
remains to pay them. Dimming reduces their unallocated assimilate first; the screen rejects a
community whose carbon income cannot pay its fixed programme. This is why the table's accounted
production scales with organism area. It does not imply that total biological production is
insensitive to light. At 10% and 20% area, allocating half the remaining assimilate to new tissue
could add another 0.86 and 1.45 Gt C/year. Those quantities remain an allocation sensitivity:
additional nutrients, structural support, shedding, growth or other destinations for that mass
have not been established, and they are excluded from the headline estimate.

The deep atmosphere offers places for a vertically differentiated community. In the inherited
atmosphere, 10 km has a mean temperature of 14.2°C and pressure near one atmosphere; 20 km has
5.1°C and 0.81 atmosphere, with mean freezing near 25.7 km. The revised 301 K climate will
change these fields. Upper, exposed forms can gather strong light, while lower reefs can
retain water, soils and nutrients in pools and support shade phototrophs on their fringes.
Microbes, grazers, flyers and decomposers can occupy that structure without every organism
requiring its own sunny horizontal patch. The 20-km shell has only 2.3% more collecting area
than the surface, so these extra habitats share almost the same incoming solar resource.

For this purpose the area variable C is the sum of projected organism areas over all heights,
divided by the Moon's area. Four levels averaging 5% each have C = 0.20, just as one community
with the same total projected area does. Independent opaque projected footprints give the
intercepted fraction F = h[1 − exp(−C/h)], where h is the fraction of the Moon over which the
organisms occur. At h = 1, 10% and 20% area intercept about 9.5% and 18.1% of the model's light
column. Additional layers share what the layers above transmit.

The main table places the most light-demanding giants above the lesser aerophytes and the
reefs. Each class receives the photons transmitted through the classes above it; class means
sum exactly to the intercepted photon budget. Within each class, the mean dose and the dose
at its bottom plane are checked against maintenance, hydrogen production, tissue renewal and
adult replacement. The 10–20% stratified cases pass both checks, including the sensitivity
concentrating organisms over 71% of the Moon. The 30% case passes both when spread globally;
concentrating it over 71% leaves the bottom of its giant layer marginally short. A randomly mixed community can pass at its mean
while its lower giants fail. At 50% globally distributed area the stratified mixture passes
at class means but the bottom of its giant layer fails; the screen does not credit that as a fully
viable fixed vertical community. These are tests of specified traits and placement, not ceilings
on every possible aerial community. Angular light, reflected light, finite patches, shadow
duration and motion remain unresolved.

The light multiplier is 0.982, the previous comparison's midpoint provisional revised sunlight
relative to its 0.950 baseline. It is a small sensitivity while final climate tuning remains
open. The 5 kg C/m²/year gross-fixation case in the older organism study is a different
physiological scenario, not a measured gain from climbing into the sky. Dense vegetation aloft
would itself change absorption, reflection, humidity, clouds and rain. None of those climate
feedbacks is computed here.

The surface underneath matters to the additional-world result. The table's low end assumes
complete loss of production beneath the intercepted area, with representative underlying
terrain. Its high end distributes the reduction in light across the surface habitats. Rain-fed
land remains water-limited until its reduced light potential falls below that limit; wet margins
and productive coasts respond more directly to light, while nutrient-poor waters respond weakly.
These are alternative direct-light scenarios, not guaranteed bounds once climate, patch location
and ecological responses are included. Concentrating organisms over highly productive forests
can cost more surface production than these averages; over dry ground or nutrient-poor water,
less. Clouds and diffuse scattering can soften shadows, but no uncomputed global cloud bonus
is added.

Water is a matter of access and timing. Cool photosynthetic surfaces with a CO₂-concentrating
mechanism require 25–92 kg of water per projected m² per year in the inherited humid and
moderately humid cases, at 1.83 kg C/m²/year of net leaf assimilation. These rates are scaled
by area without a reduction for mutual shading. The concentrating trait's energy and its
joint physiology with the carbon-budget organism have not been recomputed. Across the new
communities the water requirement amounts to:

| Projected area | Water-efficient throughput, Gt water/year, approximately km³/year | P resupply before atmospheric fallout, tonnes/year | Ideal interception of the shared fallout analogue, tonnes P/year |
|---|---:|---:|---:|
| 10% | 96–349 | 52,000–625,000 | 2,500–101,000 |
| 20% | 192–697 | 103,000–1,251,000 | 4,800–193,000 |

Water returns to the atmosphere through transpiration; this is circulation through the organisms,
not irreversible water consumption. The warm-leaf C3 reference needs roughly five to eighteen
times as much as the efficient cases. Small reef modules can buffer only about 0.04 kg/m² of
water swing with their inherited pressure-trim allowance, so they need frequent collection.
Mean rain alone cannot demonstrate survival between storms. Sorbent surfaces, cloud collectors,
routes through humid air and stable wet microhabitats are the mechanisms the existing study
proposes; their operation together remains to be tested.

Phosphorus must return upward after litter, emigrants and whole bodies fall. The table carries
the existing forest-floor-like retention scenarios, using their losses **before** crediting
rain and dust. The former sky-need figures had already subtracted a per-footprint atmospheric
input; using them here and crediting shared fallout again would count that supply twice.
The P range also retains its inherited hazard range of 0.001–0.005/year and trait/harvest
choices, while the carbon table uses 0.003/year. Fallout is an Earth analogue with unknown
lunar rates, and the stated interception assumes perfect capture of the intersected flux.
It is an optimistic geometric comparison, not an established supply. It closes only some
favourable corners. Larger communities generally need biological transport by returning
flyers and other migrants, managed nutrient return, or both. These transfers draw from the
Moon's existing nutrient reservoirs. Nitrogen fixation, its energy cost, and iron and
molybdenum supply still lack equivalent closed aerial accounts.

Suspended phototrophic microbes could contribute, but the volume of warm air alone does not
establish a large extra crop. The separate microbial calculation holds the column's optical
depth fixed, so doubling its physical depth halves the concentration and leaves its total
biomass unchanged. At optical depth 0.1, 1–10 µm dry particles and the reference 3.6-day
division time with favourable conditions for 20% of the cycle give 0.008–0.080 Gt C/year of
cell synthesis. These duty and division rates are scenario inputs; the inherited cloud-water
incubation parameter does not establish an airborne phototroph's rate. The 30-day removal case loses cells faster than they reproduce; the 100-day
case can persist in this equation. An optimistic engineered one-day division time and 50%
favourable duty raises synthesis ninefold, to 0.072–0.720 Gt C/year, with its water, nutrient
and physiological requirements unresolved. Optical depth 0.1 already removes 9.5% from the
simple light column. These microbial quantities are not added to the aerophyte table: both
populations would need to share the photons and nutrients, and the microbial model has not
subtracted its displacement of surface production.

Earth provides precedents for viable airborne phototrophs, including cyanobacteria and green
algae recovered from cloud water, and for microbial activity responding to cloud moisture.
Those observations support exploring aerial communities; they do not measure the productivity
of a stable lunar aeroplankton population. See
[Dillon et al.](https://journals.asm.org/doi/10.1128/AEM.01850-20) and
[Péguilhan et al.](https://bg.copernicus.org/articles/22/1257/2025/).

The strongest present route to a productive deep sky is thus extensive aerophyte communities
with retained water and nutrients, differentiated by height and light demand. A contribution
of a few additional Gt C/year is compatible with the reference budgets. Establishing how much
of that becomes an actual capacity requires a population and resource calculation along the
new climate's routes, including shadows, dry gaps, nutrient return and the community's effects
on its own atmosphere.

Reproduce with `python -m research.studies.aerial_productivity.run`; check with
`python -m pytest research/studies/aerial_productivity`. The tests cover photon conservation,
carbon closure, viability gates, surface loss, optical column mass and regeneration. They
establish the arithmetic, not the organisms' empirical viability. The [check record](checks.json)
records 13 passing focused tests, substantive review and the full repository check
(2835 Python tests passed, 20 skipped; JavaScript, calendar, layer, provenance and
ensemble checks passed). These checks validate implementation and records, not biological capacity.
