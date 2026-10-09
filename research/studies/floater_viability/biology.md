# An organism must replace its gas, tissues and descendants

The floater has a possible route to a positive mass and carbon balance. A complete organism still needs gas-tight
growth, maintenance, a tolerable night, reproduction and a supported wet community in the same design. The pure
functions in [biology.py](biology.py) keep those requirements separate. The root runner records their illustrative
outputs; [biology_sources.json](biology_sources.json) records the experimental evidence and actual access limits.

## Separate the tissues before calculating respiration

An inert extracellular lifting membrane can be produced by living secretory tissue without every gram of that
membrane remaining metabolically active. Carbon allocation should distinguish the membrane, active photosynthetic
and repair tissue, stored carbohydrate, reproductive tissue, and associated organisms. Water adds lift burden
without adding respiring carbon. A live cell layer that secretes or repairs a barrier has its own mass and nutrient
requirements. Stomata or other open gas-exchange pathways must communicate with outside air without opening the
hydrogen volume.

The vascular-plant reference is deliberately demanding. Van Iersel and Seymour's whole-vinca measurements give
estimated maintenance of 0.031 kg glucose per kg dry tissue per day and total construction substrate of 1.39 kg
glucose per kg new dry material. Their chamber conditions and model-based partition of respiration matter.
[Source](https://journals.ashs.org/downloadpdf/view/journals/jashs/125/6/article-p702.pdf).
Lettuce measurements give another active-tissue comparison: maintenance estimates of 0.031–0.039 kg glucose/kg/day,
with construction respiration depending on the fitting method.
[Source](https://doi.org/10.1046/j.0016-8025.2003.01067.x).

The sensitivities include lower maintenance rates of 0.001, 0.003 and 0.01. These test the benefit required from
lower metabolic activity or different tissue composition. They are not demonstrated floater values. Applying an
active-leaf coefficient to the whole acellular membrane would exaggerate upkeep; omitting its periodic replacement
would underestimate it. The relevant lifetime of a barrier in humid, biologically active air remains unknown.

## GPP and NPP are different budgets

`gross_carbon_budget` starts with gross carbon fixation, subtracts host maintenance, and then converts the remaining
assimilate into tissue using a construction cost that includes growth respiration. With tissue carbon fraction
0.45 and glucose fraction 0.4, the chosen construction input costs 0.556 kg assimilate C per kg new dry tissue;
0.45 kg is retained and 0.106 kg is respired. NPP is the retained carbon. A maintenance deficit is reported rather
than being concealed by a zero-growth result.

For an imposed GPP of 3 kg C/m²/year, 1 kg living dry tissue/m², equal day and night, and night maintenance at one
quarter of its daytime rate:

| Glucose maintenance per living dry mass per day | Annual maintenance | Resulting NPP | Minimum starch for one dark half-cycle |
|---:|---:|---:|---:|
| 0.001 kg/kg/day | 0.0913 kg C/m² | 2.35 kg C/m²/year | 3.32 g/m² |
| 0.003 kg/kg/day | 0.274 kg C/m² | 2.21 kg C/m²/year | 9.97 g/m² |
| 0.01 kg/kg/day | 0.913 kg C/m² | 1.69 kg C/m²/year | 33.2 g/m² |
| 0.031 kg/kg/day | 2.83 kg C/m² | 0.137 kg C/m²/year | 103 g/m² |

These illustrative GPP and suppression assumptions are not observations. At the last rate, 10 kg of living dry
tissue per projected square metre would require 28.3 kg C/m²/year for maintenance alone. A design carrying large
biomass therefore needs a small active fraction, low maintenance, greater productive surface, or higher measured
gross assimilation. Stacking surfaces does not multiply the available incident light.

When the input is the earlier assumed **1 kg C/m²/year NPP**, `net_surplus_budget` allocates that already net budget
among barrier and tissue replacement, consumers, hydrogen costs, buds and harvest. Baseline host respiration is
not charged a second time. At the inherited 18 MJ/kg dry biomass and 45% carbon, an additional assimilate-equivalent
cost of 1 W/m² uses 0.789 kg C/m²/year. Raw hydrogen chemical output is not automatically this cost: pathway losses,
gas processing and delivery must be supplied separately. The calculation is an energetic opportunity screen.

## Night storage has a floor

`dark_storage` uses a constant demand from living tissue and dependent organisms, with an optional additional
decay term for the reserve itself. For usable carbon reserve S, the model is `dS/dt = -F - kS`. Its minimum starting
reserve is `F T` when k is zero and the requested final reserve is zero. With k positive it is
`F (exp(kT)-1)/k`. A specified final reserve adds its decay-adjusted amount. The living structure is preserved;
respiration cannot fall merely because essential tissue has been consumed.

The table assumes no community respiration, reserve decay or inaccessible reserve. Each adds demand. Reserve is
additional payload, with its own tissue packaging and hydration supplied to the mass model. Host dormancy does
not automatically suppress algae, grazers, decomposers or hydrogen-producing symbionts. Small retained-water
habitats also need oxygen exchange through the night; the existing [oxygen box](../../../biosphere/long_night.py)
provides a separate dissolved-oxygen screening framework. A carbon surplus does not guarantee oxygen sufficiency.

Night reserve inventory and annual respiration must not be added as two independent annual carbon losses.
Reserve consumption supplies the respiration already counted in an annual carbon ledger. Initial construction of
the reserve, replacement of decayed reserve and reserve assigned to offspring are distinct flows.

## Hydrogen is a biological pathway with competitors

Three experimentally supported mechanisms matter. Sulfur-deprived green algae can enter a hydrogen-producing
state after oxygen production falls; that physiological switch requires recovery and resource accounting
([Melis et al.](https://doi.org/10.1104/pp.122.1.127)). Pulsed-light cultures can maintain photobiological hydrogen
production, but their productive state cannot be credited with normal growth at the same time
([Kosourov et al. 2018](https://doi.org/10.1039/C8EE00054A)). Isotope experiments established direct water
biophotolysis in a regime where the carbon-fixation cycle was inactive
([Kosourov et al. 2020](https://doi.org/10.1073/pnas.2009210117)). These results justify a hydrogen organ or symbiont
as a research direction, with its light, substrates and oxygen management counted.

Oxygen sensitivity is not an absolute prohibition. Air-grown algae have been shown to maintain protected microoxic
regions with hydrogenase activity. Hydrogen represented a small electron sink in those conditions
([Liran et al.](https://doi.org/10.1104/pp.16.01063)). A floater would still need sufficiently high sustained output
and gas capture; gross culture evolution is not the same as delivery into a nearly pure hydrogen envelope.

`photosynthetic_hydrogen_allocation` assigns an exclusive fraction of area or productive time to hydrogen.
An illustrative demand of 1 W/m² chemical hydrogen, 100 W/m² matched incident light and 2% conversion requires
half the area-time budget. The remaining half fixes carbon and must support the maintenance of the entire host.
The 0.5–2% efficiency sweep is hypothetical; laboratory pulsed-light or PAR efficiencies must not be multiplied
by unmatched annual all-solar irradiance. Gas pumping, separation and oxygen control still require budgets.

Fermentation instead consumes pre-existing organic matter. A regulated Clostridium experiment produced about
2.1 mol H2 per mol glucose, with repeated-culture yields spanning 1.63–2.32 and mixed rather than pure hydrogen
gas ([Masset et al.](https://doi.org/10.1186/1754-6834-5-35)). The corresponding gross glucose feed is
**38.5–54.8 kg per kg H2**; the ideal four-mol acetate-route comparison requires 22.3 kg. Carbon remaining in acids
or CO2 has a separate fate. Recovering it can improve system accounting, but the same substrate cannot also feed
consumers or growth without that route being closed. Product inhibition and a hydrogen-rich collection volume
make reactor yields conditional.

## Attached buds and lifetime replacement

An attached bud can receive carbon, nutrients and gas while too small to carry itself. `bud_attachment_budget`
checks both its individual margin and the combined parent–bud margin. The full growth trajectory must pass the
combined test, followed by independent lift, trim, night survival and carbon balance at release. A bud's own
photosynthesis can help, but is excluded from the simple parental allocation bound unless separately modeled.

Multicellular construction distributes secretion and repair; multiple gas compartments localise some failures.
Both require extra tissue or partitions. Dividing one spherical gas volume into two equal spherical volumes
increases exterior area by a factor of 2^(1/3), about 26%, before adding seals and organs. Fission cannot retain the
parent's envelope-area advantage without paying that cost. Tiny propagules can disperse, remain dormant or use a
surface stage; a buoyant-only life cycle must establish how they reach a self-supporting juvenile size.

`reproduction_budget` divides carbon left after adult obligations by each offspring's constructed dry tissue,
reserve and hydrogen-acquisition cost. It reports the minimum build time and the expected number of recruits
reaching adulthood over the parent's reproductive life. Construction respiration is charged only when the input
is post-maintenance assimilate; it is already paid when the input is NPP. Recruitment probability and reproductive
lifetime remain explicit unknowns. At least one expected mature replacement is necessary and does not establish
population stability under variable storms or correlated mortality.

The associated wet community shares every one of these budgets. Its food, water, night respiration and nutrient
inventory accompany the host. Buds remove phosphorus from their parent; population-wide recycling or capture must
replace losses. Carbohydrate reserves and nearly pure cellulose barriers have different P contents from active
cells, so a single Redfield ratio applied to the entire dry organism is an avoidable error. The existing resources
P-cycle model should receive separately stated inventories, tissue turnover, offspring export, harvest and
recovery fractions. Fe-dependent hydrogenases and other metabolic cofactors add trace-element requirements.

The next decisive experiment is a coupled trait set: a barrier's retained strength and gas leakage under realistic
wetting, paired with measured repair demand, productive tissue mass, hydrogen delivery, night survival and a
complete juvenile-to-adult route. Passing separate optimistic inequalities cannot establish that combined biology.
