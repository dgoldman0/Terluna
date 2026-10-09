# The organism's carbon: upkeep, gas, night and offspring

**A floater can carry a positive carbon balance when most of its mass is inert, its living tissue is light and its
upkeep is low: at 1 kg of living dry tissue per square metre and 0.003 kg glucose per kg a day of maintenance, 3 kg C
of gross fixation leaves 2.2 kg C of net growth.** The same tissue at a whole-plant rate of 0.031 kg/kg a day keeps
only 0.14 kg C, and 10 kg of living tissue per square metre at that rate needs 28 kg C a year for upkeep alone. Light
living tissue on a long-lived inert envelope is therefore the design. [biology.py](biology.py) keeps each budget
separate; [biology_sources.json](biology_sources.json) records the experiments and how far they were read.

## Living tissue and inert membrane

**An inert extracellular membrane, secreted and repaired by a thin living layer, carries no respiration of its own.**
Carbon allocation separates the membrane, active photosynthetic and repair tissue, stored carbohydrate, reproductive
tissue and the community the floater carries. Water adds lift burden without respiring. Stomata or other gas-exchange
paths open to the outside air while the hydrogen volume stays closed.

**Whole-plant measurements set the demanding end.** Van Iersel and Seymour's whole-vinca gas exchange gives a
maintenance of 0.031 kg glucose per kg dry tissue a day and a construction cost of 1.39 kg glucose per kg of new dry
material, at 22/18 °C and 14/10-hour days ([source](https://journals.ashs.org/downloadpdf/view/journals/jashs/125/6/article-p702.pdf)).
Lettuce gives 0.031–0.039 kg/kg a day, with construction respiration depending on the fit
([source](https://doi.org/10.1046/j.0016-8025.2003.01067.x)). The sensitivities 0.001, 0.003 and 0.01 test how much
less costly tissue or suppressed metabolism would buy.

## Gross and net production

**`gross_carbon_budget` subtracts maintenance from gross fixation and builds tissue with the rest.** With 45% carbon in
dry tissue and 40% in glucose, building a kilogram of dry tissue takes 0.556 kg of assimilate carbon: 0.45 kg kept and
0.106 kg respired. Net production is the carbon kept, and a maintenance deficit is reported as such. For 3 kg C of
gross fixation, 1 kg of living dry tissue, equal day and night and night maintenance at a quarter of the day's rate:

| Maintenance, kg glucose per kg living dry mass a day | Upkeep a year | Net production | Least starch for one dark half-cycle |
|---:|---:|---:|---:|
| 0.001 | 0.0913 kg C/m² | 2.35 kg C/m²/yr | 3.32 g/m² |
| 0.003 | 0.274 kg C/m² | 2.21 kg C/m²/yr | 9.97 g/m² |
| 0.01 | 0.913 kg C/m² | 1.69 kg C/m²/yr | 33.2 g/m² |
| 0.031 | 2.83 kg C/m² | 0.137 kg C/m²/yr | 103 g/m² |

Stacking surfaces does not multiply the light that reaches them.

**`net_surplus_budget` divides a net production among barrier and tissue renewal, consumers, hydrogen, buds and
harvest.** At 18 MJ/kg dry and 45% carbon, each 1 W/m² of extra assimilate cost uses 0.789 kg C/m² a year. Host
respiration is already inside net production and is charged once.

## The night store

**The least starting reserve for a dark interval T at a constant demand F is F T, and F(e^{kT} − 1)/k when the reserve
itself decays at rate k.** `dark_storage` solves dS/dt = −F − kS with living structure kept intact; a final reserve adds
its decay-adjusted amount. The reserve is extra payload with its own packaging and water. Host dormancy leaves the
algae, grazers, decomposers and hydrogen symbionts with their own demands, and small water habitats need oxygen
through the night (the [oxygen box](../../../biosphere/long_night.py) gives a framework). The night store supplies
the respiration already in the annual ledger; building, decay and the reserve given to offspring are the separate
flows. Shorter nights aloft cut the store in proportion ([sailing](sailing.md)).

## Hydrogen from light and from sugar

**Algae make hydrogen from light in states that compete with carbon fixation.** Sulfur-deprived green algae enter a
hydrogen-producing state after oxygen evolution falls, a switch with recovery and resource costs
([Melis et al.](https://doi.org/10.1104/pp.122.1.127)). Pulsed light sustains photoproduction, with a reported early
light-to-hydrogen peak near 1.6–1.7% in cells that did not grow
([Kosourov et al. 2018](https://doi.org/10.1039/C8EE00054A)). Isotopes show water oxidation by photosystem II as the
electron source when carbon fixation is inactive ([Kosourov et al. 2020](https://doi.org/10.1073/pnas.2009210117)).
Air-grown cells keep microoxic niches with active hydrogenase, though hydrogen took under 1% of the electron flow
([Liran et al.](https://doi.org/10.1104/pp.16.01063)).

**Photolytic hydrogen takes area and time from carbon fixation one for one.** `photosynthetic_hydrogen_allocation`
gives hydrogen an exclusive share of area-time: 1 W/m² of hydrogen from 100 W/m² of matched light at 2% takes half
the budget, and the other half pays the whole organism's upkeep. The 0.5–2% sweep is hypothetical, and laboratory
efficiencies apply only on the same spectrum and duty.

**Fermentation commits 22–55 kg of sugar per kilogram of hydrogen.** A regulated *Clostridium* culture made about
2.1 mol H₂ per mol glucose, 1.63–2.32 in repeated batches, in mixed biogas of 62–65% hydrogen
([Masset et al.](https://doi.org/10.1186/1754-6834-5-35)): 38.5–54.8 kg of glucose per kg H₂, against 22.3 kg at the
four-mole acetate ceiling. The carbon left in acids and CO₂ has its own fate, and the same sugar cannot also feed
growth. The [growth](growth.md) note turns these costs into ages.

## Buds, offspring and replacement

**A bud can draw carbon, nutrients and gas from its parent before it lifts itself.** `bud_attachment_budget` checks
the bud's own margin and the parent and bud together; release follows when the bud has its own lift, trim, night
reserve and carbon balance. Dividing one gas volume into two equal spheres adds 26% to the outer area, 2^(1/3), before
seals and organs.

**`reproduction_budget` divides the carbon left after the adult's obligations by each offspring's cost.** It returns
the least time per bud and the expected recruits reaching adulthood over the parent's reproductive life; at least one
is needed. Construction respiration is charged once, at the stage the input names. The community aboard shares all
these budgets. Buds take phosphorus from their parent, and carbohydrate and cellulose hold little phosphorus while
active cells hold more, so each part needs its own content ([growth](growth.md) gives the phosphorus guesses).
Iron-dependent hydrogenases add trace-element needs.

The decisive experiment is one trait set measured together: a barrier's retained strength and gas loss when wet, its
repair demand, the productive tissue's mass, hydrogen delivery, night survival and a full juvenile-to-adult route.
