# Aquatic biomes and the environmental basis

4 November 2026. Basis for the [aquatic-food study](aquatic_novelty.md), read from
main `b7eddacd674270ba8365f82a2d1bd5c3001ca883`. Exact source hashes are in
[aquatic_basis.json](aquatic_basis.json). Values below are inherited model results
or design choices, not new observations or fishery capacities.

## Baseline and version boundaries

| Input | Working basis | Consequence for food research |
|---|---|---|
| Climate design | 301 K target; latest adopted design years at this main revision give global surface mean 300.84 K, `A28R_301k_design_l15_snow_s934` years 10–19 ([plan](../../plan.md)) | Global surface temperature is not a local sea temperature |
| Sea thermal/oxygen study | [sea_oxygen](../sea_oxygen/README.md), 18 October: sea surface median 29.6 °C, approximate basin means 25–32.5 °C, on interpolated target-climate forcing | Useful warm-water habitat screen; do not relabel these columns as runs of the latest GCM years |
| Geography and forthcoming changes | Current main's atlas; the author has chosen an approximately 8° reorientation and deferred coastline changes from water loading to atlas work ([decisions](../../decisions.md)) | Named sites are provisional. Record reference frame, source revision and scenario; refresh locations before detailed stocking maps |
| Fresh seas | [sea_chemistry](../sea_chemistry/README.md): maintained water near 0.14 g/kg and pH 8.3 | Hard fresh water, chemically different from both chloride seawater and concentrated soda lakes |
| Essential minerals | [Aquatic plan](../aquatic_biogeochemistry/plan.md): the original managed sink route leaves no magnesium; a retained magnesium floor is unresolved | The old salt total is not proof of a complete physiological medium. Specify Ca, Mg, Na, K, chloride, sulfate, alkalinity and diet mineral provision before trials |
| Oxygen | Basin-specific column results under alternative mixing/trophic assumptions | Mean oxygen, saturation concentration and the minimum available to a growing fish are distinct quantities |
| Light and gravity | [shared constants](../../../shared/constants.json), [lunar ecology](../lunar_cycle_ecology/README.md), [sea appearance](../sea_appearance/README.md) | Use the actual 29.53059-day cycle and site-dependent underwater light. Low gravity changes hydrostatic gradients and settling; it does not remove swimming drag |
| Nutrients and harvest | [Sea phosphorus](../../../resources/sea_phosphorus.md), existing [food accounts](README.md) | Productive local fisheries must fit retained nutrient and oxygen budgets; vast water area alone supplies no harvest estimate |

The cold-water study exists on **domain/hydrosphere**, not in this main checkout.
Its [README at the inspected revision](https://github.com/dgoldman0/Terluna/blob/6ac19391f2728f9a2899f5f4492248653575db67/research/studies/cold_waters/README.md)
was read as a supplementary result: deep highland lakes near 21 °C, shallow high
ponds near 17 °C late in the night, with roughly 15 °C possible when unresolved
cold hours are considered. It uses the newer design years but has lapse-rate,
local-atmosphere and temporal-resolution limits. It does not justify cold-water
fish in every deep basin or establish a universal minimum temperature. No
hydrosphere commits are merged by this study.

## Basin oxygen: preserve both results and uncertainty

These are **saved-decade average oxygen concentrations in the deepest tenth of
the model's typical column**, in nutrient-poor water. They are not habitat-wide
minima, species tolerance thresholds or abundance predictions.

| Basin/habitat | Osborn central case, mg/L | Relevant bound or limitation | Food-system use |
|---|---:|---|---|
| Procellarum | 5.57 | Alternative energetic-mixing case averages at most 1.54 and was not settled | Surface pelagic programme is distinct from conditional deep-slope fisheries |
| Ingenii | 2.49 | Alternative mixing gives anoxic deep water | Surface/upper-water habitats first; no automatic deep predator allocation |
| Nubium | 5.93 | Alternative energetic form not run here | Oxygenated slope candidate pending profiles and extremes |
| Humboldtianum | 5.87 | Alternative energetic form not run here | Cooler warm-water basin comparison, not a salmon analogue |
| Fecunditatis | 5.22 | Alternative energetic form not run here | Shallow/coastal and conditional slope comparison |
| Western Procellarum shallows | 6.30 | Local productivity can increase bottom oxygen demand | Leading carbonate/prey-community test setting |
| Smythii | 0.55 | Even the central nutrient-poor deep case is severely oxygen-limited | Warm upper-water and shore study; exclude the deep typical column from initial food-fish habitats |

Values and qualifications are from [sea_oxygen](../sea_oxygen/README.md), including
the unsettled Procellarum alternative. Cooler basin ends can have different oxygen
than the typical cells. The column model can also produce short surface oxygen
crashes during overturn; a model that mixes an entire basin simultaneously is
not a resolved three-dimensional event. A food organism needs sufficient oxygen
for feeding, growth and reproduction, not merely survival through a mean value.

## Biomes to trophic roles and dishes

| Biome | Basal resources and community to evaluate | Candidate culinary niche | Principal gate |
|---|---|---|---|
| Nutrient-poor pelagic sea, productive edges | Phytoplankton and microbial loop → zooplankton → small fish; predators only after the prey budget | H1 cyclic rich small fish; later larger fillet fish; cooked whole fish and qualified raw portions | Underwater light, actual prey continuity, reproductive thermal limits, trophic export |
| Clear carbonate coast, tufa/marl structures and shellfish beds | Attached algae, biofilms, invertebrates, molluscs and grazers; structure may include mineral substrate rather than animal-built reef | H3 pale fish, edible crustaceans, selected bivalves, local grilled and steamed foods | Essential ions, habitat persistence, suspended food, night oxygen and biodeposition |
| Lake/river margins and vegetated deltas | Detritus, insects, plants, invertebrates and small fish, with controlled nutrient loads | Catfish, eel where recruitment is solved, rich cooked fillets, broth; feeder source for other habitats only if sustainable | Off-flavor, flood/recruitment cycles, sediment contaminants and oxygen minima |
| Moderately mineral-concentrating lake | Chemically characterized mixed food web with moderate carbonate/bicarbonate, evaluated separately from chloride salinity | H2 sweet/savory fish or prawns; springy cooked flesh; possible qualified raw finfish | Independent pH/alkalinity/ion effects, acclimation and full life cycle |
| Extreme soda lake | Specialist microbial producers, grazers and alkali-adapted fish only where measured tolerances permit | Mechanistic comparison, possible small specialist foods after evidence | Small body size, ion exclusion, nitrogen metabolism, conservation and source ecology; no generic shellfish farm |
| Highland lake and night-cooled pond | The modeled warm-to-mildly-cool exposure, with lake depth and shading specified | Regional texture/diet comparisons; a refuge from the hottest shore conditions | Entire cycle including warm days, not the coldest pond reading alone |
| Controlled hatchery/grow-out/finishing | Known broodstock, feeds and water treatment; saline larval system where needed | Reliable culinary reference, lineage comparisons, qualified raw/cooked products | Complete nutrient/feed/energy accounting, isolation of treatments and monitoring |

The first field map should be a **habitat eligibility map**, with unmeasured cells
left unknown. For each cell or site record: coordinate frame and scenario;
temperature range by depth and phase; oxygen duration below life-stage criteria;
light spectrum/PAR and attenuation; pH, alkalinity and individual ions;
hydraulic residence and connections; primary/prey production; substrate and
sediment load; harvest access and chilling delay. A single water-body mean cannot
stand in for these exposures.

## Budgets that keep the portfolio coherent

Maintain four linked accounts for each site and farm:

1. **Carbon and food:** production, respiration, prey grazing, fish maintenance,
   reproduction, edible harvest and all imported feed. Lipid content changes the
   energy density of prey; it does not create additional carbon or energy.
2. **Oxygen and wastes:** community respiration throughout the light cycle,
   atmosphere/water exchange, circulation, sediment oxygen demand, methane and
   nitrogen by-products. Filter feeding shifts particles into feces and
   pseudofeces as well as animal tissue.
3. **Nutrients and minerals:** phosphorus and nitrogen in all feed, living biomass,
   food, shell/nacre, feces, sediment and return flows; calcium and magnesium
   tracked separately. Shells returned or dissolved are not permanent mineral
   removal, and a harvested pearl mussel cannot also remain a living filter.
4. **Recruitment:** broodstock, spawning/larval habitat, host fish when required,
   nursery survival and return of the breeding population. Edible adults in
   suitable water do not demonstrate a self-sustaining population.

The earlier fish-equivalent phosphorus ceilings and broad aquaculture accounts
are screening constraints, not quotas that may be added together for every
candidate. New local feed or phosphorus delivery changes the inflow account and
the oxygen/greenhouse-gas burden. The author’s nutrient-poor open-water direction
remains the basis; productive margins and contained culture are separate cases.
