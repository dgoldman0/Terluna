# The sky's food web

**Aerophytes are the sky's foundation species: they gather the light of open air into bodies an animal can land on,
graze and nest in.** The aerial screen found that filter-feeding the open sky barely pays even at 1–4.5 mg of food per
m³; the sky's resident aeroplankton, spread through its warm air, comes to 10⁻⁵–10⁻⁴ mg/m³, so a filter feeder there
recovers 2 × 10⁻⁷ to 10⁻⁴ of the energy it spends. An aerophyte holds its kilogram of green tissue per square metre on
a surface, where a grazer walks to its food. Kelp, coral and trees play the same part in their worlds. Numbers are
*screen* values from [food_web.py](food_web.py) (`food_web` in the [product](results/sky_ecology.json)) unless marked
*derived*; the literature values and design guesses are named once each in `RANGES`, with their sources in
[sources.json](sources.json).

## Production by class

**A growing aerophyte builds 0.43–1.09 kg C of biomass per square metre a year at the reference light, one that holds
its size 0.35–1.05 kg C, and each gives its tenants 0.1 kg C of sugar besides.** The coupled aerophyte model, called
here for each class at 10 km, gives the carbon ledger; this component splits its spending into what becomes biomass and
what does not. Of the reference's 2 kg C/m² of gross fixation, maintenance takes 0.27, night gas 0.14–0.17 and the
photolytic replacement of permeated gas 0.03–0.04; the rest renews half the green tissue a year (0.225 kg C/m²), renews
structure (0.13–0.78 kg C/m², most of it in the giants' tendons and the reefs' larger modules), gives tenants 0.1 kg C/m²
of sugar and leaves a growth surplus of 0.2–1.0 kg C/m². Building new bodies costs gas as well as tissue: at the
reference light gas is 42–68% of that cost by photolysis and 76–90% by fermentation, so only part of the surplus
becomes biomass. The sugar is counted apart from the biomass throughout.

**Two states of the sky bound its production.** In a growing sky every body spends its surplus on new bodies: at the
reference light lesser bodies add 2.1–10% to their standing biomass a year, reefs 0.4–8.4% and giants 0.3–1.1%, and in
canopy-like light up to 20%, 18% and 4%. In a mature sky the bodies hold their size and number: deaths at the adult
hazard of 0.1–0.5% a year (the aerophyte growth note) are replaced, and the surplus beyond that pays for tenants' gas and
food, domatia and sailing ([tenants](tenants.md)). A real sky lies between, its young growing and its old holding.

| Class (bodies) | Growing, reference light | Mature, reference light | Growing, canopy-like light | At 1% cover, reference light (growing; mature) |
|---|---:|---:|---:|---:|
| Lesser aerophytes, G02 (round, 50–150 m) | 0.43–0.72 kg C/m²/yr | 0.35–0.42 | 0.66–1.13 | 162–274; 134–159 Mt C/yr |
| Round giants, G03 (500 m, 1 km; 2 km in canopy light) | 0.79–1.00 | 0.70–0.93 | 1.34–2.45 | 299–380; 267–353 Mt C/yr |
| Sky reefs, G03 (rafts of 60–150 m modules) | 0.54–1.09 | 0.46–1.05 | 0.79–1.71 | 204–413; 176–397 Mt C/yr |

Growing low ends take fermented gas and high ends photolytic gas; mature ends take the hazard's ends. Over 1% of the
Moon the sky builds 162–413 Mt C of biomass a year growing and 134–397 Mt mature at the reference light, up to 928 Mt for
growing giants in canopy-like light, and gives its tenants 38 Mt C of sugar. The aerial study's central 1,000 g C/m² a
year of biomass sits at the top of the giants' and the reefs' ranges.

**The cover runs from 0.1% to 2% of the Moon.** One per cent is the aerial study's central case. The sky is shared: the
provisioning branch's five billion people aloft live in towns that shade 1.8% of the Moon (*derived*, 0.36% a billion
in its ways-of-living product), so a sky as busy as that holds aerophytes on a similar share, and a sparse wild sky a
tenth of a per cent. Every sky-wide number scales with the cover in proportion.

## Grazing, transfer and the levels of the web

**Giants take grazing as trees do: grazers taking a tropical canopy's share of the green tissue cost them 3.5–12% of
their growth surplus.** Tropical forest canopies lose 12–19% of their leaf production to herbivores (Metcalfe et al.
2014), the share applied here to the renewed green tissue; forests and shrublands lose 7.9% of their whole above-ground
production (Cebrian 2004, reproduced by Metcalfe et al. 2014) and terrestrial plants a median 18% (Cyr & Pace 1993, as
quoted by Burghardt & Schmitz 2015). A 60-m reef module regrows the grazed tissue for 3.8–6.0% of its surplus, a 100-m
lesser body for 3.5–5.5% and a 1-km round giant for 7.6–12%. Metcalfe et al. found 79% of the foliar tissue eaten by
invertebrates returned to the soil as frass, moults and bodies and 21% respired.

**The web's levels, each range at the end that bounds production:**

| Per m² of a reef a year | Low | High | Basis |
|---|---:|---:|---|
| Green tissue grazed | 27 g C | 43 g C | 12–19% of 225 g C (Metcalfe et al. 2014) |
| The host's sugar for its tenants | 100 g C | 100 g C | the reference organism's food to tenants |
| Assimilation of tissue | 30% | 50% | design guess anchored by Metcalfe et al.'s 21% respired |
| Assimilation of sugar | 90% | 90% | design guess |
| Production efficiency | 20% | 40% | poikilotherm populations about 20% (Golley 1968, in van der Meer 2021) |
| Level 2 production at the carbon limit | 20 g C | 45 g C | |
| Level 2 production at the phosphorus limit | 1.5 g C | 18 g C | food 1–3 g P/kg, 50–80% of it built in, consumers 0.97–0.62% P |
| Transfer to level 3 | 4.3% | 16% | 10.13 ± 5.81% (Pauly & Christensen 1995) |
| Invertebrate production/biomass | 4.9 a year | 2.1 a year | 0.65 M^−0.37 (Banse & Mosher 1980), 1–10 mg dry adults |

**The tenants that eat the host are short of phosphorus.** Level 2 is the grazers on the host's tissue and the guards on
its sugar. The sugar carries no phosphorus, and green tissue at 1–3 g P per kg holds too little for insect bodies at
0.62–0.97% P (Woods et al. 2004). On host food alone level 2 reaches 7.9–41% of its carbon limit at the bounding edges;
over all 64 corners of the six ranges that set it, 3.7–89%, with a median of 18% and 56 corners below half. Reaching the
carbon limit takes 0.41–0.70 g P per m² a year more at the edges (0.03–1.6 g over the corners), brought by the droppings
of roosting flyers, caught dust and pollen. The flyers that bring the sea's phosphorus up feed the whole web.

**At 1% cover the sky's grazers and guards produce 0.59–7.0 Mt C a year and hold 0.12–3.4 Mt C; its invertebrate
predators hold 0.005–0.53 Mt C.** With the phosphorus that lifts it to its carbon limit, level 2 would produce 7.4–17 Mt
C. Earth's terrestrial arthropods hold about 200 Mt C and its wild birds about 2 Mt C (Bar-On et al. 2018). The
aerophytes' own green tissue at 1 kg dry/m² over 1% of the Moon is 0.17 Gt C (*derived*), so the web on top of it is
lean, as a forest's is.

The phosphorus-limited consumer ledger reports the same tissue production as the food-web totals. Assimilated
carbon that the carbon-only efficiency would have put into tissue remains an explicit **unallocated assimilate**
pool after the phosphorus limit. Its fate could involve changed respiration, excretion, storage or intake; this
screen does not decide that physiology or credit the pool to predators or litter. Accounting closure therefore
includes that unresolved pool and is not a complete physiological model.

## How many great flyers the sky feeds

**A lunar albatross of 1.7 t burns 185 MJ a day, and a kilometre of its flight costs 0.45–0.93 MJ, a sixth of the same
animal's on Earth.** The people product carries Earth's largest flyers six times larger in every length, as the rule of
similarity allows: the wandering albatross to 18 m across and 1.7–2.1 t, Argentavis to 42 m and 15 t, Quetzalcoatlus to
60–66 m and 44–55 t. Their field metabolic rates follow the bird allometry FMR = 10.5 M^0.681 kJ/day with M in grams
(Nagy 2005): 185–210 MJ/day at the albatross scale, 824–838 at the Argentavis scale and 1,684–1,961 at the
Quetzalcoatlus scale. A species can sit at 30–170% of the prediction, and the counts below carry that spread. The
allometry puts each flyer at 39 times its Earth original. The physics brackets it: at the same speed and hours, flight
work scales as 216 times the mass at a sixth of the weight per kilogram (36 times) and upkeep as mass to the 0.72 (48
times), 40–45 times in all for flight shares of 30–70% (*derived*). Any allometric exponent between those two gives a
ratio between them, so the agreement checks consistency, and the extrapolation from tens of kilograms to tonnes stays
open. A day's FMR pays for 190–410 km of flight at the albatross scale and 66–150 km at the Quetzalcoatlus scale, at
glide ratios of 15–25 (a design guess between the sky fleet's paraglider at 10 and its sailplane at 50) and 20–25%
muscle efficiency.

**At 1% cover the sky feeds 1.8–4.8 million albatross-scale grazers, or 200,000–530,000 of the Quetzalcoatlus scale: a
quarter to nearly twice the biomass of Earth's wild birds.**

| At 1% cover | As grazers of green tissue | The same over the FMR spread | As predators on the grazers |
|---|---:|---:|---:|
| Albatross scale, 1.7 t | 1.8–4.8 million, 0.47–1.2 Mt C | 1.1–16 million | 28,000–1.6 million, 0.007–0.43 Mt C |
| Argentavis scale, 15 t | 410,000–1.1 million, 0.94–2.5 Mt C | 240,000–3.6 million | 6,200–370,000 |
| Quetzalcoatlus scale, 44 t | 200,000–530,000, 1.3–3.5 Mt C | 120,000–1.8 million | 3,100–180,000 |

The albatross-scale grazers hold 23–62% of Earth's wild birds' 2 Mt C and the Quetzalcoatlus-scale ones 66–174%
(*derived*). These are the numbers a food supply feeds at their FMR, each class taken alone. Grazers take the
herbivores' 30–50% assimilation; predators take 10–50% of the grazers' production (a design guess) at 80% assimilation
(Otero et al.'s 2018 value for seabirds) and carry the phosphorus limit above.

**Flying between bodies is cheap at 1% cover.** Bodies of 100 m stand 890 m apart, 1-km giants 8.9 km and 2-km reefs
18 km. One flight between neighbours costs an albatross-scale flyer 0.21–0.45% of its day's energy among 100-m bodies,
2.1–4.5% among giants and 4.3–8.9% among reefs; a Quetzalcoatlus-scale flyer 0.6–25%. On Earth each would cost six
times as much.

## The eaten base

**A 10-m aerophyte carries 0.36 of the reference organism's living tissue if it stores no water, close to the
two-fifths estimate.** The gross lift per projected square metre of a sphere is (2/3)DΔρ, 7.32 kg/m² at 10 m and 10 km.
With the reference skin and its ballonet wall (0.98 kg/m² dry) and the 30% reserve, the gas gate leaves room for 0.36 kg
of living dry tissue per m² at 90% water (0.38 without the wall). With the storm gust, the colony's 100 Pa of spare
pressure and half a kilogram of water it carries 0.28. Small bodies need no fire-barrier skin: a burning one is lost
whole, which a short-lived body accepts.

**Small aerophytes are costly food: 5–8% of the carbon a small body costs is living tissue a grazer can eat.** Building
an 8–40 m body costs 2.0–8.8 kg C per m² at the reference light with photolytic gas: 59–66% goes into its gas, 19–28%
into its regenerated-cellulose envelope, 6–8% into construction respiration, under 1% into its starch store and 5.0–8.1%
into living tissue; with fermented gas the tissue is 1.6–2.5%. The envelope goes to the decomposers with the litter, and
cellulose-digesting grazers, as termites are with their gut symbionts, could take a part of it; the gas is lost when a
body is eaten. A grazer eating a giant's green tissue gets 81% of what the host spent on it.

**The small bodies turn over fast, and their edible production is a third to two-fifths of a giant's renewed green
tissue.** An 8–15 m population can lose 40–78% of its numbers a year to grazers at the reference light, 64–121% in
canopy-like light with photolytic gas (12–75% with fermented gas), and hold steady. Its edible production at 8–15 m is
0.077–0.094 kg C per m² a year at the reference light, 0.34–0.42 of the 0.225 kg C/m² of green tissue a long-lived body
renews (*derived*), and 0.12–0.15 kg C in canopy-like light. Grazers take 30% (macrophytes) to 79% (algae) of it, the
medians of Cyr & Pace (1993) for aquatic producers, and on it level 2 reaches 7.5–59% of its carbon limit at the
bounding edges (median 22% over the corners), short of phosphorus as on a reef. The base works best small: at 30–40 m turnover falls to 8–30% a year.

## Sky snow

**Litter falls as sky snow and lands within 6–12 km of its track from 10 km and 41–81 km from 20 km, along the wind.**
Shed film and tissue at 0.05–0.2 kg/m² fall flat at 0.34–0.67 m/s at 10 km and 0.36–0.73 m/s at 20 km, slowing as the
air grows denser below (the sky-ship product's density profile); they take 4.3–8.6 hours from 10 km and 8.3–17 hours
from 20 km. Frass pellets of 1–3 mm fall at 1.9–3.3 m/s at 10 km. In the design run's area-weighted mean wind, which
blows west below about 3 km and east above, flakes drift 5.9–12 km east from 10 km and 41–81 km from 20 km, and pellets
1.2–2.1 km and 8.3–14 km. Through each latitude's own mean winds the slowest flakes drift at most 20 km east and 20 km
north from 10 km, and 117 km east and 18 km north from 20 km. Across the wind every drift is far smaller than the wind
product's 171-km cell, so litter lands beneath the band of the route it was shed on.

**Without pools to catch it, a reef's litter is 0.45–0.46 kg C per m² a year, 171–176 Mt C over 1% of the Moon: 4.5–4.6 g
C per m² of ground a year on average.** The ground grows 2.8–3.2 kg C/m² a year in the canopy model's clear sky (the people
product), so sky snow is a small carbon subsidy; its phosphorus matters more ([phosphorus](phosphorus.md)), and the
routes put three quarters of it on land.

## Checks

[test_food_web.py](test_food_web.py) checks:
- that the production partition closes the carbon ledger to 10⁻¹², and that the growing and mature biomass follow from
  it with the tenants' sugar left out;
- that each consumer level closes (consumption = egestion + assimilation; assimilation = production + respiration) and
  that litter counts every flow;
- that the bounding edges bound every corner of the phosphorus limit;
- that the eaten base feeds on living tissue alone and its turnover and shares close;
- Banse and Mosher's slope;
- that the allometric ratio lies between the similarity terms for the people product's own spans;
- flight at a sixth of Earth's cost per kilometre, and the glide ratios within the sky fleet's gliders;
- the FMR spread carried into the counts;
- the small body exactly at the gas gate, and the 10-m check against (2/3)DΔρ;
- the filter return;
- fall speeds through the density profile, drift through the winds and the cell read from the wind product;
- the 1% cover's area against the aerial study's.
