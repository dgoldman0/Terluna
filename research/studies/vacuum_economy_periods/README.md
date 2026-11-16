# The vacuum economy, period by period

A model of the vacuum buildout's economy from 2041 to 2241 whose foundation is the technology of each period. Every year the settlement builds and renews with the engineering of that year, chosen by cost among the options the technology then offers: its towns and cities come from their design studies, its power from the cheapest supply at each kind of site, its machines from its own factories as far as they can make them, its water from the near-Earth sources nearest first, and its money from what Earth's customers pay against what it pays Earth. It is the third approach to this economy. The coupled study in [vacuum_economy](../vacuum_economy/README.md) is closed up as it stands; the compact model in [vacuum_economy_compact](../vacuum_economy_compact/README.md) is kept and not extended.

This is the first complete version, of 16 November 2026. Its parts were designed together before it was built, and five of the design's choices are still open; each is set to the option recommended that day and carried as a switch, and their effects are measured below.

```
python3 -m research.studies.vacuum_economy_periods.run                  # a table of every case (about a minute)
python3 -m research.studies.vacuum_economy_periods.run --case central:central:constant_growth
python3 -m research.studies.vacuum_economy_periods.run --write          # results/periods.json, about 30 seconds
```

A case on a ramp of residents takes about a seventh of a second; growth as a result searches for its rate in about two seconds.

## What it shows

On the programme's constant-growth ramp, 12,500 residents at the opening in 2049 growing 3.3% a year to the planned 6.47 million in 2241, in trillions of 2041 dollars:

| Demand | Technology | Funding needed | Of it, by the opening | Funders' net present value at 3%, to 2241 | Stays repaid from year |
|---|---|---|---|---|---|
| Low | held at 2041 · slow · central · fast | 938 · 43.1 · 2.64 · 0.93 | 1.56 · 1.43 · 1.05 · 0.83 | −9.86 · −2.43 · −0.42 · 0.05 | never in any |
| Central | held at 2041 · slow · central · fast | 988 · 2.87 · 1.54 · 1.21 | 1.55 · 1.42 · 1.05 · 1.21 | −9.22 · −0.84 · 2.80 · 4.12 | never · 74 · 34 · 25 |
| High | held at 2041 · slow · central · fast | 738 · 4.79 · 3.15 · 2.45 | 1.59 · 1.46 · 1.13 · 2.45 | −1.59 · 16.3 · 28.9 · 12.6 | never · 49 · 27 · 18 |

The technology of the period decides whether the settlement pays, more than the demand does. With 2041's engineering held for two centuries no scenario pays: a settlement of millions at that technology pays Earth about 7 trillion dollars a year by 2241, for tugs and their water, supplies landed at 2041's price, equipment and work done from Earth, against receipts of 0.7 trillion. On the central path the same settlement needs 1.54 trillion dollars of funding, most of it by the opening, and has repaid it by year 34 (2075), with 2.8 trillion of net present value to its funders at 3% (82.6 trillion undiscounted, 22.8 at 1%, −0.40 at 7%).

In the central case the engineering changes as the roadmap says it can. Power comes from fission through the 2040s to 2060s, as the power lookup found megawatt reactors the lightest supported supply away from the poles, and from beamed stations from about 2070, as stations approach the 2020s gigawatt designs; the polar grid is never chosen, since its lines grow with the square of their length to sites 600 to 2,700 km away and the poles supply less than half a gigawatt on the central path. A kilowatt-hour costs 3.7 dollars in 2041, 2.7 cents in 2091 and half a cent in 2241; a kilogram landed from Earth 2,800 dollars, then 184 and 47. The settlement makes 98% of its equipment by year 30 (2071), after an imported seed and a ramp of twenty years, and buys from Earth mostly the supplies its farms and homes consume, about 450,000 tonnes a year by 2241. Its tugs fly chemical propellant from the sources' water throughout, on legs that lengthen from 0.8 km/s at the opening to 3.0 km/s by 2241 as the nearest water is drawn; the towns' lakes are much of what they carry.

By 2241 the central case receives 1.87 trillion dollars a year (computing 0.78, visits 0.65, sport and entertainment 0.28, hosted work 0.13, propellant and material 0.03) and pays Earth 1.14 trillion (what partners and Earth's services keep 0.59, the share of residents' pay spent on Earth 0.43, supplies 0.13). It uses 162 GW of firm power, carries 1.6 million visitors a year and launches 36 times a day from Earth, needing 172 pads at today's highest authorised rate. In the high scenario the checks show what its demand implies: 1,036 GW and 779 pads with central technology, and 8,000 GW and 7,000 pads with fast technology, where the residents cannot staff its lines in most years and they are cut to the hands available. The high scenario compounds its upper elasticities for two centuries and these sizes say it does not hold.

The front-loaded ramp, the founding system of 250,000 by year 33, costs more in every case: in the central case 6.02 trillion dollars of funding, a net present value of −0.29 trillion, repaid from year 65, because it builds earlier with dearer engineering for residents whose upkeep the demand does not yet pay.

### Residents as a result

When the residents are a result, the settlement grows from the opening as fast as it can while its funders still earn the programme's 3% (a net present value at 3% to year 220 that is not negative), with its own work staffed and its sources' water not exhausted, and holds at the planned 6.47 million once it gets there. Residents cost the account more than they bring, since the receipts follow Earth's demand, so this bounds their growth.

| Case | Highest steady growth | Residents in 2091 · 2241 | Reaches the plan | Funding needed |
|---|---|---|---|---|
| Central demand, central technology | 7.95% a year | 310,000 · 6.47 million | 2131 | 1.82 |
| Central demand, fast technology | 12.2% a year | 1.57 million · 6.47 million | 2106 | 2.08 |
| Central demand, slow technology | 2.23% a year | 32,000 · 0.86 million | not reached | 2.31 |
| Low demand, fast technology | 3.61% a year | 55,000 · 6.47 million | 2226 | 4.34 |
| Low demand, central or slow technology | none | 12,500 · 12,500 | not reached | 1.15 · 1.68 |

On central demand and technology the settlement could reach the plan seventy years before the ramp does and still pay its funders their 3%. The rule of the coupled study, the fastest steady rate period by period, spends the same allowance on its first years: the central case grows 30% a year to 2053 and 20% to 2059 and then hardly at all, holding about 112,600 residents for the rest of the period.

## The five open choices, measured

Each of these is carried as a switch, set to the recommendation of 16 November; the cases named `central:central:constant_growth:<switch>` in the product measure them one at a time.

1. **What the money counts.** The account with Earth counts the share of residents' pay spent on Earth, 10%. Counting the whole payroll in dollars instead, as the coupled study did, raises the central case's funding from 1.54 to 77.4 trillion dollars and leaves it never repaid: by 2241 3.6 million workers at Earth's projected wages cost 4 trillion dollars a year against receipts of 1.87. Raising the share spent on Earth to 20% lowers the net present value from 2.80 to 2.72 trillion.
2. **Births and fares.** With births at 1% a year, 4.5 million residents arrive from Earth over the period and their fares cost 185 billion dollars; with none, 6.5 million and 256 billion. Most arrive when a one-way passage costs tens of thousands of dollars (1.4 million in 2049, 55,000 in 2141), and either way the net present value moves by 0.003 trillion. Migrants paying their own fares moves it by 0.03.
3. **Towns or cities.** With the four caves of the founding system the only ones assumed, newcomers after it go to towns whatever the rule; with twenty caves the cities take 800,000 more residents and the net present value moves by 0.008 trillion.
4. **The shape of growth as a result.** The steady rule and the fastest rule are compared above; the fastest rule's early spending leaves a settlement a fiftieth of the steady rule's.
5. **The financing in U.** Not carried; money is in 2041 dollars. In the coupled study credit in U changed the funding by at most 0.25 trillion dollars.

The choices that move the result most are the first, which changes the funding fifty-fold, and the fourth, which changes how many people the same money carries.

## How it is built

**Time and stocks.** The model steps year by year from 2041 (year 0) to 2261 (year 220); the founding deployment is built in years 0 to 7 and the settlement opens in 2049 with 12,500 residents. Capital wears out geometrically, as national accounts treat it: each year a share of each stock equal to one over its life retires and is replaced with that year's engineering. Cohorts retiring all at once would echo their first build for two centuries.

**Technology.** [period_technology.json](period_technology.json) holds 28 quantities as paths through 2241 in a slow, a central and a fast case, with a static case that holds 2041, and six technologies that arrive in dated windows: fusion, a grid from the polar ridges, a mass driver on the Moon, nuclear-electric tugs, chips made on the Moon and a belt of arrays around the equator. Each value carries its basis: the lookup and source key it was read from on the research drive (engineering, local content and transport of 15 November; food, power and launch, and manufacturing, automation and labour of 16 November), or the study's own projection. The lookups of 16 November changed the roadmap: output per job grows 0.7%, 1.5% or 3% a year (4, 20 or 369 times by 2241), megawatt fission falls from 75 to 12-20 kg/kWe, the industry that builds its own machines doubles in 7, 2 or 1.4 years once established, launch to low orbit floors at 13 dollars a kilogram, a lunar mass driver needs 0.5-0.8 kWh a kilogram, and beamed stations start at the lunar designs' 12 W/kg.

**Settlements from their designs.** The founding system is the designed one: four lava-tube cities of 50,000 and a hundred crater towns of 500 ([lava_tube_cities](../lava_tube_cities/buildout.md), [vacuum_crater_buildout](../vacuum_crater_buildout/first_hundred.md)). After it, newcomers go to cities while caves remain and otherwise to towns of 750 on the 8,886 sites of the crater programme's initial quota ([crater_gardens](../../../habitation/crater_gardens/README.md)), the size at which those sites hold the plan. [designs.py](designs.py) reads each design's own quantities (regolith handled, bulk products, metal, pressure components, equipment, the first air and water, floor and loads) and works them out for the year built: farms follow the food path, structures lighten, lamps and crops change the food power. With 2041's technology the founding system returns its design's quantities. Towns built after the founding system put their share of the region's farms and light industry in halls on the surface under three metres of regolith.

**Power.** [engineering.py](engineering.py) prices a firm kilowatt-year at each kind of site from arrays with storage through the 354-hour night, beamed stations, fission, fusion once it arrives, the polar grid once it arrives and within the poles' supply, and the equatorial belt; capacity is added each year from the cheapest at the cities' and the towns' representative distances from a pole.

**Industry and labour.** Six lines (regolith handling, bulk products, metal, manufactured equipment, water at the sources and propellant at the depot) each have a plant mass per tonne a year, an energy, a life and jobs that fall with output per job. Plant is sized a quarter above the work asked of it; the bulk lines are built to the year's work, the factories from a seed bought from Earth and then from their own output no faster than their doubling, and what they cannot yet make is bought finished. Hands go to the settlements' services first, the lines of sale next (sized on the year before's balance) and the industry last; the industry's work beyond the hands left is done from Earth by teleoperation, at an Earth wage times the 1.74 that a 3-second delay cost in a laboratory trial.

**Transport and sources.** Launch from Earth at the year's price, landed at the year's multiple; passengers of 270 kg with a tonne of cabin. Lunar material lifted to the depot by rocket or, once it arrives, by mass driver. Water from the source survey's bodies nearest first ([reach.py](../../../resources/source_access/reach.py)), by the tug that costs least among chemical, electric on water and nuclear-electric, sized by the rocket equation on the leg of the water being mined; the towns' lakes and the settlements' first stocks fill over ten years.

**Sales.** The coupled study's scenarios of demand at their own date ([scenarios.json](../vacuum_economy/scenarios.json)): sport and entertainment, hosted work, computing, propellant and material at the depot, precious metals from the feed, and visits, priced by the coupled study's rule from this model's costs. A line is served in a year when its receipts exceed its annualized costs. Visits and depot sales saturate at bounds the study states, since one elasticity at every price would otherwise sell without limit as costs fall.

**Money.** The account with Earth, in 2041 dollars: equipment, supplies and hardware bought there and landed, fares, what partners and Earth's services keep, teleoperation, and the share of residents' pay spent on Earth, against Earth's customers' payments. Work on the Moon costs hands, power and plant, which the model counts physically; every choice by cost values the hands it takes at the wage.

## What it does not yet do

The rates of work, the plant masses, most paths after their sourced starting values, the mapping of Earth-orbit beaming designs onto lunar orbit, the bounds on visits and depot sales and the share of pay spent on Earth are the study's own, stated with their bases in [period_inputs.json](period_inputs.json) and [period_technology.json](period_technology.json). The settlement's own economy beyond its needs is not modelled: in the central case only about a tenth of the workforce is needed by 2241 for the work the model counts, and the rest of the residents' work goes to one another, unseen by the account. Carbon and nutrients still come from Earth, though the sources' feed holds them. Power is chosen at one representative distance from a pole for towns and one for cities. Lines are served or dropped year by year on their annualized cost, without the coupled study's search for later starts, and a line cut for want of hands is sized on the year before's balance. The industry builds its own plant on a margin of a quarter; the sources' and the depot's own power is carried in their plant. The high scenario's results rest on demand the size checks show cannot be met.

## Evidence

The model runs on stated paths and stated inputs, and its results hold only on them; it qualifies no plant, market or price. The tests check that each path holds its stated values and interpolates, and that the static case holds 2041; that with 2041's technology the founding system returns its design's quantities; that improving technology lowers a resident's food power and built mass; that the power supply chosen is the cheapest offered and that the polar grid and fusion wait for their arrival; that a chemical tug draws more water the farther it flies and stops closing; that the account adds up in every year, with visits within their bound and lines served within their hands; that better technology needs less funding; that the growth found for residents as a result is the highest the funders allow; and that the saved product is current with its sources.
