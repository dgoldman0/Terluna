# A compact model of the vacuum economy

A second approach to the economy of the vacuum buildout, beside the coupled study in [vacuum_economy](../vacuum_economy/README.md), which stays as it is. That study sizes an industry of hundreds of equipment cohorts year by year on coefficients carried from the founding design, most of them without a source, and holds that design's engineering for 200 years. This one starts from the questions the economy has to answer and keeps to the quantities that decide them:

- what each resident needs and what a resident's share of the settlement costs to build and to run;
- what the settlement buys from Earth and what it draws from its sources;
- what it sells, at each date, and at what price;
- whether the whole pays its way, how much funding it needs and when it is repaid;
- whether the whole is of a sensible size: its power, its imports, the launches and pads they need, its visitors.

Every one of them moves with the technology of the year. [technology.json](technology.json) gives each quantity of the stack (food under lamps and from electricity, power by array, storage, fission, fusion and beaming, launch, transport, processing, local manufacture, output per job, water recovery, built mass and lifetimes) as a path from 2041 to 2241 in a slow, a central and a fast case, with a static case that holds 2041 throughout. The starting values come from the lookups of 15 and 16 November 2026 where those found one; the paths after them are the study's own projections, and the first of them is being sourced layer by layer, food first.

```
python3 -m research.studies.vacuum_economy_compact.compact [--ramp NAME] [--technology CASE]   # a table by year
python3 -m research.studies.vacuum_economy_compact.compact --write                            # results/compact.json
```

The whole product, every ramp, scenario of demand and technology case, takes about a second and a half.

## How it is built

**Needs.** A resident's needs start from the founding design's own computation for 250,000 residents in four lava-tube cities and a hundred crater towns ([lava_tube_cities](../lava_tube_cities/README.md), `results/buildout.json`): about 25 kW of average power, of which food under lamps is 14.5 kW; 318 t of material made on the Moon and 17 t bought from Earth; 3.4 t of water a year. Food follows the lamps' efficacy, the food crops yield per unit of light, and the share of the diet grown from electricity by hydrogen-oxidising microbes at their efficiency; cooling follows the heat; water drawn follows the recovery; the farm floor follows the crops; the mass bought from Earth follows the share of new equipment still imported.

**Supply and costs.** Firm power is the cheapest of the options the technology offers in each year: arrays with storage through the 354-hour night, beaming from arrays in lunar orbit (the founding design's ten stations with two carrying the load), fission, and fusion once it has arrived. A kilogram from Earth lands at the launch price times a landing multiple that falls as transport improves. Local material and imported equipment carry the staged model's unit values, which have no source, falling with the price of equipment; the residents' operators follow output per job at a wage that rises with it.

**Sales.** The scenarios of demand of [scenarios.json](../vacuum_economy/scenarios.json), judged at each year's date, with the coupled study's world: sport and entertainment, hosted work and computing as stated; propellant and material at the depot and visits priced from this model's costs by the same rule as the coupled study (cost times e/(e−1) where buyers answer price strongly, the ceiling or the anchor otherwise). Visits are held to a share of the world's people a year, the study's own bound, since one elasticity at every price would otherwise carry more visitors than there are people.

**Account.** The founding build of the opening population over the eight years before the opening; each year after, building for new residents, running for all, and the lines' receipts and own costs; the funding needed, the net present value at 3% and the payback year.

## What it shows so far

On constant growth to 6.47 million residents, in trillions of 2041 dollars:

| Demand | Technology | Funding needed | Net present value at 3% | Payback year |
|---|---|---|---|---|
| Low | static · slow · central · fast | 2,232 · 187 · 3.9 · 0.6 | −22.7 · −5.1 · −0.4 · 0.3 | none · none · 60 · 32 |
| Central | static · slow · central · fast | 2,200 · 112 · 0.9 · 0.6 | −21.6 · −3.5 · 3.1 · 4.8 | none · none · 29 · 17 |
| High | static · slow · central · fast | 1,345 · 1.8 · 0.8 · 0.6 | −5.9 · 19.2 · 39.7 · 56.8 | 53 · 43 · 17 · 10 |

The technology of the period decides the outcome more than the demand does. Held at 2041, a new resident costs about 92 million dollars of 2041 to build and 5 million a year to run, most of it landing what is bought from Earth at about 1,800 dollars a kilogram, and no scenario pays. On the central path the same resident costs about 1.3 million to build and 93,000 a year to run by 2100, and even the low demand nearly pays. The front-loaded ramp, which builds sooner at dearer prices, needs about three times the funding of constant growth in the central case. In the high scenario at year 200 the checks show what that demand implies: 10 million visitors a year, 95 launches a day from 457 pads with central technology, and 3,700 pads with fast technology.

## What it does not yet do

The unit values of local work and imported equipment, the landing multiples and the bound on visits are the study's own; the power choice has no grid from the polar ridges and prices its options on those values; metals sold to Earth are not yet in the sales; the number of residents follows the ramps, not yet what the account sustains. The settlement's needs per resident are the founding design's, scaled with its mix of cities and towns.

## Evidence

The model runs on stated paths and stated unit values, and its results hold only on them. The tests check the paths' interpolation, that the static case returns the founding design's needs, that improving technology lowers a resident's needs and costs, that visits stay within their bound, and that the account adds up.
