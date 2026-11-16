# The economy of the vacuum buildout

This study follows the economy of the settlements beyond Earth through the vacuum period: from the opening of production in year 8 of the programme (2049) to the planned population of 6.47 million residents in year 200, and for twenty years past it. Its unit of account is U, the general economic unit of the whole planned system beyond Earth before the Open Moon, the array and the outer Solar System's buildout. For each case it draws up what customers order in each year, sizes the industry that fills those orders, finances that industry, and reports what the funders, the residents and the sources each give and receive.

The study couples parts that live in their own lanes and reads each where it is kept. The cohort model of industry is [provisioning/vacuum_industry](../../../provisioning/vacuum_industry/README.md), with its measures of material recovery and equipment refurbishment. The financing model is [enclosed_returns/finance.py](../enclosed_returns/finance.py). The run in this folder joins them for one case at a time.

## The programme definition

[programme.json](programme.json) states the calendar, the money and the population once, and [programme.py](programme.py) reads it for every other part. Year zero is 2041; production opens in year 8 with 12,500 residents; the planned population of 6,468,289 residents, 0.1% of the full buildout, stands in year 200; and every amount is in constant dollars of 2041. Three ramps carry the population from the opening to the plan:

| Ramp | Shape | Residents in year 33 | Residents in year 100 |
|---|---|---|---|
| `constant_growth` | 3.3% a year throughout | 28,200 | 249,600 |
| `front_loaded` | 12.7% a year until the founding system of 250,000 residents stands in year 33, then 2.0% a year | 250,000 | 922,100 |
| `first_pass` | the staged studies' schedule, slow at both ends | 16,100 | 222,200 |

The first two are the programme's cases. The third is kept so that the staged studies' figures can be reproduced. The run refuses a case whose order book disagrees with the definition on the residents of any year, and it refuses to start if the staged inputs name a different calendar, population, horizon or discount rate.

Two further ramps state a rule in place of a schedule. On `resident_led` and `resident_led_beyond_plan` the number of residents is a result: the settlement grows as fast as its own production and financing sustain, held at the planned population from the year it reaches it on the first and not held before completion on the second. [residents.py](residents.py) finds their schedules, as the section on the residents as a result describes, and the definition refuses a solved ramp whose residents it is not given.

## The run

```
python3 -m research.studies.vacuum_economy.run [--case NAME ...] [--jobs N]
```

The run takes the cases listed in [cases.json](cases.json) and writes [results/economy.json](results/economy.json), 2.2 MB for 68 cases. It reads its inputs through [inputs.py](inputs.py), as residents.py does, and a case on a solved ramp reads its schedule from [results/residents.json](results/residents.json), which must be produced first and is refused where any file it was computed from has changed. A first-pass case takes about ten seconds on one thread, and about twenty where the staffing of local manufacture is solved; a scenario case sizes the industry some twenty times and takes five to ten minutes. `--jobs` solves that many cases at once. Each record holds the funding at the opening; for each of three financing arrangements (equity alone, funded credit, and credit in U) and each horizon (years 100, 200 and 220) the receipts, the operating cash, the payouts, the equity added after the opening, the rate of return and the net present values at 0%, 1%, 3% and 7%; for each of seven periods the receipts by sector and the costs by category, with the share of the costs paid to Earth; the asteroid feed, the lunar regolith, the water, the nitrogen and the fresh metal drawn by each horizon; the traffic of each horizon's year, as launches from Earth with the pads they imply and the round trips of the ferries, the lunar landers and the tugs; the jobs of that year with the share of the machinery's working parts made in place; and what the mined feed yields by each horizon beside what the system uses and sells of it. A scenario case adds the account of each line of sale: whether it is served and from which year, its worth, its receipts, cost and jobs by period, and its quantity and price at each horizon. A case may change engineering coefficients, each named by its path, to a value or to a schedule by year. The product names its producer and the digest of every file it was computed from.

## The first pass

The scenario `first_pass` is the demand of the staged studies of 12 to 14 November 2026: their nine export portfolios at their strongest market, with real demand growing 0.5% a year, with material recovery and equipment refurbishment applied. On the `first_pass` ramp, through this entry point, every figure of the [staged appraisal](../enclosed_returns/revised.md) is reproduced exactly; the table repeats its headline, in trillions of 2041 dollars with credit in U and a 3% discount rate, and adds two columns the staged write-up does not carry.

| Portfolio | Equity added by year 200 | Operating cash a year, year 200 | Payouts by year 200 | Funders' net present value, year 200 | Costs paid to Earth, years 151–200 | Asteroid feed by year 200 |
|---|---|---|---|---|---|---|
| Raw-price offer | 48.08 | −1.18 | 0.00 | −3.89 | 44% | 134 trillion kg |
| Commodity offer | 31.05 | −0.19 | 0.00 | −3.83 | 48% | 139 trillion kg |
| Specialty products | 12.98 | 0.72 | 1.20 | −3.85 | 64% | 193 trillion kg |
| Premium lunar products | 8.27 | −1.53 | 0.00 | −3.87 | 64% | 213 trillion kg |
| Propellant and logistics | 0.00 | 0.75 | 119.93 | 8.79 | 44% | 146 trillion kg |
| Computing and digital | 58.38 | −1.23 | 0.00 | −3.99 | 49% | 139 trillion kg |
| Worldwide services | 48.08 | −1.18 | 0.00 | −3.89 | 44% | 134 trillion kg |
| Diversified | 0.00 | 1.54 | 136.45 | 7.35 | 63% | 199 trillion kg |
| Earlier manufacturing mix | 0.00 | 3.26 | 220.56 | 9.08 | 57% | 178 trillion kg |

The two added columns show where the first pass leaves the system's own production out of its projections. The model makes structure locally and buys the working parts of every new machine from Earth throughout its 220 years, at one price for each kilogram and one price for its launch, 200 dollars of 2026 a kilogram to low Earth orbit; refurbishment returns part of the retired machinery to service. As the industry grows, the purchases grow with it: in the diversified case they are 4% of all costs in years 8 to 25 and 63% in years 151 to 200. The asteroid feed is the quantity of rock that yields the water the system needs, which exceeds the rock its nitrogen and its metal need in every year. In the diversified case 94% of that water through year 200 is propellant: 54% for the tugs that carry the feed, 21% for the ferries between the Moon and orbit, 18% for the ferries between Earth orbit and the Moon, and 2% sold to customers. The residents and their guests use the remaining 6%. The source study's distribution of what the near-Earth objects hold is compared with the draw after the run and changes nothing in it.

Five further choices of the first pass bear on the same point. Its visitors pay for their nights and no fare, so the cost of a trip is absent from the number who come. Its sales are capped by four budgets of buyers on Earth that are set beside the model and grow by the stated rate whatever the system can offer. Its materials earn receipts when a capsule returns them to Earth, and its only customers in space are missions that buy propellant and its handling. Its returns credit no coproducts. And its population follows the staged schedule in every case.

The first pass counts what the system sells in founding systems. Its order book takes the export mix of the founding system of 250,000 residents and offers 5% of it at the opening and 4% more in each later year, whatever the industry of that year could make. It staffs and costs its premium products and its logistics per founding scale, takes its visitors as the founding figure times a ratio of residents, and adds service facilities as identical units of 350 staff, each of which brings 65 billion dollars a year. The sizing of equipment was corrected on 13 November to follow the work to be done; the order book was left as it was. The three scenarios below replace it, and the first pass stays as the reference against which they are measured.

## The first pass on the programme's ramps

The same nine portfolios run on the programme's two ramps, with nothing changed but the residents of each year. Every case solves and passes its capacity screens. Each cell gives the staged schedule, constant growth and the front-loaded ramp in that order, in trillions of 2041 dollars with credit in U and a 3% discount rate, and the feed in trillions of kilograms.

| Portfolio | Equity added by year 200 | Payouts by year 200 | Funders' net present value, year 200 | Asteroid feed by year 200 |
|---|---|---|---|---|
| Raw-price offer | 48.08 · 7.35 · 47.50 | 0.00 · 0.00 · 0.00 | −3.89 · −3.70 · −3.92 | 134 · 104 · 149 |
| Commodity offer | 31.05 · 0.00 · 32.67 | 0.00 · 4.87 · 0.00 | −3.83 · −0.99 · −3.87 | 139 · 109 · 154 |
| Specialty products | 12.98 · 0.00 · 14.82 | 1.20 · 32.94 · 4.97 | −3.85 · 3.74 · −3.89 | 193 · 163 · 207 |
| Premium lunar products | 8.27 · 0.00 · 7.09 | 0.00 · 2.62 · 0.00 | −3.87 · −2.24 · −3.92 | 213 · 184 · 228 |
| Propellant and logistics | 0.00 · 0.00 · 0.00 | 119.93 · 162.02 · 121.14 | 8.79 · 10.46 · 7.23 | 146 · 116 · 161 |
| Computing and digital | 58.38 · 17.68 · 57.91 | 0.00 · 0.00 · 0.00 | −3.99 · −3.78 · −4.02 | 139 · 110 · 154 |
| Worldwide services | 48.08 · 7.35 · 47.50 | 0.00 · 0.00 · 0.00 | −3.89 · −3.70 · −3.92 | 134 · 104 · 149 |
| Diversified | 0.00 · 0.00 · 0.00 | 136.45 · 181.59 · 138.41 | 7.35 · 9.04 · 5.74 | 199 · 170 · 214 |
| Earlier manufacturing mix | 0.00 · 0.00 · 0.00 | 220.56 · 263.93 · 222.08 | 9.08 · 10.77 · 7.45 | 178 · 148 · 192 |

Constant growth improves every portfolio on the staged schedule. The staged schedule is slow at both ends and steep between: it grows 5.1% a year at its steepest, in year 106, and houses two million new residents in years 101 to 150, where constant growth houses one million. That middle stretch concentrates construction and, twenty years on, the replacement of what it built, in decades when the receipts are the same on either schedule. Spread evenly, the same build needs 7 trillion dollars of added equity in the services portfolio where the staged schedule needs 48, needs none in the commodity, specialty and premium portfolios, and draws about 30 trillion kg less feed in each portfolio.

The front-loaded ramp houses 250,000 residents by year 33, sixteen times the staged schedule's 16,100 in that year. In the first pass's demand the early services are limited by the residents who staff them, so receipts rise with the residents: in the diversified portfolio the receipts of years 26 to 50 are 55 trillion dollars against 25, and the costs 41 against 15. Through year 100 the ramp adds about 72 trillion dollars of costs and 66 trillion of receipts to that portfolio, the costs arriving first, and the funders' net present value at year 200 falls by about 1.6 trillion dollars in each of the three portfolios that finance themselves. The equity the others need is within 2 trillion dollars of the staged schedule's, and the feed is about 15 trillion kg higher in each portfolio, from the larger number of resident-years.

## Where the Earth leg's propellant comes from

Every tonne the settlement imports crosses from low Earth orbit to the depot on a ferry. In the staged model the depot fuels both of the ferry's burns from asteroid water, at 6.2 tonnes of propellant for each tonne of cargo, and the fleet of tugs and mines that supplies it is built before production opens. The [industry model](../../../provisioning/vacuum_industry/README.md) now carries two other arrangements: `split`, in which the outward burn's 1.8 tonnes of propellant for each tonne of cargo are bought on Earth and launched with it, and a schedule that fuels from Earth until year 20 and from the depot afterwards. The table runs the diversified and services portfolios on constant growth at three launch prices held through the programme. Each cell gives the depot throughout, `split` throughout and `split` until year 20 in that order; money is in trillions of 2041 dollars with credit in U and a 3% discount rate.

| Portfolio | Launch price, dollars of 2026 a kg | Funding at the opening | Funders' net present value, year 200 | Equity added by year 200 | Launches from Earth in year 200 |
|---|---|---|---|---|---|
| Diversified | 200 | 4.24 · 1.99 · 1.99 | 9.04 · 7.91 · 10.80 | 0.00 · 0.00 · 0.00 | 61,017 · 124,900 · 61,017 |
| Diversified | 100 | 3.96 · 1.69 · 1.69 | 10.22 · 11.67 · 12.04 | 0.00 · 0.00 · 0.00 | 61,017 · 124,900 · 61,017 |
| Diversified | 50 | 3.82 · 1.54 · 1.54 | 10.80 · 12.81 · 12.66 | 0.00 · 0.00 · 0.00 | 61,017 · 124,900 · 61,017 |
| Worldwide services | 200 | 4.09 · 1.91 · 1.91 | −3.70 · −2.06 · −1.72 | 7.35 · 97.87 · 7.77 | 19,975 · 55,931 · 19,975 |
| Worldwide services | 100 | 3.82 · 1.62 · 1.62 | 2.42 · −1.46 · 4.21 | 0.00 · 5.15 · 0.00 | 19,975 · 55,931 · 19,975 |
| Worldwide services | 50 | 3.69 · 1.48 · 1.48 | 5.45 · 6.75 · 7.27 | 0.00 · 0.00 · 0.00 | 19,975 · 55,931 · 19,975 |

Fuelling the first years from Earth halves the funding the programme needs before it earns anything, from 4.2 to 2.0 trillion dollars at the staged launch price, because the founding deployment no longer includes the mines and tugs that would fuel the ferries. Moving to the depot's propellant at year 20 then keeps the mature years as they were, and the funders' net present value rises by 1.8 to 2.0 trillion dollars in both portfolios at every price tested. Staying on Earth's propellant throughout is the better course only at 50 dollars a kilogram in the diversified portfolio, and by 0.15 trillion: it doubles the launches, to 125,000 a year at year 200, and at 200 dollars a kilogram it costs the services portfolio 98 trillion dollars of added equity.

The launch price itself moves the result more than the arrangement does. At 100 dollars a kilogram in place of 200 the services portfolio finances itself, and its funders' net present value rises from −3.7 to 2.4 trillion dollars before any change of arrangement.

The traffic these cases imply is part of the result. In the diversified portfolio with the depot's propellant, year 200 has 61,000 launches from Earth, which is 305 pads at the model's 200 launches a pad a year, with 259,000 lunar landings, 25,700 ferry trips and 10,800 returns of tugs from the sources. In year 100 it has 5,200 launches and 13,700 landings.

## Working parts made in place

The staged model buys the working parts of every new machine from Earth in every year. The [industry model](../../../provisioning/vacuum_industry/README.md) now lets the industry make a share of them, as far as its residents can staff the work. In these cases the industry becomes able to make 78% of working parts over the thirty years from the opening, which leaves 5% of all equipment mass imported, the central value of the published studies. The labour of a locally made part is this study's assumption, since no source gives it, and runs at 50%, 25% and 10% of the part's import price. Both portfolios are on constant growth, at the staged launch price and Earth-leg arrangement; money is in trillions of 2041 dollars with credit in U and a 3% discount rate.

| Portfolio | Working parts | Share made in place, orders of years 100 · 200 | Funders' net present value, year 200 | Equity added by year 200 | Costs paid to Earth, years 151–200 | Launches from Earth in year 200 |
|---|---|---|---|---|---|---|
| Diversified | All bought, as staged | 0% · 0% | 9.04 | 0.00 | 62% | 61,017 |
| Diversified | Labour at 50% of the import price | 10% · 38% | 9.44 | 0.00 | 58% | 50,727 |
| Diversified | Labour at 25% | 20% · 75% | 9.95 | 0.00 | 56% | 45,111 |
| Diversified | Labour at 10% | 51% · 78% | 11.27 | 0.00 | 51% | 44,339 |
| Worldwide services | All bought, as staged | 0% · 0% | −3.70 | 7.35 | 42% | 19,975 |
| Worldwide services | Labour at 50% of the import price | 33% · 78% | 5.80 | 0.00 | 29% | 7,673 |
| Worldwide services | Labour at 25% | 65% · 78% | 7.36 | 0.00 | 17% | 7,599 |
| Worldwide services | Labour at 10% | 78% · 78% | 7.71 | 0.00 | 17% | 7,598 |

The residents are what limits the share. With labour at half the import price, the diversified portfolio's people can staff a tenth of their working parts in year 100, with nine tenths of the workforce employed, and 38% in year 200. At a quarter of the import price they staff a fifth and then three quarters; at a tenth, half and then the whole able share. The services portfolio employs fewer people in its other work and makes more of its machinery: a third of its working parts in year 100 even at the dearest labour.

That portfolio is the one the change transforms. Bought machinery, its replacement and their launch are 42% of what the services portfolio spends in years 151 to 200, and with every part bought it needs 7 trillion dollars of added equity on constant growth and has a net present value of −3.7 trillion. Making its working parts in place, it finances itself at every labour cost tried, with a net present value of 5.8 to 7.7 trillion dollars, and its launches in year 200 fall from 20,000 to about 7,600. The diversified portfolio gains 0.4 to 2.2 trillion dollars, and its payments to Earth stay above half of its costs. About a third of those payments is the cost of selling material to Earth, which local manufacture leaves as it was: operating the capsules that return the goods costs 82 trillion dollars in years 151 to 200, and relaunching them is a third of all the mass launched. The rest is machinery that its residents, busier than the services portfolio's, staff more slowly.

## Demand judged at its own date

[scenarios.json](scenarios.json) states three scenarios of demand, `low`, `central` and `high`, and [demand.py](demand.py) turns one into the orders of every year. Each value in the file carries its basis: sourced, with the lookup it was read from, or constructed, where no opened source gives a figure and the study states its own. The lookups are on the research drive.

The customers live in the world of the year in question. Its population follows the United Nations' medium projection, 9.66 billion in 2050 and 10.18 billion in 2100, and is held there after 2100, where the projection ends. Income per person grows 1.6% a year to 2049, as projected, and 0.5%, 1.0% or 1.6% a year after. World output follows the two from a projected 163 trillion dollars of 2015 in 2050. No quantity is a multiple of the founding system, and none depends on how many residents the settlement has.

| Line of sale | What sets it | Low | Central | High |
|---|---|---|---|---|
| Sport and entertainment | Followers at maturity, as a share of the world's people | 0.92% | 3.7% | 10% |
| | Receipts per follower a year (dollars of 2024), rising with income | 45 | 45 | 130 |
| | Growth of receipts from the opening portfolio of 1.2 billion a year | 6% | 10% | 19% |
| | Share others keep, production and administration included | 85% | 65% | 43% |
| Hosted research and professional work | Share of world output (billions of 2041 dollars in 2049) | 0.0015% (5) | 0.005% (18) | 0.02% (71) |
| | Share others keep; what the settlement keeps per resident job (millions) | 85%; 1 | 72%; 2 | 65%; 5 |
| Computing | Capacity wanted in 2049 at the price of the same service on Earth | 1 GW | 5 GW | 20 GW |
| Propellant at the depot | Tonnes a year in 2049 at 500 dollars a kilogram | 60 | 450 | 1,640 |
| | Elasticity to price; growth of activity in space to 2100 | 0.6; 3% | 1.0; 5% | 1.4; 9% |
| Structural material at the depot | A tonne for each tonne of propellant, on the same terms | | | |
| Precious metals to Earth | Supply as a share of the world's platinum-group production | none | 50% | 200% |
| | Cost of treating a kilogram of feed (dollars of 2041) | 1.00 | 0.10 | 0.01 |
| Visits | Elasticity to the price of a visit; to income | 0.5; 1.0 | 1.25; 1.5 | 2.0; 2.0 |
| Launch from Earth | Dollars of 2026 a kilogram to low orbit | 273 | 100 | 33 |
| The Earth leg's propellant | Launched from Earth until year 20, from the depot after | | | |
| Working parts made in place | Share the industry is able to make, and the years it takes | 57% in 30 | 78% in 20 | 96% in 12 |
| | Resident labour as a share of the part's import price | 50% | 25% | 10% |
| Technology from the opening | Fall in the price of machinery's working parts for each doubling of those acquired | 5% | 10% | 15% |
| | Growth a year of what a job does, the wage rising with it | 0.5% | 1.0% | 1.6% |

A follower is a person who watches the competitions or programmes regularly, and it is the one audience count in the scenarios. The opened evidence gives 45 dollars a year for a regular viewer of a leading world circuit, against 2 to 4 for a claimed fan and 130 for a paying subscriber. The [humanities study](../../../humanities/economics/README.md) counts distinct fans, of whom one in five pays; its portfolio of 1.2 billion U a year from 140 million distinct fans is the opening point here, with U read as a dollar of 2041. Receipts grow from it at the rate found for properties that succeeded until they meet the mature following, and follow population and income after.

Propellant and material are handed over at the depot. Their price has a ceiling: 80% of what the same kilogram costs launched from Earth to the depot in that year, with the propellant that carries it. The anchor for visits is a published model of 2022, about 1,370 travellers a year to orbit at 4 million dollars a seat, carried forward with the number of wealthy buyers. Each visitor brings 270 kg launched from Earth, a tonne of cabin on the ferry and the lander each way, and 14 nights.

Three lines rest on values the study constructed, and the tables below should be read with that in mind. No opened source gives a market for structural material in space, and the line is tied to propellant tonne for tonne. No source gives the computing capacity customers would contract, and the first pass's three figures are carried; with technology improving, computing becomes the central scenario's largest line by year 200, so this construction weighs on the central results. The size of the hosted work is set around the one published scenario of lunar surface activity to 2050. The high scenario compounds its upper elasticities for two centuries.

Technology improves from the opening in every scenario. The price of machinery's working parts falls by 5%, 10% or 15% for each doubling of the working parts acquired since the founding deployment: the rates published for satellite programmes of 1 to 10, 11 to 50 and more than 50 units, where across more than a hundred manufacturing processes the fall runs from 46% to a rise of 8%. An order pays the price of its own year. What a job does grows at each scenario's growth of income per person on Earth, a link the study constructs because no source gives the productivity of an industry beyond Earth; the wage rises with it, so the same work costs the same and needs fewer of the residents' hands. The launch price keeps each scenario's one value for every year, and chemical freight, the masses each flow needs, the lifetimes, the settlement's equipment per resident and the prices of operating materials keep the founding design's. The case `<scenario>:constant_growth:frozen_technology` runs each scenario with the founding design throughout.

## What the system sells, at what price and how far

[mix.py](mix.py) turns the demand into orders the system can fill and is better off filling. It first sizes the industry for the residents alone and for every line of sale, and cuts every line alike in the years the residents could not staff them all. It then accounts each line by sizing the industry again without it: the difference is what the line costs, with the equipment, the freight and the propellant it calls for, and the jobs it needs.

The system sets the price of propellant, material and visits. Where buyers answer a price cut with more than a matching rise in quantity, the price that brings the seller most is the cost of a unit times e/(e−1), where e is the elasticity; where they do not, the seller holds the ceiling, or for visits the price at which the travel model is anchored. A line that costs more than it brings over the programme, at a 3% discount, is tried with a later start, at the first year of each later period, and dropped where no start pays. The residents' free hands then go to the lines kept, year by year, in the order of what a job in each brings, and to making the working parts of machinery after that. The industry is sized for those orders, and each served line is accounted once more against them.

The prices are the seller's. The case `central:constant_growth:cost_plus` prices the three lines at cost with a margin of half, which carries more travellers and sells more propellant for less profit.

## The three scenarios

Each scenario runs on both population ramps, with its technology improving from the opening. Money is in trillions of 2041 dollars a year or in total, with credit in U and a 3% discount rate; the first pass's diversified portfolio on constant growth stands beside them.

| Scenario and ramp | Receipts a year, year 100 · 200 | Costs a year, years 51–100 · 151–200 | Costs paid to Earth, years 151–200 | Funding at the opening | Equity added by year 200 | Funders' net present value, year 100 · 200 |
|---|---|---|---|---|---|---|
| First pass, diversified, constant growth | 3.78 · 13.29 | 1.86 · 10.45 | 62% | 4.24 | 0.00 | 6.78 · 9.04 |
| Low, constant growth | 0.06 · 0.15 | 0.11 · 1.79 | 60% | 1.86 | 103.89 | −1.99 · −2.81 |
| Low, front-loaded | 0.09 · 0.15 | 0.32 · 1.80 | 62% | 1.98 | 138.10 | −6.77 · −8.05 |
| Central, constant growth | 0.58 · 1.87 | 0.28 · 1.32 | 31% | 1.64 | 0.00 | −0.48 · −0.08 |
| Central, front-loaded | 0.61 · 1.87 | 0.34 · 1.35 | 32% | 1.73 | 2.76 | −2.03 · −1.82 |
| High, constant growth | 7.44 · 153.59 | 1.97 · 26.77 | 15% | 1.84 | 0.00 | 11.16 · 37.94 |
| High, front-loaded | 10.73 · 153.65 | 2.18 · 26.80 | 15% | 1.94 | 2.75 | 22.22 · 52.92 |

Judged at its own date, the central demand brings 0.58 trillion dollars a year in year 100 and 1.87 trillion in year 200, about a seventh of the first pass's receipts in each of the two years. On constant growth the system covers what it spends in every period after its first: 83 billion dollars a year of receipts against 90 of costs in years 8 to 25, then 227 against 193 in years 26 to 50, 449 against 279 in years 51 to 100, 803 against 514 in years 101 to 150 and 1,435 against 1,323 in years 151 to 200, while the residents go from 250,000 to 6.47 million. The planned settlement needs no equity after its opening funding of 1.64 trillion dollars, and the funders, paid 27 trillion by year 200, come within 0.08 trillion of recovering the opening in present value at 3%. The front-loaded ramp spends earlier on residents the central demand does not yet need: it adds 2.8 trillion dollars of equity, and its net present value is 1.7 trillion lower.

The low scenario's receipts, 0.15 trillion dollars a year at year 200, do not cover the planned settlement's second century, which needs 104 trillion dollars of added equity on constant growth. The high scenario finances the whole programme with a net present value of 38 trillion dollars on constant growth; 126 trillion of its 154 trillion dollars of receipts in year 200 are 238 million visits a year at half a million dollars each, the product of its upper elasticities to price and income held for two centuries. The first pass's receipts lie between central's and high's from year 51 on; in the first fifty years, when the first pass already sells a share of a founding system's exports, they exceed high's.

The lines of sale of the central scenario on constant growth, with prices and costs in 2041 dollars:

| Line | Year 100: quantity · price · receipts a year | Year 200: quantity · price · receipts a year | Cost of a unit, years 51–100 · 151–200 | Worth at 3% |
|---|---|---|---|---|
| Visits | 51,000 visitors · 2.2 million · 113 billion | 327,400 · 1.7 million · 541 billion | 522,000 · 339,000 | 980 billion |
| Computing | 13.2 GW · 22 a watt-year · 290 billion | 35.8 GW · 22 · 784 billion | 14.9 · 7.9 a watt-year | 384 billion |
| Sport and entertainment | 377 million followers · 276 · 104 billion | 377 million · 747 · 282 billion | 147 · 388 | 298 billion |
| Hosted research and professional work | 2,610 resident jobs · 17.8 million · 47 billion | 2,610 · 48.3 million · 126 billion | 12.1 million · 29.1 million | 80 billion |
| Propellant at the depot | 13,200 t · 480 a kilogram · 6 billion | 35,800 t · 480 · 17 billion | 15 · 9 a kilogram | 57 billion |
| Material at the depot | 13,200 t · 480 a kilogram · 6 billion | 35,800 t · 480 · 17 billion | 25 · 17 a kilogram | 55 billion |
| Precious metals to Earth | 269 t · 55,700 a kilogram · 15 billion | 2,290 t · 44,500 · 102 billion | 46,300 · 46,200 a kilogram | 24 billion |

Computing is the line that technology opens. At the price of the same service on Earth, 22 dollars a watt-year, it does not repay its hardware, its launch and its power with the founding design's prices, and the founding-design cases drop it in every scenario; with the price of its working parts falling as the industry acquires more of them, its cost falls from 15 dollars a watt-year in years 51 to 100 to 8 in years 151 to 200, and it becomes the largest line in year 200, 784 billion dollars of the 1.87 trillion. The first pass sold 8% of its year-200 receipts in computing. Propellant made from the sources' water costs the system 9 to 15 dollars a kilogram at the depot once the industry is running, against 600 for a kilogram launched from Earth at a launch price of 100 dollars a kilogram to low orbit; that gap is the margin on every tonne sold, and the published demand it is sold into is small, 13,200 tonnes in year 100. A visit costs the system 0.3 to 0.5 million dollars to carry and house, and with buyers who answer price as the central scenario states, the seller's price is four to five times that. The hosted work brings the same receipts with 2,610 residents' jobs where the founding design needed 6,500 in year 100 and 17,600 in year 200, because each job does more.

## The Earth leg, local manufacture, the launch price and the price of travel

Four cases change one thing in the central scenario on constant growth. They repeat on this demand the comparisons first made on the first pass's.

| Case | What changes | Funding at the opening | Equity added by year 200 | Costs paid to Earth, years 151–200 | Funders' net present value, year 200 |
|---|---|---|---|---|---|
| Central, constant growth | | 1.64 | 0.00 | 31% | −0.08 |
| `depot_throughout` | The ferries are fuelled at the depot from the opening, as in the staged model | 3.60 | 0.00 | 31% | −2.00 |
| `parts_bought` | Every working part is bought from Earth | 1.64 | 45.68 | 58% | −1.70 |
| `launch_200` | Launch costs 200 dollars a kilogram, as in the staged model | 1.63 | 0.00 | 38% | −0.56 |
| `cost_plus` | Propellant, material and visits are priced at cost and a half | 1.82 | 0.13 | 32% | −0.67 |

Fuelling the ferries' outward burn from Earth for the first twenty years takes 2.0 trillion dollars off the funding at the opening. The founding deployment is imported whole, its freight burns propellant at the depot, and the staged arrangement buys on Earth the tugs that fetch that propellant's water before anything is produced: for the residents alone, 570 tugs of 1,200 tonnes each, 1.6 trillion dollars of a founding cost of 2.3 trillion. Making working parts in place leaves the opening as it is and takes 46 trillion dollars off the equity the second century needs, because the growing settlement then buys 31% of what it spends from Earth, against 58% with every part bought. A launch price of 200 dollars a kilogram needs no equity either and costs the funders 0.5 trillion dollars of net present value. Pricing travel at cost and a half carries 229,000 visitors in year 100 against 51,000, at 0.67 million dollars against 2.2, and sells 258,000 tonnes of propellant against 13,200; it needs 0.13 trillion dollars of equity and costs the funders 0.6 trillion of net present value.

## What improving technology changes

The cases `<scenario>:constant_growth:frozen_technology` hold each scenario's technology at the founding design throughout. Beside each stands the same scenario with its technology improving. Money is in trillions of 2041 dollars, with credit in U and a 3% discount rate.

| Scenario and technology | Receipts a year, year 100 · 200 | Costs a year, years 51–100 · 151–200 | Costs paid to Earth, years 151–200 | Funding at the opening | Equity added by year 200 | Funders' net present value, year 100 · 200 | Working parts made in place, orders of years 100 · 200 |
|---|---|---|---|---|---|---|---|
| Low, founding design | 0.06 · 0.10 | 0.14 · 2.48 | 66% | 1.86 | 151.84 | −2.18 · −3.41 | 37% · 57% |
| Low, improving | 0.06 · 0.15 | 0.11 · 1.79 | 60% | 1.86 | 103.89 | −1.99 · −2.81 | 57% · 57% |
| Central, founding design | 0.28 · 1.02 | 0.17 · 1.49 | 39% | 1.34 | 37.15 | −1.21 · −1.40 | 61% · 78% |
| Central, improving | 0.58 · 1.87 | 0.28 · 1.32 | 31% | 1.64 | 0.00 | −0.48 · −0.08 | 78% · 78% |
| High, founding design | 4.22 · 56.65 | 1.04 · 16.79 | 50% | 1.51 | 0.00 | 8.07 · 17.17 | 3% · 10% |
| High, improving | 7.44 · 153.59 | 1.97 · 26.77 | 15% | 1.84 | 0.00 | 11.16 · 37.94 | 12% · 96% |

The two improvements act on different parts of the account. Cheaper working parts lower what the growing industry buys and open the line of computing, which pays from the opening in the central and high scenarios and from year 101 in the low one; the opening funding rises by 0.3 trillion dollars in central and high because the computing served from year 8 is part of the first orders. More productive work frees the residents' hands. In year 200 the central settlement's industry and lines employ 63,000 of its residents against 341,000 with the founding design, 2% of its workforce against 13%, and the high scenario's employ 551,000 against 2.5 million, 21% of its workforce against 96%; with those hands the high scenario makes 96% of its working parts in place by year 200 against 10%, and pays 15% of its costs to Earth against 50%.

Together they carry the central scenario's planned settlement without equity after the opening, where the founding design needed 37 trillion dollars: receipts nearly double, from computing and from visits that cost less to carry and house, and costs in years 151 to 200 fall from 1.49 to 1.32 trillion dollars a year even with computing's hardware added. In years 201 to 220 the planned 6.47 million residents cost 1.67 trillion dollars a year against receipts of 2.11, where with the founding design they cost 1.69 against 1.21. The low scenario still needs 104 trillion dollars of equity; its customers buy too little for the improvement it states to cover a settlement of millions. The high scenario's net present value at year 200 more than doubles, to 38 trillion dollars, and what it draws from its sources grows faster still: on constant growth it draws 383 billion tonnes of contained water by year 200 against 59 with the founding design, and it passes the median of the source survey's wider near-Earth screen in year 134 against 188, because its sales of propellant and of visits, whose ferries burn the sources' water, grow faster.

Technology as stated here makes machinery cheaper and work more productive, and leaves the freight as it was: chemical, on tugs that take three years for a round trip and burn more propellant than they deliver. It brings no farther source within reach. Every case now records its water against two pools of the source survey: the wider near-Earth screen, 39 billion tonnes of contained water at its median, and the broader near-Earth population without Don Quixote, 961 billion tonnes at its median, whose use would need transport that reaches beyond the window of the programme's sources (semimajor axis to 1.8 AU, eccentricity to 0.5, inclination to 30 degrees). Reaching those sources is a question of freight: of tugs chosen by what their freight costs, electric among them, and of how far cheaper freight widens the window.

## What the mined feed makes available

Every case records what the feed mined at the sources yields, beside what the system uses and sells of it. In the central scenario on constant growth the sources yield 84.3 trillion kg of feed by year 200, and water sets that quantity in every year.

| From the feed by year 200 | Recoverable | Used by the system | Sold | Left |
|---|---|---|---|---|
| Water | 5.90 trillion kg | 5.90 trillion kg | 2.9 billion kg of it as propellant | none |
| Metal | 4.05 trillion kg | 0.13 trillion kg | 2.9 billion kg as material at the depot | 3.92 trillion kg, 97% |
| Nitrogen | 121 billion kg | 23 billion kg | | 98 billion kg |
| Platinum-group metals and gold | 348,000 t contained | | 133,000 t | 215,000 t in the tailings |

The system mines for water, uses all of it and sells a two-thousandth of it. The metal and the nitrogen that the same rock holds are worked only as far as the system's own structures and air need them, and no scenario has a buyer for the rest. Nothing is held back: the protection gates of the source-access study reserve no quantity yet.

## The residents as a result

The programme's ramps say how many residents are housed in each year, and the run shows whether the system pays for them. [residents.py](residents.py) asks the question the other way round: how many residents each scenario of demand can carry, and when. From the opening's 12,500 residents in year 8 the settlement grows through twelve periods in turn (years 9–12, 13–18, 19–25, 26–37, 38–50, 51–65, 66–80, 81–100, 101–125, 126–150, 151–175 and 176–200), each at the highest steady rate at which, with no growth after that period, three conditions hold to year 220:

- **Money.** The system needs no equity after its opening funding, with credit in U, while it keeps its surplus for growth. The opening funding may be no larger than the constant-growth ramp's: the settlement's construction is ordered up to eight years ahead, so the housing of the first newcomers is bought before opening, and growing faster than the plan in those years would need more money from the funders.
- **Hands.** In every year the residents' own work, for themselves, for their newcomers' housing and supply and for the industry these need, fits within the workforce. The lines of sale take the hands left, and the making of working parts takes what the lines leave.
- **Income.** What the system holds does not fall over years 201 to 220, so the last population is carried by what the system earns and not by what it saved.

```
python3 -m research.studies.vacuum_economy.residents [--scenario NAME ...] [--jobs N]
```

A period's rate is found to within 5% of itself in six to eight trials, each of which sizes the industry twice or more and finances it, in 25 seconds to a little over a minute; the three scenarios take about an hour and ten minutes on three processes, and the search keeps its progress in research/runs so that an interrupted solution resumes. The lines of sale, their prices and their first years are those [mix.py](mix.py) chooses on the constant-growth ramp. The run then chooses them again on each solved schedule and finances it under the study's ordinary rules, which pay the funders whatever the remaining schedule does not need; that is the check on the search. The ramp `resident_led_beyond_plan` is the schedule found, and `resident_led` is the same schedule held at the planned population from the year it reaches it.

The solved schedules beside the constant-growth ramp, with each scenario's technology improving, financed by the run with credit in U; money in trillions of 2041 dollars at a 3% discount rate. The low scenario never reaches the plan, so its two solved cases are one schedule.

| Scenario and ramp | Residents in year 50 · 100 · 200 | 250,000 residents | 6.47 million residents | Receipts a year, year 100 · 200 | Equity added by year 200 | Funders' net present value, year 100 · 200 |
|---|---|---|---|---|---|---|
| Low, constant growth | 49,000 · 249,600 · 6.47 million | year 101 | year 200 | 0.06 · 0.15 | 103.89 | −1.99 · −2.81 |
| Low, residents as a result | 38,100 · 79,800 · 194,300 | not reached | not reached | 0.06 · 0.15 | 0.01 | −1.69 · −1.69 |
| Central, constant growth | 49,000 · 249,600 · 6.47 million | year 101 | year 200 | 0.58 · 1.87 | 0.00 | −0.48 · −0.08 |
| Central, residents as a result, held at the plan | 70,500 · 1.57 million · 6.47 million | year 73 | year 197 | 0.61 · 1.87 | 0.00 | −1.48 · −1.46 |
| Central, residents as a result, beyond the plan | 70,500 · 1.57 million · 6.86 million | year 73 | year 197 | 0.61 · 1.87 | 0.00 | −1.48 · −1.47 |
| High, constant growth | 49,000 · 249,600 · 6.47 million | year 101 | year 200 | 7.44 · 153.59 | 0.00 | 11.16 · 37.94 |
| High, residents as a result, held at the plan | 539,600 · 6.47 million · 6.47 million | year 46 | year 62 | 10.72 · 154.00 | 5.50 | 14.26 · 44.60 |
| High, residents as a result, beyond the plan | 539,600 · 56.7 million · 959 million | year 46 | year 62 | 10.65 · 156.15 | 6.59 | −0.04 · 1.30 |

The central scenario now carries the plan. It reaches the founding system's 250,000 residents in year 73, a million in year 92 and the planned 6.47 million in year 197; with the founding design it carried 2.15 million by year 200. Its growth is uneven. It grows 4.9% and 7.0% a year to year 18 on the opening funding and almost not at all in years 19 to 37; in years 8 to 25 the computing it serves from the opening costs 1.33 trillion dollars and brings 0.87. It then grows 9.0% a year in years 38 to 50, 2.5% in years 51 to 65, 12.0% in years 66 to 80 and 5.3% in years 81 to 100, and 0.9% to 1.8% a year after. The low scenario carries 194,300 residents by year 200, 3% of the plan, against 141,300 with the founding design: about 5% a year to year 18, almost nothing for the twenty years after, and then between almost nothing and 4.6% a year as its receipts allow. In both, the opening funding sets the pace of the first period and money sets it after. The residents' hands bind in no period, where with the founding design they held central's growth in years 38 to 50. Financed by the run under the ordinary rules, with the sales chosen again on each schedule, the solved schedules of both scenarios need at most a hundredth of a trillion dollars of equity after the opening. The funders do not recover the opening on these schedules: they lose 1.69 trillion dollars in present value in the low scenario and 1.46 in the central one, where the constant-growth ramp, slower in its first century, loses them only 0.08.

The receipts barely move with the residents. The central scenario sells 1.87 trillion dollars a year in year 200 on its solved schedule and on constant growth alike, because its customers' demand is filled before its hands run out. With improving technology the planned settlement pays its way after the plan is reached: in years 201 to 220 it spends 1.67 to 1.72 trillion dollars a year against receipts of 2.11, where with the founding design it spent 1.69 against 1.21. A scenario carries as many residents as the surplus of its sales pays the upkeep of.

The high scenario reaches the planned population in year 62, and the same money would carry 959 million residents by year 200, 148 times the plan; with the founding design it reached the plan in year 78 and would carry 264 million. Those schedules rest on the sales the system chooses on the constant-growth ramp. On its own schedules the system serves the same seven lines from year 8, and prices them from the costs of a larger industry: a settlement of 6.47 million by year 62 makes visits and propellant more cheaply, so the system sells more of them at lower prices. In year 100 it sells 9.5 million visits at 0.54 million dollars each, against 2.8 million at 1.0 million on the ramp, and propellant at 23 dollars a kilogram against 58. Sold at those prices, the lines pay over the programme at 3%, which is how mix.py judges them, but cost more than they bring in their first years; the run needs 5.5 trillion dollars of equity before year 100 on the schedule held at the plan and 6.6 trillion beyond it, where the search, with the ramp's prices, needed none. With the founding design the gap was 29 to 31 trillion. In the high scenario the rule that grows the residents and the rule that chooses the sales still disagree, by less.

The high scenario's growth beyond the plan outruns the water of its sources sooner with improving technology. On the run's sales its schedule draws the median contained water of the wider near-Earth screen, 39 billion tonnes, by year 88, when it has 43 million residents, and passes the screen's 5th and 95th percentiles in years 79 and 104; it passes the broader near-Earth population's 5th percentile, 504 billion tonnes, in year 163 and its median, 961 billion tonnes, in year 185, with 621 million residents. Held at the plan, the high scenario passes the wider screen's median in year 114, and on constant growth in year 134. Low and central draw 8 billion tonnes by year 200 on constant growth; on their solved schedules low draws 0.7 billion tonnes and central 13.8 billion, a third of the wider screen's median. Every case records these comparisons. The pools' water is the survey's conditional range before recovery and before any limit on access, and it does not stop growth here, since the sales draw on the same water and the study does not yet divide it between them and the residents.

The rule grows each period as fast as it can given the periods before it, so its schedules are uneven, as central's pause in years 19 to 37 and its bursts of 9% and 12% a year after it show. A rule that spread growth evenly, or one that sought the most residents by a stated year, would draw other schedules from the same accounts.

## What these results rest on, and what is not yet in them

The scenario values are stated one by one in scenarios.json, and several are the study's own. On the planned ramps the residents follow the plan whatever the demand, so a scenario that cannot pay for them shows as added equity; on the solved ramps they follow what the system's money and hands sustain, with the sales chosen on the constant-growth ramp, and the cost of each resident's housing and supply is the founding system's, per head. Technology improves in two ways only, the price of working parts with the working parts acquired and what a job does with the years, at rates the scenarios state; the productivity rates follow income on Earth by the study's construction, and the two improvements are not yet run apart. All freight is chemical, on tugs that take three years for a round trip and burn more propellant than they deliver, and the Earth leg changes arrangement in year 20 because that year was tried, with no search for the best. Power is beamed from the constellation the staged model sizes. The residents' hands go to the lines of sale before the making of working parts; in the high scenario the lines leave that manufacture 12% of working parts in year 100, and 96% in year 200 once productivity has freed the hands. A line is served from the opening or from the first year of a later period, or dropped. Each of these is a step still to take, and each moves the figures above.

## Evidence

The records are conditional accounts. The industry is a numerically closed sizing of equipment, freight and power for the orders it is given; the orders are scenarios of demand; the financing follows stated rules. A converged case shows that the stated machinery fills the stated orders within the model, and the qualification of the sources, the plant, the markets and the prices lies outside it. The reproduction of the staged appraisal is exact equality of 22 measures at each of 729 checkpoints, checked for all 81 staged order books on 15 November 2026; the tests repeat it for the nine cases on the staged schedule against the saved product and solve one case afresh. The scenario cases are checked against the saved product: each serves only lines that pay and that its residents can staff, its receipts are the receipts of its lines, and its resource account uses no more than the feed yields. They are not solved afresh in the tests, since one takes five to ten minutes. When technology began to improve, on 16 November 2026, the three cases that hold it at the founding design reproduced the earlier constant-growth records value for value, and the 49 first-pass cases kept every value; the learning index and the growth of productivity are tested in the industry model, with the quantities they leave as they were. The solved schedules are checked against [results/residents.json](results/residents.json): each follows the rates of its periods, opens on the opening population and never falls, the schedule held at the plan holds from the year the plan is reached, and the product names the digest of every file it was computed from. The search itself is checked on a stand-in trial whose highest sustained rates are known. A schedule is not sought afresh in the tests, since the three take about an hour and ten minutes.
