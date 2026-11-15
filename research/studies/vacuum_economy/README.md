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

## The run

```
python3 -m research.studies.vacuum_economy.run [--case NAME ...]
```

The run takes the cases listed in [cases.json](cases.json) and writes [results/economy.json](results/economy.json), 1.2 MB for 49 cases. A case takes about ten seconds on one thread, and about twenty where the staffing of local manufacture is solved. Each record holds the funding at the opening; for each of three financing arrangements (equity alone, funded credit, and credit in U) and each horizon (years 100, 200 and 220) the receipts, the operating cash, the payouts, the equity added after the opening, the rate of return and the net present values at 0%, 1%, 3% and 7%; for each of seven periods the receipts by sector and the costs by category, with the share of the costs paid to Earth; the asteroid feed, the lunar regolith, the water, the nitrogen and the fresh metal drawn by each horizon; the traffic of each horizon's year, as launches from Earth with the pads they imply and the round trips of the ferries, the lunar landers and the tugs; and the jobs of that year with the share of the machinery's working parts made in place. A case may change engineering coefficients, each named by its path, to a value or to a schedule by year. The product names its producer and the digest of every file it was computed from.

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

## Evidence

The records are conditional accounts. The industry is a numerically closed sizing of equipment, freight and power for the orders it is given; the orders are scenarios of demand; the financing follows stated rules. A converged case shows that the stated machinery fills the stated orders within the model, and the qualification of the sources, the plant, the markets and the prices lies outside it. The reproduction of the staged appraisal is exact equality of 22 measures at each of 729 checkpoints, checked for all 81 staged order books on 15 November 2026; the tests repeat it for the nine cases on the staged schedule against the saved product and solve one case afresh.
