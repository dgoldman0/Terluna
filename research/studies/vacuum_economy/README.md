# The economy of the vacuum buildout

This study follows the economy of the settlements beyond Earth through the vacuum period: from the opening of production in year 8 of the programme (2049) to the planned population of 6.47 million residents in year 200, and for twenty years past it. Its unit of account is U, the general economic unit of the whole planned system beyond Earth before the Open Moon, the array and the outer Solar System's buildout. For each case it draws up what customers order in each year, sizes the industry that fills those orders, finances that industry, and reports what the funders, the residents and the sources each give and receive.

The study couples parts that live in their own lanes and reads each where it is kept. The cohort model of industry is [provisioning/enclosed_industry/model.py](../../../provisioning/enclosed_industry/model.py), with the recovery and refurbishment measures of [reclamation.py](../../../provisioning/enclosed_industry/reclamation.py) beside it. The financing model is [enclosed_returns/finance.py](../enclosed_returns/finance.py). The run in this folder joins them for one case at a time.

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

The run takes the cases listed in [cases.json](cases.json) and writes [results/economy.json](results/economy.json), 205 kB for nine cases. A case takes about ten seconds on one thread. Each record holds the funding at the opening; for each of three financing arrangements (equity alone, funded credit, and credit in U) and each horizon (years 100, 200 and 220) the receipts, the operating cash, the payouts, the equity added after the opening, the rate of return and the net present values at 0%, 1%, 3% and 7%; for each of seven periods the receipts by sector and the costs by category, with the share of the costs paid to Earth; and the asteroid feed, the lunar regolith, the water, the nitrogen and the fresh metal drawn by each horizon. The product names its producer and the digest of every file it was computed from.

## The first pass

The scenario `first_pass` is the demand of the staged studies of 12 to 14 November 2026: their nine export portfolios at their strongest market, with real demand growing 0.5% a year, on the `first_pass` ramp, with material recovery and equipment refurbishment applied. Through this entry point every figure of the [staged appraisal](../enclosed_returns/revised.md) is reproduced exactly; the table repeats its headline, in trillions of 2041 dollars with credit in U and a 3% discount rate, and adds two columns the staged write-up does not carry.

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

## Evidence

The records are conditional accounts. The industry is a numerically closed sizing of equipment, freight and power for the orders it is given; the orders are scenarios of demand; the financing follows stated rules. A converged case shows that the stated machinery fills the stated orders within the model, and the qualification of the sources, the plant, the markets and the prices lies outside it. The reproduction of the staged appraisal is exact equality of 22 measures at each of 729 checkpoints, checked for all 81 staged order books on 15 November 2026; the tests repeat it for the nine listed cases against the saved product and solve one case afresh.
