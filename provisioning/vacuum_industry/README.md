# The industry of the vacuum buildout

This model sizes the industry that houses the residents of the vacuum buildout and fills its customers' orders. It takes a schedule, year by year, of residents and of goods, propellant, computing, services and visitors' beds to supply, and it returns the equipment that schedule requires, the dated orders that put the equipment in place, the flows of material, freight, propellant and power, and their costs. The coupled [study of the buildout's economy](../../research/studies/vacuum_economy/README.md) supplies the schedule and reads the result.

## How it works

Residents set the pressure enclosures, the ordinary equipment and the utilities of the settlement. Everything else follows from work to be done: the mines, the process trains, the fabrication cells, the power plant, the radiators and the vehicles are each sized from the flow they must carry in that year, with a quarter of capacity spare.

Each asset has a structure, which the industry makes from lunar or asteroid material, and working parts, which it buys from Earth. A component is ordered two years before it enters service, half of its manufacture and freight falling in each year. It then survives by age on a Weibull curve of shape three whose mean is the component's stated life, and it retires without resale. Equipment that outlives its need stays available; equipment that is no longer needed is not replaced.

The flows close on themselves. Imported parts need freight; freight needs propellant; propellant, water, nitrogen and metal need feed from the near-Earth sources and regolith from the Moon; mining and processing need power and heat rejection; and all of these need equipment, whose orders return to the freight. [model.py](model.py) iterates this loop until the equipment required in every year stops rising, and it reports a year as feasible when the equipment in place covers the requirement and the jobs fit within 40% of the residents.

The founding deployment is sized once, for the eight years before production opens in year 8, and its equipment is imported whole. Its cost and its feed are reported apart from the operating years.

## Recovery and refurbishment

[recovery.py](recovery.py) holds two measures that the model calls at three points: after the ledger of cohorts, after a year's physical flows, and after a year's equipment requirement. [recovery.json](recovery.json) states their coefficients.

Closure raises the recovery of domestic and industrial water from 98% to 99.5%, of crop water from 99.5% to 99.7% and of fish-purge water from 95% to 98%. Household goods are collected at 98% and accepted at 97%, and retired structure at 99% and 99%. The sorting and treatment equipment and their energy are sized with the rest of the industry.

Refurbishment collects 80% of ten kinds of retired machine, among them the asteroid tugs, the landers, the mining and process trains and the fabrication cells. Inspection a year later accepts 60% of those collected. An accepted machine returns to service with new parts for 20% of its mass, in a workshop at the place it works, and the import it replaces falls by the other 80%. The workshops, their energy, consumables and labour are sized and paid for in the same balance.

[industry.py](industry.py) solves the model for one book of orders under a named set of measures (`none`, `closure`, `refurbishment` or `combined`) and returns the industry year by year together with an account of the loops: each water stream with its makeup, the retired structure and goods with what was collected, accepted and used, the machines withdrawn for repair, and the material drawn from the sources and the Moon.

## Coefficients

[engineering.json](engineering.json) holds the coefficients: the asset templates with their lives and structural shares, the yields and energies of the lunar and asteroid processes, the vehicles, and the prices. They come from the founding industrial study (commit `146a5582`) through its exporter, and they are today the same bytes the staged studies of 12 to 14 November 2026 read. Prices are in dollars of 2026 and are multiplied once by 1.03¹⁵ = 1.558 to give constant dollars of 2041.

The coefficients hold several things fixed that bear on what the model can say about a growing industry. Launch to low Earth orbit costs 200 dollars a kilogram, imported hardware 500, spacecraft 1,500 and imported materials 150, in every year. The working share of each asset is imported in every year. All freight is chemical, on hydrogen and oxygen from asteroid water at a specific impulse of 430 s: tugs of 10,000 t payload make a round trip to a near-Earth source every three years at 2.5 km/s each way, landers of 100 t carry cargo between the depot and the lunar surface at 2.8 km/s each way, and ferries of 200 t link low Earth orbit with the depot at 3.8 km/s each way, braking with their engines on the return. The surface's power is beamed from a constellation in orbit that is sized with the industry. The source rock yields water at a grade of 10% with 70% recovered, and the feed is the larger of the rock its water, its nitrogen and its metal each need.

## Evidence

The model's results are conditional engineering requirements. Convergence shows that the stated equipment carries the stated flows under the stated coefficients; it qualifies no plant, source, trajectory or price. The lifetime curve is a scheduling assumption with the stated mean. The recovery and refurbishment coefficients are hypotheses: the cited practice supports the mechanisms, and the industrial fractions, costs and restored service lives await demonstration.

This model is the staged studies' cohort model ([provisioning/enclosed_industry](../enclosed_industry/README.md)) with its recovery experiments made part of it and its statements written one to a line. For four of the staged order books under five sets of measures it returns the staged model's records value for value, and the coupled study's nine first-pass cases keep every figure. The staged model stays in its folder, unchanged, for the staged studies that read it.

## Tests

`provisioning/vacuum_industry/tests` checks the freight against the founding study's controls, the lifetime curve's mean, the dating of orders, the feed as the largest of its three needs, the conservation of mass and money through refurbishment, the closure of every year's water and material balances on a small settlement, and that a run with no measures is the plain model.
