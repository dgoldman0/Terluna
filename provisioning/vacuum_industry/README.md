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

The coefficients hold several things fixed that bear on what the model can say about a growing industry. Launch to low Earth orbit costs 200 dollars a kilogram, imported hardware 500, spacecraft 1,500 and imported materials 150, in every year, unless the technology described below lowers the price of working parts. The working share of each asset is imported in every year. All freight is chemical, on hydrogen and oxygen from asteroid water at a specific impulse of 430 s: tugs of 10,000 t payload make a round trip to a near-Earth source every three years at 2.5 km/s each way, landers of 100 t carry cargo between the depot and the lunar surface at 2.8 km/s each way, and ferries of 200 t link low Earth orbit with the depot at 3.8 km/s each way, braking with their engines on the return and, unless the arrangement of the Earth leg says otherwise, fuelled for both burns at the depot. The surface's power is beamed from a constellation in orbit that is sized with the industry. The source rock yields water at a grade of 10% with 70% recovered, and the feed is the larger of the rock its water, its nitrogen and its metal each need.

## The Earth leg and coefficients that change with the year

Two coefficients may be given as a schedule, `{"by_year": [[year, value], ...]}`, in place of one value: the launch price `industrial.transport.LEO_USD2026_kg`, interpolated in a straight line between the listed years and held outside them, and the arrangement `industrial.transport.earth_leg`, which steps at each listed year. The founding deployment reads a schedule at year 4, the middle of its years.

`earth_leg` names where the ferry between low Earth orbit and the depot gets its propellant. Under `space`, the staged model's arrangement and the default, the depot supplies both burns from asteroid water: each returning ferry carries down the propellant of the next outward burn, so a tonne of cargo sent outward draws 6.2 tonnes of propellant from the depot. Under `split` the outward burn's propellant is bought on Earth and launched with the cargo, 1.8 tonnes for each tonne of cargo, and the depot fuels the return alone with 0.3 tonnes. Under `earth` both burns' propellant is launched, 2.8 tonnes for each tonne of cargo, as in the founding study's first arrangement. These figures are for cargo outward with an empty return. The ferry of the founding study carries tanks for the `space` arrangement; under the others it keeps its fixed mass and carries tanks for what it holds. The model reports the propellant launched for the Earth leg, counts it in the mass launched from Earth and prices it at the launch price plus its purchase on Earth.

## Working parts made in place

`economy.production.local_working_parts` states how far the industry makes the working parts of its own machinery. `share` is the share of a newly ordered machine's working parts that the industry is able to make in the order's year. The staged model's value is 0, and every working part is bought from Earth. The share made is fabricated from local metal as structure is, with the fabrication cells, the energy and the feed that takes, and it costs resident labour worth `labour_fraction_of_import_price` of the part's import price, paid at the industry's wage and counted among its jobs. The capsules that return goods to Earth are bought there in every case. Both coefficients may be schedules by year.

The residents limit the share. [industry.py](industry.py) solves the industry first with every working part bought, counts in each year the workforce left after the industry's and the settlement's other jobs, and staffs the able share as far as those hands reach in the order's year and the next. It then solves again and lowers the staffing of any year whose jobs still exceed the workforce. A caller may state the staffing itself, as `staffed`, and the industry is then solved once.

The published studies of industry beyond Earth give the imported remainder by mass: 10%, 5% or 1% of all equipment, reached within 12 to 30 years of the first plant, with microelectronics, bearings, instruments, motors and sensors the last to be made in place. Working parts are 23% of the founding industry's mass, so a remainder of 5% of equipment is a share of 78% of working parts. All of those figures are from paper studies. None gives a cost, and the labour fraction is this model's own assumption, to be carried as a range.

## Technology that improves with the years

`economy.technology` states how the industry's technology improves after the founding. Both of its rates are nought unless stated, which keeps the founding design throughout and the staged studies' records with it.

`learning_per_doubling` lowers the price of machinery's working parts by that fraction for each doubling of the working parts acquired. The founding deployment's working parts are the first unit, every later order adds its working parts to the total whether bought from Earth or made in place, and an order pays the price of its own year. Orders placed before `from_year`, the opening of production unless stated, keep the founding price. Structure, launch, operating materials and the charges for provision at the sources keep their stated prices. The quantities the model sizes do not depend on prices, so it sizes the industry once and then prices the same orders at the learned index.

`labour_productivity_growth` is the yearly growth, from `from_year`, of what a job does. Every count of jobs in the model falls by the accumulated factor: the industry's staff, the settlement's operators, the factory's and the computing staff, the makers of working parts, the visitors' staff and the workshops of refurbishment. The wage rises by the same factor, so the payroll for the same work stays as it was and the work needs fewer of the residents' hands.

## Goods handed over in space, and travellers

Four kinds of order go no further than the depot, and the staged order books state none of them. `space_source_goods` is material made at the sources from the metal of the mined feed and carried to the depot by the tugs. `space_goods` is goods made in the lunar factory and lifted to the depot by the landers; the factory and its staff are sized for these and the goods for Earth together. Neither kind is packed for return, rides a ferry to Earth or needs a capsule. `passenger_launch_kg_year` is the mass launched from Earth for travellers, charged at the launch price. `passenger_cabin_kg_year` is the cabins that carry them, which ride the ferries of the Earth leg and the landers in both directions.

## Evidence

The model's results are conditional engineering requirements. Convergence shows that the stated equipment carries the stated flows under the stated coefficients; it qualifies no plant, source, trajectory or price. The lifetime curve is a scheduling assumption with the stated mean. The recovery and refurbishment coefficients are hypotheses: the cited practice supports the mechanisms, and the industrial fractions, costs and restored service lives await demonstration.

This model is the staged studies' cohort model ([provisioning/enclosed_industry](../enclosed_industry/README.md)) with its recovery experiments made part of it and its statements written one to a line. With the staged arrangement of the Earth leg, the staged launch price and every working part bought, it returns the staged model's records value for value, beside the two outputs added since, for the staged order books compared: four books under five sets of measures before those options were added, and one book under the five sets after each. The coupled study's nine first-pass cases keep every figure. The staged model stays in its folder, unchanged, for the staged studies that read it.

## Tests

`provisioning/vacuum_industry/tests` checks the freight against the founding study's controls, the three arrangements of the Earth leg and the schedules by year, the goods handed over at the depot and the travellers, the locally made share of working parts with its labour, its staffing and its part in refurbishment, the learning index and the growth of labour productivity with the quantities they leave as they were, the lifetime curve's mean, the dating of orders, the feed as the largest of its three needs, the conservation of mass and money through refurbishment, the closure of every year's water and material balances on a small settlement, and that a run with no measures is the plain model.
