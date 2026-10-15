# Geography and hydrology

**Current condition:** global LOLA topography (4 and 16 pixels per degree) and
the GRAIL GL0420A gravity field are fetched and hash-checked
([inputs.json](inputs.json), [fetch_inputs.py](fetch_inputs.py)). Heights are
referred to the geoid, and hydrostatic water storage is computed. The IAU lunar
nomenclature is fetched and hash-checked the same way and names the atlas's seas,
islands and targets. A first estimate of rivers and rain-fed lakes routes the
climate run's runoff over the terrain. The monthly tide of every sea is computed
from the lunar orbit and librations, with each sea's volume and the water's own
attraction, and checked against a dynamic solution. Erosion, groundwater and
climate-coupled shorelines are the next layer.

| File | Holds |
|---|---|
| [topography.py](topography.py) | LOLA radius grids; the geoid from degrees 2–200 of GL0420A plus Earth's static tide (×(1 + k2)) and rotation; heights above that equipotential; exact spherical cell areas |
| [basins.py](basins.py) | Level curves (area and volume below a common water level), separate seas, and the depression tree: every closed depression's spill level and capacity, from a union-find merge by height |
| [water_inventory.py](water_inventory.py) | Writes [results/](results/): the level table and the major basin joins as water rises, with schema, hashes and evidence statement |
| [atlas.py](atlas.py) | The atlas at the scenario water share: seas and lakes with depths and IAU water names, islands with summits, mare flooding and a comparison of candidate shares ([results/atlas.json](results/atlas.json); grid product in `products/`) |
| [coastal_grid.py](coastal_grid.py) | Native 16 px/deg LOLA terrain for the eastern Smythii wave sector, with the geoid evaluated directly and the atlas water level retained; writes the ignored `terluna.geography.coastal-grid/1` product ([wave experiment](../climate/waves/coastal.md)) |
| [shore_terrain.py](shore_terrain.py) | Native 256 px/deg LOLA terrain, about 118 m spacing, for the eastern Smythii shore and its surrounding sector; pinned archive rows and label in [shore_inputs.json](shore_inputs.json), direct GRAIL geoid, four-connected sea mask, and bulk products on the research drive ([shore exposure study](../climate/waves/shore.md)) |
| [coast_terrain.py](coast_terrain.py) | Native 256 px/deg LOLA terrain round the sea-appearance study's four rendering coasts: pinned row slices of three more LOLA tiles in [coast_inputs.json](coast_inputs.json) plus the Smythii rows, heights above the degree-200 geoid and the atlas sea level, resampled to ±72 km local grids (±69 km at Smythii) on the research drive |
| [shore_profiles.py](shore_profiles.py) | Cross-shore profiles along the local shallower-water normal at the Smythii stations; a versioned product separates interpolated 118 m rock information from explicit beach slopes for the [breaking and run-up study](../climate/waves/runup.md) |
| [nomenclature.py](nomenclature.py) | IAU lunar feature names from the USGS gazetteer archive (read by a small dBase reader built on the standard library) |
| [drainage.py](drainage.py) | Rivers and rain-fed lakes above sea level: runoff from the climate run routed over the 16 px/deg terrain, with fill-and-spill lakes set by each depression's water balance ([results/drainage.json](results/drainage.json); grid product in `products/`). `--climatology` takes another run's climate product, `--out` writes elsewhere |
| [groundwater.py](groundwater.py) | A first estimate of the water the porous crust takes up, from GRAIL's porosity, beside the seas and lakes ([results/groundwater.json](results/groundwater.json)) |
| [groundwater_uptake.py](groundwater_uptake.py) | How fast the crust takes the water: heat and water down columns under the seas, lakes and land, the cold airless crust thawing from above ([results/groundwater_uptake.json](results/groundwater_uptake.json); sources in [groundwater_sources.json](groundwater_sources.json)) |
| [groundwater_table.py](groundwater_table.py) | The water table on land: steady sideways flow over the atlas at 15 km with the seas and lakes held, springs where it reaches the ground, the mature and young Moon ([results/groundwater_table.json](results/groundwater_table.json)) |
| [lunar_ephemeris.py](lunar_ephemeris.py) | The Earth's and the Sun's directions over the Moon and their distances: the truncated ELP-2000/82 lunar theory with the optical librations, checked against JPL Horizons ([tides_ephemeris_check.csv](tides_ephemeris_check.csv)) |
| [Calendar ephemeris check](CALENDAR_EPHEMERIS.md) | Independent pinned JPL comparisons across 2000–2500 TT, with horizon timing and conditional clear-sky brightness error propagation |
| [tides.py](tides.py) | The monthly tide of every sea, 2026–2045: equilibrium height with each sea's volume held and the water's self-attraction, constituents, ranges and stations ([results/tides.json](results/tides.json); grid product in `products/`) |
| [tide_dynamics.py](tide_dynamics.py) | Linear shallow-water response of a sea to one tidal line, with the Moon's rotation and bed friction, and the sea's seiches |

## What the topography allows

These are hydrostatic storage geometries, not predicted shorelines. The global
layer is the water volume spread evenly over the whole Moon.

| Water-covered share of the Moon | Common level above the geoid | Water (global layer) | Water mass | Mean depth |
|---|---|---|---|---|
| 10% | −2.69 km | 95 m | 3.6×10¹⁸ kg | 0.95 km |
| 20% | −2.05 km | 188 m | 7.1×10¹⁸ kg | 0.94 km |
| 30% | −1.56 km | 311 m | 1.2×10¹⁹ kg | 1.04 km |
| 40% | −1.10 km | 470 m | 1.8×10¹⁹ kg | 1.18 km |
| 50% | −0.60 km | 696 m | 2.6×10¹⁹ kg | 1.39 km |
| 71% (Earth's ocean share) | +0.73 km | 1514 m | 5.7×10¹⁹ kg | 2.13 km |

For scale, the 1.2-atm atmosphere is about 2.8×10¹⁸ kg.

The nearside maria form one connected lowland:

- **About −2.5 km:** Imbrium, Frigoris and Serenitatis join.
- **By −1.7 km:** Procellarum, Nectaris, Fecunditatis, Nubium–Humorum and
  Crisium join. That sea covers 12% of the Moon and holds 71 m.
- **At −0.9 km:** it reaches 23% and 244 m.

South Pole–Aitken fills separately on the far side: 3.8% and 59 m at −2.7 km;
8% and 230 m at −0.3 km. At −0.3 km it joins the nearside ocean, and above the
geoid level the Moon has one world ocean.

## Atlas at the selected water share

`python -m geography.atlas` fills the Moon to one common level covering the scenario share in
[shared/scenarios/water.json](../shared/scenarios/water.json), 28%. The product
([results/atlas.json](results/atlas.json), schema `terluna.geography.atlas/1`) lists the separate seas and lakes, the
islands, how much of each mare floods, and the comparison of shares from 25% to 35% behind the choice. The
grid product `products/atlas_28pct_4ppd.npz` (kept out of Git, regenerated in seconds) holds heights, water and land
labels and 16 px/deg polar caps. The [atlas sheets](../visualization/atlas/) draw it, and the
[conservation study](../research/studies/conservation/) places heritage sites and scientific targets on it.

At 28%:

- **Sea level and water.** Sea level is −1,654 m above the geoid. The standing water equals a 282 m global layer, 1.07×10¹⁹ kg.
- **Land.** The main land mass covers 71.1% of the Moon and islands 0.9%.
- **The Near-side Sea** covers 13.9% of the Moon, with a mean depth of 648 m and a greatest depth of 3.75 km. It joins Oceanus Procellarum and Maria Imbrium, Serenitatis, Crisium, Nectaris, Fecunditatis, Nubium, Humorum, Cognitum and Frigoris.
- **The South Pole–Aitken Sea** covers 6.1%, with a mean depth of 2.1 km and a greatest depth of 6.7 km.
- **Smaller waters.**
  - Smythii–Marginis and Humboldtianum seas, about 0.75% each;
  - Orientale and Moscoviense lakes;
  - crater lakes, about 5.6% of the Moon in total.
- **Maria.**
  - Fully flooded: Imbrium, Serenitatis, Crisium, Nectaris, Nubium and Humorum.
  - Mostly flooded: Procellarum, Frigoris and Fecunditatis.
  - About one-third flooded: Tranquillitatis, in its southwest arm.
  - Still land: Vaporum and Insularum.
- **Islands.** The largest are Caucasus Island (68,000 km², summit 4.2 km) and Jura Island (58,000 km²), which encloses the flooded Sinus Iridum.

| Share | Sea level | Near-side water | Near-side maria joined to Imbrium |
|---|---|---|---|
| 25% | −1,802 m | 32% | 4 of 12 |
| 27% | −1,702 m | 35% | 9 of 12 |
| 28% | −1,654 m | 37% | 9 of 12 |
| 30% | −1,557 m | 40% | 9 of 12 |
| 35% | −1,328 m | 48% | 9 of 12 |

The atlas uses 4 px/deg for global geometry, and its levels match the 16 px/deg level table within 10 m. Poleward
of about 80°, the cylindrical grids disagree by hundreds of metres for small craters, so polar work moves to LOLA's
polar stereographic grids.

## Rivers and rain-fed lakes: a first estimate

`python -m geography.drainage` routes water over LOLA at 16 px/deg (1.9 km) with the atlas's sea level. It writes
[results/drainage.json](results/drainage.json) (schema `terluna.geography.drainage/1`) and the grid product
`products/drainage_28pct_16ppd.npz`, which is kept out of Git and takes about 20 seconds to regenerate.

- **Routing.** Each cell drains to the steepest of its eight neighbours. The floor of each closed depression spills
  through the lowest pass on its way to the sea. A priority flood over the graph of catchments and their shared
  passes finds those passes.
- **Water balance.** Runoff is precipitation less evaporation over land, where positive, from climate run
  A28_dim5_moon (model years 20–29, T21): the 28% seas with the lakes as water, PlaSim's cloud scheme corrected for
  lunar gravity, the author's chosen shield, 5% dimmer than the titania stack, and the GCM's Sun, tilt and sunlight
  split corrected (see climate/gcm). A depression fills and overflows when its inflow and the rain on its lake exceed the
  lake's open-water evaporation (the run's evaporation over sea at that latitude). Otherwise it keeps the smaller
  lake whose evaporation balances its inflow.

Results at 28%:

- **Lakes.** 180,000 closed depressions lie above sea level. 27,900 of them hold rain-fed lakes, 2,500 of those
  without an outlet. Together they cover 10.3% of the Moon on top of the 28% sea and hold 3.1 million km³, three
  tenths of the seas' volume.
- **Where.** Lakes cover 24% of the surface within 15° of the equator, 12% between 15° and 30°, 5% between 30°
  and 45° and 0.1% between 45° and 60°, where the corrected climate hardly rains (22, 13, 10 and 4% with the
  climate before its correction). The largest fill far-side
  basins: Hertzsprung (200,000 km², up to 4.7 km deep, its surface 5.1 km above sea level) and Korolev (141,000 km²,
  up to 6.7 km deep, 7.8 km above sea level).
- **Rivers.** 229,000 m³/s reaches the seas, about 1.1 times Earth's runoff per unit of land. The largest river
  enters the Orientale sea with 57,400 m³/s from a basin of 833,000 km².
- **Lakes and climate agree.** The first estimate took its runoff from run A (25% water, PlaSim's clouds as they
  were): 39,400 lakes over 12.0% of the Moon. Run A28 was set up with those lakes as water; its runoff gave
  12.9%, the chosen dimmer shield's climate 11.6%, and the same climate with the GCM's Sun, tilt and sunlight
  split corrected gives the 10.3% above. At the GCM's T21 resolution these lakes differ from A28's coastline by
  111 of 2,048 cells (68 water cells become land and 43 land cells water) and lower the water from 40.1% to
  38.3%. Whether to rebuild the GCM's coastline with them, which needs a cold start, is open (climate/gcm).

With a fixed water inventory, water held in lakes above sea level comes out of the seas. The sea level then falls
below −1,654 m unless the lakes' water is delivered on top, or groundwater and infiltration keep the basins
drier. [tests/test_drainage.py](tests/test_drainage.py) checks
mass conservation, fill-and-spill and the closed-lake balance on synthetic terrain.

The [terrain study](../research/studies/terrain_through_time/README.md) fills these lakes and closed hollows with the
land's eroded rock: lakes under 10 km² last a median 4,600 years at the central erosion and the largest 1.5 million,
and each lake's and hollow's fill times are in its grids for the atlas round.

## Groundwater: a first estimate

The crust under the seas and land is porous. GRAIL's gravity gives the highland
crust a bulk density of 2550 kg/m³ and an average porosity of 12% to depths of
at least a few kilometres (Wieczorek et al. 2013), about 4% at 20 km (as
summarised by Wiggins et al. 2022). `python -m geography.groundwater` fits an
exponential profile through both (13.5% at the surface, falling off over
16.4 km) and counts the pore water once the crust is saturated to a given depth
([results/groundwater.json](results/groundwater.json)):

| Crust saturated to | Pore water (global layer) | Seas + lakes + crust | Against the seas alone |
|---|---|---|---|
| none | 0 | 1.38×10¹⁹ kg (364 m) | 1.29 |
| 1 km | 131 m | 1.88×10¹⁹ kg | 1.76 |
| 4 km | 480 m | 3.20×10¹⁹ kg | 2.99 |
| 10 km | 1,013 m | 5.22×10¹⁹ kg | 4.88 |
| 20 km | 1,563 m | 7.31×10¹⁹ kg | 6.83 |

The seas are a 282 m global layer and the lakes 82 m. So the crust, not the
seas, may set how much water has to be delivered for the seas to stay at 28%:
saturating even its upper kilometre takes half as much again as the seas hold,
and the whole porous crust five to seven times as much. Seas laid on dry crust
drain into it until the pores beneath are full; how fast, and how deep the water
reaches, depend on a permeability unknown within orders of magnitude, so the
saturation depth is a question of time this estimate does not settle; the next
section follows it. It leaves out the maria's less porous basalt fill, closed
pores, water bound into new minerals and the water table's shape on land.

## Groundwater: how fast the crust takes the water

`python -m geography.groundwater_uptake` writes [results/groundwater_uptake.json](results/groundwater_uptake.json)
(schema `terluna.geography.groundwater-uptake/1`) in a few seconds; [groundwater_sources.json](groundwater_sources.json)
records the published values and how far each was read. It follows heat and water down one-dimensional columns under
the seas' floors, the rain-fed lakes' floors and the land in six latitude bands, from the moment the planned climate
and the water are in place.

- **The airless crust is cold near the surface and hot below.** Its shattered megaregolith, 2-3 km thick under the
  highlands and 1 km under the maria, conducts heat at about 0.2 W/m/K (Warren & Rasmussen 1987), so the Moon's
  13-17 mW/m² of internal heat warms it by 65-85 K a kilometre from about 255 K at the equator, below the month's
  heat wave, to some 440 K at the highland megaregolith's base. Rock below freezing reaches 240 m down at the equator,
  630 m at 45-60°, about 1.5 km under the polar highlands and 5 km under the polar seas' mare floors. Water entering
  it freezes and seals it, as in permafrost.
- **The planned climate thaws it from above.** Ground at 289-294 K and floors near 290 K warm the cold layer through
  in 530-560 years at the equator, 2,100-2,200 years at 30-45°, 5,600-7,900 at 45-60°, 18,000-21,000 at 60-75° and
  56,000-250,000 near the poles (central values).
- **Once thawed, the crust takes water fast.** The water runs down into the warm crust below under the seas' and
  lakes' head. Moon-wide the crust holds a 103 m global layer within 1,000 years (14-296 m across the cases), 514 m
  within 10,000 (198-840) and 1,075 m within 100,000 (528-1,264), against the seas' 282 m: 0.36, 1.8 and 3.8 times
  the seas. A crust with no cold layer would hold 237 m by 1,000 years, so the cold halves the first millennium's
  uptake and delays it by centuries. Over ten million years it nears the 1,560 m its pores hold to 20 km.
- **Holding the seas takes a top-up of thousands of cubic kilometres a year.** At the central values the crust takes
  1,500-4,100 km³ a year through the first thousand years and 1,700 km³ a year to 10,000, falling to 240 by 100,000
  and about 10 by a million; the fast case's first millennium takes 12,000 a year. The rivers bring the seas
  7,200 km³ a year. Without the top-up the seas would fall 17 m in a century and 370 m within a thousand years; a
  fall of 50 m splits the Near-side Sea, whose maria join about 50 m below its planned level.
- **On land the crust drinks the rain that soaks in.** In the wet tropics it takes about 60% of the recharge in the
  first thousand years and 20% to 10,000; between 15° and 30° half to all of it for 10,000 years; between 30° and 45°
  nearly all of it for 10,000 years and half to 100,000; on the dry lands all of it for a million. Springs and baseflow
  stay weak while it fills, the water table's next stage.
- **The floors' sediment slows it little.** Under kilometres of head a few metres of clay pass water faster than the
  deep crust takes it, and the deep crust's low permeability sets the pace; the tightest seal costs the slow case a
  quarter of its uptake by 10,000 years.

| Moon-wide water held in the crust, global layer | 100 years | 1,000 | 10,000 | 100,000 | 1 million |
|---|---|---|---|---|---|
| Slow | 3 m | 14 m | 198 m | 528 m | 991 m |
| Central | 5 m | 103 m | 514 m | 1,075 m | 1,347 m |
| Fast | 7 m | 296 m | 840 m | 1,264 m | 1,347 m |
| Central with a crust that never freezes the water | 105 m | 237 m | 555 m | 1,087 m | 1,347 m |
| Central with an airless megaregolith of 0.5 W/m/K | 4 m | 17 m | 294 m | 884 m | 1,342 m |
| The seas, for scale | 282 m | 282 m | 282 m | 282 m | 282 m |

| Top-up that holds the seas, km³ a year | Slow | Central | Fast |
|---|---|---|---|
| First 10 years | 3,100 | 4,100 | 5,700 |
| 10-100 years | 1,100 | 1,500 | 2,300 |
| 100-1,000 years | 440 | 4,100 | 12,200 |
| 1,000-10,000 years | 780 | 1,700 | 2,300 |
| 10,000-100,000 years | 140 | 240 | 180 |
| 100,000 years to 1 million | 20 | 11 | 4 |

**How it is set up.** Heat moves by implicit conduction: the megaregolith's conductivity rises from its airless
0.2 W/m/K to 0.3-0.8 once air fills its cracks and 1.0-2.0 once water does (*design guesses*), the bedrock below
conducts 1.4 W/m/K, and water raises the heat capacity. The surface holds the ground's temperature from the climate
run or the floors' water near the coldest sea surface, and the internal heat enters from below at 30 km. Water
infiltrates on Green and Ampt's terms through the crust's permeability, log k = -14 - 3.2 log z at a kilometre's depth
(Manning & Ingebritsen 1999) an order either way, capped at volcanic rock's 10⁻¹² m² near the surface (Gleeson et
al. 2011), at lunar gravity and the water's temperature. It is ponded under the seas and lakes at their mean depths
and fed on land by 10-30% of the rain (Earth's mean 17%; Moeck et al. 2020), at most the rain less evaporation and
at least 1% of the rain on dry land (Scanlon et al. 2006). It stops where the rock is below freezing, and the floors'
fine sediment builds up as a seal at 0.08-0.36 mm a year under the seas and 0.15-2.4 mm under the lakes (*design
guesses* from the sea chemistry and terrain studies' totals). Porosity follows the profile above to 20 km, where the
pores close. The slow and fast cases set every value toward slower or faster uptake together. The same formulation
on a warm Martian crust (Shadab et al. 2025) reaches a water table 2-5 km down in 50-200 years; for the same rock the
Moon's gravity makes it 2.3 times slower.

**What it leaves out.** The build's own warming of the ground before the water arrives, which would thaw the low
latitudes' cold layer sooner; the delivery's heat; water flowing sideways, from the seas under the land; the heat the
inflowing water carries and the circulation it would set up through the hot megaregolith, with warm springs; air
trapped ahead of the water; the maria's less porous basalt fill, for which a map of the maria's extent is needed;
clays from weathering closing pores; and the crust's response to the water's load.

## Groundwater: the water table on land

`python -m geography.groundwater_table` writes [results/groundwater_table.json](results/groundwater_table.json)
(schema `terluna.geography.groundwater-table/1`) in about a minute and a half; its grids go to
`research/runs/resources/groundwater/groundwater_table_2ppd.npz` (ignored by Git) for the atlas round. It solves the
groundwater's sideways flow over the atlas averaged to 15 km, about the resolution of Earth's global groundwater
models (de Graaf et al. 2015; Reinecke et al. 2019), for the mature Moon once the deep crust has filled and for the
young Moon at 1,000 and 10,000 years.

- **In the mature Moon groundwater comes out across a third of the land.** The water table stands at the ground across
  32% of the land (25-39% across the cases), and groundwater comes out there as springs, seeps, wetlands and baseflow:
  three quarters of the land within 15° of the equator, two fifths at 15-30°, a fifth at 30-45°, 4% at 45-60° and almost
  none nearer the poles. Elsewhere it lies deep, a median 364 m below the ground (170-635), deeper than 200 m under half
  the land and some 800 m under the dry high latitudes.
- **A quarter of the rivers' flow is groundwater, and half of the water that comes off the land.** Of the 1,870 km³
  a year that recharges the crust (1,070-2,490), 92% comes out at the surface, 6% seeps into the lakes and 2% into the
  seas: a quarter of the rivers' 7,240 km³ a year (15-34%), where Earth's rivers draw 59% of their flow from
  groundwater (Xie et al. 2024). Nearly half of the Moon's river flow is rain on its lakes, net of their evaporation,
  and the lakes even out what passes through them; of the water that comes off the land itself, groundwater is about
  half ([the springs' load](../research/studies/springs/README.md)).
- **The lakes gain groundwater, except at the dry poles.** On the 15 km grid 886 of 1,061 lakes, 96% of the lake area,
  take groundwater in, a median 5 cm a year. Of the 47 lakes losing more than 10 cm a year, 46 lie poleward of 75°,
  perched over deep water tables in the dry polar lands; they hold only if their beds seal, and the narrow polar cells
  overstate their rates.
- **The young Moon is drier, and unevenly.** At 1,000 years the crust still takes three fifths of the tropics' recharge
  and half of 15-30°'s: groundwater comes out across 23% of the land, 63% of the tropics, a third of 15-30° and 5% of
  30-45°. At 10,000 years the tropics are wetter, 73%, and the middle latitudes drier, a tenth of 15-30° and 2.5% of
  30-45°, since there the crust takes nearly all the recharge for ten thousand years. With the fastest-filling crust
  the first thousand years' water table lies a kilometre down at the median and the lakes lose 44 km³ a year to it.
- **The seas trade little.** They take 35 km³ a year of groundwater (21-61); near low coasts they feed the land's
  crust less than 1 km³ a year.
- **Crater hollows hold ponds where groundwater comes out.** Across three quarters of the mature tropics and two fifths
  of 15-30° the water table reaches the ground within each 15 km cell, so its hollows hold water; poleward of 45° they
  mostly stay dry. Under 2% of the land has its water table within 10 m of the ground without reaching it: at this
  resolution a cell either seeps somewhere or keeps its water table far below its mean ground.

| Mature Moon, central | Recharge | Land where it comes out | Median depth to water | Groundwater out at the ground |
|---|---|---|---|---|
| 0-15° | 264 mm a year | 76% | 0 m | 1,264 km³ a year |
| 15-30° | 77 mm | 43% | 60 m | 371 km³ |
| 30-45° | 19 mm | 18% | 295 m | 79 km³ |
| 45-60° | 0.9 mm | 4% | 556 m | 3 km³ |
| 60-75° | 0.2 mm | 0.5% | 808 m | 0.3 km³ |
| 75-90° | 0.4 mm | 0.2% | 808 m | 5 km³ |

| Land where groundwater comes out, and the median depth to water | 1,000 years | 10,000 years | Mature |
|---|---|---|---|
| Slow (permeability 10⁻¹⁵ m² at 1 km, 10% of the rain recharging) | 35%, 389 m | 32%, 419 m | 39%, 170 m |
| Central | 23%, 687 m | 20%, 812 m | 32%, 364 m |
| Fast (10⁻¹³ m², 30%) | 4%, 1,047 m | 16%, 1,010 m | 25%, 635 m |

The slow crust carries less water sideways and takes less down, so its land is the wettest.

**How it is set up.** The seas and the drainage product's lakes, cells by majority, hold their levels. Recharge on each
cell follows stage 1's rule from the design climate. The transmissivity is stage 1's permeability integrated below the
water table: 5×10⁻⁴ m²/s with the water table at the ground at the central values, 3×10⁻⁵ at 500 m and 1.5×10⁻⁶ at
2 km. A smoothed drain lets water out where the water table reaches the ground. Newton's method with pseudo-transient
continuation solves each state to a root-mean-square imbalance of 10⁻⁴ of the recharge, the mature state from a water
table at the ground and the young states from the mature. The young states take stage 1's share of each band's
recharge that the deep crust still takes.

**What it leaves out.** Relief within a 15 km cell, where the valleys seep and the ridges stay dry; the atlas round's
districts need finer regional solutions, which this solver did not reach on the 7.6 km grid. Evaporation from shallow
water tables and wetlands, which would lower the seepage in the dry belt. The lakes' beds, left open, so their leakage
is at its most. Flow through the deep crust; the water table's own changes through time, the young states being
steady at each time; and the maria's basalt.

### Shores by the fresh seas

15 October 2026, a *screen* of the water table under the shores, which the 15 km solver leaves out; the author asked
for the shorelines and the seas' freshness to be weighed more. The shore study at tens of metres is
[resources' question 16](../resources/questions.md).

- **The sea feeds the shore's water table.** Under a salt sea, fresh groundwater floats on seawater as a lens reaching
  about 40 times as far below sea level as its water table stands above it (the Ghyben–Herzberg relation, set by the
  two densities and so the same at any gravity), and plants that draw on it pull salt water up. Under the held fresh
  seas the crust is fresh to its base, the sea itself supplies the water table at the shore, and the spray carries no
  salt.
- **Lunar gravity lifts water higher in the soil and moves it more slowly through the ground.** In Gardner's
  exponential soil, the height above the water table from which a soil still lifts a flux q to the roots is
  ln(1 + εK/q)/(αε), with ε the ratio of lunar to Earth gravity (0.166) and the conductivity K and Gardner's α taken
  in Earth's terms. A sandy loam (K 10⁻⁵ m/s, α 2 m⁻¹) lifts 3 mm a day from 11.7 m on the Moon and 2.8 m on Earth,
  and a finer soil (α 1 m⁻¹) from 23 m and 5.7 m. The same rock conducts a sixth as well, so the sea's water moves
  sideways through it six times more slowly; stage 2 already integrates at lunar gravity.
- **Shore belts.** With 3–5 m of roots, plants draw on the water table from 15–25 m above the sea, against 6–10 m on
  Earth, so on steep shores their belts are two to three times as wide. On flat shores the ground's supply sets the
  width: a belt drawing E from the sea lowers the water table at its inland edge by EW²/2T, so at a metre a year and
  stage 2's central 5×10⁻⁴ m²/s, with 10–20 m to spare, it is 0.56–0.79 km wide, and 1.8–2.5 km at ten times the
  transmissivity, as in sandy ground. A freshwater shore on Earth, six times as transmissive with 4 m to spare, gives
  0.87 km.
- **Small islands.** The sea lies under the whole of a small island. Wooded throughout, a long island of half-width a
  lowers its middle's water table by Ea²/2T: 8 m for an island a kilometre across at the central transmissivity and
  32 m for one two kilometres across, a tenth of that in sandy ground. Low islands up to about a kilometre across in
  ordinary ground, or a few kilometres in sand, can stay wooded from shore to shore on groundwater alone.
- **The sea's load stays where the draw ends.** A metre a year at 0.14 g/kg leaves 140 g/m² a year. Most of it is
  calcite, magnesium minerals and opal, a hardpan that cements about 18 cm of soil at 30% porosity each millennium.
  Sodium and potassium are 2% of weathering's load ([sea chemistry](../research/studies/sea_chemistry/README.md)), so
  salt in the harmful sense builds slowly. On the Okavango Delta's islands, whose water is fresher than the held seas
  but carries more sodium beside its calcium, trees drawing the delta's water leave wooded rims round bare, salty
  centres (McCarthy and Ellery 1994; *recalled*). Rain, the monthly tide (3.7 m typical on the Near-side Sea) or
  flushing by people carries the load away.
- **Jura Island.** Measured at LOLA's 16 px/deg on research/jura-island, only 2% of the island lies within 25 m of the
  sea and its first 1–2 km inland stand a median 112 m up, so its belts are a rim of 1,100–2,600 km², 2–4.5% of it.

## The monthly tide

Earth raises a tide on the Moon 12.9 m high at the sub-Earth point. Its permanent
part already shapes the atlas's equipotential. The orbit's eccentricity and the
evection change the Earth's distance through each month, and the optical
librations swing the sub-Earth point by about 7° in longitude and latitude, so
the bulge swells, shrinks and rocks. The Sun adds 7 cm. `python -m geography.tides`
(about five minutes) writes [results/tides.json](results/tides.json) (schema
`terluna.geography.tides/1`) and the grid product `products/tides_28pct_4ppd.npz`,
kept out of Git.

- **Ephemeris.** [lunar_ephemeris.py](lunar_ephemeris.py) takes the truncated
  ELP-2000/82 lunar theory and the optical librations as tabulated by Meeus
  (*Astronomical Algorithms*, 1998, chapters 47 and 53). Against JPL Horizons
  (DE441) over 2026–2045, the sub-Earth and sub-solar points agree within
  0.05°, the physical librations it leaves out, and the Earth's distance within
  11 km.
- **Height.** The degree-2 and degree-3 potential of the Earth and the Sun, every
  three hours from 2026 to 2045, through the 18.6-year nodal cycle, on the atlas's
  quarter-degree nodes. The solid Moon's own tide enters through the Love numbers
  k2 = 0.02405 (GRAIL) and h2 = 0.0371 (LOLA). Each sea keeps its volume. The
  water's self-attraction, solved for all seas together through spherical
  harmonics to degree 60, adds 9% in the largest seas.

Ranges below are the highest minus the lowest level within an anomalistic month
(27.55 days), the median over the record's 251 months, weighted by area over
each sea:

| Sea | Typical monthly range | 95th percentile of the sea's area | Largest in 19 years |
|---|---:|---:|---:|
| Nearside | 3.7 m | 6.0 m | 9.8 m |
| South Pole–Aitken | 2.1 m | 3.5 m | 6.7 m |
| Humboldtianum | 0.85 m | 1.9 m | 3.9 m |
| Smythii–Marginis | 0.55 m | 1.4 m | 2.8 m |
| Moscoviense | 0.40 m | 0.61 m | 1.1 m |
| Orientale | 0.37 m | 0.81 m | 1.5 m |

The nearside sea's typical monthly range is 3.3–3.5 m at Procellarum, Imbrium
and Frigoris, 4.1 m in Serenitatis, 5.0 m in Tranquillitatis and 5.5–6.0 m in
Fecunditatis, Nubium and Humorum. Mare Ingenii's is 3.1 m. Mare Smythii's is
0.23 m and the eastern Smythii headland of the wave studies 0.78 m (1.0 m at
most): there the monthly tide moves the waterline 40–80 m on 1:50 and 1:100
beaches.

The forcing is a monthly pair: the anomalistic month M′ (27.55 days, from the
eccentricity and the libration in longitude) and the draconic month F (27.21 days,
the libration in latitude), each up to 2.2 m. The evection 2D − M′ (31.8 days)
and 2D (14.77 days) reach 0.44 m. The fortnightly lines M′ + F, 2M′ and 2F
(13.6–13.8 days) reach 0.13–0.31 m, and 2D + M′ (9.61 days) 0.09 m. Twelve lines
hold 95% of the forcing's variance.

[tide_dynamics.py](tide_dynamics.py) solves the linear shallow-water equations
for the eight leading lines on each main sea's quarter-degree cells, with the
Moon's rotation and bed friction for 1 cm/s flows. It covers the 94–99% of each
sea joined through shared cell faces.

- **Smythii–Marginis** follows its equilibrium within 1.4%. Its seiches take
  40, 21 and 17 hours.
- **The South Pole–Aitken sea** follows within 6% in RMS. A sub-basin by Lacus
  Tenebrarum, near the pole, joins it through narrow passages and, at
  quarter-degree cells, exchanges water over 15 days.
- **The nearside sea's monthly lines** follow within 7% in RMS.
- **Mare Fecunditatis** joins the nearside sea through a strait near 40° E, 7° S,
  60–180 m deep and a few cells wide. It exchanges water with the main sea every
  9.5 days: 13.0 at half-degree cells and 9.1 at whole degrees, so the strait's
  geometry sets the period. Inside Fecunditatis the fortnightly lines rise to
  about twice their equilibrium height. The 9.61-day line, close to the exchange
  period, resonates: 39 times its equilibrium height with bed friction for
  1 cm/s flows, 5.4 times with ten times that friction and 0.6 times with a
  hundred times. A strait moving a metre of Fecunditatis's level flows at
  about 1 m/s, the hundredfold case. The strait's true geometry and nonlinear
  friction set Fecunditatis's tide.

The self-attraction holds the solid Moon rigid under the water's load and uses
the 1,025 kg/m³ water density of the wave studies. Fresh water would lower the
9% self-attraction by 2.5% of itself. Degree 3 also uses a rigid Moon and adds
at most 6 cm.

## Limits

- **Datums:** LOLA (mean-Earth frame) and GL0420A (principal-axis frame) differ
  by about 0.02°, below the geoid's resolution. A degree-100 geoid differs from
  degree 200 by 5 m RMS.
- **Resolution:** the depression tree uses the 4 pixel/degree grid, where small
  craters are unresolved; level curves use 16 pixels/degree.
- **Not modelled:** crustal loading by the water (subsidence of order 10–30% of
  the depth under large seas), groundwater beyond the first estimate above, polar
  ice and the water held in the air and soil (the design case's air holds about
  180 kg/m², 0.06% of the seas).

## Next work

- How fast, and how deep, the crust takes up water, which sets the inventory to deliver (see "Groundwater").
- The straits of Mare Fecunditatis and of the Tenebrarum sub-basin at 16 px/deg, with nonlinear bed friction, which
  set those basins' fortnightly and 9.6-day tides; the tide's water levels in the wave and run-up calculations.
- Groundwater storage in the porous crust (GRAIL: about 12%) and crustal loading, which set how much water
  to deliver for 28% cover and where the shore falls.
- LOLA polar stereographic grids for the polar cold traps.
- Land–sea masks for chosen inventories are available:
  `python -m geography.water_inventory --masks 0.25 0.35 --resolution 1.0` writes
  `products/land_sea_{25,35}pct_1deg.npz` (schema `terluna.geography.land-sea/1`:
  water fraction, mean ground height and mean water depth per cell; not
  committed).
- Polar cold-trapping of water.
- Crustal loading.

Couple candidate geographies to [climate](../climate/) and the water and element
accounts in [biosphere](../biosphere/). Keep several candidate geographies, so
the core can describe a whole world rather than one landscape. Older atlas,
mist-highland and inland-sea descriptions stay hypotheses until their water and
climate budgets are established.
