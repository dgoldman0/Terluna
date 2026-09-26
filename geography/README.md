# Geography and hydrology

**Current condition:** global LOLA topography (4 and 16 pixels per degree) and
the GRAIL GL0420A gravity field are fetched and hash-checked
([inputs.json](inputs.json), [fetch_inputs.py](fetch_inputs.py)). Heights are
referred to the geoid, and hydrostatic water storage is computed. The IAU lunar
nomenclature is fetched and hash-checked the same way and names the atlas's seas,
islands and targets. A first estimate of rivers and rain-fed lakes routes the
climate run's runoff over the terrain. Erosion, groundwater and climate-coupled
shorelines are the next layer.

| File | Holds |
|---|---|
| [topography.py](topography.py) | LOLA radius grids; the geoid from degrees 2–200 of GL0420A plus Earth's static tide (×(1 + k2)) and rotation; heights above that equipotential; exact spherical cell areas |
| [basins.py](basins.py) | Level curves (area and volume below a common water level), separate seas, and the depression tree: every closed depression's spill level and capacity, from a union-find merge by height |
| [water_inventory.py](water_inventory.py) | Writes [results/](results/): the level table and the major basin joins as water rises, with schema, hashes and evidence statement |
| [atlas.py](atlas.py) | The atlas at the scenario water share: seas and lakes with depths and IAU water names, islands with summits, mare flooding and a comparison of candidate shares ([results/atlas.json](results/atlas.json); grid product in `products/`) |
| [nomenclature.py](nomenclature.py) | IAU lunar feature names from the USGS gazetteer archive (read by a small dBase reader built on the standard library) |
| [drainage.py](drainage.py) | Rivers and rain-fed lakes above sea level: runoff from the climate run routed over the 16 px/deg terrain, with fill-and-spill lakes set by each depression's water balance ([results/drainage.json](results/drainage.json); grid product in `products/`). `--climatology` takes another run's climate product, `--out` writes elsewhere |
| [groundwater.py](groundwater.py) | A first estimate of the water the porous crust takes up, from GRAIL's porosity, beside the seas and lakes ([results/groundwater.json](results/groundwater.json)) |

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
- **Water balance.** Runoff is precipitation less evaporation over land, where positive, from climate run A28
  (model years 20–29, T21): the 28% seas with the lakes as water, and PlaSim's cloud scheme corrected for lunar
  gravity (see climate/gcm). A depression fills and overflows when its inflow and the rain on its lake exceed the
  lake's open-water evaporation (the run's evaporation over sea at that latitude). Otherwise it keeps the smaller
  lake whose evaporation balances its inflow.

Results at 28%:

- **Lakes.** 180,000 closed depressions lie above sea level. 47,300 of them hold rain-fed lakes, 15,100 of those
  without an outlet. Together they cover 12.9% of the Moon on top of the 28% sea and hold 3.8 million km³, about
  a third of the seas' volume.
- **Where.** Lakes cover 17–19% of the equatorial belt and 0.1% of the dry polar regions. The largest fill far-side
  basins: Hertzsprung (200,000 km², up to 4.7 km deep, its surface 5.1 km above sea level) and Korolev (141,000 km²,
  up to 6.7 km deep, 7.8 km above sea level).
- **Rivers.** 370,000 m³/s reaches the seas, about 1.7 times Earth's runoff per unit of land. The largest river
  enters the Orientale sea with 67,100 m³/s from a basin of 833,000 km².
- **Lakes and climate agree.** The first estimate took its runoff from run A (25% water, PlaSim's clouds as they
  were): 39,400 lakes over 12.0% of the Moon. Run A28 was set up with those lakes as water; its runoff, 22%
  larger, gives the lakes above. At the GCM's T21 resolution they would flip 84 of 2,048 cells (51 to water, 33
  to land) and raise the water from 40.1% to 40.9%, worth about 0.1 K, so no further round is needed.

With a fixed water inventory, water held in lakes above sea level comes out of the seas. The sea level then falls
below −1,654 m unless the lakes' water is delivered on top, or groundwater and infiltration keep the basins
drier. [tests/test_drainage.py](tests/test_drainage.py) checks
mass conservation, fill-and-spill and the closed-lake balance on synthetic terrain.

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
| none | 0 | 1.45×10¹⁹ kg (382 m) | 1.35 |
| 1 km | 131 m | 1.94×10¹⁹ kg | 1.82 |
| 4 km | 480 m | 3.27×10¹⁹ kg | 3.05 |
| 10 km | 1,013 m | 5.29×10¹⁹ kg | 4.94 |
| 20 km | 1,563 m | 7.38×10¹⁹ kg | 6.89 |

The seas are a 282 m global layer and the lakes 99 m. So the crust, not the
seas, may set how much water has to be delivered for the seas to stay at 28%:
saturating even its upper kilometre takes half as much again as the seas hold,
and the whole porous crust five to seven times as much. Seas laid on dry crust
drain into it until the pores beneath are full; how fast, and how deep the water
reaches, depend on a permeability unknown within orders of magnitude, so the
saturation depth is a question of time this estimate does not settle. It leaves
out the maria's less porous basalt fill, closed pores, water bound into new
minerals and the water table's shape on land.

## Limits

- **Datums:** LOLA (mean-Earth frame) and GL0420A (principal-axis frame) differ
  by about 0.02°, below the geoid's resolution. A degree-100 geoid differs from
  degree 200 by 5 m RMS.
- **Resolution:** the depression tree uses the 4 pixel/degree grid, where small
  craters are unresolved; level curves use 16 pixels/degree.
- **Not modelled:** crustal loading by the water (subsidence of order 10–30% of
  the depth under large seas), groundwater beyond the first estimate above, polar
  ice and the water held in the air and soil (the design case's air holds about
  260 kg/m², 0.1% of the seas).

## Next work

- How fast, and how deep, the crust takes up water, which sets the inventory to deliver (see "Groundwater").
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
