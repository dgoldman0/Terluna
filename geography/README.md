# Geography and hydrology

**Current condition:** global LOLA topography (4 and 16 pixels per degree) and
the GRAIL GL0420A gravity field are fetched and hash-checked
([inputs.json](inputs.json), [fetch_inputs.py](fetch_inputs.py)). Heights are
referred to the geoid, and hydrostatic water storage is computed. The IAU lunar
nomenclature is fetched and hash-checked the same way and names the atlas's seas,
islands and targets. Drainage, runoff, erosion and climate-coupled shorelines are
the next layer.

| File | Holds |
|---|---|
| [topography.py](topography.py) | LOLA radius grids; the geoid from degrees 2–200 of GL0420A plus Earth's static tide (×(1 + k2)) and rotation; heights above that equipotential; exact spherical cell areas |
| [basins.py](basins.py) | Level curves (area and volume below a common water level), separate seas, and the depression tree: every closed depression's spill level and capacity, from a union-find merge by height |
| [water_inventory.py](water_inventory.py) | Writes [results/](results/): the level table and the major basin joins as water rises, with schema, hashes and evidence statement |
| [atlas.py](atlas.py) | The atlas at the scenario water share: seas and lakes with depths and IAU water names, islands with summits, mare flooding and a comparison of candidate shares ([results/atlas.json](results/atlas.json); grid product in `products/`) |
| [nomenclature.py](nomenclature.py) | IAU lunar feature names from the USGS gazetteer archive (read by a small dBase reader built on the standard library) |

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

## Limits

- **Datums:** LOLA (mean-Earth frame) and GL0420A (principal-axis frame) differ
  by about 0.02°, below the geoid's resolution. A degree-100 geoid differs from
  degree 200 by 5 m RMS.
- **Resolution:** the depression tree uses the 4 pixel/degree grid, where small
  craters are unresolved; level curves use 16 pixels/degree.
- **Not modelled:** crustal loading by the water (subsidence of order 10–30% of
  the depth under large seas), groundwater, polar ice and the water held in the
  air and soil.

## Next work

- Fill-and-spill driven by rainfall over catchments. Where water actually stands
  depends on precipitation minus evaporation, which needs the climate model.
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
