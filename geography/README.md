# Geography and hydrology

**Current condition:** global LOLA topography (4 and 16 pixels per degree) and
the GRAIL GL0420A gravity field are fetched and hash-checked
([inputs.json](inputs.json), [fetch_inputs.py](fetch_inputs.py)). Heights are
referred to the geoid, and hydrostatic water storage is computed. There is no
drainage, runoff, erosion or climate-coupled shoreline model yet.

| File | Holds |
|---|---|
| [topography.py](topography.py) | LOLA radius grids; the geoid from degrees 2–200 of GL0420A plus Earth's static tide (×(1 + k2)) and rotation; heights above that equipotential; exact spherical cell areas |
| [basins.py](basins.py) | Level curves (area and volume below a common water level), separate seas, and the depression tree: every closed depression's spill level and capacity, from a union-find merge by height |
| [water_inventory.py](water_inventory.py) | Writes [results/](results/): the level table and the major basin joins as water rises, with schema, hashes and evidence statement |

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
