# The springs' load

14 October 2026. Groundwater's third stage, at the author's word: how long the water stays underground, how much rock it
dissolves on the way, how hard the springs are and where they build travertine, and what the groundwater carries of the
weathering and carbon the carbon cycle counts. Stage 1 ([how fast the crust takes the
water](../../../geography/README.md#groundwater-how-fast-the-crust-takes-the-water)) and stage 2 ([the water table on
land](../../../geography/README.md#groundwater-the-water-table-on-land)) are geography's; this study couples them with
the carbon cycle's chemistry, so it sits among the studies
([register](../../decisions.md#the-terrain-through-time-13-october-2026)).

Two runners write two products:

- `python -m research.studies.springs.run` re-solves stage 2's nine states (the mature Moon and the young Moon at 1,000
  and 10,000 years, for stage 1's slow, central and fast crust), follows each land cell's water from rain to outflow
  and mixes it along the flow, then works out its chemistry for the central case and three variants. It writes
  [results/springs.json](results/springs.json) (schema `terluna.research.springs/1`) in about eighteen minutes with
  2.3 GB; its grids, each spring's median age, hardness, pH and saturation in air, go to
  `research/runs/resources/groundwater/springs_2ppd.npz` (ignored by Git) for the atlas round.
- `python -m research.studies.springs.carbon` reads that product with stage 1 and the carbon cycle and writes
  [results/springs_carbon.json](results/springs_carbon.json) (schema `terluna.research.springs-carbon/1`) in seconds:
  the carbon the filling crust takes down, and the weathering and calcite against the carbon cycle's figures. It is
  kept apart so that the carbon cycle's rebuilds reach it in seconds.

[sources.json](sources.json) records each published value and how far it was read.

## What it finds

- **The springs' water is young.** In the mature Moon it reaches the springs a median 76 years after it fell as rain
  (76–101 across the cases): 58% of it within a century, 92% within 1,000 years and 2% older than 10,000. The tropics'
  springs run 66 years, 15–30°'s 107, 30–45°'s 200 and the few at 45–60° 7,600. The crust's permeability falls with
  depth, so most of the water crosses its top few hundred metres while the deep crust's old water reaches the springs a
  little at a time. In the young Moon the springs' water is a median 177 years old at 1,000 years and 85 at 10,000.
- **The springs carry about the soil water's load, because the rock below the soil runs out of CO2.** They carry 1.6
  meq/kg of alkalinity in the mature Moon, 1.25 times the soil water's 1.26, and 1.4–1.6 in the young Moon, 0.87–1.02
  times its soil water's. Below the soil the water has only the CO2 it took up there, 0.2–0.9 mmol/kg: the rock turns it
  into alkalinity and, past saturation, lays the calcium in the rock as calcite, about half of it in the mature Moon and
  nearly all of it in the young Moon, whose freshly wetted rock reacts 0.8 meq/kg a year at 1,000 years against 0.012 at
  a million.
- **They are moderately hard, and alkaline magnesium waters in the young Moon.** The springs hold 65–80 mg/L as CaCO3
  across the wet bands. The mature Moon's are calcium and magnesium bicarbonate waters near pH 7.1 with about fifteen
  times the air's CO2. The young Moon's are magnesium bicarbonate waters at pH 9.4–9.5 with 1–3 mg/L of calcium and less
  CO2 than the air, as the groundwater of Norway's Egersund anorthosites runs to pH 9.9 with calcium down to 1.5 mg/L
  (Banks & Frengstad 2006); where they come out they take CO2 from the air.
- **The springs build almost no travertine.** In air the mature Moon's springs come to calcite saturation (saturation
  index 0.02) and the young Moon's stay below it (−0.7 to −1.1); Earth's streams build tufa from an index of 0.75, 5.6
  times saturation. Twice saturation, where designed builders are taken to work, is reached by under 1% of the springs'
  flow (6% with the slow rock). With the rock above the water table as rich in CO2 as the soil, tufa forms only at the
  young Moon's dry-belt springs, a sixth of 30–45°'s spring flow at 1,000 years. The calcium goes into the rock instead:
  72 Mt of calcite a year in the mature Moon (34–181 across the variants) and 108–124 Mt in the young Moon, cementing
  the aquifers. The Moon's travertine is the young soils' rivers' ([terrain
  study](../terrain_through_time/README.md#travertine)), and lasting barriers in the mature Moon rest on builders.
- **The carbon cycle counts more weathering than the land exports.** The weathering that reaches the waters is 0.065 Gt
  C a year in the mature Moon, 59% of the carbon cycle's catchment total (56–61% across the cases), and 20–49% of it in
  the young Moon, least where the crust fills fastest. The carbon cycle applies the catchment law to the whole river
  flow, and nearly half of the rivers' 7,240 km³ a year is rain on the lakes, net of their evaporation, which touches no
  rock; the groundwater also keeps part of its calcium in the rock. Groundwater carries 55% of what the land exports
  (33–71% across the cases), and 8–54% in the young Moon, as the filling crust takes more or less of the recharge. The
  rock's calcite, 0.009 Gt C a year in the mature Moon, is 3.3 times the carbon cycle's soil carbonate.
- **The filling crust takes carbon down.** The water the crust draws from the land, the lakes and the seas carries 0.08
  Gt C a year down with it over 100–1,000 years and 0.05 over 1,000–10,000, 17% and 26% of the carbon cycle's import
  then, falling to 0.007 by 100,000 years; most of it is the land's water with the soil's CO2. Earth's flux of dissolved
  carbon into groundwater is about 0.2 Gt C a year (Kessler & Harvey 2001). The floors also take 4–14% of the alkalinity
  that weathering brings the waters over those spans.
- **Warm springs are rare.** 1.8% of the mature Moon's spring water has been below 0.93 km, where the wetted crust is 10
  K warmer than the ground (0.6–4.9% across the cases), and 0.16% below 2.8 km, where it is 30 K warmer. The few springs
  poleward of 45° carry the most deep water, a sixth to a third of theirs.

| Mature Moon, central | Spring flow | Median age | Soil CO2 | Alkalinity | Calcium | Magnesium | Hardness | pH | CO2 | Saturation in air |
|---|---|---|---|---|---|---|---|---|---|---|
| 0-15° | 1,263 km³ a year | 66 years | 2.6% | 1.61 meq/kg | 0.41 mmol/kg | 0.37 mmol/kg | 78 mg/L | 7.07 | 8,100 μatm | 0.06 |
| 15-30° | 371 km³ | 107 years | 2.0% | 1.54 | 0.33 | 0.41 | 74 mg/L | 7.26 | 4,900 μatm | −0.06 |
| 30-45° | 79 km³ | 198 years | 0.54% | 1.34 | 0.24 | 0.41 | 65 mg/L | 7.67 | 1,700 μatm | −0.34 |
| 45-60° | 3 km³ | 7,600 years | air's | 1.04 | 0.13 | 0.37 | 50 mg/L | 9.16 | 37 μatm | −0.79 |

The polar band's 4 km³ a year of spring flow is lake water passing through the banks of the polar lakes, which the
narrow polar cells exaggerate (stage 2).

| Springs, Moon-wide, central crust | Median age | Alkalinity | Calcium | Magnesium | Hardness | pH | Saturation in air | Calcite laid in the rock |
|---|---|---|---|---|---|---|---|---|
| Mature, central rock | 76 years | 1.58 meq/kg | 0.38 mmol/kg | 0.38 mmol/kg | 77 mg/L | 7.12 | 0.02 | 72 Mt a year |
| Mature, slow rock | 76 years | 1.52 | 0.43 | 0.30 | 74 mg/L | 6.80 | 0.04 | 34 Mt |
| Mature, fast rock | 76 years | 1.37 | 0.11 | 0.54 | 65 mg/L | 9.20 | −0.64 | 181 Mt |
| Mature, open rock | 76 years | 1.67 | 0.39 | 0.42 | 81 mg/L | 7.15 | 0.06 | 93 Mt |
| 10,000 years, central rock | 85 years | 1.44 | 0.08 | 0.60 | 68 mg/L | 9.36 | −0.74 | 124 Mt |
| 1,000 years, central rock | 177 years | 1.63 | 0.03 | 0.74 | 77 mg/L | 9.50 | −1.06 | 108 Mt |

| The filling crust's carbon, central | 0-10 years | 10-100 | 100-1,000 | 1,000-10,000 | 10,000-100,000 |
|---|---|---|---|---|---|
| Water it takes | 4,110 km³ a year | 1,510 | 4,140 | 1,730 | 236 |
| Carbon it takes down | 0.059 Gt C a year | 0.030 | 0.082 | 0.049 | 0.007 |
| Against the carbon cycle's import | 5% | 2% | 17% | 26% | 4% |

**Against Earth's groundwaters.** Basaltic aquifers worldwide carry mostly 1–10 meq/L of alkalinity at pH 6–8 (Romano
et al. 2024), and the Columbia Plateau's 2.8–3.3 meq/L at pH 7.6–7.9 after about 3,800 years (Steinkampf & Hearn
1996). The Moon's springs lie at the low end because its crust has no volcanic CO2 and its water is younger. Like
Derry et al.'s (2026) basaltic systems, where about 70% of the calcium leaves the water as carbonate in the deeper
vadose zone, groundwater and riparian zones, the Moon's groundwater lays much of its calcium underground. Beyond about
ten years of residence weathering comes under thermodynamic control (Maher 2010), and the Moon's springs, decades to
centuries old, sit at the limit the CO2 they carry sets.

**What it means for the other studies.** The terrain study takes springs at 1.5–3 times the rivers' load as a design
guess; the springs carry 0.9–1.3 times the soil water's. The carbon cycle's weathering sink would fall to about 0.6 of
its present value in the mature Moon and about 0.4 in the young Moon, and the filling crust adds a sink of up to a
quarter of its import over the first 10,000 years; carrying both into it, and the floors' share of the seas' alkalinity
into the sea chemistry, waits for the author.

## How it is set up

**Transit.** Stage 2's heads, recharge and seepage are re-solved with its own code, and the flows between cells come
from its transmissivities. The soil, a 1.5 m root zone of loam (*design guess* for its depth; Clapp & Hornberger's
loam), and the rock between the soil and the water table drain at unit gradient, each holding the water at which its
conductivity equals the recharge (Campbell 1974), so the water crosses them in plug flow. Below the water table each
cell's spread of ages follows a generalized Vogel model: with uniform recharge R in Dupuit flow, the share of the
outflow older than t is the share of the flow passing below the depth z that the streamline of age t reaches, with
t(z) = (1/R) ∫ φ / G dz from the water table, whatever the length of the flow path. With constant properties it gives
Vogel's exponential ages with mean φH/R (Haitjema 1995); with stage 1's permeability falling as depth^-3.2 below
about 240 m and its porosity falling over 16 km, it gives the power-law tail Cardenas & Jiang (2010) find for basins
whose permeability falls with depth. Each land cell's water, its recharge with any water the seas or lakes lose to it,
is mixed with the flows into it, and the mixture leaves as seepage and as flow to its neighbours, the seas and the
lakes; the mixing is solved as one sparse system over the steady flow. In the young states no water is older than the
Moon's water at that age.

**Chemistry.** The carbonate system at each cell's ground temperature, with the carbon cycle's constants (Plummer &
Busenberg 1982) and Davies activities, open to a gas or sealed. The soil water carries the carbon cycle's catchment load
(Dessert et al. 2003 at the land's air, 1.26 meq/kg, spread by temperature as the terrain study's routing spreads it) at
the young soils' factor of its central corner, with the highland soil's cations (65% calcium, 33% magnesium and 2%
sodium and potassium by equivalent), open to the soil's CO2 from Brook et al. (1983) on the climate run's evaporation:
2.6% in the tropics, 0.5% at 30–45° and the air's 480 μatm in the dry high latitudes. Below the soil the water is sealed
from gas. It keeps dissolving the rock at the wetted rock's rate, which falls as White & Brantley's (2003) field rates
do with the time since the water came: 0.52 meq/kg a year at 1,000 years, 0.13 at 10,000 and 0.008 at a million in the
central case at 15 C (0.80, 0.20 and 0.012 at the tropics' ground), on 100 m²/kg of reactive surface (Maher &
Chamberlain's default) over the 18 kg of rock a kilogram of pore water wets. It consumes the CO2 it carries, lays
calcite in the rock wherever it passes saturation, and stops at pH 9.5 (*design guess*, within the 9.4–9.9 that the
Columbia Plateau's basalt and Norway's anorthosite groundwaters reach). The variants: the slow rock reacts at the
Columbia Plateau's rate carried back to the Moon's younger rock; the fast rock has ten times the central surface; the
open rock keeps the rock above the water table open to CO2 at the soil's level, as Earth's deep unsaturated zones can be
(2% in a capillary fringe at 29 m, nearly 4% in karst boreholes to 60 m), with its alkalinity held below the Columbia
Plateau's highest, 6.7 meq/L.

**Springs and rivers.** The springs' water, mixed before it meets the air, is brought to the air's CO2 and judged
against the terrain study's thresholds: a saturation index of 0.75 for Earth's streams' tufa (Merz-Preiß & Riding
1999), three times saturation with pH above 8 for Plitvice's barriers, twice saturation for designed builders. The
rivers carry the quick flow off the land at the soil water's load, the groundwater that comes back up at its own, and
the rain on the lakes, net of their evaporation, at none.

## What it leaves out

The chemistry of the secondary clays and zeolites, which set where basaltic groundwaters stop (Derry et al. 2026); a
single pH stands in for them. Magnesium carbonates and silicates, which would cap the magnesium the dry bands' springs
gather. Evaporation from springs and shallow water tables, which concentrates their load in the dry belt. The maria's
basalt, which awaits a map of the maria's extent. Organic acids and the engineered biology's own effect on the soil's
CO2 and on weathering; deep CO2 sources, of which the Moon has none. Relief within a 15 km cell: valley springs fed by
short paths are younger and softer than the cell's mixture. The lakes' and seas' own chemistry, which is the carbon
cycle's.

## Files

| File | What it does |
|---|---|
| [run.py](run.py) | Transit times, the water's path and chemistry, the mixing to the springs, lakes and seas, and the summaries |
| [carbon.py](carbon.py) | The filling crust's carbon and the weathering and calcite against the carbon cycle |
| [test_run.py](test_run.py) | The chemistry against the carbon cycle's solver, the transit times against Vogel's exponential case and Campbell's water content, the mixing's conservation, the products' binding and agreement with stage 2, and an opt-in regeneration (`TERLUNA_SPRINGS_REGENERATE=1`, about twenty minutes) |
| [sources.json](sources.json) | Published values and how far each was read |
