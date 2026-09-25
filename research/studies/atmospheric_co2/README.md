# Atmospheric CO2 on the Open Moon

25 September 2026. What CO2 level the Open Moon's plants need, how fast the land
would take CO2 out of the air, and what a settled Moon has to return to keep it.

This is a literature synthesis with small calculations, not a carbon-cycle
model. `python -m research.studies.atmospheric_co2.run` writes
[results/atmospheric_co2.json](results/atmospheric_co2.json) (schema `terluna.research.atmospheric-co2/1`,
with its evidence statement and reading rule). [sources.json](sources.json)
records each source and how far it was checked. The atmosphere's response to CO2
comes from the atmosphere domain's results, cited below.

## What plants need

The design air holds O2 at Earth's partial pressure (21.2 kPa) in 1.2 atm, so
plants respond to the CO2 partial pressure in pascals, and Earth experiments
carry over in Pa. At 121.6 kPa, 1 Pa is 8.2 ppm; the design's 400 ppm is 48.6 Pa
(Earth today is about 41 Pa).

| Level | Partial pressure | ppm at 1.2 atm | Basis |
|---|---|---|---|
| Hard floor | about 20 Pa | about 165 | Glacial Earth: C3 crops roughly halved, C4 favoured |
| Diverse ecosystems | 28–30 Pa | about 240 | Pre-industrial Earth with a margin for heat and pressure |
| Productive C3 farming | about 35 Pa | about 290 | At the highest farmed ground |
| Comfortable | 43–60 Pa | 350–490 | Today's Earth rate after the pressure penalty, up to where gains flatten |

- **No sharp threshold.** C3 photosynthesis falls smoothly with CO2. The leaf
  model here (Farquhar–von Caemmerer–Berry with Bernacchi et al. 2001 kinetics,
  light-saturated, internal CO2 at 0.7 of ambient) gives 0.37, 0.67, 1.18 and
  1.30 times today's rate at 18, 28, 49 and 60 Pa at 25 °C. A meta-analysis of
  halved CO2 finds photosynthesis −38% and biomass −47% (Temme et al. 2013).
- **Warmth deepens every loss.** At 35 °C the model's 18 Pa rate falls to 0.19
  of today's, and experiments find the same (Cowling & Sage 1998). That matters
  on a Moon averaging about 26 °C.
- **The higher total pressure costs a little.** Gas diffuses more slowly in
  denser air. At a fixed stomatal conductance, 1.2 atm lowers the rate by 5–7%,
  and 43.5 Pa restores Earth's rate at 41 Pa.
- **Height matters.** The partial pressure falls about 9% per 5 km of height,
  and the lakes and farms of the high far side sit several kilometres up.
- **Humans set no floor.** The occupational limits (5,000 ppm) are ten times
  above the plant range.

The design's 49 Pa is comfortable; nothing here asks for a change.

## Climate and the upper air barely notice

Across this range CO2 is a weak climate lever. Each doubling changes the
Moon's clear-sky balance by about 4 W/m², as on Earth: −2.0 W/m² at 280 ppm and
−5.4 W/m² at 150 ppm against 400 ppm
(`atmosphere/radiative_convective/results/inverse_climate.csv`, the
`carbon_dioxide` sweep). With the GCM's 1.0–1.6 W/m² per K
(`climate/gcm/compare.py --settle`), that is 1.3–2 K and 3.5–5.5 K of cooling.
Less CO2 warms the upper air, its main coolant: the exobase behind the titania
film rises from 192–193 K to 197–198 K at 280 ppm and 206–208 K at 150 ppm, far
under the 250 K target ([atmosphere/middle_atmosphere](../../../atmosphere/middle_atmosphere/)).

## Build-up: imported

An Earth-like land biosphere (15–16 kg C per m² in plants and the top metre of
soil) on the Moon's 72% land holds 1.5–1.6×10¹⁵ kg CO2, 87–92% of the
atmosphere's whole 400-ppm inventory (1.73×10¹⁵ kg, 35.6 Gt per Pa). Regrowing
forest takes up its share in about 50 years; soils take centuries. The CO2 for
this is imported with the biosphere.

## The settled Moon

Once forests and soils stop growing, photosynthesis is balanced by respiration
and decay, apart from a small leak into peat and sediments. The geological loop
does not close by itself:

- **Weathering is a one-way sink.** Silicate rock takes CO2 into bicarbonate;
  carbonate forming in the seas returns half. With the basalt law of Dessert et
  al. (2003) at 26 °C and the rivers' runoff (349 mm/yr from the geography
  domain's first estimate), the net sink is 0.36 Gt CO2 a year, and up to about
  3.6 Gt if the fine, glassy regolith weathers ten times faster than catchment
  rock. Unreturned, that takes the air from 49 to 28 Pa in about 200–2,100 years.
- **Nothing returns it naturally.** Earth's volcanoes and plate tectonics
  return about 0.3–0.4 Gt a year; the Moon has neither.
- **The sink slows but lasts.** Weathering slows as soils form, but the regolith
  is 4–15 m deep, and a quarter of a metre of it (at about 0.19 kg CO2 per kg,
  from Apollo soils' CaO and MgO) would bind the whole inventory.

So the settled Moon needs one managed step in place of a volcano: recover the
carbonate sediments and heat them to drive the CO2 back out, as a lime kiln
does. That takes 128 GW of heat per Gt of CO2 a year in theory, about twice
that in practice, a sliver of the Moon's sunlight. The lime has to go into
construction or storage; spread back on soil it would take the CO2 up again.
Carbon then becomes one strand of the mineral cycle, with the calcium,
magnesium, potassium and phosphorus that weathering frees for plants.

## Limits

- Rates are Earth's: catchments, inventories and experiments. Lunar regolith has
  no soil structure or nitrogen yet, so the biological sink may start slower.
- The leaf model is textbook C3 at light saturation; C4 plants and whole
  canopies are not modelled.
- The ocean's carbonate chemistry and uptake are not included; its volume and
  alkalinity are open.
- The runoff comes from the first drainage estimate, which used climate run A.

## Next work

A box model of the air, plants, soils, ocean and carbonate sediments, with
weathering that slows as soils develop and a managed return, run together with
the water inventory and the mineral cycle. It would size the return step and
show how it changes over the Moon's first thousands of years.
