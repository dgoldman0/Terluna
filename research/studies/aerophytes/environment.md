# Water and heat for aerophytes

**A wet aerophyte at 10 km with the reference traits transpires 1.25 kg of water per square metre a day, and a closed
leaf in full sun runs 14 K above the air.** The 10–20 km band is warm (14.2 °C at
10 km, 5.1 °C at 20 km) and above freezing, which suits wet tissue; above 25.7 km the mean air freezes. Water flow,
leaf heat, rain loads and the long night's heat set this component's requirements; [water](water.md) finds where the
air supplies the water.

[environment.py](environment.py) returns the environmental component; `INPUT_FILES` names the committed products it
reads. The tests check the diffusion units, water and energy balances, limits and the inherited air.

## The air

The [sky-ship product](../sky_ships/results/sky_ships.json) gives the design GCM's three-day means; heights are above
the model's sea level and the values are global profile means.

| Height | Mean temperature | Pressure | CO₂ partial pressure at 400 ppm |
|---|---:|---:|---:|
| 10 km | 14.2 °C | 0.988 atm | 40.04 Pa |
| 20 km | 5.1 °C | 0.813 atm | 32.95 Pa |
| 25 km | 0.7 °C | 0.735 atm | 29.79 Pa |
| 30 km | −3.8 °C | 0.665 atm | 26.95 Pa |
| 40 km | −13.9 °C | 0.542 atm | 21.97 Pa |

The mean freezing height is 25.7 km (the ecology product). At the ground the canopy model's 400 ppm is 48.6 Pa of CO₂
on the Moon against 40.5 Pa on Earth. The 1,000 and 2,000 ppm cases here are physiological sensitivities.

The equatorial CM1 ring's cloud fractions are 3.4% by day and 0.7% by night at 10 km, 2.4% and 0.97% at 20 km, and
0.27% and 0.011% at 40 km; its column water stays in the air 97.6 days. PlaSim's cloud water serves its radiation and
is rained out by its convection scheme, so this study takes cloud liquid from CM1. The ground's PAR is 1,555 µmol/m²/s
with the Sun overhead, and the plant model's ground twilight stays above 1 µmol/m²/s for about 58 hours after
equatorial sunset; aloft, the light is the illumination solver's to compute.

## Transpiration

**Per square metre of one-sided leaf, water leaves at E = g [e_s(T_leaf) − RH e_s(T_air)]/P, with
g = (1/g_leaf + 1/g_boundary)⁻¹.** Conductances are in mol/m²/s; saturation follows the updated Buck equation for
liquid water (NCAR EOL). Condensation is flagged and counted only where a collector takes it.

**Closed leaves still lose water.** Minimum conductances of nine European tree species are 0.8–4.8 mmol/m²/s at 20 °C
on two-sided area (Wang et al. 2024), 1.6–9.6 mmol/m²/s on one side; active conductance near 0.1 mol/m²/s is
measured by Gauthier et al. (2018). At 10 km and 60% RH, an open leaf (0.1) 5 K warm with boundary conductance 0.5
loses 1.624 kg/m² a day and 46.1 W/m² of latent heat; a closed leaf (0.0016) at air temperature in 90% air loses
0.0040 kg/m² a day. Half the cycle in each state averages 0.814 kg/m² of leaf a day.

| Leaf area per projected area | Water per lunar cycle | Days from 9 kg/m² tissue water at 20% usable | Tissue store for a whole dry cycle |
|---|---:|---:|---:|
| 1 | 24.0 kg/m² | 2.21 days | 120 kg/m² |
| 3 | 72.1 kg/m² | 0.74 days | 361 kg/m² |

The 20% usable share is a scenario; tissue water cannot all be spent.

## Water for carbon

**A C3 leaf spends E/A = 1.6 VPD/(x_CO₂ P (1 − C_i/C_a)) moles of water per mole of net carbon.**
`water_for_assimilation` takes the year's net leaf uptake per projected area. At 10 km, 60% RH, a leaf 5 K warm and
C_i/C_a 0.7, a gross fixation of 2 kg C/m² with 70–90% net at the leaf needs:

| CO₂ | Water a year | Water per lunar cycle |
|---|---:|---:|
| 400 ppm | 351–451 kg/m² | 28.3–36.4 kg/m² |
| 1,000 ppm | 140–180 kg/m² | 11.3–14.6 kg/m² |
| 2,000 ppm | 70–90 kg/m² | 5.67–7.29 kg/m² |

Five kilograms of gross fixation multiply these by 2.5. For 1 kg C/m² of net production, net leaf uptake twice that
and a leaf 5 K warm, raising humidity from 60% to 90% and 98% lowers the need from 501 to 307 and 255 kg/m² a year:
the warm leaf still sees a vapour-pressure deficit. Cool leaves, humid layers and CO₂-concentrating mechanisms cut it
much further ([water](water.md)).

## Droplets, condensation and rain

**Droplets come only with relative flow: 0.1 g/m³ at 3 m/s and 30% capture gives 0.324 kg/m² an hour in fog.** At 10%
or 1% fog duty that is 0.778 or 0.0778 kg/m² a day. A body drifting with a uniform wind has no flow through it;
settling, sailing and shear supply it. Droplet water was condensed elsewhere, so catching it releases no latent heat.

**Condensing vapour needs a surface below the dewpoint and the latent heat removed.** With 10 or 50 W/m² of net
cooling the energy ceilings are 0.353 and 1.76 kg/m² a day. Sorbents and their heat are costed in [water](water.md).

**The ring's storms rain 3.71 mm/h at the median, 18.13 at the 90th percentile and 68.52 at the largest sampled
event.** Fully caught, the largest adds 6.85 kg/m² in six minutes and 68.5 kg/m² in an hour unless drained; storms
last 6 hours at the median, 15 at the 90th percentile and 57 at most in the three-hourly tracks. [Storms](storms.md)
gives the shedding the body needs.

## Heat

**In full sun a closed leaf runs hotter and an open leaf runs wetter.** The two-face leaf balance at 10 km, with air
and infrared surroundings at 14.2 °C, 60% RH, 300 W/m² absorbed and 10 W/m²/K of sensible exchange:

| Leaf | Leaf temperature | Water loss | Latent cooling |
|---|---:|---:|---:|
| Closed, 0.0016 | 28.31 °C | 0.0714 kg/m²/day | 2.02 W/m² |
| Open, 0.1 | 24.76 °C | 2.789 kg/m²/day | 79.08 W/m² |

At 600 W/m² and a weak 2 W/m²/K exchange the closed leaf reaches 54.3 °C and the open one 40.3 °C while losing
8.45 kg/m² a day. The infrared surroundings are a stated temperature; the leaf never faces the cold of space through
this dense air. The exchange coefficients 2–30 W/m²/K and the absorbed 100–600 W/m² are scenarios.

**Exposed wet tissue at 40 km needs 199–728 W/m² of heat to stay at 5 °C against air at −13.9 °C**, at 2–30 W/m²/K of
exchange; warm tissue there belongs inside insulated bodies. Liquid water's sensible heat is finite: 9 kg/m² over a
10 K range holds 0.376 MJ/m², which a steady 1 W/m² deficit spends in 4.35 days and 10 W/m² in 0.435 days; a 14.77-day
interval at those deficits would need 30.5 or 305 kg/m² of water. Bodies in the 10–20 km band sit in warm air and
lose little heat; shorter nights aloft help further ([sailing](sailing.md)).

Sources: [Wang et al. 2024](https://doi.org/10.1093/treephys/tpae027), [Gauthier et al. 2018](https://doi.org/10.1104/pp.16.00741),
[NCAR EOL saturation formulations](https://doi.org/10.26023/S2SX-K278); fog-catcher context in the
[source register](environment_sources.json).
