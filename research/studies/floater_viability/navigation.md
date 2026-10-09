# Holding the Sun: wind and the cost of motion

**Holding the Sun costs nothing when a matching wind carries the floater and up to 90 W/m² when it must push a sphere
through calm air at the equator.** The westward ground speed that holds mean solar time is 4.303 cos(latitude) m/s at
10 km. In calm air a sphere (C_d 0.47) driven at that speed by muscle-like actuators at 25% efficiency spends
90.18 W/m² at the equator and 11.27 W/m² at 60°; the reference's whole gross fixation is worth 2.535 W/m². The answer
is to ride winds, sail and pick latitudes, which [sailing](sailing.md) maps on the project's winds. The author's
direction of 9 October 2026 is that colonies may follow the slow daylight and that much floating habitation may roam.

## Ground motion and airspeed

**The air-relative velocity is u = v_ground − w_air.** East and north are positive, so holding the Sun needs ground
velocity (−v_☉, 0). A westward wind of 4.303 m/s at the equator matches it; an eastward wind of the same speed needs
8.606 m/s of airspeed; a matching zonal wind with 1 m/s northward needs a 1 m/s correction. A broad daytime corridor can
replace exact phase lock.

## Drag work

**Pushing a body at airspeed u costs (1/2)ρ C_d (A_D/A_C) u³/η per square metre of collecting area.** The drag uses one
consistent reference area (NASA drag equation), and the sphere's coefficient depends on Reynolds number and roughness
(NASA sphere drag). C_d 0.47 is the sphere case and 0.05 the streamlined one, both at area ratio 1. The efficiency of
0.25 bundles chemical-to-actuator and propulsive losses for a biological mover; 0.7 is an engineering electrical
chain. At 10 km, density 1.204 kg/m³:

| Airspeed | Useful work, sphere | Chemical input at η 0.25 | Carbon a year |
|---:|---:|---:|---:|
| 0.5 m/s | 0.0354 W/m² | 0.1415 W/m² | 0.1116 kg C/m² |
| 1 m/s | 0.2829 W/m² | 1.1318 W/m² | 0.8929 kg C/m² |
| 2 m/s | 2.2635 W/m² | 9.0541 W/m² | 7.1431 kg C/m² |

The carbon equivalent uses 18 MJ/kg dry biomass at 45% carbon, 40 MJ/kg C. A steady 2 m/s through still air costs more
than the reference fixes; 1 m/s takes most of a favourable surplus. Drag work and collecting power both scale with area,
so size alone does not change these per-area costs.

**The wind must match within a small circle in wind-vector space.** Spending a margin Q on drag allows a mismatch of
u_max = (2ηQ/(ρ C_d A_D/A_C))^(1/3):

| Margin per collecting m² | Sphere, η 0.25 | Streamlined, η 0.25 | Sphere, η 0.7 | Streamlined, η 0.7 |
|---:|---:|---:|---:|---:|
| 0.1 W | 0.445 m/s | 0.940 m/s | 0.628 m/s | 1.325 m/s |
| 0.5 W | 0.762 m/s | 1.607 m/s | 1.073 m/s | 2.266 m/s |
| 1 W | 0.960 m/s | 2.025 m/s | 1.352 m/s | 2.854 m/s |
| 2 W | 1.209 m/s | 2.552 m/s | 1.704 m/s | 3.596 m/s |

## Calm air and the worth of daylight

In calm air the ground speed is the airspeed. At η 0.25:

| Latitude | Westward speed at 10 km | Sphere | Streamlined |
|---:|---:|---:|---:|
| 0° | 4.303 m/s | 90.18 W/m² | 9.594 W/m² |
| 30° | 3.727 m/s | 58.57 W/m² | 6.231 W/m² |
| 60° | 2.152 m/s | 11.27 W/m² | 1.199 W/m² |
| 70° | 1.472 m/s | 3.608 W/m² | 0.3838 W/m² |
| 80° | 0.747 m/s | 0.4722 W/m² | 0.05023 W/m² |
| 85° | 0.375 m/s | 0.05970 W/m² | 0.00635 W/m² |

High latitudes cost little to hold but give less light: the ground's equinox-noon PAR is 334.0 W/m² at the equator,
133.3 at 60°, 35.05 at 80° and 16.08 at 85° ([photoperiod](photoperiod.md)). Following the Sun pays only if the extra
fixation ΔG exceeds the motion and any other added costs; `sunlight_advantage` makes that comparison once the light
gain is known.

## Steering with winds

**Wind-aware altitude changes can steer a balloon for weeks.** Loon's engineered superpressure balloons navigated a
39-day Pacific flight with wind information, using opposing wind layers (Bellemare et al. 2020). Two steady layers w₁
and w₂ give any average on the segment between them: westward layers of 2 and 6 m/s average 4.303 m/s with 57.58% of
the time in the faster one, and a shared 1 m/s northward component stays. `two_layer_advection` finds the nearest
average; altitude changes cost trim, gas and time ([trim cycle](trim_cycle.md)).

**Short bursts of thrust cost more than steady correction.** A mean correction u delivered over a duty d needs u/d of
airspeed while active, so the cubic drag cost rises as 1/d²: 0.5 m/s costs 0.1415 W/m² steadily, 0.5659 W/m² at 50%
duty and 14.147 W/m² at 10%.

## Checks and sources

[test_navigation.py](test_navigation.py) checks vector invariance, force-times-speed work, cubic and area scaling, the
inverse power bound, time and displacement conservation, the carbon accounting and the product's provenance.

Sources: [NASA drag equation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/);
[NASA sphere drag](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-of-a-sphere/);
[Colozza & Dolce 2005](https://ntrs.nasa.gov/api/citations/20050080709/downloads/20050080709.pdf), whose airship
assessment couples solar collection, propulsion, aerodynamics and the wind environment;
[Bellemare et al. 2020](https://www.nature.com/articles/s41586-020-2939-8). Access and scope are in
[navigation_sources.json](navigation_sources.json).
