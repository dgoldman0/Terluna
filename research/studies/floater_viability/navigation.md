# Following the slow Sun: wind access and locomotion costs

Author direction, 2026-10-09: colonies may follow the slowly moving daylight and much floating habitation may roam. This is a useful possibility. A wind carrying a colony along the desired route needs no steady horizontal propulsion in the air frame. Whether such winds are reachable, persistent and compatible with water, temperature, collisions and storms is unresolved. The inherited atmosphere provides useful densities and scalar wind statistics, not a navigable vector trajectory.

This module is a conditional calculation supporting the coupled floater study. It supplies no evolved propulsor, sensory system, flight controller or meteorological route. The twilight and solar geometry are owned by [photoperiod.md](photoperiod.md); the gas, carbon, water and mechanics ledgers remain separate costs that the integrated runner must pay.

## Ground motion and motion through air

For mean synodic period \(P\), lunar radius \(R_M\), height \(h\), and latitude \(\phi\), the westward ground speed preserving local solar longitude is

\[
v_\odot={2\pi(R_M+h)\over P}\cos\phi.
\]

At 10 km this is **4.303 cos(latitude) m/s**. The function comes from the shared-constant photoperiod module. East and north are positive, so desired horizontal ground velocity is \((-v_\odot,0)\). The required air-relative velocity is the vector

\[
\mathbf u=\mathbf v_{\rm ground}-\mathbf w_{\rm air}.
\]

A westward wind of 4.303 m/s at the equator matches the mean solar motion; an eastward wind of the same magnitude requires 8.606 m/s relative motion. A matching zonal component with 1 m/s northward wind still requires a 1 m/s correction. Neither a global median wind magnitude nor a percentile establishes these directions, their persistence or their accessible altitude.

Exact phase lock is sufficient for constant solar elevation in the mean fixed-declination model. A colony may instead tolerate drift inside a broad daytime corridor. Slower-than-Sun westward motion can lengthen both the next day **and the eventual night**; it does not by itself eliminate darkness. The companion photoperiod model quantifies this distinction. Real illumination also changes with declination, clouds, terrain, altitude and the atmosphere.

## Drag work and photosynthetic accounting

Let \(A_D\) be the frontal area used with drag coefficient \(C_D\), \(A_C\) the projected collecting footprint, and \(\eta\) the overall conversion from supplied energy to useful drag work. Then

\[
{D\over A_C}=\tfrac12\rho C_D{A_D\over A_C}u^2,\qquad
{P_{\rm useful}\over A_C}=\tfrac12\rho C_D{A_D\over A_C}u^3,\qquad
P_{\rm input}=P_{\rm useful}/\eta.
\]

NASA's drag equation requires a consistent reference area. Its sphere discussion also warns that Reynolds number and roughness change the coefficient. Here **Cd 0.47 is a sphere sensitivity**, not a validated universal value. **Cd 0.05 is a hypothetical streamlined comparison**. Both use area ratio 1 for the table; a real elongated, leaf-bearing body must be re-evaluated with its frontal area, collecting geometry, skin friction, attachments and orientation. More internal leaf surface does not multiply externally intercepted sunlight.

The biological conversion efficiency **0.25 is assumed**, combining chemical-to-actuator and propulsive losses. It is not a measured organism capability. The engineering comparison **0.7 refers to a supplied electrical/mechanical energy chain**, not a biological conversion demonstrated at that value. These numbers omit propulsor construction and upkeep, induced wake losses beyond the imposed efficiency, acceleration, turns, vertical travel and storm maneuvers. Buoyancy itself need not require continuous propulsive lift.

In the inherited 10 km air, density is **1.204 kg/m³**. Per collecting m²:

| Continuous relative speed | Sphere useful mechanical work | Sphere chemical input, η0.25 | Annual motion carbon equivalent |
|---:|---:|---:|---:|
| 0 m/s | 0 W | 0 W | 0 kg C |
| 0.5 m/s | 0.0354 W | 0.1415 W | 0.1116 kg C |
| 1 m/s | 0.2829 W | 1.1318 W | 0.8929 kg C |
| 2 m/s | 2.2635 W | 9.0541 W | 7.1431 kg C |

The carbon equivalent uses the study's **18 MJ/kg dry biomass and 0.45 carbon fraction**, or 40 MJ/kg C. This is a common accounting basis, not a measured muscle-substrate yield. The reference gross fixation of 2 kg C/m²/year corresponds to only **2.535 W/m²** on this basis, before maintenance, gas renewal, construction, consumers and reproduction. Thus sustained 2 m/s spherical locomotion cannot be paid from that reference photosynthesis. Even 1 m/s can exhaust a favorable remaining construction budget. Apply motion cost to the remaining assimilate once; do not subtract maintenance a second time or spend the same surplus on both locomotion and growth.

For a margin of 1 W/m² already available after other obligations, sphere correction at 0.5 m/s leaves 0.859 W/m²; 1 m/s leaves −0.132 W/m². Doubling radius does not fix this particular per-area shortage when coefficient, area ratio and productivity are unchanged: collecting power and drag work both scale as area. Size still changes Reynolds number, structure, gas inventory and control mass.

## How close must the wind match be?

Spending an available continuous input margin \(Q\) gives

\[
u_{\max}=\left({2\eta Q\over\rho C_D(A_D/A_C)}\right)^{1/3}.
\]

This is the radius of an acceptable region in **wind-vector space**, not a tolerance on wind magnitude alone. In the 10 km air, with area ratio 1:

| Available input per collecting m² | Sphere, biological η0.25 | Streamlined hypothesis, biological η0.25 | Sphere, engineering η0.7 | Streamlined hypothesis, engineering η0.7 |
|---:|---:|---:|---:|---:|
| 0.1 W | 0.445 m/s | 0.940 m/s | 0.628 m/s | 1.325 m/s |
| 0.5 W | 0.762 m/s | 1.607 m/s | 1.073 m/s | 2.266 m/s |
| 1 W | 0.960 m/s | 2.025 m/s | 1.352 m/s | 2.854 m/s |
| 2 W | 1.209 m/s | 2.552 m/s | 1.704 m/s | 3.596 m/s |

The engineering columns require their own electrical supply; they cannot be directly debited from a biological carbon budget without its conversion efficiency and hardware. All columns spend the entire specified margin on propulsion and leave no new margin for contingencies.

## Calm-air comparison and the value of daylight

Only in this explicitly **calm-air benchmark** is solar-following ground speed equal to airspeed. The chemical input below uses η0.25 and does not credit any additional fixation from daylight:

| Latitude | Required westward speed at 10 km | Sphere chemical input | Streamlined hypothesis chemical input |
|---:|---:|---:|---:|
| 0° | 4.303 m/s | 90.18 W/m² | 9.594 W/m² |
| 30° | 3.727 m/s | 58.57 W/m² | 6.231 W/m² |
| 60° | 2.152 m/s | 11.27 W/m² | 1.199 W/m² |
| 70° | 1.472 m/s | 3.608 W/m² | 0.3838 W/m² |
| 80° | 0.747 m/s | 0.4722 W/m² | 0.05023 W/m² |
| 85° | 0.375 m/s | 0.05970 W/m² | 0.00635 W/m² |

High latitude reduces travel speed but changes illumination. The inherited **ground, horizontal-plane, equinox-noon** PAR comparator is 334.0 W/m² at latitude 0°, 133.3 at 60°, 35.05 at 80°, and 16.08 at 85°. These are not aerial absorbing-surface fluxes or bolometric fluxes, and are not converted here into a floater's carbon output. A colony's orientation, direct/diffuse light collection and light-response curve matter. A universal 200 W/m² illumination budget across all latitudes would be unjustified.

For a realized additional chemical fixation rate \(\Delta G\), solar tracking has an energetic advantage only if

\[
\Delta G>P_{\rm locomotion}+\Delta P_{\rm other}.
\]

`sunlight_advantage` makes this comparison without manufacturing a light-to-carbon conversion. Extra costs can include increased daytime maintenance, water supply and cooling; benefits can include a smaller dark reserve and less dark hydrogen-production substrate. Continuous illumination need not double fixation: light saturation, carbon supply, water and biological rhythms remain. Existing respiration is already paid in the baseline ledger. Twilight can shorten or remove optically dark intervals without supplying enough light for net carbon gain; the companion photoperiod analysis distinguishes the quantities.

## Steering with winds is an open physical route

The Bellemare et al. 2020 primary balloon study includes an actual **39-day Pacific experiment** with engineered superpressure balloons. Publisher descriptions show navigation using wind information and opposing wind sheets. This is evidence that autonomous wind-aware navigation can work in Earth's stratosphere. It supplies neither the lunar wind field nor a biological controller or free altitude changes. The study's complete methods were not available through the retrieved publisher page; that access limit is recorded in the source register.

For an optimistic local screen, two accessible constant wind vectors \(\mathbf w_1,\mathbf w_2\) allow an average

\[
\overline{\mathbf w}=(1-f)\mathbf w_1+f\mathbf w_2,\quad 0\le f\le1.
\]

`two_layer_advection` finds the point on that velocity segment closest to the desired ground vector. Imposed westward layers at 2 and 6 m/s can average 4.303 m/s with **57.58%** of time in the faster layer and no steady horizontal propulsive work in this idealization. If both layers also carry a 1 m/s northward component, that residual cannot be canceled by switching between them. The example is not extracted from the inherited scalar wind statistics.

Altitude changes have finite time, alter buoyancy/trim and may require pumping, gas production, gas venting, thermal control or ballast exchange. Higher and lower layers can differ in light, humidity, liquid water, freezing and shear. Whether recurrent vertical steering is economical therefore requires a coupled trajectory and control mechanism; it cannot be inferred from the gross gravitational energy of the whole displaced body alone. Wind shear and buoyancy cycles may supply useful mechanical energy, but no extraction mechanism or energy credit is assumed here. An isolated floater comoving with uniform wind has no free relative airflow to power a turbine or collect fog by sweep-through.

Occasional powered collision avoidance, storm escape or correction is compatible with roaming. For an imposed powered speed and duty, mean expenditure is duty times the active cubic power, but this does not prove solar tracking between maneuvers. In constant wind, supplying a constant mean correction \(u\) with powered fraction \(d\), coasting comoving during the rest, requires active airspeed \(u/d\). Even ignoring accelerations, mean cubic drag cost rises by **1/d²** over continuous gentle correction. A mean 0.5 m/s correction costs 0.1415 W/m² continuously, 0.5659 W/m² at 50% duty and 14.147 W/m² at 10% duty. Broad permitted drift, useful changing winds and altitude steering can change the optimal strategy; arbitrary duty reduction cannot make a persistent adverse wind free.

## What this calculation establishes

Sun-following is kinematically slow, and wind-matched passive travel is a legitimate favorable case. Continuous whole-body propulsion across adverse winds can consume more energy than a photosynthetic colony fixes. The unresolved question is which wind corridors and biological/control mechanisms allow low residual airspeed while preserving water, carbon, gas, storm and collision budgets. Scalar wind statistics cannot decide that question in either direction. No new GCM or CFD run was made.

`navigation.py` exposes relative velocity, drag power, wind tolerance, carbon-equivalent motion budget, conditional solar-following cases, a two-layer advection bound, pulse-correction accounting and sunlight advantage. `evaluate()` reads the versioned sky-ships product and records input hashes. `test_navigation.py` checks vector invariance, force-times-speed work, cubic and area scaling, inverse power bounds, time/displacement conservation, carbon accounting and product provenance. Full integrated-run validation belongs to the parent study runner.

Primary and official sources: [NASA drag equation](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-equation/); [NASA sphere drag](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/drag-of-a-sphere/); [Colozza and Dolce 2005](https://ntrs.nasa.gov/api/citations/20050080709/downloads/20050080709.pdf); [Bellemare et al. 2020](https://www.nature.com/articles/s41586-020-2939-8). Their access conditions, scope and exclusions are preserved in [navigation_sources.json](navigation_sources.json). Colozza and Dolce's engineering assessment motivates joint solar, propulsion, aerodynamic and environmental sizing; none of its numerical hardware performance is transferred to the biological cases.
