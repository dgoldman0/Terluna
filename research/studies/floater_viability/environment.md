# Water and heat for photosynthetic floaters

The committed atmosphere supports a warm lower-sky investigation. It also makes continuous water replacement,
solar leaf heating, rain shedding and cold-altitude survival explicit requirements. Atmospheric water vapour is
not an inventory of immediately usable liquid water. These calculations do not establish a viable organism.

`environment.evaluate()` returns the environmental component of the parent study. `INPUT_FILES` identifies the
committed products it reads. Functions accept explicit units and independently specified conditions; tests check
diffusion units, water and energy conservation, limiting cases and the inherited air. No new GCM campaign,
trajectory simulation, external spectral restoration or weather-field generation was performed.

## Which environment is actually available

The [sky-ship product](../sky_ships/results/sky_ships.json) supplies the following air, based on the design GCM's
three-day output means. Heights are above model sea level. These are global profile means, not instantaneous
weather along a floater's path.

| Height | Mean temperature | Pressure | Approximate CO2 pressure at inherited400ppm |
|---|---:|---:|---:|
| 10km | 14.2°C | 0.988atm | 40.04Pa |
| 20km | 5.1°C | 0.813atm | 32.95Pa |
| 25km | 0.7°C | 0.735atm | 29.79Pa |
| 30km | −3.8°C | 0.665atm | 26.95Pa |
| 40km | −13.9°C | 0.542atm | 21.97Pa |

The ecology product's mean freezing height is25.7km. Adequate gas lift at40km therefore does not establish
unprotected, continuously active wet photosynthetic tissues there. Cold-tolerant dormant tissues, thaw cycles or
protected internal habitats are different candidate designs. The favourable temperature screen points first to
roughly10–20km, while water access, storm exposure and solar heating still require evaluation.

The canopy product records48.636Pa CO2 on the Moon and40.53Pa on Earth: both use400ppm, with the difference due to
surface pressure. The atmospheric-CO2 study and the latest upper-air cooling calculations also use400ppm. The
1000/2000ppm cases here are physiological sensitivities, not adopted compositions or climate calculations.

The equatorial CM1 ring's Eulerian cloud fractions at10km are3.4% by day and0.7% by night; at20km2.4% and0.97%;
at40km0.27% and0.011%. Cloud presence does not specify liquid content, droplet size or the water a moving organism
can capture. Cloud fraction is not a parcel's wet duty. The ring's water-vapour residence time of97.6days is a
column water-budget measure, not a cloud lifetime or a floater's hydration guarantee. Surface relative humidity
cannot be assigned to the10–40km band. The calculations therefore expose60/90/98% RH as assumptions.

PlaSim's cloud water is diagnosed for radiation from an Earth-derived formula; its condensation/convection scheme
rains out excess vapour. The implementation in `climate/gcm/cloud_water.py` does not supply a prognostic reservoir
of collectable cloud liquid. It is not used as such here.

Existing ground optics give1,555µmol PAR/m²/s at overhead Sun. The plant calculation's ground twilight remains
above1µmol/m²/s for roughly58hours after sunset at its equatorial reference. Neither number is an aerial light
field. Altitude changes the direct horizon, scattering, clouds and infrared surroundings. This component takes
no automatic production or night-storage credit from altitude.

## A wet body needs a water flux

Per square metre of one-sided leaf area, the dilute-vapour diffusion screen uses

\[
g=(g_{leaf}^{-1}+g_{boundary}^{-1})^{-1},\qquad
E=g\,[e_s(T_{leaf})-RH\,e_s(T_{air})]/P.
\]

Conductances are mol/m²/s, vapour pressures and total pressure are Pa, and E is mol water/m²/s. Negative E is
flagged as potential condensation and is not credited as biological water uptake. Multiply by0.01801528kg/mol
and86,400s/day for water flow. The associated latent cooling uses2.45MJ/kg. Saturation follows the updated Buck
liquid-water equation documented by [NCAR EOL](https://doi.org/10.26023/S2SX-K278).

The [Wang et al.2024 measurements](https://doi.org/10.1093/treephys/tpae027) place minimum conductance of nine
European tree species at0.8–4.8mmol/m²/s at20°C. Their denominator is **total two-sided area**. This model doubles
those endpoints to1.6–9.6mmol/m²/s on a flat leaf's one-sided area. These Earth measurements bracket a sensitivity;
they are not established lunar traits. Closure still loses water, and ordinary night conductance can exceed the
fully closed minimum. Active conductances of0.05/0.1/0.3mol/m²/s and boundary conductance0.5 are explicit scenarios;
an experimental active value near0.1 is reported by [Gauthier et al.2018](https://doi.org/10.1104/pp.16.00741).

At10km,60%RH and leaf temperature19.2°C, g_leaf=0.1 with g_boundary=0.5 gives **1.624kg/m²leaf/day** and46.1W/m²
latent cooling. A closed leaf at14.2°C and90%RH with g_leaf=0.0016 loses0.0040kg/m²leaf/day. Assuming those two
states each occupy half the cycle, the mean is0.814kg/m²leaf/day. This illustrative cycle has no measured weather
trajectory or carbon-performance claim.

| One-sided leaf area / projected body area | Water lost per lunar cycle | Autonomy from9kg/m² tissue water,20% usable | Total tissue store needed for a whole dry cycle at20% usable |
|---|---:|---:|---:|
| 1 | 24.0kg/m² | 2.21days | 120kg/m² |
| 3 | 72.1kg/m² | 0.74days | 361kg/m² |

The20% expendable fraction is a scenario. A separate freely drainable reservoir could have a larger usable share;
its water and container both consume lift. Exhausting all tissue water is not a survival strategy. Stomatal
closure can prolong hydration while reducing carbon acquisition. Water retained in a9kg/m² standing stock cannot
be repeatedly counted as new supply.

## Couple water to actual carbon uptake

For a stomata-dominated C3 diffusion case,

\[
E/A=1.6\,{VPD\over x_{CO2}P(1-C_i/C_a)}.
\]

`water_for_assimilation` accepts annual **net leaf CO2 uptake**, per projected organism area. If the parent carbon
model provides GPP, subtract leaf respiration first. Whole-organism NPP also subtracts other maintenance and
growth costs; substituting NPP directly would undercount the CO2 transported through stomata. The output records
water per mole of net leaf C, so the parent can use its own carbon result without rerunning this environment.

For the10km,60%RH,leaf+5K case above, Ci/Ca=0.7 and annual GPP2kgC/m², a70–90% net-leaf fraction requires:

| CO2 mixing sensitivity | Water per year | Water per lunar cycle |
|---|---:|---:|
| 400ppm | 351–451kg/m² | 28.3–36.4kg/m² |
| 1000ppm | 140–180kg/m² | 11.3–14.6kg/m² |
| 2000ppm | 70–90kg/m² | 5.67–7.29kg/m² |

The5kgGPP sensitivity multiplies these flows by2.5. These estimates omit cuticular and night losses. The higher
CO2 reductions hold Ci/Ca, VPD and carbon uptake fixed; they do not prove the same productivity, nutrient supply
or climate at the richer composition. CAM, C4, microbial carbon pathways and external wet biofilms need their
own physiology. CO2 and water diffuse through boundary layers differently, so the stomatal ratio is an explicit
approximation rather than a complete leaf gas-exchange solver.

Humidity alone does not remove solar-heating costs. For1kgNPP/m²/year, net-leaf uptake twice NPP and a leaf5K above
the10km air, increasing RH from60% to90% and98% changes water demand from501 to307 and255kg/m²/year. The warmed
leaf still has a vapour-pressure deficit relative to nearly saturated ambient air.

## Harvest droplets, reject condensation heat, shed rain

Droplet interception is LWC × relative speed × capture efficiency × collector area. At0.1g liquid/m³,3m/s and30%
capture, each collector square metre gathers0.324kg/hour while in fog. At an assumed10% fog duty that is
0.778kg/day. With1% duty it is0.0778kg/day. Capture area relative to leaf/body area matters. The0.05/0.1/0.5g/m³
grid is a set of analog-informed sensitivities; the high value is an explicit extension. The source register
separates terrestrial field collectors from modeled tree interception and records incomplete source access.

A body drifting exactly with a uniform wind has no through-flow from that wind. Droplet settling, controlled
relative motion and wind shear can supply different encounter velocities, each with a load and energy account.
Following clouds may increase collection duty, but requires a viable trajectory and may expose the organism to
storms, shade and freezing. Droplet water was already condensed elsewhere; interception does not itself release
the full vapour condensation heat again. Fog depletion and changes to rain must nevertheless enter the shared
water budget.

Condensing **vapour** is different. A surface must be colder than the dewpoint and reject the latent heat. Even
with10 or50W/m² of *net* cooling available, the energy ceilings are0.353 or1.76kg/m²/day; actual mass transfer can
be lower. These are not guaranteed daily yields. Absorbing vapour through hygroscopic material also requires
water-activity, heat and regeneration accounts.

The current committed CM1 ring records sampled storm rain maxima of3.71mm/hour at the median,18.13 at p90 and
68.52 at its largest sampled event. Those are a **surface-rain stress proxy**, not precipitation flux measured
onto a floater at altitude. At full interception the last value adds6.85kg/m² in six minutes or68.5kg/m² in one
hour unless drained. Routing water inward and shedding excess need distinct capacities. Icing, hail, retained
surface films and nutrient washout remain additional loads. Storm tracks last6hours at the median,15 at p90 and
57 maximum at three-hour sampling; they do not establish a safe floater route or avoidance schedule.

## The solar and night heat limits

The two-face leaf screen balances prescribed absorbed radiation with sensible heat, evaporation and infrared
exchange with an effective environmental radiation temperature. Its2/10/30W/m²/K sensible coefficients include
both faces and are scenarios, not consequences of the kilometre-scale body size. Its IR surroundings are explicitly
specified; a leaf in dense air is never silently assigned a view of3K space. Absorbed100/300/600W/m² spans thermal
requirements without claiming an altitude-resolved solar spectrum. Cloud radiation, body conduction and actual
leaf orientation require a resolved model.

At10km with air and effective IR surroundings both14.2°C,60%RH and300W/m² absorbed, the10W/m²/K case gives:

| Leaf state | Equilibrium leaf temperature | Water loss | Latent cooling |
|---|---:|---:|---:|
| Minimum conductance0.0016 | 28.31°C | 0.0714kg/m²/day | 2.02W/m² |
| Active conductance0.1 | 24.76°C | 2.789kg/m²/day | 79.08W/m² |

Closing stomata saves water but warms the leaf. Under600W/m² and the weak2W/m²/K sensible exchange, the minimum-
conductance leaf reaches54.3°C, while the active leaf reaches40.3°C and loses8.45kg/m²/day. Conductances here are
held fixed; heat injury, stomatal response and temperature-dependent cuticular loss are not modeled. Thus these
are balance solutions and thermal demands, not a claim that the tissues remain functional.

At40km, keeping exposed closed-stomata leaf tissue at5°C with IR surroundings at the mean−13.9°C requires199,
350 or728W/m² under the2/10/30W/m²/K sensible cases. This makes sustained warm exposed tissue through the night
demanding. Protected, insulated internal liquid stores and freeze-tolerant dormant tissues require different
geometries and models; they are not excluded by the exposed-leaf calculation.

Water's sensible heat alone is also finite. With a10K usable cooling range,9kg/m² stores0.376MJ/m². A persistent
1W/m² net heat deficit consumes that in4.35days;10W/m² takes0.435day. A14.77day interval at1 or10W/m² would require
30.5 or305kg water/m² for that10K range. These constant-deficit sensitivities are not a prediction of night weather:
a body can equilibrate with warm surrounding air, and insulation can reduce loss. They show why thermal mass
alone cannot justify warm internal water through arbitrary cold nights. Latent freezing heat is not credited as
continued liquid-habitat function.

The next coupled test must follow water, carbon, heat and lift along actual weather histories: leaf area and
conductance, liquid uptake and mobile stores, rain drainage, radiation, relative-flow collection and freezing
strategy. These requirements favour a manageable warm-altitude design with explicit water replenishment over
assuming that atmospheric depth and long-lived vapour alone sustain wet aerial production.
