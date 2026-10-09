# Photosynthetic floaters: coupled requirements and failure modes

**Large buoyant organisms remain physically worth investigating, but the earlier lift-only sizes do not establish
that they can sustain themselves.** This study adds pressure-bearing structure, gas permeability, wet tissues,
water retention, carbon allocation, hydrogen production, dark reserves, reproduction, motion and changing buoyancy.
It also incorporates the author's correction: Terluna's persistent twilight and slowly moving Sun give organisms
several possible light histories. A stationary 14.77-day dark interval is one stress case, not a universal habitat.

At 10 km, a conditional spherical reference passes mass, photon-allocation and positive-carbon gates at tested
diameters of **40, 60, 100 and 150 m**. That is a set of tested points under one imposed trait combination, not a
certified size range. Joint material properties, water supply and trim, weather survival and a complete life cycle
remain open. No species, standing population or harvest capacity is selected by this work.

Run `python -m research.studies.floater_viability.run`; check with
`python -m pytest research/studies/floater_viability`. The versioned
[product](results/floater_viability.json) binds eight producer files, six source registers and inherited inputs by
SHA-256 and records named shared constants. [Checks](checks.json) distinguish focused tests from repository gates.
The renderer is [visualization/floater_viability/plot.py](../../../visualization/floater_viability/plot.py);
generated plots stay outside Git.

## What such an organism could look like

The working biological concept is a large gas-supported body with a mostly inert, replaceable envelope, distributed
photosynthetic tissue and small, protected liquid habitats. Thin wet films, enclosed pockets, wicks and many small
reservoirs would distribute water more gently than one open pond. Microbes, algae, small grazers and decomposers
could occupy these retained-water habitats. Their food, water, gas exchange and nutrient return still have to fit
the host's accounts; a food web cannot be inferred from the envelope's volume.

Possible lobes and tension fibres divide membrane spans. Small gas cells could limit the consequences of a puncture,
but extra partitions add mass, permeability area and developmental complexity. The current sphere and ideal lobed
models do not solve that multicellular geometry. A shaded, protected gas barrier behind a tougher outer surface has
historical material precedents, but its additional layers and their biological renewal need explicit designs.

An attractive first concept is therefore a light aerial organism with small retained-water communities, modest
active biomass, a slow replacement rate and routes through suitable light and humidity. A forest, soil bed or large
open lake aloft demands a substantially different load and water budget. The model's free-water allowance of
5 kg/m² projected area is just 5 mm averaged over that area; concentrating it into pockets does not create more water.

## Environment and units

The study reads the existing sky-ship air product. Heights are above sea level. Organism-area quantities use its
projected area, `A = pi R²`; membrane area and gas-contact area are recorded separately. Diameter is `2R`.

| Height | Air density | Mean temperature | Net supported mass per m³ of lifting mixture |
|---|---:|---:|---:|
| 10 km | 1.204 kg/m³ | 14.2 °C | 1.098 kg/m³ |
| 20 km | 1.024 kg/m³ | 5.1 °C | 0.933 kg/m³ |
| 40 km | 0.734 kg/m³ | −13.9 °C | 0.669 kg/m³ |

These are inherited mean states, not proven regional flight conditions. Mean freezing occurs near 25.7 km;
high-latitude light and a globally averaged warm layer cannot simply be combined into a demonstrated warm route.
The adopted atmosphere is approximately 400 ppm CO2. Higher surface partial pressure does not mean that the
10-km habitat has 1,000 or 2,000 ppm; those concentrations are separate sensitivities.

The lifting mixture convention is 98% hydrogen and 2% ambient air by mole fraction. The inherited rounded density
gives hydrogen inventory about 1.05% below ideal `x p V / RT` at 10 km. Supported mass, inventory, trim and daughter
inflation consistently use the density convention; the discrepancy is exposed in every feasible case.

## Size helps lift but also increases structural load

Gross supported mass grows as volume, but a tall hydrogen envelope also develops a pressure difference with height:

\[
\frac{d\Delta p}{dz}=(\rho_{air}-\rho_{gas})g,
\qquad t_{sphere}=\frac{\Delta p R}{2\sigma_{allow}}.
\]

At 10 km, the hydrostatic contribution at the crown grows by about **3.53 Pa per metre of radius**. Using the crown
load over a uniform sphere is a conservative pressure-sizing screen, not a solved membrane shape. Once that head
dominates, the required thickness grows as `R²` and skin mass as `R⁴`. Bigger therefore eventually becomes heavier
per unit lift. Gust loads, partition loads, seams, defects, fatigue and biological attachment all matter.

At 1 MPa allowable strength and 1,500 kg/m³ material density, a simple sphere without other pressure or barrier
loads has a maximum wet-payload allowance of about 50.7 kg/m² using all lift, or 12.7 kg/m² if only half is used.
The earlier 100-kg/m² wet-biomass example cannot fit that architecture at any radius. Lobes reduce local film spans;
their tendons still carry global forces. The study pays tendon mass but assumes their strength separately.

Measured biopolymer/cuticle coupon strengths make MPa-scale studies reasonable. The assembled long-lived wet
envelope's allowable strength is not known. The reference's **5 MPa** is an imposed requirement, not a validated
biological construction. See [mechanics](mechanics.md) and its [source register](mechanics_sources.json).

## There is a real biological gas-barrier precedent

Historical airships used processed intestinal collagen membranes, called goldbeater's skin, as hydrogen gas-cell
barriers. NACA's 1922 report describes finished fabrics around 130–150 g/m² and low hydrogen transmission, with
cotton, glue, glycerol and varnish among the construction materials. Many cells were protected inside a hull.
This is a macroscopic biogenic barrier precedent; living growth, repair and exposed weathering were not demonstrated.

Other primary studies report low hydrogen permeability in cellulose-derived films. The 2022 cellophane result
used here is about **0.173 Barrer**; a historical film gives a lower sensitivity around 0.039 Barrer only if its
reported litres are interpreted at standard conditions. The original litre reference and test humidity are not
established. Recent microfibrillated-cellulose values around 10–12 Barrer show that the material name alone does
not specify performance. Nanocellulose fuel-cell experiments give useful hydrated-material evidence, but their
80 °C, 95% RH device results cannot be transplanted directly into ambient lunar permeability.

Hydrogen leaks in response to its **partial-pressure difference**, close to one atmosphere here. The small total
mechanical overpressure is not the permeation driving pressure. Pinholes, seams and wet-state changes can overwhelm
an intact-film estimate. Making a barrier thicker reduces intact permeation and adds mass and replacement carbon.
Growing larger lengthens the gas-inventory turnover time, but does not reduce replacement power per projected area
when barrier thickness and gas-contact-area ratio are fixed.

The product compares measured-condition leads and explicit sensitivities. It does not combine their best values
into a claimed existing material. [Envelope analysis](envelope.md), [sources and access scope](envelope_sources.json).

## Reference requirements and joint budgets

The illustrative reference combines 5 MPa allowable structure, 25-µm minimum structural film, an additional
100-µm barrier at 0.173 Barrer, 1,500 kg/m³ material density, one great-circle gas/air partition, 20 Pa bottom
overpressure and a 10 m/s relative-gust pressure case. Full-capacity pressure sizing is paid before nominal trim.
The whole exterior plus partition is treated as gas-contact area, an explicit conservative-area sensitivity.

Per projected m² it carries 1 kg active host dry matter, 0.1 kg community dry matter, 90% water in those tissues,
5 kg free water, and hydrated dark reserves. The structural material is treated as 90% dry and mostly inert.
At most 70% of full supported lift is used, leaving 30% capacity for requirements not yet designed. This spare
capacity is not an actuator and is not counted as an upward force at nominal neutral buoyancy.

The carbon case starts with 2 kg C/m²/year gross photosynthesis, maintenance of 0.003 kg glucose/kg host dry
matter/day, quarter-rate dark maintenance, 0.1 kg C/m²/year consumer food, annual replacement of 10% of inert
structure and 50% of active host tissue. Construction costs 1.39 kg glucose/kg new dry material, an explicitly
transferred biological analogue. Community renewal must fit its aggregate food allowance; it is not a resolved
trophic model. Higher active-tissue maintenance, material turnover or wet loads can reverse the conclusion.

Hydrogen uses an assumed **1% full-solar-to-H2 LHV efficiency** and supplied cycle-mean irradiance of 200 W/m².
Daytime gas production diverts collecting area from carbon fixation. Dark fermentation uses 2.1 mol H2 per mol
glucose, with gross substrate committed and no unproven coproduct recovery credited. Measured fermentation also
produces mixed gases; purification, oxygen management and the organs' mass remain unresolved. Laboratory pulsed
photosynthetic-H2 efficiencies have different spectral and temporal denominators and do not validate the 1% input.

| Material sensitivity, fixed half-cycle dark stress | Diameter | Fraction of full lifting volume needed | Remaining construction assimilate | Parent-funded equal-daughter allocation time |
|---|---:|---:|---:|---:|
| 0.173 Barrer, 100 µm | 40 m | 0.639 | 1.029 kg C/m²/year | 7.04 years |
| 0.173 Barrer, 100 µm | 100 m | 0.336 | 0.735 kg C/m²/year | 16.14 years |
| 0.039 Barrer, 25 µm | 40 m | 0.619 | 1.076 kg C/m²/year | 6.32 years |
| 10 Barrer, 100 µm | 40 m | 0.747 | −9.668 kg C/m²/year | No surplus |
| 10 Barrer, 1,000 µm | 100 m | 0.432 | −0.503 kg C/m²/year | No surplus |

The 40-m reference supports about 18.70 kg/m² (23.5 tonnes in total); the 100-m one supports 24.59 kg/m²
(193 tonnes). Most of the volume is gas. These masses include the listed stocks but exclude unresolved organs,
control systems and their loads. The larger reference has better fractional lift reserve yet a larger structural
and daughter-construction bill. With 10 kg/m² of active dry host tissue, none of the tested sizes passes the three
joint mass, photon and carbon gates. Inert, dormant and actively respiring dry mass cannot be interchanged.

The daughter time pays construction and initial reserves plus extra hydrogen photons diverted from the parent's
fixed collecting area. It assumes no daughter photosynthesis, developmental upkeep, mortality or attachment cost;
it is a conditional allocation time, not a universal minimum or population doubling time. Recruitment and lifetime
must yield at least one successful replacement per adult. A large body that merely remains aloft is not yet a species.
See [biology](biology.md), [sources](biology_sources.json) and the explicit sensitivities in the product.

## Twilight and moving with the Sun

The inherited clear-ground optical calculation reaches PAR of at least 1 µmol photons/m²/s throughout the equinoctial cycle
poleward of about 60.4°; allowing the project's solar declination range shifts the all-season geometric threshold
to about 61.9°. That PAR threshold is an optical definition used in this comparison. It is not a universal
photosynthetic compensation point. Clouds, elevated-observer spectra and local temperature need their own treatment.

At 60° latitude, an inherited ground-canopy example spends only about 34 hours per cycle below that optical
threshold, yet its whole-plant carbon deficit lasts about 307–309 hours. The canopy's physiology is not assigned to
floaters; it demonstrates why continuous twilight and positive net carbon are different conditions. Its diffuse
spectral proxy gives PAR around 5.44 W/m² with the Sun 10° below the horizon and 1.11 W/m² at −20°. Lux alone does
not set hydrogen-production power. Cached spectral transport needed for direct 10/20-km PAR evaluation is absent,
so the study makes no aerial-light extrapolation or new climate run.

The westward ground speed that maintains mean local solar time is

\[
v_\odot=\frac{2\pi(R_{Moon}+h)\cos\phi}{T_{synodic}}.
\]

At 10 km it is **4.30 m/s at the equator, 2.15 m/s at 60°, and 0.747 m/s at 80°**. This is a slow moving target.
An organism carried by a matching wind needs no steady horizontal propulsion to keep pace. Real trajectories must
match the wind vector, including north/south drift. Selecting among altitude winds could help when the desired
velocity lies between accessible winds; changing altitude and surviving those layers still have costs.

For the spherical drag comparison, chemical input per projected area is `rho Cd v_rel³ / (2 eta)` with `Cd=0.47`
and `eta=0.25`. A **0.5 m/s** mismatch costs **0.141 W/m²**, **1 m/s** costs **1.13 W/m²**, and **2 m/s** costs
**9.05 W/m²**. A streamlined `Cd=0.05` case is a separate architectural sensitivity. The ground speed above is
only an airspeed in explicitly calm air. Sparse corrective activity cannot cheaply cancel a steady mismatch:
concentrating the same displacement into duty fraction `d` raises mean cubic-drag cost as `1/d²`.

Drifting west slightly too slowly can lengthen both daylight and darkness; it is not automatically Sun-following.
The product evaluates fixed-half-cycle, shorter-deficit and continuously useful-light cases separately, keeping
annual GPP and mean irradiance fixed to isolate timing. With these fixed inputs the reference's continuous-light
case has slightly less surplus (0.984 versus 1.029 kg C/m²/year), because the imposed active maintenance replaces
dark idling. This does not predict that Sun-following reduces real productivity: increased usable light and its
photosynthetic response have not been solved. [Photoperiod](photoperiod.md) and [navigation](navigation.md).

## Water supply and changing buoyancy are major remaining gates

At 10 km, 60% relative humidity, leaves 5 K warmer than air and 400 ppm CO2, the reference diffusion calculation
uses about **458 kg water/m²/year** for its net leaf assimilation. That is about 37 kg per lunar cycle, compared
with 5 kg/m² free-water storage. Higher humidity, cooler leaves, different carbon physiology and recapture could
reduce losses. Those traits must also preserve photosynthesis and heat balance. Night cuticular loss and metabolic
water have separate entries; the model does not pretend this diffusion estimate is a measured annual water demand.

Fog capture needs liquid-water concentration, residence time and relative droplet flow. Free drift in a uniform
wind supplies no horizontal collector flow. Existing Eulerian cloud fractions cannot establish a floater's water
route. Rain is both a source and a transient overload: retaining an extra 1 kg/m² gives a 2.38 m/s residual
terminal-descent diagnostic in the chosen drag model, before shedding or control. This is not a storm trajectory.

Water changes lift requirements. The signed trim ledger combines water intake, evaporation, drainage, metabolic
water, dry matter, leakage, production and venting before any control response. At 10 km, one kilogram of the
reference hydrogen inventory corresponds to about **13.40 kg of supported mass**. If an organism loses 1 kg water
per projected m² and holds altitude by venting gas, then restores both after collecting water, each daily cycle
requires 0.0746 kg H2/m², equivalent to **103.6 W/m² of chemical hydrogen output** averaged over that day. That
strategy overwhelms the reference's 2 W/m² full-surface H2 production capacity. It is not a tax on all gross
evaporation: simultaneous intake can balance it, and the product tests that cancellation explicitly.

Retaining gas and compressing it has a different cost. From a neutral supported load of 20 kg/m², accommodating
1 kg/m² loss by shrinking the retained gas volume at fixed altitude raises its pressure by about **4.78 kPa**.
A fixed outer hull with compressed air admitted to a common-pressure ballonet is another arrangement and has a
different pressure increment. Both require additional pressure-bearing structure or storage; a flexible
equal-pressure air/gas partition alone does not remove buoyancy while retaining all hydrogen.

Allowing altitude changes, timing uptake, suppressing evaporation and exchanging stored water are possible design
routes. A soft pressure-balanced gas bag with gas temperature following air has nearly altitude-independent lift
until its expansion limit; losing water need not produce a small, self-correcting climb. A trajectory/controller
model is needed. Night leakage-only ballast numbers likewise hold other payload mass constant: evaporation and
respiration may offset that loss, while rain and uptake may increase it. Metabolic product fates must be tracked;
carbon disappearance alone does not determine total mass change. [Environment](environment.md),
[trim cycle](trim_cycle.md), and their conservation tests make these boundaries explicit.

## What would decide the question next

1. **A joint wet membrane test:** hydrogen and oxygen/air permeability, allowable sustained stress, water sorption,
   seams, fatigue, puncture/repair and replacement rate for one assembled living-compatible construction. The
   central requirement is joint performance at its operating humidity and temperature.
2. **A coupled water and flight trajectory:** light spectrum, clouds, precipitation, temperature, wind vector and
   shear along accessible routes; actual gas geometry, water buffering, altitude control and propulsion. This is
   where the author's twilight and Sun-following direction can be tested quantitatively.
3. **A complete host metabolism and development:** shared photons for carbon and H2, gas separation, oxygen safety,
   active versus inert tissue, nutrient return, reserve turnover, daughter growth and successful recruitment.
4. **Population and ecosystem closure:** storm and collision loss rates, disease, predation, surface exchange,
   phosphorus return and shared light with underlying habitats. No harvest is credited before that account closes.

The current result narrows an architecture and exposes falsifiable requirements. It supports further work on
lightweight floaters with retained-water communities, while leaving their sustainable size and abundance open.
