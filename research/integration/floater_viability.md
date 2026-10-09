# Photosynthetic floaters: coupled requirements and remaining tests

The 9 October 2026 follow-up evaluates large aerial organisms as a research direction, including the author's proposals that colonies can follow the slow Sun and exploit prolonged twilight. It finds conditional combinations that pass static lift, gas-replacement and carbon calculations. It does **not** establish a self-sustaining species, a safe size range or harvest available to human settlements.

The model, data product, source registers and tests belong to [`domain/ecology`, `research/studies/floater_viability`](https://github.com/dgoldman0/Terluna/tree/domain/ecology/research/studies/floater_viability). This main-branch record summarizes that work and connects it to the [earlier aerial biosphere and population evaluation](aerial_biosphere_population.md). It imports no model into another domain and changes no population allocation. Branch checkpoint: [`f1f8383`](https://github.com/dgoldman0/Terluna/commit/f1f8383bd1aa1ba60a301b5372a037c946537360). All 101 focused floater tests pass; repository check limits are recorded below.

## Twilight and motion change the assumed night

A stationary colony does not universally face 14.77 days of biological darkness. The inherited clear-ground twilight approximation keeps PAR above 1 µmol/m²/s for about 58 hours after equatorial sunset. At latitude 60°, its plant product has only 33.7 hours per cycle below that threshold, yet roughly 309 hours when photosynthesis fails to cover the reference plant's maintenance. Optical darkness, carbon deficit and a hydrogen pathway's usable illumination are different quantities.

Applied geometrically, the same ground twilight table predicts no interval below PAR 1 poleward of about 60.37° at zero declination, or 61.91° across the stated declination sensitivity. This is a **ground-reference “eternal twilight” condition**, not guaranteed year-round productivity. At 70° equinox midnight, the reference is only about 5.14 µmol/m²/s. Clouds, aerial light transport, biological action spectra and maintenance remain unresolved. Visual lux cannot be treated as hydrogen-production power.

At 10 km, maintaining mean solar longitude requires westward **ground speed of 4.303 cos(latitude) m/s**. Air-relative motion is the desired ground vector minus the local wind vector. A matching wind therefore permits steady horizontal drift without propulsion; an equally fast wind in the opposite direction does not. Selecting accessible wind layers is a legitimate control possibility. Terrestrial engineered balloons provide an actual navigation precedent, including a 39-day controlled Pacific experiment, but the inherited lunar wind magnitudes do not establish vector routes, controllable shear or persistence.

For a spherical sensitivity with Cd 0.47 and assumed overall chemical-to-drag-work efficiency 0.25, continuous correction of **1 m/s costs 1.13 W/m²** of projected collecting area; 2 m/s costs 9.05 W/m². The reference gross fixation of 2 kg C/m²/year contains only 2.535 W/m² on the study's biomass-energy basis, before other obligations. Available chemical margins of 0.1, 0.5, 1 and 2 W/m² tolerate vector mismatches of approximately 0.45, 0.76, 0.96 and 1.21 m/s. Streamlining can help, but its assumed Cd 0.05 is unvalidated for a leaf-bearing colony.

Higher latitude reduces travel speed and changes illumination: the ground equinox-noon PAR comparator falls from 334 W/m² at the equator to 133 at 60° and 35 at 80°. These are neither aerial nor bolometric fluxes. Sun-following can reduce reserve and dark gas-production burdens, but extra fixation must pay motion and changed water, heat and maintenance costs. Drifting westward more slowly than the Sun can lengthen the eventual night as well as the day. Low maneuver duty alone does not guarantee continuous daylight.

## Buoyancy needs a pressure-bearing body

The retained sky-ship atmosphere gives mean temperatures of 14.2°C at 10 km and 5.1°C at 20 km, making the lower sky the first warm-tissue comparison. Mean air at 40 km is −13.9°C; lifting capacity there does not establish exposed wet physiology.

Large envelopes are not governed by volume-to-area lift scaling alone. For a full gas body, hydrostatic pressure contrast rises approximately as diameter times density contrast times gravity. Spherical membrane thickness then grows approximately with radius squared at fixed allowable strength. Small bodies cannot carry sufficient wet payload; sufficiently large homogeneous spheres spend excessive lift on their own pressure structure. Minimum film thickness, gusts, seams, barriers and partitions further constrain the window.

In the sampled reference sweep, a 5 MPa assumed assembled allowable and a 70% maximum lift-utilization criterion admit the mass condition at diameters 40–500 m among tested points, while 30 m and 1,000 m fail. At 1 MPa none of the sampled sizes passes that same reserve criterion. These are architecture-dependent sampled conditions, not certified continuous size ranges. Lobed films and tendons can change the scaling, while retaining global loads, attachment and stability requirements. Coupon tensile strength does not establish a giant membrane's long-term allowable.

Neutral trim uses actual gas, ambient-air ballast, structure and wet mass. Unused lifting capacity is a reserve, not an upward force continually acting on the body. A soft, pressure-matched gas bag also has no automatic stable altitude merely because ambient density decreases: expansion approximately offsets that density change until geometry or another control changes the balance.

## Biological barriers have real precedents

Hydrogen containment by biologically derived materials is experimentally grounded. A 2022 regenerated-cellulose measurement reports about 0.173 Barrer at room temperature, with test humidity unrecovered. The official 1934 cellophane study reports 0.1 L H₂/m²/day for approximately 25.4 µm material at 25°C, while the reporting-liter convention and gas-test humidity remain unverified. Historical goldbeater's-skin gas cells used processed intestinal membranes with cotton, adhesive and varnish. They demonstrate useful biological material, not a growing exposed living envelope.

Humid cellulose fuel-cell tests provide another important comparator: a 30 µm membrane at 80°C and 95% RH has measurable hydrogen crossover; a thinner spray membrane suffered cracks and pinholes. Such observations prevent assigning dry-film performance universally to wet material. Oxygen or water-vapor permeability cannot substitute for hydrogen measurements.

Leakage depends on hydrogen **partial-pressure** difference, potentially near one bar, rather than the much smaller mechanical gauge pressure. An extracellular, largely inactive barrier could reduce metabolic upkeep, but it still needs secretory/repair tissue, durable joints and turnover. No measured material presently supplies the entire assumed package of wet strength, low gas permeability, synthesis, repair and lifetime. Gas purity, oxygen/nitrogen ingress and delivery of biologically produced hydrogen remain open requirements.

## What the coupled reference cases actually pay

The runner combines structure, an explicit partition, hydrated host and community, free water, reserves, gas renewal, maintenance, consumers and tissue replacement. The following **diameters**, not radii, use the half-cycle-darkness stress case; that light regime is retained as a comparison rather than assigned to every location.

| Reference quantity at 10 km | 40 m diameter | 100 m diameter |
|---|---:|---:|
| Actual supported mass per projected m² | 18.70 kg | 24.59 kg |
| Required lifting-mixture volume fraction | 63.9% | 33.6% |
| Actual hydrogen inventory | 1.75 t | 14.41 t |
| Remaining assimilate after current obligations | 1.029 kg C/m²/year | 0.735 kg C/m²/year |
| Conditional parent-funded daughter allocation time | 7.04 years | 16.14 years |

The reference combines an assumed 5 MPa allowable, 100 µm barrier, low cellulose-like permeability, 1 kg living dry host/m², 90% tissue water, 5 kg free water/m², low maintenance, and imposed 200 W/m² cycle-mean illumination with 1% solar-to-hydrogen conversion. It is a requirements tuple, not an observed phenotype. Its separate hydrogen and carbon allocations cannot spend the same photons twice. Short-deficit and wind-matched continuous-light cases change timing while holding specified production inputs fixed; they do not derive latitude-dependent productivity.

Positive remaining assimilate has not yet paid future construction respiration, locomotion or every unresolved organ. For example, continuous spherical correction at 1 m/s consumes another 0.893 kg C/m²/year on the declared energy basis, enough to reverse some otherwise favorable cases. A hypothetical daughter allocation time excludes a demonstrated growing attachment, juvenile support, organs, mortality and a complete independent life cycle.

## Water changes both physiology and flight trim

The reference carbon-uptake calculation at 10 km, 400 ppm CO₂, 60% RH and leaves 5 K above air requires roughly **458 kg water/m²/year**, before additional cuticular/night losses. Cloud fraction is not collectable liquid or a colony's wet duty. Fog interception requires relative encounters; exact drift in uniform wind supplies no sweep-through. Condensing vapor requires latent-heat rejection. Rain must be admitted or shed within finite mass and drainage limits.

The new trim ledger nets uptake, evaporation, drainage, food and metabolic products against gas leakage and production. Opposing simultaneous flows can cancel. For the retained mixture, each kilogram of hydrogen provides about 13.40 kg of net support. If a prescribed daily cycle loses **1 kg water/m²**, holds altitude by venting hydrogen, then replaces water and gas, it vents and later replaces 0.0746 kg H₂/m². Repeating that strategy requires **103.6 W/m² of hydrogen chemical output** before conversion losses. This is a conditional **vent/refill strategy cost**, not a tax on every kilogram evaporated; matched water uptake, retained gas or other controls may avoid it.

Retaining gas by compression instead needs real pressure structure and actuation. An illustrative 1 kg/m² payload loss raises pressure by several kilopascals in the modeled arrangements, exceeding a 400 Pa starting allowance. Ideal pumping work can be modest while the structural requirement remains substantial. Variable altitude, thermal buoyancy, finite storage and synchronized intake therefore deserve evaluation, with their actual mass, heat and water boundaries. Static lift and annual carbon closure alone cannot establish controlled flight.

## Evidence checkpoint

The discriminating next work is an assembled wet barrier under representative tension and gas gradients, followed by a coupled water/gas/pressure cycle and a reachable weather-and-light trajectory. Reproduction must include growth to independent lift, inherited nutrients, repair, survival and recruitment sufficient to replace losses. Correlated storms and consumers share these budgets.

The climate programme remains paused; no new GCM or CFD campaign underlies this study. Component conservation tests and independent reviews establish implementation consistency. They do not certify material synthesis, a viable organism, ecological abundance or population capacity. Detailed assumptions, primary sources and access limits remain in the branch-owned component reports and source registers.

## Checks and publication

All **101 focused floater tests passed**, including exact product regeneration, named constants and all producer/input
hashes. Independent review covered the component physics, coupled allocation, water trim, twilight units and this
synthesis. The plot was regenerated from the final product and visually checked; its source is tracked and generated
images remain outside Git. The ecology [check record](https://github.com/dgoldman0/Terluna/blob/f1f8383bd1aa1ba60a301b5372a037c946537360/research/studies/floater_viability/checks.json)
records the full scope.

Required `make -k check` runs completed before committing: ecology **1,355 passed, 111 skipped, 30 failed**; main
**1,221 passed, 111 skipped, 30 failed**. The failing Python test identifiers match the preceding checks: 29 need absent
ignored external inputs or generated products, and one process-lookup test fails in this execution environment.
Layer checks, provenance/baseline checks and ensemble integrity passed on both branches. Calendar checks passed
13 with one skipped per branch because optical atlas/page assets were absent. Immersion bakes completed; the
Makefile's Node 24 directory-test invocation still fails and stops the subsequent JavaScript subcommands. Full
`make check` is therefore not green. No unrelated environment or model changes were made to suppress those failures.

The detailed work is published on `domain/ecology`, with this summary on main. Provisioning and resources retain
the previously published cross-domain evaluations; no new capacity, resource demand or operating route is adopted.
