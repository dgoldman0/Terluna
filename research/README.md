# Terluna research

The hub for the project's shared research: what exists in each domain, the
cross-domain studies, the findings, the plan, the condition register and the
provenance of imported material. It supports the five-paper
[Open Moon ensemble](../ensemble/). **Constructing and Sustaining an Open Moon is
the next full manuscript**; the core's opening landscape does not define the
limits of the research.

Models live in their domain folders. This folder holds what spans domains:

| File | Holds |
|---|---|
| [findings.md](findings.md) | Equations, bounds, source limits and numerical examples from the environment screens |
| [plan.md](plan.md) | The core-first work order |
| [status.json](status.json) | Condition, open questions and next task for each topic |
| [decisions.md](decisions.md) | The author's decisions, where each is recorded, and whether it stands |
| [provenance.json](provenance.json), [archive_status.md](archive_status.md) | Every original archive member and its disposition; unresolved recovery |
| [check.py](check.py), [checks.json](checks.json) | Byte integrity of the pinned imports and reproduction of the historical baseline |
| [baselines/feasibility](baselines/feasibility/) | The September multi-domain feasibility script and its full report, kept intact |
| [studies/](studies/) | Cross-domain studies, each with its runner, results and write-up |
| [integration/](integration/README.md) | The joint integrations: of 8 October 2026, the five lines merged, the conflicts and repairs, the checks, and Codex's review that planned it; and [of 12 October 2026](integration/domains.md), the three domains merged with four refreshes of main's studies |

## Domains

| Topic | Current calculation | Principal open condition |
|---|---|---|
| [Atmosphere](../atmosphere/) | Solved conduction/advection/Jeans column, band ledger, lower air, column soundings, line-by-line radiative–convective column with shield-filtered sunlight, correlated-k thermal scheme, global-mean radiative–photochemical equilibrium (ozone, stratosphere, surface UV) per shield with solved balance temperatures, non-LTE CO2 cooling bounds, exobase with idealised-filter and titania-film heating and atomic-oxygen escape, line-by-line radiation benchmarks for a GCM, fluorinated trace-gas forcing from measured bands, the exobase with less CO2, the loss response to protection (screening); the thermal column's infrared cooling and traced heating; [electricity](../atmosphere/electricity/): charging, ionization and conductivity, Moon-wide lightning and a tower's exposure in the storms | 3-D middle-atmosphere transport and its effect on the exobase; light bypassing the shield film through aperture gaps |
| [Climate](../climate/) | Conservative periodic latitude–longitude thermal screen; canopy and forest-patch flow; clear-sky line-by-line budget from the atmosphere domain; single-column month-long day and night over land and sea; GCM boundary files, experiment plan and a restartable ExoPlaSim runner (design case near 299.6 K with PlaSim's cloud scheme corrected for lunar gravity; water and cloud sensitivities); [the seas' waves](../climate/waves/): two lunar cycles of SWAN waves over the nearside sea and Smythii–Marginis from the GCM's surface stress (mean Hs 1.39 and 0.88 m), coastal nests to 118 m terrain, and SWASH breaking and run-up | Climate programme paused on 2026-09-30 (decisions register): CM1 and the GCM disagree on warmth and humidity over land, so comfort is carried as a range; cloud amount, dynamics and terrain. Waves: further lunar cycles and the South Pole–Aitken sea, the monthly tide's water level at the coasts, a finer basin grid, sediment on regolith shores The atmospheric-electricity study reads its CM1 storms (Studies below) |
| [Biosphere](../biosphere/) | Carbon reserve theorem/model, bounded oxygen box, tree and forest-patch mechanics, clear-sky canopy photosynthesis by wavelength under the Moon's and Earth's skies, the whole plant's carbon through the lunar night with twilight, fruit designed for the lunar day on that plant (program, set time, size and share); [ecology](../biosphere/ecology/README.md)'s first screen of the ultraviolet by action spectrum, settling at lunar gravity, and methane and N2O without OH | Measured traits, complete life cycles and ecological interactions; whether leaves survive weeks of darkness and use two weeks of continuous light; whether a fruit's program can be compressed to one lunar day |
| [Illumination](../illumination/) | Spectral clear-sky solver, A1–A3 light-transport references, site geometry and earthlight, bright stars, surface spectra, a shielded spherical sky on the solved columns and regional evening cloud radiance with evolving views; the dated light calendar (2000–2500), Earthlight at the measured albedo, rough-water reflection and the water's colour, the ring fleet's light in the night ([fleet_light](../illumination/fleet_light/model.py)) and haze ([aerosol](../illumination/aerosol/haze.py)) | Regional profiles, three-dimensional cloud structure, aerosols, refraction, polarization and the earthlight spectrum |
| [Protection](../protection/) | Original component model and its full September design report, plus spectral coverage audit; the stored film's transmission from hard X-rays to the far ultraviolet; DE440 dynamics, the ring fleet's films and the finite-tile geometry behind the [solar shield study](studies/solar_shield_array/README.md) | Aperture leakage (gaps, pinholes, edges), clean operation, particle transport and lifetime resources |
| [Engineering](../engineering/) | Original transport/renewal accounts plus plume-heat sensitivity; tall lattice towers sized for gravity and wind, and wind devices screened for output, load and fliers; buoyant, winged and rotor flyers sized for any gravity and air; the shield's power, storage, electromagnetic, large-coil and heat-recovery accounts | Complete industrial network and safe source-to-use routes |
| [Geography](../geography/) | LOLA topography above the GRAIL geoid; hydrostatic water storage (level curves, basin joins); the atlas at the selected 28% water share, named from the IAU gazetteer; a first estimate of rivers and rain-fed lakes from the climate run's runoff; explicit cross-shore profiles for the wave studies; [every sea's monthly tide](../geography/README.md#the-monthly-tide) (3.7 m typical range over the nearside sea) with a dynamic check | How the water inventory divides between seas and rain-fed lakes (lakes hold about a third of the seas' volume in the first estimate); groundwater; Mare Fecunditatis's strait, which sets its fortnightly tides |
| [Resources](../resources/) | [Phosphorus return requirements and conservative finite reservoirs](../resources/phosphorus.md), with shared population material comparisons and no selected allocation; first screen of the inventories, the heat of delivery, the hydrogen left by shipping oxygen as water, oxygen taken by weathering rock, the shield's film, flooded ore, carbon beyond the biosphere and fresh seas, with an evaluation of the repository and the literature | The nitrogen's source under the conservation principles; a redox account of every import; where deliveries' energy goes; the shield's material life cycle and recovery fraction; carbon, oxygen and element budgets through time |
| [Provisioning](../provisioning/) | First screen of the heat ledger with the ring fleet's own infrared, land per person for food, the metropolis's night, tides and rivers, transit at lunar gravity and fusion fuel, with an evaluation of the repository and the literature | The fleet's infrared in the climate and shield design; the Moon's population and demand; the power system through the night and power delivered from orbit; the food system; freight between the surface and orbit; water, sanitation and work |
| [Habitation](../habitation/) | Design concepts; a first optical-comfort screen with surface luminance, eye illumination, shadow contrast and Earth controls; the [summit metropolis brief](../habitation/summit_metropolis/README.md), with the port, its fire screen, and the sky ships and flyers sized in the studies below | Practical inhabited capacity; visual comfort in resolved scenes and human validation |

Rendering that displays these results is in [visualization](../visualization/);
the explorable experience, which reads them as baked products, is in
[immersion](../immersion/). Constants and scenarios shared by every lane are in
[shared](../shared/).

## Studies

| Study | What it couples | Boundary |
|---|---|---|
| [environment_screens](studies/environment_screens/) | Thermal column, band ledger, climate screen, carbon and oxygen requirements, protection plume heat | Conditional consequences of selected equations; the spectrum/material/chemistry interface is unfilled |
| [megaforest_wind](studies/megaforest_wind/README.md) | Lower air (atmosphere), canopy flow (climate), tree mechanics (biosphere): 480 static load envelopes | Material and soil traits hypothetical; wind climate, gusts and evolution unmodelled |
| [conservation](studies/conservation/README.md) | Geography (the atlas), atmosphere (design pressure), climate (albedo, polar temperature), engineering (enclosures, sourcing) and a written heritage evaluation: the two preservation programmes, gates, register, scientific zones and community settings | Heritage evaluation is a first pass awaiting the author's review; numbers are hydrostatic geometry and screening estimates |
| [forest_patch](studies/forest_patch/README.md) | Three-dimensional patch airflow, tree load sharing, compliant foundations, shared soil and damage | Wind climatology, gusts, nonlinear failure and evolved morphology open |
| [atmospheric_co2](studies/atmospheric_co2/README.md) | Biology (plant CO2 needs, a C3 leaf model), weathering (a basalt law) and atmosphere (CO2 forcing, exobase): the CO2 level plants need, drawdown by biosphere build-up and weathering, and the return a settled Moon needs | Literature synthesis with Earth rates; its box model is carbon_cycle |
| [carbon_cycle](studies/carbon_cycle/README.md) | Resources (weathering, imports, kilns), ecology (plants, soils, sky, lakes, methane), provisioning (the food blend, people's stocks) and geography (waters, runoff): the carbon cycle as a box model from the seeding to year 10,000, the managed return that holds 400 ppm, the build's import and the oxygen burial leaves | Box model on Earth rates and the project's products; waters well mixed, design guesses labelled |
| [sea_chemistry](studies/sea_chemistry/README.md) | Resources (weathering, the seas and lakes, phosphorus), the carbon cycle's chemistry, geography (the waters) and conservation (erosion): the seas' make-up from the water's delivery through a billion years, fresh, brackish, salted or with phosphorus, left alone or held, and when sediment and their own deposits fill them | Box model on textbook constants and design-guess thresholds; the brines past an ionic strength of 0.5 by mass balance |
| [terrain_through_time](studies/terrain_through_time/README.md) | Geography (the drainage routing and relief), resources (the seas' deposits, dredging), conservation (erosion), the carbon cycle's chemistry and ecology's builders: how long the lakes, hollows, small craters and relief last without upkeep, what keeping them costs, and where travertine and dam builders keep pools going | Routing on the 1.9 km grid with Earth's erosion, trapping, creep, tufa and beaver-dam relations; design guesses labelled; features below the grid by statistics |
| [springs](studies/springs/README.md) | Geography (stages 1 and 2 of groundwater), the carbon cycle's chemistry and weathering, and the terrain study's travertine: how long groundwater stays underground, what it dissolves, how hard the springs are and where they build travertine, and the weathering and carbon the groundwater carries | Steady flow and a generalized Vogel transit model on the 15 km grid; carbonate chemistry with field weathering rates, design guesses labelled; no lunar groundwater measured |
| [protection_architecture](studies/protection_architecture/README.md) | Atmosphere (the loss response), climate (dimming), protection (hardware), engineering (supply) and the decisions register: the requirements for the optical shield and charged-particle protection, and what each loss budget asks of them | Screening loss model with the exosphere's losses; the loss budget, the ultraviolet cut-off and the protected radius are the author's to set; hardware, holding and supply are still to design |
| [solar_shield_array](studies/solar_shield_array/README.md) | Protection, engineering, illumination and habitation: holding, cycling and formation dynamics of the shield fleet with the full ephemeris. The [relative-orbit plan](studies/solar_shield_array/relative_orbits.md) (2026-10-06) traces the cycle's ~345 km return error to the depth stack's orbital-energy spread and the hour-13.8 clearance stop to the pattern folding through its orbital plane | No accepted design yet. The [integrated comparison](studies/solar_shield_array/integrated_comparison.md) (2026-10-06) judges candidates on the shield requirements, net-positive habitat electricity, resources, stability and safety, with every outflow traced by direction. A zoned held screen needs 32.8 TW, but as designed its exhaust exceeds the proposed outflow allowance; with the planned magnetosphere holding off its ions it would still need gridded ion thrusters and capture of most of its unionized gas, and its propellant fails S6. The ring fleet, kept with photon forces, narrows to radial-facing tiles near 20,000 km. Edge-on turns would collide with the next ring, and a filter centred on each tile with reflectivity trim holds the attitude. Nested rings need one sail loading, which puts annulus tiles near 13 g/m². The planes shift the whole pattern about 1,600 km each way over a year, which 23% more rings cover, and the edge strips need steering or more oversizing: about 27–33 million tiles and 45–53 Gt. A 589-tile patch with every tile propagated stays clear and covers every interior receiver ray for four orbits. Light annulus films stop the extreme and far ultraviolet but pass solar X-rays, whose effect is the atmosphere domain's to settle |
| [inhabited_volume](studies/inhabited_volume/README.md) | Habitation, provisioning, flight engineering and ecology: shared surface/aerial population comparisons, wet district mass, buoyancy, projected area, wind-holding power, water and whole-population food/heat | Requirement screens on committed air and service products; no complete architecture, ecological occupancy allowance, population ceiling or adopted allocation |
| [lunar_cycle_ecology](studies/lunar_cycle_ecology/README.md) | Biosphere (plant and fruit carbon), geography (seas, lakes, runoff), atmosphere (air column, flight) and illumination (darkness): a register of design ideas for plants, animals, the night food web and the nutrient and mineral cycles, with light calculations and a recommended selection | Ideas and screening arithmetic; no ecosystem, population or nutrient-cycle model; the selection is the author's to make |
| [optical_comfort](studies/optical_comfort/README.md) | Illumination, surfaces and human vision: baseline screens plus 630 gaze cases and 24 finite matte scenes using a consistent shielded spherical sky, with an Earth control and local cloud and water-slope inputs | Conditional clear molecular atmosphere and imposed scene geometry; cloud radiance, measured surface reflectance, rough-water reflection and human-factors evidence remain open |
| [cloud_twilight](studies/cloud_twilight/README.md) | Second-cycle regional cloud survey, spectral cloud radiance and colour, an 8–80 h fixed-observer history and sheltered viewing deep into evening | One sampled cycle, explicit cloud extrusion and particle-optics scenarios; three-dimensional regional structure, terrain and cloud Earthlight remain further inputs |
| [sea_appearance](studies/sea_appearance/README.md) | Climate (the waves), geography (the tide and the seas), illumination (sky, clouds, earthlight) and biosphere (living water): how the seas look by day, through the long evening, under earthlight and in the far side's twilit night | Plan adopted 2026-10-04; lighting calendar of six coasts and the Earth over the seas done, giving dawn and dusk their photometric definition; the sea's slopes through the month, the reflection model, the waters' colours, the Earth's observed light and the seas by regime at four coasts done; clear sky, uniform column, earthlight with sunlight's colour |
| [atmospheric_electricity](studies/atmospheric_electricity/README.md) | Climate (the CM1 storms' graupel, ice and supercooled water), atmosphere (ionization and conductivity) and protection (the particles it lets through): the plan of 2026-10-02, from the storms already run and a conductivity column to an electrified CM1 storm, the global circuit and transient luminous events | Stage 1 done (2026-10-03). Stage 2 written up (2026-10-08): the electrified box's two lunar days under the lunar rules, its windows (the breakdown cap, ground strikes, charging laws, leakage, the leader's crossing), the storms' lives, the thunder and the nitrogen; one site, a 6-km grid, Earth-calibrated discharge rules. Its open questions, the global circuit and the plasma above it among them, are gathered in the README's "Where the study stands" (2026-10-08); stage 3 is next in the plan |
| [open_moon_aerosol](studies/open_moon_aerosol/README.md) | Biosphere (the ecology register's landscapes and organisms), climate (the GCM's regions, CM1's weather near the ground), geography (seas and lakes) and atmosphere (ions and conductivity): the particles the seas, forests, plains and polar lands put into the air near the ground, the conductivity and the cloud nuclei they give | Measured Earth analogues changed by lunar factors, in clean, central and loaded cases; the ground only, not the storms' levels; the seas' composition, crust cover and fuels await the domains |
| [summit_tower](studies/summit_tower/README.md) | Geography (the summit, its highland rise and high lakes), climate (the design run's winds at the site and the cloud-resolving ring), illumination (sunlight for panels), the ecology register's fliers and the engineering tower models: how the wind bears on a sky tower's height, how an open frame stands up to it, and what wind devices in the frame give and cost | Concept sizing and screening arithmetic; the models' winds raised by stated factors; no structural design or terrain-resolving wind |
| [port_fire](studies/port_fire/README.md) | The summit tower study's port (zones, air, floors, programme and winds) and the engineering fire models: smoke, heat, detectors and sprinklers at 0.16 g in the port's spaces, design fires and a materials rule, getting out, and the tower's shafts, water, sealed zone and firefighting from the air | First-order fire-engineering arithmetic; Earth correlations carried by Froude scaling, untested at partial gravity beyond small samples; lunar design fires a bracket; crowd flows and aerial evacuation assumed |
| [sky_ships](studies/sky_ships/README.md) | Climate (the design run's air and winds over the whole Moon, and the cloud-resolving ring's winds, updrafts and storms), the summit tower study's port and the engineering flight models: how large a buoyant or winged ship the air can carry at 0.16 g, the rule of similarity that sets it, and what grows with a ship | First-order sizing calibrated against LZ 129 Hindenburg and large transports; loads from Earth's airship and aircraft practice carried to lunar gravity; no ship designed |
| [sky_fleet](studies/sky_fleet/README.md) | The sky-ship study's product, the summit tower study's port and sized form, climate (the design run's air, the cloud-resolving ring's daytime mixed layer and storms) and the engineering flight models: the Open Moon's flyers from canopies, micro gliders and pedalled wings to sky boats, ferries, regional and long-haul ships, freighters and high platforms; how many the summit metropolis and its port need; how they share the air by height; the long-haul class for the crown, drawn to scale in visualization/sky-fleet | First-order sizing; reference flyers, people's power and cities' travel from published sources; the metropolis's travel a planning case; no flyer designed; the crown's berths a first arrangement |
| [infrastructure_review](studies/infrastructure_review/README.md) | The four studies above and the metropolis brief, read against main's corrected climate, waves, tides and optical studies, the atmospheric electricity branch's storms, the sea-appearance branch's twilight and the solar-shield branch's work off the Moon: what holds, what changes, storm protection, the night, the coasts, infrastructure off the Moon and the open questions; the branch's closing work | A review of committed results on each branch, the studies rerun on a temporary copy with main's corrected inputs, the corrected run's winds exported, and light arithmetic |
| [joint_synthesis](studies/joint_synthesis/README.md) | Every line of work after the joint integration: each line's assumptions about the others against their final answers, and joint calculations of the fleet's light in the night, the air under the shield with its lightning, lightning over the Moon and its nitrogen, the port's crown in the storms, haze, the ledgers with their cycle times, the regional magnets on the map and seven places through a lunar day | Built from the merged products, each calculation with its stated assumptions; the decisions it lists are the author's |
| [aerial_ecology](studies/aerial_ecology/README.md) | Ecology, resources and provisioning: aeroplankton and large photosynthetic aerial organisms together, their consumers, food harvest and a functioning phosphorus cycle | A requirements evaluation; no safe yield or species count assumed |
| [uv_biology](studies/uv_biology/README.md) | Protection, atmosphere and biology: restricted ultraviolet against vitamin D, plants and animals, with local ultraviolet provision | A fixed-column comparison; no added ultraviolet transmission adopted |
| [aerophytes](studies/aerophytes/README.md) | Ecology, the design climate's winds and storms, the air's gases and the engineering of membranes and structure: the giant aerophytes' size limits, ages, water, storms, sailing and photoperiod, and the sky reefs' living sails | Screen values and conditional requirements; no organism size or route adopted |
| [sky_ecology](studies/sky_ecology/README.md) | The aerophyte product, the design winds and resources' phosphorus: the sky's food web, tenants, guards and symbioses, and phosphorus aloft | Accounting on design guesses; no community or phosphorus cycle simulated |
| [ways_of_living](studies/ways_of_living/README.md) | Provisioning, habitation, climate, illumination and geography: eighteen ways of living on and around the Moon, each placed on the atlas with its physical situation and its costs per billion residents | Costs by setting; the distribution is the author's choice |
| [daily_life](studies/daily_life/README.md) | The ways of living through the Moon's 708.7-hour cycle: homes, movement, work, gathering, care and cultural life, from the products | The physical situation computed; how people organize it is open for the human paper |
| [food_systems](studies/food_systems/README.md) | Provisioning, ecology's diets and yields, resources' phosphorus, climate, geography and the heat ledger: eight food systems per billion residents, their fit to the light cycle, seven mixes for 4–8 billion, and the author's blend with heritage meat | Screens on published measurements; the blend's fields not yet placed band by band |
| [land_cities](studies/land_cities/README.md) | Geography (the atlas and the high ground), climate (rain and the water table) and engineering: the lava tubes and the highland cities for the ground's people | Mapped pits and rilles; the tubes' interiors unsurveyed |
| [distributions](studies/distributions/README.md) | Every study of people: four distributions of 25–30 billion across Earth, the Moon, the array and the colonies, the author's C′, the wild share, traffic, phosphorus and synergies | Placement corners on design guesses; C′ placed with the power-heavy mix |
| [living_light](studies/living_light/README.md) | Ecology's glowing organisms and night stores, illumination's twilight and Earthlight, and people at night: what living light gives, what people need outdoors at night, and the sky kept dark | Screen values; no glowing organism built |
| [domain_synthesis](studies/domain_synthesis/README.md) | Resources, ecology and provisioning after the joint integration of 12 October: heat, land, methane and N2O, phosphorus and the sky's layers at the author's choices | Arithmetic on screen products; the gases' warming in the deep column comes next |
| [light_limitation](studies/light_limitation/README.md) | The canopy's light response, the food study's land by band, the GCM's soil water, the waters' and the sky's production: where and when light limits growth, and what the planned dimming costs the biosphere through light and, from Earth's models, through the rain | A screen on committed products and Earth's waters; the Moon's rain response is the GCM's, in gases_glow_dimming |
| [gases_glow_dimming](studies/gases_glow_dimming/README.md) | Climate's GCM with the gases, the fleet's glow and the dimming, the atmosphere domain's correlated-k radiation, the light study's growth: whether the dimming holds the design climate, what the rain does, and what the glow's infrared needs | Three 20-year runs of one model at T21 and a clear-sky spectral check; the held point interpolated |
| [jura_island](studies/jura_island/README.md) | Geography (the atlas and LOLA's shores), the CRM's 3-D box on the corrected design case and the groundwater's shore reach: Jura Island's layout, its climate and its shore belts | One box of two lunar days on placeholder land; the rerun waits on the GCM's longwave and the fresh shores |
| [climate_targets](studies/climate_targets/README.md) | The GCM's settled warm Moon and design case, the gases, the fleet's glow, the beamed power and the water top-up, light and rain: what each target climate takes and the heat beyond the sunlight, after the joint integration of 17 October | Arithmetic on two runs of different physics, interpolated by the ground's temperature; no target is run |
| [lake_oxygen](studies/lake_oxygen/README.md) | Climate (the GCM's forcing over three crater lakes and the waters' surface), geography (the rain-fed lakes and their depths), resources (the held fresh seas) and ecology (oxygen demand by trophic state): the oxygen the seas and lakes hold at saturation in each target climate, and a lake model through the lunar month for stratification, mixing and dissolved oxygen in deep crater lakes and productive shallows | A one-dimensional lake model on one GCM's 3-day means at T21, corrected after the review of 18 October, the stirring below the night's mixing set by Earth's measured share of the wind's energy with its range and the mixing's form (Osborn's rule; Shih and Bouffard's energetic regime in mixing_forms) as brackets, and the resort interpolated between two physics; oxygen demand as bounded scenarios from Earth's lakes |
| [lake_anoxia](studies/lake_anoxia/README.md) | The waters' oxygen study's columns, the domain synthesis's methane ledger, the deep column's gas forcing, the sea chemistry and the land's phosphorus: the lakes' methane and what it does to the air and the dimming, the bubble point, phosphorus from anoxic floors and the remedies, sized | Arithmetic on the lake model's columns with Earth's rates, at the wind's share and its bracket by Osborn's rule, the energetic regime's on the 0–15° lakes; each overturn's release bounded by two rules |
| [lake_kinds](studies/lake_kinds/README.md) | Ecology's lake model at the resort's 301 K, geography's rain-fed lakes on LOLA, groundwater stage 2, the sea chemistry's soda path, the springs and the anoxic depths' remedies: every lake sorted by kind with what each kind needs | Lake-model columns at 301 K on a first-estimate lake set, Osborn's form, with the energetic regime's 352 m lake and giant beside them; the ventilated depth interpolated between the 100 m and 352 m rungs |
| [sea_oxygen](studies/sea_oxygen/README.md) | Geography (LOLA's seas and the atlas), climate (the GCM's weather at the seas' cells moved to 301 K by the lake kinds' construction), resources (the held fresh seas) and ecology (lake_oxygen's model and trophic states): the seas' basins, sills and hypsography, and the deep water's mixing, renewal and oxygen at each basin's typical and coldest cells, with the shelves' and bays' nights | One-dimensional columns on one GCM's 3-day means at T21 bounding a basin's circulation; Osborn's rule throughout, Shih and Bouffard's form at the two largest basins |
| [lake_methane](studies/lake_methane/README.md) | The corrected lake model's saved columns (lake_anoxia's targets, the lake kinds' 301 K ladder and giants, mixing_forms, the seas' basins), lake_anoxia's air and the sea chemistry's acceptors: the lakes' and seas' methane at 301 K by depth and kind, the air it holds, and what aerobic, overturn and anaerobic methanotrophy can cut | An offline methane tracer on saved columns, reproducing lake_anoxia's methane within 4%; Earth's methanotroph rates as design variables; Osborn's form with Shih and Bouffard's beside it |

The forest-patch folder also holds a [reviewed checkpoint](studies/forest_patch/reviewed/README.md)
from a separate lineage (57 tests, 90-m spacing, 20-m/s reservoir flow), whose
source is kept at commit `363b161` and on `checkpoint/forest-patch-review-57-tests`.
The two lineages differ in their assumptions. In both, the dome was an imposed
candidate: neither simulates an evolved canopy shape or establishes an optimum.

Run everything from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 python -m pytest
OPENBLAS_NUM_THREADS=1 python -m research.studies.environment_screens.run   # results in studies/environment_screens/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.megaforest_wind.run       # results in research/runs/megaforest_wind
OPENBLAS_NUM_THREADS=1 python -m research.studies.forest_patch.run --quick  # results in research/runs/forest_patch
python -m research.studies.atmospheric_co2.run                             # results in studies/atmospheric_co2/results
python -m research.studies.carbon_cycle.run                                # results in studies/carbon_cycle/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.protection_architecture.run # results in studies/protection_architecture/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.lunar_cycle_ecology.run     # results in studies/lunar_cycle_ecology/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.optical_comfort.run
OPENBLAS_NUM_THREADS=1 python -m research.studies.optical_comfort.directional  # results in studies/optical_comfort/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.cloud_twilight.run      # results in studies/cloud_twilight/results
python -m research.studies.open_moon_aerosol.run                           # results in studies/open_moon_aerosol/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.summit_tower.run            # results in studies/summit_tower/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.port_fire.run               # results in studies/port_fire/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.sky_ships.run               # results in studies/sky_ships/results
OPENBLAS_NUM_THREADS=1 python -m research.studies.sky_fleet.run               # results in studies/sky_fleet/results
python -m research.studies.joint_synthesis.night_sky                          # and the synthesis's other steps; its README lists them
python -m research.studies.aerophytes.run                                    # results in studies/aerophytes/results
python -m research.studies.sky_ecology.run --output research/studies/sky_ecology/results/sky_ecology.json
python -m research.studies.ways_of_living.run                                # results in studies/ways_of_living/results
python -m research.studies.daily_life.run                                    # results in studies/daily_life/results
python -m research.studies.food_systems.run                                  # then food_systems.blend for the blend
python -m research.studies.land_cities.run                                   # results in studies/land_cities/results
python -m research.studies.distributions.run                                 # results in studies/distributions/results
python -m research.studies.living_light.run --output research/studies/living_light/results/living_light.json
python -m research.studies.domain_synthesis.ledgers                        # results in studies/domain_synthesis/results
python -m research.studies.climate_targets.run                             # after its extract, with the climate venv
python -m research.studies.lake_oxygen.run                                 # after its lakes and its extract (climate venv)
python -m research.studies.lake_oxygen.mixing_forms                        # after lake_oxygen.run; the energetic regime
python -m research.studies.lake_anoxia.profiles                            # then lake_anoxia.run
python -m research.studies.lake_kinds.ladder                               # then lake_kinds.run; about 2.3 h on two cores
python -m research.studies.sea_oxygen.run --columns                        # after sea_oxygen.seas and its extract
python -m research.studies.lake_methane.run --work                         # after the lake and sea columns; then run
```

For the environment screens, supply `--protection-archive /path/to/Lunar_Protection_Model.zip`
to inspect the original optical inputs directly. Source-access and validation
status are in [sources.json](studies/environment_screens/sources.json) and
[checks.json](studies/environment_screens/checks.json).

## Historical baseline reproduction

```sh
python -m pip install -r research/requirements.txt
python research/check.py
# Restore the exact five protection inputs from the original project ZIP:
python protection/fetch_inputs.py --archive /path/to/Lunar_Protection_Model.zip
python research/check.py --require-inputs
```

Upstream retrieval is also available through `python protection/fetch_inputs.py --download`;
the exact expected hash and size must match. With missing protection data the
checker reports `BLOCKED` explicitly; it never substitutes synthetic inputs.

The [checked snapshot](checks.json) records an actual run using the original
archived inputs. Future checks write to ignored `research/runs/`. Numerical
reproduction and file integrity do not grant physical validation, full-source
admission, or manuscript clearance.

The September feasibility implementation is a historical multi-domain script;
its single intact copy lives in [baselines/feasibility](baselines/feasibility/),
and topic folders hold its original compact reference tables. Protection has its
own intact implementation. Model code and CSV tables are verbatim; the two JSON
reference files select original fields with unchanged values. New calculations
should state assumptions, track conservation and residuals, test convergence and
say which conclusion their outputs can change.


The 9 October [aerial-biosphere, UV and population follow-up](integration/aerial_biosphere_population.md) connects
new evaluations on the ecology, provisioning and resources branches. They share population comparison cases,
include roaming aerial habitation, and distinguish adopted directions from still-conditional capacities.

The [aerophyte viability follow-up](integration/floater_viability.md) summarizes the ecology branch's coupled
structure, membrane, carbon, water/trim and life-cycle evaluation, including persistent twilight and wind-assisted
Sun-following. These are conditional requirements; no sustainable organism size or new population capacity is adopted.

The ecology branch's [aerial-food-web](studies/aerial_ecology/README.md) and
[restricted-UV](studies/uv_biology/README.md) follow-ups evaluate the 9 October cross-domain discussion.
They retain explicit evidence limits and use [shared population comparison cases](../shared/scenarios/population.json)
with the provisioning and resources branches; no allocation or additional UV transmission is selected.

The [aerophyte viability study](studies/aerophytes/README.md) deepens the aerial work with coupled mass,
pressure, permeability, carbon, water/trim, reproduction and navigation accounts. Persistent twilight and following
the slow Sun are explicit regimes. Positive budgets remain conditional requirements, with no selected organism size.
Its [living-sail note](studies/aerophytes/sail_biology.md) gives a sky reef's sails distributed sensing and
movement, staged growth, a renewable suspension, other life histories and other uses, and computes signal travel,
fluid head, local torque and renewal; one sail tested for stability is next.

The [sky-ecology study](studies/sky_ecology/README.md) gives a sky of aerophytes its food web, tenant habitat and
phosphorus as screen values, with sudden loads timed against a 25-minute buoyancy period. Its [guards and symbioses
analysis](studies/sky_ecology/symbioses.md) compares partners that share a reef's defence and computes the live mass,
response distance and shared diet each needs; whether a mixture defends and persists is open.

