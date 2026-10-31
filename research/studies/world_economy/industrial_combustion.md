# Industrial production in the oxygen–nitrogen atmosphere

The existing fire work supports a first industrial assessment. It supplies useful transport calculations and evidence that material behaviour can change substantially with gravity. It does not yet establish a production process, material acceptance envelope or financial benefit at **121.59 kPa, 17.5% oxygen and lunar gravity**. The extra nitrogen nevertheless has a credible, material-specific benefit: it can make sustained flame propagation harder while leaving oxygen partial pressure near Earth's. That benefit warrants investigation in polymer and composite fabrication; its size must be established separately from the effects of low gravity.

This assessment is for industrial planning within the three-world economy. It retains the 1.2 atm lock and uses 301 K for matched gas calculations. It evaluates production materials, process atmospheres and equipment consequences. The [calculation](industrial_combustion.py), [inputs](industrial_combustion_inputs.json), [results](results/industrial_combustion.json) and [source ledger](industrial_combustion_sources.json) distinguish conservation identities, experiments under other conditions and proposed applications. No fire-loss allowance or world-product total changes.

## What the earlier work establishes

| Existing work | Useful coverage | Boundary for industrial use |
|---|---|---|
| [Archived settlement assessment](../../../archive/july_2026_knowledge_bundle/terluna_moon_terraforming_bundle/07_HUMAN_SETTLEMENT_HEALTH_AND_CULTURE.md#7-fire-safety) | Recognizes altered buoyancy, smoke transport and the need for partial-gravity research | Conceptual precautions under an older 80 kPa / 24% O₂ example; no qualified production process |
| [Engineering fire models](../../../engineering/fire/README.md) | Earth correlations, explicit air properties, Froude scaling, plume/ceiling flow, compartment heat, detection proxies and ventilation | Numerical reproduction of correlations; radiation, soot, fuel pyrolysis, flame chemistry and spray suppression are not established at lunar gravity |
| [Central port study](../port_fire/README.md) | Applies those models to rooms, halls, shafts, water systems and evacuation; brackets design fires | Earlier tower geometry and roughly 0.74–0.95 atm port zones. These are neither the current surface production environment nor a factory fire study |
| [Ecology synthesis](../../../biosphere/ecology/README.md#climate-biomes-and-fire) | Distinguishes ignition frequency from uncertain natural-fuel propagation | No natural-fuel measurements at the combined 17.5% O₂, 1.2 atm and lunar gravity; no industrial biomass or dust qualification |

Three inherited interpretations need particular care. First, the port's peak heat-release multiplier of about **0.65** derives from a correlation between normalized mass loss and oxygen concentration at a fire's base. The underlying experiments involved ventilation and combustion-product mixing, supplemented by smaller fuel tests. The port also doubles time to reach 1 MW as a separate slower-growth assumption. Neither establishes a universal reduction for surface factories, and the Earth design fire is a comparison case rather than a proven worst case. The original equation and its conditions have now been inspected directly. [Peatross and Beyler, 1997](https://publications.iafss.org/publications/fss/5/403/view/fss_5-403.pdf).

Second, the port's **23.5% oxygen screening proposal** adds the largest observed shift in a small material sample to the design oxygen fraction. It can prioritize specimens for testing. It cannot guarantee that every assembly passing a terrestrial test will self-extinguish under lunar process conditions. Third, its **13 K ceiling-jet rise** represents a smoke-detection surrogate; it does not calculate particle transport into a particular optical detector. Its sprinkler response calculation establishes modeled activation, with subsequent fire control assumed. These distinctions are explicit in the model audit; the historical numerical product remains unchanged.

## Evidence for the nitrogen benefit

NASA's compiled material data extend to **119.3 kPa for some specimens**, close to the design's 121.59 kPa. Across the central pressure range, the authors found relatively weak pressure dependence of oxygen-percentage self-extinguishment thresholds for 17 of 22 tested materials. Maintaining oxygen partial pressure while lowering its percentage can therefore improve propagation behaviour. The result concerns particular upward-flame tests; the small remaining pressure extrapolation does not resolve the larger gravity, geometry or airflow differences. [Harper et al., 2016, methods, findings and Table 2](https://ntrs.nasa.gov/api/citations/20160001047/downloads/20160001047.pdf).

NASA's 2024 configuration handbook separately recognizes oxygen and diluent partial pressures. Its December 2024 revision also addresses assemblies that fail even though their constituent materials pass individually, and explicitly leaves partial-gravity effects largely unresolved. This is directly relevant to composite layups, bonded joints and cable assemblies manufactured on the Moon. [JSC 29353D, sections 5.1–5.2 and revision history](https://standards.nasa.gov/system/files/tmp/JSC%2029353D%20Flammability%20Configuration%20Analysis%20for%20Spacecraft%20Applications.pdf).

Gravity supplies a competing effect. Cast-PMMA experiments across 0.04–1 g found greatest flammability near lunar gravity, at an average 57.6 kPa and 13.2–15.2% oxygen. Those conditions differ from the design atmosphere, but they rule out treating partial gravity as a monotonic improvement between Earth and freefall. [Olson et al., 2024, publisher abstract and methods excerpt](https://doi.org/10.1016/j.firesaf.2024.104267).

Flow is part of that competition. ISS experiments on PMMA rods found that low forced flow could sustain propagation at oxygen fractions where the corresponding terrestrial flame failed; spread also changed nonmonotonically with flow. Removing buoyancy can remove cooling as well as oxygen transport. A factory's ventilation, moving web or extraction hood can consequently change material behaviour. [Link et al., 2018, methods and results](https://pmc.ncbi.nlm.nih.gov/articles/PMC5760683/). Ignition needs a separate assessment: subatmospheric PMMA experiments found a nonmonotonic ignition-delay response to pressure, including at constant oxygen partial pressure. [McAllister et al., 2010, abstract](https://research.fs.usda.gov/download/treesearch/35607.pdf).

The evidence supports **material-specific propagation resistance from the atmosphere**. Its net industrial effect at lunar gravity remains unmeasured.

## Quantities that can already be calculated

The comparison below uses rounded dry binary mixtures, 21% O₂ at 1 atm versus 17.5% at 1.2 atm, both at 301 K. Argon, humidity and trace gases are omitted. The [existing atmosphere properties](../../../atmosphere/radiative_convective/thermodynamics.py) supply molecular masses and heat capacities near 300 K.

| Quantity | Lunar value relative to the reference | Industrial meaning |
|---|---:|---|
| Oxygen partial pressure | 1.000; both 21.278 kPa | Equal oxygen inventory in an equal gas volume at equal temperature |
| Nitrogen per mole of oxygen | 1.253 | About 25.3% more nitrogen accompanies a specified oxygen feed |
| Gas density | 1.194 | Mass-flow and volume-flow settings cannot be interchanged |
| Heat capacity per gas volume near 301 K | 1.200 | More heat absorption for a small imposed temperature rise; this is not a flame-temperature calculation |
| Oxygen per kilogram of air | 0.837 | A fixed air-mass feed delivers about 16.3% less oxygen |
| Air volume for complete combustion of a fixed fuel amount, at local ambient pressure and equal temperature | 1.000 | A volumetric oxygen-supply duty does not gain the fixed-heat fan saving of the earlier cooling example |

Extra nitrogen distributes released heat over more gas and changes transport and reaction conditions. Equal oxygen partial pressure still leaves oxygen available for surface oxidation and dissolution. Consequently, lower oxygen percentage alone does not establish an oxidation-resistant annealing atmosphere, an oxygen-free polymer-curing environment or protection for reactive powders. Such processes need their own oxygen specification and enclosure.

Two counterfactuals must remain distinct. Comparing **1 atm / 21%** with **1.2 atm / 17.5%** isolates normoxic nitrogen addition. Comparing **1 atm / 17.5%** with **1.2 atm / 17.5%** adds oxygen as well as nitrogen; that is a different process question and need not improve flammability.

## Production opportunities worth pursuing

**Selected polymer and composite products are the strongest material opportunity.** Candidate lines include cable insulation, bonded panels, fibre composites, conveyor components, reusable moulds and temporary fabrication fixtures. If an assembly has adequate propagation resistance in the design air at its actual airflow and temperature, its manufacturer might use less flame-retardant loading or a thinner barrier, retain more ductility, tolerate a wider recycled-feedstock range, or reduce rejected output. Each is a separate hypothesis: lower additive loading can alter strength and processing, and a passing surface can still conceal a failing adhesive or core. The test must compare the finished assembly and its production configuration, including exposed edges and process residues.

This opportunity is most interesting when combined with existing low-gravity fabrication candidates. A large, lightly supported composite can benefit from low gravitational stress while its selected matrix or surface treatment benefits from additional diluent. The combined product must outperform an optimized terrestrial process or an enclosed/rotating array factory. A result qualified only for lunar air would serve domestic production first; Earth and array buyers need qualification in their own atmospheres. Qualification for the buyer's environment is part of the export cost.

**Large assembly and repair operations could gain flexibility if testing confirms a lower tendency for incidental flames to propagate.** Long bonded seams, extensive temporary tooling and partly completed insulation create exposures unlike finished products. Ambient dilution could reduce the amount of specially conditioned workspace needed around a particular qualified operation. The economic gain would be avoided enclosure, handling or fire-related interruption per accepted unit. Solvent removal, worker exposure, process cleanliness and quality can still require an enclosure even when propagation improves, so the benefit will differ sharply between operations.

**Processes needing modest additional nitrogen dilution may use less supplied gas.** This is a calculable opportunity once the required endpoint is known, but lower starting oxygen percentage does not guarantee lower nitrogen consumption. At higher pressure, each vessel-volume exchange contains more moles. For ideal continuous purging of a perfectly mixed vessel,

\[
n_{N_2}=\frac{PV}{RT}\ln\!\left(\frac{y_{O_2,0}}{y_{O_2,f}}\right).
\]

The following endpoints are mathematical examples, **not inerting limits or operating recommendations**. A material-specific endpoint must first be established at the actual pressure, temperature and process condition.

| Endpoint, 1 m³ vessel at 301 K | Reference nitrogen | Lunar nitrogen | Change in supplied mass |
|---|---:|---:|---:|
| Both end at 10% O₂ | 0.841 kg | 0.762 kg | −9.5% |
| Both end at 1% O₂ | 3.453 kg | 3.895 kg | +12.8% |
| Both end at 1 kPa oxygen partial pressure | 3.468 kg | 4.162 kg | +20.0% |

The first two rows deliver different final oxygen partial pressures across worlds. The last delivers equal oxygen partial pressure. Real purges also depend on leaks, desorption, dead volumes and whether flow displaces or mixes the original gas. The savings above concern gas quantity alone; separation, compression and recovery remain separate costs.

**Controlled combustion may benefit where peak temperature is the constraint.** Added nitrogen could moderate local heating in an appropriately designed burner or thermal-treatment process, potentially reducing hot spots or the external dilution needed. Lower peak temperature can reduce thermal NOx, but nitrogen availability, oxygen and hot-gas residence time also matter. There is no quantified emissions credit here. [EPA process-heater mechanism](https://archive.epa.gov/ttn/ozone/web/pdf/procheat.pdf). This candidate is weaker than the material opportunity: Earth and array equipment can deliberately dilute or recirculate gases, while an electrically heated process may avoid combustion entirely.

## Changes that can increase production costs

**High-temperature combustion pays for the extra nitrogen.** With complete stoichiometric methane combustion as an illustrative balance, nitrogen rises from 7.524 to 9.429 mol per mole of fuel. Wet exhaust rises from 10.524 to 12.429 mol, **18.1% more gas molecules** to heat or process. Fuel chemistry is unchanged. Furnace heat recovery can reduce the added sensible-heat loss; hot processes may favour oxygen enrichment or electrical heating. The established reason for industrial oxygen enrichment is precisely to avoid heating and exhausting unnecessary nitrogen. [DOE process-heating sourcebook, Tip Sheet 3](https://www.energy.gov/sites/default/files/2014/05/f15/process_heating_sourcebook2.pdf).

More exhaust moles do not mean 18.1% larger ducts. At the respective ambient pressures and the **same exhaust temperature**, the calculated wet-volume ratio is **0.984**, because higher pressure nearly offsets the added moles. This ratio is independent of the chosen common temperature under the ideal-gas assumption, provided water remains vapour. Actual temperatures, heat transfer, pressure losses and condensation need a furnace model.

**Controls and gas analyses need the correct composition.** In the same complete-combustion methane example with 10% excess oxygen, the dry exhaust contains **2.11% O₂ on Earth and 1.73% on the Moon**. Both have the same fractional oxygen excess. An Earth-calibrated conversion from measured oxygen percentage to excess air would give the wrong answer. Instrument-specific pressure compensation and fuel-flow calibration are separate requirements. Similarly, unchanged contaminant mass can appear as fewer parts per million when more nitrogen accompanies the products; emissions comparisons need mass per accepted output or energy service.

Dry stoichiometric CO₂ fraction falls from **11.73% to 9.59%**, while its partial pressure at the respective ambient pressures falls by only about 1.9%. Thus the separation duty cannot be inferred from concentration alone. Exhaust molar throughput, partial pressure, solvent or membrane behaviour, heat integration and available oxygen must be compared together. Nitrogen removal changes the entire capture process, as the [NETL oxy-combustion account](https://netl.doe.gov/carbon-management/combustion) explains.

**Powders require their own assessment.** Slow settling can improve suspension processes and also prolong combustible clouds or transport fines into filters. Primary dust data show strong effects of size, turbulence, temperature, oxygen and initial pressure; coal dust in the cited experiments still propagated well below 17.5% oxygen. A polymer's upward-flame result therefore supplies no dust-explosion limit. Nor should less visible flame be taken as freedom from hot deposits or incomplete oxidation. [Cashdollar, 2000, sections 4.4–4.7](https://stacks.cdc.gov/view/cdc/9139).

A useful thermodynamic check prevents another overinterpretation. In an ideal sealed constant-volume gas with unchanged mole count and constant heat-capacity ratio, adding heat gives \(\Delta P=(\gamma-1)Q/V\). Additional inert gas lowers temperature rise but does not, by itself, lower this pressure rise for the same released heat. Real reactive mixtures change mole count, heat capacity, losses and reaction completeness. This identity is not a blast prediction; it shows why thermal dilution alone cannot supply an explosion-pressure discount.

Finally, combustible vapours and polymer parts are only some industrial hazards. Smouldering deposits, self-heating bulk material, battery thermal runaway, reactive chemistry and oxidizers carried within a product require their own mechanisms and controls. The port's gas-phase buoyancy model does not assess them. Equipment designed for intentional flames also needs ignition, stability and burnout checks: a change that helps a finished polymer self-extinguish can hinder a burner or leave process contaminants unoxidized.

## A practical next comparison

The most useful first material trial is a **bonded panel or cable/insulation assembly** already proposed for lunar production. Keep material lot, thickness, support, exposed edges and heating history fixed. Compare 1 atm / 21% and 1.2 atm / 17.5% at the same gravity, then repeat the matched atmospheres at lunar gravity with still air and representative process flow. A 1 atm / 17.5% control distinguishes nitrogen addition at equal oxygen partial pressure from simply reducing oxygen availability. Vary process heat and airflow deliberately; short drop tests alone cannot establish long thermal soak, full assembly behaviour or production-scale consequences.

Record ignition response, propagation, sustained burning/extinction, heat release, smoke and residue alongside the actual product's strength, quality and accepted yield. For dust-bearing lines, assess the process dust and hot deposits separately. For a burner candidate, compare equal useful heat and product quality, including fuel, nitrogen/oxygen conditioning, fan work, heat recovery and exhaust treatment. The array comparison should use its chosen enclosed process atmosphere and optimized freefall or rotating equipment.

A positive result would justify a specific material substitution, reduced conditioning duty or measured yield improvement. It would then enter the existing economic accounts through accepted output, operating cost or asset losses, with the atmosphere's construction and maintenance already counted. Until those comparisons are available, the supported result is a ranked industrial opportunity and a set of calculable equipment changes, without an assumed general manufacturing or fire-safety premium.

Reproduce the reference calculation without replacing the committed product:

```sh
PYTHONDONTWRITEBYTECODE=1 python -m research.studies.world_economy.industrial_combustion
```

The default output is ignored `research/runs/world_economy/industrial_combustion.json`. Verification and source-reading limits are recorded in [checks.json](checks.json), stage 11.
