# Hydrogen containment and inflation work

Biologically derived material can contain hydrogen: manufactured cellulose sheets and historical animal-membrane gas cells establish that much. They do not establish a living, growing, weather-exposed envelope with durable seams, controlled gas composition and a positive whole-organism energy balance. This calculation separates measured transport from those unresolved requirements. The coupled study must also pass the mechanical, water, nutrient, night-storage and ecological screens.

`envelope.py` provides pure scalar calculations; `envelope_sources.json` records the source, access level, conditions and missing evidence. No measurement here is relabelled as a validated lunar membrane. In particular, a gas-barrier film can be extracellular and metabolically inactive while supported or renewed by living tissue; it need not be a single exposed layer of living cells. That architecture remains a proposal requiring measurements.

## Measurements and evidence limits

| Comparator | Actual H₂ result | Conditions and what was accessed |
|---|---|---|
| Regenerated cellulose, Kurata et al. 2022 | 5.78 × 10⁻¹⁷ mol m/(m² s Pa), about 0.173 Barrer | Room temperature; publisher abstract read. Test RH and thickness not recovered. |
| Polypropylene, same experiment | 2.24 × 10⁻¹⁵ SI, about 6.69 Barrer | Engineered comparator; identical access limitation. |
| Cellophane, NBS RP750, 1934 | 0.1 L H₂/m²/day; approximately 25.4 µm, 51 g/m², 16% glycerol | 25°C, hydrogen opposite air at 30 mm water gauge. Full official PDF and table checked. RH of the gas test and the reporting-liter reference state remain unverified. |
| Goldbeater's skin composite, NACA TM172, 1922 | Finished fabric 130–150 g/m²; only a qualitative few L H₂/m²/day | Official full report; processed intestinal membranes, cotton, adhesive and varnish. Protected gas cells are not an exposed living skin. No precise permeability inferred. |
| Sulfonated cellulose nanofibers, Bayer et al. 2021 | 30 µm membrane crossover 0.03–0.3 mA/cm² | Full paper: humid fuel-cell configuration, 80°C, 95% RH; H₂/N₂ crossover measurement. This is a measured device flux, not a room-temperature material coefficient. |
| MFC/LMFC, 2025 | 10.7, 12.1 and 11.9 Barrer | Primary publisher indexed excerpt only. Room temperature stated; H₂-test humidity and full methods unverified. Conditional sensitivity. |
| Chitin nanofibers, Wu thesis 2014, Table 2.3 | H₂ 0.024 Barrer; separately O₂ 0.006 and N₂ 0.0034 Barrer | Official repository table excerpt only; temperature, RH and thickness not verified. Conditional sensitivity, not wet-flight evidence. |

Stable primary links: [Kurata](https://doi.org/10.2472/jsms.71.903), [NBS RP750](https://nvlpubs.nist.gov/nistpubs/jres/13/jresv13n6p879_A2b.pdf), [NACA TM172](https://ntrs.nasa.gov/api/citations/19930083385/downloads/19930083385.pdf), [Bayer](https://doi.org/10.1007/s10570-020-03593-w), [MFC/LMFC](https://doi.org/10.1021/acssuschemeng.5c05172), [Wu thesis](https://repository.gatech.edu/bitstreams/485965a3-039e-4ad2-a174-4c2bdb76b99e/download). The source register includes the related chitin journal article and access failures.

The NBS 65% RH condition concerns weighing, so it is not used to label its H₂ test humid. Different cellulose preparations differ substantially; their range is not a measured universal wet/dry multiplier. [Wu and Yuan 2002](https://doi.org/10.1016/S0376-7388(02)00037-6) did test dry and water-swollen cellulose with H₂, but the accessible abstract did not supply a numerical wet coefficient. Water-vapor or O₂ barrier performance is never substituted for H₂ permeability. Wet-strength and natural-cuticle measurements belong to [the mechanics assessment](mechanics.md); neither supplies the missing transport coefficient for an assembled living membrane.

## Model and integration contract

For each gas species, intact-film flux is

\[
j_i=P_i\frac{x_{i,\mathrm{in}}p_{\mathrm{in}}-x_{i,\mathrm{out}}p_{\mathrm{out}}}{t}.
\]

`permeation_flux` returns signed outward mol/(m² s). One Barrer is represented by the conventional rounded conversion 3.35 × 10⁻¹⁶ mol m/(m² s Pa). The H₂ gradient can approach one bar while the mechanical gauge pressure is only 100 Pa. Using the latter would underpredict loss by roughly a thousandfold. `layered_permeance` adds ideal serial resistances, Σ(t/P); pinholes, seams and edge bypass are additional paths, not resistance in series.

`hydrogen_budget` accepts actual lifting **mixture** volume, gas-contacting barrier area and projected collecting area independently. It calculates H₂ inventory as xPV/(RT), then flux × barrier area. The coupled geometry may bound gas-contacting area between a sphere containing the trimmed lifting mixture and the full exterior area. Actual cap shape, gas partitions and their masses require explicit treatment. A gas fraction does not automatically imply the same area fraction. If a structural film also provides the barrier, count its mass once; a separate coating or partition adds its own mass.

Returned `initial_turnover_days` means initial inventory divided by initial loss. Sustaining that loss rate assumes replenishment maintains gas composition and pressure. It is not an emptying time for an unmaintained flexible balloon. The ideal calculation excludes seams, pores, punctures, material weathering and purification losses; these omissions can only make intact-film containment look easier.

Water splitting has a reference reversible work requirement of **237.141 kJ/mol H₂**, with **285.830 kJ/mol** total enthalpy at 298.15 K, one bar ([NIST-JANAF](https://janaf.nist.gov/tables/H-063.html)). The difference is reversible heat input. Neither should be silently identified with the approximately 242 kJ/mol lower heating value. The model's `gibbs_efficiency` is reference Gibbs work divided by supplied useful work: input = nΔG/η, once. It does not include a second photosynthetic conversion automatically. Initial inflation, maintenance work, water feed and oxygen coproduct are returned separately. Stoichiometry is one mole water and half a mole O₂ per mole H₂; the coproduct is vented or otherwise handled separately, not mixed into lifting gas.

## Thickness and energy sensitivity

The reproducible `reference_thickness_scan()` uses a full sphere, 100 kPa total pressure, 298.15 K, 98% H₂, η = 0.5, and fixed measured/reference permeability while thickness changes. This is a conditional comparison at common pressure, not a correction of every material to common temperature or hydration. Values below are **input W per m² projected area**; minimum Gibbs work is half as large. Sphere area/projected area is four.

| Comparator | 1 µm | 25 µm | 100 µm | 1,000 µm |
|---|---:|---:|---:|---:|
| Kurata cellophane | 10.75 | 0.430 | 0.107 | 0.0107 |
| Polypropylene | 416.5 | 16.66 | 4.165 | 0.416 |
| Chitin, incomplete conditions | 1.495 | 0.0598 | 0.0149 | 0.00149 |
| MFC 10.7 Barrer, incomplete conditions | 666.4 | 26.66 | 6.664 | 0.666 |
| LMFC 12.1 Barrer, incomplete conditions | 753.6 | 30.14 | 7.536 | 0.754 |
| Hypothetical 100 Barrer stress case | 6,228 | 249.1 | 62.28 | 6.228 |
| Film mass at hypothetical 1,500 kg/m³, kg/m² skin | 0.0015 | 0.0375 | 0.150 | 1.50 |

The scan also contains 10 µm cases. Thickening by 1,000 reduces intact-film leakage by 1,000 but adds 1,000 times the film mass at fixed density. Fabrication defects and thickness-dependent morphology can break that inverse relation. At fixed P and pressure, replacement power per projected area does not improve merely by increasing sphere radius; inventory per projected area grows in proportion to radius.

The prior aerial-ecology NPP scenario, 1,000 g C/m²/year with dry carbon fraction 0.45 and dry energy 18 MJ/kg, corresponds to **1.2675 W/m² projected** of net biomass energy. This is a comparison ceiling, not spare power: maintenance, reproduction, food harvest, locomotion, nutrient capture and repair compete for it. With a 100 µm barrier, several measured-material sensitivities already exceed that entire energy flow. The low-permeability cases do not fail this isolated leakage screen, but require compatible hydration, mechanical construction and organism-level evidence.

For the humid 30 µm S-CNF device, Faraday conversion of crossover current gives **1.55 × 10⁻⁶–1.55 × 10⁻⁵ mol H₂/(m² s)**. Applied over a spherical skin, that is **1.475–14.747 W/m² projected** minimum Gibbs work, or **2.949–29.494 W/m²** at η = 0.5. It is an 80°C humid comparator, not predicted lunar loss. The same study's 8 µm spray membrane reached 449 mA/cm² with cracks/pinholes; reducing thickness did not preserve ideal barrier behavior.

For historical liter-based transmission, retain the reported units. If a reporting liter were STP at 273.15 K/101325 Pa, each 1 L/(m² day) would imply approximately **0.98 W/m² projected** replacement input at η = 0.5. This conversion is a declared convention, not an inferred test condition or a numerical interpretation of the word few.

## Initial inflation remains substantial

For a fully gas-filled **100 m radius** sphere at the reference pressure, temperature and purity, H₂ inventory is **333.84 tonnes**, water consumed during inflation is **2,983.23 tonnes**, and minimum inflation work is **39.269 TJ**. At η = 0.5, input is **78.538 TJ**, or **2.500 GJ/m² projected**. Filling within one year requires **79.22 W/m² projected**, before leakage and construction. A one-bar atmosphere gives large lifting-gas inventory even at small mechanical overpressure.

That input equals approximately **62.5 years** of the prior NPP energy flow over the *final* projected area if all net biomass energy could be redirected to H₂. This is an energy-equivalent duration, not a growth trajectory: a smaller juvenile has a different collecting area, support load and gas requirement. External energy, staged inflation, low-density architectures, direct photochemical H₂ production and vegetative growth are separate alternatives requiring explicit budgets. Inflation energy and gas stock scale with actual filled mixture volume; root's trim model supplies that volume. Low leakage does not remove the inflation investment.

## Composition and open measurements

`fixed_volume_exchange` analytically integrates each measured species permeability in an isothermal **rigid** volume without production, reaction or purge. It permits inward O₂ and N₂ diffusion even while total internal pressure exceeds ambient. It approaches each external partial pressure independently. It is tested as a diagnostic, not used to claim a flexible envelope's flight duration or a flammability threshold. Numerical O₂/N₂ evolution for a proposed biological envelope is deferred until matched, condition-resolved measurements exist; the incomplete chitin table is not sufficient to establish wet flight behavior.

The decisive experiment is an assembled barrier tested for H₂/O₂/N₂ transport and water balance while under representative wet tension, temperature cycling, irradiation and biological exposure. Coupon strength, a water barrier and isolated low H₂ permeability do not together prove a viable assembled organism. The study needs growth-compatible seams, puncture repair, aging rates, gas production/purification costs and environmental survival. No actual closed ecology or autonomous giant biological floater is established by these calculations.

Run the focused checks from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q research/studies/floater_viability/test_envelope.py
```

The checks cover partial-pressure rather than gauge-pressure diffusion, film resistance, energy accounting, wet-device current conversion, independent gas volume/area, stoichiometry and exact limiting cases of species exchange. They establish arithmetic and model behavior; they do not validate the proposed material.
