# Holding hydrogen

**Biologically derived films hold hydrogen well enough for a floater: regenerated cellulose permeates 0.173 Barrer,
and goldbeater's skin lined the gas cells of historical airships.** At 100 µm the cellulose film loses hydrogen at
0.107 W/m² of replacement input per square metre of floater, a small share of the 1.27 W/m² of net production in the
aerial-ecology screen. The driving pressure is hydrogen's partial-pressure difference, close to one atmosphere, so a
film's loss is set by its permeability and thickness, whatever the small mechanical overpressure. Filling the gas is
the large cost: a full 100-m-radius sphere holds 334 t of hydrogen, 62.5 years of that net production over its area
([growth](growth.md) gives the ages that follow).

[envelope.py](envelope.py) holds the scalar calculations; [envelope_sources.json](envelope_sources.json) records each
measurement's conditions and how far it was read.

## Measurements

| Comparator | Hydrogen result | Conditions and what was read |
|---|---|---|
| Regenerated cellulose (Kurata et al. 2022) | 5.78 × 10⁻¹⁷ mol m/(m² s Pa), about 0.173 Barrer | Room temperature; publisher abstract; test humidity and thickness not recovered |
| Polypropylene, same experiment | 2.24 × 10⁻¹⁵ SI, about 6.69 Barrer | Engineered comparator; same access |
| Cellophane (NBS RP750, Sager 1934) | 0.1 L/m²/day; about 25.4 µm, 51 g/m², 16% glycerol | 25 °C, hydrogen against air at 30 mm water gauge; full official PDF and table read; the litre's reference state and the gas-test humidity unstated |
| Goldbeater's skin fabric (Chollet 1922, NACA TM 172) | Finished fabric 130–150 g/m²; a few litres of H₂ per m² a day | Full report; processed intestinal membranes with cotton, glue and varnish, used as protected gas cells inside hulls |
| Sulfonated cellulose nanofibre membrane (Bayer et al. 2021) | 30 µm membrane, crossover current 0.03–0.3 mA/cm² | Full paper; humid fuel-cell configuration at 80 °C and 95% RH, H₂/N₂ crossover |
| Microfibrillated cellulose, lignin-rich MFC (2025) | 10.7, 12.1 and 11.9 Barrer | Publisher excerpt; room temperature; humidity and methods unverified |
| Chitin nanofibres (Wu thesis 2014, table 2.3) | H₂ 0.024 Barrer; O₂ 0.006 and N₂ 0.0034 Barrer | Repository table excerpt; temperature, humidity and thickness unverified |

Links: [Kurata](https://doi.org/10.2472/jsms.71.903), [NBS RP750](https://nvlpubs.nist.gov/nistpubs/jres/13/jresv13n6p879_A2b.pdf),
[NACA TM172](https://ntrs.nasa.gov/api/citations/19930083385/downloads/19930083385.pdf),
[Bayer](https://doi.org/10.1007/s10570-020-03593-w), [MFC/LMFC](https://doi.org/10.1021/acssuschemeng.5c05172),
[Wu thesis](https://repository.gatech.edu/bitstreams/485965a3-039e-4ad2-a174-4c2bdb76b99e/download). Wu & Yuan (2002)
tested dry and water-swollen cellulose with hydrogen; the accessible abstract gave no wet coefficient, so none is
used. Cellulose preparations differ widely, and the material's name does not fix its permeability: the MFC films
permeate about 60 times more than the regenerated cellulose.

## The model

**Each species permeates as j_i = P_i (x_in p_in − x_out p_out)/t.** `permeation_flux` returns the signed outward flux
in mol/(m² s), with one Barrer taken as 3.35 × 10⁻¹⁶ mol m/(m² s Pa). The hydrogen gradient approaches one bar while
the gauge pressure is about 100 Pa, so using the gauge pressure would understate the loss about a thousandfold.
`layered_permeance` adds layer resistances Σ(t/P); pinholes, seams and edges are parallel paths.

**`hydrogen_budget` takes the lifting mixture's volume, the gas-contact area and the projected area separately.** It
counts the inventory as xPV/(RT) and the loss as flux times contact area; `initial_turnover_days` is the inventory
over the first-day loss, the cycle time of the gas at that loss rate. If one film is both structure and barrier its
mass is counted once.

**Splitting water takes 237.141 kJ per mole of hydrogen of reversible work and 285.830 kJ of enthalpy at 298.15 K and
one bar (NIST-JANAF).** `gibbs_efficiency` divides the Gibbs work by the supplied work once; the difference from the
enthalpy is reversible heat. Each mole of hydrogen consumes a mole of water and yields half a mole of oxygen, vented
or used apart from the lifting gas.

## Thickness and replacement power

**Permeation input falls in proportion to thickness, and film mass rises in proportion.** The reproducible
`reference_thickness_scan()` uses a full sphere, 100 kPa, 298.15 K, 98% hydrogen and η = 0.5, at each material's own
measured permeability. Input in W per m² of projected area (the minimum Gibbs work is half; sphere area is four times
the projection):

| Comparator | 1 µm | 25 µm | 100 µm | 1,000 µm |
|---|---:|---:|---:|---:|
| Regenerated cellulose | 10.75 | 0.430 | 0.107 | 0.0107 |
| Polypropylene | 416.5 | 16.66 | 4.165 | 0.416 |
| Chitin, conditions incomplete | 1.495 | 0.0598 | 0.0149 | 0.00149 |
| MFC 10.7 Barrer, conditions incomplete | 666.4 | 26.66 | 6.664 | 0.666 |
| LMFC 12.1 Barrer, conditions incomplete | 753.6 | 30.14 | 7.536 | 0.754 |
| 100 Barrer sensitivity | 6,228 | 249.1 | 62.28 | 6.228 |
| Film mass at 1,500 kg/m³, kg per m² of skin | 0.0015 | 0.0375 | 0.150 | 1.50 |

The scan also holds 10 µm cases. Defects and thickness-dependent structure can break the inverse law in thin films:
the same study's 8 µm spray membrane cracked and passed 449 mA/cm². At fixed permeability and pressure, replacement
power per projected area is the same at every sphere radius, while the inventory per projected area grows with
radius.

**The aerial-ecology production of 1,000 g C/m²/yr is 1.2675 W/m² of net biomass energy** (0.45 carbon fraction,
18 MJ/kg dry). The low-permeability films replace their gas for a small share of it; at 100 µm, several of the
measured high-permeability films would take all of it.

**The humid S-CNF device's crossover is 1.55 × 10⁻⁶–1.55 × 10⁻⁵ mol/(m² s) of hydrogen,** 1.475–14.747 W/m² projected
of minimum Gibbs work over a sphere, or 2.949–29.494 W/m² at η = 0.5. For historical litre-based figures, one litre
per m² a day at STP would need about 0.98 W/m² projected at η = 0.5.

## Inflation

**A full 100-m-radius sphere holds 333.84 t of hydrogen, made from 2,983.23 t of water with 39.269 TJ of minimum
work.** At η = 0.5 that is 78.538 TJ, or 2.500 GJ per m² of projection; filling in one year would take 79.22 W/m²
before any permeation or construction. The input equals about 62.5 years of the aerial-ecology production over the
final projected area. A growing body fills in stages from a smaller collecting area, and its trimmed gas is a fraction
of the full volume; [growth](growth.md) integrates that path.

## Composition and the decisive measurement

`fixed_volume_exchange` integrates each species' exchange in a rigid, isothermal volume with no production or purge;
oxygen and nitrogen enter even while the total pressure inside exceeds the outside. It is a tested diagnostic.
Oxygen's slow ingress and its microbial scrubbing are costed in [storms](storms.md).

The decisive experiment is an assembled barrier tested for hydrogen, oxygen and nitrogen transport and water balance
under wet tension, temperature cycling, irradiation and biological exposure, with its seams, punctures and repairs.

## Checks

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q research/studies/floater_viability/test_envelope.py
```

The tests cover partial-pressure diffusion, film resistance, energy accounting, the crossover-current conversion,
independent gas volume and area, stoichiometry and the exact limits of species exchange.
