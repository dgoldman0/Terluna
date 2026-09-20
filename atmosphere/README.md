# Atmosphere

Two distinct reduced-model lineages are available. Both retain their own assumptions.

| Material | Condition | Appropriate use |
|---|---|---|
| [April simulator](sim/sim1.py) and [report](sim/report.md) | Existing repository history; source inspected in the preceding audit, not rerun in this organization pass | Sensitivity to prescribed cooling floor/law, diffusive structure, parametric photochemistry and adjustable Jeans/Parker interpolation |
| [September implementation](../research/baselines/feasibility/model.py) | Executable prescribed-temperature hydrostatic/Jeans baseline; imported without code changes | Spherical inventory and structure, molecular losses and imposed atomic-composition stress cases |
| [Reference tables](reference/) | Original September outputs | Compare identical assumptions; failed static closures remain in the tables |

The September script also owns the original climate and industrial calculations linked from their topic folders. Its [model notes](../research/baselines/feasibility/README.md) state the limits and rerun command.

## Interpretation and conflicts

The April solver's default 240 K radiative floor, assumed hot-state calibration, photochemical loss scales and adjustable outflow coupling remain inputs. Its temperature result is conditional on those choices. The September baseline instead prescribes a cold middle/upper profile and uses a different composition and exobase treatment. Recency and a temperature solver alone do not rank the two models by physical validity.

The existing 1.2 atm, 17.5% O2 mixture and literal Earthlike proportions are separate scenarios. A universal 250/260 K retention cliff is superseded by the structure/composition-dependent work. The approximately 250 K target has not been established as the outcome of a realizable filter.

## Next executable work

Reconcile pressure, composition, collision partners, homopause and temperature structure. Connect wavelength-dependent protection to heating and chemical source terms; test an energy-balanced reduced column with alternative upper boundaries. Track water supply to escaping regions separately. New results require conservation, convergence and benchmark checks.

No global climate, self-consistent photochemistry/fluid-kinetic solution or radiation-dose simulation is present in this folder. [Research plan](../research/plan.md) prioritizes the core's whole-world requirements.
