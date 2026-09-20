# Atmosphere

Three reduced-model lineages are available. Historical code and reference tables remain unchanged.

| Material | Condition | Appropriate use |
|---|---|---|
| [April simulator](sim/sim1.py) and [report](sim/report.md) | Prescribed cooling floor/law, parameterized photochemistry and adjustable outflow | Historical sensitivity baseline |
| [September implementation](../research/baselines/feasibility/model.py) and [reference tables](reference/) | Prescribed temperature and atomic fractions | Spherical inventory and conditional escape comparison |
| [Thermal column](thermal_column.py) | Numerically solved conduction, advected energy, molecular Jeans boundary and exobase | Conditional response to deposited sensible heat and a prescribed lower atmosphere |
| [Spectral interface](spectral_interface.py) | Explicit band energy ledger with missing-data rejection | Connect supplied irradiance and response parameters to surface-normalized heat; audit input coverage |

## New calculation and remaining closure

The thermal column solves its upper temperature without a fixed exobase floor. Its lower temperature/pressure, conductivity proxy, collision cross-section, molecular mixing ratio and heating distribution remain inputs. Atoms, ions, explicit radiative cooling, species chemistry, eddy diffusion, acceleration and tides remain outside this limiting model. Kinetic and tidal approximation warnings accompany outputs. Numerical convergence is not a test of physical stability.

The original solar table starts at 202 nm. The inherited optical implementation leaves a response interval between its X-ray treatment below about 24.8 nm and titania optical constants starting near 120.18 nm. A historical five-band solar benchmark is now transcribed separately, with all atmospheric/filter responses left unspecified. This recovers a coarse energy reference while preserving the unresolved optical and chemical requirements.

Read [findings](../research/findings.md) for equations, bounds, source limitations and numerical examples. The new [232-case output](../research/results/environment_screens/molecular_columns.csv) is separate from the byte-pinned historical tables. Molecular losses do not constitute total escape or a billion-year lifetime estimate.

## Run and next task

From repository root:

```sh
python -m unittest discover -s tests -v
python research/run_environment_screens.py
```

Next: acquire defensible EUV material response and resolved irradiance; replace assumed absorption/heating fractions with species-resolved absorption, chemistry and cooling; then test atomic transport and kinetic boundaries. Reconcile a calibrated lower-atmospheric profile with this upper model before inferring integrated habitability.
