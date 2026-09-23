# Climate

The original [scalar thermal-response table](reference/thermal_response.csv) and [September baseline](../research/baselines/feasibility/model.py) are preserved.

[cycle.py](cycle.py) adds a finite-volume latitude–longitude energy-balance model. It uses an equal-area spherical grid, moving sunlight, prescribed land/water heat capacities and a conservative diffusive heat-transport operator. An exact linear interval propagator and a periodic-boundary solve remove arbitrary initial-temperature/spin-up choices.

## Interpretation

This is an uncalibrated thermal screen. Synthetic dry, concentrated-water and distributed-water maps are compared. Albedo, atmospheric participation, transport and a linear outgoing-radiation law are inputs. The global mean follows analytically from those inputs; a temperate mean is not independent evidence that a 1.2-atm lunar atmosphere produces it.

Moisture, clouds, winds, latent heat, freezing feedbacks, actual terrain and photochemistry remain outside the calculation. Outputs reaching freezing or far outside the reference temperature range are flagged. The reported 273–313 K fraction is a chosen grid-time diagnostic, not a habitability or liquid-water certification.

The [54-case output](../research/studies/environment_screens/results/climate_budget.csv) closes its discrete global energy balance and periodic condition. Timestep tests pass, while a 6x24 to 12x48 spatial comparison changes some coarse-cell averaged temperatures by about 2 K; spatial resolution uncertainty remains. [Findings](../research/findings.md) explain the exact checks and boundaries.

## Next step

Replace the assumed outgoing-radiation law with a calibrated vertical radiative-convective calculation, then add water/ice treatment and measured terrain. The present temperature/irradiance interface can supply biological requirement calculations, with all atmospheric and trait assumptions retained.

Run `python -m pytest climate` and `python -m research.studies.environment_screens.run` from the repository root. This folder supplies no GCM weather prediction.
