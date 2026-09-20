# Climate

**Current material:** an executable scalar thermal-response calculation inside the [September feasibility model](../research/baselines/feasibility/model.py), with its original [18-case table](reference/thermal_response.csv).

It solves the periodic response of `C dT/dt + B T = F1 cos(omega t)` using imposed participating atmospheric heat capacity, water mixed-layer depth, a linear outgoing-radiation slope and the slow solar cycle. Temperature amplitude is a conditional response, not a geographical climate prediction.

**Missing:** a radiative-convective column, latitude-longitude heat transport, dynamically predicted winds, humidity, clouds, rainfall, persistent fog and storm statistics. None is represented as implemented.

## Next work

Combine a small set of gas/spectrum/water scenarios with time-dependent land/water energy balances. Check periodic solutions against the current analytic case before adding transport. Compare heat-storage and transport timescales and test sensitivity to cloud/albedo assumptions. Couple to measured terrain and basin alternatives from [geography](../geography/).

A specialist GCM is a setup-and-benchmark candidate, not an installed or verified capability of this repository. Regional climate descriptions in seeds and historical training data remain proposals.

Run `python research/baselines/feasibility/model.py --out research/runs/feasibility` from the repository root to reproduce the existing diagnostic with the other original cases. See [status](../research/status.json).
