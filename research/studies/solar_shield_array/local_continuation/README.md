# Local continuation of the expanded coupled cycle search

This package accompanies `../expanded_cycle.md`. It preserves an **unresolved
search**, not a feasible shielding cycle. The physical comparator remains the
zero-translation slower-turn prefix. The last trial and local restoration
models are diagnostic. Read the final report and `restart.json` before choosing
a candidate or increasing the budget.

## Environment and compact inputs

Run from the repository root on `research/solar-shield-habitat-array`.
Use Python 3.12, NumPy 2.3.5, SciPy 1.17.0, Shapely 2.1.2, jplephem 2.24
and a C++17 compiler. `requirements.txt` records the actually inspected versions.
The runner limits its address space to 2 GiB and uses one numerical thread.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r research/studies/solar_shield_array/local_continuation/requirements.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m research.studies.solar_shield_array.expanded_run inspect
```

`ephemeris.npz` contains geometric DE440 samples at 1,800-second spacing for
six days around the scenario epoch. They are derived from the pinned NASA/JPL
DE440s kernel, using the repository's original sampler and interpolator.
`inputs.json` includes the kernel URL, SHA-256, citation and epoch; the compact
file has its own hash. NASA/JPL Solar System Dynamics and NAIF receive credit.
No network retrieval or ignored kernel is required for ordinary replay.
The committed `joint_seed.json` and checked `coupled_restart.json` retain the
original exact states, hardware allocation and previous failed search.
`basis.npz` contains the residual-directed acquisition patterns and nominal
terminal centre references. The original sources and spectral products are
unchanged. The seed is deterministic; no random training or search is used.

The complete executed member histories are stored as two byte parts to fit
the upload limit. `inspect` verifies each part's hash. They contain the exact
original NPZ bytes; read them without writing another file:

```python
import hashlib, io, json
from pathlib import Path
import numpy as np
p = Path('research/studies/solar_shield_array/local_continuation')
h = json.loads((p/'restart.json').read_text())['executed_history']
data = b''.join((p/name).read_bytes() for name in h['parts'])
assert len(data) == h['size_bytes']
assert hashlib.sha256(data).hexdigest() == h['sha256']
history = np.load(io.BytesIO(data), allow_pickle=False)
print({name: history[name].shape for name in history.files})
```

## Replay and resume

The representative smoke command starts at the saved actual hour-6 state
and integrates 360 seconds. It is a loading/execution check, not validation.

```sh
python -m research.studies.solar_shield_array.expanded_run replay \
  --run-dir /tmp/terluna-cycle-smoke \
  --candidate research/studies/solar_shield_array/local_continuation/last_trial.json \
  --case smoke --end-hour 6.1 --cpu-budget 60
```

For a complete diagnostic replay and independent original-epoch validation:

```sh
python -m research.studies.solar_shield_array.expanded_run replay \
  --run-dir research/runs/solar_shield_array/local_cycle \
  --candidate research/studies/solar_shield_array/local_continuation/last_trial.json \
  --case retained --cpu-budget 1200
python -m research.studies.solar_shield_array.expanded_run validate \
  --run-dir research/runs/solar_shield_array/local_cycle \
  --candidate research/studies/solar_shield_array/local_continuation/last_trial.json \
  --case retained --cpu-budget 1200
```

For the larger search, allocate a cumulative CPU budget explicitly. Each
stage adds to the run-directory ledger; increasing `--cpu-budget` raises that
directory's cumulative ceiling. The Jacobian stage skips completed columns.
Completed columns survive interruption. A fresh baseline starts from the
original physical state and reproduces the anchor; it does not require old
Work trajectories. The saved `optimizer_model.npz` preserves the measured
local Jacobian and contacts for inspection. Rebuild after changing the anchor,
basis, forces, ephemeris, mass, cadence or contact branches.

To continue the saved local SQP model immediately, without regenerating any
ignored trajectories, use its exact coefficients and Jacobian:

```sh
python -m research.studies.solar_shield_array.expanded_run solve \
  --run-dir /tmp/terluna-saved-model \
  --resume-model research/studies/solar_shield_array/local_continuation/optimizer_model.npz \
  --radius 0.1 --case resumed --cpu-budget 60
```

This computes a proposal only. The documented coupled replay and validation
remain required; the original physical comparator is the retained state.

```sh
python -m research.studies.solar_shield_array.expanded_run benchmark \
  --run-dir research/runs/solar_shield_array/local_cycle --cpu-budget 12000
python -m research.studies.solar_shield_array.expanded_run sensitivity \
  --run-dir research/runs/solar_shield_array/local_cycle \
  --cpu-budget 12000 --reserve 1500 --columns 6 7 0 10 11 12 13 14 15 17 8 9 1 2 3 4 5
python -m research.studies.solar_shield_array.expanded_run solve \
  --run-dir research/runs/solar_shield_array/local_cycle --cpu-budget 12000 --case local
python -m research.studies.solar_shield_array.expanded_run replay \
  --run-dir research/runs/solar_shield_array/local_cycle --cpu-budget 12000 --case local
python -m research.studies.solar_shield_array.expanded_run assess \
  --run-dir research/runs/solar_shield_array/local_cycle --cpu-budget 12000 --case local
```

These commands describe one relinearization and trial. A larger iterative
optimizer must explicitly adopt/replay the new anchor, rebuild contacts and
Jacobian and compare actual/predicted merit again. The current runner never
silently treats a rejected trial as the new state. Arrival timing is an
available basis direction but changes a fixed nominal centre target; include
centre/receiver-path freedom before interpreting a failed timing solve broadly.

## Acceptance constraints that remain

Initial and next six-hour services must protect the complete moving 20 km
receiver for the finite Sun, with 9.89 km clear apertures, original seams,
physical orientation and unchanged EUV/ionizing spectral requirements. The
receiver stays inside four lunar radii. Sampled endpoint restoration is only
a proposal constraint; whole-service polygon coverage and missed-ray searches
remain mandatory acceptance tests.

Every square must retain 100 m clearance after numerical allowance, using
adaptive between-sample checks. Independently refine source/area quadrature,
integration step and tolerances. Check individual translation, rate, torque
couple and installed power limits. Energy ceilings remain 51.859 TJ at twelve
hours and 100 TJ through the next service; the latter is a search cutoff.
Any virtual control or numerical restoration slack must be zero. A trajectory
after its first failed clearance event receives no physical credit.

The finite-sequence repeatability gate is 50 m, 0.01 m/s and 1e-4 rad at
matched service and departure phases. Further actual-date cycles must continue
from the resulting state without reset, meeting coverage, clearance, resource
and bounded-defect gates. The nonautonomous ephemeris precludes inferring
indefinite recurrence from one endpoint match.

Actuator capacity and electrical supply are separate. No generator, storage,
bus or collector placement is installed here. Additional actuator mass at
300 W/kg, supply mass, optical interception and changed photon forces require
an original-epoch replay before hardware acceptance. Preserve the 238.35 TW
held-array comparison and 23.835 TW operating target. The 17.29-million-tile
number remains unvalidated; cadence, full-target coverage, inventory and
handover require separate fleet work.

The bounded Work run does not start the larger local campaign.
