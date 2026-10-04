# Solar shield and habitat array

The author directed this work on 4 October 2026: develop the combined solar
shield, power system and habitat fleet, and investigate the holding problem
with a full Sun–Earth–Moon ephemeris. The branch starts at main commit
`2f2e0a104a3a1ea962a6dead4179f522762f997e`.

[design.md](design.md) records the proposed architecture before dynamics work.
The [full-ephemeris study](report.md) compares the September held screen with
conserved solar momentum, variable sunward distance, transverse motion,
eclipses and propagated control. The selected monthly trajectory reduces the
50-g/m² array's mean holding power from about 309 to 238 TW under the reference
propulsion assumptions. It still consumes about 371,000 kg/s of propellant.

[Photogravitational cycling](cycling.md) now has propagated point-tile
candidates at the same 50-g/m² base mass. A 15,000-km retrograde seed with
filtering during service and feathering between passes stays bounded for
three years and averages about 15.7% useful projected area. A separate
optimized 8.9-day solar-sail arc closes within 17 m on independent replay.
These establish orbital motion and an inventory comparison.
[tacking.md](tacking.md) retains the initial literature review.

The [simultaneous fleet study](fleet.md) now tests five common-epoch populations
with multiple planes and radii, finite square shadows, solar-disk quadrature,
seams, swept collision and mutual-shadow exclusions, and Earth-safe ideal
reflections. Continuous coverage remains unachieved. A 3,072-member large-square
geometry probe averages 61.76% interception before exclusions; its selected
19-member subset averages 0.78% and reaches zero. Those 1,000 km squares need
extended-body dynamics; the 10 km replay also fails the 50 m phase-accuracy
gate. No sufficient inventory or collector/habitat power is demonstrated.
Four lunar radii remains primary, with three radii as a comparison.

The study couples protection, engineering, illumination and habitation. Its
runners, compact numerical results and interpretation belong together here;
dynamics components belong in the protection domain. Computation files and
external ephemeris kernels go in ignored `research/runs/`. The historical
protection implementation and its imported results remain byte-pinned.

## Reproduce the study

```sh
python -m pip install -r research/studies/solar_shield_array/requirements.txt
python -m protection.dynamics.ephemeris --download
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.run
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m research.studies.solar_shield_array.validate
python -m research.studies.solar_shield_array.publish
python -m research.studies.solar_shield_array.tacking_screen
python -m research.studies.solar_shield_array.cycling_replay
python -m pytest protection/dynamics research/studies/solar_shield_array
python visualization/solar-shield-array/plot.py
python visualization/solar-shield-array/cycling_plot.py
```

The runner writes raw results to `research/runs/solar_shield_array/`.
`publish.py` checks source identity, relevant constants, time/area convergence
and the matching validation product before exporting compact
[holding](results/holding.json) and [validation](results/validation.json)
snapshots. The plot and its provenance are generated locally under
`visualization/solar-shield-array/results/`. `--quick` provides the centre-force
pilot without optimization, nodal-span propagation or the final exporter.

The final comparisons preserve a four-lunar-radius protected region and the
selected climate spectrum. Optical control is an optimistic lower bound;
10-g/m² films require a replacement optical stack. The 19-year test and ideal
feedback propagation establish numerical behavior for the tested families.
Final holding architecture, actual filter/PV forces, safe plumes, tile seams,
habitat capacity and long-term resource closure remain open.

## Checks

`make check` was run before this design-only commit: 545 Python tests passed,
55 skipped and 13 failed. The failures concern absent solar-spectrum inputs,
CM1 build configuration, the ExoPlaSim process lookup and ring-comfort tests.
No executable model was changed in that initial commit. Current targeted
and repository-wide check status, including the fleet phase-accuracy failure,
and the numerical
verification are recorded in [checks.json](checks.json) and
[validation.json](results/validation.json).

The cycling [report](cycling.md) gives the full run matrix, convergence limits,
source/constant provenance and reproduction commands. Its annual useful-area
statistic is numerically well resolved; orbit phase and perturbed tile
positions drift by kilometres, so metre-scale seam control is not established.
