# Waves arriving along the Smythii shore

The coastal model follows changing waves through a complete lunar solar cycle,
after a full preceding cycle of spin-up. It extends the selected steady
states in the [shore exposure study](shore.md) into a time history at the
western face, eastern face and southern shore of the eastern Smythii headland.
The [result product](results/coastal_history.json) records wave heights,
periods, incoming energy and every threshold interval. Selected spectra drive
the [individual-wave and run-up calculations](runup.md).

## What reaches each coast

The western face receives the largest sustained exposure in this simulated
cycle. The eastern face becomes active in a later weather phase. The southern
shore has substantial metre-scale episodes and lower shoreward power than
the western face. All durations below use the second cycle alone.

| Station | Peak Hs | Mean period at that peak | Hours with Hs ≥ 0.5 m | Hours with Hs ≥ 1 m | Longest ≥1 m span within the cycle |
|---|---:|---:|---:|---:|---:|
| Western face | 1.80 m | 12.24 s | 385.8 | 142.8 | 68.2 h |
| Eastern face | 0.97 m | 8.68 s | 77.7 | 0.0 | 0.0 h |
| Southern shore | 1.53 m | 12.42 s | 384.3 | 96.3 | 44.3 h |
| Local offshore reference | 1.95 m | 12.53 s | 622.8 | 203.3 | 70.2 h |

The western and offshore longest spans begin with waves already above the
threshold at the reporting boundary. Their complete episodes extend into
the preceding cycle. The eastern face's longest stretch above 0.5 m lasts
49.8 hours, starting 14.19 Earth days into the second cycle. Its peak follows
at day 14.71, about nine days after the principal western and southern peak.

The clearest new metre-scale episode begins at these locations:

| Location | First upward 1 m crossing, Earth days into cycle | Duration above 1 m in that episode |
|---|---:|---:|
| Local offshore reference | 5.208 | 50.9 h |
| Western face | 5.260 | 47.8 h |
| Southern shore | 5.305 | 44.3 h |

The western crossing follows the local offshore crossing by 1.24 hours;
the southern crossing follows by 2.34 hours. These are observed onset
differences within the calculation. Changing direction, local wind input
and the different exposure of each location accompany propagation.

| Coast | Maximum incoming power | Mean incoming power | Incoming energy through the cycle |
|---|---:|---:|---:|
| Western face | 546.5 W/m | 91.6 W/m | 233.8 MJ/m |
| Eastern face | 93.0 W/m | 7.3 W/m | 18.7 MJ/m |
| Southern shore | 243.4 W/m | 35.8 W/m | 91.2 MJ/m |

The western mean is 12.5 times the eastern mean and 2.56 times the southern
mean. This describes strong geographic and directional concentration of
wave energy. The units refer to a metre of the local depth contour at each
station; an alongshore total requires its spatially varying transport.

The [time-history figure](../../visualization/waves/results/coastal_history.png)
and [exposure maps](../../visualization/waves/results/coastal_history_maps.png)
show how those episodes sit within the complete cycle. The coastal stations
remain in 14–49 m water. The [profile study](runup.md) carries selected incoming
spectra onward through breaking to the moving waterline.

## Clock, geography and forcing

All three spectral models begin at hour zero and run through hour 1,419.
The reported window is **708.73416–1,417.46832 Earth hours**, the second
29.53059-day solar cycle. The first cycle supplies spin-up throughout the
basin, regional and coastal models. Report endpoints are interpolated to the
exact cycle boundaries.

| Scale | Propagation grid | Boundary and forcing |
|---|---:|---|
| Smythii–Marginis basin | About 30.3 km | The existing 150 s continuous-cycle calculation, replayed to retain hourly directional spectra |
| Eastern regional sea | About 1.90 km | Four nearby wet basin spectra at each wet perimeter point, with the established distance bounds |
| Headland and surrounding coast | About 474 m | Hourly regional directional spectra, with the same bounded wet-node interpolation |

Both finer models use the existing native 118 m LOLA/GRAIL terrain, sampled
at their propagation nodes. The coastal box covers about 30 km square.
The western, eastern and southern stations occupy water 22.38, 49.04 and
13.55 m deep respectively. They remain the same physical nodes across the
earlier spatial refinements.

The local offshore reference is at 92.84570°E, 2.34570°N, in about 1,271 m
of water, roughly 12 km west of the western-face station. Its occurrence
statistics describe this coastal reference. The earlier basin report uses
the more distant reference at 89.125°E, 1.125°S.

The atmospheric source is the recovered global three-hourly surface-stress
history. The existing coupling maps that stress into SWAN's drag law; its
15-minute forcing file is carried through both nests. This retains the
atmospheric model's resolved spatial and temporal variations. Local coastal
gusts and topographic wind circulations require finer atmospheric forcing.

Hourly boundary spectra are interpolated linearly in time by SWAN. The
initial regional and coastal timestep is 900 s. Station diagnostics are
written every 15 minutes, directional spectra every hour and full fields
every six hours. An additional final field brackets the exact reporting
endpoint. The map maxima therefore describe six-hourly samples; the station
maxima use the denser station history.

The propagation, wind input, gravity and seawater assumptions follow the
earlier studies. Water level is fixed, and the spectral model omits currents,
reflection, diffraction and sediment evolution. The SWASH profile model
adds reflected waves, bulk breaking and wetting and drying within its
explicit one-dimensional geometry.

## Reading arrival times and duration

For each stated height or incoming-power threshold, the analysis finds the
crossings of the piecewise-linear sampled history. It reports every interval,
the total duration, the longest continuous episode and whether the cycle
boundary cuts through an episode. These describe this particular simulated
cycle. Additional cycles are needed to measure variability between months.

An offshore-to-coast onset difference includes propagation, local growth and
changes of direction. The offshore reference and each shore also have
different exposure. A threshold crossing is therefore an event onset at a
stated location; identifying a pure travel time requires a controlled signal
or explicit tracing of its spectrum.

Incoming power integrates the directional spectrum times group velocity and
the positive projection toward shallower water, multiplied by water density
and gravity. Its units are watts per metre of local depth contour. Opposing
directions are retained separately. The offshore reference uses net transport
magnitude and carries that label in the product.

## Numerical evidence

The basin replay compares every retained wave value with the established
continuous-cycle calculation. Prepared boundaries retain the input spectrum,
upstream manifest, terrain and remapping hashes. Completed runs retain every
input and required output hash, their producer source and restart checkpoints.
The independent spectral integral is compared with SWAN's transport output
through the reporting cycle.

The replay preserves all 4,076,820 retained basin values exactly. A separate
restart control reproduces all 624 saved station values over its three-hour
continuation. Two-, four- and eight-thread coastal previews also produce
identical station and spectral output bytes.

Two coastal restarts halve the timestep from 900 to 450 s. Each comparison
allows 24 hours of adjustment after the saved checkpoint. The 720–912 h
window covers the beginning and end of the principal new metre-scale
episode. Its onset shifts by 0.47 minutes offshore, 2.60 minutes at the
western face and 5.20 minutes at the southern shore. All four stations pass
the stated height, period and transport bounds in this window.

The later 840–1,097 h comparison includes the changing exposure toward the
eastern face. Its 95th-percentile height differences are 1.26% west, 4.54%
east, 1.48% south and 1.47% offshore. Time above 1 m changes by 19.8 minutes
west, 0.26 minutes south and 0.70 minutes offshore across this window.
The eastern energy-transport vector difference reaches 13.83% at the 95th percentile,
exceeding the 10% bound. Restricting the paired eastern records to Hs ≥0.5 m
gives 8.12%; the complete-window failure remains in the
[check product](results/coastal_history_checks.json). Its upward 0.5 m crossing
shifts by 1.77 minutes; the comparison ends while this episode is still active.
These are selected coastal-window checks. Upstream and complete-cycle
refinements remain further work.

The evolving nests use one iteration per timestep, SWAN's default transient
setting. In the second cycle, the coastal log reports its 99% iteration
accuracy criterion unmet at all 2,835 recorded steps; the regional run meets
it at 68 steps. The timestep comparisons above provide a separate measure
of the resulting station sensitivity. Iterations per step and upstream
temporal resolution remain additional refinement dimensions. The
[SWAN numerical manual](https://swanmodel.sourceforge.io/online_doc/swanuse/node29.html)
defines this distinction between an iteration limit and an accuracy criterion.

The same window retains 2,424 coastal and 1,936 regional warnings for at
least one boundary point whose computed total Hs differs from the supplied
total by at least 10%. SWAN prescribes incoming directions while outgoing
waves evolve inside the domain, so this total-height comparison also includes
that outgoing field. The
[model guidance](https://swanmodel.sourceforge.io/online_doc/swanuse/node5.html)
describes this boundary behavior. The warnings remain in the run records;
evolving boundary-position and iteration controls are further work alongside
the earlier stationary boundary controls.

A controlled narrowband wave-energy pulse crosses a flat 10 km domain with
wind input and dissipation removed. The analytic group velocity is 1.45725 m/s.
At 8 km, the expected centroid delay is 1.52494 hours. The centroid agrees
within 0.0041 minutes at the 900 s timestep, while numerical spreading moves
the rising crossing of half the peak energy about 11.0 minutes early. Reducing the timestep
to 450 and 150 s reduces that rising-crossing bias to 6.59 and 3.55 minutes.
The pulse-height ratios at 8 km are 0.956, 0.972 and 0.984 respectively.
This separates a well-preserved mean delay from the broader event front.

The independent directional transport integral agrees with SWAN within
0.552% throughout the reporting cycle. Of the eligible hourly spectra,
37.2% exceed the 1% high-frequency edge bound. Their resolved heights range
from 0.10 to 0.75 m; every sampled metre-scale spectrum passes this particular
edge check. The failures occur 266 times at the west, 553 at the east, 214
at the south and 15 at the local offshore reference. This leaves the smaller
seas and their periods especially sensitive to wider frequency coverage.

The earlier spatial study remains part of the uncertainty: its 237→118 m
comparison has a 2.84% 95th-percentile height difference and a 39.09% maximum
at local outliers. The 474 m evolving calculation retains its own geographic
resolution. Weak, short waves also retain the established upper-frequency
coverage problem. Numerical checks and lunar physical calibration have
separate evidence states.

## Reproduction and storage

The runners preserve the existing checkout and the external wave-data link.
All spectra, binary checkpoints, full maps and model logs go to
`/media/projectspace/terluna-research/wave-runs/shore_history/` through
`research/runs/waves`. Compact results and reports stay in the repository.

```bash
python -m climate.waves.shore_history parent --executable <coupled-air-swan>
python -m climate.waves.coastal_history region --executable <coupled-air-swan>
python -m climate.waves.coastal_history coast --executable <coupled-air-swan>
python -m climate.waves.history_controls --executable <coupled-air-swan>
python -m climate.waves.coastal_history_checks
python -m visualization.waves.coastal_history
```

A selected-window timestep check restarts from a saved coastal checkpoint:

```bash
python -m climate.waves.coastal_history coast --executable <coupled-air-swan> \
  --step 450 --start <checkpoint-hour> --end <end-hour> \
  --initial <coastal-case>/h<checkpoint-hour>.hot --name <refined-case>
python -m climate.waves.coastal_history_checks --case <coastal-case>/product.json \
  --compare <refined-case>/product.json --start <comparison-start> --end <end-hour> \
  --output climate/waves/results/coastal_history_checks.json
```

The comparison begins at least 24 hours after the checkpoint and stays within
the second cycle. The regional boundary and wind history remain the same in
the paired runs. The completed cases restart at hours 696 and 816 and compare
720–912 and 840–1,097 h respectively. Add `--append` to retain both windows in
the comparison product; the analysis recomputes the existing entries from
their recorded inputs.

The [verification record](shore_history_verification.json) collects input and
producer hashes, numerical outcomes, rendered-array checks and repository
validation for this study.

The parent uses two CPUs; each finer spectral run uses eight, with a 2 GiB
address-space cap. The renderer checks its displayed arrays against the
hashed domain product and writes image/provenance pairs outside Git.
