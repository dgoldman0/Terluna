# Lighting-calendar geometry, 2000–2500

The existing [lunar ephemeris](lunar_ephemeris.py) has been compared independently
with NASA/JPL Horizons across **2000-01-01 00:00 to 2500-12-31 23:30 TT**. No
ephemeris coefficients were changed or fitted. The saved experiment supports a
2000–2500 calendar with the measured angular accuracy below. These are sampled
maxima, not rigorous bounds at every intervening instant.

[calendar_ephemeris_check.py](calendar_ephemeris_check.py) fetches and reproduces
the check. [The CSV fixture](calendar_ephemeris_check.csv) pins the reference
values, [input provenance](calendar_ephemeris_inputs.json) records the exact API
queries, returned frames and response hashes, and [the result](calendar_ephemeris_check.json)
contains all residual statistics and event comparisons. The fixture occupies
about 0.52 MB and requires no network for routine checks.

## Sampling and reference conventions

The 5,021 rows contain 5,010 unique epochs:

- A 73-day grid across 2000–2500, with a separate end-of-2500 direction sample.
  Geometric Sun/Moon and Earth/Moon ranges are checked at 2,507 grid epochs.
- Fourteen 24-hour windows, sampled every 30 minutes around a sunrise and sunset
  in 2000, 2050, 2100, 2200, 2300, 2400 and 2500. Their sites cover both
  hemispheres, both near/far sides, and latitudes from the equator to ±80°.
- Every half hour from March 1 through April 3, 2030, for twilight-brightness
  threshold comparisons at the six coast presets and five additional locations.
- Two one-minute follow-up windows at western Procellarum, where the 0.1-lux
  threshold nearly touches the nightly brightness minimum.

Horizons uses **DE441** for the positions and **MOON_ME**, the DE421
mean-Earth/polar-axis cartographic standard, for the lunar orientation. The
reference is east-longitude positive. The fetcher rejects an observer response
that substitutes IAU_MOON. Horizons' modern principal-axis solution and its
standard mean-Earth cartographic frame have distinct roles; MOON_ME is the
applicable frame for these longitude/latitude comparisons.

Every input and output epoch is **TT**. Horizons performs its internal TT-to-TDB
conversion; the periodic difference is under two milliseconds. UTC inputs need
32.184 seconds plus the applicable TAI−UTC offset. Future leap seconds are
unspecified, so extending a fixed offset into the future is an explicit civil
calendar convention. It does not reduce the TT interval checked here.

Horizons observer quantities 14 and 15 are apparent sub-observer and sub-solar
coordinates as seen from Earth's center, including their documented light-time
corrections. The lightweight ephemeris provides instantaneous directions. Its
reported residuals therefore include this convention difference together with
truncation and omitted physical libration. The separate vector requests use
`VEC_CORR=NONE` at the same TT epochs, so range residuals compare geometric
distances consistently.

## Results

| Quantity | Largest sampled difference |
|---|---:|
| Sun direction, full angular separation | 0.059158° |
| Sun longitude / latitude components | 0.050588° / 0.041724° |
| Earth direction, full angular separation | 0.055875° |
| Earth longitude / latitude components | 0.042838° / 0.042439° |
| Earth distance | 13.88 km |
| Sun distance | 12,276.63 km |
| Earth distance contribution to inverse-square flux | 0.0075% |
| Sun distance contribution to inverse-square flux | 0.0163% |

The fourteen selected Sun-center horizon crossings differ by at most **7.48
minutes**. Halving the reference temporal resolution from 30 to 60 minutes
changes their interpolated event times by at most 0.121 seconds. Those events
use a level spherical horizon and Moon-center source directions. Surface
parallax, finite source radius, terrain and atmospheric refraction are separate
effects; near-tangent crossings can have substantially larger timing residuals.

## Propagation into the saved lighting curve

The checker reads the committed sea-appearance product's
`sky.ground_lux_by_sun_elevation`, pins its hash, and linearly interpolates lux
against elevation. It applies both ephemerides' solar directions to the same
curve at eleven locations. This isolates the brightness difference caused by
geometry within that conditional clear-sky model. It does not test the
radiative-transfer solution, regional atmospheric variation, or Earthlight.

| Reference solar elevation | Largest sampled relative illuminance difference |
|---|---:|
| Above 0° | 0.50% |
| −6° to 0° | 0.75% |
| −30° to −6° | 1.25% |
| −60° to −30° | 1.43% |
| −90° to −60° | 1.49% |

The March 2030 sample contains 73 crossings of 100, 10, 1 or 0.1 lux across those
locations, with the same crossing counts in both geometries. At western
Procellarum (28.875° N, 73.875° W), the 0.1-lux level is barely crossed. The
lightweight geometry predicts darkening **55.75 minutes early** and brightening
**60.52 minutes late** on March 10. The one-minute versus two-minute reference
interpolation changes either crossing by less than 0.15 seconds. The longest
remaining reference interpolation sensitivity is 5.80 seconds.

Small brightness residuals can therefore produce an hour of threshold-timing
shift when the brightness curve nearly touches a threshold. A uniform
minute-level accuracy claim would exceed this evidence. All threshold times in
the report are TT, and their definition uses this precise exported photometry
curve; another interpolator or atmosphere can produce different crossings.

## Reproduction

From the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m geography.calendar_ephemeris_check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest geography/tests/test_calendar_ephemeris.py -q
```

To explicitly refresh the reference fixture from JPL, append `--fetch` to the
first command. Requests are sequential and their total is 21. New references
are saved only after every request succeeds. A changed saved CSV fails the
offline hash check until explicitly re-pinned through a successful fetch.

Primary reference documentation:

- [NASA/JPL Horizons manual](https://ssd.jpl.nasa.gov/horizons/manual.html),
  especially the lunar mean-Earth frame, time scales and quantities 14/15.
- [NASA/JPL Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html), including
  `TIME_TYPE`, observer quantities and geometric vector tables.

This run used one compute thread and only the targeted tests above. Full
repository validation belongs to the integrating change.
