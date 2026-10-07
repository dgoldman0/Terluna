# Date-based light calendar

This domain supplies the [single-page calendar](../../visualization/light-calendar/).
It combines the date ephemeris with the solved design atmosphere and separate
solar and terrestrial source spectra. The browser needs no application backend:
the expensive scattering solve is performed when producing the data.

## Products and geometry

`results/astronomy.json` (`terluna.illumination.calendar-astronomy/1`) exports
the compact lunar theory, shared constants, scenario coast coordinates, a water
mask and Python geometry samples. `geometry.py` subtracts the surface observer
from the lunar-centred source vectors before calculating elevations, azimuths,
angular sizes, distances, Earth's phase and the orientation of its bright limb.
This includes Earth's roughly quarter-degree surface parallax near the limb.

The supported Gregorian calendar is **2000–2500 inclusive**. The orbital input
is TT. The browser's UTC conversion uses historical leap seconds through 2017,
then assumes TT−UTC = 69.184 seconds; it cannot predict future civil-time changes.
The [JPL comparison](../../geography/CALENDAR_EPHEMERIS.md) contains 5,021 rows
at 5,010 unique epochs, with sparse coverage across the full interval and denser
event windows. Maximum sampled angular differences are 0.0592° for the Sun and
0.0559° for Earth. Tested ordinary centre-horizon crossings differ by up to
7.48 minutes. Near-grazing brightness crossings can differ by about an hour.
These sampled comparisons do not establish a uniform error bound at every date.

## Spectral light and finite disks

`results/transfer.json` (`terluna.illumination.calendar-transfer/1`) stores direct
and diffuse **horizontal illuminance in lux** per source elevation, and Earth's
response by phase. It also carries band-integrated source irradiance and the
apparent disk-reflectance spectrum at each phase. It uses the solved 1.2-atm molecular column, titania-stack
solar transmission multiplied by 0.95, and the established spherical scalar
Rayleigh operator with a Lambertian ground albedo of 0.1. Its spectral reduction
retains line absorption and targets a one-percent training budget. The
independent sampled beam check reaches a 1.11% lux difference; this is distinct
from an error bound for the complete scattering model. Earth reweights the
reduced channels using its calibrated irradiance in 10-nm bands.

`transport.py` reuses the sky solver's grids, quadrature, shell paths and
scattering operator. Its point-source driver retains scattering orders until
the increment criterion is met across the full −90° to +90° range, including
deep night. Source and input hashes, numerical settings and reduction checks
travel with the product. The original finite-Sun solver and its saved products
remain separate historical calculations.

The reader integrates a uniform finite Sun and a phase-dependent Earth disk
over that response. Earth uses [Glenar et al.'s spectral shape](https://doi.org/10.1016/j.icarus.2018.12.025)
at [Robinson et al.'s observed visual phase brightness](https://arxiv.org/abs/2507.22258),
with the existing model's continuous Lambert-shaped continuation beyond 144°.
Earthlight is not passed through the solar shield. The spatial disk weights
follow Lambert incidence and are normalized, so the empirical phase law is
applied once. The model resolves the illuminated hemisphere, angular extent,
bright-limb orientation and horizon clipping; it does not resolve Earth's
individual clouds or continents. Each source element uses parallel rays in
the horizontally uniform atmosphere, with topocentric distances and angles
at the observer.

The shared browser reader is `evaluator.mjs`; the visualization packages it
without reimplementing the model. Linear interpolation is in source elevation
and phase, before finite-disk quadrature. Distances rescale the sources, and
Sunlight and Earthlight are added in linear units. There is no stellar or
airglow floor. The logarithmic chart is only a display choice.

## Conditions and exclusions

These are conditional clear-sky predictions for a spherical, level surface,
using one atmospheric column at every site and date. Terrain obstruction,
observer altitude, refraction, aerosols, local weather, stars and airglow are
additional inputs. Earth is a mean source: its daily weather and surface
changes are not forecasts. Spectral phase extrapolation and the unresolved
disk pattern contribute uncertainty separate from the ephemeris comparison.

When Earth overlaps the solar disk, the reader sets `eclipse=true`. Its fluxes
remain the **unobscured reference**, since partial and total eclipse shadow
transport through this atmosphere has not been solved. The app, calendar and
CSV retain that distinction. Rise/set events refer to a geometric disk edge,
not a refraction-corrected or terrain-dependent visible horizon.

## Reproduce

The offline producer uses `research/requirements.txt` plus
`illumination/sky/requirements.txt`; the published app needs neither Python
nor these packages. Restore the byte-pinned atmospheric spectroscopy and solar spectrum through
their fetchers and the Earthshine table through `illumination.earthlight.fetch_inputs`.
The latter is [Mendeley Data v2](https://doi.org/10.17632/xfjm6nmh3m.2),
Glenar, Stubbs, Schwieterman, Robinson and Livengood, licensed
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Inputs are ignored; hashes are checked rather than silently updated.

```sh
python -m illumination.calendar.portable
OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=4 python -m illumination.calendar.build --compute --export
OPENBLAS_NUM_THREADS=1 python -m illumination.calendar.validate
python visualization/light-calendar/build.py
```

The heavy optics and point-source solution live in ignored
`research/runs/light_calendar/`. Only compact products are committed. Existing
caches must match their recorded producer/input identities. To rerun a changed
solver, use a fresh `--directory`; to export an unchanged completed solution,
use `--export` alone. `validate.py` records brightness propagation, finite-disk
quadrature comparisons and Python golden light values for the JavaScript tests.

The initial export converged in 182 scattering orders over 114 channels. Its
last retained contribution is at most 7.68×10⁻⁷ of solar diffuse surface lux.
This stopping measure is distinct from total numerical or physical accuracy.
The saved audit has 66 portable light cases. At 2,507 JPL epochs and eleven
locations, changing both directions and ranges changes total light by at most
1.2375% among non-eclipse comparisons whose reference is at least 0.1 lux.

The order-12/order-24 disk comparison changes total light above that threshold
by at most 0.00159%; it tests disk quadrature, not the atmospheric grid or source
model. The historical sea-calendar curve differs by up to 1.816% in daylight
and 2.335% in deep night. The largest daylight change mainly comes from replacing
its coarse log interpolation of the direct beam. Deep-night differences include
the retained scattering orders and convolving the piecewise-linear elevation
response over the finite disk. The audit retains both curves and their distinct
producer settings.

```sh
python -m pytest illumination/calendar geography/tests/test_calendar_ephemeris.py
node --test illumination/calendar/tests/*.test.mjs
```
