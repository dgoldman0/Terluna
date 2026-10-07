# Terluna light calendar

A static single-page almanac for the Open Moon, built on the
[date-based illumination products](../../illumination/calendar/). Its source
lives on `study/sea-appearance` alongside the research it displays.

Choose a coast or map point and a date from 2000–2500. The page provides surface
illuminance, separate direct and scattered sunlight/Earthlight, topocentric
positions and Earth's illuminated face. A monthly logarithmic graph, next
rise/set and brightness events, and daily summaries make the long lunar day
readable on an ordinary Earth calendar. The Earthlight switch recalculates the
brightness curve and threshold events. Links retain the selected viewpoint;
CSV exports include units, coordinates, time scale and eclipse flags.

The layout takes the location → current conditions → graph → calendar sequence
as its guide. The map, source balance and sky geometry address this world's
additional light sources. It uses local HTML, CSS, SVG, canvas and ES modules;
there are no runtime external libraries, tracking calls or font downloads.

## Build and serve

```sh
python visualization/light-calendar/build.py
python -m http.server 8080 --directory visualization/light-calendar/dist
```

Open `http://localhost:8080`. An HTTP origin is required for module workers;
opening `index.html` directly from disk is unsupported. `--out PATH` selects
another output directory. The build refuses missing, incompatible or stale
domain products and writes `build.json` with source/product hashes. Generated
`dist/` is ignored. Deploy that directory to any static host; no server-side
calculation or database is required.

The worker samples brightness every 30 minutes, refines bracketed events to
sub-second numerical precision, refines sampled turning points to catch grazing
double crossings, and reports sampled daily extrema. Event
**accuracy** remains limited by the orbital and lighting models; sub-second
root finding is not a claim of sub-second physical accuracy. The source model's evidence, time-scale
assumptions and eclipse limitation are available in the page's model details.

The physics reader is packaged unchanged from
`illumination/calendar/evaluator.mjs`. Model changes belong in that domain;
this folder owns controls, displays, exports and static packaging.

## Initial build verification

The calendar's 19 Python tests and nine JavaScript tests pass, including the
independent light-reader cases, grazing crossings, Earthlight switching, leap
days and the end of the supported calendar. Packaging checks the product and
source hashes. The repository layer check also passes.

`make check` was run before the commit: 705 Python tests passed, with four
failures from a missing `climate/gcm/products/moon_gcm_configuration.json` and
the execution environment's process lookup. A separate legacy JavaScript run
passed 176 tests and found one regenerated immersion PNG byte-size fixture
mismatch. The historical provenance check reports `PASS_WITH_BLOCKED_INPUTS`;
the ensemble setup validator passes. These are separate from the new calendar
checks. Browser-based visual and interaction QA was unavailable in this build
environment; the packaged modules and calculation paths were checked directly.
