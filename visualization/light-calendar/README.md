# Terluna light calendar

Two static pages built on the [date-based illumination products](../../illumination/calendar/):

| Page | For | What it shows |
|---|---|---|
| `index.html`, **Open Moon Skies** | Anyone curious about the Open Moon | The sky at one of six places now and through the coming month: the condition in plain words, one sentence about the Sun, Earth's phase and place, a comparison with light on Earth, and the next sunrises, sunsets, twilights and Earth phases |
| `almanac.html`, **Light almanac** | Research and checking | Any coast or map point and any date from 2000 to 2500: surface illuminance split into direct and scattered sunlight and Earthlight, topocentric positions, a monthly logarithmic graph, rise/set and brightness crossings, daily summaries and CSV export |

Both run in the browser from static files: HTML, CSS, canvas and ES modules with
a module worker. There is no server-side calculation, database, external library,
tracking call or font download.

## Build and serve

```sh
python visualization/light-calendar/build.py
python -m http.server 8080 --directory visualization/light-calendar/dist
```

Open `http://localhost:8080` for Open Moon Skies and `http://localhost:8080/almanac.html`
for the almanac. Module workers need an HTTP origin; the page opened from disk
says so and shows these two commands. `--out PATH` selects another output
directory. Deploy `dist/` to any static host.

The build checks the calendar products' schemas, producer hashes and constants
and refuses stale or missing ones. [assets.py](assets.py) then writes the sky
page's display assets into `dist/assets/` from checked products:

| Asset | Source | Check |
|---|---|---|
| `sky-moon.bin` | The solved spherical sky atlas of the 1.2-atm Open Moon ([illumination/sky](../../illumination/sky/)), log luminance and chromaticity on the atlas grid | SHA-256 recorded in `solved_sky.json` |
| `assets.json` | The atlas grid, the direct beam's colour by Sun height, and the light under the solver's Earth control sky, used for "as bright as Earth…" | Same product |
| `earth.jpg` | NASA Blue Marble surface and cloud maps, composited at 1024 × 512 | SHA-256 in [nubium-earth-inputs.json](../sea-appearance/nubium-earth-inputs.json); fetched if absent |
| `places.json` | Six places: four calendar coasts and two Apollo sites from the [conservation register](../../research/studies/conservation/) at the [geography atlas](../../geography/README.md#atlas-at-the-selected-water-share) sea level | Register and atlas sea levels must agree |
| `stars.json` | Yale Bright Star Catalogue stars to magnitude 5.5 ([illumination/stars](../../illumination/stars/)) | Product hash |

`build.json` records the hash of every source file, product and asset. The sky
atlas lives in the optical comfort study's run storage. When a source is
missing, the almanac still builds and the build reports which source the sky page
needs; `--require-assets` makes that a failure.

## Open Moon Skies

The page opens on the Ocean of Storms at the present moment, in the visitor's
own time zone. Choosing a place changes the whole page. Play runs through the
coming month at eight hours a second, the bar beside it shows the month's
conditions in colour, and each upcoming moment can be tapped to see its sky.

**Computed.** Every light level, Sun and Earth position, phase and event comes
from the calendar products through [evaluator.mjs](../../illumination/calendar/evaluator.mjs);
the worker [explorer-worker.mjs](explorer-worker.mjs) finds sunrises, sunsets,
the 2.98-lux twilight thresholds, Earth's phases and Earth's eclipses of the Sun.
The sky's colour and brightness in every direction come from the solved atlas,
scaled to the calendar's scattered light for the Sun and for Earth. Earth shows
its real rotation: [sky-frame.mjs](sky-frame.mjs) turns Earth's surface and the
stars into the observer's sky, and its sub-observer and sub-solar points agree
with the saved JPL Horizons record within 0.05°. The Earth disk carries the
calendar's direct Earthlight; stars are dimmed by the solved column's Rayleigh
depth.

**Words.** [words.mjs](words.mjs) names five conditions: daytime while any of
the Sun is up; golden twilight while the ground keeps at least 2.98 lux, the
light at the end of Earth's civil twilight in the same solver; Earthlit night
when Earth gives more light than the Sun's twilight; night down to 0.05 lux; and
dark night below that. Light is compared with the Sun height on Earth that gives
the same horizontal light under the solver's Earth control sky, and nights with
the brightest full Moon on Earth, about 0.3 lux ([Kyba, Mohar and Posch 2017](https://doi.org/10.1093/astrogeo/atx025)).

**Display choices** ([sky-render.mjs](sky-render.mjs)). The picture is a
panorama with azimuth across and elevation up at one scale, framed on the Sun,
its twilight glow or Earth. Its exposure follows the average sky in view the way
an eye adapts, so days are bright, twilights glow and nights stay dark; colour
fades toward a moonlit blue at night. Earth is exposed on its own, as a
photographer exposes the Moon, and reddens near the horizon with the direct
beam's colour. The foreground is open sea at coasts and seas and level land at
the Descartes Highlands, with the ground albedo of the light model. The sea
reflects the sky from facets tilted about 7°, and the Sun and Earth glitter on it
from Cox–Munk facets with a ripple texture. The Sun's glare is a bloom.

By day the solved sky is a pale cyan-blue overhead and cream toward the horizon,
and the direct Sun is golden; in twilight the sky toward the Sun turns gold and
orange. These colours are the atlas's own.

## Light almanac

Choose a coast or map point and a date. The page provides surface illuminance,
separate direct and scattered sunlight and Earthlight, topocentric positions and
Earth's illuminated face. A monthly logarithmic graph, next rise/set and
brightness events, and daily summaries make the long lunar day readable on an
ordinary Earth calendar. The Earthlight switch recalculates the brightness curve
and threshold events. Links retain the selected viewpoint; CSV exports include
units, coordinates, time scale and eclipse flags.

The worker samples brightness every 30 minutes, refines bracketed events to
sub-second numerical precision, refines sampled turning points to catch grazing
double crossings, and reports sampled daily extrema. Event accuracy follows the
orbital and lighting models; the root finding resolves crossings well below that.
The model's evidence, time-scale assumptions and eclipse treatment are in the
page's model details. The physics reader is packaged unchanged from
`illumination/calendar/evaluator.mjs`; model changes belong in that domain.

## Checks

```sh
python visualization/light-calendar/build.py
node --test visualization/light-calendar/tests/*.test.mjs
```

[calendar.test.mjs](tests/calendar.test.mjs) checks grazing crossings, event
thresholds, the Earthlight switch, leap days and the end of the calendar.
[sky.test.mjs](tests/sky.test.mjs) checks Earth's orientation against JPL
Horizons, finds 2026's two eclipses of the Sun by Earth over the Sea of Showers
on the dates of that year's lunar eclipses (3 March and 28 August), and checks
the twilight sentences, the brightness comparisons and the six places. Both pages
were reviewed in a browser at desktop and phone widths, without horizontal
scrolling. The first public page (commit 6745ee6) stays in history as a failed
attempt.
