# Open Moon — A month of light

> **Repository layout:** this directory tracks the reproducible application source,
> methods, third-party notices and numerical check records. The generated offline
> HTML, `.npz` sky atlases, CSV outputs and PNG textures/screenshots are excluded by
> `.gitignore`. They are included in the complete full-cycle distribution supplied
> with this project; a clean checkout builds them using the commands below.
> `SOURCE_BUNDLE_MANIFEST.json` preserves the original distribution's file hashes,
> including its original README. `repository_import.json` records this import.

A self-contained 360° full-solar-cycle comparison viewer, extending the earlier
*The long light* experiment. Open `Open_Moon_Full_Cycle.html` in a current browser.
The application and its numerical data are embedded. No web server or external
JavaScript libraries are required. WebGL 2 provides the full visual renderer;
a software renderer uses the same clear-sky maps with simplified cloud layers.

## Use

The initial state is **lunar noon with clear air**, with the Sun at +90°.
Drag to look around; scroll to zoom. `Overhead` reveals the overhead Sun;
`West / bay` frames the cove and `East / forest` looks the other way.
`Follow Sun` keeps the solar direction in view, or the horizon when it is below it.
`H` hides the controls, and `Escape` restores them. Arrow keys move one simulated
hour through the selected clock. The timeline runs from noon to the next noon.

`Moon`, `Earth` and `Compare` share the same solar phase and chosen conditions.
The lunar clock takes 29.53059 Earth days; the Earth clock takes 24 hours.
Playback speed changes the solar clock only. Cloud, water and rain motion have a
separate animation clock. `Conditions` exposes weather inputs, animation, loop
control and graphics quality. High quality raises rendering resolution and, on
WebGL, the cloud integration sample count.

The six weather presets are Clear, Scattered cloud, Overcast, Haze, Fog and Rain.
Their coverage, optical-depth parameter, visibility and movement are prescribed.
They do not specify expected lunar weather frequencies, cloud microphysics,
rainfall, wind statistics or a coupled climate solution.

## Evidence boundary

The clear-sky calculation has been extended from +12° to +90° and recomputed for
all three optical profiles. It includes spherical straight-ray spectral scalar
molecular transport, repeated scattering, ozone absorption, a finite solar disk
and a neutral reflecting ground. At lunar noon the ozone-bearing scenario gives
about **94.2 klx** (60.5 klx direct plus 33.7 klx diffuse). These are conditional
outputs of the chosen optical profiles.

The scene is an illustrative wooded cove, with local tangent geometry. Its
geography and organisms are visual reference choices. It is neither a DEM nor a
validated ecosystem. The sky solver retains the relevant planetary radius.
Scenery, cloud shading, cloud shadows, haze/fog, wetness and rain use additional
reduced display approximations. Readouts retain the **clear-sky baseline**;
there is no claim that weather-altered lux has been scientifically validated.

The light source is the Sun. Earthlight, starlight, airglow and artificial
illumination remain outside the calculation. The night extension below −35°
continues to use the coarser diffuse-source grid. A shared camera-style automatic
exposure is available; it is not an eye adaptation or night-vision model.

## Actual checks

- 41 unit tests pass (21 retained numerical tests plus 20 new full-cycle/data tests).
- Six new noon radiance comparisons use 100,000 independently traced Monte Carlo
  paths each. Selected differences range up to about 0.8%; some exceed two
  Monte Carlo standard errors. These are limited numerical comparisons, not a
  whole-sky error bound or empirical lunar-atmosphere validation.
- The exact vertex and fragment shaders compile and link in Mesa OpenGL ES 3.2,
  with rendered clear/cloudy/fog/rain scenes inspected separately.
- The browser controls and embedded-data loading are exercised in Chromium's
  software fallback. This environment disables WebGL and local-file navigation;
  the test inserts the complete HTML into the page. No external network request
  is needed. No claim of hardware-browser WebGL testing is made.

The preceding study reported a global lunar energy accounting residual up to
4.6%. This remains unresolved. Its records are retained in `validation/legacy_v1/`.

## Rebuild

```sh
python -m pip install -r requirements.txt
mkdir -p work
NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 python calculate_full_cycle.py
NUMBA_NUM_THREADS=4 python make_landscape.py --width 3072
python -m unittest test_solver test_full_cycle -v
python build_cycle.py
```

To check the new noon samples: `python check_noon.py`. To run browser controls:
`python browser_cycle_check.py` (requires Playwright and Chromium). The optional
`egl_scene.py` / `render_preview.py` checks require Mesa EGL. The generated maps,
geometry textures, source, methods, checks and third-party notices are included.
Large intermediate source-field caches can be regenerated and are omitted from
the compact source archive. No GitHub Actions workflow is involved.
