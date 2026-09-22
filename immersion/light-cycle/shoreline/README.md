# Open Moon — atmospheric cloud studies / revision 04

## Accuracy status and next milestone

The revision-04 atmospheric column is numerically checked against idealized
reference calculations. Its prescribed cloud shapes, retained condensate and
approximate scattering remain separate physical uncertainties; overall image
accuracy is unbounded. See [Cloud accuracy assessment](CLOUD_ACCURACY.md) for the
limitation register and the measured 128/8 integration errors, and
[A1: next accuracy milestone](CLOUD_NEXT_ACCURACY_MILESTONE.md) for the ordered
work, deliverables and proposed acceptance gates. All A1 tasks remain planned.
The [audit summary](validation/clouds/accuracy-audit-summary.json) retains the
recorded model samples and comparison metrics with original-artifact provenance.
[Publication checks](PUBLICATION_04.json) distinguish fresh build/unit checks from
historical browser and GPU evidence.

Based on published shoreline commit
`032765c51b1aeb140c732e8f484e47efda15cf4d`, using the same verified Three.js r186
modules, spectral clear-sky data, terrain and surface assets.

## Explore

Open `Open_Moon_Shoreline.html`. Choose **Show controls → Conditions**, then an
**Vertical atmosphere study**: Coastal fog, Elevated clouds or Deep convection.
**Look into the sky** raises the camera. The world selector compares the Moon and
Earth; the Sun slider changes illumination while holding the selected sounding.
The panel displays cloud bounds, temperature profiles, scale height, liquid/ice
water paths and optical depth. **Export column** saves the inputs, numerical
profile and assumptions. **M1 noon** returns to the published clear-noon camera
and lighting. **Resume motion** advects the selected cloud field with its wind
shear, holding the sounding fixed.

Cloud support comes from condensation, buoyancy and the selected inversion/cap.
The calculated lunar ceilings range from about 123 m for the fog study to 59 km
for the deep convective sounding. These are prescribed experiments, not climate
predictions. See `CLOUD_METHODS.md` for the thermodynamics, optical assumptions,
rendering budgets and open work.

The **authored weather replay** remains separately selectable. It retains the
preceding rain episode and its older cloud geometry. Atmospheric-column studies
hold surface rainfall at zero; diagnosing cloud condensate does not supply a
precipitation forecast. The TSL renderer lab remains a clear-sky experiment and
does not yet render the new column clouds.

## Build

The source package contains the selected original r186 modules, their license,
reference surface images and spectral sky atlas. Build offline with Python:

```sh
python build_landscape.py
```

The builder checks dependency and surface-image hashes and produces the WebGL
viewer plus `Open_Moon_Renderer_Lab.html`. The package's `build.py` is only a
convenience wrapper. The repository's existing source-first builder is preserved
by the revision patch.

For a fresh repository checkout, restore the ignored dependencies and surfaces:

```sh
python tools/import_three.py /path/to/three.js-r186.zip
python -m pip install -r requirements-surfaces.txt
python tools/generate_surfaces.py
python build_landscape.py
```

The importer selects four release modules and their MIT license. Examples,
manual assets and font files remain excluded. A hash mismatch requires an
investigation rather than a replacement reference hash.

## Checks

The numerical/runtime suite uses Node 22.16.0, including ES-module loading:

```sh
node --test --test-concurrency=1 tests/core.test.cjs tests/m1.test.cjs \
  tests/landscape.test.cjs tests/ponds.test.cjs \
  tests/weather-column.test.cjs tests/cloud-cache.test.cjs
```

`CLOUD_VALIDATION.json` records actual results and exact build hashes. The new
browser harness treats context-loss messages as failures, including late events.
It supports column scenarios, world/Sun selection, GPU field readback, and a
pixel-identical M1 reset check:

```sh
python -m pip install -r tests/requirements.txt
python tests/capture_clouds.py --regime fair --world moon --sun 45 \
  --probe --reset --output validation/clouds/local-fair.png
```

The harness uses headed Chromium and SwiftShader with an in-memory HTML load
and a harness-only digest bridge where the restricted origin lacks Web Crypto.
Its checks exercise real WebGL shaders. They are software-renderer diagnostics,
not consumer-GPU performance measurements. The `--cloud-quality` option isolates
sky-cache resolution while retaining the specified scene quality; its value is
recorded in the resulting state.

Economy is the automatic software-renderer default. Higher-resolution
Balanced software tests still lost their contexts, including a larger sky cache
with Economy scene quality; those failures are retained.
Completed short scenarios do not establish long-duration stability. Native
WebGPU cloud execution and consumer-GPU speed remain unmeasured.

## Apply the source patch

`revision04.patch` contains ordinary source, tests, documentation and numerical
records, restricted to `immersion/light-cycle/shoreline`. From the repository
root, check it before application:

```sh
git apply --check /path/to/revision04.patch
git apply /path/to/revision04.patch
```

Review the diff and run the checks above. Generated viewers, dependency bundles,
screenshots and source archives remain outside the repository patch. The source
package has no automatic commit/push behavior. `CLOUD_SOURCE_MANIFEST.json`
records baseline identities and payload hashes for future publication.

Earlier revision-03 documentation and evidence are retained as historical
material in the original package/repository. The current checkpoint's evidence
is the separate `CLOUD_VALIDATION.json` and `validation/clouds` directory.
