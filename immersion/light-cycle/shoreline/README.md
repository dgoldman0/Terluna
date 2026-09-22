# Open Moon — shoreline r186 / revision 03

This revision extends recovered repository commit
`caa4caa2b2855e9ed8f03597317a05d8224949e6`. It includes the r186 WebGL viewer,
volume-solved connected ponds and wet banks, bathymetry-limited waves, CPU terrain
update reductions, and a separate node-renderer experiment. Native WebGPU remains
an unmeasured experimental path; the lab's actual backend is displayed and exported.

## Open and build

The delivered HTML viewers are self-contained. Open the WebGL viewer directly in
a current browser. In Conditions, choose **After rain**, then **Inspect ponds**
to visit the largest water-filled basin. **M1 noon** restores the reference state.

The source package includes the selected r186 modules and reference material PNGs.
Building needs Python's standard library and performs no network download:

```sh
python build_landscape.py
```

This produces `Open_Moon_Shoreline.html` and `Open_Moon_Renderer_Lab.html`.
The lab offers Automatic, Native WebGPU and TSL/WebGL2 selection, a repeatable
route benchmark and a measurement export. To provide a localhost secure context:

```sh
python -m http.server 8000 --bind 127.0.0.1
# Open http://127.0.0.1:8000/Open_Moon_Renderer_Lab.html
```

The native experiment requires a compatible WebGPU adapter with float32-filterable
textures. A fallback is identified as TSL/WebGL2. Cloud transport, animated rain,
exact leaf lighting and sway remain parity gaps. Water replay in the lab retains
clear-sky lighting; it is separate from a cloud/weather rendering comparison.

## Reproduce dependencies and tests

The archive importer selects four original modules and their MIT license, with
all inputs verified before writing. Examples, unrelated assets and fonts are excluded:

```sh
python tools/import_three.py /path/to/three.js-r186.zip
# Regenerate surface images only when needed:
python -m pip install -r requirements-surfaces.txt
python tools/generate_surfaces.py
```

Preserve the reference surface manifest; investigate a mismatch rather than
accepting changed hashes. Node 22.16.0 was used for the unit suite, including its
ES-module support:

```sh
node --test --test-concurrency=1 tests/core.test.cjs tests/m1.test.cjs tests/landscape.test.cjs tests/ponds.test.cjs
```

`REVISION_03.md` explains the model, calibration change and remaining gates.
`REVISION_VALIDATION.json` indexes measured results. Earlier LANDSCAPE_* documents
and records describe the preceding checkpoint. The source package's small
`build.py` wrapper is a convenience; repository application preserves the existing
source-first repository builder and its delegation to `build_landscape.py`.

## Apply to the recovered repository

The source package includes a dry-run-first helper and an ordinary source patch.
The helper checks the recovered baseline and every payload digest, rejects
conflicting edits or symlinks and leaves unrelated files alone:

```sh
python apply_revision.py /path/to/Terluna
python apply_revision.py /path/to/Terluna --apply --copy-assets
```

The optional asset copy supplies ignored local dependency modules and material
images for immediate building. The helper does not commit, push, change branches,
replace the repository's original builder, or add a source archive to Git.
Review the resulting ordinary-file diff before committing.
