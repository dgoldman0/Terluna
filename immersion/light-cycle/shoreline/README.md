# Open Moon — shoreline landscape/ecology checkpoint

This checkpoint recovers the reviewed landscape source based on M1 commit
`4c7affae780d3ff7ef5581e88aac7ae32ada7505`. All twelve runtime modules and the
landscape page template retain their exact review-package bytes. The renderer
remains Three.js r180 / WebGL. Visual acceptance and continuous SwiftShader
traversal remain open.

## Build and test

The repository stores ordinary editable source, the material generator and its
reference manifest. Generate the two surface PNGs locally before the first build:

```sh
python -m pip install -r requirements-surfaces.txt
python tools/generate_surfaces.py
git diff --exit-code -- assets/surfaces/manifest.json
python build.py
```

The manifest must remain unchanged. A mismatch requires investigation; preserve
the reference hashes. Both generated PNGs were reproduced byte-for-byte during
recovery. They are ignored by Git, along with the generated viewer and optional
vendor bundle. No source ZIP or encoded source archive replaces the editable files.

For the exact self-contained review build and the complete unit suite, place the
reviewed r180 bundle at `vendor/three.cjs`. Its Git blob identity is
`ca4833532c363b72477b2e8a6f47cc0e2fc7b09a`; the builder checks it. The saved landscape
source package includes that bundle. Without it, the viewer retains the existing
pinned, hash-checked first-open dependency loader.

```sh
python build.py
node --test --test-concurrency=1 tests/core.test.cjs tests/m1.test.cjs tests/landscape.test.cjs
```

With the pinned vendor and reference maps, the generated HTML has SHA-256
`2bfe13b25a155cf24383c39a2392320bf7dad6d6a697651164515f77943cc3be` and 14,218,498 bytes.
The original repository `build.py` is retained with a small delegation to
`build_landscape.py`; the M1 template and inherited atmosphere stay unchanged.

## Model and evidence

`landscape.js` owns shared elevation, drainage, substrate and habitat fields.
`terrain.js` supplies moving, stitched terrain grids. `ecology.js` supplies stable
plant placement and connected branching. `surface-water.js` supplies spatial
water stores; `materials.js` consumes the shared fields and coordinated maps.
`LANDSCAPE_METHODS.md` defines the assumptions and model boundaries.

`LANDSCAPE_VALIDATION.json` is the unchanged original review record. Its
`remote_changes_performed: false` describes that review's time, before recovery.
Its individual run and image paths refer to the saved source/review package.
The failed walking log and verified optical log are also tracked here; other raw
review logs and screenshots remain in that package. `validation/recovery.json`
records this recovery's new build, regeneration and 87-test checks separately.
The browser harness and both GPU probes are ordinary source under `tests/`.

Continuous SwiftShader walking loses its graphics context even in Economy.
Consumer-GPU performance remains unmeasured. Materials and foliage remain
procedural, terrain updates can hitch, and visible connected puddles, material
height displacement and resolved shoreline run-up still need work.

## Next revision

See `RENDERER_ROADMAP.md` and `LANDSCAPE_HANDOFF.md`: first migrate the unchanged
WebGL scene to r186, then improve model-to-render integration and evaluate a
separate WebGPU/TSL implementation against the same state. The recovery commit
performs none of that renderer migration.
