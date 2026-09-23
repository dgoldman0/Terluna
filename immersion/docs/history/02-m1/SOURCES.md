# M1 source and dependency record

## Repository and inherited data

Base: `dgoldman0/Terluna`, commit
`4877d8e6c020de51e237deb377955442778917da`,
`immersion/light-cycle/shoreline`.

The seven baseline runtime JavaScript files match their Git blob identities.
A recovered source archive supplies local build assets; earlier packaging
metadata in that archive is superseded by this M1 record. The inherited packed
`data/atmosphere.json` is unchanged. Its own provenance retains original NPZ
hashes and the preceding atmospheric calculation's references.

The earlier source references are preserved in
`docs/checkpoint-sources.md`, and the earlier methods in
`docs/checkpoint-methods.md`. Their historical tests apply to that checkpoint.

## Dependency verified in this work

Three.js r180 / package 0.180.0, provided as `three.cjs` by the user and embedded
without modification. Byte count: 2,006,834. Git blob SHA-1:
`ca4833532c363b72477b2e8a6f47cc0e2fc7b09a`.
The original MIT notice is in `LICENSES.txt` and the viewer.

The bundle's actual MeshStandardMaterial, PMREM and ACES tone-mapping shaders
were executed in the integration and gray-card tests. The implementation stays
on r180. No renderer-version migration or claim about the latest release is
made by this build.

## New implementation inputs

Terrain, fracture-shaped rocks, vegetation, material detail and acoustic code
are procedural project code. New optical coefficients, wind-wave amplitudes,
visibility and surface wetting bands are selected inputs documented in
`METHODS.md`. They have no implied measured lunar validation.

No newly downloaded photographs, scanned materials, commercial models, recorded
sound assets or externally generated illustrative image is used for the runtime
or validation screenshots. Screenshots come from the actual bundled renderer.
