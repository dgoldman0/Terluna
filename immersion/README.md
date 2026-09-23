# Immersion

Interactive experiences for exploring how an Open Moon could look and feel.

Each experience lives in its own subfolder, with its source, assumptions, methods
and check records. Perceptual exploration and visual illustration have separate
evidence status from the repository's theoretical, environmental and engineering
feasibility research.

## Experiences

- [Shoreline](light-cycle/shoreline/README.md): the walkable development scene
  (terrain, vegetation, water, clear-sky atlas, column clouds).

Science that grew inside the experience now lives in its domain: the clear-sky
solver in [illumination/sky](../illumination/sky/), the atmospheric column in
[atmosphere/column](../atmosphere/column/), and the A1–A3 light-transport
references in [illumination/references](../illumination/references/). Tools that
render science faithfully (the month-of-light viewer, accuracy labs and the B1
reference renderer) are in [visualization](../visualization/). Past checkpoint
records are in [docs/history](docs/history/).

## Development

Preserve the distinction between calculated clear-sky illumination, reduced
weather/display approximations, and illustrative scenery. Keep assumptions and
known numerical limitations close to the controls and results. Simulation clocks,
weather animation clocks and display exposure represent separate choices.

Shared scientific models remain in their existing subject folders. Immersion
experiments may reference those models, while claims about feasibility or physical
validation require the corresponding research evidence.

Generated distributions, numerical atlases and image assets stay outside Git when
they can be rebuilt from the tracked source. Each experience documents its build
and preserves provenance for distributed snapshots.
