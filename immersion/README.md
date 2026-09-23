# Immersion

Interactive experiences for exploring how an Open Moon could look and feel.

Each experience lives in its own subfolder, with its source, assumptions, methods
and check records. Perceptual exploration and visual illustration have separate
evidence status from the repository's theoretical, environmental and engineering
feasibility research.

## Experiences

- [Shoreline](light-cycle/shoreline/README.md): a walkable development scene on an
  authored cove, with the clear-sky atlas, terrain, vegetation and water.

The fixed-viewpoint "month of light" sky viewer is scientific visualization of
the illumination domain's clear-sky atlases; it now lives in
[visualization/month-of-light](../visualization/month-of-light/), and its solver
in [illumination/sky](../illumination/sky/).

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
