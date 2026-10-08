# Sea scenes

A spectral path tracer for the [sea appearance study](../../../research/studies/sea_appearance/README.md)'s
four coasts. It draws a photographic frame from an eye 2 m above the water at each of the study's 24 moments,
with explicit waves from the hour's spectra, the coast from LOLA heights at the moment's tide, the Sun, the
Earth and the stars. The frames are radiance estimates in CIE XYZ (cd/m²) with their sampling error. The
[scene descriptions](../../../research/studies/sea_appearance/scenes.md) say what each frame holds.

**State:**
- The renderer is built and its parts are tested.
- One 480 by 270 test frame (western Procellarum, the Sun 4° up) ran end to end.
- The full frames, 1920 by 1080, wait for a quieter machine.
- Before them, the display curve needs revising: it rolls highlights off toward white keeping their colour,
  so a blown-out Sun stays an orange dot where a photograph burns it to white.

## How it extends B1

It keeps [B1](../B1_METHODS.md)'s coordinates (a tangent plane over a spherical Moon), its Fresnel water, its
finite solar disk with its own normalisation, heightfield tracing with a grid of height bounds, and sampling
errors kept with every frame. It departs from B1 in four ways:

- **The sky and the air.** They come from the solved spherical sky that the study's results use
  ([illumination/sky/aerial.py](../../../illumination/sky/aerial.py)): the sky's radiance in every
  direction, and for each surface the light the air adds in front of it and the transmittance to it.
  B1 traces every path through the air by Monte Carlo instead, with three limits for this work:
  - it has no earthlight;
  - it cannot reach the deep twilight with usable noise;
  - it uses a different atmosphere from the study's.
- **The sea is explicit.** Waves from the hour's spectra ([realization.py](../../../illumination/water_surface/realization.py))
  run on four tiles from 2,048 m to 1.43 m:
  - the largest tile is traced as geometry over the curved sea;
  - every tile's slopes are filtered over each pixel's footprint (LEAN mapping);
  - the slopes left inside the footprint form a Gaussian distribution of facets, reflecting as the study's
    model does.
- **The water has its own light**, from the study's water-column optics, where B1 absorbed everything that
  entered it.
- **The scene is real.** LOLA terrain at the tide's level, with a placeholder bare-soil colour, replaces
  B1's authored coast. The Earth is a second light source, drawn with its phase, and the bright stars are
  a separate layer.

## Files

- [kernels.py](kernels.py): numba kernels.
  - Heightfield tracing over the curved Moon, and filtered slope statistics.
  - Visible-facet sampling for Gaussian slopes with a mean (Heitz and d'Eon 2014).
  - The Sun's glitter by the Cox–Munk formula, and the Earth's by multiple importance sampling.
  - The air in front of each surface, and the render loop.
- [scene.py](scene.py): a frame's inputs from the domains and the study. The cameras, views and wave spectra
  come from [scenes.py](../../../research/studies/sea_appearance/scenes.py).
- [stars.py](stars.py): the star layer, from [illumination/stars/sky.py](../../../illumination/stars/sky.py).
- [render.py](render.py): the driver. It renders blocks of rows, saves each as it completes, and resumes
  after an interruption.
- [test_kernels.py](test_kernels.py): seven checks.
  - Sky reflection within 1% of the study's quadrature, down to 0.5° grazing.
  - Visible-facet sampling with a tilted mean.
  - Both samplings of the Earth's glitter against the disk formula.
  - Heightfield tracing against brute-force marching.

```sh
NUMBA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 <numba environment>/bin/python \
    visualization/reference-renderer/seas/render.py --coasts "W Procellarum" --moments "low Sun"
<numba environment>/bin/python -m pytest visualization/reference-renderer/seas/test_kernels.py
```

Frames go to `research/runs/sea_appearance/renderings/`. They need:
- the coast grids from `python -m geography.coast_terrain`;
- the wave runs on the research drive;
- the solved sky in `research/runs/optical_comfort/spherical/`.
