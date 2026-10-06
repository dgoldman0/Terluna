# Selected sea illustrations

These selected sea images are preserved in Git at the author's
request. Routine rendered outputs and intermediate generation passes remain in
the ignored `../results/illustrations/` or study `research/runs/` directories.
The final images contain no labels; scene identity, prompts, source hashes and
reviews live in these sidecars.


The current Nubium references are the
[2048 × 3072 calculated wide view](nubium-midnight-visibility-view.png) and the
[1024 × 1024 Earth-centred view](nubium-earth-visibility-view.png), with
[display method and checks](nubium-visibility-display.json),
[review](nubium-visibility-review.json) and
[image-generation prompts](nubium-visibility-generation.json).
They incorporate observed Earth structure, the solved atmospheric dimming and
an explicit CIE99 optical approximation. The visibility study supports detectable
interior detail; the display colour and adaptation remain approximations.

The [first new generated pass](nubium-midnight-visibility-first-generated.png)
and [sky-correction pass](nubium-midnight-visibility-generated.png) are retained
as illustration candidates at their native 1024 × 1536. The second pass reduces
sky mottling but both alter the measured cloud pattern. Their planet-luminance
rank correlations with the calculated interior are only 0.46 and 0.32 after
matching canvas size. They are not replacements for the quantitative reference.
The prior white-disk `nubium-midnight-human-view.png` is preserved below for
history and superseded as an appearance reference.

- [IC-1: Ingenii at midday](ic-1-ingenii-noon.png), with
  [prompts and provenance](ic-1-ingenii-noon.json) and
  [review](ic-1-final-review.json).
- [WP-6: Earthlight over western Procellarum](wp-6-procellarum-earthlight.png),
  with [prompts and provenance](wp-6-procellarum-earthlight.json),
  [review](wp-6-final-review.json),
  [composition-reference record](wp-6-calculated-guide.json) and
  [wave-reference record](wp-6-wave-structure-guide.json).

- [Nubium midnight: 2048 × 3072 calculated view](nubium-midnight-human-view.png),
  with the [human-vision method and render record](nubium-midnight-human-view.json),
  [review](nubium-midnight-review.json) and
  [generation prompts and candidate provenance](nubium-midnight.json).
  The [generated candidate](nubium-midnight-generated-candidate.png) remains
  1024 × 1536 and differs from the calculated geometry; it is not the larger
  physical reference. The [first Nubium attempt](nubium-midnight-first-attempt.png)
  is preserved as a rejected photographic treatment, not a naked-eye view.

The earlier Nubium calculated view applies a documented human-vision display operator
to scene radiance. Its glare, colour response and adaptation are approximations;
the featureless bright Earth in that wide view is partly display clipping, not
proof that an observer cannot resolve the planet. The 2048 × 3072 file should be
opened at full size to inspect its larger canvas. A thumbnail shrinks every
feature. The source script can also generate an interactive matched-field size
comparison and an Earth-centred sensitivity view; see the
[workflow](../README.md#nubium-midnight-and-human-vision).

These are illustrations guided by the study, with explicit remaining differences
recorded in their reviews. They are not measured landscapes or validated spectral
renderings. Visible land detail is illustrative; the omission of vegetation is
the optical study's placeholder, not a mature ecosystem design.

## Nobili lake: thin-cloud night

The Nobili view looks west from 0 N, 75.6278 E in the selected flooded atlas.
Earth's geometry, the native terrain silhouette and the night phase are
constrained by the [scene product](../../../research/studies/sea_appearance/results/nobili-night.json).
Cloud weather and lake wave age remain conditional, with their evidence limits
recorded beside the renderings. The generated image is a photographic
illustration; the calculated reference is retained separately.

- `nobili-night-first-generated.png`: first built-in imagegen pass, preserved
  immediately as requested; native output 1448 by 1086.
- `nobili-night-first-prompt.txt`: its complete prompt.
- `nobili-night.png`: selected thin-cloud illustration after three built-in passes; native 1448 by 1086, without artificial upscaling.
- `nobili-night-prompt.txt`: the selected pass's complete prompt.
- `nobili-night-reference.png`: 2560 by 1920 calculated cloud-conditioned guide.
- `nobili-night-clear-reference.png`: calculated clear-sky control.
- `nobili-earth-reference.png`: the separately calculated Earth detail.
- `nobili-night-display.json` and `nobili-night-clear-display.json`: radiance/display records.
- `nobili-night-image-checks.json`: reproducible image-space geometry and texture screens.
- `nobili-night-review.json`: quantitative image checks, visual review,
  source hashes, generation history and remaining limits.

NASA LRO/LOLA and GRAIL supply terrain; NASA/JPL supplies orientation;
NASA/GSFC Blue Marble supplies historical Earth texture. Credits and modelling
notes stay outside the images.

The generated image preserves the broad geometry, while changing fine wave and
Earth texture. Its lower-limb radius fit is 1.4% above nominal, with a small
placement shift. These checks do not validate exact naked-eye appearance. The
cloud field produces a diffuse glow; distinct high clouds require a different
weather sample. Run `python visualization/sea-appearance/nobili_review.py --out
/tmp/nobili-image-checks.json` from the repository root to repeat the selected
image comparison.

## Two cloud illustrations

- [Fecunditatis: nearly full Earth and high clouds](fecunditatis-earth-clouds.png),
  with the [2048 × 2560 calculated reference](fecunditatis-earth-clouds-reference.png)
  and [Earth detail reference](fecunditatis-earth-clouds-earth-reference.png).
- [Eastern Smythii: high clouds at twilight, Earth below the horizon](smythii-twilight-clean.png),
  with the [2560 × 1920 calculated reference](smythii-twilight-clouds-reference.png)
  and [current texture-repair review](smythii-twilight-clean-review.json).
  The [original illustration](smythii-twilight-clouds.png) is preserved.

Both are natural photographic illustrations without text inside the images.
The built-in imagegen tool supplies cloud microtexture and water detail;
calculated references constrain geometry, cloud opacity, illumination and waves.
Native generated sizes are 1122 × 1402 and 1448 × 1086, respectively. They have
not been artificially upscaled. Their first generated images are saved as
`fecunditatis-earth-clouds-first-generated.png` and
`smythii-twilight-clouds-first-generated.png`.

[Initial generation history and full prompt files](cloud-scenes-generation.json),
[initial review and remaining differences](cloud-scenes-review.json),
[reference-copy hashes](cloud-scenes-reference-copies.json), and the per-scene
`*-image-checks.json` files keep the evidence outside the pictures. The two
`*-cloud-opacity.png` guides show upper-cloud opacity, not visible cloud colour.
The clear references and display records are also retained.

These are conditional scenes from the saved equatorial experiment. The cloud
altitudes near 54 and 60 km are not cloud-tower thicknesses or atmospheric scale
heights. Future weather, cross-ring structure, fine cloud forms and exact human
colour appearance are not predicted. The generated images retain residual
composition and brightness differences from the references, itemized in the
review. The Nobili hill-framed night image remains a separate selected scene.

Credits: NASA LRO/LOLA, GRAIL and JPL provide the admitted terrain/datum and dated
Earth geometry; NASA/GSFC Blue Marble supplies historical Earth texture. Source
URLs, rights notes and hashes remain in the study and Earth-input manifests.

## Smythii depth-cue study

The [separate depth-cue illustration](smythii-twilight-depth.png) revises the
Smythii photograph with more legible overlapping veils and gradual recession.
The improvement is modest; the image remains a single view, with unresolved
cloud morphology supplied by image generation. It does not establish a 3-D
cloud reconstruction or supersede the calculated reference. The original
selected illustration is unchanged.

The built-in imagegen tool produced one native 1448 × 1086 edit. That first
pass is saved without further image processing. The [complete prompt](smythii-twilight-depth-prompt.txt),
[review and provenance](smythii-twilight-depth-review.json), and
[image-space comparisons](smythii-twilight-depth-image-checks.json) record the
inputs, measured brightness differences and remaining limits. No text or labels
appear inside the image.

## Smythii texture repair

The [current illustration](smythii-twilight-clean.png) reduces the repeated
cloud curls, embossed edges and harsh water outlines in the earlier images.
Three targeted edits retained too much of that texture; the selected fourth
attempt reconstructs the cloud and water detail directly from the calculated
lighting reference. The broad central shadow, warm flanks and offshore framing
remain. Fine cloud forms and wave detail are illustrative.

The built-in imagegen tool produced a native 1448 × 1086 image, saved without
post-processing or upscaling. The [first repair pass](smythii-twilight-clean-first-generated.png)
is also preserved in Git; the original and depth-cue images remain intact.
The complete prompt set is [initial cleanup](smythii-twilight-clean-prompt.txt),
[cloud repair](smythii-twilight-clean-cloud-prompt.txt),
[water repair](smythii-twilight-clean-water-prompt.txt), and
[selected reconstruction](smythii-twilight-clean-rebuild-prompt.txt).

The [review and provenance](smythii-twilight-clean-review.json) and
[image-space comparisons](smythii-twilight-clean-image-checks.json) record the
selection and its limits. The horizon is about two generated pixels above the
reference; upper-sky, middle-cloud and water median display luminances are about
27%, 14% and 25% lower. Colour also differs. This is a texture repair, not a
new physical solution or a radiometrically equivalent image. No labels appear
inside it.
