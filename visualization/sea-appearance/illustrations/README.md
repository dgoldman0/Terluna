# Selected sea illustrations

These selected sea images are preserved in Git at the author's
request. Routine rendered outputs and intermediate generation passes remain in
the ignored `../results/illustrations/` directory. The final images contain no
labels; scene identity, prompts, source hashes and reviews live in these sidecars.


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
