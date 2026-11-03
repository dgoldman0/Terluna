# Dance renderer

Renders a dance motion (the segment poses of `humanities/dance/motion_io.py`) as a skinned dancer on a lit stage,
in Blender 4.0.2 with EEVEE (legacy). The body mesh is Blender Studio's CC0 realistic female base mesh (Human Base
Meshes v1.4.1); the skeleton is `humanities/dance/body/body.json`, which the physics owns.

## Running it

From the worktree root (any Python 3 works; the launcher only calls Blender):

    python visualization/dance/run_render.py MOTION.npz --out OUT.mp4 \
        [--garment CLOTH.npz ...] [--preset draft|test|detail|video|final] [--stills 2.0,4.4,6.5] [--stills-preset final] [--stills-only] \
        [--frames 1:120] [--cam-az DEG] [--cam-dist M] [--cam-aim M] [--cam-fstop N] [--rebuild] [--save-scene SCENE.blend] \
        [--bodice R,G,B] [--bodice-look matte|satin] [--yoke R,G,B,A] [--sheer hashed|blend] [--samples N] \
        [--bodice-shell STYLE[,NECKLINE]] [--trim R,G,B] [--light stage|soft]

- Presets: `draft` is 854x480 with 24 samples and sampled soft shadows, the default for whole videos until a complete
  video is near (from the per-pixel cost below, about 45% of the 720p frame time; not yet timed); check stills render
  at `final` by default (`--stills-preset`), so detail is judged at full size while videos stay small, and the
  1080p presets are kept for a final draft. `test` is 1280x720 with 24 samples and sampled soft shadows; `video` is
  1920x1080 with 32 samples and filtered shadows; `final` is 1920x1080 with 80 samples and soft shadows, meant for
  stills. `--samples` overrides a preset's sample count. Motion blur is off because EEVEE's motion blur and an 8-pixel
  volumetric tile both returned black frames in background renders with Blender 4.0.2, on the integrated GPU (below).
- The GPU. On a laptop with hybrid graphics Blender renders through Mesa on the integrated GPU unless OpenGL is
  offloaded to the NVIDIA GPU, which run_render.py does whenever the NVIDIA driver is loaded
  (`__NV_PRIME_RENDER_OFFLOAD=1`, `__GLX_VENDOR_LIBRARY_NAME=nvidia`). Until 2 November every render here ran on the
  integrated GPU (an Intel UHD GT1): stills at 1080p or at 80 samples came back black or striped, frames took 5 to 10 s,
  and blended transparency crashed in the Mesa driver. On the NVIDIA GPU a 1080p still at 80 samples renders cleanly
  in about 6 s after 8 s of scene setup, and 512 samples at 720p in about 20 s in all. The frame timings in this
  README before that date are the integrated GPU's.
- Sheer cloth is drawn with hashed transparency (`--sheer hashed`, the default), which is dithered and smooths out only
  with many samples: at 24 samples it read as grain; at 200-256 samples it is smooth (a measure of pixel noise in the
  skirt fell from 0.032 at 24 samples to 0.023 at 200 and 0.022 at 512, the rest being the folds' own detail).
  `--sheer blend` alpha-blends instead, smooth at any sample count, but layers sort by object, the folds lose their
  shading and the body shows through every layer, so the design stills use hashed transparency at 256 samples.
- `--stills` writes PNGs beside the video as `<stem>_tSS.SS.png`; `--stills-only` skips the video. `detail` is the
  final quality at 720p.
- The camera stands at a fixed bearing (`--cam-az`) and distance (`--cam-dist`) and follows the smoothed centre of
  mass. `--cam-aim` sets how far above the centre of mass it aims (0.06 m by default) and raises the camera with its
  aim, so a close still is taken level; `--cam-fstop` sets the aperture (2.8 by default). The bodice stills of
  3 November used `--cam-aim 0.36 --cam-dist 2.1 --cam-fstop 8`, which frames the standing dancer from the hips to the
  head with the whole torso in focus.
- A garment file may carry its own render settings (meta `render`: look, sRGB colour, alpha), which the design
  stills use for provisional colours; cloth gets one level of subdivision smoothing at render time. Opaque cloth
  (satin, silk, wool, velvet, jersey) is given 1.5 mm of thickness with a Solidify modifier. A sheer layer stays a
  single surface seen from both sides: solidified, it was two veils 1.5 mm apart, and stacked sheer layers read as
  opaque white. Satin has a silk satin's soft sheen (roughness 0.38, specular 0.45, anisotropy 0.3); until
  2 November it was glossier (0.22, 0.7, 0.6) and read as rubber or plastic.
- Design stills can dress the leotard zone and the skin above it. `--bodice r,g,b` makes the leotard a bodice of that
  sRGB colour, in the look `--bodice-look` names: `matte` (the default), a matte silk with a soft sheen in its own
  colour that sits with chiffon skirts, or `satin`, the glossy satin of the first design stills, which read as a
  rubber catsuit beside chiffon. `--yoke r,g,b,a` lays a sheer chiffon yoke of that sRGB colour over the skin from
  the bodice's neckline to a jewel neckline at the base of the neck, front and back, and over the shoulders to the
  shoulder point. Its strength `a` (0 to 1) is the share of the skin the chiffon veils seen face on. The veil
  thickens toward the silhouette, where light crosses the sheet obliquely (transmission `(1 - a) ** (1 / cos)`), and
  doubles in the hems at its neckline and armholes. Ivory chiffon over this dancer's skin reads as a pale grey veil
  at 0.5 and shows the skin clearly at 0.3. Without `--yoke` the skin is bare as before. The build stores the
  yoke's zone as the vertex fields `yoke` and `yoke_hem`, and the render adds the chiffon to the body material
  (`set_yoke`), so its look can change without a rebuild; a dancer built before 2 November lacks the zone and needs
  `--rebuild`.
- `--bodice-shell STYLE[,NECKLINE]` builds the bodice as a garment over the torso instead of a colour on the leotard
  zone, which follows every contour of the skin and reads as moulded on the body. The shell is a copy of the body's
  torso from 1.00 m up to just above the neckline (arm, head and neck weights under 0.3), skinned by the body's own
  armature, so it moves as the body does and costs no simulation. It spans the body's hollows as boning and a lined
  shell do: Taubin's smoothing (steps of 0.5 and -0.53) takes out relief a few centimetres across without shrinking
  the torso or the bust, the shell lies at most 1.5 mm inside that smoothed torso and up to 15 mm outside the body,
  an outward-only fill then closes any hollow left, such as the cleavage and the spine's groove, and the whole is
  lifted 3 mm and given 1.5 mm of thickness and one level of subdivision. The copy's custom split normals are
  cleared: the first shells kept them and shaded as the body itself, the abdominal muscles showing whatever the
  shell's shape. NECKLINE is `sweetheart` (the default: 1.212 m between the breasts, arching to 1.252 m over each),
  `vee` (from 1.165 m at the centre front) or `straight` (level at 1.235 m), each 1.205 m across the back; the shader
  cuts it by alpha clipping from a per-vertex field and pipes it with 8 mm of satin in the `--trim` colour (gold by
  default), and the body's leotard is cleared above it so the skin or the yoke shows there. STYLE sets the surface
  as a bump on the rest position: `structured`, princess seams curving from the bust to the waist and down the back;
  `ruched`, soft gathers rising from the centre front's waist to the sides, about 17 mm apart, their spacing and
  depth varied by two noise fields (a regular 11 mm sine read as stripes); `pleated`, knife pleats fanning from 7 cm
  below the waist. The colour is `--bodice`'s (ivory without it). The shell does not crease or slide as cloth would.
- `--garment` takes a cloth sequence: `faces` (F, 3) int and `verts` (N, V, 3) in the dance world frame at the motion's
  rate (or at the file's own `fps`). The renderer resamples it to 30 fps, writes a PC2 cache into `pc2/` beside the
  output and plays it with a Mesh Cache modifier. `test_sheet.py MOTION.npz OUT.npz` makes a kinematic test skirt for
  checking this path.
- The rigged dancer is built once into `/media/projectspace/terluna-research/humanities/dance-runs/renders/dancer_rigged.blend`
  (about 40 s) and reused; `--rebuild` remakes it and rewrites `landmark_check.json` (and `proposal_landmarks.json`
  when a landmark differs from body.json by more than 2 cm) in the render agent's folder.
- Blender runs with `--threads 3`; EEVEE renders on the GPU (about 0.3 GB). The frame time is mostly GPU work per pixel
  and sample (subsurface, reflections, haze, shadow maps), so 1080p costs about 2.2 times 720p. A virtual-memory limit
  (`ulimit -v`) fails under the NVIDIA driver, which reserves far more address space than it uses (a 14 GB limit
  stopped a three-still render with std::bad_alloc while Blender held about 0.3 GB); cap real memory instead, above
  the scene's peak and without swap, for example with `systemd-run --user --scope -p MemoryMax=10G -p MemorySwapMax=0`.
  With the fourth costume round's fourteen garment files (each tier's bound hem a file of its own) Blender's resident
  memory peaked at 7.7-8.5 GB for a still at 720p or 1080p, and under an 8 GB cap the memory controller stopped the
  canopy still; what in the scene takes that much is not yet measured.

## Pipeline

`render.py build` (Blender side):

1. Imports the collection "Body Female - Realistic" (body and both eyes), turns it to face +x with its left on +y,
   centres it and puts the soles on z = 0.
2. Measures the mesh's joint landmarks and compares them with body.json (rules in `measure_landmarks`).
3. Builds a fit armature at the mesh's own joints (A-pose), takes Blender's heat weights on the 10.6 k-vertex base mesh,
   and cleans them by rule (below), then applies two multires levels (169 k vertices; the subdivision interpolates the
   weights).
4. Reposes the mesh into body.json's rest pose: each limb bone is rotated onto body.json's rest direction, stretched
   along its axis to body.json's length and placed on its joint centre; the trunk, neck and head keep the mesh's shape.
5. Builds the render armature from body.json's joint centres alone (21 parentless bones, one per segment) and binds the
   reposed mesh to it with the same vertex groups.
6. Adds the costume zones (leotard, slippers, yoke), the face and hair, the eyes and the materials, and saves the
   rigged dancer.

`render.py render` opens the rigged dancer, keys every bone from the motion, builds the stage, lights and camera, adds
garments and renders.

Driving: each bone follows its segment's world pose. With the bone's rest matrix `[R_r | h]` and the segment's pose
`(q, p)`, the bone's basis is `rest^-1 [R(q) R_r | p]`, keyed every render frame (30 fps, linear interpolation between
the motion's 60 fps samples, normalised lerp for the orientations). The root moves with the pelvis segment.

## Designed by hand and generated by rule

Designed by hand (choices in the code, judged on renders):
- the costume: a camisole-cut leotard in deep wine as a smooth vertex-field zone of the body material (high-cut leg
  line, scoop neck and back, sleeveless), 1.8 mm proud of the skin, with the fabric bridging hollows (cleavage, navel,
  spine groove); pale pink satin slippers below 5.5 cm; a yoke zone for the sheer yoke of the design stills,
  running from under the leotard to a jewel neckline at the base of the neck and over the shoulders to the shoulder
  point, its armholes on the leotard's armhole line (summed arm weight 0.32), with a hem band 4 to 6 mm wide inside
  its neckline and armholes. The neckline follows the neck, taken as an elliptic cylinder of semi-axes 6.0 cm front
  to back and 7.0 cm across about an axis leaning forward with it (x = 0.012 + 0.4 (z - 1.36) m), so it meets the skin
  at the jugular notch (z 1.337 m), at the sides of the neck's base (1.39 m) and just below the seventh vertebra
  (1.39 m);
- the hair: a dark base cap 1.5 to 6 mm off the scalp that thins out at the hairline, about 13.7 k hair-curve strands
  combed over it to a classical bun at the crown, 2.6 k strands wound round the bun, all rigid with the head bone;
- the face: lips (the base mesh's lip face sets) and brows (an arch above each eye) as colour fields;
- materials: subsurface skin (warm medium brown, red-weighted scattering radius), knit leotard, satin slippers, hair,
  eyes with iris and pupil from the angle to the eye's axis; at render time, the bodice's matte silk and satin looks
  and the yoke's chiffon (a BSDF in 0.35 of the fabric's colour plus 0.25 of the skin's albedo, for light passed
  through the threads to the skin and back, with a soft sheen, mixed over the skin). The yoke has no weave bump: a
  0.5 mm weave, below a pixel at these framings, made 80-sample stills come back black on the integrated GPU;
- the stage: matte black Marley (roughness 0.52 to 0.68), a curved cyclorama with an emitted dusk wash (amber horizon
  to deep indigo), a black ground row at its foot, light haze;
- the lighting (`--light soft`, the default since 3 November): booms at shin, mid and head height on both sides (a pale
  warm tint stage left, a pale cool one stage right) at 0.8 of the stage rig's power, the shin-level pair at 0.4 so the
  legs light like the rest of her, a soft front key of 950 W, a soft fill from the camera's side, a high back rim of
  2.1 kW and a top pool that follows the dancer, graded with AgX's base contrast. `--light stage`, the earlier rig, has
  the booms in warm amber stage left and steel blue stage right (2.2-2.7 kW each), a key several times dimmer than a
  boom, a 4.2 kW rim and AgX's medium-high contrast; under it the skin turns orange or blue-grey with the side it faces
  and the face goes dark when it turns from the booms to the camera, which reads as jarring changes of skin tone between
  frames;
- the camera: 50 mm at f/2.8 following the smoothed centre of mass (across the floor with sigma 0.8 s, in height with
  sigma 0.35 s, aimed 6 cm above it). Phrases in place: 5.5 m away at a 30 deg bearing from the dancer's front, 1.1 m
  high. Travelling phrases (more than 1 m of travel): 6.8 m away, square to the travel so the line reads in profile.
  The lighting rig turns with the camera bearing and travels with the dancer; the cyclorama is centred on the path.

Generated by rule (repeatable from the scripts):
- landmarks and the landmark check;
- the weights: Blender's heat weights, then left/right masks at the midline, arm groups removed from the legs and lower
  trunk, the upper arm allowed onto the torso only at the deltoid cap, thigh, shank, foot, toes and clavicle bounded in
  height, then Blender's own smooth (4 passes), clean, limit to four and normalise. The head and neck are weighted by
  height up the neck, from its base (the notch in front, the seventh vertebra behind) to the jaw line and the skull's
  base: thorax to neck over the lower half, neck to head over the upper half, and everything above the jaw line (face,
  jaw, mouth) rigid with the head. The artist's head set also covers the neck, and the jaw and mouth are a set of their
  own; until 2 November the rules held the neck rigid with the head and left the jaw half on the neck, so the mouth
  opened and the throat creased at its base whenever the head moved on the trunk;
- the skinning stack: leotard thickness (Displace), linear-blend Armature, a dual-quaternion Armature mixed in over the
  hips and lower trunk (`dq_mix`), Corrective Smooth (length-weighted, 8 iterations, masked off the head, hands and toes);
- the hair strands' paths (each step toward the bun in the scalp's tangent plane, snapped back to the cap);
- the camera path (Gaussian-smoothed centre of mass, sigma 0.8 s) and the top pool's path.

## Files

- `render.py`: the Blender side (build and render).
- `run_render.py`: the launcher.
- `test_sheet.py`: the kinematic test skirt for the garment path.
- `dress_form.py`: measures the built dancer's skin against body.json's capsules and writes `dress_form.json`
  (schema `terluna.visualization.dance-dress-form/1`), the collision body for dressed renders' cloth. body.json's
  capsules lie inside the rendered skin by up to 5 cm at the thighs, and the rendered waist sits 5 cm forward of the
  abdomen capsule, so cloth simulated on them passed through the skin and the waistband sank into the waist. The
  form splits each limb capsule into three tapered pieces reaching the skin (95th percentile of its distance from the
  axis near each node, over skin within the capsule's span; the hands keep their capsules inside the sleeves' cuffs), turns the trunk's upright capsules
  into left and right pairs spanning its cross-section, and gives the waist's outline for the waistband. It is a
  stand-in until the body model and the rendered mesh are fitted as one body (question 1). Rebuild it after
  `--rebuild`: `blender -b DANCER.blend --python visualization/dance/dress_form.py -- --out visualization/dance/dress_form.json`.

The rigged dancer, renders, assembled scenes and PC2 caches live on the drive under `dance-runs/renders/`; the landmark
check and its proposal live in `agents/render/`.
