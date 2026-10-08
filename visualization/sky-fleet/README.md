# Sky fleet

Drawings for the [sky-fleet study](../../research/studies/sky_fleet/README.md), built from its product
(`results/sky_fleet.json`, schema `terluna.research.sky-fleet/1`) and the summit tower study's
(`results/summit_tower.json`):

- **The long-haul classes at the crown.** Each of the four liner sizes berthed at the crown in the flight band, from
  one camera from the south-east so the sizes compare directly; each class's berth plan to scale, with its
  dimensions, north, the prevailing wind and a scale bar; the four classes in elevation beside the Hindenburg.
- **How the flyers share the air.** A section to scale through the summit metropolis: the ground, the tower, each
  class's height band, the daytime mixed layer and the storm tops.

From the products: the frame's height and width profile, the flight band's lower edge and the disks' heights, each
class's hull length and diameter, fin span, berths, levels, pitch and arm lengths, the height bands, the ground along
the section (the geography atlas, 7.6 km cells), the mixed layer and storm tops (the cloud-resolving ring). Drawn as
a design proposal at plausible sizes: the round diagrid's members, the sealed disks and their core, the glazed terminal
drum, the arms as deep trusses with enclosed gangways, and each ship's 36-sided hull, fins, window bands, engine pods
and propellers. The hull follows s^0.45 (1 - s)^0.8 along its length, which holds the flight models' volume
(0.475 L D^2), and the fins have the model's area. Nothing here is structural design.

| File | Makes |
|---|---|
| [crown.py](crown.py) | The massing model in Blender (Eevee): one perspective per class (`view_<length>.png`) and the side elevation (`lineup.png`, with `lineup_labels.json`) |
| [sheet.py](sheet.py) | The berth plans (`plans/plan_<length>.png`), a sheet of the four classes with their numbers (`long-haul-classes.png`) and a page of the same with the plans (`long-haul-classes.html`) |
| [figures.py](figures.py) | The section of the air over the metropolis (`heights.png`) |

```sh
blender -b --factory-startup --python visualization/sky-fleet/crown.py -- visualization/sky-fleet/build/renders
python3 visualization/sky-fleet/sheet.py
python3 visualization/sky-fleet/figures.py
```

Eevee renders on the graphics card through OpenGL, so the model needs a display. Inside the VS Code snap, run Blender
with a clean environment (`env -i PATH=/usr/bin:/bin HOME=$HOME DISPLAY=$DISPLAY XAUTHORITY=$XAUTHORITY
XDG_RUNTIME_DIR=$XDG_RUNTIME_DIR blender ...`). This machine's Blender (4.0.2) has no GPU path tracing and no
denoiser, which is why the model uses Eevee. Naming views after the output folder renders only those; `--small`
renders at a quarter of the size. The model builds in a few seconds, each view renders in about 40 seconds, and the
whole set uses under 1 GB. Everything written goes to `build/`, which Git ignores.

**Faithfulness.** The sheet and page print every dimension and number from the product; the plans draw the ships
at their berths' positions, pitch and arm lengths from the product; the renders place the same berths. The background
is a drawing's neutral sky, not the Open Moon's calculated sky.
