# Summit port

Renderings of the central port on the summit, built from the
[summit tower study](../../research/studies/summit_tower/README.md)'s results (`results/summit_tower.json`,
schema `terluna.research.summit-tower/1`):
- **From the study.** The frame's height, its width at every level, a floor band every 200 m, the storey area, the
  zones and their air, and the programme by zone.
- **A design proposal.** Everything else is drawn in so the port can be seen working, at plausible sizes:
  - six leg towers, each a lattice column with a lift-and-stair core;
  - a four-storey wing every 200 m, wrapped round two opposite legs and stepping one leg round per band;
  - berth fingers every 60 m on each wing's outer rim, and roof gardens;
  - gate interchanges at the leg feet, joined by rail;
  - a sealed long-haul terminal at the crown with eight docking arms.

  Ships, trees and people are simple stand-ins. Nothing here is structural design.

| File | Makes |
|---|---|
| [model.py](model.py) | The architectural model of the hexagonal frame in Blender, rendered with Workbench from six viewpoints (the whole tower, inside looking up, a segment, a wing at human scale, a leg foot, the crown), and the screen positions of labelled points |
| [page.py](page.py) | A page of the renders with numbered markers, and the stacking plan |
| [drawings.py](drawings.py) | First-pass SVG sheets: the stacking plan used by the page, and elevations, plans, band details and the crown of both frames |

```sh
blender -b --factory-startup --python visualization/summit-port/model.py -- visualization/summit-port/build/renders
python3 visualization/summit-port/page.py        # build/summit-port-blueprints.html
python3 visualization/summit-port/drawings.py    # build/sheets/*.svg
```

The model renders on the graphics card through OpenGL (Workbench), so it needs a display. Inside the VS Code snap,
run Blender with a clean environment (`env -i PATH=/usr/bin:/bin HOME=$HOME DISPLAY=$DISPLAY XAUTHORITY=$XAUTHORITY
XDG_RUNTIME_DIR=$XDG_RUNTIME_DIR blender ...`), or the snap's libraries break it. Naming views after the output
folder renders only those, and the labels of the others are kept. The model builds in about 20 seconds and uses under
1 GB. Everything written goes to `build/`, which Git ignores.

Open items, as of 2026-09-27:
- **Only the hexagonal frame is modelled.** The round diagrid should carry the same architecture on its six lift
  spines, and be rendered from the same viewpoints, so the two frames can be compared at the scale of a building.
- **The crown's docking arms are not credible structure.** They cantilever 0.9–1.3 km at about 90 times their depth,
  from the tower's narrowest point.
- **The long-haul berths are too few unless ships are very large.** At 30,000–50,000 long-haul passengers an hour,
  2,000-passenger ships need about 15–25 at berth, and 500-passenger ships need 55–95. There is no sky-ship design yet.
- **No faithfulness checks.** Nothing tests the renders yet against the study's numbers.
