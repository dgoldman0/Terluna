# The summit port and metropolis

The explorable map of the Moon's central sky port (a 24 km tower on the far side's summit) and the metropolis round
it, for the Unreal game. The research behind it: [the summit port and metropolis decisions](../../../research/decisions.md#the-summit-port-and-metropolis),
the [metropolis brief](../../../habitation/summit_metropolis/README.md) and the [tower study](../../../research/studies/summit_tower/README.md),
whose hand sizing (`form.py`) sets every level and member size.

## What is here

- `ring0/design.py`: ring 0's design as explicit data. Ring 0 stands on the transfer ring, 3.06 km above the summit,
  150 m wide and 20.5 km round. It has six gates where the legs arrive, six quarters, and Quarter 1 (the Market
  Quarter) laid out place by place at hall level and on its roof park.
- `ring0/draw.py`: its drawings, in `ring0/drawings/`:
  - a key plan;
  - Quarter 1's hall-level and roof plans, unrolled along the ring;
  - a cross-section through the Grand Market.

  Rerun with `python3 immersion/world/summit/ring0/draw.py`. The SVGs are the record; the PNG previews are not
  committed.

Kinds, in this lane's terms:
- **Informed:** the tower's structure, its levels and sizes, and the climate at each level.
- **Artistic:** the architecture and landscape, as a design proposal. The trees and plants are placeholders, since
  the Open Moon's plants are not yet designed.

## Status

This is the zoning stage of ring 0's design: where each place goes, for review before modelling. Each building and
garden is then modelled part by part in Blender, and its model and modelling steps are committed here.
