# Flight

Concept-level models of flyers on any world, from buoyant ships to wings and rotors. They size a flyer's lift,
structure and power for a gravity and an air, at the first order used to compare concepts, sizes and worlds.

| File | Holds |
|---|---|
| [buoyant.py](buoyant.py) | Lifting gases and the net lift per cubic metre (purity, superheat); a streamlined hull's volume and area; Lamb's apparent-mass coefficients and Munk's moment; the midship bending moments from a gust (Woodward's envelope) and from a deflated gas cell; `rigid`, a Zeppelin-type hull sized by the groups of the airship weight statements (longitudinals, frames, empennage, cover, cells, power plant, fixed equipment), with its structure split into weight-driven, aerodynamic and areal shares; `pressure_hull`, a non-rigid envelope's pressure, crown tension and fabric use; the power to fly at a speed |
| [winged.py](winged.py) | The speed a wing needs to carry its mass; Pratt's gust load factor; the ideal spar caps of an elliptically loaded, straight-tapered wing and a wing mass with a non-optimum factor and an areal term; the power to fly at a lift-to-drag ratio and to hover |

## The rule of similarity

Gravity cancels in buoyant lift and returns in the loads. Every load on a flyer is one of three kinds:

- **Weight-driven.** The flyer's own weight, and bending from an uneven spread of lift and weight. The share of the
  flyer's mass its structure needs for these grows as g L ρ_m / σ.
- **Aerodynamic.** Gusts and manoeuvres. At a given speed, gust and air density, the share they need is the same at
  any size and gravity for a buoyant hull, and for a wing of a given loading in pascals.
- **Areal.** Covers, gas cells and minimum gauges. Their share falls as the flyer grows.

So under a sixth of Earth's gravity, a flyer scaled up by g_earth / g_moon = 6.04 in every length, flown at the same
speed in air of the same density, carries the same stresses from its weight and the same from the same gusts, with
the same Froude number: the same flyer, six times longer and 220 times heavier. The tests check that both models
keep this rule.

## Assumptions

- **Rigid hulls.** The longitudinals carry the ultimate midship moment as a thin-walled tube of the hull's radius
  and average 0.6 of that section along the hull. The gust moment is Woodward's envelope, 0.10 (U/V) q Vol^(2/3) L
  (1976), and the static one 0.016 of weight times length with a cell deflated. The frames carry the ring force of
  the cells' pressure head and a share of the hull's bending. Each group has its own factor on its ideal amount; the
  [sky-ship study](../../research/studies/sky_ships/README.md) fits them to four rigid airships' weight statements.
- **Pressure hulls.** The envelope's pressure is the larger of a margin on the dynamic pressure and what keeps its
  longitudinal tension above the bending stress at limit load (FAA-P-8110-2, §4.43(a)); the crown carries the gas
  head as well; the fabric carries the hoop tension four times over (§4.43(b)).
- **Wings.** Elliptic loading on a straight-tapered, unswept planform; the spar caps' depth a fixed share of the
  thickness; a non-optimum factor and an areal term, which the study fits to ten transports' wing groups.
- **Left out.** Joints, fatigue, flutter, dynamic gust response, ground handling, landing gear, sweep and the
  details of fins and control surfaces, except as the fits absorb them.

[tests/test_flight.py](../tests/test_flight.py) checks:
- buoyant lift without gravity in it, for pure and impure gases and with superheat;
- Lamb's coefficients for spheroids of fineness 4, 6 and 10, and the sphere's one half;
- Woodward's 3,950,000 lb ft for the Shenandoah at 91 ft/s in a 35 ft/s gust;
- the rule of similarity for rigid hulls and for wings, to nine places;
- the crown's tension in a pressure hull; the elliptic loading's planform factor (1/32 for a rectangular wing);
- Pratt's alleviation factor and the gust's extra lift, which gravity leaves alone; hover power as gravity to 1.5.
