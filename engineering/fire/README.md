# Fire

Fire engineering for buildings under any gravity and air. The models carry Earth's published correlations to
other gravity by Froude scaling, and hold the design-fire data and the people-movement model a fire study needs.
They are first-order tools of the kind used for performance-based design on Earth. They are not validated at
partial gravity, where no fire larger than a small sample has burned.

| File | Holds |
|---|---|
| [smoke.py](smoke.py) | `Air` (gravity, pressure, temperature, composition) and `froude`, the factors that carry an Earth correlation to other gravity and air. Heskestad's flame height with the air's oxygen; NFPA 92's plume; the ceiling jet of a t-squared fire (Heskestad and Delichatsios, NFPA 72 Annex B) and detector and sprinkler response (RTI); a two-zone smoke-filling model, and the exhaust or natural vents that hold smoke above a clear height; MQH's hot layer and flashover threshold; the ventilation-controlled heat release (EN 1991-1-2) with the air supply scaled to gravity, and the supply wind gives through failed openings; stack pressure; the pressure of a water column |
| [design_fires.py](design_fires.py) | EN 1991-1-2 Annex E's occupancies (growth class, peak heat release per m², fire load). Peatross and Beyler's burning rate against oxygen. The measured fall in materials' oxygen limits at lunar gravity (Ferkul and Olson 2011). A slower lunar bracket for a design fire |
| [egress.py](egress.py) | The hydraulic model of people's movement (SFPE: speed against density, 1.3 persons per second per metre of door), the required safe egress time, and the throughput of clearing a floor by air |

What the scaling gives, at lunar gravity in Earth-like air:
- **An equivalent fire.** In a space of the same size, a fire behaves as an Earth fire about 2.5 times larger, played
  about 2.5 times more slowly (more where the air is thinner).
- **Plumes.** They take in about 0.55 times Earth's air for the same fire, so smoke is hotter and less diluted.
- **Flames.** They stand about 1.4 times taller, and more in air with less oxygen.
- **Rooms.** They flash over at a fire about 0.64 times as large. A room's openings feed a fully developed fire about
  0.4 times as much air.

What it does not carry:
- **The fire itself.** How fast a fire grows and how much heat it releases per unit area are inputs.
- **Physics outside the scaling.** Radiation, soot, heat lost to walls, the change to laminar flow at low Grashof
  numbers, and sprinkler sprays do not scale with Froude.

[tests/test_fire.py](../tests/test_fire.py) checks the models in Earth air and in other gravity:
- **Earth air.** In Earth air they reproduce:
  - Heskestad's flame height, 0.235 Q^(2/5) − 1.02 D;
  - NFPA 72's corrected ceiling-jet constants, 0.861, 0.146 and 0.242;
  - MQH's flashover, 610 (h_k A_T A_o √H_o)^½;
  - EN 1991-1-2's ventilation limit, 1.4 MW per m^(5/2).
- **Other gravity.**
  - The ceiling jet and the smoke-filling model follow Froude scaling exactly.
  - Plume entrainment goes as (g ρ²/T)^(1/3).
  - Flame height goes as g^(−1/5).
  - The flashover threshold goes as g^(1/4).
- **Filling.** The zone model's smoke layer comes down later than NFPA 92's first indication of smoke, by less than
  a factor of 2.5.
- **Sensors, people and fires.** A sensing element follows a step exactly. The hydraulic model gives 1.19 m/s and
  1.3 persons/s/m. The design-fire helpers keep the published values.

The [port fire study](../../research/studies/port_fire/README.md) applies them to the summit's central port.
