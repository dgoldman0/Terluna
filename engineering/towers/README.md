# Towers

Concept-level models for very tall structures, on the Moon or on Earth. They size
a frame and screen devices; they are not structural design.

| File | Holds |
|---|---|
| [lattice.py](lattice.py) | A tapered square lattice tower of four corner legs and braced faces, with enclosed floor bands spaced up it, sized from the top down so each leg carries the weight and the wind's moment above it at an allowable stress. It returns the frame's mass, its sway under a wind, its first natural frequency (Rayleigh, static-deflection shape), a whole-frame buckling factor (Rayleigh–Ritz, same shape) and the footing load, for any gravity and air-density profile, and picks the lightest frame that meets the sway, buckling and base-width checks |
| [wind_energy.py](wind_energy.py) | Weibull wind distributions fitted to percentiles, the wind's power, the mean output, thrust and parked drag of rotor and bladeless devices, and Band's single-crossing collision probability for a flier passing through a rotor |

The tower model's assumptions:

- **Loads.** Member buckling is one factor on the strength (0.75). One safety factor (1.6) covers gravity and wind together. The wind is a peak gust, taken as uniform up the tower, with a dynamic factor of 1.2 on its pressure. The lattice's force coefficient is TIA-222-G's square-tower formula, Cf = 4.0e² − 5.9e + 4.0, times 0.57 − 0.14e + 0.86e² − 0.24e³ (at most 1) for round members in subcritical flow, and times 1 + 0.75e (at most 1.2) for wind along a diagonal. The floor bands' facades have a drag coefficient of 1.4.
- **Bracing.** It is 35% of the legs' mass, or what the wind's shear needs if that is more.
- **Left out.** Joints, fatigue, member-level vortex shedding, construction stages, ice and the ground's own mechanics.

[tests/test_towers.py](../tests/test_towers.py) checks the model against textbook cases:
- a column of constant stress under a top load grows as the exponential of height over its self-weight length;
- a uniform cantilever's sway under a uniform load is wH⁴/8EI;
- Rayleigh gives its first frequency within 0.4% of the exact value;
- Rayleigh–Ritz gives Greenhill's self-weight buckling load within 2%.

The [summit tower study](../../research/studies/summit_tower/README.md) applies both models to a sky tower on the Moon's highest ground.
