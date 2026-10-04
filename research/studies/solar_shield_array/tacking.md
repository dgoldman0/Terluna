# Photogravitational tacking and a cycling shield fleet

The author proposed photogravitational tacking on 4 October 2026 as a route
to lower holding costs at the existing mass per area. The first
[numerical cycling study](cycling.md) now tests that direction at **50 g/m²**:
free orbital motion, a physically constrained ideal specular sail, lunar and
Earth eclipses, returning arcs and multi-year service-and-coast propagation.
It establishes candidate tile motion and measured shadow-service fractions.
Fleet coverage and handovers remain unsolved.

The sections below record the original literature review and inventory screen.
Their assumed duty fractions are superseded for the computed candidates by
[results/cycling.json](results/cycling.json); they remain useful algebraic comparisons.

## Why the existing result leaves this route open

The [held-aperture study](report.md) computes the force required to follow
chosen trajectories and minimizes the remaining electrical thrust. It
already permits an optimistic mixture of passive outgoing photon directions
at every instant. Alternating optical settings on one of those same paths
cannot beat its instantaneous momentum envelope. The chosen moving path's
0.375-mm/s² sunward floor belongs to that path.

Orbital sailing changes the path as well as the optical setting. In an ideal
central solar gravity field, sail acceleration changes specific orbital
energy and angular momentum according to

\[
\dot{\mathcal E}=\mathbf v\cdot\mathbf a_{\rm sail},\qquad
\dot{\mathbf h}=\mathbf r\times\mathbf a_{\rm sail}.
\]

A transverse optical force can oppose orbital velocity and lower orbital
energy even while its radial component points away from the Sun. Gravity
curves the resulting motion. In the actual Sun–Earth–Moon problem, the moving
primaries also exchange energy and angular momentum with the spacecraft;
the central-field equations above are explanatory, not the cislunar model.
The return trajectory and control law must be solved with the ephemeris.

The previous search used three synodic harmonics, 20,000–180,000 km sunward
distance and at most 2,000 km of transverse centre motion. Its aperture cost
enlarges a complete disk when its centre moves sideways. It did not optimize
free orbital states, long return arcs, many phased trajectories or exchanging
tiles in the aperture. Its 238-TW result therefore does not establish a
minimum for an orbiting fleet.

## Relevant primary work

The literature usually describes this work through solar-sail periodic or
quasi-periodic orbits, resonances, and attitude-controlled station keeping.

| Source | What it establishes | Relation to this study |
|---|---|---|
| Heiligers, Macdonald & Parker (2016), [Extension of Earth–Moon libration point orbits with solar sail propulsion](https://doi.org/10.1007/s10509-016-2783-3) | Differential correction and continuation produce sail-assisted Lyapunov, halo, vertical, Earth-centred and distant retrograde orbit families. Stability varies with the family and sail acceleration; stronger sailing can destabilize a previously stable orbit. Earth/Moon shadowing is omitted. | Supplies candidate orbit families and a continuation method. Neither the orbit families nor their stability analysis impose continuous lunar UV coverage. |
| Gao et al. (2023), [Station-keeping for Earth–Moon solar-sail resonant libration point orbits](https://doi.org/10.1016/j.cnsns.2023.107274), [author-hosted text](https://www.maia.ub.es/~gerard/2023/CNSNS-2023.pdf) | Sail attitude feedback reshapes unstable and central modes in the Earth–Moon quasi-bicircular problem. Long-term numerical control and orbit-determination errors are examined. | Shows how optical control can maintain an orbit by correcting its dynamics. The model is a quasi-bicircular approximation; coverage and actual filter optics need separate constraints. |
| Chujo (2024), [Quasi-periodic orbits of small solar sails with time-varying attitude around Earth–Moon libration points](https://doi.org/10.1007/s42064-023-0186-0) | Time-varying attitude produces orbit families that are transferred from linear and nonlinear circular models to DE435 ephemeris dynamics by multiple shooting. An L2 example uses A/m = 2 m²/kg and remains nearby for about 2,000 days. | Demonstrates an appropriate progression to an ephemeris solution. Earth visibility around L2 is a different geometrical constraint from maintaining a solar shadow on the Moon. |
| Estrada et al. (2025), [Solar Sail Orbits in Cislunar Space for Lunar Surface Surveying](https://ntrs.nasa.gov/citations/20250006966) | A 0.056-mm/s² sail maintains an L1 halo-like trajectory for 12 revolutions in a bicircular model; the attitude schedule includes long coasting intervals. Eclipses are omitted. A proposed returning excursion could only be solved with an unphysical sunward sail force. | Supports alternating sailing and coasting, and makes full-cycle force admissibility essential: the outward trip alone cannot establish a recurring service orbit. |
| Bellinazzi, Matonti & Romano (2026), [Instantaneous Photo-Gravitational Equilibria and Solar-Sail Dynamics in the Sun–Earth–Moon System](https://iris.polito.it/handle/11583/3015641) | The workshop abstract describes time-dependent admissible equilibrium regions and attitude-controlled periodic motion for Earth planetary sunshades in a photogravitational bicircular four-body model. | Directly relevant sunshade research, but the poster is restricted. The abstract supplies no lunar coverage demonstration or reproducible design to adopt. |

The full-text dynamics, methods and results relevant here were inspected for
the first four sources. This is a targeted technical reading, separate from
a full scholarly review. The 2026 conference work was available as an abstract
only. No published trajectory has been implemented as a Terluna solution.

## Force authority at the current mass

At 1 AU, an unloaded Sun-facing 50-g/m² patch has an ideal full-reflection
anti-sunward acceleration of **0.182 mm/s²**. The relaxed photon envelope's
maximum transverse component is **0.0908 mm/s²**. While protecting the central
climate window, the available 13.795% optical fraction reduces those bounds
to **0.0251 and 0.0125 mm/s²**. Payload and power hardware reduce all of them.
Actual UV absorption, power extraction and permitted outgoing directions
still need to be applied.

These are relevant acceleration scales for orbital control. The 2016 study
uses 0.215 mm/s² as its illustrative characteristic acceleration, and the
2025 NASA example uses 0.056 mm/s². Comparing magnitudes establishes a reason
to search at the present mass; it does not demonstrate that our optical
states can follow either orbit.

A cycling central-window tile may use a larger portion of the spectrum when
its redirected rays cannot change lunar illumination. That makes optical
authority depend on the tile's location and service phase. The ray geometry
must identify those intervals. The present model's central spectral limit
continues to apply whenever its rays can reach the surface.

## Inventory in place of expended propellant

The moving aperture's base optical mass is **8.73×10¹² kg**. Its reference
propulsion expends **1.17×10¹³ kg per Julian year**. An expanded fleet can
therefore be worth examining even if only a small share is useful in the
shadow corridor at a time.

Let δ denote the average fraction of total deployed physical area that is
usable projected filtering area in the corridor. It includes orientation
and useful optical state, not simply the fraction of tile centres nearby.
For comparison with the selected aperture, perfect placement gives the
necessary inventory scale

\[
A_{\rm deployed}\geq A_{\rm aperture}/\delta,\qquad
M_{\rm optical}\geq\sigma A_{\rm aperture}/\delta.
\]

For a varying required aperture, the time-averaged required area replaces
the selected fixed area. This area-time bound says nothing about whether
usable tiles can occupy every needed ray at each instant. Spatial overlap,
seams and failed handovers can make a fleet inadequate despite excess area.

| Assumed usable projected fraction | Deployed area multiple | Base optical inventory | Same mass in reference propellant years |
|---:|---:|---:|---:|
| 50% | 2 | 1.75×10¹³ kg | 1.49 |
| 10% | 10 | 8.73×10¹³ kg | 7.46 |
| 1% | 100 | 8.73×10¹⁴ kg | 74.6 |

These are illustrative fractions, not orbital results. The final column
compares masses only. It is not a construction, energy or economic break-even
time. Power equipment, control, replacements, failure reserves, habitats and
deployment are excluded. Recurring sail corrections can still consume
propellant, so the fleet's measured residual impulse belongs in the eventual
lifetime comparison.

[tacking_screen.py](tacking_screen.py) reproduces these algebraic quantities
from the committed holding product and named shared constants. Its
[data product](results/tacking_screen.json) records the input and source
identities. It contains no candidate trajectory or achieved duty fraction.

## Study sequence and remaining acceptance conditions

1. Generate seeds from distant retrograde and synodic-resonant libration
   orbit families, then use multiple shooting or direct collocation to vary
   free states and optical controls. Also retain a broader single-aperture
   trajectory search as a comparison. Keep the 50-g/m² base mass fixed.
2. Continue candidates into the existing DE440 force model. Require bounded
   return states over repeated service intervals, admissible photon momentum
   on every leg, solar visibility with both Earth and Moon blocking, and
   collision and Earth-exclusion clearances. The present eclipse function
   models Earth alone; general orbital trajectories need Moon eclipses too.
   A recurring full-ephemeris trajectory need not repeat exactly each month.
3. Optimize the population, phases and handovers alongside the orbits. Test
   time-dependent, finite-Sun UV coverage across the four-lunar-radius region
   and the climate spectrum at the lunar surface. Use retarded moving-tile
   intersections, local seams and projected area. Do not infer coverage from
   an orbit plot or aggregate area alone.
4. Report achieved useful-area fractions, deployed mass, residual propulsion,
   electricity and blackout recovery. Propagate perturbations with achievable
   optical settings and control rates. Extend promising annual solutions over
   the nodal cycle before comparing lifetime inventories.

Returning point-tile paths and measured service intervals are now available in
[cycling.md](cycling.md). Continuous fleet coverage is the next acceptance
condition. Solar-only orbital control becomes a shield solution when both
have been demonstrated. Habitat wheels and collection hubs can use their
own candidate trajectories; attaching their mass to a tile requires a new
control and inventory calculation.
