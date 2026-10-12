# Mechanics of biological aerophytes

**A gas body's skin has a useful size window, and lobed films with tendons widen it to kilometre radii in the
arithmetic.** The pressure across the skin grows with the height of the gas column, so a plain sphere's skin mass
grows as R⁴ while its lift grows as R³. Lobes shrink the film's local span and leave the global load to tendons,
whose mass grows more slowly. The full account of joints, gusts, damage and creep that brings these limits in is in
[structure](structure.md); this note holds the underlying screens.

[mechanics.py](mechanics.py) holds the pure functions the coupled study uses. Air comes from the
[sky-ship product](../sky_ships/results/sky_ships.json) and gravity is GM/(R_Moon + h)² from shared constants. The gas
table's lifting mixture is 98% hydrogen and 2% ambient air by volume. The [source register](mechanics_sources.json)
records each source and how far it was read.

## Pressure grows across a tall gas body

**At 10 km the pressure difference across the skin grows by 3.53 Pa per metre of radius: 35 Pa at 10 m, 353 Pa at
100 m and 3.53 kPa at 1 km.** With Δρ the air density less the lifting mixture's, hydrostatics gives

\[
\frac{d\Delta p}{dz}=\Delta\rho\, g,
\qquad \Delta p_{\rm crown}=p_{\rm bottom}+2r\Delta\rho g.
\]

The 10-km air has density 1.204 kg/m³, mixture lift 1.098 kg/m³ and gravity 1.60568 m/s². A balloon open at its base
still carries this head at its crown. NASA's balloon design work identifies the same size-dependent pressure and the
advantage of local lobes with global load tapes (Said et al. 2004).

A gust allowance adds C_p ρ_air v²/2; with C_p 1 and a 10 m/s relative gust it adds 60.2 Pa beside a 10 Pa base
pressure. The storm-based design gusts that replace this allowance in the coupled cases are in
[structure](structure.md). An independent isothermal check with exponential hydrostatics and exact lunar gravity
agrees with the local head to within 0.6% up to 1 km radius. A body 20 km tall would span a large part of the air's
scale heights (52 km for pressure and 62 km for density at 10 km) and need the real temperature profile.

## A plain sphere has a size window

**A uniform sphere's skin must be t = pR/(2σ) thick, so above the size where the head dominates, its mass grows as R⁴
and payload falls.** Here σ is the allowable stress of the assembled, hydrated and aged membrane, seams included. The
screen sizes the whole skin to the crown pressure, which is conservative. With base and gust pressure p₀, a separate
barrier of areal mass μ_b and a share f of full lift given to all carried mass, the wet payload per projected area is

\[
B(r)=a r-b r^2-c,
\quad a=\frac{4 f\Delta\rho}{3}-\frac{2\rho_s p_0}{\sigma},
\quad b=\frac{4\rho_s\Delta\rho g}{\sigma},
\quad c=4\mu_b.
\]

It peaks at r = a/(2b), and a payload B* fits only if a² − 4b(c + B*) ≥ 0, between the two roots. A minimum
thickness can only lower the payload.

**At 1 MPa a plain sphere carries at most 50.7 kg/m² of wet payload with all its lift, or 12.7 kg/m² with half.** The
skin density is 1,500 kg/m³ with no added pressure or barrier. A heavy design with 10 kg of dry biomass per square
metre at 90% water carries 100 kg/m² of water and tissue, and fits at no size at this strength. With the 70.2 Pa
allowance and a 0.05 kg/m² barrier:

| Assembled allowable stress | Share of lift used | Peak wet payload | Radii for 10 kg/m² | Radii for 100 kg/m² |
|---|---:|---:|---:|---:|
| 0.1 MPa | 0.5 | below zero | none | none |
| 1 MPa | 0.5 | 6.22 kg/m² at 24.6 m | none | none |
| 1 MPa | 1.0 | 36.9 kg/m² at 59.2 m | 8.8–109.7 m | none |
| 10 MPa | 0.5 | 119 kg/m² at 336 m | 14.7–657 m | 201–471 m |

A body using its full lift has nothing left for rain or growth; a share below one is held by a partly filled gas
bladder with an air ballonet, at neutral buoyancy.

## Lobes shrink the film's span; tendons carry the body

**Lobes cut the film's thickness to t = p r_l/σ for local radius r_l, and tendons of total mass π²ρ_t p R³/σ_t carry the
global load.** The model uses N ≈ πR/r_l meridional tendons of length πR sharing the hemisphere force pπR², a 1.1 film
area multiplier and the spherical-equivalent volume. More lobes thin the film without changing the tendons' mass.
Carrying the crown pressure along each tendon overstates the exact equatorial load by up to 20%
([structure](structure.md)). Node fittings, sleeves, partitions and redundancy are costed in the structure note.

With a 1 MPa film, 100 MPa tendons, 1,500 kg/m³ materials, a 50 µm film floor, the same barrier and allowance, half
the lift and lobes of 1 m radius (R/5 below 5 m):

| Body radius | Bounding pressure | Sphere film | Sphere wet payload | Lobed wet payload with tendons |
|---|---:|---:|---:|---:|
| 0.1 m | 70.6 Pa | 0.050 mm | −0.43 kg/m² | −0.48 kg/m² |
| 1 m | 73.7 Pa | 0.050 mm | 0.23 kg/m² | 0.18 kg/m² |
| 10 m | 105 Pa | 0.527 mm | 3.96 kg/m² | 6.35 kg/m² |
| 30 m | 176 Pa | 2.64 mm | 5.92 kg/m² | 20.3 kg/m² |
| 100 m | 423 Pa | 21.1 mm | −53.8 kg/m² | 68.2 kg/m² |
| 300 m | 1.13 kPa | 169 mm | −796 kg/m² | 196 kg/m² |
| 1 km | 3.60 kPa | 1.80 m | −10,057 kg/m² | 539 kg/m² |
| 3 km | 10.65 kPa | 15.97 m | −93,640 kg/m² | 620 kg/m², panels too thick for the thin-panel screen |
| 10 km | 35.33 kPa | 176.7 m | below zero; thin-shell screen fails | below zero |

A thin-film result is flagged when thickness exceeds 1% of the curvature radius. Braided balloon tendons fall below
the direct scaling of their fibres' strength, and slow loading and age matter (Sterling 2004); full balloon analyses
include viscoelastic film, thermal cycling and tendon-sleeve detail (Wakefield & Bown 2017).

## What biological materials show

**Wet biological coupons reach 14–33 MPa.** Isolated tomato fruit cuticle, 3 × 9 mm strips about 7.2 µm thick tested
at 10–45 °C, 40% and 80% RH and in buffer with 20-minute load steps, held about 21 MPa wet below its thermal
transition and 14 MPa above it; warmth and hydration lowered strength and stiffness (Matas et al. 2005). A
nanocellulose–gelatin–tannin film, 54.2 ± 1.3 µm dry and taking up about 25% water, held 33 MPa wet after an hour's
immersion, tested on 20 × 3 mm strips at a 10 mm gauge and 1 mm/min (Kriechbaum & Bergström 2020). The structure note
carries these through joints and sustained load to the 5.8 MPa assembled skin. The 0.1–10 MPa film and 100 MPa tendon
cases here are sensitivities, and 1,500 kg/m³ is a hydrated-density sensitivity.

## Neutral trim

**A gas bladder with an ambient-air ballonet floats neutral when f_gas = M/((ρ_a − ρ_g)V).** M is all carried mass,
including wet tissue, films, partitions, fittings and stores. The total M + f_gas ρ_g V + (1 − f_gas)ρ_a V equals
ρ_a V, which the code checks. If f_gas exceeds one the external volume cannot float that mass. `partition_mass`
carries the partition sheet's mass; a spherical inner bladder's least area is 4πr² f_gas^(2/3), and the outer envelope
bounds the gas-contact area from above. Pressure sizing keeps the full-height gas body as the reserve-expansion case.

**A soft bladder has no restoring height; a full one does.** In a slack bladder with gas and air at one temperature,
V = nRT/P and the net lift n(M_air − M_gas) is the same at every height. A full, closed, fixed-volume body gains
restoring buoyancy as the air thins above it, and gas heating then changes its pressure and skin stress. A fixed outer
hull with air ballonets changes its total mass by valves and pumping. Gas production, temperature, pressure, reserve
expansion, ballast and flight are therefore one coupled problem ([trim cycle](trim_cycle.md), [storms](storms.md)).

## Rain, punctures and tears

**Each millimetre of retained rain adds 1 kg/m², and uncompensated loads of 1, 5 and 10 kg/m² sink a body at 2.38,
5.33 and 7.53 m/s.** The still-air balance v = sqrt(2ΔMg/(ρ_a C_D A)) with C_D 0.47 at 10 km shows why drainage and
trim timing matter even for drifting organisms; storm rain rates are in [storms](storms.md).

**A puncture's first flow is Q = C_d A sqrt(2Δp/ρ_gas): at the crown of a 100-m-radius body a 1-mm hole passes
3.37 m³ of mixture a day, a 1-cm hole 337 m³ and a 10-cm hole 33,677 m³.** `initial_hole_leak` applies this with
C_d 0.6 for overpressures up to 10% of ambient, the low-compressibility range; faster flows need compressible and
choked-flow models (NASA Glenn). A zero overpressure drives no flow through the hole, while diffusion continues. Tear
growth, pressure redistribution and repair set the consequence of large openings; balloon experiments treat flaws and
repairs under biaxial load (Portanova 1989), and the arrested-tear times are in [structure](structure.md).

## Checks

Fifteen tests in [test_mechanics.py](test_mechanics.py) cover the hydrostatic and exact-profile heads, the hemisphere
force balance, sphere and lobe scaling, the separate barrier mass, the analytic window's roots, neutral gas–air
closure, the puncture flow and the rain-load drag balance.

Sources: [Said et al. 2004](https://ntrs.nasa.gov/citations/20040079786),
[Sterling 2004](https://ntrs.nasa.gov/api/citations/20040040072/downloads/20040040072.pdf),
[Wakefield & Bown 2017](https://ntrs.nasa.gov/citations/20180007750),
[Matas et al. 2005](https://doi.org/10.3732/ajb.92.3.462),
[Kriechbaum & Bergström 2020](https://doi.org/10.1021/acs.biomac.9b01655),
[Portanova 1989](https://ntrs.nasa.gov/citations/19890020453),
[NASA Glenn mass flow](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/).
