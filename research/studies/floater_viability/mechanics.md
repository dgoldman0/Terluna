# Mechanical feasibility of biological floaters

The earlier aerial-ecology screen identified a minimum radius from a chosen skin mass. A load-bearing gas envelope
also becomes harder to build as its vertical extent grows. Its upper surface experiences a hydrostatic pressure
excess even when gas and air pressures match at the bottom. Larger floaters therefore gain buoyant mass and incur
larger pressure, structural, damage and control requirements together.

A homogeneous membrane and a lobed body with strong tendons have materially different feasible ranges. The present
calculations reject some combinations across **all radii within a stated model**, while leaving others mechanically
possible in a first screen. They do not establish a biological maximum size or certify an organism. Materials,
membranes, attachment points, repeated storms, growth and repair remain coupled research requirements.

[mechanics.py](mechanics.py) contains pure functions used by the combined study. The atmospheric inputs are the
committed [sky-ship product](../sky_ships/results/sky_ships.json); gravity is `GM/(Rmoon+h)²` from shared constants.
The [source register](mechanics_sources.json) records primary evidence and its limits. The existing gas table assumes
98% lifting species and 2% ambient air by volume; mixture mass and actual hydrogen mass must be distinguished.

## Pressure increases across a large gas body

For a body of radius r in a locally uniform atmosphere, let Δρ be ambient density minus lifting-gas density.
Hydrostatics gives

\[
\frac{d\Delta p}{dz}=\Delta\rho g,
\qquad \Delta p_{\rm crown}=p_{\rm bottom}+2r\Delta\rho g.
\]

At 10 km, the retained product has air density 1.204 kg/m³, hydrogen-mixture lifting allowance 1.098 kg/m³ and gravity
1.60568 m/s². The hydrostatic head is 3.526 Pa per metre of radius: 35.3 Pa at r=10 m, 353 Pa at r=100 m and 3.53 kPa
at r=1 km. Calling a balloon “zero pressure” at its lower opening does not make its upper membrane stress-free.
NASA balloon design work identifies the same size-dependent pressure requirement and the advantage of local lobes
and global load tapes ([Said et al. 2004](https://ntrs.nasa.gov/citations/20040079786)).

A gust/shear sensitivity adds `Cp × ρair v²/2`. The sample calculations use Cp=1 and a 10 m/s relative gust, adding
60.2 Pa, beside a 10 Pa bottom pressure. These are imposed loads, not a sampled storm return level. Co-moving drift
reduces mean horizontal relative wind; it cannot remove velocity differences across a large deformable body or
sudden storm outflows. Pressure coefficients can be negative and spatially variable, with local collapse and
wrinkling possible; adding one positive pressure allowance does not solve that fluid–structure problem.

The code includes an independent isothermal check using exponential hydrostatics and exact varying lunar gravity.
For kilometre scales it differs little from the local-density head in this example. For very large bodies, the
body's height becomes appreciable relative to the density scale height, approximately 52 km here. A 20-km-tall body
cannot use one parcel's temperature, density or weather without qualification. The exact isothermal check is still
an idealized consistency calculation; Terluna's actual nonisothermal profiles must replace it in a full design.

## A homogeneous spherical skin has a finite useful size window

A thin spherical membrane under uniform differential pressure p carries tension per unit length `p r/2`, hence
required thickness `t = p r/(2 σallow)`. Here σallow means allowable stress of the assembled, hydrated, aged membrane,
including seams and acceptable defects. It is not a short coupon's breaking stress.

Sizing the **entire** skin to the crown pressure is a deliberately conservative local-pressure screen. It does not
prove a real spherical body can transmit its nonuniform gravity and payload loads in that shape. A complete model
needs an equilibrium shape, distributed masses, biaxial constitutive behavior, local attachment details and stability.
The skin mass is `4πr² ρs t`; separate gas barriers or biological coverings add their own areal mass. If one membrane
performs both structural and barrier roles, its mass is entered once.

With bottom pressure and gust allowance combined into p₀, additional barrier areal mass μb, and fraction f of full
lift allocated to **all actual supported mass**, the remaining wet payload per projected area is

\[
B(r)=a r-b r^2-c,
\quad a=\frac{4 f\Delta\rho}{3}-\frac{2\rho_s p_0}{\sigma_{\rm allow}},
\quad b=\frac{4\rho_s\Delta\rho g}{\sigma_{\rm allow}},
\quad c=4\mu_b.
\]

For a>0, its optimum occurs at `r=a/(2b)`. A requested wet payload B* is possible in this screen only if
`a²−4b(c+B*) ≥ 0`; the roots delimit the radius interval. The code also enforces a chosen minimum skin thickness
when evaluating individual bodies. The analytic interval assumes stress sets thickness; a thickness floor can only
reduce the available payload. The equations are a restriction of this architecture, not a general theorem forbidding
other shapes or materials.

For example, with no added pressure or barrier, material density 1,500 kg/m³ and allowable stress 1 MPa, the maximum
wet payload is about 50.7 kg/m² when f=1, or 12.7 kg/m² when f=0.5. The old heavy example of 10 kg dry biomass/m² at
90% water requires 100 kg wet/m², so that homogeneous 1 MPa screen cannot support it at any radius. Merely enlarging
the organism does not rescue that case.

Adding the sample 70.2 Pa pressure allowance and 0.05 kg/m² separate barrier gives:

| Assembled allowable stress | Loaded fraction f | Peak wet payload per projected area | Radius interval for 10 kg wet/m² | Radius interval for 100 kg wet/m² |
|---|---:|---:|---:|---:|
| 0.1 MPa | 0.5 | No positive payload | None | None |
| 1 MPa | 0.5 | 6.22 kg/m² at 24.6 m | None | None |
| 1 MPa | 1.0 | 36.9 kg/m² at 59.2 m | 8.8–109.7 m | None |
| 10 MPa | 0.5 | 119 kg/m² at 336 m | 14.7–657 m | 201–471 m |

The densities and allowable stresses are explicit sensitivities, not measured properties of a candidate organism.
Full loading f=1 leaves no lift margin for rain or future growth. Lower f means a partly filled gas bladder or another
neutral-trim arrangement; it does not mean the organism flies indefinitely with unopposed excess upward force.

## Lobes reduce local film span; tendons still carry the whole body

A lobed membrane can make local curvature radius rₗ far smaller than body radius r. Approximating each film panel
as cylindrical gives local hoop thickness `tfilm = p rₗ/σfilm`. Strong meridional tendons then react the global
pressure loads. The idealized model uses `N≈πr/rₗ` tendons of length πr, sharing a bounding hemispherical force
`F=pπr²`. Each tendon has section `F/(Nσtendon)`. Therefore total tendon mass is

\[
M_{\rm tendons}=\pi^2\rho_t p r^3/\sigma_t.
\]

Increasing the number of lobes reduces panel span but does not eliminate this global force or tendon mass. As
hydrostatic pressure grows with body size, the tendon mass fraction also grows. An attachment-force allowance can
be passed separately, but a real load-path solution must prevent both omission and double counting of payload load.
The present force estimate excludes node fittings, sleeves, partitions, tendon redundancy and load sharing after
damage. It uses spherical-equivalent volume and a stated 1.1 film-area multiplier, rather than claiming a solved
pumpkin shape.

The following example scans radii from 0.1 m to 10 km. It uses film allowance 1 MPa, tendon allowance 100 MPa,
material densities 1,500 kg/m³, 50 μm minimum film thickness, the same additional barrier and pressure allowance,
and f=0.5. The lobe radius is 1 m, reduced to r/5 for bodies smaller than 5 m.

| Body radius | Bounding pressure | Sphere film thickness | Sphere wet payload | Lobed wet payload, including ideal tendons |
|---|---:|---:|---:|---:|
| 0.1 m | 70.6 Pa | 0.050 mm | −0.43 kg/m² | −0.48 kg/m² |
| 1 m | 73.7 Pa | 0.050 mm | 0.23 kg/m² | 0.18 kg/m² |
| 10 m | 105 Pa | 0.527 mm | 3.96 kg/m² | 6.35 kg/m² |
| 30 m | 176 Pa | 2.64 mm | 5.92 kg/m² | 20.3 kg/m² |
| 100 m | 423 Pa | 21.1 mm | −53.8 kg/m² | 68.2 kg/m² |
| 300 m | 1.13 kPa | 169 mm | −796 kg/m² | 196 kg/m² |
| 1 km | 3.60 kPa | 1.80 m | −10,057 kg/m² | 539 kg/m² |
| 3 km | 10.65 kPa | 15.97 m | −93,640 kg/m² | 620 kg/m², panel approximation fails |
| 10 km | 35.33 kPa | 176.7 m | Negative; thin-shell approximation fails | Negative; panel approximation fails |

Negative payload rejects the case. Positive payload is mass arithmetic only. Thin-film formulas are marked outside
their screening range when thickness exceeds 1% of the relevant curvature radius; that is an explicit modeling
cutoff, not a physical failure threshold. At 3 km the lobed calculation already needs panels too thick for its stated
thin-panel approximation despite positive arithmetic payload. Many smaller positive cases can still fail on creep,
attachments, gas retention, nutrition, control or storms.

NASA's experience provides a relevant distinction: effective braided-tendon strength can fall well below direct
scaling of constituent-fiber strength. Slow loading and material age matter
([Sterling 2004](https://ntrs.nasa.gov/api/citations/20040040072/downloads/20040040072.pdf)). Likewise, modern analyses
include viscoelastic film behavior, thermal cycling and local tendon-sleeve details
([Wakefield and Bown 2017](https://ntrs.nasa.gov/citations/20180007750)).

## What biological material measurements actually establish

Two primary studies support investigating biological or biologically derived structural skins while leaving the
assembled-envelope requirement open.

[Matas et al. 2005](https://doi.org/10.3732/ajb.92.3.462) tested isolated tomato cuticle strips, 3×9 mm and approximately
7.2 μm thick. Samples experienced 10–45 °C and 40% RH, 80% RH or immersion, with 20-minute loading steps. Hydration
and warmth reduced strength; the reported wet strength was approximately 21 MPa below the thermal transition and
14 MPa above it. These are excised, screened, millimetre specimens, not a large living sheet or long-term creep limit.

[Kriechbaum and Bergström 2020](https://doi.org/10.1021/acs.biomac.9b01655) produced cellulose-nanofibril/gelatin/tannin
films with 33 MPa wet tensile strength. Testing used 20×3 mm specimens, a 10 mm gauge, 1 mm/min extension and one hour
of prior immersion. The hybrid's dry thickness was 54.2±1.3 μm and water uptake approximately 25%. This is an engineered
nonliving film, with no demonstrated hydrogen barrier, giant assembled membrane or lifelong self-repair.

Consequently, neither coupon value is adopted as a design allowance. The 0.1–10 MPa assembled-film range asks how
much performance survives wetting, aging, sustained load, joints, flaws and biological growth. The 100 MPa tendon
case asks for a different directional material and load architecture. An effective hydrated density of 1,500 kg/m³
is an imposed mass sensitivity. Measured wet thickness, water uptake and mechanical stress must be kept consistent;
dry density cannot be paired silently with wet geometry.

## Neutral trim and altitude control

Let V be the fixed external volume, fgas its lifting-mixture volume fraction, ρg the mixture density and M the actual
supported mass, including living wet biomass, films, partitions, fittings and emergency stores. A gas bladder plus
ambient-air ballonet is neutral when

\[
f_{\rm gas}=\frac{M}{(\rho_a-\rho_g)V}.
\]

The total mass is `M + fgas ρg V + (1−fgas) ρa V = ρa V`. The code verifies that closure and reports actual mixture and
hydrogen mass, separately from full-capacity inventory. If fgas>1, this architecture cannot trim neutral with that
external volume. Actual partitions and ballonet fabric need mass; `partition_mass` makes the extra sheet budget
explicit. A spherical inner bladder's minimum area is `4πr² fgas^(2/3)`; an outer-envelope-sized barrier supplies an
alternative leakage-area bound. Neither ideal geometry establishes a growth or partition design. Pressure sizing above retains a full-height gas body as a conservative reserve-expansion case even when nominal trim uses less gas. A particular inner bladder may have a smaller head, but its interfaces, geometry and variable filling would need their own load solution.

Altitude stability also depends on the architecture. In the ideal soft, ambient-pressure case with gas and air at the
same temperature, `V=nRT/P` and `ρair V=nMair`; net gas lifting mass is approximately `n(Mair−Mgas)`, independent of
altitude. Such a slack bladder has no automatic restoring density level merely because atmospheric density falls
upward. A fully inflated, fixed-volume, closed-gas body can have restoring buoyancy as ambient density changes, but
then gas heating changes pressure and membrane stress. A fixed outer volume with air-exchanging ballonets instead
changes its total mass through valves and pumping. Gas production, temperatures, pressure control, reserve expansion,
air ballast and flight behavior must therefore be modeled together.

## Rain, damage and loss of gas

One millimetre of retained water over projected area adds 1 kg/m². If this remains uncompensated after drainage and
trim, a still-air steady-drag screen gives

\[
v=\sqrt{2\Delta M g/(\rho_a C_D A)}.
\]

At 10 km with Cd=0.47, residual loads of 1, 5 and 10 kg/m² correspond to downward speeds 2.38, 5.33 and 7.53 m/s.
This is a force-balance diagnostic, not a descent trajectory, impact prediction or storm-safety result. Acceleration,
added air mass, deformation, turbulent currents and actual drag all remain. It shows why reliable drainage, ballast
exchange and the timing of control are central requirements even for gently drifting organisms.

A small fixed puncture and a propagating tear are different failure modes. `initial_hole_leak` calculates only initial
pressure-driven continuum outflow, `Q=Cd A sqrt(2Δp/ρgas)`, with Cd=0.6 as a sensitivity. It rejects pressure excess above
10% of ambient, beyond its low-compressibility screen; compressible flow and choking require a different calculation
([NASA Glenn](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/)). Zero overpressure
means no advective flux in this equation, not no hydrogen diffusion.

For the 100 m radius example at its static crown pressure, a 1 mm diameter circular hole initially passes about
3.37 m³/day of mixture, a 1 cm hole 337 m³/day, and a 10 cm hole 33,677 m³/day. Those rates cannot be extrapolated
unchanged through an entire collapse. Crack propagation, pressure redistribution, rain intrusion and repair can
dominate the consequence. Existing balloon experiments explicitly treat manufacturing flaws and repairs under
biaxial loading ([Portanova 1989](https://ntrs.nasa.gov/citations/19890020453)); the floater has no measured tear-arrest
or repair performance yet.

## Mechanical conclusion and next discriminating evidence

Large photosynthetic floaters remain worth pursuing, especially bodies whose local gas barrier and global support
fibres perform different jobs. A viable candidate must close the full wet mass and pressure account at its chosen
size. The earlier lift-only radius was a starting requirement; the present size windows prevent “make it larger”
from being used as a universal answer.

The next discriminating measurements are assembled wet-membrane biaxial strength and creep, gas permeability after
strain and damage, tendon/junction durability, and repair throughput under pressure. The next model is a compartmented
body with solved shape, actual tissue placement, thermal/pressure cycles, drainage and neutral trim. No safe species,
maximum biological size, population density or storm tolerance is selected by these screens.

Fifteen lightweight tests pass, covering hydrostatic and exact-profile limits, hemisphere force balance, sphere and
lobe scaling, separate barrier mass, analytic size-window roots, neutral gas/air closure, finite puncture flow and
rain-load drag balance. Full study integration and repository checks belong to the combined runner.
