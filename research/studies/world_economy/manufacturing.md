# Partial-gravity manufacturing: three specific hypotheses

31 October 2026. This investigation examines candidate process advantages for enclosed lunar factories. It supports the [early-export account](early_exports.md) and supplies reproducible physical calculations for design and economic review. Values below are calculations unless identified as published observations; the candidate factories and their commercial advantages remain unvalidated.

There are defensible reasons to test partial gravity. A process may need enough gravity to remove gas, bring grains into contact or establish a liquid surface, while excessive gravity produces segregation, distortion or larger rotating loads. Lunar gravity can occupy that compromise. The strongest candidates identified here are filled castings with a constrained cure schedule, large spin-cast optical blanks, and selected liquid-phase-sintered parts. Their physical cases differ; none supports a general multiplier for lunar manufacturing productivity.

| Priority | Candidate product | Specific possible lunar advantage | Decisive unresolved issue |
|---|---|---|---|
| 1: most direct experiment | Particle-filled resin castings; possibly filled glasses or glass-ceramics | Gas can rise before immobilization while dense particles remain acceptably distributed during a required processing interval | The actual chemistry and geometry must prevent cheaper Earth controls from reproducing the same result |
| 2: clearest distinct load scaling | Large cured/spin-cast parabolic blanks for orbital mirrors or concentrators | The same passive figure forms at 40.7% of terrestrial angular speed and 16.6% of its centrifugal acceleration | Slower leveling, finishing, transport and competing orbital fabrication may outweigh reduced loads |
| 3: exploratory materials programme | Selected liquid-phase-sintered refractory or composite components | Some gravity could assist skeletal contact and pore movement while reducing segregation | Published Earth/microgravity results do not locate an optimum near lunar gravity |

## 1. Filled castings with an independently constrained cure or thermal cycle

Consider an isolated spherical filler particle and gas bubble in a dilute Newtonian liquid. For an effectively immobile bubble interface, both follow the rigid-sphere Stokes expression:

\[
v_b=\frac{2\rho_f a_b^2 g}{9\mu},\qquad
v_p=\frac{2(\rho_p-\rho_f)a_p^2g}{9\mu}.
\]

Here radii are in metres, viscosity in Pa s, density in kg/m³ and speed in m/s. If bubbles must travel depth \(H\) within a required mobile interval \(T\), while particles may drift at most \(\delta\), the permitted acceleration is

\[
\frac{9\mu H}{2\rho_fa_b^2T}\le g\le
\frac{9\mu\delta}{2(\rho_p-\rho_f)a_p^2T}.
\]

Both bounds have units m/s². The interval exists when
\(H/\delta\le[\rho_f/(\rho_p-\rho_f)](a_b/a_p)^2\).

Nassar et al.'s Table 4 gives theoretical settling inputs for SiC/epoxy: liquid density 1,170 kg/m³, viscosity 0.029 Pa s, particle **diameter** 25 μm and particle density 3,210 kg/m³. The table does not identify viscosity temperature. Their Kemapoxy 150 JM preparation used 4 min mixing at 500 rpm, a 90°C heater for 1 min, and an 8 h pressing period. Room-temperature observations found about 45 min pot life and almost complete settling by 35 min. Thus these inputs do not establish a 100 s gel schedule. [Nassar et al., Table 4 and §§2.1, 4.1](https://www.mdpi.com/2073-4360/14/2/333); [retrieved full PDF](https://pdfs.semanticscholar.org/7409/1e04de5d8aa1aabc5ea3aa2d15c440766411.pdf).

Combining those four table entries with **assumed** 100 μm bubble radius, 10 mm depth, 100 s interval and 0.5 mm settling tolerance gives **1.1154–2.0471 m/s²**. Lunar particle drift is 0.3967 mm and bubble clearance 68.67 s; Earth's drift is 2.3953 mm. Earth also has a feasible interval, 11.37–20.87 s, if the material can be immobilized then. This arithmetic checks a parameter scale rather than demonstrating that the cited resin realizes the proposed window.

The process dimensions and timing were chosen to demonstrate a possible window; they do not identify an optimum from an existing manufacturing recipe.

| Result for the specified 100 s cycle | Moon | Earth | Ideal zero gravity |
|---|---:|---:|---:|
| Bubble travel time over 10 mm | 68.67 s | 11.37 s | No buoyant clearance |
| Particle displacement | 0.397 mm | 2.395 mm | Zero gravitational settling |
| Meets both stipulated tolerances | Yes | No | No |

This establishes a **fixed-process** window. It does not establish superiority over an optimized terrestrial process. At constant viscosity, Earth can reproduce the lunar displacements by actually immobilizing the casting after 16.56 s, or by using 6.04 times the viscosity for 100 s. Merely inspecting earlier does not stop subsequent settling. A useful lunar product therefore needs an independently constrained filling, reaction, diffusion or thermal history that prevents those substitutions without an unacceptable property or cost penalty.

For evolving viscosity, displacements scale with \(g\int dt/\mu(t)\), provided particle sizes and interface factors remain fixed. This preserves the same degeneracy. Jackson and colleagues experimentally validated an integral treatment of bubble motion during changing glass temperatures; the calculation needs the viscosity history, not a convenient mean viscosity. Their work does not validate the assumed 100 s process. [Jackson et al., 2022](https://doi.org/10.1111/jace.18680).

The creeping-flow check is favorable: bubble Reynolds numbers are 0.001175 on the Moon and 0.007094 on Earth. With assumed surface tension 0.03 N/m, the diameter-based bubble Bond numbers are 0.00253 and 0.0153, consistent with nearly spherical bubbles. This surface tension is not a measured property of the proposed formulation. A clean mobile bubble may rise up to 1.5 times the rigid-interface result; particle hindrance, walls or yield stress can instead slow or arrest it. The example fails on the Moon if bubble mobility falls below 68.7% of the assumed value. Particle-rich experiments demonstrate transitions from individual bubble movement to contact-network deformation and fracture-like gas paths, so isolated-sphere extrapolation to heavily filled compounds needs measurement. [Oppenheimer et al., 2015](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2015.00060/full).

Bubble arrival at the surface also does not guarantee venting through a skin or foam film. A bubble born at the full 10 mm depth must appear by 31.33 s in this example. Late volatile formation creates a real chemical clock, but it can make weaker gravity worse. Direct resin-transfer-moulding observations show coupling among volatile release, cure, viscosity and pressure, while also demonstrating terrestrial remedies. [Anders et al., author manuscript](https://bpb-us-w1.wpmucdn.com/sites.usc.edu/dist/f/80/files/2020/05/J230.pdf), DOI 10.1016/j.compositesa.2016.02.024.

Potential products include wear-resistant filled parts and thermally functional encapsulation compounds; the latter is an application hypothesis, not a result demonstrated by the cited wear study. Local basalt/glass fillers would need their own particle-size, chemistry, wetting and performance qualification. Imported SiC/resin cannot be relabelled as lunar feedstock.

**Discriminating test:** identify a real product and required schedule first; measure viscosity, bubble birth/size, particle/agglomerate size, actual venting and final concentration/void gradients. Compare identical schedules and independently optimized schedules at zero, lunar and Earth gravity. Earth controls must include feasible faster cure, vacuum/pressure degassing, rheology modifiers, thinner travel paths and imposed motion. A lunar advantage is falsified if those controls deliver equal quality at lower total cost, or if measured mobility eliminates the acceleration window.

## 2. Large spin-cast or cured parabolic blanks

A rotating liquid with constant vertical gravity has
\[
z(r)=\frac{\omega^2r^2}{2g}=\frac{r^2}{4f},\qquad
\omega=\sqrt{\frac{g}{2f}}.
\]
For the same aperture and focal length, the Moon therefore needs 0.407 times Earth's angular speed. Radial acceleration, rotational kinetic energy for the same rotating inertia, and the centrifugal stress scale \(\rho\omega^2R^2\) fall to 0.166 of Earth's. These ratios describe particular loads, not the entire furnace price or energy consumption. Unlike a rescaled settling hold, a slower Earth spin changes the desired focal length.

For an **illustrative** 10 m diameter, 20 m focal-length blank, both liquid surfaces have 0.3125 m edge-to-centre sag. The required rotation is **1.924 rpm on the Moon versus 4.728 rpm on Earth**; edge speeds are 1.008 versus 2.476 m/s. Lower gravitational support loads and slower rotation could help very large, delicate rotating tooling. Existing terrestrial spin-casting already produces 8.4 m lightweight honeycomb mirror blanks, followed by precision finishing, so scale or delivered cost must supply the reason to relocate. [University of Arizona Mirror Lab](https://astro.arizona.edu/research-areas/design-testing-and-fabrication-large-optics).

This architecture has a useful lower-gravity limit: at zero axial gravity a freely rotating liquid does not make the same gravity-supported paraboloid. Published lunar liquid-mirror studies supply a related physical precedent, including component demonstrations, but their proposed 20–100 m observatories are not manufactured solid exports or completed lunar instruments. [Angel et al., 2008](https://arxiv.org/abs/0806.2241).

Gravity must remain strong enough relative to capillarity at the intended length scale. With illustrative liquid density 2,500 kg/m³ and surface tension 0.3 N/m, the lunar capillary length is 8.60 mm, against 3.50 mm on Earth; the 5 m radius Bond number is approximately 338,000. The bulk is therefore in the gravity-dominated regime, while a larger edge zone needs correction or trimming. These chosen properties do not qualify a particular optical melt.

The tradeoff is slower leveling: in the gravity-dominated, long-wave thin-film approximation the relaxation rate scales with \(\rho g h^3k^2/(3\mu)\). Holding thickness, wavelength and viscosity fixed makes lunar leveling about six times slower. Annealing shrinkage, residual stress, polishing, coating and launch handling can dominate final quality. A blank figured under lunar loading must be designed and verified for its eventual operating support; an Earth telescope would impose greater self-weight, while orbital use has a different unloading distortion. A molten surface is not automatically an optical-quality exported part.

**Discriminating test:** compare realistic maximum furnace/tool loads and final optical figure for the same blank at both gravities, including leveling and annealing time. Compare orbital centrifuges, actively shaped liquids and capillary fabrication as separate production architectures. NASA's FLUTE programme explores large liquid optics in microgravity; its proposed full-scale instrument is not yet evidence that a factory competitor is mature. [NASA FLUTE](https://www.nasa.gov/science-research/astrophysics/what-is-the-fluidic-telescope/). The lunar case fails commercially if savings in rotating/support structure are smaller than extra finishing and delivery costs. High-value orbital buyers are more natural than Earth importers, but their demand is finite.

## 3. Liquid-phase sintering with limited settling and sufficient contact

This is a credible experiment rather than a quantified operating window. Johnson and German's report describes 77 microgravity tungsten-alloy samples alongside terrestrial comparisons. Grain settling on Earth creates contacts that can strengthen the skeleton but also segregates material; microgravity samples can develop large pores and greater distortion. This evidence concerns particular liquid-phase-sintered alloys and thermal histories. It does not mean all densification is impossible without gravity, nor establish lunar gravity as optimal. The author report was read through indexed primary-text passages because direct PDF retrieval failed. [Johnson and German](https://citeseerx.ist.psu.edu/document?doi=e64f1d0986323a4a14527fc409bc07c8574305d1&repid=rep1&type=pdf). NASA's GEDS account independently describes pore coalescence/swelling and the dependence on liquid fraction, grain size, height and interfacial geometry. [NASA GEDS](https://science.nasa.gov/biological-physical/investigations/gravitation-effects-on-distortion-in-sintering/).

A fixed sintering cycle sets diffusion, neck formation and coarsening. The hypothesized optimum is enough gravity-assisted rearrangement and pore movement during that cycle, with acceptable segregation and shape strain. A dilute settling screen would be
\[
g\int\frac{\Delta\rho\,d(t)^2h(\phi,t)}{18\mu_l(t)}dt\le\delta,
\]
where \(h\) is a measured hindrance factor. Effective grain overburden scales as \(\Delta\rho\phi_sgH\), and bulk shape load as \(\rho_{bulk}gH\). Comparisons require evolving contact strength, wetting, liquid fraction and pore connectivity. Pore escape and granular rearrangement are not independent Stokes spheres, so inventing lower/upper lunar thresholds from this schematic would create false precision.

**Discriminating test:** select an alloy and cycle known to segregate at 1g and distort or retain pores at near-zero gravity; map density, composition and dimensional fidelity across intermediate acceleration. Include pressure-assisted sintering, a presintered skeleton and adjusted liquid fraction as controls. Only an intermediate minimum in combined defects, surviving those controls, supports the proposed advantage. Translating a tungsten-alloy result to regolith glass-ceramics requires a separate materials programme.

## Boundaries for the early-export account

Ordinary soldering is a useful adverse control. NASA's plated-through-hole tests found mean flux-core void area of 11.0% near zero gravity, 7.75% at nominal lunar gravity and 3.21% at 1g; the lunar group had only 16 samples and the authors cautioned against statistical conclusions from the small partial-gravity groups. Capillary gap-filling and buoyant gas removal suggest a possible geometry-specific compromise, but these data do not establish one. [Watson et al., Table 2](https://ntrs.nasa.gov/api/citations/20060013420/downloads/20060013420.pdf?attachment=true). Generic claims that lunar solder, foam or printed regolith is superior are unsupported by this analysis.

A rotating orbital factory can reproduce lunar acceleration. For example, a 10 m centrifuge radius requires about 3.85 rpm for 1.624 m/s²; a larger radius reduces acceleration gradients. Its support, containment and material handling have costs, while lunar factories have dust, transport and site costs. The economic hypothesis is therefore a combination of a useful acceleration, suitable nearby feedstock, stationary rock-supported facilities and buyer access. Partial gravity alone supplies no defensible monopoly rent.

The repository does not specify an early factory's commissioned equipment, throughput, yield or customer contracts for these products. The proposed manufacturing revenue row must remain an independent capacity/price scenario; these physical hypotheses neither validate its assumptions nor justify adding a gravity premium. A real process advantage should first appear as measured yield, equipment or operating-cost improvement, then enter the territorial production account once. Moon-sited factories produce Moon value added; an orbital factory remains array production even if lunar residents own it.

## Reproduction and evidence

[partial_gravity.py](partial_gravity.py) computes the SiC-parameter casting window, Earth time/viscosity controls, a longer constant-viscosity sensitivity and the spin-cast geometry using the shared gravity constants. [partial_gravity_inputs.json](partial_gravity_inputs.json) separates tabulated material inputs from assumed process constraints; [results/partial_gravity.json](results/partial_gravity.json) binds their identities and producer. Run `python -m research.studies.world_economy.partial_gravity` for an ignored-run result.

A 45-minute hold at the same constant viscosity fails the stated particle tolerance at both lunar and Earth gravity. This sensitivity establishes why the assumed 100-second interval matters; it does not simulate the source resin's evolving viscosity or establish that all other accelerations fail. The dilute displacement model omits container walls and loses physical applicability after particles reach a boundary.

The [experimental evidence audit](manufacturing_sources.md) records additional direct manufacturing comparisons and their access scopes. [sources.json](sources.json) records targeted reading for this investigation, including incomplete retrievals; [checks.json](checks.json) records numerical and repository verification. No materials experiment, completed optical design, manufacturing cost comparison or empirical lunar validation was performed.
