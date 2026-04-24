# Report: Retention of a Shielded, Terraformed Lunar Atmosphere

## Core conclusion

A dense lunar atmosphere with Earthlike N₂/O₂/Ar composition can land in the intended retention range when the Moon has both a solar-wind shield and a strong high-energy spectral shield. The modeled atmosphere can cycle on timescales of hundreds of millions of years, and in colder upper-atmosphere cases it can reach or exceed a billion years.

The viable operating window is narrow and physically demanding. The central requirement is a cold, weakly escaping nitrogen exosphere. The atmosphere survives longest when the upper atmosphere stays near **230–250 K**, high-energy solar leakage remains around **0.1% or lower**, and nitrogen escape remains closer to a Jeans-tail regime than to a full planetary wind.

The strongest modeled cases produce loss rates around **10–50 kg/s**, corresponding to atmospheric cycling times of **2–8 billion years**.

The more nominal cases produce loss rates around **70–200 kg/s**, corresponding to cycling times of roughly **0.45–1.3 billion years**.

The marginal-but-maintainable cases produce loss rates around **200–400 kg/s**, corresponding to cycling times of roughly **200–500 million years**.

The major physical hinge is nitrogen. Oxygen and argon remain comparatively well retained. Nitrogen forms the high-altitude, light-enough, weakly bound tail that controls the maintenance burden.

---

## What “cycling time” means here

The modeled lunar atmosphere has a surface pressure of **1.2 atm** and an Earthlike bulk composition. Its mass is estimated from hydrostatic column weight:

[
M_{\rm atm} \approx \frac{P,4\pi R_{\rm Moon}^{2}}{g_{\rm Moon}}
]

For the Moon, this gives:

[
M_{\rm atm} \approx 2.8\times10^{18}\ {\rm kg}
]

The cycling time is the time required for maintenance top-off to replace one full atmosphere:

[
\tau = \frac{M_{\rm atm}}{\dot M}
]

This framing converts escape rates into intuitive long-term maintenance requirements.

| Average loss rate | Full atmosphere cycling time |            Annual top-off |
| ----------------: | ---------------------------: | ------------------------: |
|           10 kg/s |           ~8.9 billion years | ~0.32 million tonnes/year |
|           30 kg/s |           ~3.0 billion years | ~0.95 million tonnes/year |
|           90 kg/s |           ~1.0 billion years |  ~2.8 million tonnes/year |
|          125 kg/s |           ~720 million years |  ~3.9 million tonnes/year |
|          300 kg/s |           ~300 million years |  ~9.5 million tonnes/year |
|          800 kg/s |           ~110 million years |   ~25 million tonnes/year |

This is one of the most important quantitative anchors. A billion-year atmosphere permits losses near **90 kg/s**. A 300-million-year atmosphere permits losses near **300 kg/s**. These are large industrial flows, yet they remain small compared with the total undertaking of terraforming and maintaining a lunar biosphere-scale environment.

---

## System architecture assumed by the model

The modeled system includes three major parts.

The first is the atmosphere itself: a dense, Earthlike mixture dominated by nitrogen and oxygen, with argon as a minor constituent.

The second is a solar-wind shield. This is assumed to suppress ion pickup, sputtering, direct plasma erosion, and most solar-wind-driven atmospheric stripping. This assumption is foundational. A dense lunar atmosphere without strong plasma protection faces a much harsher retention problem.

The third is a high-energy spectral shield. This shield attenuates EUV, FUV, and soft X-rays while allowing most visible and near-infrared sunlight through. The goal is to reduce upper-atmosphere heating, dissociation, ionization, and hot-atom production without turning the lunar surface into a dark world.

This spectral shield is powerful because the relevant high-energy solar bands carry only a tiny fraction of total solar power. The structure remains Moon-scale in area, but its radiative power burden is far smaller than a visible-light sunshade. Blocking high-energy photons attacks the part of sunlight that drives exosphere heating rather than the part that drives ordinary illumination and climate.

The model therefore evaluates a Moon that receives normal visible sunlight, has plasma shielding, and has strongly reduced thermosphere-heating radiation.

---

## Model development

The model began with simple atmospheric mass accounting and progressed toward a reduced aeronomy/escape model. Each added layer tested whether the viability window survived stronger physical constraints.

### Atmospheric mass and pressure

The first step was the mass of a 1.2 atm lunar atmosphere. The Moon’s weak gravity means a given surface pressure requires less atmospheric mass than Earth would require over Earth’s surface, though the lunar atmosphere is still enormous in absolute terms.

The result was approximately:

[
2.8\times10^{18}\ {\rm kg}
]

This number sets the maintenance budget. A loss of tens to hundreds of kg/s corresponds to geological cycling times rather than rapid collapse.

### Hydrostatic vertical structure

The next step modeled a spherical hydrostatic atmosphere above the Moon. Hydrostatic balance was solved with lunar gravity decreasing with altitude:

[
g(r)=\frac{GM_{\rm Moon}}{r^2}
]

The model used:

[
\frac{dP}{dz}=-\rho g(z)
]

with:

[
P=nkT
]

The key result from this layer is that a lunar atmosphere is vertically huge compared with Earth’s. For N₂/O₂ near Earthlike temperatures, lunar scale height is around **50–55 km**, compared with Earth’s ~8.5 km. The atmosphere thins slowly. That pushes the exobase to high altitude, where escape speed is lower.

The Moon’s retention problem is therefore controlled by high-altitude structure, not simply by surface pressure.

### Exobase determination

The model defines the exobase where mean free path becomes comparable to scale height:

[
\ell \sim H
]

The mean free path is:

[
\ell = \frac{1}{\sqrt{2}n\sigma}
]

The scale height is:

[
H=\frac{kT}{mg}
]

This makes the exobase an outcome of temperature, composition, density, gravity, and collision cross section. In a dense lunar atmosphere, the exobase can sit thousands of kilometers above the surface. That is a crucial result because escape speed falls with altitude.

### Jeans escape

The first escape model used classical Jeans escape. The Jeans parameter is:

[
\lambda = \frac{GMm}{rkT}
]

The Jeans flux has the approximate form:

[
\Phi_J =
\frac{n v_{\rm th}}{2\sqrt{\pi}}
(1+\lambda)e^{-\lambda}
]

with:

[
v_{\rm th}=\sqrt{\frac{2kT}{m}}
]

This stage showed that N₂ escape is strongly temperature-sensitive. A cold exobase gives slow loss. A modestly warmer exobase gives much faster loss.

At upper-atmosphere temperatures near **250 K**, Jeans-only N₂ loss can sit in the tens of kg/s range. At temperatures approaching **270 K**, the atmosphere becomes much leakier. At still warmer temperatures, the hydrostatic assumption begins to approach a weak-outflow regime.

### Diffusive separation

The model then added diffusive separation above a homopause. Below the homopause, the atmosphere remains well mixed. Above it, species follow individual scale heights.

This upgrade made nitrogen clearly dominant in the escape budget.

O₂ is heavier than N₂ and sits lower. Argon is heavier still. N₂ extends upward farther, reaches lower escape-speed regions, and dominates loss. With diffusive separation included, the maintenance burden becomes mostly a nitrogen top-off problem.

This result was stable across a range of homopause altitudes. The precise homopause value mattered less than the upper-atmosphere temperature.

### Parametric photochemistry

A first-order photochemistry module added budget terms for:

* H escape from trace upper-atmosphere water,
* hot atomic N escape,
* hot atomic O escape.

This module was intentionally parametric. It estimates escape budget pressure rather than solving a full reaction network.

The result was encouraging. In shielded cases, photochemical escape contributed a few to a few tens of kg/s. That matters, yet it usually remained smaller than N₂ thermal or transitional escape.

Hydrogen escape is important for water inventory over long timeframes, but in the dry-upper-atmosphere cases it remains manageable. Hot N and O escape become significant when FUV leakage rises, but they fit inside the atmospheric budget when leakage is near 0.1% and the upper atmosphere stays cold.

### Upper-atmosphere energy balance

The model then stopped assigning the upper-atmosphere temperature directly and instead solved for it from a simplified energy budget.

The budget includes:

* residual high-energy solar heating,
* radiative cooling,
* optional thermal coupling to the lower atmosphere,
* escape cooling.

This upgrade changed the interpretation. Strong high-energy shielding matters, but once high-energy heating is suppressed enough, lower-atmosphere heat coupling can become the limiting source of upper-atmosphere warmth. In that regime, stronger spectral shielding provides diminishing returns unless the atmosphere can radiatively cool or thermally isolate the exobase.

This means a real design must focus on upper-atmosphere thermal architecture, not only photon blocking. The atmosphere needs mechanisms that keep the thermosphere cold: radiative cooling, limited upward heat transport, possibly composition choices that support infrared cooling, and global circulation that avoids warm exobase regions dominating escape.

### Transitional outflow

The final major upgrade added a non-Jeans escape check.

The model computes both classical Jeans escape and a Parker-wind-style outflow estimate. It then blends between them using an outflow coupling parameter. This is necessary because many favorable N₂ cases sit at Jeans parameters around:

[
\lambda_{\rm N2} \sim 6\text{–}10
]

That range is marginal. Classical Jeans escape gives the low-loss end. A fully organized Parker-like wind gives the high-loss end. The real lunar atmosphere would likely sit somewhere between these limits depending on collisionality, heat deposition, cooling, composition, and flow geometry.

This upgrade sharpened the conclusion. The concept works well when N₂ escape remains weakly or moderately enhanced above Jeans. It misses the billion-year target when N₂ behaves close to a full Parker wind.

---

## Main numerical results

### Outflow sensitivity at a 250 K upper atmosphere

With a 250 K upper floor and 0.1% high-energy leakage, the model gives the following range:

| Escape behavior              | Total loss |       Cycling time |
| ---------------------------- | ---------: | -----------------: |
| Jeans only                   |   ~50 kg/s | ~1.8 billion years |
| Weak outflow enhancement     |   ~72 kg/s | ~1.3 billion years |
| Moderate outflow enhancement |  ~123 kg/s | ~730 million years |
| Strong outflow enhancement   |  ~270 kg/s | ~330 million years |
| Full Parker-like upper bound |  ~783 kg/s | ~110 million years |

This table is the most direct statement of uncertainty. The same thermal atmosphere can be long-lived or merely medium-lived depending on how the N₂ exobase actually escapes.

Weak-to-moderate outflow supports the target. Full Parker-like outflow pushes the atmosphere into rapid engineered cycling.

### Temperature sensitivity with moderate outflow

With moderate outflow coupling and 0.1% leakage:

| Upper-atmosphere floor | Total loss |       Cycling time |
| ---------------------: | ---------: | -----------------: |
|                  230 K |   ~11 kg/s |   ~8 billion years |
|                  240 K |   ~31 kg/s | ~2.9 billion years |
|                  250 K |  ~123 kg/s | ~730 million years |

This is the strongest temperature result. A 10 K change has a large effect. A 20 K change has an enormous effect.

The billion-year target becomes comfortable near **230–240 K**. It remains plausible near **250 K** if outflow is modest. It becomes increasingly difficult above **255–260 K**.

### Lower-atmosphere thermal coupling

With thermal coupling from the lower atmosphere included, the upper atmosphere warms. In a representative coupling case:

| Radiative floor | Solved upper T | Total loss |       Cycling time |
| --------------: | -------------: | ---------: | -----------------: |
|           230 K |         ~247 K |   ~77 kg/s | ~1.2 billion years |
|           240 K |         ~250 K |  ~125 kg/s | ~720 million years |
|           250 K |         ~253 K |  ~190 kg/s | ~470 million years |

This result is probably closer to the practical design problem. The thermosphere needs both reduced high-energy heating and favorable heat balance. A cold radiative floor helps. Strong upward coupling warms the exobase and reduces retention time.

### High-energy leakage sensitivity

With lower-atmosphere coupling, moderate outflow, and a 240 K radiative floor:

| Effective high-energy leakage | Solved upper T | Total loss |       Cycling time |
| ----------------------------: | -------------: | ---------: | -----------------: |
|                         0.01% |         ~249 K |  ~109 kg/s | ~830 million years |
|                         0.03% |         ~250 K |  ~112 kg/s | ~800 million years |
|                          0.1% |         ~250 K |  ~125 kg/s | ~720 million years |
|                          0.3% |         ~252 K |  ~162 kg/s | ~560 million years |
|                          1.0% |         ~256 K |  ~315 kg/s | ~290 million years |

This shows that high-energy shielding remains essential, yet the design becomes thermally limited after leakage falls below roughly 0.1%. In that regime, better cooling and weaker upward heat coupling matter as much as additional spectral attenuation.

---

## Strengths of the model

The model captures the main atmospheric-retention structure of the problem. It connects surface pressure, lunar gravity, composition, exobase altitude, thermal escape, photochemical loss, and top-off burden in one framework.

It captures spherical geometry. This matters because the exobase can sit thousands of kilometers above the lunar surface, where gravity and escape speed differ strongly from surface values.

It captures diffusive separation. This is essential. Treating the atmosphere as bulk-mixed underestimates the importance of nitrogen’s high-altitude tail.

It captures the dominant role of upper-atmosphere temperature. The model repeatedly shows that thermal state drives escape more strongly than small changes in surface pressure or homopause altitude.

It captures the correct budget logic. The relevant question is total loss rate compared with the atmosphere’s mass. The model translates kg/s into cycling times and annual top-off masses.

It captures the importance of high-energy solar shielding. EUV/FUV/soft-X-ray suppression directly reduces heating and photochemical loss, and it does so without requiring broadband sunlight reduction.

It captures the main uncertainty through an explicit outflow parameter. Rather than hiding the Jeans-to-wind transition, the model exposes it. That makes the feasibility judgment clearer.

It distinguishes between different gases. N₂, O₂, and Ar behave differently, and the model shows that nitrogen is the controlling bulk species.

---

## Limitations of the model

The model is one-dimensional and globally averaged. A real lunar atmosphere would have day-night structure, circulation, wave activity, weather, and local exobase variation. These effects can redistribute heat and change where escape is strongest.

The upper-atmosphere energy balance is simplified. It uses reduced terms for heating, cooling, lower-atmosphere coupling, and escape cooling. A higher-fidelity model would calculate radiative transfer, species-specific heating, molecular conduction, eddy diffusion, and infrared line cooling.

The photochemistry module is parametric. It budgets H, hot N, and hot O escape, while a full model would solve altitude-dependent reaction networks for N₂, O₂, O, N, NOx, O₃, H₂O, H, OH, and related species.

The transitional outflow model is a bracket. It interpolates between Jeans escape and Parker-like outflow. The real solution requires solving the momentum, continuity, and energy equations through the exobase transition.

The solar-wind shield is treated as effective. A complete system model would include shield leakage, storm-time compression, cusp losses, reconnection-like escape pathways, and residual sputtering.

Surface sinks are still budget terms rather than modeled reservoirs. O₂ reaction with regolith, nitrogen fixation, adsorption, cold trapping, impact gardening, and volatile burial can matter over geological time if they reach tens of kg/s equivalent.

The current atmosphere is Earthlike by assumption. Alternate compositions could improve retention. A heavier buffer gas mixture would lower scale height and reduce nitrogen loss, while changing climate, breathing mixture, and biospheric constraints.

---

## How well the model resolves feasibility

The model resolves the first-order feasibility question well enough to classify the concept.

A shielded lunar atmosphere is not automatically removed by thermal escape. A dense Moon atmosphere can sit in a maintenance range compatible with hundreds of millions to a billion years of cycling.

The model identifies the feasibility window:

* solar wind stripping suppressed,
* high-energy heating suppressed by roughly 99.7–99.9%,
* upper atmosphere held near 230–250 K,
* nitrogen escape weak-to-moderately enhanced above Jeans,
* photochemical and surface sinks kept to tens of kg/s,
* maintenance top-off available at million-tonne-per-year scale.

Inside that window, the system reaches the target.

The model also identifies the failure path:

* upper atmosphere warms above roughly 255–260 K,
* N₂ escape approaches full outflow behavior,
* lower-atmosphere thermal coupling overwhelms cooling,
* high-energy leakage rises toward ~1% or more,
* nonthermal and surface sinks add hundreds of kg/s.

In those cases, the atmosphere still exists as an engineered atmosphere, but its cycling time shifts toward 100–300 million years or below. The billion-year goal becomes difficult.

The model’s current resolution is therefore strong at the architectural level and medium at the predictive level. It can say what conditions are required. It can estimate loss ranges. It can rank sensitivities. It cannot yet produce a single definitive loss number for a final design.

---

## Conditions for a successful long-lived atmosphere

### Solar-wind shield

The solar-wind shield must suppress plasma erosion strongly. This includes ion pickup, sputtering, direct stripping, and storm-time loss. The atmospheric retention result assumes those losses remain small compared with the thermal and photochemical budget.

A practical target would be residual plasma-driven losses well below **10–30 kg/s** for billion-year ambitions, or below **100 kg/s** for a hundreds-of-millions-year system.

### High-energy spectral shield

The high-energy shield should suppress thermosphere-heating EUV/FUV/soft-X-ray input to around:

[
f_{\rm leak} \sim 10^{-3}
]

or lower in effective heating terms.

That corresponds to roughly:

[
99.7\text{–}99.9%
]

effective reduction of the radiation that matters for upper-atmosphere heating and photochemistry.

The shield can allow visible and near-infrared sunlight through. The goal is spectral selectivity, not darkness.

### Cold upper atmosphere

The thermosphere/exosphere should stay near:

[
230\text{–}250\ {\rm K}
]

for strong long-term retention.

The billion-year target is most comfortable when the upper atmosphere is near **230–240 K**. It remains plausible near **250 K** when outflow is weak or moderate. It becomes marginal above **255–260 K**.

This temperature requirement is the central physical design requirement.

### Controlled nitrogen outflow

N₂ escape should remain closer to Jeans escape than to full Parker-like outflow. Moderate enhancement above Jeans can be tolerated. Full Parker-like flow drives losses too high for the billion-year target in most modeled cases.

This condition depends on the detailed exobase structure, collisionality, heating altitude, and cooling rates.

### Low upper-atmosphere water leakage

Water vapor reaching the upper atmosphere should be limited. H escape can become a water-loss concern if the upper atmosphere is wet and FUV leakage is significant.

A dry stratosphere/mesosphere analogue helps. Cold traps, atmospheric circulation, and water-management chemistry matter.

### Manageable surface sinks

Oxygen and nitrogen sequestration into the lunar surface should remain budgeted. A dense atmosphere over fresh lunar regolith will alter the surface chemically. Over long timescales, even kg/s-scale sinks matter.

The atmosphere is easiest to maintain if the surface reaches chemical passivation or if regolith reactions are included in the top-off plan.

---

## Practical maintenance burden

The viable regimes imply top-off rates on the order of millions of tonnes per year.

Examples:

| Total loss |            Annual top-off |
| ---------: | ------------------------: |
|    30 kg/s | ~0.95 million tonnes/year |
|    90 kg/s |  ~2.8 million tonnes/year |
|   125 kg/s |  ~3.9 million tonnes/year |
|   300 kg/s |  ~9.5 million tonnes/year |

For a civilization capable of lunar terraforming, this is manageable compared with the initial atmospheric import, construction of shields, climate management, and biosphere support. The maintenance problem is large in contemporary industrial terms and modest in planetary-engineering terms.

The replacement material is mostly nitrogen. Oxygen and argon are thermally retained much better. Nitrogen supply and recycling therefore become strategic resources.

---

## Design implications

The atmospheric composition should be treated as a design variable. Earthlike composition is a demanding baseline because N₂ is light enough to escape efficiently from the lunar exobase. A modified breathing atmosphere with lower nitrogen fraction, higher argon fraction, or another inert buffer could improve retention. Human habitability and biosphere compatibility would constrain that trade.

The high-energy shield should be optimized for upper-atmosphere outcomes rather than raw attenuation alone. The most important quantity is effective deposited heat in the thermosphere, followed by photodissociation and photoionization rates.

The atmosphere should include radiative cooling pathways in the upper atmosphere. Trace species that cool efficiently could materially improve retention if they avoid causing additional photochemical escape. CO₂, NO, O₃, and related species may matter here.

The lower atmosphere should avoid excessive upward thermal coupling. Strong eddy mixing into the thermosphere can warm the exobase. A stable upper atmosphere with controlled mixing favors retention.

The solar-wind shield should be evaluated by residual kg/s loss, not simply by magnetic field strength or standoff distance. The atmosphere’s retention budget cares about mass flux.

The system should be designed around nitrogen conservation. Since N₂ dominates loss, monitoring and active nitrogen replenishment become the core maintenance function.

---

## Remaining work that matters most

The next high-value model is a one-dimensional transonic escape solver for N₂. It should solve continuity, momentum, and energy equations from the collisional upper atmosphere through the exobase region. This would determine whether the real solution resembles Jeans escape, weak outflow, or Parker-like flow.

The second high-value model is a real thermosphere energy-balance model. It should include altitude-dependent heating, radiative cooling, conduction, eddy diffusion, and composition-dependent thermal structure. This would determine whether 230–250 K is physically maintainable.

The third high-value model is a basic N/O/H photochemical network. It should solve production and loss rates for N₂, O₂, O, N, NOx, O₃, H₂O, H, and OH under the shielded spectrum.

The fourth high-value model is a surface-sink budget. It should estimate O₂ uptake, nitrogen fixation, volatile cold trapping, adsorption, and impact gardening over million- to billion-year intervals.

These additions would turn the current screening model into a credible preliminary retention model.

---

## Feasibility statement

The shielded lunar atmosphere concept remains viable under defined conditions.

The atmosphere reaches the desired **hundreds-of-millions-year** class across a reasonably broad range of cold, shielded cases.

The atmosphere reaches the **billion-year** class when the upper atmosphere is colder, outflow is weak-to-moderate, and nonthermal losses stay within a tens-of-kg/s budget.

The atmosphere reaches the **multi-billion-year** class in the coldest modeled cases, especially near a 230–240 K upper atmosphere with weak outflow.

The controlling requirement is clear:

[
T_{\rm upper} \lesssim 240\text{–}250\ {\rm K}
]

with N₂ escape staying below roughly:

[
\dot M_{\rm N2} \sim 30\text{–}100\ {\rm kg/s}
]

for billion-year ambitions.

A terraformed Moon with a solar-wind shield and strong EUV/FUV/soft-X-ray filtering can therefore maintain a dense atmosphere on geological timescales in the modeled parameter space. The design challenge is not atmospheric mass alone. It is the creation of a cold, stable, weakly escaping nitrogen exosphere.
