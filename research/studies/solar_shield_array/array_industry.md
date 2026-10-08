# Heat, computing and industry of the array

On 8 October 2026 the author kept thermodynamic computing as a research
direction for the solar shield and habitat array, for the array's own control
and for scientific and AI services. The author also asked whether the heat the
array already handles could support it without adding much to its burden. The
same day the author made this the last work on the solar-shield branch. One
later branch, `research/array-industry`, is planned for heat reuse, general
computing and industry of the array, with the shield's reflectivity trim and
momentum store and their power in its plan
([decisions.md](../../decisions.md)). This note gathers the first
investigation:
- what thermodynamic computers are and what they can do;
- a screen of where the array's heat arises and what it is worth
  ([array_heat.py](array_heat.py));
- the shield's own loads;
- the planned branch.

## What the investigation finds

**The tiles' own heat repays no device.** The ring fleet's film absorbs
65–83 PW of sunlight. A tile facing the Sun runs near 175 K, with 0.1–0.4 mK
between its faces. Carnot work across every film in the fleet would come to a
quarter of a terawatt. The best thermoradiative diode measured so far, in
mercury cadmium telluride, would give 27 TW only if it covered all
2.8 billion km². The light the tiles pass and redirect is the resource.

**The heat worth reusing comes from what the array does.**
- Collectors give off 1.9–3.4 W of heat for each watt of electricity, at
  300–315 K. Thermoradiative converters on their backs, at the best measured
  performance, would add two to four parts in a hundred thousand.
- The film plant that remakes the fleet's film, 2.3–4.6 Gt a year, runs at
  0.15–2.9 TW near 2,000 K, the best-grade heat in the system.
  Thermophotovoltaic cells lining it could return 9–24% of its energy.
- Computing rejects its heat near 330 K, at the bottom of the cascade, on
  830 km² of radiator panel per terawatt.

**Where heat is released matters more than recovering it.** A terawatt used on
the Moon adds 0.026 W/m² to its global heat budget, so the 100 TW reference
case would add about Earth's present anthropogenic forcing. Released in orbit
at the ring radius, 0.19% of it reaches the Moon. Energy-intensive industry and
computing belong in orbit.

**Thermodynamic computers are processors for sampling, inference and some
linear algebra.** They take their randomness from thermal noise and their
energy from electricity, like any computer, so they run on waste heat only
through a heat engine. Their advantages over GPUs, ten- to ten-thousand-fold,
are modelled, and the machines measured so far are small. For the array they
suit the sampling-heavy parts of the fleet's estimation and scheduling, and the
science, inference and generative services the author named. They are to be
judged by useful results per joule on matched tasks.

**The shield's own loads depend on hardware the planned branch chooses.**
- Trim devices that need power to hold their state would draw 0.7–13 PW over
  the fleet, and devices switched and left 0.006–16 TW.
- The momentum store's rotors weigh 0.1–1.7 Gt across the fleet and lose
  4–155 GW, depending on their size and speed.
- Holding the stack's eccentricities takes up to a fifth of the outermost
  rings' eccentricity drive, about the whole of the trim's reach there.

## Thermodynamic computing

**How it works.** A thermodynamic computer lets a physical system settle under
its own fluctuations and reads the answer from the states it visits. Its
couplings encode the problem. Coupled electrical oscillators driven by noise
settle into a Gaussian distribution whose covariance is the inverse of their
coupling matrix, so sampling their voltages inverts the matrix
([Melanson and colleagues](https://doi.org/10.1038/s41467-025-59011-x);
[Aifer and colleagues](https://doi.org/10.1038/s44335-024-00014-0)).
Probabilistic bits flip at random with biases set by their neighbours and so
sample a Boltzmann distribution, which carries optimization and inference
problems. The noise explores the problem's states, where a digital computer
pays for every pseudo-random number. The energy comes from the power supply:
it programs the couplings, holds the noise at the level the problem needs,
reads out samples and resets. In present designs the readout's
analogue-to-digital converters take most of it, 40–200 µW against 0.2–3 µW for
the computing core in a simulated chip. Normal Computing's board ran on noise
from a random-number generator, since ambient noise was too weak to drive it
([Quanta, 2026](https://www.quantamagazine.org/thermodynamic-computers-go-with-the-energy-flow-20260715/)).
Computing driven by heat flowing between two temperatures exists in theory
([Lipka-Bartosik and colleagues](https://doi.org/10.1126/sciadv.adm8792)). Any
of these machines can run on waste heat only through a heat engine, whose
electricity is Carnot-limited like any other.

**What it does well and badly.** These machines sample Gaussian and Boltzmann
distributions. They serve Bayesian inference and uncertainty estimates, Ising
and combinatorial optimization, energy-based and generative models, and linear
algebra done by sampling. Their accuracy is statistical: the error falls as
the square root of the number of samples, so the time to a given accuracy grows
as its inverse square. They do deterministic high-precision arithmetic and
general code badly, and every one built so far is a co-processor to a digital
host. Their size is limited by the connections between cells, the cost of
readout and the calibration of each cell.

**Where the hardware stands.** No computer that draws on physical noise has yet
been measured beating a GPU end to end on the same task; its reported
advantages are modelled. FPGAs that emulate probabilistic bits digitally have
been measured beating GPU code on sparse Ising problems:
- **Normal Computing.** The measured board has eight cells. It sampled Gaussian
  distributions and inverted an 8×8 matrix to 25–31% error in 12,000 samples. A
  model of larger chips crosses over a GPU near dimension 3,000 and gains an
  order of magnitude at 10,000. The CN101 chip, taped out in August 2025,
  targets up to a thousandfold; no silicon results were published by October
  2026.
- **Extropic.** Its X0 test chip draws its randomness from transistor noise at
  hundreds of attojoules to a few femtojoules per sample, at tens of megahertz.
  A hardware model of denoising thermodynamic models matches GPU models on a
  small image benchmark at about a ten-thousandth of the energy per sample,
  against GPU energy estimated from operation counts. The Z1 chip of 269,568
  probabilistic bits, under 1 W, has been taped out, with early access planned
  for 2027.
- **Probabilistic bits.** On FPGAs with a few thousand of them, measured
  machines reach 140–185 flips per nanosecond. One ran 5–18 times faster than
  optimized GPU and TPU code; another took 0.05 nJ per flip, against 22 nJ
  quoted for a V100 GPU on the same problems. An integrated magnetic chip is
  reported to run 96,000 spins at under 40 fJ each.

By the Landauer bound, erasing a bit costs at least kT ln 2,
2.9×10⁻²¹ J at 300 K. Present logic signals take about 700 aJ, some
200,000 times that ([Frank and Conte](https://www.sandia.gov/app/uploads/sites/210/2022/06/FrankConte-HotCarbon22-v4SAND.pdf)),
and probabilistic devices sit as far above it. The search found no radiation
test of a thermodynamic or probabilistic computer. Magnetic tunnel junctions,
the devices of some probabilistic bits, came through 147 kGy of gamma rays
unchanged ([Montoya and colleagues](https://doi.org/10.1038/s41598-020-67257-2)).

**Measuring their computing.** There is no accepted measure, and each group
uses its own:
- probabilistic bits: flips per second and joules per flip;
- optimizers: time to a solution;
- samplers: samples per joule at a stated fidelity;
- generative hardware: joules per sample or per token.

FLOPs do not carry over, and published comparisons often model one side
([Stroev and colleagues](https://doi.org/10.1038/s42005-024-01870-9)). For the
array the useful measure is results per joule on a matched task: samples at a
stated fidelity, an optimization at a stated quality, an inference at a stated
accuracy. It has to count the digital host, readout, calibration and the time
between independent samples on both sides.

**How much computing a terawatt buys.** Present accelerators do about
1.4×10¹² operations per joule at 16-bit precision
([NVIDIA H100](https://www.nvidia.com/en-us/data-center/h100/)), so a terawatt
gives about 1.4×10²⁴ operations a second. Each terawatt also needs
2,450–3,670 km² of collectors at 20–30% and 830 km² of radiator panel at
330 K ([below](#where-the-arrays-heat-arises)).

**For the array.** The fleet's control is an estimation and scheduling problem
over 28 million tiles: their states with uncertainty, keeping and service
plans, failures. Its sampling-heavy parts suit thermodynamic co-processors
beside radiation-tolerant digital controllers. The same machines could serve
the science, inference and generative work the author named. Like any computer
they draw electricity and need radiators.

## Where the array's heat arises

[array_heat.py](array_heat.py) follows each path in closed form from the
comparison's products and the sourced assumptions it lists.

| Heat | Power | Temperature | What recovering it could give |
|---|---:|---:|---|
| Sunlight the tiles absorb | 65–83 PW | 173–176 K facing the Sun | 0.25 TW at Carnot across every film |
| Collectors for 100–1,000 TW delivered | 1.9–3.4 W per watt delivered | 300–315 K | 2–35 GW from measured thermoradiative diodes on their backs |
| Film plant, 73–147 t of glass a second | 0.15–2.9 TW | near 2,000 K | 0.01–0.7 TW from thermophotovoltaics, 9–24% of its energy |
| Computing | per terawatt | near 330 K | warmth for habitats beside it |

**The tiles.** The ring fleet's 2.8 billion km² of film intercepts
1,900–2,400 PW of sunlight over its orbits. That is eight to ten times what the
aperture intercepts, since each ring's tiles circle the whole orbit. The coated
annulus film absorbs 3.0% of the sunlight on it and the window stack 4.5%, the
protection domain's figures, 65–83 PW in all. The protection domain's silica
and titania tables give the films' faces 0.37–0.44 (annulus) and 0.44–0.68
(window) of a blackbody's thermal emission. A tile facing the Sun therefore
runs at 173–176 K, and colder through the rest of its orbit. The shaded face's
share of the heat crosses the film with 0.1–0.4 mK between the faces, so
Carnot work across every film in the fleet would come to 0.25 TW. The best
thermoradiative diode yet measured gives 9.4 mW/m² at room temperature facing a
150 K surface
([Radchenkov and colleagues](https://doi.org/10.1063/5.0265431)). It would give
27 TW only if it covered all 2.8 billion km². It is made of mercury cadmium
telluride, and at 175 K it would give far less. Ideal limits for a 300 K
emitter facing space are 48–153 W/m²
([Buddhiraju and colleagues](https://doi.org/10.1073/pnas.1717595115)), and a
175 K tile has a ninth of a 300 K body's emission to work with. The author's
"water from pumice" holds for the tiles: the light they pass and redirect is
the resource, and the heat they absorb is too cold and too even to repay a
device.

**The collectors.** Converting light at 20–30%, a collector gives off
1.9–3.4 W of heat for every watt of electricity, at a multijunction array's
absorptance of 0.88. Facing the Sun and radiating from both faces it runs at
300–315 K; geostationary arrays measure 41–46 °C. For the 100–1,000 TW
electricity cases that is 190–3,400 TW of heat, the largest stream after the
tiles' own.
[Landis](https://ntrs.nasa.gov/citations/20220003553) proposes thermoradiative
converters on an array's shaded face as a bottoming cycle. At the best
measured 9.4 mW/m² they would add 2–35 GW to 100–1,000 TW, two to four parts
in a hundred thousand. Ideal limits allow 50–150 W/m² from a 300 K face,
against the array's 270–410 W/m² of electricity. The array then runs warmer,
losing about 0.23% of its power per degree, behind a front that has to emit
little in the infrared.

**The film plant.** Remaking the fleet's film every 10–20 years means
2.3–4.6 Gt of glass a year, 73–147 t a second. Raising silica from room
temperature into the melt takes 2.0–2.3 MJ/kg at 2,000–2,200 K, and evaporating
it takes at least 13–15 MJ/kg (NIST-JANAF). At 2–20 MJ/kg the plant runs at
0.15–2.9 TW near 2,000 K, where fused silica softens and forms glass. This is
the best-grade heat in the system. Thermophotovoltaic cells have converted
36–41% of the radiation from emitters at 2,100–2,700 K, with the cells held at
25 °C ([LaPotin and colleagues](https://doi.org/10.1038/s41586-022-04473-y);
[Tervo and colleagues](https://doi.org/10.1016/j.joule.2022.10.002)). Suppose
cells line the plant's hot surfaces, capture 30–60% of its heat and convert
30–40% of that:
- they return 9–24% of the plant's energy, 0.01–0.7 TW;
- they leave themselves 0.03–1.2 TW to reject near 330 K, on 20–1,000 km² of
  panel.

The glass itself makes a poor source, since silica is transparent across the
cells' band. The cells face the refractory surfaces that hold the melt.

**Computing.** Processors reject their heat near 330 K. There a two-faced
panel with an emissivity of 0.9 rejects 1,210 W/m², so a terawatt of computing
needs:
- 830 km² of panel, weighing 2.5–11.6 Mt at 3–14 kg/m²;
- 2,450–3,670 km² of collectors to power it.

The 3–14 kg/m² spans NASA's goal to the ISS's panels. An engine working off
this heat into a 250 K radiator could take at most 24%, on a radiator three
times the size. Habitats beside the processors could use some of the heat for
warmth.

All of these paths begin as sunlight. Collectors take their electricity from
it before any of it becomes heat, so recovering their leftover heat adds little
to what they already make. The plant's heat is worth recovering because the
plant has to reach 2,000 K anyway.

## Where heat is released

All of this heat leaves by radiation, and where it is released decides what it
does to the Moon. A terawatt used on the Moon adds 0.026 W/m² to its
global-mean heat budget, 0.009% of the 293 W/m² of sunlight the climate design
gives it. The 100 TW reference case would add 2.6 W/m², close to the 2.7 W/m²
of Earth's present anthropogenic forcing
([IPCC AR6](https://www.ipcc.ch/report/ar6/wg1/)). A terawatt radiated
evenly from the ring radius, 20,000 km, reaches the Moon in proportion to the
solid angle the Moon covers there: 0.19%, about 5×10⁻⁵ W/m². Industry and
computing in orbit therefore keep their heat out of the Moon's climate.
Radiators turned edge-on to the Moon and kept out of the aperture's sight
lines put even less on it. Power delivered to the Moon ends as heat there,
counted in the Open Moon's allocation, one of the three power ledgers the
engineering domain keeps apart.

## The shield's own loads

The comparison left the shield's own loads, the reflectivity trim and the
momentum store, unsized for power. They depend on hardware not yet chosen, so
the screen bounds the classes the planned branch has to choose between.

**Trim.** The trim changes the reflectivity of the redirected band by up to a
quarter, from devices over 25–89% of each tile. Across the fleet that is
0.7–2.5 billion km², so each milliwatt per square metre the devices draw is
0.7–2.5 TW. The device classes compare as follows:
- **Held films.** Liquid-crystal films of the kind IKAROS flew reflect
  specularly only while powered, at 1–5 W/m²
  ([Mori and colleagues](https://www.issfd.org/ISSFD_2009/AOCSI/Mori.pdf);
  [Ma and colleagues](https://doi.org/10.1002/adom.201600668)). Over the
  fleet's devices that is 0.7–13 PW.
- **Electrochromic films with memory.** These take about 10–80 J/m² for a
  quarter's change in reflection; window products take 0.2–2 kJ/m² for a full
  tint ([Villa and colleagues](https://doi.org/10.1016/j.enbuild.2024.114293)).
  Switched 2–20 times an orbit, they draw 0.06–16 TW.
- **Display-type bistable media.** These take about 1–5 J/m² per update
  ([Lin and colleagues](https://doi.org/10.3390/mi15091076)), 0.006–1 TW.

Photovoltaics on the tile could power any of them, over under 0.003% of its
area for the switched devices. The held films would need petawatts and up to
2.3% of every tile in photovoltaics, which rules them out. The switching energy
and the number of switches an orbit are the figures to design for.

**Its duty.** The trim changes the band's reflection, which gives 2R of a
film's normal pressure 2R + A: about three-quarters for the annulus film and
two-thirds for the window stack. A quarter of the band therefore moves 17–20%
of the push. Holding the stack's eccentricities, which the radius layout
relies on, takes a steady share of each ring's eccentricity drive
([held_eccentricity.json](results/held_eccentricity.json)):
- 5% over the stack;
- up to 19.5% at the outermost rings below the Moon's orbit plane at 0.6 km
  clearance, and 26% at 1 km;
- more than a tenth on 334 rings at 0.6 km.

At the outermost rings holding alone would take the whole trim, which also
holds the tiles' attitude. The tilts have to carry part of it; the parked
keeping study settles how much.

**The momentum store.** The store's orbit-mean demand lies about the tile's
roll and pitch axes, 2,450 and 3,650 N·m, with almost none about its normal.
Its rotors therefore spin about axes in the tile's plane and reach out of it
toward the next rings, at most a few hundred metres. A rotor of radius r at rim
speed v holds momentum H with a mass of H/(v r) and an energy of H v/(2r). For
the 1.9×10⁹ N·m·s each tile needs:

| Rotor radius | Rim speed | Rotor mass per tile | Energy per tile | Rotors in the fleet | Losses, cycled once an orbit |
|---:|---:|---:|---:|---:|---:|
| 100 m | 300 m/s | 62 t (3.8%) | 2.8 GJ | 1.7 Gt | 15–47 GW |
| 100 m | 1,000 m/s | 19 t (1.1%) | 9.3 GJ | 0.5 Gt | 52–155 GW |
| 400 m | 300 m/s | 16 t (0.9%) | 0.7 GJ | 0.4 Gt | 4–12 GW |
| 400 m | 1,000 m/s | 4.7 t (0.3%) | 2.3 GJ | 0.13 Gt | 13–39 GW |

The percentages are of a mean tile's 1,640 t, and the losses take 85–95% of the
energy back each cycle. A NASA study of a 3.2 km solar power satellite found
about the same peak momentum, 2×10⁹ N·m·s. It proposed five to seven
space-built wheels of 350 m radius and 6 t of aluminium each
([Wie and Roithmayr](https://ntrs.nasa.gov/citations/20010071579)). Compact
wheels cannot do it: at the 17.5 N·m·s per kilogram of the ISS's gyroscopes, a
tile's store would weigh 65 times the tile, and 8 times at the 150 N·m·s/kg
foreseen for advanced flywheels. Wide rotors are the lightest and the
cheapest to run. Even the narrowest and fastest lose 0.2 TW across the fleet,
small beside its other flows, while their mass, 0.1–1.7 Gt, belongs in the
46 Gt inventory.

## The planned branch

`research/array-industry` takes up heat reuse, general computing and industry
of the shield and habitat array, with the shield's reflectivity trim and
momentum store and their power. It opens when the author chooses; on
8 October too many branches were open off a stale main. It starts from:
- this note, [array_heat.json](results/array_heat.json) and
  [held_eccentricity.json](results/held_eccentricity.json);
- the [integrated comparison](integrated_comparison.md) and its ledger: the
  ring fleet's 28 million tiles and 46 Gt, its film upkeep, the electricity
  cases and the gates, among them net-positive electricity for habitat;
- the reconstructed industrial architecture of 19 September
  ([engineering/reference/industrial_architecture](../../../engineering/reference/industrial_architecture/README.md)):
  its power ledgers, heat rejection, source industry, manufacturing closure
  and matter stream;
- the engineering domain's [network model and supply ledger](../../../engineering/network/README.md),
  and the three power ledgers it keeps apart (civilization-wide power, the Open
  Moon's allocation, and where heat is released);
- the September feasibility report's call to close where industrial heat is
  placed ([research/baselines/feasibility](../../baselines/feasibility/report.md)).

Its questions, in order:
1. **The shield's own loads as hardware.**
   - A trim device that holds its state without power and switches a few times
     an orbit. Its electrodes, wiring and sensors have to be made of common
     elements: a milligram per square metre of anything is 2.8 Mt across the
     fleet.
   - The momentum store's rotors within the clearance to the next ring, with
     their mass and losses.
   - How attitude, keeping and eccentricity holding share the trim and the
     tilts.
   - Then the ring fleet's electricity gate.
2. **The film plant.**
   - Where it orbits, and its chain from worn tile to new: recovery, remelting
     or evaporation, drawing or deposition, coating, and assembly into 10 km
     tiles.
   - Its 73–147 t/s and 0.15–2.9 TW, and its heat cascade.
   - The logistics of 4,000–8,000 tile exchanges a day, failed tiles included.
3. **The power and heat ledgers.** The array's collection, conversion and
   delivery, with the heat of each step placed in orbit or on the Moon.
4. **Computing.**
   - Demand for the fleet's own control and for science and AI services.
   - Placement beside collectors and industry, with radiators edge-on to the
     Moon and outside the shield's light path.
   - The radiation environment at the ring radius, outside Earth's
     magnetosphere most of each month.
   - Conventional and thermodynamic processors compared by results per joule
     on matched tasks.
5. **Industry beyond the shield.** Which energy-intensive industries belong in
   orbit and which on the Moon, judged by where their heat goes and by their
   need for vacuum and energy. They have to meet the outflow requirement S7 and
   the six traffic-safety rules.

The rest of the solar-shield work stays parked on this branch
([the integrated comparison](integrated_comparison.md#where-the-comparison-stands-8-october)):
- the keeping in the regressing frame, with the eccentricity holding and the
  along-track keeping;
- the atmosphere domain's refinements;
- what waits on inputs this machine lacks;
- a full fleet's coverage, handovers and the shadow on Earth.

## Sources and running

[array_industry_sources.json](array_industry_sources.json) lists every source,
what it was used for and how much of it was read. The first investigation
(Codex, 8 October) read sections of six. Research agents re-read those the
same day, corrected the register and found the rest.

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.held_eccentricity  # about 6 minutes on 2 workers
python -m research.studies.solar_shield_array.array_heat         # seconds
python -m pytest engineering/test_heat_recovery.py research/studies/solar_shield_array/test_array_heat.py research/studies/solar_shield_array/test_held_eccentricity.py
```

[held_eccentricity.py](held_eccentricity.py) rebuilds the stack's coupled
radius layout from [ring_layout.json](results/ring_layout.json) and writes
[results/held_eccentricity.json](results/held_eccentricity.json).
[array_heat.py](array_heat.py) reads it with the integrated ledger, the
attitude schemes and the protection domain's films, and writes
[results/array_heat.json](results/array_heat.json), in closed form. Both are
screens: they run no dynamics and model no device or plant.
