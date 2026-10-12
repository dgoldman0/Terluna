# Computing in the array: supply against demand

The author called computing the valuable resource that turns the array's otherwise wasted energy into use, and asked
what it is for: the array's own control, the Moon's management, heavy simulation, and the far-out tasks people
consider, set against what the array can supply ([register](../research/decisions.md#solar-shield-and-habitat-array),
9 October). This note answers with a screen, [computing.py](computing.py), which writes
[results/computing.json](results/computing.json) from the shield study's
[array heat](../research/studies/solar_shield_array/results/array_heat.json), the
[joint synthesis](../research/studies/joint_synthesis/README.md#6-ledgers)'s night half, the fleet-light model and the
[first screen](README.md#heat-and-energy). Numbers marked *screen* come from that product. Literature values are cited
(Author year) and listed, with how far each was read, in [computing_sources.json](computing_sources.json). "Today"
means the NVIDIA H100's 1.4×10¹² FLOP per joule at 16 bits.

## Findings

1. **The light is far larger than any computing the array can build.** The night half intercepts 851 PW that the Moon
   never needs. A thousand terawatts of computing at 20–30% conversion takes 0.4–0.6% of it, and 10,000 TW takes
   4–6%. The light the night half's films reflect, about 14 PW, would itself power 2,800–4,200 TW (*screen*).
2. **Mass binds first.** A terawatt of computing needs 7–34 Mt in orbit: radiators of 2.5–11.6 Mt at 330 K (the array
   study), collectors of 1–10 Mt at the shield study's 100–1,000 W/kg, and processors of 3.5–12.8 Mt. The processors
   range from a 100 t launch container holding about 40 MW, with 1 kg of shielding per kW (Feilden et al. 2024), to a
   DGX H100 server (NVIDIA 2024). So 1,000 TW weighs 7–34 Gt, 15–74% of the fleet's 46 Gt. Renewing its processors on
   a GPU's five-year life (Ho et al. 2023) takes 0.7–2.6 Gt a year, the scale of the film plant's 2.3–4.6 Gt. At
   10,000 TW the computing outweighs the fleet 1.5–7.4 times (*screen*).
3. **Heat placement binds next, and distance answers it.** Orbital computing releases its own heat and its
   collectors', 2.9–4.4 W per watt computed, and from the ring radius 0.19% of it reaches the Moon. From there
   460–1,100 TW of computing warms the Moon by 0.1 K and 4,600–11,000 TW by 1 K. At 10,000 TW the Moon receives
   1.5–2.2 W/m², 26–39% of the fleet's glow at its low end. A face turned to the Moon doubles the share and a radiator
   edge-on to it cuts it. The Earth–Moon L1 and L2 points cut it 9.5 times, 100,000 km 25 times and the Earth–Sun hubs
   5,600 times. A terawatt used on the Moon adds 0.016–0.026 K, so the Moon keeps the computing that must be close: a
   hundred times today's 5.8 W of data centre per person (IEA 2025; UN 2024) is 0.58 TW per billion people,
   0.01–0.015 K (*screen*).
4. **A terawatt today is 1.4×10²⁴ FLOP/s, 94 times the world's AI computing.** The world held about 15 million
   H100-equivalents in January 2026, drawing over 10 GW (Epoch AI 2026), and its data centres drew 47 GW in 2024
   (IEA 2025) (*screen*).
5. **Energy per operation has about 200-fold to go in CMOS and 4,700–470,000-fold to the irreversible floor.** A
   16-bit operation costs 2.3×10⁸ times kT ln 2 at 330 K today. CMOS tops out near 2.9×10¹⁴ FLOP/J (Ho et al. 2023),
   12–21 years away at past doubling times. With Ho et al.'s 480–48,000 switches per operation, each erasing at least
   kT ln 2 (Landauer 1961; Bérut et al. 2012), the floor at 330 K is 6.6×10¹⁵–6.6×10¹⁷ FLOP/J. Reversible, adiabatic
   logic goes below it by running slower (Bennett 1973; Frank et al. 2020). Over the first thousand years the planning
   range runs from the CMOS limit to the floor (*screen*). Thermodynamic computers raise results per joule on
   sampling, by a modelled 10–10,000 times
   ([array_industry.md](../research/studies/solar_shield_array/array_industry.md)). They suit the fleet's estimation
   and scheduling and generative services.
6. **Radiators should run hot.** At the floor, cooling a radiator buys operations per joule as 1/T and loses operations
   per square metre as T³. At 150 K a terawatt needs 19,400 km² of radiator, each square metre doing 9% of the computing
   it does at 330 K; at 500 K it needs 157 km², each doing 3.5 times as much (*screen*). With light to spare and mass to
   save, the array is best served by radiators as hot as the logic allows.
7. **The array's control and the Moon's management take megawatts to terawatts down to 10 m.** A 1 TFLOP/s spaceflight
   processor on each of 28 million tiles (Microchip 2024) draws 20 MW today; the momentum store loses 4–155 GW. A
   coupled weather and climate twin of the whole Moon, run as a forecast service of 50 members at 20 times real time,
   takes 6–12 MW at 1 km, 6–12 GW at 100 m and 9–15 TW at 10 m. These scale COSMO's measured 596 MWh per simulated
   year (Schulthess et al. 2019) (*screen*).
8. **Micromanagement at metre scale fills the array, and following cells one by one exceeds it.** At 1 m the twin takes
   38–44 TW in real time and 38,000–44,000 TW as a forecast service, 15–26% of the night half's light, or 180–210 TW at
   the CMOS limit. A gram of soil holds up to 10¹⁰ cells (Torsvik & Øvreås 2002). At the 4D minimal cell's 300 A100
   GPU-hours per 105-minute cycle (Thornburg et al. 2026), following each of them takes 350 TW. Earth's 4–6×10³⁰
   prokaryotes (Whitman et al. 1998) would take 360–550 million Suns (*screen*).
9. **Science of today's kind takes megawatts, and molecular-scale biology terawatt-years.** Correlating FarView's
   100,000 dipoles (Polidan et al. 2024) takes 1 MW and the SKA's science processing 0.2 MW. Lattice QCD's largest
   proposed ensemble (Boyle et al. 2022) is 31 s of a terawatt. An all-atom minimal cell through one cycle,
   3.8×10³¹ FLOP, is 312 days at 1 TW (*screen*).
10. **AI for 20–30 billion people matches the array today and shrinks to terawatts at the limits.** At 1–100
    brain-equivalents of 10¹⁵ FLOP/s per person (Carlsmith 2020) it takes 14–2,100 TW today, up to 1.3% of the night
    half's light, and 0.07–10 TW at the CMOS limit. Counted as 100–10,000 tokens a second from 10¹²-parameter models
    (Kaplan et al. 2020), it takes 3–430 TW; a person uses about 0.06 tokens a second today (Alphabet 2025; UN 2024).
    One training run at the latency wall, 3×10³⁰–10³² FLOP (Sevilla et al. 2024), takes 0.25–8 TW over 100 days.
    A brain-equivalent draws 714 W today and 3.4 W at the CMOS limit, so the [population](population.md) note's
    0.1–1 kW of computing a person is 0.14–1.4 brain-equivalents today and 29–290 at the limit (*screen*).
11. **The far-out tasks span the array and pass it.** A digital population equal to ours, at 10¹⁵–10¹⁸ FLOP/s a mind,
    takes 14–21,000 TW today. A trillion ems at 1,000 times human speed (Hanson 2016) take 2.8–4,200 times the night
    half's light, and 1.5–150,000 TW at the floor. One brain emulated at the metabolome level (Sandberg & Bostrom 2008)
    takes 7 TW. An ancestor simulation (Bostrom 2003) runs in 8 days to 23 years on 1,000 TW today. Re-running the
    computation of evolution's nervous systems, about 10⁴¹ FLOP (Cotra 2020), takes 10,800 years on 1,000 TW at the
    CMOS limit. A Matrioshka brain (Bradbury 1999) draws 3.5×10⁸ times the night half's light (*screen*).
12. **Latency divides the work.** The 0.13 s round trip to the ring radius fits inside the 0.21 s people leave between
    turns in conversation (Stivers et al. 2009), so the Moon's interactive AI, management and forecasts can run in the
    array. Millisecond control stays on the Moon, and Earth's 2.56 s round trip suits training, science and batch
    services (*screen*).

## Supply

Collectors fly on platforms of their own in light the Moon never needs, and radiators face outward or stand edge-on to
the Moon. A platform for one terawatt spans 57–67 km. Ranges pair 30% conversion and the light hardware with 20% and the
heavy hardware.

| Computing power | 1 TW | 10 TW | 100 TW | 1,000 TW | 10,000 TW |
|---|---:|---:|---:|---:|---:|
| FLOP/s today | 1.4×10²⁴ | 1.4×10²⁵ | 1.4×10²⁶ | 1.4×10²⁷ | 1.4×10²⁸ |
| FLOP/s at the CMOS limit | 2.9×10²⁶ | 2.9×10²⁷ | 2.9×10²⁸ | 2.9×10²⁹ | 2.9×10³⁰ |
| Share of the night half's 851 PW | 0.0004–0.0006% | 0.004–0.006% | 0.04–0.06% | 0.4–0.6% | 3.9–5.9% |
| Collectors and radiators | 3,300–4,500 km² | 33,000–45,000 km² | 0.33–0.45 million km² | 3.3–4.5 million km² | 33–45 million km² |
| Mass in orbit | 7–34 Mt | 70–340 Mt | 0.7–3.4 Gt | 7–34 Gt | 70–340 Gt |
| Share of the fleet's 46 Gt | 0.015–0.07% | 0.15–0.7% | 1.5–7.4% | 15–74% | 150–740% |
| Processor renewal at a five-year life | 0.7–2.6 Mt a year | 7–26 Mt a year | 0.07–0.26 Gt a year | 0.7–2.6 Gt a year | 7–26 Gt a year |
| Heat on the Moon from the ring radius | 1.5–2.2×10⁻⁴ W/m² | 0.0015–0.0022 W/m² | 0.015–0.022 W/m² | 0.15–0.22 W/m² | 1.5–2.2 W/m² |
| Its warming | 0.0001–0.0002 K | 0.001–0.002 K | 0.009–0.022 K | 0.09–0.22 K | 0.9–2.2 K |
| The same power used on the Moon | 0.016–0.026 K | 0.16–0.26 K | 1.6–2.6 K | 16–26 K | 160–260 K |
| Brain-equivalents of 10¹⁵ FLOP/s, today | 1.4×10⁹ | 1.4×10¹⁰ | 1.4×10¹¹ | 1.4×10¹² | 1.4×10¹³ |

## Demand against supply

**Continuous services.** Power is at today's efficiency unless named. Warming is for the computing and its collectors
at the ring radius. Weather and cell rows carry measured energy on P100 and A100 GPUs, scaled to the H100 by 32-bit
FLOP per joule.

| Task | Computing, FLOP/s | Today | At the CMOS limit | Share of the night half's light, today | Warming from the ring radius | Runs |
|---|---:|---:|---:|---:|---:|---|
| The fleet's control: a 1 TFLOP/s processor on each tile (Microchip 2024) | 2.8×10¹⁹ | 20 MW | 0.1 MW | about 10⁻⁸% | below 10⁻⁸ K | on the tiles |
| The Moon's weather and climate twin, 1 km, forecast service (Schulthess et al. 2019) | 0.8–1.6×10¹⁹ | 6–12 MW | 0.03–0.06 MW | below 10⁻⁸% | below 10⁻⁸ K | array |
| The same at 100 m | 0.9–1.7×10²² | 6–12 GW | 29–57 MW | below 10⁻⁵% | below 10⁻⁵ K | array |
| The same at 10 m | 1.3–2.1×10²⁵ | 9–15 TW | 43–70 GW | 0.004–0.009% | 0.001–0.003 K | array |
| The same at 1 m, real time | 5.3–6.1×10²⁵ | 38–44 TW | 0.18–0.21 TW | 0.015–0.026% | 0.003–0.01 K | array |
| The same at 1 m, forecast service | 5.3–6.1×10²⁸ | 38,000–44,000 TW | 180–210 TW | 15–26% | 3.5–9.6 K | far from the Moon |
| Cell by cell, one gram of soil (Torsvik & Øvreås 2002; Thornburg et al. 2026) | 4.9×10²⁶ | 350 TW | 1.7 TW | 0.14–0.21% | 0.03–0.08 K | anywhere |
| Cell by cell, Earth's prokaryotes (Whitman et al. 1998) | 2–3×10⁴⁷ | 360–550 million Suns | | | | beyond |
| Radio correlation, FarView's 100,000 dipoles of one polarization over 5–40 MHz (NASA, FarView) | 1.4×10¹⁸ | 1 MW | 5 kW | below 10⁻⁹% | below 10⁻⁹ K | by the telescope |
| Ten million elements, full / by FFT (Tegmark & Zaldarriaga 2009) | 1.4×10²² / 4×10¹⁶ | 10 GW / 29 kW | 48 MW / 0.1 kW | below 10⁻⁵% | below 10⁻⁵ K | by the telescope |
| The SKA's science data processors (SKAO) | 2.7×10¹⁷ | 0.19 MW | 0.9 kW | about 10⁻⁸% | below 10⁻¹⁰ K | anywhere |
| Folding@home at its 2020 peak (Zimmerman et al. 2021) | 1.0×10¹⁸ | 0.72 MW | 3.4 kW | below 10⁻⁷% | below 10⁻⁹ K | anywhere |
| The world's AI computing, January 2026 (Epoch AI 2026) | 1.5×10²² | 11 GW | 51 MW | below 10⁻⁵% | below 10⁻⁵ K | Earth |
| One training run at the latency wall over 100 days (Sevilla et al. 2024) | 3.5×10²³–1.2×10²⁵ | 0.25–8.3 TW | 1.2–39 GW | 0.0001–0.005% | 0.00002–0.002 K | anywhere, on one platform |
| AI services for 20–30 billion people, 1–100 brain-equivalents each (Carlsmith 2020) | 2×10²⁵–3×10²⁷ | 14–2,100 TW | 0.07–10 TW | 0.006–1.3% | 0.001–0.47 K | the array, for the Moon's people |
| A digital population equal to ours, 10¹⁵–10¹⁸ FLOP/s a mind (Carlsmith 2020; Sandberg & Bostrom 2008) | 2×10²⁵–3×10²⁸ | 14–21,000 TW | 0.07–100 TW | 0.006–13% | 0.001–4.7 K | array, or farther out |
| A trillion ems at 1,000 times human speed (Hanson 2016) | 10³⁰–10³³ | 7×10⁵–7×10⁸ TW | 3,400–3.4×10⁶ TW | 2.8–4,200 night halves | | beyond the array |
| One emulated brain, spiking network to electrophysiology (Sandberg & Bostrom 2008) | 10¹⁸–10²² | 0.71 MW–7.1 GW | 3.4 kW–34 MW | | | anywhere |
| One emulated brain, metabolome | 10²⁵ | 7.1 TW | 34 GW | 0.003–0.004% | 0.0007–0.002 K | anywhere |
| One emulated brain, states of protein complexes | 10²⁷ | 710 TW | 3.4 TW | 0.28–0.42% | 0.07–0.16 K | anywhere |
| One emulated brain, distribution of complexes | 10³⁰ | 7.1×10⁵ TW | 3,400 TW | 2.8–4.2 night halves | | beyond the array |
| One emulated brain, single molecules | 10⁴³ | 18,700 Suns | | | | beyond |
| A planetary-mass computer, the earlier "Jupiter brain" (Bostrom 2003; Bradbury 1999) | 10⁴² | 1,900 Suns | | | | beyond |
| A Matrioshka brain (Bradbury 1999) | 3×10⁴² | on 0.8 Suns at its own 10¹⁶ operations per joule | | | | beyond |

**Campaigns.** Times are for the whole task at the stated power.

| Campaign | Operations | 1 TW today | 1,000 TW today | 1,000 TW at the CMOS limit | 10,000 TW at the floor |
|---|---:|---:|---:|---:|---:|
| AlphaFold's training, 128 TPU v3 cores for about 11 days (Jumper et al. 2021) | 7.5×10²¹ | 5 ms | | | |
| Frontier-E, 4 trillion particles (Frontiere et al. 2025) | 3.1×10²³ | 0.2 s | | | |
| The largest proposed lattice QCD ensemble (Boyle et al. 2022) | 4.3×10²⁵ | 31 s | | | |
| The largest AI training run of 2025, Grok 3 (Epoch AI 2025b) | 4.6×10²⁶ | 5.5 min | 0.3 s | | |
| An all-atom minimal cell through one cycle (Stevens et al. 2023; Sandberg & Bostrom 2008) | 3.8×10³¹ | 312 days | 7.5 hours | 2 min | 0.006–0.6 s |
| An ancestor simulation (Bostrom 2003) | 10³³–10³⁶ | 23–23,000 years | 8 days–23 years | 1 hour–39 days | 0.15 s–4 hours |
| The evolution anchor (Cotra 2020) | 10⁴¹ | 2.3 billion years | 2.3 million years | 10,800 years | 0.5–48 years |
| Counting through 2²⁵⁶ keys at kT ln 2 a step at 2.7 K, after Schneier 1996 | 1.2×10⁷⁷ | 3×10⁵⁴ J, the Sun's output for 2.5×10²⁰ years | | | |

The far end has its own limits. A kilogram computes at most 5.4×10⁵⁰ operations a second (Margolus & Levitin 1998;
Lloyd 2000) and 1.36×10⁵⁰ bits a second by Bremermann's bound (Bremermann 1962). Ordinary matter computing with its
nuclear spins manages about 10⁴⁰ a second per kilogram (Lloyd 2000). The universe has performed at most 10¹²⁰
operations (Lloyd 2002), and waiting for a colder universe multiplies what can be computed by 10³⁰ (Sandberg,
Armstrong & Ćirković 2017). A year of the Sun's output at 3.2 K steps a counter through 2¹⁸⁷ states (Schneier 1996).
Against these, today's computing in orbit gives 4×10¹³–2×10¹⁴ FLOP/s per kilogram of radiators, collectors and
processors (*screen*). Bradbury moved from planet-sized "Jupiter brains" to nested shells because a solid computer
that size cannot shed a star's power without melting and spends much of its energy moving coolant (Bradbury 1999): the
radiator problem of finding 6 at the scale of a star.

## Where the work runs

- **On the tiles and on the Moon:** control loops that close in milliseconds, such as tile keeping, sky-boat
  separation, grid protection and robots.
- **In the array, within the 0.13 s round trip:** the Moon's AI services and digital minds that meet its people, its
  weather and ecosystem twins, and the fleet's estimation and scheduling. Learned forecast models run a 10-day global
  forecast in under a minute (Lam et al. 2023), and the physics twin trains and checks them.
- **Farther out or anywhere:** training, science campaigns, populations of minds that live in the computer, and
  computing beyond about 1,000 TW, whose heat from the ring radius would pass 0.1 K on the Moon; from 100,000 km the
  same heat lands 25 times weaker.
- **For Earth, at a 2.56 s round trip:** batch services, models and science. Text for Earth's 10.3 billion people at
  100–10,000 tokens a second each needs 33 Tb/s–3.3 Pb/s, 165–16,500 links of TBIRD's 200 Gb/s (Schieler et al.
  2023). Video at 10 Mb/s each needs about 515,000 of them (*screen*). A laser link has spanned the Earth–Moon
  distance at 622 Mb/s (Boroson & Robinson 2014).
- **By the telescope:** a low-frequency radio array, which the Open Moon's ionosphere moves into deep space behind a
  shield of its own ([conservation](../research/studies/conservation/README.md)), with its correlator beside it.
  The correlation runs at reduced precision on tensor cores (Romein 2021).

## Open questions, in order of what they settle

1. **How much computing the civilization wants.** This means AI services per person, whether digital minds count among
   the 20–30 billion, and what share of the array's light the people choose to turn into computing. It sets the scale
   between 1 and 10,000 TW, and with it the mass and the heat. It belongs with the [population](population.md) work and
   the human paper.
2. **Processor mass and life at the ring radius.** These are kilograms per kilowatt in orbit, wear in the radiation
   outside Earth's magnetosphere for most of each month, and the renewal they force. At 1,000 TW renewal matches the
   film plant. The answer settles whether that scale can be built beside the fleet.
3. **Where large computing sits.** The choices are the ring radius, the Earth–Moon L1 and L2 points, 100,000 km or the
   Earth–Sun hubs, with radiators edge-on to the Moon and the collectors' reflections kept off its night as the
   fleet's are. The planned climate needs a share of its heat budget set for orbital computing: at the ring radius
   0.1 K comes at 460–1,100 TW.
4. **The path in energy per operation.** It runs from CMOS to adiabatic and reversible logic, through hot radiators and
   logic that tolerates them, and through thermodynamic co-processors judged by results per joule on matched tasks.
   It settles how much computing each tonne buys.
5. **The management twin's design.** It needs resolution chosen by place, metre scale over cities, farms and flight
   corridors, with learned emulators trained on the physics. That decides whether management stays at gigawatts or
   reaches the array's scale.
6. **Services to Earth, and the array as a prototype of Earth's own** (open in the register). Bandwidth and the 2.56 s
   round trip decide what can be exported.
7. **Collectors as dimmers.** Collectors in the aperture take light bound for the Moon, so they dim it while feeding
   computing. The dimmer's 586 TW of optical power would give 117–176 TW
   ([integrated comparison](../research/studies/solar_shield_array/integrated_comparison.md#electricity)). This belongs
   with the shield's work on the glow.

## How the numbers were made

`python -m provisioning.computing` rewrites the product in under a second, and
`python -m pytest provisioning/tests/test_computing.py` checks that it binds its inputs and that its arithmetic holds.
It reads `array_heat.json` for collectors, radiators, heat placement and the fleet, the joint ledger for the night half,
the fleet-light product for the night half's reflections, and the first screen for the warming per W/m² and the latency.
The screen is closed-form arithmetic. The weather twin scales COSMO's measured energy per cell and step to the Moon's
area:
- each tenfold refinement costs at least a thousandfold, a hundredfold in columns and tenfold in steps;
- 180–360 levels, with the lowest kilometre refined to the grid;
- a coupling factor of 1.2 for land, ocean, waves and ice.

These choices are assumptions for scale:
- the per-person AI use;
- the weather twin's levels and its forecast service of 50 members at 20 times real time;
- the 20-state, 1,000-member fleet filter;
- the 100-day training run.

The bounds apply the 16-bit ratios to every task. Ideal energy grows as the square of precision (Ho et al. 2023), so
from the H100's own units the gain to the CMOS limit is 210 at 16 bits, 770 at 32 bits and 190–380 at 64 bits
(*screen*).
