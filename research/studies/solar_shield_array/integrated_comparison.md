# Integrated comparison of the shield candidates

On 6 October 2026 the author set the acceptance for the whole shield system
([research/decisions.md](../../decisions.md)). It has to meet the shield
requirements, give reasonable net-positive electricity for habitat, balance
its overall resource requirements, and be stable and safe. Every mass flow is
judged by its direction as well as its size: whatever a design expels must
stay away from the Moon's cold, dense exobase. This comparison replaces the
decision gate of the [relative-orbit plan](relative_orbits.md), which tested the
ring fleet's recurrence cost against a 23.835 TW target. The 23.835 TW figure
and the 100 TW civilization case remain as comparison references, with higher
demand as sensitivities.

A candidate has to pass every gate: the shield requirements, stability, safety
and net-positive electricity for habitat after its own loads. Candidates that
pass are weighed on resources: hardware, fresh feedstock, recycled throughput
and delivery. [integrated_ledger.py](integrated_ledger.py) gathers what the
products establish, adds the arithmetic that follows from stated assumptions,
marks what no model yet supplies and writes
[results/integrated_ledger.json](results/integrated_ledger.json).
[exhaust_isolation.py](exhaust_isolation.py) follows the held screen's exhaust
and writes [results/exhaust_isolation.json](results/exhaust_isolation.json).

## The candidates

**The held screen with a zoned aperture** holds independent tiles on the
selected moving trajectory, 67,700–105,500 km sunward of the Moon over each
month and 90,700 km on average
([the plan's zoned-aperture section](relative_orbits.md#a-zoned-aperture-for-the-held-screen)).
The window, 2,245 km in radius, carries the climate stack at 26 g/m²; the
annulus out to 7,454 km carries a 5 g/m² UV film that passes visible light.
Its area is about 1.75 million panel-equivalents of 10 km (1.78 million of
9.89 km clear aperture), with no packing layout yet.

**The ring fleet** flies nested ring bundles on natural orbits about the Moon,
each ring one orbit of shingled tiles, kept with photon forces. No full-fleet
run exists. Its inventory below follows from the bundle's geometry: one ring
for every 9.29 km height step across the aperture, and one tile for every
8.5 km of ring, with every tile of a ring that crosses the window carrying the
window's areal mass. The count comes before any oversizing for the strip
pattern's annual motion.

## The gates

Each state is an assessment of the cited evidence for the author: met,
conditional (met if a stated condition holds), open (no model yet decides it)
or not met. The ledger product carries the evidence behind each.

| Gate | Held screen, zoned | Ring fleet |
|---|---|---|
| UV transmission (O1, O8) | Open: the 5 g/m² annulus film is unproven, because the 24.8–120 nm opacity proof rests on 22 g/m² of silica; continuous area with no seams | Open: projected overlaps stay positive in the kept bundle; receiver-ray coverage, handovers and the annulus film are unproven |
| Protected radius (O2) | Met: the 7,454 km aperture covers four lunar radii with the finite Sun | Open: the strip pattern's sunward crossing drifts 2,000–8,500 km over a year at 15,000 km until the frozen common orbit is found |
| Window spectrum (E2, O4) | Conditional on the 26 g/m² window and the dimmer's form | The same, on the rings that cross the window |
| Earth's shadow and rejected light (O5, O7) | Open: the optical bound lets redirected light leave toward the Moon; the shadow on Earth is unmapped | Open: night-side tiles reflect back past the Moon, with light the screen has already filtered (to check); the shadow on Earth is unmapped |
| Stability | Conditional on continuous thrust (75.7 TW installed); an unpowered tile drifts about 2,700 km in its first day | Open: the kept bundle's edges first come within 150 m at day 8.1; the interior stays clear for 23 days in a one-tile-per-ring model |
| Outflows (S7, proposed) | Not met at the 45° cant, by the direct plumes and by the slow gas alone (below) | Open: photon keeping releases nothing, pending its check; electric keeping would release 0.4 t/s or more at or inside the planned magnetosphere |
| Holding cheaper than the atmosphere it saves (S6) | Not met: 51,000 kg/s against 300–80,000 kg/s of unshielded loss | Conditional on photon keeping |
| Net-positive electricity for habitat | Open: its own load is 32.8 TW mean | Open: its own load is the attitude energy, not yet computed |

## Where the held screen's exhaust goes

Holding a screen upstream of the Moon needs a sunward push that sunlight
cannot give, so all of its reaction mass leaves into the Moon's hemisphere. At
the zoned design point that is 51,000 kg/s carrying 23 TW of jet power. The
proposed requirement S7
([requirements.md](../protection_architecture/requirements.md)) keeps
outflows out of the escape region. The September report's energy diagnostic
([protection/report.md](../../../protection/report.md), section 8) lets a
tenth of the deposited energy do escape work against 0.94 MJ/kg, so each kg/s
of loss budget allows 9.4 MW deposited there: 4.1×10⁻⁷ of the jet power for a
1 kg/s budget and 4.1×10⁻⁶ for 10 kg/s.

**The direct path.** Each tile's propulsive need is split between two
thrusters canted from the resultant, in the plane that keeps both plumes off
the Moon's side. Sampled daily over the year at 96 nodes, and weighted by each
node's propellant flow, this is the share of fast ions that leaves on straight
lines into the protected sphere, whose half-angle from the screen is
3.8–5.9°:

| Cant | Propellant and power | Gridded ion (NEXT) | Hall (BPT-4000) |
|---|---:|---:|---:|
| 45° | 1 | 6.2×10⁻⁶ | 1.6×10⁻⁴ |
| 60° | 1.41 | 1.6×10⁻⁶ | 6.7×10⁻⁵ |
| 75° | 2.73 | 4.6×10⁻⁷ | 2.6×10⁻⁵ |

At 45° the worst node and date reach 7×10⁻⁵ for the gridded ion plume and
1.3×10⁻³ for the Hall plume. Canting each pair in the plane through the Moon
instead points one plume of a pair into the protected sphere at times and
sends 1–2% of the fast ions there. Hall plumes miss the allowance for 10 kg/s
by 6–40 times at every cant. A gridded ion plume meets it from about 60° and
comes to the 1 kg/s allowance at 75°, at 1.4 and 2.7 times the 45° design's
propellant and power. The plume profiles are read from published figures
(Young and others, 2019, for NEXT; Goebel and Katz, 2008, for the BPT-4000) to
about ±30%, and these shares count only fast ions on straight lines.

**The slow gas.** Flight-class Hall thrusters reach about 95% mass
utilization, so at least about 5% of the propellant leaves unionized, at
thermal speed; at 5–15% that is 2,550–7,650 kg/s here. Test particles left
every node every two days for a month, with each node's own velocity, a cosine
law about each thruster axis and the effusive speeds of oxygen molecules or
xenon atoms at 500 K, and flew for 60 days under DE440s point-mass gravity:

| Share of the unionized gas | Oxygen molecules | Xenon |
|---|---:|---:|
| Enters the protected sphere within 10 days | 1.2% | 4.6% |
| Enters it within 60 days | 2.2% | 8.1% |
| Passes within 10 Earth radii of Earth | 12.4% | 3.6% |
| Leaves beyond 3 million km | 30.4% | 11.5% |
| Still in flight at 60 days | 55.1% | 76.7% |
| Delivered, with a 30-day neutral lifetime | 1.5% | 5.5% |
| Delivered, kg/s | 38–114 | 140–420 |
| Kinetic power on arrival, MW | 30–91 | 97–291 |

The gas arrives at about 1.2 km/s at four lunar radii and gains speed as it
falls further in. Halving the integration step changed the fate of none of the
first launch's 2,304 particles. Neutral lifetimes of 10 to 100 days give
delivered shares of 1.0–1.9% for oxygen and 3.5–7.1% for xenon. The slow gas
alone deposits 3–30 times the 1 kg/s allowance, and at
the upper end three times the 10 kg/s allowance. The September report
proposed capturing the neutral gas locally; it would have to stop nearly all of
it. A further 4–12% of the gas passes near Earth, which S7 also asks to trace.

**The charged exhaust.** The other 85–95% of the flow is ions. 51 t/s is about
10⁵ times the solar wind's own mass flux through the aperture, so the screen
acts as a comet upstream of the Moon. The wind loads, slows and carries the
ions downstream over the Moon, energising picked-up oxygen to tens of keV, and
the [loss response](../../../atmosphere/loss_response/README.md) credits each
returning pickup ion with 1–10 sputtered molecules. No model here follows that
plasma. It stays open, and it adds to both paths above.

**What follows.** On present evidence the held screen's exhaust breaks the
proposed S7 on both paths computed here, before the plasma path is counted.
Meeting it would take gridded ion thrusters canted 60–75°, capture of nearly
all unionized propellant, and a demonstration that the solar wind returns
almost none of the charged exhaust. The held screen also fails S6: its
51,000 kg/s of propellant exceeds the atmosphere's unshielded loss of
300–80,000 kg/s everywhere but the top of that range.

## Resources

Optical masses use 26 g/m² wherever a tile serves the window and 5 g/m²
elsewhere. Film replacement follows O8's 10–20 years. Worn film can be
recovered and remade; propellant is expelled for good.

| | Held screen, zoned | Ring fleet at 15,000 km | at 19,000 km | at 20,000 km |
|---|---:|---:|---:|---:|
| Tiles of 10 km | 1.75 million (equivalents) | 16.8 million | 21.3 million | 22.4 million |
| Tile area over the aperture's | 1 | 10.8 | 13.7 | 14.4 |
| Optical mass | 1.21 Gt | 17.4 Gt | 22.2 Gt | 23.5 Gt |
| Propellant | 51,000 kg/s, 1.6 Gt a year | none with photon keeping | none | none |
| Film replacement | 0.06–0.12 Gt a year | 0.9–1.7 Gt a year | 1.1–2.2 Gt a year | 1.2–2.3 Gt a year |
| Thrust power | 32.8 TW mean, 75.7 TW installed | none; attitude energy open | none | none |

The held screen burns 1.08 kg of propellant per held kilogram each year, about
2% of the Moon's mass over a billion years, and every kilogram added to it,
whether collector or habitat, burns the same. The ring fleet has about ten to
thirteen times the held screen's tiles to make and replace, as recycled
throughput. Done by 30 km/s electric thrust, its keeping and steering would
release 0.4 t/s and 0.06–5.7 t/s of exhaust in lunar orbit, which is why
photon keeping was adopted. The kept run's keeping demand peaks at 4.6% of the
sail acceleration that the filter's redirected band gives a 5 g/m² tile and
24% of a 26 g/m² tile's; the steering of the steady drift needs 10⁻⁷–10⁻⁵ m/s².
Whether a sail force that always points away from the Sun can follow that
demand around each orbit is the next check.

## Electricity

The light is far larger than any demand case: the aperture intercepts about
240 PW. The dimmer's 586 TW of optical power would give 117–176 TW at 20–30%
conversion, and opaque collectors on 1% of the annulus 432–648 TW, at the cost
of about 1% dimming on Earth for a few hours at eclipse-season new moons unless
they turn edge-on then. What collection costs depends on where the collectors
ride.

| Demand | Held screen: generation with its own load | Collectors at 20–30% | Share of the annulus | Propellant if they ride the held screen |
|---:|---:|---:|---:|---:|
| 100 TW | 133 TW | 0.33–0.49 million km² | 0.2–0.3% | 11–17 kg/s per g/m², plus 70–104 kg/s for the light they absorb |
| 300 TW | 333 TW | 0.82–1.22 million km² | 0.5–0.8% | 28–42 kg/s per g/m², plus 174–262 kg/s |
| 1,000 TW | 1,033 TW | 2.5–3.8 million km² | 1.6–2.4% | 87–130 kg/s per g/m², plus 541–812 kg/s |

At 30 g/m² of collector that is 0.4–0.6 t/s more propellant at 100 TW and
3.1–4.7 t/s at 1,000 TW, all of it into the Moon's hemisphere. The ring fleet
needs the demand's collectors alone, 0.24–0.37 million km² for 100 TW, on its
tiles or flying free, with no holding cost; its own attitude load is still
open. Device yield, transfer losses, storage through eclipses and reserves
decide the power delivered to habitat and industry, and none is modelled yet.

## Inhabited infrastructure

Habitats need their own structure, radiation shielding and life support,
none budgeted here. They belong on natural orbits in either design: on the held
screen a habitat would burn about its own mass in propellant every year.

## What the comparison shows so far

The held screen meets the protected radius on its geometry, but on present
evidence its exhaust breaks the proposed outflow requirement on both paths
computed here, and its propellant fails S6. The ring fleet passes no gate yet
and fails none: its open items are photon keeping, receiver-ray coverage with
independently propagated tiles, the edge rings, the frozen common orbit and its
oversizing, and the attitude energy. Both share the annulus film's unproven UV
transmission and the unmapped shadow on Earth.

## Next steps

These follow the order in the [plan](relative_orbits.md#plan):

4. Check photon keeping for the ring fleet against the kept run's demand.
5. Compute the annulus film's UV transmission in the protection domain.
6. Validate the bundle with independently propagated tiles, finite squares and
   receiver rays.
7. Redo the keeping in the regressing frame, treat the edge rings and find the
   frozen common orbit, as far as the ledger keeps the fleet in contention.

The held screen's plasma path needs a model this repository lacks; a bound
from published comet and ion-release studies comes first.

## Running

```sh
OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.exhaust_isolation
python -m research.studies.solar_shield_array.integrated_ledger
python -m pytest research/studies/solar_shield_array/test_exhaust_isolation.py research/studies/solar_shield_array/test_integrated_ledger.py
```

The exhaust study takes about 27 CPU minutes, almost all of it the slow-gas
particles, and keeps one checkpoint per launch date in
`research/runs/solar_shield_array/exhaust_isolation/`. The ledger takes a
second.

## Sources

- [protection/report.md](../../../protection/report.md), section 8: exhaust
  isolation and the energy diagnostic.
- [atmosphere/loss_response](../../../atmosphere/loss_response/README.md):
  returning pickup ions sputter 1–10 molecules each.
- J. A. Young, T. S. Matlock, M. Nakles and M. W. Crofton, *Far Field Plume
  Distribution and Divergence for NEXT: DART Mission*, AIAA SciTech Forum 2019
  ([NTRS 20190004942](https://ntrs.nasa.gov/citations/20190004942)), figure 5.
- D. M. Goebel and I. Katz, *Fundamentals of Electric Propulsion: Ion and Hall
  Thrusters*, JPL Space Science and Technology Series, 2008: chapter 7 (mass
  utilization of SPT thrusters) and chapter 8, figures 8-6 and 8-12 (plume
  current density).
- A. N. Heays, A. D. Bosman and E. F. van Dishoeck, *Photodissociation and
  photoionisation of atoms and molecules of astrophysical interest*, A&A 602,
  A105 (2017), tables 18–19.
