# Calculation provenance and assumptions

## 1. Kardashev / Sagan interpolation

The current repository itself contains the formula

`K = (log10(P) - 6) / 10`, with `P` in watts.

For K=1.17:

`P = 10^(10*1.17 + 6) = 5.011872336e17 W = 501.187 PW`.

**Status:** repository-encoded calculation. In the recovered September 19 discussion, this was a civilization-scale controllable-power scenario, not a measured project load.

## 2. D-T fusion fuel scale

Recomputed here from 17.6 MeV per D-T reaction and the reacting D+T mass.

At 501.187 PW continuously for one Julian year:

- ideal fusion-energy conversion: ~4.685×10^10 kg/year reacting D+T;
- at 50% net conversion: ~9.370×10^10 kg/year reacting D+T.

The historical discussion rounded these to ~46 Mt/year and ~90 Mt/year.

**Status:** recovered design-history number, independently recomputed in this bundle. This is reacting fuel, excluding breeding, processing, inventory, downtime, plant losses beyond the chosen net factor, and non-electric uses.

## 3. Fission fuel scale

Using 200 MeV per fission and 235 atomic mass units as a simple scale gives ~1.93×10^11 kg/year of actually fissioned heavy nuclei at 501.187 PW. The historical discussion rounded this to ~2×10^8 tonnes/year.

**Status:** recovered comparison scale, independently recomputed. This is fissioned mass, not mined/enriched fuel throughput.

## 4. Solar collection area

The September baseline uses a convenient net solar-electric output of 300 W/m² near 1 AU.

At that net output:

- 10 PW → 33.33 million km²;
- 100 PW → 333.33 million km²;
- 500 PW → 1.667 billion km².

At 0.1 AU, inverse-square sunlight is ~100× the 1-AU value. Holding the same conversion fraction and ignoring thermal/material limits would reduce a 500-PW collector area to ~16.67 million km². This 0.1-AU figure is an ideal geometric scaling, not a collector design.

The September 19 discussion used this to motivate near-Sun generation plus power transmission/beaming as one mature-system option.

## 5. Solar-luminosity share

Taking solar luminosity as 3.828×10^26 W, 501.187 PW is ~1.31×10^-9 of the Sun's luminosity.

This shows that the scenario is large industrially while remaining an extremely small fraction of stellar output.

## 6. Radiator scale

For one-sided deep-space radiators with emissivity 0.9:

`A = P / (epsilon * sigma * T^4)`.

For 501.187 PW of rejected heat:

- 350 K → 654 million km²;
- 600 K → 75.8 million km²;
- 1000 K → 9.82 million km²;
- 1500 K → 1.94 million km².

The historical discussion rounded these to ~650, ~75, ~10, and ~2 million km². The September baseline independently gives the same scaling at 10 PW: ~13 million km² at 350 K and ~1.5 million km² at 600 K.

These are ideal view-to-deep-space areas. Real systems incur view-factor, sunlight, structure, pumping, temperature-gradient, conversion, and reliability penalties.

## 7. Construction flow and pipeline inventory

The September baseline gives ~2.67×10^8 kg/s nitrogen flow for a 300-year N2 build and ~3.32×10^8 kg/s total atmospheric flow. The recovered September 19 discussion used a representative bulk flow of ~2.7×10^8 kg/s.

At 2.7×10^8 kg/s and a six-year transit time:

`M_pipeline ≈ 5.11×10^16 kg`.

The September report independently describes a six-year N2 route as placing roughly 5×10^16 kg in transit.

## 8. Packet-rate examples

At 2.7×10^8 kg/s:

- 10^6 kg packets → 270/s;
- 10^8 kg packets → 2.7/s;
- 10^9 kg packets → one every 3.70 s;
- 10^12 kg packets → one every 61.7 min;
- 10^13 kg packets → one every 10.3 h;
- 10^14 kg packets → one every 4.29 d;
- 10^15 kg packets → one every 42.9 d.

The larger values were throughput illustrations, not preferred safety choices. The recovered safety discussion favored roughly 10^6–10^8 kg packets for mature high-density traffic because individual consequences remain more bounded.

## 9. Packet kinetic-energy examples

At 10 km/s:

- 10^6 kg → 5×10^13 J (~12 kt TNT);
- 10^8 kg → 5×10^15 J (~1.20 Mt TNT);
- 10^9 kg → 5×10^16 J (~11.95 Mt TNT);
- 10^12 kg → 5×10^19 J (~11,950 Mt TNT).

This is why the historical safety architecture treated packet size, intercept geometry, remote capture and passive-failure trajectories as first-class design variables.

## 10. Final braking heat

The September baseline explicitly calculates that at the ~300-year nitrogen flow:

- dissipating 3 km/s locally → ~1.2 PW, ~32 W/m² globally averaged;
- dissipating 10 km/s locally → ~13.4 PW, ~352 W/m² globally averaged.

The recovered architecture therefore moves major braking/capture away from the final atmosphere and meters material inward after capture.
