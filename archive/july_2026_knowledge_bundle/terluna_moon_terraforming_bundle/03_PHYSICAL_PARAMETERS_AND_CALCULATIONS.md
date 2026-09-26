# Physical Parameters and Screening Calculations

## 1. Purpose and limits

This file makes the Terluna reference design numerically explicit. The calculations are order-of-magnitude checks intended to expose the controlling variables; they are not substitutes for a general circulation model, kinetic escape model, combustion program, ecological model, or final engineering design.

The reference case uses an **80 kPa** dry surface atmosphere, **285 K** reference temperature, and **3 × 10^17 kg** water inventory. Those values are reconstructed modeling choices rather than recovered canon.

## 2. Fixed lunar quantities

| Quantity | Value | Use |
|---|---:|---|
| Mean radius | 1,737.4 km | Surface area and orbital scales |
| Surface area | 3.793 × 10^13 m² | Atmosphere, water and land inventories |
| Surface gravity | 1.624 m/s² | Pressure-column relation, lapse rate and structures |
| Escape speed | about 2.38 km/s | Atmospheric escape and transport screening |
| Vacuum circular speed near the surface | about 1.68 km/s | Space-access reference before atmosphere |
| Synodic solar cycle | 29.530588 Earth days | Sunrise-to-sunrise environmental cycle |
| Illumination half-cycle | 14.765294 Earth days | Approximate bright and dark halves near the equator |
| Solar angular rate | 0.50795 degrees/hour | Twilight and moving-terminator geometry |
| Equatorial terminator ground speed | about 4.3 m/s | Mobile ecology, weather and operations |

The Moon is synchronously rotating relative to Earth but still rotates once per sidereal month. Atmospheric dynamics therefore feel a rotation rate about one twenty-seventh of Earth’s, while the Sun crosses the sky on the 29.53-day synodic period [S01][S02].

## 3. Atmospheric mass

For a hydrostatic global atmosphere with surface pressure `P`, surface area `A`, and surface gravity `g`, the first-order total mass is:

`M_atm = P A / g`

At 80,000 Pa:

`M_atm ≈ 80,000 × 3.793 × 10^13 / 1.624 ≈ 1.87 × 10^18 kg`

### Pressure sensitivity

| Surface pressure | Atmosphere mass | Main interpretation |
|---:|---:|---|
| 5 kPa | 1.17 × 10^17 kg | Thin global experimental atmosphere |
| 10 kPa | 2.34 × 10^17 kg | Climate and escape plateau; not open-air breathable |
| 20 kPa | 4.67 × 10^17 kg | Very low-pressure habitat design space |
| 40 kPa | 9.34 × 10^17 kg | Intermediate plateau |
| 60 kPa | 1.40 × 10^18 kg | Lower edge of reference open-world range |
| 80 kPa | 1.87 × 10^18 kg | Reference case |
| 101.325 kPa | 2.37 × 10^18 kg | Earth sea-level pressure on the Moon |

Every additional kilopascal requires about **2.34 × 10^16 kg** of gas. Pressure optimization therefore has civilization-scale consequences.

## 4. Reference gas inventory

The reference dry mole fractions are 75% N2, 24% O2, 0.9% Ar, 800 ppm CO2, and 200 ppm other trace gases. The resulting mean molar mass is about 29.09 g/mol.

| Gas | Partial pressure | Approximate mass | Function |
|---|---:|---:|---|
| Nitrogen | 60.0 kPa | 1.35 × 10^18 kg | Inert buffer, nitrogen cycle and pressure reserve |
| Oxygen | 19.2 kPa | 4.93 × 10^17 kg | Respiration, oxidation and ozone precursor |
| Argon | 0.72 kPa | 2.31 × 10^16 kg | Optional inert buffer |
| Carbon dioxide | 0.064 kPa | 2.26 × 10^15 kg | Photosynthesis and greenhouse control |
| Other trace gases | 0.016 kPa | about 3.7 × 10^14 kg | Water vapor, ozone, methane and controlled species |

Mole fraction and mass fraction are different. Oxygen is 24% by molecule count but about 26.4% by mass; nitrogen is 75% by molecule count and about 72.2% by mass.

The reference oxygen partial pressure is close to Earth sea-level oxygen availability, but this does not establish fire safety. Ignition, flame spread and material compatibility depend on oxygen fraction, oxygen partial pressure, total pressure, humidity, flow and low-gravity combustion behavior.

## 5. Column mass and shielding

The column mass above one square meter is:

`m_column = P / g ≈ 49,261 kg/m² at 80 kPa`

This is about 4.8 times the mass column of Earth’s sea-level atmosphere because the same pressure requires more overlying mass in lower gravity. That large column is potentially valuable for cosmic-ray, solar-particle, ultraviolet and meteoroid protection, although shielding quality must be computed with particle-transport models because secondary radiation can matter [S10].

## 6. Scale height and vertical depth

For an isothermal atmosphere:

`H = R T / (M g)`

At 285 K and mean molar mass 29.09 g/mol, `H ≈ 50.2 km`. The reference dry adiabatic lapse rate is approximately:

`Γ_d = g / c_p ≈ 1.62 K/km`

Both values are radically different from Earth’s because lunar gravity is weak. A convective atmosphere can therefore be tens of kilometers deep without exhausting its temperature margin.

### Isothermal pressure screen

| Altitude | Screening pressure from 80 kPa surface |
|---:|---:|
| 10 km | 65.6 kPa |
| 20 km | 53.7 kPa |
| 35 km | 39.8 kPa |
| 40 km | 36.1 kPa |
| 45 km | 32.6 kPa |
| 100 km | 10.9 kPa |
| 200 km | 1.49 kPa |
| 300 km | 0.202 kPa |
| 500 km | 0.00375 kPa, or 3.75 Pa |

These numbers deliberately over-simplify the real atmosphere. Temperature, changing gravity, condensation, photochemistry, molecular diffusion and a hot thermosphere all alter the profile. The calculation is still enough to show that present-day 50–100 km low lunar orbits would lie deep inside a thick atmosphere.

## 7. Twilight geometry

The Sun moves across the lunar sky at about `360° / 708.734 h = 0.50795°/h`. A five-hour practical transition corresponds to about **2.54 degrees** of solar elevation change.

Terluna therefore defines the canonical five hours as the operational brightening or dimming interval around local sunrise or sunset. Very faint astronomical scattering can begin earlier, and mountains can delay or advance local illumination. The five-hour transition is part of the 708.7-hour cycle; it is not extra time added to the bright and dark halves.

A schematic equatorial cycle is:

- Practical dawn: 5 h.
- Full daylight: about 349.37 h.
- Practical dusk: 5 h.
- Full night: about 349.37 h.

## 8. Hydrosphere scales

Global equivalent depth is:

`GED = M_water / (ρ_water A)`

| Water mass | Global equivalent depth | Mean depth over 10% of surface |
|---:|---:|---:|
| 1 × 10^17 kg | 2.64 m | 26.4 m |
| 3 × 10^17 kg | 7.91 m | 79.1 m |
| 1 × 10^18 kg | 26.36 m | 263.6 m |

The 3 × 10^17 kg reference inventory is equivalent to a pure-water sphere about 83 km in diameter. That comparison illustrates the transport problem without implying delivery as one impactor.

If all reference water were manufactured using lunar oxygen, the imported hydrogen alone would be about **3.36 × 10^16 kg**. Importing water or hydrated material may be operationally easier, while indigenous polar ice should first be treated as a limited scientific and strategic resource [S03][S07][S08].

## 9. Oxygen production from regolith

Lunar regolith is roughly 40–45% oxygen by mass bound in oxides, and oxygen-extraction processes have been demonstrated at laboratory scale [S04]. If feedstock contains 45% oxygen and a planetary process recovers 75% of it, producing the reference 4.93 × 10^17 kg atmospheric oxygen requires processing:

`M_regolith = M_O2 / (0.45 × 0.75) ≈ 1.46 × 10^18 kg`

At a bulk density of 1,600 kg/m³, that equals a global average layer near 24 m, although real mining would be concentrated into selected industrial provinces. The operation yields metals or metal alloys as coproducts, so the atmosphere program and construction economy should be designed together.

## 10. Oxygen-production energy

Using an illustrative system energy intensity of 20–50 MJ per kilogram of product oxygen gives:

- Total process energy: about **9.9 × 10^24 to 2.5 × 10^25 J**.
- Average power over 1,000 years: about **313 to 782 TW**.

This excludes much of the mine development, comminution, transport, volatile import, nitrogen processing, water circulation, habitat construction and loss replacement. It is a scale marker, not a cost estimate.

## 11. Planetary thermal-inertia check

At a screening Bond albedo of 0.30, Terluna would absorb roughly **9 petawatts** averaged over the illuminated disk and time. Over one 14.765-day half-cycle, the absorbed energy is of order **1.1 × 10^22 J**.

The heat capacity of 3 × 10^17 kg of liquid water is about 1.25 × 10^21 J/K, so a 10 K participating temperature swing stores about **1.25 × 10^22 J**. The 80 kPa atmosphere has a comparable 10 K sensible-heat scale of about **1.9 × 10^22 J**.

This is one reason the long night is not automatically fatal: the reference atmosphere and hydrosphere contain enough thermal capacity in principle to buffer energy on the same order as a half-cycle’s absorbed sunlight. Distribution, clouds, ice, infrared cooling and heat-transfer rates determine whether that capacity is actually usable.

## 12. Flight scaling

For the same aircraft, atmospheric density, wing area and lift coefficient, stall speed scales with the square root of weight. Reducing gravitational acceleration to 0.16 g gives an idealized stall-speed factor of:

`sqrt(0.16) ≈ 0.40`

A vehicle that stalls at 25 m/s on Earth could, ignoring redesign and atmospheric differences, stall near 10 m/s on Terluna. Aircraft can trade this benefit for larger payload, smaller wings, slower flight or higher operating altitude.

Buoyant lift in kilograms does not gain the same one-sixth factor because both buoyant force and payload weight scale with gravity. Airships still benefit from lower structural loads and a deep atmosphere.

## 13. Artificial-gravity geometry

For a rotating habitat, `a = ω²r`. Representative 1 g rotations are:

| Radius | Rotation rate for 1 g |
|---:|---:|
| 50 m | 4.23 rpm |
| 100 m | 2.99 rpm |
| 250 m | 1.89 rpm |
| 500 m | 1.34 rpm |

Large radii reduce rotation rate and head-to-foot gravity gradients. Terluna settlements should reserve structural corridors for rotating sleep quarters, clinics, maternity facilities and exercise habitats until partial-gravity thresholds are known [S14][S15].

## 14. Imported-material momentum and impact energy

A volatile program must never treat a giant comet impact as an ordinary delivery method. A 10^15 kg packet arriving at 2 km/s carries `2 × 10^21 J`; a 10^18 kg body at the same speed carries `2 × 10^24 J`.

Resources must be subdivided, processed and decelerated through high-orbit depots, electromagnetic systems, tethers, propulsion or controlled aerocapture after a suitable atmosphere exists. Momentum disposal is part of the mass budget.

## 15. Mass fraction of the Moon

The 80 kPa atmosphere is only about **2.5 × 10^-5** of lunar mass. It does not meaningfully alter the Moon’s orbit or surface gravity. The engineering difficulty comes from acquiring, processing and retaining the material, not from its effect on lunar bulk dynamics.

## 16. Reproducibility

The script `models/terluna_screening_calculations.py` reproduces core formulas. The CSV files in `data/` preserve pressure, altitude, water, mass and energy sensitivities. All outputs should be replaced by higher-fidelity models as the concept develops.
