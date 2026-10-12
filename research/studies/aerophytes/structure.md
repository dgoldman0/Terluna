# Structure beyond the ideal arithmetic

**A single round aerophyte survives from about 50 m to 700 m across at 10 km with cellulose-class tendons, to 1 km
with 300 MPa fibre and to 3 km in canopy-like light; a colony of 60–150 m modules has no structural size limit.**
The ranges come from the coupled cases in [run.py](run.py) after joints, sustained load, storm gusts, the spare
superpressure for trim, wet fire-barrier skins, gas cells and tendon renewal are paid. The low edge is set by mass:
a sphere's gas column is two thirds of its diameter, and below about 45 m that column cannot lift the living tissue,
the water store and the wet skin with a 30% reserve. The high edge of a round body is set by carbon: its tendons'
mass per square metre grows with its size, and renewing them eats its photosynthetic surplus. A colony's modules
each pass the gates at their own size, and the ties that hold the colony together weigh the same per square metre at
any width.

| Body at 10 km (20 km) | Reference light, GPP 2 kg C/m²/yr | Carbon runs out at | Canopy-like light, GPP 5 kg C/m²/yr | Carbon runs out at |
|---|---:|---:|---:|---:|
| Round, collagen tendons (5.3 MPa sustained) | 100–200 m (none) | 324 m | 100–200 m (none) | 750 m |
| Round, cellulose-class fibre, 100 MPa (31 MPa at a century) | 50–700 m (60–700 m) | 856 m | 50–1,500 m (60–2,000 m) | 1.89 km |
| Round, 300 MPa fibre (92 MPa at a century) | 50–1,000 m (60–1,500 m) | 1.44 km | 50–3,000 m (60–3,000 m) | 3.22 km |
| Colony module, one-layer raft | 60–150 m modules (80–150 m) | 172 m | 60–300 m and beyond (80–300 m) | beyond 300 m |

All numbers are *screen* values from the product (`size_limits`) unless marked *derived*. The ranges are grid
diameters that pass; the carbon column finds where the surplus reaches zero by bisection on the model. Collagen bodies
above 200 m fail the lift reserve before carbon, because their tendons are heavy. The gates are: actual mass within
70% of the full lift, a positive carbon surplus after renewal, the daytime hydrogen share within the area-time budget,
and staying aloft after losing one gas cell by dropping free water. The joints, gas cells, wet skins and gusts below
enter every case.

**Every hull carries the spare superpressure its water strategy needs.** A round body banks storm rain within its
trim ([water](water.md)), so its hull carries 300 Pa above its base, head and gust; a colony module drinks as it
loses and carries 100 Pa (both design guesses, 100–1,000 Pa and 100–300 Pa). With the round bodies' 300 Pa, colony
modules would pass at 60–100 m only (`colony_module_spare_300_pa`).

## Joints and sustained load set the allowable stress

**Joints and a decade of sustained load cut a 20 MPa coupon to a 5.8 MPa assembled skin.** The chain is coupon
strength × joint efficiency × duration-of-load factor ÷ design factor: 20 MPa × 0.7 × 0.62 ÷ 1.5 = 5.79 MPa. The
reference's imposed 5 MPa therefore follows from the measured wet coupon strengths of 14–33 MPa (tomato cuticle and
nanocellulose–gelatin films; [mechanics](mechanics.md)).

**Wood carries 62% of its short-term strength for ten years and 58% for a century.** The Madison curve fitted to
126 long-time bending tests on clear Douglas fir is y = 108.4/x^0.04635 + 18.3, with x the seconds under load and y
the percentage of the standard five-minute strength (Wood 1951). It gives 67% at one year, 62% at ten, 58% at a
hundred and 56% at three hundred years, falling slowly toward 18.3%. Wood is the measured biological analogue for
a cellulose fibre composite, so the cellulose-class tendons and the film use this curve at their renewal interval.

**Collagen creeps to rupture at a tenth of its dynamic strength.** Wallaby tail tendon ruptures in creep at 10 MPa
and above, with time to rupture falling exponentially with stress between 20 and 80 MPa, against a dynamic yield
near 144 MPa (Wang & Ker 1995). Hind-limb tendons held at their muscles' maximum isometric force ruptured after
about 4.2 hours on average (Ker et al. 2000). A collagen tendon held below 10 MPa, through a joint at 0.8 and a 1.5
design factor, works at 5.3 MPa, which makes it no stronger in long service than the film.

The joint efficiencies (0.5–0.9 for film seams, 0.6–0.9 for tendon terminations), renewal intervals (1–30 years for
film, 30–300 years for tendons) and the 1.5 design factor are design guesses. The product carries their whole grid
(`structure.film_allowables`, `structure.tendon_allowables`).

## Tendons take a share of lift that grows with size

**The tendons of a round body take a share of its lift equal to (3π/4)ρp_b/(Δρσ) + (5π/4)ρgR/σ: the radius term
alone reaches 10% at 325 m radius for a 31 MPa fibre.** The upper hemisphere's pressure resultant across the equator
is πR²(p_b + ΔρgR) + (2/3)πR³Δρg, from the hydrostatic excess p(z) = p_b + Δρg(z − z_bottom). Meridional tendons of
length πR carry it at their allowable stress, so their mass is πRρF/σ and their share of the gross lift
(4/3)πR³Δρ is the sum above. The radius term grows in proportion to radius over the fibre's breaking length σ/(ρg):
12.8 km for the cellulose-class fibre at lunar gravity, 38 km for the 300 MPa fibre and 2.2 km for collagen. With
the round hulls' 445 Pa of base, spare and gust pressure, the first term is 4.7% of the lift for the 31 MPa fibre,
and the share reaches 10% at 174 m radius, 25% at 661 m and 50% at 1.47 km; for the 300 MPa fibre, 823 m, 2.29 km
and 4.72 km. For collagen the first term alone is 27%, so collagen tendons take more than a quarter of the lift at
any size, 39% in the coupled 100-m body (`round_tendon_ceilings`). The coupled cases use the existing lobed model,
which carries the crown pressure along the whole tendon and so sits above this exact resultant by up to 20%.

**End fittings, organ suspension and seams are small items.** Converging tendons put a hoop force F/2π into a crown
or base ring; with a ring a twentieth of the body's radius, the two rings weigh 3.2% of the tendons (2r_f/πR). Organs
and payload hung on fibre at 28 MPa cost 0.9% of their own mass per 100 m of load path at lunar gravity, because
the fibre's breaking length is 11.6 km. Seam overlaps of 2–5 cm on 1–4 m gores add 0.5–5% to the film. The coupled
round bodies carry 3.2% fittings and 2.5% seams.

## Gusts against a drifting body

**An aerophyte drifting at 10 km meets storm gusts of 14.4 m/s, a dynamic pressure of 125 Pa; at 20 km, 17.0 m/s and
148 Pa.** A drifting body rides the large-scale flow (the GCM's 50-year level of 3-day means is 8.7 m/s at 10 km), so
what it feels are the storms' sharper changes. The horizontal design gust is the jump from the CM1 ring's median
to its highest 3-hourly wind at the height (1.6 to 11.9 m/s at 10 km), and the vertical one the ring's strongest
updraft (9.8 m/s), which a body holding its height meets head on. Each takes the 1.4 factor the summit tower and
sky-ship studies apply for gusts sharper than the 6-km columns resolve.

| Height | Median, highest wind | Strongest updraft | Design gust | Design pressure | Centuries case (× 1.2) |
|---|---:|---:|---:|---:|---:|
| 1 km | 2.8, 15.6 m/s | 1.9 m/s | 18.0 m/s | 224 Pa | 21.6 m/s, 323 Pa |
| 5 km | 2.5, 13.0 m/s | 8.2 m/s | 14.8 m/s | 143 Pa | 17.8 m/s, 206 Pa |
| 10 km | 1.6, 11.9 m/s | 9.8 m/s | 14.4 m/s | 125 Pa | 17.3 m/s, 180 Pa |
| 20 km | 3.4, 15.6 m/s | 11.3 m/s | 17.0 m/s | 148 Pa | 20.4 m/s, 214 Pa |
| 40 km | 10.1, 23.0 m/s | 15.5 m/s | 21.7 m/s | 172 Pa | 26.0 m/s, 248 Pa |

The 1.2 centuries factor is the summit tower's decades factor; a centuries return level needs storm records longer
than one lunar day of one ring.

**Bodies larger than about 10 m meet a gust front's full jump.** A neutrally buoyant sphere with added mass 0.5
sheds half its relative wind in τ = 8(1 + C_a)R/(3C_d u₀): 0.6 s at 1 m radius, 12 s at 20 m, 59 s at 100 m and
590 s at 1 km for the 14.4 m/s jump. Streamlining (C_d 0.05) lengthens each tenfold.

**Gusts load small bodies hardest and the steady head dominates giants.** With the round hull's 320 Pa of base and
spare pressure counted as steady, the gust's share of the design membrane pressure is 27% at 5 m radius, 24% at 20 m,
16% at 100 m and 3% at 1 km. Small bodies therefore meet the largest stress cycles, and giants carry a steady load
that creep governs.

## Damage, tear arrest and repair

**An arrested 4-m tear empties one gas cell in minutes in a small body and in hours in a giant.** A slit of length
L under tension opens to about 0.05–0.2 L² (a design guess) and passes C_d A sqrt(2Δp/ρ) at the crown pressure (the
[mechanics](mechanics.md) orifice law). With the tear stopped at a 4-m ripstop or lobe boundary and the opening at
0.05 L², one cell empties in 2 minutes in a 50-m body (8 cells), 31 minutes in a 200-m body (10), 2.7 hours in a
500-m body (20) and 8 hours in a 1-km body (22–52). A colony module of 60–150 m empties in 0.6–4.8 hours. Larger
cells empty slowly through a fixed tear because their volume grows as R³ while the jet grows only as the square root
of the crown pressure. Healing faster than these times saves the cell; otherwise the body drops water for it.

**Losing one gas cell is survivable when the free water outweighs its lift.** With 4 kg/m² of free water above the
store's 1 kg/m² floor, a body needs at least supported mass ÷ 4 cells, and never fewer than eight; each body takes the
fewest that serve: for the cellulose-class body 8 up to 150 m, 10 at 200 m, 20 at 500 m and 52 at 1 km; for the 300
MPa body 8 up to 300 m, 12 at 500 m, 22 at 1 km, 64 at 2 km and 128 at 3 km. The cells are slack gas cells inside the
pressure-bearing hull, as in a rigid airship: a lost cell's volume fills with ballonet air at the hull's pressure, so
no wall carries the crown pressure and each wall stays a light partition (0.05 kg/m² plus its barrier). Walls between
cells hold hydrogen on both sides and lose none; only the outer skin and the gas–air partition pass hydrogen out.
Walls that had to hold the crown pressure across a lost neighbour would weigh tens of kilograms per square metre in a
500-m body and shrink every range. Ripstop fibres at 2–6 m spacing and the lobe boundaries arrest tears; the spacing
is a design guess.

**The skin holds the water a burning neighbour cell would dry.** A torn cell's hydrogen burns in the air outside
([storms](storms.md)); the skin of the next cell, two cell radii away, must stay wet through the fire's radiant heat.
Each coupled body carries the larger of 0.5 kg of water per square metre of skin and that need, as mass: for cellulose
and 300 MPa bodies 0.5 kg/m² up to 300 m, 0.5–0.65 at 500 m, 0.7–1.3 at 1 km and 1.4–2.8 at 2 km, which is 2.2–12 kg
per projected square metre; collagen bodies of 150–200 m carry 0.62–0.75 kg/m² and colony modules of 60–150 m 0.62–1.0
kg/m². The dose model radiates the burning cell's heat from its centre. Skin that adjoins the burning cell, one cell
radius away, would take four times the dose; carrying that water moves the smallest round body to 60 m and the
smallest colony module to 100 m (`round_fire_at_one_radius`, `colony_module_fire_at_one_radius`), and the fire
distance is a design range of one to two cell radii.

## Creep, fatigue and renewal over centuries

**A giant lives for centuries on young parts.** Living renewal sets the duration each part must carry its load: the
film is renewed every 10 years (62% of its short-term strength) and the tendons and colony ties every century (58%).
Renewal costs carbon in proportion to the renewed mass: the cellulose-class 500-m body renews 47 kg/m² of tendon
every century for 0.26 kg C/m² a year and keeps 0.58 kg C/m² of surplus; at 1 km its 157 kg/m² of tendon cost 0.87 kg
C/m² a year, the largest single item, and the surplus is gone. Hydrogen replacement costs a steady 0.2 kg C/m² a year
at every size, because only the outer skin and one partition lose gas.

**Storm gusts load a body 1,800–180,000 times a year, 0.5–54 million times in 300 years.** Fifty to five hundred
storm hours a year at gust periods of 10–100 seconds give the count (design guesses). Small, gust-dominated bodies
need fatigue-tolerant films and frequent renewal; a giant's stress is 97% steady, so creep-rupture at its renewal
interval sets its allowable.

## Colonies and flat bodies escape the head

**A module 36 m across minimises skin per gas volume, and a colony of such modules carries the same skin share at any
size.** A spherical module's skin per gas volume is 3ρ_s p₀/(2σ) + 3ρ_sΔρg r/σ + 3μ_b/r, with μ_b the barrier's areal
mass and the stress-bearing film at least its 25 µm floor; its minimum lies at r* = sqrt(μ_bσ/(ρ_sΔρg)). With the 5.8
MPa film and a 100 µm barrier (0.15 kg/m²), r* is 18.1 m and the skin takes 8.8% of the lift at the module's 120 Pa of
base and spare pressure (13.2% with the 125 Pa gust); the stress film there is above its floor. A 25 µm barrier moves
r* to 9.1 m. N modules of radius r carry the skin share of one module, whatever N is, because each module's
hydrostatic head is set by its own height. The coupled modules are larger, 60–150 m, because their payload of tissue,
water and wet skin needs a taller gas column; they carry the round bodies' 2.5% seam and ripstop allowance.

**The ties that hold a raft together weigh the same per square metre at any width: 0.37–0.92 kg/m² for 60–150 m
modules.** A one-layer raft of touching spheres on a hexagonal grid covers 1.103 times the modules' discs, so each
module lifts the payload of its hexagonal cell. The windward edge meets the gust pressure over the raft's thickness,
and the ties pass that force back through the raft: ρq h/σ per square metre, with h the module diameter, q the 125 Pa
design gust and σ the 31 MPa century allowable of tendon fibre. The ties renew once a century with the tendons. At the
centuries gust a 40-m layer with 28 MPa ties needs 0.39 kg/m² (`structure.colony_network`).

**A quilt of gas cells under one skin, 30–60 m thick, weighs 2.1–4.6 kg/m².** Its faces bulge between cell walls
20–40 m apart and the walls tie the faces together, each wall carrying p × w per metre of length. Like the raft,
its strength needs are set by its thickness and cell size, whatever its width. Wide, thin bodies are the giants'
structural form: lift comes from area times thickness, and the head from thickness alone.

## Sources and checks

The [source register](structure_sources.json) records Wood (1951), Wang & Ker (1995) and Ker et al. (2000), with
what was read. The winds are the [GCM wind product](../../../climate/results/gcm/global_winds_A28_dim5_moon.json) and
the [CM1 equatorial ring](../../../climate/results/crm/ring_ring_equator.json); the air is the
[sky-ship product](../sky_ships/results/sky_ships.json). [test_structure.py](test_structure.py) checks the Madison
curve, the joint chain, the equatorial resultant against a numerical hemisphere integral, the tendon-share law and
its inverse, the module optimum against neighbouring radii, the raft and quilt bookkeeping, the gust relaxation
against an integrated equation of motion, and the gusts against the ring's own numbers;
[test_run.py](test_run.py) checks that walls between gas cells lose no hydrogen, that skin water and joints enter the
mass once, that every hull carries its spare superpressure, and that the cells and skin water answer damage and fire.
An independent read-only review on 9 October 2026 found the gas-cell walls counted as hydrogen-loss area, the spare
superpressure and the fire-barrier water left out of the hull and mass, and the colony ties renewed too often; a
second found the cell count settling above the fewest, the colony modules without the round bodies' seams, the
carbon stops interpolated low and these tables without the spare pressure. The coupled cases carry all eight
corrections.
