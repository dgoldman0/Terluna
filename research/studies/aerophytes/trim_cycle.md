# Water, carbon and lifting-gas trim

**Venting gas to answer a daily water swing is ruinous, and keeping the gas means carrying the swing as
superpressure.** Losing 1 kg/m² of water and holding height by venting throws away 0.0746 kg of hydrogen per square
metre; making it back each day takes 104 W/m² of hydrogen output, against the reference's 2 W/m² of production
capacity. Keeping the gas instead raises its pressure by 4.78 kPa per kg/m² of swing at a 20 kg/m² load, or 2.77 kPa in
a fixed hull with an air ballonet. Intake must therefore track loss closely, and the [water](water.md) note gives the
swing each body's spare pressure absorbs. [trim_cycle.py](trim_cycle.py) uses the
[mechanics](mechanics.md), [environment](environment.md), [envelope](envelope.md) and [biology](biology.md)
components and the [sky-ship air](../sky_ships/results/sky_ships.json); chemical energy uses a rounded 120 MJ/kg for
hydrogen.

## The mass ledger

**The upward residual changes by ΔR = ΔL − ΔM**, where M is all carried solids and liquids and L the gas's displaced
air less its own mass. The gas's weight is already inside L. The ledger takes humidity uptake, rain, intake,
evaporation, drainage, dumping, dry-mass and metabolic-water changes, hydrogen input, permeation and venting, and nets
opposing flows before choosing a control. It counts the water held, whatever a tank's capacity.

**One kilogram of hydrogen in the lifting mixture supports b = (ρ_a − ρ_g)/(ρ_g − (1 − x)ρ_a) = 13.4033 kg at 10 km**
(ρ_a 1.204, ρ_g 0.106 kg/m³, x 0.98). Losing 1 kg of water and venting to hold height removes 0.074608 kg of hydrogen;
in a fixed hull, ambient air replaces the vented volume and the masses close.

**After a mass change ΔM and a natural gas loss q, the hydrogen correction is ΔM/b + q.** With b = 10, evaporating 5 kg
while losing 0.1 kg of hydrogen to permeation calls for 0.4 kg of extra venting; regaining the 5 kg while another 0.1
kg permeates calls for 0.6 kg of new hydrogen. Over the closed cycle, production equals permeation plus venting, which
the tests check.

## Water stores and repeated venting

**A store of 5 kg/m² with a 1 kg/m² floor and 9 kg/m² capacity holds 4 kg/m² for routine dumping;** a steady 1 kg/m² a
day of net loss reaches the floor in four days, and a week needs 3 kg/m² more. Tissue water stays out of the store.

Losing the stated water and regaining it once a day at 10 km, venting on the loss and making gas on the refill
(permeation left out):

| Daily swing | Hydrogen vented and replaced | Mean hydrogen output | Fermented feed carbon a year |
|---|---:|---:|---:|
| 0.01 kg/m² | 0.000746 kg/m² | 1.04 W/m² | 4.64 kg C/m² |
| 0.1 kg/m² | 0.00746 kg/m² | 10.4 W/m² | 46.4 kg C/m² |
| 1 kg/m² | 0.0746 kg/m² | 103.6 W/m² | 464 kg C/m² |
| 5 kg/m² | 0.373 kg/m² | 518 W/m² | 2,319 kg C/m² |

The fermentation column uses 2.1 mol H₂ per mol glucose, as gross diverted carbon. At 200 W/m² and 1% conversion the
first row already takes 52% of the collecting area-time; the 1 kg/m² row takes 52 times what there is. Hydrogen made
for trim is charged once, through `biology.photosynthetic_hydrogen_allocation` or as fermentation feed. Even a 0.01
kg/m² daily swing met this way costs more than the aerial-ecology screen's 1 kg C/m² of net production.

## Keeping the gas

**Shrinking the gas to shed lift at constant gas mass needs δV = δM/ρ_a and raises its pressure.** From a neutral
20 kg/m² at 10 km the gas volume is 18.215 m³/m²; losing 1 kg/m² takes 0.831 m³/m² of contraction and 4.783 kPa more
pressure, so a structure sized for 400 Pa would need about 5.18 kPa. The reversible work on the gas is 85.1 kJ/m², of
which 1.96 kJ/m² is above what the ambient air supplies. A freely moving equal-pressure partition cannot hold this
state; a pressure-bearing one, or a shrinking outer bladder, can.

**A fixed hull can instead pump in air equal to the lost payload, at a shared rising pressure
Δp = δM p_ambient/(ρ_a V_outer).** With 30 m³/m² of hull and the same 20 kg/m², a 1 kg/m² loss raises the pressure by
2.772 kPa at 1.14 kJ/m² of isothermal pumping work. More unused hull lowers the rise in proportion, and the hull's own
mass must be carried; this is the giants' advantage in [water](water.md), where V/A reaches 1.3 km.

## Night ballast and metabolic mass

**Dumping b × the night's hydrogen loss in water offsets that loss alone.** Evaporation, dew, rain and respiration
change the balance too, and the water dumped must later be regained. Respiring starch, C₆H₁₀O₅ + 6 O₂ → 6 CO₂ +
5 H₂O, uses about 2.25 kg of dry reserve per kg of carbon and makes 1.25 kg of water: payload falls by about 1 kg if the
water stays and 2.25 kg if it leaves. Food passed between organisms aboard changes neither mass nor lift. Splitting
water for hydrogen uses water and releases oxygen; fermenting it commits organic matter and makes CO₂ and acids. A
calculation that restores both water and gas includes the makeup water and those products.

Variable altitude, timed rain intake, matched gas production, gas storage, aerodynamic forces and thermal buoyancy are
the routes that avoid venting; [water](water.md) and [storms](storms.md) cost the ones that matter most.
