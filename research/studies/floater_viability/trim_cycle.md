# Water, carbon and lifting-gas trim

An annual gas-leakage balance does not establish neutral flight. Rain, fog uptake,
evaporation, drainage, stored food and metabolic products change the mass that the
gas must support. In some otherwise favourable barrier cases, routine water
turnover requires much more control than intact-film leakage. This is a coupled
mass and pressure requirements calculation, not a predicted trajectory or an
established biological flight mechanism.

`trim_cycle.py` uses the existing [mechanics](mechanics.md),
[water and thermal model](environment.md), [gas barrier assessment](envelope.md)
and [organism carbon budget](biology.md). It adds no empirical trait measurements.
Atmospheric density, pressure and mixture lift come from the retained
[sky-ship product](../sky_ships/results/sky_ships.json). Chemical energy uses a
declared rounded H₂ lower heating value of 120 MJ/kg; fermentation uses the
separately documented measured reactor-yield sensitivity. These are alternative
accounting routes, not two costs charged to the same hydrogen.

## The actual mass ledger

At a fixed ambient state, define supported mass `M` as all actual solids, liquids
and payload, excluding the lifting mixture and ambient-air ballonet contents.
Define `L` as the gas's displaced-air mass minus its own mass. Then

\[
\Delta R=\Delta L-\Delta M,
\]

where positive residual `R` tends upward. Hydrogen weight is already included in
`L`; subtracting the escaping hydrogen mass a second time gives a wrong answer.
The ledger explicitly accepts humidity uptake, rain, intake, evaporation,
drainage, water dumping, dry-mass change, metabolic-water change, hydrogen input,
leakage and deliberate venting. Opposing simultaneous flows are netted before
control is selected. A tank's unfilled capacity does not count as actual water.

For a mixture containing H₂ mole fraction `x` and otherwise ambient air at the
same pressure and temperature,

\[
b=\frac{\rho_a-\rho_g}{\rho_g-(1-x)\rho_a}
\]

is net supported kg per kg of H₂. This constant-composition volume-change screen
assumes purity control and suitable expansion or contraction; it is not the
constant-volume gas-exchange model. At the retained 10 km state, `ρa=1.204`,
`ρg=0.106 kg/m³`, `x=0.98` and `b=13.4033 kg/kg`. Consequently, losing 1 kg of
water while holding altitude with gas venting removes **0.074608 kg H₂**. In a
fixed outer hull, the removed gas volume is replaced with real ambient air; gas
loss, air admission and water loss close the total-mass balance.

After a supported-mass change `ΔM` and natural hydrogen leak `q`, the additional
signed hydrogen correction is `ΔM/b+q`: positive requires production or stored
gas; negative requires additional venting under this particular control choice.
For example, with `b=10`, evaporation of 5 kg and leakage of 0.1 kg H₂ requires
0.4 kg additional venting, not 0.5 kg venting plus replacement of the leak. Later
recovering the 5 kg water while another 0.1 kg H₂ leaks requires 0.6 kg H₂.
Over that closed payload/gas cycle, production equals natural leakage plus
actual venting. The tests check this identity independently.

## Water inventory and the price of repeated venting

`water_inventory_band` represents actual free water with explicit lower and
upper bounds. With 5 kg/m² present, a 1 kg/m² minimum and 9 kg/m² tank capacity,
only 4 kg/m² is available for routine dumping or net evaporation. A 1 kg/m²/day
loss reaches the minimum in four days. Seven days requires an additional
3 kg/m² of supply; the function reports the deficit rather than silently
clamping the water budget. Tissue hydration and retained-community water cannot
be borrowed without an independently justified physiological allowance.

The following conditional cycles lose the stated water and then regain it once
per day at 10 km. They hold altitude by venting on the loss and making new gas on
the refill; natural leakage is omitted so the table isolates trim. All quantities
are per projected m². Continuous matched uptake and loss would not incur this
gross-cycle trim burden.

| Daily water swing | H₂ vented and later replaced | Mean H₂ chemical output required | Fermentative gross feed carbon, if repeated all year |
|---|---:|---:|---:|
| 0.01 kg/m² | 0.000746 kg/m² | 1.04 W/m² | 4.64 kg C/m²/year |
| 0.1 kg/m² | 0.00746 kg/m² | 10.4 W/m² | 46.4 kg C/m²/year |
| 1 kg/m² | 0.0746 kg/m² | 103.6 W/m² | 464 kg C/m²/year |
| 5 kg/m² | 0.373 kg/m² | 518 W/m² | 2,319 kg C/m²/year |

The fermentation column uses 2.1 mol H₂/mol glucose and does not credit residual
organic-acid recovery. It is gross diverted substrate carbon, not a claim that
all of it leaves as CO₂. Product recycling would require its own energy and mass
ledger. At an assumed cycle-mean irradiance of 200 W/m² and 1% full-solar-to-H₂
efficiency, the first row already needs 52% of the collecting area/time; the
1 kg row needs about 52 times the available area/time. Those efficiency and light
inputs are illustrative requirements, not measured floater performance.

For integration, pass chemical power from `hydrogen_refill_budget` into
`biology.photosynthetic_hydrogen_allocation`, or charge the reported annual
fermentation feed against available assimilate. Do not count that same hydrogen
again as a biomass-energy-equivalent charge. The remaining carbon must still fund
existing maintenance, turnover, consumers and reproduction; if starting from
NPP, baseline host respiration has already been paid. Under the explicit
no-coproduct-recovery fermentation screen, even a 0.01 kg/m² daily swing needs
more than the earlier 1 kg C/m²/year baseline NPP.

These are prescribed cycles, not forecasts. Their high cost motivates retaining
gas, adjusting altitude, smoothing water collection, changing body temperature
or avoiding large uncontrolled inventory swings. They do not prove all possible
floaters impossible.

## Retaining hydrogen still needs a pressure-bearing mechanism

Two distinct compression arrangements are provided. Both conserve the existing
gas molecules and require actual actuation and heat transfer.

**Contracting gas displacement with an ambient-pressure ballonet or shrinking
outer bladder.** At fixed gas mass, lowering net support by payload loss `δM`
requires `δV=δM/ρa`, because gas mass is constant. Using
`δM/(ρa−ρg)` would incorrectly assume that gas escaped. Isothermal pressure then
obeys `p1=p0 V0/(V0−δV)`. If a fixed hull admits ambient-pressure air around the
contracted bladder, the internal gas/air partition must bear the pressure
difference. An ordinary freely moving partition at equal pressure cannot do this.

For a neutral initial supported mass of 20 kg/m² at 10 km, initial gas volume is
18.215 m³/m². Losing 1 kg/m² requires contraction by 0.831 m³/m² and raises gas
pressure by **4.783 kPa**. If the earlier structure was sized for 400 Pa, retaining
the same allowance for the original loads requires approximately 5.183 kPa in
this additive screen. The unchanged 400 Pa structure does not pass this new
pressure requirement. Partition geometry, wrinkling, load redistribution and
the changed hydrostatic head still require a structural solution.

Reversible work done on the compressed gas is `p0 V0 ln(V0/V1)`. Ambient pressure
supplies `pambient (V0−V1)` of that work as volume contracts, so the minimum
additional work above ambient is their difference. The example gives 85.1 kJ/m²
on the gas and 1.96 kJ/m² above ambient. The latter is an ideal mechanical lower
bound, not demonstrated biological pump performance or a hardware mass budget.

**Fixed hull with gas and air compartments at a shared rising pressure.** Pumping
an actual air mass equal to the lost payload preserves total mass and external
volume while retaining all H₂. Here
`Δp=δM pambient/(ρa Vouter)`, and the gas volume decreases in inverse proportion
to the common pressure. Ballonet air becomes denser; assuming unchanged ambient
density and equal pressure simultaneously would be inconsistent. With outer
volume 30 m³/m² and the same 20 kg/m² supported mass, a 1 kg/m² loss raises
pressure by 2.772 kPa and has ideal isothermal air-pumping work 1.14 kJ/m².
This still exceeds the example's original 400 Pa budget. Enlarging unused hull
capacity can reduce the pressure rise, but the larger hull's mass and structure
must be supplied; capacity does not act as a free control force.

Neither arrangement proves that an organism has suitable valves, a compressor,
gas storage, a pressure-bearing partition or adequate heat rejection. Reversing
the water change requires expansion or air release with compatible geometry.
Stored high-pressure gas has its own vessel mass, heat and pressure costs.

## Night ballast and carbon mass need the same boundary

Dumping water of mass `b × night H₂ loss` offsets only that specified gas-loss
term. It does not establish night-long neutral flight if evaporation, dew,
rain or respiration also occur. The whole ledger must be reevaluated against
the actual free-water band. On later gas restoration, recovered water or another
real payload increase is needed to restore the prior neutral state.

Carbon accounting does not by itself supply payload change. For example, the
ideal starch reaction `C6H10O5 + 6 O2 → 6 CO2 + 5 H2O` consumes about 2.25 kg dry
reserve and forms 1.25 kg water per kg respired carbon. If all product water is
retained, payload falls by about 1 kg; if it escapes, payload falls by about
2.25 kg. The example ignores solutes and uses nominal atomic-mass ratios; the
ledger accepts the actual dry and metabolic-water changes separately. Food
transfer between organisms aboard changes neither whole-floater mass nor lift
until material enters or leaves its boundary.

Hydrogen production also has products and feedstocks. Water splitting consumes
water and releases oxygen; fermentation commits stored organic matter and forms
CO₂ and retained or expelled acids. Refill costs alone therefore do not close
physiology. A calculation that restores both water and gas inventories must
include the makeup water and product fates, or state that these are prescribed
external balances. Gas production from internal water cannot simultaneously
count that water as unchanged ballast.

Variable altitude, staged rain intake, synchronized gas production, finite gas
storage, aerodynamic forces and thermal buoyancy may avoid the vent/refill
strategy. The present outputs diagnose their required mass/pressure functions;
they do not solve a weather-resolved flight cycle. Passing a static lift or
annual carbon gate remains insufficient without that coupled trajectory.
