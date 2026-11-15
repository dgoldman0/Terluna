# Reclamation in an enclosed industrial economy

The [reclamation experiment](reclamation.py) follows water and materials through use, recovery, rejection and replacement, then recalculates the equipment, electricity, freight and mining needed to sustain those flows. It retains the enclosed population and commercial orderbook, the chemical freight routes and the source grades of the [published physical baseline](README.md). Its 36 cases cover the services, propellant/logistics and diversified portfolios under strong markets, with either 0.5% or 1% annual real market growth. No Open Moon, planetary array or outer-system buildout enters the calculation.

The recovery and refurbishment coefficients are engineering scenarios. The numerical results show what those specified capabilities would change; they do not establish demonstrated industrial yields, qualified repair procedures or accessible source reserves. The [source record](reclamation_sources.json) separates the primary evidence from the assumptions in [the experiment inputs](reclamation_inputs.json).

## Closure already present

At the inherited normalization of 250,000 residents, the water model processes 9.131 billion kg/year of domestic water, 1.826 billion kg/year of industrial water, 29.311 billion kg/year of crop water and 0.486 billion kg/year of fish-system purge. Recovery is already credited at 98%, 98%, 99.5% and 95%, respectively. Their combined losses are 389.992 million kg/year. A separate landscape/agriculture allowance adds 456.563 million kg/year, giving the inherited 846.554 million kg/year of operational water makeup. Visitors receive the same per-guest-year requirements. Initial settlement, visitor and construction stocks are additional dated requirements.

The separate allowance originates in [the maintenance study](../../research/studies/vacuum_crater_buildout/maintenance.md): 5 kg per resident-day for garden, lake and food-production losses beyond its domestic/process loop. The later settlement calculation retains this allowance while also adding explicit crop and fish losses. Its allocation is therefore unresolved, with possible overlap. The main experiment retains all 5 kg; a separate endpoint removes that allowance alone as a maximum overlap or complete-capture credit. That endpoint preserves every explicit water-loop loss and does not establish that the allowance is entirely redundant.

The baseline also credits 90% material return for ordinary household goods and 95% usable return from retired local structures. Structural bulk and metal have carried inventories; recovered material first supplies subsequent fabrication, and only the remaining demand requires virgin feed. Full inherited material-processing energy remains charged. Active machinery components, including imported tug equipment, receive no baseline refurbishment credit. Recovering structural material thus does not eliminate the much larger cost of purchasing replacement active systems.

## Prospective comparisons

| Variant | Changes from the inherited baseline |
|---|---|
| `baseline` | Exact control, including all existing recovery and the landscape/agriculture water allowance. |
| `closure` | Domestic and industrial water recovery rise to 99.5%, crop recovery to 99.7%, and fish-purge recovery to 98%. Structural collection and quality acceptance are each 99%, giving 98.01% usable return. Goods collection of 98% and quality acceptance of 97% give 95.06% usable return. Extra sorting and water-treatment equipment and energy are charged. |
| `refurbishment` | Selected mechanical active components receive 80% collection and 60% acceptance after inspection. New imported parts replace 20% of each accepted unit's mass; the displaced old parts receive no further credit. Water and structural recovery retain baseline rates. |
| `combined` | Both closure and refurbishment changes, including their costs and physical feedback. |
| `tug_life_30` | The source tug's active mean life rises from 20 to 30 years. Other component lives and recovery remain unchanged. This durability target is separate from refurbishment. |
| `combined_without_residual_allowance` | The combined experiment with only the legacy 5 kg/person/day water allowance removed. This is a credit ceiling for overlap or complete capture, not a replacement baseline. |

The selected refurbishment assets are source tugs, lunar landers, source mining/metal/fabrication trains, and lunar excavators, haulers and bulk/metal/fabrication trains. Repair occurs at each asset's inherited destination: source facilities service the source fleet and plants, depot facilities service lunar landers, and lunar facilities service surface equipment. This is a co-located workshop design assumption. No retired source hardware is moved to a lunar workshop without freight.

A component retiring in year *t* enters inspection no earlier than year *t+1*. Accepted stock can then supply a two-year work order, returning to service no earlier than year *t+3*. Components retain their original expected service life after successful refurbishment; that reset, including repeated later refurbishment, requires engineering qualification. Accepted stored units do not age while awaiting an order in this scenario. Parts that are not collected, fail inspection or are displaced by new repair parts remain uncredited. The model tracks the pending inspection stock, accepted inventory and work in progress separately.

For an order requiring 100 kg, a usable old unit contributes 80 kg of retained parts and requires 20 kg of newly imported parts. The procurement saving is the retained mass, valued at the inherited component price per kilogram. The 20 kg of discarded old material is recorded as waste. Imported repair parts are assumed to have that same unit purchase price; a premium for specialized spares is not calibrated. Inspection and remanufacturing consume 1 and 10 kWh per kilogram processed, respectively, and import consumables equal to 0.5% of their combined processed mass. These are explicit design coefficients. Direct labour handles 100,000 kg per worker-year of inspection or 20,000 kg per worker-year of remanufacturing; it is added to the inherited equipment-support staffing proxy and paid once.

Workshops use a stated 50,000 kg equipment module processing 500 kg/hour, with the inherited 80% process availability and 25% capacity margin. Their active and structural mean lives are 15 and 30 years. Additional water treatment uses 1,000 kg modules at 1,000 kg/hour, with 10-year active and 30-year structural lives, and 0.02 kWh per kilogram of gross processed water. Additional structural sorting uses 0.02 kWh/kg bulk or 1 kWh/kg metal; its equipment is sized from inherited bulk/fabrication-train intensities. Added recovery of household goods also receives sorting throughput and 1 kWh/kg of extra recovery work. These equipment and productivity assumptions are prospective, rather than estimates fitted to lunar operating data.

The fixed-point solver includes those facilities' own replacements, imported consumables, power generation, heat rejection, freight and fuel. Full bulk/metal conversion energy remains in the inherited account, so recovered material never becomes free processed feed. Additional direct repair workers enter the workforce screen. The source chemistry/utility equipment intensity also applies to additional source process power as a conservative inherited proxy; it is not a detailed repair-shop engineering design.

## Physical results

All 36 calculations converge and pass the dated equipment and resident-workforce screens. The six baseline controls reproduce the frozen physical accounts exactly. The following services comparison retains 0.5% real annual market growth; operating figures include each variant's applicable process hardware, repair consumables and direct repair payroll.

| Services variant | Source feed through year 200, trillion kg | Year 200 source feed, trillion kg/year | Year 200 physical operation, trillion USD2041/year |
|---|---:|---:|---:|
| `baseline` | 164.593 | 3.008 | 4.128 |
| `closure` | 136.481 | 2.386 | 3.608 |
| `refurbishment` | 161.734 | 2.945 | 3.409 |
| `combined` | 134.073 | 2.335 | 3.029 |
| `tug_life_30` | 162.168 | 2.965 | 3.592 |
| `combined_without_residual_allowance` | 81.222 | 1.170 | 2.318 |

Across the six orderbooks, the combined scenario reduces cumulative source feed through year 200 by 12.10–18.54% and physical operating outlays in year 200 by 9.99–26.64%. The services and diversified portfolios differ substantially because their ongoing production, exports and replacement bills differ. These are physical cash comparisons before commercial receipts, financing or the humanities capital-cost allowance.

| Portfolio and real market growth | Baseline cumulative source feed through year 200, trillion kg | Combined cumulative source feed, trillion kg | Baseline year 200 operation, trillion USD2041/year | Combined year 200 operation, trillion USD2041/year |
|---|---:|---:|---:|---:|
| propellant/logistics (0.5%) | 176.774 | 146.038 | 4.236 | 3.107 |
| services (0.5%) | 164.593 | 134.073 | 4.128 | 3.029 |
| diversified (0.5%) | 231.014 | 199.388 | 9.735 | 8.356 |
| propellant/logistics (1%) | 187.781 | 156.862 | 4.410 | 3.235 |
| services (1%) | 164.721 | 134.194 | 4.135 | 3.034 |
| diversified (1%) | 264.624 | 232.616 | 16.812 | 15.132 |

For these six orderbooks, the baseline source-feed range is 3.956–6.102 trillion kg through year 100, 164.593–264.624 trillion kg through year 200, and 225.267–385.078 trillion kg through year 220, including seed demand once. The combined experiment changes those ranges to 3.380–5.497, 134.073–232.616 and 181.276–338.394 trillion kg. The zero-legacy-allowance endpoint gives 2.394–4.512, 81.222–179.780 and 105.126–262.265 trillion kg. These are selected-orderbook ranges rather than statistical confidence intervals.

The water-source account distinguishes settlement makeup, depot demand and source-ferry fuel. Chemical exhaust has zero recovery credit in every case. Reducing imported replacement mass also reduces its transport fuel and the mines, power equipment and vehicles needed to produce and deliver that fuel. This feedback is why a hardware-repair intervention can reduce virgin source-feed demand even when it recycles no engine exhaust. The three-year source-route convention, annual-rate cargo arrival and unresolved mission scheduling remain inherited conditions.

The calculation reports through year 220. It retains the original solver's 40-year constant-demand continuation to size near-boundary requirements, but does not publish a qualified year-300 industrial path. A resource consumer may separately show how long a specified stock lasts at the reported final annual demand; that is an analytical depletion comparison, not an extension of the physical buildout.

## Evidence and interpretation

NASA's [water-recovery report](https://www.nasa.gov/missions/station/iss-research/nasa-achieves-water-recovery-milestone-on-international-space-station/) describes an ISS demonstration reaching 98% total recovery using wastewater, cabin moisture, urine processing and brine treatment. It supports the recovery mechanism and the need for filtration, quality checks and maintained equipment. It does not validate the proposed higher recovery for lunar farms, industry or cities.

NASA's [manufacturing portfolio plan](https://ntrs.nasa.gov/api/citations/20250004020/downloads/ISM%20Portfolio%20Plan%20Baseline%20for%20Public%20Release%20R1c.pdf), sections 1.1.3–1.1.4, distinguishes promising recycling processes from qualified hardware production. It records the Refabricator's incomplete orbital demonstration and a terrestrial metal-chip-to-wire demonstration. Those findings motivate collection and quality losses in the experiment; they supply no aggregate lunar-fleet acceptance fraction. The [CLEAR research abstract](https://ntrs.nasa.gov/citations/20100033105) explains why component-level repair might reduce the logistics burden of replacing complete units, without establishing the experiment's repair success or service-life reset.

The inherited air account represents atmospheric leakage and new air stocks. Food, carbon and metabolic oxygen are not resolved sufficiently to add a second oxygen-recovery credit. NASA's [SCOR description](https://www.nasa.gov/spacecraft-oxygen-recovery-scor/) concerns recovery from exhaled carbon dioxide and retention of hydrogen otherwise lost in methane, a different flow from leaked atmosphere or rocket exhaust. Nutrients, manufacturing scrap, active-hardware downcycling and element-specific alloy quality remain further model extensions. The reported balances close the represented loops; they do not claim a completely closed civilization.

## Product and reproduction

The [compressed product](results/reclamation.json.gz) has schema `terluna.provisioning.reclamation/1`. Each case has a top-level `variant` and `baseline_case`, a `physical_case` with the existing annual cash, equipment and workforce fields, and a distinct `annual` loop audit for years 0–220. The nested physical case retains the original orderbook name and legacy scenario metadata for compatibility; the top-level variant identifies the experimental physics. Loop masses are subsets or circulations within the physical account; adding recovered mass to virgin demand would double count them. All money remains constant USD2041, with no new U conversion or financial-return appraisal.

Seed feed appears once across years 0–7. `seed_audit.physical_annual_average` and `late_preproduction_physical_total` retain the original average deployment and the separately charged years 6/7 procurement. Their sum supplies exact seed water, nitrogen and virgin-metal requirements under that inherited convention. Detailed dated seed circulation and retirement balances remain null where the producer has no executable deployment sequence. The preproduction average does not establish a mine commissioning chronology.

```sh
PYTHONDONTWRITEBYTECODE=1 python3 -m provisioning.enclosed_industry.reclamation
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest provisioning/enclosed_industry/tests/test_reclamation.py -q
```

The runner uses one numerical worker. Cached cases require identical producer, complete case and input identities. The six baseline controls must reproduce the frozen product exactly; a mismatch aborts production. Tests also check delayed availability, imported repair-parts substitution, rejected mass, source-location preservation, finite stocks, seed inclusion, cash totals and reproduction of a complete combined case. These numerical checks establish the implementation's accounting behavior separately from empirical validation.

[The verification record](reclamation-checks.md) gives the selected check scope, source-reading scope and final prose-only provenance refresh.
