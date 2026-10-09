# Biological precedents for holding gas

**Living things already make gas, hold it at great pressure and wall it in with crystals: fish fill their swim
bladders with oxygen at over a hundred atmospheres, guanine plates make the bladder wall a hundred times tighter, and
the Portuguese man-of-war secretes its own float gas.** Each precedent below was checked against primary sources
(the [register](gas_biology_sources.json) says how far each was read), and each carries into a floater number in
[gas_biology.py](gas_biology.py).

## The swim bladder's gas gland and rete mirabile

**Deep-sea fish secrete oxygen into their swim bladders against pressures above 100 atmospheres.** Fish with gas-
filled swim bladders live at least as deep as 266 atmospheres; the bladder gas is mostly oxygen, with nitrogen at
2–15% of the oxygen, in fish caught at 320–1,300 m (Scholander & van Dam 1954). The rete mirabile's capillaries are
7–13 mm long in deep-sea species and lengthen generally with depth (Scholander 1954).

**The gas gland acidifies the blood, and the rete multiplies the effect by countercurrent exchange.** Gas-gland
cells turn about 80% of the glucose they take up into lactate and make CO₂ in the pentose phosphate shunt; the acid
switches on the Root effect, which lowers the blood's oxygen capacity, and lactate salts gases out. Back-diffusion
from the venous to the arterial capillaries of the rete multiplies this single effect (Kuhn et al. 1963, as reviewed
by Pelster 2015). In the eel's rete about 58,000 arterial and 44,000 venous capillaries, 1–2 µm apart over a few
millimetres, raise the oxygen and CO₂ partial pressures seven- to eightfold (Pelster 2023). The same countercurrent
keeps secreted gas in the bladder while blood keeps circulating (Scholander & van Dam 1954, summarising Jacobs 1930).

**For a floater the rete is a refinement.** If the whole reference transpiration stream, 1.25 kg/m² a day, wetted the
gas cells and left saturated with hydrogen at 98 kPa, it would carry off 3.9% of the hydrogen that permeation loses
(Henry's-law solubility 7.7–7.8 × 10⁻⁶ mol m⁻³ Pa⁻¹ at 25 °C; Sander's compilation). A countercurrent exchanger
returning 90% of it cuts that to 0.4%, and the thrifty water needs cut it further, to 0.2–0.8% before any exchanger
(*screen*).

## Guanine plates in the bladder wall

**Removing the silvery guanine layer makes a swim-bladder wall about a hundred times more permeable.** The intact wall
is far less permeable to carbon dioxide, oxygen and nitrogen than ordinary connective tissue, and overlapping
crystals account for it; two deep-sea eels carried about ten times the guanine per unit area of the conger eel
(Denton, Liddicoat & Taylor 1972). The low oxygen conductance lies in the wall's middle layer and comes from a low
diffusion coefficient caused by multiple layers of very thin (about 0.02 µm) and broad (up to 100 µm) crystals that
occupy only a small part of the tissue (Lapennas & Schmidt-Nielsen 1977). Sardine swim-bladder plates are 2–50 µm
wide and under 20 nm thick (Pinsk et al. 2022): aspect ratios of 100 to 5,000 (*derived*). At 3,000 m the oxygen
pressure difference across a swim-bladder wall is close to 300 atmospheres (Wittenberg et al. 1980).

**Aligned plates of that shape could cut a floater barrier's hydrogen permeation 10–50 times.** For impermeable plates
at volume share φ and aspect ratio α, Nielsen's tortuosity model gives P/P₀ = (1 − φ)/(1 + (α/2)φ) (Nielsen 1967, in
the form given by Idris et al. 2022). Plates at 2% of the barrier's volume with aspect ratios of 1,000–5,000 cut
permeability to 1.9–8.9% of the plain film's, and 4% at 5,000 gives the hundredfold cut that stripping the guanine
layer reverses. A plated barrier a tenth to a fiftieth as thick holds hydrogen as well as the reference's 100 µm,
since the thickness for equal loss scales with P/P₀. It saves most of the barrier's 0.15 kg/m², lift that the body
can spend on tissue and water ([structure](structure.md)).

## The Portuguese man-of-war's gas gland

**The man-of-war secretes carbon monoxide into its float from serine.** Its float gas holds 0.5–13% carbon monoxide
and 15–20% oxygen with negligible CO₂; the gas gland (pneumadena) forms the carbon monoxide from L-serine, and the
gland carries a strikingly large concentration of folic acid. Secretion inflates the float, and the carbon monoxide
is later replaced by air through diffusion and exchange (Wittenberg 1960; Munro et al. 2019 describe the gland's
aeriform cells). It is a colonial animal that makes its own lifting gas in a dedicated tissue and sails with it, the
closest living model for a floater's gas organ.

## Siphonophores: colonies of specialised parts

**Siphonophores build one body from many specialised zooids budded from growth zones.** Each zooid is homologous to
a solitary animal; a colony carries a gas-filled float (the pneumatophore), swimming bells (nectophores), feeding
polyps (gastrozooids), palpons, reproductive gonophores on gonodendra, and bracts, and new zooids arise as probuds in
the growth zones that also lengthen the stem (Dunn & Wagner 2006). The man-of-war is a siphonophore that uses its
float as a sail (Munro et al. 2019). A floater colony can follow the plan: float modules, collector modules with
fibre fringes, gas-making modules, sail and tether modules, and reproductive modules that bud juveniles, all added
at growth zones as the colony grows and regrown where storms or fire take them.

## Goldbeater's skin

**Processed animal membrane held airships' hydrogen.** Finished goldbeater's-skin fabric of processed intestinal
membranes, crossed in layers with cotton, glue and varnish, weighed 130–150 g/m² and passed a few litres of hydrogen
per square metre a day (Chollet 1922, in [envelope](envelope.md)). It shows that a collagen membrane makes a
practical hydrogen barrier, and it is the precedent for harvesting floater membrane for sky towns' gas cells
([resources](resources.md)).

## Sources and checks

[gas_biology_sources.json](gas_biology_sources.json) records each source with what was read. The tests
([test_gas_biology_resources.py](test_gas_biology_resources.py)) check the plate model's limits, its hundredfold case,
the equal-loss thickness and the dissolved-loss scaling.
