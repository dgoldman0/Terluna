# Floaters as a resource and a hazard to people

**Floaters offer people membrane, fibre, shade and habitat before food, and their hazards are their tethers, their
falls and the storms that carry them.** A 2-km colony renewing its skin every 10–30 years sheds 38,000–570,000 m² of
gas film a year at 10–50% harvest, enough for half to nearly eight sets of gas cells for the smallest sky ships. As
food, a floater is a garnish: harvesting 5–20% of its renewed living tissue feeds 380–1,500 people a year from a 2-km
colony. Numbers are *screen* values (`resources` in the product); harvest shares and renewal intervals are design
guesses.

## Membrane for gas cells

**A floater's shed skin is gas-cell film, as goldbeater's skin once lined airship cells.** The sky-ship study's hulls
of 245 m to 3 km carry 73,000 m² to 11 million m² of gas-cell film in sixteen cells. A 200-m round floater, renewing
its 138,000 m² envelope every 10–30 years and yielding 10–50% of it, gives 460–6,900 m² a year (a tenth of the
smallest ship's cells at most); a 2-km colony gives 38,000–570,000 m², 0.5–7.8 of the smallest ship's cells or a
twentieth of the largest's. Goldbeater's skin fabric weighed 130–150 g/m² (Chollet 1922); a floater's barrier film of
0.15–0.19 kg/m² is of the same order. Processing a living film into fabric, and its permeability after harvest, are
open.

## Food, fibre and habitat

**Food.** Harvesting 5–20% of the living tissue a floater renews each year (0.55 kg dry/m²) gives 0.15–0.6
person-years of food from a 40-m module, 3.8–15 from a 200-m body and 380–1,500 from a 2-km colony, at the aerial
ecology's 4,000 kcal/kg and 2,500 kcal a day. The floaters' own communities and consumers take their share first.

**Fibre.** A round giant's tendons are its largest material store: 157 kg/m² in a 1-km cellulose-class body, about
123,000 t, and 188 kg/m² in a 2-km body of 300 MPa fibre, about 590,000 t, both growing in canopy-like light. They come
to people when a giant dies and falls, or as the century-renewed tendons it sheds, about 1% of the store a year.

**Shade.** A floater narrower than the Sun's disk seen across its height casts no full shadow: below 93 m across at
10 km, or 186 m at 20 km, it only dims the ground (a 40-m body at 10 km removes 9% of the light at its shadow's
centre). A 200-m body at 10 km casts an umbra 107 m wide that removes the direct beam, 49% of the noon light at the
ground, leaving the bright sky; a 2-km body's umbra is 1.9 km wide. Floater shadows drift across the land like
cloud shadows.

**Habitat.** Floaters are islands of retained water and roosts in the sky: small wet communities, the collectors and
grazers the [aerial ecology](../aerial_ecology/README.md) favours around floaters, and the migrants that would bring
them phosphorus ([growth](growth.md)). People can visit and study them; a giant's lift margin could carry a station,
with its weight counted in the trim.

## Hazards

**Tethers are the hidden hazard.** Sailing tethers of 5.6–10.3 km hang weighted wings from floaters at 10 km down to 5
km or to the lowest layer near 0.7 km, 1 mm to 39 cm across ([sailing](sailing.md)). The longest reach into the sky
boats' band (0.3–1.5 km), with the wing in its middle at 0.74 km, and stay above the wing band (30–300 m). Sky lanes
need floater corridors kept clear of tethers, or tethers marked and tracked.

**Storms carry floaters out of their band.** A soft floater at 95–64% fill lifted by an updraft rises 3–31 km before
it is full ([storms](storms.md)), so storm-entrained floaters can reach the regional ships' 22–27 km and, rarely, the
liners' 35–45 km. Floaters belong below 22 km, and ships in the flight band already keep out of storms.

**A deflated floater falls slowly and drapes.** An envelope of 3–25 kg/m² falling flat sinks at 2.3–6.7 m/s near the
ground, faster aloft where the air is less dense, and takes 24–134 minutes to come down from 10–20 km. It lands
like a tent of hundreds of tonnes (a 200-m body's dry structure is 220–460 t of the 820–1,090 t it carries): a hazard
to towns and lines below, and a hydrogen fire risk while gas still escapes.

**Floaters and the summit tower.** The tower's crown at 35.6 km launches upward leaders in some 100–120 episodes a
year, mostly while storms are near (the joint synthesis); floaters keep clear of it in storms, and its port's traffic
meets their tethers.

## Sources and checks

[resources_sources.json](resources_sources.json) lists the inputs; the sky-ship hulls come from the
[sky-ship product](../sky_ships/results/sky_ships.json) and goldbeater's skin from [envelope](envelope.md).
[test_gas_biology_resources.py](test_gas_biology_resources.py) checks the shadow geometry, the gas-cell area, the
food arithmetic, the falling envelope's force balance and its fall through the density profile.
