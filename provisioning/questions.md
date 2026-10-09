# Provisioning's open questions

In order of what each would settle for the whole Open Moon. Each names what decides it and how it could be answered
on this laptop. The findings behind them are in the [README](README.md).

1. **The ring fleet's infrared under the planned climate.** The screen puts 5.7–9.9 W/m² of the films' own infrared
   on the Moon, which the design climate leaves out. The author decided on 9 October that the investigated climate is
   the planned one, so the glow is answered with added solar protection, by redirecting it away from the Moon, or
   both.
   - **What decides it:** how the stack's layers pass heat among themselves and to the Moon, with the films' thermal
     transmission and reflection, which the two-layer estimate leaves out; where in the column the 175 K emission is
     absorbed and how much it warms the ground; whether it fills part of the upper air's CO2 cooling band, which
     the loss response relies on; and how the array's industry and collectors change the films' absorption, the
     dimming and where their heat leaves (README, finding 1).
   - **The means to size:** films that absorb less, with the ultraviolet reflected before the titania absorbs it; a
     heat mirror on Moon-facing faces, on the layer nearest the Moon where the stack allows it, adding little
     absorption (a sparse metal mesh or metasurface the first candidate); strong emitters on the outer faces; and
     added dimming for what remains.
   - **How:** the ring stack's radiative transfer with the protection domain's optical tables (2–125 µm) on the
     fleet-light model's geometry, minutes to hours; a downward band from the fleet into the radiative–convective
     column; then the added dimming that holds the planned climate's energy budget. The climate programme stays
     paused.
   - **Whose:** the shield's work on its films.
2. **The distribution of people among the ways of living.** The [shared comparison cases](../shared/scenarios/population.json)
   set lunar surface and aerial settlements beside Earth and the orbital habitats. The
   [inhabited-volume screen](../research/studies/inhabited_volume/README.md) makes lift, wet mass, floor area,
   projected footprint, holding power and domestic recovery explicit. The
   [ways-of-living study](../research/studies/ways_of_living/README.md) gives eighteen settings with their physical
   situation, their costs per billion and the cases' totals by mix: every mix carries every case at its own split,
   with food, the dark night and the sky towns' hydrogen the largest costs.
   - **What decides it:** the author's choice among the costs, after how people live
     ([daily life](../research/studies/daily_life/README.md)), steps five and six of the approved order.
   - **What remains to size:**
     - town and hull designs, with their structure, gas cells and wet mass;
     - drifting districts' routes, their steering against the gathering the design winds show (onto one track within
       months near the ground, over the poles within years at 40–60 km), their storm routing and rendezvous
       ([climate/gcm/zonal_winds.md](../climate/gcm/zonal_winds.md));
     - rain at the sky towns' height;
     - crop water in the plant model;
     - nutrient return;
     - the dark-night requirement's value;
     - health over a lifetime at a sixth of Earth's gravity.
3. **The power system through the night.** Fusion is a planning assumption and a conditional technology. The
   alternatives:
   - a grid from the day side, with 9–19% losses over 4,000–5,500 km;
   - storage at eight times Earth's pumped hydro for the metropolis alone;
   - power from orbit, including microwaves from the fleet's night half.

   Beaming keeps the night dark but closes 4,000–7,000 km² of receivers to people and sky ships for 200 GW.

   **How:** a network model of sources, storage, lines and reserves through the lunar cycle, and the summit city's
   power plan that its studies ask for.
4. **The food system.**
   - Diets, and staple crops beyond the day fruit.
   - Protein from microbial production on firm power.
   - Water for farms outside the tropics.
   - Irrigation beyond the rain: at the high end of the diet range ten billion lunar residents irrigate 2.2 million
     km² with 8,700 km³ a year from the fresh seas, 1.5 TW of lifting ([ways of living](../research/studies/ways_of_living/README.md), finding 2).
   - The monthly harvest and its stores.
   - The land within reach of each city.

   With the ecology work it must also make no methane and keep ethylene low enough for cereals.

   **How:** land and energy per diet on the canopy model's yields and the GCM's rain, beside the ecology work's
   water model.
5. **Freight between the surface and orbit.** Surface mass drivers fail in the air, so freight to orbit has to start
   from altitude, a tether or an elevator:
   - the high platforms near 70 km and the crown at 35.6 km are candidate starts;
   - an elevator meets the ring stack;
   - a launch rule in the spirit of S7 is needed, alongside the resources domain's rule for arrivals.

   Orbital industry then draws on space resources.

   **How:** an ascent study through the Open Moon's column from several release heights, and the upper air's
   tolerance for exhaust from the loss response.
6. **Water, sanitation and the return of nutrients.**
   - Treatment that carries the whole job of disinfection.
   - Sewage as a main phosphorus return path.
   - Methane captured from waste.
   - Water systems that ride out the storms.

   **How:** per-person flows from the literature, set against the metropolis and the ecology work's nutrient budget.
7. **Building materials of the far side.** The aluminium, glass, ceramic and calcium-silicate materials the dry
   highland gives; where the steel comes from; the metropolis's 4–34 Gt of stocks and their embodied energy. The
   resources domain's geochemical map feeds this.
8. **Transit at a sixth of gravity.**
   - Rail and metro designed around 0.17–0.25 m/s² for standing riders, or seated travel.
   - Walking and running at 0.16 g: people switch to running at 1.42 m/s in lunar gravity (De Witt et al. 2014), and
     braking on a dry path holds about 1 m/s², so stopping distances are six times Earth's
     ([daily life](../research/studies/daily_life/README.md)).
   - The sky boats' traffic system, which the sky-fleet study left open.
9. **Computing.** Supply and demand are in [computing](computing.md): the light outruns any computing the array can
   build, mass binds first and heat placement next. Its open questions, in order: how much computing the civilization
   wants, digital minds included; processor mass and life at the ring radius; where large computing sits and its
   share of the planned climate's heat budget; the path in energy per operation; the management twin's resolution by
   place; services to Earth; and collectors that dim the Moon while feeding computing. Each is weighed by what it does
   to the Moon's heat and the glow (README, finding 1). With the planned `research/array-industry` branch.
10. **Work, maintenance and the economy.** Maintenance runs 2–4% of replacement value a year; who does it, how much is
    automated, and how work is shared belong to the human paper. Provisioning supplies the hours and energy.
11. **Lamps and the dark night.** The summit's night has about 190 hours that need lamps. At the lunar design's
    lighting each billion people light 0.79–1.42 million km² of sky, so 6–15 billion lunar residents light 20–91% of
    the dry land ([ways of living](../research/studies/ways_of_living/README.md), finding 4). Sleep needs at most 1 lux
    indoors (Brown et al. 2022), which shutters give anywhere, so the people's side of the night requirement lies in
    the outdoor dark sky ([daily life](../research/studies/daily_life/README.md)). Their energy, the skyglow of a city
    of 100 million in this scattering air and every other setting's lighting are weighed against the dark-night
    requirement once its value is set.
12. **Health and life support.**
    - Vitamin D and other biological needs under the filtered spectrum, using diet/adaptation and controlled local UV; a narrow natural UV-B option remains a separately evaluated candidate under the cold-atmosphere requirements.
    - Artificial-gravity facilities, carried open from the July canon. A room at Earth's gravity on the surface is a
      55 m wheel at 4 rpm with its floor banked 80.5°, and facilities for 1–8 hours a day take 420–3,330 km² of
      rotating floor per billion residents ([ways of living](../research/studies/ways_of_living/README.md)). How much
      time at Earth's gravity a lifetime needs is unknown (Clément et al. 2015). No human has lived in partial gravity
      beyond 75 hours; the evidence for children, older residents and return to Earth is gathered in
      [daily life](../research/studies/daily_life/README.md), section 5.
    - Life support for enclosed places, such as the undersea community.
