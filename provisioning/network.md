# The Moon's network

The author asked for a continual network across the whole Moon, connected to the array, and chose a mesh of aerostats
drifting near 30 km with a small fleet of fixed posts
([register](../research/decisions.md#automation-the-network-and-food-12-october-2026)). This brief holds the
estimates behind that choice. They come from the project's products and published figures (*derived*); no network
model has been run.

## Reach

- **Absorption.** Below about 10 GHz the air absorbs a few thousandths to a hundredth of a decibel per kilometre
  (ITU-R P.676 at sea level, scaled to the Moon's air), so a 300 km link loses a few decibels at most to it.
- **The horizon.** The Moon's curvature sets the reach. From a height h the horizon lies √(2Rh + h²) away: 264 km
  from 20 km and 324 km from 30 km, against 505 and 619 km on Earth. Two nodes at 30 km see each other up to about
  650 km apart.
- **Refraction.** Near the ground the refractivity falls about 9–13 N-units per kilometre. The dry air gives 5.3: its
  density scale height is 59.5 km, from the 52 km pressure scale height and the design climate's 0.7 K/km lapse (people
  and land, table 1, in ecology). Water vapour gives 4–7: the CM1 rings' 130–240 mm column over a
  surface dewpoint near 15 °C puts its scale height at 10–19 km. The effective radius grows by 1.6–2.2%, so the radio
  horizon lies about 1% beyond the geometric one, where Earth's air adds about 15% (the 4/3-Earth rule). A trapping
  layer needs a gradient below −576 N-units per kilometre, against Earth's −157.
- **The ionosphere.** The shield passes little of the ultraviolet and X-rays that build Earth's ionosphere, and a
  rough scaling says no layer forms that reflects shortwave. Nothing in the project computes it: the electricity work
  stops at 149 km, and the conservation study counts a new ionosphere among the reasons the far side's radio quiet
  ends. The network relays everything past the horizon either way.
- **Millimetre waves and light.** At 30 km the air is still about 60% as dense as at the ground, where Earth's at
  20 km is 7%, and the column holds 5–10 times Earth's water vapour. A 300 km link at 70–90 GHz there loses tens of
  decibels, against about one at Earth's 20 km. Near-infrared lasers above the storms carry the links between nodes:
  Rayleigh scattering at 1.55 µm costs about 0.15 dB over 300 km at that height.

## Height and wind

- **The sky's layers.** Storms top out at 22 km at the median in the CM1 equatorial ring; sky reefs drift at
  10–20 km (sky ecology); aerial life is sparse above 26 km, where it is cold and dry; and the
  flight band runs from 35 to 45 km (people and land, table 1).
- **The winds** ([zonal_winds.md](../climate/gcm/zonal_winds.md)). From 20 km up the wind blows east in at least
  98% of 3-day means, at a mean 3.5 m/s at 20 km and 10.7 m/s at 40 km. At 30 km the drift towards the equator or
  the poles vanishes, and routes that start about 300 km apart keep a median of 117 km to their nearest neighbour
  after ten years. A node drifting there lives a solar day of 12.7 days on the equator (median) and 4.5–8.2 days at
  40–60° (the 30–45 km band), so its nights last 2–6 days against a fixed place's 15.
- **Lift.** Hydrogen lifts 1.31 kg a cubic metre at the ground and 0.63–0.74 kg in the flight band (people and
  land), about 0.8 kg at 30 km: ten times its lift at Earth's 20 km. A sphere about 11 m across, with a skin of
  0.15 kg/m², carries half a tonne.

## The layers of the network

- **The drifting mesh near 30 km**, above the storms and the reefs and below the flight band. With cells of 100 km
  radius, where a node stands at least 15° above the horizon, about 1,500 nodes cover the Moon, and spares fill the
  gaps the drift opens. The mesh flows east round the Moon in about three weeks at the equator. Nodes keep their
  spacing along the track by changing height, the wind growing by about 0.4 m/s each kilometre, and steer across it
  on a few watts: 1 m/s through the air, about 90 km a day, takes about 2 W for a half-tonne node. They need power
  only for what they carry, and every node passes the service ports once a lap.
- **Fixed posts at 30 km**, where a place needs steady capacity or a long watch: the shared sky tracks, the edges of
  the wild, and storms and fires while they last. Holding a half-tonne node against about 6.5 m/s takes about
  0.5 kW, which sunlight covers by day. The 15-day night needs a store, batteries at about 40% of a 5-tonne node's
  lift, or power beamed from the fleet ([questions](questions.md), item 3). Posts near 20 km would hold on about
  0.1 kW, storable for the night, but they would stand at the top of the reefs' layer, where reefs a kilometre or
  two across drift through at about 3.5 m/s, and inside many storm tops.
- **Tall structures.** The summit port's crown, at 35.6 km, is a fixed post.
- **The ground.** Fibre links the cities and towns, and small towers serve settled ground and the parks.
- **The sky towns** carry their own networks and are the mesh's large hubs; drifting low, they gather onto shared
  routes ([zonal_winds.md](../climate/gcm/zonal_winds.md)).
- **The night.** Every layer can fly dark, with infrared beacons only, so the dark sky stays dark.

## The array and space

- **The array.** The ring fleet at 19,600–21,900 km is 0.13 s away by round trip, which keeps millisecond control
  on the Moon, the Moon's AI services in the array and batch work for Earth at 2.56 s ([computing](computing.md)). The
  stack's rings, at 11.3–12.6 lunar radii and tilted up to 24°, keep 8–13 million tiles above every horizon, the
  highest within 20–25° of it at the poles ([power](power.md), *screen*). Nodes reach it by laser from above the
  weather.
- **Position and time.** The fleet tracks its tiles' orbits closely to keep them clear, and it sits about as far from
  its users as GPS satellites are from theirs. Broadcasting those orbits with time gives the Moon a navigation and
  clock service, with the aerostats filling in the geometry.
- **Imaging.** The air holds 7.1 Earth columns, so its vertical Rayleigh optical depth is about 0.7 at 550 nm and
  0.1 at 850 nm. Imaging from the array works best in the near-infrared, where a 10 m telescope resolves about 2 m.
- **No low orbits.** The air stays collisional out to its exobase at 2.5–3.8 lunar radii
  ([loss response](../atmosphere/loss_response/README.md)), so the band where Earth flies its low satellites lies
  inside the Moon's air, and the conservation study replaces low lunar orbit with high orbits. The mesh does the work
  of Earth's low satellites, and the array that of its geostationary ones.
- **The far side's radio quiet** ends with settlement; the conservation study gives radio astronomy a radio array in
  deep space behind an engineered shield
  ([D3](../research/studies/conservation/README.md#4-decisions)).

## What a node carries

- Radio for people and machines below 6 GHz, and laser links to its neighbours and the array.
- A weather station.
- Sensors for ethylene, methane and N2O, which this air keeps.
- Radar for flyers, aerophytes and towns, which the towns' right of way on the shared tracks relies on.
- Cameras looking down, resolving about 16 cm from 30 km with a 20 cm telescope at 850 nm, and thermal infrared for the
  long night.
- Samplers for the pollen and spores that stay aloft for weeks.

## Earth's precedent

Loon flew superpressure balloons at about 20 km, each serving about 5,000 km², an area 80 km across, for flights of
up to 312 days; its balloons linked over more than 100 km and in 2018 relayed data across about 1,000 km through
seven balloons. For people below, a node's cell is set by how high it stands above the horizon, so cells on the Moon
match Earth's at the same height; the Moon's curvature shortens the links between nodes.

## What remains

- A constellation model on the design winds: spacing control, coverage statistics and spares.
- The fixed posts' night power, from a store or a beam from the fleet.
- The ionosphere under the shield, with the electricity work's stage 3.
- The nodes themselves: superpressure through the drifting day, envelopes and the service ports.
- The network's own energy and spectrum.

## Sources

- ITU-R Recommendation P.676, Attenuation by atmospheric gases and related effects.
- Loon LLC ([Wikipedia](https://en.wikipedia.org/wiki/Loon_LLC)); the 312-day flight
  ([New Atlas, 2020](https://newatlas.com/telecommunications/loon-balloon-flight-record-312-days/)); the 1,000 km relay
  ([Slashdot, 2018][r]).
- The project's products named above, read as their reading rules give them.

[r]: https://tech.slashdot.org/story/18/09/13/183210/alphabets-loon-balloons-just-beamed-the-internet-across-620-miles
