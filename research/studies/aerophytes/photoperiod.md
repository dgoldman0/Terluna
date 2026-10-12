# Twilight and the slow Sun

**The lunar day gives aerophytes several light histories: poleward of about 60° the clear night never falls below
1 µmol/m²/s of PAR, and holding a fixed solar hour takes only 4.3 cos(latitude) m/s of westward motion.** Light and
food are different things, though: at 60° the ground reference plant spends 34 hours a cycle below that PAR level and
runs a carbon deficit for 308 hours. [photoperiod.py](photoperiod.py) does the mean-motion geometry and the spectral
bookkeeping on committed products; [sailing](sailing.md) finds which winds carry aerophytes where.

## Moving with the Sun

**Holding mean local solar time needs a westward ground speed of v = 2π(R_Moon + h) cos φ/T_synodic: at 20 km,
4.33 m/s at the equator, 2.16 m/s at 60° and 0.75 m/s at 80°.** The period and radius are shared constants. The
navigation component turns ground speed into airspeed against the wind.

| Latitude | Westward speed at 20 km | Equinox noon Sun | Ground PAR at that height | Ground 202–1000 nm power |
|---|---:|---:|---:|---:|
| 0° | 4.328 m/s | 90° | 334.0 W/m² | 512.8 W/m² |
| 30° | 3.748 m/s | 60° | 276.6 W/m² | 425.5 W/m² |
| 60° | 2.164 m/s | 30° | 133.3 W/m² | 206.5 W/m² |
| 70° | 1.480 m/s | 20° | 81.0 W/m² | 125.8 W/m² |
| 80° | 0.752 m/s | 10° | 35.1 W/m² | 54.2 W/m² |
| 85° | 0.377 m/s | 5° | 16.1 W/m² | 24.7 W/m² |

The power columns are the surface-light product's ground values for the design shield and ground albedo 0.1, over
202–1000 nm; the older daylight solution understates low Sun against the spherical atlas the plant model uses. The
solar elevation is sin a = sin φ sin δ + cos φ cos δ cos H; with the lunar equator's 1.5424° tilt, the noon Sun at
60° spans 28.46–31.54° and at 80°, 8.46–11.54°.

**A speed error drifts the solar hour by −Δv/[(R + h) cos φ].** At 20 km a 1 m/s error crosses a ±30° hour window
from its centre in 10.65 days at the equator, 5.33 days at 60° and 1.85 days at 80°; doubling the error halves it.

**Moving west slower than the Sun lengthens both day and night.** At 60° and 20 km a steady 2 m/s westward drift gives
a 389.9-day relative solar period, with 194.95 days of geometric day and as long a night at zero declination.
Matching the Sun's angular pace holds the chosen hour; matching it in dim twilight holds the twilight. Moving east
shortens both, the effect the eastward winds aloft give ([sailing](sailing.md)).

## Three meanings of night

**PAR darkness, a plant's carbon deficit and visual darkness are different intervals.** The plant product's darkness
is PAR below 1 µmol/m²/s; its maintenance deficit is the time its photosynthesis does not cover its upkeep.

| Site at 60° | Strategy | PAR below 1 | Carbon deficit | Carbon reserve | Leaf / stem–root respiration |
|---|---|---:|---:|---:|---:|
| Near side | Earth-like | 33.714 h | 308.718 h | 34.081 g C/m² | 74.041 / 35.789 g C/m² |
| Far side | Earth-like | 33.714 h | 307.487 h | 34.344 g C/m² | 71.795 / 34.592 g C/m² |
| Near side | Idling at a quarter | 33.714 h | 225.048 h | 5.538 g C/m² | 52.337 / 25.394 g C/m² |
| Far side | Idling at a quarter | 33.714 h | 225.171 h | 5.629 g C/m² | 50.044 / 24.173 g C/m² |

Hours and carbon are per synodic cycle per square metre of ground, for the plant model's LAI 5 stand and its C3,
tissue-temperature and respiration assumptions. Near and far PAR darkness match because the plant model's light has
no Earthlight; their differences come from the temperature traces. The product's historical field
`respiration_in_the_dark` integrates respiration while the Sun is below the geometric horizon, and this component
names it `sun_below_horizon_respiration_g_c_m2_cycle`.

**Even the stationary equatorial plant sees 238 hours of PAR darkness, well short of the 354-hour geometric night.** At
30° it is 217 hours. The equatorial twilight keeps PAR above 50 for 9.60 hours after sunset, above 10 for 31.51 hours
and above 1 for 57.96 hours.

**PAR stays at 1 µmol/m²/s down to a solar depression of 29.63°, so the clear night never drops below it poleward of
60.37°** (61.91° allowing the ±1.5424° declination). At 70° the equinox midnight Sun is 20° down and the ground gets
5.14 µmol/m²/s. The analytic duty gives 33.81 hours below 1 at 60° against the plant model's sampled 33.71; the
product keeps the plant model's value.

## Twilight energy

**The ground's twilight carries 0.21523 W of PAR and 0.28941 W over 202–1000 nm per µmol of PAR,** using the plant
producer's 1° diffuse-black-ground spectrum scaled by its stored twilight PAR (`twilight_energy()`).

| Sun | PAR photons | PAR energy | 202–1000 nm energy |
|---|---:|---:|---:|
| −30° | 0.907 µmol/m²/s | 0.195 W/m² | 0.262 W/m² |
| −20° | 5.137 µmol/m²/s | 1.106 W/m² | 1.487 W/m² |
| −15° | 11.794 µmol/m²/s | 2.539 W/m² | 3.413 W/m² |
| −10° | 25.288 µmol/m²/s | 5.443 W/m² | 7.319 W/m² |
| −8° | 33.431 µmol/m²/s | 7.196 W/m² | 9.675 W/m² |

The calendar product's solved ground illuminance is 24.60 lux at −30°, 180.83 at −20° and 1,192.91 at −10°, a
visual quantity kept separate. At an assumed 1% conversion on the same band, the −10° row would supply 0.0732 W/m² of
hydrogen, so geometric night need not stop hydrogen production; a gas load of 0.07 W/m² needs 7 W/m² incident at 1%.
Each candidate organism integrates its own compensation, upkeep and hydrogen response along its light history; the
table stops at −35°, and aloft the illumination solver gives the twilit sky.

Sources are the committed products listed in [photoperiod_sources.json](photoperiod_sources.json).
