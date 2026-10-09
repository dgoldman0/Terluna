# Twilight and colonies that follow the Sun

The long lunar cycle does not imply14.77days of biological darkness for every colony.
The stored ground calculation already has prolonged diffuse twilight, especially at
high latitude. Motion at a few metres per second can also hold a colony near a chosen
solar hour. Both possibilities belong in the viability calculation. Their value depends
on the organism's light requirements and the motion available along an actual route.

`photoperiod.evaluate()` provides the new component. It reads committed products and
performs mean-motion geometry and spectral bookkeeping. No weather campaign, optical
rerun or new physiological model was performed. The environment component's half-cycle
water/thermal examples remain sensitivities; they are not a universal dark-duration rule.

## Moving with the Sun

At latitude φ and altitude h, mean Sun-following requires westward ground speed

\[
v_\odot={2\pi(R_{Moon}+h)\cos\phi\over T_{synodic}}.
\]

The period and radius come from shared constants. This is ground-relative motion.
The [navigation component](navigation.md) must account for the vector difference
between that motion and the surrounding wind, then pay the propulsion cost.

| Latitude | Westward speed at20km | Equinox noon solar height | Ground PAR at that height | Ground202–1000nm power |
|---|---:|---:|---:|---:|
| 0° | 4.328m/s | 90° | 334.0W/m² | 512.8W/m² |
| 30° | 3.748m/s | 60° | 276.6W/m² | 425.5W/m² |
| 60° | 2.164m/s | 30° | 133.3W/m² | 206.5W/m² |
| 70° | 1.480m/s | 20° | 81.0W/m² | 125.8W/m² |
| 80° | 0.752m/s | 10° | 35.1W/m² | 54.2W/m² |
| 85° | 0.377m/s | 5° | 16.1W/m² | 24.7W/m² |

The power columns are **ground references**, not irradiance at20km. They use the
surface-light product's design shield and neutral ground albedo0.1. The optical band
stops at1000nm, so it is not total solar power. The older daylight product also
underestimates low-Sun light relative to the spherical atlas used by the plant model.
It is not justified to hold a200W/m² mean input constant across every latitude.

The local-horizontal solar elevation follows

\[
\sin a=\sin\phi\sin\delta+\cos\phi\cos\delta\cos H.
\]

Noon has H=0. A colony following another hour angle has less light; holding latitude
and solar hour also leaves seasonal declination changes. The shared lunar inclination
is1.5424°, and the output includes declinations0 and±1.5424°. At60° latitude the
noon elevation is28.46–31.54°. At80° it is8.46–11.54°. Those changes matter for
small light margins. Altitude does not add directly to elevation above the local
horizontal; its effect on horizon visibility and light transport requires a separate
calculation. The kinematics neglect real ephemeris speed variations.

A speed mismatch Δv changes solar hour angle at−Δv/[(R+h)cosφ]. `drift_window`
reports time from the centre of a±30° hour-angle window to its edge. At20km a1m/s
mismatch takes10.65days at the equator,5.33days at60° and1.85days at80°.
Doubling the error halves that residence. These are angular windows, not cloud
residence times or physiological thresholds.

**Steady westward drift does not necessarily remove night.** If it is slower than
Sun-following, it lengthens both the day and the night. At60° and20km, a constant
2m/s westward ground motion has a389.9day relative solar period:194.95days of
geometric day and194.95days of geometric night at zero declination. That is a
mean-motion extrapolation, not a prediction that a wind persists for a year.
Only matching the solar angular motion holds the chosen hour in this calculation;
matching while already in dim twilight also keeps the colony in dim twilight.
`passive_wind_windows` includes0 and±1/2/5m/s westward wind scenarios without
claiming that any is a controllable lunar route.

## Three different meanings of night

The primary plant product defines darkness as PAR below1µmol/m²/s. Its maintenance
deficit is the interval during which photosynthesis does not cover that particular
ground plant's upkeep. Visual darkness is different again. Neither threshold defines
the duty of a microbial or engineered hydrogen-production pathway.

| Site at60° | Strategy | PAR below1 | Maintenance deficit | Carbon reserve | Leaf / stem-root respiration |
|---|---|---:|---:|---:|---:|
| Near side | Earth-like | 33.714h | 308.718h | 34.081gC/m² | 74.041 /35.789gC/m² |
| Far side | Earth-like | 33.714h | 307.487h | 34.344gC/m² | 71.795 /34.592gC/m² |
| Near side | Quarter idle | 33.714h | 225.048h | 5.538gC/m² | 52.337 /25.394gC/m² |
| Far side | Quarter idle | 33.714h | 225.171h | 5.629gC/m² | 50.044 /24.173gC/m² |

Hours and carbon amounts are per synodic cycle, area is ground footprint, and the
stand has LAI5. These are its own C3/tissue-temperature/respiration assumptions,
not a floater's compensation threshold. Near/far PAR darkness is identical because
`plant.py` does not include Earthlight; its side differences come from temperature
traces. Its historical `respiration_in_the_dark` field actually integrates respiration
while the Sun is below the geometric horizon. This component names that quantity
`sun_below_horizon_respiration_g_c_m2_cycle` to prevent confusion with PAR darkness.

At the equator the stored PAR-dark interval is237.721h; at30° it is217.296h.
The twilight table keeps PAR above50 for9.60h after equatorial sunset, above10 for
31.51h and above1 for57.96h. Thus even the stationary equatorial plant does not
experience354h below its PAR-dark threshold.

Interpolating that same clear-ground twilight table gives PAR=1 at a solar elevation
of−29.63°. In the ideal zero-declination geometry, latitudes poleward of60.37°
never fall below that level; allowing the±1.5424° declination sensitivity moves
the all-season geometric bound to about61.91°. This does not establish year-round
carbon surplus. At70° latitude, the equinox midnight Sun is−20° down and the
table gives only5.14µmol/m²/s. Clouds and a colony's altitude remain unresolved.
The continuous analytic threshold calculation gives33.81h below1 at60°, versus
the plant model's sampled33.71h; the inherited value is preserved rather than replaced.

## Twilight energy and trait-specific requirements

The current calendar product predicts ground point-source solar diffuse illuminance
of24.60lux at−30°,180.83lux at−20°, and1192.91lux at−10°. Those are visual
quantities from the solved column. The product does not expose the corresponding
ground spectral energy, and lux cannot be converted to H2 production without a
spectrum and action response. This study leaves the calendar lux separate.

For a reproducible **inherited approximation**, `twilight_energy()` follows the plant
producer's spectral assumption: the1° diffuse-black-ground spectrum, scaled by its
stored twilight PAR. It is not a newly solved twilight spectrum. Per PAR micromole,
that shape carries0.21523W in PAR and0.28941W over202–1000nm:

| Sun elevation | PAR photons | PAR energy |202–1000nm energy |
|---|---:|---:|---:|
| −30° | 0.907µmol/m²/s | 0.195W/m² | 0.262W/m² |
| −20° | 5.137µmol/m²/s | 1.106W/m² | 1.487W/m² |
| −15° | 11.794µmol/m²/s | 2.539W/m² | 3.413W/m² |
| −10° | 25.288µmol/m²/s | 5.443W/m² | 7.319W/m² |
| −8° | 33.431µmol/m²/s | 7.196W/m² | 9.675W/m² |

A gas-production load of0.07W/m² requires7W/m² incident at1% conversion,
provided efficiency, spectral band and collector area have the same definition.
Under a deliberately imposed1% efficiency on this202–1000nm twilight proxy,
the−10° row would supply0.0732W/m². That demonstrates why geometric night cannot
automatically force H2 production to zero. It does not transfer a full-solar laboratory
efficiency to a new spectrum, nor show that enough remaining energy is available
for organism maintenance, growth and locomotion. Those uses must share the budget.

The parent should integrate each candidate's own C3 compensation, maintenance and
H2 flux response along its light history. Low-barrier gas replacement can have a
different useful-light interval from net plant growth. Flux remains nonzero below a
chosen operating threshold. The table stops at−35°; deeper values are unresolved
here rather than silently set to zero. The old daylight and twilight approximations
must not be spliced into a supposedly continuous, validated altitude light curve.

Sun-following, persistent diffuse light, periods of reduced metabolism and reserves
can therefore be combined as explicit candidate strategies. Their joint route must
still close carbon, water, heat, lift gas and motive-power budgets under the same
weather and seasonal illumination history.
