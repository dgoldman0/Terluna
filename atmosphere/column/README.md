# Atmospheric column

A deterministic one-dimensional column model: pressure-coordinate soundings with
spherical hydrostatic heights, dry and liquid/ice pseudoadiabatic parcels,
diagnosed cloud support under a prescribed inversion or cap, and retained
condensate with a mixed-phase partition. `atmospheric-profile.js` continues a
column hydrostatically to 600 km as the shared thermodynamic state used by the
light-transport references. Method, equations and sources: [METHODS.md](METHODS.md).

The soundings (fair-weather cloud, deep convection, coastal fog; Moon and Earth)
are selected experiments, not forecasts or climatology, and the retained
condensate is an optical closure, not a precipitation budget.

## Who uses it

- [illumination/references](../../illumination/references/) exports the A1 shared
  profiles from it (`tools/a1/export_profiles.cjs`).
- [atmosphere/lower_air.py](../lower_air.py) reads those exports through
  `A1Profile` (for example in the megaforest wind study).
- The shoreline experience draws its column clouds from it. That runtime use will
  be replaced by baked column products, so the experience reads results rather
  than running atmospheric code.

## Run

```sh
node --test tests/weather-column.test.cjs
```

The model is JavaScript because it was written for the browser. It uses
g = 1.62 m/s², while the repository's other lunar models use GM/R² = 1.6242 m/s².
A Python implementation checked against these tests would make it a first-class
part of the atmosphere domain.
