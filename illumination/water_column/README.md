# Water column

The colour of seawater and the light under it. [model.py](model.py) adds the
absorption, scattering and backscattering of a water's constituents and gives its
remote-sensing reflectance and the diffuse attenuation of daylight:

- pure seawater: absorption from Pope and Fry (1997) and Kou et al. (1993),
  scattering from Zhang et al. (2009) at 20 °C and 38.4 g/kg;
- phytoplankton and the particles that vary with them: absorption from Bricaud et
  al. (1998), scattering and backscattering from Loisel and Morel (1998) and
  Morel and Maritorena (2001);
- coloured dissolved organic matter, with the mean spectral slope of Babin et al.
  (2003);
- suspended regolith fines: the single-scattering albedo of the finest fraction
  (under 10 µm) of Apollo soils, from the Lunar Soil Characterization Consortium's
  reflectance through Hapke's isotropic model and equivalent slab, carried into
  water, for grains of a stated size and density;
- reflectance by Lee et al. (2002), attenuation by Lee et al. (2005).

At Earth's chlorophyll levels of 0.1–3 mg/m³ its reflectance ratios give back
the chlorophyll within a factor of two through NASA's OC4 algorithm, and its
attenuation at 490 nm in clear water follows Morel's relation within a quarter.

The tables are fetched and hash-checked, never committed:

```sh
python -m illumination.water_column.fetch_inputs --download
```

[inputs.json](inputs.json) records their sources and credits. What the Open Moon's
seas carry is the biosphere's
[design guess](../../biosphere/living_water/waters.json).
