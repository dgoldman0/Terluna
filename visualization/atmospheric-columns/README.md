# Atmospheric columns

A page for inspecting the atmosphere domain's column soundings: for each preset
(elevated fair-weather cloud, deep convection, coastal fog), the Open Moon column
beside the Earth column.

- Temperature of the air and of the lifted parcel, relative humidity, and retained
  liquid and ice condensate against altitude, with the diagnosed cloud layer shaded.
- A readout of every field at any altitude you point to.
- The model's diagnostics: scale height, lapse rate, condensation level, cloud base
  and top, buoyancy integrals, water paths, optical depth and the hydrostatic check.
- A download of any column at full resolution.

It displays [atmosphere/column](../../atmosphere/column/)'s published product,
`atmosphere/column/products/columns.json` (schema `terluna.atmosphere.column-set/1`).
Profiles are drawn from every row, and the readout samples rows by the product's
own rule, linear in altitude. The page shows the product's evidence statement,
model hash and the product's sha256. The soundings are selected experiments, not
forecasts or climatology.

```sh
node build.mjs                    # writes Open_Moon_Atmospheric_Columns.html
node --test tests/                # the page shows the product unchanged
```

`build.mjs` runs the domain's exporter first if the product is missing. The built
page is self-contained and ignored by Git.

The immersion draws clouds from the same product; this page is where its soundings
and diagnostics are read.
