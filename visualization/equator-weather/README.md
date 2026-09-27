# Equator weather

A page for reading what the cloud-resolving model gives for the Open Moon's
equator through the month-long day, for the chosen design (28% seas with lakes,
the 5% dimmer shield):

- the lunar day by local time over land and sea: air and humidity at 2 m,
  rain and cloud, with a readout at any local time, and a table from midnight
  through the day;
- the heaviest storm of the lunar day as a height–distance section, with the
  flight band marked;
- where and when it rained all the way round the ring, with local noon;
- the circulation that follows the Sun: wind along the equator, rising and
  sinking, and cloud, by local time and height;
- winds by height, from 10 m to 70 km, and the storm and water statistics.

It displays [climate/crm](../../climate/crm/)'s products,
`climate/results/crm/ring_ring.json` (schema `terluna.climate.crm-ring/1`) and
`climate/results/crm/ring_ring_fields.npz` (`terluna.climate.crm-ring-fields/1`),
both written by `climate/crm/ring_analysis.py` and kept in Git. Every number in
the page's text is filled from them, and the page states their evidence, reading
rules and sha256. The ring is two-dimensional, flat, with a placeholder land
surface and its upper air held to the GCM; the page says so.

```sh
python3 visualization/equator-weather/build.py         # writes Open_Moon_Equator_Weather.html
python3 -m pytest visualization/equator-weather/tests   # the page shows the products unchanged
```

The built page is self-contained and ignored by Git. Its tests run on their own,
like the other Python viewers' here.
