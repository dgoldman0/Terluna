"""How the runs' year files hold their 3-day records, and the humidity means to read from them.

The runs write one netCDF file a model year (MOST.<year>.nc), 120 records of 3-day means in low-I/O mode (NLOWIO 1,
NSTPW 143): PlaSim adds each output field to its mean every step and writes the means at the end of each interval
(outmod.f90, outaccu, outgp and outsp). Temperature, surface pressure, the winds, cloud, precipitable water and the
fluxes are true means, and the gridpoint diagnostics (codes 21-41) are single steps at each interval's end.

Specific humidity is not a mean. With spectral humidity (NQSPEC 1) the field PlaSim adds for code 133, sqout, is
recomputed only at output steps (plasim.f90, `if (mod(nstep,nafter)==0 .and. nqspec == 1)`), so each record of hus is
the humidity at the end of the interval before, held through the interval. sqout is not in the restart either, so each
year file's first record holds about 2% of the vapour. Each interval's mean is taken here as the average of the
humidity at its start and its end, the record and the next. Over the fifteen-level baseline's second year that brings
the column's vapour within 5.7 kg/m2 of the precipitable-water means (rms over cells and records), against 10.4 for the
record itself, with no bias (14.2 to 8.2 in the Kuo run). A year file's first record, which is empty, and its last,
whose end lies in the next file, have no such mean, so readers take records 1 to n-2 of every field. A run that writes
the humidity's true mean (hus_mean, code 233) has every record.
"""
from __future__ import annotations

import numpy as np

HUS_MEAN = 'hus_mean'         # code 233: the gridpoint humidity PlaSim accumulates every step (aadq), as a mean


def records(d) -> slice:
    """The records of a year file (an open netCDF4 dataset) whose humidity means are known: every record where the run
    writes hus_mean, otherwise the second to the second-last."""
    return slice(None) if HUS_MEAN in d.variables else slice(1, -1)


def interval_means(hus) -> np.ndarray:
    """Interval means of specific humidity from records that each hold the humidity at their interval's start (time
    first): the average of each record and the next, for records 1 to n-2."""
    hus = np.asarray(hus, float)
    return 0.5 * (hus[1:-1] + hus[2:])


def humidity(d) -> np.ndarray:
    """3-day means of specific humidity (kg/kg) over the records records(d) of a year file, (time, lev, lat, lon)."""
    if HUS_MEAN in d.variables:
        return np.asarray(d[HUS_MEAN][:], float)
    return interval_means(d['hus'][:])


def fields(d, names) -> dict:
    """The named fields of a year file over records(d), each a float array with time first; 'hus' is the humidity's
    means (humidity) and fields without a time axis come whole."""
    keep = records(d)
    out = {}
    for name in names:
        if name == 'hus':
            out[name] = humidity(d)
        elif d[name].dimensions and d[name].dimensions[0] == 'time':
            out[name] = np.asarray(d[name][:], float)[keep]
        else:
            out[name] = np.asarray(d[name][:], float)
    return out
