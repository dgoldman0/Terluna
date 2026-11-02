"""The simplified Betts-Miller scheme's relaxation of the buoyant layer alone (sbm.column's relax_from_lfc): off, the
scheme is Frierson's exactly; on, the levels below the parcel's level of free convection keep their state, the
buoyant layer alone relaxes and heats, and the column keeps its moist enthalpy and its water."""
from __future__ import annotations

import numpy as np

from climate.tests.test_sbm_dilution import C, DSIGMA, PS, run, same, sounding


def deep_columns():
    """Soundings on which the scheme convects deep, raining, with a buoyant layer above an unbuoyant lower part."""
    out = []
    for moist, warm in ((1.15, 2.0), (1.3, 3.0), (1.2, 1.0)):
        t, q = sounding(moist, warm)
        for kw in (dict(rh=0.5), dict(rh=0.7, dilution=2e-4), dict(rh=0.83, dilution=2e-4, tau=18 * 3600.0)):
            o = run(t, q, **kw)
            if o['top'] is not None and not o['shallow'] and o['rain'].sum() > 0.0:
                out.append((t, q, kw))
    return out


def test_off_is_the_scheme_as_it_was():
    for moist, warm in ((1.0, 0.0), (1.15, 2.0), (0.4, 0.0)):
        t, q = sounding(moist, warm)
        for kw in (dict(rh=0.7), dict(rh=0.7, dilution=2e-4, deep_pa=20000.0, shallow_cap=True)):
            same(run(t, q, relax_from_lfc=False, **kw), run(t, q, **kw))


def test_the_layers_below_the_level_of_free_convection_keep_their_state():
    checked = 0
    for t, q, kw in deep_columns():
        o = run(t, q, relax_from_lfc=True, **kw)
        tp = o['reference_t']
        lfc = next(k for k in range(t.size - 2, -1, -1) if tp[k] > t[k])
        if o['top'] is None or o['shallow']:
            continue
        checked += 1
        below = np.arange(t.size) > lfc
        assert np.all(o['dtdt'][below] == 0.0) and np.all(o['dqdt'][below] == 0.0) and np.all(o['rain'][below] == 0.0)
        assert np.any(o['dtdt'][~below] != 0.0)
    assert checked >= 5, f'{checked} deep test columns'


def test_the_column_keeps_its_moist_enthalpy_and_water():
    checked = 0
    for t, q, kw in deep_columns():
        o = run(t, q, relax_from_lfc=True, **kw)
        if o['top'] is None or o['shallow']:
            continue
        checked += 1
        m = PS * DSIGMA / C['ga']
        cp = C['acpd'] * (1.0 + C['ADV'] * q)
        rain = o['rain'].sum()
        assert rain > 0.0
        assert abs((-o['dqdt'] * m).sum() - rain) <= 1e-12 * rain
        assert abs((cp * o['dtdt'] * m).sum() - C['ALV'] * rain) <= 1e-9 * C['ALV'] * rain
    assert checked >= 5, f'{checked} deep test columns'

