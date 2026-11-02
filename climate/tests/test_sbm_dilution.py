"""The simplified Betts-Miller scheme's diluted parcel (sbm.diluted_parcel, sbm.column's dilution): with no dilution the
scheme is the undiluted one exactly; the parcel is undiluted below its condensation level and cooler above it; a dry
entraining plume in isothermal air follows its analytic solution; its steps converge; and a diluted column keeps its
moist enthalpy and water."""
from __future__ import annotations

import math

import numpy as np
import pytest

from climate.gcm import sbm

SIGMA = np.array([0.038, 0.119, 0.211, 0.317, 0.437, 0.567, 0.699, 0.823, 0.924, 0.983])
EDGES = [0.0]
for s in SIGMA:
    EDGES.append(2.0 * s - EDGES[-1])
DSIGMA = np.diff(EDGES)
PS, STEP, AKAP = 121490.0, 1800.0, 287.05 / 1004.6
C = sbm.PLASIM


def sounding(moist=1.0, warm_ground=0.0):
    """The Near-side Sea's mean column of test_sbm (gcm_site, A28R_301k_cold_moist years 0-5), top first, its vapour
    scaled by moist and its lowest level warmed by warm_ground K."""
    t = np.array([165.97, 194.74, 221.45, 243.81, 260.84, 272.61, 281.35, 289.48, 295.52, 299.16])
    q = np.array([1e-8, 2.7e-6, 6.02e-5, 4.808e-4, 1.3599e-3, 2.9967e-3, 4.7614e-3, 6.7683e-3, 1.13541e-2, 1.32844e-2])
    return t + warm_ground * (SIGMA > 0.95), q * moist


def run(t, q, **kw):
    return sbm.column(t, q, SIGMA, DSIGMA, PS, step=STEP, akap=AKAP, **kw)


def same(a, b):
    """Two outputs of sbm.column identical in every field."""
    assert a.keys() == b.keys()
    for k in a:
        if isinstance(a[k], np.ndarray):
            assert np.array_equal(a[k], b[k]), k
        else:
            assert a[k] == b[k], k


def test_zero_dilution_is_the_undiluted_scheme_exactly():
    p = PS * SIGMA
    for moist, warm in ((1.0, 0.0), (1.15, 0.0), (1.15, 2.0), (0.4, 0.0)):
        t, q = sounding(moist, warm)
        assert np.array_equal(sbm.diluted_parcel(t, q, SIGMA, PS, 0.0, AKAP), sbm.parcel(t, q, SIGMA, p, 0.0, AKAP))
        # the diluted path itself, with its mixing held above the top level, crosses every layer as parcel() does
        assert np.array_equal(sbm.diluted_parcel(t, q, SIGMA, PS, 2e-4, AKAP, base=0.5 * p[0]),
                              sbm.parcel(t, q, SIGMA, p, 0.0, AKAP))
        for kw in (dict(rh=0.5), dict(rh=0.7), dict(rh=0.5, deep_pa=20000.0, shallow_cap=True),
                   dict(rh=0.7, deep_pa=20000.0, shallow_cap=True, cin_limit=10.0)):
            same(run(t, q, dilution=0.0, **kw), run(t, q, **kw))


def test_the_diluted_parcel_is_undiluted_below_its_condensation_level_and_cooler_above():
    t, q = sounding(1.15, 2.0)
    p = PS * SIGMA
    base = sbm.condensation_level(t[-1], q[-1], p[-1], p[0], AKAP)
    plain = sbm.parcel(t, q, SIGMA, p, 0.0, AKAP)
    # the cloud of the undiluted parcel, from its condensation level to its top, where the air it takes in is cooler
    cloud = (p < base) & (np.arange(SIGMA.size) >= run(t, q, rh=0.5)['top'])
    tops = []
    for rate in (1e-4, 2e-4, 3.3e-4):
        tp = sbm.diluted_parcel(t, q, SIGMA, PS, rate, AKAP)
        below = p > base
        assert below.any() and cloud.sum() >= 4
        assert np.array_equal(tp[below], plain[below])
        assert np.all(tp[cloud] < plain[cloud])
        out = run(t, q, rh=0.5, dilution=rate)
        tops.append(len(SIGMA) if out['top'] is None else out['top'])
    # stronger mixing cools the parcel more, so its top lies no higher; at the first level above the condensation
    # level, where the saturated parcel takes in cooler air, it cools with every rise in the rate
    assert tops == sorted(tops)
    first = int(np.flatnonzero(p < base)[-1])
    cool = [sbm.diluted_parcel(t, q, SIGMA, PS, r, AKAP)[first] for r in (0.0, 1e-4, 2e-4, 3.3e-4)]
    assert cool == sorted(cool, reverse=True) and len(set(cool)) == 4


def test_a_dry_entraining_plume_in_isothermal_air_follows_the_analytic_solution():
    """Dry air at 250 K everywhere, the parcel mixing from the lowest level: in x = ln(sigma_lowest/sigma), the dry
    lift gives dT/dx = -kappa T and the mixing -eps H (T - Te), H = R Te/g, so T = Tinf + (Te - Tinf) exp(-(kappa +
    eps H) x) with Tinf = eps H Te/(kappa + eps H)."""
    te, rate = 250.0, 2e-4
    t, q = np.full(SIGMA.size, te), np.zeros(SIGMA.size)
    h = C['gascon'] * te / C['ga']
    x = np.log(SIGMA[-1] / SIGMA)
    tinf = rate * h * te / (AKAP + rate * h)
    exact = tinf + (te - tinf) * np.exp(-(AKAP + rate * h) * x)
    errors = []
    for max_mix in (0.02, 0.01, 0.005):
        tp = sbm.diluted_parcel(t, q, SIGMA, PS, rate, AKAP, base=PS, max_mix=max_mix)
        errors.append(float(np.abs(tp - exact).max()))
    # first order in the step, and within 0.05 K at the finest
    assert all(1.8 < a / b < 2.2 for a, b in zip(errors, errors[1:])), errors
    assert errors[-1] < 0.05


def test_the_diluted_parcels_steps_converge():
    t, q = sounding(1.15, 2.0)
    limit = sbm.diluted_parcel(t, q, SIGMA, PS, 2e-4, AKAP, max_mix=0.001)
    errors = [float(np.abs(sbm.diluted_parcel(t, q, SIGMA, PS, 2e-4, AKAP, max_mix=m) - limit)[2:].max())
              for m in (0.2, 0.1, 0.05)]
    assert errors[0] > errors[1] > errors[2]
    assert errors[2] < 0.1                                  # the default steps: within 0.1 K below sigma 0.2


def test_a_diluted_column_rains_its_vapour_and_keeps_its_enthalpy():
    t, q = sounding(1.15, 2.0)
    out = run(t, q, rh=0.5, dilution=2e-4)
    assert out['top'] is not None and not out['shallow']
    m = PS * DSIGMA / C['ga']
    cp = C['acpd'] * (1.0 + C['ADV'] * q)
    rain = out['rain'].sum()
    assert rain > 0.0 and rain == pytest.approx(-(out['dqdt'] * m).sum(), rel=1e-12)
    enthalpy = float(((cp * out['dtdt'] + C['ALV'] * out['dqdt']) * m).sum())
    assert abs(enthalpy) < 1e-9 * C['ALV'] * rain
    assert np.array_equal(out['reference_t'], sbm.diluted_parcel(t, q, SIGMA, PS, 2e-4, AKAP))


def test_the_diluted_parcel_and_the_mixing_option_do_not_combine():
    t, q = sounding()
    with pytest.raises(ValueError):
        run(t, q, eps=2.4, dilution=2e-4)
    assert math.isclose(sbm.MAX_MIX, 0.05)
