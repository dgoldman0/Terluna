"""The spectral plume model (climate/gcm/plume_spectrum.py): with one plume it is the second redesign's reference
(sbm_redesign2.zbp_reference) and with no entrainment the finely integrated undiluted parcel; every plume keeps its
moist static energy budget; the environment it sets stands at the plume neutral at each level, with the plumes still
rising warmer and those stopped colder; tops taken from a mass-flux profile start at its peak and end where it ends;
its deficit below the undiluted parcel grows with entrainment; 16 steps a layer agree with 64 within 0.2 K; and one
column costs what the module's account states."""
from __future__ import annotations

import time

import numpy as np
import pytest

from climate.gcm import plume_spectrum as psm
from climate.gcm import sbm
from climate.gcm import sbm_redesign2 as sr2
from climate.tests import test_sbm_closure as tsc      # PlaSim's fifteen levels, the design run's constants, soundings

SIGMA, AKAP, C = tsc.SIGMA, tsc.AKAP, tsc.C
# A spread like the readings' central one (CM1's rates by height at 0.70 of saturation): rates from 2.4e-4 down to
# 2.9e-5 per metre, uniform in log, and CM1's afternoon cloud tops above the condensation level (metres, shares).
TOPS = psm.tops_closure([0.0, 9282.7, 15351.0, 21076.0, 47076.0, 69076.0], [0.0, 0.41334, 0.5, 0.9, 0.99, 1.0])
SPREAD = psm.log_spectrum(2.85e-5, 2.4e-4)


def columns():
    out = []
    for name in tsc.SOUNDINGS:
        for moist, warm in ((1.0, 0.0), (1.1, 2.0)):
            t, q, ps = tsc.sounding(name, moist, warm)
            out.append((name, t, q, ps))
    return out


def base_of(t, q, ps):
    p = ps * SIGMA
    return sbm.condensation_level(t[-1], q[-1], p[-1], p[0], AKAP, C)


def fine(t, q, ps, nsub=16):
    return sbm.diluted_parcel(t, q, SIGMA, ps, 0.0, AKAP, C, max_mix=sbm.MAX_MIX, min_steps=nsub)


@pytest.mark.parametrize('name,t,q,ps', columns())
def test_one_plume_is_the_second_redesigns_reference(name, t, q, ps):
    for eps in (0.25, 0.5, 1.0):
        ref = sr2.zbp_reference(t[-1], q[-1], SIGMA, ps, C, AKAP, base_of(t, q, ps), eps=eps, rh=0.7)
        o = psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, [eps], family='inverse_z', rh=0.7)
        assert np.max(np.abs(o['t_env'] - ref)) < 1e-9, (name, eps)


@pytest.mark.parametrize('name,t,q,ps', columns())
def test_no_entrainment_gives_the_finely_integrated_undiluted_parcel(name, t, q, ps):
    one = psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, [0.0])
    spread = psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, np.zeros(8), tops=TOPS)
    assert np.max(np.abs(one['t_env'] - fine(t, q, ps))) < 1e-9
    assert np.max(np.abs(spread['t_env'] - fine(t, q, ps))) < 1e-9


@pytest.mark.parametrize('name,t,q,ps', columns())
def test_each_plume_keeps_its_moist_static_energy_budget(name, t, q, ps):
    for params, family in ((SPREAD[0], 'constant'), (np.array([2.0, 1.0, 0.5]), 'inverse_z')):
        o = psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, params, family=family, tops=TOPS)
        assert np.max(np.abs(o['budget_residual'])) < 1e-6                 # rounding, J/kg
        assert np.max(np.abs(o['budget_lift'])) < 5.0                      # the dry lift's third-order term, J/kg
        assert np.all(o['budget_entrainment'] < -1000.0)                   # entrainment lowers each plume's energy


@pytest.mark.parametrize('name,t,q,ps', columns())
def test_the_environment_stands_at_the_plume_neutral_at_each_level(name, t, q, ps):
    o = psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, SPREAD[0], tops=TOPS)
    assert o['sorting'] < 1e-9                         # still rising: no colder than it; stopped: no warmer
    assert np.all(np.diff(o['top_z']) >= 0.0)          # the less entraining plumes stop higher
    assert np.all(np.diff(o['mass_flux']) >= 0.0)      # top first: the spread's mass flux falls with height
    above = o['z_env'] > o['z_base']
    ok = above & np.isfinite(o['rh_romps'])
    assert np.any(ok) and np.all((o['rh_romps'][ok] > 0.0) & (o['rh_romps'][ok] < 1.0))
    one = psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, SPREAD[0][:1], tops=None)
    assert np.all(o['t_env'] >= one['t_env'] - 1e-9)   # warmer than its most entraining plume alone


def test_tops_from_a_mass_flux_profile_start_at_its_peak_and_end_where_it_ends():
    z = np.array([5e3, 10e3, 20e3, 30e3, 40e3, 50e3, 60e3, 70e3])
    m = np.array([0.002, 0.010, 0.020, 0.016, 0.018, 0.012, 0.004, 0.0])
    tops = psm.tops_from_mass_flux(z, m, base=6e3)
    assert tops(10e3) == 0.0 and tops(14e3) == 0.0                     # nothing stops below the peak (20 km)
    assert tops(24e3) == pytest.approx(0.2) and tops(34e3) == pytest.approx(0.2)   # kept nondecreasing
    assert tops(44e3) == pytest.approx(0.4) and tops(64e3) == 1.0
    t, q, ps = tsc.sounding('L')
    o = psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, SPREAD[0], tops=tops)
    assert o['sorting'] < 1e-9 and np.all(np.diff(o['top_z']) >= 0.0)


@pytest.mark.parametrize('name,t,q,ps', columns())
def test_the_deficit_grows_with_entrainment(name, t, q, ps):
    yard = fine(t, q, ps)
    above = SIGMA * ps < base_of(t, q, ps)
    above[-1] = False
    last = None
    for eps in (0.0, 0.25e-4, 0.5e-4, 1.0e-4, 2.0e-4):
        d = yard - psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, [eps])['t_env']
        if last is not None:
            assert np.all(d[above] >= last[above]) and np.mean(d[above]) > np.mean(last[above]) + 0.1
        last = d
    last = None
    for scale in (0.5, 1.0, 2.0):
        d = yard - psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, SPREAD[0] * scale, tops=TOPS)['t_env']
        if last is not None:
            assert np.all(d[above] >= last[above]) and np.mean(d[above]) > np.mean(last[above]) + 0.1
        last = d


@pytest.mark.parametrize('name,t,q,ps', columns()[::2])
def test_sixteen_steps_a_layer_agree_with_sixty_four_within_a_fifth_of_a_kelvin(name, t, q, ps):
    for params, family in ((SPREAD[0], 'constant'), (np.array([2.0, 1.0, 0.5]), 'inverse_z')):
        a = psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, params, family=family, tops=TOPS, nsub=16)
        b = psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, params, family=family, tops=TOPS, nsub=64)
        assert np.max(np.abs(a['t_env'] - b['t_env'])) < 0.2


def test_a_column_of_thirty_two_plumes_costs_a_few_hundredths_of_a_second():
    """One fifteen-level column with the 32-plume spread in this pure-Python copy: about 0.03-0.04 s on the author's
    laptop beside other work (4 November 2026); the bound leaves room for a loaded machine."""
    t, q, ps = tsc.sounding('L')
    psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, SPREAD[0], tops=TOPS)
    start = time.perf_counter()
    for _ in range(5):
        o = psm.equilibrium(t[-1], q[-1], SIGMA, ps, C, AKAP, SPREAD[0], tops=TOPS)
    cost = (time.perf_counter() - start) / 5
    assert o['steps'] <= 16 * (SIGMA.size - 1)
    assert cost < 0.5
